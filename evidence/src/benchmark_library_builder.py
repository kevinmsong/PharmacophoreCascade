"""
Build a retrospective GLP-1R active-vs-decoy benchmark library.

This module keeps the two benchmark-construction steps separate:
1. Curate in-domain direct small-molecule agonists from ChEMBL.
2. Match property-similar ZINC decoys from the local screening pool.
"""
from __future__ import annotations

from dataclasses import dataclass
import gzip
from pathlib import Path
from typing import Iterable, Iterator, Sequence

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

CHEMBL_GLP1R_TARGET = "CHEMBL1784"
CHEMBL_ACTIVITY_URL = "https://www.ebi.ac.uk/chembl/api/data/activity.json"

DIRECT_AGONIST_TERMS = (
    "agonist activity",
    "glp-1r-mediated agonist activity",
    "stimulation of camp",
    "camp accumulation",
)
EXCLUDED_ASSAY_TERMS = (
    "allosteric",
    "pam",
    "antagonist mode",
    "ec20",
    "glp-1 assay with glp1 endogenous agonist",
)


@dataclass(frozen=True)
class MoleculeDescriptor:
    ligand_id: str
    smiles: str
    canonical_smiles: str
    murcko_scaffold: str
    mw: float
    logp: float
    hbd: int
    hba: int
    tpsa: float
    rotatable_bonds: int
    formal_charge: int
    fingerprint: object


def _standardize_smiles(smiles: str) -> Chem.Mol | None:
    mol = Chem.MolFromSmiles(str(smiles).strip())
    if mol is None:
        return None
    return Chem.RemoveHs(mol)


def describe_molecule(ligand_id: str, smiles: str) -> MoleculeDescriptor | None:
    mol = _standardize_smiles(smiles)
    if mol is None:
        return None
    return MoleculeDescriptor(
        ligand_id=str(ligand_id),
        smiles=str(smiles),
        canonical_smiles=Chem.MolToSmiles(mol),
        murcko_scaffold=MurckoScaffold.MurckoScaffoldSmiles(mol=mol),
        mw=float(Descriptors.MolWt(mol)),
        logp=float(Crippen.MolLogP(mol)),
        hbd=int(Lipinski.NumHDonors(mol)),
        hba=int(Lipinski.NumHAcceptors(mol)),
        tpsa=float(rdMolDescriptors.CalcTPSA(mol)),
        rotatable_bonds=int(rdMolDescriptors.CalcNumRotatableBonds(mol)),
        formal_charge=int(Chem.GetFormalCharge(mol)),
        fingerprint=AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048),
    )


def descriptor_to_record(desc: MoleculeDescriptor) -> dict[str, object]:
    return {
        "ligand_id": desc.ligand_id,
        "smiles": desc.smiles,
        "canonical_smiles": desc.canonical_smiles,
        "murcko_scaffold": desc.murcko_scaffold,
        "mw": desc.mw,
        "logp": desc.logp,
        "hbd": desc.hbd,
        "hba": desc.hba,
        "tpsa": desc.tpsa,
        "rotatable_bonds": desc.rotatable_bonds,
        "formal_charge": desc.formal_charge,
    }


def is_in_stage0_domain(desc: MoleculeDescriptor) -> bool:
    return (
        200.0 <= desc.mw <= 500.0
        and -1.0 <= desc.logp <= 5.0
        and desc.hbd <= 5
        and desc.hba <= 10
    )


def fetch_chembl_activities(
    target_chembl_id: str,
    standard_type: str = "EC50",
    limit: int = 1000,
) -> pd.DataFrame:
    """Fetch ChEMBL activity rows for any target / standard-type pair."""
    import requests

    records: list[dict[str, object]] = []
    offset = 0
    while True:
        response = requests.get(
            CHEMBL_ACTIVITY_URL,
            params={
                "target_chembl_id": target_chembl_id,
                "standard_type": standard_type,
                "relation__exact": "=",
                "limit": limit,
                "offset": offset,
            },
            timeout=60,
        )
        response.raise_for_status()
        payload = response.json()
        activities = payload.get("activities", [])
        if not activities:
            break
        records.extend(activities)
        offset += len(activities)
        if offset >= int(payload.get("page_meta", {}).get("total_count", 0)):
            break
    return pd.DataFrame(records)


def fetch_glp1r_chembl_activities(
    standard_type: str = "EC50",
    limit: int = 1000,
) -> pd.DataFrame:
    """GLP-1R-specific wrapper retained for backward compatibility."""
    return fetch_chembl_activities(CHEMBL_GLP1R_TARGET, standard_type, limit)


def curate_in_domain_actives(
    target_chembl_id: str,
    *,
    standard_types: Sequence[str] = ("EC50",),
    default_potency_type: str = "EC50",
    require_assay_terms: bool = True,
    assay_terms: Sequence[str] = DIRECT_AGONIST_TERMS,
    excluded_terms: Sequence[str] = EXCLUDED_ASSAY_TERMS,
    scaffold_limit: int = 10,
    potency_threshold_nM: float = 15000.0,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Curate in-domain small-molecule actives for an arbitrary ChEMBL target.

    Generalizes the GLP-1R curation: one or more ChEMBL ``standard_type`` queries
    are pooled, then filtered by excluded-assay terms, an optional required
    assay-term whitelist (e.g. agonist context for peptide-GPCRs; disabled for
    interface-inhibitor targets such as MDM2-p53), a potency threshold, and the
    Stage-0 physicochemical design space, before canonical-SMILES and scaffold
    deduplication.
    """
    frames = [fetch_chembl_activities(target_chembl_id, st) for st in standard_types]
    frames = [frame for frame in frames if not frame.empty]
    raw = (
        pd.concat(frames, ignore_index=True)
        if frames
        else pd.DataFrame(columns=["canonical_smiles"])
    )
    curated_rows: list[dict[str, object]] = []
    excluded_rows: list[dict[str, object]] = []

    for row in raw.to_dict(orient="records"):
        assay_desc = str(row.get("assay_description") or "")
        assay_desc_l = assay_desc.lower()
        chembl_id = str(row.get("molecule_chembl_id") or "")
        display_id = chembl_id or f"activity_{row.get('activity_id')}"
        smiles = str(row.get("canonical_smiles") or "").strip()
        potency_raw = row.get("standard_value")
        try:
            potency_value = float(potency_raw)
        except (TypeError, ValueError):
            potency_value = np.nan

        if not smiles:
            excluded_rows.append(
                {
                    "ligand_id": display_id,
                    "source": "ChEMBL",
                    "potency_type": row.get("standard_type") or default_potency_type,
                    "potency_value": potency_value,
                    "reason": "missing_smiles",
                }
            )
            continue

        descriptor = describe_molecule(display_id, smiles)
        if descriptor is None:
            excluded_rows.append(
                {
                    "ligand_id": display_id,
                    "source": "ChEMBL",
                    "potency_type": row.get("standard_type") or default_potency_type,
                    "potency_value": potency_value,
                    "reason": "invalid_smiles",
                }
            )
            continue

        if any(term in assay_desc_l for term in excluded_terms):
            excluded_rows.append(
                {
                    "ligand_id": display_id,
                    "canonical_smiles": descriptor.canonical_smiles,
                    "source": "ChEMBL",
                    "potency_type": row.get("standard_type") or default_potency_type,
                    "potency_value": potency_value,
                    "reason": "excluded_assay_context",
                    "assay_description": assay_desc,
                }
            )
            continue

        if require_assay_terms and not any(term in assay_desc_l for term in assay_terms):
            excluded_rows.append(
                {
                    "ligand_id": display_id,
                    "canonical_smiles": descriptor.canonical_smiles,
                    "source": "ChEMBL",
                    "potency_type": row.get("standard_type") or default_potency_type,
                    "potency_value": potency_value,
                    "reason": "not_direct_agonist_assay",
                    "assay_description": assay_desc,
                }
            )
            continue

        if not np.isfinite(potency_value) or potency_value > potency_threshold_nM:
            excluded_rows.append(
                {
                    "ligand_id": display_id,
                    "canonical_smiles": descriptor.canonical_smiles,
                    "source": "ChEMBL",
                    "potency_type": row.get("standard_type") or default_potency_type,
                    "potency_value": potency_value,
                    "reason": "potency_threshold",
                    "assay_description": assay_desc,
                }
            )
            continue

        if not is_in_stage0_domain(descriptor):
            excluded_rows.append(
                {
                    "ligand_id": display_id,
                    "canonical_smiles": descriptor.canonical_smiles,
                    "source": "ChEMBL",
                    "potency_type": row.get("standard_type") or default_potency_type,
                    "potency_value": potency_value,
                    "reason": "out_of_stage0_domain",
                    "assay_description": assay_desc,
                    **descriptor_to_record(descriptor),
                }
            )
            continue

        curated_rows.append(
            {
                "ligand_id": display_id,
                "chembl_id": chembl_id,
                "source": f"ChEMBL:{row.get('document_chembl_id') or 'unknown_document'}",
                "pref_name": row.get("molecule_pref_name") or chembl_id,
                "potency_value": potency_value,
                "potency_type": row.get("standard_type") or default_potency_type,
                "document_year": row.get("document_year"),
                "assay_description": assay_desc,
                "in_domain": 1,
                **descriptor_to_record(descriptor),
            }
        )

    curated_df = pd.DataFrame(curated_rows)
    if curated_df.empty:
        return curated_df, pd.DataFrame(excluded_rows)

    curated_df = (
        curated_df.sort_values(["canonical_smiles", "potency_value", "document_year"])
        .drop_duplicates("canonical_smiles", keep="first")
        .sort_values(["potency_value", "document_year", "ligand_id"])
        .reset_index(drop=True)
    )

    if scaffold_limit > 0:
        included = []
        excluded_duplicates = []
        seen_scaffolds: set[str] = set()
        for row in curated_df.to_dict(orient="records"):
            scaffold = str(row.get("murcko_scaffold") or "")
            if scaffold in seen_scaffolds:
                excluded_duplicates.append(
                    {
                        "ligand_id": row["ligand_id"],
                        "canonical_smiles": row["canonical_smiles"],
                        "source": row["source"],
                        "potency_type": row["potency_type"],
                        "potency_value": row["potency_value"],
                        "reason": "scaffold_deduplicated",
                        "murcko_scaffold": scaffold,
                    }
                )
                continue
            included.append(row)
            seen_scaffolds.add(scaffold)
            if len(included) >= scaffold_limit:
                break
        curated_df = pd.DataFrame(included)
        excluded_rows.extend(excluded_duplicates)

    curated_df = curated_df.reset_index(drop=True)
    curated_df["label"] = "active"
    curated_df["chemotype"] = [
        f"scaffold_{index+1}" for index in range(len(curated_df))
    ]
    return curated_df, pd.DataFrame(excluded_rows).reset_index(drop=True)


def curate_in_domain_glp1r_actives(
    scaffold_limit: int = 10,
    potency_threshold_nM: float = 15000.0,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """GLP-1R curation (EC50 direct agonists). Backward-compatible wrapper."""
    return curate_in_domain_actives(
        CHEMBL_GLP1R_TARGET,
        standard_types=("EC50",),
        default_potency_type="EC50",
        require_assay_terms=True,
        assay_terms=DIRECT_AGONIST_TERMS,
        excluded_terms=EXCLUDED_ASSAY_TERMS,
        scaffold_limit=scaffold_limit,
        potency_threshold_nM=potency_threshold_nM,
    )


def iter_zinc_smiles(zinc_paths: Sequence[Path]) -> Iterator[tuple[str, str]]:
    for path in zinc_paths:
        with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                fields = line.split()
                if len(fields) < 2:
                    continue
                smiles, ligand_id = fields[0], fields[1]
                yield ligand_id, smiles


def load_zinc_candidate_pool(
    zinc_paths: Sequence[Path],
    max_candidates: int | None = None,
) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    seen_ids: set[str] = set()
    for ligand_id, smiles in iter_zinc_smiles(zinc_paths):
        if ligand_id in seen_ids:
            continue
        descriptor = describe_molecule(ligand_id, smiles)
        if descriptor is None or not is_in_stage0_domain(descriptor):
            continue
        seen_ids.add(ligand_id)
        records.append(descriptor_to_record(descriptor))
        if max_candidates and len(records) >= max_candidates:
            break
    return pd.DataFrame(records)


def _property_distance(active: pd.Series, candidate: pd.Series) -> float:
    scales = {
        "mw": 100.0,
        "logp": 1.5,
        "hba": 2.0,
        "hbd": 1.0,
        "tpsa": 25.0,
        "rotatable_bonds": 2.0,
    }
    distance = 0.0
    for field, scale in scales.items():
        distance += abs(float(candidate[field]) - float(active[field])) / scale
    return float(distance)


def match_decoys_from_pool(
    actives_df: pd.DataFrame,
    candidate_pool_df: pd.DataFrame,
    decoys_per_active: int = 30,
    tanimoto_threshold: float = 0.35,
    candidate_cap_per_active: int = 300,
    random_state: int | None = None,
) -> pd.DataFrame:
    """Assign unique ZINC decoys to actives using property matching and global dissimilarity.

    With ``random_state=None`` the strictly nearest property-matched decoys are
    chosen (deterministic; reproduces the published library). With an integer
    ``random_state`` the nearest neighbours are shuffled before assignment, so
    each seed produces an independent but equally property-matched decoy set for
    the multiple-decoy-set robustness analysis.
    """
    if actives_df.empty or candidate_pool_df.empty:
        return pd.DataFrame()

    active_descriptors = {
        row["ligand_id"]: describe_molecule(row["ligand_id"], row["smiles"])
        for _, row in actives_df.iterrows()
    }
    active_descriptors = {k: v for k, v in active_descriptors.items() if v is not None}
    active_rows = actives_df.set_index("ligand_id")

    candidate_records: dict[str, list[dict[str, object]]] = {ligand_id: [] for ligand_id in active_descriptors}

    for row in candidate_pool_df.to_dict(orient="records"):
        ligand_id = str(row["ligand_id"])
        smiles = str(row["smiles"])
        descriptor = describe_molecule(ligand_id, smiles)
        if descriptor is None or not is_in_stage0_domain(descriptor):
            continue

        max_similarity = 0.0
        for active_desc in active_descriptors.values():
            max_similarity = max(
                max_similarity,
                float(DataStructs.TanimotoSimilarity(descriptor.fingerprint, active_desc.fingerprint)),
            )
            if max_similarity >= tanimoto_threshold:
                break
        if max_similarity >= tanimoto_threshold:
            continue

        candidate_row = pd.Series(descriptor_to_record(descriptor))
        for active_id, active_desc in active_descriptors.items():
            if descriptor.formal_charge != active_desc.formal_charge:
                continue
            if descriptor.murcko_scaffold == active_desc.murcko_scaffold:
                continue
            distance = _property_distance(active_rows.loc[active_id], candidate_row)
            bucket = candidate_records[active_id]
            bucket.append(
                {
                    **descriptor_to_record(descriptor),
                    "matched_active_ligand_id": active_id,
                    "property_distance": distance,
                    "max_active_tanimoto": max_similarity,
                }
            )
            bucket.sort(key=lambda rec: (rec["property_distance"], rec["ligand_id"]))
            if len(bucket) > candidate_cap_per_active:
                del bucket[candidate_cap_per_active:]

    if random_state is not None:
        rng = np.random.default_rng(random_state)
        shuffle_pool = max(decoys_per_active * 4, decoys_per_active + 40)
        for active_id, bucket in candidate_records.items():
            head = bucket[:shuffle_pool]
            rng.shuffle(head)
            candidate_records[active_id] = head + bucket[shuffle_pool:]

    assigned_rows: list[dict[str, object]] = []
    used_decoys: set[str] = set()
    active_ids = list(active_descriptors)
    for round_index in range(decoys_per_active):
        for active_id in active_ids:
            bucket = candidate_records[active_id]
            chosen = None
            while bucket:
                candidate = bucket.pop(0)
                if candidate["ligand_id"] not in used_decoys:
                    chosen = candidate
                    break
            if chosen is None:
                continue
            used_decoys.add(str(chosen["ligand_id"]))
            chosen["label"] = "decoy"
            chosen["source"] = "local_zinc_pool"
            chosen["potency_value"] = np.nan
            chosen["potency_type"] = ""
            chosen["in_domain"] = 1
            chosen["chemotype"] = f"matched_to_{active_id}"
            assigned_rows.append(chosen)

    return pd.DataFrame(assigned_rows).reset_index(drop=True)


def build_benchmark_library(
    actives_df: pd.DataFrame,
    decoys_df: pd.DataFrame,
) -> pd.DataFrame:
    keep_columns = [
        "ligand_id",
        "smiles",
        "label",
        "source",
        "potency_value",
        "potency_type",
        "murcko_scaffold",
        "in_domain",
        "canonical_smiles",
        "chemotype",
        "mw",
        "logp",
        "hbd",
        "hba",
        "tpsa",
        "rotatable_bonds",
        "formal_charge",
        "matched_active_ligand_id",
        "property_distance",
        "max_active_tanimoto",
        "pref_name",
        "chembl_id",
        "document_year",
        "assay_description",
    ]
    merged = pd.concat([actives_df, decoys_df], ignore_index=True, sort=False)
    for column in keep_columns:
        if column not in merged.columns:
            merged[column] = np.nan
    return merged[keep_columns].copy()
