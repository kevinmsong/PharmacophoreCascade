#!/usr/bin/env python3
"""A readable reference implementation of Stage 0-2 cascade scoring.

This module exists so that a third party can check that the equations printed
in the manuscript are the equations the screen actually computes. It is a
plain, unoptimised reimplementation written directly from those equations and
from the published parameter values -- it shares no code with the production
screening engine, and it is not intended to be fast.

What is implemented
-------------------
* **Stage 0** (eq. 1): standardization and the property envelope; structural
  alerts are recorded without exclusion in the current configuration.
* **Stage 1** (eq. 2): the hotspot-weighted fraction ``H`` and the three-part
  gate.
* **Stage 2** (eq. 4): the typed pair-hash catalogues and the recall,
  precision, and F1 overlap ``O_pair``.
* **Cascade score** (eq. 3): ``100 x (0.4 H + 0.6 O_pair)``.

What is not implemented
-----------------------
Stage 3 and the terminal native branch both require 3D conformer generation
and alignment. They are deterministic given a conformer set but depend on
RDKit's embedding, so reproducing them exactly from a script is not
meaningful; the released per-molecule Stage-3 and native scores support
reanalysis of those stages instead.

Usage
-----
    python reference_implementation/stage_scoring_reference.py \
        --pharmacophore maps/pharmacophore_rigorous.json \
        --smiles "CC(C)(C)Nc1nc2cc(Cl)c(Cl)cc2nc1S(C)(=O)=O"
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")

# --------------------------------------------------------------------------
# Published parameters (Methods, Table of cascade settings)
# --------------------------------------------------------------------------
FEATURE_TYPES = ("anion", "cation", "donor", "acceptor", "aromatic", "hydrophobe")

#: Per-family caps on how many ligand features enter the pair catalogue.
BASE_TYPE_CAPS = {
    "anion": 2, "cation": 2, "donor": 4, "acceptor": 6, "aromatic": 6, "hydrophobe": 6,
}

#: Receptor feature family -> the ligand family that complements it.
RECEPTOR_TO_QUERY_FAMILY = {
    "positive": "anion", "negative": "cation", "hbd": "acceptor",
    "hba": "donor", "aromatic": "aromatic", "hydrophobic": "hydrophobe",
}

#: Chemically related families that earn partial credit.
QUERY_FAMILY_COMPATIBILITY = {
    "anion": ("acceptor",), "cation": ("donor",), "donor": (), "acceptor": (),
    "aromatic": (), "hydrophobe": ("aromatic",),
}

COMPATIBLE_MATCH_CREDIT = 0.5  # gamma in eq. (2); also the pair-expansion credit
CASCADE_WEIGHTS = {"hotspot": 0.4, "pair_hash": 0.6}  # eq. (3)

#: precision_5bin: bond-path and through-space bin edges (Methods, Stage 2).
GRAPH_EDGES = (2.0, 3.0, 5.0, 7.0)
RECEPTOR_EDGES = (3.5, 5.5, 7.5, 10.0)

#: Stage-0 property envelope (eq. 1). Rotatable-bond count is reported by the
#: screen but is not a gate condition.
PROPERTY_BOUNDS = {"mw": (None, 500.0), "logp": (-1.0, 5.0), "hbd": (None, 5), "hba": (None, 10)}

REACTIVE_SMARTS = {
    "azide": "[N-]=[N+]=N",
    "acyl_halide": "[CX3](=[OX1])[F,Cl,Br,I]",
    "isocyanate": "N=C=O",
    "isothiocyanate": "N=C=S",
    "sulfonyl_halide": "[SX4](=[OX1])(=[OX1])[F,Cl,Br,I]",
    "epoxide": "[OX2r3]1[#6r3][#6r3]1",
    "aldehyde": "[CX3H1](=O)[#6]",
    "michael_acceptor": "[C,c]=[C,c]-[C](=O)[O,N,S]",
}

_params = FilterCatalogParams()
for _c in (FilterCatalogParams.FilterCatalogs.PAINS_A,
           FilterCatalogParams.FilterCatalogs.PAINS_B,
           FilterCatalogParams.FilterCatalogs.PAINS_C):
    _params.AddCatalog(_c)
PAINS_CATALOG = FilterCatalog(_params)
REACTIVE_PATTERNS = {k: Chem.MolFromSmarts(v) for k, v in REACTIVE_SMARTS.items()}
_LARGEST_FRAGMENT = rdMolStandardize.LargestFragmentChooser()


def compute_type_caps(n_features: int) -> dict[str, int]:
    """Scale the per-family caps with pharmacophore size."""
    base_total = sum(BASE_TYPE_CAPS.values())
    if n_features <= base_total:
        return dict(BASE_TYPE_CAPS)
    scale = n_features / base_total
    return {k: max(v, round(v * scale)) for k, v in BASE_TYPE_CAPS.items()}


# --------------------------------------------------------------------------
# Stage 0 - eq. (1)
# --------------------------------------------------------------------------
def standardize(smiles: str) -> Chem.Mol | None:
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None
    try:
        return _LARGEST_FRAGMENT.choose(mol)
    except Exception:
        return None


def stage0(smiles: str, exclude_alerts: bool = False) -> dict:
    """Evaluate retained property limits; report alerts without excluding by default."""
    mol = standardize(smiles)
    if mol is None:
        return {"pass": False, "reason": "unparseable", "mol": None}

    props = {
        "mw": Descriptors.MolWt(mol),
        "logp": Crippen.MolLogP(mol),
        "hbd": Descriptors.NumHDonors(mol),
        "hba": Descriptors.NumHAcceptors(mol),
        "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
    }
    violations = [
        f"{k}={props[k]:.2f}"
        for k, (lo, hi) in PROPERTY_BOUNDS.items()
        if (lo is not None and props[k] < lo) or (hi is not None and props[k] > hi)
    ]
    if violations:
        return {"pass": False, "reason": "property: " + ", ".join(violations),
                "mol": mol, "properties": props}

    # The property gate is applied first; only survivors are alert-screened.
    alerts = []
    hit = PAINS_CATALOG.GetFirstMatch(mol)
    if hit is not None:
        alerts.append(f"PAINS:{hit.GetDescription()}")
    for label, patt in REACTIVE_PATTERNS.items():
        if patt is not None and mol.HasSubstructMatch(patt):
            alerts.append(f"reactive:{label}")
    if alerts and exclude_alerts:
        return {"pass": False, "reason": "alert: " + "; ".join(alerts),
                "mol": mol, "properties": props, "alerts": alerts}

    return {"pass": True, "reason": "", "mol": mol, "properties": props, "alerts": alerts}


# --------------------------------------------------------------------------
# Typed ligand features
# --------------------------------------------------------------------------
def atom_family(atom: Chem.Atom) -> str | None:
    n, q = atom.GetAtomicNum(), atom.GetFormalCharge()
    if n == 7 and q > 0:
        return "cation"
    if n == 8 and q < 0:
        return "anion"
    if n in (7, 8) and atom.GetTotalNumHs() > 0:
        return "donor"
    if n in (7, 8):
        return "acceptor"
    if atom.GetIsAromatic():
        return "aromatic"
    if n in (6, 16):
        return "hydrophobe"
    return None


def atom_priority(atom: Chem.Atom, family: str) -> float:
    """Priority score used to rank candidate features within a family."""
    n, q, deg = atom.GetAtomicNum(), atom.GetFormalCharge(), atom.GetDegree()
    in_ring, aromatic = int(atom.IsInRing()), int(atom.GetIsAromatic())
    hetero = sum(1 for nb in atom.GetNeighbors() if nb.GetAtomicNum() not in (1, 6))
    if family == "cation":
        return 200 + 50 * max(q, 0) + 10 * hetero + 5 * deg
    if family == "anion":
        return 200 + 50 * max(-q, 0) + 10 * hetero + 5 * deg
    if family == "donor":
        return (120 + 12 * atom.GetTotalNumHs() + 8 * hetero + 5 * in_ring
                + (8 if n == 7 else 6 if n == 8 else 0))
    if family == "acceptor":
        return (120 + 8 * hetero + 4 * in_ring
                + (10 if n == 8 else 8 if n == 7 else 0) - 2 * atom.GetTotalNumHs())
    if family == "aromatic":
        return 110 + 15 * aromatic + 8 * in_ring + 2 * hetero
    return 100 + 8 * in_ring + 3 * deg - 2 * hetero


def _aromatic_rings(mol: Chem.Mol) -> dict[int, set[int]]:
    memb: dict[int, set[int]] = defaultdict(set)
    for ring_i, atoms in enumerate(mol.GetRingInfo().AtomRings()):
        if all(mol.GetAtomWithIdx(a).GetIsAromatic() for a in atoms):
            for a in atoms:
                memb[int(a)].add(ring_i)
    return memb


def _redundant(mol, node, selected, rings) -> bool:
    """Suppress features that duplicate an already-selected nearby feature."""
    fam, idx = node["family"], node["atom_idx"]
    node_rings = rings.get(idx, set())
    for sel in selected:
        if sel["family"] != fam:
            continue
        sidx = sel["atom_idx"]
        if fam == "aromatic":
            if node_rings & rings.get(sidx, set()):
                return True
            continue
        d = len(Chem.GetShortestPath(mol, idx, sidx)) - 1
        if fam in {"anion", "cation", "donor", "acceptor"} and d <= 2:
            return True
        if fam == "hydrophobe":
            if d <= 1:
                return True
            if node_rings and node_rings & rings.get(sidx, set()):
                return True
    return False


def typed_features(mol: Chem.Mol) -> list[dict]:
    """Typed, deduplicated ligand pharmacophoric features, highest priority first."""
    raw = []
    for atom in mol.GetAtoms():
        fam = atom_family(atom)
        if fam is None:
            continue
        raw.append({"atom_idx": atom.GetIdx(), "family": fam,
                    "priority": atom_priority(atom, fam)})
    raw.sort(key=lambda d: (-d["priority"], d["atom_idx"]))

    rings = _aromatic_rings(mol)
    kept: list[dict] = []
    for node in raw:
        if not _redundant(mol, node, kept, rings):
            kept.append(node)
    return kept


# --------------------------------------------------------------------------
# Stage 1 - eq. (2)
# --------------------------------------------------------------------------
def stage1(features: list[dict], hotspots: list[dict], required_groups: tuple[str, ...],
           min_exact: int = 3, min_groups: int = 2) -> dict:
    """Hotspot-weighted fraction H and the three-part Stage-1 gate.

    Each receptor hotspot consumes at most one ligand feature: an exact
    family match takes full weight, otherwise a compatible family takes
    ``gamma = 0.5``. Features are consumed highest-priority first.
    """
    available: dict[str, list[dict]] = defaultdict(list)
    for f in features:
        available[f["family"]].append(f)
    for fam in available:
        available[fam].sort(key=lambda d: (-d["priority"], d["atom_idx"]))

    def take(families: tuple[str, ...]) -> dict | None:
        best, best_fam = None, None
        for fam in families:
            if available.get(fam):
                cand = available[fam][0]
                if best is None or cand["priority"] > best["priority"]:
                    best, best_fam = cand, fam
        if best is not None:
            available[best_fam].pop(0)
        return best

    total_weight = sum(float(h["weight"]) for h in hotspots)
    matched_weight = 0.0
    exact_matches = compatible_matches = 0
    exact_priority_residues: set[int] = set()
    exact_priority_groups: set[str] = set()

    for h in hotspots:
        fam = h["query_family"]
        w = float(h["weight"])
        hit = take((fam,))
        if hit is not None:
            exact_matches += 1
            matched_weight += w
            if h.get("priority_feature"):
                exact_priority_residues.add(int(h.get("resnum") or 0))
                if h.get("query_group"):
                    exact_priority_groups.add(str(h["query_group"]))
            continue
        hit = take(tuple(QUERY_FAMILY_COMPATIBILITY.get(fam, ())))
        if hit is not None:
            compatible_matches += 1
            matched_weight += w * COMPATIBLE_MATCH_CREDIT

    H = matched_weight / max(total_weight, 1e-9)
    groups_ok = (not required_groups) or bool(exact_priority_groups & set(required_groups))
    return {
        "H": H,
        "hotspot_weighted_pct": 100.0 * H,
        "exact_matches": exact_matches,
        "compatible_matches": compatible_matches,
        "exact_residue_count": len(exact_priority_residues),
        "group_count": len(exact_priority_groups),
        "required_groups_pass": groups_ok,
        "gate_pass": (len(exact_priority_residues) >= min_exact
                      and len(exact_priority_groups) >= min_groups
                      and groups_ok),
    }


# --------------------------------------------------------------------------
# Stage 2 - eq. (4)
# --------------------------------------------------------------------------
def _bucket(distance: float, edges: tuple[float, ...]) -> int:
    for i, e in enumerate(edges):
        if distance <= e:
            return i
    return len(edges)


def _pair_key(a: str, b: str, bucket: int) -> tuple:
    return tuple(sorted((a, b))) + (bucket,)


def _expanded_keys(fa: str, fb: str, bucket: int) -> tuple:
    """Compatible pair keys a receptor pair can match, with their credit."""
    out: dict[tuple, float] = {}
    for ca in (fa,) + QUERY_FAMILY_COMPATIBILITY.get(fa, ()):
        credit_a = 1.0 if ca == fa else COMPATIBLE_MATCH_CREDIT
        for cb in (fb,) + QUERY_FAMILY_COMPATIBILITY.get(fb, ()):
            credit_b = 1.0 if cb == fb else COMPATIBLE_MATCH_CREDIT
            key = _pair_key(ca, cb, bucket)
            out[key] = max(out.get(key, 0.0), credit_a * credit_b)
    return tuple(sorted(out.items(), key=lambda kv: (-kv[1], kv[0])))


def ligand_pair_catalogue(mol: Chem.Mol, features: list[dict],
                          caps: dict[str, int]) -> tuple[dict, float]:
    """Bin every retained feature pair by its shortest bond path."""
    by_fam: dict[str, list[dict]] = defaultdict(list)
    for f in features:
        by_fam[f["family"]].append(f)
    selected = []
    for fam in FEATURE_TYPES:
        selected.extend(by_fam.get(fam, [])[: caps[fam]])
    if len(selected) < 2:
        return {}, 0.0

    selected.sort(key=lambda d: (d["atom_idx"], d["family"]))
    dmat = np.asarray(Chem.GetDistanceMatrix(mol), dtype=float)
    max_priority = max(f["priority"] for f in selected)
    pairs: dict[tuple, dict] = defaultdict(lambda: {"count": 0, "weight_sum": 0.0})
    total = 0.0
    for i, left in enumerate(selected):
        for right in selected[i + 1:]:
            d = dmat[left["atom_idx"], right["atom_idx"]]
            if d <= 0:
                continue
            key = _pair_key(left["family"], right["family"], _bucket(d, GRAPH_EDGES))
            w = 0.5 * (left["priority"] / max_priority + right["priority"] / max_priority)
            pairs[key]["count"] += 1
            pairs[key]["weight_sum"] += w
            total += w
    return dict(pairs), total


def receptor_pair_query(hotspots: list[dict]) -> tuple[dict, float]:
    """Bin every receptor feature pair by its through-space distance."""
    pos = np.asarray([[h["x"], h["y"], h["z"]] for h in hotspots], dtype=float)
    pairs: dict[tuple, dict] = defaultdict(lambda: {"count": 0, "weight_sum": 0.0})
    total = 0.0
    for i, left in enumerate(hotspots):
        for j in range(i + 1, len(hotspots)):
            right = hotspots[j]
            d = float(np.linalg.norm(pos[i] - pos[j]))
            key = _pair_key(left["query_family"], right["query_family"],
                            _bucket(d, RECEPTOR_EDGES))
            w = 0.5 * (float(left["weight"]) + float(right["weight"]))
            pairs[key]["count"] += 1
            pairs[key]["weight_sum"] += w
            total += w

    finalized = {}
    for (fa, fb, bucket), stats in pairs.items():
        finalized[(fa, fb, bucket)] = {
            "count": stats["count"],
            "weight_sum": stats["weight_sum"],
            "candidate_keys": _expanded_keys(fa, fb, bucket),
        }
    return finalized, total


def stage2(ligand_pairs: dict, total_ligand_w: float,
           query_pairs: dict, total_query_w: float) -> dict:
    """Pair-hash recall, precision, and the F1 overlap O_pair.

    Each receptor pair bucket is credited by its best-matching ligand bucket
    (and vice versa for precision), scaled by the fraction of the bucket's
    count that is covered and by the family-compatibility credit. Recall and
    precision are then combined as their harmonic mean.
    """
    if not ligand_pairs or not query_pairs or total_query_w <= 0 or total_ligand_w <= 0:
        return {"pair_hash_overlap": 0.0, "pair_hash_recall": 0.0, "pair_hash_precision": 0.0}

    matched_query_w = 0.0
    reverse: dict[tuple, list[tuple[tuple, float]]] = defaultdict(list)
    for qkey, qstats in query_pairs.items():
        best = 0.0
        qc, qw = qstats["count"], qstats["weight_sum"]
        for ckey, credit in qstats["candidate_keys"]:
            reverse[ckey].append((qkey, credit))
            lstats = ligand_pairs.get(ckey)
            if not lstats or lstats["count"] <= 0:
                continue
            frac = min(lstats["count"], qc) / max(qc, 1)
            best = max(best, qw * frac * credit)
        matched_query_w += best

    matched_ligand_w = 0.0
    for lkey, lstats in ligand_pairs.items():
        lc = lstats["count"]
        if lc <= 0:
            continue
        lw = lstats["weight_sum"]
        best = 0.0
        for qkey, credit in reverse.get(lkey, []):
            qc = query_pairs[qkey]["count"]
            frac = min(qc, lc) / max(lc, 1)
            best = max(best, lw * frac * credit)
        matched_ligand_w += best

    recall = matched_query_w / max(total_query_w, 1e-9)
    precision = matched_ligand_w / max(total_ligand_w, 1e-9)
    overlap = (2 * precision * recall / max(precision + recall, 1e-9)
               if precision > 0 and recall > 0 else 0.0)
    return {"pair_hash_overlap": overlap, "pair_hash_recall": recall,
            "pair_hash_precision": precision}


# --------------------------------------------------------------------------
# Pharmacophore loading and the full Stage 0-2 evaluation
# --------------------------------------------------------------------------
NON_CURATED_PRIORITY = 9999


def _query_priority(f: dict) -> int:
    """Rank used to order receptor features: curated first, then native-supported."""
    curated = int(f.get("curated_priority") or 0)
    if curated > 0:
        return curated
    if f.get("native_supported"):
        return 1000 + int(f.get("native_support_rank") or 9999)
    return NON_CURATED_PRIORITY


def _sort_key(f: dict) -> tuple:
    return (-float(f["weight"]), _query_priority(f), str(f["type"]),
            int(f.get("resnum") or 0), str(f.get("atom") or ""))


def _collapse(features: list[dict]) -> list[dict]:
    """Keep one feature per (residue, complementary family), the heaviest."""
    best: dict[tuple, dict] = {}
    for f in features:
        fam = RECEPTOR_TO_QUERY_FAMILY.get(f["type"])
        if fam is None:
            continue
        key = (int(f.get("resnum") or 0), fam)
        if key not in best or _sort_key(f) < _sort_key(best[key]):
            best[key] = f
    return sorted(best.values(), key=_sort_key)


def select_query_features(features: list[dict], n_keep: int,
                          native_support_max_residues: int = 8) -> list[dict]:
    """Select the receptor features that seed Stages 1 and 2.

    Curated contact features come first, ordered by their curated priority;
    if fewer than ``n_keep`` remain, up to ``native_support_max_residues``
    native-supported residues are added. The screen runs this with backfill
    disabled, so the selection can be smaller than ``n_keep`` -- for the
    curated GLP-1R model it is, and that is the behaviour reproduced here.
    """
    collapsed = _collapse(features)
    curated = [f for f in collapsed if int(f.get("curated_priority") or 0) > 0]
    curated.sort(key=lambda f: (_query_priority(f), -float(f["weight"]), str(f["type"]),
                                int(f.get("resnum") or 0), str(f.get("atom") or "")))
    selected = curated[:n_keep]
    seen = {(int(f.get("resnum") or 0), str(f["type"]), str(f.get("atom") or ""))
            for f in selected}

    if native_support_max_residues > 0 and len(selected) < n_keep:
        supported = [f for f in collapsed
                     if int(f.get("curated_priority") or 0) <= 0 and f.get("native_supported")]
        supported.sort(key=lambda f: (int(f.get("native_support_rank") or 9999),
                                      -float(f["weight"]), str(f["type"]),
                                      int(f.get("resnum") or 0), str(f.get("atom") or "")))
        reserved: set[int] = set()
        for f in supported:
            if len(selected) >= n_keep:
                break
            resnum = int(f.get("resnum") or 0)
            key = (resnum, str(f["type"]), str(f.get("atom") or ""))
            if resnum in reserved or key in seen:
                continue
            selected.append(f)
            reserved.add(resnum)
            seen.add(key)
            if len(reserved) >= native_support_max_residues:
                break
    return selected


def load_hotspots(path: Path, top_hotspots: int = 25, pair_features: int = 24,
                  native_support_max_residues: int = 8) -> dict:
    """Build the Stage-1 hotspot set and the Stage-2 pair query."""
    data = json.loads(Path(path).read_text())
    features = data["features"]

    def to_query(f: dict) -> dict:
        return {
            "query_family": RECEPTOR_TO_QUERY_FAMILY[f["type"]],
            "weight": float(f["weight"]),
            "x": f["x"], "y": f["y"], "z": f["z"],
            "resnum": f.get("resnum"),
            "query_group": f.get("curated_group"),
            "priority_feature": _query_priority(f) != NON_CURATED_PRIORITY,
        }

    hotspots = [to_query(f) for f in
                select_query_features(features, top_hotspots, native_support_max_residues)]
    pair_query = [to_query(f) for f in
                  select_query_features(features, pair_features, native_support_max_residues)]
    return {"hotspots": hotspots, "pair_query": pair_query, "n_features": len(features)}


def evaluate(smiles: str, pharmacophore: dict,
             required_groups: tuple[str, ...] = ("ECD anchoring", "Upper TMD activation pocket")) -> dict:
    """Run Stages 0-2 on one molecule and return every published quantity."""
    s0 = stage0(smiles)
    if not s0["pass"]:
        return {"smiles": smiles, "stage0_pass": False, "stage0_reason": s0["reason"],
                "cascade_score_pct": 0.0}

    mol = s0["mol"]
    feats = typed_features(mol)
    if not feats:
        return {"smiles": smiles, "stage0_pass": True, "stage1_gate_pass": False,
                "stage0_reason": "no typed features", "cascade_score_pct": 0.0}

    s1 = stage1(feats, pharmacophore["hotspots"], required_groups)
    # Explicit fixed caps used by the headline and retrospective reruns.
    caps = dict(BASE_TYPE_CAPS)
    lig_pairs, lig_w = ligand_pair_catalogue(mol, feats, caps)
    q_pairs, q_w = receptor_pair_query(pharmacophore["pair_query"])
    s2 = stage2(lig_pairs, lig_w, q_pairs, q_w)

    cascade = 100.0 * (CASCADE_WEIGHTS["hotspot"] * s1["H"]
                       + CASCADE_WEIGHTS["pair_hash"] * s2["pair_hash_overlap"])
    return {
        "smiles": smiles,
        "stage0_pass": True,
        "stage0_reason": "",
        "n_typed_features": len(feats),
        **{k: v for k, v in s1.items()},
        "pair_hash_overlap_pct": 100.0 * s2["pair_hash_overlap"],
        "pair_hash_recall_pct": 100.0 * s2["pair_hash_recall"],
        "pair_hash_precision_pct": 100.0 * s2["pair_hash_precision"],
        "cascade_score_pct": cascade,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pharmacophore", default="maps/pharmacophore_rigorous.json")
    ap.add_argument("--smiles", nargs="*", help="SMILES to score")
    ap.add_argument("--examples", action="store_true",
                    help="score the bundled example set instead")
    args = ap.parse_args()

    root = Path(__file__).resolve().parent.parent
    ph = load_hotspots(root / args.pharmacophore)

    if args.examples:
        import csv

        rows = list(csv.DictReader(
            (Path(__file__).parent / "examples" / "example_molecules.csv").open()))
        smiles_list = [(r["ligand_id"], r["smiles"]) for r in rows]
    else:
        smiles_list = [(f"mol{i + 1}", s) for i, s in enumerate(args.smiles or [])]

    if not smiles_list:
        ap.error("pass --smiles or --examples")

    for name, smi in smiles_list:
        r = evaluate(smi, ph)
        print(f"\n{name}")
        if not r["stage0_pass"]:
            print(f"  Stage 0: FAIL ({r['stage0_reason']})")
            continue
        print(f"  Stage 0: pass")
        print(f"  Stage 1: H = {r['H']:.4f} ({r['hotspot_weighted_pct']:.2f}%), "
              f"exact residues {r['exact_residue_count']}, groups {r['group_count']}, "
              f"gate {'pass' if r['gate_pass'] else 'FAIL'}")
        print(f"  Stage 2: recall {r['pair_hash_recall_pct']:.2f}%, "
              f"precision {r['pair_hash_precision_pct']:.2f}%, "
              f"O_pair {r['pair_hash_overlap_pct']:.2f}%")
        print(f"  Cascade score: {r['cascade_score_pct']:.3f}")


if __name__ == "__main__":
    main()
