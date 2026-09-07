#!/usr/bin/env python3
"""Per-active attrition accounting for the production-constrained screen.

For every in-domain active in every system, this reports *which* stage removed
it and *why*, so a reader can judge whether a loss reflects a tunable threshold
(a structural alert, a property value just outside the envelope, a cascade
score just below the shortlist cut) or a genuine mismatch between the
topological gate and that target's chemistry.

Inputs are the released production artifacts
``evidence/outputs/production_<system>/stage012_evaluation.csv``, which record,
for every molecule embedded in the 30,000-molecule background, its label, the
gate it failed (or ``candidate`` if it reached the shortlist), and its cascade
score percentile.

The Stage-0 gate is reimplemented here from its published definition rather
than imported from the production engine, so this script is independently
runnable; ``--check`` verifies that the reimplementation reproduces the
engine's own per-molecule ``gate`` column.

Outputs (to ``evidence/outputs/attrition/``):
  * ``active_attrition_per_molecule.csv`` - one row per active
  * ``active_attrition_summary.csv``      - counts by system x stage x reason
  * LaTeX tables for the manuscript and appendix
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "evidence" / "outputs" / "attrition"

SYSTEMS = {
    "glp1r": "GLP-1R",
    "ghsr": "GHSR",
    "ntsr1": "NTSR1",
    "mdm2": "MDM2-p53",
}

# Stage-0 property gate exactly as implemented in the screening engine's
# ``passes_property_gate``. Note that rotatable-bond count is computed and
# reported by the engine but is NOT part of the gate.
PROPERTY_BOUNDS = {
    "MW": (None, 500.0),
    "LogP": (-1.0, 5.0),
    "HBD": (None, 5),
    "HBA": (None, 10),
}

# Reactive-group alerts applied alongside the PAINS A/B/C catalogs in the
# engine's strict chemistry gate.
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
for _cat in (
    FilterCatalogParams.FilterCatalogs.PAINS_A,
    FilterCatalogParams.FilterCatalogs.PAINS_B,
    FilterCatalogParams.FilterCatalogs.PAINS_C,
):
    _params.AddCatalog(_cat)
PAINS_CATALOG = FilterCatalog(_params)
REACTIVE_PATTERNS = {k: Chem.MolFromSmarts(v) for k, v in REACTIVE_SMARTS.items()}
_LARGEST_FRAGMENT = rdMolStandardize.LargestFragmentChooser()


def standardize(smiles: str) -> Chem.Mol | None:
    """Parse and standardize as the engine does: keep the largest fragment.

    Salt stripping matters here: several ChEMBL actives are supplied as salts
    whose parent cation only exceeds the LogP bound once the counterion is
    removed, so properties must be computed post-standardization.
    """
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None
    try:
        mol = _LARGEST_FRAGMENT.choose(mol)
    except Exception:
        return None
    return mol


def chemistry_alerts(mol: Chem.Mol) -> list[str]:
    """PAINS and reactive-group alerts, matching the engine's strict gate."""
    alerts: list[str] = []
    match = PAINS_CATALOG.GetFirstMatch(mol)
    if match is not None:
        alerts.append(f"PAINS:{match.GetDescription()}")
    for label, patt in REACTIVE_PATTERNS.items():
        if patt is not None and mol.HasSubstructMatch(patt):
            alerts.append(f"reactive:{label}")
    return alerts


def property_violations(mol: Chem.Mol) -> list[str]:
    """Stage-0 property-envelope violations, as ``name=value`` strings."""
    values = {
        "MW": Descriptors.MolWt(mol),
        "LogP": Crippen.MolLogP(mol),
        "HBD": Descriptors.NumHDonors(mol),
        "HBA": Descriptors.NumHAcceptors(mol),
    }
    out = []
    for name, (lo, hi) in PROPERTY_BOUNDS.items():
        v = values[name]
        if (lo is not None and v < lo) or (hi is not None and v > hi):
            out.append(f"{name}={v:.1f}" if name in ("MW", "LogP") else f"{name}={int(v)}")
    return out


def classify(row: pd.Series, shortlist_pct: float) -> dict[str, object]:
    """Assign one active a stage of loss and a human-readable reason."""
    mol = standardize(row["smiles"])
    if mol is None:
        return {"stage_lost": "Stage 0", "reason_class": "unparseable", "reason_detail": "invalid SMILES"}

    props = property_violations(mol)
    alerts = chemistry_alerts(mol)
    gate = str(row["gate"])

    if gate == "stage0_fail":
        # Property envelope is checked before the structural-alert screen, so
        # report a property violation in preference to an alert when both hold.
        if props:
            return {
                "stage_lost": "Stage 0",
                "reason_class": "property envelope",
                "reason_detail": ", ".join(props),
            }
        if alerts:
            return {
                "stage_lost": "Stage 0",
                "reason_class": "structural alert",
                "reason_detail": "; ".join(alerts[:2]),
            }
        return {"stage_lost": "Stage 0", "reason_class": "standardization", "reason_detail": "failed standardization"}

    if gate == "stage1_fail":
        return {
            "stage_lost": "Stage 1",
            "reason_class": "hotspot gate",
            "reason_detail": "too few exact priority-residue matches or curated groups",
        }

    # Reached the candidate pool: lost (if at all) at the shortlist cut.
    pct = float(row["cascade_score_pct"])
    if row.get("_in_shortlist", False):
        return {"stage_lost": "-", "reason_class": "retained", "reason_detail": "reached native scoring"}
    return {
        "stage_lost": "Stage-3 shortlist",
        "reason_class": "below shortlist cut",
        "reason_detail": f"cascade score {pct:.1f}, below the top-{shortlist_pct:.0f}% cut",
    }


def analyze_system(key: str, shortlist_fraction: float = 0.05) -> pd.DataFrame:
    path = ROOT / "evidence" / "outputs" / f"production_{key}" / "stage012_evaluation.csv"
    df = pd.read_csv(path)

    # Reproduce the shortlist: the top ``shortlist_fraction`` of Stage-1-pass
    # candidates by cascade score, which is what the production harness applies.
    cand = df[df.gate == "candidate"].sort_values("cascade_score_pct", ascending=False)
    k = int(len(cand) * shortlist_fraction + 0.999999)
    shortlisted = set(cand.head(k)["ligand_id"])

    actives = df[df.label == "active"].copy()
    actives["_in_shortlist"] = actives["ligand_id"].isin(shortlisted)

    rows = []
    for _, r in actives.iterrows():
        info = classify(r, shortlist_fraction * 100)
        rows.append(
            {
                "system": SYSTEMS[key],
                "ligand_id": r["ligand_id"],
                "smiles": r["smiles"],
                "cascade_score_pct": r["cascade_score_pct"],
                "engine_gate": r["gate"],
                **info,
            }
        )
    return pd.DataFrame(rows)


def check_against_engine(per_mol: pd.DataFrame) -> None:
    """Confirm the reimplemented Stage-0 gate agrees with the engine."""
    s0 = per_mol[per_mol.engine_gate == "stage0_fail"]
    explained = s0[s0.reason_class.isin({"property envelope", "structural alert", "unparseable"})]
    print(f"\n[check] engine Stage-0 active failures: {len(s0)}")
    print(f"[check] reproduced by the published gate definition: {len(explained)}")
    if len(s0):
        print(f"[check] agreement: {100 * len(explained) / len(s0):.1f}%")
    unexplained = s0[~s0.index.isin(explained.index)]
    for _, r in unexplained.iterrows():
        print(f"[check]   unexplained: {r.system} {r.ligand_id}")


def latex_summary(summary: pd.DataFrame) -> str:
    """Main-text table: where each system's actives are lost, and to what."""
    order = ["Stage 0", "Stage 1", "Stage-3 shortlist", "-"]
    lines = [
        r"\begin{table}[!t]",
        r"\caption{Stage and cause of in-domain active loss under production constraints.",
        r"Each system's actives were embedded in a 30\,000-molecule ZINC background and screened",
        r"with the million-compound-run settings. Entries give the number of actives removed at",
        r"each stage, split by cause. Structural-alert losses are strict-mode PAINS and",
        r"reactive-group rejections, which are a tunable gate setting rather than a failure of",
        r"topological matching.}",
        r"\label{tab:attrition}",
        r"\centering",
        r"\begin{tabular}{llrrrr}",
        r"\hline",
        r"Stage & Cause & GLP-1R & GHSR & NTSR1 & MDM2 \\",
        r"\hline",
    ]
    for stage in order[:-1]:
        sub = summary[summary.stage_lost == stage]
        for cause in sorted(sub.reason_class.unique()):
            row = sub[sub.reason_class == cause]
            cells = []
            for sysname in ["GLP-1R", "GHSR", "NTSR1", "MDM2-p53"]:
                n = int(row[row.system == sysname]["n"].sum())
                cells.append(str(n) if n else "--")
            lines.append(f"{stage} & {cause} & " + " & ".join(cells) + r" \\")
    lines.append(r"\hline")
    ret = summary[summary.stage_lost == "-"]
    cells = []
    for sysname in ["GLP-1R", "GHSR", "NTSR1", "MDM2-p53"]:
        n = int(ret[ret.system == sysname]["n"].sum())
        cells.append(str(n))
    lines.append(r"\textbf{Reached native scoring} & & " + " & ".join(cells) + r" \\")
    lines += [r"\hline", r"\end{tabular}", r"\end{table}", ""]
    return "\n".join(lines)


def tex_escape(text: str) -> str:
    """Escape the LaTeX specials that appear in generated cause strings.

    Unescaped ``%`` is the dangerous one: it comments out the rest of the line,
    including the row terminator, and silently merges table rows.
    """
    out = str(text)
    for char, repl in (("\\", r"\textbackslash{}"), ("&", r"\&"), ("%", r"\%"),
                       ("$", r"\$"), ("#", r"\#"), ("_", r"\_"),
                       ("{", r"\{"), ("}", r"\}")):
        out = out.replace(char, repl)
    return out


def latex_per_molecule(per_mol: pd.DataFrame) -> str:
    """Appendix long table: every lost active, with its cause."""
    lost = per_mol[per_mol.stage_lost != "-"].copy()
    # The cause column wraps rather than truncating: a reviewer asked why each
    # active was lost, and a clipped reason answers only half the question.
    header = r"System & Active & Stage lost & Cause \\"
    lines = [
        r"{\footnotesize",
        r"\begin{longtable}{@{}llp{0.16\textwidth}p{0.38\textwidth}@{}}",
        r"\caption{Every in-domain active removed before native scoring, with the stage that",
        r"removed it and the specific cause.}",
        r"\label{tab:attrition_full}\\",
        r"\hline",
        header,
        r"\hline",
        r"\endfirsthead",
        r"\hline",
        header,
        r"\hline",
        r"\endhead",
    ]
    for sysname in ["GLP-1R", "GHSR", "NTSR1", "MDM2-p53"]:
        for _, r in lost[lost.system == sysname].iterrows():
            lines.append(
                f"{sysname} & \\texttt{{{r.ligand_id}}} & {r.stage_lost} & "
                f"{tex_escape(str(r.reason_detail))} \\\\"
            )
    lines += [r"\hline", r"\end{longtable}", r"}", ""]
    return "\n".join(lines)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--shortlist-fraction", type=float, default=0.05)
    ap.add_argument("--check", action="store_true", help="verify the gate reimplementation")
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    per_mol = pd.concat(
        [analyze_system(k, args.shortlist_fraction) for k in SYSTEMS], ignore_index=True
    )
    per_mol.to_csv(OUT / "active_attrition_per_molecule.csv", index=False)

    summary = (
        per_mol.groupby(["system", "stage_lost", "reason_class"]).size().reset_index(name="n")
    )
    summary.to_csv(OUT / "active_attrition_summary.csv", index=False)

    (OUT / "attrition_table.tex").write_text(latex_summary(summary), encoding="utf-8")
    (OUT / "attrition_full_table.tex").write_text(latex_per_molecule(per_mol), encoding="utf-8")

    print("Per-system active attrition (actives removed / total):")
    for sysname in ["GLP-1R", "GHSR", "NTSR1", "MDM2-p53"]:
        sub = per_mol[per_mol.system == sysname]
        kept = int((sub.stage_lost == "-").sum())
        print(f"  {sysname:10s} retained {kept}/{len(sub)}")
        for stage in ["Stage 0", "Stage 1", "Stage-3 shortlist"]:
            ss = sub[sub.stage_lost == stage]
            if len(ss):
                causes = ss.reason_class.value_counts().to_dict()
                print(f"      lost at {stage:18s} {len(ss):3d}  {causes}")

    if args.check:
        check_against_engine(per_mol)
    print(f"\nWrote outputs to {OUT}")


if __name__ == "__main__":
    # Current entry point uses measured final survival under non-excluding
    # alerts. The functions above are retained for historical source audits.
    from analyze_alert_disabled_attrition import main as current_main
    current_main()
