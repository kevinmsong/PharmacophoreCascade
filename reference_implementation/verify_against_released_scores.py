#!/usr/bin/env python3
"""Check the reference implementation against the screen's released scores.

This answers the question the reference implementation exists for: do the
equations printed in the manuscript, implemented independently from that
description alone, reproduce what the production screen actually computed?

The comparison uses the released per-molecule GLP-1R benchmark table, which
records each molecule's status under the full cascade. Stage-0 and Stage-1
outcomes are directly comparable; the continuous Stage-2 quantities are
compared against the headline million-compound run's audit table, which
publishes ``hotspot_weighted_pct``, ``pair_hash_overlap_pct``, and
``cascade_score_pct`` per molecule.

Whatever agreement is obtained is reported as-is, including any disagreement.

    python reference_implementation/verify_against_released_scores.py
    python reference_implementation/verify_against_released_scores.py --n 500
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from stage_scoring_reference import evaluate, load_hotspots  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
BENCH = ROOT / "evidence" / "data" / "machine_readable" / "glp1r_benchmark_scored.csv"
AUDIT = ROOT / "results" / "screening_full_1M_topological_hashed.csv"
PHARM = ROOT / "maps" / "pharmacophore_rigorous.json"


def check_stage_status(ph: dict) -> None:
    """Stage-0 / Stage-1 pass-fail agreement on the released benchmark."""
    df = pd.read_csv(BENCH)
    status = df["full_cascade_status"].astype(str)

    rows = []
    for _, r in df.iterrows():
        res = evaluate(str(r["canonical_smiles"]), ph)
        engine = str(r["full_cascade_status"])
        # The engine records "property_filtered"/"chemistry_filtered" for
        # Stage 0 and "hotspot_filtered" for the Stage-1 gate.
        ref_stage0 = res["stage0_pass"]
        ref_gate = bool(res.get("gate_pass", False))
        rows.append({
            "ligand_id": r["ligand_id"],
            "engine_status": engine,
            "ref_stage0_pass": ref_stage0,
            "ref_stage1_gate_pass": ref_gate,
        })
    out = pd.DataFrame(rows)

    engine_stage0_fail = status.isin({"property_filtered", "chemistry_filtered", "invalid"})
    engine_stage1_fail = status.eq("hotspot_filtered")

    s0_agree = (out["ref_stage0_pass"] != engine_stage0_fail.to_numpy()).mean()
    reached = ~engine_stage0_fail.to_numpy()
    s1_agree = (
        out.loc[reached, "ref_stage1_gate_pass"].to_numpy()
        != engine_stage1_fail.to_numpy()[reached]
    ).mean()

    print("=== Stage-0 / Stage-1 status agreement (GLP-1R benchmark, n = "
          f"{len(out)}) ===")
    print(f"  Stage-0 pass/fail agreement : {100 * s0_agree:.1f}%")
    print(f"  Stage-1 gate agreement      : {100 * s1_agree:.1f}% "
          f"(over {int(reached.sum())} Stage-0 survivors)")
    print("  engine status counts        :",
          dict(status.value_counts()))
    return out


def check_continuous(ph: dict, n: int, seed: int) -> None:
    """Numerical agreement on H, O_pair, and the cascade score."""
    if not AUDIT.exists():
        print(f"\n[skip] {AUDIT.name} not present; continuous check needs the "
              "headline run audit table.")
        return

    audit = pd.read_csv(
        AUDIT,
        usecols=["zinc_id", "smiles", "hotspot_weighted_pct",
                 "pair_hash_overlap_pct", "cascade_score_pct"],
    )
    sample = audit.sample(n=min(n, len(audit)), random_state=seed)

    recs = []
    for _, r in sample.iterrows():
        res = evaluate(str(r["smiles"]), ph)
        if not res["stage0_pass"]:
            continue
        recs.append({
            "zinc_id": r["zinc_id"],
            "engine_H": r["hotspot_weighted_pct"],
            "ref_H": res["hotspot_weighted_pct"],
            "engine_O": r["pair_hash_overlap_pct"],
            "ref_O": res["pair_hash_overlap_pct"],
            "engine_cascade": r["cascade_score_pct"],
            "ref_cascade": res["cascade_score_pct"],
        })
    cmp = pd.DataFrame(recs)
    if cmp.empty:
        print("\n[warn] no molecules survived Stage 0 in the sample.")
        return

    print(f"\n=== Continuous-score agreement (headline run, n = {len(cmp)}) ===")
    for label, a, b in (("hotspot fraction H (%)", "engine_H", "ref_H"),
                        ("pair overlap O_pair (%)", "engine_O", "ref_O"),
                        ("cascade score (%)", "engine_cascade", "ref_cascade")):
        d = cmp[b] - cmp[a]
        r = np.corrcoef(cmp[a], cmp[b])[0, 1]
        print(f"  {label:26s} Pearson r = {r:.4f}  "
              f"mean abs diff = {np.abs(d).mean():.3f}  max = {np.abs(d).max():.3f}")

    # The cascade score must equal 0.4 H + 0.6 O_pair on the engine's own
    # columns; this confirms eq. (3) as printed.
    recomputed = 0.4 * cmp["engine_H"] + 0.6 * cmp["engine_O"]
    print(f"\n  eq. (3) check on the engine's own columns: "
          f"max |0.4H + 0.6*O_pair - cascade| = "
          f"{np.abs(recomputed - cmp['engine_cascade']).max():.6f}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, default=300, help="molecules for the continuous check")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    ph = load_hotspots(PHARM)
    print(f"Loaded receptor pharmacophore: {ph['n_features']} features, "
          f"{len(ph['hotspots'])} Stage-1 hotspots, "
          f"{len(ph['pair_query'])} Stage-2 pair-query features.\n")
    check_stage_status(ph)
    check_continuous(ph, args.n, args.seed)


if __name__ == "__main__":
    main()
