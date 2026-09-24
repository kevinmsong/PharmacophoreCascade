#!/usr/bin/env python3
"""Re-derive every figure the letter's Point 4 cites, from the run outputs.

Point 4 gained a lot of specific numbers when the superseded alert policy was
measured. Each is checked here against the per-molecule survival records and
the Stage 0-2 evaluations rather than against the prose that quotes them.
"""
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ANNOTATE = ROOT / "evidence/outputs/absolute_floor"
STRICT = ROOT / "evidence/outputs/strict_alert_check"


def survival(root, system, policy, stage):
    rows = csv.DictReader(open(root / system / policy / "per_molecule_survival.csv"))
    return sum(1 for r in rows if r["label"] == "active" and r[stage] == "True")


def final_ids(root, system, policy="percentage"):
    rows = csv.DictReader(open(root / system / policy / "per_molecule_survival.csv"))
    return {r["ligand_id"] for r in rows
            if r["label"] == "active" and r["final_ranked"] == "True"}


def alert_excluded(system):
    rows = csv.DictReader(open(STRICT / system / "stage012_evaluation.csv"))
    return sum(1 for r in rows
               if r["label"] == "active" and r["status"] == "chemistry_filtered")


def candidates(root, system):
    summary = json.load(open(root / system / "policy_summary.json"))
    entry = next(p for p in summary if p["policy"] == "percentage")
    return entry["n_candidates"], entry["n_shortlist"]


def main() -> None:
    checks = []

    def expect(label, got, want):
        checks.append((label, got, want, got == want))

    # NTSR1: the headline confirmation.
    expect("NTSR1 strict final actives (5% only)",
           survival(STRICT, "ntsr1", "percentage", "final_ranked"), 15)
    expect("NTSR1 annotate final actives (5% only)",
           survival(ANNOTATE, "ntsr1", "percentage", "final_ranked"), 18)
    expect("NTSR1 actives alert-excluded at Stage 0", alert_excluded("ntsr1"), 17)
    expect("NTSR1 strict Stage-0 survivors",
           survival(STRICT, "ntsr1", "percentage", "stage0_pass"), 31)
    expect("NTSR1 annotate Stage-0 survivors",
           survival(ANNOTATE, "ntsr1", "percentage", "stage0_pass"), 48)
    expect("NTSR1 strict Stage-0-2 survivors",
           survival(STRICT, "ntsr1", "percentage", "stage012_pass"), 23)
    expect("NTSR1 annotate Stage-0-2 survivors",
           survival(ANNOTATE, "ntsr1", "percentage", "stage012_pass"), 35)
    n_s, k_s = candidates(STRICT, "ntsr1")
    n_a, k_a = candidates(ANNOTATE, "ntsr1")
    expect("NTSR1 strict candidate pool", n_s, 384)
    expect("NTSR1 annotate candidate pool", n_a, 398)
    expect("NTSR1 shortlist size, both policies", (k_s, k_a), (20, 20))
    expect("NTSR1 strict retention under the 1,000 minimum",
           survival(STRICT, "ntsr1", "floor1000", "final_ranked"), 23)

    # The 18 -> 12 -> 15 decomposition.
    warn, strict = final_ids(ANNOTATE, "ntsr1"), final_ids(STRICT, "ntsr1")
    flagged = {r["ligand_id"] for r in
               csv.DictReader(open(STRICT / "ntsr1" / "stage012_evaluation.csv"))
               if r["status"] == "chemistry_filtered"}
    expect("NTSR1 kept under both", len(warn & strict), 12)
    expect("NTSR1 lost under strict", len(warn - strict), 6)
    expect("NTSR1 lost that are alert-flagged", len((warn - strict) & flagged), 6)
    expect("NTSR1 promoted under strict", len(strict - warn), 3)
    expect("NTSR1 promoted that are alert-flagged", len((strict - warn) & flagged), 0)

    # GHSR: the control.
    expect("GHSR strict final actives (5% only)",
           survival(STRICT, "ghsr", "percentage", "final_ranked"), 9)
    expect("GHSR annotate final actives (5% only)",
           survival(ANNOTATE, "ghsr", "percentage", "final_ranked"), 9)
    expect("GHSR final actives are the same molecules",
           final_ids(STRICT, "ghsr") == final_ids(ANNOTATE, "ghsr"), True)
    expect("GHSR actives alert-excluded at Stage 0", alert_excluded("ghsr"), 9)
    expect("GHSR strict Stage-0 survivors",
           survival(STRICT, "ghsr", "percentage", "stage0_pass"), 41)
    expect("GHSR annotate Stage-0 survivors",
           survival(ANNOTATE, "ghsr", "percentage", "stage0_pass"), 50)
    expect("GHSR strict Stage-0-2 survivors",
           survival(STRICT, "ghsr", "percentage", "stage012_pass"), 30)
    expect("GHSR annotate Stage-0-2 survivors",
           survival(ANNOTATE, "ghsr", "percentage", "stage012_pass"), 39)
    expect("GHSR shortlist size", candidates(STRICT, "ghsr")[1], 14)

    bad = 0
    for label, got, want, ok in checks:
        if not ok:
            bad += 1
            print(f"  MISMATCH {label}: got {got!r}, letter says {want!r}")
    print(f"{len(checks) - bad}/{len(checks)} Point 4 figures re-derived and consistent")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
