#!/usr/bin/env python3
"""Recompute the decoy-replicate benchmarks on the corrected replicate libraries.

The replicate libraries shipped with the first revision still held 10 actives per
system. They were built before GHSR, NTSR1, and MDM2--p53 were expanded to 50
actives and were never regenerated, so the replicate analysis and the primary
benchmarks described different active sets. ``build_decoy_replicates.py`` has
since rebuilt them from the current libraries; this script reruns the benchmarks
on top of them.

Settings mirror ``run_alert_disabled_suite.benchmark(system, workers, seed=n)``
exactly: structural alerts recorded rather than excluding, the 0.25 hotspot
cascade weight, fixed feature caps, no native-only arm, and no bootstrap inside
the per-replicate runs, because the replicate spread is what the analysis
reports.

The first-revision replicate outputs, their execution manifests, their configs,
and their identifier-ordered backups are moved aside before anything is written,
so the previous round stays reproducible.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
META = ROOT / "evidence/outputs/alert_disabled_revision"
OUT_BASE = ROOT / "evidence/outputs/decoy_replicates"
ARCHIVE = ROOT / "evidence/outputs/decoy_replicates_10active_archive"

# (base benchmark config, native rerank config, required hotspot groups)
SPECS = {
    "glp1r": ("benchmark_glp1r_full.yaml", None, None),
    "ghsr": ("benchmark_ghsr_full.yaml", "study_ghsr_benchmark.yaml", ""),
    "ntsr1": ("benchmark_ntsr1_full.yaml", "study_ntsr1_benchmark.yaml", ""),
    "mdm2": ("benchmark_mdm2_full.yaml", "study_mdm2_benchmark.yaml", ""),
}
SEEDS = [1, 2, 3, 4, 5]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def archive_prior(system: str, seed: int) -> None:
    """Move the first-revision artifacts for one replicate job out of the way."""
    job = f"{system}_seed{seed}"
    dest = ARCHIVE / job
    dest.mkdir(parents=True, exist_ok=True)

    run_dir = OUT_BASE / job
    if run_dir.exists():
        target = dest / "outputs"
        if not target.exists():
            run_dir.rename(target)

    for name in (f"{job}.json", f"{job}.log"):
        path = META / name
        if path.exists() and not (dest / name).exists():
            path.rename(dest / name)

    cfg = META / "configs" / f"{job}.yaml"
    if cfg.exists() and not (dest / f"{job}.yaml").exists():
        cfg.rename(dest / f"{job}.yaml")

    # The tie-aware step reads this as the identifier-ordered baseline; a stale
    # copy from the previous round would silently become the comparison point.
    backup = META / "identifier_order_evaluation" / job
    if backup.exists():
        target = dest / "identifier_order_evaluation"
        if not target.exists():
            backup.rename(target)


def build_config(system: str, seed: int, workers: int) -> Path:
    base, _native, _groups = SPECS[system]
    cfg = yaml.safe_load((ROOT / "evidence/configs" / base).read_text())
    job = f"{system}_seed{seed}"
    out = OUT_BASE / job

    settings = cfg["benchmark"]
    settings["library_csv"] = f"evidence/data/replicates/{system}_decoyset{seed}.csv"
    settings["include_native_only"] = False
    settings["bootstrap_iterations"] = 0
    settings.update(
        chemistry_gate_mode="warn_only",
        workers=workers,
        native_workers=workers,
        cascade_hotspot_weight=0.25,
        feature_cap_mode="fixed",
        output_dir=str((out / "benchmark_external").relative_to(ROOT)),
    )

    folder = META / "configs"
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / f"{job}.yaml"
    path.write_text(yaml.safe_dump(cfg, sort_keys=False))
    return path


#: Written into every manifest this script produces. A manifest without it
#: belongs to the first revision, whose replicate libraries held 10 actives,
#: and must be archived rather than treated as a completed run to resume from.
ROUND_MARKER = "decoy replicates rebuilt on the current active sets"


def run(system: str, seed: int, workers: int) -> None:
    job = f"{system}_seed{seed}"
    manifest = META / f"{job}.json"
    if manifest.exists():
        record = json.loads(manifest.read_text())
        if record.get("round") == ROUND_MARKER:
            if record.get("status") == "completed":
                print(f"Already completed: {job}", flush=True)
                return
            raise RuntimeError(f"Inspect unfinished job before restart: {manifest}")
        print(f"Archiving first-revision run: {job}", flush=True)

    archive_prior(system, seed)
    cfg_path = build_config(system, seed, workers)
    settings = yaml.safe_load(cfg_path.read_text())["benchmark"]
    out = OUT_BASE / job

    base, native, groups = SPECS[system]
    native_path = (ROOT / "GLP1_top_ligand_analysis/configs"
                   / (native or "study_top1000_full_1M_native_tuned_no_docking.yaml"))

    command = [sys.executable, "-u", "evidence/run_target_benchmark.py",
               "--benchmark-config", str(cfg_path), "--output-dir", str(out),
               "--native-rerank-config", str(native_path)]
    if groups is not None:
        command += ["--required-groups", groups]

    inputs = [cfg_path, ROOT / settings["library_csv"],
              ROOT / settings["pharmacophore_path"], native_path]
    record = {
        "job": job,
        "round": ROUND_MARKER,
        "status": "running",
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "command": command,
        "settings": {
            "chemistry_gate_mode": "warn_only",
            "native_pains_filter": False,
            "workers": workers,
            "cascade_hotspot_weight": 0.25,
            "feature_caps": "fixed",
            "required_groups": groups,
            "variant": None,
            "decoy_seed": seed,
            "replicate_actives": "current benchmark active set",
        },
        "cached_scores_used": False,
        "input_sha256": {str(p.relative_to(ROOT)): sha(p) for p in inputs},
    }
    META.mkdir(parents=True, exist_ok=True)
    manifest.write_text(json.dumps(record, indent=2))

    started = time.perf_counter()
    print(f"Starting {job}", flush=True)
    with open(META / f"{job}.log", "w", encoding="utf-8") as log:
        result = subprocess.run(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
    record.update(
        status="completed" if result.returncode == 0 else "failed",
        returncode=result.returncode,
        elapsed_seconds=time.perf_counter() - started,
        completed_utc=datetime.now(timezone.utc).isoformat(),
    )
    manifest.write_text(json.dumps(record, indent=2))
    if result.returncode:
        raise RuntimeError(f"{job} failed; see {META / (job + '.log')}")
    print(f"Completed {job} in {record['elapsed_seconds'] / 60:.1f} min", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--systems", nargs="+", default=list(SPECS), choices=list(SPECS))
    parser.add_argument("--seeds", nargs="+", type=int, default=SEEDS)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()

    for system in args.systems:
        for seed in args.seeds:
            run(system, seed, args.workers)


if __name__ == "__main__":
    main()
