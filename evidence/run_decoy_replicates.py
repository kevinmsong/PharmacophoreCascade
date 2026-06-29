#!/usr/bin/env python3
"""
Benchmark every decoy replicate and report metric mean +/- SD per system.

For each system x replicate-seed, runs the external benchmark (point metrics
only; native-only and bootstrap disabled for speed) as a separate process so
the per-system required-group gate and native-reference config are applied
correctly, then aggregates ROC-AUC, PR-AUC, BEDROC, EF1%, and top-10 recovery
across replicates -> evidence/outputs/revision/decoy_replicate_robustness.csv.
"""
from __future__ import annotations

import argparse
import copy
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd
import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

# per-system: base benchmark config, native-reference study config, required-groups gate
SYSTEMS = {
    "glp1r": {
        "base": HERE / "configs" / "benchmark.yaml",
        "native_cfg": None,          # default GLP-1R study config
        "required_groups": None,     # keep default ECD/TMD gate
    },
    "ghsr": {
        "base": HERE / "configs" / "benchmark_ghsr.yaml",
        "native_cfg": "GLP1_top_ligand_analysis/configs/study_ghsr_benchmark.yaml",
        "required_groups": "",
    },
    "ntsr1": {
        "base": HERE / "configs" / "benchmark_ntsr1.yaml",
        "native_cfg": "GLP1_top_ligand_analysis/configs/study_ntsr1_benchmark.yaml",
        "required_groups": "",
    },
    "mdm2": {
        "base": HERE / "configs" / "benchmark_mdm2.yaml",
        "native_cfg": "GLP1_top_ligand_analysis/configs/study_mdm2_benchmark.yaml",
        "required_groups": "",
    },
}
METRICS = ["roc_auc", "pr_auc", "bedroc", "ef_1pct", "top10_recovery"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--systems", nargs="+", default=list(SYSTEMS), choices=list(SYSTEMS))
    ap.add_argument("--seeds", nargs="+", type=int, default=[1, 2, 3, 4, 5])
    ap.add_argument("--replicate-dir", default=str(HERE / "data" / "replicates"))
    args = ap.parse_args()
    rep_dir = Path(args.replicate_dir)
    out_base = ROOT / "evidence" / "outputs" / "decoy_replicates"

    records = []
    for system in args.systems:
        spec = SYSTEMS[system]
        base_cfg = yaml.safe_load(open(spec["base"]))
        for seed in args.seeds:
            lib = rep_dir / f"{system}_decoyset{seed}.csv"
            if not lib.exists():
                print(f"[skip] missing {lib}")
                continue
            run_out = out_base / f"{system}_seed{seed}"
            cfg = copy.deepcopy(base_cfg)
            cfg["benchmark"]["library_csv"] = f"evidence/data/replicates/{system}_decoyset{seed}.csv"
            cfg["benchmark"]["output_dir"] = f"evidence/outputs/decoy_replicates/{system}_seed{seed}/benchmark_external"
            cfg["benchmark"]["include_native_only"] = False
            cfg["benchmark"]["bootstrap_iterations"] = 0
            with tempfile.NamedTemporaryFile("w", suffix=".yaml", delete=False) as tf:
                yaml.safe_dump(cfg, tf, sort_keys=False)
                cfg_path = tf.name
            cmd = [sys.executable, str(HERE / "run_target_benchmark.py"),
                   "--benchmark-config", cfg_path, "--output-dir", str(run_out)]
            if spec["native_cfg"]:
                cmd += ["--native-rerank-config", spec["native_cfg"]]
            if spec["required_groups"] is not None:
                cmd += ["--required-groups", spec["required_groups"]]
            print(f"  running {system} seed {seed} ...", flush=True)
            subprocess.run(cmd, check=True, cwd=str(ROOT), capture_output=True, text=True)
            summ = pd.read_csv(run_out / "benchmark_summary.csv").set_index("method")
            for method in summ.index:
                rec = {"system": system, "seed": seed, "method": method}
                for m in METRICS:
                    rec[m] = float(summ.loc[method, m]) if m in summ.columns else float("nan")
                records.append(rec)

    raw = pd.DataFrame(records)
    out_dir = ROOT / "evidence" / "outputs" / "revision"
    out_dir.mkdir(parents=True, exist_ok=True)
    raw.to_csv(out_dir / "decoy_replicate_raw.csv", index=False)
    agg = (
        raw.groupby(["system", "method"])[METRICS]
        .agg(["mean", "std"])
        .round(4)
    )
    agg.columns = [f"{m}_{stat}" for m, stat in agg.columns]
    agg = agg.reset_index()
    agg.to_csv(out_dir / "decoy_replicate_robustness.csv", index=False)
    print("\n=== decoy-replicate robustness (full_cascade) ===")
    fc = agg[agg["method"] == "full_cascade"]
    print(fc.to_string(index=False))


if __name__ == "__main__":
    main()
