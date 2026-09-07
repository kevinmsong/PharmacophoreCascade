#!/usr/bin/env python3
"""
Run the active-vs-decoy cascade benchmark for an arbitrary target or ablation.

This is a thin, reusable driver around ``src.benchmark.run_benchmark`` for the
external-benchmark path. It supports:

* the GLP-1R no-manual-curation ablation (automated pharmacophore + empty
  required-group gate), and
* additional peptide-receptor case studies (new library + automated interface
  pharmacophore + per-target native-reference complex).

It writes ``benchmark_summary.csv`` (point metrics + bootstrap CIs),
``benchmark_deltas.csv`` (paired deltas vs the full cascade) and
``benchmark_pairwise.csv`` to ``--output-dir``.

Usage
-----
    python evidence/run_target_benchmark.py \
        --benchmark-config evidence/configs/benchmark_glp1r_automated.yaml \
        --required-groups "" \
        --output-dir evidence/outputs/benchmark_glp1r_automated
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--benchmark-config", required=True, help="Benchmark YAML (external-benchmark schema).")
    ap.add_argument(
        "--native-rerank-config",
        default=None,
        help="Study YAML (path relative to repo root) providing the native-branch "
        "reference complex/chains. Defaults to the GLP-1R study config.",
    )
    ap.add_argument(
        "--required-groups",
        default=None,
        help='Value for CASCADIA_HOTSPOT_REQUIRED_GROUPS (";"-separated). Pass "" to '
        "disable the required-group gate for automated/generic pharmacophores.",
    )
    ap.add_argument("--output-dir", default=None, help="Where to write summary CSVs.")
    args = ap.parse_args()

    # Must be set BEFORE the screening engine is imported (read at module load).
    if args.required_groups is not None:
        os.environ["CASCADIA_HOTSPOT_REQUIRED_GROUPS"] = args.required_groups

    import run_optimized_1M_topological_hashed_screening as screening
    if args.required_groups is not None:
        screening.HOTSPOT_REQUIRED_GROUPS=tuple(x for x in args.required_groups.split(';') if x)
    screening.TYPE_CAPS=screening._BASE_TYPE_CAPS.copy()

    if args.native_rerank_config:
        screening.DEFAULT_NATIVE_RERANK_CONFIG = Path(args.native_rerank_config)

    from src.benchmark import run_benchmark

    cfg_path = Path(args.benchmark_config)
    with cfg_path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle) or {}
    requested_workers=int(config.get('benchmark',{}).get('workers',1))
    effective_workers=requested_workers
    allocation_reason='requested configuration'
    # Only the active revision queue opts into this runtime allocation profile.
    # Scoring parameters, per-ligand seeds, and original config files stay fixed.
    profile=ROOT/'evidence/outputs/alert_disabled_revision/runtime_allocation.json'
    revision_configs=(profile.parent/'configs').resolve()
    if profile.exists() and cfg_path.resolve().is_relative_to(revision_configs):
        allocation=json.loads(profile.read_text())
        headline=ROOT/'results/absolute_floor_1000/rerun_provenance.json'
        if json.loads(headline.read_text()).get('status')=='completed':
            dock=ROOT/'evidence/outputs/docking_top10/rerun_manifest.json'
            promoted=profile.parent/'headline_promotion.json'
            dock_done=(dock.exists() and promoted.exists() and
                json.loads(dock.read_text()).get('ligand_ids')==json.loads(promoted.read_text()).get('new_top10'))
            effective_workers=min(int(allocation['after_docking_workers'] if dock_done else allocation['after_headline_workers']),os.cpu_count() or 1)
            allocation_reason='headline and docking completed' if dock_done else 'headline completed; reserve cores for docking'
            config['benchmark'].update(workers=effective_workers,native_workers=effective_workers)
        else:
            production=[profile.parent/f'production_{system}.json' for system in ['glp1r','ghsr','ntsr1','mdm2']]
            if all(p.exists() and json.loads(p.read_text()).get('status')=='completed' for p in production):
                effective_workers=min(int(allocation.get('after_production_workers',requested_workers)),os.cpu_count() or 1)
                allocation_reason='production lane completed; headline still running'
                config['benchmark'].update(workers=effective_workers,native_workers=effective_workers)
    actual_out=Path(args.output_dir) if args.output_dir else ROOT/config['benchmark']['output_dir']
    actual_out.mkdir(parents=True,exist_ok=True)
    execution={'status':'running','requested_workers':requested_workers,'effective_workers':effective_workers,
        'allocation_reason':allocation_reason,'input_config_sha256':hashlib.sha256(cfg_path.read_bytes()).hexdigest(),
        'driver_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'effective_configuration':config,'chemistry_gate_mode':config['benchmark']['chemistry_gate_mode']}
    execution_path=actual_out/'benchmark_execution.json'
    execution_path.write_text(json.dumps(execution,indent=2))
    # Future revision jobs may use cores released between native invocations.
    # Only the worker count changes; the validated per-ligand scoring functions
    # and molecular parameters stay unchanged. Record each actual allocation.
    if profile.exists() and cfg_path.resolve().is_relative_to(revision_configs):
        from evidence import parallel_native
        original_install=parallel_native.install
        execution['native_worker_events']=[]
        def install_for_current_phase(workers):
            chosen=workers
            reason='job allocation; headline still running'
            current=json.loads(profile.read_text())
            head=ROOT/'results/absolute_floor_1000/rerun_provenance.json'
            if json.loads(head.read_text()).get('status')=='completed':
                dock=ROOT/'evidence/outputs/docking_top10/rerun_manifest.json'
                promotion=profile.parent/'headline_promotion.json'
                finished=(dock.exists() and promotion.exists() and
                    json.loads(dock.read_text()).get('ligand_ids')==json.loads(promotion.read_text()).get('new_top10'))
                chosen=min(int(current['after_docking_workers'] if finished else current['after_headline_workers']),os.cpu_count() or 1)
                reason='headline and docking completed' if finished else 'headline completed; reserve cores for docking'
            execution['native_worker_events'].append(dict(invocation=len(execution['native_worker_events'])+1,
                requested_workers=workers,actual_workers=chosen,reason=reason,
                started_utc=datetime.now(timezone.utc).isoformat()))
            execution_path.write_text(json.dumps(execution,indent=2))
            print(f'[run_target_benchmark] native invocation uses {chosen} workers: {reason}',flush=True)
            original_install(chosen)
        parallel_native.install=install_for_current_phase
    started=time.perf_counter()
    # Preserve actual archived benchmark coefficients, explicitly recorded in
    # revised configs. The old top-level cascade_weights described a separate
    # sensitivity sweep and was not consumed by this runner.
    weight=float(config.get('benchmark',{}).get('cascade_hotspot_weight',.25))
    screening.CASCADE_STAGE_WEIGHTS={'hotspot':weight,'pair_hash':1-weight}

    print(f"[run_target_benchmark] config={cfg_path}")
    print(f"[run_target_benchmark] pharmacophore={config.get('benchmark', {}).get('pharmacophore_path')}")
    print(f"[run_target_benchmark] required_groups={tuple(screening.HOTSPOT_REQUIRED_GROUPS)}")
    print(f"[run_target_benchmark] native_rerank_config={screening.DEFAULT_NATIVE_RERANK_CONFIG}")

    # The external-benchmark path does not consult the loaded `tables` argument.
    bench = run_benchmark({}, config)

    if args.output_dir:
        output_dir = Path(args.output_dir)
    else:
        rel = config.get("benchmark", {}).get("output_dir", "evidence/outputs/benchmark_target")
        output_dir = (ROOT / rel).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if not bench.summary_df.empty:
        bench.summary_df.to_csv(output_dir / "benchmark_summary.csv", index=False)
    if not bench.delta_df.empty:
        bench.delta_df.to_csv(output_dir / "benchmark_deltas.csv", index=False)
    if not bench.pairwise_df.empty:
        bench.pairwise_df.to_csv(output_dir / "benchmark_pairwise.csv", index=False)

    print(f"[run_target_benchmark] wrote summary to {output_dir}")
    cols = [c for c in ["method", "roc_auc", "pr_auc", "ef_1pct", "bedroc", "top10_recovery"] if c in bench.summary_df.columns]
    if cols:
        print(bench.summary_df[cols].to_string(index=False))
    execution.update(status='completed',elapsed_seconds=time.perf_counter()-started)
    execution_path.write_text(json.dumps(execution,indent=2))


if __name__ == "__main__":
    main()
