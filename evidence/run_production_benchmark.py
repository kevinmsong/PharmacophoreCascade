#!/usr/bin/env python3
"""
Production-constrained benchmark with stage-by-stage active/decoy survival.

Embeds a system's 10 actives + 300 matched decoys in a large ZINC background and
runs the real production pipeline (the same Stage 0-2 gate, the 5% Stage-3
shortlist, and the native-pool quotas/scaffold caps used in the million-compound
screen), then reports how many actives and decoys survive Stage 0, Stage 1,
Stage 2, the Stage-3 shortlist, native-pool selection, and the final ranking.

Run per system (separate process so the required-group gate + native reference
are correct):

    python evidence/run_production_benchmark.py \
        --benchmark-config evidence/configs/benchmark.yaml \
        --background-size 50000 --output-dir evidence/outputs/production_glp1r \
        [--native-rerank-config <study.yaml>] [--required-groups ""]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(HERE))


def _gate_pass_status(status: str) -> str:
    """Map an engine status to the last stage a molecule cleared."""
    s = str(status)
    if s in ("invalid", "property_filtered"):
        return "stage0_fail"
    if s == "chemistry_filtered":
        return "stage0_fail"
    if s == "hotspot_filtered":
        return "stage1_fail"
    if s in ("candidate", "warn_only_flagged"):
        return "candidate"
    return "other"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--benchmark-config", required=True)
    ap.add_argument("--native-rerank-config", default=None)
    ap.add_argument("--required-groups", default=None)
    ap.add_argument("--background-size", type=int, default=50000)
    ap.add_argument("--shortlist-fraction", type=float, default=0.05)
    ap.add_argument("--max-zinc-candidates", type=int, default=400000)
    ap.add_argument("--zinc-dir", default=str(ROOT / "tmp" / "zinc"))
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--seed", type=int, default=20260601)
    args = ap.parse_args()

    if args.required_groups is not None:
        os.environ["CASCADIA_HOTSPOT_REQUIRED_GROUPS"] = args.required_groups
    import run_optimized_1M_topological_hashed_screening as screening
    if args.native_rerank_config:
        screening.DEFAULT_NATIVE_RERANK_CONFIG = Path(args.native_rerank_config)
    from src.benchmark_library_builder import iter_zinc_smiles

    cfg = yaml.safe_load(open(args.benchmark_config))
    bench = cfg["benchmark"]
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ---- labeled set ----
    lib = pd.read_csv(ROOT / bench["library_csv"])
    labeled = lib[["ligand_id", "smiles", "label"]].copy()
    labeled_ids = set(labeled["ligand_id"].astype(str))
    n_act = int((labeled["label"] == "active").sum())
    print(f"labeled: {len(labeled)} ({n_act} actives)")

    # ---- background from ZINC (distinct from labeled decoys) ----
    rng = np.random.default_rng(args.seed)
    zinc_paths = sorted(Path(args.zinc_dir).glob("*.smi.gz"))
    seen = set()
    bg_rows = []
    target = args.background_size
    # reservoir-ish: skip a random stride to diversify tranches
    stride = max(1, args.max_zinc_candidates // max(target, 1))
    for i, (zid, smi) in enumerate(iter_zinc_smiles(zinc_paths)):
        if i % stride != 0:
            continue
        if zid in labeled_ids or zid in seen:
            continue
        seen.add(zid)
        bg_rows.append({"ligand_id": zid, "smiles": smi, "label": "background"})
        if len(bg_rows) >= target:
            break
    background = pd.DataFrame(bg_rows)
    combined = pd.concat([labeled, background], ignore_index=True)
    print(f"combined library: {len(combined)} (background {len(background)})")

    # ---- Stage 0-2 ----
    pharm = screening.load_pharmacophore(str(ROOT / bench["pharmacophore_path"]))
    views = screening.precompute_pharmacophore_views(
        pharm,
        top_hotspots=int(bench.get("top_hotspots", 25)),
        pair_features=int(bench.get("pair_features", 24)),
        rerank_query_features=int(bench.get("rerank_query_features", 28)),
        pair_hash_mode=str(bench.get("pair_hash_mode", "precision_5bin")),
        native_support_max_residues=int(bench.get("native_support_max_residues", 0)),
    )
    rows = []
    for _, r in combined.iterrows():
        record, status, _ = screening.evaluate_topological_candidate_with_timings(
            r["smiles"], r["ligand_id"], views,
            hotspot_min_exact=int(bench.get("hotspot_min_exact", 3)),
            hotspot_min_groups=int(bench.get("hotspot_min_groups", 2)),
            pair_hash_mode=str(bench.get("pair_hash_mode", "precision_5bin")),
            chemistry_gate_mode=str(bench.get("chemistry_gate_mode", "strict")),
        )
        payload = {"ligand_id": r["ligand_id"], "smiles": r["smiles"], "label": r["label"],
                   "status": status, "gate": _gate_pass_status(status),
                   "cascade_score_pct": float((record or {}).get("cascade_score_pct", 0.0) or 0.0)}
        rows.append(payload)
    ev = pd.DataFrame(rows)
    ev.to_csv(out_dir / "stage012_evaluation.csv", index=False)

    candidates = ev[ev["gate"] == "candidate"].copy()
    # Stage-3 shortlist = top fraction of candidates by cascade score
    n_short = max(1, int(np.ceil(len(candidates) * args.shortlist_fraction)))
    shortlist = candidates.sort_values(["cascade_score_pct", "ligand_id"], ascending=[False, True]).head(n_short).copy()
    short_ids = set(shortlist["ligand_id"].astype(str))

    # ---- Stage 3 rerank on the shortlist ----
    pseudo = []
    for rank, (_, r) in enumerate(shortlist.iterrows(), start=1):
        mol = screening.standardize_screening_molecule(r["smiles"])
        props = screening.calculate_properties(mol) if mol is not None else {}
        pseudo.append({
            "zinc_id": r["ligand_id"], "smiles": r["smiles"],
            "mw": props.get("mw", 0.0), "logp": props.get("logp", 0.0),
            "hbd": props.get("hbd", 0), "hba": props.get("hba", 0),
            "rotatable_bonds": props.get("rotatable_bonds", 0), "lipinski": 1,
            "hotspot_weighted_pct": 0.0, "hotspot_bits_matched": 0, "hotspot_exact_matches": 0,
            "hotspot_compatible_matches": 0, "hotspot_exact_residue_count": 0,
            "hotspot_group_count": 0, "hotspot_required_groups_pass": 0,
            "pair_hash_overlap_pct": 0.0, "pair_hash_recall_pct": 0.0, "pair_hash_precision_pct": 0.0,
            "cascade_score_pct": float(r["cascade_score_pct"]), "shortlist_rank": rank,
        })
    reranked = screening.rerank_shortlist_parallel(
        pseudo, views, workers=int(bench.get("workers", 4)),
        rerank_conformers=int(bench.get("rerank_conformers", 16)),
        tolerance=float(bench.get("rerank_tolerance", 2.75)),
        pair_tolerance=float(bench.get("rerank_pair_tolerance", 2.75)),
        rerank_score_mode=str(bench.get("rerank_score_mode", "stage3_only")),
        report_interval=500,
    )
    stage3 = pd.DataFrame(reranked, columns=screening.final_output_columns())
    stage3 = screening.sort_stage3_final_hits(stage3)
    stage3.insert(0, "stage3_screen_rank", np.arange(1, len(stage3) + 1, dtype=int))

    # ---- Native terminal rerank with production caps (1M-run settings) ----
    pa = {k: out_dir / f"prod_{k}.csv" for k in
          ["native_scored", "native_final", "shortlist", "full", "top100", "run_summary"]}
    pa["native_bundle_root"] = out_dir / "prod_native_bundle"
    pa["pharmacophore_plot"] = out_dir / "prod_pharm.png"
    pa["top20_plot"] = out_dir / "prod_top20.png"
    pa["distributions_plot"] = out_dir / "prod_dist.png"
    prod_args = SimpleNamespace(
        native_candidate_pool_k=20000, native_selection_top_k=5000,
        native_final_top_k=1000, native_rerank_max_per_scaffold=8,
        native_rerank_pair_tolerance=3.0, final_rank_mode="native_first",
    )
    native_final, _ = screening.run_native_terminal_rerank(stage3_df=stage3, output_paths=pa, args=prod_args)
    native_pool_ids = set()
    sel_csv = out_dir / "prod_native_bundle" / "input" / "native_selected_input.csv"
    if sel_csv.exists():
        native_pool_ids = set(pd.read_csv(sel_csv)["zinc_id"].astype(str))
    final_ids = set(native_final["zinc_id"].astype(str)) if not native_final.empty else set()

    # ---- stage-by-stage survival ----
    def surv(mask_ids: set) -> dict:
        sub = labeled[labeled["ligand_id"].astype(str).isin(mask_ids)]
        return {"actives": int((sub["label"] == "active").sum()),
                "decoys": int((sub["label"] == "decoy").sum())}
    stage0_ids = set(ev[ev["gate"].isin(["candidate", "stage1_fail"])]["ligand_id"].astype(str))  # passed Stage 0
    stage1_ids = set(candidates["ligand_id"].astype(str))  # passed hotspot gate
    table = []
    total_bg = len(background)
    for name, ids in [("input", set(combined["ligand_id"].astype(str))),
                      ("stage0_pass", stage0_ids), ("stage1_pass", stage1_ids),
                      ("stage3_shortlist", short_ids), ("native_pool", native_pool_ids),
                      ("final_ranked", final_ids)]:
        s = surv(ids)
        # background survivors among the labeled-excluded set
        bg_surv = int(background["ligand_id"].astype(str).isin(ids).sum())
        table.append({"stage": name, "actives": s["actives"], "decoys": s["decoys"],
                      "background": bg_surv,
                      "active_survival_pct": round(100 * s["actives"] / max(n_act, 1), 1)})
    surv_df = pd.DataFrame(table)
    surv_df.to_csv(out_dir / "stage_survival.csv", index=False)

    # final enrichment: ranks of actives in the native_final ranking over the combined pool
    fr = native_final.copy()
    if not fr.empty:
        fr["label"] = fr["zinc_id"].astype(str).map(
            labeled.set_index(labeled["ligand_id"].astype(str))["label"].to_dict()).fillna("background")
        act_final_ranks = sorted(fr.index[fr["label"] == "active"].tolist())
    else:
        act_final_ranks = []
    summary = {
        "background_size": total_bg, "shortlist_fraction": args.shortlist_fraction,
        "n_candidates": int(len(candidates)), "n_shortlist": int(n_short),
        "actives_in_final": int((surv_df.set_index("stage").loc["final_ranked", "actives"])),
        "stage_survival": table,
    }
    (out_dir / "production_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print("\n=== stage-by-stage survival ===")
    print(surv_df.to_string(index=False))


if __name__ == "__main__":
    main()
