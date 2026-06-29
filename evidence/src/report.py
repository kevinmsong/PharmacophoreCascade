"""
Report emitter: CSV + Markdown + LaTeX + PDF plots.

Entry point: emit_report(benchmark_result, ablation_result, claims_result, output_dir)
"""
from __future__ import annotations

import logging
import textwrap
from pathlib import Path
from typing import Dict, Optional

import matplotlib
matplotlib.use("Agg")  # non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd

from .benchmark import BenchmarkResult
from .ablation import AblationResult
from .claims import claims_to_summary_df

logger = logging.getLogger(__name__)

METHOD_DISPLAY_NAMES = {
    "full_cascade": "Full cascade",
    "stage3_only": "Stage-3 only",
    "standard_3d_pharmacophore": "Std. 3D pharm.",
    "active_state_docking": "Active-state docking",
    "state_preference_docking_sensitivity": "State-preference docking",
}

COLUMN_DISPLAY_NAMES = {
    "method": "Method",
    "method_a": "Method A",
    "method_b": "Method B",
    "metric": "Metric",
    "n_ranked": "Ranked",
    "ef_1pct": "EF1%",
    "ef_5pct": "EF5%",
    "roc_auc": "ROC-AUC",
    "pr_auc": "PR-AUC",
    "bedroc": "BEDROC",
    "top10_recovery": "Act@10",
    "top25_recovery": "Act@25",
    "top50_recovery": "Act@50",
    "scaffold_recovery_top10": "Scaf@10",
    "scaffold_recovery_top25": "Scaf@25",
    "scaffold_recovery_top50": "Scaf@50",
    "kendall_tau": "Kendall tau",
    "spearman_r": "Spearman r",
    "delta_mean": "Delta mean",
    "delta_ci_low": "95% CI low",
    "delta_ci_high": "95% CI high",
}


# ---------------------------------------------------------------------------
# CSV outputs
# ---------------------------------------------------------------------------

def _save_csv(df: pd.DataFrame, path: Path, label: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    logger.info("Wrote %s: %s (%d rows)", label, path.name, len(df))


# ---------------------------------------------------------------------------
# LaTeX table helper
# ---------------------------------------------------------------------------

def _df_to_latex(
    df: pd.DataFrame,
    caption: str,
    label: str,
    *,
    resize_to_textwidth: bool = False,
    size_command: str = r"\small",
    tabcolsep_pt: Optional[float] = None,
) -> str:
    """Convert a DataFrame to a LaTeX table string (booktabs style)."""
    pretty_df = df.copy()
    for col in ("method", "method_a", "method_b"):
        if col in pretty_df.columns:
            pretty_df[col] = pretty_df[col].map(
                lambda v: METHOD_DISPLAY_NAMES.get(v, str(v).replace("_", " "))
            )
    pretty_df = pretty_df.rename(
        columns={c: COLUMN_DISPLAY_NAMES.get(c, c) for c in pretty_df.columns}
    )

    col_fmt = "l" + "r" * (len(df.columns) - 1)
    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        rf"\caption{{{caption}}}",
        rf"\label{{{label}}}",
    ]
    if size_command:
        lines.insert(2, size_command)
    if tabcolsep_pt is not None:
        lines.append(rf"\setlength{{\tabcolsep}}{{{tabcolsep_pt}pt}}")
    if resize_to_textwidth:
        lines.append(r"\resizebox{\textwidth}{!}{%")
    lines += [
        rf"\begin{{tabular}}{{{col_fmt}}}",
        r"\toprule",
    ]

    # Header
    header = " & ".join(str(c).replace("%", r"\%").replace("_", r"\_") for c in pretty_df.columns) + r" \\"
    lines.append(header)
    lines.append(r"\midrule")

    # Rows
    for _, row in pretty_df.iterrows():
        def _fmt(v):
            if isinstance(v, float):
                if np.isnan(v):
                    return "--"
                return f"{v:.3f}"
            return str(v).replace("%", r"\%").replace("_", r"\_")

        lines.append(" & ".join(_fmt(v) for v in row) + r" \\")

    lines += [r"\bottomrule", r"\end{tabular}"]
    if resize_to_textwidth:
        lines.append(r"}")
    lines.append(r"\end{table}")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def _plot_rank_corr_heatmap(pairwise_df: pd.DataFrame, ax: plt.Axes) -> None:
    """Heatmap of Spearman r between all method pairs."""
    if pairwise_df.empty:
        ax.text(0.5, 0.5, "No pairwise data", ha="center", va="center")
        return

    methods = sorted(set(pairwise_df["method_a"].tolist() + pairwise_df["method_b"].tolist()))
    n = len(methods)
    mat = np.full((n, n), np.nan)
    np.fill_diagonal(mat, 1.0)

    idx = {m: i for i, m in enumerate(methods)}
    for _, row in pairwise_df.iterrows():
        i, j = idx[row["method_a"]], idx[row["method_b"]]
        v = row.get("spearman_r", np.nan)
        mat[i, j] = v
        mat[j, i] = v

    im = ax.imshow(mat, vmin=-1, vmax=1, cmap="RdYlGn")
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    short = [m.replace("_", "\n") for m in methods]
    ax.set_xticks(range(n)); ax.set_xticklabels(short, fontsize=7, rotation=45, ha="right")
    ax.set_yticks(range(n)); ax.set_yticklabels(short, fontsize=7)
    ax.set_title("Spearman r (method pairs)", fontsize=9)

    for i in range(n):
        for j in range(n):
            v = mat[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=6)


def _plot_jaccard_heatmap(pairwise_df: pd.DataFrame, ax: plt.Axes, k: int = 100) -> None:
    """Heatmap of Jaccard top-k between method pairs."""
    col = f"jaccard_top{k}"
    if pairwise_df.empty or col not in pairwise_df.columns:
        ax.text(0.5, 0.5, f"No Jaccard top-{k} data", ha="center", va="center")
        return

    methods = sorted(set(pairwise_df["method_a"].tolist() + pairwise_df["method_b"].tolist()))
    n = len(methods)
    mat = np.full((n, n), np.nan)
    np.fill_diagonal(mat, 1.0)

    idx = {m: i for i, m in enumerate(methods)}
    for _, row in pairwise_df.iterrows():
        i, j = idx[row["method_a"]], idx[row["method_b"]]
        v = row.get(col, np.nan)
        mat[i, j] = v
        mat[j, i] = v

    im = ax.imshow(mat, vmin=0, vmax=1, cmap="Blues")
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    short = [m.replace("_", "\n") for m in methods]
    ax.set_xticks(range(n)); ax.set_xticklabels(short, fontsize=7, rotation=45, ha="right")
    ax.set_yticks(range(n)); ax.set_yticklabels(short, fontsize=7)
    ax.set_title(f"Jaccard top-{k} (method pairs)", fontsize=9)

    for i in range(n):
        for j in range(n):
            v = mat[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=6)


def _plot_cascade_vs_native_scatter(
    tables: Dict[str, pd.DataFrame],
    ax: plt.Axes,
) -> None:
    """Scatter: Stage-3 weighted coverage vs native weighted coverage."""
    if "rank_comparison" not in tables:
        ax.text(0.5, 0.5, "rank_comparison table not loaded", ha="center", va="center")
        return

    rc = tables["rank_comparison"]
    x = rc["weighted_coverage_pct"].values
    y = rc["native_weighted_coverage_pct"].values

    ax.scatter(x, y, alpha=0.25, s=6, color="steelblue", rasterized=True)

    # Regression line
    m, b, r, p, _ = __import__("scipy.stats", fromlist=["linregress"]).linregress(x, y)
    x_line = np.linspace(x.min(), x.max(), 100)
    ax.plot(x_line, m * x_line + b, color="crimson", linewidth=1.5,
            label=f"r={r:.3f}, p={p:.2g}")
    ax.set_xlabel("Stage-3 weighted coverage (%)", fontsize=8)
    ax.set_ylabel("Native weighted coverage (%)", fontsize=8)
    ax.set_title("Stage-3 vs Native Coverage", fontsize=9)
    ax.legend(fontsize=7)


def _plot_rank_shift_hist(tables: Dict[str, pd.DataFrame], ax: plt.Axes) -> None:
    """Histogram of absolute rank shifts between stage-3 and native rank."""
    if "rank_comparison" not in tables:
        ax.text(0.5, 0.5, "rank_comparison table not loaded", ha="center", va="center")
        return

    rc = tables["rank_comparison"]
    shifts = rc["abs_rank_shift"].dropna()
    ax.hist(shifts, bins=50, color="steelblue", edgecolor="white", alpha=0.8)
    med = shifts.median()
    ax.axvline(med, color="crimson", linewidth=1.5, linestyle="--",
               label=f"median={med:.0f}")
    ax.set_xlabel("Absolute rank shift (screen→native)", fontsize=8)
    ax.set_ylabel("Count", fontsize=8)
    ax.set_title("Rank Shift Distribution", fontsize=9)
    ax.legend(fontsize=7)


def _plot_external_metric_bars(
    summary_df: pd.DataFrame,
    ax: plt.Axes,
    metric: str,
    title: str,
) -> None:
    if summary_df.empty or metric not in summary_df.columns:
        ax.text(0.5, 0.5, f"No {metric} data", ha="center", va="center")
        return

    plot_df = summary_df.copy()
    plot_df = plot_df.loc[plot_df[metric].notna()].copy()
    if plot_df.empty:
        ax.text(0.5, 0.5, f"No {metric} data", ha="center", va="center")
        return

    methods = plot_df["method"].tolist()
    values = plot_df[metric].astype(float).values
    yerr = None
    low_col = f"{metric}_ci_low"
    high_col = f"{metric}_ci_high"
    if low_col in plot_df.columns and high_col in plot_df.columns:
        low = plot_df[low_col].astype(float).values
        high = plot_df[high_col].astype(float).values
        yerr = np.vstack([np.clip(values - low, 0, None), np.clip(high - values, 0, None)])

    x = np.arange(len(methods))
    ax.bar(x, values, yerr=yerr, color="steelblue", edgecolor="white", alpha=0.9, capsize=3)
    ax.set_xticks(x)
    ax.set_xticklabels(
        [METHOD_DISPLAY_NAMES.get(m, m.replace("_", " ")) for m in methods],
        rotation=0,
        ha="center",
        fontsize=14,
    )
    ax.set_title(title, fontsize=17, pad=10)
    ax.set_ylabel(metric, fontsize=14)
    ax.tick_params(axis="y", labelsize=13)
    ax.set_ylim(bottom=0)


def make_plots(
    bench: BenchmarkResult,
    tables: Dict[str, pd.DataFrame],
    output_dir: Path,
) -> None:
    """Generate and save benchmark_plots.pdf."""
    fig, axes = plt.subplots(2, 2, figsize=(15, 11.5))
    fig.suptitle("GLP-1R Cascade Evidence Package: Benchmark Summary", fontsize=18, y=0.99)

    if bench.library_df is not None:
        _plot_external_metric_bars(bench.summary_df, axes[0, 0], "roc_auc", "ROC-AUC")
        _plot_external_metric_bars(bench.summary_df, axes[0, 1], "pr_auc", "PR-AUC")
        _plot_external_metric_bars(bench.summary_df, axes[1, 0], "ef_1pct", "EF1%")
        _plot_external_metric_bars(bench.summary_df, axes[1, 1], "bedroc", "BEDROC (alpha=20)")
    else:
        _plot_rank_corr_heatmap(bench.pairwise_df, axes[0, 0])
        _plot_jaccard_heatmap(bench.pairwise_df, axes[0, 1], k=100)
        _plot_cascade_vs_native_scatter(tables, axes[1, 0])
        _plot_rank_shift_hist(tables, axes[1, 1])

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    out_path = output_dir / "benchmark_plots.pdf"
    fig.savefig(out_path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    logger.info("Wrote benchmark_plots.pdf")


def make_ablation_plots(ablation: AblationResult, output_dir: Path) -> None:
    """Sensitivity heatmap: cascade weights vs Jaccard top-100."""
    if ablation.sensitivity_df.empty:
        return

    weight_df = ablation.sensitivity_df[
        ablation.sensitivity_df.get("sweep_type", pd.Series()) == "cascade_weights"
    ] if "sweep_type" in ablation.sensitivity_df.columns else pd.DataFrame()

    if weight_df.empty:
        return

    fig, ax = plt.subplots(figsize=(7, 4))
    jaccard_col = "jaccard_top100" if "jaccard_top100" in weight_df.columns else None
    if jaccard_col:
        labels = [f"({r['hotspot_weight']:.1f},{r['pair_hash_weight']:.1f})"
                  for _, r in weight_df.iterrows()]
        vals = weight_df[jaccard_col].values
        ax.bar(labels, vals, color="steelblue")
        ax.axhline(1.0, color="crimson", linestyle="--", linewidth=1)
        ax.set_ylim(0, 1.05)
        ax.set_xlabel("(hotspot_w, pair_hash_w)", fontsize=8)
        ax.set_ylabel("Jaccard top-100 vs baseline", fontsize=8)
        ax.set_title("Cascade Weight Sensitivity", fontsize=9)
        plt.xticks(rotation=30, fontsize=7)

    out_path = output_dir / "ablation_sensitivity.pdf"
    fig.savefig(out_path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    logger.info("Wrote ablation_sensitivity.pdf")


# ---------------------------------------------------------------------------
# Markdown report
# ---------------------------------------------------------------------------

def _verdict_icon(verdict: str) -> str:
    return {"supports": "[SUPPORTS]", "contradicts": "[CONTRADICTS]", "inconclusive": "[INCONCLUSIVE]"}.get(
        verdict, verdict
    )


def make_markdown_report(
    bench: BenchmarkResult,
    ablation: AblationResult,
    claims_results: Dict[str, dict],
    output_dir: Path,
) -> None:
    lines = [
        "# GLP-1R Evidence Package: Benchmark Report\n",
        "## 1. Baseline Method Comparison\n",
    ]

    if not bench.summary_df.empty:
        lines.append(bench.summary_df.to_markdown(index=False, floatfmt=".3f"))
        lines.append("\n")

    lines.append("\n## 2. Pairwise Rank Correlation\n")
    if not bench.pairwise_df.empty:
        cols = ["method_a", "method_b", "n_shared", "kendall_tau", "spearman_r"]
        jaccard_cols = [c for c in bench.pairwise_df.columns if c.startswith("jaccard")]
        lines.append(bench.pairwise_df[cols + jaccard_cols].to_markdown(index=False, floatfmt=".3f"))
        lines.append("\n")

    if not bench.delta_df.empty:
        lines.append("\n## 3. Paired Bootstrap Deltas vs Full Cascade\n")
        lines.append(bench.delta_df.to_markdown(index=False, floatfmt=".3f"))
        lines.append("\n")

    lines.append("\n## 4. Ablation Results\n")
    if not ablation.ablation_df.empty:
        lines.append(ablation.ablation_df.to_markdown(index=False, floatfmt=".3f"))
        lines.append("\n")
    else:
        lines.append("_No ablations run._\n")

    lines.append("\n## 5. Sensitivity Sweeps\n")
    if not ablation.sensitivity_df.empty:
        lines.append(ablation.sensitivity_df.to_markdown(index=False, floatfmt=".3f"))
        lines.append("\n")
    else:
        lines.append("_No sweeps run._\n")

    lines.append("\n## 6. Manuscript Claim Evidence\n\n")
    lines.append("| Claim | Verdict | Summary |\n|---|---|---|\n")
    for label, res in claims_results.items():
        v = _verdict_icon(res["verdict"])
        narr = res["narrative"].replace("|", r"\|")
        lines.append(f"| {label} | {v} | {narr} |\n")

    lines.append("\n### Claim Detail Tables\n")
    for label, res in claims_results.items():
        lines.append(f"\n#### {label}\n")
        edf = res["evidence_df"]
        if not edf.empty:
            lines.append(edf.to_markdown(index=False, floatfmt=".4f"))
        lines.append("\n")

    out_path = output_dir / "benchmark_report.md"
    out_path.write_text("\n".join(lines), encoding="utf-8")
    logger.info("Wrote benchmark_report.md")


# ---------------------------------------------------------------------------
# LaTeX tables
# ---------------------------------------------------------------------------

def make_latex_tables(
    bench: BenchmarkResult,
    ablation: AblationResult,
    output_dir: Path,
) -> None:
    parts = []
    summary_tex = None
    summary_core_tex = None
    summary_recovery_tex = None
    pairwise_tex = None
    delta_tex = None
    ablation_tex = None

    if not bench.summary_df.empty:
        display_cols = [
            "method",
            "n_ranked",
            "ef_1pct",
            "ef_5pct",
            "roc_auc",
            "pr_auc",
            "bedroc",
            "top10_recovery",
            "top25_recovery",
            "top50_recovery",
            "scaffold_recovery_top10",
            "scaffold_recovery_top25",
            "scaffold_recovery_top50",
        ]
        display_cols = [c for c in display_cols if c in bench.summary_df.columns]
        summary_caption = (
            "Retrospective external benchmark performance on the GLP-1R active-vs-decoy set. "
            "Metrics summarize discrimination and early enrichment for the full cascade and the simpler baseline methods. "
            "BEDROC used alpha=20."
        )
        summary_tex = _df_to_latex(
            bench.summary_df[display_cols],
            caption=summary_caption,
            label="tab:benchmark_summary",
            resize_to_textwidth=True,
        )

        core_cols = [
            "method",
            "n_ranked",
            "ef_1pct",
            "ef_5pct",
            "roc_auc",
            "pr_auc",
            "bedroc",
        ]
        recovery_cols = [
            "method",
            "top10_recovery",
            "top25_recovery",
            "top50_recovery",
            "scaffold_recovery_top10",
            "scaffold_recovery_top25",
            "scaffold_recovery_top50",
        ]
        core_cols = [c for c in core_cols if c in bench.summary_df.columns]
        recovery_cols = [c for c in recovery_cols if c in bench.summary_df.columns]
        core_df = bench.summary_df[core_cols].copy()
        if "method" in core_df.columns:
            core_df["method"] = core_df["method"].replace(
                {"standard_3d_pharmacophore": "Standard 3D pharmacophore"}
            )
        summary_core_tex = _df_to_latex(
            core_df,
            caption=summary_caption,
            label="tab:benchmark_summary_core",
            size_command=r"\normalsize",
            tabcolsep_pt=4,
        )
        recovery_df = bench.summary_df[recovery_cols].copy()
        if "method" in recovery_df.columns:
            recovery_df["method"] = recovery_df["method"].replace(
                {"standard_3d_pharmacophore": "Standard 3D pharmacophore"}
            )
        summary_recovery_tex = _df_to_latex(
            recovery_df,
            caption="Early active and scaffold recovery on the retrospective GLP-1R benchmark. Values denote the fraction of benchmark actives or active Murcko scaffolds recovered within the top-k ranked molecules.",
            label="tab:benchmark_summary_recovery",
            size_command=r"\normalsize",
            tabcolsep_pt=4,
        )
        parts.extend([summary_core_tex, summary_recovery_tex])

    if not bench.pairwise_df.empty:
        display = ["method_a", "method_b", "kendall_tau", "spearman_r"]
        display += [c for c in bench.pairwise_df.columns if c.startswith("jaccard_top")]
        display = [c for c in display if c in bench.pairwise_df.columns]
        pairwise_tex = _df_to_latex(
            bench.pairwise_df[display],
            caption="Pairwise rank correlation and Jaccard overlap between ranking methods.",
            label="tab:benchmark_pairwise",
        )
        parts.append(pairwise_tex)

    if not bench.delta_df.empty:
        delta_display = [
            "method",
            "metric",
            "delta_mean",
            "delta_ci_low",
            "delta_ci_high",
        ]
        delta_display = [c for c in delta_display if c in bench.delta_df.columns]
        delta_tex = _df_to_latex(
            bench.delta_df[delta_display],
            caption="Paired bootstrap deltas versus the full cascade.",
            label="tab:benchmark_deltas",
        )
        parts.append(delta_tex)

    if not ablation.ablation_df.empty:
        ablation_tex = _df_to_latex(
            ablation.ablation_df,
            caption="Rank stability under component ablations.",
            label="tab:ablation_results",
        )
        parts.append(ablation_tex)

    out_path = output_dir / "benchmark_table.tex"
    out_path.write_text("\n\n".join(parts), encoding="utf-8")
    logger.info("Wrote benchmark_table.tex")
    if summary_tex is not None:
        (output_dir / "benchmark_summary_table.tex").write_text(summary_tex, encoding="utf-8")
        logger.info("Wrote benchmark_summary_table.tex")
    if summary_core_tex is not None:
        (output_dir / "benchmark_summary_core_table.tex").write_text(summary_core_tex, encoding="utf-8")
        logger.info("Wrote benchmark_summary_core_table.tex")
    if summary_recovery_tex is not None:
        (output_dir / "benchmark_summary_recovery_table.tex").write_text(summary_recovery_tex, encoding="utf-8")
        logger.info("Wrote benchmark_summary_recovery_table.tex")
    if pairwise_tex is not None:
        (output_dir / "benchmark_pairwise_table.tex").write_text(pairwise_tex, encoding="utf-8")
        logger.info("Wrote benchmark_pairwise_table.tex")
    if delta_tex is not None:
        (output_dir / "benchmark_delta_table.tex").write_text(delta_tex, encoding="utf-8")
        logger.info("Wrote benchmark_delta_table.tex")
    if ablation_tex is not None:
        (output_dir / "ablation_results_table.tex").write_text(ablation_tex, encoding="utf-8")
        logger.info("Wrote ablation_results_table.tex")


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def emit_report(
    bench: BenchmarkResult,
    ablation: AblationResult,
    claims_results: Dict[str, dict],
    tables: Dict[str, pd.DataFrame],
    output_dir: Path,
) -> None:
    """Write all output files to output_dir."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # CSVs
    if not bench.summary_df.empty:
        _save_csv(bench.summary_df, output_dir / "benchmark_summary.csv", "benchmark_summary")
    if not bench.delta_df.empty:
        _save_csv(bench.delta_df, output_dir / "benchmark_deltas.csv", "benchmark_deltas")
    if not bench.pairwise_df.empty:
        _save_csv(bench.pairwise_df, output_dir / "benchmark_pairwise.csv", "benchmark_pairwise")
    if not ablation.ablation_df.empty:
        _save_csv(ablation.ablation_df, output_dir / "ablation_results.csv", "ablation_results")
    if not ablation.sensitivity_df.empty:
        _save_csv(ablation.sensitivity_df, output_dir / "sensitivity_results.csv", "sensitivity_results")

    claims_summary = claims_to_summary_df(claims_results)
    _save_csv(claims_summary, output_dir / "claims_summary.csv", "claims_summary")

    # Per-claim evidence tables
    for label, res in claims_results.items():
        edf = res["evidence_df"]
        if not edf.empty:
            _save_csv(edf, output_dir / f"{label}_evidence.csv", label)

    # Markdown
    make_markdown_report(bench, ablation, claims_results, output_dir)

    # LaTeX
    make_latex_tables(bench, ablation, output_dir)

    # Plots
    try:
        make_plots(bench, tables, output_dir)
        make_ablation_plots(ablation, output_dir)
    except Exception as exc:
        logger.warning("Plot generation failed: %s", exc)
