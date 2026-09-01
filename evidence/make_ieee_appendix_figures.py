#!/usr/bin/env python3
"""Build the appendix figures for the IEEE Access submission.

Same style contract as the main-text figures (``evidence/figstyle.py``):
Okabe-Ito palette, redundant encoding, 600 dpi PNG plus vector PDF.

Figures produced
----------------
figA1_upstream_vs_native   upstream screen metrics against native coverage
figA2_coverage_profile     distribution of native weighted coverage
figA3_top_ligand_overlay   best-scoring native overlay for the top ligand
figA4_top20_structures     structures of the top 20 native-ranked ligands

The previous versions of A1 and A3 relied on hue alone, including a red/green
contrast that is lost under the common forms of color-vision deficiency, and
A1 plotted 4,997 heavily overlapping points as raw scatter. Both are rebuilt
here: A1 as a density hexbin, A3 on the Okabe-Ito palette with marker shapes
carrying the feature family.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

import figstyle as fs

ROOT = Path(__file__).resolve().parent.parent

#: Output folder per journal target; ``main`` rebinds OUT from --target.
TARGET_DIRS = {
    "ieee": ROOT / "IEEE_Access_submission",
    "acs": ROOT / "ACS_Omega_resubmission",
}
OUT = TARGET_DIRS["ieee"]
BUNDLE = ROOT / "results" / "screening_full_1M_topological_hashed_native_terminal_bundle" / "analysis"
SCREEN_CSV = ROOT / "results" / "screening_full_1M_topological_hashed.csv"
FINAL_CSV = ROOT / "results" / "top_1000_glp1_mimetics_full_1M_topological_hashed_native_final.csv"


def _joined() -> pd.DataFrame:
    """Native-scored ligands with their upstream screen metrics.

    ``native_weighted_coverage_pct`` is the terminal native score; the file's
    ``weighted_coverage_pct`` is the Stage-3 geometric value and must not be
    substituted for it.
    """
    return pd.read_csv(
        BUNDLE / "ligand_best_native_mapping_summary.csv",
        usecols=[
            "ligand_id", "zinc_id", "native_weighted_coverage_pct", "stage3_screen_rank",
            "pair_hash_overlap_pct", "pair_hash_recall_pct",
            "hotspot_compatible_matches", "cascade_score_pct",
        ],
    )


def figA1_upstream_vs_native() -> None:
    """Upstream screen metrics against terminal native coverage.

    Drawn as a density hexbin: with 4,997 ligands on a narrow coverage range,
    a raw scatter mostly plots points on top of each other and hides where the
    mass of the distribution actually sits.
    """
    df = _joined()
    panels = [
        ("pair_hash_overlap_pct", "Pair-hash overlap (%)", 0.275, 0.250),
        ("pair_hash_recall_pct", "Pair-hash recall (%)", 0.292, 0.260),
        ("hotspot_compatible_matches", "Hotspot-compatible matches", 0.306, 0.221),
        ("cascade_score_pct", "Cascade score (%)", 0.113, 0.106),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(fs.PAGE_WIDTH, 4.4))
    for ax, (col, xlabel, pear, spear), letter in zip(axes.ravel(), panels, "abcd"):
        sub = df[[col, "native_weighted_coverage_pct"]].dropna()
        hb = ax.hexbin(
            sub[col], sub["native_weighted_coverage_pct"],
            gridsize=34, cmap="Blues", mincnt=1, linewidths=0.0,
        )
        # Least-squares trend, drawn in a color that reads against the blue map.
        m, b = np.polyfit(sub[col], sub["native_weighted_coverage_pct"], 1)
        xs = np.linspace(sub[col].min(), sub[col].max(), 50)
        ax.plot(xs, m * xs + b, color=fs.OKABE_ITO["vermillion"], linewidth=1.2,
                linestyle="--", zorder=4)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Native weighted coverage (%)")
        ax.set_title(f"Pearson $r$ = {pear:.3f},  Spearman $\\rho$ = {spear:.3f}", fontsize=7)
        cb = fig.colorbar(hb, ax=ax, pad=0.02, fraction=0.045)
        cb.set_label("ligands", fontsize=6.4)
        cb.ax.tick_params(labelsize=6)
        fs.panel_label(ax, f"({letter})", dx=-0.20, dy=1.16)

    fig.tight_layout(h_pad=1.5, w_pad=1.4)
    fs.save(fig, OUT, "figA1_upstream_vs_native")
    print("wrote figA1_upstream_vs_native")


def figA2_coverage_profile() -> None:
    """Distribution of native weighted coverage over the scored ligands."""
    native = pd.read_csv(
        BUNDLE / "ligand_best_native_mapping_summary.csv",
        usecols=["native_weighted_coverage_pct", "matched_reference_features"],
    )
    cov = native["native_weighted_coverage_pct"].to_numpy()

    fig, (ax, axb) = plt.subplots(1, 2, figsize=(fs.PAGE_WIDTH, 2.3))

    ax.hist(cov, bins=60, color=fs.OKABE_ITO["blue"], edgecolor="white", linewidth=0.3)
    for q, style, lab in ((np.median(cov), "-", "median"),
                          (np.percentile(cov, 25), ":", "IQR"),
                          (np.percentile(cov, 75), ":", None)):
        ax.axvline(q, color=fs.OKABE_ITO["vermillion"], linewidth=0.9, linestyle=style,
                   label=lab, zorder=3)
    ax.set_xlabel("Native weighted coverage (%)")
    ax.set_ylabel("Ligands")
    ax.legend(fontsize=6.4)
    ax.set_title(
        f"median {np.median(cov):.2f}%   IQR {np.percentile(cov, 25):.2f}–"
        f"{np.percentile(cov, 75):.2f}%   $n$ = {len(cov):,}",
        fontsize=6.8,
    )
    fs.panel_label(ax, "(a)", dx=-0.17)

    order = np.sort(cov)[::-1]
    axb.plot(np.arange(1, len(order) + 1), order, color=fs.OKABE_ITO["blue"], linewidth=1.2)
    axb.axvline(1000, color=fs.OKABE_ITO["orange"], linewidth=0.9, linestyle="--")
    axb.annotate("reported top 1,000", xy=(1000, order[999]), xytext=(1400, order[0] - 1.2),
                 fontsize=6.2, color="#3A3A3A",
                 arrowprops=dict(arrowstyle="->", lw=0.6, color="#4D4D4D"))
    axb.set_xlabel("Native rank")
    axb.set_ylabel("Native weighted coverage (%)")
    fs.panel_label(axb, "(b)", dx=-0.17)

    fig.tight_layout(w_pad=1.5)
    fs.save(fig, OUT, "figA2_coverage_profile")
    print("wrote figA2_coverage_profile")


def figA3_top_ligand_overlay() -> None:
    """Best-scoring native overlay for the top-ranked ligand.

    The retained GLP-1 peptide-contact reference features are shown with the
    ligand's matched feature centroids. Families are encoded by Okabe-Ito
    color and by marker shape, replacing the earlier red/green scheme.
    """
    ref = pd.read_csv(BUNDLE / "native_reference_features.csv")
    matches = pd.read_csv(BUNDLE / "microstate_native_feature_matches.csv")
    final = pd.read_csv(FINAL_CSV)
    top = final.iloc[0]
    top_id = top["zinc_id"]

    # Only the microstate that produced the ligand's reported best score.
    sub = matches[
        (matches["ligand_id"] == top_id)
        & (matches["microstate_id"] == top["native_microstate_id"])
    ]
    matched_ids = set(sub["reference_feature_id"])

    fam_style = {
        "polar": (fs.OKABE_ITO["orange"], "s", "Polar"),
        "anion": (fs.OKABE_ITO["blue"], "v", "Anion"),
        "cation": (fs.OKABE_ITO["purple"], "^", "Cation"),
        "hydrophobe": (fs.OKABE_ITO["green"], "o", "Hydrophobe"),
        "aromatic": (fs.OKABE_ITO["sky"], "h", "Aromatic"),
    }

    xyz = ref[["x", "y", "z"]].to_numpy()
    centered = xyz - xyz.mean(axis=0)
    _, _, vt = np.linalg.svd(centered, full_matrices=False)
    proj = centered @ vt[:2].T

    fig, ax = plt.subplots(figsize=(fs.COL_WIDTH, 3.0))
    for i, r in ref.reset_index(drop=True).iterrows():
        color, marker, _ = fam_style.get(
            r["family"], (fs.OKABE_ITO["grey"], "o", r["family"])
        )
        is_matched = r["feature_id"] in matched_ids
        ax.scatter(
            proj[i, 0], proj[i, 1],
            s=18 + 16 * r["weight"],
            c=color if is_matched else "white",
            marker=marker,
            edgecolors=color,
            linewidths=1.1 if is_matched else 0.7,
            alpha=1.0 if is_matched else 0.55,
            zorder=3 if is_matched else 2,
        )
        if is_matched:
            ax.annotate(str(r["residue_label"]), (proj[i, 0], proj[i, 1]),
                        textcoords="offset points", xytext=(5, 3.5),
                        fontsize=5.6, fontweight="bold", zorder=4)

    ax.set_xlabel("Principal axis 1 (Å)")
    ax.set_ylabel("Principal axis 2 (Å)")
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_title(f"{top_id}: {len(matched_ids)} of {len(ref)} reference features matched",
                 fontsize=6.8)
    handles = [
        Line2D([], [], marker=mk, color="none", markerfacecolor=c, markeredgecolor=c,
               markersize=4.5, label=lab)
        for c, mk, lab in fam_style.values()
    ] + [
        Line2D([], [], marker="o", color="none", markerfacecolor="white",
               markeredgecolor=fs.OKABE_ITO["grey"], markersize=4.5, label="Unmatched"),
    ]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.20),
              ncol=3, fontsize=5.9)
    fig.tight_layout()
    fs.save(fig, OUT, "figA3_top_ligand_overlay")
    print("wrote figA3_top_ligand_overlay")


def figA4_top20_structures() -> None:
    """Structure grid for the top 20 native-ranked ligands."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Draw
    from rdkit.Chem.Draw import rdMolDraw2D

    RDLogger.DisableLog("rdApp.*")
    final = pd.read_csv(FINAL_CSV).head(20)
    mols, legends = [], []
    for i, r in final.iterrows():
        m = Chem.MolFromSmiles(str(r["smiles"]))
        if m is None:
            continue
        Chem.rdDepictor.Compute2DCoords(m)
        mols.append(m)
        cov = r.get("weighted_coverage_pct", np.nan)
        legends.append(f"{i + 1}. {r['zinc_id']}\n{cov:.2f}% coverage")

    drawer = rdMolDraw2D.MolDraw2DCairo(2400, 3000, 480, 600)
    opts = drawer.drawOptions()
    opts.legendFontSize = 26
    opts.bondLineWidth = 3
    drawer.DrawMolecules(mols, legends=legends)
    drawer.FinishDrawing()
    png = drawer.GetDrawingText()
    (OUT / "figA4_top20_structures.png").write_bytes(png)

    # Stamp 600 dpi so the file declares the resolution IEEE Access asks for.
    from PIL import Image
    import io

    img = Image.open(io.BytesIO(png))
    img.save(OUT / "figA4_top20_structures.png", dpi=(600, 600))
    print("wrote figA4_top20_structures")


def main() -> None:
    global OUT
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--target", choices=sorted(TARGET_DIRS), default="ieee",
                    help="journal to size and write the figures for")
    args = ap.parse_args()

    OUT = TARGET_DIRS[args.target]
    fs.use_target(args.target)
    fs.apply()
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"target {args.target}: {OUT.name}, page width {fs.PAGE_WIDTH} in")
    figA1_upstream_vs_native()
    figA2_coverage_profile()
    figA3_top_ligand_overlay()
    figA4_top20_structures()


if __name__ == "__main__":
    main()
