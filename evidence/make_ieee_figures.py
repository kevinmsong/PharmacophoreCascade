#!/usr/bin/env python3
"""Build every main-text figure for the IEEE Access submission.

All figures share ``evidence/figstyle.py``: one Okabe-Ito colorblind-safe
palette, one type scale, IEEE Access column widths, 600 dpi PNG plus vector
PDF. Series are separated by hatch and marker as well as hue so the figures
survive grayscale reproduction.

Figures produced
----------------
fig1_cascade              staged architecture, annotated input -> output
fig2_pharmacophore        receptor-side GLP-1R pharmacophore
fig3_glp1r_benchmark      GLP-1R enrichment, four methods with bootstrap CIs
fig4_cross_system         cross-system enrichment, cascade vs native-only
fig6_efficiency           efficiency-retention trade-off of the shortlist
fig5_attrition            where in-domain actives are lost, and to what
fig7_docking              orthogonal docking with the paired Wilcoxon test

Run ``analyze_active_attrition.py``, ``analyze_efficiency_retention.py``, and
``analyze_equivalence.py`` first; this script only renders their outputs.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

import figstyle as fs

ROOT = Path(__file__).resolve().parent.parent
EV = ROOT / "evidence" / "outputs"

#: Output folder per journal target; ``main`` rebinds OUT from --target.
TARGET_DIRS = {
    "ieee": ROOT / "IEEE_Access_submission",
    "acs": ROOT / "ACS_Omega_resubmission",
}
OUT = TARGET_DIRS["ieee"]

SYSTEM_ORDER = ["GLP-1R", "GHSR", "NTSR1", "MDM2-p53"]

#: Axis labels use the en dash the manuscript sets for the protein-protein pair;
#: SYSTEM_ORDER keeps the hyphen because it also indexes the released tables.
DISPLAY = {"MDM2-p53": "MDM2–p53"}


def shown(system: str) -> str:
    return DISPLAY.get(system, system)


# --------------------------------------------------------------------------
# Figure 1 - staged architecture
# --------------------------------------------------------------------------
def fig1_cascade() -> None:
    """Funnel schematic with an explicit input -> output line per stage.

    The reviewer asked for each stage's inputs and outputs to be stated in one
    sentence, so the schematic carries them rather than deferring to the text.
    """
    screen = fs.OKABE_ITO["blue"]
    native = fs.OKABE_ITO["green"]
    edge = fs.OKABE_ITO["black"]

    # (title, input -> output line, surviving count, color, y, box width)
    rows = [
        ("Input library", "1,000,000 ZINC compounds (tranches H17–H20)", "1,000,000", "#4D4D4D", 0.955, 0.50),
        ("Stage 0 · Standardization and property limits",
         "in: SMILES  →  out: property-qualified molecules; alerts are annotations", "956,240", screen, 0.855, 0.80),
        ("Stage 1 · Receptor hotspot compatibility",
         "in: typed ligand features  →  out: hotspot fraction $H$, pass/fail gate", "948,971", screen, 0.755, 0.80),
        ("Stage 2 · Typed pair-hash comparison",
         "in: feature pairs  →  out: pair overlap $O_\\mathrm{pair}$ for the cascade score", "no gate", screen, 0.655, 0.78),
        ("Cascade-score shortlist (5%; minimum 1,000)",
         "in: ranked survivors  →  out: the only molecules given 3D treatment", "47,812", screen, 0.555, 0.72),
        ("Stage 3 · Conformer-level geometric rerank",
         "in: 16 ETKDG conformers each  →  out: anchor coverage $\\times$ geometry penalty", "47,689", screen, 0.455, 0.76),
        ("Native branch · diversified candidate pool",
         "in: Stage-3 rank, hotspot breadth, native support  →  out: capped pool", "20,000 → 5,000", native, 0.330, 0.76),
        ("Native peptide-contact rescoring",
         "in: prepared microstates  →  out: native weighted coverage $C_\\mathrm{native}$", "4,997 scored", native, 0.230, 0.72),
        ("Final ranked output",
         "in: native scores  →  out: ranked peptide-mimicking candidates", "1,000", edge, 0.130, 0.56),
    ]

    # The ACS revision is tied to the completed fresh execution, not hand-entered
    # archived counts. Retain the historical rendering for the other journal.
    fresh = ROOT / "results/absolute_floor_1000/screening_full_1M_floor1000_run_summary.json"
    if OUT == TARGET_DIRS["acs"] and fresh.exists():
        run = json.loads(fresh.read_text())
        counts = run["counts"]
        native_counts = run["native_rerank"]["counts"]
        measured = [f"{counts['total_scanned']:,}", f"{counts['property_pass']:,}",
                    f"{counts['hotspot_pass']:,}", "no gate", f"{counts['shortlist_size']:,}",
                    f"{counts['final_hits']:,}",
                    f"{native_counts['candidate_pool_size']:,} → {native_counts['selected_ligands']:,}",
                    f"{native_counts['successful_best_ligands']:,} scored", f"{native_counts['final_rows']:,}"]
        rows = [(title, sub, count, color, y, width)
                for (title, sub, _, color, y, width), count in zip(rows, measured)]
    elif OUT != TARGET_DIRS["acs"]:
        rows[4] = ("Cascade-score shortlist (top 5%)", *rows[4][1:])

    cx = 0.42
    fig, ax = plt.subplots(figsize=(fs.PAGE_WIDTH, 6.0))
    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(0.06, 1.0)
    ax.axis("off")
    ax.grid(False)

    for title, sub, count, color, y, width in rows:
        x0 = cx - width / 2
        ax.add_patch(
            FancyBboxPatch(
                (x0, y - 0.036),
                width,
                0.072,
                boxstyle="round,pad=0.004,rounding_size=0.012",
                linewidth=0.8,
                edgecolor=color,
                facecolor=color + "1A",
            )
        )
        ax.text(cx, y + 0.014, title, ha="center", va="center", fontsize=7.4, fontweight="bold")
        ax.text(cx, y - 0.018, sub, ha="center", va="center", fontsize=6.0, color="#3A3A3A")
        ax.annotate(
            count,
            xy=(x0 + width + 0.014, y),
            ha="left",
            va="center",
            fontsize=7.2,
            color=color,
            fontweight="bold",
        )

    ys = [r[4] for r in rows]
    for i in range(len(ys) - 1):
        col = native if i >= 5 else screen
        ax.add_patch(
            FancyArrowPatch(
                (cx, ys[i] - 0.036),
                (cx, ys[i + 1] + 0.036),
                arrowstyle="-|>",
                mutation_scale=8,
                linewidth=0.9,
                color=col,
            )
        )

    # Bracket marking the region evaluated in a single pass over the library.
    bx = 0.935
    ax.plot([bx, bx + 0.022, bx + 0.022, bx], [0.891, 0.891, 0.619, 0.619],
            color="#7F7F7F", linewidth=0.7)
    ax.text(bx + 0.030, 0.755, "single pass\nover the library", rotation=90, ha="left",
            va="center", fontsize=6.2, color="#4D4D4D")

    fs.save(fig, OUT, "fig1_cascade")
    print("wrote fig1_cascade")


# --------------------------------------------------------------------------
# Figure 2 - receptor-side pharmacophore
# --------------------------------------------------------------------------
def fig2_pharmacophore() -> None:
    """Receptor-side GLP-1R pharmacophore.

    Replaces the previous 3D scatter, which relied on a red/green legend and
    buried its residue labels behind overlapping markers. Here the 102 features
    are projected onto their own two principal axes, families are encoded by
    Okabe-Ito color *and* marker shape, and only the 21 curated contact
    features carry residue labels.
    """
    data = json.loads((ROOT / "maps" / "pharmacophore_rigorous.json").read_text())
    feats = pd.DataFrame(data["features"])

    xyz = feats[["x", "y", "z"]].to_numpy()
    centered = xyz - xyz.mean(axis=0)
    # Principal axes: the site is elongated along the peptide-binding axis, so
    # PC1/PC2 give the most informative flat view.
    _, _, vt = np.linalg.svd(centered, full_matrices=False)
    proj = centered @ vt[:2].T

    family_style = {
        "hydrophobic": (fs.OKABE_ITO["grey"], "o", "Hydrophobe"),
        "aromatic": (fs.OKABE_ITO["green"], "h", "Aromatic ring"),
        "negative": (fs.OKABE_ITO["vermillion"], "v", "Anion"),
        "positive": (fs.OKABE_ITO["blue"], "^", "Cation"),
        "hbd": (fs.OKABE_ITO["sky"], "D", "H-bond donor"),
        "hba": (fs.OKABE_ITO["orange"], "s", "H-bond acceptor"),
    }
    group_edge = {
        "ECD anchoring": fs.OKABE_ITO["purple"],
        "Upper TMD activation pocket": fs.OKABE_ITO["black"],
        "ECL1 support": fs.OKABE_ITO["yellow"],
    }

    fig, (ax, axb) = plt.subplots(
        1, 2, figsize=(fs.PAGE_WIDTH, 3.1), gridspec_kw={"width_ratios": [2.15, 1]}
    )

    for fam, (color, marker, _label) in family_style.items():
        m = feats["type"] == fam
        if not m.any():
            continue
        ax.scatter(
            proj[m.to_numpy(), 0],
            proj[m.to_numpy(), 1],
            s=8 + 5.5 * feats.loc[m, "weight"],
            c=color,
            marker=marker,
            linewidths=0.4,
            edgecolors="white",
            alpha=0.9,
            zorder=2,
        )

    # Place residue labels using their rendered extents, so labels do not
    # collide with each other or spill outside the plotting area.
    curated = feats[feats["curated_group"].notna()]
    seen = set()
    labels = []
    for idx, r in curated.iterrows():
        p = proj[feats.index.get_loc(idx)]
        ax.scatter(p[0],p[1],s=26+5.5*r["weight"],facecolors="none",
                   edgecolors=group_edge[r["curated_group"]],linewidths=1.0,zorder=3)
        tag = f"{r['resname']}{int(r['resnum'])}"
        if tag in seen: continue
        seen.add(tag)
        labels.append(ax.annotate(tag,xy=p,xytext=(7,7),textcoords='offset points',
            fontsize=6.0,color='#1A1A1A',zorder=4,
            bbox=dict(facecolor='white',edgecolor='none',alpha=.90,pad=.25),
            arrowprops=dict(arrowstyle='-',color='#8A8A8A',lw=.35)))

    ax.set_xlabel("Principal axis 1 (Å)")
    ax.set_ylabel("Principal axis 2 (Å)")
    ax.set_aspect("equal", adjustable="datalim")
    fs.panel_label(ax, "(a)", dx=-0.11)

    handles = [
        Line2D([], [], marker=mk, color="none", markerfacecolor=c, markeredgecolor="white",
               markersize=4.5, label=lab)
        for c, mk, lab in family_style.values()
    ] + [
        Line2D([], [], marker="o", color="none", markerfacecolor="none", markeredgecolor=c,
               markersize=6, label=g)
        for g, c in group_edge.items()
    ]
    # Legend below the axes so it never sits on top of the feature cloud.
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.19),
              ncol=3, fontsize=5.9, handletextpad=0.3, columnspacing=0.9)

    # Panel b: where the model's weight actually sits.
    feats["group"] = feats["curated_group"].fillna("Uncurated interface")
    agg = feats.groupby("group").agg(n=("weight", "size"), w=("weight", "sum"))
    order = ["ECD anchoring", "Upper TMD activation pocket", "ECL1 support", "Uncurated interface"]
    agg = agg.reindex([o for o in order if o in agg.index])
    colors = [group_edge.get(g, fs.OKABE_ITO["grey"]) for g in agg.index]
    ypos = np.arange(len(agg))[::-1]
    axb.barh(ypos, agg["w"], color=colors, edgecolor="white", linewidth=0.5, height=0.62)
    for y, (n, w) in zip(ypos, agg[["n", "w"]].to_numpy()):
        axb.text(w + 3, y, f"{int(n)} feat.", va="center", fontsize=6.2, color="#3A3A3A")
    axb.set_yticks(ypos)
    axb.set_yticklabels([g.replace("Upper TMD activation pocket", "Upper TMD\nactivation pocket")
                         .replace("Uncurated interface", "Uncurated\ninterface") for g in agg.index],
                        fontsize=6.4)
    axb.set_xlabel("Summed feature weight")
    axb.set_xlim(0, float(agg["w"].max()) * 1.32)
    axb.grid(axis="y", visible=False)
    fs.panel_label(axb, "(b)", dx=-0.42)

    fig.tight_layout(w_pad=1.6)
    from matplotlib.text import Text
    from matplotlib.transforms import Bbox
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    placed=[Bbox.from_bounds(x-3,y-3,6,6) for x,y in ax.transData.transform(proj)]
    bounds=ax.get_window_extent(renderer)
    offsets=[(dx,dy) for dy in [7,-7,17,-17,27,-27,37,-37] for dx in [7,-7,20,-20,35,-35]]
    for label in labels:
        choices=[]
        for dx,dy in offsets:
            label.set_position((dx,dy));label.set_ha('left' if dx>0 else 'right');label.set_va('bottom' if dy>0 else 'top')
            box=Text.get_window_extent(label,renderer).expanded(1.12,1.18)
            overlap=sum(max(0,min(box.x1,b.x1)-max(box.x0,b.x0))*max(0,min(box.y1,b.y1)-max(box.y0,b.y0)) for b in placed)
            outside=max(0,bounds.x0-box.x0)+max(0,box.x1-bounds.x1)+max(0,bounds.y0-box.y0)+max(0,box.y1-bounds.y1)
            choices.append((overlap*100+outside*1000+dx*dx+dy*dy,dx,dy,box))
        _,dx,dy,box=min(choices,key=lambda item:item[0])
        label.set_position((dx,dy));label.set_ha('left' if dx>0 else 'right');label.set_va('bottom' if dy>0 else 'top')
        placed.append(box)
    fs.save(fig, OUT, "fig2_pharmacophore")
    print("wrote fig2_pharmacophore")


# --------------------------------------------------------------------------
# Figure 3 - GLP-1R benchmark, four methods
# --------------------------------------------------------------------------
def _summary(path: Path) -> pd.DataFrame:
    return pd.read_csv(path).set_index("method")


def fig3_glp1r_benchmark() -> None:
    """GLP-1R enrichment for all four methods, including native-only.

    Native-only was previously absent from this figure even though it is the
    comparison the paper turns on, so it is shown here alongside the cascade.
    """
    summ = _summary(EV / "benchmark_glp1r_full" / "benchmark_summary.csv")
    methods = ["full_cascade", "native_only", "stage3_only", "standard_3d_pharmacophore"]
    metrics = [
        ("roc_auc", "ROC-AUC", None),
        ("pr_auc", "PR-AUC", None),
        ("ef_1pct", "EF1%", None),
        ("bedroc", "BEDROC ($\\alpha=20$)", None),
    ]

    fig, axes = plt.subplots(1, 4, figsize=(fs.PAGE_WIDTH, 2.15))
    x = np.arange(len(methods))
    for ax, (key, title, _), letter in zip(axes, metrics, "abcd"):
        vals = [summ.loc[m, key] for m in methods]
        lo = [summ.loc[m, key] - summ.loc[m, f"{key}_ci_low"] for m in methods]
        hi = [summ.loc[m, f"{key}_ci_high"] - summ.loc[m, key] for m in methods]
        for i, m in enumerate(methods):
            st = fs.METHOD_STYLE[m]
            ax.bar(x[i], vals[i], width=0.68, color=st["color"], hatch=st["hatch"],
                   edgecolor="white", linewidth=0.6, zorder=2)
        ax.errorbar(x, vals, yerr=[lo, hi], fmt="none", ecolor="#333333",
                    elinewidth=0.7, capsize=2, capthick=0.7, zorder=3)
        ax.set_title(title)
        ax.set_xticks(x)
        ax.set_xticklabels([fs.METHOD_STYLE[m]["label"] for m in methods],
                           rotation=38, ha="right", fontsize=6.4)
        ax.grid(axis="x", visible=False)
        fs.panel_label(ax, f"({letter})", dx=-0.30, dy=1.13)

    fig.tight_layout(w_pad=1.1)
    fs.save(fig, OUT, "fig3_glp1r_benchmark")
    print("wrote fig3_glp1r_benchmark")


# --------------------------------------------------------------------------
# Figure 4 - cross-system enrichment
# --------------------------------------------------------------------------
def fig4_cross_system() -> None:
    """Cascade, native-only, and the conventional baseline on every system.

    The panel makes the paper's central negative result visible: on MDM2-p53
    the native-only baseline is far ahead of the full cascade.
    """
    paths = {
        "GLP-1R": EV / "benchmark_glp1r_full" / "benchmark_summary.csv",
        "GHSR": EV / "benchmark_ghsr_full" / "benchmark_summary.csv",
        "NTSR1": EV / "benchmark_ntsr1_full" / "benchmark_summary.csv",
        "MDM2-p53": EV / "benchmark_mdm2_full" / "benchmark_summary.csv",
    }
    summaries = {k: _summary(p) for k, p in paths.items() if p.exists()}
    methods = ["full_cascade", "native_only", "standard_3d_pharmacophore"]
    metrics = [("roc_auc", "ROC-AUC"), ("ef_1pct", "EF1%"), ("bedroc", "BEDROC ($\\alpha=20$)")]

    fig, axes = plt.subplots(1, 3, figsize=(fs.PAGE_WIDTH, 2.35))
    systems = [s for s in SYSTEM_ORDER if s in summaries]
    x = np.arange(len(systems))
    width = 0.26

    for ax, (key, title), letter in zip(axes, metrics, "abc"):
        for j, m in enumerate(methods):
            st = fs.METHOD_STYLE[m]
            vals, lo, hi = [], [], []
            for s in systems:
                sm = summaries[s]
                v = sm.loc[m, key] if m in sm.index else np.nan
                vals.append(v)
                lo.append(v - sm.loc[m, f"{key}_ci_low"] if m in sm.index else 0)
                hi.append(sm.loc[m, f"{key}_ci_high"] - v if m in sm.index else 0)
            pos = x + (j - 1) * width
            ax.bar(pos, vals, width=width, color=st["color"], hatch=st["hatch"],
                   edgecolor="white", linewidth=0.5, label=st["label"], zorder=2)
            ax.errorbar(pos, vals, yerr=[lo, hi], fmt="none", ecolor="#333333",
                        elinewidth=0.6, capsize=1.5, capthick=0.6, zorder=3)
            # A zero-height bar is indistinguishable from a missing series, so
            # say so in print: EF1% is genuinely 0 for MDM2-p53 single-pass 3D.
            for p, v in zip(pos, vals):
                if v == 0:
                    ax.annotate("0", xy=(p, 0), xytext=(0, 2.5),
                                textcoords="offset points", ha="center", va="bottom",
                                fontsize=6.0, color=st["color"], fontweight="bold",
                                zorder=4)
        ax.set_title(title)
        ax.set_xticks(x)
        ax.set_xticklabels([shown(s) for s in systems], rotation=28, ha="right",
                           fontsize=6.6)
        ax.grid(axis="x", visible=False)
        fs.panel_label(ax, f"({letter})", dx=-0.22, dy=1.13)

    # Interpretation is written from the completed comparisons in the text;
    # an annotation tied to the old strict-filter result would be misleading.
    axes[0].set_ylim(0, 1.05)
    handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor=fs.METHOD_STYLE[m]["color"],
                      hatch=fs.METHOD_STYLE[m]["hatch"], edgecolor="white",
                      linewidth=0.5, label=fs.METHOD_STYLE[m]["label"])
        for m in methods
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=6.6,
               bbox_to_anchor=(0.5, -0.055))
    fig.tight_layout(w_pad=1.3)
    fs.save(fig, OUT, "fig4_cross_system")
    print("wrote fig4_cross_system")


# --------------------------------------------------------------------------
# Figure 5 - efficiency-retention trade-off
# --------------------------------------------------------------------------
def fig6_efficiency() -> None:
    """Actives retained against native-scoring calls as the shortlist varies.

    Reading: a favorable trade is a curve that reaches its plateau far to the
    left of the production marker. GHSR and NTSR1 do not, which is the honest
    cost of a fixed-percentage shortlist on a small candidate pool.
    """
    sweep = pd.read_csv(EV / "efficiency" / "shortlist_sweep.csv")

    fig, (ax, axb) = plt.subplots(1, 2, figsize=(fs.PAGE_WIDTH, 2.5))

    for s in SYSTEM_ORDER:
        sub = sweep[sweep.system == s].sort_values("n_native_scored")
        st = fs.SYSTEM_STYLE[s]
        ax.plot(sub.n_native_scored, sub.active_retention * 100, color=st["color"],
                marker="none", linewidth=1.3, label=s, zorder=2)
        prod = sub.iloc[(sub.shortlist_fraction - 0.05).abs().argmin()]
        ax.plot(prod.n_native_scored, prod.active_retention * 100, marker=st["marker"],
                color=st["color"], markersize=5, markeredgecolor="white",
                markeredgewidth=0.6, zorder=4)

    ax.set_xscale("log")
    ax.set_xlabel("Molecules reaching native scoring (log scale)")
    ax.set_ylabel("In-domain actives retained (%)")
    ax.set_ylim(0, 100)
    ax.legend(loc="lower right", fontsize=6.4)
    ax.text(0.03, 0.96, "markers = production 5% setting", transform=ax.transAxes,
            fontsize=6.2, va="top", color="#4D4D4D")
    fs.panel_label(ax, "(a)", dx=-0.17)

    # Panel b: retention at the production setting vs the permissive setting.
    ypos = np.arange(len(SYSTEM_ORDER))[::-1]
    prod_vals, full_vals = [], []
    for s in SYSTEM_ORDER:
        sub = sweep[sweep.system == s]
        prod_vals.append(sub.iloc[(sub.shortlist_fraction - 0.05).abs().argmin()].active_retention * 100)
        full_vals.append(sub.iloc[sub.shortlist_fraction.argmax()].active_retention * 100)
    axb.barh(ypos + 0.19, full_vals, height=0.34, color=fs.OKABE_ITO["sky"],
             edgecolor="white", linewidth=0.5, label="Permissive (100%)", zorder=2)
    axb.barh(ypos - 0.19, prod_vals, height=0.34, color=fs.OKABE_ITO["blue"],
             hatch="//", edgecolor="white", linewidth=0.5, label="Production (5%)", zorder=2)
    for y, (p, f_) in zip(ypos, zip(prod_vals, full_vals)):
        if f_ - p > 3:
            axb.annotate(f"+{f_ - p:.0f} pts", xy=(f_ + 1.5, y + 0.19), fontsize=6.0,
                         va="center", color="#3A3A3A")
    axb.set_yticks(ypos)
    axb.set_yticklabels(SYSTEM_ORDER, fontsize=6.8)
    axb.set_xlabel("In-domain actives retained (%)")
    axb.set_xlim(0, 118)
    axb.grid(axis="y", visible=False)
    axb.legend(loc="upper center", bbox_to_anchor=(0.5, -0.20), ncol=2, fontsize=6.2)
    fs.panel_label(axb, "(b)", dx=-0.26)

    fig.tight_layout(w_pad=1.6)
    fs.save(fig, OUT, "fig6_efficiency")
    print("wrote fig6_efficiency")


# --------------------------------------------------------------------------
# Figure 6 - active attrition
# --------------------------------------------------------------------------
def fig5_attrition() -> None:
    """Where each system's in-domain actives are lost, and to what cause."""
    att = pd.read_csv(EV / "attrition" / "active_attrition_per_molecule.csv")

    causes = [
        ("structural alert", fs.OKABE_ITO["vermillion"], "//", "Stage 0: structural alert"),
        ("property envelope", fs.OKABE_ITO["orange"], "\\\\", "Stage 0: property envelope"),
        ("hotspot gate", fs.OKABE_ITO["purple"], "xx", "Stage 1: hotspot gate"),
        ("below shortlist cut", fs.OKABE_ITO["sky"], "..", "Stage-3 shortlist cut"),
        ("unparseable", fs.OKABE_ITO["black"], "++", "Stage 0: invalid input"),
        ("3D scoring failure", fs.OKABE_ITO["vermillion"], "oo", "Stage 3: no valid score"),
        ("native pool or scaffold cap", fs.OKABE_ITO["grey"], "//", "Native selection"),
        ("native preparation or scoring failure", fs.OKABE_ITO["yellow"], "xx", "Native preparation/scoring"),
        ("below final rank limit", fs.OKABE_ITO["blue"], "..", "Final top-1,000 cutoff"),
        ("retained", fs.OKABE_ITO["green"], "", "In final ranking"),
    ]
    causes=[c for c in causes if (att.reason_class==c[0]).any()]
    if not set(att.reason_class).issubset({c[0] for c in causes}):
        raise ValueError('Unrepresented attrition cause')

    fig, ax = plt.subplots(figsize=(fs.PAGE_WIDTH, 2.8))
    ypos = np.arange(len(SYSTEM_ORDER))[::-1]

    for y, s in zip(ypos, SYSTEM_ORDER):
        sub = att[att.system == s]
        total = len(sub)
        left = 0.0
        for cause, color, hatch, _label in causes:
            n = int((sub.reason_class == cause).sum())
            if not n:
                continue
            pct = 100 * n / total
            ax.barh(y, pct, left=left, height=0.6, color=color, hatch=hatch,
                    edgecolor="white", linewidth=0.6, zorder=2)
            if pct >= 3:
                ax.text(left + pct / 2, y, str(n), ha="center", va="center",
                        fontsize=7, color="black", fontweight="bold", zorder=3,
                        bbox=dict(facecolor='white',edgecolor='none',alpha=.85,pad=.7))
            left += pct

    ax.set_yticks(ypos)
    ax.set_yticklabels(
        [f"{shown(s)}\n($n$={len(att[att.system == s])})" for s in SYSTEM_ORDER],
        fontsize=6.6,
    )
    ax.set_xlabel("In-domain actives (%)")
    ax.set_xlim(0, 100)
    ax.grid(axis="y", visible=False)
    handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor=c, hatch=h, edgecolor="white", linewidth=0.5, label=lab)
        for _cause, c, h, lab in causes
    ]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.30),
              ncol=2, fontsize=6.0)
    fig.tight_layout()
    fs.save(fig, OUT, "fig5_attrition")
    print("wrote fig5_attrition")


# --------------------------------------------------------------------------
# Figure 7 - orthogonal docking
# --------------------------------------------------------------------------
def fig7_docking() -> None:
    """Redocking of the current top-10 ligands against both receptor states.

    Built from ``docking_top10_summary.csv``, the run that matches the reported
    top-10 ranking. The superseded figure was drawn from an earlier run whose
    ligand set no longer matched the manuscript's table.
    """
    df = pd.read_csv(EV / "revision" / "docking_top10_summary.csv").sort_values("final_rank")
    stats = json.loads((EV / "equivalence" / "docking_wilcoxon.json").read_text())

    fig, (ax, axb) = plt.subplots(
        1, 2, figsize=(fs.PAGE_WIDTH, 2.5), gridspec_kw={"width_ratios": [2.1, 1]}
    )

    ypos = np.arange(len(df))[::-1]
    for y, (_, r) in zip(ypos, df.iterrows()):
        ax.plot([r.best_inactive, r.best_active], [y, y], color="#9A9A9A",
                linewidth=0.9, zorder=1, solid_capstyle="round")
    ax.scatter(df.best_inactive, ypos, s=22, color=fs.OKABE_ITO["orange"], marker="s",
               edgecolors="white", linewidths=0.5, label="Best inactive-state", zorder=3)
    ax.scatter(df.best_active, ypos, s=22, color=fs.OKABE_ITO["blue"], marker="o",
               edgecolors="white", linewidths=0.5, label="Best active-state", zorder=3)

    ax.set_yticks(ypos)
    ax.set_yticklabels([f"{int(r.final_rank)}. {r.ligand_id}" for _, r in df.iterrows()],
                       fontsize=6.0)
    ax.set_xlabel("AutoDock Vina affinity (kcal mol$^{-1}$)")
    ax.invert_xaxis()
    ax.grid(axis="y", visible=False)
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.04), ncol=2,
              borderaxespad=0, fontsize=6.2)
    fs.panel_label(ax, "(a)", dx=-0.33)

    # Panel b: the paired differences the significance test is run on.
    d = df["active_pref"].to_numpy()
    axb.axvline(0, color="#9A9A9A", linewidth=0.7, zorder=1)
    axb.scatter(d, ypos, s=22, color=fs.OKABE_ITO["purple"], marker="D",
                edgecolors="white", linewidths=0.5, zorder=3)
    for y,value in zip(ypos,d):
        if not np.isfinite(value):axb.text(0,y,'no pair',ha='center',va='center',fontsize=6,color='#555555')
    axb.axvline(float(stats['median_delta']), color=fs.OKABE_ITO["purple"], linewidth=0.8,
                linestyle="--", zorder=2)
    axb.set_yticks(ypos)
    axb.set_yticklabels([])
    axb.set_xlabel("$\\Delta$ = active $-$ inactive\n(kcal mol$^{-1}$)")
    axb.grid(axis="y", visible=False)
    # Statistics go in the panel title, where nothing can collide with them.
    axb.set_title(
        f"median {stats['median_delta']:.2f}; "
        f"{stats['n_favoring_active']}/{stats['n']} favor active\n"
        f"Wilcoxon $p$ = {stats['p_two_sided']:.4f}",
        fontsize=6.2,
    )
    fs.panel_label(axb, "(b)", dx=-0.18, dy=1.22)

    fig.tight_layout(w_pad=1.2)
    fs.save(fig, OUT, "fig7_docking")
    print("wrote fig7_docking")


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
    fig1_cascade()
    fig2_pharmacophore()
    fig3_glp1r_benchmark()
    fig4_cross_system()
    if args.target == "acs" and (EV / "absolute_floor/policy_comparison.csv").exists():
        from analyze_absolute_floor import figures
        figures(pd.read_csv(EV / "absolute_floor/policy_comparison.csv"),
                pd.read_csv(EV / "absolute_floor/all_stage_survival.csv"))
    else:
        fig6_efficiency()
    fig5_attrition()
    fig7_docking()


if __name__ == "__main__":
    main()
