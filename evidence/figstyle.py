#!/usr/bin/env python3
"""Shared publication figure style for the IEEE Access submission.

Every figure in the manuscript is rendered through this module so that the
whole set reads as one system: a single colorblind-safe palette, one type
scale, IEEE Access column widths, and 600 dpi raster output alongside vector
PDF.

Design rules enforced here
--------------------------
* **Palette.** Okabe-Ito, the standard eight-color qualitative palette that
  stays distinguishable under deuteranopia, protanopia, and tritanopia
  (Okabe & Ito, "Color Universal Design", 2008). Red/green pairs are never
  used to carry meaning on their own.
* **Redundant encoding.** Series are separated by hatch, marker, and direct
  label in addition to hue, so every figure survives grayscale printing.
* **Sizing.** IEEE Access is a two-column journal: ``COL_WIDTH`` (3.5 in) for
  single-column figures and ``PAGE_WIDTH`` (7.16 in) for double-column ones.
  Type is set so that body labels land near 8 pt at final printed size.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------- dimensions
# IEEE Access text block: 3.5 in single column, 7.16 in across both columns.
COL_WIDTH = 3.5
PAGE_WIDTH = 7.16

#: Text-block widths per journal, as (full width, half width) in inches.
#: ACS Omega here is single-column letter paper with 1 in margins, so the full
#: text block is 6.5 in and a half-width panel is 3.25 in.
TARGETS = {
    "ieee": (7.16, 3.5),
    "acs": (6.5, 3.25),
}


def use_target(name: str) -> None:
    """Set the page and column widths for the journal being built.

    Figure code reads ``fs.PAGE_WIDTH`` and ``fs.COL_WIDTH`` at call time, so
    calling this before rendering is enough to resize the whole set.
    """
    global PAGE_WIDTH, COL_WIDTH
    try:
        PAGE_WIDTH, COL_WIDTH = TARGETS[name]
    except KeyError:
        raise SystemExit(f"unknown target {name!r}; expected one of {sorted(TARGETS)}")

# ------------------------------------------------------------------ palette
# Okabe-Ito. Keys are semantic so call sites never hard-code a hex value.
OKABE_ITO = {
    "black": "#000000",
    "orange": "#E69F00",
    "sky": "#56B4E9",
    "green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "purple": "#CC79A7",
    "grey": "#7F7F7F",
}

#: Ordered palette for categorical series, most-distinct first.
CATEGORICAL = [
    OKABE_ITO["blue"],
    OKABE_ITO["orange"],
    OKABE_ITO["green"],
    OKABE_ITO["purple"],
    OKABE_ITO["sky"],
    OKABE_ITO["vermillion"],
]

#: Fixed color, hatch, and marker per screening method, so a method looks the
#: same in every figure it appears in.
METHOD_STYLE = {
    "full_cascade": {
        "label": "Full cascade",
        "color": OKABE_ITO["blue"],
        "hatch": "",
        "marker": "o",
    },
    "native_only": {
        "label": "Native-only",
        "color": OKABE_ITO["orange"],
        "hatch": "//",
        "marker": "s",
    },
    "stage3_only": {
        "label": "Stage-3 only",
        "color": OKABE_ITO["sky"],
        "hatch": "..",
        "marker": "^",
    },
    "standard_3d_pharmacophore": {
        "label": "Single-pass 3D",
        "color": OKABE_ITO["grey"],
        "hatch": "xx",
        "marker": "D",
    },
}

#: Fixed color and hatch per cascade stage, shared by the funnel, the
#: attrition waterfall, and the retention curves.
STAGE_STYLE = {
    "stage0": {"label": "Stage 0 (property/chemistry)", "color": OKABE_ITO["vermillion"], "hatch": "//"},
    "stage1": {"label": "Stage 1 (hotspot gate)", "color": OKABE_ITO["orange"], "hatch": "\\\\"},
    "stage3": {"label": "Stage-3 shortlist", "color": OKABE_ITO["sky"], "hatch": ".."},
    "native": {"label": "Native scoring", "color": OKABE_ITO["green"], "hatch": ""},
    "retained": {"label": "Retained", "color": OKABE_ITO["blue"], "hatch": ""},
}

#: Fixed color per benchmark system.
SYSTEM_STYLE = {
    "GLP-1R": {"color": OKABE_ITO["blue"], "marker": "o"},
    "GHSR": {"color": OKABE_ITO["orange"], "marker": "s"},
    "NTSR1": {"color": OKABE_ITO["green"], "marker": "^"},
    "MDM2-p53": {"color": OKABE_ITO["purple"], "marker": "D"},
}


def apply() -> None:
    """Install the manuscript-wide matplotlib style."""
    plt.rcParams.update(
        {
            # Times-compatible stack, matching the IEEEtran body font.
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Nimbus Roman", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "font.size": 8,
            "axes.titlesize": 8.5,
            "axes.labelsize": 8,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "legend.fontsize": 7.5,
            "figure.titlesize": 9,
            # Restrained frame: no top/right spines, hairline axes.
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.6,
            "axes.titlepad": 4.0,
            "axes.labelpad": 2.5,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            # Grid sits behind the data and stays faint.
            "axes.grid": True,
            "axes.axisbelow": True,
            "grid.color": "#D9D9D9",
            "grid.linewidth": 0.4,
            "grid.alpha": 0.9,
            "legend.frameon": False,
            "legend.handlelength": 1.4,
            "legend.handletextpad": 0.5,
            "legend.columnspacing": 1.1,
            "lines.linewidth": 1.2,
            "lines.markersize": 3.5,
            "patch.linewidth": 0.5,
            "hatch.linewidth": 0.5,
            "figure.dpi": 150,
            "savefig.dpi": 600,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.02,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def save(fig, out_dir: Path, stem: str) -> None:
    """Write ``stem`` as both a 600 dpi PNG and a vector PDF.

    IEEE Access accepts either; the PDF is preferred for line art and the PNG
    documents the 600 dpi raster requirement.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / f"{stem}.pdf")
    fig.savefig(out_dir / f"{stem}.png", dpi=600)
    plt.close(fig)


def panel_label(ax, letter: str, dx: float = -0.16, dy: float = 1.06) -> None:
    """Place a bold (a)/(b)/(c) panel label in axes coordinates."""
    ax.text(
        dx,
        dy,
        letter,
        transform=ax.transAxes,
        fontsize=9,
        fontweight="bold",
        va="top",
        ha="left",
    )
