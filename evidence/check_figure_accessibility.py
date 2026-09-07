#!/usr/bin/env python3
"""Verify the submission figures meet the resolution and color requirements.

Three checks:

1. **Resolution.** Every raster figure must declare 600 dpi.
2. **Palette.** Every color a figure paints in quantity must be close to a
   color in the Okabe-Ito colorblind-safe palette (or to a neutral). This
   catches a stray matplotlib default slipping into an otherwise safe figure.
3. **Grayscale legibility.** Converting to luminance must preserve contrast:
   if two series are distinguished by hue alone, the grayscale version
   collapses. We approximate this by requiring that the figure's luminance
   histogram keeps a wide spread, which a hue-only encoding does not.

The palette check simulates the three common forms of color-vision
deficiency and reports the minimum separation between the palette colors
that survive under each, which is the quantity that actually matters.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parent.parent
SUB = ROOT / "IEEE_Access_submission"

OKABE_ITO = {
    "black": (0, 0, 0), "orange": (230, 159, 0), "sky": (86, 180, 233),
    "green": (0, 158, 115), "yellow": (240, 228, 66), "blue": (0, 114, 178),
    "vermillion": (213, 94, 0), "purple": (204, 121, 167), "grey": (127, 127, 127),
}

# Brettel/Vienot-style linear approximations of dichromatic vision.
CVD = {
    "protanopia": np.array([[0.152, 1.053, -0.205], [0.115, 0.786, 0.099],
                            [-0.004, -0.048, 1.052]]),
    "deuteranopia": np.array([[0.367, 0.861, -0.228], [0.280, 0.673, 0.047],
                              [-0.012, 0.043, 0.969]]),
    "tritanopia": np.array([[1.256, -0.077, -0.179], [-0.078, 0.931, 0.148],
                            [0.005, 0.691, 0.304]]),
}


def simulate(rgb: np.ndarray, kind: str) -> np.ndarray:
    return np.clip(rgb @ CVD[kind].T, 0, 255)


def palette_separation() -> None:
    """Minimum pairwise distance among palette colors under each CVD type."""
    names = [k for k in OKABE_ITO if k != "black"]
    base = np.array([OKABE_ITO[n] for n in names], dtype=float)
    print("Okabe-Ito separation (min pairwise RGB distance):")
    for kind in ["normal", *CVD]:
        pts = base if kind == "normal" else simulate(base, kind)
        d = np.linalg.norm(pts[:, None, :] - pts[None, :, :], axis=-1)
        np.fill_diagonal(d, np.inf)
        i, j = np.unravel_index(np.argmin(d), d.shape)
        print(f"  {kind:14s} min = {d.min():6.1f}  ({names[i]} vs {names[j]})")


#: Figures that legitimately paint colors outside the categorical palette,
#: with the reason. Both remain colorblind-safe.
PALETTE_EXEMPT = {
    "figA1_upstream_vs_native.png":
        "sequential Blues density map - single hue varying in lightness only, "
        "which is safe under every form of color-vision deficiency",
    "figA4_top20_structures.png":
        "RDKit structure drawings retain explicit element labels; the current "
        "ACS exports use a monochrome atom palette",
}


def check_figure(path: Path) -> dict:
    im = Image.open(path)
    dpi = im.info.get("dpi", (0, 0))[0]
    rgb = np.asarray(im.convert("RGB"), dtype=float)
    flat = rgb.reshape(-1, 3)

    # Ignore near-white background when judging color usage.
    ink = flat[flat.sum(axis=1) < 720]
    if len(ink) == 0:
        return {"file": path.name, "dpi": dpi, "off_palette": 0.0, "grey_spread": 0.0}

    palette = np.array(list(OKABE_ITO.values()), dtype=float)
    sample = ink[np.random.default_rng(0).choice(len(ink), min(60000, len(ink)),
                                                 replace=False)]
    # Neutrals (any grey) are always acceptable.
    is_neutral = np.ptp(sample, axis=1) < 40
    d = np.linalg.norm(sample[:, None, :] - palette[None, :, :], axis=-1).min(axis=1)
    off = float((~is_neutral & (d > 90)).mean())

    lum = sample @ np.array([0.299, 0.587, 0.114])
    spread = float(np.percentile(lum, 95) - np.percentile(lum, 5))
    return {"file": path.name, "dpi": dpi, "off_palette": off, "grey_spread": spread}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dir", default=str(SUB))
    args = ap.parse_args()

    palette_separation()
    print("\nPer-figure checks (off-palette = fraction of ink far from the palette;")
    print("grey spread = luminance range surviving grayscale conversion):")
    print(f"  {'file':<34s} {'dpi':>6s} {'off-palette':>12s} {'grey spread':>12s}  status")

    ok = True
    paths = sorted(Path(args.dir).glob("fig[0-9A]*.png"))
    toc = Path(args.dir) / "toc_graphic.png"
    if toc.exists():
        paths.append(toc)
    if not paths:
        raise SystemExit("No publication figures found")
    for p in paths:
        r = check_figure(p)
        flags = []
        if abs(r["dpi"] - 600) > 1:
            flags.append("DPI")
        if r["off_palette"] > 0.05 and r["file"] not in PALETTE_EXEMPT:
            flags.append("PALETTE")
        if r["grey_spread"] < 40:
            flags.append("GREYSCALE")
        exempt = r["file"] in PALETTE_EXEMPT and r["off_palette"] > 0.05
        status = "CHECK: " + ",".join(flags) if flags else ("ok (exempt)" if exempt else "ok")
        ok &= not flags
        print(f"  {r['file']:<34s} {r['dpi']:>6.0f} {r['off_palette']:>11.1%} "
              f"{r['grey_spread']:>12.1f}  {status}")

    if PALETTE_EXEMPT:
        print("\nDocumented palette exemptions:")
        for name, why in PALETTE_EXEMPT.items():
            print(f"  {name}: {why}")

    print("\nAll figures pass." if ok else "\nSome figures need review.")
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
