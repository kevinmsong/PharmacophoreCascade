#!/usr/bin/env python3
"""Regenerate the ACS Omega graphical abstract / Table-of-Contents graphic
(toc_graphic.png and toc_graphic.pdf), sized to the ACS TOC box (3.25 in wide
x 1.75 in tall).

The layout itself lives in make_floor_toc_graphic; this module is the entry
point the build scripts call.
"""
from make_floor_toc_graphic import main


if __name__ == "__main__":
    main()
