#!/usr/bin/env python
"""Figure 4.  The measurement floor is two-sided; the rainfall loss is not.

Two dry interferograms of one burst differ by a median 12.4% of the coherence
of the one nominated as control, which is of the same size as the loss a modest
rainfall event causes.  Magnitude alone therefore cannot separate the two.
Direction can: the dry contrast straddles the one-to-one line and the rainfall
contrast sits below it.

    (a) two DRY interferograms of one burst, plotted against each other
    (b) a rainfall pair against its dry control

Both panels carry identical axes so that the displacement in (b) is read
against the scatter in (a) rather than against itself.

Reads the deposited tables, so `--root` is the deposit directory (or omit it
and set RCOH_DATA).

Usage:
    python -m rainfall_coherence.figures.fig_one_sidedness --root $RCOH_DATA
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .. import io

plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans",
                     "axes.linewidth": 0.7, "axes.spines.top": False,
                     "axes.spines.right": False, "savefig.dpi": 400,
                     "savefig.bbox": "tight", "legend.frameon": False})

RAIN, DRY = "#9E2A2B", "#335C81"
FLOOR = 0.077          # coherence of an uncorrelated burst mean
HEADROOM = FLOOR + 0.15  # the condition imposed on every control


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args()

    v = io.validation(a.root)
    u = io.units(a.root)

    # one of the two dry interferograms is nominated control at random; the
    # other takes the place the event pair occupies in (b), so that the two
    # panels are read on the same construction.
    v = v.assign(g_other=np.where(np.isclose(v.ctrl, v.g_a), v.g_b, v.g_a))

    fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.8), constrained_layout=True)

    # ---- (a) dry against dry -------------------------------------------
    A = ax[0]
    A.scatter(v.ctrl, v.g_other, s=5, color=DRY, alpha=0.35, lw=0,
              rasterized=True)
    A.set_xlabel("coherence of the nominated dry control")
    A.set_ylabel("coherence of the second dry pair")
    below_a = (v.red > 0).mean() * 100
    A.set_title("(a)  two dry interferograms of one burst\n"
                f"n = {len(v):,} units,  {below_a:.1f}% below the line,  "
                f"median {v.red.median():+.2f}%",
                fontsize=8, loc="left")

    # ---- (b) rainfall against its control ------------------------------
    B = ax[1]
    B.scatter(u.g_ctrl, u.g_case, s=2, color=RAIN, alpha=0.20, lw=0,
              rasterized=True)
    B.set_xlabel("coherence of the dry control")
    B.set_ylabel("coherence of the rainfall pair")
    red = (1 - u.g_case / u.g_ctrl) * 100
    below = (u.g_case < u.g_ctrl).mean() * 100
    B.set_title("(b)  a rainfall pair against its dry control\n"
                f"n = {len(u):,} units,  {below:.1f}% below the line,  "
                f"median {red.median():+.1f}%",
                fontsize=8, loc="left")

    # the nominated control must stand 0.15 clear of the floor; its partner
    # is left free, so each cloud has a wall on x and none on y
    for X in (A, B):
        X.plot([0, 1], [0, 1], color="0.3", lw=0.9, zorder=3)
        X.axvline(HEADROOM, color="0.45", lw=0.8, ls="--")
        X.axhline(FLOOR, ls=":", color="0.55", lw=0.8)
        X.set_xlim(0, 1)
        X.set_ylim(0, 1)
        X.set_aspect("equal")
    A.text(HEADROOM + 0.015, 0.97, "headroom condition\non the control",
           color="0.4", fontsize=6.5, va="top")
    B.text(0.98, FLOOR + 0.02, "noise floor 0.077", color="0.45",
           ha="right", fontsize=6.5)

    outdir = a.out or Path("figs")
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "fig4_one_sidedness"
    fig.savefig(out.with_suffix(".png"))
    fig.savefig(out.with_suffix(".pdf"))
    print(f"(a) n={len(v):,}  below the line {below_a:.1f}%  "
          f"median {v.red.median():+.2f}%")
    print(f"(b) n={len(u):,}  below the line {below:.1f}%  "
          f"median {red.median():+.2f}%")
    print("wrote", out.with_suffix(".png"))


if __name__ == "__main__":
    main()
