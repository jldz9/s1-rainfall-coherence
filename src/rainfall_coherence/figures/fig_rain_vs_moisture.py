#!/usr/bin/env python
"""Rainfall, soil water, and what the control pair absorbs.

(a) and (b) put the two predictors of the loss on an equal footing.  (c) is the
check that the control earns its place: control pairs carry soil-water changes
nearly as large as the cases do, those changes depress the control's own
coherence, and yet they leave the ratio untouched, because they move numerator
and denominator together.
"""
import argparse
from pathlib import Path
import glob
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .. import io

plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans",
                     "axes.linewidth": 0.7, "axes.spines.top": False,
                     "axes.spines.right": False, "savefig.dpi": 420,
                     "savefig.bbox": "tight", "legend.frameon": False})
CASE, CTRL, THIRD = "#9E2A2B", "#335C81", "#1B7A6B"
rng = np.random.default_rng(0)


def cboot(df, col, B=1500):
    """Resample bursts, not units."""
    bs = df.burst.unique()
    idx = {b: g[col].to_numpy() for b, g in df.groupby("burst")}
    return np.percentile(
        [np.median(np.concatenate([idx[b] for b in rng.choice(bs, len(bs), True)]))
         for _ in range(B)], [2.5, 97.5])


def binned(df, var, ycol, inner, nmin=40):
    """Open-ended bins: every unit is plotted, and thin bins merge outwards."""
    edges = [-np.inf] + list(inner) + [np.inf]
    while True:
        cnt = pd.cut(df[var], edges).value_counts(sort=False).to_numpy()
        small = np.where(cnt < nmin)[0]
        if not len(small) or len(edges) <= 3:
            break
        k = small[0]
        edges.pop(1 if k == 0 else (len(edges) - 2 if k == len(cnt) - 1 else k))
    q = pd.cut(df[var], edges)
    rows = []
    for _, g in df.groupby(q, observed=True):
        lo, hi = cboot(g, ycol)
        rows.append((g[var].median(), g[ycol].median(), lo, hi, len(g)))
    b = pd.DataFrame(rows, columns=["x", "y", "lo", "hi", "n"])
    b.attrs["kept"] = int(b.n.sum()); b.attrs["total"] = len(df)
    return b

RAIN_E = [3, 6, 10, 15, 20, 25, 30, 40, 50, 70, 100]
DTH_E = [-0.10, -0.05, -0.02, 0, 0.02, 0.05, 0.08, 0.12, 0.18, 0.25]
ADJ_E = [-0.20, -0.12, -0.06, -0.02, 0.02, 0.06, 0.12, 0.20, 0.30]


def panel(a, b, colour, marker, xlab, ylab, title):
    a.errorbar(b.x, b.y, yerr=[b.y - b.lo, b.hi - b.y], fmt=marker, color=colour,
               ms=6, lw=1.3, capsize=2.5, zorder=3)
    a.plot(b.x, b.y, "-", color=colour, lw=1.1, zorder=2)
    a.set_xlabel(xlab, fontsize=11)
    a.set_ylabel(ylab, fontsize=11)
    a.set_title(title, fontsize=12.5, loc="left")
    a.tick_params(labelsize=10)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args()
    out = a.out or Path("figs")
    out.mkdir(parents=True, exist_ok=True)

    # io.soil_water orients the change to the event; see its docstring
    j = io.soil_water(a.root)
    j["dth_case"] = j.dtheta
    m = j                      # one sample for all three panels
    j["dth_ctrl"] = j.dtheta_ctrl
    j["adth_ctrl"] = j.dth_ctrl.abs()

    j["dth_adj"] = j.dth_case - j.dth_ctrl
    fig, axg = plt.subplots(2, 2, figsize=(11.0, 8.4), constrained_layout=True,
                            sharey=True)
    ax = axg.ravel()

    specs = [("rain_mm", RAIN_E, CASE, "o", "24 h accumulation (mm)",
              "(a)  the rainfall"),
             ("dth_case", DTH_E, CTRL, "s",
              "$\\Delta\\theta$ across the case pair (m$^3$ m$^{-3}$)",
              "(b)  soil water, case pair"),
             ("dth_ctrl", DTH_E, THIRD, "D",
              "$\\Delta\\theta$ across the control pair (m$^3$ m$^{-3}$)",
              "(c)  soil water, control pair"),
             ("dth_adj", ADJ_E, "#7A5195", "^",
              "$\\Delta\\theta_{case} - \\Delta\\theta_{control}$ (m$^3$ m$^{-3}$)",
              "(d)  the difference of the two")]

    tabs = []
    for a, (var, edges, colour, mk, xlab, title) in zip(ax, specs):
        b = binned(j, var, "reduction_pct", edges)
        panel(a, b, colour, mk, xlab, "", title)
        if var != "rain_mm":
            a.axvline(0, color="0.8", lw=1.0, zorder=1)
        a.annotate(f"span {b.y.max() - b.y.min():.1f} pp\n"
                   f"{b.attrs['kept']:,} of {b.attrs['total']:,} units",
                   (0.04, 0.95), xycoords="axes fraction", fontsize=10.5,
                   va="top", color="0.3")
        tabs.append((title[1], b))
    for a in (ax[0], ax[2]):
        a.set_ylabel("coherence reduction (%)", fontsize=11)
    ax[0].set_xscale("log")
    ax[0].set_xticks([1, 2, 5, 10, 20, 50, 100, 200])
    ax[0].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax[0].get_xaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())

    fig.savefig(out / "fig8_rain_vs_moisture.png")
    fig.savefig(out / "fig8_rain_vs_moisture.pdf")
    plt.close(fig)
    print(f"wrote {out}/fig8_rain_vs_moisture.png / .pdf")
    for nm, b in tabs:
        print(f"  ({nm}) span {b.y.max() - b.y.min():5.1f} pp   "
              f"plotted {b.attrs['kept']:,} of {b.attrs['total']:,}")


if __name__ == "__main__":
    main()
