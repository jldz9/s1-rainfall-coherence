#!/usr/bin/env python
"""Event and control coherence considered separately, against rainfall.

The ratio reported throughout could in principle move because the event pair
falls or because the control rises.  Plotting the two separately settles which:
the control is nearly flat across the whole range of accumulation while the
event pair descends, so the trend resides in the numerator.

Bands, statistic and interval follow Table 2 -- medians, with a cluster
bootstrap over bursts -- so that this figure and that table are the same
arithmetic applied to the two coherences rather than to their ratio.

Usage:
    python fig_r12_dose.py --root <wx root>
"""
from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .. import io

plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans",
                     "axes.linewidth": 0.7, "axes.spines.top": False,
                     "axes.spines.right": False, "savefig.dpi": 420,
                     "savefig.bbox": "tight", "legend.frameon": False})
RAIN, CTRL = "#9E2A2B", "#335C81"
FLOOR = 0.077
EDGES = [1, 3, 6, 10, 15, 20, 25, 30, 40, 50, 70, 100, 320]   # Table 2
rng = np.random.default_rng(0)


def cboot(d: pd.DataFrame, col: str, B: int = 2000):
    """Resample bursts, not units: units within a burst are not independent."""
    idx = {b: g[col].to_numpy() for b, g in d.groupby("burst")}
    bs = list(idx)
    o = [np.median(np.concatenate([idx[b] for b in rng.choice(bs, len(bs), True)]))
         for _ in range(B)]
    return np.percentile(o, [2.5, 97.5])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args()
    f = io.units(a.root)
    outdir = a.out or Path("figs")
    outdir.mkdir(parents=True, exist_ok=True)

    rows = []
    for lo, hi in zip(EDGES[:-1], EDGES[1:]):
        s = f[(f.rain_mm >= lo) & (f.rain_mm < hi)]
        if len(s) < 20:
            continue
        clo, chi = cboot(s, "g_ctrl")
        elo, ehi = cboot(s, "g_case")
        rows.append(dict(x=np.sqrt(lo * hi), n=len(s),
                         c=s.g_ctrl.median(), clo=clo, chi=chi,
                         e=s.g_case.median(), elo=elo, ehi=ehi))
    t = pd.DataFrame(rows)

    fig, ax = plt.subplots(figsize=(5.2, 4.3), constrained_layout=True)
    ax.errorbar(t.x, t.c, yerr=[t.c - t.clo, t.chi - t.c], fmt="s-", color=CTRL,
                ms=4, lw=1.1, capsize=2, label="control pair")
    ax.errorbar(t.x, t.e, yerr=[t.e - t.elo, t.ehi - t.e], fmt="o-", color=RAIN,
                ms=4, lw=1.1, capsize=2, label="event pair")
    ax.axhline(FLOOR, ls=":", color="0.55", lw=0.9)
    ax.text(t.x.iloc[-1], FLOOR + 0.008, f"noise floor {FLOOR}", color="0.45",
            ha="right", va="bottom", fontsize=7.5)

    # labels go in clear space, not beside the last point where the error bars sit
    rc, re = t.c.max() - t.c.min(), t.e.max() - t.e.min()
    ax.text(0.97, 0.88, f"control varies by {rc:.3f}", transform=ax.transAxes,
            ha="right", va="center", color=CTRL, fontsize=7.5)
    ax.text(0.97, 0.26, f"event varies by {re:.3f}", transform=ax.transAxes,
            ha="right", va="center", color=RAIN, fontsize=7.5)

    ax.set_xscale("log")
    ax.set_xlabel("24 h accumulation at the event acquisition (mm)")
    ax.set_ylabel("coherence")
    ax.set_ylim(0, 0.62)
    ax.legend(loc="lower left")
    out = outdir / "fig7_case_control.png"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".pdf"))
    print(f"control {t.c.iloc[0]:.3f} -> {t.c.iloc[-1]:.3f}  range {rc:.3f}")
    print(f"event   {t.e.iloc[0]:.3f} -> {t.e.iloc[-1]:.3f}  range {re:.3f}")
    print(t[["x", "n", "c", "e"]].round(3).to_string(index=False))
    print("wrote", out)


if __name__ == "__main__":
    main()
