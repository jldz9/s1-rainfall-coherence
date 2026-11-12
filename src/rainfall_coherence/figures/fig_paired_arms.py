#!/usr/bin/env python
"""The two contrasts taken against the SAME control interferogram.

The null arm was built by staging a second dry pair beside a control that a
rain unit already used.  So for those bursts both contrasts share a denominator:

    rain pair  vs control      and      second dry pair  vs that same control

Nothing differs between them except whether the numerator's interferogram spans
rain — not the burst, not the season, not the control, not the processing run.
Any burst-level property cancels exactly, which no unpaired comparison can
claim.
"""
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
RAIN, DRY, WITHIN = "#9E2A2B", "#335C81", "#1B7A6B"
rng = np.random.default_rng(0)


def cboot(df, col, B=2000):
    """Resample bursts, not units: units within a burst are not independent."""
    bs = df.burst.unique()
    idx = {b: g[col].to_numpy() for b, g in df.groupby("burst")}
    o = [np.median(np.concatenate([idx[b] for b in rng.choice(bs, len(bs), True)]))
         for _ in range(B)]
    return np.percentile(o, [2.5, 97.5])


def by_dose(df, col, edges, labels):
    g = pd.cut(df.rain_mm, edges, labels=labels, right=False)
    rows = []
    for k in labels:
        s = df[g == k]
        lo, hi = cboot(s, col)
        rows.append((s[col].median(), lo, hi, len(s)))
    return pd.DataFrame(rows, columns=["m", "lo", "hi", "n"])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args()
    out = a.out or Path("figs")
    out.mkdir(parents=True, exist_ok=True)
    j = io.shared_control(a.root)
    w = io.within_burst(a.root)      # already carries loss_wet, loss_dry, diff

    fig, ax = plt.subplots(1, 3, figsize=(16.5, 4.9), constrained_layout=True)

    # ---- (a) the two contrasts as distributions ------------------------
    A = ax[0]
    bins = np.linspace(-60, 90, 62)
    A.hist(j.red_dry, bins=bins, density=True, color=DRY, alpha=0.40, lw=0,
           label=f"second dry pair (n = {len(j):,})", zorder=2)
    A.hist(j.red_rain, bins=bins, density=True, histtype="step", color=RAIN,
           lw=1.8, label=f"rain pair (n = {len(j):,})", zorder=3)
    for v, c in ((j.red_dry.median(), DRY), (j.red_rain.median(), RAIN)):
        A.axvline(v, color=c, lw=1.6, zorder=4)
    A.annotate(f"{j.red_dry.median():+.1f}%", (j.red_dry.median(), 0.56),
               xycoords=("data", "axes fraction"), color=DRY, fontsize=10.5,
               ha="right", xytext=(-3, 0), textcoords="offset points")
    A.annotate(f"{j.red_rain.median():+.1f}%", (j.red_rain.median(), 0.68),
               xycoords=("data", "axes fraction"), color=RAIN, fontsize=10.5,
               ha="left", xytext=(3, 0), textcoords="offset points")
    A.set_xlabel("coherence reduction against the shared control (%)", fontsize=11)
    A.set_ylabel("probability density", fontsize=11)
    A.set_title("(a)  both measured against the same control", fontsize=12.5,
                loc="left")
    A.tick_params(labelsize=10)
    A.legend(loc="upper right", fontsize=10)

    # ---- (b) the paired difference, by dose ----------------------------
    B = ax[1]
    j["diff"] = j.red_rain - j.red_dry
    lab = ["< 5", "5–15", "15–30", "≥ 30"]
    y = np.arange(len(lab))
    t = by_dose(j, "diff", [0, 5, 15, 30, 1e9], lab)
    B.axvline(0, color="0.65", lw=1.1, zorder=1)
    B.hlines(y, t.lo, t.hi, color=RAIN, lw=2.8, alpha=0.5, zorder=2)
    B.plot(t.m, y, "o", color=RAIN, ms=8, zorder=3)
    for i, r in enumerate(t.itertuples()):
        B.annotate(f"{r.m:.1f} pp", (0.99, i), xycoords=("axes fraction", "data"),
                   fontsize=11, color="0.3", va="center", ha="right")
    B.set_yticks(y)
    B.set_yticklabels([f"{k} mm\nn = {int(r):,}" for k, r in zip(lab, t.n)],
                      fontsize=10.5)
    B.set_xlabel("extra loss attributable to the rain (pp)", fontsize=11)
    B.set_xlim(0, 46)
    B.set_ylim(-0.6, len(lab) - 0.4)
    B.set_title("(b)  paired difference, by accumulation", fontsize=12.5,
                loc="left")
    B.tick_params(axis="x", labelsize=10)

    # ---- (c) wetted against unwetted ground of one interferogram -------
    C = ax[2]
    labw = ["< 2", "2–5", "5–15", "≥ 15"]
    yw = np.arange(len(labw))
    u = by_dose(w, "diff", [0, 2, 5, 15, 1e9], labw)
    C.axvline(0, color="0.65", lw=1.1, zorder=1)
    C.hlines(yw, u.lo, u.hi, color=WITHIN, lw=2.8, alpha=0.5, zorder=2)
    C.plot(u.m, yw, "o", color=WITHIN, ms=8, zorder=3)
    for i, r in enumerate(u.itertuples()):
        C.annotate(f"{r.m:.1f} pp", (0.99, i), xycoords=("axes fraction", "data"),
                   fontsize=11, color="0.3", va="center", ha="right")
    C.set_yticks(yw)
    C.set_yticklabels([f"{k} mm\nn = {int(r):,}" for k, r in zip(labw, u.n)],
                      fontsize=10.5)
    C.set_xlabel("extra loss over the wetted ground (pp)", fontsize=11)
    C.set_xlim(0, 46)
    C.set_ylim(-0.6, len(labw) - 0.4)
    C.set_title("(c)  within one interferogram, by accumulation",
                fontsize=12.5, loc="left")
    C.tick_params(axis="x", labelsize=10)

    fig.savefig(out / "fig3_paired_arms.png")
    fig.savefig(out / "fig3_paired_arms.pdf")
    plt.close(fig)
    d = j.red_rain - j.red_dry
    print(f"wrote {out}/fig3_paired_arms.png / .pdf")
    print(f"  rain {j.red_rain.median():+.2f}%   dry {j.red_dry.median():+.2f}%   "
          f"difference {d.median():+.2f} pp, rain larger in {100*(d>0).mean():.1f}%")
    print(t.round(2).to_string(index=False))
    print(f"  within-burst: n = {len(w):,} over {w.burst.nunique():,} bursts, "
          f"median {w['diff'].median():+.2f} pp, "
          f"wet lower in {100 * (w['diff'] > 0).mean():.1f}%")
    print(u.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
