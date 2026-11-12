#!/usr/bin/env python
"""The dose response, with the logarithmic fit drawn over it.

The horizontal axis is logarithmic, so a logarithmic response is a straight
line: the fit is there to let a reader check the form of the relationship
rather than take the variance comparison on trust.  It is a median regression,
because the plotted points are medians and the distribution is skewed enough
that a least-squares line would sit below every one of them.
"""
import argparse
from pathlib import Path
import numpy as np, pandas as pd
import statsmodels.api as sm
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .. import io

plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans",
                     "axes.linewidth": 0.7, "axes.spines.top": False,
                     "axes.spines.right": False, "savefig.dpi": 420,
                     "savefig.bbox": "tight", "legend.frameon": False})
RAIN, GREY = "#9E2A2B", "0.45"
EDGES = [1, 3, 6, 10, 15, 20, 25, 30, 40, 50, 70, 100, 320]
rng = np.random.default_rng(0)


def cboot(df, col="reduction_pct", B=2000):
    bs = df.burst.unique()
    idx = {b: g[col].to_numpy() for b, g in df.groupby("burst")}
    return np.percentile(
        [np.median(np.concatenate([idx[b] for b in rng.choice(bs, len(bs), True)]))
         for _ in range(B)], [2.5, 97.5])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--linear", action="store_true",
                    help="plot the accumulation on a linear axis, so that the "
                         "logarithmic fit shows as the concave curve it is")
    a = ap.parse_args()
    out = a.out or Path("figs")
    out.mkdir(parents=True, exist_ok=True)

    u = io.units(a.root)
    v = io.validation(a.root)

    rows = []
    for lo, hi in zip(EDGES[:-1], EDGES[1:]):
        g = u[(u.rain_mm >= lo) & (u.rain_mm < hi)]
        if len(g) < 40:
            continue
        l, h = cboot(g)
        rows.append((g.rain_mm.median(), g.reduction_pct.median(), l, h, len(g), lo))
    b = pd.DataFrame(rows, columns=["x", "y", "lo", "hi", "n", "edge"])

    # The measured zero: both members dry, the same arithmetic applied.
    # Its own generator, so the interval does not depend on how many bands
    # were bootstrapped before it, and 50,000 draws, since 2,000 leave some
    # 0.1 of Monte Carlo noise upon the bounds.  The arm holds one unit per
    # burst, so resampling units and resampling bursts are the same thing.
    nul = v.red.median()
    nrng = np.random.default_rng(0)
    r = v.red.to_numpy()
    nlo, nhi = np.percentile(
        [np.median(nrng.choice(r, r.size, True)) for _ in range(50_000)],
        [2.5, 97.5])

    # median regression, since the points are medians
    X = sm.add_constant(np.log2(u.rain_mm.to_numpy()))
    q = sm.QuantReg(u.reduction_pct.to_numpy(), X).fit(q=0.5)
    c0, slope = q.params
    ols = sm.OLS(u.reduction_pct.to_numpy(), X).fit()

    fig, ax = plt.subplots(figsize=(5.6, 4.6), constrained_layout=True)
    ax.axhspan(nlo, nhi, color="0.85", lw=0, zorder=1)
    ax.axhline(nul, color="0.35", ls="--", lw=1.1, zorder=2)
    ax.annotate(f"dry against dry, {nul:+.2f}% [{nlo:+.2f}, {nhi:+.2f}]",
                (0.97, nul), xycoords=("axes fraction", "data"), fontsize=8.5,
                ha="right", va="bottom", color="0.35",
                xytext=(0, 4), textcoords="offset points")

    xs = np.linspace(np.log2(u.rain_mm.min()), np.log2(u.rain_mm.max()), 100)
    ax.plot(2 ** xs, c0 + slope * xs, "-", color=RAIN, lw=1.5, alpha=0.7,
            zorder=3, label=f"median fit, {slope:.2f} pp per doubling")

    light = b.edge < 15
    ax.errorbar(b.x[light], b.y[light], yerr=[(b.y - b.lo)[light], (b.hi - b.y)[light]],
                fmt="o", mfc="white", color=RAIN, ms=6, lw=1.3, capsize=2.5,
                zorder=5, label="sampled (< 15 mm)")
    ax.errorbar(b.x[~light], b.y[~light], yerr=[(b.y - b.lo)[~light], (b.hi - b.y)[~light]],
                fmt="o", color=RAIN, ms=6, lw=1.3, capsize=2.5, zorder=5,
                label="every event (≥ 15 mm)")

    for r in b.itertuples():
        ax.annotate(f"{r.n:,}", (r.x, -0.075), xycoords=("data", "axes fraction"),
                    fontsize=7, rotation=90, ha="center", va="top", color=GREY)
    ax.annotate("n", (0.0, -0.075), xycoords="axes fraction", fontsize=7,
                ha="right", va="top", color=GREY)

    if not a.linear:
        ax.set_xscale("log")
        ax.set_xticks([1, 2, 5, 10, 20, 50, 100, 200])
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.get_xaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())
    else:
        ax.set_xlim(0, 135)
    ax.set_ylim(-3, 1.08 * b.hi.max())      # the rejected fit leaves the frame
    ax.set_xlabel("24 h accumulation at the event acquisition (mm)", fontsize=10,
                  labelpad=34)
    ax.set_ylabel("coherence reduction (%)", fontsize=10)
    ax.tick_params(labelsize=9)
    ax.legend(loc="upper left", fontsize=8.5)

    stem = "fig_dose_linear" if a.linear else "fig6_dose"
    fig.savefig(out / f"{stem}.png")
    fig.savefig(out / f"{stem}.pdf")
    plt.close(fig)
    print(f"wrote {out}/{stem}.png / .pdf")
    print(f"  median fit {slope:.3f} pp/doubling (intercept {c0:.2f}); "
          f"OLS {ols.params[1]:.3f}, R2 {ols.rsquared:.4f}")
    print(f"  validation arm {nul:+.2f}% [{nlo:+.2f}, {nhi:+.2f}], n = {len(v):,}")
    print(b.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
