#!/usr/bin/env python
"""Perpendicular baseline: real, monotonic, and unable to bias the result.

The design fixes the temporal baseline at twelve days but leaves the orbital
separation free, so a case pair and its control differ in perpendicular
baseline.  That difference does move the coherence.  What protects the
measurement is not that the effect is negligible but that it is symmetric:
the case is no more likely than the control to carry the wider orbit, so the
loss it imparts cancels in the mean.
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
RAIN, DRY = "#9E2A2B", "#335C81"
rng = np.random.default_rng(0)


def cboot(df, col, B=2000):
    """Resample bursts, not units."""
    bs = df.burst.unique()
    idx = {b: g[col].to_numpy() for b, g in df.groupby("burst")}
    o = [np.median(np.concatenate([idx[b] for b in rng.choice(bs, len(bs), True)]))
         for _ in range(B)]
    return np.percentile(o, [2.5, 97.5])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args()
    out = a.out or Path("figs")
    out.mkdir(parents=True, exist_ok=True)

    b = io.baselines(a.root).dropna(subset=["bperp"])
    u = io.units(a.root)
    p = b.pivot_table(index="unit_id", columns="role", values="bperp", aggfunc="first")
    m = (u.merge(p[["case", "control"]], left_on="unit_id", right_index=True,
                 how="inner").dropna(subset=["case", "control"]))
    m["dbp"] = m["case"].abs() - m["control"].abs()

    fig, ax = plt.subplots(1, 2, figsize=(11.0, 4.3), constrained_layout=True)

    # ---- (a) the difference is centred on zero -------------------------
    q25, q75 = m.dbp.quantile([0.25, 0.75])
    A = ax[0]
    A.hist(m.dbp.clip(-250, 250), bins=np.linspace(-250, 250, 66), color=DRY,
           alpha=0.45, lw=0, zorder=2)
    A.axvspan(q25, q75, color=DRY, alpha=0.13, lw=0, zorder=1)
    A.axvline(0, color="0.5", lw=1.0, zorder=3)
    A.axvline(m.dbp.mean(), color=RAIN, lw=1.8, zorder=4)
    A.annotate(f"mean {m.dbp.mean():+.1f} m,  sd {m.dbp.std():.0f} m\n"
               f"half the units lie within\n{q25:.0f} to {q75:+.0f} m (shaded)",
               (0.03, 0.95), xycoords="axes fraction", fontsize=10.5,
               va="top", color="0.25")
    A.set_xlabel("case minus control in |perpendicular baseline| (m)", fontsize=11)
    A.set_ylabel(f"units (n = {len(m):,})", fontsize=11)
    A.set_title("(a)  the orbital separation is no wider for the case",
                fontsize=12.5, loc="left")
    A.tick_params(labelsize=10)

    # ---- (b) it moves the loss, but hardly ever by much -----------------
    B = ax[1]
    q = pd.qcut(m.dbp, 10, labels=False)
    rows = []
    for k, g in m.groupby(q):
        lo, hi = cboot(g, "reduction_pct")
        rows.append((g.dbp.median(), g.reduction_pct.median(), lo, hi, len(g)))
    t = pd.DataFrame(rows, columns=["x", "y", "lo", "hi", "n"])
    lvl = m.loc[m.dbp.abs() <= 20, "reduction_pct"].median()
    fit = np.polyfit(m.dbp, m.reduction_pct, 1)

    B.axvspan(q25, q75, color=DRY, alpha=0.13, lw=0, zorder=1)
    B.axhline(lvl, color=RAIN, lw=1.3, ls="--", zorder=2)
    B.annotate(f"level set by the rainfall, {lvl:.1f}%", (0.985, lvl),
               xycoords=("axes fraction", "data"), fontsize=10, color=RAIN,
               ha="right", va="top", xytext=(0, -5), textcoords="offset points")
    B.errorbar(t.x, t.y, yerr=[t.y - t.lo, t.hi - t.y], fmt="o", color=DRY,
               ms=6.5, lw=1.4, capsize=2.5, zorder=4)
    B.plot(t.x, t.y, "-", color=DRY, lw=1.1, zorder=3)
    B.annotate(f"{fit[0] * 100:.2f} pp per 100 m, but\n"
               f"only {fit[0] * (q75 - q25):.1f} pp across the shaded half\n"
               f"and {fit[0] * m.dbp.mean():.2f} pp of bias at the mean",
               (0.03, 0.97), xycoords="axes fraction", fontsize=10.5,
               va="top", color="0.25")
    B.set_xlabel("case minus control in |perpendicular baseline| (m)", fontsize=11)
    B.set_ylabel("total coherence reduction (%)", fontsize=11)
    B.set_ylim(0, 38)
    B.set_title("(b)  median total reduction by decile of that difference",
                fontsize=12.5, loc="left")
    B.tick_params(labelsize=10)

    fig.savefig(out / "fig5_baseline.png")
    fig.savefig(out / "fig5_baseline.pdf")
    plt.close(fig)
    print(f"wrote {out}/fig5_baseline.png / .pdf")
    print(f"  n = {len(m):,} units over {m.burst.nunique():,} bursts")
    print(f"  dbp mean {m.dbp.mean():+.2f} m, median {m.dbp.median():+.1f} m, "
          f"sd {m.dbp.std():.1f} m")
    print(f"  slope {fit[0] * 100:.2f} pp per 100 m -> "
          f"{fit[0] * m.dbp.mean():.3f} pp of bias at the mean")
    print(t.round(2).to_string(index=False))
    print(f"  central half {q25:.0f} to {q75:+.0f} m -> "
          f"{fit[0] * (q75 - q25):.2f} pp across it; rainfall level {lvl:.2f}%")


if __name__ == "__main__":
    main()
