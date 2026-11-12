#!/usr/bin/env python
"""How the land surface conditions the coherence lost to rain.

Land cover and climate are confounded in any global sample: shrubland is mostly
arid, moss is mostly polar. Comparing cover classes across bursts therefore
measures climate as much as it measures the surface. The coherence here is
resolved per cover class within each burst, so the comparison can instead be
made inside a single burst, where the storm, the climate, the incidence angle
and the processing are all held fixed by construction.

Both estimates are reported: the naive between-burst one, which is what a
cover-only stratification would give, and the within-burst one, which is the
effect of the surface alone.
"""
import argparse, numpy as np, pandas as pd
from scipy import stats
from pathlib import Path

MIN_PX = 20_000
FLOOR = 0.077
MIN_HEADROOM = 0.15
REF = "grass"               # most evenly spread across climates
KEEP = ["tree", "grass", "crops", "bare", "shrub", "built", "moss"]


def load(root: Path) -> pd.DataFrame:
    v = pd.read_parquet(root / "plan" / "analysis_v4.parquet")
    v = v[v.event.str.startswith("rain")].copy()
    v = v[(v.n_px >= MIN_PX) & (v.g_ctrl - FLOOR >= MIN_HEADROOM)]
    v = v[v.cover.isin(KEEP)].copy()
    v["red"] = (1 - v.g_case / v.g_ctrl) * 100
    v["rain"] = np.log2(1 + v.case_24h.clip(lower=0))
    # a burst can contribute more than one row per cover class; collapse them
    # into a single pixel-weighted value so each burst-cover pair appears once
    def collapse(g):
        w = g.n_px.to_numpy(float)
        return pd.Series({"red": float(np.average(g.red, weights=w)),
                          "rain": float(g.rain.iloc[0]),
                          "koppen": g.koppen.iloc[0],
                          "n_px": float(w.sum())})
    v = (v.groupby(["unit_id", "cover"], as_index=False)
           .apply(collapse, include_groups=False).reset_index(drop=True))
    return v


def naive(v: pd.DataFrame) -> None:
    print("Between-burst comparison (confounded with climate)")
    print(f"  {'cover':<10}{'n':>6}{'reduction':>12}{'se':>7}")
    for c, g in sorted(v.groupby("cover"), key=lambda kv: -kv[1].red.mean()):
        se = g.red.std() / np.sqrt(len(g))
        print(f"  {c:<10}{len(g):6d}{g.red.mean():11.1f}%{se:7.1f}")


def paired(v: pd.DataFrame) -> pd.DataFrame:
    """Difference against the reference cover inside the same burst."""
    ref = v[v.cover == REF].set_index("unit_id")["red"]
    rows = []
    for c in KEEP:
        if c == REF:
            continue
        g = v[v.cover == c].set_index("unit_id")
        both = g.index.intersection(ref.index)
        if len(both) < 20:
            continue
        diff = g.loc[both, "red"] - ref.loc[both]
        t, p = stats.ttest_1samp(diff, 0)
        rows.append({"cover": c, "n_pairs": len(both), "diff": diff.mean(),
                     "se": diff.std() / np.sqrt(len(diff)), "p": p,
                     "rain": g.loc[both, "rain"].to_numpy(),
                     "d": diff.to_numpy()})
    return pd.DataFrame(rows).sort_values("diff", ascending=False)


def fixed_effects(v: pd.DataFrame) -> None:
    """Cover effect with a free intercept per burst, SEs clustered by burst."""
    v = v[v.groupby("unit_id").unit_id.transform("size") >= 2].copy()
    y = v.red.to_numpy(float)
    lab = v.cover.astype(str)
    levels = [c for c in KEEP if c != REF and (lab == c).sum() >= 20]
    D = np.column_stack([(lab == c).to_numpy(float) for c in levels])
    # within-burst demeaning removes the burst intercepts
    g = v.unit_id.to_numpy()
    def demean(M):
        M = np.atleast_2d(M.T).T.astype(float)
        out = M.copy()
        df = pd.DataFrame(M); df["g"] = g
        return (df.groupby("g").transform(lambda s: s - s.mean())).to_numpy()
    Xw, yw = demean(D), demean(y).ravel()
    b, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
    resid = yw - Xw @ b
    XtXi = np.linalg.pinv(Xw.T @ Xw)
    meat = np.zeros_like(XtXi)
    for u in np.unique(g):
        m = g == u
        s = Xw[m].T @ resid[m]
        meat += np.outer(s, s)
    V = XtXi @ meat @ XtXi
    se = np.sqrt(np.diag(V))
    nclu = len(np.unique(g))
    print(f"\nWithin-burst comparison, {nclu} bursts, reference = {REF}")
    print(f"  {'cover':<10}{'extra loss':>13}{'se':>7}{'p':>10}")
    order = np.argsort(-b)
    for i in order:
        t = b[i] / se[i]
        p = 2 * stats.t.sf(abs(t), nclu - 1)
        star = "***" if p < 1e-3 else "**" if p < 0.01 else "*" if p < 0.05 else ""
        print(f"  {levels[i]:<10}{b[i]:+12.1f}%{se[i]:7.1f}{p:10.2e} {star}")


def by_koppen(v: pd.DataFrame) -> None:
    print("\nIs the surface effect the same in every climate?")
    ref = v[v.cover == REF].set_index("unit_id")["red"]
    kop = v.drop_duplicates("unit_id").set_index("unit_id")["koppen"]
    for c in KEEP:
        if c == REF:
            continue
        g = v[v.cover == c].set_index("unit_id")
        both = g.index.intersection(ref.index)
        if len(both) < 60:
            continue
        diff = (g.loc[both, "red"] - ref.loc[both]).to_frame("d")
        diff["k"] = kop.loc[both].to_numpy()
        parts = [(k, s.d) for k, s in diff.groupby("k") if len(s) >= 15]
        if len(parts) < 2:
            continue
        F, p = stats.f_oneway(*[s for _, s in parts])
        txt = "  ".join(f"{k}:{s.mean():+5.1f}(n={len(s)})" for k, s in parts)
        star = "***" if p < 1e-3 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
        print(f"  {c:<8s} vs {REF}:  {txt}   between-climate p={p:.3f} {star}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    a = ap.parse_args()
    v = load(a.root)
    print(f"{len(v)} cover-rows over {v.unit_id.nunique()} rain bursts "
          f"(n_px >= {MIN_PX:,}, headroom >= {MIN_HEADROOM})\n")
    naive(v)
    pr = paired(v)
    print(f"\nPaired against {REF} inside the same burst")
    print(f"  {'cover':<10}{'pairs':>7}{'extra loss':>13}{'se':>7}{'p':>10}")
    for _, r in pr.iterrows():
        star = "***" if r.p < 1e-3 else "**" if r.p < 0.01 else "*" if r.p < 0.05 else ""
        print(f"  {r.cover:<10}{r.n_pairs:7d}{r['diff']:+12.1f}%{r.se:7.1f}{r.p:10.2e} {star}")
    fixed_effects(v)
    by_koppen(v)

    print("\nDoes the surface effect grow with the size of the storm?")
    for _, r in pr.iterrows():
        rho, p = stats.spearmanr(r.rain, r.d)
        print(f"  {r.cover:<10} rho={rho:+.3f}  p={p:.3f}")


if __name__ == "__main__":
    main()
