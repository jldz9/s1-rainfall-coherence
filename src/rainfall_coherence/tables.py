"""The tables of the paper, each reproduced from the deposited data."""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import io
from .config import DOSE_EDGES
from .stats import cluster_bootstrap, log_fit


def dose(root=None) -> pd.DataFrame:
    """Table 3.  Coherence reduction by accumulation band.

    Reductions are medians; standard errors resample bursts.
    """
    u = io.units(root)
    rows = []
    for lo, hi in zip(DOSE_EDGES[:-1], DOSE_EDGES[1:]):
        s = u[(u.rain_mm >= lo) & (u.rain_mm < hi)]
        if len(s) < 20:
            continue
        med, blo, bhi = cluster_bootstrap(s.reduction_pct, s.burst)
        rows.append(
            dict(band=f"{lo}-{hi}", units=len(s), bursts=s.burst.nunique(),
                 reduction=round(med, 1), lo=round(blo, 2), hi=round(bhi, 2))
        )
    return pd.DataFrame(rows)


def landcover_contrasts(root=None, reference: str = "grass") -> pd.DataFrame:
    """Table 4.  Extra loss of each cover against grassland, within a unit.

    Collapse to one pixel-weighted value per unit and cover, then pair against
    the reference inside the same unit.

    The interval resamples bursts, as everywhere else in the paper.  Pairing
    removes the burst from the contrast but not the clustering: units of one
    burst still share its terrain, and sd/sqrt(n) over units runs 1.2 to 1.8
    times too narrow depending on the cover.

    No further filter is applied, and none should be.  The headroom condition
    belongs to the unit and was imposed when the analysed sample was built;
    imposing it again here removes the dim covers, moves snow and ice from
    -12.6 to -4.8, and turns the wetland contrast from -1.6 to +1.6.
    """
    d = io.landcover(root)
    g = (
        d.groupby(["unit_id", "cover"])
        .apply(lambda x: np.average(x.reduction_pct, weights=x.n_px),
               include_groups=False)
        .rename("reduction_pct")
        .reset_index()
    )
    burst = io.units(root).set_index("unit_id").burst
    ref = g[g.cover == reference].set_index("unit_id").reduction_pct

    rows = []
    for cover in sorted(g.cover.unique()):
        if cover == reference:
            continue
        s = g[g.cover == cover].set_index("unit_id").reduction_pct
        shared = s.index.intersection(ref.index)
        if len(shared) < 20:
            continue
        diff = s.loc[shared] - ref.loc[shared]
        b = burst.reindex(diff.index)
        mean, lo, hi, draws = cluster_bootstrap(
            diff, b, stat=np.mean, n=4000, return_draws=True)
        rows.append(
            dict(cover=cover, units=len(shared), bursts=int(b.nunique()),
                 extra_loss=round(mean, 1),
                 se=round(float(draws.std(ddof=1)), 2),
                 lo=round(lo, 2), hi=round(hi, 2))
        )
    return pd.DataFrame(rows).sort_values("extra_loss", ascending=False)


def climate_and_season(root=None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Table 5.  Reduction by Koppen class and by hemisphere-corrected season."""
    u = io.units(root)
    k = (u.groupby("koppen")
           .agg(units=("unit_id", "size"),
                reduction=("reduction_pct", "mean"),
                median_rain=("rain_mm", "median"))
           .round(1).reset_index())
    s = (u.groupby("season")
           .agg(units=("unit_id", "size"),
                reduction=("reduction_pct", "mean"),
                median_rain=("rain_mm", "median"))
           .round(1).sort_values("reduction", ascending=False).reset_index())
    return k, s


def dose_relationship(root=None) -> dict:
    """Equation (3), with the comparison of functional forms it rests on."""
    u = io.units(root)
    y, x = u.reduction_pct.values, u.rain_mm.values
    intercept, slope, r2, zero = log_fit(x, y)

    def _r2(xx):
        X = np.column_stack([np.ones(len(xx)), xx])
        b, *_ = np.linalg.lstsq(X, y, rcond=None)
        res = y - X @ b
        return 1 - (res**2).sum() / ((y - y.mean()) ** 2).sum()

    heavy = u[u.rain_mm >= 15]
    return {
        "intercept": round(intercept, 2),
        "slope": round(slope, 2),
        "r2_log": round(r2, 4),
        "r2_linear": round(_r2(x), 4),
        "r2_sqrt": round(_r2(np.sqrt(x)), 4),
        "zero_mm": round(zero, 3),
        "n": len(u),
        "slope_above_15mm": round(log_fit(heavy.rain_mm, heavy.reduction_pct)[1], 2),
        "n_above_15mm": len(heavy),
    }


def measurement_floor(root=None) -> dict:
    """The width of the dry-against-dry distribution, by cover and by climate."""
    v = io.validation(root)
    u = io.units(root)
    cover = (u.assign(b=u.burst.str.lower())
               .groupby("b").dom_cover.agg(lambda s: s.mode().iat[0]))
    v = v.assign(b=v.burst.str.lower())
    v["cover"] = v.b.map(cover)
    v["abs_pc"] = v.red.abs()
    by_cover = v.dropna(subset=["cover"]).groupby("cover").abs_pc.agg(["median", "size"])
    return {
        "median_pc": round(v.abs_pc.median(), 2),
        "sd_pc": round(v.red.std(), 1),
        "n": len(v),
        "by_cover": by_cover[by_cover["size"] >= 100]["median"].round(1).to_dict(),
        "by_koppen": v.groupby("koppen").abs_pc.median().round(1).to_dict(),
    }


def population_estimate(root=None) -> dict:
    """The headline reduction weighted to the census of rainfall events.

    The processed sample holds 58.4% of its units above 15 mm where the census
    holds 12.0%, so the unweighted mean describes the units rather than the
    rainfall.  Weighting each unit by the events its band represents gives the
    loss a rainfall event costs.
    """
    import numpy as np

    from . import io
    from .stats import phase_sigma, phase_to_displacement_mm

    u = io.units(root)
    w = io.population_weights(root).to_numpy()
    y = u.reduction_pct.to_numpy()
    ok = np.isfinite(w)
    w, y, u = w[ok], y[ok], u[ok]

    def wmean(v):
        return float((v * w).sum() / w.sum())

    def wmedian(v):
        o = np.argsort(v)
        c = np.cumsum(w[o]) / w.sum()
        return float(np.interp(0.5, c, v[o]))

    idx = {b: g.index.to_numpy() for b, g in u.reset_index(drop=True).groupby("burst")}
    keys, rng, draws = list(idx), np.random.default_rng(0), []
    for _ in range(2000):
        k = np.concatenate([idx[keys[i]] for i in rng.integers(0, len(keys), len(keys))])
        draws.append((y[k] * w[k]).sum() / w[k].sum())

    d = np.rad2deg(phase_sigma(u.g_case.to_numpy()) - phase_sigma(u.g_ctrl.to_numpy()))
    return {
        "weighted_mean_pct": round(wmean(y), 2),
        "ci": (round(float(np.percentile(draws, 2.5)), 2),
               round(float(np.percentile(draws, 97.5)), 2)),
        "unweighted_mean_pct": round(float(y.mean()), 2),
        "below_control_pct": round(wmean((y > 0).astype(float)) * 100, 1),
        "weighted_median_pct": round(wmedian(y), 2),
        "median_accumulation_mm": round(wmedian(u.rain_mm.to_numpy()), 1),
        "mean_accumulation_mm": round(wmean(u.rain_mm.to_numpy()), 1),
        "phase_sigma_deg": round(wmean(d), 1),
        "displacement_mm": round(float(phase_to_displacement_mm(np.deg2rad(wmean(d)))), 2),
        "effective_n": int(w.sum() ** 2 / (w ** 2).sum()),
    }
