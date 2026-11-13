"""Statistics shared by every table and figure.

The one thing that matters here is that units are not independent.  They are
clustered within bursts, a mean of nine to each and as many as forty-six, and
units of one burst share its terrain and its scattering.  Treating them as
independent understates the uncertainty by a factor of 1.7, so every interval
in the paper resamples bursts rather than units.
"""
from __future__ import annotations

from typing import Callable

import numpy as np
import pandas as pd


def reduction(g_case, g_ctrl):
    """The percentage by which a case pair falls below its control.

    Reported as a ratio rather than a difference because coherence magnitude is
    specific to the processing chain.  A ratio survives a change of scale
    between chains; an absolute difference in coherence does not.
    """
    return (1.0 - np.asarray(g_case) / np.asarray(g_ctrl)) * 100.0


def cluster_bootstrap(
    values: pd.Series | np.ndarray,
    groups: pd.Series | np.ndarray,
    stat: Callable[[np.ndarray], float] = np.median,
    n: int = 2000,
    seed: int = 0,
    interval: tuple[float, float] = (2.5, 97.5),
    return_draws: bool = False,
):
    """Resample whole bursts, recompute `stat`, return (point, low, high).

    With `return_draws` the bootstrap distribution is appended as a fourth
    element, so a caller wanting a standard error can take its spread rather
    than assume the units independent.

    Resampling units instead would understate the spread by roughly 1.7.
    """
    v = np.asarray(values, dtype=float)
    g = np.asarray(groups)
    ok = np.isfinite(v)
    v, g = v[ok], g[ok]
    if v.size == 0:
        empty = (np.nan, np.nan, np.nan)
        return (*empty, np.array([])) if return_draws else empty

    by_group: dict = {}
    for key, idx in pd.Series(np.arange(v.size)).groupby(g):
        by_group[key] = idx.to_numpy()
    keys = list(by_group)

    rng = np.random.default_rng(seed)
    draws = np.empty(n)
    for i in range(n):
        pick = rng.integers(0, len(keys), len(keys))
        draws[i] = stat(v[np.concatenate([by_group[keys[k]] for k in pick])])
    lo, hi = np.percentile(draws, interval)
    point = float(stat(v))
    if return_draws:
        return point, float(lo), float(hi), draws
    return point, float(lo), float(hi)


def paired_difference(a, b, groups, **kw):
    """Difference taken WITHIN each unit, then summarised.

    Not the gap between two summaries.  A median is not additive, so the median
    of the differences is not the difference of the medians, and only the
    former respects the pairing.
    """
    d = np.asarray(a, dtype=float) - np.asarray(b, dtype=float)
    return cluster_bootstrap(d, groups, **kw)


def phase_sigma(gamma, looks: int | None = None) -> np.ndarray:
    """Phase standard deviation in radians, from the Cramer-Rao bound.

    Reported alongside coherence because it maps onto displacement uncertainty
    and is unbounded where gamma is floor-censored.
    """
    from .config import LOOKS_DISPLACEMENT

    L = LOOKS_DISPLACEMENT if looks is None else looks
    g = np.asarray(gamma, dtype=float)
    return np.sqrt(1.0 - g**2) / (g * np.sqrt(2.0 * L))


def phase_to_displacement_mm(sigma_rad) -> np.ndarray:
    """Line-of-sight displacement noise, in millimetres."""
    from .config import WAVELENGTH_MM

    return np.asarray(sigma_rad, dtype=float) * WAVELENGTH_MM / (4.0 * np.pi)


def log_fit(x, y):
    """Least squares of y on log2(x).  Returns (intercept, slope, r2, zero_mm).

    The dose relationship of Eq. (3).  `zero_mm` is the accumulation at which
    the fitted line crosses zero loss.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y) & (x > 0)
    lx, ly = np.log2(x[ok]), y[ok]
    slope, intercept = np.polyfit(lx, ly, 1)
    pred = intercept + slope * lx
    r2 = 1.0 - ((ly - pred) ** 2).sum() / ((ly - ly.mean()) ** 2).sum()
    return float(intercept), float(slope), float(r2), float(2 ** (-intercept / slope))
