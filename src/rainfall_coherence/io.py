"""Loading the deposited tables.

Point `RCOH_DATA` at the directory holding the Zenodo deposit, or pass `root`
explicitly.  Every loader returns a plain DataFrame with the derived columns
the analysis uses, so a caller never recomputes a reduction by hand.
"""
from __future__ import annotations

import os
from pathlib import Path

import pandas as pd

from .config import NOISE_FLOOR
from .stats import reduction

FILES = {
    "manifest": "00_unit_manifest.csv",
    "census": "09_event_census.parquet",
    "units": "01_units_analysed.parquet",
    "landcover": "02_units_by_landcover.parquet",
    "validation": "03_validation_arm.parquet",
    "shared": "04_shared_control_arm.parquet",
    "masked": "05_rain_masked.parquet",
    "baselines": "06_perpendicular_baselines.parquet",
    "withinburst": "07_within_burst_wet_dry.parquet",
    "soilwater": "08_soil_water.parquet",
}


def data_root(root: str | Path | None = None) -> Path:
    """Where the deposited tables live.

    An explicit `root` wins, then `RCOH_DATA`, then the `data/` directory that
    ships with the repository.  A clone therefore needs no configuration; set
    RCOH_DATA only to point at a copy held elsewhere.
    """
    if root is not None:
        return Path(root)
    env = os.environ.get("RCOH_DATA")
    if env:
        return Path(env)
    for parent in Path(__file__).resolve().parents:
        bundled = parent / "data"
        if (bundled / FILES["units"]).exists():
            return bundled
    raise RuntimeError(
        "No data directory. Pass root=, set RCOH_DATA, or run from a clone "
        "that still has its data/ directory."
    )


def _read(name: str, root=None) -> pd.DataFrame:
    p = data_root(root) / FILES[name]
    if not p.exists():
        raise FileNotFoundError(f"{p} not found; see README for the deposit.")
    return pd.read_csv(p) if p.suffix == ".csv" else pd.read_parquet(p)


def units(root=None) -> pd.DataFrame:
    """The 21,010 analysed units, with `reduction_pct` and `headroom`."""
    d = _read("units", root)
    if "reduction_pct" not in d:
        d["reduction_pct"] = reduction(d.g_case, d.g_ctrl)
    if "headroom" not in d:
        d["headroom"] = d.g_ctrl - NOISE_FLOOR
    return d


def landcover(root=None, analysed_only: bool = True) -> pd.DataFrame:
    """Coherence resolved by cover class within each unit.

    `analysed_only` restricts to units in the analysed sample, which is the
    only filter this table should receive.  Do not filter it on headroom: that
    condition belongs to the unit and was applied when `units` was built.
    """
    d = _read("landcover", root)
    d = d[(d.g_case > 0) & (d.g_ctrl > 0)].copy()
    if analysed_only:
        d = d[d.unit_id.isin(set(units(root).unit_id))]
    d["reduction_pct"] = reduction(d.g_case, d.g_ctrl)
    return d


def validation(root=None) -> pd.DataFrame:
    """Dry-against-dry pairs; `red` is the signed contrast, `ctrl` its denominator."""
    return _read("validation", root)


def shared_control(root=None) -> pd.DataFrame:
    """Rainfall and dry pairs measured against one shared interferogram."""
    return _read("shared", root)


def masked(root=None, analysed_only: bool = True) -> pd.DataFrame:
    """Coherence and accumulation over the wetted pixels, and over all of them."""
    d = _read("masked", root)
    if analysed_only:
        d = d[d.unit_id.isin(set(units(root).unit_id))]
    # a unit without a wetted population cannot be compared both ways
    d = d.dropna(subset=["g_case_wet", "g_ctrl_wet"]).copy()
    d["loss_all"] = reduction(d.g_case_all, d.g_ctrl_all)
    d["loss_wet"] = reduction(d.g_case_wet, d.g_ctrl_wet)
    return d


def within_burst(root=None) -> pd.DataFrame:
    """Wetted against unwetted ground inside one interferogram.

    `diff` is the within-unit difference of the two losses, which is the
    quantity of interest; the gap between their medians is not.
    """
    d = _read("withinburst", root).dropna(
        subset=["g_case_wet", "g_ctrl_wet", "g_case_dry", "g_ctrl_dry"]
    ).copy()
    d["loss_wet"] = reduction(d.g_case_wet, d.g_ctrl_wet)
    d["loss_dry"] = reduction(d.g_case_dry, d.g_ctrl_dry)
    d["diff"] = d.loss_wet - d.loss_dry
    return d


def soil_water(root=None) -> pd.DataFrame:
    """ERA5-Land theta at the four acquisitions, with the change oriented.

    `dtheta` is theta at the wet acquisition minus theta at the dry one.  The
    orientation matters: taking `sec - ref` without it inverts the sign for
    roughly half the sample and cancels the signal.
    """
    import numpy as np

    t = _read("soilwater", root).dropna().drop_duplicates("unit_id")
    # the unit table carries ctrl_ref/ctrl_sec as DATES; rename to avoid a clash
    t = t.rename(columns={"ctrl_ref": "th_ctrl_ref", "ctrl_sec": "th_ctrl_sec"})
    j = units(root).merge(t, on="unit_id", how="inner")
    event_first = (j.t_event == j.t_ref).to_numpy()
    j["dtheta"] = np.where(
        event_first, j.case_ref - j.case_sec, j.case_sec - j.case_ref
    )
    j["dtheta_ctrl"] = j.th_ctrl_sec - j.th_ctrl_ref
    j["theta_dry"] = np.where(event_first, j.case_sec, j.case_ref)
    return j


def baselines(root=None) -> pd.DataFrame:
    """Perpendicular baseline per interferogram, case and control by `role`."""
    return _read("baselines", root)


def census(root=None) -> pd.DataFrame:
    """The 166,117 rainfall events the units were drawn from.

    Needed for the population weights: the sampling is stratified, every event
    above 15 mm having been processed where its acquisitions allowed while
    lighter events were sampled, so the processed units over-represent heavy
    rain roughly fivefold.
    """
    return _read("census", root)


def population_weights(root=None, edges=None) -> pd.Series:
    """Weight per unit: census events in its band over units processed in it.

    Multiply a per-unit quantity by these before averaging across
    accumulations, and the result describes a rainfall event rather than the
    sample that was processed.  The weights sum to the census total.

    Do NOT apply them to a quantity already conditioned on the accumulation --
    a band of Table 3, the fit of Eq. (3), a cover or climate contrast, or
    either validation arm.  Those are measured at their own dose and weighting
    only adds noise.
    """
    import pandas as pd

    from .config import DOSE_EDGES

    e = DOSE_EDGES if edges is None else edges
    e = list(e[:-1]) + [float("inf")]
    u = units(root)
    c = census(root)
    cb = pd.cut(c.mm, e, right=False)
    ub = pd.cut(u.rain_mm, e, right=False)
    w = c.groupby(cb, observed=True).size() / u.groupby(ub, observed=True).size()
    return ub.map(w).astype(float)
