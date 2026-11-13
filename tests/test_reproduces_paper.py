"""Every published number this package can reach, checked against the paper.

Needs RCOH_DATA pointing at the deposit; skipped otherwise.
"""
import os
import pytest

pytest.importorskip("pandas")
from rainfall_coherence import io, tables, stats  # noqa: E402

pytestmark = pytest.mark.skipif(
    not os.environ.get("RCOH_DATA"), reason="RCOH_DATA not set"
)


def test_sample_size():
    u = io.units()
    assert len(u) == 21010
    assert u.burst.nunique() == 2329


def test_sample_mean_is_unweighted():
    """The mean over the processed units; the paper's headline is weighted."""
    assert round(io.units().reduction_pct.mean(), 2) == 19.21


def test_equation_3():
    f = tables.dose_relationship()
    assert (f["intercept"], f["slope"]) == (0.82, 5.27)
    assert round(f["r2_log"] * 100, 2) == 10.32
    assert round(f["r2_linear"] * 100, 2) == 7.64
    assert f["slope_above_15mm"] == 6.46


def test_table4_landcover():
    t = tables.landcover_contrasts().set_index("cover")
    assert t.loc["tree", "units"] == 20291
    assert t.loc["tree", "bursts"] == 2128
    assert t.loc["tree", "extra_loss"] == 3.4
    assert t.loc["built", "extra_loss"] == -5.4
    # the sign that a stray headroom filter reverses
    assert t.loc["wetland", "extra_loss"] < 0
    # the interval must resample bursts.  sd/sqrt(n) over units gives 0.11
    # here, and every interval in the table would be half again too tight.
    assert t.loc["tree", "se"] == 0.19
    assert t.loc["snowice", "se"] == 1.25


def test_measurement_floor():
    f = tables.measurement_floor()
    assert f["n"] == 1907
    assert round(f["median_pc"], 1) == 12.4
    assert f["by_cover"]["bare"] < f["by_cover"]["tree"]


def test_within_burst():
    w = io.within_burst()
    assert len(w) == 4167
    assert round(w["diff"].median(), 2) == 1.64


def test_masking_is_negligible():
    m = io.masked()
    assert len(m) == 20944
    assert round(m.loss_wet.mean() - m.loss_all.mean(), 2) == 0.13


def test_soil_water_is_oriented():
    s = io.soil_water()
    from scipy.stats import spearmanr
    rho = spearmanr(s.dtheta, s.reduction_pct).statistic
    assert rho > 0.25, "unoriented delta-theta collapses this to about -0.09"


def test_cluster_bootstrap_is_wider_than_naive():
    u = io.units()
    _, lo, hi = stats.cluster_bootstrap(u.reduction_pct, u.burst, stat=lambda v: v.mean())
    _, nlo, nhi = stats.cluster_bootstrap(
        u.reduction_pct, u.index, stat=lambda v: v.mean(), n=400
    )
    assert (hi - lo) > 1.4 * (nhi - nlo)


def test_population_weighting():
    """The weighted headline, and that weighting is doing real work."""
    e = tables.population_estimate()
    assert e["weighted_mean_pct"] == 11.76
    assert e["unweighted_mean_pct"] == 19.21
    assert e["displacement_mm"] == 0.88
    # the sample over-represents heavy rain, so weighting must lower the mean
    assert e["weighted_mean_pct"] < e["unweighted_mean_pct"]
    # and cost precision
    assert e["effective_n"] < 21010


def test_weights_sum_to_the_census():
    w = io.population_weights()
    assert abs(w.sum() - len(io.census())) < 1
