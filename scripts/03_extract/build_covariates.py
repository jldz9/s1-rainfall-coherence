#!/usr/bin/env python
"""build_covariates.py -- the surface state, not just the forcing.

Rainfall is what the sky does. What the radar sees is the surface, and the
quantity that moves the backscatter is the change in soil moisture, which is not
proportional to the rain that caused it. Forty millimetres onto desert soil
moves the dielectric constant a long way; the same forty onto saturated tropical
ground runs off and moves it very little.

That one mechanism would explain three results the analysis currently reports as
unrelated: arid B losing the most coherence despite receiving the least rain, the
dose-response being concave rather than linear, and the within-burst wet/dry
enhancement being weak -- rainfall is the forcing and varies sharply in space,
while the moisture state is smooth.

So this collects, for the case and control acquisitions of every unit:

    theta      volumetric soil water, top layer, at the acquisition
    d_theta    the change between the pair's two acquisitions
    ndvi       vegetation, which intercepts rain and decorrelates on its own
    sand/clay  texture, which sets how much water a soil takes and holds

If d_theta predicts the loss better than rainfall does, the paper has a
mechanism rather than a correlation, and later sampling should target moisture
state rather than rainfall thresholds. If it does not, the dielectric story is
wrong and something else -- vegetation motion, erosion, surface reworking -- is
doing the work.

Usage:
    python build_covariates.py --root <wx root> --shard N --of M
"""
from __future__ import annotations

import argparse
import datetime as dt
import logging
import time
from pathlib import Path

import ee
import pandas as pd

log = logging.getLogger("covar")
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
HALF_LON, HALF_LAT = 0.40, 0.12


def sample(lon, lat, t_ref, t_sec):
    reg = ee.Geometry.Rectangle([lon - HALF_LON, lat - HALF_LAT,
                                 lon + HALF_LON, lat + HALF_LAT])
    era = ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR").select(
        ["volumetric_soil_water_layer_1", "volumetric_soil_water_layer_2"])

    def theta(t):
        d = ee.Date(t.strftime("%Y-%m-%d"))
        return era.filterDate(d, d.advance(1, "day")).first()

    th_r, th_s = theta(t_ref), theta(t_sec)
    # NDVI nearest the second acquisition; 16-day composite, so a window is used
    ndvi = (ee.ImageCollection("MODIS/061/MOD13A2").select("NDVI")
            .filterDate((t_sec - dt.timedelta(days=24)).isoformat(),
                        (t_sec + dt.timedelta(days=24)).isoformat())
            .mean().multiply(0.0001).rename("ndvi"))
    # texture, top 5 cm, as a static property of the site
    soil = (ee.Image("projects/soilgrids-isric/sand_mean").select("sand_0-5cm_mean")
            .rename("sand")
            .addBands(ee.Image("projects/soilgrids-isric/clay_mean")
                      .select("clay_0-5cm_mean").rename("clay")))
    img = (th_r.rename(["th1_ref", "th2_ref"])
           .addBands(th_s.rename(["th1_sec", "th2_sec"]))
           .addBands(ndvi).addBands(soil))
    return img.reduceRegion(ee.Reducer.mean(), reg, 1000,
                            maxPixels=int(1e8), bestEffort=True).getInfo()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--of", type=int, default=1)
    ap.add_argument("--project", default="ifgqual")
    a = ap.parse_args()
    ee.Initialize(project=a.project)
    R = a.root

    man = pd.read_parquet(R / "plan" / "manifest_verified.parquet")
    man["unit_id"] = [f"{b}_{pd.Timestamp(s).strftime('%Y%m%d')}"
                      for b, s in zip(man.burst_id_jpl, man.case_sec)]
    man = man.drop_duplicates("unit_id").set_index("unit_id")
    cas = pd.read_parquet(R / "plan" / "cascade_staged.parquet").drop_duplicates("unit_id")
    prof = pd.read_parquet(R / "static" / "burst_profile.parquet").set_index("burst_id_jpl")
    cas = cas[cas.unit_id.isin(man.index)].iloc[a.shard::a.of].reset_index(drop=True)
    log.info("units in shard: %s", f"{len(cas):,}")

    rows, t0 = [], time.time()
    for i, u in enumerate(cas.itertuples(), 1):
        m = man.loc[u.unit_id]
        lon = float(prof.at[u.burst_id_jpl, "lon_c"])
        lat = float(prof.at[u.burst_id_jpl, "lat_c"])
        t_cs = pd.Timestamp(m.case_sec).to_pydatetime().replace(tzinfo=None)
        t_cr = pd.Timestamp(m.case_ref).to_pydatetime().replace(tzinfo=None)
        d3s = str(u.new_ctrl_sec); d3r = str(u.new_ctrl_ref)
        t_ks = t_cs.replace(year=int(d3s[:4]), month=int(d3s[4:6]), day=int(d3s[6:]))
        t_kr = t_cs.replace(year=int(d3r[:4]), month=int(d3r[4:6]), day=int(d3r[6:]))
        rec = {"unit_id": u.unit_id, "event": u.event}
        for role, (tr, ts) in (("case", (t_cr, t_cs)), ("ctrl", (t_kr, t_ks))):
            try:
                r = sample(lon, lat, tr, ts)
            except Exception as e:                                # noqa: BLE001
                log.warning("  %s %s: %s", u.unit_id, role, str(e)[:80])
                continue
            for k, v in r.items():
                rec[f"{role}_{k}"] = v
        rows.append(rec)
        if i % 25 == 0:
            el = time.time() - t0
            log.info("  %s/%s  %.1f s/unit  eta %.0f min", f"{i:,}", f"{len(cas):,}",
                     el / i, (len(cas) - i) * el / i / 60)
    out = R / "plan" / f"covar_shard{a.shard:02d}.parquet"
    pd.DataFrame(rows).to_parquet(out, index=False)
    log.info("wrote %s rows -> %s", f"{len(rows):,}", out)


if __name__ == "__main__":
    main()
