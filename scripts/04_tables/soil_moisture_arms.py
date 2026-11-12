#!/usr/bin/env python
"""Volumetric soil water at the four acquisitions of every unit in an arm.

build_covariates.py predates the release table: it reconstructs each unit from
the manifest and the cascade plan, and spends two Earth Engine round trips per
unit doing it.  The release table already carries the centre of the burst and
all four dates, so the reconstruction is unnecessary and the calls can be
batched, which is what makes the new arms affordable at all.

Only ERA5-Land is pulled here.  NDVI and soil texture are static or
slowly-varying and were collected once for the earlier arms; the moisture is
the quantity the analysis differences across a pair.

Usage:
    python soil_moisture_arms.py --root <wx root> --arms M12,P12,Q12,H12 \
        --shard N --of M
"""
from __future__ import annotations

import argparse, logging, time
from pathlib import Path
import pandas as pd

log = logging.getLogger("theta")
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
HALF_LON, HALF_LAT = 0.40, 0.12
CHUNK = 10


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--arms", default="M12,P12,Q12,H12")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--of", type=int, default=1)
    ap.add_argument("--project", default="ifgqual")
    ap.add_argument("--tag", default="arms",
                    help="output goes to plan/theta_<tag>_shard<N>.parquet")
    ap.add_argument("--skip-done", action="store_true",
                    help="drop units already retrieved by any earlier run")
    a = ap.parse_args()

    import ee
    ee.Initialize(project=a.project)
    era = ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR").select(
        "volumetric_soil_water_layer_1")

    u = pd.read_parquet(a.root / "release" / "rain_coherence_units.parquet")
    u = u[u.arm.isin(a.arms.split(",")) & (u.headroom >= 0.15)]
    cols = ["unit_id", "lon", "lat", "t_ref", "t_sec", "ctrl_ref", "ctrl_sec"]
    u = u[cols].dropna()
    if a.skip_done:
        done = set()
        for f in sorted((a.root / "plan").glob("theta_*_shard*.parquet")):
            done |= set(pd.read_parquet(f, columns=["unit_id"]).unit_id)
        before = len(u)
        u = u[~u.unit_id.isin(done)]
        log.info("already retrieved: %s of %s", f"{before - len(u):,}", f"{before:,}")
    u = u.iloc[a.shard::a.of].reset_index(drop=True)
    log.info("units in shard: %s", f"{len(u):,}")

    def one(r):
        """Four daily means at one burst centre, as a single band-stacked read."""
        box = ee.Geometry.Rectangle(
            [r.lon - HALF_LON, r.lat - HALF_LAT, r.lon + HALF_LON, r.lat + HALF_LAT],
            None, False)

        def day(t, name):
            d = ee.Date(pd.Timestamp(t).strftime("%Y-%m-%d"))
            return era.filterDate(d, d.advance(1, "day")).first().rename(name)

        img = ee.Image.cat(day(r.t_ref, "case_ref"), day(r.t_sec, "case_sec"),
                           day(r.ctrl_ref, "ctrl_ref"), day(r.ctrl_sec, "ctrl_sec"))
        return img.reduceRegion(ee.Reducer.mean(), box, 1000,
                                maxPixels=int(1e8), bestEffort=True)

    rows, t0 = [], time.time()
    for i0 in range(0, len(u), CHUNK):
        ch = u.iloc[i0:i0 + CHUNK]
        try:
            vals = ee.List([one(r) for r in ch.itertuples()]).getInfo()
        except Exception as e:                                    # noqa: BLE001
            log.warning("chunk %d failed: %s", i0, str(e)[:110])
            time.sleep(5)
            continue
        for r, v in zip(ch.itertuples(), vals):
            v = v or {}
            rows.append({"unit_id": r.unit_id,
                         **{k: v.get(k) for k in
                            ("case_ref", "case_sec", "ctrl_ref", "ctrl_sec")}})
        if i0 and i0 % 500 == 0:
            el = time.time() - t0
            log.info("  %s / %s   eta %.0f min", f"{i0:,}", f"{len(u):,}",
                     (len(u) - i0) * el / i0 / 60)

    out = a.root / "plan" / f"theta_{a.tag}_shard{a.shard:03d}.parquet"
    pd.DataFrame(rows).to_parquet(out, index=False)
    log.info("wrote %s (%s rows)", out, f"{len(rows):,}")


if __name__ == "__main__":
    main()
