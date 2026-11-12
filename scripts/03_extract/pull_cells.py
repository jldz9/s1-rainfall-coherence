#!/usr/bin/env python
"""Per-cell soil moisture and rainfall over each burst, rather than one mean.

A burst spans about five cells of the 0.1 degree grid that both ERA5-Land and
IMERG are served on. Reducing those to a single mean, as the first pull did,
discards the only information that would let moisture and rainfall act as
within-burst contrasts, which is how land cover carries its effect: the
coherence raster is per pixel, so a covariate that is one number per burst can
only ever explain differences between bursts, never inside one.

The grids are kept at their native resolution rather than resampled onto the
coherence grid. Interpolating a 10 km field onto a 10 m one invents structure
that is not in the measurement, and costs four orders of magnitude in storage.
"""
import argparse, logging, time
from pathlib import Path
import numpy as np, pandas as pd

HALF = 0.10          # degrees; roughly a burst
log = logging.getLogger("cells")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--of", type=int, default=1)
    ap.add_argument("--project", default="ifgqual")
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    import ee
    ee.Initialize(project=a.project)
    era = ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR").select(
        "volumetric_soil_water_layer_1")
    imerg = ee.ImageCollection("NASA/GPM_L3/IMERG_V07").select("precipitation")

    u = pd.read_parquet(a.root / "release" / "rain_coherence_units.parquet")
    u = u.iloc[a.shard::a.of].reset_index(drop=True)
    log.info("units in shard: %s", f"{len(u):,}")

    def box(lon, lat):
        return ee.Geometry.Rectangle([lon - HALF, lat - HALF, lon + HALF, lat + HALF],
                                     None, False)

    def theta_grid(lon, lat, ts):
        d = ee.Date(pd.Timestamp(ts).strftime("%Y-%m-%d"))
        im = era.filterDate(d, d.advance(1, "day")).first()
        return im.sampleRectangle(box(lon, lat), defaultValue=-9999)

    # a reduced collection carries no projection of its own, so sampleRectangle
    # would fall back to the default 1 degree grid and return a single cell.
    # The projection of a member image is reimposed to keep IMERG's 0.1 degree.
    IMERG_PROJ = ee.ImageCollection("NASA/GPM_L3/IMERG_V07").first() \
                   .select("precipitation").projection()

    def rain_grid(lon, lat, ts):
        end = ee.Date(pd.Timestamp(ts).strftime("%Y-%m-%dT%H:%M:%S"))
        acc = imerg.filterDate(end.advance(-24, "hour"), end).sum().multiply(0.5)
        return (acc.rename("precipitation").reproject(IMERG_PROJ)
                   .sampleRectangle(box(lon, lat), defaultValue=-9999))

    rows, t0, CH = [], time.time(), 12
    for i0 in range(0, len(u), CH):
        sub = u.iloc[i0:i0 + CH]
        feats = []
        for r in sub.itertuples():
            wet = r.t_sec if r.rain_at_sec else r.t_ref
            dry = r.t_ref if r.rain_at_sec else r.t_sec
            feats.append(theta_grid(r.lon, r.lat, r.ctrl_ref).set("uid", r.unit_id).set("k", "th_base"))
            feats.append(theta_grid(r.lon, r.lat, wet).set("uid", r.unit_id).set("k", "th_wet"))
            feats.append(theta_grid(r.lon, r.lat, dry).set("uid", r.unit_id).set("k", "th_dry"))
            feats.append(rain_grid(r.lon, r.lat, wet).set("uid", r.unit_id).set("k", "rain"))
        try:
            got = ee.FeatureCollection(feats).getInfo()["features"]
        except Exception as e:                                  # noqa: BLE001
            log.warning("chunk %d failed: %s", i0, str(e)[:110]); continue
        acc = {}
        for f in got:
            p = f["properties"]; key = p.get("k"); uid = p.get("uid")
            band = "precipitation" if key == "rain" else "volumetric_soil_water_layer_1"
            grid = p.get(band)
            if grid is None:
                continue
            arr = np.array(grid, dtype="float32")
            arr[arr <= -9998] = np.nan
            e = acc.setdefault(uid, {"unit_id": uid})
            e[f"{key}_mean"] = float(np.nanmean(arr)) if np.isfinite(arr).any() else np.nan
            e[f"{key}_sd"] = float(np.nanstd(arr)) if np.isfinite(arr).any() else np.nan
            e[f"{key}_min"] = float(np.nanmin(arr)) if np.isfinite(arr).any() else np.nan
            e[f"{key}_max"] = float(np.nanmax(arr)) if np.isfinite(arr).any() else np.nan
            e[f"{key}_n"] = int(np.isfinite(arr).sum())
            e[f"{key}_cells"] = arr.ravel().tolist()
            e[f"{key}_shape"] = list(arr.shape)
        rows += list(acc.values())
        if (i0 // CH) % 10 == 0:
            done = min(i0 + CH, len(u))
            log.info("  %d/%d  %.2f s/unit", done, len(u), (time.time() - t0) / max(done, 1))

    d = pd.DataFrame(rows)
    dst = a.root / "plan" / f"cells_shard{a.shard:03d}.parquet"
    d.to_parquet(dst, index=False)
    log.info("wrote %s rows -> %s", f"{len(d):,}", dst)


if __name__ == "__main__":
    main()
