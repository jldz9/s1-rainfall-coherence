#!/usr/bin/env python
"""build_events.py -- global daily extreme-weather rasters (stage A of step 2).

Two-stage by design
-------------------
The campaign needs to know whether extreme weather was happening AT the minute a
burst was acquired.  Asking that directly -- sample IMERG/ERA5 at each of ~68 M
burst-acquisitions -- is ~200,000 Earth Engine calls and a fast route to a quota
ban.  Almost all would return "nothing happening": extreme weather is rare in
space-time (on 2023-09-10, 10,480 of 5,227,200 cells exceeded 20 mm/hr, 0.2%).

  stage A (here)  one global daily raster per variable at native 0.1 deg,
                  downloaded once.  Candidates are then found LOCALLY with
                  numpy, at ANY threshold, as often as we like.
  stage B         only for surviving candidates, read the exact half-hourly
                  IMERG / hourly ERA5-Land image at the acquisition timestamp.

Keeping stage A dense rather than pre-thresholded is deliberate: the whole study
turns on where the "extreme" bar sits (the parent project's antecedent test gave
-0.037 at a 4 mm bar and -0.030 at a 90 mm bar), so the threshold must stay a
free parameter.  ~13 GB buys that.

Why pre-aggregated daily collections
------------------------------------
The cost driver is TEMPORAL REDUCTION, not transfer.  Measured cold-cache on the
north strip alone:

    IMERG half-hourly  .max()+.sum() over 48 images   22.7 s   10.4 MB
    ERA5-Land hourly   .min()/.max() over 24 images   10.7 s   15.6 MB
    ERA5-Land DAILY_AGGR  (pre-reduced, 1 image)       1.2 s   15.6 MB
    CHIRPS DAILY          (already daily, 1 image)      1.5 s   20.8 MB

Same pixels, ~10x cheaper.  Stage A therefore screens on daily products and
leaves sub-daily intensity to stage B, which is where it actually matters.

Beware timing this with a warm cache: an early version measured 6 s/day on a date
Earth Engine had already cached from a failed attempt.  The true cold cost was
70 s/day -- the estimate was wrong by 12x.  Always time a date never requested.

Screening variables (int16, scaled; nodata = -32768)
----------------------------------------------------
  era5_precip     total_precipitation_sum     mm x 10   (caps at 3276 mm)
  era5_snowfall   snowfall_sum                mm we x 100
  era5_snowdepth  snow_depth_max              cm   (mm would overflow:
                                                   global max is 33.3 m)
  era5_snowcov_lo snow_cover_min              % x 10
  era5_snowcov_hi snow_cover_max              % x 10
  era5_tmin       skin_temperature_min        K x 10 - 2000
  era5_tmax       skin_temperature_max        K x 10 - 2000
  chirps_precip   CHIRPS daily precipitation  mm x 10   (50S-50N only)

tmin/tmax straddling 273.15 K is the freeze-thaw flag, which is why both are
kept rather than a precomputed boolean.  snow_cover min AND max are kept because
the parent project's strongest per-cover signal was snow-cover CHANGE (bare and
sparse vegetation went -0.032 -> -0.121 monotonically with it), not snow presence.

CHIRPS spans 50S-50N and is nodata outside; ERA5-Land covers the rest.  ERA5-Land
is land-only (~30% of cells valid), which suits us -- the burst DB is land-
buffered too.

BOTH screening products UNDER-DETECT extremes, and badly.  Storm Daniel over
Derna on 2023-09-10, measured here:

    CHIRPS daily            3.3 mm     gauge-blended, and Libya has ~no gauges
    ERA5 total_precip      68.8 mm     reanalysis smooths deep convection
    IMERG daily accum     206.7 mm     (stage B; closest to the ~400 mm reported)

So stage A is a SCREEN, never a measurement: take the max over the available
products, set the bar LOW, and let stage B's IMERG at the exact half hour do the
rejecting.  A magnitude quoted from stage A would be wrong by more than an order
of magnitude in exactly the gauge-sparse regions the campaign cares about.

Sizing traps, both hit for real
-------------------------------
getDownloadURL refuses > 50 MB.  Seven int16 bands over a 721-row strip is
54.5 MB, so the globe is cut into THREE strips (~484 rows, ~24 MB) not two.
And `scale` is advisory: getDownloadURL snaps to the collection's native grid
regardless of what is asked, so the returned shape must always be verified.

Usage:
    python build_events.py --out .../events --start 2014-10-01 --end 2026-09-21
"""
from __future__ import annotations

import argparse
import datetime as dt
import io
import logging
import random
import time
from pathlib import Path

import ee
import numpy as np
import requests

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("events")

SCALE = 11132.0                                  # 0.1 deg -- advisory, see above
STRIPS = [(85, 37), (37, -11), (-11, -60)]       # ~484 rows, ~24 MB each
BACKOFF = [5, 15, 45, 120, 300]
FATAL = ("must be less than or equal", "Invalid argument", "no bands")
NODATA = -32768

ERA5_BANDS = [
    ("total_precipitation_sum",  10000, 0,    "era5_precip"),
    ("snowfall_sum",            100000, 0,    "era5_snowfall"),
    ("snow_depth_max",             100, 0,    "era5_snowdepth"),
    ("snow_cover_min",              10, 0,    "era5_snowcov_lo"),
    ("snow_cover_max",              10, 0,    "era5_snowcov_hi"),
    ("skin_temperature_min",        10, 2000, "era5_tmin"),
    ("skin_temperature_max",        10, 2000, "era5_tmax"),
]


def era5_image(day: dt.date) -> ee.Image:
    d0, d1 = day.isoformat(), (day + dt.timedelta(days=1)).isoformat()
    im = ee.Image(ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR")
                  .filterDate(d0, d1).first())
    return ee.Image.cat([
        im.select(src).multiply(mul).subtract(off).toInt16().rename(dst)
        for src, mul, off, dst in ERA5_BANDS
    ]).toInt16().unmask(NODATA)


def chirps_image(day: dt.date) -> ee.Image:
    d0, d1 = day.isoformat(), (day + dt.timedelta(days=1)).isoformat()
    c = ee.ImageCollection("UCSB-CHG/CHIRPS/DAILY").filterDate(d0, d1)
    return ee.Image(c.first()).multiply(10).toInt16().rename("chirps_precip").unmask(NODATA)


def _download(img: ee.Image, north: float, south: float) -> np.ndarray:
    url = img.getDownloadURL({
        "region": ee.Geometry.BBox(-180, south, 180, north),
        "scale": SCALE, "format": "NPY", "crs": "EPSG:4326",
    })
    r = requests.get(url, timeout=900)
    r.raise_for_status()
    return np.load(io.BytesIO(r.content), allow_pickle=False)


def fetch_day(day: dt.date, out: Path, sleep: float) -> bool:
    fp = out / f"{day:%Y}" / f"{day:%Y-%m-%d}.npz"
    if fp.exists():
        return True
    fp.parent.mkdir(parents=True, exist_ok=True)

    for attempt, wait in enumerate([0] + BACKOFF):
        if wait:
            j = wait * (0.5 + random.random())
            log.warning("  retry %d for %s in %.0fs", attempt, day, j)
            time.sleep(j)
        try:
            e5, ch, has_chirps = [], [], True
            for north, south in STRIPS:
                e5.append(_download(era5_image(day), north, south))
                time.sleep(sleep)
                if has_chirps:
                    try:
                        ch.append(_download(chirps_image(day), north, south))
                        time.sleep(sleep)
                    except Exception as ce:   # CHIRPS lags ~3 weeks near real time
                        log.warning("  %s CHIRPS unavailable (%s): %s",
                                    day, type(ce).__name__, str(ce)[:120])
                        has_chirps = False
            bands = {n: np.vstack([p[n] for p in e5]).astype("int16")
                     for n in e5[0].dtype.names}
            if has_chirps and len(ch) == len(STRIPS):
                bands["chirps_precip"] = np.vstack(
                    [p["chirps_precip"] for p in ch]).astype("int16")
            np.savez_compressed(fp, **bands)
            return True
        except Exception as e:                                        # noqa: BLE001
            msg = str(e)
            log.warning("  %s failed (%s): %s", day, type(e).__name__, msg[:150])
            if any(f in msg for f in FATAL):
                log.error("  %s: not retryable, skipping", day)
                return False
    return False


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--start", required=True)
    ap.add_argument("--end", required=True)
    ap.add_argument("--project", default="ifgqual")
    ap.add_argument("--sleep", type=float, default=0.3)
    a = ap.parse_args()

    ee.Initialize(project=a.project)
    d, last = dt.date.fromisoformat(a.start), dt.date.fromisoformat(a.end)
    total = (last - d).days + 1
    n_ok = n_bad = 0
    t0 = time.time()
    while d <= last:
        n_ok += (ok := fetch_day(d, a.out, a.sleep))
        n_bad += (not ok)
        done = n_ok + n_bad
        if done % 25 == 0:
            el = (time.time() - t0) / 60
            log.info("%s | ok %d bad %d | %.1f min | eta %.1f h",
                     d, n_ok, n_bad, el, el / done * (total - done) / 60)
        d += dt.timedelta(days=1)
    log.info("done ok=%d bad=%d in %.1f min", n_ok, n_bad, (time.time() - t0) / 60)


if __name__ == "__main__":
    main()
