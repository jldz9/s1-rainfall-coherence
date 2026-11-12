#!/usr/bin/env python
"""build_burst_grid.py -- the static Sentinel-1 burst grid in lon/lat.

The OPERA burst DB (opera-burst-bbox-only.sqlite3, v0.10.0) holds 1,127,661
land bursts, but each bbox is stored in its OWN UTM zone, which is useless for a
global spatial join.  This reprojects every bbox to EPSG:4326 once and writes a
flat parquet that later steps can load into memory (~1.1 M rows, ~60 MB).

Bursts, not scenes, are the unit the campaign marks as "affected": a burst is
~20 x 25 km, close to the scale of a convective rain cell or a flood reach,
whereas a 250 x 250 km IW scene footprint would count an event hundreds of km
away as a hit.

Track number and subswath are parsed out of burst_id_jpl (t001_000001_iw1) so
the scene -> burst join can be restricted to the scene's own path before any
geometry is touched -- that turns a 1.1 M-row intersect into a ~4,000-row one.

Caveat: the DB is NOT land-filtered.  Its metadata carries land_buffer_deg = 0.3
but land_optimized = 0, so the optimisation was never applied and all 1,127,661
bursts are present, ocean included.  Measured against the Koppen raster, 735,622
of them (66%) sit on open water -- which matches Earth's ocean fraction, so this
is the full global grid, not a land subset.  Filter on a valid Koppen class (or
land_fraction from profile_bursts.py) before treating a burst count as land.

Usage:
    python build_burst_grid.py --out /smith-scratch/PROJECTS/InSAR/extreme_wx/timeline
"""
from __future__ import annotations

import argparse
import logging
import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd
from pyproj import Transformer

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("burstgrid")

DB_CANDIDATES = [
    Path.home() / "data" / "s1-burst-db" / "opera-burst-bbox-only.sqlite3",
    Path.home() / "dev" / "InSARHub" / "p56" / "s1-burst-db" / "opera-burst-bbox-only.sqlite3",
]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--db", type=Path, default=None)
    a = ap.parse_args()

    db = a.db or next((p for p in DB_CANDIDATES if p.exists()), None)
    if db is None:
        raise SystemExit(f"burst DB not found; looked in {DB_CANDIDATES}")
    log.info("burst DB %s", db)

    con = sqlite3.connect(db)
    df = pd.read_sql("SELECT burst_id_jpl, epsg, xmin, ymin, xmax, ymax FROM burst_id_map", con)
    con.close()
    log.info("read %s bursts", f"{len(df):,}")

    # t001_000001_iw1 -> track 1, subswath 1
    df["track"] = df.burst_id_jpl.str.slice(1, 4).astype("int16")
    df["iw"] = df.burst_id_jpl.str.slice(-1).astype("int8")

    lon_min = np.empty(len(df)); lon_max = np.empty(len(df))
    lat_min = np.empty(len(df)); lat_max = np.empty(len(df))
    for epsg, idx in df.groupby("epsg").groups.items():
        i = df.index.get_indexer(idx)
        tr = Transformer.from_crs(int(epsg), 4326, always_xy=True)
        sub = df.iloc[i]
        # all four corners: a UTM bbox is not a lon/lat bbox, so take the hull
        xs = np.concatenate([sub.xmin, sub.xmax, sub.xmin, sub.xmax])
        ys = np.concatenate([sub.ymin, sub.ymin, sub.ymax, sub.ymax])
        lo, la = tr.transform(xs, ys)
        lo = lo.reshape(4, -1); la = la.reshape(4, -1)
        lon_min[i], lon_max[i] = lo.min(0), lo.max(0)
        lat_min[i], lat_max[i] = la.min(0), la.max(0)

    df["lon_min"], df["lon_max"] = lon_min, lon_max
    df["lat_min"], df["lat_max"] = lat_min, lat_max
    df["lon_c"] = (df.lon_min + df.lon_max) / 2
    df["lat_c"] = (df.lat_min + df.lat_max) / 2

    # A burst straddling the antimeridian gets a bbox spanning ~360 deg; flag it
    # rather than silently letting it intersect every event on Earth.
    df["dateline"] = (df.lon_max - df.lon_min) > 180
    log.info("dateline-straddling bursts: %d", int(df.dateline.sum()))

    a.out.mkdir(parents=True, exist_ok=True)
    fp = a.out / "burst_grid.parquet"
    df.drop(columns=["epsg", "xmin", "ymin", "xmax", "ymax"]).to_parquet(fp, index=False)
    log.info("wrote %s  (%.1f MB)", fp, fp.stat().st_size / 1e6)
    log.info("tracks %d..%d | mean burst %.1f x %.1f km",
             df.track.min(), df.track.max(),
             float(((df.lon_max - df.lon_min) * 111 * np.cos(np.radians(df.lat_c)))[~df.dateline].mean()),
             float(((df.lat_max - df.lat_min) * 111).mean()))


if __name__ == "__main__":
    main()
