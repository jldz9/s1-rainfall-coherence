#!/usr/bin/env python
"""mark_affected.py -- join the S1 acquisition timeline to the weather rasters.

Step 2 of the campaign, and the reason stages A and B are split: everything here
is LOCAL.  No API, no quota, no network.  Re-run it with a different threshold as
often as you like -- that is exactly what keeping stage A dense bought.

    scenes (timeline/scenes/*.parquet)      ~2.5 M IW SLC acquisitions
      x burst grid (timeline/burst_grid.parquet)  1,127,661 land bursts
      x daily weather (events/YYYY/*.npz)   0.1 deg global, 8 bands
      -> one row per BURST-ACQUISITION with the weather over that burst

Burst level, not scene level
---------------------------
A burst is ~80 x 20 km; an IW scene footprint is ~250 x 250 km.  Marking a whole
scene "affected" because a storm clipped one corner would put hundreds of km of
dry ground into the case group and bias every Delta-gamma toward zero.

Scene -> burst without an API call
----------------------------------
The OPERA burst grid is static, so covering bursts are found geometrically:
restrict to the scene's own track (pathNumber == track), then bracket by latitude
with searchsorted, then test longitude overlap.  That turns a 1.1 M-row intersect
into ~40 candidates per scene.

Burst acquisition time
----------------------
Interpolated along-track from the scene's startTime/stopTime by latitude, honouring
flight direction.  The whole scene spans ~33 s, so this is precise to seconds --
far finer than IMERG's 30 min or ERA5's hour, which is the point made in
build_timeline.py for pulling scene-level metadata rather than 68 M burst records.

Sampling geometry, and a caveat
-------------------------------
burst_grid bboxes are axis-aligned hulls of ROTATED footprints plus the burst
DB's 5 km margin, so they over-cover by ~1.5x (mean 99 x 64 km against a true
~80 x 20 km burst).  Using them directly would inflate the affected rate.  So the
weather is sampled over a window centred on the burst centroid and sized to the
TRUE burst extent, and both the centroid value and the window max are written --
the centroid is conservative, the max is inclusive, and the choice between them
is a downstream decision rather than one baked in here.

KNOWN BIAS -- this join OVER-ASSIGNS bursts to scenes by ~45%
-------------------------------------------------------------
Measured on S1A_IW_SLC__1SDV_20141003T165117_..._79FE: this code assigns 46
bursts (15/15/16 across IW1/2/3), but the scene spans 165117-165146 = 29 s and
the IW burst cycle is ~2.75 s, so the truth is ~10-11 per subswath, ~31 total.

The cause is geometric slack, all of it deliberate elsewhere: the burst DB adds a
5 km margin, the lon/lat bbox is the hull of a ROTATED footprint (~99 x 64 km for
a true ~80 x 20 km burst), and bursts_for_scene() pads latitude by another 0.5 deg.

This is left in place because the join is a CANDIDATE GENERATOR, not a manifest.
Burst availability has to be confirmed against ASF SLC-BURST before processing
anyway -- the parent project's find_bursts.py already does exactly that -- and
confirmation is cheap once the candidate list is down to the few thousand bursts
actually selected.  Phantom candidates cost nothing but are never processed.

Two consequences to keep in mind:
  * absolute burst-acquisition COUNTS from this step are ~45% high;
  * affected RATES are roughly unbiased, since numerator and denominator inflate
    together -- but only roughly, so do not quote them to two decimals.

To remove the bias properly, store the scene footprint POLYGON in
build_timeline.py (asf_search returns it; only the bbox is kept today) and
require the burst centroid to fall inside it.

Outputs one parquet per month:
    affected/YYYY-MM.parquet
      burst_id_jpl, track, iw, scene, acq_time, lon_c, lat_c,
      <band>_c   value at the burst centroid
      <band>_x   max over the burst window

Usage:
    python mark_affected.py --root /smith-scratch/PROJECTS/InSAR/extreme_wx
    python mark_affected.py --root ... --months 2014-10,2014-11
"""
from __future__ import annotations

import argparse
import datetime as dt
import logging
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("affected")

# The stage-A grid: -180..180, 85..-60, 0.1 deg, north-up.
GRID_N, GRID_W, GRID_RES = 85.0, -180.0, 0.1
NODATA = -32768

# True burst half-extent in degrees, used for the sampling window.  A burst is
# ~80 km across-track and ~20 km along-track; at 0.1 deg (~11 km) that is a few
# cells.  Deliberately NOT the burst_grid bbox -- see the module docstring.
HALF_LON_KM, HALF_LAT_KM = 40.0, 10.0


def load_bursts(root: Path) -> dict[int, pd.DataFrame]:
    g = pd.read_parquet(root / "timeline" / "burst_grid.parquet")
    g = g[~g.dateline].copy()          # 6,188 of 1.13 M; they would match everything
    out = {}
    for trk, sub in g.groupby("track"):
        sub = sub.sort_values("lat_c").reset_index(drop=True)
        out[int(trk)] = sub
    log.info("burst grid: %s bursts on %d tracks", f"{len(g):,}", len(out))
    return out


def day_raster(root: Path, day: dt.date) -> dict[str, np.ndarray] | None:
    fp = root / "events" / f"{day:%Y}" / f"{day:%Y-%m-%d}.npz"
    if not fp.exists():
        return None
    with np.load(fp) as z:
        return {k: z[k] for k in z.files}


def sample(bands: dict[str, np.ndarray], lon: np.ndarray, lat: np.ndarray) -> dict[str, np.ndarray]:
    """Centroid value and window max for every band, at every burst.

    Fully vectorised: bursts are grouped by their window half-width (which varies
    with latitude, because 0.1 deg of longitude is 11 km at the equator and 5 km
    at 60 N) and each group is gathered in one fancy-indexed take.  The obvious
    per-burst Python loop costs ~1.5 h over the full archive; this is seconds.

    NODATA is int16 min, so it can never win a max unless the whole window is
    nodata -- which is the answer we want there anyway.
    """
    any_band = next(iter(bands.values()))
    ny, nx = any_band.shape
    r = np.clip(((GRID_N - lat) / GRID_RES).astype(int), 0, ny - 1)
    c = np.clip(((lon - GRID_W) / GRID_RES).astype(int), 0, nx - 1)

    hr = max(int(HALF_LAT_KM / 111.0 / GRID_RES), 1)
    hc_f = HALF_LON_KM / (111.0 * np.maximum(np.cos(np.radians(lat)), 0.05)) / GRID_RES
    hc = np.clip(hc_f.astype(int), 1, 60)

    out = {f"{n}_c": a[r, c] for n, a in bands.items()}
    for n in bands:
        out[f"{n}_x"] = np.empty(r.size, dtype="int16")

    dr = np.arange(-hr, hr + 1)
    for w in np.unique(hc):
        sel = np.nonzero(hc == w)[0]
        dc = np.arange(-w, w + 1)
        rows = np.clip(r[sel, None, None] + dr[None, :, None], 0, ny - 1)
        # longitude wraps: a burst near the antimeridian must not clip to the edge
        cols = (c[sel, None, None] + dc[None, None, :]) % nx
        for n, a in bands.items():
            out[f"{n}_x"][sel] = a[rows, cols].max(axis=(1, 2))
    return out


def bursts_for_scene(bt: pd.DataFrame, row) -> np.ndarray:
    """Indices of bursts on this track whose bbox overlaps the scene bbox."""
    lo = np.searchsorted(bt.lat_c.values, row.lat_min - 0.5, side="left")
    hi = np.searchsorted(bt.lat_c.values, row.lat_max + 0.5, side="right")
    if hi <= lo:
        return np.empty(0, dtype=int)
    sl = slice(lo, hi)
    ok = ((bt.lon_max.values[sl] >= row.lon_min) & (bt.lon_min.values[sl] <= row.lon_max)
          & (bt.lat_max.values[sl] >= row.lat_min) & (bt.lat_min.values[sl] <= row.lat_max))
    return np.nonzero(ok)[0] + lo


def process_month(root: Path, month: str, tracks: dict[int, pd.DataFrame], out_dir: Path) -> int:
    fp = root / "timeline" / "scenes" / f"{month}.parquet"
    if not fp.exists():
        log.warning("%s: no timeline", month); return 0
    dst = out_dir / f"{month}.parquet"
    if dst.exists():
        log.info("%s skip (exists)", month); return 0

    sc = pd.read_parquet(fp)
    if sc.empty or "lon_min" not in sc:
        log.warning("%s: empty timeline", month); return 0
    sc["day"] = pd.to_datetime(sc.startTime).dt.date

    recs = []
    for day, chunk in sc.groupby("day"):
        bands = day_raster(root, day)
        if bands is None:
            continue
        for row in chunk.itertuples():
            bt = tracks.get(int(row.pathNumber))
            if bt is None:
                continue
            idx = bursts_for_scene(bt, row)
            if idx.size == 0:
                continue
            sub = bt.iloc[idx]
            # along-track time: interpolate by latitude between start and stop
            t0 = pd.Timestamp(row.startTime); t1 = pd.Timestamp(row.stopTime)
            span = row.lat_max - row.lat_min
            frac = ((sub.lat_c.values - row.lat_min) / span) if span > 1e-6 else np.full(len(sub), .5)
            if str(row.flightDirection).upper().startswith("DESC"):
                frac = 1.0 - frac
            t = t0 + (t1 - t0) * np.clip(frac, 0, 1)
            d = {"burst_id_jpl": sub.burst_id_jpl.values, "track": sub.track.values,
                 "iw": sub.iw.values, "scene": row.sceneName, "acq_time": t,
                 "lon_c": sub.lon_c.values, "lat_c": sub.lat_c.values}
            d |= sample(bands, sub.lon_c.values, sub.lat_c.values)
            recs.append(pd.DataFrame(d))

    if not recs:
        log.warning("%s: nothing joined (weather rasters missing?)", month); return 0
    df = pd.concat(recs, ignore_index=True)
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_parquet(dst, index=False)
    log.info("%s  %s burst-acquisitions -> %s", month, f"{len(df):,}", dst.name)
    return len(df)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--months", default=None, help="comma list, else every month present")
    a = ap.parse_args()

    tracks = load_bursts(a.root)
    months = (a.months.split(",") if a.months else
              sorted(p.stem for p in (a.root / "timeline" / "scenes").glob("*.parquet")))
    out_dir = a.root / "affected"
    total = 0
    for m in months:
        total += process_month(a.root, m, tracks, out_dir)
    log.info("total %s burst-acquisitions", f"{total:,}")


if __name__ == "__main__":
    main()
