#!/usr/bin/env python
"""profile_bursts.py -- land cover and Koppen class for every Sentinel-1 burst.

Turns the affected table into the (event x land cover x Koppen) cell counts that
decide whether this study is a paper or a pilot.

Two stratifiers, two very different roles
-----------------------------------------
LAND COVER carries the response.  It varies WITHIN a burst -- a single burst in
the parent project spanned grassland, tree cover, shrubland and cropland -- so an
event effect can be estimated as a within-burst contrast, immune to the site
effects that dominate everything else (between-site spread was 0.076-0.101 for
vegetated covers, larger than the effects themselves).  Hence FRACTIONS, not a
majority label: a burst is 80 x 20 km and is never one class.

KOPPEN carries the exposure.  In the parent project's 24 sites it had 1-3 sites
per class and ICC = 1.000 -- a pure site identifier, and dropping it IMPROVED
out-of-site skill.  It only becomes estimable here because ~400k bursts spread
across the globe give many independent bursts per class.  Recorded as dominant
class plus PURITY, because the parent project's site audit found real classes as
low as 0.38 pure and that silently moved three sites into the wrong climate.

Sampling
--------
Koppen (Beck et al. 2018, 1 km) is sampled on a 9 x 5 grid spanning the burst's
TRUE extent (~80 x 20 km), not the burst_grid bbox -- that bbox is the hull of a
rotated footprint plus a 5 km margin and over-covers by ~1.5x.  Dominant class and
purity come from those 45 points.

Land cover is the mean of the 0.1 deg fraction raster over the same true extent.
Both rasters share the 0.1 deg / -180 / 85N origin, so one index formula serves
(the weather rasters run 2 rows further south, which never matters here).

Output: static/burst_profile.parquet
    burst_id_jpl, track, iw, lon_c, lat_c,
    koppen, koppen_name, koppen_group, koppen_purity,
    lc_<class> x10 fractions, land_fraction

Usage:
    python profile_bursts.py --root /smith-scratch/PROJECTS/InSAR/extreme_wx
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("profile")

KOPPEN_TIF = Path.home() / ".cache" / "ifgqual" / "climate" / "Beck_KG_V1_present_0p0083.tif"

# Beck et al. (2018) legend, codes 1-30.
KG = ["", "Af", "Am", "Aw", "BWh", "BWk", "BSh", "BSk", "Csa", "Csb", "Csc",
      "Cwa", "Cwb", "Cwc", "Cfa", "Cfb", "Cfc", "Dsa", "Dsb", "Dsc", "Dsd",
      "Dwa", "Dwb", "Dwc", "Dwd", "Dfa", "Dfb", "Dfc", "Dfd", "ET", "EF"]

GRID_N, GRID_W, GRID_RES = 85.0, -180.0, 0.1     # land-cover raster
HALF_LON_KM, HALF_LAT_KM = 40.0, 10.0            # TRUE burst half-extent
NODATA = -32768
CHUNK = 200_000


def koppen_sample(lon: np.ndarray, lat: np.ndarray, kg: np.ndarray,
                  tr) -> tuple[np.ndarray, np.ndarray]:
    """Dominant Koppen code and purity from a 9 x 5 grid over the true burst."""
    ny, nx = kg.shape
    dlat = np.linspace(-HALF_LAT_KM, HALF_LAT_KM, 5) / 111.0
    dlon_km = np.linspace(-HALF_LON_KM, HALF_LON_KM, 9)
    cosl = np.maximum(np.cos(np.radians(lat)), 0.05)

    counts = np.zeros((lon.size, 31), dtype="int16")
    for dy in dlat:
        for dx_km in dlon_km:
            la = lat + dy
            lo = lon + dx_km / (111.0 * cosl)
            r = np.clip(((90.0 - la) / 0.008333333).astype(int), 0, ny - 1)
            c = ((lo + 180.0) / 0.008333333).astype(int) % nx
            v = kg[r, c]
            counts[np.arange(lon.size), np.clip(v, 0, 30)] += 1
    counts[:, 0] = 0                       # code 0 is ocean/no-data, never dominant
    tot = counts.sum(1)
    dom = counts.argmax(1).astype("int8")
    pur = np.where(tot > 0, counts.max(1) / np.maximum(tot, 1), 0.0)
    return dom, pur.astype("float32")


def lc_sample(lon: np.ndarray, lat: np.ndarray,
              lc: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    """Mean cover fraction over the true burst extent."""
    any_b = next(iter(lc.values()))
    ny, nx = any_b.shape
    r = np.clip(((GRID_N - lat) / GRID_RES).astype(int), 0, ny - 1)
    c = ((lon - GRID_W) / GRID_RES).astype(int) % nx
    hr = max(int(HALF_LAT_KM / 111.0 / GRID_RES), 1)
    hc = np.clip((HALF_LON_KM / (111.0 * np.maximum(np.cos(np.radians(lat)), 0.05))
                  / GRID_RES).astype(int), 1, 60)
    dr = np.arange(-hr, hr + 1)
    out = {k: np.full(r.size, NODATA, dtype="int16") for k in lc}
    for w in np.unique(hc):
        sel = np.nonzero(hc == w)[0]
        dc = np.arange(-w, w + 1)
        rows = np.clip(r[sel, None, None] + dr[None, :, None], 0, ny - 1)
        cols = (c[sel, None, None] + dc[None, None, :]) % nx
        for k, arr in lc.items():
            win = arr[rows, cols].astype("float32")
            win[win <= NODATA + 1] = np.nan
            with np.errstate(invalid="ignore"):
                m = np.nanmean(win, axis=(1, 2))
            out[k][sel] = np.where(np.isnan(m), NODATA, m).astype("int16")
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, type=Path)
    a = ap.parse_args()

    g = pd.read_parquet(a.root / "timeline" / "burst_grid.parquet")
    g = g[~g.dateline].reset_index(drop=True)
    log.info("bursts: %s", f"{len(g):,}")

    with np.load(a.root / "static" / "landcover_0p1.npz") as z:
        lc = {k: z[k] for k in z.files}
    log.info("land cover: %d bands %s", len(lc), next(iter(lc.values())).shape)

    with rasterio.open(KOPPEN_TIF) as src:
        kg = src.read(1)
        tr = src.transform
    log.info("koppen: %s %s", kg.shape, kg.dtype)

    parts = []
    for i0 in range(0, len(g), CHUNK):
        sub = g.iloc[i0:i0 + CHUNK]
        lon = sub.lon_c.to_numpy(); lat = sub.lat_c.to_numpy()
        dom, pur = koppen_sample(lon, lat, kg, tr)
        d = {"burst_id_jpl": sub.burst_id_jpl.values, "track": sub.track.values,
             "iw": sub.iw.values, "lon_c": lon, "lat_c": lat,
             "koppen": dom, "koppen_purity": pur}
        d |= {f"lc_{k}": v for k, v in lc_sample(lon, lat, lc).items()}
        parts.append(pd.DataFrame(d))
        log.info("  %s / %s", f"{min(i0+CHUNK, len(g)):,}", f"{len(g):,}")

    df = pd.concat(parts, ignore_index=True)
    df["koppen_name"] = [KG[c] if 0 < c <= 30 else "" for c in df.koppen]
    df["koppen_group"] = df.koppen_name.str[:1]
    frac = df[[c for c in df.columns if c.startswith("lc_")]].replace(NODATA, np.nan)
    water = frac[["lc_water_permanent", "lc_water_seasonal"]].sum(1, min_count=1)
    df["land_fraction"] = (1000.0 - water.fillna(0)) / 1000.0

    out = a.root / "static" / "burst_profile.parquet"
    df.to_parquet(out, index=False)
    log.info("wrote %s  (%.1f MB)", out, out.stat().st_size / 1e6)
    log.info("koppen groups: %s", df.koppen_group.value_counts().to_dict())


if __name__ == "__main__":
    main()
