#!/usr/bin/env python
"""Coherence loss restricted to the ground that actually received the rain.

The dose is an IMERG accumulation averaged over the burst, while the coherence
is averaged over every land pixel in it.  Where the storm covered the burst the
two describe the same ground, and above about 7.5 mm they do: coverage above
1 mm is 1.00 and the footprint mean equals the wet-cell mean to two decimals.
Below 5 mm they do not.  In the 1-2 mm band the rain reaches less than half the
burst, so the reported loss is an average over ground much of which stayed dry,
and the measurement understates what the rain did where it fell.

This re-reads the stored coherence rasters and takes the same paired difference
over the wet cells only.  The rain field is fetched onto the coherence grid
exactly as the land cover is -- same CRS, same transform, same dimensions --
and cached beside it, so a cell is wet or dry at 40 m resolution rather than
being resampled.

Only units below the coverage threshold need this; above it the answer is
already correct, and the script leaves them alone.

Usage:
    python extract_rainmasked.py --root <wx root> --shard N --of M [--max-mm 5]
"""
from __future__ import annotations

import argparse, io, logging, time
from pathlib import Path
import numpy as np, pandas as pd, requests, rasterio

log = logging.getLogger("rainmask")
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
WET_MM = 1.0            # a cell counts as wetted at or above this 24 h total.
# 1 mm is the ETCCDI wet-day threshold (Zhang et al. 2011): CDD counts days
# below it, CWD counts days at or above it, PRCPTOT sums the latter.  Using the
# standard keeps the split comparable with the climate-extremes literature.
# The same value, independently, is the criterion deciding whether an
# acquisition is dry enough to enter a control pair (build_campaign.py).


def rain_mask_for(group: Path, ref_tif: Path, t_event, ee_project: str):
    """The IMERG 24 h field on the coherence grid, as a boolean wet mask."""
    cache = group / f"rainmask_{pd.Timestamp(t_event).strftime('%Y%m%d')}.npy"
    if cache.exists():
        return np.load(cache)
    import ee
    if not getattr(rain_mask_for, "_init", False):
        ee.Initialize(project=ee_project)
        rain_mask_for._init = True
    with rasterio.open(ref_tif) as s:
        H, W = s.shape
        T = s.transform
        crs = s.crs.to_string()
    d1 = ee.Date(pd.Timestamp(t_event).strftime("%Y-%m-%dT%H:%M:%S"))
    acc = (ee.ImageCollection("NASA/GPM_L3/IMERG_V07").select("precipitation")
           .filterDate(d1.advance(-24, "hour"), d1).sum().multiply(0.5))
    # nearest-neighbour onto the fine grid: an IMERG cell is a block, not a
    # surface to be interpolated, and interpolating would blur the wet edge
    img = acc.gt(WET_MM).rename("wet").reproject(crs=crs, crsTransform=[T.a, T.b, T.c, T.d, T.e, T.f])
    for attempt in range(4):
        try:
            url = img.getDownloadURL({
                "crs": crs, "crs_transform": [T.a, T.b, T.c, T.d, T.e, T.f],
                "dimensions": f"{W}x{H}", "format": "NPY"})
            r = requests.get(url, timeout=600)
            r.raise_for_status()
            a = np.load(io.BytesIO(r.content), allow_pickle=False)
            m = a[a.dtype.names[0]].astype(bool)
            if m.shape != (H, W):
                log.warning("%s: mask shape %s != %s", group.name, m.shape, (H, W))
                return None
            np.save(cache, m)
            return m
        except Exception as e:                                   # noqa: BLE001
            log.warning("  %s rainmask attempt %d: %s", group.name, attempt, str(e)[:110])
            time.sleep(5 * (attempt + 1))
    return None


def cor_path(root: Path, group: str, ref: str, sec: str) -> Path | None:
    p = root / "proc" / "burst" / group / "pairs" / f"{ref}_{sec}" / "interferograms"
    f = p / f"{ref}_{sec}.int.cor.tif"
    return f if f.exists() else None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--of", type=int, default=1)
    ap.add_argument("--max-mm", type=float, default=5.0,
                    help="only units at or below this accumulation need masking; "
                         "pass a large value to mask the whole range")
    ap.add_argument("--tag", default="",
                    help="suffix for the output shard, so a full-range run does "
                         "not overwrite the light-rain one")
    ap.add_argument("--skip-done", action="store_true",
                    help="skip units any earlier rainmasked shard already holds, "
                         "so a rerun spends no Earth Engine quota on them")
    ap.add_argument("--project", default="ifgqual")
    a = ap.parse_args()
    R = a.root

    d = pd.read_parquet(R / "plan" / "units_with_coverage.parquet")
    d = d[d.rain_mm <= a.max_mm].copy()
    if a.skip_done:
        done = set()
        for f in sorted((R / "plan").glob("rainmasked*_shard*.parquet")) + \
                 sorted((R / "plan").glob("rainmasked_shard*.parquet")):
            try:
                done |= set(pd.read_parquet(f, columns=["unit_id"]).unit_id)
            except Exception:                                    # noqa: BLE001
                pass
        before = len(d)
        d = d[~d.unit_id.isin(done)].copy()
        log.info("skip-done: %s already masked, %s remain of %s",
                 f"{len(done):,}", f"{len(d):,}", f"{before:,}")
    # the group directory is not carried on the unit; recover it from the burst
    bycase = {}
    for g in (R / "proc" / "burst").glob("*"):
        if (g / ".done").exists():
            bycase.setdefault(g.name.rsplit("_", 1)[0].lower(), []).append(g)
    d = d.iloc[a.shard::a.of].reset_index(drop=True)
    log.info("units in shard: %s (at or below %s mm)", f"{len(d):,}", a.max_mm)

    out = R / "plan" / f"rainmasked{a.tag}_shard{a.shard:03d}.parquet"
    rows, t0 = [], time.time()
    for i, r in enumerate(d.itertuples(), 1):
        groups = bycase.get(r.burst, [])
        hit = None
        for g in groups:
            c_case = cor_path(R, g.name,
                              pd.Timestamp(r.t_ref).strftime("%Y%m%d"),
                              pd.Timestamp(r.t_sec).strftime("%Y%m%d"))
            c_ctrl = cor_path(R, g.name,
                              pd.Timestamp(r.ctrl_ref).strftime("%Y%m%d"),
                              pd.Timestamp(r.ctrl_sec).strftime("%Y%m%d"))
            if c_case and c_ctrl:
                hit = (g, c_case, c_ctrl)
                break
        if hit is None:
            continue
        g, c_case, c_ctrl = hit
        m = rain_mask_for(g, c_case, r.t_event, a.project)
        if m is None or not m.any():
            continue
        try:
            with rasterio.open(c_case) as s:
                gc = s.read(1).astype("float32")
            with rasterio.open(c_ctrl) as s:
                gk = s.read(1).astype("float32")
        except Exception as e:                                   # noqa: BLE001
            log.warning("  %s: %s", r.unit_id, str(e)[:90])
            continue
        if gc.shape != m.shape or gk.shape != m.shape:
            continue
        lc = g / "worldcover.npy"
        land = np.load(lc) != 80 if lc.exists() else np.ones_like(m, bool)
        ok = m & land & np.isfinite(gc) & np.isfinite(gk) & (gc > 0) & (gk > 0)
        allpx = land & np.isfinite(gc) & np.isfinite(gk) & (gc > 0) & (gk > 0)
        if ok.sum() < 200:
            continue
        rows.append({"unit_id": r.unit_id,
                     "g_case_wet": float(gc[ok].mean()), "g_ctrl_wet": float(gk[ok].mean()),
                     "g_case_all": float(gc[allpx].mean()), "g_ctrl_all": float(gk[allpx].mean()),
                     "n_px_wet": int(ok.sum()), "n_px_all": int(allpx.sum())})
        if i % 50 == 0:
            el = time.time() - t0
            log.info("  %s / %s   kept %s   eta %.0f min", f"{i:,}", f"{len(d):,}",
                     f"{len(rows):,}", (len(d) - i) * el / i / 60)
        if rows and i % 100 == 0:
            # checkpoint: a walltime kill must not discard hours of retrieval
            tmp = out.with_suffix(".tmp")
            pd.DataFrame(rows).to_parquet(tmp, index=False)
            tmp.replace(out)

    pd.DataFrame(rows).to_parquet(out, index=False)
    log.info("wrote %s (%s rows)", out, f"{len(rows):,}")


if __name__ == "__main__":
    main()
