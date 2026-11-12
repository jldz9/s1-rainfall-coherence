#!/usr/bin/env python
"""Coherence and accumulation, both restricted to the ground that was rained on.

The masked run that preceded this one took the coherence over the wetted cells
but left the accumulation as the burst mean, which measures the response over
one set of ground and the dose over another.  Below 3 mm the two differ by a
factor of 1.65; above it they do not differ at all, the storm covering the
burst.  This retrieves the IMERG field itself rather than a threshold of it, so
that the dose can be averaged over exactly the pixels whose coherence entered
the response.

Everything is taken over one pixel set per unit -- wetted, land, and finite and
positive in both interferograms -- and the all-pixel counterparts are carried
alongside, so masked and unmasked remain a within-unit comparison.

Usage:
    python extract_rainfull.py --root <wx root> --shard N --of M [--skip-done]
"""
from __future__ import annotations

import argparse, io, logging, time
from pathlib import Path
import numpy as np, pandas as pd, requests, rasterio

from extract_rainmasked import cor_path                          # noqa: E402

log = logging.getLogger("rainfull")
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
WET_MM = 1.0            # a cell counts as wetted at or above this 24 h total.
# 1 mm is the ETCCDI wet-day threshold (Zhang et al. 2011): CDD counts days
# below it, CWD counts days at or above it, PRCPTOT sums the latter.  Using the
# standard keeps the split comparable with the climate-extremes literature.
# The same value, independently, is the criterion deciding whether an
# acquisition is dry enough to enter a control pair (build_campaign.py).


def rain_acc_for(group: Path, ref_tif: Path, t_event, ee_project: str):
    """The IMERG 24 h accumulation on the coherence grid, in millimetres.

    Cached as float16: the field is a step function over 0.1 deg blocks, so the
    fine grid holds very few distinct values and half precision costs nothing
    a rainfall retrieval could resolve.
    """
    cache = group / f"rainacc_{pd.Timestamp(t_event).strftime('%Y%m%d')}.npy"
    if cache.exists():
        try:
            return np.load(cache).astype("float32")
        except Exception:                                        # noqa: BLE001
            pass
    import ee
    if not getattr(rain_acc_for, "_init", False):
        ee.Initialize(project=ee_project)
        rain_acc_for._init = True
    with rasterio.open(ref_tif) as s:
        H, W = s.shape
        T = s.transform
        crs = s.crs.to_string()
    d1 = ee.Date(pd.Timestamp(t_event).strftime("%Y-%m-%dT%H:%M:%S"))
    acc = (ee.ImageCollection("NASA/GPM_L3/IMERG_V07").select("precipitation")
           .filterDate(d1.advance(-24, "hour"), d1).sum().multiply(0.5))
    # nearest neighbour, exactly as the boolean mask was taken: an IMERG cell is
    # a block, and interpolating it would invent a gradient across the wet edge
    img = acc.rename("acc").reproject(
        crs=crs, crsTransform=[T.a, T.b, T.c, T.d, T.e, T.f])
    for attempt in range(4):
        try:
            url = img.getDownloadURL({
                "crs": crs, "crs_transform": [T.a, T.b, T.c, T.d, T.e, T.f],
                "dimensions": f"{W}x{H}", "format": "NPY"})
            r = requests.get(url, timeout=600)
            r.raise_for_status()
            a = np.load(io.BytesIO(r.content), allow_pickle=False)
            m = a[a.dtype.names[0]].astype("float32")
            if m.shape != (H, W):
                log.warning("%s: acc shape %s != %s", group.name, m.shape, (H, W))
                return None
            np.save(cache, m.astype("float16"))
            return m
        except Exception as e:                                   # noqa: BLE001
            log.warning("  %s rainacc attempt %d: %s", group.name, attempt, str(e)[:110])
            time.sleep(5 * (attempt + 1))
    return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--of", type=int, default=1)
    ap.add_argument("--tag", default="")
    ap.add_argument("--skip-done", action="store_true")
    ap.add_argument("--wet-mm", type=float, default=WET_MM,
                    help="24 h accumulation at or above which a cell is wetted")
    ap.add_argument("--project", default="ifgqual")
    a = ap.parse_args()
    R = a.root

    d = pd.read_parquet(R / "plan" / "units_with_coverage.parquet")
    if a.skip_done:
        done = set()
        for f in sorted((R / "plan").glob(f"rainfull{a.tag}_shard*.parquet")):
            try:
                done |= set(pd.read_parquet(f, columns=["unit_id"]).unit_id)
            except Exception:                                    # noqa: BLE001
                pass
        before = len(d)
        d = d[~d.unit_id.isin(done)].copy()
        log.info("skip-done: %s held, %s remain of %s",
                 f"{len(done):,}", f"{len(d):,}", f"{before:,}")

    bycase = {}
    for g in (R / "proc" / "burst").glob("*"):
        if (g / ".done").exists():
            bycase.setdefault(g.name.rsplit("_", 1)[0].lower(), []).append(g)
    d = d.iloc[a.shard::a.of].reset_index(drop=True)
    log.info("units in shard: %s", f"{len(d):,}")

    out = R / "plan" / f"rainfull{a.tag}_shard{a.shard:03d}.parquet"
    rows, t0 = [], time.time()
    for i, r in enumerate(d.itertuples(), 1):
        hit = None
        for g in bycase.get(r.burst, []):
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
        acc = rain_acc_for(g, c_case, r.t_event, a.project)
        if acc is None:
            continue
        try:
            with rasterio.open(c_case) as s:
                gc = s.read(1).astype("float32")
            with rasterio.open(c_ctrl) as s:
                gk = s.read(1).astype("float32")
        except Exception as e:                                   # noqa: BLE001
            log.warning("  %s: %s", r.unit_id, str(e)[:90])
            continue
        if gc.shape != acc.shape or gk.shape != acc.shape:
            continue
        lcf = g / "worldcover.npy"
        land = np.load(lcf) != 80 if lcf.exists() else np.ones(acc.shape, bool)
        base = land & np.isfinite(gc) & np.isfinite(gk) & (gc > 0) & (gk > 0)
        wet = base & (acc >= a.wet_mm)
        if base.sum() < 200:
            continue
        row = {"unit_id": r.unit_id, "burst": r.burst, "wet_mm": a.wet_mm,
               "rain_bursted": float(r.rain_mm),
               "g_case_all": float(gc[base].mean()), "g_ctrl_all": float(gk[base].mean()),
               "rain_all": float(acc[base].mean()), "n_px_all": int(base.sum())}
        if wet.sum() >= 200:
            row |= {"g_case_wet": float(gc[wet].mean()),
                    "g_ctrl_wet": float(gk[wet].mean()),
                    "rain_wet": float(acc[wet].mean()),
                    "rain_wet_max": float(acc[wet].max()),
                    "n_px_wet": int(wet.sum())}
        else:
            row |= {"g_case_wet": np.nan, "g_ctrl_wet": np.nan, "rain_wet": np.nan,
                    "rain_wet_max": np.nan, "n_px_wet": int(wet.sum())}
        rows.append(row)
        if i % 50 == 0:
            el = time.time() - t0
            log.info("  %s / %s   kept %s   eta %.0f min", f"{i:,}", f"{len(d):,}",
                     f"{len(rows):,}", (len(d) - i) * el / i / 60)
        if rows and i % 100 == 0:
            tmp = out.with_suffix(".tmp")
            pd.DataFrame(rows).to_parquet(tmp, index=False)
            tmp.replace(out)

    pd.DataFrame(rows).to_parquet(out, index=False)
    log.info("wrote %s (%s rows)", out, f"{len(rows):,}")


if __name__ == "__main__":
    main()
