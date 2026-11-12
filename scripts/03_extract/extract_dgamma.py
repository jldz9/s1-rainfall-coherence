#!/usr/bin/env python
"""extract_dgamma.py -- per-pixel, per-land-cover coherence loss for each unit.

The payoff step: turns 5,792 coherence rasters into the
(event x land cover x Koppen) matrix the campaign was built to fill.

Why PER PIXEL and PAIRED
------------------------
Case and control are the same burst, and COMPASS geocodes a burst to the same
fixed UTM grid whatever the date -- verified identical (shape, transform, CRS)
on 120 of 120 units checked.  So the difference can be taken pixel by pixel
rather than between two burst means.

That matters twice over.  It removes the static scene from the estimate, so what
is left is the event rather than the terrain; and it is the only way land cover
can act as a WITHIN-BURST contrast, which is the whole reason cover carries the
response while Koppen carries the exposure.  A burst-level label cannot do this:
a burst is 80 x 20 km and never one cover class.

Land cover comes from ESA WorldCover 10 m, resampled by Earth Engine onto the
coherence grid itself (crs_transform + dimensions, never `scale`, which is
advisory and snaps).  ~1.1 s per burst, cached as worldcover.npy.

WATER IS MASKED. The probe burst was 1,049,290 water pixels against 15,921 of
grass; unmasked, a coastal burst reports the coherence of the sea.

Output: one tidy row per (unit x cover):
    unit_id, event, cover, koppen, magnitude, btemp_days, n_px,
    g_case, g_ctrl, dgamma (mean per-pixel case-control),
    dgamma_med (median of the same, robust to fully-decorrelated pixels)

Usage:
    python extract_dgamma.py --root ... [--groups-from FILE] [--shard I --of N]
"""
from __future__ import annotations

import argparse
import io
import logging
import time
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
import requests

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("dgamma")

WC = {10: "tree", 20: "shrub", 30: "grass", 40: "crops", 50: "built",
      60: "bare", 70: "snowice", 80: "water", 90: "wetland", 95: "mangrove",
      100: "moss"}
DROP = {80}                 # water: not a land surface, and it can dominate a burst
MIN_PX = 200                # below this a per-cover mean is noise


def worldcover_for(group: Path, ref_tif: Path, ee_project: str) -> np.ndarray | None:
    cache = group / "worldcover.npy"
    if cache.exists():
        return np.load(cache)
    import ee
    if not getattr(worldcover_for, "_init", False):
        ee.Initialize(project=ee_project)
        worldcover_for._init = True
    with rasterio.open(ref_tif) as s:
        H, W = s.shape
        T = s.transform
        crs = s.crs.to_string()
    img = ee.Image(ee.ImageCollection("ESA/WorldCover/v200").first()).select("Map")
    for attempt in range(4):
        try:
            url = img.getDownloadURL({
                "crs": crs, "crs_transform": [T.a, T.b, T.c, T.d, T.e, T.f],
                "dimensions": f"{W}x{H}", "format": "NPY"})
            r = requests.get(url, timeout=600)
            r.raise_for_status()
            a = np.load(io.BytesIO(r.content), allow_pickle=False)
            lc = a[a.dtype.names[0]].astype("uint8")
            if lc.shape != (H, W):
                log.warning("%s: worldcover shape %s != %s", group.name, lc.shape, (H, W))
                return None
            np.save(cache, lc)
            return lc
        except Exception as e:                                   # noqa: BLE001
            log.warning("  %s worldcover attempt %d: %s", group.name, attempt, str(e)[:110])
            time.sleep(5 * (attempt + 1))
    return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--out", default=None, type=Path)
    ap.add_argument("--project", default="ifgqual")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--of", type=int, default=1)
    # proc/burst holds every campaign ever run, so a run that means one of them
    # must say which; the default keeps the old behaviour of taking them all
    ap.add_argument("--pattern", default="*",
                    help="glob over group directory names, e.g. '*_r12*'")
    a = ap.parse_args()

    groups = sorted(p for p in (a.root / "proc" / "burst").glob(a.pattern)
                    if p.is_dir() and (p / ".done").exists()
                    and not (p / ".parked").exists())
    groups = groups[a.shard::a.of]
    log.info("groups in this shard: %s", f"{len(groups):,}")

    rows = []
    for i, G in enumerate(groups, 1):
        try:
            pl = pd.read_csv(G / "pair_list.csv")
        except Exception:                                        # noqa: BLE001
            continue
        lc = None
        for uid, sub in pl.groupby("unit_id"):
            roles = set(sub.role)
            if "case" not in roles:
                continue
            # control  = the original match, which IMERG shows was often NOT dry
            #            (36% were wetter than their own case)
            # control2 = the v1 re-selection: driest by IMERG, calendar gap
            #            uncapped, which reaches 11 years back
            # control3 = the v4 re-selection: same search, gap capped at 4 years
            # All available arms are emitted, tagged by `contrast`, so each
            # selection rule reads as a before/after on the SAME units rather
            # than being asserted.
            contrasts = [("v1", "control")]
            if "control2" in roles:
                contrasts.append(("v2", "control2"))
            if "control3" in roles:
                contrasts.append(("v4", "control3"))
            # control_sy = the same-year arm: the same burst a few weeks away,
            # which the seasonal series shows reproduces a weather-free pair
            # about 3.5x more faithfully than any control a year or more away
            # (11.2% of baseline against 25.1%).
            if "control_sy" in roles:
                contrasts.append(("sy", "control_sy"))
            # control_cas = the cascade: the nearest pair in time that is
            # actually dry. One rule, applied to every unit, which is what the
            # paper reports; the other arms remain as the record of how it was
            # arrived at.
            if "control_cas" in roles:
                contrasts.append(("cas", "control_cas"))

            for tag, ctrl_role in contrasts:
                arr, ok, ref = {}, True, None
                for r in sub[sub.role.isin(["case", ctrl_role])].itertuples():
                    f = G / "pairs" / f"{r.ref_date}_{r.sec_date}" / "interferograms" / \
                        f"{r.ref_date}_{r.sec_date}.int.cor.tif"
                    if not f.exists():
                        ok = False
                        break
                    with rasterio.open(f) as src:
                        arr["case" if r.role == "case" else "control"] = \
                            src.read(1).astype("float32")
                    ref = f
                if not ok or len(arr) != 2 or arr["case"].shape != arr["control"].shape:
                    continue
                if lc is None:
                    lc = worldcover_for(G, ref, a.project)
                if lc is None or lc.shape != arr["case"].shape:
                    break
                c, k = arr["case"], arr["control"]
                good = (np.isfinite(c) & np.isfinite(k) & (c > 0) & (c <= 1)
                        & (k > 0) & (k <= 1))
                m0 = sub.iloc[0]
                for code, name in WC.items():
                    if code in DROP:
                        continue
                    sel = good & (lc == code)
                    n = int(sel.sum())
                    if n < MIN_PX:
                        continue
                    cc, kk = c[sel], k[sel]
                    rows.append(dict(
                        unit_id=uid, contrast=tag, group=G.name, event=m0.event,
                        cover=name, koppen=m0.koppen,
                        magnitude=float(m0.magnitude),
                        btemp_days=int(m0.btemp_days), n_px=n,
                        g_case=float(cc.mean()), g_ctrl=float(kk.mean()),
                        dgamma=float((cc - kk).mean()),
                        dgamma_med=float(np.median(cc - kk))))
        if i % 100 == 0:
            log.info("  %s / %s groups | %s rows", f"{i:,}", f"{len(groups):,}", f"{len(rows):,}")

    out = a.out or (a.root / "plan" / f"dgamma_shard{a.shard:02d}.parquet")
    pd.DataFrame(rows).to_parquet(out, index=False)
    log.info("wrote %s rows -> %s", f"{len(rows):,}", out)


if __name__ == "__main__":
    main()
