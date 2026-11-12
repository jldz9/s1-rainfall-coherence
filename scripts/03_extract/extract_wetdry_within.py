#!/usr/bin/env python
"""Coherence loss on rained-on ground against un-rained-on ground in one scene.

Every comparison in the campaign so far sets one interferogram against another,
and must argue that the two are alike.  This one does not.  An IMERG cell is
about 11 km and a burst is 80 by 20 km, so a storm covering part of a burst
leaves wet and dry ground inside the same image -- same orbit, same date, same
processing, same control pair.  Nothing needs matching because nothing differs.

The loss is formed separately over the two populations, each against its own
part of the control, so that whatever static difference exists between the
ground that catches rain and the ground beside it is divided out before the two
are compared.

This can only be done where the storm covers part of the burst and not all of
it.  Above about 7.5 mm coverage reaches 1.00 and no dry ground remains, so the
test is confined to light rainfall by construction -- which is the regime the
dose curve makes the strongest claim about.

Usage:
    python extract_wetdry_within.py --root <wx root> --shard N --of M
"""
from __future__ import annotations

import argparse, logging, sys, time
from pathlib import Path
import numpy as np, pandas as pd, rasterio

sys.path.insert(0, str(Path(__file__).resolve().parent))
from extract_rainmasked import rain_mask_for, cor_path    # noqa: E402

log = logging.getLogger("wetdry")
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
MIN_PX = 200
PX_M2 = 1740.0          # 42 m posting of the coherence grid


def cell_px(lat: float) -> float:
    """One IMERG cell, 0.1 deg square, in coherence pixels at this latitude.

    The rain field knows wet from dry at 11 km and no finer, so a population
    smaller than one cell is the burst outline clipping through a cell rather
    than an observation of ground that did or did not receive rain.
    """
    km2 = (0.1 * 110.574) * (0.1 * 111.320 * np.cos(np.radians(lat)))
    return km2 * 1e6 / PX_M2
PCT = (10, 25, 50, 75, 90)


def stats(gc, gk, sel, tag):
    cc, kk = gc[sel], gk[sel]
    d = cc - kk
    out = {f"g_case_{tag}": float(cc.mean()), f"g_ctrl_{tag}": float(kk.mean()),
           f"n_px_{tag}": int(sel.sum()), f"dmed_{tag}": float(np.median(d))}
    for p, v in zip(PCT, np.percentile(d, PCT)):
        out[f"d{p}_{tag}"] = float(v)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--of", type=int, default=1)
    ap.add_argument("--lo", type=float, default=0.2)
    ap.add_argument("--hi", type=float, default=0.8)
    ap.add_argument("--project", default="ifgqual")
    ap.add_argument("--min-cells", type=float, default=1.0,
                    help="each population must span at least this many IMERG "
                         "cells; 1.0 is one full cell, the resolution of the "
                         "rain field")
    ap.add_argument("--tag", default="",
                    help="suffix for the output, so a rerun with different "
                         "coverage limits does not overwrite an earlier one")
    a = ap.parse_args()
    R = a.root

    d = pd.read_parquet(R / "plan" / "units_final.parquet")
    d = d[d.frac_gt1.between(a.lo, a.hi)].copy()
    # skip the rasters of units that cannot clear the cell rule, from the
    # footprint fraction alone; the exact test is made on the pixels below
    need = a.min_cells * d.lat.map(cell_px)
    d = d[np.minimum(d.frac_gt1, 1 - d.frac_gt1) * d.n_px >= 0.8 * need].copy()
    bycase = {}
    for g in (R / "proc" / "burst").glob("*"):
        if (g / ".done").exists():
            bycase.setdefault(g.name.rsplit("_", 1)[0].lower(), []).append(g)
    d = d.iloc[a.shard::a.of].reset_index(drop=True)
    log.info("units in shard: %s (coverage %.1f-%.1f, >= %.1f IMERG cell each side)",
             f"{len(d):,}", a.lo, a.hi, a.min_cells)

    rows, t0 = [], time.time()
    for i, r in enumerate(d.itertuples(), 1):
        hit = None
        for g in bycase.get(r.burst, []):
            cc = cor_path(R, g.name, pd.Timestamp(r.t_ref).strftime("%Y%m%d"),
                          pd.Timestamp(r.t_sec).strftime("%Y%m%d"))
            ck = cor_path(R, g.name, pd.Timestamp(r.ctrl_ref).strftime("%Y%m%d"),
                          pd.Timestamp(r.ctrl_sec).strftime("%Y%m%d"))
            if cc and ck:
                hit = (g, cc, ck)
                break
        if hit is None:
            continue
        g, c_case, c_ctrl = hit
        m = rain_mask_for(g, c_case, r.t_event, a.project)
        if m is None or not m.any() or m.all():
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
        lcf = g / "worldcover.npy"
        lc = np.load(lcf) if lcf.exists() else None
        land = (lc != 80) if lc is not None else np.ones_like(m, bool)
        ok = land & np.isfinite(gc) & np.isfinite(gk) & (gc > 0) & (gc <= 1) \
             & (gk > 0) & (gk <= 1)
        wet, dry = ok & m, ok & ~m
        need = max(MIN_PX, a.min_cells * cell_px(r.lat))
        if wet.sum() < need or dry.sum() < need:
            continue
        row = {"unit_id": r.unit_id, "burst": r.burst, "rain_mm": r.rain_mm,
               "frac_gt1": r.frac_gt1, "koppen": r.koppen}
        row.update(stats(gc, gk, wet, "wet"))
        row.update(stats(gc, gk, dry, "dry"))
        # does the cover composition differ between the two populations?
        if lc is not None:
            for tag, sel in (("wet", wet), ("dry", dry)):
                vals, cnt = np.unique(lc[sel], return_counts=True)
                row[f"domcover_{tag}"] = int(vals[cnt.argmax()])
                row[f"domfrac_{tag}"] = float(cnt.max() / cnt.sum())
        rows.append(row)
        if i % 50 == 0:
            el = time.time() - t0
            log.info("  %s / %s  kept %s  eta %.0f min", f"{i:,}", f"{len(d):,}",
                     f"{len(rows):,}", (len(d) - i) * el / i / 60)

    out = R / "plan" / f"wetdry{a.tag}_shard{a.shard:03d}.parquet"
    pd.DataFrame(rows).to_parquet(out, index=False)
    log.info("wrote %s (%s rows)", out, f"{len(rows):,}")


if __name__ == "__main__":
    main()
