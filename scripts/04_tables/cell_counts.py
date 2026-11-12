#!/usr/bin/env python
"""cell_counts.py -- (event x land cover x Koppen) case/control availability.

The number that decides the campaign.  A cell is usable only if enough distinct
bursts in it have BOTH an event acquisition and a quiet acquisition of the SAME
burst in another year -- the case-control unit.  Counting events is not enough;
an event with no matched control is not a measurement.

Cover assignment
----------------
A burst is 80 x 20 km and is never one cover class, so profile_bursts.py stores
FRACTIONS.  Here a burst is assigned to a cover only when that cover exceeds
`--cover-purity` (default 0.5) of its area; otherwise it is "mixed" and reported
separately rather than silently forced into an argmax bin.  Mixed bursts are not
waste -- they are where the within-burst cover contrast lives -- but they cannot
label a stratum.

Koppen is reported at GROUP level (A/B/C/D/E), not the 30-class legend: 7 events
x 30 classes x 9 covers is 1,890 cells and guarantees empties.  Groups keep the
climate contrast the parent project's data hinted at (rain on arid bare ground
-0.123 against -0.013 temperate) while leaving cells populated.

Usage:
    python cell_counts.py --root /smith-scratch/PROJECTS/InSAR/extreme_wx
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from cover import COVERS, NODATA, assign_cover


def event_masks(d: pd.DataFrame) -> dict[str, np.ndarray]:
    def v(col, scale):
        if col not in d:
            return np.full(len(d), np.nan, "float32")
        x = d[col].to_numpy().astype("float32")
        x[x <= NODATA] = np.nan
        return x / scale
    rain = np.fmax(v("era5_precip_x", 10), v("chirps_precip_x", 10))
    snowf = v("era5_snowfall_x", 100)
    dsnow = v("era5_snowcov_hi_x", 10) - v("era5_snowcov_lo_x", 10)
    tmin = v("era5_tmin_c", 1) / 10 + 200.0
    tmax = v("era5_tmax_c", 1) / 10 + 200.0
    return {
        "rain>50mm":      np.nan_to_num(rain) > 50,
        "rain>100mm":     np.nan_to_num(rain) > 100,
        "snowfall>10mm":  np.nan_to_num(snowf) > 10,
        "snowswing>50%":  np.nan_to_num(dsnow) > 50,
        "freeze-thaw":    (np.nan_to_num(tmin, nan=999) < 273.15)
                          & (np.nan_to_num(tmax, nan=-999) > 273.15),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--cover-purity", type=float, default=0.5)
    a = ap.parse_args()

    prof = pd.read_parquet(a.root / "static" / "burst_profile.parquet")
    prof = prof[prof.koppen > 0].copy()
    prof["cover"] = assign_cover(prof, a.cover_purity)
    pmap = prof.set_index("burst_id_jpl")[["cover", "koppen_group"]].to_dict("index")
    print(f"land bursts profiled: {len(prof):,}  "
          f"(labelled, not mixed: {100*(prof.cover != 'mixed').mean():.0f}%)")

    files = sorted((a.root / "affected").glob("*.parquet"))
    print(f"affected months: {len(files)}  ({files[0].stem} .. {files[-1].stem})")

    ev_years: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    quiet_years: dict[str, int] = defaultdict(int)
    cols = ["burst_id_jpl", "acq_time", "era5_precip_x", "chirps_precip_x",
            "era5_snowfall_x", "era5_snowcov_lo_x", "era5_snowcov_hi_x",
            "era5_tmin_c", "era5_tmax_c"]
    for f in files:
        d = pd.read_parquet(f, columns=cols)
        yr = (pd.to_datetime(d.acq_time).dt.year.to_numpy() - 2014).astype("int8")
        bid = d.burst_id_jpl.to_numpy()
        m = event_masks(d)
        anyev = np.zeros(len(d), bool)
        for name, mask in m.items():
            anyev |= mask
            tgt = ev_years[name]
            for b, y in zip(bid[mask], yr[mask]):
                tgt[b] |= 1 << int(y)
        for b, y in zip(bid[~anyev], yr[~anyev]):
            quiet_years[b] |= 1 << int(y)
        del d

    groups = ["A", "B", "C", "D", "E"]
    for name in m:
        tbl = defaultdict(int)
        for b, ym in ev_years[name].items():
            p = pmap.get(b)
            if p is None or not (quiet_years.get(b, 0) & ~ym):
                continue
            tbl[(p["cover"], p["koppen_group"])] += 1
        covers = [c for c in COVERS + ["mixed"]
                  if any(tbl.get((c, g), 0) for g in groups)]
        if not covers:
            continue
        print(f"\n=== {name} — bursts with BOTH case and control ===")
        print(f"{'cover':8s}" + "".join(f"{g:>9s}" for g in groups) + f"{'total':>10s}")
        for c in covers:
            row = [tbl.get((c, g), 0) for g in groups]
            print(f"{c:8s}" + "".join(f"{v:>9,}" for v in row) + f"{sum(row):>10,}")
        tot_row = [sum(tbl.get((c, g), 0) for c in covers) for g in groups]
        print(f"{'TOTAL':8s}" + "".join(f"{v:>9,}" for v in tot_row) + f"{sum(tot_row):>10,}")


if __name__ == "__main__":
    main()
