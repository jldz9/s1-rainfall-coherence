#!/usr/bin/env python
"""stage_groups.py -- materialise each processing group as a parent-style site.

The parent project's ISCE3_Burst pipeline is per SITE, and a site is a directory
holding burst_ids.txt, acq_dates.txt, pair_list.csv and site_manifest.json.  This
writes one such directory per GROUP, so make_pair_workdirs.py, the insarhub
config and the whole estimator chain are reused UNCHANGED rather than
reimplemented.  That reuse is the point: the estimator must be byte-identical
across every pair in the campaign or the coherences are not comparable.

Layout (deliberately mirrors the parent so --workdir just works):

    proc/burst/<group_id>/
        burst_ids.txt      fullBurstIDs, one per line
        acq_dates.txt      YYYYMMDD, deduplicated and sorted
        pair_list.csv      two rows per unit: the CASE pair and its CONTROL
        site_manifest.json aoi + provenance

Each unit contributes TWO rows to pair_list.csv, tagged `role` = case / control.
They are deliberately adjacent and share a unit_id, because the analysis is the
DIFFERENCE between them -- if one is ever dropped the other is worthless, and a
silently half-processed unit is worse than a missing one.

bperp_m is written as NaN, not 0.0.  The parent obtains it from plan_pairs.py
using orbit state vectors; nothing here knows it, and a 0.0 would read as a
measured zero baseline.  It is not used by make_pair_workdirs.py -- only
btemp_days is -- so NaN is honest and harmless.

Usage:
    python stage_groups.py --root /smith-scratch/PROJECTS/InSAR/extreme_wx
"""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("stage")

SEASONS = {12: "DJF", 1: "DJF", 2: "DJF", 3: "MAM", 4: "MAM", 5: "MAM",
           6: "JJA", 7: "JJA", 8: "JJA", 9: "SON", 10: "SON", 11: "SON"}


def window_label(bt: int) -> str:
    for lo, hi in ((0, 12), (12, 24), (24, 36), (36, 48)):
        if lo < bt <= hi:
            return f"{lo}-{hi}"
    return ">48"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--out", default=None, type=Path)
    a = ap.parse_args()

    plan = a.root / "plan"
    man = pd.read_parquet(plan / "manifest_grouped.parquet")
    avail = json.loads((plan / "burst_availability.json").read_text())
    gran = {b: v.get("granules", {}) for b, v in avail.items()}
    pol_of = {b: (v.get("pol_used") or "VV") for b, v in avail.items()}
    out_root = (a.out or (a.root / "proc")) / "burst"
    out_root.mkdir(parents=True, exist_ok=True)

    prof = pd.read_parquet(a.root / "static" / "burst_profile.parquet").set_index("burst_id_jpl")

    n_groups = n_pairs = 0
    missing = 0
    # Split each group by POLARISATION.  A 2 deg cell can straddle the VV/HH
    # acquisition-mode boundary, and one downloader call takes a single
    # --polarization, so a mixed directory would silently fetch half its bursts.
    # Coherence is not comparable across polarisations either, so they must not
    # share a stratum.
    man = man.assign(pol=man.burst_id_jpl.map(lambda b: pol_of.get(b, "VV")))
    staged = []
    for (gid0, pol), sub in man.groupby(["group_id", "pol"]):
        # UPPERCASE: make_pair_workdirs.py upper-cases its --site before looking
        # for the directory (the parent's sites were YAZD, ERZURUM...), so a
        # lowercase dir is simply not found.  Keeping the on-disk name uppercase
        # is cheaper than patching a script the campaign deliberately reuses.
        gid = (gid0 if pol == "VV" else f"{gid0}_{pol}").upper()
        rows, dates, bursts = [], set(), set()
        for u in sub.itertuples():
            g = gran.get(u.burst_id_jpl, {})
            legs = {"case": (u.case_ref, u.case_sec), "control": (u.ctrl_ref, u.ctrl_sec)}
            resolved = {}
            for role, (r, s) in legs.items():
                rd, sd = pd.Timestamp(r).strftime("%Y-%m-%d"), pd.Timestamp(s).strftime("%Y-%m-%d")
                if rd not in g or sd not in g:
                    resolved = None
                    break
                resolved[role] = (g[rd], g[sd], rd, sd)
            if resolved is None:
                missing += 1
                continue                    # never stage half a unit
            uid = f"{u.burst_id_jpl}_{pd.Timestamp(u.case_sec).strftime('%Y%m%d')}"
            for role, (rs, ss, rd, sd) in resolved.items():
                rdc, sdc = rd.replace("-", ""), sd.replace("-", "")
                rows.append(dict(
                    ref_scene=rs, sec_scene=ss, ref_date=rdc, sec_date=sdc,
                    btemp_days=int(u.baseline_days), bperp_m=float("nan"),
                    season=SEASONS[int(sdc[4:6])], window=window_label(int(u.baseline_days)),
                    unit_id=uid, role=role, event=u.event, cover=u.cover,
                    koppen=u.koppen, magnitude=float(u.magnitude),
                    burst_id_jpl=u.burst_id_jpl))
                dates.update({rdc, sdc})
            bursts.add(u.burst_id_jpl)

        if not rows:
            continue
        d = out_root / gid
        d.mkdir(parents=True, exist_ok=True)
        fb = sorted(b[1:].upper() for b in bursts)
        (d / "burst_ids.txt").write_text("\n".join(fb) + "\n")
        (d / "acq_dates.txt").write_text("\n".join(sorted(dates)) + "\n")
        pd.DataFrame(rows).to_csv(d / "pair_list.csv", index=False)
        lon = prof.loc[list(bursts), "lon_c"]
        lat = prof.loc[list(bursts), "lat_c"]
        # s1_path is the relative orbit as a PLAIN INTEGER (57, not 057).  The
        # downloader's --stacks takes PATH:BURST_ID tokens and rejects a bare
        # burst id outright, so this is not optional metadata.
        track = int(sorted(bursts)[0][1:4])
        json.dump({"group_id": gid, "polarization": pol, "s1_path": track,
                   "aoi": [float(lon.min() - 0.5), float(lat.min() - 0.5),
                           float(lon.max() + 0.5), float(lat.max() + 0.5)],
                   "bursts": fb, "n_units": len(rows) // 2,
                   "task_id": int(sub.task_id.iloc[0])},
                  open(d / "site_manifest.json", "w"), indent=1)
        staged.append({"group_id": gid, "polarization": pol,
                       "task_id": int(sub.task_id.iloc[0]),
                       "n_units": len(rows) // 2})
        n_groups += 1
        n_pairs += len(rows)

    log.info("staged %s groups | %s pairs (%s units) | %s units dropped "
             "for a missing acquisition",
             f"{n_groups:,}", f"{n_pairs:,}", f"{n_pairs // 2:,}", f"{missing:,}")
    st = pd.DataFrame(staged)
    st.to_parquet(plan / "staged_groups.parquet", index=False)
    log.info("polarisation split: %s", st.polarization.value_counts().to_dict())
    log.info("root: %s  (index: %s)", out_root, plan / "staged_groups.parquet")


if __name__ == "__main__":
    main()
