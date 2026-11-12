#!/usr/bin/env python
"""build_timeline.py -- the Sentinel-1 IW SLC acquisition timeline, 2014 -> today.

Step 1 of the extreme-weather coherence campaign.  Pulls scene-level metadata
for every Sentinel-1 IW SLC ever acquired and writes one parquet per month.

Why SCENE level and not BURST level
-----------------------------------
The archive is ~2.5 M IW SLC scenes, i.e. ~68 M bursts.  Every burst in a scene
is acquired inside the scene's own ~33 s (startTime 23:57:45 -> stopTime
23:58:18 on a measured sample), so burst time differs from scene time by less
than the sampling interval of ANY weather product we can match against -- IMERG
is 30 min, ERA5-Land is 60 min.  Pulling 68 M burst records would cost 27x the
API traffic and buy no temporal precision.

Burst geometry comes from the static OPERA burst DB instead (1,127,661 land
bursts, ~51 MB, already on disk), so the scene -> burst expansion is a local
spatial join.  See build_burst_grid.py.

Being a good API citizen
------------------------
ASF publishes no hard rate limit, so the ceiling here is self-imposed:

  * ONE request in flight at a time -- never parallelised;
  * day-sized chunks (~600-1,100 records, measured 1.9 s each) so no single
    request is huge and any failure costs one day, not one month;
  * a configurable pause between requests (default 0.5 s -> ~2.4 s/day);
  * exponential backoff with jitter on failure (5/15/45/120 s, 5 tries);
  * a day that still fails is recorded in failures.json and SKIPPED, so one bad
    day never aborts a 3 h run.

Total: ~4,350 days at ~2.4 s = ~2.9 h wall, ~2.5 M records.

Resume
------
Output is one parquet per month.  A month whose file exists is skipped whole,
so re-running after an interruption costs at most the partial month (~30
requests).  Re-running a finished archive exits in seconds.

Usage:
    python build_timeline.py --out /smith-scratch/PROJECTS/InSAR/extreme_wx/timeline
    python build_timeline.py --start 2016-01 --end 2016-03 --out ... --dry-run
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import logging
import random
import sys
import time
from pathlib import Path

import asf_search as asf
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("timeline")

# S1A launched 2014-04-03; routine IW tasking ramps through that year.  A test
# query over the Balkans in 2014-05 returned 0 bursts, so the usable record
# really begins later -- we still pull from 2014-10 and let the data say where.
DEFAULT_START = "2014-10"

KEEP = ["sceneName", "startTime", "stopTime", "pathNumber", "orbit",
        "flightDirection", "polarization", "platform", "centerLat", "centerLon",
        "frameNumber", "fileID"]

BACKOFF = [5, 15, 45, 120, 300]


def months(start: str, end: str):
    y, m = (int(x) for x in start.split("-"))
    ey, em = (int(x) for x in end.split("-"))
    while (y, m) <= (ey, em):
        yield y, m
        y, m = (y + 1, 1) if m == 12 else (y, m + 1)


def fetch_day(day: dt.date, sleep: float) -> list[dict] | None:
    """One day of global IW SLC metadata, with backoff.

    Returns the records, or None if every retry was exhausted.  An empty list is
    a genuine answer (days early in the mission really do have no acquisitions)
    and must not be confused with a failure.
    """
    for attempt, wait in enumerate([0] + BACKOFF):
        if wait:
            jitter = wait * (0.5 + random.random())
            log.warning("  retry %d for %s in %.0fs", attempt, day, jitter)
            time.sleep(jitter)
        try:
            res = asf.search(
                platform=asf.PLATFORM.SENTINEL1,
                processingLevel="SLC",
                beamMode="IW",
                start=f"{day}T00:00:00Z",
                end=f"{day + dt.timedelta(days=1)}T00:00:00Z",
                maxResults=60000,
            )
            time.sleep(sleep)
            out = []
            for r in res:
                p = r.properties
                rec = {k: p.get(k) for k in KEEP}
                geom = r.geometry or {}
                coords = geom.get("coordinates") or []
                # footprint as WKT-ish bbox; full polygon kept for the burst join
                if coords and coords[0]:
                    lons = [c[0] for c in coords[0]]
                    lats = [c[1] for c in coords[0]]
                    rec |= {"lon_min": min(lons), "lon_max": max(lons),
                            "lat_min": min(lats), "lat_max": max(lats)}
                out.append(rec)
            return out
        except Exception as e:                      # noqa: BLE001 - any API fault
            log.warning("  %s failed (%s): %s", day, type(e).__name__, str(e)[:120])
    return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--start", default=DEFAULT_START)
    ap.add_argument("--end", default=dt.date.today().strftime("%Y-%m"))
    ap.add_argument("--sleep", type=float, default=0.5,
                    help="pause between requests, seconds (politeness knob)")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    scenes_dir = a.out / "scenes"
    scenes_dir.mkdir(parents=True, exist_ok=True)
    fail_path = a.out / "failures.json"
    failures = json.loads(fail_path.read_text()) if fail_path.exists() else []

    todo = list(months(a.start, a.end))
    log.info("months %s .. %s  (%d)", a.start, a.end, len(todo))
    t0, grand = time.time(), 0

    for y, m in todo:
        fp = scenes_dir / f"{y:04d}-{m:02d}.parquet"
        if fp.exists():
            log.info("%04d-%02d  skip (exists)", y, m)
            continue
        d = dt.date(y, m, 1)
        last = dt.date(y + (m == 12), (m % 12) + 1, 1)
        rows = []
        while d < last and d <= dt.date.today():
            got = fetch_day(d, a.sleep)
            if got is None:
                failures.append(str(d))
            else:
                rows.extend(got)
            d += dt.timedelta(days=1)
        if a.dry_run:
            log.info("%04d-%02d  %6d scenes (dry-run, not written)", y, m, len(rows))
            continue
        if rows:
            df = pd.DataFrame(rows)
            df["startTime"] = pd.to_datetime(df["startTime"], format="mixed", utc=True)
            df["stopTime"] = pd.to_datetime(df["stopTime"], format="mixed", utc=True)
            df.to_parquet(fp, index=False)
        else:
            pd.DataFrame(columns=KEEP).to_parquet(fp, index=False)
        grand += len(rows)
        fail_path.write_text(json.dumps(failures, indent=1))
        el = time.time() - t0
        log.info("%04d-%02d  %6d scenes | total %8d | %.1f min elapsed",
                 y, m, len(rows), grand, el / 60)

    log.info("done: %d scenes, %d failed days, %.1f min", grand, len(failures),
             (time.time() - t0) / 60)


if __name__ == "__main__":
    sys.exit(main())
