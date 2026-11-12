#!/usr/bin/env python
"""verify_availability.py -- confirm every planned acquisition really exists.

The gate between the candidate manifest and processing.  Until this runs, the
manifest is a plan drawn on the OPERA burst grid and the weather rasters; neither
knows whether ASF actually holds an SLC-BURST product for a given burst on a given
date.

This is also where the join's deliberate over-assignment gets filtered out.
mark_affected.py assigns ~45% more bursts to a scene than were acquired in it
(46 against ~31 measured on one scene), because the burst bbox is the hull of a
rotated footprint plus a 5 km margin.  That was left in on purpose -- it is a
candidate generator -- and this is the step that was always going to correct it.

Query shape
-----------
ASF returns no frameNumber on BURST products, so a burst stack is keyed by
fullBurstID, which is burst_id_jpl with the leading "t" dropped and upper-cased:

    t018_038640_iw2  ->  018_038640_IW2      (841 products, 437 dates, measured)

One query per DISTINCT burst returns its whole history, so the four dates a unit
needs are then checked locally.  3,121 bursts, not 13,112 burst-dates.

Politeness: one request in flight, a pause between, exponential backoff with
jitter, and a per-burst cache so a re-run costs nothing.  Compute nodes were
confirmed to reach ASF (job 603918), so this can run under SLURM.

Outputs
-------
    plan/burst_availability.parquet   burst -> available dates, polarizations
    plan/manifest_verified.parquet    units whose four acquisitions all exist

Usage:
    python verify_availability.py --root /smith-scratch/PROJECTS/InSAR/extreme_wx
"""
from __future__ import annotations

import argparse
import json
import logging
import random
import time
from pathlib import Path

import asf_search as asf
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("verify")

BACKOFF = [5, 15, 45, 120]
CHECKPOINT_EVERY = 100

# Retry only on a COMPLETELY empty answer.
#
# An earlier version set this to 40 on the reasoning that "the plan only picked
# bursts carrying both a case and a control, so a short history means a truncated
# response".  That reasoning is wrong: a case plus a control is FOUR
# acquisitions, not forty, so a sparsely-tasked burst legitimately has a short
# history.  Measured: t012_023918_iw1 returns 13 dates on every one of five
# retries -- stable, therefore real -- and 7 of the first 500 bursts have fewer
# than 40 dates.  Each one burned the full 5/15/45/120 s backoff for nothing,
# which is what took the job from 2.6 s/burst to 13 s/burst.
#
# Truncation needs no floor of its own: a truncated response simply fails the
# four-date lookup later and drops the unit.  That errs toward losing a unit
# rather than inventing one, which is the right direction.
MIN_PLAUSIBLE_DATES = 1


def fbid(burst_id_jpl: str) -> str:
    """t018_038640_iw2 -> 018_038640_IW2"""
    return burst_id_jpl[1:].upper()


def fetch_burst(b: str, sleep: float) -> dict | None:
    for attempt, wait in enumerate([0] + BACKOFF):
        if wait:
            j = wait * (0.5 + random.random())
            log.warning("  retry %d for %s in %.0fs", attempt, b, j)
            time.sleep(j)
        try:
            r = asf.search(platform=asf.PLATFORM.SENTINEL1, processingLevel="BURST",
                           fullBurstID=fbid(b), maxResults=3000)
            # Keep the GRANULE NAME per date, not just the date.  The parent
            # pipeline's pair_list.csv is keyed by ASF granule
            # (S1_120480_IW2_20241227T141922_VV_4379-BURST), so discarding names
            # here would force a second full pass over ASF later.
            #
            # Group by polarization rather than filtering to VV.  SLC-BURST
            # granules are single-polarisation, and Sentinel-1 images the high
            # latitudes in HH+HV, not VV+VH: burst 010_020042_IW3 has 983
            # products, every one of them HH or HV, and none VV.  An earlier
            # version filtered to VV, got an empty result, read that as a
            # truncated response and burned the whole 5/15/45/120 s backoff
            # chain on it -- 74 bursts had done that before it was caught, and
            # 382 of the planned units sit above 60 degrees.
            by_pol: dict[str, dict[str, str]] = {}
            for p in r:
                pr = p.properties
                pol = (pr.get("polarization") or "").upper()
                if not pol:
                    continue
                day = str(pr["startTime"])[:10]
                by_pol.setdefault(pol, {}).setdefault(day, pr.get("fileID") or pr.get("sceneName"))
            pols = sorted(by_pol)

            # Prefer VV; fall back to HH where the mission flies HH.  Never mix
            # the two inside one burst -- a case/control pair must be same-pol,
            # since coherence magnitude is not comparable across polarisations.
            pol_used = "VV" if "VV" in by_pol else ("HH" if "HH" in by_pol else
                                                    (pols[0] if pols else ""))
            by_date = by_pol.get(pol_used, {})
            dates = set(by_date)

            # Retry only on a genuinely EMPTY or truncated answer, judged on the
            # whole product list rather than on one polarisation.
            all_days = {d for m in by_pol.values() for d in m}
            if len(all_days) < MIN_PLAUSIBLE_DATES:
                log.warning("  %s: ASF returned nothing at all, retrying", b)
                continue
            return {"burst_id_jpl": b, "n_products": len(r),
                    "dates": sorted(dates), "pols": pols,
                    "pol_used": pol_used, "granules": by_date}
        except Exception as e:                                   # noqa: BLE001
            log.warning("  %s failed (%s): %s", b, type(e).__name__, str(e)[:120])
        finally:
            time.sleep(sleep)
    return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--sleep", type=float, default=0.3)
    a = ap.parse_args()

    man = pd.read_parquet(a.root / "plan" / "manifest.parquet")
    bursts = sorted(man.burst_id_jpl.unique())
    log.info("units %s | distinct bursts %s", f"{len(man):,}", f"{len(bursts):,}")

    cache_fp = a.root / "plan" / "burst_availability.json"
    cache: dict = json.loads(cache_fp.read_text()) if cache_fp.exists() else {}
    # Drop cached FAILURES so they are retried.  A failure is cached like any
    # other result, which is right for resuming across a restart but wrong after
    # an outage: the CMR outage on 2026-09-23 wrote 8 bursts as n_products=-1,
    # and without this they would be skipped forever and their units silently
    # dropped from the campaign.  A real "this burst does not exist" costs one
    # cheap re-query to rediscover; a false failure costs units.
    stale = [b for b, v in cache.items() if v.get("n_products", -1) <= 0]
    for b in stale:
        del cache[b]
    log.info("cached bursts: %s  (dropped %d cached failures for retry)",
             f"{len(cache):,}", len(stale))

    t0 = time.time()
    todo = [b for b in bursts if b not in cache]
    for n, b in enumerate(todo, 1):
        got = fetch_burst(b, a.sleep)
        cache[b] = got or {"burst_id_jpl": b, "n_products": -1, "dates": [],
                           "pols": [], "pol_used": "", "granules": {}}
        if n % CHECKPOINT_EVERY == 0:
            cache_fp.write_text(json.dumps(cache))
            el = (time.time() - t0) / 60
            log.info("%s / %s | %.1f min | eta %.1f h",
                     f"{n:,}", f"{len(todo):,}", el, el / n * (len(todo) - n) / 60)
    cache_fp.write_text(json.dumps(cache))

    av = pd.DataFrame([{"burst_id_jpl": k, "n_products": v["n_products"],
                        "n_dates": len(v["dates"]), "pols": ",".join(v["pols"])}
                       for k, v in cache.items()])
    av.to_parquet(a.root / "plan" / "burst_availability.parquet", index=False)

    have = {k: set(v["dates"]) for k, v in cache.items()}
    cols = ["case_ref", "case_sec", "ctrl_ref", "ctrl_sec"]
    ok = pd.Series(True, index=man.index)
    per_col = {}
    for c in cols:
        d = pd.to_datetime(man[c]).dt.strftime("%Y-%m-%d")
        hit = [dt in have.get(b, set()) for b, dt in zip(man.burst_id_jpl, d)]
        per_col[c] = sum(hit)
        ok &= pd.Series(hit, index=man.index)

    man["available"] = ok
    man[ok].to_parquet(a.root / "plan" / "manifest_verified.parquet", index=False)

    log.info("per-acquisition hit rate:")
    for c in cols:
        log.info("   %-9s %s / %s  (%.1f%%)", c, f"{per_col[c]:,}", f"{len(man):,}",
                 100 * per_col[c] / len(man))
    from collections import Counter
    log.info("polarisation used per burst: %s",
             dict(Counter(v.get("pol_used", "") for v in cache.values())))
    log.info("units with ALL FOUR acquisitions: %s / %s (%.1f%%)",
             f"{int(ok.sum()):,}", f"{len(man):,}", 100 * ok.mean())
    log.info("interferograms to process: %s", f"{2*int(ok.sum()):,}")
    if int(ok.sum()):
        print("\nverified units by event x cover:")
        print(pd.crosstab(man[ok].cover, man[ok].event).to_string())


if __name__ == "__main__":
    main()
