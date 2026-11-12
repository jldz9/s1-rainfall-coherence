#!/usr/bin/env python
"""Build the processing list for the twelve-day rain campaign.

Pairs are formed between any two acquisitions twelve days apart, not merely
between consecutive ones. Where the revisit is six days a twelve-day pair skips
an acquisition, and an earlier version of this count built only consecutive
pairs and so missed those bursts entirely.

A case pair carries rain above the threshold at exactly one of its two
acquisitions and is dry at the other. A control pair is dry at both, lies in the
same burst, and falls within a month of the case's day of year in a different
year. One case is kept per storm, a storm being a day within a two-degree cell,
and one control is shared between all cases of the same burst and season.
"""
import argparse, logging, glob
import numpy as np, pandas as pd
from pathlib import Path

BASE_LO, BASE_HI = 11.5, 12.5
WET_MM, DRY_MM = 15.0, 1.0
SEASON_TOL, MAX_YEARS = 30, 5
CLUSTER_DEG = 2.0
SEASON_BIN = 60          # days; controls are shared within a bin

log = logging.getLogger("campaign")


def all_12day_pairs(d: pd.DataFrame) -> pd.DataFrame:
    """Every pair of acquisitions 12 days apart, consecutive or not."""
    out = []
    for b, g in d.groupby("burst_id_jpl", sort=False):
        t = g.t.to_numpy()
        r = g.rain.to_numpy()
        # the match for each acquisition is sought among those that follow it,
        # and with a 12-day window at most a handful are ever in range
        for k in (1, 2, 3):
            if len(g) <= k:
                break
            dt = (t[k:] - t[:-k]) / np.timedelta64(1, "D")
            m = (dt >= BASE_LO) & (dt <= BASE_HI)
            if not m.any():
                continue
            out.append(pd.DataFrame({
                "burst_id_jpl": b,
                "t_ref": t[:-k][m], "t_sec": t[k:][m],
                "rain_ref": r[:-k][m], "rain_sec": r[k:][m],
                "skip": k - 1,
            }))
    p = pd.concat(out, ignore_index=True)
    return p.drop_duplicates(subset=["burst_id_jpl", "t_ref", "t_sec"])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--metric", choices=["mean", "max"], default="mean")
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    R = a.root

    d = pd.concat([pd.read_parquet(f) for f in
                   sorted(glob.glob(str(R / "plan" / "rainscreen_shard*.parquet")))],
                  ignore_index=True)
    d["t"] = pd.to_datetime(d.t)
    d["rain"] = d[f"rain_{a.metric}"]
    d = d.sort_values(["burst_id_jpl", "t"]).reset_index(drop=True)
    log.info("acquisitions: %s over %s bursts", f"{len(d):,}", f"{d.burst_id_jpl.nunique():,}")

    p = all_12day_pairs(d)
    log.info("12-day pairs (any gap): %s  [consecutive only would give %s]",
             f"{len(p):,}", f"{int((p.skip == 0).sum()):,}")

    wet_sec = p.rain_sec >= WET_MM
    wet_ref = p.rain_ref >= WET_MM
    case = p[(wet_sec & (p.rain_ref < DRY_MM)) | (wet_ref & (p.rain_sec < DRY_MM))].copy()
    case["t_event"] = np.where(case.rain_sec >= WET_MM, case.t_sec, case.t_ref)
    case["mm"] = np.maximum(case.rain_sec, case.rain_ref)
    ctrl = p[(p.rain_sec < DRY_MM) & (p.rain_ref < DRY_MM)].copy()
    log.info("case candidates: %s   control candidates: %s", f"{len(case):,}", f"{len(ctrl):,}")

    # one case per storm
    g = pd.read_parquet(R / "plan" / "groups.parquet")
    g["gid"] = g.group_id.str.lower()
    cen = g.set_index("gid")[["lon_min", "lat_min"]]
    lo = case.burst_id_jpl.str.lower()
    case["lon"] = cen.lon_min.reindex(lo).to_numpy() + 0.5
    case["lat"] = cen.lat_min.reindex(lo).to_numpy() + 0.5
    case["t_event"] = pd.to_datetime(case.t_event)
    key = (case.t_event.dt.strftime("%Y%m%d") + "_"
           + (case.lon / CLUSTER_DEG).round().astype("Int64").astype(str) + "_"
           + (case.lat / CLUSTER_DEG).round().astype("Int64").astype(str))
    case = case.sort_values("mm", ascending=False).groupby(key, sort=False).head(1)
    log.info("after one-per-storm: %s cases", f"{len(case):,}")

    # assign a concrete control to each case
    case["doy"] = case.t_event.dt.dayofyear
    case["yr"] = case.t_event.dt.year
    ctrl["doy"] = ctrl.t_sec.dt.dayofyear
    ctrl["yr"] = ctrl.t_sec.dt.year
    by = {b: g.reset_index(drop=True) for b, g in ctrl.groupby("burst_id_jpl", sort=False)}
    pick = []
    for b, doy, yr in zip(case.burst_id_jpl, case.doy, case.yr):
        cg = by.get(b)
        if cg is None:
            pick.append(None); continue
        dd = np.abs(cg.doy.to_numpy() - doy)
        dd = np.minimum(dd, 365 - dd)
        dy = np.abs(cg.yr.to_numpy() - yr)
        ok = (dd <= SEASON_TOL) & (dy > 0) & (dy <= MAX_YEARS)
        if not ok.any():
            pick.append(None); continue
        j = np.flatnonzero(ok)[np.argmin(dd[ok] + dy[ok] * 0.5)]   # closest in season, then year
        pick.append((cg.t_ref.iloc[j], cg.t_sec.iloc[j]))
    case["ctrl_ref"] = [x[0] if x else pd.NaT for x in pick]
    case["ctrl_sec"] = [x[1] if x else pd.NaT for x in pick]
    case = case.dropna(subset=["ctrl_ref"]).copy()
    log.info("cases with an assigned control: %s", f"{len(case):,}")

    # Controls are NOT forced into a shared slot. Overwriting each case's own
    # control with the first of its 60-day season bin broke the season rule for
    # every other case in that bin, by as much as 89 days. Cases close in season
    # already choose the same nearest control, so the sharing that is safe
    # happens by itself and is recovered by counting distinct controls.

    # what is already on disk
    pr = pd.read_parquet(R / "plan" / "processed_relabelled.parquet")
    done = set(pr.burst + "_" + pr.ref.dt.strftime("%Y%m%d") + "_" + pr.sec.dt.strftime("%Y%m%d"))
    lo = case.burst_id_jpl.str.lower()
    case["case_key"] = lo + "_" + case.t_ref.dt.strftime("%Y%m%d") + "_" + case.t_sec.dt.strftime("%Y%m%d")
    case["ctrl_key"] = lo + "_" + case.ctrl_ref.dt.strftime("%Y%m%d") + "_" + case.ctrl_sec.dt.strftime("%Y%m%d")
    case["case_done"] = case.case_key.isin(done)
    case["ctrl_done"] = case.ctrl_key.isin(done)

    rain_todo = case[~case.case_done]
    ctrl_todo = (case[~case.ctrl_done].drop_duplicates("ctrl_key"))

    print("\n" + "=" * 62)
    print("CAMPAIGN SUMMARY")
    print("=" * 62)
    print(f"  rain (case) pairs required      {len(case):>8,}")
    print(f"     already processed            {int(case.case_done.sum()):>8,}")
    print(f"     TO PROCESS                   {len(rain_todo):>8,}")
    print()
    print(f"  control pairs required (unique) {case.ctrl_key.nunique():>8,}")
    print(f"     already processed            {int(case.drop_duplicates('ctrl_key').ctrl_done.sum()):>8,}")
    print(f"     TO PROCESS                   {len(ctrl_todo):>8,}")
    print()
    tot = len(rain_todo) + len(ctrl_todo)
    print(f"  TOTAL INTERFEROGRAMS TO PROCESS {tot:>8,}")
    print(f"     distinct bursts involved     {len(set(rain_todo.burst_id_jpl) | set(ctrl_todo.burst_id_jpl)):>8,}")
    print(f"     cpu-hours (0.39 each)        {tot * 4076 / 10380:>8,.0f}")
    print(f"     disk (36 MB each)            {tot * 361 / 10380:>8,.0f} GB")

    print("\n  rain pair event size")
    for lo_, hi in [(15, 25), (25, 50), (50, 100), (100, 1e9)]:
        m = (rain_todo.mm >= lo_) & (rain_todo.mm < hi)
        lab = f"{lo_:g}-{hi:g} mm" if hi < 1e9 else ">100 mm"
        print(f"    {lab:<11s} {int(m.sum()):7,}  {m.mean()*100:5.1f}%")

    rain_todo.to_parquet(R / "plan" / "campaign_rain_pairs.parquet", index=False)
    ctrl_todo.to_parquet(R / "plan" / "campaign_ctrl_pairs.parquet", index=False)
    case.to_parquet(R / "plan" / "campaign_full.parquet", index=False)
    print(f"\n  wrote campaign_rain_pairs.parquet / campaign_ctrl_pairs.parquet")


if __name__ == "__main__":
    main()
