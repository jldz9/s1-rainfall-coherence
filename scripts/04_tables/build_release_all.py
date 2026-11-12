#!/usr/bin/env python
"""Consolidate every campaign arm into one analysis-ready table.

build_release.py joined the two arms that existed when it was written. The
low-dose arms added since carry the same unit schema, so the arms are listed
here instead of named in the code, and an arm missing a column simply leaves it
empty rather than being excluded.

Soil moisture is present only for the arms whose Earth Engine pull has run; the
moisture analysis drops the rest rather than the whole arm being held back.
"""
import argparse, logging
from pathlib import Path
import numpy as np, pandas as pd

FLOOR = 0.077
log = logging.getLogger("release")

ARMS = [("r12_units_full.parquet",  "campaign_full.parquet", "main"),
        ("L12_units_full.parquet",  "lowdose_full.parquet",  "low-dose"),
        ("M12_units_full.parquet",       "lowdose2_full.parquet", "low-dose-2"),
        ("P12_units_full.parquet",       "lowdose3_full.parquet", "low-dose-3"),
        ("Q12_units_full.parquet",       "lowdose4_full.parquet", "proportional"),
        ("H12_units_full.parquet",       "heavy_full.parquet",    "heavy-top-up")]


def season(doy: float, lat: float) -> str:
    if not (np.isfinite(doy) and np.isfinite(lat)):
        return ""
    d = doy if lat >= 0 else (doy + 182.6) % 365.25
    return ["DJF", "MAM", "JJA", "SON"][int(((d + 30) % 365.25) // 91.3)]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    R = a.root

    parts = []
    for units, plan, arm in ARMS:
        f = R / "plan" / units
        if not f.exists():
            log.warning("arm %s absent (%s)", arm, units); continue
        u = pd.read_parquet(f).copy()
        p = pd.read_parquet(R / "plan" / plan).copy()
        p["t_event"] = pd.to_datetime(p.t_event)
        p["unit_id"] = p.burst_id_jpl + "_" + p.t_event.dt.strftime("%Y%m%d")
        p = p.drop_duplicates("unit_id").set_index("unit_id")
        for c in ("t_ref", "t_sec", "ctrl_ref", "ctrl_sec", "lon", "lat", "doy"):
            if c in p.columns and c not in u.columns:
                u[c] = p[c].reindex(u.unit_id).to_numpy()
        u["t_event"] = p.t_event.reindex(u.unit_id).to_numpy()
        u["arm"] = arm
        parts.append(u)
        log.info("%-14s %7s units", arm, f"{len(u):,}")
    d = pd.concat(parts, ignore_index=True)
    d = d.rename(columns={"mm": "rain_mm", "th_baseline": "theta_baseline",
                          "red": "red_pct"})

    # a burst may appear in more than one arm if the sampling overlapped; the
    # event is the same event, so the first is kept rather than counted twice
    dup = d.unit_id.duplicated().sum()
    if dup:
        log.warning("%s units appear in more than one arm; keeping first", f"{dup:,}")
        d = d.drop_duplicates("unit_id")

    # the event is at whichever acquisition carries the rain, so any quantity
    # differenced across the pair must be oriented to it rather than to sec-ref
    d["rain_at_sec"] = (pd.to_datetime(d.t_event).dt.normalize()
                        == pd.to_datetime(d.t_sec).dt.normalize())
    if {"th_case_ref", "th_case_sec"} <= set(d.columns):
        d["theta_wet"] = np.where(d.rain_at_sec, d.th_case_sec, d.th_case_ref)
        d["theta_dry"] = np.where(d.rain_at_sec, d.th_case_ref, d.th_case_sec)
        d["d_theta"] = d.theta_wet - d.theta_dry
    if "doy" not in d.columns or d.doy.isna().all():
        d["doy"] = pd.to_datetime(d.t_event).dt.dayofyear
    d["season"] = [season(x, y) for x, y in zip(d.doy, d.lat)]
    d["ratio"] = d.g_case / d.g_ctrl
    d["reduction_pct"] = (1 - d.ratio) * 100
    d["headroom"] = d.g_ctrl - FLOOR

    keep = ["unit_id", "arm", "burst", "koppen", "season", "lon", "lat",
            "t_event", "t_ref", "t_sec", "ctrl_ref", "ctrl_sec", "rain_at_sec",
            "rain_mm", "g_case", "g_ctrl", "ratio", "reduction_pct", "headroom",
            "n_px", "dom_cover", "theta_baseline", "theta_wet", "theta_dry", "d_theta"]
    d = d[[c for c in keep if c in d.columns]]
    out = R / "release" / "rain_coherence_units.parquet"
    out.parent.mkdir(exist_ok=True)
    d.to_parquet(out, index=False)
    d.to_csv(out.with_suffix(".csv"), index=False)
    log.info("wrote %s (%s units, %s columns)", out, f"{len(d):,}", len(d.columns))

    covs = []
    for pre in ("r12", "L12", "M12", "P12", "Q12", "H12"):
        f = R / "plan" / f"{pre}_cover.parquet"
        if f.exists():
            covs.append(pd.read_parquet(f))
    cov = pd.concat(covs, ignore_index=True).drop_duplicates(["unit_id", "cover"])
    cov.to_parquet(R / "release" / "rain_coherence_by_cover.parquet", index=False)
    log.info("wrote by-cover table (%s rows)", f"{len(cov):,}")

    ok = d[d.headroom >= 0.15]
    print(f"\n  total units      {len(d):,}")
    print(f"  headroom pass    {len(ok):,} ({len(ok)/len(d)*100:.1f}%)")
    print(f"  mean reduction   {ok.reduction_pct.mean():.2f}%")
    print(f"  with moisture    {int(ok.d_theta.notna().sum()):,}" if "d_theta" in ok else "")


if __name__ == "__main__":
    main()
