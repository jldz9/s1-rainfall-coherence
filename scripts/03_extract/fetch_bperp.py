#!/usr/bin/env python
"""Perpendicular baselines for every processed pair, from the ASF baseline service.

ASF computes baselines for SLC scenes but not for BURST products, which is why
the bperp_m column of pair_list.csv was never filled. The parent scene is
recovered from the burst granule name, whose trailing four hex characters are
the scene's product identifier, and one baseline query per stack then returns
B_perp for every acquisition relative to a common reference. The baseline of a
pair is the difference of two such values.

Baselines are scene-level, evaluated at the scene centre; within one burst the
variation is negligible for the present purpose.
"""
import argparse, json, logging, re, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd
import requests

API = "https://api.daac.asf.alaska.edu/services/search"
BURST = re.compile(r"S1_(\d{6})_(IW\d)_(\d{8})T\d{6}_VV_([0-9A-F]{4})-BURST")
log = logging.getLogger("bperp")


def seed_scene(sess, track, lon, lat, date):
    """One parent SLC scene for this burst, used to seed the stack query."""
    r = sess.get(f"{API}/param", params={
        "platform": "SENTINEL-1", "processingLevel": "SLC",
        "relativeOrbit": int(track),
        "start": f"{date[:4]}-{date[4:6]}-{date[6:]}T00:00:00Z",
        "end": f"{date[:4]}-{date[4:6]}-{date[6:]}T23:59:59Z",
        "intersectsWith": f"POINT({lon:.4f} {lat:.4f})",
        "output": "jsonlite"}, timeout=90)
    r.raise_for_status()
    res = r.json().get("results", [])
    return res[0]["granuleName"] if res else None


def stack(sess, ref):
    """date -> perpendicular baseline, for the whole stack of that scene."""
    r = sess.get(f"{API}/baseline", params={"reference": ref, "output": "jsonlite"},
                 timeout=180)
    r.raise_for_status()
    out = {}
    for x in r.json().get("results", []):
        b = x.get("perpendicularBaseline")
        if b is None:
            continue
        out[x["granuleName"][17:25]] = float(b)      # stopTime date of the scene
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")

    P = pd.read_parquet(a.root / "plan" / "pair_scenes.parquet")
    g = pd.read_parquet(a.root / "plan" / "groups.parquet")
    g["b"] = g.group_id.str.lower()
    cen = g.set_index("b")[["lon_min", "lat_min"]]

    # one representative acquisition and centre per burst
    B = (P.groupby("burst_id_jpl")
           .agg(group=("group", "first"), date=("ref_date", "first"))
           .reset_index())
    B["track"] = B.burst_id_jpl.str.extract(r"^t(\d{3})")[0].astype(int)
    B["lon"] = cen.lon_min.reindex(B.burst_id_jpl).to_numpy() + 0.5
    B["lat"] = cen.lat_min.reindex(B.burst_id_jpl).to_numpy() + 0.5
    B = B.dropna(subset=["lon"])
    log.info("bursts to resolve: %s", f"{len(B):,}")

    out = a.root / "plan" / "bperp_by_burst.jsonl"
    done = set()
    if out.exists():
        for line in out.open():
            try: done.add(json.loads(line)["burst"])
            except Exception: pass
        log.info("resuming, already done: %s", f"{len(done):,}")
    B = B[~B.burst_id_jpl.isin(done)]

    fh = out.open("a")
    sess = requests.Session()
    n = [0]

    def one(r):
        for attempt in range(3):
            try:
                ref = seed_scene(sess, r.track, r.lon, r.lat, str(r.date))
                if not ref:
                    return {"burst": r.burst_id_jpl, "err": "no parent scene"}
                return {"burst": r.burst_id_jpl, "ref": ref, "bperp": stack(sess, ref)}
            except Exception as e:
                if attempt == 2:
                    return {"burst": r.burst_id_jpl, "err": str(e)[:120]}
                time.sleep(3 * (attempt + 1))

    with ThreadPoolExecutor(a.workers) as ex:
        for rec in ex.map(one, [r for r in B.itertuples()]):
            fh.write(json.dumps(rec) + "\n"); fh.flush()
            n[0] += 1
            if n[0] % 100 == 0:
                log.info("%s / %s", f"{n[0]:,}", f"{len(B):,}")
    fh.close()
    log.info("done: %s records", f"{n[0]:,}")


if __name__ == "__main__":
    main()
