#!/usr/bin/env python
"""merge_and_analyse.py -- fold the extraction shards together and report every
control arm separately.

Run after the extraction array finishes. Merging was previously done by hand at
a prompt, which is how plan/dgamma_all.parquet came to hold two contrasts that
the analysis then averaged together without anyone noticing.

Each contrast present is analysed on its own and its output kept, so the
comparison between control-selection rules is a paired result on the same
units rather than a claim.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd

FINAL = Path("/home/C838053462/dev/extreme_wx/scripts/final_analysis.py")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--keep-shards", action="store_true")
    a = ap.parse_args()
    R = a.root
    plan = R / "plan"

    shards = sorted(plan.glob("dgamma_shard*.parquet"))
    if not shards:
        raise SystemExit("no dgamma_shard*.parquet -- run the extraction first")
    parts = [pd.read_parquet(s) for s in shards]
    d = pd.concat(parts, ignore_index=True)
    print(f"merged {len(shards)} shards -> {len(d):,} cover-rows, "
          f"{d.unit_id.nunique():,} units")

    # A unit must not appear twice under one contrast: shards partition the
    # groups, so a duplicate means two shards saw the same group -- which
    # happened to the campaign itself when two submitters overlapped.
    dup = d.duplicated(subset=["unit_id", "contrast", "cover"]).sum()
    if dup:
        print(f"WARNING: {dup:,} duplicated unit/contrast/cover rows; keeping first")
        d = d.drop_duplicates(subset=["unit_id", "contrast", "cover"], keep="first")

    print("\nrows by contrast:")
    print(d.contrast.value_counts().to_string())
    n_arms = d.groupby("unit_id").contrast.nunique()
    print(f"units with more than one arm: {(n_arms > 1).sum():,} of {len(n_arms):,}")

    out = plan / "dgamma_all.parquet"
    d.to_parquet(out, index=False)
    print(f"wrote {out}")

    if not a.keep_shards:
        for s in shards:
            s.unlink()
        print(f"removed {len(shards)} shard files")

    rc = 0
    for c in sorted(d.contrast.unique()):
        print(f"\n{'#' * 30} contrast {c} {'#' * 30}")
        r = subprocess.run([sys.executable, str(FINAL), "--root", str(R),
                            "--contrast", c])
        rc |= r.returncode
    sys.exit(rc)


if __name__ == "__main__":
    main()
