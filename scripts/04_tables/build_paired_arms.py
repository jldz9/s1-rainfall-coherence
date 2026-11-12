#!/usr/bin/env python
"""The rainfall and validation contrasts taken against one shared control.

A burst visited by a rainfall event holds a control interferogram, and where the
validation arm paired a second dry interferogram against that same control, the
two contrasts can be formed with a common denominator.  Everything belonging to
the burst then cancels in the subtraction rather than merely matching: not the
terrain, not the look geometry, not the season, not the control itself.

Two points of construction matter.  The control appears twice, once in each
arm's extraction, and the two measurements differ slightly because each arm
masks its own pixels; rows where they disagree by 1% or more are not the same
interferogram and are dropped.  For the rest the rainfall arm's own value is
used as the single denominator of both ratios, which keeps the rainfall contrast
within one extraction and places the cross-run difference in the dry contrast,
where it works against the result rather than for it.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np, pandas as pd

TOL = 1.0          # per cent disagreement above which the control is not shared


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--out", default="paired_rain_vs_null.parquet")
    a = ap.parse_args()
    R = a.root

    u = pd.read_parquet(R / "plan" / "units_final.parquet")
    v = pd.read_parquet(R / "plan" / "validation_arm_final.parquet")

    # which of the validation unit's two members is the rainfall unit's control
    # is not recorded — the case and control roles were assigned at random when
    # the arm was built — so both are tested and the closer one taken
    j = u.merge(v[["burst", "g_a", "g_b"]].drop_duplicates("burst"),
                on="burst", how="inner")
    n_join = len(j)

    da = (j.g_ctrl - j.g_a).abs() / j.g_ctrl * 100
    db = (j.g_ctrl - j.g_b).abs() / j.g_ctrl * 100
    a_is_ctrl = da <= db
    j["g_shared"] = np.where(a_is_ctrl, j.g_a, j.g_b)
    j["g_otherdry"] = np.where(a_is_ctrl, j.g_b, j.g_a)
    j = j[np.minimum(da, db) < TOL].copy()

    # one denominator for both ratios
    j["red_rain"] = (1 - j.g_case / j.g_ctrl) * 100
    j["red_dry"] = (1 - j.g_otherdry / j.g_ctrl) * 100
    j["paired_diff"] = j.red_rain - j.red_dry

    cols = ["unit_id", "burst", "rain_mm", "g_case", "g_ctrl", "g_otherdry",
            "g_shared", "red_rain", "red_dry", "paired_diff"]
    j[cols].to_parquet(R / "plan" / a.out, index=False)

    d = j.paired_diff
    rng = np.random.default_rng(0)
    grp = {b: x.values for b, x in d.groupby(j.burst)}
    ks = list(grp)
    bs = [np.median(np.concatenate([grp[k] for k in rng.choice(ks, len(ks), True)]))
          for _ in range(4000)]
    lo, hi = np.percentile(bs, [2.5, 97.5])

    print(f"joined on burst            {n_join:,}")
    print(f"control shared to <{TOL}%    {len(j):,} over {j.burst.nunique():,} bursts")
    print(f"  rainfall contrast        {j.red_rain.median():.2f}%")
    print(f"  dry contrast             {j.red_dry.median():.2f}%")
    print(f"  paired difference        {d.median():.2f} pp  [{lo:.2f}, {hi:.2f}]")
    print(f"  rainfall lower           {(j.red_rain > j.red_dry).mean()*100:.1f}%")
    print("  by accumulation:")
    for lo_, hi_, lab in ((0, 5, "<5 mm"), (5, 15, "5-15 mm"),
                          (15, 30, "15-30 mm"), (30, 1e9, ">30 mm")):
        m = j[(j.rain_mm >= lo_) & (j.rain_mm < hi_)]
        if len(m):
            print(f"    {lab:>9}  n={len(m):5,}  {m.paired_diff.median():5.1f} pp")
    print(f"wrote {R / 'plan' / a.out}")


if __name__ == "__main__":
    main()
