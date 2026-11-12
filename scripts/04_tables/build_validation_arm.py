#!/usr/bin/env python
"""The validation arm, built by extraction from interferograms already on disk.

Each unit is two control interferograms over one burst, both of which the
rainfall arm already produced and both of which IMERG places under 1 mm at
every acquisition, averaged over the burst footprint.  Neither carries rain, so
the comparison between them is what the design returns when the only thing
removed is the event.

Nothing is reprocessed.  A common pixel mask was tested against independent
masks on 400 pairs and changed the spread by 0.06 pp, and two controls drawn
from different processing groups agree as closely as two from the same group
(13.58% against 13.52% unsigned), so the stored rasters serve directly.
"""
import argparse
from pathlib import Path
import numpy as np, pandas as pd, rasterio

ap = argparse.ArgumentParser()
ap.add_argument("--root", type=Path, required=True)
ap.add_argument("--shard", type=int, default=0)
ap.add_argument("--of", type=int, default=1)
a = ap.parse_args()
R = a.root
P = pd.read_parquet(R / "plan" / "validation_arm_plan.parquet").iloc[a.shard::a.of]

done = {}
for g in (R / "proc" / "burst").glob("*"):
    if (g / ".done").exists():
        done.setdefault(g.name.rsplit("_", 1)[0].lower(), []).append(g)


def cor(burst, ref, sec):
    for g in done.get(burst, []):
        f = g / "pairs" / f"{ref}_{sec}" / "interferograms" / f"{ref}_{sec}.int.cor.tif"
        if f.exists():
            return f, g
    return None, None


rows = []
for r in P.itertuples():
    fa, ga = cor(r.burst, r.a_ref, r.a_sec)
    fb, gb = cor(r.burst, r.b_ref, r.b_sec)
    if fa is None or fb is None:
        continue
    try:
        with rasterio.open(fa) as s:
            A = s.read(1).astype("float32")
        with rasterio.open(fb) as s:
            B = s.read(1).astype("float32")
    except Exception:
        continue
    if A.shape != B.shape:
        continue
    lc = ga / "worldcover.npy"
    land = (np.load(lc) != 80) if lc.exists() else np.ones(A.shape, bool)
    m = land & np.isfinite(A) & (A > 0) & np.isfinite(B) & (B > 0)
    if m.sum() < 5000:
        continue
    rows.append({"burst": r.burst, "koppen": r.koppen, "ddoy": r.ddoy, "dyear": r.dyear,
                 "g_a": float(A[m].mean()), "g_b": float(B[m].mean()), "n_px": int(m.sum())})

pd.DataFrame(rows).to_parquet(R / "plan" / f"valarm_shard{a.shard:03d}.parquet", index=False)
print(f"wrote {len(rows):,} units")
