# 03 — Extract

Measurement over the coherence rasters, and retrieval of everything else.

| script | produces |
| --- | --- |
| `extract_dgamma.py` | coherence resolved by land cover class |
| `extract_rainmasked.py` | coherence over rain-wetted pixels |
| `extract_rainfull.py` | coherence **and** accumulation over one identical pixel set |
| `extract_wetdry_within.py` | wetted against unwetted ground inside one interferogram |
| `build_covariates.py` | ERA5-Land soil water, NDVI, soil texture |
| `pull_cells.py` | soil water at all four acquisitions |
| `fetch_bperp.py` | perpendicular baselines from the ASF service |

All take `--root <store>` and shard with `--shard N --of M` for a SLURM array.

---

## The within-burst contrast

`extract_wetdry_within.py` is the arm no confound touches. One rain mask,
derived from the event day, is applied to **both** interferograms:

```python
m = rain_mask_for(group, c_case, t_event)      # one mask, from the event
ok = land & finite(gc) & finite(gk) & (gc > 0) & (gk > 0)
wet, dry = ok & m, ok & ~m
stats(gc, gk, wet, "wet")     # both rasters, same selection
stats(gc, gk, dry, "dry")     # both rasters, same selection
```

Each population gets its own control baseline from the same control
interferogram, which is what removes the "that patch is just rougher"
explanation. The two populations also lie beneath one sky at one instant, so
any tropospheric delay is common to them and cancels.

!!! note "Why only 4,167 units qualify"

    A burst contributes only while a storm covers part of it. Units wholly wet
    or wholly dry are discarded, and each population must span at least one
    whole IMERG cell so that the comparison is of ground rather than of the
    burst outline clipping a cell corner.

---

## Masking both sides

`extract_rainfull.py` exists because masking the coherence while leaving the
accumulation as a burst mean measures the response over one patch of ground
and the dose over another. It retrieves the IMERG field itself rather than a
threshold of it, caches it as float16 beside each burst, and averages
coherence and accumulation over **one identical pixel set**.

??? note "Nearest neighbour, deliberately"

    The rain field is resampled to the coherence grid by nearest neighbour. An
    IMERG cell is a 0.1° block; interpolating it would invent a gradient
    across the wet edge that does not exist.
