# From scratch

Rebuilding the measurement itself: retrieve the imagery, form the
interferograms, measure the coherence, assemble the tables, then run the
analysis.

!!! warning "What this costs"

    The published study processed 56,930 interferograms across 31,580 burst
    groups. End to end this is weeks of cluster time and roughly 300 GB. If
    you want to check the statistics rather than the interferometry, use
    [the deposit](deposit.md) instead.

---

## What you need

| | |
| --- | --- |
| Earth Engine | an authenticated cloud project |
| cluster | SLURM, or adapt the submitters |
| disk | ~300 GB for the burst store |
| [InSARHub](https://github.com/jldz9/InSARHub) | forms the interferograms; installed separately |
| this package | `pip install -e ".[extract]"` |

```bash
earthengine authenticate
pip install -e ".[extract]"
export WX_ROOT=/path/to/the/store      # where everything is written
```

---

## The four phases

```
01_acquire    what exists, where it rained, which bursts were affected
     |
02_process    choose the units, stage them, form the interferograms
     |
03_extract    measure coherence over the rasters, pull the covariates
     |
04_tables     assemble the tables the analysis reads
     |
          then: With the provided data
```

### 01 — Acquire

Establish the burst grid, the acquisition timeline and the rainfall census,
then join them.

```bash
python scripts/01_acquire/build_burst_grid.py     --root $WX_ROOT
python scripts/01_acquire/build_timeline.py       --root $WX_ROOT
python scripts/01_acquire/build_events.py         --root $WX_ROOT --project <ee-project>
python scripts/01_acquire/mark_affected.py        --root $WX_ROOT
python scripts/01_acquire/verify_availability.py  --root $WX_ROOT
python scripts/01_acquire/profile_bursts.py       --root $WX_ROOT --project <ee-project>
```

Produces the event census: 166,117 events over 2,764 burst stacks, of which
19,918 reach 15 mm. [Details](acquire.md)

!!! important "Screening is not measurement"

    `build_events.py` screens candidates with CHIRPS and ERA5-Land, then
    measures every accumulation with IMERG. The screening products
    under-detect convective extremes by an order of magnitude. Never quote a
    Stage A value as a rainfall amount.

### 02 — Process

Select the units, stage them, and form the interferograms.

```bash
python scripts/02_process/build_campaign.py   --root $WX_ROOT
python scripts/02_process/stage_groups.py     --root $WX_ROOT
python scripts/02_process/wx_submit_daemon.py --root $WX_ROOT
```

The daemon submits a SLURM job per group and keeps a sliding window inside the
queue limits. This is the phase that takes weeks. [Details](process.md)

??? tip "Skipping phase 01"

    `00_unit_manifest.csv` in the deposit already names the burst and the four
    acquisition dates of every unit. Feed it to `stage_groups.py` and you skip
    the census and the selection entirely, reproducing the interferometry
    against exactly the published sample.

### 03 — Extract

Measure over the resulting rasters and retrieve everything else. All of these
shard with `--shard N --of M` for a SLURM array.

```bash
python scripts/03_extract/extract_dgamma.py        --root $WX_ROOT --shard 0 --of 20
python scripts/03_extract/extract_wetdry_within.py --root $WX_ROOT --shard 0 --of 20
python scripts/03_extract/extract_rainfull.py      --root $WX_ROOT --shard 0 --of 20
python scripts/03_extract/build_covariates.py      --root $WX_ROOT --shard 0 --of 20
python scripts/03_extract/pull_cells.py            --root $WX_ROOT
python scripts/03_extract/fetch_bperp.py           --root $WX_ROOT
```

[Details](extract.md)

??? note "Checkpointing"

    Each writes its parquet every 100 units, so a walltime kill cannot discard
    hours of Earth Engine retrieval. Resume with `--skip-done`.

### 04 — Tables

Assemble the tables that the analysis reads.

```bash
python scripts/04_tables/build_paired_arms.py    --root $WX_ROOT
python scripts/04_tables/build_validation_arm.py --root $WX_ROOT
python scripts/04_tables/merge_and_analyse.py    --root $WX_ROOT
python scripts/04_tables/soil_moisture_arms.py   --root $WX_ROOT
```

[Details](tables.md)

### Then the analysis

The outputs of phase 04 are the deposit. Point the package at them and
continue exactly as in [With the provided data](deposit.md):

```bash
export RCOH_DATA=$WX_ROOT/plan
rcoh fit
```

---

## Checking you landed in the same place

```bash
pytest
```

Eleven tests tied to published numbers. If your rebuild reproduces them,
the chain matched; if `test_equation_3` passes but the headline differs, the
difference is in unit selection rather than in the interferometry.

!!! warning "Coherence magnitude is not portable"

    A different window, multilook or estimator will not reproduce the absolute
    coherences, and a chain applying a Goldstein filter will report smaller
    reductions. The parameters are fixed in [02 Process](process.md) and every
    one of them matters. Compare the ordering and the form of the response
    before comparing magnitudes.
