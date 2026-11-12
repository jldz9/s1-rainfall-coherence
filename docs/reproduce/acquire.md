# 01 — Acquire

Establishing what exists before anything is processed: the burst grid, the
Sentinel-1 acquisition timeline, the rainfall census, and the join between
them.

| script | what it does |
| --- | --- |
| `build_burst_grid.py` | the static OPERA burst grid, 1,127,661 bursts, ~392,000 on land |
| `build_timeline.py` | the Sentinel-1 IW SLC acquisition timeline per burst |
| `build_events.py` | the rainfall census |
| `mark_affected.py` | joins the timeline to the census: which acquisitions carry rain |
| `verify_availability.py` | confirms every planned acquisition is actually retrievable |
| `profile_bursts.py` | land cover and Köppen class per burst |

---

## The two-stage rainfall census

!!! important "Screening is not measurement"

    `build_events.py` runs two stages, and conflating them is the commonest
    serious error in event-based remote sensing.

**Stage A** screens candidates from daily CHIRPS within 50° of the equator and
ERA5-Land beyond it, taking the larger of the two. Both badly under-detect
convective extremes:

| product | Derna, 10 September 2023 |
| --- | ---: |
| CHIRPS daily | 3.3 mm |
| ERA5 total precipitation | 68.8 mm |
| **IMERG daily accumulation** | **206.7 mm** |

An order of magnitude, in exactly the gauge-sparse regions the campaign cares
about. Stage A therefore sets its bar far below any threshold used in the
analysis and serves only to find candidates.

**Stage B** measures with GPM IMERG V07 half-hourly on its native 0.1° grid,
accumulated over the 24 h preceding each acquisition and averaged over the
burst footprint. **Every accumulation reported in the paper is a Stage B
value.**

---

## Event definition

An event is a single day within a two-degree cell, so a storm crossing several
bursts is counted once. An event acquisition must reach **1 mm**, the ETCCDI
wet-day threshold.

Between 3 October 2014 and 23 September 2026 this gives **166,117 events**
across 2,764 burst stacks, of which 19,918 carry 15 mm or more.

!!! note "One threshold, three jobs"

    1 mm is the floor an event must reach, the criterion deciding whether an
    acquisition is dry enough to enter a control pair, and the spatial mask
    separating wetted from dry ground within a burst. It is defined once, in
    `config.WET_MM`.
