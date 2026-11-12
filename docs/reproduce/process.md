# 02 — Process

Choosing the units and forming the interferograms.

| script | what it does |
| --- | --- |
| `build_campaign.py` | unit selection: case and control pairing, the dry rules |
| `stage_groups.py` | materialises each processing group for the cluster |
| `submit_managed.py` | one SLURM job per group |
| `wx_submit_daemon.py` | sliding-window submitter, respecting the queue limits |

---

## How a unit is built

For each qualifying event, `build_campaign.py` forms the **case** pair, a
twelve-day interferogram straddling the rain, then searches the same burst
stack for a **control**: a twelve-day pair, dry at both acquisitions, falling
in the same season, and between one and five years from the event. Where
several qualify, the one closest to the event's point in the year is taken.

```python
WET_MM, DRY_MM = 15.0, 1.0   # see the warning below
```

!!! warning "`WET_MM` in this script is not the wet-day threshold"

    The name is misleading and is kept only because it is what was run. In
    `build_campaign.py`, written for the heavy arm:

    | constant | value | what it actually is |
    | --- | ---: | --- |
    | `WET_MM` | 15.0 | the accumulation above which **every** available event was processed; below it events were sampled |
    | `DRY_MM` | 1.0 | the criterion deciding whether an acquisition is dry enough for a control |

    The paper's wet-day threshold is **1 mm**, the ETCCDI definition, and that
    is also the event floor. 15 mm is a statement about **coverage**, not
    about rainfall. The package names them unambiguously as
    `config.WET_DAY_MM`, `config.DRY_ACQ_MM` and
    `config.FULL_PROCESSING_MM`.

Four conditions then reduce 30,593 formed units to the 21,010 analysed:

| condition | removes |
| --- | ---: |
| headroom: the control too near the noise floor | 4,825 |
| no liquid fraction returned by IMERG | 112 |
| snow or mixed phase, below nine tenths liquid | 2,641 |
| a non-event acquisition wetter than 1 mm | 2,005 |

---

## The processing chain

!!! warning "Fixed for the whole study"

    | parameter | value |
    | --- | --- |
    | backend | ISCE3 / dolphin, driven through InSARHub |
    | DEM | Copernicus GLO-30, buffered by 0.5° |
    | geocoding | OPERA burst grid, 10 m range × 20 m azimuth |
    | multilook | 4 range × 2 azimuth, giving a **40 m square pixel** |
    | ionosphere | JPL total electron content |
    | SHP selection | generalised likelihood ratio test, significance 0.001 |
    | estimator | phase-only |
    | window | 11 × 11 Gaussian, σ = one third the window width |
    | window on the ground | 11 × 40 m = **440 m square** |
    | adaptive phase filter | **none** |
    | noise floor | 0.077 |
    | temporal baseline | twelve days, every pair processed alone |
    | polarisation | VV |

    The absence of a Goldstein filter is deliberate. It raises coherence most
    where coherence is lowest, which is where the rainfall cases sit, so
    filtering would compress the contrast under measurement by an amount that
    differs between chains.

---

## Cluster notes

??? note "Queue limits"

    The submitters assume a per-user cap on queued array tasks. Raise the
    throttle on a running array rather than cancelling and resubmitting:

    ```bash
    scontrol update jobid=<id> MinCPUsNode=1 MinMemoryNode=4096
    ```

??? note "Checkpointing"

    The extraction scripts write their parquet every 100 units so a walltime
    kill cannot discard hours of Earth Engine retrieval. Resume with
    `--skip-done`.
