# API

```python
from rainfall_coherence import config, io, stats, tables
```

---

## `config`

Constants fixed for the whole study, defined once.

| name | value | |
| --- | ---: | --- |
| `NOISE_FLOOR` | 0.077 | coherence where no correlation survives |
| `MIN_HEADROOM` | 0.15 | how far a control must stand above it |
| `WET_DAY_MM` | 1.0 | the ETCCDI wet-day threshold |
| `EVENT_FLOOR_MM` | 1.0 | alias: the smallest accumulation counted as an event |
| `DRY_ACQ_MM` | 1.0 | alias: dry enough to enter a control pair |
| `FULL_PROCESSING_MM` | 15.0 | above this every event was processed, below it sampled |
| `DOSE_EDGES` | | the accumulation bands of Table 3 |
| `LOOKS_DISPLACEMENT` | 4 | multilooking for the phase conversion |
| `WAVELENGTH_MM` | 55.5 | Sentinel-1 C-band |

---

## `io`

Every loader takes `root=`, else reads `RCOH_DATA`, else falls back to the
`data/` directory bundled with the repository. All attach derived columns.

| function | returns |
| --- | --- |
| `units()` | the analysed sample, with `reduction_pct` and `headroom` |
| `landcover(analysed_only=True)` | cover-resolved; **no headroom filter** |
| `validation()` | dry-against-dry pairs |
| `shared_control()` | the shared-denominator arm |
| `masked(analysed_only=True)` | with `loss_all` and `loss_wet` |
| `within_burst()` | with `loss_wet`, `loss_dry` and `diff` |
| `soil_water()` | with `dtheta` oriented to the event |
| `baselines()` | perpendicular baselines |
| `census()` | the 166,117 events the units were drawn from |
| `population_weights()` | census events per band over units processed |

---

## `stats`

| function | |
| --- | --- |
| `reduction(g_case, g_ctrl)` | the percentage loss |
| `cluster_bootstrap(values, groups, stat, n, seed)` | **groups must be the burst** |
| `paired_difference(a, b, groups)` | within-unit, then summarised |
| `phase_sigma(gamma, looks=4)` | Cramér–Rao phase standard deviation |
| `phase_to_displacement_mm(sigma)` | line-of-sight noise in mm |
| `log_fit(x, y)` | `(intercept, slope, r2, zero_mm)` |

```python
stats.cluster_bootstrap(u.reduction_pct, u.burst)
# (point, low, high) resampling bursts, never units
```

---

## `tables`

| function | produces |
| --- | --- |
| `dose()` | Table 3 |
| `landcover_contrasts(reference="grass")` | Table 4 |
| `climate_and_season()` | Table 5, as two frames |
| `dose_relationship()` | Equation 3 and the form comparison |
| `measurement_floor()` | the floor, by cover and by climate |
| `population_estimate()` | the headline weighted to the event census |

---

## CLI

```
rcoh {dose,landcover,climate,fit,floor,population} [--data DIR]
```
