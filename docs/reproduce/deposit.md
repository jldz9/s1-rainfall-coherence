# With the provided data

Every table, every fitted relationship and every figure in the paper
regenerates from the 25 MB deposit in under a minute. No imagery, no Earth
Engine account, no cluster.

```bash
```

---

## From the command line

=== "Equation 3"

    ```bash
    rcoh fit
    ```

    ```
      intercept            0.82
      slope                5.27
      r2_log               0.1032
      r2_linear            0.0764
      r2_sqrt              0.0991
      zero_mm              0.898
      n                    21010
      slope_above_15mm     6.46
      n_above_15mm         12272
    ```

    The three $R^2$ values are the comparison of functional forms reported in
    the Results: a logarithm explains 10.32% of the variance against 7.64% for
    a straight line and 9.91% for a square root.

=== "Table 3"

    ```bash
    rcoh dose
    ```

    Reduction by accumulation band, medians with cluster bootstrap intervals
    over bursts.

=== "Table 4"

    ```bash
    rcoh landcover
    ```

    ```
       cover  units  extra_loss   se
        tree  20291         3.4 0.11
       shrub   8765         2.3 0.17
     wetland   8962        -1.6 0.24
       crops  17142        -3.0 0.12
    mangrove    887        -4.0 0.91
        bare  18606        -5.1 0.14
       built  17537        -5.4 0.11
        moss   2526        -6.9 0.54
     snowice   1753       -12.5 1.00
    ```

=== "Table 5"

    ```bash
    rcoh climate
    ```

    Reduction by Köppen class and by hemisphere-corrected season.

=== "The floor"

    ```bash
    rcoh floor
    ```

    The dry-against-dry width, 12.4% at the median, broken out by cover and by
    climate. Cover spans 4.0% on bare ground to 14.6% under trees; climate
    spans only 10.1 to 13.8, and does not follow humidity.

---

## From Python

```python
from rainfall_coherence import io, tables, stats

tables.dose()                   # Table 3
tables.landcover_contrasts()    # Table 4
tables.climate_and_season()     # Table 5
tables.dose_relationship()      # Equation 3
tables.measurement_floor()      # the floor, by cover and climate
```

The headline, weighted to the event census:

```python
tables.population_estimate()
# {'weighted_mean_pct': 11.76, 'ci': (11.17, 12.38),
#  'unweighted_mean_pct': 19.21, 'effective_n': 10873, ...}
```

The sample over-represents heavy rain roughly fivefold by design, so anything
averaged across accumulations is weighted first. The unweighted mean of the
processed units is available for comparison:

```python
u = io.units()
stats.cluster_bootstrap(u.reduction_pct, u.burst, stat=lambda v: v.mean())
# (19.21, 18.63, 19.82) -- the sample, not the rainfall
```

The validation arm, which fixes the zero point:

```python
v = io.validation()
stats.cluster_bootstrap(v.red, v.burst)     # (0.17, -0.69, 1.07)
```

The within-burst contrast, the one measurement immune to every confound
because both populations lie in the same image under the same sky:

```python
w = io.within_burst()
stats.cluster_bootstrap(w["diff"], w.burst) # (1.64, 1.33, 1.95)
```

---

## The figures

Each lands under the number it carries in the paper, as `figN_*.png` and
`figN_*.pdf`.

| figure | command | what it shows |
| --- | --- | --- |
| **3** | `python -m rainfall_coherence.figures.fig_paired_arms` | the shared-control and within-burst contrasts |
| **4** | `python -m rainfall_coherence.figures.fig_one_sidedness` | the floor falls either way, the rainfall loss does not |
| **5** | `python -m rainfall_coherence.figures.fig_baseline` | the perpendicular baseline affects but cannot bias |
| **6** | `python -m rainfall_coherence.figures.fig_dose` | the dose–response |
| **7** | `python -m rainfall_coherence.figures.fig_r12_dose` | event and control coherence separately |
| **8** | `python -m rainfall_coherence.figures.fig_rain_vs_moisture` | rainfall against soil water |

All of them at once:

```bash
for f in fig_paired_arms fig_one_sidedness fig_baseline fig_dose \
         fig_r12_dose fig_rain_vs_moisture; do
    python -m rainfall_coherence.figures.$f
done
```

??? note "Options"

    `fig_dose` takes `--linear` to draw the same data on a linear axis, which
    is how the concavity was checked. `fig_baseline` and `fig_paired_arms`
    take further arguments documented in their module docstrings.

!!! note "Figures 1 and 2"

    Figure 1 is a schematic, drawn by hand. The generator for Figure 2 is
    among the scripts being reconstructed; see the repository README.

---

## Verifying everything at once

```bash
pytest
```

Eleven tests and twenty-nine assertions, each tied to a published number: the sample size, the
headline, Equation 3 and its form comparison, five rows of Table 4, the
floor, the within-burst difference, the masking sensitivity, and that the
soil-water change is oriented. A pass means this copy of the code reproduces
the published analysis exactly. The suite runs against the bundled `data/`;
`RCOH_DATA` need not be set.

!!! tip "One test is a tripwire"

    `test_table3_landcover` asserts that the wetland contrast is negative. It
    is negative in the published table and becomes positive if anyone
    reinstates the headroom filter on the land-cover rows. See
    [Pitfalls](pitfalls.md).
