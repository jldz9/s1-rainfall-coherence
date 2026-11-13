# Pitfalls

Three mistakes made during the analysis, each of which changed a published
number before it was caught. The package guards against all three, but anyone
working with the tables directly can still walk into them.

---

## 1. The headroom condition belongs to the unit

!!! danger "Never filter the land-cover table on headroom"

    ```python
    # WRONG
    d = io.landcover()
    d = d[d.g_ctrl - 0.077 >= 0.15]     # do not do this
    ```

A control must stand 0.15 above the noise floor for a loss to be resolvable.
That condition is a property of the **unit**, applied once to the burst mean
when the analysed sample was built.

Applying it again to the rows of `02_units_by_landcover` asks a different
question — whether each *cover class* is bright — and conditions on the
denominator of the very ratio being measured. It removes the dim covers
preferentially:

| cover | median $\gamma_{ctrl}$ | correct | with the stray filter | rows dropped |
| --- | ---: | ---: | ---: | ---: |
| built-up | 0.619 | −5.4 | −5.5 | 1% |
| tree | 0.388 | +3.4 | +4.0 | 11% |
| wetland | 0.357 | **−1.6** | **+1.6** | 21% |
| snow and ice | 0.371 | **−12.6** | **−6.2** | 30% |
| mangrove | 0.276 | −4.1 | +0.1 | 37% |

The wetland contrast reverses sign. Snow and ice halves. The two covers quoted
in the abstract barely move, which is why the error survived three attempts to
find it.

`io.landcover()` does not apply the filter, and `tables.landcover_contrasts()`
relies on that.

---

## 2. The soil water change must be oriented to the event

!!! danger "Δθ is wet minus dry, not second minus first"

    ```python
    # WRONG
    dtheta = theta_sec - theta_ref
    ```

The campaign allows the rain at either acquisition of a pair. Taking
`sec - ref` inverts the sign for roughly half the sample, the two halves
cancel, and the correlation with coherence loss collapses:

| | ρ with the reduction |
| --- | ---: |
| oriented, θ(wet) − θ(dry) | **+0.288** |
| unoriented, θ(sec) − θ(ref) | −0.086 |

The wet acquisition is the one equal to `t_event`. `io.soil_water()` orients
it and exposes `dtheta`.

---

## 3. Resample bursts, never units

!!! danger "Units within a burst are not independent"

    ```python
    # WRONG
    bootstrap over rows
    ```

The 21,010 units occupy only 2,329 bursts, an average of nine to each and as
many as forty-six, and units of one burst share its terrain and its scattering
behaviour. Resampling units rather than bursts understates the spread by a
factor of **1.7**.

```python
stats.cluster_bootstrap(u.reduction_pct, u.burst)   # groups = burst
```

Every interval in the paper is computed this way.

!!! note "Pairing does not remove the clustering"

    Table 4 pairs each cover against grassland inside the same burst, which
    removes the burst from the contrast but not from the sampling: its units
    still share one terrain. Taken as $sd/\sqrt{n}$ over units its standard
    errors ran 1.2 to 1.8 times too narrow, worst for tree cover. They
    resample bursts like every other interval in the paper.

---

## 4. Weight before averaging across accumulations

!!! danger "The sample over-represents heavy rain by design"

    ```python
    # WRONG, if you want a population figure
    io.units().reduction_pct.mean()      # 19.21%
    ```

Every event above 15 mm was processed where its acquisitions allowed, while
lighter events were sampled for precision. So heavy rain is 58.4% of the units
against 12.0% of the events recorded, and the unweighted mean describes the
sample rather than the rainfall.

```python
tables.population_estimate()["weighted_mean_pct"]    # 11.76%
```

| | unweighted | weighted |
| --- | ---: | ---: |
| mean loss | 19.21% | **11.76%** |
| displacement noise | 1.43 mm | **0.88 mm** |
| falling below control | 80.7% | 72.4% |
| median accumulation | 16.5 mm | **3.9 mm** |

!!! note "Do not weight a conditional quantity"

    Table 3's bands, Eq. (3), the cover and climate contrasts and both
    validation arms are already conditioned on dose or paired within units.
    Weighting them adds noise without removing bias, and would change Eq. (3)
    from `0.82 + 5.27 log₂R` to `1.46 + 4.98 log₂R` for no gain.

The weighted headline is also the one that reconciles with the curve: the
weighted median accumulation is 3.9 mm, for which Eq. (3) gives 11.2%.

---

## Two more worth knowing

**A median is not additive.** The median of the within-unit differences is not
the difference of the two medians, and only the former respects the pairing.
In the within-burst arm the wetted population falls 9.25% and the unwetted
6.05%, a gap of 3.20 — but the median paired difference is **1.64**. Use
`stats.paired_difference`.

**Coherence magnitude is not portable.** These values come from a phase-only
estimator with an 11 × 11 Gaussian window on a 40 m pixel, 440 m square on the
ground, and a floor of 0.077. A chain with a different window, a different
multilook, a maximum-likelihood estimator or a Goldstein filter will not
reproduce them, and will report different reductions. Compare the ordering and
the form of the response, not the magnitudes.

!!! note "One parameter is unresolved"

    The manuscript states the estimator's effective number of looks as 169 in
    one place and 133 in another. Both reach the 0.077 floor, by different
    conventions: $1/\sqrt{169}=0.0769$ is the rms of the estimate at zero
    coherence, $\sqrt{\pi}/(2\sqrt{133})=0.0768$ its mean. Neither matches
    the 102 an 11 × 11 Gaussian with $\sigma=w/3$ would give on independent
    pixels, which is probably because the SHP test varies the sample per
    pixel. The floor itself is empirical and unaffected.
