# s1-rainfall-coherence

Code for *A global dose–response for rainfall and Sentinel-1 interferometric
coherence*, Li et al.

Coherence loss is measured over 21,010 rainfall events, each matched against a
dry interferogram of the same burst, baseline and season. A rainfall event
costs **11.8%** of coherence once the sample is weighted to the census it was
drawn from. The relationship is logarithmic, `L = 0.82 + 5.27 log₂R`, and it
has no threshold: one to three millimetres already costs 6.0%.

---

## Install

```bash
git clone https://github.com/jldz9/s1-rainfall-coherence
cd s1-rainfall-coherence
pip install -e .            # add [extract] for the retrieval scripts
```

## Data

**The tables ship with this repository**, in `data/` (27 MB, ten files). A
clone needs no download and no configuration. The same tables are archived at
**doi:10.5281/zenodo.XXXXXXX** under CC-BY 4.0.

```bash
rcoh population      # works straight after the clone
```

To run against a copy held elsewhere, point `RCOH_DATA` at it:

```bash
export RCOH_DATA=/path/to/deposit
```

## Use

```bash
rcoh dose        # Table 3, reduction by accumulation band
rcoh landcover   # Table 4, extra loss by cover against grassland
rcoh climate     # Table 5, by Köppen class and by season
rcoh population  # the headline, weighted to the event census
rcoh fit         # Equation 3, with the comparison of functional forms
rcoh floor       # the dry-against-dry measurement floor
```

```python
from rainfall_coherence import io, tables, stats

u = io.units()                       # the 21,010 analysed units
tables.dose_relationship()           # {'slope': 5.27, 'r2_log': 0.1032, ...}

# every interval resamples bursts, never units
stats.cluster_bootstrap(u.reduction_pct, u.burst)
```

---

## Layout

```
data/                       the ten deposited tables, 27 MB; nothing to download

src/rainfall_coherence/     analysis of the deposited tables; no raw data needed
    config.py     constants fixed for the study: the floor, the headroom
                  condition, the wet-day threshold, the bands of Table 3
    io.py         loaders for the deposited tables, with derived columns
    stats.py      cluster bootstrap, paired differences, the dose fit
    tables.py     Tables 2, 3 and 4, and Equation 3
    figures/      one module per figure

scripts/                    the pipeline, in the order it was run
    01_acquire/   the burst grid, the Sentinel-1 timeline, the rainfall
                  census from IMERG, and which bursts were affected
    02_process/   unit selection, staging, and the SLURM jobs that form
                  the interferograms
    03_extract/   measurement over the resulting coherence rasters, and
                  retrieval of soil water, baselines and the rain masks
    04_tables/    assembly of the deposited tables
```

Phases 01 to 03 need the burst store and the `[extract]` extras; phase 04 and
the package itself run on the deposit alone.

---

## Four things that will bite you

**Weight before averaging across accumulations.** Every event above 15 mm was
processed while lighter events were sampled, so the units over-represent heavy
rain roughly fivefold. `io.units().reduction_pct.mean()` gives 19.21%, which
describes the sample; `tables.population_estimate()` gives **11.76%**, which
describes the rainfall. Do not weight a quantity already conditioned on dose —
Table 3's bands, Eq. (3), the cover and climate contrasts, both validation
arms.


**Never filter the land-cover table on headroom.** The condition
`g_ctrl - 0.077 >= 0.15` belongs to the *unit* and is applied when the analysed
sample is built. Applying it again to the rows of `02_units_by_landcover`
removes the dim covers, moves snow and ice from −12.6 to −4.8, and reverses the
sign of the wetland contrast from −1.6 to +1.6. `io.landcover()` does not apply
it; `tables.landcover_contrasts()` relies on that.

**Orient the soil water change to the event.** Δθ must be θ at the wet
acquisition minus θ at the dry one. Taking `sec - ref` inverts the sign for
roughly half the sample and the correlation collapses from +0.288 to −0.086.
`io.soil_water()` orients it.

**Resample bursts, not units.** Units are clustered within bursts, a mean of
nine to each and up to forty-six. Resampling units understates the spread by a
factor of 1.7. `stats.cluster_bootstrap` takes the burst as its group.

---

## Reproducing the interferograms

They are not deposited — 31,580 burst groups at roughly 75 MB each. They are
regenerable from `00_unit_manifest.csv`, which gives the burst identifier and
the four acquisition dates of every unit:

1. retrieve the bursts from the Alaska Satellite Facility
2. form twelve-day pairs with ISCE3/dolphin at 40 m posting
3. estimate coherence with a phase-only estimator, 11 × 11 Gaussian window of
   standard deviation one third the window width, and no adaptive phase filter
4. run the scripts in `scripts/`

Coherence magnitude is specific to this chain. A different window or estimator
will not reproduce the absolute values, and a chain applying a Goldstein filter
will report smaller reductions.

---

## Not yet included

Two scripts from the published analysis are being reconstructed: the generator
for Figure 2, and the builder of Table 4 as originally run.
`tables.landcover_contrasts()` reproduces Table 4 exactly from the deposited
data, and `fig_r12_dose` was rewritten against the bands of Table 3.

The interferograms themselves are formed by InSARHub with an ISCE3/dolphin
backend, which is a separate package; `02_process/` drives it rather than
containing it.

## Licence

MIT. The data deposit is CC-BY 4.0. Please cite the paper and the deposit.
