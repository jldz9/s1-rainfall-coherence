# Rainfall and Sentinel-1 interferometric coherence — the data

Data supporting *A global dose–response for rainfall and Sentinel-1
interferometric coherence*, Li et al.

Every table, figure and number in the paper is produced from the ten files in
this directory by the code at
https://github.com/jldz9/s1-rainfall-coherence. The same files are archived at
doi:10.5281/zenodo.XXXXXXX.

Licence: CC-BY 4.0. Please cite the paper and the deposit.

---

## What is and is not here

The deposit holds **derived** products. The Sentinel-1 bursts they were computed
from are Copernicus open data and are not duplicated here; they may be retrieved
from the Alaska Satellite Facility (https://search.asf.alaska.edu) using
`00_unit_manifest.csv`, which gives the burst identifier and the four
acquisition dates of every unit. With that manifest and the processing
parameters in the paper's 'Interferometric processing' section, every
interferogram in the study can be rebuilt.

Coherence rasters are not deposited. There are 31,580 burst groups at roughly
75 MB each, some 300 GB in all, and all of it is regenerable from the manifest.

---

## Files, and what each one is for

| file | rows | what it holds | what it produces |
| --- | ---: | --- | --- |
| `00_unit_manifest.csv` | 21,010 | burst and the four acquisition dates per unit | nothing directly — the recipe for rebuilding the interferograms |
| `01_units_analysed.parquet` | 21,010 | the analysed sample, one row per unit | Table 3, Eq. (3), the 11.8% headline, Figs. 4b, 5, 6, 7 |
| `02_units_by_landcover.parquet` | 221,229 | coherence by cover class within each unit | Table 4, the paired cover contrasts |
| `03_validation_arm.parquet` | 1,907 | two **dry** interferograms of one burst | the zero point (+0.17%) and the floor (12.4%); Fig. 4a |
| `04_shared_control_arm.parquet` | 4,936 | rainfall and dry pairs against **one shared** control | Fig. 3a, +18.70% against +0.28% |
| `05_rain_masked.parquet` | 23,015 | coherence over wetted pixels only, beside the all-pixel value | the masking check: 0.13 points |
| `06_perpendicular_baselines.parquet` | 89,069 | perpendicular baseline per interferogram | Fig. 5, the baseline cannot bias the result |
| `07_within_burst_wet_dry.parquet` | 4,167 | wetted against unwetted ground **inside one image** | Fig. 3c |
| `08_soil_water.parquet` | 25,768 | ERA5-Land soil water at all four acquisitions | Fig. 8, the dielectric and mechanical pathways |
| `09_event_census.parquet` | 166,117 | every rainfall event screened, processed or not | the population weights — the difference between 19.2% and **11.8%** |

### How they fit together

`01` is the centre. Of the twelve analyses — six tables and six figures —
eleven read it, Figure 3 being the only one that does not. Everything else
either extends it or defends it.

**Three tables defend the design**, each removing a different objection.
`03` asks what the method returns when there is no rain to find. `04` removes
the matching by giving both arms the same denominator. `07` removes the second
interferogram altogether by comparing two populations of pixels within one
image.

**Three extend it along a covariate** — `02` by land cover, `06` by baseline,
`08` by soil water. `05` is a sensitivity check.

`09` is the one most easily skipped and the most consequential. The sampling is
stratified: every event of 15 mm or more was processed where its acquisitions
allowed, while lighter events were sampled, so the units over-represent heavy
rain roughly fivefold. Any quantity averaged **across** accumulations must be
weighted by this census first. Quantities already conditioned on the
accumulation must not be.

---

## Definitions

A **unit** is one rainfall event expressed as a comparison between two
interferograms of the same burst stack. The **case** pair spans the event, one
of its acquisitions having received the rain and the other, twelve days away,
less than 1 mm. The **control** pair is dry at both acquisitions and is drawn
from the same burst, the same twelve-day baseline and the same season.

The **reduction** is `(1 - g_case / g_ctrl) * 100`, in percentage points. It is
reported as a ratio because coherence magnitude is specific to the processing
chain; see the paper's Limitations.

---

## Column reference

### 01_units_analysed.parquet

| column | meaning |
| --- | --- |
| `unit_id` | burst and event date |
| `burst` | OPERA burst identifier; the unit of clustering for all intervals |
| `g_case`, `g_ctrl` | burst-mean coherence of the event and control pairs |
| `reduction_pct` | `(1 - g_case/g_ctrl) * 100`; the paper's headline quantity |
| `dgamma` | `g_case - g_ctrl`, in coherence units |
| `rain_mm` | IMERG 24 h accumulation before the event acquisition, burst mean |
| `r0`, `r1` | accumulation at the first and second acquisition of the case pair |
| `headroom` | `g_ctrl - 0.077`; units below 0.15 are excluded |
| `koppen`, `season`, `dom_cover` | climate class, hemisphere-corrected season, dominant cover |
| `frac_gt1`, `frac_gt5`, `frac_gt15` | fraction of the burst exceeding 1, 5 and 15 mm |
| `liquid_pct` | liquid fraction of the precipitation, **0–100**; units below 90 are excluded |
| `t_event`, `t_ref`, `t_sec`, `ctrl_ref`, `ctrl_sec` | the four acquisition dates |
| `lon`, `lat`, `doy`, `btemp_days`, `n_px` | centre, day of year, baseline, valid pixels |

### 02_units_by_landcover.parquet

One row per unit and cover class. `cover` is the WorldCover class, `n_px` its
pixel count within the burst, `g_case`/`g_ctrl` its coherences.

**Do not filter this table on headroom.** Headroom is a property of the unit,
already applied when `01_units_analysed` was built. Re-applying it per cover row
discards the dim classes and reverses the sign of the wetland contrast.

### 03_validation_arm.parquet

`g_a`, `g_b` are two dry interferograms of the same burst. `ctrl` is **one of
those two, nominated at random** as the denominator — not a third
interferogram — and `red` is the signed contrast of the other against it.
`ddoy` and `dyear` give their separation in day of year and in years.

Plotting `g_a` against `g_b` is not the comparison `red` describes. Put `ctrl`
on the control axis to reproduce the paper's figures and its 50.7%.

### 04_shared_control_arm.parquet

`g_shared` is one dry interferogram serving as the denominator for both the
rainfall pair (`red_rain`) and a second dry pair (`red_dry`). `paired_diff` is
their within-unit difference, which is the quantity of interest; the gap
between the two medians is not, a median not being additive.

### 05_rain_masked.parquet

`_all` columns are averaged over every valid land pixel, `_wet` over those in
cells that received 1 mm or more. `rain_wet` is the accumulation over the same
wetted pixels, so response and dose share one support.

### 06_perpendicular_baselines.parquet

One row per interferogram, keyed by `unit_id` and `role` (`case` or `control`),
with `bperp` in metres beside the scene identifiers and dates it was computed
from. Pivot on `role` for the case-minus-control difference per unit. Some rows
carry no baseline; drop them rather than filling. Dates are stored as strings.

### 07_within_burst_wet_dry.parquet

Wetted and unwetted populations of one interferogram, each with its own control
baseline. `d10`–`d90` are percentiles of the per-pixel difference.

### 08_soil_water.parquet

ERA5-Land `volumetric_soil_water_layer_1` (0–7 cm), daily mean on each
acquisition date, in `case_ref`, `case_sec`, `ctrl_ref`, `ctrl_sec`. These hold
soil water **values**, not dates, despite the names; join to
`01_units_analysed` on `unit_id` for the dates.

Orient the change to the event as `theta(wet) - theta(dry)`; the wet
acquisition is the one whose date equals `t_event`. Taking `sec - ref` without
orienting cancels the sign for roughly half the sample and collapses the
correlation from +0.288 to −0.086.

### 09_event_census.parquet

Every rainfall event the screening identified, processed or not, with its burst
and accumulation `mm`. Bin it on the same edges as the units, divide the census
count in each band by the units processed in that band, and those are the
weights. They sum to 166,117 and give an effective sample of 10,873.

---

## Products used

| product | use |
| --- | --- |
| Sentinel-1 SLC, OPERA burst grid | interferograms |
| GPM IMERG V07 half-hourly | rainfall; all accumulations reported |
| CHIRPS daily, ERA5-Land daily | candidate screening only, never a magnitude |
| ERA5-Land hourly | soil water, snow, temperature |
| ESA WorldCover 2021 | land cover |
| Köppen-Geiger 1980–2016 | climate class |
| MODIS MOD13A2 | NDVI |
