# Data dictionary

Nine files, 25 MB, deposited at **doi:10.5281/zenodo.XXXXXXX** under CC-BY 4.0.

---

## Definitions

A **unit** is one rainfall event expressed as a comparison between two
interferograms of the same burst stack.

| term | meaning |
| --- | --- |
| case | the twelve-day pair spanning the event |
| control | a dry pair, same burst, same baseline, same season |
| reduction | `(1 - g_case / g_ctrl) * 100`, in percentage points |
| headroom | `g_ctrl - 0.077`; a unit needs 0.15 to be analysed |

---

## `00_unit_manifest.csv` — 21,010 rows

The file that makes the study rebuildable from the public archive.

| column | meaning |
| --- | --- |
| `unit_id` | burst and event date |
| `burst` | OPERA burst identifier |
| `t_event` | the day the rain fell |
| `t_ref`, `t_sec` | the two acquisitions of the case pair |
| `ctrl_ref`, `ctrl_sec` | the two acquisitions of the control pair |
| `lon`, `lat` | burst centre |

---

## `01_units_analysed.parquet` — 21,010 rows

| column | meaning |
| --- | --- |
| `g_case`, `g_ctrl` | burst-mean coherence of the event and control pairs |
| `reduction_pct` | the headline quantity |
| `dgamma` | `g_case - g_ctrl`, in coherence units |
| `rain_mm` | IMERG 24 h accumulation, burst mean |
| `r0`, `r1` | accumulation at the first and second case acquisition |
| `headroom` | `g_ctrl - 0.077` |
| `koppen`, `season`, `dom_cover` | climate, hemisphere-corrected season, cover |
| `frac_gt1`, `frac_gt5`, `frac_gt15` | fraction of the burst above 1, 5, 15 mm |
| `liquid_pct` | liquid fraction; below 90% is excluded |
| `btemp_days`, `n_px`, `doy` | baseline, valid pixels, day of year |

## `02_units_by_landcover.parquet` — 221,229 rows

One row per unit and cover class. `cover` is the WorldCover class, `n_px` its
pixel count, `g_case`/`g_ctrl` its coherences.

!!! danger "Do not filter on headroom"

    See [Pitfalls](../reproduce/pitfalls.md). It reverses the wetland sign.

## `03_validation_arm.parquet` — 1,907 rows

`g_a`, `g_b` are two dry interferograms of one burst; `red` is the signed
contrast with the numerator assigned at random; `ctrl` is the nominated
denominator; `ddoy` and `dyear` their separation in day of year and in years.

## `04_shared_control_arm.parquet` — 4,936 rows

`g_shared` serves as the denominator for both the rainfall pair (`red_rain`)
and a second dry pair (`red_dry`); `paired_diff` is their within-unit
difference.

## `05_rain_masked.parquet` — 23,015 rows

`_all` columns average over every valid land pixel, `_wet` over those in cells
receiving 1 mm or more. `rain_wet` is the accumulation over the same wetted
pixels, so response and dose share one support.

## `06_perpendicular_baselines.parquet` — 89,069 rows

One row per interferogram, `role` distinguishing case from control.

## `07_within_burst_wet_dry.parquet` — 4,167 rows

Wetted and unwetted populations of one interferogram, each with its own
control baseline. `d10`–`d90` are percentiles of the per-pixel difference.

## `08_soil_water.parquet` — 25,768 rows

ERA5-Land `volumetric_soil_water_layer_1` (0–7 cm), daily mean on each
acquisition date.

!!! danger "Orient the change to the event"

    Δθ is θ(wet) − θ(dry). The wet acquisition is the one equal to `t_event`.
    Taking `sec - ref` cancels the sign for half the sample.

---

## `09_event_census.parquet` — 166,117 rows

Every rainfall event the units were drawn from, with its burst, date and IMERG
accumulation. Needed for the population weights, since the sampling was
stratified upon the accumulation.

---

## Products

| product | use |
| --- | --- |
| Sentinel-1 SLC, OPERA burst grid | interferograms |
| GPM IMERG V07 half-hourly | rainfall; every accumulation reported |
| CHIRPS daily, ERA5-Land daily | candidate screening only, never a magnitude |
| ERA5-Land hourly | soil water, snow, temperature |
| ESA WorldCover 2021 | land cover |
| Köppen-Geiger 1980–2016 | climate class |
| MODIS MOD13A2 | NDVI |
