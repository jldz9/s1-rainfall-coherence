# Getting the data

There is nothing to get. The nine tables and the manifest the analysis reads,
27 MB in all and licensed CC-BY 4.0, ship with the repository in `data/`.

```bash
git clone https://github.com/jldz9/s1-rainfall-coherence
cd s1-rainfall-coherence
pip install -e .
rcoh population          # no download, no configuration
```

The loaders find `data/` by walking up from the installed package, so this
works from any working directory in an editable install.

!!! info "Reading a copy held elsewhere"

    Set `RCOH_DATA` to override the bundled tables, or pass `root=` to any
    loader. The same tables are archived at doi:10.5281/zenodo.XXXXXXX, which
    is the citable form; the copy here exists so that a reviewer need not
    fetch anything.

    ```bash
    export RCOH_DATA=/path/to/deposit
    ```

---

## What is in the deposit

| file | rows | what it holds |
| --- | ---: | --- |
| `00_unit_manifest.csv` | 21,010 | burst and the four acquisition dates of every unit |
| `01_units_analysed.parquet` | 21,010 | the analysed sample, one row per unit |
| `02_units_by_landcover.parquet` | 221,229 | coherence resolved by cover class within each unit |
| `03_validation_arm.parquet` | 1,907 | dry-against-dry pairs; the zero point and the floor |
| `04_shared_control_arm.parquet` | 4,936 | rainfall and dry pairs against one shared control |
| `05_rain_masked.parquet` | 23,015 | coherence and rainfall over wetted pixels only |
| `06_perpendicular_baselines.parquet` | 89,069 | perpendicular baseline per interferogram |
| `07_within_burst_wet_dry.parquet` | 4,167 | wetted against unwetted ground in one image |
| `08_soil_water.parquet` | 25,768 | ERA5-Land soil water at all four acquisitions |
| `09_event_census.parquet` | 166,117 | every rainfall event screened, for the population weights |

Column-by-column definitions are in the [data dictionary](../reference/data.md).

---

## What is not in the deposit

!!! info "The interferograms"

    31,580 burst groups at roughly 75 MB each, some 300 GB. They are not
    deposited because they are regenerable, and because the Sentinel-1 bursts
    they derive from are Copernicus open data already archived at the Alaska
    Satellite Facility.

`00_unit_manifest.csv` is the bridge. It gives the burst identifier and the
four acquisition dates of every unit, so the whole store can be rebuilt:

```
unit_id                   burst            t_event     t_ref       t_sec       ctrl_ref    ctrl_sec
t001_000033_iw1_20150602  t001_000033_iw1  2015-06-02  2015-05-21  2015-06-02  2017-05-10  2017-05-22
```

See the [Pipeline overview](../reproduce/scratch.md) for the rebuild.

---

## Loading

```python
from rainfall_coherence import io

u = io.units()            # 21,010 rows, with reduction_pct and headroom
c = io.landcover()        # cover-resolved, restricted to analysed units
v = io.validation()       # dry-against-dry
s = io.soil_water()       # with dtheta oriented to the event
m = io.masked()           # loss_all and loss_wet side by side
w = io.within_burst()     # with the within-unit difference
b = io.baselines()        # perpendicular baselines
```

Each loader attaches the derived columns the analysis needs, so a caller never
recomputes a reduction by hand or merges the soil-water table unoriented.

!!! warning "Do not re-filter"

    `io.landcover()` deliberately does not apply the headroom condition, and
    `io.soil_water()` orients the change to the event. Both are described in
    [Pitfalls](../reproduce/pitfalls.md); getting either wrong changes published results.
