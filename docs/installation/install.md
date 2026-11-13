# Installation

The package splits into two halves. The analysis half reads the deposited
tables and needs nothing but the scientific Python stack. The extraction half
reads coherence rasters and talks to Earth Engine, and is only needed if you
intend to rebuild from the archive.

---

=== "Analysis only"

    Everything needed to reproduce the tables, the figures and Equation 3 from
    the data deposit.

    ??? note "A fresh environment"

        ```bash
        conda create -n rcoh python=3.12
        conda activate rcoh
        ```

    ```bash
    git clone https://github.com/jldz9/s1-rainfall-coherence
    cd s1-rainfall-coherence
    pip install -e .
    ```

    Pulls `numpy`, `pandas`, `pyarrow`, `scipy` and `matplotlib`.

=== "With extraction"

    Adds raster reading and Earth Engine, for phases 01 to 03 of the pipeline.

    ```bash
    pip install -e ".[extract]"
    ```

    ??? warning "rasterio and GDAL"

        `rasterio` ships manylinux wheels with GDAL bundled, so pip alone
        usually suffices. If it fails to build, install GDAL from conda first:

        ```bash
        conda install -c conda-forge gdal
        pip install -e ".[extract]"
        ```

---

## Verifying the install

```bash
rcoh --help
```

With the data in place (see [Getting the data](data.md)):

```bash
# nothing to set: data/ ships with the clone
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

If those match, the install reproduces Equation 3 of the paper.

### Checking every published number

```bash
pip install pytest
pytest
```

Eleven tests and twenty-nine assertions, each tied to a published number: the sample size, the
headline reduction, Equation 3 and its comparison of functional forms, five
rows of Table 4, the measurement floor, the within-burst difference, the
masking sensitivity, and that the soil water change is oriented to the event.
A pass means this install reproduces the published analysis exactly.

The suite runs against the bundled `data/`; set `RCOH_DATA` only to test a copy held elsewhere.

---

## Earth Engine

Only needed for `01_acquire/` and parts of `03_extract/`. Authenticate once:

```bash
earthengine authenticate
```

??? note "Which scripts need it"

    `build_events.py`, `build_covariates.py`, `extract_rainmasked.py`,
    `extract_rainfull.py` and `pull_cells.py`. Each takes a `--project`
    argument naming your Earth Engine cloud project.

---

## Requirements

| | |
| --- | --- |
| Python | ≥ 3.11 |
| analysis | numpy, pandas, pyarrow, scipy, matplotlib |
| extraction | rasterio, requests, earthengine-api |
| disk, analysis only | 25 MB for the deposit |
| disk, full rebuild | roughly 300 GB for the burst store |
