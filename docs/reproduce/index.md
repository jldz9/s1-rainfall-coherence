# Reproduce

Two routes to the same numbers. Pick by whether you want to check the analysis
or rebuild the measurement.

=== "With the provided data"

    **A minute, on a laptop.** Download the 25 MB deposit and every table,
    every fitted relationship and every figure in the paper regenerates.

    ```bash
    pip install -e .
    rcoh fit                 # Equation 3
    rcoh dose                # Table 3
    rcoh landcover           # Table 4
    rcoh climate             # Table 5
    rcoh floor               # the measurement floor
    ```

    And the figures:

    ```bash
    python -m rainfall_coherence.figures.fig_dose          
    python -m rainfall_coherence.figures.fig_paired_arms   
    python -m rainfall_coherence.figures.fig_one_sidedness 
    python -m rainfall_coherence.figures.fig_baseline      
    python -m rainfall_coherence.figures.fig_r12_dose      
    python -m rainfall_coherence.figures.fig_rain_vs_moisture
    ```

    This route starts from the measured coherences. It verifies the
    statistics, the fits and the figures, and it cannot verify the
    interferometry.

    :material-arrow-right: [Full walkthrough](deposit.md)

=== "From scratch"

    **Weeks, on a cluster.** Retrieve the Sentinel-1 bursts, form the
    interferograms, measure the coherence, and rebuild the tables the other
    route starts from.

    ```
    01_acquire    the burst grid, the S1 timeline, the rainfall census
         |
    02_process    unit selection, then ~57,000 interferograms on SLURM
         |
    03_extract    coherence over the rasters, soil water, baselines
         |
    04_tables     the tables that are deposited
         |
    (then the route on the left)
    ```

    Needs an Earth Engine account, a SLURM cluster, roughly 300 GB, and
    [InSARHub](https://github.com/jldz9/InSARHub) for the interferometry.

    :material-arrow-right: [Full walkthrough](scratch.md)

---

## Which route verifies what

| | with the deposit | from scratch |
| --- | :---: | :---: |
| the fitted relationships | ✅ | ✅ |
| the tables and figures | ✅ | ✅ |
| the cluster bootstrap and the intervals | ✅ | ✅ |
| that the coherences were measured correctly | ❌ | ✅ |
| that the units were selected correctly | ❌ | ✅ |
| that the rainfall was retrieved correctly | ❌ | ✅ |
| time | a minute | weeks |
| hardware | a laptop | a cluster |

!!! tip "A middle route"

    `00_unit_manifest.csv` in the deposit names the burst and the four
    acquisition dates of every unit, so you can rebuild the interferograms
    without repeating the event census or the unit selection. Start at
    [02 Process](process.md).

---

## Before either route

Read [Pitfalls](pitfalls.md). Three mistakes in handling these tables each
changed a published number during the analysis, and two of them are silent —
the result looks plausible and is wrong.
