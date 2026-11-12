# s1-rainfall-coherence

Code and data for *A global dose–response for rainfall and Sentinel-1
interferometric coherence*, Li et al.

Coherence loss was measured over **21,010 rainfall events**, each set against a
dry interferogram of the same burst, the same twelve-day baseline and the same
season. The relationship is logarithmic and has no threshold.

$$L = 0.82 + 5.27\,\log_2 R$$

for $L$ the percentage by which the coherence of the event pair falls below its
control and $R$ the 24 h accumulation in millimetres.

---

## What the study found

| | |
| --- | --- |
| mean reduction against a dry control | **11.8%**, 95% interval 11.2 to 12.4 |
| per doubling of the accumulation | **5.27 percentage points** |
| at 1–3 mm, where practice assumes nothing happens | **6.0%** |
| zero point, measured not assumed | **+0.17%** on dry-against-dry pairs |
| tree cover against grassland, same storm | **+3.4 points** |
| built-up against grassland | **−5.4 points** |

The headline is weighted to the census of 166,117 events the units were drawn
from, since every event above 15 mm was processed while lighter events were
sampled. Unweighted, the processed units give 19.2%, which describes the sample
rather than the rainfall.

The effect begins far below the accumulations usually treated as extreme. Half
of the maximum loss is reached by 10 mm, so routine light rainfall is a
systematic contributor to temporal decorrelation and is not currently screened
for in operational processing.

---

## Two ways to use this repository

=== "Reproduce the analysis"

    The 27 MB of tables ship with the repository. No download, no raw
    imagery, no Earth Engine account, no cluster.

    ```bash
    pip install -e .
    rcoh fit
    ```

    See [Getting the data](installation/data.md) and
    [Reproducing the paper](reproduce/deposit.md).

=== "Rebuild from the archive"

    Start from `00_unit_manifest.csv`, pull the Sentinel-1 bursts from ASF,
    form the interferograms, and run the full pipeline.

    This needs a cluster, an Earth Engine account and some weeks. See the
    [Pipeline overview](reproduce/scratch.md).

---

## Design in one paragraph

A **unit** is one rainfall event expressed as a comparison between two
interferograms of the same burst stack. The **case** pair spans the event, one
acquisition having received the rain and the other, twelve days away, less than
1 mm. The **control** pair is dry at both acquisitions and is drawn from the
same burst, the same baseline and the same season, so that the ground and the
interval are held in common and only the weather differs. Their ratio is the
unit's measurement of the loss.

!!! note "Why a ratio and not a difference"

    Coherence magnitude is specific to the processing chain. A ratio survives a
    change of scale between chains; an absolute difference in coherence does
    not. Neither survives a change of noise floor, so the reductions here are
    comparable with another study in their ordering and form, not in their
    magnitudes.

---

## Citing

Please cite the paper and the data deposit. `CITATION.cff` in the repository
root carries the machine-readable form.
