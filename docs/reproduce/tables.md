# 04 — Tables

Assembly of the deposited tables from the extraction outputs.

| script | writes |
| --- | --- |
| `build_paired_arms.py` | `units_final`, `validation_arm`, `paired_rain_vs_null` |
| `build_validation_arm.py` | the dry-against-dry pairs |
| `merge_and_analyse.py` | `dgamma_all`, the cover-resolved table |
| `soil_moisture_arms.py` | the soil-water join |
| `build_release_all.py` | the release tables |
| `landcover_rain.py` | the within-burst cover contrasts |
| `cell_counts.py` | per-cell diagnostics |

---

## The three validation arms

The paper carries three independent constructions, all built here, and they
count different objects. This is the commonest source of confusion:

| arm | n | counts |
| --- | ---: | --- |
| validation arm | 1,907 | **pairs**, one per burst |
| shared control | 4,936 | **rainfall units**, 1,832 bursts |
| unrestricted substitution | 64,074 | **pairings**, 1,751 bursts |

The validation arm takes one pair per burst because pairs from a burst share
ground and often an interferogram; seven controls would yield twenty-one pairs
resting on seven measurements, and counting them all would narrow the
bootstrap interval beyond what the evidence supports.

The shared-control arm is larger because several rainfall units can lean on
the same control — a median of two, and as many as fourteen.

---

## Once these have run

Everything downstream is in the package rather than in scripts:

```python
from rainfall_coherence import tables
tables.dose()
tables.landcover_contrasts()
tables.dose_relationship()
```

See [Reproducing the paper](deposit.md).
