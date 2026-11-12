"""Constants fixed for the whole study.

These were duplicated across a dozen analysis scripts, which is how the noise
floor came to be applied in a place it does not belong.  They live here now and
are imported everywhere.
"""
from __future__ import annotations

# ---------------------------------------------------------------- estimator --
NOISE_FLOOR = 0.077
"""Coherence returned by the estimator where no correlation survives.

A phase-only estimator over an 11 x 11 Gaussian window of standard deviation
one third the window width.  The value is empirical, confirmed against the
distribution of observed coherences.  It is NOT portable: a chain with a
different window or a maximum-likelihood estimator has a different floor.
"""

MIN_HEADROOM = 0.15
"""A control must stand this far above the floor for a loss to be resolvable.

Equivalently a control coherence of 0.227 or more, which leaves room for a
reduction of two thirds before the case reaches the floor.

This is a property of the UNIT, applied once to the burst mean when the
analysed sample is built.  Applying it again to the rows of the land-cover
table discards the dim covers and reverses the sign of the wetland contrast.
See `tables.landcover_contrasts`.
"""

LOOKS_DISPLACEMENT = 4
"""Multilooking assumed when converting coherence to phase standard deviation.

A convention for a displacement product, not the number of looks used to
estimate the coherence itself.
"""

WAVELENGTH_MM = 55.5
"""Sentinel-1 C-band wavelength, for converting phase noise to displacement."""

# ----------------------------------------------------------------- rainfall --
WET_DAY_MM = 1.0
"""Threshold separating a wet day from a dry one, in millimetres over 24 h.

The ETCCDI wet-day definition (Zhang et al. 2011).  One value serves three
distinct purposes, which is why the aliases below exist:

* the accumulation an event must reach to enter the study at all,
* the criterion deciding whether an acquisition is dry enough for a control,
* the spatial mask separating wetted ground from dry within a burst.
"""

EVENT_FLOOR_MM = WET_DAY_MM
"""The smallest accumulation counted as an event.  Same value, first role."""

DRY_ACQ_MM = WET_DAY_MM
"""An acquisition below this is dry enough to enter a control pair."""

FULL_PROCESSING_MM = 15.0
"""Above this every available event was processed; below it events were sampled.

This is a statement about *coverage*, not about rainfall.  Of the 19,918 events
reaching 15 mm, 17,436 were processed, the rest failing on the availability of
one of their four acquisitions.  Below 15 mm events are plentiful -- 146,199 of
them -- so the number processed in each band was governed by the precision
sought rather than by the supply, and 13,157 were drawn at random within
logarithmically spaced bands.

Beware `build_campaign.py`, which names this same 15.0 `WET_MM`.  That script
was written for the heavy arm, where 15 mm was the threshold a case acquisition
had to reach; it is not the wet-day threshold of the paper and it is not an
event floor.
"""

# retained for backward compatibility with scripts that import WET_MM
WET_MM = WET_DAY_MM
HEAVY_MM = FULL_PROCESSING_MM

DOSE_EDGES = [1, 3, 6, 10, 15, 20, 25, 30, 40, 50, 70, 100, 320]
"""Accumulation bands of Table 2, in millimetres.

Used by the dose table and by the figures, so that a figure and the table
beside it are the same arithmetic.
"""

# ------------------------------------------------------------------- figure --
COLOUR = {
    "rain": "#9E2A2B",
    "control": "#335C81",
    "within": "#1B7A6B",
    "grey": "0.55",
}

PLOT_STYLE = {
    "font.size": 8,
    "font.family": "DejaVu Sans",
    "axes.linewidth": 0.7,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "savefig.dpi": 420,
    "savefig.bbox": "tight",
    "legend.frameon": False,
}
