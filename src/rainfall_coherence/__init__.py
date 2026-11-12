"""Rainfall and Sentinel-1 interferometric coherence.

Code supporting "A global dose-response for rainfall and Sentinel-1
interferometric coherence", Li et al.

    import os
    os.environ["RCOH_DATA"] = "/path/to/the/zenodo/deposit"

    from rainfall_coherence import io, tables
    tables.dose()                   # Table 2
    tables.landcover_contrasts()    # Table 3
    tables.dose_relationship()      # Equation (3)
"""
from . import config, io, stats, tables   # noqa: F401

__version__ = "1.0.0"
__all__ = ["config", "io", "stats", "tables", "figures"]
