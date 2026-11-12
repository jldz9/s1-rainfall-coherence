"""Command line entry point: `rcoh <table|fit|floor> [--data DIR]`."""
from __future__ import annotations

import argparse

from . import tables


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="rcoh", description=__doc__)
    p.add_argument("what", choices=["dose", "landcover", "climate", "fit", "floor",
                                    "population"])
    p.add_argument("--data", default=None, help="the deposit directory")
    a = p.parse_args(argv)

    if a.what == "dose":
        print(tables.dose(a.data).to_string(index=False))
    elif a.what == "landcover":
        print(tables.landcover_contrasts(a.data).to_string(index=False))
    elif a.what == "climate":
        k, s = tables.climate_and_season(a.data)
        print(k.to_string(index=False), "\n")
        print(s.to_string(index=False))
    elif a.what == "fit":
        for k, v in tables.dose_relationship(a.data).items():
            print(f"  {k:20} {v}")
    elif a.what == "floor":
        for k, v in tables.measurement_floor(a.data).items():
            print(f"  {k:20} {v}")
    elif a.what == "population":
        for k, v in tables.population_estimate(a.data).items():
            print(f"  {k:24} {v}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
