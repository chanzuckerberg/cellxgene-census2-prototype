"""Read the whole obs of one organism from Census 1 and print one JSON record.

Each run is a new process. The read is `cellxgene_census.get_obs`, with its defaults:
every column, and the TileDB context that cellxgene_census makes. The record splits
the time into the read to one Arrow table and the conversion to pandas, as get_obs does.

    python census1_obs.py homo_sapiens --census-version 2025-11-08
"""

import argparse
import importlib.metadata
import json
import platform
import resource
import time


def cpu() -> float:
    u = resource.getrusage(resource.RUSAGE_SELF)
    return u.ru_utime + u.ru_stime


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("organism")
    p.add_argument("--census-version", default="2025-11-08")
    p.add_argument("--repeat", type=int, default=0)
    args = p.parse_args()

    import cellxgene_census

    t0 = time.monotonic()
    census = cellxgene_census.open_soma(census_version=args.census_version)
    obs = census["census_data"][args.organism].obs
    t1, c1 = time.monotonic(), cpu()
    table = obs.read().concat()
    t2, c2 = time.monotonic(), cpu()
    frame = table.to_pandas()
    t3, c3 = time.monotonic(), cpu()
    peak = (
        resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20
    )  # Linux KiB -> GiB
    record = {
        "config": "census1",
        "census_version": args.census_version,
        "organism": args.organism,
        "repeat": args.repeat,
        "open_seconds": round(t1 - t0, 3),
        "read_arrow_seconds": round(t2 - t1, 3),
        "to_pandas_seconds": round(t3 - t2, 3),
        "query_seconds": round(t3 - t1, 3),
        "cpu_seconds": round(c3 - c1, 3),
        "read_arrow_cpu_seconds": round(c2 - c1, 3),
        "rows": len(frame),
        "columns": len(frame.columns),
        "arrow_gib": round(table.nbytes / 2**30, 2),
        "pandas_gib": round(frame.memory_usage(deep=False).sum() / 2**30, 2),
        "peak_rss_gib": round(peak, 2),
        "machine": platform.machine(),
        "versions": {
            m: importlib.metadata.version(m)
            for m in ("cellxgene-census", "tiledbsoma", "pandas", "pyarrow")
        },
    }
    print(json.dumps(record), flush=True)
    census.close()


if __name__ == "__main__":
    main()
