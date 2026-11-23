"""Split the cost of one store open into its parts, one store after another.

Usage: python opencost.py CATALOG [--consolidated]

For each store that a read of the catalog opens (the per-dataset stores, or with
--consolidated the stores that serve the datasets), measure in this process:

- client: the CPU time of `obstore.store.from_url` for the store
- parse: the CPU time of `zarr.open_group` over a memory store that already holds
  the root `zarr.json`, so no request is made
- sdk: the wall and CPU time of the SDK's open (`catalog.open_store`)
- lean: the wall and CPU time of an open with `zarr_format=3` through one shared
  client for the bucket

and the size of the root `zarr.json`. Print one JSON object with the medians and
the totals.
"""

import json
import statistics
import sys
import time

import obstore
import obstore.store
import pandas as pd
import zarr
import zarr.storage
from census2_sdk import catalog as _catalog
from zarr.core.buffer import default_buffer_prototype


def cpu() -> float:
    return time.process_time()


def uris_of(path: str, consolidated: bool) -> list[str]:
    rows = _catalog.load(path)
    uris = [str(u) for u in rows["store_uri"]]
    if consolidated:
        uris = [
            str(c) if pd.notna(c) else u
            for c, u in zip(rows["consolidated_store_uri"], uris, strict=True)
        ]
    return list(dict.fromkeys(uris))


def main() -> None:
    path = sys.argv[1]
    consolidated = "--consolidated" in sys.argv
    uris = uris_of(path, consolidated)
    bucket = uris[0].split("://", 1)[1].split("/", 1)[0]
    shared = zarr.storage.ObjectStore(
        obstore.store.from_url(f"s3://{bucket}"), read_only=True
    )
    parts: dict[str, list[float]] = {
        k: []
        for k in (
            "client_cpu",
            "parse_cpu",
            "sdk_wall",
            "sdk_cpu",
            "lean_wall",
            "lean_cpu",
            "json_kb",
        )
    }
    # One warm-up open of each kind, so imports and the first TLS setup do not count.
    _catalog.open_store(uris[0])
    key0 = uris[0].split(bucket + "/", 1)[1].rstrip("/")
    zarr.open_group(shared, path=key0, mode="r", zarr_format=3)
    for uri in uris:
        key = uri.split(bucket + "/", 1)[1].rstrip("/")
        c0 = cpu()
        obstore.store.from_url(uri)
        parts["client_cpu"].append(cpu() - c0)

        raw = bytes(obstore.get(shared.store, f"{key}/zarr.json").bytes())
        parts["json_kb"].append(len(raw) / 1024)
        memory = zarr.storage.MemoryStore(
            {"zarr.json": default_buffer_prototype().buffer.from_bytes(raw)}
        )
        c0 = cpu()
        zarr.open_group(memory, mode="r", zarr_format=3)
        parts["parse_cpu"].append(cpu() - c0)

        t0, c0 = time.perf_counter(), cpu()
        _catalog.open_store(uri)
        parts["sdk_wall"].append(time.perf_counter() - t0)
        parts["sdk_cpu"].append(cpu() - c0)

        t0, c0 = time.perf_counter(), cpu()
        zarr.open_group(shared, path=key, mode="r", zarr_format=3)
        parts["lean_wall"].append(time.perf_counter() - t0)
        parts["lean_cpu"].append(cpu() - c0)
    out = {"catalog": path, "consolidated": consolidated, "stores": len(uris)}
    for name, values in parts.items():
        scale = 1 if name == "json_kb" else 1000
        unit = "" if name == "json_kb" else "_ms"
        out[f"{name}{unit}_median"] = round(statistics.median(values) * scale, 2)
        out[f"{name}{unit}_total"] = round(sum(values) * scale, 1)
    print(json.dumps(out))


if __name__ == "__main__":
    main()
