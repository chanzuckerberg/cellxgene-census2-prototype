#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.10,<3.14"  # tiledbsoma has no 3.14 wheel
# dependencies = ["cellxgene-census"]
# ///
"""Sample datasets from a Census 1 release and write an ingestion manifest.

The manifest names the dataset versions that the Census 1 release holds, not the
current Discover versions. Stores made from it therefore hold the same source
data as the release, and a query against them can be compared with the same
query against Census 1.

The sample is uniform over the datasets that contribute cells to one organism's
experiment in the release. The records are in selection order, so each prefix of
the manifest is also a uniform sample.

A schema migration gives every dataset a new version, so the versions in an older
release are usually not current. Discover keeps the older h5ads in its asset
bucket, under their own dataset_version_id, and this script gets their metadata
from the curation API. It then asserts that each file size agrees with the copy
of the h5ad that the release keeps.

Unlike tools/cxg-manifest.py, this script needs cellxgene-census, so run it with
uv, which reads the dependency from the header above.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import random
import sys
import urllib.request
from pathlib import Path

import cellxgene_census
import pyarrow.compute as pc


def load_cxg_manifest():
    """Import tools/cxg-manifest.py, so the records have the same shape."""
    path = Path(__file__).with_name("cxg-manifest.py")
    spec = importlib.util.spec_from_file_location("cxg_manifest", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def newest_lts() -> str:
    versions = cellxgene_census.get_census_version_directory(lts=True)
    return max(versions)


def organism_datasets(census, organism: str) -> set[str]:
    """Get the IDs of the datasets that have cells in the organism's experiment."""
    obs = census["census_data"][organism].obs
    ids: set[str] = set()
    for tbl in obs.read(column_names=["dataset_id"]):
        ids.update(pc.unique(tbl["dataset_id"]).to_pylist())
    return ids


def http_get_json(url: str, user_agent: str) -> object:
    req = urllib.request.Request(url, headers={"User-Agent": user_agent})
    with urllib.request.urlopen(req, timeout=120) as resp:
        return json.load(resp)


def http_content_length(url: str, user_agent: str) -> int:
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": user_agent})
    with urllib.request.urlopen(req, timeout=120) as resp:
        return int(resp.headers["Content-Length"])


def dataset_version(api_base: str, row: dict, user_agent: str) -> dict:
    """Get the curation API metadata of the dataset version that the release holds."""
    url = f"{api_base}/curation/v1/datasets/{row['dataset_id']}/versions"
    versions = http_get_json(url, user_agent)
    for v in versions:
        if v["dataset_version_id"] == row["dataset_version_id"]:
            return v
    raise RuntimeError(
        f"dataset {row['dataset_id']}: Discover has no version "
        f"{row['dataset_version_id']}"
    )


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--census-version",
        default=None,
        help="Census 1 release, for example 2025-11-08 (default: the newest LTS)",
    )
    p.add_argument(
        "--organism",
        default="homo_sapiens",
        help="Census experiment name (default homo_sapiens)",
    )
    p.add_argument("--limit", type=int, default=100, help="0 means all")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--out", default="manifest.jsonl")
    args = p.parse_args()

    cxg = load_cxg_manifest()
    version = args.census_version or newest_lts()
    h5ads = cellxgene_census.get_census_version_description(version)["h5ads"]
    bucket = h5ads["uri"].removeprefix("s3://").split("/", 1)[0]
    h5ads_url = (
        f"https://{bucket}.s3.{h5ads['s3_region']}.amazonaws.com{h5ads['relative_uri']}"
    )

    with cellxgene_census.open_soma(census_version=version) as census:
        datasets = census["census_info"]["datasets"].read().concat().to_pylist()
        members = organism_datasets(census, args.organism)
    print(
        f"sample: Census {version} holds {len(datasets)} datasets, "
        f"{len(members)} in {args.organism}",
        file=sys.stderr,
    )

    by_id = {d["dataset_id"]: d for d in datasets}
    population = sorted(members)
    limit = args.limit or len(population)
    chosen = random.Random(args.seed).sample(population, limit)

    records = []
    for dataset_id in chosen:
        row = by_id[dataset_id]
        record = cxg.to_record(
            dataset_version(cxg.API_BASE, row, cxg.USER_AGENT), cxg.SOURCE_BUCKET
        )
        census_size = http_content_length(
            h5ads_url + row["dataset_h5ad_path"], cxg.USER_AGENT
        )
        if census_size != record["filesize"]:
            raise RuntimeError(
                f"dataset {dataset_id}: Discover gives {record['filesize']} bytes "
                f"for version {row['dataset_version_id']}, and Census {version} "
                f"keeps {census_size}"
            )
        records.append(record)

    with open(args.out, "w") as fh:
        fh.writelines(json.dumps(r) + "\n" for r in records)

    total = sum(r["filesize"] for r in records)
    print(
        f"sample: wrote {len(records)} records to {args.out} "
        f"({total / 1024**3:.2f} GiB, seed={args.seed})",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
