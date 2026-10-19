#!/usr/bin/env python3
"""Build the prototype catalog from a prefix of stores: one Parquet row for each store.

The catalog builder is not part of the SDK. The SDK reads a catalog, and the Census 2
team publishes one. This script is the prototype's stand-in for that builder. It reads
a prefix that holds `<dataset_version_id>.zarr` stores. Each store's summary record, in
the attributes of its root group, gives the dataset's identifiers and the gene reference
version. The rest of the row comes from the store's `obs` and `uns`. The builder reads
no receipt.

    python census2-sdk/tools/build_catalog.py s3://bucket/prefix/ catalog.parquet

A remote prefix needs the `s3` extra and credentials in the environment.
"""

import pathlib
import sys
from typing import Any

import pandas as pd
import zarr
from anndata.io import read_elem

from census2_sdk.catalog import INDEXED, is_remote, open_store
from census2_sdk.gene_reference import TAXON_TO_ORGANISM

# Assays that the curation schema lists as spatial.
SPATIAL_ASSAYS = frozenset({"EFO:0010961", "EFO:0030062"})

# The summary record that `census2-dataset-writer` writes into the attributes of
# each store's root group: `attrs[SUMMARY_NAMESPACE][SUMMARY_KEY]`. The SDK does
# not depend on the writer, so the names are repeated here. The writer's README
# describes the record.
SUMMARY_NAMESPACE = "cellxgene-census2"
SUMMARY_KEY = "summary"
SUMMARY_VERSION = 1


def _stores(prefix: str) -> list[str]:
    """Names of the `.zarr` stores directly under `prefix`."""
    if is_remote(prefix):
        import obstore
        import obstore.store

        store = obstore.store.from_url(prefix)
        listing = obstore.list_with_delimiter(store, "")
        names = [p.rstrip("/") for p in listing["common_prefixes"]]
    else:
        names = [p.name for p in pathlib.Path(prefix).iterdir() if p.is_dir()]
    return sorted(n for n in names if n.endswith(".zarr"))


def summary_record(group: zarr.Group, name: str) -> dict[str, Any]:
    """The summary record of the store `name`. Refuse a store that has no usable one.

    The writer adds the record last, so a store without one is partial. A catalog
    that dropped such a store quietly would drop its cells from every query.
    """
    namespace = group.attrs.get(SUMMARY_NAMESPACE)
    record: dict[str, Any] | None = (
        namespace.get(SUMMARY_KEY) if isinstance(namespace, dict) else None
    )
    if record is None:
        raise ValueError(f"{name}: no summary record; the store is partial")
    if record.get("version") != SUMMARY_VERSION:
        raise ValueError(
            f"{name}: summary record version {record.get('version')!r}, "
            f"and this builder reads version {SUMMARY_VERSION}"
        )
    if f"{record.get('dataset_version_id')}.zarr" != name:
        raise ValueError(
            f"{name}: the summary record names dataset version "
            f"{record.get('dataset_version_id')!r}"
        )
    return record


def build(prefix: str, out: str | pathlib.Path | None = None) -> pd.DataFrame:
    """One catalog row for each store under `prefix`.

    The store's summary record supplies `dataset_id`, `dataset_version_id`, and
    `gene_reference_version`. The store's `obs` and `uns` supply the rest.
    """
    prefix = prefix.rstrip("/")
    rows: list[dict[str, Any]] = []
    for name in _stores(prefix):
        store_uri = f"{prefix}/{name}"
        group = open_store(store_uri)
        record = summary_record(group, name)
        rows.append(
            {
                "dataset_id": str(record["dataset_id"]),
                "dataset_version_id": str(record["dataset_version_id"]),
                "store_uri": store_uri,
                "gene_reference_version": record.get("gene_reference_version"),
                **summarize(group),
            }
        )
    frame = pd.DataFrame(rows)
    if out is not None:
        frame.to_parquet(out, index=False)
    return frame


def summarize(group: zarr.Group) -> dict[str, Any]:
    """The catalog fields that one store supplies."""
    obs = group["obs"]
    assert isinstance(obs, zarr.Group)
    row: dict[str, Any] = {}
    n_obs: int | None = None
    # ponytail: list<string> columns; bitmaps over a dictionary (design 2.2) when the
    # catalog size matters.
    for column in INDEXED:
        if column not in obs:
            row[column] = None
            continue
        values = pd.Series(read_elem(obs[column]))
        n_obs = len(values)
        # A null is not a value. astype(str) would turn it into the string "nan".
        row[column] = sorted(values.dropna().astype(str).unique())
    if "is_primary_data" in obs:
        primary = pd.Series(read_elem(obs["is_primary_data"])).astype(bool)
        n_obs = len(primary)
        row["is_primary_data"] = (
            "all_true"
            if primary.all()
            else "all_false"
            if not primary.any()
            else "mixed"
        )
    else:
        row["is_primary_data"] = "mixed"
    if n_obs is None:
        n_obs = len(pd.Series(read_elem(obs[str(obs.attrs["_index"])])))
    row["n_obs"] = n_obs

    uns = group["uns"]
    assert isinstance(uns, zarr.Group)
    row["curation_schema_version"] = (
        str(read_elem(uns["schema_version"])) if "schema_version" in uns else None
    )
    taxon = str(read_elem(uns["organism_ontology_term_id"]))
    row["organism"] = TAXON_TO_ORGANISM.get(taxon, taxon)
    assays = row.get("assay_ontology_term_id") or []
    row["spatial"] = bool(SPATIAL_ASSAYS & set(assays))
    return row


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(__doc__, file=sys.stderr)
        return 2
    frame = build(argv[0], argv[1])
    print(f"{len(frame)} rows -> {argv[1]}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
