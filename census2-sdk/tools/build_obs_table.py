#!/usr/bin/env python3
"""Register the `obs` fragments of a catalog as one Iceberg table, in place.

The writer puts one Parquet fragment next to each store (`census2_dataset_writer.
store_format`). This tool defines the release-wide `obs` table of one organism as the
union of those fragments. It does not copy or rewrite a fragment: each fragment becomes
an Iceberg data file at its own location, through pyiceberg's `add_files`. The tool
writes only the Iceberg metadata, under the table location:

    <table location>/metadata/<n>-<uuid>.metadata.json   one file for each commit
    <table location>/metadata/version-hint.text          the name of the current one
    <table location>/metadata/snap-*.avro, *-m0.avro     the manifest list and manifests

    python census2-sdk/tools/build_obs_table.py catalog.parquet s3://bucket/census/obs/homo_sapiens

The tool needs no Iceberg catalog service. It opens a SQLite catalog in a temporary
directory, registers the current metadata file there when the table exists, and commits
one change: it adds the fragments of new catalog rows and deletes the data files of the
rows that the catalog no longer holds. A deletion drops whole files, because the table
is partitioned by `dataset_id` and each fragment holds one dataset. The commit writes a
new metadata file. The tool then writes `version-hint.text`, which is the commit point
for a reader of the live table. A release pins the metadata file itself.

The snapshot summary records the SHA-256 digest of the catalog file in
`census2.catalog-sha256`. The SDK reads the table only with that catalog, so the table
holds a fragment for each catalog row and nothing else.

The catalog row of each store gives the fragment's size and digest (`build_catalog.py`
copies them from the store's fragment record). The tool checks the size of each new
fragment, and with `--verify-digest` also its digest. A remote location needs
credentials in the environment and `AWS_REGION`.
"""

import argparse
import hashlib
import io
import pathlib
import sys
import tempfile
import time
from typing import Any, cast

import census2_dataset_writer
import pandas as pd
import pyarrow.fs as pafs  # type: ignore[import-untyped]  # pyarrow has no types
import pyarrow.parquet as pq  # type: ignore[import-untyped]
from census2_dataset_writer import fragment_location
from pyiceberg.catalog import Catalog
from pyiceberg.catalog.sql import SqlCatalog
from pyiceberg.partitioning import PartitionField, PartitionSpec
from pyiceberg.schema import assign_fresh_schema_ids
from pyiceberg.table import Table
from pyiceberg.transforms import IdentityTransform

from census2_sdk import expressions, obs_table
from census2_sdk.catalog import is_remote

NAMESPACE = "census2"

# The table properties. A string column's bounds are truncated to 16 characters, which
# is the Iceberg default; `observation_joinid` gets counts only, because its bounds
# prune nothing and would only make each manifest entry larger.
PROPERTIES = {
    "write.metadata.metrics.default": "truncate(16)",
    "write.metadata.metrics.column.observation_joinid": "counts",
    "format-version": "2",
}


def filesystem(location: str) -> tuple[Any, str]:
    """A pyarrow filesystem for `location`, and the path in it."""
    if is_remote(location):
        return cast("tuple[Any, str]", pafs.FileSystem.from_uri(location))
    return pafs.LocalFileSystem(), str(pathlib.Path(location).resolve())


def uri(location: str) -> str:
    """`location` as pyiceberg takes it: a URL, with `file://` for a local path."""
    if is_remote(location):
        return location.rstrip("/")
    return "file://" + str(pathlib.Path(location).resolve())


def fragments(catalog: pd.DataFrame, organism: str) -> pd.DataFrame:
    """The catalog rows of `organism`, with the location of each fragment.

    Refuse a row without a fragment record: the table must cover every catalog row,
    or a query through the table would drop the cells of that dataset.
    """
    rows = catalog[catalog["organism"] == organism].copy()
    if "obs_fragment_bytes" not in rows.columns:
        raise ValueError("the catalog has no obs fragment columns; rebuild it")
    missing = rows[rows["obs_fragment_bytes"].isna()]
    if len(missing):
        raise ValueError(f"no obs fragment for {sorted(missing['dataset_id'])}")
    # The fragment is next to the store, unless the row names it (`bench.fragments`).
    named = rows.get(
        "obs_fragment_uri", pd.Series([None] * len(rows), index=rows.index)
    )
    rows["fragment_uri"] = [
        uri(str(f) if isinstance(f, str) else fragment_location(str(s)))
        for s, f in zip(rows["store_uri"], named, strict=True)
    ]
    return rows


def check_fragments(rows: pd.DataFrame, verify_digest: bool) -> None:
    """Each fragment exists with the size, and optionally the digest, of its record."""
    row: Any
    for row in rows.itertuples():
        fs, path = filesystem(str(row.fragment_uri).removeprefix("file://"))
        info = fs.get_file_info(path)
        if info.type != pafs.FileType.File or info.size != int(row.obs_fragment_bytes):
            raise ValueError(
                f"{row.fragment_uri}: size {info.size}, the record gives "
                f"{int(row.obs_fragment_bytes)}"
            )
        if verify_digest:
            with fs.open_input_stream(path) as fh:
                digest = hashlib.sha256(fh.read()).hexdigest()
            if digest != row.obs_fragment_sha256:
                raise ValueError(f"{row.fragment_uri}: digest differs from its record")


def open_catalog(tmp: str, location: str) -> Catalog:
    """A SQLite Iceberg catalog in `tmp`, with the table's namespace."""
    catalog = SqlCatalog(
        "census2-build",
        uri=f"sqlite:///{tmp}/catalog.db",
        warehouse=location,
    )
    catalog.create_namespace(NAMESPACE)
    return catalog


def create(catalog: Catalog, name: str, location: str, first_fragment: str) -> Table:
    """A new table with the schema of `first_fragment`, partitioned by `dataset_id`."""
    fs, path = filesystem(first_fragment.removeprefix("file://"))
    # Field IDs from 1, in column order. The partition field names its source by ID.
    schema = assign_fresh_schema_ids(
        Catalog._convert_schema_if_needed(pq.read_schema(path, filesystem=fs))
    )
    source = schema.find_field(census2_dataset_writer.DATASET_ID_COLUMN)
    spec = PartitionSpec(
        PartitionField(
            source_id=source.field_id,
            field_id=1000,
            transform=IdentityTransform(),
            name=census2_dataset_writer.DATASET_ID_COLUMN,
        )
    )
    return catalog.create_table(
        (NAMESPACE, name),
        schema=schema,
        location=location,
        partition_spec=spec,
        properties=PROPERTIES,
    )


def build(
    catalog_path: str,
    location: str,
    organism: str = "homo_sapiens",
    *,
    verify_digest: bool = False,
    tag: str | None = None,
) -> dict[str, Any]:
    """Commit the fragments of the catalog's `organism` rows to the table at `location`.

    Return a summary: the metadata file, the counts of added, deleted, and kept
    fragments, and the elapsed time of each step.
    """
    t0 = time.monotonic()
    content = obs_table.read_bytes(catalog_path)
    catalog_sha256 = hashlib.sha256(content).hexdigest()
    rows = fragments(pd.read_parquet(io.BytesIO(content)), organism)
    location = uri(location)
    current = obs_table.current_metadata(location)

    with tempfile.TemporaryDirectory() as tmp:
        iceberg = open_catalog(tmp, location)
        table: Table | None = None
        existing: dict[str, str] = {}
        if current is None:
            if not len(rows):
                raise ValueError(f"the catalog has no {organism} rows")
        else:
            table = iceberg.register_table((NAMESPACE, organism), current)
            # Each data file is one fragment, so its partition is its dataset_id.
            existing = {
                task.file.file_path: str(task.file.partition[0])
                for task in table.scan().plan_files()
            }
        added = rows[~rows["fragment_uri"].isin(set(existing))]
        removed = sorted(set(existing) - set(rows["fragment_uri"]))
        t1 = time.monotonic()
        check_fragments(added, verify_digest)
        t2 = time.monotonic()
        if table is None:
            table = create(
                iceberg, organism, location, str(rows["fragment_uri"].iloc[0])
            )

        properties = {"census2.catalog-sha256": catalog_sha256}
        with table.transaction() as tx:
            if removed:
                gone = sorted({existing[path] for path in removed})
                column = census2_dataset_writer.DATASET_ID_COLUMN
                tx.delete(expressions.is_in(column, gone), properties)
            if len(added):
                tx.add_files(list(added["fragment_uri"]), properties)
            if not removed and not len(added):
                # Nothing changed, but the table must name this catalog.
                tx.set_properties({"census2.catalog-sha256": catalog_sha256})
        table = iceberg.load_table((NAMESPACE, organism))
        if tag is not None and table.current_snapshot() is not None:
            snapshot = table.current_snapshot()
            assert snapshot is not None
            table.manage_snapshots().create_tag(snapshot.snapshot_id, tag).commit()
            table = iceberg.load_table((NAMESPACE, organism))
        metadata = table.metadata_location
        t3 = time.monotonic()
    obs_table.write_version_hint(location, metadata)
    return {
        "metadata": metadata,
        "catalog_sha256": catalog_sha256,
        "fragments": len(rows),
        "added": len(added),
        "deleted": len(removed),
        "plan_seconds": round(t1 - t0, 3),
        "check_seconds": round(t2 - t1, 3),
        "commit_seconds": round(t3 - t2, 3),
    }


def main(argv: list[str]) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("catalog", help="the catalog Parquet file, a path or a URL")
    p.add_argument("location", help="the table location, a path or a URL")
    p.add_argument("--organism", default="homo_sapiens")
    p.add_argument("--verify-digest", action="store_true")
    p.add_argument("--tag", default=None, help="tag the new snapshot, e.g. a release")
    args = p.parse_args(argv)
    summary = build(
        args.catalog,
        args.location,
        args.organism,
        verify_digest=args.verify_digest,
        tag=args.tag,
    )
    for key, value in summary.items():
        print(f"{key}: {value}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
