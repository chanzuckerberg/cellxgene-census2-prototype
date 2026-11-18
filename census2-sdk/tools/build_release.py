#!/usr/bin/env python3
"""Build a release from the live catalog: one directory that pins the catalog.

The release builder is not part of the SDK. The SDK reads a release; the Census 2 team
publishes one (design 3.3). This script is the prototype stand-in for that builder.

    python census2-sdk/tools/build_release.py catalog.parquet releases/ 2026-09-01

The builder writes the directory `<releases>/<name>/`, on the local file system:

    catalog.parquet               the live catalog content, unchanged
    presence/<organism>.parquet   dataset_id, dataset_version_id, measured
    manifest.json                 the last write

The builder reads the live catalog one time and writes that content unchanged. It
derives each other file from the same content. The live catalog and the stores can be
remote. A remote location needs the `s3` extra and credentials in the environment.

Each presence file has one row for each store of the organism. `measured` lists the
positions in `var` of the genes that the dataset measured, from the `measured` column
of the store's `var`. The builder reads `var` and no expression data.

The builder writes `manifest.json` last, so a directory with a manifest is a complete
release. `census2_sdk.release` defines the manifest. Each release is a new directory,
and the builder refuses a directory that exists. A build that stopped leaves a
directory without a manifest. Remove that directory before the next build.

This version makes no dictionary file, no check of a store against its inventory, and
no entry in the release index. The live catalog has no dictionary, and a store has no
inventory.
"""

import hashlib
import io
import json
import os
import pathlib
import sys

import census2_dataset_writer
import numpy as np
import pandas as pd
import pyarrow as pa  # type: ignore[import-untyped]  # pyarrow has no types
import pyarrow.parquet as pq  # type: ignore[import-untyped]
import zarr
from anndata.io import read_elem

from census2_sdk import catalog, release

MEASURED = census2_dataset_writer.MEASURED_COLUMN


def _require_name(value: str, kind: str) -> None:
    """Refuse a value that cannot be the name of a file or a directory."""
    if not value or "/" in value or value in (".", ".."):
        raise ValueError(f"the {kind} {value!r} is not a file name")


def _write(directory: pathlib.Path, path: str, data: bytes) -> str:
    """Write one file of the release. Return its SHA-256 for the manifest."""
    target = directory / path
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(data)
    return hashlib.sha256(data).hexdigest()


def _parquet(table: pa.Table) -> bytes:
    content = io.BytesIO()
    pq.write_table(table, content)
    return content.getvalue()


def reference_versions(frame: pd.DataFrame) -> dict[str, str]:
    """The gene reference version of each organism of the catalog (design 1.5).

    A release holds one version for each organism. The builder refuses an organism
    whose rows record two versions, or a row that records none.
    """
    versions: dict[str, str] = {}
    for organism, rows in frame.groupby("organism", sort=True):
        recorded = rows["gene_reference_version"]
        found = sorted({str(v) for v in recorded.dropna()})
        if bool(recorded.isna().any()) or len(found) != 1:
            raise ValueError(
                f"{organism}: a release holds one gene reference version for each "
                f"organism, and the catalog rows record {found}"
            )
        versions[str(organism)] = found[0]
    return versions


def measured_genes(store_uri: str, bucket_clients: catalog.Clients) -> np.ndarray:
    """The `measured` column of the `var` of one store: one flag for each gene."""
    var = catalog.open_store(store_uri, bucket_clients)["var"]
    if not isinstance(var, zarr.Group):
        raise TypeError(f"{store_uri}: var must be an encoded dataframe")
    if MEASURED not in var:
        raise ValueError(f"{store_uri}: var holds no {MEASURED} column")
    return np.asarray(read_elem(var[MEASURED]), dtype=bool)


def presence(organism: str, rows: pd.DataFrame, measured: list[np.ndarray]) -> pa.Table:
    """The presence table of one organism: the measured genes of each store.

    `measured` holds the flags of each row of `rows`, from `measured_genes`. The
    stores of one organism have one `var`, so the positions of one store are valid
    in each other store. The builder refuses stores whose `var` differ in length.
    """
    lengths = sorted({len(flags) for flags in measured})
    if len(lengths) != 1:
        raise ValueError(f"{organism}: the stores have var of {lengths} genes")
    positions = [np.flatnonzero(flags) for flags in measured]
    return pa.table(
        {
            "dataset_id": pa.array(rows["dataset_id"].astype(str), type=pa.string()),
            "dataset_version_id": pa.array(
                rows["dataset_version_id"].astype(str), type=pa.string()
            ),
            "measured": pa.array(positions, type=pa.list_(pa.int32())),
        }
    )


def build(
    live_catalog: str | os.PathLike[str], releases: str | os.PathLike[str], name: str
) -> release.Manifest:
    """Write the release `name` of the live catalog under `releases`.

    Return the manifest. The steps follow design 3.3: the catalog snapshot, the
    presence file of each organism, and then the manifest.
    """
    _require_name(name, "release name")
    directory = pathlib.Path(releases) / name
    if directory.exists():
        raise FileExistsError(f"{directory} exists, and a release is a new directory")

    # One read, so that each file of the release comes from one catalog version.
    content = catalog.read(live_catalog).getvalue()
    frame = catalog.load(io.BytesIO(content))
    if frame.empty:
        raise ValueError(f"{live_catalog}: the live catalog has no rows")
    unnamed = int(frame["organism"].isna().sum())
    if unnamed:
        # A group step drops these rows, and the presence files then lack their stores.
        raise ValueError(f"{live_catalog}: {unnamed} rows record no organism")
    versions = reference_versions(frame)
    for name_of_organism in versions:
        _require_name(name_of_organism, "organism")

    # Two builds of one name cannot both make the directory.
    directory.mkdir(parents=True)
    files = {release.CATALOG_NAME: _write(directory, release.CATALOG_NAME, content)}
    bucket_clients = catalog.clients(frame["store_uri"].astype(str))
    for organism, rows in frame.groupby("organism", sort=True):
        measured = [
            measured_genes(uri, bucket_clients) for uri in rows["store_uri"].astype(str)
        ]
        path = f"presence/{organism}.parquet"
        table = presence(str(organism), rows, measured)
        files[path] = _write(directory, path, _parquet(table))

    manifest = release.validate_manifest(
        {
            "version": release.MANIFEST_VERSION,
            "release": name,
            "gene_reference_version": versions,
            "stores": [
                {"dataset_version_id": str(v)} for v in frame["dataset_version_id"]
            ],
            "files": files,
        }
    )
    # A rename is atomic, so `manifest.json` is complete when it exists.
    partial = directory / f"{release.MANIFEST_NAME}.partial"
    partial.write_text(json.dumps(manifest, indent=1) + "\n", encoding="utf-8")
    partial.replace(directory / release.MANIFEST_NAME)
    return manifest


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print(__doc__, file=sys.stderr)
        return 2
    manifest = build(argv[0], argv[1], argv[2])
    directory = pathlib.Path(argv[1]) / argv[2]
    print(f"{len(manifest['stores'])} stores -> {directory}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
