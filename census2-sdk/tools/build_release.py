#!/usr/bin/env python3
"""Build a release from the live catalog: one directory that pins the catalog.

The release builder is not part of the SDK. The SDK reads a release; the Census 2 team
publishes one (design 3.3). This script is the prototype stand-in for that builder.

    python census2-sdk/tools/build_release.py catalog.parquet releases/ 2026-09-01

The builder writes the directory `<releases>/<name>/`, on the local file system:

    catalog.parquet               the live catalog content, unchanged
    presence/<organism>.parquet   dataset_id, dataset_version_id, measured
    aggregates/var.parquet        organism, feature_id, nnz, n_measured_obs
    manifest.json                 the last write

The builder reads the live catalog one time and writes that content unchanged. It
derives each other file from the same content. The live catalog and the stores can be
remote. A remote location needs the `s3` extra and credentials in the environment.

Each presence file has one row for each store of the organism. `measured` lists the
positions in `var` of the genes that the dataset measured, from the `measured` column
of the store's `var`. The builder reads `var` and no expression data.

A position has one meaning in each store of an organism only when the stores hold the
same genes in the same order. Thus the builder compares the ordered gene identifiers
of each store with those of the first store of the organism. It stops at a store that
differs. For this check, the builder reads the index of `var` of each store, which is
one string for each gene (78,894 for a human store). It reads the index in the same
open that reads `measured`.

`aggregates/var.parquet` has one row for each gene of each organism. `nnz` is the
number of cells with a nonzero count of the gene. The builder adds the per-gene nnz
files that the writer makes next to the stores (design 3.4). `n_measured_obs` is the
number of cells of the datasets that measured the gene. The builder adds the arrays by
position, which the check of the gene identifiers makes valid. A store without its
file stops the build.

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
import posixpath
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
AGGREGATES = "aggregates/var.parquet"


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
    whose rows record two versions, or a row that records none. `open_manifest` of
    the SDK applies the same rule to the snapshot of a release.
    """
    versions: dict[str, str] = {}
    for organism, rows in frame.groupby("organism", sort=True):
        recorded = catalog.recorded(rows, "gene_reference_version")
        found = sorted({str(v) for v in recorded if v is not None})
        if None in recorded or len(found) != 1:
            raise ValueError(
                f"{organism}: a release holds one gene reference version for each "
                f"organism, and the catalog rows record {found}"
            )
        versions[str(organism)] = found[0]
    return versions


def read_var(
    store_uri: str, bucket_clients: catalog.Clients
) -> tuple[np.ndarray, np.ndarray]:
    """The gene identifiers of the `var` of one store, in order, and its `measured` flags.

    The function opens the store one time. It reads the index of `var`, which is one
    string for each gene, and the `measured` column, which is one flag for each gene.
    """
    var = catalog.open_store(store_uri, bucket_clients)["var"]
    if not isinstance(var, zarr.Group):
        raise TypeError(f"{store_uri}: var must be an encoded dataframe")
    if MEASURED not in var:
        raise ValueError(f"{store_uri}: var holds no {MEASURED} column")
    feature_ids = np.asarray(read_elem(var[str(var.attrs["_index"])])).astype(str)
    return feature_ids, np.asarray(read_elem(var[MEASURED]), dtype=bool)


def organism_var(
    organism: str, uris: list[str], bucket_clients: catalog.Clients
) -> tuple[np.ndarray, list[np.ndarray]]:
    """The gene identifiers of one organism, and the `measured` flags of each store.

    The stores of one organism must hold the same genes in the same order. A position
    in `var` then names one gene in each store (design 1.5). One gene reference
    version does not prove that order. Thus the function compares the
    identifiers of each store with those of the first store, and it refuses a store
    that differs. It keeps the identifiers of the first store only.
    """
    feature_ids, flags = read_var(uris[0], bucket_clients)
    measured = [flags]
    for uri in uris[1:]:
        other, flags = read_var(uri, bucket_clients)
        if not np.array_equal(other, feature_ids):
            raise ValueError(
                f"{uri}: the genes of var differ from {uris[0]}, the first store of "
                f"{organism}, in identifier or in order"
            )
        measured.append(flags)
    return feature_ids, measured


def presence(rows: pd.DataFrame, measured: list[np.ndarray]) -> pa.Table:
    """The presence table of one organism: the measured genes of each store.

    `measured` holds the flags of each row of `rows`, from `organism_var`. The
    stores hold the same genes in the same order, so a position names one gene in
    each store.
    """
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


def gene_nnz(
    store_uri: str, dataset_version_id: str, measured: np.ndarray, n_obs: int
) -> np.ndarray:
    """The array of the per-gene nnz file next to one store (design 3.4).

    The array must be `int64` with one value for each gene of `var`. A gene that the
    dataset did not measure must have zero, and no value can be above `n_obs`.
    """
    name = census2_dataset_writer.gene_nnz_name(dataset_version_id)
    uri = posixpath.join(posixpath.dirname(store_uri.rstrip("/")), name)
    try:
        content = catalog.read(uri)
    except FileNotFoundError as e:
        raise FileNotFoundError(
            f"{uri}: the store {store_uri} has no per-gene nnz file"
        ) from e
    try:
        # This call reads the `.npy` format only, and no pickle and no archive.
        nnz = np.lib.format.read_array(content, allow_pickle=False)
    except (ValueError, EOFError) as e:
        raise ValueError(f"{uri} is not a .npy array: {e}") from e
    if nnz.dtype != np.int64 or nnz.shape != measured.shape:
        raise ValueError(
            f"{uri} holds {nnz.dtype} of shape {nnz.shape}, and the file must hold "
            f"int64 for the {len(measured)} genes of var"
        )
    if nnz[~measured].any():
        raise ValueError(f"{uri} counts a gene that the dataset did not measure")
    if nnz.min(initial=0) < 0 or nnz.max(initial=0) > n_obs:
        raise ValueError(f"{uri} holds a count outside 0 to {n_obs}, the n_obs")
    return nnz


def aggregates(
    organism: str,
    rows: pd.DataFrame,
    feature_ids: np.ndarray,
    measured: list[np.ndarray],
) -> pa.Table:
    """The per-gene aggregates of one organism, one row for each gene of `var`.

    `nnz` is the sum of the per-gene nnz files of the stores. `n_measured_obs` is the
    sum of `n_obs` over the datasets that measured the gene. `feature_ids` and
    `measured` come from `organism_var`, and `measured` holds the flags of each row
    of `rows`.
    """
    # Each sum is at most the cells of the organism, so this limit excludes an overflow.
    if sum(int(n) for n in rows["n_obs"]) > np.iinfo(np.int64).max:
        raise ValueError(f"{organism}: the n_obs of the catalog rows exceed int64")
    nnz = np.zeros(len(feature_ids), dtype=np.int64)
    n_measured_obs = np.zeros(len(feature_ids), dtype=np.int64)
    for uri, version_id, n_obs, flags in zip(
        rows["store_uri"],
        rows["dataset_version_id"],
        rows["n_obs"],
        measured,
        strict=True,
    ):
        nnz += gene_nnz(str(uri), str(version_id), flags, int(n_obs))
        n_measured_obs += flags * int(n_obs)
    return pa.table(
        {
            "organism": pa.array([organism] * len(feature_ids), type=pa.string()),
            "feature_id": pa.array(feature_ids.tolist(), type=pa.string()),
            "nnz": pa.array(nnz, type=pa.int64()),
            "n_measured_obs": pa.array(n_measured_obs, type=pa.int64()),
        }
    )


def build(
    live_catalog: str | os.PathLike[str], releases: str | os.PathLike[str], name: str
) -> release.Manifest:
    """Write the release `name` of the live catalog under `releases`.

    Return the manifest. The steps follow design 3.3: the catalog snapshot, the
    presence file and the aggregates of each organism, and then the manifest.
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
    repeated = sorted(
        {str(d) for d in frame["dataset_id"][frame["dataset_id"].duplicated()]}
    )
    if repeated:
        # The SDK refuses a source that holds a dataset more than one time.
        raise ValueError(f"{live_catalog}: more than one row names {repeated}")
    versions = reference_versions(frame)
    for name_of_organism in versions:
        _require_name(name_of_organism, "organism")

    # Two builds of one name cannot both make the directory.
    directory.mkdir(parents=True)
    files = {release.CATALOG_NAME: _write(directory, release.CATALOG_NAME, content)}
    bucket_clients = catalog.clients(frame["store_uri"].astype(str))
    tables = []
    for organism, rows in frame.groupby("organism", sort=True):
        uris = [str(uri) for uri in rows["store_uri"]]
        feature_ids, measured = organism_var(str(organism), uris, bucket_clients)
        path = f"presence/{organism}.parquet"
        files[path] = _write(directory, path, _parquet(presence(rows, measured)))
        tables.append(aggregates(str(organism), rows, feature_ids, measured))
    files[AGGREGATES] = _write(
        directory, AGGREGATES, _parquet(pa.concat_tables(tables))
    )

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
