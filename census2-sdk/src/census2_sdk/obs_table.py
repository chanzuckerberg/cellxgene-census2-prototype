"""The release-wide `obs` table: one Iceberg table over the `obs` fragments of the stores.

The writer puts one Parquet fragment next to each store, with one row for each cell: the
store's `obs` columns, the store row ordinal (`store_row`), and the range of the cell's
values in `X` (`x_start`, `x_end`). `tools/build_obs_table.py` registers the fragments
of a catalog as the data files of one Iceberg table for each organism, without a copy.
It writes the table metadata under the table location and names the current metadata
file in `metadata/version-hint.text`.

A query scans the table once. The scan returns the cell key, the store row ordinals,
and the `obs` columns of the result. The SDK then opens only the stores that hold a
selected cell, and reads the `X` ranges of those cells. It does not read a store's
`obs`. See `Census.get_anndata` and `Census.get_obs`.

The table is a static table: the SDK reads a metadata file and needs no Iceberg catalog
service. `open_table` takes a metadata file, or a table location and its version hint.
The snapshot summary of the table names the catalog that it was built from, and
`check_catalog` refuses another catalog. Without that check, a catalog row whose
fragment is not in the table would lose its cells from every result.

pyiceberg is the optional `iceberg` extra. pyiceberg ships type information; pyarrow
does not.
"""

import hashlib
import io
import logging
import os
import pathlib
import time
from dataclasses import dataclass, field
from typing import Any

import pyarrow as pa  # type: ignore[import-untyped]  # pyarrow has no types
from pyiceberg.expressions import BooleanExpression
from pyiceberg.io.pyarrow import ArrowScan
from pyiceberg.table import StaticTable

from census2_sdk import catalog as _catalog

# The columns that a fragment adds to the store's obs columns. They agree with
# census2_dataset_writer.store_format, which the SDK does not import.
STORE_ROW = "store_row"
X_START = "x_start"
X_END = "x_end"
ROW_COLUMNS = (STORE_ROW, X_START, X_END)

# The snapshot summary property that names the catalog of the table.
CATALOG_SHA256 = "census2.catalog-sha256"

HINT = "metadata/version-hint.text"

# pyiceberg logs a warning for each file that it reads with a dictionary-encoded
# column. The SDK asks for those columns, so the warning carries no information.
logging.getLogger("pyiceberg.io.pyarrow").addFilter(
    lambda record: "does not have a dictionary type" not in record.getMessage()
)


def _properties() -> dict[str, str]:
    """The pyiceberg file IO properties: the region, from the environment."""
    region = os.environ.get("AWS_REGION") or os.environ.get("AWS_DEFAULT_REGION")
    return {"s3.region": region} if region else {}


def read_bytes(uri: str | os.PathLike[str]) -> bytes:
    """The content of one file, local or remote."""
    return _catalog.read(uri).getvalue()


def _location(location: str) -> str:
    """A table location as a URL. A local path becomes a `file://` URL."""
    if _catalog.is_remote(location) or location.startswith("file://"):
        return location.rstrip("/")
    return "file://" + str(pathlib.Path(location).resolve())


def current_metadata(location: str) -> str | None:
    """The current metadata file of the table at `location`, or None when it has none."""
    location = _location(location)
    hint = f"{location}/{HINT}"
    try:
        content = read_bytes(hint.removeprefix("file://")).decode().strip()
    except (FileNotFoundError, OSError):
        return None
    if not content:
        return None
    return f"{location}/metadata/{content}"


def write_version_hint(location: str, metadata: str) -> None:
    """Name `metadata` as the current metadata file of the table at `location`.

    The hint holds the file name. One small object write is the commit point for a
    reader of the live table.
    """
    location = _location(location)
    name = metadata.rstrip("/").rsplit("/", 1)[1]
    if f"{location}/metadata/{name}" != metadata:
        raise ValueError(f"{metadata} is not in {location}/metadata")
    if _catalog.is_remote(location):
        _catalog._require_obstore()
        import obstore
        import obstore.store

        obstore.put(
            obstore.store.from_url(f"{location}/metadata"),
            "version-hint.text",
            name.encode(),
        )
        return
    path = pathlib.Path(location.removeprefix("file://")) / HINT
    path.write_text(name)


def open_table(uri: str) -> StaticTable:
    """The table of a metadata file, or of a table location through its version hint."""
    if uri.endswith(".metadata.json"):
        metadata = uri
    else:
        found = current_metadata(uri)
        if found is None:
            raise FileNotFoundError(f"no {HINT} under {uri}")
        metadata = found
    if not _catalog.is_remote(metadata) and not metadata.startswith("file://"):
        metadata = "file://" + str(pathlib.Path(metadata).resolve())
    return StaticTable.from_metadata(metadata, properties=_properties())


def catalog_sha256(content: io.BytesIO) -> str:
    return hashlib.sha256(content.getvalue()).hexdigest()


def check_catalog(table: StaticTable, sha256: str, source: str) -> None:
    """Refuse a table that was not built from the catalog with this digest."""
    snapshot = table.current_snapshot()
    recorded = None
    if snapshot is not None and snapshot.summary is not None:
        recorded = snapshot.summary.additional_properties.get(CATALOG_SHA256)
    if recorded is None:
        recorded = table.properties.get(CATALOG_SHA256)
    if recorded != sha256:
        raise ValueError(
            f"the obs table {table.metadata_location} was built from another catalog "
            f"than {source}; build it again with tools/build_obs_table.py"
        )


@dataclass
class ScanStats:
    """What one scan read: the plan, the data files, and the rows."""

    plan_seconds: float = 0.0
    scan_seconds: float = 0.0
    files_total: int = 0
    files_planned: int = 0
    rows: int = 0
    # The sizes of the metadata files that the plan read: the metadata file, the
    # manifest list, and the manifests.
    metadata_bytes: int = 0
    manifests: int = 0
    # The compressed size of the projected columns of the planned files, from the
    # column sizes in the manifests. An estimate of the data that the scan reads: the
    # scan also reads each footer, and it skips the row groups that its statistics
    # exclude.
    data_bytes: int = 0
    extra: dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> dict[str, Any]:
        out = {k: v for k, v in self.__dict__.items() if k != "extra"}
        return {**out, **self.extra}


def scan(
    table: StaticTable,
    row_filter: BooleanExpression,
    selected: tuple[str, ...],
    *,
    dictionary: tuple[str, ...] = (),
) -> tuple[pa.Table, ScanStats]:
    """The rows of `table` that satisfy `row_filter`, with the `selected` columns.

    `dictionary` names the string columns to return dictionary-encoded, which is what a
    categorical column of pandas needs. pyiceberg reads the manifests and the data files
    with its thread pool; `PYICEBERG_MAX_WORKERS` sets its size.
    """
    stats = ScanStats()
    t0 = time.monotonic()
    data_scan = table.scan(row_filter=row_filter, selected_fields=selected)
    tasks = list(data_scan.plan_files())
    t1 = time.monotonic()
    arrow = ArrowScan(
        table.metadata,
        table.io,
        data_scan.projection(),
        data_scan.row_filter,
        dictionary_columns=dictionary,
    ).to_table(tasks)
    t2 = time.monotonic()

    snapshot = table.current_snapshot()
    stats.plan_seconds = round(t1 - t0, 3)
    stats.scan_seconds = round(t2 - t1, 3)
    stats.files_planned = len(tasks)
    stats.rows = int(arrow.num_rows)
    if snapshot is not None:
        summary = snapshot.summary.additional_properties if snapshot.summary else {}
        stats.files_total = int(summary.get("total-data-files", 0))
        manifests = snapshot.manifests(table.io)
        stats.manifests = len(manifests)
        stats.metadata_bytes = sum(int(m.manifest_length) for m in manifests)
    projected = {f.field_id for f in data_scan.projection().fields}
    stats.data_bytes = sum(
        int(size)
        for task in tasks
        for field_id, size in (task.file.column_sizes or {}).items()
        if field_id in projected
    )
    return arrow, stats
