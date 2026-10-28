#!/usr/bin/env python3
"""Build the prototype catalog from the Iteration 1 stores: one Parquet row for each dataset.

The catalog builder is not part of the SDK. The SDK reads a catalog; the Census 2 team
publishes one (design 2.3, 3.3). This script is the Iteration 1 stand-in for that
builder. It reads one or more prefixes that hold `<dataset_version_id>.zarr` stores.
Each store's summary record, in the attributes of its root group, gives the dataset's
identifiers and the gene reference version. The rest of the row comes from the store's
`obs` and `uns`. The builder reads no receipt.

Each row also holds the cell counts of its store, in the `cell_counts` column
(tracker 5.9). The column has one entry for each combination of the `COUNTED` values
and `is_primary_data`, with the number of cells in `n_cells`. The builder counts the
cells in the same pass that reads the values for the catalog row.

    python census2-sdk/tools/build_catalog.py s3://bucket/prefix/ catalog.parquet

`--compacted` names a prefix of compacted stores, each a second copy of several
datasets (see `census2_dataset_writer.compaction`). The builder checks each compacted
store against the per-dataset stores, and then writes the store and the rows of each
member into the member's row (`COMPACTED`). A row keeps its per-dataset `store_uri`,
which stays the primary copy. The builder refuses the catalog when a check fails.

    python census2-sdk/tools/build_catalog.py s3://bucket/stores/ catalog.parquet \\
        --compacted s3://bucket/compacted/

A remote prefix needs the `s3` extra and credentials in the environment.
"""

import argparse
import pathlib
import sys
from collections.abc import Sequence
from concurrent.futures import ThreadPoolExecutor
from typing import Any, cast

import census2_dataset_writer
import numpy as np
import pandas as pd
import pyarrow as pa  # type: ignore[import-untyped]  # pyarrow has no types
import pyarrow.parquet as pq  # type: ignore[import-untyped]
import zarr
from anndata.io import read_elem
from census2_dataset_writer import compaction

from census2_sdk.catalog import (
    COMPACTED,
    COUNTED,
    COUNTS,
    INDEXED,
    is_remote,
    open_store,
)
from census2_sdk.gene_reference import TAXON_TO_ORGANISM

# The declared type keeps a column that no store holds as a string in every catalog.
COUNT_ENTRY = pa.struct(
    [pa.field(c, pa.string()) for c in COUNTED]
    + [pa.field("is_primary_data", pa.bool_()), pa.field("n_cells", pa.int64())]
)
# Declared for the same reason: a catalog without compacted stores has the columns,
# with the same types, and every value null.
COMPACTED_TYPES = dict(
    zip(COMPACTED, (pa.string(), pa.int64(), pa.int64()), strict=True)
)


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


def summary_record(
    group: zarr.Group, name: str
) -> census2_dataset_writer.SummaryRecord:
    """The summary record of the store `name`. Refuse a store that has no usable one.

    `census2-dataset-writer` defines the record, reads it, and checks its fields.
    A catalog that dropped a store without a record quietly would drop its cells
    from every query.
    """
    try:
        record = census2_dataset_writer.read_summary(group)
    except census2_dataset_writer.SummaryRecordError as e:
        raise ValueError(f"{name}: {e}") from e
    if record is None:
        raise ValueError(f"{name}: no summary record")
    if f"{record['dataset_version_id']}.zarr" != name:
        raise ValueError(
            f"{name}: the summary record names dataset version "
            f"{record['dataset_version_id']!r}"
        )
    return record


def _row(store_uri: str) -> tuple[dict[str, Any], dict[str, Any]]:
    """The catalog row of one per-dataset store, and the identity of its cells.

    The identity is what a compacted store must match: the digest of the
    `observation_joinid` values and the stored values of `X`.
    """
    name = store_uri.rsplit("/", 1)[1]
    group = open_store(store_uri)
    record = summary_record(group, name)
    row, cells = summarize(group)
    identity: dict[str, Any] = {"obs_digest": None, "nnz": None}
    obs = group["obs"]
    assert isinstance(obs, zarr.Group)
    if "observation_joinid" in obs:
        joinids = pd.Series(read_elem(obs["observation_joinid"])).astype(str)
        identity["obs_digest"] = census2_dataset_writer.obs_digest(joinids)
    x = group["X"]
    if isinstance(x, zarr.Group) and "indptr" in x:
        indptr = cast("zarr.Array[Any]", x["indptr"])
        identity["nnz"] = int(np.asarray(indptr[-1:])[0])
    return {
        "dataset_id": record["dataset_id"],
        "dataset_version_id": record["dataset_version_id"],
        "store_uri": store_uri,
        "gene_reference_version": record["gene_reference_version"],
        **row,
        COUNTS: cells.to_dict("records"),
    }, identity


def build(
    prefixes: Sequence[str],
    out: str | pathlib.Path | None = None,
    *,
    compacted: Sequence[str] = (),
    workers: int = 8,
) -> pd.DataFrame:
    """One catalog row for each per-dataset store under `prefixes`.

    The store's summary record supplies `dataset_id`, `dataset_version_id`, and
    `gene_reference_version`. The store's `obs` and `uns` supply the rest, the
    cell counts included. The compacted stores under `compacted` supply the
    `COMPACTED` columns; see `compacted_rows`. `workers` stores are read at the
    same time.
    """
    uris = [
        f"{prefix.rstrip('/')}/{name}"
        for prefix in prefixes
        for name in _stores(prefix.rstrip("/"))
    ]
    with ThreadPoolExecutor(max_workers=workers) as pool:
        results = list(pool.map(_row, uris))
    frame = pd.DataFrame([row for row, _ in results])
    if len(frame) and frame["dataset_version_id"].duplicated().any():
        dupes = frame[frame["dataset_version_id"].duplicated(keep=False)]
        raise ValueError(
            f"a dataset version occurs twice: {sorted(dupes['store_uri'])}"
        )
    identities = {
        row["dataset_version_id"]: {**row, **identity} for row, identity in results
    }
    serving = compacted_rows(compacted, identities, workers=workers)
    for column, kind in zip(COMPACTED, ("object", "Int64", "Int64"), strict=True):
        values = [
            serving.get(v, (None, None, None))
            for v in frame.get("dataset_version_id", [])
        ]
        position = COMPACTED.index(column)
        frame[column] = pd.Series([v[position] for v in values], dtype=kind)
    if out is not None:
        table = pa.Table.from_pandas(
            frame.drop(columns=[COUNTS, *COMPACTED], errors="ignore"),
            preserve_index=False,
        )
        for column, kind in COMPACTED_TYPES.items():
            cells = frame[column].astype(object).where(frame[column].notna(), None)
            table = table.append_column(column, pa.array(list(cells), type=kind))
        counts = pa.array(list(frame.get(COUNTS, [])), type=pa.list_(COUNT_ENTRY))
        pq.write_table(table.append_column(COUNTS, counts), str(out))
    return frame


def compacted_rows(
    prefixes: Sequence[str],
    rows: dict[str, dict[str, Any]],
    *,
    workers: int = 8,
) -> dict[str, tuple[str, int, int]]:
    """The compacted store and the rows of each member, by dataset version.

    `rows` holds the catalog row and the cell identity of each per-dataset store.
    Each compacted store must agree with them; see `check_compacted`. A dataset
    version may be in one compacted store only.
    """
    uris = [
        f"{prefix.rstrip('/')}/{name}"
        for prefix in prefixes
        for name in _stores(prefix.rstrip("/"))
    ]
    with ThreadPoolExecutor(max_workers=workers) as pool:
        records = list(pool.map(lambda uri: check_compacted(uri, rows), uris))
    serving: dict[str, tuple[str, int, int]] = {}
    for uri, record in zip(uris, records, strict=True):
        for member in record["members"]:
            version = member["dataset_version_id"]
            if version in serving:
                raise ValueError(
                    f"{version} is in two compacted stores: {serving[version][0]}, {uri}"
                )
            serving[version] = (uri, member["row_start"], member["row_stop"])
    return serving


def check_compacted(
    uri: str, rows: dict[str, dict[str, Any]]
) -> census2_dataset_writer.CompactionRecord:
    """Check one compacted store against the per-dataset stores. Return its record.

    Each member must match a per-dataset store: the same `dataset_id`, number of
    cells, gene reference, `obs_digest`, and `nnz`. The members must compose: one
    organism, one gene reference, and one spatial value. The compacted store must
    hold what its record says: the member's `dataset_id` and `observation_joinid`
    values in its rows, and its stored values of `X`. These checks read two `obs`
    columns and some values of `X/indptr`, not the matrix.
    """
    name = uri.rsplit("/", 1)[1]
    group = open_store(uri)
    try:
        record = compaction.read_compaction(group)
    except census2_dataset_writer.SummaryRecordError as e:
        raise ValueError(f"{name}: {e}") from e
    if record is None:
        raise ValueError(f"{name}: no compaction record")
    members = record["members"]
    matched = []
    for member in members:
        version = member["dataset_version_id"]
        row = rows.get(version)
        if row is None:
            raise ValueError(f"{name}: member {version} has no per-dataset store")
        n_obs = member["row_stop"] - member["row_start"]
        expected = {
            "dataset_id": row["dataset_id"],
            "n_obs": row["n_obs"],
            "gene_reference_version": row["gene_reference_version"],
            "obs_digest": row["obs_digest"],
            "nnz": row["nnz"],
        }
        found = {
            "dataset_id": member["dataset_id"],
            "n_obs": n_obs,
            "gene_reference_version": record["gene_reference_version"],
            "obs_digest": member["obs_digest"],
            "nnz": member["nnz"],
        }
        differ = sorted(k for k in expected if expected[k] != found[k])
        if differ:
            raise ValueError(
                f"{name}: member {version} differs from its per-dataset store in {differ}"
            )
        matched.append(row)
    for field in ("organism", "gene_reference_version", "spatial"):
        values = {str(r[field]) for r in matched}
        if len(values) > 1:
            raise ValueError(f"{name}: the members differ in {field}: {sorted(values)}")

    # The compacted store holds what its record says.
    x = group["X"]
    obs = group["obs"]
    assert isinstance(x, zarr.Group) and isinstance(obs, zarr.Group)
    n_rows = int(cast("list[int]", x.attrs["shape"])[0])
    if n_rows != members[-1]["row_stop"]:
        raise ValueError(
            f"{name}: {n_rows} rows; the record lists {members[-1]['row_stop']}"
        )
    dataset_ids = pd.Series(read_elem(obs["dataset_id"])).astype(str).to_numpy()
    joinids = pd.Series(read_elem(obs["observation_joinid"])).astype(str)
    bounds = [m["row_start"] for m in members] + [n_rows]
    indptr = np.asarray(
        cast("zarr.Array[Any]", x["indptr"]).get_orthogonal_selection(np.array(bounds))
    )
    for k, member in enumerate(members):
        start, stop = member["row_start"], member["row_stop"]
        if set(dataset_ids[start:stop]) != {member["dataset_id"]}:
            raise ValueError(f"{name}: rows {start}:{stop} hold another dataset_id")
        if (
            census2_dataset_writer.obs_digest(joinids[start:stop])
            != member["obs_digest"]
        ):
            raise ValueError(f"{name}: rows {start}:{stop} hold other cells")
        if int(indptr[k + 1] - indptr[k]) != member["nnz"]:
            raise ValueError(f"{name}: rows {start}:{stop} hold other values of X")
    return record


def summarize(group: zarr.Group) -> tuple[dict[str, Any], pd.DataFrame]:
    """The catalog fields that one store supplies, and its cell counts."""
    obs = group["obs"]
    assert isinstance(obs, zarr.Group)
    row: dict[str, Any] = {}
    counted: dict[str, pd.Series] = {}
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
        if column in COUNTED:
            counted[column] = values
    primary: pd.Series | None = None
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
    # The writer's list of spatial assays, which is the Census 1 list.
    row["spatial"] = bool(census2_dataset_writer.assays.SPATIAL & set(assays))
    return row, cell_counts(counted, primary, n_obs)


def cell_counts(
    counted: dict[str, pd.Series], primary: pd.Series | None, n_obs: int
) -> pd.DataFrame:
    """The cell counts of one store, with the number of cells in `n_cells`.

    The table has one row for each combination of the `COUNTED` values and
    `is_primary_data`. A column that the store does not hold is null in every
    row. For a store without `is_primary_data`, every row has `is_primary_data`
    false. The rows keep a null as a separate value, so the rows of a store add
    up to its `n_obs`.
    """
    frame = pd.DataFrame(
        {
            # Categorical keys make the group step fast on stores of millions of cells.
            c: _categorical(counted[c])
            if c in counted
            else pd.Categorical([None] * n_obs)
            for c in COUNTED
        }
    )
    frame["is_primary_data"] = (
        primary.to_numpy() if primary is not None else [False] * n_obs
    )
    sizes = frame.groupby(list(frame.columns), observed=True, dropna=False).size()
    table = sizes.rename("n_cells").reset_index()
    for c in COUNTED:
        table[c] = table[c].astype(object).where(table[c].notna(), None)
    return table


def _categorical(values: pd.Series) -> pd.Categorical:
    """The values as a categorical of strings. A null stays null."""
    cat = pd.Categorical(values)
    return cat.rename_categories(str)


def main(argv: list[str]) -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "paths", nargs="+", help="prefixes of per-dataset stores, then the catalog"
    )
    p.add_argument(
        "--compacted",
        action="append",
        default=[],
        help="a prefix of compacted stores; give it once for each prefix",
    )
    p.add_argument(
        "--workers", type=int, default=8, help="stores read at the same time"
    )
    args = p.parse_args(argv)
    if len(args.paths) < 2:
        p.error("give one or more store prefixes, then the catalog")
    *prefixes, out = args.paths
    frame = build(prefixes, out, compacted=args.compacted, workers=args.workers)
    served = int(frame[COMPACTED[0]].notna().sum())
    stores = int(frame[COMPACTED[0]].nunique())
    print(f"{len(frame)} rows -> {out}; {served} read from {stores} compacted stores")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
