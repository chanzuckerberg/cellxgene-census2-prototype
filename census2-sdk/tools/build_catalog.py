#!/usr/bin/env python3
"""Build the prototype catalog from the stores: one Parquet row for each store.

The catalog builder is not part of the SDK. The SDK reads a catalog; the Census 2 team
publishes one (design 2.3, 3.3). This script is the prototype stand-in for that
builder. It reads a prefix that holds `<dataset_version_id>.zarr` stores. Each store's
summary record, in the attributes of its root group, gives the dataset's identifiers and
the gene reference version. The rest of the row comes from the store's `obs` and `uns`.
The builder reads no receipt.

Each row also holds the cell counts of its store, in the `cell_counts` column. The
column has one entry for each combination of the `COUNTED` values and
`is_primary_data`, with the number of cells in `n_cells`. The builder counts the cells
in the same pass that reads the values for the catalog row.

    python census2-sdk/tools/build_catalog.py s3://bucket/prefix/ catalog.parquet

A remote prefix needs the `s3` extra and credentials in the environment.
"""

import pathlib
import sys
from typing import Any

import census2_dataset_writer
import pandas as pd
import pyarrow as pa  # type: ignore[import-untyped]  # pyarrow has no types
import pyarrow.parquet as pq  # type: ignore[import-untyped]
import zarr
from anndata.io import read_elem

from census2_sdk.catalog import COUNTED, COUNTS, INDEXED, is_remote, open_store
from census2_sdk.gene_reference import TAXON_TO_ORGANISM

# The declared type keeps a column that no store holds as a string in every catalog.
COUNT_ENTRY = pa.struct(
    [pa.field(c, pa.string()) for c in COUNTED]
    + [pa.field("is_primary_data", pa.bool_()), pa.field("n_cells", pa.int64())]
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


def build(prefix: str, out: str | pathlib.Path | None = None) -> pd.DataFrame:
    """One catalog row for each store under `prefix`.

    The store's summary record supplies `dataset_id`, `dataset_version_id`, and
    `gene_reference_version`. The store's `obs` and `uns` supply the rest, the
    cell counts included.
    """
    prefix = prefix.rstrip("/")
    rows: list[dict[str, Any]] = []
    for name in _stores(prefix):
        store_uri = f"{prefix}/{name}"
        group = open_store(store_uri)
        record = summary_record(group, name)
        row, cells = summarize(group)
        rows.append(
            {
                "dataset_id": record["dataset_id"],
                "dataset_version_id": record["dataset_version_id"],
                "store_uri": store_uri,
                "gene_reference_version": record["gene_reference_version"],
                **row,
                COUNTS: cells.to_dict("records"),
            }
        )
    frame = pd.DataFrame(rows)
    if out is not None:
        table = pa.Table.from_pandas(
            frame.drop(columns=COUNTS, errors="ignore"), preserve_index=False
        )
        counts = pa.array(list(frame.get(COUNTS, [])), type=pa.list_(COUNT_ENTRY))
        pq.write_table(table.append_column(COUNTS, counts), str(out))
    return frame


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
    if len(argv) != 2:
        print(__doc__, file=sys.stderr)
        return 2
    frame = build(argv[0], argv[1])
    print(f"{len(frame)} rows -> {argv[1]}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
