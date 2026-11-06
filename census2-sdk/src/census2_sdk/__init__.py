"""Query Census 2 stores through a catalog and return one AnnData.

The SDK runs on the user's machine. It reads the catalog and excludes the datasets that
cannot match. It then opens the remaining stores and concatenates the selected cells by
gene position (design 4.1 to 4.4).

    import census2_sdk

    census = census2_sdk.open_census(sources=["catalog.parquet"])
    adata = census.get_anndata("homo_sapiens", obs_filter={"tissue": "lung"})

`obs_filter` maps a catalog column to one value or to a list of values. The SDK selects
a cell when every column holds one of its values. The columns are the 18 indexed
columns of design 2.2, `dataset_id`, and `is_primary_data` with a boolean value. A
filter grammar is open (design section 8).

A source can also have an `obs` table: one Iceberg table over the `obs` fragments that
the writer puts next to the stores (`census2_sdk.obs_table`). With it, a query scans the
table once, opens only the stores that hold a selected cell, and reads only the `X`
values of those cells. The filter can then also be an Iceberg expression, a string such
as `"tissue IN ('lung') AND (cell_type = 'macrophage' OR cell_type = 'monocyte')"`; see
`census2_sdk.expressions`. `get_obs` returns the `obs` of the selected cells from the
table alone.

    census = census2_sdk.open_census(["catalog.parquet"], obs_tables=["obs/homo_sapiens"])

With `obs_engine="duckdb"`, DuckDB reads the fragments of the catalog's datasets
directly, without a table (`census2_sdk.obs_sql`). The filter can then also be a SQL
expression, such as `"tissue = 'lung' AND nnz >= 500"`. `get_obs` and `count_obs`
then open no store.

    census = census2_sdk.open_census(["catalog.parquet"], obs_engine="duckdb")

The result carries the counts in `X` and the pinned gene reference as `var`. Its `obs`
has the columns that all selected stores hold. `uns["census2"]` has the sources, the
number of catalog rows for the organism, and the number of stores opened. The cell
key `(dataset_id, observation_joinid)` is in `obs`, and `obs_names` is
`<dataset_id>_<observation_joinid>`. When no cell matches, `get_anndata` returns None.

The stores hold only the counts. `layers=["normalized"]` adds the Census 1
normalized values in `layers["normalized"]`, calculated on read.

The writer projects each store onto the reference. Thus every selected store has
the same genes in `var`. The SDK checks that every store has the same gene identifiers
in the same order, and concatenates by position. A store's `var` also marks the
dataset's assayed genes in a `measured` column, and the result leaves that column out,
because the result spans datasets.

`count_cells` counts the selected cells and opens no store. It sums the cell counts
that each catalog row holds. Its filter takes the columns of
`catalog.COUNTED` and `dataset_id`.

    n = census.count_cells("homo_sapiens", {"cell_type": "neutrophil"})

Each store's `obs` holds `dataset_id`, which the writer adds. The SDK refuses a
store whose `dataset_id` differs from the dataset of its catalog row.
"""

import os
from collections.abc import Mapping, Sequence
from concurrent.futures import ThreadPoolExecutor
from typing import TYPE_CHECKING, Any, Literal, cast

import anndata
import numpy as np
import pandas as pd
import scipy.sparse as sp
import zarr
from anndata.io import read_elem, sparse_dataset

from census2_sdk import catalog as _catalog
from census2_sdk import fragments
from census2_sdk import normalize as _normalize
from census2_sdk.catalog import CanonicalFilter, Filter, members, open_store, unrecorded
from census2_sdk.gene_reference import MEASURED_COLUMN, ORGANISMS

if TYPE_CHECKING:  # pyiceberg is the optional `iceberg` extra
    from pyiceberg.expressions import BooleanExpression
    from pyiceberg.table import StaticTable

    from census2_sdk.expressions import AnyFilter
    from census2_sdk.obs_sql import Engine, Options

__all__ = [
    "Census",
    "SourceConflict",
    "open_census",
]

PathLike = str | os.PathLike[str]


class SourceConflict(ValueError):
    """The sources or the selected stores do not compose (design 4.3, 4.4)."""


class Census:
    """`Census` holds a set of opened sources.

    Each call to `get_anndata` or `count_cells` runs one query over them.
    """

    def __init__(
        self,
        sources: Sequence[PathLike],
        obs_tables: Sequence[PathLike | None] | None = None,
        obs_engine: "str | Options | None" = None,
    ) -> None:
        self.sources = [str(s) for s in sources]
        self._counts: pd.DataFrame | None = None
        # One read per catalog, so that the rows and the counts come from one version.
        self._contents = [_catalog.read(source) for source in self.sources]
        frames = []
        for source, content in zip(self.sources, self._contents, strict=True):
            frame = _catalog.load(content)
            frame["source"] = source
            frames.append(frame)
        self.catalog = pd.concat(frames, ignore_index=True)
        dupes = self.catalog[self.catalog.duplicated("dataset_id", keep=False)]
        if len(dupes):
            listing = dupes.groupby("dataset_id")["source"].agg(sorted).to_dict()
            raise SourceConflict(
                f"a dataset_id occurs in more than one source: {listing}"
            )
        self._tables: dict[str, StaticTable] = {}
        # The scan statistics of the last query through a table, also when it
        # selected no cell and returned None.
        self.last_scan: list[dict[str, Any]] = []
        if obs_tables is not None:
            self._open_tables([None if t is None else str(t) for t in obs_tables])
        self._sql: Engine | None = None
        if obs_engine is not None:
            self._open_engine(obs_engine)

    def _open_engine(self, obs_engine: "str | Options") -> None:
        """Open DuckDB over the fragments of every catalog row (`census2_sdk.obs_sql`)."""
        from census2_sdk import obs_sql

        if self._tables:
            raise ValueError("a census takes obs_tables or obs_engine, not both")
        if isinstance(obs_engine, str):
            if obs_engine != "duckdb":
                raise ValueError(f"no obs engine {obs_engine!r}; the one is 'duckdb'")
            obs_engine = obs_sql.Options()
        self._sql = obs_sql.Engine(self.catalog, obs_engine)

    @property
    def sql_engine(self) -> "Engine":
        """The DuckDB engine of a census opened with `obs_engine`."""
        if self._sql is None:
            raise ValueError("this census has no obs engine")
        return self._sql

    def _open_tables(self, uris: list[str | None]) -> None:
        """Open the `obs` table of each source that has one, and check its catalog."""
        from census2_sdk import obs_table

        if len(uris) != len(self.sources):
            raise ValueError("obs_tables needs one entry for each source")
        for source, content, uri in zip(
            self.sources, self._contents, uris, strict=True
        ):
            if uri is None:
                continue
            table = obs_table.open_table(uri)
            obs_table.check_catalog(table, obs_table.catalog_sha256(content), source)
            self._tables[source] = table

    def get_anndata(
        self,
        organism: str,
        obs_filter: "AnyFilter | None" = None,
        *,
        layers: Sequence[str] = (),
        max_workers: int = 1,
    ) -> anndata.AnnData | None:
        """Return the selected cells of `organism` as one AnnData.

        Return None when no cell matches, whether or not the catalog excluded
        every store.

        `layers` names the layers to add next to the counts. The one layer is
        `normalized`, which the SDK calculates from the counts with the Census 1
        algorithm, without its rounding (`census2_sdk.normalize`).

        `max_workers` is the largest number of stores that the SDK reads at the
        same time, and it must be a positive integer. The result is the same for
        every allowed value: the stores keep their catalog order.

        When the census has an `obs` table, the query goes through it; see
        `_get_anndata_from_tables`. The filter can then be an Iceberg expression.
        """
        unknown = sorted(set(layers) - {"normalized"})
        if unknown:
            raise ValueError(f"no layer {unknown}; the one layer is 'normalized'")
        if self._sql is not None:
            from census2_sdk import obs_sql

            return obs_sql.get_anndata(self, organism, obs_filter, layers, max_workers)
        if self._tables:
            return self._get_anndata_from_tables(
                organism, obs_filter, layers, max_workers
            )
        if obs_filter is not None and not isinstance(obs_filter, Mapping):
            raise ValueError("an expression filter needs an obs table")
        canonical = _catalog.canonicalize(obs_filter)
        rows = _catalog.exclude(self.catalog, organism, canonical)
        _check_composition(rows, ORGANISMS[organism]["version"])

        records = cast("list[dict[str, Any]]", rows.to_dict("records"))
        # ponytail: serial by default; raise the default after the store-count benchmark.
        with ThreadPoolExecutor(max_workers=max_workers) as pool:
            parts = list(
                pool.map(lambda row: _read_store(row, canonical, layers), records)
            )
        uris = list(rows["store_uri"].astype(str))
        _check_var(parts, uris)
        if not sum(part.n_obs for part in parts):
            return None
        # Zero-row parts stay in, so the obs schema is the columns that all
        # opened stores hold.
        result = anndata.concat(parts, join="inner")
        if not result.obs_names.is_unique:
            raise SourceConflict("the cell key is not unique across the result")
        # Each part carries only the var index. The stores record one reference
        # version and one index, so the first store supplies the var columns,
        # except the assayed genes of its own dataset.
        var = cast("pd.DataFrame", read_elem(open_store(uris[0])["var"]))
        result.var = var.drop(columns=MEASURED_COLUMN, errors="ignore")
        result.uns["census2"] = {
            "sources": self.sources,
            "stores_in_catalog": int((self.catalog["organism"] == organism).sum()),
            "stores_opened": len(rows),
        }
        return result

    def _candidates(self, organism: str, expr: "BooleanExpression") -> pd.DataFrame:
        """The catalog rows of `organism` that can hold a cell that satisfies `expr`."""
        from census2_sdk import expressions

        scoped = self.catalog[self.catalog["organism"] == organism]
        records = cast("list[dict[str, Any]]", scoped.to_dict("records"))
        keep = [expressions.can_match(expr, row) for row in records]
        return scoped[np.array(keep, dtype=bool)]

    def _scan(
        self,
        rows: pd.DataFrame,
        expr: "BooleanExpression",
        columns: Sequence[str] | None,
    ) -> tuple[dict[str, Any], list[dict[str, Any]]]:
        """Scan the `obs` table of each source for the cells of `rows` that satisfy `expr`.

        Return the selected rows of each dataset, as an Arrow table in store row order,
        and the statistics of each scan. `columns` names the `obs` columns to return;
        None returns all. The cell key and the row columns are always in the result.
        """
        import pyarrow as pa  # type: ignore[import-untyped]
        import pyarrow.compute as pc  # type: ignore[import-untyped]

        from census2_sdk import expressions, obs_table

        cells: dict[str, Any] = {}
        stats: list[dict[str, Any]] = []
        for source, table in self._tables.items():
            mine = rows[rows["source"] == source]
            if not len(mine):
                continue
            records = cast("list[dict[str, Any]]", mine.to_dict("records"))
            ids = {str(d) for d in mine["dataset_id"]}
            row_filter = expressions.with_members(expr, records)
            names = [f.name for f in table.schema().fields]
            obs_names = [n for n in names if n not in fragments.ROW_COLUMNS]
            wanted = obs_names if columns is None else list(columns)
            key = ["dataset_id", "observation_joinid"]
            selected = tuple(dict.fromkeys([*key, *wanted, *fragments.ROW_COLUMNS]))
            # Read the label columns dictionary-encoded, for pandas categoricals. Not
            # the cell key, whose values are unique, and not a column of a prefix
            # term: pyarrow has no starts_with kernel for a dictionary column.
            prefixed = expressions.prefix_columns(row_filter)
            strings = tuple(
                f.name
                for f in table.schema().fields
                if f.name in selected
                and str(f.field_type) == "string"
                and f.name not in {"observation_joinid", *prefixed}
            )
            if SCAN_KEY_FIRST and expressions.columns(expr) - {"dataset_id"}:
                # A first scan of the key alone finds the datasets with a match, so
                # that the second scan reads every column of those fragments only.
                keys, first = obs_table.scan(
                    table, row_filter, ("dataset_id",), datasets=ids
                )
                stats.append({"source": source, "step": "key", **first.as_dict()})
                ids = set(keys.column("dataset_id").to_pylist())
                if not ids:
                    continue
            arrow, scan_stats = obs_table.scan(
                table, row_filter, selected, dictionary=strings, datasets=ids
            )
            stats.append({"source": source, **scan_stats.as_dict()})
            if not arrow.num_rows:
                continue
            # Group the rows by dataset, in store row order. Arrow sorts no
            # dictionary column, so sort the integer codes with numpy, and only
            # when the scan did not already return each file's rows together.
            ids = arrow.column("dataset_id").combine_chunks()
            if not pa.types.is_dictionary(ids.type):
                ids = pc.dictionary_encode(ids)
            codes = ids.indices.to_numpy(zero_copy_only=False)
            rows_of = arrow.column("store_row").to_numpy()
            change = np.flatnonzero(codes[1:] != codes[:-1]) + 1
            grouped = len(change) + 1 == len(np.unique(codes)) and not bool(
                np.any((np.diff(rows_of) <= 0) & (codes[1:] == codes[:-1]))
            )
            if not grouped:
                order = np.lexsort((rows_of, codes))
                arrow, codes = arrow.take(order), codes[order]
                change = np.flatnonzero(codes[1:] != codes[:-1]) + 1
            names = ids.dictionary.to_pylist()
            starts = [0, *change.tolist()]
            stops = [*change.tolist(), len(codes)]
            for start, stop in zip(starts, stops, strict=True):
                cells[str(names[codes[start]])] = arrow.slice(start, stop - start)
        return cells, stats

    def _get_anndata_from_tables(
        self,
        organism: str,
        obs_filter: "AnyFilter | None",
        layers: Sequence[str],
        max_workers: int,
    ) -> anndata.AnnData | None:
        """`get_anndata` through the `obs` tables.

        The catalog excludes the datasets that cannot match (`expressions.can_match`).
        One scan of each table then selects the cells of the remaining datasets. The
        SDK opens only the stores that hold a selected cell, in catalog order, and reads
        the `X` values of the selected cells. The result `obs` comes from the scan. A
        source without a table goes through its stores, as without tables, and takes
        a mapping filter only.

        The result equals the result through the stores, value for value. A categorical
        column of the result holds the labels that the result uses, where the result
        through the stores holds the labels of every store that it opened.
        """
        from census2_sdk import expressions

        expr = expressions.from_filter(obs_filter)
        rows = self._candidates(organism, expr)
        _check_composition(rows, ORGANISMS[organism]["version"])
        plain = rows[~rows["source"].isin(list(self._tables))]
        canonical: CanonicalFilter = {}
        if len(plain):
            if obs_filter is not None and not isinstance(obs_filter, Mapping):
                raise ValueError(
                    "an expression filter needs an obs table for every source"
                )
            canonical = _catalog.canonicalize(obs_filter)
        cells, stats = self._scan(rows, expr, None)
        self.last_scan = stats

        def part(row: dict[str, Any]) -> anndata.AnnData | None:
            if row["source"] not in self._tables:
                return _read_store(row, canonical, layers)
            selected = cells.get(str(row["dataset_id"]))
            return None if selected is None else _read_cells(row, selected, layers)

        records = cast("list[dict[str, Any]]", rows.to_dict("records"))
        with ThreadPoolExecutor(max_workers=max_workers) as pool:
            found = list(pool.map(part, records))
        opened = [(p, r) for p, r in zip(found, records, strict=True) if p is not None]
        parts = [p for p, _ in opened]
        uris = [str(r["store_uri"]) for _, r in opened]
        if not sum(part.n_obs for part in parts):
            return None
        # The first opened store supplies var. A part from the table holds gene
        # positions only, so give it the reference genes when the number agrees.
        var = cast("pd.DataFrame", read_elem(open_store(uris[0])["var"]))
        for (p, r), uri in zip(opened, uris, strict=True):
            if r["source"] not in self._tables:
                continue
            if p.n_vars != len(var):
                raise SourceConflict(
                    f"{uri}: X has {p.n_vars} genes, and {uris[0]} has {len(var)}; "
                    "the stores are not on one gene reference"
                )
            p.var_names = cast("Sequence[str]", var.index)
        _check_var(parts, uris)
        result = anndata.concat(parts, join="inner")
        if not result.obs_names.is_unique:
            raise SourceConflict("the cell key is not unique across the result")
        result.var = var.drop(columns=MEASURED_COLUMN, errors="ignore")
        result.uns["census2"] = {
            "sources": self.sources,
            "stores_in_catalog": int((self.catalog["organism"] == organism).sum()),
            "catalog_candidates": len(rows),
            "stores_opened": len(parts),
            "obs_table": stats,
        }
        return result

    def get_obs(
        self,
        organism: str,
        obs_filter: "AnyFilter | None" = None,
        *,
        columns: Sequence[str] | None = None,
        dtype_backend: Literal["numpy", "pyarrow"] = "numpy",
    ) -> pd.DataFrame:
        """The `obs` of the selected cells of `organism`, from the `obs` tables alone.

        The SDK opens no store. `columns` names the `obs` columns to return; None
        returns all. The frame always holds the cell key and the row columns of
        `fragments.ROW_COLUMNS`, and the columns keep the order of the table schema.
        Its index is `<dataset_id>_<observation_joinid>`. Every source needs an `obs`
        table, or the census needs an `obs_engine`.

        With `dtype_backend="numpy"`, a label column is a categorical of the labels
        that the result uses, sorted, as in `get_anndata`. With `"pyarrow"`, each
        column is a `pd.ArrowDtype` column over the scan's Arrow data, which pandas
        does not copy: a label column is then an Arrow dictionary column.

        `frame.attrs["census2"]` holds the statistics of the scan, and the seconds
        of the scan (`table_seconds`) and of the frame (`build_seconds`).
        """
        import time

        import pyarrow as pa
        import pyarrow.compute as pc

        from census2_sdk import expressions

        if self._sql is not None:
            from census2_sdk import obs_sql

            return obs_sql.get_obs(self, organism, obs_filter, columns, dtype_backend)
        if set(self.sources) - set(self._tables):
            raise ValueError("get_obs needs an obs table for every source")
        if dtype_backend not in ("numpy", "pyarrow"):
            raise ValueError(
                f"dtype_backend is 'numpy' or 'pyarrow', not {dtype_backend!r}"
            )
        t0 = time.monotonic()
        expr = expressions.from_filter(obs_filter)
        rows = self._candidates(organism, expr)
        cells, stats = self._scan(rows, expr, columns)
        order = [cells[str(d)] for d in rows["dataset_id"] if str(d) in cells]
        t1 = time.monotonic()
        # One Arrow table first, so that each categorical has one set of labels.
        frame = pd.DataFrame()
        if order:
            arrow = pa.concat_tables(order).unify_dictionaries()
            if dtype_backend == "pyarrow":
                frame = cast(
                    "pd.DataFrame", arrow.to_pandas(types_mapper=pd.ArrowDtype)
                )
            else:
                frame = _obs_frame(arrow, keep_rows=True)
            # The cell key, joined in Arrow rather than one Python string at a time.
            key = pc.binary_join_element_wise(
                pc.cast(arrow.column("dataset_id"), pa.string()),
                arrow.column("observation_joinid"),
                "_",
            )
            frame.index = pd.Index(key.to_pandas())
        frame.attrs["census2"] = {
            "catalog_candidates": len(rows),
            "obs_table": stats,
            "table_seconds": round(t1 - t0, 3),
            "build_seconds": round(time.monotonic() - t1, 3),
        }
        return frame

    def count_obs(
        self,
        organism: str,
        obs_filter: "AnyFilter | None" = None,
        *,
        by: Sequence[str] = (),
    ) -> pd.DataFrame:
        """The number of selected cells of `organism`, from the `obs` tables alone.

        The count takes any filter and any `obs` column in `by`, for example
        `donor_id`, which the catalog does not count. The frame has the `by` columns
        and `n_cells`, one row for each group, in no fixed order. Every source needs
        an `obs` table, or the census needs an `obs_engine`.
        """
        if self._sql is not None:
            from census2_sdk import obs_sql

            return obs_sql.count_obs(self, organism, obs_filter, by)

        import pyarrow as pa

        from census2_sdk import expressions, obs_table

        if set(self.sources) - set(self._tables):
            raise ValueError("count_obs needs an obs table for every source")
        expr = expressions.from_filter(obs_filter)
        rows = self._candidates(organism, expr)
        parts, stats = [], []
        for source, table in self._tables.items():
            mine = rows[rows["source"] == source]
            if not len(mine):
                continue
            records = cast("list[dict[str, Any]]", mine.to_dict("records"))
            arrow, scan_stats = obs_table.scan(
                table,
                expressions.with_members(expr, records),
                tuple(dict.fromkeys(["dataset_id", *by])),
                datasets={str(d) for d in mine["dataset_id"]},
            )
            stats.append({"source": source, **scan_stats.as_dict()})
            parts.append(arrow.group_by(list(by)).aggregate([([], "count_all")]))
        if parts:
            counts = (
                pa.concat_tables(parts)
                .group_by(list(by))
                .aggregate([("count_all", "sum")])
            )
            frame = cast("pd.DataFrame", counts.to_pandas())
            frame = frame.rename(columns={"count_all_sum": "n_cells"})
        else:
            frame = pd.DataFrame({**{c: [] for c in by}, "n_cells": []})
        frame.attrs["census2"] = {"catalog_candidates": len(rows), "obs_table": stats}
        return frame[[*by, "n_cells"]]

    def count_cells(
        self,
        organism: str,
        obs_filter: Filter | None = None,
        *,
        primary_only: bool = False,
    ) -> int:
        """Return the number of cells of `organism` that match `obs_filter`.

        The SDK sums the cell counts in the catalog rows and opens no store. A
        cell matches by the rule of `get_anndata`, and the catalog excludes the
        same datasets. Thus, without `primary_only`, the count equals the `n_obs`
        of the result of `get_anndata` when there is one. The count does not
        apply the checks of design 4.3, because a count joins no matrices.

        The filter takes the columns of `catalog.COUNTED` and `dataset_id`, and a
        filter on another column raises ValueError. With `primary_only`, the count
        includes only the cells whose `is_primary_data` is true. A catalog without
        the `cell_counts` column raises ValueError.
        """
        canonical = _catalog.canonicalize(obs_filter)
        uncounted = sorted(
            set(canonical) - {"dataset_id", _catalog.PRIMARY, *_catalog.COUNTED}
        )
        if uncounted:
            raise ValueError(f"the catalog counts no {uncounted}; use get_anndata")
        counts = self._cell_counts()
        candidates = _catalog.exclude(self.catalog, organism, canonical)
        rows = counts[
            counts["dataset_version_id"].isin(candidates["dataset_version_id"])
        ]
        if primary_only:
            rows = rows[rows["is_primary_data"].astype(bool)]
        keep = np.ones(len(rows), dtype=bool)
        for column, wanted in canonical.items():
            if column == "dataset_id":
                continue  # the catalog answered it
            # Test each distinct value once; a null's code -1 indexes the last False.
            codes, values = pd.factorize(rows[column])
            by_value = [_hit(v, column, wanted) for v in values]
            keep &= np.array(by_value + [False], dtype=bool)[codes]
        return int(rows["n_cells"].to_numpy()[keep].sum())

    def _cell_counts(self) -> pd.DataFrame:
        """The cell count entries of all sources, from the catalog content."""
        if self._counts is None:
            frames = []
            for source, content in zip(self.sources, self._contents, strict=True):
                try:
                    frames.append(_catalog.load_counts(content))
                except ValueError as e:
                    raise ValueError(f"{source}: {e}") from e
            self._counts = pd.concat(frames, ignore_index=True)
        return self._counts


def open_census(
    sources: Sequence[PathLike],
    obs_tables: Sequence[PathLike | None] | None = None,
    obs_engine: "str | Options | None" = None,
) -> Census:
    """Open catalogs, and the `obs` table of each source that has one.

    `obs_tables` has one entry for each source: a table location or metadata file, or
    None for a source without a table. See `census2_sdk.obs_table`.

    `obs_engine` is `"duckdb"` or a `census2_sdk.obs_sql.Options`: DuckDB then reads
    the `obs` fragments of the catalog rows directly, and the filter can be SQL. See
    `census2_sdk.obs_sql`. A census takes `obs_tables` or `obs_engine`, not both.
    """
    return Census(sources, obs_tables, obs_engine)


def _check_composition(rows: pd.DataFrame, pinned: str) -> None:
    """Design 4.3: the checks before concatenation. Name the conflicting values.

    Every row must record the `gene_reference_version` that this SDK pins for the
    organism. The SDK does not remap a store onto another reference. A row that
    records none has no summary record that names the reference of its store. A
    catalog without the column records none.

    Spatial and non-spatial stores do not compose. A row that records no
    `spatial` value cannot prove that its store composes with another, so a
    result of two or more stores refuses it. One store always composes with
    itself. A catalog without the column records no value.
    """
    recorded = _recorded(rows, "gene_reference_version")
    missing = sorted(
        str(d) for d, v in zip(rows["dataset_id"], recorded, strict=True) if v is None
    )
    if missing:
        raise SourceConflict(f"no gene_reference_version recorded for {missing}")
    versions = sorted({str(v) for v in recorded})
    if versions and versions != [pinned]:
        raise SourceConflict(
            f"gene_reference_version {versions} differs from the pinned {pinned!r}"
        )
    if len(rows) < 2:
        return
    spatial = _recorded(rows, "spatial")
    unknown = [d for d, v in zip(rows["dataset_id"], spatial, strict=True) if v is None]
    if unknown:
        raise SourceConflict(f"no spatial value is recorded for {sorted(unknown)}")
    if len({bool(v) for v in spatial}) > 1:
        raise SourceConflict("spatial and non-spatial stores in one result")


def _recorded(rows: pd.DataFrame, column: str) -> list[Any]:
    """One entry per row: the recorded value, or None when the row records none.

    A catalog without the column records none. `unrecorded` gives the rule.
    """
    if column not in rows.columns:
        return [None] * len(rows)
    return [None if unrecorded(v) else v for v in rows[column]]


def _mask(obs: pd.DataFrame, obs_filter: CanonicalFilter) -> np.ndarray:
    """The cells of one store that match `obs_filter`, by the same rule as the catalog."""
    mask = np.ones(len(obs), dtype=bool)
    for column, wanted in obs_filter.items():
        if column == "dataset_id":
            continue  # the catalog answered it; the store holds the dataset whole
        if column not in obs.columns:
            return np.zeros(len(obs), dtype=bool)  # the catalog kept an unproven store
        col = obs[column]
        if column == _catalog.PRIMARY:
            hit = np.isin(col.to_numpy().astype(bool), [w == "True" for w in wanted])
        elif isinstance(col.dtype, pd.CategoricalDtype):
            # Apply the rule once for each category, then index by the codes. A
            # null has code -1, which indexes the False at the end.
            by_category = [_hit(c, column, wanted) for c in col.cat.categories]
            hit = np.array(by_category + [False], dtype=bool)[col.cat.codes.to_numpy()]
        else:
            hit = np.fromiter((_hit(v, column, wanted) for v in col), bool, len(col))
        mask &= hit
    return mask


def _hit(value: Any, column: str, wanted: frozenset[str]) -> bool:
    """True when one cell's value holds a wanted member. A null holds none."""
    return not pd.isna(value) and bool(wanted & set(members(str(value), column)))


def _read_store(
    row: dict[str, Any], obs_filter: CanonicalFilter, layers: Sequence[str] = ()
) -> anndata.AnnData:
    """The selected cells of one store. `X` holds the counts on the reference genes.

    `var` holds only the index, the `feature_id`s, which is all that `_check_var`
    compares. `get_anndata` reads the var columns from one store. `layers` names
    the layers to calculate from the counts.
    """
    group = open_store(row["store_uri"])
    obs_group = group["obs"]
    if not isinstance(obs_group, zarr.Group):
        raise TypeError(f"{row['store_uri']}: obs must be an encoded dataframe")
    # To read obs faster, the SDK decodes the two per-cell string columns only
    # where the result needs them. It does not read the index (the source cell
    # names), because the result replaces the index with the cell key. It reads
    # the cell key after the filter, for the kept cells only.
    columns = cast("list[str]", obs_group.attrs["column-order"])
    key = "observation_joinid"
    obs = pd.DataFrame(
        {name: read_elem(obs_group[name]) for name in columns if name != key}
    )
    dataset_id = str(row["dataset_id"])
    _check_dataset_id(obs, dataset_id, str(row["store_uri"]))
    keep = np.flatnonzero(_mask(obs, obs_filter))
    obs = obs.iloc[keep].copy()
    if key in columns:
        obs.insert(columns.index(key), key, _read_column_rows(obs_group[key], keep))
    _check_cell_key(obs, str(row["store_uri"]))
    var_group = group["var"]
    if not isinstance(var_group, zarr.Group):
        raise TypeError(f"{row['store_uri']}: var must be an encoded dataframe")
    index = read_elem(var_group[str(var_group.attrs["_index"])])
    var = pd.DataFrame(index=pd.Index(np.asarray(index).astype(str)))
    x = _rows(group["X"], keep)
    computed = {}
    if "normalized" in layers:
        assay = obs["assay_ontology_term_id"].astype(str)
        full_gene = assay.isin(_normalize.FULL_GENE_ASSAYS).to_numpy()
        length = None
        if full_gene.any():
            length = np.asarray(read_elem(var_group["feature_length"]), dtype=np.int64)
        computed["normalized"] = _normalize.normalized(x, full_gene, length)

    obs.index = pd.Index(dataset_id + "_" + obs["observation_joinid"].astype(str))
    return anndata.AnnData(X=x, obs=obs, var=var, layers=computed)


# Read the X values of the selected cells from the ranges in the obs table. False reads
# the whole X/indptr of each store first, as a read without the table does. The
# benchmark sets it to compare the two.
X_FROM_RANGES = True

# Scan the key columns first, then every column of the fragments that hold a match.
# The second scan then reads no fragment of a dataset that the catalog kept and that
# holds no match. The benchmark sets it.
SCAN_KEY_FIRST = False


def _read_cells(
    row: dict[str, Any], cells: Any, layers: Sequence[str] = ()
) -> anndata.AnnData:
    """The selected cells of one store, from their rows in the `obs` table.

    `cells` is an Arrow table of the store's selected rows, in store row order. The
    result `obs` comes from it, so the SDK does not read the store's `obs`. `X` holds
    the values of the selected rows, from their ranges in `X/data` and `X/indices`.
    `layers` names the layers to calculate from the counts.

    The SDK does not read the store's `var` index. The catalog row records the
    pinned `gene_reference_version`, and the writer projects every store onto that
    reference, so `var` holds positions only: `_get_anndata_from_tables` puts the
    reference genes in after it checks the number of genes.
    """
    store_uri = str(row["store_uri"])
    group = open_store(store_uri)
    x = group["X"]
    if not isinstance(x, zarr.Group):
        raise TypeError(f"{store_uri}: X must be an encoded sparse matrix")
    shape = tuple(int(n) for n in cast("list[int]", x.attrs["shape"]))
    n_obs = row.get("n_obs")
    if not unrecorded(n_obs) and shape[0] != int(cast("int", n_obs)):
        raise ValueError(
            f"{store_uri}: X has {shape[0]} rows; the catalog gives {n_obs}"
        )
    obs = _obs_frame(cells)
    _check_cell_key(obs, store_uri)
    if X_FROM_RANGES:
        starts = cells.column(fragments.X_START).to_numpy()
        ends = cells.column(fragments.X_END).to_numpy()
        matrix = _rows_by_range(x, starts, ends, shape[1])
    else:
        matrix = _rows(x, cells.column(fragments.STORE_ROW).to_numpy())
    computed = {}
    if "normalized" in layers:
        assay = obs["assay_ontology_term_id"].astype(str)
        full_gene = assay.isin(_normalize.FULL_GENE_ASSAYS).to_numpy()
        length = None
        if full_gene.any():
            var_group = group["var"]
            if not isinstance(var_group, zarr.Group):
                raise TypeError(f"{store_uri}: var must be an encoded dataframe")
            length = np.asarray(read_elem(var_group["feature_length"]), dtype=np.int64)
        computed["normalized"] = _normalize.normalized(matrix, full_gene, length)
    obs.index = pd.Index(
        str(row["dataset_id"]) + "_" + obs["observation_joinid"].astype(str)
    )
    var = pd.DataFrame(index=pd.RangeIndex(shape[1]).astype(str))
    return anndata.AnnData(X=matrix, obs=obs, var=var, layers=computed)


def _obs_frame(cells: Any, *, keep_rows: bool = False) -> pd.DataFrame:
    """The `obs` columns of an Arrow table of selected rows, as the store holds them.

    A dictionary-encoded column becomes a categorical whose categories are the labels
    that the rows use, sorted. The row columns are in the frame with `keep_rows` only.
    """
    import pyarrow as pa

    names = [
        n for n in cells.column_names if keep_rows or n not in fragments.ROW_COLUMNS
    ]
    columns: dict[str, Any] = {}
    for name in names:
        column = cells.column(name)
        if pa.types.is_dictionary(column.type):
            columns[name] = _categorical(column.combine_chunks())
        else:
            columns[name] = column.to_pandas()
    return pd.DataFrame(columns)


def _categorical(array: Any) -> pd.Categorical:
    """A dictionary-encoded Arrow array as a categorical of the labels that it uses, sorted.

    The work is on the integer codes, with one Python string for each label and not
    for each cell. A null becomes the code -1.
    """
    labels = [str(v) for v in array.dictionary.to_pylist()]
    indices = array.indices.fill_null(0).to_numpy(zero_copy_only=False)
    used = np.flatnonzero(np.bincount(indices, minlength=len(labels)) > 0)
    order = sorted(used.tolist(), key=labels.__getitem__)
    remap = np.full(len(labels), -1, dtype=np.int32)
    remap[order] = np.arange(len(order), dtype=np.int32)
    codes = remap[indices]
    if array.null_count:
        codes[array.is_null().to_numpy(zero_copy_only=False)] = -1
    categories = pd.Index([labels[i] for i in order])
    # The codes are valid by construction, so pandas need not check each one.
    return cast(
        "pd.Categorical",
        pd.Categorical.from_codes(cast("Any", codes), categories, validate=False),
    )


def _rows_by_range(
    node: zarr.Group, starts: np.ndarray, ends: np.ndarray, n_vars: int
) -> sp.csr_matrix:
    """The rows of a CSR matrix whose values are at `[starts[i], ends[i])`.

    The ranges come from the `obs` table, so the read does not touch `X/indptr`. One
    selection of the value positions reads the touched chunks of `data` and `indices`
    with zarr's concurrency.
    """
    starts = np.asarray(starts, dtype=np.int64)
    ends = np.asarray(ends, dtype=np.int64)
    data = node["data"]
    indices = node["indices"]
    assert isinstance(data, zarr.Array) and isinstance(indices, zarr.Array)
    lengths = ends - starts
    indptr = np.concatenate([[0], np.cumsum(lengths)]).astype(np.int64)
    if not len(starts):
        return sp.csr_matrix((0, n_vars), dtype=data.dtype)
    if bool((starts[1:] == ends[:-1]).all()):
        # One contiguous run: one slice, and no position array.
        selection: Any = slice(int(starts[0]), int(ends[-1]))
    else:
        selection = np.repeat(starts - indptr[:-1], lengths) + np.arange(
            int(indptr[-1]), dtype=np.int64
        )
    values = np.asarray(data[selection])
    columns = np.asarray(indices[selection])
    return sp.csr_matrix((values, columns, indptr), shape=(len(starts), n_vars))


def _check_dataset_id(obs: pd.DataFrame, dataset_id: str, store_uri: str) -> None:
    """The store's `dataset_id` column names the dataset of its catalog row."""
    if "dataset_id" not in obs.columns:
        raise ValueError(f"{store_uri}: obs has no dataset_id column")
    # unique() works on the categorical codes, not on one string per cell.
    found = {str(v) for v in obs["dataset_id"].unique()}
    if found != {dataset_id}:
        raise ValueError(
            f"{store_uri}: obs dataset_id holds {sorted(found)}, and the catalog row "
            f"names {dataset_id}"
        )


def _read_column_rows(node: zarr.Array[Any] | zarr.Group, rows: np.ndarray) -> Any:
    """The values of one `obs` column at the row positions `rows`.

    For a nullable string column, the SDK makes strings for these rows only. The
    dtype is the one that `anndata` reads for the same column.
    """
    if node.attrs.get("encoding-type") != "nullable-string-array":
        return pd.array(cast("Any", read_elem(node)))[rows]
    if not isinstance(node, zarr.Group):
        raise TypeError(f"{node.store_path} must be an encoded group")
    na_value = np.nan if node.attrs.get("na-value", "NA") == "NaN" else pd.NA
    dtype = pd.StringDtype(na_value=na_value)
    if not rows.size:
        return pd.array([], dtype=dtype)
    values, mask = node["values"], node["mask"]
    if not (isinstance(values, zarr.Array) and isinstance(mask, zarr.Array)):
        raise TypeError(f"{node.store_path} must hold a values and a mask array")
    result: Any = pd.array(np.asarray(values.oindex[rows]), dtype=dtype)
    result[np.asarray(mask.oindex[rows])] = pd.NA
    return result


def _check_cell_key(obs: pd.DataFrame, store_uri: str) -> None:
    """Design 1.6: the kept cells' `observation_joinid` is present, non-null, and unique.

    The writer checks every cell of the store.
    """
    if "observation_joinid" not in obs.columns:
        raise ValueError(f"{store_uri}: obs has no observation_joinid column")
    joinid = obs["observation_joinid"]
    if bool(joinid.isna().any()) or not joinid.astype(str).is_unique:
        raise ValueError(f"{store_uri}: observation_joinid is null or not unique")


def _check_var(parts: list[anndata.AnnData], store_uris: list[str]) -> None:
    """Concatenation by position needs one `var`: every store's genes, in one order."""
    if not parts:
        return
    first = parts[0].var_names
    for part, store_uri in zip(parts[1:], store_uris[1:], strict=True):
        if not part.var_names.equals(first):
            raise SourceConflict(
                f"{store_uri}: var differs from {store_uris[0]}; the stores are not "
                "on one gene reference"
            )


def _rows(node: zarr.Array[Any] | zarr.Group, keep: np.ndarray) -> sp.csr_matrix:
    """The selected rows of a CSR matrix, read without loading the rest."""
    if not isinstance(node, zarr.Group):
        raise TypeError(f"{node.store_path} must be an encoded sparse matrix")
    dataset = sparse_dataset(node)
    if not len(keep):
        return sp.csr_matrix((0, dataset.shape[1]), dtype=dataset.dtype)
    return sp.csr_matrix(cast("sp.csr_matrix", dataset[keep]))
