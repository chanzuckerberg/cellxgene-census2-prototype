"""Query Census 2 stores through a catalog and return one AnnData.

The SDK runs on the user's machine. It reads the catalog and excludes the datasets that
cannot match. It then opens the remaining stores and concatenates the selected cells by
gene position (design 4.1 to 4.4).

    import census2_sdk

    census = census2_sdk.open_census(sources=["catalog.parquet"])
    hs = census.homo_sapiens
    adata = census.get_anndata(hs.obs.cell_type.isin(["neutrophil"]) & (hs.obs.nnz >= 500))
    adata = census.get_anndata(hs.where({"tissue": "lung"}))
    adata = census.get_anndata(hs)  # every cell of the organism

A query selects the cells of one organism. The columns of `hs.obs` have pandas names
and give pyarrow expressions (`census2_sdk.query`). pyarrow is the filter language and
the evaluator, with SQL's three-valued logic for nulls. `hs.where` takes a mapping of a
catalog column to one value or to a list of values, and selects a cell when every
column holds one of its values.

A query also records the columns that it reads, and, from its top-level `&` of `==`,
`isin`, and `has_any` terms on catalog columns, the allowed values of each column (its
hints). The catalog excludes a dataset whose value sets prove that no cell matches the
hints. Then pyarrow evaluates the query on the combinations of each remaining dataset's
value sets (`census2_sdk.value_sets`): it excludes a dataset where no cell can match, and
proves a dataset whole where every cell matches. For a whole dataset, the SDK reads no
filter column of the store and keeps all of its cells. For every other candidate store,
the SDK reads only the query's columns, as codes for a categorical column, and
evaluates the expression with pyarrow (`census2_sdk.obs_read`).

The result carries the counts in `X` and the pinned gene reference as `var`. Its `obs`
has the columns that all stores with a selected cell hold. A categorical column has
only the categories that the selected cells use. `uns["census2"]` has the sources, the
number of catalog rows for the organism, the number of stores opened, and the number of
them that the catalog proved whole. The cell key `(dataset_id, observation_joinid)` is
in `obs`, and `obs_names` is `<dataset_id>_<observation_joinid>`. When no cell matches,
`get_anndata` returns None.

The stores hold only the counts. `layers=["normalized"]` adds the Census 1
normalized values in `layers["normalized"]`, calculated on read.

The writer projects each store onto the reference. Thus every selected store has
the same genes in `var`. The SDK checks that every store has the same gene identifiers
in the same order, and concatenates by position. A store's `var` also marks the
dataset's assayed genes in a `measured` column, and the result leaves that column out,
because the result spans datasets.

`count_cells` counts the selected cells and opens no store. It sums the cell counts
that each catalog row holds, so it takes only a query that is exactly its hints, on the
columns of `catalog.COUNTED`, `dataset_id`, and `is_primary_data`.

    n = census.count_cells(hs.obs.cell_type == "neutrophil")

Each store's `obs` holds `dataset_id`, which the writer adds. The SDK refuses a
store whose `dataset_id` differs from the dataset of its catalog row.
"""

import os
from collections.abc import Sequence
from concurrent.futures import ThreadPoolExecutor
from typing import Any, cast

import anndata
import numpy as np
import pandas as pd
import scipy.sparse as sp
import zarr
from anndata.io import read_elem, sparse_dataset

from census2_sdk import catalog as _catalog
from census2_sdk import normalize as _normalize
from census2_sdk import obs_read
from census2_sdk import value_sets as _value_sets
from census2_sdk.catalog import Filter, members, open_store, unrecorded
from census2_sdk.gene_reference import MEASURED_COLUMN, ORGANISMS
from census2_sdk.obs_columns import Obs
from census2_sdk.query import HasQuery, Query, QueryLike, as_query

__all__ = [
    "Census",
    "Organism",
    "Query",
    "SourceConflict",
    "open_census",
]

PathLike = str | os.PathLike[str]


class SourceConflict(ValueError):
    """The sources or the selected stores do not compose (design 4.3, 4.4)."""


class Organism(HasQuery):
    """The `obs` columns of one organism, and the query of all of its cells.

    `census.homo_sapiens` gives one. Use it as a query to select every cell, build a
    query from the typed columns of `obs`, or use `where` for a mapping.
    """

    def __init__(self, name: str) -> None:
        if name not in ORGANISMS:
            raise ValueError(f"no organism {name!r}; use one of {sorted(ORGANISMS)}")
        self.name = name
        self.obs = Obs(name)

    def query(self) -> Query:
        """The query of every cell of the organism."""
        return Query(self.name, None)

    def where(self, obs_filter: Filter) -> Query:
        """The query of a mapping of a catalog column to one value or to a list.

        A cell matches when every column holds one of its values. For the
        multi-valued `disease` columns, a member of the cell's value matches.
        """
        query = self.query()
        for column, values in _catalog.canonicalize(obs_filter).items():
            labels = sorted(values)
            col = getattr(self.obs, column)
            if column in _catalog.MULTIVALUED:
                term = col.has_any(labels)
            elif len(labels) == 1:
                term = col == labels[0]
            else:
                term = col.isin(labels)
            query = query & term
        return query

    def __repr__(self) -> str:
        return f"<Organism {self.name}>"


class Census:
    """`Census` holds a set of opened sources.

    Each call to `get_anndata` or `count_cells` runs one query over them.
    """

    def __init__(self, sources: Sequence[PathLike]) -> None:
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

    @property
    def homo_sapiens(self) -> Organism:
        return Organism("homo_sapiens")

    @property
    def mus_musculus(self) -> Organism:
        return Organism("mus_musculus")

    def _candidates(self, query: QueryLike) -> pd.DataFrame:
        """The catalog rows that can hold a selected cell, with a `whole` column.

        Two tests run, and a row must pass both. First the query's hints exclude a
        row (`catalog.excludes`), which costs little. Then pyarrow evaluates the
        query on every combination of each remaining row's value sets
        (`census2_sdk.value_sets`), which also tests `|` and `~`. `whole` is true when
        every combination matches; for a row whose table passes the size limit, when
        the query is exactly its hints and the value sets lie inside them
        (`catalog.covers`). Both proofs of `whole` rely on the store
        format's rule that the indexed columns hold no null.
        """
        q = as_query(query)
        catalog_hints = {
            column: values
            for column, values in q.hints.items()
            if column == "dataset_id" or column in _catalog.INDEXED
        }
        rows = _catalog.exclude(self.catalog, q.organism, catalog_hints)
        records = cast("list[dict[str, Any]]", rows.to_dict("records"))
        keep, whole, tested = _value_sets.test(records, q)
        hinted = np.array(
            [
                q.expression is None or (q.exact and _catalog.covers(row, q.hints))
                for row in records
            ],
            dtype=bool,
        )
        whole = np.where(tested, whole, hinted)
        return rows[keep].assign(whole=whole[keep])

    def get_anndata(
        self,
        query: QueryLike,
        *,
        layers: Sequence[str] = (),
        max_workers: int = 1,
    ) -> anndata.AnnData | None:
        """Return the cells that `query` selects as one AnnData.

        Return None when no cell matches, whether or not the catalog excluded
        every store.

        `layers` names the layers to add next to the counts. The one layer is
        `normalized`, which the SDK calculates from the counts with the Census 1
        algorithm, without its rounding (`census2_sdk.normalize`).

        `max_workers` is the largest number of stores that the SDK reads at the
        same time, and it must be a positive integer. The result is the same for
        every allowed value: the stores keep their catalog order.
        """
        unknown = sorted(set(layers) - {"normalized"})
        if unknown:
            raise ValueError(f"no layer {unknown}; the one layer is 'normalized'")
        q = as_query(query)
        rows = self._candidates(q)
        _check_composition(rows, ORGANISMS[q.organism]["version"])
        records = cast("list[dict[str, Any]]", rows.to_dict("records"))
        uris = [str(u) for u in rows["store_uri"]]
        bucket_clients = _catalog.clients(uris)
        # The default is one store at a time. Raise it after the store-count benchmark.
        with ThreadPoolExecutor(max_workers=max_workers) as pool:
            parts = list(
                pool.map(
                    lambda row: _read_store(row, q, layers, bucket_clients), records
                )
            )
        kept = [(p, u) for p, u in zip(parts, uris, strict=True) if p is not None]
        if not kept:
            return None
        _check_var([p for p, _ in kept], [u for _, u in kept])
        result = anndata.concat([p for p, _ in kept], join="inner")
        if not result.obs_names.is_unique:
            raise SourceConflict("the cell key is not unique across the result")
        # Each part carries only the var index. The stores record one reference
        # version and one index, so the first store supplies the var columns,
        # except the assayed genes of its own dataset.
        first = kept[0][1]
        var = cast("pd.DataFrame", read_elem(open_store(first, bucket_clients)["var"]))
        result.var = var.drop(columns=MEASURED_COLUMN, errors="ignore")
        result.uns["census2"] = {
            "sources": self.sources,
            "stores_in_catalog": int((self.catalog["organism"] == q.organism).sum()),
            "stores_opened": len(rows),
            "stores_whole": int(rows["whole"].sum()),
        }
        return result

    def _get_obs(
        self,
        query: QueryLike,
        *,
        columns: Sequence[str] | None = None,
        max_workers: int = 1,
    ) -> pd.DataFrame:
        """The `obs` of the cells that `query` selects, with no read of `X`.

        The frame is the `obs` of `get_anndata`, with only `columns` when given. A
        query that selects no cell gives a frame with no rows. The measure tools use
        it.
        """
        q = as_query(query)
        rows = self._candidates(q)
        records = cast("list[dict[str, Any]]", rows.to_dict("records"))
        bucket_clients = _catalog.clients(str(u) for u in rows["store_uri"])
        with ThreadPoolExecutor(max_workers=max_workers) as pool:
            parts = [
                part
                for part in pool.map(
                    lambda row: _read_obs(row, q, columns, bucket_clients), records
                )
                if part is not None
            ]
        if not parts:
            return pd.DataFrame(columns=list(columns or []))
        # anndata unites the categories of each categorical column.
        frames = [anndata.AnnData(obs=frame) for frame, _, _ in parts]
        return cast("pd.DataFrame", anndata.concat(frames, join="inner").obs)

    def count_cells(self, query: QueryLike) -> int:
        """Return the number of cells that `query` selects.

        The SDK sums the cell counts in the catalog rows and opens no store. A
        cell matches by the rule of `get_anndata`. Thus the count equals the `n_obs`
        of the result of `get_anndata` when there is one. The count does not apply the checks of
        design 4.3, because a count joins no matrices.

        The query must be exactly its hints (`Query.exact`): a `&` of `==`, `isin`,
        and `has_any` terms on the columns of `catalog.COUNTED` and `dataset_id`, and
        of `is_primary_data`. Any other query raises ValueError. A catalog without
        the `cell_counts` column raises ValueError.
        """
        q = as_query(query)
        if not q.exact:
            raise ValueError(
                "count_cells sums the catalog's cell counts, so it takes only a & of "
                "==, isin, and has_any terms and is_primary_data; use get_anndata "
                f"for {q.text}"
            )
        counted = {"dataset_id", "is_primary_data", *_catalog.COUNTED}
        uncounted = sorted(set(q.hints) - counted)
        if uncounted:
            raise ValueError(f"the catalog counts no {uncounted}; use get_anndata")
        canonical = {c: v for c, v in q.hints.items() if c != "is_primary_data"}
        counts = self._cell_counts()
        candidates = _catalog.exclude(self.catalog, q.organism, canonical)
        rows = counts[
            counts["dataset_version_id"].isin(candidates["dataset_version_id"])
        ]
        primary = q.hints.get("is_primary_data")
        if primary is not None:
            rows = rows[rows["is_primary_data"].astype(bool).isin(primary)]
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


def open_census(sources: Sequence[PathLike]) -> Census:
    """Open catalogs."""
    return Census(sources)


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


def _hit(value: Any, column: str, wanted: frozenset[str]) -> bool:
    """True when one cell's value holds a wanted member. A null holds none."""
    return not pd.isna(value) and bool(wanted & set(members(str(value), column)))


def _n_obs(group: zarr.Group) -> int:
    """The number of cells of a store, from the metadata of `X`."""
    x = group["X"]
    if isinstance(x, zarr.Array):
        return int(x.shape[0])
    return int(cast("list[int]", x.attrs["shape"])[0])


def _listed(row: dict[str, Any]) -> dict[str, list[str]]:
    """The catalog's sorted value set of each indexed column that the row records."""
    return {
        name: sorted(str(v) for v in row[name])
        for name in _catalog.INDEXED
        if not unrecorded(row.get(name))
    }


def _select(
    row: dict[str, Any], q: Query, bucket_clients: _catalog.Clients | None
) -> tuple[zarr.Group, np.ndarray, dict[str, obs_read.ColumnData]] | None:
    """Open a store and find its selected cells. Return None when it selects none.

    A row that the catalog proved whole keeps every cell and reads no column. For
    the other rows, the SDK reads the query's columns and evaluates the expression.
    A column that the query names and the store does not hold raises `ValueError`.
    The columns that it read go back with the positions of the kept cells, for the
    result.
    """
    uri = str(row["store_uri"])
    group = open_store(uri, bucket_clients)
    obs = group["obs"]
    if not isinstance(obs, zarr.Group):
        raise TypeError(f"{uri}: obs must be an encoded dataframe")
    n = _n_obs(group)
    if row["whole"]:
        return (group, np.arange(n), {}) if n else None
    held = set(cast("list[str]", obs.attrs["column-order"]))
    names = sorted(q.columns)
    missing = [c for c in names if c not in held]
    if missing:
        raise ValueError(f"{uri}: obs has no column {missing}, which {q.text} reads")
    data = obs_read.read_columns(obs, names, _listed(row))
    try:
        keep = obs_read.evaluate(q.expression, data, n)
    except (TypeError, ValueError, NotImplementedError) as error:
        raise ValueError(f"{uri}: pyarrow cannot evaluate {q.text}: {error}") from error
    if not keep.size:
        return None
    return group, keep, data


def _obs_frame(
    row: dict[str, Any],
    obs: zarr.Group,
    keep: np.ndarray,
    data: dict[str, obs_read.ColumnData],
    columns: Sequence[str] | None = None,
) -> pd.DataFrame:
    """The `obs` of the kept cells, with only `columns` when given.

    The SDK reads each other column as codes, and takes the categories of a catalog
    column from the row. It does not read the index (the source cell names),
    because the result replaces the index with the cell key. It reads the cell key
    for the kept cells only.
    """
    uri = str(row["store_uri"])
    order = cast("list[str]", obs.attrs["column-order"])
    key = "observation_joinid"
    if "dataset_id" not in order:
        raise ValueError(f"{uri}: obs has no dataset_id column")
    if key not in order:
        raise ValueError(f"{uri}: obs has no observation_joinid column")
    wanted = [c for c in order if columns is None or c in columns or c == key]
    names = [c for c in order if c != key and (c in wanted or c == "dataset_id")]
    unread = [c for c in names if c not in data]
    data = {**data, **obs_read.read_columns(obs, unread, _listed(row))}
    dataset_id = str(row["dataset_id"])
    _check_dataset_id(data["dataset_id"], dataset_id, uri)
    frame = pd.DataFrame({c: data[c].take(keep) for c in names})
    frame[key] = obs_read.take_rows(obs[key], keep)
    _check_cell_key(frame, uri)
    frame.index = pd.Index(dataset_id + "_" + frame[key].astype(str))
    return frame[wanted]


def _read_obs(
    row: dict[str, Any],
    q: Query,
    columns: Sequence[str] | None,
    bucket_clients: _catalog.Clients | None,
) -> tuple[pd.DataFrame, zarr.Group, np.ndarray] | None:
    """The `obs` of the selected cells of one store, or None when it selects none."""
    selected = _select(row, q, bucket_clients)
    if selected is None:
        return None
    group, keep, data = selected
    obs = group["obs"]
    assert isinstance(obs, zarr.Group)  # _select checked it
    return _obs_frame(row, obs, keep, data, columns), group, keep


def _read_store(
    row: dict[str, Any],
    q: Query,
    layers: Sequence[str] = (),
    bucket_clients: _catalog.Clients | None = None,
) -> anndata.AnnData | None:
    """The selected cells of one store, or None when it selects none.

    `X` holds the counts on the reference genes. `var` holds only the index, the
    `feature_id`s, which is all that `_check_var` compares. `get_anndata` reads the
    var columns from one store. `layers` names the layers to calculate from the
    counts. `bucket_clients` holds the obstore clients of the query
    (`catalog.clients`). A store that selects no cell reads nothing more.
    """
    read = _read_obs(row, q, None, bucket_clients)
    if read is None:
        return None
    obs, group, keep = read
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
    return anndata.AnnData(X=x, obs=obs, var=var, layers=computed)


def _check_dataset_id(
    column: obs_read.ColumnData, dataset_id: str, store_uri: str
) -> None:
    """The store's `dataset_id` column names the dataset of its catalog row.

    The store format gives a categorical column no unused category, so the store's
    categories are the values of every cell, and a code of -1 is a null.
    """
    n = len(column.codes if isinstance(column, obs_read.Codes) else column.values)
    if isinstance(column, obs_read.Codes):
        found = {str(c) for c in column.categories}
        if (column.codes < 0).any():
            found.add("nan")
    else:
        found = {str(v) for v in pd.Series(column.take(np.arange(n))).unique()}
    if found != {dataset_id}:
        raise ValueError(
            f"{store_uri}: obs dataset_id holds {sorted(found)}, and the catalog row "
            f"names {dataset_id}"
        )


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


# These limits are fixed. Tune them when stores change their row order or chunk size.
MAX_RUNS = 256
RUN_THREADS = 16  # range reads in flight for one store


def _rows(node: zarr.Array[Any] | zarr.Group, keep: np.ndarray) -> sp.csr_matrix:
    """The rows `keep` of a CSR matrix, read without the other rows.

    The plain call is `sparse_dataset(node)[keep]`. For row positions, anndata makes one
    array position for each stored value of each kept row. zarr then selects the values
    by position. The CPU time and the memory of that read increase with the stored values.

    The stored values of adjacent kept rows are one span of `data` and of `indices`.
    Thus this function reads one slice of each array for each run of adjacent kept rows.
    Above `MAX_RUNS` runs, the SDK uses the plain call, which makes fewer requests. A
    matrix in another encoding, such as CSC, also uses the plain call.
    """
    if not isinstance(node, zarr.Group):
        raise TypeError(f"{node.store_path} must be an encoded sparse matrix")
    arrays = [node[name] for name in ("data", "indices", "indptr")]
    if not all(isinstance(a, zarr.Array) for a in arrays):
        raise TypeError(f"{node.store_path} must hold data, indices, and indptr arrays")
    data, indices, indptr_node = cast("list[zarr.Array[Any]]", arrays)
    n_vars = int(cast("list[int]", node.attrs["shape"])[1])
    if not len(keep):
        return sp.csr_matrix((0, n_vars), dtype=data.dtype)
    if node.attrs.get("encoding-type") != "csr_matrix":
        return sp.csr_matrix(cast("sp.csr_matrix", sparse_dataset(node)[keep]))

    # A new run starts at each kept row that does not follow the kept row before it.
    breaks = np.flatnonzero(np.diff(keep) != 1) + 1
    if len(breaks) + 1 > MAX_RUNS:
        return sp.csr_matrix(cast("sp.csr_matrix", sparse_dataset(node)[keep]))
    runs = np.split(keep, breaks)

    indptr = np.asarray(indptr_node[:], dtype=np.int64)
    spans = [slice(int(indptr[run[0]]), int(indptr[run[-1] + 1])) for run in runs]

    # A zarr slice read blocks the thread that calls it, so one thread reads one span at
    # a time. This pool keeps as many as RUN_THREADS spans of the store in flight.
    # `get_anndata` calls this function from one thread for each store that it reads.
    # Thus a query has no more than max_workers * RUN_THREADS span reads in flight.
    # On 103 stores, the pool decreases the wall time of a query by 5 to 60 percent.
    # It adds as much as 10 percent of CPU time.
    with ThreadPoolExecutor(min(RUN_THREADS, 2 * len(spans))) as pool:
        # `map` submits each read at once, so both arrays are in flight before the wait.
        pending = [pool.map(array.__getitem__, spans) for array in (data, indices)]
        kept_data, kept_indices = (
            np.concatenate([np.asarray(part) for part in parts]) for parts in pending
        )

    row_nnz = np.diff(indptr)[keep]
    return sp.csr_matrix(
        (kept_data, kept_indices, np.r_[0, np.cumsum(row_nnz)]),
        shape=(len(keep), n_vars),
    )
