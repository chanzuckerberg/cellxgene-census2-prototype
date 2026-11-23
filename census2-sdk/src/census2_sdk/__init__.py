"""Query Census 2 stores through a catalog and return one AnnData.

The SDK runs on the user's machine. It reads the catalog and excludes the datasets that
cannot match. It then opens the remaining stores and concatenates the selected cells by
gene position (design 4.1 to 4.4).

    import census2_sdk

    census = census2_sdk.open_census(sources=["catalog.parquet"])
    adata = census.get_anndata("homo_sapiens", obs_filter={"tissue": "lung"})

`obs_filter` maps a catalog column to one value or to a list of values. The SDK selects
a cell when every column holds one of its values. The columns are the 18 indexed
columns of design 2.2 and `dataset_id`. A filter grammar is open (design section 8).

The result carries the counts in `X` and the pinned gene reference as `var`. Its `obs`
has the columns that all selected stores hold. `uns["census2"]` has the sources, the
number of catalog rows for the organism, the number of datasets read, and the number
of stores opened. The cell key `(dataset_id, observation_joinid)` is in `obs`, and
`obs_names` is `<dataset_id>_<observation_joinid>`. When no cell matches,
`get_anndata` returns None.

The stores hold only the counts. `layers=["normalized"]` adds the Census 1
normalized values in `layers["normalized"]`, calculated on read.

The writer projects each store onto the reference. Thus every selected store has
the same genes in `var`. The SDK checks that every store has the same gene identifiers
in the same order, and concatenates by position. A store's `var` also marks the
dataset's assayed genes in a `measured` column, and the result leaves that column out,
because the result spans datasets.

`get_obs` returns the `obs` of the selected cells as a DataFrame. It reads no `X`,
and `column_names` limits the columns that it reads.

    obs = census.get_obs("homo_sapiens", {"tissue": "lung"}, column_names=["sex"])

`get_var` returns the `var` of an organism, which is the pinned gene reference as the
stores hold it. It reads one store.

    var = census.get_var("homo_sapiens")

`get_presence_matrix` returns the genes that each dataset of an organism measured, as
a boolean DataFrame of `dataset_id` by `feature_id`. It opens every store of the
organism.

    presence = census.get_presence_matrix("homo_sapiens")

`count_cells` counts the selected cells and opens no store. It sums the cell counts
that each catalog row holds. Its filter takes the columns of
`catalog.COUNTED` and `dataset_id`.

    n = census.count_cells("homo_sapiens", {"cell_type": "neutrophil"})

Each store's `obs` holds `dataset_id`, which the writer adds. The SDK refuses a
store whose `dataset_id` differs from the dataset of its catalog row.

A catalog row can also name a consolidated store, a second copy of several small
datasets in one store, with the rows of its dataset there (`consolidated_store_uri`,
`consolidated_row_start`, `consolidated_row_stop`). The SDK then reads the dataset from
the consolidated store, and opens each store once for all of its selected datasets.
Exclusion stays per dataset, so a consolidated store opens only when one of its
datasets can match. The SDK reads only the rows of the selected datasets, and the
result is the result of the per-dataset stores. When the selected datasets hold at
least `WHOLE_STORE_SHARE` of the rows of a consolidated store, the SDK reads the `obs`
columns of that store whole, once, and cuts each dataset's rows from memory.
`open_census(..., consolidated=False)` reads the per-dataset stores, which stay the
primary copy.
"""

import os
from collections.abc import Sequence
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from typing import Any, cast

import anndata
import numpy as np
import pandas as pd
import pyarrow as pa  # type: ignore[import-untyped]  # pyarrow has no types
import scipy.sparse as sp
import zarr
from anndata.io import read_elem, sparse_dataset

from census2_sdk import catalog as _catalog
from census2_sdk import normalize as _normalize
from census2_sdk.catalog import (
    CONSOLIDATED,
    CanonicalFilter,
    Filter,
    members,
    open_store,
    unrecorded,
)
from census2_sdk.gene_reference import MEASURED_COLUMN, ORGANISMS

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

    Each call to `get_anndata`, `get_obs`, `get_var`, `get_presence_matrix`, or
    `count_cells` runs one query over them.
    """

    def __init__(
        self, sources: Sequence[PathLike], *, consolidated: bool = True
    ) -> None:
        self.sources = [str(s) for s in sources]
        # Read a dataset from the consolidated store that its catalog row names, if any.
        self.consolidated = consolidated
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

    def get_anndata(
        self,
        organism: str,
        obs_filter: Filter | None = None,
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

        `max_workers` is the largest number of reads that the SDK runs at the
        same time, and it must be a positive integer. The SDK first opens each
        store, then reads each selected dataset, and both steps use up to
        `max_workers` threads. The datasets of one consolidated store thus read at
        the same time as well. The result is the same for every allowed value:
        the datasets keep their catalog order.
        """
        unknown = sorted(set(layers) - {"normalized"})
        if unknown:
            raise ValueError(f"no layer {unknown}; the one layer is 'normalized'")
        canonical = _catalog.canonicalize(obs_filter)
        rows = _catalog.exclude(self.catalog, organism, canonical)
        _check_composition(rows, ORGANISMS[organism]["version"])

        reads = _reads(rows, self.consolidated)
        # The default is one store at a time. Raise it after the store-count benchmark.
        with ThreadPoolExecutor(max_workers=max_workers) as pool:
            opened = _open_stores(reads, pool, with_var=True)
            uris = list(opened)
            parts = list(
                pool.map(
                    lambda r: _read_rows(opened[r.uri], r, canonical, layers), reads
                )
            )
        _check_var(parts, [r.uri for r in reads])
        if not sum(part.n_obs for part in parts):
            return None
        # Zero-row parts stay in, so the obs schema is the columns that all
        # opened stores hold.
        result = anndata.concat(parts, join="inner")
        if not result.obs_names.is_unique:
            raise SourceConflict("the cell key is not unique across the result")
        # Each part carries only the var index. The stores record one reference
        # version and one index, so the first store supplies the var columns,
        # except the assayed genes of its own dataset. A consolidated store holds
        # those in `varm`, so its `var` has no such column.
        var = cast("pd.DataFrame", read_elem(opened[uris[0]].group["var"]))
        result.var = var.drop(columns=MEASURED_COLUMN, errors="ignore")
        result.uns["census2"] = {
            "sources": self.sources,
            "datasets_in_catalog": int((self.catalog["organism"] == organism).sum()),
            "datasets_read": len(reads),
            "stores_opened": len(uris),
        }
        return result

    def get_obs(
        self,
        organism: str,
        obs_filter: Filter | None = None,
        *,
        column_names: Sequence[str] | None = None,
        max_workers: int = 1,
    ) -> pd.DataFrame | None:
        """Return the `obs` of the selected cells of `organism`, without a read of `X`.

        The frame equals the `obs` of the result of `get_anndata` for the same
        filter. It has the same cells and the same index, and the SDK applies the
        same checks. Return None when no cell matches.

        `column_names` names the columns of the result, in that order. From each
        store, the SDK then reads only these columns, the filter columns,
        `dataset_id`, and `observation_joinid`. A store without a named column
        raises ValueError. Without `column_names`, the result has the columns that
        all opened stores hold.

        `max_workers` has the same function as in `get_anndata`. A dataset that a
        catalog row places in a consolidated store reads from there, as in
        `get_anndata`.
        """
        canonical = _catalog.canonicalize(obs_filter)
        rows = _catalog.exclude(self.catalog, organism, canonical)
        _check_composition(rows, ORGANISMS[organism]["version"])

        reads = _reads(rows, self.consolidated)

        def read(r: _Read) -> anndata.AnnData:
            obs = _read_obs(opened[r.uri], r, canonical, column_names)[0]
            return anndata.AnnData(obs=obs)

        key = "observation_joinid"
        columns = (
            None
            if column_names is None
            else {*column_names, *canonical, "dataset_id", key}
        )
        with ThreadPoolExecutor(max_workers=max_workers) as pool:
            opened = _open_stores(reads, pool, with_var=False, columns=columns)
            parts = list(pool.map(read, reads))
        if not sum(part.n_obs for part in parts):
            return None
        # `anndata.concat` gives the frame the dtypes of the `get_anndata` result.
        result = cast("pd.DataFrame", anndata.concat(parts, join="inner").obs)
        if not result.index.is_unique:
            raise SourceConflict("the cell key is not unique across the result")
        return result

    def get_var(self, organism: str) -> pd.DataFrame:
        """Return the `var` of `organism`: the pinned gene reference, as the stores hold it.

        The frame is the `var` of a `get_anndata` result. Every store of one organism
        holds the same reference (design 1.5), so the SDK reads the `var` of one
        store. The frame leaves out the `measured` column, which describes one dataset.

        Every catalog row of the organism must record the gene reference version
        that the SDK pins. Sources that hold no store of the organism raise
        ValueError.
        """
        rows = _catalog.exclude(self.catalog, organism, None)
        if not len(rows):
            raise ValueError(f"the sources hold no store of {organism}")
        _check_reference(rows, ORGANISMS[organism]["version"])
        uri = str(rows["store_uri"].iloc[0])
        return _read_var(uri, _catalog.clients([uri]))

    def get_presence_matrix(
        self, organism: str, *, max_workers: int = 1
    ) -> pd.DataFrame:
        """Return the genes that each dataset of `organism` measured.

        The result is a boolean DataFrame. It has one row for each store of the
        organism, in catalog order, and the index is the `dataset_id`. It has one
        column for each gene of `get_var`, in reference order, and the column
        labels are the `feature_id`s. A value is True when the dataset measured
        the gene.

        The call opens every store of the organism, and it reads only the
        `measured` column of `var` from each. `get_var` supplies the genes, so the
        first store opens one more time. `max_workers` is the largest number of
        stores that the SDK reads at the same time.

        The SDK applies the checks of `get_var`. It also checks that each store
        holds one `measured` value for each gene of the reference. It does not
        compare the gene identifiers of each store.
        """
        genes = self.get_var(organism).index
        rows = _catalog.exclude(self.catalog, organism, None)
        uris = list(rows["store_uri"].astype(str))
        bucket_clients = _catalog.clients(uris)
        with ThreadPoolExecutor(max_workers=max_workers) as pool:
            measured = list(
                pool.map(
                    lambda uri: _read_measured(uri, len(genes), bucket_clients), uris
                )
            )
        datasets = list(rows["dataset_id"].astype(str))
        return pd.DataFrame(np.stack(measured), index=datasets, columns=genes)

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
        uncounted = sorted(set(canonical) - {"dataset_id", *_catalog.COUNTED})
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


def open_census(sources: Sequence[PathLike], *, consolidated: bool = True) -> Census:
    """Open catalogs.

    With `consolidated` false, the SDK reads every dataset from its per-dataset store,
    and ignores the consolidated stores that the catalog names.
    """
    return Census(sources, consolidated=consolidated)


def _check_composition(rows: pd.DataFrame, pinned: str) -> None:
    """Design 4.3: the checks before concatenation. Name the conflicting values.

    Every row must record the pinned gene reference (`_check_reference`).

    Spatial and non-spatial stores do not compose. A row that records no
    `spatial` value cannot prove that its store composes with another, so a
    result of two or more stores refuses it. One store always composes with
    itself. A catalog without the column records no value.
    """
    _check_reference(rows, pinned)
    if len(rows) < 2:
        return
    spatial = _recorded(rows, "spatial")
    unknown = [d for d, v in zip(rows["dataset_id"], spatial, strict=True) if v is None]
    if unknown:
        raise SourceConflict(f"no spatial value is recorded for {sorted(unknown)}")
    if len({bool(v) for v in spatial}) > 1:
        raise SourceConflict("spatial and non-spatial stores in one result")


def _check_reference(rows: pd.DataFrame, pinned: str) -> None:
    """Every row records the `gene_reference_version` that this SDK pins for the organism.

    The SDK does not remap a store onto another reference. A row that records
    none has no summary record that names the reference of its store. A catalog
    without the column records none.
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
        if isinstance(col.dtype, pd.CategoricalDtype):
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


@dataclass(frozen=True)
class _Read:
    """One selected dataset, and where the SDK reads it.

    `rows` is the dataset's row range in a consolidated store, or None for a
    per-dataset store, which holds the dataset whole.
    """

    row: dict[str, Any]
    uri: str
    rows: tuple[int, int] | None


@dataclass
class _Opened:
    """One open store, and what all of its datasets share.

    `var` holds only the index, or is None when the query reads no `X`. `n_obs` is
    the number of rows of the store. `preloaded` holds the `obs` columns of a
    consolidated store that the open read whole (`_read_whole`).
    `categories` holds the categories of each categorical `obs` column of a
    consolidated store, read once for all the datasets of the store.
    """

    uri: str
    group: zarr.Group
    var: pd.DataFrame | None
    n_obs: int
    categories: dict[str, pd.Index] = field(default_factory=dict)
    preloaded: dict[str, Any] = field(default_factory=dict)


def _reads(rows: pd.DataFrame, consolidated: bool) -> list[_Read]:
    """The read of each selected dataset, in catalog order.

    A row that names a consolidated store is read there, unless `consolidated` is false.
    A catalog without the columns names none.
    """
    records = cast("list[dict[str, Any]]", rows.to_dict("records"))
    use = consolidated and all(c in rows.columns for c in CONSOLIDATED)
    reads = []
    for row in records:
        uri = row.get(CONSOLIDATED[0]) if use else None
        if unrecorded(uri):
            reads.append(_Read(row, str(row["store_uri"]), None))
        else:
            span = (int(row[CONSOLIDATED[1]]), int(row[CONSOLIDATED[2]]))
            reads.append(_Read(row, str(uri), span))
    return reads


def _group(parent: zarr.Group, name: str, uri: str, kind: str) -> zarr.Group:
    """The child group `name`, or a TypeError that names the store and the node."""
    node = parent[name]
    if not isinstance(node, zarr.Group):
        raise TypeError(f"{uri}: {name} must be {kind}")
    return node


# A consolidated store reads its `obs` columns whole, once, when the selected datasets
# hold at least this share of its rows. A member read then cuts its rows from memory.
# Otherwise each member reads its own rows, and each partial read of a sharded array
# also reads the shard index. The share is a first guess, not a measured optimum.
WHOLE_STORE_SHARE = 0.5


def _open_stores(
    reads: list[_Read],
    pool: ThreadPoolExecutor,
    *,
    with_var: bool,
    columns: set[str] | None = None,
) -> dict[str, _Opened]:
    """Open each store of `reads` once, in order, with the threads of `pool`.

    The stores of one bucket share one obstore client (`catalog.clients`). `with_var`
    reads the var index of each store, which `get_anndata` compares. `columns` names
    the `obs` columns that the query reads, or None for all of them.
    """
    uris = list(dict.fromkeys(r.uri for r in reads))
    selected: dict[str, int] = {}
    for r in reads:
        if r.rows is not None:
            selected[r.uri] = selected.get(r.uri, 0) + r.rows[1] - r.rows[0]
    bucket_clients = _catalog.clients(uris)

    def one(uri: str) -> _Opened:
        return _open(
            uri,
            uri in selected,
            bucket_clients,
            with_var=with_var,
            selected_rows=selected.get(uri, 0),
            columns=columns,
        )

    return dict(zip(uris, pool.map(one, uris), strict=True))


def _open(
    uri: str,
    consolidated: bool,
    bucket_clients: _catalog.Clients | None = None,
    *,
    with_var: bool = True,
    selected_rows: int = 0,
    columns: set[str] | None = None,
) -> _Opened:
    """Open a store, and read what its datasets share: the var index, and categories.

    `var` holds only the index, the `feature_id`s, which is all that `_check_var`
    compares. `get_anndata` reads the var columns from one store. The number of
    rows comes from the attributes of `X`, which the open already read. The
    categories of a consolidated store's `obs` columns serve each member's read.
    `bucket_clients` holds the obstore clients of the query (`catalog.clients`).

    When the selected datasets of a consolidated store hold `selected_rows` of its
    rows, at least `WHOLE_STORE_SHARE` of them, the open also reads the `obs`
    columns that the query needs (`columns`, or all) whole, into `preloaded`.
    """
    group = open_store(uri, bucket_clients)
    var = None
    if with_var:
        var_group = _group(group, "var", uri, "an encoded dataframe")
        index = read_elem(var_group[str(var_group.attrs["_index"])])
        var = pd.DataFrame(index=pd.Index(np.asarray(index).astype(str)))
    x = _group(group, "X", uri, "an encoded sparse matrix")
    opened = _Opened(uri, group, var, int(cast("list[int]", x.attrs["shape"])[0]))
    if not consolidated:
        return opened
    _check_record(opened)
    obs = _group(group, "obs", uri, "an encoded dataframe")
    names = [str(name) for name in cast("list[str]", obs.attrs["column-order"])]
    categorical = [
        name
        for name in names
        if isinstance(obs[name], zarr.Group) and _encoding(obs[name]) == "categorical"
    ]

    def categories(name: str) -> pd.Index:
        node = _group(obs, name, uri, "an encoded categorical")
        return cast("pd.Index", pd.Index(np.asarray(read_elem(node["categories"]))))

    whole = [
        name
        for name in names
        if selected_rows >= WHOLE_STORE_SHARE * opened.n_obs
        and (columns is None or name in columns)
    ]
    # Each column is a read of its own, so read them at the same time, as `_obs_rows`
    # reads the columns of a dataset.
    threads = max(1, min(len(categorical) + len(whole), OBS_THREADS))
    with ThreadPoolExecutor(max_workers=threads) as pool:
        pending = pool.map(categories, categorical)
        preloaded = pool.map(lambda name: _read_whole(obs[name]), whole)
        opened.categories = dict(zip(categorical, pending, strict=True))
        opened.preloaded = dict(zip(whole, preloaded, strict=True))
    return opened


def _read_whole(node: zarr.Array[Any] | zarr.Group) -> Any:
    """One `obs` column read whole, in the form that `_column_rows` cuts rows from.

    That is the array of an array column, the codes of a categorical column, and the
    values and the mask of a nullable string column. Another encoding reads with
    `read_elem`.
    """
    kind = _encoding(node)
    if isinstance(node, zarr.Array) and kind in ("array", "string-array"):
        return node[:]
    if isinstance(node, zarr.Group) and kind == "categorical":
        return np.asarray(cast("zarr.Array[Any]", node["codes"])[:])
    if isinstance(node, zarr.Group) and kind == "nullable-string-array":
        values, mask = node["values"], node["mask"]
        if not (isinstance(values, zarr.Array) and isinstance(mask, zarr.Array)):
            raise TypeError(f"{node.store_path} must hold a values and a mask array")
        return np.asarray(values[:]), np.asarray(mask[:])
    return read_elem(node)


def _encoding(node: zarr.Array[Any] | zarr.Group) -> str:
    return str(node.attrs.get("encoding-type", ""))


def _check_record(opened: _Opened) -> None:
    """A consolidated store holds a consolidation record, in its root attributes."""
    namespace = opened.group.attrs.get("cellxgene-census2")
    if not isinstance(namespace, dict) or "consolidation" not in namespace:
        raise ValueError(
            f"{opened.uri}: no consolidation record; not a consolidated store"
        )


def _member_rows(opened: _Opened, read: _Read) -> tuple[int, int]:
    """The rows of the dataset in a consolidated store, which its record must list."""
    namespace = cast("dict[str, Any]", opened.group.attrs["cellxgene-census2"])
    listed = cast("list[dict[str, Any]]", namespace["consolidation"].get("members", []))
    version = str(read.row["dataset_version_id"])
    spans = [
        (m.get("row_start"), m.get("row_stop"))
        for m in listed
        if m.get("dataset_version_id") == version
    ]
    if read.rows is None or spans != [read.rows]:
        raise ValueError(
            f"{opened.uri}: the consolidation record does not list {version} at rows "
            f"{read.rows}"
        )
    return read.rows


def _read_obs(
    opened: _Opened,
    read: _Read,
    obs_filter: CanonicalFilter,
    column_names: Sequence[str] | None = None,
) -> tuple[pd.DataFrame, np.ndarray, int]:
    """The `obs` of the selected cells of one dataset, their rows, and the first row.

    The read takes only the rows of the dataset: all of a per-dataset store, or the
    member's rows of a consolidated store. The row positions count from the first
    row of the dataset. The index is the cell key, `<dataset_id>_<observation_joinid>`,
    and it has no name. `column_names` names the columns of the frame. Without it,
    the frame has every column of the store.
    """
    dataset_id = str(read.row["dataset_id"])
    if read.rows is None:
        where, (start, stop) = opened.uri, (0, opened.n_obs)
    else:
        where, (start, stop) = f"{opened.uri}[{dataset_id}]", _member_rows(opened, read)
    obs_group = _group(opened.group, "obs", opened.uri, "an encoded dataframe")
    # To read obs faster, the SDK decodes the per-cell string columns only where
    # the result needs them. It does not read the index (the source cell names, or
    # the cell key in a consolidated store), because the result sets the cell key. It
    # reads the cell key after the filter, for the kept cells only.
    columns = [str(c) for c in cast("list[str]", obs_group.attrs["column-order"])]
    key = "observation_joinid"
    named = columns if column_names is None else list(dict.fromkeys(column_names))
    absent = sorted(set(named) - set(columns))
    if absent:
        raise ValueError(f"{where}: obs has no column {absent}")
    # The filter and the dataset check need their columns also when the frame
    # leaves them out.
    needed = {*named, *obs_filter, "dataset_id"} - {key}
    names = [c for c in columns if c in needed]
    obs = _obs_rows(opened, obs_group, names, read.row, start, stop)
    _check_dataset_id(obs, dataset_id, where)
    keep = np.flatnonzero(_mask(obs, obs_filter))
    obs = obs.iloc[keep].copy()
    if key in columns:
        obs[key] = _read_column_rows(
            obs_group[key], start + keep, opened.preloaded.get(key)
        )
    _check_cell_key(obs, where)
    # pandas gives the index the name of the key column. `write_h5ad` refuses an
    # index name that is also a column with other values, so the index has no name.
    obs.index = pd.Index(dataset_id + "_" + obs[key].astype(str)).rename(None)
    return obs[named], keep, start


def _read_rows(
    opened: _Opened,
    read: _Read,
    obs_filter: CanonicalFilter,
    layers: Sequence[str] = (),
) -> anndata.AnnData:
    """The selected cells of one dataset. `X` holds the counts on the reference genes.

    The read takes only the rows of the dataset from `obs` and `X` (`_read_obs`).
    `layers` names the layers to calculate from the counts.
    """
    obs, keep, start = _read_obs(opened, read, obs_filter)
    x = _rows(opened.group["X"], keep, start)
    computed = {}
    if "normalized" in layers:
        assay = obs["assay_ontology_term_id"].astype(str)
        full_gene = assay.isin(_normalize.FULL_GENE_ASSAYS).to_numpy()
        length = None
        if full_gene.any():
            var_group = _group(opened.group, "var", opened.uri, "an encoded dataframe")
            length = np.asarray(read_elem(var_group["feature_length"]), dtype=np.int64)
        computed["normalized"] = _normalize.normalized(x, full_gene, length)

    return anndata.AnnData(X=x, obs=obs, var=opened.var, layers=computed)


# The obs columns of one dataset that the SDK reads at the same time.
OBS_THREADS = 8


def _obs_rows(
    opened: _Opened,
    group: zarr.Group,
    columns: list[str],
    row: dict[str, Any],
    start: int,
    stop: int,
) -> pd.DataFrame:
    """Rows `start` to `stop` of the `obs` columns `columns`, without the other rows.

    `row` is the catalog row of the dataset. A read of all the rows of a store reads
    each column whole with `_read_obs_column`, which gives the column of `read_elem`.
    A read of part of the rows gives each column the type that `read_elem` gives
    it, and a categorical column keeps only the categories that its rows use. A
    per-dataset store holds no unused category, because the writer writes none, so
    a member of a consolidated store gets the categories of its per-dataset store.
    """
    whole = (start, stop) == (0, opened.n_obs)

    def one(name: str) -> Any:
        if whole:
            return _read_obs_column(group[name], row.get(name), name)
        return _column_rows(opened, group, name, start, stop)

    # Each column is a read of its own, with a round trip or two, so read them at the
    # same time. The pool waits on no other task, so it cannot deadlock the pool of
    # `get_anndata`.
    with ThreadPoolExecutor(max_workers=max(1, min(len(columns), OBS_THREADS))) as pool:
        values = list(pool.map(one, columns))
    return pd.DataFrame(dict(zip(columns, values, strict=True)), columns=columns)


def _column_rows(
    opened: _Opened, group: zarr.Group, name: str, start: int, stop: int
) -> Any:
    """Rows `start` to `stop` of one `obs` column, with the type of `read_elem`.

    The encodings that anndata writes for the store's column types read in part:
    categorical, array, string array, and nullable string array. Another encoding
    reads whole, and then the rows are cut from it. A column that the open read
    whole (`opened.preloaded`) cuts the rows from memory.
    """
    node = group[name]
    kind = _encoding(node)
    held = opened.preloaded.get(name)
    if isinstance(node, zarr.Array) and kind in ("array", "string-array"):
        return node[start:stop] if held is None else held[start:stop]
    if isinstance(node, zarr.Group) and kind == "categorical":
        source = cast("zarr.Array[Any]", node["codes"]) if held is None else held
        codes = np.asarray(source[start:stop])
        categories = opened.categories[name]
        used = np.unique(codes[codes >= 0])
        # Renumber the codes over the used categories. A null keeps the code -1,
        # which indexes the -1 at the end.
        renumber = np.full(len(categories) + 1, -1, dtype=np.int64)
        renumber[used] = np.arange(len(used))
        return pd.Categorical.from_codes(
            renumber[codes],
            categories=categories[used],
            ordered=bool(node.attrs.get("ordered", False)),
        )
    if isinstance(node, zarr.Group) and kind == "nullable-string-array":
        return _read_column_rows(node, np.arange(start, stop), held)
    # `read_elem` is typed for every element kind; here it is an array.
    return cast("Any", read_elem(node) if held is None else held)[start:stop]


def _read_var(store_uri: str, bucket_clients: _catalog.Clients) -> pd.DataFrame:
    """The `var` of one store, without the assayed genes of its own dataset."""
    var = cast("pd.DataFrame", read_elem(open_store(store_uri, bucket_clients)["var"]))
    return var.drop(columns=MEASURED_COLUMN, errors="ignore")


def _read_measured(
    store_uri: str, n_genes: int, bucket_clients: _catalog.Clients
) -> np.ndarray:
    """The `measured` column of one store: True for each gene that the dataset assayed."""
    var_group = open_store(store_uri, bucket_clients)["var"]
    if not isinstance(var_group, zarr.Group):
        raise TypeError(f"{store_uri}: var must be an encoded dataframe")
    if MEASURED_COLUMN not in var_group:
        raise ValueError(f"{store_uri}: var has no {MEASURED_COLUMN} column")
    node = var_group[MEASURED_COLUMN]
    if not isinstance(node, zarr.Array) or node.dtype != np.bool_:
        raise TypeError(f"{store_uri}: var {MEASURED_COLUMN} must be a boolean array")
    if node.shape != (n_genes,):
        raise ValueError(
            f"{store_uri}: var {MEASURED_COLUMN} holds {node.size} values, and "
            f"the reference has {n_genes} genes"
        )
    return np.asarray(node[:])


def _read_obs_column(node: zarr.Array[Any] | zarr.Group, listed: Any, name: str) -> Any:
    """One `obs` column of a store. `listed` is the catalog row's entry for `name`.

    This function is a faster form of `read_elem(node)` and returns the same column.
    For a listed categorical column, the function reads only the codes and takes the
    categories from the row. The writer stores the sorted values of the column as its
    categories, and the catalog builder checks that order.

    A column keeps `read_elem` when the store has a different number of categories
    than the row lists. `dataset_id` keeps `read_elem`, because `_check_dataset_id`
    compares the store with the row.
    """
    if (
        name not in _catalog.INDEXED
        or unrecorded(listed)
        or not isinstance(node, zarr.Group)
        or node.attrs.get("encoding-type") != "categorical"
    ):
        return read_elem(node)
    codes, categories = node["codes"], node["categories"]
    values = sorted(str(v) for v in listed)
    if not (
        isinstance(codes, zarr.Array)
        and isinstance(categories, zarr.Array)
        and categories.shape == (len(values),)
    ):
        return read_elem(node)
    return pd.Categorical.from_codes(
        cast("Any", np.asarray(codes[:])),
        categories=pd.Index(values),
        ordered=bool(node.attrs.get("ordered", False)),
    )


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


def _read_column_rows(
    node: zarr.Array[Any] | zarr.Group, rows: np.ndarray, held: Any = None
) -> Any:
    """The values of one `obs` column at the row positions `rows`.

    For a nullable string column, the SDK makes strings for these rows only. The
    dtype is the one that `anndata` reads for the same column. The SDK builds the
    pandas array through pyarrow, which is faster than the pandas conversion.
    `held` is the values and the mask of a nullable string column that the open
    read whole (`_read_whole`), if any.
    """
    if node.attrs.get("encoding-type") != "nullable-string-array":
        return pd.array(cast("Any", read_elem(node)))[rows]
    if not isinstance(node, zarr.Group):
        raise TypeError(f"{node.store_path} must be an encoded group")
    na_value = np.nan if node.attrs.get("na-value", "NA") == "NaN" else pd.NA
    dtype = pd.StringDtype(na_value=na_value)
    if not rows.size:
        return pd.array([], dtype=dtype)
    # One run of rows reads as a slice, which zarr reads faster than a selection.
    contiguous = rows[-1] - rows[0] + 1 == rows.size
    selection: Any = slice(int(rows[0]), int(rows[-1]) + 1) if contiguous else rows
    if held is not None:
        strings, nulls = held[0][selection], held[1][selection]
    else:
        values, mask = node["values"], node["mask"]
        if not (isinstance(values, zarr.Array) and isinstance(mask, zarr.Array)):
            raise TypeError(f"{node.store_path} must hold a values and a mask array")
        strings = np.asarray(values.oindex[selection])
        nulls = np.asarray(mask.oindex[selection])
    try:
        arrow = pa.array(strings, mask=nulls)
    except pa.ArrowNotImplementedError:
        # pyarrow before 26 does not read the NumPy string dtype.
        arrow = pa.array(strings.astype(object), mask=nulls)
    return pd.array(arrow, dtype=dtype)


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


def _rows(
    node: zarr.Array[Any] | zarr.Group, keep: np.ndarray, start: int = 0
) -> sp.csr_matrix:
    """The rows `keep` of a CSR matrix, read without the other rows.

    `keep` holds row numbers from `start`. A member of a consolidated store starts
    at its row in the store.

    The plain call is `sparse_dataset(node)[start + keep]`. anndata then reads all of
    `indptr`, which in a consolidated store holds the rows of every member. For row
    positions, it makes one array position for each stored value of each kept row.
    zarr then selects the values by position. The CPU time and the memory of that
    read increase with the stored values.

    This function reads `indptr` from the first to the last kept row only. The stored
    values of adjacent kept rows are one span of `data` and of `indices`. Thus this
    function reads one slice of each array for each run of adjacent kept rows. Above
    `MAX_RUNS` runs, it selects the values by position, as the plain call does, which
    makes fewer requests. A matrix in another encoding, such as CSC, uses the plain
    call.
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
        return sp.csr_matrix(cast("sp.csr_matrix", sparse_dataset(node)[start + keep]))

    first, last = start + int(keep[0]), start + int(keep[-1])
    indptr = np.asarray(indptr_node[first : last + 2], dtype=np.int64)
    local = keep - keep[0]  # the kept rows, as positions in `indptr`
    begins = indptr[local]
    row_nnz = indptr[local + 1] - begins

    # A new run starts at each kept row that does not follow the kept row before it.
    breaks = np.flatnonzero(np.diff(keep) != 1) + 1
    if len(breaks) + 1 > MAX_RUNS:
        # One position for each stored value of each kept row, as the plain call makes.
        offsets = np.cumsum(row_nnz) - row_nnz
        positions = np.arange(int(row_nnz.sum())) + np.repeat(begins - offsets, row_nnz)
        kept_data, kept_indices = (np.asarray(a[positions]) for a in (data, indices))
    else:
        runs = np.split(local, breaks)
        spans = [slice(int(indptr[run[0]]), int(indptr[run[-1] + 1])) for run in runs]
        # A zarr slice read blocks the thread that calls it, so one thread reads one
        # span at a time. This pool keeps as many as RUN_THREADS spans of the store in
        # flight. `get_anndata` calls this function from one thread for each dataset
        # that it reads. Thus a query has no more than max_workers * RUN_THREADS span
        # reads in flight. On 103 stores, the pool decreases the wall time of a query
        # by 5 to 60 percent. It adds as much as 10 percent of CPU time.
        with ThreadPoolExecutor(min(RUN_THREADS, 2 * len(spans))) as pool:
            # `map` submits each read at once, so both arrays are in flight before the
            # wait.
            pending = [pool.map(array.__getitem__, spans) for array in (data, indices)]
            kept_data, kept_indices = (
                np.concatenate([np.asarray(part) for part in parts])
                for parts in pending
            )

    return sp.csr_matrix(
        (kept_data, kept_indices, np.r_[0, np.cumsum(row_nnz)]),
        shape=(len(keep), n_vars),
    )
