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
result is the result of the per-dataset stores. `open_census(..., consolidated=False)`
reads the per-dataset stores, which stay the primary copy.
"""

import functools
import os
import threading
import warnings
from collections.abc import Sequence
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from typing import Any, cast

import anndata
import numpy as np
import pandas as pd
import scipy.sparse as sp
import zarr
from anndata.io import read_elem

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
from census2_sdk.obs_arrow import (
    OBS_THREADS,
    CategoryCache,
    ColumnData,
    encoding,
    evaluate,
    read_column,
    read_columns,
    take,
)
from census2_sdk.obs_arrow import read_rows as _read_column_rows
from census2_sdk.obs_columns import Obs
from census2_sdk.obs_filter import Expr, FilterError, FilterLike, is_const
from census2_sdk.obs_filter import columns as columns_of
from census2_sdk.obs_plan import Plan, ValueSets, resolve, suggestions, value_sets

__all__ = [
    "Census",
    "FilterError",
    "Organism",
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

    @functools.cached_property
    def homo_sapiens(self) -> "Organism":
        """The human datasets, with the typed filter builder in `.obs`."""
        return Organism(self, "homo_sapiens")

    @functools.cached_property
    def mus_musculus(self) -> "Organism":
        """The mouse datasets, with the typed filter builder in `.obs`."""
        return Organism(self, "mus_musculus")

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
        # Each store opens once, for all of its selected datasets.
        uris = list(dict.fromkeys(r.uri for r in reads))
        consolidated = {r.uri for r in reads if r.rows is not None}
        # ponytail: serial by default; raise the default after the store-count benchmark.
        with ThreadPoolExecutor(max_workers=max_workers) as pool:
            opened = dict(
                zip(
                    uris,
                    pool.map(lambda u: _open(u, u in consolidated), uris),
                    strict=True,
                )
            )
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

    `var` holds only the index. `n_obs` is the number of rows of the store.
    `categories` holds the categories of each categorical `obs` column of a
    consolidated store, read once for all the datasets of the store.
    """

    uri: str
    group: zarr.Group
    var: pd.DataFrame
    n_obs: int
    categories: dict[str, pd.Index] = field(default_factory=dict)
    # The categories of the builder path, read once for each column when needed.
    cache: CategoryCache | None = None


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


def _open(uri: str, consolidated: bool) -> _Opened:
    """Open a store, and read what its datasets share: the var index, and categories.

    `var` holds only the index, the `feature_id`s, which is all that `_check_var`
    compares. `get_anndata` reads the var columns from one store. The number of
    rows comes from the attributes of `X`, which the open already read. The
    categories of a consolidated store's `obs` columns serve each member's read.
    """
    group = open_store(uri)
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

    # Each column is a read of its own, so read them at the same time, as `_obs_rows`
    # reads the columns of a dataset.
    threads = max(1, min(len(categorical), OBS_THREADS))
    with ThreadPoolExecutor(max_workers=threads) as pool:
        opened.categories = dict(
            zip(categorical, pool.map(categories, categorical), strict=True)
        )
    return opened


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


def _read_rows(
    opened: _Opened,
    read: _Read,
    obs_filter: CanonicalFilter,
    layers: Sequence[str] = (),
) -> anndata.AnnData:
    """The selected cells of one dataset. `X` holds the counts on the reference genes.

    The read takes only the rows of the dataset from `obs` and `X`: all of a
    per-dataset store, or the member's rows of a consolidated store. `layers` names
    the layers to calculate from the counts.
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
    obs = _obs_rows(opened, obs_group, [c for c in columns if c != key], start, stop)
    _check_dataset_id(obs, dataset_id, where)
    keep = np.flatnonzero(_mask(obs, obs_filter))
    obs = obs.iloc[keep].copy()
    if key in columns:
        kept = _read_column_rows(obs_group[key], start + keep)
        obs.insert(columns.index(key), key, kept)
    return _part(opened, where, dataset_id, obs, keep, start, layers)


def _part(
    opened: "_Opened",
    where: str,
    dataset_id: str,
    obs: pd.DataFrame,
    keep: np.ndarray,
    start: int,
    layers: Sequence[str],
) -> anndata.AnnData:
    """The AnnData of the kept cells of one dataset: `obs`, the counts, and `layers`.

    `keep` holds the kept rows from `start`, and `obs` their rows. Both query paths
    end here, so their results agree.
    """
    group, uri, var = opened.group, opened.uri, opened.var
    _check_cell_key(obs, where)
    x = _rows(group["X"], keep, start)
    computed = {}
    if "normalized" in layers:
        assay = obs["assay_ontology_term_id"].astype(str)
        full_gene = assay.isin(_normalize.FULL_GENE_ASSAYS).to_numpy()
        length = None
        if full_gene.any():
            var_group = _group(group, "var", uri, "an encoded dataframe")
            length = np.asarray(read_elem(var_group["feature_length"]), dtype=np.int64)
        computed["normalized"] = _normalize.normalized(x, full_gene, length)

    obs.index = pd.Index(dataset_id + "_" + obs["observation_joinid"].astype(str))
    return anndata.AnnData(X=x, obs=obs, var=var, layers=computed)


def _obs_rows(
    opened: _Opened, group: zarr.Group, columns: list[str], start: int, stop: int
) -> pd.DataFrame:
    """Rows `start` to `stop` of the `obs` columns `columns`, without the other rows.

    A read of all the rows of a store reads each column whole with `read_elem`.
    A read of part of the rows gives each column the type that `read_elem` gives
    it, and a categorical column keeps only the categories that its rows use. A
    per-dataset store holds no unused category, because the writer writes none, so
    a member of a consolidated store gets the categories of its per-dataset store.
    """
    whole = (start, stop) == (0, opened.n_obs)

    def one(name: str) -> Any:
        if whole:
            return read_elem(group[name])
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
    reads whole, and then the rows are cut from it.
    """
    node = group[name]
    kind = _encoding(node)
    if isinstance(node, zarr.Array) and kind in ("array", "string-array"):
        return node[start:stop]
    if isinstance(node, zarr.Group) and kind == "categorical":
        codes = np.asarray(cast("zarr.Array[Any]", node["codes"])[start:stop])
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
        return _read_column_rows(node, np.arange(start, stop))
    # `read_elem` is typed for every element kind; here it is an array.
    return cast("Any", read_elem(node))[start:stop]


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


def _rows(
    node: zarr.Array[Any] | zarr.Group, keep: np.ndarray, start: int = 0
) -> sp.csr_matrix:
    """The selected rows of a CSR matrix, read without loading the rest.

    `keep` holds row numbers from `start`. The read takes `indptr` from the first
    to the last selected row, then the `data` and `indices` of the selected rows.
    When the selected rows follow one another, that is one slice. Otherwise it is
    one selection of every selected stored value, as anndata's row read makes.
    """
    if not isinstance(node, zarr.Group):
        raise TypeError(f"{node.store_path} must be an encoded sparse matrix")
    arrays = [node[name] for name in ("data", "indices", "indptr")]
    if not all(isinstance(a, zarr.Array) for a in arrays):
        raise TypeError(f"{node.store_path} must hold data, indices, and indptr arrays")
    data_array, indices_array, indptr_array = cast("list[zarr.Array[Any]]", arrays)
    n_var = int(cast("list[int]", node.attrs["shape"])[1])
    if not len(keep):
        return sp.csr_matrix((0, n_var), dtype=data_array.dtype)
    first, last = start + int(keep[0]), start + int(keep[-1])
    indptr = np.asarray(indptr_array[first : last + 2])
    local = keep - keep[0]
    begins, lengths = indptr[local], indptr[local + 1] - indptr[local]
    selection: slice | np.ndarray
    if len(keep) == last - first + 1:
        selection = slice(int(indptr[0]), int(indptr[-1]))
    else:
        offsets = np.cumsum(lengths) - lengths
        selection = np.arange(int(lengths.sum())) + np.repeat(begins - offsets, lengths)
    data = np.asarray(data_array[selection])
    indices = np.asarray(indices_array[selection])
    new_indptr = np.concatenate([[0], np.cumsum(lengths)])
    return sp.csr_matrix((data, indices, new_indptr), shape=(len(keep), n_var))


# The typed filter builder: one namespace for each organism.


class Organism:
    """The queries of one organism, with a typed filter: `census.homo_sapiens`.

    `obs` holds the typed columns of `census2_sdk.obs_columns`. Each query resolves
    its filter once against the catalog (`plan`), and then binds it to the `obs` of
    each store that can hold a match (`census2_sdk.obs_arrow`). A filter of another
    organism raises `FilterError`.
    """

    def __init__(self, census: "Census", name: str) -> None:
        if name not in ORGANISMS:
            raise ValueError(f"no organism {name!r}; use one of {sorted(ORGANISMS)}")
        self.census = census
        self.name = name
        self.obs = Obs(name)
        self._known: ValueSets | None = None
        # The counts of the last query: catalog candidates, stores opened, reads.
        self.last_query: dict[str, int] = {}

    def __repr__(self) -> str:
        return f"<Organism {self.name} of {len(self.census.sources)} sources>"

    def plan(self, obs_filter: FilterLike | None = None) -> Plan:
        """Resolve `obs_filter` against the catalog: see `census2_sdk.obs_plan`."""
        if self._known is None:
            catalog = self.census.catalog
            self._known = value_sets(catalog[catalog["organism"] == self.name])
        return resolve(obs_filter, self.name, self.census.catalog, self._known)

    def get_anndata(
        self,
        obs_filter: FilterLike | None = None,
        *,
        layers: Sequence[str] = (),
        max_workers: int = 1,
    ) -> anndata.AnnData | None:
        """The selected cells as one AnnData, as `Census.get_anndata` gives them.

        For a filter that a mapping filter can express, the result equals the result
        of `Census.get_anndata`, value for value, except that each categorical column
        holds only the labels that the result uses. Return None when no cell matches.
        """
        unknown = sorted(set(layers) - {"normalized"})
        if unknown:
            raise ValueError(f"no layer {unknown}; the one layer is 'normalized'")
        plan = self.plan(obs_filter)
        rows = pd.DataFrame([c.row for c in plan.candidates])
        if len(rows):
            _check_composition(rows, ORGANISMS[self.name]["version"])
        scan = self._scan(plan, "anndata", max_workers, layers=layers)
        parts = [p for p in scan.values if isinstance(p, anndata.AnnData)]
        if not parts or not sum(p.n_obs for p in parts):
            return None
        _check_var(parts, [str(p.uns.pop("_uri")) for p in parts])
        result = anndata.concat(parts, join="inner")
        if not result.obs_names.is_unique:
            raise SourceConflict("the cell key is not unique across the result")
        first = scan.opened[scan.uris[0]]
        var = cast("pd.DataFrame", read_elem(first.group["var"]))
        result.var = var.drop(columns=MEASURED_COLUMN, errors="ignore")
        result.uns["census2"] = {
            "sources": self.census.sources,
            "datasets_in_catalog": plan.catalog_rows,
            "datasets_read": len(parts),
            "stores_opened": len(scan.uris),
        }
        return result

    def get_obs(
        self,
        obs_filter: FilterLike | None = None,
        *,
        columns: Sequence[str] | None = None,
        max_workers: int = 1,
    ) -> pd.DataFrame:
        """The `obs` of the selected cells, without `X`, indexed by the cell key.

        `columns` names the columns of the result, all of them by default. A
        categorical column holds the labels that the result uses.
        """
        plan = self.plan(obs_filter)
        scan = self._scan(plan, "obs", max_workers, columns=columns)
        frames = [f for f in scan.values if isinstance(f, pd.DataFrame)]
        return _concat_frames(frames)

    def count(
        self, obs_filter: FilterLike | None = None, *, max_workers: int = 1
    ) -> int:
        """The number of selected cells.

        A dataset whose residual filter is True counts from its catalog row, with no
        store open. Another dataset reads only the columns of its residual filter.
        """
        plan = self.plan(obs_filter)
        return sum(int(n) for n in self._scan(plan, "count", max_workers).values)

    def count_by(
        self,
        by: str | Sequence[str],
        obs_filter: FilterLike | None = None,
        *,
        max_workers: int = 1,
    ) -> pd.Series:
        """The number of selected cells for each combination of the `by` columns."""
        names = [by] if isinstance(by, str) else list(by)
        plan = self.plan(obs_filter)
        scan = self._scan(plan, "count_by", max_workers, columns=names)
        parts = [p for p in scan.values if isinstance(p, pd.DataFrame)]
        if not parts:
            return pd.Series(dtype=np.int64, name="n_cells")
        frame = pd.concat(parts, ignore_index=True)
        return frame.groupby(names, sort=True)["n_cells"].sum().rename("n_cells")

    def _scan(
        self,
        plan: Plan,
        mode: str,
        max_workers: int,
        *,
        layers: Sequence[str] = (),
        columns: Sequence[str] | None = None,
    ) -> "_Scan":
        """Read each candidate dataset by `mode`: anndata, obs, count, or count_by."""
        units = _units(plan, self.census.consolidated)
        values: list[Any] = [None] * len(plan.candidates)
        need_open = []
        for unit in units:
            if mode == "count" and all(
                is_const(m.residual, True) for m in unit.members
            ):
                for m in unit.members:
                    whole = int(m.read.row["n_obs"])
                    values[m.index] = m.stop - m.start if m.read.rows else whole
            else:
                need_open.append(unit)
        uris = list(dict.fromkeys(u.uri for u in need_open))
        consolidated = {u.uri for u in need_open if u.consolidated}
        seen: set[str] = set()
        missing: dict[str, list[str]] = {}
        lock = threading.Lock()

        def scan(unit: _Unit, store: _Opened) -> None:
            present, lacking, out = _scan_unit(unit, store, mode, layers, columns)
            with lock:
                seen.update(present)
                for name in lacking:
                    missing.setdefault(name, []).append(unit.uri)
            for index, value in out:
                values[index] = value

        need_var = mode == "anndata"
        with ThreadPoolExecutor(max_workers=max_workers) as pool:
            opened = dict(
                zip(
                    uris,
                    pool.map(lambda u: _store(u, u in consolidated, need_var), uris),
                    strict=True,
                )
            )
            list(pool.map(lambda u: scan(u, opened[u.uri]), need_open))
        self.last_query = {
            "candidates": len(plan.candidates),
            "stores_opened": len(uris),
            "reads": len(need_open),
        }
        _check_missing(missing, seen, len(need_open))
        return _Scan(values, uris, opened)


@dataclass(frozen=True)
class _Scan:
    """The values of a scan, one for each candidate, and the stores that it opened."""

    values: list[Any]
    uris: list[str]
    opened: dict[str, "_Opened"]


def _check_missing(missing: dict[str, list[str]], seen: set[str], reads: int) -> None:
    """Refuse a filter column that no read store holds, and warn for a column that
    some stores lack. Such a store gives no cells, as in `_mask`.
    """
    for name, uris in sorted(missing.items()):
        if name not in seen:
            near = suggestions(name, frozenset(seen))
            hint = f"; did you mean {' or '.join(map(repr, near))}?" if near else ""
            raise FilterError(f"no store holds the obs column {name!r}{hint}")
        stores = sorted(set(uris))
        shown = ", ".join(stores[:3]) + (", ..." if len(stores) > 3 else "")
        warnings.warn(
            f"{len(stores)} of the stores read have no obs column {name!r}, so they "
            f"give no cells: {shown}",
            UserWarning,
            stacklevel=4,
        )


@dataclass(frozen=True)
class _Member:
    """One candidate dataset of a read: its rows in the store, and its filter."""

    index: int
    read: "_Read"
    residual: Expr
    # The rows of a consolidated member; 0 and 0 for a per-dataset store.
    start: int
    stop: int


@dataclass(frozen=True)
class _Unit:
    """One read: a per-dataset store, or adjacent members of a consolidated store."""

    uri: str
    consolidated: bool
    members: tuple[_Member, ...]


def _units(plan: Plan, consolidated: bool) -> list[_Unit]:
    """The reads of the candidates. Adjacent members of a consolidated store read once."""
    if not plan.candidates:
        return []
    rows = pd.DataFrame([c.row for c in plan.candidates])
    reads = _reads(rows, consolidated)
    by_uri: dict[str, list[_Member]] = {}
    units: list[_Unit] = []
    for index, (candidate, read) in enumerate(zip(plan.candidates, reads, strict=True)):
        start, stop = read.rows if read.rows is not None else (0, 0)
        member = _Member(index, read, candidate.residual, start, stop)
        by_uri.setdefault(read.uri, []).append(member)
    for uri, listed in by_uri.items():
        if listed[0].read.rows is None:
            units.extend(_Unit(uri, False, (m,)) for m in listed)
            continue
        listed.sort(key=lambda m: m.start)
        run = [listed[0]]
        for m in listed[1:]:
            if m.start == run[-1].stop:
                run.append(m)
            else:
                units.append(_Unit(uri, True, tuple(run)))
                run = [m]
        units.append(_Unit(uri, True, tuple(run)))
    return units


def _store(uri: str, consolidated: bool, need_var: bool) -> "_Opened":
    """Open a store with the lean open, and read the var index when `need_var`."""
    group = open_store(uri, lean=True)
    var = pd.DataFrame()
    if need_var:
        var_group = _group(group, "var", uri, "an encoded dataframe")
        index = read_elem(var_group[str(var_group.attrs["_index"])])
        var = pd.DataFrame(index=pd.Index(np.asarray(index).astype(str)))
    x = _group(group, "X", uri, "an encoded sparse matrix")
    opened = _Opened(uri, group, var, int(cast("list[int]", x.attrs["shape"])[0]))
    if consolidated:
        _check_record(opened)
    opened.cache = CategoryCache(_group(group, "obs", uri, "an encoded dataframe"))
    return opened


def _scan_unit(
    unit: _Unit,
    opened: "_Opened",
    mode: str,
    layers: Sequence[str],
    columns: Sequence[str] | None,
) -> tuple[set[str], set[str], list[tuple[int, Any]]]:
    """Bind and evaluate the filter of each member of one read, then read by `mode`.

    Return the columns of the store, the filter columns that it lacks, and one value
    for each member.
    """
    obs = _group(opened.group, "obs", opened.uri, "an encoded dataframe")
    order = [str(c) for c in cast("list[str]", obs.attrs["column-order"])]
    present = set(order)
    spans = [
        (0, opened.n_obs) if m.read.rows is None else _member_rows(opened, m.read)
        for m in unit.members
    ]
    start, stop = spans[0][0], spans[-1][1]
    wanted: set[str] = set()
    for m in unit.members:
        wanted |= columns_of(m.residual)
    lacking = wanted - present
    extra = set(columns or ()) & present if mode == "count_by" else set()
    cache = opened.cache
    data = read_columns(obs, sorted((wanted & present) | extra), start, stop, cache)
    out: list[tuple[int, Any]] = []
    for m, (lo, hi) in zip(unit.members, spans, strict=True):
        dataset_id = str(m.read.row["dataset_id"])
        where = opened.uri if m.read.rows is None else f"{opened.uri}[{dataset_id}]"
        names = columns_of(m.residual)
        if names & lacking:
            keep = np.zeros(0, dtype=np.int64)  # the rule of `_mask`
        else:
            local = {k: data[k].part(lo - start, hi - start) for k in names}
            keep = evaluate(m.residual, local, hi - lo)
        if mode == "count":
            out.append((m.index, len(keep)))
        elif mode == "count_by":
            span = (lo - start, hi - start)
            out.append((m.index, _group_counts(data, columns or (), span, keep)))
        elif mode == "obs":
            names_out = order if columns is None else list(columns)
            frame = _kept_obs(obs, names_out, data, lo, start, keep, cache, where)
            if "dataset_id" in frame.columns and len(frame):
                _check_dataset_id(frame, dataset_id, where)
            frame.index = pd.Index(_cell_keys(obs, dataset_id, lo + keep))
            out.append((m.index, frame))
        else:
            frame = _kept_obs(obs, order, data, lo, start, keep, cache, where)
            if len(frame):
                _check_dataset_id(frame, dataset_id, where)
            part = _part(opened, where, dataset_id, frame, keep, lo, layers)
            part.uns["_uri"] = opened.uri
            out.append((m.index, part))
    return present, lacking, out


def _cell_keys(obs: zarr.Group, dataset_id: str, rows: np.ndarray) -> np.ndarray:
    """The cell keys of the rows `rows`: `<dataset_id>_<observation_joinid>`."""
    joinid = pd.Series(_read_column_rows(obs["observation_joinid"], rows))
    return np.asarray(dataset_id + "_" + joinid.astype(str), dtype=object)


def _kept_obs(
    obs: zarr.Group,
    names: Sequence[str],
    data: dict[str, ColumnData],
    lo: int,
    start: int,
    keep: np.ndarray,
    cache: CategoryCache | None,
    where: str,
) -> pd.DataFrame:
    """The columns `names` of the kept rows of one member, as `read_elem` types them.

    A column that the filter read is cut from that read. Another categorical or numeric
    column reads the rows from the first to the last kept row; a column of another
    encoding, such as the cell key, reads the kept rows only.
    """
    order = [str(c) for c in cast("list[str]", obs.attrs["column-order"])]
    unknown = [n for n in names if n not in order]
    if unknown:
        raise FilterError(f"{where}: obs has no column {unknown}")
    rows = lo + keep  # store rows
    first, last = (int(rows[0]), int(rows[-1]) + 1) if len(rows) else (lo, lo)

    def one(name: str) -> Any:
        node = obs[name]
        kind = encoding(node)[0]
        if kind not in ("categorical", "array") or name == "observation_joinid":
            return _read_column_rows(node, rows)
        found = data.get(name)
        if found is not None:
            return take(found, rows - start)
        cats = None
        if cache is not None and kind == "categorical":
            cats = functools.partial(cache, name)
        return take(read_column(node, first, last, cats), rows - first)

    with ThreadPoolExecutor(max_workers=max(1, min(len(names), OBS_THREADS))) as pool:
        values = list(pool.map(one, names))
    return pd.DataFrame(dict(zip(names, values, strict=True)), columns=list(names))


def _group_counts(
    data: dict[str, ColumnData],
    by: Sequence[str],
    span: tuple[int, int],
    keep: np.ndarray,
) -> pd.DataFrame:
    """The kept cells of one member, counted by the `by` columns."""
    missing = [b for b in by if b not in data]
    if missing:
        raise FilterError(f"obs has no column {missing} to count by")
    frame = pd.DataFrame({b: take(data[b].part(*span), keep) for b in by})
    counts = frame.groupby(list(by), observed=True, sort=False).size()
    out = counts.rename("n_cells").reset_index()
    for b in by:
        if isinstance(out[b].dtype, pd.CategoricalDtype):
            out[b] = out[b].astype(str)
    return out


def _concat_frames(frames: list[pd.DataFrame]) -> pd.DataFrame:
    """Concatenate result frames; a categorical column unions its categories."""
    if not frames:
        return pd.DataFrame()
    if len(frames) == 1:
        return frames[0]
    columns = list(frames[0].columns)
    out: dict[str, Any] = {}
    for name in columns:
        parts = [f[name] for f in frames]
        if all(isinstance(p.dtype, pd.CategoricalDtype) for p in parts):
            out[name] = pd.api.types.union_categoricals(
                [p.array for p in parts], sort_categories=True
            )
        else:
            out[name] = pd.concat(parts, ignore_index=True).array
    index = pd.Index(np.concatenate([f.index.to_numpy() for f in frames]))
    return pd.DataFrame(out, columns=columns, index=index)
