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

A catalog row can also name a compacted store, a second copy of several small
datasets in one store, with the rows of its dataset there (`compacted_store_uri`,
`compacted_row_start`, `compacted_row_stop`). The SDK then reads the dataset from
the compacted store, and opens each store once for all of its selected datasets.
Exclusion stays per dataset, so a compacted store opens only when one of its
datasets can match. The SDK reads only the rows of the selected datasets, and the
result is the result of the per-dataset stores. `open_census(..., compacted=False)`
reads the per-dataset stores, which stay the primary copy.
"""

import os
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
    COMPACTED,
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

    Each call to `get_anndata` or `count_cells` runs one query over them.
    """

    def __init__(self, sources: Sequence[PathLike], *, compacted: bool = True) -> None:
        self.sources = [str(s) for s in sources]
        # Read a dataset from the compacted store that its catalog row names, if any.
        self.compacted = compacted
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
        `max_workers` threads. The datasets of one compacted store thus read at
        the same time as well. The result is the same for every allowed value:
        the datasets keep their catalog order.
        """
        unknown = sorted(set(layers) - {"normalized"})
        if unknown:
            raise ValueError(f"no layer {unknown}; the one layer is 'normalized'")
        canonical = _catalog.canonicalize(obs_filter)
        rows = _catalog.exclude(self.catalog, organism, canonical)
        _check_composition(rows, ORGANISMS[organism]["version"])

        reads = _reads(rows, self.compacted)
        # Each store opens once, for all of its selected datasets.
        uris = list(dict.fromkeys(r.uri for r in reads))
        compacted = {r.uri for r in reads if r.rows is not None}
        # ponytail: serial by default; raise the default after the store-count benchmark.
        with ThreadPoolExecutor(max_workers=max_workers) as pool:
            opened = dict(
                zip(
                    uris,
                    pool.map(lambda u: _open(u, u in compacted), uris),
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
        # except the assayed genes of its own dataset. A compacted store holds
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


def open_census(sources: Sequence[PathLike], *, compacted: bool = True) -> Census:
    """Open catalogs.

    With `compacted` false, the SDK reads every dataset from its per-dataset store,
    and ignores the compacted stores that the catalog names.
    """
    return Census(sources, compacted=compacted)


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

    `rows` is the dataset's row range in a compacted store, or None for a
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
    compacted store, read once for all the datasets of the store.
    """

    uri: str
    group: zarr.Group
    var: pd.DataFrame
    n_obs: int
    categories: dict[str, pd.Index] = field(default_factory=dict)


def _reads(rows: pd.DataFrame, compacted: bool) -> list[_Read]:
    """The read of each selected dataset, in catalog order.

    A row that names a compacted store is read there, unless `compacted` is false.
    A catalog without the columns names none.
    """
    records = cast("list[dict[str, Any]]", rows.to_dict("records"))
    use = compacted and all(c in rows.columns for c in COMPACTED)
    reads = []
    for row in records:
        uri = row.get(COMPACTED[0]) if use else None
        if unrecorded(uri):
            reads.append(_Read(row, str(row["store_uri"]), None))
        else:
            span = (int(row[COMPACTED[1]]), int(row[COMPACTED[2]]))
            reads.append(_Read(row, str(uri), span))
    return reads


def _group(parent: zarr.Group, name: str, uri: str, kind: str) -> zarr.Group:
    """The child group `name`, or a TypeError that names the store and the node."""
    node = parent[name]
    if not isinstance(node, zarr.Group):
        raise TypeError(f"{uri}: {name} must be {kind}")
    return node


def _open(uri: str, compacted: bool) -> _Opened:
    """Open a store, and read what its datasets share: the var index, and categories.

    `var` holds only the index, the `feature_id`s, which is all that `_check_var`
    compares. `get_anndata` reads the var columns from one store. The number of
    rows comes from the attributes of `X`, which the open already read. The
    categories of a compacted store's `obs` columns serve each member's read.
    """
    group = open_store(uri)
    var_group = _group(group, "var", uri, "an encoded dataframe")
    index = read_elem(var_group[str(var_group.attrs["_index"])])
    var = pd.DataFrame(index=pd.Index(np.asarray(index).astype(str)))
    x = _group(group, "X", uri, "an encoded sparse matrix")
    opened = _Opened(uri, group, var, int(cast("list[int]", x.attrs["shape"])[0]))
    if not compacted:
        return opened
    _check_record(opened)
    obs = _group(group, "obs", uri, "an encoded dataframe")
    for name in cast("list[str]", obs.attrs["column-order"]):
        node = obs[str(name)]
        if isinstance(node, zarr.Group) and _encoding(node) == "categorical":
            categories = np.asarray(read_elem(node["categories"]))
            opened.categories[str(name)] = pd.Index(categories)
    return opened


def _encoding(node: zarr.Array[Any] | zarr.Group) -> str:
    return str(node.attrs.get("encoding-type", ""))


def _check_record(opened: _Opened) -> None:
    """A compacted store holds a compaction record, in its root attributes."""
    namespace = opened.group.attrs.get("cellxgene-census2")
    if not isinstance(namespace, dict) or "compaction" not in namespace:
        raise ValueError(f"{opened.uri}: no compaction record; not a compacted store")


def _member_rows(opened: _Opened, read: _Read) -> tuple[int, int]:
    """The rows of the dataset in a compacted store, which its record must list."""
    namespace = cast("dict[str, Any]", opened.group.attrs["cellxgene-census2"])
    listed = cast("list[dict[str, Any]]", namespace["compaction"].get("members", []))
    version = str(read.row["dataset_version_id"])
    spans = [
        (m.get("row_start"), m.get("row_stop"))
        for m in listed
        if m.get("dataset_version_id") == version
    ]
    if read.rows is None or spans != [read.rows]:
        raise ValueError(
            f"{opened.uri}: the compaction record does not list {version} at rows "
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
    per-dataset store, or the member's rows of a compacted store. `layers` names
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
    # the cell key in a compacted store), because the result sets the cell key. It
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
    _check_cell_key(obs, where)
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

    obs.index = pd.Index(dataset_id + "_" + obs["observation_joinid"].astype(str))
    return anndata.AnnData(X=x, obs=obs, var=opened.var, layers=computed)


# The obs columns of one dataset that the SDK reads at the same time.
OBS_THREADS = 8


def _obs_rows(
    opened: _Opened, group: zarr.Group, columns: list[str], start: int, stop: int
) -> pd.DataFrame:
    """Rows `start` to `stop` of the `obs` columns `columns`, without the other rows.

    A read of all the rows of a store reads each column whole with `read_elem`.
    A read of part of the rows gives each column the type that `read_elem` gives
    it, and a categorical column keeps only the categories that its rows use. A
    per-dataset store holds no unused category, because the writer writes none, so
    a member of a compacted store gets the categories of its per-dataset store.
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
    # One run of rows reads as a slice, which zarr reads faster than a selection.
    contiguous = rows[-1] - rows[0] + 1 == rows.size
    selection: Any = slice(int(rows[0]), int(rows[-1]) + 1) if contiguous else rows
    result: Any = pd.array(np.asarray(values.oindex[selection]), dtype=dtype)
    result[np.asarray(mask.oindex[selection])] = pd.NA
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
