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
number of catalog rows for the organism, and the number of stores opened. The cell
key `(dataset_id, observation_joinid)` is in `obs`, and `obs_names` is
`<dataset_id>_<observation_joinid>`. When no cell matches, `get_anndata` returns None.

The writer projects each store onto the reference. Thus every selected store has
the same `var`. The SDK checks that every store has the same gene identifiers in the
same order, and concatenates by position.

`count_cells` counts the selected cells and opens no store. It sums the cell counts
that each catalog row holds (tracker 5.9). Its filter takes the columns of
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
from census2_sdk.catalog import (
    COMPACTED,
    Filter,
    Normalized,
    members,
    open_store,
    unrecorded,
)
from census2_sdk.gene_reference import ORGANISMS

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
        max_workers: int = 1,
    ) -> anndata.AnnData | None:
        """Return the selected cells of `organism` as one AnnData.

        Return None when no cell matches, whether or not the catalog excluded
        every store.

        `max_workers` is the largest number of reads that the SDK runs at the
        same time, and it must be a positive integer. The SDK first opens each
        store, then reads each selected dataset, and both steps use up to
        `max_workers` threads. The datasets of one compacted store thus read at
        the same time as well. The result is the same for every allowed value:
        the datasets keep their catalog order.
        """
        normalized = _catalog.normalize(obs_filter)
        rows = _catalog.exclude(self.catalog, organism, normalized)
        _check_composition(rows, ORGANISMS[organism]["version"])

        reads = _reads(rows, self.compacted)
        # Each store opens once, for all of its selected datasets.
        uris = list(dict.fromkeys(r.uri for r in reads))
        compacted = {r.uri for r in reads if r.rows is not None}
        # ponytail: serial by default; raise the default from the 6.4 results.
        with ThreadPoolExecutor(max_workers=max_workers) as pool:
            opened = dict(
                zip(
                    uris,
                    pool.map(lambda u: _open(u, u in compacted), uris),
                    strict=True,
                )
            )
            parts = list(
                pool.map(lambda r: _read_rows(opened[r.uri], r, normalized), reads)
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
        # version and one index, so the first store supplies the var columns.
        result.var = cast("pd.DataFrame", read_elem(opened[uris[0]].group["var"]))
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
        normalized = _catalog.normalize(obs_filter)
        uncounted = sorted(set(normalized) - {"dataset_id", *_catalog.COUNTED})
        if uncounted:
            raise ValueError(f"the catalog counts no {uncounted}; use get_anndata")
        counts = self._cell_counts()
        candidates = _catalog.exclude(self.catalog, organism, normalized)
        rows = counts[
            counts["dataset_version_id"].isin(candidates["dataset_version_id"])
        ]
        if primary_only:
            rows = rows[rows["is_primary_data"].astype(bool)]
        keep = np.ones(len(rows), dtype=bool)
        for column, wanted in normalized.items():
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


def _mask(obs: pd.DataFrame, obs_filter: Normalized) -> np.ndarray:
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
    `categories` holds the categories of each categorical `obs` column, read once
    for all the datasets of the store.
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


def _open(uri: str, compacted: bool) -> _Opened:
    """Open a store, and read what its datasets share: the var index and categories.

    `var` holds only the index, the `feature_id`s, which is all that `_check_var`
    compares. `get_anndata` reads the var columns from one store. The number of
    rows comes from the attributes of `X`, which the open already read.
    """
    group = open_store(uri)
    var_group = group["var"]
    assert isinstance(var_group, zarr.Group), "var must be an encoded dataframe"
    index = read_elem(var_group[str(var_group.attrs["_index"])])
    var = pd.DataFrame(index=pd.Index(np.asarray(index).astype(str)))
    x = group["X"]
    assert isinstance(x, zarr.Group), "X must be an encoded sparse matrix"
    opened = _Opened(uri, group, var, int(cast("list[int]", x.attrs["shape"])[0]))
    if compacted:
        _check_record(opened)
    obs = group["obs"]
    assert isinstance(obs, zarr.Group), "obs must be an encoded dataframe"
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
    members = cast("list[dict[str, Any]]", namespace["compaction"].get("members", []))
    version = str(read.row["dataset_version_id"])
    spans = [
        (m.get("row_start"), m.get("row_stop"))
        for m in members
        if m.get("dataset_version_id") == version
    ]
    if spans != [read.rows]:
        raise ValueError(
            f"{opened.uri}: the compaction record does not list {version} at rows "
            f"{read.rows}"
        )
    assert read.rows is not None
    return read.rows


def _read_rows(opened: _Opened, read: _Read, obs_filter: Normalized) -> anndata.AnnData:
    """The selected cells of one dataset. `X` holds the counts on the reference genes.

    The read takes only the rows of the dataset from `obs` and `X`: all of a
    per-dataset store, or the member's rows of a compacted store.
    """
    dataset_id = str(read.row["dataset_id"])
    if read.rows is None:
        where, (start, stop) = opened.uri, (0, opened.n_obs)
    else:
        where, (start, stop) = f"{opened.uri}[{dataset_id}]", _member_rows(opened, read)
    obs = _obs_rows(opened, start, stop)
    _check_cell_key(obs, where)
    _check_dataset_id(obs, dataset_id, where)
    keep = np.flatnonzero(_mask(obs, obs_filter))
    obs = obs.iloc[keep].copy()
    x = opened.group["X"]
    assert isinstance(x, zarr.Group), "X must be an encoded sparse matrix"
    matrix = _x_rows(x, start, keep)

    obs.index = pd.Index(dataset_id + "_" + obs["observation_joinid"].astype(str))
    return anndata.AnnData(X=matrix, obs=obs, var=opened.var)


def _obs_rows(opened: _Opened, start: int, stop: int) -> pd.DataFrame:
    """Rows `start` to `stop` of a store's `obs`, read without the other rows.

    Each column has the type that `read_elem` gives it. A categorical column keeps
    only the categories that its rows use. In a per-dataset store that is every
    category, because the writer writes no unused one. In a compacted store it
    gives the categories of the per-dataset store. The index is not read, because
    the caller sets the cell key.
    """
    group = opened.group["obs"]
    assert isinstance(group, zarr.Group), "obs must be an encoded dataframe"
    columns = [str(c) for c in cast("list[str]", group.attrs["column-order"])]
    data = {name: _column_rows(opened, group, name, start, stop) for name in columns}
    return pd.DataFrame(data, columns=columns)


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
        values = np.asarray(cast("zarr.Array[Any]", node["values"])[start:stop])
        mask = np.asarray(cast("zarr.Array[Any]", node["mask"])[start:stop])
        # As anndata reads this encoding.
        na_value: Any = np.nan if node.attrs.get("na-value", "NA") == "NaN" else pd.NA
        dtype = pd.StringDtype(na_value=na_value)
        strings = pd.array(
            values.astype(np.dtypes.StringDType(na_object=na_value)), dtype=dtype
        )
        strings[mask] = pd.NA
        return strings
    # `read_elem` is typed for every element kind; here it is an array.
    return cast("Any", read_elem(node))[start:stop]


def _x_rows(x: zarr.Group, start: int, keep: np.ndarray) -> sp.csr_matrix:
    """The selected rows of a CSR matrix, read without loading the rest.

    `keep` holds row numbers from `start`. The read takes `indptr` from the first
    to the last selected row, then the `data` and `indices` of the selected rows.
    When the selected rows follow one another, that is one slice.
    """
    data_array = cast("zarr.Array[Any]", x["data"])
    indices_array = cast("zarr.Array[Any]", x["indices"])
    n_var = int(cast("list[int]", x.attrs["shape"])[1])
    if not len(keep):
        return sp.csr_matrix((0, n_var), dtype=data_array.dtype)
    first, last = start + int(keep[0]), start + int(keep[-1])
    indptr = np.asarray(cast("zarr.Array[Any]", x["indptr"])[first : last + 2])
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
    """Design 1.6: `observation_joinid` is present, non-null, and unique in a store."""
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
