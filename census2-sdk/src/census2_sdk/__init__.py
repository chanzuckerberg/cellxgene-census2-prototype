"""Query Census 2 stores through a catalog and return one AnnData.

The SDK runs on the user's machine. It reads the catalog, excludes the datasets that
cannot match, opens the remaining stores, and concatenates the selected cells by gene
position (design 4.1 to 4.4).

    import census2_sdk

    census = census2_sdk.open_census(sources=["catalog.parquet"])
    adata = census.get_anndata("homo_sapiens", obs_filter={"tissue": "lung"})

`obs_filter` maps a catalog column to one value or to a list of values. The SDK selects
a cell when every column holds one of its values. The columns are the 18 indexed
columns of design 2.2 and `dataset_id`. A filter grammar is open (design section 8).

The result carries counts in `X`, the pinned gene reference as `var`, the `obs` columns
that all selected stores hold, and `uns["census2"]` with the sources and the store
counts. The cell key `(dataset_id, observation_joinid)` is in `obs`, and `obs_names`
is `<dataset_id>_<observation_joinid>`. When no cell matches, `get_anndata` returns
None.

The writer projects each store onto the reference. Thus every selected store has
the same `var`. The SDK checks this, and concatenates by position.

`count_cells` counts the selected cells and opens no store. It sums the cell count
table that the catalog builder writes next to each catalog (tracker 5.9). Its filter
takes the columns of `catalog.COUNTED` and `dataset_id`.

    n = census.count_cells("homo_sapiens", {"cell_type": "neutrophil"})

Iteration 1 scaffolding, marked `ponytail:` in the code: the SDK reads `dataset_id`
from the catalog because the stores do not hold the column.
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
from census2_sdk.catalog import Filter, Normalized, members, open_store, unrecorded
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
    """A set of sources, opened. `get_anndata` and `count_cells` run one query over them."""

    def __init__(self, sources: Sequence[PathLike]) -> None:
        self.sources = [str(s) for s in sources]
        self._counts: pd.DataFrame | None = None
        frames = []
        for source in self.sources:
            frame = _catalog.load(source)
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

        `max_workers` is the largest number of stores that the SDK reads at the
        same time, and it must be a positive integer. The result is the same for
        every allowed value: the stores keep their catalog order.
        """
        normalized = _catalog.normalize(obs_filter)
        rows = _catalog.exclude(self.catalog, organism, normalized)
        _check_composition(rows, ORGANISMS[organism]["version"])

        records = cast("list[dict[str, Any]]", rows.to_dict("records"))
        # ponytail: serial by default; raise the default from the 6.4 results.
        with ThreadPoolExecutor(max_workers=max_workers) as pool:
            parts = list(pool.map(lambda row: _read_store(row, normalized), records))
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
        # version and one index, so the first store supplies the var columns.
        result.var = cast("pd.DataFrame", read_elem(open_store(uris[0])["var"]))
        result.uns["census2"] = {
            "sources": self.sources,
            "stores_in_catalog": int((self.catalog["organism"] == organism).sum()),
            "stores_opened": len(rows),
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

        The count comes from the cell count tables, and no store opens. A cell
        matches by the rule of `get_anndata`, and the catalog excludes the same
        datasets. When `get_anndata` returns a result, the count thus equals its
        `n_obs`. The count does not apply the checks before concatenation: a
        count joins no matrices. The filter takes the columns of
        `catalog.COUNTED` and `dataset_id`; a filter on another column raises
        ValueError. With `primary_only`, only the cells with `is_primary_data`
        count.
        """
        normalized = _catalog.normalize(obs_filter)
        uncounted = sorted(set(normalized) - {"dataset_id", *_catalog.COUNTED})
        if uncounted:
            raise ValueError(
                f"the cell count table holds no {uncounted}; use get_anndata"
            )
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
        """The cell count tables of all sources, read once.

        Refuse a table whose total for a dataset differs from the catalog
        `n_obs`. A table from another build gives a wrong count and no error.
        """
        if self._counts is None:
            frames = []
            for source in self.sources:
                uri = _catalog.cell_counts_uri(source)
                try:
                    frames.append(_catalog.load(uri))
                except FileNotFoundError as e:
                    raise FileNotFoundError(
                        f"{source}: no cell count table at {uri}; "
                        "build the catalog again with tools/build_catalog.py"
                    ) from e
            counts = pd.concat(frames, ignore_index=True)
            totals = counts.groupby("dataset_version_id")["n_cells"].sum()
            wrong = sorted(
                str(v)
                for v, n in zip(
                    self.catalog["dataset_version_id"],
                    self.catalog["n_obs"],
                    strict=True,
                )
                if int(totals.get(v, 0)) != int(n)
            )
            if wrong:
                raise SourceConflict(
                    f"the cell count table differs from the catalog n_obs for {wrong}"
                )
            self._counts = counts
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


def _read_store(row: dict[str, Any], obs_filter: Normalized) -> anndata.AnnData:
    """The selected cells of one store. `X` holds the counts on the reference genes.

    `var` holds only the index, the `feature_id`s, which is all that `_check_var`
    compares. `get_anndata` reads the var columns from one store.
    """
    group = open_store(row["store_uri"])
    obs = cast("pd.DataFrame", read_elem(group["obs"]))
    _check_cell_key(obs, str(row["store_uri"]))
    keep = np.flatnonzero(_mask(obs, obs_filter))
    obs = obs.iloc[keep].copy()
    var_group = group["var"]
    assert isinstance(var_group, zarr.Group), "var must be an encoded dataframe"
    index = read_elem(var_group[str(var_group.attrs["_index"])])
    var = pd.DataFrame(index=pd.Index(np.asarray(index).astype(str)))
    x = _rows(group["X"], keep)

    # ponytail: the stores hold no dataset_id column; it comes from the writer
    # once disposition row 13 lands in harmonize.
    dataset_id = str(row["dataset_id"])
    obs.insert(
        0,
        "dataset_id",
        pd.Series(pd.Categorical([dataset_id] * len(obs)), index=obs.index),
    )
    obs.index = pd.Index(dataset_id + "_" + obs["observation_joinid"].astype(str))
    return anndata.AnnData(X=x, obs=obs, var=var)


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


def _rows(node: zarr.Array[Any] | zarr.Group, keep: np.ndarray) -> sp.csr_matrix:
    """The selected rows of a CSR matrix, read without loading the rest."""
    assert isinstance(node, zarr.Group), "X must be an encoded sparse matrix"
    dataset = sparse_dataset(node)
    if not len(keep):
        return sp.csr_matrix((0, dataset.shape[1]), dtype=dataset.dtype)
    return sp.csr_matrix(cast("sp.csr_matrix", dataset[keep]))
