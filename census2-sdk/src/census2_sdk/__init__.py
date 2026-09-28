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
is `<dataset_id>_<observation_joinid>`. When the catalog excludes every store, the
SDK opens no store, and `var` has no genes.

The writer projects each store onto the reference. Thus every selected store has
the same `var`. The SDK checks this, and concatenates by position.

Iteration 1 scaffolding, marked `ponytail:` in the code: the SDK reads `dataset_id`
from the catalog because the stores do not hold the column.
"""

import os
from collections.abc import Sequence
from typing import Any, cast

import anndata
import numpy as np
import pandas as pd
import scipy.sparse as sp
import zarr
from anndata.io import read_elem, sparse_dataset

from census2_sdk import catalog as _catalog
from census2_sdk.catalog import Filter, Normalized, members, open_store
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
    """A set of sources, opened. `get_anndata` runs one query over them."""

    def __init__(self, sources: Sequence[PathLike]) -> None:
        self.sources = [str(s) for s in sources]
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
        self, organism: str, obs_filter: Filter | None = None
    ) -> anndata.AnnData:
        """Return the selected cells of `organism` as one AnnData."""
        normalized = _catalog.normalize(obs_filter)
        rows = _catalog.exclude(self.catalog, organism, normalized)
        _check_composition(rows, ORGANISMS[organism]["version"])

        parts = [
            _read_store(row, normalized)
            for row in cast("list[dict[str, Any]]", rows.to_dict("records"))
        ]
        _check_var(parts, list(rows["store_uri"].astype(str)))
        if parts:
            # Zero-row parts stay in, so the obs schema is the columns that all
            # opened stores hold even when no cell matched.
            result = anndata.concat(parts, join="inner", merge="same")
            if not result.obs_names.is_unique:
                raise SourceConflict("the cell key is not unique across the result")
        else:
            empty = pd.DataFrame(
                {
                    "dataset_id": pd.Series([], dtype=str),
                    "observation_joinid": pd.Series([], dtype=str),
                }
            )
            result = anndata.AnnData(
                X=sp.csr_matrix((0, 0), dtype=np.float32), obs=empty
            )
        result.uns["census2"] = {
            "sources": self.sources,
            "stores_in_catalog": int((self.catalog["organism"] == organism).sum()),
            "stores_opened": len(rows),
        }
        return result


def open_census(sources: Sequence[PathLike]) -> Census:
    """Open catalogs."""
    return Census(sources)


def _check_composition(rows: pd.DataFrame, pinned: str) -> None:
    """The checks before concatenation. Name the conflicting values.

    Every row must record the `gene_reference_version` that this SDK pins for the
    organism. The SDK does not remap a store onto another reference, and a row
    that records none is from a store that the writer did not project.
    """
    if not len(rows):
        return
    recorded = rows["gene_reference_version"]
    if bool(recorded.isna().any()):
        missing = sorted(rows.loc[recorded.isna(), "dataset_id"].astype(str))
        raise SourceConflict(f"no gene_reference_version recorded for {missing}")
    versions = sorted(set(recorded.astype(str)))
    if versions != [pinned]:
        raise SourceConflict(
            f"gene_reference_version {versions} differs from the pinned {pinned!r}"
        )
    spatial = sorted(set(rows["spatial"].astype(bool)))
    if len(spatial) > 1:
        raise SourceConflict("spatial and non-spatial stores in one result")


def _mask(obs: pd.DataFrame, obs_filter: Normalized) -> np.ndarray:
    """The cells of one store that match `obs_filter`, by the same rule as the catalog."""
    mask = np.ones(len(obs), dtype=bool)
    for column, wanted in obs_filter.items():
        if column == "dataset_id":
            continue  # the catalog answered it; the store holds the dataset whole
        if column not in obs.columns:
            return np.zeros(len(obs), dtype=bool)  # the catalog kept an unproven store
        # A null matches no value, and it is not the string "nan".
        hit = [
            not pd.isna(value) and bool(wanted & set(members(str(value), column)))
            for value in obs[column]
        ]
        mask &= np.asarray(hit, dtype=bool)
    return mask


def _read_store(row: dict[str, Any], obs_filter: Normalized) -> anndata.AnnData:
    """The selected cells of one store. `X` holds the counts on the reference genes."""
    group = open_store(row["store_uri"])
    obs = cast("pd.DataFrame", read_elem(group["obs"]))
    _check_cell_key(obs, str(row["store_uri"]))
    keep = np.flatnonzero(_mask(obs, obs_filter))
    obs = obs.iloc[keep].copy()
    var = cast("pd.DataFrame", read_elem(group["var"]))
    x = _rows(group["X"], keep)

    # ponytail: the stores hold no dataset_id column until the writer stamps it
    # into obs.
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
