"""Compaction: write a second copy of several per-dataset stores as one store.

A compacted store holds the cells of several datasets. It is a materialization:
the per-dataset stores stay the primary copy, and stay available. A catalog row
can then tell the SDK to read a dataset's rows from the compacted store, so a
query opens one store for several datasets. `compact` writes a compacted store,
and `verify` compares one with its per-dataset stores.

    record = compaction.compact(
        ["s3://bucket/stores/v1.zarr", "s3://bucket/stores/v2.zarr"],
        "s3://bucket/compacted/c-0123456789abcdef.zarr",
    )

Layout
------
A compacted store is an AnnData-Zarr store like a per-dataset store, so a plain
anndata reader can open it. Its members follow one another along `obs`, each in
the row order of its per-dataset store. The compaction record in the root
attributes lists the members with their row ranges (see
`store_format.CompactionRecord`). A reader of one member reads only that row
range of `obs` and `X`. These differences from a per-dataset store make that
read efficient and the store usable on its own:

- Each `obs` array and `X/indptr` has fixed chunks of `OBS_CHUNK_ROWS` rows in
  shards of `OBS_SHARD_ROWS` rows. A per-dataset store keeps anndata's
  default, which is often one chunk for a column, because a reader reads its
  `obs` whole.
- `X/data` and `X/indices` have fixed chunks of `X_CHUNK_VALUES` values in shards
  of `X_SHARD_VALUES` values.
- Each categorical column holds the union of the members' categories: the
  observed values, sorted, so no category is unused. A reader of one member
  removes the categories that its rows do not use, and so gets the categories
  of the per-dataset store.
- The `obs` index is the cell key, `<dataset_id>_<observation_joinid>`, which is
  unique across the corpus. The index of a per-dataset store is the author's
  cell name, which can occur in more than one dataset.
- `var` is the gene reference without `MEASURED_COLUMN`, because that column
  differs between datasets. `varm["measured"]` holds it instead: a boolean
  matrix of genes x members, with column `k` for member `k` of the record. It is
  the members' part of a presence matrix of genes and datasets.
- `uns` holds the keys that all members share with one value. The compaction
  record holds the `uns` of each member.

Compatibility
-------------
The members must compose as the SDK composes stores: one organism, one gene
reference, and one `var` apart from `MEASURED_COLUMN`; one `obs` column order;
the same `X` value type; all spatial or all non-spatial. `compact` refuses other
members. One further rule is the caller's: the members must share one access
grant, because the compacted store has one. Compact only public datasets for now.

Memory
------
`compact` holds the `obs` of all members in memory, and copies `X` in slices of
`X_SHARD_VALUES` values. A shard is written once, when it is full. Memory thus
follows the cell count of the group, not its stored values. The caller bounds
the group, as the compaction planner does.
"""

from __future__ import annotations

import os
from collections.abc import Callable, Iterator, Sequence
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Any, cast

import numpy as np
import pandas as pd
import zarr
from anndata.io import read_elem, write_elem

from . import as_store, assays, read_summary
from .store_format import (
    COMPACTION_KEY,
    COMPACTION_VERSION,
    DATASET_ID_COLUMN,
    MEASURED_COLUMN,
    SUMMARY_NAMESPACE,
    CompactedMember,
    CompactionRecord,
    obs_digest,
    validate_compaction,
)

# Rows per chunk and per shard of each `obs` array and of `X/indptr`. A read of one
# member reads the chunks that its rows overlap. The values are the layout that the
# design record gives for `obs`.
OBS_CHUNK_ROWS = 8192
OBS_SHARD_ROWS = 65536
# Values per chunk and per shard of `X/data` and `X/indices`. One shard of
# `X/indices` is 64 MiB before compression.
X_CHUNK_VALUES = 2**17
X_SHARD_VALUES = 2**23

JOINID_COLUMN = "observation_joinid"


class CompactionError(ValueError):
    """The members do not compose into one store."""


@dataclass
class _Member:
    """One per-dataset store, opened, with the parts that compaction reads first."""

    uri: str
    group: zarr.Group
    dataset_id: str
    dataset_version_id: str
    gene_reference_version: str
    obs: pd.DataFrame
    var: pd.DataFrame
    uns: dict[str, str]
    n_obs: int
    nnz: int


def _open(uri: str) -> _Member:
    """Open one per-dataset store and read its summary record, `obs`, `var`, and `uns`."""
    group = zarr.open_group(as_store(uri), mode="r")
    summary = read_summary(group)
    if summary is None:
        raise CompactionError(f"{uri}: no summary record; not a per-dataset store")
    x = group["X"]
    if not isinstance(x, zarr.Group) or x.attrs.get("encoding-type") != "csr_matrix":
        raise CompactionError(f"{uri}: X is not a CSR matrix")
    uns = cast("dict[str, Any]", read_elem(group["uns"]))
    obs = cast("pd.DataFrame", read_elem(group["obs"]))
    return _Member(
        uri=uri,
        group=group,
        dataset_id=summary["dataset_id"],
        dataset_version_id=summary["dataset_version_id"],
        gene_reference_version=summary["gene_reference_version"],
        obs=obs,
        var=cast("pd.DataFrame", read_elem(group["var"])),
        uns={str(k): str(v) for k, v in uns.items()},
        n_obs=len(obs),
        nnz=int(np.asarray(_x_part(group, "indptr")[-1:])[0]),
    )


def _check_members(members: Sequence[_Member]) -> None:
    """Refuse members that do not compose into one store. Name the conflict."""
    first = members[0]
    versions = [m.dataset_version_id for m in members]
    if len(set(versions)) != len(versions):
        raise CompactionError(f"a dataset version occurs twice: {sorted(versions)}")
    for field in ("dataset_id", "gene_reference_version"):
        values = {getattr(m, field) for m in members}
        if field == "dataset_id" and len(values) != len(members):
            raise CompactionError("a dataset_id occurs twice")
        if field == "gene_reference_version" and len(values) > 1:
            raise CompactionError(f"gene_reference_version differs: {sorted(values)}")
    organisms = {m.uns.get("organism_ontology_term_id") for m in members}
    if len(organisms) > 1:
        raise CompactionError(f"organism differs: {sorted(map(str, organisms))}")
    spatial = {_spatial(m.obs) for m in members}
    if len(spatial) > 1:
        raise CompactionError("spatial and non-spatial datasets in one group")
    columns = list(first.obs.columns)
    for m in members:
        if list(m.obs.columns) != columns:
            raise CompactionError(f"{m.uri}: obs columns differ from {first.uri}")
        if MEASURED_COLUMN not in m.var.columns:
            raise CompactionError(f"{m.uri}: var has no {MEASURED_COLUMN} column")
        if not _reference(m.var).equals(_reference(first.var)):
            raise CompactionError(f"{m.uri}: var differs from {first.uri}")
        if JOINID_COLUMN not in m.obs.columns or DATASET_ID_COLUMN not in m.obs.columns:
            raise CompactionError(f"{m.uri}: obs lacks {JOINID_COLUMN} or dataset_id")
        joinids = m.obs[JOINID_COLUMN]
        if bool(joinids.isna().any()) or not joinids.is_unique:
            raise CompactionError(f"{m.uri}: {JOINID_COLUMN} is null or not unique")
        found = {str(v) for v in m.obs[DATASET_ID_COLUMN].unique()}
        if found != {m.dataset_id}:
            raise CompactionError(f"{m.uri}: obs dataset_id holds {sorted(found)}")
    dtypes = {str(_x_part(m.group, "data").dtype) for m in members}
    if len(dtypes) > 1:
        raise CompactionError(f"X value types differ: {sorted(dtypes)}")


def _reference(var: pd.DataFrame) -> pd.DataFrame:
    """The part of a store's `var` that all stores of one organism share."""
    return var.drop(columns=MEASURED_COLUMN, errors="ignore")


def _spatial(obs: pd.DataFrame) -> bool:
    """Whether the dataset is spatial, by the rule of the catalog builder."""
    if "assay_ontology_term_id" not in obs.columns:
        return False
    terms = {str(t) for t in obs["assay_ontology_term_id"].unique()}
    return bool(assays.SPATIAL & terms)


def union_obs(frames: Sequence[pd.DataFrame]) -> pd.DataFrame:
    """The `obs` frames one after another, each categorical over the union of categories.

    A categorical column gets the sorted union of the members' categories, as
    strings. Its codes are remapped, so no string is made for each cell. Another
    column keeps its type. The index is the cell key,
    `<dataset_id>_<observation_joinid>`.
    """
    columns = list(frames[0].columns)
    data: dict[str, Any] = {}
    for name in columns:
        parts = [f[name] for f in frames]
        if all(isinstance(p.dtype, pd.CategoricalDtype) for p in parts):
            categories = pd.Index(
                sorted({str(c) for p in parts for c in p.cat.categories})
            )
            codes = [
                np.where(
                    p.cat.codes.to_numpy() < 0,
                    -1,
                    categories.get_indexer(p.cat.categories.astype(str))[
                        p.cat.codes.to_numpy()
                    ],
                )
                for p in parts
            ]
            data[name] = pd.Categorical.from_codes(
                np.concatenate(codes), categories=categories
            )
        else:
            data[name] = pd.concat(parts, ignore_index=True)
    dataset_ids = np.concatenate(
        [f[DATASET_ID_COLUMN].astype(str).to_numpy(dtype=object) for f in frames]
    )
    joinids = np.concatenate(
        [f[JOINID_COLUMN].astype(str).to_numpy(dtype=object) for f in frames]
    )
    index = pd.Index(dataset_ids + "_" + joinids)
    # `.array` keeps an extension type, such as the `str` type of pandas 3, and a
    # frame built from arrays does not align them on their index.
    return pd.DataFrame(
        {k: v.array if isinstance(v, pd.Series) else v for k, v in data.items()},
        index=index,
    )


class _ShardWriter:
    """Writes a 1-D array in order, one whole shard at a time.

    A shard written once, whole, needs no read of the shard before the write. The
    last shard is written by `close`, also once.
    """

    def __init__(self, array: zarr.Array[Any]) -> None:
        self.array = array
        self.shard = (array.shards or array.chunks)[0]
        self.parts: list[np.ndarray] = []
        self.buffered = 0
        self.position = 0

    def append(self, values: np.ndarray) -> None:
        self.parts.append(values)
        self.buffered += len(values)
        while self.buffered >= self.shard:
            self._write(self.shard)

    def close(self) -> None:
        if self.buffered:
            self._write(self.buffered)
        if self.position != self.array.shape[0]:
            raise CompactionError(
                f"wrote {self.position} values of {self.array.shape[0]}"
            )

    def _write(self, n: int) -> None:
        buffer = np.concatenate(self.parts) if len(self.parts) > 1 else self.parts[0]
        self.array[self.position : self.position + n] = buffer[:n]
        self.position += n
        rest = buffer[n:]
        self.parts = [rest] if len(rest) else []
        self.buffered = len(rest)


def _x_part(group: zarr.Group, name: str) -> zarr.Array[Any]:
    """One array of the CSR matrix `X` of a store: `data`, `indices`, or `indptr`."""
    x = group["X"]
    assert isinstance(x, zarr.Group), "X must be an encoded sparse matrix"
    array = x[name]
    assert isinstance(array, zarr.Array), f"X/{name} must be an array"
    return array


def _slices(n: int, step: int) -> Iterator[slice]:
    return (slice(i, min(i + step, n)) for i in range(0, n, step))


def _copy_x(members: Sequence[_Member], x: zarr.Group, n_var: int) -> np.ndarray:
    """Copy the `X` of each member into `x`, in member order. Return the new indptr."""
    n_obs = sum(m.n_obs for m in members)
    nnz = sum(m.nnz for m in members)
    first = _x_part(members[0].group, "data")
    first_indices = _x_part(members[0].group, "indices")
    x.attrs.update(
        {
            "encoding-type": "csr_matrix",
            "encoding-version": "0.1.0",
            "shape": [n_obs, n_var],
        }
    )
    writers = {
        name: _ShardWriter(
            x.create_array(
                name,
                shape=(nnz,),
                dtype=source.dtype,
                chunks=(X_CHUNK_VALUES,),
                shards=(X_SHARD_VALUES,),
            )
        )
        for name, source in (("data", first), ("indices", first_indices))
    }
    indptr = np.zeros(n_obs + 1, dtype=np.int64)
    row, offset = 0, 0
    for m in members:
        member_indptr = np.asarray(_x_part(m.group, "indptr")[:])
        n = len(member_indptr) - 1
        indptr[row + 1 : row + n + 1] = member_indptr[1:] + offset
        for name, writer in writers.items():
            array = _x_part(m.group, name)
            for part in _slices(m.nnz, X_SHARD_VALUES):
                writer.append(np.asarray(array[part]))
        row += n
        offset += m.nnz
    for writer in writers.values():
        writer.close()
    x.create_array(
        "indptr",
        shape=indptr.shape,
        dtype=indptr.dtype,
        chunks=(OBS_CHUNK_ROWS,),
        shards=(OBS_SHARD_ROWS,),
    )[:] = indptr
    return indptr


@contextmanager
def _threads(threads: int | None) -> Iterator[None]:
    """Bound zarr's codec thread pool, as `write` bounds dask's. See `write`."""
    if not threads:
        yield
        return
    with zarr.config.set({"threading.max_workers": threads}):
        yield


def compact(
    sources: Sequence[str | os.PathLike[str]],
    dest: str | os.PathLike[str],
    *,
    threads: int | None = None,
) -> CompactionRecord:
    """Write the per-dataset stores `sources` as one compacted store at `dest`.

    The members follow the order of `sources`. Return the compaction record, which
    is also the last write to the store. Raise `CompactionError` when the members
    do not compose; see the module docstring. `threads` bounds the threads that
    read the members and that encode chunks.
    """
    if not sources:
        raise CompactionError("a compacted store needs one or more members")
    uris = [str(s) for s in sources]
    with ThreadPoolExecutor(max_workers=threads or 8) as pool:
        members = list(pool.map(_open, uris))
    _check_members(members)
    # The assayed genes of each member, one column each. The check found the rest
    # of every `var` equal, so keep one copy: a group of many members otherwise
    # holds one copy of the reference for each.
    measured = np.column_stack(
        [m.var[MEASURED_COLUMN].to_numpy(dtype=bool) for m in members]
    )
    var = _reference(members[0].var)
    for m in members:
        m.var = var

    record: CompactionRecord = {
        "version": COMPACTION_VERSION,
        "gene_reference_version": members[0].gene_reference_version,
        "members": [],
    }
    row = 0
    for m in members:
        joinids = m.obs[JOINID_COLUMN].astype(str)
        entry: CompactedMember = {
            "dataset_id": m.dataset_id,
            "dataset_version_id": m.dataset_version_id,
            "row_start": row,
            "row_stop": row + len(m.obs),
            "nnz": m.nnz,
            "obs_digest": obs_digest(joinids),
            "uns": m.uns,
        }
        record["members"].append(entry)
        row += len(m.obs)
    validate_compaction(record)

    obs = union_obs([m.obs for m in members])
    shared = {
        k: v
        for k, v in members[0].uns.items()
        if all(m.uns.get(k) == v for m in members)
    }
    store = as_store(dest)
    with _threads(threads):
        root = zarr.open_group(store, mode="w", zarr_format=3)
        root.attrs.update({"encoding-type": "anndata", "encoding-version": "0.1.0"})
        obs_layout = {"chunks": (OBS_CHUNK_ROWS,), "shards": (OBS_SHARD_ROWS,)}
        write_elem(root, "obs", obs, dataset_kwargs=obs_layout)
        del obs
        for m in members:  # free the frames before the copy of X
            m.obs = m.obs.iloc[:0]
        write_elem(root, "var", var)
        write_elem(root, "uns", shared)
        # One chunk: the matrix is small, and a reader of one member reads it whole.
        write_elem(
            root,
            "varm",
            {MEASURED_COLUMN: measured},
            dataset_kwargs={"chunks": measured.shape},
        )
        for name in ("layers", "obsm", "obsp", "varp"):
            write_elem(root, name, {})
        _copy_x(members, root.create_group("X"), len(var))
        zarr.consolidate_metadata(store)
    # The record is the last write, as the summary record is in `write`. The root
    # group opens with its consolidated metadata, and the update keeps it.
    root = zarr.open_group(store, mode="r+")
    root.update_attributes({SUMMARY_NAMESPACE: {COMPACTION_KEY: dict(record)}})
    return record


def read_compaction(
    dest: str | os.PathLike[str] | zarr.Group,
) -> CompactionRecord | None:
    """The compaction record of a store, or None when the store has none.

    A record that is present but malformed raises `SummaryRecordError`.
    """
    root = (
        dest
        if isinstance(dest, zarr.Group)
        else zarr.open_group(as_store(dest), mode="r")
    )
    namespace = root.attrs.get(SUMMARY_NAMESPACE)
    if not isinstance(namespace, dict) or COMPACTION_KEY not in namespace:
        return None
    return validate_compaction(namespace[COMPACTION_KEY])


def verify(
    dest: str | os.PathLike[str],
    sources: Sequence[str | os.PathLike[str]],
    *,
    progress: Callable[[str], None] | None = None,
) -> list[str]:
    """Compare a compacted store with its per-dataset stores. Return the differences.

    For each member, the compacted rows must hold the per-dataset store's `obs`
    values and its `X`, value for value: the same `indptr` steps, `indices`, and
    `data`. A categorical value is compared as a string. The store's `var` must be
    the member's `var` without `MEASURED_COLUMN`, and the member's column of
    `varm["measured"]` must be that column. An empty list means that the stores
    agree. The check reads both copies whole, so it costs about as much as the
    compaction.
    """
    root = zarr.open_group(as_store(dest), mode="r")
    record = read_compaction(root)
    if record is None:
        return [f"{dest}: no compaction record"]
    if len(record["members"]) != len(sources):
        return [f"{len(record['members'])} members, {len(sources)} sources"]
    failures: list[str] = []
    obs = cast("pd.DataFrame", read_elem(root["obs"]))
    var = cast("pd.DataFrame", read_elem(root["var"]))
    varm = cast("dict[str, Any]", read_elem(root["varm"]))
    measured = np.asarray(varm.get(MEASURED_COLUMN, np.zeros((len(var), 0), bool)))
    if measured.shape != (len(var), len(sources)) or measured.dtype != bool:
        failures.append(
            f"varm {MEASURED_COLUMN} is {measured.dtype} {measured.shape}; "
            f"want bool {(len(var), len(sources))}"
        )
    indptr = np.asarray(_x_part(root, "indptr")[:])
    for k, (member, source) in enumerate(zip(record["members"], sources, strict=True)):
        if progress:
            progress(f"verify {member['dataset_version_id']}")
        name = member["dataset_version_id"]
        start, stop = member["row_start"], member["row_stop"]
        group = zarr.open_group(as_store(source), mode="r")
        summary = read_summary(group)
        if summary is None or summary["dataset_version_id"] != name:
            failures.append(f"{name}: {source} holds another dataset version")
            continue
        want_var = cast("pd.DataFrame", read_elem(group["var"]))
        if not _reference(want_var).equals(var):
            failures.append(f"{name}: var differs")
        assayed = want_var.get(MEASURED_COLUMN)
        if (
            assayed is None
            or k >= measured.shape[1]
            or not np.array_equal(assayed.to_numpy(dtype=bool), measured[:, k])
        ):
            failures.append(f"{name}: varm {MEASURED_COLUMN} differs")
        want = cast("pd.DataFrame", read_elem(group["obs"]))
        got = obs.iloc[start:stop]
        if len(want) != len(got) or list(want.columns) != list(got.columns):
            failures.append(f"{name}: obs shape or columns differ")
            continue
        for column in want.columns:
            a, b = want[column], got[column]
            if isinstance(a.dtype, pd.CategoricalDtype):
                a, b = a.astype(str), b.astype(str)
            if not np.array_equal(a.to_numpy(), b.to_numpy()):
                failures.append(f"{name}: obs {column} differs")
        want_indptr = np.asarray(_x_part(group, "indptr")[:])
        got_indptr = indptr[start : stop + 1]
        if not np.array_equal(np.diff(want_indptr), np.diff(got_indptr)):
            failures.append(f"{name}: X indptr differs")
            continue
        base = int(got_indptr[0])
        for part_name in ("data", "indices"):
            want_values = _x_part(group, part_name)
            got_values = _x_part(root, part_name)
            for part in _slices(member["nnz"], X_SHARD_VALUES):
                shifted = slice(base + part.start, base + part.stop)
                if not np.array_equal(want_values[part], got_values[shifted]):
                    failures.append(f"{name}: X {part_name} differs at {part.start}")
                    break
    return failures
