"""The pinned gene references, and the projection of a dataset's genes onto them.

The reference is the gene file that the CELLxGENE dataset schema publishes for each
organism in `single-cell-curation`, pinned by repository tag. Its row order is the
gene position. Every store of one organism has the reference as its `var`, in this
order, so a reader concatenates stores by position (design 1.5).

`project` keeps the counts matrix: `raw.X` when the dataset has `raw`, else `X`. It
drops the features that are not genes of the dataset's organism, and moves each
remaining column to its reference position. A gene that the dataset did not measure
is a zero column. Thus the store does not record which genes the dataset measured.

The writer never fetches a reference. The caller stages the files into `refdata`,
and `load_reference` checks each file against its pinned SHA-256 digest.
`bulk-ingest-nf/tools/stage-refdata.py` stages them from the `url` of each entry.
"""

from __future__ import annotations

import gzip
import hashlib
import io
import os
import pathlib
from collections.abc import Mapping
from typing import Any, TypedDict, cast

import anndata
import dask.array as da
import numpy as np
import pandas as pd
import scipy.sparse as sp

# The pin: one single-cell-curation tag. Every organism's reference file and its
# annotation version come from this tag. A new tag is a reference migration
# (design 3.2), and it changes the store format.
CURATION_TAG = "v7.0.0"
_BASE = (
    "https://raw.githubusercontent.com/chanzuckerberg/single-cell-curation/"
    f"{CURATION_TAG}/cellxgene_schema_cli/cellxgene_schema/gencode_files"
)

# The columns of a reference file, which has no header row.
COLUMNS = [
    "feature_id",
    "feature_name",
    "feature_version",
    "feature_length",
    "feature_type",
]


class Reference(TypedDict):
    """One organism's reference file, as the configuration names it."""

    file: str  # resolved against `refdata`
    version: str  # recorded in the store as `uns["gene_reference_version"]`
    sha256: str
    url: str  # where `stage-refdata.py` fetches the file; the writer never does


# The default `gene_references` configuration, keyed by organism ontology term.
# gene_info.yml at the tag declares GENCODE v48 and M37, and the files have these
# digests.
REFERENCES: dict[str, Reference] = {
    "NCBITaxon:9606": {
        "file": "genes_homo_sapiens.csv.gz",
        "version": "GENCODE v48",
        "sha256": "2ef2138de00aa2b6def4f8779195e40d98bcb4d00556400b40b9cc8d51b4e047",
        "url": f"{_BASE}/genes_homo_sapiens.csv.gz",
    },
    "NCBITaxon:10090": {
        "file": "genes_mus_musculus.csv.gz",
        "version": "GENCODE M37",
        "sha256": "a6c827698cb3130df51e2394eae8dfc247c18a6c59a54d55da494891d166fef4",
        "url": f"{_BASE}/genes_mus_musculus.csv.gz",
    },
}


class GeneReferenceMismatch(ValueError):
    """The dataset holds a gene that is not in the pinned reference.

    The writer refuses the dataset. It does not map genes across reference versions.
    """


def load_reference(
    spec: Reference, refdata: str | os.PathLike[str] | None
) -> pd.DataFrame:
    """The reference as a `var`: indexed by `feature_id`, in reference order.

    Reads `spec["file"]` from `refdata` and checks its SHA-256 digest before
    anything decompresses it.
    """
    if refdata is None:
        raise ValueError(
            f"gene projection needs a refdata directory that holds {spec['file']}"
        )
    path = pathlib.Path(refdata) / spec["file"]
    if not path.exists():
        raise FileNotFoundError(f"{spec['file']} is not staged in refdata ({refdata})")
    data = path.read_bytes()
    digest = hashlib.sha256(data).hexdigest()
    if digest != spec["sha256"]:
        raise ValueError(
            f"{path} has digest {digest}; the configuration says {spec['sha256']}"
        )
    frame = pd.read_csv(io.BytesIO(gzip.decompress(data)), header=None, names=COLUMNS)
    if not frame["feature_id"].is_unique:
        raise ValueError(f"{path} has duplicate feature_ids")
    frame.index = pd.Index(frame["feature_id"].astype(str).to_numpy())
    return frame


def locate(
    var: pd.DataFrame, taxon: str, reference: pd.DataFrame, version: str
) -> tuple[np.ndarray, np.ndarray]:
    """The `var` rows to keep, and the reference position of each kept row.

    `var` is indexed by `feature_id`. When `var` has `feature_biotype` and
    `feature_reference`, the function drops a non-gene feature, such as an ERCC
    spike-in, and a gene of another organism. The Census 1 builder drops the same
    features (disposition row 12). A kept gene that is not in the reference raises
    `GeneReferenceMismatch`.
    """
    ids = pd.Index(var.index.astype(str))
    if not ids.is_unique:
        raise ValueError("var has duplicate feature_ids")

    kept = np.ones(len(var), dtype=bool)
    if "feature_biotype" in var.columns:
        kept &= (var["feature_biotype"].astype(str) == "gene").to_numpy()
    if "feature_reference" in var.columns:
        kept &= (var["feature_reference"].astype(str) == taxon).to_numpy()

    positions = reference.index.get_indexer(ids[kept])
    missing = positions < 0
    if bool(missing.any()):
        sample = list(ids[kept][missing][:5])
        raise GeneReferenceMismatch(
            f"{int(missing.sum())} of {int(kept.sum())} feature_ids are not in the "
            f"pinned reference {version} ({CURATION_TAG}). First missing: {sample}"
        )
    return kept, positions.astype(np.int64)


def _remap(
    block: Any, columns: np.ndarray, positions: np.ndarray, n_genes: int
) -> sp.csr_matrix:
    """One row block, with the kept columns moved to their reference positions."""
    csr = sp.csr_matrix(block)[:, columns]
    out = sp.csr_matrix(
        (csr.data, positions[csr.indices], csr.indptr), shape=(csr.shape[0], n_genes)
    )
    # The reference order is not the source order. sum_duplicates sorts the
    # indices of each row, and merges nothing because the feature_ids are unique.
    out.sum_duplicates()
    return out


def project_matrix(
    x: Any, kept: np.ndarray, positions: np.ndarray, n_genes: int
) -> da.Array | sp.csr_matrix:
    """`x` with its kept columns moved to reference positions. The result is CSR.

    A dask array stays lazy. Each row block is remapped when anndata writes the
    store, and it holds at most the stored values of its source block. A dense
    source becomes sparse, because most reference genes are zero in one dataset.
    """
    columns = np.flatnonzero(kept)
    if not isinstance(x, da.Array):
        return _remap(x, columns, positions, n_genes)
    if len(x.chunks[1]) != 1:
        raise ValueError("the matrix must be one block wide to remap its columns")
    # dask types `map_blocks` for a function of the blocks alone. This function
    # takes more arguments, so the call goes through Any, and a cast names the
    # result.
    map_blocks: Any = x.map_blocks
    return cast(
        "da.Array",
        map_blocks(
            _remap,
            columns,
            positions,
            n_genes,
            chunks=(x.chunks[0], (n_genes,)),
            dtype=x.dtype,
            meta=sp.csr_matrix((0, 0), dtype=x.dtype),
        ),
    )


def project(
    adata: anndata.AnnData,
    references: Mapping[str, Reference],
    refdata: str | os.PathLike[str] | None,
) -> anndata.AnnData:
    """Return a new AnnData whose `X` is the counts on the reference genes.

    `X` is `raw.X` when `raw` exists, else `X` (disposition row 8), and the result
    has no `raw`. `var` is the reference. `uns` gains `gene_reference_version`.
    `obs`, `uns`, `obsm`, and `obsp` pass through. `varm`, `varp`, and `layers` are
    indexed by the source genes, so they are dropped.
    """
    taxon = adata.uns.get("organism_ontology_term_id")
    if taxon is None:
        raise ValueError(
            "uns holds no organism_ontology_term_id, so the gene reference is unknown"
        )
    spec = references.get(str(taxon))
    if spec is None:
        raise ValueError(f"no gene reference is configured for {taxon}")
    reference = load_reference(spec, refdata)

    counts = adata.raw if adata.raw is not None else adata
    kept, positions = locate(counts.var, str(taxon), reference, spec["version"])
    x = project_matrix(counts.X, kept, positions, len(reference))

    uns = dict(adata.uns)
    uns["gene_reference_version"] = spec["version"]
    obsm: dict[str, Any] = {key: adata.obsm[key] for key in adata.obsm}
    obsp: dict[str, Any] = {key: adata.obsp[key] for key in adata.obsp}
    # `AnnData.obs` is typed `DataFrame | Dataset2D`, the second being what
    # `read_lazy` produces. `build` assembles obs with `read_elem`.
    return anndata.AnnData(
        X=x,
        obs=cast("pd.DataFrame", adata.obs),
        var=reference,
        uns=uns,
        obsm=obsm,
        obsp=obsp,
    )
