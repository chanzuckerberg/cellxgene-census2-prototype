"""The pinned gene reference, and the projection of a store's `var` onto it.

The reference is the gene file that the CELLxGENE dataset schema publishes for each
organism in `single-cell-curation`, pinned by repository tag. Its row order is the gene
position, so stores of one organism concatenate by position (design 1.5).

The SDK projects `var` because the Iteration 1 stores carry the source `var`. When the
writer projects `var` at ingestion, `project` reduces to a check that the store's
`feature_id` order equals the reference.
"""

import gzip
import hashlib
import io
import os
import pathlib
import tempfile
import urllib.request
from dataclasses import dataclass
from typing import TypedDict

import numpy as np
import pandas as pd

# The pin: one single-cell-curation tag. Every organism's reference file and its
# annotation version come from this tag. A new tag is a reference migration
# (design 3.2).
CURATION_TAG = "v7.0.0"
_BASE = (
    "https://raw.githubusercontent.com/chanzuckerberg/single-cell-curation/"
    f"{CURATION_TAG}/cellxgene_schema_cli/cellxgene_schema/gencode_files"
)
DEFAULT_REFDATA = pathlib.Path.home() / ".cache" / "census2-sdk" / CURATION_TAG

COLUMNS = [
    "feature_id",
    "feature_name",
    "feature_version",
    "feature_length",
    "feature_type",
]


class Organism(TypedDict):
    taxon: str
    version: str
    n_genes: int
    sha256: str


# Recorded at pinning time: gene_info.yml at the tag declares GENCODE v48 and M37,
# and the files have these digests and row counts.
ORGANISMS: dict[str, Organism] = {
    "homo_sapiens": {
        "taxon": "NCBITaxon:9606",
        "version": "GENCODE v48",
        "n_genes": 78894,
        "sha256": "2ef2138de00aa2b6def4f8779195e40d98bcb4d00556400b40b9cc8d51b4e047",
    },
    "mus_musculus": {
        "taxon": "NCBITaxon:10090",
        "version": "GENCODE M37",
        "n_genes": 78317,
        "sha256": "a6c827698cb3130df51e2394eae8dfc247c18a6c59a54d55da494891d166fef4",
    },
}
TAXON_TO_ORGANISM = {spec["taxon"]: name for name, spec in ORGANISMS.items()}


class ReferencePinMismatch(Exception):
    """The store's genes are not from the pinned reference.

    The SDK refuses the store. It does not remap genes during a query (design 4.3).
    """


@dataclass(frozen=True)
class Projection:
    """A store's `var` mapped onto the pinned reference."""

    positions: np.ndarray  # reference position for each kept var row, in var order
    kept: np.ndarray  # boolean over the input var rows
    gene_reference_version: str


def load_reference(
    organism: str, refdata: str | os.PathLike[str] | None = None
) -> pd.DataFrame:
    """The pinned reference table. Row order is the gene position.

    `refdata` is a directory that holds `genes_<organism>.csv.gz`, staged by the
    caller and read as given, so a test can supply a small reference. With `None`,
    the SDK downloads the file once into the user cache and checks its SHA-256 digest
    and row count against the pin on every load.
    """
    spec = ORGANISMS[organism]
    if refdata is not None:
        return _read(pathlib.Path(refdata) / f"genes_{organism}.csv.gz")
    path = DEFAULT_REFDATA / f"genes_{organism}.csv.gz"
    if not path.exists():
        _download(f"{_BASE}/genes_{organism}.csv.gz", path, spec)
    _verify(path, spec)
    frame = _read(path)
    if len(frame) != spec["n_genes"]:
        path.unlink()
        raise ValueError(
            f"{path.name} has {len(frame)} rows; the pin says {spec['n_genes']}"
        )
    return frame


def _download(url: str, path: pathlib.Path, spec: Organism) -> None:
    """Fetch into a unique temporary file, check the digest, then replace atomically."""
    DEFAULT_REFDATA.mkdir(parents=True, exist_ok=True)
    descriptor, name = tempfile.mkstemp(dir=DEFAULT_REFDATA, suffix=".part")
    os.close(descriptor)
    partial = pathlib.Path(name)
    try:
        urllib.request.urlretrieve(url, partial)
        _verify(partial, spec)
        os.replace(partial, path)
    finally:
        partial.unlink(missing_ok=True)


def _verify(path: pathlib.Path, spec: Organism) -> None:
    """Check the file's SHA-256 against the pin before anything decompresses it."""
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != spec["sha256"]:
        path.unlink()
        raise ValueError(
            f"{path.name} has digest {digest}; the pin {CURATION_TAG} says {spec['sha256']}"
        )


def _read(path: pathlib.Path) -> pd.DataFrame:
    frame = pd.read_csv(
        io.BytesIO(gzip.decompress(path.read_bytes())), header=None, names=COLUMNS
    )
    if not frame["feature_id"].is_unique:
        raise ValueError(f"{path} has duplicate feature_ids")
    return frame


def reference_var(reference: pd.DataFrame) -> pd.DataFrame:
    """The reference as the `var` of a query result: indexed by `feature_id`."""
    var = reference.copy()
    var.index = pd.Index(var["feature_id"].astype(str).to_numpy())
    return var


def project(var: pd.DataFrame, organism: str, reference: pd.DataFrame) -> Projection:
    """Map a store's `var` onto `reference`.

    `var` is indexed by `feature_id`. `feature_biotype` and `feature_reference` are
    applied when present: a non-gene feature or a feature of another organism is
    dropped, as the Census 1 builder drops them (disposition row 12). A gene id that
    is absent from the reference raises `ReferencePinMismatch`.
    """
    spec = ORGANISMS[organism]
    ids = pd.Index(var.index.astype(str))
    if not ids.is_unique:
        raise ValueError("the store's var has duplicate feature_ids")

    kept = np.ones(len(var), dtype=bool)
    if "feature_biotype" in var.columns:
        kept &= (var["feature_biotype"].astype(str) == "gene").to_numpy()
    if "feature_reference" in var.columns:
        kept &= (var["feature_reference"].astype(str) == spec["taxon"]).to_numpy()

    position = pd.Series(
        np.arange(len(reference)), index=reference["feature_id"].astype(str)
    )
    looked_up = position.reindex(ids[kept])
    missing = looked_up.isna()
    if bool(missing.any()):
        sample = list(looked_up.index[missing][:5])
        raise ReferencePinMismatch(
            f"{int(missing.sum())} of {int(kept.sum())} feature_ids are not in the pinned "
            f"reference {spec['version']} ({CURATION_TAG}). First missing: {sample}"
        )
    return Projection(
        positions=looked_up.to_numpy(dtype=np.int64),
        kept=kept,
        gene_reference_version=spec["version"],
    )
