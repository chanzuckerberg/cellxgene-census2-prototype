"""The Census 1 normalized values, calculated from the counts on read (tracker 4.13).

Census 1 stores `X["normalized"]`. A Census 2 store holds only the counts, and
`get_anndata(..., X_layers=["normalized"])` calculates the same values with the
Census 1 algorithm (`experiment_builder.py` of `cellxgene_census_builder`):

1. The nonzero counts as float32. Census 1 drops stored zeros first. A cell of a
   full-gene assay divides each count by the gene's `feature_length`.
2. Each value divided by its row sum, accumulated in float64.
3. Each value rounded to 15 mantissa bits, half to even.
4. A value that rounds to zero replaced by the smallest float32 subnormal, so the
   layer has the sparsity of the nonzero counts.

The SDK reads whole rows, so it sums each row itself, in reference gene order.
Census 1 sums in the dataset's gene order. Integer counts give the same sum in
any order, but in a rare full-gene row the order changes the rounded value. The
row-sum column that the writer adds in source order removes that difference, and
a read of part of a row needs it (tracker 4.13).
"""

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp

# Census 1 `FULL_GENE_ASSAY` (`globals.py`): the assays with a gene-length correction.
FULL_GENE_ASSAYS = frozenset(
    {
        "EFO:0003755",  # FL-cDNA
        "EFO:0008747",  # FRISCR
        "EFO:0008763",  # Hi-SCL
        "EFO:0008797",  # MATQ-seq
        "EFO:0008877",  # Quartz-seq
        "EFO:0008930",  # Smart-seq
        "EFO:0008931",  # Smart-seq2
        "EFO:0008956",  # SUPeR-seq
        "EFO:0009999",  # G&T-Seq
        "EFO:0010004",  # SCRB-seq
        "EFO:0010006",  # scM&T-seq
        "EFO:0010022",  # Smart-3Seq
        "EFO:0010058",  # Fluidigm C1-based SMARTer library preparation
        "EFO:0010184",  # Smart-like
        "EFO:0022396",  # TruSeq
        "EFO:0022488",  # Smart-seq3
        "EFO:0022839",  # STORM-seq
        "EFO:0030031",  # SCOPE-chip
        "EFO:0030060",  # mCT-seq
        "EFO:0030061",  # mcSCRB-seq
        "EFO:0700016",  # Smart-seq v4
    }
)

KEEP_BITS = 15
SIGMA = np.finfo(np.float32).smallest_subnormal


def normalized(
    counts: sp.csr_matrix,
    full_gene: npt.NDArray[np.bool_],
    feature_length: npt.NDArray[np.integer] | None,
) -> sp.csr_matrix:
    """The Census 1 normalized values of `counts`, at its nonzero entries.

    `full_gene` marks the rows of a full-gene assay. `feature_length` is the
    length of each column's gene, and it is needed only when a row is marked.
    """
    counts = counts.copy()
    counts.eliminate_zeros()  # as Census 1 does before it normalizes
    rows = np.repeat(np.arange(counts.shape[0]), np.diff(counts.indptr))
    data = counts.data.astype(np.float32)
    if full_gene.any():
        assert feature_length is not None, "a full-gene assay needs feature_length"
        corrected = full_gene[rows]
        data = np.where(corrected, data / feature_length[counts.indices], data)
        data = data.astype(np.float32)
    # ponytail: reference order; rare full-gene rows differ (module docstring).
    row_sum = np.bincount(rows, weights=data, minlength=counts.shape[0])
    values = (data / row_sum[rows]).astype(np.float32)
    _round_half_to_even(values, KEEP_BITS)
    values[values == 0] = SIGMA
    return sp.csr_matrix((values, counts.indices, counts.indptr), shape=counts.shape)


def _round_half_to_even(a: npt.NDArray[np.float32], keepbits: int) -> None:
    """Round float32 values in place to `keepbits` mantissa bits, half to even.

    The bit operations of Census 1 `_roundHalfToEven`, on the int32 view.
    """
    maskbits = 23 - keepbits
    b = a.view(np.int32)
    b += ((b >> maskbits) & 1) + np.int32((1 << (maskbits - 1)) - 1)
    b &= np.int32(-(1 << maskbits))
