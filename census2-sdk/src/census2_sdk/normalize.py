"""The Census 1 normalized values, calculated from the counts on read.

Census 1 stores `X["normalized"]`. A Census 2 store holds only the counts, and
`get_anndata(..., layers=["normalized"])` calculates the values with the
Census 1 algorithm (`experiment_builder.py` of `cellxgene_census_builder`):

1. The nonzero counts as float32. Census 1 drops stored zeros first. A cell of a
   full-gene assay divides each count by the gene's `feature_length`.
2. Each value divided by its row sum, accumulated in float64.

Census 1 also rounds each value to 15 mantissa bits, so that its stored layer
compresses better. The SDK keeps the values in memory and does not round them, so
a value differs from Census 1 by less than 2**-15 relative. This bound holds while
each value is a normal float32. A value below 2**-126 needs a row sum more than
2**126 times the smallest value of its row.

The SDK reads whole rows, so it sums each row itself, in reference gene order.
Census 1 sums in the dataset's gene order. Integer counts give the same sum in
any order, but in a rare full-gene row the order changes the last float32 bit.
A read of part of a row needs a row-sum column in the store, summed in source gene
order. The writer does not write that column yet.
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


def normalized(
    counts: sp.csr_matrix,
    full_gene: npt.NDArray[np.bool_],
    feature_length: npt.NDArray[np.integer] | None,
) -> sp.csr_matrix:
    """The Census 1 normalized values of `counts`, unrounded, at its nonzero entries.

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
    # ponytail: reference order; a rare full-gene row differs in the last bit.
    row_sum = np.bincount(rows, weights=data, minlength=counts.shape[0])
    values = (data / row_sum[rows]).astype(np.float32)
    return sp.csr_matrix((values, counts.indices, counts.indptr), shape=counts.shape)
