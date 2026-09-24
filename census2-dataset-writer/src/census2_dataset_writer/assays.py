"""Assay ontology terms that Census 2 admits.

The default is the CELLxGENE Census 1 list. It merges `RNA_SEQ` and
`ALLOWED_SPATIAL_ASSAYS` from the `chanzuckerberg/cellxgene-census` repository,
at `tools/cellxgene_census_builder/src/cellxgene_census_builder/build_soma/globals.py`.
Census 1 keeps the two apart because it writes them to separate experiments.
Census 2 admits a dataset from either, so this module merges them.

The list is a starting point, not a decision. It is what Census 1 covers today,
which keeps the two corpora comparable while Census 2 is a prototype. Each term
carries its label, so a reader does not need a second lookup.
"""

from __future__ import annotations

# RNA sequencing assays. Census 1 `RNA_SEQ`.
RNA_SEQ: frozenset[str] = frozenset(
    {
        "EFO:0003755",  # FL-cDNA
        "EFO:0008640",  # 3'T-fill
        "EFO:0008641",  # 3’-end-seq
        "EFO:0008643",  # 3′-Seq
        "EFO:0008661",  # Bru-Seq
        "EFO:0008669",  # CAGEscan
        "EFO:0008673",  # CapSeq
        "EFO:0008675",  # CaptureSeq
        "EFO:0008679",  # CEL-seq
        "EFO:0008694",  # ClickSeq
        "EFO:0008697",  # cP-RNA-Seq
        "EFO:0008703",  # CytoSeq
        "EFO:0008708",  # DeepCAGE
        "EFO:0008710",  # Digital RNA
        "EFO:0008718",  # DP-Seq
        "EFO:0008720",  # DroNc-seq
        "EFO:0008722",  # Drop-seq
        "EFO:0008735",  # FACS-seq
        "EFO:0008747",  # FRISCR
        "EFO:0008748",  # FRT-Seq
        "EFO:0008752",  # GMUCT 1.0
        "EFO:0008753",  # GMUCT 2.0
        "EFO:0008756",  # GRO-CAP
        "EFO:0008763",  # Hi-SCL
        "EFO:0008780",  # inDrop
        "EFO:0008796",  # MARS-seq
        "EFO:0008797",  # MATQ-seq
        "EFO:0008824",  # NanoCAGE
        "EFO:0008825",  # Nanogrid RNA-Seq
        "EFO:0008826",  # NET-Seq
        "EFO:0008850",  # PAS-Seq
        "EFO:0008859",  # PEAT
        "EFO:0008863",  # PLATE-Seq
        "EFO:0008868",  # PRO-cap
        "EFO:0008869",  # PRO-seq
        "EFO:0008877",  # Quartz-seq
        "EFO:0008896",  # RNA-Seq
        "EFO:0008897",  # RNAtag-Seq
        "EFO:0008898",  # RNET-seq
        "EFO:0008903",  # SC3-seq
        "EFO:0008919",  # Seq-Well
        "EFO:0008929",  # SMA
        "EFO:0008930",  # Smart-seq
        "EFO:0008931",  # Smart-seq2
        "EFO:0008937",  # snDrop-seq
        "EFO:0008941",  # sNuc-Seq
        "EFO:0008945",  # SPET-seq
        "EFO:0008953",  # STRT-seq
        "EFO:0008954",  # STRT-seq-2i
        "EFO:0008956",  # SUPeR-seq
        "EFO:0008962",  # TARDIS
        "EFO:0008966",  # TCR Chain Paring
        "EFO:0008967",  # TCR-LA-MC PCR
        "EFO:0008972",  # TL-seq
        "EFO:0008974",  # Tomo-Seq
        "EFO:0008975",  # TRAP-Seq
        "EFO:0008978",  # TSS Sequencing
        "EFO:0008980",  # UMI Method
        "EFO:0009309",  # Div-Seq
        "EFO:0009899",  # 10x 3' v2
        "EFO:0009900",  # 10x 5' v2
        "EFO:0009901",  # 10x 3' v1
        "EFO:0009919",  # SPLiT-seq
        "EFO:0009922",  # 10x 3' v3
        "EFO:0009991",  # Nuc-Seq
        "EFO:0009999",  # G&T-Seq
        "EFO:0010003",  # RASL-seq
        "EFO:0010004",  # SCRB-seq
        "EFO:0010005",  # DR-Seq
        "EFO:0010006",  # scM&T-seq
        "EFO:0010007",  # scTrio-seq
        "EFO:0010010",  # CEL-seq2
        "EFO:0010022",  # Smart-3Seq
        "EFO:0010034",  # Cappable-Seq
        "EFO:0010041",  # Nascent-Seq
        "EFO:0010058",  # Fluidigm C1-based SMARTer library preparation
        "EFO:0010184",  # Smart-like
        "EFO:0010550",  # sci-RNA-seq
        "EFO:0011025",  # 10x 5' v1
        "EFO:0022396",  # TruSeq
        "EFO:0022488",  # Smart-seq3
        "EFO:0022490",  # ScaleBio single cell RNA sequencing
        "EFO:0022600",  # Parse Evercode Whole Transcriptome v1
        "EFO:0022601",  # Parse Evercode Whole Transcriptome v2
        "EFO:0022602",  # Parse Evercode Whole Transcriptome v3
        "EFO:0022604",  # 10x 3' v4
        "EFO:0022605",  # 10x 5' v3
        "EFO:0022606",  # 10x gene expression flex
        "EFO:0022839",  # STORM-seq
        "EFO:0022845",  # modified STRT-seq
        "EFO:0022846",  # 5' STRT-seq
        "EFO:0022962",  # SHARE-seq
        "EFO:0030001",  # sci-CAR
        "EFO:0030002",  # microwell-seq
        "EFO:0030003",  # 10x 3' transcription profiling
        "EFO:0030004",  # 10x 5' transcription profiling
        "EFO:0030019",  # Seq-Well S3
        "EFO:0030021",  # Nx1-seq
        "EFO:0030026",  # sci-Plex
        "EFO:0030028",  # sci-RNA-seq3
        "EFO:0030030",  # Quant-seq
        "EFO:0030031",  # SCOPE-chip
        "EFO:0030059",  # 10x multiome
        "EFO:0030060",  # mCT-seq
        "EFO:0030061",  # mcSCRB-seq
        "EFO:0030074",  # SORT-seq
        "EFO:0700003",  # BD Rhapsody Whole Transcriptome Analysis
        "EFO:0700004",  # BD Rhapsody Targeted mRNA
        "EFO:0700010",  # TruDrop
        "EFO:0700011",  # GEXSCOPE technology
        "EFO:0700016",  # Smart-seq v4
        "EFO:0900000",  # particle-templated instant partition sequencing
        "EFO:0900001",  # Asteria scRNA-seq kit
        "EFO:0900002",  # HIVE CLX Single-Cell RNAseq Solution
    }
)

# Spatial assays. Census 1 `ALLOWED_SPATIAL_ASSAYS`.
SPATIAL: frozenset[str] = frozenset(
    {
        "EFO:0022857",  # Visium Spatial Gene Expression V1
        "EFO:0030062",  # Slide-seqV2
    }
)

# The default `assay_allowlist`. A dataset is admitted when it holds at least
# one cell on a term in this set. See `census2_dataset_writer.admit`.
DEFAULT_ALLOWLIST: frozenset[str] = RNA_SEQ | SPATIAL
