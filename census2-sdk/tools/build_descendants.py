#!/usr/bin/env python3
"""Write the descendants table of the CxG ontologies: one Parquet row for each pair.

The SQL filter language selects the descendants of a term with a semi-join on this
table (`census2_sdk.obs_sql`, "Ontology closure"):

    cell_type_ontology_term_id IN (
        SELECT descendant FROM descendants WHERE ancestor = 'CL:0000084')

The columns are `ontology` (the prefix of the descendant), `ancestor`, and
`descendant`. Each term is its own descendant. The rows come from
`cellxgene_ontology_guide`, at the ontology versions of one CxG schema version; the
writer pins the same package and version for `tissue_general`. A release would
publish the table next to its catalog, so that a query and the stores see one version
of each ontology. The rows are sorted by `ancestor`, so that the Parquet statistics
of a row group exclude most row groups for one ancestor.

    python census2-sdk/tools/build_descendants.py descendants.parquet
"""

import argparse
import pathlib

import pyarrow as pa  # type: ignore[import-untyped]  # pyarrow has no types
import pyarrow.parquet as pq  # type: ignore[import-untyped]
from cellxgene_ontology_guide.ontology_parser import (  # type: ignore[import-untyped]
    OntologyParser,
)

# The ontologies of the obs columns with an ontology term, without NCBITaxon, whose
# closure has more than 30 million pairs and which no obs column holds per cell.
ONTOLOGIES = ("CL", "UBERON", "MONDO", "EFO", "HsapDv", "MmusDv", "PATO", "HANCESTRO")
SCHEMA_VERSION = "7.0.0"


def build(ontologies: tuple[str, ...], schema_version: str) -> pa.Table:
    parser = OntologyParser(f"v{schema_version}")
    ancestors: list[str] = []
    descendants: list[str] = []
    names: list[str] = []
    for name in ontologies:
        for term in parser.cxg_schema.ontology(name):
            for ancestor in parser.get_term_ancestors(term, include_self=True):
                ancestors.append(ancestor)
                descendants.append(term)
                names.append(name)
    table = pa.table(
        {"ontology": names, "ancestor": ancestors, "descendant": descendants}
    )
    return table.sort_by([("ancestor", "ascending"), ("descendant", "ascending")])


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("out", type=pathlib.Path)
    p.add_argument("--ontologies", default=",".join(ONTOLOGIES))
    p.add_argument("--schema-version", default=SCHEMA_VERSION)
    args = p.parse_args()
    table = build(tuple(args.ontologies.split(",")), args.schema_version)
    pq.write_table(table, args.out, compression="zstd", row_group_size=65_536)
    size = args.out.stat().st_size
    print(f"{table.num_rows} pairs, {size / 2**20:.1f} MiB -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
