"""The cost of no proof engine: catalog candidates and column reads, two builders.

`run` tests each query of the earlier filter spike's corpus against a catalog, over
the SDK that this interpreter imports: this branch (`--api spike`) or the earlier
filter spike (`--api pr43`, a checkout of `spike/consolidation-filters`). The corpus
is the named queries of `bench/filters.py`, the bound queries of the harness corpus
for `homo_sapiens` (each SDK's own `harness.filter_language_builder.builder_of`), and
its probes. It reads only the catalog. For each query it records:

- `candidates`: the catalog rows that the SDK opens;
- `whole`: the candidates whose cells all match, so that the SDK reads no filter
  column of the store;
- `column_reads`: the filter columns that the SDK reads, summed over the candidates.
  This branch reads all of the query's columns from each candidate that is not
  whole. The earlier spike reads the columns of each candidate's residual filter.

`compare` joins the two records, and counts the queries where this branch opens
more stores or reads more columns.

    python cost.py run --api spike --catalog catalog-103.parquet \\
        --tools census2-sdk/tools --filters census2-sdk/tools/bench/filters.py \\
        --out cost-spike.json
    python cost.py compare cost-spike.json cost-pr43.json --out cost.json
"""

import argparse
import importlib.util
import json
import pathlib
import sys
from typing import Any


def spike_cost(census: Any, built: Any) -> dict[str, Any]:
    from census2_sdk.query import as_query

    q = as_query(census.homo_sapiens if built is None else built)
    rows = census.candidates(q)
    whole = int(rows["whole"].sum())
    return {
        "candidates": len(rows),
        "whole": whole,
        "columns": len(q.columns),
        "column_reads": (len(rows) - whole) * len(q.columns),
    }


def pr43_cost(census: Any, built: Any) -> dict[str, Any]:
    from census2_sdk.obs_filter import columns, is_const

    plan = census.homo_sapiens.plan(built)
    reads = [
        0 if is_const(c.residual, True) else len(columns(c.residual))
        for c in plan.candidates
    ]
    return {
        "candidates": len(plan.candidates),
        "whole": sum(1 for c in plan.candidates if is_const(c.residual, True)),
        "columns": len(plan.columns),
        "column_reads": sum(reads),
    }


def load_filters(path: pathlib.Path) -> Any:
    """`bench/filters.py` of this branch, by its path, beside another branch's tools."""
    spec = importlib.util.spec_from_file_location("spike_filters", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules["spike_filters"] = module
    spec.loader.exec_module(module)
    return module


def run(args: argparse.Namespace) -> None:
    sys.path.insert(0, str(args.tools))
    import warnings

    import census2_sdk
    from harness import filter_language_builder as flb
    from harness.corpus import expand, load_values

    census = census2_sdk.open_census([args.catalog])
    cost = spike_cost if args.api == "spike" else pr43_cost
    obs = census.homo_sapiens.obs
    out: dict[str, Any] = {}
    for name, query in load_filters(args.filters).QUERIES.items():
        out[f"filters:{name}"] = cost(census, query.build(obs))
    for bound in expand(load_values()):
        if bound.experiment != "homo_sapiens" or not bound.obs_filter:
            continue
        if bound.accessor in ("get_var", "get_presence_matrix"):
            continue
        built, _ = flb.builder_of(str(bound.obs_filter), obs)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                out[f"corpus:{bound.label}"] = cost(census, built)
        except ValueError as error:  # the earlier spike refuses an unknown label
            out[f"corpus:{bound.label}"] = {"refused": str(error)[:200]}
    for kind, build, _ in flb.PROBES:
        if build is not None:
            out[f"probe:{kind}"] = cost(census, build(obs))
    args.out.write_text(json.dumps(out, indent=1) + "\n")
    print(f"{len(out)} queries")


def compare(args: argparse.Namespace) -> None:
    ours = json.loads(args.spike.read_text())
    theirs = json.loads(args.pr43.read_text())
    rows = []
    for name, a in ours.items():
        b = theirs.get(name)
        if b is None or "refused" in a or "refused" in b:
            continue
        rows.append(
            {
                "query": name,
                "candidates": a["candidates"],
                "pr43_candidates": b["candidates"],
                "whole": a["whole"],
                "pr43_whole": b["whole"],
                "column_reads": a["column_reads"],
                "pr43_column_reads": b["column_reads"],
            }
        )
    more_stores = [r for r in rows if r["candidates"] > r["pr43_candidates"]]
    more_reads = [r for r in rows if r["column_reads"] > r["pr43_column_reads"]]
    fewer = [r for r in rows if r["candidates"] < r["pr43_candidates"]]
    summary = {
        "queries": len(rows),
        "more_stores": len(more_stores),
        "fewer_stores": len(fewer),
        "more_column_reads": len(more_reads),
        "extra_stores": sum(
            r["candidates"] - r["pr43_candidates"] for r in more_stores
        ),
        "extra_column_reads": sum(
            r["column_reads"] - r["pr43_column_reads"] for r in more_reads
        ),
    }
    for r in rows:
        mark = "*" if r in more_stores or r in more_reads else " "
        print(
            f"{mark} {r['query'][:60]:<60} {r['candidates']:>4}/{r['pr43_candidates']:<4} "
            f"whole {r['whole']:>3}/{r['pr43_whole']:<3} "
            f"reads {r['column_reads']:>4}/{r['pr43_column_reads']:<4}"
        )
    print(json.dumps(summary))
    if args.out is not None:
        args.out.write_text(
            json.dumps({"summary": summary, "queries": rows}, indent=1) + "\n"
        )


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = p.add_subparsers(dest="command", required=True)
    r = sub.add_parser("run")
    r.add_argument("--api", choices=["spike", "pr43"], required=True)
    r.add_argument("--catalog", required=True)
    r.add_argument("--tools", type=pathlib.Path, required=True, help="harness")
    r.add_argument(
        "--filters", type=pathlib.Path, required=True, help="bench/filters.py"
    )
    r.add_argument("--out", type=pathlib.Path, required=True)
    c = sub.add_parser("compare")
    c.add_argument("spike", type=pathlib.Path)
    c.add_argument("pr43", type=pathlib.Path)
    c.add_argument("--out", type=pathlib.Path, default=None)
    args = p.parse_args()
    if args.command == "run":
        run(args)
    else:
        compare(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
