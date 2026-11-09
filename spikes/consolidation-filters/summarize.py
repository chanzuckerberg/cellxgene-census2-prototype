"""Summarize the measures of the filter builder spike as markdown tables.

The inputs are the JSON-lines files under `census2-sdk/tools/results/consolidation_filters/`:
`runs-ec2.jsonl` (`bench.matrix` of this branch: `main`, the mapping path, and the builder over
per-dataset, consolidated, and the DuckDB spike's stores), `duckdb-ec2.jsonl` (`bench.matrix` of
the DuckDB spike, via `sql`), and `bind-ec2.jsonl` (`bench.bind`). Each figure is the median of
the repeats.

    python summarize.py census2-sdk/tools/results/consolidation_filters
"""

import argparse
import json
import pathlib
import statistics
from collections import defaultdict
from typing import Any

# The DuckDB spike's shapes, by the name of the query here.
DUCKDB_SHAPES = {
    "single": "one_dataset",
    "cell_type": "neutrophil",
    "corpus_obs": "all_obs",
    "obs:wmg": "wmg",
    "count_donor:neutrophil": "count_b_donor",
    "count_donor:all_obs": "count_all_donor",
    "count:wmg:tissue_general_ontology_term_id+cell_type_ontology_term_id": "wmg_count",
}


def load(path: pathlib.Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def duckdb_name(shape: str) -> str:
    if shape in DUCKDB_SHAPES:
        return DUCKDB_SHAPES[shape]
    return shape.split(":", 1)[1] if ":" in shape else shape


def medians(rows: list[dict[str, Any]], key: Any) -> dict[Any, dict[str, float]]:
    groups: dict[Any, list[dict[str, Any]]] = defaultdict(list)
    for r in rows:
        if r.get("status") == "ok":
            groups[key(r)].append(r)
    out = {}
    for k, rs in groups.items():
        out[k] = {
            f: statistics.median(float(r[f]) for r in rs)
            for f in (
                "query_seconds",
                "cpu_seconds",
                "requests",
                "mb_payload",
                "peak_rss_mb",
            )
            if all(f in r for r in rs)
        }
        out[k]["cells"] = rs[0].get("cells", 0)
        out[k]["stores_opened"] = rs[0].get("stores_opened", 0)
        out[k]["n"] = len(rs)
    return out


def fmt(value: float | None, digits: int = 2) -> str:
    if value is None:
        return "–"
    if digits == 0 or value >= 100:
        return f"{value:,.0f}"
    return f"{value:.{digits}f}"


def runs_table(
    runs: list[dict[str, Any]], duck: list[dict[str, Any]], workers: int
) -> str:
    mine = medians(
        [r for r in runs if r.get("workers") == workers],
        lambda r: (r["config"], r["query"]),
    )
    theirs = medians(
        [r for r in duck if r.get("workers") == workers and r.get("via") == "sql"],
        lambda r: ("duckdb", duckdb_name(str(r["shape"]))),
    )
    allm = {**mine, **theirs}
    configs = [
        "main",
        "main_lean",
        "mapping_pd",
        "builder_pd",
        "builder_cons",
        "builder_ice",
        "duckdb",
    ]
    queries = list(dict.fromkeys(q for (_, q) in sorted(mine, key=lambda k: k[1])))
    lines = [
        f"Seconds at {workers} workers (median), then requests, CPU s, peak RSS GB:",
        "",
        "| query | cells | " + " | ".join(configs) + " |",
        "|---|---|" + "---|" * len(configs),
    ]
    for q in queries:
        cells = next((allm[(c, q)]["cells"] for c in configs if (c, q) in allm), 0)
        cols = [fmt(allm.get((c, q), {}).get("query_seconds")) for c in configs]
        lines.append(f"| {q} | {cells:,} | " + " | ".join(cols) + " |")
    lines += [
        "",
        "| query | " + " | ".join(configs) + " |",
        "|---|" + "---|" * len(configs),
    ]
    for q in queries:
        cols = []
        for c in configs:
            m = allm.get((c, q))
            if m is None:
                cols.append("–")
                continue
            cols.append(
                f"{fmt(m.get('requests'), 0)}, {fmt(m.get('cpu_seconds'), 1)}, "
                f"{fmt(m.get('peak_rss_mb', 0) / 1024, 1)}"
            )
        lines.append(f"| {q} | " + " | ".join(cols) + " |")
    return "\n".join(lines)


def per_store(stores: dict[str, list[dict[str, Any]]], f: Any) -> list[float]:
    """For each store, the median of `f` over its repeats."""
    return [statistics.median(f(r) for r in v) for v in stores.values()]


def bind_table(rows: list[dict[str, Any]]) -> str:
    by: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for r in rows:
        by[(r["layout"], r["query"])].append(r)
    lines = [
        "Per store, median over stores of the median over repeats, in ms; sums over stores in s:",
        "",
        "| layout | query | stores | columns | main read_elem(obs) ms | bind ms (read + eval) | eval ms | main CPU ms | bind CPU ms | main req | bind req | sum main s | sum bind s |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for (layout, query), rs in sorted(by.items()):
        stores: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for r in rs:
            stores[r["store"]].append(r)

        main = per_store(stores, lambda r: r["main"]["ms"])
        bind = per_store(stores, lambda r: r["bind_ms"])
        ev = per_store(stores, lambda r: r["eval"]["ms"])
        main_cpu = per_store(stores, lambda r: r["main"]["cpu_ms"])
        bind_cpu = per_store(stores, lambda r: r["bind_cpu_ms"])
        main_req = per_store(stores, lambda r: r["main"]["requests"])
        bind_req = per_store(stores, lambda r: r["read"]["requests"])
        lines.append(
            f"| {layout} | {query} | {len(stores)} | {len(rs[0]['columns'])} | "
            f"{fmt(statistics.median(main), 1)} | {fmt(statistics.median(bind), 1)} | "
            f"{fmt(statistics.median(ev), 2)} | {fmt(statistics.median(main_cpu), 1)} | "
            f"{fmt(statistics.median(bind_cpu), 1)} | {fmt(statistics.median(main_req), 0)} | "
            f"{fmt(statistics.median(bind_req), 0)} | {fmt(sum(main) / 1e3, 2)} | "
            f"{fmt(sum(bind) / 1e3, 2)} |"
        )
    return "\n".join(lines)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("results", type=pathlib.Path)
    args = p.parse_args()
    runs = load(args.results / "runs-ec2.jsonl")
    duck = load(args.results / "duckdb-ec2.jsonl")
    bind = load(args.results / "bind-ec2.jsonl")
    for workers in (8, 1, 32):
        print(runs_table(runs, duck, workers))
        print()
    print(bind_table(bind))
    laptop = load(args.results / "runs-laptop.jsonl")
    if laptop:
        print()
        print(runs_table(laptop, [], 8))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
