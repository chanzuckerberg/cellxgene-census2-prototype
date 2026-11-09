"""Summarize the measures of the filter builder spike as the markdown tables of the notes.

The inputs are the files under `census2-sdk/tools/results/consolidation_filters/`:

- `runs-ec2.jsonl`: `bench.matrix` of this branch: `main`, the mapping path, and the
  builder over per-dataset stores, consolidated stores, and the DuckDB spike's stores.
- `runs-ec2-obs.jsonl`: a later batch of the builder's `get_obs` queries (`wmg`,
  `all_obs`), after a change that reads the cell key once. Its records replace those of
  `runs-ec2.jsonl` for these queries.
- `duckdb-ec2.jsonl`: `bench.matrix` of the DuckDB spike, via `sql`.
- `bind-ec2.jsonl`: `bench.bind`.
- `runs-laptop.jsonl`: a few queries from a laptop.
- `queries-duckdb.json`: the names of the DuckDB spike's filters.

Each figure is the median of the repeats.

    python summarize.py census2-sdk/tools/results/consolidation_filters
"""

import argparse
import json
import pathlib
import statistics
from collections import defaultdict
from typing import Any

# The DuckDB spike's records name a shape and hold the filter, not its name. These
# shapes have one query each; a `query`, `obs`, or `count` record takes the name whose
# SQL it holds, from `queries-duckdb.json`.
DUCKDB_SHAPES = {
    "single": "one_dataset",
    "cell_type": "neutrophil",
    "corpus_obs": "all_obs",
}
# The queries that the benchmark counts (`count_by` in `bench.filters`).
COUNTED = {
    "count_b_donor",
    "count_all_donor",
    "wmg_count",
    "compound",
    "has_any",
    "between",
    "not_in",
    "prefix",
    "substring",
    "regex",
}
SQL: dict[str, list[str]] = {}
REPLACED = {"wmg", "all_obs"}

SELECTIONS = [
    ("one_dataset", "(a) one dataset"),
    ("single_member", "one dataset, a member"),
    ("single_alone", "one dataset, alone"),
    ("neutrophil", "(b) neutrophils"),
    ("lung_macrophage", "lung macrophages (mapping form)"),
    ("lung_macrophage_primary", "(c) lung macrophages, primary"),
    ("fp_none", "(d1) false positive, no cell"),
    ("fp_some", "(d2) false positive, 2 of 4"),
    ("false_positive", "false positive of the consolidation spike"),
    ("tissue_general", "(f) `tissue_general == 'pancreas'`"),
    ("numeric", "(b) & `nnz >= 1000`"),
    ("membership", "COVID-19 neutrophils (`has_any`)"),
    ("or_columns", "`\\|` across columns"),
    ("has_all", "`has_all`, the joined value"),
    ("numeric_range", "(b) & `nnz >= 2000` & `raw_sum < 20000`"),
    ("isna", "`raw_variance_nnz.isna()`"),
]
OBS_ONLY = [
    ("count_b", "count (b)"),
    ("count_b_donor", "count (b) by donor"),
    ("count_all_donor", "count all by donor"),
    ("wmg_count", "(g) WMG, count by tissue and cell type"),
    ("compound", "compound `\\|` and `~`, by donor"),
    ("not_in", "`NOT IN`, by donor"),
    ("has_any", "`has_any` of two terms, by donor"),
    ("between", "`between`, by donor"),
    ("prefix", "`.str.startswith`, by donor"),
    ("substring", "`.str.contains`, by donor"),
    ("regex", "`.str.match`, by donor"),
    ("wmg", "WMG `get_obs` (all columns)"),
    ("all_obs", "(e) `get_obs` of all cells"),
]


def load(path: pathlib.Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def duckdb_name(record: dict[str, Any]) -> str:
    """The name of the query of a record of the DuckDB spike's benchmark."""
    shape = str(record["shape"])
    if shape in DUCKDB_SHAPES:
        return DUCKDB_SHAPES[shape]
    names = SQL.get(str(record.get("query")), [])
    if shape == "obs":
        names = [n for n in names if n == "wmg"]
    elif shape == "count":
        names = [n for n in names if n in COUNTED]
    else:
        names = [n for n in names if n not in COUNTED and n != "wmg"]
    return names[0] if names else f"{shape}?"


Medians = dict[tuple[str, str, int], dict[str, float]]


def medians(rows: list[dict[str, Any]], key: Any) -> Medians:
    groups: dict[Any, list[dict[str, Any]]] = defaultdict(list)
    for r in rows:
        if r.get("status") == "ok":
            groups[key(r)].append(r)
    out = {}
    for k, rs in groups.items():
        fields = (
            "query_seconds",
            "cpu_seconds",
            "requests",
            "mb_payload",
            "peak_rss_mb",
        )
        out[k] = {
            f: statistics.median(float(r[f]) for r in rs)
            for f in fields
            if all(f in r for r in rs)
        }
        out[k]["cells"] = rs[0].get("cells", 0)
        out[k]["stores_opened"] = rs[0].get("stores_opened", 0)
    return out


def fmt(value: float | None, digits: int = 2) -> str:
    if value is None:
        return "–"
    if digits == 0 or value >= 100:
        return f"{value:,.0f}"
    return f"{value:.{digits}f}"


def all_medians(results: pathlib.Path, runs_name: str = "runs-ec2.jsonl") -> Medians:
    runs = load(results / runs_name)
    later = load(results / "runs-ec2-obs.jsonl")
    if later:
        runs = [r for r in runs if r.get("query") not in REPLACED] + later
    mine = medians(runs, lambda r: (r["config"], r["query"], r["workers"]))
    duck = [r for r in load(results / "duckdb-ec2.jsonl") if r.get("via") == "sql"]
    theirs = medians(duck, lambda r: ("duckdb", duckdb_name(r), r["workers"]))
    return {**mine, **theirs}


def table(
    m: Medians,
    rows: list[tuple[str, str]],
    configs: list[str],
    workers: int = 8,
    field: str = "query_seconds",
) -> str:
    lines = [
        "| Query | Cells | " + " | ".join(configs) + " |",
        "|---|---|" + "---|" * len(configs),
    ]
    for name, label in rows:
        found = [m[(c, name, workers)] for c in configs if (c, name, workers) in m]
        if not found:
            continue
        cells = int(found[0]["cells"])
        cols = []
        for c in configs:
            v = m.get((c, name, workers), {}).get(field)
            cols.append(fmt(v))
        lines.append(f"| {label} | {cells:,} | " + " | ".join(cols) + " |")
    return "\n".join(lines)


def cost_table(m: Medians, rows: list[tuple[str, str]], configs: list[str]) -> str:
    lines = [
        "| Query | " + " | ".join(configs) + " |",
        "|---|" + "---|" * len(configs),
    ]
    for name, label in rows:
        cols = []
        for c in configs:
            v = m.get((c, name, 8))
            if v is None:
                cols.append("–")
                continue
            cols.append(
                f"{fmt(v.get('requests'), 0)} / {fmt(v.get('cpu_seconds'), 1)} / "
                f"{fmt(v.get('peak_rss_mb', 0) / 1024, 1)}"
            )
        lines.append(f"| {label} | " + " | ".join(cols) + " |")
    return "\n".join(lines)


def workers_table(m: Medians, configs: list[str]) -> str:
    lines = [
        "| Query | Workers | " + " | ".join(configs) + " |",
        "|---|---|" + "---|" * len(configs),
    ]
    for name, label in [("fp_none", "(d1)"), ("neutrophil", "(b)")]:
        for w in (1, 8, 32):
            cols = [fmt(m.get((c, name, w), {}).get("query_seconds")) for c in configs]
            lines.append(f"| {label} | {w} | " + " | ".join(cols) + " |")
    return "\n".join(lines)


def per_store(stores: dict[str, list[dict[str, Any]]], f: Any) -> list[float]:
    """For each store, the median of `f` over its repeats."""
    return [statistics.median(f(r) for r in v) for v in stores.values()]


def bind_table(rows: list[dict[str, Any]]) -> str:
    by: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for r in rows:
        by[(r["layout"], r["query"])].append(r)
    header = (
        "| Stores | Query | Columns | `read_elem(obs)` ms | Builder ms: read + eval | "
        "of it eval | CPU ms: `read_elem` / builder | Requests: `read_elem` / builder | "
        "Sum over stores, s |"
    )
    lines = [header, "|---|---|---|---|---|---|---|---|---|"]
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
        med = statistics.median
        lines.append(
            f"| {len(stores)} {layout.replace('_', '-')} | {query} | "
            f"{len(rs[0]['columns'])} | {fmt(med(main), 1)} | {fmt(med(bind), 1)} | "
            f"{fmt(med(ev), 2)} | {fmt(med(main_cpu), 1)} / {fmt(med(bind_cpu), 1)} | "
            f"{fmt(med(main_req), 0)} / {fmt(med(bind_req), 0)} | "
            f"{fmt(sum(main) / 1e3, 1)} / {fmt(sum(bind) / 1e3, 2)} |"
        )
    return "\n".join(lines)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("results", type=pathlib.Path)
    args = p.parse_args()
    queries = args.results / "queries-duckdb.json"
    if queries.exists():
        for name, form in json.loads(queries.read_text()).items():
            SQL.setdefault(str(form["sql"]), []).append(name)
    m = all_medians(args.results)
    sel = [
        "main",
        "main_lean",
        "mapping_pd",
        "builder_pd",
        "builder_cons",
        "builder_ice",
        "duckdb",
    ]
    obs = ["builder_pd", "builder_cons", "builder_ice", "duckdb"]
    print("Selections with X, seconds at 8 workers:\n")
    print(table(m, SELECTIONS, sel))
    print("\nobs-only questions, seconds at 8 workers:\n")
    print(table(m, OBS_ONLY, obs))
    print("\nRequests / CPU s / peak RSS GiB at 8 workers:\n")
    rows = [r for r in SELECTIONS + OBS_ONLY if r[0] in {
        "one_dataset", "neutrophil", "lung_macrophage", "fp_none", "fp_some",
        "membership", "count_b_donor", "count_all_donor", "wmg_count", "wmg", "all_obs",
    }]  # fmt: skip
    print(
        cost_table(
            m, rows, ["main", "builder_pd", "builder_cons", "builder_ice", "duckdb"]
        )
    )
    print("\nWorkers, seconds:\n")
    print(
        workers_table(
            m, ["main", "builder_pd", "builder_cons", "builder_ice", "duckdb"]
        )
    )
    bind = load(args.results / "bind-ec2.jsonl")
    if bind:
        print("\nBind and evaluate per store, against `read_elem(obs)`:\n")
        print(bind_table(bind))
    laptop = load(args.results / "runs-laptop.jsonl")
    if laptop:
        ml = medians(laptop, lambda r: (r["config"], r["query"], r["workers"]))
        print("\nLaptop, seconds at 8 workers, median of 2:\n")
        rows = [r for r in SELECTIONS + OBS_ONLY if any(k[1] == r[0] for k in ml)]
        print(table(ml, rows, ["main", "builder_pd", "builder_cons"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
