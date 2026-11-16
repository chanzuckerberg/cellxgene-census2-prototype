"""Run the end-to-end `get_anndata` queries over the SDK of one branch, or a matrix.

`run` runs one named query of `bench/filters.py` in this process, over the SDK that
this interpreter imports, and prints one JSON record. Three query APIs:

- `mapping`: `census.get_anndata("homo_sapiens", mapping)`, as on `main`.
- `pr43`: `census.homo_sapiens.get_anndata(filter)`, the builder of the earlier
  filter spike (`spike/consolidation-filters`).
- `spike`: `census.get_anndata(query)`, this branch.

The record holds the wall time and the CPU time of the query (user and system, all
threads), each read at the zarr store interface and its bytes, the peak RSS of the
process, the cells, and the stores that the query opened. The count wraps
`zarr.storage.ObjectStore`, which each branch uses for an S3 store, so it counts the
reads of every branch the same way. Each `get` is one request; a `get_partial_values`
counts one request for each range. A watchdog stops a run whose resident memory
passes `--max-rss-gb`, with the status `memory_limit`, before the host runs out.

`matrix` runs `run` in a new process for each arm, query, worker setting, and repeat,
with the arms interleaved in each repeat, and appends each record to `--out`. It skips
each run that `--out` already holds, so a matrix that stopped can start again.

    python e2e.py matrix --catalog catalog-103.parquet --tools census2-sdk/tools \\
        --arm main=venv-main/bin/python:mapping --arm spike=venv-spike/bin/python:spike \\
        --queries liver,neutrophil --workers 8,32 --repeats 3 --where ec2 --out e2e.jsonl
"""

import argparse
import json
import os
import pathlib
import resource
import subprocess
import sys
import threading
import time
from typing import Any, Generic, TypeVar

import zarr.storage


class Counter:
    def __init__(self) -> None:
        self.lock = threading.Lock()
        self.requests = 0
        self.nbytes = 0

    def add(self, buffers: list[Any]) -> None:
        with self.lock:
            self.requests += len(buffers)
            self.nbytes += sum(len(b) for b in buffers if b is not None)


COUNTER = Counter()
T = TypeVar("T")


class CountingObjectStore(zarr.storage.ObjectStore[T], Generic[T]):
    """`zarr.storage.ObjectStore`, generic as it is, because a branch subscripts it."""

    async def get(self, key: str, prototype: Any, byte_range: Any = None) -> Any:
        buffer = await super().get(key, prototype, byte_range)
        COUNTER.add([buffer])
        return buffer

    async def get_partial_values(self, prototype: Any, key_ranges: Any) -> Any:
        buffers = await super().get_partial_values(prototype, key_ranges)
        COUNTER.add(list(buffers))
        return buffers


def usage() -> tuple[float, float]:
    u = resource.getrusage(resource.RUSAGE_SELF)
    return u.ru_utime, u.ru_stime


def rss_mb() -> float:
    """The resident memory of this process now, from /proc on Linux."""
    try:
        for line in pathlib.Path("/proc/self/status").read_text().splitlines():
            if line.startswith("VmRSS:"):
                return int(line.split()[1]) / 1024
    except OSError:
        pass
    return 0.0


def guard_memory(limit_mb: float, record: dict[str, Any]) -> None:
    """Stop this process when its memory passes `limit_mb`, before the host runs out."""

    def watch() -> None:
        while True:
            now = rss_mb()
            if now > limit_mb:
                over = {**record, "status": "memory_limit", "rss_mb": round(now, 1)}
                print(json.dumps(over), flush=True)
                os._exit(3)
            time.sleep(0.2)

    threading.Thread(target=watch, daemon=True).start()


def peak_rss_mb() -> float:
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return round(peak / (2**20 if sys.platform == "darwin" else 2**10), 1)


def run(args: argparse.Namespace) -> dict[str, Any]:
    zarr.storage.ObjectStore = CountingObjectStore  # type: ignore[misc]
    guard_memory(args.max_rss_gb * 1024, {"query": args.query, "api": args.api})
    sys.path.insert(0, str(args.tools))
    import census2_sdk
    from bench.filters import QUERIES

    query = QUERIES[args.query]
    census = census2_sdk.open_census([args.catalog])
    record: dict[str, Any] = {"query": args.query, "api": args.api}
    user, system = usage()
    t0 = time.monotonic()
    if args.api == "mapping":
        if query.mapping is None:
            raise SystemExit(f"{args.query} has no mapping form")
        adata = census.get_anndata(
            "homo_sapiens", query.mapping, max_workers=args.workers
        )
        opened = len(
            census2_sdk.catalog.exclude(census.catalog, "homo_sapiens", query.mapping)
        )
        whole = None
    elif args.api == "pr43":
        hs = census.homo_sapiens
        adata = hs.get_anndata(query.build(hs.obs), max_workers=args.workers)
        opened, whole = hs.last_query.get("stores_opened"), None
    else:
        hs = census.homo_sapiens
        built = query.build(hs.obs)
        q = hs if built is None else built
        adata = census.get_anndata(q, max_workers=args.workers)
        rows = census._candidates(q)
        opened, whole = len(rows), int(rows["whole"].sum())
    seconds = time.monotonic() - t0
    end_user, end_system = usage()
    record.update(
        workers=args.workers,
        cells=0 if adata is None else int(adata.n_obs),
        stores_opened=opened,
        stores_whole=whole,
        seconds=round(seconds, 3),
        cpu_seconds=round(end_user - user + end_system - system, 3),
        requests=COUNTER.requests,
        mb=round(COUNTER.nbytes / 2**20, 2),
        peak_rss_mb=peak_rss_mb(),
        status="ok",
    )
    return record


def matrix(args: argparse.Namespace) -> None:
    sys.path.insert(0, str(args.tools))
    from bench.filters import QUERIES

    arms = []
    for spec in args.arm:
        name, rest = spec.split("=", 1)
        python, api = rest.rsplit(":", 1)
        arms.append((name, python, api))
    names = args.queries.split(",")
    script = pathlib.Path(__file__).resolve()
    # A matrix that stopped starts again where it stopped: skip each finished run.
    finished = set()
    if args.out.exists():
        for line in args.out.read_text().splitlines():
            r = json.loads(line)
            if r.get("status") == "ok":
                finished.add((r["query"], r["arm"], r["workers"], r["repeat"]))
    for repeat in range(args.repeats):
        for workers in (int(w) for w in args.workers.split(",")):
            for name in names:
                for arm, python, api in arms:
                    if api == "mapping" and QUERIES[name].mapping is None:
                        continue
                    if (name, arm, workers, repeat) in finished:
                        continue
                    command = [
                        python, str(script), "run", "--api", api,
                        "--catalog", str(args.catalog), "--tools", str(args.tools),
                        "--query", name, "--workers", str(workers),
                    ]  # fmt: skip
                    try:
                        done = subprocess.run(
                            command,
                            check=False,
                            capture_output=True,
                            text=True,
                            timeout=args.timeout,
                        )
                        lines = done.stdout.strip().splitlines()
                        if lines and '"memory_limit"' in lines[-1]:
                            record = json.loads(lines[-1])
                        elif done.returncode or not lines:
                            record = {
                                "query": name,
                                "status": "error",
                                "stderr": done.stderr[-2000:],
                            }
                        else:
                            record = json.loads(lines[-1])
                    except subprocess.TimeoutExpired:
                        record = {"query": name, "status": "timeout"}
                    record.update(
                        arm=arm, repeat=repeat, workers=workers, where=args.where
                    )
                    with args.out.open("a") as fh:
                        fh.write(json.dumps(record) + "\n")
                    print(
                        f"r{repeat} w{workers:<3} {name:<14} {arm:<6} "
                        f"{record.get('seconds', record['status'])}",
                        flush=True,
                    )


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = p.add_subparsers(dest="command", required=True)
    r = sub.add_parser("run")
    r.add_argument("--api", choices=["mapping", "pr43", "spike"], required=True)
    r.add_argument("--catalog", required=True)
    r.add_argument("--tools", type=pathlib.Path, required=True)
    r.add_argument("--query", required=True)
    r.add_argument("--workers", type=int, default=8)
    r.add_argument("--max-rss-gb", type=float, default=48)
    m = sub.add_parser("matrix")
    m.add_argument("--catalog", required=True)
    m.add_argument("--tools", type=pathlib.Path, required=True)
    m.add_argument("--arm", action="append", required=True, help="name=python:api")
    m.add_argument("--queries", required=True)
    m.add_argument("--workers", default="8,32")
    m.add_argument("--repeats", type=int, default=3)
    m.add_argument("--timeout", type=float, default=1800)
    m.add_argument("--where", required=True)
    m.add_argument("--out", type=pathlib.Path, required=True)
    args = p.parse_args()
    if args.command == "run":
        print(json.dumps(run(args)), flush=True)
    else:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        matrix(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
