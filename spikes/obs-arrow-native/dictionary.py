"""Time the builder's pyarrow terms on a dictionary array, against a lookup by codes.

One column of 2 million cells, with 300 or 30,000 categories of about 30 characters,
a tenth of them joined values of two members (`a || b`). One thread
(`pyarrow.set_cpu_count(1)`), best of 7 repeats, in milliseconds. The measures:

- `decode`: the cast of the dictionary array to strings, one for each cell.
- `==`, `isin`: the builder's expression for each, on the dictionary array.
  `isin` is `is_in` in `if_else(is_valid(x), ..., null)`, for the three-valued logic.
- `has_any`: the builder's anchored RE2 expression over the cast strings.
- `project(isin & ==)`: an Acero projection of a query of two terms over two columns,
  as `obs_read.evaluate` runs it.
- `lut ==`, `lut isin`, `lut has_any`: the earlier filter spike's evaluation: a test
  of each category in Python, then a NumPy lookup by the codes.
- `dict has_any`: `match_substring_regex` on the dictionary's values, then `take` by
  the indices: what a dictionary kernel of pyarrow could do.

    python dictionary.py --where ec2 --out dictionary.jsonl
"""

import argparse
import json
import pathlib
import platform
import re
import sys
import time
from collections.abc import Callable
from typing import Any

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
from pyarrow import acero

N = 2_000_000


def best(fn: Callable[[], Any], repeats: int = 7) -> float:
    times = []
    for _ in range(repeats):
        t0 = time.perf_counter()
        fn()
        times.append(time.perf_counter() - t0)
    return round(min(times) * 1e3, 2)


def labels(k: int, rng: np.random.Generator) -> list[str]:
    words = [f"term{i:05d} of a tissue" for i in range(k)]
    out = []
    for i, w in enumerate(words):
        out.append(f"{w} || {words[(i * 7 + 3) % k]}" if i % 10 == 0 else w)
    return sorted(set(out))


def measure(k: int) -> dict[str, float]:
    sys.path.insert(0, str(pathlib.Path(__file__).parents[2] / "census2-sdk" / "src"))
    from census2_sdk import obs_read
    from census2_sdk.query import members_pattern

    rng = np.random.default_rng(0)
    cats = labels(k, rng)
    codes = rng.integers(0, len(cats), N).astype(np.int16 if k > 127 else np.int8)
    col = obs_read.Codes(codes, cats, False)
    other = obs_read.Codes(
        rng.integers(0, 3, N).astype(np.int8), ["a", "b", "c"], False
    )
    array = col.arrow()
    picked = cats[: max(1, len(cats) // 100)]
    one = cats[len(cats) // 2]
    member = "term00003 of a tissue"
    pattern = members_pattern([member])
    f = pc.field("x")
    isin = pc.if_else(
        f.is_valid(), f.isin(pa.array(picked)), pa.scalar(None, pa.bool_())
    )
    table = pa.table({"x": array, "y": other.arrow()})
    query = isin & (pc.field("y") == "b")

    def project(expression: Any) -> Any:
        return acero.Declaration.from_sequence(
            [
                acero.Declaration("table_source", acero.TableSourceNodeOptions(table)),
                acero.Declaration(
                    "project", acero.ProjectNodeOptions([expression], ["m"])
                ),
            ]
        ).to_table(use_threads=False)

    rx = re.compile(r"(?:^| \|\| )\s*(?:" + re.escape(member) + r")\s*(?: \|\| |$)")
    wanted = set(picked)

    def lut(test: Callable[[str], bool]) -> np.ndarray:
        table_ = np.fromiter((test(c) for c in cats), bool, len(cats))
        return np.append(table_, False)[codes]

    return {
        "decode": best(lambda: array.cast(pa.string())),
        "==": best(lambda: pc.equal(array, pa.scalar(one))),
        "isin": best(lambda: project(isin)),
        "has_any": best(
            lambda: pc.match_substring_regex(array.cast(pa.string()), pattern=pattern)
        ),
        "project(isin & ==)": best(lambda: project(query)),
        "lut ==": best(lambda: lut(lambda c: c == one)),
        "lut isin": best(lambda: lut(lambda c: c in wanted)),
        "lut has_any": best(lambda: lut(lambda c: rx.search(c) is not None)),
        "dict has_any": best(
            lambda: pc.match_substring_regex(array.dictionary, pattern=pattern).take(
                array.indices
            )
        ),
    }


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--where", required=True)
    p.add_argument("--out", type=pathlib.Path, required=True)
    args = p.parse_args()
    pa.set_cpu_count(1)
    with args.out.open("a") as fh:
        for k in (300, 30_000):
            for name, ms in measure(k).items():
                record = {
                    "where": args.where,
                    "machine": platform.machine(),
                    "pyarrow": pa.__version__,
                    "cells": N,
                    "categories": k,
                    "measure": name,
                    "ms": ms,
                }
                fh.write(json.dumps(record) + "\n")
                print(f"{k:>6} {name:<20} {ms:>8}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
