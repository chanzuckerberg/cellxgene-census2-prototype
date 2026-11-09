"""Does pyarrow evaluate `==` and `is_in` on a dictionary array over its dictionary?

The script times, on one thread (`pyarrow.set_cpu_count(1)`), best of 7, for a
dictionary array of 2 million cells:

- the decode of the dictionary array to strings (`cast`)
- `equal` and `is_in` on the dictionary array, and on the decoded strings
- `starts_with` on the dictionary array (pyarrow 26 has no kernel) and on the strings
- the same terms through an expression, with `Table.filter` and with an Acero
  `project` node
- the evaluation by codes: the term on the dictionary, then `take` by the indices,
  or a NumPy lookup table indexed by the codes

When `equal` on the dictionary array costs about the decode plus `equal` on the
strings, pyarrow decodes. One JSON line is printed for each measure.

    python dictionary.py --out dictionary.jsonl
"""

import argparse
import json
import pathlib
import platform
import time
from collections.abc import Callable
from typing import Any

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
from pyarrow import acero

N = 2_000_000


def best(fn: Callable[[], Any], reps: int = 7) -> float:
    out = float("inf")
    for _ in range(reps):
        t0 = time.perf_counter()
        fn()
        out = min(out, time.perf_counter() - t0)
    return round(out * 1e3, 3)


def project(table: Any, expression: Any) -> Any:
    plan = acero.Declaration.from_sequence(
        [
            acero.Declaration("table_source", acero.TableSourceNodeOptions(table)),
            acero.Declaration("project", acero.ProjectNodeOptions([expression], ["m"])),
        ]
    )
    return plan.to_table(use_threads=False)["m"]


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--out", type=pathlib.Path, default=None)
    p.add_argument("--where", default="laptop")
    args = p.parse_args()
    pa.set_cpu_count(1)
    pa.set_io_thread_count(1)
    rng = np.random.default_rng(0)
    records = []
    for n_categories in (300, 30_000):
        labels = [f"cell type number {i:06d}" for i in range(n_categories)]
        categories = pa.array(labels)
        codes = rng.integers(0, n_categories, N).astype(np.int32)
        dictionary = pa.DictionaryArray.from_arrays(pa.array(codes), categories)
        strings = dictionary.cast(pa.string())
        want = labels[42]
        values = pa.array([labels[42], labels[100], labels[200]])
        table = pa.table({"c": dictionary})
        table_strings = pa.table({"c": strings})
        lut = np.zeros(n_categories + 1, dtype=bool)
        lut[[42, 100, 200]] = True

        def starts_with_dictionary() -> Any:
            return pc.starts_with(dictionary, "cell type number 00")

        try:
            starts_with_dictionary()
            dictionary_kernel = True
        except pa.ArrowNotImplementedError:
            dictionary_kernel = False
        measures: dict[str, Callable[[], Any]] = {
            "decode": lambda: dictionary.cast(pa.string()),
            "equal, dictionary": lambda: pc.equal(dictionary, want),
            "equal, strings": lambda: pc.equal(strings, want),
            "is_in, dictionary": lambda: pc.is_in(dictionary, value_set=values),
            "is_in, strings": lambda: pc.is_in(strings, value_set=values),
            "starts_with, strings": lambda: pc.starts_with(
                strings, "cell type number 00"
            ),
            "Table.filter(==), dictionary": lambda: table.filter(pc.field("c") == want),
            "Table.filter(==), strings": lambda: table_strings.filter(
                pc.field("c") == want
            ),
            "Acero project(==), dictionary": lambda: project(
                table, pc.field("c") == want
            ),
            "Acero project(isin), dictionary": lambda: project(
                table, pc.field("c").isin(values)
            ),
            "by codes: is_in on the dictionary, take": lambda: pc.take(
                pc.is_in(categories, value_set=values), dictionary.indices
            ),
            "by codes: NumPy lookup table": lambda: lut[codes],
        }
        for name, fn in measures.items():
            record = {
                "categories": n_categories,
                "cells": N,
                "measure": name,
                "ms": best(fn),
                "starts_with_dictionary_kernel": dictionary_kernel,
                "pyarrow": pa.__version__,
                "machine": platform.machine(),
                "where": args.where,
            }
            records.append(record)
            print(json.dumps(record), flush=True)
    if args.out is not None:
        args.out.write_text("".join(json.dumps(r) + "\n" for r in records))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
