"""Run the query corpus through the Census 2 SDK and compare with Census 1 (task 6.2).

For each bound query, the driver loads the recorded Census 1 result, answers the query
through `sdk_adapter`, and runs the comparator. It prints one status line for
each query (the harness README lists the five statuses). It writes every finding to a
JSON record, with the task that owns it, and exits 1 on a failure.

Run from `census2-sdk/tools/` with the stores' credentials in the environment:

    python -m harness.compare_stacks CATALOG --oracle harness/oracle_snapshot \\
        --out results/sdk_vs_census1.json

The SDK's `var` is the whole pinned reference, and the Census 1 `var` is the genes of
its build. The comparator compares no X value while the gene sets differ, so the driver
restricts both sides to the shared genes. A Census 1 gene that the SDK does not hold
fails the query.

A finding in the declared list of expected differences (`EXPECTED` and `FLOAT32_SUMS`,
task 4.9) stays on record and does not fail the query. Census 1 is the left side.
"""

import argparse
import collections
import hashlib
import json
import pathlib
import re
import sys
import time
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse

import census2_sdk
from harness import sdk_adapter
from harness.canonical import (
    Arrays,
    decode_strings,
    frame_arrays,
    layer_arrays,
    load_result,
)
from harness.comparator import Finding, compare
from harness.corpus import BoundQuery, expand, load_values

HERE = pathlib.Path(__file__).parent
MANIFEST = HERE / "snapshot_manifest.json"

# The task that owns each kind of finding: (where pattern, task, reason). First match wins.
OWNERS: tuple[tuple[str, str, str], ...] = (
    (
        r"^obs\.(raw_sum|nnz|raw_mean_nnz|raw_variance_nnz|n_measured_vars)$",
        "4.5",
        "per-cell statistics that the writer computes",
    ),
    (r"^obs\.tissue_general", "4.5", "tissue_general, which the writer derives"),
    (r"^X\.normalized$", "4.13", "the normalized layer, which 4.13 decides"),
    (
        r"^var\.(nnz|n_measured_obs)$",
        "4.14",
        "corpus aggregates that the release builder sums from the per-dataset files",
    ),
    (
        r"^var\.feature_version$",
        "4.2",
        "a column of the pinned reference that Census 1 drops",
    ),
    (
        r"^var genes, SDK only$",
        "4.2",
        "genes of the pinned reference that the Census 1 build does not hold",
    ),
    (r"^obs\.", "4.4", "the obs columns that a store keeps"),
)
COLUMN_RULES = {r"^obs\."}  # these rules apply to column-set findings only
SDK_ONLY_GENES = "var genes, SDK only"
# The task that adds what an unsupported query needs: (reason pattern, task).
UNSUPPORTED_OWNERS: tuple[tuple[str, str], ...] = (
    (r"has no (get_obs|get_var|get_presence_matrix) accessor", "6.7"),
    (r"takes equality and membership only", "6.7"),
)

# The declared expected differences: (kind, where pattern, detail pattern, reason).
EXPECTED: tuple[tuple[str, str, str, str], ...] = (
    (
        "row-set",
        rf"^{SDK_ONLY_GENES}$",
        "",
        "the SDK's var is the whole pinned reference",
    ),
    (
        "column-set",
        r"^var\.feature_version$",
        r"^present on right only$",
        "the SDK's var keeps this column of the pinned reference",
    ),
    (
        "column-set",
        r"^var\.(nnz|n_measured_obs)$",
        r"^present on left only$",
        "the release builder adds the corpus aggregates in Iteration 3 (4.14)",
    ),
)
# Census 1 asks scipy for a float64 sum, but scipy < 1.18 adds float32 counts in float32.
# Census 1 thus rounds a row sum above 2**24 (largest relative difference 1.7e-7). The
# ingest round-trip check verifies the SDK sum (README).
FLOAT32_SUMS = {"obs.raw_sum", "obs.raw_mean_nnz"}
FLOAT32_RTOL = 1e-6


def _differ(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    """The cells where two numeric columns differ. NaN equals NaN."""
    left, right = left.astype(np.float64), right.astype(np.float64)
    same = (left == right) | (np.isnan(left) & np.isnan(right))
    return np.asarray(~same, dtype=bool)


def _float32_rounding(census1: Arrays, sdk: Arrays, column: str) -> bool:
    """True if `column` differs only where Census 1 rounds the row sum (harness README)."""
    needed = {"obs.raw_sum"} | ({"obs.nnz"} if column == "obs.raw_mean_nnz" else set())
    if not needed <= census1.keys() & sdk.keys():
        return False
    differ = _differ(census1[column], sdk[column])
    left, right = (side[column].astype(np.float64)[differ] for side in (census1, sdk))
    sum1, sum2 = (
        side["obs.raw_sum"].astype(np.float64)[differ] for side in (census1, sdk)
    )
    ok = (
        np.all(_differ(sum1, sum2))
        and np.all(sum2 > 2**24)
        and np.allclose(left, right, rtol=FLOAT32_RTOL, atol=0)
    )
    if ok and column == "obs.raw_mean_nnz":
        nnz1, nnz2 = (
            side["obs.nnz"].astype(np.float64)[differ] for side in (census1, sdk)
        )
        ok = (
            np.array_equal(nnz1, nnz2)
            and np.allclose(left, sum1 / nnz1, rtol=1e-12, atol=0)
            and np.allclose(right, sum2 / nnz2, rtol=1e-12, atol=0)
        )
    return bool(ok)


def expected(finding: Finding, census1: Arrays, sdk: Arrays) -> str | None:
    """The reason that `finding` is an expected difference, or None."""
    for kind, where, detail, reason in EXPECTED:
        if (
            finding.kind == kind
            and re.search(where, finding.where)
            and re.search(detail, finding.detail)
        ):
            return reason
    if (
        finding.kind == "value"
        and finding.where in FLOAT32_SUMS
        and _float32_rounding(census1, sdk, finding.where)
    ):
        return "Census 1 rounds a float32 row sum above 2**24"
    return None


def owner(finding: Finding, row_columns: set[str]) -> tuple[str, str]:
    """The task that owns `finding`. A row-set finding takes the owner of its filter."""
    where = finding.where
    if finding.kind in {"row-set", "strict-subset"} and where == "obs":
        for column in sorted(row_columns):
            if column.startswith("tissue_general"):
                return "4.5", "the filter names tissue_general, which the stores lack"
        return "6.2", "row sets differ; triage"
    for pattern, task, reason in OWNERS:
        if pattern in COLUMN_RULES and finding.kind != "column-set":
            continue
        if re.search(pattern, where):
            return task, reason
    return "6.2", "triage"


def drop_soma_joinid(arrays: Arrays) -> Arrays:
    """Census 1 rows carry `soma_joinid`, which a Census 2 result does not have."""
    return {k: v for k, v in arrays.items() if k.split(".")[1] != "soma_joinid"}


def _var_frame(arrays: Arrays) -> pd.DataFrame:
    columns: dict[str, Any] = {}
    for key, value in arrays.items():
        parts = key.split(".")
        if parts[0] != "var" or parts[-1] == "bytes":
            continue
        if parts[-1] == "offsets":
            name = ".".join(parts[1:-1])
            columns[name] = decode_strings(value, arrays[f"var.{name}.bytes"])
        else:
            columns[".".join(parts[1:])] = value
    return pd.DataFrame(columns)


def restrict_genes(arrays: Arrays, keep: set[str]) -> Arrays:
    """`arrays` with `var` and every layer restricted to the genes in `keep`."""
    var = _var_frame(arrays)
    mask = var["feature_id"].isin(keep).to_numpy()
    out = {k: v for k, v in arrays.items() if not k.startswith(("var.", "X."))}
    out.update(frame_arrays("var", var[mask].reset_index(drop=True)))
    for name in sorted({k.split(".")[1] for k in arrays if k.startswith("X.")}):
        rows, cols = (int(n) for n in arrays[f"X.{name}.shape"])
        matrix = sparse.csr_matrix(
            (
                arrays[f"X.{name}.data"],
                arrays[f"X.{name}.indices"],
                arrays[f"X.{name}.indptr"],
            ),
            shape=(rows, cols),
        )
        out.update(layer_arrays(name, matrix[:, np.flatnonzero(mask)].tocsr()))
    return out


def _obs_rows(arrays: Arrays) -> int:
    offsets = arrays.get("obs.dataset_id.offsets")
    return 0 if offsets is None else len(offsets) - 1


def _row_columns(query: BoundQuery, pins: dict[str, Any]) -> set[str]:
    """The obs columns that the filter names. `answer` accepted it, so this cannot raise."""
    if not query.obs_filter:
        return set()
    return set(sdk_adapter.translate(query.obs_filter, pins))


def _finding_record(
    finding: Finding, row_columns: set[str], expected: str | None = None
) -> dict[str, Any]:
    who, why = owner(finding, row_columns)
    record = {
        "kind": finding.kind,
        "where": finding.where,
        "detail": finding.detail,
        "owner": who,
        "reason": why,
    }
    if expected:
        record["expected"] = expected
    return record


def run_one(
    query: BoundQuery,
    oracle: pathlib.Path,
    census: census2_sdk.Census,
    stores: set[str],
    values: dict[str, Any],
    samples: dict[str, Any],
) -> dict[str, Any]:
    record: dict[str, Any] = {"query": query.label}
    sample = samples[f"{query.modality}/{query.experiment}"]["dataset_ids"]
    if not stores & set(sample):
        record.update(
            status="skip", reason="the catalog holds no store of the clamp sample"
        )
        return record
    result_dir = oracle / query.label.replace("/", "__")
    meta = json.loads((result_dir / "manifest.json").read_text())
    clamp = meta.get("clamp", [])
    if clamp and sorted(clamp) != sorted(sample):
        raise SystemExit(
            f"{query.label}: the recorded clamp differs from snapshot_manifest.json"
        )
    record["oracle"] = {"corpus_version": meta["env"]["corpus_version"], "clamp": clamp}
    pins = values[f"{query.modality}/{query.experiment}"]
    t0 = time.monotonic()
    try:
        sdk, notes = sdk_adapter.answer(query, census, clamp, pins)
    except sdk_adapter.Unsupported as reason:
        who = next(
            (task for p, task in UNSUPPORTED_OWNERS if re.search(p, str(reason))), "6.2"
        )
        record.update(status="unsupported", reason=str(reason), owner=who)
        return record
    except sdk_adapter.NoCells as reason:
        rows = _obs_rows(load_result(result_dir))
        found = Finding("row-set", "obs", f"{reason}; Census 1 selected {rows}")
        record.update(
            seconds=round(time.monotonic() - t0, 1),
            notes=[],
            cells={"census1": rows, "sdk": 0},
            findings=[_finding_record(found, _row_columns(query, pins))]
            if rows
            else [],
            status="fail" if rows else "vacuous",
        )
        return record
    record["seconds"] = round(time.monotonic() - t0, 1)
    record["notes"] = notes
    census1 = drop_soma_joinid(load_result(result_dir))
    record["cells"] = {"census1": _obs_rows(census1), "sdk": _obs_rows(sdk)}

    findings: list[Finding] = []
    genes_c1 = set(_var_frame(census1)["feature_id"])
    genes_sdk = set(_var_frame(sdk)["feature_id"])
    shared = genes_c1 & genes_sdk
    for where, genes in (
        (SDK_ONLY_GENES, genes_sdk - genes_c1),
        ("var genes, Census 1 only", genes_c1 - genes_sdk),
    ):
        if genes:
            detail = f"{len(genes)} genes; X compared on {len(shared)} shared genes"
            findings.append(Finding("row-set", where, detail))
    if genes_c1 != genes_sdk:
        census1, sdk = restrict_genes(census1, shared), restrict_genes(sdk, shared)
    findings += compare(census1, sdk)

    row_columns = _row_columns(query, pins)
    record["findings"] = [
        _finding_record(f, row_columns, expected(f, census1, sdk)) for f in findings
    ]
    blocking = [f for f in record["findings"] if "expected" not in f]
    if record["cells"]["census1"] == 0 and record["cells"]["sdk"] == 0:
        record["status"] = "vacuous"  # no cell to compare; the findings stay on record
    else:
        record["status"] = "fail" if blocking else "pass"
    return record


def summary(records: list[dict[str, Any]]) -> dict[str, Any]:
    status = collections.Counter(r["status"] for r in records)
    by_owner: dict[str, collections.Counter[str]] = collections.defaultdict(
        collections.Counter
    )
    for record in records:
        for f in record.get("findings", []):
            tag = " (expected)" if "expected" in f else ""
            by_owner[f["owner"]][f"{f['kind']} {f['where']}{tag}"] += 1
    return {
        "status": dict(status),
        "by_owner": {k: dict(sorted(v.items())) for k, v in sorted(by_owner.items())},
    }


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("catalog", help="the catalog of the stores under test")
    p.add_argument("--oracle", type=pathlib.Path, default=HERE / "oracle_snapshot")
    p.add_argument("--out", type=pathlib.Path, default=None)
    p.add_argument("--only", default="")
    args = p.parse_args()

    census = census2_sdk.open_census([args.catalog])
    stores = set(census.catalog["dataset_id"])
    values = load_values()
    samples = json.loads(MANIFEST.read_text())
    records = []
    for query in [q for q in expand(values) if args.only in q.label]:
        record = run_one(query, args.oracle, census, stores, values, samples)
        records.append(record)
        line = f"{record['status']:<11} {query.label}"
        if record["status"] in {"skip", "unsupported"}:
            line += f": {record['reason']}"
            if "owner" in record:
                line += f"  (owner {record['owner']})"
        print(line)
        for f in record.get("findings", [])[:5]:
            tag = ", expected" if "expected" in f else ""
            print(
                f"    [{f['kind']}] {f['where']}: {f['detail']}  (owner {f['owner']}{tag})"
            )
    catalog = {
        "file": pathlib.Path(args.catalog).name,
        "sha256": None
        if census2_sdk.catalog.is_remote(args.catalog)
        else hashlib.sha256(pathlib.Path(args.catalog).read_bytes()).hexdigest(),
        "stores": len(census.catalog),
    }
    report: dict[str, Any] = {
        "catalog": catalog,
        "summary": summary(records),
        "queries": records,
    }
    print(json.dumps(report["summary"], indent=1))
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(report, indent=1) + "\n")
    return 1 if report["summary"]["status"].get("fail") else 0


if __name__ == "__main__":
    sys.exit(main())
