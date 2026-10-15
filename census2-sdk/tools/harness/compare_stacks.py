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
restricts both sides to the shared genes. The genes that only the SDK holds are an
expected finding. A Census 1 gene that the SDK does not hold fails the query.
"""

import argparse
import collections
import hashlib
import json
import os
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
    (r"^X\.normalized$", "4.5", "the normalized layer that the writer computes"),
    (
        r"^var\.(nnz|n_measured_obs)$",
        "4.5",
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
# The SDK's var is the whole pinned reference, so its extra genes do not block a query.
EXPECTED = "var genes, SDK only"


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


def _finding_record(finding: Finding, row_columns: set[str]) -> dict[str, Any]:
    who, why = owner(finding, row_columns)
    return {
        "kind": finding.kind,
        "where": finding.where,
        "detail": finding.detail,
        "owner": who,
        "reason": why,
    }


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
        record.update(status="unsupported", reason=str(reason))
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
        (EXPECTED, genes_sdk - genes_c1),
        ("var genes, Census 1 only", genes_c1 - genes_sdk),
    ):
        if genes:
            detail = f"{len(genes)} genes; X compared on {len(shared)} shared genes"
            findings.append(Finding("row-set", where, detail))
    if genes_c1 != genes_sdk:
        census1, sdk = restrict_genes(census1, shared), restrict_genes(sdk, shared)
    findings += compare(census1, sdk)

    row_columns = _row_columns(query, pins)
    record["findings"] = [_finding_record(f, row_columns) for f in findings]
    blocking = [f for f in record["findings"] if f["where"] != EXPECTED]
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
            by_owner[f["owner"]][f"{f['kind']} {f['where']}"] += 1
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
        print(line)
        for f in record.get("findings", [])[:5]:
            print(
                f"    [{f['kind']}] {f['where']}: {f['detail']}  (owner {f['owner']})"
            )
    uris = sorted(census.catalog["store_uri"].astype(str))
    catalog = {
        "file": pathlib.Path(args.catalog).name,
        "sha256": None
        if census2_sdk.catalog.is_remote(args.catalog)
        else hashlib.sha256(pathlib.Path(args.catalog).read_bytes()).hexdigest(),
        "stores": len(uris),
        "store_prefix": os.path.commonprefix(uris),
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
