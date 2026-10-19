"""The oracle runner (task 2.2): record the answers of the pinned Census 1 stack.

The runner executes each corpus query through the Census 1 SDK surface (`get_anndata`,
`get_obs`, `get_var`, `get_presence_matrix`) against the pinned release, and it
serializes each result in the canonical form (`canonical.py`). Two runs of one query
against one release give the same `result_digest`. `--verify` checks that.

The top-level SDK accessors do not reach the spatial experiments. A spatial
`get_anndata` query runs through `experiment.axis_query` on the same census handle,
the code path `get_anndata` reaches internally (evidence 11.2). Spatial `obs`, `var`,
and presence queries read the SOMA objects directly. Each result records which path ran.

The snapshot is bounded. An obs-filtered query is clamped to the N smallest datasets
of each experiment (`snapshot_manifest.json`, built once from each experiment's obs).
Differential validation is meaningful at one to five datasets (evidence 11.1), and a
full-corpus X read belongs to the benchmark rig. `var` and presence queries are small
and run unclamped.

Run from `census2-sdk/tools/` in the oracle environment (see README.md):

    python -m harness.oracle_runner build-manifest   # once; writes snapshot_manifest.json
    python -m harness.oracle_runner run [--only substr] [--out DIR]
    python -m harness.oracle_runner run --verify     # second pass, compares result digests
"""

import argparse
import json
import pathlib
import platform
import sys
import time

import cellxgene_census
import tiledbsoma as soma

from harness.canonical import (
    canonicalize_anndata,
    frame_arrays,
    layer_arrays,
    result_digest,
    write_result,
)
from harness.corpus import (
    CORPUS_VERSION,
    clamp_filter,
    expand,
    load_values,
    parse_picks,
)

HERE = pathlib.Path(__file__).parent
MANIFEST = HERE / "snapshot_manifest.json"
DATASETS_PER_EXPERIMENT = 3
ORGANISM_LABEL = {
    "homo_sapiens": "Homo sapiens",
    "mus_musculus": "Mus musculus",
    "macaca_mulatta": "Macaca mulatta",
    "callithrix_jacchus": "Callithrix jacchus",
    "pan_troglodytes": "Pan troglodytes",
}


def open_census(values):
    return cellxgene_census.open_soma(uri=values["_source"]["uri"])


def env_record():
    import anndata
    import numpy

    return {
        "cellxgene_census": cellxgene_census.__version__,
        "tiledbsoma": soma.__version__,
        "anndata": anndata.__version__,
        "numpy": numpy.__version__,
        "python": platform.python_version(),
        "corpus_version": CORPUS_VERSION,
    }


def build_manifest(values, picks=None):
    """The N smallest datasets of each experiment, from the experiment's own obs.

    `picks` overrides the sample of an experiment key (the v1.1 human sample, README).
    """
    picks = picks or {}
    manifest = {"datasets_per_experiment": DATASETS_PER_EXPERIMENT}
    with open_census(values) as census:
        for modality in ("census_data", "census_spatial_sequencing"):
            for experiment in census[modality]:
                exp = census[modality][experiment]
                col = (
                    exp.obs.read(column_names=["dataset_id"])
                    .concat()
                    .column("dataset_id")
                    .to_pandas()
                )
                counts = col.value_counts().sort_values(kind="stable")
                key = f"{modality}/{experiment}"
                if key in picks:
                    chosen = sorted(picks[key])
                    missing = [d for d in chosen if d not in counts.index]
                    if missing:
                        raise SystemExit(
                            f"{key}: picked datasets not in the release: {missing}"
                        )
                    cells = int(counts[chosen].sum())
                else:
                    chosen = sorted(
                        str(d) for d in counts.index[:DATASETS_PER_EXPERIMENT]
                    )
                    cells = int(counts.iloc[:DATASETS_PER_EXPERIMENT].sum())
                manifest[key] = {"dataset_ids": chosen, "cells": cells}
                print(f"{key}: {chosen} ({cells:,} cells)")
    unknown = sorted(set(picks) - set(manifest))
    if unknown:
        raise SystemExit(f"--pick names experiments the release lacks: {unknown}")
    MANIFEST.write_text(json.dumps(manifest, indent=1, sort_keys=True) + "\n")
    print(f"wrote {MANIFEST}")


def run_query(census, q, sample):
    """Execute one bound query. Return the canonical arrays and the path that ran."""
    label = ORGANISM_LABEL[q.experiment]
    clamped = q.key != "presence_matrix" and q.accessor != "get_var"
    obs_filter = (
        clamp_filter(q.obs_filter, sample["dataset_ids"]) if clamped else q.obs_filter
    )

    if q.accessor == "get_presence_matrix":
        if q.modality == "census_data":
            pm = cellxgene_census.get_presence_matrix(census, label, "RNA")
            path = "sdk"
        else:  # spatial: the same array, read directly
            exp = census[q.modality][q.experiment]
            pm = (
                exp.ms["RNA"]["feature_dataset_presence_matrix"]
                .read()
                .coos()
                .concat()
                .to_scipy()
            )
            path = "direct"
        return layer_arrays("presence", pm.tocsr()), path

    if q.accessor == "get_var":
        if q.modality == "census_data":
            var = cellxgene_census.get_var(census, label, value_filter=q.var_filter)
            path = "sdk"
        else:
            kw = {"value_filter": q.var_filter} if q.var_filter else {}
            var = (
                census[q.modality][q.experiment]
                .ms["RNA"]
                .var.read(**kw)
                .concat()
                .to_pandas()
            )
            path = "direct"
        var = var.sort_values("feature_id", kind="stable").reset_index(drop=True)
        return frame_arrays(
            "var", var.drop(columns=["soma_joinid"], errors="ignore")
        ), path

    if q.accessor == "get_obs":
        cols = list(dict.fromkeys(("dataset_id", "observation_joinid") + q.obs_columns))
        if q.modality == "census_data":
            obs = cellxgene_census.get_obs(
                census, label, value_filter=obs_filter, column_names=cols
            )
            path = "sdk"
        else:
            obs = (
                census[q.modality][q.experiment]
                .obs.read(value_filter=obs_filter, column_names=cols)
                .concat()
                .to_pandas()
            )
            path = "direct"
        order = obs.sort_values(
            ["dataset_id", "observation_joinid"],
            key=lambda s: s.astype(str),
            kind="stable",
        ).reset_index(drop=True)
        return frame_arrays(
            "obs", order.drop(columns=["soma_joinid"], errors="ignore")
        ), path

    first, *rest = q.layers
    if q.modality == "census_data":
        adata = cellxgene_census.get_anndata(
            census,
            label,
            X_name=first,
            X_layers=rest,
            obs_value_filter=obs_filter,
            var_value_filter=q.var_filter,
        )
        path = "sdk"
    else:
        exp = census[q.modality][q.experiment]
        var_q = (
            soma.AxisQuery(value_filter=q.var_filter)
            if q.var_filter
            else soma.AxisQuery()
        )
        with exp.axis_query(
            "RNA", obs_query=soma.AxisQuery(value_filter=obs_filter), var_query=var_q
        ) as aq:
            adata = aq.to_anndata(X_name=first, X_layers=rest)
        path = "axis_query"
    adata.uns["x_layer_name"] = first
    return canonicalize_anndata(adata), path


def main():
    p = argparse.ArgumentParser()
    p.add_argument("command", choices=["build-manifest", "run"])
    p.add_argument("--only", default="")
    p.add_argument(
        "--pick",
        action="append",
        default=[],
        metavar="EXPERIMENT=ID,ID,...",
        help="build-manifest: use these datasets for one experiment",
    )
    p.add_argument("--out", type=pathlib.Path, default=HERE / "oracle_snapshot")
    p.add_argument(
        "--verify",
        action="store_true",
        help="compare digests against an existing snapshot",
    )
    args = p.parse_args()
    if args.pick and args.command != "build-manifest":
        p.error("--pick applies to build-manifest only")
    try:
        picks = parse_picks(args.pick, DATASETS_PER_EXPERIMENT)
    except ValueError as error:
        p.error(str(error))

    values = load_values()
    if args.command == "build-manifest":
        build_manifest(values, picks)
        return

    manifest = json.loads(MANIFEST.read_text())
    queries = [q for q in expand(values) if args.only in q.label]
    env = env_record()
    print(
        f"oracle run: corpus v{CORPUS_VERSION}, release {values['_source']['release']}, "
        f"{len(queries)} queries, sdk {env['cellxgene_census']}"
    )

    mismatches = 0
    with open_census(values) as census:
        for q in queries:
            sample = manifest[f"{q.modality}/{q.experiment}"]
            t0 = time.monotonic()
            arrays, path = run_query(census, q, sample)
            out_dir = args.out / q.label.replace("/", "__")
            if args.verify:
                old = json.loads((out_dir / "manifest.json").read_text())[
                    "result_digest"
                ]
                new = result_digest(arrays)
                same = "identical" if new == old else "DIGEST MISMATCH"
                mismatches += new != old
                print(f"  {same}  {q.label}  {time.monotonic() - t0:.0f}s")
            else:
                meta = {
                    "query": q.label,
                    "path": path,
                    "env": env,
                    "source": values["_source"],
                    "clamp": sample["dataset_ids"],
                }
                digest = write_result(out_dir, arrays, meta)
                print(f"  {digest[:12]}  {q.label}  {time.monotonic() - t0:.0f}s")
    if args.verify:
        print(f"{'PASS' if not mismatches else 'FAIL'}: {mismatches} digest mismatches")
        sys.exit(1 if mismatches else 0)


if __name__ == "__main__":
    main()
