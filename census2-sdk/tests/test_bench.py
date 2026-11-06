"""Tests for the store-count benchmark tools. No network, no AWS."""

import asyncio
import os
import threading
from pathlib import Path
from typing import Any

import anndata
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
import zarr.storage
from bench import catalogs, counting, matrix, run
from zarr.abc.store import RangeByteRequest
from zarr.core.buffer import default_buffer_prototype


def source_catalog() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "dataset_id": ["a", "b"],
            "store_uri": ["s3://bkt/source/v-a.zarr", "s3://bkt/source/v-b.zarr/"],
            "cell_type": [["T cell", "B cell"], ["T cell"]],
        }
    )


def test_expand_gives_unique_ids_and_copy_uris() -> None:
    rows = catalogs.expand(source_catalog(), "s3://bkt/stores/", 2)
    assert list(rows["dataset_id"]) == ["a", "b", "a-c01", "b-c01", "a-c02", "b-c02"]
    assert rows["dataset_id"].is_unique
    assert list(rows["copy"]) == [0, 0, 1, 1, 2, 2]
    assert list(rows["store_uri"][2:4]) == [
        "s3://bkt/stores/c01/v-a.zarr",
        "s3://bkt/stores/c01/v-b.zarr",
    ]
    # A smaller catalog is a prefix of a larger one.
    assert list(rows.head(3)["dataset_id"]) == ["a", "b", "a-c01"]


def test_pick_cell_type_is_nearest_share() -> None:
    rows = catalogs.expand(source_catalog(), "s3://bkt/stores", 1)
    assert run.pick_cell_type(rows, share=0.5) == "B cell"
    assert run.pick_cell_type(rows, share=1.0) == "T cell"


def test_cell_type_share_counts_the_source_rows_only() -> None:
    # An extra copy of "a" gives B cell a share of 4/7 over all rows. The source
    # rows give it a share of 1/2.
    rows = catalogs.expand(source_catalog(), "s3://bkt/stores", 2)
    uneven = pd.concat([rows, rows[rows["dataset_id"] == "a-c01"]], ignore_index=True)
    assert list(run.source_rows(uneven)["dataset_id"]) == ["a", "b"]
    assert run.pick_cell_type(run.source_rows(uneven), share=0.5) == "B cell"


def reset_counter() -> None:
    counting.COUNTER.__init__()  # type: ignore[misc]


def test_get_ranges_is_counted() -> None:
    store = counting.CountingStore(zarr.storage.MemoryStore())
    buffer = default_buffer_prototype().buffer.from_bytes(bytes(range(100)))

    async def ranges() -> list[Any]:
        await store.set("k", buffer)
        requests = [RangeByteRequest(0, 10), RangeByteRequest(50, 70)]
        prototype = default_buffer_prototype()
        return [g async for g in store.get_ranges("k", requests, prototype=prototype)]

    reset_counter()
    groups = asyncio.run(ranges())
    assert counting.COUNTER.requests["zarr"] == len(groups) >= 1
    assert counting.COUNTER.nbytes["zarr"] == 30


@pytest.fixture(scope="module")
def stores(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Two small stores in the writer's form, and a catalog over them."""
    root = tmp_path_factory.mktemp("bench")
    genes = [f"ENSG{j:011d}" for j in range(8)]
    rows = []
    for k, n_obs in enumerate((30, 20)):
        x = sp.random(n_obs, len(genes), density=0.3, format="csr", random_state=k)
        obs = pd.DataFrame(
            {
                "dataset_id": pd.Categorical([f"ds-{k}"] * n_obs),
                "observation_joinid": [f"j{i}" for i in range(n_obs)],
            },
            index=[str(i) for i in range(n_obs)],
        )
        obs["cell_type"] = pd.Categorical(["T cell"] * n_obs)
        obs["donor_id"] = pd.Categorical([f"donor-{k}"] * n_obs)
        var = pd.DataFrame({"feature_id": genes}, index=genes)
        uri = root / f"v-{k}.zarr"
        anndata.AnnData(X=x.astype(np.float32), obs=obs, var=var).write_zarr(uri)
        rows.append(
            {
                "dataset_id": f"ds-{k}",
                "store_uri": str(uri),
                "organism": "homo_sapiens",
                "spatial": False,
                "gene_reference_version": "GENCODE v48",
                "cell_type": ["T cell"],
            }
        )
    pd.DataFrame(rows).to_parquet(root / "catalog.parquet")
    return root


@pytest.mark.parametrize("shape", run.SHAPES)
@pytest.mark.parametrize("workers", [1, 2])
def test_every_shape_runs_and_counts(stores: Path, shape: str, workers: int) -> None:
    reset_counter()
    query = {"cell_type": "T cell"} if shape in ("query", "count") else None
    catalog = str(stores / "catalog.parquet")
    record = run.run(catalog, shape, workers, None, query=query)
    assert record["status"] == "ok"
    assert record["catalog_stores"] == 2
    expected = {"single": (1, 30), "open": (2, 0)}.get(shape, (2, 50))
    assert (record["stores_opened"], record["cells"]) == expected
    assert record["requests"] > 0
    assert record["mb_payload"] > 0
    assert record["cpu_seconds"] >= 0 and record["cpu_open_seconds"] >= 0


def test_memory_limit_stops_the_run(
    stores: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # The run is a new process, so it needs the import path of the tools and the SDK.
    here = Path(__file__).resolve().parents[1]
    monkeypatch.setenv("PYTHONPATH", f"{here / 'tools'}{os.pathsep}{here / 'src'}")
    record = matrix.one(
        stores / "catalog.parquet", "open", 1, 120, "T cell", max_rss_gb=0.001
    )
    assert record["status"] == "memory_limit"
    assert record["peak_rss_mb"] > record["limit_mb"]
    assert record["shape"] == "open"


def test_only_the_first_record_is_published(
    capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    # The guard and the normal end of a run both publish. The first one wins.
    monkeypatch.setattr(run, "PUBLISHED", threading.Event())
    assert run.publish({"status": "memory_limit"}) is True
    assert run.publish({"status": "ok"}) is False
    assert capsys.readouterr().out.splitlines() == ['{"status": "memory_limit"}']


def test_a_run_over_the_limit_between_checks_is_not_ok() -> None:
    record = {"status": "ok", "peak_rss_mb": 2048.0}
    over, code = run.final_record(record, limit_mb=1024.0)
    assert over["status"] == "memory_limit"
    assert over["limit_mb"] == 1024.0
    assert code == run.MEMORY_LIMIT_EXIT
    assert run.final_record(record, limit_mb=4096.0) == (record, None)
    assert run.final_record(record, limit_mb=None) == (record, None)


@pytest.fixture(scope="module")
def table(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Two stores with their obs fragments, the catalog, and the obs table."""
    import build_catalog
    import build_obs_table
    import census2_dataset_writer

    root = tmp_path_factory.mktemp("bench-table")
    prefix = root / "out"
    genes = [f"ENSG{j:011d}" for j in range(8)]
    for k, n_obs in enumerate((30, 20)):
        x = sp.random(n_obs, len(genes), density=0.3, format="csr", random_state=k)
        obs = pd.DataFrame(
            {
                "dataset_id": pd.Categorical([f"ds-{k}"] * n_obs),
                "cell_type": pd.Categorical(["T cell", "B cell"] * (n_obs // 2)),
                "donor_id": pd.Categorical([f"donor-{k}"] * n_obs),
                "observation_joinid": [f"j{i}" for i in range(n_obs)],
            },
            index=[str(i) for i in range(n_obs)],
        )
        var = pd.DataFrame({"feature_id": genes}, index=genes)
        adata = anndata.AnnData(X=sp.csr_matrix(x, dtype=np.float32), obs=obs, var=var)
        adata.uns = {"organism_ontology_term_id": "NCBITaxon:9606"}
        path = prefix / f"v-{k}.zarr"
        census2_dataset_writer.write(
            adata,
            path,
            summary={
                "version": census2_dataset_writer.SUMMARY_VERSION,
                "dataset_id": f"ds-{k}",
                "dataset_version_id": f"v-{k}",
                "gene_reference_version": "GENCODE v48",
            },
            fragment=census2_dataset_writer.fragment_location(str(path)),
        )
    build_catalog.build(str(prefix), root / "catalog.parquet")
    build_obs_table.build(str(root / "catalog.parquet"), str(root / "obs"))
    return root


@pytest.mark.parametrize(
    ("shape", "expected"),
    [("query", (2, 25)), ("count", (0, 25)), ("corpus_obs", (0, 50))],
)
def test_the_table_shapes_run(
    table: Path, shape: str, expected: tuple[int, int]
) -> None:
    reset_counter()
    record = run.run(
        str(table / "catalog.parquet"),
        shape,
        2,
        None,
        query={"cell_type": "T cell"},
        obs_table=str(table / "obs"),
    )
    assert record["status"] == "ok" and record["via"] == "table"
    assert (record["stores_opened"], record["cells"]) == expected
    assert record["files_planned"] == 2
    assert record["files_total"] == 2
    assert record["manifests"] == 1


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("cell_type = 'T cell'", [True, False, False, True]),
        ("cell_type != 'T cell'", [False, True, False, False]),
        ("cell_type IN ('B cell')", [False, True, False, False]),
        ("cell_type NOT IN ('B cell')", [True, False, False, True]),
        ("cell_type IS NULL", [False, False, True, False]),
        ("nnz >= 500 AND cell_type != 'B cell'", [True, False, False, False]),
        ("is_primary_data = true OR nnz < 100", [True, True, False, True]),
    ],
)
def test_evaluate_follows_the_scan_rules(text: str, expected: list[bool]) -> None:
    """A null label matches no comparison, as in the scan of the table."""
    from pyiceberg.expressions import parser

    obs = pd.DataFrame(
        {
            "cell_type": pd.Categorical(["T cell", "B cell", None, "T cell"]),
            "nnz": np.array([600, 700, 800, 50], dtype=np.int64),
            "is_primary_data": np.array([True, True, False, False]),
        }
    )
    assert list(run.evaluate(parser.parse(text), obs)) == expected


def test_pruning_reports_the_files_and_row_groups_that_can_match(table: Path) -> None:
    from bench import pruning

    from census2_sdk import catalog, obs_table

    frame = catalog.load(catalog.read(table / "catalog.parquet"))
    opened = obs_table.open_table(str(table / "obs"))
    out = pruning.prune(frame, opened, "cell_type = 'B cell'")
    assert (out["files"], out["files_after_bounds"], out["rows"]) == (2, 2, 25)
    assert out["row_groups"] == out["row_groups_after_stats"] == 2
    # A label outside the bounds of both files: no file, no row group.
    out = pruning.prune(frame, opened, "cell_type = 'Z cell'")
    assert (out["files_after_bounds"], out["row_groups"], out["rows"]) == (0, 0, 0)


def test_arrow_counting_counts_the_table_reads(
    table: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Every file that pyiceberg opens is one HEAD, and each read of it one GET."""
    from pyiceberg.io.pyarrow import PyArrowFileIO

    monkeypatch.setattr(PyArrowFileIO, "_initialize_fs", PyArrowFileIO._initialize_fs)
    counting.install_arrow_counting()
    reset_counter()
    record = run.run(
        str(table / "catalog.parquet"),
        "count",
        2,
        None,
        query="cell_type = 'B cell'",
        obs_table=str(table / "obs"),
    )
    assert record["cells"] == 25
    requests, heads = record["requests_by_kind"], record["heads_by_kind"]
    # Two fragments: one HEAD and at least one GET each.
    assert heads["fragment"] == 2 and requests["fragment"] >= 4
    # The metadata file, the manifest list, and one manifest, with the same rule.
    assert heads["table_metadata"] >= 3 and requests["table_metadata"] >= 6
    assert record["mb_by_kind"]["fragment"] > 0


@pytest.mark.parametrize(
    ("query", "by", "narrow"),
    [
        ({"cell_type": "T cell"}, ["donor_id"], False),
        ("cell_type != 'B cell' AND donor_id IN ('donor-1')", ["cell_type"], True),
    ],
)
def test_a_count_from_the_stores_equals_the_count_from_the_table(
    table: Path, query: Any, by: list[str], narrow: bool
) -> None:
    catalog = str(table / "catalog.parquet")
    with_table = run.run(
        catalog, "count", 2, None, query=query, obs_table=str(table / "obs"), by=by
    )
    if isinstance(query, str):
        # The store path takes an expression for a count, with its own catalog test.
        import census2_sdk

        census = census2_sdk.open_census([catalog])
        groups, _, cells = run.count_from_stores(census, query, by, 2, narrow=narrow)
    else:
        record = run.run(catalog, "count", 2, None, query=query, by=by, narrow=narrow)
        groups, cells = record["groups"], record["cells"]
    assert (groups, cells) == (with_table["groups"], with_table["cells"])
