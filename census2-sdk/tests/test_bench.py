"""Tests for the store-count benchmark tools (task 6.4). No network, no AWS."""

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
from bench import catalogs, matrix, run
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


def test_get_ranges_is_counted() -> None:
    store = run.CountingStore(zarr.storage.MemoryStore())
    buffer = default_buffer_prototype().buffer.from_bytes(bytes(range(100)))

    async def ranges() -> list[Any]:
        await store.set("k", buffer)
        requests = [RangeByteRequest(0, 10), RangeByteRequest(50, 70)]
        prototype = default_buffer_prototype()
        return [g async for g in store.get_ranges("k", requests, prototype=prototype)]

    run.COUNTER.requests = run.COUNTER.nbytes = 0
    groups = asyncio.run(ranges())
    assert run.COUNTER.requests == len(groups) >= 1
    assert run.COUNTER.nbytes == 30


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
    run.COUNTER.requests = run.COUNTER.nbytes = 0
    obs_filter = {"cell_type": "T cell"} if shape == "filter" else None
    record = run.run(
        str(stores / "catalog.parquet"), shape, workers, None, obs_filter=obs_filter
    )
    assert record["status"] == "ok"
    assert record["catalog_stores"] == 2
    expected = {"single": (1, 30), "open": (2, 0)}.get(shape, (2, 50))
    assert (record["stores_opened"], record["cells"]) == expected
    assert record["requests"] > 0
    assert record["mb_payload"] > 0


def test_serving_uris_open_a_compacted_store_once() -> None:
    rows = pd.DataFrame(
        {
            "store_uri": ["a.zarr", "b.zarr", "c.zarr"],
            "compacted_store_uri": ["ab.zarr", "ab.zarr", None],
        }
    )
    assert run.serving_uris(rows) == ["ab.zarr", "c.zarr"]
    assert run.serving_uris(rows, per_dataset=True) == ["a.zarr", "b.zarr", "c.zarr"]
    assert run.serving_uris(rows.drop(columns="compacted_store_uri")) == [
        "a.zarr",
        "b.zarr",
        "c.zarr",
    ]


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
