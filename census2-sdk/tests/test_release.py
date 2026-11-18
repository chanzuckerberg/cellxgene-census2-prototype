"""Tests for a release: the manifest as a source, and the release builder.

No network, no AWS. The fixture writes three small stores and their live catalog, as
`test_sdk` does, and the per-gene nnz file next to each store, as the ingest task
does. The tests of the manifest make a release directory by hand: a copy of
the catalog and a manifest. Thus they do not depend on the release builder. The tests
of the release builder run `tools/build_release.py` on the fixture.
"""

import hashlib
import io
import json
import shutil
from pathlib import Path
from typing import Any, cast

import anndata
import build_catalog
import build_release
import census2_dataset_writer
import numpy as np
import pandas as pd
import pytest
from test_sdk import (
    REFERENCE_IDS,
    STORES,
    X,
    make_store,
    one_store_catalog,
    query,
    same,
    set_summary,
    write_one,
)

import census2_sdk
from census2_sdk import catalog, release

NAME = "2026-09-01"


def nnz_file(prefix: Path, dataset_id: str) -> Path:
    """The location of the per-gene nnz file of one fixture store."""
    return prefix / census2_dataset_writer.gene_nnz_name(f"v-{dataset_id}")


def write_store(prefix: Path, adata: anndata.AnnData, dataset_id: str) -> None:
    """One store, and its per-gene nnz file next to it."""
    write_one(prefix, adata, dataset_id)
    nnz = census2_dataset_writer.gene_nnz(prefix / f"v-{dataset_id}.zarr")
    np.save(nnz_file(prefix, dataset_id), nnz)


@pytest.fixture
def live(tmp_path: Path) -> dict[str, Path]:
    """Three stores under a prefix, and the live catalog of that prefix."""
    prefix = tmp_path / "out"
    prefix.mkdir()
    for seed, (dataset_id, spec) in enumerate(STORES.items()):
        write_store(prefix, make_store(dataset_id, *spec, seed=seed), dataset_id)
    build_catalog.build(str(prefix), tmp_path / "current.parquet")
    return {"prefix": prefix, "catalog": tmp_path / "current.parquet"}


def good_manifest(snapshot: Path) -> dict[str, Any]:
    """The manifest of a release whose catalog snapshot is the file `snapshot`."""
    rows = catalog.load(snapshot)
    return {
        "version": release.MANIFEST_VERSION,
        "release": NAME,
        "gene_reference_version": {"homo_sapiens": "GENCODE v48"},
        "stores": [{"dataset_version_id": v} for v in rows["dataset_version_id"]],
        "files": {
            release.CATALOG_NAME: hashlib.sha256(snapshot.read_bytes()).hexdigest()
        },
    }


def make_release(live: dict[str, Path], **changes: Any) -> Path:
    """A release directory of the live catalog. Return the location of its manifest.

    `changes` replaces fields of the manifest, so that a test can write a malformed one.
    """
    directory = live["catalog"].parent / "releases" / NAME
    directory.mkdir(parents=True)
    snapshot = directory / release.CATALOG_NAME
    shutil.copy(live["catalog"], snapshot)
    manifest = directory / release.MANIFEST_NAME
    manifest.write_text(json.dumps(good_manifest(snapshot) | changes))
    return manifest


def test_a_manifest_gives_the_results_of_the_live_catalog(
    live: dict[str, Path],
) -> None:
    manifest = make_release(live)
    pinned = census2_sdk.open_census([manifest])
    current = census2_sdk.open_census([live["catalog"]])

    result = query(pinned, {"tissue": "lung"})
    same(result, query(current, {"tissue": "lung"}))
    for obs_filter in (None, {"tissue": "lung"}, {"disease": "flu"}):
        want = current.count_cells("homo_sapiens", obs_filter)
        assert pinned.count_cells("homo_sapiens", obs_filter) == want
    assert pinned.count_cells("homo_sapiens") == 75

    digest = hashlib.sha256(manifest.read_bytes()).hexdigest()
    assert result.uns["census2"]["releases"] == {NAME: digest}
    assert result.uns["census2"]["sources"] == [str(manifest)]
    assert result.uns["census2"]["stores_in_catalog"] == 3
    # A catalog is not a release.
    assert query(current, {"tissue": "lung"}).uns["census2"]["releases"] == {}


def test_a_manifest_and_a_catalog_compose(
    live: dict[str, Path], tmp_path: Path
) -> None:
    manifest = make_release(live)
    lab = tmp_path / "lab"
    lab.mkdir()
    adata = make_store("ds-lab", 8, REFERENCE_IDS[:6], ["lung"], ["normal"], 21)
    census = census2_sdk.open_census(
        [manifest, one_store_catalog(lab, adata, "ds-lab")]
    )

    assert census.count_cells("homo_sapiens") == 75 + 8
    result = query(census, {})
    assert result.n_obs == 75 + 8
    assert set(result.obs["dataset_id"]) == {*STORES, "ds-lab"}
    assert result.uns["census2"]["stores_opened"] == 4
    assert list(result.uns["census2"]["releases"]) == [NAME]


def test_a_dataset_in_a_release_and_in_a_catalog_is_refused(
    live: dict[str, Path],
) -> None:
    manifest = make_release(live)
    with pytest.raises(census2_sdk.SourceConflict, match="more than one source"):
        census2_sdk.open_census([manifest, live["catalog"]])


def test_a_release_in_two_sources_is_refused(live: dict[str, Path]) -> None:
    manifest = make_release(live)
    with pytest.raises(census2_sdk.SourceConflict, match=f"release '{NAME}'"):
        census2_sdk.open_census([manifest, manifest])


def test_validate_manifest_returns_a_good_manifest(live: dict[str, Path]) -> None:
    manifest = good_manifest(live["catalog"])
    assert release.validate_manifest(manifest) == manifest


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        ({"version": 2}, "version 2"),
        ({"version": "1"}, "version '1'"),
        ({"version": True}, "version True"),
        ({"release": ""}, "release"),
        ({"gene_reference_version": {"homo_sapiens": 48}}, "gene_reference_version"),
        ({"stores": {"dataset_version_id": "v-ds-a"}}, "stores"),
        ({"stores": [{"dataset_version_id": "v-ds-a", "n": 1}]}, "stores"),
        ({"stores": ["v-ds-a"]}, "stores"),
        ({"files": {release.CATALOG_NAME: "abc"}}, "SHA-256"),
        ({"files": {}}, f"no {release.CATALOG_NAME}"),
        ({"writer_version": 1}, r"unknown \['writer_version'\]"),
    ],
)
def test_a_malformed_manifest_is_refused(
    live: dict[str, Path], changes: dict[str, Any], message: str
) -> None:
    manifest = make_release(live, **changes)
    with pytest.raises(release.ReleaseError, match=message) as error:
        census2_sdk.open_census([manifest])
    assert str(manifest) in str(error.value)


@pytest.mark.parametrize("content", [b"not a catalog", b"[1, 2]", b'{"version": 1}'])
def test_a_source_that_is_no_catalog_and_no_manifest_is_refused(
    tmp_path: Path, content: bytes
) -> None:
    source = tmp_path / "source"
    source.write_bytes(content)
    with pytest.raises(release.ReleaseError, match=str(source)):
        census2_sdk.open_census([source])


def test_a_release_without_its_snapshot_is_refused(live: dict[str, Path]) -> None:
    manifest = make_release(live)
    (manifest.parent / release.CATALOG_NAME).unlink()
    with pytest.raises(release.ReleaseError, match="holds no catalog.parquet") as error:
        census2_sdk.open_census([manifest])
    assert str(manifest) in str(error.value)


def test_a_changed_snapshot_is_refused(live: dict[str, Path], tmp_path: Path) -> None:
    manifest = make_release(live)
    # The catalog builder writes a catalog of one other store over the snapshot.
    lab = tmp_path / "lab"
    lab.mkdir()
    adata = make_store("ds-lab", 8, REFERENCE_IDS[:6], ["lung"], ["normal"], 21)
    shutil.copy(
        one_store_catalog(lab, adata, "ds-lab"), manifest.parent / release.CATALOG_NAME
    )
    with pytest.raises(release.ReleaseError, match="SHA-256") as error:
        census2_sdk.open_census([manifest])
    assert str(manifest) in str(error.value)


def parquet(frame: pd.DataFrame) -> bytes:
    content = io.BytesIO()
    frame.to_parquet(content)
    return content.getvalue()


@pytest.mark.parametrize(
    "content", [b"PAR1 and no table", parquet(pd.DataFrame({"n_obs": [1]}))]
)
def test_a_pinned_snapshot_that_is_no_catalog_is_refused(
    live: dict[str, Path], content: bytes
) -> None:
    """The manifest pins the content, so the digest is right and the content is not."""
    manifest = make_release(live)
    (manifest.parent / release.CATALOG_NAME).write_bytes(content)
    fields = json.loads(manifest.read_text())
    fields["files"][release.CATALOG_NAME] = hashlib.sha256(content).hexdigest()
    manifest.write_text(json.dumps(fields))
    with pytest.raises(release.ReleaseError, match="is not a catalog") as error:
        census2_sdk.open_census([manifest])
    assert str(manifest) in str(error.value)


def test_a_snapshot_with_other_stores_is_refused(live: dict[str, Path]) -> None:
    manifest = make_release(live, stores=[{"dataset_version_id": "v-ds-a"}])
    with pytest.raises(release.ReleaseError, match="not the stores") as error:
        census2_sdk.open_census([manifest])
    assert str(manifest) in str(error.value)


def test_an_absent_source_is_not_found(tmp_path: Path) -> None:
    """An absent manifest fails as an absent catalog does."""
    for name in (release.MANIFEST_NAME, "catalog.parquet"):
        with pytest.raises(FileNotFoundError):
            census2_sdk.open_census([tmp_path / name])


def test_a_missing_store_fails_the_query(live: dict[str, Path]) -> None:
    manifest = make_release(live)
    store = live["prefix"] / "v-ds-b.zarr"
    shutil.rmtree(store)
    # `open_census` opens no store.
    census = census2_sdk.open_census([manifest])
    assert census.count_cells("homo_sapiens") == 75
    with pytest.raises(FileNotFoundError) as error:
        census.get_anndata("homo_sapiens")
    for part in (str(store), "ds-b", str(manifest)):
        assert part in str(error.value)
    # A query that excludes the dataset does not open its store.
    assert query(census, {"dataset_id": "ds-a"}).n_obs == 30


# --- the release builder -------------------------------------------------------


def add_store(live: dict[str, Path], adata: anndata.AnnData, dataset_id: str) -> None:
    """Write one more store, and build the live catalog again."""
    write_store(live["prefix"], adata, dataset_id)
    build_catalog.build(str(live["prefix"]), live["catalog"])


def built(live: dict[str, Path]) -> Path:
    """Build the release `NAME` of the live catalog. Return its directory."""
    releases = live["catalog"].parent / "releases"
    build_release.build(live["catalog"], releases, NAME)
    return releases / NAME


def test_the_snapshot_is_the_live_catalog_content(live: dict[str, Path]) -> None:
    directory = built(live)
    snapshot = (directory / release.CATALOG_NAME).read_bytes()
    assert snapshot == live["catalog"].read_bytes()
    files = {str(p.relative_to(directory)) for p in directory.rglob("*") if p.is_file()}
    assert files == {
        release.CATALOG_NAME,
        "presence/homo_sapiens.parquet",
        "aggregates/var.parquet",
        release.MANIFEST_NAME,
    }


def test_the_manifest_pins_the_stores_and_the_files(live: dict[str, Path]) -> None:
    directory = built(live)
    manifest = json.loads((directory / release.MANIFEST_NAME).read_text())
    assert release.validate_manifest(manifest) == manifest
    assert manifest["version"] == 1
    assert manifest["release"] == NAME
    assert manifest["gene_reference_version"] == {"homo_sapiens": "GENCODE v48"}
    assert manifest["stores"] == [
        {"dataset_version_id": f"v-{dataset_id}"} for dataset_id in STORES
    ]
    # The digests are those of the other files of the release, and of no other file.
    others = {
        str(p.relative_to(directory)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in directory.rglob("*")
        if p.is_file() and p.name != release.MANIFEST_NAME
    }
    assert manifest["files"] == others


def test_the_presence_rows_are_the_measured_columns(live: dict[str, Path]) -> None:
    directory = built(live)
    table = pd.read_parquet(directory / "presence" / "homo_sapiens.parquet")
    assert list(table.columns) == ["dataset_id", "dataset_version_id", "measured"]
    assert list(table["dataset_id"]) == list(STORES)
    assert list(table["dataset_version_id"]) == [f"v-{d}" for d in STORES]
    for dataset_id, positions in zip(
        table["dataset_id"], table["measured"], strict=True
    ):
        assert [REFERENCE_IDS[j] for j in positions] == STORES[dataset_id][1]


def test_the_aggregates_add_the_fixture_matrices(live: dict[str, Path]) -> None:
    directory = built(live)
    table = pd.read_parquet(directory / "aggregates" / "var.parquet")
    assert list(table.columns) == ["organism", "feature_id", "nnz", "n_measured_obs"]
    assert set(table["organism"]) == {"homo_sapiens"}
    assert list(table["feature_id"]) == REFERENCE_IDS
    # The expected values come from the matrices that the fixture wrote.
    nnz = np.zeros(len(REFERENCE_IDS), dtype=np.int64)
    n_measured_obs = np.zeros(len(REFERENCE_IDS), dtype=np.int64)
    for seed, (dataset_id, spec) in enumerate(STORES.items()):
        counts = X(make_store(dataset_id, *spec, seed=seed)).toarray()
        nnz += np.count_nonzero(counts, axis=0)
        n_measured_obs += np.isin(REFERENCE_IDS, spec[1]) * spec[0]
    assert table["nnz"].dtype == table["n_measured_obs"].dtype == np.int64
    assert list(table["nnz"]) == list(nnz)
    assert list(table["n_measured_obs"]) == list(n_measured_obs)
    # Each store measured a part of the genes, and no gene has more cells than that.
    assert 0 < nnz.sum() < n_measured_obs.sum()
    assert set(n_measured_obs) > {75}


def test_a_store_without_its_nnz_file_stops_the_build(live: dict[str, Path]) -> None:
    missing = nnz_file(live["prefix"], "ds-b")
    missing.unlink()
    with pytest.raises(FileNotFoundError, match=missing.name) as error:
        built(live)
    assert "v-ds-b.zarr" in str(error.value)
    directory = live["catalog"].parent / "releases" / NAME
    assert not (directory / release.MANIFEST_NAME).exists()


def wrong_nnz(change: str) -> Any:
    """The nnz array of the store `ds-a`, which measured genes 0 to 11, with one fault."""
    nnz = np.zeros(len(REFERENCE_IDS), dtype=np.int64)
    nnz[:12] = 3
    if change == "length":
        return nnz[:-1]
    if change == "dtype":
        return nnz.astype(np.int32)
    if change == "unmeasured":
        nnz[15] = 1
    if change == "above n_obs":
        nnz[0] = 31
    if change == "negative":
        nnz[0] = -1
    return nnz


@pytest.mark.parametrize(
    ("change", "message"),
    [
        ("length", r"shape \(19,\)"),
        ("dtype", "int32"),
        ("unmeasured", "did not measure"),
        ("above n_obs", "outside 0 to 30"),
        ("negative", "outside 0 to 30"),
    ],
)
def test_a_wrong_nnz_file_stops_the_build(
    live: dict[str, Path], change: str, message: str
) -> None:
    path = nnz_file(live["prefix"], "ds-a")
    np.save(path, wrong_nnz(change))
    with pytest.raises(ValueError, match=message) as error:
        built(live)
    assert path.name in str(error.value)
    directory = live["catalog"].parent / "releases" / NAME
    assert not (directory / release.MANIFEST_NAME).exists()


@pytest.mark.parametrize("content", [b"", b"not an array", b"PK\x03\x04 no archive"])
def test_an_nnz_file_that_is_no_array_stops_the_build(
    live: dict[str, Path], content: bytes
) -> None:
    path = nnz_file(live["prefix"], "ds-a")
    path.write_bytes(content)
    with pytest.raises(ValueError, match=f"{path.name} is not a .npy array"):
        built(live)


def test_the_builder_refuses_cell_counts_that_exceed_int64(
    live: dict[str, Path],
) -> None:
    """Each count is in range, and the sum of two of them is not."""
    frame = pd.read_parquet(live["catalog"])
    frame["n_obs"] = 2**62
    frame.to_parquet(live["catalog"])
    with pytest.raises(ValueError, match="homo_sapiens.*exceed int64"):
        built(live)


def test_each_organism_has_its_presence_file(live: dict[str, Path]) -> None:
    """A mouse store has another `var`, and its positions are in that `var`."""
    mouse_ids = [f"ENSMUSG{j:011d}" for j in range(7)]
    adata = make_store(
        "ds-m", 5, mouse_ids[2:5], ["lung"], ["normal"], 31, var_ids=mouse_ids
    )
    adata.uns["organism_ontology_term_id"] = "NCBITaxon:10090"
    write_store(live["prefix"], adata, "ds-m")
    record = dict(
        census2_dataset_writer.read_summary(live["prefix"] / "v-ds-m.zarr") or {}
    )
    set_summary(
        live["prefix"] / "v-ds-m.zarr",
        record | {"gene_reference_version": "GENCODE M37"},
    )
    build_catalog.build(str(live["prefix"]), live["catalog"])

    directory = built(live)

    manifest = json.loads((directory / release.MANIFEST_NAME).read_text())
    assert manifest["gene_reference_version"] == {
        "homo_sapiens": "GENCODE v48",
        "mus_musculus": "GENCODE M37",
    }
    mouse = pd.read_parquet(directory / "presence" / "mus_musculus.parquet")
    assert list(mouse["dataset_id"]) == ["ds-m"]
    assert list(mouse["measured"][0]) == [2, 3, 4]
    human = pd.read_parquet(directory / "presence" / "homo_sapiens.parquet")
    assert list(human["dataset_id"]) == list(STORES)
    # The aggregates hold the genes of each organism, from the stores of that organism.
    table = pd.read_parquet(directory / "aggregates" / "var.parquet")
    assert list(table["organism"]) == ["homo_sapiens"] * 20 + ["mus_musculus"] * 7
    genes = table[table["organism"] == "mus_musculus"]
    assert list(genes["feature_id"]) == mouse_ids
    assert list(genes["nnz"]) == list(np.count_nonzero(X(adata).toarray(), axis=0))
    assert list(genes["n_measured_obs"]) == [0, 0, 5, 5, 5, 0, 0]
    assert list(table["feature_id"][:20]) == REFERENCE_IDS
    census = census2_sdk.open_census([directory / release.MANIFEST_NAME])
    assert census.count_cells("mus_musculus") == 5
    assert census.count_cells("homo_sapiens") == 75


def test_a_later_write_does_not_change_the_release(live: dict[str, Path]) -> None:
    manifest = built(live) / release.MANIFEST_NAME
    before = query(census2_sdk.open_census([manifest]), {"tissue": "lung"})
    # The live catalog matches the release until the next write.
    same(before, query(census2_sdk.open_census([live["catalog"]]), {"tissue": "lung"}))

    adata = make_store("ds-new", 9, REFERENCE_IDS[:4], ["lung"], ["normal"], 41)
    add_store(live, adata, "ds-new")

    pinned = census2_sdk.open_census([manifest])
    after = query(pinned, {"tissue": "lung"})
    same(before, after)
    assert after.n_obs == before.n_obs
    assert after.uns["census2"]["stores_in_catalog"] == 3
    assert pinned.count_cells("homo_sapiens") == 75
    current = census2_sdk.open_census([live["catalog"]])
    assert current.count_cells("homo_sapiens") == 75 + 9


@pytest.mark.parametrize(
    ("owner", "step", "written"),
    [
        # The build stops after the snapshot, and then at the rename of the manifest.
        (build_release, "measured_genes", [release.CATALOG_NAME]),
        (
            Path,
            "replace",
            [
                release.CATALOG_NAME,
                "presence/homo_sapiens.parquet",
                "aggregates/var.parquet",
            ],
        ),
    ],
)
def test_a_stopped_build_leaves_no_manifest(
    live: dict[str, Path],
    monkeypatch: pytest.MonkeyPatch,
    owner: Any,
    step: str,
    written: list[str],
) -> None:
    def stop(*args: Any) -> Any:
        raise RuntimeError("stop")

    with monkeypatch.context() as patch:
        patch.setattr(owner, step, stop)
        with pytest.raises(RuntimeError, match="stop"):
            built(live)
    directory = live["catalog"].parent / "releases" / NAME
    for path in written:
        assert (directory / path).exists()
    assert not (directory / release.MANIFEST_NAME).exists()
    with pytest.raises(FileNotFoundError):
        census2_sdk.open_census([directory / release.MANIFEST_NAME])
    # A release is a new directory, so the next build needs the old one removed.
    with pytest.raises(FileExistsError, match=NAME):
        built(live)
    shutil.rmtree(directory)
    assert built(live) == directory
    assert (
        census2_sdk.open_census([directory / release.MANIFEST_NAME]).count_cells(
            "homo_sapiens"
        )
        == 75
    )


@pytest.mark.parametrize("rows", [[0], [0, 1, 2]])
def test_the_builder_refuses_a_row_without_an_organism(
    live: dict[str, Path], rows: list[int]
) -> None:
    frame = pd.read_parquet(live["catalog"])
    frame.loc[rows, "organism"] = None
    frame.to_parquet(live["catalog"])
    with pytest.raises(ValueError, match=f"{len(rows)} rows record no organism"):
        built(live)
    assert not (live["catalog"].parent / "releases").exists()


def test_the_builder_refuses_a_directory_that_holds_a_release(
    live: dict[str, Path],
) -> None:
    directory = built(live)
    manifest = (directory / release.MANIFEST_NAME).read_bytes()
    adata = make_store("ds-new", 9, REFERENCE_IDS[:4], ["lung"], ["normal"], 41)
    add_store(live, adata, "ds-new")
    with pytest.raises(FileExistsError, match=NAME):
        built(live)
    assert (directory / release.MANIFEST_NAME).read_bytes() == manifest
    assert (directory / release.CATALOG_NAME).read_bytes() != live[
        "catalog"
    ].read_bytes()


@pytest.mark.parametrize("name", ["", ".", "..", "a/b"])
def test_the_builder_refuses_a_name_that_is_no_directory_name(
    live: dict[str, Path], name: str
) -> None:
    releases = live["catalog"].parent / "releases"
    with pytest.raises(ValueError, match="release name"):
        build_release.build(live["catalog"], releases, name)
    assert not releases.exists()


def test_the_builder_refuses_two_reference_versions_of_one_organism(
    live: dict[str, Path],
) -> None:
    store = live["prefix"] / "v-ds-b.zarr"
    record = dict(census2_dataset_writer.read_summary(store) or {})
    set_summary(store, record | {"gene_reference_version": "GENCODE v47"})
    build_catalog.build(str(live["prefix"]), live["catalog"])
    with pytest.raises(ValueError, match=r"homo_sapiens.*GENCODE v47"):
        built(live)
    assert not (live["catalog"].parent / "releases").exists()


def test_the_builder_refuses_stores_with_another_var_length(
    live: dict[str, Path],
) -> None:
    adata = make_store(
        "ds-short",
        4,
        REFERENCE_IDS[:3],
        ["lung"],
        ["normal"],
        51,
        var_ids=REFERENCE_IDS[:7],
    )
    add_store(live, adata, "ds-short")
    with pytest.raises(ValueError, match=r"homo_sapiens.*\[7, 20\]"):
        built(live)
    directory = live["catalog"].parent / "releases" / NAME
    assert not (directory / release.MANIFEST_NAME).exists()


def test_the_builder_refuses_a_store_without_measured(live: dict[str, Path]) -> None:
    adata = make_store("ds-old", 4, REFERENCE_IDS[:3], ["lung"], ["normal"], 61)
    adata.var = cast("pd.DataFrame", adata.var).drop(columns="measured")
    add_store(live, adata, "ds-old")
    with pytest.raises(ValueError, match=r"v-ds-old.zarr.*measured"):
        built(live)


def test_the_command_builds_a_release(
    live: dict[str, Path], capsys: pytest.CaptureFixture[str]
) -> None:
    releases = live["catalog"].parent / "releases"
    assert build_release.main([str(live["catalog"]), str(releases), NAME]) == 0
    assert "3 stores" in capsys.readouterr().out
    assert (releases / NAME / release.MANIFEST_NAME).exists()
    assert build_release.main([str(live["catalog"])]) == 2
