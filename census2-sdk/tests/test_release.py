"""Tests for a release manifest as a source. No network, no AWS.

The fixture writes three small stores and their live catalog, as `test_sdk` does. Each
test then makes a release directory by hand: a copy of the catalog and a manifest.
"""

import hashlib
import io
import json
import shutil
from pathlib import Path
from typing import Any

import build_catalog
import pandas as pd
import pytest
from test_sdk import (
    REFERENCE_IDS,
    STORES,
    make_store,
    one_store_catalog,
    query,
    same,
    write_one,
)

import census2_sdk
from census2_sdk import catalog, release

NAME = "2026-09-01"


@pytest.fixture
def live(tmp_path: Path) -> dict[str, Path]:
    """Three stores under a prefix, and the live catalog of that prefix."""
    prefix = tmp_path / "out"
    prefix.mkdir()
    for seed, (dataset_id, spec) in enumerate(STORES.items()):
        write_one(prefix, make_store(dataset_id, *spec, seed=seed), dataset_id)
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
