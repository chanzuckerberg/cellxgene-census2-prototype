"""The release manifest: the file that pins the catalog and the stores of one release.

A release is a directory that the release builder writes (design 3.3). The directory
holds a snapshot of the live catalog in `CATALOG_NAME`, and the manifest in
`MANIFEST_NAME`. The release builder writes the manifest last, so a directory with a
manifest is a complete release.

`open_census` accepts a manifest as a source (design 4.1). `open_manifest` reads the
catalog snapshot from the directory of the manifest and checks it against the manifest.
The SDK then uses the snapshot as it uses a catalog.

A release holds one gene reference version for each organism (design 1.5). Thus each
row of the snapshot must record the version that the manifest gives for its organism.
The SDK refuses a release with a row that records no version, and the release builder
refuses a catalog with such a row.

The manifest is UTF-8 JSON, and `Manifest` declares its fields. `validate_manifest`
checks them. The design record reserves the format of the manifest (design section 8),
so these fields are the format of `MANIFEST_VERSION`.

The manifest pins the identity of each store and not its content. It holds no content
digest for a store and no `writer_version`, because no store and no catalog row records
either value.
"""

import hashlib
import io
import json
import pathlib
import re
from typing import Any, TypedDict, cast

from census2_sdk import catalog

# The format of the manifest. Raise it when a field changes meaning or goes away.
# The SDK refuses a version that it does not know.
MANIFEST_VERSION = 1
MANIFEST_NAME = "manifest.json"
CATALOG_NAME = "catalog.parquet"

_SHA256 = re.compile(r"[0-9a-f]{64}")


class StoreEntry(TypedDict):
    """One store of a release."""

    dataset_version_id: str


class Manifest(TypedDict):
    """The fields of a release manifest."""

    version: int
    # The name of the release, which is the name of its directory.
    release: str
    # The gene reference version of each organism, such as "GENCODE v48" (design 1.5).
    gene_reference_version: dict[str, str]
    # One entry for each row of the catalog snapshot.
    stores: list[StoreEntry]
    # The SHA-256 of each other file of the release, in hex. The key is the path of
    # the file in the release directory.
    files: dict[str, str]


class ReleaseError(ValueError):
    """A release that the SDK cannot use.

    The manifest is malformed, or the catalog snapshot is absent or differs from the
    manifest.
    """


def _name(value: object) -> bool:
    return isinstance(value, str) and bool(value)


def validate_manifest(obj: object) -> Manifest:
    """Return `obj` as a manifest, or raise `ReleaseError` with the cause.

    A manifest comes from a file, which is external input, so this function checks
    its type at run time. The manifest must have exactly the fields of `Manifest`,
    and `version` must be `MANIFEST_VERSION`. `files` must hold `CATALOG_NAME`.
    """
    if not isinstance(obj, dict):
        raise ReleaseError(
            f"a release manifest is a JSON object, not {type(obj).__name__}"
        )
    fields = set(Manifest.__annotations__)
    if set(obj) != fields:
        missing, extra = sorted(fields - set(obj)), sorted(set(obj) - fields)
        raise ReleaseError(f"manifest fields: missing {missing}, unknown {extra}")
    version = obj["version"]
    if type(version) is not int or version != MANIFEST_VERSION:
        raise ReleaseError(
            f"manifest version {version!r}; this SDK reads {MANIFEST_VERSION}"
        )
    if not _name(obj["release"]):
        raise ReleaseError(f"manifest release is {obj['release']!r}, not a name")
    versions = obj["gene_reference_version"]
    if not isinstance(versions, dict) or not all(
        _name(organism) and _name(reference) for organism, reference in versions.items()
    ):
        raise ReleaseError(
            "manifest gene_reference_version must map each organism to a version"
        )
    stores = obj["stores"]
    if not isinstance(stores, list) or not all(
        isinstance(store, dict)
        and set(store) == {"dataset_version_id"}
        and _name(store["dataset_version_id"])
        for store in stores
    ):
        raise ReleaseError(
            "manifest stores must be a list of objects with one dataset_version_id each"
        )
    files = obj["files"]
    if not isinstance(files, dict) or not all(
        _name(path) and isinstance(digest, str) and _SHA256.fullmatch(digest)
        for path, digest in files.items()
    ):
        raise ReleaseError("manifest files must map each path to a SHA-256 in hex")
    if CATALOG_NAME not in files:
        raise ReleaseError(f"manifest files holds no {CATALOG_NAME}")
    return cast("Manifest", obj)


def _sibling(source: str, name: str) -> str:
    """The location of the file `name` in the directory of `source`."""
    if catalog.is_remote(source):
        return f"{source.rpartition('/')[0]}/{name}"
    return str(pathlib.Path(source).with_name(name))


def _strings(values: list[Any]) -> list[str | None]:
    return [None if value is None else str(value) for value in values]


def open_manifest(source: str, raw: bytes) -> tuple[Manifest, io.BytesIO]:
    """The manifest in `raw`, and the content of the catalog snapshot that it pins.

    `raw` is the content of `source`. The snapshot is `CATALOG_NAME` in the directory
    of `source`. Its SHA-256 must be the digest in the manifest. It must be a catalog,
    and its rows must name the stores of the manifest. The rows of each organism must
    record the one gene reference version that the manifest gives for that organism.
    The manifest must name no other organism. Each failure raises `ReleaseError`, and
    the message names `source`.
    """
    try:
        manifest = validate_manifest(json.loads(raw))
    except ValueError as e:
        raise ReleaseError(
            f"{source} is not a catalog and not a release manifest: {e}"
        ) from e
    try:
        snapshot = catalog.read(_sibling(source, CATALOG_NAME))
    except FileNotFoundError as e:
        raise ReleaseError(f"{source}: the release holds no {CATALOG_NAME}") from e
    pinned = manifest["files"][CATALOG_NAME]
    digest = hashlib.sha256(snapshot.getvalue()).hexdigest()
    if digest != pinned:
        raise ReleaseError(
            f"{source}: {CATALOG_NAME} has the SHA-256 {digest}, "
            f"and the manifest pins {pinned}"
        )
    try:
        frame = catalog.load(snapshot)
        rows = sorted(str(v) for v in frame["dataset_version_id"])
    except (KeyError, ValueError, OSError) as e:
        # pyarrow raises a ValueError or an OSError for content that is not Parquet.
        raise ReleaseError(f"{source}: {CATALOG_NAME} is not a catalog: {e!r}") from e
    if rows != sorted(store["dataset_version_id"] for store in manifest["stores"]):
        raise ReleaseError(
            f"{source}: the rows of {CATALOG_NAME} are not the stores of the manifest"
        )
    # A row that records no organism or no version gives None, which no manifest holds.
    # Each other value is compared as a string, which is how the release builder
    # writes it into the manifest.
    found: dict[str | None, set[str | None]] = {}
    for organism, version in zip(
        _strings(catalog.recorded(frame, "organism")),
        _strings(catalog.recorded(frame, "gene_reference_version")),
        strict=True,
    ):
        found.setdefault(organism, set()).add(version)
    versions = manifest["gene_reference_version"]
    if found != {organism: {version} for organism, version in versions.items()}:
        listing = {o: sorted(found[o], key=str) for o in sorted(found, key=str)}
        raise ReleaseError(
            f"{source}: the rows of {CATALOG_NAME} record the gene reference versions "
            f"{listing}, and the manifest gives {versions}"
        )
    return manifest, snapshot
