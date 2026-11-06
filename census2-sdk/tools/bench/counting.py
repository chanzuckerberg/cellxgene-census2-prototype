"""Count the object reads of one benchmark process, for every kind of object.

A query reads four kinds of object, through three clients:

- `zarr`: the store objects, through zarr and obstore. `CountingStore` wraps the zarr
  store and counts at the zarr store interface. One group of byte ranges that zarr
  fetches together counts as one request.
- `catalog` and `table_metadata`: the catalog file and the table's version hint,
  through `census2_sdk.catalog.read` and obstore. One read is one GET.
- `table_metadata` and `fragment`: the Iceberg metadata file, the manifest list, the
  manifests, and the `obs` fragments, through pyiceberg's `PyArrowFileIO` and the
  pyarrow S3 file system. `install_arrow_counting` makes pyiceberg open each file
  through `CountingHandler`, a `pyarrow.fs.FileSystemHandler` over the real file
  system. An open is one HEAD, because pyarrow's S3 file system asks for the size of
  each file that it opens. Each read of the open file is one GET, because pyarrow
  reads a range of an S3 object with one GET.

The file-system wrapper puts Python in each read of a fragment, and reads of one file
take turns. Thus `--count-arrow` can change the times of the table path, and the
benchmark measures times without it (see `run.py`).
"""

import threading
from collections.abc import AsyncIterator, Iterable, Sequence
from typing import Any

import pyarrow as pa  # type: ignore[import-untyped]  # pyarrow has no types
import pyarrow.fs as pafs  # type: ignore[import-untyped]
import zarr
import zarr.storage
from zarr.abc.store import ByteRequest, Store
from zarr.core.buffer import Buffer, BufferPrototype

KINDS = ("zarr", "fragment", "table_metadata", "catalog")


class Counter:
    """The requests and the bytes of each kind of object in this process."""

    def __init__(self) -> None:
        self.lock = threading.Lock()
        self.requests = dict.fromkeys(KINDS, 0)
        self.heads = dict.fromkeys(KINDS, 0)
        self.nbytes = dict.fromkeys(KINDS, 0)

    def add(self, kind: str, requests: int, nbytes: int, *, head: bool = False) -> None:
        with self.lock:
            self.requests[kind] += requests
            self.nbytes[kind] += nbytes
            if head:
                self.heads[kind] += requests

    def snapshot(self) -> dict[str, Any]:
        """The counts so far: the total, and the requests and MB of each kind."""
        with self.lock:
            return {
                "requests": sum(self.requests.values()),
                "requests_by_kind": dict(self.requests),
                "heads_by_kind": dict(self.heads),
                "mb_by_kind": {k: round(v / 2**20, 3) for k, v in self.nbytes.items()},
            }


COUNTER = Counter()


def kind_of(path: str) -> str:
    """The kind of a file that is not a store object, from its name."""
    if path.endswith(".obs.parquet"):
        return "fragment"
    if "/metadata/" in path or path.endswith((".avro", ".metadata.json")):
        return "table_metadata"
    return "catalog"


class CountingStore(zarr.storage.WrapperStore[Store]):
    """A zarr store that counts each read and its data in `COUNTER`.

    A sharded array reads parts of a shard through `get_ranges`. The wrapped store
    fetches each group of ranges together, so one group counts as one request. The
    bytes are the data that zarr keeps, without the gaps that a group fetch also reads.
    """

    async def get(
        self,
        key: str,
        prototype: BufferPrototype,
        byte_range: ByteRequest | None = None,
    ) -> Buffer | None:
        buffer = await self._store.get(key, prototype, byte_range)
        COUNTER.add("zarr", 1, 0 if buffer is None else len(buffer))
        return buffer

    async def get_partial_values(
        self,
        prototype: BufferPrototype,
        key_ranges: Iterable[tuple[str, ByteRequest | None]],
    ) -> list[Buffer | None]:
        buffers = await self._store.get_partial_values(prototype, key_ranges)
        sizes = [0 if b is None else len(b) for b in buffers]
        COUNTER.add("zarr", len(sizes), sum(sizes))
        return buffers

    async def get_ranges(
        self,
        key: str,
        byte_ranges: Sequence[ByteRequest | None],
        *,
        prototype: BufferPrototype,
        max_concurrency: int | None = None,
        max_gap_bytes: int | None = None,
        max_coalesced_bytes: int | None = None,
    ) -> AsyncIterator[Sequence[tuple[int, Buffer | None]]]:
        groups = super().get_ranges(
            key,
            byte_ranges,
            prototype=prototype,
            max_concurrency=max_concurrency,
            max_gap_bytes=max_gap_bytes,
            max_coalesced_bytes=max_coalesced_bytes,
        )
        async for group in groups:
            size = sum(0 if b is None else len(b) for _, b in group)
            COUNTER.add("zarr", 1, size)
            yield group


class CountingReader:
    """A Python file over a pyarrow random-access file. Each read counts one GET."""

    def __init__(self, native: Any, kind: str) -> None:
        self.native = native
        self.kind = kind
        self.position = 0
        self.closed = False

    def readable(self) -> bool:
        return True

    def seekable(self) -> bool:
        return True

    def writable(self) -> bool:
        return False

    def size(self) -> int:
        return int(self.native.size())

    def tell(self) -> int:
        return self.position

    def seek(self, position: int, whence: int = 0) -> int:
        base = {0: 0, 1: self.position, 2: self.size()}[whence]
        self.position = base + position
        return self.position

    def read(self, nbytes: int | None = None) -> bytes:
        remaining = self.size() - self.position
        n = remaining if nbytes is None or nbytes < 0 else min(nbytes, remaining)
        if n <= 0:
            return b""
        data = self.native.read_at(n, self.position)
        self.position += len(data)
        COUNTER.add(self.kind, 1, len(data))
        return bytes(data)

    def close(self) -> None:
        self.closed = True
        self.native.close()


class CountingHandler(pafs.FileSystemHandler):  # type: ignore[misc]
    """A pyarrow file system that counts the HEAD and the GETs of each file it reads."""

    def __init__(self, fs: Any) -> None:
        self.fs = fs

    def __eq__(self, other: object) -> bool:
        return isinstance(other, CountingHandler) and self.fs == other.fs

    def __ne__(self, other: object) -> bool:
        return not self == other

    def get_type_name(self) -> str:
        return "counting-" + str(self.fs.type_name)

    def normalize_path(self, path: str) -> str:
        return str(self.fs.normalize_path(path))

    def get_file_info(self, paths: list[str]) -> Any:
        for path in paths:
            COUNTER.add(kind_of(path), 1, 0, head=True)
        return self.fs.get_file_info(paths)

    def get_file_info_selector(self, selector: Any) -> Any:
        return self.fs.get_file_info(selector)

    def _open(self, path: str) -> Any:
        kind = kind_of(path)
        native = self.fs.open_input_file(path)
        COUNTER.add(kind, 1, 0, head=True)
        return pa.PythonFile(CountingReader(native, kind), mode="r")

    def open_input_file(self, path: str) -> Any:
        return self._open(path)

    def open_input_stream(self, path: str) -> Any:
        return self._open(path)

    def create_dir(self, path: str, recursive: bool) -> None:
        self.fs.create_dir(path, recursive=recursive)

    def delete_dir(self, path: str) -> None:
        self.fs.delete_dir(path)

    def delete_dir_contents(self, path: str, missing_dir_ok: bool = False) -> None:
        self.fs.delete_dir_contents(path, missing_dir_ok=missing_dir_ok)

    def delete_root_dir_contents(self) -> None:
        self.fs.delete_dir_contents("", accept_root_dir=True)

    def delete_file(self, path: str) -> None:
        self.fs.delete_file(path)

    def move(self, src: str, dest: str) -> None:
        self.fs.move(src, dest)

    def copy_file(self, src: str, dest: str) -> None:
        self.fs.copy_file(src, dest)

    def open_output_stream(self, path: str, metadata: Any = None) -> Any:
        return self.fs.open_output_stream(path, metadata=metadata)

    def open_append_stream(self, path: str, metadata: Any = None) -> Any:
        return self.fs.open_append_stream(path, metadata=metadata)


def install_arrow_counting() -> None:
    """Make pyiceberg open its files through `CountingHandler`. Call before a table opens."""
    from pyiceberg.io.pyarrow import PyArrowFileIO

    original = PyArrowFileIO._initialize_fs

    def counted(self: Any, scheme: str, netloc: str | None = None) -> Any:
        return pafs.PyFileSystem(CountingHandler(original(self, scheme, netloc)))

    PyArrowFileIO._initialize_fs = counted  # type: ignore[method-assign]


def install_read_counting() -> None:
    """Count the reads of `census2_sdk.catalog.read`: the catalog and the version hint."""
    from census2_sdk import catalog as _catalog

    original = _catalog.read

    def counted(uri: Any) -> Any:
        content = original(uri)
        COUNTER.add(kind_of(str(uri)), 1, len(content.getvalue()))
        return content

    # The SDK modules call `catalog.read` through the module, so one patch counts all.
    _catalog.read = counted
