"""Run `bench.run`, and add the CPU time and the S3 requests of pyarrow to its record.

`bench.run` counts the reads of the stores at the zarr store interface. On the table
path, pyiceberg reads the fragments through pyarrow's S3 file system, which that count
does not see. This wrapper takes the same arguments as `bench.run`, and adds:

- `cpu_seconds`: the user and system CPU time of the process during `bench.run.run`.
- With `--count-arrow-requests`: `arrow_requests`, the requests of pyarrow's S3 client
  by method, and `arrow_mb`, the `Content-Length` of their responses. The count comes
  from the debug log of the AWS SDK for C++, which pyarrow uses. That log is slow and
  holds credentials, so it goes to a private temporary file, which the wrapper reads
  and deletes. Use the counts of such a run, and the times of a run without it.

    python -m bench.counted catalog-103.parquet --shape query --query '{...}' \\
        --obs-table s3://my-bucket/obs-103 --workers 8 --count-arrow-requests
"""

import os
import re
import resource
import sys
import tempfile
from typing import Any

from bench import run as bench_run

_REQUEST = re.compile(r"\(HeaderOut\) (GET|HEAD|PUT|POST|DELETE) ")
_LENGTH = re.compile(r"\(HeaderIn\) Content-Length: (\d+)", re.IGNORECASE)


def cpu_seconds() -> float:
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return usage.ru_utime + usage.ru_stime


def count_log(path: str) -> tuple[dict[str, int], float]:
    """The requests by method in an AWS SDK debug log, and the MB of the responses."""
    counts: dict[str, int] = {}
    nbytes = 0
    with open(path, errors="replace") as fh:
        for line in fh:
            found = _REQUEST.search(line)
            if found:
                counts[found.group(1)] = counts.get(found.group(1), 0) + 1
                continue
            length = _LENGTH.search(line)
            if length:
                nbytes += int(length.group(1))
    return counts, round(nbytes / 2**20, 2)


def main() -> int:
    argv = sys.argv[1:]
    count = "--count-arrow-requests" in argv
    sys.argv = [sys.argv[0], *(a for a in argv if a != "--count-arrow-requests")]
    log = None
    saved = None
    if count:
        handle, log = tempfile.mkstemp(prefix="aws-sdk-", suffix=".log")
        os.chmod(log, 0o600)
        saved = os.dup(2)
        os.dup2(handle, 2)
        os.close(handle)
        import pyarrow.fs as pafs  # type: ignore[import-untyped]

        pafs.initialize_s3(pafs.S3LogLevel.Debug)

    timed: dict[str, float] = {}
    inner_run, inner_publish = bench_run.run, bench_run.publish

    def run(*args: Any, **kwargs: Any) -> dict[str, Any]:
        start = cpu_seconds()
        record = inner_run(*args, **kwargs)
        timed["cpu"] = cpu_seconds() - start
        return record

    def publish(record: dict[str, Any], exit_code: int | None = None) -> bool:
        """Add the CPU time and the counts to the terminal record of `bench.run`.

        The memory guard publishes with an exit code and stops the process at once,
        so its record goes out as it is.
        """
        if exit_code is None:
            record = dict(record)
            if "cpu" in timed:
                record["cpu_seconds"] = round(timed["cpu"], 2)
            if log is not None:
                counts, mb = count_log(log)
                record.update(arrow_requests=counts, arrow_mb=mb, count_run=True)
                total = sum(counts.values())
                record["requests_total"] = int(record.get("requests", 0)) + total
        return inner_publish(record, exit_code)

    bench_run.run = run
    bench_run.publish = publish
    try:
        return bench_run.main()
    finally:
        if saved is not None:
            os.dup2(saved, 2)
            os.close(saved)
        if log is not None:
            os.unlink(log)


if __name__ == "__main__":
    raise SystemExit(main())
