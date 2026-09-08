# Dependencies, plus the census2-store package.
#
# Two kinds of code, with two different loops. The pipeline's own scripts ship as
# Nextflow module binaries, so editing them needs no rebuild. The data model does
# live in this image, so a change to census2-store means a rebuild, and for
# HealthOmics a push under a new tag.
#
# Build from the repository root: the context has to include the package.
#   docker build -t census2-bulk-ingest:dev -f bulk-ingest-nf/docker/Dockerfile .
#
# HealthOmics does not support ARM containers. Build linux/amd64 for ECR:
#   docker buildx build --platform linux/amd64 \
#       -t census2-bulk-ingest:dev -f bulk-ingest-nf/docker/Dockerfile .
FROM python:3.12-slim

RUN apt-get update \
 && apt-get install -y --no-install-recommends curl unzip ca-certificates procps \
 && rm -rf /var/lib/apt/lists/*

# AWS CLI v2. All S3 traffic goes through this; nothing in the Python stack speaks S3.
RUN set -eux; \
    case "$(uname -m)" in \
      x86_64)          A=x86_64 ;; \
      aarch64|arm64)   A=aarch64 ;; \
      *) echo "unsupported arch $(uname -m)" >&2; exit 1 ;; \
    esac; \
    curl -fsSL "https://awscli.amazonaws.com/awscli-exe-linux-${A}.zip" -o /tmp/awscliv2.zip; \
    unzip -q /tmp/awscliv2.zip -d /tmp; \
    /tmp/aws/install; \
    rm -rf /tmp/awscliv2.zip /tmp/aws; \
    aws --version

# Dependencies first, in their own layer, resolved from the package's
# pyproject.toml so that they are declared in exactly one place. The layer
# matters: these are ~400 MB against ~1 MB for the package, so editing the
# package must not invalidate them.
#
# pip cannot install only a project's dependencies, so this stages the metadata
# beside an empty module, installs that, and replaces it below. The alternative
# is a second dependency list to keep in step, which is what this removes.
#
# The image therefore takes whatever satisfies the package's ranges at build
# time. The range that matters is anndata <0.14, and every receipt records the
# versions that wrote a store. Generate a lock file if a rebuild ever has to
# resolve identically.
COPY census2-store/pyproject.toml census2-store/README.md /src/census2-store/
RUN mkdir -p /src/census2-store/src/census2_store \
 && touch /src/census2-store/src/census2_store/__init__.py \
 && pip install --no-cache-dir '/src/census2-store[s3]'

# The real source, replacing the empty module. --no-deps because everything it
# needs is installed above, so editing the package rebuilds one small layer and
# downloads nothing.
COPY census2-store /src/census2-store
RUN pip install --no-cache-dir --no-deps --force-reinstall /src/census2-store

# Fail fast if the import surface the pipeline relies on is missing. Keep this
# list matched to what census2_store and cxg-ingest.py import.
RUN python -c "\
import census2_store, anndata, zarr, h5py, scipy, numpy, pandas, dask, obstore, obstore.store; \
import importlib.metadata as md; \
from anndata.experimental import read_elem_lazy; \
from anndata.io import read_elem; \
from zarr.storage import ObjectStore; \
print('census2-store', md.version('census2-store'), '| anndata', md.version('anndata'), \
      '| zarr', md.version('zarr'))"

WORKDIR /work
