# Dependencies, plus the census2-dataset-writer package and the census2-sdk
# package, which defines the obs schema that the writer imports.
#
# The pipeline's own scripts ship as Nextflow module binaries, so an edit to
# them needs no rebuild. The data model lives in this image. A change to
# census2-dataset-writer or census2-sdk therefore needs a rebuild, and for
# HealthOmics a push under a new tag.
#
# Build from the repository root: the context has to include the packages.
#   docker build -t census2-bulk-ingest:dev -f bulk-ingest-nf/docker/Dockerfile .
#
# HealthOmics does not support ARM containers. Build linux/amd64 for ECR:
#   docker buildx build --platform linux/amd64 \
#       -t census2-bulk-ingest:dev -f bulk-ingest-nf/docker/Dockerfile .
FROM python:3.12-slim

# aria2 is for tools/stage-refdata.py, which CI runs in this image. No task uses it.
RUN apt-get update \
 && apt-get install -y --no-install-recommends curl unzip ca-certificates procps aria2 \
 && rm -rf /var/lib/apt/lists/*

# AWS CLI v2. The task fetches each source h5ad with it. Writes to object storage
# go through obstore instead.
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

# Dependencies first, in their own layer, resolved from the packages'
# pyproject.toml files so that they are declared in exactly one place. The layer
# matters: these are ~400 MB against ~1 MB for the packages, so editing a
# package must not invalidate them.
#
# pip cannot install only a project's dependencies, so this stages the metadata
# of each package beside an empty module, installs that, and replaces it below.
# The writer names census2-sdk as a dependency, and pip takes it from the path
# given next to the writer's. The dependency set therefore stays declared in one
# place.
#
# The image therefore takes whatever satisfies the package's ranges at build
# time. The range that matters is anndata <0.14, and every receipt records the
# versions that wrote a store. Generate a lock file if a rebuild ever has to
# resolve identically.
COPY census2-sdk/pyproject.toml census2-sdk/README.md /src/census2-sdk/
COPY census2-dataset-writer/pyproject.toml census2-dataset-writer/README.md \
     /src/census2-dataset-writer/
RUN mkdir -p /src/census2-sdk/src/census2_sdk \
             /src/census2-dataset-writer/src/census2_dataset_writer \
 && touch /src/census2-sdk/src/census2_sdk/__init__.py \
          /src/census2-dataset-writer/src/census2_dataset_writer/__init__.py \
 && pip install --no-cache-dir /src/census2-sdk '/src/census2-dataset-writer[s3]'

# The real sources, replacing the empty modules. --no-deps because everything
# they need is installed above, so editing a package rebuilds one small layer and
# downloads nothing.
COPY census2-sdk /src/census2-sdk
COPY census2-dataset-writer /src/census2-dataset-writer
RUN pip install --no-cache-dir --no-deps --force-reinstall \
        /src/census2-sdk /src/census2-dataset-writer

# Fail fast if the import surface the pipeline relies on is missing. Keep this
# list matched to what census2_dataset_writer and cxg-ingest.py import.
RUN python -c "\
import census2_dataset_writer, census2_sdk.schema, anndata, zarr, h5py, scipy, numpy, pandas, dask; \
import obstore, obstore.store; \
import cellxgene_ontology_guide; \
import importlib.metadata as md; \
from anndata.experimental import read_elem_lazy; \
from anndata.io import read_elem; \
from zarr.storage import ObjectStore; \
print('census2-dataset-writer', md.version('census2-dataset-writer'), \
      '| census2-sdk', md.version('census2-sdk'), \
      '| anndata', md.version('anndata'), '| zarr', md.version('zarr'))"

WORKDIR /work
