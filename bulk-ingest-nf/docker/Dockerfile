# Runtime only. No project code lives here.
#
# The pipeline's own scripts ship as Nextflow module binaries, so day-to-day edits
# need no image rebuild. Rebuild this image only when a dependency changes.
#
# HealthOmics does not support ARM containers. Build linux/amd64 for ECR:
#   docker buildx build --platform linux/amd64 -t census2-bulk-ingest:dev -f docker/Dockerfile .
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

COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt && rm /tmp/requirements.txt

# Fail fast if the import surface the converter relies on is missing. Keep this
# list matched to what census2_store.py and cxg-ingest import.
RUN python -c "\
import anndata, zarr, h5py, scipy, numpy, dask, obstore, obstore.store; \
from anndata.experimental import read_elem_lazy; \
from anndata.io import read_elem; \
from zarr.storage import ObjectStore; \
print('anndata', anndata.__version__, '| zarr', zarr.__version__)"

WORKDIR /work
