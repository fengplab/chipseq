# Allo 1.2.0 with the standalone Keras 3 package that TensorFlow >= 2.16 needs for `tensorflow.keras`.
# The bioconda/biocontainers image (allo:1.2.0--pyhdfd78af_0) lacks `keras` and fails with
#   ModuleNotFoundError: No module named 'tensorflow.keras'
#
# Build:   docker build -t allo:1.2.0-keras containers/allo
# Use:     nextflow run ... --with_allo --allo_container allo:1.2.0-keras
FROM python:3.12-slim

LABEL org.opencontainers.image.description="Allo 1.2.0 (seqcode/allo) with TensorFlow 2.21 + Keras 3 for nf-core/chipseq"

# procps provides `ps`, which Nextflow needs for task metrics
RUN apt-get update \
    && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*

# Versions verified together on test data (same allocation results as a working reference install)
RUN pip install --no-cache-dir \
        "bio-allo==1.2.0" \
        "tensorflow==2.21.0" \
        "keras==3.15.1" \
    && allo --version \
    && python -c "import tensorflow.keras, Allo; print('Allo import OK')"

ENV TF_CPP_MIN_LOG_LEVEL=2
