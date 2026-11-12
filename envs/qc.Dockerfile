FROM condaforge/mambaforge:latest@sha256:050753d138b6708128c9bb45055a9d943f06363f32f63a96f7dd5f304ada52a6

# Proxy settings (passed from build args)
ARG HTTP_PROXY
ARG HTTPS_PROXY
ARG http_proxy
ARG https_proxy
ARG NO_PROXY
ARG no_proxy
ENV HTTP_PROXY=${HTTP_PROXY} \
    HTTPS_PROXY=${HTTPS_PROXY} \
    http_proxy=${http_proxy} \
    https_proxy=${https_proxy} \
    NO_PROXY=${NO_PROXY} \
    no_proxy=${no_proxy}

RUN mamba create -n nanopore_qc -c bioconda -c conda-forge \
    nanoplot=1.48.0 \
    porechop=0.2.4 \
    filtlong=0.3.1 \
    multiqc=1.25.1 \
    seqkit=2.14.0 \
    -y && \
    mamba clean --all -y

ENV PATH="/opt/conda/envs/nanopore_qc/bin:$PATH"

WORKDIR /data

CMD ["bash"]
