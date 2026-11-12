FROM condaforge/mambaforge:latest

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
    nanoplot \
    porechop \
    filtlong \
    multiqc \
    seqkit \
    -y && \
    mamba clean --all -y

ENV PATH="/opt/conda/envs/nanopore_qc/bin:$PATH"

WORKDIR /data

COPY qc/run_pipeline.sh /opt/pipeline/run_pipeline.sh
RUN chmod +x /opt/pipeline/run_pipeline.sh

CMD ["bash"]
