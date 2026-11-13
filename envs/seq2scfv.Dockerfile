# seq2scfv with IgBLAST, seqkit and Clustal Omega, adapted from the Dockerfile in Mikolaj Kocikowski's fork of seq2scfv.
FROM ubuntu:22.04@sha256:74b12b11ffa550fba1929ad4d8d16533c34f6612728469671afba1b03b66881c

SHELL ["/bin/bash", "-c"]

ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
ENV PATH="/usr/miniconda/bin:$PATH"

ARG DEBIAN_FRONTEND=noninteractive
ARG VERSION_SEQKIT="2.3.1"
ARG VERSION_IGBLAST="1.20.0"
ARG VERSION_PANDAS="2.2.2"
ARG VERSION_BIOPYTHON="1.84"
ARG VERSION_CLUSTALO="1.2.4"
ARG VERSION_BIOAWK="1.0"
ARG VERSION_PYTHON="3.10"

RUN apt-get update -y && apt-get install -y unzip wget perl python3 make build-essential libbz2-dev zlib1g-dev liblzma-dev lzma liblzma-dev git libxml2

RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda3.sh
RUN bash miniconda3.sh -b -u -p /usr/miniconda
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda config --set plugins.auto_accept_tos yes

RUN conda install -y python=${VERSION_PYTHON} \
                     pandas=${VERSION_PANDAS} \
                     biopython=${VERSION_BIOPYTHON}

RUN conda install -y -c bioconda seqkit=${VERSION_SEQKIT} \
                                 clustalo=${VERSION_CLUSTALO} \
                                 bioawk=${VERSION_BIOAWK}

# the last AntPack version that does not require a licence key
RUN pip3 install logomaker==0.8.7 \
                 antpack==0.3.8.5

ARG PATH_IGBLAST="/igblast"
ENV IGDATA=${PATH_IGBLAST}/ncbi-igblast-${VERSION_IGBLAST}
WORKDIR ${PATH_IGBLAST}
RUN wget --no-check-certificate https://ftp.ncbi.nih.gov/blast/executables/igblast/release/${VERSION_IGBLAST}/ncbi-igblast-1.20.0-x64-linux.tar.gz &&\
    tar zxf ncbi-igblast-${VERSION_IGBLAST}-x64-linux.tar.gz &&\
    ln -t /usr/local/bin ncbi-igblast-${VERSION_IGBLAST}/bin/* &&\
    rm ncbi-igblast-${VERSION_IGBLAST}-x64-linux.tar.gz

ARG SEQ2SCFV_REPO
ARG SEQ2SCFV_COMMIT
RUN git clone ${SEQ2SCFV_REPO} /home/seq2scfv && git -C /home/seq2scfv checkout ${SEQ2SCFV_COMMIT}
RUN cp /home/seq2scfv/bin/* /usr/local/bin/
RUN chmod +x /usr/local/bin/*
