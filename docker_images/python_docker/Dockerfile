FROM ubuntu:jammy

SHELL ["/bin/bash", "-c"]

RUN apt-get update && \
    apt-get install -y \
    build-essential \
    wget \
    git \
    python3 \
    python3-pip

# Install all software under /opt/software:
RUN mkdir -p /opt/software 
RUN mkdir -p /opt/software/resources
ADD requirements.txt /opt/software/requirements.txt

# Download the prior files into this Docker image:
RUN wget https://granddb.s3.amazonaws.com/optPANDA/ppi/ppi_complete.txt -O /opt/software/resources/tissues_ppi.tsv
RUN wget https://granddb.s3.amazonaws.com/optPANDA/motifs/regMatQval005.csv -O /opt/software/resources/tissues_motif.tsv

# Install the NetZoo python library
RUN wget https://github.com/netZoo/netZooPy/archive/refs/tags/0.9.14.tar.gz \
    -P /opt/software && \
    tar xzvf /opt/software/0.9.14.tar.gz -C /opt/software && \
    pip3 install -e /opt/software/netZooPy-0.9.14/