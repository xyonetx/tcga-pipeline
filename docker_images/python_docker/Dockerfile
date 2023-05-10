FROM ubuntu:jammy

SHELL ["/bin/bash", "-c"]

RUN apt-get update && \
    apt-get install -y \
    build-essential \
    git \
    curl \
    python3 \
    python3-pip


# Install all software under /opt/software:
RUN mkdir -p /opt/software 
RUN mkdir -p /opt/software/resources
RUN mkdir -p /opt/software/scripts

# Download the prior files into this Docker image:
RUN curl https://granddb.s3.amazonaws.com/optPANDA/ppi/ppi_complete.txt -o /opt/software/resources/tissues_ppi.tsv
RUN curl https://granddb.s3.amazonaws.com/optPANDA/motifs/Hugo_motifCellLine.txt -o /opt/software/resources/tissues_motif.tsv
RUN curl ftp://ftp.broadinstitute.org://pub/gsea/msigdb/human/annotations/Human_Ensembl_Gene_ID_MSigDB.v2023.1.Hs.chip -o /opt/software/resources/Human_Ensembl_Gene_ID_MSigDB.v2023.1.Hs.chip

# Add the various python requirements and script files
ADD requirements.txt /opt/software/

# Install python deps
RUN pip3 install -r /opt/software/requirements.txt

# Install the NetZoo python library
RUN cd /opt/software && \
    curl -sL https://github.com/netZoo/netZooPy/archive/refs/tags/0.9.14.tar.gz | tar -xz && \
    pip3 install -e /opt/software/netZooPy-0.9.14/

ADD extract_tcga_type.py /opt/software/scripts/
ADD run_panda.py /opt/software/scripts/
ADD run_lioness.py /opt/software/scripts/
ADD map_ensg_to_symbol.py /opt/software/scripts/
ADD expression_filter.py /opt/software/scripts/
ADD determine_lioness_splits.py /opt/software/scripts/
ADD merge_slices.py /opt/software/scripts/
ADD lioness_stats.py /opt/software/scripts/