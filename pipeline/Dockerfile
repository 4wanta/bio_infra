FROM mambaorg/micromamba:1.5.6

# awscli を追加。conda-forge チャンネルから最新版を取得します
RUN micromamba install -y -n base -c bioconda -c conda-forge \
    fastp=0.23.4 \
    bwa=0.7.17 \
    samtools=1.15.1 \
    picard=3.1.1 \
    homer=4.11 \
    macs3=3.0.1 \
    awscli \
    && micromamba clean --all --yes

USER root
# Nextflowの動作に必要な procps も忘れずに
RUN apt-get update && apt-get install -y procps && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

ENV PATH="/opt/conda/bin:${PATH}"
CMD ["/bin/bash"]
