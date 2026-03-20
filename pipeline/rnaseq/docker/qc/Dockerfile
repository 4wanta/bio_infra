FROM mambaorg/micromamba:1.5.6

RUN micromamba install -y -n base -c bioconda -c conda-forge \
    fastqc=0.12.1 \
    && micromamba clean --all --yes

USER root
RUN apt-get update && apt-get install -y procps && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

ENV PATH="/opt/conda/bin:${PATH}"
CMD ["/bin/bash"]
