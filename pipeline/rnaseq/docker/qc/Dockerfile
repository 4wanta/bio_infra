FROM mambaorg/micromamba:1.5.6

RUN micromamba install -y -n base -c bioconda -c conda-forge \
    fastqc=0.12.1 \
    trimmomatic \
    && micromamba clean --all --yes

USER root
RUN apt-get update && apt-get install -y procps default-jre-headless curl unzip && rm -rf /var/lib/apt/lists/* \
    && curl -fsSL "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o /tmp/awscliv2.zip \
    && unzip -q /tmp/awscliv2.zip -d /tmp \
    && /tmp/aws/install \
    && rm -rf /tmp/aws /tmp/awscliv2.zip
USER $MAMBA_USER

ENV PATH="/opt/conda/bin:${PATH}"
CMD ["/bin/bash"]
