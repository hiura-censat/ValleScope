# Dockerfile (full, with Snakemake)
FROM mambaorg/micromamba:1.5.8

ENV MAMBA_ROOT_PREFIX=/home/mambauser/micromamba \
    MAMBA_NO_BANNER=1

SHELL ["/bin/bash", "-lc"]

USER root
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
      git build-essential zlib1g-dev libcurl4-openssl-dev tabix curl wget && \
    rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/lh3/dna-nn /usr/local/bin/dna-nn \
    && cd /usr/local/bin/dna-nn \
    && make \
    && ln -s /usr/local/bin/dna-nn/dna-brnn /usr/local/bin/dna-brnn

RUN curl -fsSL https://deb.nodesource.com/setup_lts.x | bash - \
    && apt-get install -y nodejs \
    && npm install -g @jbrowse/cli@3.6.5

USER mambauser
COPY --chown=mambauser:users envs/conda-linux-64.lock /tmp/conda.lock
RUN micromamba create -y -n base -f /tmp/conda.lock && micromamba clean -a -y

WORKDIR /opt/app
COPY --chown=mambauser:users . .
RUN micromamba run -n base python -m pip install .

ENTRYPOINT ["bash", "-lc", "sleep infinity"]

