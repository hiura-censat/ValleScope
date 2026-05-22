# Dockerfile (full, with Snakemake)
FROM mambaorg/micromamba:1.5.8

ENV MAMBA_ROOT_PREFIX=/home/mambauser/micromamba \
    MAMBA_NO_BANNER=1

SHELL ["/bin/bash", "-lc"]

USER root
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
      git build-essential zlib1g-dev libcurl4-openssl-dev tabix cmake curl wget ca-certificates && \
    rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/seryrzu/unialigner.git /opt/unialigner \
    && cd /opt/unialigner \
    && git checkout c5a1eecab7bd17485a0fe3422684409c3e884f31 \
    && cd /opt/unialigner/tandem_aligner \
    && grep -q '^#include <algorithm>' src/tools/suffix_array/suffix_array.hpp \
       || sed -i '1i#include <algorithm>' src/tools/suffix_array/suffix_array.hpp \
    && make -j"$(nproc)" \
    && test -x /opt/unialigner/tandem_aligner/build/bin/tandem_aligner \
    && ln -sf /opt/unialigner/tandem_aligner/build/bin/tandem_aligner /usr/local/bin/tandem_aligner \
    && command -v tandem_aligner

RUN git clone https://github.com/lh3/dna-nn /usr/local/bin/dna-nn \
    && cd /usr/local/bin/dna-nn \
    && make \
    && ln -sf /usr/local/bin/dna-nn/dna-brnn /usr/local/bin/dna-brnn

RUN curl -fsSL https://deb.nodesource.com/setup_lts.x | bash - \
    && apt-get update \
    && apt-get install -y --no-install-recommends nodejs \
    && npm install -g @jbrowse/cli@3.6.5 \
    && rm -rf /var/lib/apt/lists/*

USER mambauser
COPY --chown=mambauser:users envs/conda-linux-64.lock /tmp/conda.lock
RUN micromamba create -y -n base -f /tmp/conda.lock && micromamba clean -a -y

WORKDIR /opt/app
COPY --chown=mambauser:users . .
RUN micromamba run -n base python -m pip install .

ENTRYPOINT ["bash", "-lc", "sleep infinity"]
