ARG CODE_BASE="/opt/metagraph"

FROM ubuntu:22.04 AS metagraph_dev_env

# contains all dependencies to build metagraph. Can also be used during development by mounting the code base and
# build dir on the host (this is done in `make build-metagraph env=docker`)

RUN export DEBIAN_FRONTEND="noninteractive" && apt-get update && apt-get install -y \
    autoconf \
    automake \
    ccache \
    g++-11 \
    git \
    cmake \
    make \
    libboost-all-dev \
    libbz2-dev \
    libdeflate-dev \
    libjemalloc-dev \
    liblzma-dev \
    libzstd-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt
ENV LD_LIBRARY_PATH=/usr/local/lib
ENV CC /usr/bin/gcc-11
ENV CXX /usr/bin/g++-11

RUN mkdir -p /opt/metagraph/build_docker /opt/ccache_docker
RUN chmod o+rwx /opt/metagraph /opt/ccache_docker

ENV CCACHE_DIR=/opt/ccache_docker

FROM metagraph_dev_env as metagraph_bin
ARG CODE_BASE

COPY . ${CODE_BASE}

WORKDIR ${CODE_BASE}
RUN make build-sdsl-lite \
    && make build-metagraph alphabet=DNA \
    && make build-metagraph alphabet=DNA5 \
    && make build-metagraph alphabet=Protein

FROM ubuntu:24.04
ARG CODE_BASE

# the image used in production. It contains a basic runtime environment for metagraph without build tools along with
# the metagraph binary and python API code. This image is published on github's container registry (`ghcr.io/ratschlab/metagraph`).

RUN apt-get update && apt-get install -y \
    libdeflate-dev \
    libzstd-dev \
    libjemalloc2 \
    python3 \
    python3-pip \
    python3-venv \
    && rm -rf /var/lib/apt/lists/*

COPY --from=metagraph_bin ${CODE_BASE}/metagraph/build/metagraph_* /usr/local/bin/

RUN ln -s /usr/local/bin/metagraph_DNA /usr/local/bin/metagraph

RUN mkdir ${CODE_BASE}
COPY . ${CODE_BASE}

RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:${PATH}"

RUN pip3 install ${CODE_BASE}/metagraph/api/python
RUN pip3 install ${CODE_BASE}/metagraph/workflows

# check that it runs fine
RUN metagraph --version \
    && metagraph_DNA --version \
    && metagraph_DNA5 --version \
    && metagraph_Protein --version

ENTRYPOINT []
