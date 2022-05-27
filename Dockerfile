ARG CODE_BASE="/opt/metagraph"

FROM ubuntu:20.04 AS metagraph_dev_env

# contains all dependencies to build metagraph. Can also be used during development by mounting the code base and
# build dir on the host (this is done in `make build-metagraph env=docker`)

RUN export DEBIAN_FRONTEND="noninteractive" && apt-get update && apt-get install -y \
    autoconf \
    automake \
    binutils-dev \
    ccache \
    curl \
    g++-10 \
    git \
    libboost-all-dev \
    libbrotli-dev \
    libbz2-dev \
    libdouble-conversion-dev \
    libevent-dev \
    libgflags-dev \
    libgoogle-glog-dev \
    libiberty-dev \
    libjemalloc-dev \
    liblz4-dev \
    liblzma-dev \
    libsnappy-dev \
    libssl-dev \
    libtool \
    libunwind-dev \
    make \
    pkg-config \
    python3 \
    python3-pip \
    python3-venv \
    vim \
    wget \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install parameterized==0.7.1 "cmake>=3.19"

WORKDIR /opt
ENV LD_LIBRARY_PATH=/usr/local/lib
ENV CC /usr/bin/gcc-10
ENV CXX /usr/bin/g++-10

### Compile our own libcurl (compiling from source to avoid linking issues in static builds)
RUN git clone https://github.com/curl/curl.git \
    && mkdir curl/_build \
    && cd curl/_build \
    && cmake -DBUILD_SHARED_LIBS=on .. \
    && make -j $(nproc) \
    && make install \
    && cmake -DBUILD_SHARED_LIBS=off .. \
    && make -j $(nproc) \
    && make install \
    && cd /opt && rm -rf curl

### Installing htslib (compiling from source to avoid linking issues in static builds)
ARG HTSLIB_VERSION=1.10.2
RUN wget https://github.com/samtools/htslib/releases/download/$HTSLIB_VERSION/htslib-$HTSLIB_VERSION.tar.bz2 \
    && tar -vxjf htslib-$HTSLIB_VERSION.tar.bz2 \
    && cd htslib-$HTSLIB_VERSION \
    && make -j $(nproc) \
    && make install \
    && cd /opt && rm -rf htslib-$HTSLIB_VERSION htslib-${HTSLIB_VERSION}.tar.bz2

RUN mkdir -p /opt/metagraph/build_docker /opt/ccache_docker
RUN chmod o+rwx /opt/metagraph /opt/ccache_docker

ENV CCACHE_DIR=/opt/ccache_docker

FROM metagraph_dev_env as metagraph_bin
ARG CODE_BASE

COPY . ${CODE_BASE}

WORKDIR ${CODE_BASE}
RUN make build-sdsl-lite \
    && make build-metagraph alphabet=Protein \
    && make build-metagraph alphabet=DNA

FROM ubuntu:20.04
ARG CODE_BASE

# the image used in production. It contains a basic runtime environment for metagraph without build tools along with
# the metagraph binary and python API code. This image is published on github's container registry (`ghcr.io/ratschlab/metagraph`).

RUN apt-get update && apt-get install -y \
    libatomic1 \
    libcurl4-nss-dev \
    libgomp1 \
    libhts-dev \
    libjemalloc2 \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/* && pip3 install -U "pip>=22.1"

COPY --from=metagraph_bin ${CODE_BASE}/metagraph/build/metagraph_* /usr/local/bin/

RUN ln -s /usr/local/bin/metagraph_DNA /usr/local/bin/metagraph

RUN mkdir ${CODE_BASE}
COPY . ${CODE_BASE}

RUN pip3 install ${CODE_BASE}/metagraph/api/python
RUN pip3 install ${CODE_BASE}/metagraph/workflows

ENTRYPOINT ["metagraph"]
