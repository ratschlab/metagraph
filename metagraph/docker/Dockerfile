FROM ubuntu:20.04

RUN apt-get update && apt-get install -y \
    libatomic1 \
    libgomp1 \
    libcurl4-nss-dev \
    libhts-dev \
    libjemalloc2 \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

COPY docker/bin_tmp/metagraph_docker /usr/local/bin/metagraph

RUN mkdir /opt/metagraph
COPY api/python /opt/metagraph/api/python

RUN pip3 install /opt/metagraph/api/python

USER nobody

ENTRYPOINT ["metagraph"]
