FROM ubuntu:20.04

RUN apt-get update && apt-get install -y \
    libatomic1 \
    libgomp1 \
    libcurl4-nss-dev \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# htslib
ENV LD_LIBRARY_PATH=/usr/local/lib

COPY docker/bin_tmp/metagraph_docker /usr/local/bin/metagraph

COPY docker/bin_tmp/libhts_docker.so /usr/local/lib/libhts.so
RUN ln -s /usr/local/lib/libhts.so /usr/local/lib/libhts.so.3

RUN mkdir /opt/metagraph
COPY api/python /opt/metagraph/api/python

RUN pip3 install /opt/metagraph/api/python

USER nobody

ENTRYPOINT ["metagraph"]
