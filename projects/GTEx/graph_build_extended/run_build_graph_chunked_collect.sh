#!/bin/bash

set -e

if [ -z "$2" ]
then
    echo "Usage: $0 <K> <level>"
    exit 1
fi
K="$1"
level="$2"

basedir=/cluster/work/grlab/projects/metagenome/data
graphdir=${basedir}/gtex/graphs/output_k${K}_trimmed_clean_graph_extended_chunked
metagraph=/cluster/home/akahles/git/software/metagraph_current/metagraph/build/metagraph

/usr/bin/time -v ${metagraph} concatenate -v -p 10 --mode canonical --len-suffix $level --clear-dummy -i ${graphdir}/graph_k${K} -o ${graphdir}/graph_k${K}

