#!/bin/bash

set -e

if [ -z "$2" ]
then
    echo "Usage: $0 <K> <level>"
    exit 1
fi
K="$1"
level="$2"

#mem=300000
#threads=2
#pmem=$(($mem / $threads))
basedir=/cluster/work/grlab/projects/metagenome/data
graphdir=${basedir}/metagut/graphs/output_k${K}_graph_chunked
metagraph=/cluster/home/akahles/git/software/metagraph_current3/metagraph/build/metagraph

/usr/bin/time -v ${metagraph} concatenate -v -p 12 --mode canonical --len-suffix $level --clear-dummy -i ${graphdir}/graph_merged_k${K} -o ${graphdir}/graph_merged_k${K}

