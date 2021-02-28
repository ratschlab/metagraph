#!/bin/bash

set -e

## source paths
# set metagraph
# set basedir
. ../paths.sh

if [ -z "$2" ]
then
    echo "Usage: $0 <K> <level>"
    exit 1
fi
K="$1"
level="$2"

mem=300000
threads=2
pmem=$(($mem / $threads))
graphdir=${basedir}/gtex/graphs/output_k${K}_trimmed_clean_graph_chunked

/usr/bin/time -v ${metagraph} concatenate -v --mode canonical --len-suffix $level -i ${graphdir}/graph_merged_k${K} -o ${graphdir}/graph_merged_k${K}

