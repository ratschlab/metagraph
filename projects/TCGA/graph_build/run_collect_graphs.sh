#!/bin/bash

set -e

if [ -z "$3" ]
then
    echo "Usage: $0 <K> <level> <tissue>"
    exit 1
fi
K="$1"
level="$2"
tissue="$3"

#mem=300000
#threads=2
#pmem=$(($mem / $threads))
basedir=/cluster/work/grlab/projects/metagenome/data
graphdir=${basedir}/tcga/output_k${K}_trimmed_clean_graph_chunked/${tissue}
metagraph=/cluster/home/akahles/git/software/metagraph_current_dev/metagraph/build/metagraph

/usr/bin/time -v ${metagraph} concatenate -v -p 10 --mode canonical --len-suffix $level --clear-dummy -i ${graphdir}/graph_merged_k${K} -o ${graphdir}/graph_merged_k${K}

