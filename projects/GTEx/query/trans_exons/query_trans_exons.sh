#!/bin/bash

set -e 

K=41
threads=16
### source paths
# set basedir
# set metagraph
source ../../paths.sh

query_dir=${basedir}/gtex/queries

graph=${basedir}/gtex/output_k${K}_trimmed_clean_graph_chunked/graph_merged_k${K}.dbg
query=${query_dir}/gencode.v30.trans_exons.fa
$metagraph align -v -p ${threads} --fwd-and-reverse -i ${graph} ${query} | tee ${query_dir}/trans_exons/gencode.v30.trans_exons.result.txt
