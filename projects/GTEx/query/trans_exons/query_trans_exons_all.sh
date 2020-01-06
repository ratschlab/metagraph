#!/bin/bash

set -e 

K=41
threads=16
### source paths
# set basedir
# set metagraph
source ../../paths.sh

query_dir=${basedir}/gtex/queries
result_dir=${query_dir}/results
mkdir -p $result_dir

graph=${basedir}/gtex/graphs/output_k${K}_trimmed_clean_graph_chunked/graph_merged_k${K}.dbg
query=${query_dir}/trans_exons/gencode.v30.trans_exons_all.fa
$metagraph align -v -p ${threads} --fwd-and-reverse -i ${graph} ${query} | tee ${query_dir}/trans_exons/gencode.v30.trans_exons_all.result.txt
