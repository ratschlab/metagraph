#!/bin/bash

set -e 

K=41
threads=4
### source paths
# set basedir
# set metagraph
source ../../paths.sh

query_dir=${basedir}/gtex/queries
result_dir=${query_dir}/results
mkdir -p $result_dir

graph=${basedir}/gtex/graphs/output_k${K}_trimmed_clean_graph_chunked/graph_merged_k${K}.dbg
annotation=${basedir}/gtex/graphs/output_k${K}_trimmed_clean.samples.brwt.annodbg
query=${query_dir}/trans_exons/gencode.v30.trans_exons_all.result.hits.fa
$metagraph query -v -p ${threads} --discovery-fraction 0.8 --query-mode matches -i ${graph} -a ${annotation} ${query} | tee ${query_dir}/trans_exons/gencode.v30.trans_exons_all_hits.result_samples.txt

