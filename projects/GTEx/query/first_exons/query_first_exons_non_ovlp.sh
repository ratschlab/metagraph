#!/bin/bash

set -e 

K=41
threads=16
### source paths
# set basedir
# set metagraph
source ../../paths.sh

query_dir=${basedir}/gtex/queries

graph=${basedir}/gtex/graphs/output_k${K}_trimmed_clean_graph_chunked/graph_merged_k${K}.dbg
annotation=${basedir}/gtex/output_k${K}_trimmed_clean.abundance.brwt.annodbg
query=${query_dir}/gencode.v30.first_exons_non_ovlp.fa

$metagraph query -v -p ${threads} --discovery-fraction 0.0 --query-mode matches -i ${graph} -a ${annotation} ${query} | tee ${query_dir}/gencode.v30.first_exons_non_ovlp.result.txt
