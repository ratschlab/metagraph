#!/bin/bash

set -e 

K=41
threads=16
### source paths
# basedir = 
# metagraph =
source ../paths.sh

query_dir=${basedir}/gtex/queries
result_dir=${query_dir}/results
mkdir -p $result_dir

$metagraph query -p ${threads} --discovery-fraction 0.0 --count-labels -i ${basedir}/gtex/output_k${K}_trimmed_clean_graph_chunked/graph_merged_k${K}.dbg -a ${basedir}/gtex/output_k${K}_trimmed_clean.abundance.brwt.annodbg ${query_dir}/gencode.v30.all_junctions.fa | gzip > ${query_dir}/gencode.v30.all_junctions.result.txt.gz
