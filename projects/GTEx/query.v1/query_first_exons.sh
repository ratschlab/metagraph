#!/bin/bash

set -e 

K=41
threads=30
metagraph=$(pwd)/../../../metagraph/build/metagengraph

basedir=/cluster/work/grlab/projects/metagenome/data/gtex
query_dir=${basedir}/queries

$metagraph query -p ${threads} --discovery-fraction 0.0 --count-labels -i ${basedir}/output_k${K}_merged/all_merged_k${K}.dbg -a ${basedir}/output_k${K}_merged.abundances.brwt.annodbg ${query_dir}/gencode.v30.first_exons.fa > ${query_dir}/gencode.v30.first_exons.result.txt
