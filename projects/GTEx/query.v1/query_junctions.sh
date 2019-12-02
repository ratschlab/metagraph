#!/bin/bash

set -e 

K=41
threads=20
threads=1
#metagraph=/cluster/home/akahles/git/software/metagraph_brwt_rrr/metagraph/build/metagengraph
metagraph=/cluster/home/akahles/git/software/metagraph_clean2/metagraph/build/metagraph

basedir=/cluster/work/grlab/projects/metagenome/data/gtex
query_dir=${basedir}/queries

$metagraph query -p ${threads} --discovery-fraction 0.0 --count-labels -i ${basedir}/output_k${K}_merged/all_merged_k${K}.dbg -a ${basedir}/output_k${K}_merged.abundances.brwt.annodbg ${query_dir}/gencode.v30.all_junctions.fa > ${query_dir}/gencode.v30.all_junctions.result.txt
