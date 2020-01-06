#!/bin/bash

set -e

if [ -z "$1" ]
then
    echo "Usage: $0 <k>"
    exit 1
fi
K=$1

basedir=/cluster/work/grlab/projects/metagenome/data
annobase=${basedir}/metasub/graphs/output_k${K}_cleaned_graph_annotation_clean.collect.relabeled
graph=${basedir}/metasub/graphs/output_k${K}_cleaned_graph_chunked/graph_merged_k${K}.dbg
outdir=${basedir}/metasub/queries/
mkdir -p $outdir
query=/cluster/work/grlab/projects/metagenome/raw_data/AMR/mustard/all_ard.fna

metagraph=/cluster/home/akahles/git/software/metagraph_current/metagraph/build/metagraph

threads=8
mem=200000
pmem=$((${mem} / ${threads}))

log=${outdir}/amr_mustard_all_ard_k${K}.log
out=${outdir}/amr_mustard_all_ard_k${K}.tsv
echo "/usr/bin/time -v $metagraph align -i $graph --fwd-and-reverse -o ${out} -p ${threads} ${query} 2>&1" | bsub -J ms_q_amr${K} -oo ${log} -We 24:00 -n $threads -M ${mem} -R "rusage[mem=${pmem}]" -R "span[hosts=1]":
