#!/bin/bash

set -e

if [ -z "$1" ]
then
    echo "Usage: $0 <k>"
    exit 1
fi
K=$1

basedir=/cluster/work/grlab/projects/metagenome/data
annotation=${basedir}/metasub/graphs/output_k${K}_cleaned_graph_annotation_clean.collect.relabeled.brwt.annodbg
graph=${basedir}/metasub/graphs/output_k${K}_cleaned_graph_chunked/graph_merged_k${K}.dbg
outdir=${basedir}/metasub/queries_fast
mkdir -p $outdir
query=/cluster/work/grlab/projects/metagenome/raw_data/plasmids/plasmids_db_v2/All_plasmids.fna

metagraph=/cluster/home/akahles/git/software/metagraph_current_dev/metagraph/build/metagraph

threads=2
mem=700000
pmem=$((${mem} / ${threads}))

log=${outdir}/plasmids_db_v2_k${K}.fast.log
out=${outdir}/plasmids_db_v2_k${K}.fast.tsv
echo "/usr/bin/time -v $metagraph query --align --suppress-unlabeled -i $graph -a ${annotation} -p ${threads} --query-mode matches --min-kmers-fraction-label 0.1 ${query} > $out" | bsub -J ms_q_plsm${K} -oo ${log} -We 24:00 -n $threads -M ${mem} -R "rusage[mem=${pmem}]" -R "span[hosts=1]"
