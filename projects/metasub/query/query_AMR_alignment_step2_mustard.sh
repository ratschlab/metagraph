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
outdir=${basedir}/metasub/queries/
mkdir -p $outdir

metagraph=/cluster/home/akahles/git/software/metagraph_current/metagraph/build/metagraph

threads=8
mem=600000
pmem=$((${mem} / ${threads}))

log=${outdir}/amr_mustard_all_ard_k${K}.align_step2.log
out=${outdir}/amr_mustard_all_ard_k${K}.align_step2.tsv
query=${outdir}/amr_mustard_all_ard_k${K}.fa
if [ ! -f ${query} ]
then
    python tsv2fasta.py ${query%fa}tsv > $query
fi
echo "/usr/bin/time -v $metagraph query -i $graph -a ${annotation} -p ${threads} --query-mode matches --min-kmers-fraction-label 0.1 ${query} > $out" | bsub -J ms_q_amr${K} -oo ${log} -We 24:00 -n $threads -M ${mem} -R "rusage[mem=${pmem}]" -R "span[hosts=1]"
