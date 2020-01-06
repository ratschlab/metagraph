#!/bin/bash

set -e

K=41

mem=50000
threads=6
pmem=$(($mem / $threads))
basedir=/cluster/work/grlab/projects/metagenome/data
outdir=${basedir}/gtex/graphs/output_k${K}_trimmed_clean_graph_chunked
mkdir -p $outdir
metagraph=/cluster/home/akahles/git/software/metagraph_current/metagraph/build/metagraph

echo "/usr/bin/time -v ${metagraph} transform -v --to-fasta --parallel ${threads} -o ${outdir}/graph_merged_k${K} ${outdir}/graph_merged_k${K}.dbg 2>&1 | tee ${outdir}/graph_to_fasta.k${K}.log" | bsub -J gt2fa -oo ${outdir}/grapp_to_fasta.k${K}.lsf.log -We 8:00 -n $threads -M $mem -R "rusage[mem=${pmem}]" -R "span[hosts=1]"

