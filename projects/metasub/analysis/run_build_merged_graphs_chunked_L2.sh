#!/bin/bash

set -e

K=19

mem=320000
threads=6
pmem=$(($mem / $threads))
basedir=/cluster/work/grlab/projects/metagenome/data
seqdir=${basedir}/metasub/graphs/output_k${K}_cleaned
outdir=${basedir}/metasub/graphs/output_k${K}_cleaned_graph_chunked
mkdir -p $outdir
metagraph=/cluster/home/akahles/git/software/metagraph_clean/metagraph/build/metagraph

for F in {\\\$,A,C,G,T}{\\\$,A,C,G,T}
do
    echo "/usr/bin/time -v ${metagraph} build -v --parallel ${threads} --mode canonical -k ${K} --mem-cap-gb $((${mem} / 2000)) -o ${outdir}/graph_merged_k${K} --suffix $F ${seqdir}/*.fasta.gz 2>&1 | tee ${outdir}/build_k${K}_$F.log" | bsub -J ms_g_${F} -o ${outdir}/build_k${K}_$F.lsf.log -We 8:00 -n $threads -M $mem -R "rusage[mem=${pmem}]" -R "span[hosts=1]"
done

