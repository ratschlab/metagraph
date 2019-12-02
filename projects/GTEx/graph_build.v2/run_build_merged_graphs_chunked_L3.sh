#!/bin/bash

set -e

K=41

mem=250000
threads=6
pmem=$(($mem / $threads))
basedir=/cluster/work/grlab/projects/metagenome/data
seqdir=${basedir}/gtex/output_k${K}_clean
outdir=${basedir}/gtex/output_k${K}_clean_graph_chunked
mkdir -p $outdir
metagraph=/cluster/home/akahles/git/software/metagraph_current/metagraph/build/metagraph

for F in {\\\$,A,C,G,T}{\\\$,A,C,G,T}{\\\$,A,C,G,T}
do
    echo "/usr/bin/time -v ${metagraph} build -v --parallel ${threads} --canonical -k ${K} --mem-cap-gb $((${mem} / 2000)) -o ${outdir}/graph_merged_k${K} --suffix $F ${seqdir}/*.fasta.gz 2>&1 | tee ${outdir}/build_k${K}_$F.log" | bsub -J gt_g_${F} -o ${outdir}/build_k${K}_$F.lsf.log -We 8:00 -n $threads -M $mem -R "rusage[mem=${pmem}]" -R "span[hosts=1]"
done

