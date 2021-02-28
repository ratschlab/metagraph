#!/bin/bash

set -e

K=31

mem=320000
threads=6
threads=12
pmem=$(($mem / $threads))
basedir=/cluster/work/grlab/projects/metagenome/data
seqdir=${basedir}/metagut/graphs/output_k${K}
outdir=${basedir}/metagut/graphs/output_k${K}_graph_chunked
mkdir -p $outdir
metagraph=/cluster/home/akahles/git/software/metagraph_current3/metagraph/build/metagraph

for F in {\\\$,A,C,G,T}{\\\$,A,C,G,T}{\\\$,A,C,G,T}
do
    FF=$(eval "echo $F")
    if [ -f "${outdir}/graph_merged_k${K}.${FF}.dbg.chunk" ]
    then
        echo chunk $FF exists
    else
        echo "/usr/bin/time -v ${metagraph} build -v --parallel ${threads} --mode canonical -k ${K} --mem-cap-gb $((${mem} / 2000)) -o ${outdir}/graph_merged_k${K} --suffix $F ${seqdir}/*.fasta.gz 2>&1 | tee ${outdir}/build_k${K}_$F.log" | bsub -J mg_g_${F} -oo ${outdir}/build_k${K}_$F.lsf.log -We 36:00 -n $threads -M $mem -R "rusage[mem=${pmem}]" -R "span[hosts=1]"
    fi
done

