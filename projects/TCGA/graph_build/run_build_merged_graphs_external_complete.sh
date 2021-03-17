#!/bin/bash

set -e

K=31

#mem=250000
mem=200000
threads=16
#threads=8
pmem=$(($mem / $threads))
basedir=/cluster/work/grlab/projects/metagenome/data
metagraph=/cluster/home/akahles/git/software/metagraph_current_dev/metagraph/build/metagraph

seqdir=${basedir}/tcga/output_k${K}_trimmed_clean_graph
outdir=${basedir}/tcga/
mkdir -p $outdir

tmpdir=${outdir}/tmp
mkdir -p $tmpdir

if [ -f "${outdir}/graph_merged_complete_k${K}.dbg" ]
then
    echo complete graph already exists
else
    echo "ls -1 ${seqdir}/*.fasta.gz | /usr/bin/time -v ${metagraph} build -v --parallel ${threads} --mode canonical -k ${K} --mem-cap-gb $((${mem} / 2000)) --container vector_disk --tmp-dir $tmpdir --disk-cap-gb 500 -o ${outdir}/graph_merged_complete_k${K} 2>&1 | tee ${outdir}/build_complete_k${K}.log" | bsub -G ms_raets -J tc_g_cplt -oo ${outdir}/build_complete_k${K}.lsf.log -We 12:00 -n $threads -M $mem -R "rusage[mem=${pmem}]" -R "span[hosts=1]"
fi
