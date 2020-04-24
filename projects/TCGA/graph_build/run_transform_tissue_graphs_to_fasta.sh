#!/bin/bash

set -e

if [ -z "$1" ]
then
    echo "Usage: $0 <tissue>"
    exit 1
fi
tissue=$1

K=31

mem=5000
threads=1
pmem=$(($mem / $threads))
basedir=/cluster/work/grlab/projects/metagenome/data
metagraph=/cluster/home/akahles/git/software/metagraph_current_dev/metagraph/build/metagraph

if [ $tissue  == "all" ]
then
    tissues=$(find ${basedir}/tcga/output_k31_trimmed_clean -maxdepth 1 -mindepth 1 -type d | rev | cut -f 1 -d '/' | rev | tr '\n' ' ')
else 
    tissues=$tissue
fi

for tissue in $tissues
do
    seqdir=${basedir}/tcga/output_k${K}_trimmed_clean/${tissue}
    graphdir=${basedir}/tcga/output_k${K}_trimmed_clean_graph
    outdir=$graphdir

    if [ ! -f "${outdir}/graph_merged_${tissue}_k${K}.dbg" ]
    then
        echo $tissue graph does not exist yet
    else
        echo submitting $tissue
        echo "/usr/bin/time -v ${metagraph} transform -v --to-fasta --primary-kmers -o ${outdir}/graph_merged_${tissue}_k${K} ${outdir}/graph_merged_${tissue}_k${K}.dbg 2>&1 | tee ${outdir}/dump_fasta_${tissue}_k${K}.log" | bsub -G ms_raets -J tc_g_${tissue} -oo ${outdir}/dump_fasta_${tissue}_k${K}.lsf.log -We 12:00 -n $threads -M $mem -R "rusage[mem=${pmem}]" -R "span[hosts=1]"
    fi
done
