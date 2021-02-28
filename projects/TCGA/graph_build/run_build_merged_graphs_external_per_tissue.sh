#!/bin/bash

set -e

if [ -z "$1" ]
then
    echo "Usage: $0 <tissue>"
    exit 1
fi
tissue=$1

K=31

#mem=250000
mem=200000
threads=16
#threads=8
pmem=$(($mem / $threads))
basedir=/cluster/work/grlab/projects/metagenome/data

if [ $tissue  == "all" ]
then
    tissues=$(find ${basedir}/tcga/output_k31_trimmed_clean -maxdepth 1 -mindepth 1 -type d | rev | cut -f 1 -d '/' | rev | tr '\n' ' ')
else 
    tissues=$tissue
fi

for tissue in $tissues
do
    seqdir=${basedir}/tcga/output_k${K}_trimmed_clean/${tissue}
    outdir=${basedir}/tcga/output_k${K}_trimmed_clean_graph
    mkdir -p $outdir
    metagraph=/cluster/home/akahles/git/software/metagraph_current_dev/metagraph/build/metagraph

    tmpdir=${outdir}/tmp/${tissue}
    mkdir -p $tmpdir

    if [ ! -f ${seqdir}_all_files.txt ]
    then
        find ${seqdir} -name \*.fasta.gz > ${seqdir}_all_files.txt
    fi

    if [ -f "${outdir}/graph_merged_${tissue}_k${K}.dbg" ]
    then
        echo $tissue graph already exists
    else
        echo submitting $tissue
        echo "cat ${seqdir}_all_files.txt | /usr/bin/time -v ${metagraph} build -v --parallel ${threads} --mode canonical -k ${K} --mem-cap-gb $((${mem} / 2000)) --container vector_disk --tmp-dir $tmpdir --disk-cap-gb 500 -o ${outdir}/graph_merged_${tissue}_k${K} 2>&1 | tee ${outdir}/build_${tissue}_k${K}.log" | bsub -G ms_raets -J tc_g_${tissue} -oo ${outdir}/build_${tissue}_k${K}.lsf.log -We 12:00 -n $threads -M $mem -R "rusage[mem=${pmem}]" -R "span[hosts=1]"
    fi
done
