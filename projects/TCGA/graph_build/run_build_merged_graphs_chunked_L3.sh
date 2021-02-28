#!/bin/bash

set -e

if [ -z "$1" ]
then
    echo "Usage: $0 <tissue>"
    exit 1
fi
tissue=$1

K=31

mem=250000
mem=200000
threads=18
threads=8
pmem=$(($mem / $threads))
basedir=/cluster/work/grlab/projects/metagenome/data
seqdir=${basedir}/tcga/output_k${K}_trimmed_clean/${tissue}
outdir=${basedir}/tcga/output_k${K}_trimmed_clean_graph_chunked/${tissue}
mkdir -p $outdir
metagraph=/cluster/home/akahles/git/software/metagraph_current_dev/metagraph/build/metagraph

if [ ! -f ${seqdir}_all_files.txt ]
then
    find ${seqdir} -name \*.fasta.gz > ${seqdir}_all_files.txt
fi

for F in {\\\$,A,C,G,T}{\\\$,A,C,G,T}{\\\$,A,C,G,T}
do
    ### get rid of all combinations that are not allowed
    if [ -z $(echo $F | grep -v -E '[ACGT]\\\$(\\\$)+$' | grep -v -E '(\\\$)+[ACGT]+(\\\$)+$' | grep -v -E '[ACGT]+\\\$(\\\$)*[ACGT]+') ]
    then
        continue
    fi

    FF=$(eval "echo $F")
    if [ -f "${outdir}/graph_merged_k${K}.${FF}.dbg.chunk" ]
    then
        echo chunk $FF exists
    else
        echo "submitting $FF"
        echo "cat ${seqdir}_all_files.txt | /usr/bin/time -v ${metagraph} build -v --parallel ${threads} --mode canonical -k ${K} --mem-cap-gb $((${mem} / 2000)) -o ${outdir}/graph_merged_k${K} --suffix $F 2>&1 | tee ${outdir}/build_k${K}_$F.log" | bsub -J tc_g_${FF} -oo ${outdir}/build_k${K}_${FF}.lsf.log -We 12:00 -n $threads -M $mem -R "rusage[mem=${pmem}]" -R "span[hosts=1]"
    fi
done

