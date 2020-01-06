#!/bin/bash

set -e

if [ -z "$1" ]
then
    echo "Usage: $0 <K>"
    exit 1
fi
K=$1

mem=100000
threads=1
pmem=$(($mem / $threads))
basedir=/cluster/work/grlab/projects/metagenome/data
graphdir=${basedir}/metasub/graphs/output_k${K}
outdir=${basedir}/metasub/graphs/output_k${K}_cleaned
mkdir -p $outdir
metagraph=/cluster/home/akahles/git/software/metagraph_current_dev/metagraph/build/metagraph

metadata=$(pwd)/../complete_metadata_extended.clean.v2.csv

while IFS=',' read -r -a line || [[ -n "$line" ]]
do
    uuid=${line[1]}
    if [ "$uuid" == "uuid" ]
    then
        continue
    fi
    if [ ! -f ${graphdir}/${uuid}.dbg ]
    then
        echo "graph incomplete from $uuid"
        continue
    fi
    if [ -f ${outdir}/${uuid}.clean.fasta.gz ]
    then
        echo "$uuid complete"
        continue
    fi

    echo "/usr/bin/time -v $metagraph clean --prune-tips $(($K * 2)) --prune-unitigs 0 --fallback 3 -o ${outdir}/${uuid}.clean -v ${graphdir}/${uuid}.dbg" | bsub -J ms_cl_k${K} -o ${outdir}/${uuid}.clean.lsf.log -We 48:00 -n 1 -M $mem -R "rusage[mem=${pmem}]" -R "span[hosts=1]"
done < <(cat $metadata)
