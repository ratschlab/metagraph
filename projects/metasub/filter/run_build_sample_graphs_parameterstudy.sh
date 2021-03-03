#!/bin/bash

set -e

K=19

mem=50000
threads=2
pmem=$(($mem / $threads))
basedir=/cluster/work/grlab/projects/metagenome/data
outdir=${basedir}/metasub/graphs_metagraph/output_k${K}
mkdir -p $outdir
metagraph=/cluster/home/akahles/git/software/metagraph_clean/metagraph/build/metagraph
metadata="complete_metadata_extended.clean.v2.pstudy.csv"

while IFS=',' read -r -a line || [[ -n "$line" ]]
do
    uuid=${line[1]}
    fq1=$(echo ${line[40]} | cut -f 1 -d ':')
    fq2=$(echo ${line[40]} | cut -f 2 -d ':')

    outfname=${outdir}/${uuid}.dbg
    logfname=${outdir}/${uuid}.lsf.log
    if [ -f ${outfname} ]
    then
        echo "$uuid complete"
        continue
    fi

    echo "/usr/bin/time -v $metagraph build -p $threads -k $K -o ${outfname%.dbg} --mode canonical --count-kmers $fq1 $fq2" | bsub -J mg_k${K} -o ${logfname} -We 8:00 -n $threads -M $mem -R "rusage[mem=${pmem}]" -R "span[hosts=1]"
done < <(cat ${metadata})
