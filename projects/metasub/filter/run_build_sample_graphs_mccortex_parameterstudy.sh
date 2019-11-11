#!/bin/bash

set -e

K=19

mem=50000
threads=2
pmem=$(($mem / $threads))
basedir=/cluster/work/grlab/projects/metagenome/data
outdir=${basedir}/metasub/graphs_mccortex/output_k${K}
mkdir -p $outdir
mccortex=/cluster/home/akahles/git/software/mccortex/bin/mccortex31
metadata="complete_metadata_extended.clean.v2.pstudy.csv"

while IFS=',' read -r -a line || [[ -n "$line" ]]
do
    uuid=${line[1]}
    #fq1=$(echo ${line[40]} | cut -f 1 -d ':')
    #fq2=$(echo ${line[40]} | cut -f 2 -d ':')
    fqs=${line[40]}

    outfname=${outdir}/${uuid}.ctx
    logfname=${outdir}/${uuid}.lsf.log
    if [ -f ${outfname} ]
    then
        echo "$uuid complete"
        continue
    fi

    echo "/usr/bin/time -v $mccortex build -m $((${mem} / 2000))G -k $K -t $threads --sample $uuid --seq2 $fqs ${outfname}" | bsub -J mc_k${K} -o ${logfname} -We 8:00 -n $threads -M $mem -R "rusage[mem=${pmem}]" -R "span[hosts=1]"
done < <(cat ${metadata})
