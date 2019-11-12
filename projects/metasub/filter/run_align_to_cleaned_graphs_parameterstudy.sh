#!/bin/bash

set -e

K=19

mem=10000
threads=1
pmem=$(($mem / $threads))
basedir=/cluster/work/grlab/projects/metagenome/data
graphdir=${basedir}/metasub/graphs_metagraph/output_k19_cleaned_pstudy
outdir=${basedir}/metasub/graphs_metagraph/output_k19_cleaned_pstudy_realign
mkdir -p $outdir
metagraph=/cluster/home/akahles/git/software/metagraph_clean/metagraph/build/metagraph
metadata="complete_metadata_extended.clean.v2.pstudy.csv"

while IFS=',' read -r -a line || [[ -n "$line" ]]
do
    uuid=${line[1]}
    fq1=$(echo ${line[40]} | cut -f 1 -d ':')
    fq2=$(echo ${line[40]} | cut -f 2 -d ':')
    for graph in ${graphdir}/${uuid}*.clean_s1.dbg
    do
        graphbase=$(basename $graph)
        outtmp=${outdir}/${graphbase%.dbg}.align.result
        outfname=${outdir}/${graphbase%.dbg}.align.stats
        logfname=${outdir}/${graphbase%dbg}.align.lsf.log
        if [ -f ${outfname} ]
        then
            echo "$uuid $graph complete"
            continue
        fi

        echo "/usr/bin/time -v $metagraph align -i ${graph} --query-presence -r $fq1 $fq2 > ${outtmp} && grep -w 0 ${outtmp} | wc -l > ${outfname} && grep -w 1 ${outtmp} | wc -l >> ${outfname} && rm ${outtmp}" | bsub -J mg_k${K} -o ${logfname} -We 8:00 -n $threads -M $mem -R "rusage[mem=${pmem}]" -R "span[hosts=1]"
    done
done < <(cat ${metadata})
