#!/bin/bash

set -e

K=31

mem=10000
threads=1
pmem=$(($mem / $threads))

basedir=/cluster/work/grlab/projects/metagenome/data
fastqdir=${basedir}/MetaGut/nobackup/human_gut_sra
graphdir=${basedir}/metagut/parameterstudy/output_k${K}_cleaned
outdir=${basedir}/metagut/parameterstudy/output_k${K}_cleaned_realign
mkdir -p $outdir
metagraph=/cluster/home/akahles/git/software/metagraph_current/metagraph/build/metagraph
metadata="samples_parameterstudy.tsv"

index="10 45 167 288 496 1022 1340 1670 1988 2130"

for i in $index
do
    graph=$(sed -n ${i}p $metadata)
    uuid=$(basename $graph | cut -f 1 -d '.')
    fq=$(basename $graph)
    fq=${fastqdir}/${fq%.dbg}
    for graph in ${graphdir}/${uuid}*.clean_s1.dbg
    do
        graphbase=$(basename $graph)
        outtmp=${outdir}/${graphbase%.dbg}.align.result
        outfname=${outdir}/${graphbase%.dbg}.align.stats
        logfname=${outdir}/${graphbase%.dbg}.align.lsf.log
        if [ -f ${outfname} ]
        then
            echo "$uuid $graph complete"
            continue
        fi

        echo "/usr/bin/time -v $metagraph align -i ${graph} --query-presence $fq > ${outtmp} && grep -w 0 ${outtmp} | wc -l > ${outfname} && grep -w 1 ${outtmp} | wc -l >> ${outfname} && rm ${outtmp}" | bsub -J mg_k${K} -o ${logfname} -We 8:00 -n $threads -M $mem -R "rusage[mem=${pmem}]" -R "span[hosts=1]"
    done
done
