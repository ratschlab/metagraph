#!/bin/bash

set -e

K=41

mem=10000
mem=15000
threads=1
pmem=$(($mem / $threads))

fastqdir=/cluster/work/grlab/SHARED_DATA_PROPOSAL_OTHERWISE_DELETE/GTEx/extract/dbGaP-9608/fastq
#fastqdir=/cluster/work/grlab/projects/metagenome/data/gtex/parameterstudy_largest_trimmed/fastq_trimmed
basedir=/cluster/work/grlab/projects/metagenome/data
#graphdir=${basedir}/gtex/parameterstudy/output_k${K}_cleaned
#outdir=${basedir}/gtex/parameterstudy/output_k${K}_cleaned_realign
graphdir=${basedir}/gtex/parameterstudy_largest_trimmed/output_k${K}_cleaned
outdir=${basedir}/gtex/parameterstudy_largest_trimmed/output_k${K}_cleaned_realign
mkdir -p $outdir
metagraph=/cluster/home/akahles/git/software/metagraph_current2/metagraph/build/metagraph
#metadata="samples_parameterstudy.tsv"
metadata="samples_parameterstudy.largest.tsv"

while IFS=',' read -r -a line || [[ -n "$line" ]]
do
    uuid=${line[16]}
    fq1=${fastqdir}/${uuid}_1.fastq.gz
    fq2=${fastqdir}/${uuid}_2.fastq.gz
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

        echo "/usr/bin/time -v $metagraph align -i ${graph} --query-presence $fq1 $fq2 > ${outtmp} && grep -w 0 ${outtmp} | wc -l > ${outfname} && grep -w 1 ${outtmp} | wc -l >> ${outfname} && rm ${outtmp}" | bsub -J mg_k${K} -o ${logfname} -We 8:00 -n $threads -M $mem -R "rusage[mem=${pmem}]" -R "span[hosts=1]"
    done
done < <(cat ${metadata} | tr $'\t' ',')
