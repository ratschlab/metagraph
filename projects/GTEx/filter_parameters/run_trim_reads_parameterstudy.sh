#!/bin/bash

set -e

mem=10000
mem=15000
threads=1
pmem=$(($mem / $threads))

fastqdir=/cluster/work/grlab/SHARED_DATA_PROPOSAL_OTHERWISE_DELETE/GTEx/extract/dbGaP-9608/fastq
basedir=/cluster/work/grlab/projects/metagenome/data
outdir=${basedir}/gtex/parameterstudy_largest_trimmed/fastq_trimmed
mkdir -p $outdir
#metadata="samples_parameterstudy.tsv"
metadata="samples_parameterstudy.largest.tsv"
q=0.05

while IFS=',' read -r -a line || [[ -n "$line" ]]
do
    uuid=${line[16]}
    fq1=${fastqdir}/${uuid}_1.fastq.gz
    fq2=${fastqdir}/${uuid}_2.fastq.gz
    fq1_out=${outdir}/${uuid}_1.trim_${q}.fastq.gz
    fq2_out=${outdir}/${uuid}_2.trim_${q}.fastq.gz

    if [ ! -f ${fq1_out} ]
    then
        logfname=${fq1_out%.fastq.gz}.lsf.log
        echo "module load seqtk; seqtk trimfq -q $q $fq1 | gzip > ${fq1_out}" | bsub -J trim -o $logfname -We 8:00 -n $threads -M $mem -R "rusage[mem=${pmem}]" -R "span[hosts=1]"
    fi
    if [ ! -f ${fq2_out} ]
    then
        logfname=${fq2_out%.fastq.gz}.lsf.log
        echo "module load seqtk; seqtk trimfq -q $q $fq2 | gzip > ${fq2_out}" | bsub -J trim -o $logfname -We 8:00 -n $threads -M $mem -R "rusage[mem=${pmem}]" -R "span[hosts=1]"
    fi
done < <(cat ${metadata} | tr $'\t' ',')
