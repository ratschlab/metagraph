#!/bin/bash

set -e 

if [ -z "$1" ]
then
    echo "Usage: $0 <tissue>"
    exit 1
fi
tissue=$1

K=31
#mem=32000
mem=75000
#mem=100000
threads=2
pmem=$(($mem / $threads))

fastqdir=/cluster/work/grlab/projects/metagenome/raw_data/tcga/data/${tissue}
outdir=/cluster/work/grlab/projects/metagenome/data/tcga/output_k${K}_trimmed/${tissue}
mkdir -p $outdir
metadata=/cluster/work/grlab/projects/metagenome/raw_data/tcga/metadata/gdc_manifest.2020-02-22.${tissue}.txt
#fnames_dir=${outdir}/fnames_kmc
#mkdir -p $fnames_dir

while IFS=',' read -r -a line || [[ -n "$line" ]]
do
    uuid=${line[0]}
    fname=${line[1]}
    if [ "$uuid" == "id" ]
    then 
        continue
    fi  

    if [ -f "${outdir}/${uuid}.k${K}.kmc_suf" ]
    then
        echo ${uuid} already done
        continue
    fi
    if [ ! -d ${fastqdir}/${uuid} ]
    then
        continue
    fi
    fq1=$(find ${fastqdir}/${uuid} -name \*sorted.r1.fq.gz)
    fq2=$(find ${fastqdir}/${uuid} -name \*sorted.r2.fq.gz)
    if [ -z "${fq1}" -a -z "${fq2}" ]
    then
        echo No fastq file for $uuid
        continue
    fi
    echo "/usr/bin/time -v $(pwd)/count_kmers.sh $K $uuid $outdir $fq1 $fq2" | bsub -G ms_raets -M $mem -n $threads -R "rusage[mem=${pmem}]" -We 12:00 -J kmc_tcga -oo ${outdir}/${uuid}.k${K}.lsf.log
done < <(cat $metadata | tr $'\t' ',')
