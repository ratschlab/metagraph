#!/bin/bash

set -e

if [ -z "$1" ]
then
    echo "Usage: $0 <tissue>"
    exit 1
fi
tissue=$1

basedir=/cluster/work/grlab/projects/metagenome/raw_data/tcga/data/${tissue}/
for fname in $(find $basedir -name \*.bai)
do
    log=${fname%.bai}.lsf.log
    bamname=${fname%.bai}.bam
    echo "$(pwd)/process_bam_fastq.sh $bamname" | bsub -G ms_raets -oo $log -W 48:00 -n 1 -R "rusage[mem=4000]" -J cmprs_bam  
    #echo "$(pwd)/process_bam_fastq.sh $bamname" | bsub -oo $log -W 6:30 -n 1 -R "rusage[mem=4000]" -J cmprs_bam  
done
