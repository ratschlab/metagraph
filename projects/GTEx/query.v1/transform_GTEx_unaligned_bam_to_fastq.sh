#!/bin/bash

set -e

samtools=/cluster/work/grlab/share/modules/packages/samtools/1.6/bin/samtools
bam2fastq=/cluster/home/akahles/git/projects/2013/PanCancerTCGA/rerun2017/sample_processing_rna/bam2fastq_single.py

datadir=/cluster/work/grlab/projects/GTEx/rna/results/alignments

for fname in $(find ${datadir} -name \*unaligned.sorted.bam)
do
    if [ ! -f ${fname%bam}r1.fq.gz ]
    then
        echo "python $bam2fastq ${fname} && gzip ${fname%bam}r1.fq" | bsub -M 4000 -We 8:00 -n 3 -R "rusage[mem=4000]" -R "span[hosts=1]" -J b2fq -o /dev/null
    else
        echo $fname done
    fi
done

