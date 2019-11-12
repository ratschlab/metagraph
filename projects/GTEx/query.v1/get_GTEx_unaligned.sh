#!/bin/bash

set -e

module load samtools

for fname in /cluster/work/grlab/projects/GTEx/rna/results/alignments/*.bam
do
    outfname=${fname%.bam}.unaligned.sorted.bam
    if [ ! -f ${outfname} ]
    then
        samtools view -h -f4 $fname | samtools sort -n - | samtools view -bS -o ${outfname} -
    fi
done
