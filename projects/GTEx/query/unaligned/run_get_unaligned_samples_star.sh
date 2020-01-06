#!/bin/bash

set -e

basedir=/cluster/work/grlab/projects/metagenome/data/gtex/
aligndir=${basedir}/alignments_star_small/results/alignments

mem=1000
threads=1
pmem=1000

for fname in ${aligndir}/*.bam
do
    outfname=${fname%.bam}.unaligned.txt

    if [ ! -f ${outfname} ]
    then
        echo "module load samtools; samtools view -F 256 ${fname} | python $(pwd)/get_unaligned_samples_star.py > $outfname" | bsub -M $mem -n $threads -R "rusage[mem=${pmem}]" -We 20:00 -n 1 -J sum_gtex -o /dev/null
    fi
done

