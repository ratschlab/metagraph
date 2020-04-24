#!/bin/bash

set -e

for fname in /cluster/work/grlab/projects/metagenome/data/gtex/alignments_star_small/results/alignments/*.bam
do
    echo processing $fname
    outfname=${fname%.bam}.unaligned.sorted.bam
    if [ ! -f ${outfname} ]
    then
        echo "module load samtools; samtools view -h -f4 $fname | samtools sort -n - | samtools view -bS -o ${outfname} -" | bsub -M 4000 -We 8:00 -n 1 -R "rusage[mem=4000]" -R "span[hosts=1]" -J b2fq -o /dev/null

    fi
done
