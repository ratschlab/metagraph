#!/usr/bin/env bash

for dir in /cluster/home/mikhaika/metagenome/benchmark_TCGA/nobackup/*; do
  for x in $dir/*.bam; do
    if [ -f $x ] && [ ! -f ${x%.bam}.fa.gz ]; then
      bsub -J "process_bam_$dir" -W 72:00 -n 1 -R "rusage[mem=1000]" "./process_bam.sh $x"
    fi
  done
done
