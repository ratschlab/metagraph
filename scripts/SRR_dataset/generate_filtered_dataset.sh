#!/usr/bin/bash

home="/cluster/home/mikhaika"

exe="$home/metagengraph"


for file in $(ls /cluster/work/grlab/projects/metagenome/benchmark_kingsford/data_fasta/SRR*.fasta.gz); do
  #file="/cluster/work/grlab/projects/metagenome/benchmark_kingsford/data_fasta/$x.fasta.gz"
  #ls -lha $file
  x=$(basename ${file%.*.*})
  args="filter -v -k 20 --generate-fastq --filter-abund 3 $file"
  job="/usr/bin/time -v $exe $args"
  #$job 2>&1 | tee $log_dir/log_annotate_$x
  bsub -J annotate$x -W 12:00 -n 1 -R "rusage[mem=1000]" "$job"
done
