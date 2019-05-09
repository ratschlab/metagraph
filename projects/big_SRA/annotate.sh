#!/usr/bin/bash

home="/cluster/home/mikhaika"

exe="$home/metagengraph_DNA"

log_dir="$home/big_graph/gut_amplicon_k31_8_4"
graph="$home/big_graph/gut_amplicon_k31_8_4/graph_gut_AMPLICON_8_4_k31.dbg"
out_dir="$home/big_graph/gut_amplicon_k31_8_4/annotation"


for file in $(cat files_AMPLICON.txt); do
  #file="/cluster/work/grlab/projects/metagenome/benchmark_kingsford/data_fasta/$x.fasta.gz"
  #ls -lha $file
  x=$(basename ${file%.*.*})
  args="annotate -v --anno-filename -i $graph -o $out_dir/$x --filter-k 20 --min-count 9 --filter-thres 4 -p 11 $file"
  job="/usr/bin/time -v $exe $args"
  #$job 2>&1 | tee $log_dir/log_annotate_$x
  bsub -J annotate$x -W 12:00 -n 5 -R "rusage[mem=2500]" -o /dev/null "$job 2>&1 | tee $log_dir/log_annotate_$x"
done
