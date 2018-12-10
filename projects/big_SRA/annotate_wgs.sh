#!/usr/bin/bash

home="/cluster/home/mikhaika"

exe="$home/metagengraph_DNA"

dir="$home/big_graph/gut_WGS_k31_30_4"
graph="$dir/graph_gut_WGS_30_4_k31.dbg"
log_dir="$dir/logs"
out_dir="$dir/annotation"


for file in $(cat files_WGS_sorted.txt); do
  #file="/cluster/work/grlab/projects/metagenome/benchmark_kingsford/data_fasta/$x.fasta.gz"
  #ls -lha $file
  if [ -f "$out_dir/$x.color.annodbg" ]; then
    continue
  fi
  x=$(basename ${file%.*.*})
  args="annotate -v --anno-filename -i $graph -o $out_dir/$x --filter-k 20 --filter-abund 30 --filter-thres 4 -p 11 $file"
  job="/usr/bin/time -v $exe $args"
  #$job 2>&1 | tee $log_dir/log_annotate_$x
  bsub -J annotate$x -W 48:00 -n 5 -R "rusage[mem=12000]" -o /dev/null "$job 2>&1 | tee $log_dir/log_annotate_$x"
done
