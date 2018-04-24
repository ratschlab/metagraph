#!/usr/bin/bash

home="/cluster/home/mikhaika"

exe="$home/metagengraph"

log_dir="$home/big_graph/graph_SRR_k20_1"
graph="$home/big_graph/graph_SRR_k20_1/graph_SRR_1.dbg"
out_dir="$home/big_graph/graph_SRR_k20_1/annotation"


for file in $(ls /cluster/work/grlab/projects/metagenome/benchmark_kingsford/data_fasta/SRR*.gz); do
  #file="/cluster/work/grlab/projects/metagenome/benchmark_kingsford/data_fasta/$x.fasta.gz"
  #ls -lha $file
  x=$(basename ${file%.*.*})
  args="annotate -v --anno-filename -i $graph -o $out_dir/$x --noise-freq 1 -p 10 $file"
  job="/usr/bin/time -v $exe $args"
  #$job 2>&1 | tee $log_dir/log_annotate_$x
  bsub -J annotate$x -W 12:00 -n 9 -R "rusage[mem=2500]" "$job 2>&1 | tee $log_dir/log_annotate_$x"
done
