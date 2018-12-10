#!/usr/bin/bash

home="/cluster/home/mikhaika"

exe="$home/metagengraph_DNA"

dir="$home/big_graph/metasub_k31_1_4"
graph="$dir/graph_metasub_1_4_k31.dbg"
log_dir="$dir/logs"
out_dir="$dir/annotation"


for file in $(cat metasub_cured.txt | tail -n +11); do
  #file="/cluster/work/grlab/projects/metagenome/benchmark_kingsford/data_fasta/$x.fasta.gz"
  #ls -lha $file
  if [ -f "$out_dir/$x.color.annodbg" ]; then
    continue
  fi
  x=$(basename ${file%.*.*})
  args="annotate -v --anno-filename -i $graph -o $out_dir/$x --filter-k 20 --filter-abund 1 --filter-thres 4 -p 25 $file"
  job="/usr/bin/time -v $exe $args"
  #$job 2>&1 | tee $log_dir/log_annotate_$x
  bsub -J annotate$x -W 48:00 -n 12 -R "rusage[mem=12000]" -o /dev/null "$job 2>&1 | tee $log_dir/log_annotate_$x"
done


# bsub -J "annotate_metasub[1-2809]%600" -W 96:00 -n 15 -R "rusage[mem=12000]" -o /dev/null "file=\"\$(sed -n \${LSB_JOBINDEX}p metasub_cured.txt)\"; x=\$(basename \${file%.*.*}); gtime -v ~/projects2014-metagenome/metagraph/build/metagengraph_DNA annotate -v --anno-filename -i ~/big_graph/metasub_k31_1_4/graph_metasub_1_4_k31.dbg -o ~/big_graph/metasub_k31_1_4/annotation/\$x --filter-k 20 --filter-abund 1 --filter-thres 4 -p 30 \$file 2>&1 | tee ~/big_graph/metasub_k31_1_4/logs/log_annotate_\$x"
