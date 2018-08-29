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
  args="annotate -v --anno-filename -i $graph -o $out_dir/$x --filter-abund 1 -p 10 $file"
  job="/usr/bin/time -v $exe $args"
  #$job 2>&1 | tee $log_dir/log_annotate_$x
  bsub -J annotate$x -W 12:00 -n 9 -R "rusage[mem=2500]" "$job 2>&1 | tee $log_dir/log_annotate_$x"
done



# bsub -J "annotate_kingsford[1-2586]%600" -W 96:00 -n 6 -R "rusage[mem=3000]" -o /dev/null "file=\"\$(sed -n \${LSB_JOBINDEX}p ~/projects2014-metagenome/scripts/kingsford/kingsford_filtered_list.txt)\"; x=\$(basename \${file%.*.*}); gtime -v ~/projects2014-metagenome/metagraph/build/metagengraph_DNA annotate -v --anno-filename -i ~/big_graph/kingsford/graph_k19.dbg -o ~/big_graph/kingsford/annotation/\$x -p 13 \$file 2>&1 | tee ~/big_graph/kingsford/logs/log_annotate_\$x"
