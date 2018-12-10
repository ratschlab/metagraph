#!/usr/bin/bash

DATA_DIR="/cluster/work/grlab/projects/metagenome/benchmark_kingsford/data_fasta"

ls -1S "$DATA_DIR"/*.fasta.gz > kingsford_list.txt
num_files=$(cat kingsford_list.txt | wc -l)

bsub -J "filter[1-$num_files]%500" -W 72:00 -n 10 -R "rusage[mem=3400]" "./filtering.sh \"\$(sed -n \${LSB_JOBINDEX}p kingsford_list.txt)\""
