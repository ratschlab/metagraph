#!/usr/bin/bash

for F in {\\\$,A,C,G,T,N}; do bsub -J assemble$F -W 8:00 -n 30 -R "rusage[mem=15000]" "ls -1a /cluster/work/grlab/projects/metagenome/benchmark_kingsford/data_fasta/SRR*.gz | /usr/bin/time -v ~/projects2014-metagenome/metagraph/build/metagraph build -v --parallel 30 -k 20 --mem-cap-gb 350 --suffix $F --min-count 3 -o ~/big_graph/graph_SRR_2 2>&1 | tee ~/big_graph/log_assemble_reads_2.txt.30_$F"; done
