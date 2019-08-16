#!/usr/bin/bash

bsub -J "trim[1-20000]%500" -W 72:00 -n 5 -R "rusage[mem=2400]" "../filtering.sh /cluster/work/grlab/projects/metagenome/benchmark_human_metagenome/nobackup/human_gut_sra/\"\$(sed -n \${LSB_JOBINDEX}p files_sorted_desc_1.txt | sed 's|fastq.gz|trimfq_0.02.fastq.gz|g')\" 8 4"

