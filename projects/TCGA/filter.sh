#!/usr/bin/bash

ls -1Sa ~/metagenome/benchmark_TCGA/nobackup/*/*/*.fa.gz > files.txt

bsub -J "filter[1-2030]%500" -W 72:00 -n 9 -R "rusage[mem=3400]" "./filtering.sh \"\$(sed -n \${LSB_JOBINDEX}p files.txt)\" 10 0"
