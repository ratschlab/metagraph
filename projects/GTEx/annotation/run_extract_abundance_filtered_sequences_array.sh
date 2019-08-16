#!/bin/bash

set -e

mem=150000
K=41
basedir=/cluster/work/grlab/projects/metagenome/data
outdir=${basedir}/gtex/output_k${K}_sequences_per_abundance
mkdir -p $outdir
filelist=${basedir}/gtex/gtex_kmc_k${K}.txt

echo  "kmc_file=\"\$(sed -n \${LSB_JOBINDEX}p ${filelist})\"; \
$(pwd)/extract_abundance_filtered_sequence.sh \
\${kmc_file} ${outdir} $K" | bsub -J contigs_k${K}[1-$(wc -l $filelist | cut -f1 -d ' ')]%400 -o seq_ext_abundances_k${K}.lsf -We 20:00 -n 1 -R "rusage[mem=${mem}]" -R "span[hosts=1]"
