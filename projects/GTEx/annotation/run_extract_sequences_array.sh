#!/bin/bash

set -e

mem=30000
K=41
basedir=/cluster/work/grlab/projects/metagenome/data
graphdir=${basedir}/gtex/output_k${K}
outdir=${basedir}/gtex/output_k${K}_sequences
mkdir -p $outdir
filelist=${basedir}/gtex/gtex_kmc_k${K}.txt

echo  "kmc_file=\"\$(sed -n \${LSB_JOBINDEX}p ${filelist})\"; \
      x=\"\$(echo \$kmc_file | xargs -n 1 basename)\"; \
      ../../../metagraph/build/metagraph transform --to-fasta \
      -o ${outdir}/\${x%.kmc_suf}.sequences \
      ${graphdir}/\${x%.kmc_suf}.dbg" | bsub -J contigs_k${K}[1-$(wc -l $filelist | cut -f1 -d ' ')]%400 -o seq_extraction_k${K}.lsf -We 20:00 -n 1 -R "rusage[mem=${mem}]" -R "span[hosts=1]"
