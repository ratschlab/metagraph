#!/bin/bash

set -e

mem=30000
K=19
basedir=/cluster/work/grlab/projects/metagenome/data
graphdir=${basedir}/metasub/graphs/output_k${K}
outdir=${basedir}/metasub/contigs/output_k${K}
mkdir -p $outdir
filelist=${graphdir}/metasub_kmc_k${K}.txt

echo  "kmc_file=\"\$(sed -n \${LSB_JOBINDEX}p ${filelist})\"; \
      x=\"\$(echo \$kmc_file | xargs -n 1 basename)\"; \
      ../../metagraph/build/metagraph assemble \
      -o ${outdir}/\${x%.kmc_suf}.sequences \
      ${graphdir}/\${x%.kmc_suf}.dbg" | bsub -J contigs_k${K}[1-4173]%400 -o seq_extraction_${K}.lsf -We 20:00 -n 1 -R "rusage[mem=${mem}]" -R "span[hosts=1]"
