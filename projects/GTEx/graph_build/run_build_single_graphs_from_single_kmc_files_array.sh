#!/bin/bash

set -e

#mem=75000
mem=50000
K=41
threads=4
basedir=/cluster/work/grlab/projects/metagenome/data
outdir=${basedir}/gtex/output_k${K}
mkdir -p $outdir
filelist=${basedir}/gtex/gtex_kmc_k${K}.txt
ls -1 ${outdir}/*.kmc_suf > $filelist

echo "kmc_file=\"\$(sed -n \${LSB_JOBINDEX}p ${filelist})\"; \
      x=\"\$(echo \$kmc_file | xargs -n 1 basename)\"; \
      if [ -f ${outdir}/\${x%.kmc_suf}.dbg ]; then exit 0; fi; \
      ../../../metagraph/build/metagengraph build -p $threads -k $K -o ${outdir}/\${x%.kmc_suf} --kmc \${kmc_file}" | bsub -J contigs_k${K}[1-$(wc -l $filelist | cut -f1 -d ' ')]%400 -o graph_build_gtex_k${K}.lsf -We 20:00 -n $threads -R "rusage[mem=${mem}]" -R "span[hosts=1]"
