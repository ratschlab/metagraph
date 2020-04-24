#!/bin/bash

set -e

### source paths
# basedir=
# metagraph=
. ../paths.sh

if [ -z "$1" ]
then
    echo "Usage: $0 <K>"
    exit 1
fi
K=$1

mem=30000
graphdir=${basedir}/gtex/graphs/output_k${K}
outdir=${basedir}/gtex/graphs/output_k${K}_sequences
mkdir -p $outdir
filelist=${basedir}/gtex/graphs/gtex_kmc_k${K}.txt

echo  "kmc_file=\"\$(sed -n \${LSB_JOBINDEX}p ${filelist})\"; \
      x=\"\$(echo \$kmc_file | xargs -n 1 basename)\"; \
      $metagraph transform --to-fasta \
      -o ${outdir}/\${x%.kmc_suf}.sequences \
      ${graphdir}/\${x%.kmc_suf}.dbg" | bsub -J contigs_k${K}[1-$(wc -l $filelist | cut -f1 -d ' ')]%400 -o seq_extraction_k${K}.lsf -We 20:00 -n 1 -R "rusage[mem=${mem}]" -R "span[hosts=1]"
