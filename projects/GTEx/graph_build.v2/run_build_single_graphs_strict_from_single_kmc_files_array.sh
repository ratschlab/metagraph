#!/bin/bash

set -e

## source paths
# metagraph=
# basedir=
. ../paths.sh

if [ -z "$1" ]
then
    echo "Usage: $0 <K>"
    exit 1
fi
K=$1

#mem=75000
mem=50000
threads=4

countdir=${basedir}/gtex/output_k${K}
outdir=${basedir}/gtex/output_k${K}_clean_strict
mkdir -p $outdir
filelist=${basedir}/gtex/gtex_kmc_k${K}.txt
if [ ! -f ${filelist} ]
then
    ls -1 ${countdir}/*.kmc_suf > $filelist
fi

echo "kmc_file=\"\$(sed -n \${LSB_JOBINDEX}p ${filelist})\"; \
      x=\"\$(basename \$kmc_file)\"; \
      if [ -f ${outdir}/\${x%.kmc_suf}.dbg ]; then exit 0; fi; \
      outbase=\"${outdir}/\${x%.kmc_suf}\"; \
      ${metagraph} build -v -p $threads -k $K -o \${outbase} --canonical --count-kmers \${kmc_file}; \
      ${metagraph} clean -v --prune-tips $(($K * 3)) --prune-unitigs 0 --fallback 2 -o \${outbase}.clean \${outbase}.dbg; \
      rm \${outbase}.dbg* \${outbase}.edgemask" #| bsub -J contigs_k${K}[2-$(wc -l $filelist | cut -f1 -d ' ')]%400 -o graph_clean_gtex_k${K}.lsf -We 20:00 -n $threads -R "rusage[mem=${mem}]" -R "span[hosts=1]"
