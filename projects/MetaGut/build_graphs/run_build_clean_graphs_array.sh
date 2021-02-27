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

mem=50000
#mem=80000
#mem=150000
mem=350000
threads=4
pmem=$((${mem} / ${threads}))

outdir=${basedir}/metagut/graphs/output_k${K}
logdir=${basedir}/metagut/graphs/output_k${K}_logs
mkdir -p $outdir
mkdir -p $logdir

metadata=$(pwd)/../files_WGS_sorted.txt
N=$(wc -l $metadata | cut -f 1 -d ' ')
N1=4000
N2=20648
N1=3000
N2=3999

echo "fq1=\"\$(sed -n \${LSB_JOBINDEX}p ${metadata})\"; \
      uuid=\"\$(basename \${fq1} | cut -f 1 -d ',')\"; \
      if [ -f ${outdir}/\${uuid}.clean.fasta.gz ]; then exit 0; fi; \
      if [ ! -f ${outdir}/\${uuid}.dbg ]; then ${metagraph} build -v -p $threads -k $K -o ${outdir}/\${uuid} --mode canonical --count-kmers \${fq1}; fi; \
      if [ ! -f ${outdir}/\${uuid}.clean.fasta.gz ]; then ${metagraph} clean --prune-tips $(($K * 2)) --prune-unitigs 0 --fallback 3 -o ${outdir}/\${uuid}.clean ${outdir}/\${uuid}.dbg; fi; \
      rm ${outdir}/\${uuid}.dbg* ${outdir}/\${uuid}.edgemask" | bsub -J contigs_k${K}[${N1}-${N2}]%20 -o ${logdir}/graph_build%I.lsf.log -We 20:00 -n $threads -M $mem -R "rusage[mem=${pmem}]" -R "span[hosts=1]"
