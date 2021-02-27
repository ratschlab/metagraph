#!/bin/bash

set -e

K=31

mem=40000
threads=4
pmem=$(($mem / $threads))
basedir=/cluster/work/grlab/projects/metagenome/data
outdir=${basedir}/metagut/graphs/output_k${K}
logdir=${basedir}/metagut/graphs/output_k${K}_logs
mkdir -p $outdir
metagraph=/cluster/home/akahles/git/software/metagraph_current/metagraph/build/metagraph
mkdir -p $logdir

metadata=$(pwd)/../files_WGS_sorted.txt
N=$(wc -l $metadata | cut -f 1 -d ' ')
N=10000

echo "fq1=\"\$(sed -n \${LSB_JOBINDEX}p ${metadata})\"; \:
uuid=\"\$(basename \${fq1} | cut -f 1 -d ',')\"; \
if [ -f ${outdir}/\${uuid}.dbg ]; then exit 0; fi; \
${metagraph} build -v -p $threads -k $K -o ${outdir}/\${uuid} --mode canonical --count-kmers \${fq1}" | bsub -J contigs_k${K}[5000-$N]%400 -o ${logdir}/graph_build%I.lsf.log -We 20:00 -n $threads -M $mem -R "rusage[mem=${pmem}]" -R "span[hosts=1]"
