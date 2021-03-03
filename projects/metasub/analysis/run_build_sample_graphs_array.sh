#!/bin/bash

set -e

K=41

mem=40000
mem=700000
threads=4
pmem=$(($mem / $threads))
basedir=/cluster/work/grlab/projects/metagenome/data
outdir=${basedir}/metasub/graphs/output_k${K}
mkdir -p $outdir
metagraph=/cluster/home/akahles/git/software/metagraph_clean/metagraph/build/metagraph

N=4236

metadata=$(pwd)/../complete_metadata_extended.clean.v2.csv
echo "fq1=\"\$(sed -n \${LSB_JOBINDEX}p ${metadata} | cut -f 41 -d ',' | cut -f 1 -d ':')\"; \
fq2=\"\$(sed -n \${LSB_JOBINDEX}p ${metadata} | cut -f 41 -d ',' | cut -f 2 -d ':')\"; \
uuid=\"\$(sed -n \${LSB_JOBINDEX}p ${metadata} | cut -f 2 -d ',')\"; \
if [ -f ${outdir}/\${uuid}.dbg ]; then exit 0; fi; \
${metagraph} build -p $threads -k $K -o ${outdir}/\${uuid} --mode canonical --count-kmers \${fq1} \${fq2}" | bsub -J contigs_k${K}[2-$(($N + 1))]%400 -o ${outdir}.graph_build.lsf -We 20:00 -n $threads -M $mem -R "rusage[mem=${pmem}]" -R "span[hosts=1]"
