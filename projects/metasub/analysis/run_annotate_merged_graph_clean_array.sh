#!/bin/bash

set -e

if [ -z "$1" ]
then
    echo "Usage: $0 <K>"
    exit 1
fi
K=$1

mem=150000
threads=6
pmem=$(($mem / $threads))
basedir=/cluster/work/grlab/projects/metagenome/data
graphdir=${basedir}/metasub/graphs/output_k${K}_cleaned_graph_chunked
outdir=${basedir}/metasub/graphs/output_k${K}_cleaned_graph_annotation
seqdir=${basedir}/metasub/graphs/output_k${K}_cleaned
mkdir -p $outdir
metagraph=/cluster/home/akahles/git/software/metagraph_current_dev/metagraph/build/metagraph

N=4236

metadata=$(pwd)/../complete_metadata_extended.clean.v2.csv
echo "uuid=\"\$(sed -n \${LSB_JOBINDEX}p ${metadata} | cut -f 2 -d ',')\"; \
if [ -f ${outdir}/\${uuid}.column.annodbg ]; then exit 0; fi; \
${metagraph} annotate -p $threads -k $K -i ${graphdir}/graph_merged_k${K}.dbg -o ${outdir}/\${uuid} --fwd-and-reverse --anno-label ${seqdir}/\${uuid}.clean.fasta.gz" | bsub -J ms_anno_k${K}[2-$(($N + 1))]%400 -o ${outdir}.annotate.lsf -We 20:00 -n $threads -M $mem -R "rusage[mem=${pmem}]" -R "span[hosts=1]"
