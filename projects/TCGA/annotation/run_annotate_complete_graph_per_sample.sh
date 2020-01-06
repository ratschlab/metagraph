#!/bin/bash

set -e

### source paths 
# metagraph=
# basedir=
. ../paths.sh

if [ -z "$1" ]
then
    echo "Usage: $0 <K>"
    exit 1
fi
K="$1"

mem=10000
threads=1
pmem=$(($mem / $threads))

graph=${basedir}/tcga/graph_merged_complete_k${K}.dbg
seqdir=${basedir}/tcga/output_k${K}_trimmed_clean
outdir=${basedir}/tcga/graph_merged_complete_k${K}_annotation_per_sample
mkdir -p $outdir

for fname in $(find ${seqdir} -name \*.k${K}.clean.0.000000.0.200000.fasta.gz)
do
    sample=$(basename $fname | cut -f 1 -d '.')
    logfile=${outdir}/graph_merged.${sample}.lsf.log
    files=$(ls -1 ${fname%.0.000000.0.200000.fasta.gz}*.fasta.gz | tr '\n' ' ')
    echo "/usr/bin/time -v $metagraph annotate -v -i ${graph} --outfile-base ${outdir}/graph_merged.${sample} --anno-label ${sample} $files" | bsub -J ms_anno_k${K} -We 72:00 -n $threads -M $mem -R "rusage[mem=${pmem}]" -R "span[hosts=1]" -oo $logfile
done
