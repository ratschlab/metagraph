#!/bin/bash

set -e

### source paths 
# metagraph=
# basedir=
. ../paths.sh

if [ -z "$2" ]
then
    echo "Usage: $0 <K> <tissue>"
    exit 1
fi
K="$1"
tissue=$2

mem=10000
threads=1
pmem=$(($mem / $threads))

graph=${basedir}/tcga/output_k${K}_trimmed_clean_graph_chunked/${tissue}/graph_merged_k${K}.dbg
seqdir=${basedir}/tcga/output_k${K}_trimmed_clean/${tissue}
outdir=${basedir}/tcga/output_k${K}_trimmed_annotation_per_sample/${tissue}
mkdir -p $outdir

for fname in $(find ${seqdir} -name \*.k${K}.clean.0.000000.0.200000.fasta.gz)
do
    sample=$(basename $fname | cut -f 1 -d '.')
    logfile=${outdir}/graph_merged.${sample}.lsf.log
    echo "/usr/bin/time -v $metagraph annotate -v -i ${graph} --outfile-base ${outdir}/graph_merged.${sample} --anno-label ${sample} ${seqdir}/${sample}*.fasta.gz" | bsub -J ms_anno_k${K} -We 72:00 -n $threads -M $mem -R "rusage[mem=${pmem}]" -R "span[hosts=1]" -oo $logfile
done
