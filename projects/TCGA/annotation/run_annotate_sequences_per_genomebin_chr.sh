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

mem=100000
threads=20
pmem=$(($mem / $threads))

graph=${basedir}/tcga/output_k${K}_trimmed_clean_graph_chunked/${tissue}/graph_merged_k${K}.dbg
outdir=${basedir}/tcga/output_k${K}_trimmed_clean.${tissue}.genomebins_per_chr
mkdir -p $outdir
for seq in ${basedir}/tcga/annotation/hg38_tiles_per_chr/*.fa
do
    seqbase=$(basename $seq)
    outbase=${outdir}/${seqbase%.fa}
    logfile=${outdir}/${seqbase%.fa}.lsf.log
    echo "/usr/bin/time -v $metagraph annotate -v -i ${graph} --outfile-base ${outbase} --parallel $threads --anno-header ${seq}" | bsub -J tc_anno_k${K} -We 72:00 -n $threads -M $mem -R "rusage[mem=${pmem}]" -R "span[hosts=1]" -oo $logfile
done
