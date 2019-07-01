#!/bin/bash

set -e

K=41

basedir=/cluster/work/grlab/projects/metagenome/data
filelist=${basedir}/gtex/output_k${K}_sequences_per_abundance.files.txt
outdir=${basedir}/gtex/output_k${K}_merged_annotated_per_abundance
mkdir -p $outdir
graph=${basedir}/gtex/output_k${K}_merged/all_merged_k${K}.dbg
mem=15000

chunk=1
if [ ! -z "$1" ]
then
    chunk=$1
fi

chunksize=5000
from=$(( (${chunk} - 1) * ${chunksize} + 1))
to=$(( ${chunk} * ${chunksize}))
files=$(wc -l $filelist | cut -f1 -d ' ')
to=$((${to}<${files}?${to}:${files}))

echo from $from
echo to $to
echo files $files

logfile=build_column_anno_GTEx_per_abundance.k${K}.lsf
echo "fname=\"\$(sed -n \${LSB_JOBINDEX}p ${filelist})\"; \
sample=\$(basename \$fname | cut -f 1-3 -d '.'); \
outbase=${outdir}/all_merged.\${sample}; \
fbase=\$(basename \${fname}); \
fbase=\${fbase%.sequences.fasta.gz}; \
if [ ! -f \${outbase}.column.annodbg ]; then $(pwd)/../../../metagraph/build/metagengraph annotate -i $graph -o \${outbase} --anno-label \${fbase} \${fname}; fi" | bsub -M $mem -n 1 -R "rusage[mem=${mem}]" -We 24:00 -J annoGT${K}[${from}-${to}]%400 -o $logfile
