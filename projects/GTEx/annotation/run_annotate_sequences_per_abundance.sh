#!/bin/bash

set -e

K=41

basedir=/cluster/work/grlab/projects/metagenome/data
seqdir=${basedir}/gtex/output_k${K}_sequences_per_abundance
outdir=${basedir}/gtex/output_k${K}_merged_annotated_per_abundance
mkdir -p $outdir
graph=${basedir}/gtex/output_k${K}_merged/all_merged_k${K}.dbg
mem=15000

for fname in ${seqdir}/*.fasta.gz
do
    sample=$(basename $fname | cut -f 1-3 -d '.')
    outbase=${outdir}/all_merged.${sample}
    logfile=${outbase}.cluster.log
    fbase=$(basename $fname)
    fbase=${fbase%.sequences.fasta.gz}
    if [ ! -f ${outbase}.column.annodbg ]
    then
        echo "/usr/bin/time -v $(pwd)/../../../metagraph/build/metagraph annotate -i $graph -o $outbase --anno-label $fbase $fname" | bsub -M $mem -n 1 -R "rusage[mem=${mem}]" -We 24:00 -J annoGTEx -o $logfile
    else
        echo "$sample already processed"
    fi
done
