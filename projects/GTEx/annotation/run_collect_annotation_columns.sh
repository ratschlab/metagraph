#!/bin/bash

set -e

mem=500000

K=41
basedir=/cluster/work/grlab/projects/metagenome/data/gtex
datadir=${basedir}/output_k${K}_merged_annotated

outbase=${basedir}/output_k${K}_merged.anno
logfile==${basedir}/output_k${K}_merged.column.lsf
if [ ! -f ${outbase}.column.annodbg ]
then
    echo "/usr/bin/time -v $(pwd)/../../../metagraph/build/metagengraph merge_anno -o $outbase ${datadir}/all_merged_k${K}.*.column.annodbg" | bsub -M $mem -n 1 -R "rusage[mem=${mem}]" -We 48:00 -J annoGTEx -o $logfile
fi
