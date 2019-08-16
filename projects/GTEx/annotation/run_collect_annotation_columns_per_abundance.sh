#!/bin/bash

set -e

mem=500000

K=41
basedir=/cluster/work/grlab/projects/metagenome/data/gtex
datadir=${basedir}/output_k${K}_merged_annotated_per_abundance

for level in $(seq 0 4)
do
    outbase=${basedir}/output_k${K}_merged.exp_level${level}
    logfile==${basedir}/output_k${K}_merged.exp_level${level}.column.lsf
    if [ ! -f ${outbase}.column.annodbg ]
    then
        echo "/usr/bin/time -v $(pwd)/../../metagraph/build/metagraph merge_anno -o $outbase ${datadir}/*.k${K}.exp_level${level}.column.annodbg" | bsub -M $mem -n 1 -R "rusage[mem=${mem}]" -We 48:00 -J annoGTEx -o $logfile
    fi
done
