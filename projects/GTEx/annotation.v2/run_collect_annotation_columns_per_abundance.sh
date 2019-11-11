#!/bin/bash

set -e

mem=750000
mem=900000
mem=1200000

K=41
basedir=/cluster/work/grlab/projects/metagenome/data/gtex
datadir=${basedir}/output_k${K}_annotation_per_abundance_quantiles

metagraph=/cluster/home/akahles/git/software/metagraph_current/metagraph/build/metagraph

for level in 0 #$(seq 0 4)
do
    outbase=${basedir}/output_k${K}_merged.exp_level${level}
    lsffile=${basedir}/output_k${K}_merged.exp_level${level}.column.lsf.log
    logfile=${basedir}/output_k${K}_merged.exp_level${level}.column.log
    if [ ! -f ${outbase}.column.annodbg ]
    then
        echo "/usr/bin/time -v $metagraph merge_anno -o $outbase ${datadir}/*.k${K}.exp_level${level}.sequences.fasta.gz.column.annodbg 2>&1 | tee $logfile" | bsub -M $mem -n 1 -R "rusage[mem=${mem}]" -We 48:00 -J annoGTEx -o $lsffile
    fi
done
