#!/bin/bash

set -e

K=41

basedir=/cluster/work/grlab/projects/metagenome/data/gtex
#seqdir=${basedir}/parameterstudy_largest/output_k${K}_cleaned
seqdir=${basedir}/parameterstudy_largest_trimmed/output_k${K}_cleaned

for fname in ${seqdir}/*.clean_s2.fasta.gz
do
    if [ ! -f ${fname%.fasta.gz}.stats ]
    then
        echo "python $(pwd)/gen_fasta_stats.py ${fname}" | bsub -J ms_pst -o /dev/null -We 2:00 -n 1 -M 5000 -R "rusage[mem=5000]" 
    else
        echo "stats for $fname completed"
    fi
done
