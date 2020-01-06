#!/bin/bash

set -e

basedir=/cluster/work/grlab/projects/metagenome/data/gtex/
aligndir=${basedir}/alignments_unaligned

mem=1000
threads=1
pmem=1000

for fname in ${aligndir}/*.result.txt.gz
do
    echo "python $(pwd)/analyse_results.query_unaligned2.py $fname" | bsub -M $mem -n $threads -R "rusage[mem=${pmem}]" -We 20:00 -n 1 -J plot_gtex -o /dev/null
done

