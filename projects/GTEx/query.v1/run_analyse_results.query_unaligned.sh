#!/bin/bash

set -e

basedir=/cluster/work/grlab/projects/metagenome/data/gtex/queries

source activate python36

for fname in ${basedir}/query_unlabeled.*.result.txt
do 
    echo processing $fname
    python analyse_results.query_unaligned.py $fname
done
