#!/bin/bash

set -e

basedir=/cluster/work/grlab/projects/metagenome/data/gtex/queries

for fname in ${basedir}/query_unaligned.*.result.txt
do 
    echo processing $fname
    python analyse_results.query_unaligned.py $fname
done
