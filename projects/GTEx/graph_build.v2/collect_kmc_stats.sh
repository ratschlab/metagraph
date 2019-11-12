#!/bin/bash

set -e

mem=10000

basedir=/cluster/work/grlab/projects/metagenome/data/gtex
K=41

script=/cluster/home/akahles/git/software/metagraph_current3/metagraph/build/run_experiments

for fname in ${basedir}/output_k${K}/*.kmc_suf
do
    fbase=${fname%.kmc_suf}
    outfname=${fbase}.kmc_stats
    if [ ! -f ${outfname} ]
    then
        echo "${script} estimate_abundance_threshold ${fname} > $outfname || rm ${outfname}" | bsub -J kmc_stat -o /dev/null -We 2:00 -n 1 -M $mem -R "rusage[mem=${mem}]" -R "span[hosts=1]"
    fi
done
