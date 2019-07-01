#!/bin/bash

set -e

K=41
basedir=/cluster/work/grlab/projects/metagenome/data/gtex
outdir=${basedir}/output_k${K}_merged_eval
mkdir -p $outdir

mem=120000
threads=8
pmem=$(($mem / $threads))

for c in 10 20 30 40 50 100 200 #300 400 500
do
    for t in 3 5 10
    do
        if [ -f ${outdir}/merged.${c}.t${t}.kmc_suf ]
        then
            continue
        fi
        echo "$(pwd)/merge_kmc_dbs.sh $c $t" | bsub -M $mem -n $threads -R "rusage[mem=${pmem}]" -We 12:00 -J kmc_merge -o ${outdir}/merged.${c}.t${t}.cluster.log

    done
done
