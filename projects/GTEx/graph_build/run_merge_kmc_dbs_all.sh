#!/bin/bash

set -e

K=41
basedir=/cluster/work/grlab/projects/metagenome/data/gtex
datadir=${basedir}/output_k${K} 
outdir=${basedir}/output_k${K}_merged
mkdir -p $outdir

mem=120000
threads=8
pmem=$(($mem / $threads))

chunk=100
thresh=20

total_lines=$(ls -1 ${datadir}/*.kmc_suf | wc -l)
for FROM in $(seq 1 $chunk $(($total_lines)) )
do
    TO=$((${FROM} + ${chunk}))
    if [ -f ${outdir}/merged.${FROM}_${TO}.t${thresh}.kmc_suf ]
    then
        continue
    fi
    echo "$(pwd)/merge_kmc_dbs.sh $chunk $thresh $FROM" | bsub -M $mem -n $threads -R "rusage[mem=${pmem}]" -We 12:00 -J kmc_merge -o ${outdir}/merged.${FROM}_${TO}.t${thresh}.cluster.log
done
