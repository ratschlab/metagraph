#!/bin/bash

set -e

#mem=75000
mem=25000
K=41
threads=4
basedir=/cluster/work/grlab/projects/metagenome/data
outdir=${basedir}/gtex/output_k${K}_merged_eval
mkdir -p $outdir

for c in 10 20 30 40 50 100 200 #300 400 500
do
    for t in 3 5 10
    do  
        if [ -f ${outdir}/merged.${c}.t${t}.dbg ]
        then
            continue
        fi  
        echo "$(pwd)/../../metagraph/build/metagraph build -p $threads -k $K -o ${outdir}/merged.${c}.t${t} --kmc ${outdir}/merged.${c}.t${t}.kmc_suf" | bsub -J kmc2dbg -o ${outdir}/merged.${c}.t${t}.kmc2dbg.log -We 20:00 -n $threads -R "rusage[mem=${mem}]" -R "span[hosts=1]"
    done
done
