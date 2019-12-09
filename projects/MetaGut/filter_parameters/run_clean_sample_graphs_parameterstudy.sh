#!/bin/bash

set -e

K=31

mem=15000
#mem=50000
#mem=100000
threads=1
pmem=$(($mem / $threads))
basedir=/cluster/work/grlab/projects/metagenome/data
graphdir=${basedir}/metagut/graphs/output_k${K}
outdir=${basedir}/metagut/parameterstudy/output_k${K}_cleaned
mkdir -p $outdir
metagraph=/cluster/home/akahles/git/software/metagraph_current/metagraph/build/metagraph

samples="samples_parameterstudy.tsv"
if [ ! -f "$samples" ]
then
    ls -1 ${graphdir}/*.dbg > $samples
fi
index="10 45 167 288 496 1022 1340 1670 1988 2130"

tipfactor="1 2 3 4 5 6"
fallback="2 4 6 8 10 12"

for i in $index
do
    graph=$(sed -n ${i}p $samples)
    uuid=$(basename $graph | cut -f 1 -d '.')
    for tf in $tipfactor
    do
        for fb in $fallback
        do
            outfname=${outdir}/${uuid}.tf_${tf}.fb_${fb}.clean_s1.fasta.gz
            logfname=${outdir}/${uuid}.tf_${tf}.fb_${fb}.clean_s1.lsf.log
            if [ -f ${outfname} ]
            then
                echo "$uuid complete"
                continue
            fi

            echo "/usr/bin/time -v $metagraph clean --prune-tips $(($K * $tf)) --prune-unitigs 0 --fallback $fb -o ${outfname%.fasta.gz} -v ${graph}" | bsub -J mg_clean -o ${logfname} -We 8:00 -n 1 -M $mem -R "rusage[mem=${pmem}]" -R "span[hosts=1]"
        done
    done
done
