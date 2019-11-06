#!/bin/bash

set -e

K=19

#mem=10000
mem=50000
#mem=100000
threads=1
pmem=$(($mem / $threads))
basedir=/cluster/work/grlab/projects/metagenome/data
graphdir=${basedir}/metasub/graphs_metagraph/output_k${K}
outdir=${basedir}/metasub/graphs_metagraph/output_k${K}_cleaned_pstudy
mkdir -p $outdir
metagraph=/cluster/home/akahles/git/software/metagraph_current/metagraph/build/metagraph

tipfactor="1 2 3 4"
fallback="2 3 4 5 6 8 10 12"

for graph in $(ls -1S ${graphdir}/*.dbg | head -n 1000 | tail -n 10)
do
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

            echo "/usr/bin/time -v $metagraph clean --prune-tips $(($K * $tf)) --prune-unitigs 0 --fallback $fb -o ${outfname%.fasta.gz} -v ${graphdir}/${uuid}.dbg" | bsub -J ms_clean -o ${logfname} -We 8:00 -n 1 -M $mem -R "rusage[mem=${pmem}]" -R "span[hosts=1]"
        done
    done
done
