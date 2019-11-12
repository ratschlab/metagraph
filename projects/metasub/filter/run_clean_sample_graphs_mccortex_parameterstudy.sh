#!/bin/bash

set -e

K=19

#mem=10000
mem=80000
#mem=100000
threads=2
pmem=$(($mem / $threads))
basedir=/cluster/work/grlab/projects/metagenome/data
graphdir=${basedir}/metasub/graphs_mccortex/output_k${K}
outdir=${basedir}/metasub/graphs_mccortex/output_k${K}_cleaned_pstudy
mkdir -p $outdir
mccortex=/cluster/home/akahles/git/software/mccortex/bin/mccortex31

tipfactor="1 2 3 4"
fallback="2 4 5 6 8 10 12"

for graph in ${graphdir}/*.ctx
do
    uuid=$(basename $graph | cut -f 1 -d '.')
    for tf in $tipfactor
    do
        for fb in $fallback
        do
            outgraph=${outdir}/${uuid}.tf_${tf}.fb_${fb}.clean.ctx
            outseq=${outdir}/${uuid}.tf_${tf}.fb_${fb}.clean.fasta.gz
            logfname=${outdir}/${uuid}.tf_${tf}.fb_${fb}.clean.lsf.log
            if [ -f ${outseq} ]
            then
                echo "$uuid complete"
                continue
            fi

            echo "/usr/bin/time -v $mccortex clean -m $((${mem} / 2000))G --tips=$(($K * $tf)) --fallback $fb -o ${outgraph} ${graph} && $mccortex unitigs -m $((${mem} / 2000))G -n 800M ${outgraph} | gzip > ${outseq}" | bsub -J mc_cl_k${K} -o ${logfname} -We 8:00 -n ${threads} -M $mem -R "rusage[mem=${pmem}]" -R "span[hosts=1]"
        done
    done
done
