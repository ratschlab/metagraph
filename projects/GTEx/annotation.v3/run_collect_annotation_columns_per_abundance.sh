#!/bin/bash

set -e

mem=450000

K=41
basedir=/cluster/work/grlab/projects/metagenome/data/gtex
datadir=${basedir}/output_k${K}_trimmed_annotation_per_abundance_quantiles

metagraph=/cluster/home/akahles/git/software/metagraph_current/metagraph/build/metagraph

filelist=${datadir}.files.txt
find ${datadir} -name \*.k${K}.clean\*annodbg > $filelist
#for level in 0 #$(seq 0 4)
#do
    outbase=${basedir}/output_k${K}_trimmed_clean.abundance
    lsffile=${basedir}/output_k${K}_trimmed_clean.abundance.column.lsf.log
    logfile=${basedir}/output_k${K}_trimmed_clean.abundance.column.log
    if [ ! -f ${outbase}.column.annodbg ]
    then
        echo "cat $filelist | /usr/bin/time -v $metagraph merge_anno -o $outbase 2>&1 | tee $logfile" | bsub -M $mem -n 1 -R "rusage[mem=${mem}]" -We 48:00 -J annoGTEx -o $lsffile
        #echo "/usr/bin/time -v $metagraph merge_anno -o $outbase ${datadir}/*.k${K}.clean*annodbg 2>&1 | tee $logfile" | bsub -M $mem -n 1 -R "rusage[mem=${mem}]" -We 48:00 -J annoGTEx -o $lsffile
    fi
#done
