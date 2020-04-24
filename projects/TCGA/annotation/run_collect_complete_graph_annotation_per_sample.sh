#!/bin/bash

set -e

mem=450000
K=31
basedir=/cluster/work/grlab/projects/metagenome/data/tcga
datadir=${basedir}/graph_merged_complete_k${K}_annotation_per_sample

metagraph=/cluster/home/akahles/git/software/metagraph_current_dev/metagraph/build/metagraph

filelist=${datadir}.files.txt
find ${datadir} -name \*.annodbg > $filelist

outbase=${basedir}/graph_merged_complete_k${K}.samples
lsffile=${basedir}/graph_merged_complete_k${K}.samples.column.lsf.log
logfile=${basedir}/graph_merged_complete_k${K}.samples.column.log
if [ ! -f ${outbase}.column.annodbg ]
then
    echo "cat $filelist | /usr/bin/time -v $metagraph merge_anno -o $outbase 2>&1 | tee $logfile" | bsub -G ms_raets -M $mem -n 1 -R "rusage[mem=${mem}]" -We 48:00 -J annoTCGA -o $lsffile
fi
