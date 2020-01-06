#!/bin/bash

set -e

if [ -z "$1" ]
then
    echo "Usage: $0 <tissue>"
    exit 1
fi
tissue=$1

mem=350000
K=31
basedir=/cluster/work/grlab/projects/metagenome/data/tcga
datadir=${basedir}/output_k${K}_trimmed_annotation_per_sample/${tissue}

metagraph=/cluster/home/akahles/git/software/metagraph_current_dev/metagraph/build/metagraph

filelist=${datadir}.files.txt
find ${datadir} -name \*.annodbg > $filelist

outbase=${basedir}/output_k${K}_trimmed_clean.samples
lsffile=${basedir}/output_k${K}_trimmed_clean.samples.column.lsf.log
logfile=${basedir}/output_k${K}_trimmed_clean.samples.column.log
if [ ! -f ${outbase}.column.annodbg ]
then
    echo "cat $filelist | /usr/bin/time -v $metagraph merge_anno -o $outbase 2>&1 | tee $logfile" | bsub -M $mem -n 1 -R "rusage[mem=${mem}]" -We 48:00 -J annoTCGA -o $lsffile
fi
