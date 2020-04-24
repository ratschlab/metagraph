#!/bin/bash

set -e

mem=950000

K=41
basedir=/cluster/work/grlab/projects/metagenome/data/gtex
datadir=${basedir}/output_k${K}_trimmed_extended_annotation

metagraph=/cluster/home/akahles/git/software/metagraph_current/metagraph/build/metagraph

filelist=${datadir}.files.txt
find ${datadir} -name \*annodbg > $filelist
outbase=${basedir}/output_k${K}_trimmed_extended.abundance
lsffile=${basedir}/output_k${K}_trimmed_extended.abundance.column.lsf.log
logfile=${basedir}/output_k${K}_trimmed_extended.abundance.column.log
if [ ! -f ${outbase}.column.annodbg ]
then
    echo "cat $filelist | /usr/bin/time -v $metagraph merge_anno -o $outbase 2>&1 | tee $logfile" | bsub -M $mem -n 1 -R "rusage[mem=${mem}]" -We 48:00 -J annoGTEx -o $lsffile
fi
