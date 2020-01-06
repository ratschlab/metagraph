#!/bin/bash

set -e

### source paths
# set basedir
# set metagraph
. ../paths.sh

mem=450000

K=41
datadir=${basedir}/gtex/graphs/output_k${K}_trimmed_annotation_per_abundance_quantiles

filelist=${datadir}.files.txt
find ${datadir} -name \*.k${K}.clean\*annodbg > $filelist

outbase=${basedir}/gtex/graphs/output_k${K}_trimmed_clean.abundance
lsffile=${basedir}/gtex/graphs/output_k${K}_trimmed_clean.abundance.column.lsf.log
logfile=${basedir}/gtex/graphs/output_k${K}_trimmed_clean.abundance.column.log
if [ ! -f ${outbase}.column.annodbg ]
then
    echo "cat $filelist | /usr/bin/time -v $metagraph merge_anno -o $outbase 2>&1 | tee $logfile" | bsub -M $mem -n 1 -R "rusage[mem=${mem}]" -We 48:00 -J annoGTEx -o $lsffile
fi
