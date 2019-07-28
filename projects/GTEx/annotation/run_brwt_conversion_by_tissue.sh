#!/bin/bash

set -e 

K=41
mem=300000
threads=16
pmem=$(($mem / $threads))

annograph=/cluster/home/akahles/git/software/genome_graph_annotation/build/annograph

basedir=/cluster/work/grlab/projects/metagenome/data/gtex
annodir=${basedir}/output_k${K}_merged_annotated_tissue
outdir=${basedir}/output_k${K}_merged_annotated_tissue
mkdir -p $outdir
all_metadata=${basedir}/metadata/SraRunTable_20180218.txt

for tissue in $(tail -n+2 $all_metadata | cut -f 24 | sort -u | tr ' ' '_')
do
    if [ "$tissue" == "histological_type" ]
    then
        continue
    fi
    logfile=${outdir}/output_k${K}_merged.${tissue}.brwt.lsf
    if [ ! -f "${annodir}/output_k${K}_merged.${tissue}.column.annodbg" ]
    then
        continue
    fi
    if [ "$1" == "local" ]
    then
        echo processing $tissue, writing log to $logfile
        /usr/bin/time -v $annograph transform_anno -v -o ${outdir}/output_k${K}_merged.${tissue} --anno-type brwt --greedy ${annodir}/output_k${K}_merged.${tissue}.column.annodbg -p $threads > $logfile 2>&1 
    else
        echo "/usr/bin/time -v $annograph transform_anno -v -o ${outdir}/output_k${K}_merged.${tissue} --anno-type brwt --greedy ${annodir}/output_k${K}_merged.${tissue}.column.annodbg -p $threads > $logfile 2>&1" | bsub -M $mem -n 1 -R "rusage[mem=${mem}]" -We 24:00 -J annoGTEx -o /dev/null
    fi
done
