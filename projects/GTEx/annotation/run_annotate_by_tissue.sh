#!/bin/bash

set -e 

K=41
mem=100000
threads=1
pmem=$(($mem / $threads))

basedir=/cluster/work/grlab/projects/metagenome/data/gtex
annodir=${basedir}/output_k${K}_merged_annotated
outdir=${basedir}/output_k${K}_merged_annotated_tissue
mkdir -p $outdir
all_metadata=${basedir}/metadata/SraRunTable_20180218.txt

for tissue in $(tail -n+2 $all_metadata | cut -f 24 | sort -u | tr ' ' '_')
do
    files=""
    while IFS=',' read -r -a line || [[ -n "$line" ]]
    do
        assay=${line[0]}
        pair=${line[8]}
        uuid=${line[16]}
        if [ "$assay" == "Assay_Type" ]
        then 
            continue
        fi  
        if [ "$assay" != "RNA-Seq" ]
        then 
            continue
        fi  
        if [ "$pair" != "PAIRED" ]
        then
            continue
        fi  

        files="$files ${annodir}/all_merged_k${K}.${uuid}.column.annodbg"
    done < <(cat ${basedir}/metadata/by_tissue/SraRunTable_20180218.${tissue}.txt | tr $'\t' ',')
    
    logfile=${outdir}/output_k${K}_merged.${tissue}.column.lsf
    if [ -f "${outdir}/output_k${K}_merged.${tissue}.column.annodbg" ]
    then
        continue
    fi
    if [ "$1" == "local" ]
    then
        echo processing $tissue, writing log to $logfile
        /usr/bin/time -v $(pwd)/../../../metagraph/build/metagengraph merge_anno -v -o ${outdir}/output_k${K}_merged.${tissue} $files > $logfile 2>&1
    else
        echo "/usr/bin/time -v $(pwd)/../../../metagraph/build/metagengraph merge_anno -v -o ${outdir}/output_k${K}_merged.${tissue} $files > $logfile 2>&1" | bsub -M $mem -n 1 -R "rusage[mem=${mem}]" -We 4:00 -J annoGTEx -o /dev/null
    fi
done
