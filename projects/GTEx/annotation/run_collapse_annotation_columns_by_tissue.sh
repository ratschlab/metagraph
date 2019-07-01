#!/bin/bash

set -e 

K=41
mem=2500000
threads=1
pmem=$(($mem / $threads))

basedir=/cluster/work/grlab/projects/metagenome/data/gtex
seqdir=${basedir}/output_k${K}_sequences
outdir=${basedir}/output_k${K}_merged_annotated_collapsed
mkdir -p $outdir
all_metadata=${basedir}/metadata/SraRunTable_20180218.txt
label_file=${outdir}/label_map.txt

### generate label file
if [ ! -f ${label_file} ]
then
    while IFS=',' read -r -a line || [[ -n "$line" ]]
    do
        assay=${line[0]}
        pair=${line[8]}
        uuid=${line[16]}
        tissue=$(echo ${line[23]} | tr ' ' '_')
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

        echo "${seqdir}/${uuid}.k${K}.sequences.fasta.gz $tissue" >> $label_file
    done < <(cat ${all_metadata} | tr $'\t' ',')
fi

logfile=${outdir}/output_k${K}_merged.collapse_anno.lsf
if [ "$1" == "local" ]
then
    /usr/bin/time -v $(pwd)/../../../metagraph/build/metagraph transform_anno -v -o ${outdir}/output_k${K}_merged.collapsed --rename-cols $label_file ${basedir}/output_k${K}_merged.anno.column.annodbg > $logfile 2>&1
else
    echo "/usr/bin/time -v $(pwd)/../../../metagraph/build/metagraph transform_anno -v -o ${outdir}/output_k${K}_merged.collapsed --rename-cols $label_file ${basedir}/output_k${K}_merged.anno.column.annodbg > $logfile 2>&1" | bsub -M $mem -n 1 -R "rusage[mem=${mem}]" -We 4:00 -J annoGTEx -o /dev/null
fi
