#!/bin/bash

set -e

fastqdir=/cluster/work/grlab/SHARED_DATA_PROPOSAL_OTHERWISE_DELETE/GTEx/extract/dbGaP-9608/fastq
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
    fq1=${fastqdir}/${uuid}_1.fastq.gz
    fq2=${fastqdir}/${uuid}_2.fastq.gz
    for fq in $fq1 $fq2
    do 
        out=${fq}.fqchk.txt.gz
        if [ ! -f ${out} ]
        then
            echo "module load seqtk; seqtk fqchk -q0 $fq | gzip > $out" | bsub -o /dev/null -M 2000 -W 1:00 -n 1 -J stats
        fi
    done
done < <(cat /cluster/work/grlab/projects/GTEx/metadata/SraRunTable_20180218.txt | tr $'\t' ',')
