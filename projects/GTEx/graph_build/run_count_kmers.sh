#!/bin/bash

set -e 

K=41
mem=75000
threads=2
pmem=$(($mem / $threads))

fastqdir=/cluster/work/grlab/SHARED_DATA_PROPOSAL_OTHERWISE_DELETE/GTEx/extract/dbGaP-9608/fastq
outdir=/cluster/work/grlab/projects/metagenome/data/gtex/graphs/output_k${K}_trimmed
mkdir -p $outdir
metadata=/cluster/work/grlab/projects/GTEx/metadata/SraRunTable_20180218.txt

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

    if [ -f "${outdir}/${uuid}.k${K}.kmc_suf" ]
    then
        echo ${uuid} already done
        continue
    fi
    fq1=${fastqdir}/${uuid}_1.fastq.gz
    fq2=${fastqdir}/${uuid}_2.fastq.gz
    if [ -z "${fq1}" -o -z "${fq2}" ]
    then
        echo No fastq file for $uuid
        continue
    fi
    echo "/usr/bin/time -v $(pwd)/count_kmers.sh $K $uuid $outdir $fq1 $fq2" | bsub -M $mem -n $threads -R "rusage[mem=${pmem}]" -We 12:00 -n 1 -J kmc_gtex -o ${outdir}/${uuid}.k${K}.cluster.log
done < <(cat $metadata | tr $'\t' ',')
