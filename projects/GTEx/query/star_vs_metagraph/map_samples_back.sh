#!/bin/bash

set -e

K=41
threads=8
mem=80000
pmem=$(($mem / $threads))

metadata=/cluster/work/grlab/projects/GTEx/metadata/SraRunTable_20180218.txt
metagraph=/cluster/home/akahles/git/software/metagraph_current/metagraph/build/metagraph
fastqdir=/cluster/work/grlab/SHARED_DATA_PROPOSAL_OTHERWISE_DELETE/GTEx/extract/dbGaP-9608/fastq
basedir=/cluster/work/grlab/projects/metagenome/data/gtex
outdir=${basedir}/align_samples
mkdir -p $outdir

graph=${basedir}/output_k${K}_trimmed_clean_graph_chunked/graph_merged_k${K}.dbg

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

    outfile=${outdir}/remap.${uuid}.k${K}.txt.gz
    donefile=${outdir}/remap.${uuid}.k${K}.done
    logfile=${outdir}/remap.${uuid}.k${K}.lsf.log
    if [ -f ${donefile} ]
    then
        echo "$uuid already done"
        continue
    fi

    fq1=${fastqdir}/${uuid}_1.fastq.gz
    fq2=${fastqdir}/${uuid}_2.fastq.gz
    if [ -z "${fq1}" -o -z "${fq2}" ]
    then
        echo No fastq file for $uuid
        continue
    fi
    echo "(/usr/bin/time -v $metagraph align --count-kmers -i $graph -p $threads $fq1 $fq2 | cut -f 2 -d':' | tr -d ' ' | cut -f 1 -d '/' | gzip > $outfile) && touch $donefile" | bsub -M $mem -n $threads -R "rusage[mem=${pmem}]" -We 48:00 -n 1 -J map_gtex -oo $logfile  -R "select[model==XeonGold_6150]"
done < <(shuf -n 50 --random-source $metadata $metadata | tr $'\t' ',')
