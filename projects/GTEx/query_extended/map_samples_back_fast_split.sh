#!/bin/bash

set -e

K=41
threads=8
mem=350000
pmem=$(($mem / $threads))

metadata=/cluster/work/grlab/projects/GTEx/metadata/SraRunTable_20180218.txt
metagraph=/cluster/home/akahles/git/software/metagraph_current/metagraph/build/metagraph
fastqdir=/cluster/work/grlab/SHARED_DATA_PROPOSAL_OTHERWISE_DELETE/GTEx/extract/dbGaP-9608/fastq
basedir=/cluster/work/grlab/projects/metagenome/data/gtex
outdir=${basedir}/align_samples_extended_fast
mkdir -p $outdir
tmpdir=${outdir}/tmp
mkdir -p $tmpdir

graph=${basedir}/output_k${K}_trimmed_clean_graph_extended_chunked/graph_k${K}.dbg
anno=${basedir}/output_k${K}_trimmed_extended.samples.brwt.annodbg

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
    if [ -z "${fq1}" -o -z "${fq2}" ]
    then
        echo No fastq file for $uuid
        continue
    fi
    ### create splits
    splitdir=${tmpdir}/${uuid}
    mkdir -p $splitdir
    if [ ! -f ${splitdir}/splits.done ]
    then
        echo generating splits for $fq1
        zcat $fq1 | split -l 1000000 --additional-suffix .fastq - ${splitdir}/${uuid}_1.
        for fname in ${splitdir}/${uuid}_1*fastq
        do
            echo compressing $fname
            gzip $fname
        done
        echo generating splits for $fq2
        zcat $fq2 | split -l 1000000 --additional-suffix .fastq - ${splitdir}/${uuid}_2.
        for fname in ${splitdir}/${uuid}_2*fastq
        do
            echo compressing $fname
            gzip $fname
        done
        touch ${splitdir}/splits.done
    fi

    for fname in ${splitdir}/*.fastq.gz
    do
        fbase=$(basename $fname)
        fbase=${fbase%.fastq.gz}
        resultdir=${outdir}/results_split/${uuid}
        mkdir -p $resultdir
        outfile=${resultdir}/remap.${fbase}.k${K}.txt.gz
        donefile=${resultdir}/remap.${fbase}.k${K}.done
        logfile=${resultdir}/remap.${fbase}.k${K}.lsf.log
        if [ -f ${donefile} ]
        then
            echo "$fbase already done"
            continue
        fi
        echo "(/usr/bin/time -v $metagraph query --discovery-fraction 0.0 --query-mode matches -v -i $graph -a $anno -p $threads $fname | gzip > $outfile) && touch $donefile" | bsub -M $mem -n $threads -R "rusage[mem=${pmem}]" -We 48:00 -n 1 -J map_gtex -oo $logfile  #-R "select[model==XeonGold_6150]"
    done
    exit
done < <(shuf -n 10 --random-source $metadata $metadata | tr $'\t' ',')
