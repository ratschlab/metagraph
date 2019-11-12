#!/bin/bash

set -e 

K=41
#mem=32000
mem=100000
threads=4
pmem=$(($mem / $threads))

fastqdir=/cluster/work/grlab/projects/metagenome/data/gtex/parameterstudy_largest_trimmed/fastq_trimmed
outdir=/cluster/work/grlab/projects/metagenome/data/gtex/parameterstudy_largest_trimmed/fastq_trimmed_kmc_k${K}
mkdir -p $outdir
fnames_dir=${outdir}/fnames_kmc
mkdir -p $fnames_dir

metadata="samples_parameterstudy.largest.tsv"

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
    if [ -z "${fastqdir}/${uuid}_1.fastq.gz" ]
    then
        echo No fastq file for $uuid
        continue
    fi
    fnames_file=${fnames_dir}/${uuid}.txt
    files=$(ls -1 ${fastqdir}/${uuid}*fastq.gz)
    echo $files | tr ' ' $'\n' > $fnames_file
    echo "$(pwd)/count_kmers.sh $K $fnames_file $outdir $uuid" | bsub -M $mem -n $threads -R "rusage[mem=${pmem}]" -We 12:00 -n 1 -J kmc_gtex -o ${outdir}/${uuid}.k${K}.cluster.log
done < <(cat ${metadata} | tr $'\t' ',')
