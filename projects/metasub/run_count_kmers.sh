#!/bin/bash

set -e 

K=19
mem=32000
threads=4
pmem=$(($mem / $threads))

outdir=/cluster/work/grlab/projects/metagenome/data/metasub/kmc_counts/output_k${K}
mkdir -p $outdir
fnames_dir=${outdir}/fnames_kmc
mkdir -p $fnames_dir

while IFS='' read -r line || [[ -n "$line" ]]
do
    uuid="$(echo $line | cut -f 2 -d ',')"
    if [ "$uuid" == "uuid" ]
    then
        continue
    fi
    if [ -f "${outdir}/${uuid}.${K}.kmc_suf" ]
    then
        echo ${uuid} already done
        continue
    fi
    fnames_file=${fnames_dir}/${uuid}.txt
    echo $line | cut -f 41 -d ',' | tr ':' '\n' > $fnames_file
    echo "$(pwd)/count_kmers.sh $K $fnames_file $outdir $uuid" | bsub -M $mem -n $threads -R "rusage[mem=${pmem}]" -We 12:00 -n 1 -J kmc_metasub -o ${outdir}/${uuid}.k${K}.cluster.log
done < complete_metadata_extended.clean.csv
