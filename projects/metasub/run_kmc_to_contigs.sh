#!/bin/bash

set -e 

K=15
mem=25000
threads=4
pmem=$(($mem / $threads))

outdir=/cluster/work/grlab/projects/metagenome/results/metasub_wasabi
mkdir -p $outdir

while IFS='' read -r line || [[ -n "$line" ]]
do
    uuid="$(echo $line | cut -f 2 -d ',')"
    if [ "$uuid" == "uuid" ]
    then
        continue
    fi
    kmc_file=${outdir}/${uuid}.kmc_suf
    if [ ! -f "${kmc_file}" ]
    then
        echo "KMC not complete for $uuid"
        continue
    fi
    outfile=${outdir}/${uuid}.k${K}.contigs.fasta.gz
    if [ ! -f ${outfile} ]
    then
        echo "$(pwd)/../../metagraph/scripts/kmc_to_contigs_2step.sh $kmc_file $K" | bsub -M $mem -n $threads -R "rusage[mem=${pmem}]" -We 12:00 -n 1 -J kmc_msub1 -o ${outdir}/${uuid}.kmc2contS1.k${K}.cluster.log 
    fi
done < complete_metadata_extended.clean.csv
