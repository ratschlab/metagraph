#!/bin/bash

set -e

basedir=/cluster/work/grlab/projects/metagenome/data/gtex/
aligndir=${basedir}/align_samples_extended_fast

K=41

mem=1000
threads=1
pmem=1000

for fname in ${aligndir}/*.k${K}.txt.gz
do
    if [ ! -f ${fname%.txt.gz}.done ]
    then
        echo $fname incomplete
        continue
    fi

    outfname=${fname%.txt.gz}.summary.txt

    if [ ! -f ${outfname} ]
    then
        echo "python $(pwd)/summarize_mapped_samples_metagraph.py $fname > $outfname" | bsub -M $mem -n $threads -R "rusage[mem=${pmem}]" -We 20:00 -n 1 -J sum_gtex -o /dev/null
    else
        echo "${outfname} already present"
    fi
done

