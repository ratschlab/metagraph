#!/bin/bash

set -e

K=41

basedir=/cluster/work/grlab/projects/metagenome/data
outdir=${basedir}/metasub/graphs/output_k${K}
mkdir -p $outdir
metagraph=/cluster/home/akahles/git/software/metagraph_clean/metagraph/build/metagraph

metadata=$(pwd)/../complete_metadata_extended.clean.v2.csv
N=$(wc -l $metadata | cut -f 1 -d ' ')

for i in $(seq 1 $N)
do
    uuid="$(sed -n ${i}p ${metadata} | cut -f 2 -d ',')"
    if [ ! -f ${outdir}/${uuid}.dbg ]
    then
        echo $i
    fi
done
