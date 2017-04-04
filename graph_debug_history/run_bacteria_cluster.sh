#!/bin/bash

set -e

k=27

basedir=/cluster/project/raetsch/lab/07/shared/data/HMP/genbank/Bacteria/Bacterial_genomes/
outbase=/cluster/project/grlab/projects/metagenome/bacteria/results/${k}
mkdir -p $outbase

graphtool=/cluster/project/grlab/home/akahles/git/projects/2014/metagenome/metagraph/metagraph

for dd in $(ls -1 $basedir)
do
    if [ "$dd" != "${dd%gz}" ]
    then
        continue
    fi
    outdir=${outbase}/$dd
    logfile=${outbase}/${dd}.log
    donefile="${outbase}/${dd}.done"

    if [ ! -f ${donefile} ]
    then
        echo "time $graphtool build -v -k ${k} -O $outdir ${basedir}$dd/*.fna && touch $donefile" | bsub -M 3000 -J metag -W 12:00 -o $logfile -n 1 -R "rusage[mem=3000]" -R "span[hosts=1]" 
    fi
done
