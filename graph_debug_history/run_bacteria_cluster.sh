#!/bin/bash

set -e

k=27

basedir=/cbio/shared/data/HMP/genbank/Bacteria/Bacterial_genomes/
outbase=/cbio/grlab/projects/metagenome/bacteria/results/${k}
mkdir -p $outbase

graphtool=/cbio/grlab/home/akahles/git/projects/2014/metagenome/seqan_graph/metagraph

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
        echo "time $graphtool -v -k ${k} -O $outdir ${basedir}$dd/*.fna && touch $donefile" | qsub -l nodes=1:ppn=1,mem=3G,vmem=3G,pmem=3G,walltime=4:00:00 -N meta -j oe -o $logfile
    fi
done
