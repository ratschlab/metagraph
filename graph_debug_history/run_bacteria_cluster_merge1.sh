#!/bin/bash

set -e

k=27
level=3

basedir=/cbio/grlab/projects/metagenome/bacteria/results/${k}
outbase=/cbio/grlab/projects/metagenome/bacteria/results/${k}_merge1
mkdir -p $outbase

graphtool=/cbio/grlab/home/akahles/git/projects/2014/metagenome/seqan_graph/metagraph
cnt=0
total=0
for dd in $(ls -1 $basedir/*.done)
do
    
    outdir=${outbase}/${dd%.done}
    if [ "${cnt}" == "0" ]
    then
        mergelist=${dd%.done}
    else
        mergelist="${mergelist},${dd%.done}"
    fi
    cnt=$(($cnt + 1))
    total=$(($total + 1))
    if [ "${cnt}" == "${level}" ]
    then
        cnt=0
    else
        continue
    fi
    logfile=${outbase}/merge_${total}.log
    donefile="${outbase}/merge_${total}.done"
    outfile="${outbase}/merge_${total}"
    if [ ! -f ${donefile} ]
    then
        echo "time $graphtool -v -m $mergelist -O $outfile DUMMY && touch $donefile" | qsub -l nodes=1:ppn=1,mem=4G,vmem=4G,pmem=4G,walltime=24:00:00 -N meta -j oe -o $logfile
    fi
done
