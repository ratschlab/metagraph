#!/bin/bash

set -e

k=27
level=2
stage=4
mem=5G

basedir=/cbio/grlab/projects/metagenome/bacteria/results/${k}_merge$(($stage - 1))
outbase=/cbio/grlab/projects/metagenome/bacteria/results/${k}_merge${stage}
mkdir -p $outbase

graphtool=/cbio/grlab/home/akahles/git/projects/2014/metagenome/seqan_graph/metagraph2
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
        echo "time $graphtool -v -m $mergelist -O $outfile DUMMY && touch $donefile" | qsub -l nodes=1:ppn=1,mem=${mem},vmem=${mem},pmem=${mem},walltime=72:00:00 -N meta_${stage} -j oe -o $logfile
    else
        echo "$donefile already exists"
    fi
done
