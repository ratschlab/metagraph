#!/bin/bash

set -e

k=27
level=2
stage=8
mem=1500
cpus=24

basedir=/cluster/project/grlab/projects/metagenome/bacteria/results/${k}_merge$(($stage - 1))
outbase=/cluster/project/grlab/projects/metagenome/bacteria/results/${k}_merge${stage}
mkdir -p $outbase

graphtool=/cluster/project/grlab/home/akahles/git/projects/2014/metagenome/metagraph/metagraph
cnt=0
total=0
for dd in $(ls -1 $basedir/*.done)
do
    
    outdir=${outbase}/${dd%.done}
    if [ "${cnt}" == "0" ]
    then
        mergelist=${dd%.done}
    else
        mergelist="${mergelist} ${dd%.done}"
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
        echo "time $graphtool merge -v --parallel $cpus --bins-per-thread 2 -O $outfile $mergelist && touch $donefile" | bsub -M $((${mem} * ${cpus})) -J metag -W 72:00 -o $logfile -n $cpus -R "rusage[mem=${mem}]" -R "span[hosts=1]" 
    else
        echo "$donefile already exists"
    fi
done

if [ "${cnt}" != "0" ]
then
    logfile=${outbase}/merge_${total}.log
    donefile="${outbase}/merge_${total}.done"
    outfile="${outbase}/merge_${total}"
    if [ ! -f ${donefile} ]
    then
        echo "time $graphtool merge -v --parallel $cpus --bins-per-thread 2 -O $outfile $mergelist && touch $donefile" | bsub -M $((${mem} * ${cpus})) -J metag -W 72:00 -o $logfile -n $cpus -R "rusage[mem=${mem}]" -R "span[hosts=1]" 
    else
        echo "$donefile already exists"
    fi
fi
