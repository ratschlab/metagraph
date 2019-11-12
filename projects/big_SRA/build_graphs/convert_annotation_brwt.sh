#!/bin/bash

set -e

if [ -z "$1" ]
then
    echo "Usage: $0 <k>"
    exit 1
fi
K=$1

basedir=/cluster/work/grlab/projects/metagenome/data
graphbase=${basedir}/metagut/graphs/output_k${K}_annotation_clean.collect.relabeled
metagraph=/cluster/home/akahles/git/software/metagraph_current/metagraph/build/metagraph
#mem=1500000
mem=1600000
threads=24
pmem=$(($mem / $threads))

echo "/usr/bin/time -v $metagraph transform_anno -v -o ${graphbase} --anno-type brwt --parallel-nodes $((${threads} / 2)) --greedy ${graphbase}.column.annodbg -p ${threads} 2>&1 | tee ${graphbase}.brwt.log" | bsub -J mg_brwt_k${K} -oo ${graphbase}.lsf.log -W 120:00 -n $threads -M ${mem} -R "rusage[mem=${pmem}]" -R "span[hosts=1] -oo ${graphbase}.brwt.lsf.log" 
#/usr/bin/time -v $metagraph transform_anno -v -o ${graphbase} --anno-type brwt --greedy ${graphbase}.columns.annodbg -p ${threads} 2>&1 | tee conversion_to_brwt.k${K}.log 
