#!/bin/bash

set -e

if [ -z "$1" ]
then
    echo "Usage: $0 <K>"
    exit 1
fi
K="$1"

mem=450000
threads=1
pmem=$(($mem / $threads))

basedir=/cluster/work/grlab/projects/metagenome/data
annodir=${basedir}/metasub/graphs/output_k${K}_cleaned_graph_annotation
metagraph=/cluster/home/akahles/git/software/metagraph_current_dev/metagraph/build/metagraph

echo "/usr/bin/time -v ${metagraph} merge_anno -v -o ${annodir}.collect ${annodir}/*.annodbg 2>&1 | tee ${annodir}.collect.log" | bsub -J ms_anno_k${K} -We 24:00 -n $threads -M $mem -R "rusage[mem=${pmem}]" -R "span[hosts=1]" -oo ${annodir}.collect.lsf.log

