#!/bin/bash

set -e

if [ -z "$1" ]
then
    echo "Usage: $0 <K>"
    exit 1
fi
K="$1"

#mem=300000
#threads=2
#pmem=$(($mem / $threads))
basedir=/cluster/work/grlab/projects/metagenome/data
annodir=${basedir}/metagut/graphs/output_k${K}_annotation_clean
metagraph=/cluster/home/akahles/git/software/metagraph_current3/metagraph/build/metagraph

/usr/bin/time -v ${metagraph} merge_anno -v -o ${annodir}.collect ${annodir}/*.annodbg | tee > ${annodir}.collect.log


