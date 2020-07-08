#!/bin/bash

set -e

if [ -z "$1" ]
then
    echo "Usage: $0 <k>"
    exit 1
fi
K=$1

basedir=/cluster/work/grlab/projects/metagenome/data
graphbase=${basedir}/metasub/graphs/output_k${K}_cleaned_graph_annotation.collect.relabeled
metagraph=/cluster/home/akahles/git/software/metagraph_current_dev/metagraph/build/metagraph
mem=700000
threads=20
pmem=$(($mem / $threads))

echo "/usr/bin/time -v $metagraph relax_brwt -v -o ${graphbase} --relax-arity 15 -p ${threads} ${graphbase}.brwt.annodbg 2>&1 | tee ${graphbase}.relax_brwt.log" | bsub -G ms_raets -J rl_brwt_k${K} -oo ${graphbase}.relax_brwt.lsf.log -W 120:00 -n $threads -M ${mem} -R "rusage[mem=${pmem}]" -R "span[hosts=1]" 
