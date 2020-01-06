#!/bin/bash

set -e

if [ -z "$1" ]
then
    echo "Usage: $0 <tissue>"
    exit 1
fi
tissue=$1

basedir=/cluster/work/grlab/projects/metagenome/data/tcga
metagraph=/cluster/home/akahles/git/software/metagraph_current_dev/metagraph/build/metagraph
K=31
mem=2000000
threads=48
pmem=$((${mem} / ${threads}))

outbase=${basedir}/output_k${K}_trimmed_clean.${tissue}.genomebins

### transform
echo "/usr/bin/time -v $metagraph transform_anno -v -o ${outbase} --parallel-nodes 2 --anno-type brwt --greedy ${outbase}.column.annodbg -p $threads 2<&1 | tee ${outbase}.brwt.log" | bsub -J convert_metasub_to_brwt -oo ${outbase}.brwt.lsf.log -We 24:00 -n $threads -R "rusage[mem=${pmem}]" -M $mem -R "span[hosts=1]"
