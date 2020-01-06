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
mem=200000
threads=16
pmem=$((${mem} / ${threads}))

### transform
echo "/usr/bin/time -v $metagraph transform_anno -v -o ${basedir}/output_k${K}_trimmed_clean.${tissue}.samples --parallel-nodes $(($threads / 2)) --anno-type brwt --greedy ${basedir}/output_k${K}_trimmed_clean.${tissue}.samples.column.annodbg -p $threads 2<&1 | tee ${basedir}/output_k${K}_trimmed_clean.${tissue}.samples.brwt.log" | bsub -J convert_metasub_to_brwt -oo ${basedir}/output_k${K}_trimmed_clean.${tissue}.samples.brwt.lsf.log -We 48:00 -n $threads -R "rusage[mem=${pmem}]" -M $mem -R "span[hosts=1]"
