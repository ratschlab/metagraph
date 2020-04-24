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
mem=600000
threads=32
pmem=$((${mem} / ${threads}))

outbase=${basedir}/output_k${K}_trimmed_clean.${tissue}.genomebins_per_chr

### transform
for fname in ${outbase}/*.column.annodbg
do
    fbase=${fname%.column.annodbg}
    echo "/usr/bin/time -v $metagraph transform_anno -v -o ${fbase} --parallel-nodes 2 --anno-type brwt --greedy ${fname} -p $threads 2<&1 | tee ${fbase}.brwt.log" | bsub -J convert_metasub_to_brwt -oo ${fbase}.brwt.lsf.log -We 24:00 -n $threads -R "rusage[mem=${pmem}]" -M $mem -R "span[hosts=1]"
done
