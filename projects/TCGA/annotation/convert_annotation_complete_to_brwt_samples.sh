#!/bin/bash

set -e

basedir=/cluster/work/grlab/projects/metagenome/data/tcga
metagraph=/cluster/home/akahles/git/software/metagraph_current_dev/metagraph/build/metagraph
K=31
mem=600000
threads=16
pmem=$((${mem} / ${threads}))

### transform
echo "/usr/bin/time -v $metagraph transform_anno -v -o ${basedir}/graph_merged_complete_k${K}.samples --parallel-nodes $(($threads / 2)) --anno-type brwt --greedy ${basedir}/graph_merged_complete_k${K}.samples.column.annodbg -p $threads 2<&1 | tee ${basedir}/graph_merged_complete_k${K}.samples.brwt.log" | bsub -G ms_raets -J convert_metasub_to_brwt -oo ${basedir}/graph_merged_complete_k${K}.samples.brwt.lsf.log -We 48:00 -n $threads -R "rusage[mem=${pmem}]" -M $mem -R "span[hosts=1]"
