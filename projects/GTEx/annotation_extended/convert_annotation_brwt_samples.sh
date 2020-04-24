#!/bin/bash

set -e

basedir=/cluster/work/grlab/projects/metagenome/data
metagraph=/cluster/home/akahles/git/software/metagraph_current/metagraph/build/metagraph
K=41
mem=1200000
threads=24
pmem=$((${mem} / ${threads}))

### transform
echo "/usr/bin/time -v $metagraph transform_anno -v -o ${basedir}/gtex/output_k${K}_trimmed_extended.samples --parallel-nodes $(($threads / 6)) --anno-type brwt --greedy ${basedir}/gtex/output_k${K}_trimmed_extended.samples.column.annodbg -p $threads 2<&1 | tee ${basedir}/gtex/conversion_to_brwt_extended_samples.k${K}.log" | bsub -J convert_metasub_to_brwt -oo ${basedir}/gtex/conversion_to_brwt_extended_samples.k${K}.lsf.log -We 150:00 -n $threads -R "rusage[mem=${pmem}]" -M $mem -R "span[hosts=1]"
