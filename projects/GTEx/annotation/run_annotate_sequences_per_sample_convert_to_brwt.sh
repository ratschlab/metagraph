#!/bin/bash

set -e

### source paths
# set basedir
# set metagraph
. ../paths.sh

K=41
mem=600000
threads=24
pmem=$((${mem} / ${threads}))

### transform
echo "/usr/bin/time -v $metagraph transform_anno -v -o ${basedir}/gtex/graphs/output_k${K}_trimmed_clean.samples --parallel-nodes $(($threads / 2)) --anno-type brwt --greedy ${basedir}/gtex/graphs/output_k${K}_trimmed_clean.samples.column.annodbg -p $threads 2<&1 | tee ${basedir}/gtex/graphs/conversion_to_brwt_samples.k${K}.log" | bsub -J convert_metasub_to_brwt -oo ${basedir}/gtex/graphs/conversion_to_brwt_samples.k${K}.lsf.log -We 150:00 -n $threads -R "rusage[mem=${pmem}]" -M $mem -R "span[hosts=1]"
