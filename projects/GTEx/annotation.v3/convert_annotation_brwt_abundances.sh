#!/bin/bash

set -e

basedir=/cluster/work/grlab/projects/metagenome/data
metagraph=/cluster/home/akahles/git/software/metagraph_current3/metagraph/build/metagraph
K=41
mem=900000
threads=30
pmem=$((${mem} / ${threads}))

#echo "/usr/bin/time -v $metagraph transform_anno -v -o ${basedir}/gtex/output_k${K}_merged.anno --anno-type brwt --greedy ${basedir}/gtex/output_k${K}_merged.anno.column.annodbg -p $threads 2>&1" #| bsub -J convert_metasub_to_brwt_pm -oo conversion_to_brwt.k${K}.lsf -W 150:00 -n $threads -R "rusage[mem=${mem}]" -R "span[hosts=1]" 

### transform
echo "/usr/bin/time -v $metagraph transform_anno -v -o ${basedir}/gtex/output_k${K}_trimmed_clean.abundance --anno-type brwt --greedy ${basedir}/gtex/output_k${K}_trimmed_clean.abundance.column.annodbg -p $threads 2<&1 | tee ${basedir}/gtex/conversion_to_brwt_abundance.k${K}.log" | bsub -J convert_metasub_to_brwt -oo ${basedir}/gtex/conversion_to_brwt_abundance.k${K}.lsf.log -We 150:00 -n $threads -R "rusage[mem=${pmem}]" -M $mem -R "span[hosts=1]"
