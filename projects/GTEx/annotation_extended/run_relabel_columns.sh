#!/bin/bash

set -e

if [ -z "$1" ]
then
    echo "Usage: $0 <k>"
    exit 1
fi
K="$1"

mem=910000
threads=24
pmem=$(($mem / $threads))

basedir=/cluster/work/grlab/projects/metagenome/data
anno_in=${basedir}/gtex/output_k${K}_trimmed_extended.abundance.column.annodbg
anno_out=${basedir}/gtex/output_k${K}_trimmed_extended.samples.column.annodbg

metagraph=/cluster/home/akahles/git/software/metagraph_current/metagraph/build/metagraph

### generate label map
$metagraph stats -a ${anno_in} --print-col-names | grep -v INFO | cut -f 5 -d ' ' > ${anno_in}.labels
cat ${anno_in}.labels | rev | cut -f 1 -d '/' | rev | cut -f 1 -d '.' > k${K}_labels_new
paste ${anno_in}.labels k${K}_labels_new > k${K}_labels_map

echo "cd $(pwd); $metagraph transform_anno -v --outfile-base ${anno_out} --parallel $threads $anno_in --rename-cols k${K}_labels_map 2>&1" | bsub -J ms_anno_k${K} -We 72:00 -n $threads -M $mem -R "rusage[mem=${pmem}]" -R "span[hosts=1]" -oo ${anno_out}.lsf.log
 
