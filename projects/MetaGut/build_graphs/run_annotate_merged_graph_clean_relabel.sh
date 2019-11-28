#!/bin/bash

set -e

if [ -z "$1" ]
then
    echo "Usage: $0 <k>"
    exit 1
fi
K="$1"

mem=1200000
threads=1
pmem=$(($mem / $threads))

basedir=/cluster/work/grlab/projects/metagenome/data
anno_in=${basedir}/metagut/graphs/output_k${K}_annotation_clean.collect.column.annodbg
anno_out=${basedir}/metagut/graphs/output_k${K}_annotation_clean.collect.relabeled

metagraph=/cluster/home/akahles/git/software/metagraph_current/metagraph/build/metagraph

### generate label map
$metagraph stats -a ${anno_in} --print-col-names | grep -v INFO > ${anno_in}.labels
rev ${anno_in}.labels | cut -f 1 -d '/' | rev | cut -f 1 -d '.' > k${K}_labels_new
paste ${anno_in}.labels k${K}_labels_new > k${K}_labels_map

echo "cd $(pwd); $metagraph transform_anno -v --outfile-base ${anno_out} --parallel $threads $anno_in --rename-cols k${K}_labels_map 2>&1" | bsub -J mg_anno_k${K} -We 24:00 -n $threads -M $mem -R "rusage[mem=${pmem}]" -R "span[hosts=1]" -oo ${anno_out}.lsf.log
 
