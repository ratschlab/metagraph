#!/bin/bash

set -e

K=19

#mem=240000
#threads=18
#mem=500000
threads=24
#pmem=$(($mem / $threads))
basedir=/cluster/work/grlab/projects/metagenome/data
anno_in=${basedir}/metasub/graphs/output_k${K}_cleaned_graph_annotation_clean.collect.column.annodbg
anno_out=${basedir}/metasub/graphs/output_k${K}_cleaned_graph_annotation_clean.collect.relabeled
label_map=../complete_metadata_extended.clean.v2.relabel_map.tsv

metagraph=/cluster/home/akahles/git/software/metagraph_current/metagraph/build/metagraph

$metagraph transform_anno -v -i ${graphdir}/graph_merged_k${K}.dbg --outfile-base ${anno_out} --parallel $threads $anno_in --rename-cols ${label_map} 2>&1 | tee $(pwd)/relabel_annotation_clean.k${K}.log
