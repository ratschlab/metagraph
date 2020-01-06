#!/bin/bash

set -e

K=19

threads=24
basedir=/cluster/work/grlab/projects/metagenome/data
anno_in=${basedir}/metasub/graphs/output_k${K}_cleaned_graph_annotation.collect.column.annodbg
anno_out=${basedir}/metasub/graphs/output_k${K}_cleaned_graph_annotation.collect.relabeled

metagraph=/cluster/home/akahles/git/software/metagraph_current_dev/metagraph/build/metagraph

label_map=../complete_metadata_extended.clean.v2.relabel_map.tsv

$metagraph transform_anno -v -i ${graphdir}/graph_merged_k${K}.dbg --outfile-base ${anno_out} --parallel $threads $anno_in --rename-cols ${label_map} 2>&1 | tee $(pwd)/relabel_annotation_clean.k${K}.log
