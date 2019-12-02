#!/bin/bash

set -e 

K=41
threads=2
metagraph=/cluster/home/akahles/git/software/metagraph_brwt_rrr/metagraph/build/metagengraph

basedir=/cluster/work/grlab/projects/metagenome/data/gtex
query_dir=${basedir}/queries
datadir=/cluster/work/grlab/projects/GTEx/rna/results/alignments

for fname in ${datadir}/*.all.unaligned.sorted.r1.fq.gz
do
    sample=$(basename $fname | cut -f 1 -d '.')
    outfile=${query_dir}/query_unlabeled.${sample}.result.txt
    if [ ! -f ${outfile} ]
    then
        $metagraph query -p ${threads} --discovery-fraction 1.0 -i ${basedir}/output_k${K}_merged/all_merged_k${K}.dbg -a ${basedir}/output_k${K}_merged.anno.brwt.annodbg ${fname} > $outfile
    fi
done
