#!/bin/bash

set -e 

## source paths
# set metagraph
# set basedir
. ../../paths.sh

K=41
threads=2

query_dir=${basedir}/gtex/queries
datadir=/cluster/work/grlab/projects/GTEx/rna/results/alignments

for fname in ${datadir}/*.all.unaligned.sorted.r1.fq.gz
do
    sample=$(basename $fname | cut -f 1 -d '.')
    outfile=${query_dir}/query_unaligned.${sample}.result.txt
    if [ ! -f ${outfile} ]
    then
        graph=${basedir}/gtex/graphs/output_k${K}_merged/all_merged_k${K}.dbg
        annotation=${basedir}/gtex/graphs/output_k${K}_merged.anno.brwt.annodbg
        $metagraph query -p ${threads} --min-kmers-fraction-label 1.0 -i ${graph} -a ${annotation} ${fname} > $outfile
    fi
done
