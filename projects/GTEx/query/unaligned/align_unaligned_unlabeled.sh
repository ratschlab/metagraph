#!/bin/bash

set -e 

## source paths
# set metagraph
# set basedir
. ../../paths.sh

K=41
threads=6
mem=40000
pmem=$(($mem / $threads))

datadir=/cluster/work/grlab/projects/GTEx/rna2/results/alignments
query_dir=${basedir}/gtex/queries/unaligned_unlabeled

for fname in ${datadir}/*.all.unaligned.sorted.r1.fq.gz
do
    sample=$(basename $fname | cut -f 1 -d '.')
    outfile=${query_dir}/align_unaligned_unlabeled.${sample}.result.txt.gz
    logfile=${query_dir}/align_unaligned_unlabeled.${sample}.lsf.log
    if [ ! -f ${outfile} ]
    then
        graph=${basedir}/gtex/graphs/output_k${K}_trimmed_clean_graph_chunked/graph_merged_k${K}.dbg
        echo "/usr/bin/time -v $metagraph align -p ${threads} -i ${graph} ${fname} | gzip > $outfile" | bsub -M $mem -We 24:00 -n ${threads} -R "rusage[mem=${pmem}]" -R "span[hosts=1]" -J aln_gtex -oo $logfile
    fi
done
