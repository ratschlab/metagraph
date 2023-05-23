#!/bin/bash

set -e 

K=41
threads=8
mem=300000
pmem=$(($mem / $threads))
metagraph=/cluster/home/akahles/git/software/metagraph_current_dev/metagraph/build/metagraph

basedir=/cluster/work/grlab/projects/metagenome/data/gtex
query_dir=${basedir}/alignments_unaligned
datadir=/cluster/work/grlab/projects/metagenome/data/gtex/alignments_star_small/results/alignments
mkdir -p $query_dir

for fname in ${datadir}/*.all.unaligned.sorted.r1.fq.gz
do
    sample=$(basename $fname | cut -f 1 -d '.')
    outfile=${query_dir}/map_unaligned.${sample}.result.txt.gz
    logfile=${query_dir}/map_unaligned.${sample}.lsf.log
    if [ ! -f ${outfile} ]
    then
        echo "/usr/bin/time -v $metagraph query --discovery-fraction 0.0 --count-labels -v -p ${threads} -i ${basedir}/output_k${K}_trimmed_clean_graph_extended_chunked/graph_k${K}.dbg -a ${basedir}/output_k${K}_trimmed_extended.samples.brwt.annodbg ${fname} | gzip > $outfile" | bsub -M $mem -We 24:00 -n ${threads} -R "rusage[mem=${pmem}]" -R "span[hosts=1]" -J aln_gtex -oo $logfile
    fi
done
