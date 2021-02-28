#!/bin/bash

set -e

K=41

#mem=10000
#mem=50000
mem=100000
#mem=150000
threads=1
pmem=$(($mem / $threads))
basedir=/cluster/work/grlab/projects/metagenome/data
#graphdir=${basedir}/gtex/parameterstudy/output_k${K}_cleaned
#outdir=${basedir}/gtex/parameterstudy/output_k${K}_cleaned
graphdir=${basedir}/gtex/parameterstudy_largest/output_k${K}_cleaned
outdir=${basedir}/gtex/parameterstudy_largest/output_k${K}_cleaned
#graphdir=${basedir}/gtex/parameterstudy_largest_trimmed/output_k${K}_cleaned
#outdir=${basedir}/gtex/parameterstudy_largest_trimmed/output_k${K}_cleaned
mkdir -p $outdir
metagraph=/cluster/home/akahles/git/software/metagraph_current2/metagraph/build/metagraph

for inseq in ${graphdir}/*.clean_s1.fasta.gz
do
    outgraph=${inseq%.clean_s1.fasta.gz}.clean_s1.dbg
    outseq=${inseq%.clean_s1.fasta.gz}.clean_s2.fasta.gz
    logfname=${inseq%.clean_s1.fasta.gz}.clean_s2.lsf.log
    if [ -f ${outseq} ]
    then
        echo "$inseq complete"
        continue
    fi

    echo "/usr/bin/time -v $metagraph build -k $K -o ${outgraph%.dbg} --mode canonical -v ${inseq} && $metagraph assemble --unitigs -o ${outseq%.fasta.gz} $outgraph" | bsub -J ms_clean -o ${logfname} -We 8:00 -n 1 -M $mem -R "rusage[mem=${pmem}]" -R "span[hosts=1]"
done
