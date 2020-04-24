#!/bin/bash

set -e

K=41

mem=250000
threads=6
pmem=$(($mem / $threads))
genome=/cluster/work/grlab/projects/metagenome/raw_data/ref_genomes/hg38/GRCh38.primary_assembly.genome.fa.bgz
basedir=/cluster/work/grlab/projects/metagenome/data
outdir=${basedir}/gnomad/release_3.0
mkdir -p $outdir
metagraph=/cluster/home/akahles/git/software/metagraph_current/metagraph/build/metagraph

for chr in $(seq 1 22) X Y
do 
    vcf=/cluster/work/grlab/projects/metagenome/raw_data/gnomad/release_3.0/gnomad.genomes.r3.0.sites.chr${chr}.vcf.bgz
    echo "/usr/bin/time -v ${metagraph} build --reference ${genome} -v -k ${K} --mem-cap-gb $((${mem} / 1200)) --parallel ${threads} -o ${outdir}/graph_gnomad_chr${chr}_k${K} ${vcf} 2>&1 | tee ${outdir}/gnomad_graph_chr${chr}.k${K}.log" | bsub -J gnmd_${chr} -oo ${outdir}/gnomad_graph_chr${chr}.k${K}.lsf.log -We 12:00 -n $threads -M $mem -R "rusage[mem=${pmem}]" -R "span[hosts=1]"
done
