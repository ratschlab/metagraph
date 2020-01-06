#!/bin/bash

set -e

## source paths
# set metagraph
# set basedir
. ../paths.sh

K=41

mem=30000
threads=6
pmem=$(($mem / $threads))
genome=/cluster/work/grlab/projects/metagenome/raw_data/ref_genomes/hg38/GRCh38.primary_assembly.genome.fa.bgz
outdir=${basedir}/gnomad/release_3.0
mkdir -p $outdir

for chr in $(seq 1 22) X Y
do 
    echo "/usr/bin/time -v ${metagraph} transform --to-fasta -v --parallel ${threads} -o ${outdir}/graph_gnomad_chr${chr}_k${K} ${outdir}/graph_gnomad_chr${chr}_k${K}.dbg  2>&1 | tee ${outdir}/gnomad_graph_fasta_chr${chr}.k${K}.log" | bsub -J gnmd_${chr} -oo ${outdir}/gnomad_graph_fasta_chr${chr}.k${K}.lsf.log -We 12:00 -n $threads -M $mem -R "rusage[mem=${pmem}]" -R "span[hosts=1]"
done
