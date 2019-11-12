#!/bin/bash

set -e

K=41
kmcdir=/cluster/work/grlab/projects/metagenome/data/gtex/output_k${K}
fnames_dir=${kmcdir}/fnames_kmc
outdir=/cluster/work/grlab/projects/metagenome/data/gtex/output_k${K}_ci1
mkdir -p $outdir

for uuid in SRR1098336 SRR1352481
do
    fnames_file=${fnames_dir}/${uuid}.txt
    $(pwd)/count_kmers_debug.sh $K $fnames_file $outdir $uuid
done
