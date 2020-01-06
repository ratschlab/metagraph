#!/bin/bash

set -e

anno=/cluster/work/grlab/projects/metagenome/raw_data/ref_genomes/gencode_v30/gencode.v30.annotation.gtf
outdir=/cluster/work/grlab/projects/metagenome/raw_data/ref_genomes/gencode_v30/

. "/cluster/home/akahles/anaconda3/etc/profile.d/conda.sh"
conda activate spladder
cd /cluster/home/akahles/git/software/spladder_p3
python -m spladder.spladder build -a $anno -o $outdir -b $anno --no-extract-ase --no-insert-ni --no-insert-es --no-insert-ir --filter-overlap-exons -v  
