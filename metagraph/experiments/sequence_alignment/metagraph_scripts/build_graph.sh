#!/bin/bash

if [ $# -ne 2 ]; then
    echo "usage: ./build_graph <k> <graph>"
    echo "graph could be: succinct, hash, hashstr, bitmap"
    echo "this script should be in same directory as metagraph binary"
    exit
fi

K=$1
GRAPH_TYPE=$2
DATA_FILE="/cluster/work/grlab/projects/metagenome/raw_data/ref_genomes/hg19_hs37d5/genome.fa"

VCF_DIR="/cluster/work/grlab/projects/metagenome/data/alignment/1000G/vcf_per_chr"
VCF_FILES=$(find $VCF_DIR -name "*.vcf.gz")

GRAPH_FILE_NAME="hg19_hs37d5_${GRAPH_TYPE}_k_${K}"

echo "Building graph $GRAPH_FILE_NAME"
bsub -J "build_graph" -R"rusage[mem=80000]" -q normal.24h -n 1 ./metagengraph build -k $K -o $GRAPH_FILE_NAME --mem-cap-gb 80 --graph $GRAPH_TYPE $DATA_FILE
