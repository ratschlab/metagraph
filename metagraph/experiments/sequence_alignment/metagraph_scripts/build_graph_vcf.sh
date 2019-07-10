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

GRAPH_FILE_NAME="hg19_hs37d5_${GRAPH_TYPE}_k_${K}_w_vcf_all"

echo "Building graph $GRAPH_FILE_NAME"
#echo "with reference file $DATA_FILE and vcf files $VCF_FILES" 
bsub -J "build_graph" -R"rusage[mem=200000]" -q normal.24h -n 1 ./metagengraph build -k $K -o $GRAPH_FILE_NAME --mem-cap-gb 150 --graph $GRAPH_TYPE --reference $DATA_FILE $DATA_FILE $VCF_FILES
