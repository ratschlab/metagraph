#!/bin/bash

# 150bp Illumina read simulation workflow

COVERAGE=$1
REFERENCE=$2
OUT_PREFIX=$3

READ_LEN=150
MEAN_FRAGMENT_LENGTH=300
SD_FRAGMENT_LENGTH=10
SEED=42

#generate
art_illumina \
    -ss HSXt \
    -sam \
    -i $REFERENCE \
    -l $READ_LEN \
    -m $MEAN_FRAGMENT_LENGTH \
    -s $SD_FRAGMENT_LENGTH \
    -f $COVERAGE \
    -rs $SEED \
    -p \
    -o $OUT_PREFIX

echo "Interleaving paired read files"
paste -d'\n' <(cat $OUT_PREFIX\1.fq | paste - - - -) \
             <(cat $OUT_PREFIX\2.fq | paste - - - -) |
    awk '{print $1; print $2; print $3; print $4}' > $OUT_PREFIX.fq

#report
echo "Extracting ground truth data"
paste <(paste -d'\n' <(grep -v "^[#@]" $OUT_PREFIX\1.aln | paste - - -) \
                     <(grep -v "^[#@]" $OUT_PREFIX\2.aln | paste - - -)) \
      <(grep -v "^@" $OUT_PREFIX.sam) |
    awk '{print $10,$4,$12,$5}' > $OUT_PREFIX.truth.txt
