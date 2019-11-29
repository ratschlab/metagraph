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
TRUTH=$(paste <(paste -d'\n' <(grep -v "^[#@]" $OUT_PREFIX\1.aln | paste - - -) \
                             <(grep -v "^[#@]" $OUT_PREFIX\2.aln | paste - - -)) \
              <(grep -v "^@" $OUT_PREFIX.sam) |
              awk '{gsub(/-/, "", $5); print $9,$10,$4,$12,$5,$16}')

METAGRAPH_EXP=$PWD/run_experiments
paste <(echo "$TRUTH") \
      <($METAGRAPH_EXP evaluate_alignment score_cigar dna non-unit \
          <(echo "$TRUTH" | cut -d" " -f3,4,5,6)) > $OUT_PREFIX.truth.txt
