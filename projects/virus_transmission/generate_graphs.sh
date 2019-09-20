#!/bin/bash

METAGRAPH=~/metagraph/projects2014-metagenome/metagraph/build/metagraph
FLAGS="-p 10 --verbose"
K=$1
NAME=$2

for a in /cluster/work/grlab/projects/metagenome/raw_data/huber_virology/$NAME\_*.fastq.gz; do
    echo $a
    SNAME=$(basename $a .fastq.gz)
    $METAGRAPH build -k $K --canonical -o k$K/$SNAME --count-kmers $a $FLAGS
    $METAGRAPH clean \
        --prune-tips $((2*K)) \
        --prune-unitigs 0 \
        --fallback 2 \
        --to-fasta \
        -o k$K/$SNAME.contigs k$K/$SNAME.dbg $FLAGS
done

$METAGRAPH build -k $K k$K/$NAME\_*.contigs.fasta.gz -o k$K/$NAME $FLAGS
$METAGRAPH assemble --unitigs $FLAGS k$K/$NAME.dbg -o k$K/$NAME.unitigs

$METAGRAPH annotate -i k$K/$NAME.dbg k$K/$NAME\_*.contigs.fasta.gz -o k$K/$NAME --anno-filename $FLAGS

for a in 0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 1.00; do
    $METAGRAPH assemble -a k$K/$NAME.column.annodbg --unitigs $FLAGS k$K/$NAME.dbg \
        -o k$K/$NAME.diff.$a \
        --label-mask-in k$K/$NAME\_D_0w_urine_raw.contigs.fasta.gz \
        --label-mask-in k$K/$NAME\_R_4-6w_urine_raw.contigs.fasta.gz \
        --label-mask-in k$K/$NAME\_R_52w_urine_raw.contigs.fasta.gz \
        --label-mask-out k$K/$NAME\_R_0w_urine_raw.contigs.fasta.gz \
        --label-mask-out-fraction $a
done

