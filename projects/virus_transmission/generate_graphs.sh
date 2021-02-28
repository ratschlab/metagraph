#!/bin/bash

## Clean raw reads, construct and annotate joint graph, then extract differential unitigs
# USAGE
# ./generate_graphs.sh ${graph k} ${sample prefix}

METAGRAPH=~/metagraph/projects2014-metagenome/metagraph/build/metagraph
FLAGS="-p 10 --verbose"
K=$1
NAME=$2

for a in /cluster/work/grlab/projects/metagenome/raw_data/huber_virology/$NAME\_*.fastq.gz; do
    echo $a
    SNAME=$(basename $a .fastq.gz)
    $METAGRAPH build -k $K --mode canonical -o k$K/$SNAME --count-kmers $a $FLAGS
    $METAGRAPH clean \
        --prune-tips $((2*K)) \
        --prune-unitigs 0 \
        --fallback 3 \
        --to-fasta \
        --header "k$K/$SNAME.unitig." \
        -o k$K/$SNAME.contigs k$K/$SNAME.dbg $FLAGS
done

$METAGRAPH build -k $K --mode primary k$K/$NAME\_*.contigs.fasta.gz -o k$K/$NAME $FLAGS
$METAGRAPH assemble --unitigs $FLAGS k$K/$NAME.dbg -o k$K/$NAME.unitigs

$METAGRAPH annotate -i k$K/$NAME.dbg k$K/$NAME\_*.contigs.fasta.gz -o k$K/$NAME --anno-filename $FLAGS

for a in 0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45; do
    $METAGRAPH assemble -a k$K/$NAME.column.annodbg --unitigs $FLAGS k$K/$NAME.dbg \
        -o k$K/$NAME.diff.unitig.$a \
        --header "k$K/$NAME.diff.unitig.$a." \
        --enumerate \
        --label-mask-in k$K/$NAME\_D_0w_urine_raw.contigs.fasta.gz \
        --label-mask-in k$K/$NAME\_R_4-6w_urine_raw.contigs.fasta.gz \
        --label-mask-in k$K/$NAME\_R_52w_urine_raw.contigs.fasta.gz \
        --label-mask-out k$K/$NAME\_R_0w_urine_raw.contigs.fasta.gz \
        --label-other-fraction 0.0 \
        --label-mask-out-fraction $a
done

