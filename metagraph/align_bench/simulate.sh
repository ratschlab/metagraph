#!/bin/bash
# Read simulation recipe (run from data/). Seeds are fixed.
set -e
B=../tools/venv/bin/badread
S=../metagraph/metagraph/align_bench
python3 $S/simulate_illumina.py genomes 2000 150 0.005 100 1 > reads/illumina150.fa
NOADAPT='--start_adapter_seq "" --end_adapter_seq "" --start_adapter 0,0 --end_adapter 0,0 --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0'
# contig_map.tsv: contig, genome, contig length (see README)
run() { # name length identity quantity nreads minlen maxlen seed
  eval $B simulate --reference reads/allref.fa --quantity $4 --length $2 --identity $3 --seed $8 $NOADAPT > reads/$1.fq 2> reads/$1.badread.log
  python3 $S/badread_to_fasta.py reads/contig_map.tsv $6 $7 $5 $(( $5 / 20 )) $8 < reads/$1.fq > reads/$1.fa
}
run ont400_e3 400,40 97,99,1 600000 1000 300 500 11 &
run ont400_e5 400,40 95,98,1.5 600000 1000 300 500 12 &
run ont400_e10 400,40 90,94,2 600000 1000 300 500 13 &
run ontlong_e3 3000,1000 97,99,1 1500000 300 1000 5000 21 &
run ontlong_e5 3000,1000 95,98,1.5 1500000 300 1000 5000 22 &
run ontlong_e10 3000,1000 90,94,2 1500000 300 1000 5000 23 &
wait
grep -c '>' reads/*.fa
