#!/bin/bash
# Rerun after the redone engine honours max_nodes_per_seq_char and
# max_num_seeds_per_locus (align_refactor ac3e22615): all labeled redone arms
# at CLI defaults (budget 5.0, seed cap 1000), redone with the node budget
# unlimited on CANONICAL, and the default engine on CANONICAL as a same-time
# latency reference. Output under runs2/.
# usage: run_caps.sh BENCH_BINARY [PARALLEL_JOBS]
BENCH=$1; J=${2:-4}
S=$(dirname "$0")
SETS="illumina150 ont400_e3 ont400_e5 ont400_e10 ontlong_e3 ontlong_e5 ontlong_e10"
declare -A SEED=( [fullk]="--align-min-seed-length 31" [s19]="--align-min-seed-length 19"
                  [s19me0]="--align-min-seed-length 19 --align-min-exact-match 0"
                  [s19inf]="--align-min-seed-length 19 --align-max-nodes-per-seq-char inf"
                  [s19me0inf]="--align-min-seed-length 19 --align-min-exact-match 0 --align-max-nodes-per-seq-char inf" )
arm() { echo "$1 labeled_$2/$3 $4 -i data/graphs/$2.dbg -a data/graphs/$2.column.annodbg ${SEED[$3]}"; }
jobs() {
  for a in fullk s19 s19me0; do for r in $SETS; do arm redone canonical $a $r; arm default canonical $a $r; done; done
  for a in s19inf s19me0inf; do for r in $SETS; do arm redone canonical $a $r; done; done
  for g in basic primary; do for a in fullk s19 s19me0; do for r in $SETS; do arm redone $g $a $r; done; done; done
}
jobs | xargs -P $J -L 1 bash -c 'mkdir -p runs2/$1; python3 '"$S"'/run_arm.py '"$BENCH"' $0 10 data/reads/$2.fa runs2/$1/$0.$2 -- "${@:3}" || echo "FAILED $0 $1 $2"'
echo ALL_DONE
