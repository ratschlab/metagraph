#!/bin/bash
# Full comparison matrix. Run from the directory holding data/ and runs/.
# usage: run_all.sh BENCH_BINARY [PARALLEL_JOBS]
BENCH=$1; J=${2:-4}
S=$(dirname "$0")
SETS="illumina150 ont400_e3 ont400_e5 ont400_e10 ontlong_e3 ontlong_e5 ontlong_e10"
declare -A SEED=( [fullk]="--align-min-seed-length 31" [s19]="--align-min-seed-length 19"
                  [s19me0]="--align-min-seed-length 19 --align-min-exact-match 0" )
labeled() {
  for a in fullk s19 s19me0; do for e in default redone; do for r in $SETS; do
    echo "$e labeled_$1/$a $r -i data/graphs/$1.dbg -a data/graphs/$1.column.annodbg ${SEED[$a]}"
  done; done; done
}
jobs() {
  for g in canonical basic; do labeled $g; done
  for a in fullk s19; do for e in default redone; do for r in $SETS; do
    echo "$e unlabeled_canonical/$a $r -i data/graphs/canonical.dbg ${SEED[$a]}"
  done; done; done
  labeled primary
}
jobs | xargs -P $J -L 1 bash -c 'mkdir -p runs/$1; python3 '"$S"'/run_arm.py '"$BENCH"' $0 10 data/reads/$2.fa runs/$1/$0.$2 -- "${@:3}" || echo "FAILED $0 $1 $2"'
echo ALL_DONE
