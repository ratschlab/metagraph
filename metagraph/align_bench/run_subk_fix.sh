#!/bin/bash
# Rerun after align_refactor 51cf39976 (redone's sub-k reverse-complement pass
# on CanonicalDBG): redone, labeled, seed 19 on PRIMARY and BASIC, and on
# CANONICAL as a control. Output under runs3/.
# usage: run_subk_fix.sh BENCH_BINARY [PARALLEL_JOBS]
BENCH=$1; J=${2:-4}
S=$(dirname "$0")
SETS="illumina150 ont400_e3 ont400_e5 ont400_e10 ontlong_e3 ontlong_e5 ontlong_e10"
for g in primary basic canonical; do for r in $SETS; do
  echo "redone labeled_$g/s19 $r -i data/graphs/$g.dbg -a data/graphs/$g.column.annodbg --align-min-seed-length 19"
done; done | xargs -P $J -L 1 bash -c 'mkdir -p runs3/$1; python3 '"$S"'/run_arm.py '"$BENCH"' $0 10 data/reads/$2.fa runs3/$1/$0.$2 -- "${@:3}" || echo "FAILED $0 $1 $2"'
echo ALL_DONE
