# Aligner benchmark harness

Benchmark-only tooling (branch `bench/align-harness`). Nothing here is built by
default.

`align_bench` (CMake target, `EXCLUDE_FROM_ALL`) aligns each read with one
aligner object reused across reads, the way a language binding calls it. For
each read it records latency and enforces a watchdog. A crash is pinned to one
read, because output is flushed per read. It also has a seed fan-out mode.
Usage is in the source header.

The scripts are the same as on `bench/align-engine-compare`, where the data
recipe is documented (`align_bench/README.md` there):
- `run_arm.py` runs one arm, resuming past timeouts and crashes;
- `metrics.py` and `summarize.py` score the arms;
- `simulate_illumina.py`, `badread_to_fasta.py` and `simulate.sh` simulate the reads.

Example:

    python3 align_bench/run_arm.py build/align_bench default 10 reads.fa out/default.reads -- \
        -i graph.dbg -a graph.column.annodbg --align-min-seed-length 19
