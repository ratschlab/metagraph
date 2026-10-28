# Results (2026-10-28, laptop, 4 concurrent single-threaded arms)

- `summary_bca751a37.tsv`: the full matrix, run on align_refactor `bca751a37` plus
  this bench commit. Redone has no work caps at this base. The default-engine
  numbers do not depend on the base. The rebased binary (`fe8744554`) gave
  identical non-latency output on 4/4 re-checked arms.
- `summary_caps_e9677f22b.tsv`: rerun on `e9677f22b`, where redone honours
  `max_nodes_per_seq_char` and `max_num_seeds_per_locus`. It covers all labeled
  redone arms at CLI defaults, redone with an unlimited node budget
  (`s19inf`, `s19me0inf`), and the default engine on CANONICAL at the same
  time as a latency reference.
- `fanout_s19_k31_canonical.tsv`: real k-mers that end with (`suffix`) and
  start with (`prefix`, forward resolution) each 19-mer read window, CANONICAL
  graph. Windows with at least one match.

Column definitions are in `metrics.py` and `fanout_summary.py`. Seeds:
- `fullk`: `--align-min-seed-length 31`;
- `s19`: 19;
- `s19me0`: 19 plus `--align-min-exact-match 0`;
- `*inf`: plus `--align-max-nodes-per-seq-char inf`.
