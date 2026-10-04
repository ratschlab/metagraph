# Review request: attempts, deadlines, identity (8a98759a..ce949da5), and efficiency

**What I want from you:** two things on branch `gr/labeled-traversal`:
- **A. Correctness.** Review everything since your stage-3 review: commits `abc120fe`…`ce949da5`. Read their messages; each states what it claims.
- **B. Efficiency.** Critique and extend our efficiency plan, which is based on measurements against a deployed refseq33m server.

The normative text is `docs/SPEC-labeled-traversal-core.md` (§5, §6.8, §7.0, §7.3, §7.5.2, §10.3). The design records are `docs/DESIGN-traverse-graphlet.md` §19–§23. MGT v1 is unchanged.

## A. Correctness: what is new since 8a98759a

1. **Your stage-2 recheck and stage-3 review fixes** (`abc120fe`, `73a6bdde`):
   - floats priced at their widest canonical width;
   - every returned seed row charged before its comparison;
   - failed work stops state their largest charge;
   - trace-validation coordinates observed;
   - the zero-cache warm hang, the recorder's cache/cost split and the radius-0 false stop fixed;
   - extra labels reachable through a chain of switches;
   - continuation branch budgets with `reset_branches`;
   - per-walk `left_out`;
   - the spooled-receipt preflight;
   - the reworded "at least" and work statements.
2. **Stage 4, backend half** (`abc120fe`; DESIGN §22): `attempt_id` / `budget_id` / `locus_id`, plus:
   - a 409 on a reused id;
   - `POST /traverse/cancel` (200; 404; 404 with `tombstone: true`; 429);
   - `GET /traverse/attempt/{id}`;
   - the `usage` block;
   - an attempt bound the server enforces itself, including while the response is built and compressed (503 at the bound);
   - the stop when the client is gone (socket and TCP state polled at checkpoints);
   - `feature_level`. `schema_version` stays the request schema.
3. **Pass 5** (`87dbf983`, `c1605d39`, `ce949da5`; DESIGN §23; `feature_level` 3):
   - **`not_after_ms`:** 409 "expired" when the request starts too late.
   - **Per-graph identity on multi-graph servers:** CSV `manifest_path` / `index_ns`, checked at start-up against every file the pair loads; `capabilities?graph=&graph_path=`. An ambiguous `graph_path` is refused.
   - **The server-wide `GET /capabilities`.**
   - **Chunked deadlines:** a read the deadline may fall into is decoded in time-sized chunks in walk order; a read far from it stays one piece, with identical bytes.
   - **The delivery window:** each seed's text is written as soon as it is built; the walk stops early enough to build and compress what it walked; traversal routes use gzip level 1.
   - **Library:** `standalone_text`, and the J envelope's usage reduced to totals plus the seed's own entry.
   - **Tooling:** `index_manifest.py --server-csv`.

The evidence: unit filter 353, annotation 564, integration Traverse 75 and api 171 + 20, Python 620. Requests without the new features are byte-identical to 0a880475, decompressed:
- mini_refseq, 1,240 requests;
- UHGG, 1,649;
- a multi-graph list, 1,314;
- 1 ms chunks, 1,205.

On SRA, 8 differences were explained: cold-cache artefacts, `walk_until_ms`, and one reserve cut (see Known).

**Priorities for A**
1. **The lease guarantee, end to end.** Can an attempt outlive `not_after_ms + clock_skew_allowance_ms + bound_ms` (plus the stated uninterruptible overrun) on a server with the attempt contract? Places to look:
   - a request queued before the header read;
   - a slow upload;
   - a proxy holding the connection;
   - the cancel/finish race;
   - the 429 path;
   - restart, which gives a new `server_instance`.
2. **Chunked deadlines.** Find a request that keeps neither its bytes, when far from the deadline, nor its budget plus one chunk, when stopped. Also check: counters (`rows_requested`, `direct_reads`, `keys_mapped`) across chunk targets; the derivation window; the lookahead; annotate and trace modes; both arms.
3. **The delivery reserve.** Is the estimate (account / ratio + finished seeds' bytes, rates learnt per detail) sound? Does a 503 still happen where a 200 partial was possible? Is every statement true (`walk_until_ms`, `attempt_deadline`)?
4. **Identity.** Find a way for two different indexes on one multi-graph server to state the same `index_fp`, or for one pair to load a file its manifest does not cover. Check the `--server-csv` / `--digests` tool for mismatches it misses.
5. **The tombstone and 409 semantics** as a ledger would use them (DESIGN §22): is "release only on a finished state or a tombstoned cancel answer from the same server_instance" safe?
6. **MGT v1 and the library.** `standalone_text` must be byte-equal to `dump(from_response())`. Check the reduced-usage rule, the reserved names, and whether old saved documents stay canonical.

## B. Efficiency: what we measured, and what we plan

We ran 176 paced requests against the staging server `https://metagraph.ethz.ch:8080/api/refseq33m-experimental`. It runs 8a98759a on RefSeq 97: k=31, a basic graph of 627 G nodes; RowDiff<BRWT> with coordinates; 85,375 taxid columns; ~32.9 M accession headers. We then wrote a rerunnable benchmark (`scripts/traversal/bench_traverse.py`, `bench_library.py`; baseline in `~/.cache/metagraph-traverse-bench/staging-baseline`, window offset 2500). Findings:

- **Decode dominates:** 83–99.7% of server time on wide loci. Serialisation takes at most 6.8 ms; label-state work at most 2%; mapping coordinates through CoordToHeader 3–7%.
- **First touch:** 13–383 ms per row on wide loci; repeats are 9–83× faster; rows stay warm for minutes and go cold within about 20 minutes.
- **Row width:** a warm batched read costs 0.04 ms per row for one header and 10.6 ms for 16S (27.9k headers and 11.3k columns per k-mer). Cost follows columns more than coordinates. Column labels cost as much as header labels, because `TupleRowDiff::get_rows` builds every coordinate and then drops it.
- **The row-diff path is decoded again in every fetch call.** On warm 23S rows, 1 k-mer takes 25–86 ms, 3 k-mers 22–88 ms, 30 k-mers 72–138 ms. Duplicates are removed only within a call (row_diff.cpp:84,114), and LabelQuery caches only final rows. A walk on a branched locus fetches a few keys per level, so it pays the path again at every level: 27–41 ms per row against 4.6 ms batched.
- **The lookahead** read 1,047 rows to use 32. It ignored the remaining radius and ran outside the deadline (the deadline part is fixed by pass 5's chunking).
- **/resolve:** discovery decodes every row twice; header discovery counts into a `std::map` (0.87–1.83 s on 23S); `max_labels` bounds the output, not the work.
- **The Python library** is efficient: median call 1 ms or less, parsing at 10.7 MB/s. Hot spots on 1000-label graphlets:
  - merge partitions as arrays (as sets: claims 7.5×);
  - per-label scans in `routes` / `label_walks` (an index: 10.5–11×);
  - lookup by name (a dict: 370×);
  - rebuilt walk spellings (memoised: to_fasta 6×).

**Our plan, in order:**
1. **The cheap fixes, next:**
   - a per-request, byte-bounded cache of decoded row-diff path rows reused across fetch calls (estimated 5–19× per warm row);
   - the lookahead capped at the remaining radius;
   - `/resolve` without the double decode, and a flat sort instead of `std::map`;
   - a header-range binary search for coordinate mapping;
   - the library speed-ups;
   - raising the delivery reserve's starting estimates.
2. **Stage 3c:** decode only the selected labels' columns and header ranges (a prototype was exact in 966,495 checks; 25–110 µs per row against 180–980 µs on large binary BRWTs).
3. **Stage 3b:** budget-aware reads for the other formats, and `/resolve` budgets.
4. **Stage L:** budgets for the library's local operations.

**What I want from you on B:**
- Is the plan right, and in the right order?
- What is missing? For example:
  - BRWT-level batching or prefetching for first touch;
  - anchor density or the row-diff path length;
  - an mmap/page-cache strategy;
  - decoding rows straight into the permitted-label filter;
  - avoiding coordinate materialisation in constrain mode;
  - parallel decode inside one request;
  - server-side caching across requests, and its implications for determinism and for the budgets.
- For path reuse: how to keep budgets deterministic (charges batch-independent) while the actual decode shares paths? Today a returned row is charged its full dependency path ("deterministic logical work").
- For 3c: the crossover rule (about 64 columns) and the charge model (8 per dependency row + 1 per BRWT probe + 1 per selected coordinate).

Feel free to probe the staging server, but it is shared: strictly serial, at least 2 s between requests, a few dozen at most. The benchmark (`bench_traverse.py --window-offset 2500`) is paced.

## Known; do not re-report

- **The delivery reserve's starting estimates are cautious.** Ratio 20 for JSON and 40 for graphlets, until a server has measured a detail; measured values run 34.5–1,344. So the first large attempt per detail on a fresh server can be cut early (SRA: 42 MB partial instead of a complete result). It is stated, and we will raise the defaults.
- **Not chunked:**
  - `/resolve` (stage 3b);
  - a seed's k-mer mapping;
  - one row's decode (`max_uninterruptible_ms` null until 3c);
  - a cancel during a whole read far from the deadline, which is seen after the read.
- **Stated rate limits:** a read decoded whole overruns if its rows are more than 64× slower than any the request read before; the rest of a split read overruns if its rows are more than 4× slower than its chunk's.
- **No cap on `not_after_ms`.**
- **A name with a single pair ignores a wrong `graph_path`.**
- **Linux/GCC compilation is not verified here:** the explicit-instantiation access to the connection, and the TCP_INFO branch.
- **Open owner decisions:**
  - a short walk for a derived seed that runs out of time;
  - displaying the merge parent with the most labels;
  - `compare()` above 64 labels;
  - server-side maxima for memory and work (accepted, low priority).

## How to verify

```bash
cd metagraph/build && make -j metagraph unit_tests       # AppleClang: -DCMAKE_CXX_FLAGS=-Wno-error=unused-result
./unit_tests --gtest_filter='*Walker*:*LabelOracle*:*Resolve*:*Trie*:MiniRefSeq*:*Graphlet*:Stage2Review*:*CoordToHeader*:*HeaderIndex*:*Attempt*'   # 353
./unit_tests --gtest_filter='*RowDiff*:*BRWT*:*TupleRowDiff*:*IntRowDiff*:*Annotat*'           # 564
test_venv/bin/python ../integration_tests/main.py --test_filter='*Traverse*'                    # 75
cd .. && python3 -m unittest discover -s api/python/tests -p 'test_traverse_*.py'               # 620 (11 skipped)
python3 scripts/traversal/graphlet_fixtures.py --check && python3 scripts/traversal/graphlet_fixtures.py --from-cli build/metagraph --check
python3 scripts/traversal/make_codec_vectors.py --check && python3 docs/source/_graphlet_examples/run_examples.py
python3 scripts/traversal/bench_traverse.py --help-long                                          # the benchmark
```

## How to report

For each correctness finding, give:
- severity;
- whether it breaks the lease guarantee, byte identity or MGT v1;
- `file:line`;
- a concrete failure;
- your proposed fix.

Mark each one as confirmed, suspected or design disagreement.

For each efficiency suggestion, give:
- the expected gain, and on which workload;
- the evidence (measured or reasoned);
- the cost and risk;
- the interaction with determinism and with the budgets.

End with:
- a correctness verdict (sound, sound after named fixes, not sound);
- your ranking of the efficiency work.
