# Review request: levels 4 and 5 (e7d99e0c..dcc0cebd)

**What I want from you:** three things on branch `gr/labeled-traversal`:
- **A. Re-check the blockers of your last review.** Findings 1–4 broke the lease-release and identity guarantees. Are they closed?
- **B. Correctness of everything new since e7d99e0c.** Four commits:
  - level 4: `c18280d2` (C++), `f667d775` (Python);
  - level 5: `7cdde4c2` (C++), `dcc0cebd` (Python).
  Read their messages; each states what it claims.
- **C. Efficiency.** Critique the measurements and the plan for the next stages.

The normative text is `docs/SPEC-labeled-traversal-core.md`. The design record is `docs/DESIGN-traverse-graphlet.md`; §24 covers your review and level 5. MGT v1 is unchanged.

Staging (`https://metagraph.ethz.ch:8080/api/refseq33m-experimental`) runs **level 4** (`f667d775`), with `max_time_ms` 20000 and `max_query_bp` 10000. Level 5 is not deployed yet.

## A. Your findings, and what we did

| # | Finding | Fix (level 5) | Regression test |
|---|---|---|---|
| 1 | A tombstone expires; a half-uploaded request then runs | See the next list | Your half-upload probe, with and without `not_after_ms` in the cancel |
| 2 | Symlinked main files hide different sidecars | One loader inventory per listed spelling: sidecars resolved next to the spelling. Bundles deduplicated by the complete inventory's realpaths; every bundle validated against its manifest | Your A/B symlink case: start-up is refused |
| 3 | `.bloom` missing from manifests | The inventory covers every file the loaders open. `metagraph traverse --index-inventory` and `index_manifest.py --inventory` share it, with a cross-check test | Your bloom case: refused with the old manifest; the generated manifest lists `.bloom` |
| 4 | Negative retention wraps | Retention seconds, count and tombstone cap parsed as bounded integers; anything else refuses to start. Retention 0 answers 429 `no_suppression` | Your probe plus eight values, as start-up refusals |
| 5 | Budgeted lookahead chunks globally sorted rows | Chunks are chosen in first-occurrence (walk) order and sorted within each chunk. Charges are order-independent | Your pacing probe, as decoder-work assertions |
| 6 | Delivery checks not every 64 KiB | Writes and assembly copy in bounded pieces with a check between them. A single large token's preparation is measured into `observed_max_uninterruptible_ms` | Your 16 MiB probe: one check per piece |
| 7 | Malformed costs raise `TypeError` | Validated before any membership test: field-qualified `ValueError`; `bad_argument` from the tools | Your probe plus a fuzz over malformed shapes |
| — | "never started" | Now "cannot start subsequently" | — |

**How finding 1 was fixed:**
- `POST /traverse/cancel` takes `not_after_ms`.
- Every tombstone answer states:
  - `suppressed_until_ms`;
  - `not_after_ms`;
  - `covers_admission`;
  - `covers_admission_reason`: `no_not_after_ms` or `beyond_tombstone_max`.
- How long a tombstone is held:
  - until `not_after_ms` + `clock_skew_allowance_ms`;
  - never less than `retention_s`;
  - capped by `--traverse-attempt-tombstone-max-s` (default 86400);
  - kept on both the steady clock and the wall clock.
- Repeats and refused copies only extend a tombstone. A dispatch-time 409 for a tombstoned id extends it to that request's own `not_after_ms`.
- New request field `expect_server_instance`, answered with 409 `instance_mismatch`. It protects against a replay after a restart.
- `attempts.release_rule` in the capabilities is normative. Release early only when **all** hold:
  - the answer has `covers_admission: true`;
  - it comes from the same `server_instance`;
  - the attempt was sent with exactly that `not_after_ms` **and** with `expect_server_instance`.
  Otherwise, release only on a finished state, or after `not_after_ms` + skew + `bound_ms`.

**Priorities for A**
1. **The release rule, end to end.** Can an attempt run after an answer with `covers_admission: true`? Try:
   - slow and half uploads;
   - queued requests;
   - repeated cancels with different `not_after_ms`;
   - the count limit and 429;
   - retention 0;
   - the cap;
   - wall-clock steps;
   - a restart with and without `expect_server_instance`;
   - proxy replays.
   Is every sentence in `attempts` true?
2. **Identity.** Can two different indexes on one multi-graph server state the same `index_fp` / `index_meta_fp`? Can a pair load a file its manifest does not cover? Read the loaders yourself (graph, annotation, CoordToHeader, header index, row-diff files) and compare them with the inventory. Check hard links, relative and absolute spellings, and a sidecar present for one spelling only.

## B. What is new besides your findings

**Level 4 (`c18280d2`, `f667d775`)**
- **Row-diff path cache** (`row_diff_cache.hpp`, `--traverse-path-cache-mb`, default 128 MiB per request):
  - a later fetch stops at a cached row instead of decoding to the anchor;
  - budget-aware decode stops only at rows that carry their whole path's aggregates, so `RowCost` and admission stay cache-independent.
- **Lookahead:** capped at the remaining radius, and paced by the walk deadline at every depth.
- **/resolve:** one decode for discovery, and flat label counts (same ranking).
- **Coordinate mapping:** an adaptive run mapper through `CoordToHeader::sequence_range`.
- **Delivery reserve:** starting estimates 30/50, `stop_ms` = chunk_target + 950.
- **Server maxima:** `--traverse-max-memory-mb` and `--traverse-max-work-units` clamp budgets (off by default).
- **Library:** speed-ups (sets, a runs-by-label index, a name dict, memoised spellings) and **stage L** local budgets (DESIGN §21 L1–L9, off by default).

**Level 5, besides your findings (`7cdde4c2`, `dcc0cebd`)**
- **The first-touch regression of level 4 on staging.** A 1 s budget took 9.9 s on first touch.
  - 9.3 s of that was the 33M-header reverse index, which the first header-named request built inside its own budget. The server now builds it while the index loads.
  - Separately, the path cache copied every row of every path. On wide tuple rows that made a one-row-per-call walk take 30 s instead of 3.6 s (synthetic index with 8,000 labels and 24,000 coordinates per row).
  - Fix: a retention rule (requested rows, 8 successors, a checkpoint every 16 rows, narrow rows) plus flat storage of tuple rows. Copying on that walk fell from 78,472 rows / 50 GB to 345 rows / 0.22 GB.
- **Timing (only with `output.timing`):**
  - `seed_phase_ms`, `seed_fetch_ms`, `label_resolve_ms`;
  - `path_cache {...}`;
  - `deadline {longest_piece {ms, kind, rows, coordinates}, stopped_by, stop_after_deadline_ms}`.
- **The library's lazy indexes:** small single-label calls are back near ce949da5.
- **Benchmark:** `bench_traverse.py` compares reach and consumed work for deadline-limited requests.

**Priorities for B**
1. **Path-cache determinism.**
   - Is a returned row's logical charge exactly its standalone (uncached) charge, whatever the cache held?
   - Are outputs and budgeted stops identical across cache capacities (0, tiny, default), `batch_kmers`, chunk targets, and cold or warm state?
   - Does the retention rule ever change a result?
   - Is the cache's share of a memory budget sound?
2. **The deadline record.** Is `longest_piece` honest about what it measures? Is any uninterruptible piece missing from its kinds?
3. **The level-4 items.**
   - Does the lookahead cap change only the stated counters?
   - Does /resolve keep its count and tie ordering?
   - Is the run mapper exact?
   - Are the maxima clamped and echoed correctly, including `requested: "unlimited"`?
   - Is the reserve's calibration text true?
4. **Stage L.**
   - Are its units deterministic across cold and warm state?
   - Are there no partial exports, `parsed: false` on a stopped parse, `compare` answering "unknown", and cursors that rebuild a list exactly?
   - Is `ops.compare_cost` an honest estimate?
5. **The header index built at load.** It costs start-up time and about 1–2 GB on every server with coordinates. Is that the right trade-off, or should it be a flag?

**Evidence**
- **Tests:** unit 383, annotation 571, integration Traverse 86 and api 171 + 20, Python 738 (11 skipped), fixtures, codec vectors, 40 docs examples.
- **GCC 13.5:** `-O3 -Werror` on all 112 `src/` files. Staging builds with GCC 13.
- **Byte identity against f667d775** (decompressed): 2,470 unbudgeted requests, 588 in multi-graph mode and 1,968 budgeted ones. With timing off they are identical; with timing on, only the new timing keys appear.
- **Byte identity, level 4 against ce949da5:** 4,264 unbudgeted and 1,176 budgeted requests differ only in the stated changes.
- **Staging benchmarks** (`~/.cache/metagraph-traverse-bench/`):
  - `staging-baseline` (8a98759a), `staging-after-pass5` (e7d99e0c), `staging-level4` (f667d775);
  - comparisons: `staging-after-pass5/compare_vs_baseline.md`, `staging-level4/compare_vs_pass5.md`.
- **Wide-index reproduction:** `~/.cache/metagraph-traverse-bench/level5-evidence/` (`final_table.txt` and the generator scripts).

## C. Efficiency: measurements and plan

**Level 4 against level 3 on staging (refseq33m, same panel, K=2500)**

| What | Level 3 | Level 4 |
|---|---|---|
| 23S named walk, warm | 15.1 ms/row | 1.5 ms/row |
| Walks cut off by their budget | — | reach about 2× more bases |
| Lookahead suite | — | 4–12× faster |
| /resolve header discovery on 23S | 784 ms | 171 ms |
| Async e2e: cold /resolve of 573 bp | 10.3 s | 1.19 s |

At level 4, first-touch reads regressed (fixed in level 5, as above). Explicit /resolve and single calls were 5–40% slower. Wide rows still cost several ms per row warm.

**What we learned about the deployment**
- The container is capped at 128 GiB with `--mmap --madv-random`.
- The index is 648 GiB: graph 338 GiB plus annotation 310 GiB.
- The page cache counts against the container's limit, so most of the index cannot stay warm. That is very likely the main source of first-touch cost. The host has 1.1 TB of RAM.
- `--madv-random` turns off read-ahead.

**Plan for the next stages** (`docs/PLAN-traverse-next-stages.md`, approved by the owner; your comments on 3c and 3b were folded in):

| Level | Stage |
|---|---|
| 6 | Record coordinates (DESIGN §18) |
| 7 | 3c-core: selected-label decoding |
| 8 | 3c-ii: checkpoints inside reads, so `max_uninterruptible_ms` can become finite |
| 9 | 3b-1: /resolve budgets |
| 10 | 3b-2: the remaining annotation formats |

**What I want from you on C**
- Is the plan right?
  - Is the order right?
  - Is the 3c charge model right (decision X-R3: 8 per dependency row; BRWT probes counted only below a nonzero root rank; stored selected tuples; stored coordinates in the canonical window, cancelled ones included; a bounded deep term; lazy shifting)?
  - Is the 64-column crossover right?
- Ops levers that likely matter more than code: the container memory limit, and `--madv-random` against read-ahead. How would you measure them without disturbing a shared host?
- Anything else in the plan that looks wrong.

## Known; do not re-report
- A same-size replacement of a manifest file is accepted at start-up with an unchanged `index_fp`. Only `index_manifest.py --verify` catches it; this is stated.
- `unit_tests` is not GCC `-Werror` clean in 3 test files; this was already so at the base. `src/` is clean.
- Stage L accounts are conservative, up to about 12× the peak on some operations. Budgeted overhead is 0–16% on millisecond calls. The time calibration was done on a loaded machine.
- /resolve costs +24% / +48% CPU on queries that repeat wide rows beyond 256 MiB (synthetic index), and up to +14% on 3 kb header discovery.
- Held finished attempts are bounded by `tombstone_max_s`, not by `retention_count`.
- The attribution of the staging 9.3 s to the header index is not yet confirmed on staging, because level 5 is not deployed.
- The server ignores SIGTERM; `docker stop` waits its full timeout.
- Owner decisions since your last review:
  - a derived seed that runs out of time returns a partial result if anything was derived;
  - the displayed walk at merges follows the parent carried by the most labels;
  - `compare()` above 64 labels is explained in the MCP layer.
  These land with level 6.

## How to verify

```bash
cd metagraph/build && make -j metagraph unit_tests       # AppleClang: -DCMAKE_CXX_FLAGS=-Wno-error=unused-result
./unit_tests --gtest_filter='*Walker*:*LabelOracle*:*Resolve*:*Trie*:MiniRefSeq*:*Graphlet*:Stage2Review*:*CoordToHeader*:*HeaderIndex*:*Attempt*:*PathCache*'   # 383
./unit_tests --gtest_filter='*RowDiff*:*BRWT*:*TupleRowDiff*:*IntRowDiff*:*Annotat*'            # 571
test_venv/bin/python ../integration_tests/main.py --test_filter='*Traverse*' </dev/null        # 86
cd .. && python3 -m unittest discover -s api/python/tests -p 'test_traverse_*.py'                # 738 (11 skipped)
python3 scripts/traversal/graphlet_fixtures.py --check && python3 scripts/traversal/graphlet_fixtures.py --from-cli build/metagraph --check
python3 scripts/traversal/make_codec_vectors.py --check && python3 docs/source/_graphlet_examples/run_examples.py
./build/metagraph traverse --index-inventory --json -i G.dbg -a A.annodbg                       # the loader inventory
```

Probe staging only gently: it is shared. Requests strictly serial, at least 2 s apart, a few dozen at most. Never between 01:00 and the end of an announced overnight run.

## How to report

For each correctness finding, give:
- severity;
- whether it breaks the lease guarantee, index identity, byte identity or MGT v1;
- `file:line`;
- a concrete failure;
- your proposed fix.

Mark each one as confirmed, suspected or design disagreement.

For each efficiency suggestion, give:
- the expected gain and on which workload;
- the evidence;
- the cost and risk;
- the interaction with determinism and with the budgets.

End with:
- a correctness verdict (sound, sound after named fixes, not sound);
- your ranking of the plan.
