# Review request: stage 2 of the resource contract, and the fixes to your implementation review

**What I want from you:** an adversarial review of two things on branch `gr/labeled-traversal` (read the commit
messages: each states what it claims):

- **Part A — stage 2** of the resource contract you approved (`docs/DESIGN-traverse-graphlet.md` v5.2, §14 and
  §14.1): commits `d4ee8f53`…`4d332729`.
- **Part B — the fixes** to your implementation review of `cb9eac2d` and to the 16 bugs a new real-index test suite
  found: commits `f10b1fc4`…`6e8891b2`. See the section *Part B* below.

## Part A: stage 2

The spec text is
`docs/SPEC-labeled-traversal-core.md` §5 (the request budgets), §6.8, §6.10, §7.0 (the `walk_domain` and
`memory_bound_soft` rows), §7.5 (MGT `Q` and the end reason `Y`, now written) and §10.3 (the capabilities).

Stage 2 claims three things. Please try to break each one:

1. **Atomicity.** Every head is planned, admitted and then committed.
   - The plan mutates nothing.
   - The admission reserves what the plan allocates. It also reserves what ending, merging and delivering those
     objects will cost later.
   - The commit replays the old mutation order and cannot fail.
   - A head that is not admitted is censored like a head beyond a cap: end reason `resource_limit` (`Y`), and its
     level does not count toward `complete_to_bp`.
2. **Bounds that hold, on the model.**
   - `bounds.max_memory_mb` charges a deterministic per-object model, never measured RSS. The model covers the
     heads' label states, segments, runs, events, recorded sets, fixed cache allotments and the output in the
     requested detail.
   - `bounds.max_work_units` charges a weighted sum of the work counters.
   - Both budgets, and the deadline, are checked before every head and at least every `W` = 65,536 units. That
     includes inside a derivation and in the seed phase.
   - What cannot be charged in advance is stated as `memory_bound_soft`, with the observed excess: the rows a
     level decodes, cache working-set excess and annotate dictionary growth.
3. **Byte identity when no budget is set.** With budgets unset, every result is byte-identical to the
   pre-stage-2 binary, with one intended exception: per-arm attribution of the work counters. This fixes your
   implementation-review finding 6: a left-arm re-minimisation was not stated as `greedy_losses`. The
   `ambiguous_split` fixture's `pair_evaluations` move between its arms.

   Lazy paths are part of this claim. `PathResult` keeps only its leaf, and chains are rebuilt from parent
   pointers by the JSON writer, `spell_path` and `write_continuation`. So a comb-shaped trie no longer finalises
   in quadratic memory.

The evidence I have:

| Check | Result |
|---|---|
| Golden CLI corpus (181 requests × 4 details) | 724/724 identical |
| Random requests (400 × 3 details) | 1,200/1,200 identical |
| Denial sweeps over all admission indices (~65k partial results checked) | 0 failures |
| Debug unit filter | 250/250 (257/257 after Part B, Release and Debug) |
| Integration suite | 38/38 (39/39 after Part B) |
| Python suite | 170 OK (519 OK, 11 skipped, after Part B) |
| Fixture checks | up to date |

An internal adversarial review (two lenses: atomicity, identity-and-bounds) found 7 issues. All are fixed and each
has a `Stage2Review.*` test in `tests/cli/test_graphlet_codec.cpp`.

### Priorities

1. **Atomicity** (`src/graph/traversal/walker.cpp`: `process_item`, `process_item_annotate`, `plan_outcomes`,
   `commit_item`, `commit_entries`, `merge_level`, `run_level`, `trip`, `finalize`).
   - Find a mutation of `SeedResult` or `ArmState` that happens before admission, or a commit step that can
     fail or allocate beyond its reservation.
   - Find a partial result whose runs, events, splits, paths, label summary and limitations disagree. The test
     hook `WalkerHooks::deny` injects a refusal at the N-th admission; `Graphlet.DenialLeavesAConsistentPrefix`
     uses it.
   - Look hardest at switches (the source run is ended in the commit, not the plan), splits (`begin_split` /
     `taken_run` clones), merges (reservations settled in `merge_level`), hairpins, trace and annotate mode, and
     the left arm.
2. **Byte identity.** Find an unbudgeted request whose output at `4d332729` differs from `aa5e3222` beyond the
   attribution change. Part B intentionally changes split clones' `route_bp`, annotate continuation labels and
   non-UTF-8 names; nothing else. Candidates:
   - event order at equal `at_bp`;
   - run ids and path ids after the lazy-path change;
   - counters;
   - continuations under `sequences: false`;
   - time stops. The deadline is now also checked inside a level, at depth > 0. That changes only results that
     hit the time budget, which were outside the determinism contract already. Is that acceptable?
3. **Honesty of the bounds.**
   - Is the memory model an upper bound on what the walker actually retains, with each object counted twice
     for vector growth and events three times for the stable sort? Measure it if you can, for example with
     jemalloc stats around a budgeted walk.
   - Does the reserve always let a partial result finalise and serialise?
   - Is the `memory_bound_soft` observed value an honest statement of the excess?
   - Can a hostile request make the walker hold much more than budget plus soft excess? Ideas: many labels,
     wide first levels, a huge derived seed set, annotate mode with `max_labels_per_node: unlimited`, or many
     seeds.
4. **The guarantee rule** (design v5.2 §14: a stated limitation is never omitted, and no outcome dimension is
   non-`complete` without the limitation that explains it).
   - A seed whose depth-0 state or seed phase exceeds its budget now fails: `walks: failed`, no arms, a
     `resource_stop`, and a **seed-level** `walk_domain` that names a lever. The levers are
     `lower_max_seed_labels`, `name_fewer_labels`, `lower_max_labels_per_node` and `shorten_seed`. Is that the
     right shape? The alternative was cutting the permitted set.
   - Find a budgeted result with a limitation missing, or with `outcome` disagreeing with its limitations.
5. **Determinism and scope.**
   - Budgets are per seed (a request of n seeds can hold n × budget), so a seed's stop does not depend on its
     position in the request.
   - Cache allotments are fixed fractions of the budget (label cache limit/4, lookahead limit/16, each capped
     at 64 MiB), so admission does not depend on `batch_kmers`.
   - `Graphlet.BudgetsAreDeterministic` checks batch size and seed order. Find a dependency it misses.
6. **Delivery.**
   - A memory stop depends on the detail by design: the output is charged in the requested detail. So the same
     budget stops detail full earlier than the graphlet; see the comb fixtures in
     `api/python/tests/data/traverse/documents/budgets/`.
   - Work is charged on the walk only, so a graphlet still rebuilds the detail-full result of the same walk.
     Consequence: neither work nor the deadline bounds the serialisation of a detail-full response. Only
     `max_memory_mb` bounds it. Is that acceptable, or is a deadline check inside the serialiser needed?
   - `DeliveryCostsBoundTheOutput` checks the per-object delivery costs against the real serialisers. Is the
     6× (JSON) / 10× (MGT) usefulness margin reasonable?

### Specific questions

- **No trip at the radius.** A budget does not trip at level start when the level is at the radius. Such a trip
  would mark the arm truncated with `complete_to_bp` equal to the radius, which breaks the `finalize` invariant.
  Is it correct that the level then completes even though its admission was not charged first?
- **Hook denial without a budget** is reported as a memory stop with limit `unlimited` (MGT `u`), so that it stays
  representable in `K` and `Q`. That only matters for C++ callers. Acceptable?
- **Annotate dictionary compaction.** After a budget stop or a refused admission, the annotate dictionary is
  compacted to the labels the result records. Caps without a budget keep the pre-stage-2 behaviour, so orphan
  labels remain possible there, for byte identity. Should that change in a later, versioned step?
- **Wholesale eviction.** `LabelQuery` and `LabelRecorder` now take an optional byte bound (set only under a
  memory budget) and still evict wholesale. Under a tight budget, could that make work explode, with the same
  rows fetched again level after level? The work units would charge it, but is the behaviour sensible?
- **Header-index cache** (`d4ee8f53`): it is keyed by address plus a fingerprint (column count, sizes and the
  corner headers). Can two distinct `CoordToHeader` objects at the same address share a fingerprint and still
  differ in ways that matter? A server loads one per process; the failure was seen only in tests.

## Part B: the fixes

### Your implementation review of `cb9eac2d`

| Finding | Fix | Commit | Tests |
|---|---|---|---|
| 1 UTF-8 replacement continues under another label | the replacement path is removed; such a seed is refused per seed (cause `unrepresentable_label_name`, label named by column/seq_id); the library refuses to build `next_request()` / `as_seed()` for names shared within the retrieval or containing U+FFFD (`UnverifiableLabelName`) | `c1617aaa`, `232d624c` | `Graphlet.NamesThatAreNotUtf8AreRefusedPerSeed`, integration `test_t37_*`, `TestFinding1*` |
| 2 store RAM accounting | `memory_bytes()` is a deep, deduplicating count (within 11 % of tracemalloc on every fixture); the store re-charges a model when its caches grow | `232d624c` | `TestFinding2*` |
| 3 annotate merges drop routes | claims, routes and label_walks follow the §5.1 union rule through any merge parent, route support kept apart from displayed support | `232d624c` | `TestFinding3*` |
| 4 comparisons beyond the cut | summaries derived from runs and intervals clipped to `depth_used` first | `232d624c` | `TestFinding4*` |
| 5 missing sequences | claims / walks / prefix_subset are `unknown` with a reason naming `output.sequences` | `232d624c` | `TestFinding5*` |
| 6 left-arm greedy losses | per-arm attribution of every work counter | `28c51c31` (stage 2) | `Walker.CountersAreAttributedToTheirArm` |
| 7 continuation skips identity | continue and replay check identity before sending and before storing; one-sided manifest → `index_mismatch`; none on either side → `identity.verified: false` | `232d624c` | `TestFinding7*` |
| 8 RANGES at `UINT64_MAX` | overflow-safe ordering and adjacency in encoder and decoder | `c1617aaa` | `GraphletCodec.RangesAtTheTopOfTheIdSpace` |
| requested: one time-stopped SeedResult both ways | serialised in detail full and graphlet from one object; the full result is rebuilt from the body | `c1617aaa` | `Graphlet.TimeStopSerialisesBothWaysFromOneResult` |

Run against the previous library, every finding's test class fails.

One place where we read a finding differently: for **finding 5**, a claims comparison without bases, matched by
label and interval, was arguably sound as `qualified`. We followed your proposal (`unknown`).

The library also gained a few fixes from writing its documentation:
- `Graphlet.splits()`;
- `exact` per arm (a clean arm no longer reads inexact because of the other arm);
- a claims filter that counts merge-entered claims as `merge_entered`;
- walks-mode comparison under a label filter drops walks the filter emptied;
- `python_requires >= 3.10`.

### The real-index suite

`api/python/tests/real/` drives the library against live servers on SRA (primary), UHGG and mini_refseq:
- 581 seed × strategy cells;
- an oracle built on `/resolve` alone;
- spellers that are independent of the library.

Latest run: 9,436 tests OK. 16 small retrievals are committed as offline CI fixtures
(`test_traverse_real_offline.py`, 270 tests).

It found 16 bugs, all fixed (`f10b1fc4`, `232d624c`). Two are in the walker and change unbudgeted output on
purpose:
- **Split clones and `route_bp`:** a run cloned at a split now inherits its source's `route_bp`. In the real cache
  that changed 976 runs.
- **Annotate continuation labels:** they now come from the tail's n − k + 1 k-mer nodes only. 10,031 continuations
  changed, each to a superset of the old list.

The rest are in the library:
- annotate claims and routes;
- compare cost: was quadratic, now linear;
- memory accounting;
- page-wise `graphlet_walks`;
- parse-error hygiene: huge integers, lone surrogates, echoed tokens;
- MCP tool argument and ceiling handling.

Please check that these two walker changes are right.

### Design questions where the code had to choose

The design text (v5.2 §5/§6) is now behind the code on these. Say which choice you would freeze:

1. **Annotate claims.** Each claim is a maximal end of a label's routes under the union rule. Its displayed
   support starts at the last merge its witness route enters through a non-first parent. `routes()` uses one
   witness per (label, maximal end), and the first carrying parent in `parents` order wins ties.
2. **`R.route_bp` against spec §7.1.** The walker stamps the *earliest* non-first-parent merge on the lineage's
   route. Spec §7.1 reads `[route_bp, to_bp)` as displayed support, which overstates support after two such
   merges. The library derives `evidence_from` from the partitions (the latest merge) and is unaffected. Example:
   `sra_hub__limit1`, left arm, run 1150 (stamp 61, merges at 113 and 104). Which definition should the spec use?
3. **`next_request()` refusal.** Any name containing U+FFFD is refused, with no opt-out. Should fixed servers
   advertise "names are verbatim" so the library can relax this?
4. **One-sided manifest.** It is now `index_mismatch` for continue and replay, which is stricter than before.
5. **Filter keys.** `graphlet_walks` still reports `filtered.route_only` for merge-entered walks; only the claims
   filter says `merge_entered`.
6. **MCP ceilings.** The `max_bytes` floor is 64. export, save and load take no `max_bytes`. The dry-run
   continuation returns the full normalised request under a 16 KB ceiling with `ceiling_bytes` stated.
7. **Local work budgets.** compare, routes and export have no work budget; an interrupted comparison would return
   `unknown`. Should §14.1 list this as a staged item with a stated limitation?

## Known limitations — do not re-report

- **Stage 3**, the budget-aware row-diff decoder (`memory_bound_soft` stays stated), and **stage 4**, the
  request-level ledger, are not built.
- `/resolve` and `/select` still pass label names through the JSON writer. This has the same root cause as
  finding 1 and is latent: it needs non-UTF-8 headers or file names. It is queued.
- The shared golden vectors (`codec_vectors.tsv`) have no RANGES case at `UINT64_MAX`. The C++ gtest and a
  Python parity test cover it.
- The Python suite ran on 3.14 only. `python_requires` is now `>=3.10`, and no 3.10 interpreter was available
  here.
- Record coordinates (design §18) are designed but not implemented.

## How to verify

```bash
cd metagraph/build && make -j metagraph unit_tests       # AppleClang: configure with -DCMAKE_CXX_FLAGS=-Wno-error=unused-result
./unit_tests --gtest_filter='*Walker*:*LabelOracle*:*Resolve*:*Trie*:MiniRefSeq*:*Graphlet*:Stage2Review*'   # 257 tests (cwd = build/)
test_venv/bin/python ../integration_tests/main.py --test_filter='*Traverse*'                    # 39 tests
cd .. && python3 -m unittest discover -s api/python/tests -p 'test_traverse_*.py'               # 519 tests (11 skipped)
python3 scripts/traversal/make_codec_vectors.py --check
python3 docs/source/_graphlet_examples/run_examples.py                                          # 37 examples
python3 scripts/traversal/graphlet_fixtures.py --check && python3 scripts/traversal/graphlet_fixtures.py --from-cli build/metagraph --check
# byte identity: build aa5e3222 in a worktree and compare full JSON and MGT on the same requests
```

A budgeted request, for example:

```json
"bounds": {"max_extension_bp": 240, "max_memory_mb": 1}
"bounds": {"max_extension_bp": 100, "max_work_units": 1500}
```

## How to report

For each finding, give:
- severity (blocker / major / minor);
- whether it breaks byte identity of unbudgeted results, or forces a change to MGT v1 (frozen) (yes / no);
- `file:line`;
- a concrete failure (request → inconsistent partial result, a bound exceeded beyond the stated soft excess, a
  missing limitation, a crash);
- the fix you propose.

Separate *confirmed* (reproduced), *suspected* (say which experiment settles it) and *design disagreement*. End
with an explicit verdict for each part: sound, sound after named fixes, or not sound; and your choice for each
design question.
