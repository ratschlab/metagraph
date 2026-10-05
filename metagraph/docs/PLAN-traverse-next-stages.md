# Next stages: program plan — record coordinates, stage 3c, stage 3b

**Date:** 2026-10-04. **Base:** f667d775 (feature_level 4, deployed on staging) plus the level-5 fix batch
(R1–R10, P7–P10), which is still running and must be committed first.
**Status:** plan only. No repository file was changed.

**Sources**
- The three stage plans:
  - `plans/coordinates/plan.md`
  - `plans/stage3c/plan.md`
  - `plans/stage3b/plan.md`
- Their critiques, all checked against the code at f667d775.
- DESIGN §18 and §21 (your decisions C1–C12, R2–R7, B1–B16).
- The external review of 5801aea1, part B.
- Staging measurements at level 4.
- The fix batch task list (`wf/review-pass5-fixes.js`).

**How to read this.**
- The stage plans remain the detailed reference: file lists, function names, test names.
- This document overrides them wherever it records a revision.
- Each stage section ends with a table that maps every confirmed critique issue to its resolution.

---

## 0. At a glance

| # | Stage | Level | Effort (workflow h) | Elapsed (h) | Starts when |
|---|---|---|---|---|---|
| 0 | Fix batch R1–R10, P7–P10 (running now) | **5** | — | — | — |
| 1 | Record coordinates | 6 | 38–48 | 29–36 | fix batch committed as level 5 |
| 2 | Stage 3c-core: selected-label decoding | 7 | 45–55 | 43–52 | coordinates committed |
| 3 | Stage 3c-ii: checkpoints inside reads | 8 | 15–20 | 15–20 | 3c-core committed |
| 4 | Stage 3b-1: /resolve budgets and level rows | 9 | 26–35 | 26–35 | 3c-ii committed |
| 5 | Stage 3b-2: the remaining annotation formats | 10 | 25–34 | 18–26 | 3b-1's /resolve part committed (its library part can start earlier in a worktree) |
| | External review fix passes (one per level) | | 32–51 | 32–51 | after each review returns |
| | **Total** | | **181–243** | **about 165–220** | |

**The one urgent decision (P1).**
- The fix batch must commit as **level 5**.
- Level 4 is on staging: release `refseq97-f667d775…` reports `feature_level 4`.
- The batch's own task text says the fixes "go into the same level 4" because level 4 is "NOT deployed anywhere". That premise is false.
- The working tree still has `kTraverseFeatureLevel = 4`.
- The batch changes capabilities: tombstone and release fields, uninterruptible observations, calibration text. A second level-4 capability set would break the monotonic rule (§21).
- Recommendation: the main session sets 5 in the commit step. That means `traverse.hpp`, the SPEC §10.3 mapping line and the DESIGN §24 text. A quick capabilities check follows. The running workflow need not be stopped.

**Recommended order:** coordinates → 3c-core → 3c-ii → 3b-1 (/resolve budgets first) → 3b-2. This is your §21 order. Two refinements:
- 3c's in-read checkpoints become their own level.
- /resolve budgets come first inside 3b.

---

## 1. Order and rationale

### 1.1 Recommendation: keep coordinates → 3c → 3b

**Why coordinates first**
- It is the smallest and least risky stage. It is opt-in, additive, and adds no annotation reads or work units.
- It delivers a contract you have fully decided (C1–C12). The library and the agents are waiting for it (C7, C8).
- The refseq33m wide-row cost does not grow: coordinates reuse rows already fetched.
- It gives 3c an end-to-end positional oracle. With coordinates on, a wrong window shift in 3c shows up as wrong reported positions, not only as different trace decisions.
- Its byte-identity set with coordinates on becomes part of the 3c and 3b gates.

**Why 3c before 3b**
- The external review ranks selected-column decoding ahead of "/resolve budgets and remaining format coverage".
- 3b reuses three things from 3c:
  - `DecodeControl` and `INTERRUPTED`, for /resolve windows and the new readers;
  - the selected decoder, for /resolve explicit labels;
  - the price definition (R3).

  Building /resolve budgets first would mean retrofitting all three later.
- 3c is the main performance win on staging. Walks with one header label are time-bound at 5 s with 1.5–9.4 ms per row (level 4).

**Why split 3c into core and checkpoints (P3)**
- 3c-core can be deployed and measured on staging without waiting for the checkpoint work.
- The checkpoint work touches about ten decode call sites. Its routing choice (D8) needs its own A/B/A data.
- 3b's new readers can then poll `DecodeControl` from their first line.
- The cost is one more feature level. Levels are cheap integers.

**Why /resolve budgets first inside 3b**
- The external review advises it.
- /resolve budgets and walker change (c) serve the production row-diff indexes (refseq33m, UHGG, SRA).
- The remaining formats have no known production user.

### 1.2 Alternatives considered

| Order | For | Against | Verdict |
|---|---|---|---|
| A. coordinates → 3c → 3b (yours) | Smallest first; contract delivered early; positional oracle guards 3c; 3b reuses 3c | refseq33m walks stay slow until 3c lands, about 30–36 elapsed hours later | **Recommended** |
| B. 3c → coordinates → 3b | Faster refseq33m walks about 35 h sooner | Delays coordinates by about 60 h; the 3c gate lacks coordinate requests; saves no effort; riskiest stage first | Not recommended |
| C. coordinates → /resolve budgets → 3c → rest of 3b | Closes the unbudgeted /resolve exposure sooner (23S header discovery 5.7 ms per k-mer; 16S discovery 28.9 s cold) | The /resolve engine is built without DecodeControl or selected reads, and 3c then retrofits both; goes against the review's ranking | Not recommended. Mitigate /resolve now through launch flags (ops note O1) |

---

## 2. Stage plans

Invariants that every stage keeps:
- MGT v1 frozen.
- Your guarantee rule: never weaken silently, state limitations, be conservative.
- Opt-out requests byte-identical (decompressed) to the previous level, apart from the `feature_level` digit and timing.
- Budget charges independent of batch, cache, chunking and timing.
- The attempts' lease bound unaffected.
- `feature_level` monotonic.
- -Werror on AppleClang and GCC (Docker gcc:12).
- Python parity for every new field.

### 2.1 Record coordinates (level 6)

**Goal**
- With `output.coordinates: true`, every seed result reports where each run's sequence lies.
- Positions are within the record for header labels and global for column labels.
- Every cut and every incompleteness is stated.

**Design: unchanged core** (`plans/coordinates/plan.md` §2)
- Occurrences are the lineage's live coordinate chains, copied in `Walker::end_run`. Every trace run ends there.
- The interval rule:
  - right arm: [c1+k−L, c1+k);
  - left arm: [c1, c1+L).
  - It was verified on 3,631 real runs.
- Per-run coordinates go in a side table, `ArmResult::run_coordinates`. `LabelRun`, `Entry` and `Item` stay unchanged, because the memory model charges `sizeof(LabelRun)`.
- Memory is charged once, when a run is created, at min(chains, cap) occurrences.
- There is no work charge: coordinates were charged when their row was fetched.
- The delivery reserve prices the coordinate share of the account at a fixed 12 account bytes per text byte.
- Python:
  - the subclasses `CoordClaim`, `CoordWalk` and `CoordLabelWalk`;
  - parsed coordinates kept in `g.cache`;
  - eager strict validation (C9).
- Capability detection uses `GET /traverse/capabilities?graph=` (available since level 3).

**Design: revisions after the critique**
1. **Switch-entered runs can understate occurrences.**
   - A run that switches into a label whose own lineage is still live carries only that label's chains that continue from before the switch. Chains that start at the switch node are missing.
   - This contradicts the §18.2 sentence "Runs entered by a switch carry the new label's own coordinates from their start".
   - Fix: `commit_entries` marks such a run in the side table.
   - The response states it (decision C-N2, recommended option a):
     - the run gets `lower_bound: true`;
     - the block gets `complete: false` and a count `runs_lower_bound`.
   - The walk itself is not changed, which would break opt-out trace identity.
   - It only arises with a switch cost (the default is `forbid`).
2. **`drop_coordinates` must work through continuations.**
   - The library removes `max_coordinate_occurrences` whenever the merged `output.coordinates` is false or absent.
   - This applies in `_next_request` and in `build_request`.
   - The server stays strict and refuses a cap without `coordinates: true` with 400.
3. **The reserve's ratio sample must not leak into the server pool.**
   - The exact coordinate text length is computed from digit counts and fixed punctuation.
   - The sample becomes (account − coordinate account)/(text − coordinate text).
   - A test checks that an opt-in attempt leaves `measured_account_per_text_byte` unchanged for later opt-out attempts.
4. **The D4 gate is redefined to measure the right regime.**
   - Each trace cell runs at budgets of 50% and 75% of its own opt-out peak.
   - A column-kind fixture with at least 16 chains per run is added.
   - An offline refseq33m projection is added: R × (3,280 + min(chains, 16) × 872) against the walk account.
   - The run-cost multiplier is restated as 4–20×, depending on chains per run.
5. **M2 gets a dedicated wide coordinate fixture.**
   - The same few-kb sequence is put in at least 5,000 records under one column, with sparse anchors.
   - A header-kind variant uses `CoordToHeader`.
   - R10's many-label fixture is used only for the case of many entries with one chain each.
6. **`compare()` must not clip coordinates.**
   - `_run_claim(…, with_coords=False)` is the default.
   - Coordinates are attached only in the public paths that price them: `_claims_plain`, `_claims` and `_walk_objects`.
   - A test checks that `compare()` and `compare_cost()` are identical with and without the block.
7. **The true output bound is stated.**
   - `--traverse-max-memory-mb` defaults to 0, so it does not bound the block by default.
   - The capabilities rule and the docs state the real bound: occurrences ≤ the coordinates read, each charged 1 work unit at fetch.
   - Ops note O2: set `--traverse-max-memory-mb` on staging and production.
8. **The capabilities diff is stated correctly.**
   - Both GET routes change in:
     - `feature_level`;
     - the `attempts.delivery_reserve.rule` text;
     - a new numeric `attempts.delivery_reserve.coordinate_account_per_text_byte: 12`, so a ledger can reproduce the reserve.
   - `/traverse/capabilities` also gains the `coordinates` block.
   - Per-request capabilities change only in the digit.
9. **The frozen-contract text is aligned with C12.**
   - DESIGN §14.1 and SPEC §7.5 currently freeze limitation kinds.
   - Both are amended: limitation kinds are an open `[a-z_]` token set that readers keep.
   - The precedent is `memory_bound_soft` (28c51c31), not §22.
10. **A coordinate fixture is added for T6.**
    - `WalkerTest.LabelEndCoordinates` (T6) has no k-mer coordinates.
    - Increment 1 adds a BASIC coordinate variant of the T6 fixture, in a column version and a header version.

**Contract additions**
- **Request** (`strategy.output`):
  - `coordinates`: bool, default false, echoed only when true.
  - `max_coordinate_occurrences`: integer ≥ 1 or "unlimited", default 16.
    - Refused with 400 without `coordinates: true`.
    - Accepted and inert under `support: kmer`.
- **Response**, per seed in every detail:
  - either the `coordinates` block, with kind, k, max_occurrences and complete; `seed[]`; and per-arm runs carrying run, label, from_bp, to_bp, occurrences, `occurrences_total` (only when cut), `chains_ended` (only when > 0) and **`lower_bound` (only when true; new)**;
  - or `coordinates: null` with `coordinates_reason`.
- **Limitation** `coordinates`, in no outcome class (C2):
  - knob `output.max_coordinate_occurrences`;
  - `lists_cut`.
- **Resource-stop action:** `drop_coordinates`.
- **MGT free tokens:** K kind `coordinates` with the extra field `lists_cut`, and the action `drop_coordinates` (C12). Nothing else in the body changes.
- **Capabilities:** as in revision 8.
- **Feature level:** 6.
- **Python:**
  - `Coordinates`, `RunCoordinates`, `SeedOccurrences` and the `Coord*` subclasses;
  - `Graphlet.coordinates` and `coordinates_reason`;
  - `to_fasta(coordinates=None)`;
  - `build_request(coordinates=, max_coordinate_occurrences=)`, `supports_coordinates()` and `traverse(coordinates='auto')`;
  - MCP claims rows gain `coordinates`, `coordinates_total`, `coordinates_truncated`, `coordinates_kind` and `coordinates_lower_bound`;
  - `budget.W_COORD = 2`.

**Increments** (each one can be committed)

| # | Increment | Done when | h |
|---|---|---|---|
| C0 | Rebase onto level 5; re-read the touch points (R10 moved the seed phase; R9(b) changed the reserve); add a `--coordinates` mode to `byteid3.py`; build the level-5 base binary | Rebased; base binary built | 0.5–1 |
| C1 | Walker: side tables, seed occurrences, `chains_ended`, lower-bound marking, run-creation charges, two-argument progress, `validate_strategy`; T6 coordinate fixture (column and header); `WalkerCoordinates.*` | Release and Debug green; budgeted opt-out harness on mini_refseq shows 0 differences | 7–9 |
| C2 | Request and JSON layer: parse and echo, `coordinates_json`, the limitation, null reasons in the four failed-seed builders, `delivery_costs(…, coordinates)`, `drop_coordinates`; `GraphletCoordinates.*`; integration tests | Opt-out harness (a) and opt-in harness (b) pass on mini_refseq, UHGG, SRA and /resolve, unbudgeted and budgeted. **The wire contract is frozen here** | 5–6.5 |
| C3 | Reserve split with exact coordinate text; capabilities block and attempts field; level 6 | The capabilities diff equals revision 8 exactly | 3–4 |
| C4 | Real-index tests: `MiniRefSeq.CoordinatesAgainstTheSourceRecords`, now with switch cells (cost 0.5 and 1, loss_budget 2, limit 0 and 2; `next/coords-plan/swreq_*`), and the path-cache equality test; M1 at the cells' stop points; M2 on the dedicated fixture; gcc:12 | Committed with the measurement table and the D4 and D5 data | 5–6.5 |
| C5 | Python library (in parallel with C3–C4 once C2 is frozen): `coords.py`, subclasses, parser, ops with `with_coords` gating, export, client auto-on, cap stripping, MCP, stage-L charges | Golden gate unchanged; full Python suite green | 8–9 |
| C6 | Real suite (`trace_coords`, repeat seeds, column cell, positional oracle, conformance, offline snapshots), CLI fixtures, `bench_traverse.py --coordinates` | Committed | 3–4 |
| C7 | Docs: SPEC §5, §7.0, §7.1, §7.5 (freeze sentence), §10.3, §11.3, §13; DESIGN §14.1 amendment, §25 and v5.8; graphlets.rst | Committed | 2–2.5 |
| C8 | Adversarial review with three lenses (guarantees and byte identity; budgets, reserve and determinism; library, stage L and MCP), fixes, final verify | Ready for your review request | 4–5 |

**Tests and verification**
- **Unit:**
  - `WalkerCoordinates.*`: exact intervals; partition at a split; switch-entered runs; lower-bound marking; zero-length runs; cap cuts; column kind; `chains_ended`; work independence; memory stop not deeper; attempt stops.
  - `GraphletCoordinates.*`: the 400 rule; null reasons; echo; same block in every detail; K record only when cut; stripped response equals opt-out; delivery-cost bound.
  - `GraphletAttempt.ReserveCountsCoordinateText`, plus a no-pool-pollution test.
- **Positional oracle (mini_refseq):**
  - The bases at every [s, e) must equal the spelled run.
  - The occurrence count must equal the string count for unmarked runs, and be ≤ it for lower-bound runs.
- **Byte identity** (byteid3 against the level-5 binary):
  - 588 mini_refseq, 792 UHGG and 804 SRA /traverse requests, plus 688 /resolve.
  - Unbudgeted, and budgeted at 4 MiB, 300k units and 16 MiB.
  - Only the digit and timing are masked. Zero differences.
  - With coordinates on: the stripped response equals opt-out, and the MGT body differs by at most one K record.
- **Determinism:** batch_kmers 1 and 64; path cache off, tiny and default; chunk targets 1 and 50 ms. Stops must be equal at the budget stop points. These cases join R4's and R6's matrices.
- **Builds:** AppleClang -Werror Release and Debug, and the Docker gcc:12 compile of `src/`.

**Measurements**

| | Where | What decides |
|---|---|---|
| M1 | mini_refseq trace cells at their own stop points; switch cells; column cell | D4 gate: depth at the stop with and without coordinates; frequency of lower-bound runs (C-N2) |
| M2 | Dedicated wide coordinate fixture, cap 16 vs "unlimited" | D5: the reserve estimate must be ≥ the real text; text, account, build time, walk_until |
| M3 | Offline: 12 recorded staging trace responses (done) and the refseq33m projection | About 441 B per seed today; the multiplier on taxid columns |
| M4 | Staging after deploy (you): `bench_traverse.py --coordinates` | Elapsed within noise; ec_23S_dV shows 7 occurrences per run (7 rRNA operons); cuts only with taxid columns |

**Risks**
- The fix batch is editing `walker.cpp`, `traverse.cpp`, `traverse_attempts.*` and the timing block. Start only from its commit.
- The layout trap: any future per-run field in `LabelRun`, `Entry` or `Item` changes every budgeted account. Guards: a comment on `LabelRun` and the budgeted harness.
- Budgeted trace requests that ask for coordinates stop shallower (4–20× per run). This is stated, `drop_coordinates` is the lever, and the D4 gate applies.
- On refseq33m, taxid columns cut at cap 16 routinely.
- "unlimited" without a memory maximum gives large responses (ops note O2).
- Assumptions the index cannot reveal: `--fwd-and-reverse` builds and duplicate headers. They are documented, and the oracle catches them on test indexes.

**Critique issues and how each is resolved**

| # | Severity | Issue | Resolution |
|---|---|---|---|
| 1 | high | Switch-entered runs understate occurrences while `complete` stays true | Revision 1; decision C-N2; switch oracle cells |
| 2 | medium | `drop_coordinates` gives 400 through `next_request`, `deepen` and `traverse_continue` | Revision 2; decision C-N6; three tests |
| 3 | medium | The reserve sample leaks into the server pool | Revision 3; test |
| 4 | medium | The D4 gate is vacuous at 4 MiB | Revision 4 |
| 5 | medium | The M2 fixture has one chain per label | Revision 5 (+1–2 h in C4) |
| 6 | medium | `compare()` would clip coordinates without pricing them | Revision 6; test |
| 7 | low | Feature level branch "stays at 4" | P1: the fix batch is 5, coordinates 6 |
| 8 | low | D6 premise false | Revision 7; decision C-N5; ops note O2 |
| 9 | low | Wrong capabilities done-check | Revision 8 |
| 10 | low | D3 case against option (b) overstated | Corrected comparison in decision C-N3; recommendation unchanged, with the reason |
| 11 | low | §14.1 and SPEC §7.5 freeze text | Revision 9; decision X-C12 |
| 12 | low | T6b has no coordinate fixture | Revision 10 (C1 raised to 7–9 h) |

---

### 2.2 Stage 3c-core: selected-label decoding (level 7)

**Goal**
- A constrain-mode read on a row-diff annotation decodes only the selected columns.
- For header labels it decodes only the coordinate windows of their sequences.
- Hits are identical and unbudgeted bodies byte-identical.
- The price is deterministic and conservative.

**Design: unchanged core** (`plans/stage3c/plan.md` §3)
- **Exact decoder** `IRowDiff::decode_selected<RowT>`:
  - windows shift with depth: [lo+i+1, hi+i+1) at depth i;
  - a clear marker is detected from the stored tuple's full size;
  - presence is never XORed from bits;
  - a per-row slack σ keeps cached projections valid stops.
- The same decoder serves binary rd_direct (R5).
- Hits are built by the unchanged `hits_from_*` functions on the restricted rows, so they are identical.
- A projection cache sits next to the level-4 full-row cache.
- A seed in selected mode never touches the full-row cache.
- The crossover is deterministic: 64 distinct columns, set by `--traverse-selected-max-columns` (R4). It never depends on timing or cache state.
- A new knob `strategy.annotation.decode`: auto | full | selected, echoed only when given. "selected" where unsupported is a 400.
- `access_path` is unchanged (R2). The mode goes into timing and capabilities.
- The price goes into the existing `RowCost` fields, so admission code stays unchanged.

**Design: revisions after the critique**
1. **The price (decision X-R3).**
   - 8 per dependency row.
   - Plus, per path row:
     - 1 per rank of the canonical pruned descent, **counted only below a nonzero root rank** (an empty stored row costs 0, as in the full price);
     - 1 per stored selected tuple;
     - 1 per stored coordinate in the canonical window, cancelled ones included.
   - Plus a bounded deep term.
   - A measured price-ratio table (selected/full) is published for mini, UHGG, SRA and the wide fixture.
   - The claim that "budgeted requests get a lower price" is dropped. It is replaced by an acceptance criterion: no index class loses more than 10% median budgeted reach against level 6.
2. **Shift work is bounded.**
   - The decoder shifts lazily. Values stay in the anchor's frame with a per-row offset, are trimmed by binary search, and are merged only at non-empty diffs.
   - Per-step work is then bounded by the windowed stored coordinates that the price already counts.
   - Test: a long path of empty diffs under a large anchor tuple stays within the stated bound.
   - Fallback if lazy shifting proves impractical: a second path aggregate Σ canon(x) × depth(x), added to the price.
3. **The slack D is not a path bound.**
   - The builder bounds row-diff paths only to about twice `max_path_length`.
   - Proposal: D = 256, stated in capabilities as a slack, not a path bound.
   - The deep term becomes Σ(whole stored selected tuple sizes) for paths with L−1 > D. It never exceeds the full read's coordinates for those columns, and it is bounded.
   - Increment 3c.0 measures the real maximum path length on mini, UHGG and SRA. A timing counter (max dependency_rows) reports it on staging.
4. **Every response stays unchanged apart from the digit.**
   - `selected_decode` goes only into `probe_json` (GET /traverse/capabilities).
   - `capabilities_to_json` is written into every response and stays untouched.
   - The identity criterion masks only the digit and timing.
5. **The opt-out lever is defined precisely.**
   - `columns_over_crossover` = (decode == auto) ∧ supported ∧ max_columns > 0 ∧ distinct_columns > max_columns.
   - Byte-identity tests show that budgeted stops under `decode: "full"` and under `max_columns 0` match level 6 exactly.
6. **One cache bound for all members.**
   - Each member's room is bound minus its siblings' bytes, or the accounting moves into `RowDiffPathCache`.
   - Test: `bytes()` ≤ bound with mixed full and selected reads, with and without a memory budget.
7. **The validation query and the walk query are reported separately.**
   - Timing reports each one's mode and column count.
   - D6's reuse benefit is stated as conditional: it holds only when the two selections are equal.
   - A stop names the lever of the query that stopped.
   - Not adopted: serving S' ⊆ S from a projection made under S. The per-selection price aggregates (pruned-tree probes with shared ancestors) do not decompose per column, so prices would depend on the cache.
8. **/resolve stays on full rows in 3c.** `resolve.cpp` passes `Decode::FULL` explicitly. 3b-1 turns selected decoding on for /resolve explicit labels as a planned work item (decision 3c-N7).
9. **Model sizes are pinned.** `static_assert` on `sizeof(SelVisit)` and on the projection entry sizes, or named constants. This is the level-4 `Visit` precedent.
10. **No half-built interruption.**
    - 3c-core gives the new decoder a `DecodeControl*` parameter that is always null.
    - `DecodeStatus::INTERRUPTED` and all real polls arrive together in 3c-ii, with every call site handled.
11. **Default-on scope (decision 3c-N1).**
    - On for tuple (coordinate) row-diff indexes.
    - For binary rd_direct, on only if M3's price ratio and A/B/A pass. Otherwise opt-in through `decode: "selected"` until measured.
12. **Staging expectations are split by run type.**
    - Warm runs: at least 10× lower ms per row.
    - First-touch runs: measured with page faults per row. "Reach the radius" is not expected for `.first` runs.
    - The claim that 3c removes the R10 first-touch regression is conditional on R10's attribution.

**Contract additions**
- **Request:** `strategy.annotation.decode` (auto | full | selected), echoed only when given; 400 for "selected" where unsupported or in annotate mode.
- **Flag:** `--traverse-selected-max-columns N` (default 64; 0 means never chosen by auto) on the server and the CLI.
- **Timing only** (R2):
  - per query (validation, walk): decode mode and `selected_columns`;
  - counters: `selected_rows`, `selected_visits`, `selected_probes`, `selected_coords`, `selected_cache_hits`, max dependency_rows.
- **GET /traverse/capabilities:**
  - `selected_decode {supported, max_columns, window_slack: 256, applies_to, rule, price, echo}`;
  - amended `work_bound` and `decode_cache.rule` texts.
- **GET /capabilities:** `features += "selected_decode"`.
- **Budgeted requests**, intended and stated: `usage.work_units`, stop amounts and stop positions change for constrain-mode seeds in selected mode. Decode and work stops offer `name_fewer_labels` or `lower_max_seed_labels`.
- **MGT:** nothing new. The actions reuse existing tokens; the messages are free text.
- **Feature level:** 7.
- **Python:**
  - `decode` passes through, and `next_request` carries it.
  - `capabilities()` passes the block through.
  - The bench reads the new timing fields and tolerates their absence.

**Increments**

| # | Increment | Done when | h |
|---|---|---|---|
| 3c.0 | Rebase onto level 6; read R6's header-range outcome and R10's final cache API; add permitted-range filtering for full reads if R6 did not; measure max path lengths | Byte identity; `coords_mapped` lower (timing only); path-length table | 2–3 |
| 3c.1 | Library primitives, not used yet: `BRWT::prune` and `probe`, `ColumnMajor::probe`, `TupleCSCMatrix::tuple_extent` and `column_values`, `Selection` | `BRWTProbeColumns`, `ColumnMajorProbeColumns` and `TupleRowDiffWindow` pass; no behaviour change | 5–6 |
| 3c.2 | Decoder, not wired: trace, dmax, valid stops, σ, lazy shift, reconstruction, price (X-R3), bounded deep term, demand, refusal, projection cache with the shared bound, `admit_as_uncached` mirror, static_asserts | `RowDiffSelectedDecode.*` green, including the denial sweep, jemalloc, the long-empty-path bound and batch/cache independence; mini_refseq exhaustive check passes at library level (expected 966,495 checks, 0 mismatches) | 14–17 |
| 3c.3 | Wiring: LabelOracle and LabelQuery modes, walker, timing per query, stop actions with the lever rule, knob and flag, probe-only capabilities, `Decode::FULL` at resolve.cpp, level 7 | `LabelOracleSelected`, `WalkerSelectedDecode`, `MiniRefSeqSelected.WalksIdentical` and `TestTraverseSelectedDecode` pass; identity run: ≥ 1,000 unbudgeted requests identical; budgeted differences only in the listed classes; opt-out tests | 10–12 |
| 3c.4 | Docs (SPEC §6.8, §7.0, §8.2, §8.4, §10.3, §11.3, §13; DESIGN §26, v5.9), Python parity tests, bench `selected` suite | Python suite green; old baselines still read | 4–5 |
| 3c.6 | Measurements M1–M3 and the price table; adversarial review with four lenses (exactness; budgets and determinism; cache and R10; statements); fixes; gcc:12 | Ready for your review request | 10–12 |

The synthetic index builds (WT-16S, WT-Ecoli, WT-long, each under 5 GB) take 3–6 h of wall time. They run in the background during 3c.1–3c.2, at low priority.

**Tests and verification**
- Unit and real-index tests as in the stage plan (`RowDiffSelectedDecode.*`, `TupleRowDiffWindow.*`, `LabelOracleSelected.*`, `MiniRefSeqSelected.*`, `WalkerSelectedDecode.*`), plus:
  - `ProjectionCacheSharesOneBound`;
  - `LongEmptyPathWorkWithinPrice`;
  - `DeepTermBoundedByWholeTuples`;
  - `OptOutLeverMatchesLevel6` (decode full, max_columns 0, budgeted).
- **The coordinates gate is part of 3c's gate:**
  - the positional oracle with `output.coordinates: true`;
  - the opt-in byte-identity set.
- **Byte identity:** ≥ 1,000 unbudgeted requests against the level-6 binary (588 mini_refseq, 396–792 UHGG, SRA summary), decompressed, with timing off and on (timing masked).
- **Budgeted:**
  - `decode: "full"` is identical to level 6;
  - with auto, a masking script requires every difference to fall into a listed class.
- **Determinism:** budgeted selected requests are identical across batch_kmers {1, 2, 7, 64, 1000}, chunk targets, cache 0/1 MiB/128 MiB, cold and warm, seed permutations, and one seed per request.

**Measurements**

| | Where | What |
|---|---|---|
| M1 | Synthetic wide tuple indexes | ms per row, selected vs full (alone, batch 64, in a walk; cache on and off); crossover sweep 1–1024 columns and 1/16/256 headers per column; demand vs jemalloc peak; price ratio |
| M2 | mini_refseq | Exhaustive correctness; ≤ 10% slower on narrow rows; price ratio |
| M3 | Local UHGG and SRA | rd_direct byte identity and A/B/A elapsed; price ratio (decides 3c-N1 for binary indexes) |
| M5 | Staging after deploy (you) | Warm named and trace walks on ec_rplJ, ec_23S_dV, pa_gyrB and ndm1 at ≥ 10× lower ms per row; first touch with page faults; `selected` suite; crossover probes |

**Risks**
- The window and slack logic is subtle. Mitigations: the invariant proof, targeted tests, the exhaustive mini check and an exactness review lens.
- Exactness rests on four builder invariants: sorted tuples; clear markers only as size-0 tuples on non-anchor rows; full anchors; SHIFT = 1. A debug verifier checks them.
- Prices depend on server flags (max_columns, D). A ledger reads them from capabilities.
- Paths longer than D re-decode more often and get a positive deep term. Stated.
- R10 may change the cache API; rebase on its final form. Comparisons are always against level 6.
- Cold graph I/O remains on first touch. Stated.

**Critique issues and how each is resolved**

| # | Severity | Issue | Resolution |
|---|---|---|---|
| 1 | high | INTERRUPTED would become a false DECODE stop | Moved to 3c-ii, where every call site is handled before the REFUSED branch, with a poll-ordinal sweep test (§2.3) |
| 2 | high | The in-read bound is overclaimed: hit building and the 256-row fallback are not covered | 3c-ii polls hit building, `map_coords` and recorder naming; the ops bound is advertised only for checkpointed paths (§2.3) |
| 3 | high | The price can exceed level 5 on sparse-diff indexes | Revision 1; decision X-R3; ratio table; integration assertion limited to the wide fixture |
| 4 | medium | `selected_decode` in every response breaks byte identity | Revision 4 |
| 5 | medium | D = 128 is not a path bound; the deep term is unbounded | Revision 3; decision 3c-N9 |
| 6 | medium | Shift work is not counted | Revision 2 |
| 7 | medium | D8 would reroute almost every full read | Separate increment in 3c-ii, behind its own flag, default off; `decode: "full"` stays on the level-6 path (decision 3c-N8) |
| 8 | low-med | Cache bound is per member | Revision 6 |
| 9 | low-med | Validation vs walk selections | Revision 7 (subset reuse rejected, reason given) |
| 10 | low-med | /resolve ownership unclear | Revision 8; 3b-1 work item |
| 11 | low-med | Unsupported M5 expectations | Revision 12 |
| 12 | low-med | Effort optimistic | Re-estimated; external review pass budgeted; index builds as wall time; 3c.5 split out |
| 13 | low | `sizeof(SelVisit)` not pinned | Revision 9 |
| 14 | low | Opt-out lever undefined | Revision 5 |

---

### 2.3 Stage 3c-ii: checkpoints inside reads (level 8)

**Goal**
- A /traverse read on a row-diff format can be stopped within a stated number of index operations, through a deadline, a cancel or an attempt stop.
- Nothing is advertised that does not hold.

**Design**
- `DecodeControl {stop, interval = 4096 ops, poll(n)}`.
  - An op is one trace step, rank or select.
  - n coordinates count as ⌈n/64⌉ ops.
- `DecodeStatus::INTERRUPTED` behaves like a refusal: no side effects, budget restored.
- **Every call site handles INTERRUPTED before its REFUSED-and-retry branch.**
  - The sites at f667d775 are `label_oracle.cpp` 972, 1168, 1241, 1664, 1879, 1895, 1916 and 2035; `walker.cpp` 1396; `tuple_row_diff.hpp` 143. They are re-located after the rebase.
  - Each maps to `FetchRefusal::INTERRUPTED` → `ReadPacing::interrupted` / `paced_trip`, charging the units of the keys already built.
  - There is no `switch` over `DecodeStatus` anywhere, so `-Wswitch` finds nothing. The site list is the checklist.
- **Polls:**
  - in the selected decoder;
  - in `decode_budgeted` (trace, `BRWT::descend_row` per node, `row_tuples` per tuple, `xor_into` between columns);
  - in hit building (`hits_from_tuples`, `hits_budgeted`, `map_coords`, `LabelRecorder` naming), every 4096 coordinates. The hit-cache insert comes after, so there is no side effect.
- **Routing of full reads (D8)**, behind a new server flag that defaults to off:
  - Paced full reads on row-diff go through `decode_budgeted` with an unlimited budget.
  - Almost every read counts as paced, because `time_budget_ms` defaults to 30 s. The flag is therefore off until A/B/A passes on every bench suite.
  - `decode: "full"` stays on the level-7 path.
- **`deadline_check`, on both capabilities routes:**
  - `checkpoint_ops: 4096`;
  - `checkpointed`: the {format, path} pairs actually polled. Selected reads always; full reads only with the flag on; never the 256-row fallback;
  - `max_uninterruptible_ops` for those pairs;
  - `max_uninterruptible_ms` stays null (decision 3c-N5);
  - the rule text is amended on top of R5's.

**Contract additions:** the `deadline_check` fields above; the new flag; feature level 8. No response body changes.

**Increments**

| # | Increment | Done when | h |
|---|---|---|---|
| 3c.5a | INTERRUPTED at every site; polls in the decoders and in hit building; poll-ordinal sweep test across fetch, warm, recorder, derivation window and validation | No sweep point yields a DECODE or DEMAND refusal; the walk stops where a stop between chunks would; ops between polls ≤ 4096 on a wide fixture, hit building included; byte identity far from the deadline | 8–10 |
| 3c.5b | D8 routing behind the flag; A/B/A on walks, annotate, budgets.derived, lookahead and batch, on UHGG, SRA, mini and the wide fixture | Flag-on is no slower than level 7 beyond noise, or it stays off and this is stated | 4–6 |
| 3c.5c | `deadline_check` fields and text, docs (DESIGN §27, v5.10), review, verify, gcc:12 | Ready for your review request | 3–4 |

**Measurements:**
- local: ops between polls; A/B/A overhead;
- staging (you): R8's per-seed longest-piece record, before and after.

**Risks**
- Poll overhead on narrow rows. It is measured, and the interval can be tuned.
- The flag may stay off. Then full reads remain unbounded per piece, which is stated.

---

### 2.4 Stage 3b-1: /resolve budgets and level rows (level 9)

**Goal**
- /resolve accepts time, memory and work budgets.
- A stopped /resolve is exactly the resolve of a query prefix (B7).
- Under a memory budget, a /traverse level's rows count in the account (B5, B6, B15).

**Design: unchanged core** (`plans/stage3b/plan.md` §2)
- Strict optional bounds: `time_budget_ms`, `max_memory_mb`, `max_work_units`.
- One `ProfileBuilder` serves the budgeted and unbudgeted paths. The unbudgeted output is byte-identical (decision 3b-N2).
- Fields appear only under a budget: `outcome`, `resolved`, limitations (`memory_bound_soft`, `query_prefix`), `resource_stop`.
- `--resolve-max-time-ms`, default 0 (B11). No server memory default (B12).
- Walker change (c): `level_held_` joins `accounted()`. Each head's rows are released after the head. A `LEVEL_ROWS` statement without `more_selective_seed`. The lookahead key list is reserved.

**Design: revisions after the critique**
1. **x does not depend on the window.**
   - k-mer j is admitted against its own row plus the profile and the reservations only.
   - The window's later rows are droppable, like the walker's lookahead (admitted as uncached). When j does not fit, they are released and fetched again.
   - Without a memory budget, windows are sized by pacing and by work, as traverse does.
   - Test: `ResolveBudget.PrefixIndependentOfWindow` (x equal for W ∈ {1, 7, 64} under memory, work and time budgets, with an injected clock).
2. **Work is charged for every decoded row.**
   - Under a work budget each fetch is sized like `fetch_chunk`: max(1, units_left / (8 + widest row so far)), growing 1, 2, 4….
   - Every row a fetch decodes is charged in `resource_stop.used` and `resolved.work_units`. This does not touch the prefix property, because those fields are excluded from it.
   - The overrun bound is stated in `resolve_budget_rule`.
   - The work formula differs from the Oct-3 SPEC text (decision X-B8b).
3. **B10 follows the approved wording.**
   - Explicit intervals wholly inside the prefix are selected from it.
   - An interval that ends past x fails the selection. The shape is decision X-B10.
   - `selection` is excluded from the prefix comparison only in that failure case.
   - Tests: `PrefixProperty` is fixed to match; new `ExplicitSeedsInsidePrefixAreSelected`.
4. **The 3c hand-offs are explicit work items.**
   - Explicit-label /resolve reads through `decode_selected`, with the X-R3 price terms, the mode echoed in timing, and byte identity.
   - /resolve windows poll `DecodeControl`. INTERRUPTED becomes a time stop at an exact prefix.
5. **Select_seeds is bounded.**
   - It is restructured byte-identically: it sorts label ids instead of copying names.
   - Reservations are derived from the loops:
     - runs: R × (1 + max_seeds) × 2 × sizeof(Run);
     - seed sequence: ×3;
     - covering labels per seed.
   - If the work charge is not a loop count, it is stated as nominal.
   - `ModelBoundsTheHeap` covers labels above the cap and `max_seeds` ≥ 10 with overlapping runs.
6. **Discovery's delivery reservation is bounded.**
   - min(labels_met, max_labels) × (fixed JSON + longest escaped name seen) + Σ runs × run JSON bound. It is monotonic and needs no heap.
   - The B9 "smaller than level 4" claim is made only after measurement on the ec_23S_dV, pa_gyrB and rplL loci (wide fixture) and on mini_refseq.
7. **Limits plumbing.**
   - `ResolveLimits` gains `chunk_target_ms` and `path_cache_bytes`, beside `max_time_ms` and `max_labels`.
   - The server and the CLI pass them.
   - Test: a budgeted /resolve uses the cache and is paced.
8. **Time-only budgets use the default decode with pacing**, as traverse does. Memory and work budgets use the budget-aware path.
9. **`--resolve-max-time-ms` when set:**
   - It lowers larger time budgets and sets omitted ones (decision 3b-N6).
   - The docs and the ops note state that this changes the bytes, and the engine, of every unbudgeted /resolve.
   - Turning it on in production is gated on a staging comparison against level 8, with a revert path.
10. **What stays uninterruptible is stated.**
    - Listed in `deadline_check.rule` and `resolve_budget_rule`: label resolution and the header-lookup build; per-window k-mer mapping; `finish()`; `select_seeds`; response building.
    - `time_budget_ms` bounds the reading of the profile only.
    - Timing records the longest piece by kind.
11. **The discover clamp is dropped.**
    - `--resolve-max-labels` caps only the explicit label list (400 when over), as you asked.
    - `clamped` appears only for `bounds.time_budget_ms` (decision X-B8a).
12. **D2 gates are local.**
    - Local replays of /resolve on mini_refseq, UHGG, SRA and the wide fixture.
    - Elapsed ≤ 1.02× level 8 and peak memory ≤ level 8, for both 3b-N2 options.
    - Staging is a confirmation after deploy, with a fallback switch.
13. **Walker change (c) leaves the walk-until unchanged.**
    - `progress()` receives accounted() − `level_held_`, or the level is fully released before progress.
    - `usage.per_seed[].memory.peak_admitted_bytes` now includes the level's rows. This is stated in the SPEC §10.3 level-9 line and in the diff classifier (decision 3b-N16).
    - An attempts test checks that the walk-until is unchanged for a memory-budgeted row-diff walk with no level stop.

**Contract additions**
- **/resolve request:** `bounds.time_budget_ms`, `bounds.max_memory_mb`, `bounds.max_work_units`.
- **/resolve response under a budget:**
  - `outcome {profile, selection?}`;
  - `resolved {kmers, query_kmers, bp, query_bp, work_units}`;
  - limitations `memory_bound_soft` and `query_prefix`;
  - `resource_stop`;
  - `clamped` (time budget only);
  - timing additions.
- **JSON-only tokens:**
  - limitation `query_prefix`;
  - scope `request`;
  - phase `profile`;
  - actions `shorten_query`, `resolve_remainder`, `name_labels`;
  - outcome values `complete`, `partial`, `failed`.
- **/traverse:**
  - no new field or token;
  - `LEVEL_ROWS` is free text;
  - the memory-budgeted `peak_admitted_bytes` change is stated.
- **Flags:** `--resolve-max-time-ms` and `--resolve-max-labels`, both default 0.
- **GET /traverse/capabilities only:**
  - `resolve_budgets`;
  - `resolve_max_time_ms`;
  - `resolve_max_labels`;
  - `resolve_budget_rule`;
  - the `deadline_check.rule` text.
- **MGT:** nothing.
- **Feature level:** 9.
- **Python:**
  - `client.resolve(bounds=)`;
  - MCP `traverse_resolve` budget arguments;
  - outcome, resolved and stop placed first in results, never cut by `max_bytes`.

**Increments**

| # | Increment | Done when | h |
|---|---|---|---|
| 3b.0 | Reference binary (level 8); bench baselines; /resolve peak-memory harness; read R6's final resolve.cpp | Identity replay of the reference against itself (noise calibration) | 1–2 |
| 3b.1 | /resolve budgets with revisions 1–12; Python and MCP parity; SPEC §4 and §10.3 | Unbudgeted identity on ≥ 1,000 requests (/resolve included); `ResolveProfileBuilder.*` and `ResolveBudget.*` (prefix, window independence, determinism, denial sweep, heap model, delivery bound, client gone); D2 local gates; gcc:12 | 18–24 |
| 3b.2 | Walker change (c) with revision 13; fixture regeneration with a classified diff; ASan | `GraphletStage3bLevel.*` under ASan; UHGG measurement E (constrain unchanged, annotate earlier); identity for results with no observed level excess; attempts walk-until test | 7–9 |

**Measurements**
- **Local:**
  - mini_refseq queries of 10 kbp, 100 kbp and 1 Mbp: explicit, discover column, discover header and trace; x per budget; overhead per budget instead of a fixed 1.10× target; peak memory against level 8;
  - UHGG and SRA: discover column;
  - wide fixture: discover header and trace;
  - a refseq33m projection table.
- **Staging (you):**
  - `resolve.ec_23S_dV.{explicit, discover_header}.{mem64, work2M, time500}`;
  - `resolve.r16S_univ.discover_header.time5s`;
  - 10 kbp queries;
  - a client-side prefix resubmission check.

**Risks**
- Conservative reservations stop /resolve early. Mitigations: the bounded delivery term, droppable window rows and the measured x per budget.
- R6 and R9(c) rewrite the same /resolve code. 3b.1 rebases on their final form.
- Platform-dependent model constants (libc++ vs libstdc++). Tested per platform.
- Walker change (c) churns `memory_bound_soft` text in every memory-budgeted row-diff response.

---

### 2.5 Stage 3b-2: remaining annotation formats (level 10)

**Goal:** budget-aware reads for 13 more formats (B1). The four legacy formats stay stated as uncovered.

**Design: unchanged core** (`plans/stage3b/plan.md` §2, part 3)
- A `BudgetedReader` factory (dynamic_cast chain).
- `decode_single_rows<RowT>`.
- ColumnMajor as charged column-wise pairs (B2).
- BRWT single-row descent only under a memory budget (B3).
- Disk readers charged as a measured constant `kDiskReaderBytes` (B14).
- An `IntRowDiff` count decode.
- Budgeted DIRECT replaces the assert at `label_oracle.cpp:968`.
- Work-only budgets keep the default decode on plain formats (B4).
- Format-aware statements.
- `annotation_reads` on GET /traverse/capabilities only (B16).

**Design: revisions after the critique**
1. Every new reader polls `DecodeControl` and honours INTERRUPTED. That covers single-row reads, ColumnMajor pairs, disk fetchers and the IntRowDiff decode.
   - Each format gets a poll-ordinal sweep test.
   - Each covered format joins 3c-ii's `checkpointed` list.
2. **B3 is kept as you decided.** The production measurement is possible: the local UHGG and SRA `row_diff_brwt` indexes contain production BRWTs.
   - Measured: single-row descent against batched slice.
   - Separately for anchor rows and diff rows, for consecutive and random rows, and for 16S/23S loci.
   - The ratio is reported.
3. **B15 and `direct_reads`.**
   - Once budgeted DIRECT is active on plain formats (for example `column_coord`, the CLI default), skipping the lookahead changes `annotation.direct_reads`.
   - The B15 SPEC text states this, as level 4 did.
   - The B15 test asserts `direct_reads` separately for DIRECT formats.
4. **Selected reads on plain formats are not in this program.** 3b's DIRECT path uses `get()`, as planned. 3c's probe primitives are named correctly in the docs (`BRWT::probe`, `ColumnMajor::probe`, `TupleCSCMatrix::tuple_extent`) for a later stage.

**Contract additions**
- `annotation_reads {budget_aware, kind}` on GET /traverse/capabilities.
- The `work_bound` text: dependency rows only on row-diff.
- Feature level 10.
- Memory-budgeted /traverse on the 13 formats is charged. Work-only budgets on plain formats are unchanged.

**Increments**

| # | Increment | Done when | h |
|---|---|---|---|
| 3b.3 | Library formats (annotation code only), with polls, sweeps and the B3 measurement; in a worktree, may start during 3b.2 | Disassembly identity of the default decode functions; runtime ratio 0.98–1.02; jemalloc bounds on macOS and in Docker; gcc:12 | 13–17 |
| 3b.4 | Integration: `reader_`, the `decode_charged_` rule (B4), budgeted DIRECT, statements, `column_coord` fixtures, new format builders and probes | Unbudgeted and work-only identity; budgeted diff classified; `work_stop` fixture identical | 9–12 |
| 3b.5 | Capabilities, docs (DESIGN §29, v5.12), final gates | Ready for your review request | 3–5 |

**Risks**
- Fixture churn on `column_coord`.
- The disk constant depends on libc++/libstdc++ and sdsl internals. If jemalloc cannot run in Docker, the fallback is `mallinfo2`, which is weaker and stated (decision 3b-N13).
- No traversal fixture builds several of these formats today. New builders are needed.

**Critique issues for 3b (both parts) and how each is resolved**

| # | Severity | Issue | Resolution |
|---|---|---|---|
| 1 | high | The window rule starves the profile (x = 0) | 3b-1 revision 1 |
| 2 | high | Rows decoded past x are not charged | 3b-1 revision 2; decision X-B8b |
| 3 | high | The prefix property contradicts B10 | 3b-1 revision 3; decision X-B10 |
| 4 | high | 3c hand-offs unowned | 3b-1 revision 4; 3b-2 revisions 1 and 4 |
| 5 | medium | Walker change (c) moves `peak_admitted_bytes` and the walk-until | 3b-1 revision 13; decision 3b-N16 |
| 6 | medium | D8 waived B3 on a false premise | 3b-2 revision 2 (B3 kept; no decision needed) |
| 7 | medium | D9 over-reserves about 25× | 3b-1 revision 6 (supersedes old D9) |
| 8 | medium | The selection reservation does not bound `select_seeds` | 3b-1 revision 5 |
| 9 | medium | D6 plus the ops note silently move /resolve to the new engine | 3b-1 revisions 8 and 9 |
| 10 | medium | `chunk_target_ms` and `path_cache_bytes` not plumbed | 3b-1 revision 7 |
| 11 | medium | Level numbering conflict | P1; levels relative to previous + 1 |
| 12 | low | `deadline_check` text overclaims | 3b-1 revision 10 |
| 13 | low | D2 gates need staging before merge | 3b-1 revision 12 |
| 14 | low | B15 and `direct_reads` | 3b-2 revision 3 |
| 15 | low | Discover clamp is a scope addition | 3b-1 revision 11 (dropped) |
| 16 | low | Effort and overhead target | Re-estimated; per-budget overhead measurements |

---

## 3. Cross-stage interactions

### 3.1 With the level-5 fix batch (it lands first; every stage rebases on its commit)

| Item | Coordinates | 3c | 3b |
|---|---|---|---|
| R1 tombstones, `expect_server_instance` | Merge point in `client.build_request` parameters | — | — |
| R2 bounded-integer parser | — | `--traverse-selected-max-columns` uses it | `--resolve-*` flags use it |
| R3 identity inventory | Header positions need the `.seqs` sidecar; without it the kind is column | No new file loaded | Disk formats reopen covered files only |
| R4 walk-order chunk grouping | Coordinates must be identical; add opt-in requests to its identity test | Selected reads inherit it; extend its decoder-count test | Budgeted DIRECT and /resolve first-occurrence order build on it |
| R5 delivery checkpoints | Block built with `delivery_tick()` per element | 3c-ii rebases on its `deadline_check` text | /resolve responses use the same writer |
| R6 efficiency audit | Add opt-in requests to its cache/batch/chunk matrix | Its header-range outcome decides 3c.0; matrix extended to selected reads | Its /resolve rewrite is the base of `ProfileBuilder`; matrix extended to budgeted /resolve and (c) |
| R7 SPEC, DESIGN §24 | DESIGN §25 | §26 (core), §27 (checkpoints) | §28 (3b-1), §29 (3b-2) |
| R8 per-seed overrun record | Merge in the timing block of `seed_result_to_json` | Record gains the read mode; 3c-ii shortens pieces | /resolve timing reuses the longest-piece kinds |
| R9(a) injected clock | — | Reused for checkpoint tests | Template for /resolve time-stop tests |
| R9(b) per-detail ratios | The reserve split builds on `ratio_locked` and `set_delivery_detail` | — | — |
| R9(c) explicit /resolve speed | — | — | Re-measure explicit /resolve against level 8 |
| R10 selective path caching, seed phase | Rebase on its seed-phase edits in `validate_seed`/`run()`; its fixture only for many entries | Projections live in its cache; adopt its API and copy-volume counters; compare against level 6, never level 4 | /resolve windows are its successor-rule case; its counters go into /resolve timing |
| P7 continuation validation | Merge in `_next_request` (cap stripping) | — | — |
| P8 memo audit | Coordinate index built once per parse | — | — |
| P9 bench reach and work columns | `--coordinates` columns beside them | Staging comparison by certified reach and consumed work | Resolve cases report `resolved.kmers` and work |
| P10 index thresholds | Lookups by run id; no per-call index | — | — |

### 3.2 Path cache and R10
- **Coordinates:**
  - Coordinate charges belong to the account, not to the cache allotment.
  - Responses must be byte-equal with the cache on and off.
- **3c:**
  - Projection members join `RowDiffPathCache` under one shared bound.
  - Selected seeds never touch the full-row cache. Full reads still need R10's fix.
- **3b:**
  - Budgeted /resolve uses the cache under traverse's allotment rules. That needs `path_cache_bytes` plumbed.
  - `RowDiff<RowFlat/RowSparse/RowDisk>` gain `PathAggregates`. IntRowDiff gets no cache.
- **Invariant everywhere:** `RowCost` and admission never depend on cache contents. R6's matrix tests this for every stage.

### 3.3 Attempts, deadlines and the delivery reserve
- **The lease bound formula is untouched by every stage.**
- **Coordinates:**
  - The reserve gains a conservative coordinate term.
  - Its ratio samples exclude the coordinate text exactly, so opt-out walk-until is unchanged.
- **3c-core:** the walk-until and the reserve are unchanged.
- **3c-ii:**
  - INTERRUPTED feeds the existing `paced_trip`, `cancelled` and `attempt_deadline` stops.
  - Attempt bounds hold more tightly.
- **3b-1:**
  - Walker change (c) must not move the walk-until; `progress()` excludes `level_held_`.
  - /resolve is not an attempt (decision 3b-N11).

### 3.4 Library (`metagraph.traverse`)
- **Merge points:**
  - `client.build_request`: R1's fields, the coordinates fields, the `decode` passthrough.
  - `_next_request`: P7's validation, coordinates cap stripping.
- **Stage L** charges only where coordinates exist (`W_COORD`). `WORK_MODEL` stays 1, so the golden gate over `committed.json.gz` is unchanged.
- **Every stage keeps Python parity** in the same level:
  - coordinates: the full library;
  - 3c: passthrough and bench;
  - 3b: `resolve(bounds=)` and MCP.

### 3.5 Shared files and merge order
- **Hotspots:**
  - `walker.cpp` admission: coordinates charges, then 3b's (c);
  - `traverse.cpp`: coordinates JSON, 3c actions and timing, 3b /resolve;
  - `server.cpp probe_json`: all stages add blocks; 3c-ii and 3b both edit `deadline_check`;
  - `label_oracle.cpp`: 3c dispatch, 3c-ii sites, 3b `reader_` and DIRECT;
  - SPEC §10.3: every stage.
- Strict sequencing removes conflicts. Each stage starts from the previous stage's commit.
- **DESIGN numbering is fixed now:**
  - §24 (v5.7) fix batch;
  - §25 (v5.8) coordinates;
  - §26 (v5.9) 3c-core;
  - §27 (v5.10) 3c-ii;
  - §28 (v5.11) 3b-1;
  - §29 (v5.12) 3b-2.
- The stage plans each claimed §25. This numbering resolves that.

### 3.6 Shared gates (they grow per stage)
- The byte-identity harness (`p6/cpp/byteid3.py`) gains one mode per stage:
  - `--coordinates`;
  - decode modes;
  - /resolve budgets.

  Each stage reruns all earlier modes against the previous level's binary.
- The determinism matrix (R6) gains each stage's budgeted cases.
- The coordinates positional oracle, including the switch cells, is part of every later gate.
- The canonical capabilities diff is checked per stage against an expected list.
- Docker gcc:12 compiles every increment that touches `src/`.

---

## 4. Owner decisions

You can answer "accept all recommendations" and change only the items you disagree with.

The IDs keep the stage plans' numbers. A missing number (3c-N2; 3b-N3, N4, N5, N8, N15) is an item that moved to §4.3 as a change to an existing decision, or that the critique resolved (the old 3b D8 is now "B3 kept").

### 4.1 Urgent (before the fix batch commits)

| ID | Question | Options | Recommendation |
|---|---|---|---|
| **P1** | Feature level of the fix batch | (a) 5; (b) stay at 4, as its task text says | **(a).** Level 4 is deployed on staging. Then 6 coordinates, 7 3c-core, 8 3c-ii, 9 3b-1, 10 3b-2, each assigned as previous + 1 at commit. The main session sets 5 in the commit step |

### 4.2 New questions

**Program**

| ID | Question | Recommendation |
|---|---|---|
| P2 | Order | Coordinates → 3c → 3b (§1); alternatives B and C rejected |
| P3 | Split 3c into core (7) and checkpoints (8)? | Yes |
| P4 | Workflow shape (§5.2) | Chunks of 1–3 increments; the main session commits between chunks; worktrees only for 3b.3 and reference binaries |

**Record coordinates**

| ID | Question | Options | Recommendation |
|---|---|---|---|
| C-N1 | `chains_ended` definition (clarifies C5) | As drafted: counts chains that continue on no followed path of the lineage; chains partitioned among sibling clones are not counted; clones inherit the prefix's count | Adopt |
| C-N2 | Switch-entered runs into a still-live label (see also X-18.2) | (a) mark: per-run `lower_bound: true`, block `complete: false`, `runs_lower_bound`, capabilities sentence; MGT unchanged; (b) exact: carry an output-only chain set in a side table (+4–6 h, more memory charge); (c) a separate limitation kind | **(a)** now. Revisit (b) if M1 shows such runs on real requests |
| C-N3 | Python API shape | (a) `Coord*` subclasses; (b) a base-class slot listed in `ADDITIVE_FIELDS` (golden digests stay unchanged) | **(a).** Corrected comparison: (b) keeps digests but moves stage-L finite-budget stop points by 8 B per row for graphlets without coordinates, unless `_CLAIM_BYTES` is compensated. (a) costs two record types per kind, and `CoordClaim != Claim` under dataclass equality (no code path compares them) |
| C-N4 | Delivery reserve | (a) a fixed bound of 12 for the coordinate share, with ratio samples using the exact coordinate text; (b) let the measured ratio absorb it | **(a)** |
| C-N5 | Server cap on occurrences | (a) none in v1, true bound stated, ops note O2; (b) a `--traverse-max-coordinate-occurrences` flag; (c) clamp "unlimited" when no memory maximum is set | **(a)**; revisit after M4 |
| C-N6 | A cap sent with `coordinates: false` or absent | (a) the server keeps the 400; the library strips the cap; (b) the server accepts it as inert | **(a)** |
| C-N7 | A cap under `support: kmer` with `coordinates: true` | Accept as inert (the null reason explains), or 400 | Accept as inert |
| C-N8 | Coordinates in GFA export | No in v1 | No |
| C-N9 | `compare()` with coordinates | A note "coordinates are not compared (v1)"; no clipping in compare paths | Adopt |

**Stage 3c**

| ID | Question | Options | Recommendation |
|---|---|---|---|
| 3c-N1 | Default on, and the knob | (a) on (`decode: auto`) with opt-outs per request and per server; (b) opt-in; (c) on for tuple row-diff, binary rd_direct opt-in until M3 passes | **(c)**, using the new knob `strategy.annotation.decode` |
| 3c-N3 | Crossover | (a) R4 literal, 64 distinct columns; (b) also skip indexes with ≤ max_columns columns in total; (c) a cost estimate | **(a)**; tune from the sweep and staging data |
| 3c-N4 | rd_direct echo | `access_path` stays "rows" | Yes |
| 3c-N5 | What to advertise after checkpoints | (a) `max_uninterruptible_ms` null, `max_uninterruptible_ops` only for checkpointed paths; (b) a finite ms value per server; (c) a finite ms value with `--traverse-assume-resident` | **(a)**; (c) later if staging measures the tail (refseq33m is not resident) |
| 3c-N6 | Projection cache in the seed phase without a memory budget | On / off | On; the benefit is conditional on equal selections |
| 3c-N7 | /resolve explicit labels | (a) 3c keeps `Decode::FULL`, 3b-1 turns selected on; (b) 3c turns it on for unbudgeted /resolve | **(a)** |
| 3c-N8 | Routing paced full reads through `decode_budgeted` | Separate flag, default off, enabled after A/B/A on every suite | Adopt |
| 3c-N9 | Slack D | 128; 256; measured | **256**, stated as a slack, with the bounded deep term |

**Stage 3b**

| ID | Question | Options | Recommendation |
|---|---|---|---|
| 3b-N1 | Order inside 3b | /resolve → (c) → formats | Adopt |
| 3b-N2 | Profile engine | (a) shared `ProfileBuilder`, unbudgeted loop as R6 leaves it; (b) unbudgeted also streams | **(a)**; (b) if the local gates pass |
| 3b-N6 | `--resolve-max-time-ms` when set | (a) sets omitted and lowers larger budgets, echoed; (b) clamps only given budgets | **(a)**, with the bytes-change statement and gated production use |
| 3b-N7 | Window rule | x independent of W, droppable window rows | Adopt (replaces the old D7) |
| 3b-N9 | Discovery delivery reservation | The bounded formula (revision 6) | Adopt (replaces the old D9) |
| 3b-N10 | Selection reservation | Restructure `select_seeds` byte-identically and derive the reservations | Adopt |
| 3b-N11 | /resolve as a ledger attempt | In 3b or later | Later |
| 3b-N12 | Path cache in budgeted /resolve | On, under traverse's rules | On |
| 3b-N13 | B14 validation on Linux | (a) jemalloc in Docker gcc:12; (b) a stated margin with `mallinfo2` | **(a)**; (b) as fallback, stated as weaker |
| 3b-N14 | Levels | Two bumps (9, 10) or one | Two |
| 3b-N16 | `peak_admitted_bytes` after (c) | Accept the stated change (more honest) with the walk-until kept unchanged | Accept |

### 4.3 Proposed changes to existing decisions

| ID | Existing decision | Proposed change | Why | Recommendation |
|---|---|---|---|---|
| **X-18.2** | DESIGN §18.2: "Runs entered by a switch carry the new label's own coordinates from their start" | Amend: except when the label's own lineage is still live; such runs are lower bounds, marked per C-N2 | The walker keeps only continuing chains there. Making it exact costs extra tracking (C-N2 b) | Amend, with the stated exception |
| **X-C8** | C8: the library requests coordinates for trace when supported | Under a request memory budget, auto-on only if the D4 gate passes (depth at the stop within 10%); otherwise off, and the tool result says how to ask | A run's coordinate account is 4–20× the run's own; budgeted fetches would stop much shallower on taxid columns | Adopt the gate |
| **X-C12** | C12 adds the K kind `coordinates`; DESIGN §14.1 and SPEC §7.5 freeze limitation kinds | Amend both texts: limitation kinds are an open `[a-z_]` token set that readers keep (precedent `memory_bound_soft`) | Keeps the frozen-contract text consistent with C12; no wire change | Adopt |
| **X-R3** | R3: 8 per dependency row, 1 per BRWT probe, 1 per selected coordinate | Probes counted only below a nonzero root rank; tuple lookups counted; selected coordinates = stored coordinates in the canonical window, cancelled ones included; a bounded deep term; lazy shifting so the price bounds the work | R3 taken literally prices selected reads above full reads on sparse-diff indexes. The review asks to count examined, shifted and cancelled coordinates | Adopt, and publish the measured ratio table |
| **X-B8a** | B8: /resolve fields "as proposed" | Add `clamped` (time budget only), `outcome.selection` and `resolved.work_units` (consumption on success); drop the planner's `discover.max_labels` clamp | Lets the ledger see consumption without a stop; keeps `--resolve-max-labels` to the scope you asked for | Adopt |
| **X-B8b** | B8 work formula from the Oct-3 proposal (8 per k-mer, entries and coordinates, dependency rows, "1 per run update") | 8 per in-graph k-mer + entries + coordinates + dependency units (+ the X-R3 terms for selected explicit reads) + a selection reservation per qualifying run; every decoded row is charged in `used`, including rows past x | Traverse's stage-2 F3 rule (charge what is decoded); a deterministic upper bound | Adopt |
| **X-B10** | B10: "an explicit `select.policy` that stops fails the request" (approved text: "a stop before an explicit interval's end fails the request") | Trigger as approved. Shape: (a) HTTP 200, prefix profile kept, `outcome.selection: "failed"`, no `selection`; (b) an error status | (a) keeps the useful prefix; (b) is the literal reading | **(a)** |

### 4.4 Kept as decided (despite proposals to change)
- **B3:** the production measurement is kept. It runs on the local UHGG and SRA BRWTs (3b-2 revision 2).
- **C2:** no outcome class for now. Revisit after the first external review of level 6.
- **R6:** no local refseq33m. Staging measurements are run by you after each deploy.

### 4.5 Ops notes (launch configuration, not contract)
- **O1:** set `--resolve-max-query-bp` explicitly on staging and production now. It exists at level 4 and defaults to 0. Unbudgeted /resolve has no time bound before level 9.
- **O2:** set `--traverse-max-memory-mb` on staging and production before coordinates are used with "unlimited". It defaults to 0.
- **O3:** set `--resolve-max-time-ms` only after the level-9 staging comparison (3b-1 revision 9).

---

## 5. Effort, calendar and workflow shape

### 5.1 Effort per increment (workflow hours, including internal review and verification)

| Stage | Increments | Hours | Elapsed |
|---|---|---|---|
| Coordinates | C0 0.5–1 · C1 7–9 · C2 5–6.5 · C3 3–4 · C4 5–6.5 · C5 8–9 · C6 3–4 · C7 2–2.5 · C8 4–5 | 38–48 | 29–36 (C5–C6 in parallel with C3, C4 and C7) |
| 3c-core | 3c.0 2–3 · 3c.1 5–6 · 3c.2 14–17 · 3c.3 10–12 · 3c.4 4–5 · 3c.6 10–12 | 45–55 | 43–52 (3c.4 overlaps 3c.6), plus 3–6 h wall for index builds in the background |
| 3c-ii | 3c.5a 8–10 · 3c.5b 4–6 · 3c.5c 3–4 | 15–20 | 15–20 |
| 3b-1 | 3b.0 1–2 · 3b.1 18–24 · 3b.2 7–9 | 26–35 | 26–35 |
| 3b-2 | 3b.3 13–17 · 3b.4 9–12 · 3b.5 3–5 | 25–34 | 18–26 (3b.3 starts during 3b.2 in a worktree) |
| External review fix passes | coordinates 6–10 · 3c-core 10–15 · 3c-ii 4–6 · 3b-1 6–10 · 3b-2 6–10 | 32–51 | 32–51 |
| **Total** | | **181–243** | **about 165–220** |

The stage plans estimated 32–39, 51 and 40–55 hours. The increases come from:
- the confirmed critique items;
- the external review loop, which no plan budgeted;
- the 3c split.

### 5.2 Workflow shape
- **One workflow per chunk of 1–3 increments.** Each one ends in verification. The workflows never commit. The main session commits after each chunk, so every committed state is verified.
- **Shape of each chunk** (the proven fix-batch pattern):
  - writers;
  - 2–4 read-only adversarial reviewers with lenses;
  - a fixer;
  - a verifier, which runs every suite with exact counts, byte identity against the previous level's binary (≥ 1,000 requests: mini_refseq, UHGG, SRA, /resolve), the determinism matrix, the canonical capabilities diff and Docker gcc:12.
- **Parallel writers only with disjoint file ownership in one tree** (C++ vs Python). Never two C++ writers in one tree: they share build directories.
- **Worktrees** only for:
  - reference binaries, as now;
  - 3b.3, the annotation-only formats, with their own build directory.

  At most one extra heavy build at a time, because the machine is shared.
- **Course corrections:** stop, edit and resume (no messages to running workflow agents).

**Chunks**

| Chunk | Content | Commit |
|---|---|---|
| W1 | C0–C2, C++ | "coordinates: walker and wire contract" (no level bump yet) |
| W2 | C3, C4, C7 (C++ writer) in parallel with C5–C6 (Python writer); then C8 | Level 6 |
| W3 | 3c.0–3c.2 | Decoder, not wired, no output change |
| W4 | 3c.3, 3c.4, 3c.6 | Level 7 |
| W5 | 3c.5a–c | Level 8 |
| W6 | 3b.0–3b.2 (3b.3 starts in a worktree during 3b.2) | Level 9 |
| W7 | 3b.3 (finish), 3b.4, 3b.5 | Level 10 |

- After each level: you deploy to staging, run the bench, and request the external review.
- The next chunk may start its library-only increments while a review is out.
- The review fix pass is applied before the next chunk's wiring increment.

### 5.3 Calendar
- **Assumptions:**
  - workflows run about 16 hours a day, unattended, including overnight;
  - your turnaround per level (deploy, staging bench, external review) is 1–2 days and overlaps the next chunk.
- **Pure workflow time:** about 10–14 days.
- **Milestones**, counted from the fix batch's commit:

| Milestone | Approximately |
|---|---|
| Level 6 (coordinates) | day 2–3 |
| Level 7 (3c-core) | day 5–7 |
| Level 8 (3c-ii) | day 7–9 |
| Level 9 (3b-1) | day 10–12 |
| Level 10 (3b-2) | day 12–15 |

- Each milestone includes the previous level's external review fix pass.
- **Total:** about 3–4 calendar weeks to level 10, including the review fix passes and the waits for your reviews.
- **The largest uncertainties:**
  - 3c.2: window and slack exactness, and the memory model against jemalloc;
  - 3c.5b: whether routed full reads stay as fast;
  - 3b.1: the /resolve charge model on wide rows.
