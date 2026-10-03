# Review request: stage 3 of the resource contract (budget-aware row-diff decoding)

**What I want from you:** an adversarial review of stage 3 of the resource contract
(`docs/DESIGN-traverse-graphlet.md` §14, §14.1, §19.1). Branch `gr/labeled-traversal`; the stage-3 commits follow
`d8d3ef86` (read their messages: each states what it claims). The spec carries the normative text:
`docs/SPEC-labeled-traversal-core.md` §6.8 (the head-admission and decode paragraphs), §7.0 (the `memory_bound_soft`
and `walk_domain` rows, `resource_stop`) and §7.5.2 (the Q record).

**The promise** (§14): the annotation reads used by traversal go through a decode path that charges each row-diff
dependency row and coordinate tuple as it is materialised. When the request's remaining budget would be exceeded,
it aborts the read without side effects, and the step stops with `resource_stop.phase = annotation_decode`. The
suggested actions are a more selective seed or a label-constrained query; a smaller radius would not help.

## What was built

1. **An opt-in decode path in the annotation library.**
   - `IRowDiff::decode_rows` / `decode_row_tuples` take a `DecodeBudget` (`src/annotation/binary_matrix/base/
     decode_budget.hpp`, `row_diff/row_diff_budgeted.cpp`). They are single-threaded and return a status instead of
     throwing. Throwing from the existing OpenMP loops kills the process, so the budget-aware path is a separate
     single-threaded version.
   - Every buffer is charged *before* it is allocated:
     - the trace to the anchors, through a small open-addressing table with modelled growth, replacing hopscotch,
       whose growth no model bounds;
     - each dependency row, read by a new depth-first single-row BRWT/ColumnMajor descent;
     - each coordinate tuple, sized from its delimiters before it is copied;
     - the reconstruction buffers. The XOR step reserves exactly; results are checked equal to `add_diff` on 9.3M
       mini_refseq and 3M UHGG rows.
   - A read that does not fit is refused whole.
   - The existing `get_rows` / `get_row_tuples` / `call_rows` / `get_rd_ids` and the BRWT slice code are unchanged:
     all 492 decode- and query-path functions compile to the same instructions as before.
2. **Per-row costs independent of batching.** Each returned row carries two costs:
   - its dependency units (8 per dependency row plus its stored entries), charged as work;
   - its *standalone demand*: what reading it alone, with its dependency rows, holds at most.

   Every key a fetch returns is admitted against its standalone demand, cached or not. A run of keys that does not
   fit is retried in halves, and only a key that does not fit alone stops the walk. So where a walk stops does not
   depend on `annotation.batch_kmers` or the lookahead.
3. **Consumers.** LabelQuery and LabelRecorder have budgeted fetches; the walker's level fetch, seed phase,
   derivation and annotate roots use them. The stops:

   | Situation | Result |
   |---|---|
   | A refused read in a level | stops the seed like a refused head at the level's first head, with phase `annotation_decode` |
   | A refused read in the seed phase or at an annotate root | fails the seed |
   | Annotate mode: the row fits but its new dictionary labels do not | phase `traversal`, stating the labels' number and bytes |
   | A level whose own key and successor lists leave nothing for its read | phase `traversal`, naming the lists |

4. **`memory_bound_soft` stays, and on budget-aware formats it now names only what is still uncharged:**
   - a level's key and successor lists, and its fetched rows, until its heads are processed;
   - the seed phase's intersection and hits;
   - a failed seed's depth-0 dictionary;
   - the index-wide header lookup (not observed);
   - the dictionary's first table and growth transient (about 1.5 KB; not observed);
   - a failed seed's echo.

   Other formats keep stage 2's statement.
5. **MGT v1 is unchanged.** The phase `annotation_decode` and the actions `more_selective_seed` and
   `label_constrained_query` are free `[a-z_]` tokens.

## The evidence

| Check | Result |
|---|---|
| Byte identity with budgets unset: 2,060 CLI requests (fixtures over 15 annotation formats, mini_refseq, UHGG, SRA), on HEAD and new, Release and Debug | 0 differences |
| Server routes (`/search`, `/align`, `/resolve`, `/traverse`, `/stats`, `/column_labels`) on mini_refseq and UHGG; concurrent vs sequential, 16 clients | 0 differences |
| Index building (every `transform_anno` type) and `metagraph query` | identical output |
| Default decode speed, Release, new/HEAD | 0.997–1.000 (UHGG query, single thread) |
| Allocator vs charged model (jemalloc thread peak) on SRA, UHGG and mini_refseq rows | ≤ 1.0 |
| Every charge refused in turn (about 11,000 charges over 15 row-diff cases; binary and coordinate; constrain, annotate, merge, derived, trace) | always a consistent prefix |
| `batch_kmers` 1/3/64/1000 sweeps | no dependence of the stop |
| Unit tests | traversal filter 299 (281 before), annotation suites 564 (557 before) |
| Integration `*Traverse*` | 50 |
| Python | 577 OK, 11 skipped |

For the consistent prefixes: a refused read gave walks equal to the unrefused walk up to `complete_to_bp`, valid
bins, a codec round trip and a stated stop.

The last verification found one defect: under a memory budget a level stop's *message* depended on `batch_kmers`
(13 of 408 UHGG cases). Whether the refused row's exact demand was known depended on what the lookahead had read.
Fixed: a level stop now states only the batch-independent fact, that the row needed more than the bytes left.

**Expected consequence.** Budgeted walks are never deeper than before, and some stop earlier: on memory-only
budgets the median ratio is 1.00, the 10th percentile 0.55–1.00, and a few go to 0 bp. Before, the decoder's
buffers sat outside the account and were reported only as soft excess; now they are charged.

## Questions

1. **Work weight per dependency row (open owner decision).** Each returned row is charged its full dependency path
   (8 units per dependency row plus 1 per stored entry), so that charges do not depend on batching. Shared
   dependencies are therefore counted many times: 10–50× the earlier units for the same walk. Options:
   - keep it (conservative);
   - (a) 1 unit per dependency row, no per-entry charge;
   - (b) each distinct dependency row once per seed, with the per-seed set charged to memory.

   Which is right for a work bound that should track decode effort?
2. **Lower-bound statements.** Some depth-0 failures from unread roots or unnamed labels state "need at least
   X MiB", a lower bound rather than the exact need. Raising the budget to X may still not suffice. Acceptable?
3. **Statement texts.** §14.1 lists effect templates as frozen, while K stores free text. Are the new texts of
   the label stop, the level-lists stop, the root and seed failures, and the work wording ("cut W units past the
   budget") acceptable as free text?
4. **Capabilities.** `/traverse/capabilities` `work_bound` does not mention the dependency weights. Changing it
   would change the capabilities document. Should it?
5. **Stage 3b.** Formats without row-diff (BRWT, ColumnMajor, TupleCSC alone), the disk and flat row formats,
   IntRowDiff, and `/resolve` still decode without a budget and keep stage 2's `memory_bound_soft`. A level's rows
   and lists are stated, not charged. Is that the right boundary for this stage?

## Known and queued — do not re-report

- **Your recheck of `278a53dd..d8d3ef86`** (P1–P7) is queued for the next pass, together with your design answers:
  - P1 float widths in delivery pricing;
  - P2 seed fetch charging and failed-seed `largest_charge`;
  - P2 trace-validation coordinate copies;
  - P3 the spooled receipt and multi-walk `left_out`;
  - P3 the continuation-equivalence wording;
  - the continuation branch budget, extra-label reachability and the `seed_id` wording.

  Stage 3 changed the seed-phase code that P2 touches, so P2 is fixed on top of this stage.
- **Queued from the owner:**
  - stop a traversal when the client connection is gone;
  - `request_id` with `POST /traverse/cancel`.

  A request-level total budget is deferred.
- **Pending an owner decision:** per-graph identity (manifests) for multi-graph servers.
- **Not part of this stage:** stage 4 (the request ledger) and record coordinates (§18).

## How to verify

```bash
cd metagraph/build && make -j metagraph unit_tests       # AppleClang: configure with -DCMAKE_CXX_FLAGS=-Wno-error=unused-result
./unit_tests --gtest_filter='*Walker*:*LabelOracle*:*Resolve*:*Trie*:MiniRefSeq*:*Graphlet*:Stage2Review*:*CoordToHeader*:*HeaderIndex*'   # 299
./unit_tests --gtest_filter='*RowDiff*:*BRWT*:*TupleRowDiff*:*IntRowDiff*:*Annotat*'           # 564
test_venv/bin/python ../integration_tests/main.py --test_filter='*Traverse*'                    # 50
cd .. && python3 -m unittest discover -s api/python/tests -p 'test_traverse_*.py'               # 577 (11 skipped)
python3 scripts/traversal/graphlet_fixtures.py --check && python3 scripts/traversal/graphlet_fixtures.py --from-cli build/metagraph --check
python3 scripts/traversal/make_codec_vectors.py --check && python3 docs/source/_graphlet_examples/run_examples.py
```

The new tests:
- C++: `RowDiffBudgetedDecode.*`, `LabelOracleBudgeted*`, `GraphletStage3Decode.*`, `GraphletStage3Review.*`,
  `MiniRefSeq.BudgetAwareReadsDoNotDependOnBatchKmers`.
- Integration: `test_stage3_*`.
- Python: `test_traverse_stage3_decode.py`, with the `budgets/decode_stop` fixture.

## How to report

For each finding, give:
- severity;
- whether it changes the MGT v1 format, or unbudgeted output;
- `file:line`;
- a concrete failure;
- your proposed fix.

Mark each one as confirmed, suspected or design disagreement. End with a verdict (stage 3 sound, sound after
named fixes, or not sound) and an answer to each question above.
