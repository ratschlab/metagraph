# Review request: recheck of the fixes to your stage-2 and Part-B reviews

**What I want from you:** a recheck of the fixes to your last two reviews of `gr/labeled-traversal`:
- the delayed review of `4d332729` plus its working tree (findings 1–8);
- the review of clean `278a53dd` (findings N1–N4 and your answers to the design questions).

The fixes are the commits after `278a53dd`; read their messages. `docs/DESIGN-traverse-graphlet.md` v5.3 §19 records the resource contract as implemented and the library decisions we adopted from your answers. The spec carries the normative text: §5, §6.8, §7.0, §7.1 and §10.3.

Please try to break each fix, and say for each finding whether it is fixed, partially fixed or not fixed. **MGT v1 is unchanged.** Without budgets, server output is byte-identical to `278a53dd`: 3,148 runs over two random corpora and mini_refseq, in all four details, on Debug and Release, after normalising elapsed time.

## Your findings and what we did

| Finding | Fix | Reproducer, before → after |
|---|---|---|
| 1 escaping beyond the delivery reserve | Every per-object delivery cost is a derived worst case: jsoncpp's tree plus three copies of the text, or four copies of a graphlet body. Names are charged at their real escaped length, once per place they are written. A seed whose depth-0 state cannot hold that fails at depth 0 with the size it needs. | control: 2,163,840 B delivered under 2 MiB with soft 0 → seed failed at depth 0 ("need 8 MiB"), response 3.9 KB. percent and percent1: same shape. |
| 2 work overrun beyond W | Work is charged when a fetch returns a row: 8 units per key, 1 per entry and per coordinate. It is compared after every charge. Fetch calls are sized from the budget left and start from one key per level. **Every work stop states the most its seed charged between two comparisons (`largest_charge`).** The promise is that number, no longer W. | work: used 125,044 → 25,008, with largest charge 25,008 stated. |
| 3 interrupted fetch hides memory | What a fetch holds beyond the account is recorded before any check after it can stop the walk. This covers every fetch call, the roots' rows, the lookahead's warming, the seed phase's label validation and the depth-0 dictionary. | soft: observed 0 → 744 MiB, a model price that errs high. |
| 4 continuation restores spent loss | The loss budget is reduced by the **largest** terminal loss among the continued labels and stated as `{original, effective, largest_terminal_loss}`. | continuation.py: next budget 3 → 1. The continued walk ends at 120 bp against 150 uninterrupted, which is conservative. |
| 5 switched continuation invalid | `labels.extra` is rebuilt around the new seed labels. Continuations that cannot share one request raise `IncompatibleContinuations`, with `next_requests()` as the fallback. | duplicate-extra.py: backend 400 → accepted. |
| 6 later merges cause false differences | The DAG is restricted to the comparison depth. Merges at or beyond the depth are left out and named; where the cut cannot be exact, the result is `qualified`. | repro.py: equal False → True in all four modes. |
| 7 header-cache fingerprint | The header index belongs to the `CoordToHeader` object (`find_header`). It is built on first lookup, dropped by `load()`, and never carried to copies. The address cache is deleted. | probe_header_cache: X → seq 1 (named Y) → X → seq 2. |
| 8 failed derivations omit `memory_bound_soft` | It is stated on every failure path under a memory budget, with the value carried in `SeedDerivationError`. | probe_failed_memory: absent → present. |
| N1 false time stop at the radius | A stop is created only when a head below the radius was censored. | radius_zero_request.json: `resource_stop`/Q → none, outcome complete. |
| N2 binding errors exceed `max_bytes` | The supplied ceiling is resolved before the error is fitted. | 77 B and 74 B → 64 B. |
| N3 live oracle ignores `index_fp` | Available manifests are compared; a missing manifest stays unverifiable. | `live_matches_cache` True → False. |
| N4 spec §7.1 on `route_bp` | Spec §7.1, the MGT R legend and the `LabelRun::route_bp` comment now say: R = the earliest non-first-parent merge, T = the latest entry, displayed support = `[evidence_from, to_bp)`. A run closed by a merge of three or more parents may carry `route_bp == to_bp`. | text only; route-stamp.py unchanged (61 vs 113). |

**Your design answers, as implemented** (design §19.2):
- Annotate claims: union rule, witness chosen by stored parent order, documented as one witness.
- R/T/`evidence_from` as you proposed.
- The U+FFFD refusal is kept.
- A one-sided manifest gives `index_unverifiable` ("not proven different"), with the opt-in `allow_unverified_index` (default off).
- `merge_entered` is used in the walks filter as well as the claims filter.
- Successful export/save/load always keep their receipt: optional fields are dropped first.
- Local operations are documented as having no work or allocation budget; listed as stage L in §14.1.

**Part A answers, as implemented:**
- The radius exemption is kept.
- An injected denial says it was injected, and its only action is `continue_from_leaves`.
- The spec states that a memory bound is not a wall-clock bound, and that delivery is not deadline-checked.

## What our own adversarial review found after the first fixes (all fixed)

Three internal reviewers attacked the first fixes. They confirmed 11 more problems:

- **F1** — failed results echoed index names without charging them. Now:
  - under a memory budget, an echoed header is cut to 256 bytes plus `(N bytes; column C, sequence S)`;
  - the request's own `seed_id` is not cut; its excess is stated in `memory_bound_soft`.
- **F2** — the label validation failed the seed before recording what it held.
- **F3** — rows that were decoded but never consumed went uncharged. We kept your definition, "rows decoded" (design §14), instead of counting only consumed rows.
- **F4** — both roots were charged before the first comparison. This is now stated as one charge, by design: a stop between the two roots would deliver no result complete to 0 bp.
- **F5** — trace coordinates were charged as one unchecked sum. They are now part of the row's charge.
- **F6** — the depth-0 dictionary was not recorded before its admission.
- **F7** — some charges had no stated size. `largest_charge` now covers every kind.
- **A** — a merge exactly at the comparison depth caused false differences.
- **B** — malformed continuation overrides raised exceptions instead of returning errors.
- **C** — continuation receipts for walks with many labels were refused at the default ceiling.
- **D** — labels alive at the leaf but not seeded were dropped silently. They are now named in `left_out`.

## Questions

1. **Work bound.** The guarantee is now "overrun ≤ the stated largest single charge", and on a wide index that number can far exceed W: 70,000 records on one k-mer give 70,008 units, and with both arms 140,016. Until stage 3 makes the decoder interruptible, is a stated bound acceptable as the contract? Is "rows decoded, charged at fetch" the right definition?
2. **Delivery reserve.** The fixed part of the delivery model is about 450 KiB: zlib state, the envelope, and the maximum number of limitations. As a result, a 5,000-character control-character header already fails at 1 MiB in detail full. That is conservative but coarse. Should the fixed part shrink, for example by pricing limitations as they arise?
3. **Continuations.** Should a continuation also reduce `max_label_branches` by the branches already used, the way it now reduces the loss budget?
4. **Extra labels.** Should the server enforce a one-switch reachability rule for extra labels? Today the library keeps every extra label the original request named.
5. **`seed_id` echo.** It is stated in `memory_bound_soft`, not cut. Is that the right choice, given that it is the caller's own input?

## Known limitations — do not re-report

- **Stage 3** (the budget-aware, interruptible row-diff decoder) and **stage 4** (the request-level ledger) are not built. `memory_bound_soft` stays stated until stage 3 lands.
- Rows the lookahead decodes without a fetch asking for them are not charged. Annotate header rows are charged by their distinct labels, not by their coordinates. Spec §6.8 says both.
- `/resolve` and `/select` still pass label names through the JSON writer unchanged.
- No 3.10 interpreter was available; the Python suite ran on 3.14 only.
- Record coordinates (design §18) are not implemented.

## How to verify

```bash
cd metagraph/build && make -j metagraph unit_tests       # AppleClang: configure with -DCMAKE_CXX_FLAGS=-Wno-error=unused-result
./unit_tests --gtest_filter='*Walker*:*LabelOracle*:*Resolve*:*Trie*:MiniRefSeq*:*Graphlet*:Stage2Review*:*CoordToHeader*:*HeaderIndex*'   # 281 tests
test_venv/bin/python ../integration_tests/main.py --test_filter='*Traverse*'                    # 46 tests
cd .. && python3 -m unittest discover -s api/python/tests -p 'test_traverse_*.py'               # 571 tests (11 skipped)
python3 scripts/traversal/graphlet_fixtures.py --check && python3 scripts/traversal/graphlet_fixtures.py --from-cli build/metagraph --check
python3 scripts/traversal/make_codec_vectors.py --check && python3 docs/source/_graphlet_examples/run_examples.py
python3 /tmp/metagraph-stage2-resource-repros.py --case all --binary build/metagraph           # your resource cases
```

The new tests:
- C++: `Stage2ReviewRound3.*`, `LabelOracleHeaderIndex.IsBoundToItsCoordToHeader`, the extended `Graphlet.DeliveryCostsBoundTheOutput`, and `WorkStopsStateTheirLargestCharge`.
- Integration: `test_stage2_*`.
- Python: `api/python/tests/test_traverse_review3.py`, with committed responses in `api/python/tests/data/traverse/review3/`.

## How to report

For each finding, give:
- severity;
- whether it changes the MGT v1 format, or unbudgeted output;
- `file:line`;
- a concrete failure;
- your proposed fix.

Mark each one as confirmed, suspected or design disagreement. End with a verdict: stage 2 and the library fixes sound, sound after named fixes, or not sound. Give an answer to each question above.
