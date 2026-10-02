# Review request (round 2): label-consistent graph traversal in MetaGraph

**What I want from you:** a second adversarial review of the traversal feature on branch `gr/labeled-traversal`,
now **committed and pushed** (11 commits on top of upstream `1de02cec`, head `01a5230a`). Round 1 found a
blocker (trace support fabricated direct support across a merge), a crash and eleven further defects; all are
fixed and listed in §6 so you do not re-report them. Since then the feature grew a verification layer and a
number of real-index results, and that is what this round is about.

**Priorities, in order:**

1. **Correctness of the verification itself.** The new `annotate` mode produces an exhaustive, label-free trie
   that is used as an *oracle* for the label-constrained walk, and a string-level model of the records is used as
   a third party. If the oracle, the model, or the way they are compared is wrong, every green test built on them
   is worthless. Attack that first (§4 C13–C21, §5.1).
2. **Correctness of the new surface** — `exhaustive`, `complete_to_bp`, the per-node label cap, derived seed
   labels (`per_hit`), switching with a budget (§5.2).
3. **Efficiency, with the measurements I have** (§5.3) — what decides whether this works on `refseq33m`
   (33 M labels, 627 G nodes, coordinate annotation), and whether the three complexity fixes are exact.
4. **Evidence honesty**, again (§5.4) — any output that claims more biological support than the evidence warrants.
5. **Design critique and API** (§5.5–§5.6).

Please be blunt. If the verification approach is circular, say so and say what would not be.

---

## 1. What the feature does (one page; unchanged in substance from round 1)

MetaGraph indexes DNA as a de Bruijn graph (nodes = k-mers) plus an **annotation matrix** mapping each k-mer to
the **labels** (samples, genomes, accessions…) that contain it. Search answers *"which labels contain my query?"*;
it cannot answer *"what sequence lies **around** my query in those labels?"*. This feature adds that:

```
query sequence
  └─ resolve  → which labels support which contiguous k-mer intervals of the query   (no traversal)
       └─ select  → freeze an explicit seed list (sequence + label set), policy recorded
            └─ traverse → extend each seed outward through the graph, constrained by labels
```

The traversal walks outward from the seed and at every step keeps only the labels that still support the path;
a step is admissible while at least one permitted label survives. The state is a **set** of labels (one seed
carried by 1,000 genomes is one walk whose set splits where the genomes diverge), with a **label-change cost**
plus a per-arm **loss budget** (`forbid` = one label must carry the whole path; `constant 1` with budget 1 = at
most one switch), per-**lineage** branch limits and quorums, per-path edge reuse, hairpin skipping and seed
re-entry as the structural rule, and optional reconvergence merging into a DAG.

**New since round 1:**

- **Derived seed labels (`per_hit`).** `seeds[].labels` may be omitted: the permitted set is then the labels that
  carry every k-mer of the seed (the intersection of the rows seed validation reads anyway — no extra I/O), capped
  at `max_seed_labels` and reported as such. This is the design note's default and the user's primary use:
  "traverse consistently labelled continuations from how the match is labelled, without naming labels".
- **`labels.mode: annotate`.** Admissibility becomes purely structural (edge reuse, hairpins, seed re-entry);
  the labels present at every node are *recorded* (full rows, a `LabelRecorder`) instead of filtered. No permitted
  set, no label state, no switching, no quorum, no branch limit — the code path shares the level loop and
  `check_structure()` with the constrained mode and **nothing** of the label machinery.
- **`exhaustive: true`.** Pins the knobs that would prune (unlimited branches and splits, `on_reconverge: keep`,
  no quorum, `on_overflow: stop`) and *refuses* a request that sets them otherwise, so "the exhaustive trie" cannot
  silently be a beam-pruned one.
- **`complete_to_bp`.** Exploration is level-synchronous; a size cap trips between two heads of one level; that
  level is marked partial and excluded. Every arm reports the depth up to which *every* admissible walk is
  present, and `status: complete ⇔ complete_to_bp == max_extension_bp`. The walk rule the claim quantifies over
  is echoed as `walk_rule`.
- **Per-node label cap.** Recorded lists are cut at `labels.max_labels_per_node` (default 64) but every list
  carries its true count, nodes_truncated is reported per arm, and the oracle refuses to prove anything with cut
  lists.
- **Verification layer** (§3): the oracle contract, a records model, 25 edge/non-edge cases, and the contract on
  a real index.

Design note (project owner): `docs/DESIGN-labeled-traversal-endpoint.md`. Implementation spec (mine):
`docs/SPEC-labeled-traversal-core.md` (1105 lines; §6.9 the oracle, §6.10 the completeness boundary, §11.3 the
test list T1–T36). **Where code and spec disagree, that is a finding.**

---

## 2. How to see the code

```bash
git fetch origin gr/labeled-traversal && git checkout gr/labeled-traversal
git log --oneline 1de02cec..HEAD        # 11 commits; read the messages, they state what each one claims
git diff --stat 1de02cec HEAD
```

| Path | Lines | Role |
|---|---:|---|
| `src/graph/traversal/traversal_types.{hpp,cpp}` | 95 / 74 | enums: regimes, label kinds, support kinds, end reasons |
| `src/graph/traversal/label_oracle.{hpp,cpp}` | 275 / 690 | per-request context: regime, oriented node → row key, batched membership (`LabelQuery`), full-row recording (`LabelRecorder`) |
| `src/graph/traversal/resolve.{hpp,cpp}` | 155 / 579 | `resolve` (support profile) and `select` (seed freezing) |
| `src/graph/traversal/walker.{hpp,cpp}` | 549 / 3037 | the traversal — **the core**; `walker.cpp:48–100` is the implementation-notes block |
| `src/cli/traverse.{hpp,cpp}` | 105 / 1229 | strict JSON parsing/validation, serialisation, CLI entry |
| `tests/graph/traversal/test_trie_oracle.hpp` | 667 | **the oracle contract**: E from recorded sets vs claims(A); tuned-run subset; merged-route soundness |
| `tests/graph/traversal/test_trie_reference.hpp` | 586 | **the records model**: walk rule over strings; per-label walks; §6.3 recurrence; trace |
| `tests/graph/traversal/test_trie_checks.hpp` | 485 | the three-way checks every fixture runs |
| `tests/graph/traversal/test_trie.cpp` | 1298 | annotate/exhaustive smoke tests, cap boundary, cut lists, `TrieOracle` suite |
| `tests/graph/traversal/test_trie_cases.cpp` | 986 | 25 edge/non-edge cases, caps sweep, random graphs |
| `tests/graph/traversal/test_walker.cpp` | 2517 | traversal semantics, typed over graph × annotation × mode |
| `tests/graph/traversal/test_label_oracle.cpp` | 485 | regime/key mapping, access paths, eviction |
| `tests/graph/traversal/test_resolve.cpp` | 409 | support profiles, selection policies, scale |
| `tests/graph/traversal/test_mini_refseq.cpp` | 723 | **against a real index** in the `refseq33m` format, incl. the contract vs the 42 source records |
| `integration_tests/test_traverse.py` | 935 | CLI + HTTP end to end (23 tests) |
| `scripts/traversal/build_mini_refseq.sh` | — | builds the real-format fixture from 42 NCBI accessions (~25 s) |

Modified pre-existing files (173 lines): `src/cli/server.cpp` (+106: routes, shard addressing),
`src/cli/config/config.{hpp,cpp}` (+60), `src/main.cpp` (+4), `integration_tests/main.py` (6).

**Suggested reading order for this round:** spec §6.9–§6.10 and §7.4 → `test_trie_oracle.hpp` (read the header
comment: *what* is compared and why claims, not leaves) → `test_trie_reference.hpp` (the model; is it the rule?)
→ `test_trie_checks.hpp` → `walker.cpp` `process_item_annotate`, `check_structure`, `cap_check`, `trip`,
`stop_frontier`, `bounded`, `finalize` (the `complete_to_bp` computation) and `validate_strategy` →
`test_trie_cases.cpp` → `test_mini_refseq.cpp` (last test).

---

## 3. The verification layer — what it claims to establish

For a seed S, radius R and permitted set P, all under the `exhaustive` preset:

```
T = trie(S, R, annotate)                 structural, labels recorded, no label logic
E = { (w, l) : l ∈ P, w a walk of T on every node of which l is present, maximal per label }
A = trie(S, R, constrain, P, forbid)     the label-constrained walker
contract:  claims(A) == E                a claim is a (walk, label) pair — a label's END position, so a label
                                         ending inside a walk other labels continue is checked at its own end
```

plus, for any tuned run (branch limits, quorums, beam, caps): every tuned leaf is a prefix of an exhaustive leaf
under each of its labels, and every (walk, label) of E the tuned run does not reach has a *recorded* reason where
the run left it (a branch event naming the base not taken, a quorum text, a cap, the beam). Nothing drops
silently.

**The third party.** T and A share the graph, the annotation and `check_structure()`, so their agreement cannot
catch an error common to both. `test_trie_reference.hpp` therefore models the index from the records a fixture is
built from and nothing else: per-label packed k-mer sets, and the walk rule of §6.10 re-implemented over strings —
the structural walks with the labels present at every node (compared with T's recorded sets and path reasons),
each label's maximal walks with the reason each ends (compared with claims(A)), the §6.3 recurrence for a constant
switch cost and a budget (compared with constant-cost runs: leaves, losses, switch events), and the records' own
continuations for `support: trace` (compared with trace runs). Every fixture is also checked for: the run over P
equals the union of the single-label runs; a permitted set derived from the seed equals the explicit list; the
recurrence at budget 0 equals the forbid leaves.

**On a real index** (`MiniRefSeq.TrieContractAgainstTheSourceRecords`): the model is built from the 42 source
FASTA records (9.3 Mbp, 2.8 s) of the fixture — k = 31, basic, **unmasked** `DBGSuccinct` (so `'$'` dummy
k-mers are live), `RowDiff<BRWT>` with coordinates, header (accession) labels, exactly the `refseq33m` shape — and
the contract holds on the whole blaNDM-1 seed (19 carriers; 24 structural walks left / 2 right at 300 bp) and on
three 150 bp windows cut from records, under k-mer and under trace support, at 300 bp (all parties) and 1000 bp
(claims vs records).

**What this does and does not establish.** It shows the constrained walker returns exactly the label-consistent
walks the records admit *under the stated walk rule*, with the stated reasons. It does **not** show the walk rule
is the right finite set of walks to return (see §5.5.1 — the homopolymer consequence is the clearest example), and
it does not cover `on_reconverge: merge` except through the route-soundness check.

---

## 4. The claims to attack

C1–C12 from round 1 still stand (restated briefly; C1 and C3 were qualified after round 1). C13–C21 are new.

| # | Claim | Where |
|---|---|---|
| C1 | For a fixed path with no branch limit and `constant` cost, the per-label loss is the exact minimum over label assignments under *switch_on: loss* (a switch only from a label lost at that step); under a `table` cost exactness is subject to `max_switch_sources` (the preset forces it unlimited). | spec §6.3; `derive()` |
| C2 | Branching is counted per lineage; a renamed lineage cannot evade its branch limit. | `process_item` exclusion loop |
| C3 | Under `forbid` with `on_reconverge: keep`, the result over P equals the union of the single-label walks. | `check_single_label_union` on every fixture |
| C4 | Each (k+1)-mer edge is used at most once **per path** (canonical in canonical regimes, orientation kept). | `check_structure`, `edge_key` |
| C5 | Block precedence `rejoined_seed > edge_reuse_rc > edge_reuse`; seed re-entry keyed on the target node (either strand). | `block_rank` |
| C6 | A hairpin step (self-RC (k+1)-mer; for even k also into/out of a self-RC node) is skipped by default and never counts toward ambiguity. | `is_hairpin` |
| C7 | `trace` requires consecutive coordinates (+1 right, −1 left) and is only defined for BASIC graphs; `trace` + `merge` is rejected. | `validate_seed`, `derive` |
| C8 | Results are invariant under `batch_kmers` in both modes (annotate assigns label ids at consumption). | `AnnotateIsDeterministicAndBatchInvariant` |
| C9 | `result` is byte-identical across runs and request splits; per-seed annotation counters are deltas. | `Determinism` |
| C10 | `complete` ⇔ the requested domain was fully explored. | — |
| C11 | Strict JSON validation; nothing silently weakened (`annotation.access` must be `"auto"`). | `Strict` |
| C12 | One physical shard per traversal; ambiguity is a routing 400, not a parse error. | `resolve_traverse_index` |
| **C13** | **Annotate mode shares no label logic with constrain mode**: `process_item_annotate` never touches `derive`, `Entry`, `State`, the cost, the budget, quorums or branch limits; the only shared pieces are the level loop, `check_structure`, the caps, merging and the segment DAG. If any label decision leaks in, the oracle is circular. | `walker.cpp` |
| **C14** | **`complete_to_bp` is sound**: for every n ≤ `complete_to_bp`, every walk of n bases from the seed boundary obeying the walk rule is in `segments`, and every path ending before it ended for its reported semantic reason. The partial level is marked, not rolled back, and excluded. Heads already at the radius when a seed-level cap trips are ended `max_extension_bp`, not censored. | `cap_check`, `trip`, `stop_frontier`, `finalize`; `TrieCasesCaps`, `TrippedCapReportsTheCompleteDepth` |
| **C15** | Under the `exhaustive` preset nothing prunes except the level boundary; the preset refuses conflicting knobs naming the required value, including `max_switch_sources` under a `table` cost, and a derived set cut by `max_seed_labels`. | `validate_strategy`, `ExhaustiveRejectsConflictingKnobs`, `ExhaustiveRefusesACutDerivedSet` |
| **C16** | A derived permitted set is exactly the labels carrying every seed k-mer, in `(column, seq_id)` order, costs no extra annotation reads, is validated in O(sum of hits) not O(L²·M), honours the time budget, and yields results identical to naming those labels. | `validate_seed`; `check_derived` on every fixture; `ManyDerivedLabelsMatchTheExplicitList` |
| **C17** | Every recorded label list carries its true count (`labels_total` / `labels_distinct`), including blocked/hairpin events and the root; a cut anywhere increments `nodes_truncated`; live-label counts over cut lists are stamped inexact. | `bounded`, `distinct_labels`; `CutRecordedListsAreReportedAndBreakTheOracle`, `CutEventListsAndTheRootTotalAreReported` |
| **C18** | The oracle E in `test_trie_oracle.hpp` reads only segment sequences, the recorded sets and the path/segment structure, and compares **claims** (label ends), so that a label ending inside a continuing walk is checked at its own end. | `expected_leaves`, `constrained_claims`, `OracleComparesPerLabelClaimsNotLeaves` |
| **C19** | The records model implements the walk rule exactly, and "a label's maximal walks" (per-label DFS over its k-mers) equals "the maximal per-label prefixes of the structural walks" (what E computes from T). | `test_trie_reference.hpp` `WalkRule::steps`, `detail::dfs`, `label_trie` |
| **C20** | `hits_from_tuples` is O(coordinates) per node (a slot map, not a scan of the hits built so far); `max_support` selection is O(R log R) per round with a single materialisation; `LabelQuery::fetch` keeps a call's whole working set resident across eviction. | `label_oracle.cpp`, `resolve.cpp`; `MaxSupportScalesWithManyLabels`, `EvictionKeepsTheCurrentBatchAnswerable` |
| **C21** | The trace model's end reasons: at a record's end, `dead_end` when the node has no successor, `record_end` when a successor carries the label with non-continuing coordinates, `label_lost` otherwise; structural blocks still apply under trace. | `trace_trie` vs `TrieCasesTrace.*` |

---

## 5. Specific questions

### 5.1 Is the verification circular, and is the model the rule?

1. **C13.** Read `process_item_annotate` and everything it calls. Does any label-derived quantity influence which
   successor is followed, which path ends, or where? (`labels_at`, `record_present`, `bounded` are meant to be
   write-only.) Is the `continuation` of an annotate path, whose labels are "recorded on every node of the tail",
   a label decision in disguise?
2. **C18/C19.** `expected_leaves` walks each structural path and, per permitted label, takes the longest prefix on
   which the label is present at every node, then keeps the maximal prefixes per label by a lexicographic
   successor test. (a) Is "maximal per label" exactly the walker's label end? Think of a label present on a walk's
   nodes 0..j, absent at j+1, present again at j+2 on the *same* walk. (b) The lexicographic test claims "every
   extension of s directly follows s"; is that true with the empty walk and with walks that are prefixes of
   several others? (c) The records model's `label_trie` is a per-label DFS with the per-path used-edge set; I
   argue it equals the maximal-prefix definition because the prefix shares the path and hence the used set. Is
   there a case where a label's own DFS takes a step the structural DFS also takes but *blocks differently*?
3. **The walk rule in the model** (`WalkRule::steps`): hairpin → skipped; seed node → `rejoined_seed`; used key →
   `edge_reuse_rc` if the orientation differs else `edge_reuse`; `stop_reason` picks the strongest block among the
   label's passing candidates, else `label_lost` if none passes, else `dead_end`. Compare with
   `Walker::check_structure` and the end-reason selection in `process_item` (`walker.cpp` ~2300–2380). Where do
   they differ? In particular: a label whose only continuation is a *blocked* step while *another* label's
   continuation is followed; a candidate that is both a seed node and a used edge; two blocked candidates of
   different ranks.
4. **Recorded sets vs the model's `labels_at`.** On canonical and primary graphs the model inserts both
   orientations of every record k-mer and the recorder reads the row of the base node; on PRIMARY the row key is
   `get_base_node`. The three-way agreement held in all three modes on every fixture — but is there a k-mer for
   which the two disagree *and* no fixture exercises it (e.g. a k-mer whose RC is also in a different record)?
5. **What the real-index test proves.** The model for `mini_refseq` is built from the FASTA the index was built
   from, by my own builder script. If the builder drops or alters sequence (N handling, line endings, lowercase),
   the model and the index would agree on the wrong thing. The records are ACGT-only and the index builder is
   `scripts/traversal/build_mini_refseq.sh`; please check that the records the test loads are the ones the index
   was built from and that no k-mer could be silently absent from one side.
6. **`check_tuned_subset`.** Clause 2 follows each exhaustive claim through the tuned trie and accepts a `cut`
   reason (branch limit, loss budget, beam, cap) or a recorded branch event where the run leaves the claim. Is
   "a branch event naming the base not taken" too permissive — could a run drop a walk for a *wrong* reason that
   happens to be recorded?
7. **A cautionary tale about the JSON-level oracle.** The integration tests' `_walks` / `_structural_walks`
   compare *leaves* on the *right* arm. When I ported them to a script for the SRA run I (a) compared leaves to
   per-label claims and (b) reversed the left arm's flank as one string instead of per segment (a left-arm
   segment stores its bases in natural orientation, so walking order is `reverse(seg)` concatenated root→leaf —
   see `spell_path`). Both mistakes produce a "DIFFER" that looks like a walker bug and is not; with claims and
   the right orientation the SRA result is 114 = 114 and 2 = 2. The C++ oracle does this correctly, but the two
   conventions (claims not leaves; per-segment reversal) are easy to get wrong in any new consumer — is the
   output format itself inviting this, and should the left arm be serialised in walking order instead?

### 5.2 Correctness of the new surface

1. **C14, the boundary.** `trip()` sets `boundary` to the first unexpanded head's depth; the beam sets it to the
   pruned heads' depth; `MAX_STEPS` (a seed-level cap) stops the *other* arm at its own frontier depth. With both
   arms alternating by depth, can arm B's boundary be set by a cap that tripped in arm A at a depth B has already
   completed? Can `complete_to_bp` ever be *larger* than the last complete level? (The sweep in `TrieCasesCaps`
   found no case, but it is one fixture.)
2. **Partial level marked, not rolled back.** Children created before the cap tripped stay in the output ended
   with the cap reason, above `complete_to_bp`. Is anything downstream (`label_summary.reach_bp`, `direct_bp`,
   `continuation`) computed over them in a way that reads as evidence?
3. **The per-node cap and E.** With lists cut, `expected_leaves` sets `cut` and the test refuses to compare. But
   `constrain` mode also uses `bounded()` for its *reported* lists while its *decisions* use the full state. Is
   there any place where a cut list feeds a decision?
4. **Derived labels (C16).** The derivation intersects hits across the seed's k-mers starting from the cheapest row
   of the first batch. Under `trace`, the derived set must also be coordinate-consecutive over the seed. Check the
   compaction of per-k-mer coordinates (`TraceCoordinatesSurviveCompaction`) and the ambiguous-header refusal: can
   an accession present in two columns be derived under `column` kind and traverse under the wrong column?
5. **Switching (C1).** `OneSwitch` and `TwoSwitches` pin the recurrence against the records model: a target enters
   from the cheapest *lost* source at source loss + cost within the budget, stays when cheaper or equal, and a
   lost label whose only switch would exceed the budget ends `loss_budget` with the needed value. Is "switch only
   from a lost source" (`switch_on: loss`) what the spec says, and does `switch_on: any` have any reference at
   all? (It does not; is that acceptable for a knob that exists?)
6. **Even k.** The hairpin rule for even k treats a self-RC *node* as a hairpin in itself, so a record containing a
   palindromic k-mer cannot be reconstructed through it (`EvenKPalindromicNode` pins the walk stopping one base
   short). Is that the right rule, or should a palindromic node be enterable and only its *RC retrace* blocked?
7. **`switch_on: loss` is per successor, and that makes one switch + branch limit 0 useless.** Measured on the
   SRA index (§5.3, 24 random positions, 10 kb): `constant 1` / budget 1 with the default limit 0 reached a
   median of **57 bp** against 224 bp under `forbid`. Mechanism: at a plain divergence (A continues on P, B on Q)
   A is *lost on Q*, so it is an eligible source there and the lineage of A continues on Q renamed to B while it
   also continues on P as A — an ambiguous branch (C2, correctly), so limit 0 ends both labels at every
   divergence between two carriers. The user's intent ("exact continuation within the same label, allow one
   switch") is a switch **only when the label can continue on no successor at all**. Should `switch_on: loss` be
   defined globally per step (lost on every successor), with the current behaviour available as a third value?
   I believe this is a semantics bug against the design note, not a tuning matter.

### 5.3 Efficiency — measurements and the open problem

Measured, all single-threaded on a laptop (Apple M-series, Release build):

| What | Index | Numbers |
|---|---|---|
| 3 kb two-arm walk, 19 labels, merge on | mini_refseq (9 Mbp, RowDiff<BRWT>+coords, 9 columns) | 0.14 s incl. resolve; 24,877 rows requested (~1 row per step per arm) |
| exhaustive constrained, 19 labels, R = 300 | mini_refseq | 0.33 s both arms (26 leaves); structural annotate run 0.11 s |
| exhaustive constrained, 19 labels, R = 1000, k-mer / trace | mini_refseq | 0.49 s / 0.43 s |
| the records model, 42 records, 9.28 Mbp | — | 2.8 s to build (packed k-mers, hopscotch; the same with `std::hash` on packed k-mers took 80 s — identity hashing of 2-bit k-mers in an open-addressing table) |
| 2000 bp exact flank recovery, single label | UHGG public index (9.68 G nodes, 4,644 labels, RowDiff<BRWT>, no coords) | ~1.4 annotation rows per step; 6 s wall dominated by index load; 6.6 GB RSS |
| 16S window, 8 carriers, multi-label == union of 8 single-label walks | UHGG | 10 leaves; every leaf exact in its MGnify genome and absent from the others |
| switching at a real contig end | UHGG | forbid → 0 bp `label_lost`; constant 1 / budget 1 → 500 bp via two labels; `route_bp = 0` label verified exact over 800 bp, `route_bp = 451` label's agreement ends at exactly 419 |
| **Cost at random positions**, 12 random 50-mers (from gut genomes, fully present in the graph) | SRA, as below | `/search` 50 bp: **2 ms** (1–3). Depth 100: default walk **15 ms** (5–63), exhaustive constrained trie **17 ms** (3–595), label-free trie extraction **20 ms** (3–1708); ≈ 0.05–0.1 ms per node entered; 23/24 arms complete at 100 bp. The tail is the repeat: 616 labels on the seed → 114 constrained / 1139 structural leaves, the 1000-live-path cap at 51 bp |
| **Label-consistent walks to 10 kb**, 24 random positions, derived labels (median 6) | SRA, as below | strict (`forbid`, limit 0): flank median **224 bp**, 10/48 arms ≥ 1 kb, max 7.5 kb, none reaches 10 kb; ends `branch` 205 / `label_lost` 94 / `dead_end` 9; 44 ms median (max 884). Limit 1: median 428 bp; limit 2: 520 bp (13 → 14 arms ≥ 1 kb). `min_successor_labels 2`: 74 bp (too few labels for a quorum). One switch + limit 0: **57 bp** (see §5.2.7). Exhaustive capped at 200 live paths: cap hit on 27/48 arms, 24 `edge_reuse` blocks in total, no `edge_reuse_rc`, no `rejoined_seed` — loops are rare at random positions and, under a branch limit, are met as `branch` at the junction before they are entered |
| **PRIMARY regime on a real index**: a full-length 16S (1396 bp window of PZ326290.1) | SRA `sra_random_100studies` (42.4 G nodes, k = 31, primary `DBGSuccinct` small + `RowDiff<BRWT>` 19.6 GB, no coords; 12 GB RSS, 20 s load, server mode) | `resolve`: 1.0 s / 0.4 s per 16S query, ≥ 300 labels (discover cap hit), 6 carriers of the whole window (one study). Exhaustive constrained trie, derived labels, R = 300: **1000 live paths at 237 bp (right) / 291 bp (left)** — 1116 / 1437 leaves, 44.7 k steps, 124.6 k rows requested, **5.4 s**; the structural trie (2836 labels recorded, up to 827 at one node) trips the same cap at 109 / 46 bp, 5.0 s. **The JSON-level oracle at the common complete depth agrees on both arms: 114 = 114 claims (right, 10 of them labels ending inside walks), 2 = 2 (left).** One switch allowed: 398 / 282 switch events, 19 `loss_budget` ends. Default strategy (branch limit 0, merge on), R = 3000: one leaf per arm, **61 bp right / 52 bp left, every label ends `branch`** — 0.2 s; each of the 6 samples alone gives the same 61 / 52 (one 18): the 16S flank is ambiguous *within* each metagenome within ~60 bp |

The complexity audit (asymptotics in `walker.cpp`, `label_oracle.cpp`, `resolve.cpp`) found six superlinear
terms. Fixed and tested (C20): `hits_from_tuples` (O(coords × hits) → O(coords) per node — on the coordinate path
`refseq33m` uses) and `max_support` (re-materialising per candidate → once per round). Landed in `01a5230a`
("Remove three superlinear terms from the walker, byte-identical"): the edge-reuse check keyed by
`(edge, segment)` and probed over the current path's own ancestors when that is fewer than the edge's uses (was
O(P²) per level in `keep` mode, which the preset forces; the new `EdgeReuseProbes` test fails on the old code,
62,074 probes against a bound of 32,220), the reconvergence state union as a two-pointer merge (was O(σa·σb)),
and an explicit bound and counters on re-minimisation rounds at an ambiguous node (`edge_reuse_probes`,
`reminimisation_rounds`, `max_reminimisation_rounds` under `counters`). Each was required to leave every result
byte-identical; the three-way suite is the regression net.

1. **The term that decides `refseq33m`.** Every label there is a header label, so every step reads a *tuple* row
   from `TupleRowDiff<TupleCSC<BRWT>>`: the full row and all coordinate lists, of which ~20 labels are wanted.
   Row population W, not the walker, is the cost. Round 1 showed the proposed boolean `rd_direct` (XOR of
   `get(row, col)` along the row-diff path) is **invalid for tuple matrices** (membership is "tuple non-empty",
   not a bit XOR). What *is* the right selected-column extraction for tuple row-diff — decode only the permitted
   columns' tuples along the anchor path? Is it implementable without touching the on-disk format, and what does
   it cost per step?
2. **Annotate mode on a wide locus.** The oracle reads full rows by design and the spec calls it "a verification
   tool at small radius". On `refseq33m` a node can carry 10⁴ labels; is there any use of annotate mode that
   remains sensible there, or should the API refuse it above a label count?
3. **The perf fixes.** For each, is the new code exact (same results, same reasons, same event order)? The brief
   is in the commit messages; the edge-reuse change in particular must preserve C5's precedence and the `blocked`
   events.
4. **Caches.** `LabelQuery`'s cache is bounded by entries, not bytes (each entry carries per-label coordinate
   lists); fetch keeps the working set resident. At |P| = 10³ with coordinates, what is the realistic footprint?
5. Where would you expect the first wall on `refseq33m` now — and what single measurement would you run first?

### 5.4 Evidence honesty

1. `LabelEnd::route_bp` (0 = the label carries this path's spelled bases; > 0 = merged in at that depth, supports
   the leaf and the shared tail only) and `LabelArmSummary::direct_bp` (a label-consistent *route* of that length
   exists, not a contiguous occurrence) were corrected in wording after round 2. Is every place a label is named
   now unambiguous about which of the three kinds of support it has — k-mer route, trace, or merged-in?
2. In annotate mode `label_summary` is computed from the recorded sets (`direct_bp` = continuous presence along
   *some* route, exact when no list was cut, a lower bound otherwise). Is "lower bound otherwise" actually a lower
   bound, or can a cut list make it an *over*-estimate?
3. The derived set is the seed's *carrier* set, not one hit. Under `forbid` this is `path_common` over the carriers.
   Can a consumer mistake a derived-set result for evidence about one particular record?
4. `Event::text` qualifiers (`superseded`, `minority`, `below_min_labels`, `split_limit`, `hairpin`,
   `switch_sources`) still carry distinctions that are not `EndReason` values. Round 1 flagged this as lossy; it
   was kept. Is there a concrete misreading a consumer would make?

### 5.5 Design critique

1. **The walk rule.** "No (k+1)-mer twice per path" makes the set of walks finite but has consequences the tests
   now pin as *known*: a homopolymer A^(k+3) yields the walks A^(k−1)·Z, A^k·Z, A^(k+1)·Z and **never the
   record's own** A^(k+3)·Z — under k-mer *and* trace support (the coordinates force the record's path onto the
   loop edge a second time). Likewise a tandem repeat is unrolled at most once. Is this the right rule for a tool
   whose output is "candidate sequence"? Alternatives: no *node* twice (stricter, fewer walks), a per-path edge
   multiplicity bound, or letting `trace` override edge reuse when coordinates continue (then trace walks could
   be longer than k-mer walks, and `check_trace`'s prefix property would have to go).
2. **Annotate as oracle vs annotate as product.** It was built as a verification tool; the user also wants it as
   a way to "extract the whole trie from a match up to a size, independent of labels but with labels annotated".
   Are the size bounds (`max_steps`, `max_output_bp`, `max_paths`, `max_live_paths`, time) plus `complete_to_bp`
   the right contract for that, or is a bound on the *number of nodes entered per level* missing?
3. **Breadth-first as the only order under the preset.** The completeness boundary depends on level-synchronous
   exploration; `lowest_loss_first` and `most_supported_first` are refused with `exhaustive`. Is that the right
   coupling?
4. **Merging is still on by default** outside the preset, and still the thing that weakens per-path evidence
   (`route_bp`). Round 1 asked whether it should be opt-in; I have not changed it. Argue it either way.
5. **Not built, by choice:** an `hll` cost model (the user's formulation needs 0 or 1 switches, not similarity);
   discovered switch targets (hopping into a label that does not carry the seed); a per-label `label_walks` view
   (the union-of-single-label check covers the semantics; the view would be convenience). Which of these would
   you build first, if any?

### 5.6 API and usability

1. Given a result, can an agent tell why a trie stayed shallow or exploded, and what to change? Evidence offered:
   `complete_to_bp`, `cap_trigger`, `labels_per_node`, `growth` bins, `branch_events`, `needed_budgets`,
   `walk_rule`. What is missing for the loop "probe → read → retune → extend"?
2. The derived-label default: is returning the *carrier set* the right default when a user says "traverse from
   how the match is labelled", or should the API require the user to pick one hit's labels?
3. Trace support is BASIC-only (coordinates carry no strand in canonical indexes) and rejects merging. Is that a
   reasonable surface, or a trap?

---

## 6. Already found and fixed — do not re-report, but do check the fixes

From round 1 and its follow-up:

1. **Blocker:** trace support merged coordinate sets across a reconvergence and could report direct support no
   single occurrence justifies → `trace` + `merge` is rejected with an actionable error.
2. `LabelQuery` eviction mid-call → "Couldn't find key" crash → fetch keeps the call's working set; always inserts.
3. `max_support` did not clip runs between rounds → fixed, tested (`MaxSupportClipsRunsBetweenRounds`).
4. Multi-graph requests rejected with `unknown field 'graph'` → routing fields declared; regression test.
5. `resolve` accepted and ignored `bounds.time_budget_ms` → rejected.
6. Per-request annotation counters broke C9 → per-seed deltas.
7. Negative `kmer_interval`, `annotation.access` other than `auto` → validated / rejected.
8. C1 qualified by `max_switch_sources`; C3 restricted to `on_reconverge: keep`.
9. Doc wording: "contiguously in its own sequence" overstated k-mer evidence; `route_bp` comments were inverted
   (the unsupported part is the prefix *before* `route_bp`); the spec's exact read-out requires
   `route_bp: 0` as well as `entered_by: seed` / `from_bp: 0`.
10. Derived labels: O(L²·M) seed validation → one data-major sweep; no clock check during derivation → time budget
    honoured, plus `--traverse-max-seed-labels`; `table` cost refused with derived labels (a by-name table would
    silently mis-apply).
11. From the review of the trie work: the oracle compared leaves, not claims → claims; event label lists and the
    root carried no totals → every list has its true count; `exhaustive` could be subverted by
    `max_switch_sources` under a `table` cost, by a derived set cut by `max_seed_labels`, and by the 64-label
    annotate default → refused / reported; the per-node cap was enforced after the coordinate work it should
    bound → before; `walk_rule` understated the enforced rule (even-k hairpins, hashed edges) → stated;
    live-label counts over cut lists were stamped exact → inexact.

---

## 7. Known limitations — please do not spend time re-reporting these

- No `hll` cost model; no discovered switch targets; no tip/bubble windows; no GFA export; no `label_walks` view.
- `trace` is BASIC-only and rejects merging; with *column* labels it cannot detect a record boundary that has
  consecutive global coordinates (header labels can).
- The homopolymer / tandem-repeat consequence of the walk rule (§5.5.1) is pinned, not solved.
- `discover` in `resolve` does not scale (it truncates output, not work) — the supported path is labels from
  `/search`; `resolve` is for a known label list.
- Caches are bounded by entry count, not bytes. The header → (column, seq_id) index is keyed on a raw
  `CoordToHeader*`.
- No measurement on a `refseq33m` shard. The PRIMARY regime is measured on the SRA index (§5.3) but not
  validated against records there (SRA reads have no "record" to compare with); its real-index evidence is the
  JSON-level oracle agreement and the unit suite's canonical/primary cells.
- `switch_on: any` has no independent reference; `on_reconverge: merge` is checked only through route soundness.

---

## 8. How to verify anything you suspect

```bash
# Release and Debug (asserts on). AppleClang needs the -Wno-error flag; GCC does not.
cd metagraph/build       && cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS=-Wno-error=unused-result .. \
                         && make -j metagraph unit_tests
cd metagraph/build_debug && cmake -DCMAKE_BUILD_TYPE=Debug   -DCMAKE_CXX_FLAGS=-Wno-error=unused-result .. \
                         && make -j unit_tests
# re-run cmake after ADDING a file: the source globs are not CONFIGURE_DEPENDS

./unit_tests --gtest_filter='*Walker*:*LabelOracle*:*Resolve*:*Trie*:MiniRefSeq*'   # 195 tests, ~12 s / ~15 s
./unit_tests --gtest_filter='*TrieCases*'                                            # the 47 edge/non-edge cases
TRIE_TIMINGS=1 ./unit_tests --gtest_filter='MiniRefSeq.TrieContract*'                # the real-index contract, phase timings

# the real-format fixture (~25 s, downloads 42 NCBI records; MiniRefSeq* skips without it; run tests with cwd = build/)
../scripts/traversal/build_mini_refseq.sh ./mini_refseq ./metagraph ../scripts/traversal/mini_refseq_accessions.tsv

# CLI and HTTP routes end to end
./test_venv/bin/python ../integration_tests/main.py --test_filter='*Traverse*'      # 23 tests, ~13 s
```

Status at `01a5230a`: **195** traversal unit tests pass in Release and Debug, **23** integration tests pass. The
oracle contract held on every cell (4 graph × annotation pairs × 3 modes × 5 permitted sets × 2 arms on the
original fixture; 2 pairs × 3 modes × 2 arms on each of the 25 cases; 16 random fixtures at k = 7; the real
index). Every disagreement found while building the suite was on the test side (my literal expectations), never
in the walker — **which is exactly the pattern a circular oracle would also produce. Please distrust it.**

---

## 9. How to report

For each finding: **severity** (blocker / major / minor), `file:line`, a **concrete failure scenario** (inputs or
graph shape → the wrong output), and a suggested fix. Separate **confirmed** (traced), **suspected** (say which
experiment would settle it) and **design disagreement**. If the honest summary is "the verification is sound and
the open problem is tuple-row extraction at scale", say that; if you find the oracle circular anywhere, that is
the finding I most want.
