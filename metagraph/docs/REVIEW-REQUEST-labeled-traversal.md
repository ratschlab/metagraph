# Review request: label-consistent graph traversal in MetaGraph

**What I want from you:** an adversarial technical review of a new, uncommitted feature. Please try to *break*
it rather than confirm it. I care most about defects that would make the output **claim more biological support
than the evidence warrants**, and about whether the design will survive the scale it is aimed at.

You have access to the repository. Nothing is committed, so see §2 for how to find the code.

**Priorities, in order:**

1. **Correctness** — real defects in the traversal semantics.
2. **Evidence honesty** — can any output overstate biological support? (One such bug was already found; §5.2.)
3. **Efficiency at scale** — will this work on a 33 M-label, 627 G-node index, or only on toys?
4. **Robustness** — large disk-backed annotations, unmasked graphs, concurrent server use, hostile input.
5. **Design critique** — is the formulation right at all? See §6 for the alternatives I rejected.
6. **API / usability** — can an LLM agent actually drive and tune this from the JSON it gets back?

Please be blunt. If the whole approach is wrong, say so and say what you would do instead.

---

## 1. What the feature does, in one page

MetaGraph indexes DNA as a de Bruijn graph (nodes = k-mers) plus an **annotation matrix** mapping each k-mer to
the **labels** (samples, genomes, accessions…) that contain it. Existing search answers *"which labels contain my
query?"*. It cannot answer *"what sequence lies **around** my query in those labels?"* — the flanks are not in the
query, so they cannot be returned by a lookup.

This feature adds that, as a three-step pipeline:

```
query sequence
  └─ resolve  → which labels support which contiguous k-mer intervals of the query   (no traversal)
       └─ select  → freeze an explicit seed list (sequence + label set), policy recorded
            └─ traverse → extend each seed outward through the graph, constrained by labels
```

The traversal walks outward from the seed, and at every step keeps only the labels that still support the path.
A step is admissible while at least one permitted label survives. The central generalisation: traversal state
carries a **set** of labels, not one, so one seed carried by 1,000 labels is **one** walk whose label set splits
where the genomes diverge — rather than 1,000 walks.

On top of that set are two knobs that subsume the modes in the design note:

- a **label-change cost** `cost(ℓ → ℓ′) ≥ 0` plus a per-arm **loss budget**: with cost `forbid` a single label must
  carry the whole path (strict, the default); with a finite cost the path may switch labels, paying for it, which
  is the "close enough labels" idea;
- **per-lineage branch limits** and quorum thresholds, so breadth is a property of the *strategy*, not just the
  graph: a tight strategy prunes most structural branches and runs deep with a small result, a loose one fans out
  and hits a cap.

Design note (goals, written by the project owner, deliberately cautious about what graph evidence proves):
`docs/DESIGN-labeled-traversal-endpoint.md`. Implementation spec (mine, the contract this code is written
against): `docs/SPEC-labeled-traversal-core.md` (785 lines). **Where code and spec disagree, that is a finding.**

---

## 2. How to see the code

Branch `gr/labeled-traversal`, head is still upstream `1de02cec` — **there are no commits**. So:

```bash
git status --short                 # 5 modified files, the rest untracked
git diff                           # ONLY shows the 5 modified files
```

Everything else is **untracked** and will not appear in a diff. The complete feature:

| Path | Lines | Role |
|---|---:|---|
| `src/graph/traversal/traversal_types.{hpp,cpp}` | 95 / 74 | shared enums: regimes, label kinds, support kinds, end reasons |
| `src/graph/traversal/label_oracle.{hpp,cpp}` | 204 / 525 | per-request context: graph regime, oriented-node → annotation-row mapping, batched label membership |
| `src/graph/traversal/resolve.{hpp,cpp}` | 152 / 548 | `resolve` (support profile) and `select` (seed freezing) |
| `src/graph/traversal/walker.{hpp,cpp}` | 330 / 2089 | the traversal itself — **the core, review this hardest** |
| `src/cli/traverse.{hpp,cpp}` | 90 / 944 | strict JSON parsing/validation, serialisation, CLI entry point |
| `tests/graph/traversal/test_walker.cpp` | 2085 | traversal semantics, typed over graph × annotation × mode |
| `tests/graph/traversal/test_label_oracle.cpp` | 345 | regime/key mapping, access-path agreement |
| `tests/graph/traversal/test_resolve.cpp` | 374 | support profiles, selection policies, scale |
| `tests/graph/traversal/test_mini_refseq.cpp` | 376 | **against a real index** in the public `refseq33m` format |
| `integration_tests/test_traverse.py` | 319 | CLI + HTTP routes end to end |
| `scripts/traversal/build_mini_refseq.sh` | — | builds the real-format fixture from 42 NCBI accessions (~25 s) |

Modified existing files (151 inserted lines total): `src/cli/server.cpp` (+101: three routes and a shard-
addressing helper), `src/cli/config/config.{hpp,cpp}` (+43: the `traverse` command and its flags),
`src/main.cpp` (+4: dispatch), `integration_tests/main.py` (server-based module must not be chunked).

**Suggested reading order:** spec §1–§9 → `traversal_types.hpp` (the enum contract: note `EndReason` has no
`switched`/`superseded`/`minority`/`hairpin` values — those are `Event::text` qualifiers) → `walker.hpp` (the
output contract; every field is an evidence claim) → `walker.cpp`, whose **lines 42–78 are an
implementation-notes block listing the known deviations** → `label_oracle.cpp` → `resolve.cpp` →
`test_walker.cpp`. The `MiniRefSeq*` tests must run with `cwd = build/` (they locate the fixture by relative
path and skip if it is absent).

---

## 3. Background you need (MetaGraph specifics)

- **Label kinds.** A *column label* is an annotation column (on `refseq33m`, an NCBI taxid). A *header label* is
  one indexed sequence inside a column — an accession — recovered through a `CoordToHeader` (`.seqs`) sidecar that
  maps a column's k-mer coordinates to the FASTA record they came from. **On `refseq33m` the useful labels are
  header labels**, so the per-accession path is the production path, not an edge case.
- **Three graph regimes**, and oriented traversal node ids are *not* annotation row keys:
  - `BASIC`: key = node.
  - `PRIMARY`: the graph is wrapped in `CanonicalDBG`; ids above an offset are reverse-complement orientations;
    key = `get_base_node(node)`.
  - `CANONICAL` (native, both strands stored): no arithmetic mapping; the key is `map_to_nodes()` of the spelled
    k-mer (for `DBGSSHash`, `min(n, reverse_complement(n))`).
  Passing an oriented id straight to the annotation indexes out of range on `PRIMARY` and **silently returns an
  empty row** on `CANONICAL` — a false negative with no error. `CanonicalDBG::get_mode()` returns `CANONICAL` for
  primary graphs, so the regime is detected by `dynamic_cast`, not by mode.
- **Annotation representations** differ in what a membership test costs: `ColumnMajor`, `RowFlat`, `BRWT` support a
  single-cell `get(row, col)`; `RowDiff` variants do **not** — a row is reconstructed by walking a successor path
  to an anchor (soft cap ~100 steps) and XOR-ing diffs. Coordinates live in `TupleRowDiff<TupleCSC<BRWT>>`, read
  with `get_row_tuples`. The design note is explicit that whole-column extraction must never be an online path.
- **Dummy k-mers.** `DBGSuccinct` without a mask exposes sentinel `'$'` neighbours; `refseq33m` ships no
  `.edgemask`, so this is the production case. Degree functions count them, so branch decisions must come from
  the character-filtered, label-filtered neighbour list.
- **Conventions in the output.** Support runs are **0-based half-open intervals of k-mer starts**; a run `[a,b)`
  covers bases `[a, b+k−1)`. Extension is measured in **bases outward from the seed boundary** per arm, counted
  outside the seed; step *j* adds outward base *j−1*. Left-arm sequences are reported in natural orientation.
- **What the evidence means.** `within_label` support proves every k-mer of the path carries that label. It does
  **not** prove the path is a contiguous molecule: one label can hold many contigs, repeat copies or strains, so a
  graph path through them is a *candidate* sequence. The design note is emphatic about this, and it is the crux of
  §5.2 below.

**The target index** (public, verified via its own `/stats` and S3 listing): `refseq33m` — k=31, **basic** mode
`DBGSuccinct`, **forward strand only**, 626,753,663,468 nodes (362 GB), annotation
`RowDiff<BRWT>` **with k-mer coordinates** (333 GB), 32,881,371 labels, columns = taxids, accessions via a
0.64 GB `.seqs` sidecar, no dummy mask.

---

## 4. The claims to attack

These are the load-bearing claims. Each is meant to be decidable by reading the code. **If any is false, that is
the finding I most want.**

| # | Claim | Where |
|---|---|---|
| C1 | For a fixed path with no branch limit, the per-label loss is the exact minimum over label assignments. Branch-limit pruning is explicitly greedy, not a constrained optimum. | spec §6.3; `walker.cpp` `derive()` |
| C2 | Branching is counted per **lineage**, not per label name: with a finite switch cost a label can be renamed, and a renamed lineage must not evade its branch limit. | spec §6.4; the ambiguity/exclusion loop in `process_item` |
| C3 | With `forbid` and `max_label_branches = 0`, the result equals the union of independent single-label walks, and each arm has at most \|S\| leaves. | `test_walker.cpp` `ReferenceWalker` |
| C4 | Every path uses each directed edge at most once, **per path**, with merged lineages' edge sets unioned (conservative). Canonical regimes treat `u→v` and `rc(v)→rc(u)` as one edge. | spec §6.6; `edge_key()`, `is_ancestor_or_self()` |
| C5 | Precedence when a step is blocked: `rejoined_seed` > `edge_reuse_rc` > `edge_reuse`; seed re-entry is keyed on the target **node**, not the edge. | `block_rank()` |
| C6 | A hairpin step (a self-reverse-complementary (k+1)-mer) is skipped by default and never counts toward ambiguity — otherwise every label would stop at any RC-palindromic (k−1)-mer. | spec §6.5 |
| C7 | `trace` support requires strictly consecutive coordinates in one indexed sequence: `+1` on the right arm, `−1` on the left. | `derive()` coordinate filter |
| C8 | Results do not depend on `annotation.batch_kmers` — prefetching only warms a cache and can never change which steps are taken. | spec §6.8; `prefetch()`; `Determinism` test |
| C9 | The `result` object is byte-identical across runs, seed permutations, and one-seed-per-request splits. Wall-clock/cache counters live in a separate `timing` object excluded from that contract. | spec §6.8 |
| C10 | `complete` means the requested domain was fully explored; `truncated`/`pruned` mean it was not. A capped result is never evidence of absence. Seed-level caps stop both arms, per-arm caps stop one. | spec §6.7 |
| C11 | Strict JSON validation: unknown keys, bad types, out-of-range values and unreachable `extra` labels are rejected naming the field; nothing is silently weakened. | `traverse.cpp` `Strict` |
| C12 | A traversal addresses exactly one physical (graph, annotation) shard; an index name covering several is rejected unless disambiguated. | `server.cpp` `resolve_traverse_index()` |

---

## 5. Specific questions

### 5.1 Correctness

1. **C2, the re-minimisation loop.** When an ambiguous lineage is excluded, every successor's state is recomputed
   and the loop repeats. Can this loop fail to terminate, or terminate with a state that depends on iteration
   order? Is the bound (`round <= |state|`) actually sufficient?
2. **C1 vs C2 interaction.** `derive()` fixes one predecessor per target *before* branch-limit exclusion runs. Can
   that discard a cheaper assignment that would have satisfied the branch limit — and if so, is the greedy
   behaviour the spec admits, or something worse?
3. **C4 after a merge.** Segments gain multiple parents at a reconvergence and `is_ancestor_or_self` walks all
   branches. Is the ancestor marking still valid when segments are created *during* the same level as a merge?
   Can an edge be reused because a lineage's ancestry was not yet linked?
4. **C7 on the left arm.** Coordinates decrease leftward. Check the `coord == 0` boundary and the case where a
   label has several coordinates at one node (a repeat) — can a trace jump between copies and still be reported
   as one run?
5. **Regime mapping.** For native `CANONICAL` graphs the key needs the spelled k-mer. Is the spelling always the
   oriented k-mer at that node, including on the left arm and after a hairpin? A wrong spelling yields a silently
   empty row, i.e. a false "label lost".
6. **`resolve`'s trace runs** (`resolve.cpp`): the run/break construction is subtle. Does a coordinate jump always
   produce exactly one break and two runs, and is `kmers_supported` (k-mer presence) correctly distinguished from
   trace continuity?
7. **`select`'s `max_support` sweep**: it claims that for a fixed left end the best interval ends at the smallest
   run end ≥ `a + min_kmers`, with candidate left ends being run starts and `end − min_kmers`. Is that exact, or
   is there a case where a non-endpoint left end wins?

### 5.2 Evidence honesty — the highest-value area

I found exactly one real bug this way, and I want you to look for more of the same shape.

**What happened.** I extended a blaNDM-1 gene 3 kb through the real fixture and then checked every recovered
flank against the source RefSeq record it was attributed to. 29 of 30 were exact substrings. The one exception
was a path that passed a **reconvergence** at 1239 bp: two routes met, the merged node's label set became the
union, and an accession that arrived via the *other* route was listed at a leaf whose spelled sequence it does not
contain. The output was asserting that an accession carries a 3 kb flank it does not carry.

**The fix** was to expose `LabelEnd::route_bp` (already computed internally for continuations): `0` means the label
travelled this path's own bases; non-zero means its route diverged there, so it supports the leaf *node* but not
the spelled bases. `test_mini_refseq.cpp` now asserts the substring property for all `route_bp == 0` labels and
requires the merged case to occur.

**Questions:**

1. Are there other paths by which a label can be reported against sequence it does not support? Candidates to
   check: `LabelRun` intervals after a merge; `LabelArmSummary::direct_bp` vs `reach_bp`; `Segment::labels_start`
   / `labels_end` on a merged segment; `Continuation::labels`.
2. Is `direct_bp` (advertised as "the longest run entered from the seed", the only per-sample evidence) actually
   immune to the merge problem?
3. After a switch, a path is supported by a *chain* of labels, not one sample. Is that distinguishable in the
   output in every place a label is named, or can a switched-in label be mistaken for a direct carrier?
4. Does any end reason overstate knowledge? Specifically `LABEL_LOST` ("no continuation carries it") vs merely
   "we did not look" — e.g. when a successor was skipped as a hairpin or blocked by edge reuse.
5. Several semantic distinctions are carried in `Event::text` qualifiers (`"superseded"`, `"minority"`,
   `"below_min_labels"`, `"split_limit"`, `"hairpin"`) rather than distinct `EndReason` values. Is that a
   lossy encoding a consumer can misread?

### 5.2b Open leads I have not resolved

An internal pass over this code produced the following specific leads. I verified and fixed three of them
(§5.2c); **the rest are open and are good starting points** — each names a place and a decidable question.

1. **Determinism of the `annotation` counters.** They are per-request (per `LabelOracle`), not per seed
   (`walker.cpp` ~:2049, serialised `traverse.cpp` ~:632). For a request with seeds [A, B], does seed B report
   the same counters as when submitted alone? If not, claim **C9** is false as stated, or the `annotation` block
   belongs in `timing` (which is excluded from the contract).
2. **`max_switch_sources` truncates by loss, not by switch cost** (`derive`, `walker.cpp` ~:1047). With a `table`
   cost, |σ| = 70 and the cap at 64, if the only finite-cost source into the surviving target ranks 70th, is it
   found? If not, **C1's exactness fails for table costs** — and does `ReferenceWalker`/the loss oracle ever run
   with |σ| above the cap?
3. **Caps are evaluated before `merge_level`** (`walker.cpp` ~:1708 vs ~:1523). Can an arm be truncated for
   exceeding `max_live_paths` by heads that were about to merge away? And is the growth-bin invariant
   (`max live paths ≤ max_live_paths`) trivially true because it is recorded *post*-merge while the cap is
   checked *pre*-merge?
4. **`select` MAX_SUPPORT in later rounds** (`resolve.cpp` ~:344–424). Candidate right ends are restricted to run
   ends and any interval overlapping an already-taken one is dropped wholesale. With one label supporting
   `[0,1000)` and `max_seeds: 3`, how many seeds come back, and does the "one binary search per candidate left
   end" optimality argument still hold once `taken` is non-empty?
5. **`merge_overlapping` (default true) can shrink a long seed** (`resolve.cpp` ~:458). With candidates
   `A=[0,800)` and `B=[700,760)` under `longest_first`, does the frozen seed collapse to `[700,760)`? That would
   directly defeat the "seeds as long as possible" requirement.
6. **Trace semantics** (`resolve.cpp` ~:226, `walker.cpp` ~:1560). Are coordinate sets *unioned* across
   occurrences, so a repeat lets a trace hop between copies inside one record? Is `trace_break` really a record
   boundary, given that `num_kmers_in_sequence()` is never consulted — i.e. is a true record end ever detected?
7. **Hairpin vs seed re-entry precedence** (`walker.cpp` ~:1586, `is_hairpin` ~:763). For a circular molecule in a
   PRIMARY index whose closing step is self-reverse-complementary with `hairpins: skip`, is the reason
   `rejoined_seed` (with `overlap_bp`) or a plain dead end? And for even k, does a palindromic source k-mer make
   *every* outgoing step a hairpin, killing all labels at once?
8. **`resolve` is quadratic in the label count** (`resolve.cpp` ~:207–263 linear-scans each node's hit list per
   label). With `discover: {max_labels: 1000}` on a 100 kb query where most k-mers carry hundreds of labels, what
   is the wall time and peak RSS — and is there anything that aborts it? (`bounds.time_budget_ms` is parsed and
   discarded in `resolve`.)
9. **Internal errors return HTTP 400** (`server_utils.cpp` `process_request`, unchanged). Spec §10.3 requires 500
   for non-client errors. What status does a `CanonicalDBG` primary-graph inconsistency produce today, and does
   the CLI still write its per-request error JSON on a non-`InvalidRequest` throw?
10. **The header→(column, seq_id) index is cached forever, keyed on a raw `CoordToHeader*`**
    (`label_oracle.cpp` ~:157–204). Is that identity-safe if one sidecar is destroyed and another allocated at the
    same address, and how large does this static grow for a `refseq33m`-scale `.seqs`?

### 5.2c Already found and fixed — do not re-report, but do check the fix

1. **Multi-graph mode was dead on arrival.** The server consumes `graph` / `graph_path` to pick a shard, but the
   strict parser never declared them, so the unknown-field check rejected **every** multi-graph request with
   `unknown field 'graph'`. Invisible to the suite because the integration tests run a single-graph server. Fixed
   by declaring both as known routing fields; a regression test now asserts the error is a routing error and not
   a parse error. **Claim C12 was untested when written — please check it properly.**
2. **`direct_bp` and the merge.** I first clamped `direct_bp` at the merge depth and a test correctly rejected it:
   a reconvergence joins paths *at the same node*, so bases after the merge are identical on every route in, and
   the label really does carry them in its own sequence. The resolution is that `direct_bp` is per-label along its
   **own** route, while `LabelEnd::route_bp` is the per-path qualifier. Please check that this distinction is
   actually sound and consistently applied (`LabelRun::route_bp` is now also reported).
3. **`annotation.access` was parsed and silently ignored**; now rejected unless `"auto"`, since the spec promises
   nothing is silently weakened.

### 5.3 Efficiency at `refseq33m` scale — I am most unsure here

Measured on the 9 Mbp fixture: a 3 kb two-arm traversal from a seed with 19 labels took **0.14 s** and requested
**24,877 annotation rows**. That is roughly one row per step per arm, which is the expected shape — but the
fixture has **9 columns**, and `refseq33m` has **33 M labels and a 333 GB coordinate annotation**.

1. **The tuples path is the production path.** Header (accession) labels require coordinates, so every
   `refseq33m` traversal uses `get_row_tuples` on `TupleRowDiff<TupleCSC<BRWT>>` — full row *and tuple*
   reconstruction per node, where a node's row may carry thousands of labels each with coordinate lists, of which
   we want ~20. Is that viable at 10³–10⁴ steps per request, or fatal? What is the realistic per-step cost?
2. **The selected-column accessor was specified but not implemented.** The spec (§8.2) proposes `rd_direct`:
   for `RowDiff<M>` where `M` supports `get(row,col)`, compute membership as the XOR of `diffs().get(r, c)` along
   the row-diff path for only the permitted columns, decoding incrementally along an unbranched run. Only
   `direct` / `rows` / `tuples` exist today. Is `rd_direct` the right fix, is it correct for *tuple* matrices
   (where membership is "tuple non-empty"), and is there a better approach I am missing?
3. **Cache growth.** `LabelQuery`'s cache is bounded by entry count (default 10⁶), not bytes, and each entry holds
   a per-label coordinate list. At large permitted sets, is this a memory blowup? Is full eviction on overflow the
   right policy?
4. **Per-step graph cost.** On `PRIMARY`, each first visit adds O(k) reverse-complement index work; the left arm
   on `DBGSuccinct` uses a `NodeFirstCache`. Are the spelling hints actually passed everywhere they should be, or
   does some path fall back to `get_node_sequence()` per step?
5. Where would you expect the first wall to be hit on `refseq33m` — annotation I/O, row reconstruction, label-set
   churn, or output size? What measurement would you run first?

### 5.4 Robustness

1. Concurrency: all annotation reads are const and the loaded index is shared, with per-request caches. Is
   anything actually shared-mutable? (Note a per-request `CanonicalDBG` clone is made precisely to avoid mutating
   the shared wrapper.)
2. The HTTP handlers run on the server's io threads, whose default count is **1**, so one long traversal blocks
   `/search`. Documented, not fixed. How would you bound it?
3. `process_request` maps **every** `std::exception` to HTTP 400, so an internal failure looks like a client
   error. Pre-existing, but my `InvalidRequest` relies on it. Worth changing?
4. Hostile or pathological input: a seed of exactly `k`; a seed that is one long homopolymer; a 1 Mbp query to
   `resolve`; `max_extension_bp` of 10⁹; 10⁴ seeds in one request; a label name containing a tab or comma (seed
   ids are built from length-prefixed label strings — is that injection-proof?).
5. Disk-backed annotations (`RowDisk`) open a file view per `get_rows` call. Does the batching avoid per-node
   calls on those representations?

### 5.5 Design critique

Please challenge the formulation itself, not just the code.

1. Is **label set + change cost + loss budget** the right generalisation of the design note's discrete modes
   (`within_label`, `path_common`, `any_of`, `sequence_trace`)? It makes them special cases — but does it make the
   evidence harder to interpret than four explicit modes would be?
2. **Reconvergence merging** is on by default: it is what keeps a population-scale walk from producing one path
   per label through every bubble, but it is exactly what caused the §5.2 bug. Is merging-by-default right, given
   that it weakens per-path evidence? Should it be opt-in?
3. `max_support` selection maximises the *number* of labels covering a block, so on blaNDM it prefers the 231-k-mer
   prefix shared by 25 carriers over the 783-k-mer gene carried in full by 19. The project owner wants seeds "as
   long as possible", so `longest_first` exists too. Is a single objective with a length floor better than two
   policies?
4. Deliberately **not** implemented: tip windows and bubble windows (skip a short dead-end or a closed bubble
   without spending a branch). Without them, a SNP bubble inside one sample stops that label at the first bubble
   under `forbid`, and raising the branch limit grows the result exponentially. Is that the right thing to build
   next, or is there a cheaper mechanism?
5. The intended consumer is an LLM agent iterating: probe cheaply, read the diagnostics, retune, extend. Does the
   contract support that loop, or does it presume a human reading JSON?

### 5.6 API and usability

1. Given a result, can a consumer actually tell **why** a traversal stayed shallow or blew up, and which knob to
   turn? The evidence offered is `growth` (per-bin live paths, distinct live labels, splits, reconvergences,
   blocked repeats, label ends by reason), `branch_events`, `needed_budgets` (what budget *would* have kept a
   label), and `cap_trigger`. What is missing?
2. Output size: bounded in bases (`max_output_bp`) and events, but label lists scale with the permitted set. At
   \|S\| = 1,000 is the JSON still usable? `output.detail: summary|tree|full` exists — is `summary` the right cheap
   probe?
3. Is `continuation` (a leaf's tail sequence + the labels covering it, re-submittable as a seed) a sound basis for
   iterative deepening, given that it resets edge-reuse and branch state?
4. Are the knob names and defaults in spec §5 the ones you would expose?

---

## 6. Known limitations — please do not spend time re-reporting these

- **No `hll` cost model.** The "close enough labels" cost is `forbid` / `constant` / `table` only. The planned HLL
  variant (from branch `origin/hm/aln_alt_label_change`: `log2(min(|B|, |A|+|B|−|A∪B|)) − log2|B|` over per-column
  HyperLogLog sketches) needs a `libcount` submodule and a `.hll` sidecar, and must clamp the estimate so sketch
  error cannot make a cost negative. **Do tell me if the cost model's shape is wrong.**
- **No tip/bubble windows** (fields exist, stay zero). See §5.5.4.
- **No GFA export** of the explored subgraph.
- `max_steps` is enforced per path head, not at single-step granularity inside a chunk.
- Beam overflow ranks by support only; there is no `queue` ranking option.
- Trace support with *column* labels follows consecutive column coordinates and does not detect a record boundary
  with consecutive global coordinates; header labels use local coordinates and do detect it.
- No measurement on a real ATB or `refseq33m` shard yet — only the 9 Mbp fixture.
- Untested: reverse-complement seed symmetry, a brute-force loss oracle, access-path invariance as an explicit
  test, and several later spec cases.

---

## 7. How to verify anything you suspect

```bash
# Release and Debug (asserts on). AppleClang needs the -Wno-error flag; GCC does not.
cd metagraph/build       && cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS=-Wno-error=unused-result .. \
                         && make -j metagraph unit_tests
cd metagraph/build_debug && cmake -DCMAKE_BUILD_TYPE=Debug   -DCMAKE_CXX_FLAGS=-Wno-error=unused-result .. \
                         && make -j unit_tests
# re-run cmake after ADDING a file: the source globs are not CONFIGURE_DEPENDS

./unit_tests --gtest_filter='*Walker*:*LabelOracle*:*Resolve*:MiniRefSeq*'   # 104 tests
./unit_tests                                                                  # full suite, 3767 tests

# the real-format fixture (~25 s, downloads 42 NCBI records; MiniRefSeq* skips without it)
../scripts/traversal/build_mini_refseq.sh ./mini_refseq ./metagraph \
    ../scripts/traversal/mini_refseq_accessions.tsv

# CLI and HTTP routes end to end
./test_venv/bin/python ../integration_tests/main.py --test_filter='*Traverse*'   # 12 tests
```

Current status: **104** traversal unit tests pass in Release and Debug; the **full 3767-test suite** passes in
Release (3708 in Debug before the last change); **12** integration tests pass. A mutation check was run on the
walker — reverting its semantic fixes fails 16 tests — so the suite has some teeth, but §5.2 shows it had a real
blind spot, so **please distrust the tests as evidence of correctness.**

---

## 8. How to report

For each finding: **severity** (blocker / major / minor), `file:line`, a **concrete failure scenario** (inputs or
graph shape → the wrong output), and a suggested fix. Please separate:

- **confirmed** — you traced the code and are sure;
- **suspected** — needs a run or a fixture to settle (say which experiment would settle it);
- **design disagreement** — the code does what it says, but it should say something else.

If you think the honest summary is "this is sound but unproven at scale", say that too — I would rather hear it
than a list of minor nits. And if a claim in §4 is simply unverifiable from the code as written, that itself is a
finding: it means the contract is not expressed where a reader can check it.
