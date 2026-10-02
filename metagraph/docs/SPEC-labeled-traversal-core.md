# Spec: label-consistent sequence extension in the MetaGraph core (`resolve` → `select` → `traverse`)

**Status:** draft v0.2 (v0.1 revised after four independent reviews) · **Branch:** `gr/labeled-traversal` (from `master` @ `1de02cec`)
**Scope:** C++ core library, `metagraph traverse` CLI, core server routes `POST /resolve`, `POST /traverse`, `GET /traverse/capabilities`.
**Companion:** `DESIGN-labeled-traversal-endpoint.md` (roadmap). This spec implements a prioritized subset of its
Phase A/B and deliberately generalizes some choices; deviations are listed in §13.
**Target index for testing:** the public `refseq33m` index — k = 31, BASIC DBGSuccinct (forward strand only, dummy
k-mers not masked), `RowDiff<BRWT>` annotation **with k-mer coordinates**, columns = NCBI taxids, labels reported as
accession headers through a `CoordToHeader` (`.seqs`) sidecar. A small index in exactly that format
(`scripts/traversal/build_mini_refseq.sh`) is the integration fixture. Source references are to this checkout
(paths relative to `metagraph/`).

---

## 1. Goals and non-goals

**Goals**

- **G1 — resolve.** For a (long) query sequence, report where it is present in the graph and which labels support
  which contiguous k-mer intervals; derive candidate seeds. No traversal.
- **G2 — select.** Turn candidates into an explicit, frozen seed list with a recorded, deterministic selection
  policy. Hit selection is separate from traversal.
- **G3 — traverse.** Extend each frozen seed outward in both directions under a declarative *strategy* that
  controls which labels may support the extension, what evidence a label must provide (k-mer presence or a
  consecutive coordinate trace), and how much label inconsistency is tolerated (a label-change cost and a per-arm
  loss budget).
- **G4 — tunable breadth.** Breadth is a function of the strategy, not only of the graph: a tight strategy prunes
  most structural branches and extends deep with a small result; a loose one fans out and hits size caps early.
  Every result reports the evidence an agent needs to tune it, including what a different setting would have done.
- **G5 — honest and bounded.** Deterministic, exact semantics; explicit per-label and per-path stop reasons;
  resource caps distinguished from semantic stops; no graph-sized allocation, whole-column extraction or global
  unitig scan per request.

**Non-goals for this branch:** async service, MCP tools, search-derived manifests (`hit_set_id`), path sampling,
resumable frontiers, motif or feature stops, multi-hit joint paths, bubble/allele *calling* (bounded bubble
*skipping* is in scope, §6.5), discovery of labels not named in the request during a walk, abundance/count
weighting, protein and DNA5 builds (code compiles for them; validation is DNA-only in v1).

## 2. Terminology

| Term | Meaning |
|---|---|
| node | Oriented graph node id as returned by `map_to_nodes_sequentially` on the walk graph (`sequence_graph.hpp:264`). `npos = 0`. |
| annotation key | Node id under which annotation is stored; annotation row = key − 1 (`annotated_dbg.hpp:50`). |
| regime | `basic`: key = node. `primary`: the walk graph is a `CanonicalDBG` wrapper over a primary graph (how PRIMARY indexes are loaded, `load_annotated_graph.cpp:126`), key = `get_base_node(node)`. `canonical`: a native graph storing both strands, key = the canonical node of the spelled k-mer. Detection: `dynamic_cast<const CanonicalDBG*>` → primary; else `get_mode() == CANONICAL` → canonical (DBGSSHash: key = `min(n, reverse_complement(n))`; others: `map_to_nodes` over the spelled k-mer); else basic. Any other wrapper is rejected as unsupported. |
| column label | An annotation column (`/column_labels`), e.g. a taxid on refseq33m. |
| header label | One indexed sequence inside a column, `(column, seq_id)`, named by its FASTA header (an accession on refseq33m). Requires a `CoordToHeader`. |
| label | A column label or a header label. Each is resolved once per request to a `LabelRef`. |
| support kind | `kmer`: the label carries every k-mer of the path. `trace`: the label additionally has consecutive k-mer coordinates in one indexed sequence along the path (trace-consistent). `trace` requires coordinates and is defined for the `basic` regime (coordinates in canonical/primary indexes carry no strand). |
| seed | An exact sequence (≥ k bp) whose every k-mer is in the graph, plus a non-empty set of seed labels each supporting every seed k-mer (under the strategy's support kind) — either given by the caller or derived from the seed itself (§6.1). |
| arm | The right (outgoing from the last seed k-mer) or left (incoming to the first seed k-mer) extension of one seed. |
| step | Entering one node on one path; adds exactly one base outside the seed. Step j (j ≥ 1) adds outward base index j−1. |
| extension DAG | Per arm, the explored paths as a DAG of segments: a split creates children, a reconvergence (§6.5) joins them. Without merging it is a tree. |
| segment | A maximal run of steps between splits/joins. |
| path | A root-to-leaf route through the DAG. |
| lineage | One label's continuation along a path, possibly renamed by switches (§6.3). |
| label state σ | For a position on a path: `label → entry(loss, branches, run, predecessor)` of labels still supporting it. |

## 3. Pipeline and why hit selection is separate

In a de Bruijn graph every k-mer is one node, so a sequence maps to exactly one oriented node path. "The sequence
as pointer" is well defined as a *location*. What multiplies is

1. **labels** — thousands of samples/contigs/accessions may carry the same path;
2. **partial support** — a label may support only some blocks of the query;
3. **alignment** — with `align=true` the hit sequence is a graph path that differs from the query (the caller passes
   that aligned sequence to `resolve`);
4. **occurrences inside a label** — several positions of one accession (visible with coordinates; `trace` support
   keeps them apart).

So a hit is effectively `(graph sequence, supported interval, label)` and choosing hits means choosing intervals
and labels:

```text
query sequence (+ optional candidate labels from a previous /search)
   └─ resolve ─▶ support profile + seed candidates                              [no traversal]
         └─ select(policy) ─▶ frozen seed list (seed_id, sequence, labels, population)   [recorded policy]
               └─ traverse(strategy) ─▶ per-seed extension DAGs + diagnostics   [never searches]
```

Because the traversal state carries a *set* of labels (§6.3), one seed with |S| = 1,240 labels is one traversal.
Labels are partitioned among successors where they diverge and re-joined where paths reconverge (§6.5).
Selection matters for bounding |S| and for choosing representative or complementary subsets, not for multiplying
walks.

## 4. `resolve` and `select`

### 4.1 Request

```json
{
  "sequence": "ACGT…",
  "labels": ["573", "NZ_STEQ01000045.1"],
  "discover": {"max_labels": 1000, "kind": "header"},
  "support": "kmer",
  "min_block_kmers": 1,
  "select": {"policy": "max_support", "max_seeds": 10, "min_block_bp": 31,
             "max_labels_per_seed": 1000, "label_order": "hash", "sample_seed": 1,
             "merge_overlapping": true},
  "run_format": "intervals",
  "bounds": {"max_query_bp": 100000, "time_budget_ms": 30000}
}
```

- Exactly one of `labels` (explicit, ≥ 1; column or header names) or `discover` must be present.
  `discover.kind` is `column` (default) or `header` (requires a `CoordToHeader`).
- `select` is optional; without it only the profile and candidates are returned.
- Unknown keys, wrong types and out-of-range values are rejected (strict validation, as in §5).

### 4.2 Semantics

Let L = |sequence|, k = graph k, k-mer starts i ∈ [0, L−k+1).

- `key(i)` = annotation key of k-mer i (`map_to_nodes`; strand-agnostic in primary/canonical regimes, forward
  strand only in basic). `in_graph(i) = key(i) ≠ npos`.
- Column label ℓ: `support(ℓ, i) = in_graph(i) ∧ ℓ ∈ row(key(i) − 1)`. Header label (c, s): in_graph(i) and some
  coordinate of column c at that row maps to sequence s (`CoordToHeader::map_single_coord`). With
  `support: "trace"` a header label additionally requires that consecutive supported k-mers have consecutive local
  coordinates (+1 per step); a coordinate jump ends the run (reported as `trace_break`).
- Three states per (k-mer, label): *absent from graph*, *in graph without ℓ*, *in graph with ℓ*. Labels not
  queried are never reported as unsupported.
- **Runs** are 0-based half-open intervals of k-mer starts `[a, b)`; base interval `[a, b + k − 1)`.
- **Discovery**: labels are collected from full rows (column kind) or row tuples (header kind) of all in-graph
  k-mers. If more than `max_labels` distinct labels occur, the top `max_labels` by supported k-mer count are kept
  (tie: column id, then seq_id) and `labels_truncated: {kept, total, min_kept_kmers, max_dropped_kmers,
  dropped_full_length}` is reported.
- **Seed candidates:** every maximal support run of each label with length ≥ `min_block_kmers`; labels with an
  identical run `[a, b)` are grouped. Candidate = `{kmer_interval, bp_interval, labels}`. Order: longer first,
  then more labels, then smaller `a`. Disconnected supported blocks of one label are separate candidates.

### 4.3 `select`

Deterministic, recorded in the output.

- `longest_first`: candidates in candidate order, skipping those shorter than `min_block_bp`, until `max_seeds`.
- `max_support`: choose `[a, b)` with `b − a ≥ min_block_kmers` maximizing `|{ℓ : a run of ℓ contains [a, b)}|`
  (sweep over run endpoints; tie: longer, then smaller `a`); repeat on the remaining k-mers until `max_seeds`.
  This is the right policy when many labels support near-identical intervals (run ends differ by a few k-mers).
  **It optimizes breadth, not length:** on the blaNDM fixture it returns the 231-k-mer prefix that all 25 carriers
  share rather than the 783-k-mer gene that 19 carry in full (a few SNPs break the rest). For the longest
  possible seed use `longest_first`, or raise `min_block_bp` to put a length floor under `max_support`.
- `explicit`: the caller passes `[{kmer_interval, labels}]`. The interval must lie inside one graph run; each
  label must have a run containing it (labels that do not are reported in `labels_not_covering`, and the seed is
  rejected if none remain).
- `merge_overlapping: true` (default): candidates whose intervals overlap are merged onto their common interval
  with the union of labels; `overlaps_with` lists seeds that still share k-mers.
- `max_labels_per_seed`: if a candidate has more labels, keep that many in `label_order`: `hash` (default; FNV-1a of
  label + `sample_seed`, unbiased w.r.t. column order), `column_id`, or `kmers_supported`. Dropped labels are
  reported as `labels_dropped: {count, digest, labels (if ≤ 1000)}` so the complementary subset can be run.

Output: a **frozen seed list** `[{seed_id, sequence, labels, label_population}]` where `sequence` is the block
substring, `label_population = {supporting_total, included, order, sample_seed}`, and
`seed_id = hex FNV-1a-64( release_id ∥ "\t" ∥ canon(sequence) ∥ "\t" ∥ join(len(ℓ) ":" ℓ for ℓ in bytewise-sorted labels) )`
with `canon(sequence)` = min(sequence, rc(sequence)) in primary/canonical regimes and the sequence itself in basic.
Plus counts `{candidates, eligible, selected}` and the normalized policy. The list is directly usable as
`/traverse` input.

### 4.4 Response (abridged)

```json
{
  "k": 31, "regime": "basic", "num_kmers": 4970, "support": "kmer",
  "graph_runs": [[0, 4970]],
  "labels": [{"label": "NZ_STEQ01000045.1", "kind": "header", "kmers_supported": 4782,
              "runs": [[0, 2100], [2131, 4970]]}],
  "labels_truncated": null,
  "candidates": [{"kmer_interval": [2131, 4970], "bp_interval": [2131, 5000], "labels": ["NZ_STEQ01000045.1"]}],
  "selection": {"policy": {...}, "counts": {...}, "seeds": [{"seed_id": "…", "sequence": "…", "labels": [...],
                "label_population": {...}}]},
  "capabilities": {"...": "§10.3"},
  "timing": {"elapsed_ms": 12}
}
```

`run_format: "rle"` emits runs as the existing presence-mask RLE (`encode_presence_mask`, `query.cpp:131`, moved
to a header).

## 5. `traverse` request

```json
{
  "release": "optional release id to pin",
  "seeds": [{"seed_id": "optional", "sequence": "ACGT…", "labels": ["NZ_STEQ01000045.1", "573"]},
            {"sequence": "ACGT…"}],
  "strategy": {
    "schema_version": 1,
    "direction": "both",
    "support": "kmer",
    "exhaustive": false,
    "labels": {
      "mode": "constrain",
      "max_labels_per_node": 64,
      "extra": [],
      "change_cost": {"model": "forbid"},
      "loss_budget": 0,
      "switch_on": "loss",
      "max_switch_sources": 64,
      "max_seed_labels": 1000,
      "seed_label_kind": "auto"
    },
    "branching": {
      "max_label_branches": 0,
      "min_successor_labels": 1,
      "min_successor_fraction": 0.0,
      "tip_window_bp": 0,
      "bubble_window_bp": 0,
      "on_reconverge": "merge",
      "hairpins": "skip",
      "max_splits_per_path": 64
    },
    "bounds": {
      "max_extension_bp": 5000, "min_live_labels": 1,
      "max_steps": 1000000, "max_live_paths": 1000, "max_paths": 10000,
      "max_output_bp": 2000000,
      "time_budget_ms": 30000
    },
    "frontier": {"order": "breadth_first", "beam_rank": "support", "on_overflow": "stop"},
    "output": {"detail": "full", "sequences": true, "label_lists": "delta", "label_runs": true,
               "profile_bin_bp": 100, "max_branch_events": 100, "continuation_bp": 1000, "timing": true},
    "annotation": {"access": "auto", "batch_kmers": 64}
  }
}
```

- **Strict validation.** Unknown keys, wrong types, out-of-range values and unsupported `schema_version` are
  rejected with a message naming the field; nothing is silently weakened. Defaults are filled in and the
  **normalized strategy** is echoed. If the server clamps a bound (§10.3) the echo carries
  `clamped: [{field, requested, effective}]`.
- `direction`: `both | left | right`. `support`: `kmer | trace` (`trace` rejected unless coordinates are indexed
  and the regime is basic).
- **`seeds[].labels` is optional.** Omitting it is the default and realizes the design note's
  `labels: {"permit": "per_hit"}` (design note line 183), which "explicitly means the resolved hit's annotation
  label, not all graph labels" (line 331): the permitted set is then **derived from the seed** — the labels that
  support every one of its k-mers (§6.1). An explicitly empty array is rejected naming the field (it reads as
  "nothing permitted", and the design note has empty permit sets invalid). A `change_cost` `table` is authored by
  label *name* over the request order and is therefore rejected together with a derived seed.
- `labels.max_seed_labels` (default 1000, 1 … 100 000; the server may clamp it lower, see below): cap on a
  **derived** set. It bounds the traversal state, not the
  discovery cost — the labels above it are found in the same pass either way and are reported, never dropped
  silently (§6.1, §7.1). Ignored when `labels` is given (that list is the caller’s own cap). Under
  `exhaustive` a derived set over the cap is **refused** per seed (a `SeedDerivationError` naming the knob),
  not cut: every walk of a dropped carrier would be missing from a trie that promises no walk is dropped.
- `labels.seed_label_kind`: `auto` (default) | `column` | `header` — what a *derived* label denotes. `auto` is
  `header` when the index has a `CoordToHeader` (so a derived label is an indexed sequence / accession) and
  `column` otherwise. `header` without a `CoordToHeader`, and any kind needing coordinates on an annotation
  without them, is rejected rather than downgraded. Ignored when `labels` is given (each name resolves by itself).
- `labels.extra`: additional permitted labels (switch targets). Rejected when unreachable: model `forbid`, or every
  `cost(seed label → extra) > loss_budget`. `extra` applies on top of a derived set, and an `extra` label the
  derivation also produced is rejected as a duplicate seed label.
- `labels.change_cost.model`: `forbid` (default) | `constant` (`value`) | `table`
  (`entries: [[from, to, cost], …]`, `default`: `"forbid"` or a number) | `hll` (increment I8; `scale`,
  `denominator`). Costs ≥ 0; `cost(ℓ → ℓ) = 0`; `forbid` = +∞.
- `labels.loss_budget` ≥ 0, per arm. `labels.switch_on`: `loss` (default; a switch into ℓ′ is allowed only from a
  label that does **not** continue at v — switching is forced, never free-riding) | `any`.
  `max_switch_sources` (default 64, or `"unlimited"`): for pairwise models, only the cheapest that many sources
  are considered per step. Under `exhaustive` its default is `"unlimited"` and, with a `table` cost, a number is
  **rejected**: a target label reachable only from a cut source is not entered, and if it was the only label on
  that successor the walk into it would be pruned with no path-level reason.
- `branching.*`: see §6.4–§6.5. `max_label_branches` is a **depth** along a root-to-leaf path (a lineage may reach
  up to d^b leaves, d ≤ number of successor characters), not a total.
- `bounds.*` scopes are in §6.8. `frontier.*` in §6.8. `output.detail`: `summary | tree | full` (§7).
- **`labels.mode`**: `constrain` (default; everything above) | `annotate` (§6.9: every structural successor is
  followed, the labels present are recorded, nothing is filtered). In `annotate` mode `seeds[].labels`,
  `labels.extra`, a `change_cost` other than `forbid`, a `loss_budget`, the quorum knobs, `bounds.min_live_labels`
  and `support: trace` are **rejected** naming the field (there is no permitted set for them to act on), and
  `branching.max_label_branches` / `max_splits_per_path` may only be `"unlimited"` (there is no lineage to count).
- **`labels.max_labels_per_node`** (1 … 100 000; default 64 in `constrain` mode and `max(64, max_seed_labels)` in
  `annotate` mode, so that the structural oracle can by default record what the constrain side admits as a
  derived set): cap on every per-node label list in the output — the recorded sets of `annotate` mode, the
  lists on `blocked` / `hairpin` events, the root's boundary list and the per-branch lists of the trie view in
  both modes. The true count is always reported beside a cut list (§7.4); a cut is never silent.
- **`exhaustive`** (default `false`): the preset of §6.9 — unlimited branch allowance, `on_reconverge: keep`, no
  split limit, no quorum, no beam. It changes the *defaults* of those knobs; a knob given explicitly with a
  conflicting value is **rejected** (400) rather than overridden.
- `branching.max_label_branches` and `branching.max_splits_per_path` accept an integer or the string
  `"unlimited"`, and the normalized echo uses `"unlimited"` where it applies.
- Values above are placeholders, not benchmarked defaults.

## 6. Traversal semantics

### 6.1 Seed validation

1. `|sequence| ≥ k`; every character is in `graph.alphabet()` minus `'$'` after the build's case mapping
   (DNA builds: upper-case; `N` is invalid in DNA builds). The first invalid position is reported.
2. `nodes = map_to_nodes_sequentially(walk graph, sequence)`; every node ≠ npos, otherwise the seed is rejected
   with its `graph_runs` (a partially present sequence is never treated as one seed).
3. Every seed label must resolve (else the request is rejected naming the label). Each label must support every
   seed k-mer under `support`; labels that do not are **dropped** with `{reason: seed_unsupported, runs}`. A seed
   with no remaining label is rejected. `validated_seed_id` is computed over the remaining labels; a supplied
   `seed_id` that differs from the recomputed one yields `seed_id_mismatch: true`.
4. **Derived seed labels (the `permit: per_hit` default, §5).** When `labels` is omitted, the permitted set is the
   seed’s own: the labels that support **every** k-mer of the seed, of the kind `seed_label_kind` names, in
   ascending `(column, seq_id)` order. The derived set is the **intersection of the very rows step 3 reads
   anyway** — one pass over the seed’s annotation rows, no second scan, and nothing outside the seed is ever
   looked at: deriving costs no extra annotation reads over validating a given list. (With header labels the
   coordinates of the still-live columns are mapped to sequence ids per k-mer; after the first k-mer only columns
   that still carry a candidate are touched, so the work shrinks with the intersection.) The intersection is taken
   incrementally and the seed is rejected as soon as it is empty. Rows are reconstructed in sub-batches of at most
   64 k-mers (the batch is what makes a long seed fast), so an empty intersection stops the scan within one
   sub-batch rather than at the exact k-mer; `rows_requested` counts the rows actually CONSUMED, which is what the
   derivation read, not what the batch fetched.
   The pass starts from the **cheapest row among the seed's first 64 k-mers**, not from k-mer 0: that one row is
   the only one no intersection has narrowed yet, so starting at k-mer 0 would make the peak cost an accident of
   where the caller cut the seed. The window is 64 k-mers **whatever `annotation.batch_kmers` is** (only the later
   sub-batches follow the knob): the cheapest-row choice and the guard below decide whether a seed is accepted,
   and acceptance must not depend on a fetch-size setting (§6.8). The order does not change the result (an
   intersection is commutative and the set stays ascending by `(column, seq_id)`), and a seed whose cheapest row
   in that window still has more annotation entries (columns plus k-mer coordinates — an upper bound on its
   labels) than the derivation will materialise is refused up front (a seed of exactly k bp has no intersection
   at all, so its "derived" set would be a whole annotation row). The derivation also watches
   `bounds.time_budget_ms` itself — it runs before the walk, i.e. before the only other place the clock is read.
   `max_seed_labels` then truncates the ordered set; the truncation is reported as
   `{labels_from_seed, labels_supporting_total, labels_dropped, labels_dropped_digest}` (FNV-1a-64 over the cut
   names), and truncated labels are **not** in `dropped_labels`: they support the seed, they were only not taken.
   For an **explicit** list `labels_supporting_total` is `|labels after dropping|` and the other three are
   `false / 0 / ""`, so `labels_supporting_total > |seed.labels|` means "truncated" in both cases and
   `labels_supporting_total` never reads as "nothing supports this seed".
   A derived list is **resubmittable verbatim** as an explicit one, and three cases that would break that promise
   are refused instead of echoed: a `header` set whose names are not distinct (a FASTA header is unique per
   column, not per index, and the set is deduplicated by `(column, seq_id)`, while the explicit path rejects
   duplicate names), a `header` whose name does not **resolve back** to its own `(column, seq_id)` (the explicit
   path resolves a name to the first column holding it, which may be a column that does not carry the seed), and
   the too-wide first row above. Each refusal names the header and the two ways out (`seed_label_kind: column`
   or an explicit list).
   The cap exists to bound the traversal state (§6.3 carries a label set per path position), not the discovery.
   Under `support: trace` the set is derived from k-mer **presence** and then held to step 3 unchanged, so a
   derived label without a coordinate-consecutive occurrence of the seed appears in `dropped_labels` with
   `seed_unsupported` exactly as a named one would. Derived labels take the place of the given ones everywhere:
   `label_dict`, every output label id, and the request indices of a `table` cost model (§9) — which is why a
   table, authored by name, requires an explicit list.
   Note what the derived set *is*: the seed’s **carrier set** — every label whose annotation supports the whole
   seed — not one distinguished hit. Under the default `change_cost: forbid` that is the design’s `path_common`
   over exactly those carriers (§6.3), which is the honest reading of "the resolved hit’s annotation label" when a
   sequence is carried by several records; a caller who wants `within_label` for one specific record still names
   that one label (and `/resolve` + `/select` remain the way to freeze a reproducible subset).
   **Per-seed failure.** Deriving can fail for reasons the caller cannot act on, because it named no labels at all:
   nothing carries the seed in full, the derived names are ambiguous, the first row is too wide, or the time budget
   ran out mid-derivation. Those do not fail the request: the seed's entry in `results` is
   `{"seed": {"seed_id", "length_bp", "labels_from_seed": true}, "error": "<reason>"}` — no `arms` — and the other
   seeds are traversed as usual. A client distinguishes the two shapes by the presence of `error`. Everything that
   *is* the caller's own doing still fails the whole request with HTTP 400: a malformed seed (too short, invalid
   character, not fully present in the graph), an unknown or duplicate **explicit** label, and every strategy error.
5. The seed is never rewritten.
6. The permitted universe is `P = seed labels ∪ extra`. Duplicate `seed_id`s in one request are computed once and
   reported with `duplicate_of`.

### 6.2 Arms, steps, orientation

- **Right arm**: successors of the last seed node via outgoing edges; **left arm**: predecessors of the first
  seed node via incoming edges, both in the seed's orientation. Left-arm sequences are reported in natural
  orientation.
- Implementation per regime (§8.4): `primary` uses the 3-argument `CanonicalDBG::call_*_kmers(node, spelling_hint,
  cb)` with the node's k-mer from the arm's rolling buffer; whenever the walk graph is or wraps a `DBGSuccinct`
  the request owns a `NodeFirstCache` and the left arm calls `cache.call_incoming_kmers` (basic/canonical) or the
  cache is attached to the per-request `CanonicalDBG` clone (primary).
- Successor characters `'$'` (`boss::BOSS::kSentinel`) are ignored; successors are processed in ascending
  character order. Degree functions are never used for decisions (they count dummy k-mers on unmasked
  DBGSuccinct, `dbg_succinct.cpp:557-640`).
- **Positions.** Step j adds outward base index j−1 (seed coordinate `|seed| − 1 + j` on the right arm, `−j` on the
  left). Runs are half-open `[from_bp, to_bp)` over outward base indices. `label_end at_bp = x` means the label's
  last supported step is x (whether a base x+1 exists is given by the reason). On the left arm, outward index i of
  a segment maps to `sequence[length_bp − 1 − (i − from_bp)]` of the natural-orientation segment string.

### 6.3 Label state and the loss recurrence

At the seed boundary of each arm: `σ₀ = {ℓ ↦ (loss 0, branches 0, run(ℓ, from 0, entered_by seed), pred none)}`
for every validated seed label.

Stepping from node `u` with state σ to a successor `v` with admissible labels `A(v) = labels(v) ∩ P` (under
`support: trace`, a header label is in `A(v)` only if it has a coordinate equal to its coordinate at `u` + 1 for
the right arm, − 1 for the left arm; otherwise the lineage ends with `trace_break`):

```text
for each ℓ' in A(v):
    stay   = σ[ℓ'].loss                                   if ℓ' ∈ σ else +∞
    sources = { ℓ ∈ σ : ℓ ≠ ℓ' and (switch_on == any or ℓ ∉ A(v)) }   (pairwise models: cheapest max_switch_sources)
    switch = min over sources of σ[ℓ].loss + cost(ℓ → ℓ')
    loss'  = min(stay, switch);  keep ℓ' iff loss' ≤ loss_budget
    pred(ℓ') = ℓ' on stay, else the argmin source
σ_v = { kept ℓ' ↦ (loss', branches(pred), run', pred) }
```

- Tie-break (deterministic): `(loss, stay before switch, branches of the predecessor, column id / seq_id of the
  predecessor)`.
- `run'` is the predecessor's run on `stay`; on a switch a new run record `{label ℓ', from_bp, entered_by switch,
  from ℓ, cost}` is created, linked to the predecessor's run (a persistent list, memory ∝ number of switches).
- For `constant`, `switch` is `(min loss over sources) + c` computed with first and second minimum: O(|σ| + |A|).
- `v` is admissible for the path iff `σ_v ≠ ∅`, the edge `u → v` is unused on this path (§6.6), and `v` is not a
  hairpin step when `hairpins = skip` (§6.6).

**Exactness.** For a fixed path, `max_label_branches = ∞` **and a non-binding `max_switch_sources`**, `loss'` is
the exact minimum over label assignments (costs ≥ 0, so budget pruning is safe). Two documented approximations:

- Branch-limit pruning (§6.4) is greedy — applied after each successor has fixed one predecessor per label, so it
  is not a constrained optimum. Tests assert optimality only without branch limits.
- **`max_switch_sources` truncates the source list by `(loss, branches, label)`, not by the cost into the
  particular target**, so with a `table` cost it can miss the cheapest switch: if A→E costs 2 and B→E costs 0.5
  from equal starting losses, `max_switch_sources: 1` reports 2 and `2` reports 0.5. Raising the cap restores
  exactness; the approximation is reported, not hidden.

**Special cases and their meaning.**

| Cost model / budget | Equivalent | Evidence a surviving label provides |
|---|---|---|
| `forbid`, single seed label | design `within_label` | That label carries every k-mer from the seed to here |
| `forbid`, several labels | design `path_common`, with per-label splitting | Each surviving label carries every k-mer from the seed to here |
| `constant c`, budget B | ≤ ⌊B/c⌋ label switches | A chain of labels with ≤ that many switches covers the path |
| `table` / `hll` | "close enough labels" | Each switch is priced by label similarity; total ≤ B |
| cost 0 for all pairs | design `any_of` | Every k-mer has some permitted label — exploratory only |
| any of the above with `support: trace` | design `sequence_trace` | Additionally, each run follows one indexed sequence occurrence |

A path with switches is a **candidate sequence supported by a chain of labels**, never evidence that one sample
contains it; every switch and its cost are reported, and switched-in runs are distinguished from seed runs.

### 6.4 Branching is per lineage

At node `u` with state σ, let `v₁ … v_m` be the admissible successors with states `σ₁ … σ_m` (computed by §6.3).

- A **source entry** ℓ ∈ σ is **ambiguous at u** if entries whose predecessor is ℓ occur in two or more σ_i.
- An ambiguous source with `branches + 1 ≤ max_label_branches` continues on each such successor; every derived
  entry gets `branches + 1`. Otherwise every entry derived from ℓ is removed from all σ_i and each removed target
  ℓ′ is re-minimised over the remaining sources (dropped if none is finite and within budget); ℓ ends at u with
  reason `branch`. The branch event records the successor characters and per-successor label counts.
- **Quorum.** A successor whose `|σ_i|` is below `min_successor_labels` or below `min_successor_fraction · |σ|` is
  not followed; its entries end with reason `minority` (reported like `branch`, with counts). A path whose `|σ|`
  falls below `bounds.min_live_labels` ends with `below_min_labels`. Both are semantic stops.
- A non-ambiguous source continues on its single successor **without** counting a branch, even if the path splits
  (a **divergence**: the sources of the successors are disjoint).
- Successors whose state became empty are not followed; if two or more remain, the path splits. A path that
  would split beyond `max_splits_per_path` ends with `split_limit` (semantic).

**Consequences.** With `forbid` and `max_label_branches = 0` each arm has at most |S| leaves: a structural branch
costs nothing unless some permitted label actually supports both sides. **Equality with the union of independent
single-label walks holds only under `on_reconverge: keep`**: merging unions the merged lineages' used-edge sets
(§6.6), which is conservative, so a label can be stopped by `edge_reuse` on an edge only another label's
pre-merge route had used. The reference-walker test therefore runs with merging off. Under any finite change cost, divergences
followed by re-entry via switch can multiply paths without consuming branches; `max_splits_per_path`, the quorum
knobs and reconvergence merging (§6.5) bound that. Note that under `forbid`, a bubble whose two sides are both in
one label (within-sample variation) is an ambiguous branch for that label; §6.5 is the remedy.

### 6.5 Tips, bubbles, reconvergence and hairpins

- **Tip window** (`tip_window_bp > 0`): before declaring a source ambiguous at u, each alternative is explored
  structurally (graph only, counted toward `max_steps`) for up to `tip_window_bp`; an alternative whose permitted
  continuation ends (dead end or label lost) within the window is a **tip**: it does not count toward ambiguity,
  is not followed, and a `tip` event `{at_bp, char, length_bp}` is recorded.
- **Bubble window** (`bubble_window_bp > 0`): if all ℓ-admissible alternatives reconverge at one node r within
  the window, with ℓ on every node of every side and no further branching inside, no branch is counted; the path
  continues from r along the lowest-character side; a `bubble` event `{at_bp, alleles: [{sequence, labels}],
  length_bp}` is recorded. Bubbles that fail these conditions are ordinary ambiguous branches.
- **Reconvergence** (`on_reconverge`): per arm, a map `oriented node → (segment, bp)` of first arrivals, bounded by
  `max_steps` (nothing graph-sized). `merge` (default): when a live path enters a node another path of the same
  arm holds at the same outward distance, the paths join into one DAG node with `parents: […]`; σ is the per-label
  union keeping the min `(loss, branches)` entry when a label is on both (exact Viterbi minimum on the DAG); each
  label still enters through exactly one parent, so per-label routes stay reconstructible; used-edge sets become
  the union (conservative). A `reconverge` event `{at_bp, segments, bubble_len_bp}` is emitted and counted in
  `growth`. Arrivals at different distances are not merged (reported as `revisit` events). `keep`: no merging.
- **Hairpins.** A step through a self-reverse-complementary (k+1)-mer (`u → rc(u)`; for even k also into/out of a
  palindromic k-mer) is a hairpin step. `skip` (default): not admissible, never counts toward ambiguity, a
  `hairpin` event is recorded. `follow`: followed and flagged. (In canonical/primary regimes every node whose
  (k−1)-suffix is an RC palindrome has such a successor carrying exactly its own labels; without this rule every
  label would stop there.)

### 6.6 Edge reuse, seed re-entry and cycles

- **Edge identity** is the (k+1)-mer of the step: in primary/canonical regimes canonicalized (`min(s, rc(s))`,
  with the orientation bit kept), in basic the forward (k+1)-mer. Encoding: 2-bit packed when the alphabet is
  DNA and k+1 ≤ 32, otherwise a 128-bit hash with full-string comparison on collision.
- A path uses each edge at most once. Storage: one map per seed and arm `edge_key → small vector of (segment,
  offset, orientation)`; a step is blocked iff some entry's segment is an ancestor-or-self of the current segment
  (precomputed depth and parent pointers; for merged DAG nodes, any parent). No graph-sized arrays.
- **Seed re-entry.** A step whose target is any seed node (either orientation in canonical/primary) ends the path
  at that node with `rejoined_seed{seed_offset, orientation}`, regardless of seed length; the re-entered bases are
  reported as `overlap_bp` and not counted in `extension_bp`. Interpretation: a circular molecule (offset 0 on the
  right arm / the last offset on the left arm, same orientation) or a repeat shared with the seed. Arms do not
  share used-edge sets; a circle is reported by both arms.
- Precedence when several apply: `rejoined_seed > edge_reuse_rc > edge_reuse`.
- **Blocked continuations are always reported.** Whenever a successor that would carry at least one entry of σ is
  inadmissible only because of edge use (or a skipped hairpin), a `blocked` event `{at_bp, char, reason, labels}`
  is emitted whether or not the path continues, and counted per bin as `blocked_repeat`. A flank that looks unique
  but passed a collapsed repeat is therefore visible as such.
- **Cycles.** Each (k+1)-mer is used at most once per path, so a cycle can be traversed only as far as returning
  to a used edge; a loop-then-exit continuation exists only if the exit leaves from the node where the cycle was
  entered. Cycle junctions and self-loops are ambiguous branches for labels on both sides (so with
  `max_label_branches = 0` the label stops at the junction with `branch`; with 1 the loop is unrolled once).
  One unrolling is not an estimate of copy number.
- **Homopolymers (a consequence worth knowing).** A run A^n with n ≥ k+1 is a self-loop on the node A^k, and every
  node whose (k−1)-suffix is A^(k−1) has the exit A^(k−1)·Z[0] as a successor. For a record X·A^(k+3)·Z the trie
  therefore holds A^(k−1)·Z, A^k·Z and A^(k+1)·Z — and never the record's own A^(k+3)·Z, under k-mer *or* trace
  support (the loop edge may be used once; the coordinates force the record's path onto it a second time). Pinned
  by `TrieCases.HomopolymerSelfLoop` and `TrieCasesTrace.HomopolymerLimitation`.

### 6.7 Termination reasons

**Per label** (the name leaves σ on this path; a name may re-enter later by a switch, which starts a new run):

| Reason | Meaning |
|---|---|
| `dead_end` | `u` has no structural successor in this direction (graph tip) |
| `label_lost` | successors exist but none carries the label and no finite switch from it exists |
| `loss_budget` | a switch from the label was possible but only above the budget; `needed_budget` is reported |
| `superseded` | a within-budget switch from the label existed but a cheaper predecessor was chosen for every target |
| `switched` | the lineage continued under another name (`to`, `cost`) |
| `branch` | ambiguous with `max_label_branches` exhausted |
| `minority`, `below_min_labels`, `split_limit` | quorum / split stops (§6.4) |
| `trace_break` | `support: trace`: the coordinate did not continue (+1/−1) |
| `edge_reuse`, `edge_reuse_rc`, `rejoined_seed`, `hairpin` | §6.6 / §6.5 |
| `beam_pruned`, and every resource reason | the path ended for a path-level reason; the label is **censored** with that same reason |

**Per path**: ends when no admissible successor remains (its `end_reasons` count the label reasons above), with
the semantic reason `max_extension_bp` (the requested radius; the domain is complete, the biology may continue),
or with a **resource** reason: `max_steps`, `time_budget` (per seed: end all live paths of both arms), or
`max_live_paths`, `max_paths`, `max_output_bp` (per arm: end that arm only). `beam_pruned` marks pruned paths.

Each arm reports `status: complete | truncated (a resource reason) | pruned (beam)`, `frontier_remaining:
{live_paths, live_labels, exact}` and `cap_trigger: {reason, arm, at_bp, segment, live_paths, live_labels,
exact}` when truncated. `exact` is `false` when a head counted carried a list cut by `labels.max_labels_per_node`
(`annotate` mode): `live_labels` is then a lower bound. The same flag sits on every `growth[]` bin
(`exact`); in `constrain` mode it is always `true` (the live state is never cut). `max_output_bytes` / `max_events` are **delivery** truncation (`delivery: {truncated, next_cursor}`),
separate from exploration status. A capped or pruned result is never evidence of absence.

### 6.8 Exploration order, bounds, determinism

- One priority queue per arm. Keys: `breadth_first` = `(extension_bp, path_id)`; `lowest_loss_first` =
  `(min loss in σ, extension_bp, path_id)`; `most_supported_first` = `(−|σ|, min loss, extension_bp, path_id)`,
  where in `annotate` mode |σ| is the number of labels recorded at the head node (§6.11).
  Arms alternate by `extension_bp` so both advance fairly.
- A popped path advances along its structurally unbranched run in chunks of up to `batch_kmers` nodes (labels
  batch-fetched, §8.3). **Cap and overflow decisions are made in queue-key order at single-step granularity:** a
  chunk is truncated at the step where a cap trips, so `batch_kmers` never changes which steps are taken.
- Scope table: `max_extension_bp` per path; `min_live_labels` per path; `max_steps`, `time_budget_ms` per seed;
  `max_live_paths`, `max_paths`, `max_output_bp` per arm; `max_output_bytes`, `max_events` per seed (delivery);
  request-level `time_budget_ms` (optional, `request_time_budget_ms`) marks seeds not started as `not_started`.
- **Overflow** (`on_overflow`): `stop` ends the arm's live paths with that reason; `beam` keeps the best paths by
  `beam_rank` (`support` = `(−|σ|, min loss, extension_bp, path_id)`, or `queue`) within the arm and ends the
  others with `beam_pruned`; pruned paths are listed in `beam_pruned: [{path, at_bp, n_labels}]`.
- **Determinism.** The `result` object (seeds, arms, segments, events, label_summary, growth, branch_events,
  step/enumeration/rows-mapped counters) is byte-identical across runs, across seed permutations and when the
  request is split into one seed per request, for the same graph, annotation, seeds and normalized strategy,
  provided no time budget tripped. All wall-clock values and cache-hit counters live in a separate `timing`
  object excluded from that contract (`output.timing: false` omits it). The annotation cache is per seed.

### 6.9 Label modes and the `exhaustive` preset (the trie oracle)

The walk's own output cannot validate its label filtering: it *is* the filtered result. An exhaustive trie built
**without** the label constraint, with the labels merely recorded, shares no label logic with §6.3–§6.4 and is
therefore usable as independent evidence for it. That is what the second label mode is for.

- **`labels.mode: constrain`** (default): §6.1–§6.8 as written — a successor is admissible iff some permitted
  label survives on it.
- **`labels.mode: annotate`**: admissibility is **purely structural**. Every non-`'$'` successor is followed,
  subject only to the per-path edge-reuse rule, hairpin handling and seed re-entry of §6.5–§6.6. There is no
  permitted set, no label state, no loss budget, no switching, no quorum and no branch limit: `seeds[].labels`
  must be omitted and the label machinery is rejected (§5). The labels present at every node are **recorded**
  instead — `label_dict` is filled in the order labels are first met, each segment carries the sets present
  along it as runs (§7.4), and `labels_start` / `labels_end` are the sets at the segment's entry node (the seed
  boundary for the root) and last node. A path ends with a **path-level** reason only — `dead_end`,
  `edge_reuse`, `edge_reuse_rc`, `rejoined_seed` (a successor that was only a skipped hairpin is a `dead_end`
  with its `hairpin` event), `max_extension_bp` or a resource cap — and `end_labels` / `end_reasons` are empty,
  because no lineage is tracked. Splits are never `ambiguous` (a lineage notion). `label_summary` is computed
  from the recorded sets: `reach_bp` is the furthest position a label was recorded at, `direct_bp` its
  continuous presence from the seed boundary along *some* route (exact when no list was cut, a lower bound
  otherwise), `runs` is empty.
  **Cost.** Every node costs a **full row** (or tuple row, for `header` labels: every coordinate is mapped so that
  the true count is known). Nothing narrows the read, which is why this is a verification tool at small radius
  and not a search primitive; the per-node cap `labels.max_labels_per_node` bounds the *output* and the
  dictionary, not the read. The label kind recorded is `labels.seed_label_kind` (`auto`: `header` with a
  `CoordToHeader`, `column` otherwise).
- **`exhaustive: true`** (both modes): unlimited `max_label_branches`, `on_reconverge: keep` (the result is a
  **trie** of walks, not a DAG), unlimited `max_splits_per_path`, no quorum (`min_successor_labels: 1`,
  `min_successor_fraction: 0`, `min_live_labels: 1`), `on_overflow: stop`. The preset sets the defaults of
  exactly these knobs and **refuses** a request that sets any of them to something else — asking for the
  exhaustive trie and receiving a beam-pruned result would be the failure this flag exists to prevent, so the
  conflict is a 400 naming the knob and the value the preset requires. **No pruning of any kind happens under
  the preset**: the only legitimate stop besides a semantic end is "end of the last complete level" (§6.10).

The verification contract this enables, for a seed S, radius R and permitted set P, all exhaustive:
`T = trie(S, R, annotate)`, `E = {paths p of T : for some l in P, every node of p carries l}` (computed by
filtering the raw recorded sets — no walker label logic), `A = trie(S, R, constrain, P, cost forbid)`; then
`claims(A) == E` where a claim is a (walk, label) pair — `E` holds, per permitted label, its maximal
label-consistent walks, and `claims(A)` the walker's label ends, so that a label ending inside a walk that other
labels continue is checked at its own end position (comparing leaves alone would never see it) — and for any
tuned run `leaves(tuned) ⊆ leaves(A)` with every omission carrying a reported reason. Comparisons are
restricted to `min(complete_to_bp)` of the runs compared (§6.10).

**The third party.** T and A share the graph, the annotation and the walker's structural code, so their
agreement cannot catch an error common to both (a successor the graph fails to enumerate, a dummy k-mer followed,
a label the annotation attaches to the wrong node). The tests therefore add a model that shares nothing with the
index: `tests/graph/traversal/test_trie_reference.hpp` builds, from the records a fixture is made of, the k-mer
sets per label and enumerates the walk rule of §6.10 over the strings — the structural walks with the labels
present at every node, each label's maximal walks with the reason each ends, the §6.3 recurrence for a constant
switch cost and a budget, and the records' own continuations for `support: trace`. Every fixture of
`test_trie_cases.cpp` is checked three ways (records ⇔ T ⇔ A, with end reasons), plus: the run over P is the union
of the single-label runs (the complete list of single-label walks), a permitted set derived from the seed runs
like the explicit list, and the recurrence at budget 0 gives the forbid leaves. On the real mini-refseq index the
same model is built from the 42 source FASTA records (9.3 Mbp, 2.8 s packed) and the contract holds on the whole
blaNDM seed and on windows cut from records, under k-mer and under trace support
(`MiniRefSeq.TrieContractAgainstTheSourceRecords`).

### 6.10 Size bounds and the completeness guarantee (`complete_to_bp`)

The structural trie can explode, so the bound is the crux, and **its order decides whether the result proves
anything**: a depth-first walk cut by a size cap yields some deep paths and some wholly unexplored ones — complete
for nothing — while a breadth-first walk cut by the same cap yields a trie in which every walk of length ≤ d has
been fully explored, which is a usable oracle, just a shallower one. Exploration here is level-synchronous (§6.8),
so no new order is needed; what the bound needs is a **per-level completion boundary and an honest report of it**.

- The size caps are the existing ones: `max_steps` (nodes entered, per seed), `max_output_bp` (bases emitted, per
  arm), `max_paths` (leaves, per arm), `max_live_paths` (frontier, per arm), `time_budget_ms`. Radius is what the
  caller pushes; size is what is capped.
- A cap trips **between two heads of one level**: the heads expanded before it have all their children, the heads
  after it have none. That level is **partial and does not count**. Every arm reports **`complete_to_bp`** = the
  extension depth up to which *every* admissible walk is present: for every n ≤ `complete_to_bp`, every walk of
  n bases from the seed boundary that obeys the walk rule below is in `segments`, and every path that ended
  before `complete_to_bp` ended for the reported semantic reason. Walks longer than that may be present (the
  partial level's children, honestly ended with the cap reason) but nothing is claimed about them.
- `status: complete` ⇔ `complete_to_bp == bounds.max_extension_bp`; otherwise `truncated` (or `pruned`) with
  `cap_trigger` naming the cap. Heads that have already reached the radius when the time budget or a seed-level
  cap trips are ended as `max_extension_bp`, not censored, so an arm whose every walk is present is reported
  complete; a beam never prunes heads at the radius.
- **The walk rule**, stated in every response as `walk_rule`, is what makes "all walks" a finite set in a graph
  with cycles and therefore defines what `complete_to_bp` quantifies over: a walk uses no (k+1)-mer edge twice on
  its own path (canonical (k+1)-mers in canonical/primary regimes, each use keeping its orientation), enters no
  seed node (of either strand in those regimes), and takes no hairpin step (`hairpins: skip`) or takes them
  flagged (`follow`); in `constrain` mode it is additionally followed only while some permitted label supports
  every node of it. The statement is exact for the index it is issued on: for even k a step into or out of a
  self-reverse-complementary k-mer node is a hairpin too, and when a (k+1)-mer does not pack into 64 bits
  ((k+1) · bits per symbol > 64, e.g. k ≥ 32 on DNA) edges are identified by a 128-bit FNV-1a hash of the
  (k+1)-mer, so a collision blocks a step as a reuse — the finite set of walks is then slightly smaller than
  "no (k+1)-mer twice". In a `basic` graph (one strand) no step is a hairpin and the statement says so. **Under
  `on_reconverge: merge` the certificate is weaker and says so too:** a merge unites the edge histories of the
  routes it joins (conservative), so a walk admissible on its own path can be blocked by an edge another route
  used, and `complete_to_bp` then quantifies over the walks admissible under the *united* history — fewer than
  the per-path set. The per-path claim holds only with `keep`, which the `exhaustive` preset forces and which is
  the default in `annotate` mode (the mode that promises "the whole trie"); in `constrain` mode merging stays the
  default and the response's `walk_rule` carries the qualification.
- The depth at which the structural trie stops is itself a measurement: "the unconstrained trie explodes at 40 bp
  here" is a reportable property of the locus.

### 6.11 Label-free exploration: `annotate` with a beam

A complete trie is bounded by **breadth**: every walk up to `complete_to_bp` must be present, so in read data —
where every sequencing error and every strain difference is a fork — the label-free trie dies within tens to a
few hundred bases whatever the size cap (the walk count grows geometrically with depth; raising the cap buys depth
logarithmically). What goes deep is giving up completeness: `labels.mode: annotate` **without** `exhaustive`, with
`frontier.on_overflow: beam` and `bounds.max_live_paths: N` (N = 1 is a single guided walk). This is a legitimate
tool for an agent exploring the graph *locally* before committing to labels, and the contract for it is:

- **What the result is.** A beam path is a set of **candidate bases with per-node support**, not a sequence any
  sample is claimed to contain. At every fork the beam kept some heads and dropped the rest, so the spelled walk
  may switch, at any node, from one set of carriers to another. The response says so: `status: pruned`,
  `cap_trigger.reason: beam_pruned`, `walk_rule` (which, without `exhaustive`, makes no completeness claim beyond
  `complete_to_bp`, the depth before the first pruning), and `label_mode: annotate`.
- **What the evidence is.** The labels present at every node are recorded (`segments[].label_sets`, true counts),
  so the agent *sees* the support change along the path — the carrier set shrinking, a new set taking over —
  without constraining by it. Per label, `label_summary.direct_bp` is how far that sample follows the spelled
  walk contiguously from the seed boundary (exact when no list was cut); that, not the walk's length, is the
  per-sample claim. A stretch supported by no label at all is recorded as an empty set.
- **How the beam ranks.** `order: most_supported_first` keeps the heads with the most labels *recorded at the head
  node* (in `constrain` mode: the most labels alive), ties broken by `extension_bp` and `path_id`, so the result is
  deterministic and a width-1 beam follows the majority continuation at each fork — the most-supported local
  path, not an arbitrary one. `breadth_first` keeps the earliest-created heads (no preference) and is the default
  only because it is the completeness order; for exploration pass `most_supported_first` explicitly.
- **How to use it.** Explore with a beam → read where the support changes → either constrain
  (`labels.mode: constrain`, naming the carriers seen) or re-seed from a `continuation` → extend. The design
  note's "probe cheaply, read the diagnostics, retune" loop; the beam is the cheap probe that reaches far, the
  constrained walk is the one whose output is a claim.

### 7.1 Per seed (`detail: full`)

```json
{
  "seed": {"seed_id": "…", "validated_seed_id": "…", "seed_id_mismatch": false, "length_bp": 4870,
           "num_kmers": 4840, "labels": [...], "dropped_labels": [], "label_population": {...},
           "labels_from_seed": false, "labels_supporting_total": 2, "labels_dropped": 0,
           "labels_dropped_digest": ""},
  "label_dict": [{"name": "NZ_STEQ01000045.1", "kind": "header"}, {"name": "573", "kind": "column"}],
  "arms": {
    "right": {
      "status": "complete", "frontier_remaining": {"live_paths": 0, "live_labels": 0, "exact": true},
      "segments": [
        {"id": 0, "parents": [], "from_bp": 0, "length_bp": 812, "sequence": "…",
         "labels": [0, 1], "labels_at_end": [0],
         "events": [{"at_bp": 300, "type": "switch", "from": 1, "to": 2, "cost": 0.7},
                    {"at_bp": 455, "type": "label_end", "label": 2, "reason": "label_lost",
                     "structural_successors": 1, "unpermitted_labels_present": 3}]},
        {"id": 1, "parents": [0], "from_bp": 812, "length_bp": 40, "labels": [], "labels_removed": [], "sequence": "…"}
      ],
      "splits": [{"at_bp": 812, "segment": 0, "children": [1, 2], "kind": "divergence"}],
      "joins": [], "tips": [], "bubbles": [], "blocked": [],
      "paths": [{"id": 0, "segments": [0, 1], "length_bp": 852, "end_reasons": {"dead_end": 1},
                 "end_labels": [{"label": 0, "loss": 0, "branches": 0}],
                 "continuation": null}]
    },
    "left": {"…": "same shape"}
  },
  "label_summary": [{"label": 0, "right": {"direct_bp": 852, "reach_bp": 852, "reentries": 0,
                                           "runs": [{"path": 0, "from_bp": 0, "to_bp": 852, "entered_by": "seed", "end_reason": "dead_end"}]},
                                  "left": {...}}],
  "diagnostics": {"…": "§7.2"}
}
```

- **A seed whose permitted set could not be DERIVED** (§6.1 step 4) takes the place of this object with
  `{"seed": {"seed_id", "length_bp", "labels_from_seed": true}, "error": "<reason>"}` and no `arms`: the request
  still succeeds and the other seeds are traversed. Presence of `error` is how a client tells the two apart.
- `segment.labels` = names live at the segment's **first** step; `labels_at_end` at its last; with
  `label_lists: delta` a child lists only `labels_removed` relative to its parent's end set. Label sets change only
  through events, sorted by `at_bp`.
- `sequence` holds the bases the segment adds, natural orientation. A right flank is the concatenation root →
  leaf; a left flank leaf → root. `"sequences": false` omits strings (sizes and events stay). Whole flanks are
  reconstructed client-side so shared prefixes are never duplicated.
- `label_summary` per label and arm: `direct_bp` = longest stay-run from the seed. It means *there is a graph
  route of that length on which every k-mer carries the label* — **not** that the label's sequence contains those
  bases contiguously (a label may hold several contigs or repeat copies; only `support: trace` or validation
  against the source record establishes a contiguous occurrence). `reach_bp` = end of the label's lineage
  including switches, plus `reentries` and all runs.
  Three distinct claims must not be conflated: **label-consistent route support** (`direct_bp`), **support for a
  particular displayed path** (that run's `route_bp == 0`, or only its `[route_bp, to_bp)` part), and **a
  contiguous source occurrence** (trace or source validation). So "N of M labels share this exact displayed flank
  to x bp" counts runs with `entered_by: seed`, `from_bp: 0` **and `route_bp: 0`** — the first two alone establish
  support along *some* route, not along the one shown.
- **Continuation.** Every leaf ended by `max_extension_bp`, a resource reason or `beam_pruned` carries
  `continuation: {sequence, labels, loss_used, branches_used}`: the last `max(k, bp since the last switch)` bases
  up to `continuation_bp`, natural orientation, with the labels covering that whole tail. It is valid `/traverse`
  input (a fresh seed: edge-reuse and branch state reset; the agent may reduce `loss_budget` by `loss_used`).
- `detail: summary` returns `seed`, `label_summary`, per-leaf `{length_bp, n_labels, end_reasons}`, `diagnostics`
  and `continuation`s, but no segments, events or sequences — the cheap probe to run before committing to a radius.
  `tree` adds segments without sequences.

### 7.2 Diagnostics (the tuning evidence, always present)

- `growth[arm]`: per bin of `profile_bin_bp`: max live paths, distinct live labels, live (path, label) pairs,
  divergences, ambiguous branches taken, splits, reconvergences, bubbles, tips, `blocked_repeat`, label ends by
  reason, steps.
- `branch_events[arm]`: the first `max_branch_events` in queue order **plus** the top `max_branch_events` by labels
  affected; each with position, successor characters, per-successor label counts and lookahead
  `{bp_until_end_or_window, reason}`, ambiguous labels, labels dropped by `branch`/`minority`/`loss_budget`;
  `branch_events_truncated` = total − shown.
- `needed_budget` histogram per arm (from `loss_budget` events); `cost_preview` for P under pairwise models
  (min/median/max of `cost(seed label → extra)`, units stated).
- `counters`: steps, successor enumerations, annotation access (path used, keys mapped, rows reconstructed,
  direct cell reads, tuple rows, pair evaluations, rc_index_range calls), resource-cap hits; and the work of
  the structural and branching rules — `edge_reuse_probes` (uses scanned or (edge, segment) pairs probed by
  the per-path edge-reuse check, whichever was fewer), `reminimisation_rounds` (re-derivations after the
  first at ambiguous nodes, §6.4; bounded by |σ| per node) and `max_reminimisation_rounds` (the largest at
  one node) — so that a pathological locus is visible rather than silent.
- `timing` (excluded from determinism): elapsed per phase, cache hits, and the **physical fetch counters**
  `rows_fetched`, `tuple_rows_fetched`, `coords_mapped` — prefetching along unbranched runs changes them with
  `annotation.batch_kmers` while the walk, `rows_requested` and `keys_mapped` do not (§6.8), so they are not part
  of the invariant `result`.

Consistency invariants (tested, T25): Σ_bins steps = counters.steps; Σ_bins label ends by reason = label_summary
ends by reason; splits of kind `ambiguous` = Σ_bins ambiguous branches taken; max live paths in any bin ≤
`max_live_paths`.

### 7.3 Provenance (response level)

`k`, `regime`, `alphabet`, graph and annotation types, annotation access path, `algorithm_version`,
`release` (configured id, §10.3), the normalized strategy, the seeds as validated, `capabilities` (§10.3, now
including `label_modes`), and `walk_rule` (§6.10).

### 7.4 The trie view and the completeness contract

Per arm, in both modes:

- `complete_to_bp` (§6.10) next to `status`, `frontier_remaining` and `cap_trigger`.
- `labels_per_node: {cap, max_seen, nodes_truncated}` — `cap` is `labels.max_labels_per_node`, `max_seen` the
  largest true label count met at a node, `nodes_truncated` how many recorded lists lost labels to the cap: the
  root's boundary list, every node entered, every successor not taken (listed on its `blocked` / `hairpin`
  event, the only place it appears) and every trie-view branch. Non-zero means the recorded sets are
  incomplete and must not be read as "these labels and no other". Every such list carries its true count:
  `segments[].labels_total` / `labels_truncated` (the entry node; for the root, the seed boundary),
  `label_sets[].labels_total`, `events[].labels_total` / `truncated` on `blocked` and `hairpin` events, and
  `branches[].labels_distinct`.
- `splits[]` is the trie (design note §5.1.1): `at_bp` (echoed as `prefix_bp`) is the shared prefix length from
  the seed boundary, and `branches: [{segment, char, labels_distinct, labels, labels_truncated}]`, parallel to
  `children`, gives per branch its first base and the labels at its first node — the true count and a list cut
  at the cap. `labels_before` is the count at the split node. **A branch's count is not a share of
  `labels_before`: one label may follow several branches** (an ambiguous branch in `constrain` mode; always
  possible in `annotate` mode), so the children's counts can sum to more than the parent's.
- `label_mode` on every seed result.

In `annotate` mode additionally, per segment, `label_sets: [{from_bp, to_bp, labels, labels_total, truncated}]`:
the labels present along the segment as maximal runs with an identical (capped list, count), half-open over
outward base indices so that `[from_bp, to_bp)` covers the nodes entered by steps `from_bp + 1 … to_bp`. The
root segment has no run (it adds no base); its `labels` list is the seed boundary's own set. `end_labels` and
`end_reasons` are empty, `path_reason` is always set, and a `continuation` (at the radius or a cap) lists the
labels recorded on *every* node of its tail — each of them validates as a seed label of that continuation, so the
tail stays valid `/traverse` input.

## 8. Annotation access (`LabelOracle`)

### 8.1 Oriented node → annotation key

| Regime | Key for oriented node `n` |
|---|---|
| `basic` | `n` |
| `primary` (`CanonicalDBG`) | `get_base_node(n)` (`canonical_dbg.hpp:38`), O(1) |
| `canonical`, DBGSSHash | `min(n, reverse_complement(n))` |
| `canonical`, DBGSuccinct | `min(n, node of rc(k-mer))`: map `rc(window)` with `map_to_nodes_sequentially`, reverse, pairwise min with the oriented ids (one O(k) index_range per chunk instead of per node) |
| `canonical`, hash/bitmap | `map_to_nodes(spelled window)` |

Never call `AnnotatedSequenceGraph::get_labels(node)` / `graph_to_anno_index` with an oriented id: on `primary`
it indexes out of range, on `canonical` it silently returns an empty row. `CanonicalDBG`'s empty
`runtime_error("")` on an inconsistent primary graph is caught and rethrown as an internal error naming the node.

### 8.2 Membership for a small permitted set P

Backend-specific, chosen once per request (`get_matrix()` is called once: it locks and flushes on
`ColumnCompressed`):

| Backend | Path |
|---|---|
| ColumnMajor (`.column`), RowFlat, BRWT, BinRelWT, CSC | `direct`: `GetEntrySupport::get(row, c)` for c ∈ P |
| `RowDiff<M>` with M ∈ {ColumnMajor, BRWT, RowFlat} (binary) | `rd_direct`: XOR of `diffs().get(r, c)` along the `row_diff_successor` path to the first anchor, decoded incrementally along unbranched runs (`row_P(v) = row_P(u) ⊕ diff_P(u)` when v is u's rd-successor) |
| Tuple matrices (`*_coord`) incl. `TupleRowDiff<TupleCSC<BRWT>>` (refseq33m) | `tuples`: batched `get_row_tuples(rows)`; header labels resolved by `CoordToHeader::map_single_coord`; column labels by column |
| RowDisk, RowDiffDisk, Rainbow, CSR, others | `rows`: batched `get_rows(rows)` ∩ P |

`access: auto` picks by this table; `direct`/`rows` may be forced (a 400 if unsupported — no silent fallback).
Header labels always use `tuples`. Never `get_column`, never `AnnotationBuffer`. A per-seed cache
`annotation key → interned hit set` bounded by `max_steps`. T1 asserts `direct = rd_direct = rows = tuples`
agreement where several apply.

### 8.3 Batching along unbranched runs

Walk structurally up to `batch_kmers` nodes while each node has exactly one non-`'$'` successor; map and fetch
labels for the chunk in one call; then apply the recurrence node by node and truncate at the first node where
the path ends (also at a cap, §6.8). At a structural branch, fetch all successors in one call. Unverified nodes
are never emitted. Results do not depend on `batch_kmers` (T23).

### 8.4 Graph handle, caches, threading

All reads are const and concurrent requests share the loaded `AnnotatedDBG` (as `/search` does). The request owns
a `TraversalContext` holding: the walk graph (for `primary`: a lightweight `CanonicalDBG` clone over the same
primary graph with a `NodeFirstCache` attached — the shared wrapper is never mutated), a `NodeFirstCache` for any
underlying `DBGSuccinct`, the matrix cross-casts, the `CoordToHeader` pointer, the label cache and counters.
Per-step costs to expect: `basic` DBGSuccinct right step O(σ) rank/select, left step 2 bwd with the cache;
`primary` first-visit step adds O(k) (`rc_index_range`, counted).

## 9. Label-change cost models

| Model | `cost(ℓ → ℓ')` for ℓ ≠ ℓ' |
|---|---|
| `forbid` | +∞ |
| `constant` | `value` |
| `table` | listed entries; others `default` |
| `hll` (I8) | `−scale · log2(min(D, Î) / D)` with `Î = |ℓ| + |ℓ'| − |ℓ ∪ ℓ'|_est` from per-column HyperLogLog sketches; +∞ if `Î ≤ 0`; clamped to ≥ 0. `denominator` D ∈ {`from`, `to`, `min`} |

The `hll` model ports the label-change score from `origin/hm/aln_alt_label_change`
(`AnnotationBuffer::get_label_change_score`: `log2(min(|B|, |A|+|B|−|A∪B|)) − log2|B|`). The port requires the
`libcount` submodule (or vendoring it), the CMake `hll` target, `hll_counter.{hpp,cpp}`, `hll_matrix.hpp`, and a
`transform_anno` hook writing a `.hll` sidecar that **also stores the label list** (refused on mismatch with the
loaded encoder). `LabelChangeCost` owns the sketches; nothing is attached to the graph. `hll` is a 400 when no
sidecar is loaded. Pair costs are cached per request; the argument order of the branch's chainer call
(`a_i.col, a_last->col`) is checked before fixing the default `denominator`. Header labels have no sketch; `hll`
applies to column labels only.

## 10. Interfaces

### 10.1 C++ (core library, `src/graph/traversal/`, no Config/JSON/HTTP)

```cpp
namespace mtg::graph::traversal {
class TraversalContext {            // §8.4: per-request graph handle, caches, oracle, counters
  public:
    TraversalContext(const AnnotatedDBG &, const OracleOptions &);
    const DeBruijnGraph& walk_graph() const;  Regime regime() const;
    LabelOracle& oracle();  Counters& counters();
};
SupportProfile resolve_support(TraversalContext &, std::string_view query, const ResolveOptions &);
SeedSelection  select_seeds(const SupportProfile &, std::string_view query, const SelectionPolicy &);
struct Seed { std::string seed_id, sequence; std::vector<std::string> labels; };
struct Strategy;                     // §5, validated and normalized
class LabelChangeCost;               // §9
SeedResult traverse_seed(TraversalContext &, const Seed &, const Strategy &, const LabelChangeCost &);
}
```

Helpers to add: `map_to_rows` declared in `annotated_dbg.hpp`; `encode_presence_mask` moved to a header; a label
lookup that throws `InvalidRequest` naming the label. JSON parsing/validation and serialization live in
`src/cli/traverse.{hpp,cpp}`, shared by CLI and server.

### 10.2 CLI

```text
metagraph traverse  -i GRAPH -a ANNOTATION [-p THREADS] request.json [...]
metagraph traverse --resolve -i GRAPH -a ANNOTATION request.json [...]
```

Request files are exactly the JSON bodies of §4.1/§5; one JSON document per request on stdout. One code path with
the server.

### 10.3 Server

`POST /resolve`, `POST /traverse`, `GET /traverse/capabilities` in `server.cpp`.

- **Capabilities** (also in `/stats` and every resolve/traverse response): `k`, `regime`, `alphabet`,
  `num_labels`, `has_coordinates`, `has_coord_to_header`, `cost_models_available`, access path, server maxima,
  `schema_version`, `release`.
- Single-graph mode: as `/search`. **Multi mode:** `graph` (name) is required; the server selects the
  `(graph, annotation)` pairs under that name whose labels contain every seed and `extra` label: exactly one →
  use it; several sharing one loaded graph → a union oracle (membership = OR, column ids namespaced per pair);
  otherwise 400 listing candidate pairs and how many requested labels each holds. `graph_path` is an optional
  override. The exact pair(s) are echoed. `in_ram` deployments are rejected. The 4th CSV column (`release`) is
  parsed explicitly (today silently dropped) and `--index-release ID` serves single mode.
- Caps: `--traverse-max-time-ms` (default 30 000), `--traverse-max-seeds` (64), `--traverse-max-seed-bp`
  (100 000), `--traverse-max-seed-labels` (10 000) and `--resolve-max-query-bp` (0 = unlimited). `0` means
  unlimited for each. A cap that lowers a request bound is **echoed** as `clamped` (§5); `max_seeds` and
  `max_seed_bp` are refusals (400), not clamps. The traverse caps are non-zero by DEFAULT because a request may
  name no labels at all: the server then chooses the permitted set, so the request's own bounds no longer describe
  the work it asks for, and an unconfigured deployment must not be the unlimited one. `--traverse-max-time-ms` is
  also the only bound on the derivation phase, so when a request has a derived seed a `bounds.time_budget_ms` of 0
  (which for the walk means "do not extend") is clamped UP to it rather than left unbounded — also echoed as
  `clamped`. The CLI (`metagraph traverse`), being operator-run, applies no caps. A request pinned to another
  `release` is a 400 naming both ids.
- Errors: a new `InvalidRequest` exception → 400 with the message; every other exception → 500 with context
  (small change to `process_request`). 503 while the index loads. Handlers run on the io threads: `-p ≥ 2` is
  documented as required for serving `/traverse` beside `/search`.

## 11. Tests

### 11.1 Build configurations and gates

Gate for every increment: `./unit_tests --gtest_filter="*Traversal*:*LabelOracle*:*Resolve*"` passes in `build/`
(Release) **and** in `build_debug/` (`cmake -DCMAKE_BUILD_TYPE=Debug`, asserts on). At I4 and I7 the concurrency
test T24 runs in `build_tsan/` (`-DCMAKE_BUILD_TYPE=Threads`). One ASan run before I7.

### 11.2 Fixture tiers

- **Tier O** (oracle tests T1–T3): {DBGSuccinct, DBGHashFast, DBGHashOrdered, DBGBitmap, DBGSSHash} ×
  {ColumnCompressed, RowFlat} plus DBGSuccinct × {RowDiffColumn, RowDiffDisk, BRWT, RowDiff<BRWT>,
  ColumnCoord, RowDiffCoord, RowDiffBRWTCoord}, modes basic/canonical/primary, with the expected regime per cell
  (SSHash "primary" ⇒ `canonical`; DBGHashString excluded). New instantiations added to
  `tests/annotation/test_annotated_dbg_helpers.cpp`.
- **Tier W** (walker semantics T5–T16, T21–T26): {DBGSuccinct, DBGHashFast, DBGSSHash} × ColumnCompressed × 3 modes,
  plus DBGSuccinct-primary × RowDiff<BRWT>, plus DBGSuccinct-basic × RowDiffBRWTCoord with a CoordToHeader
  (the refseq33m shape). All fixtures use odd k ≥ 11 and call `expect_topology()` before asserting results.
- **Tier U** (unmasked): a `build_anno_graph` overload with `mask_dummy_kmers = false` (ported from
  `origin/hm/aln_alt_label_change`), DBGSuccinct basic and primary × {ColumnCompressed, RowDiff<BRWT>}.
- **Tier R** (real format): a small index in exactly the `refseq33m` format, built from 42 real blaNDM-carrying
  RefSeq records by `scripts/traversal/build_mini_refseq.sh` (~9.3 Mbp, ~25 s, output `build/mini_refseq/`).
  `tests/graph/traversal/test_mini_refseq.cpp` (gtest filter `MiniRefSeq*`) runs against it and **skips when it is
  absent**, so it never blocks a normal unit-test run; `integration_tests/test_traverse.py` covers the CLI and the
  server routes on an equivalent fixture it builds itself. Later: one ATB shard for measurements.
  RowDiff anchor selection is not byte-deterministic across rebuilds, so these tests compare query results, never
  index bytes.

### 11.3 Test list

| ID | Case | Expectation | Design row |
|---|---|---|---|
| T1 | oracle agreement | `direct` = `rd_direct` = `rows` = `tuples` (where applicable) on random rows × label sets, all Tier O cells | Annotation backend |
| T2 | key mapping | primary/canonical: both orientations map to one key and see the label; basic: rc lookup gives npos or a row without the label | Reverse complements |
| T3 | resolve masks | label-present bits equal `get_top_label_signatures(seq, SIZE_MAX, 0, 0)`; T3b: hand-built query (200 bp, 15 bp foreign insert at 80, A on [0,120), B on [100,200)) gives exact runs and the three states; header labels and `trace` runs on the Coord fixture | Inside-query support |
| T4 | candidates/select | disconnected blocks separate; gapped seed rejected; T4b literal expected seeds for each policy incl. `max_support` and `merge_overlapping`; T4c discovery truncation {kept 2, total 3}; T4d RLE round-trip; `label_order: hash` stable under column permutation | Partial seed |
| T5 | linear | T5a right arm, T5b both arms: left + seed + right reconstructs the source exactly; seed once; bp counted outside the seed | Linear sequence |
| T6 | label end | T6b exact coordinates: k=11, A on W·R·X, B on R·X·Y (R 60 bp, X 25, Y 30, W 20): right arm A `label_lost` at 25, B `dead_end` at 55; left arm B `label_lost` at 0, A `dead_end` at 20; every run checked against a per-base brute-force scan | Support changes / Sequence interval labels |
| T7 | divergence | shared prefix then split; 2 leaves; no branch counted | — |
| T8 | ambiguous | A on X·P·Y and X·Q·Z, seed X: limit 0 → one leaf, `branch` at 0 with successor chars {P[0], Q[0]}; limit 1 → two leaves at |P|+|Y| and |Q|+|Z|, branches 1; exact `growth` for bin 10 | — |
| T9 | switching chain {A,B}→{A}→{A,B}→{B} | pinned events: `forbid`: A `label_lost` at end(b3), B `label_lost` at end(b1), one leaf; `constant 1`/budget 1: switch A→B at from(b3), A `switched`, leaf reaches end(b4) with B loss 1; cost 0: same path, loss 0; `switch_on: loss` vs `any` differ only when both continue | Sample-switching chain |
| T10 | diamond | A on X·P·Y·Q, B on X·P'·Y·Q', seed X {A,B}, forbid, limit 0: one divergence, 2 leaves, P-leaf ends at Q with {A} only, Q' never reached from P; with `on_reconverge: merge` the Y segment is shared with parents [P, P'] and labels {A,B} | Diamond |
| T11 | cycles | T11a simple cycle L·J·C·J·E: limit 0 → `branch` at end(L); limit 1 → exactly leaves L→J→E and L→J→C→J→E; T11b homopolymer self-loop taken once + one `blocked` event; T11c figure-eight limit 2 → enumerated leaf set; T11d circle of C bp, seed s ∈ {k, k+20}: each arm one leaf `rejoined_seed` at C − s + k − 1 with `overlap_bp` k−1 | Cycle / self-loop |
| T12 | hairpins | linear fixture containing an RC-palindromic (k−1)-mer in canonical/primary: reconstructs with no branch counted, one `hairpin` event; `follow` flags the hairpin leaf and blocks the RC retrace with `edge_reuse_rc` | Reverse complements / palindromes |
| T13 | RC seed | for `status: complete` on both arms, rc(seed) yields the same multiset of (rc(flank), label runs, end reasons) with arms swapped; ids/order may differ | Reverse complements |
| T14 | dummies (Tier U) | `'$'` never counts as a branch; tips → `dead_end`; identical results masked vs unmasked | — |
| T15 | caps | `max_extension_bp` → `status: complete`; each of max_steps/max_live_paths/max_paths/max_output_bp → that reason on every live path of its scope, `truncated`, `cap_trigger` set, other arm unaffected for per-arm caps; `time_budget_ms: 0` stops after the first step; T15b `lowest_loss_first` + `max_live_paths: 1` on T9 keeps the loss-0 path; beam survivors listed for T8 limit 3 | Caps, cancellation |
| T16 | determinism | `result` byte-identical across 3 runs, seed permutations, and one-seed-per-request splits | Batch reorder |
| T17 | validation | unknown keys, bad labels, short/invalid seeds, `extra` under forbid, `trace` without coords, `access: direct` on RowDisk → precise errors; T17b `seed_unsupported` drop | Stale anchor / missing label |
| T18 | JSON golden + CLI | golden files in `tests/data/traverse/*.json`, regenerated only with `--update-golden`; CLI end-to-end on the mini-refseq index | — |
| T19 | server | `integration_tests/test_traverse.py`: `/resolve`, `/traverse`, capabilities; clamp echoed; release mismatch 400; multi-mode label routing and ambiguity 400; 503 during load | — |
| T20 | hll | T20a cost equals §9 formula on sketch estimates, +∞ at Î ≤ 0, clamp ≥ 0; T20b on 6 columns with 10–90 % overlaps of 10⁵ rows, relative error ≤ 3·1.04/√m; label-list mismatch refused | — |
| T21 | brute-force loss oracle | 200 random graphs (k=11, 3–5 labels, fixed RNG): reported min loss per leaf = brute-force min over assignments, constant and table costs, `max_label_branches = ∞` | — |
| T22 | reference walker | naive per-label DFS implementing §6.4–6.6 for one label; under forbid/limit 0 the (label, end at_bp, reason) set equals the union over labels | — |
| T23 | invariance | `result` identical for access ∈ {direct, rd_direct, rows, tuples} where supported, `batch_kmers ∈ {1, 2, 7, 64}`, `sequences ∈ {true, false}` (minus strings) | — |
| T24 | concurrency (TSan) | 8 threads × 50 seeds on shared primary CanonicalDBG(DBGSuccinct)+ColumnCompressed and DBGSuccinct+RowDiff<BRWT>: equal to sequential; shared wrapper's extension list unchanged | — |
| T25 | diagnostics invariants | §7.2 invariants on every Tier W fixture | Aggregate drill-down |
| T26 | repeat / two molecules | A on X·R·Y and Z·R·W, seed R, limit 1 → 2 leaves both A; output states support kind, never "one molecule" | Same-label repeated motif |
| T27 | G4 population fixture | 40 synthetic 5 kb genomes sharing a 2 kb element E; 10 share right flank F1, 10 share F2, rest unique. forbid/limit 0 + merge: leaves ≤ 40, F1/F2 cohorts reach 1 kb, every label ends `branch`/`dead_end`, growth shows the junction; cost 0 + `max_live_paths 8` → truncated with `cap_trigger` at the junction; `min_successor_labels 5` → only F1/F2 followed, others `minority` | — |
| T28 | tips and bubbles | SNP bubble inside one label: limit 0 stops with `branch`; `bubble_window_bp 100` passes with a `bubble` event and no branch; 5 bp tip with `tip_window_bp 20` → `tip` event | Bubble phase |
| T29 | trace support | Coord fixture with two occurrences of R in one accession: `kmer` support follows both continuations as one label; `trace` keeps them apart and ends with `trace_break` at a coordinate jump; a cross-record boundary stops with `trace_break` even with consecutive global coordinates | Cross-record boundary |
| T30 | derived seed labels | Seed carried in full by 2 of 3 labels, `labels` omitted: derives exactly those 2, `labels_from_seed`, results identical to naming them; `max_seed_labels` truncation reports counts + digest and does NOT drop them; empty intersection rejected; `extra` on top of a derived set; coordinate fixture: header kind derives accessions, column kind the column, `trace` holds the derived set to coordinate continuity; mini-refseq: the whole blaNDM gene with no labels derives the 19 full-length carriers and traverses identically to the explicit run | Hit-labelled seed |
| T30b | derived-set cost and usability | 25 labels, one carrying only a prefix: the derived set and the same list named explicitly agree field by field, and `labels_supporting_total` equals the kept count for the explicit list too; `bounds.time_budget_ms` stops the derivation (which runs before the walk's clock check) while leaving an explicit request a truncated walk; 100 decoy records sharing only the seed's first k-mer: the intersection starts from the cheapest row of the first batch, so `coords_mapped` stays at two per k-mer instead of paying for that row; the same accession in two columns makes a derived `header` set ambiguous and is refused naming the header and `seed_label_kind`, while the `column` kind and an explicit `ACC1` still work; under `trace`, a set that halves outside the first batch (8 carriers of the first 100 k-mers, 2 of the rest) still reports the coordinates the explicit run finds, i.e. the retained per-k-mer coordinates survive being compacted; mini-refseq: derived and explicit agree on dropped labels, label summary, runs, access path and the row/key counters | Hit-labelled seed |
| T30c | per-seed derivation failure | Three seeds, the middle one carried in full by no label: HTTP 200 / exit 0 with `{"seed": {...}, "error": …}` in its place and the other two traversed; the same seed with an EXPLICIT label list still fails the whole request; `max_seed_labels` out of `[1, 100000]` is a parse error; a seed over `--traverse-max-seed-bp` and more than `--traverse-max-seeds` seeds are 400s; `max_seed_labels` above `--traverse-max-seed-labels` is clamped and reported | Hit-labelled seed |
| T31 | trie oracle (§6.9) | `test_trie.cpp`: annotate records what constrain filters; the exhaustive preset refuses conflicting knobs; a tripped cap reports `complete_to_bp` and the partial level is excluded; cut recorded lists are reported and break the oracle; `TrieOracle` (4 graph × annotation pairs, 3 modes): `claims(A) == E` for 5 permitted sets on both arms, tuned runs are prefix-subsets with a reason for every omission, merged routes are sound | — |
| T32 | records model, edge and non-edge cases | `test_trie_cases.cpp` + `test_trie_reference.hpp` + `test_trie_checks.hpp`, two graph × annotation pairs, three modes, both arms, every case three ways (records ⇔ T ⇔ A with end reasons) plus single-label union, derived == explicit and the recurrence at budget 0: linear; seed at record start / end / the whole record; radius 0, 1, 2, |R|−1, |R|, |R|+1 (a head at the radius is `max_extension_bp`, not `dead_end`); fork; unequal bubble (+ merged routes); nested bubbles; one-base tip at the seed boundary; homopolymer self-loop (3 walks, never the record's); cycle junction; circle (`rejoined_seed` on both arms); repeat in two records plus a chimera; seed twice in one record (the walk runs through the k−1 junction k-mers before `rejoined_seed`); seed spanning a bubble (the other label is dropped, nothing else changes); an unlabeled region (recorded empty, `label_lost` at its edge); a reverse-complement record (not a carrier in basic, the second fork branch in canonical/primary); hairpin skip and follow; even-k palindromic node; 100 labels (default cap cuts the lists and the oracle says so; uncapped the contract holds); unmasked DBGSuccinct ('$' never recorded); 16 random fixtures at k = 7 | Cycle / self-loop, Reverse complements |
| T33 | switching vs the §6.3 recurrence | `OneSwitch`: forbid ends A in Y; budget 1 spells P·Y·Z under B at loss 1 switched at |P·Y|, budget 0 reports `loss_budget` needed 1; `TwoSwitches`: budget 1 cuts A's walk at |P·Y·Z| needing 2, budget 2 spells P·Y·Z·W under C at loss 2, B's walk needs one; leaves, losses and switch events equal the recurrence over the records | Sample-switching chain |
| T34 | trace vs the records | `TrieCasesTrace`: k-mer support follows the four combinations of two records sharing a stretch, trace the two records; a record ending where another goes on is `record_end` under trace and continues under k-mer; a seed twice in one record gives two traces; the homopolymer limitation (§6.6) | Cross-record boundary |
| T35 | size caps vs the completeness guarantee | `TrieCasesCaps`: `max_steps` 1…120, `max_live_paths` 1–3 (stop and beam), `max_output_bp` 1…90, `max_paths` 1–2 on the nested bubbles, three modes: every tuned run is a prefix-subset with a reason per omission, its walks equal the exhaustive trie's at every depth ≤ `complete_to_bp`, and `complete_to_bp` never decreases as `max_steps` grows | Caps |
| T36 | the contract on a real index | `MiniRefSeq.TrieContractAgainstTheSourceRecords`: the string model from the 42 source records (k = 31, basic, unmasked, RowDiff<BRWT> + coordinates, header labels); the whole blaNDM seed (19 carriers; 24 structural walks left, 2 right at 300 bp) and three 150 bp windows: records ⇔ T ⇔ A, trace ⇔ records, no '$' in any output; at 1000 bp the 19 carriers' claims under k-mer and trace support equal the records' | — |

## 12. Implementation increments (each with tests, then an adversarial review)

| # | Deliverable | Tests | Est. |
|---|---|---|---|
| I1 | `TraversalContext`/`LabelOracle`: regime detection, key mapping per regime, `direct`/`rd_direct`/`rows`/`tuples`, header labels via CoordToHeader, cache, counters; Tier O/U fixture instantiations | T1, T2 | 4–5 d |
| I2 | `resolve_support` + `select_seeds` (incl. `trace` runs, discovery, policies, seed_id) | T3, T4 | 2–3 d |
| I3a | Right-arm walker, `forbid`, per-lineage branching, edge reuse (ancestor-indexed map), `rejoined_seed`, hairpin skip, stop reasons, segment DAG, events recorded at decision time, bounds | T5a, T7, T8, T10 (no merge), T11, T12, T14, T15 (non-beam), T22, T25 | 1.5 wk |
| I3b | Chunked walking with hints, `NodeFirstCache`, left arm, seed validation, reconvergence merge, tip/bubble windows, quorum | T5b, T6, T10 (merge), T13, T23, T26, T27, T28 | 1 wk |
| I3c | Thin `metagraph traverse` CLI (fields implemented so far) + `benchmarks/graph/bench_traverse.cpp` (access paths × |P| × backends) + mini-refseq build script | T18 (partial) | 3 d |
| I4 | Cost models `constant`/`table`, switching with `switch_on`, runs, `lowest_loss_first`/`most_supported_first`, beam, `trace` support, continuation | T9, T15b, T16, T21, T24, T29 | 1 wk |
| I5 | Growth profile, branch-event selection, `needed_budget`, `cost_preview`, `detail` modes, delivery bounds | T25 extended | 3 d |
| I6 | Strict JSON validator/normalizer, full CLI, golden files | T17, T18 | 3–4 d |
| I3d | Measurement on the mini-refseq index and (when available) one ATB shard: latency, RSS, rows per step, leaves; fix default bounds, `batch_kmers`, auto thresholds; recorded in `docs/` | — | 2 d |
| I8 | `hll` cost: libcount, sketches with label list, build hook | T20 | 1–1.5 wk |
| I7 | Server routes, capabilities, multi-mode label routing, release id, error classes | T19 | 4–5 d |
| I9 | Compacted explored subgraph + GFA export | (later) | 1 wk |

## 13. Deviations from the design note and open decisions

**Deviations (deliberate):**

1. `within_label` / `path_common` / `any_of` / `sequence_trace` are special cases of one mechanism (label set +
   support kind + change cost + loss budget), §6.3.
2. Branching is counted per **lineage** along a path; divergences of label sets are not branches (§6.4).
3. Edge reuse is defined on canonical (k+1)-mers in canonical regimes (§6.6); arms do not share used-edge sets.
4. Loss budget and branch limits are per arm; arms are not joined into whole paths.
5. Reconvergent paths are merged into a DAG by default (§6.5); the design's "unique-node exploration records
   cycles as edges" is realized by `blocked` events rather than by a global visited set.
6. The primary output is the extension DAG; the compacted subgraph/GFA comes later (I9).
7. Permissible merges (in-degree > 1) are not detected in v1.
8. **Delivery bounds are not implemented.** `max_output_bytes` and `max_events` and the
   `delivery: {truncated, next_cursor}` block of §6.7 do not exist; output is bounded in bases
   (`max_output_bp`) and events (`max_branch_events`) only. The request schema rejects the two unimplemented
   fields rather than accepting and ignoring them.
9. **`support: trace` cannot be combined with `on_reconverge: merge`** (rejected with an actionable error).
   A merge unites both parents' coordinate sets but keeps one parent's `(loss, branches, run, route)`, so a later
   step can continue on the discarded occurrence's coordinates while carrying the retained one's evidence —
   fabricating direct support, a loss or a run that no single occurrence justifies. Supporting it needs
   per-coordinate provenance, which is deferred.
10. **`annotation.access` is not client-selectable.** The access path follows the annotation representation and
   the label kinds (header labels require the coordinate path); the C++ API can override it, the JSON API
   rejects anything but `"auto"`. The `rd_direct` accessor of §8.2 is **not implemented** — `RowDiff` membership
   goes through full row (or tuple) reconstruction, which is the main open scaling question.
11. ~~**`seeds[].labels` is a required explicit list**, so the design note’s `labels: {"permit": "per_hit"}`
   default (design note line 183, line 331: "`per_hit` explicitly means the resolved hit’s annotation label, not
   all graph labels") had no counterpart: every request had to name the labels, which is not what a hit’s own
   label means.~~ — **resolved.** The list is optional, and omitting it derives the permitted set from the seed:
   the labels supporting every seed k-mer, which *is* `per_hit` (§5, §6.1). The derivation is the intersection of
   the rows the seed validation already reads, so it adds no annotation reads and never looks past the seed;
   `max_seed_labels` bounds the traversal state, not the discovery, and what it cuts is counted and digested.
   Two consequences are deliberate rather than accidental: (a) a `change_cost` `table` is authored by label *name*
   over the request order (§9, note 7 below), which a derived seed does not have at request time, so the
   combination is rejected naming the field instead of being mapped onto the wrong indices; (b) an `extra` label
   the derivation also produced is rejected as a duplicate seed label — the caller named a switch target that is
   already permitted, which under a derived set is worth an error rather than a silent merge.

**Implementation notes of increment I3 (`src/graph/traversal/walker.cpp`), deviations to resolve or confirm:**

1. **Tip and bubble windows are not implemented** (`tip_window_bp`, `bubble_window_bp`, the `tip` / `bubble`
   events and the `tips` / `bubbles` growth counters). The fields exist and stay zero; T28 is deferred.
2. `EndReason` has no `switched`, `superseded`, `minority`, `below_min_labels`, `split_limit`, `switch_sources`
   or `hairpin` value (the enum is the contract of `traversal_types.hpp`). A source whose lineage continues only
   under another name ends its run with `label_lost` and the continuation is documented by the `switch` event;
   `superseded` is `label_lost` with `text: "superseded"`; quorum and split-limit stops are `branch` with
   `text: "minority"`, `"below_min_labels"` or `"split_limit"`; a label whose only continuation is a skipped
   hairpin ends with `dead_end`, `text: "hairpin"`; a label cut from the switch sources by `max_switch_sources`
   (pairwise models) whose cheapest switch would have been finite ends with `label_lost`, `text:
   "switch_sources"`.
3. `loss_budget` is reported only when the switch target is **not** already carried by a cheaper predecessor at
   that successor (otherwise the label would have been superseded anyway and raising the budget would not
   continue its lineage); `needed_budget` is the cheapest such switch. Only the sources the recurrence actually
   considered count (the `max_switch_sources` cheapest by `(loss, branches, column, seq_id)`), and `superseded`
   is reported only for a switch that was within budget: an over-budget switch into a target another source
   carries is a plain `label_lost`.
4. Exploration is level-synchronous (one step per level for every live path of an arm, arms alternating by
   extension depth) rather than a single priority queue; `order` decides the processing order within a level.
   `max_steps` is decided per path head (a split that would exceed it ends that head), so a split never overshoots
   the cap but the last level may stop one head early.
5. `revisit` is reported once per stretch of revisited nodes, not once per node.
6. Runs closed by a reconvergence merge carry `ended: false` with `to_bp` set (the lineage continues in the kept
   entry's run); `reentries` counts `switch` events into the label. Runs are per path: at a split the first child
   containing the label keeps the run, other children get a clone with the same `from_bp`.
7. `label_dict` excludes dropped seed labels (they are neither permitted nor switch targets). The indices of a
   `table` cost model refer to the **request** (seed labels in the given order, then `extra`), not to the
   compacted dictionary: a caller authors the table before validation can drop a label, so dropped labels keep
   their index and entries naming them are ignored; the walker remaps the table to dictionary ids once. With a
   **derived** seed (§6.1) the derived labels occupy the seed-label positions of that same request order; since
   they are not knowable when the request is written, the JSON API rejects a `table` together with a derived seed
   rather than resolve its names against the wrong indices (§13 deviation 11).
8. `trace` support with **column** labels follows consecutive column coordinates and does not detect a record
   boundary with consecutive global coordinates (T29's cross-record case needs the `CoordToHeader` lookup).
   `dropped_labels.runs` are k-mer presence runs, not trace runs.
9. The graph is node-centric: a homopolymer `A^(k+1)` already branches at the node before `A^k` (`A^k` and
   `A^(k-1)·x` are both successors), so unrolling the self-loop once needs `max_label_branches: 2` (T11b).
10. `continuation.labels` falls back to all live labels when no run covers the whole tail (only reachable when
    the tail is the single current k-mer, which every live label supports).
11. **Hairpins, `follow` mode:** the self-RC step is followed, flagged with a `hairpin` event (`text:
    "followed"`) on the parent segment and never counts toward ambiguity, so the split it causes is a
    divergence; the RC retrace behind it is then blocked as `edge_reuse_rc` (or `rejoined_seed`).
12. **Reconvergence merge and routes:** a merged segment lists per parent the labels whose kept entry entered
    through it (`Segment::labels_via_parent`, parallel to `parents`), which makes per-label routes
    reconstructible although `paths` are spelled through first parents. An entry taken from a later parent
    carries the merge position as a route boundary: a `continuation` tail is cut at the last switch **or route
    join**, so its labels always cover the spelled tail. Under `support: trace` the live coordinates of a label
    present on both parents are united (the min `(loss, branches)` entry decides loss, branches and run).
13. **Edge identity above 64 bits** uses a 128-bit FNV-1a pair without keeping the (k+1)-mer: two distinct
    (k+1)-mers colliding on the pair are treated as the same edge (conservative: blocks), the spec's
    "full-string comparison on collision" is not implemented.
14. **Counters vs `batch_kmers`:** `successor_enumerations`, `keys_mapped` and `rows_requested` count work the
    walk consumes and are independent of `batch_kmers` (lookahead enumerations and keys are charged when the
    walk reaches the node, prefetch warms the cache without counting requests). `rows_fetched`,
    `tuple_rows_fetched` and `cache_hits` include lookahead work and belong to `timing`.
15. `seed_id` / `validated_seed_id` are computed over the case-mapped sequence (`make_seed_id` applies the
    build's case mapping), so `select_seeds` and `traverse_seed` agree for a seed submitted in either case.
16. **Annotate mode reads full rows** (§6.9): `LabelRecorder` reconstructs one whole row (or tuple row) per
    node, cached per seed like the permitted-set query. For `header` labels the true count needs the sequences
    a k-mer occurs in, which only the coordinate mapping tells; a sequence's coordinates are contiguous, so the
    recorder maps ONE coordinate per (column, sequence) and skips the rest of that sequence's range
    (`coords_mapped` is O(distinct sequences · log coordinates) per row, not O(coordinates)). The row cache is
    bounded in rows (1M) and in kept keys (64M, about 1 GB), whichever trips first. Accepted for a verification
    tool at small radius; the per-node cap is mandatory and bounds the output, not the read. In that mode
    `max_label_branches` and `max_splits_per_path` have no meaning (no lineage, no path-level split budget: a
    split limit would be a pruning the oracle must not have) and are `"unlimited"`; a number is refused.
    Annotate mode's `validated_seed_id` is the seed's id under an empty label list.
17. **The completeness boundary** (§6.10) is `ArmState::boundary`: the extension depth of the first head a cap
    left unexpanded (or the first head a beam pruned), minimised over events; `complete_to_bp` is its minimum with
    the radius. Two existing behaviours changed to keep `complete ⇔ complete_to_bp == radius` exact: a frontier
    head that has already reached the radius when the time budget or `max_steps` (other arm) trips is ended with
    `max_extension_bp` instead of being censored with the cap reason, and a beam does not prune heads sitting at
    the radius. `algorithm_version` is `traverse-0.2`.
18. **Seed validation** sweeps the k-mers once and scatters every hit into its label's run list (O(hits), as
    `resolve.cpp` does) instead of looking each label up at each k-mer; under `trace` the live coordinate sets
    of all labels advance in that same sweep. Nothing L × M is materialised. The one quadratic term left in
    the walk is `derive()` under a `table` cost: every (source, target) pair of a successor is evaluated,
    O(|σ| · |A(v)|) per successor, bounded by `max_switch_sources` except under `exhaustive`, where that bound
    is refused. It is avoidable for a sparse table (the best source for a target is the best-two-by-loss
    source unless an explicit entry beats it, so O(|A(v)| + entries) suffices); not implemented.
19. **The verification contract is label-granular** (§6.9): the test-side oracle compares (walk, label) CLAIMS
    — for every permitted label its maximal label-consistent walks — with the walker's label ends, not leaves
    with leaves: a label whose walk ends inside a walk that other labels continue is a claim of its own and
    its end position is checked.

**Open decisions (need input or data):**

1. Default bounds, `batch_kmers`, tip/bubble windows — fixed after I3d measurements.
2. `hll` denominator and scale (verify the chainer's argument order on the branch).
3. ~~Whether `discover` is needed at Logan-scale label counts~~ — **answered: `discover` does not scale and is
   largely redundant.** `max_labels` truncates the *output*, not the work: header discovery reconstructs every full
   row and performs a rank/select per coordinate of every label at every k-mer, so the cost scales with
   (query k-mers × labels per k-mer × coords per label). On a 33 M-label index a conserved k-mer can be carried by
   10⁴–10⁶ accessions. Worse, `/search` already returns a per-label presence mask over k-mer positions (and with a
   `.seqs` sidecar those labels are accessions), filtered server-side by top-N and `discovery_fraction`.
   **The supported path is therefore: get candidate labels from `/search`, pass them to `resolve` explicitly.**
   What `resolve` adds over `/search` is: exactly the labels you name (no top-N distortion, zeros reported), the
   distinction between "k-mer absent from the graph" and "present but unlabelled", coordinate-consecutive
   (`trace`) runs, and a frozen reproducible seed. `discover` should be demoted to a small-index convenience or
   reimplemented on top of the search path. Note that **deriving a seed’s labels (§6.1) is not `discover`** and
   does not inherit its scaling problem: it reads only the seed’s own rows — the ones the seed validation reads
   anyway — it intersects instead of ranking (so the candidate set only shrinks, and after the first k-mer only
   still-live columns are touched), and it stops at the first k-mer with an empty intersection.
4. Which ATB shard hosts the first non-RefSeq measurement.
