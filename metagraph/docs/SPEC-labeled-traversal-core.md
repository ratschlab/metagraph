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
    "frontier": {"order": "breadth_first", "on_overflow": "stop"},
    "output": {"detail": "full", "sequences": true,
               "profile_bin_bp": 100, "max_branch_events": 100, "continuation_bp": 1000, "timing": true},
    "annotation": {"access": "auto", "batch_kmers": 64}
  }
}
```

- **Strict validation.** Unknown keys, wrong types, out-of-range values and unsupported `schema_version` are
  rejected with a message naming the field; nothing is silently weakened. Defaults are filled in and the
  **normalized strategy** is echoed. If the server clamps a bound (§10.3) the echo carries
  `clamped: [{field, requested, effective}]`, each value in its knob's type (an integer for an integer knob such
  as `labels.max_seed_labels`, a number for `bounds.time_budget_ms`).
- **Request budgets** (`DESIGN-traverse-graphlet.md` §14, stage 2 of §14.1): `bounds.max_memory_mb` (integer,
  1 … 1 048 576) and `bounds.max_work_units` (integer ≥ 1), both optional — omitted means no budget, and an
  omitted budget is not echoed, so the response to a request without one is what it always was. Memory is the
  modelled bytes the request retains at once: the walker's state, the label caches (a fixed allotment of the
  budget, which they evict within) and the output **in the requested detail** (a graphlet costs the least,
  `detail: full` spells every path's chain). Work is charged work units (§6.8). Either budget stops the whole
  seed when a head cannot be admitted (§6.8), with `resource_limit`, a `resource_stop` and a `walk_domain`
  limitation (§7.0). Both budgets are **per seed** (the design's *locus* scope): each seed of a request is
  walked under budgets of its own, so that its result does not depend on the other seeds, and a request of n
  seeds can hold up to n times the memory budget (each seed's output is kept until the response is written;
  the server bounds n). A budget that does not hold the seed itself — the work of reading the seed (its
  validation, or the derivation of its permitted set), or the memory of its depth-0 state (the label
  dictionary and both roots, with what ending and delivering them costs in the requested detail) — fails
  that seed (`outcome.walks: failed`, §7.0) before anything per label is delivered; the other seeds are
  traversed. **What bounds the output**: work and the deadline bound the walk, the same in every detail (so a
  graphlet rebuilds the `detail: full` response of the same walk); the size of the output, and so the time to
  serialise it, is bounded by `bounds.max_memory_mb` alone, which charges each expansion's delivery in the
  requested detail (`detail: full` spells every leaf's chain: on a comb-shaped trie the output is quadratic
  in the walk, and only the memory budget bounds it).
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
  that successor the walk into it would be pruned with no path-level reason. Outside the preset every cut that
  may have changed a step is counted and stated (`counters.switch_sources_cut`, the `switch_sources` limitation,
  §7.0).
- `branching.*`: see §6.4–§6.5. `max_label_branches` is a **depth** along a root-to-leaf path (a lineage may reach
  up to d^b leaves, d ≤ number of successor characters), not a total. `tip_window_bp` and `bubble_window_bp` are
  **not implemented**: only `0` (the default, and what the echo carries) is accepted; any other value is rejected
  (400 naming the field) rather than accepted and walked as 0.
- `output.continuation_bp` (default 1000): the length of a leaf's `continuation` (§7.1). `0` means **no
  continuation sequence** (`sequence: ""`; its labels and loss are still reported); `1 … k − 1` is **rejected**
  naming k, because a continuation shorter than k is not valid `/traverse` input; from k on the continuation is a
  resubmittable seed.
- `bounds.*` scopes are in §6.8. `frontier.*` in §6.8. `output.detail`: `summary | tree | full | graphlet`
  (§7; `graphlet` is the lossless compact retrieval format of §7.5: a small JSON summary per seed with the MGT
  text embedded). `output.timing` (default `true`) adds the `timing` blocks. Both are request fields, not walker
  knobs, and the normalized echo carries them like every other field, so the echoed strategy is resubmittable
  verbatim.
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
- `branching.max_label_branches`, `branching.max_splits_per_path` and `output.max_branch_events` accept an
  integer or the string `"unlimited"`, and the normalized echo uses `"unlimited"` where it applies.
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
   `{"seed": {"seed_id", "length_bp", "labels_from_seed": true}, "outcome": {"walks": "failed", …}, "error":
   "<reason>", "limitations": [{"kind": "derivation", "cause", "knob", …}]}` — no `arms` — and the other seeds are
   traversed as usual. The `derivation` entry names the request field that would get past the cause (§7.0); a
   client distinguishes the two shapes by `outcome.walks` (or the presence of `error`). Everything that
   *is* the caller's own doing still fails the whole request with HTTP 400: a malformed seed (too short, invalid
   character, not fully present in the graph), an unknown or duplicate **explicit** label, and every strategy error.
5. The seed is never rewritten.
6. The permitted universe is `P = seed labels ∪ extra`. A seed whose `validated_seed_id` repeats an earlier seed
   of the same request is traversed again and its result carries `duplicate: true` (the `duplicate_of` reference
   and the compute-once reuse are not implemented).

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
  exactness; the approximation is reported, not hidden: every successor derivation in which a cut source had a
  finite switch into one of the successor's targets is counted (`counters.switch_sources_cut`) and stated as the
  arm's `switch_sources` limitation with `outcome.label_evidence: lower_bound` (§7.0) — also when the cut source
  goes on along another successor and so never ends with the `switch_sources` qualifier.

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

**The tip and bubble windows below are not implemented** (§13, I3 note 1): a non-zero `tip_window_bp` or
`bubble_window_bp` is rejected (400 naming the field), no `tip` / `bubble` event is ever emitted and the
`tips` / `bubbles` growth counters stay 0. The two bullets describe the intended behaviour.

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
  at that node with `rejoined_seed`, regardless of seed length (the `{seed_offset, orientation}` qualifier and an
  `overlap_bp` for the re-entered bases are **not implemented**: the label end carries the reason only, and the
  step into the seed node is not taken but reported as a `blocked` event with that reason). Interpretation: a
  circular molecule (offset 0 on the
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
or with a **resource** reason: `max_steps`, `time_budget`, `resource_limit` (a head the request's memory or work
budget did not admit; per seed: end all live paths of both arms), or `max_live_paths`, `max_paths`,
`max_output_bp` (per arm: end that arm only). `beam_pruned` marks pruned paths.

Each arm reports `status: complete | truncated (a resource reason) | pruned (beam)`, `frontier_remaining:
{live_paths, live_labels, exact}` and `cap_trigger: {reason, arm, at_bp, segment, live_paths, live_labels,
exact}` when truncated. `exact` is `false` when a head counted carried a list cut by `labels.max_labels_per_node`
(`annotate` mode): `live_labels` is then a lower bound. The same flag sits on every `growth[]` bin
(`exact`); in `constrain` mode it is always `true` (the live state is never cut). (`cap_trigger` has no `arm`
field: it sits inside its arm.) **Delivery truncation is not implemented** (§13 deviation 8): `max_output_bytes`,
`max_events` and a `delivery: {truncated, next_cursor}` block do not exist, the two fields are rejected, and every
response is delivered whole (`outcome.delivery: inline`, §7.0; spooled / paged delivery is designed in
`DESIGN-traverse-graphlet.md` §14). A capped or pruned result is never evidence of absence.

### 6.8 Exploration order, bounds, determinism

- One priority queue per arm. Keys: `breadth_first` = `(extension_bp, path_id)`; `lowest_loss_first` =
  `(min loss in σ, extension_bp, path_id)`; `most_supported_first` = `(−|σ|, min loss, extension_bp, path_id)`,
  where in `annotate` mode |σ| is the number of labels recorded at the head node (§6.11).
  Arms alternate by `extension_bp` so both advance fairly.
- A popped path advances along its structurally unbranched run in chunks of up to `batch_kmers` nodes (labels
  batch-fetched, §8.3). **Cap and overflow decisions are made in queue-key order at single-step granularity:** a
  chunk is truncated at the step where a cap trips, so `batch_kmers` never changes which steps are taken.
- Scope table: `max_extension_bp` per path; `min_live_labels` per path; `max_steps`, `time_budget_ms`,
  `max_memory_mb`, `max_work_units` per seed (the design's *locus* scope); `max_live_paths`, `max_paths`,
  `max_output_bp` per arm.
- **Head admission** (`DESIGN-traverse-graphlet.md` §14, "atomic commit per head"). Every head is processed in
  three phases: **plan** (successors, derived states, label ends, switch events, splits — computed without
  touching the result), **admit** (the plan's objects and the reservations of the heads it creates — what ending,
  merging and delivering each of them costs — against the memory budget), **commit** (cannot fail). A head that
  is not admitted is censored exactly like a head beyond a cap (`resource_limit`), and its level does not count
  toward `complete_to_bp`; every admitted prefix can be finished and delivered within what it reserved. The
  model is deterministic (fixed bytes per object, never a measurement), so a memory stop is reproducible and
  independent of `annotation.batch_kmers`. Work units are a weighted sum of what the walk consumed — successor
  enumerations (4), annotation keys requested (8) and entries returned (1), pair evaluations, refusal scans,
  edge-reuse probes, derivation scans and steps (1 each) — and the work budget and the deadline are checked
  before every head and at least every `W` = 65 536 units (within a derivation and between chunks of a level's
  annotation fetch), so neither is overrun by more than `W`. The seed phase is charged too (8 per seed k-mer
  read, 1 per annotation entry and coordinate, as the walk charges a fetched row; `account.work_seed` beside
  each arm's `work_units`), checked at the same interval: a seed phase within one interval is followed by a
  stop at the first head (a valid result complete to 0 bp), a longer one is cut within about the interval
  and fails the seed. The depth-0 state is admitted like a head: a memory budget that does not hold it fails
  the seed. Every counter belongs to the arm whose head did the work. The deadline is never checked at depth 0
  (`time_budget_ms: 0` still means one level, then the boundary check). **Not implemented** (see `DESIGN-traverse-graphlet.md`
  §14): the per-seed delivery bounds `max_output_bytes` / `max_events` (§6.7), and a request-level time budget
  (`request_time_budget_ms`) that would mark the seeds it leaves unstarted as `not_started` — every seed of a
  request is traversed, each under its own `time_budget_ms`, and the server bounds a request by
  `--traverse-max-seeds` (§10.3).
- **Overflow** (`on_overflow`): `stop` ends the arm's live paths with that reason; `beam` keeps, at the end of each
  level, the `max_live_paths` heads with the most support (`(−|σ|, min loss, path_id)`; in `annotate` mode |σ| is
  the true count recorded at the head node, §6.11), whatever `frontier.order` is, and ends the others with
  `beam_pruned` (their paths carry `path_reason: beam_pruned` and the arm is `pruned`, with the first pruning in
  `cap_trigger` and a `walk_domain` limitation, §7.0). There is no `beam_rank` knob (the request schema rejects
  it) and no separate `beam_pruned: [...]` list: the pruned paths are the `paths` with that reason.
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
  instead — `label_dict` is filled in the order labels are first met (a level's labels are named when its
  annotation is fetched, before its heads are admitted: a walk that a request budget or a refused admission
  stopped after that fetch keeps only the labels its result records, in the same order, so that no label is
  "met" that appears nowhere in the result; a cap, and a time stop without a budget, keep the dictionary as
  it was), each segment carries the sets present along it as runs (§7.4), and `labels_start` / `labels_end` are the sets at the segment's entry node (the seed
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
  arm), `max_paths` (leaves, per arm), `max_live_paths` (frontier, per arm), `time_budget_ms`, and the request
  budgets `max_memory_mb` and `max_work_units` (per seed, §6.8), which stop between two heads like any cap. Radius is what the
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
  default and the response's `walk_rule` carries the qualification. Every arm states the scope structurally as
  well, `completeness_scope: per_path | united_history` (§7.4), so that a consumer verifying a merged result
  against the exhaustive trie knows which ends it may compare: route support is per-label and scope-free, but
  a structural block (`edge_reuse`, `edge_reuse_rc`, `rejoined_seed` by block precedence) on a route that
  passed a merge is decided on the united history and need not be where the per-path trie ends the label.
- **Two boundaries, one per question.** `complete_to_bp` bounds the *walks* (what is present);
  `evidence.complete_to_bp` (§7.0, §7.2) bounds the *reasons* (why something is absent): below it every
  branch decision and refusal the walker made is reported. They are independent — the first is set by the
  size caps, the second by `output.max_branch_events` — and both are level boundaries for the same reason
  (exploration is level-synchronous, so decisions are made, and events produced, in non-decreasing depth).
  A successor not taken at depth d is guaranteed its recorded reason when d < `evidence.complete_to_bp`; at
  or beyond it the reason may be among the events not kept, and the response says so (a `branch_events`
  limitation). A head the size cap left unexpanded needs no event: it was ended with the cap's reason. The
  tuned-run checker applies exactly this rule (§6.9: such an omission counts as `unexplained_capped`, not as
  a failure, and only at or beyond the boundary).
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
  without constraining by it. Per label, `label_summary.direct_bp` is the longest stretch from the seed boundary
  on which the label is recorded at every node of **some recorded route** — existential over every route the
  result holds, the stubs the beam pruned included (a one-base stub left at a fork gives its carriers
  `direct_bp ≥ 1`), exact when no list was cut. It is label-consistent route support (§7.1), **not** support
  along the spelled walk: how far a sample follows *this* walk is read off the walk's own `label_sets`. A
  stretch supported by no label at all is recorded as an empty set.
- **How the beam ranks.** When the frontier overflows, the beam keeps the heads with the most labels *recorded
  at the head node* — the **true count**, `labels_total`, not the list cut at `labels.max_labels_per_node`
  (in `constrain` mode: the most labels alive) — ties broken by `extension_bp` and `path_id`, so the result is
  deterministic and a width-1 beam follows the majority continuation at each fork, the most-supported local
  path rather than an arbitrary one. **The beam ranks by support whatever `frontier.order` says:** `order`
  only sets the sequence in which the heads of one level are expanded (and so where a size cap trips within
  a level), it never decides which heads a beam retains; there is no earliest-created beam. `breadth_first`
  stays the default only because it is the completeness order; for exploration pass `most_supported_first`
  so that the expansion order matches the beam's.
- **How to use it.** Explore with a beam → read where the support changes → either constrain
  (`labels.mode: constrain`, naming the carriers seen) or re-seed from a `continuation` → extend. The design
  note's "probe cheaply, read the diagnostics, retune" loop; the beam is the cheap probe that reaches far, the
  constrained walk is the one whose output is a claim.

## 7. Response

### 7.0 Guarantees and stated limitations

**Principle.** Every response certifies what it covers and states every cap that limited it, with the request
field to turn: nothing is cut, ignored or weakened silently. A cap either leaves the result complete, or the result
says where the cut starts, what was cut and which knob controls it — the model is `complete_to_bp` (§6.10). A knob
that would be accepted and then change nothing (`tip_window_bp`, `bubble_window_bp`, a `continuation_bp` below k)
is rejected instead (§5).

- **The outcome, per seed result.** Every entry of `results` carries `outcome`, one value per guarantee, because
  the guarantees are independent (a complete walk can have cut diagnostics; complete walks can carry cut label
  lists):

  | axis | values | meaning |
  |---|---|---|
  | `walks` | `complete` \| `partial` \| `failed` | `complete` only when no walk-class limitation applies: every requested arm has `status: complete` (`complete_to_bp == bounds.max_extension_bp`) **per path** (keep, or merge where no merge united a history) and no carrier of the seed was left out. `partial`: a valid certified prefix whose limits are stated — a `walk_domain` (an arm was truncated or pruned), a `scope` with `observed > 0` (a merge united histories: the walks are complete for the united-history rule, not per path) or a `seed_labels` limitation (the walks only a cut carrier carries are missing). `failed`: no traversal — the permitted set could not be derived (§6.1 step 4): the result has no `arms`, an `error`, and a `derivation` limitation; or a request budget does not hold the seed itself (§5, request budgets): no `arms`, an `error`, a seed-level `walk_domain` naming the budget's knob and a `resource_stop` |
  | `branch_diagnostics` | `complete` \| `cut` | `cut`: a `branch_events` limitation — some arm's `evidence.complete` is `false`, and branch decisions and refusals at or beyond `evidence.complete_to_bp` are not reported |
  | `label_evidence` | `complete` \| `lower_bound` \| `qualified` | `lower_bound`: evidence may be *missing or understated* — recorded lists were cut (`label_lists`, `inexact_counts`: annotate mode's `label_summary` and the counts flagged `exact: false` are lower bounds), carriers of the seed were cut (`seed_labels`), a `max_switch_sources` cut may have raised a loss or missed a switch entry (`switch_sources`), or losses were re-minimised greedily (`greedy_losses`). `qualified`: something reported may be *overstated* — a column-label trace cannot see a record boundary (`trace_record_boundaries`); `qualified` wins when both apply |
  | `delivery` | `inline` | the whole result is in this response; `spooled` / `paged` are reserved for the graphlet delivery path (`DESIGN-traverse-graphlet.md` §14) |

  **Conservative by construction** (`DESIGN-traverse-graphlet.md` §14, v5.2): each axis is computed from the
  result's own `limitations` (seed level and every arm), so an axis reads `complete` exactly when no limitation
  of its class is stated — a reader never finds a limitation whose axis still says `complete`, and never a
  non-complete axis without the limitation that explains it. `server_clamp` and `derivation` belong to no class:
  what a clamp caused is stated by the `walk_domain` or `seed_labels` entry it led to. A failed derivation
  reports `{"walks": "failed", "branch_diagnostics": "complete", "label_evidence": "complete", "delivery":
  "inline"}` — nothing was walked, so nothing else was cut — except that carriers cut before a trace check
  (its `seed_labels` entry) make `label_evidence` a `lower_bound`. Each axis is backed by the per-arm fields and
  the `limitations` below, which say how far it holds and which knob would go further.
- **What is certified, per arm.** `complete_to_bp` under `walk_rule` and `completeness_scope` (§6.10): every
  admissible walk of at most that many bases is present. `evidence: {complete, complete_to_bp}`: below
  `evidence.complete_to_bp` every branch decision and every refusal the walker made is in `branch_events`
  (§7.2), so every successor not taken there carries its reason; `{"complete": true, "complete_to_bp": null}`
  when no event was dropped. How the two boundaries relate: §6.10.
- **What limited it.** `limitations`, per arm and per seed result, in every `detail` level; `[]` when nothing
  did. Each entry is `{kind, knob, limit, observed, effect}`, plus `complete_to_bp` for a limitation with a
  depth boundary. `knob` names the request field relative to `strategy` (as `strategy.clamped` does; a seed's
  own field is `seeds[].…`), `limit` is its value in this run, `observed` what the run met against it — for a
  cap, the demand that exceeded `limit` — and `effect` says in one sentence what is missing and how to get it.
  Numbers keep their knob's type: integers for integer knobs (also in `strategy.clamped`), numbers for time
  budgets. A kind is emitted only when it applies:

| `kind` | Where | Emitted when | `knob` | `observed` |
|---|---|---|---|---|
| `walk_domain` | arm; seed, on a seed a request budget failed (no `complete_to_bp`: nothing was walked) | `status` is not `complete`: the cap in `cap_trigger`, then any other cap that ended walks after it (a beam, then a step cap) | the cap's field: `bounds.max_steps`, `bounds.max_live_paths` (also a beam's width), `bounds.max_paths`, `bounds.max_output_bp`, `bounds.time_budget_ms`, `bounds.max_memory_mb`, `bounds.max_work_units` (`resource_limit`) | what the cap compared at the trip: steps of the seed, bases of the arm, leaves + live heads, live heads (beam: heads of the level), elapsed ms, the MiB admitting the head needed (rounded up), the work units used; for a later cap the walks it ended. `complete_to_bp` = the arm's |
| `branch_events` | arm | events were dropped by the cap | `output.max_branch_events` | `branch_events_total`; `complete_to_bp` = `evidence.complete_to_bp` |
| `label_lists` | arm | `labels_per_node.nodes_truncated > 0` | `labels.max_labels_per_node` | `labels_per_node.max_seen` |
| `inexact_counts` | arm | a live-label count in `frontier_remaining`, `cap_trigger` or a `growth` bin is flagged `exact: false` | `labels.max_labels_per_node` | how many counts are flagged |
| `switch_sources` | arm | a `table` cost's source list was cut while a cut source had a finite switch into a target of that successor (`counters.switch_sources_cut > 0`), or a label ended `label_lost` with the qualifier `switch_sources` — whether or not the cut source ended: one that goes on along another successor leaves no end, yet the target it was the cheapest way into was entered at a higher loss or not at all | `labels.max_switch_sources` (accepts `"unlimited"`) | `counters.switch_sources_cut`, the successor derivations the cut may have changed (an over-approximation: the kept sources may still have been the cheapest); `label_ends` = the `switch_sources` label ends. Effect: losses may be overestimated and switch entries missed |
| `scope` | arm | `completeness_scope: united_history` | `branching.on_reconverge` (`limit: "merge"`) | merges done (0: no history was united) |
| `greedy_losses` | arm | an exclusion re-ran the derivation with priced switches under a finite branch limit | `branching.max_label_branches` | re-minimisation rounds after the first |
| `trace_record_boundaries` | seed | `support: trace` with column labels | `labels.seed_label_kind` (`limit: "column"`) | column labels in the dictionary |
| `seed_labels` | seed | the derived permitted set was cut (`labels_dropped > 0`); on a failed result, carriers were cut before the trace check (`no_trace_carrier`) | `labels.max_seed_labels`, with `server_limit` when the server clamped it | `labels_supporting_total` |
| `server_clamp` | seed | an entry of `strategy.clamped` bound this seed: a lowered derived-set cap that cut its set, a lowered time budget that tripped, or a budget raised from zero (the walk ran under it) | the clamped field | the requested value (`limit` is the effective one) |
| `memory_bound_soft` | seed | `bounds.max_memory_mb` is set: the budget is enforced on the modelled state and output, but the annotation rows a level decodes (and an annotate dictionary's growth) are held before they can be charged — stage 3 of `DESIGN-traverse-graphlet.md` §14.1 charges them inside the decoder. In no outcome class | `bounds.max_memory_mb` | the largest excess over the budget seen of what was held beyond the admitted account (the decoded rows, a cache beyond its allotment, the dictionary a fetch grew), MiB rounded up (0: none); the admitted account itself, the depth-0 state included, never exceeds the budget. Also on a seed a budget failed |
| `derivation` | failed seed | the permitted set could not be derived (`outcome.walks: failed`); `cause` names why, `server_limit` is added when the server clamped the knob | per `cause`: `no_carrier` → `seeds[].sequence`; `no_trace_carrier` → `support` (`limit: "trace"`); `too_wide` → `seeds[].sequence`; `time_budget` → `bounds.time_budget_ms`; `ambiguous_header` → `labels.seed_label_kind` (`limit: "header"`); `over_seed_label_cap` (`exhaustive`) → `labels.max_seed_labels` | `no_carrier`: the seed k-mers read when no candidate was left (`limit`: the seed's k-mers); `no_trace_carrier`: the labels carrying every k-mer by presence; `too_wide`: annotation entries of the narrowest of the first 64 k-mers (`limit`: 64 · `max_seed_labels`, at least 65 536); `time_budget`: elapsed ms; `ambiguous_header`: the header; `over_seed_label_cap`: the carriers |

```json
"outcome": {"walks": "complete", "branch_diagnostics": "cut", "label_evidence": "complete", "delivery": "inline"},
"evidence": {"complete": false, "complete_to_bp": 25},
"limitations": [{"kind": "branch_events", "knob": "output.max_branch_events", "limit": 1, "observed": 3,
                 "complete_to_bp": 25, "effect": "branch decisions and refusals at or beyond complete_to_bp are
                 not reported; more were produced: raise the knob or set it to \"unlimited\""}]
```

- **A resource stop**, per seed result, when a request budget stopped the walk (and for a time stop only when the
  request set a budget, so that a response without one keeps its form): `resource_stop: {scope: "locus",
  resource: memory | work | time, phase: "traversal", requested, effective, used, remaining, actions, message}`,
  the amounts in the knob's unit (MiB — `used` rounded up, `remaining` down —, work units, ms; `requested` is the
  request's value where the server clamped it), `actions` the levers (`raise_memory_budget`, `use_graphlet`,
  `drop_sequences`, `raise_work_budget`, `raise_time_budget`, `continue_from_leaves`; on a seed a budget failed,
  which has no leaves, the levers on the seed instead of `continue_from_leaves`: `lower_max_seed_labels` for a
  derived set, `name_fewer_labels` for a named one, `lower_max_labels_per_node` in annotate mode, and
  `shorten_seed` for a work budget the seed phase spent). The MGT `Q` record carries
  it (§7.5.2). It is a property of the walk, the same in every detail; a memory stop itself depends on the
  requested detail (the output is charged), so the same request can stop at another depth in another detail.
  `phase` is always `traversal` in this stage: finalisation and serialisation are reserved at admission.
- **Not limitations:** semantic stops (`dead_end`, `label_lost`, `max_extension_bp`, …) — the requested domain
  is complete there — and `dropped_labels` (named labels that do not support the seed). A seed whose permitted
  set could not be derived is a `failed` result with a `derivation` limitation (§6.1 step 4), not a limited one.
- So an arm with `limitations: []` is exactly what its strategy defines to `bounds.max_extension_bp`: every
  walk per path, every reason, every recorded list in full, every loss as the §6.3 recurrence gives it, and a result whose `outcome` axes are all `complete` / `inline` and whose seed-level `limitations` are
  `[]` covers everything it certifies.
- **Formerly unstated gaps, now stated per response:** (1) `support: trace` with **column** labels follows consecutive
  column coordinates and cannot detect a record boundary whose global coordinates are consecutive (header labels
  detect it): stated as the seed-level limitation `trace_record_boundaries` (knob `labels.seed_label_kind`,
  observed = the number of column labels), with `outcome.label_evidence: qualified`. (2) Under `support: trace`, `dropped_labels[].runs` are k-mer presence
  runs: every dropped label carries `runs_kind: "presence"`. (3) Under a finite change cost with a finite
  `max_label_branches`, losses after a branch-limit exclusion are re-minimised greedily (§6.3, §6.4) and need not be
  optimal: stated as the arm limitation `greedy_losses` (knob `branching.max_label_branches`, observed = the
  re-minimisation rounds after the first) whenever an exclusion actually re-ran the derivation with priced switches,
  and `outcome.label_evidence` is then `lower_bound`. No known gap remains unstated.
- **Described in this spec but not implemented** (a request naming one is rejected; no response carries one):
  delivery bounds and the `delivery` block (§6.7, §13 deviation 8), a request-level time budget and `not_started`
  (§6.8), `beam_rank` and a `beam_pruned: [...]` list (§6.8), tip and bubble windows (§6.5), the `hll` cost model
  (§9), `duplicate_of` (§6.1; the response has `duplicate: true`), `seed_offset` / `orientation` / `overlap_bp` on
  `rejoined_seed` (§6.6), the `label_lists: delta` encoding (§7.1; superseded by the graphlet's delta-coded
  sets, §7.5), and in §7.2 `cost_preview`, the
  `rc_index_range` and resource-cap-hit counters and the `bubble_len_bp` of a `reconverge` event.

### 7.1 Per seed (`detail: full`)

```json
{
  "seed": {"seed_id": "…", "validated_seed_id": "…", "seed_id_mismatch": false, "length_bp": 4870,
           "num_kmers": 4840, "labels": [...], "dropped_labels": [],
           "labels_from_seed": false, "labels_supporting_total": 2, "labels_dropped": 0,
           "labels_dropped_digest": ""},
  "outcome": {"walks": "complete", "branch_diagnostics": "complete", "label_evidence": "complete",
              "delivery": "inline"},
  "label_mode": "constrain",
  "label_dict": [{"name": "NZ_STEQ01000045.1", "kind": "header"}, {"name": "573", "kind": "column"}],
  "limitations": [],
  "arms": {
    "right": {
      "status": "complete", "frontier_remaining": {"live_paths": 0, "live_labels": 0, "exact": true},
      "evidence": {"complete": true, "complete_to_bp": null},
      "limitations": [{"kind": "scope", "knob": "branching.on_reconverge", "limit": "merge", "observed": 0,
                       "effect": "…"}],
      "segments": [
        {"id": 0, "parents": [], "from_bp": 0, "length_bp": 812, "sequence": "…",
         "labels": [0, 1], "labels_at_end": [0],
         "events": [{"at_bp": 300, "type": "switch", "from": 1, "to": 2, "cost": 0.7},
                    {"at_bp": 455, "type": "label_end", "label": 2, "reason": "label_lost",
                     "structural_successors": 1}]},
        {"id": 1, "parents": [0], "from_bp": 812, "length_bp": 40, "labels": [], "labels_at_end": [],
         "sequence": "…", "events": []}
      ],
      "splits": [{"at_bp": 812, "segment": 0, "children": [1, 2], "kind": "divergence"}],
      "paths": [{"id": 0, "segments": [0, 1], "length_bp": 852, "end_reasons": {"dead_end": 1},
                 "end_labels": [{"label": 0, "loss": 0, "branches": 0, "run": 0, "route_bp": 0}]}],
      "runs": [{"label": 0, "from_bp": 0, "to_bp": 852, "entered_by": "seed", "route_bp": 0,
                "end_reason": "dead_end"}],
      "growth": [], "branch_events": [], "branch_events_truncated": 0, "needed_budgets": [],
      "counters": {"…": "§7.2"}
    },
    "left": {"…": "same shape"}
  },
  "label_summary": [{"label": 0, "right": {"direct_bp": 852, "reach_bp": 852, "reentries": 0, "runs": [0]},
                                  "left": {...}}],
  "annotation": {"access_path": "…", "keys_mapped": 0, "rows_requested": 0, "direct_reads": 0},
  "timing": {"…": "§7.2"}
}
```

- **A seed whose permitted set could not be DERIVED** (§6.1 step 4) takes the place of this object with
  `{"seed": {"seed_id", "length_bp", "labels_from_seed": true}, "outcome": {"walks": "failed", …}, "error":
  "<reason>", "limitations": [{"kind": "derivation", …}]}` and no `arms`: the request still succeeds and the other
  seeds are traversed. `outcome.walks` (or the presence of `error`) is how a client tells the two apart.
- `segment.labels` = names live at the segment's **first** step; `labels_at_end` at its last (both lists in full:
  a `label_lists: delta` encoding with `labels_removed` is not implemented). Label sets change only through
  events, sorted by `at_bp`.
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
  input (a fresh seed: edge-reuse and branch state reset; the agent may reduce `loss_budget` by `loss_used`):
  `continuation_bp` is either 0 or at least k (1 … k − 1 is rejected, §5). With `continuation_bp: 0` the sequence
  is `""` and `labels` / `loss_used` / `branches_used` describe the leaf's head node.
- `detail: summary` returns `seed`, `outcome`, `limitations`, `label_summary`, per-leaf `{length_bp, n_labels,
  end_reasons}`, the per-arm diagnostics of §7.2 and `continuation`s, but no segments, splits, runs, events or
  sequences — the cheap probe to run before committing to a radius.
  `tree` adds segments without sequences. `graphlet` returns per seed the summary of §7.5 and the whole result
  as MGT text (`graphlet`), from which the `full` result is reproduced exactly (§7.5.6).
- **Fields that make the result reconstructible** (additive; no existing field changed):
  `runs[]`: `segment` — the segment on which the run ended (its `label_end` event, if it has one, is there at
  `to_bp`) or was closed by a merge (a parent of the merged segment); `branches`, `loss` — the lineage's
  values at that point. A run can end with **no** event: a switch source whose lineage goes on only under
  other names ends `label_lost` silently, documented by the `switch` event(s) at `to_bp` (on the children when
  the source ends at a split; §13 note 20). `segments[].labels_via_parent` on a merged segment (more than one
  parent): per parent, the labels whose kept entry came through it (empty lists in `annotate` mode).
  `hairpin` events: `followed` (`hairpins: follow`, the child is present) vs skipped. `revisit` events:
  `length_bp` (the distance difference to the first arrival) and `same_distance` (a `keep`-mode join at one
  node and depth). `label_dict[]`: `column`, and `seq_id` for header labels — the `LabelRef` that joins a label
  across retrievals (a name alone can be ambiguous). `growth[].label_ends` is an object even when empty (was
  `null`).
- **Derivable fields** (the graphlet stores none of them): `children`, `paths`, `splits`, `end_labels` and
  `end_reasons` of a leaf (its labels' runs), `label_end` events (the runs), `reconverge` events (the merged
  segments), `labels_at_end`, `label_summary`, `needed_budgets` (the `loss_budget` ends; a histogram, its
  order is not information), `branch_events_truncated`, every `*_truncated` flag. **Events at one `at_bp` are
  unordered**: the list is sorted by position only.

### 7.2 Diagnostics (the tuning evidence, always present)

- `growth[arm]`: per bin of `profile_bin_bp`: max live paths, distinct live labels, live (path, label) pairs,
  divergences, ambiguous branches taken, splits, reconvergences, bubbles, tips, `blocked_repeat`, label ends by
  reason, steps.
- `branch_events[arm]`: the first `output.max_branch_events` events in processing order (default 100, or
  `"unlimited"`). Exploration is level-synchronous, so an arm produces its events in non-decreasing `at_bp`
  and the kept prefix ends at a depth: the arm's `evidence.complete_to_bp` (§7.0) is the `at_bp` of the first
  event **not** kept, and every event at a smaller depth — every ambiguity and every refusal decided there — is
  in the list (at that depth and beyond some may be missing; a `branch_events` limitation then states the cut
  with `branch_events_total` as `observed`). Each event has its position, successor characters, per-successor
  label counts, ambiguous labels, labels dropped by `branch`/`minority`/`loss_budget`, and
  `refused: [{char, cause, labels}]` — one entry per (successor, cause) the walker decided **not** to follow for
  labels whose lineage would have continued on it, with those labels (the sources alive at the node,
  ascending). The causes and when each is emitted:

  | `cause` | Emitted | Knob that must be active |
  |---|---|---|
  | `branch` | for a source ambiguous at the node (its lineage on two or more successors) over its allowance: on every successor that carried its lineage when it was excluded; it ends there with `branch` | finite `max_label_branches` |
  | `minority` | for every source on a successor whose label count (after the branch-limit exclusions) is below `min_successor_labels` or `min_successor_fraction · \|σ\|` | either quorum knob |
  | `below_min_labels` | the same, for a count below `bounds.min_live_labels` | `min_live_labels > 1` |
  | `split_limit` | for every source on every successor of a split the path may not make (it has made `max_splits_per_path` splits) | finite `max_splits_per_path` |
  | `loss_budget` | for every source on every successor where its lineage could go on only by a switch above the budget into a target nobody entered — **whether or not the source continues on another successor or ends for another reason** (only its *label end* `loss_budget` depends on that). Not for a source excluded by the branch limit: the limit took it out of every successor's source set, so none of its switches was priced; its refusals are the `branch` ones | a finite change cost |

  This is the explicit per-successor refusal evidence the tuned-run property (§6.9) is checked against:
  `dropped` and `ambiguous` alone cannot tell a refused successor from one that was followed and then lost
  from the output, and a missing branch says nothing about why it is missing. Because a refusal is a claim
  about a pruning decision, a consumer checks it against the strategy the run was made with: the cause must
  name an active knob whose condition holds at that node (the test checkers do, `refusal_problem`), and a
  `hairpin` event is a reason for a missing successor only under `hairpins: skip` — one marked `followed`
  records a child that is present. A step with a refusal always emits a branch event; `branch_events_truncated`
  = total − shown. Refusals live in the events, so a cut by `max_branch_events` cuts refusal evidence — but
  only from `evidence.complete_to_bp` on, and the response says where: below it an omission must carry its
  reason, at or beyond it a consumer cannot expect one (the checker counts such an omission as
  `unexplained_capped`). Set the knob to `"unlimited"` when the evidence for the whole trie is what is wanted.
- `needed_budgets` per arm: the `needed_budget` of every `loss_budget` end, as a list. (A `cost_preview` for P
  under pairwise models — min/median/max of `cost(seed label → extra)` — is not implemented.)
- `counters` per arm: `steps`, `successor_enumerations`, `output_bp`, `pair_evaluations`; the work of the
  structural and branching rules — `edge_reuse_probes` (uses scanned or (edge, segment) pairs probed by the
  per-path edge-reuse check, whichever was fewer), `reminimisation_rounds` (re-derivations after the first at
  ambiguous nodes, §6.4; bounded by |σ| per node), `max_reminimisation_rounds` (the largest at one node) and
  `refusal_scans` (the work of recording the refusals: per round that excludes sources, the successor state
  entries scanned once for all of them, plus one loss-budget test per (source, successor) under a constant cost
  or per target under a table — linear in |σ| + Σ|σ_v| per round under `forbid` and `constant`) — so that a
  pathological locus is visible rather than silent; and `switch_sources_cut`, the successor derivations of
  committed steps in which a `table` cost's source list was cut by `max_switch_sources` while a cut source had a
  finite switch into a target there (the `switch_sources` limitation, §7.0); and, when the request sets a budget,
  `work_units`, the arm's charged work units (§6.8; listed last, so a reader that does not know it keeps it).
  Every counter belongs to the arm whose head did the work. Annotation access is per seed
  (`annotation: {access_path, keys_mapped, rows_requested, direct_reads}`); `rc_index_range` calls and
  resource-cap hits are not counted.
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

- `complete_to_bp` (§6.10) next to `status`, `frontier_remaining` and `cap_trigger`, with
  `completeness_scope: per_path | united_history` naming what it quantifies over (`per_path` under
  `on_reconverge: keep`, `united_history` under `merge`) so that a checker can pick the termination scope it
  verifies against instead of parsing `walk_rule`. `frontier_remaining` and `cap_trigger` count the heads a cap
  cut, not heads that had already reached the radius (those end as `max_extension_bp`).
- `evidence` (the boundary of `branch_events`) and `limitations` (every cap that limited the arm, with the knob
  to turn), §7.0.
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

### 7.5 The graphlet (MGT v1)

**MGT v1 is frozen (2026-10-02).** This section is the wire contract; any change to a record, field or token is MGT v2.

`output.detail: graphlet` returns, per seed, a small JSON summary with the **whole** result embedded as one
line-based text, MGT v1 (`DESIGN-traverse-graphlet.md` §2, as implemented; where this section and the design
differ, the differences are the decisions listed in §7.5.8). The text is lossless: the `detail: full` result of the
same request is reproduced from it and its summary (§7.5.6). It is written by a streaming writer straight from the
walker's result (no JSON tree), is self-contained (seed, dictionary, provenance and orientation rule are in it). Measured with
`scripts/traversal/graphlet_measure.py --cli` on the mini-refseq fixture (the whole blaNDM-1 seed, nine strategies,
200–1000 bp): the body is 11–29 % of the compact `full` JSON (23–60 % gzipped against gzipped), much of it the
bases themselves; the largest case, annotate at 600 bp (2,020 segments, 1,012 leaves), is 210 KB against 1.90 MB
(11 %; 16 KB against 70 KB gzipped) and parses in 20 ms. **Size targets: unverified.** The design's targets for
constrained retrievals (< 15 % of the compact `full` JSON, ~7.5 % projected for the 2,500-leaf SRA trie, gzip
~33–50 %) are **not met on mini-refseq**: every constrained case is 18–29 % (gzip 48–60 %; merge 600 bp 22 %,
exhaustive 1000 bp 23 %, in which the fixed records `H S A K L O Z` are 8 % of the body and the bases ~70 %), so
the miss is not a small-locus effect of fixed costs, and the 7.5 % projection rests on SRA's ~24 bp per
segment, which has not been measured. Met: annotate (~10 % target; 11 %), the per-seed summary (< 10 KB;
1.5–2.4 KB), T40's body < ½ of `full`; per-record costs match the estimates (G ~26 B per segment, R ~27 B per
run). The figures are not to be quoted as met until `graphlet_measure.py --server <SRA>` has run, and the MGT v1
freeze is gated on that measurement; if SRA misses too, the bases (`G`) are the lever (optional 2-bit packing,
or bases only on request with spelling through `graphlet_sequence`). A seed whose permitted set could not be
derived keeps the failed shape of §7.1 (no `graphlet`).

#### 7.5.1 Lexical rules

- One record per line, terminated by LF (the only line separator: names may hold VT, FF, NEL, LS, PS raw). Fields
  are separated by one space; the first field is a capital letter, the record type. Every field is printable
  ASCII without spaces, except the **last field of `X`, `L`, `Q` and `K`** (a name, a suffix, a message, an
  effect), which is free text running to the end of the line, preserved exactly: in it `%`, LF and CR are
  written `%25 %0A %0D` (uppercase) and nothing else is escaped; any other `%` sequence or a raw LF/CR is
  invalid.
- **Integers** decimal, no sign, no leading zero. **Floats** (costs, losses, budgets, demands; all ≥ 0): the
  shortest round-trip significant digits of the double (C++ `std::to_chars(..., chars_format::scientific)`
  without a precision = Python `repr`) **expanded positionally**: zeros up to the decimal point, a `.` only before
  a non-empty fraction, no trailing zeros, `-0` → `0`, +∞ → `inf` (`1e23` → `100000000000000000000000`,
  `2.5e-07` → `0.00000025`). Booleans `0|1`.
- `*` = "absent" or "derive by this field's rule" (§7.5.3); `.` = the empty list. Lists are comma-separated.
- **RANGES**: strictly ascending ids, every maximal run of ≥ 2 consecutive ids written `a-b` (`{5,6}` is `5-6`),
  a lone id `a`; `.` = empty. A reader checks every id against the id space (the number of `L` records) before it
  expands a run, so a corrupt token of a few bytes (`0-3000000000`) fails without allocating in proportion to
  the range.
- **SETEXPR** against a base set: `RANGES` (explicit) | `!` (= the base) | `!R` | `!+A` | `!R+A` (base minus R
  plus A; R and A non-empty RANGES; an empty part is omitted, never written `.`). The writer emits the shorter in
  bytes, a tie explicit; a delta needs a base, its R must lie in the base and its A outside it (else the reader
  derived another base: an error). Bases: `G.entry` → parent[0]'s end set, for a merge the union of all its
  parents' end sets, the root none (explicit only); `G.end` → this segment's entry set; `P` → the previous `P`
  of the segment, the first `P` the entry set; `C.labels` → the leaf's end set.
- **Typed values** (`K`): `i:<integer>` (0 … 2⁶⁴−1) | `f:<float>` | `s:<string>` | `u` (unlimited). In `s:` every
  byte outside 0x21–0x7E, and `%` and `,`, is `%XX` (uppercase hex), every other byte raw; the decoded bytes must
  be UTF-8. The string `unlimited` is `s:unlimited`; the knob value "unlimited" is `u`.
- **Codes.** End reasons: `D` dead_end, `L` label_lost, `B` loss_budget, `R` branch, `U` edge_reuse,
  `V` edge_reuse_rc, `J` rejoined_seed, `T` trace_break, `X` max_extension_bp, `S` max_steps, `P` max_live_paths,
  `N` max_paths, `O` max_output_bp, `M` time_budget, `W` beam_pruned, `Y` resource_limit. A label
  end's text qualifier is a second letter: `Rm` minority, `Rb` below_min_labels, `Rs` split_limit, `Dh` hairpin,
  `Ls` superseded, `Lx` switch_sources.
- Every token has exactly one valid spelling and readers reject every other one (`1.0`, `01`, `0,1`, `%41`, a
  `!` without a base, …); the only freedoms are field-level (`*` vs an explicit value, explicit vs delta
  SETEXPR, a shorter front-coding prefix), which parse and dump in canonical form. The shared golden vectors
  `api/python/tests/data/traverse/codec_vectors.tsv` (10,433 vectors; generated by
  `scripts/traversal/make_codec_vectors.py`) are read byte-exactly by the C++ test
  (`GraphletCodec.GoldenVectors`) and the Python library.

#### 7.5.2 Records, in document order

`H S X* L* O Q? K*`, then per requested arm, left before right, `A B* V*`, then per segment in id order
`G P* E* T? C?`, then the arm's `R*`; finally `Z`. (`J`, the response envelope of a saved file, is written
by the library after `H`, never by the server.)

```
H mgt 1 <k> <regime> <alphabet> <mode c|a> <support k|t> <reconverge m|k> <cap> <continuation_bp>
  <seed_index> walk <index_ns|*> <index_fp|*> <index_meta_fp>
    cap = labels.max_labels_per_node; seed_index = the seed's position in the request; walk = the orientation
    rule (§7.5.4); the index identity of §10.3 (* = no name / no manifest)
S <validated_seed_id> <length_bp> <num_kmers> <num_seed_labels> <SEQUENCE>
    the seed as validated: case-mapped, request orientation
X <reason> <a-b,c-d,…|.> <name>
    a dropped seed label: its reason and its k-mer presence runs on the seed, as given ([a, b] pairs)
L <c|h> <column> <seq_id|*> <prefix_len> <suffix>
    one per label_dict entry, label id = ordinal. name = previous name[:prefix_len] + suffix: prefix_len counts
    BYTES of the previous name's UTF-8, cut back to a code-point boundary, so the suffix is valid UTF-8; the
    first L's previous name is "". The space before the suffix is always written (an empty suffix leaves the
    line ending in a space). seq_id is * for a column label. Names come from FASTA headers and file names,
    which need not be UTF-8: the server makes every label and dropped-label name valid UTF-8 once, before any
    output (each maximal ill-formed subsequence → U+FFFD, = Python's decode('utf-8', 'replace')), so
    detail: full, the summary and the body carry the same bytes and graphlet_bytes holds after transport
O <walks c|p|f> <branch_diagnostics c|x> <label_evidence c|l|q> <delivery i|s|p>
    the per-seed outcome (§7.0); q = qualified (DESIGN v5.2)
Q <scope> <resource> <phase> <requested VALUE|*> <effective VALUE|*> <used VALUE|*> <remaining VALUE|*>
  <actions csv|.> <message>
    a resource stop (DESIGN §14): the result's resource_stop (§7.0), written when a request budget stopped
    the walk. The four
    amounts are typed like K's values (the design left them untyped; a budget can be an integer, a float or
    "unlimited"), * = not stated; actions = the suggested next actions ([a-z_] tokens); message = free text
K <arm l|r|*> <kind> <knob> <limit VALUE> <observed VALUE> <complete_to_bp|*> <extra name=VALUE,…|.> <effect>
    one entry of `limitations` (§7.0), the JSON entry itself: arm * = seed level. Seed-level entries first,
    then the left arm's, then the right arm's, each list in its JSON order. extra = every further key of the
    entry (server_limit, label_ends, …) in byte order of name; effect = the stored sentence (free text)
A <l|r> <status c|t|p> <complete_to_bp> <scope p|u> <live_paths> <live_labels> <exact> <max_seen>
  <nodes_truncated> <counters> <cap_trigger|*> <branch_events_total> <evidence_complete_to_bp|*>
  <n_segments> <n_runs> <n_leaves> <n_splits> <n_merges> <bases>
    counters = name=value,… in the walker's order (steps, successor_enumerations, output_bp,
    pair_evaluations, edge_reuse_probes, reminimisation_rounds, max_reminimisation_rounds, refusal_scans,
    switch_sources_cut); a reader keeps a name it does not know. cap_trigger =
    reason,at_bp,segment,live_paths,live_labels,exact,demand. evidence * = no branch event dropped. The last six
    validate the body: G records, R records, T records, G records with a split, G records with > 1 parent, and
    the sum of G.length_bp
B <from_bp> <max_live_paths> <distinct_live_labels> <live_pairs> <exact> <steps> <divergences>
  <ambiguous_branches> <splits> <reconvergences> <bubbles> <tips> <blocked_repeat> <label_ends code:n,…|.>
    one per growth bin; label_ends in end-reason order. An arm may have none: when a cap trips on the left
    arm at depth 0, the right arm is stopped before it ran a level (its growth is [], every count 0)
V <at_bp> <segment> <chars|.> <labels_per_successor csv|.> <ambiguous RANGES> <dropped RANGES>
  <refused char:cause:RANGES;…|.>
    the stored branch events (the only record of a successor not taken)
G <parents csv|*> <from_bp> <length_bp> <entry SETEXPR|*> <entry_total|*> <end SETEXPR|*>
  <partition RANGES|RANGES|…|*> <split 0|1|*> <first_base|*> <bases|.|*>
    segment id = ordinal (parents are created first). split = Split::ambiguous of the split at this segment's
    end, * = it does not end in a split. bases in WALKING order (§7.5.4), . = none (length 0), * = absent
    (`sequences: false`); first_base only without bases: a split child's first base, * otherwise
P <from_bp> <to_bp> <total> <SETEXPR>
    annotate mode: one recorded LabelSetRun
E <at_bp> s <from> <to> <cost>                      switch
E <at_bp> b <char> <reason code> <total> <RANGES>   blocked successor
E <at_bp> h <char> <total> <RANGES> <f|s>           hairpin followed | skipped
E <at_bp> v <segment> <length_bp | =>               revisit: distance difference, = same_distance
E <at_bp> t <char> <length_bp>   E <at_bp> u <length_bp> <alleles>     tip / bubble: reserved
    the segment's events in stored order, without label_end (from R) and reconverge (from G)
T <path_reason code|*> <label:loss:branches:route_bp,…|.>
    this segment is a leaf; * = no path reason (semantic end); the extras list the leaf labels whose loss,
    branches or route_bp (the entry's latest merge stamp) is non-zero
C <n> <loss_used> <branches_used> <labels SETEXPR> [<sequence>]
    the leaf's continuation; the sequence is written only where it is not derived (§7.5.3)
R <segment> <label> <from_bp> <to_bp> <end> <route_bp> <from_label:cost|*> <prev_run|*>
  <structural_successors|*> <branches> <loss> [<needed_budget>]
    run id = ordinal (ArmResult::runs order). segment = LabelRun::segment. end = code[+qualifier] (ended with
    a label_end event on <segment> at to_bp) | Lw (ended label_lost silently: a switch source that went on
    under other names) | m (closed by the merge at to_bp, ended = false). from_label:cost * = entered by seed;
    structural_successors * for Lw and m; needed_budget only for code B. branches, loss = the lineage's
    values when the run ended or was closed (§7.1)
Z <line count of the document, this line included>
```

#### 7.5.3 The `*` rules (canonical form)

The writer computes each rule from the walker's own data and writes `*` **exactly when** the rule reproduces the
walker's value, the explicit value otherwise; every SETEXPR in its shortest form. So every document the server
writes is canonical, `dump(parse(x)) == x`, and every dump tests the rules. A value the grammar cannot express
(a derivation that does not hold) fails the request with a 500 instead of producing a wrong document.

- `G.parents *` = the root. `G.entry *` (constrain only): the root → `0 … num_seed_labels−1`; a merge → the union
  of its partition. Never `*` in annotate mode.
- `G.entry_total *` = |entry|.
- `G.end *`: constrain → start from the entry set and apply, ordered by position with removals before additions
  at one position, every `R` anchored on this segment with `to_bp < from_bp + length_bp` (remove its label) and
  every `E … s` on it (add its `to`); annotate → the last `P` set, or the entry set without `P`.
- `G.partition *` = no merge (|parents| ≤ 1), or annotate mode, where it means one empty list per parent.
- `T.path_reason *` = none; `R.from_label:cost *` = entered by seed; `R.prev_run *` = none.
- `C.sequence` absent = derived: right arm the last n bases of seed + natural(right flank), left arm the first n
  bases of natural(left flank) + seed (n ≤ |seed| + flank; it may cross into the seed). It is written when the
  flank's bases are absent (`sequences: false`) and n > 0.

#### 7.5.4 Orientation

`walk` in `H`: every `G` stores its bases in walking order on **both** arms, so outward index
`i ∈ [from_bp, from_bp + length_bp)` is `bases[i − from_bp]` and every position field indexes bases directly.
Natural orientation: right flank = concatenation root → leaf; left flank = reverse(concatenation root → leaf)
(= today's `sequence` fields, which hold each left segment reversed). Whole molecule: natural(left) + seed +
natural(right); seed coordinate of outward `i`: |seed| + i (right), −(i+1) (left).

#### 7.5.5 The summary (per seed)

```
{"seed": {seed_id, validated_seed_id, seed_id_mismatch, length_bp, num_kmers, labels_from_seed,
          labels_supporting_total, labels_dropped, labels_dropped_digest, num_labels, num_seed_labels},
 "label_mode", "duplicate"?, "outcome", "limitations",
 "arms": {"left"?|"right"?: {status, complete_to_bp, completeness_scope, frontier_remaining, labels_per_node,
          cap_trigger?, evidence, limitations, counters, branch_events_total,
          counts: {segments, leaves, splits, merges, runs, bases, max_bp,
                   label_ends: {reason: n},            run ends by reason, silent ends included (= growth)
                   leaves_by_reason: {path_reason | "semantic": n}}}},
 "annotation", "timing"?, "graphlet": "<MGT text>", "graphlet_bytes", "graphlet_lines"}
```

No names, per-segment, per-leaf or per-label data (they are in the text). `timing` (per seed and per response)
gains `serialize_ms`, the time spent building the response for it (JSON and text).

#### 7.5.6 Stored and derived

Stored: bases, the DAG, entry sets where not implied, merge partitions, annotate presence runs, runs with their
anchor, end token, `route_bp`, switch source and cost, `prev_run`, structural successors, terminal branches and
loss and `needed_budget`, leaf path reasons with their labels' (loss, branches, route_bp), continuations, the
switch/blocked/hairpin/revisit events, growth bins, branch events, seed, dropped labels, the dictionary,
outcome and limitations. Derived exactly as the walker's finalisation: `children`; leaves (segments with `T`);
`paths` (leaf ordinal in segment order, chain through `parents[0]`, `length_bp = from_bp + length_bp`);
`end_labels` (the `R` anchored at the leaf with `to_bp` = the path length, ascending label, with the `T` extras)
and `end_reasons`; `label_end` events (`R` with a code end, at `to_bp`, reason = qualifier text or the enum,
`needed_budget` for `B`); `reconverge` events (a `G` with > 1 parents, at `from_bp`); `splits` (single-parent
children grouped by parent, ordered by (at_bp, first child); `kind` from `G.split`; `labels_before` = |parent end
set| (constrain) or the parent's last `P` total, else its entry total (annotate); per branch: first base,
`labels_distinct` = the child's entry total, `labels` = its entry set cut to `cap`); `labels_at_end` = `G.end`;
`label_summary` (the normative pseudocode of `DESIGN-traverse-graphlet.md` §5.1); `needed_budgets`;
`branch_events_truncated`; continuation sequences; `seed.labels`; every `*_truncated` flag (total > |list|);
`dropped_labels[].runs_kind` (`presence`). `Graphlet.ReaderRebuildsTheFullResult` (C++) rebuilds the `full`
result from the text and the summary and compares it with the server's own on merged, kept, annotate, beam,
switching, capped, cut-list, hairpin, `sequences: false` and batch fixtures; the Python library
(`api/python/metagraph/traverse`) does the same on every body of T38–T41 (§11.3) and on the whole-document
fixtures `api/python/tests/data/traverse/documents/` (most generated by the CLI from tiny indexes,
`scripts/traversal/graphlet_fixtures.py --from-cli … --check`).

#### 7.5.7 Index identity

`H` carries the identity of §10.3: labels are joined across retrievals by `(kind, column, seq_id)` only when the
`index_fp` agree; a missing `index_fp` makes a join unverifiable, never equal; `index_meta_fp` can only prove two
indexes different.

#### 7.5.8 Decisions where the design was silent or stale

- `H` has three identity fields (the design's worked excerpt §2.6 still shows two); `O.label_evidence` admits
  `q` (v5.2's `qualified`, written under `trace_record_boundaries`, §7.0).
- `G.bases` `.` for a zero-length segment when sequences are on (`*` would not say whether `sequence` is `""` or
  absent); `first_base` is written only for split children (the only first bases today's JSON reports with
  `sequences: false`).
- The `G.entry` base of a merge is the union of its parents' end sets in **both** modes (§2.1 of the design;
  its §2.3 mentions parent[0] for annotate).
- `C.sequence` is written whenever it is not derivable or differs, not only "when G carries no bases".
- `K`: a JSON string "unlimited" is `u`, every other string `s:`; a JSON integer `i:`, a JSON real `f:` (a
  floating knob stays a float even when integral: `f:30000`); extras in byte order of name; seed-level
  entries first.
- `T` keeps the grammar's extras although `R` now carries the same terminal loss and branches (only
  `route_bp` differs from the run's).
- `R`: a single-letter end code when the label end has no qualifier; `Lw` and `m` as in the grammar.

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
  `schema_version`, `release`, `graphlet_format` (1: the MGT version `detail: graphlet` writes, §7.5),
  `detail_levels` (`["summary", "tree", "full", "graphlet"]`); `GET /traverse/capabilities` adds the request
  budgets it accepts (`budgets: ["max_memory_mb", "max_work_units"]`), `work_check_interval` (`W` of §6.8, in
  work units) and `memory_bound: "soft"` (until stage 3 charges annotation decoding; not in the per-request
  capabilities, which stay as they were); and the **index identity**
  (`DESIGN-traverse-graphlet.md` §3.1), also in every graphlet's `H` record:
  - `index_ns`: `--index-name NAME` (`[A-Za-z0-9._-]+`), a name for humans and routing, not identity; `null`
    when unset.
  - `index_fp`: the identity — the lowercase hex sha256 of the index bundle's **manifest** file list
    (`--index-manifest FILE`): a JSON object with `files: [{path, size, sha256}, …]` (every file the server
    loads: graph, annotation and sidecars such as `.anchors`, `.rd_succ`, `.seqs`), hashed as the lines
    `<path>\t<size>\t<sha256>\n` in ascending byte order of path; other keys (builder, inputs) are metadata. The
    build hashes the files once; the server reads the digests, checks that every loaded file is listed (by base
    name) with its size, and that an `index_fp` the manifest states is the computed one, and refuses to serve
    otherwise. `null` without a manifest: labels joined across retrievals are then unverifiable, never equal.
    `scripts/traversal/build_mini_refseq.sh` writes `<out>/annotation.relaxed.relabeled.manifest.json`
    (`MANIFEST_ONLY=1` writes it for an existing build).
  - `index_meta_fp`: FNV-1a-64 (16 hex) over k, regime, alphabet, the annotation's row count and its ordered
    column names, and whether coordinates and a `CoordToHeader` exist — always present, and only a **negative**
    check: a mismatch proves two indexes different, equality proves nothing (an index whose columns swap their
    memberships keeps it).
  Both flags describe one index (`-i`/`-a`) and are refused with a graph list; in multi-graph mode `index_ns` and
  `index_fp` are `null` and `index_meta_fp` is computed per index on first use. The CLI (`metagraph traverse`)
  takes the same two flags.
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
| T30c | per-seed derivation failure | Three seeds, the middle one carried in full by no label: HTTP 200 / exit 0 with `{"seed": {...}, "outcome": {"walks": "failed", …}, "error": …, "limitations": [derivation]}` in its place (cause `no_carrier`, knob `seeds[].sequence`) and the other two traversed with their own `outcome` (`Walker.FailedDerivationIsAStatedOutcome`, also `over_seed_label_cap` and `time_budget` with `server_limit`); the same seed with an EXPLICIT label list still fails the whole request; `max_seed_labels` out of `[1, 100000]` is a parse error; a seed over `--traverse-max-seed-bp` and more than `--traverse-max-seeds` seeds are 400s; `max_seed_labels` above `--traverse-max-seed-labels` is clamped and reported | Hit-labelled seed |
| T31 | trie oracle (§6.9) | `test_trie.cpp`: annotate records what constrain filters; the exhaustive preset refuses conflicting knobs; a tripped cap reports `complete_to_bp` and the partial level is excluded; cut recorded lists are reported and break the oracle; `TrieOracle` (4 graph × annotation pairs, 3 modes): `claims(A) == E` for 5 permitted sets on both arms, tuned runs are prefix-subsets with a reason for every omission, merged routes are sound | — |
| T32 | records model, edge and non-edge cases | `test_trie_cases.cpp` + `test_trie_reference.hpp` + `test_trie_checks.hpp`, two graph × annotation pairs, three modes, both arms, every case three ways (records ⇔ T ⇔ A with end reasons) plus single-label union, derived == explicit and the recurrence at budget 0: linear; seed at record start / end / the whole record; radius 0, 1, 2, |R|−1, |R|, |R|+1 (a head at the radius is `max_extension_bp`, not `dead_end`); fork; unequal bubble (+ merged routes); nested bubbles; one-base tip at the seed boundary; homopolymer self-loop (3 walks, never the record's); cycle junction; circle (`rejoined_seed` on both arms); repeat in two records plus a chimera; seed twice in one record (the walk runs through the k−1 junction k-mers before `rejoined_seed`); seed spanning a bubble (the other label is dropped, nothing else changes); an unlabeled region (recorded empty, `label_lost` at its edge); a reverse-complement record (not a carrier in basic, the second fork branch in canonical/primary); hairpin skip and follow; even-k palindromic node; 100 labels (default cap cuts the lists and the oracle says so; uncapped the contract holds); unmasked DBGSuccinct ('$' never recorded); 16 random fixtures at k = 7 | Cycle / self-loop, Reverse complements |
| T33 | switching vs the §6.3 recurrence | `OneSwitch`: forbid ends A in Y; budget 1 spells P·Y·Z under B at loss 1 switched at |P·Y|, budget 0 reports `loss_budget` needed 1; `TwoSwitches`: budget 1 cuts A's walk at |P·Y·Z| needing 2, budget 2 spells P·Y·Z·W under C at loss 2, B's walk needs one; leaves, losses and switch events equal the recurrence over the records | Sample-switching chain |
| T34 | trace vs the records | `TrieCasesTrace`: k-mer support follows the four combinations of two records sharing a stretch, trace the two records; a record ending where another goes on is `record_end` under trace and continues under k-mer; a seed twice in one record gives two traces; the homopolymer limitation (§6.6) | Cross-record boundary |
| T35 | size caps vs the completeness guarantee | `TrieCasesCaps`: `max_steps` 1…120, `max_live_paths` 1–3 (stop and beam), `max_output_bp` 1…90, `max_paths` 1–2 on the nested bubbles, three modes: every tuned run is a prefix-subset with a reason per omission, its walks equal the exhaustive trie's at every depth ≤ `complete_to_bp`, and `complete_to_bp` never decreases as `max_steps` grows | Caps |
| T36 | the contract on a real index | `MiniRefSeq.TrieContractAgainstTheSourceRecords`: the string model from the 42 source records (k = 31, basic, unmasked, RowDiff<BRWT> + coordinates, header labels); the whole blaNDM seed (19 carriers; 24 structural walks left, 2 right at 300 bp) and three 150 bp windows: records ⇔ T ⇔ A, trace ⇔ records, no '$' in any output; at 1000 bp the 19 carriers' claims under k-mer and trace support equal the records' | — |
| T37 | stated limitations and outcome (§7.0) | `Walker.LimitationsStateExactlyWhatLimitedTheResult` and `test_traverse_states_every_limitation`: each kind produced by its knob and absent without it, the four `outcome` axes per case; `WalkerTest.SwitchSourcesCutWithoutALabelEnd`: a `table` cost whose cut source goes on along another successor — no `switch_sources` end, E entered at 0.8 instead of 0.5, `switch_sources_cut` 1 and the limitation stated; `"unlimited"`: 0.5, nothing stated; `Walker.UnimplementedWindowsAndShortContinuationsAreRefused`: a non-zero tip / bubble window and `continuation_bp` 1 … k − 1 are 400s naming the field (and k), 0 is accepted; clamped integer knobs serialize as integers | — |
| T38 | the graphlet round trip (`DESIGN-traverse-graphlet.md` §9 T37) | `TestTraverseGraphlet.test_t37_round_trip` (CLI, 14 strategies: merge, keep, annotate exhaustive and beam, quorum, left only, `max_steps`, `max_steps: 1` in annotate mode (an arm with no growth bin), `max_labels_per_node: 1`, `sequences: false`, a radius below k, trace with column labels, a switch at a split, and a 3-seed batch with a failed derivation and a duplicate) and `test_api_graphlet` (HTTP, gzip): `detail: full` and `detail: graphlet` of one request; `from_response(r, out).to_json()` equals the full result after the §2.5 normalisation, `parse(text).dump() == text`, `is_canonical`, the `A` counts and the summary's counts describe the body, `graphlet_lines == Z`, `check_rules()` empty, label ends ≤ runs, every `Lw` run continued by a switched run at its `to_bp`. `test_t37_cli_fixtures_are_reproduced`: the whole-document fixtures are what this binary writes for their requests on their tiny indexes, byte for byte, and each still shows its counterexample (§9 v2 list). `GraphletCodec.GoldenVectors` (C++) and `test_traverse_codec` (Python): the 10,433 shared codec vectors | — |
| T39 | orientation (design T38) | `spell(left) + seed + spell(right)` is each record (acc1/acc2/acc3); the library's continuations equal the server's on both arms, contained (40 of 60 bp) and crossing into the seed (k of 10 bp); a continuation the library derives (`next_request`) is accepted and extends the walk | — |
| T40 | the oracle and the identity through the library (design T39) | `constrain.compare(annotate, mode='claims')` equal under the permitted labels and `only_in_b` = acc2 over all; exhaustive vs annotate equal in `claims`, `walks`, `labels`; tuned vs exhaustive `prefix_subset` with acc2's omission at 0 for `minority`; §3.1: two indexes over swapped memberships share `index_meta_fp`, differ in `index_fp` and never compare equal, without a manifest `unverifiable` | — |
| T41 | the graphlet protocol (design T40) | summary key sets; `output.detail/timing` echoed and resubmittable; `graphlet_format`, `detail_levels` and the index identity in capabilities and `H`; `sequences: false` (no bases, first bases on split children, continuation sequences, `splits[].char` reproduced); two runs byte-equal; body < ½ of the compact full JSON | — |

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
   events and the `tips` / `bubbles` growth counters). The fields exist and stay zero; a non-zero window is
   rejected (the JSON API with a 400 naming the field, `validate_strategy` for the C++ API) rather than accepted
   and walked as 0; T28 is deferred.
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
20. **Run anchors and terminal values** (`LabelRun::segment`, `branches`, `loss`; `DESIGN-traverse-graphlet.md`
    §4). The walker records on every run the segment on which it ended or was closed and the lineage's branch
    count and loss at that point, set where the run ends (`end_run`, from the ending entry and the head's
    segment) and where a merge closes it (on the parent the closed entry came through). Neither the events nor
    the paths determine them: clones made at a split are identical rows, a merge of three parents closes several
    runs at one depth, and branch counts are taken before quorum filtering, so an interior end's count can
    depend on a successor that left no trace. **A run can end without a `label_end` event:** a switch source
    whose lineage continues only under other names ends `label_lost` silently (the `switch` events at its
    `to_bp` document it, on the children when it ends at a split), so "every ended run has a `label_end`
    event" is false while "every `label_end` event belongs to exactly one run" holds.
    `Trie.EveryRunIsAnchoredWhereItEnds` checks all three kinds of end on the trie fixtures and dense
    switching graphs.
21. **Orientation of the graphlet** (§7.5.4): bases are stored in walking order on both arms, so positions index
    them directly; the JSON keeps the natural orientation.

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
