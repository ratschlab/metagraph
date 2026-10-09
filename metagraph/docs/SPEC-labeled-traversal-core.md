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
  "discover": {"max_labels": 1000, "kind": "header"},
  "support": "kmer",
  "min_block_kmers": 1,
  "select": {"policy": "max_support", "max_seeds": 10, "min_block_bp": 31,
             "max_labels_per_seed": 1000, "label_order": "hash", "sample_seed": 1,
             "merge_overlapping": true},
  "run_format": "intervals",
  "bounds": {"max_query_bp": 100000}
}
```

- `bounds.time_budget_ms` *(milestone 1b of `DESIGN-pattern-search.md`; a 400 before, "no deadline is
  enforced")*: optional, the request's deadline in milliseconds, with the semantics of §4.5 — the work stops
  at the budget less a finalisation reserve and the answer is then exactly the resolve of a query prefix,
  with a `stop` block; without the field nothing changes. `bounds.max_query_bp` bounds the work too: every
  in-graph k-mer's annotation row is read, so the k-mers read are the work. `discover.max_labels` caps the
  labels profiled and returned, not the rows read (a discovery reads every k-mer's whole row whatever it
  keeps).
- Exactly one of `labels` (explicit, ≥ 1; column or header names) or `discover` must be present
  (both is a 400): the example discovers; with explicit labels it would carry
  `"labels": ["573", "NZ_STEQ01000045.1"]` in place of `discover`.
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
  dropped_full_length}` is reported. *(Efficiency pass, the same answers byte for byte: the rows are decoded
  once — the support profile took its hits from the discovery's rows (tuple rows when the profile needs
  coordinates: header labels or `trace`; their columns are those of the whole rows) instead of decoding every row
  again, which made a column discover cost 1.6–2.1 times the profile of explicit labels on refseq33m —, the
  sequences' labels go into one flat hash table instead of a `std::map` node per label (0.8–1.7 µs a pair, 1.83 s
  for a header discover on 23S), and a header label's coordinates are mapped by their runs (§8.2), each sequence
  of a column counted once per k-mer. Feature level 5: the discovery accumulates the profiles while it reads the
  rows, below.)*
- **Seed candidates:** every maximal support run of each label with length ≥ `min_block_kmers`; labels with an
  identical run `[a, b)` are grouped. Candidate = `{kmer_interval, bp_interval, labels}`. Order: longer first,
  then more labels, then smaller `a`. Disconnected supported blocks of one label are separate candidates.

**A `/resolve`'s memory** *(review of pass 5, part B, and of its fixes, finding 7; feature level 5)*. A `/resolve`
decodes the present k-mers' rows in **batches** and holds one batch at a time: the first of 64 rows, each next one
sized from the widest row of the one before to about 64 MiB (`row_copy_bytes`), at most twice the previous
batch's rows and at most 4,096. A **discovery** accumulates, per label, while it reads the rows: its k-mers (the
ranking) and its profile — support runs; for `trace` the chain of coordinates and its breaks, k-mer by k-mer —
so the top labels' profiles need no second read. It keeps the row of a k-mer that occurs again later in the
query for that occurrence, at most 256 MiB of such rows (`row_copy_bytes`), and decodes a repeated row again
beyond that. An **explicit** profile primes its label query with each distinct present row once, batch by batch
(the hits are kept per k-mer, for its labels only). What a request holds is therefore one batch, the rows kept
for repeats, and per label found its accumulator (counts, runs, trace breaks, the chain's current coordinates);
a discovery on column labels indexes its accumulators by a slot of 4 bytes per column of the index. Before, the
efficiency pass decoded and held every present row until the profile was done, and the first fix of this review
kept a prefix of 256 MiB for the profile pass, which decoded every other row again in one call — the same peak,
and 40–56% more CPU on wide rows. Measured on a synthetic coordinate index of 8,000 labels whose core rows hold
24,000 coordinates (user CPU, medians of 5, interleaved, on a machine loaded by other work — about ±10% noise; the
previous build → now): on queries of 0.3–3 kb of that locus column presence, column and header trace and explicit
labels 3–17% faster and header presence from 7% faster to 14% slower, at 144–167 MB instead of 0.23–1.66 GB; on a 9.8 kb
genome repeating the locus three times column presence faster (1.20 → 0.91 s) and explicit labels within 8%, at
165–236 MB instead of 1.7–1.8 GB, and header and trace discovery 24% and 48% slower (5.72 → 7.10 s, 3.04 →
4.50 s) at about 370 MB instead of 4.8 GB — its repeats beyond the 256 MiB kept are decoded again. Every answer
was the same. The ranking (more k-mers first, ties by column, then sequence) and the profiles do not depend on
the batches or on what is kept (tested with batches of 1 to 4,096 rows, byte targets of 1 B to 64 MiB, 0 to
256 MiB kept, on a query and on one repeating two thirds of it, column and header labels, presence and trace;
a discovery's profiles equal those of the same labels given explicitly).

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
to a header). A request with `bounds.time_budget_ms` also carries `limits` and `stop` (§4.5); a request without
it carries neither.

### 4.5 The deadline (`bounds.time_budget_ms`)

*(Milestone 1b of `DESIGN-pattern-search.md`: the search service runs `/resolve` as a job and asked for a
deadline. `src/cli/traverse.cpp` `process_resolve_request`, `src/graph/traversal/resolve.cpp`.)*

- **Opt-in.** A request without `bounds.time_budget_ms` runs without a deadline and is answered byte for byte as
  before — no `limits`, no `stop`, no deadline set or read (the clock of `timing.elapsed_ms` is read as it
  always was), the explicit labels' hits one fetch —, whatever the server's cap. With it: a number of
  milliseconds **above the finalisation reserve** (250 ms, the same as `--pattern-finalize-ms`'s default, fixed
  in this build: no flag sets it; a budget not above it, a negative number or a non-number is a 400 naming the
  field), lowered to the server's `--traverse-max-time-ms` (the cap of `/traverse`'s budget, default 30 000) and
  never above 899 000 ms on the server, also when the flag is 0 or larger (the 900 s content timeout less 1 s for
  the transport; review of 2026-10-07, R1-05); the CLI applies none; a lowered budget is stated in
  `limits.clamped: [{"field": "bounds.time_budget_ms", "requested": …, "effective": …}]`. The deadline starts
  when the request's body has been parsed, as `/pattern`'s (`DESIGN-pattern-search.md` §5.3). A cap at or below
  the reserve would leave no time to work: every budgeted request would stop before its first row.
- **The work stops at the budget less the reserve.** The deadline is read where the client's connection is
  checked during the work, after it: before the first row is read; between two batches of rows (a discovery's
  pass and the explicit labels' priming: 64 rows, then up to 4,096 rows or about 64 MiB of rows, §4.2); and,
  on a direct-access annotation, before every 4,096 k-mers of the explicit labels' hits, which under a deadline
  are fetched in pieces of 4,096 k-mers (one fetch of the whole query, the path without a deadline, is a piece
  no clock read can end: it reads every k-mer's cell of every label). On the row paths the hits of the k-mers
  each priming batch completes (those before the next batch's first k-mer) are taken from the primed rows right
  after that batch, before the deadline is read, so a stop between batches keeps every row read (review of
  2026-10-07, T3-01/V1-01: a separate pass read the deadline again and cut such answers to about 4,096
  k-mers). A key's hits do not depend on the piece they are fetched in, so the profile is the one fetch's.
- **A stop answers the resolve of a query prefix** (`DESIGN-traverse-graphlet.md` §13), status
  200: exactly the profile of the sequence's first x k-mers — `num_kmers` = x, `graph_runs` cut at x, the
  labels with their support in the prefix (a discovery's labels are those met in it, ranked and truncated
  there: `labels_truncated` is the prefix's), the candidates and the selection made from them; no field
  describes the rest of the query — and the `stop` block, beside `limits`:

  ```json
  "stop": {"phase": "rows", "reason": "time", "resolved_kmers": 1984, "query_kmers": 2970,
           "resolved_bp": 2014, "query_bp": 3000, "remainder_from_bp": 1984,
           "message": "bounds.time_budget_ms 1000 ms less the finalisation reserve of 250 ms ran out while …"},
  "limits": {"time_budget_ms": 1000, "finalize_reserve_ms": 250, "clamped": []}
  ```

  x (`resolved_kmers`) is the first k-mer in the graph whose labels were not read (the k-mers before it that are
  absent from the graph are resolved: the mapping told); `resolved_bp` = x + k − 1 (0 for x = 0);
  `remainder_from_bp` = x: the sequence from that base holds exactly the k-mers not resolved, so a client
  resolves the rest with it. A run ending at x may continue past it (in the remainder's answer it starts at
  its k-mer 0). `phase` is `rows` (the rows of a discovery or of the explicit labels' priming were being read,
  or none yet) or `support` (a direct-access annotation's explicit labels' hits, or a stop before any).
  An explicit `select` whose interval ends past x, or, on a discovery, names a label the prefix did not
  discover, is **not made**: `selection: null`, the message says which seed and why (whether it holds on the
  whole query is unknown, and a 400 would blame the request for the deadline). What holds whatever the stop is a
  400 as without one, every seed checked before any selection is left unmade (review of 2026-10-07, V1-02): an
  interval empty, past the query's k-mers or not fully in the graph (the whole query's k-mers are mapped before
  the deadline is first read), and, with explicit labels (each of them profiled whatever the stop), a seed label
  that is not one of them. A request whose work completes answers `stop: null`; `limits` (the effective budget,
  the reserve, the clamps) is stated either way.
- **The answer is built within the reserve.** After the work, the loops over the labels (a discovery's ranking
  and naming, the profiles, the candidates' grouping) read the whole budget every 4,096 labels, the selection
  once it is done, the JSON every 4,096 labels, candidates or seeds, its text every 64 KiB and its compression
  every block (the server's transport check, as `/pattern`'s). Past the budget the answer is **503**
  `{"error": "resolve: the answer could not be built and written within bounds.time_budget_ms (… ms, the
  finalisation reserve of 250 ms included): nothing partial is sent", "code": "deadline"}`, never a partial
  answer; the CLI writes that body and exits 1. (A 503 with `Retry-After` and no `code` is still the loading
  index's.)
- **Not polled** — each a piece the deadline can fall into, so a stop comes up to one of them after the work
  deadline, and their overrun of the reserve is a 503: the request's parse; the mapping of the query's k-mers
  (one call, linear in the query, which `bounds.max_query_bp` and `--resolve-max-query-bp` bound); the
  explicit labels' resolution and their query's setup (linear in their bytes, and no server limit caps their
  number; an unknown label is refused before the deadline is first read, so a stop never hides one); one row
  batch's decode (a row-diff row's whole path) and its accumulation (every label of its rows); a discovery
  pass's setup (the next occurrence of every present k-mer, one hash pass) and the explicit priming's distinct
  keys; one piece of hits and its accumulation; the ranking sort of a discovery's labels and the candidates'
  sort (n log n); `select_seeds` (a `max_support` round sweeps every run); and the transport. The client's checks
  are where they were (§10.3); a reserve that these overrun is not adapted to them (the reserve is fixed).
- **Capabilities.** Both GET routes carry the `resolve` block (§10.3). The field is not a `feature_level`: every
  `/resolve` and `/traverse` response states the level, so a bump would change every answer; a client gates
  on the block. An older server refuses the field (400 naming it).
- **Not in this milestone** (stage 3b-1 of `PLAN-traverse-next-stages.md`): `/resolve`'s memory and work
  budgets, a server default deadline for a request without one (`--resolve-max-time-ms`), reads chunked inside
  a batch (`DecodeControl`), and stage 3b's `outcome` / `resource_stop` vocabulary.

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
  `detail: full` spells every path's chain); on an annotation whose reads are budget-aware (§6.8, stage 3 of
  `DESIGN-traverse-graphlet.md` §14.1) also every annotation read, charged inside the decoder. Work is charged
  work units (§6.8). Either budget stops the whole
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
  **Delivery is charged at demonstrated upper bounds, not averages**: every object at the worst case of what the
  serialisers hold of it at once — the JSON tree (jsoncpp's nodes, key copies, maps and string buffers) with
  three copies of its text (`Json::writeString` holds its stream's buffer, grown by doubling, and the copy it
  returns), or a graphlet's body four times (the writer's text, grown by doubling, and its JSON copy; then that
  copy with the escaped body three times), each field at its widest — plus zlib's deflate state, the envelope
  and the largest number of limitations a result can state (5 at the seed level, 8 per arm). Names (label
  names, dropped labels' names, the `seed_id`) are the one input whose delivered size is no fixed multiple of
  its length: JSON writes a control character as six bytes and a non-ASCII character as six or twelve, MGT a
  `%` as three, and a graphlet's body is escaped once more inside its JSON string; so each name is charged at
  its escaped length in every place the detail writes it (a seed label twice in `full`: `label_dict` and
  `seed.labels`). Numbers are the other: canonical MGT writes a float positionally (§7.5.1), so `1e-300` is 302
  characters; every MGT float — a cost, a loss, a needed budget, a time — is charged at the widest any value the
  request allows can be written (`mgt_float_width`: the smallest positive cost bounds the leading zeros, the
  loss budget plus the largest finite cost the integer digits, and the time budgets, requested and effective,
  theirs; an elapsed time is measured in clock ticks of 1 ns and is at least the budget it tripped), 24
  characters for every usual request (review of the stage-2 recheck, P1: switches of 1e-300 delivered 22 MB
  within an account of 16 MiB). `Graphlet.DeliveryCostsBoundTheOutput` checks the bounds against the
  serialisers' output, adversarial names and wide floats included, and the widths they assume (keys, nesting,
  numbers, effects). The request's own
  echo (`strategy`, with `labels.extra` and a cost table) is outside the per-seed budgets: its size is the
  request's.
  **A failed or refused seed's result is not admitted** (no walk exists to admit): it echoes the request's
  `seed_id`, priced as a delivered result prices it (the fixed part and the `seed_id` in every copy the
  serialisers hold), and whatever exceeds the budget is stated as `memory_bound_soft`'s observed excess (§7.0),
  never left unstated (a `seed_id` the depth-0 admission refused still comes back whole). A name the INDEX
  supplies that a failure echoes — the header of an `ambiguous_header` derivation, in `error` and as
  `observed`, and a column it collides with — is, under `bounds.max_memory_mb`, echoed whole only up to 256
  bytes; a longer one as its first 256 bytes (never inside a UTF-8 sequence), `...`, and `(<n> bytes; column
  <c>, sequence <s>)`, so that nothing the request does not bound can make a failed result larger than the
  budget. Without a budget it is echoed whole, as before.
  **A memory bound is not a wall-clock bound**, and the walk's deadline does not cover delivery: serialising
  and compressing an admitted result take time the deadline does not check (a hard overall deadline would need
  enforcement during delivery, a staged item, §6.8); what bounds the delivery's time is the size the memory
  budget admits — except for a request with `attempt_id`, whose duration bound the server enforces through the
  delivery too (§6.8, the attempt's bound).
- **Attempt fields** (`DESIGN-traverse-graphlet.md` §14 and §14.1: the frozen wire contract's request fields;
  stage 4, backend half — the ledger itself lives in the search service): top-level `attempt_id`, `budget_id`
  and `locus_id`, each optional, each a string matching `^[A-Za-z0-9._:-]{1,128}$`. A request with
  `attempt_id` is an **attempt** of a ledger that reserved an allowance for it: the server registers it while it
  runs; an id does not run twice at once on a server process, nor again while it is retained or held there — a
  second request with an id that is running, retained, held or tombstoned (§10.3) is refused with 409 and runs
  nothing; once its retention and hold are over (at once with `retention_s` 0) a request with the id runs again
  *(review of 2026-10-06, D3: "an id runs once per server process" was read as unconditional)* —; it can be
  cancelled by id (`POST /traverse/cancel`) and its
  state read (`GET /traverse/attempt/{attempt_id}`, also by a worker other than the one that sent it); every
  response to it carries a response-level `usage` block (§7.3), its 4xx/5xx answers after the request was read
  included; and the server enforces a duration bound on it (§6.8). `budget_id` and `locus_id` are echoed in
  `usage` and nothing else: no allowance is enforced by them. Either one without `attempt_id`, or a malformed
  id, is a 400 naming the field (without `usage`: nothing was registered). A request **without** `attempt_id`
  keeps its response byte for byte (no `usage`; the one change for it is that a client that is gone is not
  answered, §6.8). No request-level total budget exists (deferred by the owner). The CLI accepts the fields and
  states `usage`, with `bound.enforced: false` (operator-run: nothing to lease), its error output for a request
  error after the fields were read included (`{"error", "usage"}`, as the server's 400; a malformed id: no
  `usage`). `/resolve` does not accept them (400, unknown field).
- **`not_after_ms`** (pass 5; top level, optional, with or without `attempt_id`): an integer in
  [0, 2⁵³ − 1], a Unix epoch instant in milliseconds on the server's clock, after which the request
  must not be **started**. A ledger that has no answer for an attempt treats it as one that **cannot start
  subsequently** once its own clock passes `not_after_ms + attempts.clock_skew_allowance_ms` (§10.3) — an
  unanswered request may already be running, which the attempt's bound covers *(review of pass 5: "never
  started" claimed more than the check gives)* —, and as past its bound once it passes that plus the attempt's
  `bound_ms`, **apart from what it runs past its bound up to its next delivery check, a run whose length has no
  stated bound** (§6.8, the attempt's bound and the stated limits; §10.3, the release rule's clock release)
  *(review of 2026-10-06, X2: "as stopped" left that run out; the review of its fixes: it is not one
  uninterruptible step, since only the delivery checks compare the bound)*. The
  server makes the first of these safe — the start, not the run — by refusing, when the
  request's handler starts, a request whose `not_after_ms` is earlier than its clock: **409**
  `{error, state: "expired", not_after_ms, server_time_ms, attempt_id?, budget_id?, locus_id?,
  server_instance}` (the ids as given), nothing run and nothing registered — no `usage`, and a later
  `GET /traverse/attempt/{id}` answers 404 (a later copy of the request is expired too). The check is
  strict: no allowance is added (the skew is the ledger's to add), and a request whose instant has
  not passed runs whatever happens afterwards (it is a start deadline, not a run deadline: the
  attempt's bound limits the run). The refusals are made in one order, under one lock: an
  `expect_server_instance` naming another process first (`instance_mismatch`), then an id that is
  running, retained, held or tombstoned (the duplicate's 409, which carries `attempt`, whatever the
  request's `not_after_ms`; a held or tombstoned id's refusal extends that hold, §10.3), the expiry
  last — so an `expired` 409 always means that no attempt with that id existed on that
  `server_instance` when the request was judged, and no copy can start there later while its clock
  does not step back below `not_after_ms` (a release ground, §10.3). The capabilities' `not_after`
  text states this order *(review of 2026-10-06, C24)*. A fraction, a sign, a string or a value above
  2⁵³ − 1 is a 400 naming the field, without `usage`. It is echoed as `usage.not_after_ms` and in the
  attempt's state, only when given. The CLI applies the same check (the 409 body on stdout, exit
  status 1). `/resolve` does not accept it (400, unknown field).
- **`expect_server_instance`** *(review of pass 5; with `attempt_id` only)*: the `server_instance` (§10.3) the
  request is meant for, a string matching the ids' pattern `^[A-Za-z0-9._:-]{1,128}$` (as every
  `server_instance` does), compared as given. Any other value — null, a non-string, `""`, more than 128
  characters, another character — is a 400 naming the field, without `usage`: nothing runs and nothing is
  registered (refusing `""` keeps a request from being un-pinned; the capabilities' `instance` text states it,
  *review of 2026-10-06, X4*). A process whose own `server_instance` is another refuses it
  before anything runs or is registered: **409** `{error, state: "instance_mismatch", expect_server_instance,
  server_instance, attempt_id, budget_id?, locus_id?, not_after_ms?}` (checked first, ahead of a duplicate id and
  of `not_after_ms`). Tombstones (§10.3, `POST /traverse/cancel`) and the holds of finished attempts live in
  memory, so a restarted process — a new `server_instance` — holds none, and a delayed copy of a cancelled or of a
  finished request would otherwise run there; a ledger that releases capacity early on a tombstone, or on a
  finished state before `not_after_ms + clock_skew_allowance_ms + bound_ms` (§10.3, the release rule), sends its
  attempts with this field *(the finished state: review of levels 4–5, finding 6)*. Without
  `attempt_id` it is a 400 naming the field. The CLI, whose instance is its own and random, refuses any
  well-formed value (the 409 body on stdout, exit status 1; a malformed one gets the 400's error output, exit
  status 1). `/resolve` does not accept it (400, unknown field), nor does a server below feature level 5 (400,
  unknown field).
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
- `labels.extra`: additional permitted labels (switch targets). Rejected when unreachable: model `forbid`, or no
  **chain** of switches from a seed label (through the request's labels: seed and extra) enters it with a summed
  cost within `loss_budget` — the cheapest chain, found by a shortest path over the cost model (a `table`'s
  default relaxed lazily, so a pool of thousands of extra labels validates in O((labels + entries) log labels)).
  The walk enforces the cumulative loss switch by switch, so a label two switches away is a valid target
  (review of the stage-2 recheck, design answer 2: `A → B = 1`, `B → C = 1` under a budget of 2 rejected `C`).
  The rejection names every unreachable extra label (the first eight, and how many more). `extra` applies on
  top of a derived set, and an `extra` label the derivation also produced is rejected as a duplicate seed
  label.
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
- **`output.coordinates`** (boolean, default `false`) and **`output.max_coordinate_occurrences`** (an integer in
  [1, 2^64 − 2] or `"unlimited"`, default 16) *(feature level 6; opt-in, `DESIGN-traverse-graphlet.md` §18)*:
  report where each run's sequence lies (§7.1, the `coordinates` block). The cap **without**
  `coordinates: true` is refused (400 "only with output.coordinates true, where it bounds the occurrence
  lists; …"): it would change nothing (§7.0). 0 or a number above the maximum (`out of range [1,
  18446744073709551614] (or "unlimited")`), another string (`expected an integer or "unlimited"`), a negative or
  non-integral number (`expected a non-negative integer`) and a non-boolean `coordinates` (`expected a boolean`)
  are refused naming the field. Both are accepted and inert under `support: kmer` and in `annotate` mode, where
  the result states why it has none (`coordinates_reason`), so that a request switching its support stays valid.
  The echo carries both fields — the cap's default included — only when `coordinates` is true, so it stays
  resubmittable; `coordinates: false` is byte-identical to omitting it, and a request without coordinates is
  answered as before (feature level 5) apart from the `feature_level` digit.
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
   with its `graph_runs` (a partially present sequence is never treated as one seed). A request that selected its
   graph with `graphs` (§10.3, multi-graph mode: the fan-out of `/search`, which sends each seed to graphs that
   need not hold it) gets such a seed as its result instead — `{"seed": {…}, "error": "Seed is not fully present
   in the graph; graph runs: …", "not_in_graph": {"kmers": n, "kmers_present": p}, "limitations": […],
   "outcome": {"walks": "not_in_graph", …}}`, no arms, no graphlet, `p` the seed's k-mers this graph has (0: none;
   a partially present seed is answered so too) — and the other seeds are traversed; `limitations` is empty but
   for `memory_bound_soft` under a memory budget, and a request with `output.coordinates` states `coordinates:
   null` with its reason, as for a failed seed.
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
   incrementally and the seed is rejected as soon as it is empty. Rows are reconstructed in windows of 64 k-mers
   (the last one shorter; the batch is what makes a long seed fast), so an empty intersection stops the scan within
   one window rather than at the exact k-mer; `rows_requested` counts the rows actually CONSUMED, which is what the
   derivation read, not what the batch fetched.
   The pass starts from the **cheapest row among the seed's first 64 k-mers**, not from k-mer 0: that one row is
   the only one no intersection has narrowed yet, so starting at k-mer 0 would make the peak cost an accident of
   where the caller cut the seed. **Every window is 64 k-mers whatever `annotation.batch_kmers` is**: the
   cheapest-row choice and the guard below decide whether a seed is accepted, and each window is one work charge
   before its k-mers are consumed (§6.8), so its width also places the work comparisons, `largest_charge`, the
   outcome under a work budget and, after an early `no_carrier`, `work_seed` — none of which may depend on a
   fetch-size setting (§6.8). *(Feature level 6, the review of 2026-10-06, W1: on a format whose reads are not
   budget-aware the later windows followed the knob, `clamp(batch_kmers, 1, 64)` k-mers wide, and a work-budgeted
   derived seed of 195 k-mers over 992 columns walked partial at `batch_kmers` 1 and 7 and failed at 64; an
   unbudgeted `no_carrier` at k-mer 90 charged 639 units at 1 and 1,152 at 64. The fixed width costs at most 63
   rows read past an early `no_carrier`.)* The order does not change the result (an
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
   ran out before the derivation consumed its first k-mer. A time budget that runs out after 1 ≤ j < n of the
   seed's n k-mers does **not** fail the seed (D3, feature level 6; *review of 2026-10-06, X1*: this sentence
   listed every mid-derivation timeout as a failure): the seed is walked with the set of the labels carrying the
   j k-mers read — a superset of the whole seed's carriers — truncated at `complete_to_bp` 0, `walks: partial`,
   `label_evidence: qualified`, a `derivation` limitation first with `observed` = j (§7.0, the `derivation` row),
   unless that superset would fail the seed in another way, which then fails it with the time budget as before
   D3. One failure of the same shape is not a derivation's: a seed whose labels — the
   `label_dict` of either mode (an `annotate` dictionary is complete only after the walk) or a `dropped_labels`
   entry — include a name that is not valid UTF-8 (names come from FASTA headers and file names) is **refused**
   with cause `unrepresentable_label_name`, whether or not its labels were named. No output carries such a name
   verbatim, and a replaced one (U+FFFD, which earlier servers wrote) can be another label's valid name anywhere in
   the index: resubmitted as a continuation's seed label it resolved to that label and walked other bases. The
   refusal names the label by its column (and sequence id), never by its bytes. Those do not fail the request: the
   seed's entry in `results` is
   `{"seed": {"seed_id", "length_bp", "labels_from_seed": true}, "outcome": {"walks": "failed", …}, "error":
   "<reason>", "limitations": [{"kind": "derivation", "cause", "knob", …}]}` — no `arms` — and the other seeds are
   traversed as usual. The `derivation` entry names the request field that would get past the cause (§7.0); a
   client distinguishes the two shapes by `outcome.walks` (or the presence of `error`). Everything that
   *is* the caller's own doing still fails the whole request with HTTP 400: a malformed seed (too short, invalid
   character, not fully present in the graph — but with `graphs`, step 2), an unknown or duplicate **explicit**
   label, and every strategy error.
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
  model is deterministic (fixed bytes per object, never a measurement), so a memory stop, in either phase
  (`traversal`, `annotation_decode`), is reproducible and independent of `annotation.batch_kmers`. Work units
  count **deterministic logical work, not measured decode effort** (review of stage 3, answer 1: the physical
  decode counters — `rows_fetched`, `tuple_rows_fetched`, `coords_mapped`, `annotation_fetch_ms` — are in
  `timing`, apart): a weighted sum of what the walk did — successor enumerations (4), the annotation rows its fetches return
  (8 per key, 1 per entry and, under `trace`, 1 per coordinate: a recorded row its whole width; on a
  budget-aware annotation each row also its row-diff dependency rows, below), pair evaluations, refusal scans,
  edge-reuse probes, derivation scans and steps (1 each). A row is charged **when a fetch returns it** (the design's "rows decoded"), whether
  or not a head then consumes it: a level cut mid-way decoded its later rows all the same, and charging only
  consumed rows hid that decoding from the budget and from `used`. A row the lookahead decoded ahead is charged
  when a level's fetch returns it from the cache; a row the lookahead decoded that no fetch asks for (at most
  `min(batch_kmers, the radius left)` per head of a level, §8.3) is not charged; a row evicted and fetched again is
  charged again. The row-diff path cache (§8.4) changes none of this: a row is charged its whole row-diff path
  whatever the decode shared or cached.
  **The work bound is stated with its number, not as a fixed maximum**: the walk compares the work budget after
  every charge (each fetch call, a derivation's scan, each priced target, each edge-reuse probe, each
  enumeration); the previous comparison passed, so a stop exceeds the budget by at most what was charged since,
  and every work stop's message states the most its seed charged between two comparisons
  (`ResourceAccount::largest_charge`), whatever the charge was: a fetch call's rows, decoded whole with their
  coordinates (on a budget-aware annotation each with its row-diff dependency rows: near the budget a call reads
  one key, so the overrun is then one row with its whole dependency path), a node's
  label-state scan (|σ| + |A(v)|), or the seed phase with the roots' rows (both arms' roots belong to the
  result complete to 0 bp, so they are charged before the first comparison). It can exceed `W` on a wide index
  (70,000 records on one k-mer: one row of 70,008 units; with `direction: both`, both roots: 140,016): only the
  seed phase can fail, at a comparison finding it at least `W` units over budget (`W` is that threshold, not a
  ceiling on how far it ran past the budget), while the roots' rows are one charge as wide as the index makes
  them. Under a
  request budget a level's annotation is fetched in calls of at most 8,192 keys, sized so that a call holds about
  `W` units at the widest row fetched so far and, under a work budget, no more than the budget has left (down to
  one key near the budget), growing from one key at the start of each level so that rows wider than any fetched
  before are met by a small call; without a request budget the level is one logical call, as before (its counters
  depend on the batching). The deadline is read before every head, between those calls, before every time-sized
  chunk an annotation read is decoded in (*chunked deadlines*, below), and at least every `W` = 65 536 units
  otherwise. The seed phase is charged too (8 per seed k-mer read, 1 per annotation entry and coordinate, as the
  walk charges a fetched row; `account.work_seed` beside each arm's `work_units`), **every row a read returned
  charged before a comparison can fail the seed** (review of the stage-2 recheck, P2: charged row by row, a
  failure left the rest of a decoded batch uncharged): the validation reads, under a work budget, in calls grown
  from one key, doubling, each holding about `W` units at the widest row read so far (its hits, coordinates and,
  budget-aware, dependency rows), each call one charge; the derivation's window (64 k-mers whatever
  `annotation.batch_kmers` is, §6.1,
  held at once to choose its cheapest row) is one charge as soon as it is read — a window found `too_wide`
  is added to the seed's work without a comparison, so that `too_wide` stays the stated cause and the work it
  decoded is still the seed's (`work_seed`, the usage's `work_units`; review of the stage-3 fixes, P3: such a
  seed reported 0 units for 400,019 decoded). It is compared once every `W`
  units of its own work: a comparison that finds it past the budget by less than an interval lets it go on, so
  that once it ends the stop at the first head delivers a valid result complete to 0 bp (what it charged since
  the last comparison that passed, with the roots' rows charged after it, is the one further charge stated
  above, `largest_charge`); a comparison that finds it an interval or more past the budget fails the seed (so a
  failed seed phase ran past the budget by less than two intervals plus one fetch call or window; the failure's
  message states the most its seed charged between two comparisons, as a walk's work stop does). The depth-0 state is admitted like a head: a memory budget that does not hold it
  fails the seed (what its dictionary held is stated, §7.0 `memory_bound_soft`). Its account includes the
  caches' allotments of the budget (a quarter of it, at most 64 MiB, for the label cache and a sixteenth, at
  most 64 MiB, for the lookahead), which grow with the budget, so the need at this budget is not the knob value
  that holds it: a memory stop states, as its `walk_domain` observed and in a depth-0 failure's message, the
  smallest budget (whole MiB) whose account holds the need with that budget's own allotments
  (`memory_budget_holding`; review of the stage-3 fixes, P2: "need 7 MiB" at 1 MiB, then 9 at 7 and 10 at 8;
  10 held it). For an exact need, raising the budget to it holds that state; a lower bound stays "at least".
  Every counter
  belongs to the arm whose head did the work. The deadline is never checked at depth 0 (`time_budget_ms: 0`
  still means one level, then the boundary check), and **a level at the radius is never stopped**: its heads
  only end (`max_extension_bp`), within what they reserved, so neither a budget nor the deadline trips there,
  and a deadline that passes when every remaining head has reached the radius censors nothing and states no
  stop (no `resource_stop`, no `Q`: the walk is complete). A level with no annotation key to read (the radius,
  or heads that are all dead ends) is not admitted as a read either: it finishes within its heads'
  reservations (review of stage 3, F3: a completed radius-0 walk whose depth-0 state filled the budget exactly
  was stopped by its empty lists). **Not implemented** (see `DESIGN-traverse-graphlet.md`
  §14): the per-seed delivery bounds `max_output_bytes` / `max_events` (§6.7), a request-level time budget
  (`request_time_budget_ms`) that would mark the seeds it leaves unstarted as `not_started` — every seed of a
  request is traversed, each under its own `time_budget_ms`, and the server bounds a request by
  `--traverse-max-seeds` (§10.3) — and a hard overall deadline for a request without `attempt_id` (the walk's
  deadline stops at the walk; §5, request budgets). For a request with `attempt_id` the server enforces the
  attempt's bound, through the delivery (below).
- **Stops from outside the walk** (stage 4 of `DESIGN-traverse-graphlet.md` §14, backend half; the core
  library's `AttemptControl`, which knows no HTTP). The walk polls its caller where it already reads its
  budgets and its deadline — before every head, at least every W charged units, at every charge of the seed
  phase (a validation's fetch call, a derivation's window) and per k-mer of a derivation — and the request polls
  between seeds. A poll is an atomic load (the stop flag); every `poll_stride`-th poll of a walk (8, stated as
  `deadline_check.poll_stride` on both capabilities routes), and every poll between seeds, also reads the clock
  (the bound) and, at most every 100 ms, the client's socket — so a cancel reaches the walk at its next
  checkpoint, the bound and a gone client within `poll_stride`. The poll before every chunk of a paced
  annotation read (*chunked deadlines*, below) reads the clock and the socket too, so a stop that falls inside
  a read the deadline splits is seen within one chunk; a read far from its deadline is one piece, and a
  cancel or a gone client is seen at the checkpoint after it.
  - **`cancelled`** (`POST /traverse/cancel`) and **`attempt_deadline`** (the attempt reached its bound,
    below), requests with `attempt_id` only: the seed being walked stops like a budget at the next checkpoint —
    its heads censored with `resource_limit` (`Y`), `complete_to_bp` the last complete level, a `resource_stop`
    with `scope: "attempt"`, `resource: "cancelled" | "attempt_deadline"`, `phase: "traversal"`, `requested`
    = `effective` the attempt's bound and `used` its elapsed time (whole ms, rounded up), `actions:
    ["continue_from_leaves"]`, and an arm `walk_domain` naming `attempt_id` (§7.0). In its seed phase the seed
    fails — no result exists yet, not even one complete to 0 bp: no arms, an `error`, a seed-level
    `walk_domain` naming `attempt_id` and `actions: ["retry_attempt"]`. The seeds not started yet are failed per
    seed with `resource_stop.phase: "not_started"` (§7.0). No budget ran out, and the statements say so
    (`time_budget` would tell a reader to raise one). A cancel that arrives after every seed was walked discards
    nothing: the response is written (within the bound) and `usage.reason` is `completed`. MGT v1 is unchanged:
    only free tokens take new values (`Q` scope `attempt`, resources `cancelled` and `attempt_deadline`, phase
    `not_started`, action `retry_attempt`; the `K` knob `attempt_id`); the amounts are integers, so no float is
    priced wider.
  - **The client is gone** (every `/traverse`, with or without `attempt_id`; `/resolve`, which **a client
    that goes away stops within one row batch**: it is checked at its phase boundaries, between the row
    batches of the discovery and of the explicit labels' priming — 64 rows, then up to 4,096 rows or about
    64 MiB of rows (§4.2) — and every 4,096 k-mers of an explicit-label support pass (*review of 2026-10-06,
    X-EFFICIENCY-04*), with or without `bounds.time_budget_ms`, each check a peek at the socket; the one
    exception is the direct cell path — at most 16 explicit column labels on an annotation with single-cell
    reads, none of the row-diff family — without `bounds.time_budget_ms`, whose read of the query's cells is
    one piece): the client closed or
    reset its connection, or half-closed it (a client that half-closes after its request is treated as gone, as
    nginx does by default). Found by a non-blocking peek at the request's socket that consumes nothing (at most
    every 100 ms); where bytes past the request wait on the socket (a pipelined request, a trailing CRLF), which
    hide the close behind them from a peek, the kernel's TCP state is read instead (`TCP_INFO` on Linux,
    `TCP_CONNECTION_INFO` on macOS: every state past `ESTABLISHED` is gone — the peer's FIN (`CLOSE_WAIT`) or
    reset (`CLOSED`), and the server's own shutdown at its content timeout (`FIN_WAIT1`/`FIN_WAIT2` and the
    states after them), in which nothing can be delivered either; `SYN_RECV`, a TCP Fast Open connection's
    state, which this server's listener does not enable, reads connected; review of the stage-4 backend, F3:
    such a client was walked to the end and answered into a dead socket) — on another platform such a client
    is seen only once its bytes are read, i.e. not during the walk. A descriptor that holds no connection (closed,
    or no socket) reads gone. The walk is abandoned where it is (nothing finalised, everything freed), nothing is
    written and the connection is closed. A walk that outlived the HTTP server's content timeout (900 s at the
    earliest — the timer runs when an io thread is free —, after which it shuts the connection) is stopped by
    this too, where it used to compute on: also with bytes waiting, which on Linux the server's own shutdown
    keeps, so that its state (`FIN_WAIT2`) read connected until the client closed its end *(review of
    2026-10-06, U13-02: only the peer's states were read)*.
  - **The attempt's bound**, enforced for a request with `attempt_id` (the server's operator bound, derived from
    the per-seed budgets — not a request knob): `bound_ms = min(n_seeds × T + allowance_ms, 899 000)`, T the
    effective per-seed `bounds.time_budget_ms` after the server's clamp (0 when not positive), `allowance_ms` the
    server's `--traverse-attempt-allowance-ms` (default 10 000) and 899 000 the content timeout less one second,
    on the attempt's clock, which starts when the server read the request's header (time queued before that,
    every server thread busy, is not in it: §10.3). The seeds stop being walked at
    `bound_ms − max(allowance_ms / 2, reserve_ms)` (`attempt_deadline`, as above; the seeds left are
    `not_started`) — the walk-until, which moves with the reserve; the walk stops at its first poll **that reads
    the clock** after it (one in `poll_stride` — 8, `deadline_check.poll_stride` — of a walk's polls, the one
    before a head; every poll before a paced read's chunk or in the lookahead; the poll before each seed), so up
    to `poll_stride` − 1 heads after it, while a cancel is seen at the next poll of any kind —
    leaving the rest to deliver what was walked; the response is built and written under the
    bound — checked every 4096 objects of the JSON tree and of the MGT text, every 64 KiB of the JSON text,
    between compression blocks and once before it is handed to the transport — and one not ready by the bound
    is not written: 503 with `usage` (reason `deadline`). **The bound itself is compared only at those delivery
    checks** (the walk's polls compare the walk-until, below it), so past the bound the attempt runs on until
    its next delivery check (then the 503) or its handler's return: the rest of the piece it was in when the
    bound passed and, when its walk had not stopped by then, the walk up to its next poll that reads the clock,
    the stopped seed's finalisation and the building of its result up to the first delivery check — several of
    the unpolled pieces of the stated limits below, a run whose length has no stated bound; a ledger's clock
    release takes the attempt as stopped at its bound and assumes that run ended (§10.3, the release rule)
    *(review of 2026-10-06, X2; the review of its fixes: "one uninterruptible step" understated it — on the
    reviewer's 6.94 Mbp seed the k-mer mapping, the seed phase up to its first poll and the failed result's
    building ran past the bound, and GET answered `running` 2,720 ms past the clock release's instant)*.
    **The delivery reserve** (pass 5): `reserve_ms = 1.25 × ((T + E) / (compress_mbps × 1000) + E / (build_mbps
    × 1000)) + stop_ms`, T the exact bytes of the text of the seeds finished so far (the server writes each seed's result
    as compact text once it is built and assembles the response from those texts — byte for byte the text of
    the whole tree, which is freed seed by seed), E the text the seed being walked is estimated to write (its
    modelled account, published at every level's end, / the account per text byte of the requested detail; from
    feature level 6 the account's record-coordinates share C apart: E = ⌈(A − C) / ratio⌉ +
    ⌈C / coordinate_account_per_text_byte⌉ with `coordinate_account_per_text_byte` 12, a bound — every byte of
    coordinate text, whatever its digits, costs at least 12 account bytes (each part's ratio at its widest text is
    14.9 to 21.9: an occurrence's 44 bytes for 872, a run's entry 221 for 3,680, a seed label's 79 for 1,728, the
    block's skeleton with a cut list's limitation, K record and `drop_coordinates` 960 for 15,151 in a graphlet,
    the null form 106 for 1,584) —, where an occurrence's
    account is 20–73 times its text (872 bytes for 12–44) against 115–129 for a tree's other output, which the
    measured ratio would have understated 2.6–4.1 times on coordinate-heavy seeds (`DESIGN-traverse-graphlet.md`
    §18.8, M2); a ratio is
    measured as (A − C) / (the text − the text the coordinates wrote, counted exactly from digits and
    punctuation), on a rest of at least 1 MiB, so that a request with coordinates measures what the same walk
    without them would, and the server's measured ratios do not depend on whether one came first),
    with the configured values `account_per_text_byte` (30 for the JSON details, 50 for a graphlet — *calibrated
    in the efficiency pass, feature level 4*: just below the smallest account/text ratios measured on real
    responses, 33.5 (UHGG `full`) to 1,344 (SRA `summary`; 37.2 the smallest SRA JSON ratio) and 58.4 and more for
    a graphlet; they stay conservative: a server's first large attempt of a detail can still be cut early until it
    has measured that detail itself — on a warm SRA server the first `tree` and `full` attempts of a 16S beam with a
    40 s bound were cut at 3-5 s with 20/40 and with 30/50 alike, the `tree` text measuring 115-129 account units a
    byte against the 30 assumed; the earlier report of 18.2 s instead of 11.0 s on a fresh server came from a cold
    walk against an intermediate build and did not reproduce; per-detail starting ratios were considered and not
    adopted, one locus being too narrow a measurement to start a reserve from),
    `compress_mbps` (`--traverse-delivery-compress-mbps`, default 50) and `build_mbps`
    (`--traverse-delivery-build-mbps`, default 10) until the server has measured its own: then, for each, the
    slowest rate — or for the detail the smallest ratio — of its last 16 measurements on `/traverse` responses
    and seeds' results of at least 1 MiB of text, and the attempt's own seeds' where more conservative, so
    that a machine slower than assumed is seen at once and a fast one is not held to the guess (both
    capabilities routes state the measurements). On the M5 Max of the measurements the configured guesses
    alone stopped the walk of a 403 MB SRA response at 2.3 s of its 35 s walk-until (account 35.6 GB, 88 per
    text byte, estimated at 20); measured, they hold it to its own needs. The factor 1.25 is a margin for rates
    that vary between responses, and `stop_ms` the time from the walk-until to the walk's end — the walk stops
    at its first poll that reads the clock after the walk-until (after a chunk of an annotation read, a
    lookahead's poll, or the heads between two readings of the clock), then finalises the stopped seed: `delivery_reserve.stop_ms`
    (`--traverse-chunk-target-ms` + 950: SRA attempts stopped 352 ms after their walk-until on a quiet fresh
    server and 1,001 ms on a loaded one, once 1,699 ms under load on a 400 MB result — a finalisation that grows
    with the result is covered by the 1.25 margin as long as it runs faster than 4 × `build_mbps`, about 235 MB/s
    for that result against 40 MB/s needed), or the longest such time the server
    measured over its last 16 attempts that walked past their walk-until when longer (`measured_stop_ms`). Both
    capabilities routes refer to the measurements these starting estimates come from
    (`delivery_reserve.calibration`, feature level 4: *the delivery reserve: calibration*, below); only the
    walk-until of attempts (`usage.bound.walk_until_ms`, and where a reserve stops a walk) can
    differ from the previous build's. Without the margin and the stop
    time a response the model fitted exactly was a 503 about half the time (review of pass 5, F3: a walk that
    stopped 124 ms after its walk-until left the delivery 1,277 ms of the 1,401 ms reserved, 98 ms short).
    `usage.bound.walk_until_ms` states, when the walk-until stopped the walk, the walk-until in force then;
    otherwise the lowest walk-until seen: the lowest that a poll reading the clock compared with — every 8th
    poll of a walk, every poll before a paced read's chunk or in the lookahead, and the poll before each seed —
    (the floor, `bound_ms − allowance_ms / 2`, when no such poll saw the reserve above it; the one computed once
    the last seed's text is written bounds no walk); a lower walk-until in force only between two such polls (the
    reserve moves at every level's end) is not stated *(review of 2026-10-06, C20: "the lowest in force" claimed
    every value)*;
    `usage.stopped_at` is when the walk saw the stop. The traversal routes compress at zlib level
    `--traverse-compression-level` (default 1; the other routes keep 9): on 10 real responses (541 MB of JSON)
    level 1 wrote 630 MB/s against level 9's 153 MB/s (graphlet bodies 51–64 MB/s at level 9) for 1.8 times
    the bytes, and the decompressed bytes are the same at every level. **When a 503 becomes likely**: before
    pass 5, as soon as the response's text exceeded allowance/2 × the machine's building-and-compression rate
    (about 87 MB on the staging server at level 9: 5 s × 17.5 MB/s measured; a few hundred MB on the M5 Max of
    the measurements); now only when the machine builds or compresses more slowly than the rates the reserve
    uses allow for with the 1.25 margin (a server's first large responses use the configured rates; later ones
    its measured rates, the slowest of its recent responses), when the seed being walked when its walk-until
    passes writes more text than its account / `account_per_text_byte`, or when the walk ends later after its
    walk-until than `stop_ms` (a read decoded whole because it was predicted to end well before the deadline,
    a head that runs long; a longer time is measured and kept for the next attempts). The configured rates are conservative guesses not calibrated
    on the staging server (both capabilities routes state them and the measured ones; they are
    configuration). So a large response is no longer
    lost whole for want of time: the reserve stops the walk early enough to deliver what was walked, with the
    rest stated (`attempt_deadline`, `not_started`). The reserve can stop a large attempt earlier than before
    pass 5; what it did not walk is stated, as any stop.
    **Stated limits** — the uninterruptible pieces, each of which a stop waits out (a cancel is seen after
    the piece it falls into, a walk-until after the pieces up to the next poll that reads the clock, the bound
    after those up to the next delivery check, above): one chunk of an
    annotation read (at least one row; no bound on one row's decode exists before stage 3c, so
    `deadline_check.max_uninterruptible_ms` is null — and it stays null after stage 3c-ii, whose checkpoints bound
    index operations, not time: page faults and scheduling leave the time open (decision 3c-N5) —, and
    `usage.observed_max_uninterruptible_ms` states the
    longest read, chunk or head piece of the attempt's walks, below), or a read decoded whole because it was
    predicted to end well before the walk-until, whose rows were more than 64 times slower than any read before
    (*chunked deadlines*, below); the mapping of the seed's k-mers to nodes before its validation,
    which is not polled (up to `--traverse-max-seed-bp` = 100 000 k-mers on the server, unbounded in the CLI;
    review of the stage-4 backend, F7); the seed phase's own processing — resolving the named labels, their
    duplicate check, the extra labels, the seed id, the depth-0 state — linear in the bytes of the labels named,
    which no server limit caps; a head's processing between two work checks; the lookahead's chains (§8.3)
    between their polls — every 16 graph steps and before each chain's key mapping, which is one call over at
    most `annotation.batch_kmers` nodes — and its clearing once it outgrows its bound (1,000,000 entries without
    a budget; linear in its entries, at most a level's heads × `min(batch_kmers, the radius left)`); a read's
    preparation before its first chunk — its cache lookups and the ordering of its keys in the walk's order, n log
    n in its keys (a lookahead read's are up to a level's heads × `batch_kmers`: 90–260 ms for the 1 to 3.8 million
    keys of the reviewer's W3 fixture, the longest head piece there once its chains were polled); a seed's
    finalisation and summary (linear in its admitted result); the building of the response between two delivery
    checks; and the transport, which sends the written bytes after the handler returned (the content timeout
    bounds that). *(Review of 2026-10-06, W3: the chains read only the seed's own deadline; at `batch_kmers`
    60,000 a cancel was seen 1.5–11.8 s late and a walk-until passed by 0.6–6 s, an HTTP 503 at the bound, while
    `observed_max_uninterruptible_ms` said 1 ms. D12: the seed phase's processing was missing here.)*
- **The delivery reserve: calibration** (the capabilities' `attempts.delivery_reserve.calibration` refers here).
  The reserve's configured values are starting estimates, each replaced by this server's own measurements once it
  has them (the `measured_*` fields of `delivery_reserve`, above):

  | starting estimate | value | measured | why this value |
  |---|---|---|---|
  | `account_per_text_byte.json` | 30 | account bytes per text byte on real responses: 33.5 (UHGG `full`) to 1,344 (SRA `summary`), 37.2 the smallest SRA JSON ratio | just below the smallest ratio measured; one starting ratio for every JSON detail, since one locus is too narrow a measurement to start a detail's reserve from |
  | `account_per_text_byte.graphlet` | 50 | 58.4 and more | just below the smallest ratio measured |
  | `stop_ms` | `chunk_target_ms` + 950 | from the walk-until to the walk's end: 352 ms on a quiet SRA server, 1,001 ms on a loaded one, 1,699 ms once (under load, a 400 MB result) | 950 ms for the heads between two readings of the clock and the stopped seed's finalisation; the rest of a finalisation that grows with the result is covered by the margin while it runs faster than 4 × `build_mbps`; a longer stop measured replaces it (`measured_stop_ms`) |
  | `compress_mbps` | 50 | 460–670 MB/s on SRA responses at zlib level 1 | conservative |
  | `build_mbps` | 10 | 13.6–51 MB/s on SRA responses | conservative |

  They stay conservative, so a server's first large attempt of a detail can still be cut early until it has
  measured that detail on its own responses: on a warm SRA server the first `tree` and `full` attempts of a 16S
  beam (a 40 s bound) were cut at 3–5 s, the `tree` text measuring 115–129 account units a byte against the 30
  assumed.
- **Chunked deadlines** (pass 5; `LabelOracle::pacer`, `ReadPacing`). Under a deadline — the seed's
  `bounds.time_budget_ms` where the walk reads it (depth > 0, as a checkpoint does; in a derivation, a positive
  budget, as the derivation reads it after every k-mer) and, for a request with `attempt_id`, the attempt's
  walk-until (everywhere) — every annotation read of `/traverse` that the deadline may fall into is decoded in
  **chunks**, the deadline checked before each. A read is **one piece**, as before pass 5, when at the slowest
  per-row time the request has seen (any piece, any read: pessimistic, which costs at most a first chunk) it
  would take less than 1/64 of the time left — rows up to 64 times slower than any seen before still end before
  the deadline — and without a deadline. Otherwise its first chunk is at most 8 rows (it measures this read's
  own rows: another read's rows, another level's, can be two orders of magnitude cheaper), each next chunk at
  most 4 times the previous one and sized at the rate the previous chunk measured to take `chunk_target_ms`
  (`--traverse-chunk-target-ms`, default 50, on the server and in the CLI) or the time left, and the rest of
  the read is one piece once predicted at that rate to take less than a quarter of the time left. The chunks
  are taken in the walk's order (a lookahead's rows path by path), each sorted as a whole read is: on a
  row-diff annotation the rows of one call share the decoding of their row-diff paths, which a chunk of rows
  from as many paths decodes again (review of pass 5, F1: chunks of sorted rows took 0.75 ms a row against
  0.016 ms for the read whole, and a row_diff walk of 1.8 s took 30 s and was stopped by its budget). So
  splitting costs a few small chunks per read near its deadline and nothing far from it. Chunked are a level's fetch, the
  lookahead, a seed's validation — by the attempt only: a seed's own time budget never stops its validation
  (the deadline is never checked at depth 0, and a seed validated past its budget still delivers its result
  complete to 0 bp) — and a derivation's window. The lookahead is paced by the seed's time budget at every
  depth, depth 0 included, its chains' graph steps stop at it, and it is skipped once the budget has passed
  *(efficiency pass)*: run() reads that budget before depth 1 and a head's checkpoint before the other arm's next
  head, so a lookahead read past it is never consumed (before, a chain's steps and the roots' lookahead read ran
  in full: 2.2 s on a 250 ms budget on refseq33m, 1.3 s of it graph steps of 1,000-node chains). **The chains
  read the attempt's stop too** *(feature level 6, the review of 2026-10-06, W3)*: they run between two
  checkpoints and charge no units while they are enumerated, so every 16 graph steps of the level's lookahead
  (`kLookaheadPollSteps`) and before each chain's key mapping they read the seed's deadline and, for a request
  with `attempt_id`, a cancel, the walk-until and a gone client (the poll that reads the clock, as before a
  paced read's chunk); a stop ends the lookahead there — the chain being built is dropped, nothing is read —
  and the next head's checkpoint (at depth 0 run(), before depth 1) stops the walk on it: the lookahead reads
  the attempt's stops only through that clock-reading poll, which records a passed walk-until with the
  attempt, so the checkpoint's poll returns it whatever its stride *(the review of the P2 fixes: without
  pacing, `--traverse-chunk-target-ms` 0, the lookahead stopped on the walk-until without that poll, and the
  walk ran up to 7 more heads before a poll that read the clock stopped it)*. Nothing depends on
  the cache, so a walk no stop ends is unchanged; one a stop ends differs only in timing and the physical
  counters of reads it no longer makes. A chunked read **returns exactly what one read returns**: its
  counting (`annotation.rows_requested`, the cache hits), its cache's eviction and its result are decided once
  for the whole read and only the decoding is split (a budget-aware read admits every key against its
  standalone demand, so its runs may be cut anywhere), so a request that no deadline stops is byte-identical
  to the unchunked request, `annotation.direct_reads` and `keys_mapped` included, whatever the chunk sizes
  (checked against the previous build on 1,380 real requests — 588 mini_refseq, 792 UHGG — and on the 588 with
  chunks of 1 ms, and in one-row chunks, every read split, by the unit tests). A deadline that stops a read stops the walk at the read, as the check after the whole read would
  have: a level's fetch at the level's first head (`time_budget`, `cancelled` or `attempt_deadline`), the
  lookahead silently (the next head's checkpoint, which reads the same deadline, then stops the walk), a
  validation by failing the seed (the attempt's stop, §7.0), a derivation with `time_budget` after the k-mers
  consumed before its window (or the attempt's stop). The unchunked walk ran the read to its end, past the
  deadline, and could walk on (cached levels) before its next check, so a walk stopped in a read can end at a
  smaller `complete_to_bp` than the unchunked walk that overran (on the `row_diff_brwt` fan-out index at 100
  ms, 5 where the overrunning walk reached 69). Only such a walk's counters (`rows_requested`, which a stopped read does not
  add) and the time values differ from the unchunked request's — a stop by time is not deterministic anyway
  (§6.8, determinism). The rows a stopped read's chunks decoded are **charged** as
  work (decoded, though no row was returned; review of the stage-3 fixes, P3: what was decoded is in
  `usage.work_units`): to the level's arm, or to the seed phase (`work_seed`) without a comparison, as a
  window found `too_wide`. One chunk stays uninterruptible: at least one row, and no bound on one row's decode
  exists before stage 3c (selected-label decoding) — `deadline_check.max_uninterruptible_ms: null`, which stays
  null once stage 3c-ii's checkpoints bound the operations of a read (decision 3c-N5);
  `deadline_check.observed_max_uninterruptible_ms` is the longest single piece the server process saw — a read
  or chunk (a whole read far from its deadline included), a **head piece** (the walk between two readings of
  the clock for a stop, its reads excluded, the lookahead's chains between their polls included; feature level
  6, the review of 2026-10-06, W3: only reads were counted before) or a gap between two delivery checks: how
  late a cancel can be seen —, and `usage.observed_max_uninterruptible_ms` the attempt's reads and head pieces:
  observations, not bounds. Neither counts the mapping of a seed's k-mers, the seed phase's own processing, a
  derivation's steps or a seed's finalisation (each seed's `timing.deadline.longest_piece` names its longest
  piece of any kind), so neither is a margin a ledger can add (§10.3, the release rule). **Stated limits**: a read
  decoded whole because it was predicted to end before 1/64 of the time left overruns the deadline when its rows
  are more than 64 times slower than any the request read before; the rest of a split read, when its rows are
  more than 4 times slower than its chunk's. **Not chunked**: `/resolve`
  (stage 3b; its own deadline, `bounds.time_budget_ms`, is read between its row batches, §4.5), the mapping of a seed's k-mers, the seed phase's own processing (linear in the bytes of the
  labels named), a head's processing (checked every `work_check_interval` units), a lookahead chain's key
  mapping (one call, at most `batch_kmers` nodes; the chain's graph steps are polled, above), a read's
  preparation before its first chunk (its cache lookups and the ordering of its keys, n log n in them) and the
  lookahead's clearing once it outgrows its bound, a seed's finalisation and summary, the response's building
  between two delivery checks (every 64 KiB of text and every 4,096 objects of a seed's tree; one token's
  preparation in between), and the transport.
  `--traverse-chunk-target-ms 0` decodes every read in one piece, as before pass 5. Measured on the fan-out
  index of `TestTraverseWideIndex` (1,024 header labels, one lookahead read of about 66,000 rows), on a loaded
  M5 Max, budgets of 50, 100, 200 and 500 ms: with a column annotation, which the previous build overran by up to
  0.9 s, kept to within 15 ms; transformed to `row_diff` (whose depth-6 fetch of 1,024 unshared rows takes 1–2.6
  s) and `row_diff_brwt`, with and without a memory budget, which the previous build overran by up to 1.0 and
  1.9 s, kept to within 60 ms (the largest piece 106 ms); and a walk its 30 s budget does not stop returns the
  same bytes in the time of the unchunked walk (one or two chunks of 8 and 32 rows per large read). On 60 real
  UHGG requests the walks took 15,973 ms against 15,932 ms for the previous build and 16,021 ms unchunked.
  *Review of pass 5, finding 5:* the budget-aware lookahead (`LabelQuery::warm`, `LabelRecorder::warm` with a
  budget) cut its runs from the globally sorted missing keys; its runs and their pieces are taken in the walk's
  order now, each piece sorted, as the default reads take their chunks (32 independent paths, 2,048 walked rows,
  runs of 8: 31,297 decoder charges against 161,679, the walk-ordered runs alone 28,736). Only what is cached
  depends on it: a key is admitted and charged by its own costs when a level's fetch reads it. *Review of
  2026-10-06, U05-01:* such a lookahead ends at the first run the cache cannot keep beside its own runs (a cache
  of 8 keys: the walk's first run of 8 decoded and kept, 130 charges, where every run was decoded and the last
  kept; the decoding a memory budget still adds is stated under the refusals below). *R9:* the tests of
  this paragraph (`WalkerDeadlineChunks`) run on a virtual clock (`DecodePacer::test_clock_ms`) that only their slow
  reads advance, so their bounds are exact and they no longer fail under load (0 of 600 executions with 24 busy
  processes; 5 of ~400 before).
- **Delivery checks** *(review of pass 5, finding 6; feature level 5)*. The text of each seed's result and of the
  response is written under the attempt's delivery check (its client and its bound) every 64 KiB: a piece the JSON
  writer hands over, and a seed's text appended to the response, are copied in pieces up to the next check — before,
  one 16 MiB string value was appended whole, one check after it, in the writing and in the assembly alike. What
  stays uninterruptible is what the writer does between two pieces: preparing one token (escaping one string
  value, e.g. a graphlet of many MB, before it is copied). The longest time between two checks of a written text is
  measured and added to `deadline_check.observed_max_uninterruptible_ms` (not to `usage`, which states the reads
  and, from feature level 6, the head pieces). *Delivery memory (review of 2026-10-06, C9):* the response's text
  is assembled into a string of its exact size, each seed's text freed once copied into it, and the transport's
  copy is made into a buffer sized once, after the last check. The seeds' texts lived until the handler
  returned, beside the assembled text and the transport's buffer, both grown by doubling: measured on
  mini_refseq (n identical seeds, a fresh server each, jemalloc), the server's peak grew by 2.60 bytes per byte
  of text with gzip and 3.51 without (185 MB of text: 590 and 764 MB), and now by 1.52 and 2.57 (394 and 584
  MB) — not by one and two: jemalloc keeps the pages of the texts freed during the assembly for a while. No byte
  changes, and nothing of it is promised: no budget bounds a whole response (§6.8, memory_bound_soft). The
  compressor takes the text in pieces of at most 2^32 − 1 bytes, zlib's input counter, so a text of 4 GiB or
  more is compressed whole; it was cut to its size modulo 2^32 and answered 200 *(U13-04; an output change on
  every route that compresses, a level-6 correction, §10.3)*; a text below 4 GiB is compressed by the same
  calls, its bytes unchanged.
- **The deadline record** *(R8; feature level 5; in `timing` only, never in an untimed body)*. Per seed,
  `timing.deadline`: `longest_piece {ms, kind, rows, coordinates}` — the seed's longest uninterruptible piece,
  `kind` one of `read` (an annotation read in one piece), `chunk`, `rest` (the piece that ended a split read),
  `kmer_mapping` (the seed's k-mers mapped to nodes, then to keys), `coord_mapping` (a derived seed's k-mer whose
  coordinates were mapped to headers), `derivation_step` (a derived seed's k-mer otherwise), `head` (the walk
  between two readings of the clock, its reads excluded; the last one ends at the walk's stop or end), `setup`
  (the seed phase's own processing between its other pieces: validation, resolving the label names and their
  duplicate check, the extra labels, the depth-0 state — everything before the walk's first checkpoint that no
  other piece covers) or `finalisation` (from the walk's stop or end to its result), with its rows and the
  coordinates mapped to headers in it — the pieces cover the seed's time from its start to its result —, and,
  when a stop ended the walk, `stopped_by` (`time_budget`, `attempt`, `cancelled`, `memory`, `work`) and `stop_after_deadline_ms` (for
  `time_budget` and `attempt`: how long after the seed's budget, or the attempt's walk-until, the walk saw it). A
  read maps its rows' coordinates to headers inside its piece, so the rate a chunk is sized by includes the
  mapping (checked: the staging overrun of a warm tuple walk, 6,136 ms on 5,000, could not say which piece was
  late; it now can). Beside it, `timing.seed_phase_ms` (the seed phase: from the seed's start to the walk's first
  checkpoint, or to its stop or end when it reaches none — the validation or derivation, the extra labels and the
  depth-0 state, a failed seed's included), `seed_fetch_ms` (its reads, also in `annotation_fetch_ms`) and
  `label_resolve_ms` (resolving its label names), and, on a row-diff
  annotation with the path cache, `timing.path_cache {hits, stored_rows_read, rows_kept, bytes_kept, peak_bytes}`
  (§8.4): physical, timing only. *Review of levels 4–5, finding 3:* a valid request naming 2,500 headers with a
  shared 1,024-character prefix stated, on a warm server, a seed phase of 126 ms and a longest piece of 0.297 ms:
  resolving the names (3 ms) and their pairwise duplicate check (quadratic; most of the rest) ran between the
  validation's k-mer mapping and its read, where no piece was measured, and the first head's checkpoint did not
  count what came before it. The seed phase is now cut into pieces at its reads, k-mer mappings and derivation
  steps, the spans between them `setup` pieces, flushed on every way out of the seed (a failed one too);
  `seed_phase_ms` runs to the first checkpoint (before: the validation or derivation alone; it now also counts
  the depth-0 state, and `seed_fetch_ms` an annotate root's read); the last `head` piece runs to the walk's stop
  or end (before: no piece); the duplicate checks of the seed labels and of the extra labels are hash sets (the
  same refusals, in the request's order). The same request: seed phase 13.9 ms, longest piece 9.8 ms (`setup`:
  mostly the seed id over the 2,500 names, linear) against 139.1 and 1.28 (`read`) before, from the CLI. Timing
  only: no untimed byte changes.
- **Annotation reads under a budget** (stage 3 of `DESIGN-traverse-graphlet.md` §14.1). On an annotation whose
  reads are budget-aware — `RowDiff` over BRWT or ColumnMajor, with or without coordinates
  (`row_diff_brwt`, `row_diff_brwt_coord`, `row_diff`, `row_diff_coord`) — a request with a budget reads its
  annotation through the decoder's budget-aware path (`IRowDiff::decode_rows` / `decode_row_tuples`); every
  other read, and every read without a budget, is the default decode, unchanged.
  - **Work.** A returned row's work also counts its **dependency rows** — every other row on its row-diff path to
    the anchor: 8 units each and 1 per entry and coordinate they store — whichever read decoded them (a fetch,
    the lookahead, an earlier level). The weight is a property of the row, kept with it in the caches, so work
    stays independent of `annotation.batch_kmers`; it is conservative (a batched read shares dependency rows
    between keys, and the charge does not).
  - **Memory.** Under a memory budget the decoder charges every buffer before allocating it: the row-diff trace,
    each dependency row (read alone by the index matrix's single-row descent), each coordinate tuple (its size
    read from the delimiters before it is copied), the reconstruction's buffers, and what building a row's hits
    or label list holds beside it. The bytes are a deterministic model (jemalloc's size classes and the
    containers' growth), never a measurement; the tests check it against jemalloc's own peak. A read is charged
    against what the request has left — the budget minus the admitted account and what the level's fetch holds
    beyond it — and one that does not fit is **refused whole**: nothing is returned, and no cache, dictionary or
    result counter (`annotation.rows_requested`) changes; the physical counters in `timing` (`rows_fetched`,
    `tuple_rows_fetched`, `coords_mapped`, `annotation_fetch_ms`) count every decode made, a refused read's
    included. In annotate mode a key's admission also counts the dictionary labels it names first, each priced as
    the account charges a dictionary label (its entry and its delivery in the requested detail) plus a fixed 640
    bytes that bound what naming it provisionally holds (one flat table and list per call, their growth
    included); the names are charged with the dictionary on success, the naming charge stays with the level
    until its heads are processed, and a refused read names nothing.
  - **Admission of every returned row.** Every key a fetch returns is admitted against its **standalone
    demand**, a deterministic upper bound of what reading it alone holds (the decode, building its hits, the
    result), computed with the row and kept with it in the caches. A key is admitted, cached or not, exactly when
    its demand fits what is left at its position in the level (the budget minus the account and the rows the
    level's earlier keys returned), so that where a level stops does not depend on how it is cut into calls nor
    on what the lookahead cached. Keys are decoded in runs; a run that does not fit is retried in halves, and a
    key that does not fit alone is the stop.
  - **What a refusal does.** In a level's fetch it stops the seed like a refused head at the level's first head:
    every head of the level ends `resource_limit` and the arm is complete to the level's depth. The stop names
    its cause (§7.0): a row that does not fit what was left (its standalone demand, or its read alone with its
    dependency rows) stops in phase `annotation_decode`, stated as needing more than the bytes left — whether the
    exact demand was known depends on what the lookahead had read, so the statement does not give it; in annotate mode a row that fits but whose new dictionary labels do not
    stops in phase `traversal` (the dictionary is walker state), with the labels' number and bytes; a level whose
    own key and successor lists, held beside the account, leave none of the budget for its read stops in phase
    `traversal` too, before reading anything. In the seed phase — the derivation's rows, the validation's — or at
    an annotate root a row that does not fit fails the seed in phase `annotation_decode` (no result complete to 0
    bp exists yet). At an annotate root, labels that fit beside the row but not beside the depth-0 state, and a
    depth-0 state that already reaches the budget before a root's read (the other arm's root, its labels and
    reservation), fail the seed as a depth-0 state (phase `traversal`), without naming the labels or reading the
    row, its need stated as a lower bound ("at least"). The lookahead decodes within what is left and gives up
    silently; it also ends at the first run the label cache cannot keep beside the runs it cached, or cannot
    keep at all (only entries cached before it may be evicted, for its first run) *(review of 2026-10-06,
    U05-01: a lookahead larger than the cache decoded every run and evicted its own earlier runs, keeping the
    last one or two — 31-39% more rows decoded at `batch_kmers` 2,048-8,192 under 16 MiB, every result the
    same; only physical work, time and the `timing` counters change)*. **What a memory budget still costs in
    decoding (physical work, stated, not fixed):** a lookahead's first run still evicts the label cache
    wholesale when it does not fit beside what is cached, and with it the rows earlier lookaheads read ahead
    that the walk had not reached yet; the walk decodes those again. On mini_refseq, in walks the budget does
    not cut, the rows decoded exceed those of the same walk under a work budget alone (whose caches have no
    byte bound) by up to 20% at 16 MiB (39% before the fix above), 52% at 32 MiB and 32% at 64 MiB (50% and
    31% before) at `batch_kmers` 2,048-8,192, and by up to 0.9%, 4.3% and 8.0% at the default 64 (unchanged).
    The fix is not uniformly fewer decodes: of 180 walks under 16-64 MiB (three seeds, both modes, one and
    two arms, `batch_kmers` 64-8,192) it decoded fewer tuple rows in 24 (up to 47% fewer), more in 10 (up to
    1.4% more; up to 4.4% more with the stored rows the path cache read), the same in 146, every result the
    same. Evicting the oldest rows first instead (a lookahead keeping a quarter of the cache) was measured and
    not taken: it cut the measured walks' tuple rows by 12% but kept the label cache full, and the row-diff
    path cache, which holds only what the label cache leaves of their allotment (§8.4), read 2.5 times the
    stored rows. Not changed either: the structural lookahead clears itself when a level's chains (up to
    `batch_kmers` per head) would pass its allotment (`max_memory_mb` × 512 entries) and builds and reads them
    again at later levels — on mini_refseq up to 8.4 million tuple rows for walks that decode 10,000-50,000
    under a work budget alone (4-16 MiB, `batch_kmers` 2,048-8,192), every result the same. On a format whose
    reads are not budget-aware nothing is refused inside a read, and an annotate level's new dictionary labels
    are charged after its read: when they put the account over the budget the level is censored at its first
    head as at a refused read, and the stop names them — phase `traversal`, their number and bytes, the need a
    lower bound ("at least": the account with them; its heads need more), the levers of a stop by labels *(review of
    2026-10-06, U03-03: it was stated as a refused head, the dictionary's need as exact, so that raised to it the
    walk stopped at the same depth again, and without `lower_max_labels_per_node`)*. Under a request budget the
    label caches never exceed their allotments, and the
    structural lookahead is cleared before a chain would push it past its allotment (counted keys: each key the
    walk consumes once, `annotation.keys_mapped`).
  - The derivation reads 64 k-mers per window whatever `annotation.batch_kmers` is, under a budget or not (§6.1;
    feature level 6), and holds the window's rows at once (in the first it chooses the cheapest row). The window's
    distinct rows are read in runs, a run that does not fit retried in halves, every row admitted against its
    standalone demand beside the rows read before it. A window whose rows do not fit together fails the seed,
    and the statement names the refused row's seed k-mer and demand (or that its read alone did not fit), what
    the seed phase had left, and what the window's earlier rows held — conservative (each row alone may fit),
    and stated.
  Formats without the budget-aware path are read by the default decode under a budget too; what their reads hold
  is observed, not charged (§7.0, `memory_bound_soft`). **The formats with budget-aware reads are exactly the four
  above**; every other format (`column`, `column_coord`, `brwt`, `row_diff_disk`, `row_flat`, …) keeps stage 2's
  observed-only reads. **Stage 3 does not establish a hard request-wide memory bound** (review of stage 3, answer
  5): on every format `memory_bound_soft` names what is still held uncharged beside the admitted account — a
  level's keys and fetched rows until its heads are processed, the seed phase's intersection, hits and
  coordinate copies, a failed seed's depth-0 dictionary, a failed result's echo of `seed_id`, and, not observed,
  an index-wide header lookup and the label dictionary's first table — and each seed of a request has budgets
  of its own (§5), so a request of n seeds holds up to n budgets.
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
  | `walks` | `complete` \| `partial` \| `failed` | `complete` only when no walk-class limitation applies: every requested arm has `status: complete` (`complete_to_bp == bounds.max_extension_bp`) **per path** (keep, or merge where no merge united a history) and no carrier of the seed was left out. `partial`: a valid certified prefix whose limits are stated — a `walk_domain` (an arm was truncated or pruned), a `scope` with `observed > 0` (a merge united histories: the walks are complete for the united-history rule, not per path) or a `seed_labels` limitation (the walks only a cut carrier carries are missing). `failed`: no traversal — the permitted set could not be derived (§6.1 step 4): the result has no `arms`, an `error`, and a `derivation` limitation; or a request budget does not hold the seed itself (§5, request budgets): no `arms`, an `error`, a seed-level `walk_domain` naming the budget's knob and a `resource_stop`; or a stop from outside the walk reached the seed in its seed phase or before it started (§6.8: a cancel, the attempt's bound): the same shape, the `walk_domain` naming `attempt_id`, the `resource_stop` with `scope: "attempt"` (`phase: "not_started"` for a seed never started) |
  | `branch_diagnostics` | `complete` \| `cut` | `cut`: a `branch_events` limitation — some arm's `evidence.complete` is `false`, and branch decisions and refusals at or beyond `evidence.complete_to_bp` are not reported |
  | `label_evidence` | `complete` \| `lower_bound` \| `qualified` | `lower_bound`: evidence may be *missing or understated* — recorded lists were cut (`label_lists`, `inexact_counts`: annotate mode's `label_summary` and the counts flagged `exact: false` are lower bounds), carriers of the seed were cut (`seed_labels`), a `max_switch_sources` cut may have raised a loss or missed a switch entry (`switch_sources`), or losses were re-minimised greedily (`greedy_losses`). `qualified`: something reported may be *overstated* — a column-label trace cannot see a record boundary (`trace_record_boundaries`), or the permitted set was derived from part of the seed (a `derivation` on a walked result, D3, feature level 6: it may hold labels the whole seed would exclude); `qualified` wins when both apply |
  | `delivery` | `inline` | the whole result is in this response; `spooled` / `paged` are reserved for the graphlet delivery path (`DESIGN-traverse-graphlet.md` §14) |

  **Conservative by construction** (`DESIGN-traverse-graphlet.md` §14): each axis is computed from the
  result's own `limitations` (seed level and every arm), so an axis reads `complete` exactly when no limitation
  of its class is stated — a reader never finds a limitation whose axis still says `complete`, and never a
  non-complete axis without the limitation that explains it. `server_clamp`, `coordinates` (decision C2:
  record coordinates are opt-in output, whose own completeness their block states) and a failed seed's
  `derivation` belong to no class: what a clamp caused is stated by the `walk_domain` or `seed_labels` entry it
  led to. A `derivation` on a **walked** result — a permitted set derived from part of the seed (D3, below) — is
  in the `qualified` class: such a set may hold labels the whole seed would exclude. A failed derivation
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
| `walk_domain` | arm; seed, on a seed a request budget failed, or the attempt stopped in its seed phase or before it started (no `complete_to_bp`: nothing was walked) | `status` is not `complete`: the cap in `cap_trigger`, then any other cap that ended walks after it (a beam, then a step cap) | the cap's field: `bounds.max_steps`, `bounds.max_live_paths` (also a beam's width), `bounds.max_paths`, `bounds.max_output_bp`, `bounds.time_budget_ms`, `bounds.max_memory_mb`, `bounds.max_work_units` (`resource_limit`); `attempt_id` for a stop from outside the walk (§6.8: a cancel or the attempt's bound, `resource_limit`; `limit` the attempt's bound in ms, and the effect says no budget ran out) | what the cap compared at the trip: steps of the seed, bases of the arm, leaves + live heads, live heads (beam: heads of the level), elapsed ms, the smallest budget in MiB that admits the head — what admitting it needed with that budget's own caches' allotments, which grow with the knob (§6.8) —, the work units used, the attempt's elapsed ms; for a later cap the walks it ended. `complete_to_bp` = the arm's |
| `branch_events` | arm | events were dropped by the cap | `output.max_branch_events` | `branch_events_total`; `complete_to_bp` = `evidence.complete_to_bp` |
| `label_lists` | arm | `labels_per_node.nodes_truncated > 0` | `labels.max_labels_per_node` | `labels_per_node.max_seen` |
| `inexact_counts` | arm | a live-label count in `frontier_remaining`, `cap_trigger` or a `growth` bin is flagged `exact: false` | `labels.max_labels_per_node` | how many counts are flagged |
| `switch_sources` | arm | a `table` cost's source list was cut while a cut source had a finite switch into a target of that successor (`counters.switch_sources_cut > 0`), or a label ended `label_lost` with the qualifier `switch_sources` — whether or not the cut source ended: one that goes on along another successor leaves no end, yet the target it was the cheapest way into was entered at a higher loss or not at all | `labels.max_switch_sources` (accepts `"unlimited"`) | `counters.switch_sources_cut`, the successor derivations the cut may have changed (an over-approximation: the kept sources may still have been the cheapest); `label_ends` = the `switch_sources` label ends. Effect: losses may be overestimated and switch entries missed |
| `scope` | arm | `completeness_scope: united_history` | `branching.on_reconverge` (`limit: "merge"`) | merges done (0: no history was united) |
| `greedy_losses` | arm | an exclusion re-ran the derivation with priced switches under a finite branch limit | `branching.max_label_branches` | re-minimisation rounds after the first |
| `trace_record_boundaries` | seed | `support: trace` with column labels | `labels.seed_label_kind` (`limit: "column"`) | column labels in the dictionary |
| `seed_labels` | seed | the derived permitted set was cut (`labels_dropped > 0`); on a failed result, carriers were cut before the trace check (`no_trace_carrier`) | `labels.max_seed_labels`, with `server_limit` when the server clamped it | `labels_supporting_total` |
| `server_clamp` | seed | an entry of `strategy.clamped` bound this seed: a lowered derived-set cap that cut its set, a lowered time budget that tripped, a budget raised from zero (the walk ran under it), or a memory or work budget at the server's maximum (§10, R16) that stopped its walk or failed the seed | the clamped field | the requested value (`limit` is the effective one) |
| `memory_bound_soft` | seed | `bounds.max_memory_mb` is set: something is always held beyond the admitted account. **(i) On an annotation whose reads are budget-aware** (§6.8; `ResourceAccount::decode_charged`) the reads are charged inside the decoder, and the statement names only what is still uncharged: a level's key and successor lists and its fetched rows until its heads are processed, the seed phase's intersection and hits, the dictionary a failed seed's depth-0 state built, a failed result's echo of `seed_id` (named in the effect, review of the stage-2 recheck, design answer 5), an index-wide header lookup that the first request needing it builds (not observed: charging it would make a result depend on the server's history), and the label dictionary's first table (about 1.5 KB per request) with the transient copy while its list grows (not observed). **(ii) Other formats**: the budget is enforced on the modelled state and output, but the annotation rows the seed phase and each level decode (with an annotate dictionary's growth and a cache beyond its allotment), and a failed result's echo of `seed_id`, are held before they can be charged. Neither is a hard request-wide memory bound (§6.8). In no outcome class | `bounds.max_memory_mb` | the largest excess over the budget seen of what was held beyond the admitted account, MiB rounded up (0: none), observed wherever it is held and before any check after it can stop the walk, so a stop inside a fetch still states what the fetch held; the admitted account itself never exceeds the budget. Observed once a level's key and successor lists are built, after every call of a level's fetch (with the call's raw rows in (ii)) and at a refused read (with the rows read before it), at **every head admission** (the account the head needs with the level's lists and fetched rows beside it), after the roots' rows and the lookahead's warming, and over the seed phase's rows and intersection (after each sub-batch is read and again once it is consumed; each fetch call of the label validation before its charge, which can fail the seed; under `trace` the validation's peak — the hits, every label's live coordinate set and both arms' boundary coordinates — before anything can fail the seed, review of the stage-2 recheck, P2). In (ii) the decoder's own transient buffers are not counted (the materialised rows are). Also on a seed a budget failed — what its depth-0 state built before the admission refused it (the dictionary's labels with their names in every copy, the caches) is observed first — on a seed whose derivation failed (what the derivation's rows held), and on every failed or refused seed what its result's echo of the request's `seed_id` holds beyond the budget (§5) |
| `coordinates` | seed | `output.coordinates` is true and an occurrence list — a seed label's or a run's — was cut at the cap to its first occurrences by start (each cut list states `occurrences_total`; the block `complete: false`). In no outcome class (decision C2) | `output.max_coordinate_occurrences` | the largest true count of a cut list; the extra field `lists_cut` counts the lists cut |
| `derivation` | failed seed; a walked seed whose set was derived from part of it (D3) | the permitted set could not be derived (`outcome.walks: failed`); `cause` names why, `server_limit` is added when the server clamped the knob; also a seed refused for a label name that is not UTF-8 (`unrepresentable_label_name`, §6.1 step 4, either mode). **D3** *(the owner's decision, feature level 6)*: a derived seed whose time budget runs out after 1 ≤ j < n of its k-mers is delivered as a walked result, without `error`, when nothing else would fail it — the set of the labels carrying the j k-mers read (a superset of the whole seed's carriers; `labels_from_seed`, `labels_supporting_total` its size), every requested arm truncated at `complete_to_bp` 0 (`cap_trigger` `time_budget`), this limitation first with `cause: time_budget` and `observed` = j, no `num_kmers` field, `outcome.walks: partial` and `label_evidence: qualified`; j = 0, or a superset that would fail the seed in any other way, fails as before | per `cause`: `no_carrier` → `seeds[].sequence`; `no_trace_carrier` → `support` (`limit: "trace"`); `too_wide` → `seeds[].sequence`; `time_budget` → `bounds.time_budget_ms`; `ambiguous_header` → `labels.seed_label_kind` (`limit: "header"`); `over_seed_label_cap` (`exhaustive`) → `labels.max_seed_labels`; `unrepresentable_label_name` → the knob that avoids recording the label: `labels.mode` (`limit: "annotate"`) in annotate mode, `labels.seed_label_kind` (`limit: "header"`) for a derived header label, `labels.extra` (`limit`: its labels) for an extra one, else `seeds[].labels` (`limit`: the labels named or derived) | `no_carrier`: the seed k-mers read when no candidate was left (`limit`: the seed's k-mers); `no_trace_carrier`: the labels carrying every k-mer by presence; `too_wide`: annotation entries of the narrowest of the first 64 k-mers (`limit`: 64 · `max_seed_labels`, at least 65 536); `time_budget`: **two units** — on a failed seed the elapsed ms (a number), on a walked seed (D3) j, the k-mers read (an integer; n is the seed's `num_kmers`, and the effect says "after j of n k-mers"): one `kind`, `cause` and `knob`, told apart by the result's shape (walked: `arms`, no `error`), which the owner kept rather than change the wire within level 6 *(review of 2026-10-06, X1; in an MGT `K` record j is written `i:j`, an elapsed time `f:`)*; `ambiguous_header`: the header (under `bounds.max_memory_mb` one longer than 256 bytes as a bounded prefix with its length, column and sequence id, §5); `over_seed_label_cap`: the carriers; `unrepresentable_label_name`: the names that are not UTF-8 (the effect names the first by column and sequence id, never its bytes) |

```json
"outcome": {"walks": "complete", "branch_diagnostics": "cut", "label_evidence": "complete", "delivery": "inline"},
"evidence": {"complete": false, "complete_to_bp": 25},
"limitations": [{"kind": "branch_events", "knob": "output.max_branch_events", "limit": 1, "observed": 3,
                 "complete_to_bp": 25, "effect": "branch decisions and refusals at or beyond complete_to_bp are
                 not reported; more were produced: raise the knob or set it to \"unlimited\""}]
```

- **A resource stop**, per seed result, when a request budget stopped the walk (and for a time stop only when the
  request set a budget, so that a response without one keeps its form): `resource_stop: {scope: "locus",
  resource: memory | work | time, phase: "traversal" | "annotation_decode", requested, effective, used,
  remaining, actions, message}`,
  the amounts in the knob's unit (MiB — `used` rounded up, `remaining` down —, work units, ms; `requested` is the
  request's value where the server clamped it), `actions` the levers (`raise_memory_budget`, `use_graphlet`,
  `drop_sequences`, `raise_work_budget`, `raise_time_budget`, `continue_from_leaves`; `drop_coordinates` (feature
  level 6) right after `use_graphlet` and `drop_sequences` wherever those are offered — memory stops, depth-0
  failures included — when the request has `output.coordinates: true`, whose block or null form is part of the
  account; on a seed a budget failed,
  which has no leaves, the levers on the seed instead of `continue_from_leaves`: `lower_max_seed_labels` for a
  derived set, `name_fewer_labels` for a named one, `lower_max_labels_per_node` in annotate mode, and
  `shorten_seed` for a work budget the seed phase spent). A stop from outside the walk (§6.8) is stated the same
  way with `scope: "attempt"`, `resource: "cancelled" | "attempt_deadline"`, the attempt's bound (`requested`,
  `effective`) and elapsed time (`used`) in whole ms, and `actions: ["continue_from_leaves"]` on a walked seed,
  `["retry_attempt"]` on a seed failed in its seed phase (`phase: "traversal"`) or never started (`phase:
  "not_started"`); its message says when the walk stopped on the attempt's clock and that no budget ran out.
  The MGT `Q` record carries
  it (§7.5.2). It is a property of the walk, the same in every detail; a memory stop itself depends on the
  requested detail (the output is charged), so the same request can stop at another depth in another detail.
  `phase` is `annotation_decode` when a budget-aware annotation row did not fit the memory the request had left
  (§6.8; memory only), and `traversal` otherwise: finalisation and serialisation are reserved at admission. A
  memory stop names its cause in its message (review of stage 3): **a row** (`annotation_decode`) — its standalone
  demand, or that its read alone with its row-diff dependency rows did not fit, against what was left beside the
  account, the level's lists and the rows read before it (a seed-phase failure: beside the seed and what the
  seed phase held — for the derivation the window's earlier rows; a root's: beside the depth-0 state built so
  far); **the labels a row names first** (annotate mode, `traversal`): their number and bytes, the row fitting;
  on a format whose reads are not budget-aware, **the labels a level's rows named** (annotate mode,
  `traversal`), charged after the read: their number and bytes and the account they put over the budget, its
  need stated as a lower bound (*review of 2026-10-06, U03-03*; stated as a refused head before);
  **the level's own lists** (`traversal`): the key and successor lists that left none of the budget for the
  level's read. For a refused read `used` is what the walk held (its account and what the level's fetch held
  beside it), `remaining` the budget minus it, and the `walk_domain`'s `observed` what admitting the refused row
  needed beside what was held, in MiB rounded up — its demand, or, where its read alone was refused, the least
  that read was seen to need (the effect then says "at least", also for a level stop, whose demand is always
  stated as the least it needed); either way more than the budget. The levers of a
  row read fewer rows: `more_selective_seed`, and in annotate mode `label_constrained_query` (annotate mode reads
  every node's row); a level stop also offers `raise_memory_budget`, `use_graphlet`, `drop_sequences` (they
  shrink the account the read competes with) and `continue_from_leaves`; a seed-phase failure
  `raise_memory_budget` and the levers above (the detail does not change what a seed-phase read needs); an
  annotate root's the same, plus `use_graphlet`, `drop_sequences` and `lower_max_labels_per_node` when the other
  arm's root, already in the account, is what the row would fit without. A stop by labels offers
  `raise_memory_budget`, `use_graphlet`, `drop_sequences`, `lower_max_labels_per_node`, `label_constrained_query`
  and `continue_from_leaves`; a stop by the level's lists `raise_memory_budget`, `use_graphlet`,
  `drop_sequences` and `continue_from_leaves`. Every memory `walk_domain` observed made after the caches'
  allotments were charged (a head, a level's read or lists, a root's read, a depth-0 state) is the smallest budget
  that admits the demand with that budget's own allotments (§6.8), and its effect says so; a seed-phase read's,
  made before them, is its demand. A depth-0 failure whose root row or labels were not built because
  the budget was reached (§6.8) states its need, and its `walk_domain` observed, as a lower bound ("at least").
  **A lower bound is no promise** (review of stage 3, answer 2): every "at least" — a depth-0 failure's message
  and effect, a refused read's effect at the seed or at a level, a level's lists — also says that raising the
  budget to that value may still fail, because the roots, labels, rows or later state not read or built yet may
  need more.
  `lower_max_seed_labels` and `name_fewer_labels` are not offered for a
  decode stop: on a row-diff annotation every read decodes a whole row, whatever labels are named. A work
  stop's message states how far `used` can exceed the budget, with its number: at most what was charged since
  the previous comparison, and the most its seed charged between two comparisons (one indivisible charge: a
  fetch call's rows with their coordinates — on a budget-aware annotation each with its row-diff dependency rows,
  and the message gives the weights —, a label-state scan, or the roots' rows with the end of the seed phase,
  charged as one so that the result complete to 0 bp is delivered; only the seed phase can fail, at a
  comparison finding it at least `W` units over budget; §6.8); on a row-diff annotation without the budget-aware
  path it says that dependency rows are not counted. A seed the work budget failed states the same in its
  message: the comparison every `W` units that failed it, and the most its seed phase charged between two
  comparisons (review of the stage-2 recheck, P2). The texts are free text within MGT v1 (the `Q` message, the
  `K` effects): each message fits 1,024 characters and each effect 640 (`GraphletStage3Decode.StatementsFitTheirWidths`). A refusal injected by the test hooks
  (`WalkerHooks::deny`, and `WalkerHooks::deny_decode` for an annotation read, C++ callers only) is reported as a
  memory stop — with `limit` `"unlimited"` when no
  budget is set, so that it stays representable in `K` and `Q` — whose message and `walk_domain` effect say that
  the refusal was injected and that no budget caused it, and whose only action is `continue_from_leaves`.
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
  seeds are traversed. `outcome.walks` (or the presence of `error`) is how a client tells the two apart. Under
  `bounds.max_memory_mb` the limitations also state `memory_bound_soft` (§7.0), whatever the cause: the
  derivation decoded rows no admission charged, and `observed` is what they held beyond the budget.
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
  particular displayed path**, and **a contiguous source occurrence** (trace or source validation). A run's
  `route_bp` (MGT `R`) is the **earliest** merge at which its lineage entered through a non-first parent, anywhere
  on the lineage's own route; a leaf label's `route_bp` (`end_labels[]`, MGT `T`) is the **latest** such merge
  the entry came through. Neither is the start of displayed support: after a lineage entered through non-first
  parents at two merges d1 < d2, the run keeps d1 while the displayed path spells another route's bases up to d2
  (reproduced on SRA: a run with `route_bp` 61 whose displayed support starts at 113). Displayed support of a
  run is `[evidence_from, to_bp)`, with `evidence_from` derived per run and displayed path from the merge
  partitions — the latest merge on the displayed first-parent chain that the lineage entered through a
  non-first parent, at least the run's `from_bp` (`DESIGN-traverse-graphlet.md` §5.1, `Claim.evidence_from` in
  the library); at a leaf the merge part equals the leaf label's `route_bp`. `route_bp == 0` on a run means no
  such merge on its route, so the whole run is displayed support. One exception to the definition, on the
  conservative side: at a merge of three or more parents the walker stamps an entry that beat an earlier
  parent before a later parent beats it, so the run that merge then closes (code `m`, `to_bp == route_bp`)
  carries the merge although it is not on the run's route; its displayed support is still `[evidence_from,
  to_bp)` (on the offline real data: 6 of 73,575 runs, all whole-run support). A reader counting `route_bp: 0`
  runs then undercounts, never overstates; stamping only the entries that survive the whole merge would change
  `R` in unbudgeted output. So "N of M labels share this exact displayed flank to x bp" counts runs with
  `entered_by: seed`, `from_bp: 0` **and `route_bp: 0`** — the first two alone establish support along *some*
  route, not along the one shown.
- **Continuation.** Every leaf ended by `max_extension_bp`, a resource reason or `beam_pruned` carries
  `continuation: {sequence, labels, loss_used, branches_used}`: the last `max(k, bp since the last switch)` bases
  up to `continuation_bp`, natural orientation, with the labels covering that whole tail. It is valid `/traverse`
  input — a fresh seed: edge-reuse and branch state reset, and every label restarts at loss 0. `loss_used` is the
  SMALLEST terminal loss of the labels: for no route of the continuation to exceed its original budget, reduce
  `loss_budget` by the LARGEST terminal loss of the continued labels (their `T` records), as the library's
  `next_request()` does; one request carries one loss budget and no per-label starting loss, so this is
  conservative for the labels that ended at a lower loss. The branch allowance is reduced the same way (review
  of the stage-2 recheck, design answer 3: a continuation from 40 branched on to 100 where one uninterrupted
  walk stopped at 46): `next_request()` lowers `branching.max_label_branches` by the LARGEST terminal branch
  count of the continued labels (`"unlimited"` stays unlimited), states it (`branch_budget`, and a note where
  it is conservative for lineages that used fewer), and keeps the original only when asked
  (`reset_branches=True`, stated). What the guarantee covers is the cumulative loss and branch count of each
  continued lineage under unchanged costs; edge reuse is not carried over, so a continuation is a new traversal,
  not one uninterrupted walk. The continuation's labels become the seed labels: `labels.extra` must not repeat
  them, and keeps only the other permitted labels that a chain of switches from a seed label enters within the
  reduced budget (the server's rule, §5); every label left out is listed per continued walk, and the note on it
  says when one uninterrupted walk could still have entered it from a continued label at a lower loss — by a
  chain through any label of the retrieval's pool, the kept extra labels included (review of the stage-3 fixes,
  P3: chains through a kept label were not searched, understating what the continuation may miss). A label alive at the
  leaf that does not cover the whole tail (it switched in within the last bases) is not among the
  continuation's labels, so the continuation does not seed it and may lack its lineage; the library states each
  such label. A continuation's
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
  parent): per parent, the labels whose kept entry came through it (empty lists in `annotate` mode). **The
  first parent of a merge is the one carried by the most labels** *(the owner's decision R21 (4), feature level
  6)*: in `constrain` mode the labels whose lineages its head brings into the merge node, in `annotate` mode the
  fewest labels present at a node of the parent's own segment (every head at the node holds the node's own
  labels; a merged parent holds the union of its parents' labels, so after nested merges the walk displayed
  upstream can be carried by fewer, DESIGN §11); ties keep the arrival order, as every parent was before.
  **Except in `tree` and `full` detail under a memory budget**, where a merge keeps the arrival order: there
  every head reserves the delivery of its displayed segment chain, so the displayed parent would move the memory
  stop, and the rule changes the display only, never the depth a budget certifies (`graphlet` and `summary`
  charge no chain and follow the rule under every budget). The displayed walk — a path's `segments` chain and
  spelled bases, its `continuation`, and `route_bp` (what of that spelling a label does not carry, §13 note 20) —
  follows first parents, so it passes a merge through its majority parent. Which labels reach a leaf, with which
  loss and branches, does not depend on the order; on a tie of loss and branches the first parent's lineage
  continues (its run, its `labels_via_parent` entry).
  `hairpin` events: `followed` (`hairpins: follow`, the child is present) vs skipped. `revisit` events:
  `length_bp` (the distance difference to the first arrival) and `same_distance` (a `keep`-mode join at one
  node and depth). `label_dict[]`: `column`, and `seq_id` for header labels — the `LabelRef` that joins a label
  across retrievals (a name alone can be ambiguous). `growth[].label_ends` is an object even when empty (was
  `null`).
- **Record coordinates** *(feature level 6; opt-in: `output.coordinates`, §5; `DESIGN-traverse-graphlet.md`
  §18)*. With `coordinates: true` every seed result, in every detail (a graphlet's JSON summary included,
  never its MGT body), carries either `"coordinates": {block}` or `"coordinates": null` with
  `"coordinates_reason"`: `"index has no coordinates"`, `"support kmer"` (annotate mode included), `"no
  traversal"` (a failed derivation, a seed a budget failed, a seed refused for a label name, a seed never
  started) or `"partial derivation"` (a `trace` seed walked with a set derived from part of it, D3), checked in
  that order. The block (keys sorted): `kind` — `record` (every label a header: positions within its record),
  `column` (every label a column: positions in the column's k-mer index space) or `mixed` (each label's kind
  says which) —, `k`, `max_occurrences` (the cap, an integer or `"unlimited"`), `complete` (false iff a list was
  cut or a run is a lower bound), `runs_lower_bound` (only when > 0), `seed` — per seed label, ascending:
  `{label, occurrences: [[s, e], …], occurrences_total?}`, its occurrences of the whole seed (e − s =
  `length_bp`) — and `arms` — per requested arm one entry per run, in `runs` order: `{run, label, from_bp, to_bp,
  occurrences, occurrences_total?, chains_ended?, lower_bound?}`. **Intervals**: 0-based, half-open, forward
  strand; an occurrence of a run is a chain of the label's coordinates live at the run's last node (chains only
  continue or die along a run; a split hands each to at most one child), c its k-mer coordinate there and L =
  `to_bp − from_bp`: `[c + k − L, c + k)` on the right arm, `[c, c + L)` on the left (L = 0: an empty interval at
  the seed boundary); the interval holds the run's own bases. Each list keeps its first `max_occurrences`
  occurrences by start and states `occurrences_total` when cut (the `coordinates` limitation, §7.0).
  `chains_ended` (only when > 0, decision C-N1): chains live on the run's lineage that continued on no followed
  path of it before its last node — a record end, a successor not followed, a continuation a switch took over;
  chains partitioned among a split's children are not ended, and a clone inherits its prefix's count.
  `lower_bound: true` (only when true; revision 1, decision C-N2): the run was entered by a switch into a label
  whose own lineage was live there, and that label's chains starting at the switch node were left out, so
  its count is a lower bound (the walk itself is unchanged). **Column positions**: record i's k-mer j is
  `offset_i + j` with `offset_{i+1} = offset_i + len_i − k + 1`, so a column interval numbers base p of record i
  as `offset_i + p`: **a record's last k − 1 bases share their numbers with the next record's first k − 1
  positions**, an interval there is attributed to one record only with the record lengths (not stated; C10),
  and column intervals of different records can overlap (`trace_record_boundaries` states that a column trace
  can cross records). **Assumed, not detectable from the index**: one strand per record (no `--fwd-and-reverse`
  build) and unique headers. Memory: each entry is charged once when its run is created, at min(chains, cap)
  occurrences, in the requested detail (no work charge: every coordinate was charged one unit with its row);
  budgeted requests with coordinates stop no deeper than without (`drop_coordinates`). **What bounds the block**
  (the probe's `coordinates.output_bound` refers here, §10.3): its occurrences are coordinate chains the walk
  read (every coordinate charged one work unit with its row), so a seed lists at most the coordinates it read,
  and at most min(chains, `max_coordinate_occurrences`) per run and seed label; under a memory budget —
  `bounds.max_memory_mb`, or the server's `max_memory_mb`, which a request without one runs under when it is not
  0 (§10.3) — every entry and occurrence is charged when its run is created (at the count it starts with, in the
  requested detail), so the block is within the budget; under a work budget (`bounds.max_work_units`, or the
  server's `max_work_units`) the coordinates read are within it; with neither (`--traverse-max-memory-mb` and
  `--traverse-max-work-units` 0, their default) and `"unlimited"` nothing else bounds the block: a client sets
  a budget, keeps the cap or drops the coordinates.
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

**`usage`** (a request with `attempt_id` only, §5; `DESIGN-traverse-graphlet.md` §14.1: "usage in every
response" — the ledger reconciles its reservation against it): in every response to such a request — 200
(complete, partial, failed seeds, cancelled), 400 and 500 after the request was read, and 503 when the bound was
reached while the response was built — but not a 409 (a duplicate id: it would be reconciled against the other
attempt) nor a 400 for a malformed id:

```json
"usage": {"attempt_id": "a-17", "budget_id": "b-3", "locus_id": "l-9", "server_instance": "9f3c0d1e2a4b5c6d",
          "not_after_ms": 1791137060000,
          "reason": "completed", "received_at": "2026-10-03T19:17:45.540Z",
          "stopped_at": "2026-10-03T19:17:48.565Z", "elapsed_ms": 3026, "bound_ms": 40000,
          "bound": {"seeds": 1, "time_budget_ms": 30000, "allowance_ms": 10000, "walk_until_ms": 35000,
                    "capped_by": null, "enforced": true},
          "seeds": {"requested": 1, "started": 1, "finished": 1, "abandoned": 0},
          "work_units": 47720760, "observed_max_uninterruptible_ms": 37,
          "memory": {"peak_admitted_bytes": 698981842, "soft_excess_bytes": null, "held_bound_bytes": null},
          "per_seed": [{"index": 0, "outcome": "partial", "stopped_by": null, "work_units": 47720760,
                        "work_seed": 0, "peak_admitted_bytes": 698981842, "final_bytes": 698857938,
                        "soft_excess_bytes": null, "refused_bytes": null, "elapsed_ms": 3021}]}
```

- `reason`: `completed` (no stop reached the walk), `cancelled`, `deadline` (the attempt's bound, in the walk or
  while the response was built), `error` (a 400/500). `client_gone` appears only in the attempt's state (no
  response is written then).
- Instants are ISO-8601 UTC with milliseconds; durations are whole ms, rounded up. `received_at` is when the
  server read the request's header (the attempt's clock), `stopped_at` when the walk stopped (the checkpoint
  that saw the stop, or the end of the last seed's walk), `elapsed_ms` when the block was built (just before
  the response was written; the attempt's state has the final value).
- `bound_ms` is the bound the server enforces (§6.8) and `bound` its derivation: `seeds` × `time_budget_ms`
  (the effective per-seed budget) + `allowance_ms`, `capped_by: "content_timeout"` when the HTTP server's cap
  applied, `walk_until_ms` where the seeds stop being walked — `bound_ms − max(allowance_ms / 2, the delivery
  reserve)`, which moves with the reserve (pass 5, §6.8): when it stopped the walk, its value then (the walk
  stopped at its first poll that read the clock after it, `stopped_at`); otherwise the lowest walk-until seen:
  the lowest that a poll reading the clock compared with — every 8th poll of a walk, every poll before a paced
  read's chunk or in the lookahead, and the poll before each seed —, which is `bound_ms − allowance_ms / 2` when
  no such poll saw the reserve above that floor; a lower one in force only between two such polls is not stated
  *(review of 2026-10-06, C20: "the lowest in force while a seed was walked" claimed every value; pinned with
  the server's stride by `GraphletAttempt.WalkUntilStatedIsTheLowestAClockReadingPollSaw`)* —, `enforced` (false
  in the CLI). A request that failed before it was parsed states `seeds: 0`,
  `time_budget_ms: null` and the cap.
- `not_after_ms` (pass 5): the request's own, echoed only when it gave one (§5).
- `observed_max_uninterruptible_ms` (pass 5; every response with `usage`): the longest single uninterruptible
  piece of the attempt's walks — a chunk of a split read, a read decoded whole because its deadline was far
  (§6.8, chunked deadlines), or, from feature level 6, a head piece (the walk between two readings of the clock
  for a stop, its reads excluded; the lookahead's chains between their polls included: the review of 2026-10-06,
  W3) — whole ms rounded up: how late a stop could have come in them. An observation of this attempt, not a bound
  and not a margin: the mapping of a seed's k-mers, the seed phase's own processing, a derivation's steps and a
  seed's finalisation are not in it (§6.8).
- Number types: every count, byte amount and duration in `usage` is an integer (ms rounded up), except
  `bound.time_budget_ms`, a number as the request's knob is.
- `seeds`: `requested`; `started`; `finished` — the seeds whose walk ended (complete, partial or failed),
  delivered or not; `abandoned` — the seeds whose walk the client's departure cut (neither finished nor
  delivered; visible in the attempt's state, since no response is written then). A seed never started is in
  none of the last three.
- `per_seed` (in responses; the attempt's retained state keeps the totals only): `outcome` — `complete |
  partial | failed | not_started` as the result's `outcome.walks`; `stopped_by` — the resource that stopped
  it: `memory | work | time | cancelled | attempt_deadline` (null: none, or a cap); `work_units` — the seed's
  charged work (its seed phase `work_seed` and both arms, the account a work budget is compared with; computed
  whether or not a budget is set: deterministic logical work, §6.8; a derivation found `too_wide` includes the
  window it decoded); `peak_admitted_bytes` — the modelled account's peak, which under `bounds.max_memory_mb` is
  at most the budget (an account past it was refused — the seed failed at depth 0 or the walk stopped — or is
  held beyond it as the soft excess); `final_bytes` — what the seed's result holds until the response is
  written: the account when its walk ended for a walked result, the failed result as priced (its fixed part
  and its echo of `seed_id`, §5) for a failed, refused or never started one; `soft_excess_bytes` — what the
  result's `memory_bound_soft` states, in bytes (the walk's observed excess, and for a failed or never started
  result its echo of `seed_id` beyond the budget), null without `bounds.max_memory_mb`; `refused_bytes` — the
  demand the memory budget refused when it stopped or failed the seed (bytes at that budget: a head, a read —
  the least it was seen to need where that read alone was refused —, or the depth-0 state; its knob value is
  the `walk_domain`'s observed), null when no memory budget did; `elapsed_ms` — the seed's walk and its
  result's building. A seed failed in its seed phase reports what that phase consumed; a request failed by one
  seed (400) reports that seed too.
- Totals: `work_units` sums the seeds; `memory.peak_admitted_bytes` = max over seeds i of (the `final_bytes` of
  the seeds before i, whose results are held until the response is written, + i's peak, for a walked seed), and
  at least all results' `final_bytes` together; `memory.soft_excess_bytes` is the largest per-seed excess (each
  seed's over its own budget), null without a memory budget; `memory.held_bound_bytes` — an upper bound of what
  the walks and the results held at once, **as observed** (what `memory_bound_soft` lists as not observed is not
  in it), under `bounds.max_memory_mb` only: max over walked seeds i of (the `final_bytes` before i + the budget
  + i's soft excess), and at least all results together; null without a memory budget. The modelled account is
  **no bound of what was held**: without a memory budget it leaves out the label caches (up to 1,000,000 rows),
  the lookahead and a level's decoded rows (review of the stage-4 backend, F2: 0.79 MB stated for a walk that
  raised the server's RSS by about 111 MB), and under one the rows within the budget beside the account are
  in neither number — read `held_bound_bytes`, not `peak_admitted_bytes`, as the request's memory. The block
  lies outside the per-seed budgets (like the strategy echo).

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

**MGT v1 is frozen (2026-10-02).** This section is the wire contract; any change to a record, a field or the
grammar is MGT v2. **The `[a-z_]` token fields are open token sets** *(decision X-C12, amending the freeze text;
the precedent is `memory_bound_soft`, 28c51c31)*: a K record's `kind` and its extra field names, and a Q record's
`scope`, `resource`, `phase` and `actions`, can take values a later feature level adds — level 6 adds the kind
`coordinates` with the extra field `lists_cut` and the action `drop_coordinates` (C12) — and a reader keeps a
value it does not know, uninterpreted. Record coordinates (§7.1) are **not** in the body: only their limitation's K
record is, when a list was cut (`K * coordinates output.max_coordinate_occurrences i:<cap> i:<largest> *
lists_cut=i:<n> <effect>`, Z one higher); the block travels in the JSON summary (§7.5.5), with the J line of a
standalone document.

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
by the library after `H`, never by the server: the envelope with `results` reduced to this seed's summary — the
graphlet string removed — and *(pass 5)* an attempt's `usage` reduced to its totals plus this seed's `per_seed`
entry, the same rule in `save()`, `dump()`, `standalone_text()` and the store; `view` and `derived_from` are the
library's own top-level names in it (a view's selector, a derived graphlet's source), and a response that
carries either is refused by `from_response` and `standalone_text` rather than misread. A file saved by an
earlier library (feature level 2) whose `J` kept every seed's `per_seed` entry stays valid and loads as it is
(`load()`, `parse()`), but it is not in the canonical form of this rule: `is_canonical()` is false for it,
`dump(parse(text))` differs from it, and `save()` rewrites it reduced. The `J` record is the library's, not
MGT v1's: no record, field or token of the format changed.)

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
    which need not be UTF-8: a seed whose label or dropped-label names include one that is not is refused
    (§6.1 step 4, cause unrepresentable_label_name) and has no body, so every name the body carries is the
    index's own, verbatim — a replaced name (U+FFFD, as servers before this rule wrote it) could be another
    label's; the writer refuses a name that is not UTF-8 as a backstop
O <walks c|p|f> <branch_diagnostics c|x> <label_evidence c|l|q> <delivery i|s|p>
    the per-seed outcome (§7.0); q = qualified (DESIGN §14)
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
    values when the run ended or was closed (§7.1). route_bp = the run's earliest non-first-parent merge
    stamp (§7.1), not where its displayed support begins (that is derived: evidence_from); a run closed by
    a merge of three or more parents can carry that merge itself (to_bp == route_bp, §7.1)
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
 "label_mode", "duplicate"?, "outcome", "limitations", "coordinates"?, "coordinates_reason"?,
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
  `q` (the design's `qualified`, §14, written under `trace_record_boundaries`, §7.0).
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

**Coordinates through `CoordToHeader`** *(efficiency pass)*: the coordinates of one sequence of a column occupy a
contiguous range, so a tuple row's coordinates of a column (sorted) are mapped by the runs they form
(`LabelOracle::map_coords`): once two coordinates in a row fell into one sequence or adjacent ones, a coordinate in
the range of the previous one's sequence maps by a subtraction (its last coordinate: one select,
`CoordToHeader::last_coord`), one in the next sequence by one select, any other as `map_single_coord` maps it (a
rank and a select); lists of fewer than 8 coordinates are mapped by `map_single_coord`. The header label of a
sequence is looked up once per run. The recorder (annotate), which needs only the sequences, maps one coordinate
per sequence (`CoordToHeader::sequence_range`: a rank and two selects, one select fewer than before) and skips the
rest of its run; a header discovery counts a sequence once per k-mer from the run mapper. The same `(seq_id,
local)` per coordinate as `map_single_coord`, so the same hits, lists and `timing.coords_mapped` for the same rows
(one per coordinate mapped; one per sequence for the recorder; the lookahead's radius cap, §8.3, maps fewer
rows). Measured on a column of 2,000 sequences
(`BM_TraversalCoordRuns` in `benchmarks/traversal/bench_traversal_measurements.cpp`): runs of 7 coordinates in one sequence (rRNA operons) at 22 ns a
coordinate against 42–46 by `map_single_coord`, one coordinate in each of consecutive sequences (a single-copy gene)
at 32–35 against 48, scattered ones at 40–42 against 39–42; on the mini refseq rows (2–6 coordinates a column)
unchanged. A per-coordinate rank and select were 3–7% of a refseq33m request, the remaining per-row cost once
selected-label decoding (3c) shrinks the decode.

**Permitted-range filtering of header labels** *(review of pass 5, part B; feature level 5)*. A query of named
header labels needs only its own sequences' coordinates: it holds, per column, the coordinate ranges `[first,
last]` of its requested sequences (from `CoordToHeader::last_coord`), and looks each coordinate of the column up
among them by a binary search — a coordinate outside them is dropped without being mapped, one inside maps to
`(label, coordinate − first)`. The run mapper above (rank and select per run) still maps the coordinates where every
sequence matters: the recorder (annotate), a header discovery and a derivation. The same hits and lists as mapping
every coordinate (tested on mini refseq against `map_coord` for every coordinate, both read paths, with and without
coordinates); `timing.coords_mapped` still counts every coordinate of the column examined. This is a filter applied
to the reconstructed row; pushing the permitted columns and ranges into the row's decoding is stage 3c, where the
ranges must be shifted along the row-diff path (a coordinate moves by one at each step), never applied as one
absolute interval at every dependency row.

### 8.3 Batching along unbranched runs

Walk structurally up to `batch_kmers` nodes while each node has exactly one non-`'$'` successor; map and fetch
labels for the chunk in one call; then apply the recurrence node by node and truncate at the first node where
the path ends (also at a cap, §6.8). At a structural branch, fetch all successors in one call. Unverified nodes
are never emitted. Results do not depend on `batch_kmers` (T23). *(Efficiency pass:)* a chain stops at the
radius — at most `min(batch_kmers, max_extension_bp − d − 1)` nodes from a level at depth d: a head at the radius
is ended, not expanded, so neither its successors nor the rows beyond it are ever asked for (before, a chain ran
`batch_kmers` nodes whatever the radius left: 1,047 rows read for 32 consumed on refseq33m) —, and the lookahead
is skipped once the seed's deadline has passed and paced by it at every depth, depth 0 included (§6.8, chunked
deadlines: no later level runs then, so nothing it would read can be consumed). Only physical work changes
— `timing`, including its counters `rows_fetched`, `tuple_rows_fetched`, `coords_mapped` and `cache_hits`, which
count the lookahead's reads (on the mini_refseq set 250 of 588 responses: `mini_batch3__batch_refuse`
`tuple_rows_fetched` 428 → 370, `coords_mapped` 430 → 372; `rows_fetched` in 321 of the 396 UHGG responses with
timing on) — and two counters of the `annotation` block that count lookahead work: **`direct_reads`** (on an
annotation read by single cells, `access_path: direct` — a column or BRWT annotation with at most 16 named column
labels —, the lookahead's reads are single-cell reads too, so a chain that ran past the radius made reads no walk
consumed: e.g. 568 → 520 on the `merge` CLI fixture; the row-diff indexes read rows and are unaffected), and
`keys_mapped` where an unbudgeted lookahead outgrew its 1,000,000 entries and its clearing counted the keys (none
of 2,184 real requests on row-diff indexes changed).

### 8.4 Graph handle, caches, threading

All reads are const and concurrent requests share the loaded `AnnotatedDBG` (as `/search` does). The request owns
a `TraversalContext` holding: the walk graph (for `primary`: a lightweight `CanonicalDBG` clone over the same
primary graph with a `NodeFirstCache` attached — the shared wrapper is never mutated), a `NodeFirstCache` for any
underlying `DBGSuccinct`, the matrix cross-casts, the `CoordToHeader` pointer, the label cache and counters.
Per-step costs to expect: `basic` DBGSuccinct right step O(σ) rank/select, left step 2 bwd with the cache;
`primary` first-visit step adds O(k) (`rc_index_range`, counted).

The reverse header index (header → `(column, seq_id)`, for resolving a header label by name) belongs to the
`CoordToHeader` itself: built while the index loads (§10.3; feature level 5 — before, on its first lookup, inside
that request's budget), shared by every request on the loaded index, and freed with it. It is never a process-wide cache keyed by the object's address: an object made later at a freed address
could be answered from another's index, which no partial fingerprint of the headers can rule out.

**The row-diff path cache** *(efficiency pass; feature level 4)*. On a row-diff annotation every read of a
`/traverse` request — the default reads and the budget-aware ones (§6.8), a level's fetch, the lookahead, a
derivation's window, a validation — keeps rows it reconstructs (until feature level 5 every row of every path; now
those its retention rule selects, below) in a cache of at most `--traverse-path-cache-mb` MiB (default 128, 0 = off; the CLI takes
the same flag), so that a later read's row-diff path stops at a cached row instead of decoding to its anchor
again. On refseq33m a fetch of one warm 23S k-mer cost one whole path decode (66–86 ms) and each further row of
the same path 2.1 ms, so a walk fetching a few keys per level paid 27–41 ms a row against 4.6 ms read in one call;
locally, warm, the same binary with the cache off and on, the cache cut the annotation fetch time of the real
requests 3.1-fold on UHGG (198 requests; elapsed 2.0-fold) and 3.05-fold on SRA (201; elapsed 1.65-fold) at the
default `batch_kmers`, 12.9-fold (UHGG; elapsed 6.9-fold) and 7.2-fold (the SRA 16S cells; elapsed 3.05-fold) at
`batch_kmers: 1`, 2.3- to 2.8-fold under memory budgets of 256 and 8 MiB and a work budget of 10⁶ units (UHGG);
on a 16S walk with five named UHGG columns (radius 1,000, 3,572 rows) fetch calls of 1, 3 and 30 k-mers took
0.202, 0.116 and 0.077 ms a row without the cache and 0.029, 0.028 and 0.023 with it (7.1-, 4.2- and 3.4-fold).
A row's content does not
depend on how it was reached, so **what a read returns never depends on the cache**, nor does anything a budget
charges or admits: a budget-aware read whose path stops at a cached row takes the row's path aggregates from the
cache (`PathAggregates`: the whole path's length, entries, stored bytes, scratch and reconstruction peak), so every
row's `RowCost` — its work (8 per dependency row, 1 per entry) and its standalone demand — is that of its whole path,
the cache only shortening the decoding (and the transient charges of a call; a row read alone never needs more
than without the cache, so where a fetch stops is unchanged). The cache is bounded in bytes (each row as its exact
copy plus 640 B for its share of the tables) in two generations: a row is **admitted before it is copied** — its
stored size computed from the row, the shared bound read and the older entries evicted first, so that a row the
cache refuses (more than half the bound) is never copied and the copy of one it keeps is made within the bound
*(review of levels 4–5, finding 1: copying first, a cache bounded at 1 KiB allocated 4,210,688 bytes for a row of
1,048,576 columns and then refused it, its bytes, peak and kept rows all 0; and a cache full of two 1.5 MiB rows
under 4 MiB held 4.5 MiB while it copied a third)* —; rows go to the current one, which becomes the
older one when it fills half the bound (the older one is dropped then); a hit in the older one moves the row to
the current one; a generation dropped frees its table (a cleared hash map keeps its bucket array, which after
many narrow rows held up to 74.7 MiB of heap for 63.9 MiB accounted at a 64 MiB bound; review of the efficiency
pass: 62.6 MiB now), and each generation's count moves with its table at a rotation (so that a trim to a shrunk
shared bound sees the older generation too). Without a memory budget it is the request's, kept from seed to seed (only physical work depends on it).
**Under `bounds.max_memory_mb`** it is each seed's: empty and off until the depth-0 state is admitted — the seed
phase (derivation, validation) and an annotate root's read decode as before, so that what their refusals state is
as before too (a root read alone refused states a lower bound, a root whose standalone demand did not fit states
that demand; with the cache on, the left root's path could hold the right root's row and turn the one into the
other: review of the efficiency pass) —, then within what the label cache leaves of the label cache's allotment
(`min(budget / 4, 64 MiB)`, held by the account): the path cache is trimmed before the label cache grows and
re-reads the room at every insert, so that both stay within the allotment while the label cache evicts exactly as
before — no stop, admission, result counter or soft overshoot changes —, and emptied when the seed ends. **The
lookahead's reads under a memory budget** (§8.3) are admitted as without the cache: the lookahead reads ahead until
a run does not fit, and a run cut short by cached rows holds less, so with the cache it read rows no level asked
for (UHGG 16S, annotate mode, 8 MiB: 75k rows warmed against 52–55k, 13–27% more instructions than without the
cache). Each of its reads therefore also charges, from its trace until its stored rows are read, what decoding the
rows the cache spared it would have held — per anchor (paths to one anchor share their tail) the rows of its
longest cached path beyond the cached row, in the trace's containers and per-visit arrays, and that path's stored
rows (`RowDiffCache::admit_as_uncached`): exact for a chain, more than the union where paths of one anchor join
another read's (the lookahead then reads less than without the cache), less where cached paths branch. Measured
warm against the previous build (and the same build with the cache off), instructions: the 198 UHGG requests at
8 MiB 37.8 G against 71.4 (69.7), `uhgg_16s__annotate_merge` among them 2.35 against 3.48 (13–17% more than the
previous build before this rule), the 15 UHGG 16S cells at 6, 10, 12 and 16 MiB 4.7, 6.7, 5.1 and 7.9 against
9.4, 16.0, 12.8 and 12.0, the 201 SRA requests at 8 MiB 76.3 against 106.6 (105.5) and the 147 mini_refseq ones
under 4 MiB 17.1 against 19.5 (19.4); no request 5% above the previous build, every body identical. Stated in the
probe as `decode_cache: {path_cache_mb, rule}`.
`annotation_fetch_ms` and the elapsed times are what shrinks; the physical counters in `timing` count the rows the
reads return, and under a memory budget the cache can change how many that is (`rows_fetched`,
`tuple_rows_fetched`, `cache_hits`, `coords_mapped`): the lookahead's reads stop about where they stopped without
it, not exactly, and a level read the cache lets through in one piece decodes the rows after a refused row of the
piece (§6.8: a stop's position does not depend on it). Without a memory budget they are the same with and without
the cache. With several traversals at once each holds its own cache (up to `path_cache_mb` each without a memory
budget).

**Retention** *(review of feature level 4, R10; feature level 5)*. Keeping every row of every path copied each row
of a long path into the cache. On a coordinate annotation with wide rows (refseq33m 23S: about 28,000 coordinates
in 11,000 columns, one allocation per column with more than two of them) a first read, whose paths no later read
meets, then cost several times the read without the cache, and a walk reading one row per call churned the cache's
generations so that every call decoded its whole path again. A read now keeps (`RowDiffCache::keeps`): every row
it was asked for (later paths meet them: a predecessor's path runs through its successor; a repeated read is a hit);
the first 8 rows after each of them on its path (a forward walk reads these next); every row whose distance to its
anchor is a multiple of 16, anchors included (a **checkpoint**: a later path through rows of this one stops within
15 rows of a kept row instead of running to the anchor; the distance is a property of the row, so which rows are
checkpoints does not depend on the calls); and every **narrow** row, whose copy holds less than 4,096 bytes (copying
it costs little beside the descent that decoded it: on UHGG, binary rows of 10-200 bytes, the rule keeps every
row, the same stored rows read and hits as before). A call of n rows whose paths read s stored rows thus keeps at
most n × 11 + s / 16 rows that are not narrow, where it kept s; `timing.path_cache` counts what it copies
(`rows_kept`, `bytes_kept`), the stored rows read (`stored_rows_read`), the hits and the peak. A tuple row is kept
**flat** (its columns, the end of each column's coordinates, the coordinates: three buffers instead of one per
column with more than two coordinates), its bytes so counted against the bound; a hit rebuilds the row exactly and
a budget-aware read is charged its copy as before. Measured on a synthetic coordinate index (8,000 labels each
carrying three copies of a 3 kb core with SNPs; 24,000 coordinates a row; paths of up to 100 and of up to 1,000
rows), medians of 5, against the cache off and the pass-5 binary: no request slower beyond noise; isolated first
reads 1.1-2.9 times faster than with the cache off (1.6-6.2 times faster than keeping every row), walks with
`batch_kmers: 1` 6-15 times faster than with the cache off (keeping every row was 4.9 times faster than off on paths
of 100, and 7.7 times slower on paths of 1,000: 16.9 s against 2.2 s), time-bounded walks cut with the cache off
reaching 1.5-2.0 times its rows; on UHGG the same work as keeping every row. The seed phase (a validation, a
derivation) keeps reading through the cache, as decided by these measurements: the isolated first reads above are
seed phases. The rule changes what is kept only, never a row, a charge or an admission
(tested across rules, capacities, batches, chunks and cold and warm seeds, T52).

Under a memory budget the label caches get a fixed allotment of it (§5) and **evict wholesale** when they exceed
it, as without one (a walk moves forward, so recency is not worth tracking; evicting the oldest rows first was
measured in the review of 2026-10-06 and starved the path cache, which holds what the label cache leaves of
the allotment: §6.8). A row evicted and needed again is
decoded again and charged again when the fetch returns it, like every row (work counts the rows the walk's
fetches return, §6.8), so the refetching is bounded by the work budget and the deadline, and the memory it
holds beyond the allotment is observed as `memory_bound_soft`. After a budget stop or a refused admission an annotate dictionary is compacted
to the labels the result records (§6.9); a cap, and a time stop without a budget, keep it as it was.

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

`POST /resolve`, `POST /traverse`, `GET /traverse/capabilities`, `GET /capabilities`, `POST /traverse/cancel` and
`GET /traverse/attempt/{attempt_id}` in `server.cpp` (the attempts in `traverse_attempts.{hpp,cpp}`; the pattern
search's routes, `POST /pattern` and `GET /pattern/capabilities`, are `SPEC-pattern-search.md`'s).

- **Shutdown** *(level 6's commit)*: on SIGTERM or SIGINT the server exits within seconds, status 0. Before, it
  had no handler, and in a container, where it is PID 1, the kernel drops a signal with the default action: each
  `docker stop` waited its whole timeout (120 s on staging) and killed it. Now nothing new is served (a
  `/traverse` or `/resolve` finds its client gone at its first check), every traversal in flight — with or
  without `attempt_id` — is stopped at its next poll as for a client gone, so **nothing is written for it** and
  its connection is closed (a partial result would state a cause that is not true; a ledger treats the attempt
  as unanswered: `not_after_ms`, `bound_ms`), then, once they returned or after 3 s, the HTTP server stops (a
  response still being sent is cut short of its `Content-Length`, never a shorter complete one), and a handler
  that has not returned 3 s later (a `/search`, an uninterruptible read) does not hold the process. The process
  exits once the server stopped, whatever still runs — in single-index mode the server listens at once and answers
  503 while its index loads, and a signal then does not wait for the load. Before the server accepts connections
  the process exits at once; a second signal ends it immediately.
  Attempts' states, tombstones and holds live in the process and end with it (§10.3, the release rule).

- **Capabilities: which route carries which fields.** (`/stats` carries none of them: it states the graph and
  annotation statistics only, unchanged.)
  - **Every `/resolve` and `/traverse` response**, `capabilities`: `schema_version` (the request schema the server
    accepts, `strategy.schema_version`; 1), `feature_level`, `k`, `regime`, `alphabet`, `num_labels`,
    `has_coordinates`, `has_coord_to_header`, `supports_trace`, `cost_models_available`, `label_modes`,
    `direct_access`, `release`, `graphlet_format` (1: the MGT version `detail: graphlet` writes, §7.5),
    `detail_levels` (`["summary", "tree", "full", "graphlet"]`) and the index identity (`index_ns`, `index_fp`,
    `index_meta_fp`, below); `/traverse` responses also carry `algorithm_version` at the top level.
  - **`GET /traverse/capabilities`** (the probe; per graph in multi-graph mode, below): all of the above for its
    index, plus the server maxima (`max_time_ms`, `max_seeds`, `max_seed_bp`, `max_seed_labels`,
    `max_query_bp`; from feature level 4 also `max_memory_mb` and `max_work_units`, the server's maxima of the
    budgets, 0 when off, below), the request budgets it accepts (`budgets: ["max_memory_mb", "max_work_units"]`),
    `work_check_interval` (`W` of §6.8, in work units), `work_bound` (a reference to §6.8, work units: what
    `bounds.max_work_units` counts and how far a work stop can exceed it — deterministic logical work, not
    measured decode effort, with its weights, the dependency weighting of a budget-aware row included; what was
    charged since the previous comparison, one indivisible charge such as a fetch call's rows decoded whole, with
    each stop's message stating the most its seed charged between two comparisons; the seed phase failing at a
    comparison finding it at least `work_check_interval` units over budget; the physical decode counters in
    `timing`), `memory_bound: "soft"` (stage 3 charges budget-aware annotation reads, §6.8, but
    what is held beyond the admitted account is still soft, §7.0 `memory_bound_soft`: no hard request-wide memory
    bound), `attempts` (below), from feature level 6 `coordinates` (record coordinates, §7.1: `supported` —
    whether this index reports them, which is `supports_trace` —, `knob` `"output.coordinates"`, `cap_knob`
    `"output.max_coordinate_occurrences"`, `max_occurrences_default` (16), `kinds` (the block kinds this index can
    report: `record` and `mixed` need a `CoordToHeader`; `[]` when not supported), `limitation`
    `"coordinates"`, `action` `"drop_coordinates"`, `output_bound` — a reference to §7.1, what bounds the block:
    the coordinates read, min(chains, cap) per list, the memory budget that charges each entry at its run's
    creation, the work budget that charged each coordinate with its row; nothing else with neither and
    `"unlimited"` (`--traverse-max-memory-mb` is 0 by default) — and `rule`, a reference to §7.1, record
    coordinates: positions, the interval rule, the column numbering in which a record's last k − 1 bases share
    their numbers with the next record's first positions, `chains_ended`, `lower_bound`, the assumptions; only
    here: the per-request capabilities change only in their `feature_level`), `content_encodings` (`["gzip",
    "deflate"]`), `pattern` (the pattern search's block with `details`, `SPEC-pattern-search.md` §23) and, from
    feature level 3,
    `algorithm_version`, `compression_level` (the zlib level of the traversal routes' compressed bodies) and
    `deadline_check` (below), from feature level 4 `decode_cache: {path_cache_mb, rule}` (the row-diff path cache
    of the reads; `rule` a reference to §8.4); in multi-graph mode also `graph` and `graph_path`, the pair it
    describes.
  - **`GET /capabilities`** (feature level 3; the owner's server-wide document, `DESIGN-traverse-graphlet.md`
    §17.1), small and index-free, answered while a single index still loads (`ready: false`; the routes that read
    the index answer 503 then, the cancel and state routes do not):

    ```json
    {"algorithm_version": "traverse-0.2", "attempts": {"…": "as the probe's"}, "compression_level": 1,
     "content_encodings": ["gzip", "deflate"],
     "deadline_check": {"chunk_target_ms": 50, "max_uninterruptible_ms": null,
                        "observed_max_uninterruptible_ms": 37, "rule": "…"},
     "feature_level": 3, "features": ["search", "align", "resolve", "traverse", "attempts"],
     "graphs": null, "mode": "single", "ready": true, "release": "",
     "routes": {"align": "POST /align", "attempt": "GET /traverse/attempt/{attempt_id}",
                "cancel": "POST /traverse/cancel", "capabilities": "GET /capabilities",
                "column_labels": "GET /column_labels", "resolve": "POST /resolve", "search": "POST /search",
                "stats": "GET /stats", "traverse": "POST /traverse",
                "traverse_capabilities": "GET /traverse/capabilities"},
     "schema_version": 1, "server_instance": "9f3c0d1e2a4b5c6d"}
    ```

    In multi-graph mode `mode` is `"multi"`, `graphs` the sorted list of the graph list's names, `align` is in
    neither `features` nor `routes` (it answers 400 there), and `routes.traverse_capabilities` is
    `"GET /traverse/capabilities?graph={name}[&graph_path={path}]"`. `server_instance` is the attempts' (below).
    Both modes list the pattern search: `pattern` in `features`, `routes.pattern` and
    `routes.pattern_capabilities` (`"GET /pattern/capabilities"`; in multi-graph mode
    `"GET /pattern/capabilities?graph={name}[&graph_path={path}]"`), and the `pattern` block
    (`SPEC-pattern-search.md` §10, §23, §24; a multi-graph server's without a graph). Both state `in_ram` and,
    in multi-graph mode, `graph_summary` (`null` in single-graph mode; below, multi-graph mode).
  - **`deadline_check`** (both capabilities routes, feature level 3): `chunk_target_ms` (integer,
    `--traverse-chunk-target-ms`, an integer in [0, 2⁵³ − 1]; 0: reads are not chunked), `max_uninterruptible_ms`
    (null: no bound on one row's decode exists before stage 3c, and it stays null after stage 3c-ii, whose
    checkpoints bound index operations, not time — decision 3c-N5, §6.8), `observed_max_uninterruptible_ms`
    (integer: the longest single uninterruptible piece of any `/traverse` in this process's lifetime — a chunk, a
    read decoded whole far from its deadline, a head piece (from feature level 6: the walk between two readings
    of the clock for a stop, the lookahead's chains between their polls included, W3) or a gap between two
    delivery checks — an observation, which leaves out the mapping of a seed's k-mers, the seed phase's own
    processing, a derivation's steps and a seed's finalisation, §6.8), `poll_stride` (integer, 8: one in this many
    of a walk's polls — the one before a head — reads the clock for an attempt's walk-until and bound, so a
    walk-until is seen up to `poll_stride` − 1 heads late, §6.8) and `rule` (a reference to §6.8, chunked
    deadlines, which states the rule with the factors 64 and 4, the first chunk of 8 rows, the lookahead's polls
    every 16 graph steps, the steps left unpolled, what the observation counts, and how late a cancel, a
    walk-until and the attempt's bound are seen). It describes `/traverse`'s reads; `/resolve`'s deadline is the
    `resolve` block's.
  - **`resolve`** (both capabilities routes; milestone 1b of `DESIGN-pattern-search.md`, no feature level, §4.5):
    `{"time_budget": {"accepted": true, "knob": "bounds.time_budget_ms", "default": null, "max_time_ms": 30000,
    "finalize_reserve_ms": 250, "finalize_reserve_configurable": false, "check_kmers": 4096, "check_labels": 4096,
    "stop_phases": ["rows", "support"], "rule": "…"}}` — `default` null: no deadline without the field, whatever
    the cap; `max_time_ms` the cap (`--traverse-max-time-ms`, at most 899 000, and 899 000 when the flag is 0),
    `finalize_reserve_ms` the reserve inside every budget (fixed in this build), `check_kmers` the explicit labels'
    hits between two readings, `check_labels` the labels or JSON objects between two readings of the answer's
    time, and `rule` a reference to §4.5, which states the rule: opt-in, the 400 and the clamp, where the deadline
    is read, the stop block and the prefix it answers, the selection not made and the 400s that hold whatever the
    stop, the reserve and its 503, and what is not polled. On `GET /capabilities` it is index-free, so stated while a single index loads
    too. Its presence on either route is how a client knows the field is accepted (a
    server without it refuses the field, 400). **Neither `resolve` nor `pattern`** (`DESIGN-pattern-search.md`
    §7.3) bumps `feature_level`: the level is in every `/resolve` and `/traverse` response, so a bump would
    change every answer, which the milestone's byte-identity rule forbids; the probe and `GET /capabilities` of
    such a build differ from the build before exactly by these blocks (and `/capabilities`' `pattern` feature
    and route).
  - **The rules are references.** Every rule of the two GET routes is a string naming the section of this
    SPEC that states it, `"SPEC-labeled-traversal-core.md section N"` with, where the section holds several rules,
    the rule's name after a comma (a colon separating a part of it), printable ASCII, for people, not parsed: the
    two documents a service's MCP tool returns in one piece have a ceiling of 32 KiB, which the rules' prose
    nearly filled. The numbers a client computes with are fields beside them (`chunk_target_ms`, `poll_stride`,
    `work_check_interval`, `path_cache_mb`, the attempts' and the reserve's numbers). The texts are prose and
    raise no `feature_level`; nor does `deadline_check.poll_stride`, which is in no `/resolve` or `/traverse`
    response (as `resolve` and `pattern` are not).

    | field | route | the section |
    |---|---|---|
    | `work_bound` | probe | §6.8, work units |
    | `decode_cache.rule` | probe | §8.4, the row-diff path cache |
    | `coordinates.rule` | probe | §7.1, record coordinates |
    | `coordinates.output_bound` | probe | §7.1, record coordinates: what bounds the block |
    | `deadline_check.rule` | both | §6.8, chunked deadlines |
    | `resolve.time_budget.rule` | both | §4.5 |
    | `attempts.bound` | both | §6.8, the attempt's bound |
    | `attempts.not_after` | both | §5, `not_after_ms` |
    | `attempts.instance` | both | §5, `expect_server_instance` |
    | `attempts.suppression` | both | §10.3, `POST /traverse/cancel`: the tombstone |
    | `attempts.release_rule` | both | §10.3, the release rule |
    | `attempts.delivery_reserve.rule` | both | §6.8, the delivery reserve |
    | `attempts.delivery_reserve.calibration` | both | §6.8, the delivery reserve: calibration |

  - **`feature_level`** — what the server offers beyond the base contract, monotonic and only ever extended, so a
    client states a feature as `feature_level >= n` (`DESIGN-traverse-graphlet.md` §17.1: each change that adds
    capabilities fields or routes raises it by one; `schema_version` stays the request schema version). Absent or
    1: the base contract through stage 3. **2**: attempts — `attempt_id` / `budget_id` / `locus_id`, the `usage`
    block, `POST /traverse/cancel`, `GET /traverse/attempt/{attempt_id}`, the enforced attempt bound — and the stop
    when the client is gone. **3** (pass 5): `not_after_ms` with its expired 409 and
    `attempts.clock_skew_allowance_ms`; the per-graph identity (the graph list's columns 4–5, checked at start-up)
    and `GET /traverse/capabilities?graph=`; `GET /capabilities`; `algorithm_version` in the capabilities;
    `attempts.hard_cap_ms` and `attempts.allowance_ms` as an integer; `deadline_check`, the chunked deadlines and
    `usage.observed_max_uninterruptible_ms`; `compression_level` and the delivery reserve
    (`attempts.delivery_reserve`). **4** (the efficiency pass): the server's maxima of the budgets
    (`--traverse-max-memory-mb`, `--traverse-max-work-units`; the probe's `max_memory_mb` and `max_work_units`;
    the clamps below), the row-diff path cache (`--traverse-path-cache-mb`; the probe's `decode_cache`, §8.4) and
    the delivery reserve's calibrated starting estimates (`attempts.delivery_reserve.calibration`; ratios 30 and
    50, `stop_ms` = `chunk_target_ms` + 950). Feature level 4 changes no response of a request without budgets
    beyond the `feature_level` digit, `annotation.direct_reads` on a direct-access annotation whose lookahead ran
    past the radius (§8.3), the `timing` block's values — its times, and its physical counters `rows_fetched`,
    `tuple_rows_fetched`, `coords_mapped` and `cache_hits`, which count the lookahead's reads the radius cap
    removed (§8.3) — and an attempt's walk-until (the reserve's estimates): the lookahead's radius cap, the single
    decode of a discovery, the coordinate runs and the path cache change only physical work (checked against the
    previous build, T51); the same holds for budgeted requests while the budgets' maxima are unset, where the
    path cache can also change those counters (§8.4). **5** (the review of pass 5): the tombstone's suppression —
    `POST /traverse/cancel` with `not_after_ms`, `suppressed_until_ms` / `covers_admission` /
    `covers_admission_reason` in every tombstone answer, the refused copy extending the tombstone,
    `--traverse-attempt-tombstone-max-s` and `attempts.tombstone_max_s`, `cancel_fields`, the `suppression`,
    `release_rule` and `instance` texts — and `expect_server_instance` with its `instance_mismatch` 409; bounded
    retention settings (retention 0: no tombstones, 429 `no_suppression`); and from the review of these fixes a
    finished attempt's id held until its `not_after_ms` + skew (the release rule's finished state), the hold
    live through `suppressed_until_ms` inclusive, and a refused copy's suppression judged against its own
    `not_after_ms` alone; the loader dependency inventory (a manifest must cover every file the pair loads, the
    `.bloom` included — until the owner's decision #17 of 2026-10-08 took the mask and the Bloom filter out of it,
    below —, with distinct base names and no optional file of the inventory the pair does not load;
    bundles grouped by their whole inventory; `traverse --index-inventory`); the path cache's retention rule and flat tuple rows (§8.4), the
    sequence header index built while the index loads, the header labels' range filter (§8.2), a `/resolve`'s
    rows decoded in bounded batches, a discovery accumulating its profiles as it reads (§4.2), the budgeted
    lookahead's runs in the walk's order (§8.3), the
    delivery checks every 64 KiB of a large token and the delivery gaps in `observed_max_uninterruptible_ms`
    (§6.8), and in `timing` the seed phase (`seed_phase_ms`, `seed_fetch_ms`, `label_resolve_ms`), the path cache's
    physical work (`path_cache`) and the deadline record (`deadline`, §6.8). Feature level 5 changes no response of
    a request without budgets beyond the digit, the `timing` block (new fields; `coords_mapped` counts the
    coordinates a header label's filter examined, as before) and the texts of the capabilities; budgeted requests
    keep their bytes too (checked against the previous build, T52). The server states 5 from `7cdde4c2`. The fixes
    of the review of levels 4–5 (T53) keep level 5: they change no request or response field, only the
    capabilities' `attempts.release_rule` and `attempts.instance` texts (a finished state is replay-safe only for
    an attempt pinned with `expect_server_instance`, §10.3), the values in `timing` (the deadline record's
    `setup` pieces and its last `head` piece, `seed_phase_ms` and `seed_fetch_ms` to the first checkpoint, §6.8)
    and physical work (the path cache admits a row before copying it, §8.4). **Builds `7aaee760`..`67bef367`
    are not level 5** although they state it: they carry D3 (below) and accept the first `output.coordinates`
    fields, so a client keyed on `feature_level` cannot tell them from a level-5 build — `dcc0cebd` answers
    `walks: failed` to a derived seed whose time budget runs out after part of it, `67bef367` `walks: partial`
    with `label_evidence: qualified`, both stating 5. None was deployed (staging ran `32959610`), and none may
    be *(review of 2026-10-06, X1)*. **6** (record coordinates,
    `DESIGN-traverse-graphlet.md` §18): `output.coordinates` and `output.max_coordinate_occurrences` (§5) with the
    `coordinates` block or null form per seed (§7.1), the `coordinates` limitation and its K record (§7.0, §7.5),
    the action `drop_coordinates`, the probe's `coordinates` block,
    `attempts.delivery_reserve.coordinate_account_per_text_byte` with the reserve's coordinate share (§6.8). Its capabilities differ from level 5's exactly in:
    `feature_level` (both GET routes and every response), `attempts.delivery_reserve.rule` and the new
    `attempts.delivery_reserve.coordinate_account_per_text_byte` (both GET routes), `deadline_check.rule` (both:
    `max_uninterruptible_ms` stays null, no later stage promises a bound, decision 3c-N5) and the probe's
    `coordinates` block. A request without coordinates is answered as at level 5 apart from the digit, the
    `timing` block and an attempt's walk-until, with three intended changes: **D3** (a derived seed whose time
    budget ran out after part of its k-mers is delivered as a partial walk, §7.0 `derivation`), the **first
    parent of a merge** (R21 (4), §7.1: the displayed spelling, chains, continuations, `route_bp`, the order of
    `parents` / `labels_via_parent` and on a tie which lineage continues, all `on_reconverge: merge`; not in
    `tree` and `full` detail under a memory budget, where a head reserves its displayed chain's delivery and a
    merge keeps the arrival order, so no budget's stop moves — T55, DESIGN §11), and the server's
    **shutdown on SIGTERM** (below). The server states 6 from this level's commit. **The fixes of the review of
    2026-10-06** (its P2 items; level 6 was deployed nowhere, so they are folded into it and `feature_level`
    stays 6) change: the capabilities texts `attempts.bound`, `attempts.not_after`, `attempts.instance` and
    `attempts.release_rule` — what an attempt runs past its bound (up to its next delivery check: the rest of
    the piece the bound fell into, the walk up to its next poll that reads the clock, the stopped seed's
    finalisation and the building up to that check) and the clock release's assumption that it ended (X2), the
    409 refusing a copy as a finished-state source and the `expired` 409 as a release ground (LRG-R1/R2), the
    refusal order (C24), the `expect_server_instance` pattern and its 400 (X4), the hold that refused copies
    extend (C30), the content-timeout cap (C16), the walk-until stated and the poll that reads the clock (C20) —,
    `attempts.delivery_reserve.rule` (the walk stops at its first poll that reads the clock after the
    walk-until, C20) and `deadline_check.rule` (the lookahead's polls, the steps left unpolled, what the
    observation counts, how late a cancel, a walk-until and the bound are seen: W3, D12, X2); the texts of the
    attempt routes' 404/409/429 errors that state the retention (C30: the 404 of an unknown id, the 409s of a
    duplicate and of a tombstoned id, the cancel's 429 `tombstones_full`), the duplicate's 409 also no longer
    saying "an attempt runs once" (D3: "is running, or retained or held after it finished, on this server
    (…): it is not run while so"), of an `attempt_deadline` statement (the arm's `walk_domain` effect,
    `resource_stop.message`, a not-started seed's `error`) and of the 503 at the bound (C16: they name the
    content-timeout cap); and three outputs. **W1**: a derivation's windows are 64 k-mers whatever
    `annotation.batch_kmers` is (§6.1), so a request with a derived seed and `batch_kmers` below 64 can change
    where its seed phase's work is compared — its outcome, its `resource_stop` and `largest_charge` under a work
    budget, its `work_seed` and `usage.work_units` (an early `no_carrier`), its `rows_requested` — and a D3
    seed's j, and, on a format whose reads are not budget-aware, its soft memory observation under a memory
    budget, since a later window now holds 64 rows where it held `batch_kmers` (`usage.memory.soft_excess_bytes`
    and `held_bound_bytes`, the per-seed `soft_excess_bytes`, `memory_bound_soft`'s `observed` where it crosses a
    MiB: on the reviewer's 5,000 bp trace seed at 1 MiB and `batch_kmers` 1, 393,256 and 1,441,832 bytes before,
    400,800 and 1,449,376 now, the values at 64 — the more accurate count of what is held); at `batch_kmers` 64
    (the default) and above nothing changes. **W3**: the lookahead's chains read the attempt's stop and the
    seed's deadline every 16 graph steps and before each chain's key mapping (§6.8), through the poll that reads
    the clock and records a passed walk-until with the attempt (paced or not), so an attempt that a cancel or
    its walk-until stops inside a lookahead stops sooner (its timing values, `usage.stopped_at`, the physical
    counters of reads it no longer makes, and where its walk is censored) — a walk no stop ends is unchanged;
    and `observed_max_uninterruptible_ms` (usage and `deadline_check`) counts head pieces too, so it reads
    larger. **W2** changes only time (a trace
    derivation under a memory budget is linear again). Every other request is answered as by `6897db99` apart
    from timing (checked on the CLI fixtures and the bench capabilities, T58). **The fixes of the review's P3
    server items** *(the owner's approval, 2026-10-06; folded into level 6 as corrections: `feature_level` stays 6, and a level-6
    build before them — `ea285c2e`, the build deployed to staging — states what they correct, so a client keyed on
    `feature_level` cannot tell the two apart)* change three outputs of requests that ask for nothing new, every
    other untimed byte being as at `ea285c2e` (T59). **C10**: a seed an attempt never started (`resource_stop.phase:
    "not_started"`) states `seed.labels_from_seed` by the rule every result of its request states it — false in
    annotate mode, which derives no set, true for a constrain seed that names no labels — where it stated
    `seeds[].labels` empty, true for every annotate seed, beside a walked seed's false in the same response.
    **W9**: on a format whose reads are not budget-aware, an annotate level whose new dictionary labels put the
    account over `bounds.max_memory_mb` (they are charged after the read, §6.8) is stated as a stop by those labels
    — `resource_stop.message` names them and their bytes, `actions` gain `lower_max_labels_per_node` and
    `label_constrained_query`, and the arm's `walk_domain` effect says "did not admit the dictionary labels" with
    "observed: at least" (in MGT, the `Q` record's actions and message and that `K` record's effect, and so the
    result's `graphlet_bytes`) — where it was stated as a refused head with an exact need (its `observed`,
    `used` and where the walk stops are unchanged). **C26**: a response whose text is 4 GiB or more, sent
    compressed (`Accept-Encoding` gzip or deflate), is the whole text compressed, where it was a well-formed
    stream of the text's first (size modulo 2^32) bytes answered 200 — on every route that compresses (`/search`,
    `/align` and `/resolve` as well as `/traverse`: the server's one compressor); a text below 4 GiB is compressed
    by the same calls as before, its bytes unchanged. The rest changes no untimed byte:
    **C8** a seed's stop latency is measured at its walk's end only, so `measured_stop_ms` (and the walk-until of
    later attempts) no longer takes up the time a result was built when the client left during it or a writer
    refused it (10.5 s recorded for a 400 ms stop; short attempts then walked nothing); **C27** with bytes waiting,
    a socket past `ESTABLISHED` — the server's own shutdown at its content timeout included — reads gone, so on
    Linux such a walk is stopped instead of computing on (nothing could be delivered either way), and a descriptor
    that holds no connection reads gone; **C9** delivering a response
    holds less beside its text (the peak per byte of text 1.52 with gzip and 2.57 without, against 2.60 and 3.51,
    §6.8), the bytes unchanged; **W11** a budget-aware lookahead ends at the first run the cache cannot keep beside its
    own (physical work: the `timing` counters `rows_fetched`, `tuple_rows_fetched`, `coords_mapped`,
    `cache_hits`, `path_cache` and the times; fewer decodes in most walks it changes, slightly more in some, and
    what a memory budget still adds is stated in §6.8; §8.3); **W16** `/resolve`'s explicit labels are scattered
    in one pass, its profiles unchanged, and its client checked every 4,096 k-mers of it (time only).
  - **`schema_version`** is the one request schema the server accepts: 1, with no compatibility window. A future
    change of the request schema adds `request_schema_versions: [..]` to the capabilities and keeps accepting 1 for
    a stated window; a `schema_version` above 1 without that list means a server whose requests a client written
    for version 1 cannot make. Clients key features on `feature_level`, never on `schema_version`.
  - **The index identity** (`DESIGN-traverse-graphlet.md` §3.1), also in every graphlet's `H` record:
    - `index_ns`: a name for humans and routing, not identity (`[A-Za-z0-9._-]+`): `--index-name NAME` for a single
      index, the graph list's fifth column per (graph, annotation) pair; `null` when unset. A manifest's own
      `index_ns` (metadata) is not read for it; a different one is logged.
    - `index_fp`: the identity — the lowercase hex sha256 of the index bundle's **manifest** file list
      (`--index-manifest FILE` for a single index, the graph list's fourth column per pair): a JSON object with
      `files: [{path, size, sha256}, …]` (every identity file the server loads: graph, annotation and the
      sidecars it reads beside them, never the graph's derived mask and Bloom filter, below), hashed as the lines `<path>\t<size>\t<sha256>\n` in ascending byte order of path; other keys
      (builder, inputs, an entry's `digest`) are metadata. The build hashes the files once; the server reads the
      digests and checks, before loading anything, that every identity file it loads is listed (by base name) with its
      size — the **loader dependency inventory** of the pair *(review of pass 5, findings 2 and 3; one C++
      definition, `index_load_inventory`, which `metagraph traverse --index-inventory -i GRAPH -a ANNOTATION`
      prints as JSON and `scripts/traversal/index_manifest.py` mirrors, an integration test comparing the two on
      every sidecar kind)*: the graph; the annotation;
      `<graph>.anchors` and `<graph>.rd_succ` for a `.row_diff.annodbg` (required: the loader exits without them);
      `<annotation without .<type>.annodbg>.seqs` for a coordinate annotation unless `--no-coord-mapping` — each
      path derived from the LISTED spelling, as the loaders derive it (beside a symlink, not beside its target).
      *(The owner's decision #17 of 2026-10-08.)* The graph's **derived data** is **not part of `index_fp`**:
      `<graph without .dbg>.edgemask` (the dummy-edge mask, loaded when it opens) and `<graph without .dbg>.bloom`
      (the Bloom filter, loaded when it exists and the mask was loaded), both computed from the graph alone. They
      decide which counts are exact (`/pattern`: exact with the mask, bounds and an estimate without,
      `SPEC-pattern-search.md` §18) and how fast k-mers are looked up, never what an exact answer is; so adding,
      removing or rebuilding them leaves `index_fp` unchanged: a manifest made for a graph without a mask stays
      valid, with the same `index_fp`, after `metagraph transform --mask-dummy`. A manifest lists neither: one with
      an entry whose base name ends in `.edgemask` or `.bloom` (whichever graph it belongs to, loaded or not) is
      refused at start-up (server, single and multi-graph, and the CLIs), naming the entry and the rule
      (`kIndexDerivedDataRule`: "the dummy-edge mask (.edgemask) and the Bloom filter (.bloom) are derived data of
      the graph, not part of index_fp: they decide which counts are exact and how fast k-mers are looked up, never
      what an exact answer is, so adding, removing or rebuilding one leaves index_fp unchanged, and a manifest
      lists neither") and saying to write it again without the entry, which changes its `index_fp` once.
      `index_manifest.py` never writes them, refuses them with `--extra` and flags them with `--verify`.
      `traverse --index-inventory` (and `index_manifest.py --inventory`) lists them apart (`derived: [{path,
      role: graph_mask | graph_bloom, exists, loaded}]`, `derived_rule`), and the start-up log names the derived
      files loaded beside a checked manifest (`[Server] Derived data of <graph> loaded beside it, not part of
      index_fp: …`). **Limitation, stated:** no fingerprint covers derived data, so a stale or foreign derived
      file beside the graph is not detected by `index_fp`. The loader checks a mask's size and its W = `$` edges
      (`mask_invalid` on `/pattern`), and only a Bloom filter's k and mode: a Bloom filter of another graph with
      the same k and mode can hide k-mers (the pass-5 probe: `/resolve`'s `graph_runs` `[[0, 36]]` became `[]`)
      under an unchanged `index_fp`. Write derived data only with `build --mask-dummy`, `transform --mask-dummy`
      or `transform --initialize-bloom` on the graph it sits beside.
      A file of the inventory the manifest does not name is refused, naming it; a named one of another size too.
      *(Review of the pass-5 fixes.)* Loaded files are matched to entries by base name, so the entries' **base
      names must be distinct**: a manifest written for a directory of bundles sharing names (`A/graph.dbg`,
      `B/graph.dbg`, `A/annotation.seqs`, `B/annotation.seqs`) passed for each of them, and two single-index
      servers stated one `index_fp` and one `index_meta_fp` while answering with different accessions; such a
      manifest is refused, naming the shared base name (`index_manifest.py --verify` refuses it too, and writes
      bare base names). And it must **not list an optional file of the pair's inventory that the pair does not
      load** — the annotation's `.seqs`, when it is missing beside the listed spelling or with
      `--no-coord-mapping` —, so that `index_fp` identifies the loaded files (a server without the `.seqs` beside
      its symlinks, or with `--no-coord-mapping`, stated the `index_fp` of one that loads it). Files
      outside the inventory (`index_manifest.py --extra`: a `.coords`, a `.weights`, anchors beside another
      annotation type) stay allowed.
      Deliberately not in it, because the server does not open them: a column annotation's `.coords`
      (`merge_load` reads the columns only), the graph's `.weights`, an `<graph>.anchors` beside another
      annotation type (those are inside the annotation file), and the reverse header index (built in memory from
      the `.seqs`). *(Pass 5 review: before, the sidecars were not checked, and a manifest whose sidecars were
      another build's was accepted.)* — It also checks that the manifest lists no graph (`*dbg`) or
      annotation (`*.annodbg`) file other than those (a manifest describes one graph with one annotation: one
      written for a directory holding several annotations of one graph lent one fingerprint to each), and that
      an `index_fp` the manifest states is the computed one; it refuses to start otherwise. Sizes, not contents,
      are checked: a file replaced by another of the same size and name is not noticed (that is what hashing at
      build time and `index_manifest.py --verify` are for). `null` without a manifest: labels joined across
      retrievals are then unverifiable, never equal.
      `scripts/traversal/build_mini_refseq.sh` writes `<out>/annotation.relaxed.relabeled.manifest.json`
      (`MANIFEST_ONLY=1` writes it for an existing build); `scripts/traversal/index_manifest.py` writes one for any
      bundle (`-i/-a`), or one per pair of a server's graph list (`--server-csv`, below).
    - `index_meta_fp`: FNV-1a-64 (16 hex) over k, regime, alphabet, the annotation's row count and its ordered
      column names, and whether coordinates and a `CoordToHeader` exist — always present (computed per pair on
      first use in multi-graph mode), and only a **negative** check: a mismatch proves two indexes different,
      equality proves nothing (an index whose columns swap their memberships keeps it).
    `--index-name` and `--index-manifest` describe one index (`-i`/`-a`) and are refused with a graph list. The CLI
    (`metagraph traverse`) takes the same two flags.
- Single-graph mode: as `/search`; `graph` / `graph_path` in a request (or as probe parameters) are refused.
  **Multi-graph mode** (`server_query GRAPHS.csv`):
  - **The graph list** has one line per (graph, annotation) pair, `name,graph_path,annotation_path` and *(feature
    level 3)* optionally `,manifest_path` and `,manifest_path,index_ns` (an empty column means none; paths as the
    server opens them, relative to its working directory). A line of three columns is read as it always was. More
    than five columns, fewer than three, or an `index_ns` that is not `[A-Za-z0-9._-]+` refuses to start, naming
    the line (a fourth column used to be dropped silently, so a list that carried something else there — a
    release id, as an earlier version of this section described — now refuses to start rather than have it
    misread). Every listed manifest is checked against the files its pair loads exactly as `--index-manifest` is
    (sizes and the digest of its list, no re-hashing), before anything is loaded; a mismatch refuses to start. One
    index, one identity: lines naming the same index — the same files, however their paths are spelled: two lines
    are one index only when their pairs' **whole loader inventories** (its identity files: derived data does not
    split an index, so two spellings of one graph, one with a mask beside it, are one index with one identity)
    resolve to the same real paths in the same roles *(review of pass 5, finding 2: grouped by the graph's and annotation's real paths alone, two pairs of
    symlinks to one graph and annotation with different `.seqs` beside them shared one manifest check, and both
    stated A's `index_fp` while answering with different accessions)*, so each distinct bundle is validated
    against the manifest it names — must agree on what they state (an empty column states nothing and takes what another line of the
    pair states), and two different pairs never state one `index_fp` (a manifest whose files had the base names
    and sizes of two pairs, such as two annotations of one graph with swapped memberships, would let one index
    be taken for the other; copies of one index under two paths are refused too, as telling them apart would
    need hashing: list one path), else the server refuses to start. Every `/resolve` and
    `/traverse` response, and every graphlet's `H` record, states its pair's `index_ns`, `index_fp` and
    `index_meta_fp`. `scripts/traversal/index_manifest.py --server-csv GRAPHS.csv [--jobs N] [--digests FILE …]
    [--digests-only] [--match-base-names] [--out-dir DIR] [--write-csv OUT.csv] [--force]` writes those
    manifests: one per pair (to the manifest path its lines give — two different ones, or two different
    `index_ns`, are an error —, else `DIR/<name>.<annotation base>.manifest.json`, else next to the annotation),
    every distinct file hashed once by N parallel streams or taken from precomputed sha256 digests (`<sha256>
    <path>` lines, e.g. from storage checksums; a line names a file by its real path, relative paths taken from
    the digest file's own directory, and two different digests for one file are an error however spelled; by
    base name only with `--match-base-names`, only for lines whose paths name no file on this machine, only a
    base name one such line carries, never one line for two different files — bundles share names such as
    `graph.dbg`, and the server, checking sizes only, could not notice another bundle's digest —, each match
    listed; each entry records `"digest": "hashed" | "precomputed" | "precomputed-by-base-name"`), and
    `--write-csv` the list with every line's manifest column filled with its pair's manifest.
  - **Routing**: a request's `graph` (a name) is required; when the name lists one pair (the same pair listed
    twice is one) that pair is used; when its pairs share one graph with several annotations, the request is
    refused (400 naming them, with or without `graph_path`): a traversal reads one annotation, and `graph_path`
    cannot choose among them (list each annotation under a name of its own); when it lists several graphs,
    `graph_path` must name one of them (else 400 listing the graphs), and is refused like the above when that
    graph has several annotations under the name. There is **no union oracle**: a traversal reads one pair. The
    pair is not echoed (its identity is); a name with one pair ignores `graph_path`. The release id is
    server-wide (`--index-release ID`), in both modes.
  - **`graphs: [name]`** (the owner's decision of 2026-10-09): `/search`'s field, accepted by `/traverse` and
    `/resolve` as an alias of `graph` (one name: a traversal reads one graph; `graph_path` selects as above). With
    both fields, a list that is not one name, or on a single-graph server, 400. A request that selects with
    `graphs` is a fan-out's — the service sends `/traverse` to every selected chunk, without a presence check, as
    it sends `/search` — so a seed the chunk does not hold (or holds in part) is its result, `outcome.walks:
    "not_in_graph"` with `not_in_graph.kmers_present` (§6.1 step 2), not the whole request's 400; with `graph` the
    400 stays.
  - **`in_ram`** (the owner's decisions of 2026-10-09): `/traverse` and `/resolve` accept it exactly as `/search`:
    on a multi-graph server that runs on mmap, `in_ram: true` loads the selected pair into RAM for the request when
    its graph and annotation fit `--mem-cap-gb` (the load waits until that much of it is free — one pool with
    `/search`'s and `/pattern`'s loads, `/search`'s unbounded wait —, holds it while the request runs and frees it
    after; the reverse index of its headers is built in the load, as at start-up: seconds on a chunk with millions
    of records, inside `load_ms`); a pair above the cap, a server that loaded its
    graphs into RAM, and a single-graph server serve from the index they hold. The request's **budgets start after
    the load**: a `/traverse`'s time budgets, its memory and work budgets and its attempt's bound (`usage.bound`
    states the load as `load_ms`; the bound is `load + n_seeds × T + allowance`, the hard cap still counted from
    the header), a `/resolve`'s `bounds.time_budget_ms`. `timing.load_ms` (top level, with `in_ram`; a
    `/traverse` with `output.timing`) states the wait and the load, 0 when nothing was loaded. The HTTP server's
    content timeout counts both. A malformed request is refused before its load; the wait is abandoned when the
    client leaves or the server stops (nothing written), and a cancel of an attempt that waits takes effect when
    its work begins (its seeds not started). The loaded copy states the pair's identity: `index_fp` and
    `index_meta_fp` are per (graph, annotation) pair.
  - **`GET /capabilities`** of a multi-graph server (the owner's decision of 2026-10-09): `graph_summary`, so that
    a service learns every pair with one probe — `{columns_disjoint, shared_columns, pairs: [{graph, graph_path,
    annotation_path, index_ns, index_fp, k, graph_mode, available, unavailable_reason, mask, counting, traversal:
    {regime, num_labels, has_coordinates, has_coord_to_header, supports_trace}}]}`, one entry per (name, pair),
    the names in byte order, computed once at start-up (the schema is `SPEC-pattern-search.md` §24.4; about 470
    bytes per entry, which makes the route's document large on a long graph list). `columns_disjoint`:
    no column name is a column of two distinct pairs — the chunks of an index partition its samples —, which
    licenses summing a label's counts and occurrences over the pairs; `shared_columns` counts the names that are
    not (the start-up log names one). `null` on a single-graph server. Both GET routes state `in_ram`:
    `{"routes": ["pattern", "resolve", "search", "traverse"], "loads": <a multi-graph server on mmap with
    --mem-cap-gb above 0>, "mem_cap_gb": <--mem-cap-gb>, "budgets_start": "after_load"}`.
  - **`GET /traverse/capabilities?graph=<name>[&graph_path=<path>]`** (feature level 3; a 400 before): the probe of
    the pair those parameters select by the rules above (percent-decoded values), with `graph` and `graph_path`
    (the pair's graph) added; without `graph` a 400 (`Bad request: in multi-graph mode GET /traverse/capabilities
    needs ?graph=<name> (and graph_path=<path> when the name spans several graphs); GET /capabilities lists the
    graphs`); an unknown or repeated parameter is a 400. In single-graph mode `graph` and `graph_path` are a 400,
    any other parameter is ignored, as before.
- Caps: `--traverse-max-time-ms` (default 30 000; also the cap of `/resolve`'s `bounds.time_budget_ms`, which
  it lowers, never sets, and which is never above 899 000 ms on the server, §4.5), `--traverse-max-seeds` (64),
  `--traverse-max-seed-bp` (100 000), `--traverse-max-seed-labels` (10 000) and `--resolve-max-query-bp` (0 = unlimited). `0` means
  unlimited for each. A cap that lowers a request bound is **echoed** as `clamped` (§5); `max_seeds` and
  `max_seed_bp` are refusals (400), not clamps. **The budgets' maxima** *(feature level 4; R16, the owner's
  decision)*: `--traverse-max-memory-mb` (an integer in [0, 1 048 576]) and `--traverse-max-work-units` (an integer
  in [0, 2⁵³ − 1]), both 0 = off by default — a budget changes how a walk stops, so a deployment that did not
  choose one keeps the results it always gave. Set, a request's larger `bounds.max_memory_mb` /
  `bounds.max_work_units` is lowered to the maximum and an omitted one (no budget: unbounded) **set** to it, both
  echoed in `strategy.clamped` — `{"field": "bounds.max_memory_mb", "requested": 8192 | "unlimited", "effective":
  4096}` (`"unlimited"` for an omitted budget: how a knob spells no limit, which the graphlet's K and Q records
  hold, unlike null) — and in `strategy.bounds`; a smaller budget is kept, unclamped. Every request on such a
  server therefore runs under the budget the operator chose (and reads the annotation by the budget-aware path,
  §6.8). A seed such a budget stopped — its walk, or *(review of the efficiency pass)* its seed phase or depth-0
  state, failing the seed — states a `server_clamp` limitation for it (`limit` the maximum, `observed` what the
  request gave; "the request gave no budget and the server set its maximum, which this seed ran into; a request
  cannot raise it further", or that it was lowered); its `resource_stop.requested` is what the request gave, its
  `actions` offer no `raise_memory_budget` / `raise_work_budget`, and the knob's `walk_domain` entries (an arm's,
  or a failed seed's) carry `server_limit` (the maximum, as a `derivation`'s and `seed_labels`' do) and say "the
  knob is at the server's maximum (server_limit), which a request cannot raise" where they would say "raise the
  knob". A budget the request gave below the maximum keeps its lever. The probe states the maxima (`max_memory_mb`,
  `max_work_units`). The CLI applies none. The traverse caps are non-zero by DEFAULT because a request may
  name no labels at all: the server then chooses the permitted set, so the request's own bounds no longer describe
  the work it asks for, and an unconfigured deployment must not be the unlimited one. `--traverse-max-time-ms` is
  also the only bound on the derivation phase, so when a request has a derived seed a `bounds.time_budget_ms` of 0
  (which for the walk means "do not extend") is clamped UP to it rather than left unbounded — also echoed as
  `clamped`. The CLI (`metagraph traverse`), being operator-run, applies no caps. A request pinned to another
  `release` is a 400 naming both ids.
- Errors: a new `InvalidRequest` exception → 400 with the message; every other exception → 500 with context
  (small change to `process_request`). 503 while the index loads. Handlers run on the io threads: `-p ≥ 2` is
  documented as required for serving `/traverse` beside `/search`.
- **Attempts** (requests with `attempt_id`, §5; stage 4, backend half). Both capabilities routes state them
  under `attempts` (the same block): `fields` (`["attempt_id", "budget_id", "locus_id", "not_after_ms",
  "expect_server_instance"]`; the last from feature level 5), `cancel_fields` (`["attempt_id", "wait_ms",
  "not_after_ms"]`, feature level 5), `id_pattern`, `cancel`, `state`, `server_instance` (16 hex, random per
  process: a ledger tells a restarted backend, whose attempts and tombstones are gone, from an expired attempt),
  `retention_s`, `retention_count` (the finished attempts kept, and apart from them the most tombstones held at
  once), `tombstone_max_s` (feature level 5: the longest a tombstone is held to cover a `not_after_ms`, as
  applied — `max(--traverse-attempt-tombstone-max-s, retention_s)`), `suppression` (a reference to this
  section's `POST /traverse/cancel`: the tombstone), `release_rule` (to the release rule, below) and `instance`
  (to §5, `expect_server_instance`; the three from feature level 5), `allowance_ms`, `hard_cap_ms`
  (`content_timeout_s` × 1000 − 1000 = 899 000: the bound never exceeds it), `content_timeout_s` (900),
  `client_check_ms` (100), `clock_skew_allowance_ms` (`--traverse-clock-skew-ms`, default 2000: what a ledger adds
  to `not_after_ms`, §5; the server's own check adds nothing), `bound` and `not_after` (references to §6.8, the
  attempt's bound, and to §5, `not_after_ms`) and
  `delivery_reserve` (`compress_mbps`, `build_mbps`, `account_per_text_byte: {json, graphlet}` (30 and 50),
  `coordinate_account_per_text_byte` (feature level 6: 12, an integer, so that a ledger can reproduce the
  reserve, §6.8),
  `margin` (1.25), `stop_ms` (the configured time from the walk-until to the walk's end:
  `--traverse-chunk-target-ms` + 950), `calibration` (feature level 4: a reference to the measurements the
  starting estimates come from, §6.8, the delivery reserve: calibration), `measured_text_bytes`, `rate_window`, `measured_compress_mbps`, `measured_build_mbps`,
  `measured_account_per_text_byte: {summary, tree, full, graphlet}` and `measured_stop_ms` — the slowest rate,
  the smallest ratio and the longest stop time of the last `rate_window` this server measured, null until it has
  measured one — and its `rule`, a reference to §6.8, the delivery reserve). **Number types**: `allowance_ms`, `hard_cap_ms`,
  `clock_skew_allowance_ms`, `content_timeout_s`, `client_check_ms`, `retention_s`, `retention_count`,
  `measured_text_bytes`, `rate_window`, `stop_ms`, `measured_stop_ms` (or null), `tombstone_max_s`,
  `coordinate_account_per_text_byte` and, in
  `deadline_check`,
  `chunk_target_ms` and `observed_max_uninterruptible_ms` are integers (`allowance_ms` was written 10000.0 before
  feature level 3; `--traverse-attempt-allowance-ms`, `--traverse-clock-skew-ms` and
  `--traverse-chunk-target-ms` take integers in [0, 2⁵³ − 1], anything else is refused at start-up, so that a
  JSON client reads each as written and a ledger can add to it; *review of pass 5, finding 4*:
  `--traverse-attempt-retention-s` and `--traverse-attempt-tombstone-max-s` take integers in [0, 31 536 000] (a
  year) and `--traverse-attempt-retention` in [0, 10 000 000] — a negative, fractional, non-numeric or larger value
  refuses to start, naming the option and its range: `-1` had started, stated 2⁶⁴ − 1 seconds and expired every
  tombstone at once); `compress_mbps`, `build_mbps`,
  `account_per_text_byte`, `margin` and the measured rates and ratios (or null) are numbers (a rate may be
  fractional), as are the probe's `max_time_ms` and the time budgets.
  - `POST /traverse` with `attempt_id`: 200 with `usage` (complete, partial, failed seeds, cancelled); 400 with
    `{error, usage}` for a request error after the attempt was registered (a bad strategy, an unknown label or
    graph, too many seeds, an unreachable extra label, …); 400 with `{error}` alone for a malformed id, a malformed
    `not_after_ms`, a malformed `expect_server_instance` (§5: the ids' pattern), or `budget_id`, `locus_id` or
    `expect_server_instance` without `attempt_id` *(the last two: review of 2026-10-06, X4)*; **409** with
    `{error, attempt}` (the other
    attempt's state, no `usage`) when the id is running, retained, held past its retention (below), or tombstoned
    by a cancel (a tombstone's `attempt` carries its suppression, judged against this request's own `not_after_ms`
    alone — absent when it has none: then `covers_admission` is false, `no_not_after_ms` — and the refusal extends
    the tombstone, below); **409** `{error, state: "instance_mismatch", …}` when
    `expect_server_instance` names another process (§5; answered first); **409** with
    `{error, state: "expired", not_after_ms, server_time_ms, attempt_id, budget_id?, locus_id?, server_instance}`
    (no `usage`, nothing registered) when `not_after_ms` has passed (§5; the same without the ids for a request
    without `attempt_id`; the duplicate's 409 is answered first); 503 with `{error, usage}` (reason `deadline`)
    when the bound was reached while the response was built; 500 with `{error, usage}` for a non-standard
    exception; 503 while the index loads (no `usage`: the request was not read). A client that is gone gets no
    response (with or without `attempt_id`).
  - `POST /traverse/cancel` `{"attempt_id": "…", "wait_ms": 0, "not_after_ms": …}` (`wait_ms` 0 … 10 000, default
    0: how long to wait for the attempt to finish; `not_after_ms` optional, feature level 5, validated as on
    `/traverse` — the instant of the request being cancelled): **200** `{attempt_id, server_instance, cancelled:
    true, state: "stopping" | "finished", attempt}` when the attempt was asked to stop (now or before: idempotent;
    `state` `finished` when it finished within `wait_ms`); **404** `{error, attempt_id, server_instance, cancelled:
    false, state: "finished", attempt}` when it has finished; **404** `{error, attempt_id, server_instance,
    cancelled: false, state: "unknown", tombstone: true, suppressed_until_ms, not_after_ms?, covers_admission,
    covers_admission_reason?}` for an id this process does not know — the id is then **tombstoned**: a request
    arriving later with it is refused (409) and never runs here (a cancel can overtake its request: still queued,
    on the wire, its body half uploaded). *Review of pass 5, finding 1*: held for `retention_s` alone, a tombstone
    expired while the request it suppressed could still be admitted — with `retention_s` 1, a request with
    `not_after_ms` = now + 60 s was half uploaded and cancelled (404, `tombstone: true`), and its upload completed
    at 1.18 s ran. A tombstone is now held **at least `retention_s`** from the cancel and, when the cancel names the
    request's `not_after_ms`, **until this server's clock reads `not_after_ms + clock_skew_allowance_ms`**, at most
    `tombstone_max_s` from the cancel. It is held while **either** clock says it is live: the steady clock for the
    duration computed when it was set or extended, plus one millisecond (a forward step of the wall clock does not
    shorten it), and the wall clock while it reads **at most** `suppressed_until_ms`, that millisecond included (a
    backward step lengthens it). *(Review of the pass-5 fixes, finding 2: the tombstone ended when the clock read
    `suppressed_until_ms` exactly, while the strict check still admits a request in the millisecond its clock
    reads `not_after_ms`; with `--traverse-clock-skew-ms 0` a covered copy whose handler started in that
    millisecond ran, 27 of 40 tries.)* A repeated cancel, and a
    refused copy of the request (its 409 above, with its own `not_after_ms`: a proxy's replay carries the same),
    **never shorten** it: each extends it to at least `retention_s` from then and to its own `not_after_ms` +
    `clock_skew_allowance_ms` within the cap, so every later copy is refused while it could still be admitted.
    Every tombstone answer — the cancel's 404, a repeated cancel's 404, `GET /traverse/attempt`'s 404 and the
    refused request's 409 (in its `attempt`) — states `suppressed_until_ms` (the hold's last millisecond on the wall
    clock, inclusive, Unix epoch ms on this server's clock), `not_after_ms` (the one it was judged against: for the
    refused request's 409 its own alone, absent when it has none — *review of the pass-5 fixes, finding 3*: it fell
    back to a cancel's, and stated `covers_admission: true` for a copy without `not_after_ms` that ran once re-sent
    after the hold; for a cancel its own, else, as for GET, the largest named for the id; absent when none was),
    `covers_admission` (**true iff a `not_after_ms` was judged against and `suppressed_until_ms ≥ not_after_ms +
    clock_skew_allowance_ms`**: no copy of that request can start on this `server_instance` — while the tombstone
    is live it is refused, and once it is gone this server's clock has read later than `suppressed_until_ms`, so
    the copy's `not_after_ms` has passed on the clock the strict check reads) and, when false,
    `covers_admission_reason` (`no_not_after_ms` | `beyond_tombstone_max`). Tombstones are kept apart from the
    finished attempts and expire by their hold only, so later finishes never evict one early (review of the
    stage-4 backend, F6); a cancel of an unknown id is tombstoned only while fewer than `retention_count`
    tombstones and held finished attempts (below) are held, and refused beyond that, **429** `{error, attempt_id,
    server_instance, cancelled: false, state: "unknown", tombstone: false, reason: "tombstones_full"}`: the id
    was not tombstoned and nothing is promised (retry the cancel later). With `retention_s` 0 (finding 4) the
    server keeps **no tombstones**: such a cancel is refused, **429** with `reason: "no_suppression"`, and a request
    with the id arriving later runs (the capabilities' `suppression` says so). 400 for a bad body. Not refused
    while the index loads. What a clock step can do: a backward step of this server's wall clock by more than
    `clock_skew_allowance_ms` after a tombstone expired can make a copy's `not_after_ms` lie in the future again
    (exactly: the hold ends once the clock reads later than `not_after_ms` + skew, and the strict check admits a
    copy only while it reads at most `not_after_ms`; the assumption the ledger's skew allowance makes anyway).
  - **The release rule** (normative; the capabilities' `release_rule`, feature level 5; its grounds completed
    at feature level 6 by the review of 2026-10-06 — X2, and the search service's release parity, LRG-R1/R2 —
    as text only: what the server does is unchanged): a ledger may release an
    attempt's capacity **before it finished only** on a 404 with `tombstone: true` **and** `covers_admission:
    true` from the same `server_instance`, for an attempt it sent with **exactly that `not_after_ms`** and with
    `expect_server_instance` equal to that `server_instance` (a copy reaching a restarted process is then refused,
    `instance_mismatch`). An attempt sent without `not_after_ms`, or without `expect_server_instance`, is never
    released early on a tombstone. Otherwise it releases only on a finished state (`GET /traverse/attempt`, a
    cancel's 404 or 200 with `state: "finished"`, the 409 refusing a copy of the request whose `attempt` — the
    id's state as `GET` answers it, with its `attempt_id` and `server_instance` — says `state: "finished"`, or the
    response), on an `expired` 409 (below), or once its own clock passes `not_after_ms +
    clock_skew_allowance_ms + bound_ms` (the clock release, below). **A finished state** promises that the request does not run again on this
    process (and, pinned, anywhere; below): a finished attempt's id stays refused (409, its state) while it is
    retained (`retention_s`, `retention_count`) and, for an attempt sent with `not_after_ms`, until this server's
    clock reads later than `not_after_ms + clock_skew_allowance_ms`, at most `tombstone_max_s` after it finished
    **or after its latest refused copy** — past its retention if need be, **held** among the tombstones (on both
    clocks, as a tombstone; a copy refused with its own `not_after_ms`, an identical replay included, extends the
    hold to that + `clock_skew_allowance_ms`, at most `tombstone_max_s` from the copy's arrival: *review of
    2026-10-06, C30*, "at most `tombstone_max_s` after it finished" left the refused copies out). *(Review of the pass-5 fixes, finding 1: with `retention_s` 1, or
    `retention_count` 1 and one later finish, a replay of a finished request — same `attempt_id`, `not_after_ms`
    and `expect_server_instance` — arrived after the ledger had released on the finished state and ran again.)*
    **The hold is the process's, in memory:** a restarted process (a new `server_instance`) holds none, and a copy
    of the request reaching it runs there unless it was sent with `expect_server_instance` (then 409
    `instance_mismatch`). So a finished state is **replay-safe** — no copy of the request runs after it, on this
    process or a restarted one — **only for an attempt sent with `expect_server_instance`** equal to that
    `server_instance` **and with `not_after_ms`** whose `not_after_ms + clock_skew_allowance_ms` lies within
    `tombstone_max_s` of the finish: a ledger that releases on a finished state before its clock passes
    `not_after_ms + clock_skew_allowance_ms + bound_ms` pins the instance, as for a tombstone. *(Review of levels
    4–5, finding 6: a request finished with `not_after_ms` 60 s ahead, without `expect_server_instance`, replayed
    after a restart of the same endpoint before then ran again — 200, 170 work units — while the text said a
    finished request within its hold was never run again; the pinned copy was refused, `instance_mismatch`.)* A
    finished state of an attempt sent **without** `expect_server_instance` assumes that no copy of the request
    reaches a restarted process (or another server at the same address) while it could still be admitted there;
    one sent **without** `not_after_ms`, or with one beyond `tombstone_max_s` of its finish, assumes that no copy of
    the request arrives after the attempt left retention; with `retention_s` 0 nothing is kept, and a finished
    state assumes that no copy arrives after it (a pinned copy is still refused by a restarted process). Held
    attempts are not bounded by `retention_count` (they are never dropped early): their number is at most the
    finishes and refused copies within `tombstone_max_s`, so `--traverse-attempt-tombstone-max-s` and the rate of
    refused copies bound their memory (copies that keep arriving keep an id held, and held attempts count toward
    `retention_count` for the tombstones of unknown ids: 429 `tombstones_full`).
    **An `expired` 409** (`state: "expired"`) for an attempt sent with exactly that `not_after_ms` — from the
    `server_instance` it named in `expect_server_instance`, when it named one — **releases** it: the server checks
    the id's registration (running, retained, held, tombstoned) before the expiry, under one lock (§5), so no copy
    of the request was running or registered on that `server_instance` when it was judged, and none can start
    there later. It **settles nothing**: an earlier copy may have run, finished and left retention, and its usage
    is then unknown. It assumes that this server's clock does not step back below `not_after_ms` (the 409 states
    `server_time_ms`; `server_time_ms − not_after_ms` is the step it survives, and once that exceeds
    `clock_skew_allowance_ms` the step is the one assumed below) and, for an attempt sent without
    `expect_server_instance`, that no copy of the request reaches another server at the same address while it
    could still be admitted there *(LRG-R2: the search service releases on it; the rule's list omitted it)*.
    **The clock release** takes the attempt as stopped once the ledger's clock passes `not_after_ms +
    clock_skew_allowance_ms + bound_ms`: it is then past its bound **apart from what it runs past it until its
    next delivery check** — the rest of the piece the bound fell into and, when its walk had not stopped, the
    walk up to its next poll that reads the clock, the stopped seed's finalisation and the building up to that
    check (§6.8, the attempt's bound) —, **a run whose length has no stated bound**
    (`deadline_check.max_uninterruptible_ms` is null, and `observed_max_uninterruptible_ms`, which leaves out the
    mapping of a seed's k-mers, the seed phase's own processing and a seed's finalisation, is an observation, not
    a margin). A `GET /traverse/attempt` that still answers `running` or `stopping` after that instant shows that
    run still going: a ledger that must not let two attempts overlap keeps polling while it answers so, or adds
    a margin of its own. *(Review of 2026-10-06, X2: one 6.94 Mbp seed, allowance 0, skew 0, bound 1,500 ms —
    GET answered `running` up to 2,720 ms past the ledger's release instant, and the POST ended in a 503 at
    4,542 ms, while `observed_max_uninterruptible_ms` said 0–1 ms: the k-mer mapping is not observed. Its fixes
    first said "one uninterruptible step"; that run is several pieces, the review of the fixes found. The
    library's `release_verdict` holds the clock release while the last answer said `running` or `stopping`, by
    default.)*
    A 429, a cancel's 200 with `state: "stopping"`, a 409 whose `attempt` says `running` or `stopping` or is a
    tombstone, and a 409 `instance_mismatch` release nothing.
    Assumed, not checked: the ledger's clock is within `clock_skew_allowance_ms` of this server's, this server's
    clock does not step back by more than that, nothing between the ledger and the server rewrites a request's
    `attempt_id`, `not_after_ms` or `expect_server_instance`, no copy of the request reaches **another
    server** that serves the same ledger (a server knows only its own tombstones and holds;
    `expect_server_instance` refuses such a copy only where it names that server's instance), and, for the clock
    release, that the attempt's run past its bound has ended by then. The capabilities' `release_rule`
    refers here.
  - `GET /traverse/attempt/{attempt_id}`: **200** with the attempt's state, **404** `{error, attempt_id,
    state: "unknown", server_instance}` (and for a tombstoned id `tombstone: true` with its suppression, as the
    cancel states it; reading it does not extend it) once it is no longer
    retained, 400 for a malformed id (the id is matched as sent, not percent-decoded):

    ```json
    {"attempt_id": "a-17", "budget_id": "b-3", "locus_id": "l-9", "server_instance": "9f3c0d1e2a4b5c6d",
     "state": "finished", "reason": "cancelled", "stop_requested_by": "cancel",
     "received_at": "2026-10-03T19:17:48.808Z", "cancel_requested_at": "2026-10-03T19:17:49.812Z",
     "stop_requested_at": "2026-10-03T19:17:49.812Z", "stopped_at": "2026-10-03T19:17:49.812Z",
     "finished_at": "2026-10-03T19:17:49.943Z", "deadline_at": "2026-10-03T19:19:28.808Z",
     "bound_ms": 100000, "elapsed_ms": 1136,
     "seeds": {"requested": 3, "started": 1, "finished": 1, "abandoned": 0},
     "response": {"written": true, "status": 200, "bytes": 13721360}, "usage": {"…": "as in §7.3, without per_seed"}}
    ```

    The state carries `not_after_ms` at its top level when the request gave one (pass 5). `state`: `running`,
    `stopping` (a stop was requested — a cancel, the bound, a gone client — and the handler has not returned) or
    `finished` (the handler returned: the walk's memory is freed and the response, if any,
    handed to the transport, which holds its bytes until they are sent or the content timeout). `reason` (once
    finished): `completed | cancelled | deadline | client_gone | error`; a walk that completed before it saw a
    cancel stays `completed`. `stop_requested_by`: `cancel | deadline | client_gone | null` (the first reason
    wins). `stopped_at`: when the walk stopped (the checkpoint that saw the stop, or the end of the last seed);
    `response.bytes` is the body's size as written (compressed when the client asked for it), `written: false`
    when the client was gone. `seeds` as in `usage` (§7.3): a seed whose walk the client's departure cut is
    `abandoned`, not `finished`. `usage` is null until the attempt finished, and keeps the totals only (no
    `per_seed`).
  - Retention: finished attempts are kept `--traverse-attempt-retention-s` seconds (3600) and at most
    `--traverse-attempt-retention` of them (10 000, oldest dropped first; a running attempt is never dropped);
    GET answers 404 after that and the id may be used again — unless the attempt was sent with `not_after_ms`:
    its id is then **held** (GET 200, a request with it 409) until `not_after_ms + clock_skew_allowance_ms`
    within `tombstone_max_s` of its finish or of its latest refused copy (the release rule, above). Tombstones are held **at least**
    `retention_s` and, for a `not_after_ms` a cancel or a refused copy names, until it + the skew allowance
    within `tombstone_max_s` (above); a cancel of an unknown id is tombstoned only while fewer than
    `retention_count` tombstones and held attempts are held. Every attempt is logged when it is
    registered, cancelled and finished (`[Server] Attempt <id> (request <n>) finished (<reason>): walk stopped at
    <instant> (<d> ms after the stop was requested by <who>), <s> of <n> seed(s) walked[ (<a> abandoned in its
    walk)], bound <b> ms; response
    <status>, <bytes> bytes, <ms> ms`, or `no response written (the client is gone)`), a request without
    `attempt_id` whose client is gone too (`[Server] Request <n>: client gone, walk stopped at …; no response
    written`).
  - Threads: the cancel and state routes run on the same io threads as the traversals; with every thread busy
    they wait. `-p` must exceed the number of concurrent traversals (by one at least) for a cancel or a state
    query to be answered while they run; the attempt's clock starts when its header is read, so a request queued
    behind busy threads is not yet on it.
- **Deployment** (`DESIGN-traverse-graphlet.md` §17.5):
  - `-p` must exceed the number of concurrent traversals by one at least, so that the cancel, state and
    capabilities routes find an io thread while traversals run.
  - A reverse proxy in front of the server must forward `POST /traverse/cancel`, `GET /traverse/attempt/*`,
    `GET /capabilities` and `GET /traverse/capabilities` with their query strings (`?graph=…&graph_path=…`,
    percent-encoded), must close its upstream connection when its client goes (the server stops a traversal whose
    client is gone; a proxy that keeps the upstream connection open hides that), must not buffer a request's
    answer past the client's own timeout, and must keep its read timeout at or above `hard_cap_ms` plus the
    transfer of the largest response.
  - The server's clock must be synchronised (NTP) to within `attempts.clock_skew_allowance_ms` of the ledger's,
    for `not_after_ms` to mean the same instant to both; the stated allowance is configuration
    (`--traverse-clock-skew-ms`).
  - A graph list that gives manifests (column 4) refuses to start when one does not describe its pair's files:
    regenerate them with `index_manifest.py --server-csv` after an index changes. Manifests written before the
    review of pass 5 may lack files the server now checks (`<base>.seqs` of a coordinate annotation, which the
    tool used to look for under other names), or list a second annotation (`--extra`): such a list refuses to
    start, naming the file; regenerate the manifests (`index_manifest.py`, whose default bundle is now the loader
    inventory: it no longer adds a `.coords` or a `.weights` the server does not open — `--extra` adds them,
    changing the fingerprint). A manifest that lists a mask or a Bloom filter (written between the review of
    pass 5 and the owner's decision #17 of 2026-10-08) refuses to start, naming the rule: write it again without
    them (its `index_fp` changes once; after that, adding a mask leaves it unchanged).
  - **The sequence header index** (header name → `(column, seq_id)`) is built while the index loads, before the
    server answers traversal requests, and logged with its size and time (`[Server] Sequence header index built:
    N headers in S s`; 4.6 s for 33 M synthetic headers on the M5 Max): built by the first request naming a header
    instead, it cost that request seconds inside its time budget and outside every read timer (R10: staging,
    refseq33m, a first header-named walk ran 9,947 ms on a 1,000 ms budget, 689 ms of it reading). It holds a
    hash entry per header (on refseq33m, 33 M headers, on the order of 1-2 GB), which a server that serves header
    labels would have built anyway. The CLI builds it once before its requests.

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
| T30c | per-seed derivation failure | Three seeds, the middle one carried in full by no label: HTTP 200 / exit 0 with `{"seed": {...}, "outcome": {"walks": "failed", …}, "error": …, "limitations": [derivation]}` in its place (cause `no_carrier`, knob `seeds[].sequence`) and the other two traversed with their own `outcome` (`Walker.FailedDerivationIsAStatedOutcome`, also `over_seed_label_cap`; a time budget the server lowered so far that it runs out after the derivation's first k-mer is, from feature level 6, not a failure but D3's walked result — `partial/complete/qualified`, `observed` 1 k-mer, `server_limit` — and a budget spent before the first k-mer still fails the seed, `WalkerDeadlineChunks.DerivationWindowIsPaced`, T57; *review of 2026-10-06, X1*: this row named `time_budget` among the failures it asserts); the same seed with an EXPLICIT label list still fails the whole request; `max_seed_labels` out of `[1, 100000]` is a parse error; a seed over `--traverse-max-seed-bp` and more than `--traverse-max-seeds` seeds are 400s; `max_seed_labels` above `--traverse-max-seed-labels` is clamped and reported | Hit-labelled seed |
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
| T42 | the request budgets' bounds (stage 2, the re-review) | `Graphlet.DeliveryCostsBoundTheOutput`: per detail, the model is at least the demonstrated peak of the serialisers (the JSON tree from this platform's jsoncpp types with `writeString`'s three texts, the compressor's, the graphlet writer's), for every budget case and for names of control characters, `%`, quotes, tabs and DEL and 2-, 3- and 4-byte UTF-8 (seed, extra, recorded and dropped labels, the seed_id); the walker charges at least the model; the output keeps to the widths the model assumes; with a fixed four bytes per name byte the case fails. `GraphletCodec.JsonEscapedSizeIsTheWriters`. `Stage2ReviewRecheck.*`: a work stop within one row of the budget (annotate root row, a constrain level of 3,000-label rows; the message states the most charged between two comparisons), an interrupted fetch stating what it held, failed derivations (no carrier, the time budget) stating `memory_bound_soft`, a walk complete to its radius under `time_budget_ms: 0` stating no stop, an injected refusal saying so, escaped names within the budget (180,000 control characters, 200,000 `%`). `LabelOracleHeaderIndex.IsBoundToItsCoordToHeader`: the reviewer's `[A,X,Y,Z]` → `[A,Y,X,Z]` at one address resolves by the live headers, and the index is built once per `CoordToHeader`. Integration `test_stage2_*`: the reviewer's CLI cases (25,000 headers under a work budget of 1; 100,000 headers under 1 MiB; 180,000 control characters and `%` headers under 1–2 MiB) | — |
| T43 | the review of the stage-2 fixes (round 3) | `Stage2ReviewRound3.*`: a failed seed's echo — an `ambiguous_header` of 180,000 control characters echoed as a bounded prefix with its length and place under a budget (whole without one), a `seed_id` of 100,000 emoji whose failed result's excess is stated as `memory_bound_soft`; the label validation observing a row of 200,000 coordinates before the charge that fails the seed; every fetched row charged (a level's second head refused after the fetch: work ≥ rows_requested × (8 + n) at width n); work stops within their stated largest charge (both arms' roots, a row's trace coordinates, a budget sweep over 1,000-label constrain rows and label-state scans); a failed depth-0 admission stating the dictionary it held (six 200 KB names, both modes). Integration `test_stage2_*` (the same through the CLI, including a budget sweep on an index of uniform row width). Python `test_traverse_review3.TestRound3*`: compare() equal at a radius equal to a merge position (constrain 36, annotate 36 and 62) with true differences kept; malformed continuation overrides, arms and resolve labels as bounded errors; a 30-label continuation's receipt within the default ceiling, the per-label loss list cut first; a label alive at the leaf but not seeded stated in `notes` and `left_out` | — |
| T44 | the recheck of the stage-2 fixes and the review of stage 3 | `Graphlet.DeliveryCostsBoundTheOutput` with wide floats (switches of 1e-300 on alternating labels, sums of 0.1 under a loss budget of 1e300, a time budget of 1e-300: every R, E and float field within the priced width). `GraphletStage2Recheck.*`: every row a seed-phase read returned charged before its comparison (validation: exactly the rows read; derivation: the window, 400,019-unit analogue), a failed work seed stating its largest charge and the threshold wording, the trace validation's coordinate copies observed (four copies under 1 MiB), extra labels reachable by a chain accepted and unreachable ones refused by name, `switch_reach` against a fixpoint over every pair (and 20,000 labels under a dense default), `mgt_float_width` against encoded costs, sums and times. `GraphletStage3Review.RadiusOnlyLevelIsNotStopped` (radius 0 and 1, 256 budgets from the first admitted). `LabelOracleBudgeted{Query,Recorder}.ZeroCacheWarmReturns`, `.AlternatingFetchPathsKeepCosts`. Python `test_traverse_stage2_recheck`: the spooled receipt's boundary, `left_out` per walk, the branch allowance (reduced, reset, unlimited, the caller's, the tool's receipt; CLI: the continuation stops at 46 like one walk), extra labels by chains (library and CLI) | — |
| T46 | the reviews of the stage-3 fixes and of the stage-4 backend | `GraphletStage3Fixes.MemoryStopsStateTheBudgetThatHoldsThem`: `memory_budget_holding` against 20,000 random needs and budgets (the smallest whole MiB whose own allotments fit beside the need), the roots case at 1 … 12 MiB (every depth-0 failure states the same knob value in its message and `walk_domain`; raised to it the state is held, one MiB less is not), and head stops of every budget case (the cap trigger's demand is the knob value; raised to it the refused head is admitted). `.TooWideWindowIsTheSeedsWork` (column and row-diff, with and without a work budget). `GraphletAttempt.UsageStatesWhatFailedResultsHold` (a depth-0 failure and a never started seed with 2 MiB `seed_id`s under 1 MiB: each seed's `soft_excess_bytes` is what its `memory_bound_soft` states, `final_bytes` the priced failed result, `refused_bytes`, the admitted peak within the budget, the request's peak and `held_bound_bytes`; no bound without a memory budget), `.AbandonedWalksAreNotFinished`. `GraphletAttemptRegistry.TombstonesOutliveLaterFinishes` (ten finishes within `retention_s` leave the tombstone; beyond `retention_count` tombstones a 429; expiry by age). `GraphletServer.PeerClosedSeesACloseBehindWaitingBytes` (a trailing CRLF and a pipelined request, then close, half-close or reset). Integration `test_traverse_attempt_errors_state_their_usage` (CLI), `test_api_attempts` (a failed seed's usage under a memory budget, an error's usage in the library, the MCP tool's capabilities at its default ceiling from the real server), `test_a_close_behind_waiting_bytes_stops_the_walk`, `test_the_registry_over_http`, the library's `AttemptAtBound` at the bound. Python `test_traverse_stage4_review`: the note's chain through a kept extra label, the two 503s told apart (`AttemptAtBound`, `ServerInitializing`), an error's `usage`, `traverse_capabilities` at its default ceiling (2.4 KB and 12 KB) while an explicit ceiling holds, a 429 cancel returned | — |
| T45 | stage 4, backend half: attempts (cancel, client gone, usage, the bound) | `GraphletAttempt.StopAtEveryPollLeavesAConsistentPrefix`: a cancel and the attempt's bound injected at sampled polls of every budget case leave every walk up to `complete_to_bp` the unstopped walk's, a later stop never walks less, the stop is stated (`Q attempt cancelled\|attempt_deadline traversal`, a `walk_domain` naming `attempt_id`, no "raise the knob") and serialises both ways; in the seed phase the seed fails. `.NoStopIsByteIdentical` (with and without budgets), `.GoneClientAbandonsTheWalk`, `.SeedPhaseStopFailsTheSeed` (validation and derivation, through the request), `.UnstartedSeedsAreFailedWithTheStop`, `.UsageStatesWhatTheSeedsConsumed` (per seed the walk's account; with `attempt_id` the response is otherwise byte for byte the one without), `.StatementsFitTheirWidths`, `.BoundIsEnforcedOnTheAttemptsClock`, `.IdsAreValidated`, `.GoneClientAbandonsTheAttempt`. `GraphletAttemptRegistry.*`: an id runs once (running, retained, expired), cancel idempotent and acknowledged, a cancel of an unknown id tombstones it, retention by count and age, `wait_ms`, eight threads starting, cancelling, reading and finishing. `GraphletServer.PeerClosedTellsAGoneClient` (closed, half-closed, reset, a pipelined request waiting), `.CheckedWriterIsByteIdentical`. Integration `TestTraverseAPI.test_api_attempts` (usage on success, partial, a failed seed and a 400; 409 and malformed ids without usage; 404s; the library) and `TestTraverseAttempts` (a slow index: a cancel mid-walk with the client connected, a client that goes away with and without `attempt_id`, `wait_ms`, the bound enforced, retention, sixteen concurrent attempts with mixed cancels) | — |
| T47 | pass 5: `not_after_ms`, per-graph identity, capabilities | `GraphletAttempt.NotAfterMsIsParsedStrictly`, `GraphletAttemptRegistry.NotAfterMsRefusesAtStart` (strict at the instant, the exact 409 body, nothing registered, a duplicate answered first, usage keys unchanged without it), `GraphletAttemptRegistry.CapabilitiesStateIntegers`; `GraphletServer.GraphListLinesAreParsedStrictly` and `.GraphListIdentitiesAgreePerPair`; `test_api_not_after_ms` (server and CLI, with and without `attempt_id`, the library's `AttemptExpired`), `test_api_server_capabilities`; `TestTraverseMultiGraph` (identity per pair in responses and `H`, the probe per pair, the several-annotations refusal, start-up refusals for a wrong manifest, conflicting identities, six columns, a bad name, a missing manifest; a three-column list keeps nulls; review of pass 5: a pair listed twice needs no `graph_path`, a name over one graph with several annotations is refused as such, and a directory's manifest, another build's sidecar and one `index_fp` for two pairs refuse to start, while two spellings of one pair are one index); `Graphlet.IndexIdentity` (the sidecars checked, another graph or annotation in a manifest refused); `test_traverse_index_manifest.py` (batch = single mode, each file hashed once, parallel streams deterministic, precomputed digests, refusals; review of pass 5: digests by real path from the digest file's directory, never another bundle's by base name, a pair's manifest path from any of its lines, the sidecars the server loads) | §5, §10.3 |
| T48 | pass 5: chunked deadlines | `LabelOraclePacing.*` (paced reads answer as one read on the direct, rows, tuples and budget-aware paths, recorder included, with evictions, in one-row chunks; an interrupted read changes nothing and states its decoded work); `WalkerDeadlineChunks.*` (a 100 ms budget held within a chunk where the whole read overran by 150+ ms, the same `complete_to_bp`; no stop: the same result; the attempt reaches a validation and a level's read within a chunk when its walk-until is near, and after the read when it is far; a derivation's window; review of pass 5: `FarDeadlineReadsAreOnePiece` — rows that share a per-call cost, as row-diff rows share their paths, read whole far from the deadline: the same calls, result and time — and `SplitReadsStartSmall` — a read whose rows are ten times slower than earlier reads' starts with a small chunk); `LabelOraclePacing.PacerSizesChunksByTime` (one piece far from the deadline, the first chunk at most 8 rows, growth and rest rules); `MiniRefSeq.PathCacheKeepsTheResponse` (reads in one piece and in one-row chunks, budgets included); `TestTraverseWideIndex` (budgets of 50–200 ms kept within a chunk + 150 ms on the fan-out index; review of pass 5: walks no deadline stops keep their bytes in the unchunked time on `row_diff_brwt`, budgets of 50 and 200 ms kept with and without a memory budget — and the integer flags refused when negative, fractional, above 2⁵³ − 1 or malformed) | §6.8 |
| T49 | pass 5: delivery | `GraphletServer.AssembledResponseIsByteIdentical` (per-seed texts assembled = the whole tree's text, every detail, a failed seed, usage); `GraphletAttempt.DeliveryReserveMovesTheWalkUntil` (with the 1.25 margin, the configured and the measured stop time, and `usage.bound.walk_until_ms` the walk-until in force when it stopped the walk); `test_api_bodies_are_the_cli_output_at_every_encoding`; `test_wide_index_delivery_reserve_stops_the_walk`; the byte-identity harness against the previous build (1,380 real requests at chunks of 50 ms, 588 of them at 1 ms, 588 from a three-column graph list, 804 SRA requests) | §6.8, §10.3 |
| T50 | pass 5: `standalone_text` and the reduced `J` | `test_traverse_standalone.py`: byte-equal to `dump(from_response(…), envelope=True)` on every fixture result, `usage` reduced to the totals and this seed's `per_seed`, `save()` and the store write the same bytes, reserved names refused | §7.5.2 |
| T51 | the efficiency pass (feature level 4) | `RowDiffPathCache.*` (the cached default decode returns the default decode's rows call after call, under bounds that keep everything, evict often and keep nothing, and a shared bound that shrinks; the budget-aware decode with the cache gives every row the costs of its whole path and at most the held bytes, a row read alone at most its peak without the cache and within its demand, and refuses below its peak; the generations; a dropped generation frees its table, `DroppedGenerationsReleaseTheirTables`, and the counts follow the tables at a rotation and under a shrinking shared bound; the lookahead's reads admitted as without the cache, `LookaheadAdmittedAsWithoutTheCache`), `LabelOraclePathCache.SameAnswersAndCounters` (query and recorder, budgeted and not, every counter, the shared bound), `MiniRefSeq.PathCacheKeepsTheResponse` (96 responses byte-equal with and without the cache: constrain/annotate, memory budgets at their stops, work budgets, `batch_kmers` 1 and 64, an evicting cache, one-row chunks), `MiniRefSeq.PathCacheKeepsRootRefusals` (an annotate seed whose right root's read the seed_id's length moves across its refusal under 1 MiB: every response byte-equal with and without the cache, the refusals of a read alone among them), `MiniRefSeq.ServerBudgetMaxima` (also a seed failed in its seed phase by the work maximum and at an annotate root by the memory maximum: server_clamp, `requested`, no raise action, `server_limit`; a seed failed by the request's own budget keeps its lever), `CoordToHeader.SequenceRangeAgreesWithMapSingleCoord`, `LabelOracleCoordRuns.SameAsMapSingleCoord` (and the measurements in `benchmarks/traversal/`); `test_row_diff_path_cache_keeps_the_bytes`, `test_server_budget_maxima`, `test_api_server_capabilities`; the byte-identity harness against the previous build (unbudgeted: 588 mini_refseq — also from a three-column graph list and with chunks of 1 ms —, 792 UHGG and 804 SRA `/traverse` requests, 688 `/resolve` requests; budgeted: the 588 mini_refseq requests under 4 MiB and under 300,000 work units, the 792 UHGG ones under 16 MiB; differences only in attempts' `walk_until_ms` and in walks the previous build's cold first touch cut by their time budget), and with the cache on and off (201 SRA and 198 UHGG requests at the default `batch_kmers`, at 1, the UHGG ones under 256 and 8 MiB and 10⁶ work units, 147 mini_refseq requests without and under 4 MiB) | §4.2, §6.8, §8.2–8.4, §10.3 |
| T52 | the review of pass 5 and R10 (feature level 5) | `GraphletAttemptRegistry.CancelNotAfterMsHoldsTheTombstoneThroughAdmission` (the reviewer's half-upload timeline with retention 1 s on injected clocks: the cancel naming `not_after_ms` holds the tombstone to `not_after_ms + skew`, the repeat at 0.763 s states the same expiry, the copy at 1.177 s is refused, at the expiry the strict check refuses it; without `not_after_ms` `covers_admission` is false and the copy runs), `.TombstonesAreNeverShortenedAndCapped`, `.ARefusedCopyExtendsTheTombstone`, `.TombstonesSurviveWallClockSteps` (forward and backward steps), `.RetentionZeroKeepsNoTombstones`, `.ExpectServerInstanceRefusesAnotherProcess`, `.CapabilitiesStateIntegers` (the release rule's phrases, `cancel_fields`, `tombstone_max_s`, "cannot start subsequently"; what a finished state promises and assumes, the inclusive hold, the 409's own `not_after_ms`); the review of these fixes: `.FinishedAttemptsAreHeldThroughTheirNotAfterMs` (the reviewer's p1 — retention 1 s, the replay at 1.3 s — and p1b — retention_count 1, one later finish — refused until `not_after_ms` + skew inclusive, then expired; held attempts fill the tombstone table; without `not_after_ms` dropped as before; the cap from the finish and a refused copy's extension; retention 0 keeps nothing), `.TombstonesHoldThroughSuppressedUntilInclusive` (skew 0, a copy at the wall clock's `not_after_ms` after the steady hold passed is refused), `.ARefusedCopyIsJudgedByItsOwnNotAfterMs` (no `not_after_ms`: not covered, `no_not_after_ms`, and it runs once re-sent after the hold); `GraphletAttempt.IdsAreValidated` (`expect_server_instance`); `GraphletServer.InventoryTableMatchesTheLoadersTypes` (every annotation type `initialize_annotation` builds is in the inventory's table, row-diff anchors exactly for `RowDiff<ColumnMajor>`, headers exactly for a `MultiIntMatrix`), `.SymlinkedMainFilesDoNotHideTheirSidecars`, `.DeliveryChecksEvery64KiBOfALargeToken` (the reviewer's 16 MiB probe: ≥ size / 64 KiB checks writing and assembling, bytes unchanged, the gap measured); `Graphlet.IndexIdentity` (the `.bloom` and the mask in the inventory, a manifest not covering one refused naming it; entries `a/x.seqs` and `b/x.seqs` refused for their shared base name; a `.seqs` the pair does not load refused); `RowDiffPathCache.RetentionKeepsRowsAndCostsAndBoundsTheCopies` (keep-all, 4/1, 16/8, requested rows only and all-narrow rules: the default and budget-aware decodes' rows and costs unchanged; copies per call ≤ n × (successors + 3) + stored / checkpoint; keep-all copies every reconstructed row; the rule copies fewer), `.FlatTupleRowsRoundTrip`; `LabelOracleBudgeted.LookaheadRunsFollowTheWalk` (the reviewer's pacing probe: warm's decoder charges 31,297 against 161,679 for sorted runs, 28,736 for the walk-ordered runs alone; query and recorder); `MiniRefSeq.PathCacheKeepsTheResponse` (+ three retention rules, 168 comparisons), `.PathCacheWarmSeedIsTheColdSeed` (a seed after a copy of itself equals the seed alone, every budget, batch, chunk and cache capacity), `.DeadlineRecordNamesTheLongestPiece`, `.ResolveDecodesEachRowOnceInBoundedBatches` (the review of these fixes, finding 7: on the query and on one repeating two thirds of it, column and header labels, presence and trace, batches of 1 to 4,096 rows × byte targets 1 B / 20 kB / 64 MiB × 0 / 30 kB / 256 MiB kept (explicit labels keep none) — 384 runs: the profiles equal the reference and the same labels given explicitly, no read beyond its batch, the first 64 rows, each at most twice the one before, every distinct row decoded once when the repeats can be kept and always for explicit labels), `.HeaderHitsFilterRequestedRanges`; `WalkerDeadlineChunks.*` on a virtual clock (exact bounds; the deadline record of a paced stop). Integration `test_retention_settings_are_validated_at_start_up` (the reviewer's negative-retention probe and thirteen more values, the p3 refusals of the review of these fixes among them, refuse to start; the largest accepted are stated), `test_a_cancel_naming_not_after_ms_covers_a_half_uploaded_request` (the reviewer's upload probe on a raw socket), `test_expect_server_instance_refuses_another_process` (server and CLI), `test_inventory_is_mirrored_by_index_manifest_py` (masked graph with a Bloom filter, row_diff, coordinate annotations with and without `.seqs`, `.coords`, `--no-coord-mapping`, symlinks; the tables equal), `test_symlinked_main_files_do_not_hide_their_sidecars` and `test_a_loaded_bloom_filter_is_part_of_the_identity` (the reviewer's `/tmp/metagraph-pass5-identity` cases, rebuilt); the review of these fixes: `test_a_manifest_names_one_bundle_and_its_loaded_files` (the reviewer's dup probe — a directory manifest of A/ and B/ refused by single-index servers, the CLI and `--verify`, naming the shared base name — and its C/ and `--no-coord-mapping` cases refused naming the `.seqs`), `test_retention_settings_are_validated_at_start_up` (each refusal names its option's own range). Byte identity against the previous build (f667d775): the unbudgeted and budgeted real requests on mini_refseq and UHGG, decompressed | §5, §6.8, §8.2, §8.4, §10.3 |
| T53 | the review of levels 4–5 (level 5 kept) | `RowDiffPathCache.ARefusedRowIsNeverCopied` (the reviewer's probe: a 1 KiB bound, its shared room 1 KiB, a row of 1,048,576 columns — no byte allocated by the insert nor at the moment the room is read, by jemalloc's thread counters; bytes, peak and kept rows 0, `rows_refused` 1; a tuple row of 65,536 columns likewise), `.AKeptRowIsCopiedWithinTheBound` (two rows of 1.5 MiB under 4 MiB and a third: the heap the cache holds stays within the bound during the insert, jemalloc's peak; the stated peak is what it held); `WalkerSeedPhase.SetupPiecesCoverTheSeedPhase` (2,500 headers with a 1,024-character prefix on a virtual clock that resolving a name advances by 1 ms: the longest piece is `setup`, exactly 2,500 ms, and `seed_phase_ms` 2,500 ms, at radius 5 and at radius 0 (no checkpoint); 10 seed and 20 extra labels: two setup pieces, 20 ms stated; a duplicate after 2,500 names: the failed seed's phase and piece stated), `.DuplicateChecksRefuseAsBefore` (the hash sets refuse what the scans refused, in the request's order, an unknown name before a later duplicate, extra labels against kept seed labels and earlier extras, a column and its header distinct targets); `GraphletAttemptRegistry.AFinishedStateIsReplaySafeOnlyWhenPinned` (the reviewer's finish-restart-replay probe on two registries: held by its process, the unpinned replay runs at the restarted one, the pinned one is refused `instance_mismatch`; the release rule's new phrases); integration `TestTraverseSeedPhase.test_the_longest_piece_covers_the_seed_phase` (the reviewer's deadline probe, CLI and warm server: the longest piece at least the names' resolution and a sixteenth of `seed_phase_ms`), `TestTraverseAttempts.test_the_registry_over_http` (the reviewer's attempt probe on a real server restarted on the same port: 200 with the same work units unpinned, 409 `instance_mismatch` pinned). Untimed bytes against the previous build (7aaee760) | §5, §6.8, §8.4, §10.3 |
| T54 | record coordinates (feature level 6) | Walker: `WalkerCoordinates.*` (exact intervals on both arms against a per-base scan, a split's partition, switch-entered runs and their lower bounds, zero-length runs, cap cuts with true counts, column kind, `chains_ended` and its inheritance, work and walk independent of coordinates, memory stops no deeper, denials and attempt stops keeping the side table). Request and JSON: `GraphletCoordinates.*` (the 400s, null reasons in their order, the echo resubmittable, one block in every detail, the K record only where a list was cut and the stripped response equal to the opt-out one, the delivery model bounding the block at 20-digit positions, `drop_coordinates` on memory stops), `.CoordinateTextIsExactAndBoundedByItsAccount` (142 results: the text less `coordinates_text_bytes` is the stripped result's, a graphlet's K and Q tokens and its counts' digits included; the account less the coordinate share is the opt-out walk's; the share ≥ 12 × the text), `.CoordinateAccountBoundsItsText` (each part at its widest: 14.9 the least ratio), `.AttemptsMeasureTheSameRatioWithCoordinates` (the server's path on 1.3 MB seeds: the attempt's measured ratio equal with and without coordinates), `.CapabilitiesBlockFollowsTheIndex` (the probe's block per index: `supported` = `supports_trace`, `kinds` with and without a CoordToHeader and none without coordinates, the cap's default, the record-end sentence; never in the per-request capabilities); `GraphletAttempt.ReserveCountsCoordinateText`, `.CoordinatesLeaveTheServersRatioUnchanged`. Real index: `MiniRefSeq.CoordinatesAgainstTheSourceRecords` (39 cells × caps 1, 16, "unlimited": the blaNDM gene, its reverse complement, three 200-bp windows and two repeat windows under limits 0, 2, the exhaustive preset and the switch cells of every recorded switch request's seed — constant 0.5 and 1, loss budget 2, limits 0 and 2, to their 3,000 bp, a run of each seed's ending past 1,000 bp —, column labels with and without a switch to 3,000 bp: 9,498 runs, 5,508 switch-entered, 48 lower bounds (27 strictly below the string count), 10,696 occurrences whose bases are the run's, 326 lists cut, 3 column intervals in a record's last k − 1 bases read as the next record's), `.CoordinatesKeepTheResponseUnderThePathCache` (144 comparisons over budgets, batches, cache capacities and chunk targets), `.CoordinateShareIsExact` (60 results, every detail and cap: the opt-in text less the coordinates' is the opt-out text). Measurements (`benchmarks/traversal/bench_traversal_measurements.cpp`, DESIGN §18.8): `BM_TraversalCoordinatesDepthAtTheStop` (M1, the D4 gate on mini_refseq, with what completing takes and the column cells under the refseq33m projection), `BM_TraversalCoordinatesDepthByRegime` (M1 by regime: the column fixtures of `scripts/traversal/make_column_coord_fixtures.sh`, 16 chains a run, and the wide fixture's column and header labels), `BM_TraversalCoordinateReserveOnTheWideFixture` (M2, `scripts/traversal/make_wide_coord_fixture.sh`). Integration `test_api_coordinates`, `test_api_capabilities` (the probe's block, not in the per-request capabilities, `max_uninterruptible_ms` null), `test_api_server_capabilities` (`coordinate_account_per_text_byte` 12, an integer), `test_multi_probe_describes_one_pair`. Byte identity against level 5 (67bef367): opt-out requests identical apart from the digit, timing and R21 (4) | §5, §6.8, §7.0, §7.1, §7.5, §10.3 |
| T55 | the first parent of a merge (R21 (4)) | `WalkerTest.MergeIsSpelledThroughTheParentWithTheMostLabels` (three graph types × modes: the majority route (two labels) through either branch is the first parent, in constrain and annotate mode, the path spelled through it, the minority's label routed at the merge, losses unchanged); `WalkerTest.MergeKeepsTheArrivalOrderWhereTheBudgetChargesTheChain` (three graph types × modes: with the chain charged (`tree` / `full`) under a memory budget the route that arrived first is first whichever labels it carries; without the budget or without the chain charge, the majority); `MiniRefSeq.AnnotateMergeRanksParentsByTheirOwnSegment` (the annotate rule pinned as it stands on a nested merge of the real cache's `mini_win200_02__annotate_merge`: the 5-bp merged parent of 6 labels first before the 32-bp one of 5, its walk upstream through 4 labels, the path's continuation carried by label 3); the CLI fixtures `merge` (the 62 merge's majority arrived second: its parents and partition swapped, a.fa now routed at 62, b.fa at 36 only) and `merge_ties` (the same locus without c.fa: ties keep the arrival order and show `R:route_split`); the byte-identity differences against level 5 classified (DESIGN §11: every differing request is `on_reconverge: merge`, and equal once parents, partitions, chains, continuations, run tables and the end labels' run / route_bp are canonicalised; no budgeted stop moves) | §7.1 |
| T56 | shutdown on SIGTERM | Integration `TestTraverseAttempts.test_sigterm_stops_the_server_promptly` (an idle server exits 0 within 5 s; one walking a three-seed attempt stops it at its next poll — nothing written, the connection closed — and exits 0 well within the drain time), `.test_sigterm_while_the_index_loads` (a single-index server whose graph is a FIFO nobody writes, so its load never ends: `ready: false`, 503, and exit 0 within 5 s of the signal) | §10.3 |
| T57 | D3: a derivation the time budget cut after part of the seed (the owner's decision, feature level 6) | `WalkerDerive.PartialDerivationDeliversTheSetOfTheKmersRead` (virtual clock, 1 ms a row, 70 ms: j = 65 of 110 k-mers — the first window consumed, the second read whole; 70 when later windows followed `batch_kmers` 1 — the superset {A, C} against the whole seed's {A}, arms truncated at 0, under `kmer` and `trace`, no coordinates, `partial derivation`), `.PartialSetThatWouldFailTheSeedFailsAsBefore` (`exhaustive` over the cap, an extra label the superset duplicates, an ambiguous derived header, a depth-0 memory failure: each fails with the time budget "after 65 of 110 k-mers", as before D3; a cap that cuts the superset states it as a superset), `Walker.PartialSetWithAnUnrepresentableNameFailsAsBefore`, `Walker.FailedDerivationIsAStatedOutcome` (the server-lowered budget: walked, `observed` 1), `WalkerDeadlineChunks.DerivationWindowIsPaced` (j = 0 fails), `Graphlet.PartialDerivationIsPricedAndBound` (the limitation priced and bound in every detail); integration `test_traverse_partial_derivation_states_its_set`; the fixtures `coords/d3_kmer`, `d3_trace_coordinates` and `d3_memory_stop` (the library's `WALKED_QUALIFIED_CLASS`). The `derivation` limitation's `observed` for `time_budget` has two units — elapsed ms failed, j walked — stated in §7.0 (*review of 2026-10-06, X1*) | §6.1, §7.0 |
| T58 | the review of 2026-10-06, P2 items (feature level 6) | **W1** `WalkerDerive.WindowIsSixtyFourKmersWhateverBatchKmers` (the reviewer's U01-01 shape: 992 columns, a 195-k-mer derived seed, 100,000 work units — the response at `batch_kmers` 1, 7, 63, 65 and 1,000 equal to the one at 64, a work stop that fails the seed; an unbudgeted `no_carrier` at k-mer 90: `work_seed` and the response equal at 1, 7 and 64). **W2** `.CoordinateStateIsTheSumItReplaced` (16 carriers shrinking to 8 and 4, two compactions, on a column annotation and on a row-diff one: the running total recounted by `WalkerHooks::recount_derivation_state` wherever the state is counted — the soft observation under a memory budget of 4 GiB or 1 MiB, and on the row-diff annotation also `beside`, before every window's budget-aware read, under those budgets and under a work budget alone; without a budget nothing counts it, so that case only checks that the hook changes nothing — the result equal to the run without the recount; a debug build asserts it always; checked to bite with a one-byte drift after a compaction, the row-diff work-only case included) and `.TraceDerivationUnderAMemoryBudgetIsLinear` (200 kbp, 4 carriers: under a budget that never binds within three times the unbudgeted time; the base took 4.1–4.4 s against 0.3 s on the CLI; skipped in a debug build, whose assert recounts the state at every observation — the quadratic cost measured). *(The review of the P2 fixes: this row said the recount ran "without a budget"; GCC rejected an unbraced `EXPECT_EQ` in the first test, `-Werror=dangling-else`, so the changed test files are now compiled with GCC 13 too.)* **W3** `WalkerDeadlineChunks.LookaheadChainsReadTheStops` (a 20,000-k-mer chain on a clock every reading advances by 10 µs: a cancel and a walk-until at 50 ms, unpaced and paced, end the walk by 55 ms where the chain alone takes 200 ms; no head piece longer than 5 ms; `DecodePacer::max_head_ms` in the observation) and `.LookaheadHandsTheWalkUntilToTheAttempt` (the same chain under a control that reads the clock on every 8th poll and on every forced one, as the server's attempt does: unpaced and paced, the walk-until is taken by the lookahead's clock-reading poll and no checkpoint after it answers "no stop"; it fails with the earlier unpaced `ms_left` shortcut). **C20** `GraphletAttempt.WalkUntilStatedIsTheLowestAClockReadingPollSaw` (the reviewer's U11-02 driver: stride 8 states 65,000 where 58,900 was in force between clock-reading polls, stride 1 states 58,900; the `bound` text names the poll that reads the clock). **Texts** `GraphletAttemptRegistry.ReleaseTextsStateTheOverrunAndTheGrounds` (the X2, LRG-R1/R2, C16, C20, C24, C30 and X4 phrases, with and without retention; no "uninterruptible step" in `bound`, `not_after` or `release_rule`, no "first poll after" there or in `delivery_reserve.rule`). The reviewer's reproductions re-run before and after (U01-01, U01-04, U03-01). Byte identity against `6897db99` on the CLI fixtures (`graphlet_fixtures.py --from-cli --check`) and the bench capabilities: only the texts and W1's stated cases differ | §5, §6.1, §6.8, §7.0, §10.3 |
| T59 | the review of 2026-10-06, P3 server items (level-6 corrections) | **C10** `GraphletAttempt.UnstartedSeedsStateLabelsFromSeedAsTheWalked` (explicit, derived and annotate seeds cancelled after the first: the not-started results state the walked one's value, false in annotate mode); integration `test_wide_index_delivery_reserve_stops_the_walk` (both label modes) and `test_cancel_mid_walk_with_the_client_connected` (annotate: false). **W9** `GraphletStage3Review.LabelsThatDoNotFitAreNotARowStop` (on a column annotation a row naming 300 labels of 4 KB at 2 MiB: the labels' message, "observed: at least", `lower_max_labels_per_node` and `label_constrained_query`; one label per node walks past it) and `GraphletStage3Decode.StatementsFitTheirWidths` (its message at the widest values). **C8** `GraphletAttempt.StopLatencyIsMeasuredWhereTheWalkEnds` (the reviewer's timeline: 400 ms kept after a second call 10 s later, abandoned or failed; a walk that ended before its walk-until measures nothing). **C27** `GraphletServer.PeerClosedSeesTheServersOwnShutdown` (the server's own shutdown behind a CRLF and behind a pipelined request, and with nothing waiting; a pipe; a closed descriptor). **C26** `GraphletServer.CompressionTakesTheTextInPieces` (pieces of 1, 7, 32,768 and 65,537 bytes and one short of the text inflate to the whole text in both containers under the check; a text of one piece gives the old bytes; the output change itself, a text of 4 GiB or more, by the reproduction: 2^32 + 10 bytes inflated to 10 at `ea285c2e`). **C9** `GraphletServer.AssemblyReservesTheResponseOnce` (the response's capacity within 16 bytes of its size, members on either side of `results`, with and without checks; the bytes those of the tree). **W11** `LabelOracleBudgeted.LookaheadRunsFollowTheWalk` (a cache of 8 keys: the walk's first run decoded and kept, 130 charges where every run's were at least 28,736; under the byte bound one run kept, the second dropped, the warm ended) and `MiniRefSeq.LookaheadKeepsItsRunsUnderAMemoryBudget` (the review's walk under 16 MiB at `batch_kmers` 2,048, 4,096 and 8,192: the walk of the work budget alone, its tuple rows at most 25% more than that walk's — 33,015, 33,053 and 30,667 against 27,642; `ea285c2e` 36,275, 36,838 and 38,441). **W16** `ResolveCoordTest.ExplicitLabelsGiveTheDiscoverysProfiles` (301 explicit labels, presence and trace, column and row-diff annotations: the discovery's runs, trace breaks and counts). The reproductions re-run before and after: U11-01, U03-03, X-DUP-01, U13-02 (Linux in a container and macOS), U13-04 (2^32 + 10 bytes inflate whole), X-EFFICIENCY-02 (the delivery peak), U05-01 and X-EFFICIENCY-04. Byte identity against `ea285c2e`: the CLI fixtures (`graphlet_fixtures.py --from-cli --check`: up to date), the coordinate snapshots (only time values differ), 238 replays of the mini_refseq bench requests (with memory budgets at `batch_kmers` up to 8,192): no untimed byte differs; both capabilities documents identical | §6.8, §7.0, §10.3 |
| T60 | `/resolve`'s deadline (milestone 1b, §4.5; no feature level) | `ResolveCoordTest.DeadlineStopIsTheResolveOfAPrefix` (decision B7: the deadline passing at every reading of the work in turn — rows read three at a time, a query with a stretch absent from the graph and a repeat, a discovery with and without truncation and 120 explicit labels, presence and trace, a column and a row-diff coordinate annotation — gives exactly the resolve of the prefix of `resolved_kmers` k-mers, a k-mer in the graph, never decreasing; no stop: the unbudgeted profile), `ResolveTest.DeadlineStopsTheExplicitHitsPass` (the hits fetched 4,096 k-mers at a time: a stop before each piece, at the next k-mer in the graph, on the direct path and the row paths, prefixes again), `Resolve.FinishCheckIsReadInTheLoopsOverTheLabels`; `ResolveDeadline.*` on the route with a clock that steps 1 ms a reading: `.WithoutTheFieldTheAnswerIsAsBefore` (no `limits`/`stop`; a budget that does not run out gives the same answer), `.AStopAnswersThePrefixWithTheStopBlock` (stops before the rows and after 64, 192, 448, 960, 1,984 rows, each byte for byte the unbudgeted resolve of the prefix apart from `limits`, `stop`, `timing`), `.ReserveOverrunAnswers503` (reserve 0: `ResolveDeadline`, the body `{error, code: deadline}`; the delivery check before and after the budget), `.BudgetValidationAndTheCap`, `.AnExplicitSelectionPastThePrefixIsNotMade`, `.RequestErrorsAre400WhateverTheStop` (review of 2026-10-07, V1-02), `.CapabilitiesBlock`; `.WithoutTheFieldTheAnswerIsAsBefore` also asserts that a request without the field reads no clock for a deadline (T3-05). Integration `TestTraverseAPI.test_api_resolve_deadline` (server and CLI: a budget just above the reserve stops at once with the stop block, a 400 at the reserve, the clamp to `max_time_ms`, the `resolve` block on both GET routes); `test_pattern.TestPatternMini.test_resolve_deadline_is_503` (the HTTP 503 end to end, with and without gzip: no `Retry-After`, no `Content-Encoding`, the exact body; the CLI's exit 1; V1-04); `test_traverse.TestTraverseResolveRegression` (byte identity of the requests without the field against the previous build, server with and without gzip and CLI, `timing.elapsed_ms` apart: it runs only with `$METAGRAPH_BASE_BINARY`, so not in CI; X-TESTS-05) and the bench mini panel | §4.5, §10.3 |
| T61 | the capabilities' rules as references (§10.3, "The rules are references"; no feature level) | `GraphletAttemptRegistry.CapabilitiesRulesAreSpecReferences` (each attempts rule the reference of the table, whatever `retention_s` and `poll_stride`, printable ASCII); `GraphletCoordinates.CapabilitiesBlockFollowsTheIndex` and `ResolveDeadline.CapabilitiesBlock` (the coordinates' and `/resolve`'s references); `test_traverse_capabilities_references.py` (every fixture server's `/traverse/capabilities` and `/capabilities` states the table's reference in each rule field, no other text above 100 characters; each referenced section exists, holds the reference's named parts and states the rule's substance, phrase by phrase — the phrases the texts' tests checked —; a section that does not exist, a part or a rule it lacks, another section's reference and a text written out again each fail); integration `test_api_capabilities`, `test_api_resolve_deadline`, `test_api_server_capabilities` (`deadline_check.poll_stride` 8, an integer) and `test_capabilities_only_gain` (the texts compared as present, `poll_stride` the one field added) | §10.3 |

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
22. **Record coordinates are opt-in output, not in the MGT body** (feature level 6, §7.1; DESIGN §18):
    the design's field, with the owner's decisions C1–C12 and C-N1..C-N9 — in no outcome
    class (C2), a cut stated by the block (`complete`, `occurrences_total`) and its limitation, a switch-entered
    run into a still-live label marked `lower_bound` rather than walked differently (C-N2, X-18.2), column
    positions in the column's k-mer index space without the record lengths (C10: a record's last k − 1 bases
    share their numbers with the next record's first), the K kind `coordinates` and the action
    `drop_coordinates` as free tokens (C12, X-C12). The cap has no server maximum (D6: the memory budget and the
    coordinates read bound the block; set `--traverse-max-memory-mb` in deployments, ops note O2).
23. **A merge displays its majority parent** (R21 (4), §7.1): the first parent of a merged segment — whose chain
    the displayed walk, its continuation and `route_bp` follow — is the one carried by the most labels, not the
    first to arrive; ties keep the arrival order. Before level 6 it was always the first to arrive, and it still is
    in `tree` and `full` detail under a memory budget, where the displayed chain is charged and would move the
    stop (the owner may trade that for the rule everywhere, DESIGN §19).

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
