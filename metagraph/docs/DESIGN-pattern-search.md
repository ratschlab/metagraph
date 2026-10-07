# Design: pattern search — motifs shorter than k, IUPAC patterns and peptides, counted first, found completely

**Status:** v6 (2026-10-07), **approved for phased implementation** by the sixth external review (readiness,
§21), after five design reviews (v1: nine findings, §16; v2: eight and a resource policy, §17; v3: seven and
four contract corrections, §18; v4: six findings and three corrections on the predicate extension, §19; v5: six
findings and one naming decision, §20; all confirmed against the code). Written from a read-only check of the code at
`804731aa` (gr/labeled-traversal). Nothing is implemented. Line numbers are pointers, not anchors.

Companion documents: `NOTE-mismatch-search.md` (bounded mismatches: a later increment on the same engine, §12),
`SPEC-labeled-traversal-core.md` (`/resolve` and `/traverse`, whose conventions for caps, truncation statements,
determinism and budgeted annotation access this route follows), `DESIGN-traverse-graphlet.md` §14 (the guarantee
rule and the budgeted reads).

## 1. Problem

No server route answers a query shorter than k. `/search` returns nothing for it (every `AnnotatedDBG::get_*`
answers `{}` for a sequence without a full k-mer), `/resolve` refuses it ("Query shorter than k",
`resolve.cpp`), and `/traverse` rejects the seed (`walker.cpp`, `validate_seed`). An agent with a 16-nt primer,
a 20-nt spacer, a transcription-factor site written in IUPAC codes, or a 9-residue peptide has to invent a longer
seed, which the index may not contain. A query with a degenerate base (R, Y, N) finds nothing on any route: the
k-mer lookup encodes it as an out-of-alphabet symbol and gives up (`dbg_succinct.cpp:319-322`).

The navigation exists. The BOSS index answers "which nodes end with this string" by range narrowing
(`BOSS::index_range`, `tighten_range`, `boss.hpp:474-483`), the alignment seeder uses it for sub-k seeds
(`SuffixSeeder`, `aligner_seeder_methods.cpp:153-320`), and `suffix_to_prefix` (`:95-139`) already walks a range
one symbol at a time trying every real symbol: a degenerate position is exactly that. What the aligner lacks is
completeness (it keeps best hits, by design), an occurrence model that says where in a record a match lies, an
output that says which labels carry it, budgets on the annotation reads, a count that is known before any row is
decoded, and a server route.

## 2. Goals and non-goals

Goals:

1. **One route, `POST /pattern`**, that **counts first** and, when the aggregate count is admitted, finds
   **every** occurrence of a pattern in the retained, k-mer-covered sequence of the index, and reports per label:
   presence, the number of distinct graph contexts, and (on BASIC indexes with record mapping) the record,
   1-based position and strand of each occurrence.
2. **Three pattern kinds on one engine:** an exact DNA motif of any length (1 ≤ L, including L < k); an IUPAC
   pattern (a set of bases per position); a peptide over all synonymous codons, as a codon automaton.
3. **Both strands**, each oriented pattern searched completely and independently; strand stated where the index
   can state it.
4. **Every count carries its unit and its relation** (`exact`, `at_least`, `bounds`, `unknown`); nothing decoded
   is dropped silently; absence is claimed only by an answer whose retrieval completed.
5. **Reuse without changing what the aligner does.** `/align`, `/search` with `align: true`, and `metagraph align`
   keep their outputs byte for byte (§11).
6. **Bounded like `/traverse`:** annotation reads through the budget-aware, paced traversal classes, one deadline
   with a reserved finalisation window and per-shard memory accounts per request, every stop stated with its
   phase (§5).
7. **Served synchronously:** the server route is a plain request/response like `/search`; the Python client gets
   a method beside `align()`; the search service exposes it as a synchronous MCP tool and REST route, modelled on
   `traverse_resolve` (§9).
8. **Annotation predicates** (§5.6): a broad pattern narrowed by a logical condition on the labels of each
   context, evaluated before anything is materialised, with the retrieval threshold applied to what passes.

Non-goals (v1):

- Mismatches and indels (§12; `NOTE-mismatch-search.md`).
- A tolerant best-hit search for divergent variants: that is the aligner; IUPAC rows in its scoring table on
  `/align` are a separate item (§12).
- Verified placement on CANONICAL or PRIMARY indexes (§3: strand and offset are ambiguous there; v1 reports
  presence and contexts only).
- The `suffix` scope on PRIMARY indexes wrapped in `CanonicalDBG` (§4.1: a suffix of a virtual node is a prefix of
  a stored one; `any_offset` covers it).
- Loading an index inside a request (§5.3: resident indexes only).
- Protein graphs (the protein alphabet build): the engine is for DNA graphs; a peptide is searched as DNA.
- Archive-wide batch runs through the job queue: the sync path first; the queue can wrap the same route later.

## 3. Semantics

**Pattern.** A pattern is a sequence of L positions; each position allows a set of bases from {A, C, G, T}. For a
peptide the allowed set at a position depends on the bases already spelled within the current codon (§6); the
engine's one primitive is therefore `allowed(position, spelled_prefix)`. A pattern's N means {A, C, G, T}; it
never matches a record's N symbol (DNA5 builds, `alphabets.hpp:81`); the unconstrained flanks of §4.1 do match it.

| kind | input | positions |
|---|---|---|
| `dna` | a string over A, C, G, T | one base each |
| `iupac` | a string over the 15 IUPAC codes | the code's set (R = {A, G}, N = {A, C, G, T}, …) |
| `protein` | a string over the 20 amino acids (no stop, no ambiguity codes in v1) | 3 per residue; the codon automaton of the residue (§6) |

A pattern's **information** is Σ log2(4 / |set_i|) bits over its positions (for a peptide, log2(64 / codons) per
residue). For a pattern longer than k the **anchor window's** information, over positions [0, k), is stated
separately: a highly informative suffix does not make an `N…N` anchor cheap.

**Scope.** A request names what it counts and finds:
- `suffix` (L ≤ k): the k-mers whose last L symbols instantiate the pattern — every occurrence that starts at
  position ≥ k − L within its **retained island** (the stretch of a record between its start, unsupported
  symbols and pruned k-mers that the index kept, see "Covered sequence"), each once per distinct k-mer. The
  cheap, exactly countable scope; BASIC and native CANONICAL indexes only (§4.1).
- `any_offset` (L ≤ k, the default): every k-mer that contains the pattern at any offset p ∈ [0, k − L] — every
  occurrence inside a retained k-mer, at the cost of discovering the flank ranges (§4.1).
- `long` (L > k, implied): anchors instantiating positions [0, k), then completed extensions (§4.2).

A `suffix` answer states `absence_scope: suffix_only`: a label absent from it may still carry the pattern inside
the first k − L bases of a retained island.

**Covered sequence.** The index represents only the k-mers it retained: a record is split at every symbol the
build did not admit (N on a DNA4 build), a stretch shorter than k yields no k-mer, and graph cleaning may have
pruned k-mers anywhere. A **retained island** is a maximal run of consecutive k-mer start positions whose k-mers
the index retained; the base spans of two islands may overlap (with k = 3, keeping `ACG` and `GTA` of `ACGTA` but
pruning `CGT` leaves every base covered and two islands). Every completeness statement of this document is over
retained islands, never over records as deposited: with k = 5 on a DNA4 build, `AAAAANACNCCCCC` contains `AC` in
a stretch of two bases that no k-mer covers, so no route can see it; and the `long` pattern `ACGTA` above has no
retained path although every base is covered. A pattern of L ≤ k is found wherever one retained k-mer contains
it; a pattern of L > k is found only where **all** its L − k + 1 k-mers lie in one island. The capabilities state
the build's alphabet and whether the graph was cleaned.

**Graph context.** A graph context is (orientation, k-mer path, offset): a path of n = max(1, L − k + 1) k-mers
and, when L < k, the offset p at which the pattern sits inside the single k-mer. Two contexts are the same iff all
three agree. A k-mer that contains the pattern twice (`ACGAC` contains `AC` at offsets 0 and 3) gives two
contexts. A graph context is **not** a physical occurrence: a k-mer that occurs ten times across the records of a
column is one k-mer in the graph, and repeated copies of `AAA` in one record share one context. Occurrence counts
exist only with record mapping (below). A context is the unit of a result, and a result always carries the
context's k-mer and the pattern instance it matched (§7.2): the row of the annotation is identified by its k-mer,
and nothing else identifies it where no placement exists.

**Occurrence.** On a BASIC index with record mapping (a `CoordToHeader` `.seqs` file, loaded at start-up unless
`--no-coord-mapping`, `load_annotated_graph.cpp:33-60`), a placed occurrence is (column, `seq_id`, 1-based start,
strand), with `seq_id` the `CoordToHeader`'s stable record index inside the column (header text need not be
unique and is carried for display). It is derived from a context by mapping the context's first k-mer coordinate
to (`seq_id`, local position) **first** and adding the offset p **after** (§4.3). Placed occurrences are
deduplicated by their own identity, never by node id (§5.4).

**Counts and their relations.** Every count in the answer is `{value, relation, unit}` with
`relation ∈ {exact, at_least, bounds, unknown}` (`bounds` carries `lower` and `upper`) and
`unit ∈ {graph_contexts, anchors, paths, placed_occurrences, labels}`. `exact` is answered only for a count whose
discovery completed (no step, time or memory stop touched it). A stop turns the count of the interrupted phase
into `at_least`, the sum of the lower bounds of everything explored; it is `bounds` only when **every**
unexplored portion of the search — every IUPAC branch, offset and strand not yet visited — has a valid upper
bound, which in practice means that all ranges were discovered and only mask scans were interrupted (§4.1); a
range's own bounds otherwise stay diagnostics under `work`. A phase not started is `unknown`, and so is every
phase after a stop. A count is never promoted by assumption: an anchor count for `long` says
nothing about completed paths (`AAAC` can have an `AAA` anchor and no path), and a graph-context count says
nothing about placed occurrences.

**Label.** A label carries a context when it annotates every k-mer of the context's path. Per label the answer
gives `contexts` (distinct graph contexts in the requested scope), `anchors` and `paths` for `long`, and, with
record mapping and admitted retrieval, `occurrences` (placed, deduplicated, exact when retrieval completed).

**Strand.**
- On a BASIC index (k-mers stored as deposited, e.g. refseq33m), contexts of P are on the record's + strand and
  contexts of rc(P) on the − strand; a − strand occurrence of P occupies the same record interval as the + strand
  occurrence of rc(P) it is found as. Both oriented patterns are searched. A **palindromic** pattern
  (P = rc(P)) is searched once; its contexts carry `strand: "="`, count once in every total, and the per-strand
  split reports them under `both`, never under `+` and `−` separately.
- On a CANONICAL or PRIMARY index a stored k-mer may be the reverse complement of the deposited one, so neither
  the strand nor the offset is known: an offset p in the stored k-mer is deposited offset p or k − L − p. v1
  reports presence and contexts only on such indexes (`placement: "none_canonical"`); the alignment code's own
  rule is the same (`AnnotationBuffer` drops coordinates outside BASIC mode, `annotation_buffer.cpp:27-29`).

## 4. The engine

Two phases over the succinct (BOSS) graph; the second only for `long`. BOSS facts the design rests on
(`boss.hpp`): a node is a (k − 1)-mer, an edge is a k-mer (node + its last symbol `W`); the edges leaving one node
are contiguous, the group ends marked in `last`; `index_range` narrows to the edges whose source node ends with a
string of at most k − 1 symbols (:474-483), and a range is normalised to whole node groups with
`pred_last(first − 1) + 1` and `succ_last(last)` (:398-404; the seeder does this at
`aligner_seeder_methods.cpp:303`); `tighten_range(&first, &last, s)` appends one symbol (:684-685, by `rank_W`);
`rank_W(i, c)` counts edges with `W = c` up to i (:416); `W` stores an edge's symbol c either plain or **marked**
as `c + alph_size` when another edge into the same target carries the same symbol (`get_minus_k_value`, :280) —
both are k-mers ending in c; `pick_edge(node_edge, c)` selects a node's edge with symbol c (:365). On top of BOSS,
`DBGSuccinct` keeps `valid_edges_` (`dbg_succinct.hpp:124-147`, the `.edgemask` file), a `bit_vector` with rank
and select that is 0 on every dummy edge (those containing `$`, sources and sinks) and on pruned ones. **The mask
is required:** without it `in_graph` treats every edge as a k-mer (`dbg_succinct.cpp:935`), dummies included, and
no count of this route is right. A server whose graph loaded without its mask states `pattern.mask: absent` in
the capabilities and refuses the route with `mask_required`.

### 4.1 Phase 1: anchors and contexts, by range narrowing (reused from the seeder)

**The lookup rule.** Let Q be an oriented pattern and m = min(L, k).

- Positions 0 … m − 2 are matched on **node ranges**: a depth-first search that starts from the whole edge range
  and at depth d tries only the symbols `allowed(d, prefix)`, dropping a branch whose range becomes empty. This is
  `suffix_to_prefix`'s loop (`aligner_seeder_methods.cpp:95-139`) with its symbol set made a parameter; its own
  callers keep today's set, every non-sentinel symbol of the graph's alphabet (`s = 1 … alph_size − 1`, which on
  a DNA5 build includes N). A leaf is a range of nodes whose suffix instantiates positions 0 … m − 2.
- Position m − 1 is matched on **`W`**: for each allowed symbol c, the k-mers whose last m symbols instantiate
  Q[0, m) are the valid edges of the leaf range with `W ∈ {c, c + alph_size}`. For L = k the leaf is one node and
  this is `pick_edge`; for `long` it yields the anchors of §4.2.

This keeps the (k − 1)-node / final-edge distinction the existing function keeps
(`call_nodes_with_suffix_matching_longest_prefix`, `dbg_succinct.cpp:308-390`: `index_range` over at most k − 1
symbols, then `pick_edge` for the k-th). That function is left as it is; the new entry points live beside it (§10).

**Counting.** For a normalised range [first, last] and a symbol c:

```
plain   = rank_W(last, c)              − rank_W(first − 1, c)
marked  = rank_W(last, c + alph_size)  − rank_W(first − 1, c + alph_size)
invalid = (last − first + 1) − (valid.rank1(last) − valid.rank1(first − 1))
count   = plain + marked − |{ e ∈ [first, last] : valid[e] = 0 and W[e] ∈ {c, c + alph_size} }|
```

The first three lines are O(1). The last term is zero when `invalid` is zero (the common case: a node range for
a pattern suffix contains a source-dummy node only when some record prefix ends with that suffix). Otherwise it
is a **scan** of the range's invalid edges through `valid.select0`, charged to `max_steps` one step per edge and
checked against the deadline like every other loop; it is real work, not a rank. With R = plain + marked, I the
range's invalid edges, s of them scanned and t of those matching, an interrupted scan leaves the range with
`lower = max(0, R − t − (I − s))` (every unscanned invalid edge might match) and `upper = R − t`. The count of
the whole pattern is `exact` only when every range's scan completed; when a scan was interrupted but every
range of every branch, offset and strand was discovered, it is `bounds` (the sums of the ranges' lower and upper
values); when discovery itself was interrupted, it is `at_least` over what was explored, since an undiscovered
branch has no upper bound (§3). The oracle suite (§13) pins the three cases the first draft got wrong: marked
edges (a motif with one plain and two marked matching edges), a sink dummy `AC$` inside a range, and a source
dummy `$$AC` leaving a range, plus an interrupted scan.

**All offsets (`any_offset`, L < k).** The k-mers ending with Q cover every occurrence at record position
≥ k − L. An occurrence inside a record's first k − L bases is contained in the record's first k-mer at some offset
p < k − L, and that k-mer is real and annotated; what is not guaranteed is a dummy chain to it (dummy edges are
pruned at build, `boss.hpp:202-239`), so no dummy pass can find it. The engine instead enumerates **every k-mer
that contains Q at any offset**:

- offset p with p + L ≤ k − 1 (Q inside the node): the nodes whose suffix is Q·X with |X| = k − 1 − p − L, found
  by continuing the DFS from the Q leaves with **every real symbol** for |X| more steps (the seeder's loop,
  unchanged: on a DNA5 build the flank admits N, so `ACNTA` at k = 5 is found for `AC` at offset 0); **all** valid
  edges leaving those nodes are such k-mers, counted by `valid.rank1(last) − valid.rank1(first − 1)`, which also
  excludes sink dummies;
- offset p = k − L (Q is the suffix): the `W` rule above.

Discovering the flank ranges is **branching work, even for an exact DNA pattern**, and it is not bounded by the
suffix count: with k = 5 and the record `ACGTT`, the motif `AC` has zero suffix anchors and one context at offset
0. The engine therefore counts the ranges it visits (`ranges_visited`) against `max_steps`, and the answer reports
them; no closed-form bound is claimed. Per offset the counts are exact once its ranges are discovered. This is
what makes the completeness claim hold without dummies (§5.1), at the price of context counts that are not
occurrence counts, which the answer states per unit (§3).

**Orientation.** Each oriented pattern Q ∈ {P, rc(P)} is searched completely on its own and extended in its own
reading direction from its own first k-mer. No anchor id is converted between orientations. `COMPL_TAB`
(`reverse_complement.hpp`) maps every IUPAC code correctly (A↔T, C↔G, R↔Y, K↔M, B↔V, D↔H; S, W, N to themselves),
and a codon automaton is reverse-complemented by reversing the residue order and reverse-complementing each
codon (§6).

Three graph modes, told apart by `get_mode()` and the wrapper:
- **BASIC** (`DBGSuccinct`, mode BASIC): one search per oriented pattern; every scope; placement possible.
- **native CANONICAL** (`DBGSuccinct`, mode CANONICAL, both orientations stored as k-mers): one search per
  oriented pattern on the base BOSS, no mapping; every scope; presence and contexts only.
- **PRIMARY wrapped in `CanonicalDBG`** (the server wraps at load, `load_antotated_graph.cpp:119-125`): one
  orientation is stored, so an instance of Q in a deposited sequence lies in a k-mer x that is stored either as x
  or as rc(x).
  - `any_offset`: if rc(x) is stored, it contains rc(Q) at offset k − L − p, so the searches of Q and rc(Q), each
    over all offsets, together find every retained instance. Complete for presence and contexts.
  - `suffix`: Q as the suffix of x is rc(Q) as the **prefix** of the stored rc(x), and a prefix is not a BOSS
    suffix range (runtime-verified by the review: PRIMARY stores `ACGA`, the wrapper exposes `TCGT`, and `CGT`
    is found neither as stored suffix `CGT` nor as stored suffix `ACG`). v1 does not offer `suffix` on wrapped
    PRIMARY graphs: the request answers `scope_unsupported: suffix_on_primary` and names `any_offset`; a later
    prefix lookup with orientation conversion may lift this (§12).
  - `long`: phase 1 runs on the base BOSS for **Q[0, k) and for rc(Q[0, k))** — the reverse complement of the
    anchor window, not the anchor window of rc(Q): with k = 3 and Q = `AACG`, rc(Q)[0, 3) = `CGT` maps to `ACG`,
    Q's last k-mer, whereas rc(Q[0, 3)) = `GTT` maps to `AAC`, Q's first. Anchors of rc(Q[0, k)) are mapped with
    `CanonicalDBG::reverse_complement(node)`, exact for a full k-mer (the seeder's mapping,
    `aligner_seeder_methods.cpp:303-309`, used here only for full k-mers). The two probes' anchors are **united
    by wrapper node id before anything counts, admits or extends them**: a palindromic anchor window is found
    by both probes and `reverse_complement` maps a palindromic node to itself (`canonical_dbg.cpp:536-539`),
    so with k = 4 the pattern `ACGTA` would otherwise count `ACGT` twice and extend it twice; an IUPAC window
    can have some palindromic instances and some not, and the union handles each anchor on its own. Extension
    then runs on the wrapper (`CanonicalDBG::call_outgoing_kmers`) along Q. For a peptide, the automaton's state at position k (residue
    index, position in codon, codon prefix) is recomputed from the anchor's spelled sequence, which instantiates
    positions [0, k) exactly.

### 4.2 Phase 2: extension beyond k (new, small)

For `long`, each anchor (a k-mer instantiating positions [0, k)) is extended one base at a time along outgoing
edges (`DeBruijnGraph::call_outgoing_kmers`, `sequence_graph.hpp:224`), keeping only the symbols
`allowed(position, prefix)` at the next position, by depth-first search to position L. Every complete path is a
context. The search stops at the first disallowed base, so it is complete and never scores anything.

Why not the extender: `DefaultColumnExtender` keeps the best DP cell per node and query position and prunes by
`xdrop`, `rel_score_cutoff` and `num_alternative_paths` (`aligner_config.hpp`, `aligner_extender_methods.cpp`).
With an exact-only scoring it would still merge distinct paths through one node and drop alternatives by design.
For an exact, IUPAC or codon pattern the plain DFS is complete and cheaper. The extender stays what it is; the
tolerant variant is §12.

Two thresholds, two admissions (§5.2): `max_anchors` admits the **extension** (an anchor count above it answers
`withheld: anchors_above_threshold` with the exact anchor count and no paths), and `max_paths` admits
**retrieval** (paths are counted as the DFS completes them, with `candidates_examined`, the branches entered,
beside them). `max_steps` bounds the DFS together with phase 1's `ranges_visited`.

A later optimisation (§12): when the anchor window carries little information (a leading `N…N`), anchor on the
most informative k-window inside the pattern and verify in both directions (`call_incoming_kmers` for the left
part). v1 anchors on [0, k) and states the anchor window's bits.

### 4.3 Labels, coordinates and placement (reused from the traversal, not from the aligner)

**Two steps, because the request names no labels.** The traversal's `LabelQuery` reads the rows of a key set
against a **permitted label set** given to its constructor (`label_oracle.cpp:454`) and cannot discover labels;
`LabelRecorder` discovers them (annotate mode: the labels of each key, ascending, at most `max_labels_per_node`,
`truncated()` when cut, `label_oracle.hpp:698-727`) but returns no coordinates. Retrieval therefore runs as:

1. **Discovery** with a `LabelRecorder` over the admitted anchors (every k-mer of every context; for `long`,
   every k-mer of every path): the budget-aware, paced fetch (`LabelRecorder::fetch` with `DecodeBudget` and
   `ReadPacing`, :727), a per-anchor cap `max_labels_per_anchor` (default 64, as `max_labels_per_node`). An
   anchor whose row was cut is stated (`anchors_truncated`, with the cap and the row's total); its label set is
   incomplete, so in `all_or_count` the pattern's results are withheld (`withheld: anchor_labels_truncated`)
   unless the request raises the cap, and in `partial` the cut is stated per anchor.
2. **Placement** with a `LabelQuery` built on the discovered label set, `with_coords`, over the same keys
   (`LabelQuery::fetch(keys, n, DecodeBudget&, ReadPacing*)`, :541-560): the coordinates of every (anchor, label)
   pair, all or nothing per fetch, interrupted or refused fetches leaving the cache and the budget as on entry and
   saying why. Only on BASIC indexes with coordinates; elsewhere step 1 is the whole retrieval.

Both steps charge the request's annotation account and read the deadline at their run boundaries. A paced
full-row/tuple visitor that does both in one pass (extracted from the two fetch loops) is the planned
refinement (§12); v1 takes the two existing classes as they are. This is integration work the earlier drafts
priced as reuse; the estimate in §13 carries it.

**Supported backends.** The budget-aware path exists only where `LabelOracle::decode_charged()`
(`label_oracle.hpp:347`) is true: the row-diff family with the budgeted decode (`row_diff_budgeted.cpp`;
RowDiff<BRWT>, the TupleRowDiff coordinate formats). On other backends (column, column_coord, brwt, row_flat,
row_diff_disk, …) the budget answers `UNSUPPORTED` (`decode_budget.hpp:79`). The route states the backend's
support in `/capabilities` (`pattern.annotation: budgeted | unbudgeted`); the `count` mode is available on every
backend; retrieval on an unbudgeted backend is refused with `annotation_unbudgeted` unless the request sets
`allow_unbudgeted_annotation: true`, in which case the read runs without a memory account, the deadline is
checked between keys only, and the answer says so (`annotation: unbudgeted`).

The alignment's `AnnotationBuffer` is not used: it materialises complete rows with no work, memory or time bound,
and one common anchor (a conserved 16S motif) can carry a row of a hundred thousand labels with coordinates. A
row the budget refuses is reported with its anchor (`rows_refused`), never skipped silently.

**Label consistency for `long`.** A label carries a path when it annotates every k-mer of it. With coordinates
the intersection is taken on consecutive coordinates (the `utils::set_intersection(..., dist)` idiom of
`LabeledExtender::call_alignments`, `aligner_labeled.cpp:361-458`). Consecutive coordinates alone do **not**
prove one record: a column's records are concatenated without separators, so records `ACG` and `CGT` at
coordinates 0 and 1 make the path `ACGT` look consecutive, and the aligner's own test allows exactly that
(`test_aligner_labeled.cpp:320`, `CrossBoundary`). Support is therefore a property of **each label** of a path, not of the path: one label can be
verified while another only carries the constituent k-mers. Each returned label states its `support`:
- `record_verified`: the label's coordinates for the path are consecutive, the first maps to (`seq_id`, local)
  through the `CoordToHeader`, and local + n − 1 < `num_kmers_in_sequence(column, seq_id)`
  (`coord_to_header.hpp`), so the last k-mer is inside the same record. Possible only on a **BASIC** index with
  coordinates **and** record mapping: canonical coordinates cannot establish the deposited orientation (the
  traversal states the same restriction, `walker.cpp:1157`), and coordinates without a `.seqs` file cannot
  establish record bounds.
- `label_intersection`: the label annotates every k-mer of the path, nothing more; the note says why
  (`record_bounds_unknown` with coordinates but no mapping, `label_intersection_only` without coordinates,
  `strand_unknown_canonical` on a canonical index).
- `kmer` (L ≤ k): the label annotates the one k-mer; there is nothing to verify.
A context's `support` summarises its labels (`record_verified` when every returned label is, else
`label_intersection` or `mixed`). The capabilities state the best support an index can give
(`pattern.support`), and a request may demand it: `require_support: record_verified` has two effects. It
refuses the request on an index that cannot verify at all (400 with the index's best), and on an index that
can, it **excludes unverified labels from the returned projection**: a path verified in A but only intersected
in B is returned with A alone, and B is counted in `labels_excluded_unverified` so that the exclusion is
visible. Without the field every label is returned with its own `support`, and the weaker mode is stated.

**Placement (BASIC index, record mapping).** For a context with first-k-mer coordinate c in column col and offset
p: `(seq_id, local) = cth.locate(col, c)` (the global → (seq_id, local) mapping, `coord_to_header.hpp:63-71`), then
`start = local + p + 1` (1-based), `end = start + L − 1`, strand from the orientation; the occurrence carries
`column`, `seq_id` and `record` (the header, for display) separately, since one column can hold many records and
two records of equal length can place a hit at the same local coordinates. The offset is added to the
record-local position, never to the global coordinate: with k = 5 and records `ACGTA`, `CCCCC` at global 0 and 1,
the `TA` at offset 3 of the first record is local 3, not global 3. Without record mapping the answer gives
`kmer_coord` (global, 0-based) and `offset` as two fields and places nothing. The `nt_coords` / `nt_length` format
of `metagraph align --json` is reused for placed occurrences (`alignment.cpp:946-1001`; the writer moves into a
shared header, its caller's output unchanged).

## 5. Admission, budgets and determinism

### 5.1 What "found every occurrence" rests on

- Phase 1 enumerates BOSS ranges, which are exact: a valid edge is in a leaf's `W`-selected set iff its last m
  symbols instantiate Q[0, m); a node is in an offset range iff its suffix is Q·X.
- For L ≤ k, every occurrence of P inside a retained island lies inside at least one retained k-mer, and
  `any_offset` visits every k-mer containing Q (§4.1). No dummy edge is relied on. For `long`, every occurrence
  whose L − k + 1 k-mers were all retained is a path the extension walks from its first k-mer; an occurrence
  with a pruned k-mer is not in the index (§3, "Covered sequence").
- Phase 2 is an exhaustive DFS with no scoring.
- Both strands are searched explicitly.

So, within the admitted work, a label absent from an answer carries no occurrence in the covered sequence
**only when** that answer's `retrieval_complete` is true: every count `exact`, every shard's retrieval finished,
no anchor truncated, no row refused — and the answer was not filtered. With a `predicate`, or with
`output.labels: predicate_only`, the claim shrinks to what was asked: `retrieval_complete` then means that every
**selected** context is returned with its **projected** labels, and it says nothing about the contexts the
predicate excluded or the labels the projection left out; `absence_scope` carries `filtered` beside the scope.
A `count` answer names no labels and claims nothing about any label; a `suffix` answer says
`absence_scope: suffix_only`.

### 5.2 Count first, then retrieve: the modes

The route never spends annotation work on a pattern whose aggregate graph count was not admitted. Two modes,
plus an explicit opt-in for the traversal's "deliver what was built" convention:

| mode | what happens |
|---|---|
| `count` | discover and count the contexts (or anchors and paths) in the requested scope. **Without a predicate** no annotation row is decoded, and completed interval descriptors are **discarded** as they are counted, so the memory account bounds only the DFS frontier and broad patterns can be counted. **With a predicate** (§5.6) the answer carries the raw and the selected counts: descriptors are retained up to `max_predicate_contexts`, the predicate pass reads annotation under its budgets and the backend restrictions of §4.3, and still no result is returned |
| `all_or_count` (default) | the count, with the interval descriptors **retained while retrieval is still possible**: without a predicate until the raw count exceeds `max_contexts`, with a predicate until it exceeds `max_predicate_contexts` (20,000 raw contexts may select 5); past that they are discarded and counting continues within the compute budget; retrieval of **all** requested results only if the completed aggregate count — the selected count when there is a predicate — is within the retrieval threshold and the retrieval finishes on every shard; otherwise the counts alone, with the reason |
| `partial` | the count, then retrieval up to the caps with what was built delivered and every cut stated — an agent that wants the first results of a large set asks for it explicitly |

**Descriptor retention by phase.** For L ≤ k the descriptors are the contexts' intervals, kept or discarded as
the table says. For `long` the phases differ: **anchors** are always kept through the extension admission
(they are what the extension starts from) and discarded once every anchor has been extended or the admission
failed; **completed paths** are discarded at once in an unfiltered `count`, kept up to `max_paths` for ordinary
retrieval, and kept up to `max_predicate_contexts` when a predicate follows. `stop_at_threshold: true` ends the
current phase as soon as its own threshold is crossed and names it. Without a predicate the thresholds are
`anchors` (`max_anchors`), `paths` (`max_paths`) and `contexts` (`max_contexts`). With a predicate the raw
phases stop at `max_predicate_contexts` (`raw_contexts` or `raw_paths`; anchors still at `max_anchors`), and
only the `selected` phase stops at `max_contexts` or `max_paths`: a raw count above the retrieval threshold is
not a reason to stop when a filter follows. Either way the counts are `at_least` and results are withheld
(`withheld: threshold_crossed`). Descriptors live for the request only: there is no token for a later retrieval;
an agent that wants the results of an over-threshold count narrows the pattern, the scope or the strand and
asks again.

**Across shards, barriers per pattern.** Every admission is on the aggregate: every shard counts, the counts are
summed per unit with the weakest relation, and only an aggregate `exact` count within the threshold passes —
two shards with 6,000 contexts each do not pass a threshold of 10,000 separately. Without a predicate there are
two barriers: **retrieval admission** (the raw aggregate against `max_contexts`, or `max_paths` for `long`) and
**publication** (results are published only when **every** shard's retrieval finished; a shard whose retrieval
was refused or stopped withholds the pattern's results on all shards, `withheld` naming the shard and the
reason, with the counts kept). With a predicate there are three: **compute admission** (the raw aggregate
against `max_predicate_contexts`; for `long`, anchors against `max_anchors` first and raw paths against
`max_predicate_contexts`), **selection admission** (the selected aggregate against `max_contexts` or
`max_paths`), then publication. Per-shard work and memory budgets stay what §5.3 says; the barriers add no
shared budget.

| situation | answer |
|---|---|
| aggregate counting completes within the threshold; retrieval finishes on every shard | exact graph counts, all requested results, `retrieval_complete: true` |
| aggregate counting completes above the threshold | exact graph counts, results withheld (`withheld: count_above_threshold`) |
| the threshold is crossed and `stop_at_threshold` is set | `at_least` counts, results withheld (`withheld: threshold_crossed`) |
| range or path discovery reaches `max_steps` or the deadline | `at_least`, `bounds` or `unknown` counts, results withheld (`withheld: discovery_budget`) |
| graph counts complete, but a shard's annotation or output budget is exceeded, an anchor's labels are truncated, or a row is refused | exact graph counts, results withheld with the reason and the shard (`withheld: annotation_budget` with `rows_refused`, `anchor_labels_truncated`, or `output_budget`) |

An agent can then tell "many matches" (count above threshold: refine the pattern, the scope or the strand),
"expensive discovery" (discovery budget: shorten the flank scope, add information) and "expensive annotation"
(a common anchor: raise `max_labels_per_anchor` for a narrower pattern, or ask for the count only) apart, and
act on each.

### 5.3 Budgets, the deadline and loading

**The deadline.** One per request, started when the request is parsed and compared at every phase boundary and
inside every loop (the range DFS, the mask scans and the extension every 4,096 steps, the paced reads at their
run boundaries as in `/resolve` and `/traverse`, serialisation every 4,096 objects). The route keeps a
**finalisation reserve** (`--pattern-finalize-ms`, default 250 ms, stated in the capabilities): the count answer
of every pattern is maintained incrementally during the search, so that when the deadline less the reserve
passes, the route stops all work and serialises the counts it has, with `stop: {phase, shard}`,
`withheld: deadline` and the relations of §3. If serialisation itself overruns the reserve, the answer is the
explicit outcome 503 `deadline`, as a `/traverse` attempt past its bound; nothing partial is sent as if whole.

**Resident indexes only.** The route serves the graphs resident in the process (mmap or RAM, as loaded at
start-up or by an earlier `/search` with `in_ram`); it never loads an index inside a request, because a
synchronous load cannot be interrupted at the deadline. `in_ram` is refused on `/pattern` (400), and on a
multi-graph server a shard whose index is not resident answers at once with `stop: not_resident` for that shard,
counted as a stopped shard in the barriers of §5.2. The capabilities list which graphs are resident.

**Memory and work per shard.** Each shard's search gets a fixed share of `max_memory_mb` (equal shares, or the
request's `budget_split`), its own `DecodeBudget` (single-thread-owned, as the class requires), and its own
`max_steps` share, so what a shard admits does not depend on what another shard admitted in the meantime.

**What the account holds.** A `DecodeBudget` bounds the rows a fetch decodes, not what the reused classes keep:
`LabelQuery` and `LabelRecorder` copy their results into caches that are unbounded unless the caller sets a
limit (`set_max_cache_bytes`, `label_oracle.hpp:583-587`), and the traversal's walker reserves cache space and
prices its dictionaries, label maps and output into one account beside the reads (`walker.cpp:3286-3300`). The
pattern search keeps the same discipline, so that fetches admitted one by one cannot add up past the advertised
limit: the shard's share is one account that holds, over **all patterns of the request**, (1) the oracle caches,
given a fixed allotment of the share through `set_max_cache_bytes` and evicting within it; (2) the label maps and
dictionaries of discovery; (3) the retained descriptors of §5.2; (4) the deduplication state of §5.4; (5) the
contexts and labels built for the answer and the buffered response text; and (6) the `DecodeBudget` for the
reads, charged against what the rest has left. A fetch is admitted against the account's remainder, not on its
own; what any item would push over the share is a stated memory stop.

| cap | default | server flag | what is stated |
|---|---|---|---|
| `min_information_bits` | 24 (≈ 12 specified bases) | `--pattern-min-information-bits` | gates **discovery** for `iupac`, `protein`, `any_offset` and `long` (its anchor window); an exact `dna` `suffix` count is always admitted, however short: one range, a few ranks, and a scan of the range's invalid edges when it has any, charged like any step |
| `max_contexts` | 10,000 per pattern, all offsets and both strands, aggregate over shards | `--pattern-max-contexts` | the retrieval threshold of `all_or_count`; the count is always returned with its relation |
| `max_anchors` (`long`) | 1,000 | `--pattern-max-anchors` | the extension threshold; `withheld: anchors_above_threshold` |
| `max_paths` (`long`) | 1,000 | `--pattern-max-paths` | the retrieval threshold on completed paths; `candidates_examined` |
| `max_steps` | 1,000,000 range, scan or edge steps per shard | `--pattern-max-steps` | `ranges_visited`, `steps`, the phase that hit it |
| `max_labels_per_anchor` | 64 | `--pattern-max-labels-per-anchor` | `anchors_truncated` with the cap and each row's total |
| `max_predicate_contexts` (§5.6) | 100,000 per pattern | `--pattern-max-predicate-contexts` | the compute admission of a predicate pass; `withheld: predicate_above_threshold` with the unfiltered count |
| `max_predicate_work` (§5.6) | the oracle's units, default as `max_annotation_work` | `--pattern-max-predicate-work` | `selection.tested` as `at_least`, `withheld: predicate_budget` |
| `max_predicate_labels` (§5.6) | 10,000 names per request | `--pattern-max-predicate-labels` | refused above it (400) |
| `max_annotation_work` | the oracle's units, default as `/resolve`'s | `--pattern-max-annotation-work` | `rows_refused` with their anchors |
| `max_memory_mb` | 256 per request, split per shard | `--pattern-max-memory-mb` | `memory.stop` with the phase and the shard |
| `max_labels` | 1,000 | `--pattern-max-labels` | `labels` count with its relation; labels kept by (contexts desc, column asc), as `/search`'s top-N; `partial` only (`all_or_count` returns all or none) |
| `max_occurrences_per_label` | 16 | `--pattern-max-occurrences` | `occurrences` count per label with its relation; `partial` only |
| `time_budget_ms` | 5,000 | `--pattern-max-time-ms` | `stop: time` with the phase and the shard; the finalisation reserve inside it |

Why information rather than length, and why only as a gate: expected suffix anchors ≈ (k-mers in the graph) ×
2^−bits. On a 50-billion-k-mer graph a 16-mer expects about 12 per strand, a 12-mer about 3,000, a 10-mer about
48,000. That is a planning heuristic; what bounds a run is the range work and the deadline, which is why
`count` is always available and `max_steps` is the real admission for discovery.

Low-complexity patterns are not refused: they are legitimate questions, the caps bound them, and
`is_low_complexity` (`aligner_seeder_methods.cpp:22`) marks them in `notes` so the reader knows why the counts
are large.

### 5.4 Deduplication

Graph contexts are distinct by (orientation, path, offset). Placed occurrences are distinct by (column, `seq_id`,
start, strand): the up to k − L + 1 contexts of one interior occurrence collapse to one placed occurrence, and
two offsets in one k-mer stay two. Node id is never a deduplication key.

### 5.5 Determinism

Contexts are enumerated in BOSS edge order, then by offset; labels are ordered by (contexts desc, column asc);
occurrences by (`seq_id`, start, strand); shards by (graph name, `index_fp`) (§8). The same request on the same
index gives the same answer, including which results a cap cut, because every budget a shard spends is its own
fixed share. The exception is a deadline stop, which depends on the machine and is stated
(`determinism: time_limited`).

### 5.6 Annotation predicates: select contexts by a condition on their labels (extension, the owner's request)

A broad pattern with thousands of contexts becomes a small, biologically meaningful answer when only the
contexts whose annotation satisfies a condition are returned: "present in sample A or B and absent from C", "in
at least three samples of cohort A and in none of cohort B". The condition is evaluated on the decoded
annotation of each context, which costs annotation work, but far less than returning the whole answer, and it
keeps the count-first contract: the request is admitted on counts, filtered or not.

**Two operations, kept apart.**
- **Selection:** which contexts pass the predicate. Evaluated on the labels of the context's row (L ≤ k) or, for
  `long`, on the **path's shared label set** (the intersection along the path, coordinate-consistent where
  coordinates exist, §4.3), never row by row: one row carrying A and another carrying B both satisfy
  `any(A, B)`, while neither sample supports the whole path. The predicate asks about support of the context,
  so for a path it is asked of the labels that support all of it.
- **Projection:** which labels and coordinates are returned for a selected context. `output.labels:
  predicate_only` (the default) returns only the labels the predicate names, with their coordinates where
  placement exists; `all` runs the ordinary discovery of §4.3 for the selected contexts, under the per-anchor
  cap. A selected row can carry thousands of other labels; selecting it does not mean returning them. A
  selected context is returned **whether or not** its projected list is empty: under `none(A)` a context
  annotated only with B passes and comes back with its k-mer, instance, offset and strand and `labels: []`,
  which is why results are contexts and not columns (§7.2) — the k-mer is the identity and the traversal
  seed, and the agent asks for `all` when it wants to know what the context does carry.

**The predicate language (v1).** A predicate names **annotation columns** only, by their column labels (on a
header-labelled index the headers are the columns), as explicit lists of at most `max_predicate_labels`
(default 10,000) names; named cohorts are a service-side convenience that expands to lists (§9). The records a
`CoordToHeader` maps inside a column (`seq_id` on a taxid-labelled index) are not predicate terms in v1: header
text need not be unique, and the existing header lookup keeps the first duplicate, so a predicate over records
would need `(column, seq_id)` references with a defined expansion, which is a later increment (§12). A name that
is not a column is "unknown" (above).

| predicate | meaning |
|---|---|
| `any(labels)` | at least one listed label is present |
| `all(labels)` | every listed label is present |
| `none(labels)` | no listed label is present |
| `at_least(n, labels)` | at least n of the listed labels are present |
| `and`, `or`, `not` | combine predicates |

A label unknown to a shard is simply absent there (`any` false, `none` true); the answer lists unknown labels
per shard (`predicate.unknown_labels`) so that a typo is not read as absence.

**Evaluation.** The predicate names a label set, which is exactly what `LabelQuery` takes as its permitted set
(§4.3): the rows of the admitted contexts are read through `LabelQuery::fetch` restricted to the predicate's
labels, budget-aware and paced, and the predicate is evaluated on the restricted row. Labels unknown to a
shard are folded away before the query is built (`any` of unknowns is false, `none` true, `at_least(n, …)` with
fewer than n known labels false): `LabelQuery` rejects an empty permitted set, and a predicate whose labels are
all unknown needs no annotation read at all. **Access path:** the budgeted, paced reads decode rows
(`LabelQuery::decode_run`, `label_oracle.cpp:1014`) and keep the named columns; the oracle's selected-column
membership (`Access::DIRECT`) exists only on the unbudgeted path. v1 therefore states `predicate.access: rows`
for budgeted predicates, offers `columns` only through `allow_unbudgeted_annotation: true`, and prices a paced,
budgeted direct access as a later increment (§12). On the row-diff family, asking for fewer labels does not make
a row cheaper to decode. For L ≤ k, coordinates and label-name expansion are deferred until a context has
passed.

**`long`: support before selection.** For a path the predicate is asked of the labels that support the whole
path, and on a coordinate index "support" means consecutive coordinates inside one record (§4.3), which has to
be established **before** the predicate is evaluated, not after: with k = 3, a sample A holding the separate
records `ACG` and `CGT` annotates every node of the graph path `ACGT` and contains no occurrence of it, so a
row-wise `any(A)` would select the path and a row-wise `none(A)` would reject it, both wrongly. The pass runs in
two steps: (1) the membership intersection along the path, restricted to the predicate's labels, as a
preliminary filter (a label absent from any node cannot support the path); (2) for the (path, label) pairs that
survive, the coordinates are read and the consecutive-and-in-bounds check of §4.3 decides support; the
predicate is then evaluated on the **supporting** labels and the selected paths are counted. The support a predicate is evaluated on is the index's best (§4.3): `predicate.support: record_verified`
needs a BASIC index with coordinates and record mapping; on every other index the first step is all there is,
the answer states `predicate.support: label_intersection`, and every selected path's labels carry that support.
Because this changes what `none(A)` means — "no record of A contains the path" against "A does not annotate
every k-mer of it" — a request that needs the strong reading sets `require_support: record_verified`: the
predicate is then evaluated on verified labels only, unverified labels are excluded from the projection and
counted as excluded, and an index that cannot verify refuses the request rather than answering in the weaker
mode.

**Admission and counts.** The pipeline becomes: (1) count the unfiltered contexts (or anchors and paths) as in
§5.2 — the unfiltered aggregate admits the **compute** of the predicate pass against `max_predicate_contexts`
(default 100,000 per pattern; for `long`, anchors against `max_anchors` first, then raw paths): a broad pattern
is not refused for being broad when its filter is selective, and its descriptors are retained up to that cap;
(2) evaluate the predicate on the admitted contexts under `max_predicate_work` (annotation units) and the
deadline; (3) count the contexts that pass (`counts.selected`, with `counts.tested` beside it); (4) the
**retrieval** threshold — `max_contexts`, or `max_paths` for `long` — applies to the selected aggregate, and
`all_or_count` returns all selected results or the counts. A predicate pass stopped early leaves
`counts.selected` as `at_least` over the contexts tested and withholds results (`withheld: predicate_budget`).
The three shard barriers of §5.2 apply.

**What a filter does and does not do.** A predicate reduces what is returned and what is placed; it does not
reduce the BOSS backtracking that discovers the contexts, which runs before it. A `withheld: discovery_budget`
therefore calls for a more informative pattern or a narrower scope, not for a filter; a
`withheld: count_above_threshold` is where a selective predicate helps.

**Scope of the claim: `predicate.scope: shard_context`.** The predicate is evaluated per context (one row, or one
path) and per shard. It asks whether **this sequence context** is supported by the labels it names, not whether
the motif occurs anywhere in a sample: with k = 5 and the pattern `AC`, a sample A holding `TACGG` and a sample
C holding `CACCC` give two different rows, and `any(A) and none(C)` selects A's row although C does carry the
motif, in a different flank. Likewise C's support on another shard does not touch the evaluation on this one.
That is the right semantics for finding sequence variants and contexts; it does **not** answer "motifs present
in A and absent throughout C". That question is an aggregation over all contexts of the pattern — a label
carries the motif if any context carries it — which needs the unfiltered discovery of every context's labels
and a predicate over the pattern-level label sets, across shards; it is the later increment "motif-level
predicates" (§12), and the answer names its scope so that no one reads a context-level result as a motif-level
one. Absence of selected contexts is claimed only with `retrieval_complete: true` and a predicate pass that
tested every admitted context (`counts.tested` equal to the unfiltered count, `exact`). Predicates over
coordinates (strand, position within a record), over label properties the index does not carry, and quantified
over cohorts larger than the list cap are later increments too.

## 6. Peptides: the codon automaton

A residue is a set of codons; the automaton's state is (residue index, position in codon, codon prefix so far)
and `allowed(position, prefix)` returns the bases that extend the prefix to one of the residue's codons. This
admits exactly the residue's codons: no superset for Leu (TTA, TTG, CTN), Ser (TCN, AGY) or Arg (CGN, AGR), no
translation filter, and no branch through a stop codon. The standard code (table 1) is built in; `genetic_code`
is a request field for others.

Both strands: the reverse-complemented automaton — residues in reverse order, each codon replaced by its
reverse complement (GCN becomes NGC, AAR becomes YTT) — is searched as the oriented pattern rc(P); it is a
pattern like any other. A peptide of m residues is a pattern of 3m positions, so with k = 31 peptides up to 10
residues stay in the one-k-mer case; longer ones extend (§4.2), carrying the automaton state across the k
boundary (§4.1), and need `record_verified` support for a record claim.

## 7. Contract: request and answer

### 7.1 Request (`POST /pattern`)

```json
{
  "patterns": [
    {"id": "p1", "dna": "GCGTCAGGTACCGGA"},
    {"id": "p2", "iupac": "TTGACANNNNNNNNNNNNNNNNNTATAAT"},
    {"id": "p3", "protein": "MKVLAAGIVG"}
  ],
  "mode": "all_or_count",
  "scope": "any_offset",
  "strands": "both",
  "stop_at_threshold": false,
  "max_contexts": 10000, "max_anchors": 1000, "max_paths": 1000, "max_steps": 1000000,
  "max_labels_per_anchor": 64, "max_annotation_work": null, "max_memory_mb": 256, "time_budget_ms": 5000,
  "max_labels": 1000, "max_occurrences_per_label": 16,
  "allow_unbudgeted_annotation": false,
  "require_support": null,
  "predicate": {"and": [{"any": ["SRR100", "SRR101"]}, {"not": {"any": ["SRR900"]}}]},
  "max_predicate_contexts": 100000, "max_predicate_work": null,
  "graphs": ["name"],
  "output": {"occurrences": true, "labels": "predicate_only", "paths": false}
}
```

`predicate` is optional (§5.6); without it `output.labels` is `all` and the predicate caps are unused.

Exactly one of `dna`, `iupac`, `protein` per pattern; at most `--pattern-max-patterns` (default 16) per request;
`id` optional. Unknown fields are refused (400), as `/traverse` refuses them; so is `in_ram` (§5.3). The route
name avoids the `/search` prefix on purpose: the server matches `^/search` without an end anchor
(`server.cpp:768`).

### 7.2 Answer

One entry per pattern, in request order. Every count is `{value, relation, unit}`; the example abbreviates
`exact` counts as numbers.

```json
{
  "id": "p1", "kind": "dna", "length": 15, "information_bits": 30, "anchor_information_bits": null,
  "scope": "any_offset", "strands": ["+", "-"], "palindromic": false, "mode": "all_or_count",
  "placement": "record",
  "counts": {
    "contexts": {"value": 412, "relation": "exact", "unit": "graph_contexts",
                 "suffix": 37, "by_offset": {"0": 25, "1": 24, "16": 37}, "by_strand": {"+": 210, "-": 202}},
    "labels": {"value": 12, "relation": "exact", "unit": "labels"},
    "occurrences": {"value": 51, "relation": "exact", "unit": "placed_occurrences"}
  },
  "work": {"ranges_visited": 1203, "mask_scans": 0, "steps": 1203, "annotation_rows": 412, "memory_bytes": 184320},
  "retrieval_complete": true, "withheld": null, "rows_refused": [], "anchors_truncated": [],
  "selection": {"predicate": {"...": "..."}, "access": "rows",
                "tested": {"value": 412, "relation": "exact", "unit": "graph_contexts"},
                "selected": {"value": 9, "relation": "exact", "unit": "graph_contexts"},
                "unknown_labels": {"refseq33m": []}},
  "results": [
    {"graph": "refseq33m", "index_fp": "907ec9…", "release": "refseq97-804731aa…",
     "kmer": "TTGACGCGTCAGGTACCGGAACTGATCCAGT", "instance": "GCGTCAGGTACCGGA", "offset": 4, "strand": "+",
     "support": "kmer",
     "labels": [
       {"column": "562", "support": "kmer", "occurrences": [
          {"seq_id": 17, "record": "NZ_CP012345.1", "nt_coords": "18231-18245", "nt_length": 4912034},
          {"seq_id": 403, "record": "NZ_CP099999.1", "nt_coords": "77-91", "nt_length": 3120000}]},
       {"column": "573", "occurrences": [
          {"seq_id": 2, "record": "NZ_AP018001.1", "nt_coords": "1200456-1200470", "nt_length": 5300000}]}
     ]}
  ],
  "by_label": [
    {"graph": "refseq33m", "column": "562", "contexts": 3, "contexts_suffix": 1,
     "occurrences": {"value": 3, "relation": "exact", "unit": "placed_occurrences"}}
  ],
  "absence_scope": "any_offset",
  "stop": null, "determinism": "full",
  "limits": {"max_contexts": 10000, "...": "..."},
  "timing": {"elapsed_ms": 41.2, "phase1_ms": 3.1, "extension_ms": 0, "discovery_ms": 20.0, "placement_ms": 15.0},
  "notes": []
}
```

- **A result is a graph context, and it always carries its k-mer** (`kmer`, the k bases of the anchor, read with
  `get_node_sequence`; for `long`, the spelled sequence of the whole path, L bases, and `anchor_kmer`), the
  pattern instance it matched (`instance`, the L bases at `offset`, which for an IUPAC pattern or a peptide names
  the variant and the codons), its `offset` and its strand or orientation. The k-mer is the row's identity: on
  indexes without placement it is the only handle the agent gets, and every reported k-mer is a valid seed for
  `/traverse` on its own graph. `labels` lists the context's labels (all of them under the per-anchor cap, or the
  predicate's with `output.labels: predicate_only`), each with its placed occurrences where placement exists;
  `column` and `properties` are written by the same `get_label_as_json` as `/search` (`query.cpp:159`);
  `seq_id` and `record` are per occurrence.
- `by_label` is the per-label summary over the returned contexts (contexts, suffix contexts, placed occurrences
  with their relation), ordered as §5.5 orders labels; it answers "which samples carry it" while `results`
  answers "which sequences match".
- For `long`: `counts.anchors` (unit `anchors`), `counts.paths` (unit `paths`, with `candidates_examined`), and
  `support` on every returned label (`record_verified` | `label_intersection`, §4.3), summarised per context;
  `kmer` on every label of an L ≤ k result.
- `by_label` occurrence totals, and the `max_occurrences_per_label` cap, are over the **deduplicated union** of
  a label's placed occurrences across the returned contexts (§5.4), never the sum of the nested lists: the
  k − L + 1 contexts of one interior occurrence contribute one.
- `placement` is `record` (BASIC, record mapping), `global` (BASIC, coordinates without `.seqs`: occurrences
  carry `kmer_coord` and `offset`), `none` (no coordinates) or `none_canonical`.
- `retrieval_complete` is the one flag that licenses an absence claim (§5.1); `withheld` names why results are
  absent (§5.2) and the shard; `stop` names the phase and shard of a deadline, step or memory stop;
  `determinism` is `full` or `time_limited`.
- `selection` is present only with a `predicate` (§5.6): the contexts tested and selected, each with its
  relation, the access path and the labels unknown per shard; `results` then holds the selected contexts only,
  and with `output.labels: predicate_only` each result carries only the predicate's labels.
- `notes`: `strand_unknown_canonical`, `low_complexity_pattern`, `record_bounds_unknown`,
  `label_intersection_only`, `annotation_unbudgeted`, `graph_cleaned`.
- `output.paths` adds the node path (ids, scoped by `graph`) for `long` results; the k-mer and the spelled
  sequence are always present, so a client needs the ids only to address nodes directly.

A refused pattern (too little information for its scope, bad alphabet, `suffix` on a wrapped PRIMARY graph) gets
an `error` entry in its slot and the request still answers the others, as `/search` answers per sequence; a
malformed request is a 400 for the whole; a request whose finalisation reserve is overrun is a 503 `deadline`.

### 7.3 Capabilities

`GET /capabilities` lists `"pattern"` under `features` and `routes.pattern = "POST /pattern"`
(`server.cpp:1421-1480`), with a `pattern` block: the floors, caps and the finalisation reserve, the modes and
scopes (and which scope each graph mode supports), `placement` (the best this index can give), `strand_stated`,
`mask: present | absent`, `annotation: budgeted | unbudgeted`, the build alphabet and `graph_cleaned`, the
resident graphs, `records_shorter_than_k: not_indexed`, and `pattern_contract_version: 1`. A client gates on the
feature, as the search service gates `traverse` on the feature list (§9).

## 8. Multi-graph servers

`/pattern` fans out over the shards of the named graphs as `/search` does (`server.cpp:790-885`: one task per
(graph, annotation) pair on `graphs_pool`), through a helper factored out of that loop and used by both routes.
Differences from `/search`'s merge, all deliberate:

- every result carries its shard: `graph`, `index_fp` and `release`, so node ids and columns are scoped; a column
  present in two shards is two results, and `labels` counts (shard, column) pairs, stated by `per_shard: true`;
- every shard runs under its own fixed share of the memory, step and annotation budgets (§5.3), with its own
  `DecodeBudget`, so admission is deterministic; the deadline is the request's one instant;
- the barriers of §5.2: the aggregate count admits retrieval; results are published only when every shard's
  retrieval finished;
- the merged list is **sorted** after the fan-out, by (graph name, `index_fp`) and then the order of §5.5:
  `/search` appends in completion order, which this route may not;
- counts are summed per unit over shards with the weakest relation (`exact` only if every shard's is), and listed
  per shard under `by_shard`;
- no shard is loaded inside the request (§5.3).

`/align` stays single-graph (400 on multi-graph servers, `server.cpp:898-902`); nothing here changes that.

## 9. Serving the extension synchronously

Three layers, all request/response; no job, no polling, no S3.

**Server.** `POST /pattern` goes through `process_request` like `/search` (compact JSON, gzip when accepted,
`kContentTimeoutS` 900 s, `server.cpp:52`). The route's own deadline (default 5 s, cap `--pattern-max-time-ms`,
finalisation reserve inside it) keeps it far under that.

**Python client.** `GraphClientJson.pattern(patterns, **options)` and `GraphClient.pattern(...)` beside `align()`
(`api/python/metagraph/client.py:100`, `_json_seq_query` at :117): the same `requests.post` with the route's
JSON and timeout handling; `GraphClient.pattern` returns two DataFrames: `contexts`, one row per returned context (pattern, graph, k-mer,
instance, offset, strand, support, the number of its labels), which keeps a selected context whose projected
label list is empty; and `occurrences`, one row per (pattern, graph, context, column, placed occurrence) with
nullable occurrence fields where nothing is placed. The counts come as attributes with their relations, and
per-label totals are the deduplicated union of §7.2. `MultiGraphClient.pattern` fans out over hosts with the existing thread
pool.

**Search service** (its repository; the owner's call which session builds it). A synchronous MCP tool
`pattern_search(database, dna | iupac | protein, mode, scope, strands, ...)` and REST `POST /pattern/search`,
built like `traverse_resolve` (`app/mcp_server.py:5838`; REST `app/api/routes_traversal.py:766`, admission
`_admit` at :502):

- admission by a per-host semaphore with its own settings (`PATTERN_MAX_CONCURRENCY`, `PATTERN_PER_HOST`,
  `PATTERN_TIMEOUT_S`, `PATTERN_MAX_PATTERNS`), the pattern the service already uses for its synchronous tools
  (`_COMPARE_SLOTS`, `app/mcp_server.py:5516`; `TRAVERSAL_RESOLVE_*`, `app/settings.py:1292-1326`);
- the backend called from the API process on the thread pool (`anyio.to_thread.run_sync`), one host per
  database, timeout = the server's `max_time_ms` plus an allowance, as resolve does;
- served only for databases whose host lists `pattern` in `/capabilities`; others answer
  `pattern_unsupported` with the host's feature list, as `traversal_disabled` does;
- the answer passed through with the service's additions only: `database`, `hit_unit`, the untrusted-data notice
  on label and record strings, the standard error envelope; the tool's docstring explains the `withheld`
  reasons, `retrieval_complete` and `absence_scope`, and what to change for each;
- rate limits and the anonymous-call budget as for the other synchronous tools; `count` mode is the cheap default
  the tool suggests for a first look at a new pattern.

Archive-scale databases whose hosts serve several graphs are reachable through the same tool once their hosts
build the feature; the fan-out happens inside the server (§8). A batch of many patterns, or a scan over every
database, is the job queue's shape and comes later, wrapping this route.

## 10. Code placement and reuse

| piece | where | reused from | change to the source |
|---|---|---|---|
| range DFS with `allowed(depth, prefix)` | `src/graph/alignment/pattern_search.{hpp,cpp}`, `PatternSearch` | `suffix_to_prefix` (`aligner_seeder_methods.cpp:95-139`) | the template moves to `aligner_seeder_methods.hpp` and gains the symbol-set callback whose default is today's loop (every non-sentinel symbol), so the seeder's calls compile to today's behaviour on every alphabet |
| leaf counting and enumeration (`W` plain and marked, the `valid_edges_` correction with its charged scan, range normalisation) | `DBGSuccinct::count_edges_with_suffix_pattern(...)` / `call_...` beside the existing function | `call_nodes_with_suffix_matching_longest_prefix` (`dbg_succinct.cpp:308-390`), `BOSS::index_range`, `tighten_range`, `rank_W`, `pick_edge`, `pred_last`/`succ_last`, `valid_edges_` rank/select | the existing function is untouched; the new ones answer its TODO at :350 (count plus first N) |
| reverse-complement anchors on wrapped PRIMARY graphs (`long`) | `PatternSearch` | the full-k-mer mapping of `SuffixSeeder` (`aligner_seeder_methods.cpp:303-309`), `get_base_dbg_succ` | used as functions, for full k-mers only; the seeder is untouched |
| extension beyond k | `PatternSearch::extend` | none (§4.2) | new, ~150 lines |
| label discovery and coordinate retrieval with budgets and pacing | `PatternSearch::retrieve` | `LabelRecorder::fetch` then `LabelQuery::fetch` with `DecodeBudget` and `ReadPacing` (`label_oracle.hpp:491-566, 698-727`), `LabelOracle::decode_charged` | used as they are; a one-pass paced visitor is §12 |
| record mapping and bounds | `PatternSearch::place` | `CoordToHeader` (`coord_to_header.hpp`: global → (seq_id, local), `num_kmers_in_sequence`) | used as is |
| JSON | `src/cli/pattern.cpp` | the `nt_coords` writer of `Alignment::to_json` (`alignment.cpp:946-1001`), `get_label_as_json` (`query.cpp:159`) | the writers become shared helpers in a header; their callers' outputs do not change |
| codon automaton, genetic code | `src/common/seq_tools/genetic_code.{hpp,cpp}` | none in the tree | new, with the standard table |
| route, fan-out, barriers, capabilities | `src/cli/server.cpp`, `server_utils.cpp` | the `/search` fan-out loop (`server.cpp:790-885`) | factored into a helper both routes call; the merge sorted, budget-shared and barriered for this route |
| CLI | `metagraph pattern` (`src/cli/pattern.cpp`, `config.cpp`) | the `align --json` output style | new subcommand, for tests and offline use |
| Python client | `api/python/metagraph/client.py` | `align()`, `_json_seq_query` | one method per client class |

## 11. Not breaking the aligner

- No scoring, seeding or extension default changes. `initialize_aligner_config` (`align.cpp:33-69`),
  `DBGAlignerConfig`, `SuffixSeeder`, `DefaultColumnExtender`, `LabeledAligner`, `AnnotationBuffer` keep their
  code paths; the only edits inside alignment files are the template signature with a default that reproduces
  today's symbol loop on every alphabet (§10) and the move of two JSON writers into a shared header.
- **Golden gate for alignment:** before the first change, record the outputs of every
  `integration_tests/test_align.py` case (8 tests × graph representations, DNA5 included), of `test_query.py`'s
  `align` cases and of `tests/annotation/test_aligner_labeled.cpp` as fixtures; the gate compares byte for byte
  after every increment.
- The seeder's seed lists on the test graphs are captured once (a small gtest dumping `SuffixSeeder` seeds for
  the mini graphs, a DNA5 graph among them) and compared after the `allowed` change.
- `/align` and `/search` with `align: true` answer byte-identically on the fixture requests (the byte-identity
  harness of the traversal work, `scripts/traversal/bench_traverse.py --compare`, with an alignment request set).

## 12. Later increments on the same engine

- **A one-pass paced row/tuple visitor** replacing the discovery-then-placement pair (§4.3), extracted from the
  two fetch loops.
- **Motif-level predicates** (§5.6): quantifiers over all contexts of a pattern (`any_context`, `no_context`),
  across shards, on the pattern-level label sets, for "present in A and absent throughout C".
- **Paced, budgeted selected-column access** for predicates (`predicate.access: columns` under a budget), so that
  `any(A) and none(B)` on a column-major backend reads two bits instead of a row.
- **`suffix` scope on wrapped PRIMARY graphs** through a prefix lookup with orientation conversion.
- **Bounded mismatches** (`NOTE-mismatch-search.md`): automaton states (position, mismatches used) in the DFS;
  the admission policy and caps carry over.
- **Tolerant best-hit search on `/align`:** IUPAC rows in the extender's score table
  (`aligner_extender_methods.cpp:38-59`) plus `min_seed_length` per request and the labeled aligner on the server
  route. Best hits, not every carrier; a separate item.
- **Internal anchors for `long`:** anchor on the most informative k-window and verify both ways.
- **Placement on canonical indexes** for `long` from coordinate steps, shared with `/search`'s `query_coords`.
- **Queue path** in the search service for batches and archive-wide scans.

## 13. Implementation increments (each with tests, then an adversarial review)

0. **The oracle suite first.** A gtest suite of tiny graphs built from explicit records, with **two** brute-force
   oracles, because no single one validates both counting units: a **graph-walk oracle** that enumerates the
   built graph's k-mers and paths directly (every k-mer containing the pattern at each offset, every path
   spelling a long pattern), giving the expected raw graph contexts; and a **record-scan oracle** that scans the
   records' retained islands as strings, both strands, with the pattern as a regex or the codon sets, giving the
   expected placed occurrences and verified supports. At k = 3 the records `ACG` and `CGT` make the graph path
   `ACGT` and no record occurrence: the first oracle expects one context, the second none, and the engine must
   match both. Each oracle lists its expectations with the exact count per unit. It starts
   with the reviews' counterexamples: a record start whose dummy chain is pruned (`TACGAT` then `ACGA`, `AC`);
   the length-k lookup (`AAC` at k = 3); a motif with one plain and two marked matching edges; a sink dummy `AC$`
   inside a range and a source dummy leaving one; an interrupted mask scan with its bounds; a graph loaded
   without its mask (refused); the DNA5 flank (`ACNTA`, `AC` at offset 0); a DNA4 island shorter than k
   (`AAAAANACNCCCCC`, `AC`: not covered, not claimed); zero suffix anchors with a prefix context (`ACGTT`, `AC`);
   the cross-record path (`ACG`, `CGT` at coordinates 0, 1); the double occurrence in one k-mer (`ACGAC`, `AC`);
   the record-local offset (`ACGTA`, `CCCCC`, `TA`); a repeated k-mer within one record (one context, several
   occurrences); two equal-length records with hits at equal local coordinates (`seq_id` tells them apart); a
   palindrome (`ACGT`) counted once; a `suffix` query on a wrapped PRIMARY graph (`ACGA` stored, `CGT` asked:
   refused for `suffix`, found by `any_offset`); a reverse hit of a pattern longer than k on a wrapped PRIMARY
   graph (`AACG` at k = 3, anchored through rc(`AAC`)); an anchor without a completed path (`AAAC` with `AAA`
   present); two shards of 6,000 contexts against a threshold of 10,000 (withheld on the aggregate); a shard
   whose retrieval fails (withheld everywhere); an anchor with more labels than `max_labels_per_anchor`; a native
   CANONICAL graph (presence only); an unbudgeted backend (count only); a refused row; a deadline inside
   serialisation (counts answered within the reserve) and one past it (503); a palindromic anchor window on a
   wrapped PRIMARY graph (`ACGTA` at k = 4: one anchor, one extension) and an IUPAC window with some palindromic
   instances; a pruned middle k-mer (`ACG`, `GTA` kept, `CGT` pruned at k = 3: `ACGTA` not found, both bases
   covered); a path verified in one label and only intersected in another (per-label support); a request with
   `require_support: record_verified` on an index without record mapping (refused); a sequence of fetches whose
   caches would exceed the share (one memory stop, not a silent overrun). Every later increment runs against
   it; the alignment golden gate of §11 is recorded in this step too.
1. **Count-only engine and route** (`mode: count`, `scope: suffix`, exact DNA, one strand, BASIC): mask handling
   settled first (the requirement, the charged scan, the bounds), then `PatternSearch` with the range DFS and the
   `W` rule with marked edges, the exact counts, `POST /pattern` on a single-graph server answering counts, the
   finalisation reserve, capabilities. The early milestone: a server that counts correctly before any annotation
   is read.
2. **IUPAC, both strands, `any_offset`, graph modes.** The `allowed` callback; the two oriented searches with the
   palindrome rule; the flank ranges with every real symbol and `ranges_visited`; native CANONICAL; wrapped
   PRIMARY with `suffix` refused and `any_offset` complete.
3. **Retrieval: discovery, placement, JSON, admission.** `LabelRecorder` discovery with the per-anchor cap, then
   `LabelQuery` placement, both with `DecodeBudget` and pacing; backend support and the unbudgeted opt-in;
   descriptors retained only while eligible; placement through `CoordToHeader` with `seq_id`; deduplication; the
   modes, `retrieval_complete` and the `withheld` reasons; the CLI subcommand. Integration test
   `integration_tests/test_pattern.py` on the mini index: placed occurrences equal a regex scan over
   `build/mini_refseq`'s FASTA (record, 1-based position, strand), counts per unit exact.
4. **`long`.** The extension DFS, anchor-window bits, the two thresholds, per-label `support` with record bounds and `require_support`, the
   anchor mapping on wrapped PRIMARY graphs, the peptide state across the k boundary; tests with patterns
   spanning two and three k-mers, the cross-record path the bounds check must reject, the anchor without a path,
   and the reverse hit.
5. **Peptides.** The codon automaton and genetic code; tests with Leu/Ser/Arg-rich peptides and a peptide whose
   only graph instance crosses a stop codon.
5b. **Predicates** (§5.6). The predicate parser and evaluator over `LabelQuery` reads restricted to the
   predicate's labels; `output.labels: predicate_only`; the compute admission on the unfiltered count and the
   retrieval threshold on the selected count; path-level evaluation for `long`; unknown labels per shard. Tests:
   the oracle suite extended with per-row and per-path cases (one row with A and another with B on a path:
   `any(A, B)` selects no path), `at_least` on cohorts, a typo in a label name reported as unknown, a predicate
   pass stopped by its budget (`at_least`, withheld), and a selective filter on a pattern whose unfiltered count
   exceeds `max_contexts`.
6. **Multi-graph, barriers, per-shard budgets, benchmark.** The shared fan-out helper with the sorted merge,
   shard identity, per-shard budget shares and the barriers; resident-only shards; a benchmark on
   refseq33m-experimental (16-, 20-, 25-nt motifs in both scopes, a 29-nt IUPAC promoter, three peptides)
   recording `ranges_visited`, mask scans, counts per unit, rows, elapsed and memory, and an estimate for the
   Logan hosts from their k-mer counts.
7. **Clients and service.** The Python client methods; the search service's synchronous tool and REST route (its
   repository); docs (`docs/source`); the contract of §7 promoted into a SPEC section once frozen.

Estimate: increments 0–2 about 2 weeks for one engineer who knows the BOSS and traversal code (the count-only
milestone, mask handling included); 3 about 2 weeks (the two-step retrieval is integration work, not reuse);
4 about a week; 5 about a week; 5b about a week (the parser is small; the two-step support check for `long` and
the three admissions are the work); 6–7 about a week and a half, of which the service side is half. The known unknown
is the benchmark in 6: a repeat-rich motif in `any_offset` scope on a Logan host sets the real defaults for
`max_steps`, `max_contexts` and the information floor.

## 14. Open decisions for the owner

1. Route and tool names: `/pattern` and `pattern_search` (this draft), or `/motif`.
2. The default mode (`all_or_count`, this draft) and default scope (`any_offset` for L ≤ k, this draft; `suffix`
   is the cheaper alternative and always states `absence_scope`).
3. Defaults: the information floor (24 bits, discovery only), `max_contexts` (10,000), `max_anchors` and
   `max_paths` (1,000), `max_labels_per_anchor` (64), `max_steps` (10⁶ per shard), `max_memory_mb` (256 per
   request), the finalisation reserve (250 ms); all server policy.
4. Whether `label_intersection` labels are returned for `long` without `require_support`, or only counted. This
   draft returns them, each marked with its `support`.
5. Whether canonical indexes get the route at all in v1 (presence and contexts only), or answer
   `pattern_unsupported` until placement exists.
6. Whether `allow_unbudgeted_annotation` exists at all, or unbudgeted backends get `count` only.
7. The CLI subcommand: keep (tests and offline use) or server-only.
8. Which session builds the search-service side (§9), and whether it waits for the queue path.

## 15. Known limitations, stated in the answer

- Completeness is over retained, k-mer-covered regions: islands shorter than k, unsupported symbols and pruned
  k-mers are outside it; records shorter than k are not in the index at all.
- No placement and no strand on canonical and primary indexes (`placement: none_canonical`); no `suffix` scope on
  wrapped PRIMARY graphs.
- Graph contexts are not occurrences: repeated copies of a k-mer share one context; exact occurrence counts need
  record mapping and admitted retrieval.
- A `suffix` answer covers only occurrences starting at record position ≥ k − L (`absence_scope: suffix_only`).
- Absence of a label is claimed only with `retrieval_complete: true`; a `count` answer claims nothing about labels.
- A label's path for `long` without record mapping is not a proof of one record (`support: label_intersection`).
- A pattern's N never matches a record's N symbol; flanks do.
- The engine needs the succinct graph representation with its edge mask; other representations, or a graph
  without its mask, answer 400 with the reason.
- Budgeted retrieval exists only on the row-diff family with the budgeted decode; elsewhere `count` only, or
  unbudgeted by explicit opt-in; an anchor with more labels than the per-anchor cap withholds the results.
- Indexes are served resident only; a shard not resident is a stopped shard.
- Caps and budgets are per shard on multi-graph servers (`per_shard`); a column in two shards is two results;
  admission is on the aggregate and publication needs every shard.
- A deadline stop is not deterministic (`determinism: time_limited`); a finalisation overrun is a 503.
- Support is per label: `record_verified` only on a BASIC index with coordinates and record mapping; elsewhere
  `label_intersection`, and `none(A)` means only that A does not annotate every k-mer of the path.
- Predicates name annotation columns only; records inside a column are not addressable in v1.

## 16. Changes in v2 (the first review's findings)

| # | finding | change |
|---|---|---|
| 1 | dummy edges do not preserve every record start (pruned at build) | §4.1: the dummy pass is gone; every k-mer containing the pattern at any offset is enumerated from the node ranges |
| 2 | the range loop consumed one symbol too many; the count formula was wrong | §4.1: positions 0 … m − 2 on node ranges, the last on `W`; the existing k − 1 / `pick_edge` distinction kept |
| 3 | consecutive coordinates do not prove one record | §4.3: `chain_verified` needs record mapping, consecutive coordinates and the record's bounds |
| 4 | `AnnotationBuffer` has no budgets; the 5 s claim was unsupported | §4.3, §5: reads through the traversal's budgeted access; one deadline with a defined scope; refused rows stated |
| 5 | canonical placement is ambiguous in offset, not only strand | §3: v1 places only on BASIC indexes |
| 6 | offsets were added to global coordinates | §4.3: map to (record, local) first, add the offset after |
| 7 | node-only deduplication loses distinct occurrences | §3, §5.4: graph-context identity vs placed-occurrence identity |
| 8 | reverse-complement anchor conversion lacked a direction invariant | §4.1: each oriented pattern searched and extended in its own direction |
| 9 | multi-shard results lacked identity and a defined order | §8: shard identity per result; sorted merge |
| — | superset back-translation spends the search budget on false candidates | §6: the codon automaton in v1 |
| — | difficult semantics were priced as reuse | §13: the oracle suite first; estimate raised |

## 17. Changes in v3 (the second review's findings and policy)

| # | finding | change |
|---|---|---|
| 1 | `rank_W(c)` misses marked edges (`c + alph_size`); raw interval sizes include dummies | §4.1: plain + marked ranks, the `valid_edges_` correction, range normalisation |
| 2 | the canonical anchor mapped rc(Q)[0, k), which is Q's last k-mer | §4.1: search Q[0, k) and rc(Q[0, k)); native CANONICAL vs wrapped PRIMARY distinguished |
| 3 | "all four" in the flank excludes DNA5's N and would have changed the seeder | §3, §4.1, §10: flanks and the template's default use every non-sentinel symbol, as today |
| 4 | the (k − L + 1)× work bound is false; flank discovery branches even for exact DNA | §4.1: the bound is withdrawn; `ranges_visited` measured against `max_steps` |
| 5 | counts conflated units and never stated exactness | §3, §7.2: every count is `{value, relation, unit}` |
| 6 | the raw oracle methods are budget-aware but not paced; backends vary; `in_ram` waits untimed | §4.3: the paced classes; supported backends stated; `count` on every backend |
| 7 | occurrences omitted the record identity | §4.3, §7.2: record and column separately per occurrence |
| 8 | shared budgets and parallel shards break deterministic membership | §5.3, §5.5, §8: fixed per-shard budget shares and per-shard `DecodeBudget` |
| policy | partial results after caps hide what the agent needs to know | §5.2: count first; `count` and `all_or_count` modes, `partial` by opt-in; `withheld` reasons; `scope` with `absence_scope`; the information floor gates discovery only |

## 18. Changes in v4 (the third review's findings and corrections)

| # | finding | change |
|---|---|---|
| 1 | `LabelQuery` needs a permitted label set and cannot discover labels | §4.3: two steps, `LabelRecorder` discovery with a per-anchor cap (stated, withholding when cut) then `LabelQuery` placement; a one-pass visitor deferred to §12; estimate raised |
| 2 | `all_or_count` undefined across shards; exact counts alone do not justify absence | §5.2, §8: two barriers per pattern (aggregate admission, publication only after every shard's retrieval); `retrieval_complete` is the one flag that licenses absence; `count` answers claim nothing about labels |
| 3 | `suffix` on wrapped PRIMARY graphs misses matches (a virtual suffix is a stored prefix) | §3, §4.1: `suffix` refused on wrapped PRIMARY (`scope_unsupported`), `any_offset` shown complete there; a prefix lookup deferred to §12 |
| 4 | the mask may be absent; the invalid-edge scan is work; an interrupted subtraction is an upper bound | §4: the mask is required (`mask_required`); §4.1: the scan charged to `max_steps` and deadline-checked; `bounds` relation with `lower` and `upper`; §5.3 wording corrected |
| 5 | completeness claimed over whole records, not over k-mer-covered sequence | §3 "Covered sequence", §5.1, §15: every claim is over retained, k-mer-covered regions; the DNA4 island example in the oracle suite |
| 6 | count-only mode should not retain descriptors; "later retrieval" undefined | §5.2: `count` discards descriptors; `all_or_count` retains them only while eligible; descriptors live for the request only |
| 7 | loading cannot be interrupted; serialisation checks do not produce the count-only fallback | §5.3: resident indexes only, `in_ram` refused, non-resident shards stopped; a finalisation reserve with the counts maintained incrementally; 503 `deadline` as the explicit outcome |
| — | header text is not a stable identity | §3, §4.3, §5.4, §7.2: `seq_id` in the occurrence identity, header for display |
| — | `long` admission on anchors or on paths | §4.2, §5.3: `max_anchors` admits the extension, `max_paths` admits retrieval |
| — | palindromes and strand-specific counts | §3, §7.2: searched once, `strand: "="`, counted once, reported under `both` |
| — | the rc of a codon automaton | §4.1, §6: residues reversed and each codon reverse-complemented |
| — | the owner: a result is a row, and a row needs its k-mer | §3, §7.2: results are context-level and always carry `kmer`, `instance`, `offset` and strand (the spelled sequence for `long`); `by_label` is the per-label summary; `output.kmers` removed as an option |
| ext | the owner's extension request: filter rows by a logical condition on the annotation | §5.6: annotation predicates (`any`, `all`, `none`, `at_least`, `and`/`or`/`not`) evaluated through `LabelQuery` restricted to the predicate's labels; selection kept apart from projection; path-level evaluation for `long`; the compute admission on the unfiltered count, the retrieval threshold on the selected count; §7, §5.3 and increment 5b |

## 19. Changes in v5 (the fourth review: the predicate extension)

| # | finding | change |
|---|---|---|
| 1 | a negative-only filter selects contexts with no projected label, and the column-organised answer had nowhere to put them; absence was over-claimed under filtering | §5.6, §7.2: results are contexts with their k-mer and are returned with `labels: []`; §5.1: with a predicate or `predicate_only`, `retrieval_complete` covers selected contexts and projected labels only, `absence_scope: filtered` |
| 2 | path-level selection deferred the coordinates it depends on | §5.6: support before selection — membership intersection as a preliminary filter, then the consecutive-and-in-bounds check, then the predicate on the supporting labels; `predicate.support: label_intersection` stated where no coordinates exist |
| 3 | `count` promised no annotation; descriptors were discarded at `max_contexts`; admissions and barriers were undefined with a predicate | §5.2, §5.6: `count` with and without a predicate defined; descriptors retained up to `max_predicate_contexts`; compute, selection and publication barriers; the `long` rule on anchors, raw paths and selected paths |
| 4 | a range's upper bound is no bound on the whole query; the lower bound could go negative | §3, §4.1: global `at_least` unless every unexplored portion has a valid upper bound; `lower = max(0, R − t − (I − s))`, `upper = R − t` |
| 5 | context absence is not motif absence; evaluation is shard-local | §5.6: `predicate.scope: shard_context` stated, the flank example, what the semantics answer and what they do not; motif-level predicates deferred to §12 |
| 6 | selected-column access is not on the budgeted path | §5.6: budgeted predicates decode rows (`access: rows`); columns only unbudgeted by opt-in; a paced budgeted direct access deferred to §12; fewer labels do not make a row-diff row cheaper |
| — | all-unknown predicates would build an empty `LabelQuery` | §5.6: constant folding before the query is built |
| — | suffix coverage relative to retained boundaries | §3: positions counted within the retained island, not the record |
| — | 5b missing from the estimate | §13: about a week |

## 20. Changes in v6 (the fifth review)

| # | finding | change |
|---|---|---|
| 1 | verified path support needs BASIC orientation, coordinates and record mapping; support belongs on each label; the weaker mode changes what `none(A)` means | §4.3, §5.6, §7: `support` per label (`record_verified` / `label_intersection` / `kmer`) with a context summary; `pattern.support` in the capabilities; `require_support: record_verified` refuses where unavailable; `chain_verified` retired |
| 2 | a `DecodeBudget` does not bound the reused classes' caches, maps, descriptors or output | §5.3: one account per shard over all patterns — caches allotted through `set_max_cache_bytes`, label maps, descriptors, deduplication state, buffered response, the reads against the remainder — on the walker's model |
| 3 | the two PRIMARY probes double-count palindromic anchors | §4.1: anchors united by wrapper node id before counting, admission and extension; oracle cases for a palindromic and a partly palindromic window |
| 4 | descriptor retention undefined per phase for `long` | §5.2: anchors kept through extension admission; paths discarded, kept to `max_paths` or to `max_predicate_contexts` by mode; `stop_at_threshold` names the phase's threshold |
| 5 | the DataFrame dropped contexts with empty projections; occurrence totals were sums of nested lists | §9, §7.2: `contexts` and `occurrences` tables with nullable fields; totals and caps over the deduplicated union |
| 6 | coverage was defined by bases, not by consecutive retained k-mers | §3, §5.1: retained islands are maximal runs of consecutive retained k-mer starts; `long` needs every constituent k-mer; the "inside one k-mer" argument limited to L ≤ k |
| — | whether predicate names are columns or mapped headers | §5.6: columns only in v1; `(column, seq_id)` references deferred |

## 21. Readiness review (2026-10-07): three clarifications, then implementation

The sixth review found no blocker and asked for three clarifications, made here: with a predicate,
`stop_at_threshold` stops the raw phases at `max_predicate_contexts` and only the selected phase at
`max_contexts` / `max_paths` (§5.2); the oracle suite uses two oracles, a graph-walk enumerator for raw contexts
and a record scan for verified occurrences (§13, increment 0); `require_support: record_verified` excludes
unverified labels from the projection and counts them, besides refusing indexes that cannot verify (§4.3, §5.6).
Stale `chain_verified` wording was replaced by per-label `support`. Implementation starts with increments 0–2
(the oracle suite and the BASIC exact-suffix count milestone, then ambiguity, offsets and graph modes), with an
early real-index benchmark before the annotation and service layers.
