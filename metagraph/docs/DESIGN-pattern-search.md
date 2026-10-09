# Design: pattern search — motifs shorter than k, IUPAC patterns and peptides, counted first, found completely

**Status:** v6 (2026-10-07), **approved for phased implementation** by the sixth external review (readiness,
§21), after five design reviews (v1: nine findings, §16; v2: eight and a resource policy, §17; v3: seven and
four contract corrections, §18; v4: six findings and three corrections on the predicate extension, §19; v5: six
findings and one naming decision, §20; all confirmed against the code). Written from a read-only check of the code at
`804731aa` (gr/labeled-traversal); its line numbers point there, as pointers, not anchors. Implemented since on
that branch: milestone 1 (increments 0–2), milestone 1b (§13), increments 3, 4 and 5, and the owner's decisions of
2026-10-08 (§4.4, SPEC §18): graphs without the dummy-edge mask served with upper bounds and an estimate, the mask
and the Bloom filter outside the index identity, the stop `*` in peptides. `SPEC-pattern-search.md` is the
contract as built; where this design and the build differ, its §13 says how and why, and the review of 2026-10-07
corrected the sentences below that the build did not keep.

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
7. **Served synchronously by the server, as a job by the service:** the server route is a plain
   request/response like `/search`; the search service serves it as a job (submit → status → results), the way
   it wraps `/search` (§9, the owner's decision of 2026-10-07). No Python client method until someone needs one.
8. **Annotation predicates** (§5.6): a broad pattern narrowed by a logical condition on the labels of each
   context, evaluated before anything is materialised, with the retrieval threshold applied to what passes.
9. **A label-free path** (§4.3, `output.labels: none`): counting **and** extracting the matching k-mers, rows
   and paths without reading a single annotation row. Labels are what slows a request down; the k-mers alone
   answer most follow-up questions and seed `/traverse`, so the cheap path is a first-class mode, on every
   backend and every graph mode.

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
| `protein` | a string over the 20 amino acids, the ambiguity codes X (any residue), B (D/N), Z (E/Q), J (I/L) (the owner's decision #15 of 2026-10-08) and the stop `*` (a stop codon of the genetic code; decision #19) | 3 per residue; the codon automaton of the residue (§6) |

A pattern's **information** is Σ log2(4 / |set_i|) bits over its positions (for a peptide, log2(64 / codons) per
residue). For a pattern longer than k the information of each searched orientation's **anchor window** (positions
[0, k) of P, and of rc(P) when the reverse orientation is searched) is stated separately, and the least of them
gates (§5.3): a highly informative suffix does not make an `N…N` anchor cheap. Information is not cost, for any
L: a run of N inside a searched window costs about min(4^run, edges / 4^a) ranges per level, a being the
specified bases before it in that window, whatever the bits; a run at the start of a window costs nothing on a
`$ACGT` graph, where the engine skips it (§5.3).

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
makes the counts exact:** without it `in_graph` treats every edge as a k-mer (`dbg_succinct.cpp:935`), dummies
included. The mask lives in a separate `.edgemask` file beside the `.dbg` (`serialize`, `dbg_succinct.cpp:898`;
loaded if present at :802), and **neither the mini index nor refseq33m-experimental's graph has one today**
(`build/mini_refseq/graph_k31.dbg` and staging's `graph_k31.indexed.dbg` sit without it). Until the owner's
decision #16 of 2026-10-08 the route required the mask and refused a graph without it (`mask_required`); since,
some graphs have masks and some not, and a graph without one is served with upper bounds and an estimate (§4.4).
Three states, stated in the capabilities as `pattern.mask: file | built_at_load | absent` and `pattern.counting:
exact | upper_bound`:
- `file`: the `.edgemask` written by `metagraph build`, or once, offline, by a new `metagraph transform
  --mask-dummy` that calls the build's own `mask_dummy_kmers(threads, false)` (no pruning: node ids and the
  annotation stay valid) and writes the file. The cost (review of 2026-10-07, D1-02, M1-01, X-EFFICIENCY-06):
  - **memory**: an uncompressed bit vector of edges + 1 bits beside the loaded graph, alive until its compressed
    copy is built: about 78 GB for refseq33m's 6.27 × 10^11 edges;
  - **time**: the dummy-tree traversal (`BOSS::mark_source_dummy_edges` → `traverse_dummy_edges`, proportional
    to the dummy edges, at most k − 1 per record, the one part `-p` speeds up), the sink pass (O(sinks), a select
    over W = $; before the review it read every edge's W on one thread: hours on refseq33m by extrapolation), and word-level
    passes over the vector (its flip and its compression). On 2 × 10^7 edges the marking took 2 ms (0.19 s on a
    STAT and 5.0 s on a SMALL graph before the sink pass changed); at refseq33m's size it has not been run, and
    the 78 GB, loading the graph and the word-level passes are what it costs.
  - **where**: on the host, not inside the 128 GiB container, with `--mmap` so that the graph's pages are mapped
    rather than read into RAM.

  This is the one-time step for staging's exact counts, run by the owner on mex (decision #18: in a staging-only
  directory, once decisions #16 and #17 are deployed; the route answers there before it, with upper bounds). What
  it changes beyond the route:
  - every loader of the graph reads the `.edgemask` from then on, the build already deployed included: node ids,
    rows and what `/search`, `/align`, `/resolve` and `/traverse` find stay as they were, but `GET /stats`
    `graph.nodes` becomes the k-mer count, and a `.bloom` beside the graph is loaded from then on;
  - a server started with `--index-manifest` keeps its manifest and its `index_fp` (the owner's decision #17 of
    2026-10-08: the mask is derived data of the graph, not part of the identity; before it, the manifest had to
    be regenerated and `index_fp` changed). The order is: `transform --mask-dummy`, then `update.sh` (it restarts
    the container even without a new commit). No manifest is regenerated and nothing is re-hashed;
  - the mask is read when the graph is loaded: a server running when the file is written counts upper bounds
    (`counting: upper_bound`, §4.4) until it is restarted;
  - the mask is checked once when the graph is loaded for the one known defect (the owner's decision of
    2026-10-07): `metagraph extend` on a masked graph used to write a mask that marks its new dummy edges valid
    (`DBGSuccinct::add_sequence`, its TODO), and counts on such a graph could be overstated, `exact` included. A
    mask that marks a W = `$` edge valid makes the route answer `mask_invalid` (capabilities `available: false`)
    until the graph is masked again with `transform --mask-dummy --force` and the server restarted; `extend`
    now rebuilds the mask of a masked graph (§15). Beyond that check the mask is trusted as written.
- `built_at_load`: a server started with `--pattern-build-mask` builds the same mask in memory when the file is
  absent, with `--threads-each` threads on `server_query` (`-p` in the `pattern` CLI), at every start-up and
  before any route answers (the same transient vector, so for small indexes and the tests — "small" meaning few
  edges, not the SMALL state —, not for a 362 GB graph on every restart). No file changes: `/stats` changes, and
  `mask: built_at_load` tells it from a file (`index_fp` stays either way, decision #17).
- `absent`: served since the owner's decision #16 (`counting: upper_bound`, §4.4); the start-up log names the two
  ways to exact counts (`transform --mask-dummy` once, read at the next start; `--pattern-build-mask`). Before the
  decision the route answered `mask_required` here.
The test setup builds the mini index's mask with the `transform` step.

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
checked against the deadline like every other loop; it is real work, not a rank. With R = plain + marked, J the
range's invalid edges whose `W` is not the sentinel (a sink dummy never carries c), s of them scanned and t of
those matching, an interrupted scan leaves the range with `lower = max(0, R − t − (J − s))` (every unscanned such
edge might match) and `upper = R − t`; as built, the engine scans the candidates instead when they are fewer than
the invalid edges (either resolves the count exactly), and a range at an offset where a palindromic k-mer is
possible on an even-k wrapped PRIMARY graph is scanned for its palindromes (SPEC §7.6). The count of
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
  unchanged: on a DNA5 build the flank admits N, so `ACNTA` at k = 5 is meant to be found for `AC` at offset 0 —
  design intent: no DNA5 build has run the route, §15); **all** valid
  edges leaving those nodes are such k-mers, counted by `valid.rank1(last) − valid.rank1(first − 1)`, which also
  excludes sink dummies;
- offset p = k − L (Q is the suffix): the `W` rule above.

Discovering the flank ranges is **branching work, even for an exact DNA pattern**, and it is not bounded by the
suffix count: with k = 5 and the record `ACGTT`, the motif `AC` has zero suffix anchors and one context at offset
0. The engine therefore counts the ranges it visits (`ranges_visited`) against `max_steps`, and the answer reports
them; no closed-form bound is claimed. Per offset the counts are exact once discovery completed (version 1 tracks
no per-offset completeness under a discovery stop: every per-offset count is then `at_least`). This is
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
    suffix range (for the record `ACGA` at k = 4 this build's PRIMARY graph stores `TCGT`, and the wrapper
    exposes `ACGA` as the virtual k-mer; `CGA`, the suffix of the virtual k-mer, is found neither as stored
    suffix `CGA` nor as stored suffix `TCG`; PatternSearch.SuffixOnWrappedPrimary checks which orientation is
    stored and asks the virtual k-mer's suffix). v1 does not offer `suffix` on wrapped
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
context. The search stops at the first disallowed base, so it is complete and never scores anything. As built
(review GPT-3, SPEC §18), the clock is read before each anchor's extension (its spelling, k − 1 BOSS steps that no
step charges, and its search), every 64 anchors listed and before every 64th node expanded, besides the
4,096-step readings: an extension whose 1,129 steps never reached one ran past the work time on refseq33m. The
answer states `work.extension_anchors` (the anchors whose extension began) and `work.extension_branches` (the
nodes with two or more allowed successors) beside `extension_edges`.

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
part). v1 anchors each oriented pattern on its own [0, k); it states the bits of P[0, k)
(`anchor_information_bits`) and of the least informative searched window (`min_anchor_information_bits`, an
addition of the review of 2026-10-07), which the information floor gates on (SPEC §7.7).

### 4.3 Labels, coordinates and placement (reused from the traversal, not from the aligner)

**The label-free path first (`output.labels: none`).** Nothing else in this section runs for it. A request with
`output.labels: none` returns the contexts themselves — k-mer, instance, offset, strand, the graph node id and
the annotation row id (`AnnotatedDBG`'s graph-to-annotation mapping, so the row can be named later), and for
`long` the path — without reading one annotation row: counting and extraction stay graph work, bounded by
`max_steps`, the deadline (read within 4,096 steps or 64 released contexts) and `max_contexts` (`partial` retains
only what it can release); version 1 has no memory account on this path, its memory being bounded by those caps
(SPEC §7.6).
It is available on every backend (no `decode_charged` requirement) and every graph mode (on canonical indexes
the contexts come without placement, as everywhere). It is the cheap retrieval: labels cost a row decode per
context, and most follow-up questions — seed `/traverse`, inspect the variants, choose rows — need the k-mers,
not the labels. A later increment reads the labels of **given row ids** (§12), so an agent can count, extract,
choose, and only then pay for annotation on the rows it chose. `all_or_count` applies its threshold to the
contexts as usual; `retrieval_complete` then means every context is returned, and the answer claims nothing
about any label. The engine's `enumerate()` is this path; the two steps below are what `labels: all` and
`predicate_only` add on top of it.

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
verified while another only carries the constituent k-mers. As built (review GPT-3, SPEC §12.1, §18) each
(path, label) is verified once: its chains are the intersection of its k-mers' coordinate lists, each shifted by
its k-mer's position, found by a leapfrog join from the shortest list with galloping seeks; consecutive chains are
one run (a homopolymer's path is a few units of work per k-mer, not one per coordinate), each run is placed record
by record, and the runs are kept (charged to the memory account) for the output, which no longer joins again;
every seek, run and record is clocked. Each returned label states its `support`:
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

### 4.4 Graphs without a mask (the owner's decisions #16 and #17, 2026-10-08)

Some graphs have masks, some not: the route serves both. With the mask (`mask: file | built_at_load`,
`counting: exact`) nothing changes. Without it (`mask: absent`, `counting: upper_bound`) the engine
(`PatternSearch`, "Graphs without the dummy-edge mask") discovers the same ranges with the same steps, and every
place that read the mask counts candidates instead:

- A flank range counts its entries whose W is not `$` (`DBGSuccinct::count_non_sink_edges_in_range`: ranks of W,
  so a sink dummy never counts), the last position its entries with W the allowed base, plain or marked
  (`count_edges_with_symbol`); no deferred scan of masked edges exists.
- **Checked and unchecked ranges.** A source dummy is `$^j` followed by the first k − j bases of a sequence start
  no k-mer enters. A range whose nodes the search spelled whole (depth k − 1: every node symbol a pattern base or a
  flank base) cannot hold one, so its count is exact. Any other range (an offset p > 0 of `any_offset`, whose
  first p node symbols are not spelled; `suffix` scope with L < k; a window whose leading N run was skipped) is
  **unchecked**: its candidates go into the upper bound U and not into the lower bound. On an even-k wrapped
  PRIMARY graph the palindrome scans spell every candidate they check, which settles it too.
- **The count** is `exact` when nothing of it is unchecked, when U = 0 (an empty block: absence holds), or after a
  release that enumerated every candidate; otherwise `bounds` {lower, U}. Stops leave `at_least` and `unknown` as
  on a masked graph. **Tiny blocks are checked** (the owner's decision #24): when discovery and its scans complete
  with at most `--pattern-max-checked-entries` (default 50, at most 1,000, 0 off; `caps.max_checked_entries`)
  unchecked candidates in all, each is tested with `node_has_sentinel` after the scans (k − 1 steps each, phase
  `mask_scan`), so every count of the pattern is `exact` in every mode; larger blocks are not touched (no sampling).
- **The estimate** U × f, rounded and kept inside [lower, U], is stated beside every `bounds` count (the route,
  `count_json`), with the note `estimate_sampled_dummy_fraction`. f is the fraction of real k-mers among the
  graph's entries whose W is not `$`, sampled once per graph at load (`pattern::sample_real_fraction`: 10,000
  entries drawn uniformly with replacement, `std::mt19937_64` seeded with the number of edges, so the same graph
  gives the same f after every restart; a drawn entry with W = `$` is drawn again; an entry is a source dummy iff
  its node holds `$`, found by `BOSS::node_has_sentinel` within k − 1 symbols), stated with Wilson's 95% interval
  and its source (`sampled`) in the capabilities and in every answer's `index`. On the mini index: f 0.9999
  [0.999434, 0.999982], the exact f (`stats --count-dummy`: 375 source and 12 sink dummies among 8,335,760
  edges) 0.999955 inside it; 40–60 ms at load. The estimate is not a bound: a pattern at a record start, held by
  the dummies before that start, can have a true count far below it (the mini's island start
  `GATGCCGGTGAACAAC`: U = 32, estimate 32, 17 contexts).
- **Admission on U** (conservative): `all_or_count`'s threshold, `stop_at_threshold` (on the running U) and the
  extension's admission (`max_anchors`) compare the upper bound, and the note `threshold_upper_bound` says when such
  a decision went against the request while the lower bound did not cross.
- **Lists stay exact.** The release tests each unchecked candidate with `node_has_sentinel` (at most k − 1
  symbols, about the spelling the route does per context anyway, under the deadline and not charged as steps) and
  never releases a source dummy. A release that enumerated every candidate (`all_or_count` within U; `partial`
  drained; a long pattern's anchors listed before their extension) makes the counts `exact`; a cut `partial`
  release raises each lower bound to what it released. The release, the palindrome scans and the check step over
  a range's sink edges (W = `$`) with `DBGSuccinct::next_non_sink_edge`: W read at the first 16, then a jump to the
  next symbol that is not `$` by rank and select, so a run of sinks costs the same whatever its length (review
  GPT-3; 100,000 sinks had cost 2.5 ms, past a work time of 1 ms, one edge at a time).
- **Stated limitations**: without the mask `partial` keeps every unchecked range it discovers (24 bytes each, at
  most one per step) rather than about `max_contexts`; on an even-k wrapped PRIMARY graph `stop_at_threshold` can
  fire late or not at all (the running U leaves out ranges whose palindrome scan is pending; the admission after
  discovery compares the final U); mode `count` states `bounds` where a retrieval of the same request states
  `exact`.

**Identity (decision #17).** The mask (and the `.bloom` the loader reads with it) is derived data of the graph,
closely linked to it: it changes which results are exact, never the number of an exact count (an exact count is
the same with and without the mask). It is not part of `index_fp` and a manifest never lists it
(`SPEC-labeled-traversal-core.md`, "The index identity"): masking a deployed index keeps its `index_fp`. The
limitation that comes back with it: no fingerprint covers derived data, so a stale or foreign Bloom filter of the
same k and mode beside the graph is not detected by `index_fp` (a stale mask is caught only by its size and by the
`mask_invalid` check).

Tests (independent oracles): `PatternUnmasked` builds every graph twice from the same records, one with its mask
reset, and checks U against a scan of the unmasked graph's entries, the lower bound and the exact counts against a
graph walk of the masked twin and a scan of the records, the lists against the masked twin's, and the sampled f
against BOSS's own dummy-tree count; the route's `PatternRoute.Unmasked*` and `PatternMaskUnmasked.*` and the
integration's `TestPatternMini.test_unmasked_*` check the same on the real loader and the mini index.

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
| `count` | discover and count the contexts (or anchors and paths) in the requested scope. **Without a predicate** no annotation row is decoded, and completed interval descriptors are **discarded** as they are counted, so the memory account bounds only the DFS frontier and broad patterns can be counted (as built: except the ranges awaiting a deferred scan, none in the common case but on an even-k wrapped PRIMARY graph one per range at a palindrome-capable offset, up to `max_steps` of 88 bytes, SPEC §7.6). **With a predicate** (§5.6) the answer carries the raw and the selected counts: descriptors are retained up to `max_predicate_contexts`, the predicate pass reads annotation under its budgets and the backend restrictions of §4.3, and still no result is returned |
| `all_or_count` (default) | the count, with the interval descriptors **retained while retrieval is still possible**: without a predicate until the raw count exceeds `max_contexts`, with a predicate until it exceeds `max_predicate_contexts` (20,000 raw contexts may select 5); past that they are discarded and counting continues within the compute budget; retrieval of **all** requested results only if the completed aggregate count — the selected count when there is a predicate — is within the retrieval threshold and the retrieval finishes on every shard; otherwise the counts alone, with the reason |
| `partial` | the count, then retrieval up to the caps with what was built delivered and every cut stated — an agent that wants the first results of a large set asks for it explicitly |

**Projection `none`.** With `output.labels: none` the retrieval of `all_or_count` and `partial` is the
label-free extraction of §4.3: the contexts with their k-mers, offsets, strands, node and row ids, no annotation
read, no annotation budget touched. The same thresholds and barriers apply to the contexts; `rows_refused` and
`anchors_truncated` stay empty by construction.

**Descriptor retention by phase.** For L ≤ k the descriptors are the contexts' intervals, kept or discarded as
the table says. For `long` the phases differ: **anchors** are always kept through the extension admission
(they are what the extension starts from) and discarded once every anchor has been extended or the admission
failed; **completed paths** are discarded at once in an unfiltered `count`, kept up to `max_paths` for ordinary
retrieval, and kept up to `max_predicate_contexts` when a predicate follows. `stop_at_threshold: true` ends the
current phase's discovery as soon as the running lower bound of what it counted crosses its own threshold, and
names it; the deferred scans (§4.1) do not consult it (the owner's decision of 2026-10-07), so a count that
crosses the threshold only in them, or on an even-k wrapped PRIMARY graph only in the union of the two base
searches, ends `exact` with `count_above_threshold`. Without a predicate the thresholds are
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
| the running lower bound crosses the threshold in discovery and `stop_at_threshold` is set | `at_least` counts, results withheld (`withheld: threshold_crossed`) |
| range or path discovery reaches `max_steps` or the deadline | `at_least`, `bounds` or `unknown` counts, results withheld (`withheld: discovery_budget`) |
| graph counts complete, but a shard's annotation or output budget is exceeded, an anchor's labels are truncated, or a row is refused | exact graph counts, results withheld with the reason and the shard (`withheld: annotation_budget` with `rows_refused`, `anchor_labels_truncated`, or `output_budget`) |

An agent can then tell "many matches" (count above threshold: refine the pattern, the scope or the strand),
"expensive discovery" (discovery budget: shorten the flank scope, add information) and "expensive annotation"
(a common anchor: raise `max_labels_per_anchor` for a narrower pattern, or ask for the count only) apart, and
act on each.

### 5.3 Budgets, the deadline and loading

**The deadline.** One per request, started when the request is parsed and compared at the start of every
pattern and before every release, and inside every loop: the range DFS, the deferred scans and the extension
every 4,096 steps (one step count, so the boundary from discovery to the scans has no reading of its own); the
release every 4,096 descriptors prepared and edges examined, and every 64 contexts handed to the route, as
`all_or_count`'s delivery of its buffered release; the paced reads at their run boundaries as in `/resolve` and
`/traverse`, and before the labels of each context are built; serialisation every 4,096 objects. As built since
the review GPT-3 (SPEC §7.6, §18) also: the extension before each anchor, every 64 anchors listed and every 64th
node expanded (§4.2); the low-complexity diagnostic between its pieces (below); and the labelled retrieval's own
work between its reads — the occurrences of each context or path, the paths' label lists and their verification —
at least every 4,096 units, before the work. A stop there is stated (`{extension | placement | output, time}`),
never a 503 after the budget. The route keeps
a **finalisation reserve** (`--pattern-finalize-ms`, default 250 ms, stated in the capabilities) and, as built,
more: the time it estimates writing the answer built so far will take, from the bytes of its results and the
configured delivery rates (SPEC §7.6, as `/traverse`'s delivery reserve). The count answer of every pattern is
maintained incrementally during the search, so that when the deadline less that time passes, the route stops all
work, early enough for the answer it holds, and serialises it, with `stop: {phase, shard}`, `withheld: deadline`
and the relations of §3. If serialisation itself overruns the reserve, the answer is the
explicit outcome 503 `deadline`, as a `/traverse` attempt past its bound; nothing partial is sent as if whole.

**Loads as `/search` loads (`in_ram`).** The owner's decision of 2026-10-09: the route follows `/search`'s
choices ("some indexes are loaded in RAM and some not"). Every production server maps its graphs at start-up
(`--mmap`); `in_ram: true` loads the selected pair into RAM for the request when it fits `--mem-cap-gb`, with
`/search`'s reservation and wait, and the request's resource limits (time, memory, work) apply to its work only,
after the load, which is unavoidable and cannot be interrupted at the deadline; `timing.load_ms` states it (SPEC
§24.2). Without `in_ram` the route serves the index the process holds.

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
own; what any item would push over the share is a stated memory stop. (As built, increment 3: the account is a
deterministic model of the labelled retrieval of `output.labels: all`, SPEC §14.4; the label-free path has no
account, its descriptors being bounded by the caps, SPEC §7.6.)

| cap | default | server flag | what is stated |
|---|---|---|---|
| `min_information_bits` | 24 (≈ 12 specified bases) | `--pattern-min-information-bits` | gates **discovery** for patterns with an ambiguity code, `protein`, `any_offset` and `long` (every searched orientation's anchor window, §3); an exact pattern (every position one base, whatever its kind) in `suffix` scope is always admitted, however short: one range, a few ranks, and a scan of the range's invalid edges when it has any, charged like any step |
| `max_contexts` | 10,000 per pattern, all offsets and both strands, aggregate over shards | `--pattern-max-contexts` | the retrieval threshold of `all_or_count`; the count is always returned with its relation |
| `max_anchors` (`long`) | 1,000 | `--pattern-max-anchors` | the extension threshold; `withheld: anchors_above_threshold` |
| `max_paths` (`long`) | 1,000 | `--pattern-max-paths` | the retrieval threshold on completed paths; `candidates_examined` |
| `max_steps` | 100,000,000 range, scan or edge steps per shard (measured in RAM at 0.2–0.7 µs a step on the mini index and on random graphs of 0.5 and 2 billion edges, an M-series Mac: 10^8 steps take 20–75 s; unmeasured on a deployed mmapped index until the benchmark of §13 increment 6; which of this cap and the time budget stops a heavy request first depends on the machine and its load; a request cannot raise it) | `--pattern-max-steps` | `ranges_visited`, `steps`, the phase that hit it |
| `max_labels_per_anchor` | 64 | `--pattern-max-labels-per-anchor` | `anchors_truncated` with the cap and each row's total |
| `max_predicate_contexts` (§5.6) | 100,000 per pattern | `--pattern-max-predicate-contexts` | the compute admission of a predicate pass; `withheld: predicate_above_threshold` with the unfiltered count |
| `max_predicate_work` (§5.6) | the oracle's units, default as `max_annotation_work` | `--pattern-max-predicate-work` | `selection.tested` as `at_least`, `withheld: predicate_budget` |
| `max_predicate_labels` (§5.6) | 10,000 names per request | `--pattern-max-predicate-labels` | refused above it (400) |
| `max_annotation_work` | the oracle's units, default as `/resolve`'s | `--pattern-max-annotation-work` | `rows_refused` with their anchors |
| `max_memory_mb` | 256 per request, split per shard | `--pattern-max-memory-mb` | `memory.stop` with the phase and the shard |
| `max_labels` | 1,000 | `--pattern-max-labels` | `labels` count with its relation; labels kept by (contexts desc, column asc), as `/search`'s top-N; `partial` only (`all_or_count` returns all or none) |
| `max_occurrences_per_label` | 16 | `--pattern-max-occurrences` | `occurrences` count per label with its relation; `partial` only |
| `time_budget_ms` | 60,000 (the owner, 2026-10-07: "5 s is not sufficient in general"; some `/search` calls take longer today, and this is the more complex function) | `--pattern-max-time-ms`, default **600,000**: under the 900 s content timeout with room for serialisation and compression (at most 899,000 on `server_query`, the content timeout less 1 s); the service sends an explicit budget per task under it | `stop: time` with the phase and the shard; the finalisation reserve inside it. Nothing in the engine assumes a short run: the clock is read every 4,096 steps (and 64 released contexts) whatever the budget; the retained descriptors are bounded by their caps, not by time; `max_memory_mb` bounds the labelled retrieval (there is no memory account on the label-free path in version 1); `max_steps` does not grow with the budget, so a long budget is usually ended by steps first; and a long call costs the one server thread it holds |

Why information rather than length, and why only as a gate: expected suffix anchors ≈ (k-mers in the graph) ×
2^−bits. On a 50-billion-k-mer graph a 16-mer expects about 12 per strand, a 12-mer about 3,000, a 10-mer about
48,000. That is a planning heuristic; what bounds a run is the range work and the deadline, which is why
`count` is always available and `max_steps` is the real admission for discovery. The range work is set by
where a pattern's N runs sit, not by its bits (review of 2026-10-07, X-EFFICIENCY-01): a run costs about
min(4^run, edges / 4^a) ranges per level, a being the specified bases before it in the searched window, so a
40-bit pattern with a leading N^11 cost 8.6 million steps on the mini index where its 20-mer alone costs 61. As
built, a run at the start of a searched window is skipped on a `$ACGT` graph (the rest of the window is searched
with its offsets shifted: the same contexts), and a pattern's base searches run cheapest first by this
estimate, so that a budget stop leaves the cheaper orientation complete; a run inside the pattern stays a cost,
and an internal anchor (§12) is what would remove it.

Low-complexity patterns are not refused: they are legitimate questions, the caps bound them, and
`is_low_complexity` (`aligner_seeder_methods.cpp:22`) marks them in `notes` so the reader knows why the counts
are large. As built (review GPT-3, SPEC §7.8, §18): on a completed search only (never beside a stop), with sdust
read in pieces of 128 + 63 bases that stop at the first flagged one and read the clock before every piece but the
first; one sdust over a 30,000-base repeat had taken 0.8 s without a reading. A pattern of more than 191 bases
whose diagnostic the work time cuts loses the note and says `time_limited` with no stop.

### 5.4 Deduplication

Graph contexts are distinct by (orientation, path, offset). Placed occurrences are distinct by (column, `seq_id`,
start, strand): the up to k − L + 1 contexts of one interior occurrence collapse to one placed occurrence, and
two offsets in one k-mer stay two. Node id is never a deduplication key.

### 5.5 Determinism

Contexts are enumerated in BOSS edge order, then by offset; labels are ordered by (contexts desc, column asc);
occurrences by (`seq_id`, start, strand); shards by (graph name, `index_fp`) (§8). The same request on the same
index gives the same answer, including which results a cap cut, because every budget a shard spends is its own
fixed share. The exception is a deadline stop, which depends on the machine and is stated
(`determinism: time_limited`); as built, also a low-complexity diagnostic the deadline cut, whose answer is
`time_limited` with no stop and complete counts (SPEC §7.9).

### 5.6 Annotation predicates: select contexts by a condition on their labels (extension, the owner's request)

**As built (increment 5b, 2026-10-08; the contract is `SPEC-pattern-search.md` §19).** The owner's decisions
P1–P24 of 2026-10-08 change this section where it said otherwise; the text below keeps the design and is marked
where it differs:
- **Strands, consistent per context.** On a BASIC graph a k-mer's row lists the columns whose records hold it as
  deposited; a record holding the motif on its other strand annotates the reverse complement. A predicate is
  evaluated with `predicate_strands: "either"` (the default) on the labels of the context's k-mer **or** of its
  reverse complement — a label supports a context in one orientation as a whole, never by a mix of strands (the
  owner's answer to P11; for a supported path, every k-mer of the walk as spelled, or every k-mer of its reverse
  walk, SPEC §20.9). `"context"` reads the context's own row only, for stranded indexes; CANONICAL and PRIMARY graphs share
  one row and answer `"either"`. The answer states per selected result and label which orientation supported it
  (`selection_strands`: `"context"`, `"reverse_complement"`, `"both"`, or `"either"` where one row serves both).
- **The default projection stays `none`** (P2): `output.labels` omitted is `none` with a predicate too; a client
  that wants the predicate's labels sends `"predicate_only"` (the "(the default)" below is the design's, not the
  built contract's).
- **The narrowed absence claim** is `absence_filter: "predicate"` beside `absence_scope`, which keeps its closed
  values (P15).
- **Counts**: `counts.tested` and `counts.selected` with relations; a pass stopped early states `selected` as
  `bounds` [S, S + R_upper − T] where the raw count has an upper bound (P13), not `at_least`. The per-entry object
  is `selection` {`pass`, `support`, `access`}; the request-wide part is the answer's `predicate` block (normal
  form, unknown labels, `vacuous`, scope, strands); unknown labels are one list (a single graph), not per shard.
- **Units**: a membership read is charged its whole decoded row (`KeyCost::entries`, P17) under its own
  `max_predicate_work`, separate from `max_annotation_work` (P3).
- **Long patterns**: a predicate selects among **supported paths** only (`long_search: "supported_paths"`,
  increment 5s, built in round C, SPEC §20; P23, P24): the search carries each branch's support (the labels
  carrying every k-mer, or the records holding the walk whole) and prunes a branch nothing supports; the
  predicate is decided at a walk's completion on its support at the level searched, a monotone normal form also
  pruning, and `"either"` holding the walks until their mirrors' support is known (found by the search on the
  other strand, or read). The two-step raw-path design of "`long`: support before selection" below is not built;
  with `long_search: "paths"` a predicate is refused, with `"anchors"` a long pattern's selection is
  `not_started`.
- **Motif-level predicates** (§12's item, built in round C for patterns of L ≤ k, SPEC §25): with
  `predicate_scope: "motif"` the normal form is also evaluated once per pattern on the union of its contexts'
  labels (each context's set as the selection evaluates it), from the pass's own reads: exact when every context
  was tested, else Kleene's value on the labels found (definite only when the untested contexts cannot change it).
  "Present in A, absent throughout C" on one index; across chunks the client evaluates the predicate on the union
  of the chunks' labels.
- **Labels are column names only** (P6): no taxonomy, a cohort is an explicit list (expanded by the client);
  record headers are unknown names (P20).

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
  placement exists (the design's default; as built the default stays `none`, P2); `all` runs the ordinary
  discovery of §4.3 for the selected contexts, under the per-anchor
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
`counts.selected` as `bounds` [S, S + R_upper − T] where the raw count has an upper bound, `at_least` otherwise
(as built, P13), and withholds results (`withheld: predicate_budget`).
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
translation filter, and no branch through a stop codon except where the peptide has `*`. The standard code
(table 1) is built in; `genetic_code` is a request field for others. As served (increment 5, the owner's decision
#15 of 2026-10-08; SPEC §12.2): every NCBI translation table (1–6, 9–16, 21–33), the ambiguity codes X (every
codon that is not a stop), B, Z and J, and since decision #19 the stop `*`: the stop codons of the table at that
position (X still never matches one). Tables 27, 28 and 31 have no unconditional stop codon (their context stops
code their residue, decision #21): there `*` matches nothing, the peptide has no instance, and every entry of it
says so (note `no_stop_codon`; L ≤ k is answered `exact` 0 without a search). (`4596bb3b` refused `*` in the
pattern's slot, `stop_unsupported`; that code is retired, `4596bb3b` the one build that answers it.)

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
  "max_contexts": 10000, "max_anchors": 1000, "max_paths": 1000, "max_steps": 100000000,
  "max_labels_per_anchor": 64, "max_annotation_work": null, "max_memory_mb": 256, "time_budget_ms": 60000,
  "max_labels": 1000, "max_occurrences_per_label": 16,
  "allow_unbudgeted_annotation": false,
  "require_support": null,
  "predicate": {"and": [{"any": ["SRR100", "SRR101"]}, {"not": {"any": ["SRR900"]}}]},
  "max_predicate_contexts": 100000, "max_predicate_work": null,
  "graphs": ["name"],
  "output": {"occurrences": true, "labels": "predicate_only", "paths": false}
}
```

The values are illustrative (here the defaults of §5.3). `predicate` is optional (§5.6); without it
`output.labels` defaults to `all`, and the predicate caps are unused (as built, the default is `none`, frozen for
contract version 1: SPEC §1, §13).
`output.labels: none` is the label-free path (§4.3): contexts with k-mers, offsets, strands, node and row ids,
no annotation read; it is the cheapest retrieval and the one a client should ask for first.

Exactly one of `dna`, `iupac`, `protein` per pattern; at most `--pattern-max-patterns` (default 16) per request;
`id` optional. Unknown fields are refused (400), as `/traverse` refuses them; `in_ram` is accepted as `/search`
accepts it (§5.3). The route
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
     "kmer": "TTGACGCGTCAGGTACCGGAACTGATCCAGT", "instance": "GCGTCAGGTACCGGA", "offset": 5, "strand": "+",
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
  `get_node_sequence`; for `long`, a path, the new fields `sequence`, the spelled sequence of the whole path, L
  bases, and `anchor_kmer`, its anchor's k bases — never `kmer`, which keeps its contract-version-1 meaning, a
  context's k-mer; paths only for a request with `long_search: "paths"`, an opt-in: the owner's decision of
  2026-10-07, SPEC §12), the
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
- With `output.labels: none` each result carries `kmer`, `instance`, `offset`, `strand`, `node` (the graph node
  id, scoped by `graph`) and `row` (the annotation row id) and no `labels`; `by_label` is absent, and
  `counts.labels` and `counts.occurrences` are `unknown` by construction, since nothing was read.
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
(`server.cpp:1421-1480`), with a `pattern` block: the floors, caps and the finalisation reserve (with the delivery rates of the time kept back for the answer,
`delivery_mbps`; the prose rules are references to the SPEC since the owner's decision P9, so that the document
an MCP tool returns whole keeps 1 KiB under its 32 KiB ceiling), the modes, the projections (`none`,
`all`, `predicate_only`, with `none` available everywhere) and scopes (and which scope each graph mode supports), `placement` (the best this index can give), `strand_stated`,
`mask: file | built_at_load | absent` with `counting: exact | upper_bound` and the sampled `dummy_fraction` (§4.4), `annotation: budgeted | unbudgeted`, the build alphabet and `graph_cleaned`, the
resident graphs, `records_shorter_than_k: not_indexed`, and `pattern_contract_version: 1`. The full block also
has a route of its own, `GET /pattern/capabilities` (the owner's decision of 2026-10-09, SPEC §23), and
`GET /traverse/capabilities`, the document the search service's probe reads, carries the block with `details`
naming that route: in two phases, first the full block, later only the fields a client gates on, so that the
probe's document stays small while the block grows (SPEC §23.2). A client gates on the block (a contract version it implements, and `available: true`; SPEC
§1, §10.3: since the outside review of 2026-10-07 a higher version is refused, not read by presence), reads
fields by presence within it, and offers the label projections exactly when the block's `projections` list has
them (§9).

## 8. Multi-graph servers

As decided (the owner, 2026-10-09: "the same logic as for the general search", "no new recipe"; SPEC §24):
`/pattern` takes `/search`'s `graphs`, answers each selected (graph, annotation) pair as a single-graph server
answers it (its own deadline, caps and memory account: the request's budgets apply per pair), in parallel on
`graphs_pool`, and returns the answers concatenated, each tagged with its pair and `index_fp`. The requester
merges, as the search service merges `/search`'s answers of a chunked database, sending one graph per request
(PROMPT §3.1 items 2 and 3); `GET /capabilities` states every pair (`graph_summary`) and whether the pairs' columns
are disjoint, which licenses summing a label's counts over them.

The merged answer below was the earlier design; it is not built:

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

Two layers. The server route is request/response; the search service wraps the route in
a **job** (submit → status → results), as it wraps `/search` and the traversal, the owner's decision of
2026-10-07: the backend stays synchronous, the service is the asynchronous layer.

**Server.** `POST /pattern` goes through `process_request` like `/search` (compact JSON, gzip when accepted,
`kContentTimeoutS` 900 s, `server.cpp:52`). The route's own deadline (default 60 s, cap `--pattern-max-time-ms`
600 s, at most 899 s on `server_query`; finalisation reserve inside it, §5.3) stays under that wall with room for
serialisation and compression. As `/resolve` and `/traverse`, it stops a request whose client has left, or that
runs at shutdown, at its next clock reading and writes nothing, not even an error (a half-close counts as gone;
`ResponseControl::gone`, which only `/pattern` sets: `/resolve` and `/traverse` still write their errors).

**No Python client yet** (the owner, 2026-10-07: "write it lazily, when we need it"; unused code is weight to
maintain). Nothing calls one: the search service has its own HTTP client, agents reach the route through the
service's tools, and the tests call the route directly. When a user needs it, it is one method per class in
`api/python/metagraph/client.py` beside `align()` on `_json_seq_query`. The same holds for passing
`time_budget_ms` through the library's MCP `traverse_resolve`.

**Search service** (its repository; the service session builds it, org-id integrates;
`docs/PROMPT-search-service-pattern.md` is the request). A **kind of the existing search job** (the owner:
"it is essentially a search, and the results are like search too"): the same Search/Task/Result tables, the same
leaf split (one task per (database, graph chunk, pattern chunk of at most the host's `max_patterns`)), the same
semaphore, per-database caps, 1,200 s client timeout, status, results, CSV, S3, lock and retention; the worker
calls `POST /pattern` instead of `POST /search`, skips scoring, truncation and enrichment for that kind, and
stores contexts as result rows with the per-pattern summary per task. The merged view adds counts with the
relation algebra of SPEC §7.4 (`exact` 7 + `unknown` = `at_least` 7; `bounds` [3, 5] + `exact` 7 = `bounds`
[10, 12]), never sums label counts across tasks, and states `retrieval_complete` only when every contributing
task's is true; `withheld` and `stop` per task (the PROMPT, §3.1 item 2). A chunked database on the service side is many **graphs on one multi-graph
server process**, selected per task through the request's `graphs: ["{label}-{i}/{N}"]` as `/search` selects a
shard today; the job type is built and tested on refseq33m-experimental (one graph) with milestone 1, and serving
the chunked databases needs milestone 6, whose `graphs` selection keeps exactly `/search`'s names and semantics
(§8):

- admission by the queue, **shared with search and traversal** (the owner: "they use essentially the same
  resources: load an index and access it"): the service's per-database queues, `META_DB_CAPS` and its distributed
  semaphore count search, pattern and traversal calls against **one** cap per database server, and on the server
  side the one request pool (`-p`; `--threads-each` is the threads within a request) serves `/search`,
  `/pattern`, `/traverse` and `/resolve` alike — no per-route pool, no reservation. Every route, the GET
  capabilities routes included, runs on those `-p` threads: with `-p` long calls in flight the probe waits up to
  their budgets, so the cap stays below `-p` (`-p` ≥ cap + 1, as `SPEC-labeled-traversal-core.md` "Deployment"
  requires for traversals), and a full pool is not a dead host. The engine needs nothing beyond "at most cap
  concurrent calls per database server across all routes"; a request's memory is bounded with `output.labels:
  all` by `max_memory_mb`, and on the label-free path by the caps (about `max_contexts` descriptors in `partial`
  and `all_or_count`, the DFS frontier in `count`, and in every mode on an even-k PRIMARY index up to
  `max_steps` × 88 bytes of ranges awaiting their palindrome check; SPEC §7.6), so the server's worst case for
  pattern calls is cap × the larger;
- caps at admission, never truncation afterwards: `max_contexts`, `max_patterns` and `time_budget_ms` above the
  service's ceilings are refused with a 400 naming the field; an answer is passed through whole, since cutting it
  would falsify `retrieval_complete`; the backend called from the API process on the thread pool
  (`anyio.to_thread.run_sync`), one host per database, timeout = the request's `time_budget_ms` plus an allowance;
- served only for databases whose host carries the `pattern` block (on `/traverse/capabilities`, which the
  service's probe reads, as well as on `/capabilities`, §7.3) with a contract version the service knows or a
  higher one, fields read by presence; a lower or missing version answers `pattern_unsupported` with the host's
  feature list, as `traversal_disabled` does; the label projections are offered exactly when the host's
  `projections` list has them, so they switch on by themselves when a later milestone lands;
- the answer passed through with the service's additions only: `database`, the untrusted-data notice on label
  and record strings, the standard error envelope (no `hit_unit`: a context is not a hit, and the tool says what
  it is instead — one distinct k-mer of the index containing the pattern at one offset and orientation); the
  tool's docstring explains the `withheld` reasons, `retrieval_complete` and `absence_scope`, what to change for
  each, uses the service's strand vocabulary, and says that canonical and primary indexes report contexts only;
- the tool's default mode is `count` (the route's own default is `all_or_count`), `scope` defaults to
  `any_offset`, `partial` is exposed, and `labels: none` is the cheap second step: the k-mers without the
  annotation, from which the agent picks rows before asking for labels; rate limits and the anonymous-call
  budget as for the other synchronous tools; no caching of answers in v1 (`determinism: time_limited` matters
  only once an answer is cached).

Budgets exactly like search: the service sends `time_budget_ms` equal to the host's cap (600 s, §5.3) unless the
caller asks for less, and its 1,200 s client timeout (`META_CALL_TIMEOUT_SECS`) stays above it; a
job-originated call holds no client connection, so a long budget costs only the request-pool slot it occupies.

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
| Python client | `api/python/metagraph/client.py` | `align()`, `_json_seq_query` | deferred until a user needs it (§9) |

## 11. Not breaking the aligner

- No scoring, seeding or extension default changes. `initialize_aligner_config` (`align.cpp:33-69`),
  `DBGAlignerConfig`, `SuffixSeeder`, `DefaultColumnExtender`, `LabeledAligner`, `AnnotationBuffer` keep their
  code paths; the only edits inside alignment files are the template signature with a default that reproduces
  today's symbol loop on every alphabet (§10) and the move of two JSON writers into a shared header.
- **Golden gate for alignment:** before the first change, record the outputs of every
  `integration_tests/test_align.py` case (8 tests × graph representations, DNA5 included), of `test_query.py`'s
  `align` cases and of `tests/annotation/test_aligner_labeled.cpp` as fixtures; the gate compares byte for byte
  after every increment. As run (milestone 1 and the review's fixes): on a DNA4 build only — no DNA5 build was
  made —, by scripts kept outside the repository (a manual gate, not in CI).
- The seeder's seed lists on the test graphs are captured once (a small gtest dumping `SuffixSeeder` seeds for
  the mini graphs, a DNA5 graph among them) and compared after the `allowed` change. As built, the unit test
  `PatternSearch.SuffixToPrefixDefaultUnchanged` checks instead that the seeder's calls (without a symbol set)
  visit the same nodes in the same order as before, on DNA4 graphs.
- `/align` and `/search` with `align: true` answer byte-identically on the fixture requests (the byte-identity
  harness of the traversal work, `scripts/traversal/bench_traverse.py --compare`, with an alignment request set).

## 12. Later increments on the same engine

- **A one-pass paced row/tuple visitor** replacing the discovery-then-placement pair (§4.3), extracted from the
  two fetch loops.
- **Motif-level predicates** (§5.6): built for patterns of L ≤ k on one index (round C, SPEC §25,
  `predicate_scope: "motif"`); left: patterns longer than k (on their supported paths), and a server-side merge
  across shards (today the client's: the predicate on the union of the chunks' labels).
- **Supported paths** (moved into increment 5s, built in round C, SPEC §20): what remains for later is their
  **selective anchor** (anchoring each orientation on its least ambiguous k-window and extending both ways, DECISIONS
  P30: a pattern with an ambiguous end is still `anchors_above_threshold` today) and an **alignment projection**
  of the supported paths (consensus and variant columns, P27). Counting walks by dynamic programming is not planned
  (later at most: the enumeration counts millions of walks in seconds, and `supported_paths` is the count that
  answers the question).
- **Paced, budgeted selected-column access** for predicates (`predicate.access: columns` under a budget), so that
  `any(A) and none(B)` on a column-major backend reads two bits instead of a row.
- **Labels for given rows:** a request naming the row ids of contexts an earlier label-free answer returned,
  answered with their labels (and placement where it exists) under the budgets of §4.3, so that annotation is
  paid for only on the rows the agent chose.
- **`suffix` scope on wrapped PRIMARY graphs** through a prefix lookup with orientation conversion.
- **Bounded mismatches** (`NOTE-mismatch-search.md`): automaton states (position, mismatches used) in the DFS;
  the admission policy and caps carry over.
- **Tolerant best-hit search on `/align`:** IUPAC rows in the extender's score table
  (`aligner_extender_methods.cpp:38-59`) plus `min_seed_length` per request and the labeled aligner on the server
  route. Best hits, not every carrier; a separate item.
- **Internal anchors:** for `long`, anchor on the most informative k-window and verify both ways; for L ≤ k,
  start the search after an N run inside the pattern, so that the run is matched on narrow ranges (as built, a run
  at the start of a window is skipped on `$ACGT`, but a run inside costs about min(4^run, edges / 4^a) ranges per
  level, §5.3).
- **Placement on canonical indexes** for `long` from coordinate steps, shared with `/search`'s `query_coords`.
- **Queue path** in the search service for batches and archive-wide scans.

## 13. Implementation increments (each with tests, then an adversarial review)

0. **The oracle suite first.** A gtest suite of tiny graphs built from explicit records, with **two** brute-force
   oracles, because no single one validates both counting units: a **graph-walk oracle** that enumerates the
   built graph's k-mers and paths directly (every k-mer containing the pattern at each offset, every path
   spelling a long pattern), giving the expected raw graph contexts; and a **record-scan oracle** that scans the
   records' retained islands as strings, both strands, matching the pattern position by position through the
   test's own IUPAC table (the bases of each code; the reverse complement and palindromy derived from it, never
   from the engine's `Pattern`) or the codon sets, giving the expected placed occurrences and verified supports.
   (As built: the graph-walk oracle shares the valid-edge mask, `in_graph`, and the node ids with the engine; the
   record-scan oracle shares neither. The IUPAC semantics are pinned by the unit suite's own table for all 15
   codes, `IUPACTableAgainstTheOracles`, by the integration oracle for R, Y, S, W, N, and by the frozen fixtures
   for Y, M, V, W, N.) At k = 3 the records `ACG` and `CGT` make the graph path
   `ACGT` and no record occurrence: the first oracle expects one context, the second none, and the engine must
   match both. Each oracle lists its expectations with the exact count per unit. It starts
   with the reviews' counterexamples: a record start whose dummy chain is pruned (`TACGAT` then `ACGA`, `AC`);
   the length-k lookup (`AAC` at k = 3); a motif with one plain and two marked matching edges; a sink dummy `AC$`
   inside a range and a source dummy leaving one; an interrupted mask scan with its bounds; a graph loaded
   without its mask (refused until the owner's decision #16 of 2026-10-08; answered since, with bounds, its lists
   the masked twin's: `PatternUnmasked`, §4.4); the DNA5 flank (`ACNTA`, `AC` at offset 0; written under `_DNA5_GRAPH`, run once
   by hand on a DNA5 build of the engine and its tests, not in CI: CI builds DNA and Protein only); a DNA4 island shorter than k
   (`AAAAANACNCCCCC`, `AC`: not covered, not claimed); zero suffix anchors with a prefix context (`ACGTT`, `AC`);
   the cross-record path (`ACG`, `CGT` at coordinates 0, 1); the double occurrence in one k-mer (`ACGAC`, `AC`);
   the record-local offset (`ACGTA`, `CCCCC`, `TA`); a repeated k-mer within one record (one context, several
   occurrences); two equal-length records with hits at equal local coordinates (`seq_id` tells them apart); a
   palindrome (`ACGT`) counted once; a `suffix` query on a wrapped PRIMARY graph (record `ACGA`, `TCGT` stored,
   `CGA` asked, the suffix of the virtual k-mer only: refused for `suffix`, found by `any_offset` at the virtual
   node, offset k − L); a reverse hit of a pattern longer than k on a wrapped PRIMARY
   graph (`AACG` at k = 3, anchored through rc(`AAC`)); an anchor without a completed path (`AAAC` with `AAA`
   present); two shards of 6,000 contexts against a threshold of 10,000 (withheld on the aggregate); a shard
   whose retrieval fails (withheld everywhere); an anchor with more labels than `max_labels_per_anchor`; a native
   CANONICAL graph (presence only); an unbudgeted backend (count only); a refused row; the deadline around the
   work and its answer (as built: a time stop before any work, the work time passing after the work — exact,
   `determinism: full` —, the deadline read in the middle of the assembly, and past the deadline: 503); a
   palindromic anchor window on a wrapped PRIMARY graph (`ACGTA` at k = 4: one anchor, one extension) and an IUPAC window with some palindromic
   instances; a pruned middle k-mer (`ACG`, `GTA` kept, `CGT` pruned at k = 3: `ACGTA` not found, both bases
   covered); a path verified in one label and only intersected in another (per-label support); a request with
   `require_support: record_verified` on an index without record mapping (refused); a sequence of fetches whose
   caches would exceed the share (one memory stop, not a silent overrun). Every later increment runs against
   it; the alignment golden gate of §11 is recorded in this step too.
1. **Count-only engine and route** (`mode: count`, `scope: suffix`, exact DNA, one strand, BASIC): mask handling
   settled first (the requirement, the charged scan, the bounds, `metagraph transform --mask-dummy` and
   `--pattern-build-mask`, since neither the mini index nor staging's graph has a mask file), then `PatternSearch` with the range DFS and the
   `W` rule with marked edges, the exact counts, `POST /pattern` on a single-graph server answering counts, the
   finalisation reserve, capabilities. The early milestone: a server that counts correctly before any annotation
   is read.
1b. **The edge mask's tooling and `/resolve`'s deadline** (after milestone 1's contract froze, 2026-10-07): the
   mask handling of increment 1 landed as its own milestone — `metagraph transform --mask-dummy`,
   `--pattern-build-mask` and `mask: built_at_load`, the `mask_required` message naming both remedies — with no
   field of the pattern contract changed; and `/resolve` accepts `bounds.time_budget_ms` with a finalisation
   reserve, as this route's deadline, because the search service runs `/resolve` as a job
   (`SPEC-labeled-traversal-core.md` §4.5).
2. **IUPAC, both strands, `any_offset`, graph modes, and the label-free extraction.** The `allowed` callback;
   the two oriented searches with the palindrome rule; the flank ranges with every real symbol and
   `ranges_visited`; native CANONICAL; wrapped PRIMARY with `suffix` refused and `any_offset` complete; and
   `enumerate()` behind `output.labels: none`: the contexts with k-mer, instance, offset, strand, node and row
   ids returned by the route under `max_contexts`, `all_or_count` and `partial`, with no annotation read —
   tested against the graph-walk oracle (every context, in BOSS edge order) on the tiny graphs and against the
   FASTA oracle on the mini index.
3. **Retrieval: discovery, placement, JSON, admission.** `LabelRecorder` discovery with the per-anchor cap, then
   `LabelQuery` placement, both with `DecodeBudget` and pacing; backend support and the unbudgeted opt-in;
   descriptors retained only while eligible; placement through `CoordToHeader` with `seq_id`; deduplication; the
   modes, `retrieval_complete` and the `withheld` reasons; the CLI subcommand. Integration test
   `integration_tests/test_pattern.py` on the mini index: placed occurrences equal a regex scan over
   `build/mini_refseq`'s FASTA (record, 1-based position, strand), counts per unit exact.
4. **`long`.** Opt-in per request (`long_search: "paths"`; without it the anchor-only answer of version 1 stays;
   the owner's decision of 2026-10-07, SPEC §12). The extension DFS, anchor-window bits, the two thresholds, per-label `support` with record bounds and `require_support`, the
   anchor mapping on wrapped PRIMARY graphs, the peptide state across the k boundary; tests with patterns
   spanning two and three k-mers, the cross-record path the bounds check must reject, the anchor without a path,
   and the reverse hit. **Served (2026-10-08, built in parallel with 5; SPEC §12.1, §17)**: the labels of each
   path returned with their support (`label_intersection`, `record_verified`) and `require_support:
   "record_verified"` filtering to the verified ones (the owner's decision #14); an addition to contract version
   1, every answer to a request without the option unchanged.
5. **Peptides.** The codon automaton and genetic code; tests with Leu/Ser/Arg-rich peptides and a peptide whose
   only graph instance crosses a stop codon. **Served (2026-10-08; SPEC §12.2, §17)**: the 20 residues and X, B,
   Z, J, every NCBI genetic code (1 the default; the owner's decision #15); peptides longer than k through the
   paths of increment 4. **The stop `*` served too (the owner's decisions #19 and #21 of 2026-10-08; SPEC §18)**:
   a stop codon of the table; in tables 27, 28 and 31 it matches nothing, with the note `no_stop_codon`; tested
   against a six-frame translation oracle in tables 1, 2, 11, 27, 28 and 31 (`PatternPeptide.StopResidue`,
   `RandomStopPeptidesAgainstOracles`). (`4596bb3b` refused it, `stop_unsupported`, now retired.)
5c. **Masks optional (the owner's decisions #16 and #17 of 2026-10-08) — served (§4.4; SPEC §18)**: graphs
   without the dummy-edge mask answered with upper bounds and an estimate, lists exact, `mask_required` retired;
   the mask and the Bloom filter outside `index_fp`. Masked answers unchanged (checked against `4596bb3b`; the
   alignment golden gate IDENTICAL). Staging's mask (#18) follows once this is deployed: built in a staging-only
   directory and wired into `~/metagraph-refseq33m-experimental`, its `index_fp` unchanged.
5b. **Predicates** (§5.6). The predicate parser and evaluator over `LabelQuery` reads restricted to the
   predicate's labels; `output.labels: predicate_only`; the compute admission on the unfiltered count and the
   retrieval threshold on the selected count; path-level evaluation for `long`; unknown labels per shard. Tests:
   the oracle suite extended with per-row and per-path cases (one row with A and another with B on a path:
   `any(A, B)` selects no path), `at_least` on cohorts, a typo in a label name reported as unknown, a predicate
   pass stopped by its budget (`at_least`, withheld), and a selective filter on a pattern whose unfiltered count
   exceeds `max_contexts`. **Built for L ≤ k (2026-10-08, SPEC §19)**: the language (`pattern_predicate.cpp`), the
   selection pass (`pattern_selection.cpp`, strand-consistent `"either"`, whole-row units) and the route
   (`pattern.cpp`); the path-level part waits for the supported-path search (5s) and runs on it (5b-5).
5s. **Supported paths** (P23–P29). The engine's extension with a support tracker and a path sink
   (`pattern_search.cpp`); the label and record-level trackers over the request's labelled retrieval, the
   anchors' rows read whole as the permitted set, a row cache, strand-consistent support
   (`pattern_support.cpp`); the route's `long_search: "supported_paths"` (`pattern.cpp`); and 5b-5, a predicate's
   selection of the supported paths (`pattern_supported.cpp`: decided at completion, monotone pruning, `"either"`
   with the mirror walks found or read). **Built in round C (2026-10-09, SPEC §20)** with the motif-level predicate
   of L ≤ k (SPEC §25), tested against a walk oracle over the records, a recursive evaluator of the request's
   predicate and the labels of `long_search: "paths"` (`PatternSupportedRoute.*`, `PatternSupport.*`,
   `PatternMotif.*`).
6. **Multi-graph, barriers, per-shard budgets, benchmark.** The shared fan-out helper with the sorted merge,
   shard identity, per-shard budget shares and the barriers; resident-only shards; a benchmark on
   refseq33m-experimental (16-, 20-, 25-nt motifs in both scopes, a 29-nt IUPAC promoter, three peptides)
   recording `ranges_visited`, mask scans, counts per unit, rows, elapsed and memory, and an estimate for the
   Logan hosts from their k-mer counts.
7. **Service.** The search service's job kind (its repository); docs (`docs/source`); the contract of §7
   promoted into a SPEC section once frozen. The Python client methods are deferred until a user needs them (§9).

Estimate: increments 0–2 about 2 weeks for one engineer who knows the BOSS and traversal code (the count-only
milestone, mask handling included); 3 about 2 weeks (the two-step retrieval is integration work, not reuse);
4 about a week; 5 about a week; 5b about a week (the parser is small; the two-step support check for `long` and
the three admissions are the work); 6–7 about a week and a half, of which the service side is half. The known unknown
is the benchmark in 6: a repeat-rich motif in `any_offset` scope on a Logan host sets the real defaults for
`max_steps`, `max_contexts` and the information floor.

## 14. Open decisions for the owner

Decided on 2026-10-07 (the owner, with the service session's review): the route is `POST /pattern`; the service
serves it **as a job**, not as a synchronous tool, with the job's defaults `mode: count` and `scope: any_offset`
and `partial` exposed; every host carrying the block offers it, canonical indexes included; the service session
builds the job, org-id integrates, this side delivers the frozen contract, the fixture bodies and the `pattern`
block on both capabilities routes with milestone 1; the route's default `time_budget_ms` is 60 s and the cap
600 s (§5.3).

Decided on 2026-10-08 (the owner's decisions #12–#21, `OWNER_DECISIONS_2026-10-07.md`): semantic compatibility
for contract version 1 (#12); paths opt-in through `long_search: "paths"` (#13); the labels of paths returned, each
with its support, `require_support: "record_verified"` filtering (#14); peptides with X, B, Z, J and every NCBI
table (#15); **masks optional**: a graph without its mask served with upper bounds and an estimate, f sampled at
load (#16, §4.4); **the mask and the Bloom filter derived data**, outside `index_fp` (#17); staging's mask built
once this is deployed (#18, the owner's step); **the stop `*` served** as a stop codon of the table, a note where
the table has none (#19); `genetic_code_unknown` (#20); the context stops of tables 27, 28 and 31 as their residue,
never as `*` (#21). Still open:

1. Whether to keep `/pattern` or rename to `/motif` (this draft keeps `/pattern`).
2. The route's own default mode (`all_or_count`, this draft) and default scope (`any_offset` for L ≤ k; `suffix`
   is the cheaper alternative and always states `absence_scope`).
3. Defaults: the information floor (24 bits, discovery only), `max_contexts` (10,000), `max_anchors` and
   `max_paths` (1,000), `max_labels_per_anchor` (64), `max_steps` (10⁸ per shard, §5.3), `max_memory_mb` (256 per
   request), the finalisation reserve (250 ms); all server policy.
4. *(Decided, #14: returned, each marked with its `support`.)* Whether `label_intersection` labels are returned
   for `long` without `require_support`, or only counted.
5. *(Decided 2026-10-07, above: every host carrying the block, canonical indexes included.)* Whether canonical
   indexes get the route at all in v1 (presence and contexts only), or answer `pattern_unsupported` until
   placement exists.
6. Whether `allow_unbudgeted_annotation` exists at all, or unbudgeted backends get `count` only.
7. The CLI subcommand: keep (tests and offline use) or server-only.
8. The service's job kind and tool names (the service's call, in the `pattern_` family).

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
- A pattern's N never matches a record's N symbol; flanks do (DNA5 builds only, untested: below).
- The engine needs the succinct graph representation; other representations answer 400 with the reason. Its
  edge mask (the `.edgemask` file, or built at load with `--pattern-build-mask`) makes the counts exact; without
  it (the owner's decision #16) counts are `exact` where provable and otherwise `bounds` with an estimate that is
  not a bound, the thresholds admit on the upper bound (note `threshold_upper_bound`), and the lists stay exact
  (§4.4).
- Budgeted retrieval exists only on the row-diff family with the budgeted decode; elsewhere `count` only, or
  unbudgeted by explicit opt-in; an anchor with more labels than the per-anchor cap withholds the results.
- Indexes are served resident only; a shard not resident is a stopped shard.
- Caps and budgets are per shard on multi-graph servers (`per_shard`); a column in two shards is two results;
  admission is on the aggregate and publication needs every shard.
- A deadline stop is not deterministic (`determinism: time_limited`); a finalisation overrun is a 503.
- Support is per label: `record_verified` only on a BASIC index with coordinates and record mapping; elsewhere
  `label_intersection`, and `none(A)` means only that A does not annotate every k-mer of the path.
- Predicates name annotation columns only; records inside a column are not addressable in v1.

Limitations of the build, found or confirmed by the review of 2026-10-07 and stated in the SPEC:
- **DNA5** (`$ACGTN`) graphs are refused by the route and the CLI (`alphabet_untested`, capabilities
  `available: false`; the owner's decision of 2026-10-07) until a DNA5 build passes the pattern tests; the engine
  itself still supports them. No DNA5 server or CLI has been built or run, and CI builds DNA and Protein only. The oracle suite is DNA5-aware; its DNA5 branches ran once, by hand,
  on a DNA5 build of the engine and its two unit-test files (the integration of the review's fixes, 2026-10-07:
  79 of 79 PatternSearch and PatternSearchFixes tests passed, after two test expectations were corrected to the
  stated behaviour). On an odd-k wrapped PRIMARY DNA5 graph a k-mer with N at its centre between complementary
  flanks (`ACNGT`) equals its reverse complement; `CanonicalDBG` exposes it at two node ids, and the engine counts
  and releases it twice (that run confirms it: 2 + 2 for `AC` in `ACNGT` at k = 5). A leading N run is
  searched there, not skipped (§7.8 of the SPEC). These are the items a DNA5 build must settle before the
  refusal is lifted.
- **Without the mask** (§4.4; the owner's decision #16 of 2026-10-08): the estimate assumes source dummies as
  frequent among a count's candidates as in the whole graph (false for a pattern at a record start); `partial`
  keeps every unchecked range it discovers (24 bytes each, at most one per step) rather than about
  `max_contexts`; on an even-k wrapped PRIMARY graph `stop_at_threshold` can fire late or not at all. The mask and
  the Bloom filter are fingerprinted nowhere (decision #17): a foreign Bloom filter of the same k and mode can hide
  k-mers under an unchanged `index_fp`.
- **The stop `*`** (§6; decisions #19, #21): a long peptide with `*` in tables 27, 28 and 31 has no path but is
  still searched for its anchors and gated by the floor on its anchor windows; the `bad_alphabet` message of a
  protein pattern does not name `*` (kept byte for byte).
- **A mask with a valid W = `$` edge is refused** (§4; the owner's decision of 2026-10-07): `extend` on a masked
  graph used to write a mask with valid dummy edges, on which counts could be overstated, `exact` included. Such
  a mask is found once at load and answered `mask_invalid` (capabilities `available: false`) until the graph is
  masked again (`transform --mask-dummy --force`) and the server restarted; `extend` now rebuilds the mask of a
  masked graph. Beyond that check the mask is trusted as written.
- **N runs inside a pattern** cost about min(4^run, edges / 4^a) ranges per level whatever the bits (§5.3); only
  a run at the start of a searched window is skipped, and only on `$ACGT`. Internal anchors (§12) would remove
  the rest.
- **O(L) per pattern** (parsing, its bits, its palindrome test) is inside the deadline but charged no step and
  interrupted by no clock reading: many very long patterns can make a request's 503 come after its
  `time_budget_ms`. The low-complexity note reads the clock between its pieces since the review GPT-3 (§5.3).
- **The time kept back for the answer** grows with what the answer holds, at configured delivery rates (10 and
  50 MB/s, `/traverse`'s starting values), not rates measured on the host: on a fast host a request with many
  results stops its work earlier than it needed to (on an M-series Mac, a 2 s request of 15 patterns with 19,283
  contexts each stopped its work after about 0.3 s, and writing what it held took about 0.13 s). Measuring the
  rates on the server, as `/traverse`'s attempts do, is open.
- **What CI does not check**: CI runs the unit suites, the integration tests that build their own graphs
  (`TestPatternSynthetic`) and the fixture bodies against the SPEC (`TestPatternFixtureBodies`). The fixtures
  against the binary, the byte identity of `/search`, `/align` and an unbudgeted `/resolve` against the base
  build, the oracle comparisons on the mini index and the alignment golden gate (§11) run only on a developer
  machine (`$METAGRAPH_REQUIRE_GUARDS=1` makes a run without their inputs fail instead of skip); a DNA5 build
  runs in no CI job (it ran once by hand, above).

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
| 4 | the mask may be absent; the invalid-edge scan is work; an interrupted subtraction is an upper bound | §4: the mask is required (`mask_required`; optional since the owner's decision #16 of 2026-10-08, §4.4); §4.1: the scan charged to `max_steps` and deadline-checked; `bounds` relation with `lower` and `upper`; §5.3 wording corrected |
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

## 22. The owner's addition (2026-10-07): the label-free path

"There must be a fast query without touching the labels. Labels slow down a lot; preserve the functionality
without labels (count and extract rows/k-mers)." Folded in as goal 9 and `output.labels: none` (§4.3, §5.2,
§7): counting and extracting contexts — k-mer, instance, offset, strand, node and row ids, paths — with no
annotation read, on every backend and graph mode, under the graph-side budgets only; the engine's `enumerate()`
is this path and the annotation steps sit on top of it. It belongs to milestone 1 (increment 2), and a later
increment adds "labels for given rows" (§12) so that annotation is paid for only on chosen rows.

## 23. The owner's decision (2026-10-07): no Python client yet

The Python client methods (goal 7, §9, §10, increment 7) are dropped until a user needs them: nothing in the
stack calls them, and unused code is weight to maintain. The route, the service's job kind and the SPEC with its
fixtures are the serving layers.
