# Design note: a label-constrained graph-traversal endpoint (`/traverse`)

**Status:** proposal / design sketch (uncommitted draft; APIs below are proposed)
**Date:** 2026-10-02
**Context:** repo `ratschlab/metagraph` @ `1de02cec` (branch `master`), with orchestration in
`metagraph-search-service`. Source references describe the inspected checkout, not verified
capabilities of every deployed index.

**Central interaction:** an agent defines a **parametrized traversal strategy**, applies it to a
**set of graph locations, often resolved from search hits**, inspects an **aggregated result**, and
**iteratively refines the strategy**.
Individual graph walks are backend work units, not the primary agent interface. The new async job
type is **traversal**: apply a strategy to explicit graph locations and return sequences with label
support along their paths, a compact subgraph, or both. Small returned graphs are reusable analysis
objects with their own MCP inspection and extraction tools.

The objective is to recover and compare sequence context outside and through the original query:
flanks, branches, neighbouring features, within-query support and candidate sequence variants.
Target investigations include
mobile-element context around *bla*<sub>NDM</sub>, plasmid neighbourhoods of *mcr-1*, adjacency of
the caffeine-pathway genes ndmA–D, allele enumeration, and TMPRSS2–ERG junction interpretation.
These are proposed applications, not promised outcomes: a graph path, shared sample label, or pair
of recovered arms does not by itself establish
physical linkage, plasmid membership, an operon, or a fusion.

---

## 0. First working release: a small end-to-end implementation

The first release should complete the agent loop on **one verified ATB shard**, using explicit
sequence-based seeds and `within_label`. Begin with deterministic bidirectional extension that stops
at the first permissible branch, the per-path edge-reuse limit or the requested bounds. Save the
sequences and boundary evidence, group exact flanks, and expose the operation through one async job
and an MCP workflow such as `extend_hits`. Include a reproducible pilot subset through the same operation.
Do not require a durable node-ID registry, search-retention integration, general subgraph editing,
sampling, source tracing, or resumable frontiers before this works and is measured.

The minimum provenance is the exact input list/digest, normalized strategy, fixed graph/annotation
release, per-item outcomes and saved result. The larger contracts below describe the roadmap;
API v1 advertises only implemented capabilities. The source-release identifier can be configured
from the deployed artifact manifest; cross-release anchor migration is not required.

ATB is the proposed first target because its service metadata describes contig-level labels;
confirm the actual shard, annotation format and source-contig availability before running the pilot.
Validate candidate paths against source contigs, including repeats and negative controls. Logan's
run-level labels are a later, harder case: many organisms/contigs may share a label. Neither linear
walks on ATB nor rapid branching on Logan are measured guarantees. Do not assume strict
`sequence_trace` works on the large CANONICAL/PRIMARY deployments; the current annotation buffer
disables it there. Direct source-record retrieval/validation is a separate option.

Priority additions are **inside-query support diagnostics (§4.6), compressed flank tries (§5.1.1),
and annotation through existing search tools (§4.4)**, assessed with capability-dependent benchmark
tasks (§9.4). These can ship independently as their inputs become available; they do not all require
new core traversal. Metadata summaries (§5.1.2) complement them. The most valuable next traversal
objective remains two-hit reachability (§4.5), supported by measured annotation access and pilots.

Generalize to three or more hit groups only for a biological question that requires them and after
measuring state growth. Literal motif stops and richer graph inspection follow selected needs.
Contrastive context analysis (§6.2) and bounded bubble analysis (§4.7) are later directions; neither
is assumed to be a free output mode. Sampling, weighted strategies, coverage-drop heuristics and
resumption remain deferred.

## 1. What already exists (reuse, with explicit limits)

**Traversal primitives** — `src/graph/representation/base/sequence_graph.hpp`:

| Primitive | Use |
|---|---|
| Sequence-to-node mapping | Resolve the seed; preserve ordered, oriented traversal nodes rather than only canonical annotation keys |
| `adjacent_outgoing_nodes` / `adjacent_incoming_nodes` | Expand the frontier in both directions |
| `traverse` / `traverse_back` | Single-character steps |
| `call_outgoing_kmers` / `call_incoming_kmers` | Enumerate neighbours with the appended/prepended character |
| `outdegree` / `indegree` | Structural branch detection; permissible degree must be recomputed after label filtering |
| `get_node_sequence(node)` | Materialise one node's k-mer; spell longer paths with validated overlaps or extension characters |

`src/cli/query.cpp:503` already implements neighbourhood-like query-hull expansion with depth and
fork limits. It uses **DFS**, and its backward expansion has limitations for noncanonical graphs
(`:795–820`). Reuse the ideas, rather than treating this helper as a complete endpoint algorithm.

**Label-consistent extension** — `src/graph/alignment/aligner_labeled.cpp:175–300`:
`LabeledExtender::call_outgoing` propagates label intersections and checks coordinate offsets.
Its unitig handling includes deferred checks corrected by `flush()`, and it operates within an
alignment engine. Extract reusable support operations; implement and test a traversal-specific
state machine instead of copying the method as a standalone BFS filter.

**Existing positional support** — `src/cli/query.hpp:81` defines per-label `LabelSigVec` bit vectors.
Core `with_signature` serializes `signature.presence_mask` and a score (`query.cpp:221–234`), and the
async client already accepts `with_signature` (`app/metagraph/client.py:197`). Exposing useful query
support intervals may therefore need service/MCP adaptation rather than a new graph algorithm.
Verify deployed availability, mask encoding and retention through result storage/compaction. With
`align=true`, the core replaces the input with the best-alignment sequence before computing the
signature (`query.cpp:1197–1228`); original and transformed coordinate frames must remain distinct.

**Annotation access** — `src/graph/alignment/annotation_buffer.hpp`:

- `queue_path` + `fetch_queued_annotations` supports batched annotation reads.
- Cached/interned column sets reduce repeated label decoding and set storage.
- `get_labels` / `get_labels_and_coords` provides node support; mapping oriented nodes to
  annotation rows must respect canonical wrappers.
- Coordinate availability depends on graph mode as well as annotation type. The current buffer
  explicitly disables coordinates for CANONICAL/PRIMARY graphs (`annotation_buffer.cpp:27–30`).

Use a **backend-aware, batched label-membership predicate**, with encoded permitted columns, rather
than assuming every step requires all labels. The optional `GetEntrySupport::get(row, column)`
interface exists (`annotation/binary_matrix/base/binary_matrix.hpp:104`). ColumnMajor provides a
compressed-bit lookup (`column_sparse/column_major.cpp:20`); BRWT descends the relevant tree branch
(`multi_brwt/brwt.cpp:25`). These are promising for small permitted-label sets.

RowDiff does not implement that optional interface (`row_diff/row_diff.hpp:97`). Start with cached,
batched row reconstruction there, and compare against a selected-column implementation only if
measurements justify adding one. **Never use whole-column extraction as an online fallback:**
`RowDiff::get_column` currently scans all graph nodes and reconstructs rows (`:150–163`). Even
`path_common` can use targeted membership if its allowed-label universe is explicitly bounded;
full rows are needed when discovering or reporting unrestricted support. In v1, report the checked
permitted label(s), not an implied enumeration of every label carried by a returned segment.

Batch by frontier or bounded chunks; measure direct membership versus row reuse on real layouts,
varying the permitted-label count and cold/warm state. Annotation time, reconstructed rows/bytes and
memory matter more than the accessor's name. The aligner proves reusable machinery, not equal
runtime: its query/scoring constraints can prune work that unrestricted extension must explore.

**Compaction/export** — `src/cli/assemble.cpp:220–265` provides a unitig/GFA precedent.
Avoid running whole-graph extraction per request: differential assembly allocates graph-sized
masks/counts (`annotated_graph_algorithm.cpp:345–355`), and masked-graph unitig extraction can also
materialise a global mask. Build and compact only the bounded local result.

---

## 2. Backend gap and service boundary

The new core operation combines an explicit seed, user-selected annotation labels, traversal
policy, bounds, and a compact result. The existing labeled aligner is not exposed as this operation.

**Follow `/search`'s annotated-graph loading/routing pattern when adding the core route:**

| Existing route | Graph access | Multi-graph servers |
|---|---|---|
| `/search` (`server.cpp:379`) | `AnnotatedDBG` | Supported through `graphs` and `filter_graphs_from_list` |
| `/align` (`server.cpp:503`) | Plain graph | Rejects multi-graph requests |

Reuse `/search`'s routing machinery, but preserve **exact physical shard identity** in every work
item and result. A known anchor should not trigger a fresh search across its entire taxonomy scope.
Never join paths across shards merely because their node IDs or sequence strings coincide.

There are two API layers:

1. **C++ backend `POST /traverse`:** execute a bounded batch of explicit anchors on a named physical
   graph/index generation. Return compact per-anchor results and counters. No `search_id` lookup,
   Redis job lifecycle, cohort aggregation, or public polling endpoints belong here.
2. **Async service `/traversals`:** resolve/freeze a hit selection, group anchors by shard, schedule
   core requests, persist artifacts, reduce results, expose progress and refinement through MCP.

This keeps per-node traversal inside C++; the async service must not implement a graph walk through
thousands of neighbour HTTP calls.

Before graph work, use a **source-record extraction fast path** when an authoritative record and
occurrence are available and the request is simply for its observed flanks. Slice the requested
window and retain source/version/coordinate provenance. Use graph exploration for alternative or
unresolved context and topology questions. Report the execution/evidence mode explicitly; source
extraction preserves observed repeat copies, while graph unrolling follows the declared cycle policy.
Source-record boundaries and graph-imposed caps are different stopping reasons. Graph traversal
cannot recover missing sequence or establish contig joins without supporting evidence.

---

## 3. Agent contract: strategy × frozen hit set → results → refinement

### 3.1 Strategy as a versioned value

A strategy is declarative data, validated against index capabilities. It is not arbitrary code
executed on the server. Normalize defaults and store both the requested and effective strategy;
reject unsupported policies rather than silently weakening evidence requirements.

The example describes the extensible contract, not a commitment to every field in v1. Start with
inline strategies normalized and digested in the job; a separate strategy registry can wait.
Sampling, advanced support modes and resumption are omitted from the first accepted schema and
rejected if requested. Illustrative values below are **not benchmarked production defaults**:

```json
{
  "schema_version": 1,
  "direction": "both",
  "labels": {"mode": "within_label", "permit": "per_hit"},
  "exploration": {"method": "bounded_bfs"},
  "bounds": {
    "max_extension_bp_per_direction": 5000,
    "max_states": 200000,
    "max_nodes": 100000,
    "max_edges": 200000,
    "max_frontier_states": 10000,
    "max_branch_depth": 8,
    "max_output_bp": 200000,
    "max_output_bytes": 2000000,
    "time_budget_ms": 30000
  },
  "stop_when": [],
  "emit": {"sequences": "handles_and_windows", "labels": true, "coordinates": "if_available"}
}
```

The normalized strategy receives a `strategy_id`/digest. A refinement creates a new immutable
version with `parent_strategy_id`; it never changes the strategy attached to a completed job.
Maintain total budgets across the refinement campaign as well as per-anchor and per-job bounds.

### 3.2 Freeze the input population

The same operation accepts one anchor as a batch of size one for focused follow-up. **For the first
release**, accept an explicit list of seeds/labels on the configured shard. Persist the exact list,
its digest, stable per-item identities, release and normalization/selection rules with the job.
This is the minimal frozen population; a search-derived manifest service is not a prerequisite.
Any known source-search cap/scope remains attached; an unknown eligible population stays unknown.

Later the agent can select hits by `search_id` plus filters without copying sequences or graph IDs.
The async service resolves that expression into an immutable `hit_set_id` manifest containing:

- Search/result snapshot and original search scope, retention limits, truncation and completeness.
- Filter, deterministic sort/tie-break, limit, deduplication unit, and selection-sampling policy/seed.
- Exact selected result identities, query identities, immutable shard provenance and resolved anchors.
- Index/annotation generations and a content digest of the selection.
- Counts of eligible rows, selected rows, distinct samples, resolved occurrences, and unresolved hits.
  An unknown eligible total stays unknown; selection from capped search results is not an archive census.

Specify whether selection waits for a completed search or explicitly snapshots a partial result.
Repeated strategy runs use the same `hit_set_id` by default. A changed hit population creates a new
manifest and is identified separately from a strategy change. Archive the manifest so refinement
does not depend on the original search's retention lifetime.

One search row may resolve to multiple disconnected matches or source occurrences. Expose an
explicit selection policy (all within a cap, user-selected, ranked, or sampled), with excluded and
unresolved counts. Do not silently pick the first occurrence or treat rows as unique carriers.

### 3.3 Sequence-based seeds first; durable anchors later

The minimal seed is `(physical_graph_id, sequence_or_handle, query_interval, orientation, label)`,
resolved against the release pinned for the job. An exact oriented k-mer is the smallest seed; a
longer interval must map as a contiguous path with the required label on every seed node. Do not
use the first and last supported k-mers across a partially matching query as one supported interval.
Reject the interval or report disconnected supported blocks with an explicit selection policy.
A k-mer/path identifies graph sequence, not a unique physical source occurrence.

Re-map these sequence seeds inside the core; no persisted node-path registry is necessary for v1.
Extend the selected interval's left/right endpoints. Record the actual graph, annotation and source-
mapping release: a seed may still map after an update while its neighbourhood or support changes.
Successful re-mapping therefore does not detect an unchanged index. Freeze the loaded release for
the job; reject requests explicitly pinned to unavailable releases. A configured release manifest
is sufficient initially; an endpoint URL or node count is not a release identifier.

Later an opaque `anchor_id` can bind the resolved oriented path and optional source occurrence to
that physical shard/release, avoiding repeated resolution. Bare node IDs without release and
orientation are never accepted. There is no initial requirement to migrate IDs between releases.

Integration gaps for search-derived selections and richer anchors, established by source inspection:

- Core search JSON returns labels/counts/signatures/optional coordinates, not reusable oriented node
  paths (`src/cli/query.hpp:68–107`, `query.cpp:188`). `/align` also omits its node path.
- Async MCP `_compact_row` omits even the result ID (`app/mcp_server.py:2408`).
- Async aggregation rewrites child `Result.task_id` to the parent
  (`app/results/staging_merge.py:163`); retain source-shard provenance in a separate immutable field.
- Displayed coordinates can be subsampled. They are not an authoritative anchor representation.
- Current backend-cache `/stats` identity is not necessarily a graph-build identity
  (`app/backend_cache/fingerprint.py:340`). Add a build/generation identifier covering graph,
  annotations and source-record mapping; do not substitute an endpoint address or node count.

Mint richer anchors lazily when needed. Sequence-based inputs remove the need to export node paths
or retain original result IDs in the first explicit-list workflow, but do not remove exact shard
routing or input/release provenance. Unknown source shard identity must be resolved explicitly;
matching sequence strings alone do not authorize joining graphs.

### 3.4 Public async API and core execution

The first public request uses explicit seeds and an inline strategy as described above. The following
is the later search-selection form of the same operation; it is not required before the first pilot:

```json
{
  "job_type": "traversal",
  "hits": {
    "search_id": "search_...",
    "filter": {"normalized_score": {"gt": 0.9}},
    "sort": ["normalized_score:desc", "result_id:asc"],
    "limit": 500,
    "snapshot": "require_complete"
  },
  "strategy_id": "strategy_...",
  "result_view": "sequences_and_subgraph",
  "job_bounds": {"max_total_states": 20000000, "max_wall_time_ms": 600000}
}
```

Filter syntax is schematic; use the async service's validated grammar in the final implementation.
Choose exactly one input form: a `hits` selection, an existing `hit_set_id`, or explicit `anchors`
for direct exploration, including jobs without a prior search. Explicit anchors carry the graph,
generation and location/sequence information from §3.3; resolution saves the explicit-list snapshot
initially, or the richer manifest when supported, and reports ambiguity. Accept either an inline
`strategy` object (normalized and assigned an ID on submission) or a previously registered `strategy_id`, never both.
`result_view` can select `sequences`, `subgraph`, or `sequences_and_subgraph`. A sequence view is a
bounded set of supported extracted paths, not necessarily all walks through the subgraph; its path
selection method and completeness are explicit. Keep the authoritative local graph/evidence artifact
when available so a later view does not require exploring the large index again.

```text
POST /traversals                    -> job_id, hit_set_id, strategy_id, capability/selection report
GET  /traversals/{job_id}            -> progress, per-shard status, observed cost
GET  /traversals/{job_id}/results    -> bounded aggregate and artifact references
GET  /traversals/{job_id}/details    -> paged selected anchors, segments, paths or groups
POST /traversals/{job_id}/summarize  -> re-group stored evidence; no graph traversal
POST /traversals/{job_id}/cancel     -> cancel outstanding work
```

MCP tools mirror this workflow: submit a strategy over a hit set, poll, inspect aggregates/details,
and submit a refinement. The async worker sends core `/traverse` an explicit bounded batch:
`{graph_id, index_generation, anchors, effective_strategy, execution_budget}`. Deduplicate backend
work where possible while retaining the original selected-hit membership and counting unit.

## 4. Traversal semantics and stopping rules

### 4.1 Label policies and evidence

The API's traversal annotation labels are distinct from existing async `node_labels`, which route
searches through taxonomic partitions. Resolve taxonomy scope to graphs first; constrain walks with
raw annotation-column identities and report whether a label represents a sample, contig or other unit.

| Mode | Rule | Evidence obtained |
|---|---|---|
| `within_label` (MVP default) | Traverse separately for each selected anchor/label; require that label at every step | Every path k-mer occurs under that label |
| `path_common` | Propagate the intersection of permitted labels over the entire seed and each extension; stop when empty | At least one common permitted label supports all path k-mers |
| `any_of` | Each node has some permitted label; labels may differ between steps | Neighbourhood of the union of the permitted labels; candidate paths may switch labels |
| `all_of` | Every node, including the seed, carries every permitted label | All permitted labels support every path k-mer |
| `sequence_trace` | Preserve a specific source record/occurrence and successive coordinates throughout extension | Path follows an indexed source-sequence occurrence |

No unrestricted default is exposed. Empty permit sets are invalid; `per_hit` explicitly means the
resolved hit's annotation label, not all graph labels. `any_of` is an explicit exploratory mode and
must not be reported as a path present in one sample. With multiple permits, `within_label` produces
separate work items rather than joining labels.

**Even `within_label` is not proof of a contiguous molecule.** A label can contain multiple reads,
contigs, organisms or repeat occurrences. Echo both the mode and `support_kind` in every result.
A path through the permitted graph is a candidate sequence unless stronger evidence is available.

For `sequence_trace`, carry surviving `(label, source_record, seed_occurrence, strand)` states and
require consecutive coordinates in that same record. Use `.seqs` / `CoordToHeader` boundaries.
Existing coordinate alignment alone is weaker: `tests/annotation/test_aligner_labeled.cpp:319–337`
explicitly accepts a cross-record path. Support the strict mode only on graph/index combinations
whose coordinate and orientation semantics have been validated; reject unsupported requests.
Joining left and right arms also requires compatible labels/occurrences, not just a shared seed.

### 4.2 BFS state, cycles and strand

For union-filtered, radius-only neighbourhood reachability, visited oriented nodes plus
direction/distance may suffice. Branch-depth limits and motif stops introduce further path/matcher
state even in that mode. For `path_common` or trace exploration, state also includes surviving
support/occurrences.
Two arrivals at one node can have different legal continuations. Do not use a single `visited[node]`
bit or union their supports in a way that invents a recombinant path. Intern support sets; apply
state dominance only when remaining budgets and support semantics make it valid. Retain predecessor
compatibility if returning paths. A bounded graph is not an enumeration of all its possible walks.

Unique-node exploration records cycles as edges. For extracted candidate extensions, the proposed
initial cycle policy is **at most one traversal of each directed edge per candidate path**. This is
a conservative implementation of "take each loop at most once": a simple cycle can be traversed once
and left through an unused exit, while overlapping cycles that share edges may be restricted further.
Apply the policy per path, not globally across candidates; retain cycles in the compact graph.
The supplied seed is retained as matched input rather than rewritten to satisfy an extension
limit. Combined extension arms and multi-hit paths must enforce the same edge-use budget across
their newly traversed portions, rather than resetting it at every intermediate hit. Intermediate
hit intervals count as newly traversed sequence; they do not become additional exempt seeds.
The edge cap bounds graph exploration: direct source-record extraction preserves observed repeat
copies and reports that this cap was not applied, while retaining span/gap and output bounds.
Report skipped cycle continuations and unresolved repeat multiplicity; one unrolling is not an
estimate of biological copy number. Length/state/output bounds remain necessary for branching.
Sampled walks with other revisit policies are later capabilities with explicit length/visit limits.
Label filtering can create a permissible unbranched run inside a globally branching graph;
branch-based stops and limits use the **permissible** outgoing degree.

Direction is relative to the seed's mapped orientation. Preserve oriented traversal IDs; map to
base IDs only for annotation lookup. Left extension prepends sequence correctly, and reverse-
complement seeds must produce equivalent contexts after orientation normalization. Protein indexes
do not receive DNA reverse-complement semantics. Graph-relative left/right is not automatically
biological upstream/downstream without gene-strand annotation.

### 4.3 Bounds and sampling

Require per-direction extension bounds plus limits on states, nodes, edges, frontier, annotation
payload/coordinate occurrences, output bases/bytes, time and cancellation. Count bases **outside
the seed**, with overlaps counted once. Schedule both directions fairly so one does not consume all
resources first. Job-level budgets additionally limit the entire hit-set operation; use explicit
per-anchor quotas/scheduling and report anchors not started because the global budget was exhausted.

Use `bounded_bfs` first, with stable neighbour ordering and deterministic work-item ordering.
Later `sample_paths` requires `n_paths`, maximum path length, RNG seed,
cycle/duplicate policy and a named transition distribution. A random choice at each branch is not uniform
sampling of complete paths. Derive per-anchor seeds from stable identities so retries/scheduling
do not change the requested sampling policy. Abundance weighting requires appropriate per-label
counts and a defined normalization; existing search abundance parameters do not implement this
sampler; define zero/missing-count behaviour without silently substituting another distribution.
`label_diverse` similarly needs a precise diversity objective and weighting rule.

Distinguish:

- Reaching the requested radius or a declared semantic stop: the requested local domain may be
  complete, while the biological sequence may continue beyond it.
- Hitting a resource cap or cancellation: exploration of that requested domain is incomplete.
- Path sampling: nonexhaustive by design, even when all requested samples have been generated.
- Result pagination: delivery can be incomplete independently of exploration completeness.

Return stop counts/reasons per anchor/frontier and whether remaining-frontier counts are exact.
Never reinterpret a stopped/sampled/failed search as absence. Time-limited results can vary with
load; reproduce their exact content from the saved artifact, not from the strategy digest alone.

### 4.4 Semantic stops and feature annotation

Potential stops include first permissible branch, loss of required support, source-record boundary
in trace mode, and a supplied literal motif. Define OR/AND composition (start with OR), strand,
whether the seed itself is searched, and whether the stopping node/motif is included in the result.
Match motifs across segment boundaries on admissible paths, counting overlaps once. Carry matcher
state where needed; do not miss a motif because compaction or a branch split its sequence.

Recognizing an IS element, gene, or domain is a **feature-analysis operation**, not a free consequence
of BFS. Named feature stops/grouping require a versioned reference, matching method and thresholds.
Keep literal exact-motif stopping separate from biological annotation. A useful first workflow is
bounded traversal → annotate recovered sequence with established tools → aggregate → refine.
Adding a new feature analysis may reuse stored sequence without traversal, but must receive its own
analysis version and cost; it is not merely a cheap `group_by`.

**Early reuse: annotate recovered flanks by searching them.** Submit saved flank handles to existing
sequence search (for example, a suitable RefSeq index), retrieve candidate references, and verify
the matching interval against its annotation. A top hit may be a whole chromosome; its title does
not identify the recovered flank. Report match coverage/identity, coordinates and reference release,
and retain unresolved or ambiguous assignments. Where existing reference-search/annotation tools
suffice, no new backend feature detector is needed. Short or divergent features may require more
sensitive protein/profile analysis. This closes the extension → identify context → refine loop with
existing tools while keeping the annotation provenance separate from traversal support.

**Deferred idea: coverage-drop stopping.** On an index with suitable counts, a strategy could stop
at a specified relative drop. First define whether counts are per-sample read coverage, assembled-
sequence multiplicity or archive-wide support, and measure retrieval cost. Count access can require
annotation I/O; it is not automatically free. A drop is a heuristic, not proof of a chimera or contig
boundary, and a lack of a drop does not establish continuity. Keep this out of v1 pending calibration.

### 4.5 Joint paths through gene hits: candidate operons and gene neighbourhoods

A second traversal objective is **connect two or more hit groups on one bounded, consistently
supported path**. This supports questions such as: "Which selected loci contain genes A, B and C
close together, in what order and orientation, and with what intervening sequence?" It uses the same
batch-job, saved-artifact and refinement machinery as flank extension. Start with two groups after
measuring the extension primitive. Three or more groups wait for a biological question that needs
them and measured state-growth limits; they are not an automatic first-release generalization.

The two-group form can be exposed as `reach_between` or a `reaches_target` traversal objective:
return a bounded label/anchor table with `reached`, `not_reached_in_complete_scope` or `unknown`,
distance and a witness reference. Keep the witness sequence on the backend even when the initial
response is only a table. Target identity is a specified oriented hit interval, not an arbitrary
shared k-mer. Small output does not guarantee cheap exploration. Iterating the distance bound uses
the same refinement workflow; branch/cap stops cannot supply a negative reachability claim.
Examples include the supplied mcr-1/ISApl1 and blaNDM/ISAba125 hit groups, without inferring physical
adjacency or element identity merely from a graph connection.

Each group represents one requested gene/family and contains alternative mapped hits/occurrences.
Require one compatible occurrence from each requested group, not every hit returned by the gene
search. Homology discovery and protein-to-nucleotide mapping happen before this operation; the
traversal itself does not establish that a partial sequence hit is a complete gene.

The query specifies:

- Required hit groups, with an explicit ordered list or `order: any`.
- Gene-strand constraints, such as same strand, relative to the returned path; do not infer gene
  strand solely from the direction in which the graph was explored.
- Maximum total span and adjacent-hit gaps, with an explicit overlap policy. Thresholds are query
  choices, not universal operon definitions.
- Label/support mode, per-path cycle policy, maximum witnesses per locus and traversal/job budgets.

Return **one joint witness path for each reported arrangement**, including selected hit IDs and
source occurrences where known. Pairwise reachability is useful for pruning but is insufficient:
A–B and B–C can have different supporting labels, different B occurrences, or incompatible path
histories. Labels/occurrences, total span and edge-use restrictions must remain compatible across
all hits and the connecting sequence. A connected subgraph alone is not a witness for a joint path.

Each result includes the gene order, path-relative hit intervals and strands, total span, signed
adjacent gaps, common label support, evidence kind and sequence handles for the bridges/full path.
For half-open intervals ordered along the path, `gap_bp = next.start - previous.end`; negative values
indicate overlap. Call this an **intergenic gap only when full gene boundaries are established**;
otherwise it is a gap between mapped hit intervals. Identifying intervening genes requires separate
gene annotation; an unrequested stretch of sequence must not be assumed to be intergenic.
Report distances for the returned witnesses, and claim a minimum only when the search establishes
it under the declared constraints. Graph distance and source-coordinate distance remain distinct.

Use a source-coordinate fast path where available: group compatible occurrences by label/source
record, sort intervals and scan for qualifying arrangements before fetching their sequence. Otherwise,
restrict to candidate common labels/shards, seed from selective hit groups and explore the required
orientations within the span bound. Do not enumerate a Cartesian product of every gene's hits.
For two groups, unit-step BFS (or weighted search on compacted segments) can find connecting paths.
For several groups, traversal state additionally tracks ordered progress or a visited-group bitmask;
retain the support and path history needed by the chosen policies. Bound group count and states.

Report qualifying witnesses found, no qualifying path within completely evaluated **declared
constraints**, and unresolved cases separately. A cap, branch stop, missing source context or
nonexhaustive sampling cannot establish biological absence. Aggregate by arrangement and gap/span
profile, retaining member evidence and distinguishing labels, contigs and independent samples.
The agent can refine distances, order constraints or cohort selection and compare saved results.
A complementary discovery workflow starts from **one seed family**, annotates recurring neighbours
across a representative locus cohort, formulates candidate multi-gene sets, and tests those sets on
a broader held-out cohort. Compare conserved order, rearrangements, insertions and spacing; control
for duplicate/closely related loci before interpreting recurrence or enrichment.

This produces **candidate operon neighbourhoods**. Proximity, strand and conserved gene order are
useful prediction features; establishing an operon also concerns co-transcription and regulatory
organization. Add promoter/terminator, expression or curated transcription-unit evidence in a separate
analysis layer when available. See [operon prediction from genome features](https://pmc.ncbi.nlm.nih.gov/articles/PMC1800777/)
and [RegulonDB's operon/transcription-unit definitions](https://pmc.ncbi.nlm.nih.gov/articles/PMC6060643/).
A proposed MCP workflow such as `find_gene_neighbourhoods` can expose this objective with an
arrangement table and inspectable sequence evidence, without requiring generic subgraph editing.

### 4.6 Inside-query diagnostics: explain support before extending it

Given a query interval and selected labels, ask **where is the query supported, and what evidence
connects its supported parts?** This is a priority investigation mode for low-containment hits. The
query being diagnosed may match only partially; the contiguous-seed rule in §3.3 applies to anchors
used for extension/detours, not to the diagnostic query as a whole.

Implement three distinct evidence levels:

1. **Positional support profile:** reuse per-label k-mer signatures to return supported and unsupported
   runs, the longest supported run and missing-span locations. Include k, orientation, the exact
   queried sequence/digest and 0-based half-open k-mer-start intervals. Distinguish unavailable or
   unqueried data from measured missing support. If alignment transformed the sequence, expose that
   sequence and its mapping/CIGAR; do not present its mask as original-query support. This layer
   requires no BFS and does not establish physical continuity.
2. **Compatible path through the interval:** validate an exact query path or use bounded reference-
   guided alignment/anchored detours to explain unsupported spans. Return a label-compatible witness,
   its alignment to the original query, alternatives and stop reasons. Reuse aligner primitives,
   but do not equate unguided BFS with a query-explaining alignment. Reaching one endpoint alone is
   insufficient; the witness must satisfy the declared interval/matching constraints.
3. **Source-occurrence continuity:** where available, verify a specific record occurrence through the
   interval or junction. This is stronger evidence than every k-mer sharing one label or a candidate
   graph path, especially for mixtures, repeated genes and fusion-like queries.

A single substitution can disrupt multiple overlapping k-mers. Scattered gaps, one clean break, or two
supported arms can motivate hypotheses of divergence, missing sequence or a junction, but do not
automatically diagnose an allele, truncation, contig boundary or fusion. Prefer an evidence report
with targeted next steps. Distinguish a supported witness, complete failure under declared search
constraints, and unresolved/capped work. Query-only data may settle the question without traversal;
that is an intended outcome, measured explicitly in §9.4.

### 4.7 Later: bounded bubble inspection and allele evidence

For a selected biological case, inspect simple alternatives between agreed flanking anchors
**through the locus**, align/normalize the alternative sequences to a reference, and return branch
support with witness paths. Nested bubbles, indels, repeats and ambiguous anchors require additional
algorithmic work; this is not simply a different BFS emit format.

Per-bubble sample counts do not determine complete alleles: different variant sites may be supported
by different copies/strains in the same label. Named allele calls require a compatible path/source
occurrence across distinguishing sites and appropriate reference verification. Preserve unresolved
phase and multiple candidates rather than combining marginal branches. Graph-only names remain
candidate alleles; do not count them as observed phased alleles without source-continuity evidence.
Start with bounded local inspection; a general allele caller is not part of the first release. See the distinction between
partial matches and full-length allele identification in [NCBI's AMRFinderPlus methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC8208984/).

### 4.8 Sequence-space neighbours of a hit: one or two substitutions

**Question:** given the exact sequence of a hit, which closely related sequences also occur in the
index, and in which labels? Here "nearby" means sequence similarity, independently of genomic
proximity. This is a useful bounded search mode alongside context traversal; it does not require
walking neighbouring nodes in the de Bruijn graph. A similar k-mer alone does not identify the same
gene, locus or biological allele.

**Small first implementation — a single k-mer.** Resolve the hit to an explicit, oriented A/C/G/T
k-mer and freeze its graph/annotation release. For substitution-only distance, enumerate the
Hamming neighbourhood and perform exact graph-membership probes inside one backend batch:

- Exactly one substitution: `3k` candidates.
- Exactly two substitutions at distinct positions: `9 * k * (k - 1) / 2` candidates.
- For `k = 31`: **93 + 4,185 = 4,278 mutant candidates**, or **4,279 including the original**.
  The `31^2` intuition captures the position-pair scaling; unordered position pairs and the three
  alternative bases determine the actual count. Use the index's actual k; 31 is an example.

Generate position pairs once (`i < j`), with three non-reference bases at each position. Probe the
full requested set: a two-substitution neighbour can be present even when every corresponding
one-substitution intermediate is absent. A BFS restricted to present one-mutation neighbours
would miss such candidates. The original hit is a separate distance-zero result.

Use `DeBruijnGraph::kmer_to_node` (`src/graph/representation/base/sequence_graph.hpp`) for
representation-aware exact membership with a seed of length k. Filter by graph membership first,
then retrieve annotations only for present candidates, reusing the node/annotation mapping and
measured access strategy in §1. `AnnotationBuffer` queues paths: unrelated candidates must be queued
as singleton paths (then fetched together), or handled through a verified row-batch API, rather
than concatenated into an artificial path. Keep this inside one job/MCP call rather than thousands
of HTTP calls. For multiple seeds, cache/deduplicate
probes within the same release while preserving each seed-to-variant edit relation. Canonical or
reverse-complement-equivalent storage keys may share a lookup, but edits and distances must still
be reported in the declared seed orientation. Follow the selected index's strand contract.

**Proposed tool shape:** `find_sequence_neighbours(seed, max_substitutions=2, label_scope, bounds)`.
Allow either selected labels (e.g. the original hit's labels) or the full declared index scope, since
new variants may occur in labels where the original sequence is absent. Return a compact table of
variant sequences/handles, distance, 0-based substitution positions and reference/alternate bases,
label support with drill-down, release/scope and completeness. Distinguish distinct variant strings,
storage nodes, labels and source occurrences; a label count is not automatically abundance. Reject
ambiguous seed bases in the first mode rather than silently changing the enumeration semantics.

Cap seeds, candidate probes, annotation work and returned bytes separately. Report theoretical and
completed probe counts, present candidates, evaluated label scope, omitted result rows and stop
reasons. Exhaustive membership within the requested radius does not imply complete label reporting.
A negative result requires the requested membership and label scope to have been fully evaluated;
a capped/sampled result remains partial. Measure membership and annotation costs separately on the
deployed index before claiming this is inexpensive: the combinatorics alone are not a runtime test.

**Longer matches and indels are distinct modes.** For a length-L query, enumerating up to two
substitutions grows as `3L + 9L(L-1)/2`; many independent k-mer hits do not validate a whole sequence.
One practical follow-up is to extend a recovered near-neighbour seed and test its longer context
against the original hit. For a full-interval claim, require a compatible complete path and retain
the distinction from source-occurrence continuity. Indels require an edit-distance contract and
alignment/normalization; they are not included in the counts above. Existing scored/capped graph
alignment is a useful comparison, but must not be called an exhaustive edit-radius enumeration.

This operation can be implemented independently of general bubble calling (§4.7). Its immediate
value is finding local alternatives and comparing their label support; named allele identification
still requires the relevant discriminating sequence and source/reference evidence. Keep it optional
for the focused paper, and use it only if it strengthens the selected biological investigation.

---

## 5. Output: aggregate first, evidence available on demand

### 5.1 Tier 1 — reduction across the selected hit set

Report hit rows, anchors, distinct labels/samples and source occurrences separately. Declare the
counting unit and sample classification rule (for example, any supported positive occurrence versus
all evaluated occurrences). Deduplicate overlapping anchors and repeated records across shards
without discarding provenance. Partial/sampled data may establish a positive observation; lack of
a detection only supports absence within a completely evaluated requested scope.

The following **illustrative** numbers assume one selected anchor per sample. The ISCR1 group
requires the separately versioned detector; it is not inferred by the traversal primitive:

```json
{
  "hit_set_id": "hits_...",
  "strategy_id": "strategy_...",
  "analysis_id": "features_...",
  "counts": {
    "selected_rows": 412, "resolved_anchors": 412, "unique_samples": 412,
    "complete_anchors": 400, "partial_anchors": 11, "failed_anchors": 1
  },
  "feature_summary": {
    "feature": "ISCR1", "count_unit": "sample", "support_kind": "sample_kmer_support",
    "detected": 380, "not_detected_in_complete_scope": 20, "unknown": 12
  },
  "flank_groups": [
    {"group_id": "g1", "signature": "ISCR1", "n_samples": 380,
     "representative_anchor_id": "anchor_...", "representative_sequence_id": "seq_..."}
  ],
  "coordinate_coverage": {"anchors_with_source_coordinates": 0, "eligible_anchors": 412},
  "label_mode": "within_label",
  "incomplete_reasons": {"max_states": 11, "backend_failure": 1},
  "cost": {"states_expanded": 12000000, "annotation_rows_fetched": 91000,
           "elapsed_ms": 48210, "attempts": 413},
  "artifact_id": "graph_results_..."
}
```

Every group links to member anchors and supporting observations. MVP grouping can use exact
sequence hashes/lengths and topology; approximate similarity clusters require declared method,
threshold and orientation normalization. An observed representative is the default. A consensus
is a separately derived object that may not exist in any sample; never use it as evidence of an
observed path. Report prevalence only for the stated selection and evaluated scope, with original
search coverage and unknowns visible. A first agent-facing result can be "N of M evaluated labels
share this exact returned flank", with length, orientation and stop status. Different truncated
windows are not interchangeable full neighbourhoods. Pass distinct representative sequence handles
to existing comparison tools (such as `compare_references`) rather than re-reading duplicate flanks.

### 5.1.1 Compressed flank tries alongside whole-flank groups

Build a compressed prefix/radix trie from already recovered, bounded flank sequences. It exposes
shared prefix lengths and split positions: for example, a common tract followed by two context
branches. Keep exact whole-flank hashes as a complementary equality view; a trie does not merge
shared suffixes after an early substitution. This is an aggregation over saved sequences, not a
requirement to explore an unrestricted union-of-labels graph or materialize additional paths.

Compare windows from the same defined anchor boundary and orientation. Right flanks run outward
from the right boundary; left flanks are compared anchor-outward, with that display transform
recorded separately from the natural sequence-handle orientation. A compressed node stores its
sequence interval/handle, prefix length, children and member references. Report distinct samples,
occurrences and paths separately: one sample can support several branches, so child sample counts
need not sum to the parent's count.

Keep terminal reasons and denominators for observed, stopped, partial and unknown context at each
split. A short flank caused by a cap, branch stop or source-record end is not a divergent nucleotide
branch or evidence of downstream absence. Divergence positions nominate sequence-context boundaries;
calling them mobile-element boundaries or recombination hotspots requires additional evidence.
Tries are exact for the stored windows but are not guaranteed tiny: cap/page the display and retain
member drill-down, collapsed-branch counts and the authoritative sequences. Expose this as a summary
view through existing job/artifact operations rather than requiring another traversal endpoint.

### 5.1.2 Metadata-joined context summaries

Join saved flank/trie/arrangement membership to the existing hit metadata, retaining the metadata
snapshot, counting unit and join provenance. Summarize collection time, country, host or other
available fields, with missing/ambiguous values and denominators explicit. Keep collection dates
and their precision separate from submission/release dates; retain date ranges rather than treating
an interval start as an exact collection date. An earliest date is **earliest observed
among the selected dated records**, not the origin of a context or its first appearance worldwide.
Exact sequence groups, inferred features and biological source units remain separate keys. This
can be an async analysis over stored results with no new graph traversal; availability and metadata
quality still need checking. See [NCBI's metadata field definitions](https://www.ncbi.nlm.nih.gov/pathogens/pathogens_help/).

### 5.2 Tier 2 — compact graph and selected sequence views

Compact only the explored local graph into maximal permissible nonbranching segments, splitting
at support changes, anchors, trace boundaries and truncation points, or encoding those changes
explicitly. Return oriented edges with `k−1` overlaps and evidence-compatible path references.
Per-segment support alone does not authorize arbitrary concatenation across branches.

```json
{
  "anchor_id": "anchor_...", "graph_id": "graph_...", "index_generation": "generation_...",
  "k": 31, "direction_relative_to": "seed_orientation",
  "segments": [
    {"id": "s1", "sequence_id": "seq_...", "length_bp": 1243,
     "support_labels": ["label_..."], "support_kind": "sample_kmer_support",
     "source_coordinates": null}
  ],
  "edges": [{"from": "s1", "from_orientation": "+", "to": "s2",
             "to_orientation": "+", "overlap_bp": 30}],
  "path_view": "available_by_path_id",
  "exploration": {"complete_within_requested_scope": false, "sampled": false,
                  "stop_reasons": ["max_states"]},
  "delivery": {"has_more": true, "next_cursor": "cursor_..."}
}
```

This fragment illustrates the schema; paginated segment/edge references resolve through the same
artifact. Distances along a graph path are separate from source coordinates. Reconverging paths can
reach one segment at different offsets: report per-path values or ranges, not one fabricated exact
position. Coordinate-bearing results identify source record, strand and an explicit convention
(use 0-based half-open sequence intervals); raw k-mer coordinates are converted accordingly.

Store the bulk graph as a versioned, content-digested artifact with optional GFA/FASTA exports.
The digest identifies the exact saved result, including a partial result. Initial MCP output is
a bounded aggregate plus representative contiguous sequence windows. Let the agent fetch/slice
specific segments or paths through sequence handles. **Do not reduce away all sequence evidence:**
the agent may need to inspect a continuous tract to recognize a repeat or other feature. Report
window boundaries/truncation; downstream analysis tools use exact stored sequences, not model copies.

### 5.3 Sequence paths with labels along the way

For each extracted sequence, return a `path_id`, immutable `sequence_id`, ordered oriented segment
references, and an interval annotation describing label support **along that path**. Run-length
encode intervals where the label set stays constant; use a shared label dictionary instead of
repeating large strings. Mark whether an interval annotates bases or k-mer starts and define overlap
handling. The representation uses 0-based half-open intervals of path k-mer starts;
for a path of length L and k-mer size k there are L−k+1 starts. Biological feature annotations use
base intervals and remain a separate layer.

For v1 `within_label`, a single checked label can be represented as one support interval; declare
that other labels were not enumerated. Richer output keeps three distinct fields: local interval
labels, common labels supporting the entire path, and optional source-record traces. Local labels can
change even when one permitted label survives throughout. In union mode the common-label set may be empty; that must remain visible. Every path
also carries its anchor, selection/sampling method, length, stop reason and provenance.

### 5.4 Small subgraphs as first-class MCP analysis objects

Mint a `subgraph_id` for the bounded saved graph, linked to its traversal job, index generation,
anchors, strategy, support states and artifact digest. This differs from a `sequence_id`: it refers
to topology plus evidence, not just sequence bytes. A multi-shard job returns a collection of such
objects. Comparing them does not make them one connected graph.

The agent can work locally on these artifacts before deciding whether to expand the original
traversal. Proposed MCP tools, with exact naming to settle during implementation:

| Tool | Purpose | Bounds/evidence contract |
|---|---|---|
| `inspect_subgraph` | Summarize size, anchors, permissible branches/cycles, label coverage and exploration boundaries; page selected segments/edges | Small summary by default; preserve exploration versus delivery completeness |
| `extract_subgraph_paths` | Return sequence handles and interval labels for selected anchors/endpoints, or a specified bounded sampling policy | Enforce saved support/occurrence compatibility; report if requested context lies outside the saved graph; never enumerate unlimited walks |
| `find_subgraph_motifs` | Find supplied exact motifs on admissible paths, including across segment boundaries | Bound matcher/path state; return actual matching path evidence; no implied biological feature identification |
| `summarize_contexts` | Group saved paths/annotations across a traversal job by sequence signature, feature, label or length | Declare denominator and grouping/annotation version; retain member IDs and unknowns |
| `compare_contexts` (later) | Compare selected neighbourhoods, their sequence/feature differences and evidence | State comparison/alignment method; never join separate shards or infer unobserved physical linkage |

These graph tools follow the first extension workflow. In v1, reuse existing sequence slicing,
translation and verification tools on saved flank handles; general graph inspection/path extraction
is not a prerequisite. When graph analysis is added, start with inspection and bounded extraction.
Motif and comparison tools follow concrete analysis needs.
Keep a small typed tool surface rather than exposing arbitrary executable traversal code.

Local analysis still has state/output/time budgets. A derived filtered subgraph keeps its parent
digest and selection criteria and cannot claim more completeness than the parent. Paths extracted
from a compact topology must respect the saved compatibility states: new combinations of individually
valid segments are not automatically supported paths. If more context is needed, the response points
to the relevant anchor/frontier and suggests a new traversal strategy rather than fetching more of
the large graph implicitly. All local analysis cost is recorded separately from archive traversal.

---

## 6. Refinement, caching and async integration

The expected loop is:

```text
search → freeze hit set → declare strategy → submit traversal batch
       → inspect aggregate → inspect small subgraphs / extract and analyse sequences
       → refine strategy if context is missing → apply to the same hit set → compare → answer
```

Three operations have different costs and validity:

1. **Re-summarize stored observations:** change grouping/filtering without graph work, provided
   required fields are already present. A summary cannot infer branches or flanks never explored.
2. **Re-analyze stored sequences:** apply a new detector/annotation method without traversal if
   the saved sequence context is sufficient; record a new analysis version and cost.
3. **Re-traverse:** widen the radius, change path constraints, request different sampled paths, or
   investigate missing context. Resume only if a complete compatible frontier was saved; otherwise
   start a new run. A checkpoint must retain support/matcher states, frontier, visited state, budgets
   and implementation version; a compacted output graph alone is insufficient. Relaxed labels/stops
   may expose previously pruned paths, so reuse is not automatic.

Persist parent job/strategy IDs and offer a comparison over identical anchor IDs: newly resolved
cases, changed detections, remaining unknowns, and incremental cost. Changed hit sets, index builds,
support modes or feature detectors must be flagged as changed comparison conditions. Do not pool
repeated refinement runs as independent biological observations.

Use an operation-specific exact cache key containing graph/annotation generation, anchor/hit-set
digest, normalized strategy, algorithm version, effective budgets and RNG policy. Feature-analysis
and summary caches additionally include detector versions or grouping definitions. Existing search
cache reuse across `top_labels`/presence queries is unsuitable. Initially bypass it for traversal.
Keep partial artifacts with their stopping reasons; they must not satisfy a request for completed
exploration. Cache hits report both original compute and work actually incurred in the current run.

**Async service implementation points** (paths relative to `metagraph-search-service`):

- `app/metagraph/client.py`: factor shared transport, add a dedicated traversal request builder and
  parser. The current `/search` URL and parameter allow-list are search-specific; do not reuse the
  allow-list unchanged.
- `app/database_registry.py`: advertise graph generation, traversal modes, coordinate/trace support,
  label unit, counts availability and server limits per index. Separate capability from per-anchor
  evidence coverage.
- Queue a dedicated traversal job through existing per-database capacity controls. Preserve
  semaphore acquisition, retry/backoff, cancellation and stale-attempt protection explicitly.
- Partition a batch by exact shard/generation, with bounded batch sizes and fair per-anchor quotas.
  Global exhaustion reports unstarted anchors; the first jobs to finish are not a representative sample.
- Deduplicate retries by immutable work-item ID; biological counts include each selected item once,
  while cost accounts for all attempts. Retain failed, unsupported and unresolved items in denominators.
- Persist traversal jobs/manifests and graph artifacts separately from scored search-hit rows;
  reuse storage/pagination infrastructure without passing graphs through search score normalization.
- Reuse `app/seqhandles.py` for immutable sequence content. Current handles contain no graph/path
  evidence and expire; persist provenance and sequence artifacts, then remint handles as needed.

### 6.1 Pilot runs and progressive cohort expansion

A pilot applies the same strategy to a small reproducible subset of the frozen input list. A request
such as `pilot: {n_anchors: 10}` is shorthand for an explicit saved selection, not a separate traversal
algorithm. Select across available sequence diversity, label types or other relevant strata; record
selection rules and RNG seed where applicable. Keep eligible, selected, attempted and completed
counts distinct. This is **sampling loci**, separate from sampling paths through one neighbourhood.

Return useful sequence context plus distributions of extension length, branching, states, annotation
work, elapsed time, memory and stop reasons. Any full-cohort cost extrapolation is approximate and
must state assumptions: a tiny pilot can miss expensive tails, and warm caches/concurrency alter
cost. Enforce actual job budgets independently of the estimate. Prefer the smallest informative
windows initially, then widen or expand promising/unresolved cases; additional annotation on saved
artifacts should not trigger a new graph walk.

Evaluate whether the agent uses the pilot to choose an informative, affordable strategy rather than
blindly extrapolating or repeatedly exhausting caps. Retain per-attempt costs and selection history;
adaptively selected results are not an unbiased prevalence estimate of the original archive.

### 6.2 Later: contrastive contexts across two cohorts

Apply the same normalized strategy to two frozen hit sets, such as distinct environments or time
periods, and compare flank groups, trie branches or gene arrangements. Reuse saved contexts where
compatible. This is a bounded local analogue of differential-context analysis; it does not require
running whole-graph differential assembly. Begin with group frequencies, effect sizes, uncertainty
and inspectable members before adding formal enrichment tests.

Both cohorts need the same evidence rules, comparable index coverage/budgets and explicit complete/
partial/unknown denominators. Record overlap, repeated samples and lineage/study structure; define
cohorts independently of the feature being tested. Unequal truncation or missing metadata must not
become a biological difference. Statistical enrichment additionally needs a declared method,
confounder/duplicate controls and multiplicity handling for many tested branches. This is Phase C,
not merely two manifests plus an unqualified p-value; it can be useful without new core traversal.

---

## 7. Implementation phases and acceptance tests

### Phase A — measured extension plus a thin end-to-end agent loop

Build a reusable `src/graph/traversal/` operation with a small CLI/test harness, then expose the same
primitive through core `/traverse`, the async job and one MCP workflow. The release scope is §0:
one verified ATB shard/release, explicit contiguous sequence seeds, fixed-label bidirectional
extension, branch/edge-reuse/budget stops, saved flank evidence, exact grouping and pilot selection.
Initial branch-stop extension is deliberately simpler than general branching BFS. Where a cycle's
exit requires branching, stop and report it; the one-edge rule is an upper bound, not a promise to
enumerate every once-unrolled path. Preserve topology/boundaries needed for honest inspection.

Exercise batches from the first harness. Measure warm/cold latency, peak memory, annotation work,
branching and returned context; validate against source-contig truth. Compare targeted membership
with batched-row access on the actual annotation representation. Avoid graph-sized per-request
arrays, global unitig scans or a new RowDiff accessor before measurements justify it. Complete the
small agent loop before treating the remaining architecture as a release prerequisite.

In parallel, expose existing original-query support signatures where the deployed service can
preserve them. Once saved flanks exist, prioritize compressed trie summaries, search-based annotation
and metadata joins as independent service/analysis additions. Design the §9.4 benchmark at this stage
so gains from existing evidence are distinguished from gains requiring new traversal.

### Phase B — connecting hits and richer evidence inspection

Add bounded branching, two-hit reachability and reference-guided inside-query investigation.
Three-plus-hit/operon sets require a selected biological question and measured state bounds.
Introduce literal motif stops and saved-subgraph inspection/path extraction for chosen demonstrations.
Add search-derived hit selection/manifests and richer anchors when explicit lists become limiting;
keep source-shard provenance immutable. Exact grouping precedes approximate context clustering.
Reuse established feature-analysis tools on sequence handles rather than implementing biological
recognition inside the graph walker. Keep deterministic item budgets and all-attempt accounting.

### Phase C — advanced support and exploration

Add multi-label `path_common`/explicit union/intersection modes, strict source tracing only on
validated coordinate indexes, reproducible path sampling and advanced annotations as needed.
Source tracing remains unavailable on unsupported deployments; first demonstrations must not
assume it. If physical source flanks are essential, use verified source retrieval or choose a
validated coordinate index. Resumption, abundance/diversity weighting and coverage-drop stops
follow measured need and appropriate count semantics. Contrastive context analysis and bounded
bubble/allele inspection are later application-driven additions, with their own evidence and
validation requirements. No delivery-time estimate is established.

### Required acceptance cases

Apply each case when its capability is introduced; later-mode tests are not prerequisites for the
initial extension prototype.

| Case | Required behaviour |
|---|---|
| Linear sequence | Correct left/right sequence, seed included once, extension bases counted outside seed |
| Partial seed / changed release | Disconnected supported blocks are not joined; a still-matching seed does not mask a release change |
| Annotation backend | Direct-cell and batched-row predicates agree; no online whole-column/graph scan fallback |
| Pilot selection | Saved membership and selection seed reproduce; estimates remain distinct from measured full-job cost |
| Inside-query support | Masks reconstruct k-mer-start intervals in the declared original/aligned coordinate frame; gaps do not trigger automatic biological diagnoses |
| Flank trie | Stored windows and memberships reconstruct; censoring is not divergence; samples on multiple branches are not forced into a partition |
| Reference-based annotation | Feature names follow verified matching intervals; a top whole-chromosome title is not a feature call |
| Metadata/contrastive summaries | Missing dates, partial coverage, overlapping cohorts and repeated samples remain visible in denominators |
| Bubble phase | Marginal support at multiple variant sites does not fabricate a complete allele |
| Substitution neighbourhood | Exhaustive results agree with brute-force Hamming distance on a tiny graph; a present double mutant is found even if single-mutant intermediates are absent; canonical mapping preserves seed-relative edits |
| Neighbour-search budgets | Membership, label and output completeness remain distinct; capped work is not reported as absence |
| Sample-switching chain `{A,B} → {A} → {A,B} → {B}` | Union may accept; path-common rejects the complete path |
| Diamond with differing surviving labels | Retain valid continuations after reconvergence; no unsupported prefix/suffix combination |
| Same-label repeated motif / multiple molecules | Candidate graph paths remain distinguishable from occurrence-specific source traces |
| Cross-record boundary fixture | Strict source tracing stops at the record boundary even with consecutive global k-mer coordinates |
| Cycle / self-loop / overlapping cycles | Preserve cycle edges; enforce directed-edge reuse per candidate extension, including joined arms/legs; report unresolved multiplicity |
| Multi-hit joint path | Pairwise positives using incompatible labels or different middle-hit occurrences do not establish a joint witness |
| Gene-hit layout | Preserve selected occurrences, order, strand, signed gaps and total span; partial-hit gaps are not presented as intergenic distances |
| Reverse complements / palindromes | Equivalent oriented contexts; canonical annotation lookup does not erase traversal orientation |
| Support changes within an unbranched segment | Support change preserved by splitting or explicit annotation runs |
| Motif spans segment boundary | Stop/detection follows the admissible sequence with overlaps counted once |
| Caps, cancellation, failures | Explicit partial state and stop reasons; never silently counted as negative |
| Batch reorder, retry, duplicated occurrence | Stable item identity, deterministic selection, no duplicate biological counts, all attempt costs retained |
| Stale anchor / missing label / unsupported trace mode | Actionable error, no silent remapping or evidence downgrade |
| Sampling on a tiny known graph | Reproducible RNG policy and expected transition frequencies; no uniform-path claim without justification |
| Aggregate drill-down | Counts reproduce from member evidence, with declared unit and complete/partial/unknown denominators |
| Saved-subgraph analysis | Inspection/extraction makes no large-index calls; missing context remains explicit; extracted paths preserve support compatibility |
| Sequence interval labels | Label runs reconstruct every path position with correct k-mer/base coordinates and common-path support |

For the paper, compare a fixed strategy, a scripted refinement policy and an adaptive agent using
the same backend capabilities, frozen hit sets and total budgets. Measure supported context recovery,
false linkage claims, detection/rejection quality, handling of unknowns, refinement efficiency and
cost. Use independent loci/families for generalization; repeated runs characterize stochasticity.

## 8. Remaining decisions

1. Which ATB shard/release and annotation backend will host the first pilot, and can its source
   contigs be retrieved for validation? Verify actual label, coordinate and count semantics.
2. Which artifact-manifest release identifier will jobs record initially, and when are richer
   persistent anchors worth introducing?
3. Which biological question determines the first detector and required support mode?
4. What selection/deduplication unit and occurrence policy best matches that question?
5. What measured per-anchor/job budgets make the initial interactive loop practical?
6. Which result segments/windows should be included initially, and which require explicit drill-down?
7. When is saved-frontier resumption worth its state/storage cost, compared with exact reruns?


## 9. Biological demonstrations and the paper's capability claim

The intended contribution is that an agent can recover and compare relevant sequence context that
was absent from the original query, using bounded operations with inspectable evidence. MCP/chat
accuracy remains central: test whether the agent chooses useful follow-ups, reads the evidence,
handles incomplete results and reaches supported biological conclusions. Traversal is a capability
extension, but exposing an endpoint alone does not establish a discovery or a publication claim.

### 9.1 An Anthropic-like discovery workflow

The complete workflow is **sensitive family discovery → genomic anchoring → neighbourhood retrieval
→ biological analysis → cohort comparison/refinement**. MetaGraph traversal supplies the context
stage. Complementary capabilities include protein/profile homology search, protein-to-nucleotide
mapping, gene/strand annotation, repeat/domain analysis and, when justified, RNA/structure analysis,
alignments and literature verification. They operate on exact sequence handles or versioned source
records. Precomputed annotations may accelerate a census; unannotated/noncoding sequence must remain
accessible so summaries do not hide unexpected biology.

Anthropic's ART report describes a search space of about 1.94 billion protein-cluster representatives
but a neighbourhood campaign of 10,983 loci with up to 10 kb on either side. This is consistent with
representative-locus selection followed by deeper investigation, not unrestricted graph exploration
of the whole database. Nominal flank volume is about 220 Mb before truncation, excluding seed genes;
that arithmetic is not a runtime estimate. The report says 85.9% of windows hit contig ends, motivating
an evaluation of graph context without presuming it can correctly bridge those ends. Database
preparation/preannotation and agent campaign costs must be reported separately.
Source: [Anthropic technical report, Methods pp. 28–30](https://www-cdn.anthropic.com/22573675ada52a8ca8a97a1a4b4326b2f208a071.pdf).

### 9.2 Candidate demonstrations

| Investigation | Traversal/analysis question | Evidence needed |
|---|---|---|
| blaNDM context, including proposed Tn125/ISCR1 comparisons | Recover flanks and compare local element arrangements across selected loci | Source-supported sequence plus separate element annotation; context alone does not establish mobility |
| mcr-1 neighbourhoods | Compare flanks, including supplied ISApl1 hits where relevant | Validated sequence layout; plasmid membership is a separate claim |
| Caffeine-pathway genes ndmA–D and other gene sets | Find compatible joint paths, order, strand and gaps; identify recurring modules | Full-gene boundaries for intergenic gaps and separate operon evidence |
| Query-support investigation | Explain positional support and seek compatible paths through partially matching intervals | Distinguish original-query masks, graph witnesses and source continuity; gap patterns alone are not diagnoses |
| Near-hit variants | Find indexed k-mers within one or two substitutions of a selected hit and compare label support | Exhaustive declared radius/scope or explicit caps; longer context needed for locus/allele identity (§4.8) |
| Allele/bubble inspection | Identify bounded, supported sequence alternatives through a locus | Compatible evidence across variant sites, reference verification and unresolved phase; later analysis increment |
| TMPRSS2–ERG interpretation | Test whether within-query and flanking evidence improves interpretation of low-containment junction hits | Suitable human-index availability is unverified; require junction-specific validation and independent evidence; later than ATB |
| Contrastive gene contexts | Compare flank/trie branches or gene arrangements between frozen cohorts | Comparable evidence/coverage, metadata missingness, lineage/duplicate controls and declared statistical scope |
| ART-like neighbourhood investigation | Discover recurrent or unusual coding/noncoding context around a family | Sensitive upstream family search, source/repeat validation and held-out biological cases |

A fusion demonstration should test the proposed **measured failure → targeted tool/evidence →
demonstrated correction** chain, not assume it succeeds. Both arms individually resembling the
expected genes is insufficient without junction-specific support. Compare the existing search/
alignment workflow (including the available `align=true` arm where applicable), direct source-window
retrieval and the traversal-enabled workflow on matched cases. A matched-evidence condition can
separate benefit from access to new sequence from benefit of the interface/refinement strategy.

### 9.3 Validation before claims of discovery or efficiency

Start with known loci and independent source-contig truth, including repeat-rich regions, contig
ends, alternative occurrences and negative controls. Use the same anchors to compare direct source
extraction, current search/alignment and graph traversal. Measure correct flank/layout recovery,
false joins, useful additional context, repeat ambiguity and appropriate unresolved outcomes.
Scale a representative cohort from small pilots toward thousands of loci, reporting latency
(distributions, not just means), memory, annotation work, explored states, returned bytes and total
agent/compute costs including failures and retries. No source-code audit establishes speedup.

Then compare fixed strategies, scripted refinement and adaptive agents with identical underlying
capabilities and total budgets. Test evidence inspection, pilot use, unsupported-negative claims,
and detection of relevant biology. Keep truth/reference answers independent of the evaluated
pipeline; use independent loci/families/cohorts as the evaluation units and repeated runs to measure
stochasticity. An ART replay demonstrates recoverability; held-out contexts are needed to test
general investigation capability. Include lineage/diversity controls when claiming recurrence or
enrichment, rather than treating closely related samples as independent discoveries.

### 9.4 Capability-dependent tasks and selective tool use

Add an E2 task family whose answer depends on sample-specific within-query or neighbourhood evidence,
and evaluate it in the E5 tool-family ablation. Independently establish what the available evidence
can answer in each condition; do not label a task "traversal-dependent" just because traversal was
withheld when an existing tool could solve it. Include four classes:

- Existing scores/positional signatures already suffice; unnecessary traversal adds cost.
- Existing alignment or direct source retrieval supplies the needed evidence.
- Additional graph context resolves a question the preceding capabilities leave unanswered.
- Even the permitted traversal remains ambiguous/truncated; the supported answer stays unresolved.

Compare those capability levels on held-out source-validated cases with matched total budgets and
all-attempt cost accounting. Include a matched-evidence condition that supplies the same recovered
sequence to separate information access from presentation and agent-planning effects. Keep query
profiles and source truth independent of the model's proposed diagnosis; avoid answers recoverable
solely from generic knowledge or misleading benchmark hints.

Measure whether the agent identifies missing evidence, selects an appropriate tool or pilot, reads
the returned sequence/support, revises its answer correctly, and stops when evidence is sufficient.
Score justified uncertainty as well as supported answers, unsupported linkage/allele claims and
wasted traversal. The desired delta is improved supported biology per budget, not more tool calls.
This complements the existing truncated-aggregate failure mode with the converse: knowing when
additional context is actually needed.

The original MCP/QA evaluation can proceed independently of this extension. Journal positioning,
paper timelines and whether traversal becomes the centerpiece belong in
`/Users/raetsch/git/tex/2026/metagraph-agent-paper/PAPER-PLAN.md`; this design adds no venue or schedule
commitment and does not assume a successful discovery in advance.
