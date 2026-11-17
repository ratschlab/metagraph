# Pair identity in multi-graph answers

Status: proposal for review, 2026-11-17. It proposes how the staging server (branch `gr/labeled-traversal`, which
the service below runs) takes two commits of `master`, 356fe1ee3 and 82137bee7, and how its answers name a pair
once they stop stating file paths. Nothing here is implemented yet. Each decision in section 3 is a
recommendation; the project owner may overturn any of them.

## 1. Context

### 1.1 Pairs and names

A multi-graph server reads a graph list at start-up, a CSV file with one line per pair:

```text
refseq,/data/refseq/bact_other.dbg,/data/refseq/bact_other.annodbg
refseq,/data/refseq/bact_pseudo.dbg,/data/refseq/bact_pseudo.annodbg
```

A **pair** is one graph file with one annotation file. A **name** (the first column, also called the graph name)
groups pairs: a client selects names, and one name can hold several pairs. Production RefSeq is one name with 11
pairs; a human SRA name is one graph with several annotations. Each pair can also have a manifest, from which the
server computes `index_fp`, a fingerprint of the files' contents.

### 1.2 Who names pairs today, and how

On the staging server every statement of a pair is `graph` (the name) plus `graph_path` and `annotation_path`,
spelled as the graph list spells them:

| Who | Where it names a pair |
| --- | --- |
| Engine answers | the multi-graph **envelope** of `/resolve`, `/traverse` and `/pattern`, which holds one **entry** per pair: the pair's answer with its **tags** `graph`, `graph_path`, `annotation_path`, `index_fp` and `outcome`. Also the envelope's `not_reached` and `unfinished` lists; `/search` rows and their `cuts`; the request registry (`GET /requests/{id}`); the `candidates` of a selection refusal; examples in `GET /capabilities` |
| Request selectors | `pairs: [{graph, graph_path?, annotation_path?}]`, `graph` + `graph_path`, and the probe parameters `?graph=&graph_path=&annotation_path=` |
| `GET /capabilities` | `graph_summary.pairs[]`: each pair with its paths, `index_fp` and capability profile |
| The service | the API service in front of the engines: REST, and MCP (Model Context Protocol) tools for AI clients. It selects pairs by paths from its own manifest, stores traversal results, and states a public `pair_id` of its own in search answers |
| The library | the Python traversal library (`metagraph/api/python/metagraph/traverse`). It parses a `/traverse` answer into a **graphlet** document. A stored graphlet keeps its pair's entry, tags included, on its first line (the J line). The library re-asks a pair by the paths in those tags. It combines the graphlets of one seed across pairs into a **union**, which matches pairs by their tags and binds a digest to them |
| The suite | the end-to-end test suite (its own repository), which pins the key sets above and compares answers |

### 1.3 What master changed

Master serves `/traverse` on a graph list in its own way (a list of walks per seed, not the staging envelope).
Its two commits, in their own words:

- 356fe1ee3, "State chosen record places and opaque pair IDs in traversal answers": "Header labels now include
  record-local `place` (`origin`, `dir`). List-server walks identify their pair with `graph` and an opaque
  `pair_id`. The ID uses a versioned UUIDv5 encoding of the configured graph name and index paths. Raw paths are
  omitted from the response. This identifies a configured pair; content identity still comes from its deployment
  manifest. [...] Deploy the compatible API adapter before updating native servers so manifest binding can resolve
  the new IDs."
- 82137bee7, "State match-count cuts in search results": "Capped JSON query results now include `cuts` with matched
  and returned counts, completeness and cutoff tie facts [...]. List-server rows and cut records identify their
  configured pair with `graph` and an opaque `pair_id`, without raw index paths."

The id (`src/common/configured_pair_id.hpp`) is `n1-` followed by a lowercase UUIDv5 in the namespace
`d3565c73-e7e9-5bca-b6a9-64cb5c4c9f93`. The UUID's name concatenates, for the graph name, the graph path and the
annotation path, the field's UTF-8 byte length, a colon and its bytes; an absent field is `-1:`. The encoding
normalizes nothing. Master's docs: the id "is deterministic public metadata. Treat it as an identifier, never as a
credential. It is an output field. Existing request selectors keep their documented names and syntax."

The staging server already has `cuts` (one per pair and query sequence, with tie facts, a median tie-break on count
annotations and per-unit refusals). Master's `cuts` is a smaller form of the same idea. What is new for staging is
`pair_id` in place of the paths, the `matched` count, `cuts` in the CLI's JSON output, and `place`.

## 2. The problem

The owner made two decisions on 2026-11-17. First, take master's form: an answer names a pair by `graph` and
`pair_id`, and raw index paths are omitted. Second, the documents the service shows through MCP keep the host paths
in their pair tags, and the library needs no change. The library reads pairs only from those tags. It re-asks
by them, and its union matches and digests by them. An engine that no longer states paths therefore conflicts with
an unchanged library, unless something puts the paths back before the library reads a document.

## 3. The proposal

The engine speaks master's form. The service, which writes every document the MCP shows, restores the paths. The
decisions:

1. **Answers name a pair by `graph` and `pair_id`, never by paths.** This covers envelope entries, `not_reached`
   and `unfinished`, `/search` rows and `cuts`, registry pairs, selection `candidates`, and the examples in
   `GET /capabilities` (`headers_unique`, `labels_disjoint`). It is the owner's first decision and master's form.
   A single-graph server states `pair_id: null` where it stated null paths.
2. **The engine computes the id as master does**, once per configured pair at start-up, from the graph list's
   first three columns as spelled. A pair listed under two names has two ids, since the id hashes the name. The
   header and its test vectors come from master unchanged.
3. **`GET /capabilities` keeps `graph_path` and `annotation_path` beside `pair_id` in `graph_summary.pairs[]`.**
   This is the one place where the engine states the binding from id to files. It is one statement per pair, not
   per row. The service can check its own binding against it, and the library's pair discovery keeps working.
   It also adds `multigraph.pair_id: "n1"`, so a client knows which form the engine answers in before it parses
   an answer.
4. **Path selectors stay, with no end date, and `pair_id` becomes a selector too.** The forms
   `pairs: [{graph, pair_id}]` and `?graph=<name>&pair_id=<id>` are new. The library and the service build path
   selectors from what they store. A client that holds only an answer (an entry, a candidate) needs the id to name
   a pair. One object names a pair by its id or by its paths, not both. `/search` keeps `graphs` only, as today.
5. **Texts.** The selection texts that are frozen byte for byte stay unchanged, because path selectors stay and
   the texts remain true. They include `index 'refseq' spans several graphs (...); pass 'graph_path' to pick one`.
   The staging server's own pair texts name the pair as `(<name>, <pair_id>)`: not started, unfinished at the
   cutoff, loading into RAM, and a file changed on disk. Log lines keep the paths: they are the operator's view,
   and an operator acts on files.
6. **`/search` cuts.** Each cut gains `matched`: the number of labels above the request's thresholds before the
   `top_labels` cut. A unit (one query sequence on one pair) refused because its answer would pass the memory
   bound keeps the count its search computed. The CLI's `--json` output gains `cuts`, with null pair tags. The
   staging server keeps its own members (`nth_median`, `tied_same_median`, `refused`), its row order (score, then
   median abundance, then label), and `cuts` in every mode. Master omits `cuts` for `abundance_sum`. Removing
   members or modes would take away facts clients read.
7. **The service restores the paths, and the library stays unchanged.** The service computes each pair's id from
   its manifest, as master's code does. It writes `graph_path` and `annotation_path` into every document it stores,
   beside the engine's `pair_id`, and keeps selecting by paths. The library then reads the same tags as today, and
   its digests stay the same (section 5.3).
8. **`place` is not adopted.** Master's walker chooses one occurrence of the seed in each record and follows only
   that one, and `place` states it. The staging walker has no such choice. Under `support: trace` it follows every
   coordinate chain of the record, and under the default k-mer support it follows membership. A `place` would
   claim a choice the walk never made. The staging server already states every occurrence in
   `coordinates.seed[]` when `output.coordinates` is true.

## 4. The contract, before and after

The examples shorten the paths and write an `index_fp` as `3f1c...`. The ids are the real ids of the two lines
in section 1.1.

**An envelope entry** (`/resolve`, `/traverse`, `/pattern`), and `not_reached`:

```jsonc
// before
{"graph": "refseq", "graph_path": "/data/refseq/bact_other.dbg",
 "annotation_path": "/data/refseq/bact_other.annodbg", "index_fp": "3f1c...", "outcome": "answered", "...": "..."}
"not_reached": [{"graph": "refseq", "graph_path": "/data/refseq/bact_pseudo.dbg",
                 "annotation_path": "/data/refseq/bact_pseudo.annodbg"}]
// after
{"graph": "refseq", "pair_id": "n1-64b2a5df-621c-5503-96c1-4b4b95dbd8af", "index_fp": "3f1c...",
 "outcome": "answered", "...": "..."}
"not_reached": [{"graph": "refseq", "pair_id": "n1-5b5559a9-cef2-5bd7-b378-c77612f9f9fc"}]
```

**A `/search` row and its cut** on a multi-graph server:

```jsonc
// before
{"sample": "NC_003279.8", "kmer_count": 170, "graph": "refseq",
 "graph_path": "/data/refseq/bact_other.dbg", "annotation_path": "/data/refseq/bact_other.annodbg"}
{"graph": "refseq", "graph_path": "...", "annotation_path": "...", "index_fp": "3f1c...",
 "returned": 2, "complete": false, "nth_score": 139, "tied": 3, "tied_returned": 1}
// after
{"sample": "NC_003279.8", "kmer_count": 170, "graph": "refseq",
 "pair_id": "n1-64b2a5df-621c-5503-96c1-4b4b95dbd8af"}
{"graph": "refseq", "pair_id": "n1-64b2a5df-621c-5503-96c1-4b4b95dbd8af", "index_fp": "3f1c...",
 "matched": 5, "returned": 2, "complete": false, "nth_score": 139, "tied": 3, "tied_returned": 1}
```

On a single-graph server the cut states `"graph": null, "pair_id": null` and its `index_fp`.

**A selector**: both of these name the first pair; only the second is new.

```json
{"pairs": [{"graph": "refseq", "graph_path": "/data/refseq/bact_other.dbg"}]}
{"pairs": [{"graph": "refseq", "pair_id": "n1-64b2a5df-621c-5503-96c1-4b4b95dbd8af"}]}
```

`GET /traverse/capabilities?graph=refseq&pair_id=n1-64b2a5df-621c-5503-96c1-4b4b95dbd8af` probes it. An id that
no pair of the name has is refused 400 `pair_not_found`, with the name's pairs as `candidates`.

**A refusal text and its candidates:**

```jsonc
// before
{"error": "resolve: pair (refseq, /data/refseq/bact_pseudo.dbg, /data/refseq/bact_pseudo.annodbg) was not started before the request's work cutoff (...)"}
{"error": "...", "code": "pair_not_found", "candidates": [{"graph": "refseq", "graph_path": "...", "annotation_path": "..."}]}
// after
{"error": "resolve: pair (refseq, n1-5b5559a9-cef2-5bd7-b378-c77612f9f9fc) was not started before the request's work cutoff (...)"}
{"error": "...", "code": "pair_not_found", "candidates": [{"graph": "refseq", "pair_id": "n1-64b2a5df-..."}]}
```

**`GET /capabilities`:**

```jsonc
// before
"graph_summary": {"pairs": [{"graph": "refseq", "graph_path": "/data/refseq/bact_other.dbg",
                             "annotation_path": "/data/refseq/bact_other.annodbg", "index_fp": "3f1c...", "k": 31}]}
// after
"graph_summary": {"pairs": [{"graph": "refseq", "pair_id": "n1-64b2a5df-621c-5503-96c1-4b4b95dbd8af",
                             "graph_path": "/data/refseq/bact_other.dbg",
                             "annotation_path": "/data/refseq/bact_other.annodbg", "index_fp": "3f1c...", "k": 31}]},
"multigraph": {"pair_id": "n1", "selectors": {"/traverse/capabilities": ["graph", "graph_path",
               "annotation_path", "pair_id"], "...": "..."}}
```

Unchanged: single-graph answers (they name no pair), `/align` (a multi-graph server refuses it), `index_fp` and
`composite_fp` (content identity), the capability profile's `support_fp` (it excludes identity), and log lines.

## 5. Compatibility and migration

### 5.1 Deploy order

Master's instruction holds: the service's adapter first, then the engines.

1. **The service** reads both forms. A current engine states paths; a new one states `pair_id`, and its
   capabilities state `multigraph.pair_id`. The service binds each id to a pair of its manifest and checks every
   computed id against `GET /capabilities`. It writes the paths into the documents it stores and keeps sending
   path selectors. It sends a `pair_id` selector only to an engine that states `multigraph.pair_id`, because a
   current engine refuses an unknown member of a `pairs` object.
2. **The suite** accepts both forms in its key-set pins.
3. **The engines**, first on staging and then in production. A rollback to a current engine brings the paths back,
   and the adapter reads that form too.

### 5.2 Old stored documents and caches

Stored documents have paths and no `pair_id`, and the unchanged library reads them as today. Their id can be
computed from their paths. It is deterministic, so an old and a new document of the same configured pair get the
same id. The service's caches keyed on paths keep their keys if new answers are mapped to paths before caching. The
engine's `/search` result cache lives in the process, and the restart that deploys it clears it.

### 5.3 The library's digests

- A graphlet's body digest and the resume tokens of single-graphlet operations cover the records after the J line,
  not the J line, so no tag is an input.
- A union's digest covers each part's key: release, index namespace, `graph`, `graph_path`, `annotation_path`,
  `index_fp` and `index_meta_fp`, read from the J line. It ignores `pair_id`. A document whose paths the service
  restored has today's key, so its union digest, union resume tokens and export tags are unchanged.
- The hash of the deployed library (`library_sha256`, which the service reports) covers code, not documents. No
  library code changes, so it does not move.

### 5.4 What breaks, for whom

| Who | Effect | Remedy |
| --- | --- | --- |
| Clients that read `graph_path` or `annotation_path` from answers | the members are gone | read `pair_id`, and bind it through `GET /capabilities` or their own manifest |
| The service | must not forward the engine's `n1` id under its own `pair_id`, which states its public id: `p2-` and the first 32 hex digits of SHA-256 over `[database, n1]` | map `n1` to its public id; code that tells cuts of a multi-graph server from single-graph ones by `graph_path` must test `graph` or `pair_id` |
| The suite | pinned key sets (entry tags, row tags, cut keys, selectors), pairs matched to graph-list rows by paths | accept both forms; match by the id of the row |
| Direct library users on raw engine answers | for a name with several pairs, the tags carry no paths. A re-ask names the graph only and is refused as ambiguous. Without manifests the union cannot tell such pairs apart (it matches `index_fp` first, so with manifests it stays exact) | use documents the service restored, or a later library change that reads `pair_id` (section 9) |
| Requests near the `/search` memory bound | a row's tags get shorter (39 bytes instead of two paths), so the answer's price drops and a request refused before may be answered | none needed |
| Path selectors, the frozen selection texts, single-graph clients | unchanged | none |

## 6. Alternatives

| | Design | Why not |
| --- | --- | --- |
| B | the engine adds `pair_id` and keeps the paths | contradicts the owner's first decision and master's form (master tests that no path appears in walks, rows or cuts) |
| C | the library keys pairs on `pair_id` (tags, re-asks, union keys, digests) | contradicts the second decision: every union digest and union token of a document with an id changes, and the library hash moves |
| D | a server flag that keeps paths in answers for one deployment | two answer forms to test and document, and the deployed fleet would still state paths |

## 7. Risks and what the review should check

- **Encoding stability.** The id must stay byte-identical in the engine (C++) and the service (Python). Check that
  both sides run the five vectors of `tests/common/test_configured_pair_id.cpp` (Python's
  `uuid.uuid5` gives all five), and that `n1-` versions the scheme, so that a change gets a new prefix.
- **Paths spelled differently.** The id hashes the strings, not the files. A service manifest that spells a path
  through a symlink, a relative path or a trailing slash computes another id. The check against
  `GET /capabilities` (decision 3) must turn that into a loud error, not a missing pair.
- **The id hashes the name.** Renaming a graph name changes its pairs' ids; a pair under two names has two ids. Check
  that nothing keys long-lived data on an id across a rename.
- **Mixed forms.** The library compares the tag triples when it stitches and unions parts. A raw answer and a
  restored document of one pair do not compare equal. Check that the service restores every document it stores or
  hands to the library, fresh answers included.
- **Ids are public, not secret.** An id is deterministic from a name and two paths, so it confirms a guessed path.
  Check that this is acceptable where answers leave the service.
- **The old engine and the new selector.** A current engine refuses `pair_id` in a `pairs` object. Check that the
  service gates on `multigraph.pair_id`.
- **Answer prices.** The `/search` answer budget prices tags by their length. Check that the price model and its
  tests follow the new member set.
- **Divergence from master.** The staging server keeps paths in `GET /capabilities`, adds a `pair_id` selector,
  keeps its own `cuts` members and leaves out `place`. Check that none of these blocks a later merge from master.

## 8. Plan, tests and cost

| Step | Work | Days |
| --- | --- | --- |
| 1 | Merge master's two commits: master's id header and its test taken as they are; every conflicting hunk resolved to the staging code, whose `cuts` replace master's | 0.5 |
| 2 | `pair_id` in every tag writer, the selector, the capability statement, the texts and the answer prices; unit and integration tests | 2 |
| 3 | `matched`, and `cuts` in the CLI's JSON output | 0.5 |
| 4 | Re-record the contract samples (29 of 52 state paths; the `/search` ones gain `matched`) | 0.75 |
| 5 | Specifications: the envelope spec, the core traversal spec, the pattern spec, the manual pages | 0.5 |
| 6 | Library tests only: a restored document behaves exactly as a path-only one; a raw answer's limits pinned | 0.25 |
| 7 | Verification: all build targets for DNA and protein, unit and integration tests, Python tests, route and CLI comparisons against the current build, memory and timing, a GCC check, CI | 1 |
| 8 | Review and fixes | 0.5 to 1 |

About six to six and a half working days. Tests to add beyond the updated ones:

- the engine's ids equal Python's `uuid5` for the pairs of a test graph list;
- no answer contains a configured path except `GET /capabilities`;
- a `pair_id` selector names the same pair as its paths;
- a restored document gives the same union digest as a path-only one.

**Timing.** Staging changes are deployed in batches called trains. The batch in preparation, train 7, is closed, so
this change belongs in train 8. It should go after two things: the pending change to how `/resolve` holds its
answer as text, which touches the same files and one contract sample, and the service's adapter.

## 9. Open points for the owner

1. Keep the paths in `GET /capabilities` beside the ids (decision 3)? Recommended: yes, as the engine's one
   statement of the binding.
2. Switch the staging server's own pair texts to `(<name>, <pair_id>)` and keep the frozen selection texts
   (decision 5)? Recommended: yes.
3. Accept that the unchanged library cannot tell apart pairs of one name in raw answers without manifests
   (section 5.4)? Recommended: yes for now. Change the library only when a client other than the service needs it.
