# Pair identity in multi-graph answers

Status: accepted, revised after an outside review on 2026-11-17 and updated on 2026-11-19. It states how the staging
server (branch `gr/labeled-traversal`, which the service below runs) takes two commits of `master`, 356fe1ee3 and
82137bee7, and how every party names a pair once no engine response states a file path. The library's part (section
3.9) is implemented on its own, ahead of the engine, in train 7. The engine's part (steps 1 to 6 and 8 of section 8)
is being implemented for train 8. Section 2 lists the owner's decisions; every other decision in section 3 is a
recommendation the owner may overturn.

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
| Engine answers | the multi-graph **envelope** of `/resolve`, `/traverse` and `/pattern`, which holds one **entry** per pair: the pair's answer with its **tags** `graph`, `graph_path`, `annotation_path`, `index_fp` and `outcome`. Also the envelope's `not_reached` and `unfinished` lists; `/search` rows and their `cuts`; the request registry (`GET /requests/{id}`); the `candidates` of a selection refusal; refusal texts; examples in `GET /capabilities` |
| Request selectors | `pairs: [{graph, graph_path?, annotation_path?}]`, `graph` + `graph_path`, and the probe parameters `?graph=&graph_path=&annotation_path=` |
| `GET /capabilities` | `graph_summary.pairs[]`: each pair with its paths, `index_fp` and capability profile |
| The service | the API service in front of the engines: REST, and MCP (Model Context Protocol) tools for AI clients. It selects pairs by paths from its own manifest, stores traversal results as downloadable documents, and states a public `pair_id` of its own (`p2-...`) in its search answers |
| The library | the Python traversal library (`metagraph/api/python/metagraph/traverse`), which the service vendors. It parses a `/traverse` answer into a **graphlet** document. A stored graphlet keeps its pair's entry, tags included, on its first line (the J line). The library re-asks a pair by the paths in those tags. It combines the graphlets of one seed across pairs into a **union**, which matches pairs by their tags and binds a digest to them |
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
credential." This note calls it the **n1 id**.

The staging server already has `cuts` (one per pair and query sequence, with tie facts, a median tie-break on count
annotations and per-unit refusals). Master's `cuts` is a smaller form of the same idea. What is new for staging is
`pair_id` in place of the paths, the `matched` count, `cuts` in the CLI's JSON output, and `place`.

## 2. The owner's decisions

1. **Master's form.** An answer names a pair by `graph` and `pair_id`; raw index paths are omitted.
2. **No path in any engine response.** This covers answers, refusal and error texts (selection, loading, a file
   changed on disk), `GET /capabilities`, the request registry and the probes. Log lines are not responses: they keep
   the paths for operators.
3. **No path in the service's downloadable documents either.** From this change on, the service stores the engine's
   document as it comes: `pair_id` (n1), no paths. The earlier decision that these documents keep the host paths
   holds until this change lands and is retired then. Documents stored before keep their paths.

So the library changes. A new document carries `graph` and `pair_id` only, and the library continues and unions it
by its id. It still reads the path form, which from then on belongs to old documents only.

## 3. The design

### 3.1 Engine responses

Every response names a pair by `graph` and `pair_id`: envelope entries, `not_reached` and `unfinished`, `/search`
rows and `cuts`, registry pairs, selection `candidates`, the examples in `GET /capabilities` (`headers_unique`,
`labels_disjoint`), its pair list, and the texts of section 3.4. A single-graph server states `pair_id: null` where
it stated null paths. Each pair's `index_fp` stays, as content identity. `GET /capabilities` also states
`multigraph.pair_id: "n1"`, so a client knows the form before it reads an answer.

### 3.2 The id

The engine computes the n1 id as master does, once per configured pair at start-up, from the graph list's first
three columns as spelled. The header and its test vectors come from master unchanged. The id is a **configured
identity**: it names a name and two files. A pair listed under two names has two ids, since the id hashes the name.
A union's **pair identity** is content first (section 3.9): two aliases of one pair with manifests are one pair there.

### 3.3 Requests

Path selectors stay as input, with no end date: they still work, and the service and old documents use them. The
forms `pairs: [{graph, pair_id}]` and `?graph=<name>&pair_id=<id>` are new, so a client that holds only an answer
can name a pair. One `pairs` object names a pair by its id or by its paths, not both; a probe with `pair_id` and
`graph_path` or `annotation_path` is refused the same way, a plain 400 without a code. An id that no pair of the name
has is refused 400 `pair_not_found`, with the name's pairs as `candidates`. `/search` keeps `graphs` only, as today.
No response echoes a path a request sent.

### 3.4 Texts

Selection texts that clients matched byte for byte change, because they list configured paths. The codes stay; the
`candidates` state `{graph, pair_id}`. The new wording:

| Case | Before | After |
| --- | --- | --- |
| a name with pairs on several graphs | `Bad request: index '<name>' spans several graphs (<graph paths>); pass 'graph_path' to pick one` | `Bad request: index '<name>' holds several pairs (<ids>); name one by its pair_id (GET /capabilities lists them)` |
| a name with one graph and several annotations | `Bad request: index '<name>' lists one graph (<path>) with several annotations (<paths>): a traversal reads one annotation, which graph_path cannot choose; list each annotation under a name of its own` | the same new text as the row above |
| a `graph_path` on several pairs | `Bad request: graph_path '<path>' names a graph with several annotations under index '<name>' (<paths>): ...` | `Bad request: the graph_path given names several pairs of index '<name>' (<ids>); name one by its pair_id` |
| a `graph_path` not in the name | `Bad request: graph_path '<path>' is not part of index '<name>'` | `Bad request: the graph_path given is not part of index '<name>'` |
| a probe without `?graph=` | `... needs ?graph=<name> (and graph_path=<path> when the name spans several graphs); GET /capabilities lists the graphs` | `... needs ?graph=<name> (and pair_id=<id> when the name holds several pairs); GET /capabilities lists them` |
| a `pairs` object matching no pair | `request.pairs[0]: no pair of '<name>' has graph_path '<path>'` | `request.pairs[0]: no pair of '<name>' has the graph_path given` (likewise `annotation_path`, `pair_id`) |
| a `pairs` object matching several | `request.pairs[0]: '<name>' has 3 pairs on graph_path '<path>': add annotation_path` | `request.pairs[0]: '<name>' has 3 pairs on the graph_path given: add annotation_path, or name the pair by pair_id` |
| a pair not started or unfinished | `resolve: pair (<name>, <graph path>, <annotation path>) was not started ...` | `resolve: pair (<name>, <id>) was not started ...` |
| a file changed on disk | `The file '<path>' (annotation) of the pair ('<graph path>', '<annotation path>') changed on disk ...` | `The annotation file of the pair (<name>, <id>) changed on disk ...` |
| a load for a request that fails | the loader's text, which names its file | `<route>: pair (<name>, <id>) could not be loaded: ...`; the loader's text goes to the log |

The service and the suite must be told before the engines change: the service may match these texts, and the suite
pins them in its contract samples.

### 3.5 Logs

Log lines keep their paths, because an operator acts on files. Since `GET /capabilities` no longer states paths, the
server also logs one line per pair at start-up, its binding: `[Server] Pair <name> <id>: graph '<graph path>',
annotation '<annotation path>'`.

### 3.6 `/search` cuts

Each cut gains `matched`: the number of labels above the request's thresholds before the `top_labels` cut. A unit
(one query sequence on one pair) refused because its answer would pass the memory bound keeps the count its search
computed. The CLI's `--json` output gains `cuts`, with null pair tags. The staging server keeps its own members
(`nth_median`, `tied_same_median`, `refused`), its row order (score, then median abundance, then label), and `cuts`
in every mode; master omits `cuts` for `abundance_sum`. Removing members or modes would take away facts clients read.

### 3.7 Not adopted: `place`

Master's walker chooses one occurrence of the seed in each record and follows only that one; `place` states it. The
staging walker makes no such choice. Under `support: trace` it follows every coordinate chain of the record, and
under the default k-mer support it follows membership. A `place` would claim a choice the walk never made. With
`output.coordinates` true the staging server states the occurrences in `coordinates.seed[]`, up to
`output.max_coordinate_occurrences` (16 by default), with the true count when the list is cut.

### 3.8 The service

These points are the service's own view, as its maintainers stated it.

- **One meaning of `pair_id` per party.** The engine's `pair_id` is the n1 id, on its wire and in the documents the
  service stores. The service's public `pair_id` is its p2 id, live in production; renaming it would break clients.
  The service never writes a p2 id into a stored document. The join is one documented rule: **p2 = `p2-` followed by
  the first 32 hex digits of SHA-256 over the compact UTF-8 JSON `[database, n1]`**. A p2 id scopes an n1 id to one
  database. If a field turns out to be needed anyway, the traversal records of the client that holds the request's
  access token gain `native_pair_id` (n1) beside `pair_id` (p2), never anonymous search rows (section 9).
- **Binding.** The service computes the n1 id of every pair of its manifest itself (its `app/native_pairs.py`, written by
  Kulkov, which must stay byte-identical to `configured_pair_id.hpp`). On every capabilities probe and on the admin refresh it
  compares each name's set of `(graph, pair_id)` with the engine's. On a mismatch that database's capabilities become
  `unknown`: no new dispatch until the configuration is fixed, logged loudly. An answer that names an unknown id is
  refused, with or without paths.
- **Documents.** The service stores the engine's document as it comes and writes no paths into it (decision 3).
- **Selection.** The service selects by the paths of its manifest until every engine node behind an address is new.

### 3.9 The library

The library reads both forms: an old document with paths, stored before this change, and a new one with only an id.
The id form carries continuation and union from this change on. This part is done ahead of the engine, in train 7
(branch `gr/mg-t7-l1`), so that the service can vendor it before any engine changes.

- **Identity.** A document's pair id is its `pair_id`, or, for an old document without one, the n1 id computed from
  its name and paths. The library has that function, tested with master's five vectors. It reads n1 ids only and
  refuses a service p2 id.
- **Re-asking** (`pair_members`, continuations, replays of a stored request, `deepen`). A new document is re-asked
  by `pairs: [{graph, pair_id}]`. An old document is re-asked by its paths, as today; path selectors stay valid
  input on every engine. Every re-ask passes one hook: with a caller's resolver (n1 to name and paths) the id is
  resolved, checked and sent as paths, which every engine reads; without one it is sent as an id where
  `GET /capabilities` states `multigraph.pair_id: "n1"`, and refused before anything is sent otherwise.
- **Discovery** (`traverse_pairs`). It takes the pairs from `GET /capabilities` by `(graph, pair_id)` and sends
  id selectors. Record-mapping certificates (`seqs_order`) may name a pair by `(graph, pair_id)` as well as by paths.
- **Union.** The match rule keeps its order: a differing `index_meta_fp` separates, an equal `index_fp` decides,
  then the ids decide (an old document's computed id included), then the paths as today. The part key is the same
  for both forms: it holds the pair's id, stated or computed, with the paths left out wherever there is an id. An
  old document and its id-form twin therefore give one key and one digest. The cost is a one-time move: every
  union digest and union resume token changes once with this library version, a union of old documents only
  included (section 5.3).
- **Stitching** compares pairs by id, stated or computed, so an old parent and an id-form child of one pair stitch.

## 4. The contract, before and after

The examples shorten the paths and write an `index_fp` as `3f1c...`. The ids are the real ids of the two lines in
section 1.1.

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

**A `/search` row and its cut** on a multi-graph server, and a single-graph server's cut:

```jsonc
// before
{"sample": "NC_003279.8", "kmer_count": 170, "graph": "refseq",
 "graph_path": "/data/refseq/bact_other.dbg", "annotation_path": "/data/refseq/bact_other.annodbg"}
{"graph": "refseq", "graph_path": "...", "annotation_path": "...", "index_fp": "3f1c...",
 "returned": 2, "complete": false, "nth_score": 139, "tied": 3, "tied_returned": 1}
{"graph": null, "graph_path": null, "annotation_path": null, "index_fp": null, "returned": 2, "...": "..."}
// after
{"sample": "NC_003279.8", "kmer_count": 170, "graph": "refseq",
 "pair_id": "n1-64b2a5df-621c-5503-96c1-4b4b95dbd8af"}
{"graph": "refseq", "pair_id": "n1-64b2a5df-621c-5503-96c1-4b4b95dbd8af", "index_fp": "3f1c...",
 "matched": 5, "returned": 2, "complete": false, "nth_score": 139, "tied": 3, "tied_returned": 1}
{"graph": null, "pair_id": null, "index_fp": null, "matched": 5, "returned": 2, "...": "..."}
```

**A selector**: both name the first pair; only the second is new.

```json
{"pairs": [{"graph": "refseq", "graph_path": "/data/refseq/bact_other.dbg"}]}
{"pairs": [{"graph": "refseq", "pair_id": "n1-64b2a5df-621c-5503-96c1-4b4b95dbd8af"}]}
```

`GET /traverse/capabilities?graph=refseq&pair_id=n1-64b2a5df-621c-5503-96c1-4b4b95dbd8af` probes it.

**A refusal and its candidates** (`/traverse` on a name of two pairs):

```jsonc
// before
{"error": "Bad request: index 'refseq' spans several graphs (/data/refseq/bact_other.dbg, /data/refseq/bact_pseudo.dbg); pass 'graph_path' to pick one",
 "code": "traverse_single_pair", "candidates": [{"graph": "refseq", "graph_path": "...", "annotation_path": "..."}, "..."]}
// after
{"error": "Bad request: index 'refseq' holds several pairs (n1-64b2a5df-..., n1-5b5559a9-...); name one by its pair_id (GET /capabilities lists them)",
 "code": "traverse_single_pair", "candidates": [{"graph": "refseq", "pair_id": "n1-64b2a5df-..."}, "..."]}
```

**`GET /capabilities`:**

```jsonc
// before
"graph_summary": {"pairs": [{"graph": "refseq", "graph_path": "/data/refseq/bact_other.dbg",
                             "annotation_path": "/data/refseq/bact_other.annodbg", "index_fp": "3f1c...", "k": 31}]}
// after
"graph_summary": {"pairs": [{"graph": "refseq", "pair_id": "n1-64b2a5df-621c-5503-96c1-4b4b95dbd8af",
                             "index_fp": "3f1c...", "k": 31}]},
"multigraph": {"pair_id": "n1", "selectors": {"/traverse/capabilities": ["graph", "graph_path",
               "annotation_path", "pair_id"], "...": "..."}}
```

Unchanged: `/align` (a multi-graph server refuses it), `index_fp` and `composite_fp` (content identity), the
capability profile's `support_fp` (it excludes identity).

## 5. Compatibility and migration

### 5.1 Rollout order

1. **The service's adapter.** It reads both forms (a current engine states paths; a new one states `pair_id` and
   `multigraph.pair_id`), binds ids, runs the two-way check and keeps path selectors.
2. **The library**, vendored by the service, reading both forms, so that the service can continue and union the
   id-only documents a new engine gives. Its code hash (`library_sha256`, which the service reports) moves. This
   step is done in train 7, ahead of the engines.
3. **The suite**, accepting both forms in its pins.
4. **The engines**, on staging and then in production.

Path selectors stay in use until every node behind an address is new: a current engine refuses `pair_id` in a
`pairs` object as an unknown member. The service believes each production path reaches one engine process; section 9
asks for that to be confirmed. As a second guard, `expect_server_instance` and `usage.server_instance` let the
service refuse an answer from an unexpected instance. A rollback to a current engine brings the paths back, and the
adapter reads that form too.

### 5.2 Old stored documents and caches

Old documents have paths and no `pair_id`; they stay as stored. The library computes their id from the paths, so an
old and a new document of the same configured pair have the same id. The service's caches keyed on paths either map
an answer's id to its manifest pair before the lookup, or move to `(graph, pair_id)` keys. Replacing a file's contents
keeps both its paths and its id, so a cache keyed on either must also check `index_fp`. The engine's `/search`
result cache lives in the process, and the restart that deploys it clears it.

### 5.3 Digests and tokens

- Unchanged: a graphlet's body digest and the resume tokens of single-graphlet operations. They cover the records
  after the J line, so no tag is an input.
- Changed once, with the library version that reads both forms: every union digest and union resume token, and the
  `~` tags of names with several pairs, because the part key now holds the pair's id and the expected pairs are
  hashed as keys (section 3.9). A union of old documents only moves too: production RefSeq's union of old documents
  moved from `9e369aca...` to `60531fdaf...`. After that move the form makes no difference: an old document and its
  id-form twin give the same key, digest and tokens, so the service's switch to id-form documents moves nothing.
- The library's code hash moves with its code.

### 5.4 What changes, for whom

| Who | Change | Remedy |
| --- | --- | --- |
| Clients reading `graph_path` or `annotation_path` from any response, `GET /capabilities` included | the members are gone | read `pair_id`; bind it through their own configuration |
| Clients matching selection texts whole | the texts change (section 3.4) | match the code; take the new wording |
| Single-graph clients | `/search` cuts state `pair_id: null` instead of two null paths and gain `matched`; the registry's and `pair_summary`'s null tags likewise; the CLI's `--json` gains `cuts` | update strict key sets (this repository's `test_api.py`, the contract samples, the suite's cut keys) |
| Readers of the service's downloadable documents | new documents carry `pair_id`, no paths; old ones keep theirs | read `pair_id`; the library reads both |
| The service | binds ids from its manifest; two-way check; stores documents as they come; keeps p2 public | section 3.8 |
| The library | reads ids, re-asks by id, discovers by id; every union digest changes once (section 5.3) | section 3.9 |
| The suite | pinned key sets, texts and samples; it can no longer match graph-list rows to `GET /capabilities` by paths | accept both forms; match rows by their computed id |
| Requests near the `/search` memory bound | a row's tags get shorter (39 bytes instead of two paths), so the answer's price drops | none needed |

## 6. Alternatives

| Design | Why not |
| --- | --- |
| The engine adds `pair_id` and keeps the paths | contradicts master's form and the owner's decision 2 |
| Paths kept in configuration statements (`GET /capabilities`, the selection texts) as documented exceptions, or written back by the service into its stored documents | the owner decided that no response and no stored document states a path; either would have kept the library unchanged |
| The service renames its public `pair_id`, or writes p2 into stored documents | its p2 is live in production; the join rule makes a second field unnecessary |
| A server flag that keeps paths for one deployment | two answer forms to test and document, and the fleet would still state paths |

## 7. Risks and what the review should check

- **Encoding stability.** Three implementations must stay byte-identical: the engine's, the service's and the
  library's. Check that all three run the five vectors of `tests/common/test_configured_pair_id.cpp`, and that
  `n1-` versions the scheme, so that a change gets a new prefix.
- **Paths spelled differently.** The id hashes the strings, not the files. A service manifest that spells a path
  through a symlink, a relative path or a trailing slash computes another id. The two-way check must turn that into
  the loud `unknown` state, and the start-up log line lets an operator see the engine's spelling.
- **Aliases.** Renaming a graph name changes its pairs' ids; a pair under two names has two ids but one union
  identity when manifests are present. Check that nothing keys long-lived data on an id across a rename.
- **Mixed forms.** An old document (paths) and a new one (id) of one pair must have one identity in union and
  stitch. Check it with documents of both forms.
- **Ids are public, not secret.** An id is deterministic from a name and two paths, so it confirms a guessed path
  offline. Check that this is acceptable where answers leave the service.
- **The old engine and the new selector.** Check that the service stays on path selectors until every node is new.
- **Answer prices.** The `/search` answer budget prices tags by their length. Check that the price model and its
  tests follow the new member set.
- **Divergence from master.** The staging server adds a `pair_id` selector, keeps its own `cuts` members, rewords its
  texts and leaves out `place`. Check that none of these blocks a later merge from master.

## 8. Plan, tests and cost

| Step | Work | Days |
| --- | --- | --- |
| 1 | Merge master's two commits: the id header and its test taken as they are; every conflicting hunk resolved to the staging code, whose `cuts` replace master's | 0.5 |
| 2 | The engine's pair identity: tags, registry, candidates, `GET /capabilities` without paths, the selector and its conflict rules, answer prices | 2 |
| 3 | The engine's texts without paths (section 3.4), loader texts to the log, the start-up binding line | 0.75 |
| 4 | `matched`, and `cuts` in the CLI's JSON output | 0.5 |
| 5 | Re-record the contract samples (29 of 52 state paths; the texts and `GET /capabilities` change; the `/search` ones gain `matched`) | 1 |
| 6 | Specifications: the envelope spec, the core traversal spec, the pattern spec, the manual pages | 0.75 |
| 7 | The library: done in train 7, ahead of the engine (section 3.9). What is left here: its tests against this engine's answers (an id-only `GET /capabilities` with `multigraph.pair_id: "n1"`, a name of several pairs), once train 7 is merged into this branch | 0.25 |
| 8 | Verification: all build targets for DNA and protein, unit and integration tests, Python tests, route and CLI comparisons against the current build, memory and timing, a GCC check, CI | 1.25 |
| 9 | Review and fixes | 1 |

About eight working days in this repository, the library's step done separately. The service's own work (adapter,
binding, two-way check) is not counted here; it no longer writes paths into documents. Tests to add beyond the
updated ones:

- Engine: the ids equal Python's `uuid5` for the pairs of a test graph list, an alias and a non-ASCII path included;
  **no response contains a configured path**, swept over every route, probe and refusal on single- and multi-graph
  servers (selection, not reached, unfinished, a changed file, a failed load); an id selector names the pair its
  paths name; `pair_id` with a path in one `pairs` object or one probe is a plain 400; single-graph cut keys and
  `matched` with no match, `top_labels: 0`, count ties, `abundance_sum`, and a unit refused by either account.
- Library (train 7): an old and a new document of one pair have one id and one stitch identity, and are one pair
  in a union with one key and one digest; a new document is re-asked by id, an old one by its paths; discovery on an
  id-only `GET /capabilities` with a name of several pairs; two aliases with manifests are one union pair. After the
  engine: the same against this engine's answers.
- Service (in its repository): stored J lines carry the n1 id and no paths, and p2 agrees with search answers by the
  join rule; continuation, union and replay of old and new documents; an extra engine pair, a missing one and a
  misspelled path each make the database `unknown`; an answer with an unknown id is refused; a content change (new
  `index_fp`, same paths) invalidates the caches; current and new capabilities in both orders, and an answer from an
  unexpected instance refused.

**Timing.** Staging changes are deployed in batches called trains. The engine's part started on 2026-11-19, in
parallel with train 7's last fixes, on top of train 7 with its change to how `/resolve` holds its answer as text. It
lands in train 8. Train 7 deploys first, with the library's part and the outside review's fixes. The service's
adapter goes before the engines (section 5.1).

## 9. Open points

Both points gate the rollout of the engines, not their implementation.

1. **Kulkov** (the service's pair identity is his): may the service's search rows state the n1 id? Recommended: no;
   they keep the p2 id with the join rule of section 3.8, and `native_pair_id` goes into the token holder's traversal
   records (section 3.8) only if a client needs the field.
2. **Kulkov or the operator of the production proxy:** does each production path reach one engine process?
   Recommended: confirm before the engines roll out; until every node behind an address is new, the service keeps
   path selectors and checks the answering instance.
