# Feature request: pattern search in metagraph-search-service

**Status:** request, 2026-10-07. The MetaGraph side is in implementation (milestone 1 running); the wire contract
of its first milestone freezes within a day or two. This note is for scoping and the design note on the service
side; it is not a go to implement before the owner approves that note.

You are working in `~/git/services/metagraph-search-service`, the FastAPI + FastMCP service in front of the
MetaGraph `server_query` backends. Read its `CLAUDE.md` first and follow it (two Redis instances, registry changes
to Supabase and Railway, the README tool table from `scripts/gen_mcp_tool_table.py`, new tests in the CI list).

**The task.** Offer MetaGraph's new **pattern search** to agents: patterns shorter than k (a 16-nt primer, a
20-nt spacer), IUPAC patterns (a promoter written with R, Y, N), and later peptides, found **completely** in an
index and **counted before anything is read**, with a label-free path that returns the matching k-mers without
touching the annotation. First as a **synchronous tool**, like `traverse_resolve`. A queue job type for batches
and archive-wide scans is a later step and is scoped at the end.

**Process.**
1. Work out the design with the owner.
2. Record it in a design note (`dev-notes/2026-10/pattern-search/DESIGN.md`) answering §6. No implementation
   before the owner approves it.
3. Implement in increments, each with tests, gated on the backend's milestones (§2).

In §3, items marked **(required)** are the owner's decisions or follow from the MetaGraph guarantees. The rest
is a proposal from the MetaGraph side to weigh, adopt, change or reject.

## 1. What exists on the MetaGraph side

All on branch `gr/labeled-traversal` of `ratschlab/metagraph`. The normative document is
`~/git/services/metagraph/metagraph/docs/DESIGN-pattern-search.md` (v6, approved for phased implementation after six
external reviews). Read: §2 goals, §3 semantics (scope, graph context vs occurrence, counts with relations,
strand), §5.2 the count-first modes and the `withheld` reasons, §4.3 "The label-free path first", §5.6 the
annotation predicates (later), §7 the request and answer, §7.3 the capabilities block, §9 the service layer as
proposed, §13 the increments.

**The route:** `POST /pattern` on the same `server_query` process that serves `/search`, `/resolve` and
`/traverse`. Synchronous request/response, compact JSON, gzip when accepted; its own deadline (default 5 s) with
a finalisation reserve; the 900 s content timeout far above it.

**What a request says.** A list of patterns (`dna` | `iupac`, later `protein`), a `mode`
(`count` | `all_or_count` | `partial`), a `scope` (`suffix` | `any_offset`; `long` is implied for patterns longer
than k), `strands`, caps (`max_contexts`, `max_steps`, `time_budget_ms`, …), and an `output.labels` projection:
`none` (the label-free path), later `all` and `predicate_only`.

**What an answer says.** Per pattern: every count as `{value, relation, unit}` with
`relation ∈ {exact, at_least, bounds, unknown}` and `unit ∈ {graph_contexts, anchors, paths, placed_occurrences,
labels}`; `work` (ranges visited, steps); `stop`; `withheld` with its reason when results are not returned;
`retrieval_complete`; `absence_scope`; `determinism`; `notes`; and `results`, one entry per **graph context**
(a k-mer, or a path for long patterns) with `kmer`, `instance`, `offset`, `strand`, `node`, `row` and, when
labels were asked for, `labels` with placed occurrences and per-label `support`.

**Two facts a tool must not blur.** A graph context is one k-mer of the graph, not a physical occurrence: a
k-mer present in ten records of a column is one context. And a `count` answer, or a `labels: none` answer, says
nothing about any label; absence of a label is claimed only by an answer with `retrieval_complete: true` and no
filter.

**Capabilities.** The `pattern` block is carried by **both** `GET /capabilities` (under `features` and
`routes`, as the other features) and `GET /traverse/capabilities` — the document your probe already reads, which
today carries `attempts`, `coordinates` and `deadline_check` and no feature list — so one cached probe serves
both. The block: `modes`, `projections` (the list the host offers **now**: `["none"]` at milestone 1, `all` and
`predicate_only` once milestone 3 lands; gate the label projections on this list, never on a milestone number),
scopes per graph mode, the caps and floors, `placement` and `support` the index can give, `mask`,
`annotation: budgeted | unbudgeted`, `pattern_contract_version` (accept a higher version and read fields by
presence, as you did for feature level 6; refuse only a lower or a missing one). A host without the block has
no route.

## 2. When

| backend milestone | content | state |
|---|---|---|
| 1 | count (`mode: count`) and the label-free extraction (`all_or_count` / `partial` with `labels: none`): k-mers, offsets, strands, node and row ids; exact DNA and IUPAC; both strands; `suffix` and `any_offset`; single-graph servers; the capabilities block; `metagraph pattern` CLI | running now; contract freezes on its commit |
| 3 | `labels: all`: label discovery and placement (record, 1-based position, strand) on BASIC indexes with record mapping | next |
| 4 | patterns longer than k (extension), per-label `support`, `require_support` | after 3 |
| 5 / 5b | peptides (codon automaton); annotation predicates (`any`, `all`, `none`, `at_least`, `and`/`or`/`not`) | after 4 |
| 6 | multi-graph servers (per-shard budgets, barriers, shard identity per result), the real-index benchmark | after 5 |
| 7 | Python client methods; this service's tool | with you |
| fixtures | with milestone 1's freeze commit, as for level 6: the capabilities block on both routes, one answer per mode and per `withheld` reason, an error slot, from the mini index, under `api/python/tests/data/traverse/pattern/`, so your unit tests do not wait for a host | with milestone 1 |

Staging (`refseq33m-experimental`, a BASIC index with record mapping) gets the route when the owner runs
`update.sh` after milestone 1 is pushed. Until then there is no host to test against; the contract and the
capabilities block are what to build on.

## 3. Requirements and proposals

### 3.1 The synchronous tool (first)

1. **(required)** An MCP tool `pattern_search(database, dna | iupac, mode, scope, strands, max_contexts, …)` and a
   REST route `POST /pattern/search`, built like `traverse_resolve` (`app/mcp_server.py`, `app/api/routes_traversal.py`
   with its `_admit`): synchronous, no job, no S3; the backend called from the API process on the thread pool.
2. **(required)** Admission in `traverse_resolve`'s shape — a leased token per backend host with its own
   concurrency, lease and timeout settings (`PATTERN_MAX_CONCURRENCY`, `PATTERN_PER_HOST`, `PATTERN_LEASE_S`,
   `PATTERN_TIMEOUT_S`, `PATTERN_MAX_PATTERNS`) — not the in-process `_COMPARE_SLOTS` semaphore, which guards a
   local computation; this call hits a shared host. Staging and production counting separately is the known gap,
   bounded well enough for v1 by the backend's own 5 s deadline.
3. **(required)** Caps at admission, never truncation afterwards: `max_contexts`, `max_patterns` and
   `time_budget_ms` above the service's ceilings are refused with a 400 naming the field (the strategy validator's
   rule: refuses, never lowers). An answer is passed through whole; cutting it on the way out would falsify
   `retrieval_complete`. The HTTP timeout derives from the request's `time_budget_ms` plus an allowance, not from
   the host's maximum.
4. **(required)** Served only for databases whose host carries the `pattern` block (§1) with a contract version
   the service knows or a higher one; a lower or missing version, or no block, answers `pattern_unsupported` with
   the host's feature list, as `traversal_disabled` does. Label projections are offered exactly when the host's
   `projections` list has them. The probe that reads capabilities already exists (`app/traversal/probe.py`).
5. **(required)** The answer passed through with the service's additions only: `database`, the untrusted-data
   notice on label and record strings, the standard error envelope. **Never** sum per-shard label counts into one
   number, never drop `relation`, `withheld`, `retrieval_complete` or `absence_scope`, never add labels of its own.
   `hit_unit` does not apply to a context count; the tool says what a context is instead: one distinct k-mer of
   the index that contains the pattern, at one offset and orientation — not a record, not an occurrence, not a
   hit; a k-mer present in ten samples is one context.
8. **(required)** The docstring explains the three things an agent acts on: `withheld: count_above_threshold`
   (narrow the pattern, scope or strand), `withheld: discovery_budget` (a more informative pattern; a filter does
   not help), `withheld: annotation_budget` (ask for the count or `labels: none`, or a narrower pattern); that
   `count` first, then `labels: none`, is the cheap way to look at a new pattern; it uses the service's existing
   strand vocabulary, and says that canonical and primary indexes report contexts only — no strand, no position.
9. No caching of this tool's answers in v1; `determinism: time_limited` matters only once an answer is cached, and
   the design note states that.
6. Rate limits and the anonymous-call budget as for the other synchronous tools. Patterns are query data under the
   same privacy policy as sequences.
7. Row ids are opaque and valid per (host, index release, graph); the later backend request "labels for given
   rows" takes them back. The tool returns them as given.

### 3.2 Later: a queue job type

The owner asked whether the service should support "this extra type of job". The MetaGraph design says: not
first. The synchronous tool covers the interactive case, and the route answers in seconds under its own deadline.
A queue job type is the right shape for two things only: **batches** (hundreds of patterns, e.g. a primer panel)
and **archive-wide scans** across every database of a catalogue entry. When it comes, it wraps the same route:
one task per (database, host, pattern chunk), per-task completion recorded as search tasks record it, results
merged with the shard identity each result carries, large answers to S3 as searches do. What it needs from the
service that the sync tool does not: the job type registry, the merge of `withheld` and relations across tasks
(the weakest relation wins; `retrieval_complete` only when every task's is true), and a result table that keeps
contexts, not labels, as rows. Scope it in the design note; build it after milestone 6 on the backend side, when
multi-graph hosts exist.

## 4. Constraints

- Never claim absence from a `count` or `labels: none` answer; show `absence_scope` and `retrieval_complete`.
- A `determinism: time_limited` answer is not reproducible; say so where the service caches answers.
- The service does not retry a `withheld` answer with a larger cap on its own; the agent decides.
- Nothing here touches the search path or the traversal jobs.

## 5. Tests and acceptance

- Unit: capabilities gating (feature present, absent, unknown contract version); admission and timeouts; the
  pass-through keeps every field; the error envelope; the docstring's tool table regenerated.
- e2e against staging once the route is deployed: a 16-nt primer of a known refseq33m record in `count` and in
  `labels: none` (the returned k-mers contain the primer), an IUPAC promoter, a palindrome, a pattern below the
  information floor (refused), `max_contexts` 1 in `all_or_count` (withheld with the exact count) and in
  `partial` (one context, the cut stated), an unsupported host.
- Acceptance: a primer answered within 2 s on staging; every `withheld` reason surfaced verbatim; no label claim
  from a count answer in the tool's own wording.
- A new tool moves the registry count that the docs test and the e2e skill pin; update both in the same change.

## 6. Decisions (the service side's proposals of 2026-10-07, accepted by the MetaGraph side)

1. Names: `pattern_search` and `POST /pattern/search`; `pattern_` is the family name, so the later "labels for
   given rows" tool joins it. Not `/motif`.
2. Defaults: `mode: count`, `scope: any_offset` (the complete scope), with `suffix` offered as the cheap option
   and its `absence_scope` shown. `partial` is exposed: it is explicit, every cut is stated, and an agent that
   wants "the first fifty" has no other path. (The route's own default mode is `all_or_count`; the tool's is
   `count`, and the two documents now say so.)
3. Databases: every host that carries the block, canonical indexes included (contexts only). The catalogue gets
   a `pattern` field beside `traversal`, from the same probe, with the same filter.
4. Wording: "one distinct k-mer of the index that contains the pattern, at one offset and orientation. Not a
   record, not an occurrence, not a hit: a k-mer present in ten samples is one context."
5. Queue job: not first. When it comes, modelled on the traversal jobs (one document per task), not on the
   search tables, whose hit-shaped truncation and enrichment make no sense for contexts; after backend
   milestone 6.
6. Who and when: the service session writes the design note and implements through Opus agents; org-id
   integrates; the MetaGraph side supplies the frozen contract plus the fixture bodies (§2) with milestone 1's
   commit. First increment after that commit; live e2e after the owner's `update.sh`.

Still open for the owner: nothing on the MetaGraph side; the design note's approval and the start date are his.
