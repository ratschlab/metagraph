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
touching the annotation. **As a job** — submit, status, results — the way the service wraps `/search`: the
owner decided (2026-10-07) that the backend route stays synchronous and the service serves it asynchronously,
because that is what the service is for; there is no synchronous MCP tool for it.

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
`/traverse`. Synchronous request/response, compact JSON, gzip when accepted; its own deadline (default 60 s, cap
600 s, §3.2) with a finalisation reserve; the 900 s content timeout above the cap.

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
| 7 | Python client methods; this service's job type | with you; on refseq33m-experimental after backend milestone 1, on chunked databases after milestone 6 (§3.1 item 3) |
| mask | refseq33m-experimental's graph has no `.edgemask` file, and the route needs one (DESIGN §4): before the route answers on staging the owner runs `metagraph transform --mask-dummy` once on mex (offline, no change to node ids or annotation), or the server starts with `--pattern-build-mask`; until then the host states `mask: absent` and the job answers `mask_required` | before the owner's `update.sh` that enables the route |
| fixtures | with milestone 1's freeze commit, as for level 6: the capabilities block on both routes, one answer per mode and per `withheld` reason, an error slot, from the mini index, under `api/python/tests/data/traverse/pattern/`, so your unit tests do not wait for a host | with milestone 1 |

Staging (`refseq33m-experimental`, a BASIC index with record mapping) gets the route when the owner runs
`update.sh` after milestone 1 is pushed. Until then there is no host to test against; the contract and the
capabilities block are what to build on.

## 3. Requirements and proposals

### 3.1 The job (the service surface)

1. **(required)** A **kind of the existing search job** (the owner, 2026-10-07: "it is essentially a search, and
   the results are like search too"): the same Search/Task/Result tables, the same leaf split (one task per
   (database, graph chunk, pattern chunk of at most the host's `max_patterns`)), the same semaphore, per-database
   caps, 1,200 s client timeout, status, results, CSV, S3, lock and retention. The worker calls `POST /pattern`
   instead of `POST /search`, skips scoring, truncation and enrichment for that kind, and stores contexts as
   result rows with the per-pattern summary per task: a kind marker and a summary column, no new worker family,
   no new tables. No synchronous MCP tool.
2. **(required)** The merged view of a job: per count the **weakest relation wins** (`exact` only if every task's
   is); `retrieval_complete` only when every task's is true; `withheld` and `stop` carried per task with the
   graph and its reason; rows are contexts, never labels.
3. **(required)** A chunked database is many graphs on one multi-graph server process, selected per task through
   `graphs: ["{label}-{i}/{N}"]` as `/search` does (`app/download_depth.py`, `enumerate_leaf_specs`). The job type
   is built and tested on refseq33m-experimental (one graph) with backend milestone 1; serving the chunked
   databases waits for backend milestone 6, which keeps `graphs` exactly as `/search` selects a shard.
4. **(required)** Admission by the queue, **shared with search and traversal** (the owner: "traverse, pattern
   match and normal search need to share the pool inside async and inside the sync server"): the per-database
   queues, `META_DB_CAPS` and the distributed semaphore count search, pattern and traversal calls against one cap
   per database server (the cap search uses today, e.g. 15 on sra-logan-chunks). The server has one request pool
   (`-p` / `--threads-each`) for every route and reserves nothing per route; the engine needs nothing beyond
   "at most cap concurrent calls per database server across all routes", and its memory is bounded per request
   by `max_memory_mb`. The traversal jobs' separate per-host tokens move into the same semaphore (your
   follow-up).
5. **(required)** Caps at submit, never truncation afterwards: `max_contexts`, `max_patterns`, `time_budget_ms`
   above the service's ceilings are refused with a 400 naming the field (the strategy validator's rule: refuses,
   never lowers). An answer is passed through whole; cutting it on the way out would falsify
   `retrieval_complete`. The task timeout derives from the request's `time_budget_ms` plus an allowance.
6. **(required)** Served only for databases whose host carries the `pattern` block (§1) with a contract version
   the service knows or a higher one; a lower or missing version, or no block, answers `pattern_unsupported` with
   the host's feature list, as `traversal_disabled` does. Label projections are offered exactly when the host's
   `projections` list has them. The probe that reads capabilities already exists (`app/traversal/probe.py`).
7. **(required)** The answer passed through with the service's additions only: `database`, the untrusted-data
   notice on label and record strings, the standard error envelope. **Never** sum per-shard label counts into one
   number, never drop `relation`, `withheld`, `retrieval_complete` or `absence_scope`, never add labels of its own.
   `hit_unit` does not apply to a context count; the tool says what a context is instead: one distinct k-mer of
   the index that contains the pattern, at one offset and orientation — not a record, not an occurrence, not a
   hit; a k-mer present in ten samples is one context.
8. **(required)** The submit and results tools' docstrings explain the three things an agent acts on:
   `withheld: count_above_threshold` (narrow the pattern, scope or strand), `withheld: discovery_budget` (a more
   informative pattern; a filter does not help), `withheld: annotation_budget` (ask for the count or
   `labels: none`, or a narrower pattern); that `count` first, then `labels: none`, is the cheap way to look at a
   new pattern; they use the service's existing strand vocabulary, and say that canonical and primary indexes
   report contexts only — no strand, no position.
9. A stored answer keeps its `determinism`; a `time_limited` one is marked as not reproducible where the
   service shows it.
10. Rate limits and the anonymous budget as for the other jobs. Patterns are query data under the same privacy
    policy as sequences.
11. Row ids are opaque and valid per (host, index release, graph); the later backend request "labels for given
    rows" takes them back. The job returns them as given.

### 3.2 Time budget (decided 2026-10-07)

The owner: "5 s is not sufficient in general" — some `/search` calls take longer today (Logan hosts about 28 s),
and pattern search is the more complex function. So the route's **default** `time_budget_ms` is **60 s** (what
the CLI and a direct caller get) and the server cap `--pattern-max-time-ms` is **600 s**: under the 900 s content
timeout, with room for serialising and compressing a large answer. Budgets exactly like search: the service
sends `time_budget_ms` equal to the host's cap (600 s) unless the caller asks for less, and its 1,200 s client
timeout (`META_CALL_TIMEOUT_SECS`) stays above it. The deadline keeps its value under a job: it bounds a runaway
`any_offset` DFS per task, and the finalisation reserve is unchanged, so a stopped task still answers with
counts. Nothing in the engine assumes a short run: the clock is read every 4,096 steps whatever the budget; the
memory account is bounded by `max_memory_mb` and the retained descriptors by their caps, not by time;
`max_steps` (default 10⁸ per shard) is the other bound and is raised with the budget. The cost of a long call is
the request-pool slot it holds, counted by the shared per-database cap of item 4.

## 4. Constraints

- Never claim absence from a `count` or `labels: none` answer; show `absence_scope` and `retrieval_complete`.
- A `determinism: time_limited` answer is not reproducible; say so where the service caches answers.
- The service does not retry a `withheld` answer with a larger cap on its own; the agent decides.
- Nothing here touches the search path or the traversal jobs.

## 5. Tests and acceptance

- Unit: capabilities gating (block present, absent, lower, higher contract version); the per-task call and its
  timeout; the merge (weakest relation, `retrieval_complete`, per-task `withheld`); the pass-through keeps every
  field; the error envelope; the docstring's tool table regenerated. The fixture bodies of §2 feed these.
- e2e against staging once the route is deployed: a 16-nt primer of a known refseq33m record in `count` and in
  `labels: none` (the returned k-mers contain the primer), an IUPAC promoter, a palindrome, a pattern below the
  information floor (refused), `max_contexts` 1 in `all_or_count` (withheld with the exact count) and in
  `partial` (one context, the cut stated), an unsupported host, a two-host database merged.
- Acceptance: a primer job done within the service's usual job latency on staging; every `withheld` reason
  surfaced verbatim per task; no label claim from a count answer in the tools' own wording.
- A new tool moves the registry count that the docs test and the e2e skill pin; update both in the same change.

## 6. Decisions (the service side's proposals of 2026-10-07, accepted by the MetaGraph side)

1. Surface: a **job** (submit → status → results), the owner's decision; no synchronous tool. The backend route
   is `POST /pattern`; the service's job kind and tool names are the service's (the `pattern_` family, so the
   later "labels for given rows" joins it). Not `/motif`.
2. Defaults of the job: `mode: count`, `scope: any_offset` (the complete scope), with `suffix` offered as the
   cheap option and its `absence_scope` shown. `partial` is exposed: it is explicit, every cut is stated, and an
   agent that wants "the first fifty" has no other path. (The route's own default mode is `all_or_count`.)
3. Databases: every host that carries the block, canonical indexes included (contexts only). The catalogue gets
   a `pattern` field beside `traversal`, from the same probe, with the same filter.
4. Wording: "one distinct k-mer of the index that contains the pattern, at one offset and orientation. Not a
   record, not an occurrence, not a hit: a k-mer present in ten samples is one context."
5. The job is a kind of the existing search job (same tables, fan-out, semaphore, caps, timeout, status,
   results, S3, lock, retention), the owner's decision; refseq33m-experimental after backend milestone 1, chunked
   databases after milestone 6 (§3.1 items 1 and 3).
6. Who and when: the service session writes the design note and implements through Opus agents; org-id
   integrates; the MetaGraph side supplies the frozen contract, the fixture bodies (§2) and the `pattern` block
   on `/traverse/capabilities` with milestone 1's commit. First increment after that commit; live e2e after the
   owner's `update.sh`.

Open for the owner: the design note's approval and the start date.
