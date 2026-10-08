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

All on branch `gr/labeled-traversal` of `ratschlab/metagraph`. The wire contract you build on is
`~/git/services/metagraph/metagraph/docs/SPEC-pattern-search.md` (contract version 1, frozen; its §16 lists what
the review of 2026-10-07 changed). The normative design is
`~/git/services/metagraph/metagraph/docs/DESIGN-pattern-search.md` (v6, approved for phased implementation after six
external reviews). Read: §2 goals, §3 semantics (scope, graph context vs occurrence, counts with relations,
strand), §5.2 the count-first modes and the `withheld` reasons, §4.3 "The label-free path first", §5.6 the
annotation predicates (later), §7 the request and answer, §7.3 the capabilities block, §9 the service layer as
proposed, §13 the increments.

**The route:** `POST /pattern` on the same `server_query` process that serves `/search`, `/resolve` and
`/traverse`. Synchronous request/response, compact JSON, gzip when accepted; its own deadline (default 60 s, cap
600 s, §3.2) with a finalisation reserve; the 900 s content timeout above the cap.

**What a request says.** A list of patterns (`dna` | `iupac` | `protein`, the last since increment 5), a `mode`
(`count` | `all_or_count` | `partial`), a `scope` (`suffix` | `any_offset`; `long` is implied for patterns longer
than k), `strands`, caps (`max_contexts`, `max_steps`, `time_budget_ms`, …), and an `output.labels` projection:
`none` (the label-free path), `all` (since milestone 3), later `predicate_only` (5b).

**What an answer says.** Per pattern: every count as `{value, relation, unit}` with
`relation ∈ {exact, at_least, bounds, unknown}` and `unit ∈ {graph_contexts, anchors, paths, placed_occurrences,
labels}`; `work` (ranges visited, steps); `stop`; `withheld` with its reason when results are not returned;
`retrieval_complete`; `absence_scope`; `determinism`; `notes`; and `results`, one entry per **graph context**
(a k-mer, or a path for long patterns) with `kmer`, `instance`, `offset`, `strand` on BASIC hosts
(`index.strand_stated: true`) or `orientation` (forward / reverse / palindromic) on canonical and primary hosts
(counts `by_strand` or `by_orientation` likewise), `node`, `row` and, when labels were asked for, `labels` with
placed occurrences and per-label `support`.

**Two facts a tool must not blur.** A graph context is one k-mer of the graph, not a physical occurrence: a
k-mer present in ten records of a column is one context. And a `count` answer, or a `labels: none` answer, says
nothing about any label; absence of a label is claimed only by an answer with `retrieval_complete: true` and no
filter.

**Capabilities.** The `pattern` block is carried by **both** `GET /capabilities` (under `features` and
`routes`, as the other features) and `GET /traverse/capabilities` — the document your probe already reads, which
today carries `attempts`, `coordinates` and `deadline_check` and no feature list — so one cached probe serves
both. The block: `modes`, `projections` (the list the host offers **now**: `["none"]` at milestone 1,
`["none", "all"]` since milestone 3, `predicate_only` with 5b; gate the label projections on this list, never on
a milestone number),
scopes per graph mode, the caps and floors, `placement` and `support` the index can give, `mask`,
`annotation: budgeted | unbudgeted`, `pattern_contract_version` (accept only a version the service implements,
1 today; refuse a missing, malformed or other one, a higher one included, since a higher version means a field
changed meaning; within version 1 read fields by presence and pass unknown values of the extensible
enumerations through, SPEC §1). A host without the block has no route.

**Increments 4 and 5 (in the build since 2026-10-08; SPEC §12.1, §12.2, §17).** Two additions to contract
version 1, both opt-in, each gated on the capabilities block, never on a milestone number:
- **Paths of a pattern longer than k**: offer them only where `long_search` lists `"paths"`; the job then sends
  `long_search: "paths"` (and `max_paths` within `caps.max_paths`). Without it a long pattern keeps the
  anchor-only answer (`withheld: paths_later_increment`). A path result has `sequence` (the L bases),
  `anchor_kmer`, `nodes` and `rows`, never `kmer`; `counts.paths` is known with its `extension`. With labels, each
  label of a path states its `support`: `record_verified` (one record holds the whole path, with its placed
  occurrences) or `label_intersection` (every k-mer of the path carries it, no record claim). Show the two apart;
  a tool that reports "found in record X" uses `record_verified` labels only. `require_support:
  "record_verified"` lists the verified labels only (each path counts the others in
  `labels_excluded_unverified`, `null` when that number is not known: never read a `null` as 0) and is allowed
  only where the block's `support` is `record_verified` (else 400 `support_unavailable`). New values to handle: `withheld` `anchors_above_threshold`,
  `cut` `max_paths`, `stop.phase` `extension`, `stop.reason` `max_paths`, note `label_intersection_only`.
- **Peptides**: offer them only where `kinds` lists `"protein"`: a `protein` pattern over `protein_residues` (the
  20 amino acids and X, B, Z, J, and the stop `*` where the list has it), in `genetic_code` (one of
  `genetic_codes`, default `default_genetic_code`, 1; another integer is 400 `genetic_code_unknown`). The entry
  states `residues` and `genetic_code`, and `length` stays in bases (3 per residue), so a peptide of more than
  k / 3 residues is a long pattern (paths as above).

**The owner's decisions of 2026-10-08 (in the build; SPEC §18).** Additions to contract version 1:
- **Graphs without the dummy-edge mask are served.** Gate on the block's `counting`: `"exact"` (a mask) or
  `"upper_bound"` (none; the block then gives `dummy_fraction`, the sampled fraction of real k-mers with its 95%
  interval). On an `upper_bound` host a count is `exact` where the backend can prove it and otherwise `bounds`
  [`lower`, `upper`] with an `estimate` (round(upper × f)): show the estimate **as an estimate**, beside its bounds
  and the host's `dummy_fraction`, never as the count and never as a bound (a pattern at a record start can sit far
  below it). The merge algebra of §3.1 item 2 adds `lower` and `upper` and leaves estimates out; a merged estimate,
  if the service shows one, is labelled as such. The host admits on the upper bound (`all_or_count` can withhold a
  pattern whose true count fits, note `threshold_upper_bound`; `partial` lists it); its lists are exact. New notes:
  `estimate_sampled_dummy_fraction`, `threshold_upper_bound`, `no_stop_codon`. `mask_required` is retired: no
  current build answers it (an older one may: keep passing it through).
- **The stop `*`** is a residue where `protein_residues` lists it: a stop codon of `genetic_code` at that
  position (X never matches one). In tables 27, 28 and 31 it matches nothing, and the entry says so (note
  `no_stop_codon`; its contexts, or a long one's paths, `exact` 0); never read such a 0 as an absence of the
  protein. Such a pattern of at most k bases is answered without a search, so even after an earlier pattern's
  budget stop it states `exact` 0 with `stop: null` (SPEC §7.6, the one exception to "every later pattern
  `unknown`"). The slot code `stop_unsupported` (`4596bb3b`'s answer to a `*`, in that commit's fixtures) is
  retired: no current build answers it; a host still on `4596bb3b` may, so keep passing it through as a slot
  error.
- **Identity**: the mask and the Bloom filter are derived data of the graph, outside `index_fp`; a host that gains
  a mask keeps its `index_fp` (only its counting becomes `exact`).

## 2. When

| backend milestone | content | state |
|---|---|---|
| 1 | count (`mode: count`) and the label-free extraction (`all_or_count` / `partial` with `labels: none`): k-mers, offsets, strands, node and row ids; exact DNA and IUPAC; both strands; `suffix` and `any_offset`; single-graph servers; the capabilities block; `metagraph pattern` CLI | running now; contract freezes on its commit |
| 3 | `labels: all`: label discovery and placement (record, 1-based position, strand) on BASIC indexes with record mapping | in the build (SPEC §14), with fixtures |
| 4 | patterns longer than k (extension), per-label `support`, `require_support`; opt-in: only a request with `long_search: "paths"` gets paths (new fields `sequence`, `anchor_kmer`; `kmer` keeps its meaning), every other request keeps today's anchor-only answer (SPEC §12.1) | in the build (2026-10-08, SPEC §17), with fixtures (`paths*`, `support_unavailable`) |
| 5 / 5b | peptides (codon automaton); annotation predicates (`any`, `all`, `none`, `at_least`, `and`/`or`/`not`) | 5 in the build (2026-10-08, SPEC §12.2), with fixtures (`peptide*`, `genetic_code_unknown`); 5b later |
| 6 | multi-graph servers (per-shard budgets, barriers, shard identity per result), the real-index benchmark | after 5 |
| 7 | this service's job type (the backend's Python client methods are deferred until needed) | with you; on refseq33m-experimental after backend milestone 1, on chunked databases after milestone 6 (§3.1 item 3) |
| mask | refseq33m-experimental's graph has no `.edgemask` file. Since the owner's decision #16 (2026-10-08) the route answers without it: `mask: absent`, `counting: "upper_bound"`, counts `bounds` with an `estimate` where they cannot be proven, lists exact (before, it answered `mask_required`). For exact counts the owner runs `metagraph transform --mask-dummy` once on mex (decision #18: in a staging-only directory, on the host rather than in the 128 GiB container: it holds a transient bit vector of edges + 1 bits, about 78 GB, beside the graph). Node ids, rows, the annotation and `index_fp` stay (decision #17: the mask is derived data); `/stats` `graph.nodes` becomes the k-mer count and a `.bloom` beside the graph starts loading; the block then says `counting: "exact"`. `--pattern-build-mask` (the mask built in memory at every start-up) is for small indexes, not for refseq33m | the route answers from the `update.sh` that deploys it; exact counts after the mask (#18) |
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
2. **(required)** The merged view of a job: counts of one unit are added with the relation algebra of SPEC §7.4
   (`Count::operator+=`), not "the weakest relation wins": `exact` + `exact` = `exact`; `unknown` + `unknown` =
   `unknown`; `unknown` or `at_least` with anything = `at_least` of the known lower bounds (`exact` 7 +
   `unknown` = `at_least` 7); otherwise `bounds` with `lower` and `upper` summed (`bounds` [3, 5] + `exact` 7
   = `bounds` [10, 12]). Each answered entry is mapped back to its pattern's position in the job's request (the
   pattern chunk's offset plus the entry's position, SPEC §7.9), never matched by text or `id` (two patterns can
   share both); a pattern's counts are then added across its graph chunks, where a k-mer of two graphs is two
   contexts (DESIGN §8). Label counts are **never** summed across tasks (a label in two chunks would count
   twice): they stay per task. `retrieval_complete` (and any "complete" of the job) only when
   every task that contributes to that pattern answered with `retrieval_complete: true`; a failed, refused or
   missing task makes it incomplete. `withheld` and `stop` carried per task with the graph and its reason; rows
   are contexts, never labels.
3. **(required)** A chunked database is many graphs on one multi-graph server process, selected per task through
   `graphs: ["{label}-{i}/{N}"]` as `/search` does (`app/download_depth.py`, `enumerate_leaf_specs`). The job type
   is built and tested on refseq33m-experimental (one graph) with backend milestone 1; serving the chunked
   databases waits for backend milestone 6, which keeps `graphs` exactly as `/search` selects a shard.
4. **(required)** Admission by the queue, **shared with search and traversal** (the owner: "traverse, pattern
   match and normal search need to share the pool inside async and inside the sync server"): the per-database
   queues, `META_DB_CAPS` and the distributed semaphore count search, pattern and traversal calls against one cap
   per database server (the cap search uses today, e.g. 15 on sra-logan-chunks). The server has one request pool
   (`-p`) for every route, the GET capabilities routes included, and reserves nothing per route, so the cap per
   host must stay below the server's `-p` (`-p` ≥ cap + 1; a probe that waits on a full pool is not a dead host).
   The engine needs nothing beyond "at most cap concurrent calls per database server across all routes"; a
   request's memory is bounded on the label-free path by the caps (`partial` and `all_or_count` keep about
   `max_contexts` descriptors, `count` the search frontier; on an even-k primary index up to `max_steps` × 88
   bytes; SPEC §7.6) and with `labels: all` by `max_memory_mb`. `max_memory_mb` is a deterministic account
   (a model in bytes, so that where a request stops does not depend on the allocator; SPEC §14.4), not the
   process's resident memory, and reads with `allow_unbudgeted_annotation` can pass it (the label names of one
   read, and the read's own working memory, are not bounded by it): size a host's memory with headroom above
   cap × `max_memory_mb`. The traversal jobs' separate per-host tokens move into the same semaphore (your
   follow-up).
5. **(required)** Caps at submit, never truncation afterwards: `max_contexts`, `max_patterns`, `time_budget_ms`
   above the service's ceilings are refused with a 400 naming the field (the strategy validator's rule: refuses,
   never lowers). An answer is passed through whole; cutting it on the way out would falsify
   `retrieval_complete`. The task timeout derives from the request's `time_budget_ms` plus an allowance. Three
   clocks stay distinct: `time_budget_ms` is the backend's deadline for preparing the answer from the parsed
   body (the queue wait and the transport are outside it; the work stops well before it, SPEC §7.6), the
   HTTP client timeout lies above it, and the job's lifetime is the service's own.
6. **(required)** Served only for databases whose host carries the `pattern` block (§1) with a contract version
   the service implements (1) **and** `available: true` (SPEC §10.3). Any other version (missing, malformed,
   lower or higher), or no block, answers `pattern_unsupported` with the host's feature list, as
   `traversal_disabled` does. `available: false` answers unavailable with the host's `unavailable_reason` kept
   verbatim (`mask_invalid`, `alphabet_untested`, … an older build's `mask_required`, or one the service does
   not know);
   `available: null` means the index is loading: availability unknown, re-probe later, never cache it as
   unavailable. Every option is gated on its capability: label projections exactly when the host's
   `projections` list has them; on `annotation: "unbudgeted"`, `labels: all` only with the caller's explicit
   consent, sent as `allow_unbudgeted_annotation: true`; kinds and scopes by their lists; values within `caps`.
   Every request sends `output.labels` explicitly (SPEC §1; the default is `none` in contract version 1). The
   probe that reads capabilities already exists (`app/traversal/probe.py`).
7. **(required)** The answer passed through with the service's additions only: `database`, the untrusted-data
   notice on label and record strings, the standard error envelope. **Never** sum per-shard label counts into one
   number, never drop `relation`, `withheld`, `retrieval_complete` or `absence_scope`, never add labels of its own.
   `hit_unit` does not apply to a context count; the tool says what a context is instead: one distinct k-mer of
   the index that contains the pattern, at one offset and orientation — not a record, not an occurrence, not a
   hit; a k-mer present in ten samples is one context.
8. **(required)** The submit and results tools' docstrings explain the three things an agent acts on:
   `withheld: count_above_threshold` (narrow the pattern, scope or strand), `withheld: discovery_budget` (shorten
   or move an N run inside the pattern, or restrict the strands, SPEC §7.5; a filter does not help),
   `withheld: annotation_budget` (ask for the count or `labels: none`, or a narrower pattern); that `count` first,
   then `labels: none`, is the cheap way to look at a new pattern; they use the service's existing strand
   vocabulary, and say that canonical and primary indexes report contexts only — `orientation` instead of a
   strand, no position.
9. A stored answer keeps its `determinism`; a `time_limited` one is marked as not reproducible where the
   service shows it. `full` holds within one backend build and configuration: after a backend update an equal
   request may answer with other `work`, stopping points and `at_least` / `bounds` values, with every field's
   meaning kept (SPEC §1, semantic compatibility).
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
`any_offset` DFS per task, and the time kept back for writing the answer grows with what the answer holds (SPEC
§7.6), so a stopped task still answers with counts. Nothing in the engine assumes a short run: the clock is read
every 4,096 steps (and every 64 released contexts) whatever the budget; the retained descriptors are bounded by
their caps, not by time, and `max_memory_mb` bounds the labelled retrieval (as an accounted model, not the
process's memory, and not the unbudgeted reads; §3.1 item 4); `max_steps` is the host's
(`--pattern-max-steps`, default 10⁸ per request, per shard from milestone 6): a request cannot raise it (a
larger value is lowered and listed in `limits.clamped`), and it does not grow with `time_budget_ms`, so with a
600 s budget it, not the deadline, usually stops a long discovery (`stop: max_steps`, `withheld:
discovery_budget`). The cost of a long call is the request-pool slot it holds, counted by the shared
per-database cap of item 4; the server stops a `/pattern` request whose client has closed (or half-closed) its
connection at its next clock reading and writes neither an answer nor an error, so a task the service cancels frees the backend's slot
once its HTTP connection is closed. `time_budget_ms` is not the task's latency: the wait for a backend thread
comes before it, the work stops `finalize_reserve_ms` plus the estimated writing time before it (seconds, for
an answer with many results), and a finalisation that overruns it is still a 503 `deadline` (SPEC §7.6); the
task's HTTP timeout and the job's lifetime are separate, larger clocks.

## 4. Constraints

- Absence claims follow SPEC §9; show `absence_scope` and `retrieval_complete` with every answer. A count-mode or
  label-free answer establishes no absence of a label, a sample or a record, whatever its counts; an `exact` 0
  count establishes absence of its unit in its scope (`exact` 0 contexts: no k-mer of the index holds the
  pattern there); an incomplete list is not an inexact count (a `max_contexts`
  cut after a completed discovery, `labels_cut` and `occurrences_cut` keep `exact` counts), but an item missing
  from an incomplete list is not absent: one item's absence needs `retrieval_complete: true` or an `exact` 0.
- A `determinism: time_limited` answer is not reproducible; say so where the service caches answers.
- The service does not retry a `withheld` answer, or a 503 `deadline`, with a larger cap on its own; the agent
  decides.
- Nothing here touches the search path or the traversal jobs.

## 5. Tests and acceptance

- Unit: capabilities gating (block present, absent; contract version 1, missing, malformed, lower, higher;
  `available` true, false with its reason, null; `counting` exact and upper_bound, an estimate shown as one); the
  per-task call and its timeout; the merge (the algebra of
  §3.1 item 2 with its two examples, entries mapped by position across pattern chunks, label counts never
  summed, `retrieval_complete` only with every contributing task, per-task `withheld`); the pass-through keeps
  every field; the error envelope (SPEC §6: an unexpected status or body is a backend failure, never an empty
  result); the docstring's tool table regenerated. The fixture bodies of §2 feed these.
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
