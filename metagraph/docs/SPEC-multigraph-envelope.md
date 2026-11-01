# Multi-graph answers: the envelope, the search cuts, deadlines, the request registry and the selectors

Release 10, Phase 1 of the multi-graph plan (ratschlab/metagraph-suite `docs/PLAN-multigraph-implementation.md`,
section 1 is the contract; items E1 to E4, E9, E10, E12, E18, E20). This document states the engine's shapes and
rules. The samples in `metagraph/tests/data/multigraph_contract/` are their reference instances: every sample is
`{route, status, note, request?, body}`, the body as a client receives it. Where this text and a sample differ,
this text wins and the sample is a bug.

A server is in multi-graph mode when it is started with a graph list (csv); a name of the list may hold several
(graph, annotation) PAIRS. A single-graph server answers every route as before this document, except the additive
`/search` cut (section 3) and the new request fields (sections 4 and 5).

## 1. The envelope (resolve, traverse, pattern, and the capabilities routes)

On a multi-graph server `/resolve`, `/traverse` and `/pattern` answer HTTP 200 with the envelope, also for a named
graph (one entry). Its members, listed here in reading order (a JSON object: the order of its members, of an entry's
and of `deadline`'s carries no meaning; the server writes every object's members in byte order of their names, as
all its answers, and the samples list them in this reading order):

| Member | Type | Meaning |
| --- | --- | --- |
| `envelope_version` | 1 | this document's version of the envelope |
| `route` | string | `/resolve`, `/traverse`, `/pattern` |
| `pattern_contract_version` | string | `/pattern` only, as before |
| `graphs` | [string] | the selected names, in byte order |
| `answered`, `refused` | int | entries by outcome; `answered + refused` = number of entries = selected pairs |
| `not_reached` | [pair] | the pairs never started (section 4), `{graph, graph_path, annotation_path}` |
| `unfinished` | [pair] | the pairs started and cut at the work cutoff (section 4) |
| `request_id` | string or null | the client's `request_id` (resolve, pattern) or `attempt_id` (traverse); null if none |
| `server_instance` | string | this server process (the attempt ledger's value) |
| `deadline` | object | `{not_after_ms, finish_by_ms, effective_deadline_at, cutoff_at}`; nulls where not given or not applicable |
| `usage` | object | `/traverse` only, for a ledger-managed request (one with `attempt_id`): the attempt's usage, the object a single-graph traverse states at its top level |
| `answers` | [entry] | one entry per selected pair, in the index's pair order (names in byte order, each name's pairs in list order, a pair listed twice under a name once) |
| `timing`, `work` | object | the request's (`work.memory_bytes`: the answered pairs' accounts summed); `timing.elapsed_ms` is always present, an entry's own `timing` follows its request's `output.timing` |

An entry is either the pair's single-graph answer with the tags merged in, or a refusal:
- answered: `{graph, graph_path, annotation_path, index_fp, outcome: "answered", ...the single-graph answer}`.
  On `/traverse` the single-graph answer's `usage` moves to the envelope (it is the attempt's, and the attempt is
  the request's). The entry's `index_fp` equals the identity inside the answer (`capabilities.index_fp`).
- refused: `{graph, graph_path, annotation_path, index_fp, outcome: "refused", refusal: {http_status, error,
  code?, reason?, usage?}}`: what the pair alone would have been answered, or a `not_reached` / `unfinished`
  outcome (section 4).

The capabilities routes (`GET /traverse/capabilities`, `GET /pattern/capabilities`) answer without a pair
selection with the envelope over every pair (up to `max_graphs_without_selection`): `envelope_version`, `route`,
`graphs`, `answered`, `refused`, `server_instance`, `answers` (each entry the pair's capabilities block, tagged).
With a pair selection (section 6) they answer a one-entry envelope. Beyond the cap a probe without `?graph=` keeps
its 400 text of before release 10, byte for byte (section 6).

### 1.1 `/traverse` (E10)

The envelope wraps HTTP 200 answers only. Every other answer of `/traverse` stays HTTP-level and unwrapped,
exactly as before, because the service's attempt ledger reads those statuses:
- 409: a repeated `attempt_id` that is running, retained, held or tombstoned (`{error, attempt}`; also an
  `attempt_id` equal to a running or retained `request_id`, section 5, `attempt` then being that request's entry with
  `attempt_id` and `state` `running` or `finished`), and an expired `not_after_ms` (`{error, state: "expired", ...}`);
- 503 with `usage` when the attempt reaches its bound;
- 400 with `usage` for a managed request (without `usage` for an unmanaged one).
429 exists only on `POST /traverse/cancel` (an unknown id that could not be tombstoned); that route and
`GET /traverse/attempt/{attempt_id}` are unchanged and never wrapped.

A `/traverse` request runs on exactly one pair (until the graph-less traverse of a later release): a name with
several pairs, or several pairs, is refused with the candidates (section 6).

## 2. Labels absent from a pair (E9)

On a multi-graph server a request that names labels a pair does not hold is answered by that pair, not refused:
- `/resolve` `labels`: the absent ones are dropped and listed in the entry's `unknown_labels`;
- `/traverse` seed labels: listed in that seed result's `unknown_labels` and in its `dropped_labels`
  (`{label, reason: "not_in_graph"}`); a seed whose labels are all absent gets no walk, outcome
  `walks: "labels_not_in_graph"`;
- `/traverse` `strategy.labels.extra`: listed in the entry's `unknown_labels`.
A single-graph server keeps the 400 "Label not found in the annotation".

## 3. `/search` (E1)

The answer keeps its upstream shape (a list of `{seq_description, results, ...}`). Added:
- multi-graph server: every row carries `graph`, `graph_path`, `annotation_path`; the rows of one query sequence
  come in the index's pair order (not in the order the pairs finish);
- every server: each query sequence's object has `cuts`, one entry per pair searched (empty pairs included) in pair
  order: `{graph, graph_path, annotation_path, index_fp, returned, complete, nth_score, tied, tied_returned}`; on a
  single-graph server one entry with null pair tags.
  - `returned`: rows of this pair; `complete`: every match above the request's thresholds was returned (no cut):
    then `nth_score` null and `tied`, `tied_returned` 0; `complete: true` with nothing returned means nothing
    matched, and N = 0 with matches states `returned 0, complete false, nth_score null, tied 0`;
  - otherwise `nth_score` is the score of the N-th returned row, `tied` the number of matching labels with that
    score (returned or not), `tied_returned` how many of them were returned;
  - on COUNTS annotations also `nth_median` and `tied_same_median` (the labels with the N-th score and median).
    The median is read only where the N-th score is shared: a cut with `tied` 1 states `nth_median` null (and
    `tied_same_median` 1).
- The order before every engine-local cut (C5): score descending; then, on counts annotations only, the median
  abundance descending; then label ascending (byte order). The median is the LOWER median of the label's
  abundances over the query's k-mer positions where the label is present (abundance at least 1, zeros excluded),
  each position once: sorted ascending, the element at index floor((n - 1) / 2); it is read only for the rows
  it orders (those sharing a score). A counts annotation is an integer matrix that is not a coordinate matrix
  (row_diff_int_brwt, int_row_diff_disk, the integer CSC forms). On coordinate and binary annotations ties are
  ordered by label.
- The result cache key includes the selected names (with their pairs), so a pair listed under two names never
  returns rows tagged with the other name.
- `/search` takes `graphs` only (no `pairs`; it has no re-ask): `pairs` is refused 400 `pairs_unsupported` on every
  server, single-graph included, before any cached answer is looked up. On a single-graph server the field used to
  be ignored (a request it refused for `graphs` keeps that refusal).

## 4. Deadlines and the work cutoff (E4, E20)

- `not_after_ms` (integer, ms since the epoch) is the latest START time, as in SPEC-labeled-traversal-core
  section 5: `/traverse` as before; `/resolve` and `/pattern` read it too. A request whose `not_after_ms` has
  passed when its handler starts is refused 409 `{error, state: "expired", not_after_ms, server_time_ms}` and nothing
  runs. In a fan-out, a pair not started by then is `not_reached` with reason `not_after_ms`.
- `finish_by_ms` (integer, ms since the epoch) is the completion deadline, on `/traverse`, `/resolve`, `/pattern`.
  Each route's deadline becomes the minimum of its own and `finish_by_ms - clock_skew_allowance_ms`. A `/resolve`
  without `bounds.time_budget_ms` gets the server maximum as its deadline. The in_ram memory wait is not entered
  once `finish_by_ms` has passed and is given up when it passes; a load already started is not interrupted.
- The EFFECTIVE envelope deadline is min(request arrival + the route's request budget, `finish_by_ms` -
  clock_skew_allowance_ms), the request budget being `/pattern`'s `time_budget_ms`, `/resolve`'s
  `bounds.time_budget_ms` or the server maximum, the attempt's bound for `/traverse`. The WORK CUTOFF is the
  effective deadline less the route's finalisation reserve. On `/traverse` the bound already holds `finish_by_ms`
  and the transport cap, and the reserve is the delivery reserve: `effective_deadline_at` is the request's receipt
  plus `usage.bound_ms`, `cutoff_at` its receipt plus `usage.bound.walk_until_ms`. At the cutoff the fan-out stops
  admitting and waiting:
  - a pair not yet started (queued, waiting for room in the pool's bounded queue, waiting for memory) is
    `not_reached`: refusal `{http_status: 503, code: "not_reached", reason: "cutoff" | "not_after_ms" | "transport",
    error}`; a full queue never holds the handler past the cutoff or past its client's departure;
  - a pair started and not finished is `unfinished`: refusal `{http_status: 503, code: "unfinished", reason:
    "cutoff" | "loading", error, usage: {state, started_at, elapsed_ms, cpu_ms, memory_bytes}}`; it is told to stop;
    a load in progress may run past every deadline;
  - the completed pairs are written beside them. Each pair keeps its own per-pair budget, which starts when it
    starts and, when the pair runs pooled (two or more pairs), is capped at the work cutoff (the owner's decision
    of 2026-11-01): a pair whose budget would end later stops by time before the cutoff and answers what it found,
    so a time-limited answer always survives instead of racing the cutoff. The cap is stated in the pair's answer
    only when it binds, i.e. when the work stopped by time under it: `limits.time_budget_ms` is then the budget it
    ran under (the time left before the cutoff at its start) and `limits.clamped` gains its entry (`/pattern`:
    `{field: "time_budget_ms", requested: <its time_budget_ms or the default>, effective}`; `/resolve`: `{field:
    "bounds.time_budget_ms", requested: <its budget or null>, effective}`), after the `finish_by_ms` entry, which is
    stated whenever it lowers the budget, as on a single-graph server. A pair that finishes before the cutoff
    answers exactly its single-graph answer. A `/resolve` pair without a deadline (no `bounds.time_budget_ms`, no
    `finish_by_ms`) has no budget to cap and runs as its single-graph answer does (a deadline would change its
    work: the explicit labels' hits fetched in pieces); it is told to stop at the cutoff. A pair that cannot stop
    by then (a load in progress, a step that reads no clock, a pair without a deadline) is `unfinished`.
- `deadline.effective_deadline_at` and `deadline.cutoff_at` are stated only when the request gives `not_after_ms`
  or `finish_by_ms` (samples 01 and 02 give neither); otherwise they are null, the cutoff applying all the same.
- The 899 s transport cap bounds every answer.

## 5. The request registry (E4)

- `/resolve` and `/pattern` take an optional `request_id`, the same token rule as `attempt_id`
  (`[A-Za-z0-9._-]+`, length-capped), checked by the strict parsers before anything is registered. `/traverse`'s
  `attempt_id` is its request handle; no second id.
- On `/resolve` and `/pattern` the request's own 400s (the body, the selection, the fields; on `/resolve` also what
  the request and k decide of its explicit seeds: a `kmer_interval` empty or out of range, a seed label that is not
  one of `labels`) come before the registry's 409 and the expired 409: a malformed request with an expired
  `not_after_ms` gets 400. A 400 never holds a `request_id`: one that only the work can find (a seed not fully
  present in the graph, a seed label a discovery did not find, a `/pattern` the index cannot answer) frees the id when
  it is answered, so `GET /requests/{id}` answers 404 and a retry runs. `/traverse` keeps checking expiry before it
  parses the body.
- A repeated `request_id` that is running, or retained after it finished, starts no work: 409, body = the registry
  entry (below) plus `error`.
- `request_id` and `attempt_id` are one namespace, as `GET /requests/{id}` reads both: a `request_id` equal to a
  running, retained or held `attempt_id` is such a 409 (the entry being the attempt's, as `GET /requests/{id}`
  answers it), and an `attempt_id` equal to a running or retained `request_id` is `/traverse`'s 409 of a repeated
  `attempt_id` (section 1.1). The check of both registries is atomic, so two requests with one id are never both
  admitted, in either order; once its retention ended the id runs on either route. A tombstone holds an
  `attempt_id` only.
- `GET /requests/{request_id}` answers the entry: `{request_id, server_instance, route, state: "running" |
  "ended", received_at, finished_at, reason, response: {status, bytes, written} | null, retained_until, pairs:
  [{graph, graph_path, annotation_path, state: "queued" | "waiting_memory" | "loading" | "running" | "ended" |
  "never_started", outcome, usage}]}`; for `/traverse` also `attempt` (the object of
  `GET /traverse/attempt/{attempt_id}`). An entry is retained for a period and a count set by server flags, then
  the id is unknown: 404 `{error, request_id, server_instance}`. An id that breaks the token rule is 400. `reason`
  on `/resolve` and `/pattern`: `answered` (200 written), `client_gone` (nothing written), `deadline` (503),
  `error` (any other status but 400, which leaves no entry).
- A client gone, or a connection lost, stops the work as before; the entry records the end (`reason:
  "client_gone"`, with usage). That holds for an error too: on `/resolve` and `/pattern` nothing is written to a
  client that left once the request has its entry, and the entry records `client_gone`; a delivered error stays
  `error`. An answer written but lost on the way is recorded `ended` with its status and bytes; results are not
  retained.

## 6. Selectors (E12)

| Route | `graphs` | `pairs` | Note |
| --- | --- | --- | --- |
| `/search` | yes | refused (`pairs_unsupported`) | also on a single-graph server |
| `/resolve`, `/pattern` | yes | yes | any number of pairs |
| `/traverse` | `graph`, or `graphs` resolving to ONE pair | exactly one pair | several pairs: `traverse_single_pair` |
| `GET /traverse/capabilities`, `GET /pattern/capabilities` | `?graph=` | `?graph=&graph_path=&annotation_path=` | one pair per probe |

`pairs: [{graph, graph_path?, annotation_path?}]`: each object names exactly one pair of that name; `graph_path` may
be omitted when the name has one graph, `annotation_path` is required when the name's pairs share a graph. `graphs`
and `pairs` together are refused. Every selector refusal is a 400 before any pair starts, `{error, code,
candidates?}`, with `code` one of `selection_conflict`, `pair_not_found`, `pair_ambiguous`, `traverse_single_pair`,
`graph_unknown`, `pairs_unsupported`; `candidates` lists the pairs to choose from where a choice exists. The codes
are for selection outcomes only: a malformed selector (`pairs` not a non-empty array of objects, an unknown member,
a missing `graph`, a path that is no string, `graph_path` without `graph`) is a plain 400 without a section 6 code
(on `/pattern` the route's own `invalid_request`, as every `/pattern` 400). Without a selection a request runs on
every pair, up to `max_graphs_without_selection` names (that refusal stays).

The multi-graph 400 texts are those of before release 10 (8fa32e553), byte for byte: a request that server refused
is refused in the same words, a section 6 code being only an added member (a client matches these texts whole;
any later change is announced first). The names cap's three texts, `{N}` the cap and `{M}` the names listed:
- `/search`: `Bad request: requests without names (no "graphs" field) are only supported for small indexes
  (<={N} names)`;
- `/pattern`, and `/resolve` (which reads the cap since release 10) in the same words: `request.graphs: required on
  this server, which hosts {M} graph names (more than {N}; GET /capabilities lists them)`;
- a capabilities probe without `?graph=` (also one with `graph_path` or `annotation_path` but no `graph`):
  `Bad request: in multi-graph mode GET <route> needs ?graph=<name> (and graph_path=<path> when the name spans
  several graphs); GET /capabilities lists the graphs`.
The other texts kept: an unknown name on `/resolve`, `/traverse` and a probe `Bad request: unknown graph '<name>'`
(on `/search` `Request with an uninitialized graph <name>`, on `/pattern` `request.graphs: unknown graph '<name>'
(GET /capabilities lists the graphs)`), each now with `graph_unknown`; on `/resolve`, `/traverse` and the probes
without `pairs` or `annotation_path` the release-9 selection texts (`... spans several graphs (...); pass
'graph_path' to pick one`, `... with several annotations ...`, `graph_path '<path>' is not part of index
'<name>'`, `'graphs' names the one graph a traversal reads: expected [name]`, `'graph' (index name) is required
in multi-graph mode`, `give the graph as 'graph' or as 'graphs', not both`), coded `traverse_single_pair`,
`pair_ambiguous`, `pair_not_found` or `selection_conflict` where section 6 has a code. On a single-graph server
`pairs` is an unknown field, as before (`request: unknown field 'pairs'`).

## 7. Resolve discovery ties (E18)

Discovery ranks labels by k-mers supported (descending), ties by label name (byte order), then truncates at
`discover_max_labels`. `labels_truncated` gains `tied_at_boundary` (the labels with the last kept label's
support, kept or dropped) and `tied_kept` (how many of them were kept). The k-mers supported are the number the
answer states as `kmers_supported`: under `support: "trace"` the traced k-mers, otherwise the k-mers present;
the ranking, the cut and every count of `labels_truncated` (`min_kept_kmers`, `max_dropped_kmers`,
`dropped_full_length`, `tied_at_boundary`, `tied_kept`) are on that scale. Every `/traverse` capabilities block
(`GET /traverse/capabilities`, each entry of its envelope, the `capabilities` of a `/resolve` or `/traverse`
answer, single-graph servers included) states `discovery_trace_scale: "traced"`; absent, a discovery under trace
ranked by k-mers present (before release 10's post-merge fix).

## 8. Capabilities (E3)

`GET /capabilities` on a multi-graph server adds `composite_fp` (sha256 over the pairs' `index_fp` sorted as bytes,
joined by `\n`, no trailing newline; null when any pair has none; each loaded (graph_path, annotation_path) pair
counted once, so a pair listed under two names counts once), `composite_rule` (that sentence), and a
`multigraph` block advertising `envelope_version`, the envelope routes, the selectors per route,
`traverse_pairs_per_request`, the routes reading `not_after_ms`, `finish_by_ms` and `request_id`, the registry route,
`unknown_labels`, `search_cuts`, `resolve_boundary_ties` and `discovery_trace_scale` (`"traced"`, section 7).
Clients switch on these, never on dates.
