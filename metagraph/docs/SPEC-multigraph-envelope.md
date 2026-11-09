# Multi-graph answers: the envelope, the search cuts, deadlines, the request registry and the selectors

Release 10, Phase 1 of the multi-graph plan (ratschlab/metagraph-suite `docs/PLAN-multigraph-implementation.md`,
section 1 is the contract; items E1 to E4, E9, E10, E12, E18, E20). This document states the engine's shapes and
rules. The samples in `metagraph/tests/data/multigraph_contract/` are their reference instances: every sample is
`{route, status, note, request?, body}`, the body as a client receives it. Where this text and a sample differ,
this text wins and the sample is a bug.

A server is in multi-graph mode when it is started with a graph list (csv); a name of the list may hold several
(graph, annotation) PAIRS. A single-graph server answers every route as before this document, except the additive
`/search` cut (section 3), the new request fields (sections 4 and 5) and, with a cap set, the busy refusal (section
1.1, E19).

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
- 400 with `usage` for a managed request (without `usage` for an unmanaged one);
- 503 `busy` (E19, below), without `usage`.
429 exists only on `POST /traverse/cancel` (an unknown id that could not be tombstoned); that route and
`GET /traverse/attempt/{attempt_id}` are unchanged and never wrapped.

**Busy (E19).** `/traverse`, `/resolve` and `/pattern` each take an optional cap on the requests the server runs at
once (`--traverse-max-concurrent`, `--resolve-max-concurrent`, `--pattern-max-concurrent`; none by default;
`SPEC-labeled-traversal-core.md` §10.3, "Handler capacity", states what a cap counts and the threads it keeps for
the light routes). A request of a route at its cap is answered HTTP 503 with `Retry-After: 1` and the body
`{"error": <message>, "code": "busy", "server_instance": <this process>}`: no `usage`, no envelope, on a
single-graph server too. Its `server_instance` lets a ledger that sent `expect_server_instance` see a restarted
process at once, as its 409 `instance_mismatch` would tell it, rather than once a slot is free. It is refused when
its handler starts, after the 503 of a single index still loading and before its body is parsed: before every 400
and 409 and before its `attempt_id` or `request_id` is registered (section 5), so nothing of it ran, nothing of it
is held, and a re-send (the same id included) after `Retry-After` is admitted as a first copy would be. It says
nothing about other copies of its id: like the loading 503 it is no ground of the release rule; a ledger re-sends,
or releases on that rule's grounds. This is the shape a client already retries for a loading index (503,
`Retry-After`, no `usage`), told apart by its `code`. 429 stays reserved for the cancel (a 429 releases nothing).
During a shutdown a request the cap would refuse is not answered (its connection is closed, as for a client gone).
`/search` has no cap.

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
- `kmers_in_target` (since MetaGraph master #634, 2026-05-17; single-graph and multi-graph servers alike): on the
  rows of a coordinate query (`query_coords`: the rows that carry `kmer_coords`) of a pair that has a `.seqs`
  (CoordToHeader) record mapping, a non-negative integer, the target record's k-mer count (its length − k + 1,
  `CoordToHeader::num_kmers_in_sequence`), so that `kmer_count / kmers_in_target` is the fraction of the target
  the query's k-mers cover. It is omitted where the label is a whole column (a pair without a `.seqs`, or a
  server started with `--no-coord-mapping`: its `kmer_coords` are positions in the column), and never stated on
  signature or counts rows. Sample 25: one pair with a `.seqs` and one without, on one multi-graph server.
- The order before every engine-local cut (C5): score descending; then, on counts annotations only, the median
  abundance descending; then label ascending (byte order). The median is the LOWER median of the label's
  abundances over the query's k-mer positions where the label is present (abundance at least 1, zeros excluded),
  each position once: sorted ascending, the element at index floor((n - 1) / 2); it is read only for the rows
  it orders (those sharing a score). A counts annotation is an integer matrix that is not a coordinate matrix
  (row_diff_int_brwt, int_row_diff_disk, the integer CSC forms). On coordinate and binary annotations ties are
  ordered by label.
- The result cache key includes the selected names (with their pairs), so a pair listed under two names never
  returns rows tagged with the other name.
- `GET /column_labels` (H3), which tags nothing with a name, lists each distinct pair once, at its first place in
  the answers' order (the names in byte order, each name's pairs in the list's order), its labels in column
  order (section 8).
- The CLI's `metagraph query --query-mode labels` (every label above the thresholds, no counts, no cut) keeps the
  annotation's column order, not the C5 order; `/search` has no such mode, so it is not reachable over HTTP.
- `/search` takes `graphs` only (no `pairs`; it has no re-ask): `pairs` is refused 400 `pairs_unsupported` on every
  server, single-graph included, before any cached answer is looked up. On a single-graph server the field used to
  be ignored (a request it refused for `graphs` keeps that refusal).

## 4. Deadlines and the work cutoff (E4, E20)

- `not_after_ms` (integer, ms since the epoch) is the latest START time, as in SPEC-labeled-traversal-core
  section 5: `/traverse` as before; `/resolve` and `/pattern` read it too. A request whose `not_after_ms` has
  passed at its admission (section 5) is refused 409 `{error, state: "expired", not_after_ms, server_time_ms}` and
  nothing runs. In a fan-out, a pair not started by then is `not_reached` with reason `not_after_ms`.
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
- A one-pair request (the selection holds one pair) runs inline, as its single-graph route does. The pair's own
  deadline, which starts when its work starts (after an in_ram load), is capped by `finish_by_ms` -
  clock_skew_allowance_ms, never by the work cutoff; a load in progress is not interrupted; once the pair has
  started it is never `not_reached` or `unfinished` by the cutoff, it ends by its own deadline; and its envelope,
  answered or refused, is written as soon as the pair ends (a refused pair's envelope is bound by the transport cap
  only, so it can be written after `finish_by_ms`). Before it starts, an in_ram load is not started once its wait
  limit has passed (`finish_by_ms` - clock_skew_allowance_ms, else the transport cap) or its client has left:
  nothing is loaded and its entry is refused 503 `not_reached` (reason `cutoff`, `transport` or `not_after_ms` as
  above). `not_after_ms` refuses it as in a fan-out. The envelope's `deadline` states the request's instants as for
  any request.
- `deadline.effective_deadline_at` and `deadline.cutoff_at` are stated only when the request gives `not_after_ms`
  or `finish_by_ms` (samples 01 and 02 give neither); otherwise they are null, the cutoff applying all the same.
- The 899 s transport cap bounds every answer.

## 5. The request registry (E4)

- `/resolve` and `/pattern` take an optional `request_id`, the same token rule as `attempt_id`
  (`^[A-Za-z0-9._:-]{1,128}$`), checked by the strict parsers before anything is registered. `/traverse`'s
  `attempt_id` is its request handle; no second id.
- `/resolve` and `/pattern` take an optional `expect_server_instance`, as `/traverse` does
  (SPEC-labeled-traversal-core §5): the `server_instance` the request is meant for, a string matching the ids'
  pattern `^[A-Za-z0-9._:-]{1,128}$`, compared as given, **with `request_id` only**. Anything else (null, a
  non-string, `""`, more than 128 characters, another character), or the field without `request_id`, is a 400
  naming the field, one of the request's own 400s (before every 409). A server whose own `server_instance` is
  another refuses the request at its admission, before anything is registered: **409** `{error, state: "instance_mismatch",
  expect_server_instance, server_instance, request_id, not_after_ms?}`. Nothing runs and nothing is registered:
  `GET /requests/{id}` on that server answers 404 (never admitted). It is the first of the admission refusals:
  `instance_mismatch`, then the repeated id, then `expired` (section 4), the order `/traverse` has. On a
  multi-graph server it is request-level: no pair starts and no envelope is written. The CLI, whose instance is
  its own and random, refuses any value it is given (the 409 body on stdout, exit status 1), as for `/traverse`.
- A request's ADMISSION is the point where the server decides to run it: the check of `expect_server_instance`
  (`instance_mismatch`, above), of a repeated id (the 409 below) and of an expired `not_after_ms` (section 4),
  made together in that order, after which a request with an id is
  registered under it. On `/resolve` and `/pattern` it comes after the request's own 400s (below), on `/traverse`
  before its fields are parsed (its body read as JSON, its ids and `finish_by_ms` checked). A request without an
  id is admitted at the same point, without an entry.
- On `/resolve` and `/pattern` the request's own 400s (the body, the selection, the fields; on `/resolve` also what
  the request and k decide of its explicit seeds: a `kmer_interval` empty or out of range, a seed label that is not
  one of `labels`) come before every admission 409 (`instance_mismatch`, the repeated id, `expired`): a malformed request with an expired
  `not_after_ms` gets 400. The seeds are checked in their order, the work's, which checks a seed's labels before
  the next seed: under a discovery the seeds after the first that names labels are checked by the work. An
  HTTP-level 400 never holds a `request_id`: one that only the work can find (a seed not fully present in the graph,
  a seed label a discovery did not find, a fault of a seed after it, a `/pattern` the index cannot answer) frees the
  id when it is answered, so `GET /requests/{id}` answers 404 and a retry runs. On a multi-graph server such a 400
  of one pair is a refused entry (`refusal.http_status` 400) of a 200 envelope: the request was answered, its id
  is retained like any answered request's (`GET /requests/{id}` states it ended, that pair `refused`), and a
  corrected request takes a new id. A fan-out is never forgotten because one of its pairs gave 400; only a 400 of
  the whole request frees the id. `/traverse` keeps checking expiry before its fields are parsed.
- *H1 (train 3):* two refusals of the request's intake come before every check of its fields and before
  admission, on every server and alike on a single-graph and a multi-graph one (never in an envelope): a body
  longer than `--max-request-body-mb` is refused **413** `{"error", "code": "body_too_large",
  "max_request_body_mb"}` before it is read (on every route; a `Content-Length` that is no byte count is closed
  without a response, `SPEC-labeled-traversal-core.md` §10.3), and a body of `/resolve`, `/pattern`, `/traverse` or
  `POST /traverse/cancel` with more JSON values than `--max-request-json-values` **400** `{"error", "code":
  "too_many_json_values", "max_request_json_values"}` before it is parsed (members in this order;
  `SPEC-labeled-traversal-core.md` §5 and §10.3). Neither holds an id: `GET /requests/{id}` answers 404 (never
  admitted), and the same id runs once the request is sent within the limits.
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
  the id is unknown. `GET /requests/{id}` answers 404 `{error, request_id, server_instance, server_time_ms}` for an
  id this server instance does not hold, in three cases: NOT ADMITTED YET (the window between the request's receipt
  and its admission), NEVER ADMITTED (refused by an HTTP-level 400 — `too_many_json_values` included — or a 413
  `body_too_large`, also by a 400 only the work finds, which frees the id
  when it is answered, or by the expired 409 or the instance_mismatch 409; or never received here), or ENDED AND DROPPED (past its retention, by
  time or by count). An id that breaks the token rule is 400. `reason`
  on `/resolve` and `/pattern`: `answered` (200 written), `client_gone` (nothing written), `deadline` (503),
  `error` (any other status but 400, which leaves no entry).
- The 404's `server_time_ms` (integer, ms since the epoch) is the server's clock, the one every admission
  compares `not_after_ms` with, read under the request registry's lock in the same critical section as the lookup
  that found no entry. An admission of an id reads that clock under a registry's lock, after its own lookup of the
  id: `/resolve` and `/pattern` under the request registry's, `/traverse` under the attempt registry's, which
  `GET /requests/{id}` asks after the request registry. An admission of the id after that lookup therefore reads
  `server_time_ms` or a later time, and one before it is found (200) unless it ended and was dropped, or a 400
  freed it. Admission refuses a time strictly greater than `not_after_ms` (section 4), so a 404 with
  `server_time_ms` > `not_after_ms` (strictly) proves that no request with that id and that `not_after_ms` will be
  admitted on this `server_instance` after the lookup: a copy that arrives later, or a re-send, is not admitted. The
  proof holds for a `request_id` of `/resolve` and `/pattern` and for an `attempt_id` of `/traverse` alike. It
  covers the `server_instance` that answered: a 404 from another one (a restarted server) proves nothing about the
  instance the request was sent to. A request sent with `expect_server_instance` X runs on X only: every other
  instance refuses it (`instance_mismatch`) whatever its time. A 404 from X with `server_time_ms` > `not_after_ms`
  therefore proves that no copy of that request is admitted anywhere after the lookup: X refuses it `expired`,
  every other instance on its instance. A 404 from another instance proves that this instance holds none and runs
  none; whether X still runs one is for the client to know (after a restart at the same address, X is gone). It assumes that the server's wall clock does not step back to `not_after_ms`
  or earlier between the lookup and a later admission: the proof survives a backward clock step SMALLER than
  `server_time_ms − not_after_ms` (a step of exactly that brings the clock to `not_after_ms`, which admission
  accepts; a client may require the difference to exceed `clock_skew_allowance_ms`). It says nothing about
  earlier runs: an id that ran and was dropped, or that a 400 freed, also answers 404. With `server_time_ms` <=
  `not_after_ms` the 404 proves nothing; the request may still be admitted. The field's presence is the marker: it
  is absent in a 404 from an engine without this guarantee, and a client then has no proof from that 404 (its own
  clock rule applies). The expired 409's `server_time_ms` (section 4, since release 10) is unchanged: the clock
  that admission read.
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

Each pair's discovery is bounded by the server per pair (`SPEC-labeled-traversal-core.md` §4.6: the labels it
tracks, a cut stated as the entry's `labels_capped`, and the memory of its pass, a `stop` with reason `memory`), and
a client merging the pairs' discoveries reads each entry's statements by this rule: **a pair entry with
`labels_capped` and `ranking_exact: false` makes a ranking merged over pairs a lower bound at best, never exact; with
`ranking_exact: true` that pair's top list is exact, and a label it did not track supports at most
`untracked_max_kmers` k-mers in that pair. A pair entry whose `stop` has reason `memory` or `time` is the exact answer
for its first `resolved_kmers` k-mers only; the rest of the sequence starts at `remainder_from_bp`.**

## 8. Capabilities (E3)

`GET /capabilities` on a multi-graph server adds `composite_fp` (sha256 over the pairs' `index_fp` sorted as bytes,
joined by `\n`, no trailing newline; null when any pair has none; each loaded (graph_path, annotation_path) pair
counted once, so a pair listed under two names counts once), `composite_rule` (that sentence), and a
`multigraph` block advertising `envelope_version`, the envelope routes, the selectors per route,
`traverse_pairs_per_request`, the routes reading `not_after_ms`, `finish_by_ms` and `request_id`,
`expect_server_instance` (F3b: `["/traverse", "/resolve", "/pattern"]`, the routes that take the field, section 5),
the registry route,
`unknown_labels`, `search_cuts`, `resolve_boundary_ties`, `discovery_trace_scale` (`"traced"`, section 7) and
`max_request_memory_mb` (E5b: the request-wide memory bound of `/traverse`, `--traverse-max-request-memory-mb`, an
integer in MiB, 0 when off; `SPEC-labeled-traversal-core.md` §5, request budgets), which every `/traverse`
capabilities probe (each entry of its envelope) states too, beside `max_memory_mb`: a graph-less probe sees it
without reading a pair. A seed the bound leaves unstarted is a result of the pair's entry (a failed result with
`resource_stop.scope: "request"`, phase `not_started`, its `unknown_labels` stated as every result of the pair
states them); the envelope is the one-entry envelope of a named `/traverse` as with any other result.
Its `graph_summary` states `headers_unique` (section 9).
From H2 (train 3) the block also states `resolve_discovery`: `{"max_tracked_labels": n, "max_memory_mb": n,
"max_answer_labels": n}`, the bounds of a `/resolve` discovery per pair (`--resolve-max-discover-labels`,
`--resolve-max-discover-memory-mb`, `--resolve-max-discover-answer-labels`;
integers, 0 for none; `SPEC-labeled-traversal-core.md` §4.6), as the `resolve` block's `discovery` states them, and
`resolve_max_explicit_labels` (an integer, 0 = unlimited): the cap on a `/resolve`'s explicit label lists
(`--resolve-max-explicit-labels`, `SPEC-labeled-traversal-core.md` §4.1), as the `resolve` block's
`max_explicit_labels`; above it the whole request is refused 400, before any pair. And `resolve_seed_order`
(`"rank"`, as the `resolve` block's `seed_order`; `SPEC-labeled-traversal-core.md` §4.2): the order of a `/resolve`
selection's seeds, the same on every pair (one engine). Its absence states nothing: builds before 01c05323b order
by position, 01c05323b and later builds without it already rank, so a client without it takes the order from the
build (rank if and only if 01c05323b is an ancestor).
From E19 the block also states `max_concurrent`: `{"/pattern": n, "/resolve": n, "/traverse": n}`, the caps on
the requests of each route the server runs at once (`--pattern-max-concurrent`, `--resolve-max-concurrent`,
`--traverse-max-concurrent`), integers, 0 for a route without one; above a cap the request is refused 503 `busy`
(section 1.1). Every `/traverse` capabilities probe (each entry of its envelope, and a single-graph server's probe)
states the same object beside its maxima.
From H1 (train 3) the block also states `request_intake`: what the server reads of a request, the same object a
single-graph server states at the top level of `GET /capabilities`:

```json
"request_intake": {
  "max_request_body_mb": 256,
  "body_too_large": {"status": 413, "code": "body_too_large", "members": ["code", "error", "max_request_body_mb"]},
  "max_request_json_values": 262144,
  "json_values_routes": ["/pattern", "/resolve", "/traverse", "/traverse/cancel"],
  "too_many_json_values": {"status": 400, "code": "too_many_json_values", "members": ["code", "error", "max_request_json_values"]},
  "max_connections": null
}
```

`max_request_body_mb` is `--max-request-body-mb` in MiB (`null` when unlimited; the top-level member of that
name, which stays in both modes, states the same), `max_request_json_values` is `--max-request-json-values`
(`null` when unlimited), `json_values_routes` the routes whose body is counted. `members` lists each refusal's
body members as a set (sorted); the body itself writes them `error` first, then `code`, then the cap (section 5;
`SPEC-labeled-traversal-core.md` §5 and §10.3). The 413 answers a body over `max_request_body_mb` on every route,
before it is read; the 400 a body of a listed route with more values, before it is parsed. Its presence is the
marker: an engine without it closes the connection on an over-limit body and counts no values.
From H3 (train 3) `request_intake` also states `max_connections`: `--max-connections`, the most TCP connections the
server holds open at once (`null` without a cap, the default and staging's, as in the example; an integer when
the operator sets one; at it a further connection waits in the listen queue until one closes, never refused;
`SPEC-labeled-traversal-core.md` §10.3, Connections). Its presence is the
marker: an engine without it never stops accepting. The block also states `column_labels`, what `GET
/column_labels` answers, the same object a single-graph server states at the top level of `GET /capabilities`:

```json
"column_labels": {
  "max_column_labels_bytes": 16777216,
  "paging": ["limit", "offset"],
  "too_large": {"status": 400, "code": "column_labels_too_large", "members": ["code", "error", "max_column_labels_bytes", "total"]},
  "total": 651
}
```

`max_column_labels_bytes` is `--max-column-labels-bytes` (`null` when unlimited): the most bytes of the labels' JSON
strings one answer holds; the unpaged list above it is refused 400 `column_labels_too_large` (`members` as for
`request_intake`'s refusals; the body writes `error` first). `paging` names the query parameters of a page
(`?offset=&limit=`, which answers `{"labels", "next_offset", "offset", "total"}`), `total` the length of the list
(`null` while a single index loads). On a multi-graph server the list holds each distinct (graph, annotation) pair
once, at its first place in the answers' order (section 3). Its presence is the marker: an engine without it
ignores the parameters and answers the whole array (`SPEC-labeled-traversal-core.md` §10.3, `GET /column_labels`).
A single-graph server, which states its request registry as `routes.requests`, states `expect_server_instance`
(the same list) at the top level of `GET /capabilities`, beside `routes.requests` (the pattern of `headers_unique`:
one reader serves both modes), and `request_intake` there too. An engine without the key refuses the field 400 (unknown field), so clients send it
only where the key lists the route.
Clients switch on these, never on dates.

### 8.1 The capability profile of a pair (E17)

Each `graph_summary.pairs[]` entry (`SPEC-pattern-search.md` §24.4) states, besides the pair, its identity, k,
the pattern fields and `traversal`, the pair's capability profile, computed once at start-up before the server
listens:

| Member | Type | Meaning |
| --- | --- | --- |
| `index_meta_fp` | string | the pair's meta fingerprint, the value its `/traverse` and `/resolve` capabilities state (`SPEC-labeled-traversal-core.md` §10.3): one pass over the pair's column names, the pairs in parallel, logged with its time (`index_meta_fp of the N pairs (C columns) computed in S s`). A negative check only, never identity: equal means "index shape unchanged", never "file unchanged" |
| `alphabet` | string | the graph's alphabet as `/traverse` capabilities state it (a DNA build's succinct graph `$ACGT`, its hash graphs `ACGT`; a Protein build's succinct graph `$ABCDEFGHIJKLMNOPQRSTUVWYZX`) |
| `graph_representation` | string | the graph's type as `--graph` names it: `succinct`, `hash` (an ordered hash graph, packed or not: one file extension), `hashstr`, `hashfast`, `sshash` or `bitmap`; of a PRIMARY graph (wrapped in `CanonicalDBG`) the wrapped graph's |
| `annotation_form` | string | the annotation's type as `--anno-type` names it: the class the loader constructed from the file's extension, e.g. `row_diff_brwt_coord` (RefSeq r237, AllTheBacteria, UniParc), `row_diff_int_brwt` (Logan's `*_RNA` chunks), `row_diff_disk`, `column` |
| `decode_charged` | bool | whether a traversal under a memory budget charges the pair's annotation reads before they are held (`SPEC-labeled-traversal-core.md` §6.8): the traversal oracle's own predicate on the matrix. True on the row-diff forms over an index matrix with a single-row read (`row_diff`, `row_diff_coord`, `row_diff_brwt`, `row_diff_brwt_coord`, `row_diff_disk`, `row_diff_sparse`, `row_diff_flat`, `row_diff_int_brwt`) and on plain `brwt` and `column`; false on `column_coord`, `brwt_coord`, `int_brwt`, `row_diff_int_disk`, `row_diff_disk_coord`, `flat`, `row_sparse`, `row`, `bin_rel_wt`, `rbfish` and `rb_brwt`. Plain `brwt` and `column` stay `annotation: "unbudgeted"` on `/pattern` (budgeted there is a row-diff form only): the two fields together state that nuance |
| `traversal.counts` | string or null | `"presence_only"` on a count annotation (`row_diff_int_brwt`, `int_brwt`, `row_diff_int_disk`: walks read presence, `SPEC-labeled-traversal-core.md` §10.3), null on the others, so that a client without a graph sees it per pair |
| `support_fp` | string | FNV-1a-64, 16 lowercase hex digits, over the pair's capability fields (below) |

`support_fp` fingerprints what the routes can do on the pair, not its content. Two pairs with equal values agree
on the graph representation, k, regime, alphabet, annotation form, coordinates, record mapping (`.seqs`), trace,
charged reads, counts and the pattern fields (`graph_mode`, `available`, `unavailable_reason`, `mask`,
`counting`), so that a service probes the routes once per distinct `support_fp` (the chunks of one Logan build
share it) and keeps outcomes that depend on the data (a seed or a label absent from one pair) per request.
`graph_representation` stays in `support_fp` although no answer of `/traverse`, `/resolve` or `/pattern` depends on
it beyond what `alphabet` and `available` already state: pairs that run different code are probed apart (a wrong
split costs one probe, a wrong merge carries one pair's outcome to another). The
identity (`graph`, `graph_path`, `annotation_path`, `index_ns`, `index_fp`, `index_meta_fp`) and the size
(`traversal.num_labels`) are left out; the pattern block's `dummy_fraction` (sampled from the graph's data) is not
a member of the entry. The hashed text is a version line, then one line per field in this order (byte order of
the names; a member of `traversal` named `traversal.<member>`), each line ending with `\n`:

```
mgt-support 1
alphabet=<value>
annotation_form=<value>
available=<value>
counting=<value>
decode_charged=<value>
graph_mode=<value>
graph_representation=<value>
k=<value>
mask=<value>
traversal.counts=<value>
traversal.has_coord_to_header=<value>
traversal.has_coordinates=<value>
traversal.regime=<value>
traversal.supports_trace=<value>
unavailable_reason=<value>
```

`<value>` is the entry's member as compact JSON, as Python's `json.dumps` writes it: a string in double quotes
with JSON's escapes (non-ASCII as `\uXXXX`; every value is ASCII today), `true` or `false`, `null`, an integer in
decimal. A client recomputes `support_fp` from the entry alone: `'%016x' % fnv1a_64(text.encode())`, FNV-1a-64
with offset basis 14695981039346656037 and prime 1099511628211 (as `index_meta_fp`). The first line is the rule's
version: a capability field added to the profile later changes the version and the list, and with them every
pair's `support_fp`, so that values of servers of different versions differ even for one profile, while equal
values (up to FNV-1a-64 collisions) always mean one profile of one rule. A client compares `support_fp` within
one server's document or between servers of one version, and never keeps it as the identity of an index.

A single-graph server states its pair's profile at the top level of `GET /capabilities` as `pair_summary`: the
same entry, with `graph`, `graph_path` and `annotation_path` null (as its `headers_unique` tags its pair) and
`index_ns`, `index_fp` from `--index-name` and `--index-manifest`; null while the index loads (`ready: false`). A
multi-graph server has no top-level `pair_summary`. The server logs the number of distinct profiles at start-up
(multi-graph: `N capability profile(s) (support_fp)`; single-graph: `Capability profile: <representation>
<form> annotation, support_fp <fp>, index_meta_fp <fp>`).

Reference instance: sample 12, whose pairs' profile members are recorded on eleven fixture pairs of r237's
profile (succinct, k 31, basic, `row_diff_brwt_coord` with its `.seqs`, no mask) served under r237's names, their
tags rewritten to sample 12's paths: `support_fp` is r237's (a client recomputes it from each entry),
`index_meta_fp` the fixture pairs' (r237's own are read from staging's `GET /capabilities` after the deploy).

## 9. headers_unique (E22)

Whether every sequence header of a server names one record, so that a header label (a `.seqs` record name) selects
one sequence of one pair. Computed once at start-up from the pairs' header indexes (a merge of their sorted hash
arrays; equal hashes are confirmed on the header text, so a hash collision is never counted).

- Multi-graph server: `GET /capabilities` `graph_summary.headers_unique`, over every pair as listed, each once: a
  (graph_path, annotation_path) as the graph list spells it, as `composite_fp` (section 8) and
  `graph_summary.columns_disjoint` count pairs. A pair listed under several names counts once and is tagged with the
  first name in pair order; two spellings of one pair's paths are two pairs, whose headers then all repeat. There is
  no top-level member.
- Single-graph server: `GET /capabilities` `headers_unique` at the top level, over its one pair, with null pair
  tags (as its `/search` cut); null while the index loads (`ready: false`).

`{"value": true | false | null, "collisions": n, "examples": [{"header", "occurrences", "pairs": [{graph,
graph_path, annotation_path, column}...]}...]}`:

| Member | Type | Meaning |
| --- | --- | --- |
| `value` | bool or null | true: no header occurs twice; false: one does; null: a pair has no `.seqs` loaded (none beside its annotation, or the server runs with `--no-coord-mapping`): its labels are columns only, so the server cannot certify header labels |
| `collisions` | int | the number of distinct headers that occur more than once: in two pairs, in two columns of one pair, or twice in one column. With `value` null it counts over the pairs that have a `.seqs` (a repeat among them is still reported) |
| `examples` | [object] | up to 10 repeated headers, the byte-smallest ones (independent of the hash) |
| `examples[].header` | string | the header |
| `examples[].occurrences` | int | its number of records over the pairs |
| `examples[].pairs` | [object] | its first 10 occurrences in pair order, then column order, then record order: `{graph, graph_path, annotation_path, column}`, `column` the column's label, as given to annotate; a header twice in one column lists that column twice |

A header that occurs more than once is answered by its first occurrence: within a pair, the first in column order
(`find_header`); a fan-out over pairs answers each pair's own first occurrence. `value: false` tells a client that a
header label alone may not name the record it means. The server logs the certificate at start-up (`headers_unique:
true | false | null`), a repeat as a warning naming one header. The certificate is all that E22 serves: `/search`
rows keep their shape (no (pair, column, record) per row); such a fallback is later work, built only if a deployed
index's certificate comes out false.

Reference instances: sample 12 (illustrative: `false` on a multi-graph server, one header in two pairs) and the
pattern fixtures (`api/python/tests/data/traverse/pattern/`, `SPEC-pattern-search.md` §11) `capabilities`
(single-graph, `true`), `capabilities_repeated_headers` (single-graph, `false`: one header in both files of the pair,
its occurrences with null pair tags), `capabilities_representation_unsupported` (`null`: the pair has no `.seqs`) and
`capabilities_multi_graph` (`true`, in `graph_summary`).
