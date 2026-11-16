# /search partial refusal per (record, pair): design note

2026-11-14, revised 2026-11-15 (phase 1b). Lane T7-P1: a design note, no code. Base: staging `gr/labeled-traversal`
at 9726cb420 (train 5b); phase 2 stacks on train 6 lane W2's head (9e7537fed or its landed form). Owner decision
2026-11-14: option (c), a partial answer in one request; "it's fair that the service can refuse some sequences, but
it needs to be properly reported back with the search results". Fable (xhigh) reviewed the first version
(8f1baf3ae): "implement after fixes", no determinism, memory or byte-identity counterexample. The coordinator and
the service then settled the open questions (E23/T7/P1/DECISIONS.md). Section 11 lists each review finding, how I
verified it against the code and how it is settled. Remaining uncertainties are **[U1]** to **[U4]** (section 10).

## 0. Summary

- **Opt-in.** A request with `"partial_refusal": true` gets refused units reported in its answer (status 200).
  Without the field, any refused unit refuses the whole request 400 `search_answer_too_large`, as in 5b, and the
  body gains `account`. The field joins the result cache's key. Capabilities `search` state it.
- **The unit** is one FASTA record on one pair. A unit is refused when its values (its rows and per-position values)
  or its own JSON object would pass `--search-max-answer-memory-mb`. The rule is the same in both modes; the field
  only decides how a refused unit is reported.
- **With the field**, a refused unit's existing `cuts` entry gains `refused: {account, code, max_memory_mb}` (and
  `priced_mb` when `account` is `"answer"`), states `returned` 0 and `complete` false, and its rows are absent.
  Every other unit's rows and cut entry are byte-identical to today's answer. Every unit refused still gives 200.
- **The request's answer account** (all record objects together, then the text) stays per request. Passing it
  refuses the whole request 400, now with `account: "answer"`.
- **No request-wide values total** (accepted, Q1). Each record's JSON object is built as soon as its unit completes,
  and its values are freed then. No completed unit's values wait, so a total over them has nothing to bound. The
  values held are, in the model, at most (pair searches at once) x (records at once per pair) x the bound
  (section 5; model totals, not hard bounds).
- **One flag**, unchanged default 4096, no new flag. Capabilities state the concurrency factors the model total
  depends on (`-p`, `--threads-each`, the HTTP handler threads).
- **Deterministic.** Which units are refused, and the status, body and `account` of every answer, are functions of
  the request, the index and the server's flags, never of thread order (section 4).
- **Staging** (`-p 2 --threads-each 2`, at most 2 requests at once): a model total of 16 GiB of values, 8 GiB of
  answers and 0.5 GiB of bodies, 24.5 GiB in all, on the 110 GiB the container has beyond its ~18 GiB at rest.
- **Compatibility.** Every request answered today is answered byte-identically in both modes. A request refused
  today only because one pair's records together passed the bound, with each record under it, is now answered.

## 1. Problem

Train 5b (9726cb420, on staging since 2026-11-14 10:36Z) bounds what one `POST /search` holds
(SPEC-multigraph-envelope.md section 3, "The answer budget") with two kinds of account:
- one **values account per pair**: the rows of every FASTA record of the body on that pair and their per-position
  values, held until the pair's answer is built;
- one **answer account per request**: each record's JSON object, then the text with the transport's (or gzip's)
  copy and, with `--cache-size` above 0, the cache's copy.

Either account passing refuses the whole request 400. Nothing is sent for the other pairs or records. On 2026-11-14
the suite's run and the service's default shape (11 pairs, `query_coords` added automatically, `top_labels` 1000,
both strands as two records, a 300 bp slice) were refused on r237 because one eukaryotic pair's values passed 4096
MiB. Every such slice was a low-complexity chromosome end (84-100% telomere repeat). With query_coords every position
keeps each label's occurrences: about 1,500 x 270 x 8 B x 1,000 rows, roughly 3 GB per record and pair. A unique
300 bp query holds about 2 MB. The refusals are right but too coarse: one pathological (record, pair) costs the
caller the other 21 units' answers. The owner keeps refusals for pathological queries and wants them reported with
the results, in one request.

## 2. Contract (exact)

### 2.1 The request field

`partial_refusal`: a boolean, false when absent, read as /search's other boolean fields (`with_signature`,
`query_coords`). It joins the result cache's key (`make_search_cache_key`, server_cache.cpp, beside `query_mode`
and `top_labels`), so a request with it and one without never share a cached answer. /search refuses no unknown
field (only `pairs`), so a 5b engine ignores the field and answers 5b's 400. A client reads `search.unit_refused`
in capabilities (2.6) to know whether the engine reads the field.

### 2.2 The refused cut entry (with `partial_refusal`)

Section 3 of the spec already gives each record's object `cuts`, one entry per pair searched, in pair order:
`{graph, graph_path, annotation_path, index_fp, returned, complete, nth_score, tied, tied_returned}`, plus
`nth_median` and `tied_same_median` where the cut has abundances. The top-N cut is computed before the first values
charge of every path (`search_top_n` writes it, annotated_dbg.cpp:543-545, called through `top_header_matches` at
1274 and 1650, `filter_and_aggregate` at 1387, `get_top_labels` at 1766; the values charges follow at 1298, 1410,
1667, 1775 and in W2's store scratch). So a refused unit's cut is known, and the entry takes its member set from it.
The entry states what was sent, which is no row:

| Member | Refused entry | Why |
| --- | --- | --- |
| `graph`, `graph_path`, `annotation_path`, `index_fp` | the pair's, as on every entry (null tags on a single-graph server) | identifies the pair; the record is the object the entry sits in |
| `returned` | `0` | no row of the unit is sent |
| `complete` | `false` | `complete: true` with `returned` 0 states "nothing matched". A refused unit's matches were found but not sent. A client that ignores `refused` reads "matches exist that were not returned", which is true (Q3) |
| `nth_score` | `null` | no N-th row was sent |
| `tied`, `tied_returned` | `0` | as in the N = 0 state: no row was sent |
| `nth_median`, `tied_same_median` | `null`, `0`, exactly when the unit's computed cut has abundances | the member set is the one the unit would have stated (it is the same for every entry of a pair and mode) |
| `refused` | `{account, code, max_memory_mb}`, plus `priced_mb` when `account` is `"answer"` | present only on a refused entry. An answered entry has no `refused` member (not `null`), so it stays byte-identical |

`refused`:

| Member | Value |
| --- | --- |
| `account` | `"values"`: the unit's rows and per-position values passed the bound while the record was searched on the pair, or as its rows were handed over. `"answer"`: the record's object for this pair was priced past the bound before it was built (3.1, item 2) |
| `code` | `"search_answer_too_large"`, the whole-request refusal's code, so one recogniser serves both |
| `max_memory_mb` | the flag's value, an integer |
| `priced_mb` | only with `account: "answer"`: the unit's own object (its tree plus its own transient, `search_record_price`) in MiB, rounded up, the sum saturating at UINT64_MAX before the division. Never an increment over the pair's high-water mark. A function of the unit's result alone. There is none for `values`: that account stops at its first charge past the bound, which says nothing useful about the unit's need (Q5) |

No member repeats the record's header. The record is the object the entry sits in: its `seq_description` and its
position in the answer array.

### 2.3 Example: multi-graph, two records, two pairs, one unit refused

Request `{"FASTA": ">fwd\n<300 bp>\n>rc\n<300 bp>", "graphs": ["a", "b"], "query_coords": true, "top_labels": 1000,
"partial_refusal": true}`, where record `fwd` on pair `b` passes the bound. The answer is shown compact. The server
writes it as it writes every /search answer: the same writer, members in byte order of their names.

```json
[
 {"cuts":[
   {"annotation_path":"/idx/a.row_diff_brwt_coord.annodbg","complete":false,"graph":"a","graph_path":"/idx/a.dbg","index_fp":"<a's index_fp>","nth_score":270,"returned":1000,"tied":41,"tied_returned":3},
   {"annotation_path":"/idx/b.row_diff_brwt_coord.annodbg","complete":false,"graph":"b","graph_path":"/idx/b.dbg","index_fp":"<b's index_fp>","nth_score":null,"refused":{"account":"values","code":"search_answer_too_large","max_memory_mb":4096},"returned":0,"tied":0,"tied_returned":0}],
  "results":[ <the 1000 rows of pair a, exactly as today> ],
  "seq_description":"fwd"},
 {"cuts":[ <pair a's entry, as today>, <pair b's entry, as today> ],
  "results":[ <the rows of pair a, then of pair b, as today> ],
  "seq_description":"rc"}
]
```

Byte identity, stated exactly: in every record object, each answered unit's rows (a contiguous run in `results`, in
pair order) and its cut entry are the bytes of today's answer. The record-level members (`seq_description`, and with
`align` the first pair's `sequence`, `score`, `max_score`, `cigar`, `orientation`) are today's too, also when the
first pair's unit is refused: the alignment is computed before the search and is attached to the refused stub
(3.2). Only the refused unit's rows are missing and only its cut entry differs.

On a single-graph server the rule is the same with the one pair: each record of the body is its own unit, and a
refused record's object has `results: []` and one cut entry with null tags and `refused`.

### 2.4 The whole-request refusal

Status 400, `Content-Type: application/json`, no `Content-Encoding`, no `Retry-After`, compact, `error` first, then
the other members in byte order. The text is unchanged, so the service's recogniser and its displayed advice still
apply. One member is added:

```
{"error":"search: the answer would hold more than --search-max-answer-memory-mb 4096 MiB: nothing partial is sent; ask for fewer top_labels, a larger discovery_fraction, fewer or shorter FASTA records, or no per-position values (with_signature, query_counts, query_coords)","account":"answer","code":"search_answer_too_large","max_memory_mb":4096}
```

`account` is `"answer"` or `"values"`, and is a function of the request (section 4):
- **With `partial_refusal`**, only the request's answer account refuses the whole request, so `account` is always
  `"answer"`.
- **Without it**, the 400 names `"values"` if some unit's values account passes, else `"answer"` (a unit's object
  check or the request's answer account). Values come first because they are the refusal seen in practice, and
  this order lets the request stop early in that case:
  - The first unit refused on its values fixes the answer as 400 `"values"`, and every running unit stops at its
    next charge, as in 5b.
  - A unit refused by its object check, or the request's answer account passing, fixes the status (400) but not
    yet the account. From then on no object is built or charged, and the units not yet finished run their
    searches only, charging their own values accounts and freeing their values at once. If one of them passes,
    the 400 names `"values"`; once all have run, it names `"answer"`. This costs the remaining searches in that
    (rare) case. 5b had no answer-account refusal on staging; all 26 of 2026-11-14 were values.

There is **no `priced_mb`** in the 400. Several pairs charge the request's answer account at once, so the bytes held
at the refused charge depend on thread order.

### 2.5 Status and what a client concludes

| Today (9726cb420) | P1 without `partial_refusal` | P1 with `partial_refusal` |
| --- | --- | --- |
| 200 | 200, byte-identical | 200, byte-identical |
| 400, a pair's values passed | 400 `"values"` if a unit's values pass; 200 if every unit fits (the pair passed only with its records together); 400 `"answer"` if the answer account or an object check passes | 200 with the units over the bound refused (all units refused: still 200); 400 `"answer"` if the request's answer account passes |
| 400, the answer account passed | 400, `account` by the rule of 2.4 | 400 `"answer"`; or 200 with the units whose own object passes the bound refused, if what remains fits |
| any other refusal or error | unchanged (precedence in 3.3) | unchanged (precedence in 3.3) |

Why the first row holds: a unit's account is charged at most what its pair's was (W2's scratch: the unit's own, at
most the pair's `threads` x the largest), and its object's price is at most the request's answer account. So no unit
passes where its pair and its request did not. The answer account's charges are those of 9726cb420 plus W2, in
another order (3.2).

The answer is an array (the upstream shape), so there is no top-level summary. A client counts the refused entries
in `cuts`. A refused entry never means "nothing matched": the caller did not get this unit's rows. It may retry the
record alone on the pair, with smaller `top_labels`, a larger `discovery_fraction` or without per-position values.
The same request on the same server and index is refused again in the same units: the outcome is final. The service
will not store a refused cut as complete (DECISIONS, Q9).

### 2.6 Capabilities

`search` (top level on a single-graph server, `multigraph.search` on a multi-graph server, as today) becomes:

```json
"search": {
  "handler_threads": 2,
  "max_memory_mb": 4096,
  "pair_searches_at_once": 2,
  "records_at_once": 2,
  "rule": "SPEC-multigraph-envelope.md section 3, the /search answer budget",
  "too_large": {"status": 400, "code": "search_answer_too_large", "accounts": ["answer", "values"], "members": ["account", "code", "error", "max_memory_mb"]},
  "unit_refused": {"field": "partial_refusal", "code": "search_answer_too_large", "accounts": ["answer", "values"], "members": ["account", "code", "max_memory_mb", "priced_mb"]}
}
```

- `unit_refused` is new, and its presence is the marker: the engine reads `field` and refuses units alone.
  `accounts` and `members` are sets, sorted; `priced_mb` is present only with `account: "answer"`.
- `too_large` gains `accounts` and the member `account`.
- `pair_searches_at_once`: how many pair searches the server runs at once, over all requests. On a multi-graph
  server this is the graphs' pool size (`-p`): every /search pair runs on that pool, and the pool is shared with
  the fan-outs of `/resolve` and `/pattern`. On a single-graph server /search runs on its HTTP handler thread, so it
  equals `handler_threads`.
- `records_at_once`: how many records one pair search runs at once (`--threads-each`: query_fasta's pool, or the
  batched path's threads per batch, one batch at a time).
- `handler_threads`: the HTTP handler threads (`--http-threads`, else `-p`), the most requests of any route the
  server handles at once. /search has no cap of its own, so at most this many /search answers are held at once.
- None is part of any fingerprint. `max_memory_mb` stays null when the flag is 0; the other members are stated even
  then.

Pins that change: the `search` block's member set wherever it is pinned (tests/cli/test_server_envelope.cpp sample
12's expected block, integration_tests/test_traverse.py's top-level set, integration_tests/test_search_budget.py,
api/python/tests/test_pattern_fixtures.py, tests/cli/test_search_answer.cpp's block bytes), the Linux capability
fixtures (api/python/tests/data/traverse/pattern/capabilities*/answer.json), contract samples 12 and 37, and the
suite's `--expect-member` lines. Three members now depend on server flags, so every recorder and fixture generator
starts its server with `-p`, `--threads-each` and `--http-threads` set explicitly, or the pins vary with the host
(review finding 7).

## 3. Semantics

### 3.1 Accounts

`--search-max-answer-memory-mb` (MiB, default 4096, 0 = no bound) is the bound B of three kinds of charge, each made
before what it prices is built (the model of the deployed build's heap, as today, never a measurement):

1. **A unit's values account**, one per (record, pair). It is charged what today's pair account is charged, for
   that one record: its per-position values as their buffers grow, the store its column-label search gathers its
   rows in (W2's `charge_scratch`; a unit is one task, `threads` 1), and its rows as its result is handed over
   (`search_record_rows_bytes`, now made in `query_sequence`, 3.2). If it would pass B, the unit is refused with
   `account: "values"`. With `partial_refusal` the request is not marked and the other units go on; without it,
   see 2.4.
2. **A unit's object check.** Before the record's object for this pair is built, its own price (tree and
   transient, `search_record_price`, with W2's scratch for measuring it charged to the request first, as today) is
   compared with B. If the object alone would pass B, the unit is refused with `account: "answer"`, and only the
   refused stub's price is charged to the request. With one flag, an object over B would pass the request's
   account anyway. So this check never refuses what 5b answers; with `partial_refusal` it turns a whole refusal
   into a partial one.
3. **The request's answer account**, one per request, as today: each unit's object or refused stub, before it is
   built; the per-pair transient high-water mark (objects of one pair are built one at a time, 3.2); then the
   text, the transport's or gzip's copy and the cache's copy. If it would pass B, the whole request is refused 400
   (2.4). With `partial_refusal`, every unit still running stops at its next charge, as today.

Not charged, as today: the rows decoded before the top-N cut (an accepted risk, owner 2026-11-14 11:35Z;
SPEC-multigraph-envelope.md section 3), the index (mmap pages, an `in_ram` load), the body (its 256 MiB cap), the
query graph of the batched path. Section 5.1 lists what the model misses.

### 3.2 The pipeline: one refusal site per record, and the object built when the unit completes

Today `process_search_request` (cli/server.cpp) collects every record's `SeqSearchResult` of the pair, sorts them by
record id, then builds the objects. Each record's values are held until the pair's last record is done, and the
pair's one account is charged in the callback. P1 changes this as follows:

- `query_fasta` and `batched_query_fasta` take the request's `SearchBudget *` (null on the CLI) in place of the
  pair's account. Each record's task makes its unit account on its own stack and passes it down. The result's
  values never outlive the task: they detach at `shrink()`, as today, and the callback consumes the result before
  the task returns.
- **`execute_query`** (query.cpp:503-578) catches the unit-scoped `SearchUnitRefused` (a new type, not derived from
  `SearchAnswerTooLarge`) from the mode's call. The cut is a local there and complete (2.2), and the sequence is
  moved only at the return (query.cpp:577). It returns the refused stub: an empty row vector of the mode's type and
  the cut marked refused.
- **`query_sequence`** makes the handover charge (`search_record_rows_bytes`) on the unit account right after
  `execute_query` returns, unless the result is already a stub, so the sticky account never refuses twice. On
  refusal it turns the result into the stub (rows dropped, cut kept and marked). It then attaches the alignment in
  every case, so a stub keeps it. These two places are the only refusal sites of a unit's values.
- **The callback** needs no account. Under the pair's mutex it makes the object check (3.1 item 2), charges the
  request's answer account (or the stub's price), builds the object (rows tagged, cut entry) and keeps
  `(record id, object)`. When it returns, the record's result and its values are freed.
- After `query_fasta` returns, the pair's objects are sorted by record id and appended, which gives today's order
  and bytes. The pair's mutex keeps "one object of a pair at a time", the assumption behind the per-pair transient
  high-water mark (W2's `charge_search_record`).

The request-scoped `SearchAnswerTooLarge` (the request passed) is never caught by these sites; it propagates as
today. At most `pair_searches_at_once x records_at_once` unit accounts exist at once.

### 3.3 Precedence and errors

- The refusals that need no work come first, unchanged: 413 `body_too_large`, the body's JSON, the fields,
  `pairs_unsupported`, `graph_unknown`, `pair_not_found`, then a cached answer.
- **Without `partial_refusal`**, a refused unit refuses the request whatever else the work raised (another pair's
  error, a load refused `index_changed`), as in 5b. The `account` follows 2.4.
- **With it**, a refused unit is not an error. A pair whose units are all refused is answered with refused entries.
  The request's answer account passing still wins over every other error, as today. Otherwise errors are answered
  as today: the first failing pair in pair order. So, with the field, a request in which one unit is refused and
  another pair errors (for example `Annotation does not support k-mer coordinate queries`) is answered with that
  pair's error, where 5b gave its budget 400. Both are deterministic.

### 3.4 Interactions

- **Result cache.** The field is in the key (2.1). A partly refused answer is a 200 and is deterministic for the
  process (its flags, `--cache-size` and `--verbose-output` included, are fixed), so it is cached like any 200:
  STORE, HIT, duplicate waiters served the same bytes. Every 400 stays uncached, as today. The cache's copy stays
  charged to the answer account when `--cache-size` is above 0 (staging runs 0).
- **`top_labels`, `discovery_fraction`.** They set how many rows a unit keeps and which labels qualify, so they set
  its values. The advice is the 400's.
- **`in_ram`.** No change: the load is reserved (reservations), not charged. The same charges give the same
  refusals with and without it.
- **Time and the shared pool (review finding 1).** /search reads no deadline (`not_after_ms`, `finish_by_ms` and
  `request_id` belong to /traverse, /resolve and /pattern, spec section 4) and has no busy cap; the 899 s transport
  cap bounds it, as every answer.
  - In 5b, a values refusal sets `passed`, and every pair stops at its next charge (search_budget.hpp:122-129).
    Without `partial_refusal`, P1 keeps that early stop for values (2.4).
  - With `partial_refusal`, a refused unit stops only itself, and every other unit runs to completion. A partly
    refused request is therefore longer than 5b's refusal of it, and it holds the graphs' pool threads that
    /resolve and /pattern share for that time.
  - When the request's answer account passes late, the 400 comes after most units ran. This is the case of the
    telomere request (2 records x 11 pairs, query_coords, `top_labels` 1000): the suite measured 0.4-0.56 GiB of
    JSON per pair and strand, about 5 GiB in all, over the 4 GiB bound. So that request still ends in a
    whole-request 400 (`"answer"`) with the field, after most of its work.
  - The service accepts this. Its per-record fallback then sends each record alone (about 2.5 GiB of answer each),
    with vo_sauropsida's unit refused in its cut.
  - These figures are the suite's and the model's. I made no local measurement (the budget did not allow one).
- **`align`.** The alignment is computed before the search and is not charged to values until the handover. A
  refused unit keeps it (2.3).
- **The batched path** (a pair without `.seqs`, not query_coords: Logan's chunks) runs each record of a batch as one
  OpenMP iteration. The unit rule is the same.
- **`--verbose-output`** changes prices, and so which units are refused. It is a server flag, fixed per process.

### 3.5 The CLI

`metagraph query` has no budget and keeps none. The new `SearchBudget *` parameter is null there, so the CLI's
output and pymetagraph's local paths are unchanged.

## 4. Determinism

Claim: which units are refused, with which `account` and `priced_mb`, and the status and body of every answer, are
functions of (the request, the index, the server's flags). They never depend on the order in which threads, records
or pairs run or charge.

1. **A unit's values account** is charged only by its own record's work on its own pair, and only grows. Whether
   some charge passes B depends only on the total of its charges, which is a function of the unit. This holds
   whatever thread makes each charge. No member reports the charge that passed (`priced_mb` is not stated for
   values), so the result does not depend on the order of a unit's own charges.
2. **A unit's object check** compares a price computed from the unit's result alone. The same holds for
   `priced_mb`.
3. **The request's answer account** only grows. It receives a fixed set of charges:
   - one object or stub per unit, fixed by items 1 and 2;
   - the per-pair transient high-water marks, each the maximum over its pair's records of the transient and of W2's
     measuring scratch. The scratch counts for units refused by the object check too, because it is charged
     before the price is known. The maximum is the same in every order (review finding 8);
   - the text charges, which come after every unit.

   The total is therefore fixed, and every order's prefix is at most that total, so the account passes in every
   order or in none.
4. **With `partial_refusal`**, the only whole-request budget refusal is item 3, so the 400 always names `"answer"`.
5. **Without it**, the 400 names `"values"` exactly when some unit's values account would pass (item 1: a property
   of the request). Otherwise it names `"answer"`, when an object check or item 3 passes. The early stop happens
   only on a values refusal, whose name is already fixed. After an answer-side refusal, every remaining unit's
   values account still runs to its end (2.4), so whether some unit would pass is always learned, never cut off by
   thread order.
6. **No account is shared between units.** A unit's charges never read another unit's account or the request's,
   beyond "the request is already refused", which matters only when the answer is the 400 anyway.
7. **No outcome reads memory in use**, the clock or the load. The model totals of section 5 follow from the
   structure; they are not an admission rule.
8. **Errors** (3.3): a function of the field and of the pairs' own errors, as today.

The brief's request-wide values total (section 5.4) would refuse a fixed set, but it could not name its account
deterministically. That is one reason it is rejected.

## 5. Model totals and defaults

### 5.1 The model total, and what it misses

Let B = the flag (4096 MiB), P = `pair_searches_at_once`, E = `records_at_once`, H = `handler_threads`, C = /search
requests at once (at most H; the service sends at most 2 to refseq-experimental). The model totals are:
- **values**: P x E x B over all requests. A running unit's account stops at B, a unit's values are freed once its
  object is built, and at most P x E units run at once (a unit waiting for its pair's mutex still holds its thread);
- **answers**: C x B, each request's answer account;
- **bodies**: H x 256 MiB (`--max-request-body-mb`).

These are totals of the model, not bounds on memory (review finding 6). What the model misses:
- the rows decoded before the top-N cut, per running unit, unbounded (accepted risk);
- **plain modes**: rows without per-position values are built uncharged and charged once at the handover.
  `top_labels` has no cap (server_request.cpp:35), so a unit can hold more than B before it is refused;
- **responses still being written** after the handler returns (`response->send`, server_utils.cpp:1082-1087): a
  finished request's transport copy is held beside the next request's accounts;
- **the body in several copies**: the request's content string (server_utils.cpp:1061), the parsed JSON,
  `SearchRequest.fasta`, and on disk a temp file per pair;
- **gzip's output** is appended without a reserve (`compress_string`, server_utils.cpp:53-112), so incompressible
  text can take more than the charged copy;
- allocator slack. Y4 measured refused requests peaking slightly above the bound.

### 5.2 Staging (refseq-experimental, r237, 128 GiB container, about 18 GiB anonymous at rest)

Flags: `-p 2 --threads-each 2 --mmap --madv-random --cache-size 0`, no `--http-threads`, so H = 2 ([U1]).

| Model total | P1 | 9726cb420 (today, same reasoning) |
| --- | --- | --- |
| values | 2 x 2 x 4 GiB = 16 GiB | P x B = 2 x 4 GiB = 8 GiB (one pair account per pool thread) |
| answers | 2 x 4 GiB = 8 GiB | 8 GiB |
| bodies | 2 x 0.25 GiB = 0.5 GiB | 0.5 GiB |
| total | 24.5 GiB | 16.5 GiB |

That leaves about 85 GiB of the 110 GiB for what section 5.1 lists.

### 5.3 The suite's limits check, restated

The 5b check (4096 x 11 + 4096 + 256 MiB = 48.3 GiB a request, 96.5 GiB for two) assumed that every pair of a
request runs at once, independently of the pool. With P1, and from the capabilities alone:

    pair_searches_at_once x records_at_once x B  +  handler_threads x (B + 256 MiB)  <=  container - anonymous at rest

At the default B with H = 2 the check holds while P x E <= 25 (110 GiB). At P x E = 22 it gives 88 + 8.5 = 96.5
GiB, the 5b figure for two requests. Staging's P x E is 4. The default stays 4096 (owner 2026-11-14: "keep 4,096"),
now meaning 4096 MiB per unit; no unit answered today holds more (2.5).

### 5.4 Why no request-wide values total (accepted, Q1)

The brief proposed keeping today's pipeline (a pair's records held until its answer is built) and adding a request
total T, charged with each unit's final size. Its problems:
1. **T sums every unit of the request, not the units held at once**, because an account that only grows cannot
   release a pair's units. A request answered today can hold up to 11 x 4 GiB over its pairs, so any T below 44
   GiB refuses some request that is answered today, and T >= 44 GiB protects nothing.
2. **The 400 could not name its account deterministically**: T and the answer account are charged at the same time
   by different pairs.
3. **A second default** to agree.
4. **It holds more**: completed units' values wait for their pair's last record.

P1 removes the reason for T. A server with a large P x E holds more under P1 than under 5b (P x B). That is a sizing
matter, which the limits check now sees. An admission rule on memory in use would break determinism.

## 6. Code plan (phase 2, on W2's head)

- **graph/search_budget.hpp**:
  - an account's scope (the request's answer, or a unit's values);
  - a unit account that would pass marks itself refused (sticky) and throws `SearchUnitRefused{account}` without
    marking the request;
  - the request's state: open; refused with a fixed account; or, without the field only, "refused, account
    pending", in which unit values charges go on but no object is charged (2.4);
  - `charge()` checks the request's state first;
  - saturating add wherever held + charge is formed (W2's `charge_scratch` can pass UINT64_MAX; the log's sum,
    `priced_mb`);
  - the budget counts refused units for the log and keeps `values_max` per unit. W2's `threads` stays (units pass
    1).
- **graph/annotated_dbg.{hpp,cpp}**: no change in what is charged. The `SearchAccount *` they receive is the unit's.
  `SearchCut` gains the refusal (account; and the object's price for `"answer"`).
- **cli/query.{hpp,cpp}**:
  - `query_fasta` and `batched_query_fasta` take `graph::SearchBudget *` (null on the CLI), with the unit account
    on the task's stack;
  - `execute_query` catches `SearchUnitRefused` and returns the stub (the cut in scope);
  - `query_sequence` makes the handover charge and attaches the alignment;
  - `cut_to_json` writes the refused entry (2.2) from the computed cut's member set;
  - query.cpp then includes search_answer.hpp for `search_record_rows_bytes`.
- **cli/server_request.{hpp,cpp}**: `partial_refusal` read into `SearchRequest`. **cli/server_cache.cpp**:
  `make_search_cache_key` gains it.
- **cli/search_answer.{hpp,cpp}**:
  - W2's `charge_search_record` is split into measure (scratch charged), the unit object check, then charge (tree
    and transient) or the stub;
  - `search_answer_too_large_json(max_mb, account)`;
  - `search_unit_refusal_json`;
  - `search_capabilities_json(max_mb, pair_searches, records, handler_threads)`.
- **cli/server.cpp**:
  - `process_search_request` per 3.2 (no `search_results` vector; the callback builds the object under the pair's
    mutex; sort by id; append);
  - `search_budget_refusal` writes `account`;
  - the log lines (below);
  - capabilities: P (`graphs_pool`'s size, or `http_threads` on a single-graph server), E (`get_num_threads()` after
    `set_num_threads(parallel_each)`), H (`http_threads`).
- **cli/config/config.cpp**: the help text of `--search-max-answer-memory-mb` (per unit and per request; the
  field). No new flag.
- **docs**:
  - SPEC-multigraph-envelope.md section 3: the field; the cut-entry bullet gains `refused`; the budget bullet is
    rewritten (unit accounts, the object check, the request account, the 400 with `account` and its rule,
    precedence, cache); the sentence "/search has no per-pair entries ... no pair is named" goes;
  - section 8: the block of 2.6;
  - SPEC-labeled-traversal-core.md lines 4516 and 4575-4580 (the "nothing partial sent" pointer).
- **api/python/metagraph**: `GraphClient.search` and `GraphClientJson.search` take `partial_refusal` (default
  False). With it, `helpers.df_from_search_result` warns (`warnings.warn`, naming the record index,
  `seq_description` and pair) and puts the refused units in `df.attrs["refused"]`, because the DataFrame holds rows
  only.
- **Samples and fixture generators**: `-p`, `--threads-each`, `--http-threads` set explicitly (2.6). New samples
  start at 71 (train 6 uses 61-70). The changed samples are listed in LANDING for the service.
- **Log lines** (the suite parses them read-only; the exact formats go into LANDING):
  - answered: `[Server] Request <id>: answer budget: values at most <V> B a unit, answer <A> B, bound <M> MiB,
    <r> of <n> units refused` (today: `... values at most <V> B a pair, answer <A> B, bound <M> MiB`);
  - per refused unit: `[Server] Request <id>: unit refused: record <i> on <pair>: the <values|answer> account
    passed the bound of <M> MiB: held <H> B (<h> MiB) before a charge of <C> B (<c> MiB)`. `<pair>` is
    `pair (<name>, <graph>, <annotation>)` or `the graph`, and the record is named by its index, not its header.
    These are W2's fields (held bytes, refused charge), with sums saturating;
  - whole request: W2's line, its account `the answer account` or `the values account of ...` as the 400 names it.

## 7. Tests

Unit (C++; tests/annotation/test_search_budget.cpp as W2 places it, tests/cli/test_search_answer.cpp,
tests/cli/test_query.cpp):
1. Scope: a unit account that passes throws `SearchUnitRefused` and leaves the request open; other units and the
   answer go on. Once the request is refused, every charge throws as the state says. The "account pending" state
   lets values charges go on.
2. Interleavings, simulated single-threaded from a seed (review finding 8): K units with fixed charge sequences
   plus the answer account, interleaved by a seeded scheduler over 500 schedules, in both modes. The refused set,
   `priced_mb`, the status and `account` are identical in every schedule, including schedules in which the answer
   passes midway. Separately, a real-thread stress run (8 threads) checks the atomics.
3. Edges: a unit charged exactly B is answered and B + 1 B is refused; the same for the request's account; an
   object priced exactly B passes the object check and B + 1 B does not; a charge of UINT64_MAX saturates.
4. A refusal at each charge site returns the stub with the computed cut's member set: the two header paths
   (get_header_results at 1298, get_kmer_coordinates' header branch at 1667), `filter_and_aggregate` (1410), the
   signature column path at `discovery_fraction` 1 (1775), W2's store scratch, and the handover in a plain mode.
5. `query_fasta` with a budget: the stub keeps the alignment (`align`), and its siblings' `to_json` and
   `cut_to_json` bytes equal a run without a budget. No double refusal on a stub.
6. `cut_to_json` of a refused cut: the exact bytes, with and without the medians members.
7. The stub's price covers its tree (allocator-measured, as `ResolveAnswer.ThePriceCoversTheTree`).

Integration (integration_tests/test_search_budget.py, single-graph and multi-graph servers, small bounds; each case
in both modes unless it says otherwise):

8. Per mode (plain, abundance_sum, with_signature, query_counts, query_coords): two records, one over the unit
   bound and one under. With the field: 200, the over record's entry exactly, and the under record's object
   byte-identical to the bound-0 answer of a body holding only that record. Without it: 400 `"values"`.
9. With `align`, the first pair in pair order refused: the record-level members equal the bound-0 answer's.
10. Mix: 3 records x 3 pairs with refusals on a diagonal. Each answered unit's rows and cut entry are
    byte-identical to the bound-0 answer filtered to that unit, in pair order.
11. Units that each fit with a total answer that does not: 400 `"answer"` (both modes). Units refused with the
    remaining JSON passing the answer account: 400 `"answer"` with the field. Refused stubs alone passing it (many
    records, a tiny bound): 400 `"answer"`.
12. Without the field, a values refusal beside an answer-side refusal: always `"values"`, over 20 runs.
13. The object check: a verbose with_signature record whose object alone passes the bound gives `"answer"` with an
    exact `priced_mb`.
14. Every unit refused, with the field: 200, every entry refused.
15. The precedence change (3.3): a unit refused on one pair and an error on another gives that error with the field
    and 400 without it.
16. `top_labels` 0 beside a refused entry; one pair under two names (`graphs` naming both): both entries, the
    names in the tags.
17. Telomere-like: a (CCCTAA)n record and a unique record, query_coords, on a fixture whose labels hold many
    occurrences (Y3's polya25k shape): the repeat refused alone, the unique record answered.
18. Determinism: the same request 20 times on servers started with `-p 1|2|8` and `--threads-each 1|2|4`:
    identical bytes. Two different requests at once, repeated: each answer identical to its answer alone.
19. Cache (`--cache-size 1`): a partly refused answer is STOREd and the re-ask HITs with the same bytes; the same
    body without the field misses; a 400 is not cached.
20. `in_ram` gives the same bytes as mmap; bound 0 gives no refusals and the base binary's bytes
    (`METAGRAPH_BASE_BINARY`); at the default, the same.
21. Peak memory: a pair with many records whose values each fit stays near E x one record's values, not their sum
    (the values do not accumulate across a pair's records).
22. The capability block, both placements, with `-p`, `--threads-each`, `--http-threads` set.
23. The batched path (a pair without `.seqs`, plain): a unit refused alone inside a batch.
24. pymetagraph: the field passed; a warning and `attrs["refused"]` on a refused entry; nothing otherwise.

Comparisons: the route compare, base (W2's head) against the P1 head, on mini, r237form, fx11, periodic10k and
polya25k, plain and `--mmap`, without the field: no difference beyond the capability block. Samples via the shared
generator: 57-60 re-recorded (their bodies gain `account`; [U3]); new samples from 71: a multi-graph partial answer
(query_coords, one unit refused), a single-graph partial answer, every unit refused, and the 400 with each
`account`.

## 8. Cost, risks, order with W2

**Cost.** About 3.5 days: code 1.5 days (the opt-in mode and the pending-account state included), tests 1 day,
verification (builds, unit and integration DNA and Protein, route compare, samples, the GCC 13 -Werror Docker check,
CI) 1 day. Client-visible, so the text-change checklist goes in one commit (spec, help, api/python asserts, key-set
pins, Linux fixtures, samples).

**Risks.**
- R1. **Clients that send the field and ignore `cuts` lose rows silently.** Clients without the field keep 5b's 400
  (opt-in, Q2). pymetagraph warns.
- R2. **Values per pair in flight are up to E x B** (5b: B). Staging's model total is 16 GiB instead of 8 GiB; a
  large P x E needs the limits check (5.3).
- R3. **Object builds serialise per pair.** E records of one pair wait for each other's object build under the
  pair's mutex. This is cheap next to the search but not measured ([U2]).
- R4. **Samples and pins.** Samples 57-60, 12 and 37 change, and the service pinned them (1a3c46841). The service
  and suite agree first.
- R5. **The log formats change** (the suite parses them).
- R6. **Longer requests with the field (review finding 1).** A partly refused request runs all its other units,
  and an answer-account 400 can come after most units ran, holding the pool threads /resolve and /pattern share
  (3.4). Without the field, the early stop on values is kept. Against an unbounded run, refusals still shorten a
  request; against 5b's refusal, a request with the field is longer.
- R7. **The pre-cut decode stays unbounded per running unit** (accepted risk). P1 neither worsens nor bounds it.
- R8. **Master PR #765** (Kulkov, head 70c481b9f, not merged) rewrites the record fold in annotated_dbg.cpp:
  repeated rows are read once by record, the run buffers of records that can no longer be returned are freed, and
  past a fold budget (`kFoldValuesBudget`, 16 MiB in release builds at 70c481b9f; the coordinator's message said 64
  MiB) the fold keeps counts only and rereads the returned records' rows after the cut. When master is merged into
  staging, the Y4 charges (`charge_search` at about 1298, 1410, 1667, 1775) meet that code. P1 needs from that
  merge:
  - (a) the values built for the returned records, on the reread path too, take the unit's account and charge
    before each buffer grows;
  - (b) the values the fold holds before the cut stay uncharged: they are the pre-cut decode, now bounded by the
    fold budget. If they were charged, a refusal could come before the cut exists, and `execute_query`'s stub
    would need its member set from elsewhere;
  - (c) the top-N and its cut stay computed before the first charge.

**Order with train 6 lane W2.** W2 (gr/mg-t6-w2, 9e7537fed) touches the same lines: the rvalue callback,
`charge_scratch` with `threads`, `charge_search_record`, `note_passed(held, charge)`, and the refusal log line. W2
lands first; P1 stacks on its final head and adapts W2's semantics once (decided).

## 9. Decisions and remaining questions

Settled (DECISIONS.md):
- Q1: the departure is accepted.
- Q2: opt-in `partial_refusal`, in the cache key and in capabilities.
- Q3: with the field, every unit refused gives 200, and `complete` is false.
- Q4: the object check is kept.
- Q5: no `priced_mb` for values; for `"answer"`, the unit's own object; never on the 400; saturating sums.
- The 400 keeps `account`.
- Q6: the handler threads are stated.
- Q7: within the owner's rule (the budget is staging-only code).
- Q8: new samples from 71.
- Q9: caching is fine.
- Q10: the log formats go in LANDING.

Left for the coordinator:
- **L1.** The unit rule also applies without the field, so a request 5b refused only because a pair's records
  together passed the bound, with each record under it, is now answered. "Today's whole-request 400" is kept for
  every request with a refused unit, but not 5b's per-pair accounting. Recommended: one accounting for both modes.
  The alternative (keeping 5b's per-pair accounts when the field is absent) means two accounting paths, and its
  `account` member has no deterministic name.
- **L2.** The `account` rule of the 400 without the field (2.4: values first; the early stop on values kept; after
  an answer-side refusal the remaining searches run so that the name is a function of the request). The cheaper
  alternative is to name whatever passed first, which is not deterministic exactly when both would pass, as in the
  telomere request.

## 10. What I am not sure of

- **U1.** Staging's H = 2 comes from update.sh's server_query flags as read on 2026-11-10 (`-p 2 --threads-each 2
  --mmap --madv-random ... --cache-size 0`, no `--http-threads`). Not re-read; capabilities will state it.
- **U2.** The cost of building objects under the pair's mutex. Not measured.
- **U3.** Which of samples 57-60 keep status 400 at the recording bound of 1 MiB. Without the field they stay 400
  (with `account`), but which account each names needs a run.
- **U4.** Whether the service ever sends more than its two strands as records. The model total does not depend on
  it; the answer account's reachability does.

(The first version's U1, program order of a unit's charges, is gone: with no `priced_mb` for values, determinism no
longer needs it. Its U2, the medians members, is settled by finding 3: the cut is computed before every charge.)

## 11. Review findings (Fable xhigh of 8f1baf3ae) and how they are settled

Each finding was checked against 9726cb420, and against W2's 9e7537fed for the log line and the scratch.
1. **The telomere request still gets a whole-request 400, late.** Verified: search_budget.hpp:122-129 stops every
   pair after a 5b values refusal, and P1 with the field does not. Accepted: 3.4 and R6 corrected, and the
   service's fallback is described. The measurement was skipped: the figures are the suite's and the model's.
2. **Default-on.** Changed to opt-in `partial_refusal` (2.1, 2.4, 2.5, 2.6, 3.3).
3. **"Not computed" was wrong.** Verified: `search_top_n` writes the cut (annotated_dbg.cpp:543-545) before every
   values charge (1298, 1410, 1667, 1775; W2's scratch after 1387). Accepted: `execute_query` catches the refusal
   with the cut in scope (sequence moved at query.cpp:577), the member set comes from the cut (2.2, 3.2), and the
   first version's U2 is resolved.
4. **One refusal site.** Verified: the first version's handover charge in the callback would have hit the sticky
   account again for a stub. Accepted: the handover moves into `query_sequence`, `query_fasta` takes
   `SearchBudget *`, the unit account lives on the task's stack, and the callback has no account (3.2).
5. **`priced_mb`.** Verified: W2's `charge_scratch` saturates its charge at UINT64_MAX, so held + charge can
   overflow. Accepted: no `priced_mb` for values; for `"answer"`, the unit's own object; saturating sums (2.2, 6).
6. **Model totals.** Verified at server_request.cpp:35 (no `top_labels` cap), server_utils.cpp:1061 and 1082-1087
   (content copy; asynchronous send) and server_utils.cpp:53-112 (output appended without a reserve). Accepted:
   section 5 renamed and its gaps listed (5.1).
7. **Pins depend on `-p`, `--threads-each` (and now `--http-threads`).** Accepted: the generators fix them (2.6,
   6).
8. **Nits.** Accepted: the transient maximum includes the scratch of units refused by the object check (4, item 3);
   test 2 is a seeded single-threaded simulation plus a real-thread stress run (7).

The missing tests are added as tests 4, 9, 11, 12, 15, 16, 18 and 21.
