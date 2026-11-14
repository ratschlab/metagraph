# /search partial refusal per (record, pair): design note

2026-11-14. Lane T7-P1, phase 1: a design note, no code. Base: staging `gr/labeled-traversal` at 9726cb420 (train
5b). Owner decision 2026-11-14: option (c), a partial answer in one request; "it's fair that the service can refuse
some sequences, but it needs to be properly reported back with the search results". Reviewed by Fable (xhigh);
implementation (phase 2) starts after that review. Items I am not sure of are marked **[U1]** to **[U6]** and
collected in section 10; decisions that need someone else are **[Q1]** to **[Q10]** in section 9.

## 0. Summary

- The unit of refusal is one FASTA record on one pair. A unit whose values (its rows and per-position values) would
  pass `--search-max-answer-memory-mb` is refused alone. Its existing `cuts` entry gains
  `refused: {account, code, max_memory_mb, priced_mb}` and states `returned` 0; its rows are absent. Every other
  unit's rows and cut entry are byte-identical to today's answer. The status is 200.
- A unit whose own JSON object would pass the same bound is refused alone too, with `account: "answer"`.
- The request's answer account (all record objects together, then the text) stays per request. Passing it refuses
  the whole request 400 `search_answer_too_large` as today; the body gains `account: "answer"`.
- **Departure from the brief: no request-wide values total.** Today a record's values wait until its pair's answer
  is built. P1 builds each record's JSON object as soon as its unit completes and frees its values then. No
  completed unit's values are ever held, so a total over them has nothing to bound. Held values are bounded by
  the units that run at once: (pair searches at once) x (records at once per pair) x the bound. Section 5.4
  argues why the brief's total (charged from final unit sizes) would be worse.
- One flag, unchanged default 4096. No new flag. Capabilities state the two concurrency factors the memory bound
  depends on (`-p`, `--threads-each`).
- Every outcome is a function of the request, the index and the server's flags, never of thread order
  (section 4).
- Staging (`-p 2 --threads-each 2`, at most 2 requests at once): at most 16 GiB of values plus 2 x 4 GiB of
  answers plus 2 x 256 MiB of bodies, 24.5 GiB in all, on the 110 GiB the container has beyond its ~18 GiB at
  rest (section 5).
- Every request answered today is answered byte-identically. Requests refused today because one pair's values
  passed the bound are now answered, partly or fully.

## 1. Problem

Train 5b (9726cb420, on staging since 2026-11-14 10:36Z) bounds what one `POST /search` holds
(SPEC-multigraph-envelope.md section 3, "The answer budget"):
- one **values account per pair**: the rows of every FASTA record of the body on that pair and their
  per-position values, held until the pair's answer is built;
- one **answer account per request**: each record's JSON object, then the text with the transport's (or gzip's)
  copy and, with `--cache-size` above 0, the cache's copy.

Either account passing refuses the whole request 400. Nothing is sent for the other pairs or records. On
2026-11-14 the suite's run and the service's default shape (11 pairs, `query_coords` added automatically,
`top_labels` 1000, both strands as two records, a 300 bp slice) were refused on r237 because one eukaryotic pair's
values passed 4096 MiB. Every such slice was a low-complexity chromosome end (84-100% telomere repeat). With
query_coords, every position keeps each label's occurrences: about 1,500 x 270 x 8 B x 1,000 rows, roughly 3 GB per
record and pair. A unique 300 bp query holds about 2 MB. So the refusals are right, but they are too coarse: one
pathological (record, pair) costs the caller the other 21 units' answers.

The owner keeps refusals for pathological queries and wants them reported with the results, in one request.

## 2. Contract (exact)

### 2.1 The refused cut entry

Section 3 of the spec already gives each record's object `cuts`, one entry per pair searched, in pair order:
`{graph, graph_path, annotation_path, index_fp, returned, complete, nth_score, tied, tied_returned}`, plus
`nth_median` and `tied_same_median` on counts annotations. A refused unit's entry states:

| Member | Refused entry | Why |
| --- | --- | --- |
| `graph`, `graph_path`, `annotation_path`, `index_fp` | the pair's, as on every entry (null tags on a single-graph server) | identifies the pair; the record is the object the entry sits in |
| `returned` | `0` | no row of the unit is sent |
| `complete` | `false` | `complete: true` with `returned` 0 states "nothing matched", and a refusal cannot know that. `false` keeps the member's type, and a client that ignores `refused` reads "matches may exist that were not returned", which is the safe reading. Null would be more exact but changes the type for clients that read `cuts` today ([Q3]) |
| `nth_score` | `null` | there is no N-th row |
| `tied`, `tied_returned` | `0` | not computed (the work stopped), stated as in the N = 0 state |
| `nth_median`, `tied_same_median` | `null`, `0`, exactly where the pair's answered entries state these members | a refused entry has the pair's member set plus `refused` ([U2]) |
| `refused` | `{account, code, max_memory_mb, priced_mb}` | present only on a refused entry. An answered entry has no `refused` member (not `null`), so it stays byte-identical |

`refused`:

| Member | Value |
| --- | --- |
| `account` | `"values"`: the unit's rows and per-position values passed the bound while the record was searched on the pair. `"answer"`: the record's object for this pair was priced past the bound before it was built (section 3.1) |
| `code` | `"search_answer_too_large"`, the whole-request refusal's code, so one recogniser serves both |
| `max_memory_mb` | the flag's value, an integer |
| `priced_mb` | an integer: the account's bytes with the refused charge included, divided by 2^20 and rounded up. For `values` the unit's work stops at that charge, so the number is a lower bound of what the unit would hold. For `answer` it is the object's whole price (tree and transient). Deterministic per unit (section 4, [U1]) |

No member repeats the record's header. The record is the object the entry sits in: its `seq_description` and its
position in the answer array.

### 2.2 Example: multi-graph, two records, two pairs, one unit refused

Request `{"FASTA": ">fwd\n<300 bp>\n>rc\n<300 bp>", "graphs": ["a", "b"], "query_coords": true, "top_labels": 1000}`,
where record `fwd` on pair `b` passes the bound. The answer is shown compact. The server writes it as it writes
every /search answer: the same writer, members in byte order of their names.

```json
[
 {"cuts":[
   {"annotation_path":"/idx/a.row_diff_brwt_coord.annodbg","complete":false,"graph":"a","graph_path":"/idx/a.dbg","index_fp":"<a's index_fp>","nth_score":270,"returned":1000,"tied":41,"tied_returned":3},
   {"annotation_path":"/idx/b.row_diff_brwt_coord.annodbg","complete":false,"graph":"b","graph_path":"/idx/b.dbg","index_fp":"<b's index_fp>","nth_score":null,"refused":{"account":"values","code":"search_answer_too_large","max_memory_mb":4096,"priced_mb":4608},"returned":0,"tied":0,"tied_returned":0}],
  "results":[ <the 1000 rows of pair a, exactly as today> ],
  "seq_description":"fwd"},
 {"cuts":[ <pair a's entry, as today>, <pair b's entry, as today> ],
  "results":[ <the rows of pair a, then of pair b, as today> ],
  "seq_description":"rc"}
]
```

Byte identity, stated exactly: in every record object, each answered unit's rows (a contiguous run in `results`,
in pair order) and its cut entry are the bytes of today's answer. The record-level members (`seq_description`, and
with `align` the first pair's `sequence`, `score`, `max_score`, `cigar`, `orientation`) are today's too. A refused
unit keeps its alignment, which is computed before the search (section 3.2). Only the refused unit's rows are
missing and only its cut entry differs.

On a single-graph server the rule is the same with the one pair: each record of the body is its own unit, and a
refused record's object has `results: []` and one cut entry with null tags and `refused`.

### 2.3 The whole-request refusal

The request's answer account is the only whole-request budget refusal (section 3.1). Status 400, `Content-Type:
application/json`, no `Content-Encoding`, no `Retry-After`, compact, `error` first, then the other members in byte
order. The text is unchanged, so the service's recogniser and its displayed advice still apply. One member is added:

```
{"error":"search: the answer would hold more than --search-max-answer-memory-mb 4096 MiB: nothing partial is sent; ask for fewer top_labels, a larger discovery_fraction, fewer or shorter FASTA records, or no per-position values (with_signature, query_counts, query_coords)","account":"answer","code":"search_answer_too_large","max_memory_mb":4096}
```

`account` is always `"answer"` in this design. It is stated so that the body says what passed, and so that a later
request-wide account could be named without a new shape. There is **no `priced_mb`** here. Several pairs charge
the request's answer account at once, so the bytes held when the refused charge comes depend on the order of
charges, and a priced number would not be deterministic. Stating it would break "the same request gets the same
bytes". The brief asked for it only "if cheap", and it is not.

### 2.4 Status and what a client concludes

| Today (9726cb420) | P1 |
| --- | --- |
| 200 | 200, byte-identical. A unit's account is charged at most what its pair's was (W2's scratch: the unit's own, at most the pair's `threads` x the largest, see 3.1), and its object's price is at most the request's answer account, so no unit passes where its pair and request did not |
| 400, a pair's values passed | 200 with each unit over the bound refused; 200 complete if every unit fits; 400 (account `answer`) if the answer account passes once more units are answered |
| 400, the answer account passed | 400 with `account: "answer"`; or 200 with the units whose own object passes the bound refused, if what remains fits |
| any other refusal or error | unchanged (precedence in 3.3) |

A request whose every unit is refused is answered 200 with every entry refused. The rule has no exception, so
the service handles one shape ([Q3] asks whether a 400 is wanted there instead). The answer is an array (the
upstream shape), so a top-level summary would change its shape. There is none: a client counts the refused
entries it reads in `cuts`. A response header was considered and dropped, because the cache stores only the
body.

A refused entry never means "nothing matched". `refused` present: the caller did not get this unit's rows. It may
retry the record alone on the pair with smaller `top_labels`, a larger `discovery_fraction` or without per-position
values. The same request on the same server and index is refused again in the same units: the outcome is final.

### 2.5 Capabilities

`search` (top level on a single-graph server, `multigraph.search` on a multi-graph server, as today) becomes:

```json
"search": {
  "max_memory_mb": 4096,
  "pair_searches_at_once": 2,
  "records_at_once": 2,
  "rule": "SPEC-multigraph-envelope.md section 3, the /search answer budget",
  "too_large": {"status": 400, "code": "search_answer_too_large", "members": ["account", "code", "error", "max_memory_mb"]},
  "unit_refused": {"code": "search_answer_too_large", "accounts": ["answer", "values"], "members": ["account", "code", "max_memory_mb", "priced_mb"]}
}
```

- `unit_refused` is new, and its presence is the marker: an engine without it refuses the whole request (5b).
  `accounts` and `members` are sets, sorted.
- `too_large.members` gains `"account"`.
- `pair_searches_at_once`: how many pair searches the server runs at once, over all requests. On a multi-graph
  server this is the graphs' pool size (`-p`): every /search pair runs on that pool, and the pool is shared with
  the fan-outs of `/resolve` and `/pattern`. On a single-graph server /search runs on its HTTP handler thread, so
  it is the handler threads (`--http-threads`, else `-p`).
- `records_at_once`: how many records one pair search runs at once (`--threads-each`: query_fasta's pool, or the
  batched path's threads per batch, one batch at a time).
- Neither is part of any fingerprint, as today. `max_memory_mb` stays null when the flag is 0, and the two
  factors are stated even then.

Pins that change: the `search` block's member set wherever it is pinned (tests/cli/test_server_envelope.cpp sample
12's expected block, integration_tests/test_traverse.py's top-level set, integration_tests/test_search_budget.py,
api/python/tests/test_pattern_fixtures.py, tests/cli/test_search_answer.cpp's block bytes), the Linux capability
fixtures (api/python/tests/data/traverse/pattern/capabilities*/answer.json), contract samples 12 and 37, and the
suite's `--expect-member` lines (new: `search.unit_refused.*`, `search.pair_searches_at_once`,
`search.records_at_once`; changed: `search.too_large.members`).

## 3. Semantics

### 3.1 Accounts

`--search-max-answer-memory-mb` (MiB, default 4096, 0 = no bound) is the bound B of three kinds of charge, each
made before what it prices is built (the model of the deployed build's heap, as today, never a measurement):

1. **A unit's values account**, one per (record, pair). It holds what today's pair account is charged with, for
   that one record: its per-position values as their buffers grow, the store its column-label search gathers its
   rows in (W2's `charge_scratch`, which for a unit is one task, `threads` 1), and its rows as its result is
   handed over (`search_record_rows_bytes`). If it would pass B, the unit is refused with `account: "values"`.
   The request is not marked, and the other units go on.
2. **A unit's object check.** Before the record's object for this pair is built, its price (tree and transient,
   `search_record_price`, with W2's scratch for measuring it charged first as today) is compared with B. If the
   object alone would pass B, the unit is refused with `account: "answer"`, and only the refused stub's price is
   charged to the request (item 3). With one flag, an object over B would pass the request's account anyway, so
   this check only turns whole refusals into partial ones. It never refuses what 5b answers.
3. **The request's answer account**, one per request, as today: each unit's object, or its refused stub, before
   it is built; the per-pair transient high-water mark (objects of one pair are built one at a time, 3.2); then
   the text, the transport's or gzip's copy and the cache's copy. If it would pass B, the whole request is
   refused 400 (2.3), the account is marked, and every unit still running stops at its next charge, as today.

Not charged, as today: the rows decoded before the top-N cut (an accepted risk, owner 2026-11-14 11:35Z;
SPEC-multigraph-envelope.md section 3), the index (mmap pages, an `in_ram` load), the body (its 256 MiB cap), the
query graph of the batched path.

### 3.2 What changes in the pipeline: a record's object is built when its unit completes

Today `process_search_request` (cli/server.cpp) collects every record's `SeqSearchResult` of the pair, sorts them
by record id, then builds the objects. Each record's values are held until the pair's last record is done. P1
moves the build into the per-record callback. Once a unit completes, the callback:
1. charges the record's rows to the unit's account (on refusal, it turns the result into the refused stub: rows
   dropped, alignment kept, cut refused);
2. under the pair's mutex: makes the object check (3.1 item 2), charges the request's answer account, builds the
   object (rows tagged, cut entry) and keeps `(record id, object)`;
3. returns, which frees the record's result and its values.

After query_fasta returns, the pair's objects are sorted by record id and appended, which gives today's order and
bytes. The pair's mutex keeps "one object of a pair at a time", the assumption behind the per-pair transient
high-water mark (W2's `charge_search_record`). So the answer account's charges are those of 9726cb420 plus W2,
only in another order, and their total is the same.

A unit's refusal inside the search (a charge in annotated_dbg.cpp's value types or stores) throws a unit-scoped
exception (a new type, `SearchUnitRefused`, not derived from `SearchAnswerTooLarge`). `query_sequence` catches it,
because it holds both the sequence and the alignment, and returns the refused stub: an empty row vector of the
mode's type, the alignment, and a cut marked refused with its account and priced bytes. The stub flows through
the same callback. The request-scoped `SearchAnswerTooLarge` is never caught there.

Every unit's account is made when its record starts and destroyed when its object is built. At most
`pair_searches_at_once x records_at_once` unit accounts exist at once.

### 3.3 Precedence and errors

- The refusals that need no work come first, unchanged: 413 `body_too_large`, the body's JSON, the fields,
  `pairs_unsupported`, `graph_unknown`, `pair_not_found`, then a cached answer.
- A refused unit is not an error. A pair whose units are all refused is answered with refused entries.
- If the request's answer account passes, the answer is the 400 of 2.3, whatever else the work raised (another
  pair's error, a load refused `index_changed`), as today.
- Otherwise errors are answered as today: the first failing pair in pair order. One change from 5b: a request in
  which one pair's values pass and another pair errors (for example `Annotation does not support k-mer
  coordinate queries`) was the 5b budget 400. Now it is that pair's error, because the unit refusal is not a
  request refusal. Both are deterministic.

### 3.4 Interactions

- **Result cache.** A partly refused answer is a 200 and is deterministic for the process (its flags, including
  `--cache-size` and `--verbose-output`, are fixed), so it is cached like any 200: STORE, HIT, duplicate waiters
  served the same bytes. The whole-request 400 is not cached, as today. The cache's copy stays charged to the
  answer account when `--cache-size` is above 0 (staging runs 0). The cache key needs no change.
- **`top_labels`, `discovery_fraction`.** They set how many rows a unit keeps and which labels qualify, so they
  set its values. A refused entry carries no cut statistics. The advice is the 400's: fewer `top_labels`, a
  larger `discovery_fraction`, a shorter record, or no per-position values.
- **`in_ram`.** No change: the load is reserved (reservations), not charged. The same charges give the same
  refusals with and without `in_ram`.
- **Deadlines.** /search reads none (`not_after_ms`, `finish_by_ms` and `request_id` are /traverse, /resolve and
  /pattern fields, section 4 of the spec) and has no busy cap. The 899 s transport cap bounds it as every answer.
  A refused unit stops its work at the refused charge, so refusals shorten a request. Values grew at about 0.3
  GiB/s on r237 (suite, 2026-11-14), so a unit reaches 4 GiB in about 14 s.
- **`align`.** The alignment is computed before the search and is not charged to values until the handover. A
  refused unit keeps it (2.2).
- **The batched path** (a pair without `.seqs`, not query_coords: Logan's chunks) runs each record of a batch on
  one thread of the batch. The unit rule is the same.
- **`--verbose-output`** changes prices, and so which units are refused. It is a server flag, fixed per process.

### 3.5 The CLI

`metagraph query` has no budget today and keeps none. query_fasta's new per-record account factory defaults to
none, so the CLI's output and pymetagraph's local paths are unchanged.

## 4. Determinism

Claim: which units are refused, with which `account` and `priced_mb`, and whether the request is refused, are
functions of (the request, the index, the server's flags). They never depend on the order in which threads,
records or pairs run or charge.

1. **A unit's values account** is charged only by its own record's work on its own pair. The record's search
   runs on one thread: query_fasta gives each record one task, the batched path one OpenMP iteration, and
   annotated_dbg.cpp has no parallel region at its charge sites (checked at 9726cb420 and on W2's branch). The
   handover charge runs on the same thread, after the search. So the sequence of charges is fixed, and so is the
   first charge that would pass B. An account only grows, so whether some charge passes depends only on the
   total. The refused set would therefore stay deterministic even if a unit's charges one day came from several
   threads. `priced_mb` needs the fixed sequence ([U1]).
2. **A unit's object check** compares a price computed from the unit's result alone.
3. **The request's answer account** only grows. It receives a fixed set of charges: one object or stub per unit
   (fixed by items 1 and 2), the per-pair transient high-water marks (each the maximum over its pair's built
   objects, whatever their order), and the text charges, which come after every unit. The total is fixed, every
   order's prefix is at most that total, and so the account passes in every order or in none. Once it passes,
   the answer is the 400 of 2.3, whose bytes depend on nothing that ran.
4. **No account is shared between units.** A unit's charges never read another unit's account or the request's
   (beyond "has the request already passed", which only matters when the answer is the 400 anyway).
5. **No outcome reads memory in use**, the clock or the load: the concurrency bound (section 5) is a
   consequence of the structure, not an admission rule.
6. **Errors** (3.3): a unit refusal is not an error, other errors are as today, and the request's 400 wins as
   today.

The design deliberately has no request-wide values total. Section 5.4 shows that the brief's form of one stays
deterministic in what it refuses but not in which account the 400 would name.

## 5. Bounds and defaults

### 5.1 The bound

Let B = the flag (4096 MiB), P = `pair_searches_at_once`, E = `records_at_once`, C = /search requests at once (at
most the HTTP handler threads H; the service sends at most 2 to refseq-experimental). What /search holds, counted
by the model:
- **values**: at most P x E x B, over all requests. Each running unit holds at most B, a unit's values are freed
  once its object is built, and at most P x E units run at once (a unit waiting for its pair's mutex still holds
  its thread).
- **answers**: at most C x B, each request's answer account.
- **bodies**: at most H x 256 MiB (`--max-request-body-mb`).

Beside these come the uncharged pre-cut decode (per running unit, unbounded: accepted risk) and allocator slack.
The model is "not a measurement": Y4 measured refused requests peaking slightly above the bound.

### 5.2 Staging (refseq-experimental, r237, 128 GiB container, about 18 GiB anonymous at rest)

Flags: `-p 2 --threads-each 2 --mmap --madv-random --cache-size 0`, no `--http-threads`, so H = 2 (update.sh as
read on 2026-11-10; [U3]).

| | P1 | 9726cb420 (today, same reasoning) |
| --- | --- | --- |
| values | 2 x 2 x 4 GiB = 16 GiB | P x B = 2 x 4 GiB = 8 GiB (one pair account per pool thread) |
| answers | 2 x 4 GiB = 8 GiB | 8 GiB |
| bodies | 2 x 0.25 GiB = 0.5 GiB | 0.5 GiB |
| total | 24.5 GiB | 16.5 GiB |

Against 128 - 18 = 110 GiB this leaves about 85 GiB for the pre-cut decode and slack. Two concurrent worst-case
requests fit.

### 5.3 The suite's limits check, restated

The 5b check (4096 x 11 + 4096 + 256 MiB = 48.3 GiB a request, 96.5 GiB for two) assumed every pair of a request
runs at once, independently of the pool. The real bound has always depended on P, because the pool is
server-wide. With P1 the check becomes:

    P x E x B  +  C x (B + 256 MiB)  <=  container - anonymous at rest

At the default B with C = 2 this holds while P x E <= 25 (110 GiB). At P x E = 22, for example `-p 11
--threads-each 2`, it gives 88 + 8.5 = 96.5 GiB, the 5b check's own figure for two requests. Staging's P x E = 4
is far below. Hence P and E in the capabilities (2.5): the check reads them instead of assuming them. The flag's
default stays 4096 (owner 2026-11-14: "keep 4,096"). It now means 4096 MiB per unit, and no unit answered today
holds more (2.4).

### 5.4 Why no request-wide values total (the brief's design, considered and rejected)

The brief proposed keeping today's pipeline (a pair's records held until its answer is built) and adding a
request total T, charged with each unit's final size when it completes below the unit bound. Held values would
then be at most T + P x E x B. Its problems:

1. **T sums every unit of the request, not the units held at once**, because an account that may only grow
   cannot release a pair's units once its answer is built. A request answered today can hold up to 11 x 4 GiB
   over its pairs (each pair at most 4 GiB). So any T below 44 GiB refuses some request that is answered today
   (for example the measured 3.80 GiB pair beside a few others), and T >= 44 GiB protects nothing. The owner's
   5c amendment ("one 32 GiB total per request") was dropped for the same reason: it trades refusals.
2. **The 400 could not name its account deterministically.** T and the request's answer account are charged at
   the same time by different pairs. When both would pass, the first to pass depends on thread order, so
   `account` would vary between identical requests (or would have to be dropped).
3. **A second default** (a flag, or a derived total) to agree with the owner and the suite.
4. **It holds more**: every completed unit's values wait for their pair's last record, while P1 frees them.

P1 removes the reason for T: no completed unit's values wait. What remains to protect is the units in flight,
and those are bounded by P x E x B.

A server with a large P x E (a deployment with `-p 32`) would hold more under P1 than under 5b, which held only P
x B. That is an operator's sizing matter, which the limits check (5.3) now sees. An admission rule on memory in
use would break determinism and is not proposed.

## 6. Code plan (phase 2)

Base: W2's head (section 8), not 9726cb420.

- **graph/search_budget.hpp**: an account's scope (the request's answer, or a unit's values). A unit account that
  would pass keeps its held bytes and the refused charge, marks itself refused (sticky: every further charge
  throws) and throws `SearchUnitRefused{account, priced bytes}` without marking the budget. The request's
  passing keeps 5b's (and W2's) `note_passed`. `charge()` first checks the request's `passed` (and throws
  `SearchAnswerTooLarge`, as today). The budget counts refused units (for the log) and keeps `values_max` per
  unit. W2's `threads` parameter stays (unit accounts pass 1).
- **graph/annotated_dbg.{hpp,cpp}**: no change in what is charged. The `SearchAccount *` they receive is the
  unit's. `SearchCut` gains the refusal (`account`, priced bytes); or the result does (implementer's choice).
- **cli/query.{hpp,cpp}**: `query_fasta` and `batched_query_fasta` take a per-record account factory in place of
  the pair's account (null on the CLI). Each record's task makes its account, passes it to `query_sequence` and
  to the callback (callback `void(SeqSearchResult &&, graph::SearchAccount *unit)`, W2's rvalue kept).
  `query_sequence` catches `SearchUnitRefused` and returns the stub. `cut_to_json` writes the refused entry
  (2.1), with the medians members iff the pair's answered entries state them ([U2]).
- **cli/search_answer.{hpp,cpp}**: W2's `charge_search_record` is split into measure (scratch charged), the unit
  object check, then charge (tree and transient) or the stub. The stub's price is `search_record_price` with
  empty rows plus the `refused` member. `search_answer_too_large_json` gains `account`. New:
  `search_unit_refusal_json`. `search_capabilities_json(max_mb, pair_searches, records)`.
- **cli/server.cpp**: `process_search_request` per 3.2 (no `search_results` vector; the callback builds the
  object under the pair's mutex; sort by id; append). `search_budget_refusal`'s body gains `account`. Log lines
  (below). Capabilities: P (`graphs_pool`'s size, or `http_threads` on a single-graph server) and E
  (`get_num_threads()` after `set_num_threads(parallel_each)`).
- **cli/config/config.cpp**: the help text of `--search-max-answer-memory-mb` (per unit and per request, a unit
  refused alone). No new flag.
- **docs**: SPEC-multigraph-envelope.md section 3 (the cut-entry bullet gains `refused`; the budget bullet
  rewritten: unit accounts, the object check, the request account, the 400 with `account`, precedence, cache;
  the sentence "/search has no per-pair entries ... no pair is named" removed) and section 8 (the block of 2.5).
  SPEC-labeled-traversal-core.md lines 4516 and 4575-4580 (its "nothing partial sent" pointer).
- **api/python/metagraph/helpers.py**: `df_from_search_result` builds rows only, so a refused unit would vanish
  silently from the DataFrame. It warns (`warnings.warn`, naming record index, `seq_description` and pair) and
  puts the refused units in `df.attrs["refused"]`. The client's docstring says so.
- **Log lines** (the suite parses them read-only: exact formats go into the LANDING note):
  - answered: `[Server] Request <id>: answer budget: values at most <V> B a unit, answer <A> B, bound <M> MiB,
    <r> of <n> units refused` (today: `... values at most <V> B a pair, answer <A> B, bound <M> MiB`);
  - per refused unit: `[Server] Request <id>: unit refused: record <i> on <pair>: the <values|answer> account
    passed the bound of <M> MiB: held <H> B (<h> MiB) before a charge of <C> B (<c> MiB)`, where `<pair>` is
    `pair (<name>, <graph>, <annotation>)` or `the graph`. The record is named by its index, not its header;
  - whole request: W2's line unchanged (`the answer account`).

## 7. Tests

Unit (C++; tests/annotation/test_search_budget.cpp as W2 places it, tests/cli/test_search_answer.cpp,
tests/cli/test_query.cpp):
1. Scope: a unit account that passes throws `SearchUnitRefused`, leaves `passed()` false, and other units and the
   answer go on charging. Once the answer passes, every unit's next charge throws `SearchAnswerTooLarge`.
2. Interleavings: K threads x M unit accounts plus the answer account, each unit with a fixed charge sequence,
   run under 200 seeded schedules (random yields and orders). The refused set, every `priced_mb` and the
   request's outcome are identical in every schedule, including schedules in which the answer passes midway.
3. Edges: a unit charged exactly B is answered and B + 1 B is refused; the same for the request's answer account;
   an object priced exactly B passes the object check and B + 1 B does not.
4. query_fasta with the factory: the record over the bound returns the stub (rows empty; with `align` the
   alignment kept), and its siblings' `to_json` and `cut_to_json` bytes equal a run without a budget.
5. `cut_to_json` of a refused cut: the exact bytes, with and without the medians members.
6. The stub's price covers its tree (allocator-measured, as `ResolveAnswer.ThePriceCoversTheTree`).

Integration (integration_tests/test_search_budget.py, single-graph and multi-graph servers, small bounds):

7. Per mode (plain, abundance_sum, with_signature, query_counts, query_coords; one with `align`): two records,
   one over the unit bound and one under. Expect 200, the over record's entry exactly, and the under record's
   object byte-identical to the bound-0 answer of a body holding only that record.
8. Mix: 3 records x 3 pairs with refusals on a diagonal. Each answered unit's rows and cut entry are
   byte-identical to the bound-0 answer filtered to that unit, and the rows stay in pair order.
9. Units that each fit, with a total answer that does not: 400 with the exact body of 2.3.
10. The object check: a verbose with_signature record whose object alone passes the bound gives
    `account: "answer"` and an exact `priced_mb`.
11. Every unit refused: 200, every entry refused.
12. Telomere-like: a (CCCTAA)n record and a unique record, query_coords, on a fixture whose labels hold many
    occurrences (Y3's polya25k shape). The repeat is refused alone; the unique record is answered.
13. Determinism: the same request 20 times, on servers started with `-p 1|2|8` and `--threads-each 1|2|4`, gives
    identical bytes.
14. Cache (`--cache-size 1`): a partly refused answer is STOREd, and the re-ask is a HIT with the same bytes. The
    400 is not cached.
15. `in_ram` gives the same bytes as mmap.
16. Bound 0 gives no refusals and bytes identical to the base binary (`METAGRAPH_BASE_BINARY`); at the default,
    the same.
17. The capability block, both placements, with values from `-p`, `--threads-each`, `--http-threads`.
18. The batched path (a pair without `.seqs`, plain): a unit refused alone inside a batch.
19. pymetagraph: a warning and `attrs["refused"]` on a refused entry; nothing otherwise.

Comparisons: the route compare, base (W2's head) against the P1 head, on mini, r237form, fx11, periodic10k and
polya25k, plain and `--mmap`. Expect no difference beyond the capability block. Samples via the shared generator:
57-60 re-recorded ([U5]: at 1 MiB some may become 200 partial answers); new samples for a multi-graph partial
answer (query_coords, one unit refused), a single-graph partial answer and the 400 with `account`. A memory demo
(Y4's mem/ shape): the peak of a partly refused request against 5b's 400.

## 8. Cost, risks, order with W2

**Cost.** About 3 days: code 1 to 1.5 days, tests 1 day, verification (builds, unit and integration DNA and
Protein, route compare, samples, the GCC 13 -Werror Docker check, CI) 1 day. Client-visible, so it needs the
text-change checklist in one commit (spec, help, api/python asserts, key-set pins, Linux fixtures, samples) and
the service's and suite's agreement on this contract first.

**Risks.**
- R1. **Unaware clients silently lose rows.** A client that reads `results` and not `cuts` sees fewer rows, or
  `results: []`, where 5b gave a 400: pymetagraph's DataFrame client today, curl users, any service path that
  drops `cuts`. Mitigations: the capabilities marker; pymetagraph warns (section 6); the spec says a refused
  entry never means "nothing matched". The alternative is an opt-in field ([Q2]).
- R2. **Values per pair in flight are up to E x B** (5b: B). Staging holds 16 GiB server-wide instead of 8 GiB.
  Large P x E needs the limits check (5.3).
- R3. **Object builds serialise per pair.** E records of one pair wait for each other's object build under the
  pair's mutex. This is cheap next to the search, but not measured ([U4]): the route compare's timing and the
  memory demo will show it.
- R4. **`priced_mb` determinism rests on a unit's charges running in program order on one thread** ([U1]). A
  future parallel decode inside one record would make `priced_mb` vary (the refused set would not). Test 13 pins
  it.
- R5. **The medians members of a refused entry** must be computed without the work ([U2]).
- R6. **Samples 57-60, 12 and 37 change**, and the service pinned them (1a3c46841). Coordinate with the service
  and suite.
- R7. **The log formats change** (the suite parses them).
- R8. **Master.** The budget exists only on gr/labeled-traversal. When master (with Kulkov's #764 batched pre-cut
  decode, which builds values elsewhere) is merged into it, the unit account must reach those value types, as
  Y3's merge did ([Q7]).
- R9. **The pre-cut decode stays unbounded per running unit** (accepted risk). P1 neither worsens nor bounds it.

**Order with train 6 lane W2.** W2 (gr/mg-t6-w2, 9e7537fed, four commits on 9726cb420) touches the same lines:
the rvalue callback, `charge_scratch` with `threads`, `charge_search_record`, `note_passed(held, charge)`, and
the refusal log line. P1 rewrites those call sites (the pair account becomes unit accounts, and
`charge_search_record` is split). **Recommended: W2 first, then P1 stacked on W2's head** (or on the train-6
merge that contains it). P1 then keeps W2's semantics and adapts them once. Landing P1 first would force W2 to
be redone on restructured code. If both reach staging in one deploy, the suite meets one log-format change, not
two.

## 9. Open questions

- **Q1** (reviewer, coordinator): accept the departure (3.2, 5.4), that is, objects built at unit completion and
  no request-wide values total, in place of the brief's total from final unit sizes?
- **Q2** (service, owner): partial refusal on by default (recommended: the service is the client that asked for
  it, and pymetagraph warns) or opt-in with a request field such as `"partial": true`, so that unaware clients
  keep 5b's 400? /search does not refuse unknown fields (only `pairs`), so a 5b engine would ignore the field and
  either side could deploy first. The field would have to join the result cache's key.
- **Q3** (service): every unit refused gives 200 with every entry refused (recommended, uniform) or 400? And
  `complete: false` (recommended) or `null` on a refused entry?
- **Q4** (service): keep the unit object check (`account: "answer"` on a unit; recommended) or drop it, so that a
  unit is refused only for `values`?
- **Q5** (service): `priced_mb` on units (recommended; deterministic) and never on the whole-request 400 (not
  deterministic)?
- **Q6** (suite): are `pair_searches_at_once` and `records_at_once` enough for the limits check, or should the
  block also state the HTTP handler threads (the bound on requests at once, C) or the product P x E x B?
- **Q7** (coordinator, owner): the owner's 2026-11-14 11:45Z rule says search-path changes go to master first as
  a PR for Kulkov. The budget exists only on gr/labeled-traversal, so P1 can only land there. Is that within the
  rule?
- **Q8** (service, suite): samples 57-60 (pinned) change as [U5] says, and new sample numbers are needed (61 and
  up).
- **Q9** (service): is caching a partly refused answer acceptable (recommended; staging runs `--cache-size 0`
  anyway)?
- **Q10** (suite): the log formats of section 6.

## 10. What I am not sure of

- **U1.** That a unit's charges run in program order on one thread. Checked by reading only: query_fasta gives a
  record one task; the batched path one OpenMP iteration; no parallel region at annotated_dbg.cpp's charge sites,
  at 9726cb420 and on W2's branch; the handover is on the same thread. Not checked: whether any annotation
  decoder called inside a record's search spawns threads that charge (none charges today). Test 13 pins the
  result.
- **U2.** Which refused entries state `nth_median` and `tied_same_median`. `SearchCut::abundances` comes from
  whether `search_top_n` received a medians function. The header path (r237's query_coords with `.seqs`) passes
  none, and the other paths pass `medians_of_candidates`, which is empty without abundances (or without an
  abundance source). So it is a function of (pair, mode, path), not of `has_kmer_abundances()` alone. Phase 2
  computes it per pair and mode, and a test pins that a refused entry's member set equals an answered entry's of
  the same pair.
- **U3.** Staging's H = 2: update.sh's server_query flags as read on 2026-11-10 (`-p 2 --threads-each 2 --mmap
  --madv-random ... --cache-size 0`, no `--http-threads`). Not re-read today. The bound's arithmetic uses P and
  E, which are certain; H only bounds requests at once, which the service caps at 2 anyway.
- **U4.** The cost of building objects under the pair's mutex. Not measured.
- **U5.** Which of samples 57-60 become 200 partial answers at the recording bound of 1 MiB. Under P1 that is
  fixed (unit values, then the object check, then the request), but the numbers need a run.
- **U6.** Whether the service ever sends more than its two strands as records. The bound does not depend on it
  (P x E x B holds for any number of records), but the answer account's reachability does.
