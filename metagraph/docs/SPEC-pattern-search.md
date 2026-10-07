# Spec: pattern search — `POST /pattern`, contract version 1

**Status:** frozen for milestone 1 (increments 0–2 of `DESIGN-pattern-search.md`), 2026-10-07. Written from the
code as committed on `gr/labeled-traversal` at `b570800d` (route) and `63acdd7b` (engine); milestone 1b (the edge
mask: `transform --mask-dummy`, `--pattern-build-mask`, `mask: built_at_load`, the `mask_required` message naming
both remedies) added no field. **Increment 3** (2026-10-07, after `bd44e597`; `src/cli/pattern_retrieval.cpp`)
adds the projection `output.labels: "all"` to contract version 1 — additions only (§1): its request fields, answer
fields and values are in the tables below and described in §14; every answer to a request that does not ask for
it is byte for byte the milestone-1b build's; §15 lists what changed since `bd44e597`. Checked against the
fixture bodies of §11, which this build answered.
**Scope:** the server route `POST /pattern`, the `pattern` block of `GET /capabilities` and
`GET /traverse/capabilities`, and the CLI `metagraph pattern` (same answers).
**Normative design:** `DESIGN-pattern-search.md` (v6 + §22). This spec states what milestone 1 serves of it, field
by field; §13 lists where the built contract differs from the design's draft of §7 and why.
**Reader:** the search service (`PROMPT-search-service-pattern.md`), which serves this route as a job, and any
other client. Source references are to this checkout (paths relative to `metagraph/`).

---

## 1. What version 1 is

- **Frozen.** `pattern_contract_version: 1` names the request and answer of this document. A field listed here
  keeps its name, type and meaning for every server that states version 1.
- **Additions do not raise the version.** A later increment may add request fields, answer fields, values of an
  enumeration, error codes and capability fields. A field that changes meaning or type, or a value that is
  withdrawn, raises the version (`pattern_search.hpp`, `kPatternContractVersion`).
- **What a client does with that** (the rules the service's probe follows, `PROMPT-search-service-pattern.md`
  §3.1 item 6):
  - accept version 1 or higher; read fields by presence; refuse only a lower or a missing version;
  - treat an unknown value of `withheld.reason`, `cut.reason`, `stop.phase`, `stop.reason`, a slot's
    `error.code`, a refusal's `code`, a note, `placement`, `support`, `annotation`, `mask` or
    `unavailable_reason` as "not understood": pass it through, claim nothing from it;
  - gate every option on the capabilities block (§10), never on a milestone number;
  - send `output.labels` explicitly. Its default is stated in the capabilities (`default_projection`) and is
    `"none"` in version 1, also now that `"all"` is served (increment 3): a request that names no projection is
    answered as before, without reading any annotation. The design's default (`"all"`, design §7.1) would change
    what an omitted field means, so it waits for a version that says so.
- **Messages are prose.** `error` texts and slot `message`s are for people and may change; clients act on
  `code`.

What this build serves (milestone 1, and increment 3 where marked), against the design:

| design | served |
|---|---|
| modes `count`, `all_or_count`, `partial` (§5.2) | all three |
| projections `none`, `all`, `predicate_only` (§4.3, §5.6) | `none`; `all` (increment 3, §14); `predicate_only` is 400 `later_increment` |
| kinds `dna`, `iupac`, `protein` (§3) | `dna`, `iupac`; `protein` is 400 `later_increment` |
| scopes `suffix`, `any_offset`, `long` (§3) | all; `long` counts anchors only, extracts nothing (§7.7) |
| graph modes BASIC, native CANONICAL, wrapped PRIMARY (§4.1) | all three; `suffix` refused per pattern on PRIMARY |
| labels, placement, occurrences (§4.3) | read only with `output.labels: "all"` in a retrieval mode (increment 3, §14): labels on every index, placement on BASIC indexes with coordinates; otherwise none read and their counts `unknown` |
| per-label `support` for paths, `require_support` (§4.3) | later (milestone 4); a context of L ≤ k has `support: "kmer"` |
| multi-graph servers (§8) | 400 `later_increment`; the block says `multi_graph_later_increment` |
| the deadline with a finalisation reserve, 503 `deadline` (§5.3) | as designed; the annotation reads are work and stop at the work time (§14.4) |

## 2. Terms

| term | meaning |
|---|---|
| pattern | L positions, each a set of bases from {A, C, G, T}: one base (`dna`) or an IUPAC code's set (`iupac`). A pattern's N is {A, C, G, T}; it never matches a record's N symbol (design §3). |
| oriented pattern | P as given (`forward`) or its reverse complement rc(P) (`reverse`); a palindrome (P = rc(P), position by position) is one oriented pattern (`palindromic`). |
| k | the graph's k-mer length (`index.k`; 31 on refseq33m). |
| graph context | (orientation, k-mer, offset): one distinct k-mer of the index that contains the oriented pattern at that 0-based offset (L ≤ k). §7.1. |
| anchor | for L > k: a k-mer that instantiates positions [0, k) of an oriented pattern (design §4.2). Not a context. |
| retained island | a maximal run of consecutive k-mer starts of a record whose k-mers the index kept (design §3, "Covered sequence"). Every completeness statement is over retained islands. |
| information bits | Σ log2(4 / \|set_i\|) over the positions: 2 per exact base, 1 per two-base code, log2(4/3) ≈ 0.415 per three-base code, 0 per N. |
| node | the graph's id of a context's k-mer: the BOSS edge index on a DBGSuccinct, the `CanonicalDBG` wrapper id on a wrapped PRIMARY graph (§7.10). |
| row | the annotation row of the context's k-mer, named without reading it (§7.10). |
| step | one range evaluation of the discovery, or one edge examined by a mask scan (design §4.1); the unit of `max_steps`. |

## 3. Transport

- `POST /pattern` with a JSON body; the answer is compact JSON, gzip or deflate when the request's
  `Accept-Encoding` asks for it (zlib level `--traverse-compression-level`, 1 by default, so that compression
  fits the finalisation reserve). Error bodies are never compressed.
- The route runs on the server's one request pool (`-p` / `--threads-each`), shared with `/search`, `/align`,
  `/resolve` and `/traverse` (design §9). The server's content timeout is 900 s; the route's own deadline is
  `time_budget_ms` (§7.6), at most 600 s by default.
- Single-graph servers only (`server_query -i GRAPH -a ANNOTATION`). A multi-graph server (`server_query
  GRAPHS.csv`) answers every `/pattern` request with 400 `later_increment` (§6).
- `metagraph pattern -i GRAPH -a ANNOTATION [--json] REQUEST.json ...` answers each request file as the server
  would, under the same `--pattern-*` flags: the answer on stdout; a refusal's body on stdout and exit status 1.
  Answers equal the server's except `timing`.

## 4. The request

### 4.1 Top level

<!-- schema: request -->
| field | type | default | rule |
|---|---|---|---|
| `patterns` | list of pattern objects (§4.2) | required | 1 to `caps.max_patterns` (16) entries; more is 400 `invalid_request`, never cut |
| `mode` | `"count"` \| `"all_or_count"` \| `"partial"` | `"all_or_count"` | §7.5 |
| `scope` | `"suffix"` \| `"any_offset"` | `"any_offset"` | a pattern longer than k is searched as `long` whatever is named (§7.7) |
| `strands` | `"both"` \| `"forward"` \| `"reverse"` | `"both"` | `forward` searches P, `reverse` rc(P); a palindrome is searched once whatever is named (§7.3) |
| `stop_at_threshold` | boolean | `false` | stop a pattern's discovery once its count passes its threshold (§7.5) |
| `max_contexts` | integer ≥ 0 | `caps.max_contexts` (10,000) | per pattern; above the cap: lowered and listed in `limits.clamped` |
| `max_anchors` | integer ≥ 0 | `caps.max_anchors` (1,000) | per pattern, L > k: the `stop_at_threshold` threshold; lowered like `max_contexts` |
| `max_steps` | integer ≥ 1 | `caps.max_steps` (10⁸) | per **request**: the patterns spend one budget in request order (§7.6); lowered like `max_contexts` |
| `time_budget_ms` | number > `finalize_reserve_ms` | `default_time_budget_ms` (60,000) | the request's deadline (§7.6); above `caps.time_budget_ms` (600,000): lowered to it and listed |
| `output` | object (§4.3) | `{"labels": default_projection}` | what a retrieval returns |
| `max_labels_per_anchor` | integer ≥ 1 | `caps.max_labels_per_anchor` (64) | increment 3: the labels kept per row (§14.2); lowered like `max_contexts` |
| `max_annotation_work` | integer ≥ 1 | `caps.max_annotation_work` (10⁸) | increment 3: the annotation work of the request, in the oracle's units (§14.4); lowered like `max_contexts` |
| `max_memory_mb` | integer ≥ 1 | `caps.max_memory_mb` (256) | increment 3: the request's memory account (§14.4); lowered like `max_contexts` |
| `max_labels` | integer ≥ 0 | `caps.max_labels` (1,000) | increment 3, `partial` only: the labels listed per pattern (§14.5); lowered like `max_contexts` |
| `max_occurrences_per_label` | integer ≥ 0 | `caps.max_occurrences_per_label` (16) | increment 3, `partial` only: the placed occurrences listed per label (§14.5); lowered like `max_contexts` |
| `allow_unbudgeted_annotation` | boolean | `false` | increment 3: read an annotation without the budget-aware decode (§14.4) |

- An integer is a JSON number with an integral value (`5.0` is 5). A negative or fractional value is 400
  `invalid_request`.
- A field this version does not know is 400 `invalid_request` ("unknown field"), never ignored. Fields of
  later increments are refused by name with 400 `later_increment`, whatever their value, `null` included
  (§4.4); `in_ram` with 400 `resident_only`.
- The information floor is server policy (`--pattern-min-information-bits`), not a request field (§7.8).
- The six fields of increment 3 are accepted with every mode and projection; they bound only the reads of
  `output.labels: "all"` in a retrieval mode. An answer that reads no annotation (mode `count`, or
  `output.labels: "none"`) echoes none of them in `limits` (it is the milestone-1 answer) and, when the request
  named one of them or `output.labels: "all"`, carries the note `annotation_not_read` in each answered entry.

### 4.2 A pattern

<!-- schema: request_pattern -->
| field | type | rule |
|---|---|---|
| `id` | string | optional; echoed in the pattern's entry (`null` when absent). Another type is 400 `invalid_request` |
| `dna` | string over A, C, G, T (any case) | exactly one of `dna`, `iupac` |
| `iupac` | string over the 15 IUPAC codes A C G T R Y S W K M B D H V N (any case) | exactly one of `dna`, `iupac` |

- Neither or both of `dna` and `iupac`, or a value that is not a string, is 400 `invalid_request`.
- A string with another character (U, `-`, `.`, whitespace; an IUPAC code in a `dna` pattern) or an empty
  string is answered in the pattern's slot with `bad_alphabet` (§8.9); the other patterns are answered.
- No length cap: a pattern longer than k costs only its anchor window (§7.7).
- `protein` is 400 `later_increment` (design §6).

### 4.3 `output`

<!-- schema: request_output -->
| field | type | rule |
|---|---|---|
| `labels` | `"none"` \| `"all"` | `"none"`: the label-free projection (design §4.3), contexts without any annotation read. `"all"` (increment 3): every returned context with its labels, placed where the index can (§14); in mode `count` it reads nothing and is said so (note `annotation_not_read`). `"predicate_only"` is 400 `later_increment` in every mode, `count` included; any other string is 400 `invalid_request` |
| `occurrences` | boolean | with `labels: "all"`: whether the labels are placed (default `true`, the capabilities' `default_occurrences`; `false` answers `placement: "not_requested"`). With `labels: "none"` (or omitted): `false` is accepted and changes nothing, `true` is 400 `invalid_request` (occurrences are placed per label) |
| `paths` | boolean | `false` is accepted and changes nothing; `true` is 400 `later_increment` |

### 4.4 Fields of later increments

Refused by name (400 `later_increment`), whatever their value: `max_paths`, `require_support`, `predicate`,
`max_predicate_contexts`, `max_predicate_work`, `graphs`, `genetic_code`, `budget_split`;
`patterns[i].protein`; `output.labels: "predicate_only"`; `output.paths: true`. So that a request written for
a later increment is told what to wait for, not answered as if the field were absent. (Until increment 3 the
list also held `max_labels_per_anchor`, `max_annotation_work`, `max_memory_mb`, `max_labels`,
`max_occurrences_per_label`, `allow_unbudgeted_annotation`, `output.labels: "all"` and
`output.occurrences: true`; they are served now, §4.1, §4.3.)

### 4.5 Caps and defaults (server flags)

| flag | default | request field | above it |
|---|---|---|---|
| `--pattern-max-contexts` | 10,000 | `max_contexts` (default = cap) | lowered, listed in `limits.clamped` |
| `--pattern-max-anchors` | 1,000 | `max_anchors` (default = cap) | lowered, listed |
| `--pattern-max-steps` | 100,000,000 | `max_steps` (default = cap) | lowered, listed |
| `--pattern-default-time-ms` | 60,000 | `time_budget_ms` when omitted | — |
| `--pattern-max-time-ms` | 600,000 | `time_budget_ms` | lowered to the cap (not to the default), listed |
| `--pattern-finalize-ms` | 250 | — (the reserve inside every budget, §7.6) | — |
| `--pattern-min-information-bits` | 24 | — (the floor, §7.8) | — |
| `--pattern-max-patterns` | 16 | length of `patterns` | refused (400), never cut |
| `--pattern-max-labels-per-anchor` | 64 | `max_labels_per_anchor` (default = cap) | lowered, listed |
| `--pattern-max-annotation-work` | 100,000,000 | `max_annotation_work` (default = cap) | lowered, listed |
| `--pattern-max-memory-mb` | 256 | `max_memory_mb` (default = cap) | lowered, listed |
| `--pattern-max-labels` | 1,000 | `max_labels` (default = cap) | lowered, listed |
| `--pattern-max-occurrences` | 16 | `max_occurrences_per_label` (default = cap) | lowered, listed |
| `--traverse-chunk-target-ms` | 50 | — (not a cap: the annotation reads under the deadline are decoded in chunks of about this duration, as `/traverse`'s) | — |

The capabilities state every value in force (`caps`, `default_time_budget_ms`, `finalize_reserve_ms`), and every
answer echoes the effective ones (`limits`, §8.3). A clamp is never silent: `limits.clamped` lists it. The
service refuses values above its own ceilings at submit (`PROMPT-search-service-pattern.md` §3.1 item 5); the
server lowers them, so that a direct caller is not refused for a budget the host can give in part.

## 5. Order of the checks

A request is refused by the first check it fails, in this order:

1. multi-graph server: 400 `later_increment`, whatever the body;
2. the single index still loading: 503, no `code` (§6);
3. the body is not JSON: 400 `invalid_request`;
4. the graph (`PatternSearch::support`): 400 with the graph's reason (`mask_required`, …), whatever the body asks;
5. the body is not an object: 400 `invalid_request`;
6. a later-increment field (§4.4, top level) or `in_ram`, in the alphabetical order of the body's field names;
7. `patterns` (presence, list, length), then each pattern in order (`protein`, `id`, exactly one of `dna` /
   `iupac`, its type, an unknown field); a pattern's alphabet is not a refusal (§8.9);
8. `mode`, `output` (`labels`, `occurrences`, `paths`, an unknown field), `scope`, `strands`,
   `stop_at_threshold`, `max_contexts`, `max_anchors`, `max_steps`, `time_budget_ms`, then increment 3's
   `max_labels_per_anchor`, `max_annotation_work`, `max_memory_mb`, `max_labels`, `max_occurrences_per_label`,
   `allow_unbudgeted_annotation`;
9. an unknown top-level field;
10. increment 3: `output.labels: "all"` in a retrieval mode on an annotation without the budget-aware decode,
    without `allow_unbudgeted_annotation: true`: 400 `annotation_unbudgeted`.

## 6. Whole-request refusals

The body is `{"error": <message>, "code": <code>}`, except the 503 during loading, which has no `code`.

<!-- schema: refusal -->
| field | type | meaning |
|---|---|---|
| `error` | string | for people; may change |
| `code` | string | what to act on (table below) |

| status | `code` | when | what a client does |
|---|---|---|---|
| 400 | `invalid_request` | not JSON, not an object, a wrong type or value, an unknown field, an empty or too long `patterns` list, `max_steps` < 1, `time_budget_ms` ≤ the reserve, `max_labels_per_anchor`, `max_annotation_work` or `max_memory_mb` < 1, `output.occurrences: true` without `output.labels: "all"` | fix the request |
| 400 | `later_increment` | a field or value of §4.4; a multi-graph server | wait for the increment the capabilities will announce |
| 400 | `resident_only` | `in_ram`, any value: the route never loads an index inside a request (design §5.3) | drop `in_ram` |
| 400 | `mask_required` | the graph was loaded without its dummy-edge mask (`.edgemask`); without it every dummy edge would count as a k-mer | the host's operator creates the mask (`metagraph transform --mask-dummy`, or `--pattern-build-mask` at start-up; design §4); the capabilities say `mask: absent` meanwhile |
| 400 | `representation_unsupported` | not a succinct graph (nor a PRIMARY one wrapped in `CanonicalDBG`), or k < 2 | none: this host has no pattern search |
| 400 | `primary_unwrapped` | a PRIMARY graph not wrapped in `CanonicalDBG` (the server always wraps; CLI or embedding misuse) | none |
| 400 | `alphabet_unsupported` | the graph's alphabet is neither `$ACGT` nor `$ACGTN` | none |
| 400 | `annotation_unbudgeted` | increment 3: `output.labels: "all"` in a retrieval mode, and the annotation has no budget-aware decode (capabilities `annotation: "unbudgeted"`: a column, BRWT, row or disk annotation; only the row-diff family has one) | set `allow_unbudgeted_annotation: true` to read it without a memory bound on the reads (§14.4), or ask for `labels: "none"` or mode `count` |
| 503 | `deadline` | the answer could not be written by `time_budget_ms` (§7.6): nothing partial is sent | retry with a larger budget, or narrow the request |
| 503 | (none) | the single index is still loading (`Retry-After: 60`); every route answers so | retry later |

An unexpected failure (a server bug) is answered as every route answers it: 400 `{"error": …}` without a
`code` (500 `Internal server error` for a non-standard exception).

## 7. Semantics

### 7.1 Graph contexts

- A **graph context** is (orientation, k-mer, offset): one distinct k-mer of the index that contains the
  oriented pattern at that offset. Two contexts are the same iff all three agree (design §3).
- It is **not** a record, an occurrence or a hit. A k-mer present in ten records, or in ten samples, is one
  context. Repeated copies of a k-mer in one record share one context.
- A k-mer that contains the pattern twice (`ACGAC` holds `AC` at offsets 0 and 3) gives two contexts.
- An IUPAC pattern and its reverse complement can both match one instance (`NA` and `TN` both match `TA`): two
  contexts at one (k-mer, offset) that differ in orientation. An exact DNA pattern cannot.
- A context's k-mer is a valid seed for `/traverse` on the same graph (design §7.2).
- Counts of contexts are not counts of occurrences; `counts.occurrences` is `unknown` unless `output.labels:
  "all"` placed them on a record-placement index (increment 3, §14.3).

### 7.2 Scopes and `absence_scope`

| scope | requested by | counts and returns | `absence_scope` |
|---|---|---|---|
| `any_offset` | default, L ≤ k | every k-mer containing the oriented pattern at any offset p ∈ [0, k − L]: every occurrence inside a retained k-mer | `any_offset` |
| `suffix` | `scope: "suffix"`, L ≤ k | the k-mers whose last L symbols instantiate it (offset k − L only): every occurrence starting at island position ≥ k − L | `suffix_only` |
| `long` | implied by L > k | anchors (§7.7) | `long` |

- `suffix` is the cheap, exactly countable scope; `any_offset` discovers the flank ranges, branching work even for
  an exact pattern, charged to `max_steps` (design §4.1).
- `suffix` is not served on a wrapped PRIMARY graph (a virtual suffix is a stored prefix, design §4.1): such a
  pattern is answered in its slot with `scope_unsupported`; `any_offset` is complete there.
- An exact 0 count says that no k-mer of this index contains the oriented pattern(s) in that scope, over its
  retained islands. It says nothing about a label, a sample or a record as deposited (design §3, §5.1).

### 7.3 Strands, orientations and palindromes

- `strands` selects the oriented patterns: `both` searches P and rc(P), `forward` P, `reverse` rc(P). Each is
  searched completely on its own (design §4.1).
- A palindromic pattern is searched once, whatever `strands` says. Its contexts count once in every total.
- On a BASIC graph (`index.strand_stated: true`) orientations are strands: a context is on `"+"` (P), `"-"`
  (rc(P)) or `"="` (palindromic); counts are split by `by_strand` with keys `"+"`, `"-"`, `"both"` (the
  palindromic contexts, never under `+` and `-` separately).
- On a native CANONICAL graph and a wrapped PRIMARY graph a stored k-mer may be the reverse complement of the
  deposited one, so no strand is known: contexts carry `orientation` (`"forward"`, `"reverse"`, `"palindromic"`),
  counts are split by `by_orientation` with those keys, and every entry carries the note
  `strand_unknown_canonical`. Both orientations of every k-mer are in such a graph, so `forward` and `reverse`
  count alike in `any_offset` scope.
- A context's `instance` is the matched bases of the k-mer: an instance of P for `+`, `=`, `forward` and
  `palindromic`, of rc(P) for `-` and `reverse`.

### 7.4 Counts and relations

Every count is `{value, relation, unit}` (design §3):

| relation | meaning | `value` |
|---|---|---|
| `exact` | the discovery behind it completed: no step, time or threshold stop touched it | the count |
| `at_least` | a stop interrupted discovery; an undiscovered branch, offset or strand has no upper bound | the sum of the lower bounds of what was explored |
| `bounds` | every range of every branch, offset and orientation was discovered; only scans of masked edges were interrupted; `lower` ≤ true ≤ `upper` | `lower` |
| `unknown` | the phase never ran: after a stop, or not in this increment | `null` |

- `bounds` carries `lower` and `upper`; no other relation does. `bounds` stays `bounds` when `lower` = `upper`.
- Sums (a total over offsets or strands) take the weakest relation: `unknown` + `unknown` = `unknown`;
  `unknown` + anything else = `at_least`; `at_least` + anything = `at_least`; `bounds` + `bounds` or `exact` =
  `bounds`; `exact` + `exact` = `exact` (`pattern_search.hpp`, `Count::operator+=`).
- A count is never promoted by assumption: anchors say nothing about paths, contexts nothing about occurrences.
- Units: `graph_contexts`, `anchors`, `paths`, `placed_occurrences`, `labels`. `labels` and `placed_occurrences`
  are `unknown` unless `output.labels: "all"` read them in a retrieval mode (§14.6).

### 7.5 Modes, `withheld` and `cut`

| mode | results | `retrieval_complete` |
|---|---|---|
| `count` | none: the entry has no `withheld`, `returned`, `cut` or `results` | always `false` (a count returns no context) |
| `all_or_count` | every context, only when discovery completed with an `exact` total ≤ `max_contexts`; otherwise none, `withheld` says why | `true` iff every context was returned |
| `partial` | the first `max_contexts` contexts in answer order (§7.9) among those discovered, also after a step or threshold stop; `cut` says why the list is shorter than the pattern's contexts | `true` iff every context was returned |

`withheld.reason` (results absent; `returned: 0`, `results: []`):

| reason | when | the count | what to change |
|---|---|---|---|
| `count_above_threshold` | `all_or_count`: discovery completed, `exact` total > `max_contexts` | `exact` | narrow the pattern, scope or strand; or `partial` |
| `threshold_crossed` | `all_or_count`, `stop_at_threshold`, L ≤ k: discovery stopped once the count passed `max_contexts` | `at_least` | as above |
| `discovery_budget` | `all_or_count`, L ≤ k: `max_steps` reached in discovery or a mask scan | `at_least` or `bounds` | a more informative pattern or a narrower scope; a filter does not help |
| `deadline` | `all_or_count`, L ≤ k: the work time passed in discovery or a mask scan, or in the release (which is all or nothing) | as stopped | a larger `time_budget_ms`, or as for `discovery_budget` |
| `paths_later_increment` | either retrieval mode, L > k, unless the anchors are `exact` 0 (§7.7); whatever stopped the anchors is in `stop` | the anchors' | nothing yet: paths arrive with milestone 4 |
| `annotation_budget` | increment 3, `all_or_count`, `labels: "all"`: a row the memory account refused (`rows_refused`), or the reads stopped at `max_annotation_work` | `exact` | raise `max_memory_mb` or `max_annotation_work`, narrow the pattern, or `partial` |
| `anchor_labels_truncated` | increment 3, `all_or_count`, `labels: "all"`: a row carried more labels than `max_labels_per_anchor` (`anchors_truncated` lists each, with its total) | `exact` | raise `max_labels_per_anchor` to the largest total, or `partial` |
| `output_budget` | increment 3, `all_or_count`, `labels: "all"`: the memory account could not hold the answer (the contexts' descriptors or the labels and occurrences built for them) | `exact` | raise `max_memory_mb`, narrow the pattern, or `partial` |

With `labels: "all"`, `deadline` also names a time stop of the annotation reads (`stop.phase` `label_discovery`
or `placement`); the graph count is then still the discovery's.

`cut.reason` (`partial`, neither complete nor withheld; `results` holds what was released):

| reason | when |
|---|---|
| `max_contexts` | discovery completed (or stopped at the threshold) with more contexts than `max_contexts`: the first `max_contexts` in answer order |
| `max_steps` | discovery stopped at `max_steps`: the first `max_contexts` of those discovered |
| `time` | the work time passed: in discovery nothing is released (its membership would depend on the machine, `returned: 0`); in the release, what was released before |
| `max_memory` | increment 3, `labels: "all"`: the memory account held the descriptors of only the first `returned` contexts (`partial`'s descriptors take at most half of the account, §14.4; the list's length; it replaces the engine's `max_contexts` cut when both cut) |

`max_anchors` is a value of `cut.reason` reserved for a later increment's release of anchors; version 1 never
releases anchors.

### 7.6 Budget, deadline and the finalisation reserve

- **One budget per request** (design §5.3): `max_steps` and the deadline are spent by the patterns in request
  order. A stop is sticky: every later pattern answers `unknown` counts with the same `stop` (phase
  `discovery`), `work` zero, and in a retrieval mode the matching `withheld` or `cut`. `stop_at_threshold` is not
  a budget stop: it ends only its own pattern.
- **The deadline** starts when the request's body is parsed (time spent in the server's queue is not in it).
  Work stops at `time_budget_ms − finalize_reserve_ms`; the counts kept up to date during the search are then
  written in the reserve. The clock is read at every phase and pattern boundary and at least every 4,096 steps
  or released edges.
- **503 `deadline`.** The answer is assembled, serialised and compressed under the same deadline, read every
  4,096 entries and results while it is assembled, every 64 KiB of text while it is written, between
  compression blocks, and once more before it is handed to the transport. If it cannot be written by
  `time_budget_ms`, the answer is 503 `deadline` and nothing partial is sent.
- `stop` names where the work stopped: `{phase, reason}` with phase `discovery` (the range search), `mask_scan`
  (the deferred scans of masked edges, the only phase whose stop can leave `bounds`) or `extraction` (the
  release, after discovery), and reason `max_steps`, `time`, `max_contexts` or `max_anchors` (the last two:
  `stop_at_threshold`). Increment 3 (`labels: "all"`) adds the phases `label_discovery` (the first read, §14.2),
  `placement` (the second) and `output` (the descriptors and labels the account holds for the answer), and the
  reasons `max_annotation_work` and `max_memory` (§14.4); the engine's stop, when there is one, is the one
  stated.
- `work` per pattern: `ranges_visited` (range evaluations), `mask_scans` (ranges whose scan began), `steps`
  (every step charged: `ranges_visited` plus the edges scanned). The `steps` of all patterns sum to at most
  `max_steps`.

### 7.7 Patterns longer than k

- Scope `long`, whatever was requested. `anchor_information_bits` states the bits of the anchor window
  [0, k); the information floor gates on it (§7.8).
- `counts.anchors` counts the anchors (unit `anchors`, by strand or orientation); `counts.paths` is `unknown`,
  except `exact` 0 when the anchors are `exact` 0 (no path starts without an anchor: a derivation, not a
  promotion). `counts.contexts` is absent.
- Nothing is extracted: in a retrieval mode the results are withheld with `paths_later_increment`, unless the
  anchors are `exact` 0, in which case the empty answer is complete (`retrieval_complete: true`).
- The note `paths_later_increment` is on every such entry. `max_anchors` is the `stop_at_threshold` threshold.

### 7.8 The information floor

- A pattern below `caps.min_information_bits` (24 by default, about 12 specified bases) is answered in its
  slot with `information_below_floor` and costs no step. For L > k the anchor window's bits count.
- Exempt: an exact pattern (every position one base, whatever its kind) in `suffix` scope, however short: one
  range, a few ranks.
- Low-complexity patterns are not refused; an exact one that sdust flags carries the note
  `low_complexity_pattern`.

### 7.9 Ordering and determinism

- `patterns`: request order, one entry per request pattern.
- `results`: node ascending, then offset ascending, then orientation (`forward`/`+`, `reverse`/`-`,
  `palindromic`/`=`) (design §5.5, "BOSS edge order, then by offset"; the orientation separates the two
  contexts an IUPAC pattern can have at one offset).
- `strands` (entry): the orientations searched, forward before reverse; `["="]` or `["palindromic"]` for a
  palindrome.
- `notes`: `low_complexity_pattern`, `strand_unknown_canonical`, `paths_later_increment`, then increment 3's
  `annotation_unbudgeted`, `record_bounds_unknown` or `annotation_not_read`, in that order.
- `limits.clamped`: `max_contexts`, `max_anchors`, `max_steps`, `time_budget_ms`, then increment 3's
  `max_labels_per_anchor`, `max_annotation_work`, `max_memory_mb`, `max_labels`, `max_occurrences_per_label`,
  in that order.
- Increment 3: labels in `by_label` and in each result by (contexts desc, column asc) (design §5.5); a label's
  occurrences by (`seq_id`, start, strand), or (`kmer_coord`, `offset`, strand) without a record mapping.
- JSON objects (`by_offset`, `by_strand`, `by_orientation`, …) carry no order; the server writes keys sorted
  as strings (`"10"` before `"2"`).
- **Determinism.** The same request on the same index under the same caps gives the same answer, byte for byte
  apart from `timing`, including which contexts a cut kept, because every budget is spent in a fixed order. The
  exception is a time stop: the entry it touched and every answered entry after it state
  `determinism: "time_limited"`; their counts, `work` and released contexts depend on the machine. A 503
  `deadline` depends on the machine too.

### 7.10 Node and row ids

- `node` is the id of the context's k-mer in the served graph: the BOSS edge index on a DBGSuccinct (BASIC,
  CANONICAL), the wrapper id on a wrapped PRIMARY graph (a stored k-mer's id, or that plus the wrapper's offset
  for a virtual reverse complement).
- `row` is the annotation row the k-mer is annotated in, named without reading it: `node − 1` on BASIC; on
  PRIMARY the stored k-mer's (`get_base_node(node) − 1`); on native CANONICAL the canonical k-mer's (the
  annotation key of every route), so that a k-mer and its reverse complement share one row. `null` if the k-mer
  has no row in the annotation (never expected on a compatible index; stated rather than guessed).
- Both are **opaque** to clients and valid only per (host, `index.index_fp` or `index.release`, graph). A later
  increment takes row ids back ("labels for given rows", design §12); a client keeps them as given.

## 8. The answer (HTTP 200)

### 8.1 Top level

<!-- schema: answer -->
| field | type | meaning |
|---|---|---|
| `pattern_contract_version` | integer | 1 |
| `mode` | string | the request's mode, defaults applied |
| `output` | object \| null | `null` in mode `count`; in the retrieval modes the projection applied: `{"labels": "none"}`, or `{"labels": "all", "occurrences": <boolean>}` (increment 3) |
| `index` | object (§8.2) | what was searched |
| `limits` | object (§8.3) | the effective caps of this request |
| `timing` | object (§8.4) | the request's elapsed time |
| `patterns` | list of entries (§8.5) | one per request pattern, in request order |

### 8.2 `index`

<!-- schema: index -->
| field | type | meaning |
|---|---|---|
| `index_ns` | string \| null | the server's `--index-name`, `null` without it |
| `index_fp` | string \| null | the index identity (`--index-manifest`'s digest), `null` without a manifest |
| `release` | string | the server's `--index-release` (`""` when not set) |
| `k` | integer | the graph's k |
| `graph_mode` | `"basic"` \| `"canonical"` \| `"primary"` | §7.3 |
| `alphabet` | `"$ACGT"` \| `"$ACGTN"` | the BOSS alphabet, sentinel first |
| `strand_stated` | boolean | `true` on BASIC: contexts carry `strand`, counts `by_strand` |

### 8.3 `limits`

<!-- schema: limits -->
| field | type | meaning |
|---|---|---|
| `max_contexts` | integer | effective (§4.1) |
| `max_anchors` | integer | effective |
| `max_steps` | integer | effective, for the whole request |
| `time_budget_ms` | number | effective |
| `finalize_reserve_ms` | number | the reserve inside it |
| `min_information_bits` | number | the floor |
| `max_patterns` | integer | the server's |
| `stop_at_threshold` | boolean | as requested |
| `max_labels_per_anchor` | integer | increment 3, answers with `labels: "all"` only: effective (§4.1) |
| `max_annotation_work` | integer | likewise |
| `max_memory_mb` | integer | likewise |
| `max_labels` | integer | likewise |
| `max_occurrences_per_label` | integer | likewise |
| `allow_unbudgeted_annotation` | boolean | likewise, as requested |
| `clamped` | list of objects | one per request value lowered to its cap, in the order of §7.9; `[]` when none |

<!-- schema: clamped -->
| field | type | meaning |
|---|---|---|
| `field` | string | `max_contexts`, `max_anchors`, `max_steps` or `time_budget_ms`; increment 3: `max_labels_per_anchor`, `max_annotation_work`, `max_memory_mb`, `max_labels`, `max_occurrences_per_label` |
| `requested` | number | the value asked for |
| `effective` | number | the cap applied |

### 8.4 `timing`

<!-- schema: timing -->
| field | type | meaning |
|---|---|---|
| `elapsed_ms` | number | top level: from the deadline's start to the end of the search (the annotation reads included); in an entry: that pattern's search. Varies between runs |
| `label_discovery_ms` | number | increment 3, entries with `labels: "all"`: the first read (§14.2). Varies between runs |
| `placement_ms` | number | likewise, the second read (§14.3) |

### 8.5 A pattern's entry

An entry is one of three shapes: answered, refused by the engine, or refused for its alphabet (§8.9).

<!-- schema: entry -->
| field | type | present | meaning |
|---|---|---|---|
| `id` | string \| null | always | the request's `id` |
| `kind` | `"dna"` \| `"iupac"` | always | which key the pattern came in |
| `pattern` | string | answered, engine refusal | the pattern in upper case |
| `length` | integer | answered, engine refusal | L |
| `information_bits` | number | answered, engine refusal | §2 |
| `anchor_information_bits` | number \| null | answered, engine refusal | L > k: the bits of [0, k); else `null` |
| `error` | object (§8.9) | refusals only | `{code, message}` |
| `mode` | string | answered | the request's mode |
| `scope` | `"suffix"` \| `"any_offset"` \| `"long"` | answered | §7.2 |
| `strands` | list of strings | answered | the orientations searched (§7.3, §7.9) |
| `palindromic` | boolean | answered | P = rc(P) |
| `counts` | object (§8.6) | answered | |
| `work` | object (§8.7) | answered | |
| `stop` | object \| null | answered | §7.6; `null` when the pattern's work completed |
| `retrieval_complete` | boolean | answered | §7.5; the one flag that licenses an absence claim over graph contexts (never over labels) |
| `withheld` | object \| null | answered, retrieval modes | `{reason}` (§7.5) |
| `returned` | integer | answered, retrieval modes | the length of `results` |
| `cut` | object \| null | answered, retrieval modes | `{reason}` (§7.5); only in `partial` |
| `results` | list of results (§8.8) | answered, retrieval modes | the released contexts, in the order of §7.9 |
| `absence_scope` | `"suffix_only"` \| `"any_offset"` \| `"long"` | answered | §7.2 |
| `determinism` | `"full"` \| `"time_limited"` | answered | §7.9 |
| `notes` | list of strings | answered | §8.10 |
| `timing` | object (§8.4) | answered | |
| `placement` | string | answered, `labels: "all"` | increment 3: what this answer places: `record`, `global`, `none`, `none_canonical` (§14.3), or `not_requested` (`output.occurrences: false`) |
| `annotation` | `"budgeted"` \| `"unbudgeted"` | answered, `labels: "all"` | the reads' access (§14.4) |
| `by_label` | list \| null | answered, `labels: "all"` | the per-label summary over the returned contexts (§14.5); `null` when the results are withheld |
| `rows_refused` | list | answered, `labels: "all"` | the rows the memory account refused, each once (§14.4) |
| `anchors_truncated` | list | answered, `labels: "all"` | the rows cut at `max_labels_per_anchor`, each once (§14.2) |
| `labels_cut` | object \| null | answered, `labels: "all"` | `partial`: `{reason: "max_labels", returned}` when `by_label` and the results list fewer labels than were found |
| `occurrences_cut` | object \| null | answered, `labels: "all"` | `partial`: `{reason: "max_occurrences_per_label", labels}` when that many labels list fewer occurrences than they have |

### 8.6 `counts`

<!-- schema: counts -->
| field | type | present | meaning |
|---|---|---|---|
| `contexts` | contexts count | L ≤ k | unit `graph_contexts` |
| `anchors` | anchors count | L > k | unit `anchors` (§7.7) |
| `paths` | count | L > k | unit `paths`; `unknown` (`exact` 0 without anchors) |
| `labels` | count | always | unit `labels`; `unknown` unless `labels: "all"` read them (§14.6) |
| `occurrences` | count | always | unit `placed_occurrences`; `unknown` unless `labels: "all"` placed them in records (§14.6) |

<!-- schema: count -->
| field | type | meaning |
|---|---|---|
| `value` | integer \| null | §7.4; `null` iff `unknown` |
| `relation` | `"exact"` \| `"at_least"` \| `"bounds"` \| `"unknown"` | §7.4 |
| `unit` | string | §7.4 |
| `lower` | integer | `bounds` only |
| `upper` | integer | `bounds` only |

A contexts count is a count with these fields besides (every one a count of unit `graph_contexts`):

<!-- schema: contexts_count -->
| field | type | meaning |
|---|---|---|
| `suffix` | count | the contexts at offset k − L (equal to the total in `suffix` scope) |
| `by_offset` | object: offset (decimal string) → count | one key per offset of the scope (every p in [0, k − L] for `any_offset`, `k − L` for `suffix`), zeros included: a missing key never stands for zero |
| `by_strand` | object: `"+"` \| `"-"` \| `"both"` → count | BASIC graphs: one key per orientation searched |
| `by_orientation` | object: `"forward"` \| `"reverse"` \| `"palindromic"` → count | CANONICAL and PRIMARY graphs, instead of `by_strand` |

An anchors count is a count with `by_strand` or `by_orientation` besides (unit `anchors`), as above.

<!-- schema: anchors_count -->
| field | type | meaning |
|---|---|---|
| `by_strand` | object → count | BASIC graphs |
| `by_orientation` | object → count | CANONICAL and PRIMARY graphs |

### 8.7 `work`, `stop`, `withheld`, `cut`

<!-- schema: work -->
| field | type | meaning |
|---|---|---|
| `ranges_visited` | integer | range evaluations (one step each) |
| `mask_scans` | integer | ranges whose scan of masked edges began; 0 in the common case |
| `steps` | integer | every step this pattern charged |
| `annotation_rows` | integer | increment 3, `labels: "all"`: the rows this pattern's reads returned (both steps) |
| `annotation_units` | integer | likewise: the work units they cost (§14.4) |
| `memory_bytes` | integer | likewise: the request's memory account at its peak so far (the model of §14.4) |

<!-- schema: stop -->
| field | type | meaning |
|---|---|---|
| `phase` | `"discovery"` \| `"mask_scan"` \| `"extraction"`; increment 3: `"label_discovery"` \| `"placement"` \| `"output"` | §7.6 |
| `reason` | `"max_steps"` \| `"time"` \| `"max_contexts"` \| `"max_anchors"`; increment 3: `"max_annotation_work"` \| `"max_memory"` | §7.6 |

<!-- schema: reason -->
| field | type | meaning |
|---|---|---|
| `reason` | string | `withheld`: §7.5's first table; `cut`: its second. An object, so that a later increment can add the shard (design §8) |

### 8.8 A result (a context)

<!-- schema: result -->
| field | type | meaning |
|---|---|---|
| `kmer` | string | the context's k-mer (k bases), as the graph spells it |
| `instance` | string | `kmer[offset, offset + L)`: the matched bases, which for an IUPAC pattern name the variant |
| `offset` | integer | 0-based, in [0, k − L] |
| `strand` | `"+"` \| `"-"` \| `"="` | BASIC graphs (`index.strand_stated`) |
| `orientation` | `"forward"` \| `"reverse"` \| `"palindromic"` | CANONICAL and PRIMARY graphs, instead of `strand` |
| `node` | integer | §7.10 |
| `row` | integer \| null | §7.10 |
| `support` | `"kmer"` | increment 3, `labels: "all"`: the label annotates the one k-mer (design §4.3) |
| `labels_status` | string | likewise: `complete`, `truncated`, `refused`, `not_read`, `output_budget` (§14.5) |
| `labels_total` | integer \| null | likewise: the labels the row carries (its true total, also when truncated); `null` when it was not read |
| `labels` | list \| null | likewise: the context's labels (§14.5); `null` unless the row was read and held |

With `output.labels: "none"` a result has the first seven fields only: nothing of the annotation is read.

### 8.9 Per-pattern errors (the error slot)

A refused pattern is answered in its slot, inside a 200 answer; the other patterns are answered (design §7.2).
It costs no step.

<!-- schema: error -->
| field | type | meaning |
|---|---|---|
| `code` | string | table below |
| `message` | string | for people |

| `code` | when | the slot carries |
|---|---|---|
| `bad_alphabet` | a character outside the kind's alphabet, or an empty pattern (the message names the first offending 0-based position) | `id`, `kind`, `error` |
| `information_below_floor` | below the floor for its scope (§7.8) | `id`, `kind`, `pattern`, `length`, `information_bits`, `anchor_information_bits`, `error` |
| `scope_unsupported` | `scope: "suffix"` on a wrapped PRIMARY graph, L ≤ k (§7.2) | as above |

### 8.10 Notes

| note | meaning |
|---|---|
| `low_complexity_pattern` | an exact pattern that sdust flags (T = 20, W = 64, the seeder's parameters): why its counts may be large |
| `strand_unknown_canonical` | a CANONICAL or PRIMARY graph: orientations, not strands |
| `paths_later_increment` | L > k: anchors counted, paths neither extended nor extracted |
| `annotation_unbudgeted` | increment 3: the labels were read without the budget-aware decode (`allow_unbudgeted_annotation`): no memory bound on the reads themselves (what they returned is in the account), the deadline checked between chunks of keys |
| `record_bounds_unknown` | increment 3: coordinates without a record mapping (`placement: "global"`): occurrences are (`kmer_coord`, `offset`), placed in no record, not deduplicated, not counted |
| `annotation_not_read` | increment 3: the request asked for labels (`output.labels: "all"`) or named an annotation field, and this answer reads none (mode `count`, or `labels: "none"`): the fields had no effect |

## 9. What an answer licenses

- `retrieval_complete: true` (retrieval modes only): every graph context of the pattern in its scope and
  strands is in `results`, so a k-mer absent from them contains no instance there. It claims nothing about any
  label, sample or record: none were read.
- An `exact` count is the number of graph contexts (or anchors) in the scope; `exact` 0 is an absence claim
  over the index's retained k-mers in that scope (§7.2), never over labels.
- `at_least`, `bounds` and `unknown` claim what they say and no more (§7.4).
- A `count` answer, or any answer with `output.labels: "none"`, claims nothing about labels (design §5.1).
- Increment 3, `output.labels: "all"`, `retrieval_complete: true`: every graph context of the pattern in its scope
  and strands is in `results` with **all** the labels its row carries (`labels_status: "complete"` everywhere),
  each placed where `placement` is `record` or `global` — so a column absent from `by_label` carries no
  occurrence of the pattern in the index's retained k-mers in that scope, and on a `record` index every placed
  occurrence of the listed columns is in the results (`counts.occurrences` exact). Without it (withheld, cut,
  truncated, refused, stopped) the labels and occurrences listed are true but not all: their counts are
  `at_least` or `unknown`, and nothing licenses an absence.

## 10. Capabilities

### 10.1 Where the block is

- **`GET /capabilities`**: on a single-graph server, `"pattern"` in `features`, `routes.pattern =
  "POST /pattern"`, and the `pattern` block. The feature and the route are listed whether or not this graph can
  be searched (as `align` is); `pattern.available` says whether it can. A client gates on both.
- **`GET /traverse/capabilities`** (the document the service's probe reads): the same `pattern` block, so one
  cached probe serves both (design §7.3).
- While the single index loads, `/capabilities` answers with `ready: false` and the block's graph fields `null`
  (`available: null`); `/traverse/capabilities` answers 503.
- On a multi-graph server: no `pattern` feature or route; the block on `/capabilities` and on
  `/traverse/capabilities?graph=NAME` is `{pattern_contract_version: 1, available: false,
  unavailable_reason: "multi_graph_later_increment"}` and nothing else.

### 10.2 The block, field by field

<!-- schema: capabilities -->
| field | type | version 1 | meaning |
|---|---|---|---|
| `pattern_contract_version` | integer | 1 | §1 |
| `available` | boolean \| null | | `true`: `/pattern` answers on this graph; `false`: it never will (`unavailable_reason`); `null`: the index is loading |
| `unavailable_reason` | string \| null | | `mask_required`, `representation_unsupported`, `primary_unwrapped`, `alphabet_unsupported`, `multi_graph_later_increment`; `null` when available or loading |
| `modes` | list | `["count", "all_or_count", "partial"]` | §7.5 |
| `default_mode` | string | `"all_or_count"` | an omitted `mode` |
| `projections` | list | `["none", "all"]` | the `output.labels` values served **now**; gate label projections on this list (`"all"` since increment 3; `["none"]` before) |
| `default_projection` | string | `"none"` | an omitted `output.labels` (§1) |
| `projections_later_increment` | list | `["predicate_only"]` | refused with `later_increment` today (`["all", "predicate_only"]` before increment 3) |
| `default_occurrences` | boolean | `true` | increment 3: an omitted `output.occurrences` with `labels: "all"` |
| `kinds` | list | `["dna", "iupac"]` | pattern kinds served |
| `kinds_later_increment` | list | `["protein"]` | |
| `default_scope` | string | `"any_offset"` | |
| `scopes_by_graph_mode` | object | `basic`, `canonical`: `["suffix", "any_offset"]`; `primary`: `["any_offset"]` | the rule |
| `scopes` | list \| null | | this graph's requestable scopes |
| `long_patterns` | string | `"anchors_counted"` | what L > k gets (§7.7) |
| `strands` | list | `["both", "forward", "reverse"]` | |
| `default_strands` | string | `"both"` | |
| `graph_cleaned` | string | `"unknown"` | whether graph cleaning may have pruned k-mers (design §3); not known in version 1 |
| `records_shorter_than_k` | string | `"not_indexed"` | such records have no k-mer |
| `resident_only` | boolean | `true` | the route never loads an index (`in_ram` refused) |
| `caps` | object | | the maxima (§4.5): `max_contexts`, `max_anchors`, `max_steps`, `time_budget_ms`, `min_information_bits` (the floor), `max_patterns`; increment 3: `max_labels_per_anchor`, `max_annotation_work`, `max_memory_mb`, `max_labels`, `max_occurrences_per_label` |
| `default_time_budget_ms` | number | 60,000 | the budget of a request that names none, below `caps.time_budget_ms` |
| `finalize_reserve_ms` | number | 250 | §7.6 |
| `caps_rule` | string | | the clamp rule in prose |
| `graph_mode` | string \| null | | `basic`, `canonical`, `primary` |
| `k` | integer \| null | | |
| `alphabet` | string \| null | | `$ACGT` or `$ACGTN` |
| `strand_stated` | boolean \| null | | `true` on BASIC |
| `mask` | string \| null | | `file`: the `.edgemask` loaded beside the graph; `built_at_load`: built in memory at start-up (`--pattern-build-mask`, milestone 1b); `absent`: the route answers `mask_required` |
| `placement` | string \| null | | the placement `output.labels: "all"` gives on this index (`record`, `global`, `none`, `none_canonical`, §14.3; before increment 3: what a later increment could give) |
| `support` | string \| null | | the best per-label support of a path (`record_verified`, `label_intersection`; milestone 4); a context of L ≤ k has `kmer` |
| `annotation` | string \| null | | `budgeted` (row-diff with budgeted decode) or `unbudgeted`: the reads of `output.labels: "all"` (§14.4; `unbudgeted` needs `allow_unbudgeted_annotation`); `count` and `none` never read it |

The graph fields (`graph_mode` to `annotation`) are `null` while the index loads; when the graph is not
recognised (`representation_unsupported`, `primary_unwrapped`) only `k` is set; `placement`, `support` and
`annotation` are `null` if the annotation could not be described.

<!-- schema: capabilities_multi -->
| field | type | meaning |
|---|---|---|
| `pattern_contract_version` | integer | 1 |
| `available` | boolean | `false` |
| `unavailable_reason` | string | `multi_graph_later_increment` |

### 10.3 How a client gates

- No block, or a version below 1: no pattern search on this host (`pattern_unsupported` on the service).
- `available: false`: the host has the route but not for this graph; `unavailable_reason` says why (the
  `mask_required` host becomes available once its operator creates the mask).
- `available: null`: ask again after the load.
- Offer a projection, kind or scope only when its list has it; send values within `caps`.
- `output.labels: "all"`: offer it where `projections` has it; on `annotation: "unbudgeted"` it needs
  `allow_unbudgeted_annotation: true` (else 400 `annotation_unbudgeted`); `placement` says whether occurrences
  come in records (`record`), as coordinates (`global`) or not at all.

## 11. Fixtures

`api/python/tests/data/traverse/pattern/` holds one directory per situation with the body sent
(`request.json`; `null` for a GET) and the body answered (`answer.json`), generated by
`scripts/traversal/pattern_fixtures.py` from a server of this build on copies of the mini index
(`build/mini_refseq`, a BASIC index of refseq33m's format at k = 31: its graph given the `.edgemask` by
`metagraph transform --mask-dummy`, as a host gets it; served as built for `mask: absent`, and with
`--pattern-build-mask` for `mask: built_at_load`; served with `--no-coord-mapping` for `placement: global`), and a
PRIMARY index of two of its record files (a column annotation: `annotation: unbudgeted`). `index.json` states each
fixture's server, method, path and status; `README.md` says in one line what each shows. They cover the
capabilities on both routes (each `mask` value: `file`, `built_at_load`, `absent`; multi-graph; PRIMARY), each
mode, each `withheld` and `cut` reason, both scopes and every strand setting, IUPAC patterns, a palindrome, a
pattern longer than k, the three error slots, the clamps, and every whole-request refusal (two 503 bodies are
hand-made: no server produces them on demand). Increment 3 adds `labels_all` (record placement within the
threshold: 9 columns, 42 placed occurrences of the blaNDM-1 forward primer, equal to a scan of the mini's FASTA),
`labels_all_withheld` (above the threshold: nothing read), `labels_all_truncated` and
`labels_all_truncated_partial` (`max_labels_per_anchor` 2), `labels_all_partial` (the label and occurrence cuts),
`labels_all_global` (no record mapping), `labels_all_count` (mode `count`: `annotation_not_read`),
`annotation_unbudgeted` (the 400 on the PRIMARY index's column annotation) and `labels_all_unbudgeted` (the same
with `allow_unbudgeted_annotation`); `later_increment_labels` now asks for `"predicate_only"` (`"all"` is served).

- `pattern_fixtures.py --check` regenerates them and compares: bodies byte for byte, except
  `timing.elapsed_ms` (and every other `timing` value of an entry: `label_discovery_ms`, `placement_ms`),
  `server_instance` and the counts and work of a `time_limited` entry, which vary between runs. The labelled
  answers' `work.annotation_units` and `work.memory_bytes` are deterministic (§14.4) and compared.
- `api/python/tests/test_pattern_fixtures.py` validates every answer against the field lists of §§6, 8 and 10
  (the tables marked `schema`), so that this document, the fixtures and the server cannot drift apart silently.
- `integration_tests/test_pattern.py` (`TestPatternFixtures`) runs `--check` with the binary under test, so a
  build whose `/pattern` or capabilities answers differ from the stored bodies fails its own integration tests.

## 12. What changes in later increments

Contract version 1 stays; each row adds, it does not change (§1). The capabilities announce each addition before
a client may use it. (Milestone 1b, the edge mask, is in this build and its fixtures: `mask: "built_at_load"` with
`--pattern-build-mask`, and the `mask_required` message naming `transform --mask-dummy` and that flag.)

| milestone (design §13) | request | answer | capabilities |
|---|---|---|---|
| 3: labels and placement — **served in this build (§14)** | `output.labels: "all"`; `max_labels_per_anchor`, `max_annotation_work`, `max_memory_mb`, `max_labels`, `max_occurrences_per_label`, `allow_unbudgeted_annotation`; `output.occurrences` | `labels` on each result (column, support, placed occurrences `seq_id`, `record`, `strand`, `nt_coords`, `nt_length`; or `kmer_coord`, `offset`); `by_label`; `counts.labels` and `counts.occurrences` known; `rows_refused`, `anchors_truncated`, `labels_cut`, `occurrences_cut`; `withheld` reasons `annotation_budget`, `anchor_labels_truncated`, `output_budget`; notes `record_bounds_unknown`, `annotation_unbudgeted`, `annotation_not_read` (`label_intersection_only` arrives with paths, milestone 4) | `projections` gains `"all"`; `default_occurrences`; five caps; `default_projection` stays `"none"` (§1) |
| 4: patterns longer than k | `max_paths`, `require_support`, `output.paths: true` | results for L > k are paths (`kmer` = the L spelled bases, `anchor_kmer`, node path with `output.paths`); `counts.paths` known with `candidates_examined`; per-label `support`; `withheld` `anchors_above_threshold`; `paths_later_increment` no longer appears | `long_patterns` changes from `"anchors_counted"` |
| 5: peptides | `patterns[i].protein`, `genetic_code` | `kind: "protein"` (instances name the codons) | `kinds` gains `"protein"` |
| 5b: predicates | `predicate`, `max_predicate_contexts`, `max_predicate_work`, `output.labels: "predicate_only"` | `selection` (tested, selected, access, unknown labels); `withheld` `predicate_above_threshold`, `predicate_budget`; filtered answers state their narrower absence claim (design §5.1) | `projections` gains `"predicate_only"` |
| 6: multi-graph | `graphs` (as `/search` names graphs and chunks), `budget_split` | each result carries `graph`, `index_fp`, `release`; counts `by_shard` with `per_shard`; `stop` and `withheld` gain the shard; the merged order of design §8 | the block on multi-graph servers becomes available, with the resident graphs |

- What stays: every field of §8 with its type; the relations and their algebra; `retrieval_complete` as the
  only absence licence; the order of §7.9; the refusal envelope `{error, code}`.
- What a version-1 client sees on a later host: more fields and values, read by presence or passed through.

## 13. Contract deltas against the design's draft (resolved)

The design's §7 is a draft of the full contract; milestone 1 built the following, which version 1 freezes.

| design | built and frozen | why |
|---|---|---|
| §7.1: `output.labels` defaults to `all` | the default is the capabilities' `default_projection`, `"none"` in milestone 1 | `all` is not served; the owner's instruction for milestone 1. Stated in every answer (`output`) and in the capabilities; clients send it explicitly (§1) |
| §7.2: `by_offset`, `by_strand` and `suffix` as plain numbers in the example | full counts `{value, relation, unit}` | every count carries its relation (design §3); a per-offset count can be `exact` while the total is `at_least` |
| §5.2, §5.3: `withheld: <reason>`, `stop: {phase, shard}` | `withheld: {reason}`, `stop: {phase, reason}` | objects leave room for the shard (milestone 6); the stop's reason is what an agent turns |
| §7.2: no `returned`, no `cut` | `returned` and `cut: {reason}` in the retrieval modes | `partial` states every cut (design §5.2) |
| §7.2: identity per result (`graph`, `index_fp`, `release`) | one top-level `index` | single graph; per-result identity arrives with milestone 6 |
| §7.2: `placement` per entry | not in the entry; the capabilities' `placement` states what a later increment could give | nothing is placed in version 1 |
| §7.2: `work.annotation_rows`, `work.memory_bytes`; `timing.phase1_ms`, … | `work: {ranges_visited, mask_scans, steps}`; `timing: {elapsed_ms}` | no annotation read and no memory account in version 1 |
| §7.2: `strands: ["+", "-"]` on every index | strand symbols on BASIC, orientation names elsewhere | no strand is known on canonical indexes (design §3) |
| §5.5: contexts in BOSS edge order, then offset | then orientation | an IUPAC pattern and its reverse complement can share one (k-mer, offset) |
| §7.1: whole-request 400s named generically | codes `invalid_request`, `later_increment`, `resident_only`, `mask_required`, the graph reasons; 503 `deadline` | one envelope `{error, code}` |
| §4.1: `scope_unsupported: suffix_on_primary` | slot code `scope_unsupported` | the message names `any_offset` |
| §7.3: the block's fields | as §10.2, plus `available`, `unavailable_reason`, `default_*`, `*_later_increment`, `caps_rule`, `long_patterns` | a client gates on what is served now and sees what is coming |
| §4: `mask: file \| built_at_load \| absent` | as designed: `file` \| `absent` at `b570800d`, `built_at_load` added by milestone 1b | `--pattern-build-mask` landed with 1b, no field changed |
| §5.3: `max_steps` per shard | per request (one shard) | single graph; per-shard shares with milestone 6 |
| §7.2 (increment 3): a label's placed occurrences in `occurrences` | `occurrence_list` (the list) beside `occurrences` (its count `{value, relation, unit}`, as in `by_label`) | one name, one type: `occurrences` is a count everywhere in the answer |
| §7.2 (increment 3): `column` and `properties` by `get_label_as_json` | `column` is the annotation column's label as stored; no `properties` | the label is not parsed (a label with `;` or `=` stays whole); a later version may add `properties` |
| §7.2 (increment 3): `by_label` `contexts`, `contexts_suffix` as numbers | counts `{value, relation, unit: graph_contexts}` | a label's count is `at_least` when a row was cut or refused |
| §7.2 (increment 3): a result's labels are its labels | each result states `labels_status` and `labels_total` | a truncated, refused or unread row is told apart from a row without labels (`labels: null` unless read) |
| §5.3, §7.2 (increment 3): `memory.stop` with the phase | `stop: {phase, reason: max_memory}` (phases `label_discovery`, `placement`, `output`) and `work.memory_bytes` | the one `stop` object of the entry |
| §5.2 (increment 3): `max_annotation_work` gives `rows_refused` | a work stop is `stop: {phase, reason: max_annotation_work}`, the unread rows `labels_status: "not_read"`; `rows_refused` holds the rows the memory account refused, each with what it needed | a refused row is a property of the row; a work stop is not |
| §4.3 (increment 3): placement `global` with `kmer_coord` + `offset` | as designed, and nothing deduplicated or counted (`counts.occurrences` unknown) | without record bounds, `kmer_coord + offset` of two records can coincide (records are concatenated), so equal sums are not one occurrence |

## 14. Increment 3: `output.labels: "all"` (labels and placement)

Served by this build (`src/cli/pattern_retrieval.cpp`, design §4.3, §5.2–§5.5, §7.2) as additions to version 1.
`projections` lists `"all"`; a client gates on it (§10.3).

### 14.1 What is read, and when

- Only in a retrieval mode, and only for a pattern whose release the engine did not withhold: a pattern whose
  count was not admitted (`withheld` for any reason of milestone 1) reads no annotation row (design §5.2:
  `work.annotation_rows` 0). Mode `count` reads none whatever the projection (note `annotation_not_read`).
- A pattern longer than k releases nothing (`paths_later_increment`, §7.7) and reads nothing; without anchors its
  empty answer is complete, with `counts.labels` exact 0.
- The patterns are read in request order, under one memory account and one work budget for the whole request
  (§14.4); the deadline is the request's (§7.6).

### 14.2 Step 1: label discovery

- **The rows** are the distinct annotation rows of the released contexts (their `row`, §7.10), in the answer order
  of their first context: a k-mer held by two contexts (two offsets, or both orientations of an IUPAC pattern) is
  read once.
- Each row's labels are read with the traversal's `LabelRecorder` (column labels): at most
  `max_labels_per_anchor` of them, the first in ascending annotation-column id (not the names' order), and the
  row's true total beside them. A label is an annotation column (a taxid on refseq33m, a sample, or a header on an
  index annotated by header); `column` is its label as stored.
- **A row with more labels than the cap** is truncated and listed once in `anchors_truncated` with its total.
  `all_or_count` reads every row first (so that every truncated row is listed and the cap to ask for is known)
  and then withholds: `anchor_labels_truncated`. `partial` returns its contexts with `labels_status:
  "truncated"`, the kept labels and `labels_total`.

<!-- schema: anchor_truncated -->
| field | type | meaning |
|---|---|---|
| `kmer` | string | the row's k-mer (the first context's that has it) |
| `row` | integer | §7.10 |
| `cap` | integer | `max_labels_per_anchor` |
| `total` | integer | the labels the row carries |

### 14.3 Step 2: placement

- When `placement` is `record` or `global` (and `output.occurrences` is not `false`): the coordinates of every
  row read with at least one label, for the labels discovered so far, with the traversal's `LabelQuery`.
  `all_or_count` places only when every row was read completely.
- **`record`** (BASIC, coordinates, the `.seqs` record mapping): for a context of offset p and each k-mer
  coordinate c of a label (column) in its row, `(seq_id, local)` is the record mapping of c **first**, the offset
  is added **after**: start = local + p + 1 (1-based), end = start + L − 1 (`nt_coords` "start-end", the format of
  `metagraph align --json`), `nt_length` the record's length. The strand is the context's: `+` (P on the record),
  `-` (rc(P) on the record, so P on its − strand, the same bases), `=` (palindromic). With k = 5 and a column of
  the records `ACGTA`, `CCCCC`, `CCCCC` (one k-mer each, column coordinates 0, 1, 2), `GTA` at offset 2 of `ACGTA`
  is record 0, `3-5` — not the record of coordinate 0 + 2.
- **Deduplication** (design §5.4): a placed occurrence is (column, `seq_id`, start, strand). Each context lists
  the occurrences its k-mer holds, so the k − L + 1 contexts of an interior occurrence each list it; `by_label`
  and `counts.occurrences` count it once. Two offsets in one k-mer are two occurrences; a k-mer repeated in one
  record is one context with several occurrences; two records with hits at equal local coordinates are told apart
  by `seq_id` (the header need not be unique).
- **`global`** (BASIC, coordinates, no record mapping: an index without its `.seqs`, or a server started with
  `--no-coord-mapping`): each occurrence is the column coordinate of the context's k-mer and the offset, placed in
  no record. Records are concatenated in a column, so `kmer_coord + offset` of two records can coincide: nothing
  is deduplicated and nothing counted (`counts.occurrences` and every label's `occurrences` `unknown`; note
  `record_bounds_unknown`).
- **`none`** (no coordinates) and **`none_canonical`** (CANONICAL and PRIMARY indexes: a stored k-mer may be the
  reverse complement of the deposited one, so neither strand nor offset is known): labels only. **`not_requested`**:
  `output.occurrences: false`.

<!-- schema: occurrence -->
| field | type | meaning |
|---|---|---|
| `seq_id` | integer | the record's index in its column (the `.seqs`: the order the records were annotated in, FASTA order) |
| `record` | string | the record's header, for display (not an identity) |
| `strand` | `"+"` \| `"-"` \| `"="` | the context's |
| `nt_coords` | string | "start-end", 1-based, inclusive: the instance's bases in the record |
| `nt_length` | integer | the record's length (nt): its k-mers in the `.seqs` + k − 1 |

<!-- schema: occurrence_global -->
| field | type | meaning |
|---|---|---|
| `kmer_coord` | integer | the column coordinate (0-based) of the context's k-mer: its position among the column's concatenated records' k-mers |
| `offset` | integer | the pattern's offset in that k-mer (the context's) |
| `strand` | `"+"` \| `"-"` \| `"="` | the context's |

### 14.4 Budgets, the memory account and the deadline

- **Access.** `annotation: "budgeted"` when the annotation has the budget-aware decode
  (`LabelOracle::decode_charged()`: the row-diff family, e.g. `row_diff_brwt_coord` of refseq33m and the mini
  index): every read takes a `DecodeBudget` of what the memory account has left and admits each row against its
  demand, all or nothing per read (a read refused at a row is repeated for the rows before it). A row that does not
  fit alone is **refused**: listed once in `rows_refused`, `labels_status: "refused"`; `partial` reads on,
  `all_or_count` stops reading and withholds (`annotation_budget`). `annotation: "unbudgeted"` (column, BRWT, row
  and disk annotations): refused (400 `annotation_unbudgeted`) unless `allow_unbudgeted_annotation: true`; then
  the reads run without a memory bound of their own, what they return is charged to the account afterwards (a
  read it cannot hold stops the reads: `stop {phase, max_memory}`), and every entry carries the note
  `annotation_unbudgeted`. The label names such a read returns stay with the dictionary even when the rest does
  not fit, so the account can pass its maximum by the names of one read (`work.memory_bytes` shows it); while it
  is past its maximum nothing more fits: the reads stop and the later patterns' contexts are not admitted.
- **The memory account** (`max_memory_mb`, one per request, over all its patterns) is a deterministic model in
  bytes, never a measurement, so that where a request stops does not depend on the allocator: a released
  context 512 + 2k (its descriptor and its result object), a dictionary label 192 + 3 × its name's length
  (charged inside the read that names it), the rows' label lists and coordinates as the `DecodeBudget` charges
  them (held until the pattern's labels are built), a statement of a row (a `rows_refused` or
  `anchors_truncated` entry) 384 + k, a label of a result 256, a placed occurrence 256 + its record name's length
  (`global`: 192), an occurrence in a label's deduplication set 64. What the reads and the deduplication sets
  hold is freed after each pattern; the dictionary, the descriptors, the statements and the labels built stay
  (the buffered answer). The label caches of the reused classes get no allotment (nothing is cached). The order:
  - the descriptors, each charged as the engine releases its context and before its result object is built: in
    `all_or_count` they may take what the account has left, in `partial` at most half of it (the other half is
    kept for the pattern's reads and labels, so that a memory cut still returns labelled contexts). The first
    that does not fit ends the list and no result object is built after it: `all_or_count` withholds
    (`output_budget`), `partial` returns the contexts admitted (`cut: max_memory`, `stop {output, max_memory}`);
  - then the reads, each reserving one statement per row it reads before it starts (so that a refusal or a
    truncation can always be stated); when the account cannot reserve one, the reads stop: `stop
    {label_discovery | placement, max_memory}`, the rows not read `labels_status: "not_read"`, `all_or_count`
    withheld (`annotation_budget`);
  - then the labels of each context, in answer order — when a context's do not fit, `stop {output,
    max_memory}`: `all_or_count` withholds (`output_budget`), `partial` returns it and the later ones with
    `labels_status: "output_budget"`.

  `work.memory_bytes` is the account's peak so far. Apart from the forced names of an unbudgeted read (above),
  it never passes `max_memory_mb`.
- **Work** (`max_annotation_work`, the oracle's units, one budget per request): a row read costs 8, plus 1 per
  label of the row (its true total, also when truncated) in step 1, plus 1 per label and 1 per coordinate in
  step 2, plus the units of its row-diff dependency rows (8 per row and 1 per entry they store); a row read in
  both steps costs both. The reads take 1, 2, 4, … up to 512 rows, and near the budget no more than it has left
  at the widest row read so far; the budget is checked before each read, so the last one can pass it by about a
  row. The stop: `stop {label_discovery | placement, max_annotation_work}`, the rows not read `labels_status:
  "not_read"`, `all_or_count` withheld (`annotation_budget`); the request's later patterns read nothing.
- **The deadline.** The reads are work (§7.6): the work time is checked before each read and, within a read,
  between its chunks (paced at `--traverse-chunk-target-ms`, as `/traverse`'s reads; an interrupted read returns
  nothing). A time stop: `stop {label_discovery | placement, time}`, `determinism: time_limited`, `all_or_count`
  withheld (`deadline`), `partial` returns what was read; the later patterns answer as after any time stop (§7.6).
- **Determinism.** Apart from `timing` and a time stop, the labels, the cuts, the stops, `annotation_units` and
  `memory_bytes` depend only on the index and the request.

<!-- schema: row_refused -->
| field | type | meaning |
|---|---|---|
| `kmer` | string | the row's k-mer |
| `row` | integer | §7.10 |
| `phase` | `"label_discovery"` \| `"placement"` | the read that refused it |
| `reason` | `"max_memory"` | the account could not hold the row's read |
| `needed_bytes` | integer \| null | the least the row alone was seen to need (its demand, with the names it would give) |
| `available_bytes` | integer | what the account had left for it |

### 14.5 The answer

- **Each result** (a released context) gains `support: "kmer"`, `labels_status`, `labels_total` and `labels`.
  `labels_status`: `complete` (every label of the row), `truncated` (§14.2), `refused` (§14.4), `not_read` (a stop
  came first, or `all_or_count` stopped reading at a refused row), `output_budget` (read, but the account could
  not hold its labels in the answer). `labels` is `null` unless the row was read and held, and then lists the
  context's labels in label order.
- **Label order** (design §5.5): (contexts desc, column asc), over the returned contexts — `by_label`'s order.
- **`partial`'s lists.** `max_labels`: `by_label` and every result list the first `max_labels` labels of the label
  order (`labels_cut`); the counts count all of them. `max_occurrences_per_label`: each listed label lists the
  first `max_occurrences_per_label` occurrences of its deduplicated union, in (`seq_id`, start, strand) order —
  every context lists those of them it holds, possibly none (`occurrences_cut`); the label's counts stay whole.
  `all_or_count` cuts nothing: all or nothing.
- **`by_label`**: one entry per label over the returned contexts; `null` when the results are withheld, `[]` when
  none was found.

<!-- schema: label -->
| field | type | meaning |
|---|---|---|
| `column` | string | the annotation column (its label as stored) |
| `support` | `"kmer"` | the label annotates the context's one k-mer (design §4.3) |
| `occurrences` | count | `record` placement only: this label's distinct placed occurrences in this context, unit `placed_occurrences`, `exact` (also when `partial` lists fewer), `unknown` when the row's placement was refused or not reached |
| `occurrence_list` | list \| null | `record` placement: objects of the `occurrence` table; `global`: of the `occurrence_global` table; `null` when not placed; absent for `none`, `none_canonical`, `not_requested` |

<!-- schema: by_label -->
| field | type | meaning |
|---|---|---|
| `graph` | string \| null | the index's `--index-name` (`index.index_ns`), `null` without one; a shard's name with milestone 6 |
| `column` | string | the label |
| `contexts` | count | the returned contexts carrying it, unit `graph_contexts` |
| `contexts_suffix` | count | of them, those at offset k − L |
| `occurrences` | count | its deduplicated placed occurrences over the returned contexts, unit `placed_occurrences` (`unknown` without `record` placement) |

<!-- schema: labels_cut -->
| field | type | meaning |
|---|---|---|
| `reason` | `"max_labels"` | |
| `returned` | integer | the labels listed |

<!-- schema: occurrences_cut -->
| field | type | meaning |
|---|---|---|
| `reason` | `"max_occurrences_per_label"` | |
| `labels` | integer | how many labels list fewer occurrences than they have |

### 14.6 Counts and what they license

- `counts.labels`: the distinct labels of the returned contexts' rows that were read (also those a `max_labels`
  cut does not list): `exact` iff every context of the pattern was returned and every row read completely,
  `at_least` otherwise, `unknown` when the results are withheld or nothing was returned of an incomplete release.
  A complete empty release: `exact` 0.
- `counts.occurrences` (`record` placement only): the deduplicated placed occurrences summed over labels: `exact`
  iff `counts.labels` is, every label of every row was placed and the answer held them; `at_least` when some
  were placed; `unknown` otherwise, and always for `global`, `none`, `none_canonical`, `not_requested`.
- `by_label`'s `contexts` and `contexts_suffix` are `exact` iff `counts.labels` is; its `occurrences` as
  `counts.occurrences`, per label.
- `retrieval_complete` (also §7.5): the release was complete, every row was read completely, every label placed
  where placement applies, nothing cut (`labels_cut`, `occurrences_cut`, `cut`) and no annotation stop.
  `all_or_count` withholds otherwise, naming the first that applies of `deadline`, `annotation_budget` (a refused
  row, a work stop, an unbudgeted read the account could not hold), `output_budget`, `anchor_labels_truncated`;
  `rows_refused` and `anchors_truncated` stay in the withheld entry to say why. What a complete answer licenses
  is in §9.

### 14.7 Example (fixture `labels_all`, abridged)

```json
{"id": "NDM-F", "kind": "dna", "pattern": "GGTTTGGCGATCTGGTTTTC", "placement": "record",
 "annotation": "budgeted",
 "counts": {"contexts": {"value": 24, "relation": "exact", "unit": "graph_contexts", "...": "..."},
            "labels": {"value": 9, "relation": "exact", "unit": "labels"},
            "occurrences": {"value": 42, "relation": "exact", "unit": "placed_occurrences"}},
 "retrieval_complete": true, "withheld": null, "returned": 24, "rows_refused": [],
 "anchors_truncated": [], "labels_cut": null, "occurrences_cut": null,
 "results": [
   {"kmer": "CCAACGGTTTGGCGATCTGGTTTTCCGCCAG", "instance": "GGTTTGGCGATCTGGTTTTC", "offset": 5,
    "strand": "+", "node": 482091, "row": 482090, "support": "kmer", "labels_status": "complete",
    "labels_total": 9,
    "labels": [
      {"column": "1296536", "support": "kmer",
       "occurrences": {"value": 1, "relation": "exact", "unit": "placed_occurrences"},
       "occurrence_list": [{"seq_id": 0, "record": "NZ_LPPQ01000025.1", "strand": "+",
                            "nt_coords": "3682-3701", "nt_length": 11137}]},
      "... 8 more"]},
   "... 23 more"],
 "by_label": [
   {"graph": "mini_refseq", "column": "1296536",
    "contexts": {"value": 24, "relation": "exact", "unit": "graph_contexts"},
    "contexts_suffix": {"value": 2, "relation": "exact", "unit": "graph_contexts"},
    "occurrences": {"value": 2, "relation": "exact", "unit": "placed_occurrences"}},
   "... 8 more"],
 "work": {"ranges_visited": 122, "mask_scans": 0, "steps": 122, "annotation_rows": 48,
          "annotation_units": 24440, "memory_bytes": 220428},
 "timing": {"elapsed_ms": "...", "label_discovery_ms": "...", "placement_ms": "..."}}
```

## 15. What changed since `bd44e597`

Contract version 1 stays: every change below adds (§1).

- **Request.** Accepted now: `output.labels: "all"`; `output.occurrences` (with `"all"`; default `true`);
  `max_labels_per_anchor`, `max_annotation_work`, `max_memory_mb`, `max_labels`, `max_occurrences_per_label`
  (capped and clamped like the milestone-1 caps), `allow_unbudgeted_annotation`. Still `later_increment`: the
  fields of §4.4, `output.labels: "predicate_only"`, `output.paths: true`.
- **Refusals.** New: 400 `annotation_unbudgeted` (§6, check 10 of §5). Changed code: `output.occurrences: true`
  without `output.labels: "all"` is 400 `invalid_request` (it was `later_increment`); the six annotation fields
  with a wrong value are 400 `invalid_request` (they were `later_increment` whatever their value).
- **Answers to requests with `output.labels: "all"`** (new): `output.occurrences`; six `limits` fields; the entry
  fields `placement`, `annotation`, `by_label`, `rows_refused`, `anchors_truncated`, `labels_cut`,
  `occurrences_cut`; `work.annotation_rows`, `annotation_units`, `memory_bytes`; `timing.label_discovery_ms`,
  `placement_ms`; the result fields `support`, `labels_status`, `labels_total`, `labels`; `counts.labels` and
  `counts.occurrences` known; values `withheld.reason` `annotation_budget`, `anchor_labels_truncated`,
  `output_budget`; `cut.reason` `max_memory`; `stop.phase` `label_discovery`, `placement`, `output`;
  `stop.reason` `max_annotation_work`, `max_memory`; notes `annotation_unbudgeted`, `record_bounds_unknown`;
  `placement` `not_requested`.
- **Answers to the other requests.** A request valid at `bd44e597` is answered byte for byte as by the
  milestone-1b build (apart from `timing`): checked on every unchanged fixture of §11 and on a panel of 23
  count and `labels: "none"` requests (every mode, scope and strand setting, the step and threshold stops) against
  that build's binary, and again by the integration's server panel (54 requests, plain and gzip; the only
  differences were the two stated ones: mode `count` with `labels: "all"`, refused before, and the prose of the
  `predicate_only` refusal, whose code is unchanged). A request that names an annotation field or `output.labels: "all"` and reads no annotation
  (mode `count`, `labels: "none"`) gets the note `annotation_not_read`.
- **The memory account, after the review of increment 3** (§14.4): an account past its maximum (by the names
  an unbudgeted read returned) admits nothing more (before, `max - held` wrapped and every later item fitted);
  a context's descriptor is charged as the engine releases it, before its result object is built, and
  `partial`'s descriptors take at most half of the account; the statements of refused and truncated rows are
  priced (384 + k) and reserved before the read. A labelled answer's `work.memory_bytes` grows by 384 + k per
  `rows_refused` and `anchors_truncated` entry; nothing else changed in the fixtures of this increment.
- **Capabilities** (both routes): `projections` `["none", "all"]`, `projections_later_increment`
  `["predicate_only"]`, `default_occurrences`, the five caps in `caps`; `placement`, `support` and `annotation`
  now describe what is served. Nothing else changed in `/capabilities` or `/traverse/capabilities`.
- **Server flags.** `--pattern-max-labels-per-anchor` (64), `--pattern-max-annotation-work` (100,000,000),
  `--pattern-max-memory-mb` (256), `--pattern-max-labels` (1,000), `--pattern-max-occurrences` (16); the reads use
  `--traverse-chunk-target-ms` (50). `metagraph pattern` answers as the server.
- **Fixtures** (§11): new `labels_all`, `labels_all_withheld`, `labels_all_truncated`,
  `labels_all_truncated_partial`, `labels_all_partial`, `labels_all_global`, `labels_all_count`,
  `annotation_unbudgeted`, `labels_all_unbudgeted`, and the server `masked_no_map`; changed: the seven
  capabilities bodies of single-graph servers (the additions above), `later_increment_labels` (now
  `"predicate_only"`), `README.md`, `index.json`; `--check` blanks every `timing` value of an entry. Every other
  fixture is unchanged.
