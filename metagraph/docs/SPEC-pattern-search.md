# Spec: pattern search — `POST /pattern`, contract version 1

**Status:** frozen for milestone 1 (increments 0–2 of `DESIGN-pattern-search.md`), 2026-10-07. Written from the
code as committed on `gr/labeled-traversal` at `b570800d` (route) and `63acdd7b` (engine); milestone 1b (the edge
mask: `transform --mask-dummy`, `--pattern-build-mask`, `mask: built_at_load`, the `mask_required` message naming
both remedies) added no field. **Increment 3** (2026-10-07, after `bd44e597`; `src/cli/pattern_retrieval.cpp`)
adds the projection `output.labels: "all"` to contract version 1 — additions only (§1): its request fields, answer
fields and values are in the tables below and described in §14; when it was added, every answer to a request
that does not ask for it was checked byte for byte against the milestone-1b build's (§15). **The review of
2026-10-07** (milestones 1 and 1b) corrected this document sentence by sentence and changed some answers;
§16 lists every change (one new field, `min_anchor_information_bits`; two new refusal reasons, `mask_invalid`
and `alphabet_untested`; one reserved request field, `long_search`; no field changes meaning). **Increments 4 and 5**
(2026-10-08) add two opt-ins to contract version 1, additions only (§17): the paths of a pattern longer than k
for a request that sets `long_search: "paths"`, with the labels of each path and their support (§12.1), and
protein patterns, peptides searched as their codon automaton (§12.2); every answer to a request that uses
neither was checked byte for byte against the build of `44583b51` (§17). **The owner's decisions of 2026-10-08**
(#16, #17, #19, #21; §18) add to contract version 1, additions only: a graph loaded **without its dummy-edge
mask** is answered instead of refused (`mask_required` is retired): its counts are exact where the search can
prove them and otherwise `bounds` with an additive `estimate`, its lists stay exact, and the capabilities and the
answer say how the graph counts (`counting`, `dummy_fraction`); the mask (and the Bloom filter) are derived data
of the graph, outside the index identity `index_fp`; and the stop `*` is a residue of a peptide (a stop codon of
its genetic code; the slot code `stop_unsupported`, which only `4596bb3b` answered, is retired). Answers on a
masked graph were checked against the build of `4596bb3b` (panels, the stored fixtures and the alignment gate,
§18): the same apart from `timing`, the stop `*` and the capabilities' additions. **Increment 5b** (2026-10-08;
§19) adds, additions only, a request's **predicate** over the annotation columns of each graph context of a
pattern of at most k bases: the contexts it selects, read by a selection pass under its own work budget, with the
projections `"none"`, `"predicate_only"` and `"all"`; a request without a predicate was checked byte for byte
against the build of `9018f41b` (§19.15). Version 1
promises the meaning of every field and count, not identical work from build to build (§1). Checked against
the fixture bodies of §11, which this build answered.
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
- **Compatibility is semantic, not byte for byte** (the owner's decision of 2026-10-07). Between builds that
  state version 1 the meaning of every field and value is kept, and so is what each relation guarantees (§7.4:
  an `exact` count is the count, an `at_least` value a true lower bound, `bounds` a true interval, §9's licences).
  What may change between builds: the work counters (`work`, `timing`), where a budget or time stop falls, and
  therefore which valid bounded values (`at_least`, `bounds`) and which released contexts a stopped search
  states. An `exact` count of the same index and scope is the same number in every build. The same answer, byte
  for byte apart from `timing`, is promised only within one build under one effective configuration (§7.9).
- **What a client does with that** (the rules the service's probe follows, `PROMPT-search-service-pattern.md`
  §3.1 item 6):
  - accept only the contract versions it implements (a version-1 client: exactly 1); refuse a missing,
    malformed (not an integer) or unsupported version as "no pattern search on this host", a higher one
    included: a higher version means that a field changed meaning;
  - within a version, read fields by presence and tolerate fields it does not know (additions);
  - treat an unknown value of `withheld.reason`, `cut.reason`, `stop.phase`, `stop.reason`, a slot's
    `error.code`, a refusal's `code`, a note, `placement`, `support`, `annotation`, `mask`, `counting`,
    `dummy_fraction.source`, `labels_status`, `unavailable_reason`, `selection.pass`, `selection.access` or
    `selection.support` (the extensible enumerations) as "not understood": pass it through, claim nothing from
    it;
  - keep handling a code no build answers any more: `mask_required` (retired by the owner's decision #16 of
    2026-10-08, §6) is still answered by builds of version 1 before it, and the slot code `stop_unsupported`
    (a peptide holding `*`, retired by decision #19, §18) by `4596bb3b`, the one build that had it;
  - gate every option on the capabilities block (§10), never on a milestone number;
  - send `output.labels` explicitly. Its default is stated in the capabilities (`default_projection`) and is
    **frozen at `"none"` for version 1** (the owner's decision of 2026-10-07), also now that `"all"` is served
    (increment 3): a request that names no projection is answered as before, without reading any annotation.
    The design's default (`"all"`, design §7.1) would change what an omitted field means, so a server whose
    default were `"all"` would state a higher version.
- **Messages are prose.** `error` texts and slot `message`s are for people and may change; clients act on
  `code`. So are the capabilities' `caps_rule` and `protein_rule` (references to this document since the owner's
  decision P9, §18): clients act on the fields they describe.

What this build serves (milestone 1, and increment 3 where marked), against the design:

| design | served |
|---|---|
| modes `count`, `all_or_count`, `partial` (§5.2) | all three |
| projections `none`, `all`, `predicate_only` (§4.3, §5.6) | `none`; `all` (increment 3, §14); `predicate_only` with a predicate (increment 5b, §19) |
| predicates (§5.6) | increment 5b (§19): one predicate per request over the annotation columns (`any`, `all`, `none`, `at_least`, `and`, `or`, `not`), for patterns of at most k bases; a pattern longer than k keeps its anchors' answer, its selection `not_started` (the selection of supported paths is a later increment) |
| kinds `dna`, `iupac`, `protein` (§3) | all three; `protein` since increment 5 (§12.2): the 20 amino acids and X, B, Z, J, every NCBI genetic code, and the stop `*` (a stop codon of the genetic code; the owner's decision #19 of 2026-10-08) |
| scopes `suffix`, `any_offset`, `long` (§3) | all; `long` counts anchors only and extracts nothing (§7.7), unless the request sets `long_search: "paths"`: then its paths are counted and released (increment 4, §12.1) |
| graph modes BASIC, native CANONICAL, wrapped PRIMARY (§4.1) | all three; `suffix` refused per pattern on PRIMARY |
| alphabets `$ACGT`, `$ACGTN` (§3) | `$ACGT`; a `$ACGTN` (DNA5) graph is refused, `alphabet_untested` (§6, §10.2), until a DNA5 build passes the pattern tests (the owner's decision of 2026-10-07; §8.2) |
| the dummy-edge mask (§4: required) | graphs with it (`counting: "exact"`: every count of a completed discovery `exact`) and, since the owner's decision #16 of 2026-10-08, without it (`counting: "upper_bound"`: counts `exact` where provable, otherwise `bounds` with an additive `estimate`; lists exact; §7.4, §18). The mask is derived data of the graph, not part of `index_fp` (decision #17, §10.2) |
| labels, placement, occurrences (§4.3) | read only with `output.labels: "all"` (or `"predicate_only"`, §19) in a retrieval mode (increment 3, §14): labels on every index, placement on BASIC indexes with coordinates; otherwise none read and their counts `unknown`. A predicate's selection reads the rows of the contexts it tests in every mode (§19), its own labels only |
| per-label `support` for paths, `require_support` (§4.3) | served with `long_search: "paths"` (increment 4, §12.1): `label_intersection` or `record_verified` per label; a context of L ≤ k has `support: "kmer"` |
| multi-graph servers (§8) | 400 `later_increment`; the block says `multi_graph_later_increment` |
| the deadline with a finalisation reserve, 503 `deadline` (§5.3) | as designed; the annotation reads are work and stop at the work time (§14.4) |

## 2. Terms

| term | meaning |
|---|---|
| pattern | L positions, each a set of bases from {A, C, G, T}: one base (`dna`) or an IUPAC code's set (`iupac`). A pattern's N is {A, C, G, T}; it never matches a record's N symbol (design §3). A peptide (`protein`, §12.2) of m residues is a pattern of L = 3m positions whose bases follow its codon automaton: the bases allowed at a position depend on the bases already spelled in its codon. |
| oriented pattern | P as given (`forward`) or its reverse complement rc(P) (`reverse`); a palindrome (P = rc(P), position by position) is one oriented pattern (`palindromic`). |
| k | the graph's k-mer length (`index.k`; 31 on refseq33m). |
| graph context | (orientation, k-mer, offset): one distinct k-mer of the index that contains the oriented pattern at that 0-based offset (L ≤ k). §7.1. |
| anchor | for L > k: a k-mer that instantiates positions [0, k) of an oriented pattern (design §4.2). Not a context. |
| path | for L > k with `long_search: "paths"` (§12.1): a walk of n = L − k + 1 k-mers of the graph, each the next one's predecessor, that spells an instance of an oriented pattern; its first k-mer is an anchor. A path need not lie in one record. |
| retained island | a maximal run of consecutive k-mer starts of a record whose k-mers the index kept (design §3, "Covered sequence"). Every completeness statement is over retained islands. |
| information bits | Σ log2(4 / \|set_i\|) over the positions: 2 per exact base, 1 per two-base code, log2(4/3) ≈ 0.415 per three-base code, 0 per N. A peptide's: 2 per base less log2 of the distinct strings its codons spell over the positions counted, residue by residue (§12.2): log2(64 / codons) per whole residue. |
| node | the graph's id of a context's k-mer: the BOSS edge index on a DBGSuccinct, the `CanonicalDBG` wrapper id on a wrapped PRIMARY graph (§7.10). |
| row | the annotation row of the context's k-mer, named without reading it (§7.10). |
| step | one range evaluation of the discovery, or one item of a deferred scan: an edge examined among a range's masked edges or its candidates (design §4.1), or on an even-k wrapped PRIMARY graph the palindrome check of one context (§7.6); the unit of `max_steps`. |
| dummy edge | an entry of the BOSS graph that is not a k-mer of the records: a **source dummy** (a k-mer starting with `$`: the padding BOSS adds before a sequence start that no k-mer enters, `$^j` followed by the start's first k − j bases) or a **sink dummy** (its last symbol W is `$`). The dummy-edge mask (`.edgemask`) marks them; a pattern base never matches `$`, so a sink dummy never counts, and a source dummy can hold a pattern only after its `$` run. |
| upper bound U | on a graph without its mask (`counting: "upper_bound"`, §7.4): the candidate entries of a count's ranges, k-mers and source dummies alike; the true count lies between the count's `lower` and U. |
| dummy fraction f | the fraction of real k-mers among the graph's entries whose W is not `$` (the entries a pattern can count), sampled at load on a graph without its mask (`dummy_fraction`, §8.2). |
| estimate | on a graph without its mask, beside a `bounds` count: round(U × f), kept inside [`lower`, U] (§7.4). Not a bound, never `exact`. |

## 3. Transport

- `POST /pattern` with a JSON body: one RFC 8259 JSON text, with unique member names in every object, nothing
  after it (no second value, no comment, no trailing comma) and nested at most 1,000 deep; anything else is 400
  `invalid_request` (§5). The answer is compact JSON, gzip or deflate when the request's `Accept-Encoding` asks
  for it (zlib level `--traverse-compression-level`, 1 by default, to keep compression short; the time to write
  and compress the answer is kept back from the work, §7.6). Error bodies are never compressed.
- The route runs on the server's request pool (`-p`; `--threads-each` is not a pool: it sets the threads used
  within one request and at load), shared with `/search`, `/align`, `/resolve` and `/traverse` (design §9).
  Every route, `GET /capabilities` and `GET /traverse/capabilities` included, runs on these `-p` threads: with
  `-p` long calls in flight, a further call — a capabilities probe too — waits up to their `time_budget_ms`. A
  deployment therefore keeps its concurrent long calls below `-p` (`-p` ≥ cap + 1, as
  `SPEC-labeled-traversal-core.md` "Deployment" requires for traversals), and a client does not take a full pool
  for a dead host.
- The server's content timeout is 900 s; the route's own deadline is `time_budget_ms` (§7.6), at most 600 s by
  default and at most 899,000 ms on `server_query` whatever the flag (the content timeout less 1 s for the
  transport: a larger `--pattern-max-time-ms` refuses to start). `time_budget_ms` is the backend's deadline for
  preparing the answer, counted from the parsed body: not the search time (the work stops earlier, §7.6), not
  the client's end-to-end latency (the wait for a pool thread and the transport are outside it), and not the
  lifetime of a job that wraps the call. A client's HTTP timeout stays above it.
- **A request whose client has left is not answered.** A client that closed, reset or half-closed its
  connection (one that half-closes after sending its request is treated as gone, as on `/traverse`), or a
  request running when the server shuts down: the work ends at the engine's next clock reading (§7.6), the
  writing at its next check, and nothing is written. Nor is an error: a refusal (400), a 503 (the index loading,
  or `deadline`) or any other failure is not written to a client that left or during a shutdown — the server
  asks once more after the answer or the error is built (outside review GPT-2, finding 4). A 503 `deadline`
  still reaches a client that is there. A request abandoned while it waited for a thread is dropped at its first
  such reading or check.
- Single-graph servers only (`server_query -i GRAPH -a ANNOTATION`). A multi-graph server (`server_query
  GRAPHS.csv`) answers every `/pattern` request with 400 `later_increment` (§6).
- `metagraph pattern -i GRAPH -a ANNOTATION [--json] REQUEST.json ...` answers each request file as the server
  would, under the same `--pattern-*` flags: the answer on stdout; a refusal's body on stdout and exit status 1.
  A request file that fails otherwise gets the body the server's 400 without a code has (`{"error": …}`, §6),
  likewise with exit status 1, and the next file is still answered. Answers equal the server's except `timing`
  and what a time stop touched (without `--json` the CLI's indented text counts twice in the estimate of §7.6).

## 4. The request

### 4.1 Top level

<!-- schema: request -->
| field | type | default | rule |
|---|---|---|---|
| `patterns` | list of pattern objects (§4.2) | required | 1 to `caps.max_patterns` (16) entries; more is 400 `invalid_request`, never cut |
| `mode` | `"count"` \| `"all_or_count"` \| `"partial"` | `"all_or_count"` | §7.5 |
| `scope` | `"suffix"` \| `"any_offset"` | `"any_offset"` | a pattern longer than k is searched as `long` whatever is named (§7.7) |
| `strands` | `"both"` \| `"forward"` \| `"reverse"` | `"both"` | `forward` searches P, `reverse` rc(P); a palindrome is searched once whatever is named (§7.3) |
| `stop_at_threshold` | boolean | `false` | stop a pattern's discovery once its running lower bound passes its threshold (§7.5; on a graph without its mask its running upper bound, §7.4); checked in discovery only, so a pattern can still end `exact` above its threshold with no stop |
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
| `long_search` | `"anchors"` \| `"paths"` | `"anchors"` | increment 4 (§12.1): `"paths"` extends every pattern longer than k into its paths; `"anchors"` answers byte for byte as the field's absence (§7.7). A pattern of at most k bases is answered alike under both. Any other value, `null` included, is 400 `invalid_request` |
| `max_paths` | integer ≥ 0 | `caps.max_paths` (1,000) | increment 4: per pattern, the paths an `all_or_count` answer releases at most (and partial's cut, and `stop_at_threshold`'s threshold in the extension, §12.1); lowered like `max_contexts`. Accepted with any request; it acts only with `long_search: "paths"` |
| `require_support` | `"label_intersection"` \| `"record_verified"` | `"label_intersection"` | increment 4: the labels of paths listed: every label carrying the path with its support, or only those one record verifies (§12.1). An annotation field (§4.1, last bullet). `"record_verified"` with `output.occurrences: false` is 400 `invalid_request`; on an index that cannot verify, 400 `support_unavailable` (§6) |
| `genetic_code` | integer | `1` (`default_genetic_code`) | increment 5 (§12.2): the NCBI translation table the request's peptides are read in, one of the capabilities' `genetic_codes` (1–6, 9–16, 21–33); another integer is 400 `genetic_code_unknown`, a value that is not an integer 400 `invalid_request`. Accepted with any request; it acts on protein patterns only |
| `predicate` | object (§19.3) | absent | increment 5b: the condition on the annotation columns a graph context must satisfy to be selected (counted in `counts.selected`, returned); applies to every pattern of the request (§19). Its form is refused with 400 `invalid_request`, more names than `caps.max_predicate_labels` with 400 `predicate_too_large`; with `long_search: "paths"` 400 `invalid_request` |
| `max_predicate_contexts` | integer ≥ 0 | `caps.max_predicate_contexts` (100,000) | increment 5b, per pattern: the compute admission of the selection, the raw contexts it may test (§19.6); lowered like `max_contexts`. Accepted with any request; it acts only with a predicate |
| `max_predicate_work` | integer ≥ 1 | `caps.max_predicate_work` (10⁸) | increment 5b, per **request**: the work of the selection (its row reads, reverse-complement lookups and decisions) in the oracle's units (§19.9), a budget of its own beside `max_annotation_work`; lowered like `max_contexts`. Accepted with any request; it acts only with a predicate |
| `predicate_strands` | `"either"` \| `"context"` | `"either"` | increment 5b (§19.5): on a BASIC graph, `"either"`: a label is present for a context when it annotates the context's k-mer or its reverse complement (the one or the other, never a mix: a single k-mer); `"context"`: the context's own k-mer only. On CANONICAL and PRIMARY graphs one row serves both: evaluated as `"either"` whatever is asked. Accepted with any request; it acts only with a predicate |

- An integer is a JSON number with an integral value (`5.0` is 5). A negative or fractional value is 400
  `invalid_request`, except for `genetic_code`: a fractional value is 400 `invalid_request`, any integer (a
  negative one included) that is not an NCBI translation table id is 400 `genetic_code_unknown` (§6).
- A field this version does not know is 400 `invalid_request` ("unknown field"), never ignored. Fields of
  later increments are refused by name with 400 `later_increment`, whatever their value, `null` included
  (§4.4); `in_ram` with 400 `resident_only`.
- The information floor is server policy (`--pattern-min-information-bits`), not a request field (§7.8).
- The six fields of increment 3 are accepted with every mode and projection; they bound only the reads of
  `output.labels: "all"` in a retrieval mode (with a predicate, `max_memory_mb` and
  `allow_unbudgeted_annotation` also the selection's reads, in every mode, §19.9). An answer that reads no annotation (mode `count`, or
  `output.labels: "none"`, without a predicate) echoes none of them in `limits` (it is the milestone-1 answer)
  and, when the request named one of them or `output.labels: "all"`, carries the note `annotation_not_read` in
  each answered entry (a predicate answer echoes them in every mode and says `projection_not_read` instead,
  §19.10).
  `require_support` (increment 4) is such an annotation field too.
- The fields of increments 4 and 5 are opt-ins: a request that names none of them (nor a protein pattern) is
  answered as before them, byte for byte apart from `timing` (§17). `long_search` and `max_paths` are echoed in
  `limits` only in the answers to `long_search: "paths"`, `require_support` only in those that also read labels;
  `genetic_code` is stated in each protein pattern's entry. A request that names no predicate (§19) is answered
  as before increment 5b, byte for byte apart from `timing`; the three fields `max_predicate_contexts`,
  `max_predicate_work`, `predicate_strands` are then accepted, lowered to their caps (listed in `limits.clamped`)
  and echoed nowhere else.

### 4.2 A pattern

<!-- schema: request_pattern -->
| field | type | rule |
|---|---|---|
| `id` | string | optional; echoed in the pattern's entry (`null` when absent). Another type is 400 `invalid_request` |
| `dna` | string over A, C, G, T (any case) | exactly one of `dna`, `iupac`, `protein` |
| `iupac` | string over the 15 IUPAC codes A C G T R Y S W K M B D H V N (any case) | exactly one of `dna`, `iupac`, `protein` |
| `protein` | string over the 20 amino acids A C D E F G H I K L M N P Q R S T V W Y, the ambiguity codes X B Z J and the stop `*` (any case; the capabilities' `protein_residues`) | increment 5 (§12.2): a peptide, read in the request's `genetic_code`; exactly one of `dna`, `iupac`, `protein`. `*` since the owner's decision #19 of 2026-10-08: a stop codon of the genetic code at that position |

- None or more than one of `dna`, `iupac`, `protein`, or a value that is not a string, is 400 `invalid_request`.
- A string with another character (U, `-`, `.`, whitespace; an IUPAC code in a `dna` pattern; U, O or a digit
  in a `protein` pattern) or an empty string is answered in the pattern's slot with `bad_alphabet` (§8.9); the
  other patterns are answered. The stop `*` is a residue of a `protein` pattern (§12.2).
- No length cap: a pattern longer than k is charged steps for its anchor windows only (§7.7); parsing it, its
  `information_bits` and its palindrome test cost O(L) time inside the deadline that no step charges and no clock
  reading interrupts (a few ns per base; many long patterns can therefore end in a 503 after `time_budget_ms`).
  Its `low_complexity_pattern` diagnostic reads the clock (since the review GPT-3, §7.8, §18): before every piece
  of 128 bases but the first.
- A peptide's length L is in bases (3 per residue): with k = 31, peptides of up to 10 residues are within one
  k-mer; longer ones are patterns longer than k (§7.7, §12.1).

### 4.3 `output`

<!-- schema: request_output -->
| field | type | rule |
|---|---|---|
| `labels` | `"none"` \| `"all"` \| `"predicate_only"` | `"none"`: the label-free projection (design §4.3), contexts without any annotation read. `"all"` (increment 3): every returned context with its labels, placed where the index can (§14); in mode `count` it reads nothing and is said so (note `annotation_not_read`; with a predicate `projection_not_read`). `"predicate_only"` (increment 5b, §19.10): with a predicate, each selected context with the predicate's labels on its own row, placed; without a predicate 400 `invalid_request` in every mode (it was 400 `later_increment` before increment 5b); any other string is 400 `invalid_request` |
| `occurrences` | boolean | with `labels: "all"` or `"predicate_only"`: whether the labels are placed (default `true`, the capabilities' `default_occurrences`; `false` answers `placement: "not_requested"`). With `labels: "none"` (or omitted): `false` is accepted and changes nothing, `true` is 400 `invalid_request` (occurrences are placed per label) |
| `paths` | boolean | accepted with either value and changes nothing (increment 4: a path result always carries its node path, `nodes` and `rows`, §12.1; `true` was 400 `later_increment` before); another type is 400 `invalid_request` |

### 4.4 Fields of later increments

Refused by name (400 `later_increment`), whatever their value: `graphs`, `budget_split`. So that a request
written for a later increment is told what to wait for, not answered as if the field were absent.

`predicate`, `max_predicate_contexts` and `max_predicate_work` were refused so until increment 5b, which serves
them with `predicate_strands` and `output.labels: "predicate_only"` (§19); `"predicate_only"` without a
predicate is now 400 `invalid_request` (it was `later_increment`; the precedent is increment 3's
`output.occurrences: true`, §15). The selection of patterns longer than k (on supported paths, `long_search:
"supported_paths"`) is a later increment: that value is not served (400 `invalid_request`, an unknown value, as
any other), and a predicate with `long_search: "paths"` is 400 `invalid_request`.

`long_search` (the owner's decision of 2026-10-07) was reserved for the increment that extends patterns longer
than k into paths: increment 4 serves it (§12.1), with `max_paths` and `require_support`, and `output.paths:
true` is accepted. The paths are opt-in: a request that does not send `long_search: "paths"` keeps the
anchor-only answer of §7.7 (`"anchors"`, the default, is that answer). Increment 5 serves `patterns[i].protein`
and `genetic_code` (§12.2). (Until increment 3 the list also held `max_labels_per_anchor`,
`max_annotation_work`, `max_memory_mb`, `max_labels`, `max_occurrences_per_label`,
`allow_unbudgeted_annotation`, `output.labels: "all"` and `output.occurrences: true`; until increments 4 and 5,
`long_search`, `max_paths`, `require_support`, `genetic_code`, `patterns[i].protein` and `output.paths: true`.
They are served now, §4.1–§4.3.)

### 4.5 Caps and defaults (server flags)

| flag | default | request field | above it |
|---|---|---|---|
| `--pattern-max-contexts` | 10,000 | `max_contexts` (default = cap) | lowered, listed in `limits.clamped` |
| `--pattern-max-anchors` | 1,000 | `max_anchors` (default = cap) | lowered, listed |
| `--pattern-max-paths` | 1,000 | `max_paths` (default = cap; increment 4) | lowered, listed |
| `--pattern-max-steps` | 100,000,000 | `max_steps` (default = cap) | lowered, listed |
| `--pattern-default-time-ms` | 60,000 | `time_budget_ms` when omitted | — |
| `--pattern-max-time-ms` | 600,000 | `time_budget_ms` | lowered to the cap (not to the default), listed; on `server_query` the flag is at most 899,000 (§3) |
| `--pattern-finalize-ms` | 250 | — (the floor of the time kept back from the work for writing the answer, inside every budget, §7.6) | — |
| `--pattern-delivery-build-mbps` | 10 | — (the rate, MB/s, at which the answer's JSON text is assumed to be built and written, §7.6) | — |
| `--pattern-delivery-compress-mbps` | 50 | — (the rate, MB/s, at which it is assumed to be compressed, §7.6) | — |
| `--pattern-min-information-bits` | 24 | — (the floor, §7.8) | — |
| `--pattern-max-checked-entries` | 50 | — (a graph without its mask: a pattern with at most this many unchecked candidates has each tested, its counts `exact`, §7.4; 0: none; the owner's decision #24) | refused at start-up above 1,000 (or not an integer) |
| `--pattern-max-patterns` | 16 | length of `patterns` | refused (400), never cut |
| `--pattern-max-labels-per-anchor` | 64 | `max_labels_per_anchor` (default = cap) | lowered, listed |
| `--pattern-max-annotation-work` | 100,000,000 | `max_annotation_work` (default = cap) | lowered, listed |
| `--pattern-max-memory-mb` | 256 | `max_memory_mb` (default = cap) | lowered, listed |
| `--pattern-max-labels` | 1,000 | `max_labels` (default = cap) | lowered, listed |
| `--pattern-max-occurrences` | 16 | `max_occurrences_per_label` (default = cap) | lowered, listed |
| `--pattern-max-predicate-contexts` | 100,000 | `max_predicate_contexts` (default = cap; increment 5b) | lowered, listed |
| `--pattern-max-predicate-work` | 100,000,000 | `max_predicate_work` (default = cap; increment 5b) | lowered, listed; at least 1 (refused at start-up below) |
| `--pattern-max-predicate-labels` | 10,000 | — (increment 5b: the names of a predicate's lists, a name in two lists counted twice; the server's policy, `caps.max_predicate_labels`) | 400 `predicate_too_large`; the flag is refused at start-up above 1,000,000 |
| `--traverse-chunk-target-ms` | 50 | — (not a cap: the annotation reads under the deadline are decoded in chunks of about this duration, as `/traverse`'s) | — |

The capabilities state every value in force (`caps`, `default_time_budget_ms`, `finalize_reserve_ms`; the two
delivery rates in `delivery_mbps`, since the owner's decision P9, §18, in `caps_rule`'s prose before), and every
answer echoes the effective ones (`limits`, §8.3), except
`max_checked_entries`, which no request field sets and only `caps` states (§18). The delivery rates
are starting estimates, conservative on the hosts measured (§7.6); an operator who raises
`--pattern-max-contexts` or `--pattern-max-patterns`, or serves on a slow or busy host, lowers them or raises
`--pattern-finalize-ms`. A clamp is never silent: `limits.clamped` lists it. The
service refuses values above its own ceilings at submit (`PROMPT-search-service-pattern.md` §3.1 item 5); the
server lowers them, so that a direct caller is not refused for a budget the host can give in part.

## 5. Order of the checks

A request is refused by the first check it fails, in this order:

1. multi-graph server: 400 `later_increment`, whatever the body;
2. the single index still loading: 503, no `code` (§6);
3. the body is not one RFC 8259 JSON text (§3: a comment, a trailing comma, anything after the value, a
   duplicated member name, nesting deeper than 1,000): 400 `invalid_request`;
4. the graph (`PatternSearch::support`): 400 with the graph's reason (`mask_invalid`, `alphabet_untested`,
   `representation_unsupported`, …, §6), whatever the body asks; a graph without its mask passes (§7.4);
5. the body is not an object: 400 `invalid_request`;
6. a later-increment field (§4.4, top level) or `in_ram`, in the alphabetical order of the body's field names;
7. `patterns` (presence, list, length), then each pattern in order (`id`, exactly one of `dna` / `iupac` /
   `protein`, its type, an unknown field); a pattern's alphabet is not a refusal (§8.9);
8. `mode`, `output` (`labels`, `occurrences`, `paths`, an unknown field), `scope`, `strands`,
   `stop_at_threshold`, `max_contexts`, `max_anchors`, `max_steps`, `time_budget_ms`, then increment 3's
   `max_labels_per_anchor`, `max_annotation_work`, `max_memory_mb`, `max_labels`, `max_occurrences_per_label`,
   `allow_unbudgeted_annotation`, then increment 4's `long_search`, `max_paths`, `require_support` (its value,
   then `"record_verified"` with `output.occurrences: false`), then increment 5's `genetic_code` (not an
   integer: 400 `invalid_request`; an integer that is no NCBI table: 400 `genetic_code_unknown`); the peptides are
   read in the genetic code after it (their alphabet, again, is no refusal); then increment 5b's `predicate` (its
   form, the path of the first fault named: 400 `invalid_request`; then its size: 400 `predicate_too_large`),
   `max_predicate_contexts`, `max_predicate_work`, `predicate_strands`, then `output.labels: "predicate_only"`
   without a predicate and a predicate with `long_search: "paths"` (both 400 `invalid_request`);
9. an unknown top-level field;
10. increment 3: `output.labels: "all"` in a retrieval mode on an annotation without the budget-aware decode,
    without `allow_unbudgeted_annotation: true`: 400 `annotation_unbudgeted`; increment 5b: also a request
    with a predicate, in every mode;
11. increment 4: `require_support: "record_verified"` with `long_search: "paths"` and `output.labels: "all"` in a
    retrieval mode, on an index whose best support (capabilities `support`) is not `record_verified`: 400
    `support_unavailable`, whatever the patterns' lengths.

## 6. Whole-request refusals

The body is `{"error": <message>, "code": <code>}`, except the 503 during loading, which has no `code`.

<!-- schema: refusal -->
| field | type | meaning |
|---|---|---|
| `error` | string | for people; may change |
| `code` | string | what to act on (table below) |

| status | `code` | when | what a client does |
|---|---|---|---|
| 400 | `invalid_request` | not JSON (§3: a comment, a trailing comma, content after the value, a duplicated member name, nesting deeper than 1,000), not an object, a wrong type or value, an unknown field, an empty or too long `patterns` list, `max_steps` < 1, `time_budget_ms` ≤ the reserve, `max_labels_per_anchor`, `max_annotation_work` or `max_memory_mb` < 1, `output.occurrences: true` without `output.labels: "all"`; increments 4 and 5: a `long_search` or `require_support` value not listed (§4.1), `require_support: "record_verified"` with `output.occurrences: false`, a `genetic_code` that is not an integer, none or more than one of `dna` / `iupac` / `protein`; increment 5b: a predicate that breaks a rule of §19.3 (the message names the path of the first fault, `request.predicate.and[1].none[0]`; a name that is a number names the fix, "write a taxid as \"562\""), `max_predicate_contexts` < 0, `max_predicate_work` < 1, a `predicate_strands` value not listed, `output.labels: "predicate_only"` without a predicate, a predicate with `long_search: "paths"`, a `long_search` value not served (`"supported_paths"` included) | fix the request |
| 400 | `later_increment` | a field of §4.4 (`graphs`, `budget_split`); a multi-graph server | wait for the increment the capabilities will announce |
| 400 | `predicate_too_large` | increment 5b: the predicate's lists name more than `caps.max_predicate_labels` names (10,000 by default; a name in two lists counted twice); the message names the count and the cap (§19.3) | send at most `caps.max_predicate_labels` names; split the cohort |
| 400 | `support_unavailable` | increment 4: `require_support: "record_verified"` with `long_search: "paths"` and `output.labels: "all"` in a retrieval mode, on an index that cannot verify a path in one record: not BASIC, no coordinates, or no record mapping (no `.seqs`, or `--no-coord-mapping`); capabilities `support` is then not `record_verified`. The message names the index's best support and placement | ask without `require_support`: each label of a path then states its support (`label_intersection` there) |
| 400 | `genetic_code_unknown` | increment 5: `genetic_code` is an integer that is not an NCBI translation table id (1–6, 9–16, 21–33; 7 and 8 were merged into 4 and 1, 17–20 are unassigned); the message names the ids | send one of the capabilities' `genetic_codes`, or omit it (1, the standard code) |
| 400 | `resident_only` | `in_ram`, any value: the route never loads an index inside a request (design §5.3) | drop `in_ram` |
| 400 | `mask_required` | **retired** (the owner's decision #16 of 2026-10-08): no build since answers it. A graph loaded without its dummy-edge mask (`.edgemask`) is answered, its counts `exact` where provable and otherwise `bounds` with an `estimate` (`counting: "upper_bound"`, §7.4); builds of version 1 before the decision refused such a graph with this code | a client of version 1 keeps handling it (an older host): its operator creates the mask (`metagraph transform --mask-dummy` once, then a restart; or `--pattern-build-mask`), or updates the build |
| 400 | `mask_invalid` | the graph's `.edgemask` marks valid a dummy edge whose last symbol (W) is `$`, as `metagraph extend` of earlier builds wrote it on a masked graph, or a stale mask left beside a rebuilt graph; counts on such a graph could be overstated, `exact` included (the owner's decision of 2026-10-07). Checked once when the graph is loaded (the start-up log names the edges found); a mask built at load (`--pattern-build-mask`) is not checked. A graph without a mask is not checked and not refused (§7.4) | the host's operator masks the graph again (`metagraph transform --mask-dummy --force`), then restarts the server; the capabilities say `available: false`, `unavailable_reason: "mask_invalid"` meanwhile |
| 400 | `representation_unsupported` | not a succinct graph (nor a PRIMARY one wrapped in `CanonicalDBG`), or k < 2 | none: this host has no pattern search |
| 400 | `primary_unwrapped` | a PRIMARY graph not wrapped in `CanonicalDBG` (the server always wraps; CLI or embedding misuse) | none |
| 400 | `alphabet_untested` | the graph's alphabet is `$ACGTN` (a DNA5 build): no DNA5 build has passed the pattern tests yet (the owner's decision of 2026-10-07; §8.2), with its mask or without | none on this build; a later build that passes them serves it |
| 400 | `alphabet_unsupported` | the graph's alphabet is neither `$ACGT` nor `$ACGTN` | none |
| 400 | `annotation_unbudgeted` | increment 3: `output.labels: "all"` in a retrieval mode, and the annotation has no budget-aware decode (capabilities `annotation: "unbudgeted"`: a column, BRWT, row or disk annotation; only the row-diff family has one); increment 5b: a request with a predicate on such an annotation, in every mode (its selection reads rows) | set `allow_unbudgeted_annotation: true` to read it without a memory bound on the reads (§14.4; a selection then reads single cells for at most 16 labels where the annotation has direct access, `selection.access: "columns"`), or ask for `labels: "none"` or mode `count` (without a predicate) |
| 503 | `deadline` | the answer could not be written by `time_budget_ms` (§7.6): nothing partial is sent | narrow the request: fewer patterns, a smaller `max_contexts` or `max_steps`; a larger budget helps only when it lets the work end early. The operator can raise `--pattern-finalize-ms` or lower the delivery rates (§4.5) |
| 503 | (none) | the single index is still loading (`Retry-After: 60`); every route answers so | retry later |

An unexpected failure (a server bug) is answered as every route answers it: 400 `{"error": …}` without a
`code` (500 `{"error": "Internal server error"}` for a non-standard exception). One such failure is stated: a
release that disagrees with its `exact` count (§7.5; in `all_or_count`, and in `partial` when the release ends short
of both the count and `max_contexts`) fails the request this way rather than be answered as complete or as cut; it has been seen only with a mask written by `metagraph extend` on a masked graph, which
marked W = `$` edges valid and is now refused (`mask_invalid`).

**What a client does with a failure.** Every refusal of this route is a JSON body with a `code`; the bodies
without one are the 503 while the index loads (`{"error": …}`, `Retry-After: 60`), the 400 of an unexpected
failure and the 500 above. Any other status, or a body that is none of these and not a JSON answer of §8 (a
proxy's page, a truncated body), is a failure of the backend: never read as an empty result. A client retries
the loading 503 after `Retry-After`; it does not retry a 503 `deadline`, nor a `withheld` or stopped answer,
with a larger budget on its own (the caller decides, §7.5); it passes an unknown `code` through verbatim (§1).

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
| `long` | implied by L > k | anchors (§7.7); with `long_search: "paths"` also its paths (§12.1) | `long` |

- `suffix` is the cheap, exactly countable scope; `any_offset` discovers the flank ranges, branching work even for
  an exact pattern, charged to `max_steps` (design §4.1).
- `suffix` is not served on a wrapped PRIMARY graph (a virtual suffix is a stored prefix, design §4.1): such a
  pattern is answered in its slot with `scope_unsupported`; `any_offset` is complete there.
- An exact 0 count says that no k-mer of this index contains the oriented pattern(s) in that scope, over its
  retained islands, with the mask or without it (an empty block is `exact` 0 on a graph without its mask too,
  §7.4). It says nothing about a label, a sample or a record as deposited (design §3, §5.1).

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
| `at_least` | a stop interrupted discovery; an undiscovered branch, offset or strand has no upper bound | the sum of the lower bounds of what was explored (0 for a search the stop met at its first step) |
| `bounds` | every range of every branch, offset and orientation was discovered, and either only the deferred scans (§7.6) were interrupted, or (on a graph without its mask, `counting: "upper_bound"`, below) some candidates were not checked for source dummies (more of them than `caps.max_checked_entries`); `lower` ≤ true ≤ `upper` | `lower` |
| `unknown` | the phase never ran: after a stop, or not in this increment | `null` |

- `bounds` carries `lower` and `upper`; no other relation does. `bounds` stays `bounds` when `lower` = `upper`.
- Sums (a total over offsets or strands; the same algebra merges counts of one unit across the calls a client
  splits a search into) follow `pattern_search.hpp`, `Count::operator+=`: `unknown` + `unknown` = `unknown`;
  `unknown` + anything else = `at_least`, and `at_least` + anything = `at_least`, whose value is the sum of the
  known lower bounds (an `unknown` adds 0, a `bounds` its `lower`); otherwise `bounds` + `bounds` or `exact` =
  `bounds`, `lower` and `upper` summed (an `exact` adds its value to both); `exact` + `exact` = `exact`. So
  `exact` 7 + `unknown` = `at_least` 7, and `bounds` [3, 5] + `exact` 7 = `bounds` [10, 12]. Not "the weakest
  relation wins": `exact` + `unknown` is `at_least`, not `unknown`.
- A search whose discovery was entered is `at_least` even when the stop refused its very first step (`at_least`
  0); one the stop came before is `unknown`.
- After a stop in discovery every per-offset count of the interrupted search is `at_least`: version 1 does not
  track which offsets were completely discovered before the stop. The mixed answers that do occur: a per-strand
  (or per-orientation) count `exact` beside an `at_least` total, and some per-offset counts `exact` beside a
  `bounds` total (a stop in a deferred scan).
- Counts are computed by the discovery and the deferred scans, before any release, and never from it: under
  `at_least` or `bounds`, `returned` (and the results at one offset) may exceed `counts.contexts.value` (or that
  offset's). The exception is a graph without its mask (below): there a release that enumerated every candidate
  makes the counts `exact`, and a cut one raises the lower bounds to what it released.
- A count is never promoted by assumption: anchors say nothing about paths, contexts nothing about occurrences.
- Units: `graph_contexts`, `anchors`, `paths`, `placed_occurrences`, `labels`. `labels` and `placed_occurrences`
  are `unknown` unless `output.labels: "all"` read them in a retrieval mode (§14.6).

**A graph without its dummy-edge mask** (`counting: "upper_bound"` in the capabilities and in the answer's
`index`; the owner's decision #16 of 2026-10-08, §18). The mask says which BOSS entries are dummies (§2). Without
it the search discovers the same ranges with the same steps, but cannot tell a source dummy from a k-mer without
spelling it, so it counts candidates:
- **U, the upper bound**, counts every candidate entry of a count's ranges: a flank range's entries whose W is
  not `$` (a sink dummy never counts), a last position's entries with W the allowed base (plain or marked). The
  source dummies among them are included: on BASIC, CANONICAL and odd-k PRIMARY graphs U is the count a masked
  graph gives plus the source dummies that hold the pattern after their `$` run, or, for a pattern with a
  leading N run (skipped, §7.8: for `reverse` the pattern's trailing run, which leads its reverse-complement
  window), whose `$` run ends inside that leading N run. Those last dummies do not hold the pattern (N never
  matches `$`), so U minus the masked count is not the number of dummies that hold a pattern starting with N
  (`NC` at offset 1 of a k = 9 graph: masked `exact` 84; 3 dummies hold `NC` after their `$` run and 3 more,
  such as `$$CCACACA`, have their `$` under the N; U = 90).
  The skipped part is the whole leading run, except in a window of N only: its last position is searched. On
  `$ACGTN` graphs nothing is skipped. On an even-k wrapped PRIMARY graph the palindrome scans settle some
  candidates (below), and U can be lower than that sum.
- **lower** counts the ranges whose nodes the search spelled whole (all k − 1 node symbols matched against
  pattern or flank bases, so no `$` fits): those are k-mers. In practice the offset-0 contexts of a pattern of
  L ≤ k without a leading N run (its `by_offset` "0" is `exact`; a leading N run is skipped, §7.8, so its ranges
  are not spelled whole), a pattern of length k in `suffix` scope, and the anchors of a pattern longer than k
  whose searched anchor windows start with no N run (its anchors are `exact`, and its extension runs as on a
  masked graph). On an even-k wrapped PRIMARY graph the palindrome scans spell every candidate they check, which
  settles those too.
- **Few unchecked candidates are checked** (the owner's decision #24 of 2026-10-08). When discovery and the
  deferred scans completed without a stop and the pattern's unchecked candidates — the entries counted into U
  and not into `lower`, over all its orientations and offsets: on a BASIC or CANONICAL graph its total's
  `upper` − `lower`; on a wrapped PRIMARY graph an entry can enter both orientations' counts, so there they can
  be as few as half of it — number at most `caps.max_checked_entries` (the server's
  `--pattern-max-checked-entries`, 50 by default, at most 1,000; 0 checks none), each of them is tested at
  query time (its node walked back for a `$`, at most k − 1 symbols): the k-mers are counted, the source
  dummies dropped, and **every count of the pattern is `exact`** — the total, `suffix`, `by_offset`,
  `by_strand` / `by_orientation`, the anchors of a pattern longer than k — the count a masked graph gives;
  `exact` 0 when all of them were dummies (an absence claim, §7.2). Each candidate tested is charged k − 1
  steps to `max_steps` (`work.steps`; `work.mask_scans` counts the ranges checked), so the check costs at most
  `max_checked_entries` × (k − 1) steps, about 1,500 at k = 31. It runs after the deferred scans, in every mode,
  before any release, so a count and a retrieval of the same request state the same counts. A stop during it
  (`max_steps` or `time`, phase `mask_scan`) leaves every count as discovery left it, `bounds`, stated as for
  any stop (§7.5: `discovery_budget` or `deadline` in `all_or_count`, the `cut` in `partial`). A pattern with
  more unchecked candidates is answered as if the check did not exist, field for field: no candidate is
  sampled or tested (the owner: per-query checking of large blocks costs too much). Examples on the mini
  (fixture `unmasked_checked`): blaNDM-1's forward primer, `bounds` [2, 24] without the check, 22 candidates
  tested, `exact` 24; an island start, [2, 32], 30 tested, `exact` 17; its first 14 bases, [4, 72], 68
  unchecked, `bounds` with its estimate; the island start in `suffix` scope on the forward strand, one
  candidate, a dummy: `exact` 0 (`unmasked_checked_dummies`).
- A count is `exact` when nothing of it is unchecked, when U = 0 (an empty block: `exact` 0, the
  absence claim of §7.2 holds), after the check above, or after a release that enumerated every candidate
  (§7.5); otherwise `bounds` {`lower`, `upper`: U}. A stop leaves `at_least` and `unknown` as on a masked
  graph; an `at_least` value is a true lower bound (the lower parts only). The algebra above sums them alike.
- **`estimate`**: every count with relation `bounds` of a graph without its mask (`counts.contexts` with its
  `suffix`, `by_offset` and `by_strand` / `by_orientation` parts, `counts.anchors` and its parts, `counts.paths`
  and its parts) carries `estimate` = round(U × f), kept inside [`lower`, U], f the graph's dummy fraction
  (`index.dummy_fraction.value`, §8.2): what the count would be if source dummies were as frequent among its
  candidates as among all the graph's entries. It is **not a bound** and never has a relation: a pattern at the
  start of a record, which the source dummies before that start hold, can have an estimate far above its true
  count (fixture `unmasked_count`: an island start of the mini, `bounds` [2, 32], `estimate` 32, `exact` 17 with
  the mask). Each count's estimate is computed from its own bounds: a total's estimate need not be the sum of its
  parts'. An `exact`, `at_least` or `unknown` count never carries one, nor does any count of a masked graph. The
  entry then carries the note `estimate_sampled_dummy_fraction` (§8.10).
- The **thresholds** compare U (§7.5): conservative, so that nothing whose count may pass a threshold is
  admitted; the note
  `threshold_upper_bound` says when such a decision went against the request while the lower bound was within
  the threshold. A count the check made `exact` is compared as on a masked graph (`all_or_count`'s threshold,
  the extension's admission). `stop_at_threshold` still stops discovery on the running U (the check comes after
  discovery), so it can stop a pattern the check would have made `exact` within its threshold.
- The **lists stay exact**: the release tests each unchecked candidate (at most k − 1 symbols read, under the
  deadline, not charged as steps) and never releases a source dummy; every released context and path is made
  of real k-mers, the list a masked graph releases.
- Mode `count` never releases contexts, so its counts of a pattern of L ≤ k with more unchecked candidates than
  `max_checked_entries` stay `bounds` where a retrieval request with the same body states `exact` after a
  complete release; the steps charged are the same. (With at most that many, the check makes both `exact`.
  The listing of a long pattern's anchors before its extension, `long_search: "paths"`, runs in every mode.)

### 7.5 Modes, `withheld` and `cut`

| mode | results | `retrieval_complete` |
|---|---|---|
| `count` | none: the entry has no `withheld`, `returned`, `cut` or `results` | always `false` (a count returns no context) |
| `all_or_count` | every context, only when discovery completed with an `exact` total ≤ `max_contexts` — or, on a graph without its mask, a `bounds` total whose upper bound ≤ `max_contexts`: the release then enumerates every candidate, drops the source dummies and makes the counts `exact` (§7.4) — and the contexts were handed to the route before the work time passed (§7.6); otherwise none, `withheld` says why | `true` iff the answer proves every context was returned: the count `exact` and equal to `returned` |
| `partial` | the first `max_contexts` contexts in answer order (§7.9) among those discovered, also after a step or threshold stop; `cut` says why the list may be shorter than the pattern's contexts. On a graph without its mask a release that drained every candidate of a completed discovery makes the counts `exact`, and a cut one raises each `lower` to the contexts it released at that orientation and offset (§7.4) | `true` iff the answer proves every context was returned (an `exact` count equal to `returned`); after any stop of the graph search (phases `discovery`, `mask_scan`, `extraction`) it is `false` and `cut` is stated, even when the list happens to hold every context |

These are the rules of the graph release, and all of them with `output.labels: "none"`. With `output.labels:
"all"`, `retrieval_complete` also needs every label of every returned context (§14.6). An incompleteness of the
labels alone — a row truncated or refused, a stop of the annotation reads or of the output of the labels
(phases `label_discovery`, `placement`, `output`), `labels_cut`, `occurrences_cut` — sets
`retrieval_complete: false` in `partial` without a `cut` (unless the memory account also shortened the list,
`cut: max_memory`): `labels_status`, `rows_refused`, `anchors_truncated`, `labels_cut`, `occurrences_cut` and
`stop` say what is missing. `all_or_count` withholds instead (§14.6).

`withheld.reason` (results absent; `returned: 0`, `results: []`):

| reason | when | the count | what to change |
|---|---|---|---|
| `count_above_threshold` | `all_or_count`: discovery completed, `exact` total > `max_contexts`; for a path search (§12.1), the paths `exact` and more than `max_paths`. On a graph without its mask also a `bounds` total whose upper bound U > `max_contexts` (the owner's decision #16: conservative; the note `threshold_upper_bound` when its `lower` ≤ `max_contexts`, so that the true count may fit: fixture `unmasked_threshold_upper_bound`) | `exact`; without the mask `bounds` | narrow the pattern, scope or strand; or `partial` (it lists the contexts whatever U is); on a graph without its mask, also a larger `max_contexts`, up to U |
| `threshold_crossed` | `all_or_count`, `stop_at_threshold`, L ≤ k: discovery stopped once its running lower bound passed `max_contexts` (on a graph without its mask, its running upper bound: the stop can come while the true count fits, stated by the note `threshold_upper_bound` when the lower bound had not passed it). The bound lags the count by the masked edges not yet scanned and, on an even-k wrapped PRIMARY graph, by the palindromic k-mers both base searches may find (without the mask: the running upper bound leaves out the ranges whose palindrome scan is pending); the deferred scans do not consult the threshold. So the stop can come late or not at all, and a pattern above its threshold can end `exact` with `count_above_threshold` (the owner's decision of 2026-10-07: the threshold is checked in discovery only) | `at_least` | as above |
| `discovery_budget` | `all_or_count`, L ≤ k: `max_steps` reached in discovery or a deferred scan | `at_least` or `bounds` | shorten or split an N run inside the pattern, add specified bases before it, or restrict `strands` to the orientation in which more specified bases precede it: the cost is set by where N runs sit, not by the bits (§7.8); the scope hardly changes it; a filter does not help |
| `deadline` | `all_or_count`, L ≤ k: the work time passed in discovery or a deferred scan, or in the release or while its contexts were handed to the route (all or nothing: a deadline during the hand-over withholds all of them, the counts kept) | as stopped | a larger `time_budget_ms`, or as for `discovery_budget` |
| `paths_later_increment` | either retrieval mode, L > k without `long_search: "paths"`, unless the anchors are `exact` 0 (§7.7); whatever stopped the anchors is in `stop` | the anchors' | ask with `long_search: "paths"` (increment 4, §12.1); a request that does not keeps this answer |
| `anchors_above_threshold` | increment 4, `long_search: "paths"`, L > k, either retrieval mode (`partial` too): the anchors `exact` and more than `max_anchors` (on a graph without its mask also `bounds` with U > `max_anchors`, the note `threshold_upper_bound` when their `lower` ≤ `max_anchors`), so the extension was not admitted (§12.1) | the anchors' (`exact`, or without the mask `bounds`); `counts.paths` `unknown`, `extension: "not_admitted"` | raise `max_anchors`, or narrow the pattern's anchor window (its first k bases, and its last k with `strands` `both` or `reverse`) |
| `annotation_budget` | increment 3, `all_or_count`, `labels: "all"`: a row the memory account refused (`rows_refused`), or the reads stopped at `max_annotation_work` | `exact` | raise `max_memory_mb` or `max_annotation_work`, narrow the pattern, or `partial` |
| `anchor_labels_truncated` | increment 3, `all_or_count`, `labels: "all"`: a row carried more labels than `max_labels_per_anchor` (`anchors_truncated` lists each, with its total) | `exact` | raise `max_labels_per_anchor` to the largest total, or `partial` |
| `output_budget` | increment 3, `all_or_count`, `labels: "all"`: the memory account could not hold the answer (the contexts' descriptors or the labels and occurrences built for them); increment 5b: also a selected context's result object, with any projection | `exact` | raise `max_memory_mb`, narrow the pattern, or `partial` |
| `predicate_above_threshold` | increment 5b, `all_or_count` with a predicate: the raw count `exact` (without the mask: its upper bound) above `max_predicate_contexts`; nothing was read (§19.6, §19.8) | the raw count's | narrow the pattern, scope or strands, raise `max_predicate_contexts` up to its cap, or `partial` (it tests the first `max_predicate_contexts`) |
| `selected_above_threshold` | increment 5b, `all_or_count` with a predicate: the selected count `exact` and above `max_contexts` | the raw count's; `selected` `exact` | narrow the predicate or the pattern, raise `max_contexts`, or `partial` |
| `predicate_budget` | increment 5b, `all_or_count` with a predicate: the selection stopped at `max_predicate_work` (also an earlier pattern's: sticky), or the memory account could not hold it (its descriptors, its rows, a row it had to read: `rows_refused`, phase `selection`), or the predicate itself (§19.9) | the raw count's; `selected` `bounds` or `at_least` | raise `max_predicate_work` or `max_memory_mb`, narrow the pattern |

With `labels: "all"`, `deadline` also names a time stop of the annotation reads or of the output of the labels
(`stop.phase` `label_discovery`, `placement` or `output`); the graph count is then still the discovery's. With
a predicate (§19.8) `deadline` also names a time stop of the selection (`stop.phase` `selection`) or of the
results built for the selected contexts (`output`), `threshold_crossed` also a `stop_at_threshold` stop of the
raw discovery at `max_predicate_contexts` or of the selection at `max_contexts`, and `cut` reasons `max_memory`
and `time` also the selection's.

An `all_or_count` release that disagrees with its `exact` count is never published: the request fails as a
server bug (§6).

`cut.reason` (`partial`, neither complete nor withheld; `results` holds what was released):

| reason | when |
|---|---|
| `max_contexts` | discovery completed (or stopped at the threshold) with more contexts than `max_contexts`: the first `max_contexts` in answer order |
| `max_steps` | discovery or a deferred scan stopped at `max_steps`: the first `max_contexts` of those discovered |
| `time` | the work time passed: in discovery nothing is released (its membership would depend on the machine, `returned: 0`); in the release, what was released before (the clock is read every 64 contexts handed to the route, §7.6). Also the cut of a pattern already stopped by `max_steps` or the threshold whose release met the work time (§7.6: `stop` keeps the first stop) |
| `max_memory` | increment 3, `labels: "all"`: the memory account held the descriptors of only the first `returned` contexts (`partial`'s descriptors take at most half of the account, §14.4; the list's length; it replaces the engine's `max_contexts` cut when both cut) |
| `max_paths` | increment 4, `long_search: "paths"`: more paths than `max_paths` (completed, or stopped at the threshold in the extension): the first `max_paths` in answer order (§12.1) |
| `max_predicate_contexts` | increment 5b, with a predicate: the raw release held the first `max_predicate_contexts` contexts in answer order (or stopped at that threshold): the selected among them (§19.8) |
| `max_predicate_work` | increment 5b, with a predicate: the selection stopped at its work budget: the selected among the contexts it decided (§19.8) |

`max_anchors` is a value of `cut.reason` reserved for a later increment's release of anchors; version 1 never
releases anchors (a path search cut before any extension, by `stop_at_threshold` on the anchors, says
`max_anchors` with nothing returned, §12.1).

### 7.6 Budget, deadline and the finalisation reserve

- **One budget per request** (design §5.3): `max_steps` and the deadline are spent by the patterns in request
  order. A stop is sticky: every later pattern answers `unknown` counts with the same `stop` (phase
  `discovery`), `work` zero, and in a retrieval mode the matching `withheld` or `cut` (in `partial`, also after
  a `max_steps` stop, a later pattern says `cut: time` and `determinism: time_limited` when the work time has
  passed before its empty release). `stop_at_threshold` is not a budget stop: it ends only its own pattern.
  One exception, stated (§12.2, §18): a peptide without instances (its `*` read in a table without a stop codon,
  note `no_stop_codon`) of L ≤ k is answered before the budget is read, as an error slot is: `exact` 0 in every
  count, `stop` `null`, `work` zero, `determinism: "full"`, wherever it sits (fixture
  `peptide_no_stop_codon_after_stop`). Its answer needs no search, so it depends on no budget or clock.
- **The deadline** starts when the request's body is parsed (time spent in the server's queue is not in it).
  Work stops at `time_budget_ms − finalize_reserve_ms − E`, E being the estimated time to write the answer
  built so far: E = 1.25 × (B / (b × 1000) + B / (c × 1000)) ms for B bytes of compact JSON text of the results
  built (with `labels: "all"` the label objects about to be built are counted in B, and once more at the rate
  b, before they are built), b and c the delivery rates in MB/s (`--pattern-delivery-build-mbps` 10,
  `--pattern-delivery-compress-mbps` 50; the capabilities' `delivery_mbps` states the values in force). E is read with
  the work time at every reading, so it grows as results are built: about 0.15 ms per KB of results at the
  defaults, 2.8 s for 16 × 10,000 results (18.4 MB of text, which an M-series Mac wrote and gzipped in 0.3–0.45 s:
  the rates are conservative starting estimates, not measurements of the host). The answer — the counts kept up
  to date during the search and every result released — is then written in the time left, so that a request
  stopped by time still answers with its counts. `finalize_reserve_ms` is the floor that covers what E does not
  model: the work done past the work time until the next reading, assembling the counts, and the transport.
- **What `time_budget_ms` is, and is not.** It is the backend's deadline for preparing the answer, from the
  parsed body to the hand-over to the transport. It is not a search time: the work stops `finalize_reserve_ms`
  plus E before it, so a request can stop by time well before its budget (E alone is 2.8 s for 16 × 10,000
  buffered results at the default rates; a stop 1.3 s before a 2 s budget is within the rule). It is not the
  caller's latency: the wait for a request-pool
  thread (§3) and the transport are outside it. A finalisation that overruns it is still a 503 `deadline`
  (below). Three clocks are distinct: this backend deadline, the client's HTTP timeout (above it: the content
  timeout is 900 s), and the lifetime of a job that wraps the call (the service's).
- **Where the work time is read** (and, on the server, whether the client has left, §3):
  - at the start of every pattern, and before every release;
  - at every multiple of 4,096 steps charged (discovery, the deferred scans and, for L > k, the extension share
    one step count, so the boundary from discovery to the deferred scans has no reading of its own);
  - in the release: every 4,096 range descriptors it prepares and every 4,096 edges it examines;
  - every 64 contexts handed to the route — in `partial`'s release and in `all_or_count`'s delivery of its
    buffered release — since each costs the route a k-mer spelling (k − 1 graph steps) and its result object;
  - with `labels: "all"`: before each annotation read and between its chunks (§14.4), and before the labels of
    each context are built for the answer;
  - since the review GPT-3 (§18): for L > k with `long_search: "paths"`, before each anchor's extension (its
    spelling and its depth-first search), every 64 anchors listed and before every 64th node the search
    expands; after a completed exact pattern's search, before every piece but the first of its low-complexity
    diagnostic (§7.8); with `labels: "all"`, in the work between the reads — the occurrences of each context or
    path (made in pieces of at most 4,096), the paths' label lists and their verification — at least every 4,096
    units of it, before the work (§12.1, §14.4).

  Work therefore ends within one such stride after the work time. Not read: the O(L) parsing, bits and palindrome
  test of each pattern's text (§4.2), and the spelling of one anchor (k − 1 BOSS steps).
- **503 `deadline`.** The answer is assembled, serialised and compressed under the same deadline, read every
  4,096 entries and results while it is assembled, every 64 KiB of text while it is written, between
  compression blocks, and once more before it is handed to the transport. If it cannot be written by
  `time_budget_ms`, the answer is 503 `deadline` and nothing partial is sent.
- `stop` names **the first stop that touched the pattern**: `{phase, reason}` with phase `discovery` (the range
  search), `mask_scan` (the deferred scans: a range's masked edges, on an even-k wrapped PRIMARY graph the
  palindrome check of every context at an offset where a palindromic k-mer can hold the pattern, and on a graph
  without its mask the check of a pattern's few unchecked candidates, §7.4; the only phase whose stop can leave
  `bounds`) or `extraction` (the release, and `all_or_count`'s delivery of it to the
  route), and reason `max_steps`, `time`, `max_contexts` or `max_anchors` (the last two: `stop_at_threshold`).
  Increment 4 (`long_search: "paths"`) adds the phase `extension` (the depth-first extension of the anchors,
  §12.1) and the reason `max_paths` (`stop_at_threshold` in the extension).
  A later stop does not replace it: in `partial`, a pattern stopped by `max_steps` or by its threshold is still
  released, and a time stop in that release shows only as `cut: time` and `determinism: time_limited`, `stop`
  keeping the earlier reason (the owner's decision of 2026-10-07). Increment 3 (`labels: "all"`) adds the phases
  `label_discovery` (the first read, §14.2), `placement` (the second) and `output` (the descriptors and labels
  built for the answer), and the reasons `max_annotation_work` and `max_memory` (§14.4); `output` is stopped by
  `max_memory` or by `time` (§14.4); the engine's stop, when there is one, is the one stated. Every assignment of
  `stop` in this document applies only while `stop` is still `null` (first stop wins, the owner's decision of
  2026-10-07). Increment 5b (a predicate, §19.8) adds the phase `selection` (the selection pass: its
  reverse-complement lookups, its row reads and its decisions) with the reasons `max_predicate_work`,
  `max_memory`, `time` and `max_contexts` (`stop_at_threshold` on the selected count), and the reason
  `max_predicate_contexts` of the phase `discovery` (`stop_at_threshold` on the raw count with a predicate); the
  results built for the selected contexts stop in the phase `output` (`time`, `max_memory`). So a later stop shows only in what it left: `labels_status: "output_budget"` (a memory or a time
  stop of the output) or `not_read`, `cut: time`, `withheld`, while `stop` names an earlier phase. Any time stop,
  stated in `stop` or not, sets `determinism: "time_limited"` (§7.9); so does a low-complexity diagnostic the work
  time cut (§7.8), with `stop` `null`.
- `work` per pattern: `ranges_visited` (range evaluations), `mask_scans` (ranges whose deferred scan began),
  `steps` (every step charged: `ranges_visited` plus the items the deferred scans examined, plus k − 1 per
  candidate the check of §7.4 tested, plus, with `long_search: "paths"`, the outgoing edges the extension
  examined, stated as `extension_edges`). The `steps` of
  all patterns sum to at most `max_steps`. On an even-k wrapped PRIMARY graph the deferred scans check every context
  at a palindrome-capable offset, one k-mer spelling each, so an `any_offset` count there costs time and steps
  linear in those contexts (`index.graph_mode` and `index.k` tell a client so).
- **Memory.** The label-free path (mode `count`, or `output.labels: "none"`) has no memory account in version
  1; what a request holds is bounded by the caps:
  - the DFS frontier, O(k × alphabet) per base search;
  - the ranges whose count awaits a deferred scan, 88 bytes each: none in the common case; on an even-k wrapped
    PRIMARY graph one per range at a palindrome-capable offset, so up to `max_steps` of them for a short or
    degenerate pattern;
  - `all_or_count`: a 24-byte descriptor per range with contexts while the running lower bound (without the
    mask: the running upper bound) is ≤ `max_contexts` (about `max_contexts` of them), none once it is above;
  - `partial`: only the descriptors that can hold one of the first `max_contexts` contexts in answer order —
    about `max_contexts`, plus the ranges straddling the cut-off (O(k) per base search) and those whose count
    awaits a scan; between two compactions at most 4,096 or twice the number kept at the last one; none with
    `max_contexts` 0. The release indexes them with 40 bytes each. On a graph without its mask the retention
    bound counts only the contexts it is sure of, which a range not spelled whole is not: `partial` then keeps
    every discovered range not spelled whole (24 bytes each, at most one per step charged), not about
    `max_contexts` of them (a stated limitation of the owner's decision #16); `all_or_count` and the extension
    keep theirs while the running lower bound is ≤ the threshold as long as the unchecked candidates are at most
    `max_checked_entries` (the check may still make the count `exact` within it), and every mode keeps the
    unchecked ranges for the check while they are that few (at most `max_checked_entries` of them);
  - the results built for the answer: at most `max_contexts` per pattern.

  With `output.labels: "all"`, `max_memory_mb` bounds the labelled retrieval (§14.4). With a predicate the
  account exists in every mode and bounds the selection too (§19.9).

### 7.7 Patterns longer than k

- Scope `long`, whatever was requested. Each searched orientation anchors on its own window: `forward` (and a
  palindrome) on P[0, k), `reverse` on rc(P)[0, k), whose bits are those of P[L − k, L).
  `anchor_information_bits` states the bits of P[0, k), whatever the strands; `min_anchor_information_bits`
  states the bits of the least informative searched window — the lower of the two with `strands: "both"` — and
  the information floor gates on it (§7.8).
- `counts.anchors` counts the anchors (unit `anchors`, by strand or orientation); `counts.paths` is `unknown`,
  except `exact` 0 when the anchors are `exact` 0 (no path starts without an anchor: a derivation, not a
  promotion). `counts.contexts` is absent. On a graph without its mask the anchors are `exact` when every
  searched anchor window starts with no N run (spelled whole, no source dummy can hold one, §7.4), else `bounds`
  with an `estimate`; with `long_search: "paths"` and U ≤ `max_anchors` they are listed (the source dummies
  dropped, their count `exact`) and extended as on a masked graph (fixture `unmasked_paths`).
- Nothing is extracted: in a retrieval mode the results are withheld with `paths_later_increment`, unless the
  anchors are `exact` 0, in which case the empty answer is complete (`retrieval_complete: true`).
- The note `paths_later_increment` is on every such entry. `max_anchors` is the `stop_at_threshold` threshold.
- This is the answer of `long_search: "anchors"`, the default of the request field. Paths are opt-in: a request
  without `long_search: "paths"` is answered as here; with it, the anchors are extended into paths, counted and
  released as §12.1 states (increment 4).
- A peptide longer than k (more than k / 3 residues) is such a pattern: anchored on the first k bases of its
  codon automaton (a window that can cut a codon), extended through the automaton with `long_search: "paths"`
  (§12.2).

### 7.8 The information floor

- A pattern below `caps.min_information_bits` (24 by default, about 12 specified bases) is answered in its
  slot with `information_below_floor` and costs no step. For L > k every searched orientation's anchor window
  must reach the floor (§7.7; the least of them is `min_anchor_information_bits`): P + N^k is refused with
  `strands` `both` or `reverse` (the message names the reverse orientation's window), N^k + P is answered with
  `strands: "reverse"`.
- Exempt: an exact pattern (every position one base, whatever its kind; a peptide whose every residue has one
  codon in its genetic code, M and W in the standard code) in `suffix` scope, however short: one range, a few
  ranks.
- A peptide's bits are exact (§2, §12.2), also for an anchor window that cuts a codon: the floor gates peptides
  as it gates DNA.
- The floor is a planning heuristic, not a cost bound. A searched window is matched from its first position on
  the widest ranges, so the cost is set by where its N runs sit, whatever the pattern's bits: a run of N costs
  about min(4^run, edges / 4^a) ranges per level, a being the specified bases before it in that orientation's
  window. So a run at the end of a window is cheap, and a run at its start (for `forward` a leading run of P,
  for `reverse` a trailing one) would be the most expensive: on a `$ACGT` graph the engine skips it and searches
  the rest of the window with its offsets shifted, which gives the same contexts; on `$ACGTN` (a pattern's N
  never matches the graph's N) it is searched as given. A run inside the pattern stays expensive in both
  orientations when few specified bases precede it on either side; `discovery_budget` (§7.5) says what to
  change.

  A pattern's base searches run cheapest first by that estimate, so that a budget stop leaves the cheaper
  orientation complete and the other interrupted (§7.9).
- Low-complexity patterns are not refused; an exact one that sdust flags over its whole text carries the note
  `low_complexity_pattern` (for L > k the flag can come from bases outside the anchor windows) — when its answer
  states no stop. Since the review GPT-3 (§18) the diagnostic runs on a completed search only: not after a stop of
  the pattern's own nor after an earlier pattern's request-wide one, threshold stops included (a stopped answer's
  counts say what they are; the note is a hint about large counts). sdust reads the text in pieces of 128 bases,
  each overlapping the next by 63 (its window less one), and stops at the first piece it flags; the flag is the
  one over the whole text, since sdust's state at a base depends only on the 64 bases ending there. The clock is
  read before every piece but the first, so a pattern of at most 191 bases (a peptide of at most 63 residues) is
  always diagnosed; for a longer one, a work time that passes between two pieces leaves the note out and states
  `determinism: "time_limited"` with `stop` `null`, the counts complete (§7.9).

### 7.9 Ordering and determinism

- `patterns`: request order, one entry per request pattern.
- `results`: node ascending, then offset ascending, then orientation (`forward`/`+`, `reverse`/`-`,
  `palindromic`/`=`) (design §5.5, "BOSS edge order, then by offset"; the orientation separates the two
  contexts an IUPAC pattern can have at one offset). Paths (§12.1): anchor node ascending, then orientation,
  then the sequence (A < C < G < T).
- `strands` (entry): the orientations searched, forward before reverse; `["="]` or `["palindromic"]` for a
  palindrome. This is the order of the plan, not of the work: the base searches run cheapest first (§7.8), so
  after a budget stop the reverse orientation can be `exact` and the forward one `at_least`.
- `notes`: `low_complexity_pattern`, `strand_unknown_canonical`, `paths_later_increment`, then those of the
  owner's decisions of 2026-10-08 (§18), `threshold_upper_bound`, `no_stop_codon` and
  `estimate_sampled_dummy_fraction`, then increment 3's `annotation_unbudgeted`, `record_bounds_unknown` or
  `annotation_not_read`, then increment 4's `label_intersection_only`, then increment 5b's
  `predicate_constant` and `projection_not_read`, in that order.
- `limits.clamped`: `max_contexts`, `max_anchors`, `max_steps`, `time_budget_ms`, then increment 3's
  `max_labels_per_anchor`, `max_annotation_work`, `max_memory_mb`, `max_labels`, `max_occurrences_per_label`,
  then increment 4's `max_paths`, then increment 5b's `max_predicate_contexts`, `max_predicate_work`, in that
  order.
- Increment 5b: a result's `selection_labels` in the label order of the pattern's returned results (contexts
  desc, column asc), its `selection_strands` in the same order (one per label); the predicate's echo in the
  request's syntax, its names in the request's order.
- Increment 4: the labels of paths in `by_label` and in each path by (paths desc, column asc); a label's
  occurrences by (`seq_id`, start), or by `kmer_coord` without a record mapping.
- Increment 3: labels in `by_label` and in each result by (contexts desc, column asc) (design §5.5); a label's
  occurrences by (`seq_id`, start, strand), or (`kmer_coord`, `offset`, strand) without a record mapping.
- JSON objects (`by_offset`, `by_strand`, `by_orientation`, …) carry no order; the server writes keys sorted
  as strings (`"10"` before `"2"`).
- **Determinism.** The same request on the same index, answered by the same build under the same effective
  configuration (the caps and defaults in force, the floor, the delivery rates, `index.release`), gives the same
  answer, byte for byte apart from `timing`, including which contexts a cut kept, because every budget is spent
  in a fixed order. The exception is a time stop: the entry it touched (also when it shows only as `cut: time`
  or `labels_status: "output_budget"`, §7.6) and every answered entry after it state `determinism:
  "time_limited"`; their counts, `work` and released contexts depend on the machine. One `time_limited` entry
  states no stop (§7.8, §18): a completed exact pattern of more than 191 bases whose low-complexity diagnostic the
  work time cut — its counts, `work` and results are complete and deterministic, only the note's absence depends
  on the machine (a client must not read `time_limited` as "a stop is stated"). A 503 `deadline` depends on
  the machine too. `determinism: "full"` promises nothing across builds: another build of version 1 keeps the
  meaning of every field and count but may spend its budgets otherwise (§1), so its `work`, its stops and its
  `at_least` and `bounds` values can differ.

### 7.10 Node and row ids

- `node` is the id of the context's k-mer in the served graph: the BOSS edge index on a DBGSuccinct (BASIC,
  CANONICAL), the wrapper id on a wrapped PRIMARY graph (a stored k-mer's id, or that plus the wrapper's offset
  for a virtual reverse complement).
- `row` is the annotation row the k-mer is annotated in, named without reading it: `node − 1` on BASIC; on
  PRIMARY the stored k-mer's (`get_base_node(node) − 1`); on native CANONICAL the canonical k-mer's — in
  MetaGraph's sense: of the k-mer and its reverse complement the one with the smaller BOSS edge index, not the
  lexicographically smaller one — minus 1 (the annotation key of every route), so that a k-mer and its reverse
  complement share one row. `null` if the k-mer has no row in the annotation (never expected on a compatible
  index; stated rather than guessed).
- Both are **opaque** to clients and valid only per (host, `index.index_fp` or `index.release`, graph). A later
  increment takes row ids back ("labels for given rows", design §12); a client keeps them as given.

## 8. The answer (HTTP 200)

### 8.1 Top level

<!-- schema: answer -->
| field | type | meaning |
|---|---|---|
| `pattern_contract_version` | integer | 1 |
| `mode` | string | the request's mode, defaults applied |
| `output` | object \| null | `null` in mode `count`; in the retrieval modes the projection applied: `{"labels": "none"}`, or `{"labels": "all", "occurrences": <boolean>}` (increment 3), or `{"labels": "predicate_only", "occurrences": <boolean>}` (increment 5b) |
| `index` | object (§8.2) | what was searched |
| `limits` | object (§8.3) | the effective caps of this request |
| `timing` | object (§8.4) | the request's elapsed time |
| `patterns` | list of entries (§8.5) | one per request pattern, in request order |
| `predicate` | object (`predicate_block`, §19.10) | increment 5b, **exactly when the request names a predicate**: the predicate as bound to this index (its normal form, its names, the unknown ones, `vacuous`, the scope and strands of its claim) |

### 8.2 `index`

<!-- schema: index -->
| field | type | meaning |
|---|---|---|
| `index_ns` | string \| null | the server's `--index-name`, `null` without it |
| `index_fp` | string \| null | the index identity (`--index-manifest`'s digest), `null` without a manifest |
| `release` | string | the server's `--index-release` (`""` when not set) |
| `k` | integer | the graph's k |
| `graph_mode` | `"basic"` \| `"canonical"` \| `"primary"` | §7.3 |
| `alphabet` | `"$ACGT"` \| `"$ACGTN"` | the BOSS alphabet, sentinel first. This build answers `$ACGT` only: a `$ACGTN` (DNA5) graph is refused with `alphabet_untested` (§6; the owner's decision of 2026-10-07) until a DNA5 build passes the pattern tests. No DNA5 server or CLI has been built or run, and CI builds DNA and Protein only (the engine's unit suites, DNA5 branches included, passed once on a DNA5 build of the engine and its tests, by hand, 2026-10-07; not in CI). Known, for that later build: on a DNA5 wrapped PRIMARY graph with odd k, a stored k-mer equal to its reverse complement (N at its centre between complementary flanks, e.g. `ACNGT` at k = 5) is counted and released twice |
| `strand_stated` | boolean | `true` on BASIC: contexts carry `strand`, counts `by_strand` |
| `counting` | `"upper_bound"` | the owner's decision #16, **answers on a graph loaded without its dummy-edge mask only** (an answer on a masked graph has neither this field nor the next: its counting is `exact`, as the capabilities say): counts are `exact` where provable and otherwise `bounds` with an `estimate` (§7.4) |
| `dummy_fraction` | object (`dummy_fraction` below) | likewise: f, the dummy fraction the estimates rest on; the capabilities' `dummy_fraction` |

<!-- schema: dummy_fraction -->
| field | type | meaning |
|---|---|---|
| `value` | number | f = the real k-mers among the entries drawn / the entries drawn: the fraction of real k-mers among the graph's entries whose W is not `$` (the entries a pattern can count; sink dummies, W = `$`, are excluded exactly by ranks of W). 0.9999 on the mini index (exact f 0.999955: 8,335,373 of 8,335,747) |
| `interval` | [number, number] | its 95% interval (Wilson's score interval, clamped to [0, 1] and holding `value`): [0.999434, 0.999982] on the mini |
| `samples` | integer | the entries drawn: 10,000 |
| `source` | `"sampled"` | how f was obtained (an extensible enumeration, §1): drawn at load, 10,000 entries uniformly with replacement among the edges whose W is not `$`, with `std::mt19937_64` seeded by the graph's number of edges (so the same graph gives the same f in every process and after every restart); a drawn entry is real iff its source node holds no `$` (at most k − 1 symbols read). Sampled once per graph in the loading thread (40–60 ms on the mini, 188 ms on a synthetic graph of 3 × 10⁷ edges) |

### 8.3 `limits`

<!-- schema: limits -->
| field | type | meaning |
|---|---|---|
| `max_contexts` | integer | effective (§4.1) |
| `max_anchors` | integer | effective |
| `max_steps` | integer | effective, for the whole request |
| `time_budget_ms` | number | effective |
| `finalize_reserve_ms` | number | the reserve inside it: the floor of the time kept back from the work for the answer (§7.6) |
| `min_information_bits` | number | the floor |
| `max_patterns` | integer | the server's |
| `stop_at_threshold` | boolean | as requested |
| `max_labels_per_anchor` | integer | increment 3, answers with `labels: "all"` only: effective (§4.1) |
| `max_annotation_work` | integer | likewise |
| `max_memory_mb` | integer | likewise |
| `max_labels` | integer | likewise |
| `max_occurrences_per_label` | integer | likewise |
| `allow_unbudgeted_annotation` | boolean | likewise, as requested |
| `long_search` | `"paths"` | increment 4, answers to `long_search: "paths"` only: as requested |
| `max_paths` | integer | likewise: effective (§4.1) |
| `require_support` | string | likewise, and only when labels are read (`labels: "all"` in a retrieval mode): as requested, default applied |
| `max_predicate_contexts` | integer | increment 5b, answers with a predicate only: effective (§4.1). With a predicate the annotation limits above (`max_labels_per_anchor` … `allow_unbudgeted_annotation`) are echoed in every mode: the selection uses the account |
| `max_predicate_work` | integer | likewise: effective |
| `max_predicate_labels` | integer | likewise: the server's (`caps.max_predicate_labels`) |
| `predicate_strands` | string | likewise: as requested, default applied (`predicate.strands` states what was evaluated) |
| `clamped` | list of objects | one per request value lowered to its cap, in the order of §7.9; `[]` when none |

<!-- schema: clamped -->
| field | type | meaning |
|---|---|---|
| `field` | string | `max_contexts`, `max_anchors`, `max_steps` or `time_budget_ms`; increment 3: `max_labels_per_anchor`, `max_annotation_work`, `max_memory_mb`, `max_labels`, `max_occurrences_per_label`; increment 4: `max_paths`; increment 5b: `max_predicate_contexts`, `max_predicate_work` (named also in a request without a predicate) |
| `requested` | number | the value asked for |
| `effective` | number | the cap applied |

### 8.4 `timing`

<!-- schema: timing -->
| field | type | meaning |
|---|---|---|
| `elapsed_ms` | number | top level: from the deadline's start to the end of the search (the annotation reads included); in an entry: that pattern's search. Varies between runs |
| `label_discovery_ms` | number | increment 3, entries with `labels: "all"`: the first read (§14.2). Varies between runs |
| `placement_ms` | number | likewise, the second read (§14.3) |
| `extension_ms` | number | increment 4, entries of a pattern longer than k with `long_search: "paths"`: the extension (§12.1). Varies between runs |
| `label_intersection_ms` | number | review GPT-3 (§18), entries of a pattern longer than k with `long_search: "paths"` and `labels: "all"`: the intersection of the label lists of each path's rows (§12.1). Varies between runs |
| `verification_ms` | number | likewise: the verification of the labels carrying the paths — the chains' join, their record placement, the runs kept for the output (§12.1) — apart from `placement_ms`, the coordinates' reads; the loop's time also where nothing is verified (no coordinates). Varies between runs |
| `selection_ms` | number | increment 5b, every answered entry of a request with a predicate: the selection pass (its lookups, reads and decisions; 0 where it did not run). Varies between runs |

### 8.5 A pattern's entry

An entry is one of three shapes: answered, refused by the engine, or refused for its alphabet (§8.9).

<!-- schema: entry -->
| field | type | present | meaning |
|---|---|---|---|
| `id` | string \| null | always | the request's `id` |
| `kind` | `"dna"` \| `"iupac"` \| `"protein"` | always | which key the pattern came in (`protein` since increment 5) |
| `pattern` | string | answered, engine refusal | the pattern in upper case (a peptide's residues) |
| `length` | integer | answered, engine refusal | L, in bases whatever the kind (a peptide of m residues: 3m) |
| `residues` | integer | `protein` only, answered, engine refusal | increment 5: m, the peptide's residues (§12.2) |
| `genetic_code` | integer | `protein` only, answered, engine refusal | increment 5: the NCBI translation table it was read in (the request's `genetic_code`, default 1) |
| `information_bits` | number | answered, engine refusal | §2 |
| `anchor_information_bits` | number \| null | answered, engine refusal | L > k: the bits of the anchor window P[0, k), whatever the strands searched (§7.7); else `null` |
| `min_anchor_information_bits` | number \| null | answered, engine refusal | L > k: the bits of the least informative searched anchor window, the information floor's operand (§7.7, §7.8: P[0, k) for `forward` or a palindrome, rc(P)[0, k) — the bits of P[L − k, L) — for `reverse`, the lower of the two for `both`); else `null`. Added by the review of 2026-10-07 (§16) |
| `error` | object (§8.9) | refusals only | `{code, message}` |
| `mode` | string | answered | the request's mode |
| `scope` | `"suffix"` \| `"any_offset"` \| `"long"` | answered | §7.2 |
| `strands` | list of strings | answered | the orientations searched (§7.3, §7.9) |
| `palindromic` | boolean | answered | P = rc(P) |
| `counts` | object (§8.6) | answered | |
| `work` | object (§8.7) | answered | |
| `stop` | object \| null | answered | §7.6: the first stop that touched the pattern; `null` when the pattern's work completed |
| `retrieval_complete` | boolean | answered | §7.5; the one flag that licenses an absence claim from the lists: a graph context not in `results` is not in the scope, and with `output.labels: "all"` a label not in `by_label` carries no occurrence there (§9, §14.6). Without labels it licenses nothing about labels |
| `withheld` | object \| null | answered, retrieval modes | `{reason}` (§7.5) |
| `returned` | integer | answered, retrieval modes | the length of `results` |
| `cut` | object \| null | answered, retrieval modes | `{reason}` (§7.5); only in `partial` |
| `results` | list of results (§8.8) | answered, retrieval modes | the released contexts, in the order of §7.9; for a pattern longer than k with `long_search: "paths"`, the released paths (`path_result`, §12.1) |
| `absence_scope` | `"suffix_only"` \| `"any_offset"` \| `"long"` | answered | §7.2 |
| `determinism` | `"full"` \| `"time_limited"` | answered | §7.9 |
| `notes` | list of strings | answered | §8.10 |
| `timing` | object (§8.4) | answered | |
| `placement` | string | answered, `labels: "all"` (or `"predicate_only"`, increment 5b: the projection fields `placement` to `occurrences_cut` alike) | increment 3: what this answer places: `record`, `global`, `none`, `none_canonical` (§14.3), or `not_requested` (`output.occurrences: false`) |
| `annotation` | `"budgeted"` \| `"unbudgeted"` | answered, `labels: "all"` | the reads' access (§14.4) |
| `by_label` | list \| null | answered, `labels: "all"` | the per-label summary over the returned contexts (§14.5), or over the returned paths (`by_label_paths`, §12.1); `null` when the results are withheld or, in `partial`, when the memory account could not hold it (§14.4) |
| `rows_refused` | list | answered, `labels: "all"` or `"predicate_only"`, or with a predicate (any mode) | the rows the memory account refused, each once per phase (§14.4; a predicate's selection first, phase `selection`, §19.9) |
| `anchors_truncated` | list | answered, `labels: "all"` | the rows cut at `max_labels_per_anchor`, each once (§14.2) |
| `labels_cut` | object \| null | answered, `labels: "all"` | `partial`: `{reason: "max_labels", returned}` when `by_label` and the results list fewer labels than were found |
| `occurrences_cut` | object \| null | answered, `labels: "all"` | `partial`: `{reason: "max_occurrences_per_label", labels}` when that many labels list fewer occurrences than they have |
| `labels_excluded_unverified` | count | answered, `labels: "all"`, L > k, `long_search: "paths"` and `require_support: "record_verified"` | increment 4, unit `labels`: the labels carrying a returned path but verified on none, left out of `by_label` (§12.1) |
| `selection` | object (`selection`, §19.10) | answered, with a predicate | increment 5b: what this pattern's selection did (`pass`), what a label's presence means (`support`) and how the predicate's labels were read (`access`) |
| `absence_filter` | `"predicate"` | answered, with a predicate | increment 5b: the entry's absence claims are narrowed by the predicate (§19.11); `absence_scope` keeps its values |

### 8.6 `counts`

<!-- schema: counts -->
| field | type | present | meaning |
|---|---|---|---|
| `contexts` | contexts count | L ≤ k | unit `graph_contexts` |
| `anchors` | anchors count | L > k | unit `anchors` (§7.7) |
| `paths` | count | L > k | unit `paths`; `unknown` (`exact` 0 without anchors); with `long_search: "paths"` a paths count (`paths_count` below, §12.1) |
| `labels` | count | always | unit `labels`; `unknown` unless `labels: "all"` (or `"predicate_only"`: the predicate's labels on the returned contexts, §19.7) read them (§14.6); for the paths of `long_search: "paths"` with `by_support` besides (`labels_count` below) |
| `occurrences` | count | always | unit `placed_occurrences`; `unknown` unless `labels: "all"` (or `"predicate_only"`) placed them in records (§14.6) |
| `tested` | count | with a predicate | increment 5b, unit `graph_contexts` (L ≤ k; `paths` for L > k): the raw contexts the selection decided; relations §19.7 |
| `selected` | count | with a predicate | increment 5b, same unit: the contexts whose labels satisfy the predicate; relations §19.7 (never with an `estimate`) |

<!-- schema: count -->
| field | type | meaning |
|---|---|---|
| `value` | integer \| null | §7.4; `null` iff `unknown` |
| `relation` | `"exact"` \| `"at_least"` \| `"bounds"` \| `"unknown"` | §7.4 |
| `unit` | string | §7.4 |
| `lower` | integer | `bounds` only |
| `upper` | integer | `bounds` only |
| `estimate` | integer | `bounds` counts of a graph without its mask only (`counting: "upper_bound"`): round(`upper` × `index.dummy_fraction.value`), kept inside [`lower`, `upper`]; an estimate, not a bound, never `exact` (§7.4; the owner's decision #16) |

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

Increment 4: with `long_search: "paths"`, the paths count of a pattern longer than k has these fields besides
(§12.1):

<!-- schema: paths_count -->
| field | type | meaning |
|---|---|---|
| `by_strand` | object → count | BASIC graphs: the paths per orientation searched, unit `paths` |
| `by_orientation` | object → count | CANONICAL and PRIMARY graphs |
| `candidates_examined` | integer | the branches the extension entered (prefixes of k + 1 to L bases, complete paths included): work, not a count of the pattern |
| `extension` | `"no_anchors"` \| `"not_started"` \| `"not_admitted"` \| `"stopped"` \| `"completed"` | what the extension did: `completed` (relation `exact`), `no_anchors` (`exact` 0: the anchors `exact` 0), `stopped` (`at_least`: a stop in the extension), `not_started` (`unknown`: a stop before it), `not_admitted` (`unknown`: the anchors `exact` above `max_anchors`) |

and the labels count of the paths (with `labels: "all"`) has:

<!-- schema: labels_count -->
| field | type | meaning |
|---|---|---|
| `by_support` | object (`by_support` below) | the labels split by their support over the returned paths |

<!-- schema: by_support -->
| field | type | meaning |
|---|---|---|
| `record_verified` | count | unit `labels`: the labels one record verifies on at least one returned path; `unknown` without `record` placement |
| `label_intersection` | count | unit `labels`: the listed labels verified on none (carried by every k-mer of a path, no record holding it whole); the two sum to `counts.labels` when all three are `exact` |

### 8.7 `work`, `stop`, `withheld`, `cut`

<!-- schema: work -->
| field | type | meaning |
|---|---|---|
| `ranges_visited` | integer | range evaluations (one step each) |
| `mask_scans` | integer | ranges whose deferred scan began (§7.6): a range's masked edges — 0 on BASIC, CANONICAL and odd-k PRIMARY graphs unless masked edges lie among the candidates — and, on an even-k wrapped PRIMARY graph, the palindrome check of each range at a palindrome-capable offset (about one per such range). On a graph without its mask only the palindrome checks exist (each also tells a source dummy from a k-mer, §7.4), and the ranges whose few unchecked candidates were checked (§7.4, the owner's decision #24) |
| `steps` | integer | every step this pattern charged (k − 1 for each candidate the check of §7.4 tested) |
| `annotation_rows` | integer | increment 3, `labels: "all"`: the rows this pattern's reads returned (both steps) |
| `annotation_units` | integer | likewise: the work units of this pattern's reads, refused ones included (§14.4) |
| `memory_bytes` | integer | likewise: the request's memory account at its peak so far (the model of §14.4) |
| `extension_edges` | integer | increment 4, entries of a pattern longer than k with `long_search: "paths"`: the outgoing edges the extension examined, one step each (part of `steps`) |
| `extension_anchors` | integer | review GPT-3 (§18), likewise: the anchors whose extension began, each spelled once (k − 1 BOSS steps that no step charges) before its depth-first search; every anchor when the extension `completed`, 0 when it did not run (`no_anchors`, `not_started`, `not_admitted`), at most the anchors listed when a stop cut it. Not part of `steps`; deterministic unless a time stop cut the extension |
| `extension_branches` | integer | likewise: the nodes the extension expanded (anchors included) with two or more k-mers allowed at the next pattern position, where its paths fan out (each enters two candidates or more, `counts.paths.candidates_examined`); 0 when it did not run. Not part of `steps` |
| `annotation_rows_distinct` | integer | review GPT-3 (§18), every entry with `labels: "all"`: the distinct rows whose labels this pattern read (complete or truncated), each once however many contexts or path k-mers share it, the placement's second read not counted again; refused and unread rows not counted (0 where nothing was read, e.g. a withheld count). At most `annotation_rows`, which counts the reads of both steps |
| `verification_steps` | integer | likewise, entries of a pattern longer than k with `long_search: "paths"`: the verification's units of work (§12.1) — one per k-mer row looked up for a label carrying a path, per list ordered, per galloping seek, per run of chains extended and per record a run crosses; the work time is read at least every 4,096 of them. 0 where no coordinate is read. Not part of `annotation_units`; deterministic |
| `predicate_rows` | integer | increment 5b, every answered entry with a predicate: the rows the selection read (refused rows not counted); with `"predicate_only"` the projection's rows are these (its `annotation_rows` are its placement's reads only, `annotation_rows_distinct` the distinct rows of the returned contexts). With a predicate `memory_bytes` is in every answered entry (the account exists in every mode) |
| `predicate_units` | integer | likewise: the selection's work units (§19.9), refused and interrupted reads included; the request's `max_predicate_work` bounds their sum over the patterns |
| `predicate_lookups` | integer | likewise: the reverse-complement lookups made (`"either"` on a BASIC graph; k units each), one per distinct row |

<!-- schema: stop -->
| field | type | meaning |
|---|---|---|
| `phase` | `"discovery"` \| `"mask_scan"` \| `"extraction"`; increment 3: `"label_discovery"` \| `"placement"` \| `"output"`; increment 4: `"extension"`; increment 5b: `"selection"` | §7.6 |
| `reason` | `"max_steps"` \| `"time"` \| `"max_contexts"` \| `"max_anchors"`; increment 3: `"max_annotation_work"` \| `"max_memory"`; increment 4: `"max_paths"` (phase `extension` only); increment 5b: `"max_predicate_work"` (phase `selection` only), `"max_predicate_contexts"` (phase `discovery`) | §7.6 |

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
| `labels` | list \| null | likewise: the context's labels (§14.5); `null` unless the row was read and held. With `"predicate_only"` (increment 5b) the predicate's labels on its own row |
| `selection_labels` | list of strings | increment 5b, a selected context of a predicate request whose projection reads labels (`"predicate_only"`, `"all"`): the predicate's labels in the set it was evaluated on (its row's, and with `"either"` on a BASIC graph its reverse complement's), in the label order of the returned results: why it was selected; the normal form holds on them (§19.10) |
| `selection_strands` | list of strings | increment 5b, beside `selection_labels` (one per label, the same order): the orientation whose row carries the label (the owner's answer to P11). On a BASIC graph `"context"` (the row of the context's k-mer x, as `kmer` spells it), `"reverse_complement"` (the row of rc(x) only: `"either"`), `"both"` (both rows; a palindromic x, which is its own reverse complement, under `"either"`); `"either"` on CANONICAL and PRIMARY graphs (one row serves x and rc(x): no strand is known) (§19.10) |

With `output.labels: "none"` a result has the first seven fields only: nothing of the annotation is read. A
path (a result of a pattern longer than k under `long_search: "paths"`) has its own fields (`path_result`,
§12.1): never `kmer`, `node` or `row`.

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
| `bad_alphabet` | a character outside the kind's alphabet, or an empty pattern (the message names the first offending 0-based position); for `protein`, a character that is not in `protein_residues` (the stop `*` is one, §12.2; the message, kept as `4596bb3b` wrote it, lists the 20 amino acids and X B Z J and does not name `*`) | `id`, `kind`, `error` |
| `information_below_floor` | below the floor for its scope (§7.8) | `id`, `kind`, `pattern`, `length`, `information_bits`, `anchor_information_bits`, `min_anchor_information_bits`, `error`; a peptide also `residues`, `genetic_code` |
| `scope_unsupported` | `scope: "suffix"` on a wrapped PRIMARY graph, L ≤ k (§7.2) | as above |

### 8.10 Notes

| note | meaning |
|---|---|
| `low_complexity_pattern` | an exact pattern that sdust flags over its whole text (T = 20, W = 64, the seeder's parameters): a hint that its counts may be large; for L > k the flag can come from bases outside the anchor windows. Since the review GPT-3 (§18) stated only in an answer without a stop (never beside one, whatever its phase and reason, an earlier pattern's request-wide stop included), and left out of a pattern of more than 191 bases whose diagnostic the work time cut (`determinism: "time_limited"`, `stop` `null`; §7.8). Its absence is no claim that a pattern is not low-complexity |
| `strand_unknown_canonical` | a CANONICAL or PRIMARY graph: orientations, not strands |
| `paths_later_increment` | L > k without `long_search: "paths"`: anchors counted, paths neither extended nor extracted (ask with `long_search: "paths"`, §12.1) |
| `annotation_unbudgeted` | increment 3: the labels were read without the budget-aware decode (`allow_unbudgeted_annotation`): no memory bound on the reads themselves (what they returned is in the account), the deadline checked between chunks of keys. Increment 5b: also an entry whose selection read rows so, without a projection |
| `record_bounds_unknown` | increment 3: coordinates without a record mapping (`placement: "global"`): occurrences are (`kmer_coord`, `offset`), placed in no record, not deduplicated, not counted |
| `annotation_not_read` | increment 3: the request asked for labels (`output.labels: "all"`) or named an annotation field (`require_support` included), and this answer reads none (mode `count`, or `labels: "none"`): the fields had no effect |
| `label_intersection_only` | increment 4: the labels of paths where no coordinate is read (placement `none`, `none_canonical`, `not_requested`): every label is `label_intersection`, none can be verified (§12.1) |
| `threshold_upper_bound` | the owner's decision #16, a graph without its mask: a threshold was decided on a count's upper bound U against the request while its lower bound was within the threshold — `all_or_count` withheld `count_above_threshold`, `stop_at_threshold` stopped (`threshold_crossed`, or in `partial` `cut: max_contexts` / `max_anchors`), or the extension was not admitted (`anchors_above_threshold`): the true count may be within the threshold (§7.5) |
| `no_stop_codon` | the owner's decision #19: a peptide holding `*` read in a genetic code without an unconditional stop codon (tables 27, 28, 31): `*` matches nothing there, so the pattern has no instance; its contexts (L ≤ k, `exact` 0, answered without a search, §7.6) or paths (L > k) are 0 for that reason. A long one's anchors are its anchor windows', which may not reach the `*` (§12.2) |
| `estimate_sampled_dummy_fraction` | the owner's decision #16: a count of the entry carries an `estimate` (a `bounds` count of a graph without its mask): round(U × f), f sampled (`index.dummy_fraction`); not a bound (§7.4) |
| `predicate_constant` | increment 5b: the predicate's normal form is a constant (§19.4): nothing was read for the selection (`selection.pass: "constant"`) |
| `projection_not_read` | increment 5b: a predicate request named a projection field (`output.labels` `"all"` or `"predicate_only"`, `max_labels_per_anchor`, `max_annotation_work`, `max_labels`, `max_occurrences_per_label`, `require_support`) and this answer built no projection (mode `count`, or `labels: "none"`): the fields had no effect. A predicate answer never says `annotation_not_read`: its selection read rows |

## 9. What an answer licenses

- `retrieval_complete: true` (retrieval modes only): every graph context of the pattern in its scope and
  strands is in `results`, so a k-mer absent from them contains no instance there. With `output.labels:
  "none"` it claims nothing about any label, sample or record: none were read (with `"all"`, below).
- An `exact` count is the number of graph contexts (or anchors) in the scope; `exact` 0 is an absence claim
  over the index's retained k-mers in that scope (§7.2) — of contexts, not of labels. An `exact` count of
  labels or placed occurrences (`output.labels: "all"`) is that number; its `exact` 0 is an absence claim of any
  label (or placed occurrence) of the pattern in that scope.
- An incomplete list is not an inexact count: a `partial` list cut at `max_contexts` after a completed
  discovery keeps the `exact` count of contexts, and a labelled answer cut at `max_labels` or
  `max_occurrences_per_label` (`labels_cut`, `occurrences_cut`) keeps its `exact` counts of labels and
  occurrences (§14.5); but an item missing from such a list is not absent. An absence of one item (a context, a label, an
  occurrence) needs the complete list (`retrieval_complete: true`) or an `exact` 0 count of its unit.
- `at_least`, `bounds` and `unknown` claim what they say and no more (§7.4). An `estimate` claims nothing: it is
  neither a lower nor an upper bound, and no absence or presence follows from it (on a graph without its mask a
  `bounds` count with `lower` 0 may be all source dummies). A graph without its mask licenses exactly what a
  masked one does for its `exact` counts (an `exact` 0 is an absence, §7.2) and its complete lists.
- A `count` answer, or any answer with `output.labels: "none"`, establishes no absence of a label, a sample or
  a record, whatever its counts (design §5.1).
- Increment 4, `long_search: "paths"`, a pattern longer than k: `retrieval_complete: true` says that every path
  of the graph spelling an instance of the pattern (in its strands) is in `results`; an `exact` 0 count of paths
  says that no walk of the index's retained k-mers spells it. A label of a path is a claim about records only
  when its support is `record_verified` (one record holds the whole path, at the placed coordinates);
  `label_intersection` says only that every k-mer of the path carries the label (design §4.3), and its
  `occurrences` `exact` 0 that no record of the label holds the whole path. With `require_support:
  "record_verified"` a label left out (counted in `labels_excluded_unverified`) is not absent from the path's
  k-mers. With labels, `retrieval_complete: true` also says that every label carrying a returned path is listed
  with its support (§12.1).
- Increment 5b, a predicate (`absence_filter: "predicate"`): the absence claims of the entry are those of §19.11,
  per context and per index, never motif-level: `retrieval_complete: true` says that every context of the
  pattern whose label set (at `predicate.strands`) satisfies the predicate is in `results`, `selected` `exact` 0
  that none does; the raw counts keep the licences above.
- Increment 3, `output.labels: "all"`, `retrieval_complete: true`: every graph context of the pattern in its scope
  and strands is in `results` with **all** the labels its row carries (`labels_status: "complete"` everywhere),
  each placed where `placement` is `record` or `global` — so a column absent from `by_label` carries no
  occurrence of the pattern in the index's retained k-mers in that scope, and on a `record` index every placed
  occurrence of the listed columns is in the results (`counts.occurrences` exact). Without it (withheld, cut,
  truncated, refused, stopped) the labels and occurrences listed are true but not all: an omission from the lists
  establishes no absence, and each count keeps its own stated relation (it stays `exact` after an output-only cut
  such as `labels_cut` or `occurrences_cut`, and is `at_least` or `unknown` where the reads themselves stopped).

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
| `available` | boolean \| null | | `true`: `/pattern` answers on this graph, with its mask or without it (`counting`); `false`: not on this graph as loaded (`unavailable_reason`): `mask_invalid` clears when the server is restarted after the graph was masked again (`transform --mask-dummy --force`), `multi_graph_later_increment` with the increment that serves multi-graph servers, `alphabet_untested` with a build that has passed the pattern tests on DNA5 (an older build's `mask_required` when the server is restarted after the mask exists); the other reasons are permanent for this graph; `null`: the index is loading |
| `unavailable_reason` | string \| null | | `mask_invalid`, `representation_unsupported`, `primary_unwrapped`, `alphabet_untested`, `alphabet_unsupported`, `multi_graph_later_increment`; `mask_required` is retired (§6: builds of version 1 before the owner's decision #16 state it; this one never does). A later build may add others: pass an unknown one through, §1; `null` when available or loading |
| `modes` | list | `["count", "all_or_count", "partial"]` | §7.5 |
| `default_mode` | string | `"all_or_count"` | an omitted `mode` |
| `projections` | list | `["none", "all", "predicate_only"]` | the `output.labels` values served **now**; gate label projections on this list (`"all"` since increment 3, `"predicate_only"` since increment 5b, with a predicate; `["none"]` before increment 3) |
| `default_projection` | string | `"none"` | an omitted `output.labels`; frozen for version 1 (§1) |
| `projections_later_increment` | list | `[]` | refused with `later_increment` today (`["predicate_only"]` before increment 5b, `["all", "predicate_only"]` before increment 3) |
| `default_occurrences` | boolean | `true` | increment 3: an omitted `output.occurrences` with `labels: "all"` |
| `kinds` | list | `["dna", "iupac", "protein"]` | pattern kinds served (`protein` since increment 5; `["dna", "iupac"]` before) |
| `kinds_later_increment` | list | `[]` | (`["protein"]` before increment 5) |
| `protein_residues` | list of one-letter strings | the 20 amino acids then `X`, `B`, `Z`, `J`, `*` | increment 5: the residues a `protein` pattern may hold (§12.2); the stop `*` since the owner's decision #19 of 2026-10-08 (it was not among them before) |
| `genetic_codes` | list of integers | `[1, 2, 3, 4, 5, 6, 9, 10, …, 16, 21, …, 33]` | increment 5: the NCBI translation table ids `genetic_code` accepts (gc.prt version 4.6) |
| `default_genetic_code` | integer | `1` | increment 5: an omitted `genetic_code` (the standard code) |
| `protein_rule` | string | `"SPEC-pattern-search.md sections 12.2, 18"` | increment 5: a reference to the rule of how a peptide is read (§12.2: the ambiguity codes, the stop `*`, the codon automaton, the length in bases, the slot error, the context stops of tables 27, 28 and 31 and the note `no_stop_codon`; §18); since the owner's decision P9 (§18) a reference, the rule in prose before. Printable ASCII; for people, not parsed |
| `default_scope` | string | `"any_offset"` | |
| `scopes_by_graph_mode` | object | `basic`, `canonical`: `["suffix", "any_offset"]`; `primary`: `["any_offset"]` | the rule |
| `scopes` | list \| null | | this graph's requestable scopes |
| `long_patterns` | string | `"anchors_counted"` | what L > k gets without the option (§7.7); stays `"anchors_counted"` now that paths are served, since they are opt-in (`long_search`, §12.1) |
| `long_search` | list | `["anchors", "paths"]` | increment 4: the `long_search` values served; gate the paths on `"paths"` in it |
| `default_long_search` | string | `"anchors"` | increment 4: an omitted `long_search`; the paths are never switched on by a default |
| `strands` | list | `["both", "forward", "reverse"]` | |
| `default_strands` | string | `"both"` | |
| `graph_cleaned` | string | `"unknown"` | whether graph cleaning may have pruned k-mers (design §3); not known in version 1 |
| `records_shorter_than_k` | string | `"not_indexed"` | such records have no k-mer |
| `resident_only` | boolean | `true` | the route never loads an index (`in_ram` refused) |
| `caps` | object | | the maxima (§4.5): `max_contexts`, `max_anchors`, `max_steps`, `time_budget_ms`, `min_information_bits` (the floor), `max_patterns`; increment 3: `max_labels_per_anchor`, `max_annotation_work`, `max_memory_mb`, `max_labels`, `max_occurrences_per_label`; increment 4: `max_paths`; the owner's decision #24: `max_checked_entries` (no request field: the unchecked candidates a pattern on a graph without its mask may have for each to be tested, §7.4; on every server, masked or not); increment 5b: `max_predicate_contexts`, `max_predicate_work` and `max_predicate_labels` (no request field: the names a predicate may list, §19.3) |
| `default_time_budget_ms` | number | 60,000 | the budget of a request that names none, below `caps.time_budget_ms` |
| `finalize_reserve_ms` | number | 250 | §7.6 |
| `caps_rule` | string | | which caps are request fields' maxima (each named: lowered and listed in `limits.clamped` above it) and which are the server's policy (`max_patterns`, `min_information_bits`, `max_checked_entries`, `max_predicate_labels`), and a reference to the sections stating the rules (`SPEC-pattern-search.md sections 4.5, 7.4, 7.6, 12.1, 19`: the caps and defaults, the check of few unchecked candidates, the time kept back for the answer, the two admissions of `long_search: "paths"`, a predicate's admissions and budgets). Since the owner's decision P9 (§18) a reference: it stated those rules in prose before, the delivery rates among them (now `delivery_mbps`). Printable ASCII; every cap of `caps` is named in it; for people, not parsed |
| `delivery_mbps` | object (`delivery_mbps` below) | | the owner's decision P9 (§18): the rates in force of the time kept back for the answer (§7.6) |
| `predicate` | object (`capabilities_predicate` below) | | increment 5b (§19.12): the operators of a predicate served, the `predicate_strands` values and the access this index gives a selection |
| `graph_mode` | string \| null | | `basic`, `canonical`, `primary` |
| `k` | integer \| null | | |
| `alphabet` | string \| null | | `$ACGT` or `$ACGTN` (`$ACGTN`: `available: false`, `alphabet_untested`, §8.2) |
| `strand_stated` | boolean \| null | | `true` on BASIC |
| `mask` | string \| null | | `file`: the `.edgemask` loaded beside the graph (below; one that marks a W = `$` edge valid: `available: false`, `mask_invalid`); `built_at_load`: built in memory at start-up (`--pattern-build-mask`, milestone 1b); `absent`: no mask, the graph served with `counting: "upper_bound"` (the owner's decision #16; before it, `mask_required`; a `$ACGTN` graph says `alphabet_untested`, §6) |
| `counting` | `"exact"` \| `"upper_bound"` \| null | | the owner's decision #16 of 2026-10-08: how this graph counts. `exact` with a mask (`file`, `built_at_load`): every count of a completed discovery `exact`; `upper_bound` without one (`absent`): counts `exact` where provable, otherwise `bounds` with an `estimate`, lists exact (§7.4); `null` while loading or when the graph is not served (`available` not `true`). An extensible enumeration (§1). The answers on such a graph state it in `index.counting` |
| `dummy_fraction` | object (`dummy_fraction`, §8.2) \| null | | the owner's decision #16: with `counting: "upper_bound"`, f, the dummy fraction the estimates rest on (sampled at load, the same object as in each answer's `index`); `null` with a mask, while loading or when not served |
| `placement` | string \| null | | the placement `output.labels: "all"` gives on this index (`record`, `global`, `none`, `none_canonical`, §14.3; before increment 3: what a later increment could give) |
| `support` | string \| null | | the best per-label support of a path on this index (`record_verified`, `label_intersection`; served since increment 4, §12.1: `require_support: "record_verified"` needs `record_verified` here); a context of L ≤ k has `kmer` |
| `annotation` | string \| null | | `budgeted` (row-diff with budgeted decode) or `unbudgeted`: the reads of `output.labels: "all"` (§14.4; `unbudgeted` needs `allow_unbudgeted_annotation`); `count` and `none` never read it |

The graph fields (`graph_mode` to `annotation`) are `null` while the index loads; when the graph is not
recognised (`representation_unsupported`, `primary_unwrapped`) only `k` is set; `placement`, `support` and
`annotation` are `null` if the annotation could not be described.

<!-- schema: delivery_mbps -->
| field | type | meaning |
|---|---|---|
| `build` | number | MB/s at which the answer's JSON text is assumed to be built and written (`--pattern-delivery-build-mbps`, 10; a positive number, the server refuses another at start-up) |
| `compress` | number | MB/s at which it is assumed to be compressed (`--pattern-delivery-compress-mbps`, 50; likewise) |

<!-- schema: capabilities_predicate -->
| field | type | meaning |
|---|---|---|
| `operators` | list | `["any", "all", "none", "at_least", "and", "or", "not"]`: the operators a predicate may use (§19.3) |
| `strands` | list | `["either", "context"]`: the `predicate_strands` values (on CANONICAL and PRIMARY graphs both are evaluated as `"either"`, §19.5) |
| `access` | `"rows"` \| `"columns"` \| null | how a selection reads this index's annotation: `"rows"` (the budget-aware decode, or an unbudgeted annotation without direct access) or `"columns"` (an unbudgeted annotation with direct access: single cells for at most 16 labels of a predicate, rows above); `null` while the index loads or when the annotation could not be described |

**The mask, for a client** (the operator's side is design §4):
- It is read when the graph is loaded: an `.edgemask` written beside a running server changes nothing until
  the server restarts.
- Without one (`mask: absent`) the graph is served all the same since the owner's decision #16 of 2026-10-08
  (`counting: "upper_bound"`, `dummy_fraction`; §7.4): counts `exact` where provable, otherwise `bounds` with an
  `estimate`; every list exact. The mask makes every count of a completed discovery `exact`; it never changes an
  `exact` count (an `exact` count of the same index and scope is the same number with and without it).
- Once a graph has its `.edgemask`, every loader reads it, on every route: node ids, rows and what `/search`,
  `/align`, `/resolve` and `/traverse` find stay as they were, but `GET /stats` `graph.nodes` becomes the
  graph's k-mer count (the dummy edges no longer counted), and a `.bloom` beside the graph is loaded from then
  on. The index identity `index_fp` does not change: the mask and the Bloom filter are derived data of the
  graph, not part of it (the owner's decision #17 of 2026-10-08; `SPEC-labeled-traversal-core.md`, "The index
  identity"), so a manifest made without them stays valid, and one that lists them is refused at start-up.
  `mask` (`file`, `built_at_load`, `absent`) and `counting` say which counts are exact. `--pattern-build-mask`
  writes no file: `/stats` changes, and `mask: built_at_load` tells it from a file.
- `metagraph extend` on a masked graph wrote, before this build, a mask that marks the new dummy edges valid
  (this build's `extend` rebuilds the mask of a masked graph as `transform --mask-dummy` builds it), and the
  counts on such a graph could be overstated, `exact` included. Such a mask, or a stale one, is refused (the owner's decision of
  2026-10-07): a mask that marks valid a dummy edge whose W is `$` makes the graph `available: false` with
  `mask_invalid`, and every request 400 `mask_invalid`, until the graph is masked again (`metagraph transform
  --mask-dummy --force`) and the server restarted. Beyond that check a `file` mask is trusted as written. A mask
  written by `build --mask-dummy`, `transform --mask-dummy` or `--pattern-build-mask` is correct.

<!-- schema: capabilities_multi -->
| field | type | meaning |
|---|---|---|
| `pattern_contract_version` | integer | 1 |
| `available` | boolean | `false` |
| `unavailable_reason` | string | `multi_graph_later_increment` |

### 10.3 How a client gates

- Pattern search is offered on a host only when both hold: the block states a `pattern_contract_version` the
  client implements (§1), and `available: true`. No block, or a missing, malformed or unsupported version
  (lower or higher): no pattern search on this host (`pattern_unsupported` on the service).
- `available: false`: the host has the route but not for this graph; keep and show `unavailable_reason` (a
  `mask_invalid` host, or an older build's `mask_required` one, becomes available once its operator masks the
  graph and restarts the server: re-probe after a restart).
- `counting` (the owner's decision #16): `exact` hosts answer every completed count `exact`; `upper_bound`
  hosts answer `bounds` with an `estimate` where they cannot prove a count. Show an `estimate` as an estimate
  (never as the count, never as a bound) beside its `lower` and `upper`, and state the host's `dummy_fraction`
  with it; gate any use that needs exact counts (an absence claim from a count other than `exact` 0, a total
  across hosts) on `counting: "exact"` or on the count's own relation. A host whose block lacks `counting` (a
  build before the decision) counts exactly when it is available. The thresholds of an `upper_bound` host admit
  on U: a pattern whose true count fits `max_contexts` can be withheld (note `threshold_upper_bound`), and
  `partial` lists it whatever U is.
- `available: null`: availability is not known yet (the index is loading): re-probe after the load, never
  read it as `false`.
- Offer a projection, kind or scope only when its list has it; send values within `caps`; send a field of a
  later increment (§4.4) only once the capabilities announce it.
- `output.labels: "all"`: offer it where `projections` has it; on `annotation: "unbudgeted"` it needs
  `allow_unbudgeted_annotation: true` (else 400 `annotation_unbudgeted`), which a client sends only with its
  caller's explicit consent to reads without a memory bound (§14.4); `placement` says whether occurrences come
  in records (`record`), as coordinates (`global`) or not at all.
- Paths (increment 4): send `long_search: "paths"` only where `long_search` lists `"paths"`; `max_paths` within
  `caps.max_paths`; `require_support: "record_verified"` only where `support` is `record_verified` (else 400
  `support_unavailable`). A host whose capabilities lack `long_search` refuses it (400 `later_increment`, a build
before increment 4).
- Peptides (increment 5): send `protein` patterns only where `kinds` lists `"protein"`, with residues of
  `protein_residues` (the stop `*` only where that list has it), and `genetic_code` only from
  `genetic_codes`.
- Predicates (increment 5b, §19): send `predicate` only where `projections` lists `"predicate_only"` (and the
  block has `predicate`), with operators of `predicate.operators`, at most `caps.max_predicate_labels` names,
  `predicate_strands` from `predicate.strands`, `max_predicate_contexts` and `max_predicate_work` within their
  caps; on `annotation: "unbudgeted"` it needs `allow_unbudgeted_annotation: true` in every mode. A host whose
  capabilities lack `"predicate_only"` refuses a predicate (400 `later_increment`, a build before increment 5b).

## 11. Fixtures

`api/python/tests/data/traverse/pattern/` holds one directory per situation with the body sent
(`request.json`; `null` for a GET) and the body answered (`answer.json`), generated by
`scripts/traversal/pattern_fixtures.py` from a server of this build on copies of the mini index
(`build/mini_refseq`, a BASIC index of refseq33m's format at k = 31: its graph given the `.edgemask` by
`metagraph transform --mask-dummy`, as a host gets it; served as built for `mask: absent`, `counting:
"upper_bound"` (the `unmasked_*` fixtures, §18; server `unmasked` at the default `--pattern-max-checked-entries`,
server `unmasked_unchecked` with 0), and with `--pattern-build-mask` for `mask: built_at_load`;
served with `--no-coord-mapping` for `placement: global`), a
PRIMARY index of two of its record files (a column annotation: `annotation: unbudgeted`), and a hash graph of one
of its records (a graph the engine does not recognise: `representation_unsupported`). `index.json` states each
fixture's server, method, path and status; `README.md` says in one line what each shows. They cover the
capabilities on both routes (each `mask` value: `file`, `built_at_load`, `absent`; multi-graph; PRIMARY), each
mode, each `withheld` and `cut` reason, both scopes and every strand setting, IUPAC patterns, a palindrome, a
pattern longer than k, the three error slots, the clamps, the relation `bounds` with a stop in a deferred scan
(`max_steps_bounds`, `max_steps_bounds_withheld`), a `max_steps` stop whose release then met the work time
(`max_steps_then_time`: `stop` keeps the first stop, the time shows as `cut: time` and `time_limited`, §7.6), a
graph the engine does not recognise
(`*_representation_unsupported` on both capabilities routes, `representation_unsupported`), and every
whole-request refusal a server of this build gives but `mask_invalid` and `alphabet_untested` (below; the
fixture validator names them in `NO_FIXTURE` and checks its code lists against the codes the sources write;
`mask_required`, which no build since the owner's decision #16 gives, is in its `RETIRED`, with no fixture) (two 503 bodies are hand-made: no server produces them on
demand; their `request.json` is illustrative, and a live server answers it 200). Two codes have no fixture: `primary_unwrapped` (`server_query` and the CLI always wrap a PRIMARY
graph; only an embedding can reach it) and `alphabet_unsupported` (no graph of another alphabet loads in this
build). Increment 3 adds `labels_all` (record placement within the threshold: 9 columns, 42 placed occurrences
of the blaNDM-1 forward primer, equal to a scan of the mini's FASTA),
`labels_all_withheld` (above the threshold: nothing read), `labels_all_truncated` and
`labels_all_truncated_partial` (`max_labels_per_anchor` 2), `labels_all_partial` (the label and occurrence cuts),
`labels_all_work_budget` and `labels_all_work_budget_partial` (a stop at `max_annotation_work`),
`labels_all_global` (no record mapping), `labels_all_count` (mode `count`: `annotation_not_read`),
`annotation_unbudgeted` (the 400 on the PRIMARY index's column annotation) and `labels_all_unbudgeted` (the same
with `allow_unbudgeted_annotation`); `later_increment_labels` asked for `"predicate_only"` (`"all"` being served)
until increment 5b replaced it by `predicate_only_without_predicate` (below).

After the outside review of 2026-10-07 (GPT #11): `labels_all_output_budget` (a GCG repeat whose 1,828 contexts
fill the smallest account, `max_memory_mb` 1: `stop {output, max_memory}`, `withheld: output_budget`, the contexts
count `exact`), `labels_all_partial_exact_cut` (`max_labels` and `max_occurrences_per_label` cut the lists while the
counts stay `exact`) and `labels_all_mixed_slots` (refused and answered slots side by side with labels "all"; also
`bad_alphabet` and `information_floor` without labels); after its recheck `labels_all_rows_refused` (four patterns
sharing one account of 1 MB in `partial`: the last one's rows are refused, `rows_refused` and `labels_status:
"refused"`, `cut: max_memory`, its labels counted `at_least`).

Increments 4 and 5 (§17) add the paths: `paths` (label-free, with an absent 40-mer: `no_anchors`, complete),
`paths_count` (mode `count`), `paths_labels` (blaNDM-1's first 40 bases, every label `record_verified` with its
whole-path occurrences, 9 columns and 42 occurrences as for the primer; and a 51-mer joining two copies of a
repeated 31-mer of the E. coli records, a path of the graph that no record holds whole: its labels
`label_intersection`, occurrences `exact` 0), `paths_require_support` (the same with `require_support:
"record_verified"`: the chimera lists no label, `labels_excluded_unverified` 2), `paths_max_paths_partial` (`cut:
max_paths`), `paths_count_above_threshold`, `paths_stop_at_max_paths` (`stop {extension, max_paths}`,
`threshold_crossed`), `paths_anchors_above_threshold` (`not_admitted`), `paths_max_steps` (`stop {extension,
max_steps}`, the paths completed before it released, `cut: max_steps`, the sticky stop after it), `paths_global`
(no record mapping: chains, every label `label_intersection`), `paths_primary` (a PRIMARY index: orientations, note
`label_intersection_only`) and `support_unavailable` (the 400 on the index without record mapping); and the
peptides: `peptide` (NDM-1's first 10 residues within one k-mer: 4 contexts, the 9 columns and 42 placed
occurrences of blaNDM-1; with J; and 14 residues without `long_search`: anchors only), `peptide_count`
(`genetic_code: 11`), `peptide_paths` (14 residues, 42 bases, through `long_search: "paths"`: one path per strand,
`record_verified`), `peptide_bad_residue` (`bad_alphabet` for U; its `*` pattern, refused `stop_unsupported` at
`4596bb3b`, is answered since §18) and `genetic_code_unknown` (the 400 for table 7). The nine capabilities bodies of the single-graph servers gained the
fields of §17; no other stored body changed.

Situations without a stored body: `withheld: annotation_budget` for a refused row (in `all_or_count`); `by_label:
null` in `partial` (the mini's label names are too short to exhaust the account there; a real answer of a tiny
long-label index in `data/traverse/pattern_validator/by_label_null_partial` checks that the validator accepts it);
an earlier stop followed by a time stop of the output of the labels
(it depends on the machine's speed); the refusals `mask_invalid` (it needs a graph extended after masking with an
older build; the unit test `PatternMask.MaskWithAValidSentinelIsRefused` shows it) and `alphabet_untested` (a DNA5
build). The first and the third are exercised by the unit tests of `tests/cli/test_pattern_retrieval.cpp` on small graphs
(the third by `PatternRetrieval.AWorkStopThenATimeStopOfTheOutput`, on a virtual clock: a work stop in discovery
kept as the first stop, then the output's time stop).
How a client merges answers is in `PROMPT-search-service-pattern.md` §3.1 item 2, not in a fixture.

The owner's decisions of 2026-10-08 (§18) replace `mask_required` (the 400 of the unmasked server) by the
unmasked server's answers: `unmasked_count` (mode `count`: blaNDM-1's forward primer `bounds` [2, 24] with
`estimate` 24, exact 24 with the mask; `GATGCCGGTGAACAAC`, the first 16 bases of an E. coli record whose first
k-mer no k-mer enters, `bounds` [2, 32] with `estimate` 32 while the masked graph counts 17: the estimate is not
a bound; an absent primer `exact` 0; `index.counting` and `index.dummy_fraction`), `unmasked_labels_all`
(`all_or_count` with labels: every candidate enumerated, the list and the counts exact, 9 labels and 42 placed
occurrences as on the masked graph), `unmasked_threshold_upper_bound` (`max_contexts` 17: withheld
`count_above_threshold` on U = 32, note `threshold_upper_bound`; the masked graph releases the 17),
`unmasked_stop_at_threshold` (the running upper bound stops discovery: `at_least` 0, `threshold_crossed`, the
note), `unmasked_partial` (`max_contexts` 5: the masked graph's first five, `bounds` [6, 24] after the release
raised the lower bound) and `unmasked_paths` (anchors spelled whole, `exact`; the paths the masked graph's);
`unmasked_count`, `unmasked_labels_all`, `unmasked_threshold_upper_bound` and `unmasked_partial` are served with
`--pattern-max-checked-entries 0` (server `unmasked_unchecked`, their bodies unchanged), what any pattern above
the limit gets; at the default the owner's decision #24 checks those few candidates: `unmasked_checked` (the
primer `exact` 24, the island start `exact` 17, its first 14 bases, 68 unchecked, still `bounds` [4, 72] with
`estimate` 72) and `unmasked_checked_dummies` (the island start in `suffix` scope, forward: its one candidate a
source dummy, `exact` 0). The
unmasked server's two capabilities bodies say `available: true`, `counting: "upper_bound"` and its
`dummy_fraction` (f 0.9999, interval [0.999434, 0.999982], the exact f 0.999955 inside it). And the stop `*`:
`peptide_stop` (NDM-1's last 9 residues and its stop, table 1: 4 contexts whose instances end in its stop codon
TGA; the same with X for `*`: none, X never matches a stop), `peptide_no_stop_codon` (table 27: `*` matches
nothing, note `no_stop_codon`; 30 bases: `exact` 0 with no work; 45 bases: anchors counted, paths `exact` 0) and
`peptide_no_stop_codon_after_stop` (the same 30 bases after a `max_steps` stop: `exact` 0, `stop` `null`,
`determinism: "full"`, §7.6).

Increment 5b (§19) adds the predicates, on the masked copy unless named: `predicate_filter` (GCG12, 1,828
contexts, `and(any 562, none 287)`, `max_contexts` 100, `predicate_only` without occurrences: 60 selected and
returned in 1,777 rows, complete, with the predicate block and `selection_labels`), `predicate_either` and
`predicate_context` (the blaNDM-1 forward primer, `none(546)`, mode `count`: 0 with `"either"`, 12 with
`"context"`; the first also `projection_not_read`), `predicate_at_least` (`at_least` 8 of the nine columns,
`"context"`: the 12 `+` contexts), `predicate_unknown_constant` (the typo `5622`: `unknown_labels`, normal form
`false`, `pass: "constant"`, no row read, `predicate_constant`), `predicate_unknown_folded` (`and(any 562, none
5622)` folded to `any(562)`), `predicate_vacuous` (`none(287)`, partial, `max_contexts` 5: `vacuous` true, 68
selected, the first 5, `cut: max_contexts`), `predicate_budget` and `predicate_budget_partial`
(`max_predicate_work` 1, `"context"`: `tested` exact 1, `selected` `bounds` [0, 1827], `stop {selection,
max_predicate_work}`, `predicate_budget` / `cut: max_predicate_work`; the next pattern's selection
`not_started`, sticky), `predicate_above_threshold` and `predicate_partial_admission` (`max_predicate_contexts`
100: `not_admitted`, nothing read / the first 100 tested, `cut: max_predicate_contexts`),
`predicate_selected_above` (`selected_above_threshold`, 60 > 10), `predicate_stop_at_threshold` (the selected
count's early exit, `stop {selection, max_contexts}`) and `predicate_stop_at_raw` (the raw count's, `stop
{discovery, max_predicate_contexts}`, `threshold_crossed`), `predicate_only_record` (`any(562)` with
`predicate_only`: 562's 13 placed occurrences), `predicate_selection_strands` (`any(546)`, `predicate_only`: the
12 `+` contexts selected by their own row, `selection_strands: ["context"]`, the 12 `-` ones by their reverse
complement's, `["reverse_complement"]`), `predicate_all` (`none(546)`, `"context"`, labels `"all"`: the
12 `-` contexts with their 7 labels), `predicate_long_anchors` (a 40-mer: its anchors' answer, selection
`not_started`), `predicate_unmasked` (server `unmasked_unchecked`: the raw `bounds` [2, 24] admitted on U, the
release exact, 24 selected), `predicate_unbudgeted` and `predicate_unbudgeted_allowed` (the PRIMARY index: 400
`annotation_unbudgeted` in mode `count`; with the opt-in `selection.access: "columns"`, `predicate.strands:
"either"` although `"context"` was asked, note `annotation_unbudgeted`), and the refusals `predicate_invalid` (a
number as a name), `predicate_paths_refused` (with `long_search: "paths"`), `predicate_too_large` (server
`masked_small_predicate_cap`, `--pattern-max-predicate-labels 4`: five names) and
`predicate_only_without_predicate` (replacing `later_increment_labels`). The nine capabilities bodies of the
single-graph servers gained §19.12's fields; no other stored body changed.

- `pattern_fixtures.py --check` regenerates them and compares: bodies byte for byte, except
  `timing.elapsed_ms` (and every other `timing` value of an entry: `label_discovery_ms`, `placement_ms`),
  `server_instance` and the counts and work of the first `time_limited` entry of an answer (and, when the
  clock cut its release, its `returned` and `results`), which vary between runs (the later ones, stopped by the
  same budget, are compared as they are). The labelled answers'
  `work.annotation_units` and `work.memory_bytes` are deterministic (§14.4) and compared.
- `api/python/tests/test_pattern_fixtures.py` validates every answer against the field lists of §§6, 8 and 10
  (the tables marked `schema`) and the rules that tie the fields together (sums and relations, one budget per
  request, clamps, the floor, a complete list against its counts), so that this document, the fixtures and
  the server cannot drift apart silently. This validates the stored bodies, not a live server: what a build
  answers is covered only by the `--check` runs below, where they run.
- `integration_tests/test_pattern.py` (`TestPatternFixtures`) runs `--check` with the binary under test, so a
  build whose `/pattern` or capabilities answers differ from the stored bodies fails its own integration tests,
  on a machine with `build/mini_refseq` (or `$METAGRAPH_MINI_REFSEQ`). CI has neither the mini index nor a base
  binary: of these checks it runs only `TestPatternFixtureBodies` (the body validation above); the classes that
  need them are skipped there (and counted as passed by `main.py`), and `$METAGRAPH_REQUIRE_GUARDS=1` makes such
  a run fail instead of skip.

## 12. What changes in later increments

Contract version 1 stays; each row adds, it does not change (§1). The capabilities announce each addition before
a client may use it. (Milestone 1b, the edge mask, is in this build and its fixtures: `mask: "built_at_load"` with
`--pattern-build-mask`, and the `mask_required` message naming `transform --mask-dummy` and that flag; since
§18 a graph without its mask is answered and `mask_required` is retired.)

| milestone (design §13) | request | answer | capabilities |
|---|---|---|---|
| 3: labels and placement — **served in this build (§14)** | `output.labels: "all"`; `max_labels_per_anchor`, `max_annotation_work`, `max_memory_mb`, `max_labels`, `max_occurrences_per_label`, `allow_unbudgeted_annotation`; `output.occurrences` | `labels` on each result (column, support, placed occurrences `seq_id`, `record`, `strand`, `nt_coords`, `nt_length`; or `kmer_coord`, `offset`); `by_label`; `counts.labels` and `counts.occurrences` known; `rows_refused`, `anchors_truncated`, `labels_cut`, `occurrences_cut`; `withheld` reasons `annotation_budget`, `anchor_labels_truncated`, `output_budget`; notes `record_bounds_unknown`, `annotation_unbudgeted`, `annotation_not_read` (`label_intersection_only` arrives with paths, milestone 4) | `projections` gains `"all"`; `default_occurrences`; five caps; `default_projection` stays `"none"` (§1) |
| 4: patterns longer than k, **opt-in** (the owner's decision of 2026-10-07) — **served in this build (§12.1)** | `long_search: "paths"` (default `"anchors"`), `max_paths`, `require_support`; `output.paths` accepted with either value | only for a request with `long_search: "paths"`: results for L > k are paths, with the new fields `sequence` (the L spelled bases), `anchor_kmer` (the anchor's k bases) and their node path (`nodes`, `rows`) — never `kmer`, which keeps its meaning, the k-mer of a context; `counts.paths` known with its split, `candidates_examined` and `extension`; the labels of each path with their `support` (`label_intersection`, `record_verified`) and `require_support`; `withheld` `anchors_above_threshold`, `cut` `max_paths`, stop phase `extension` and reason `max_paths`; 400 `support_unavailable`. A request without it gets the answer of §7.7 as before: anchors counted, `counts.paths` `unknown` (`exact` 0 without anchors), `withheld: paths_later_increment` and its note | `long_search` `["anchors", "paths"]`, `default_long_search` `"anchors"`; `caps.max_paths`; `long_patterns` stays `"anchors_counted"` (what a request without the option gets); `support` states the best support of a path |
| 5: peptides — **served in this build (§12.2)** | `patterns[i].protein`, `genetic_code`; the stop `*` since the owner's decision #19 (§18) | `kind: "protein"` with `residues` and `genetic_code` (instances name the codons); note `no_stop_codon` (§18; the slot error `stop_unsupported`, answered by `4596bb3b` only, is retired); 400 `genetic_code_unknown` | `kinds` gains `"protein"`; `protein_residues` (with `*` since §18), `genetic_codes`, `default_genetic_code`, `protein_rule` |
| 5b: predicates, patterns of L ≤ k — **served in this build (§19)** | `predicate`, `max_predicate_contexts`, `max_predicate_work`, `predicate_strands`, `output.labels: "predicate_only"` | the top-level `predicate` block (normal form, names, known, `unknown_labels`, `vacuous`, scope, strands); per entry `selection` (`pass`, `support`, `access`), `counts.tested` and `counts.selected` (relations §19.7), `absence_filter`, `work.predicate_rows`, `predicate_units`, `predicate_lookups`, `timing.selection_ms`, `rows_refused` with phase `selection`; results the selected contexts, with `selection_labels` under a projection that reads labels; `withheld` `predicate_above_threshold`, `selected_above_threshold`, `predicate_budget`; `cut` `max_predicate_contexts`, `max_predicate_work`; `stop` phase `selection`, reasons `max_predicate_work`, `max_predicate_contexts`; notes `predicate_constant`, `projection_not_read`; 400 `predicate_too_large` | `projections` gains `"predicate_only"` (`projections_later_increment` `[]`); `caps.max_predicate_contexts`, `max_predicate_work`, `max_predicate_labels`; `predicate` {operators, strands, access} |
| 5s and 5b's L > k part: supported paths | `long_search: "supported_paths"`, `supported_paths_level`; a predicate on supported paths | `counts.supported_paths`; the selection of the supported walks (strand-consistent: a label supports a walk on one strand as a whole) | `long_search` gains `"supported_paths"` |
| 6: multi-graph | `graphs` (as `/search` names graphs and chunks), `budget_split` | each result carries `graph`, `index_fp`, `release`; counts `by_shard` with `per_shard`; `stop` and `withheld` gain the shard; the merged order of design §8 | the block on multi-graph servers becomes available, with the resident graphs |

- What stays: every field of §8 with its type and meaning; the relations and their algebra; the absence licences
  of §9; the order of §7.9; the refusal envelope `{error, code}`; the answer to a request that does not use a
  later option: the opt-ins (`output.labels: "all"`, `long_search: "paths"`, the predicates) are never switched
  on by a default.
- What a version-1 client sees on a later host: more fields and values, read by presence or passed through;
  other work counters and stopping points (§1).

### 12.1 Increment 4, served: patterns longer than k as paths (`long_search: "paths"`)

Served by this build (`src/cli/pattern.cpp`, `src/cli/pattern_retrieval.cpp`, the extension of
`src/graph/alignment/pattern_search.cpp`; design §4.2, §4.3 "Label consistency for long", owner decisions #13 and
#14), as additions to version 1. Opt-in: only a request with `long_search: "paths"` gets any of it; a request
without it, or with `"anchors"`, is answered as §7.7 says, byte for byte as before (§17). A pattern of at most k
bases is answered alike under both values (its entry is the same; only `limits` gains the echo).

**Request.** `long_search: "paths"`; `max_paths` (integer ≥ 0, default and cap `caps.max_paths`, 1,000; lowered and
listed in `limits.clamped` above the cap; accepted with any request, it acts only with `"paths"`);
`require_support` (`"label_intersection"`, the default, or `"record_verified"`, §4.1). The checks are steps 8 and 11
of §5.

**What is searched.** Each searched orientation anchors on its own window (§7.7); the anchors are counted first
(`counts.anchors`, as §7.7). The extension is admitted when the anchors are `exact` and at most `max_anchors`
(the first admission); it then extends every anchor one base at a time along the graph's outgoing edges, keeping
only the bases the pattern allows at the next position (for a peptide, its codon automaton, the codon's state
recomputed at the k boundary from the anchor's bases, §12.2), by depth-first search to L bases. Every complete
walk is a **path**: n = L − k + 1 k-mers of the graph spelling an instance of the oriented pattern. A path need not
lie in one record. Each outgoing edge examined is one step (`work.extension_edges`, inside `steps`, under the
request's one `max_steps`); the work time is read at every 4,096 steps as in discovery (§7.6) and, since the
review GPT-3 (§18), before each anchor's extension, every 64 anchors listed and before every 64th node the search
expands, so that a late clock stops the extension on time (`stop {extension, time}`, `extension: "stopped"`)
rather than after it. Beside the steps (no step charges them): `work.extension_anchors`, the anchors whose
extension began, each spelled once (k − 1 BOSS steps) before its search, and `work.extension_branches`, the nodes
where the search branched (§8.7). The paths are
released when they are `exact` and at most `max_paths` in `all_or_count` (the second admission), the first
`max_paths` in `partial`.

**`counts.paths`** (a `paths_count`, §8.6): `{value, relation, unit: "paths"}` with `by_strand` (BASIC) or
`by_orientation` (one count per orientation searched), `candidates_examined` (the branches the extension
entered, prefixes of k + 1 to L bases, complete paths included: work, not a count of the pattern) and `extension`:

| `extension` | when | relation of `counts.paths` |
|---|---|---|
| `completed` | every anchor was extended | `exact` |
| `no_anchors` | the anchors are `exact` 0 | `exact` 0 |
| `stopped` | a stop in the extension: `max_steps`, `time`, or `max_paths` with `stop_at_threshold` | `at_least` (the paths completed before the stop) |
| `not_started` | a stop before the extension (in discovery or a deferred scan) | `unknown` |
| `not_admitted` | the anchors `exact` and more than `max_anchors` | `unknown` |

The note `paths_later_increment` is absent. `timing.extension_ms` states the extension's time.

**Stops** (§7.6): the phase `extension` with reasons `max_steps`, `time` and `max_paths` (`stop_at_threshold`: the
extension stops once more than `max_paths` paths are complete); `stop_at_threshold` on the anchors stops discovery
with `max_anchors` as before. First stop wins; a stop is sticky for the later patterns as in §7.6.

**`withheld` and `cut`.** `all_or_count` releases every path or none: `anchors_above_threshold` (not admitted),
`count_above_threshold` (the paths `exact` and more than `max_paths`), `threshold_crossed` (a `stop_at_threshold`
stop, on the anchors or on the paths), `discovery_budget` (`max_steps` in any phase, the extension included),
`deadline` (the work time, anywhere), and with labels the reasons of §14.6. `partial` withholds only
`anchors_above_threshold` (nothing was extended); its `cut.reason` is `max_paths` (more paths than `max_paths`: the
first `max_paths` in answer order), `max_steps` (a stop in the extension: the paths completed before it, a prefix of
the answer order; a stop in discovery: none), `max_anchors` (`stop_at_threshold` in discovery: nothing extended,
`returned` 0), `time` or `max_memory`. `retrieval_complete` is `true` iff every path was released (the extension
`completed`, the paths `exact` and all returned; an `exact` 0 included) and, with labels, everything below was.

**A path result** (`path_result`; results in answer order: anchor node ascending, then orientation, then the
sequence, A < C < G < T):

<!-- schema: path_result -->
| field | type | meaning |
|---|---|---|
| `sequence` | string | the L bases the path spells (an instance of P for `+`, `=`, `forward`, `palindromic`; of rc(P) for `-`, `reverse`) |
| `anchor_kmer` | string | its first k bases, as the graph spells its anchor node |
| `instance` | string | the matched bases: equal to `sequence` |
| `offset` | integer | 0 |
| `strand` | `"+"` \| `"-"` \| `"="` | BASIC graphs |
| `orientation` | `"forward"` \| `"reverse"` \| `"palindromic"` | CANONICAL and PRIMARY graphs, instead of `strand` |
| `nodes` | list of integers | the n = L − k + 1 node ids of its k-mers in order, `nodes[0]` the anchor (§7.10; wrapper ids on a wrapped PRIMARY graph) |
| `rows` | list of integers \| null | the annotation row of each k-mer, as a context's `row` (§7.10); `null` for a k-mer without one |
| `support` | `"record_verified"` \| `"label_intersection"` \| `"mixed"` \| null | `labels: "all"`: its listed labels' support: all `record_verified`, none, or both; `null` when it lists no label or its labels were not built |
| `labels_status` | string | likewise: `complete`, `truncated`, `refused`, `not_read` (the worst of its rows': refused, then not_read, then truncated) or `output_budget` (read, not built) |
| `labels_total` | integer \| null | likewise: the labels on every k-mer of the path (the intersection of its rows' labels), when every row was read completely; `null` otherwise (a truncated row leaves the intersection partly known). With `require_support` the unverified ones are counted too |
| `labels` | list \| null | likewise: the labels on every k-mer of the path (with `require_support: "record_verified"` the verified ones only), in label order, each a `label` (§14.5) whose `support` is `record_verified` or `label_intersection`; `null` unless read and built |
| `labels_excluded_unverified` | integer \| null | with `require_support: "record_verified"` only: the labels of the path left out for not being verified, their true number; an integer only when every row of the path was read completely (as `labels_total`) and every label carrying it was verified or refuted; `null` otherwise: its labels not read, a row truncated (the labels cut are unknown), or its verification not done (a placement stopped or a placement read refused leaves a label neither verified nor refuted) |

Without labels a path has the first seven fields (`strand` or `orientation` once). Never `kmer`, `node`, `row`:
`kmer` keeps its meaning, the k-mer of a context.

**Support** (design §4.3, owner decision #14). A label of a path is
- `label_intersection`: the label annotates every k-mer of the path; nothing says one record holds it (a path
  may join k-mers of different records, or of different places of one record);
- `record_verified`: there is a column coordinate c of the path's first k-mer such that c + i is a coordinate of
  its i-th k-mer for every i < n, c maps (record mapping first) to (`seq_id`, local), and local + n − 1 is below
  the record's k-mer count: one record of the label holds the whole path there. A chain of consecutive coordinates
  that crosses into the next record of the column is not one. Possible only with `record` placement (BASIC,
  coordinates, the `.seqs` mapping) and `output.occurrences` not `false`; the capabilities' `support` says whether
  the index can give it.

A label of a path (`label`, §14.5): `column`, `support`; with `record` placement `occurrences` (the label's placed
occurrences of the whole path, unit `placed_occurrences`, `exact`: 0 when its coordinates show none inside one
record, then `support` is `label_intersection`; `unknown` when a row of the path was not placed) and
`occurrence_list` (`occurrence` objects: `seq_id`, `record`, `strand` (the path's), `nt_coords` "start-end" over the
L bases, 1-based, `nt_length`; ordered by (`seq_id`, start)); with `global` placement `occurrence_list` holds the
chains (`occurrence_global`: `kmer_coord` of the chain's first k-mer, `offset` 0, `strand`), nothing verified, every
label `label_intersection`; `occurrence_list` is `null` when not placed and absent for `none`, `none_canonical`,
`not_requested` (note `label_intersection_only`: every label `label_intersection`).

**`require_support: "record_verified"`** lists the verified labels only: in each path's `labels`, in `by_label`
and in the counts; each path counts the others in `labels_excluded_unverified` (an integer only when it is the
true number, else `null`, above) and the entry in `labels_excluded_unverified` (a count, unit `labels`: the labels
carrying some returned path but verified on none; `exact` or `unknown`; `exact` only when every path's is an
integer and every path was returned). It acts on the labels of paths only (entries of L ≤ k keep `support: "kmer"`). On an index
whose best support is not `record_verified` it is refused (400 `support_unavailable`) rather than answered in the
weaker mode; with `output.occurrences: false` it is 400 `invalid_request`; in an answer that reads no labels it is
stated as `annotation_not_read`.

**The entry with labels** (`labels: "all"` in a retrieval mode): `placement`, `annotation`, `rows_refused`,
`anchors_truncated` (the rows of the paths' k-mers, each named by its k-mer), `labels_cut`, `occurrences_cut` as for
contexts (§14); `by_label` lists `by_label_paths` objects, ordered by (paths desc, column asc):

<!-- schema: by_label_paths -->
| field | type | meaning |
|---|---|---|
| `graph` | string \| null | as for contexts (§14.5) |
| `column` | string | the label |
| `paths` | count | the returned paths listing it, unit `paths` |
| `paths_record_verified` | count | of them, those it verifies, unit `paths`; `unknown` without `record` placement |
| `occurrences` | count | its deduplicated placed occurrences of the whole path over the returned paths, unit `placed_occurrences`; `unknown` without `record` placement |

`counts.labels` (a `labels_count`, §8.6) carries `by_support`: `record_verified` (the labels verified on at least
one returned path; `unknown` without `record` placement) and `label_intersection` (the listed labels verified on
none); `counts.occurrences` the deduplicated `record_verified` occurrences summed over the labels.
- `counts.labels` is `exact` iff every path was returned, every row of every path read completely and every
  path's list made (with `require_support`, also every carried label verified or refuted); `at_least` otherwise
  (`unknown` when nothing was read or returned). `by_support` is `exact` iff the same holds with the verification
  complete; otherwise `record_verified` is `at_least` and `label_intersection` `unknown`. `by_label`'s `paths` as
  `counts.labels`, its `paths_record_verified` `exact` iff the verification is; its `occurrences` as
  `counts.occurrences` for contexts (§14.6).
- `all_or_count` withholds as for contexts (§14.6: `deadline`, `annotation_budget`, `output_budget`,
  `anchor_labels_truncated`); `partial` states what is missing as there.
- Note `label_intersection_only` (after increment 3's notes) for placement `none`, `none_canonical` or
  `not_requested`; `global` keeps `record_bounds_unknown`.

**Memory and the deadline** (labels: §14.4's account). A released path costs 512 + 2k + 3L + 192n bytes, charged
before its result object is built (`partial`: at most half of the account; the first that does not fit ends the
list: `cut: max_memory` and `stop {output, max_memory}`, `all_or_count` withheld `output_budget`); a path's label
list costs 32 + 8 per label (when it does not fit: `stop {output, max_memory}`, that path and the later ones
`output_budget`). Rows are read one per read, each once per step, work units as for contexts. Not in the account,
as the label-free descriptors (§7.6, "Memory"): the paths the engine retains during the extension (at most
`max_paths`, O(L) each) before their release; without labels the route has no account (bounded by `max_paths`
× O(L); a pattern has no length cap, §4.2). The work time is read before each annotation read and between its
chunks, before each path's label list (`stop {output, time}`), before each path's verification (`stop {placement,
time}`: the later paths' labels unverified) and before each path's labels are built (`stop {output, time}`); since
the review GPT-3 (§18) also inside that work, at least every 4,096 of its units and before them: in the label
lists' intersection (`stop {output, time}`), in a path's verification (`stop {placement, time}`: that path and the
later ones unverified) and in the occurrences made for a path's output, in pieces of at most 4,096 (`stop {output,
time}`).

**The verification** (review GPT-3, finding 1, §18). The label lists of a path are intersected in place, from the
shortest row list (each row's list sorted once, not copied per path; `timing.label_intersection_ms`). Each (path,
label) is then verified once: its chains are the intersection of the k-mers' coordinate lists, each shifted by
its k-mer's position, found by a leapfrog join that starts from the shortest list and seeks by galloping;
consecutive chains are extended as one run by a density check (a homopolymer's are one run, a few units per
k-mer), and each run is placed record by record (one record mapping per record it crosses: `record`, a chain
whose first k-mer is in one record and whose last is past that record's k-mers is cut, §4.3; `global`, the run
of chains as it is). Its units are `work.verification_steps` (§8.7), its time `timing.verification_ms`. The runs
of occurrences are kept for the output, charged before they are held (32 bytes for a label with occurrences, 32
per run), so the output does not join the chains again; when the account cannot hold a path's runs, `stop
{output, max_memory}`: that path and the later ones are not output (`labels_status: "output_budget"`), their
labels still verified. When the label lists stop (time or memory) at a path, the paths before it are answered
`labels_status: "output_budget"` with `labels: null` (an earlier build said `complete` with `labels: []` and
`labels_total` above 0).

**Limits echo** (answers to `long_search: "paths"` only): `limits.long_search` `"paths"`, `limits.max_paths`; and
`limits.require_support` in those that read labels.

**Capabilities**: `long_search` `["anchors", "paths"]`, `default_long_search` `"anchors"`, `caps.max_paths`
(1,000, `--pattern-max-paths`), `caps_rule` naming `max_paths` among the clamped caps (the rule is the one above; a reference since the
owner's decision P9, §18);
`long_patterns` stays `"anchors_counted"`; `support` states the best support of a path on the index.

### 12.2 Increment 5, served: peptides (`protein`)

Served by this build (`src/graph/alignment/genetic_code.{hpp,cpp}`, the codon automaton of
`src/graph/alignment/pattern_search.cpp`, the route's `protein` kind; design §6, owner decision #15) as an addition
to version 1. A request without a `protein` pattern is answered as before (§17); `genetic_code` acts on protein
patterns only.

**The kind.** `patterns[i].protein`: a peptide over the 20 amino acids A C D E F G H I K L M N P Q R S T V W Y,
the ambiguity codes **X** (any residue: every codon that is not a stop), **B** (D or N), **Z** (E or Q) and **J** (I
or L), and the stop **`*`** (the owner's decision #19 of 2026-10-08: a stop codon of the genetic code at that
position), any case (capabilities `protein_residues`). A peptide of m residues is a pattern of L = 3m bases.

**The genetic code.** `genetic_code`: the NCBI translation table (gc.prt version 4.6) the request's peptides are
read in, one of the capabilities' `genetic_codes` — 1–6, 9–16, 21–33 — default 1 (the standard code); another
integer is 400 `genetic_code_unknown`, a value that is not an integer 400 `invalid_request`. A residue's codons are
NCBI's `ncbieaa` column of the table. Tables 27, 28 and 31 list some codons both as a residue and as a stop in
context (27: TGA W; 28: TAA, TAG Q and TGA W; 31: TAA, TAG E): they code their residue here and can match as it,
which a translation that ends at them in context would not show.

**What is searched: the codon automaton.** The instances of a peptide are exactly the codon strings c1 … cm with
each ci a codon of residue i in the table: no superset (Leu TTR|CTN, Ser TCN|AGY, Arg CGN|AGR are exact, not the
per-position union), a stop codon only where the peptide has `*` (X excludes the stops too). The bases allowed at a position depend on
the bases already spelled in its codon; in discovery they are read from the range itself, in the extension from
the path spelled so far, so the automaton's state is recomputed at the k boundary from the anchor's bases.

**Strands.** As for dna (§7.3): `forward` searches P, `reverse` its reverse-complemented automaton rc(P) (the
residues in reverse order, each codon reverse-complemented: GCN becomes NGC), on BASIC strands, elsewhere
orientations. A peptide is palindromic when its codon sets equal their mirrored reverse complements: only runs of X
in tables 27, 28 and 31 (the tables without a stop codon) are, there with `*` at mirrored positions too (such a
peptide has no instance, below). In every table with stop codons a stop codon is never the reverse complement of a
residue's codon nor of a stop codon, so `*` makes no peptide palindromic.

**Counts, results and scopes** keep their meaning: offsets, `length`, scopes and instances are in bases.
`length` is 3m and `residues` m; `instance` (`kmer[offset, offset + L)`, or a path's `sequence`) names the codons
matched; `counts` count graph contexts, anchors and paths as for dna. A peptide of at most k / 3 residues (10 at
k = 31) is searched within one k-mer (`suffix`, `any_offset`); a longer one is a pattern longer than k: anchors only
(§7.7), its paths with `long_search: "paths"` (§12.1). Labels (`output.labels: "all"`) as for dna.

**Bits and the floor.** `information_bits` is exact: 2 per base less log2 of the distinct strings the peptide's
codons spell over the positions counted, residue by residue — for the whole peptide Σ log2(64 / |codons_i|) =
6m − log2(the codon strings it admits). `anchor_information_bits` and `min_anchor_information_bits` are exact for
windows that cut a codon. The floor (§7.8) gates on them as for dna; an exact peptide (every residue one codon) is
exempt in `suffix` scope.

**The stop `*`** (the owner's decisions #19 and #21 of 2026-10-08). `*` admits the codons the table's `ncbieaa`
column marks `*`: TAA, TAG and TGA in the standard code, TAA, TAG, AGA and AGG in table 2, and so on; X never
admits one. A peptide ending in `*` finds the stop codon after a coding sequence (fixture `peptide_stop`:
NDM-1's last 9 residues and `*`, table 1, 4 contexts whose instances end in TGA). Tables 27, 28 and 31 have no
unconditional stop codon: their context stops (27: TGA; 28: TAA, TAG, TGA; 31: TAA, TAG) code their residue
(decision #21, above) and are not matched by `*`, so there `*` admits no codon and a peptide holding it has no
instance. Such a peptide is answered, never refused, and never as a silent 0: every entry of it carries the note
`no_stop_codon` (§8.10). With L ≤ k it is answered `exact` 0 in every count without a search: no step, no budget
read and no information-floor refusal (§7.6: also after a sticky stop, `stop` `null`, `determinism: "full"`;
fixtures `peptide_no_stop_codon`, `peptide_no_stop_codon_after_stop`). With L > k it is searched as any other
long pattern: its anchors are its anchor windows' instances, which need not reach the `*` (they are counted, and
may be more than 0), and its paths, with `long_search: "paths"`, are `exact` 0 (the extension finds none). Its
information bits count a residue without codons as one exact codon (2 bits per base it covers), so that they stay
finite; a long one's anchor windows are gated by the floor as any other's. A pattern after it in the request
spends the budget as before.

**What is refused.**
- In the pattern's slot (§8.9): `bad_alphabet` for an empty text or a character that is not a residue (U,
  selenocysteine, and O, pyrrolysine, included; the first such is named); `information_below_floor` and
  `scope_unsupported` as for dna (the slot keeps `residues` and `genetic_code`). The slot error
  `stop_unsupported`, which `4596bb3b` gave a peptide holding `*`, is gone (§18).
- The request: `genetic_code_unknown` (above).

**Costs, stated.** A leading X is searched (the shortcut that skips a leading N run of a dna or iupac pattern,
§7.8, does not apply: X excludes the stops), and a run of X inside a window branches about as an N run of three
times its length (61 of the 64 codons in the standard code): §7.8's cost rule applies to it, whatever the bits.

**The entry** of a peptide: `kind: "protein"`, `pattern` (its residues, upper case), `length` (3m), `residues`
(m), `genetic_code` (the table used), and every other field as for dna.

**Capabilities**: `kinds` lists `"protein"` (`kinds_later_increment` is `[]`), `protein_residues` (the 24
letters and `*`), `genetic_codes`, `default_genetic_code` (1), `protein_rule` (prose).

## 13. Contract deltas against the design's draft (resolved)

The design's §7 is a draft of the full contract; milestone 1 built the following, which version 1 freezes.

| design | built and frozen | why |
|---|---|---|
| §7.1: `output.labels` defaults to `all` | the default is the capabilities' `default_projection`, frozen at `"none"` for version 1 (a default of `"all"` would raise the version) | `all` was not served in milestone 1 (the owner's instruction), and an omitted field must keep its meaning (the owner's decision of 2026-10-07). Stated in every answer (`output`) and in the capabilities; clients send it explicitly (§1) |
| §7.2: `by_offset`, `by_strand` and `suffix` as plain numbers in the example | full counts `{value, relation, unit}` | every count carries its relation (design §3): a per-strand count can be `exact` while the total is `at_least`, a per-offset count `exact` while the total is `bounds`; after a stop in discovery every per-offset count is `at_least` (§7.4) |
| §5.2, §5.3: `withheld: <reason>`, `stop: {phase, shard}` | `withheld: {reason}`, `stop: {phase, reason}` | objects leave room for the shard (milestone 6); the stop's reason is what an agent turns |
| §7.2: no `returned`, no `cut` | `returned` and `cut: {reason}` in the retrieval modes | `partial` states every cut (design §5.2) |
| §7.2: identity per result (`graph`, `index_fp`, `release`) | one top-level `index` | single graph; per-result identity arrives with milestone 6 |
| §7.2: `placement` per entry | in the entry only with `output.labels: "all"` (increment 3, §14.3); the capabilities' `placement` states what this index gives | nothing is placed without labels |
| §7.2: `work.annotation_rows`, `work.memory_bytes`; `timing.phase1_ms`, … | `work: {ranges_visited, mask_scans, steps}`, with `output.labels: "all"` also `annotation_rows`, `annotation_units`, `memory_bytes` (increment 3); `timing: {elapsed_ms}`, with labels also `label_discovery_ms`, `placement_ms` | no annotation read and no memory account on the label-free path (its memory is bounded by the caps, §7.6) |
| §7.2: `strands: ["+", "-"]` on every index | strand symbols on BASIC, orientation names elsewhere | no strand is known on canonical indexes (design §3) |
| §5.5: contexts in BOSS edge order, then offset | then orientation | an IUPAC pattern and its reverse complement can share one (k-mer, offset) |
| §7.1: whole-request 400s named generically | codes `invalid_request`, `later_increment`, `resident_only`, `mask_required` (retired, §18), the graph reasons; 503 `deadline` | one envelope `{error, code}` |
| §4: the mask is required (`mask_required`) | served without it (`counting: "upper_bound"`): counts `bounds` with an `estimate` where they cannot be proven, lists exact (§7.4) | the owner's decision #16 of 2026-10-08: some graphs have masks, some not; the mask is derived data outside `index_fp` (#17); §18 |
| §4.1: `scope_unsupported: suffix_on_primary` | slot code `scope_unsupported` | the message names `any_offset` |
| §7.3: the block's fields | as §10.2, plus `available`, `unavailable_reason`, `default_*`, `*_later_increment`, `caps_rule`, `long_patterns` | a client gates on what is served now and sees what is coming |
| §4: `mask: file \| built_at_load \| absent` | as designed: `file` \| `absent` at `b570800d`, `built_at_load` added by milestone 1b | `--pattern-build-mask` landed with 1b, no field changed |
| §5.3: `max_steps` per shard | per request (one shard) | single graph; per-shard shares with milestone 6 |
| §5.3: the floor exempts "an exact `dna` `suffix` count" | an exact pattern of either kind (every position one base: an `iupac` string without ambiguity codes too) in `suffix` scope (§7.8) | the cost is the same one range whatever the kind |
| §3, §5.3: the anchor window [0, k) of a pattern longer than k | each searched orientation's window, the least informative one gating and stated in `min_anchor_information_bits` (`anchor_information_bits` stays P[0, k)'s; §7.7) | the reverse orientation anchors on P's last k positions (review of 2026-10-07) |
| §7.2 (increment 3): a label's placed occurrences in `occurrences` | `occurrence_list` (the list) beside `occurrences` (its count `{value, relation, unit}`, as in `by_label`) | one name, one type: `occurrences` is a count everywhere in the answer |
| §7.2 (increment 3): `column` and `properties` by `get_label_as_json` | `column` is the annotation column's label as stored; no `properties` | the label is not parsed (a label with `;` or `=` stays whole); a later version may add `properties` |
| §7.2 (increment 3): `by_label` `contexts`, `contexts_suffix` as numbers | counts `{value, relation, unit: graph_contexts}` | a label's count is `at_least` when a row was cut or refused |
| §7.2 (increment 3): a result's labels are its labels | each result states `labels_status` and `labels_total` | a truncated, refused or unread row is told apart from a row without labels (`labels: null` unless read) |
| §5.3, §7.2 (increment 3): `memory.stop` with the phase | `stop: {phase, reason: max_memory}` (phases `label_discovery`, `placement`, `output`) and `work.memory_bytes` | the one `stop` object of the entry |
| §5.2 (increment 3): `max_annotation_work` gives `rows_refused` | a work stop is `stop: {phase, reason: max_annotation_work}`, the unread rows `labels_status: "not_read"`; `rows_refused` holds the rows the memory account refused, each with what it needed | a refused row is a property of the row; a work stop is not |
| §5.6 (increment 5b): a filtered answer marks `absence_scope: …filtered` | `absence_filter: "predicate"` beside `absence_scope`, which keeps its closed values (§19.11) | `absence_scope` is a closed enumeration (§1) |
| §5.6, §7.1 (increment 5b): `predicate_only` is the default projection with a predicate | the default stays `"none"`; a client sends `"predicate_only"` (owner decision P2) | `default_projection` is frozen for version 1 (§1) |
| §5.6 (increment 5b): a predicate is evaluated on the context's row | with `predicate_strands: "either"` (the default) on a BASIC graph also on its reverse complement's row, the one or the other, never a mix; `"context"` on request (owner decision P11, §19.5) | a record holding the motif on its other strand annotates the reverse complement: `none(B)` would pass contexts B carries |
| §5.6, §7.2 (increment 5b): `tested` and `selected` inside `selection`; a stopped selection `at_least` | `counts.tested` and `counts.selected` (counts with relations), `selection` {pass, support, access} per entry, the request's parts in the top-level `predicate` block; a stopped selection `bounds` [S, S + R_upper − T] where the raw count has an upper bound (P13) | they are counts; an upper bound exists |
| §5.2, §5.6 (increment 5b): raw paths kept for a selection | a predicate selects patterns longer than k only on supported paths, a later increment (`long_search: "paths"` with a predicate: 400; P24) | an enumeration cut at `max_paths` can hold no supported walk |
| §4.3 (increment 3): placement `global` with `kmer_coord` + `offset` | as designed, and nothing deduplicated or counted (`counts.occurrences` unknown) | without record bounds, `kmer_coord + offset` of two records can coincide (records are concatenated), so equal sums are not one occurrence |

## 14. Increment 3: `output.labels: "all"` (labels and placement)

Served by this build (`src/cli/pattern_retrieval.cpp`, design §4.3, §5.2–§5.5, §7.2) as additions to version 1.
`projections` lists `"all"`; a client gates on it (§10.3).

### 14.1 What is read, and when

- Only in a retrieval mode, and only for a pattern whose release the engine did not withhold: a pattern whose
  count was not admitted (`withheld` for any reason of milestone 1) reads no annotation row (design §5.2:
  `work.annotation_rows` 0). Mode `count` reads none whatever the projection (note `annotation_not_read`).
- A pattern longer than k releases nothing (`paths_later_increment`, §7.7) and reads nothing; without anchors its
  empty answer is complete, with `counts.labels` exact 0. With `long_search: "paths"` (increment 4) its released
  paths are read instead: the rows of their k-mers, each path's labels and their support (§12.1).
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
  demand, one row per read (outside review GPT-2, findings 1-2). A row that does not
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
  context 512 + 2k (its descriptor and its result object), a dictionary label 192 + 2 × its name's length
  (its name in the dictionary and in the placement's copy of it; charged inside the read that names it), the rows' label lists and coordinates as the `DecodeBudget` charges
  them (held until the pattern's labels are built), a statement of a row (a `rows_refused` or
  `anchors_truncated` entry) 384 + k, a label of a result 256 + its name's length (its copy of the name), a
  `by_label` entry 512 + its name's length (per pattern; GPT-2 finding 3: every retained copy of a name is
  charged before it is built), a placed occurrence 256 + its record name's length
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
  - then `by_label`, all its entries, before any context's labels: when it does not fit, `stop {output,
    max_memory}` (unless an earlier stop is stated), `all_or_count` withholds (`output_budget`), `partial` returns
    every context read with `labels_status: "output_budget"` and `by_label: null`;
  - then the labels of each context, in answer order — when a context's do not fit, `stop {output,
    max_memory}`: `all_or_count` withholds (`output_budget`), `partial` returns it and the later ones with
    `labels_status: "output_budget"`.

  Since the review GPT-3 (finding 4, §18), in `partial` a context's (or path's) label holds, and the account and
  the estimate E of §7.6 are charged for, only the first `max_occurrences_per_label` of its placed occurrences;
  every occurrence still enters its label's deduplication set (64 each, repeated coordinates once), so the counts
  stay exact. Once the sets are complete, each list is cut to its set's first `max_occurrences_per_label` (§14.5)
  and what the cut takes is given back to the account and to E. A set is checked against the account as it grows
  (the decision the context's final charge would make, made earlier). Before, every occurrence was held and
  charged, and E counted the text of occurrences no list would show.

  `work.memory_bytes` is the account's peak so far. Apart from the forced names of an unbudgeted read (above),
  it never passes `max_memory_mb`.
- **Work** (`max_annotation_work`, the oracle's units, one budget per request): a row read costs 8, plus 1 per
  label of the row (its true total, also when truncated) in step 1, plus 1 per label and 1 per coordinate in
  step 2, plus the units of its row-diff dependency rows (8 per row and 1 per entry they store); a row read in
  both steps costs both. The reads take one row at a time and the budget is checked before each row: a row is
  read only while the units are below `max_annotation_work`, so the reads pass it by their last row's units at
  most (`annotation_units` < `max_annotation_work` + that row's units; GPT-2 finding 2: batches of growing size
  could pass it by many rows). A refused row costs what its read decoded: a read row's units when it was read and
  then refused for its demand or its names, 8 when its read itself did not fit; it counts in `annotation_units`,
  not in `annotation_rows` (GPT-2 finding 1: refused reads were not charged). The stop: `stop {label_discovery | placement, max_annotation_work}`, the rows not read `labels_status:
  "not_read"`, `all_or_count` withheld (`annotation_budget`); the request's later patterns read nothing.
- **The deadline.** The reads and the output of the labels are work (§7.6): the work time is checked before each
  read and, within a read, between its chunks (paced at `--traverse-chunk-target-ms`, as `/traverse`'s reads; an
  interrupted read returns nothing), and before the labels of each context are built for the answer, with the
  labels about to be built counted in the estimate E of §7.6; since the review GPT-3 (§18) also while the
  occurrences of a context are made, in pieces of at most 4,096, the clock read before each piece (`stop {output,
  time}`). A time stop of the reads: `stop {label_discovery |
  placement, time}`; of the output: `stop {output, time}` — each, like every `stop` of this section, only when no
  earlier stop of the pattern is stated (first stop wins, §7.6: after a `max_steps` stop of the engine, or a
  `max_annotation_work` stop of the reads, a time stop of the output shows only as `labels_status:
  "output_budget"` and `time_limited`). Either way `determinism: time_limited`, `all_or_count`
  withheld (`deadline`), and `partial` returns what was read, a context whose labels were read but not built
  with `labels_status: "output_budget"`; the later patterns answer as after any time stop (§7.6).
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
  came first, or `all_or_count` stopped reading at a refused row), `output_budget` (read, but the answer could
  not hold its labels: the memory account, or the time to write them; `stop` says `{output, max_memory}` or
  `{output, time}` unless an earlier stop is stated, and a time stop sets `determinism: time_limited` either
  way). `labels` is
  `null` unless the row was read and held, and then lists the context's labels in label order.
- **Label order** (design §5.5): (contexts desc, column asc), over the returned contexts — `by_label`'s order.
- **`partial`'s lists.** `max_labels`: `by_label` and every result list the first `max_labels` labels of the label
  order (`labels_cut`); the counts count all of them. `max_occurrences_per_label`: each listed label lists the
  first `max_occurrences_per_label` occurrences of its deduplicated union, in (`seq_id`, start, strand) order —
  every context lists those of them it holds, possibly none (`occurrences_cut`); the label's counts stay whole.
  `all_or_count` cuts nothing: all or nothing.
- **`by_label`**: one entry per label over the returned contexts; `null` when the results are withheld or, in
  `partial`, when the memory account could not hold it (§14.4), `[]` when none was found.

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
          "annotation_rows_distinct": 24, "annotation_units": 24440, "memory_bytes": 225828},
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
- **Answers to the other requests.** This is the regression check of increment 3, made before the changes of
  §16, which changed some of these answers within the semantic compatibility of §1 (a promise of meaning, not of
  bytes, across builds). At increment 3, a request valid at `bd44e597` was answered byte for byte as by the
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

## 16. Changed after the review of 2026-10-07

The critical review of milestones 1 and 1b (2026-10-07) found no false count in a configuration a host can
reach, but budgets that did not hold, limitations left unstated and sentences of this document that the code
did not keep. The fixes below change answers; the finding ids are the review's. This document was corrected
where it overstated the build (§§1–4, 6, 7.4–7.10, 8, 10, 11, 13, 14): every such sentence now says what the code
does. Unchanged: `/search`, `/align`, `/traverse`, `/resolve` without `bounds.time_budget_ms`, `column_labels`
and `/stats` (byte for byte apart from `timing`), the alignment, and every answer not named below.

**Contract version 1 stays: additions only.**
- **New field** `min_anchor_information_bits` (every answered entry and every `information_below_floor` or
  `scope_unsupported` slot; `null` for L ≤ k): the bits of the least informative searched anchor window, the one
  the information floor now gates on (§7.7, §7.8; X-GUARANTEES-01). `anchor_information_bits` keeps its meaning,
  the bits of P[0, k), whatever the strands (the owner's decision of 2026-10-07). The two differ only with
  `strands: "reverse"`, or `"both"` when P's last k positions carry fewer bits than its first k. Every fixture
  answer with a pattern entry carries the new field; no other value of theirs changed.
- `labels_status: "output_budget"` now also covers a time stop in the output of the labels (`stop {output,
  time}`), where it meant the memory account alone (X-EFFICIENCY-04); `stop` tells the two apart (unless an
  earlier stop is stated, §7.6; `determinism: time_limited` marks the time stop either way). The owner
  accepted this as within version 1 (2026-10-07).
- **New refusal reasons** `mask_invalid` and `alphabet_untested` (values of a refusal's `code` and of
  `unavailable_reason`, both extensible, §1; the owner's decisions of 2026-10-07): see "Answers that change".
- **A reserved request field** `long_search` (§4.4): 400 `later_increment`, any value; before, 400
  `invalid_request` (an unknown field). Paths for patterns longer than k will be opt-in through it (§12), so
  that no answer to a request without it changes when they arrive.
- **The compatibility promise is stated** (§1, §7.9; the owner's decision of 2026-10-07, after the outside
  review's GPT #1): version 1 keeps the meaning of every field and value and what every relation guarantees;
  work counters, stopping points and valid bounded values may differ between builds; the byte-for-byte promise
  holds within one build and effective configuration. A client accepts only the versions it implements (§1,
  §10.3; it no longer accepts "1 or higher").

**Answers that change:**
- **DNA5 graphs** (I26: X-TESTS-01, D1-03, T1-08, X-ORACLE-02, E2-05; the owner's decision): a `$ACGTN` graph
  is refused, `available: false` with `unavailable_reason: "alphabet_untested"` and every request 400
  `alphabet_untested`, until a DNA5 build passes the pattern tests; before, it was served, untested (§8.2). No
  deployed index is DNA5.
- **A mask with a valid W = `$` edge** (I17: E1-02, E3-02; the owner's decision): refused, `available: false`
  with `mask_invalid` and every request 400 `mask_invalid`, naming `metagraph transform --mask-dummy --force`;
  before, such a mask (written by `metagraph extend` on a masked graph) was trusted and could give a false
  `exact` count (§10.2). `metagraph extend` now rebuilds the mask of a masked graph, and removes a mask left
  beside its output when the extended graph has none. Masks written by `build`, `transform` or `--pattern-build-mask` are unaffected.
- **The floor for L > k** (X-GUARANTEES-01): with `strands` `both` (the default) or `reverse`, a pattern whose
  reverse window is below the floor is refused (`information_below_floor`, the message naming the reverse
  orientation's window, `min_anchor_information_bits` its bits); before, it ran, often to `max_steps` (P + N^31
  at k = 31). With `strands: "reverse"`, N^k + P is answered; before, it was refused. The refusal of a
  low-information reverse window stays (the owner's decision of 2026-10-07).
- **N runs** (X-EFFICIENCY-01): on `$ACGT` graphs a run of N at the start of a searched window is not searched
  (§7.8): `work.steps` and `ranges_visited` of such patterns fall (N^12 + `TACAGAGCGACG` on the mini: 69 steps,
  was 20,151,974); counts and results are the same when `exact`, but under a step or time stop the `at_least`
  and `bounds` values differ (more is found per step). The base searches run cheapest first: under a stop a
  degenerate pattern can now answer its cheap orientation `exact` and the other `at_least`, where the forward one
  was `at_least` and the reverse one `unknown`; `strands` still lists the plan's order.
- **`all_or_count`, the delivery to the route** (R1-01, X-EFFICIENCY-02, C1-02, D1-06, E4-02,
  X-GUARANTEES-03): when the work time passes while the released contexts are handed to the route, the answer is
  200 with the `exact` counts, `withheld: deadline`, `returned` 0, `results: []`, `stop {extraction, time}`,
  `determinism: time_limited`; before, a complete answer past the work time or a 503.
- **`partial`, the release** (E4-03): a time stop in the release comes within 64 contexts handed to the route;
  before, within 4,096 examined edges.
- **`partial`, what is kept** (E4-01, X-EFFICIENCY-03, X-GUARANTEES-02): only the descriptors that can be
  released (§7.6, "Memory"), and the release's set-up reads the clock: requests whose set-up overran the
  deadline (a 503) now answer 200, complete or `cut: time`; the peak memory of a broad pattern falls (AC at
  `max_contexts` 1 on the masked mini, floor 0: 1,660 MB to 58 MB). Lists and work are otherwise unchanged.
- **The time kept back for the answer** (X-EFFICIENCY-04, R1-04): the work stops earlier by the estimate E of
  §7.6. A retrieval request that buffered many results and used to end in a 503 now answers 200 with a time stop;
  the pattern it stopped, and every pattern after it, can report `at_least` or `unknown` counts where the work
  would have gone on to `exact`; `partial` says `cut: time`, `all_or_count` `withheld: deadline`. With
  `labels: "all"` a time stop of the reads can come earlier, and the output of the labels can stop by time
  (`stop {output, time}`; §14.4). In mode `count`, and for answers with few results, E is a few ms.
- **`stop_at_threshold` on an even-k wrapped PRIMARY graph** (E2-02): the running lower bound adds the two base
  searches less the possible palindromes, where it took the larger: the stop fires earlier, or at all. Answers
  that were `exact` above the threshold with no stop can now be `at_least` with `stop {discovery, max_contexts}`
  and `withheld: threshold_crossed` (`cut: max_contexts` in `partial`), with less work. No other graph mode
  changes.
- **`partial` on an even-k wrapped PRIMARY graph after a discovery stop** (E2-01): the palindromic stored k-mers
  only the reverse-complement search reached are released too: more contexts, so that the list holds every
  context the counts credit unless `max_contexts` cuts it.
- **A release that disagrees with its count** (T1-02, E2-04): `all_or_count` fails the request (§6) instead of
  answering `retrieval_complete: true` with fewer results, and so does `partial` instead of reporting the short
  release as `cut: max_contexts` (the owner's decision: an internal error, not a new
  `withheld` reason); seen only with a mask written by `metagraph extend` on a masked graph, which is now
  refused (`mask_invalid`, above).
- **A client that left, or a shutdown** (R2-02, X-CONCURRENCY-01, the owner's decision): nothing is written
  (§3); before, the request ran to its stop and wrote into the closed connection. A half-close counts as gone.
  Nor is an error written (outside review GPT-2, finding 4): `{` and `{"patterns":[]}` from a half-closed client
  were answered 400; now nothing. The other routes are unchanged.
- **Request parsing** (R1-02, R1-03, R2-01): bodies the server answered 200 — with a comment, a trailing comma,
  content after the value, or a duplicated member name (whose last value won silently: `max_steps` 1 then 100000
  gave 100000, unlisted) — and bodies nested deeper than 1,000 (a 400 without a code) are now 400
  `invalid_request`. No fixture's request changes. The CLI answers every request file: one that fails
  unexpectedly gets `{"error": …}` and the run goes on (it aborted before).
- **Time caps** (R1-05): `server_query` refuses to start with `--pattern-max-time-ms` above 899,000 (§3).
- **Messages** (prose, §1): `mask_required` adds "; the mask is read when the graph is loaded: restart the server
  once the .edgemask exists" (C1-06); the 503 `deadline` message states a fractional budget as given
  ("1000.5 ms", was "1000 ms"; T3-07); a reverse window below the floor is named (above); the
  `alphabet_unsupported` message says that the alphabet is not `$ACGT` and that `$ACGTN` is not served yet
  either; `alphabet_untested` and `mask_invalid` have their own (the latter naming `transform --mask-dummy
  --force` and the restart).
- **Capabilities** (both routes): `caps_rule` rewritten (R2-04, X-EFFICIENCY-04): it names the nine clamped
  caps, says that `max_patterns` and `min_information_bits` are not request fields' maxima, and states the rule
  of §7.6 with the delivery rates in force. The `/resolve` block's `rule` (on both capabilities routes) is
  rewritten too: on the row paths the explicit labels' hits are taken from each priming batch before the
  deadline is read, so a deadline stop keeps every row read (f8919958, T3-01/V1-01), and it says which invalid
  selections are still a 400 after a deadline stop (V1-02). No field is added or removed.
- **Server flags**: `--pattern-delivery-build-mbps` (10) and `--pattern-delivery-compress-mbps` (50) are new;
  `--pattern-finalize-ms` is now the floor of the time kept back (§4.5).
- **Fixtures** (§11): new `capabilities_representation_unsupported`,
  `traverse_capabilities_representation_unsupported`, `representation_unsupported`, `max_steps_bounds`,
  `max_steps_bounds_withheld`, `max_steps_then_time` (C1-03: `stop` keeping a `max_steps` stop, the release's
  time shown as `cut: time` and `time_limited`, and the sticky stop after it); every pattern entry gains
  `min_anchor_information_bits`; changed: the capabilities bodies (`caps_rule`), `mask_required` (the message),
  `deadline` and `deadline_partial` (the heavy pattern is now `GG` + N^19 + `TTGGCGATCT`, the old one having
  become cheap; at L = 31 its `by_offset` has one key), `README.md`, `index.json`; after the outside review
  (85614d30, GPT #11) new `labels_all_output_budget`, `labels_all_partial_exact_cut` and `labels_all_mixed_slots`
  (§11); `--check` now compares every
  `time_limited` entry after an answer's first, and blanks the `returned` and `results` of that first entry when
  the clock cut its release (`max_steps_then_time`).
- **After the outside review GPT-2** (labels "all" only; count mode, labels "none" and every label-free answer
  unchanged):
  - reads take one row at a time with the work budget checked before each row, and a refused row costs what its
    read decoded (findings 1-2): under a work budget fewer rows can be read when later rows are wider, refused
    rows add to `annotation_units`, and rows that used to be refused can be `not_read` behind a work stop; near
    the memory limit which row is refused first, and its `needed_bytes`, can differ;
  - every retained copy of a label name is charged, `by_label` first (finding 3): `work.memory_bytes` grows in
    every labelled answer (the ten labelled fixtures: e.g. `labels_all` 220,428 to 225,828), contexts can become
    `output_budget` earlier, `all_or_count` can withhold `output_budget` where it completed, and `partial` can
    answer `by_label: null`;
  - nothing is written to a client that left, errors included (finding 4, §3);
  - the fixture validator knows `mask_invalid` and `alphabet_untested` and checks its code lists against the
    sources (finding 6), and accepts `by_label: null` in `partial` (its recheck);
  - a new fixture `labels_all_rows_refused` (63 in all) shows a refused row and `cut: max_memory`.
- **No answer changes** (stated for completeness): the deleted per-offset completion rule was never in effect
  (T1-04: no offset was `exact` after a discovery stop before either); the dummy-edge mask is built faster with
  the same bytes (D1-02, M1-02); an `.edgemask` that exists but cannot be opened is now named in the log, by
  the loader's warning and by the start-up note (which said "no dummy-edge mask" and named a remedy `transform`
  refuses; M1-03).

## 17. Increments 4 and 5

Increment 4 (patterns longer than k as paths, §12.1) and increment 5 (peptides, §12.2), built in parallel and
served together (2026-10-08; owner decisions #12–#15). **Contract version 1 stays: additions only** (§1, the
owner's decision #12 of semantic compatibility): no field, value or count changes its meaning or its guarantee,
and both increments are opt-ins — a request that sets neither `long_search: "paths"` nor a `protein` pattern (nor
another field below) is answered as by the build of `44583b51`, byte for byte apart from `timing`. Checked against
that build's binary on a masked copy of the mini index (`transform --mask-dummy`), with and without
`--no-coord-mapping`: every `/pattern` request of the increment-3 review panel, every stored fixture request of
those servers but the two deadline fixtures, patterns longer than k in every mode, projection and strand setting,
with step and threshold stops, and refusals that name no new field — 544 comparisons (272 requests, plain and
gzip), all identical apart from `timing`; `/search`, `/align`, `/resolve`, `/traverse`, `column_labels` and `/stats`
— 42 comparisons, all identical apart from `timing`; both capabilities routes: only the additions below. The
fixture `max_steps_then_time` depends on the clock (the base binary answers it differently from run to run) and
was compared as `pattern_fixtures.py --check` compares it. The alignment golden gate (`run_gate.sh`,
`diff_gate.sh`) is IDENTICAL (1,038 of 1,038 files). The refusal of a pattern object without a kind keeps its
message; only one that names `protein` beside another kind gets the new message naming `protein`.

**Additions, each opt-in or additive:**
- **Request fields**: `long_search` (`"anchors"` | `"paths"`, default `"anchors"`), `max_paths` (≥ 0, default and
  cap `caps.max_paths`), `require_support` (`"label_intersection"` | `"record_verified"`), `genetic_code` (an NCBI
  table id, default 1); the pattern kind `protein`; `output.paths` accepted with either value (§4.1–§4.3).
- **Refusal codes** (§6): 400 `support_unavailable` (check 11 of §5), 400 `genetic_code_unknown` (step 8).
- **Slot error** (§8.9): `stop_unsupported` (retired by §18, the stop `*` being a residue; `4596bb3b` is the
  one build that answers it).
- **Entry fields**: `residues` and `genetic_code` (protein patterns); `labels_excluded_unverified` (paths with
  `require_support: "record_verified"`); `kind` value `"protein"`.
- **Counts**: `counts.paths` of a path search with `by_strand` | `by_orientation`, `candidates_examined`,
  `extension` (`no_anchors`, `not_started`, `not_admitted`, `stopped`, `completed`); `counts.labels.by_support`
  (`record_verified`, `label_intersection`) for the labels of paths (§8.6).
- **Results**: the path result (`sequence`, `anchor_kmer`, `instance`, `offset`, `strand` | `orientation`, `nodes`,
  `rows`; with labels `support` — `record_verified`, `label_intersection`, `mixed` or null —, `labels_status`,
  `labels_total`, `labels`, `labels_excluded_unverified`); a label of a path with `support` `record_verified` or
  `label_intersection`; `by_label` of paths (`paths`, `paths_record_verified`, `occurrences`).
- **Values**: `withheld.reason` `anchors_above_threshold`; `cut.reason` `max_paths`; `stop.phase` `extension`;
  `stop.reason` `max_paths`; note `label_intersection_only`.
- **Work and timing**: `work.extension_edges`, `timing.extension_ms` (path searches).
- **Limits**: `long_search`, `max_paths` (answers to `long_search: "paths"`), `require_support` (those that read
  labels); `limits.clamped` may name `max_paths`.
- **Capabilities** (both routes): `long_search`, `default_long_search`, `caps.max_paths`, `protein_residues`,
  `genetic_codes`, `default_genetic_code`, `protein_rule`; `kinds` gains `"protein"` and `kinds_later_increment`
  becomes `[]`; `caps_rule` names `max_paths` and states the rule of `long_search`. `long_patterns` stays
  `"anchors_counted"` (what a request without the option gets) and `default_projection` stays `"none"`.
- **Server flag**: `--pattern-max-paths` (1,000; `server_query` and `pattern`).

**Requests that were refused and are now answered** (the fields were refused by name, 400 `later_increment`, any
value, §4.4 before): `long_search`, `max_paths`, `require_support`, `genetic_code`, `patterns[i].protein` and
`output.paths: true`. A wrong value of the first four is now 400 `invalid_request` (`genetic_code` that is no table:
`genetic_code_unknown`). The fields still refused by name: `predicate`, `max_predicate_contexts`,
`max_predicate_work`, `graphs`, `budget_split`, `output.labels: "predicate_only"`.

**Stated limits** (each also in §12.1 and §12.2): the paths the engine retains during the extension (at most
`max_paths`, O(L) each) are not charged to the memory account, and without labels a path search has no account
(bounded by `max_paths`); a label of a path that no record verifies is `label_intersection`, never a record claim;
`record_verified` needs a BASIC index with coordinates and its record mapping; stops `*` in peptides are refused
(`stop_unsupported`; served since §18); a leading X is searched, not skipped; the context stops of tables 27, 28
and 31 match as their residue.

**Fixtures** (§11): new `paths`, `paths_count`, `paths_labels`, `paths_require_support`,
`paths_max_paths_partial`, `paths_count_above_threshold`, `paths_stop_at_max_paths`, `paths_anchors_above_threshold`,
`paths_max_steps`, `paths_global`, `paths_primary`, `support_unavailable`, `peptide`, `peptide_count`,
`peptide_paths`, `peptide_bad_residue`, `genetic_code_unknown` (80 in all); changed: the nine capabilities bodies
of the single-graph servers (the capabilities additions above), `README.md`, `index.json`. Every other stored body
is unchanged. The validator (`test_pattern_fixtures.py`) knows the new fields, values and codes (its tables
`paths_count`, `labels_count`, `by_support`, `path_result`, `by_label_paths`), checks a peptide's bits and
instances against its own copy of NCBI's genetic codes, and refuses answers that break the rules of §12.1 and
§12.2.

## 18. After `4596bb3b`: masks optional, derived sidecars, the stop residue

The owner's decisions of 2026-10-08: **#16** (masks optional: a graph without its dummy-edge mask is answered,
its counts upper bounds with an estimate, its lists exact), **#17** (the mask and the Bloom filter are derived
data of the graph, not part of `index_fp`), **#19** (`*` in peptides is a stop codon of the chosen table; X still
never matches a stop; no `stop_unsupported`; a table without an unconditional stop codon answers with a note) and
**#21** (tables 27, 28 and 31: their context stops match as their amino acid and never as `*`); #20
(`genetic_code_unknown`) stays. Built by the engine (`pattern_search.{hpp,cpp}`, `boss.{hpp,cpp}`,
`dbg_succinct.{hpp,cpp}`), the route (`src/cli/pattern.{hpp,cpp}`, the loader) and the identity code
(`src/cli/traverse.{hpp,cpp}`, `server.cpp`, `scripts/traversal/index_manifest.py`). **Contract version 1
stays: additions only** (§1). A graph without its mask used to be refused (400 `mask_required`): answering it is
an addition, and on a masked graph every answer to a request without `*` keeps its bytes apart from `timing`.
`*` was refused (`stop_unsupported`) only by `4596bb3b`, pushed to the branch on 2026-10-08 before these
decisions were built: a client that met it reads `stop_unsupported` as a slot error (§1). Identity is not a field of this route; its rule is the traversal contract's.

**Additions:**
- **Answer, on a graph without its mask only** (§7.4, §8.2): `index.counting` (`"upper_bound"`) and
  `index.dummy_fraction` (`{value, interval, samples, source: "sampled"}`, the schema `dummy_fraction`);
  `estimate` on every `bounds` count of the graph's units (contexts and their parts, anchors, paths); the relation
  `bounds` also without a stop; the notes `estimate_sampled_dummy_fraction` and `threshold_upper_bound`.
- **Capabilities** (both routes, every single-graph server): `counting` (`"exact"`, `"upper_bound"`, or `null`
  while loading or when not served) and `dummy_fraction` (the object, or `null` with a mask); `protein_residues`
  gains `"*"`; `protein_rule` rewritten, shorter. Nothing else. The rule of the counting is stated here, not as
  prose in the capabilities: the document a service's MCP tool returns whole has a ceiling of 32 KiB
  (`api/python/metagraph/traverse/mcp_tools.py`, `CAPABILITIES_MAX_BYTES`), and the fixture servers'
  `/traverse/capabilities` is 32,490 bytes of compact JSON with the mask and 32,596 without (32,593 at
  `4596bb3b`): 172 bytes are left.
- **Request**: `*` in a `protein` pattern (§4.2, §12.2).
- **Note** `no_stop_codon` (§8.10, §12.2).
- **Values reused, their when widened** (no new value of `withheld.reason`, `cut.reason`, `stop`, `extension` or
  a refusal code): `count_above_threshold` also for a `bounds` total with U > `max_contexts`;
  `threshold_crossed` also for a stop on the running upper bound; `anchors_above_threshold` and `extension:
  "not_admitted"` also for `bounds` anchors with U > `max_anchors` (§7.5).

**Retired and removed:**
- `mask_required`, as a 400 code and as an `unavailable_reason`: retired. No configuration of this build answers
  it (its message is gone from the source); it stays in the vocabulary of version 1, since builds before the
  decision answer it, and a client keeps handling it (§1, §6, §10.2).
- `stop_unsupported` (slot code): added by `4596bb3b` (pushed on 2026-10-08, the one build that answers it),
  retired: no source of this build writes it, and a client keeps reading it as a slot error from that build (§1).
  A `protein` pattern holding `*` is answered.

**Answers that change:**
- **A graph without its mask** (`mask: absent`; as `build/mini_refseq` and refseq33m-experimental's graph stand
  today, decision #18 masking the latter later): every request was 400 `mask_required`, the capabilities
  `available: false`; now answered, the capabilities `available: true`, `counting: "upper_bound"`,
  `dummy_fraction`. Compared with the same graph masked: every `exact` count equal, every `bounds` holding the
  masked count, every complete list equal (the engine's unmasked twins, 360 random patterns; the route's panels;
  `integration_tests/test_pattern.py`). Where it is weaker, stated: `all_or_count` admits on U, so a pattern
  whose contexts fit `max_contexts` can be withheld (`unmasked_threshold_upper_bound`: 17 contexts, U = 32);
  `stop_at_threshold` stops on the running U, earlier (`partial` can return fewer contexts, even none, its cut
  stated); mode `count` states `bounds` where a retrieval of the same request states `exact`.
- **Protein patterns holding `*`**: refused `stop_unsupported` by `4596bb3b`; now answered with the stop codons of
  the table, or, in tables 27, 28 and 31, `exact` 0 (L ≤ k, no search) or searched anchors and no path (L > k),
  with the note `no_stop_codon`; or refused `information_below_floor` as any pattern with too few bits. The
  fixture `peptide_bad_residue` keeps its request; its second slot (`MELPNIMHPV*`, 33 bases) is now answered:
  scope `long`, anchors `exact` 0, paths `exact` 0.
- **A peptide without instances of L ≤ k after a sticky stop** (§7.6): `exact` 0, `stop` `null`,
  `determinism: "full"` (a new request, so no earlier answer changes).
- **Masked graphs: nothing else.** Checked against the binary of `4596bb3b`: a CLI panel of 45 requests on a
  graph of 100 transcripts at k = 15, BASIC and CANONICAL, every mode, scope and strand setting, labels `none` and
  `all`, `long_search: "paths"`, `stop_at_threshold`, `max_steps` 500 and peptides — identical apart from `timing`
  but for the `*` slots; a panel of 37 requests on a masked copy of the mini (counts, both retrieval modes, labels
  `all`, paths, thresholds, step stops, peptides, refusals) — identical apart from `timing`; every stored fixture
  of a masked server but the capabilities and `peptide_bad_residue` — unchanged (`pattern_fixtures.py --check`);
  `/traverse`, `/resolve`, `/stats` and `/traverse/capabilities` without its `pattern` block on a masked mini copy
  — identical apart from `timing` and `attempts.server_instance`; the alignment golden gate (`run_gate.sh`,
  `diff_gate.sh`) — IDENTICAL, 1,038 of 1,038 files.

**Identity (decision #17;** `SPEC-labeled-traversal-core.md`, "The index identity"**):**
- The mask (`.edgemask`) and the Bloom filter (`.bloom`) are derived data of the graph, not part of `index_fp`:
  adding, removing or rebuilding one leaves `index_fp` and `index_meta_fp` unchanged (shown on a copy of the mini:
  `transform --mask-dummy` beside a deployed graph, the same `index_fp` before and after, the stored answers and
  graphlets replayed). An exact count is the same with and without the mask (§10.2); the mask only changes which
  counts are exact.
- A manifest that lists one is refused at start-up, naming the entry and the rule; `index_manifest.py` never
  writes them (`--extra` refuses them, `--verify` flags them). Manifests written between the review of pass 5 and
  this decision that list a mask or a Bloom filter are written again without them (their `index_fp` changes
  once).
- `traverse --index-inventory` and `index_manifest.py --inventory` list the derived files apart (`derived`,
  `derived_rule`); `files` holds the identity files only.
- Limitation, stated: no fingerprint covers derived data; a foreign Bloom filter of the same k and mode can hide
  k-mers under an unchanged `index_fp`.

**Server and CLI** (no new flag): f is sampled once per graph in the loading thread, and logged ("Dummy fraction
sampled for the pattern search in … s …: f, its interval, the seed"); the start-up line of a graph without a
mask says it counts upper bounds with estimates and names the two ways to exact counts (`transform
--mask-dummy`, `--pattern-build-mask`, whose help now says it is for exact counts); `transform --mask-dummy` no
longer says that a manifest must list the mask; a server with a checked manifest logs the derived files loaded
beside the graph.

**Stated limitations** (each also where its rule is):
- The estimate is not a bound: it assumes the source dummies as frequent among a count's candidates as in the
  graph, which a pattern at a record start, held by the dummies before it, breaks (§7.4); f is sampled from
  10,000 entries, its 95% interval stated.
- Without the mask, `partial` keeps every discovered range not spelled whole (24 bytes each, at most one per step)
  instead of about `max_contexts` of them (§7.6).
- On an even-k wrapped PRIMARY graph without its mask, `stop_at_threshold` can fire late or not at all (its
  running U leaves out the ranges whose palindrome scan is pending); the admissions after discovery compare the
  final U (§7.5).
- A long peptide with `*` in tables 27, 28, 31 is still searched for its anchors and gated by the floor on its
  anchor windows, although it has no path (§12.2).
- The `bad_alphabet` message of a protein pattern does not name `*` among the residues (kept byte for byte, §8.9).
- Derived data is fingerprinted nowhere (above).
- The capabilities document has 172 bytes left under the MCP tool's ceiling (above).

**Fixtures** (§11; 88 in all): new `unmasked_count`, `unmasked_labels_all`, `unmasked_threshold_upper_bound`,
`unmasked_stop_at_threshold`, `unmasked_partial`, `unmasked_paths`, `peptide_stop`, `peptide_no_stop_codon` and
`peptide_no_stop_codon_after_stop`; removed `mask_required` (its server now answers); changed: the nine
capabilities bodies of the single-graph servers (`counting`, `dummy_fraction`, `protein_residues`,
`protein_rule`; the unmasked pair also `available: true`, `unavailable_reason: null`) and `peptide_bad_residue`
(its `*` slot, above); `README.md`, `index.json`. Every other stored body is unchanged. The validator
(`test_pattern_fixtures.py`) knows `counting`, `dummy_fraction` (its interval holding its value), `estimate` (its
formula against `index.dummy_fraction.value`, present exactly on the `bounds` counts of a graph without its
mask), the three notes and their order, the codons of `*` in every table (its own copy of NCBI's tables),
peptides without instances, the sticky-stop exception and the threshold notes; `mask_required` is in its
`RETIRED` (named in the SPEC, written by no source, held by no fixture), and each new rule refuses a mutated body
that breaks it (`test_unmasked_and_stop_rules_refuse_what_v1_never_answers`).

**Tests with independent oracles** (beside the fixtures): the engine's `PatternUnmasked` suite (graphs built in
twins, one with its mask reset: a graph walk of the masked twin, a scan of the records, a scan of the unmasked
graph's entries computing U, BOSS's own dummy-tree traversal; 360 random patterns; the sampled f against an
exact count on 30 graphs and on the mini, 375 source and 12 sink dummies, the exact f inside the interval) and
`PatternPeptide.StopResidue` and `RandomStopPeptidesAgainstOracles` (a six-frame translation oracle, tables 1, 2,
11, 27, 28 and 31); the route's `PatternRoute.Unmasked*`, `PatternMaskUnmasked.*` (f against `stats --count-dummy`) and
`PatternRoutePeptide.TheStopIsAStopCodonOfTheTable` (tables 1, 2, 11, 27, 31); the identity's `Graphlet.IndexIdentity` and
`GraphletServer.DerivedDataIsWhatTheLoaderReadsAndNotTheIdentity`; the integration's `TestPatternMini`
(`test_unmasked_*`, `test_peptide_stops_against_the_six_frames`),
`TestPatternSynthetic.test_unmasked_in_every_graph_mode` and `TestTraverseDerivedDataMini`.

**Review of this round** (one finding, documentation only): the U bullet of §7.4 said that U is the masked count
plus the source dummies holding the pattern after their `$` run. For a pattern with a leading N run, which the
engine skips (§7.8), U also counts the source dummies whose `$` run ends inside that run (a `$` under an N, so
they do not hold the pattern). No relation was wrong, only that equality: corrected, and pinned by
`PatternUnmasked.LeadingNRunCountsDummiesWithTheirSentinelUnderTheN` (a named case and the decomposition of U
against the spelled entries, BASIC and CANONICAL).

**Owner decision #24 of 2026-10-08: tiny blocks exact** (after `226bc934`; the decision #22 it corrects allowed
no per-query check). On a graph without its mask, a pattern whose discovery and deferred scans completed with at
most `max_checked_entries` unchecked candidates (§7.4: the entries in U and not in `lower`, over its orientations
and offsets; `upper` − `lower` of its total on a BASIC or CANONICAL graph) has each of them tested at query time
(`BOSS::node_has_sentinel`, k − 1 steps each, charged to `max_steps`, phase `mask_scan`), and every count of the
pattern is then `exact` — the total, `suffix`, `by_offset`, `by_strand` / `by_orientation`, the anchors of a long
pattern — in every mode, before any release, so that a count and a retrieval of the same request agree. A block
of dummies only is `exact` 0. A larger block is not touched (no in-block sampling: the owner rejected per-query
checks of large blocks as too expensive) and answers as before, field for field. Contract version 1 stays:
additions only.
- **Server and CLI:** the flag `--pattern-max-checked-entries` (default 50; any integer in [0, 1,000], refused at
  start-up beyond it; 0 checks nothing and answers exactly as `226bc934`), §4.5. Not a request field, not echoed
  in `limits`.
- **Capabilities** (both routes, every single-graph server, masked or not): `caps.max_checked_entries` and one
  sentence of `caps_rule` ("max_checked_entries: unmasked, so few unchecked candidates are tested: exact.").
  For room under the 32 KiB ceiling of the capabilities document a service's MCP tool returns whole (§10.2 and
  the 172 bytes above), `caps_rule` also drops "(long_search lists the values served)", which the field
  `long_search` itself states: the mini's `/traverse/capabilities` (`-i`/`-a` only, as
  `integration_tests/test_pattern.py` serves it) grows from 32,634 to 32,699 bytes as served, under that
  test's guard of 32,704 (the ceiling less 64); the fixture servers' (compact JSON, as above) from 32,490 and
  32,596 to 32,555 and 32,661.
- **Engine** (`pattern_search.{hpp,cpp}`): `Request::max_checked_entries` (0 by default in the engine; the route
  passes the server's), `kDefaultMaxCheckedEntries`; the unchecked ranges are kept for the check while their
  candidates are at most the limit (at most that many ranges) and freed once above; while within it, the
  retention of `all_or_count` and of the extension compares the running lower bound (the check may still make the
  count `exact` within the threshold); `stop_at_threshold` still compares the running U, so it can stop a
  pattern the check would have made `exact`. `work.mask_scans` counts the ranges checked, `work.steps` the k − 1
  per candidate.

**Answers that change** (a graph without its mask only; with the check off they are `226bc934`'s byte for byte):
- A count with 1 to 50 unchecked candidates and no stop: `exact` (the masked graph's count) where it was `bounds`
  with an `estimate` and the note `estimate_sampled_dummy_fraction`, its `work.steps` 30 per candidate more at
  k = 31; in `all_or_count` it is then admitted on its `exact` count (no `threshold_upper_bound` where U was above
  `max_contexts`); in `partial` a cut list keeps the masked graph's `exact` counts instead of raised bounds. On the
  mini: blaNDM-1's primers `exact` 24 (22 tested), the 16S V4 primers `exact` 26 and 24, the island start
  `GATGCCGGTGAACAAC` `exact` 17 (30 tested, 15 source dummies dropped), its `suffix`-scope forward count `exact` 0;
  its first 14 bases (68 unchecked) and `GCCGAATTCGGC` (209) stay `bounds`.
- A count that a complete release made `exact` before: the same counts and list, the check's steps added.
- The check's steps come out of the request's `max_steps`: a request near its step cap can stop earlier, in the
  check (`mask_scan`, the counts `bounds` as before) or at a later pattern (`unknown` where it was `at_least` 0).
- **Masked graphs: nothing** but the capabilities' `caps.max_checked_entries` and `caps_rule`. Checked against
  `bin_226bc934` on a masked copy of the mini (`transform --mask-dummy`): 166 `/pattern` requests (every stored
  request of a mini server without a time budget, and variants: three modes, both scopes, three strand settings,
  thresholds, `stop_at_threshold`, `max_steps`, labels `all`, long patterns with and without paths, peptides)
  identical apart from `timing`; and the same panel on the unmasked mini with `--pattern-max-checked-entries 0`
  identical apart from `timing`.

**Fixtures** (§11; 90 in all): new `unmasked_checked` and `unmasked_checked_dummies` (server `unmasked`, the
default); a new server `unmasked_unchecked` (`--pattern-max-checked-entries 0`) now carries `unmasked_count`,
`unmasked_labels_all`, `unmasked_threshold_upper_bound` and `unmasked_partial`, their bodies unchanged (each shows
what a pattern above the limit answers; at the default their small blocks are checked); changed: the nine
capabilities bodies of the single-graph servers (`caps.max_checked_entries`, `caps_rule`), `README.md`,
`index.json`. Every other stored body is unchanged (`--check`). The validator knows the cap (`CAPS`) and refuses a
`bounds` total without a stop whose unchecked candidates (BASIC, CANONICAL; not in `partial`) are within the
server's `max_checked_entries` (a mutated `unmasked_checked`).

**Tests:** `PatternUnmasked.TinyBlocksCheckedExact` (named cases: ACG on ACGTTGCA, `bounds` [1, 3], 2 checked:
`exact` 1, 8 steps; the limit 1 as 0; step stops in the check; `all_or_count` at 1 released; `stop_at_threshold`
unchanged; a block of one dummy, and NC's four, `exact` 0), `TinyBlocksAgainstTheMaskedTwin` (360 random cases on
twins of every mode, an oracle of the unchecked candidates from the spelled entries: at the limit E every count the
masked twin's and k − 1 steps per candidate; at E − 1 and 0 the answer without the check, field for field;
`all_or_count`, `partial`, the extension, step stops at three points of the check),
`TinyBlocksOfLongPatternsAndPeptides`; the route's `PatternRoute.UnmaskedTinyBlocksAreExact` (three modes, both limits, every mode of
request) and `PatternMaskUnmasked.CheckedEntriesFlag` (the flag through the real Config and loader, refusals);
`integration_tests/test_pattern.py` `test_unmasked_against_the_masked` (the mini: the default server against
the masked one and the server with the check off: counts, lists and labels, a pattern's few unchecked candidates
each tested at the default limit) and the CLI with the flag.

**Stated limitations:** the limit counts candidates, not the pattern's U: on a wrapped PRIMARY graph a candidate
enters both orientations' counts, so `upper` − `lower` there can be up to twice the number checked; a cut
`partial` raises `lower`, so its `upper` − `lower` says nothing of the check; `stop_at_threshold` compares U before
the check; the capabilities document keeps 5 bytes under the integration test's guard on the mini.

**Review GPT-3 and the owner's decision P9 (round fix3, 2026-10-08, after `7b8354f2`).** Five findings of the third
outside review fixed in the engine (`pattern_search.{hpp,cpp}`, `dbg_succinct.{hpp,cpp}`) and the labelled retrieval
(`pattern_retrieval.{hpp,cpp}`), their counters stated by the route (`pattern.cpp`), and room made in the
capabilities (P9 of the owner's decisions of 2026-10-08 on increments 5b and 5s). **Contract version 1 stays:
additions only.** No field changes its meaning or its guarantee; what changes is where the work stops — on time and
stated, where the base binary ran past its work time into a late answer or a 503 —, what the memory account holds,
and when the note `low_complexity_pattern` is stated. The rule each loop keeps: the work is charged and the clock
read before it (at most every 4,096 units), memory admitted before a copy, a result kept rather than computed again,
and a stop stated in the answer rather than a 503 after the budget.
- **The extension reads the clock (finding 5).** It listed, spelled and extended the anchors with the clock read
  only every 4,096 steps, and spelling an anchor (k − 1 BOSS steps) is charged no step: the 40 bp half-N pattern of
  the staging benchmark (674 anchors on refseq33m, 1,129 extension steps) never reached a reading and, cold, ran
  past its work time into a 503; on the mini, with the base binary, `GNNNCNNNTNNNANNNGNNNCNNNTNNNANNNGNNNCNNN` (331
  anchors, 615 extension edges) under a 165 ms budget and a 1 ms reserve extended for 3.4 ms without a reading and
  answered 503. The clock is now read before every anchor's extension, every 64 anchors listed and before every
  64th node expanded (§7.6, §12.1): a late clock gives `stop {extension, time}`, `extension: "stopped"`, the paths
  `at_least` (`all_or_count`: `withheld: deadline`; `partial`: `cut: time`).
- **The low-complexity diagnostic (finding 2).** One sdust over the whole text took 786 ms on ATG × 10,000 without
  a clock reading. It now reads pieces of 128 + 63 bases, stops at the first piece flagged (about 6 ms there; 30 kb
  of random sequence about 0.2 ms), reads the clock before every piece but the first, and runs on a completed
  search only (§7.8): its flag is the whole text's (`PatternSearchFixes.LowComplexityNoteAsSdustOverTheWholePattern`,
  300 texts against sdust over the whole text).
- **Runs of sink edges (finding 3; graphs without the mask).** The release, the palindrome scans and the check of
  few unchecked candidates stepped over the sink edges (W = `$`) of a range one at a time;
  `DBGSuccinct::next_non_sink_edge` reads W at the first 16 and then jumps to the next symbol that is not `$` by
  rank and select: about 5.5 µs over 100,000 sinks as over 40 (a graph of 100,000 sinks: 2.5–2.8 ms against a work
  time of 0.1–1 ms before, 0.1–0.2 ms now). Same results.
- **Verifying long paths (finding 1).** The verification of §12.1: each (path, label) once, by a leapfrog join of
  its k-mers' shifted coordinate lists with consecutive chains as one run, each run placed record by record, the
  runs kept (charged) for the output, which no longer joins the chains again; every seek, run and record clocked;
  the label lists of the paths intersected in place from the shortest (each row's list sorted once). A homopolymer
  path of 1,500 bases with one occurrence per record (28,501 occurrences, `record_verified`) was a 503 after
  2.4–3.8 s; it is answered in 7–9 ms. `(AC)^750` as a path, 500 ms: a 503 after 1.75 s, now answered in 75 ms;
  `(AC)^5000`: `stop {label_discovery, time}` at 477 ms before, `stop {placement, time}` at 241 ms now (work time 250
  ms).
- **The occurrence cap (finding 4).** §14.4: in `partial` a context's or path's label holds, and the account and the
  estimate E count, only its first `max_occurrences_per_label` occurrences; the unions take every one (the counts
  stay exact); once complete, the lists are cut to their union's first and the rest is given back. Ten patterns of
  29,998 occurrences each under a cap of 1: before, the first stopped `{output, time}` (E counted 10 MB of text no
  list would show; its occurrences `at_least` 0) and the nine after it `{discovery, time}`; now all ten are answered,
  `exact`, about 2.2 MB in the account, in 52–64 ms.
- **The counters** (§8.4, §8.7), additive fields of an entry's `work` and `timing`: `extension_anchors`,
  `extension_branches` (every entry of a path search, beside `extension_edges`), `annotation_rows_distinct` (every
  entry with `labels: "all"`, 0 where nothing was read), `verification_steps`, `label_intersection_ms`,
  `verification_ms` (entries of a path search with `labels: "all"`). The request's totals are the sums of its
  entries' and are not stated apart. For a benchmark: the anchors' uncharged spelling is about `extension_anchors`
  × (k − 1) BOSS steps, the search's work `extension_edges` (charged), `extension_branches` says whether the paths
  fan out, and `timing.extension_ms` gives the time. `scripts/traversal/bench_pattern_retrieval.py` reruns the
  benchmark of this round, the retrieval's and the engine's repros (it builds its homopolymer, dinucleotide and
  100,000-sink indexes, runs `(ATG)^10,000`, `M^10,000` and the half-N count's budget sweep on the mini, and
  compares a base binary with a candidate: untimed bodies, stops and refused runs).
- **Capabilities (P9).** `caps_rule` and `protein_rule` are references to this document (§10.2), printable ASCII (a
  section sign would be written `§`), and the delivery rates of §7.6 are numbers, `delivery_mbps` (`build`,
  `compress`). Every machine-readable field is kept. Measured as served (the body's bytes, compact JSON), the base
  binary against this build:

| server (fixtures, §11; the integration test's mini servers) | `/traverse/capabilities` before → after | left under 32,768 | under 32,768 − 1,024 |
|---|---|---|---|
| `masked` | 32,609 → 30,966 | 1,802 | 778 |
| `unmasked` | 32,728 → 31,085 | 1,683 | 659 |
| `unmasked_unchecked` | 32,727 → 31,084 | 1,684 | 660 |
| `built_at_load` | 32,618 → 30,975 | 1,793 | 769 |
| `primary` | 32,597 → 30,954 | 1,814 | 790 |
| `masked_no_map` | 32,596 → 30,953 | 1,815 | 791 |
| `hash` | 32,549 → 30,906 | 1,862 | 838 |
| `multi` (`?graph=mini_refseq`, the reduced block) | 29,164 → 29,164 | 3,604 | 2,580 |
| mini, `-i`/`-a` only: masked, unmasked, unchecked | 32,580, 32,699, 32,698 → 30,937, 31,056, 31,055 | 1,712 at least | 688 at least |

  Every single-graph document is 1,643 bytes shorter; `/capabilities` likewise (the largest, `unmasked`, 25,790
  bytes). The largest probe document keeps 659 bytes under the budget a test now holds every document to,
  32,768 − 1,024 (the room the next increment's additions need): the predicate block P9 plans for increment 5b
  (about 356 bytes) fits. This replaces the 5 bytes left under the old guard of 64 (above).

**Answers that change** (against the base binary `8f49cc88`, whose code `7b8354f2` serves, on a masked copy of the
mini and on the mini as built):
- **The counters** above, in every entry they belong to, and in the capabilities `delivery_mbps` and the two rules
  rewritten. Nothing else of a capabilities document.
- **`low_complexity_pattern` is no longer stated beside a stop**, of any phase and reason: threshold stops
  (`max_contexts`, `max_anchors`, `max_paths` with `stop_at_threshold`) and an earlier pattern's request-wide stop
  included. The base binary stated it whenever sdust flagged an exact pattern. A pattern of more than 191 bases
  whose diagnostic the work time cut loses the note and is `time_limited` with `stop` `null` and complete counts
  (§7.8, §7.9). A withheld `count_above_threshold` (no stop) keeps it.
- **Time stops on time** in the extension, in the verification and the label lists of paths, and in the
  occurrences made for the output (§7.6, §12.1, §14.4): stated stops (`{extension, time}`, `{placement, time}`,
  `{output, time}`) where the base binary ran past its work time into an answer after the deadline or a 503. The
  readings are denser, so with the same budget a stop can come earlier; `deadline` 503s are rarer.
- **`partial` with the occurrence cap cutting** (a label's union holds more than `max_occurrences_per_label`):
  `work.memory_bytes` is lower — only the listed prefix is charged (on the mini: the short panel of the identity
  check with cap 1, 225,828 → 140,940 bytes for NDM-F; `labels_all_partial` 17,721 → 17,448;
  `labels_all_partial_exact_cut` 64,361 → 61,085; NDM's 40-mer paths with cap 1, 44,518 → 37,444) — and so is the
  estimate E: stops the inflated estimate caused no longer fire. The listed occurrences, the counts and
  `occurrences_cut` are unchanged where nothing stops.
- **Paths under memory pressure:** the runs a verification keeps are charged before the output, so the output's
  `max_memory` stop can come one or more paths earlier; when the account cannot hold a path's runs, it and the later
  paths are not output (`stop {output, max_memory}`, `labels_status: "output_budget"`), their labels verified all
  the same.
- **A quirk fixed:** when the paths' label lists stopped (time or memory) at a path after the first, the paths
  before it said `labels_status: "complete"` with `labels: []` and `labels_total` above 0; they say
  `output_budget` with `labels: null`.
- **After an output stop** the unions and occurrence counts can include more of the request than the base binary's:
  they stay `at_least` (true), possibly larger.
- **Graphs without the mask:** the same results, faster over runs of sink edges.

Checked with the identity panel (`2026-10-07/pattern/tiny-identity/panel.py` and this round's cases: the
occurrence cap at 0, 1 and 2 with `max_labels` 2, labelled paths in both modes with the cap and with
`require_support`, a low-complexity pattern after an earlier pattern's `max_steps` stop, and `GCGCGCGCGCGC` (60
contexts on the mini, flagged) under `stop_at_threshold` in both modes, withheld, and complete; 184 comparisons a
graph, both capabilities routes included, on the masked copy and on the mini as built): 119 identical apart from
`timing`, 57 different by the additions only, 8 different beyond them, each one of the changes above — `memory_bytes`
with the cap cutting (5) and the note left out beside a stop (3). The same on both graphs (without the mask the
threshold stop keeps `threshold_upper_bound`).

**Fixtures** (§11; 90, none new or removed). 38 bodies changed: the nine capabilities bodies of the single-graph
servers (`delivery_mbps`, the two rules); the 13 `labels_all*` with labels read and `unmasked_labels_all`
(`annotation_rows_distinct`), `peptide` likewise; the 11 `paths*`, `unmasked_paths`, `peptide_paths` and
`peptide_no_stop_codon` (the extension's counters, and with labels the retrieval's); and `memory_bytes` of
`labels_all_partial` and `labels_all_partial_exact_cut` (the cap, above). Every other stored body is unchanged
(`--check`). No fixture shows `extension_branches` above 0 (no path of the mini branches); the tests below do. The
validator knows the new fields (`SCHEMA`: `work`, `timing`, `capabilities`, the table `delivery_mbps`) and refuses
what version 1 never answers (`test_round_fix3_rules_refuse_what_v1_never_answers`: a counter where its entry has none
or missing where it has one, a completed extension that did not begin at every anchor, an extension that did not run
with work of its own, more branchings than half the candidates, more distinct rows than rows read, verification
steps without coordinates, the note beside a stop, `time_limited` without a stop on a pattern its diagnostic reads
whole, a prose field that is no ASCII reference, a delivery rate that is not a positive number); every fixture
server's capabilities documents are held to 32,768 − 1,024 bytes as the server writes them (floats counted at 24
characters; `test_capabilities_documents_keep_a_kibibyte`, which also reads the ceiling from `mcp_tools.py`). Its
regression body `pattern_validator/by_label_null_partial` was answered again by this build (identical apart from
`timing` and `annotation_rows_distinct`). The code lists name `pattern_predicate.cpp` (increment 5b's predicate
language, built and unit-tested, not called by the route: `predicate` is still refused by name) as a source whose
new code `predicate_too_large` no answer of this build carries (`NOT_SERVED_SOURCES`).

**Tests with independent oracles:** the engine's `PatternUnmasked.LongSinkRunSkippedByRankAndSelect` (against a read
of W; the cost over 100,000 sinks under 20 times the cost over 40), `PatternSearchFixes.LowComplexity*` (sdust over
the whole text, the note after a stop, the clock) and `ExtensionReadsTheClockBeforeEveryAnchor`, and in every
`PatternSearch.Extension*` case `extension_anchors` and `extension_branches` against a path oracle; the retrieval's
`PatternPaths.RepeatsAgainstTheOracles` (graph-walk and record-scan oracles, plain and `record_verified`),
`TheCapListsTheFirstOccurrencesOfEachUnion`, `AHomopolymersChainsAreOneRun`, `TheVerificationReadsTheClock`,
`MemorySweepOverRepeats`, `TimeSweepOverRepeats`, `TheRetrievalCounters` and
`PatternRetrieval.TheCapListsTheFirstOfEachUnion`, `TheCapEstimatesOnlyWhatIsListed`, `TheOccurrencesReadTheClock`;
the route's `PatternRoute.TheCountersOfTheExtensionAndOfTheLabels` (a brute force over the records' k-mers: the
anchors, the branchings, the paths and the distinct rows of contexts and paths, in every mode, with and without
labels and paths, and the extension not admitted) and `PatternRoute.Capabilities` (the references, ASCII, every cap
named, the rates as configured); the integration's `assertPaths` (`extension_anchors` and `extension_branches`
against the graph-walk oracle `Records.branchings`, the distinct rows against the paths' k-mers,
`verification_steps` at least n × the labels carrying the paths), `assertLabelled` (the distinct rows),
`test_capabilities` and the validator's `test_capabilities_documents_keep_a_kibibyte` (every stored capabilities
document, which `pattern_fixtures.py --check` keeps equal to what the servers send, at most 32,768 − 1,024 bytes). Each engine and retrieval test fails on a mutant without its fix.

**Stated limitations:**
- Discovery still reads the clock every 4,096 steps: cold on refseq33m about 12 µs a step (4.5 million steps in
  about 60 s), so up to about 50 ms between readings.
- The reading before each anchor also asks whether the client has left (two system calls), at most `max_anchors`
  times a pattern; small next to spelling the anchor.
- The note is left out beside threshold stops too, which are exactly the answers whose counts are large; since the
  diagnostic is bounded and clocked, leaving it out after budget stops only (`max_steps`, `time`) is open for the
  owner.
- The join of a path's verification assumes a row's coordinates of one label distinct (each one k-mer position);
  the occurrences of contexts tolerate repeated coordinates (counted once).
- The 100,000-sink test compares costs (a ratio, not machine speed): a timing test still.
- The item-5 cause was confirmed on the mini, not re-measured cold on refseq33m: rerunning the staging benchmark's
  `GNCNGNTNCNGNANANANTNANGNTNTNCNCNANGNGNTN` cold should now give a stated stop.

## 19. Increment 5b: annotation predicates (`predicate`), patterns of L ≤ k

### 19.1 What is served

A request may name **one predicate** (owner decision P21): a logical condition on the annotation columns of each
graph context of a pattern of at most k bases (§7.1). The predicate **selects** contexts; the projection
(`output.labels`) decides which labels a selected context is returned with. Served by this build
(`src/cli/pattern.cpp`: the request, the per-pattern pipeline, the answer; `src/cli/pattern_predicate.cpp`: the
language, its binding and evaluation; `src/cli/pattern_selection.cpp`: the selection pass and the projection
`"predicate_only"`), as additions to version 1 (§1): a request that names no predicate is answered as before,
byte for byte apart from `timing` (§19.15). The capabilities announce it (§19.12). Design:
`DESIGN-pattern-search.md` §5.6, as corrected by the owner's decisions of 2026-10-08 (P1–P22, P24): the strands of
a BASIC graph (§19.5, P11), the default projection (P2), the relation of an interrupted selection (§19.7, P13),
the absence claim (§19.11, P15), the honest units of a read (§19.9, P17).

Not in this increment: the selection of patterns longer than k (a predicate selects among **supported paths**,
`long_search: "supported_paths"`, a later increment, P24). A pattern longer than k in a predicate request keeps the
answer of §7.7 (`long_search: "anchors"`), its selection `not_started`; with `long_search: "paths"` a predicate is
400 `invalid_request`. Labels are column names only: no taxonomy (owner decision P6; a cohort is an explicit list,
its expansion a client's business).

### 19.2 Request

The fields are §4.1's rows `predicate`, `max_predicate_contexts`, `max_predicate_work`, `predicate_strands`;
`output.labels` gains `"predicate_only"` (§4.3). With a predicate `output.labels` may be:

- `"none"` (the default, frozen, owner decision P2): the selected contexts without labels;
- `"predicate_only"`: each selected context with the predicate's labels on its **own** row, placed where the
  index can (`output.occurrences`, default `true`, as with `"all"`); the service sends it explicitly;
- `"all"` (P8): the ordinary retrieval of §14 on the selected contexts (their rows read again).

In mode `count` no projection is built (note `projection_not_read` when one was named, §19.10). With a predicate
`max_memory_mb` and `allow_unbudgeted_annotation` act in every mode (the selection reads rows under the account);
`max_labels_per_anchor`, `max_annotation_work`, `max_labels`, `max_occurrences_per_label` and
`require_support` act only on a projection.

### 19.3 The predicate

<!-- schema: predicate -->
| operator | form | meaning on a context's label set S |
|---|---|---|
| `any` | `{"any": [n1, n2, …]}` | at least one listed name is in S |
| `all` | `{"all": [n1, n2, …]}` | every listed name is in S |
| `none` | `{"none": [n1, n2, …]}` | no listed name is in S |
| `at_least` | `{"at_least": {"n": m, "labels": [n1, …]}}` | at least m of the listed names are in S |
| `and` | `{"and": [p1, p2, …]}` | every operand holds |
| `or` | `{"or": [p1, p2, …]}` | at least one operand holds |
| `not` | `{"not": p}` | p does not hold |

(P1.) Each rule is refused with 400 `invalid_request` naming the path of the first fault
(`request.predicate.and[1].none[0]: …`), unless said otherwise:

- a predicate is a JSON object with **exactly one** member, one of the seven operators;
- a name is a **non-empty JSON string**, an annotation column label as stored (case-sensitive, not trimmed; on
  refseq33m a taxid written as a string, `"562"`); a number is refused with the fix (P7: `write a taxid as
  "562"`, fixture `predicate_invalid`); a record header is not a predicate term in version 1 (P20: unknown,
  §19.4);
- a list (`any`, `all`, `none`, `at_least.labels`) has at least one name and no name twice;
- `at_least` is an object of exactly `n` and `labels`, `n` an integer with 1 ≤ `n` ≤ the list's length;
- `and` and `or` take a non-empty list of predicates, `not` one predicate; `not` may stand anywhere, the top
  level included (P5); nesting (`and`, `or`, `not`) at most 64 deep;
- **size**: the names of all lists together (a name in two lists counted twice) at most
  `caps.max_predicate_labels` (`--pattern-max-predicate-labels`, 10,000 by default, at most 1,000,000; not a
  request field, P3): above it 400 **`predicate_too_large`**, the message naming the count and the cap. Past the
  cap a name's form is still checked, its repetition in one list no longer: a body above the cap is refused for
  its size (the parse keeps at most the cap's names, whatever the body's size).

The form is checked in step 8 of §5, after `genetic_code`. Examples on the mini index (its columns are the
taxids 1296536, 158836, 287, 470, 546, 562, 573, 615, 72407): `{"and": [{"any": ["562"]}, {"none": ["287"]}]}`,
`{"at_least": {"n": 2, "labels": ["562", "573", "615", "546"]}}`, `{"or": [{"all": ["546", "615"]}, {"not":
{"any": ["287"]}}]}`.

### 19.4 Binding: unknown names, the normal form, constants

The names are resolved once per request, before the first pattern, against the index's column labels (one hash
lookup each; the work time read every 4,096 names). A name that is not a column is **unknown**: absent from every
context. Unknown names are folded away before anything is read:

| form | with unknown names |
|---|---|
| `any(L)` | the known names; none known: `false` |
| `all(L)` | an unknown name: `false` |
| `none(L)` | the known names; none known: `true` |
| `at_least(m, L)` | the known names; fewer than m known: `false` |
| `and` | `false` if an operand is `false`; `true` operands dropped; none left: `true`; one left: that operand |
| `or` | `true` if an operand is `true`; `false` operands dropped; none left: `false`; one left: that operand |
| `not` | of a constant: the other constant |

The result is the **normal form**, echoed in the answer (`predicate.normal_form`, P19) and evaluated. A normal form
that is the constant `false` or `true` is evaluated **without reading any annotation** (P18; `selection.pass:
"constant"`, note `predicate_constant`): `false` selects nothing, `true` every context. A predicate is **vacuous**
when its normal form holds on a context carrying none of its labels (`none(A)`, `not(any(A))`): it passes contexts
the index labels with other columns only, and the answer says so (`predicate.vacuous`). The bound predicate (its
labels, nodes and unknown names; model below) is charged to the request's memory account; a binding the account
cannot hold, or that the work time stops, leaves no bound predicate: every pattern's selection is then
`not_started` with `stop {selection, max_memory | time}`, and the predicate block states its `names` only (the
rest `null`).

### 19.5 What is selected

A context is selected when the normal form holds on its label set S: the predicate's labels present on the
annotation row of the context's k-mer x (§7.10) and, with `predicate_strands: "either"` (the default, owner
decision P11) on a BASIC graph, those on the row of rc(x). **Support is strand-consistent** (the owner's answer to
P11): a label supports a context in one orientation as a whole — it annotates x, or it annotates rc(x); for L ≤ k
this is a single k-mer either way, so nothing can be mixed. rc(x) has no row when no record holds it as
deposited (it then adds nothing). `"context"` reads x's row only (for stranded indexes). On CANONICAL and PRIMARY
graphs x and rc(x) share one row: S is that row's whatever is asked, and the answer states `"either"`
(`predicate.strands`).

Why `"either"` is the default: on a BASIC graph (refseq33m) a k-mer's row lists the columns whose records hold it
as deposited; a record holding the motif on its other strand annotates rc(x). On the mini the blaNDM-1 forward
primer's 12 `+` contexts carry 9 columns and its 12 `-` contexts 7 (546 and 615 hold the gene on `+` only):
`none(["546"])` selects the 12 `-` contexts with `"context"` although 546's records carry the primer, and none
with `"either"` (fixtures `predicate_context`, `predicate_either`).

**Scope `shard_context`** (`predicate.scope`): the predicate is asked of each context on its own, on this index —
"is this sequence context supported by these labels", never "does the motif occur anywhere in A and nowhere in C"
(a motif-level predicate is a later increment, design §12). With k = 5 and the pattern `AC`, a sample A holding
`TACGG` and a sample C holding `CACCC` give two contexts, and `any(A) and none(C)` selects A's although C carries
the motif in another flank.

### 19.6 The pipeline, per pattern

1. **Discovery** of the raw contexts as without a predicate, against the threshold `max_predicate_contexts`
   instead of `max_contexts`: `stop_at_threshold` stops it once its running lower bound (without the mask: its
   running upper bound) passes `max_predicate_contexts`, `stop {discovery, max_predicate_contexts}`.
2. **Compute admission**: the raw count `exact` (without the mask: its upper bound) at most
   `max_predicate_contexts`. Above it nothing is read: `all_or_count` withholds `predicate_above_threshold`,
   `count` states `selection.pass: "not_admitted"`, `partial` tests the first `max_predicate_contexts` in
   answer order (cut `max_predicate_contexts`). Without the mask an admitted release enumerates every candidate
   and drops the source dummies (§7.4): the raw count becomes `exact` (fixture `predicate_unmasked`). Each raw
   context kept for the pass is charged 64 bytes before it is kept (§19.9); the first that does not fit ends the
   pass's set (`all_or_count` reads nothing, `predicate_budget`; the other modes test what was admitted, `cut:
   max_memory`).
3. **The selection pass** (`stop` phase `selection`): (a) the rows: each tested context's row and, with
   `"either"` on a BASIC graph, its reverse complement's — the k-mer spelled, reverse-complemented and looked
   up, one lookup per distinct row, kept (the mirror of a mirror is known: no lookup); (b) the reads: one row per
   read in the order of first appearance (a context's row before its mirror's), each **restricted to the
   predicate's known labels** (rows access, budget-aware; an unbudgeted annotation with direct access single
   cells for at most 16 labels, `selection.access: "columns"`), the time and `max_predicate_work` checked before
   each; a row the account cannot hold is stated (`rows_refused`, phase `selection`) and its contexts stay
   untested (`all_or_count` stops reading there) — but an unbudgeted read (`allow_unbudgeted_annotation`), whose
   size is known only after it, that the account cannot hold stops the reads with `stop {selection,
   max_memory}` and is not stated in `rows_refused`, as §14.4's unbudgeted reads stop; (c) the decisions, in answer order as the rows are read: a
   context is **tested** once every row it needs was read, **selected** when the normal form holds on S.
   `stop_at_threshold` ends the pass once more than `max_contexts` are selected (`stop {selection,
   max_contexts}`); after any stop but time, every context whose rows were read is still decided.
4. **Selection admission**: the selected count `exact` and at most `max_contexts`. Above it `all_or_count`
   withholds `selected_above_threshold` (P14), the counts kept; `partial` returns the first `max_contexts`
   selected in answer order (cut `max_contexts`).
5. **The results** of the selected contexts, in answer order, each charged 512 + 2k bytes before its object is
   built, the clock read every 64 (`stop` phase `output`: `time`; without a projection `max_memory`, which the
   projections state as §14.4 does), then **the projection**: none; `"predicate_only"` (the pass's own rows of
   the selected contexts, their predicate labels placed by a second read of those labels' coordinates under
   `max_annotation_work`); `"all"` (the retrieval of §14).
6. **Publication**: `all_or_count` publishes every selected context or none.

| mode | with a predicate |
|---|---|
| `count` | steps 1–3 (the pass runs: a `count` with a predicate reads annotation); no results, no `withheld`, `returned`, `cut`; `retrieval_complete` false; the pass in `selection.pass` |
| `all_or_count` | steps 1–6; every selected context or none |
| `partial` | steps 1–6 on the first `max_predicate_contexts` raw contexts; the pass tests every one of them it can (P4); the first `max_contexts` selected are returned |

A normal form `false` runs discovery as mode `count` does (nothing retained, no read) and answers `selected`
`exact` 0 and, in a retrieval mode, `results: []` with `retrieval_complete: true` — also after a discovery stop: no
context can pass (P18). `true` is the unfiltered request: the engine's own thresholds (`max_contexts`) and no
read for the selection; `selected` is the raw count with its relation (`"predicate_only"` then lists no label:
the predicate names no column of the index).

### 19.7 Counts and their relations

`counts.tested` and `counts.selected` (§8.6), unit `graph_contexts` (`paths` for a pattern longer than k). With
R the raw count (`counts.contexts`), T the decisions made and S the selected among them:

| `selection.pass` | when | `tested` | `selected` |
|---|---|---|---|
| `completed` | every raw context was tested (R `exact`, T = R) | `exact` T | `exact` S |
| `stopped` | the pass stopped (work, memory, time, `stop_at_threshold`), a row was refused, the descriptors' admission cut its set, or `partial`'s release held fewer than R | `exact` T | `bounds` [S, S + R_upper − T] where R has an upper bound (`exact` R, or `bounds` up to R_upper; P13); `at_least` S where R is `at_least` |
| `not_admitted` | R above `max_predicate_contexts` (`all_or_count`, `count`) | `unknown` | `unknown` |
| `not_started` | a stop before the pass (the engine's withheld release; in `partial` an engine stop, `max_steps` or `time`, before it released a context, so that every mode answers that stop alike; a sticky `max_predicate_work` stop of an earlier pattern; a binding that stopped), or a pattern longer than k | `unknown` | `unknown`; `exact` 0 when R is `exact` 0 |
| `constant` | the normal form is `false` or `true` | R, its relation | `false`: `exact` 0; `true`: R, its relation |

- `selected` ≤ `tested` ≤ the raw count's upper bound, whatever the relations.
- The raw counts keep their meaning: what discovery found, unfiltered. A predicate never changes them, except that
  an admitted `all_or_count` release on a graph without its mask makes them `exact`, as without a predicate.
- `counts.labels` and `counts.occurrences` are over the returned contexts and the projection's labels: `unknown`
  with `"none"`; with `"predicate_only"` the predicate's labels found on the returned contexts' own rows (`exact`
  when the list and its placement are complete); with `"all"` as §14.6.
- A `selected` count carries no `estimate`, also without the mask.

### 19.8 `withheld`, `cut`, `stop`

`withheld.reason` (`all_or_count`; §7.5's table): `predicate_above_threshold` (the compute admission, nothing
read), `selected_above_threshold` (more selected than `max_contexts`, the counts kept), `predicate_budget` (the
pass stopped at `max_predicate_work`, also an earlier pattern's, sticky; or the account could not hold the
bound predicate, the descriptors, a row, a statement); and, unchanged in name: `deadline` (a time stop of the
pass, `stop {selection, time}`, or of the results, `stop {output, time}`), `threshold_crossed` (`stop_at_threshold`
on the raw count, reason `max_predicate_contexts`, or on the selected count, reason `max_contexts`),
`discovery_budget`, `output_budget` (the results' objects), and §14.6's reasons for the projection's reads.
`count_above_threshold` is not used on L ≤ k with a predicate (its raw threshold is `predicate_above_threshold`).

`cut.reason` (`partial`), the first that applies in this order: the engine's stop (`max_steps`, `time`), the
pass's stop (`max_predicate_work`, `max_memory`, `time`, `max_contexts` of `stop_at_threshold`; `max_memory` also
for refused rows and descriptors the account could not hold), the raw release's cut (`max_predicate_contexts`),
the selected list's cut (`max_contexts`); a cut of the results themselves (`time`, `max_memory`) replaces them
(the list's length is the results', as §7.5's memory cut). An incompleteness of the projection alone (a refused
placement row, a time stop of its reads) sets `retrieval_complete: false` without a cut, as §7.5 says.

`stop`: the phase **`selection`** (the pass: lookups, reads, decisions) with the reasons `max_predicate_work`,
`max_memory`, `time`, `max_contexts`; the reason **`max_predicate_contexts`** for the phase `discovery`. First stop
wins (§7.6): the engine's, then the pass's, then the results' (`output`), then the projection's. A time stop of the
pass is recorded in the request's budget: the later patterns answer as after any time stop. A `max_predicate_work`
stop is sticky for the selection of the later patterns (their discovery runs as mode `count` does, their pass
`not_started`, `withheld: predicate_budget` / `cut: max_predicate_work`), as `max_annotation_work` is for §14.4's
reads.

### 19.9 Budgets: work, memory, the deadline

**Work** (`max_predicate_work`, one budget per request, separate from `max_annotation_work`, P3), in the
oracle's units:

| item | units | when charged |
|---|---|---|
| a row read by the pass | 8 + E + D: E the decoded row's entries — **all its labels, not only the predicate's** (P17: a row-diff decode builds the whole row) — and D its row-diff dependency units (8 per dependency row, 1 per entry it stores) | gate `units < max_predicate_work` before the read; the units when it returns. A refused read: what its decode reached (counting the restricted hits, at least 8: the whole row's size is not returned for a refused key); an interrupted one at least 8; an unbudgeted read 8 + its hits (the row's size is unknown) |
| a reverse-complement lookup (`"either"`, BASIC) | k | gate before it; one per distinct row, none for a row found as another's mirror |
| deciding a context | 1, and 1 per leaf each distinct present label is listed in | after its rows were read; never refused |

The pass checks `units < max_predicate_work` before every read and every lookup: it passes the budget by at most
one read's units, as §14.4's reads do. All lookups come before the first read (step 3a), so an early exit by
`stop_at_threshold` saves reads, not lookups.

**Memory** (the request's account, `max_memory_mb`, §14.4; with a predicate it exists in every mode). The model,
beside §14.4's items:

| item | bytes (model) | when charged |
|---|---|---|
| the bound predicate | per label of the normal form 192 + 2 × its length; per name a leaf lists 8; per node 64; per unknown name 64 + its length | once, before the first pattern (its echo's text is covered by it) |
| a raw context kept for the pass | 64 | as the engine releases it, before it is kept |
| a distinct row of the pass | 96 (its entry and the index of its key) | before it is added (step 3a) |
| a row's hits | as the `DecodeBudget` models them (§14.4), restricted to the predicate's labels | by the read; held until its contexts are decided, then freed (but a listed context's own row with `"predicate_only"`) |
| the statement of a refused row | 384 + k | reserved before each read |
| a listed context's `selection_labels` and `selection_strands` | 32 + the names' lengths + 4 per label (the ids kept until the results name them); 32 + 24 per label (a byte kept per label, its string in the answer) | when it is decided |
| a label of the listed contexts' label order | 48 | likewise, once per distinct label |
| `"predicate_only"`'s dictionary | 192 + 2 × the length, per distinct predicate label on the listed rows | likewise |
| a selected context's result object | 512 + 2k (§14.4) | before it is built |

**The deadline** (§7.6): read while the names are bound (every 4,096), every 4,096 contexts and before every
lookup batch of 64, before every row read (and between its paced pieces), every 64 decisions, before the label
order of the listed contexts (stop `{output, time}`, the list then dropped; its sorts, n log n comparisons each,
are then counted as light work, so that the next light work reads the clock once they pass 4,096) and every 64
result objects built. A time stop: `stop {selection, time}` (or `{output, time}`), `determinism: "time_limited"`.
The loops that only free what the pass held (its rows, at most 2 × `max_predicate_contexts`) are not clocked: a
stop would not skip them.

The engine's own retention for the raw release (24-byte descriptors, up to `max_predicate_contexts`) is outside
the account, as the label-free path's is (§7.6). The result objects of `output.labels: "none"` are charged while
their pattern is answered and released when the next pattern's release begins.

### 19.10 The answer

**Top level.** `output` (§8.1): `{"labels": "predicate_only", "occurrences": <bool>}` beside the other values.
`predicate`, present exactly when the request named one:

<!-- schema: predicate_block -->
| field | type | meaning |
|---|---|---|
| `normal_form` | object \| boolean \| null | the predicate after folding (§19.4), in the request's syntax; `false` / `true` for a constant; `null` when the binding stopped |
| `names` | integer | the distinct names of the request's lists |
| `known` | integer \| null | of them, the index's columns (also the ones the folding dropped) |
| `unknown_labels` | list of strings \| null | the others, in the order of their first appearance in the request: a typo shows here, never as an absence |
| `vacuous` | boolean \| null | the normal form holds on a context carrying none of its labels (§19.4) |
| `scope` | `"shard_context"` | §19.5: per context, per index; not motif-level |
| `strands` | `"either"` \| `"context"` | what was evaluated (§19.5): `"either"` on CANONICAL and PRIMARY graphs whatever was asked (`limits.predicate_strands` echoes the request) |

**Entry.** Every answered entry of a predicate request gains `counts.tested`, `counts.selected` (§19.7),
`absence_filter: "predicate"` (§19.11), `rows_refused` (in every mode; the pass's first), `work.predicate_rows`,
`work.predicate_units`, `work.predicate_lookups`, `work.memory_bytes` (§8.7), `timing.selection_ms` (§8.4), and:

<!-- schema: selection -->
| field | type | meaning |
|---|---|---|
| `pass` | `"completed"` \| `"stopped"` \| `"not_admitted"` \| `"not_started"` \| `"constant"` | §19.7 |
| `support` | `"kmer"` \| `"label_intersection"` \| `"record_verified"` | what a label's presence means: on the context's k-mer (`kmer`, L ≤ k); for a pattern longer than k the index's best support of a walk (the level of the supported-path selection, a later increment) |
| `access` | `"rows"` \| `"columns"` | how the predicate's labels are read: rows decoded and restricted (budget-aware, or unbudgeted rows), or single cells (an unbudgeted annotation with direct access and at most 16 known labels, `allow_unbudgeted_annotation`) |

With `"none"` the entry has no `placement`, `annotation`, `by_label`, … (as §8.5); with `"predicate_only"` and
`"all"` it has the fields of §14. **Results** are the selected contexts, in answer order (§7.9); under a projection
that reads labels each carries `selection_labels` (§8.8, P22): the predicate's labels in the set it was evaluated
on, in label order — why it was selected (with `"either"` possibly a label of its reverse complement's row only,
which its own `labels` then lack) — and beside them `selection_strands` (§8.8; the owner's answer to P11: per
selected result and label the orientation that supported it): `"context"`, `"reverse_complement"` or `"both"` on a
BASIC graph (every label `"context"` under `predicate.strands: "context"`, which reads x's row only), `"either"` on
CANONICAL and PRIMARY graphs. Fixture `predicate_selection_strands`: `any(546)` on the blaNDM-1 forward primer
selects the 12 `+` contexts by their own row (`"context"`) and the 12 `-` ones by their reverse complement's
(`"reverse_complement"`, their own `labels` empty). With `"predicate_only"` its `labels` are the predicate's labels on its own row
(`labels_total` their number), placed as §14.3 places them; `annotation_rows` counts the placement's reads only (the
rows were read by the pass, `predicate_rows`), `annotation_rows_distinct` the distinct rows of the returned
contexts.

**`limits`** (§8.3): `max_predicate_contexts`, `max_predicate_work` (effective), `max_predicate_labels` (the
server's), `predicate_strands` (as requested), and the annotation limits of §8.3 in every mode; `clamped` may name
`max_predicate_contexts` and `max_predicate_work` (after `max_paths`).

**Notes** (§8.10, after the others in this order): `annotation_unbudgeted` (a selection that read an unbudgeted
annotation, when no projection said it), `predicate_constant`, `projection_not_read`. `annotation_not_read` is never
set on a predicate answer.

### 19.11 What a filtered answer licenses (`absence_filter: "predicate"`)

- `retrieval_complete: true` with a predicate: `selected` is `exact` and every selected context is in `results`
  with its projected labels. So a raw context of the pattern in its scope and strands that is **not** in
  `results` does not satisfy the predicate — on this index, per context (`scope: "shard_context"`), at
  `selection.support` and `predicate.strands`.
- `selected` `exact` 0: no context of the pattern in its scope and strands satisfies the predicate there.
- With `"predicate_only"` and `retrieval_complete: true`: a predicate label absent from a result's `labels` is absent
  from that context's own row; with `"either"` its `selection_labels` may still hold it, with `selection_strands`
  `"reverse_complement"` (the label carries the context's reverse complement as deposited, not the context).
- **Not licensed**: motif-level statements ("A carries the motif, C does not", §19.5); anything about the labels
  a projection did not read (`"none"`, or labels outside the predicate with `"predicate_only"`); anything about
  contexts beyond `tested` when `selected` is `bounds` or `at_least`; anything about unknown names (they are no
  columns of this index); with `"context"` on BASIC, anything about the other strand's records.
- A `count` answer with a predicate licenses its counts (`selected` `exact` 0 is the absence above), no list.
- `absence_filter: "predicate"` says that the entry's absence claims are this section's, not §9's unfiltered
  ones; the raw counts keep §9's.

### 19.12 Capabilities

Added to the `pattern` block (§10.2, both routes): `projections` lists `"predicate_only"`,
`projections_later_increment` is `[]`, `caps` gains `max_predicate_contexts` (100,000), `max_predicate_work`
(10⁸) and `max_predicate_labels` (10,000), `caps_rule` names them, and `predicate` (`capabilities_predicate`,
§10.2) states the operators, the strands and the access. A client offers predicates where `projections` lists
`"predicate_only"`, sends at most `caps.max_predicate_labels` names and `predicate_strands` from
`predicate.strands`. Measured on the mini: the largest document, the unmasked server's `/traverse/capabilities`,
grew by 285 bytes, from 31,085 to 31,370 as the server writes it, 374 under the budget of 32,768 − 1,024;
the validator's measure from above (floats counted at 24 characters) is 31,504,
240 under it (`test_capabilities_documents_keep_a_kibibyte`); the masked server's is 31,251 (31,368 from above).

### 19.13 Worked examples (the mini index, this build)

Measured with this build's server on the masked copy of the mini (the fixtures of §11; `timing` varies):

| fixture | request | answer |
|---|---|---|
| `predicate_filter` | `GCGGCGGCGGCG`, `all_or_count`, `max_contexts` 100, `and(any 562, none 287)`, `predicate_only`, occurrences false | 1,828 contexts tested in 1,777 rows (1,757 lookups, 744,365 units, ~20 ms), 60 selected and returned, complete (without the predicate: `withheld: count_above_threshold`); each with `selection_labels: ["562"]` and its own row's 562 |
| `predicate_either` / `predicate_context` | NDM-F, `count`, `none(546)` | 24 tested in 24 rows; 0 selected with `"either"` (12 lookups, 12,388 units), 12 (the `-` contexts) with `"context"` (12,004 units) |
| `predicate_budget` | GCG12, `"context"`, `max_predicate_work` 1, then NDM-F | the first row read (485 units: the budget is checked before a read), `tested` 1, `selected` `bounds` [0, 1827], `stop {selection, max_predicate_work}`, `predicate_budget`; NDM-F's discovery runs, its selection `not_started` (sticky) |
| `predicate_partial_admission` | GCG12, `partial`, `max_predicate_contexts` 100 | the first 100 raw contexts tested (103 rows), 6 selected, `bounds` [6, 1734], `cut: max_predicate_contexts` |
| `predicate_stop_at_threshold` | GCG12, `partial`, `max_contexts` 10, `stop_at_threshold` | 272 contexts decided in 269 of the 1,777 rows, 17 selected (`bounds` [17, 1573]), the first 10 returned, `stop {selection, max_contexts}` |
| `predicate_unknown_constant` | NDM-F, `any(5622)` | `unknown_labels: ["5622"]`, `normal_form: false`, `pass: "constant"`, `tested` `exact` 24, `selected` `exact` 0, no row, `predicate_constant`, complete |
| `predicate_only_record` | NDM-F, `any(562)`, `predicate_only` | 24 selected, each with 562 placed: 13 placed occurrences (`labels_all`'s 562 has 13); `selection_strands: ["both"]` (562's records hold the primer on both strands) |
| `predicate_selection_strands` | NDM-F, `any(546)`, `predicate_only`, occurrences false | 24 selected: the 12 `+` contexts by their own row (`selection_strands: ["context"]`, `labels` [546]), the 12 `-` ones by their reverse complement's (`["reverse_complement"]`, `labels` empty) |
| `predicate_unbudgeted_allowed` | the PRIMARY index, `none(573.fa)`, `"context"`, `allow_unbudgeted_annotation` | 12 rows (one per k-mer pair), `access: "columns"`, `predicate.strands: "either"`, note `annotation_unbudgeted` |

### 19.14 Stated limitations

- Context-level, per index (`scope: "shard_context"`): no motif-level predicates (design §12).
- Columns only; record headers are unknown names (P20). Cohorts larger than `caps.max_predicate_labels` are split by
  the client or need a raised cap; no taxonomy terms (P6). On refseq33m a column holds exactly its taxid's records,
  not its descendants: "in E. coli" is a list of every strain taxid.
- `access: "rows"`: a row is decoded whole to read one label of it; the units say so (P17).
- Each selected row is read again by `output.labels: "all"` (P8).
- A predicate request in mode `count` reads annotation (unlike a count without one).
- The pass costs a row per distinct row: about 1–4 ms a row on staging refseq33m (the coordinator's measurement of
  2026-10-08), so its compute admission there is the host's `--pattern-max-predicate-contexts` 10,000 (P3).
- The reverse-complement lookups all come before the first read: `stop_at_threshold` saves reads, not lookups.
- The result objects of `output.labels: "none"` are held by the account for their pattern only (§19.9).

### 19.15 What changed

**Additions only** (contract version 1): the request fields and values of §19.2, the answer fields of §19.10, the
values of §19.8, the notes `predicate_constant` and `projection_not_read`, the refusal code `predicate_too_large`,
the capabilities' additions of §19.12. **Changed for requests that were refused anyway**: `output.labels:
"predicate_only"` without a predicate is 400 `invalid_request` (was `later_increment`); `predicate`,
`max_predicate_contexts`, `max_predicate_work` are read (were `later_increment` by name), so a body naming them
with another fault can now be refused for that fault first; `caps_rule` (prose) names the three caps.

Checked with the identity panel (`round-5bA/route/panel`: every POST fixture of the single-graph mini servers but the
time budgets, and the variants of the 2026-10-07 panel and of round fix3 — 182 requests and both capabilities routes
— on the masked copy and on the mini as built, against the build of `9018f41b`): 181 of 184 identical apart from
`timing` on each graph; the three others are `/capabilities` and `/traverse/capabilities` (the additions of §19.12
and `caps_rule`, nothing else in the document) and `later_increment_labels` (the code above). Rerun on the
integrated build (`round-5bA/integrator/panel`): the same on each graph (beside them the 21 predicate requests,
which the baseline refused); 5s-2's engine panel (every POST fixture, long-pattern shapes, peptides, 294 step
budgets inside extensions) on six server configurations, its 558 requests without a predicate identical on each;
303 `/traverse`, `/resolve` and labelled `/pattern` requests identical apart from `timing` but
`/traverse/capabilities`.

**Fixtures** (§11): 25 new (`predicate_only_without_predicate` among them, replacing `later_increment_labels`), a
server variant `masked_small_predicate_cap`; the nine capabilities bodies of the single-graph servers regenerated; no other stored
body changed. The validator (`test_pattern_fixtures.py`) knows the new tables (`predicate`, `predicate_block`,
`selection`, `capabilities_predicate`) and rows, the new values, and checks the predicate block against its own
folding of the request's predicate, the relations of §19.7, the withheld and cut reasons, every result's
`selection_labels` against its own evaluation of the normal form and its `selection_strands` against its own
labels and the graph's mode, and `"predicate_only"`'s labels against the normal form's names (`test_predicate_rules_refuse_what_v1_never_answers`: mutated answers each rule refuses).

**Tests with independent oracles** (records scanned, a recursive evaluator of the request's JSON, never the
modules under test): `PatternSelection.*` (`tests/cli/test_pattern_selection.cpp`, through the route), among them
every `max_predicate_work` up to the pass's total in three modes and both strand settings, every reading of a
virtual clock, the memory account from the smallest, the descriptors' admission, refused rows, the sticky stop,
CANONICAL and PRIMARY builds, unbudgeted access, a graph without its mask; `PatternRoute.Refusals`, `RefusalOrder`
and `Capabilities`; the integration's `test_predicate_against_the_fasta` (the mini's FASTA scanned per k-mer and
its reverse complement).
