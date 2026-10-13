# Verified test and benchmark suite (/traverse, /resolve, /pattern on refseq33m-experimental)

Design: `~/.cache/metagraph-traverse-plans/2026-10-13/suite/DESIGN.md` (v2; §8 is the interface contract).

## Ground truth, manifests, checks (lane S1)

### Files
- `truth.py`: the 32 genomes of the search service's validation suite (`truth/search_service/genomes.json`,
  `manifest_v1.json`, pinned at search-service 20343460), what the index holds of each record (column, seq_id,
  column size, the 2026-10-09 last-record scan: `truth/index/`), and the record sequences fetched at their indexed
  accession.version from E-utilities (efetch FASTA, at most ~2.5 requests/s) into `truth/cache/` (git-ignored) with
  their sha256. `python3 truth.py` fetches what is missing.
- `seqtools.py`: pure sequence helpers (IUPAC, NCBI genetic codes, translation, the information bits of
  SPEC-pattern-search §2/§12.2/§21.5, seeded generators).
- `manifest_lib.py`: the request builders (every default the MCP would set otherwise is named explicitly), the MCP
  mapping of DESIGN §4 (tool and arguments, or `not_implemented` naming the field, with the async session's
  `by_design` case ids of `service_by_design.json`), the scenario constructor of §8.1.
- `areas.py`: the scenario catalogue of DESIGN §3, one generator per area: R, T, P, L, X, M, Q, C, D and N (S is
  the async session's).
- `build_manifest.py`: writes `manifest/suite_v1.json` (the full catalogue) and `manifest/suite_smoke_v1.json`
  (the scenarios tagged `smoke`), one scenario per line; `content_sha256` = sha256 of the canonical JSON of
  `scenarios` (the same seed and records give the same digest); the cross-scenario `relations` (below) one per
  line, digested in `relations_sha256`. `--no-fetch` uses the cache only.
- `checks.py`: the checks of DESIGN §8.2, pure and stdlib (it uses `suite_common.py` when present, for the
  runner's normaliser and answer digest, and falls back to the same algorithm): `evaluate(scenario, answer, arm,
  previous=None)`, `scenario_outcome(results)`, and the reference semantics of the relations,
  `evaluate_relation(rel, a_scenario, a_answer, b_scenario, b_answer, arm)`.
- `tests/test_seqtools.py`, `tests/test_checks.py`: offline (`python3 -m unittest discover -s tests`); every check
  on stored fixture answers of the mini index and on hand-made answers, passing and failing.

### Commands
    python3 truth.py                                    # fetch the record sequences (cached)
    python3 build_manifest.py                           # both manifests (seed 20261013)
    python3 -m unittest discover -s tests

### Decisions where the design is silent
- Truth taxid = the record's column in the index's `.seqs` (HSV-1 is column 10298, S. pombe 4896, Spiroplasma
  315358, not genomes.json's values). Records over 10 Mb are fetched as genomes.json's region; every coordinate
  stays a whole-record coordinate (`truth.start` 0-based on the forward strand).
- The last record of a column (the 2026-10-09 defect): the A. thaliana chloroplast NC_000932.1 (found missing) is
  used only in the scenario that documents the defect; the other last records of multi-record columns
  (NC_000908.2, NC_002754.1, NC_002163.1, NC_002528.1) were found annotated by the scan and are used; single-record
  columns (0 of 500 missing in the scan) are used.
- Windows: seeded positions, ACGT only, a DUST score below 2 (no obvious low complexity), unique on both strands of
  the fetched record, never two windows within 50 bp of each other. A check asks for every record k-mer holding
  the window (truth.flank carries the bases to recompute them), never a total count unless it is a lower bound.
- Expectations are what the deployed code answers per the SPECs and the C++ source; error texts are quoted from
  `src/` and the fixtures. Where the design text and the source differ, the source wins: duplicate resolve labels
  are refused (`Duplicate label '<name>'`), not merged; /resolve does not refuse N or IUPAC codes (the k-mers
  holding one are absent from the graph); the host has no cap on explicit resolve labels; the resolve server
  limit is 10,000 bp on staging.
- Both arms ask the same question: every /pattern body names `mode`, `max_contexts`, `time_budget_ms` (and with a
  predicate `max_predicate_contexts`, in partial with labels `max_labels`), and the defaults wherever they act
  (lane S1b, decision 6), in the body as in `mcp.arguments`: `scope`, `strands` (the MCP's plus / minus for
  forward / reverse) and `output.labels` always, `output.occurrences` with labels "all" or "predicate_only",
  `predicate_strands` with a predicate. Where they do not act they stay unnamed: `output.occurrences` with labels
  "none" is a host 400 when true and a service refusal for any value (async: pattern.occurrences_without_labels),
  `predicate_strands` without a predicate a service refusal (pattern.predicate_options_without_predicate) and a
  no-op on the host; a scenario that tests such a refusal names the field itself; explicit-label resolves name
  `bounds.time_budget_ms` 20,000 and `bounds.max_query_bp` 2,000 (500 for a discovery), as the service sends them;
  /traverse bodies are the MCP preset's merged strategy (walk: 15,000 ms, 512 MB, 1e8 work units, 5,000 bp,
  500,000 output bp, 2 label branches; trie: exhaustive, 30,000 ms, 512 MB, 1,000 live paths, 2,000 bp), with
  `output.detail` full on the direct arm (checks.py reads a graphlet body as well).
- MCP blocks may carry `by_design` (async case ids) and `requires_service: "I12"` (max_mismatches through the MCP
  needs the service's I12 deploy). Values the MCP coerces or drops (max_mismatches "1", true, null, 1.5) are
  `not_implemented` naming the field: their host refusals cannot be asked through it.
- Scenario fields beyond §8.1: `body_raw` and `headers` (lane S2's runner) for protocol cases; `truth.copies` for
  repeated seeds; `derived_from` on D scenarios (copies of other areas' fast scenarios, repeat 3).
- Tags: `smoke` (the smoke manifest), `needs-next-deploy` (the protein mismatch cap of 1 after the next restart;
  never in the smoke manifest), `code-derived` (an expectation read from the source the smoke run did not exercise),
  `machine-dependent` (time budgets: 200 or 503, stop or not), `expected-miss`, `retention-dependent` (an attempt id
  meets its own 409 within the host's retention of an hour).
- `deterministic` compares the answer's digest (suite_common.answer_digest) with `previous`, the digests of the
  earlier repeats on the same arm; the first repeat passes with its digest stated.
- What the smoke validation of 2026-10-13 changed (laptop -> public URL, 36 requests in two paced runs): short
  patterns are counted in suffix scope (a 12-mer in any_offset scope ran into the 503 deadline at 60 s three times on
  the loaded host; two such scenarios stay as `heavy`, allowing the 503); planted-mismatch windows take the
  shortest length whose ball meets fewer than ~1,000 random contexts (a 20-mer at max_mismatches 2 had 21,563
  level-2 contexts and was rightly withheld); `walk_end` passes an arm a resource stopped before the record's end
  (status truncated, a resource path reason: a stated stop) and says so, the walk so far being `walk_in_record`'s.
- Outcomes (lane S1b, the integrator's decision 1; the verification's H3 and M7). A check answers `pass`, `fail`,
  `not_exposed` (the MCP arm lacks the field), `inconclusive` (the check could not test its property on this
  answer: an HTTP status the scenario allows that is not a 200, such as a 503 `deadline`; a withheld, cut or
  stopped list that does not show the source, in contains_source_context, source_label, source_placement,
  path_spells, path_record_verified, and the absence checks source_absent, path_absent and predicate_selected;
  a walk a stated resource limit stopped before its property) or `finding` (budget_respected).
  `checks.scenario_outcome(results)` gives the execution's outcome by the precedence fail > finding >
  inconclusive > pass, and `inconclusive` when there is no check or every check is not exposed (the rule lane
  S2's comparator and report and the async MCP runner share). So an allowed 503, a withheld answer that hides the
  source, a walk stopped by a stated resource limit before its property and an MCP execution whose checks are all
  not exposed are inconclusive, never a pass. `deterministic` is left out of a scenario's first repeat (nothing to
  compare yet), so that repeat's outcome rests on the other checks.
- `budget_respected` (lane S1b, the integrator's decision 2, derived per route from the SPECs and each answer's
  own deadline record after the two verifications of 2026-10-13; no fixed allowance anywhere, since the
  capabilities state no bound on one uninterruptible piece, `deadline_check.max_uninterruptible_ms` null): every
  scenario that names a time budget and allows a 200 carries it, its args stating the rule. /pattern and /resolve:
  a 200 whose `timing.elapsed_ms` exceeds the effective `limits.time_budget_ms` broke the SPEC (past the budget
  they answer 503; the deployed /resolve sets `timing.elapsed_ms` before its last deadline check). /traverse, per
  walked seed: the deadline is the budget on the seed's clock (walker.cpp `time_exceeded`) and never stops the
  seed phase (validation or derivation, depth 0, paced against no deadline); so a time stop must come within the
  seed's own longest uninterruptible piece (`timing.deadline.longest_piece`) after max(budget, `seed_phase_ms`)
  (`stop_after_deadline_ms` is measured from the budget), and a seed phase that outlasted the budget leaves one
  level (`complete_to_bp` <= 1). A breach is the check outcome `finding`; a stop late within its stated piece, a
  seed phase longer than the budget, or a walk past the budget without a time stop (a level at the radius is never
  stopped) is a performance fact stated in the detail, for the report's timing section. The smoke run's lambda
  walk (34,201 ms under 15,000: seed phase 33,803 ms with one cold 20,556 ms coordinate-tuple read; stopped
  19,201 ms after the deadline, 398 ms after the seed phase) passes with those facts stated.
- `deterministic` is inconclusive when this answer or an earlier repeat was limited by a time budget
  (SPEC-labeled-traversal-core §6.8: determinism provided no time budget tripped; SPEC-pattern-search §7.9:
  `determinism: time_limited`); `previous` may give each earlier repeat as a digest or as
  `{"digest", "time_limited"}` (lane S2's runner passes digests today: an earlier time-limited repeat is then seen
  only through the digest difference, which this answer's own time stop explains when it has one).
- Manifests: `suite_v1.json` 587 scenarios, 627 executions per arm (content_sha256 267fb2b4...), 47 of them not
  implemented yet via MCP; `suite_smoke_v1.json` 26 scenarios, 30 executions (e9ec9564...).

### Relations (comparisons across scenarios; lane S1b, the integrator's decision 8)
The manifest's top-level `relations` list (outside `content_sha256`, digested in `relations_sha256`; the smoke
manifest keeps the relations whose two scenarios it holds). Each relation is
`{"id": "REL-<type>-nnnn", "type", "a": scenario id, "b": scenario id, "title", "args"}`. The comparator evaluates
each relation per arm on the normalised answers of a and b on that arm (the first repeat of each), with the outcome
`pass`, `fail`, `inconclusive` (the pair cannot test the claim: a non-200 answer where the type needs two answers,
no common k-mer, a count that is not exact, b's list incomplete where a result of a is missing from it) or
`not_exposed` (the MCP arm lacks a field the type reads). `checks.evaluate_relation(rel, a_scenario, a_answer,
b_scenario, b_answer, arm)` is the reference implementation (tests/test_checks.py `Relations`). Paths are dotted
(`patterns.0.counts.contexts`), `""` is the whole body; "without timing" removes every key named `timing` at any
depth.

- `same_answer` — the two requests must be answered alike. Both answers must have the same HTTP status; if it is
  not 200, the same error `code` and `text`. With `args.paths` (a list of `[path_in_a, path_in_b]`), each pair
  of values must be equal (without timing); a missing path fails. Without `paths`, the two bodies must be equal
  without timing and without the dotted paths in `args.ignore` (default none). Used for: max_mismatches 0 versus
  the same request without the field (SPEC-pattern-search §21.2), lower case versus upper case (pattern and
  resolve), and the three capabilities documents (release, server instance, deadline-check and resolve blocks of
  /capabilities and /traverse/capabilities; /capabilities' `pattern` block equal to /pattern/capabilities; the gate
  block of /traverse/capabilities a projection of /pattern/capabilities, its fields listed from the snapshot).
- `prefix_equals` — two /resolve answers of pieces of one sequence S agree on their common k-mers: a's query is
  S from k-mer `args.a_from_kmer`, b's from `args.b_from_kmer` (k-mer i starts at base i). An answer covers the
  k-mers [from, from + num_kmers) of S (num_kmers is the resolved prefix when the answer states a stop). On the
  overlap [max(fa, fb), min(fa + na, fb + nb)), the k-mers inside `graph_runs` and, for each label in
  `args.labels`, the k-mers inside its `runs` (all shifted into S's k-mer coordinates) must be the same sets.
  Inconclusive when either answer is not a 200 or the overlap is empty. Used for: a budget-stopped resolve versus
  the same sequence under the 20,000 ms budget (the stop answers the resolve of a prefix, SPEC-labeled-traversal
  §4.5), and the resolve from base 5,000 versus the whole (a continuation completes it).
- `zero_equals_absent` — a /pattern entry's exact count of an exact DNA pattern of L >= k bases searched
  `strands: "forward"` (a) and the /resolve of the same sequence (b) agree on absence: `counts.<args.count>` of
  entry `args.pattern` (contexts for L = k, paths for L > k with long_search paths) is exact 0 exactly when b's
  `graph_runs` is not `[[0, L - k + 1]]` (some k-mer of the pattern is not in the graph), and exact >= 1 exactly
  when it is (on a node-centric de Bruijn graph consecutive k-mers present spell the path). Inconclusive when the
  count is not exact, the entry is refused or either answer is not a 200. Used for random and record sequences of
  31 and 60 bp and a record 31-mer with one substitution (whatever the graph holds, the routes agree).
- `subset` — every result of a's entry `args.a_pattern` (contexts keyed by kmer, offset and strand or
  orientation; paths, `args.key: "path"`, by sequence and strand), restricted to `mismatches == args.a_level` when
  given, is among the results of b's entry `args.b_pattern` (restricted to `args.b_level` when given); and when
  both counts (`counts.contexts`, else `counts.paths`; the level's `by_mismatches` count when a level is given)
  are exact, a's is at most b's. A result of a missing from b fails when b's `retrieval_complete` is true, else
  is inconclusive. Used for: suffix scope within any_offset, an exact window within its IUPAC pattern, a partial
  list cut at 2 within the complete list, the exact answer within the level 0 of max_mismatches 1.
