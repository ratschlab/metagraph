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
  SPEC-pattern-search §2/§12.2/§21.5, seeded generators, and a port of sdust for the host's
  `low_complexity_pattern` note: checked against the compiled sdust.c on 3,000 strings and against every stored
  /pattern answer of an exact DNA pattern).
- `manifest_lib.py`: the request builders (every default the MCP would set otherwise is named explicitly), the MCP
  mapping of DESIGN §4 (tool and arguments, or `not_implemented` naming the field, with the async session's
  `by_design` case ids of `service_by_design.json`), the scenario constructor of §8.1 (it adds `budget_respected`),
  the relations, and `walk_expectations` (a header-constrained walk's exact minimum and branch depth from its
  record).
- `areas.py`: the scenario catalogue of DESIGN §3, one generator per area: R, T, P, L, X, M, Q, C, D and N (S is
  the async session's).
- `build_manifest.py`: writes `manifest/suite_v2.json` (the full catalogue) and `manifest/suite_smoke_v2.json`
  (the scenarios tagged `smoke`; `suite_v1.json` / `suite_smoke_v1.json` are lane S1's first manifests, kept), one scenario per line; `content_sha256` = sha256 of the canonical JSON of
  `scenarios` (the same seed and records give the same digest); the cross-scenario `relations` (below) one per
  line, digested in `relations_sha256`; `skipped` (planned scenarios a generator could not build, with the
  reason: nothing drops silently) and `not_applicable` (areas the host's capabilities rule out); every expectation
  validated against `checks.CHECK_ARGS`. `--no-fetch` uses the cache only.
- `checks.py`: the checks of DESIGN §8.2, pure and stdlib (it uses `suite_common.py` when present, for the
  runner's normaliser and answer digest, and falls back to the same algorithm): `evaluate(scenario, answer, arm,
  previous=None)`, `scenario_outcome(results)`, and the reference semantics of the relations,
  `evaluate_relation(rel, a_scenario, a_answer, b_scenario, b_answer, arm)`.
- `tests/test_seqtools.py`, `tests/test_checks.py`, `tests/test_manifest_lib.py`: offline (`python3 -m unittest
  discover -s tests`); every check on stored fixture answers of the mini index and on hand-made answers, each
  with a passing case and a failing case on a plausible wrong answer.

### Commands
    python3 truth.py                                    # fetch the record sequences (cached)
    python3 build_manifest.py                           # both v2 manifests (seed 20261013; ~2 min from the cache)
    python3 -m unittest discover -s tests

### Decisions where the design is silent
- Truth taxid = the record's column in the index's `.seqs` (HSV-1 is column 10298, S. pombe 4896, Spiroplasma
  315358, not genomes.json's values). Records over 10 Mb are fetched as genomes.json's region; every coordinate
  stays a whole-record coordinate (`truth.start` 0-based on the forward strand).
- The last record of a column (the 2026-10-09 defect): the A. thaliana chloroplast NC_000932.1 (found missing) is
  used only in the scenario that documents the defect; the other last records of multi-record columns
  (NC_000908.2, NC_002754.1, NC_002163.1, NC_002528.1) were found annotated by the scan (its first 200 bp) and are
  used; single-record columns (0 of 500 missing in the scan) are used. Every scenario on such a record is tagged
  `last-record-verified` and its truth carries `last_in_column` and `last_record_scan` ("ok", or "single-record
  column"), so that the report lists a failure there apart (it could be the defect, not the host's answer).
- Windows: seeded positions, ACGT only, a DUST score below 2 (no obvious low complexity), unique on both strands of
  the fetched sequence (the whole record, or for the records over 10 Mb only within the fetched region: such a
  window may occur elsewhere in its record, which the checks allow, since they ask for the source among the
  answers, never for uniqueness), never two windows within 50 bp of each other. A check asks for every record k-mer holding
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
  never in the smoke manifest), `until-next-deploy` (m = 2 peptides, answered under today's protein cap of 3 and
  refused after that restart), `code-derived` (set by the builder on every exact error text no smoke run exercised,
  and by hand on expectations read from the source: treat a first failure there as a possible suite error),
  `machine-dependent` (time budgets: 200 or 503, stop or not), `expected-miss`, `retention-dependent` (an attempt id
  meets its own 409 within the host's retention of an hour), `last-record-verified` (above), `relation` (a twin
  scenario made for a relation).
- `deterministic` compares the answer's digest (suite_common.answer_digest) with `previous`, the digests of the
  earlier repeats on the same arm; it is left out of the first repeat (below, "Outcomes").
- What the smoke validation of 2026-10-13 changed (laptop -> public URL, 36 requests in two paced runs): short
  patterns are counted in suffix scope (a 12-mer in any_offset scope ran into the 503 deadline at 60 s three times on
  the loaded host; two such scenarios stay as `heavy`, allowing the 503); planted-mismatch windows take the
  shortest length whose ball meets fewer than ~1,000 random contexts (a 20-mer at max_mismatches 2 had 21,563
  level-2 contexts and was rightly withheld); `walk_end` passed an arm a resource stopped before the record's end
  (lane S1b made that inconclusive). The three 503 deadlines of a 12-mer count at 60 s are themselves a host
  observation (SPEC-pattern-search §7.6: a request stopped by time still answers with its counts): the report
  should list them, which the P-heavy scenarios keep exercising.
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
  implemented yet via MCP; `suite_smoke_v1.json` 26 scenarios, 30 executions (e9ec9564...). Lane S1b's
  `suite_v2.json`: 679 scenarios, 719 executions per arm (content_sha256 bdf44c60...), 48 not implemented via MCP,
  22 relations (relations_sha256 cfba509f...), tags last-record-verified 187, code-derived 73, machine-dependent
  15; `suite_smoke_v2.json` 32 scenarios, 36 executions, 4 relations (8bf4e432...). Built twice from the cache:
  identical digests.
- The smoke validation of lane S1b (2026-10-13, laptop -> public URL, 36 paced requests, one in flight): 33
  scenarios that changed or are new (T-end-0001, R-junction-0001, R-budget-0002, R-rc-0001, P-scope-0011,
  P-caps-0004, P-more-0001, L-long-0002/0007, X-peptide10 of tables 2, 4, 6, X-stop of table 2, X-edge B,
  M-order-0001, M-suffix-0001, M-protein of table 4, M-refusal-0012 and its twin, Q-strands-0001,
  Q-supported-0001, N-hard-0001/0003/0011, N-realrep-0007, N-nearmiss2-0003, N-zero-0001/0002, P-scope-0003/0004,
  C-caps-0001..3) and 6 relations: 30 pass, 3 inconclusive, every relation pass. What it changed: the 16 bp
  palindromes are searched in suffix scope (in any_offset a t4 palindrome's list was rightly withheld above 10,000
  contexts after 32 s; both pass in suffix scope); the protein mismatch positives are tagged machine-dependent (a
  table-4 decapeptide at m 1 stopped by its 60 s budget in discovery: inconclusive); R-budget's budgets are 3,000
  and 8,000 ms, but the 10 kb resolve still answered 503 at budget + ~0.8 s (9.7 s for 750 ms, 8.8 s for 8,000 ms:
  an unpolled row batch overran the reserve, which SPEC §4.5 answers with a 503), so the prefix property is
  inconclusive on the loaded host and rests on the prefix_equals relation when a run sees a stop. T-end-0001 (a
  31 bp seed 60 bp from lambda's end) reached the end in 1.9 s; the A^60 max_anchors threshold, the exact junction
  end 71, the table-1-differing peptides and the 16S mirror all held.

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

### Lane S1b: what the verified review findings changed (2026-10-13)
The two Opus verifications of the Sonnet and Fable reviews (VERIFY-sonnet-findings.md, VERIFY-fable-findings.md)
named the findings to fix; the coordinator's lists merged them. Beside the outcomes, `budget_respected`, the
defaults, the relations and the tags above:
- Walks (H1). `walk_in_record` requires the longest path that follows the record to reach a minimum: for walks
  under the record's header on a whole record, `truth.walk[side].min_bp`, computed by the builder from the record
  (the first step whose (k-1)-mer context has a second successor anywhere in the record on either strand, the first
  edge or seed k-mer met again, the record's end, the radius; 31 of 45 such walks must reach the radius), a shorter
  path failing unless a stated resource limit stopped the arm (inconclusive); otherwise min(radius, the record's
  end, the truth flank), a stated branch or re-entry also accepted (other records of a column can branch where
  this record does not). A '-' truth is a seed on the reverse strand. `expect_branching` asserts a split or branch
  event at or before `truth.walk[side].branch_by_bp` (the first second successor on the record's forward strand,
  which includes the divergence of a repeat's copies); the T-repeat windows are chosen so that a copy diverges
  within the walk (M1). `walk_end` is inconclusive on a stated stop before the end; the T-end seeds are 31 bp,
  40-100 bp from the end, under staging's 60,000 ms cap (the MCP refuses a walk budget above its preset row: a
  by-design case). New checks: `traverse_stop` (T-budget: the knob's resource_stop, cap_trigger and limitation;
  a budget may fail the seed, then the walk check is inconclusive), `labels_per_node`, `label_dict_contains`
  (seed_label_kind header vs column on one seed of a multi-record column, A3), `paths_end` (the quorum); the
  traverse invariants add SPEC §7.0 (complete exactly at the radius, walks complete only with every arm complete).
- Pattern and resolve (H2, M7, L2, S1, S2, A1, A2). A stopped resolve must show its prefix fully supported by the
  record's label (a stop before the first k-mer is inconclusive; R-budget runs at 300 and 750 ms); the hard
  regions assert exact counts within 2 (k - L + 1) 4^(k - L) and `low_complexity_note` (sdust) instead of "min 0";
  A^40 and (AT)^16 are asked as patterns longer than k; `stops_stated` (P-caps-0004), `results_levels` (M-order);
  a count a budget stop left inexact is inconclusive, not fail; Q-strands' windows are 16S 20-mers whose 42 bp the
  record also holds reverse-complemented, and the mirror check fails without a mirror; palindromes are 16 bp (a
  14-mer's any_offset list is rightly withheld above the cap).
- Ground truth (M4, M2, M6). Peptide positives in every code the genomes carry (1, 2, 3, 4, 5, 6, 11, 12; 5-8
  decapeptides and 3-5 long peptides each), for codes other than 1 and 11 windows that table 1 reads differently
  first (table 12's five proteins hold no CTG in any window or shift: stated in `skipped`); an X-stop per code with
  stops; B, Z, J, X edges; M-protein positives in table 4. Planted windows from low-carrier records. R-junction
  expects the exact run ends from the records (A's run continues while B's first bases repeat A's continuation).
- Coverage (M8/M5 as verified; toward the design's ~800: 679 scenarios). Not built: P-unmasked and M-unmasked (the
  host is masked: manifest `not_applicable`), the refuted items (left-arm repeats exist, A^31 cannot be withheld),
  a body above max_request_body_mb (256 MB would load the shared host), P-modes "at the count" (the count is not
  known before the run), a deadline stop of an attempt (needs a run-time not_after_ms), max_live_paths (needs a
  branch at a known depth under a header label; the repeat windows give one, left for the next increment).
- Not changed (judged a risk to break the MCP arm): the walk preset's memory. The direct bodies send 512 MB, as the
  live service capabilities of 2026-10-13 stated; the service's code default is 256 MB and the MCP refuses an
  override above its preset row, so naming 512 in the traverse_fetch overrides would refuse every MCP walk if the
  row were 256. The run should pin the service's capabilities (run.json) and say which applied.
