# Design: the verified test and benchmark suite for /traverse, /resolve and /pattern

This is design v2 of 2026-10-13, the plan the suite under `metagraph/scripts/suite/` was built from (lanes S1 to
S5e). It was written as a working note, reviewed by the search-service (async) session, and moved here unchanged
on 2026-10-14. Where the built suite departs from it, the suite's `README.md` is authoritative; the departures as of
2026-10-14 are listed first.

## As built (2026-10-14)

- **Scenarios.** The frozen manifests hold 587 (v1), 679 (v2, v3) and 681 (v4) scenarios, not about 800. v4 by
  area: R 81, T 93, P 97, L 48, X 100, M 62, Q 32, C 12, N 135, D 21. No manifest uses the `long` time class yet
  (v4: 607 fast, 74 slow). v4 also carries relations between scenarios (`relations_sha256`;
  `checks.evaluate_relation`).
- **Area S** (service features, MCP arm only) is built in the search-service repository
  (`scripts/bench/suite_area_s.py` and `suite_area_s.json`, 14 scenarios). The MCP runner runs it on its own: it is
  not part of `manifest/suite_vN.json` or of the paired runs.
- **Check outcomes.** Besides §8.2's pass, fail and not_exposed, a check can be `inconclusive` (it could not test its
  property: an allowed 503, a withheld or cut list that hides the source, a stated resource stop) or `finding` (the
  host broke a stated bound or SPEC promise without failing the request). Executions add `error`,
  `not_implemented_mcp` and `service_by_design`. Neither inconclusive nor finding counts as pass. The rules are in
  the README ("Decisions where the design is silent").
- **`checks.evaluate`** also takes `previous` (the earlier repeats' digests and time limitation, so a time-stopped
  cold repeat followed by a complete warm one is no determinism failure) and `client` (the client wall time, which
  judges an allowed 503 deadline).
- **Comparison classes** (`compare.py`): agree, not_comparable, disagree, not_implemented_mcp, service_by_design
  (verified or not), missing, plus host states (an MCP execution stopped by a service-side limit). The by-design
  table is async's `service_by_design.json` (131 cases), vendored next to the comparator; a case can also be
  matched by the answer, not only named by the scenario.
- **Pairing.** Per scenario: MCP cold, then direct-local, then direct-remote, then MCP warm (§7 says "MCP cold ->
  direct -> MCP warm"; both direct arms sit in the middle). Timings key on the execution's `host_cache`: cold is the
  scenario's first execution on the host in the run, warm another one within 60 s, aged an older one. The host's row
  cache is not otherwise controlled.
- **Offline re-evaluation** (`reevaluate.py`, `run_suite.sh --reevaluate RUN_DIR`): a recorded run judged again with
  the current checks, without any request. It replays the MCP runner's verdict too, and that replay equals async's
  runner on the recorded executions and on mutated ones.
- **MCP authentication.** The runner can call with an API key and records `mcp_auth` (anonymous or api_key); §4
  assumes anonymous calls.
- **Commands (§8.5).** `run_suite.sh` grew options for paired or sequential runs, an MCP command or a recorded MCP
  run, resume, compare and re-evaluate; the README lists them.

---

## Original design: Verified test and benchmark suite for traversal, resolve and pattern search (design v2, 2026-10-13)

v2: async's review of v1 folded in (outcome classes, service features, MCP normalisation, overhead decomposition,
the front cut on every arm, job cleanup, vendoring paths, schedule after the service's I12 deploy) and the
interface contract of §8 added.

Owner's request (2026-10-13): (a) a comprehensive suite testing every functionality and edge case of /traverse,
/resolve and /pattern incl. mismatches, MetaGraph session leading the design, the search-service (async) session
consulted; (b) one reusable script that runs it against refseq33m-experimental with structured output, timings and
resource use, usable as a benchmark; (c) the same tests through the MCP tools: same results? MCP overhead?; (d) a
sharable report, markdown converted to PDF; repeatable with low effort; hours are fine. Owner's decisions:
single-graph, refseq33m only, the deployed staging version (update first if needed); run on mex now despite the
Logan load (timings labelled as such), rerun when mex is free; what the MCP cannot express yet is reported as "not
implemented yet"; queries must be REAL and VERIFIABLE (sequences known to be in the genomes, found back).

### 1. Principles

1. **Every scenario has a known right answer.** Positive controls are cut from records known to be in the index, at
   known positions, so the answer must name the record (its taxid column, its header, its position and strand).
   Variants are made from them with known edits (substitutions at known positions, IUPAC codes at known
   positions, translations in the record's genetic code). Negative controls cannot occur (random or shuffled
   sequences long enough to be absent, a stop codon in tables that have none). Error cases have exact expected
   codes and texts (the fixtures' texts). Everything else is checked by invariants that hold on any answer.
2. **One manifest, three arms.** A frozen, versioned scenario manifest (JSON, content_sha256) is the single
   source. Arms: `direct-local` (the runner on mex against the container's local port: no front proxy, the
   cleanest server timing), `direct-remote` (from the laptop through the public URL and nginx on dnaloc-new: the
   network and front cost), `mcp` (the same request through the staging MCP tools: the service cost). A scenario
   the MCP cannot express states "not implemented yet via MCP" in the report.
3. **Outcome classes** per scenario and arm: `pass`, `fail`, `not implemented yet via MCP` (the MCP cannot express
   the request), `service-side by design` (the service validates or reshapes before the host sees the request, so
   the arms legitimately differ: e.g. 17 patterns answered via MCP in chunks of the host's 16, max_mismatches
   refusals at submit as `pattern_out_of_bounds` / `pattern_unsupported` / `pattern_invalid` instead of the host's
   slot or 400, max_contexts above 10,000 refused, unnamed max_contexts sent as 1,000, patterns over 1,000 chars
   refused, the service's own codes for invalid resolve / traverse inputs; async supplies the mapping table from
   service codes to host expectations, and a mapped pair counts as agreement by design), `not exposed via MCP` (a
   host field, e.g. parts of `timing` / `work`, that the MCP does not surface: recorded, not compared).
4. **Checks are pure functions, shared.** One `checks` module evaluates an answer (direct or MCP, after
   normalisation) against the scenario's expectations; both runners use it, so a disagreement is never a
   difference of checking code.
5. **Measured, not guessed.** Per execution: client wall time, the server's `timing` and `work` (steps, rows,
   annotation units, `memory_bytes`), cold vs warm; container CPU and memory sampled during the run; host load.
6. **Repeatable.** One command builds nothing new (the manifest is frozen) and produces a run directory with
   `report.pdf`; `--compare RUN` shows the change against an earlier run.

### 2. Ground truth

- **Genomes.** Start from async's `genomes.json` + frozen `manifest_v1` (32 genomes across kingdoms, GC classes,
  genetic codes 1/2/3/4/5/6/11/12, organelles, phages; 160 verified proteins with record coordinates and coding
  DNA; built via the MCP NCBI tools, independently re-checked against E-utilities). Pinned copy, not a fork.
- **Record sequences** at the INDEXED accession.version (versions are immutable, so the sequence fetched is the one
  indexed), cached on disk once, sha256 per record in the manifest.
- **Exclusions** (known index facts, stated per scenario so a miss is explained, never silently skipped): the
  last record of a taxid column (the 2026-10-09 defect: often in `.seqs` but unannotated; avoid by construction,
  plus one explicit scenario that documents the defect), windows with N, windows with low complexity where the
  scenario needs a unique placement (stated: such a window may have many placements; checks then require the
  source placement to be among them, not alone), merged taxids (1.6% of columns).
- **Uniqueness.** For placement checks the window must be unique enough: the check is "the source record at its
  position is among the placements"; a separate count check only where the window is long (>= 40 bp) or the
  scenario is about repeats.

### 3. Scenario catalogue (areas, edge cases, what is verified)

IDs `AREA-nnn`; each scenario names its ground truth, its request, its expectations and its MCP mapping.
Counts are estimates for the first manifest (about 800 scenarios, ~1,100 executions per arm with repeats).

#### R: /resolve (~60)
- R-pos: windows of 31, 32, 50, 100, 200, 1,000, 2,000 bp from genome records; explicit labels = the record's
  taxid: every k-mer supported (`kmers_supported == num_kmers`), one run covering the query; with the record's
  header label (accession.version) likewise.
- R-multi: explicit labels = taxid + unrelated taxids: unrelated ones 0 or partial runs; order and dedup.
- R-junction: a query spanning two records of different taxa (concatenation): two runs, the split at the
  junction (+/- k-1).
- R-discover: discover with limit 1/20/100, kind column/header: the source taxid present when not truncated
  (the documented caveat: index order, not best-supported; checked as presence, not rank).
- R-budget: time_budget_ms small: stated stop with `resolved_kmers`, the answer equals the resolve of that prefix
  (re-resolve the prefix and compare); continuation from `remainder_from_bp` completes it.
- R-edge: exactly k bases; k-1 bases (refused); 2,001 bp (refused, `max_query_bp`); N or IUPAC (refused, exact
  text); lowercase (accepted, uppercased); unknown label (refused naming it); duplicate labels (merged); 1,001
  labels (refused); empty; discover + labels (refused); a sequence absent from the graph (0 k-mers present).

#### T: /traverse (~150)
- T-walk: seeds (31-200 bp) from genome records, label-constrained to the record's taxid, direction right / left /
  both, max_extension_bp 100 / 1,000 / 10,000: **the spelled walk must occur in the record** (substring check at the
  expected position, or the walk stops at a stated branch whose both arms carry the label: then each branch
  spelled must occur in some record of that taxid). `complete_to_bp` per arm <= the distance to the record's end.
- T-end: seeds near a record's start / end: the walk stops at the record end with the stated reason.
- T-repeat: seeds in rRNA / IS elements (multi-copy): branching stated, both arms' sequences occur in the record.
- T-labels: constrain vs annotate; quorum; max_labels_per_node; seed_label_kind header vs column; derived labels.
- T-exhaustive: trie / exhaustive with small radius: every leaf's walk occurs in a record of its label.
- T-attempts: attempt ids, ledger bounds (usage.bound, walk_until_ms), stops at the deadline with consistent prefix.
- T-coords: coordinates on walks (`trace` on this coordinate index): positions equal the record positions.
- T-budget: time_budget_ms / max_memory_mb / max_work_units small: stated stops, partial graphlets consistent.
- T-edge: seed shorter than k (refused), seed with N (refused), seed absent from graph (not_in_graph / stated),
  graph/graphs fields on a single graph (400), non-object body (400), max_seeds+1 (refused), duplicate seeds
  merged, unknown strategy key (refused with the field), bounds above the maxima (refused with field/max).

#### P: /pattern, DNA and IUPAC, patterns of at most k (~150)
- P-exact: windows of 12, 16, 20, 25, 31 bp from records: count > 0; all_or_count / partial results include the
  source k-mer at the right offset and strand; labels all: the record's taxid among the labels, the placement
  names the record at the window's position and strand (`record_verified`).
- P-scope: suffix vs any_offset (suffix ⊆ any_offset; a window at a record start: any_offset finds it, the counts'
  `absence_scope` states each); strands plus / minus / both (minus finds the reverse complement's source; a
  palindrome counted once).
- P-iupac: windows with R/Y/S/W/K/M/B/D/H/V/N at known positions: the source still found; contexts ⊇ the exact
  window's; N runs (leading N skipped; internal N cost stated).
- P-modes: count vs all_or_count vs partial on the same pattern: counts equal, list sizes as the mode says;
  max_contexts at, below, above the count: withheld count_above_threshold / cut max_contexts exactly then.
- P-negative: random 31-mers and shuffled 40-mers: exact 0, `retrieval_complete` and absence claims as the SPEC.
- P-floor: patterns below `min_information_bits` (12-mers, AT-rich, N-heavy): slot `information_below_floor` with
  the bits stated; at exactly the floor: answered.
- P-caps: max_anchors, max_steps small: stated stops; time_budget_ms small: stop time with lower counts `at_least`.
- P-unmasked: counting `upper_bound` on this host (no mask in the deployed directory? stated by capabilities):
  bounds with estimate; few-unchecked exact (max_checked_entries).
- P-edge: empty pattern list; 17 patterns (max_patterns 16 +1: refused); bad alphabet (slot `bad_alphabet`);
  lowercase; U for T; unknown field (400); duplicate ids.

#### L: /pattern, patterns longer than k (~80)
- L-paths: windows of 32, 40, 60, 100, 200 bp: `long_search: paths`: a path spelling the window exists; with labels
  all: the record `record_verified` at the window's coordinates; `require_support: record_verified`.
- L-anchors: default anchors-only answer (`withheld: paths_later_increment` absent now? as served), max_anchors
  threshold (anchors_above_threshold), max_paths cut / withheld.
- L-supported: supported_paths for record windows: supported count >= 1, the record among the supporting labels;
  level label_intersection; record_verified; budgets (max_annotation_work, max_memory_mb) with stated stops.
- L-memory: the new label-free paths admission: a large max_paths with a tight max_memory_mb: stop / cut
  max_memory as SPEC §12.1/§28.
- L-repeat: low-complexity long patterns: stated stops, `low_complexity_pattern` note rules.

#### X: peptides (~80)
- X-peptide: manifest_v1 peptide windows (10/20/30 aa) in their genetic codes: found, the record's taxid among the
  labels, record_verified at the coordinates (as async's validation, reused); shifted windows; the table-1 arm
  for non-standard codes (expected misses where a reassigned codon is used).
- X-stop: `*` in tables 1/11 (a stop codon at that position: found where the record has it); tables 27/28/31
  (`no_stop_codon`, exact 0 without search, the documented exception after a budget stop).
- X-junction: windows across one exon junction: expected 0 paths (`expected_miss_exon_junction`, never a miss).
- X-edge: X, B, Z, J residues; genetic_code unknown (400 genetic_code_unknown); peptides of 11 residues (> k/3:
  paths).

#### M: bounded mismatches (~100)
- M-planted-dna: from record windows of 15-31 bp, variants with t = 1, 2, 3 substitutions at known positions;
  query the VARIANT with max_mismatches m >= t: the source context present at level t with `substitutions` at
  exactly the planted positions (expected = the variant's base, found = the record's base, in P's coordinates;
  minus-strand results mirrored as SPEC §21.8); `by_mismatches` keys 0..m present, relations exact on a completed
  search; with m < t: the source absent from the levels (its level is t).
- M-planted-iupac: IUPAC positions + substitutions: a substitution never at N; a base inside the set is no
  substitution.
- M-planted-protein: one substituted residue (a sense codon of another residue) in a peptide window: the source at
  level 1 with `residue_substitution` (residue, expected, found, codon in P's reading direction); cap: protein 1
  on staging after the next restart (3 until then): m = 2 for a peptide -> slot `mismatches_above_cap`.
- M-order: results ordered by level, then node, offset, orientation; partial with max_contexts cutting inside a
  level.
- M-stops: max_steps small: `stop.mismatch_level`, lower levels exact, the stopped level at_least, higher unknown;
  steps_by_mismatches sums to work.steps.
- M-floor: patterns failing the floor at m but passing at m-1: slot with `max_mismatches_within_floor`.
- M-refusals: max_mismatches 4 / -1 / 1.5 / "1" / null (400 invalid_request, exact text); m > 0 with labels all,
  predicate, supported_paths (400 later_increment); L > k with m > 0 (slot mismatches_above_cap cap 0); m = 0
  byte-identical to the request without the field (except timing).
- M-unmasked: the levels' relations on the unmasked host (bounds + estimate where unprovable).

#### Q: predicates and motif (~60)
- Q-pred: windows of records of taxa A, B (A's records only hold window wA): any([A]) selects wA's contexts,
  none([A]) deselects them, all([A,B]) on a window both hold; at_least n; and/or/not; unknown labels listed in
  `unknown_labels`; normal form constant (true/false) reads nothing.
- Q-strands: predicate_strands either vs context: a window taken from the reverse strand of A's record: either
  selects it with `selection_strands: reverse_complement`, context does not.
- Q-motif: predicate_scope motif: decided_by no_instance / every_context / tested_contexts on known cases.
- Q-supported: predicates on supported paths (selection among supported walks; mirror reads under either).
- Q-budget: max_predicate_contexts / max_predicate_work small: stated stops and withheld reasons.

#### C: capabilities and protocol (~20)
- GET /capabilities, /pattern/capabilities, /traverse/capabilities: shapes, the gate block projection, sizes under
  the 32 KiB budget, release/index_fp/server_instance consistent across the three; the body limit stated.
- Protocol: gzip/deflate negotiation, compact JSON, non-object bodies, wrong content type, HTTP methods.

#### S: service features, MCP arm only (~30; async)
get_results paging to completeness (cursor, no duplicate contexts, all_or_count count == rows); cancel_search on a
running pattern job; idempotency (the same idempotency_key twice -> the same search_id); catalog(feature="pattern")
and resolve_scope for refseq33m-experimental; traverse_resolve inline (wait) vs job path giving the same `resolve`;
traverse_results views and cursor; graphlet_* views consistent with graphlet_export (graphlet_sequence(walk w) ==
the exported walk w; graphlet_walks ranks ⊆ export; graphlet_labels counts == export's); register_sequence +
sequence_id as the query (the same answer as the inline sequence); traverse_lock / traverse_stop lifecycle.

#### N: negative, near-miss and hard-region probes, across every route (~120; the owner's addition, 2026-10-13)
How the code behaves where the answer is "nothing" or the input is hostile, for every normal use case:
- **Negatives:** random sequences (31, 40, 60, 100, 200 bp; random 10-30 aa peptides) and composition-preserving
  shuffles of record windows: /pattern exact 0 with `retrieval_complete` and the stated absence scope; /resolve 0
  k-mers present (explicit labels and discovery); /traverse from a seed the graph does not hold: the stated outcome
  (`not_in_graph` / the seed's refusal as the SPEC says), never a crash or a walk.
- **Near misses derived from positives:** record windows with more substitutions than requested, so the source is
  out of reach: t = m+1 .. m+5 substitutions spread or clustered (e.g. 4-8 substitutions in a 31-mer with
  max_mismatches 0-3): the source is NOT among the contexts at any level, the levels' counts are what they are
  (other contexts may exist: never assumed 0), and the answer's relations and order stay valid; the same windows
  with t <= m: the source back at level t (pairs: reached / out of reach). Peptides with 2-4 substituted residues
  under the protein cap. Partly mutated seeds for /traverse and /resolve (a substitution every ~k/2 bases: runs
  break at the substitutions, the runs between them supported; a walk from a mutated seed refused or stated).
- **Hard regions:** low-complexity sequence from the records and synthetic: homopolymers (A^31, A^40, T^100),
  dinucleotide and short tandem repeats ((AT)^n, (CAG)^n), microsatellites, rRNA / IS / transposase windows
  (many copies), telomeric repeats in eukaryotes: /pattern answers with the `low_complexity_pattern` note where the
  SPEC states it, `count_above_threshold` / `anchors_above_threshold` withheld, max_steps / time stops stated;
  /traverse walks in repeats branch or stop with the stated reason within the budget; /resolve on repeats.
- **Every normal use case of an agent:** identify an unknown read (resolve + pattern), find a primer pair's
  binding sites (two primers, both strands), check a peptide in its genome (table of the record), extend a seed to
  its neighbourhood (traverse both arms), ask which taxa hold a motif (predicate), long marker genes (paths /
  supported paths), and their failure modes (wrong genetic code, wrong strand, too short, too degenerate).

#### D: determinism (~20 scenarios repeated 3x)
- The same request three times: identical answers apart from `timing` (and stated non-determinism flags:
  `determinism: time_limited` only where the SPEC allows).

### 4. MCP mapping (confirmed by async, 2026-10-13)

Endpoint: `https://metagraph-async-api-staging.up.railway.app/mcp`, Streamable HTTP, stateless, anonymous (raw
JSON-RPC `tools/call`, no initialize handshake; anonymous caps apply). Client: reuse `pv_mcp.py` / `mcp_e2e.py`'s
JSON-RPC client and transient-error classification (one retry on a dropped connection).

| route | MCP tool | expressible | "not implemented yet via MCP" (today) |
|---|---|---|---|
| /pattern | pattern_search (job: submit, get_status long-polls up to 240 s per call, get_results paged) | patterns, mode, scope, strands, labels, node_labels, max_contexts, time_budget_ms, long_search, max_paths, require_support, genetic_code, occurrences, predicate (+ predicate_strands, predicate_scope, max_predicate_contexts, max_predicate_work), max_labels, max_anchors, supported_paths_level | max_mismatches (service I12: built, deploys with the next bundle, then 0-3), max_steps, max_memory_mb, min_information_bits, max_annotation_work, max_occurrences_per_label, max_labels_per_anchor, stop_at_threshold, allow_unbudgeted_annotation, an explicit graphs list, anything else not listed |
| /resolve | traverse_resolve (job, inline wait <= 45 s, then traverse_status / traverse_results) | sequence / sequence_id, labels, discover (limit, kind), time_budget_ms | — |
| /traverse | traverse_fetch (presets walk / trie + overrides: bounds, direction, exhaustive, labels.max_labels_per_node, seed_label_kind); graphlet_export for the whole graphlet | T-walk, T-end, T-repeat, parts of T-labels and T-exhaustive | the other strategy fields, attempt fields, coordinates/trace on request, protocol cases |
| capabilities | traverse_capabilities, catalog | identity, limits | the raw documents |

Answer normalisation: pattern_search merges counts across the scope's graphs; refseq33m-experimental is single-graph,
so the merged view is the host answer reshaped, not recomputed: compare per pattern the counts with relations,
withheld / stop / cut, the contexts (kmer, offset, strand) or paths (sequence) from get_results, labels and
placements where both list them. traverse_resolve's `resolve` field is the host answer plus the service's stated
limits / stop. traverse_fetch: compare the graphlet from graphlet_export (walks' spelled sequences, outcomes).
MCP rows from get_results carry the service's columns (pattern_id, kmer, instance, offset, strand or orientation by
graph representation, node, row, labels_count, top_label, ...; null keys omitted), not the host's entry objects:
compare contexts by (kmer, offset, strand/orientation, node); full label lists through the REST raw include /
table_data where a check needs them. Host `timing` / `work` fields the MCP does not surface are `not exposed via
MCP`; the MCP arm records what get_status / get_diagnostics expose.
Caching: pattern and resolve are never served from the service's cache; traverse_fetch has a request cache:
`force_fresh: true` on every timed call. The host's row cache is the same for all arms (order matters: see §5).
Queue: the MCP path includes the service queue; the job timestamps (created / started / finished in get_status /
traverse_status) split queue wait from execution, so the overhead is reported both with and without the queue.

### 5. Measurements

Per execution: client wall ms (submit to final answer; for MCP also queue wait and poll overhead separately),
server `timing` (elapsed and phases) and `work` counters, answer bytes; first execution of a scenario = cold,
an immediate repeat = warm (a subset: all of D and every tenth scenario). Resource sampler: `docker stats` of
`metagraph_refseq33m_experimental` every 2 s (CPU %, memory, net/block IO) plus host load average and the
sra_all container's CPU (the Logan load), joined to executions by time. Overhead decomposition. Per MCP execution (async's runner): (a) submit acknowledgement (MCP / HTTP handling),
(b) queue wait (job created -> started), (c) the backend call as the SERVICE measured it (per task sent_at /
elapsed, GET /admin/view/search-timing/{search_id}), (d) finished -> the client sees the terminal status (long-poll
granularity), (e) result paging. Then (c) is compared with direct-remote and direct-local separately (the service
calls the host from Railway, not from the laptop, so its network path differs from direct-remote's), and
direct-remote - direct-local gives the laptop's network + front cost. Pairing MCP cold -> direct -> MCP warm: MCP
warm is compared with direct (both warm); MCP cold is the cold figure. The front proxy on dnaloc-new cuts at
~120 s for EVERY arm that goes through it (the service's base URL for refseq33m-experimental is the public one):
scenarios expected over ~100 s run direct-local only, for all arms.

### 6. Output and report

Run directory `runs/<UTC>_<release>/`: `run.json` (manifest sha, scripts commit, release, server_instance,
index_fp, capabilities snapshot, host load at start/end, arms, options), `executions.jsonl` (one line per
scenario x arm x repeat: ids, status, wall/server times, work, checks [{name, pass, detail}], normalised-answer
digest, raw answer path), `samples.jsonl`, `raw/` (gzipped answers), `jobs.jsonl` (every MCP search_id / job_id with its token, for
cleanup; the runner cancels anything still running when it aborts; anonymous jobs expire 48 h after completion), `report.md`, `report.pdf`.
Report: verdict and environment (release, load caveat); per area pass / fail / not implemented yet via MCP; every
failure with its expectation and the observed value; MCP agreement per scenario; timing tables per area (p50, p90,
max; cold vs warm; per arm); MCP and front overhead distributions; resource use per area (peak container memory and
CPU, server memory_bytes); with `--compare`, the deltas against the earlier run. pandoc + xelatex (installed).

### 7. Implementation and ownership

| part | owner | where |
|---|---|---|
| ground-truth builder (record fetch + cache, window/variant/control generation), frozen manifest | MetaGraph session (agent: Opus writes, Fable reviews) | metagraph repo, `scripts/suite/` |
| `checks` module (pure), answer normaliser for direct answers | MetaGraph session | same |
| direct runner (stdlib, both direct arms), resource sampler, orchestrator (`run_suite.sh`: ssh to mex for direct-local, local arms, merge) | MetaGraph session | same |
| MCP runner + MCP answer normaliser, the service-code mapping table, the outcome marks, area S | async session | search-service repo: `scripts/bench/suite_mcp_runner.py`, `scripts/bench/vendor_suite/` (git archive of a pinned metagraph commit, PIN.md) |
| comparator, report generator (md -> pdf), `--compare` | MetaGraph session | metagraph repo |
| first full run, review of results | both; owner reads the PDF | — |

Order: (1) async reviews this design and answers the MCP questions; (2) builder + checks + direct runner, a
20-scenario smoke run on mex; (3) MCP runner on the same 20; (4) full manifest, full run (hours) AFTER the service's I12 deploy (today), so max_mismatches is expressible via
MCP, announced to async and org-id beforehand; report; (5) the
rerun after the Logan backlog drains becomes the timing baseline.
Pacing: one request in flight per arm, never two arms at once against the host (interleaved per scenario:
MCP cold -> direct identical -> MCP warm, as traversal_mcp_bench.py pairs them), 1-2 s pauses; the service's queue
checked before and during a run (GET /admin/inspect/queue/summary); an hours-long run is announced to async and
org-id beforehand (owner rule: no deploy while staging is busy) and mex's load (the sra_all container) is recorded
next to every timing.

### 8. Interface contract (so the builder, the runners, the comparator and the report can be written in parallel)

#### 8.1 Manifest (`manifest/suite_vN.json`, frozen; `content_sha256` over the canonical JSON of `scenarios`)
```json
{
  "manifest_version": 1, "content_sha256": "...", "created_utc": "...", "index": {"database": "refseq33m-experimental",
   "release_at_build": "refseq97-...", "index_fp": "907ec986..."},
  "ground_truth": {"genomes_source": "search-service manifest_v1 @ <commit>", "records": [{"accession": "NC_...",
   "taxid": "562", "length": 4641652, "sha256": "...", "fetched_from": "efetch"}]},
  "scenarios": [{
    "id": "M-planted-dna-0007", "area": "M", "title": "...", "route": "POST /pattern",
    "request": { "...": "the exact JSON body sent to the host" },
    "mcp": {"tool": "pattern_search", "arguments": {"...": "..."}} | {"not_implemented": "max_mismatches"},
    "truth": {"record": "NC_...", "taxid": "562", "start": 123456, "strand": "+", "window": "ACGT...",
              "edits": [{"pos": 7, "from": "A", "to": "G"}], "...": "area-specific"},
    "expect": [{"check": "contains_source_context", "args": {"level": 2}}, {"check": "status", "args": {"http": 200}}],
    "timing_class": "fast | slow | long (> 100 s: direct-local only)", "repeat": 1, "tags": ["planted", "mismatch"]
  }]
}
```

#### 8.2 Checks (`checks.py`, pure, stdlib only)
`evaluate(scenario, answer, arm) -> [{"check": name, "outcome": "pass|fail|not_exposed", "detail": str}]`, where
`answer` is the NORMALISED answer: `{"http": int, "body": obj | None, "error": {"code", "text"} | None}` for the
direct arms; the MCP arm's normaliser maps its tool results into the same shape for the fields both expose (and
marks the rest `not_exposed`). Each check name is a function in `checks.py`; new checks are added there only.

#### 8.3 Executions (`executions.jsonl`, one line per scenario x arm x repeat)
```json
{"run_id": "...", "scenario": "M-planted-dna-0007", "arm": "direct-local | direct-remote | mcp", "repeat": 0,
 "cache": "cold | warm", "started_utc": "...", "wall_ms": 1234.5, "http": 200, "outcome": "pass | fail |
 not_implemented_mcp | service_by_design | error", "checks": [...], "server": {"elapsed_ms": ..., "timing": {...},
 "work": {...}}, "mcp_phases": {"submit_ms", "queue_ms", "backend_ms", "terminal_lag_ms", "paging_ms"} | null,
 "answer_digest": "sha256 of the normalised answer without timing", "raw": "raw/<scenario>.<arm>.<repeat>.json.gz",
 "job": {"search_id", "job_id"} | null}
```

#### 8.4 Samples (`samples.jsonl`): `{"t_utc", "container": "metagraph_refseq33m_experimental", "cpu_pct", "mem_bytes",
"mem_limit_bytes", "net_rx", "net_tx", "blk_read", "blk_write", "host_load1", "sra_all_cpu_pct"}` every 2 s.

#### 8.5 Commands
- `scripts/suite/build_manifest.py --out manifest/suite_vN.json` (network: E-utilities + the pinned genomes; cached).
- `scripts/suite/run_direct.py --manifest M --arm direct-local --base http://localhost:<port>/... --out RUN_DIR
  [--areas M,P] [--limit N] [--sample-docker]` (stdlib; runs on mex for direct-local).
- `scripts/suite/run_suite.sh [--areas ...] [--smoke]`: rsyncs scripts + manifest to mex, runs direct-local there with
  the sampler, pulls the run back, runs direct-remote locally, triggers (or reads) the MCP arm's results, then
  `compare.py RUN_DIR` and `report.py RUN_DIR` (md, then pandoc + xelatex to PDF).
