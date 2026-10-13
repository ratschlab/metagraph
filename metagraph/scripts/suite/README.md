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
  `scenarios` (the same seed and records give the same digest). `--no-fetch` uses the cache only.
- `checks.py`: the checks of DESIGN §8.2, pure and stdlib (it uses `suite_common.py` when present, for the
  runner's normaliser and answer digest, and falls back to the same algorithm): `evaluate(scenario, answer, arm,
  previous=None)`.
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
  predicate `max_predicate_contexts`, in partial with labels `max_labels`); explicit-label resolves name
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
- Manifests: `suite_v1.json` 587 scenarios, 627 executions per arm (content_sha256 267fb2b4...), 47 of them not
  implemented yet via MCP; `suite_smoke_v1.json` 26 scenarios, 30 executions (e9ec9564...).
