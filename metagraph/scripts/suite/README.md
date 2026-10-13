# Verified test and benchmark suite (/traverse, /resolve, /pattern on refseq33m-experimental)

Design: `~/.cache/metagraph-traverse-plans/2026-10-13/suite/DESIGN.md` (v2; §8 is the interface contract).

## Runner, sampler, orchestrator, comparator, report (lane S2)

### Files
- `suite_common.py`: shared stdlib helpers (direct answer normaliser, answer digest, server timing/work
  extraction, execution plan, statistics, jsonl I/O). Copied to mex with the runner.
- `run_direct.py`: runs a manifest against one direct arm (`direct-local` on mex, `direct-remote` from the
  laptop); writes `executions.jsonl`, `raw/`, `run.json`, `state.json`; resumable; stops cleanly on SIGINT.
- `sampler.py`: `docker stats --no-stream` of the staging container and of sra_all (CPU) plus the host load,
  every 2 s, into `samples.jsonl` (started by `run_direct.py --sample-docker`, or on its own).
- `checks_stub.py`: stand-in for lane S1's `checks.py` with the same API (HTTP status only); the runner uses
  `checks.py` whenever it exists, and `run.json` names the module that ran.
- `manifest/runner_smoke5.json`: the runner's own 5-scenario smoke manifest (status checks only).

### Run directory
`runs/<UTC>_<release-short>/` (git-ignored):

    run.json  executions.jsonl  samples.jsonl  jobs.jsonl      merged over the arms (merge_run.py)
    compare.json  report.md  report.pdf
    arms/direct-local/   run.json executions.jsonl samples.jsonl raw/ state.json manifest.json runner.log
    arms/direct-remote/  the same (its samples come from a sampler on mex during the remote phase)
    arms/mcp/            a copy of the MCP runner's run directory (--mcp-run)

Each arm writes only its own sub-directory, so it resumes on its own and a pull from mex never overwrites
another arm's lines. `run_direct.py --out X` on its own writes X/executions.jsonl in the §8.3 schema.

### Decisions where the design is silent
- Executions per scenario: its `repeat` (default 1); repeat 0 is `cold`, later repeats `warm`. A scenario with
  `repeat` 1 whose index in the whole manifest is a multiple of 10 gets one warm repeat (§5 "every tenth
  scenario"; `--warm-every 0` turns it off). The index is global, so a subset run repeats the same scenarios.
- `timing_class` long is not run on direct-remote (the front proxy cuts at ~120 s): no execution line, the
  scenario is listed under `skipped` in the arm's run.json. Client timeouts: fast 120 s, slow 600 s, long
  1800 s; direct-remote at most 125 s.
- Execution lines carry the §8.3 fields plus `area`, `route`, `timing_class`, `ended_utc`, `error`
  ({code, text} of a refused or failed request) and `client` (connect_ms, ttfb_ms, bytes_wire, bytes_body,
  content_encoding, timeout_s). `outcome` is `error` when no HTTP status arrived (http 0, error.code
  `transport` / `timeout`), when the body could not be decoded, or when `checks.evaluate` raised.
- The normalised direct answer is `{http, body, error, headers}`: `error` from the host's `{"code", "error"}`
  body; `headers` keeps content-type and content-encoding for the protocol checks.
- `answer_digest` is the sha256 of the canonical JSON (sorted keys, compact) of the normalised answer without
  `headers` and without every `timing` key; on GET */capabilities also without `measured_*` / `observed_*`
  keys (the delivery calibration and the longest uninterruptible chunk change with the server's life).
- `server` in an execution line: `elapsed_ms` = the top-level `timing.elapsed_ms`, else the sum over the
  entries (patterns / seeds); `work` sums the entries' counters (pattern `work`, traverse per-arm `counters`
  as `counters.*`, the seed's `annotation` block as `annotation.*`; nested objects as dotted keys), with
  `memory_bytes` and `max_*` as peaks.
- Optional scenario fields for protocol cases: `headers` (overrides the defaults; null removes one) and
  `body_raw` (a string sent verbatim instead of the JSON of `request`). Defaults: `Content-Type:
  application/json`, `Accept-Encoding: gzip, deflate`, `Connection: close`, one connection per request.
- run.json keeps the three capabilities documents of the run's start, release / server_instance / index_fp
  (the first non-null value in /capabilities, /traverse/capabilities, /pattern/capabilities), the host load
  at start and end, every session (a resume adds one) and `server_instance_changed` (from a closing GET
  /capabilities: did the server restart during the run). The run directory keeps a copy of the manifest.
- Resume: the same `--out` again skips the executions already in executions.jsonl (a line torn by a kill is
  ignored and its execution runs again); another manifest in the same directory is refused, and so is a
  directory whose state.json names a runner still alive on this host.
- `checks.evaluate(scenario, answer, arm, previous=...)`: `previous` is the list of answer digests
  (`suite_common.answer_digest`) of the scenario's earlier repeats on this arm, in repeat order, also across a
  resume (lane S1's determinism check compares against them, so it must digest with the same function).
- Exit codes of run_direct.py: 0 done, 130 stopped by a signal, 3 stopped at `--max-requests`, 1 failed.

### Orchestrator
`run_suite.sh [--manifest M | --smoke] [--areas ...] [--limit N] [--skip-local] [--skip-remote]
[--mcp-run DIR] [--compare RUN_DIR] [--resume RUN_DIR] [--by-design PATH] [--pause S] [--max-requests N]
[--no-pdf]`; `--help` lists the environment overrides.
1. A new run directory `runs/<UTC>_<release-short>/` (the release from GET /capabilities of the public URL),
   or the one named by `--resume` (then the run's own manifest copy unless `--manifest` is given).
   `orchestrator.json` records every invocation's options and the scripts' commit (`+dirty` when the suite
   directory has uncommitted changes).
2. direct-local: rsync of the suite's `*.py` / `*.json` to mex:~/suite-dev/scripts/ and of the manifest to
   ~/suite-dev/manifest/<run id>.json; `run_direct.py --sample-docker` started there with nohup (the run
   outlives the ssh session) into ~/suite-dev/runs/<run id>/arms/direct-local; polled every 15 s through
   state.json (an ssh that answers neither alive nor dead is retried, never read as finished, so the remote
   arm cannot start while direct-local still runs); then pulled back. A runner of an earlier invocation still
   alive on mex is attached to, not started twice. Ctrl-C sends SIGINT to the mex runner, waits for its clean
   stop, pulls, and prints the --resume command.
3. direct-remote: run_direct.py from here against the public URL, while a sampler on mex records the phase
   (stopped afterwards, at most 30 min for --smoke / 12 h otherwise if this script dies; its samples are
   appended to arms/direct-remote/samples.jsonl). `SUITE_REMOTE_SAMPLER=0` runs it without any mex contact.
4. `--mcp-run DIR`: the MCP runner's run directory (same schema) is copied to arms/mcp.
5. merge_run.py, compare.py, report.py (report.pdf unless --no-pdf).
`--smoke` defaults `--max-requests` to 40 per arm (the smoke load rule) and the manifest to
manifest/suite_smoke_v1.json (lane S1's); without `--manifest` a full run takes manifest/suite_v1.json.

### Comparator and report
- `merge_run.py RUN_DIR`: writes the top-level run.json / executions.jsonl / samples.jsonl / jobs.jsonl /
  manifest.json from `arms/*` (raw paths rewritten to `arms/<arm>/raw/...`; for a key recorded twice the
  later line wins; lines of another arm in an arm's file are dropped with a warning). Re-runnable.
- `compare.py RUN_DIR [--by-design PATH]`: merges first when `arms/` exists, then writes compare.json:
  per scenario and arm the outcome over the repeats (the common one, else error > fail > mixed) and whether
  the repeats' digests agree; direct-local vs direct-remote as `identical` (same digest), `same_checks`
  (other answer, same check outcomes) or `disagree`; the MCP class of DESIGN §1.3; overheads, cold vs warm,
  timing and resource summaries per area and arm.
- MCP classes, in this order: the manifest's `mcp.not_implemented` -> `not_implemented_mcp`; the scenario's
  `mcp.by_design` case(s) in async's table (below); the MCP runner's own `not_implemented_mcp` /
  `service_by_design` outcome; else `agree` when the checks both arms evaluated have the same outcomes (a
  check `not_exposed` on either side is not compared) and neither side alone errored, `disagree` otherwise.
  No MCP line: `missing`.
- `service_by_design.json`: async's table (124 cases), copied byte-identical from the search-service work at
  commit 22d034ba (`service_commit` inside; sha256 7ae2ee30...). A case's `agreement`
  `not_implemented_mcp` -> not implemented yet via MCP; `service_only` -> service-side by design;
  `by_design` -> service-side by design when the pair holds: the direct answer matches
  `direct_expectation` (http, code, text with digits ignored, each where stated; the code may also be a
  pattern slot's) and the MCP executions match `service_answer` (`refused_at_submit` / `relayed`: error.code,
  `text_contains` and `field` where the MCP line states them, a line without an error only when the MCP
  runner marked it service_by_design; `answered` / `clamped`: no error). A pair that does not hold is
  `disagree` with the side (direct / mcp / table for an unknown case id). With a list of cases the first
  whose pair holds decides. The MCP arm's execution lines should carry `error` {code, text, field}.
- Overheads (DESIGN §5): front = direct-remote minus direct-local per scenario (cold: repeat 0; warm: repeat
  1). The arms run one after the other, so the later one finds the host's row cache warm (smoke run: resolve
  345 ms server time on direct-local, 1.1 ms on direct-remote) and the raw wall difference mixes cache and
  network, even going negative; the client-side difference, (wall - server elapsed) remote minus the same
  locally, leaves the server's part out and is the front cost to read. Both are reported, with the server
  elapsed difference; answers without `timing` (errors, capabilities) count in the raw rows only. MCP: only executions that reached the host (outcome
  pass / fail); the MCP warm line (repeat 1, else 0) minus each direct arm's line of the same repeat (else its
  repeat 0), with and without `mcp_phases.queue_ms`, and `backend_ms` minus each direct arm; the phases'
  distributions as given; the cold MCP wall time separately.
- Resources: each execution joins the samples from its start - 1 s to its end + 2 s (docker averages CPU over
  its own window), else the nearest sample within 5 s; per area and arm the peak container CPU and memory,
  the mean sra_all CPU and host load1, and the peak `work.memory_bytes` of the answers.
- `report.py RUN_DIR [--compare EARLIER_RUN_DIR] [--no-pdf]`: report.md (ASCII-safe pipe tables: common
  signs spelled out, other non-ASCII as '?', sequences over 40 bp cut, markdown specials escaped), then
  `pandoc report.md -o report.pdf --pdf-engine=xelatex`. Verdict: INCOMPLETE when an arm did not finish,
  FAIL on any failed / errored execution or any direct / MCP disagreement, else PASS; a stub run says so in
  bold. The load caveat appears when the sra_all CPU mean exceeds 400 % or the load1 mean exceeds 48 (half
  of mex's 96 cores). compare.json is computed when missing.
