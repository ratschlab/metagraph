# Verified test and benchmark suite (/traverse, /resolve, /pattern on refseq33m-experimental)

Design: [`metagraph/docs/DESIGN-test-suite.md`](../../docs/DESIGN-test-suite.md) (v2, with how the built suite departs from it; §8 is the interface contract).

## Runner, sampler, orchestrator, comparator, report (lane S2)

### Files
- `suite_common.py`: shared stdlib helpers (direct answer normaliser, answer digest, server timing/work
  extraction, execution plan, statistics, jsonl I/O). Copied to mex with the runner.
- `run_direct.py`: runs a manifest against one direct arm (`direct-local` on mex, `direct-remote` from the
  laptop); writes `executions.jsonl`, `raw/`, `run.json`, `state.json`; resumable; stops cleanly on SIGINT.
- `sampler.py`: `docker stats --no-stream` of the staging container and of sra_all (CPU) plus the host load,
  every 2 s, into `samples.jsonl` (started by `run_direct.py --sample-docker`, by run_suite.sh, or on its own;
  it stops on a signal, a stop file, --max-seconds, a gone --parent-pid, or a --heartbeat file untouched for
  --heartbeat-max-age seconds).
- `checks_stub.py`: stand-in for lane S1's `checks.py` with the same API: it evaluates `status` as S1's check
  does (http, code, text, text_prefix, text_contains) and reports every other check `not_exposed`. The runner
  uses `checks.py` whenever it imports; it refuses the stub for a selection that expects any other check unless
  `--allow-stub`, and `run.json` names the module that ran (and the import error when auto fell back).
- `manifest/runner_smoke5.json`: the runner's own 5-scenario smoke manifest (status checks only, with the
  error code and text of the refusal); P-smoke-0001 (a /pattern count) is mapped to pattern_search, the
  other four are marked not implemented via MCP.
- `tests/fake_host.py`, `tests/fake_mcp.py`: a local stand-in for the host (HTTP/1.1 keep-alive) and for
  async's MCP runner, for the offline tests.

### Run directory
`runs/<UTC>_<release-short>/` (git-ignored):

    run.json  executions.jsonl  samples.jsonl  jobs.jsonl      merged over the arms (merge_run.py)
    compare.json  report.md  report.pdf
    arms/direct-local/   run.json executions.jsonl samples.jsonl raw/ state.json manifest.json runner.log
    arms/direct-remote/  the same, without samples (the mex samplers cover it)
    arms/mcp/            the MCP runner's executions.jsonl / jobs.jsonl (--mcp-cmd) or a copy (--mcp-run)
    samples/             the orchestrator's samplers on mex (paired-*.jsonl, direct-remote-*.jsonl)
    orchestrator.json    every invocation of run_suite.sh

Each arm writes only its own sub-directory, so it resumes on its own and a pull from mex never overwrites
another arm's lines. `run_direct.py --out X` on its own writes X/executions.jsonl in the §8.3 schema.

### Decisions where the design is silent
- Executions per scenario: its `repeat` (default 1); repeat 0 is `cache: cold`, later repeats `warm` (the
  arm's own repeat mark of §8.3; the host's cache state is `host_cache`, below). A scenario with
  `repeat` 1 whose index in the whole manifest is a multiple of 10 gets one warm repeat (§5 "every tenth
  scenario"; `--warm-every 0` turns it off). The index is global, so a subset run repeats the same scenarios.
- `timing_class` long is not run on direct-remote (the front proxy cuts at ~120 s): no execution line, the
  scenario is listed under `skipped` in the arm's run.json. Client timeouts: fast 120 s, slow 600 s, long
  1800 s; direct-remote at most 125 s.
- Execution outcome (review of lane S1; the MCP runner applies the same rule): checks may report `pass`,
  `fail`, `not_exposed`, `inconclusive` (the scenario could not test its property: an allowed 503, a withheld
  answer hiding the source, a stated resource stop) or `finding` (the host broke a stated bound or SPEC
  promise without failing the request, e.g. a walk past its time budget); any other outcome reads as a
  failure. The execution is `fail` when a check failed, else `finding` when one is a finding, else
  `inconclusive` when one is inconclusive or none could test anything (no check, or every one not_exposed),
  else `pass` (suite_common.outcome_from_checks). Over a scenario's repeats: error > fail > finding >
  inconclusive. Neither finding nor inconclusive ever counts as pass; the report lists every finding in its
  own section and counts inconclusive apart, and agreement skips a check that is inconclusive on either side.
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
  application/json`, `Accept-Encoding: gzip, deflate`, HTTP/1.1 keep-alive.
- Connection (review M1): one keep-alive connection per runner invocation, reused while the server keeps it
  open; `wall_ms` runs from sending the request to the last body byte on an established connection, and the
  set-up of a fresh connection (TCP, and TLS on the public URL) is `client.connect_ms` (null when reused;
  `client.reused`, `client.attempts`). A request on a reused connection the server had closed meanwhile is
  sent once more on a fresh one (the first attempt never reached the application), but only when it failed
  within 1 s of being sent (`STALE_RESEND_MAX_S`; review N10): a later drop may follow a request the server
  was already processing, and is recorded as a transport error. Reusing across paired steps
  would not keep the TCP window either (Linux restarts slow start after an idle second), so a fresh connection
  per paired step with the connect left out measures the same front cost.
- run.json keeps the three capabilities documents of the run's start, release / server_instance / index_fp
  (the first non-null value in /capabilities, /traverse/capabilities, /pattern/capabilities), the host load
  at start and end, every session (a resume adds one) and `server_instance_changed` (from a closing GET
  /capabilities: did the server restart during the run). The run directory keeps a copy of the manifest.
- Resume: the same `--out` again skips the executions already in executions.jsonl (a line torn by a kill is
  ignored and its execution runs again); another manifest in the same directory is refused, and so is a
  directory whose state.json names a runner still alive on this host.
- `checks.evaluate(scenario, answer, arm, previous=..., client=...)`: `previous` lists the scenario's earlier
  repeats on this arm, in repeat order, also across a resume, as `{"digest", "time_limited"}` (lane S3, the
  verification's L4: the answer digest `suite_common.answer_digest`, which lane S1's determinism check computes
  with the same function, and what made that answer time-limited, `checks.time_limited`, recorded in the
  execution line's `time_limited` and read from the raw answer for a line of before the field), so that a
  time-stopped cold repeat followed by a complete warm one is not a determinism failure; `client` is
  `{"wall_ms"}` of this execution (budget_respected judges an allowed 503 `deadline` by it), passed when the
  module's evaluate takes it. The execution's outcome is `checks.scenario_outcome(checks)` (the stub's runs:
  suite_common.outcome_from_checks).
- Exit codes of run_direct.py: 0 done (or partial: a paired step), 130 stopped by a signal, 3 stopped at
  `--max-requests`, 1 failed. `--max-requests` counts every request of the invocation, the closing GET
  /capabilities included (review L4).
- Paired steps (review H1): `--only IDS` executes just these scenarios of the selection (the three capabilities
  documents only when the arm has none yet, no closing GET); `--finish` executes nothing and takes the closing
  GET; run.json `planned` / `state` always cover the whole selection (`done` when every planned execution is
  recorded, else `partial`). `--status-line` prints a last line `RUN_DIRECT_STATUS {state, executed, requests,
  reason}` for the orchestrator.

### Orchestrator
`run_suite.sh [--manifest M | --smoke] [--paired | --sequential] [--mcp-cmd TEMPLATE | --mcp-run DIR]
[--areas ...] [--limit N] [--skip-local] [--skip-remote] [--resume RUN_DIR] [--compare RUN_DIR]
[--by-design PATH] [--pause S] [--max-requests N] [--warm-every N] [--checks auto|checks|stub] [--allow-stub]
[--no-pdf]`; `--help` lists the environment overrides. `run_suite.sh --reevaluate RUN_DIR` re-evaluates a recorded run
offline instead ("Offline re-evaluation of a recorded run", lane S4).
1. A new run directory `runs/<UTC>_<release-short>/` (the release from one GET /capabilities of the public
   URL, logged, sent only once the lock of item 3 is held and no other runner is on mex: review N4), or the
   one named by `--resume` (then the run's mode and manifest copy unless given).
   `orchestrator.json` records every invocation: options, mode, the scripts' commit (`+dirty` when the suite
   directory has uncommitted changes), the checks module, the mex clock offset (mex's `date +%s.%N` against
   the midpoint of two laptop readings; python only where date lacks %N: review N8), the MCP command TEMPLATE (never
   an expanded command or an environment value), and the loop's state (`loop`: running / done / capped /
   failed / interrupted, its reason when it stopped early, scenarios, completed, direct requests).
2. Checks (review M3; lane S3: chosen with `--checks`): `auto` (the default; `--resume` keeps the run's
   choice) runs checks.py when it is next to the script, otherwise the stub; `checks` requires checks.py; `stub`
   runs the stub even beside checks.py (the offline tests of the stub's rules choose it so). The stub is
   refused before anything runs when the selected scenarios expect any check but `status`, unless
   `--allow-stub`; the runners are passed the choice (`run_direct.py --checks`). The MCP runner judges with
   checks.py of `--checks-dir` whatever the choice. compare.py evaluates the relations with checks.py's
   `evaluate_relation` when it is next to it (`compare(..., relation_evaluator=None)` picks the built-in
   evaluation). The whole suite directory (minus runs/, tests/, __pycache__) is rsynced (with --delete) to
   mex:~/suite-dev/scripts/, so checks.py's sibling modules and data reach mex; the manifest to
   ~/suite-dev/manifest/<run id>.json.
3. Exclusion (review M4): every invocation takes `~/suite-dev/LOCK` on mex before any request to the host
   (an atomic noclobber create holding pid, host, start time, a token; the run id is written in once a new
   run is named); a lock of a process of this machine that is gone, or whose pid now belongs to a younger
   process, is removed; a live one or another machine's is refused. The lock names the machine by a stable id (review N6: `SUITE_HOST_ID`, else macOS's
   LocalHostName from scutil, else the node name; the network name, which a DHCP lease or a VPN reconnect
   changes, is recorded apart as `hostname`), so this machine's own stale lock is recognised after a
   reconnect. A stale lock is removed only while it still holds what was judged stale, the compare and the
   removal inside a mkdir mutex `~/suite-dev/LOCK.break`.
   The lock is released on exit (only if it is still this invocation's). Then any `run_direct.py` under
   ~/suite-dev (pgrep on mex) is waited for (SUITE_WAIT_RUNNER_S, 150 s) and refused if it stays, except a
   sequential --resume attaching to its own direct-local runner. An invocation refused before any arm ran
   leaves no run directory behind. All ssh and rsync calls of an invocation share its own ControlMaster
   connection (ControlPath `<SUITE_SSH_CONTROL_DIR>/%C-<pid>`, ControlPersist), closed on exit; another
   invocation's exit never closes it (review N1: a shared socket let a refused second invocation end the
   holder's master and the ssh step in flight).
4. `--paired` (the default; review H1): per scenario in manifest order, MCP cold (`--mcp-cmd`, repeat 0) ->
   direct-local (one `run_direct.py --only ID` step on mex over the shared ssh connection, its stderr into
   arms/direct-local/runner.log there) -> direct-remote (a `--only ID` step from here) -> MCP warm (each
   further repeat), one at a time; then a `--finish` step per direct arm (closing GET /capabilities, the
   arm's state over the selection). The plan (`suite_common.py paired-plan`) skips what each arm recorded,
   so --resume continues where the loop stopped; the long class runs on direct-local only (the front cuts
   every other arm); a scenario the manifest marks not implemented via MCP gets the MCP repeat 0 only (the
   runner records the mark). One sampler on mex covers the loop (`samples/paired-<stamp>.jsonl`), stopped at
   the end, after 24 h (30 min for --smoke), or when the orchestrator stops touching its heartbeat file for an
   hour (a background keeper touches it every `SUITE_HEARTBEAT_S`, 300 s, while the orchestrator lives, so one
   long step does not end the sampling: review N7). Three failing steps of one direct arm in a row stop the loop (exit 1 after the report). The MCP
   arm follows async's exit codes: the loop stops on the runner's exit 2 / 3 (a start refusal: key file,
   queue wait ran out) or 130 (interrupted), or after three MCP executions in a row that ended `error` (a
   runner-level failure) or were not recorded (`loop.reason` mcp_exit_N / mcp_interrupted / mcp_errors);
   never on `fail` (exit 1 also means a failing scenario, a result). After an execution that ended `error`
   with its job cancelled (jobs.jsonl `canceled` / `cancel_failed`, no `terminal`), the next step waits
   until the cancellation plus the scenario's time budget (`time_budget_ms` of the MCP arguments, else the
   direct request; else 60 s for /pattern, else the timing class's timeout) and 2 s, since the cancelled
   call may still run on the host (suite_common.py mcp-step); the sequential mode waits the same way
   between its cold and warm MCP passes. An ssh failure on a
   direct-local step waits until that step's runner is gone before the next arm starts (review N3): the
   check on mex answers with an explicit marker, so a failing ssh never reads as "no runner"; it is retried
   every 5 s for SUITE_WAIT_RUNNER_S. A runner still there after that, or a mex that stays silent, stops
   the loop with an error (`loop.reason` runner_still_running / mex_unreachable, exit 1 after the report),
   and then no closing step runs, so no request goes to the host while that runner may still be on it; the
   run's manifest is copied into the run directory at the start, so --resume continues it once mex answers
   (a lock that could not be released is then removed as stale).
5. `--sequential`: direct-local on mex with nohup and its own sampler (`--sample-docker`), polled through
   state.json (an ssh answering neither alive nor dead is retried, never read as finished); then
   direct-remote from here with a sampler on mex (`samples/direct-remote-<stamp>.jsonl`, at most 30 min for
   --smoke / 12 h); then the MCP arm (`--mcp-cmd`: repeat 0 of every scenario, then repeat 1) or the copy of
   a finished MCP run directory (`--mcp-run DIR` -> arms/mcp; sequential only).
6. `--mcp-cmd TEMPLATE`: async's runner (`--ids A,B --repeat N --out DIR --manifest M`; appends
   DIR/executions.jsonl and jobs.jsonl, resumable, exit 0 / 1 / 130; 3 without the admin key under
   --check-queue, 2 when the key file is unreadable). Placeholders `{ids}`, `{repeat}` (the repeat index to
   run), `{repeats}` (`{repeat}` + 1, should a runner take a count), `{out}` (RUN_DIR/arms/mcp), `{manifest}`,
   replaced with shell-quoted values; run with `bash -c` in this environment, unchanged (the runner reads its
   admin key from there); its output goes to the terminal only. A full run (not --smoke) takes the default
   template `DEFAULT_MCP_CMD` at the top of run_suite.sh, async's runner of branch feat/suite-mcp-runner:

       <search-service>/.venv/bin/python /Users/raetsch/tmp/mgss/wt-suite/scripts/bench/suite_mcp_runner.py \
         --manifest {manifest} --ids {ids} --repeat {repeat} --out {out} --check-queue \
         --queue-ignore sra-logan-chunks --checks-dir <this directory>

   so a full run is `METAGRAPH_ADMIN_KEY_FILE=/Users/raetsch/tmp/mgss/secrets/staging_admin_key
   ./run_suite.sh` (the key file is passed through, never logged; without it the script refuses to start,
   since the runner checks the service queue). `SUITE_MCP_CMD` replaces the default ('' for none), `--no-mcp`
   runs without the MCP arm, `--mcp-cmd default` asks for it in a smoke run, `--resume` keeps the run's
   template. The runner judges with checks.py of `--checks-dir` (it has no stub fallback).
7. `--max-requests N`: the direct arms' requests of this invocation in total, the closing capabilities GETs
   included (one kept back per direct arm in paired mode; sequential: direct-remote gets what direct-local
   left); the release GET is one more. `--smoke` defaults it to 40 and the manifest to
   manifest/suite_smoke_v4.json; without `--manifest` a full run takes manifest/suite_v4.json (v1-v3 stay for
   reproducing earlier runs).
   The MCP runner paces its own requests.
8. Samples (review M5): the mex samplers write ~/suite-dev/runs/<id>/samples/; every invocation pulls the
   whole directory on exit (also after Ctrl-C) and at its start, merge_run deduplicates.
9. Ctrl-C: the step in flight stops after its request (`pkill -INT` of this run's runner on mex; a local
   step gets the terminal's SIGINT), the orchestrator waits up to 120 s for the mex runner (a long request may
   take 30 min: then it says so, pulls what was recorded and leaves the rest to --resume), pulls, releases the
   lock and prints the --resume command (exit 130).
10. merge_run.py, compare.py, report.py (report.pdf unless --no-pdf). The last stdout line is the run
    directory.

### Comparator and report
- `merge_run.py RUN_DIR`: writes the top-level run.json / executions.jsonl / samples.jsonl / jobs.jsonl /
  manifest.json from `arms/*` and the orchestrator's `samples/*.jsonl` (raw paths rewritten to
  `arms/<arm>/raw/...`; for a key recorded twice the later line wins; lines of another arm in an arm's file
  are dropped with a warning; samples deduplicated by time and container). Re-runnable.
- Outcomes (lane S3): merge_run gives every execution judged by checks.py (not by the stub) the outcome of
  `checks.scenario_outcome` over its checks when that is stricter than its runner's (fail > finding > inconclusive
  > pass; error, not_implemented_mcp and service_by_design are kept), the runner's under `outcome_runner`, counted
  in run.json `outcomes_rederived`: the async MCP runner applies its own rule, which does not tell property checks
  from the generic ones.
- Host cache (review H1): merge_run gives every merged execution line `host_cache`, `since_prev_s` and
  `host_cache_after` from the order in which the run's executions of each scenario reached the host, over all
  arms (direct-local's times shifted by the mex clock offset the orchestrator measured): `cold` = the
  scenario's first execution on the host in this run, `warm` = another execution of it ended at most 60 s
  before this one started, `aged` = longer ago (cache state unknown), `not_reached` = an execution the host
  never saw (an MCP line without evidence of a backend call, i.e. neither `server.elapsed_ms` nor
  `mcp_phases.backend_ms` is a number: a refusal at submit, whatever its outcome, or a request the MCP cannot
  express; a direct request whose connection was refused), `unknown` = no start time. The MCP outcome never
  decides it (review N2): a refusal at submit may read `fail`, and an `error` or a relayed refusal may follow
  a backend call. `cache` (the arm's own
  repeat 0 / later) stays as §8.3 has it; everything timed keys on `host_cache`.
- `compare.py RUN_DIR [--by-design PATH]`: merges first when `arms/` exists, then writes compare.json:
  per scenario and arm the outcome over the repeats (the common one, else error > fail > mixed) and whether
  the repeats' digests agree; direct-local vs direct-remote as `identical` (same digest, no check or outcome
  differing other than by timing: lane S7v), `same_checks`
  (other answer, same check outcomes, at least one check compared), `not_comparable` (other answer, no check
  evaluated on both) or `disagree`, with `checks_compared`; the MCP class of DESIGN §1.3; overheads, cold vs
  warm, timing and resource summaries per area, arm and host cache.
- MCP classes, in this order: the manifest's `mcp.not_implemented` -> `not_implemented_mcp`; per MCP execution as
  async's runner matches it (since lane S5b): the scenario's `mcp.by_design` case(s) in async's table (below), else
  an answer-triggered case of the table (below); the MCP runner's own `not_implemented_mcp` / `service_by_design`
  outcome where no table is read (lane S5c, B5b-05: with a table a service_by_design mark is judged by its checks,
  and a not_implemented_mcp mark on a scenario the manifest maps to a tool, that neither mcp.not_implemented nor a
  named not_implemented_mcp case states, is `disagree`, side mcp: lane S7v); else, on the checks both arms evaluated (a check `not_exposed` on either side
  is not compared; `checks_compared` / `checks_not_exposed` recorded): `disagree` when a check differs or one
  side alone errored, `not_comparable` when no check was compared (review M2: an `agree` can no longer be
  vacuous), `agree` otherwise. No MCP line: `missing`.
- `service_by_design.json`: async's table (131 cases), copied byte-identical from the search-service work
  (`git -C /Users/raetsch/tmp/mgss/wt-suite2 show 31a2e41f:scripts/bench/vendor_suite/service_by_design.json`,
  159,522 bytes, sha256 6b55300b...; `service_commit` inside 22d034ba; lane S5: the triage's corrections, the
  splits by value and four answer-triggered cases; lane S5b: `explains_checks` on every answered / clamped case,
  `answer_match.direct_mismatch` on the four answer-triggered ones, validated at load; lane S5e: async's e427d503
  edit, traverse.duplicate_seeds explains only `([1-9][0-9]*) results for [0-9]+ seeds \(\1 distinct\)`, and
  31a2e41f's case 131 traverse.advanced_path_disabled with `answer_match.mcp_outcome` not_implemented_mcp). A
  case's `agreement` `not_implemented_mcp` -> not implemented yet via MCP; `service_only` (no host counterpart) ->
  service-side by design when the MCP executions match `service_answer`, else `disagree` side mcp; `by_design` ->
  service-side by design when the pair holds: the direct answer matches `direct_expectation` (http, code, text
  with digits ignored, each where stated; the code may also be a pattern slot's) and the MCP executions match
  `service_answer` as async's runner decides each MCP execution (lane S5b, B-02: `refused_at_submit` /
  `relayed`: error.code (alternatives written `a | b`), the text containing `text_contains` (`a | b` = either)
  and the field of the tool's error object matching `field`, each a requirement whether or not the line states
  it; `answered` / `clamped`: no error, and verified only by a line that records the service's answer, a 2xx
  status or a backend call: review N11; since lane S5b an answered / clamped case is a note, the checks are
  compared: "Lane S5b" below). A pair that does not hold is `disagree` with the side (direct / mcp / table for an
  unknown case id). The MCP executions are matched one by one as the runner matches them (lane S5b). `verified` is true only when every side the case has was checked and held; the
  report counts verified and unverified by-design pairs. The MCP arm's execution lines should carry `error`
  {code, text, field}.
- Answer-triggered cases (`answer_match` {routes, direct, mcp}, each side {http, code (or `a | b`),
  text_contains}, lane S5 adds http lists, the mcp side's field and the direct side's slot_code: "Comparator and
  re-evaluation" below; so far pattern.deadline_ends_task: the host's 503 deadline the service relays as a failed
  search, 502 search_failed with the text pattern_deadline): they apply whatever the scenario names. An MCP
  execution on one of `routes` matching the `mcp` side is service-side by design (`by` answer), verified when a
  reference direct execution of the scenario matches the `direct` side, unverified otherwise (the host answered
  within the budget on the direct arm while the MCP call did not: a timing difference under load); the
  scenario's other MCP executions are compared as usual (`other_repeats`; a disagreement there is the class,
  with `answer_case`). The MCP runner (async's vendored runner since 8ccf4243) and reevaluate.py give such an
  execution the outcome `service_by_design`.
- Overheads (DESIGN §5, review H1/M1), on pairs of host-warm executions only: front = direct-remote minus
  direct-local (each arm's first warm execution of the scenario): the wall difference, the client-side
  difference ((wall - server elapsed) remote minus the same locally) and the server difference (near 0 when
  both are warm). Wall times exclude the connection set-up (`client.connect_ms`, reported apart per arm; for
  lines of the older runner that timed from before the connect, compare subtracts connect_ms). A scenario
  whose first arm ran only its cold execution has no warm pair and is counted as such; the cold executions
  are reported apart (per arm: wall and server). MCP: only executions that reached the host (a backend call
  on record, as above); the MCP warm execution (repeat 1 preferred) minus each direct arm's first warm one, with and without
  `mcp_phases.queue_ms`, and `backend_ms` minus each direct arm; the phases of the warm execution; the MCP
  cold execution's wall time apart. Cold vs warm: an arm's host-cold execution and its next host-warm repeat.
- Relations (review of lane S1): the manifest's `relations` [{id, type, a, b, title, args}] between two
  scenarios; compare.py evaluates each per arm on the normalised answers of the first repeat of each there
  (from the raw files; both runners keep the normalised answer under `answer`), with the outcome pass, fail,
  inconclusive (the pair cannot test the claim) or not_exposed (the MCP arm lacks a field). The semantics
  (`same_answer`, `prefix_equals`, `zero_equals_absent`, `subset`) are lane S1b's README; its
  `checks.evaluate_relation` is used whenever checks.py is next to compare.py. Without it, `same_answer` is
  evaluated with the same rule and the other types read `not_evaluated`; an outcome outside the vocabulary
  reads as a failure. compare.json `relations` (per arm outcome and detail) and `summary.relations` (counts
  per arm, which implementation ran); the report has its own table, and a failed relation makes the verdict
  FAIL.
- Resources: each execution joins the samples from its start - 1 s to its end + 2 s (docker averages CPU over
  its own window), else the nearest sample within 5 s; per area and arm the peak container CPU and memory,
  the mean sra_all CPU and host load1, and the peak `work.memory_bytes` of the answers.
- `report.py RUN_DIR [--compare EARLIER_RUN_DIR] [--no-pdf]`: report.md (ASCII-safe pipe tables: common
  signs spelled out, other non-ASCII as '?', sequences over 40 bp cut, markdown specials escaped), then
  `pandoc report.md -o report.pdf --pdf-engine=xelatex`. Verdict: INCOMPLETE when a direct arm did not
  finish (state other than done) or the orchestrator's loop was cut, FAIL on any failed / errored execution,
  any direct / MCP disagreement or failed relation, FINDINGS on any finding, INCONCLUSIVE when no execution
  passed and some were inconclusive (review N5), EMPTY when nothing was recorded, UNVERIFIED when the checks
  stub ran (never PASS: review L2), else PASS; inconclusive executions are counted after it. The
  scripts commit keeps its `+dirty` marker (L1); a failing check is shown with the expectation at the same
  position among the scenario's expectations of that name (L3); the resource columns are captioned as the
  container's figures during the run, not the request's own (L5). The timing tables give the server columns
  their own n (only answers with the host's timing count: N9); `--compare` puts the host-warm wall p50 of the
  two runs side by side (like with like: N12). The load caveat appears when the sra_all CPU mean exceeds
  400 % or the load1 mean exceeds 48 (half of mex's 96 cores). compare.json is computed when missing.

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
- `build_manifest.py`: writes `manifest/suite_v3.json` (the full catalogue) and `manifest/suite_smoke_v3.json`
  (the scenarios of its SMOKE_IDS, tagged `smoke-v3`; the `smoke` tag marks the scenarios earlier smoke
  validations exercised; `suite_v1.json` / `suite_smoke_v1.json` are lane S1's first manifests and
  `suite_v2.json` / `suite_smoke_v2.json` lane S1b's, kept), one scenario per line; `content_sha256` = sha256 of the canonical JSON of
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
    python3 build_manifest.py --no-fetch                # both v3 manifests (seed 20261013; 2-6 min from the cache)
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
  `machine-dependent` (time budgets: 200 or 503, stop or not; for /pattern and /resolve the checks hold this for
  every scenario that expects a 200, tagged or not: lanes S5 and S5b, below), `expected-miss`, `retention-dependent` (an attempt id
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
- Outcomes (lane S1b, the integrator's decision 1; the verification's H3 and M7; lane S3). A check answers
  `pass`, `fail`, `not_exposed` (the MCP arm lacks the field), `inconclusive` (the check could not test its
  property on this answer: an HTTP status the scenario allows that is not a 200, such as a 503 `deadline`; a
  list that a legitimate withheld reason or cut hides the source in, in contains_source_context, source_label,
  source_placement, path_spells, path_record_verified, and the absence checks source_absent, path_absent and
  predicate_selected; a walk a legitimate stated resource limit stopped before its property: "Stated reasons"
  below) or `finding` (budget_respected). Each result also carries `property`: every check is a PROPERTY check
  (the claim the scenario was built to test) except the generic ones, `status`, `invariants`,
  `counts_consistent` and `budget_respected` (`checks.GENERIC_CHECKS`, `checks.is_property`).
  `checks.scenario_outcome(results)` gives the execution's outcome by the precedence fail > finding >
  inconclusive > pass (an outcome outside the vocabulary is a fail), and `inconclusive` when there is no check,
  every check is not exposed, or every PROPERTY check could not test anything (not exposed) although the generic
  checks passed (the verification's H3: an MCP P-exact-0001 without counts is no pass); a scenario without a
  property check (a refusal: its status, code and text are the claim) is judged by its generic checks. Both
  direct runners use it; merge_run applies it to the MCP arm's lines (above). `deterministic` is left out of a
  scenario's first repeat (nothing to compare yet), so that repeat's outcome rests on the other checks.
- `budget_respected` (lane S1b, the integrator's decision 2; lane S3 after the verification of lane S1b, X1:
  SPEC-labeled-traversal-core §6.8 and walker.cpp exactly; no fixed allowance anywhere, since the capabilities
  state no bound on one uninterruptible piece, `deadline_check.max_uninterruptible_ms` null): every scenario that
  names a time budget and allows a 200 carries it, its args stating the rule. /pattern and /resolve: a 200 whose
  `timing.elapsed_ms` exceeds the effective `limits.time_budget_ms` broke the SPEC (past the budget they answer
  503; the deployed /resolve sets `timing.elapsed_ms` before its last deadline check); an allowed 503 `deadline`
  whose client wall time (the runner's `client`) exceeds the budget by more than max(1 s, the budget) is a finding
  too (the deadline is read where the work and the delivery poll it, so the refusal comes at its next reading; the
  wall also holds the network and any wait for a request thread, hence the margin): the smoke run's
  R-budget-0002, a 503 8.9 s after a 750 ms budget, is one, its 8,000 ms twin 0.8 s late is not. /traverse, per
  walked seed: the deadline is the budget on the seed's clock (walker.cpp `time_exceeded`: elapsed >= budget),
  never read at depth 0 (`checkpoint`: depth > 0), and the depth-0 reads are unpaced. So: a time stop before
  the budget (the smallest of `elapsed_ms`, `resource_stop.used` and budget + `stop_after_deadline_ms`, which
  walker.cpp clamps at 0) is a finding; a budget that passed in the seed phase stops the walk at the first depth-1
  checkpoint, so a path longer than 1 bp is a finding (paths' lengths, not complete_to_bp: structural ends make an
  arm complete); a time stop at a checked depth (cap_trigger.at_bp >= 2) must come within one read (or chunk) and
  the head piece after it after max(budget, seed phase), a stop at depth 1 (at_bp <= 1) within the depth-0
  level's pieces of each requested arm (6: up to 3 unpaced level-fetch calls for the root's at most 4 keys, a
  lookahead chunk, two head pieces), each piece at most the seed's `timing.deadline.longest_piece` and at most the
  time since the seed phase (every such piece lies after it); a walk with no time stop may run past max(budget,
  seed phase) only by the pieces before a reading that would stop a head below the radius (1 read + 1 head when a
  path passed 1 bp, i.e. the depth-1 reading came before the deadline, else the depth-0 pieces) plus the radius
  level and the finalisation (2). A breach is the check outcome `finding`; a stop late within its pieces or a seed
  phase longer than the budget is a performance fact stated in the detail. The smoke run's lambda walk (34,201 ms
  under 15,000: seed phase 33,803 ms with one cold 20,556 ms coordinate-tuple read; stopped at depth 1, 398 ms
  after the seed phase) passes, also with its longest piece taken as the 397 ms depth-0 read; the verifier's
  broken walks (a time stop at 2 s of 15 s, a walk to the radius at 60 s without a stop, 5 s late at a checked
  depth with 300 ms pieces, 50 bp walked after a seed phase past the budget) are findings
  (tests/test_checks.py BudgetChecks).
- `deterministic` is inconclusive when this answer or an earlier repeat was limited by a time budget
  (SPEC-labeled-traversal-core §6.8: determinism provided no time budget tripped; SPEC-pattern-search §7.9:
  `determinism: time_limited`); `previous` may give each earlier repeat as a digest or as
  `{"digest", "time_limited"}` (run_direct.py passes the latter since lane S3; the async MCP runner passes digests,
  so on that arm an earlier time-limited repeat is seen only through this answer's own time stop). A time
  limitation that cannot have happened ("Stated reasons") explains no difference: fail.
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
  record's label (a stop before the first k-mer is inconclusive; R-budget runs at 3,000 and 8,000 ms); the hard
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
- Coverage (M8/M5 as verified; toward the design's ~800: 679 scenarios). Not built (lane S3 completed the list from
  the verification of lane S1b): `strategy.run_format`, `min_block_kmers`, `branching.on_reconverge`, the frontier
  knobs, Q-motif's `decided_by: tested_contexts`, the quorum knobs other than `min_live_labels`, HTTP HEAD, more
  T-attempts than the 2 scenarios, R-multi duplicate labels beyond the one R-order, partial-mode planted
  mismatch variants; and: P-unmasked and M-unmasked (the
  host is masked: manifest `not_applicable`), the refuted items (left-arm repeats exist, A^31 cannot be withheld),
  a body above max_request_body_mb (256 MB would load the shared host), P-modes "at the count" (the count is not
  known before the run), a deadline stop of an attempt (needs a run-time not_after_ms), max_live_paths (needs a
  branch at a known depth under a header label; the repeat windows give one, left for the next increment).
- Not changed (judged a risk to break the MCP arm): the walk preset's memory. The direct bodies send 512 MB, as the
  live service capabilities of 2026-10-13 stated; the service's code default is 256 MB and the MCP refuses an
  override above its preset row, so naming 512 in the traverse_fetch overrides would refuse every MCP walk if the
  row were 256. The run should pin the service's capabilities (run.json) and say which applied.

## Lane S3: the integrated suite (2026-10-13)

The suite is developed on the integrated tree from here on (base f9cdca10: lanes S1, S1b, S2, S2b, S2c). What lane
S3 changed, beside the sections above:

- Offline tests choose their checks: `run_suite.sh --checks auto|checks|stub` (above) and `compare(...,
  relation_evaluator=None)`; tests/test_orchestrator.py `test_real_checks` runs the orchestrator with checks.py.
- `budget_respected` per SPEC §6.8 and walker.cpp, and the 503 long past its budget (above).
- Property and generic checks, `scenario_outcome` in every runner (above).
- Stated reasons (the verification's X2), below.
- Residuals: N-hard-0022 asserts the repeat (`blocked_event`: each arm of the (TTAGGG)^10 seed re-enters the seed
  at 0 bp, `rejoined_seed`, SPEC §6.6, since the period closes a cycle inside the seed; the graphlet's `E b`
  records are read for the MCP arm); R-junction asserts each label's runs exactly (`resolve_runs` `runs`: a label
  with an extra run, A covering B's k-mers, fails); the N-hard all_or_count scenarios of L >= 30 expect `withheld`
  null (at most 16 contexts, never above max_contexts 100).

### Stated reasons
A check that cannot test its property because the answer states a hiding (a withheld list, a cut), a stop or a
time limitation accepts that statement only once the answer shows it can have happened; otherwise the check FAILS
(and so do the route's invariants), it is never inconclusive. Not judged where the answer lacks the quantity (a
missing cap, `work`, `timing`, an arm's `limitations` or `counters`): then the statement is taken as stated. A
graphlet result (`detail: graphlet`, and the MCP's export) carries its arms' cap_trigger, limitations and counters
as JSON beside the MGT text, and is judged like a full one (lane S4, G1 of the verification of lane S3; the MGT's
A record gives the cap trigger and counters where a result has no JSON arm).

- /pattern, `withheld.reason` (SPEC-pattern-search §7.5): `count_above_threshold` the primary count (contexts; for
  L > k paths, supported_paths) exact, or bounds with its upper bound, above max_contexts (max_paths); a lower
  bound does not show it (that is `threshold_crossed`); `threshold_crossed` stop_at_threshold set;
  `discovery_budget` a max_steps stop with the steps of all patterns at max_steps; `deadline` a time stop or
  `determinism: time_limited` (and that time reason possible, below); `paths_later_increment` L > k without
  long_search paths / supported_paths; `anchors_above_threshold` the anchors above max_anchors;
  `annotation_budget` a max_annotation_work or max_memory stop or a refused row; `anchor_labels_truncated` a
  non-empty `anchors_truncated`; `output_budget` a max_memory stop (or labels_status output_budget);
  `predicate_above_threshold` the raw count above max_predicate_contexts; `selected_above_threshold` the
  selected count above max_contexts (max_paths); `predicate_budget` a max_predicate_work or max_memory stop or a
  refused row; any other reason fails.
- /pattern, `cut.reason`: `max_contexts` the count (the selected one with a predicate) not exact within
  max_contexts, unless stopped at the threshold; `max_steps` a max_steps stop that reached max_steps; `time` a time
  stop or time limitation; `max_memory` a max_memory stop, a refused row or labels_status output_budget;
  `max_paths` the paths not exact within max_paths, unless stopped at it; `max_annotation_work`,
  `max_predicate_work`, `max_anchors` that stop; `max_predicate_contexts` the raw count not exact within it, unless
  stopped at it; any other reason fails.
- /pattern, `stop.reason` (§7.6): `time` and `determinism: time_limited` (also a low-complexity diagnostic the work
  time cut) need `timing.elapsed_ms` at or after a lower bound of the work time, time_budget_ms -
  finalize_reserve_ms - 2 E, E the §7.6 formula at the default (slowest) delivery rates 10 / 50 MB/s over the
  results in the answer plus every withheld result at a generous 2 KiB (64 KiB with labels or a predicate);
  `max_steps` the steps of all patterns at max_steps (every fixture's stop shows them equal); `max_contexts`,
  `max_anchors`, `max_paths`, `max_predicate_contexts` stop_at_threshold; `max_annotation_work` a search that reads
  annotation (labels all, supported_paths, a predicate); `max_predicate_work` a predicate; `max_memory` taken as
  stated; any other reason fails. A stop of the graph search (discovery, mask_scan, extraction, extension,
  selection) with neither a cut nor a withheld list is an invariant failure (§7.5: only these shorten a list), and
  a list with neither holds every result: a source missing from it fails (retrieval_complete false alone means
  the labels are incomplete).
- /resolve, `stop` (SPEC-labeled-traversal-core §4.5): reason `time` only, with `timing.elapsed_ms` at or after the
  budget less the reserve.
- /traverse, an arm's `cap_trigger` (§6.7, §7.0): the arm's `walk_domain` limitation of the trigger's knob
  (traverse.cpp cap_knob: max_steps, max_live_paths / beam_pruned, max_paths, max_output_bp, time_budget_ms; for
  resource_limit the knob of the seed's resource_stop, which must exist) with the observed demand past its limit
  (time: at or after it); max_output_bp the arm's `counters.output_bp`, max_steps the seed's steps within the at
  most 4 successors a head follows of the cap (walker.cpp cap_check); max_paths the arm's paths and live paths
  within that of max_paths; time_budget a stop at or after the budget (as budget_respected reads it).
- /traverse, a path's resource end (lane S4, G2 of the verification of lane S3; SPEC-labeled-traversal-core §7.0): a
  path whose `path_reason` or label end reasons name a resource reason lies on an arm the cap left short of its
  radius (`truncated`, or `pruned` by a beam) whose `cap_trigger` names a resource reason, and the path's reason is
  that trigger's or a later cap's that the arm's walk_domain limitations name ("the cap in cap_trigger, then any
  other cap that ended walks after it"); otherwise walk_in_record and walk_end fail, and so do the invariants.
- /traverse, the attempt's duration bound (SPEC §6.8; lane S3's rule of 6f945db4, narrowed by lane S4 after the
  verification of lane S3, C1-C4): a seed result whose `error` is, in full, the host's statement that the attempt's
  bound stopped it is a stated stop of THAT seed: its own checks (those with a `seed` argument) are inconclusive,
  the other seeds keep theirs, and the checks of the whole answer (invariants, budget_respected) judge the others
  and pass over it. The statements, matched in full in the deployed host's wording and the current source's (which
  adds the in_ram load term), never as a substring (the MCP tool's hint "the attempt reached the duration bound the
  server enforces for it; a retry needs a new attempt_id" is none): AttemptAtBound (traverse_attempts.cpp: the
  bound reached while the response was built or written, a 503 `deadline` that the MCP service relays as the seed's
  error), on the MCP arm only, since a direct arm receives the 503 itself; and the attempt's walk-until before the
  seed began or in its seed phase (traverse.cpp not_started_seed_to_json, walker.cpp seed_external_stop: a 200 with
  the seed failed), on any arm. The stated bound must have been reached where the answer states the attempt's
  elapsed time: AttemptAtBound's bound against `usage.elapsed_ms` and the client's wall time (`client`; the MCP
  runner passes the service's backend time), a walk-until stop at or after `usage.bound.walk_until_ms`, with a
  resource_stop (where stated) of `attempt_deadline` and the seed failed; else the seed's checks and the
  invariants fail. Both make the answer time-limited for `deterministic`. tests/test_smoke_v3_fixes.py rebuilds the
  texts from the C++ string literals, so a change of the host's wording fails a test.
- The checks that accept a stated reason: contains_source_context, source_label, source_placement, path_spells,
  path_record_verified, source_absent, path_absent, predicate_selected (a hiding), count_relation (a budget stop),
  low_complexity_note (a budget stop, a time-limited diagnostic), deterministic (a time limitation),
  resolve_runs (a resolve stop), walk_in_record and walk_end (an arm's cap trigger), and the invariants of each
  route. The verifier's examples are tests (tests/test_checks.py StatedReasons): X-peptide10-0009 (exact 2,810 of
  max_contexts 10,000) and N-realrep-0007 (exact 2) withheld as count_above_threshold, N-hard-0005 (exact 2 of
  100), a 300 ms discovery time stop of a 60,000 ms budget, a 500 bp walk "truncated by max_output_bp" under
  500,000 and one "by max_steps" under 1,000,000, a walk time stop 2 s into 15 s: each fails; the same statements
  where they can have happened stay inconclusive.

### Manifests v3
`suite_v3.json`: 679 scenarios, 719 executions per arm, 48 not implemented via MCP (content_sha256
75af21ee3122a98055f1a7fc029d77bb19030fdaada099a54faa2b6d2f5eea46), the 22 relations of v2 (relations_sha256
cfba509f111c1e1ed727035450ee983b69991e02d4933c40bb1a8da59f7c822f, unchanged); `suite_smoke_v3.json`: 11 scenarios,
11 executions, 1 relation (content_sha256 64df233674404a5c5fbd2c5dd2501d45e70cd12959db06eec118743f571bc86c,
relations_sha256 8d5fa169494d23960080e74682c99d4a6697ad8598f4f977b156546c0dafd1b4). Built twice from the cache
(`--no-fetch`, seed 20261013): identical digests, every top-level key but created_utc equal. Against v2: the same
679 ids in the same order; every budget_respected rule text is new (603 scenarios); beyond it only R-junction-0001..4
(exact runs), N-hard-0002/0005/0008/0012/0014 (withheld null), N-hard-0022 (blocked_event) and D-R-junction-0003
(the copy of R-junction-0003) changed; tags as v2's plus `smoke-v3`. The smoke selection (SMOKE_IDS): R-junction-0001,
R-budget-0001 (a 503 past its budget), T-end-0001, T-budget-0001 (a 1 ms walk budget: a stop at depth 1),
P-exact-0001 (the MCP's property exposure), X-peptide10-0009, M-refusal-0012 and M-twin-0001 (their same_answer
relation), N-hard-0005, N-hard-0022, N-realrep-0007: in paired mode with all three arms at most 9 + 4 x 11 + 4 = 57
host requests (the warm repeats of indices 0 and 10 included), fewer where the service refuses at submit
(R-budget-0001's 10 kb query, T-end-0001's budget above the walk preset's row).

## Lane S4: the last check gaps, and the offline re-evaluation of a recorded run (2026-10-13)

From the Opus verification of lane S3 (VERIFY.md; base 6f945db4):
- G1: `_arm_stop_problem` no longer exempts a graphlet result: the MCP's export and a direct `detail: graphlet` answer
  carry their arms' cap_trigger, limitations and counters as JSON beside the MGT text, and a stated cap stop there is
  judged as on a full answer ("Stated reasons"); `graphlet_arms` also reads the A record's cap trigger and counters
  for a result without JSON arms. tests/test_checks.py StatedReasons.test_mcp_shaped_stops builds the MCP shape from a
  full answer (`mcp_shaped`): the "max_output_bp at 500 of 500,000" walk fails on the MCP and on a direct arm, with
  and without its walk_domain limitation; the same stop under a max_output_bp of 500 stays inconclusive.
- G2: a path's resource end needs an arm the cap stopped, with a matching cap trigger ("Stated reasons");
  StatedReasons.test_resource_end_needs_a_cap_trigger (the verifier's E7, a later cap, walk_end).
- The attempt-bound rule of 6f945db4: per seed, the host's exact wording (both variants), the MCP arm for
  AttemptAtBound, the stated bound reached ("Stated reasons"; tests/test_smoke_v3_fixes.py AttemptBound, HostWording).
  The stopped seed's own checks are inconclusive; the invariants judge the other seeds and pass over it (6f945db4
  made them inconclusive too).
- Test isolation: run_suite.sh takes its runner pattern from SUITE_RUNNER_PATTERN (default unchanged,
  `[s]cripts/run_direct\.py`); tests/test_orchestrator.py narrows it to runners naming the test's own fake host (and
  a stand-in whose command line ends in the script), so a real run's ssh step on this machine, whose command line
  holds scripts/run_direct.py, no longer reads as the test's runner (the verification's two failures).

### Offline re-evaluation of a recorded run
`reevaluate.py RUN_DIR [--checks DIR] [--by-design PATH]` judges a recorded run's answers again with the current checks (checks.py of
DIR, default this directory), without any request; `merge_run.py RUN_DIR --reevaluated`, `compare.py RUN_DIR
--reevaluated` and `report.py RUN_DIR --reevaluated` then build the re-evaluated set, and `run_suite.sh --reevaluate
RUN_DIR [--reeval-checks DIR] [--compare EARLIER] [--by-design PATH] [--no-pdf]` runs the four (offline: before the
lock, the ssh connection or the release GET).
- Input: every arm's executions.jsonl (the last line per scenario, arm and repeat, as merge_run keeps it) and the
  raw answer of each line. Both runners keep the NORMALISED answer in the raw file under `answer` (run_direct.py;
  async's suite_mcp_runner.py writes its normaliser's answer there, beside the tool calls), so the MCP arm is
  re-evaluated like the direct ones and no normaliser runs again; on both smoke runs of lane S3 the digest of every
  raw answer equals the line's recorded answer_digest on all three arms (the answer the runner judged).
- What the checks are passed, as the runners pass it: `previous` = the scenario's earlier repeats on the arm in
  repeat order as {digest, time_limited}, both by the re-evaluating checks from the earlier raw answers (the
  recorded digest and time_limited where a raw answer is missing), so the MCP arm gets the L4 form too, which
  async's runner did not pass; `client` = {"wall_ms"}: a direct arm's recorded wall time (run_direct.py), the MCP
  arm's `mcp_phases.backend_ms` (async's runner since its 0d127f4a: the service-measured backend call, since the
  MCP wall holds the queue and the polling; none without it).
- Outcome: `checks.scenario_outcome` of the new checks, except where the runner decided it without them, which is
  kept: `error` (no answer, transport, decode, a busy deployment; but an `error` whose recorded checks are the
  checks module's own exception is judged again), `not_implemented_mcp`, `service_by_design` (async's table; also
  an MCP answer an answer-triggered case of the table matches, `--by-design PATH`, default this directory's), and
  async's host-state rule (lanes S5 and S5b: the named-case and answer-case replay, `resolve_timeout` and
  `traversal_resolve_deadline` with their plausibility rules, below; an MCP answer refused with
  `uninterruptible_overrun_exceeds_allowance` the scenario does
  not expect is inconclusive). A new exception in checks.evaluate is `error` on a direct arm and `fail` on the MCP
  arm, as the runners record it.
- Output, beside the recorded files, which all stay: `arms/<arm>/executions.reevaluated.jsonl` (the recorded line
  with checks, outcome, answer_digest and time_limited replaced, `answer_digest_recorded` where the digest changed,
  and `reevaluation` {checks_sha256 (16 hex), outcome_recorded, checks_changed, kept}), `reevaluated.json` (the
  set's mark: the checks' directory, git commit with `+dirty` when checks.py or suite_common.py differ from it,
  sha256 of both files; what each arm recorded with; per arm counts; every execution whose outcome or checks
  changed, with each changed check's before, after and detail), and from the three scripts
  `executions.reevaluated.jsonl`, `run.reevaluated.json` (`reevaluation`), `compare.reevaluated.json` (the
  relations judged by the re-evaluation's checks) and `report.reevaluated.md` / `.pdf` (title "(re-evaluated)", a
  "Re-evaluation" section after the verdict with the counts and every changed execution; never UNVERIFIED for a run
  recorded with the stub, since checks.py judged it). samples.jsonl, jobs.jsonl and manifest.json are shared (written
  by the re-evaluated merge only where missing). `--compare EARLIER` compares with the earlier run's re-evaluated set
  where it has one.
- The checks are loaded as their own module, with their directory's suite_common.py, so a re-evaluation with an
  older or newer checks.py (`--checks DIR`, e.g. an export of a commit) digests and judges as that version.
- tests/test_reevaluate.py: a synthetic recorded run over three arms (the attempt bound reached and not, a stale
  check, an evaluate exception, the kept outcomes, a 503 judged by the wall and by the backend time, a time-limited
  cold repeat), a checks.py of another directory that records what it is passed, every recorded file byte-identical
  after the re-evaluated merge, compare and report, and run_suite.sh --reevaluate with nothing reachable.

## Lane S5: the suite fixes of the verified triage (2026-10-14)

From the triage of the full run's first part (triage-1/TRIAGE.md, groups G01-G25) as its three verifications amend it
(VERIFY-T, VERIFY-RP, VERIFY-MN; where a verification amends the triage, the verification is followed).

### Checks (judged offline by re-evaluation)
- Walks through merges (G13, as VERIFY-T amends it). A merged walk (`branching.on_reconverge: merge`, the walk
  preset's default) displays each merge through its first parent (SPEC-labeled-traversal-core §7.1: the parent
  carried by the most labels; the arrival order in `tree` / `full` detail under a memory budget), so after a split
  that several records (a column's alleles, a derived set) or repeat copies carry, a displayed path may be another
  record's allele while the record's own route passes the merge through another parent. `walk_in_record` and
  `walk_end` now search the root-to-leaf routes of the arm's segment DAG through `segments[].parents`, taking ANY
  parent at a merge (`_record_routes`); a route follows the record when its whole spelled flank is the truth's
  continuation (of any copy with `any_copy`), and its end reasons are its leaf's (`_leaf_ends`: the leaf's paths[]
  entry and the label_end events on the leaf at its end). `labels_via_parent` is not followed: VERIFY-T showed that
  a label's kept entry may come through the other allele's parent (T-repeat-0001, T-walk-0036), where only the
  any-parent rule holds. The search runs forward from the roots and is pruned at the first base that departs from
  every copy, so it visits a segment at most once per copy (the children of a split differ in their first base);
  it stops after `_MAX_ROUTE_STATES` (100,000) visits, and a check that could not finish it is inconclusive, never
  failed. `walk_end` asks for a route that follows the record to its end (where the truth's flank reaches the end;
  otherwise the leaves' lengths as before), its leaf's end reasons in `reasons_in`.
- Quorum ends (G16). `paths_end` accepts a reason that the label_end event on the path's leaf at its end names: the
  deployed engine has no `below_min_labels` end reason, a quorum or split-limit stop is `branch` with the qualifier
  text in that event (SPEC §13, I3 note 2). `graphlet_arms` derives the label_end events from the R records as the
  SPEC's §7.5.6 does (a code end, on the run's segment at `to_bp`, the reason the qualifier's text through
  `_END_QUALIFIERS`, e.g. `Rb` below_min_labels, or the code's enum; `Lw` and `m` give none); on every mini-index
  fixture they equal the full result's events.
- The deadline (G21, as VERIFY-MN amends it: required, not optional; lane S5b G22 and lane S5c A2-01, A2-07). A
  POST /pattern 503 with code `deadline` is an answer SPEC-pattern-search §6 / §7.6 allow for any search whose answer
  could not be written within its time_budget_ms, so on a loaded host every searching /pattern request may get it;
  a POST /resolve one is SPEC-labeled-traversal-core §4.5's answer past the budget, for a request that names
  bounds.time_budget_ms only (the /resolve deadline is opt-in: without the field there is none, so such a 503
  fails). `checks.deadline_allowed` (the old name `pattern_deadline_allowed` an alias), on both routes: where the
  scenario's status expectation asks for a 200 and does not itself name the deadline (an http 503, the code
  `deadline` or a text naming it), such a 503 makes `status` inconclusive and joins evaluate's allowed list, so the
  body's checks read inconclusive and `budget_respected` judges it by the client's wall time (`_budget_503`, against
  min(the requested budget, the route's cap the builder states as cap_ms): before it, /resolve less its 250 ms
  reserve, a fail; long past it a finding). A scenario that expects a refusal (no 200), names the deadline itself,
  or another route or code is judged as it states.

### Comparator and re-evaluation
- The table: async's 789fb5f0 (130 cases, sha256 2fdea36b...), byte-identical, with the triage's corrections (G09,
  G10, G11, G19, G20), the splits by value (the old id keeps the value suite_v3 sends; new:
  pattern.predicate_options_noop_without_predicate, resolve.time_budget_above_reserve_below_minimum,
  resolve.discover_limit_above_100_host_answers) and the new cases pattern.caps_zero and traverse.empty_seeds;
  answer_match on pattern.deadline_ends_task, mismatch.above_kind_cap, pattern.caps_zero and traverse.empty_seeds.
  Async's runner (wt-suite2 scripts/bench/suite_mcp_runner.py, ServiceByDesign) is the reference reading of it.
- Field alternatives and placeholders (G08, G07 as VERIFY-RP has it). A case's `field` may list alternatives
  `a | b`, split as its code is (outside `<...>`), and a `<...>` placeholder matches anything, dots included, as
  async's runner (`ServiceByDesign._like`) reads it: `predicate<path of the first fault>` matches
  `predicate.at_least.n`, `bounds.<knob>` matches `bounds.a.b` too (lane S3's test said otherwise; the runner decides
  the execution's outcome with its own reading, so the comparator follows it).
- The host side of a case (G09, G10; async's additions). `direct_expectation` may state `http` as a status or a list
  (a budget that gives 200 or 503), `code` (a top-level or a pattern slot's, `a | b`), `slot_code` (a pattern
  slot's only), `text` (whole, digits ignored), `text_prefix` (as is, no digit folding, as the table's vocabulary
  states: a label-agnostic text such as "Label not found in the annotation: ") and `text_contains`; `text_source` and `behaviour` are descriptive. A
  case whose host side states nothing (G11) has its direct side judged by the scenario's own `status` check on the
  reference direct execution: pass verified, fail a disagreement (side direct), anything else unverified.
- Keys a newer table may state. A side's key the comparator neither matches nor knows as descriptive is never
  ignored when it is stated (present, not null): that side is unverified (None) with a detail naming the key; an
  answer_match whose mcp side states one never triggers. tests/test_compare_report.py checks that async's table
  states none.
- Host-state executions (G23 as VERIFY-RP amends it). `resolve_timeout` (the service's resolve client gave up before
  the host's overrun answer under load) joins `uninterruptible_overrun_exceeds_allowance` in
  `compare.MCP_HOST_STATE_CODES`, which reevaluate.py applies (inconclusive unless the scenario expects the code);
  when a named case is judged, an MCP execution that is inconclusive by such a code is skipped (it reads neither for
  nor against the case), and with none left the pair is not comparable (lane S5b, B-09), not a disagreement.
- Answer-triggered cases (B6; the order is async's runner's since lane S5b: per MCP execution the named cases
  first, an answer-triggered case only for an execution no named case matches). An answer_match case takes the MCP
  executions its mcp side matches (service side by design, `by: answer`); the named cases judge theirs against the
  direct executions no answer case explains (a direct answer of its direct side, such as the host's own 503
  deadline, says nothing about the rest's pair); a disagreement of the rest is the class, with `answer_case`. answer_match sides: direct {http (a status or a list), code (top-level or slot), slot_code (the
  codes of body.patterns[].error), text_contains}, mcp {http (one status: lane S5b, async's runner compares it by
  equality), code, text_contains, field (as above; an MCP answer that names no field does not match)}, as async's
  runner reads them (`answer_matches`); since lane S5b both apply the named cases first.
- The runner's matcher replayed (VERIFY-RP G08 point 3; the coordinator's second addendum asks the re-evaluation to
  reproduce async's own replay). reevaluate.py applies the current table to every recorded MCP answer as async's
  runner (at e334715b since lane S5b) does (`compare.named_cases`, `compare.service_answers`: the named cases with
  their service_answer, the field read from the recorded tool error as the runner's `raw_error` reads it, a case
  field against an answer naming none not matching), then the answer-triggered cases: a refused / relayed match is
  service_by_design with that case; a runner's service_by_design mark the current table does not reproduce is
  judged by the checks (lane S5b, B-06) and says so. So an execution the runner judged `fail` against a case its
  matcher could not read (G08's alternation) becomes service_by_design offline.
- A run recorded before a normaliser fix (G06: a pattern refused in its slot) is re-normalised first with async's
  `scripts/bench/suite_mcp_renormalise.py RUN --out COPY` (read-only use of its worktree), COPY put in place of a
  copy's arms/mcp; reevaluate.py itself never re-runs a normaliser.

### The v4 builder (built later by the integrator)
`build_manifest.py` now builds manifest version 4 (`manifest/suite_v4.json`, `suite_smoke_v4.json`, the smoke tag
`smoke-v4`; it refuses to write over a frozen v1-v3 manifest) with `--capabilities DIR` (default
truth/staging_capabilities/): the frozen v4 is built once the host's capabilities are refreshed after its deploy
(the snapshot of 2026-10-13 still states a protein substitution cap of 3; the deployed one is 1). The manifest's
`index` states the snapshot directory and its `mismatches_max`. What changed, by triage group:
- G03: planted mismatch windows take the shortest length whose chance estimate, the uniform ball model scaled by
  the window's composition factor where above 1 (`areas._planted_window`), stays below 1,000 (a tenth of
  max_contexts: the safety factor); a window drawn after a rejected one comes from its own random stream, so the
  area's other windows stay. `truth.chance_contexts_estimate` states it. In the provisional build M-planted-0001
  (mge, 20 -> 25), M-planted-0007 (ser, m 3, 25 -> 31) and M-planted-0013 (sso, m 3, 25 -> 31) change.
- G04: the m = 2 decapeptides only under a snapshot protein cap of at least 2, the cap slot (max_mismatches cap + 1)
  only under a cap below 3; the ids M-protein-0008..0010 are kept, a scenario not built is listed in `skipped`. The
  until-next-deploy and needs-next-deploy tags are gone (the snapshot decides).
- G05: `pattern_by_design` names mismatch.above_kind_cap above a kind's cap of the snapshot
  (`manifest_lib.set_host_capabilities`, `HOST_MISMATCH_MAX`), for patterns of at most k bases.
- G07: Q-more-0003 (at_least 3 of 2 labels) expects the refusal, 400 invalid_request naming
  request.predicate.at_least.n; `pattern_by_design` names pattern.predicate_form_other for a form fault of
  SPEC-pattern §19.3 other than a non-string name (`_predicate_form_fault`).
- G10: the builders' absent label is `manifest_lib.ABSENT_LABEL`; a resolve naming it names resolve.unknown_label.
- G12: R-edge-0019 sends an object cut before its end (no JSON reader accepts it).
- G14: the header kind of the derived-label pairs (T-labels-0004, -0006) sends labels.max_seed_labels 10,000, the
  server's and the service's maximum: the offline truth cannot count a window's header carriers (it knows the
  genomes' own records, not which RefSeq records hold a window), and the eco window's set held 6,181;
  labels.max_seed_labels is an override traverse_fetch sends.
- G15: T-coords use the trace preset (`manifest_lib.TRACE_PRESET`, the service's: on_reconverge keep).
- G17: T-attempts-0002 asks for usage on the 200 only (the duplicate's 409 has none: inconclusive).
- G18: T-budget-0003's max_work_units is 300,000 (about twice the 151,383 units its seed's validation took, below
  the run's 1,000 bp header walks); the first run should confirm a partial walk.
- G19, G02: `traverse_by_design` names traverse.empty_seeds for `seeds: []`, `pattern_by_design` pattern.caps_zero
  for max_contexts, max_anchors or max_paths 0.
- G21: the M-planted, M-protein and N-nearmiss /pattern builders allow [200, 503] and are tagged machine-dependent
  (gen_D then picks other M scenarios for determinism).
- Splits by value (async's 789fb5f0): pattern.predicate_options_without_predicate for predicate_scope "motif"
  without a predicate, pattern.predicate_options_noop_without_predicate for predicate_strands or "context" alone;
  resolve.time_budget_below_minimum at <= 250 ms, resolve.time_budget_above_reserve_below_minimum at 251..499;
  resolve.discover_limit_above_100 for max_labels below 1, resolve.discover_limit_above_100_host_answers above 100.
A provisional v4 with the current snapshot (lane S5's directory, not manifest/) differs from suite_v3 in 99
scenario ids (3 added, 4 dropped, 92 changed; every one mapped to its group), and with the run's capabilities
(protein cap 1) only in M-protein-0008/0009 (not built) and M-protein-0010 (built). suite_v3.json and
suite_smoke_v3.json are unchanged.

### Validation on the recorded run (offline)
A copy of the full run's first part (542 scenarios, 2,237 executions; the run of 2026-10-13 on b7214d29 under the
Logan load), its MCP arm re-normalised with async's suite_mcp_renormalise.py (864 of 1,045 executions rebuilt,
exactly the 12 of G06 changed), re-evaluated with this branch (`run_suite.sh --reevaluate`), against the recorded
re-evaluation reeval-B (checks of 4ff45cc5): 90 execution outcomes and 35 MCP classes change, every one in a verified
group. Fails: direct-local 32 -> 8, direct-remote 28 -> 9 (G21: 29 direct 503 deadlines now inconclusive; G13, G16),
MCP 49 -> 2 (G06 12, the table's cases replayed G02/G05/G08/G19 18, G23 1, G13/G16 16); MCP disagreements 42 -> 15,
unverified by-design pairs 8 -> 0. What is left waits for the v4 manifest (G04, G07, G10, G12, G14, G15, G17) or a
re-run of the MCP arm on the fixed service (G01), or is the load (G21, G22, G24).

## Lane S5b: the verified review of lane S5 and the triage's second part (2026-10-14)

From the review of lane S5 (review-s5/FINDINGS-A.md A-01..A-12, FINDINGS-B.md B-01..B-11) as VERIFY-REVIEW amends it,
the second part of the triage (triage-2 G27, G28, G22) as VERIFY-2 amends it, and the coordinator's addendum of the
semantics agreed with async (explains_checks, direct_mismatch, the table's validation, the host states).

### Checks
- A premature deadline (A-01). A 503 `deadline` arrives only when the answer cannot be written by the budget, and the
  deadline starts at the parsed body, inside the client's interval: `_budget_503` (budget_respected on an allowed
  503) fails one whose client wall time is below the budget less 1 ms (/resolve: below the budget less its fixed
  finalisation reserve of 250 ms, SPEC-labeled-traversal-core §4.5). A /pattern scenario where the deadline can be
  allowed carries budget_respected (v3: 426), as do the 85 budgeted /resolve scenarios that ask for a 200.
- The /resolve deadline (G22 as VERIFY-2 item 3 has it). `checks.deadline_allowed` (formerly
  `pattern_deadline_allowed`, kept as an alias) allows a POST /resolve 503 `deadline` as it allows the /pattern one
  (SPEC-core §4.5: past the budget the answer is a 503, never partial), on the direct arms (the MCP relays it as
  `traversal_resolve_deadline`: a host state, below).
- Determinism after a deadline (A-03 / G27, VERIFY-2 item 1, variant B). `time_limited` states a deadline answer
  first: the host's 503 `deadline` (/pattern, /resolve), the service's relay of a /pattern deadline (502
  `search_failed` whose text starts `pattern_deadline`) and of a /resolve one (503 `traversal_resolve_deadline`).
  `deterministic` compares a repeat only with the earlier repeats no time budget limited: none left is inconclusive
  ("no earlier repeat without a time limitation"), and a time-limited repeat is inconclusive unless its stated
  limitation cannot have happened. A repeat after an earlier time-limited one is still compared with the other
  complete ones (variant A, which only counted the deadline, masked a real difference between r1 and r2). The
  MCP runner passes plain digests (no time_limited): async passes {digest, time_limited} as run_direct.py does.
- Routes carry the record's label (A-04). Where the segments state their labels (`segments[].labels` of `detail:
  full`; a graphlet's G records' entry sets, which `graphlet_arms` now derives as SPEC §7.5.3 states, the library's
  derive.py being the reference: G.entry against parent[0]'s end set or a merge's union, `*` the seed labels at the
  root and a merge's partition union, G.end with the runs ending before the segment's last node removed and the
  switch targets added, annotate's last P set), the route search passes only segments carrying the record's own
  label (its header or its taxid's column by the label dictionary, a graphlet's L records), and a segment whose list
  the cap cut (`labels_truncated`) passes. Without the record in the dictionary (G14's cut) or without segment
  labels the search stays sequence-only, and the check's detail says so (`[sequence only: ...]`). On the mini-index
  fixtures and on the 714 recorded MCP graphlet segments of the full run whose direct-local twin has the same shape,
  the derived sets equal the full result's labels.
- The route search is memoised (A-05) on (segment, the seed copies still followed, whether the stated event was seen):
  every route to a segment has the same length and the copies it still follows fix its flank, so each state is
  searched once; `_MAX_ROUTE_STATES` (100,000 states) is a safeguard. A route counted "beyond the truth flank" now
  follows the flank to its end before it runs past it.
- A walk blocked by its own seed or edges (G28, VERIFY-2 item 2: the strict rule; the triage's loose rule accepted a
  relabelled reason, a continuation off the record and a fabricated block). An arm with no route that follows the
  record over its whole length still follows it when (1) a route spells the record for exactly its first d bases,
  where (2) d is the first depth at which the record's own continuation is inadmissible by SPEC-core §6.6,
  computed from the truth's flank and the request's seed (`_first_block`: rejoined_seed, a seed k-mer, either
  orientation in a canonical/primary regime; else edge_reuse_rc / edge_reuse, the step's (k+1)-mer or its reverse
  complement already used on the walk), (3) one of its segments states the blocked event at at_bp d with the
  record's base and that reason, and (4) past d it spells the record from a collapsed repeat's exit (an occurrence
  of the junction's k-1 bases in the flank followed by another base) for every base the flank holds. A block before
  the truth's min_bp is a truth gap (inconclusive), so is a continuation past the flank; past the record's end with
  the record as the only label it fails. A following route that ends at the computed block before the truth's
  min_bp, with that reason, is the same truth gap. The regime is the answer's capabilities or the graphlet's H
  record (basic on this index).


### The MCP verdict, the comparator and the merge (async's runner at e334715b is the reference)
- The table: async's e334715b (`git -C /Users/raetsch/tmp/mgss/wt-suite2 show
  e334715b:scripts/bench/vendor_suite/service_by_design.json`, 154,008 bytes, sha256 0fcf7e7d..., 130 cases):
  `explains_checks` on every answered / clamped case (a check name, or {check, detail_matches}: only a failure whose
  detail the regex matches completely), `answer_match.direct_mismatch` (disagree | timing) on the four answer-
  triggered cases, `answer_match.wall_rule` (budget_minus_reserve) on pattern.deadline_ends_task, text_contains
  alternatives on pattern.predicate_form_other. `compare.load_by_design` refuses what async's runner refuses at load
  (`table_problems`: an answer_match without routes, an unknown key of it or of a side, an MCP side with an http
  list, a slot_code or no code, an unknown direct_mismatch or wall_rule, an explains_checks that is not a list of
  valid entries or sits on a refused / relayed case): `TableError`, compare.py and reevaluate.py exit 4.
- The verdict (reevaluate.outcome_of replays the runner's `judge`): the cases the scenario names whose
  service_answer the MCP answer matches (in the order named; several answered ones apply together), else the first
  answer-triggered case whose MCP side matches and whose wall rule holds or cannot be decided. A refused / relayed /
  answer-triggered match is service_by_design, `finding` when a check is a finding (a by-design mark never erases
  one). An answered / clamped match never decides before the checks (B-01): F = the failing checks (fail, or an
  outcome outside the vocabulary), E = those of F with outcome fail that an entry of the union of the matched cases'
  explains_checks explains, U = F - E: fail when U is non-empty, else r = the checks' rule without E: finding, or
  service_by_design when E is non-empty, else r. The line records `service_by_design_case` (the first matched case
  explaining an item of E, else the first matched) and `explained_checks` (E) whenever a case matched; the outcome a
  runner's mark the current table does not reproduce is the checks' (B-06, said in `reevaluation.note`); '' = no
  table (then a runner's mark is kept).
- Walls and plausibility (suite_common, each threshold one function): the service-side wall only, never the
  client's (`service_wall`: mcp_phases.backend_ms; for a pattern search the deadline-ended task's one failed HTTP try
  from the diagnostic log's held-slot line; for a traversal job the last settled attempt's span, then the job's
  finished_at - started_at, from the raw record's phase_sources). pattern.deadline_ends_task (R1) applies when that
  wall >= the request's time_budget_ms (the MCP arguments', else the request's) - the finalisation reserve its
  text states (else 250 ms); a known wall below holds the case back (matching continues, the checks judge).
  resolve_timeout (R2): wall >= max(budget + reserve + 10,000 ms, the stated "within N s") - 1,000 ms;
  traversal_resolve_deadline (G22): wall >= budget (the tool error's, else the MCP arguments') - reserve; the R20
  lockout has no wall rule. Rule (c): an unknown threshold or wall makes the verdict inconclusive with
  verdict_detail "no threshold known: plausibility undecidable" / "no service-side timing: plausibility
  undecidable" (the case, or the host-state code, recorded); a finding wins on the wall-rule path, not on the
  host-state path (as async's runner; lane S5b's addendum-11 divergence is dropped since lane S5c, B5b-03). Lane
  S5c's additions (the client's wall as an upper bound, a discovery's deadline, the stated resolve_timeout phase):
  below.
- `inconclusive_reason` on every inconclusive line (rule (c)'s text, the host-state code, or "checks"); compare.json
  sums it per arm (`summary.inconclusive_by_reason`, from the line or, for a line recorded without it, from the
  raw record), the report's verdict section lists it, so a run whose deadlines are undecidable does not read clean.
- merge_run leaves a line the host-state rule or rule (c) made inconclusive (its inconclusive_reason; else the raw
  record's wall), a by-design decision (refused / relayed / answer-triggered, a finding included) and the checks an
  answered case explains out of its stricter outcome.
- compare.mcp_agreement matches each MCP execution as the runner does: an answer-triggered match judges its direct
  side (`answer_class`; direct_mismatch B-03: "disagree": any direct HTTP answer outside the direct side disagrees,
  "timing": a direct 200 or 503 deadline is unverified, another answer disagrees; a transport error or no direct
  answer is unverified; without the key unverified), a refused / relayed named match its pair (`by_design_class`,
  against the direct executions no answer-triggered case explains), and the executions an answered case matched or
  none did are compared check by check (`check_agreement`), the explained checks apart: their difference alone is
  service_by_design (`by: explains_checks`, verified by the case's direct side), any other difference disagrees, and
  the case is a note (`by_design_notes`). Host-state and undecidable executions read neither way; all of them is
  not_comparable (B-09). A disagreement anywhere is the class, except when the scenario names a not_implemented_mcp
  case (resolve.not_expressible.bounds_max_query_bp: the difference is the expressibility gap, as the first named
  case not disagreeing decided before); a by-design class is verified only when every other part agrees (B-04).
  Known keys are checked before an unknown one makes a side unverified (B-05). `--reevaluated` reads the table the
  re-evaluation replayed and warns on another (B-10).
- The direct arms' class `timing` (VERIFY-2 item 5, all its conditions): the outcomes of a repeat differ only
  because one side's own answer states a time limitation (a stop by time, a 503 deadline), that side reads
  inconclusive, the other pass or inconclusive, no check both evaluated differs, and the stopped answer is
  consistent with the complete one (/pattern: since lane S5c the ranges the population counts admit overlap and
  the released results lie among a retrieval_complete side's, below; /traverse: every stopped path a prefix of a
  complete one; a 503: no content); or (G27) a repeat inconclusive only
  by its determinism check after a time-limited earlier repeat, with the other side's answer the same. Counted
  apart in compare.json (`timing`, `timing_scenarios`) and the report.

### The v4 builder
- A-02: the M-planted, M-protein and N-nearmiss /pattern scenarios ask for a 200 again (status [200, 503] passed any
  503: a proxy page, `loading`); they keep the machine-dependent tag, and checks.deadline_allowed admits exactly the
  SPEC's 503 `deadline` in place of the 200.
- A-06: `_planted_window` tries, after the first candidate length, only longer ones, shortest first, so the last
  one, taken whatever its estimate, is the longest (no current window changes).
- A-07: `_within_kind_caps` leaves out of areas M and N every scenario whose max_mismatches lies above the
  snapshot's cap of one of its patterns' kinds (a pattern of at most k bases), listed in `skipped` with its id
  consumed, unless it expects the slot refusal (mismatches_above_cap) or a refusal of the request; D picks come from
  what is left.
- A-08: build_manifest.py refuses as an output any suite_v*.json / suite_smoke_v*.json of manifest/ that exists or
  that git tracks, by whatever path (realpath, samefile: a symbolic link to manifest/ included), and --areas with
  the default outputs (a partial build needs files of its own).
- A-09: the manifest's index records `capabilities_sha256`, the sha256 of each capabilities document, beside the
  snapshot's directory. The snapshot of 6f945db4 is committed in truth/staging_capabilities/ (protein cap 1).
- A-10: T-attempts-0002 checks the duplicate's 409 shape (answer_fields with `when_http: 409`: attempt.state and
  attempt.usage present, no response-level usage; on a 200 that check passes, and the usage checks judge it) and
  states `run_unique: ["attempt_id"]`: run_direct.py appends "-<run id>" (lane S5c: "-<run id>-<arm>", A2-06) to
  the request's attempt_id (the run directory's name, which both direct arms share, so within a run the first arm
  took the 200 path and the second met the retained id's 409; the next run, a new id), cut to the id pattern's 128
  characters, and records the request it sent; the manifest keeps the base value, so its digest stays the same for every run. The MCP arm does
  not run the scenario (not implemented via MCP).
- A-12: (a) gen_D adds a complete m = 2 pick: the forward planted window at max_mismatches 2 with the lowest chance
  estimate (D-M-planted-0005, a 31-mer at level 1); (b) a scenario for each split case of async's table without
  one: R-edge-0023 (discover.max_labels 150: resolve.discover_limit_above_100_host_answers), R-edge-0024
  (time_budget_ms 400: resolve.time_budget_above_reserve_below_minimum), Q-motif-0005 (predicate_strands without a
  predicate: pattern.predicate_options_noop_without_predicate), tagged split-case (never a D pick).
A provisional v4 with the snapshot of 6f945db4 (lane S5b's directory, not manifest/) differs from suite_v3 in 96
scenario ids (7 added, 5 dropped, 84 changed; every one mapped to its group) and from lane S5's build with the
run's capabilities (the same protein cap 1) in the 4 added scenarios, T-attempts-0002 and the 64 statuses of A-02.
suite_v3.json and suite_smoke_v3.json are unchanged; manifest/suite_v4.json is built once async's table is final.

### Validation on the recorded full run (offline)
A copy of the full run (679 scenarios, 2,900 executions; 2026-10-13 on b7214d29 under the Logan load), its MCP arm
re-normalised with async's suite_mcp_renormalise.py (1,128 of 1,330 rebuilt, 16 changed: G06), re-evaluated with
lane S5's scripts (5d5dd2fa) and with this branch (daeb497d), against the recorded re-evaluation reeval-final-B
(lane S4's checks, MCP arm not re-normalised). Fails: direct-local 40 -> 11 -> 8, direct-remote 31 -> 10 -> 9, MCP
66 -> 7 -> 2 (Q-more-0003 x2: v3 expects a 200, G07); MCP service_by_design 320 -> 344 -> 134, pass 834 -> 864 ->
1,071, inconclusive 46 -> 51 -> 59. S5 -> S5b, every change in a group: B-01 (210 MCP lines an answered case had
marked: 204 pass, 6 inconclusive: C-caps-0001/0003, N-usecase-0005), G27 (D-M-planted-0015 r1/r2 on direct-local and
MCP), G28 (N-realrep-0006 on all three arms), rule (c) (N-realrep-0002 MCP r0: no service-side timing); MCP classes:
104 by-design pairs -> agree (B-01), R-edge-0008 and N-realrep-0002 disagree -> agree (B-01 P4, rule (c)),
N-nearmiss-0028 by design verified -> unverified (B-04); direct classes: 11 disagreements -> timing (VERIFY-2 item 5:
9 deadlines or time stops, T-walk-0004's prefix, D-M-planted-0015 with G27's extension), T-attempts-0002 left (G17).
The MCP arm reproduces async's replay of its runner at e334715b (inconclusive by reason: checks 55, resolve_timeout
3, no service-side timing 1); its fails differ by the 16 re-normalised lines (G06) and the 4 lines G27 / G28 fix.

## Lane S5c: the verified review of lane S5b (2026-10-14)

From the review of lane S5b (review-s5b/FINDINGS-A.md A2-01..A2-07, FINDINGS-B.md B5b-01..B5b-10) as its
VERIFY-REVIEW amends them, and the coordinator's addendum 1 to lane S5c (async's exact rules; the runner's commit
follows, and each rule is one function here, aligned when it arrives: they landed as 62ee4214, 2a8acbaf and
f94b079d, and lane S5d below aligns the mirror with them).

### Checks
- The /resolve deadline is opt-in (A2-01, SPEC-labeled-traversal-core §4.5 l.224; capabilities
  resolve.time_budget.default null). `deadline_allowed` admits a POST /resolve 503 `deadline` only for a request
  naming bounds.time_budget_ms; on a budget-less one (the discoveries R-discover-0001..0006, R-edge-0023,
  N-negative-0014/0016, N-usecase-0010/0012) the status fails, saying why, and the other checks read "no answer".
- The effective budget (A2-03). `_budget_503` judges an allowed 503 against min(the requested budget, the
  budget_respected argument `cap_ms`), the route's cap of the snapshot the builder states where a request asks more
  (/resolve resolve.time_budget.max_time_ms 60,000, /pattern caps.time_budget_ms 600,000): R-edge-0021 (70,000 ms)
  answered 503 at 60,100 ms is inconclusive, no longer a fail; a manifest without cap_ms is judged as before.
- A block whose repeat exit lies past the truth's flank (A2-04). When the computed block at d has no exit inside the
  flank, a route that spells the record's first d bases, states the block (at_bp d, the record's base, the computed
  reason) and continues is inconclusive ("... which is not in the truth's flank: not verifiable"); one that ends at
  the block passes; without the stated block the arm still fails; another copy that decides wins.
- Duplicate seeds (B5b-07). The invariants check's count message states the distinct seeds when the request holds
  identical ones (the service's identity, planner.py normalize_seeds: the sequence whitespace-free and upper-cased,
  a FASTA header dropped; the labels trimmed and de-duplicated in order; no seed_id): exactly
  `"%d results for %d seeds (%d distinct)"`, e.g. "1 results for 2 seeds (1 distinct)"; without duplicates
  `"%d results for %d seeds"` as before. Async's regex can then explain exactly a merge:
  `([1-9][0-9]*) results for [0-9]+ seeds \(\1 distinct\)`. Until async's table has it, the table's current entry
  no longer matches T-edge-0008's detail (its MCP line reads fail).

### The MCP verdict (suite_common, reevaluate)
- A discovery's relayed deadline (A2-01, addendum 1 item 1): `discovery_without_budget`: a traverse_resolve whose
  arguments name `discover` and no time_budget_ms gets no deadline from the service, so its
  traversal_resolve_deadline is implausible (the checks judge: status fails), verdict_detail "a /resolve deadline on
  a discovery without a budget: no deadline was sent". An explicit-labels resolve without a caller budget stays
  "budget unknown": the stated budget's threshold, else undecidable.
- The client's wall as an upper bound (B5b-01, addendum 1 item 2): where the service-side wall is unknown,
  `client_wall_bound` = t1 of the call that returned the answer (the last timed submit, terminal get_status or
  traverse_results call of the raw record) - t0 of the execution's FIRST submit call (a busy retry's first call may
  have started the work: V-03). client wall + 1,000 ms below the threshold: provably premature, R3 implausible (the
  checks judge), R1 the case does not apply; else undecidable as before. Skipped when the runner called the MCP with
  an API key, and the detail says so: PROVISIONAL field `mcp_auth` "api_key" on the line or the raw record (the
  e334715b runner is anonymous and records nothing); async names the field with its commit. (Superseded by lane
  S5d: rule 2 applies only to a recorded `mcp_auth` "anonymous", with the runner's call choice and texts.)
- resolve_timeout (addendum 1 item 4): a refusal that states `timeout_phase` (read | write | deadline) and `bound_s`
  (the tool error's, else the answer error's) is a host state at a wall >= bound_s x 1,000 - 1,000 ms, any phase,
  the phase named in the detail; a legacy refusal (no stated phase): min(10,000 ms, max(A, B)) - 1,000 ms (R2's A
  and B, capped by the 10 s write phase, V-04). The recorded ~30.4 s lines stay host states.
- B5b-03: an undecidable host state is inconclusive whatever the checks say (lane S5b's addendum-11 divergence
  dropped: the replay equals the runner again). B5b-10: with no table, a recorded `finding` whose line names a
  service_by_design_case is kept, as the runner's mark is.

### Comparator and merge
- The timing class's content rule (B5b-02, VERIFY-REVIEW): each population count is the range of true values it
  admits (`count_range`: exact v [v, v], at_least v [v, inf), bounds [lower, upper], unknown anything; a relation
  outside exact / at_least / bounds / unknown, or a bound that is no number, a contradiction) and a stopped count is
  consistent with the complete side's at the same path when the ranges overlap; `tested` (work done, §19.7) and
  `estimate` are no population counts; the stopped entry's released results must lie among the complete entry's
  when that states retrieval_complete. at_most is gone.
- direct_mismatch over every direct repeat (B5b-04): `answer_class` classifies each repeat of BOTH direct arms; a
  disagreeing one wins over a match.
- Runner marks (B5b-05): with no table an unmatched execution the runner marked service_by_design is an unverified
  by-design part of its own, the others compared by their checks (a disagreement wins); with a table such a mark is
  judged by its checks.
- Explained check items (B5b-06): E per ITEM as the runner keys it: `compare.explained_items` (a name explains every
  failing item of the check, {check, detail_matches} the items it matches); `check_agreement` sets the explained
  items apart with the reference's items of the same check and occurrence; the re-evaluation records
  `explained_items` (indices into checks) beside `explained_checks`, and merge_run leaves out those items only.
- Host state vs undecidable wall rule (B5b-09): `host_state_line` holds only for a host-state code; a relayed
  /pattern deadline whose wall rule is undecidable is not_comparable "by undecidable plausibility";
  `inconclusive_reason` replays the answer-triggered wall rule for a line recorded without the field.
- The planned 409 (A2-06): `direct_agreement` reads a repeat whose arms answered 200 and 409 on a run_unique scenario
  that allows 409, neither side failing, as a planned duplicate (`planned_duplicates`), not a disagreement.

### The v4 builder
- A2-02: R-budget-0001/0002, P-caps-0004, P-heavy-0001/0002, N-usecase-0008 and R-twin-0002..0005 ask for a 200 (no
  builder states [200, 503] any more); R-twin-0002..0005 are tagged machine-dependent. All ten carry a budget, so
  the SPEC's 503 deadline stays allowed.
- A2-03: `set_host_capabilities` takes the routes' time-budget caps from the snapshot; budget_respected states
  `cap_ms` where a request asks more (R-edge-0021).
- A2-05: with --areas, an --out / --smoke-out whose realpath is a suite_v*.json / suite_smoke_v*.json of manifest/ is
  refused, by any spelling (`manifest_path`).
- A2-06: run_direct appends "-<run id>-<arm>" to a run_unique field, so each arm runs its own attempt;
  T-attempts-0002 repeats twice (r1 meets its arm's retained 409 and exercises the 409 shape check).

### Validation on the recorded full run (offline)
A new copy of the full run (lane-s5c/reeval-S5c: full-run-1-final-copy with lane S5b's re-normalised MCP arm, inputs
byte-identical to lane S5b's reeval-S5b-2) re-evaluated with this branch: every outcome and class equals lane
S5b's but T-edge-0008's two MCP lines (B5b-07): the invariants detail now reads "1 results for 2 seeds (1
distinct)", which the table's current traverse.duplicate_seeds regex does not match, so they read fail and the
pair disagree (service_by_design 134 -> 132, MCP fail 2 -> 4, MCP disagree 15 -> 16) until async's table carries
the tightened entry; with it (a scratch copy of the table) both lines read service_by_design again, E = the one
invariants item. Line fields: explained_items on the 344 MCP lines a case matched (B5b-06); the 3 recorded
resolve_timeout lines (~30.4 s) stay host states, their verdict_detail naming the legacy rule (item 4); nothing
else changed (no recorded budget-less /resolve deadline, no bounds count in a timing pair, the one undecidable line,
N-realrep-0002 r0, stays undecidable: its client wall 22,522 ms is not below 19,750 ms). The direct T-attempts-0002
disagreement stays in this v3 run (its 409 fails v3's answer_fields); A2-06 applies to runs of v4.
A provisional v4 (snapshot 6f945db4, lane-s5c/build/, not manifest/) has content_sha256 3a7417d6... (relations
cfba509f..., unchanged; smoke 8e944aeb...), deterministic under PYTHONHASHSEED 1 and 777; against lane S5b's
provisional (695acdd7...) it changes 12 scenarios and nothing else: the 10 statuses of A2-02 (R-twin-0002..0005
also tagged machine-dependent), R-edge-0021's cap_ms (A2-03) and T-attempts-0002's repeat 2 (A2-06), plus the
conventions' run_unique text.

## Lane S5d: the verified review of lane S5c, aligned with async's landed runner (2026-10-14)

From review-s5c/FINDINGS.md (S5c-01..05 and three notes), each verified first (lane-s5d/verify/VERIFY.md: every
finding CONFIRMED on its failing input against this branch at b17832df and async's runner at f94b079d). The
reference is async's runner on feat/suite-mcp-runner-vendored at f94b079d (scripts/bench/suite_mcp_runner.py):
e334715b plus 62ee4214 (rule 1, a /resolve deadline on an unbudgeted discovery), 2a8acbaf (rule 2, an anonymous
client wall bounds an unknown service-side wall) and f94b079d (rule 4, a resolve_timeout judged against the bound
it states). Its table is unchanged (0fcf7e7d..., byte-identical here).

### What the runner records (and the mirror reads)
- `mcp_auth` (2a8acbaf; every line and the raw record): "anonymous" | "api_key" (suite_mcp_client.SuiteMcp.mcp_auth:
  "api_key" when a tools/call would present a key). `client_wall_ms` (pattern_search / traverse_resolve lines) and
  the raw record's `client_wall` {ms, first_submit_t0, first_submit_call, answer_t1, answer_call, answer_tool}.
- `timeout_terms` (f94b079d; a resolve_timeout line and every raw record): {timeout_phase, bound_s, elapsed_ms,
  source}, null when the error states none (the legacy bound).
- verdict_detail texts: rule 2 "client wall X ms + 1 s slack < threshold Y ms: provably premature" (a wall-ruled
  case prefixed "<case> (wall_rule budget_minus_reserve): "); with an API key "no service-side timing: plausibility
  undecidable; rule 2 skipped: the MCP calls present an API key, ..." (inconclusive_reason without the note);
  rule 4 "the stated bound B s (phase) - 1000 ms" and "min(the write cap 10000 ms, ...) - 1000 ms".

### The mirror (suite_common, reevaluate, compare)
- S5c-01, rule 2 tri-state: `mcp_auth_of` reads the line's `mcp_auth`, else the raw record's; `client_wall_premature`
  applies the bound only for "anonymous", skips it with the runner's note for "api_key", and decides nothing for any
  other value or none (undecidable, as the runner). Every record before 2a8acbaf (the recorded full run) states
  none: `reevaluate.py --mcp-auth anonymous` (run_suite.sh --reeval-mcp-auth) states the arm's auth for such lines,
  recorded in the mark (`mcp_auth`: {assumed, recorded counts, assumed_lines}) and on each line it applied to
  (reevaluation.mcp_auth_assumed, which compare reads too); a recorded mcp_auth always wins.
- S5c-02, rule 4: `timeout_statement` as the runner's: the first of the tool error, the answer error, the job's
  stored error (the last traverse_results resolve_error.detail, any code) that states a phase in read | write |
  deadline AND a positive numeric bound_s in that one object; anything else is legacy.
- S5c-03, the client wall: `client_wall` as the runner's: the execution's FIRST submit call (its t0 unknown: no
  wall), the answer call the last timed traverse_results else the last traverse_resolve (traverse_resolve), the
  FIRST terminal get_status else the last pattern_search (pattern_search), ms to 0.1; the job wall to 0.1 from the
  runner's ISO reading. The runner's recorded client wall (raw `client_wall`, else line `client_wall_ms`) and
  `timeout_terms` are read first and cross-checked against the raw record: a mismatch is a warning (the line's
  reevaluation.warnings, the mark's warnings), the recorded value is used (it is what the runner judged with).
  verdict_detail and inconclusive_reason are composed as the runner's judge composes them; rule 1 reads the MCP
  arguments of any tool, as the runner's unbudgeted_discovery (only traverse_resolve has `discover`).
- S5c-04: run_direct's run_unique id never cuts the arm: a run id too long for the 128 characters is cut and ends in
  "." + 8 hex of its sha256 (two long run ids stay apart); a base that leaves no room is refused (ValueError).
- S5c-05: compare computes an answered case's explained items with the table in hand, as the runner computes E per
  execution; a line that records another list (a re-evaluation with another table) is noted under
  `explained_items_recorded`, not trusted. merge_run keeps reading a line's own explained_items (its own verdict).
- Note 3: a planned 409 split on one repeat now says which arm's state was not the plan (repeat 0: an id retained
  before the first send, e.g. a re-run after --resume; a later repeat: an attempt not retained, e.g. a host restart).
- Note 1 (no change): `deadline_allowed` keys on the request's bounds.time_budget_ms because the snapshot's
  resolve.time_budget.default is null (the deadline is opt-in); a service default would need that default read from
  the snapshot. Note 2 (no change): `_budget_503` on a /pattern request without time_budget_ms stays not testable;
  the v4 scenarios without one expect a 400.

### Equality with the runner and validation on the recorded full run (offline)
- Equality (lane-s5d/equality/): async's judge at f94b079d (a read-only copy, imported) and reevaluate.outcome_of on
  the same input compared on outcome, case, explained checks, explained items, verdict_detail and
  inconclusive_reason: every recorded MCP line of the full-run copy (1,282 lines with an answer), as recorded (no
  mcp_auth) and with the auth stated anonymous: 0 differences, and the re-evaluated lines as written equal the
  lane's verdict; 90,000 mutated executions (70,000 general, 20,000 timing answers mostly without a service-side
  wall; mcp_auth absent / anonymous / api_key / another value on the line, the raw record or as the assumption;
  submits with and without t0, busy retries, several terminal get_status, timed and untimed traverse_results;
  stated, unstated and invalid timeout phases in the tool, answer and job errors; the runner's recorded client wall
  and timeout_terms): 0 differences. The fuzz found one before the fix: rule 1 read only traverse_resolve's
  arguments, the runner any tool's.
- A fresh copy of the full run (lane-s5d/reeval-S5d: full-run-1-final-copy with lane S5b's re-normalised MCP arm,
  inputs byte-identical to lane S5c's) re-evaluated with this branch: every outcome, MCP class and direct class
  equals lane S5c's; the only line-field change is the verdict_detail of the 3 recorded resolve_timeout lines (~30.4
  s, still host states), now in the runner's words ("min(the write cap 10000 ms, ...) - 1000 ms"). The run states
  no mcp_auth (it predates 2a8acbaf), so rule 2 decides nothing; re-evaluated with --mcp-auth anonymous
  (reeval-S5d-anon) nothing changes either: no recorded line has an unknown service-side wall and a client wall
  provably below its threshold (the one undecidable line, N-realrep-0002 r0, has 22,522 ms >= 19,750 ms).
- The provisional v4 (lane-s5d/build/, snapshot 6f945db4) is unchanged: content_sha256 3a7417d6..., relations
  cfba509f..., smoke 8e944aeb..., deterministic under PYTHONHASHSEED 1 and 777, 0 scenarios added, dropped or
  changed against lane S5c's (run_unique's hashing is run_direct's, not the manifest's).

## The frozen v4 (2026-10-14)

`manifest/suite_v4.json` and `manifest/suite_smoke_v4.json`, built with `SUITE_TRUTH_CACHE=<lane S1's record cache>
build_manifest.py --no-fetch` from the capabilities snapshot of release refseq97-6f945db4 (truth/staging_capabilities/,
protein mismatch cap 1; each document's sha256 in the manifest's index): 681 scenarios, content_sha256
3a7417d60e4eb5d0e54433c76139df7ed1221ae35bcccf814b7fca2cf913f113, relations cfba509f111c1e1ed727035450ee983b69991e02d4933c40bb1a8da59f7c822f;
smoke 11 scenarios, content_sha256 8e944aeb792da696a544236156c3ea8c8b435f7a53548689ba520cff7b7e0d6e. It equals the provisional v4 of
lanes S5c/S5d. run_suite.sh takes v4 by default (`--smoke`: suite_smoke_v4.json); build_manifest.py refuses to overwrite either.


## Lane S5e: the mirror of async's final runner 31a2e41f (2026-10-14)

The reference is async's runner on feat/suite-mcp-runner-vendored at 31a2e41f (scripts/bench/suite_mcp_runner.py):
f94b079d plus e427d503 (the suite re-vendored at 3daea97a; the tightened traverse.duplicate_seeds regex), f52846fb
(rule 4's read floor), f386a3c0 (`explained_items` on the line and the raw record), 901893b7
(`answer_match.mcp_outcome`) and 31a2e41f (case 131, traverse.advanced_path_disabled).

- The table: async's 31a2e41f (`git -C /Users/raetsch/tmp/mgss/wt-suite2 show
  31a2e41f:scripts/bench/vendor_suite/service_by_design.json`, 159,522 bytes, sha256 6b55300b..., 131 cases),
  byte-identical here (f386a3c0's, 154,062 bytes, e6330470..., 130 cases, was the lane's first copy, a96c9303).
  Against e334715b's (0fcf7e7d...): traverse.duplicate_seeds' explains_checks entry is {"check": "invariants",
  "detail_matches": "([1-9][0-9]*) results for [0-9]+ seeds \\(\\1 distinct\\)"} (e427d503; fullmatch: as many results
  as distinct seeds, the checks' wording of lane S5c, B5b-07; the vocabulary's example "1 results for 2 seeds (1
  distinct)"); the vocabulary documents answer_match.mcp_outcome; case 131 traverse.advanced_path_disabled
  (31a2e41f): refused_at_submit at stage mcp_local, agreement not_implemented_mcp, direct_expectation http 200,
  answer_match {routes ["POST /traverse"], direct {http 200}, mcp {http 400, code preset_unavailable, text_contains
  "is on the advanced strategy path, which is not enabled on this deployment", field "overrides.<knob path>"},
  direct_mismatch disagree, mcp_outcome not_implemented_mcp}. Async's vendored checks.py and suite_common.py at
  31a2e41f equal this branch's at 3daea97a (the runner keeps its own rules; it does not call the mirror).
- Rule 4's read floor (f52846fb; `suite_common.read_bound_below_floor`, called by host_state_verdict before the
  least wall, as the runner's host_state_plausible calls it): a resolve_timeout whose statement (the runner's
  recorded `timeout_terms`, else timeout_statement) names timeout_phase "read" with a known budget (the MCP
  arguments' time_budget_ms, else the tool error's) must state bound_s x 1000 >= budget + reserve (the tool error's
  finalize_reserve_ms, else 250) + 10,000 - 1000 ms: the service's client_timeout_s gives a budgeted /resolve a read
  timeout of at least budget + reserve + 10 s. Below it the service sized its own read timeout too short: no host
  state whatever the wall (a known, a long or an unknown one; rule 2 is not consulted), the checks judge,
  verdict_detail "service read bound {b} s below its own floor {f} s" (both %g seconds, the runner's `:g`), no
  inconclusive_reason. Write and deadline phases keep bound_s x 1000 - 1000 ms (they fire earlier by design); no
  budget known and the legacy path (no stated phase) are unchanged. compare (host_state_line) and merge_run
  (host_state_inconclusive) read the same rule.
- `explained_items` (f386a3c0): the runner records on a matched case's line, right after explained_checks, and in
  the raw record the indices into the line's `checks` of the results with outcome fail that an entry of the union of
  the matched cases' explains_checks explains (explained_item_indices, in check order: the set its
  explained_outcome's E holds); [] whenever E is empty (a refused / relayed / answer-triggered case, a rule (c)
  undecidable match, an answered case explaining nothing); absent when no case matched. That is the field and the
  semantics lane S5c gave the re-evaluation (B5b-06: `compare.explained_items` over `compare.union_of` of the
  matched cases; the same entry_explains, the same union, the same order), and the re-evaluated line places it
  after explained_checks as the runner does. reevaluate.py replaces the recorded list by the replay's and notes a
  recorded one it does not reproduce (`reevaluation.explained_items_recorded`); compare computes the items with the
  table in hand and notes a differing recorded list (`explained_items_recorded`), never trusting it; merge_run
  leaves out a line's explained_items (the runner's or the re-evaluation's), and only for a line without the field
  (a runner before f386a3c0) every failing item of its explained_checks.
- Where the runner records rule 2's and rule 4's terms (its writer, Runner.execute / main at f386a3c0), and where the
  mirror reads them: `mcp_auth` on every line and in every raw record (the client's SuiteMcp.mcp_auth, "anonymous"
  | "api_key"; null only from a client without one), and in mcp_runs.jsonl / run.json (the latest invocation's at
  the top, each one's under `invocations`): `mcp_auth_of` reads the line's, else the raw record's; the
  re-evaluation's --mcp-auth assumption now applies only to a record with the key in neither (a run before
  2a8acbaf; `mcp_auth_recorded`): a recorded null stays unknown, as the runner judged it (before, the assumption
  overrode it). The mark records run.json's values (`mcp_auth.run_json`), and a line whose recorded value is none of
  them is a warning. `client_wall_ms` on a pattern_search / traverse_resolve line (null when an instant is unknown)
  and the raw record's `client_wall` ({ms, first_submit_t0, first_submit_call, answer_t1, answer_call,
  answer_tool}, null for another tool): `client_wall_ms` reads the raw record's first, else the line's, as before.
  `timeout_terms` ({timeout_phase, bound_s, elapsed_ms, source}) on a resolve_timeout line (null when the error
  states none) and in every raw record (null unless a resolve_timeout states them): `timeout_terms` reads the raw
  record's first, else the line's, as before; rule 4's bound and its read floor read it. explained_checks /
  explained_items / service_by_design_case on a matched line and in the raw record, verdict_detail and
  inconclusive_reason on the line where set (the raw record: null): read from the line, as before.
- `answer_match.mcp_outcome` (901893b7; `compare.match_mark`, `compare.not_implemented_detail`): "service_by_design"
  (the default when absent) or "not_implemented_mcp". The load validation (`compare.table_problems`, as the runner's
  ServiceByDesign.problems) admits it as an answer_match key with one of the two values and refuses another value, the
  key on an answered / clamped case and the key at a case's top level. The FIRST matched case decides (named or
  answer-triggered): with not_implemented_mcp the execution reads not_implemented_mcp (the class a manifest
  mcp.not_implemented gets), the case recorded (service_by_design_case; explained_checks and explained_items []),
  verdict_detail "not implemented via MCP on this deployment (case {case}, answer_match.mcp_outcome
  not_implemented_mcp): the service refused {field} at submit: {text!r}; on a deployment that enables it the request
  is accepted and the case does not apply" (the runner's MCP_NOT_IMPLEMENTED_DETAIL; field the tool error's `field`
  when it states the key, else the answer error's, none "the request"; the wall-rule notes first, joined "; "); rule
  (c) and the host-state rules do not apply; a known wall below a wall rule still keeps that case from matching; a
  check that is a finding makes it finding. reevaluate.py replays a runner's not_implemented_mcp mark that names its
  case with the table (the manifest's mcp.not_implemented, which names none, stays the runner's decision): the case
  reproduces it, a table without it judges the checks (`reevaluation.note`), no table keeps the mark with the
  runner's own case fields (as a kept service_by_design mark and a by-design finding now keep theirs). compare counts
  such executions under not_implemented_mcp (`by` mcp_outcome, the case, the repeats; the pair as its match would
  read it, answer_class or by_design_class, under `pair`), never service_by_design, ahead of any other part.

### Equality with the runner and validation on the recorded full run (offline)
- Equality (lane-s5e/equality/): async's judge at 31a2e41f (a read-only copy, imported; its table byte-identical to
  this branch's) and reevaluate.outcome_of on the same input compared on outcome, case, explained checks, explained
  items (the runner's Verdict.explained_items, cross-checked against its explained_outcome's E), verdict_detail and
  inconclusive_reason, the explained fields compared as the writers record them (present only on a matched line):
  every recorded MCP line of two fresh full-run copies (1,282 lines with an answer), as recorded and with the auth
  stated anonymous: 0 differences, the re-evaluated lines as written equal the lane's verdict and place
  explained_items after explained_checks; 130,000 mutated executions (seeds 31-35 general, 41-42 timing answers, 51-53
  rule 4's read floor: stated read / write / deadline bounds just below, at and above budget + reserve + 10 s - 1 s
  for the arguments' or the tool error's budget and the stated or default reserve, known, long and unknown walls;
  61-63 mcp_outcome: the real table and, for half of them, synthetic gated cases prepended, i.e. an advanced-path
  refusal, a wall-ruled pattern deadline, a resolve_timeout and a /resolve deadline code, a stated default mark and
  a gate without a field, with field variants on the answer and the tool error, named gates, findings and held-back
  walls; besides lane S5d's mutations, a recorded null mcp_auth beside an assumption and failing items of the cases
  an answered case explains): 0 differences (3,555 executions decided by the read floor, 5,065 with explained items,
  3,708 by mcp_outcome not_implemented_mcp). The harness detects each of four planted defects (no read floor, the
  first explained item only, no mcp_outcome, the detail's field from the answer error only).
- A fresh copy of the full v3 run (lane-s5e/reeval-S5e and reeval-S5e-anon with --mcp-auth anonymous: full-run-1-final-
  copy with lane S5b's re-normalised MCP arm, inputs byte-identical to lane S5d's) re-evaluated with this branch and
  the 31a2e41f table: against lane S5d only T-edge-0008's two MCP lines change, fail -> service_by_design
  (traverse.duplicate_seeds, E = invariants item 1, "1 results for 2 seeds (1 distinct)"), the pair disagree ->
  service_by_design verified (MCP fail 4 -> 2, service_by_design 132 -> 134; classes disagree 16 -> 15,
  by_design_verified 59 -> 60); every other outcome, check list, line field and class is unchanged, the
  anonymous-auth copy equals the recorded one, and the re-evaluation with f386a3c0's table equals the one with
  31a2e41f's (the v3 run sends no advanced-path key and no recorded resolve_timeout states a phase).
- v4's T-labels-0004 / -0006 (lane-s5e/synthetic/run-tlabels: the frozen v4 scenarios, the MCP arm refused as
  staging refuses them, 400 preset_unavailable on overrides.labels.max_seed_labels; r0 recorded as a runner before
  901893b7 (fail), r1 as the runner at 31a2e41f records it; the direct arms 200): every MCP line re-evaluates to
  not_implemented_mcp with case traverse.advanced_path_disabled, explained [] / [], the runner's verdict_detail and
  no inconclusive_reason; compare's class not_implemented_mcp (by mcp_outcome, the pair verified by the direct
  200s); both judges agree on all four lines.
- 251 tests (SUITE_TRUTH_CACHE=lane S5's truth cache) pass.

## Lane S6c: a strict deploy gate, the run's plan, the manifest identity and an honest report (2026-10-14)
Lane G3 of the fix phase after the outside review of the suite (GPT, verified by Opus in VERIFY-A / VERIFY-B; the owner
decided on 2026-10-14 that the deploy gate is strict: every inconclusive, missing run or failure blocks it).

### Manifest identity (finding 9: the frozen content hash did not bind the relations)
- `suite_common.manifest_identity(manifest)` is the identity of a test plan: content_sha256 of `scenarios` and
  relations_sha256 of `relations` (of [] when it has none). The relations are judged across the scenarios' answers,
  so the same scenarios with other relations are another plan. `verify_manifest` fails when a stated relations_sha256
  does not match the relations (`manifest_problems` names each stated hash its content does not match); a manifest of
  before relations_sha256 verifies on its scenarios. `suite_common.py manifest-sha FILE --write` writes both hashes.
- run_direct.py records `computed_relations_sha256` (and the stated `relations_sha256`) in run.json's manifest info,
  warns on every stated hash that does not verify, and a resume refuses other relations under the same scenarios'
  hash (an arm of before the field: the relations of the manifest copy it kept decide).
- merge_run.py writes run.json `manifest_identity`: the run manifest's hashes and the stated ones it does not match,
  each arm's identity (from the manifest copy the arm kept, else its run.json) and `differs`, one sentence per arm
  whose content or relations hash is not the run manifest's (or that does not state one). In an exploratory
  run it is a warning; `merge_run.py --strict-identity` (run_suite.sh --gate) writes the files and exits 3, and
  gate.py fails on it. The report's environment table shows both hashes and names a mismatch.
- The async runner vendors suite_common.py, but its load_manifests compares its own scenarios hash and never calls
  verify_manifest: a re-vendor changes nothing at runtime; the runner's identity check is its own (the gate checks
  the arm's manifest copy).

### The sequential mode's MCP repeats
- The sequential mode ran the MCP arm's cold pass and the warm pass of repeat 1 only, so a sequential run never
  recorded MCP repeat 2 of the repeat-3 (D) scenarios that the paired plan gives (suite_common.mcp_repeats; found by
  VERIFY-A). It runs one pass per warm repeat index of the paired plan (1, then 2 for those scenarios), each followed
  by the paired loop's rules between steps (mcp_after: a cancelled job's wait; a refusal or three errors stop the
  further passes). Exit codes are unchanged; the report's plan reconciliation shows a pass that did not run.

### The run's plan and the strict deploy gate (finding 1: completeness and gate success were not enforced)
- **The plan.** run_suite.sh writes RUN_DIR/plan.json before any arm runs (`gate.py plan`): per invocation its
  selection (manifest identity, areas, limit, ids, warm_every, arms, mode) and every (scenario, arm, repeat) it plans,
  the paired loop's plan (plan_executions for a direct arm, mcp_repeats for the MCP arm; the sequential mode runs the
  same). A resume with the same selection adds nothing, one with another selection adds its plan: the run's plan is
  the union, so a resume never narrows what the run must record (an invocation adding no execution adds no entry). A
  resume without a selection of its own (--ids, --areas, --limit, --gate) runs every scenario of that union
  (`gate.py resume-ids`), so a plain resume after a narrower targeted one finishes the run, and keeps the last
  invocation's warm_every unless given (before lane S6c a resume without --areas widened the selection to the whole
  manifest and took warm_every 10). A session of orchestrator.json without a plan.json entry (a run started before
  plan.json) is reconstructed from its selection fields and joins the union; a run with no plan.json at all is
  reconstructed from orchestrator.json alone (persisted before its arms ran too). `--ids FILE` selects scenarios by id (an id the manifest lacks is refused); the ids reach both
  direct runners, the MCP steps and the plan.
- **The report** reconciles the plan with the records (`gate.missing_planned`): a planned execution no arm recorded,
  the MCP arm's included, makes the verdict INCOMPLETE ("N planned execution(s) not recorded (mcp n; plan of ...)").
  Before, only the direct arms' state was read, the MCP arm was excluded and compare iterated the observed scenarios,
  so a missing MCP execution or a scenario absent from every arm read PASS. A complete selective MCP-only run still
  reads PASS (it fulfils its plan).
- **The policy** (`gate_policy.py MANIFEST [--routes] [--arms] [--require-agreement] [--declared FILE]
  [--advanced-path] --out POLICY`, versioned `policy_version` 1) selects by ROUTE: POST /resolve, POST /traverse and
  GET /traverse/capabilities by default, i.e. 222 scenarios of v4 (81 R + 93 T + 41 N + 6 D + C-caps-0001; areas R/T
  alone miss the 47 N and D scenarios on those routes). The scenarios an arm cannot express are explicit exclusions
  with rule and reason, derived from the manifest and the by-design table, never from an answer a run received:
  `manifest_not_implemented` (22), `table_not_expressible` (a named case whose table agreement is
  not_implemented_mcp: 26, 7 of them beyond the manifest's: R-budget-0001/0002, R-edge-0003, R-twin-0002..0005,
  10 kb queries the service refuses at its own bounds.max_query_bp; a case that cannot change the request is moot,
  `not_expressible_moot`: the query-cap case on a query within the service's cap of 2,000 bp, 500 for a discovery,
  whose direct request states no bound, so R-edge-0012/0013, 100 bp refusals whose pair compare verifies, stay) and
  `advanced_path` (traverse_fetch overrides outside the service's basic list `_TRAVERSE_BASIC_OVERRIDES`, a preset
  other than walk / trie, a strategy object: T-labels-0004/0006). It lists the executions each arm must record (v4,
  derived exclusions only, without --declared: MCP 191 scenarios, 388 executions; the shipped two-arm policy with the
  declared file: MCP 189 / 384, direct-local 220 / 253), the relations required per arm (both scenarios selected and not excluded there: v4 REL-same_answer-0003 on both arms,
  REL-prefix_equals-0004 on direct-local;
  the prefix_equals relations of the excluded R-budget / R-twin and the N-zero / C-caps relations with a partner
  outside the selection are listed as excluded), the property checks not required (none unless declared, below), and
  pins content_sha256, relations_sha256, index_fp and the by-design table's sha256. gate/policy_v4_rt_mcp.json is
  the v4 deploy gate's policy, two arms with the owner's declared entries (below; `gate_policy.py
  manifest/suite_v4.json --name v4-rt-mcp --arms mcp,direct-local --declared gate/declared_v4_rt.json --out
  gate/policy_v4_rt_mcp.json`; a test regenerates it and compares). A v5 policy is one command on the v5 manifest,
  which prints the counts (selected by route and area, exclusions by rule with each reason, required executions per
  arm, relations required and excluded); `--ids-out FILE` writes the ids the policy's arms run (v4: 222), for
  run_suite.sh --ids or async's runner standalone (`gate.py policy-ids POLICY` prints them comma-separated).
- **The gate** (`gate.py RUN_DIR POLICY [--reevaluated]`: exit 0 PASS, 1 FAIL, 2 usage; RUN_DIR/gate.json) fails,
  listing every problem by kind, on: identity (either manifest hash not the policy's, a stated hash not verified, an
  arm or a plan on another manifest, another index_fp, another by-design table); missing (a planned or required
  execution not recorded, or a required one the run did not plan); outcome (on the policy's arms and selection: fail,
  error, finding, ANY inconclusive, not_implemented_mcp on a scenario not excluded, service_by_design by a case the
  scenario does not name, a by-design pair not verified when a direct arm ran it, the MCP side not verified when none
  did); property (a property check the scenario expects not passed on a passing or explained execution, e.g. not
  exposed; the line's outcome is not trusted: every check the scenario expects must have a result on each passing
  line, determinism excepted on the first repeat only, the status must pass and no other generic check may fail or
  be inconclusive there, so a pass line without its status result or a refusal passing with no check fails);
  disagreement (MCP vs direct when the MCP arm is gated, direct-local vs direct-remote when both are); agreement (with
  --require-agreement: an MCP class other than agree or a verified pair, a required property not compared); relation
  (a required relation not pass on a gated arm); configuration (below); nothing passed. An MCP-only policy needs no
  direct arm, but it cannot verify every by-design pair: a refused or relayed case is verified on its MCP side (a note
  says so), while an answered or clamped case that explains checks (the service answers what the host refuses:
  R-edge-0009, T-edge-0007, T-edge-0008) is verified only against a direct execution, so an MCP-only policy excludes
  those scenarios (`needs_direct_arm`, below) and the deploy gate runs two arms.
- **Configuration evidence** (agreed with async's runner on 2026-10-14), with the MCP arm gated: lane S6b's
  `compare.gate_configuration_problems` is the definition (a recorded live document stating service_settings, an
  empty `{}` or a failed read included; every referenced path equal live and in the table's blocks; no MCP line whose
  request_evidence read a reference from the table or found it unknown), and the gate adds: every invocation of the
  MCP runner recorded a successful read (run.json `invocations[]`, one problem for the run naming how many did not)
  and all reads are one document (a run spanning two deployments fails); the runner's own recorded mismatches
  (run.json, a line's request_evidence_mismatches); the live service_settings.advanced_enabled equals the policy's
  advanced_path where its advanced_path exclusions assume it (false) or the policy enables it; a request_evidence
  that is not {case: {path: live | table | unknown}} (compare's reader would raise on it; the gate states it). The
  per-line evidence not read live is one problem per reference path and source, naming its cases and executions
  (v4 replayed: 11 paths, the nine preset-row budgets of traverse.bounds_budget_above_preset_row on 14 executions,
  service_settings.advanced_enabled on T-labels-0004/0006's 4, service_settings.resolve.budget_ceiling_ms on
  R-edge-0021's 2); the gate checks that compare's per-line count and its aggregation agree.
- **run_suite.sh --gate POLICY**: without --ids / --areas / --limit the selection is the policy's (less what is
  excluded on every arm the invocation runs); merge_run.py --strict-identity, then compare, then gate.py, then the
  report (a "Deploy gate" section); the exit status is the gate's. `--reevaluate RUN_DIR --gate POLICY` judges the
  re-evaluated set (gate.reevaluated.json) and exits with the gate's status. Without --gate nothing changes: the
  exit status is the orchestration's (a FAIL report exits 0). gate.json records the sha256 of what it judged (the
  merged lines, the comparison, the manifest, plan.json); a report whose files differ, e.g. after a later resume
  without --gate, shows "GATE verdict STALE" with the files that changed, never the old PASS.

### A re-evaluation of re-normalised answers says so (finding 10)
- reevaluate.py judges the answers a run recorded and never re-runs a normaliser; a normaliser fix (G06: a pattern
  refused in its slot) reaches a recorded MCP arm only through async's suite_mcp_renormalise.py, which rebuilds the
  answers from the recorded tool payloads into a copy. Re-evaluating the original v3 arm gives 337 changed records,
  the re-normalised one the S5e totals (353; MCP 1071 / 2 / 16; 542 agree / 15 disagree): the difference is the 16
  re-normalised G06 lines (VERIFY-A), and nothing in reevaluated.json or the report said so.
- merge_run.py records what a re-normalisation left in an arm (`renormalisation(arm_dir)`: renormalised.json's
  source, time, normaliser sha256, digest rule, counts, changed and skipped executions and, where it records them, the
  renormaliser's commit and runner / client / input hashes; the lines' `renormalised` marks; the raw records the
  lines name that keep the recorded answer as `answer_recorded`) in run.json `renormalisation`; a plain merge warns
  that those lines' outcomes are the recorded ones, stale until a re-evaluation; the merge of the re-evaluated set
  copies it into run.reevaluated.json `reevaluation.renormalised`. The report's verdict names the arm, and its
  Re-evaluation section states by what and when it was re-normalised, how many answers were rebuilt and changed
  (which) and not rebuilt (by kind). (reevaluate.py, lane G2's, may call merge_run.renormalisation to put the same
  into reevaluated.json itself.)
- The Re-evaluation section states the changed records and the changed outcomes apart: "353 execution record(s)
  changed: 334 changed their outcome, 19 only their checks" (the review request's "353 execution outcomes" was
  imprecise; the report's columns were not).
- The documented sequence (offline; the recorded run is never written; COPY, RUN, SUITE are paths):

      rsync -a RUN/ COPY/
      /Users/raetsch/git/services/metagraph-search-service/.venv/bin/python -I -B \
          /Users/raetsch/tmp/mgss/wt-suite/scripts/bench/suite_mcp_renormalise.py COPY/arms/mcp \
          --out COPY/mcp-renormalised --checks-dir SUITE
      mv COPY/arms/mcp COPY/mcp.recorded        # outside arms/: merge_run reads every directory under arms/
      mv COPY/mcp-renormalised COPY/arms/mcp
      SUITE/run_suite.sh --reevaluate COPY --no-pdf --reeval-mcp-auth anonymous   # a run before async's 2a8acbaf

  async's renormaliser (feat/suite-mcp-runner-vendored) rebuilds what its payloads allow and copies the rest
  (not implemented, documents not kept, runner errors), marking every line; its outcome fields stay stale until the
  re-evaluation, which the report then states as above.

### The report's wording
- The front cost (Overheads) is taken on pairs of host-warm executions: "both executions meet the recent-repeat
  heuristic (another execution of the scenario ended on the host within 60 s before); cache residency and foreign load
  were not controlled". It said both arms answered from a warm row cache, so the difference was not a cache effect: the
  host_cache value is a recency heuristic (suite_common.assign_host_cache), not a measured residency (VERIFY-A).
- The resource figures are labelled contextual telemetry: the section is "Resource telemetry per area (contextual)",
  and the load caveat's figures say so; they describe the shared staging container, not the suite's requests.
- The MCP agreement is split by how much of a scenario's claim was compared (`report.agreement_depth`, VERIFY-B's
  classes, on property checks; budget_respected, which the MCP never sees, is generic and does not count): fully
  checked, partial (naming the property checks not compared), partial with a non-pass aggregate, on failure (every
  property compared, a non-pass aggregate: the arms fail alike), generic checks only, refusal; the verdict line gives
  these counts beside not comparable, disagree, by design verified / unverified, not implemented, the scenarios
  without an MCP execution and the planned MCP executions not recorded; an agreeing scenario whose reference fails a
  check the MCP does not expose is named apart. v3-S5e: 542 agree = 480 / 26 / 16 / 1 (Q-more-0003) / 8 / 11, with
  T-labels-0004 agreeing while both direct arms fail label_dict_contains; v4: 534 = 479 / 24 / 14 / 0 / 7 / 10.
- The Re-evaluation section counts the changed records and the changed outcomes apart (above).

### Validation on the recorded runs (offline) and what the async runner sees
- Fresh copies of the v3 S5e set (lane-s5e/reeval-S5e: the re-normalised MCP arm) and of the complete v4 run
  (v4-run-copies/20261014T061053Z_refseq97-6f945db4), re-evaluated at fe563f38 and with this lane (`run_suite.sh
  --reevaluate COPY --no-pdf --reeval-mcp-auth anonymous`): 0 differences in outcome, checks, case, explained checks
  and items, verdict_detail, inconclusive_reason, time_limited and outcome_runner over 2,900 and 2,913 lines, 0 MCP or
  direct class changes, 0 relation changes, the same summaries and verdict lines, no merge warning added; the plan gap
  is 0 for both (the as-run v4 copy misses one direct-remote execution, which its verdict already showed). The v3
  report now states its re-normalised MCP arm (1,128 of 1,330 rebuilt, 16 changed: the G06 lines) and "353 records:
  334 outcomes, 19 checks only".
- With c73f0f12's MCP-only policy, `run_suite.sh --reevaluate v4-copy --gate POLICY` exits 1, GATE FAIL with 29
  problems (the two-arm policy's verdict is in the next section): the
  MCP arm's run.json records no service_capabilities (a runner of before the record); 23 inconclusive MCP executions
  (R-pos-0027 r0 resolve_timeout; N-hard-0021 r0, N-realrep-0002 r0, N-realrep-0025 r0 traversal_resolve_deadline;
  checks: T-walk-0004 r0, -0006 r0 r1, -0012 r0 r1, -0014 r0, -0018 r0 r1, -0022 r0, -0036 r0 r1, T-repeat-0002 r0,
  T-repeat-0004 r0, T-labels-0003 r0, T-budget-0003 r0 r1, N-hard-0022 r0, C-caps-0001 r0 r1); 5 required properties
  not exposed via MCP (label_dict_contains on T-labels-0003 r1 and T-labels-0005 r0 r1; answer_fields on T-edge-0008
  r0 r1). Nothing is missing, no identity or relation problem, no disagreement, every by-design pair verified. Judged
  on its recorded lines (not re-evaluated) it fails on 31: the recorded comparison's older by-design table, the
  configuration, the 23 inconclusive, T-edge-0008's two recorded fails (traverse.duplicate_seeds, the old rule) and its
  MCP disagreement, the 3 label_dict_contains; the as-run copy adds P-iupac-0001 direct-remote r0, planned and not
  recorded (32).
- The async runner vendors suite_common.py; its load_manifests compares the scenarios' hash only and never calls
  verify_manifest, so a re-vendor changes nothing at runtime (manifest_identity / manifest_problems are there to use).
  No verdict rule of the mirror changes. gate.py reads the runner's agreed records (run.json service_capabilities,
  the lines' request_evidence) and renormalised.json's commit / runner / client / input hashes where it records them.

### The two-arm gate and the owner's declared list (the verification of lane S6c, VERIFY-G3)
- The MCP-only policy could pass no MCP-only run: the answered pairs R-edge-0009 (resolve.duplicate_labels),
  T-edge-0007 (traverse.seeds_above_host_max_seeds) and T-edge-0008 (traverse.duplicate_seeds) are compared check by
  check and verified only against a direct execution, C-caps-0001 is inconclusive via MCP (its one property check,
  capabilities_doc, is not exposed: the MCP answers the service's view, without `release`), and label_dict_contains
  (T-labels-0003/0005) and answer_fields (T-edge-0008) are never exposed via MCP. Design decision (coordinator and
  design lead, 2026-10-14): **the deploy gate runs two arms, mcp and direct-local (paired)**. direct-local verifies
  the answered pairs and lets the gate catch an MCP vs direct disagreement; the properties the MCP never exposes stay
  required on direct-local and are declared not required on the MCP arm.
- `gate_policy.py --declared FILE` reads the owner's entries, a versioned JSON file beside the policy
  (gate/declared_v4_rt.json: `declared_version` 1, the manifest's content_sha256, `exclusions` [{scenario, arms,
  repeats (optional), reason}], `not_required` [{scenario, arms, check, reason}], `relations_excluded` [{relation,
  arms, reason}]). Each entry needs a reason, a selected scenario (a manifest relation) and (not required) a check the
  scenario expects; a file for another manifest is refused. An exclusion has the rule `declared` (with `repeats`
  only those repeats); describe() prints every declared entry with its reason. v4, what the MCP surface cannot show:
  C-caps-0001 excluded on the MCP arm; label_dict_contains on T-labels-0003 and T-labels-0005 and answer_fields on
  T-edge-0008 not required on the MCP arm. v4 suite errors that block every v4 gate on a quiet host too (triage-v4,
  to be removed with the v5 manifest): T-labels-0004 on direct-local (V02, its modelled memory budget), T-budget-0003
  on both arms (G18, its work units below the seed's cost), T-attempts-0002 repeat 1 on direct-local (V04, the planned
  409 duplicate), REL-prefix_equals-0001 / -0003 on direct-local (V05: a stated stop answers the prefix, never
  testable on this host) and REL-prefix_equals-0002 on direct-local (its own reason: eco's full 10 kb, R-twin-0002,
  never completed on this host, so the full-twin profile it compares against does not exist). REL-prefix_equals-0004
  (twin half against full twin, no stated stop) stays required on direct-local: it passes on direct-remote in v3 and
  v4. The owner reviews the list (lane G3's DECLARED.md) before the policy is pinned.
- `needs_direct_arm`: a policy without a direct arm excludes on the MCP arm every scenario naming an answered or
  clamped case with explains_checks (from the manifest and the table: the service answers in place of the host's
  answer), so a clean MCP-only run passes; v4 MCP-only without --declared: 35 exclusions (R-edge-0009, T-edge-0007,
  T-edge-0008, C-caps-0001 besides the 31), 187 scenarios, 380 executions; with the declared file: 36 exclusions
  (T-budget-0003 besides), 186 scenarios, 378 executions.
- The advanced_path rule reads its list of basic override paths from the by-design table (the request predicate
  `overrides.` not_in [...] of traverse.advanced_path_disabled, lane S6b), so the policy and the comparator share one
  list; gate_policy.MCP_BASIC_OVERRIDES (the service's `_TRAVERSE_BASIC_OVERRIDES`) is only the fallback of a table
  without that case. gate/policy_v4_rt_mcp.json is regenerated after merging lanes S6g1 and S6b: the same
  selection, exclusions, required executions and relations; the pin is the merged table (sha256 ab4100c9).
- v4's two-arm policy: 222 selected; MCP 189 scenarios / 384 executions (31 derived exclusions, the declared
  C-caps-0001 and T-budget-0003), direct-local 220 / 253; REL-same_answer-0003 required on both arms,
  REL-prefix_equals-0004 on direct-local; REL-prefix_equals-0001..0003 excluded on both (on MCP derived, on
  direct-local declared).
- The v4 copy (re-evaluated) under it, on this branch with lanes S6g1 and S6b merged: `run_suite.sh --reevaluate v4
  --gate gate/policy_v4_rt_mcp.json` exits 1 with 45 problems, 0 structural (with a caveat: 16 of the
  load-dependent executions, V06's four 10 kb walks on both arms (12) and V05's four cold direct-local r0 503s, may
  not complete on a quiet host either; the triage files them under ENVIRONMENT and coverage, the v5 manifest fixes
  both). Configuration and evidence 13 (v4
  predates the runner's records): no live service_settings recorded; 1,333 of 1,333 MCP invocations without a read;
  11 reference paths the replay read from the table (the nine preset-row budgets on 14 executions,
  service_settings.advanced_enabled on 4, service_settings.resolve.budget_ceiling_ms on 2). Load-dependent 32,
  which the strict rule keeps blocking: REL-prefix_equals-0004 inconclusive on direct-local (R-twin-0004's cold r0
  503), MCP 19 inconclusive (R-pos-0027 r0 resolve_timeout; N-hard-0021
  r0, N-realrep-0002 r0, N-realrep-0025 r0 traversal_resolve_deadline; time-limited checks on T-walk-0004 r0, -0006
  r0 r1, -0012 r0 r1, -0014 r0, -0018 r0 r1, -0022 r0, -0036 r0 r1, T-repeat-0002 r0, T-repeat-0004 r0, T-labels-0003
  r0, N-hard-0022 r0) and direct-local 12 (503 deadlines on R-pos-0027, R-budget-0002, R-twin-0002, R-twin-0004 r0
  and R-budget-0001 r0's finding, a 503 past its budget, triage G22; time-limited walks on T-walk-0006, -0012, -0018,
  -0036, T-repeat-0004, T-labels-0003, T-coords-0003 r0). A clean, complete run built from the same copy (those
  executions as a quiet host gives them: MCP host states replaced by the scenario's passing repeat, direct-local 503s
  by direct-remote's answer to the same request; the not-exposed properties and the suite errors left as recorded;
  the table's service block recorded as the live document and the copy re-evaluated against it, so its references
  read live) passes this policy (exit 0) and the MCP-only policy derived from the same declared file (exit 0);
  without the declared entries both fail on exactly the entries above.

## Lane S6g1: checks that accept wrong answers (2026-10-14)

The verified findings 3-6 of the outside review of 2026-10-14 (VERIFY-C), the check part of the v4 triage's V03
and V05 (VERIFY-V4): each rule below closes a one-property corruption of a valid recorded answer that the checks
of fe563f38 accepted. The runner vendors checks.py: it takes these rules by re-vendoring.

- /traverse, a structural end (SPEC-labeled-traversal-core §6.4, §6.6, §6.7; `_structural_end_problems`): a path's
  `branch`, quorum, split-limit or re-entry end is accepted only where the answer shows it can have happened, in
  `walk_in_record` (a route along the record ending short of its minimum at such an end) and in the invariants (every
  path of every requested arm), on the full and the graphlet form alike. A structural reason the path counts has a
  label_end event of it on its leaf at the path's end (judged where the leaf states its events). `branch`:
  max_label_branches finite, the lineage's branches at the end equal to it (an ambiguous source ends `branch` once
  branches + 1 exceeds the allowance, and no lineage holds more), and two or more structural successors where the
  event states them; without a quorum knob also a route to the leaf through at least that many ambiguous splits (each
  branch taken is then a split; a split whose kind is not stated counts as ambiguous). Under a quorum knob
  (min_successor_labels >= 2, min_successor_fraction > 0 or min_live_labels > 1) a branch can be taken without a
  split: the walker raises the lineage's branches on every admissible successor before the quorum drops the minority
  ones, and a single followed successor records no split (walker.cpp at fe563f38), so the split count is not required
  there. `minority` (also `branch` with that text): min_successor_labels >= 2 or min_successor_fraction > 0;
  `below_min_labels`: min_live_labels > 1; `split_limit`: max_splits_per_path finite and a route through that many
  splits; `edge_reuse`, `edge_reuse_rc`, `rejoined_seed`: a blocked event of that reason on the leaf at the path's
  end ("blocked continuations are always reported"). The knobs are the walk's: a stated clamp's effective value, else
  what the arm's requests declare (a direct arm: the scenario's request; the MCP arm: the overrides of the scenario's
  MCP arguments, which it sends beside its preset, then the scenario's direct request, whose walk the preset
  mirrors), else the default of that request's exhaustive preset or annotate mode ("unlimited" branches and splits,
  no quorum: §5, §6.9), else the echo's, else the default (the echo's preset or mode, else §5's). The declaration
  comes before the echo, so a service that drops or rewrites an override, a preset's knob or a knob the request
  implies, and echoes what it walked, is held to the scenario, not excused by its own echo. On the MCP arm an echoed
  knob that differs from the declared one fails the invariants and walk_in_record by name: "the service walked with
  branching.max_label_branches 1, the scenario declares 2: a preset mismatch" ("... the scenario's MCP overrides
  state 2: an override not applied"); a deliberate preset change needs the manifest to follow it. A leaf with an end
  that cannot have happened shows none of its ends: walk_in_record fails "a stated stop that cannot have happened".
  The graphlet reader takes the lineage's values from the R records (structural_successors on the label_end events,
  `end_labels` with the branches at each leaf), as the full result states them. The arm's route split table is built
  once per arm (a 4,096-leaf trie: 0.05 s for the invariants, 5.4 s when it was built per path). Closed: GPT's
  zero-extension mutation of T-derived-0001 (both arms a 0 bp root, every label `branch` with branches 0 under
  max_label_branches 2; the dead_end control fails as before), branches 0 or 1 of 2, branches 2 on a route without
  ambiguous splits, an echo of max_label_branches 0 under a request of 2, a quorum or split-limit end without its
  knob, a re-entry without its blocked event, a graphlet's R branches or G split flags zeroed. Measured on every
  recorded /traverse answer of v3 and v4 (647, all arms: 447 branch, 16 below_min_labels, 32 edge_reuse, 12
  rejoined_seed label ends) and on the engine's stored traverse fixtures (tests: every stored answer states possible
  structural ends): none violates a rule.
- /pattern, a stated reason tied to the request and the phase that allow it (SPEC-pattern-search §7.4-§7.7, §12.1,
  §14, §19.6-§19.8, §20.5, §20.9; `_pattern_reason_problems`, which the invariants and every check that turns a
  missing result into inconclusive read: contains_source_context, source_absent, source_label, source_placement,
  path_spells, path_record_verified, path_absent, predicate_selected, predicate_either_mirror). Besides "Stated
  reasons": a request reads annotation only with labels all in a retrieval mode (§14.1: mode count reads none
  whatever the projection), a supported-path search or a predicate (in every mode), and has a memory account only
  with labels all in a retrieval mode, a predicate or a path search (§7.6 "Memory": mode count is the label-free
  path); so `rows_refused` and `annotation_budget` need the reads, `output_budget`, `cut: max_memory` and a
  `max_memory` stop the account; `anchor_labels_truncated` needs labels all; `predicate_above_threshold`,
  `selected_above_threshold`, `predicate_budget`, `cut: max_predicate_contexts` and a `max_predicate_contexts` stop a
  predicate; `count_above_threshold` is not stated with a predicate on L <= k, unless its normal form is the constant
  `true` (`selection.pass: "constant"` with `predicate.normal_form: true`: it runs the unfiltered request under the
  engine's own max_contexts, §19.6; `false` selects nothing and withholds nothing); `anchors_above_threshold`, `cut:
  max_paths` and a `max_paths` stop need a path search; `cut: max_annotation_work` a supported-path search (labels
  all states that stop without a cut, §7.5); a `max_anchors` stop a pattern longer than k (with or without a path
  search, §7.7), a `max_contexts` stop one of L <= k or a predicate; `paths_later_increment` is not stated with the
  anchors `exact` 0 (that empty answer is complete); only `all_or_count` withholds, but for `paths_later_increment`
  and `anchors_above_threshold`. A step or threshold stop leaves the count it interrupted not `exact` (§7.4):
  `discovery_budget` with an exact primary count fails, and so does a stop_at_threshold stop whose count
  (max_contexts: the raw contexts in discovery, the selected ones in the selection; max_predicate_contexts: the raw
  contexts, on supported paths the held walks' paths; max_anchors: the anchors; max_paths: the paths, with a
  predicate the selected ones) is exact; `threshold_crossed` needs such a stop. On a masked BASIC or CANONICAL graph
  a discovery stop of max_contexts (no predicate) or max_predicate_contexts (a predicate) on L <= k states its
  contexts `at_least` above the threshold: the engine stops once the running lower bound, the sum of the counted
  ranges' lower bounds, passes it, and the halted count states that sum (pattern_search.cpp at fe563f38:
  threshold_crossed compares running_lower on a masked graph; add_item adds each range's lower(); tally and to_count
  make the halted offsets `at_least` of exactly those lowers, the scans not having run; between mismatch levels the
  bound is recomputed before the comparison; ranges counted after the halt only raise it). Not judged without the
  mask (the running upper bound decides; fixture unmasked_stop_at_threshold states at_least 0) nor on a PRIMARY graph
  (an even-k wrapped union's lower bound need not be the stated value, §18 E2-02). Closed: GPT's P-modes-0005 with
  every count exact 0 and its threshold stop kept, at_least 2 under max_contexts 3, the stop removed; P-scope-0001
  (label-free) with its list removed behind annotation_budget (a refused row or a work stop), output_budget,
  anchor_labels_truncated, predicate_budget, or discovery_budget with the count left exact. 0 of the 3,567 recorded
  entries of v3 and v4 (all arms) and none of the engine's stored pattern fixtures states a reason these rules
  refuse.
- /pattern, a result listed twice (SPEC-pattern-search §7.1, §12.1, §21.4; the invariants): a context is one k-mer at
  an offset in an orientation at one level, (node, kmer, offset, strand), a path its nodes in an orientation, (nodes,
  sequence, strand); each is listed once. The order check accepts equal keys and counts_consistent compares the
  list's length, so a copy in place of a distinct result passed both. Closed: P-modes-0002 with its second context
  replaced by a copy of the first, L-long-0006's path list [+, +] (its '-' path dropped), P-scope-0001 with its count
  and partitions doubled and every context listed twice (with every exact count anchored by a complete list, a
  doubled count can no longer be padded). 0 duplicates among the 358,157 recorded results of v3 and v4. A duplicate
  page boundary of the MCP's paging is the runner's to test (its normaliser).
- /pattern, path_record_verified places the record (SPEC-pattern-search §14.3, §20.6): the occurrence of the
  record at the window's coordinates and strand states the index's seq_id and the record's length (as
  source_placement asks of a context), and every occurrence the label verifies lies in its own record (1 <= start
  <= end <= nt_length; source_placement applies the same rule to its record's occurrence). Closed: L-verified-0001
  and X-peptide-long-0001 with the record's occurrence at seq_id 999999999 and/or nt_length 1, nt_length + 1000, or
  another occurrence outside its record. All 14,298 recorded occurrences of v3 and v4 lie in their records, and
  every one of the truth's record states its seq_id and length. (L-long-0006's `path_spells` for its '+' path is a
  manifest change, left to the v5 lane.)
- /traverse, walk_in_record attributes a cap to the record's route only at the stop frontier (the v4 triage's V03,
  confirmed by VERIFY-V4): an arm `truncated` by a resource cap stopped the longest route along the record only where
  that route reaches the arm's `complete_to_bp` (or the arm states none, or the route's own leaf states the resource
  end); a route that ended before it ended for its own stated reason (SPEC-labeled-traversal-core §6.10: "every path
  that ended before complete_to_bp ended for the reported semantic reason") and is judged by it, so a `label_lost`
  short of the minimum fails instead of reading "stopped by time_budget". Where the record's next base after that
  route is no base of the graph (an IUPAC code: no k-mer holds it), the route ends there legitimately and the truth's
  minimum did not count it: inconclusive, "a truth gap" (the manifest's cut of the flank at such a base is the v5
  lane's). Measured: T-walk-0018's left arm (the record's 'S' 4,899 bp from the seed, its route ending label_lost at
  4,898 bp below complete_to_bp 5,250 / 7,378 / 8,748) on v4 L0, R0, M1 and v3 R0, M1: the detail now names the gap
  instead of the time budget, the outcome (inconclusive, its right arm time-stopped) is unchanged; the same answers
  against the record with an A there fail. No other recorded verdict or detail moves. A leaf that states a resource
  end below complete_to_bp cannot have had it (the same §6.10 sentence: a cap ends walks at the level it trips, at or
  past complete_to_bp): walk_in_record and the invariants fail it (it read inconclusive); every one of the 346
  recorded resource-ended paths of v3 and v4 and the 62 of the engine's fixtures ends at or past its arm's
  complete_to_bp.
- /resolve, a stop the scenario is built to state (the v4 triage's V05 as VERIFY-V4 amended it): R-budget's
  property is the stated stop (resolve_runs `prefix_if_stopped`: SPEC-labeled-traversal-core §4.5, a stop answers
  the resolve of a prefix). A complete answer without a stop holds the full profile but never tested the prefix:
  resolve_runs reads INCONCLUSIVE, "no stop happened: the prefix property was not tested" (a complete profile that
  is wrong still fails), where it passed. The relation prefix_equals whose a-scenario asks `prefix_if_stopped`
  (REL-prefix_equals-0001 / -0003) reads the same when a states no stop and the two profiles agree (they still fail
  when they differ). A twin's `unless_stopped` (R-twin-*: the reference profile) and the continuation relations
  (-0002 / -0004) keep their complete-profile tests: there the complete profile is the property. Recorded runs:
  every direct R-budget answer of v3 and v4 is a 503 and every MCP one a 400 (the 2,000 bp limit), so no outcome
  moves; R-twin-0004 R0's complete answer judged as R-budget-0002's reads inconclusive (it passed), and so does
  REL-prefix_equals-0003 on it.

### Validation (offline)
- Mutation probes (lane-g1/scratch/probes/mutations.py: valid recorded answers of v3 and v4, their one-property
  corruptions, GPT's probe answers, legitimately incomplete answers): 75 probes; at fe563f38 46 of them read
  other than they must (every corruption above passed or read inconclusive), here all 75 do: every corruption
  fails, every valid answer passes, every legitimately incomplete one stays inconclusive. The lane's verification
  (its probes, lane-g1/verify/probes/my_mut.py, 43 boundary corruptions and legitimate answers) reads as it must
  on 42; the 43rd judged T-quorum-0001's MCP answer against a direct request changed to min_live_labels 1, which
  the MCP arm, held to its own overrides (min_live_labels 2), rightly passes; the same change of the MCP overrides
  fails. Its two legitimate answers the lane's first rules failed (a branch end under a quorum knob after a branch
  taken without a split; count_above_threshold with a constant `true` predicate) pass. Its follow-up probes
  (lane-g1/verify/fu/probes/fu_mut.py) read as they must: T-walk-0017's MCP answer with the walk preset's
  max_label_branches lowered to 1 or 0 and echoed fails as a preset mismatch, an MCP override the service did not
  apply fails as such, and `selection.pass: "constant"` beside a normal form other than `true` claims no exemption.
- Re-evaluation of the recorded runs (v3: lane S5e's reeval-S5e copy; v4: the complete run copy; `run_suite.sh
  --reevaluate COPY --no-pdf --reeval-mcp-auth anonymous`, against the same copies re-evaluated at fe563f38): 0
  outcome changes and 0 check outcome changes of 2,900 (v3) and 2,913 (v4) executions, 0 relation changes, compare
  summaries identical; detail changes only on T-walk-0018's walk_in_record (v4 L0, R0, M1; v3 R0, M1), as stated
  above.

## Lane S6b: exception matching, timing plausibility and the comparator (lane G2, 2026-10-14)

The verified findings 2 and 7 of the outside review of the suite (VERIFY-B), the triage's V01, the compare decision
agreed with the async session, and the amendments of the lane's own verification (VERIFY-G2). The async session's
runner implements the same matching and plausibility rules from the lane's TABLE-DIFF.md; both read
service_by_design.json at sha256 ab4100c9...

### By-design cases state the request they are about (finding 2)
- A case may state `request`, the machine-readable form of its scenario_match: a non-empty list of predicates that
  must all hold on the MCP arguments the manifest scenario states (never the wire) for the case to apply, named or
  answer-triggered alike (compare.request_holds). Before, an answer alone decided: a ">= 1" refusal of a cap the
  request set to 10,000 (or never sent) read as pattern.caps_zero, a refusal of a basic traverse override as the
  advanced path (every T scenario not_implemented_mcp), a budget refusal of max_extension_bp as the preset row; a
  regression of the service became a verified by-design pair and the gate saw no failure.
- Grammar: `{"any": [[p, ...], ...]}`, or one subject (`path` into the arguments; `field`, the refused field's name
  less a prefix; `service` / `host`, a reference) and one operator (present, absent, equals, in, not_in, gt, ge, lt,
  le); a value may be a reference `{"service": path}` / `{"host": path}`. Pinned so that two implementations read it
  alike: a path is keys joined by ".", `[n]` / `[*]` after a key, `{field}` a whole key (the refused field, the tool
  error's else the answer error's, less `field_prefix`: pattern_search's arguments are flat, its errors name caps.*),
  validated at load, and a substituted field that is no path reads absent; a list at a path is one value unless `[*]`
  steps into it (never into an object's values); null and absent are one (equals null and null in an in / not_in
  list are refused at load); equality is JSON equality (true is not 1, also inside lists; 250 and 250.0 are one); an
  absent path fails every comparison but not_in; a refused field equal to its prefix names no key; an unknown
  reference fails its predicate; `<...>` field placeholders match one or more characters, never none.
- References: the host's /capabilities (`host`) and the service's GET /traverse/capabilities root (`service`:
  service_settings.*, presets.<preset>.budgets.<knob>, service_maxima.*, backends.<db>.resolve.time_budget.*), resolved
  from the run's recorded live document (the MCP runner's run.json `service_capabilities.document`; the host's
  /capabilities the direct arms recorded), else the table's own blocks `host_capabilities` / `service_settings`, else
  unknown (suite_common.References; the order agreed with async's runner). A path both state with different values
  is a mismatch, judged with the live value. Each MCP line records `request_evidence` {case: {reference: live |
  table | unknown}} (and `request_evidence_mismatches`); compare.json's `references` and a warning line in the report
  name the source of each namespace and every mismatch; reevaluated.json records the same. Lane G3's gate calls
  compare.gate_configuration_problems (no live service_settings recorded, a mismatch, a line that read the table or
  found a reference unknown).
- A case its request holds back does not apply (matching continues; the line's verdict_detail notes once "<case>:
  request predicate <predicate> does not hold (field <f>)", also for a case both named and answer-triggered).
- Load validation (both implementations refuse the table): a malformed request, a reference the table's own blocks do
  not state, a missing request where one is required (an answer-triggered case without a wall rule; a case whose MCP
  field names several alternatives), and (compare.py's command line, reevaluate.py: the run's manifest) a scenario
  that names a case its own arguments do not satisfy for any field the case allows.
- The table: the two blocks (staging's values: protein mismatch cap 1, reserves 250, /resolve max 60,000 ms; the
  service settings of async's 7a65e913; the preset rows and service maxima as staging's traverse_capabilities stated
  them in the v4 run, walk 15,000 ms / 512 MB / 1e8); thirteen cases state their request (caps_zero,
  advanced_path_disabled, empty_seeds, above_kind_cap, the preset row, the service-maximum knobs, five pattern refusals
  with field alternatives, resolve.time_budget_above_host_max, resolve.time_budget_unnamed); caps_zero, the preset row
  and the service-maximum case name their fields one by one instead of `caps.<cap>` / `bounds.<knob>` (the
  service-maximum case: the knobs presets._check_knobs bounds by service_maxima that a basic request can reach,
  bounds.max_extension_bp / max_steps / max_paths / max_output_bp / max_live_paths and labels.max_labels_per_node);
  the advanced path never reads a refused "overrides." or "overrides.labels"; the stale "20000 on staging" reads
  60,000; tighter texts where the recorded refusals state them. Every frozen manifest satisfies the requests it names.

### Timing plausibility from the request, the host and the service (finding 7)
- The bounds a rule compares with come from the request, the host's capabilities and the service's settings (the same
  references); the bounds a failing answer states are checked against them before any is used (a stated phase's
  bound that passes sets the least wall, bound - 1 s; the legacy rule's stated timeout can only raise a term the
  write cap then caps): the /resolve budget is the request's, an explicit-labels resolve's without one the service's
  default (60,000 ms); the client timeout t = max(the mode's base 20 s / 120 s, budget + the host's reserve + 10 s).
- resolve_timeout: a stated read bound or deadline window below t - 1 s ("service read bound ... below its own floor
  ...", "service deadline window ... below its own floor ..."), a write bound other than min(10 s, t) within 1 ms
  ("service write bound ... is not its write cap ...") is no host state whatever the wall; a tool error stating a
  reserve above the host's or a budget other than the one sent ("stated finalisation reserve ... above the host's ...",
  "the tool error's time_budget_ms ... is not the request's ...") is no host state (resolve_timeout and the relayed
  /resolve deadline). The relayed deadline's threshold is budget - the host's reserve. A legacy record (no stated
  phase) without a budget takes the client timeout as a term (suite_common.legacy_unbudgeted_client_timeout_term:
  threshold 9,000 ms instead of the answer's "within <N> s" - 1 s or rule (c); pending the async session's agreement).
- pattern.deadline_ends_task's wall rule: the request's budget - 1 s (the engine answers its deadline only once the
  budget has elapsed) on the failing attempt's own wall (the deadline-ended task's one HTTP try; a traversal job's last
  settled attempt), never mcp_phases.backend_ms (a sum over a search's tasks or a job's attempts); a reserve the answer
  states above the host's keeps the case from applying. Rule 2 and rule (c) are unchanged.
- Three tests asserted answer-trusting bounds (a 5 s write bound, a 15 s deadline window on a 20 s budget, a tool
  error's 5,000 ms budget for a 20,000 ms request) and are corrected, with the read floor's and the explicit-labels
  deadline's expectations that followed the tool error's budget and the legacy rule's without a budget.

### Comparator
- A not_implemented_mcp match (answer_match.mcp_outcome) whose pair is a direct-side disagreement with
  direct_mismatch "disagree" makes the scenario's class `disagree` (side direct, the case and the direct execution
  named, the MCP executions' not_implemented_mcp stated in the detail the report shows); the MCP lines keep
  not_implemented_mcp; "timing" stays unverified (the agreed decision of 2026-10-14). Each repeat is classified on its
  own: the host's ordinary answer stays the reference beside a not_implemented_mcp group, and an MCP fault on another
  repeat is the scenario's class, never hidden under not_implemented_mcp.
- A named not-expressible case (resolve.not_expressible.bounds_max_query_bp) reads not_implemented_mcp only when every
  disagreement is the direct side of a case the scenario names; a check-level disagreement, an MCP side or a fault on
  any repeat stays a disagreement.
- A named case's direct side is judged on every repeat of both direct arms, as an answer-triggered case's is (the
  arm and repeat of a contradiction named, a contradiction by another status before one by a 503 / 429; a transport
  error reads neither way). For a case whose host answers 200, a 503 (deadline or busy) or 429 on a direct execution
  is the host's load and is set aside when another direct execution holds the side; it contradicts only when none
  does (a load-caused 503 on one arm would otherwise turn a verified pair into a disagreement that blocks a deploy
  under the strict gate). A 400 side keeps every 503 a contradiction: R-budget-0001 / -0002 and R-twin-0002 stay
  not_implemented_mcp by table.
- The direct arms' timing content rule tests every route of the stopped walk (any parent at a merge) as a prefix of a
  route of the complete walk, not the displayed first-parent paths (T-labels-0003 of the v4 run, SPEC-labeled-
  traversal-core §7.1), bounded by the checks' route-state limit.

### Validation on the recorded runs (offline)
- Copies of v3 (the S5e input set) and v4 re-evaluated (run_suite.sh --reevaluate --no-pdf --reeval-mcp-auth
  anonymous) against fe563f38: no outcome, matched case, check or MCP class moves; v4 T-labels-0003's direct pair reads
  timing (was the run's only direct disagreement). The deadline_ends_task lines state the new threshold (v3 6, v4 35
  lines); 22 lines per run gain request_evidence (T-labels-0004 / -0006 stay not_implemented_mcp, the setting read from
  the table's block). The verification's amendments move no verdict; seven compare disagreement texts now name the
  direct arm and repeat ("direct-local repeat 0: ...": R-budget-0001 / -0002, R-edge-0003, R-twin-0002 / -0003 /
  -0004 / -0005, in both runs), R-twin-0004's naming direct-remote's 200 before direct-local's 503 deadline.
- 263 tests pass (1 skipped).

## Lane S7v: the verdict findings of the two-day review (2026-10-14)
The verdict items of the two-day review (~/.cache/metagraph-traverse-reviews/2026-10-14-ultracode/REPORT.md) on
compare.py, reevaluate.py, report.py and suite_common.py's verdict parts, against 14552938 (lanes G1-G3 merged). Each
rule below changed because a recorded or probed input read wrong; tests/test_verdict.py holds one test per rule, each
failing on 14552938.

### The direct arms' class with equal digests
- `identical` was decided from equal digests before the checks and outcomes were read, so two byte-identical answers
  whose checks of the answer differed read identical and were left out of `direct_summary.disagreeing`. The order is
  disagree, identical, timing, not_comparable, same_checks, with one exception on a repeat whose answers are the same:
  a difference only in the execution checks (`compare.EXECUTION_CHECKS`: budget_respected, which reads the client's
  wall, holding the network and any wait for a request thread; deterministic, which reads the arm's earlier
  repeats) is no disagreement of the answers. It is recorded per repeat as `execution_differs` [{repeat, checks:
  [{check, local, remote}]}] (`compare.execution_only_difference`), the class stays identical, and the arm's own
  finding or failure counts in its totals (FINDINGS / FAIL through that execution; the strict gate blocks it there).
  Any other difference with equal digests (a check of the answer, e.g. arms judged by different checks.py; an outcome
  difference with no check differing; an outcome outside the vocabulary) disagrees. compare's direct summary lists
  `execution_differs`; the report's direct section names each with its checks.
- The verification of the lane (VERIFY-S-VERDICT (1)) refuted the first form of this rule, which made v4's
  R-budget-0002 and R-twin-0002 (byte-identical 503 deadline bodies; direct-local's budget_respected inconclusive,
  direct-remote's a finding: its client wall 9.5 s / 25.5 s past the budget) direct disagreements and v4's verdict
  FAIL: with both arms late the same pair read identical, so a worse host read better. They are identical with
  execution_differs budget_respected inconclusive / finding, and v4 reads FINDINGS as at 14552938.

### The MCP arm's reference, per repeat and check item
- check_agreement compared every MCP execution with one whole direct arm (direct-local whenever it had a line) and
  skipped every item that arm left inconclusive or not_exposed, so a direct-local 503 under load hid an MCP failure
  direct-remote had decided (class agree, or not_comparable when every item was skipped). Now the reference item of
  each check and occurrence is direct-local's unless direct-local left it undecided while the MCP item and
  direct-remote's item of the same repeat (else its first) are decided: then direct-remote's is taken
  (`compare.reference_checks`). The scenario's entry records `reference_alternative` and `reference_fallbacks`
  [{repeat, check, occurrence, arm, arm_repeat}] (only when an item fell back); `reference_outcome` is that of the
  lines the comparison used (a direct-local execution every compared item of which fell back is left out, the
  direct-remote executions that gave an item are in; without a fallback it is unchanged). The report's agreement
  depth reads the same items. The by-design pairs were already judged on every repeat of both direct arms (lane G2).
- Moves: v4 N-usecase-0008 and R-pos-0027 not_comparable -> agree (direct-local r0 a 503 deadline, direct-remote r0
  decides 6 / 5 items), T-repeat-0004 one more item compared (walk_in_record from direct-remote); v3
  M-planted-0017 not_comparable -> agree, T-walk-0004 one more item compared, and six v3 scenarios already
  `disagree` by an MCP runner_error now list the error line's failing items against direct-remote's
  (M-planted-0016, M-protein-0002 / -0007, N-nearmiss-0006 / -0007 / -0013). gate.py's uncompared_properties (require_agreement only; the shipped policy does
  not set it) still reads the reference arm alone, so it is stricter than compare where an item fell back.

### The host cache after a relayed deadline
- suite_common.reached_host read an MCP execution as reaching the host only with a number in server.elapsed_ms or
  mcp_phases.backend_ms; a relayed host deadline (502 search_failed "pattern_deadline: ...", 503
  traversal_resolve_deadline, 504 resolve_timeout, uninterruptible_overrun_exceeds_allowance) states neither,
  although the host worked 30-60 s, so it read not_reached and the next direct-local execution was labelled cold with
  since_prev_s null (its server time a median 2.48x the warm remote's, against 24.9x for a truly cold one). Such an
  answer (`relayed_host_deadline`: those codes, the relayed /pattern deadline's text, or a time_limited naming a relay)
  reaches the host now. `backend_timed` keeps the old test where a backend timing is what matters: the MCP overheads
  (a relayed deadline's wall is the deadline's) and whether a line records the service's answer.
- assign_host_cache adds `host_cache_after_interrupted: true` (present only then) to an execution whose previous host
  execution of the scenario ended at the host's deadline without its answer (`interrupted_on_host`: a direct 503
  deadline or a relayed one); compare's front overhead names such arms (`after_interrupted`), its summary counts
  `warm_pairs_after_interrupted`, and the report's Overheads section states the count (warm by the recent-repeat
  heuristic, the cache possibly only partly filled).
- Moves (timing only; no outcome, check, class or verdict): v4 39 MCP lines not_reached -> cold (27) / warm (12), the
  27 direct-local r0 after them cold -> warm with since_prev_s; 70 lines gain host_cache_after_interrupted; warm
  pairs 527 -> 537 (32 after an interrupted execution). v3: 10 MCP lines, 8 direct-local r0 cold -> warm.

### The text_contains grammar of an answer_match side (for async: the runner's load validation)
- A list-valued text_contains on an answer_match side raised TypeError in compare and reevaluate, and the load
  validation did not refuse it. The matcher (`compare.answer_text_contains`): on the direct side a string, or a list
  whose fragments must all be within one text (as a direct_expectation's text_contains); on the MCP side one string
  read as the runner reads it, str(text_contains) within the error's text (a list never matches there).
- The load validation (`compare.answer_value_problems`, in table_problems) refuses a side's code, text_contains,
  field or slot_code of another type: a non-empty string, on the direct side also a non-empty list of non-empty
  strings; on the MCP side one string (the MCP side's code stays the "states no code" rule). The shipped table
  passes, so nothing recorded moves. The runner reads the MCP side with str() and validates no type: it should
  refuse the same values (MIRROR, to agree with async).

### A runner's rule (c) without a table (for async: the no-table handling)
- With `--by-design ''` (no table), a runner's rule (c) inconclusive on a by-design case (inconclusive, its
  service_by_design_case, inconclusive_reason an undecidable plausibility) was not kept: reevaluate judged it by its
  checks (fail) and compare compared it check by check (disagree), while a runner's service_by_design mark and a
  by-design finding are kept in the same situation. `compare.runner_rule_c_by_design` names such a line; without a
  table reevaluate keeps it inconclusive ("decided by the runner's by-design case X by rule (c) (...); no table to
  replay it", the case, the runner's record and its inconclusive_reason kept) and compare sets it aside as
  undecidable (not_comparable when nothing else is judged; another repeat still decides). With a table nothing
  changes (the case is replayed); the recorded runs are re-evaluated with the table, so nothing moves.

### The runner's not_implemented_mcp mark (for async: how a new mark is announced)
- An MCP execution the runner marked not_implemented_mcp that no case reproduced made the scenario's class
  not_implemented_mcp ("by mcp runner") whatever the table said, while the comparator section above said a table
  judges such a mark. With a table, a mark on a scenario the manifest maps to a tool (`compare.mcp_mapped`) that
  neither mcp.not_implemented nor a named case of agreement not_implemented_mcp states, and that no case with
  answer_match.mcp_outcome matched, is `disagree`, side mcp, by mcp runner, the repeats named
  (`compare.unmatched_not_implemented_mark`; the other repeats are compared as usual). With no table, with such a
  named case, or on a scenario without an MCP tool, the mark stays the runner's class. The recorded runs have no such
  mark (their not_implemented_mcp lines are the manifest's or traverse.advanced_path_disabled's), so nothing moves. A
  new runner mark reaches the comparator through the table (answer_match.mcp_outcome) or the manifest.

### Outcomes outside the vocabulary (for async: the shared vocabulary)
- compare's totals counted only its known outcome keys, merge's stricter rule skips an unknown value and nothing
  checked a line's outcome against suite_common.OUTCOMES, so six MCP lines with `skipped_by_runner` (or an MCP-only
  run whose lines read so beside failing checks) gave PASS, against this README's "any other outcome reads as a
  failure". `suite_common.outcome_known` is the shared test (OUTCOMES: pass, fail, finding, inconclusive,
  not_implemented_mcp, service_by_design, error). compare: totals and per-area counts gain `other`, the summary lists
  `outcomes_outside_vocabulary` [{scenario, arm, repeat, outcome}], an MCP execution with such an outcome is a
  disagreement of its own (side mcp, by outcome; never compared check by check), a direct pair with one disagrees,
  and aggregate_outcome ranks it worst. reevaluate keeps such an outcome (the checks never replace it) with a
  warning. report: any `other` makes the verdict FAIL ("N with an outcome outside the vocabulary"), the arm's line
  and a list name them. gate.py already failed such a line ("outcome '...'"). merge_run (lane S-runner) and the
  paired loop's mcp_after (run_suite.sh) still read the value as recorded; their side is lane S-runner's.
- Recorded runs: every outcome is in the vocabulary; the totals gain `other: 0` and the summary an empty list.

### The load caveat without a measurement
- load_caveat read a missing sra_all figure as 0 and printed "mex was not under the Logan load" when every sample
  had failed (docker stats: a missing or renamed container, the daemon down; the sampler's `error`, its figures
  null), with load1 moderate or null. load_figures counts the samples stating each figure and the failed ones (the
  first error); with no sra_all figure, or more than 10 % failed samples (LOAD_SAMPLE_ERROR_SHARE), the verdict
  section says "**Load not measured:** ..." with the figures and the failures, never the positive sentence; a busy
  host keeps the Load caveat (with the failures named); the resource section names the failed samples. The recorded
  runs have no failed sample (v3 12,756, v4 8,693, all with sra_all): their sentence is unchanged.

### A field the MCP drops, and a part of the run that tested nothing
- A field the MCP answer lacks makes its check `not_exposed`; an execution whose property checks all read so is
  inconclusive, the pair agrees on the checks both arms evaluated, and the verdict was INCOMPLETE only when nothing
  passed at all: a service regression dropping `counts` or `results` read PASS with a footnote (the review's case:
  every MCP property check not exposed -> 6 inconclusive, class agree, PASS). Three layers judge it now:
  - **The strict deploy gate (lane G3)** already fails such a run: every inconclusive execution blocks, and a
    property check the scenario expects that is not passed on a passing execution (not exposed included) is a
    "property" problem unless the policy's declared file says it is not required on that arm (v4: label_dict_contains
    on T-labels-0003 / -0005 and answer_fields on T-edge-0008 via MCP). The declared list is the gate's record of what
    the MCP surface does not expose; anything else not exposed blocks.
  - **The report (exploratory runs)**: `report.untested_groups` names every area and route x arm whose executions
    include an inconclusive one and none passed (executions only not implemented via MCP or by design are marks,
    counted elsewhere). Such a group keeps the verdict from PASS: "INCONCLUSIVE (no execution passed in N area / route
    x arm group(s): ...)" after FAIL and FINDINGS; under another verdict a line names them. The agreement split of
    lane G3 (Agreement: generic checks only, partial) still names the scenarios whose properties were not compared.
    Consequence: v3 and v4 name area C on mcp (C-caps-0001 inconclusive via MCP: capabilities_doc not exposed, the
    10 others not implemented via MCP) and the routes GET /capabilities and GET /traverse/capabilities on mcp, so a
    clean exploratory run on v4 with the MCP arm reads INCONCLUSIVE, not PASS (the gate excludes C-caps-0001 on the
    MCP arm by the owner's declared entry); their recorded verdicts are FAIL and do not move.
  - **The comparator, with a baseline**: a scenario may state `mcp.exposes` (check names the MCP arm decided in a
    baseline run: `compare.mcp_exposure(rows)` builds it from a run's lines, for the v5 builder); a listed check that
    reads not_exposed on an MCP execution compared check by check is `disagree` (side mcp, by exposes, the checks and
    repeats named; `compare.exposure_drops`), so a dropped field fails the report and the gate's disagreement rule.
    No frozen manifest states it, so nothing recorded moves.

### A planned MCP arm that recorded nothing (the report side of a missing MCP execution)
- Lane G3's plan reconciliation already makes the verdict INCOMPLETE for a planned MCP execution no arm recorded, in
  both modes: run_suite.sh writes plan.json before any arm, and a run without it is reconstructed from
  orchestrator.json (proved here for the sequential mode with the real run_suite.sh and fake_mcp's exit3: "12
  planned execution(s) not recorded (mcp 12; plan of plan.json)"; the exit status and the loop record are lane
  S-runner's). Left: the Agreement section of such a run said "No MCP arm in this run."; it says "No MCP execution
  was recorded although the run planned N (see the Verdict: INCOMPLETE)" whenever the plan holds MCP executions.

### Validation on the recorded runs (offline)
- Fresh copies of v3 (lane-s5e/reeval-S5e) and v4 (v4-run-copies/20261014T061053Z_refseq97-6f945db4) re-evaluated
  with 14552938 and with this lane (`run_suite.sh --reevaluate COPY --no-pdf --reeval-mcp-auth anonymous`): 0 changes
  of outcome, checks, case, explained checks / items, verdict_detail, inconclusive_reason, time_limited or
  outcome_runner over 2,900 and 2,913 lines; the line changes are the host-cache labels (v3 18 / v4 66 host_cache,
  8 / 27 since_prev_s and host_cache_after, 54 / 70 lines gaining host_cache_after_interrupted).
- Class moves, each the point of its rule: v4 R-budget-0002, R-twin-0002 direct identical -> disagree; v4
  N-usecase-0008, R-pos-0027 and v3 M-planted-0017 MCP not_comparable -> agree; v4 T-repeat-0004 and v3
  T-walk-0004 one more item compared; six v3 scenarios already disagreeing by an MCP runner_error list their
  differing items. Summary: v4 MCP 534 agree / 4 not comparable -> 536 / 2, direct disagree 0 -> 2; v3 MCP 542 / 5
  -> 543 / 4. No relation moves.
- Verdicts: v4 FINDINGS -> FAIL (the 2 direct disagreements above); v3 FAIL unchanged. Both gain the line naming area
  C and the two capabilities routes on mcp where nothing passed. The deploy gate on the v4 copy
  (`--gate gate/policy_v4_rt_mcp.json`, two arms mcp + direct-local): 45 problems before and after, the same list.

### The host's references in an MCP-only run, a service document's identity, the wall rule's name (agreed with async)
Agreed with the async session on 2026-10-14 (the coordinator's item). The runner records, in run.json
service_capabilities and in each invocations[].service_capabilities, {document (GET /traverse/capabilities), host_document,
sha256, read_at, base, database, mismatches, error}; host_document has the shape of the host's /capabilities document:
"pattern" the /pattern/capabilities/<db> capabilities.raw, "resolve" {"time_budget": document.backends.<db>.budgets.
resolve_time_budget} (the host's block as received, not backends.<db>.resolve.time_budget, the service's terms).
- `compare.recorded_documents` (and so bind_run, the references every rule reads, reevaluate's replay and
  gate_configuration_problems): host = a direct arm's recorded /capabilities when the run has one, else the MCP
  runner's service_capabilities.host_document, read as live with origin "arms/mcp/run.json
  service_capabilities.host_document", else the table's host_capabilities block, else unknown; a path the live document
  lacks reads the table's block, else unknown (suite_common.References). An MCP-only run's host references, and the
  gate's live / table mismatch check on them, read the service's view of the host instead of the table.
- A service document's identity (`suite_common.service_document_sha256`, `service_document_view`): sha256 of
  canonical_json({service_settings, service_maxima, presets, backends: {<db>: {budgets, resolve, effective_maxima,
  limits}}, host_document}), a key the document lacks null in the view; probe times and counters elsewhere in the
  document change between the invocations of one deployment and would trip the gate's "spans more than one
  deployment". gate.py's `_document_sha` (edited for this function only) prefers the recorded sha256, else this view;
  compare's recorded_documents states the same. The runner vendors suite_common.py and can call the function itself.
- service_by_design.json, the vocabulary's answer_match.wall_rule entry only: the name budget_minus_reserve is kept,
  its threshold is B - 1 s (B the request's time_budget_ms) on the failing attempt's own wall, no finalisation
  reserve subtracted. The table's sha256 is d1e040aac81007713e2250a3dd765244e44026c8bcb03834fd484f143ccebc61 (was
  ab4100c9...); gate/policy_v4_rt_mcp.json pins the old one and is regenerated at integration
  (test_gate's test_shipped_v4_policy_is_regenerated_exactly fails until then, on that sha only).
- Agreed, no change: request_evidence_mismatches as the field name; the legacy no-budget term.
- Recorded runs: v3 and v4 record no service_capabilities, so nothing they read moves; their comparison and report
  name the table's new sha256.
