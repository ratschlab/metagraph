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
  1800 s, or the request's stated time budget + 30 s where that is more (suite_common.client_timeout_s; v6's
  /pattern budget of 110 s gets 140 s); direct-remote at most 125 s. A request whose stated budget reaches the
  proxy's cut less its tolerance (v6's 110,000 ms /pattern budget) and gets the proxy's own 504 at the cut (no host
  body) on direct-remote reads inconclusive on its status, as an allowed 503 deadline does, and is time-limited
  where the checks and the comparison look for one (checks.ingress_cut_allowed, checks.time_limited: the host's
  answer may still come within the window a late 503 reads inconclusive, so the 504 shows nothing of the host); a
  504 on a budget that ends before that, or on another arm, stays a status fail.
- Execution outcome (review of lane S1; the MCP runner applies the same rule): checks may report `pass`,
  `fail`, `not_exposed`, `inconclusive` (the scenario could not test its property: an allowed 503, a withheld
  answer hiding the source, a stated resource stop) or `finding` (a stated bound was broken on the clock the
  check reads: a 200 past its budget by the host's own elapsed time, a walk stopped outside its stated rule on
  the seed's clock, an allowed 503 deadline more than max(1 s, the budget) late on the client's clock; apart from
  fail whatever the HTTP status); any other outcome reads as a failure. The execution is `fail` when a check failed, else `finding` when one is a finding, else
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
  repeats on this arm that the host answered (`suite_common.feeds_previous`: an HTTP status and an outcome other
  than error; lane S7r), in repeat order, also across a resume, as `{"digest", "time_limited"}` (lane S3, the
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
   arm follows async's exit codes: the loop stops at once on the runner's exit 2 / 3 (a start refusal: key file,
   queue wait ran out), 4 (the by-design table refused), 5 (the run's snapshot refused) or 130 (interrupted),
   or after three MCP executions in a row that ended `error` (a
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
   `{checks_dir}` (the run's snapshot RUN_DIR/scripts, lane S7r), replaced with shell-quoted values; run with
   `bash -c` in this environment, unchanged (the runner reads its admin key from there); its output goes to the
   terminal and to arms/mcp/runner.log (lane S7r). A full run (not --smoke) takes the default
   template `DEFAULT_MCP_CMD` at the top of run_suite.sh, async's runner of branch feat/suite-mcp-runner:

       <search-service>/.venv/bin/python /Users/raetsch/tmp/mgss/wt-suite/scripts/bench/suite_mcp_runner.py \
         --manifest {manifest} --ids {ids} --repeat {repeat} --out {out} --check-queue \
         --queue-ignore sra-logan-chunks --checks-dir {checks_dir}

   so a full run is `METAGRAPH_ADMIN_KEY_FILE=/Users/raetsch/tmp/mgss/secrets/staging_admin_key
   ./run_suite.sh` (the key file is passed through, never logged; without it the script refuses to start,
   since the runner checks the service queue). `SUITE_MCP_CMD` replaces the default ('' for none), `--no-mcp`
   runs without the MCP arm, `--mcp-cmd default` asks for it in a smoke run, `--resume` keeps the run's
   template. The runner judges with checks.py of `--checks-dir`, the run's snapshot (it has no stub fallback).
7. `--max-requests N`: the direct arms' requests of this invocation in total, the closing capabilities GETs
   included (one kept back per direct arm in paired mode; sequential: direct-remote gets what direct-local
   left); the release GET is one more. `--smoke` defaults it to 40 and the manifest to
   manifest/suite_smoke_v5.json; without `--manifest` a full run takes manifest/suite_v5.json (the frozen v5, below;
   v1-v4 stay for reproducing earlier runs).
   The MCP runner paces its own requests.
8. Samples (review M5): the mex samplers write ~/suite-dev/runs/<id>/samples/; every invocation pulls the
   whole directory on exit (also after Ctrl-C) and at its start, merge_run deduplicates.
9. Ctrl-C: the step in flight stops after its request (`pkill -INT` of this run's runner on mex; a local
   step gets the terminal's SIGINT), the orchestrator waits up to SUITE_WAIT_RUNNER_S (150 s) for the mex runner (a long request may
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
  evaluated with the same rule, outcome and detail (compare.relation_same_answer, a stdlib copy that
  tests/test_compare_report.py `test_same_answer_built_in_is_the_reference` holds to checks.py's on every edge
  case) and the other types read `not_evaluated`; an outcome outside the vocabulary
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
  400 % or the load1 mean exceeds 48 (half of mex's 96 cores). compare.json is computed when missing, by the
  report and the gate as compare.py computes it (compare.load_or_compute: the table validated against the run
  manifest's scenarios; a table compare.py refuses makes report.py and doc/make_doc.py exit 4 as compare.py does
  ("make_doc: refusing: ...", L5 of the verification of SD-017), and is the gate's "compare.py refused" problem,
  never a cached comparison; SD-017).

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
  limitation that cannot have happened ("Stated reasons") explains no difference: fail. A /traverse seed is
  time-limited by a time or attempt_deadline resource stop, a time_budget cap trigger, the attempt's bound
  (AttemptAtBound), or any limitation of `bounds.time_budget_ms` on the result (a derivation the budget ended,
  failed or partial) or on an arm (a time cap that cut walks after another trigger); scripts/traversal's
  bench_traverse.time_bound reads a seed alike, and tests/test_checks.py `TimeLimitedSeeds` holds the two
  together (SD-036).
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
  of values must be equal (without timing; a path is read on the body, each step an object key or a list index, a
  negative one counting from the end, `*` no step); a missing path fails (not_exposed on the MCP arm, as a 200
  without a body). Without `paths`, the two bodies must be equal without timing and without the dotted paths in
  `args.ignore` (default none; an ignore path names object keys and looks through a list on its way, so
  `patterns.stop` drops `stop` in every pattern; the builder refuses an ignore path with an index, `*` or empty
  step, which would match nothing, and no manifest uses `args.ignore`). Values are compared as JSON values, as
  suite_common.json_equal: numbers by value (1 equals 1.0), a boolean only the same boolean (true is not 1, false
  not 0), in checks.py and compare.py's built-in copy together (the owner's decisions on SD-011; the stdlib copy and
  its parity test stay). Used for: max_mismatches 0 versus
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
  finalize_reserve_ms - 2 E (the reserve the answer states: the host's is configurable, --pattern-finalize-ms, and
  checks.py has no host document to check it against; else 250; see "Open items"), E the §7.6 formula at the default (slowest)
  delivery rates 10 / 50 MB/s over the
  results in the answer plus every withheld result at a generous 2 KiB (64 KiB with labels or a predicate);
  `max_steps` the steps of all patterns at max_steps (every fixture's stop shows them equal); `max_contexts`,
  `max_anchors`, `max_paths`, `max_predicate_contexts` stop_at_threshold; `max_annotation_work` a search that reads
  annotation (labels all, supported_paths, a predicate); `max_predicate_work` a predicate; `max_memory` taken as
  stated; any other reason fails. A stop of the graph search (discovery, mask_scan, extraction, extension,
  selection) with neither a cut nor a withheld list is an invariant failure (§7.5: only these shorten a list), and
  a list with neither holds every result: a source missing from it fails (retrieval_complete false alone means
  the labels are incomplete).
- /resolve, `stop` (SPEC-labeled-traversal-core §4.5): reason `time` only, with `timing.elapsed_ms` at or after the
  budget less the reserve. The reserve is the host's 250 ms (not configurable), never the one the answer states:
  `limits`, stated with every budget, must state finalize_reserve_ms 250, so invariants fails a /resolve answer
  whose `limits` states another reserve (100 as well as 900), a non-number or none, on every arm, and a time stop
  stating one fails the checks that judge the stop (L1 of the verification of SD-021). Every budget check reads the time budget
  the same way (`_effective_budget`): the answer's echo when it is a number, else the request's, so an echoed null
  never skips a time check (SD-021, SD-286, SD-296).
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
  min(the requested budget, the route's cap the builder states as cap_ms): before it a fail (lane S7c: on /resolve
  too, not before it less the 250 ms reserve); long past it a finding). A scenario that expects a refusal (no 200),
  names the deadline itself, or another route or code is judged as it states.

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
  SPEC-pattern §19.3 other than a name fault (`_predicate_form_fault`; a null, boolean, list or object as a name
  names no case: see below), and pattern.predicate_name_not_string for a number in a name list only (not where a
  predicate or a list belongs). tests/test_manifest_lib.py holds both to the shared corpus
  metagraph/tests/data/pattern_predicate_forms.json: the host's verdicts of tests/cli/test_pattern_predicate.cpp
  and the service's case per predicate, 'unknown' where the service's own reading is not stated (SD-031; the C++
  test TheSharedCorpusIsTheHostsVerdict reads the corpus on gr/sd-server). The service owner's statement of 2026-10-17
  (its app/pattern/validate.py) reads the rows the builder had only inferred and those two (`service_reading`,
  `service_source`): `{"and": [5]}` and n 2.0 are refused with pattern.predicate_form_other's text; a null, boolean or
  list name with "a name is a string, the column label as the index stores it" (no taxid hint), and a present null
  predicate is answered as no predicate, answers no case of the table states (`table_gap`, service_case null: the
  predictor names no case for them); a predicate that is no object fails the service's pydantic type, not its form
  check, so those two rows are `unknown` again.
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
  503) fails one whose client wall time is below the budget less 1 ms, on /resolve too (lane S7c: the reserve ends
  the work, not the answer; SPEC-labeled-traversal-core §4.5). A /pattern scenario where the deadline can be
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
  labels the search stays sequence-only, and the check's detail says so (`[sequence only: ...]`), except where the
  request names the record's label (its header or its taxid) as a seed label (`seeds[i].labels`; on the MCP arm
  the MCP arguments' seeds): a label dictionary (`label_dict`, a graphlet's L records, also none) that lacks it FAILS
  walk_in_record and walk_end on every arm, naming the label and the record ("seed 0: the request names the
  record's label NC_X.1 (record NC_X.1), and the answer's label dictionary (0 label(s)) lacks it"; the owner's
  decision on L8 of the re-verification: a well-formed body with no L records and num_seed_labels 0 passed on the
  sequence alone). A result that states no dictionary at all keeps the rule before. On the mini-index
  fixtures and on the 714 recorded MCP graphlet segments of the full run whose direct-local twin has the same shape,
  the derived sets equal the full result's labels. The MGT tokens are read as strictly as SPEC §7.5.1 and the
  library (_codec) read them, and tests/test_checks.py `MgtCodecVectors` holds checks.py's readers to the golden
  vectors (api/python/tests/data/traverse/codec_vectors.tsv): RANGES ascending and collapsed with one spelling per
  id (no `01`, `+1`, `0-0`, `0,1`), every id 64-bit and below the number of L records (none in a body without L
  records); a SETEXPR delta with no empty part and whose R lies in its base and A outside it; the free-text escapes
  %25, %0A and %0D only; an L record's front-coding prefix an integer within the previous name, ending on a
  character boundary; an integer field (entry_total, an L record's seq_id) with no sign or leading zero; the S
  record's num_seed_labels at most the number of L records (constrain mode, as the library's parser), and a `*`
  entry only where §7.5.3 derives it (the root, a merge's partition, constrain mode). A body that breaks one of
  these FAILS walk_in_record and walk_end on every arm (the owner's decision on SD-014): the detail names the
  record, the token and the rule (`the graphlet body breaks SPEC-labeled-traversal-core §7.5: the G record of left
  segment 0, end '00': RANGES item '00': an id is a decimal integer without sign or leading zero, ...`), and the
  walk is never judged on the sequence alone for it (the silent sequence-only fallback before). The other checks
  read such a body's structure as before. Every recorded body of the v5 and v6 baselines reads.
- The route search is memoised (A-05) on (segment, the seed copies still followed, whether the stated event was seen):
  every route to a segment has the same length and the copies it still follows fix its flank, so each state is
  searched once; `_MAX_ROUTE_STATES` (100,000 states) is a safeguard. A route counted "beyond the truth flank" now
  follows the flank to its end before it runs past it.
- A walk blocked by its own seed or edges (G28, VERIFY-2 item 2: the strict rule; the triage's loose rule accepted a
  relabelled reason, a continuation off the record and a fabricated block). An arm with no route that follows the
  record over its whole length still follows it when (1) a route spells the record for exactly its first d bases,
  where (2) d is the first depth at which the record's own continuation is inadmissible by SPEC-core §6.6,
  computed from the truth's flank and the request's seed (`_first_block`, which decides each step as the oracle's
  Walker.block and the engine do, tests/test_oracle_traverse.py FirstBlockParity: in a canonical/primary regime a
  hairpin step under `hairpins: skip` first (SPEC §6.5: a palindromic (k+1)-mer); rejoined_seed, a seed k-mer,
  either orientation in a canonical/primary regime; else edge_reuse_rc / edge_reuse, the step's (k+1)-mer or its
  reverse complement already used on the walk, a palindromic one reused being edge_reuse), (3) one of its segments
  states the blocked event at at_bp d with the
  record's base and that reason, and (4) past d it spells the record from a collapsed repeat's exit (an occurrence
  of the junction's k-1 bases in the flank followed by another base) for every base the flank holds. A block before
  the truth's min_bp is a truth gap (inconclusive), so is a continuation past the flank; past the record's end with
  the record as the only label it fails. A following route that ends at the computed block before the truth's
  min_bp, with that reason, is the same truth gap; a skipped hairpin only where the answer states it (L3 of the
  verification of SD-025): the path's `dead_end` with the label's `hairpin` end (SPEC §13 note 2; a graphlet's
  `Dh`) or the skipped `hairpin` event at that depth with the record's base on the route's leaf (a graphlet's
  `E <at_bp> h <char> ... s` record, which `graphlet_arms` now reads); a bare `dead_end` there is a walk that
  stopped short. The hairpin policy is the arm's (`_strategy_knob` with the arm: on the MCP arm its overrides,
  then the direct request, then the default skip). The regime is the answer's capabilities or the graphlet's H
  record (basic on this index, where no step is a hairpin and the change of SD-025 decides nothing).


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
  wall >= the request's time_budget_ms (the MCP arguments', else the request's) - 1,000 ms (the host answers its
  deadline only once the budget has elapsed; a reserve its text states above the host's holds the case back); a
  known wall below holds the case back (matching continues, the checks judge).
  resolve_timeout (R2): wall >= max(budget + reserve + 10,000 ms, the stated "within N s") - 1,000 ms;
  traversal_resolve_deadline (G22): wall >= the budget the service sent - 1,000 ms, the same floor
  (suite_common.deadline_least; SD-020: it was budget - reserve, the premise 84e0c7a98 refuted for the direct
  arm; tests/test_common.py DeadlineFloor holds the two rules to one floor; async's runner keeps its own copy of
  G22's budget - reserve until it follows, so for a wall in [B - 1000, B - 250) its recorded fail and
  reevaluate's inconclusive differ, both blocking); the R20
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
  consumed, unless it expects the slot refusal (mismatches_above_cap) or a refusal of the request (a status that
  admits no 200: `[200, 400]` admits the answer); D picks come from what is left. The builder side reads a
  scenario's allowed statuses through one reader, `manifest_lib.allowed_statuses` (build_manifest, areas,
  budget_expectation), and build_manifest.validate refuses a scenario with more than one status expectation or an
  expectation whose args is no object, so the judging side's readers (the last entry, every entry) read the same
  list (SD-094; checks.py and compare.py keep their own loops: the checks are pinned per run).
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
  service_by_design_case is kept, as the runner's mark is; the comparison reads it alike (SD-022,
  compare.runner_by_design_finding, which reevaluate calls too): an unverified by-design part of the runner's, as its
  service_by_design mark, where it compared the line check by check (a disagreement).

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
  (host_state_inconclusive) read a line recorded without the runner's inconclusive_reason alike
  (suite_common.host_state_of_record: the raw record's error and tool error, the run's References,
  suite_common.run_references: its recorded live documents over the table's blocks; SD-035: the merge read the
  table alone and the line's own error); on the reason path they ask different questions on purpose: the merge
  takes any host-state reason or rule (c) as the runner's decision (an answer-triggered case's undecidable wall
  included), the comparison counts a host-state repeat only with a host-state code (B5b-09).
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
  run_suite.sh --ids or async's runner standalone (`gate.py policy-ids POLICY` prints them comma-separated). A
  hand-typed list with spaces selects as the clean one: suite_common.select_scenarios strips ids as it strips
  areas, and build_manifest --areas 'P, T' builds both areas, empties left out (SD-369; it dropped T without a word).
- **The gate** (`gate.py RUN_DIR POLICY [--reevaluated]`: exit 0 PASS, 1 FAIL, 2 usage; RUN_DIR/gate.json) fails,
  listing every problem by kind, on: identity (either manifest hash not the policy's, a stated hash not verified, an
  arm or a plan on another manifest, another index_fp, another by-design table); missing (a planned or required
  execution not recorded, or a required one the run did not plan); outcome (on the policy's arms and selection: fail,
  error, finding, ANY inconclusive, not_implemented_mcp on a scenario not excluded, service_by_design by a case the
  scenario does not name in the manifest (suite_common.named_case_ids: mcp.by_design, service_by_design, sbd:<id>
  tags, as compare and the MCP runner read them; SD-030; a table case's scenario_ids, which compare.named_cases and
  the runner also read, is not a manifest naming and the gate does not accept it, by the owner's decision: compare
  flags such a match (mcp_agreement.named_by_table_only, summary.mcp_agreement.named_by_table_only_scenarios) and
  the report lists it, "By design only through the table's scenario_ids"; build_manifest.validate refuses an
  mcp.by_design that is no list), a by-design
  pair not verified when a direct arm ran it, the MCP side not verified when none did); property (a property
  check the scenario expects not passed on a passing or explained execution, e.g. not exposed; the line's outcome
  is not trusted: every check the scenario expects must have a result on each passing line, determinism excepted
  on the first repeat only, the status must pass and no other generic check may fail or
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
  list, and they match it with one rule, suite_common.override_path_listed (SD-095: `x.*` lists the bare section `x`
  too, as the table's text and the service read `{"bounds": 5}`; the gate called a bare section advanced). The
  module's own copy of the service's list (MCP_BASIC_OVERRIDES) is gone, and the table's case is required: a
  scenario's overrides with a table that states no list (no table, `--by-design ''`, included) are refused
  (gate_policy.NoBasicOverrides; gate_policy.py exits 2 naming --by-design); the policies pin the real table, which
  states it. gate/policy_v4_rt_mcp.json is regenerated after merging lanes S6g1 and S6b: the same
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
  <= end <= nt_length; source_placement applies the same rule to its record's occurrence; a label's
  occurrence_list null is "not placed", SPEC-pattern-search §14.5, which source_placement fails on every arm unless
  another source context places the record, while an absent list is a field the answer lacks, not exposed on the
  MCP arm: SD-377). Closed: L-verified-0001
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
  is a mismatch, judged with the live value. An MCP line records `request_evidence` {case: {reference: live |
  table | unknown}} for every case whose request was evaluated and names at least one {"host"|"service"} reference
  (section "Suite integration S7i"; left out when empty), and `request_evidence_mismatches`; compare.json's
  `references` and a warning line in the report
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
  /resolve deadline). The relayed deadline's threshold is the budget sent - 1 s (deadline_least, as the pattern
  wall rule's). A legacy record (no stated
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

## Lane S7r: the runner and the orchestration (the two-day review's S-runner lane, 2026-10-14)
The findings of the two-day review (2026-10-14-ultracode REPORT.md) on run_direct.py, run_suite.sh and merge_run.py,
decided against 14552938 (lanes G1-G3 merged): the resume selection's widening and the narrowed verdict were fixed
by lane S6c; everything below was open.

### An unanswered repeat is no earlier repeat (P3-run-measure-report-01, P6-end-to-end-06)
- A request that got no HTTP status (a transport failure, a client timeout) is recorded with the digest of its
  stand-in answer (`suite_common.transport_failure`). run_direct passed it in `previous`, in the session and on a
  resume, so the next repeat, which the host answered, failed `deterministic` ("differs from 1 earlier repeat(s)"):
  a false determinism failure charged to the host, and a direct or MCP disagreement. `suite_common.feeds_previous`
  is the rule: an earlier repeat feeds `previous` only with an HTTP status and an outcome other than `error` (an
  undecodable body or a runner failure has no digest of the host's answer either). A repeat after an unanswered one
  compares with the answered ones only; with none, `deterministic` has nothing to compare (left out, as on the first
  repeat). The rule is async's BENCH-01 for its runner and is reevaluate.py's to apply as well (a mirror rule: the
  same helper).

### Ctrl-C of a paired run over ssh (P3-run-measure-report-02)
- A paired direct-local step's stdout is its ssh channel. Ctrl-C ends the ssh client; the runner on mex survives on
  its request, gets the orchestrator's `pkill -INT`, records its line, and its next stdout write (the progress line,
  or with --quiet the status line) met the closed pipe: BrokenPipeError, state `failed`, a "runner failed:
  Traceback" warning that stayed in the arm's run.json (and the report) after a successful resume, exit 120 instead of
  130. Every stdout write of run_direct goes through `say`: the first write that finds the reader gone points stdout
  at /dev/null, and the run ends as its stop says (interrupted, exit 130, no warning). Tests: run_direct with its
  stdout reader closed before SIGINT (with and without --quiet), and Ctrl-C of a paired run over tests/fake_ssh.py
  followed by a resume.

### A failed step is recorded (P3-run-measure-report-05)
- A direct step of the paired loop that exited non-zero only counted towards the stop after three in a row: its
  scenario still counted as completed and the loop recorded `done N/N` (the v4 log: `python3: No such file or
  directory`, exit 127 on P-iupac-0001, then "paired: done, 681/681"). The loop record (orchestrator.json
  sessions[-1].loop, `loop_record`) counts a scenario as completed only when every step of it succeeded, lists
  `failed_steps` as ARM:SCENARIO:EXIT, and a loop that ran through with a failed step ends `incomplete` (reason
  `failed_steps`; the report's verdict reads INCOMPLETE through it) with exit status 1, like a failed loop. No
  retry: a step that failed is resumed like any other (`--resume`). A loop stopped by the request cap stays `capped`
  (exit 0, a planned stop that --resume continues) also when a step failed before the cap; its failed_steps are
  recorded and the report reads INCOMPLETE.

### The orchestration records what its arms lost (P6-end-to-end-01, the run_suite side)
- Lane S6c's plan reconciliation already reads a planned execution no arm recorded as INCOMPLETE (the report side
  of the finding: the sequential mode with the MCP cold pass refused, exit 3, reads INCOMPLETE at 14552938). The
  orchestration itself did not say so: the sequential mode wrote no loop record and exited 0 whatever its arms did,
  and the MCP runner's output went to the terminal only, so a runner crash (a traceback, no line) left nothing in
  the run directory.
- The sequential mode writes the loop record too (mode `sequential`, no per-scenario counts): `failed` with the
  MCP rule's reason when the MCP runner refused (exit 2 / 3), was interrupted (130) or ended three executions in a
  row in error or unrecorded; `capped` when an arm stopped at --max-requests; `incomplete` (reason failed_steps)
  when direct-local did not finish (`direct-local:-:STATE`) or direct-remote exited otherwise
  (`direct-remote:-:EXIT`); else `done`. failed and incomplete exit 1 in both modes.
- The MCP runner's output (stdout and stderr) goes to the terminal and to arms/mcp/runner.log, each step framed by
  its repeat and ids and its exit status (the template's expansion, which may name files, is not written; `tee -i`:
  on Ctrl-C the runner's last lines still reach the log, where a tee killed by the same Ctrl-C would leave the runner
  a closed pipe, a BrokenPipeError and exit 1 instead of 130); the
  direct-remote runner's stderr goes to arms/direct-remote/runner.log (its last lines shown when the step fails), as
  direct-local's on mex.
- mcp_after names every MCP execution a runner step did not record (`mcp_missing`, SCENARIO:rREPEAT, the first 50;
  `mcp_missing_count`); in the paired loop such a scenario is not completed, and a loop that ran through with one
  ends `incomplete` (reason mcp_missing) instead of `done`. Tests: tests/fake_mcp.py `crash` (a traceback on stderr,
  no line; FAKE_MCP_CRASH_ID limits it to one scenario) in the paired loop, and the sequential mode with the cold
  pass refused.

### A resume waits for its own run's runner (P3-run-measure-report-06)
- Ctrl-C of a paired run leaves the mex runner of the step in flight on its request (a long-class request may take
  up to 30 min); the interrupt message promised that "--resume waits for it or attaches to it", but only the
  sequential mode attaches, and a paired resume refused after SUITE_WAIT_RUNNER_S (150 s) with "... or --resume its
  run", circular while that runner was its own. A resume now first waits for its own run's runners (the mex runners
  whose command line names its --run-id, when no other run's runner is there) for up to SUITE_WAIT_OWN_RUNNER_S (1860
  s: the long class's client timeout and a minute), with a log line every minute, then goes on; only then the
  check for any runner refuses, naming the pids (another run's runner: wait for it or stop it with kill -INT, which
  records its request in flight; its own past the wait: stop it with kill -INT, then --resume again). The interrupt
  message is the mode's: the sequential mode's resume attaches, the paired mode's waits.

### A resume continues the run as it was started (P3-run-measure-report-04, P6-end-to-end-05)
- Lane S6c fixed the selection: a plain resume runs the union of every invocation's plan (it no longer widens to the
  whole manifest), and a narrowed resume's report reads INCOMPLETE against that union (the probe of 14552938: a
  `--limit 1` run resumed plainly adds no scenario; a capped run resumed with `--areas P --limit 1` reads
  `INCOMPLETE (... 12 planned execution(s) not recorded ...)`). The other options still came from the resume's
  command line: a `--skip-local` run resumed plainly ran direct-local too (the arms widened), a smoke run lost its
  request cap of 40 and its 1800 s sampler bound, the pause fell back to 1 s, --allow-stub had to be given again,
  and a gated run resumed plainly was not gated (it exited 0 on a FAIL).
- A resume keeps the options it does not give: the last invocation's mode, --checks, MCP template, --warm-every (as
  before) and --pause, and the run's FIRST invocation's flags --skip-local, --skip-remote, --allow-stub, --smoke and
  --gate (the policy must still exist). A resume that names a flag, or its negation (--local, --remote,
  --no-allow-stub, --no-gate; a flag and its negation exclude each other), applies it to that invocation only and
  records it in its session: one `--resume --skip-local` does not drop direct-local from the later resumes (taken
  from the last invocation, it did, and the run could never finish), and `--resume --no-gate` leaves one invocation
  ungated. --max-requests is not kept: it caps one invocation (a capped run is resumed to finish it); a smoke run's
  default cap of 40 applies to each invocation, its resumes' included (a larger --max-requests lifts it). The resume
  logs "the run's own options: ..." and each session of orchestrator.json records `given` (the options on its
  command line) and `restored` (those taken from the run).
- A resume that names its own selection (--ids, --areas, --limit, --gate) and leaves out scenarios of the run's
  plan says so ("this invocation's selection leaves out scenarios of the run's plan ...") and records
  `narrower_than_plan` {selected, planned, left_out} in its session; its loop and its arms' states describe its own
  selection, the report's verdict the run's plan (INCOMPLETE while any planned execution is missing).
- The tests no longer pass the run's own options again on resume (--warm-every, --skip-local, --skip-remote,
  --no-mcp), so they exercise the restore.

### The release at the run's end (P6-end-to-end-07)
- run_direct kept only `server_instance` from its closing GET /capabilities, and a resume set `run.release` only
  where it was absent: a deploy during the run, or between an interrupt and its resume, read at most as "the server
  restarted", and every execution was attributed to the first release. Each session of an arm's run.json records
  `release_start` / `index_fp_start` (the start's capabilities documents) and `release_end` / `index_fp_end` (the
  closing GET /capabilities; staging's /capabilities states no index_fp, so index_fp_end is null there and the
  release, which names the index build and the commit, carries the identity); run.json `releases` lists every
  release the arm saw in order, `release_changed` whether there was more than one. `run.release` stays the run's
  first.
- merge_run.release_identity reads both direct arms' sessions in time order: more than one release or index_fp is a
  warning naming the first and each later one with the arm, session and time, and run.json `release_identity` /
  `mixed_releases` state it. merge_run --strict-identity (run_suite.sh --gate) exits 3 on it (identity_errors) and
  gate.py fails on it (kind identity, one line in gate.py). The report shows the warning; its verdict is report.py's
  (lane S-verdict's file): it does not yet read `mixed_releases`.
- A paired step (--only) takes no capabilities GET, so the release is seen at the arm's first step and at each
  invocation's closing step: a deploy during an invocation shows at its end, one between an interrupt and the resume
  at the resume's end. In paired mode the first invocation of a run interrupted before its closing step has no end
  release; its resume's closing step has one.

### The checks are pinned per run (P6-end-to-end-03)
- Each paired direct-remote step was a fresh process importing checks.py and suite_common.py from this directory,
  the MCP runner re-read `--checks-dir` (this directory) on every invocation, direct-local used the copy rsynced to
  mex once at the start, run.json recorded only `checks_module: "checks"`, and merge_run compared manifests only: a
  checks change landing during a five-hour run (7161cbcc was committed inside the v4 window) would have judged the
  arms by different rules, silently.
- **The run's snapshot.** The run's first invocation copies what runs and judges the arms
  (`suite_common.SNAPSHOT_FILES`: run_direct.py, sampler.py, suite_common.py, checks.py, checks_stub.py,
  service_by_design.json) into RUN_DIR/scripts, with RUN_DIR/scripts/SNAPSHOT.json (created_utc, source,
  scripts_commit, sessions_before, the sha256 of each file). direct-remote runs `RUN_DIR/scripts/run_direct.py`, mex
  gets the snapshot (rsync of RUN_DIR/scripts instead of this whole directory), and the MCP template's new
  placeholder `{checks_dir}` is the snapshot (the default template passes `--checks-dir {checks_dir}`; a restored
  template of before with `--checks-dir` this directory is rewritten to it; a template without `{checks_dir}` is
  warned about at the start). A resume keeps the snapshot and says which files of this directory differ from it; the
  first resume of a run older than the snapshot takes one then (`sessions_before` the earlier invocations). Each
  session of orchestrator.json records `snapshot` {state created / kept, changed, live_differs, mcp_checks_dir}.
- **What judged each arm.** run_direct records per session `checks_sha256`, `suite_common_sha256` and `checks_dir`
  (of the modules it imported; the stub's file under --checks stub), and the latest at run.json's top. async's runner
  records per invocation (run.json `invocations[]`; it keeps `checks_module` out of them) `checks_dir` and, since
  its 1dccc344, `checks_sha256` and `suite_common_sha256` (the sha256 of checks.py and suite_common.py it imported);
  merge_run reads exactly these per-invocation keys (a runner without `invocations`, as tests/fake_mcp.py, states its
  one at run.json's top) and never `checks_module`.
- **merge_run.checks_identity** (recorded sets; a re-evaluated set is judged with one checks module, the
  re-evaluation's) compares them with the snapshot: a direct session whose checks or suite_common sha256 is not the
  snapshot's (the stub's for a stub arm), an MCP invocation whose recorded checks directory is not RUN_ID/scripts (a
  copy of the run elsewhere keeps that tail) or whose stated checks or suite_common sha256 is not the snapshot's, and
  a snapshot file changed after it was taken are `differs`; a session or invocation stating neither is `unverified`
  (sessions that started before the snapshot are counted apart); one sentence per arm and kind, since a paired run
  has a session per step. Both are warnings and run.json `checks_identity`; with a snapshot both are errors for
  --strict-identity (exit 3) and for gate.py (kind identity). Without a snapshot (a run before lane S7r) only arms
  stating different sha256 are named. The outcomes merge_run makes stricter with checks.scenario_outcome are the
  snapshot's checks' for a recorded set; compare, the gate and the report judge with them too (section "Suite
  integration S7i").
- **A gated resume of a run started before the snapshot existed always fails the checks identity.** Its first
  resume takes the snapshot, and the sessions of its earlier invocations cannot state a checks sha256: they are
  `unverified`, an error with a snapshot. This is the strict gate's intent (what judged those executions is not
  recorded); the remedy is a fresh run. (A re-evaluated set is judged with one checks module, the re-evaluation's,
  and the check does not apply to it.)
- **`--mcp-run DIR`** (a finished MCP run copied into arms/mcp) was not judged with the run's snapshot: its
  invocations state another checks directory (or none), which the merge names and the gate fails on.
- compare.py, report.py and gate.py's property reading still load this directory's checks.py and table for a
  recorded set (lane S-verdict's files; noted for the integrator).

### An outcome outside the vocabulary is never a success on the runner side (P6-end-to-end-02, hand-off of lane S7v)
- merge_run.stricter_outcomes skipped a line whose outcome is outside the shared vocabulary (it kept it as
  recorded), and the paired loop's mcp_after reset its error counter on one: six MCP lines with the outcome
  `skipped_by_runner` merged as such, and three in a row did not stop the loop. `suite_common.outcome_known` is the
  shared test (defined identically to lane S7v's, gr/suite-s-verdict 9a12a6c9). merge_run reads such a line (or one
  without an outcome) as `error`, `outcome_runner` keeping the recorded one, with a warning per arm and run.json
  `outcomes_unknown` {arm: {recorded outcome: lines}}; mcp_step_result reports it as `unknown`, which mcp_after
  counts like `error` towards the stop after three in a row and logs. (The comparator, the report and the gate read
  the recorded outcome back from outcome_runner: section "Suite integration S7i".)

### Validation (offline)
- Copies of v3 (the S5e input set) and v4 re-evaluated (`run_suite.sh --reevaluate COPY --no-pdf --reeval-mcp-auth
  anonymous`) with 14552938 and with this lane: no outcome, check outcome, direct or MCP class, verdict line or
  merged warning moves (v3 FAIL, v4 FINDINGS, as before); both merged sets state one release and no checks identity
  error (a re-evaluated set is judged with the re-evaluation's checks). The gated re-evaluation of v4 with
  gate/policy_v4_rt_mcp.json reads the same 45 problems (GATE FAIL) before and after. The recorded v4 merge (not
  re-evaluated) gives the same 2913 outcomes and host-cache labels; its checks identity has no snapshot and states
  the sessions without a sha256 as unverified, which is no warning and no error for a run older than the snapshot.
  The re-evaluation does not run run_direct: reevaluate.py applies suite_common.feeds_previous since the
  integration (section "Suite integration S7i").
- 371 tests pass (1 skipped; 345 at 14552938), the follow-up of the lane's verification included (first-session
  flags and their negations, `tee -i`, the per-invocation sha256 keys); its re-evaluations of v3 / v4 and the gated
  v4 again move nothing.

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
  direct-remote executions that gave an item are in; without a fallback it is unchanged). No item falls back to an
  execution without an answer (outcome error, or no HTTP status: a transport failure fails its items): the
  verification of the lane (VERIFY-S-VERDICT (3)) found that an MCP pass beside a direct-local 503 and a
  direct-remote transport failure read disagree, and an MCP runner error in the same place agree. The report's
  agreement depth reads the same items. The by-design pairs were already judged on every repeat of both direct arms (lane G2).
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
  reaches the host now, except a resolve_timeout in the write phase (the line's timeout_terms, else its error's
  timeout_phase: the service may never have delivered the request; the recorded ones are all read timeouts). `backend_timed` keeps the old test where a backend timing is what matters: the MCP overheads
  (a relayed deadline's wall is the deadline's) and whether a line records the service's answer.
- assign_host_cache adds `host_cache_after_interrupted: true` (present only then) to an execution whose previous host
  execution of the scenario ended at the host's deadline without its answer (`interrupted_on_host`: a direct 503
  deadline or a relayed one). The timing statistics key on `compare.timing_cache`: such a warm execution reads
  `warm_after_interrupted` ("warm after a cut" in the timing table), counted apart from the warm ones, so the front
  overhead's warm pairs, the MCP overheads and the cold vs warm pairs leave it out (v4: the 10 direct-local ones
  have a median server time 2.48x the other warm ones'); compare's summary counts
  `warm_pairs_left_out_after_interrupted` (scenarios whose warm pair it would have been) apart from the unpaired, and
  the report's Overheads section states the count (VERIFY-S-VERDICT (4)).
- Moves (timing only; no outcome, check, class or verdict): v4 39 MCP lines not_reached -> cold (27) / warm (12), the
  27 direct-local r0 after them cold -> warm with since_prev_s; 70 lines gain host_cache_after_interrupted; 32
  scenarios are left out of the front overhead's warm pairs after a cut (wall pairs 552 -> 546, server-time pairs
  527 -> 527). v3: 10 MCP lines, 8 direct-local r0 cold -> warm; 18 scenarios left out after a cut (wall pairs
  562 -> 551, server-time pairs 530 -> 530).

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
  judges such a mark. With a table, a mark on a scenario the manifest maps to a tool or a flow (`compare.mcp_mapped`:
  mcp.tool, or mcp.flow, the runner's area S flows) that
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
  and a list name them. gate.py already failed such a line ("outcome '...'"). merge_run and the paired loop's
  mcp_after read it as error (lane S7r); the comparison states the recorded one (section "Suite integration S7i").
- Recorded runs: every outcome is in the vocabulary; the totals gain `other: 0` and the summary an empty list.

### The load caveat without a measurement
- load_caveat read a missing sra_all figure as 0 and printed "mex was not under the Logan load" when every sample
  had failed (docker stats: a missing or renamed container, the daemon down; the sampler's `error`, its figures
  null), with load1 moderate or null. load_figures counts the samples stating each figure and the failed ones (the
  first error); with no sra_all figure, or more than 10 % failed samples (LOAD_SAMPLE_ERROR_SHARE), the verdict
  section says "**Load not measured:** ..." with the figures and the failures, never the positive sentence; a busy
  host keeps the Load caveat (with the failures named); the resource section names the failed samples. The recorded
  runs have no failed sample (v3 12,756, v4 8,693, all with sra_all): their sentence is unchanged. A NaN or inf
  figure is no reading either (suite_common.finite_or_none, the bench load guard's rule of c6306ec1a; the sampler
  records one as none): it never makes a mean NaN nor reads as a clean "not under the Logan load" (SD-018).

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
    An MCP execution that is inconclusive only because checks read not_exposed that the scenario's NAMED by-design
    case explains (`report.stated_mcp_limit`: the case it matched, among the cases the scenario names, whose
    explains_checks state every check item of the line that did not pass, each not_exposed, by
    suite_common.states_limit: an entry {check, detail_matches} only an item whose detail it matches, as the verdict
    pipeline reads it, a bare name every item of its check, the report's exception to item_explained (SD-023); with
    the table compare judged the run with) is a limit the table states: it does not make a group untested, and the
    verdict section lists it
    ("Inconclusive via MCP by a limit the by-design table states for the scenario (not counted as untested): ...").
    v3 and v4: C-caps-0001 (capabilities.traverse_capabilities_document) and C-caps-0003
    (capabilities.server_capabilities_document), both explaining capabilities_doc, which the MCP answer never
    carries; so no untested group remains and a clean exploratory run with the MCP arm can read PASS. A field the MCP
    drops elsewhere (no named case explaining it) still reads untested. The verification of the lane (VERIFY-S-VERDICT
    (2)) found the first form of the rule counting these by-design limits as untested; the gate is unchanged (its
    declared file excludes C-caps-0001 on the MCP arm).
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
- Class moves, each the point of its rule: v4 N-usecase-0008, R-pos-0027 and v3 M-planted-0017 MCP not_comparable
  -> agree; v4 T-repeat-0004 and v3 T-walk-0004 one more item compared; six v3 scenarios already disagreeing by an
  MCP runner_error list their differing items. No direct class moves: v4 R-budget-0002 and R-twin-0002 stay
  identical and record execution_differs (budget_respected inconclusive on direct-local, a finding on
  direct-remote). Summary: v4 MCP 534 agree / 4 not comparable -> 536 / 2, no direct disagreement; v3 MCP 542 / 5
  -> 543 / 4. No relation moves.
- Verdicts: v4 FINDINGS (4) and v3 FAIL, both unchanged. Neither gains an untested group: both state C-caps-0001 and
  C-caps-0003 on mcp as limits the by-design table states (capabilities.traverse_capabilities_document,
  capabilities.server_capabilities_document). The deploy gate on the v4 copy (`--gate gate/policy_v4_rt_mcp.json`,
  two arms mcp + direct-local, the policy regenerated for the table's sha256): 45 problems before and after, the
  same list.

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
- A service document's identity, ONE definition both sides compute (agreed with async on 2026-10-15: the runner
  calls the vendored `suite_common.service_document_sha256` and keeps no view of its own): sha256 of
  canonical_json(`service_document_view`(document, host_document, database)) = {service_settings, service_maxima,
  presets, each where the document states it; backends: {<the record's database>: {budgets, resolve,
  effective_maxima, limits}, each where stated}; host_document where recorded}. A key the document lacks is absent,
  never null; another database's entry is no part of it; a record without its database keeps no backend entry; a
  read of the host document alone (the traverse read failed) has that document's identity; an empty view (nothing
  read) has none. Probe times, cached flags, messages and counters elsewhere in the document change between the
  invocations of one deployment and would trip the gate's "spans more than one deployment". gate.py's
  `_document_sha` (edited for this function only) prefers the recorded sha256, else this; compare's
  recorded_documents states the same. On every shape the verification listed (no service_settings, a second
  database, a backend entry lacking a key, a failed pattern or traverse read) it equals async's 26fd3e7c view over
  the same host_document.
- The host's identity (`suite_common.HOST_IDENTITY_KEYS`: release, index_fp, index_meta_fp, agreed with async on
  2026-10-15): host_document = {pattern: the /pattern/capabilities/<db> capabilities.raw, resolve: {time_budget:
  backends.<db>.budgets.resolve_time_budget}, release, index_fp, index_meta_fp: from the traverse document's
  backends.<db>.identity}, each where a read states it (`suite_common.host_document_of`, for the runner too); so a
  host redeploy with unchanged limits is another deployment by the identity. A paired run's host namespace reads the
  same paths: the direct arm's /capabilities (release at its top level) with index_fp and index_meta_fp from the
  identity the arms record (their /traverse/capabilities, else the run.json's index_fp: `compare.direct_host_document`).
  A host_document without these keys (a runner before them) still resolves: the paths it states read live, the
  identity falls back to the table's block, else unknown.
- service_by_design.json, the vocabulary's answer_match.wall_rule entry only: the name budget_minus_reserve is kept,
  its threshold is B - 1 s (B the request's time_budget_ms) on the failing attempt's own wall, no finalisation
  reserve subtracted. The table's sha256 is d1e040aac81007713e2250a3dd765244e44026c8bcb03834fd484f143ccebc61 (was
  ab4100c9...); gate/policy_v4_rt_mcp.json pinned the old one and is regenerated at the integration (section
  "Suite integration S7i").
- Agreed, no change: request_evidence_mismatches as the field name; the legacy no-budget term.
- Recorded runs: v3 and v4 record no service_capabilities, so nothing they read moves; their comparison and report
  name the table's new sha256.

## Lane S7c: checks that the two-day review found too lenient or too strict (2026-10-14)

The suite findings of the review of 2026-10-14 (ultracode REPORT §3, P1-checks-01..07 and P5-v4-content-01) on
checks.py, judged against 14552938. The deploy gate is strict, so a false FAIL or a false inconclusive blocks as a
false pass misleads: each stricter rule below is measured on every recorded answer of v3 and v4 (all arms) and on the
engine's stored pattern fixtures. The runner vendors checks.py: it takes these rules by re-vendoring.

- /traverse, walk_end (the finding that any resource cap excuses a record route that ended far below it is closed in
  walk_in_record by lane S6g1): a path that states a resource end below its arm's complete_to_bp cannot have had it
  (SPEC-labeled-traversal-core §6.10), so walk_end fails "not reached: a stated stop that cannot have happened" where
  it read inconclusive (the invariants already failed that execution).
- /pattern, count_relation (SPEC-pattern-search §7.4: at_least means a stop interrupted the count): a count that is
  at_least, unknown or bounds beside a stop of any reason that can leave it so (every budget of §7.6, §14.4, §19.8,
  §20.5: max_steps, time, max_annotation_work, max_memory, max_predicate_work; every stop_at_threshold stop) cannot
  show an exact expectation or a minimum it lies below: inconclusive (the stop's own reason rules first, as before).
  A graph count (contexts, anchors, paths, supported_paths) is not excused by a stop of the labels' reads or the
  output (label_discovery, placement, output: §7.5, "the graph count is then still the discovery's"). Before, only
  max_steps and time counted, in any phase, so a supported-path search stopped by max_annotation_work or max_memory
  with supported_paths at_least 0 passed min 1, and paths at_least 3 after {extension, max_memory} failed relation
  exact. Without such a stop, a graph count at_least below the minimum fails: nothing can have left it inexact (it
  passed). The inexactness the SPEC states without a stop stays inconclusive: supported_paths after a predicate's
  pruning and the paths of a supported-path search after any branch pruned before L (§20.4, §20.9), bounds on a
  graph without its mask (§7.4); the counts with relation rules of their own (labels, occurrences, tested, selected:
  §12.1, §14.6, §19.7) read inconclusive at_least below a minimum (they passed) and keep their other rules. Bounds
  that straddle a minimum are inconclusive (they passed); an expected exact value an inexact count excludes (at_least
  above it, bounds around other values) fails also after a stop, one it does not exclude is inconclusive after a
  stop (it failed). Measured: every inexact count of v3 and v4 without a stop is a supported-path search's after a
  predicate's pruning (the L-supported scenarios, at_least 1 >= their min 1) or a selection's (bounds, pass
  stopped), or a label count; every inexact count of v3, v4 and the stored fixtures judged against min = its value
  + 1 and relation exact fails nowhere the checks of 14552938 did not (stress_count: 0 new failures; of 1,154
  such hypothetical expectations that move, 676 earlier failures and 478 passes read inconclusive). No recorded
  count_relation outcome moves; two inconclusive details (M-planted-0016 / -0017 MCP) name the relation expected.
- /resolve, resolve_discover_contains and the invariants test the discovery's ranking (SPEC-labeled-traversal-core
  §4.2: the top max_labels by supported k-mers, ties by column id, then seq_id; `labels_truncated` {kept, total,
  min_kept_kmers, max_dropped_kmers, dropped_full_length}). A truncated list without the source label passed as
  "index order, not best-supported", a misreading of the SPEC: R-discover-0001 / -0004 (limit 1, always truncated)
  could never fail. The source record's label supports every in-graph k-mer of a query taken from it (graph_runs
  summed; under a stop, the prefix's), so a truncated list may leave it out only in a tie at that support:
  min_kept_kmers and max_dropped_kmers equal to it and dropped_full_length >= 1. That answer passes (the tie's column
  ids are not in the answer; inconclusive would block the strict gate on a correct host every run, since limit 1
  with a seven-way tie always hides the source): the recorded R-discover-0001 (713600 kept, 70 of 70, 7 dropped at
  full length) passes with a detail that names the tie. A kept label below full support, a dropped label above a
  kept one or no full-length drop fails (the index-order host of the finding: 670 kept with 3 of 70, max_dropped
  70). labels_truncated without its fields is not_exposed on the MCP arm and fails on a direct one. The invariants of
  a discovery: at most max_labels labels, ranked by kmers_supported (not under support trace, whose count is of trace
  runs), none above the in-graph k-mers; a truncation states kept = the list's length < total, max_dropped_kmers <=
  min_kept_kmers = the last kept label's count, min_kept_kmers within the in-graph k-mers, and full-length drops
  (at most total - kept) only with max_dropped_kmers = the in-graph k-mers. Measured on every recorded discovery of
  v3 and v4 (R-discover, N-usecase-0010 / -0012, R-edge-0023 with 150 kept of 2,753, N-negative): none violates a
  rule; R-discover-0001's eight passes keep their outcome with the new detail. DESIGN-test-suite.md's R-discover
  line is corrected.
- /resolve, the truncation's size (the lane's verification): only more than max_labels labels are truncated, to
  exactly max_labels (SPEC-labeled-traversal-core §4.2; the request's max_labels, else the host's default 1000, with
  no server cap). A list cut to another size (max_labels 100, kept 1 of 18 at 70 of 70: R-discover-0003's shape,
  where a correct host lists all 18) passed resolve_discover_contains and the invariants as a tie; both now fail it
  (kept != max_labels, or total <= max_labels). The tie pass rests on the builder's premise that the source supports
  every in-graph k-mer: the answer cannot show the dropped source's own support. A discovery stopped before it read
  an in-graph k-mer (§4.5: the resolve of a prefix without one) has no profile to rank: resolve_discover_contains
  reads it inconclusive (it failed "not discovered in a complete list of 0 labels"; latent, no discovery scenario
  has a budget). Every recorded discovery of v3 and v4 states max_labels and keeps exactly that many when
  truncated: no outcome or detail moves.
- /pattern, low_complexity_note and the invariants read the stop's phase (SPEC-pattern-search §7.8, §8.10, §20.5):
  the note is the engine's, diagnosed after a search that completed or stopped at a threshold
  (pattern_search.cpp decides it on engine.stop()), so only a budget stop of the search in an engine phase
  (discovery, mask_scan, extraction, extension; max_steps, time, max_annotation_work, max_memory; a stop without a
  phase counts as the engine's) leaves it out. A stop of a later phase (label_discovery, placement, output, a
  predicate's selection) leaves the note as the engine decided it: the note follows sdust there (present and
  expected passes, absent and expected fails). There low_complexity_note failed a present note beside {output,
  time}, {output, max_memory}, {label_discovery, max_annotation_work} or {placement, max_memory} and read a missing
  one inconclusive; the invariants, whose vocabulary was max_steps and time, failed it beside {output, time}. A
  supported-path search's selection stop can be its own refusal of a held walk inside the search
  (§20.9) or come after it: a present note passes only where sdust flags the pattern, a missing one is inconclusive.
  The invariants fail a note beside every engine-phase budget stop (before: max_steps and time only; now also the
  supported-path search's {extension, max_annotation_work | max_memory} and the path extension's {extension,
  max_memory}, which the engine counts as budget stops). Measured: no recorded answer of v3 or v4 and no stored
  fixture states the note beside any stop, and no exact pattern sdust flags meets a later-phase stop there (14
  fixture and 11 recorded later-phase stops, none of a flagged pattern); no recorded outcome or detail moves.
- /resolve, resolve_runs on a stopped answer without prefix_if_stopped / unless_stopped: a stop answers the resolve
  of a prefix of resolved_kmers k-mers (SPEC-labeled-traversal-core §4.5), so an expectation about the whole query
  was judged on the prefix and passed ("3 of 3 k-mers" for `full`), against the suite's rule that a stop before the
  property is never a pass. Each argument is now decided only where the prefix decides it: it fails where the
  prefix contradicts it (`full`: a hole in the prefix; `covers`: an uncovered part inside the prefix; `run_end_in`:
  the run from 0 ending inside the prefix outside the range, or reaching the prefix's end past it; `run_start_in`:
  no run to the prefix's end starting in the range when the range lies inside the prefix; `runs`: the prefix's runs
  other than the expected ones cut at the prefix; `kmers` / `max_kmers` / `min_kmers`: more support than allowed,
  or too few k-mers left to reach the count), passes where the prefix proves it (a covered interval or a run end
  inside the prefix, a minimum already reached), and is inconclusive otherwise (a stop before the first k-mer
  shows nothing). 68 `full` expectations of v4 and 31 others are at risk only beyond 19.75 s of work on at most
  2,000 bp: no recorded /resolve answer of v3 or v4 states a stop, so no outcome moves.
- /resolve and /pattern, budget_respected on an allowed 503 `deadline` (`_budget_503`): the premature floor is the
  budget itself on both routes (1 ms slack). It was the budget less /resolve's 250 ms finalisation reserve, so a 503
  up to 250 ms before the budget read inconclusive; the host throws ResolveDeadline only once the elapsed time
  reaches the budget (Deadline::respond_expired, elapsed >= budget; the reserve ends the work earlier, not the
  answer), its clock starting inside the client's interval, so a correct host cannot answer one there. Budget
  1,000 ms: wall 900 (any wall below 999 ms) fails, where it read inconclusive; 999 and 1,000 are inconclusive (the
  1 ms slack); 600 fails as before. Every recorded 503 deadline of v3 and v4 (R-budget, R-twin, R-pos-0027, the
  /pattern ones) came at or after its budget (the smallest /resolve margin: R-twin-0004's 20,776 ms for 20,000, v3
  direct-local): no outcome moves. README lines of lanes S5 and S5b corrected.
- /traverse, traverse_stop with cap_reason: some requested arm of the walked seed must be truncated (or pruned) by
  that cap and every other one either so or `complete`; it required every arm truncated by the cap, a false FAIL
  under direction both, where an arm whose walks all ended for a semantic reason before the cap stays complete
  (SPEC-labeled-traversal-core §6.7), for a per-arm cap (max_output_bp, max_paths, max_live_paths) and a per-seed
  stop (max_steps, time_budget, resource_limit) alike. Both arms complete fails ("no requested arm truncated by
  ..."), an arm truncated by another cap fails as before. Recorded: every traverse_stop answer of v3 and v4
  (T-budget-0001..0004, T-caps-0001) has each requested arm truncated by the cap; no outcome moves.

### Validation (offline)
- The findings' failing inputs with their controls (lane-s-checks/scratch/probes/findings.py, 37 probes): at
  14552938 15 read as they must, here all 37. The tests that state these rules (five new, four extended or
  corrected, among them test_error_answer_skips_other_checks, which asserted the 250 ms floor) fail at 14552938
  (9 of 77 in test_checks) and pass here; the verification's amendments (the truncation's size, a discovery stopped
  before an in-graph k-mer) fail at the lane's first head f07a2fd3 too.
- Every recorded answer of v3 and v4 (5,716 answers, all arms) judged by the checks of 14552938 and these, check by
  check: no outcome moves. Every inexact count of v3, v4 and the stored pattern fixtures judged against a minimum
  one above it and relation exact (stress_count.py): no failure the checks of 14552938 did not give.
- Re-evaluation of the recorded runs (v3: lane S5e's reeval-S5e copy; v4: the complete run copy; `run_suite.sh
  --reevaluate COPY --no-pdf --reeval-mcp-auth anonymous`, against the same copies re-evaluated at 14552938): 0
  outcome changes and 0 check outcome changes of 2,900 (v3) and 2,913 (v4) executions, 0 relation changes (66 per
  run), compare.reevaluated.json identical (summary and every scenario entry), the report identical but its
  timestamps and the checks' provenance; the v4 deploy gate (gate/policy_v4_rt_mcp.json on the re-evaluated set)
  states the same 45 problems. Detail changes only: R-discover-0001's resolve_discover_contains (direct-local and
  direct-remote r0, MCP r0 and r1, both runs) names the tie, and count_relation's inconclusive on M-planted-0017
  (v3 MCP r0) and M-planted-0016 (v4 MCP r0) names the relation expected. The verification's amendments move
  nothing more (re-evaluated again at their head: the same five detail changes per run, the same 45 gate problems).
- 350 tests pass (1 skipped).

## Suite integration S7i: lanes S7r, S7v and S7c merged, and the rules that span them (2026-10-15)
Lanes S7r (runner), S7v (verdict) and S7c (checks) merged onto 81c620ef in that order; the rules below join what
one lane's verification left to another lane's file.

### One reading of an outcome outside the vocabulary
- Two mechanisms met at the merge: merge_run (S7r) read an outcome outside the vocabulary as `error` before the
  comparison ran, so the comparator's `other` bucket, its by-outcome MCP disagreement and the report's "outside the
  vocabulary" sentence (S7v) never fired on a merged set; the run read FAIL as N errored executions, the recorded
  outcome visible only in a merge warning.
- One test, `suite_common.outcome_known`, and one place that rewrites the line, merge_run (`error`, the recorded
  outcome kept under `outcome_runner`). `suite_common.recorded_outcome` reads it back: for a line whose outcome is
  error and whose outcome_runner is present and outside the vocabulary, that recorded outcome, else the line's
  outcome. compare's totals, per-area counts, arm views, the direct pair's and the MCP arm's by-outcome
  disagreement and `outcomes_outside_vocabulary` (with `read_as: error` where the merge rewrote it) read
  recorded_outcome, so a merged and an unmerged set read alike; such a line counts under `other` and not under
  error as well (the counts partition the executions). The report's Failures section lists it as "outcome '...'
  outside the vocabulary (read as error)"; the gate's outcome problem names it ("error (outcome '...' outside the
  vocabulary)"). The paired loop counts it as an error (`unknown` from mcp_step_result, lane S7r). FAIL either way.
- An outcome stored as a JSON list or object (no runner writes one) made compare raise TypeError (an unhashable
  value in aggregate_outcome's set). `suite_common.outcome_value` reads it as its canonical JSON text, outside the
  vocabulary: recorded_outcome returns it so, and aggregate_outcome, the report's untested groups and the gate's
  outcome problem read outcomes through it; such a line is error / FAIL like any other outcome outside the
  vocabulary, merged or not (VERIFY-SUITE-INT follow-up 2).

### The re-evaluation's earlier repeats (P3-run-measure-report-01 in reevaluate.py)
- run_direct.py (lane S7r) and async's runner pass in `previous` only the earlier repeats the host answered (an HTTP
  status) and the runner did not end in error (`suite_common.feeds_previous`): such a repeat's digest is its
  stand-in's, and the next repeat's `deterministic` would fail although the host answered nothing. reevaluate.py
  still passed every earlier repeat with an answer, so a run recorded with an errored repeat and then re-evaluated
  got back the false determinism failure the recording avoided. It applies feeds_previous to the recorded line, as
  the runners do (test_reevaluate's test_unanswered_repeat_feeds_no_previous).

### The comparison judges a recorded run with the run's snapshot (P6-end-to-end-03, the comparison stage)
- Lane S7r pinned the arms to the run's snapshot (RUN_DIR/scripts), but compare.py loaded checks.py next to itself
  for checks.evaluate_relation and the /traverse content rule and read this directory's service_by_design.json; the
  gate's and the report's property reading (`is_property`) did the same. A checks change landing in the tree during
  a run reached the recorded verdict silently: the arms all `pass` with the snapshot, checks_identity clean, the
  verdict PASS -> FAIL ("2 failed relation(s)", VERIFY-S-RUNNER §1).
- `compare.judging_checks(run_dir, reevaluated)` names the checks.py that judges a set beyond its arms: a recorded
  set's snapshot when RUN_DIR/scripts holds checks.py (source `snapshot`; `snapshot_changed` lists its files whose
  sha256 is no longer SNAPSHOT.json's), a re-evaluated set's re-evaluation checks (`reevaluation`), else this
  directory's (`live`: a run older than the snapshot, as v3 and v4). compare's relations and content rule, gate.py's
  and report.py's property reading load it. `compare.table_for` reads, when no --by-design is given, a recorded
  set's snapshot table (`snapshot_table`: the one its MCP arm read beside {checks_dir}); a re-evaluated set keeps the
  table its re-evaluation replayed. run_suite.sh passes `--by-design "$SNAP/service_by_design.json"` to compare
  itself unless --by-design names another table.
- compare.json `judged_with` {source, dir, checks_sha256, snapshot_changed, by_design: snapshot | live | other |
  null}; a path inside the run directory (the snapshot's `dir`, by_design_table's `path`) is stored relative to it
  (`compare.run_relative`; report.report_cases reads it against RUN_DIR), so a report regenerated after the run
  directory moved reads the same table and states the same MCP limits (VERIFY-SUITE-INT follow-up 3; with the
  absolute path the moved run lost its table and read INCONCLUSIVE). The report's Environment table states it ("comparison judged with": the run's snapshot, the
  re-evaluation's checks, or this tree's checks.py with the remark that a change of the tree since the run reaches
  the comparison) and the table's origin beside its sha256. merge_run's stricter outcomes already used the snapshot.
- The snapshot's checks.py imports the snapshot's suite_common.py in the verdict pipeline too
  (`suite_common.load_checks_bound`, reevaluate.load_checks' technique: the caller's suite_common module set aside
  and the checks directory first on sys.path while checks.py loads, both restored after), in
  compare.load_checks_module (the content rule, the relations, the gate's and the report's property reading),
  compare.load_relation_evaluator and merge_run.load_scenario_outcome: a change of a suite_common helper the checks
  call (answer_digest, normalise_direct) landing in the tree after the run cannot move its verdict (VERIFY-SUITE-INT
  follow-up 1). The pipeline itself (merge_run, compare, gate, report) keeps this directory's suite_common, as it
  keeps its own code: compare.json judged_with records `suite_common_sha256` (the checks'), `pipeline_suite_common_
  sha256` and `live_suite_common_differs`; when they differ the report's "comparison judged with" row adds a note
  and the gate a note, never a problem (re-judging an old run with a newer pipeline is intended).
- tests/test_run_snapshot.py: ComparisonSnapshotTest (judging_checks / table_for; the gate's and the report's
  property reading and the report's statement) and RunSuiteSnapshotTest.test_a_checks_change_during_the_run_does_
  not_reach_the_comparison (a copy of the suite, one same_answer relation, an evaluate_relation failing every
  relation appended to the copy's checks.py after the first request: REL-1 passes on both direct arms, PASS).

### The report reads the run side's records (mixed releases, failed steps, lost MCP executions, the sequential loop)
- Lane S7r records the host's release at each session's start and end (merge_run.release_identity: run.json
  `release_identity` {releases, index_fps, changed}, `mixed_releases`) and the loop's `failed_steps` and
  `mcp_missing`; the report read none of them: a run over two deployments read PASS without --gate, and a sequential
  loop record (no per-scenario counts) printed "loop done, - of - scenario(s)".
- A run with mixed releases is never PASS: after FAIL, FINDINGS and the INCONCLUSIVE of a run where nothing passed
  or a part tested nothing, its verdict is "INCONCLUSIVE (the run spans more than one deployment of the host: <each
  change, naming the arm, session and time that saw the new release or index_fp> (releases seen: ...); no failure,
  no disagreement)"; with another verdict a line below it states the same. The gate already failed it (identity).
- Below the verdict, the loop's failed steps (ARM:SCENARIO:EXIT) and the MCP executions the runner did not record
  (SCENARIO:rREPEAT) are listed whatever the loop's state (a capped loop keeps its failed steps). The Environment
  table's mode line states the last invocation's loop as recorded (`report.loop_text`): state (reason), the paired
  loop's completed of planned scenarios, the direct requests, the failed steps and the MCP executions not recorded; a
  sequential loop prints its own counts and no "- of -".

### A resume says which arm the run's first flags skip again (VERIFY-S-RUNNER follow-up, minor)
- A run's flags are its first invocation's (lane S7r): a run started with --skip-local and resumed once with
  --local holds direct-local executions in its plan; if that resume is cut, the next plain resume skips direct-local
  again and the run stays INCOMPLETE, while the log listed only "skip-local" among the run's own options.
- A resume that restores --skip-local / --skip-remote from the first invocation now logs it on a line of its own:
  "the run's first invocation gave --skip-remote, so this resume skips direct-remote too", with, when the run's plan
  (gate.run_plan) holds executions of that arm the arm's records lack, "the run's plan holds N direct-remote
  execution(s), M of them not recorded, which stay so (the run reads INCOMPLETE) until a resume with --remote", else
  "(--remote runs it)". The options line no longer prints a stray space before the comma of a flag without a value.

### What the gate reads from the MCP runner, and request_evidence (agreed with async)
- `request_evidence` (an MCP line's; TABLE-DIFF §2 as amended): {case: {reference path: live | table | unknown}} for
  every case whose request was EVALUATED and names at least one {"host"|"service": ...} reference: a named candidate
  whose service_answer matched, an answer-triggered candidate whose MCP side matched. Every reference that request
  names is resolved and recorded (no short-circuit), the path without its namespace. A case whose request names no
  reference records nothing, and a line with nothing recorded leaves the field out (reevaluate.py writes it only
  when not empty). `request_evidence_mismatches` [{namespace, ref, live, table, case}] beside it on a mismatch.
- The gate's inputs from the MCP runner (async's suite_mcp_runner.py), in arms/mcp/run.json at its top (the latest
  invocation) and per invocation (`invocations[]`):
  - `service_capabilities` {document (GET /traverse/capabilities), host_document (the host's /capabilities as the
    service reads it: pattern, resolve.time_budget, release, index_fp, index_meta_fp where stated), sha256 (the
    agreed suite_common.service_document_sha256, lane S7v), read_at, base, database, mismatches, error}:
    configuration_problems (a
    successful read per invocation, one deployment over all reads by the sha256, the live service_settings and
    advanced_path, the runner's own mismatches) and compare.recorded_documents (an MCP-only run's host references);
  - `checks_dir`, `checks_sha256`, `suite_common_sha256` per invocation (the runner records them since its
    1dccc344, taking the reference code from the suite_common beside --checks-dir, the run's snapshot): merge_run's
    checks_identity, an error for --strict-identity and the gate where the run has a snapshot;
  - `mcp_auth` per invocation (rule 2 of the plausibility rules) and the lines' `request_evidence` above.

### The v4 gate policy pins the table of lane S7v
- gate/policy_v4_rt_mcp.json is regenerated (`gate_policy.py manifest/suite_v4.json --arms mcp,direct-local
  --declared gate/declared_v4_rt.json --name v4-rt-mcp --out gate/policy_v4_rt_mcp.json`) for the table's sha256
  d1e040aac81007713e2250a3dd765244e44026c8bcb03834fd484f143ccebc61 (lane S7v's wall_rule entry; was ab4100c9...).
  Only by_design_table.sha256 and generated_utc differ: selection, exclusions, required executions, properties,
  relations, relations_excluded, counts, declared and manifest are equal (222 selected, 35 excluded, required
  direct-local 220 scenarios / 253 executions and mcp 189 / 384, relations 2 required and 24 excluded).
  test_gate's test_shipped_v4_policy_is_regenerated_exactly, failing since the merge of lane S7v, passes.

### Validation of the integration (offline)
- 446 tests pass (1 skipped; 345 at 81c620ef = 14552938's suite; 431 after the three merges, with two failures:
  the vocabulary test of the two mechanisms and the stale policy pin; 441 before the verification's three
  follow-ups). Each integration rule has a test failing before it.
- Copies of v3 (lane-s5e/reeval-S5e) and v4 (20261014T061053Z_refseq97-6f945db4) re-evaluated at 81c620ef and at
  the integrated head (`run_suite.sh --reevaluate COPY --no-pdf --reeval-mcp-auth anonymous`): 0 changes of outcome,
  outcome_runner, check outcome, arm outcome or relation over 2,900 / 2,913 executions; the class moves are lane
  S7v's (v3 M-planted-0017, v4 N-usecase-0008 and R-pos-0027 MCP not_comparable -> agree), the detail changes lanes
  S7v's (fallback items, execution_differs, host-cache labels) and S7c's (R-discover-0001, M-planted-0016 / -0017);
  v3 FAIL and v4 FINDINGS as before. Against lane S7v's own re-evaluation at d8cf15ff only lane S7c's five detail
  changes per run differ: the integration rules move nothing recorded. The gate on the v4 copy (the regenerated
  policy): 45 problems, the same set as at 81c620ef with the policy of its table. The verification's three
  follow-ups (the snapshot's suite_common bound, list / dict outcomes, run-relative paths) move nothing on v3 / v4
  either (re-evaluated again: identical to the integration's head, the gate the same 45 problems).

## Lane S8c: checks that GPT's second review found too lenient or too strict (2026-10-15)

The checks.py findings of GPT's second review of the suite at c11db9dd (B1, B2, B3, then B5, B6) as their Opus
verification corrected them, and its side note on an MCP walk without its graphlet. The gate is strict: each stricter
rule is measured on every recorded answer of v3 and v4 (all arms), each looser one keeps the failing neighbour it
exists for. The runner vendors checks.py: it takes these rules by re-vendoring (none needs runner logic).

- /traverse, a clamp stands only for what the engine clamps (B1). The engine states a `strategy.clamped` entry for
  four fields only (traverse.cpp, the clamps of a request): bounds.time_budget_ms above the server's max_time_ms (or
  a budget <= 0 raised to it where a seed derives its labels: that budget also bounds the derivation),
  labels.max_seed_labels above its cap, bounds.max_memory_mb and bounds.max_work_units above the server's maxima (or
  omitted and set to them, `requested` "unlimited"). A structural knob out of range is refused, never clamped. The
  checks held a walk to any clamp's effective value and exempted a clamped knob from the MCP preset comparison, so a
  zero-extension branch answer with a fabricated `branching.max_label_branches` clamp to 0 (with or without the echo
  following it, with or without `requested`) passed every check on both arms. A clamp now stands for its knob only
  where it is one the engine can state (_clamped_knob), and the /traverse invariants fail every other entry
  (_clamp_problems): a field outside the four, no `requested`, an effective value that is no positive cap, "unlimited"
  on a field that is not a budget, an entry that is no lowering (beyond the two engine raises above), a field stated
  twice, a `strategy.clamped` that is not a list. Where it is known the entry also agrees with the request (on a direct
  arm, whose request is what the runner sent: the declared value, "unlimited" for an omitted budget; the MCP arm's
  request is the service's, so its `requested` is not compared), with the cap (the scenario's own `clamped`
  expectation of the field, from the build's capabilities snapshot; the checker sees no capabilities document) and
  with the strategy echo, which states the effective value the walk ran under. The legitimate clamps keep passing:
  T-edge-0012's time budget (70,000 -> 60,000) on both direct arms, an omitted work budget set to the maximum, the
  derivation's raise of a zero budget.
- /traverse, a structural end needs its leaf's events (B2). A path that counts a structural end (branch, a quorum or
  split-limit end, a re-entry) is judged by its leaf's `events` (the label_end event of the end, the blocked event of a
  re-entry). The rule that a stated end has its label_end event read the list only where it was one, so an answer
  whose leaf had no `events` (absent, null, an object) skipped it and the branch rules too: the zero-extension branch
  answer with `events` deleted passed every check (a string crashed walk_in_record). The engine writes the list on
  every segment of a tree or full detail (traverse.cpp arm_to_json; segment events are never truncated), so its
  absence is an output or normaliser defect: _structural_end_problems raises Missing for such a leaf, which fails on a
  direct arm and is not exposed on the MCP arm (evaluate's rule), never a pass. walk_in_record treats such a leaf like
  one whose end cannot have happened (it shows none of its ends) and states the missing field only where nothing
  failed; the invariants likewise raise it only when they found no other problem, so a failure is never hidden behind
  it. Events are read through _seg_events (a list's objects; nothing from another type): no crash. The graphlet
  parser gives every segment an `events` list (an MGT states every event, SPEC §7.5; it listed them only on segments
  with an E or R record), so the MCP arm's parsed walks cannot lack it. An end that is no structural one needs no
  events. The unit fixture `walk(..., reason="branch")` stated no events on its leaf and passed through the hole; it
  states the label_end event as the engine writes it.
- /traverse, a walk without its tree is not shown (the verification's side note). An MCP answer whose result lost
  its graphlet (or an arm in summary detail) has arms without segments and paths: walk_in_record failed it "no path
  (status complete)", a wrong diagnosis (it blocks the gate either way). Such an arm is now a missing field:
  not exposed on the MCP arm, a fail on a direct one (a summary detail where the scenario asks for a walk), stated
  only where nothing failed; expect_branching does not call a split missing from an arm it cannot see. Measured on
  v4: of the 123 recorded MCP answers with a graphlet, with the graphlets removed, walk_in_record read fail on 100
  that passed; now no check turns to fail.
- /resolve, a stop's time is judged on a clock that holds the deadline's (B3). The engine starts a budgeted
  resolve's deadline once the body is parsed and its timing.elapsed_ms clock later, after the request's fields and
  the LabelOracle (traverse.cpp process_resolve_request), so elapsed_ms is no lower bound of the deadline's time: a
  valid stop at the work deadline (budget less the reserve) states less than that whenever the gap is long (a stall:
  microseconds of work, so only an OS pause makes it long; the verification found no recorded resolve stop). Three
  sites judged a stop premature by elapsed_ms (the invariants, resolve_runs on any stop, and deterministic's "time
  limitation that cannot have happened"). They now judge it only against the client's wall, whose interval holds the
  deadline's: a direct runner's send to last byte, the MCP arm's service-measured backend exchange
  (mcp_phases.backend_ms, the service's monotonic interval around the backend call); a wall below the budget less the
  reserve less 1 ms fails, no wall leaves the time unjudged (the prefix arithmetic still holds). evaluate passes its
  `client` to these three checks (_CLIENT_CHECKS); its signature is unchanged, and both runners and reevaluate.py
  already supply it. /traverse and /pattern read their stops and their elapsed times on one clock and keep their
  rules; the upper bound of budget_respected on a /resolve 200 stays on elapsed_ms (lenient by the gap, never a false
  fail). A same-origin engine field (the deadline's elapsed time) would let the rule judge without the transport in
  the wall; that is an engine change, proposed by the verification, not made here.
- /pattern, a cut after an earlier stop (B5). `stop` names the first stop in time (SPEC-pattern-search §7.6: every
  stop is assigned only while `stop` is null), the cut the first that applies in §19.8's order (the engine's stop,
  the pass's, the raw release's max_predicate_contexts, the list's max_contexts; the output's cut replaces them:
  pattern.cpp put_selection). The checks required cut max_predicate_work to equal the stop, so a partial predicate
  search whose raw discovery stopped at max_predicate_contexts (stop_at_threshold) and whose pass then spent its work
  failed, and four more legitimate compositions with it (a later pattern finding the request's selection work spent,
  sticky, §19.8, in partial and in all_or_count; a list's memory cut after an earlier extension or discovery stop:
  put_plain_cut, apply_labels). The evidence is the selection work the answer shows: work.predicate_units, per
  pattern, summed over the patterns so far (one budget per request, §19.9: the pass stops once the units reach
  max_predicate_work). cut max_predicate_work now needs a predicate and a stop; beside another stop it fails where the
  work is shown below the budget, or where this pattern's own pass ran after an engine stop of an engine phase
  (max_steps, time in discovery, mask_scan, extraction: the engine's cut comes first). A max_predicate_work stop with
  the work shown below the budget fails too (it passed). withheld predicate_budget is also the sticky budget's. A
  memory cut stands beside any stated stop: it replaces every other cut and states its own stop only when none is
  stated; it stays lenient there, since the cut states no phase of its own (no cut.phase), so a memory cut that never
  happened beside an earlier stop is not told apart (the verification proposes cut.phase as an additive engine field).
  Pending the engine owner's decision: a selection budget an earlier pattern spent beside this pattern's own engine
  stop (max_steps, time). §19.8's order puts the engine's cut first (cut max_steps / time); its sticky sentence and the
  route (a fresh extraction without the engine's cut in the pass's `before` branch) give cut max_predicate_work in
  partial and withheld predicate_budget in all_or_count. Both are accepted until §19.8 or the route states one.
  Measured: the 3,751 answered pattern entries of the 175 stored fixtures and v3 / v4 (all arms) give the same
  problems as before (none); the verification's legitimate compositions pass, its impossible neighbours fail.
- /pattern, a stop's phase states its reason and is one the request runs (B6). The checks read a stop's reason
  against the request but never its phase, so the recorded P-modes-0005 answer (a label-free threshold stop in
  discovery) passed its whole frozen scenario with the phase moved to placement, label_discovery, output, selection,
  extension, extraction or mask_scan, and a labels-all annotation-work stop passed in discovery; the phase also gates
  other rules (_threshold_value_problem judges only discovery), so a wrong phase silently turned them off. Each stop
  is now checked last in _pattern_reason_problems against SPEC-pattern-search §8.7's table (the engine's StopPhase,
  every set_stop of the route): discovery max_steps, time, max_contexts, max_anchors, max_predicate_contexts;
  mask_scan max_steps, time; extraction time; extension max_steps, time, max_paths, max_memory, max_annotation_work,
  max_predicate_contexts; label_discovery and placement max_annotation_work, max_memory, time; output max_memory,
  time; selection max_predicate_work, max_memory, time, max_contexts, max_annotation_work. And against the request:
  extraction needs a release (a retrieval mode or a predicate), extension a path search of a pattern longer than k,
  label_discovery and placement a retrieval mode reading labels (all, predicate_only), output a retrieval mode that
  builds more than the engine's release (labels, a predicate, a path search) or a motif's evaluation, selection a
  predicate; max_contexts in discovery not beside a predicate's pass (put_selection renames it
  max_predicate_contexts), in selection only on L <= k; max_predicate_contexts in discovery only on L <= k, in
  extension only for supported paths; max_annotation_work in extension and selection only for supported paths. A stop
  without a phase stays the schema's. Measured: the 26 impossible phase edits of the verification's cases fail (23
  passed), the legitimate neighbours (placement/max_annotation_work, mask_scan/max_steps, placement/max_memory,
  extraction/time) and every stored and recorded stop pass. The unit test of the low-complexity note beside a later
  phase stated a selection stop without a predicate and an extension stop without a path search; it states them on
  requests that run those phases.
- /traverse, an annotate walk's re-entry end needs its blocked event (the follow-up's R1, B2's sibling). In annotate
  mode the walker ends a path at a re-entry by its path_reason alone (edge_reuse, edge_reuse_rc, rejoined_seed): no
  lineage, so end_reasons is empty (walker.cpp, finish_path with the block's reason; the blocked event is pushed at
  that step). The B2 rule keyed on end_reasons and never saw such an end, while _leaf_ends counts the path_reason
  among the ends that explain a short walk: the recorded T-labels-0001 answer, its right arm cut at 50 bp with
  path_reason edge_reuse, passed on both arms with its leaf's events absent or empty. A structural path_reason now
  needs the leaf's events list (Missing without it, as a counted end), and a re-entry path_reason its blocked event at
  the path's end (in any mode: "blocked continuations are always reported"). The cut with events [] fails on both
  arms, the cut without the list fails on a direct arm (the MCP arm: see the events-list rule below), the cut that
  states its blocked event and the recorded answer pass. A re-entry that states its blocked event is still accepted
  where the truth has no exact minimum (in constrain mode too), the verification's remaining design limit.
- /resolve, a clamp is the route's (the follow-up of B1 for /resolve). The route clamps one field: a time budget above
  the server's resolve.time_budget.max_time_ms is lowered to it and stated in limits.clamped (`requested` the value
  as given, `effective` the cap), and limits.time_budget_ms states the budget the work ran under (traverse.cpp
  process_resolve_request; graph/traversal/resolve.cpp clamps nothing; a budget not above the finalisation reserve is
  refused). The resolve invariants now fail what the route cannot state, as the /traverse rule does: another field,
  no `requested`, an effective value that is no positive cap, an entry that is no lowering, a field twice, a
  `limits.clamped` that is not a list; where known, a `requested` other than the direct arm's request (the MCP arm's
  request is the service's), an effective value other than the cap (the scenario's `clamped` expectation, or
  budget_respected's cap_ms) or than limits.time_budget_ms. And, on a direct arm, a limits.time_budget_ms below the
  request's budget without that clamp: the stop rule reads that budget, so a budget lowered unstated would excuse a
  premature stop (a stop within 1,100 ms of a 3,000 ms request passed beside a stated budget of 1,000). The recorded
  clamp (R-edge-0021, 70,000 -> 60,000, both direct arms) and every recorded budget pass.
- /traverse, every segment states its events list (the stricter rule the B2 item left optional). The engine writes
  `events`, a list of event objects, on every segment of a tree or full detail (traverse.cpp arm_to_json; never
  truncated), and graphlet_arms gives every parsed segment its list. The /traverse invariants now fail a segment
  without the list (absent, null, another type, or entries that are not objects) on every arm, whatever ends its paths
  claim. On a direct arm this only adds a clear diagnosis beside the structural-end rules; on the MCP arm it closes
  the verification's theoretical gap: a JSON segment there without its events, on a path walk_in_record does not
  read, gave not_exposed invariants beside a passing walk, which the gate tolerates for a generic check. Measured:
  every recorded JSON segment of v3 / v4 (3,912) and of the 41 stored full documents (1,428) states its list; no
  recorded or stored verdict moves. The unit fixture walk_answer stated no list and states an empty one.
- /traverse, a budget echoed otherwise than requested needs its clamp (the /resolve rule's counterpart). The engine
  echoes the request's time budget, seed-label limit, memory and work budgets as sent unless it clamped them (and then
  states the clamp), and the walk's checks read the echo: a time stop 2 s into a 15 s request passed beside an echo of
  a 2 s budget. On a direct arm, a clampable knob the request sent whose strategy echo differs from it, lower or
  higher, without a strategy.clamped entry of that field now fails the invariants (_clamp_problems; an invalid entry
  of the field is reported by its own rule). The engine's clamps (with their effective values echoed) and its two
  raises (a zero budget raised with its clamp, an omitted budget set to the maximum: not sent, so not compared) pass.
  The MCP arm's request is the service's (its presets and its own bounds): not compared. Measured: 0 changes on v3 /
  v4, the v4 gate the same 45 problems, the 82 stored /traverse documents unchanged.

### Validation (offline)
- 455 tests pass (1 skipped, the record-cache builder test; 446 at c11db9dd, 451 before the follow-up). Each item's
  new or revised test fails on c11db9dd's checks.py: the structural-clamp rejection, the missing-events and
  missing-tree tests, the resolve stop's pre-timer and premature-wall cases, the compound cut and the phase tests (the
  low-complexity-note test is a correction and passes on both). The follow-up's tests fail on the commit before
  theirs: the annotate re-entry, the /resolve clamp, the events-list rule (the B2 and R1 tests' MCP invariants), the
  /traverse budget echoed otherwise (fails at aa5e5f48);
  the test of the discovery max_contexts stop beside a predicate's pass covers a rule already there and fails on
  c11db9dd and on a mutant without that rule, which no other test caught.
- Copies of v3 (lane-s5e/reeval-S5e) and v4 (doc/v4-c11db9dd) re-evaluated after each item and at the head against
  c11db9dd (`run_suite.sh --reevaluate COPY --no-pdf --reeval-mcp-auth anonymous`, v4 with `--gate
  gate/policy_v4_rt_mcp.json`): 0 changes of outcome, outcome_runner, check outcome or check detail over 2,900 /
  2,913 executions (every arm and the merged set); compare.reevaluated.json identical but for the table's path and
  the re-evaluation mark; v3 FAIL and v4 FINDINGS as before; the v4 gate 45 problems (configuration 13, outcome 31,
  relation 1), the same set.
- The 3,751 answered /pattern entries of the 175 stored fixtures and v3 / v4 (all arms): _pattern_reason_problems
  unchanged (none on either side); the 175 fixtures' invariants and counts_consistent on both arms: unchanged. The
  verification's 50 B5 / B6 cases (legitimate compositions, single-stop controls, impossible neighbours, the full
  P-modes-0005 scenario with each phase edit): all as intended. Of v4's 123 MCP traverse answers with a graphlet,
  stripped of it, no check turns to fail (walk_in_record failed 100).
- The follow-up (R1, the /resolve clamp, the events-list rule), re-evaluated after each item and against c11db9dd: 0
  changes on v3 / v4, compare.reevaluated.json identical but for the table's path and the re-evaluation mark, the v4
  gate the same 45 problems. The 82 stored /traverse documents (41 full, 41 graphlet; invariants on both arms) and the
  175 pattern fixtures: unchanged, all /traverse documents pass. The verification's T-labels-0001 probe (the recorded
  answer cut short by path_reason edge_reuse): with the leaf's events absent or empty it fails its whole scenario on
  both arms; with its blocked event, and the recorded answer itself, it passes.

## Lane S8g: the gate's evidence, the deadline window and the gated headline (2026-10-15)
GPT's second review of the suite (~/.cache/metagraph-traverse-reviews/2026-10-15-gpt-suite-2/REVIEW.md) findings B4,
B7 and B8, as their Opus verification (verify/VERIFY-G.md) corrected and widened them, against c11db9dd.

### What a GATE PASS needs as evidence of identity (B4)
A GATE PASS asserts which checks judged the run and which index answered it. It accepted missing evidence for both,
and a recorded hash it never recomputed:
- **Both sha256 on every MCP invocation.** merge_run.checks_identity read an MCP invocation as verified when it
  named a checks directory that is the run's snapshot, or stated only checks_sha256: a runner that loaded checks.py
  from the snapshot and imported another suite_common (the case suite_common_sha256 exists for) passed. With a
  snapshot every invocation now states both checks_sha256 and suite_common_sha256, each the snapshot's (checks_stub.py's
  for an arm whose checks_module names the stub, as for the direct arms); one missing either is `unverified` ("do not
  state both checks_sha256 and suite_common_sha256 (... missing on N; a checks directory does not stand in for
  them)"), an error for --strict-identity and the gate. The checks directory is informational; one naming another
  directory without a checks sha256 still `differs`. async's runner records both since its 1dccc344.
- **A recorded set without a snapshot blocks.** A run before lane S7r (v3, v4) has no RUN_DIR/scripts/SNAPSHOT.json,
  so nothing pins the checks its arms were judged with; checks_identity's `unpinned` names it ("the run has no
  snapshot (scripts/SNAPSHOT.json): the checks that judged its arms are not pinned (re-evaluate it to gate its
  answers)"), an error of the recorded set for --strict-identity and the gate (no merge warning: such a run stays
  readable as before). Its re-evaluated set is judged with the re-evaluation's checks (reevaluated.json) and is not
  blocked on that ground.
- **Every policy arm observed the pinned index.** The gate compared only a top-level `index_fp` of an arm's run.json,
  which async's runner never writes, and noted "not compared" when no arm stated one: the MCP arm's observed index
  (service_capabilities.host_document.index_fp, the service's backends.<db>.identity) was never compared, and a
  conflicting one passed. gate.host_identity_problems now reads merge_run.observed_host_identity per arm (a direct
  arm's run.json and its sessions' `<key>_start` / `<key>_end`; the MCP arm's host_document of run.json's record and
  every invocation's) for each of suite_common.HOST_IDENTITY_KEYS the policy pins (index_fp; release and
  index_meta_fp where a policy pins them): a policy arm that observed none ("arms/mcp states no observed index_fp;
  the policy pins ... (service_capabilities.host_document)") and any arm that observed another value ("arms/X
  observed index_fp A, B; the policy pins P") are identity problems. A direct arm outside the policy need not
  observe it. The tests' fixtures give the MCP arm no top-level index_fp (test_compare_report.write_arm), as the
  runner does.
- **Recomputed document identities.** The gate's one-document test and compare.recorded_documents read a record's
  stated sha256 in place of its identity: a changed document under a stale hash was one document, and the
  advanced_path premise, read from the first record only, never saw the change. Each record's identity is
  recomputed (suite_common.service_document_sha256) and only recomputed values are compared; a recorded sha256 the
  record does not reproduce is a configuration problem ("<where>: the recorded capabilities sha256 X is not its
  record's Y") when the run's snapshot's suite_common is the pipeline's (one definition computed both). Where the
  snapshot's suite_common is another, the recorded sha256 is checked with that file's own
  service_document_sha256, the definition the runner computed it with (gate._snapshot_document_sha: the file still
  SNAPSHOT.json's and loading; the problem then says "(the run's snapshot's suite_common.service_document_sha256)");
  only where that cannot be read (no snapshot suite_common, one changed since the snapshot, one that does not load or
  has no such function) is it a note, compared recomputed with this one. compare.recorded_documents states the
  recorded one as `sha256_recorded`. tests/fake_mcp.py records the real sha256.
- On the v4 copy re-evaluated (the shipped policy): 46 problems, the 45 of c11db9dd plus one identity problem, "mcp:
  arms/mcp states no observed index_fp; the policy pins 907ec986e19ecb04 (service_capabilities.host_document)" (v4's
  MCP runner recorded no service_capabilities). The recorded v4 set gains the unpinned problem as well.

### A shortened deadline window is placed, not failed (B7)
suite_common.stated_bound_problem failed every deadline-phase resolve_timeout whose window was below the read floor
t - 1 s. The service takes the exchange's absolute deadline (now + B, B = exchange_budget_s: t plus the connect,
write and pool caps, each at most t; 55.25 s for t = 30.25 s) before it acquires the host's token (at most A =
service_settings.resolve.acquire_timeout_ms, 2 s) and claims the attempt, and states the window that is left when
the exchange starts: a slow claim (a database pool wait, a stalled worker) legitimately leaves less than the floor,
and that load event read as a service failure. Both outcomes block the strict gate; the diagnosis was wrong.
- **One definition: suite_common.deadline_window_verdict(scenario, stated, evidence, refs)** -> (state, detail),
  which async's runner takes through shared_rule with its own evidence record. In order:
  1. a recorded exchange_budget_s other than B (beyond 1 ms): not_plausible, "the attempt's exchange budget X s is
     not t + the phase caps B s";
  2. bound_s above B: not_plausible, "service deadline window b s above its exchange budget B s";
  3. the statement's `pre_exchange_ms` (from the deadline's start, before the acquire and the claim, to the
     exchange's start) and `acquire_ms` (the acquire's own time), which the service states on every resolve_timeout
     from the async session's next round (each rounded to 0.1 ms by the service, absent when unknown), are
     preferred: acquire_ms above pre_exchange_ms is not_plausible ("the stated acquire a ms exceeds the p ms before
     the exchange it is part of"); bound_s x 1000 + pre_exchange_ms must be B x 1000 within 2 ms
     (DEADLINE_TERMS_TOLERANCE_MS: bound_s is rounded to 1 ms, the origin terms to 0.1 ms; 2 ms covers both), else
     not_plausible ("service deadline window b s ended earlier than its terms allow: exchange budget B s - p ms
     before the exchange = w s", or "... is longer than its terms allow: ..."), at any bound;
  4. at or above the read floor t - 1 s: no objection, as before;
  5. below it the failing attempt's span (from its grant, the attempt's dispatched_at, to its settling; the service
     grants at most max(A, acquire_ms) after the deadline's start and settles after the deadline) places the
     origin: a span below B - max(A, acquire_ms) - 1 s is not_plausible ("the deadline began before the acquire:
     attempt span X ms (attempt) < exchange budget B s - acquire A s - 1 s = Y ms"); otherwise no objection (load);
     no span: no objection when step 3 placed the origin, else undecidable with rule (c)'s class
     UNDECIDABLE_NO_ORIGIN ("no attempt span places the deadline's origin: plausibility undecidable; deadline window
     b s below its floor f s"), which merge_run and compare read as an undecidable inconclusive.
  No objection leaves rule 4's wall rule (the stated bound - 1 s): inconclusive as a host state.
- **The evidence: suite_common.deadline_evidence(line, doc)** reads the runner's record `deadline_evidence`
  {exchange_budget_s, timeout_s, acquire_timeout_s, granted_at, settled_at, span_ms, span_source "attempt" | "job" |
  null} (the raw record's, else the line's), cross-checked with deadline_evidence_from_record (a different span_ms
  or exchange_budget_s is a warning; the recorded one is read); without one it is built from the raw record: the
  last settled attempt of the runner's GET /admin/traversal/jobs/{id} read (reserved.bound_ms / 1000, reserved
  .timeout_s, dispatched_at, settled_at), else phase_sources (the attempts' last span_ms; the job's finished_at -
  started_at, "job", longer with retries, so never a false fail) and the job plan (traverse_results overview.plan
  max_bound_s, timeout_s); acquire_timeout_s is the service's. Never mcp_phases.backend_ms (the exchange alone).
- timeout_statement carries pre_exchange_ms / acquire_ms from the one error object it reads (non-negative numbers,
  absent otherwise); timeout_terms takes them from the record's error when a runner's recorded statement lacks them
  and agrees otherwise. TERMS gains connect_timeout_s, pool_timeout_s and acquire_timeout_ms (the live
  service_settings.resolve, else the table's: 10 s, 5 s, 2000 ms).
- No recorded resolve timeout of v3 or v4 states a phase (all legacy records): 0 recorded verdict changes.
- **What async's runner mirrors in the same round**: stated_bound_problem's deadline branch replaced by
  shared_rule("deadline_window_verdict")(scenario, stated, evidence, refs) (its write branch unchanged; the
  undecidable state answers as rule (c) with shared_rule("UNDECIDABLE")'s UNDECIDABLE_NO_ORIGIN as
  inconclusive_reason); its timeout_statement carrying pre_exchange_ms and acquire_ms; the record
  `deadline_evidence` (DEADLINE_EVIDENCE_KEYS) beside `timeout_terms` on the raw record and the line, and
  `dispatched_at` / `settled_at` kept on phase_sources.attempts rows.

### The gated headline, and UNVERIFIED for what passed nothing (B8)
A run whose MCP selection was all marks (not implemented via MCP, service-side by design) headlined PASS ("no failed
execution ..."), and a gated report showed the ordinary verdict before the gate's: a reader skimming it took the
absence of a failure for verification.
- **A gated report leads with the gate.** build_report renders the deploy gate's section (GATE PASS, GATE FAIL with
  its problems, or GATE verdict STALE) before the verdict whenever the run has a gate verdict, and the ordinary
  verdict then reads "**Exploratory verdict (not the deploy decision, which is the gate's above): ...**".
- **Nothing passed is never PASS.** A run with executions and no pass that no earlier rule caught (INCOMPLETE, FAIL,
  FINDINGS, INCONCLUSIVE with inconclusive executions, untested groups, mixed releases) reads "UNVERIFIED (no
  execution passed: 1 not implemented via MCP, 2 service-side by design; nothing was shown to hold; no failure, no
  disagreement)"; an arm with executions and no pass reads "UNVERIFIED (no execution passed on mcp: ...; nothing
  was shown to hold there; ...)", whatever another arm passed. Only verdicts that read PASS move; a complete
  selective MCP-only run with passes still reads PASS, and marks-only area or route groups keep the earlier rule.
- **The gate's nothing_passed is per policy arm** ("no execution of arm mcp passed in the policy's selection"): a
  two-arm policy whose MCP selection is entirely excluded verified nothing via MCP and no longer passes on
  direct-local alone (the shipped policy, 384 required MCP executions, cannot reach it).
- v3 FAIL and v4 FINDINGS are unchanged (every arm has passes); the v4 copy gated with the shipped policy leads with
  GATE FAIL.

### Property-level marks of what the MCP cannot show (manifest v5; owner decision 2026-10-15)
From v5 the manifest marks what the MCP cannot show by its true class: `mcp.not_implemented` (the whole scenario, as
before), `mcp.not_implemented_checks` [{check, reason}] (a missing MCP feature at property level, e.g. the label
dictionary of T-labels-0003 / -0005 / -0006; the rest of the scenario stays tested via MCP) and
`mcp.by_design_checks` [{check, case}] (a deliberate service behaviour at property level, e.g. T-edge-0008's
answer_fields under traverse.duplicate_seeds). Lane V5 adds the fields; the names are read in one place,
suite_common.MCP_CHECK_MARKS. (Integration S8i stage 2: gate_policy derives no exemption from them, and T-edge-0008 is
modelled through the table, not a mark.)
- **What a mark sets aside** (suite_common.checks_marked_aside, the one rule the outcome through unmarked_checks,
  the comparison, the gate and the report apply): on the MCP arm, a marked check recorded reading not_exposed, the
  reading the mark declares. A marked check that reads anything else (pass, fail, finding, inconclusive) is judged as
  any check: a mark hides no result. A marked check with no result on the execution is a missing result as before
  the marks (the runner drops nothing: it judges on unmarked_checks).
- **Outcome.** reevaluate.py's outcome rule (and an answered case's verdict) and merge_run.stricter_outcomes judge an
  MCP execution on its checks less the set-aside items (read on the whole execution, before an answered case's
  explained items are taken out), so a marked check never makes it inconclusive ("every property not exposed") and
  so never untested in the report.
- **Every property check marked: a scenario-level mark.** When the marks cover every property check of an MCP
  execution (suite_common.every_property_marked: at least one property check, each marked and read not_exposed; a
  check result's own `property`, else its name outside suite_common.GENERIC_CHECKS), nothing of the scenario's claim
  is tested via MCP: that is what mcp.not_implemented / mcp.by_design state for the whole scenario, and the generic
  checks alone (status, invariants) are no pass. Nothing is set aside then (checks_marked_aside is {}); unmarked_checks
  returns the marked items read inconclusive, so either outcome rule gives inconclusive (fail or finding where another
  check says so). The execution's inconclusive_reason is suite_common.SCENARIO_LEVEL_MARK ("every property check is
  marked; this is a scenario-level mark") with suite_common.scenario_level_mark_detail(names) as verdict_detail: in
  the re-evaluation, in merge_run.stricter_outcomes where it makes a recorded pass inconclusive, and in compare's
  reading of a line recorded without a reason. The gate blocks it as any inconclusive (the problem names the reason),
  and a passing MCP line that has it all the same is a property problem. The comparison and the report treat the
  marked items as unmarked not_exposed ones (skipped, partial depth); the report's mark listing names those
  executions "not set aside, every property check is marked; this is a scenario-level mark". A wrong
  all-properties mark is so visible and gate-blocking; lane V5 turns such a mark into a scenario-level one in the
  manifest (no v5 scenario after integration S8i stage 2, which models T-edge-0008 through the table). The rule is
  per execution, as
  checks.scenario_outcome's "no property passed" rule: a later repeat whose determinism check passed has a passing
  property.
- **Gate.** A marked check the MCP read not_exposed is no required-property problem, and with require_agreement
  no property "not compared with the reference" (gate.uncompared_properties leaves it out). One the MCP did not read
  on a passing execution is "the expected check ... has no result" (a property problem), and uncompared.
- **Comparison.** compare.check_agreement compares a set-aside item with nothing (removed from both sides by check
  and occurrence): it is neither counted as skipped (checks_not_exposed) nor unmatched, and the scenario's
  mcp_agreement lists it as `marked_checks` [{repeat, check, mark}] with `checks_marked`, summed in the summary;
  report.agreement_depth neither compares nor skips it (a reference failure of it is still named as hidden).
- **Report.** The Agreement section lists the marks per scenario under "Not implemented yet via MCP (property)" and
  "Service-side by design (property)" (reason or case, the MCP executions it was set aside on, those where every
  property check is marked, any other reading "judged as usual"), the agreement paragraph counts the set-aside
  items, and the verdict section has one line counting the marked checks by class.
- A manifest without the fields (v1-v4) is judged exactly as before: none of these keys or lines appears, and a
  not_exposed check that is not marked keeps its handling (skipped, partial depth, the gate's property problem
  unless the policy lists it as not required).
- **async's runner mirrors** the outcome part (the runner contract): its outcome rule judges an MCP execution on
  shared_rule("unmarked_checks")(scenario, "mcp", checks), the whole execution's checks, before an answered case's
  explained items are taken out, and drops no item; the marks covering every property check need no other code
  there, unmarked_checks returning the marked items read inconclusive. An execution that ends inconclusive by the
  checks (no rule (c), no host-state class) while shared_rule("every_property_marked")(scenario, "mcp", checks)
  names checks records inconclusive_reason shared_rule("SCENARIO_LEVEL_MARK") and verdict_detail
  shared_rule("scenario_level_mark_detail")(those names) (joined with "; " after a wall rule's notes). Until the
  runner mirrors it, a recorded line it made inconclusive by a marked check stays so in the recorded set (merge_run
  only makes outcomes stricter) and is corrected by a re-evaluation, and a recorded pass whose every property check
  is marked is made inconclusive with that reason by merge_run.

### The by-design table: the discovery cap 1000 (async fb66f1dd, 2026-10-15)
The service's discovery cap goes from 100 to 1000 and its default from 20 to 50 with async's bundle
(fix/owner-1015-pattern-discover fb66f1dd), which lands before v5; the host has no discovery cap (default 1000,
minimum 1). service_by_design.json (sha256 39f87b83e0b0004ace6e062dddf1e95eafd6595f6bbdba6387c10719bec6b228, 133
cases; the note is the lane's TABLE-DIFF.md):
- **Renamed cases** (the ids lane V5's manifest names): resolve.discover_limit_above_cap (max_labels / limit below 1:
  both refuse; request [discover equals true, limit lt 1]) and resolve.discover_limit_above_cap_host_answers (above
  the cap: the service refuses, the host answers; request [discover equals true, limit gt {"service":
  "backends.refseq33m-experimental.resolve.discover.max_labels"}], the path the live GET /traverse/capabilities
  states the cap at). Both state the refusal "`limit` must be an integer in 1..1000" (max_limit 1000), the route's
  "discover.max_labels must be an integer in 1..1000, got <value>" and the default 50; `source` lines at fb66f1dd
  (`source_commit`).
- **Aliases** for the manifests v1-v4: resolve.discover_limit_above_100 and
  resolve.discover_limit_above_100_host_answers stay, `alias_of` the renamed case, without its request (v4's
  R-edge-0023 sends 150, which does not exceed 1000) and
  with the text any cap's refusal contains ("`limit` must be an integer in 1.."), so v1-v4 load and their recorded
  runs read as before. On a deployment with the cap 1000 the service answers R-edge-0023's 150 and its pair reads
  unverified until v5 sends a limit above the cap.
- **The block** states backends.refseq33m-experimental.resolve.discover {max_labels 1000, default_max_labels 50} (a
  reference must be stated by the table's own block; the live document wins, a deployment before fb66f1dd reading
  as a configuration mismatch 100 live / 1000 in the table).
- gate/policy_v4_rt_mcp.json is regenerated for the table's sha256 (only by_design_table and generated_utc differ).
  The builder manifest_lib.py still names the old ids (lane V5's).

### Validation of lane S8g (offline)
- 464 tests (1 skipped; 446 at c11db9dd, 3 skipped there in an export without the C++ source). Every rule change has
  a test failing at c11db9dd: test_gate IdentityEvidenceTest, GatedHeadlineTest and PropertyMarksTest,
  test_run_snapshot's per-invocation sha256 test, test_reevaluate test_resolve_timeout_deadline_window,
  test_compare_report test_discovery_cap_cases; the tests that encoded the old rules are updated (a checks
  directory as the MCP arm's evidence, a recorded sha256 winning, the 15 s deadline window without a span failing,
  the table's 131 cases). Two orchestrator tests that count run_direct.py processes on the machine failed once while
  another lane's tests ran beside them and pass alone.
- Copies of v3 (lane-s5e/reeval-S5e) and v4 (doc/v4-c11db9dd) re-evaluated at c11db9dd and at this lane's head
  (`run_suite.sh --reevaluate COPY --no-pdf --reeval-mcp-auth anonymous`, v4 with `--gate gate/policy_v4_rt_mcp.json`,
  the regenerated policy at the head): 0 changes of outcome, outcome_runner, inconclusive_reason, case, explained
  items, verdict_detail or check outcome over 2,900 / 2,913 executions, 0 changes of a direct or MCP class, a
  by-design verification or a relation; v3 FAIL and v4 FINDINGS as before (v4's report leads with its gate now,
  the FINDINGS line labelled exploratory). The gate on the v4 copy: 46 problems, c11db9dd's 45 and "arms/mcp
  states no observed index_fp; the policy pins 907ec986e19ecb04 (service_capabilities.host_document)". The recorded
  v4 set (not re-evaluated) is blocked as unpinned as well.

## Suite integration S8i: lanes S8c, S8g and the bench merged (2026-10-15)
Lanes S8c (checks.py), S8g (the gate's evidence, the deadline window, the gated headline, the property-level marks,
the discovery-cap table) and the bench panel (bench/, its tests, the orchestrator tests' runner isolation)
cherry-picked onto c11db9dd in that order, 43 commits with their messages and authorship. They share only this README
(S8c and S8g each appended a section after S7i's: both kept, S8c first). Two items join what the S8g verification
left (VERIFY-G round H):

### A scenario-level mark needs every property check the scenario expects marked
every_property_marked read only the property items an MCP line records: a runner that dropped an unmarked property
check left the marked not_exposed one alone, and the execution was named "every property check is marked; this is a
scenario-level mark" (merge_run made the recorded pass inconclusive with that reason) instead of the missing result it
is. It now also needs every property check the scenario's `expect` names (suite_common.expected_property_checks: the
names outside GENERIC_CHECKS) to be marked, determinism aside (suite_common.LATER_REPEATS_ONLY, the gate's list: a
first repeat has no result of it, so the per-execution rule above stands). A dropped unmarked property is so no
scenario-level mark: the marked reading stays set aside, the merge keeps the recorded pass, and the gate names "the
expected check ... has no result on a passing execution" (a property problem). With every expected property marked
and one of them dropped, the scenario-level mark stands and the gate names the dropped one's missing result beside
it. checks_marked_aside and unmarked_checks pass their scenario through, compare.check_agreement takes the
scenario (`scenario`), and the report's agreement depth gives it its own; with marks and no scenario only the
recorded checks are read.

### A snapshot file that is not SNAPSHOT.json's is never executed
merge_run's stricter outcomes, the comparison (its relations and /traverse content rule) and the gate's and the
report's property reading loaded a run snapshot's checks.py, bound to its suite_common.py, even after either changed
against SNAPSHOT.json: the verdict blocked on "the run's snapshot changed after it was taken", but the changed code
ran. suite_common.load_checks_bound now refuses, before executing anything, a directory holding SNAPSHOT.json whose
checks.py or suite_common.py (EXECUTED_FILES) has another sha256 than the recorded one, or whose SNAPSHOT.json states
no files (snapshot_refusal; SnapshotChanged). The merge then keeps the runners' outcomes and warns ("checks.
scenario_outcome not applied, the runners' outcomes stand: ..."); the comparison uses the built-in relations and no
content rule, and warns ("the comparison does not execute the checks that judged the set: ..."); the gate's property
reading falls back to GENERIC_CHECKS and the gate records an identity problem ("the property reading refuses the
checks that judged the set (snapshot): ..."), beside checks_identity's "snapshot changed" for a recorded set;
reevaluate.py stops with the reason. A directory without SNAPSHOT.json (this tree, a re-evaluation's checks) loads as
before; gate._snapshot_document_sha checks the snapshot's suite_common.py against SNAPSHOT.json itself.

### For async: what the runner mirrors after this integration
Re-vendor checks.py (S8c: no runner logic), suite_common.py and service_by_design.json (sha256 39f87b83...; after
stage 2 d4b8636c..., below), then:
- **The deadline window** (S8g): the deadline branch of stated_bound_problem is
  shared_rule("deadline_window_verdict")(scenario, stated, evidence, refs); its undecidable state is rule (c) with
  UNDECIDABLE_NO_ORIGIN (in shared_rule("UNDECIDABLE")); timeout_statement carries pre_exchange_ms / acquire_ms
  (DEADLINE_ORIGIN_KEYS) from the one error object; deadline_evidence (DEADLINE_EVIDENCE_KEYS) is recorded on the raw
  record and the line; the phase_sources.attempts rows keep dispatched_at / settled_at.
- **The marks** (S8g and its round H, and the expected properties above): the outcome rule judges an MCP execution
  on shared_rule("unmarked_checks")(scenario, "mcp", checks), the manifest's scenario (never None: its `expect` is
  read) and the whole execution's checks before an answered case's explained items are taken out, dropping no item.
  An execution inconclusive by the checks (no rule (c), no host-state class) for which
  shared_rule("every_property_marked")(scenario, "mcp", checks) names checks records inconclusive_reason
  shared_rule("SCENARIO_LEVEL_MARK") and verdict_detail shared_rule("scenario_level_mark_detail")(names), after a
  wall rule's notes, "; "-joined. Names: GENERIC_CHECKS, LATER_REPEATS_ONLY, SCENARIO_LEVEL_MARK,
  expected_property_checks, every_property_marked, scenario_level_mark_detail, check_keys, checks_marked_aside,
  marked_aside, unmarked_checks.
- **The loader** (as corrected by the verification, VERIFY-G "Integration stage 1" §5): before the snapshot's
  suite_common.py is first executed (use_shared_common, load_checks, load_common), the runner calls its vendored,
  trusted copy, pinned_common().snapshot_refusal({checks_dir}); a non-empty refusal stops the invocation with the
  message (exit 5) and sends nothing. shared_rule("load_checks_bound") alone is the snapshot's own copy, which a
  changed file need not contain (SnapshotChanged, EXECUTED_FILES, snapshot_refusal).
- Nothing for S8g's identity evidence and headline (the runner records both sha256, host_document.index_fp and the
  record's sha256 already) or for the bench (it drives both runners unchanged: run_direct.py --only, the MCP
  runner's --ids with --repeat 0).

### Validation of the integration (offline)
- 527 tests pass (1 skipped, the record-cache builder test): 446 at c11db9dd, + 9 of S8c, + 20 of S8g, + 50 of
  the bench (its 49 tests and the orchestrator tests' isolation test) = 525 after the cherry-picks (passing there
  too), + 1 for each item above; the orchestrator's offline runs (run_suite.sh, run_direct.py) and the bench's
  offline tests are among them. The expected-properties test (test_gate
  test_a_dropped_unmarked_property_is_no_scenario_level_mark) fails before its item on the merge's outcome
  (inconclusive, not pass); the snapshot test (test_run_snapshot
  test_a_snapshot_file_changed_after_it_was_taken_is_never_executed) before its item, where the merge, the
  comparison and the gate each executed a changed checks.py or suite_common.py once.
  test_the_merge_judges_with_the_snapshot records its snapshot change in SNAPSHOT.json (it asserted that a
  changed snapshot is executed).
- Copies of v3 (lane-s5e/reeval-S5e) and v4 (doc/v4-c11db9dd) re-evaluated at the integrated head (`run_suite.sh
  --reevaluate COPY --no-pdf --reeval-mcp-auth anonymous`, v4 with `--gate gate/policy_v4_rt_mcp.json`): against
  c11db9dd and against S8g's head every execution line of the 2,900 / 2,913 is equal in every field and check item
  (detail included) but the checks' sha256; 0 changes of a direct or MCP class, a by-design verification, a relation
  or the summary; v3 FAIL and v4 FINDINGS (exploratory, under its gate) as before. The gate on the v4 copy: 46
  problems, problem for problem S8g's (c11db9dd's 45 and "arms/mcp states no observed index_fp; the policy pins
  907ec986e19ecb04 (service_capabilities.host_document)"), its notes and counts equal (judged 637, passed 550).

## Lane V5: the v5 builder, decisive under load (Phase A: a dry v5, not frozen, not the default) (2026-10-15)

`build_manifest.py` builds manifest version 5 (`suite_v5.json` / `suite_smoke_v5.json`, smoke tag `smoke-v5`); the
frozen v1-v4 stay (1387455f rebuilds v4 byte for byte from the record cache it was built from). Phase A builds v5
dry (lane V5's scratch, content_sha256 ff7afbd7..., relations_sha256 39660cfa..., with the owner's amendments and
VERIFY-V5A's fixes and its re-verification's of 2026-10-15); run_suite.sh keeps v4 as its default until the freeze (Phase B: freeze, smoke,
default, the gate policy). What changed, by finding:

### Ground truth and the builder
- **Reverse negatives k-mer by k-mer** (the outside review's finding 8, VERIFY-D; P4-builder-truth-02).
  `Ctx.rc_negative_window` draws until no k-mer of the reverse query is in the whole record
  (`manifest_lib.reverse_kmers_in_record`, which refuses a fetched region); whole-window uniqueness let about 1 % of
  accepted windows keep one. R-rc's second window is E. coli's (its reverse k-mers are in the graph through its many
  records, so the header check discriminates; v4's mja window was absent from the graph). T-strand and the
  mutated-seed refusals (N-nearmiss-0035/-0038/-0041) accept either seed refusal (`seed '<32>': `): both follow
  from the record alone; which one appears is the index's.
- **Chance bounds** (VERIFY-D). Every global negative control (truth random, shuffled_from or random_peptide)
  states `truth.chance_bound` (the union bound over the graph's k-mers; a shuffled query also its
  composition-matched worst case; `areas.chance_bound` states the model per route), and the builder refuses one
  above 1e-4. N-negative-0018's random peptide has 12 residues with paths (6.8e-8; 10 residues met 5.6e-4).
- **index-empirical**: claims the records cannot certify (X-junction's graph-wide absence, R-junction's columns,
  which seed refusal T-strand gets, the telomere warm-up's human column, a region record's mutated seed, the
  full-support ties R-discover-0001 / -0004's limits lie above) carry
  `truth.index_empirical` {index_fp, claim}: a failure on another index asks to re-validate the expectation.
- **Flanks end at a base outside ACGT** (triage-v4 V03): truth_at / flank_for cut the flank there and state it
  (`flank_cut` {at, base}); `walk_expectations` stops a header walk's minimum where the next step would take such a
  base (`first: non_acgt`) and counts no k-mer holding one as a successor. Only T-walk-0018 changes (its IUPAC S).
- **Exact runs on headers, a tolerance on columns** (P5-v4-content-03): R-junction asks the records' headers for
  their exact runs (`areas._record_runs`: the k-mers each whole record holds) and the taxid columns for one run that
  reaches into the junction, never past it (`run_end_in` / `run_start_in`, at most 100 k-mers); the single-substitution
  resolves ask the header's exact runs and the column's cover of both sides only (pae's column holds 1,762 records).
- **Titles tested** (finding 6A, P5-v4-content-04): L-long-0006 asks path_spells on '+', L-long-0005 on '-';
  N-usecase-0005 counts both primers forward. **The re-entry** (P5-v4-content-02): N-realrep-0003/-0006 expect
  `blocked_event rejoined_seed` at 0 bp (a period of at most 40 - k bases inside the 40 bp seed).
- **Nothing drops silently** (P5-v4-content-05): the X-table1 codes not reached (2, 5) or without a reassigned
  codon (12), the N area's own-table pair and M-protein's table-12 slot are `skipped` entries.
- **gen_D's pool** (P4-builder-truth-03) leaves out run-unique and retention-dependent scenarios, those that repeat
  already, those whose status admits more than a 200, warm-ups and gate-coverage scenarios (the picks stay v4's).

### Decisive under load (items 1 and 5; DECISIVE-v5.md in the lane directory lists the 32 v4 problems)
- **/resolve**: explicit-label resolves name the host's cap, 60,000 ms (`RESOLVE_TIME_BUDGET_MS`; the service's own
  default), not 20,000 (cold first calls of 21-37 s under the load); the builder refuses a snapshot with a lower cap;
  `resolve.time_budget_above_host_max` names only a budget above the snapshot's cap.
- **The stated stop** (V05 with VERIFY-V4's amendment): R-budget-0001 resolves eco's 5 kb second half (eco's
  10 kb never completed here) under 480 ms (the owner, 2026-10-15, replacing 700), R-budget-0002 pae's 10 kb under
  420 ms: the work deadline (budget - 250 ms) inside the warm work, after each one's complete twin (R-twin-0002 /
  -0004 at 60,000 ms, the warm-up; truth.stop_design, with the batch model's stop window of warm work,
  `stop_window_warm_ms`: eco 284-2,480 ms, assuming the unpolled k-mer lookup ends within the 230 ms deadline; pae
  209-21,800 ms). A stop still depends on batch times; the owner keeps both in the deploy gate, required on
  direct-local (2026-10-15), one repeat each (a repeat decides on its own); a complete answer reads inconclusive
  (lane G1's rule). A stated-stop design, like a run-unique scenario expecting its 200, decides on one execution:
  the builder keeps it off the paired plan's warm-repeat index (`keep_off_warm_repeat` moves the nearest movable
  scenario of its area from before it to right after it; `warm_repeat_problems` refuses the rest).
- **Walks** (V06, G24, as the owner amended them on 2026-10-15): the three 1 kb header walks a cold first call left
  inconclusive in v4 (T-walk-0004, -0014, -0022; the MCP's walk preset caps the budget at 15 s, cold rates 36-58
  bp/s) walk 300 bp, the other nine keep 1 kb. T-walk-0006 / -0012 / -0018 / -0036, the gate's walks, ask the host's
  60,000 ms (refused by the MCP above its preset, by design, so direct-local's repeat 0 is the cold call) and walk
  1,000 / 500 / 5,000 / 1,000 bp, what a cold first call reached in about 30 s (VERIFY-V5A: the frontiers of v3's and
  v4's successive 15 s calls chained); each has a long twin right after it on the same window,
  T-longwalk-0001..0004, the v4 test (10,000 bp) under the host's 60,000 ms instead of v4's 15,000 (which stopped all
  32 executions of v3 and v4 by time), tagged `budget_stress` (truth.budget_stress: whether it completes depends on
  the load) and kept in the full suite. T-repeat 400 bp, eco's derived-label pair 50
  bp, T-coords 250 bp. Every window is drawn as in v4 (the class's margin), so nothing else moves.
- **Warm-ups** (tag `warm-up`, the top-level `warms` naming the target; `build_manifest.place_warm_ups` keeps each
  right before its target, and `suite_common.select_scenarios` takes it into every selection that takes the target:
  ids, areas, limit, a resume's ids): a resolve of the rows the target reads (under its header: tuples; under its
  column: rows), itself a positive check. A walk's warm-up reads what the walk must reach to decide walk_in_record:
  the seed and each arm's exact minimum plus 45 bp (`minimum_warm_up`; VERIFY-V5A: a cold multi-copy region costs
  seconds per row). R-warm-0001..0004 before T-repeat-0001..0004 (298 / 295 / 320 / 743 bp; T-repeat-0002's 295
  bp split in two resolves with their own 60,000 ms budget, `areas.SPLIT_WARM_IDS`: R-warm-0002 the seed and its
  left 46 bp, R-warm-0008 the seed's last 30 bp and its right 149 bp, since its seed rows alone took 17.6-22 s cold
  in v3 / v4), R-warm-0005 before
  T-labels-0003, R-warm-0006 (the telomere repeat under the human column) before N-hard-0022, R-warm-0007 (the seed
  and 100 bp each side under column 10710) before T-walk-0009, whose seed row took 15.7 s cold in v3; R-twin-0002 /
  -0004 before R-budget. Ids are fixed per target (`areas.WARM_IDS`). This assumes each arm runs scenarios in manifest
  order (the paired loop and the async runner do). The builder refuses a run-unique scenario expecting its 200 on the
  paired plan's warm-repeat index (`warm_repeat_problems`: T-attempts-0002 would meet its own 409), and a stated-stop
  design there too (`keep_off_warm_repeat` places both off it).
- **The resume rule of a warm-up**: after an interruption between a warm-up and its target, a resume runs the warm-up
  again right before the target on every arm that runs the target (`suite_common.paired_plan` `rewarm`, decided
  before anything of the invocation runs; `suite_common.resume_rewarms` inside a runner's scope). The paired loop
  names the warm-up with the target in that arm's step (`run_direct.py --only W,T`; the MCP runner `--ids W,T`, on
  the target's first MCP step left), the sequential mode puts it into the target's first pass; a runner applies the
  rule only within its scope (`--only`, the ids of one MCP call), so a step of the target alone never re-runs its
  warm-up. The re-run is a line with the warm-up's id and repeat (run_direct.py: repeat 0; the MCP runner: the step's
  repeat), `rerun_for_resume` true and `rerun_target`, on the line and in its raw record (run_direct.py: a raw file of
  its own); the earlier execution stays. `merge_run.py` never lets it replace its key: it goes under the target's
  `warm_up_reruns` (the first execution of the target on that arm after it, with `before_s`) and into run.json
  `reruns_for_resume`; `reevaluate.py` carries it unjudged, compare, the report (a table beside the timing) and the
  gate judge the original, so it is never a second outcome, a duplicate, a missing or an extra execution. A warm-up
  the MCP arm does not send (v6's P-warm, the long class: direct-local only, it warms the host for its target's MCP
  cold call, which the paired loop runs first) runs again on direct-local right before the target's first MCP step
  left (paired_plan `rewarm` `before_mcp`, the plan line's ninth field; `run_direct.py --rewarm-only --only W,T`
  runs the re-run alone, whatever the target has left there), and not a second time in the target's direct-local
  step; where that step comes before the first MCP step left (a warm repeat), its own re-run serves (round 2's
  verification, L1: the re-run inside the target's direct-local step came after its MCP cold call). A run whose
  snapshot's run_direct.py predates `--rewarm-only` does not run the warm-up again, and the log says so (no v6 run has
  such a snapshot; round 2's re-verification, L-a). The re-run is timing
  evidence of the target's next direct-local execution, if any.
- **T-budget-0003** (G18): T-attempts-0002's eco seed, whose seed phase took 148,723 work units in v4
  (usage.per_seed[0].work_seed), under 300,000: partial walks only. **T-labels-0004** (V02): 4,096 MB for its 6,181
  derived headers (counts hard-coded from the v4 answers, `derived_labels_v4`); T-labels-0006 drops the override and
  walks via MCP. **T-end** (P4-builder-truth-01): the preset's 15,000 ms, so the MCP walks the record-end stop.
  **T-attempts** (V04): the 200 (T-attempts-0002) and the retained 409 (T-attempts-0003) are two scenarios.

### The gate's coverage (item 4, P5-v4-content-06)
- Added (MCP-expressible, on the resolve / traverse routes; tag `gate-coverage`): T-caps-0002 / -0003, the trie's
  `max_live_paths` 2 and `max_paths` 1 forced on a 3-copy M. genitalium 31-mer whose header-constrained trie holds 3
  heads by 23 / 28 bp (`manifest_lib.trie_heads`, from the record); T-caps-0004, `max_live_paths` 1,001 above the
  service's maximum (the host walks it, the MCP refuses it: traverse.bounds_knob_above_service_maximum); T-quorum-0002,
  the permit side of `min_live_labels` (header and taxid of the record at 2); T-seeds-0002, a seed of two labels of
  which the host drops one (`seed_unsupported`, labels_supporting_total 1) with a client `seed_id` (echoed,
  `seed_id_mismatch`). The stop of a running job (traverse_stop) is the async session's area S.
- The service's discovery cap is 1,000 from the v5 start (50 when none is given; the host has none): every discovery
  names `max_labels` (at most 1,000 but the by-design pair) and `kind` on both arms, and the builder refuses another
  (`SERVICE_DISCOVER_MAX_LABELS`). R-edge-0023 (150) is answered on every arm (the source listed or in a full-support
  tie); R-edge-0025 (1,001) is refused by the MCP and answered by the host
  (resolve.discover_limit_above_cap_host_answers); R-edge-0022 (0), refused on both sides, names
  resolve.discover_limit_above_cap (the table's id for the refusal below 1, resolve.discover_limit_above_100 its alias
  for v1-v4). The service states its cap at backends.<db>.resolve.discover.max_labels of /traverse/capabilities. R-discover-0001 / -0004 ask v4's full-support tie plus one (9, 3), so the
  source must be listed.
- **SPEC request paths not exercised, by decision** (each one stays out of the manifest for its stated reason):
  - /resolve `select` (policies, `merge_overlapping`, `max_labels_per_seed`), `run_format rle`, `min_block_kmers` /
    `min_block_bp`: not expressible via MCP (traverse_resolve has no such argument; VERIFY-D); a direct-arm gap,
    covered by the engine's ResolveTest fixtures; candidates for the full suite's next increment, outside the
    MCP gate.
  - /traverse label switching (`labels.extra`, `change_cost`, `loss_budget`, `switch_on`, `max_switch_sources`),
    the quorum knobs `min_successor_labels` / `min_successor_fraction`, `max_splits_per_path`, `hairpins`, the
    frontier (`order`, `on_overflow`, beam), `branching.max_label_branches` 0 / 1, `output.detail` summary / tree,
    `output.sequences`, profile bins, continuation: the MCP's advanced strategy path is off on staging
    (preset_unavailable) and `output.*` is not overridable; the engine's synthetic fixtures (SPEC §11.3 T15, T27,
    T35) cover them.
  - `bubble_window_bp` / `tip_window_bp`: not served features (only 0 is accepted); T-edge-0015 tests tip's refusal.
  - The attempt routes (GET /traverse/attempt, POST /traverse/cancel, `expect_server_instance`): not expressible via
    MCP beyond traverse_stop (area S); T-attempts-0001..0003 cover the attempt fields on the direct arms.
  - /pattern `max_labels_per_anchor`, `max_occurrences_per_label`, `occurrences: false`, a `max_labels` cut,
    `allow_unbudgeted_annotation`, `stop_at_threshold` beyond P-modes-0005, Q's motif `tested_contexts`,
    `predicate.vacuous`, `predicate_too_large`: the /pattern route is not in the resolve / traverse gate; several
    are not expressible via MCP (PATTERN_NOT_EXPRESSIBLE); left to the full suite's next increment.

### mcp.exposes, the MCP classes in the manifest and the declared gate entries (items 9 and 6b)
- `build_manifest.py --exposure-from RUN --exposure-out FILE` writes a run's exposure baseline
  (truth/mcp_exposure/v4-20261014T061053Z.json: the v4 run's MCP lines re-evaluated at 1387455f); the build fills
  `mcp.exposes` with the checks decided (pass, fail, finding) on EVERY recorded MCP execution of the scenario whose
  status passed, only where its MCP request is the baseline's (traverse_resolve's time budget aside) and restricted to
  its checks; the manifest states the baseline (`mcp_exposure_baseline`). 512 v5 scenarios carry a list.
- What the MCP cannot show is classed in the manifest (the owner, 2026-10-15), validated by the builder (a check the
  scenario expects, a reason, a case of the table, a scenario the MCP runs): `mcp.not_implemented_checks`
  [{check, reason}] for label_dict_contains of T-labels-0003..0006 (not implemented yet via MCP, the walk tested via
  MCP); C-caps-0001 `mcp.not_implemented` as a whole (the MCP's capabilities tool answers the service's view);
  `mcp.by_design_checks` [{check, case}] for a property check the service answers otherwise by design.
  `gate_policy.py` derived the MCP exemptions (rules mcp_not_implemented_check, mcp_by_design_check; removed in
  integration S8i stage 2: gate.py applies the marks); gate/declared_v5_rt.json declares nothing (the v4 file stays
  v4's); every v4 suite error of its §2 is fixed. A missing MCP feature is marked at the narrowest scope, and when
  the property-level marks would cover every property check of a scenario that scope is the scenario
  (`build_manifest.scenario_level_marks`, VERIFY-G (c): such a scenario would pass via MCP on its generic checks
  alone): `mcp.not_implemented` with the marks' reasons ("service-side by design: ..."; the manifest has no
  whole-scenario by-design form), and the builder refuses a property-level list that covers them all. T-edge-0008,
  whose only property check, answer_fields, is service-side by design (traverse.duplicate_seeds), was that case
  until integration S8i stage 2 modelled it through the table (run via MCP).
- Every `mcp.by_design` id is a case of the bundled service_by_design.json, else the build fails (a case the table
  lacks explains nothing). The table is lane F-gate's 422644b7 (sha256 39f87b83..., the renamed discovery cases),
  bundled here with its v4 policy so the integration meets the same blobs.
- `gate_policy.py` excludes a scenario tagged `budget_stress` on every arm (rule budget_stress, the manifest's
  truth.budget_stress.reason, printed by describe()): no declared entry; the owner approved it for T-longwalk-0001..0004
  (2026-10-15).

### Validation (offline, Phase A)
- The dry v5 against v4 (scenario by scenario; content ff7afbd7..., relations 39660cfa...): 700 scenarios (19
  added: R-warm-0001..0008, T-longwalk-0001..0004, T-attempts-0003, T-caps-0002..0004, T-quorum-0002, T-seeds-0002,
  R-edge-0025; none dropped); 175 changed beyond
  mcp.exposes and the smoke tag, each to one finding above; 106 unchanged; 400 only gained mcp.exposes or the smoke
  tag; gen_D's picks are v4's.
- Every truth re-checked against the cached records (633 windows, 620 flanks, all ACGT, 15 copies; 46 walk minima,
  12 header runs, 2 trie head counts, 34 chance bounds, 3 reverse negatives recomputed; 10 warm-ups placed, holding
  their seed and minima, T-repeat-0002's split pair included; 14 discoveries explicit): 0 problems.
- The v4 answers judged by v5's expectations where the request is v4's (554 scenarios): every new check passes on
  every arm and repeat; N-usecase-0005's MCP executions turn inconclusive -> pass (a property the MCP exposes);
  v4's T-attempts-0002 repeat-1 409 is T-attempts-0003's answer in v5 (pass there); R-edge-0023's v4 MCP answers
  are the service's old refusal above 100, which its cap of 1,000 retires.
- The v5 gate policy derives from the dry build (no declared entry): 241 selected, 37 excluded (by rule,
  one exclusion may carry several: manifest_not_implemented 25, table_not_expressible 27, budget_stress 4,
  advanced_path 1), MCP 204 scenarios / 414
  executions, direct-local 237 / 273, 4 property checks not required on MCP (derived); REL-same_answer-0003 and
  REL-prefix_equals-0001..0004 required.

## Suite integration S8i, stage 2: lane V5 merged, and the rules that join it (2026-10-15)
Lane V5 (gr/suite-v5 at 0c165be5, 45 commits on 1387455f) cherry-picked onto fafd20ba in order, with its messages and
authorship: 44 commits; 5f941931 (lane S8g's table commit 422644b7, already here) came out empty and was skipped. Only
this README conflicted (V5's section kept after S8i's, without its restated older S7i line). The items below join
what VERIFY-V5A's final round (the trial merge) and the check of async's runner found between V5 and S8c / S8g.

### One mechanism for the MCP marks
gate_policy.py derived a `not_required` MCP entry for every mcp.not_implemented_checks / mcp.by_design_checks entry,
while gate.py applies the same marks itself and sets aside only a marked check the MCP reads not_exposed
(suite_common.checks_marked_aside). Merged, the exemption silenced a marked check's missing result, failure or
inconclusive on MCP (test_gate's test_a_marked_check_without_a_result_is_a_missing_one failed on the merge). The
policy derives no property exemption now; the marks live in the manifest and gate.py applies them. The v5 policy has
no property check not required.

### T-edge-0008: service-side by design through the table, run via MCP
The owner classes T-edge-0008's merged duplicate seed as service-side by design. Lane V5 marked the scenario
mcp.not_implemented, which dropped v4's MCP coverage of the service's seed merge. The merge makes two readings on the
MCP arm: invariants fails "1 results for 2 seeds (1 distinct)" (which traverse.duplicate_seeds already explained) and
answer_fields reads not_exposed, "field not in the answer: results.1.duplicate" (checks.evaluate's rule for a missing
field on the MCP arm). An explains_checks entry explained only failing items, so the table now states
{"check": "answer_fields", "detail_matches": "field not in the answer: results\\.1\\.duplicate"}, and the rule is:
- **a {check, detail_matches} entry can explain a not_exposed reading** (its detail matched completely: a field the
  service's path removed from the answer);
- **a bare check name never explains a not_exposed reading** (it explains every failing item of that check, as
  before).
One definition, a shared rule of suite_common: item_explained (with EXPLAINABLE_OUTCOMES, fail and not_exposed), used
by explained_outcome (E: an execution whose only failing items E holds reads service_by_design; any other failure
still fails) and explained_items (the line's explained_items), which compare.py names, reevaluate.py records,
merge_run.py takes out of a line before the checks' rule (a recorded explained item of either outcome) and async's
runner calls through shared_rule. **The gate takes the explanation from the table, not from the line**
(VERIFY-G stage 2: a line's own claim let a not_exposed answer_fields on another path, an r1 claim beside a
legitimate r0, or a direct-arm line carrying the case's fields pass): gate.table_explanations reads compare's
table-derived explained items and case (mcp_agreement explained_items / by_design_notes), and the property loop skips
an item of either outcome only on the MCP arm, only where the table explains it under a case the scenario names; a
line whose recorded explained_items (or explained_checks, where it records no items: compare's
explained_checks_recorded) differ from the table's is an agreement problem. The agreement step requires the pair
verified, as for every answered case. Any other failure still fails. T-edge-0008 keeps mcp.by_design
(from its request) and carries no mark, so the MCP arm runs it again; the builder's scenario_level_marks and the
SCENARIO_LEVEL_MARK backstop stay for marks. The table's sha256 is
d4b8636c748438ddbed1073a22379bcf39ae148ad90e879575dbdac38cd64f44 (133 cases); the v4 policy is regenerated for it
(only by_design_table.sha256 and generated_utc differ), and v4's own declared entry for T-edge-0008 stays v4's. v4's
recorded T-edge-0008 MCP executions, re-evaluated under the v5 scenario with this table: service_by_design
(traverse.duplicate_seeds, E = invariants and answer_fields, items 1 and 2), the pair verified against direct-local;
a policy requiring T-edge-0008 on both arms raises no problem for it (with E = invariants alone, as a runner without
item_explained records it, the gate now names the line's other explanation as an agreement problem on each repeat).

### The loop stops at once on the MCP runner's exit 4 and 5
run_suite.sh's mcp_after stopped at once only on exit 2 / 3 / 130. Exit 4 (the runner refused the by-design table) and
5 (it refused the run's snapshot) read as unrecorded executions: two more MCP steps that recorded nothing and the
direct arms' steps between them ran before the loop stopped as mcp_errors. Both now stop it after the first step
(`loop.reason` mcp_exit_4 / mcp_exit_5), in the paired and the sequential mode (tests/fake_mcp.py: modes exit4,
exit5). The bench (bench/run_bench.py, a port made before this fix) stops alike on 2, 3, 4 and 5 (its session's
reason mcp_exit_N) from suite_common.MCP_START_REFUSALS, which a test holds to mcp_after's case list (SD-029).

### The gate checks which checks computed the set's comparison
A compare.json computed while the run's snapshot was changed (the refusal's built-in fallback, or other checks), the
snapshot restored before the gate, passed: every other identity check reads the restored files; so did a
compare.reevaluated.json computed while the re-evaluation's checks were changed. gate.compared_with_checks_problems
applies one rule to both sets: the comparison's judged_with names the checks that judged the set as its source (the
run's snapshot for a recorded set; the re-evaluation's checks for a re-evaluated one) and states the checks.py and
suite_common.py sha256 recorded for them (SNAPSHOT.json; reevaluated.json's checks); otherwise an identity problem asks
to recompute the comparison. A recorded set without a snapshot is judged as before. A re-evaluated set is evidence by
its checks' hashes, not its directory: compared with this tree's checks (its re-evaluation's directory gone), it
stands when their sha256 are the ones reevaluated.json records; and the files the comparison and the gate's property
reading read now (the re-evaluation's directory, or this tree) must be those, else "re-evaluate the set" (a directory
edited after its comparison no longer steers the property reading). The advice is the step that clears the problem.
The comparison also states executions_sha256, the sha256 of the merged lines it read (executions.jsonl,
executions.reevaluated.jsonl); the gate blocks on identity when the lines it judges are other ones (edited after the
comparison: "recompute the comparison").

### A broken snapshot record or a missing executed file is refused, never a crash
snapshot_refusal read SNAPSHOT.json with read_json, so an empty, torn or unreadable record raised (async's runner
exited 1 with a traceback instead of refusing with exit 5; reevaluate crashed; the merge swallowed it without its
warning), and it skipped a recorded file it could not read (a deleted suite_common.py loaded another one).
suite_common.snapshot_record reads the record without raising and names the cause (cannot be read, empty, no JSON, no
files); snapshot_refusal refuses for it and for a recorded executed file that is missing or unreadable, so
load_checks_bound and its callers (merge, compare, the gate's property reading, reevaluate) treat both as a changed
file. snapshot_state never raises: an unreadable record is no snapshot and changed (`SNAPSHOT.json`), and the
orchestrator's snapshot step stops the run (exit 2) instead of taking a new snapshot over it. The bench reads it the
same way: bench/analyze.py's comparison refuses a run whose record is empty, torn or unreadable through
snapshot_refusal (the record only), naming the cause per side, and bench/run_bench.py's resume says the record is
unreadable and why ("has no snapshot" only where there is none).

### The dry v5 and its policy on the merged tree (not frozen)
- The builder on the merged tree (lane V5's command: `SUITE_TRUTH_CACHE=<the record cache v4 was built from>
  build_manifest.py --no-fetch --out ... --smoke-out ...`): content
  72181c4a275e6bb94a905e6d43eaad4e179a6b436dd5e00175e3936ea3aa1875, relations 39660cfa... (unchanged), 700 scenarios,
  742 executions per arm, 51 MCP not implemented, 512 with mcp.exposes; the smoke selection abc0596f... (14). Against
  VERIFY-V5A's ff7afbd7... only T-edge-0008's mcp changed (its tool and arguments back, exposes [budget_respected,
  invariants, status], no not_implemented). The offline truth re-check: 0 problems.
- gate/declared_v5_rt.json names that content (it declares nothing). The V5 branch committed no dry manifest and no
  v5 policy; the policy is generated at the freeze, on this tree.
- `gate_policy.py <dry v5> --name v5-rt-mcp --arms mcp,direct-local --declared gate/declared_v5_rt.json`: 241
  selected, 36 excluded (manifest_not_implemented 24, table_not_expressible 27, budget_stress 4, advanced_path 1),
  MCP 205 scenarios / 416 executions, direct-local 237 / 273, 5 relations required (21 excluded), no property check
  not required. Against VERIFY-V5A's dry policy: T-edge-0008 runs on MCP (+1 scenario, +2 executions), no
  exemptions.

### For async: what the runner adds for stage 2
Re-vendor suite_common.py and service_by_design.json (d4b8636c...); checks.py is S8c's. Besides the S8i list above:
- **Explained items** (T-edge-0008): the runner calls the shared rule, never a copy: shared_rule("explained_outcome")
  (checks, the matched cases' entries, its outcome rule) for the outcome and E, shared_rule("explained_items") for the
  line's explained_items (explained_item_indices), shared_rule("item_explained") (entry, item) wherever it picks the
  recorded case (EXPLAINABLE_OUTCOMES: fail, not_exposed). A {check, detail_matches} entry can explain a not_exposed
  reading; a bare check name never does. explained_items are indices over the FULL recorded checks of the line (in
  the order recorded, before anything is taken out), while the outcome rule judges shared_rule("unmarked_checks")
  (scenario, "mcp", checks), as reevaluate does. The gate takes an execution's explanation from the table (compare's
  table-derived explained items and case, MCP arm only), never from the line: a line whose recorded explained_items,
  or explained_checks where it records no items, are not the table's is an agreement problem ("a runner's explanation
  is not trusted"), so a runner without the shared rule (E = invariants alone on T-edge-0008) blocks the gate.
- **T-edge-0008 runs via MCP** (traverse_fetch with its two identical seeds; the service merges them).
- **Warm-ups** (lane V5), shared rules, never a copy: shared_rule("select_scenarios") (manifest, areas, limit, ids)
  pulls a target's warm-ups (`warms`) into every selection that takes the target, each arm in manifest order;
  shared_rule("resume_rewarms") (scope, done) is the resume rule, applied only to the ids of one call (`--ids W,T`),
  never to a pulled warm-up; the re-run is a line with the warm-up's id and the step's repeat and the fields
  shared_rule("RERUN_FIELDS") names (`rerun_for_resume` true, `rerun_target`), on the line and in its own raw record
  named by shared_rule("rerun_raw_name") (scenario, arm, repeat, started); the earlier execution stays
  (is_rerun / judged_lines read it as timing evidence).
- **The snapshot refusal** (the S8i loader item, with this stage's never-raising refusal): before the snapshot's
  suite_common.py is first executed, pinned_common().snapshot_refusal({checks_dir}) refuses an unreadable, empty or
  torn SNAPSHOT.json and a recorded executed file that is missing, and the invocation exits 5 with the message
  (new: snapshot_record). run_suite.sh stops the loop at once on 4 and 5.
- The deadline's origin terms: the runner states them as the service rounds them (0.1 ms; bound_s 1 ms); the 2 ms
  tolerance covers both (wording corrected in S8g's section).

### Validation of stage 2 (offline)
- 585 tests pass (none skipped, with SUITE_TRUTH_CACHE, the record cache v4 was built from, so lane V5's builder
  tests run): 527 at fafd20ba + 53 of lane V5 (441 -> 495 on its branch, one of them 422644b7's table test, already
  here) = 580 after the cherry-picks (one failing there: the marks test above), + 2 (the explained not_exposed
  item: the rule, the gate) + 1 (exit 4 / 5) + 1 (the comparison's checks) + 1 (the broken snapshot record); V5's
  policy test is rewritten in place (no exemption derived).
  Each new test fails on its parent commit.
- Copies of v3 (lane-s5e/reeval-S5e) and v4 (doc/v4-c11db9dd) re-evaluated at the merged head, against fafd20ba's
  re-evaluation: the only change is T-edge-0008's two MCP lines on each (explained_checks invariants ->
  invariants, answer_fields; explained_items [1] -> [1, 2]); outcome, case and every check item unchanged, no class,
  verification or relation moves; T-edge-0008's checks_not_exposed 2 -> 0 (summary 854 -> 852 on v3, 828 -> 826 on
  v4). v3 FAIL; the v4 gate 46 problems, the same as fafd20ba's in order, its notes and counts equal.
- The dry v5's 554 scenarios whose request is v4's on a fresh v4 copy (2,373 executions), against the v4 copy re-
  evaluated with the v4 manifest: VERIFY-V5A's moves only (path_spells, blocked_event and count_relation pass;
  N-usecase-0005's MCP inconclusive -> pass; R-edge-0022 / -0023's MCP service_by_design -> fail under v4's cap of 100;
  T-attempts-0002 r1 is T-attempts-0003's 409).
- After the follow-ups (the shared explanation rule, the re-evaluated set's comparison check, the bench's snapshot
  record): 589 tests pass (585 + 1 + 1 + 2, each failing on its parent commit); the v3 / v4 re-evaluations equal the
  ones above line for line (every field and check item; the comparisons' scenarios and summaries byte-equal), the v4
  gate the same 46 problems, the v5-on-v4 moves the same.
- After the verification's follow-ups (the gate's table-derived explanation, a re-evaluated set's hashes, the lines'
  sha256): 592 tests pass (589 + 1 + 1 + 1, each failing on its parent commit; the failure-only control of the
  explained-item test now reads two agreement problems, the line's E not being the table's, and no property problem);
  the v3 / v4 re-evaluations equal the ones above line for line (compare gains executions_sha256), the v4 gate the same
  46 problems, the v5-on-v4 moves the same, and the dry v5 rebuilt is 72181c4a... again.

## The frozen v5 (2026-10-15)

The owner approved the freeze of v5 on 2026-10-15 at content 72181c4a..., the dry v5 built on the merged tree of
integration S8i stage 2 and verified (VERIFY-V5A "Stage 2", VERIFY-G "Stage 2" and its follow-ups).
- `manifest/suite_v5.json` and `manifest/suite_smoke_v5.json`, built on this tree with `SUITE_TRUTH_CACHE=<the record
  cache v4 was built from> build_manifest.py --no-fetch` from the capabilities snapshot of release refseq97-6f945db4
  (truth/staging_capabilities/, each document's sha256 in the manifest's index): 700 scenarios, 742 executions per
  arm, content_sha256 72181c4a275e6bb94a905e6d43eaad4e179a6b436dd5e00175e3936ea3aa1875, relations
  39660cfa0fb98cc47ac923fad3aa33419a84a3d92d662154cf2d4164752f4a86; smoke 14 scenarios (14 executions), content_sha256
  abc0596fef7977dc0fa7d51bf7b8d0194e6427ed726f0a9acf12605116546e8d. Both equal the verified dry build but for
  created_utc. The by-design table they name is d4b8636c... (133 cases).
- `gate/policy_v5_rt_mcp.json` (sha256 23ba28f5...), generated on this tree with `gate_policy.py
  manifest/suite_v5.json --name v5-rt-mcp --arms mcp,direct-local --declared gate/declared_v5_rt.json --out
  gate/policy_v5_rt_mcp.json`: 241 scenarios selected; 36 excluded (manifest_not_implemented 24,
  table_not_expressible 27, budget_stress 4, advanced_path 1; one exclusion may carry several rules); required on
  MCP 205 scenarios / 416 executions, on direct-local 237 / 273; relations REL-same_answer-0003 (both arms) and
  REL-prefix_equals-0001..0004 (direct-local) required, 21 excluded; no property check not required (gate.py
  applies the manifest's MCP marks). gate/declared_v5_rt.json names the frozen content and declares nothing.
- run_suite.sh takes v5 by default (`--smoke`: suite_smoke_v5.json), as the frozen v4 was made its default;
  build_manifest.py refuses to overwrite either; v1-v4 stay for reproducing earlier runs. test_manifest_lib's
  frozen-manifest test covers v5, and test_gate pins the policy (regenerated exactly, its figures, the contents).
- Before the first v5 run (VERIFY-V5A phase B): the service at async's fb66f1dd or later (the discovery cap 1,000,
  for R-edge-0022 / -0023 / -0025); async's runner re-vendored (suite_common, the table d4b8636c...) and calling the
  shared rules of the runner list above (the explanation, the warm-ups, the snapshot refusal); R-budget-0001 / -0002
  block whenever no stop happens (their windows rest on eco's 230 ms lookup).

## Small synthetic-graph oracle (lane O, phase 1, 2026-10-15)

The staging suite checks that known examples are found; `oracle/` checks that **every** answer is found and
every count is exact. On `/pattern` no result may be duplicated. On `/traverse`, since the verification's fixes
(H1), no segment, leaf, run or event may be listed twice. (GPT's review of 2026-10-14, "What is missing" items
2-3 and the P1 row "independent small-graph oracle".) Design: `ORACLE-DESIGN.md` in the lane directory
(`~/.cache/metagraph-traverse-plans/2026-10-15/lane-oracle/`).

### Files
- `oracle/fixtures.py`: a few tiny records with their features built in (a SNP bubble, a shared element and
  its branches, a reverse-complemented block, a period-7 repeat and a homopolymer, even palindromes, planted
  neighbours at one to three substitutions, a record split by an N, records shorter than k first and last in
  their columns, columns of several records), the same bytes on every run.
- `oracle/oracle.py`: the brute-force enumerator: every k-mer, offset and orientation tested position by
  position (own IUPAC table, reverse complement, NCBI codes 1, 2, 4, 6, 11, 27, 31). It imports `math` only
  (`tests/test_oracle_core.py` asserts it): nothing of the engine, `checks.py`, `seqtools.py` or `truth.py`.
- `oracle/compare.py`: the checks (families `shape`, `notes`, `counts`, `results`, `release`, `labels`,
  `stop`, `resolve`; outcomes pass / fail / inconclusive).
- `oracle/mutations.py`, `oracle/render.py`: one-property corruptions with the family that must catch each,
  and answers rendered from the oracle for the offline triplets.
- `oracle/requests_gen.py`: the deterministic requests (`--seed`).
- `oracle/indexes.py`, `oracle/server.py`, `oracle/report.py`, `oracle/run_oracle.py`: the build, the
  localhost server, the report and the command.

### Run

    python3 metagraph/scripts/suite/oracle/run_oracle.py --metagraph metagraph/build/metagraph --out DIR

It builds nine configurations under `DIR/build/` with the engine's CLI (`basic_k11` with the mask and a
`row_diff_brwt_coord` annotation with coordinates and `.seqs`; the same graph without its mask at
`--pattern-max-checked-entries` 50 and 0; `canonical_k11` and `primary_k11` with column annotations;
`basic_k15`; and at an even k, where k-mers equal to their reverse complement exist, `basic_k10`,
`canonical_k10`, `primary_k10`). For each it starts `server_query` on a free 127.0.0.1 port (floor 0 bits, ample caps), sends about
100 requests (60 pattern requests of up to 9 patterns, 12 budget probes, 30 `/resolve` on BASIC), compares
every answer, applies every applicable corruption to every answer that passed, and stops the server. Output:
`DIR/oracle_report.json`, `DIR/oracle_report.md` (scope, per-config outcomes, every disagreement with the
request, the oracle's set and the engine's answer, the corruptions the comparison missed),
`DIR/raw/<config>/<id>.json.gz` (each request and answer), `DIR/logs/`. `DIR` must be new or empty; nothing is
deleted. Exit 0 when every check passed, 1 when one failed, 2 on a set-up error. About four minutes.

Offline (no server): `python3 -m unittest tests.test_oracle_core tests.test_oracle_compare tests.test_oracle_run`
(about 10 s): the oracle's definitions on hand-checked cases (several are the SPEC's own examples), the triplets
on rendered answers in the three graph modes, the review's probes, the generator and the report.

### Readings the oracle takes where the SPEC leaves a choice
- `seq_id` (SPEC 14.3 "the record's index in its column (the `.seqs`: ... FASTA order)"): the index among the
  column's records of at least k bases. The `.seqs` builder skips a shorter record
  (`records_shorter_than_k: "not_indexed"`), so a record shorter than k shifts the `seq_id` of every later
  record of its column against its FASTA index, and its header is no `/resolve` label (400).
- Label order "(contexts desc, column asc)": the column's name. The fixture's column ids are not in name order,
  so an order by id would show as a disagreement.
- `/resolve` trace support: a set of live coordinates, each continued by +1; a k-mer that continues none
  starts a new run. A column label's coordinates concatenate the column's indexed records, so a chain can cross
  a record boundary (`trace_record_boundaries`).
- U of a graph without its mask (SPEC 7.4): the masked count plus the contexts in source dummies (`'$' * j +
  x[:k - j]` for every node x no k-mer enters), a dummy holding the oriented pattern where its `$` run ends at
  or before the skipped leading N run (P's trailing run for `reverse`; a window of N only keeps its last
  position; none for a peptide). With substitutions the check of few unchecked candidates makes levels exact
  one by one (21.6), so a partition's `upper` is its truth plus the dummies of the levels still in `bounds`.

### Validation of phase 1 (2026-10-15)

*Correction (verification of 2026-10-15, M1).* These runs used the main tree's binary
`metagraph/build/metagraph_DNA`, linked 2026-10-13 16:14. It is an engine older than 81c620ef: it predates 14
engine commits. The earlier wording "against the 81c620ef build" was wrong. The verification re-ran every request
against the 81c620ef Release build (worktree `eng-int`, HEAD 81c620efd, `src/` clean) and got the same answers,
byte for byte apart from timing, so the numbers below stand. Every report now records the binary's sha256 and its
build tree's HEAD.

- Seed 1 at the default sizes (`runs/final2_seed1` in the lane directory) and seeds 2-6 at twice the size
  (`runs/final2_breadth_seed2` … `6`): 8,528 requests, 32,988 answers over the nine configurations; 30,116
  answers equal the oracle (complete sets, exact counts, order, labels and placements, `/resolve` runs and
  candidates), 2,872 budget-probe answers are inconclusive on completeness with every invariant passing, 0
  disagree. 314,662 corruptions of the passing answers, of the stated mutation set (34 kinds), every one caught.
  That is not every possible one-property corruption: the verification found corruption classes outside the
  set that the comparison then missed. They were closed afterwards (see "After the verification").
- No candidate engine defect in phase 1's scope. The disagreements of the development runs were all the oracle's
  and are fixed with regressions: the two readings above (`seq_id`, trace); U per level with substitutions; U
  of the orientation a `max_steps` stop left complete beside an `at_least` total.
- The checker gaps the runs published and that are now closed: doubled `bounds` counts (before the U model);
  `retrieval_complete` on an `all_or_count` withheld on U; an offset moved within a PRIMARY `partial` cut.

### Phase 2: `/traverse` (2026-10-15)

`oracle/traverse_oracle.py` computes the walks of a seed from the records' text (SPEC-labeled-traversal-core
6.9-6.10): the walk rule (edges as (k+1)-mers, canonical with their orientation on CANONICAL/PRIMARY;
`rejoined_seed`, `edge_reuse_rc`, `edge_reuse`; hairpins), annotate's exhaustive trie with the labels at every
node, and a reference walker of 6.3-6.4 under `on_reconverge: keep`. The walker models `forbid` and `constant`
costs, `switch_on`, the branch allowance, the quorum, `min_live_labels`, the split limit and trace support. Its
claims are the label ends. A second, simpler computation of E (single-label maximal walks) must equal them
under exhaustive `forbid`. `oracle/traverse_fixtures.py` holds one record set per rule, each with a permitting
and a forbidding setting, served as `trav_basic_k11`, `trav_canonical_k11` and `trav_primary_k11`. The
phase-1 configurations get random seeds.

`traverse_compare.py` compares three kinds of arm:
- **Complete arm:** held to the oracle as complete sets (walks, leaves, claims with branches and loss, end
  labels, label end events, node label sets, blocked and hairpin events, seed and dropped labels).
- **Truncated arm** (cap probes): judged up to `complete_to_bp` only.
- **Merge arm:** judged on its invariants: soundness, the `scope` limitation, and the label on every node of a
  run's displayed support.

`tree` and `summary` are compared against `full`. `run_oracle.py --routes traverse` runs this alone. Offline:
`python3 -m unittest tests.test_oracle_traverse`, which includes the triplets on 53 recorded engine answers
in `tests/data/`.

Readings and notes:
- **End reasons** (SPEC 13, implementation note I3-2): a run's `end_reason` is the served enum. `switched` and
  `superseded` are `label_lost`; `minority`, `below_min_labels` and `split_limit` are `branch`; a skipped
  hairpin is `dead_end`. The `label_end` event states 6.7's reason as its qualifier, and a switched source has
  none. The oracle compares runs against the served enum and events against 6.7.
- **Two hairpin readings, fitted to the engine at SPEC gaps** (both in SPEC-NOTES.md in the lane directory):
  - *A skipped hairpin before `rejoined_seed`.* 6.6 ranks `rejoined_seed > edge_reuse_rc > edge_reuse` and does
    not rank `hairpin`. At even k on PRIMARY, a step into a seed node that is a palindromic k-mer is served as a
    hairpin. The oracle decides a skipped hairpin first, on every graph mode.
  - *The branch count on a followed hairpin step.* 13 I3 note 11 ("never counts toward ambiguity, so the split
    it causes is a divergence") can be read as "the hairpin entry never counts a branch", against 6.4's "every
    derived entry gets branches + 1". The oracle decides ambiguity over the non-hairpin successors. When a
    source is ambiguous, every entry derived from it gets branches + 1, the hairpin's included. The first of
    two disagreement rounds was SPEC-based (88ac1337); this second one was fitted to the engine (2618936f).
- **Radius.** `bounds.max_extension_bp: 0` is refused (400, range [1, …]), which section 5 does not state.
- **Trace with merging** is refused (400), as SPEC 13 item 9 states.

#### Validation of phase 2 (2026-10-15; the build: see the correction under phase 1's validation)
- Six runs on the code of `2618936f`: seed 1 at the default sizes (`runs/p2final3_seed1` in the lane directory)
  and seeds 2-6 at twice the size (`runs/p2final3_breadth_seed2` … `6`), every route, all twelve configurations.
  - 16,058 requests: 7,128 `/pattern`, 1,400 `/resolve`, 7,530 `/traverse`.
  - 39,441 answers equal the oracle; 4,172 are inconclusive (budget and cap probes, `merge` runs) with every
    invariant passing; 0 disagree.
  - `/traverse` alone: 9,305 seed results and 1,320 tree/summary comparisons.
  - 405,728 corruptions of the stated mutation set, 91,066 of them of `/traverse` answers, every one caught.
    Not every possible corruption: the verification's injections found classes the set lacked (below).
- **No candidate engine defect in phase 2's scope.** The disagreements of the development and first final runs
  were all the oracle's, each fixed with a regression test:
  - followed hairpin steps were not recorded as events;
  - "would carry" under trace must follow the coordinate chain;
  - the `merge` invariant now checks displayed support (§7.1);
  - a followed hairpin step never counts toward ambiguity (SPEC 13, I3 note 11), but an ambiguous label counts the
    branch on that step too.
- The suite's offline tests: 512, all passing (`python3 -m unittest discover -s tests`, about 6 minutes).

### After the verification of 2026-10-15

The verification (Opus; the lane directory's `verify/VERIFY-ORACLE.md`) found the oracle independent and its
"0 disagreements" reproducible on the real 81c620ef build. It also found five places where a wrong engine answer
would pass, and gaps in coverage and wording. Fixed, each with a regression test that fails without it:
- **H1:** `/traverse` is compared with multiplicity and order (a run, leaf, event or subtree listed twice;
  events sorted by `at_bp`; split branches in character order).
- **H2:** label ends under a finite cost follow 6.7 with 13 I3 note 3. The re-minimisation after a branch
  exclusion honours `switch_on` and the stay. Extras unreachable within the budget are a predicted 400. The
  probes T1-T4 are pinned with their 81c620ef answers.
- **H3:** on a graph without its mask, a completed count is exact or bounds; a withhold on U needs U >
  `max_contexts`, with the `threshold_upper_bound` note exactly when `lower` ≤ `max_contexts`.
- **H4:** newly compared: `needed_budget`, every switch's source, target and cost, `row` on CANONICAL and
  PRIMARY, and the PRIMARY `partial` cut against its ample twin request (M5).
- **H5:** a merge arm that merged nothing is held to the keep walker strictly; a merge states its
  `reconverge` events.
- **M1:** provenance (the binary's sha256 and build HEAD in every report).
- **M2:** `/resolve` trace queries that tell trace from presence in every seed.
- **M3:** varied constant-cost `/traverse` requests in every seed.
- **L1, L2:** sticky later patterns state unknown counts and zero work; the note rule above.

Still judged weakly or not at all:
- **L3:** the `lower` of a bounds count, held only to ≤ truth.
- **L4:** the column-coordinate layout of records holding an N, untested; no fixture chain crosses records.
- **Merge arms that did merge** (a handful per run), after their first merge: the united history's label
  states, run ends, events and needed budgets, and which walks it continues, stay invariant-only. A tie for
  the majority parent is not judged either. See "Round 2 of the verification" below.
- **The row of a CANONICAL k-mer listed without its reverse complement**, when that k-mer is not the
  canonical one, is only bounded: `row < node − 1`. Its row is the reverse complement's node − 1. The oracle
  knows the colex order but not absolute node ids (dummy edges shift them), so it can predict that node only
  when the reverse complement is listed too.

#### Validation of the verification round (2026-10-15, against 81c620ef)
- Six runs on the code of `f5a5e38d` against the 81c620ef Release build. The binary was the `eng-int` worktree's
  `metagraph_DNA`, sha256 `63b22686…`, HEAD 81c620efd, `src/` clean; each report records this.
  - Seed 1 ran at the default sizes (`runs/v3final_seed1` in the lane directory).
  - Seeds 2-6 ran at twice the size (`runs/v3final_breadth_seed2` … `6`).
  - Every route and all twelve configurations.
- 16,142 requests: 7,128 `/pattern`, 1,484 `/resolve`, 7,530 `/traverse`.
- 43,961 answers:
  - 40,469 equal the oracle;
  - 3,488 are inconclusive (budget and cap probes, merge arms that merged), with every invariant passing;
  - **4 disagree, and all 4 are oracle errors in the invariants of merged arms**, described below.
- 437,309 corruptions of the stated mutation set, 122,059 of them of `/traverse` answers, every one caught.
- Coverage per seed of what the verification found missing:
  - 14-15 `/resolve` trace queries whose trace answer differs from k-mer presence;
  - 120-246 constant-cost `/traverse` requests over 25-33 (value, budget, `switch_on`) settings, with
    37-81 switches, 9-16 `superseded`, 8-22 `label_lost` and 3 `loss_budget` label ends;
  - 5-13 predicted refusals of unreachable extras;
  - 13-23 PRIMARY twins;
  - 9-26 withholds on graphs without their mask;
  - 96-217 merge arms that merged nothing, held to the keep walker;
  - 3-5 merged arms, each with its `reconverge` event.
  - The three `loss_budget` ends per seed all come from the fixed rule case. The random constant-cost
    requests have not yet produced one.
- **No candidate engine defect.** Per the coordinator's rule, the 4 disagreements are classified, not fitted.
  The oracle code is unchanged.
  - *A displayed walk across a merge is not a label's walk.* Three answers:
    - `T-canonical_k10-0051` left (seed 2);
    - `T-canonical_k10-0046` right (seed 4);
    - `T-canonical_k11-0033` left (seed 6).

    `traverse.*.sound` requires every displayed walk of a merged arm to be a walk of the keep walker. SPEC 7.1
    says otherwise: "The displayed walk — a path's `segments` chain and spelled bases, its `continuation`, and
    `route_bp` (what of that spelling a label does not carry …) — follows first parents".

    In each case the leaf's only label entered the merge through a non-first parent, and the leaf's `route_bp`
    is that merge. Each extra walk equals a keep walk from the merge on, and that keep walk reaches the merge
    node through the other parent (`runs/v3final_triage/triage_sound.py` in the lane directory).

    *Correction:* the script's "graph walk" test (`walk_node(...) is not None`) holds for every string. The
    extra walks are graph walks by the verifier's real check (every node a k-mer of the index, no edge
    twice), which is now the oracle's. The relaxation first proposed here would also have accepted a display
    through the minority parent, and a label said to enter through a parent that never carried it. The
    verifier's exact version replaced it (below).
  - *A switch-entered run's first node.* One answer, `T-canonical_k10-0076` right (seed 2).
    `traverse.*.displayed_support` requires a run's label on nodes `from_bp … to_bp`. Two runs `[22, 23)`,
    entered by a switch and closed by the merge at 23, lack their label on node 22. That is the switch
    source's node.

    SPEC 6.2: "Step j adds outward base index j−1 … Runs are half-open `[from_bp, to_bp)` over outward base
    indices". So a run covers the nodes entered by steps `from_bp + 1 … to_bp`; 7.4's `label_sets` say the
    same in words. The keep walker gives switch-entered claims the same `from_bp`. The check now starts at
    `max(from_bp + 1, route start)` (below).

#### Round 2 of the verification (2026-10-15, against 81c620ef)
The verifier (the lane directory's `verify/VERIFY-ORACLE.md`, "Round 2") confirmed all four disagreements as
oracle errors, and confirmed that the nine fixes close H1-H5 and M1-M5. It ruled the exact form of both fixes.
- **(a) Merged arms** (`traverse_compare._merge_invariants`), as SPEC 6.5 and 7.1 state them:
  - *Displays.* A displayed walk that is not a keep walk must be a real graph walk: every node a k-mer of
    the index, and no edge twice. It must equal a keep walk from its own latest merge on, at the same
    node there. A walk with no merge before its end must be a keep walk.
  - *Merge parents.* Every merge's parents end at one node.
  - *`labels_via_parent`.* It partitions the merged labels, each under a parent whose `labels_at_end`
    holds it.
  - *The majority parent.* `parents[0]` brings the most labels, unless a memory budget applies. Ties keep
    the arrival order, which the answer does not show, so they are not judged.
  - *Own routes.* Each leaf label's own route is a keep walk the label holds. A lineage switched twice
    before a merge is marked unjudged; none occurred.

  Tests: the three real answers pass. The display through the minority parent fails, and a label said to
  enter through a parent that never carried it fails.
- **(b) Displayed support** starts at `max(from_bp + 1, route start)`. Test: `T-canonical_k10-0076`.
- **Merged arms before their first merge** are held to the keep walker exactly: walks, label ends and their
  events, switches, blocked and hairpin events. After it, the united history stays invariant-only. Three
  merge-structure corruptions joined the mutation set: `merge_minority_first`, `merge_via_wrong_parent` and
  `merge_early_reason`.
- **`loss_budget` at random.** Transition-chain requests use random seeds, values, budgets, `switch_on` and
  branching, with two switches ahead.
- **Low items:**
  - the provenance test builds its own git tree, and with it found that `build_src_uncommitted` was always 0
    (a pathspec bug, now fixed; the eng-int tree is clean);
  - the "graph walk" correction above;
  - the bounded CANONICAL row is now stated.

Validation: six runs on the code of `6ec3d280` against the same 81c620ef build. Seed 1 ran at the default
sizes and seeds 2-6 at twice the size (`runs/v4final_*`), over every route and all twelve configurations.
- 16,406 requests: 7,128 `/pattern`, 1,484 `/resolve`, 7,794 `/traverse`.
- 44,291 answers: 40,799 equal the oracle and 3,492 are inconclusive with every invariant passing.
  **0 disagree.**
- The four answers of the first final runs were re-requested. Apart from timing they are identical, and they
  now pass.
- 440,277 corruptions of the stated mutation set, 125,027 of them of `/traverse` answers, every one caught.
  55 of them are merged-arm structure corruptions.
- Per seed:
  - 23-37 `loss_budget` label ends with their `needed_budget` (3 before);
  - 3-5 merged arms, 24 in all, every one passing the exact rules;
  - 103-229 merge arms that merged nothing, held strictly;
  - 14-15 trace queries telling trace from presence;
  - 13-23 PRIMARY twins.
- **No candidate engine defect.**

#### Round 3 of the verification (2026-10-15, against 81c620ef)
The verifier ruled that the oracle can merge and that "no candidate engine defect" stands
(`verify/VERIFY-ORACLE.md`, "Round 3"). Closed since:
- **The entry a merge keeps** (SPEC 6.5: the min (loss, branches) entry of a label on both parents).
  - A run the merge closed must not be cheaper than its label's continuing run. Loss is constant along a
    run, so the loss part is exact.
  - Every run, not only a leaf's end labels, must have an own route that is a keep walk its label holds.
  - This catches the verifier's G13 (a costlier loss-1 entry kept) and F3 on `T-canonical_k10-0076` (a label
    claimed through a parent that carries it only after a switch, ending before a leaf).
  - The mutation set gains `merge_costlier_entry`.
- **Provenance.**
  - `build_head` is the tree's HEAD now. The binary prints no commit (`--version` is the release number), so
    a binary linked before the last commit touching `metagraph/src` is flagged `stale` ("STALE BINARY" in
    the report). The main tree's binary that phases 1 and 2 used first is flagged this way; the eng-int
    build is not.
  - Every git call runs with `--no-optional-locks`, so the engine tree's index is never rewritten.
  - *The limit of judging by time* (verification round 4, R4.3). Only "linked before the last commit touching
    `metagraph/src`" is caught. Two kinds of build are not:
    - a build made from an older checkout after that commit, such as a branch switched back, or a tree that
      gets the commit later by a pull or rebase whose committer time precedes the link;
    - a binary linked before uncommitted `src` edits made afterwards.

    `build_src_uncommitted` reports the second case only while the edits stay uncommitted. Only the sha256
    pins the binary itself. The engine prints no commit that would let the report name the commit it was
    built from.
- **End reasons after a merge.**
  - `dead_end` only at a graph tip.
  - `label_lost` only where no successor keeps the label. A successor holding it counts only if the label was
    re-entered there by a switch.
  - Both hold for every keep-oracle claim of the v4final requests.
  - They catch the verifier's last E3b miss (`T-primary_k11-0023`, `label_lost` at 32 after the merge at 21,
    reported as `dead_end`). The mutation set gains `merge_late_dead_end`.

Validation: six runs on the code of `bafc074a` against the same build. The report records it as not stale,
with 0 uncommitted `src` files. Seed 1 ran at the default sizes and seeds 2-6 at twice the size
(`runs/v5final_*`), over every route and all twelve configurations.
- 16,406 requests and 44,291 answers: 40,799 equal the oracle and 3,492 are inconclusive with every invariant
  passing. **0 disagree.**
- Every one of the 16,509 answer files (twins included) equals round 2's apart from timing.
- 440,280 corruptions of the stated mutation set, every one caught. 125,030 are of `/traverse` answers, and 58
  of those are merged-arm ones.
- **No candidate engine defect.**

Still open, stated:
- L3 (the `lower` of a bounds count);
- L4 (records holding an N);
- the united history's other label states after a merge;
- ties of the majority parent;
- a server-side traverse memory budget, which the suite does not use, and which would make `first_parent`
  misjudge;
- `loss_budget` comes from one chain fixture.

#### Round 4 of the verification (2026-10-15, the last; against 81c620ef)
The verdict was "merge it" (`verify/VERIFY-ORACLE.md`, "Round 4"). Its three follow-ups are closed:
- **Re-entry after a merge.**
  - A `label_lost` end whose successor still holds the label is legal only if a switch re-entered the label
    there, and the re-entering run is strictly cheaper than the lost one (6.3: min(stay, switch); a tie goes
    to the stay).
  - `switch_costs`: every switch of a merged arm is at the cost model's price, and there are none under
    forbid.
  - The verifier's W2 (0076, re-entry at loss 1 against the stay's 0) and W1 (0051, a switch under forbid)
    fail.
- **Run reasons.**
  - `run_reasons`: a run's `end_reason` is its `label_end` event's reason as served (13, I3-2).
  - A `dead_end` run ends at a graph tip or under the hairpin qualifier.
  - The seed-2 E3b/G8 miss on 0076 (a superseded run reported `dead_end`) fails. The mutation set gains
    `merge_late_run_reason`.
- **The limit of time-based staleness** is stated under round 3 above.

Validation: six runs on the code of `dc9aecf9` against the same build (not stale). Seed 1 ran at the default
sizes and seeds 2-6 at twice the size (`runs/v6final_*`).
- 16,406 requests and 44,291 answers: 40,799 equal the oracle and 3,492 are inconclusive with every invariant
  passing. **0 disagree.**
- All 16,509 answer files equal round 3's apart from timing.
- 440,283 of 440,283 corruptions caught, 61 of them merged-arm ones.
- **No candidate engine defect.**

## Suite integration S8i, stage 3: the oracle merged (2026-10-15)
Lane O's oracle (gr/suite-oracle at ef0dd05d, verified final in VERIFY-ORACLE.md "Round 5") cherry-picked onto the
frozen v5 (a76e13b1) in order, 39 commits with their messages and authorship: oracle/, tests/test_oracle_*.py,
tests/data/oracle_*_answers.json and the section above, which stays its own top-level section after the frozen v5
(the only conflict, at the first README commit; nothing of either side dropped). The oracle imports only the standard
library and its own package, and no suite module, the bench or run_suite.sh imports it (test_oracle_core checks the
oracle's side); its tests import the package, so oracle.compare / oracle.report never shadow the pipeline's modules.
Validation: 702 tests pass (593 + the oracle's 109); the oracle's 109 pass inside the git worktree and from an export
outside any git tree (the provenance tests build a git tree of their own); the frozen v5 manifests and both v5 gate
files are byte-equal to a76e13b1's; the v3 / v4 re-evaluations equal the previous head's line for line, the v4 gate
the same 46 problems. No server was contacted; the oracle's evidence runs against 81c620ef are in its lane directory.

## Manifest v6: the /pattern budgets at 110,000 ms, and gating POST /pattern (2026-10-15 / -16)

Round 1 (2026-10-15, this section up to "Round 2") set the budgets at 120,000 ms; the owner's round 2 (2026-10-16,
"Round 2: gating POST /pattern" below) set them at 110,000 ms and readied the gate for POST /pattern. Where the two
differ (the budget, the client timeouts, the ingress, the manifests' and policies' contents and shas), round 2's
figures are the current ones.

The owner's decisions (2026-10-15), after the accepted baseline run of the frozen v5 (20261015T161813Z_refseq97-81c620ef,
the deploy gate's baseline under gate/policy_v5_rt_mcp.json): every inconclusive execution of that run was a /pattern
scenario at its 60 s budget (P-heavy-0001, M-protein-0001 / -0003, the M-planted ones, ...: a 503 deadline at ~60.2 s
or a time-limited 200). (1) The /pattern budgets go from 60,000 to 120,000 ms, /pattern only (the host's /pattern cap
is 600,000 ms; /resolve and /traverse stay at 60,000, their host maximum), as manifest v6 = v5 with exactly that
change; the gate's scenarios stay byte-identical, so the v5 acceptance carries over and deploy 1 is not delayed.
(2) The report marks every printed execution time above 10 s yellow and above 55 s red.

### The manifests
- manifest_lib.PATTERN_TIME_BUDGET_MS is 120,000 (the default of pattern_request, named on both arms). Of the 455
  POST /pattern scenarios (457 with the two GETs of /pattern/capabilities), 449 named the default and get it in
  their request and, where the MCP runs them (429), in their MCP arguments; the deliberate budgets stay (P-caps-0004's 300 ms, M-protein-0010's 30,000 ms; P-edge-0012 and the
  C-protocol scenarios send no JSON object). Nothing else of any scenario and no relation changed: build_manifest.py's
  SMOKE_TAG stays `smoke-v5` (SMOKE_SELECTION_VERSION: the tag names the smoke selection, which v6 keeps), and the
  MCP exposure baseline it reads is truth/mcp_exposure/v4-20261014T061053Z-budget-free.json, the v4 baseline with
  its request digests rebased from manifest/suite_v4.json to the budget-free rule (`build_manifest.py
  --exposure-rebase truth/mcp_exposure/v4-20261014T061053Z.json --exposure-source-manifest manifest/suite_v4.json
  --exposure-out ...`; mcp_request_sha256 leaves pattern_search's time_budget_ms out as it does traverse_resolve's:
  a budget decides how far the answer gets, complete or limited by time, never which fields it carries, and a
  larger budget lets no fewer checks decide, so a list recorded under 60,000 ms is a lower bound under 120,000 ms;
  EXPOSURE_RULE's shorter wording, "a budget shapes no field of the answer", is embedded in the frozen manifests and
  stays; the exposes lists are v4's byte for byte, and every scenario of v6 has v5's list, 512 with one). Regenerating the baseline from the v4 run with --exposure-from is no longer byte-identical to
  the committed file (the run's lines were re-evaluated after it: Q-strands-0001's list), which is why the digests
  are rebased rather than the baseline regenerated.
- `manifest/suite_v6.json` and `manifest/suite_smoke_v6.json`, built on this tree with `SUITE_TRUTH_CACHE=<the record
  cache v4 was built from> build_manifest.py --no-fetch` (deterministic: two builds, one content): manifest_version 6,
  700 scenarios, 742 executions per arm, the same counts as v5, content_sha256
  e416fdc53241fd6e8da5114a0bba97c9a12253d949390495b60d232c14f8ec4e, relations
  39660cfa0fb98cc47ac923fad3aa33419a84a3d92d662154cf2d4164752f4a86 (v5's); smoke 14 scenarios, content_sha256
  183295a285c6da48aa44e3a47b580d574ab8678376a5242ff36eff7be085587e. The field-by-field diff v5 -> v6 over the 700
  scenarios: `request.time_budget_ms` 60000 -> 120000 on 449 scenarios and `mcp.arguments.time_budget_ms` 60000 ->
  120000 on 429, all POST /pattern, and nothing else; 251 scenarios byte-identical, the 22 relations identical.
- `gate/policy_v6_rt_mcp.json`, generated with `gate_policy.py manifest/suite_v6.json --name v6-rt-mcp --arms
  mcp,direct-local --declared gate/declared_v6_rt.json --out gate/policy_v6_rt_mcp.json` (gate/declared_v6_rt.json
  names v6's content and declares nothing, as v5's): it differs from policy_v5_rt_mcp.json in the pins alone (the
  manifest's file, content and version, the name, the declared file and its sha256, generated_utc) and equals it with
  them masked: 241 selected / 36 excluded, required on MCP 205 scenarios / 416 executions, on direct-local 237 / 273,
  the 5 relations, no property exemption, the by-design table d4b8636c..., no /pattern scenario. The 237 scenarios it
  selects are byte-identical to v5's (their checks included), and so are the 5 relations: the v5 acceptance carries
  over (test_builder_v6, test_gate).
- run_suite.sh keeps v5 as its default (`--smoke`: suite_smoke_v5.json) until the owner approves the v6 freeze; the
  switch is the one line of its two defaults (`MANIFEST="$HERE/manifest/suite_smoke_v5.json"` / `suite_v5.json` ->
  `_v6`), or `--manifest manifest/suite_v6.json` per run.

### The client timeouts a /pattern execution meets
A 120 s budget answers its 503 deadline at ~120.25 s at the latest, so every client timeout must exceed it with
margin, or an inconclusive turns into an error:
- the direct runner's timeout per timing class was fast 120 s (slow 600 s, long 1800 s), exactly the budget:
  suite_common.client_timeout_s gives a request the larger of its class's timeout and its stated time budget +
  CLIENT_BUDGET_MARGIN_S (30 s), and run_direct.timeout_for uses it (the 397 fast /pattern scenarios get 150 s on
  direct-local, the 52 slow ones keep 600 s; every 60 s /resolve and /traverse budget keeps its class); direct-remote
  stays capped at REMOTE_TIMEOUT_CAP_S (125 s), where the front proxy answers first (below). The 30 s margin is
  narrower than the window in which checks._budget_503 still reads a late 503 as inconclusive (up to max(1 s, the
  budget) past it, 120 s for a 120 s budget; v5's fast timeout of 120 s covered that whole window for a 60 s
  budget): a 503 delivered between 30 s and 120 s past a 120 s budget is a client timeout (an error) on
  direct-local, not an inconclusive. The lateness seen so far is at most 0.3 s, so the margin stays at 30 s
  (a longer timeout only delays a run with a hung request);
- the orchestrator's wait for the mex runner after Ctrl-C is SUITE_WAIT_RUNNER_S (150 s) instead of a fixed 120 s;
  the waits after a cancelled MCP job (suite_common.request_bound_s) and the host-cache window (end to start, 60 s)
  follow the budget as before;
- the MCP runner (the async session's suite_mcp_runner.py, not in this tree) needs no change: its deadlines are
  300 s (fast) / 900 s (slow), each get_status long-poll 60 s with a 30 s HTTP margin, the submit 120 s (the
  pattern_search submit returns at once); the service accepts a budget up to min(1,170,000 ms, the host's cap) and
  its host-call timeout is META_CALL_TIMEOUT_SECS (1,200 s). The run's snapshot of suite_common.py is the one it
  executes, so the new rule reaches it; a re-vendoring is due only for the shared rules, which did not change.
- **The front proxy: direct-remote handled, the MCP arm open (the owner decides).** The public URL's ingress (nginx
  at https://metagraph.ethz.ch:8080) cuts a request at ~120 s (the long class is not run on direct-remote for it; the
  by-design table's note on pattern.time_budget_above_host_cap: it cuts "every arm through it"). A /pattern scenario
  that runs to its 120 s deadline answers its 503 at ~120.25 s, after the cut, and a time-limited 200 at ~119.8 to
  120.0 s races it: on direct-remote the proxy's 504 arrives instead of the host's answer. The checks read that 504
  as what it is, the ingress's cut before the host's answer could arrive, which shows nothing of the host:
  checks.ingress_cut_allowed makes the status inconclusive (the allowed 503 deadline's case) for a request whose
  stated budget reaches the cut (FRONT_INGRESS_CUT_MS 120,000 less FRONT_INGRESS_CUT_TOLERANCE_MS 10,000: the cut is
  measured, not stated) when the 504 arrives no earlier than that and carries no host body, and the other checks are
  not testable on it; a 504 on a request whose budget ends before the cut (every /resolve and /traverse scenario,
  v5's /pattern ones) stays the status fail it was, since the host had to answer within the budget (test_checks,
  test_builder_v6: nothing the deploy gate selects can meet the rule, so a gate run is judged as under v5). So on
  direct-remote a v6 /pattern execution that runs to its deadline is inconclusive, as under v5, never a fail; the
  deploy gate (mcp and direct-local, no /pattern scenario) is untouched; direct-local, on mex against the
  container's port, meets no proxy and gets the host's answer. **Still open: the MCP arm.** The service runs on
  Railway and reaches refseq33m-experimental through the same ingress as far as its code shows: its settings model
  that ingress with a cut of 110 s (TRAVERSAL_INGRESS_CUT_S_BY_HOST, "measured only for the dev ingress, nginx in
  front of refseq33m-experimental": the measured ~120 s less a margin) and fit traversal batches to it
  (app/traversal/planner.py ingress_batch_cap), while the pattern path sends the request's budget as it is. A 120 s
  pattern search through it ends in the ingress's 504, which the service retries after its 60 s backoff
  (PatternBackendError.retryable for any 5xx; META_RETRY_MAX 2): about 120 + 60 + 120 s, past the MCP runner's fast
  deadline of 300 s, so the execution ends in `error` (and a second 504 would end the job as a failed search, a
  status fail), where v5 had an inconclusive; three MCP errors in a row stop the paired loop (run_suite.sh
  MCP_ERRORS), and v5's deadline-bound MCP executions ran back to back (M-planted-0007 to -0009, M-protein-0001 to
  -0007), so a full v6 run can stop in area M. Nothing in this tree changes that: the runner and the service are the
  async session's (read-only here) and the ingress is infrastructure. The owner's options: a /pattern budget of at
  most ~110,000 ms (the service's own figure for that ingress; PATTERN_TIME_BUDGET_MS is the one line, then the
  manifests, the policy and their pins are rebuilt as above; the gate identity is unaffected), an ingress read
  timeout of at least ~150 s with REMOTE_TIMEOUT_CAP_S raised with it, or v6's /pattern scenarios on direct-local
  only. Not verified live (no host was contacted).

### The report's colour marks
report.py marks every printed execution time (client wall, server elapsed, MCP phases; the tables and the lines of
the findings, failures, timing, warm-up re-runs, overheads and the change against an earlier run) above
TIME_MARK_YELLOW_MS (10 s) yellow and above TIME_MARK_RED_MS (55 s) red (fmt_time, by the value as printed: fmt_ms
rounds 100 ms and above to whole milliseconds, so 10,000.4 ms prints 10000, not above 10 s, unmarked); differences
between arms and connection set-ups stay plain. The mark is a span with a class (`t-yellow` / `t-red`) and the tint inline
(`background-color:#ffe866` / `#ff8a80`): a Markdown viewer that renders HTML shows the colour, a terminal the plain
number, and on the md -> pdf path pandoc reads the span natively and report_colour.lua (`--lua-filter`, passed by
to_pdf) turns the class into a `\colorbox` of the tint with `\fboxsep` 1pt (the default 3pt opened a visible gap
before a comma after the box), the colours defined in header-includes (xcolor, loaded by pandoc's template). The
text stays black and legible on both tints. One legend sentence follows the front matter. table() sizes and aligns
columns by a cell's visible text. Times inside a check's free-text detail (a finding's or failure's detail line)
were not marked in round 1; round 2 marks them, and the gate's problems, by the printed value
(report.mark_quoted_times; below).

### Validation (offline; no host was contacted)
- 720 tests pass, none skipped (with SUITE_TRUTH_CACHE, the record cache v4 was built from): 702 at 07dfa4875 + 18
  (test_builder_v6 8: the versions and digests, the diff is the budgets alone, the gate's scenarios and relations
  are v5's, the smoke selection, the rebased exposure baseline, the request digest, every request's client timeout
  against its budget, the frozen v6 is what the builder makes; test_report_colour 7; the client-timeout rule in
  test_common and test_runner; the v6 policy in test_gate). The v6 build is deterministic (two builds, one content);
  the field-by-field diff v5 -> v6 and the identity proof are in the lane's notes (V6-DIFF.md, v6_diff.py). The
  frozen v5 manifests and both v5 gate files are untouched (their shas as above), run_suite.sh's default is v5.

### Round 2: gating POST /pattern (2026-10-16)
The owner's decisions of 2026-10-16, after the readiness study for gating /pattern (the v5 baseline judged offline
under a dry /pattern policy: GATE FAIL, 103 problems, none an engine or service defect: 44 time at the 60 s budget,
58 fields the MCP's /pattern answer does not carry, REL-subset-0007; and four checks gaps that could pass a wrong
answer):
1. **The /pattern budget is 110,000 ms, not 120,000** (manifest_lib.PATTERN_TIME_BUDGET_MS): the public ingress at
   https://metagraph.ethz.ch:8080 cuts a request at ~120 s, and the staging service reaches the host through it (its
   own traversal cut setting is 110 s), so a 120 s search would end in the ingress's 504 on the MCP arm, which the
   service retries past the runner's deadline (VERIFY-V6-2 H1). 449 /pattern requests and 427 MCP argument sets name
   it (P-scope-0004 / -0006, marked at scenario level, send none); the fast class's client timeout for them is 140 s. checks.ingress_cut_allowed keeps its constants: the cut
   120,000 ms less the tolerance 10,000 ms covers the 110,000 ms budget inclusively, and the late-503 window (up to
   max(1 s, the budget) past it) reaches past the cut measured a tolerance late (test_builder_v6). On direct-remote
   only (the arm behind the ingress), its 504 is a time limitation where the checks and the comparison look for one
   (checks.time_limited with the scenario, client and arm; run_direct and reevaluate record it so; compare's timing
   rule reads it; VERIFY-V6-2 L1).
2. **The v6 gate also covers POST /pattern** (arms mcp and direct-local) if the v6 baseline shows it decisive:
   gate/policy_v6_rtp_mcp.json, generated with the rt policy's inputs and the routes plus POST /pattern
   (`gate_policy.py manifest/suite_v6.json --name v6-rtp-mcp --arms mcp,direct-local --declared
   gate/declared_v6_rt.json --routes "POST /resolve,POST /traverse,GET /traverse/capabilities,POST /pattern" --out
   gate/policy_v6_rtp_mcp.json`): 721 selected (480 /pattern), 94 excluded (budget_stress 35, manifest_not_implemented
   50, table_not_expressible 47, advanced_path 1, marks_cover_repeat 1), required on mcp 628 scenarios / 1,275
   executions, on direct-local 686 / 793, 19 relations. Its resolve / traverse part is gate/policy_v6_rt_mcp.json's,
   which still equals policy_v5_rt_mcp.json but for its pins (the deploy gate of today, unchanged).
3. **Heavy /pattern scenarios.** Warm-ups (areas.PATTERN_WARM_UPS, P-warm-0001..0025, by the suite's mechanism:
   `warms`, place_warm_ups, as R-warm for the walks) right before the 25 /pattern scenarios whose MCP cold call
   (repeat 0, which the paired loop runs first) ended in the host's 503 deadline relayed by the service or
   time-limited in the v5 baseline (M-planted-0007/8/9/13..18, M-protein-0001..0007, N-nearmiss-0006/13/20/27/29/30/32,
   N-usecase-0008, P-heavy-0001): each is its target's request, timing_class long (direct-local only: it warms the
   host for the MCP's first call), with the generic checks, its target's v5 figures in its truth, and tagged
   budget_stress (its cold run is by design). In 23 of them a later execution (the warm MCP repeat or a direct arm)
   completed; M-protein-0002 and P-heavy-0001 completed none, warm runs included, and are budget_stress themselves
   (P-warm-0024 / -0025, added by round 2's verification, L4: the owner's rule names every timed-out cold call; the
   warm runs of P-heavy-0001 hit the budget too, so its warm-up may change nothing). A selection that takes a target
   takes its warm-up (an rtp gate run runs 19 of them, 22 direct-local executions; the budget_stress targets' are not
   selected). After a resume between a warm-up and its target, the warm-up runs again on direct-local right before
   the target's first MCP step left (the resume rule's before_mcp, above). The warm repeat is by manifest index, so
   the 25 warm-ups would shift the deploy gate's warm repeats: build_manifest.keep_gate_warm_repeats moves free
   /pattern scenarios of area N, and P-twin-0001 (right after the R-twin block, which no free scenario of area R
   precedes: R-twin-0005 would gain a warm repeat), so that every scenario on the gate's routes keeps its v5
   warm-repeat status, and gate_warm_problems refuses a build that changes one (the 237 gate scenarios, their order
   and their required executions are v5's). Outside the gates, against v5, the shift moves the direct warm repeat of
   61 /pattern scenarios and of C-caps-0003 (GET /capabilities), and P-scope-0004 / -0006 lose their MCP warm repeat
   by their scenario-level mark. budget_stress (outside the gate, still run and reported; areas.PATTERN_BUDGET_STRESS with
   each one's reason from the baseline): M-protein-0001, -0002, -0003, -0007, P-heavy-0001 and N-nearmiss-0030.
4. **What the MCP cannot show is labelled by its true class at the narrowest scope.** Verified against the
   service's pattern_search answer, the runner's normaliser (scripts/bench/suite_mcp_normalise.py at 9d76035e) and
   the baseline's recorded MCP answers: the entry's absence_scope (the service states its own scope), strands and
   pattern text are not carried, so the entry_fields checks that read them are marked not implemented yet via MCP
   (mcp.not_implemented_checks, 22 scenarios), scenario-level where that is the only property check (P-scope-0004,
   -0006: mcp.not_implemented). The host's top-level predicate block and timing.load_ms ARE carried, on every entry
   (by_database[].predicate / .load_ms): the checks read them there on the MCP arm (checks.MCP_PATTERN_ENTRY_FIELDS;
   Q-pred-0011 / -0012, P-edge-0009), so no MCP runner change is needed. D-P-scope-0004's repeat 0 tests nothing via
   MCP (its other property check, deterministic, has no result on a first repeat): gate_policy's rule
   marks_cover_repeat excludes that MCP repeat alone. entry_fields compares every field the answer carries before it
   reports a missing one (round 2's verification, L3): P-scope-0001 / -0002's entries carry `scope` but not the
   host's absence_scope, so `scope` is held on the MCP arm too, and a wrong carried field fails (a mark sets aside only
   a check that reads not_exposed). The scenario-level marks keep the MCP arm from running P-scope-0004 / -0006, so
   REL-subset-0002 / -0003 (suffix scope's contexts among any_offset's; via MCP in the v5 baseline they passed, 60 of
   640 and 4 of 55) are required on direct-local only: keeping them on the MCP arm needs the MCP to run a fully marked
   scenario whose executions the policy excludes, which the S8i rule refuses today (the owner's decision; L2 of the
   verification). A deliberate service behaviour is service-side by design: REL-subset-0007 is set aside on mcp by a
   declared entry naming the table's case pattern.max_contexts_above_10000 (the service refuses P-modes-0004's
   max_contexts 10,001 at submit where the host clamps); it stays required on direct-local.
5. **The report** marks "not implemented yet via MCP" red like times above 55 s (red means needs attention) wherever
   it counts or lists it with something in it, and the verdict gets a red attention line counting the TESTS per kind,
   only the nonzero ones: failures, errors and findings first, with the disagreements between the arms and the failed
   relations (round 2's verification, L3 of the checks: they fail a run whose executions all pass, as the exact-count
   rule can), then the time-budget issues (an execution hit its budget: a 503 deadline, a time-limited answer, a stop
   by time, the ingress's cut; or took more than 55 s), then what is not implemented via MCP (whole executions and
   property marks). Times quoted in a check's text and in the gate's problems (whose table also gains each
   execution's wall) are marked too; every mark follows the printed value. The time-budget kind counts tests that
   stop at their budget by design as well: the warm-ups' cold runs and the heavy budget_stress scenarios (the line
   states how many of its tests are tagged budget_stress, still counting them) and the budget tests (P-caps-0004's
   300 ms, R-budget-0001 / -0002, T-budget-0001), so it is not 0 on a healthy run; counting them apart or leaving them
   out is the owner's decision (L5 of the manifest's verification, L4 of the checks').

The checks fixes the readiness study asked for, each with a test that fails on the code before it: on the MCP arm a
/pattern time stop is held to the work time by the entries' own timing (their elapsed_ms and the reads beside the
search, never the parts of a search), an answer without any timing reads inconclusive, never pass (MCP-T); a stop
codon never substitutes a residue, and X and * hold none (SD-012, SPEC §21.3); after a sticky budget stop every later
pattern states the same stop, unknown counts and zero work (SD-013's class, SPEC §7.6); an exact /pattern count that
differs between two access paths is a disagreement (compare.exact_count_differences; GPT's "doubled count"), stated
with its paths and values in both Agreement tables, and blamed on the odd arm (compare.exact_count_side: the MCP only
when its count is one no other arm states while two others agree, a direct arm when the MCP states another direct
arm's value, undetermined with two arms alone; round 2's verification, M1 and L2); the ingress's 504 as above; and a direct pair that both state a time limitation (P-heavy-0001) differs by timing, not as
a disagreement.

The manifests: manifest/suite_v6.json 725 scenarios (v5's 700 + the 25 warm-ups), 767 executions per arm, content
49907dbb95ca0a373199d40cabe888b4c136f376295a1fc526a742eae8231cd1, relations 39660cfa... (v5's); suite_smoke_v6.json
content ee731e491d7e1168665acb439de6a705b66c0bf3c3cbacae7007d2611277f7dd. The frozen v5 files are untouched.

#### Validation (offline; no host was contacted; the run directory was only read)
- The v5 baseline re-evaluated with this tree's checks (a copy, round2/r5/reeval) against the copy judged with
  07dfa4875's: 0 of the 237 rt gate scenarios' executions change verdict or any check, and the rt gate with
  policy_v5_rt_mcp.json passes as before (689 judged, 635 passed). Of /pattern, no execution's verdict changes; 6 MCP
  checks go from not_exposed to pass (answer_fields of P-edge-0009, Q-pred-0011, Q-pred-0012, r0 and r1), and
  P-heavy-0001's direct pair reads timing instead of disagree. No exact count differs across arms (0 false alarms).
- The dry /pattern gate on the v5 baseline (round2/r5/rtp: v5's requests with v6's marks and tags, the dry policy
  generated as the rtp one): 31 problems, 1 the emulation's identity and 30 time at the 60 s budget, all on warm-up
  targets (10 relayed MCP deadlines, their 10 unverified by-design pairs, 9 time-limited MCP cold calls, 1 direct-local
  time-limited cold run of M-protein-0006); no other problem. The readiness study's 58 not-exposed problems,
  REL-subset-0007 and D-P-scope-0004 are gone.
- 740 tests pass, none skipped, with SUITE_TRUTH_CACHE (722 after round 1, 734 at the start of round 2's last
  attempt with two failing; the new ones: the per-entry fields, the warm-ups and keep_gate_warm_repeats, the rtp
  policy, marks_cover_repeat, the report's attention marks and the gate table's wall). Each checks fix's test fails
  on afb46ec11's code. The v6 build is deterministic (test_frozen_v6_is_what_the_builder_makes).
- Round 2's verification fixes (2026-10-16; the Opus verifiers' round2/verify-checks and round2/verify-manifest;
  evidence in round2/fix): the v5 baseline re-evaluated with this tree's code changes no check outcome or detail of
  its 2,989 executions and no comparison entry against the verified b90cfa738 (the rt gate passes, 689 / 635); the
  dry /pattern gate gives the same 31 problems; a probe with doubled exact counts (P-exact-0003 MCP, -0005
  direct-remote, -0006 direct-local) adds exactly three MCP disagreements, side mcp for -0003 and side direct for the
  others, with their counts in both Agreement tables and the attention line; two v6 builds are byte-identical; the
  rt gate's 729 executions and policy_v6_rt_mcp.json (v5's but for its pins) do not move. 745 tests pass, none
  skipped; each new or changed test fails on b90cfa738's code.

## The results document (doc/, 2026-10-15)

One PDF per run for the MetaGraph developers: for every test what it checks, what was expected, what each of the
three arms answered and how long it took, compact enough that the 681 tests of v4 fit a readable document, every
test findable by id, every imperfect result explained in a detail block. Design: `~/.cache/metagraph-traverse-plans/
2026-10-15/doc/DESIGN-DOC.md` (the sample of chapter R was reviewed first). The generator lives in `doc/`, the PDF is
never committed.

    doc/make_doc.py RUN_DIR --out DIR [--chapters R,T,...] [--sample] [--no-pdf] [--recorded]
                    [--units FILE] [--prose-dir DIR] [--checks-dir DIR] [--by-design PATH]
                    [--policy FILE] [--also-policy FILE ...] [--declared FILE] [--identity FILE]
                    [--xelatex PATH] [--name doc] [--scaffold-prose]

After the owner's review of the chapter R sample (round 2, 2026-10-15) the document has this shape: the cover carries the
run verdict and the deploy gate as two boxes of equal weight (a gated run: the gate verdict and its counts; an ungated
run: "Deploy gate: not evaluated" and why; a gate judged offline by gate.py for the document says so; a partial build
states "Chapter R sample; whole-run summary shown for context"); Part I has four sections (I.1 decision summary: the
scope, the verdict lines, the principal findings, failures and disagreements of the whole run, the inconclusive
executions by reason and the properties still untested, computed per unit as the property checks no execution of a
member passed on any arm, with the reasons per arm, the prose of `overview.json` printed in slots between the computed
parts; I.2 methods; I.3 reading guide with one real matrix annotated by numbered markers, the legend, the marks, the
cell grammar incl. what the slash means per table type, and the record key of the ground truth's abbreviations; I.4
provenance and coverage as a compact table, the recorded outcomes per arm and the limitation and exemption counts
linked to the appendices); each chapter leads with its conclusion (`summary`), the summary block, what it leaves
untested (`coverage_gap` plus the computed lines) and only then the background (`intro`); every unit table carries a
legend line saying what `/` means there (repeats; median / maximum; the edge tables' MCP column: the recorded status,
the service's code and the class of the outcome, derived from the evidence: refused before host execution, not
expressible through MCP, answered after service normalization, host refusal relayed through the service, with an
explanation line under a normalised row); a replicates cell reads `not exercised: n` instead of `ok 0/0` and a
coverage line in words replaces "exceptions: none"; a unit whose MCP column is wholly not exercised carries one line
under its table with the short reason and the pointer to appendix B (`mcp_notes.units` of the chapter prose, else the
cases' short reasons); a detail block holds the expected line, an adjudication line only for what the rows do not say
(a failure, an error, a host-state override, a disagreement, a planned execution not recorded, and one phrase for an MCP
arm the service decided: its class and status, the short reason unless the unit's note under the table has it, the
case in B), then one row per arm and repeat (final outcome with its qualifier: allowed 503, late 503, host state, by the
checks; HTTP and code, wall, requested budget, property tested?; an arm standing ○ or ◆ reads so in its rows), a
finding's evidence per execution with the mechanism and the rule once in a line set apart, the overridden checks of a
host-state execution labelled "intermediate", every error text once, the classes, the relations and the supporting
artifacts (the raw files and the executions file's line numbers); blocks under about 28 estimated lines are unbreakable
boxes, longer ones repeat the test id on continuation, headings ask for space below them; the closing section opens with
observations filled from `closing.json` templates (sample sizes included), the caveats, a paired-latency plot per arm
pair and the per-phase box plots (pgfplots, from `plot_*.csv` the generator writes; nothing stacked), then the tables,
each section led by its method, observation and limits paragraphs (`closing.json` `sections`, the figures as
placeholders); appendix D holds the short glossary with the version-pinned README reference and the full record key,
F the environment in full, the rules of judging, the verdict and the gate (`overview.json` `appendix_F.rules`), the
vocabulary and the gate's problems (every scenario-arm pair with its repeats), G the MCP limits and the declared
exemptions (F and G with their lead-ins from `overview.json`).

After the verification of round 2 (VERIFY-R2, 2026-10-15; round 3) the facts the document states about the gate and
the by-design table come from the records that decide them, never from a count that looks alike: what blocks the gate
is the gate's own outcome problems (`sections.gate_outcome_counts`: the executions on a gated arm, in the policy's
selection and not declared exempt; I.1 and I.4 say how many of the inconclusive executions those are, or that no gate
was judged), and the by-design table the document prints is the one the re-evaluation replayed (`reader.Run._load_table`
finds it by the sha256 of `reevaluated.json by_design` among the replayed path, the suite's copy and the run's snapshot;
another table is read only with the comparator's warning kept and printed in I.4 and F with both shas). The closing
section's figures follow one population each and say which (the pairs within 1 s count |difference|; the scenarios a
pairing leaves out are counted for that pairing; the phase box plot reads the same per-scenario population as C.3's
table and states the zero paging times it leaves out; C.3's cold-wall population against C.4's whole cold column); a
family's prose may count the MCP classes of its members through placeholders the generator fills from the same
classifier as the cells (`Context.family_values`: `{n_mcp_refused_before_host}` ...), so that prose and tables cannot
drift. Part I is four pages: the methods hold the arms, the order, the ground truth and the timing definitions (the
owner's list) and the rules of judging, the verdict and the gate moved to appendix F.

After the re-verification of round 3 (round 4, 2026-10-15) every figure of the run that a chapter's prose states is a
placeholder the generator fills from the record, and every sentence whose truth depends on the run is a conditional,
so that a later run's document never prints the sample's numbers (`doc/facts.py`): the chapter prose reads
`facts.chapter_facts` (the counts by standing and by agreement, the summary's timing figures, the findings with their
overruns, the claims each family states in its `claim` sorted into what the direct arms showed, the footnoted tests
with one phrase per arm, the gate's outcome problems among them; `fam[<family>][...]` a family's values, `x[<scenario>]
[dl|dr|mcp][wall|phrase|detail|class|...]` a scenario's), a family's prose `facts.family_facts` (its standing counts, the
MCP classes with their ids and variants, the members' phrases, the findings, the 503 deadlines past their budgets, the
stated stops, the relations, the owner's declared exemptions, the expected runs); `[[flag ? when true | when false]]`,
`[[!flag ? ...]]` and `[[x[R-edge-0003][mcp][class] = refused_before_host ? ...]]` choose a text, nested innermost first;
`{n:s}` is the plural s, `{n:v}` the verb's s, `{n:word}` the number in words. A placeholder or flag the run's values
lack stays visible and `make_doc.py` warns. The literal numbers a prose file still carries are the constants of the SPEC,
of the manifest (lengths, budgets, expected runs) and of the service, each listed under the file's `_assumes` with its
source; `tests/test_doc.py` fails on any other number, on a placeholder or flag the generator does not compute, and
when chapter R's prose rendered on the fake run prints one of the v4 sample's figures. The compact record key
(abbreviation and organism, one paragraph) sits on each chapter page before its first family, where the matrices
need it (in I.3 it cost Part I a fifth page), every unit table whose variants name a record links to it and to
appendix D's full key from its legend, and appendix D prints the key before the glossary. The phase box plot's tick
labels are two lines (the phase, its n). The closing section says which population each n is (per scenario the
host-warm MCP execution with a backend timing; the cold calls without a backend timing counted apart).

After the owner's second review (round 5, 2026-10-15; OWNER-REVIEW-SAMPLE-2.md, the design frozen) the paired
scatterplots of C.2 draw a plotting area of 62 by 55 mm each (`latex.paired_plot`, `scale only axis`); an MCP
agreement's depth is judged over the repeats compare actually compared (`reader.depth_over_compared_repeats`: the
executions inconclusive by a host state, `mcp_agreement.host_state_repeats`, are left out before the depth is read),
so a detail block says "fully checked on r1; r0 was an inconclusive host-state outcome" (`blocks.depth_text`), the
chapter line counts the fully checked pairs compared over fewer repeats than recorded and the chapter prose names
them (`n_mcp_full_fewer_repeats`), and the reading guide and the glossary define the depth over the compared
repeats; no timing sentence states a cause or a bound: the owner's sentence "Latency varies substantially across
records; these correctness controls do not isolate the effect of query length" replaces the record-versus-length
reading, and "observed under concurrent host load" every "upper bound"; the host-cold MCP calls without a backend
timing are split by what ended them from the record (`sections.cold_mcp_untimed_split`, the recorded time
limitation or the error code); a filled number is joined to its unit by a non-breaking space (`latex.join_units`);
a variant label prints its thousands separators (`units.variant_text`, "2,000 bp (pfa)"); the sentences of
chapter R that asserted the run (every call refused at submit, the twins' answers the same via MCP, the source
outranked in a tie, host-warm on direct-local) are conditionals on computed facts (`mcp_all_refused_before_host`,
`rel[...][mcp]`, `x[...][dl][host_word]`; a unit's note under its table reads `facts.unit_facts`); the R-budget
prose names its footnoted members once and lets the coverage gap and the detail blocks hold the evidence. The test
`test_two_fixtures_change_conclusions` renders chapter R on two fixtures built from small data (the fake run, and a
copy with a stated stop, a late 503 and a clean direct side) and asserts that the counts, the verdict branch, the
coverage statements and the explanatory prose all change.

After the re-verification of round 5 (VERIFY-R4.md, the final patch, 2026-10-15) the closing section states no
timing cause: each of its limits names what the run measures and what it does not isolate (the service's queue,
the host's load, the cache, the loop's order, the answers' sizes; `doc/prose/closing.json`), and the
interpretive clauses of C.1 are conditionals on flags computed from the run (`sections.closing_values`:
`backend_share_small`, the median backend overhead below half of the median MCP overhead;
`paired_seconds_to_milliseconds`, the paired medians as printed; `cold_warm_mcp_close`, the MCP cold-warm ratio
within a quarter of one), so that a clause prints only when the figures make it true (`sections.observations`
fills through `facts.fill`); the number of a count's verb and noun follows the count in R.json ("is a finding" /
"are findings", "the only test" / "one of the two tests"); and `latex.join_units` also covers the findings list
of I.1, the summary lines under the tables and appendix B's quoted texts. `test_closing_two_fixtures_change_conclusions`
renders the closing section on the perturbed fixture and on a copy whose summary and walls flip each condition,
and `test_timing_wordings` bans the removed wordings ("point to", "reflects the", "rather than a slower host",
"at this load" and the others) in every prose file and in the generator.

- Reads the run directory as report.py does (manifest.json, the re-evaluated executions and compare.json when the
  run has them, else the recorded ones, run.json, reevaluated.json, samples.jsonl, gate.json, arms/mcp/run.json),
  the by-design table compare judged with, the owner-declared file of `gate/` whose manifest sha matches, and the
  checks module (docstrings for the check catalogue). Stdlib only, offline, no clock in the output.
- Writes `DIR/doc.tex`, one `chapter_<X>.tex` per chapter and `doc.sty`, then runs xelatex twice (a third time when
  the .aux changed; `/Library/TeX/texbin/xelatex` or PATH; latexmk is not used). Fonts from TeX Live: Libertinus
  Serif, Fira Sans / Mono, DejaVu Sans for the glyphs. Exit 0 built, 1 LaTeX failed (the log's first error echoed),
  2 usage or missing inputs (each named), 4 a comparison to compute (no compare.json) with a by-design table
  compare.py refuses (as compare.py and report.py; `reader.TableRefused`).
- The reading model (`doc/units.py`): a scenario belongs to one unit (one builder template, one property, varied
  over one dimension), a unit to one family (the id prefix), a family to one chapter (the area). Five layouts:
  `series` (arms as rows, variants as columns), `replicates` (aggregates per group), `refusals` (glyph only),
  `list` (the refusals shape with time cells), `pairs` (twins, relations, D repeats). The metadata comes from the
  manifest's `doc` blocks (v5) or else from the mapping `doc/units/suite_v<version>.json` (v4); a scenario in
  neither is an error, never a silent drop.
- Prose (`doc/prose/<chapter>.json`, `overview.json`, `glossary.json`, `checks.json`, `closing.json`, `cases.json`):
  keyed on family and unit ids, never on scenario ids, so a new scenario of an existing unit is covered for free.
  `<chapter>.json` holds the chapter's `summary` (the conclusion, printed first), `coverage_gap`, the `intro`
  paragraphs (the background), one title and paragraph per family and one property line per unit (a string, or
  {property, note}) and, under `mcp_notes.units`, the one line under a unit table whose MCP column is wholly not
  exercised; `overview.json` Part I's four sections (keyed decision / methods / reading / provenance, each with
  `paragraphs` and the slots the generator prints between its computed parts: I.1 `verdict`, `findings`,
  `untested`, `closing`; I.4 `coverage`, `exemptions`; the wording of the glyph legend and the marks and the
  annotated matrix's marker texts under `reading`; the lead-ins `appendix_F` and `appendix_G`) and the placeholders
  it uses (`{database}`, `{k}`, `{n_scenarios}`, `{gate_text}`, `{inconclusive_by_reason}` ...:
  `sections.Context.placeholders` fills them from the run, as plain text); `closing.json` the closing section's
  intro, the observation templates (filled by `sections.closing_values`, whose vocabulary the file declares), the
  caveats, the plot captions and, per table section (`sections.C3` to `C7`), the method, observation and limits
  paragraphs with the run's figures as placeholders; `cases.json` the short local reason per by-design case the
  tables, the detail blocks and appendix B print (one clause: the generator adds the class in words, the status and
  code and the pointer to appendix B); `checks.json` the check catalogue of appendix A (per check kind: title,
  asserts, inconclusive, fails, args, the SPEC references, per-route rules, and the `common` rules evaluate()
  applies before any check and the host-state adjudication after them; it is preferred over the docstring in
  checks.py); `glossary.json` the short entries of appendix D ({term, definition, see}, printed in alphabetical
  order, with the `lead` that names this README's glossary section and the commit that pins it, and the record
  key's lead). Inline markup: `code` and *emphasis*. `--scaffold-prose` writes the placeholder entries a chapter
  file lacks (existing text untouched). Chapter R, Part I, the glossary, the closing section, the case reasons
  and R's check kinds are written (2026-10-15, round 2); the other chapters hold placeholders marked "(to be
  written)". A chapter's prose names the run's figures as placeholders and its run-dependent sentences as
  conditionals (round 4; `doc/facts.py` computes the values, each file's `_assumes` lists the constants it may
  carry as text). `tests/test_doc.py` fails when a family, a unit or a check kind has no entry at all, when a written
  chapter (`WRITTEN_CHAPTERS`) still carries a placeholder (its summary, coverage gap, intro, families, units and
  unit notes, the overview's sections and appendix lead-ins, the glossary, the closing section's every key, the
  case reasons), when a catalogue entry's generic flag disagrees with `checks.is_property`, when the prose uses
  a placeholder the generator does not fill (or declares one it does not), when a prose string carries a number
  that is neither a placeholder nor an `_assumes` constant, or when chapter R's prose on the fake run prints a
  v4 figure. A key whose prose is awaited can be named in the test's `ROUND2_PENDING` and tolerates a placeholder
  until it lands (empty since round 2 landed).
- The cell grammar (`doc/cells.py`): glyph (pass, fail, error, finding, inconclusive, not implemented via MCP,
  service-side by design) + marks (footnote digit, disagreement, timing / same checks / not comparable, cold /
  aged, mixed repeats, unverified) + wall times in repeat order with the unit once, the server time in brackets,
  the MCP backend time on a third line. A test that is not clean on some arm gets a footnote mark and a detail
  block right after its unit's table (`doc/blocks.py`: expected vs obtained, each check's outcome and detail, the
  classes, the relations, the SPEC rule of a finding). The verdict is report.py's own (`report.verdict_text`,
  extracted so both read the same rule).
- Tests (`tests/test_doc.py`, offline): the v4 mapping complete and valid, the prose complete and the written
  chapters free of placeholders, every check kind described, the generator on the fake run of
  test_compare_report (the .tex, the index, the verdict, the escaping; the PDF when xelatex is installed), the
  cell grammar, the expectation formatters on chapter R's real expectations; and, for the review's changes, the
  four MCP outcome classes on fixtures and on the fake run, the adjudication and the property-tested reading, the
  detail block's shape (adjudication first, one row per arm and repeat, the artifacts' line numbers), the coverage
  wording, the pagination helpers (short and long boxes, the headings' space), the plots' data files (deterministic
  across two builds), the closing templates' placeholders, the untested-properties computation, the timing legends,
  the record key, the glossary pin, the gate line on the cover and the sample note; and, for round 3, the gate's
  outcome counts and clause, the replayed by-design table found by its sha (and the warning when another is read),
  the closing figures' populations, the evidence clause's parenthesis balance, the finding line's single mechanism,
  the rows' qualifiers and the one-phrase MCP arm, the family values from the classifier, the repeats in appendix F
  and the rules in appendix F; and, for round 4, the prose filling (conditionals, specs, the index syntax, the
  spacing, the unresolved recorded), no v4 literal on the fake run, every placeholder of R.json resolvable on v4,
  every literal number assumed, the facts on the fake run, the record key on the chapter page and in the legends,
  appendix D's order, the box plot's labels, the untimed cold count and the uncut B note; and, for round 5, the
  plots' size, the depth over the compared repeats (the reader, the block line, the chapter line, the prose), the
  timing wordings (no "upper bound", no record-versus-length cause), the untimed cold split, the non-breaking
  number-unit space, the variant labels, the unit note's facts, and the two fixtures whose conclusions differ;
  and, for the final patch, the closing section's conditionals on two fixtures, the banned timing wordings, the
  count grammar and the joined units outside the prose.
- The sample of 2026-10-15 (Part I, chapter R, the closing section, the appendices chapter R needs, the index)
  was built from the v4 run re-evaluated with c11db9dd: `doc/make_doc.py RUN_DIR --out DIR --sample`; its notes
  for the owner's review are `~/.cache/metagraph-traverse-plans/2026-10-15/doc/sample/SAMPLE-NOTES.md` (round 1)
  and `.../doc/sample-r2/SAMPLE-NOTES.md` (round 2, after the owner's review, point by point; round 3, after the
  verification of round 2, rebuilt in the same directory), `.../doc/sample-r3/SAMPLE-NOTES.md` (round 4, after
  the re-verification of round 3), `.../doc/sample-r4/SAMPLE-NOTES.md` (round 5, after the owner's second
  review) and `.../doc/sample-final/SAMPLE-NOTES.md` (the final patch, after the re-verification of round 5,
  `sample-r4/VERIFY-R4.md`). The template is final (the 14-unit grouping of R stays; pagination may evolve);
  rounds 4 and 5 and the final patch changed no table cell, only the prose's sources and wording, the record key's places, the plots' size and
  labels, the depth line and the variant labels' separators.

### The rebase onto the suite at 07dfa487 (F-gate, the frozen v5, the oracle; 2026-10-15)

`gr/suite-doc` (8 commits on c11db9dd) rebased onto 07dfa487. Two conflicts: `report.py`, where the doc lane had
extracted `verdict_text` and F-gate had made a gated report lead with the gate's verdict and added the UNVERIFIED
branches (a run, or an arm, with executions and no pass), resolved so that `verdict_text` holds the whole rule
(`no_pass_arms` among its facts) and `verdict_section` prints its headline through `verdict_headline` (the
exploratory label of a gated report, `exploratory_label`); and this README (both sections kept, this one last).
The document states the report's verdict and gate line by the same code as `report.py`: the cover's run verdict is
`verdict_text`'s, labelled exploratory (`report.exploratory_label("beside")`) when the document shows a gate verdict;
the cover's gate line is the run's `gate.json` when recorded and current (its words `report.gate_headline_parts`,
the ones the report's gate section prints), else, when `gate.stale_inputs` names files the run changed after the
verdict (the report's "GATE verdict STALE"), the stale verdict is named and the gate is judged offline by `gate.py`
with the policy pinning the manifest, as a run without `gate.json` is; else "not evaluated". The generator reads
the newer run fields without breaking v4: a warm-up's re-run for a resume (`merge_run`: never a line of the merged
set, `warm_up_reruns` on its target's line and `run.json` `reruns_for_resume`) is a row of I.4 and a table of
appendix F, a line in the target's detail block and its raw file in the artifacts, never an outcome; the
comparison's `executions_sha256` is held against the executions file the document reads (I.4 and appendix F: "the
lines this document reads", or the sha256 of other lines, or "states no executions sha256" for a comparison before
2026-10-15); the check items the by-design table explains (compare's `mcp_agreement.explained_items`, the gate's
reading, never the line's own claim) are grouped as explained in the detail block; the v5 manifest's property-level
marks (`suite_common.mcp_check_marks`: `mcp.not_implemented_checks`, `mcp.by_design_checks`) read as "not
implemented via MCP (the manifest's property-level mark: reason)" or "service-side by design (… case)" for that
check, in the detail block's checks line, the "property tested?" cell (the marked items counted apart) and the
untested-property lines' per-arm reasons, never as untested or as a failed check; an execution inconclusive by a
scenario-level mark (`suite_common.SCENARIO_LEVEL_MARK`: the marks cover every property check) reads as such in
the rows, the adjudication and the MCP class (`marked`), never as a host state; a check read `not_exposed` carries
the inconclusive glyph and the words, never the failure glyph; a by-design case's `explains_checks` entry that is an
object ({check, detail_matches}, the d4b8636c table) prints in appendix B; the glossary has the property-level mark.
The mapping file is picked by the run's manifest content sha among `doc/units/*.json` (`units.default_mapping_path`),
else by the version name. The chapter R sample rebuilt from the v4 run (`sample-rebased/`, no option; `sample-rebased-
table/` with `--by-design` naming the table the re-evaluation replayed, d1e040aa, extracted from git) is byte-identical
to `sample-final/` in `chapter_R.tex`, `doc.sty` and the three `plot_*.csv`; `doc.tex` differs in the cover (the
exploratory label; the offline gate 45 -> 50 problems: 5 identity, of which one is the F-gate table problem the
rebase expected (the comparison read d1e040aa, the regenerated v4 policy pins d4b8636c), one the MCP arm's
unobserved index_fp (the v4 runner recorded no `service_capabilities.host_document`), two the re-evaluation's
checks directory (the main tree, now at 07dfa487) no longer holding the checks.py and suite_common.py the
re-evaluation records, one the comparison stating no executions sha256), the I.4 and appendix F rows above, the
glossary entry and the README pin, and, without `--by-design`, the by-design table (the replayed d1e040aa table is
no longer in the tree: the live d4b8636c one prints with the WARNING, appendix B two cases longer).

### The v5 unit mapping (`doc/units/suite_v5.json`, 2026-10-15)

The frozen v5 manifest (content 72181c4a…, 700 scenarios) has no per-scenario `doc` blocks and cannot get them, so
its mapping is a file keyed to the manifest's content sha, as v4's is, built by `units_v5_build.py` (the lane's
notes directory) on the v4 builder's rules: a v4 id whose title did not change keeps its committed v4 label (645
of 681), a changed title gets its label from the new one (R-budget's stops at 480 / 420 ms after their complete
twins at 60 s, the continuation from base 2,500 / 5,000, R-discover's limits 9 and 3, R-rc's eco, R-edge-0023,
T-walk's groups by the new extensions, N-negative-0018's 12 residues, the T-labels / T-strand label cuts). The 107
v4 units are kept with their members and variants; two units are v5's: `R-warm` (family R-warm, replicates: the
eight warm-ups of T-repeat, T-walk, T-labels and N-hard targets) under chapter R, whose 14 units of the owner's
grouping stay as they are, and `T-longwalk` (family T-longwalk, pairs: the four gate walks T-walk-0006 / 0012 /
0018 / 0036 each followed by its 10 kb budget_stress twin T-longwalk-0001..0004), the only move of v4 ids between
units. R-edge-0025 joins `R-edge-accepted`, T-attempts-0003 `T-attempts`, T-caps-0002..0004 `T-caps`,
T-quorum-0002 `T-quorum`, T-seeds-0002 `T-seeds`; the 22 relations are placed as in v4. 109 units in 78 families:
37 replicates (339 scenarios), 23 series (97), 22 list (78), 20 pairs (125), 7 refusals (61); every series at most
8 variants. The maintenance tests run over both mappings (complete, valid, keyed to the manifest's sha, the moves and
the v5-only placements named) and over both for the prose rule; the prose files carry the mapping's drafts for the
v5-only family and units (`R.json` R-warm, named in `ROUND2_PENDING` until the v5 prose round; `T.json`
T-longwalk among T's pending prose). `tests/fake_doc_v5.py` builds a FAKE v5-shaped run from the frozen manifest
(every scenario on the three arms with synthetic passing answers, the MCP marks read not_exposed, one warm-up run
again for a resume; nothing of it is evidence about a host), and `V5ScaffoldTest` builds the full document from it
without a PDF: every v5 id once in the index, the ten chapters, the two v5 units, the marked checks set aside, the
re-run as timing evidence. The prose round for v5 must cover, besides the two new units: R.json's `summary`,
`coverage_gap`, the R-budget family prose and unit line and `mcp_notes.units.R-budget`, which assume v4's figures
(10 kb under 3,000 / 8,000 ms, the full twins at 20 s, the continuation from base 5,000).

The full document (2026-10-15): the nine other chapters' prose (T, P, L, X, M, Q, C, N, D) follows chapter R's rules
exactly (a chapter page that leads with its conclusion and its coverage gap, family prose plus one property line per
unit, every run figure a placeholder or a conditional on computed facts, no timing cause the run does not isolate,
"ok" only for what was tested, a catalogue entry for every check kind and a short reason for every by-design case
named), written for the union of the v4 and the v5 mapping with `[[manifest_version = 5 ? ... | ...]]` for a sentence
true in one manifest alone. The facts their prose reads beyond chapter R's (`doc/facts.py`): `manifest_version` at
every level; per arm the first execution's `time_limited` text and its checks by name
(`x[<sid>][dl][checks][walk_in_record][detail]`); the time limitations that left a 200 inconclusive, in the three
readings the chapters use (chapter T's time cuts of a walk, `n_time_cut_tests`, `time_cut_ids`, `time_cut_text`,
and the tests no arm showed the property on; chapter N's and Q's `n_time_limited` / `time_limited_text`; chapter M's
/pattern lists the deadline withheld, `n_time_limited_200_tests` / `time_limited_200_ids`); `direct_all_pass`; and at
the chapter level `unit[<unit id>][...]`, a unit's own counts (chapter M's substituted residues apart from the cap
slot of the same family). A 200 a time limitation touched but that passed is not counted: the claim it carries was
shown. `tests/test_doc.py` covers every chapter: no placeholder, no literal run figure, no banned timing wording;
every placeholder and flag resolvable on the v4 and the v5 mapping (a scenario named outside a family one mapping
alone holds is one both manifests hold); every chapter's prose rendered on the fake run prints none of the v4
figures; chapters M and T, like R, rendered on two fixtures whose data differ change their conclusions, their coverage
gaps and their family prose (`FullDocumentTest`); every check kind of both manifests in the catalogue, every by-design
case named with its short reason. A check kind without an `expect_text` formatter prints `name(args)` with break
opportunities after ( , = _ (`blocks._generic`), so the expectation column of a list table wraps instead of running
into the next column. The v4 document (every chapter, the replayed d1e040aa by-design table, two builds
byte-identical; 164 pages at that build, the current build's page count in its NOTES.md) is at
`~/.cache/metagraph-traverse-plans/2026-10-15/doc/full-v4/` with its NOTES.md.

### The tints of the printed execution times (2026-10-15, corrected after VERIFY-COLOURS-V5 on 2026-10-16)

The owner's decision of 2026-10-15: every printed execution time above 10 s is set on a yellow tint and above 55 s
on a red one (`cells.TIME_TINT_YELLOW_MS` / `TIME_TINT_RED_MS`, strictly above; the document's own constants, read
by nothing else; doc.sty's `timeyellow` FFE866 / `timered` FF8A80 behind black text, `\tyellow{}` / `\tred{}` with
`\fboxsep` 0pt, so a tinted figure sets exactly as wide as a plain one and the pagination does not move). A wall, a
server time, a backend time or an MCP phase, alone or as a median, p90 or maximum, in the cells, the replicates
aggregates, the summary and timing lines, the detail blocks' rows, phases and re-run lines, the closing tables
(`cells.fmt_exec`, `fmt_times`; the closing tables' `_stats_cells`). A difference (`fmt_delta`: a delta, an overhead,
how late a 503 was), a budget, a cap, a connection set-up and a distance in time carry none. The figures quoted in
the prose and inside a check's text carry their tint too (VERIFY-COLOURS-V5 item 1, 2026-10-16): a fact the prose
files fill in prints its time through `cells.fmt_exec_text`, which carries the tint as a marker (`latex.TINT_OPEN` /
`TINT_CLOSE`, private-use code points no prose or run text uses) that `latex.esc` turns into the macro after
escaping (`set_tints`; a marker a cut left open is closed), so the walls, server and backend times of
`x[<sid>][<arm>]`, the p50 / p90 of the chapter and closing facts and the walls of the time-limited lists are tinted
in the sentences that quote them; an execution time quoted inside a check's or the runner's text, budget_respected's
`client wall N ms` and the MCP runner's host-state `... after N ms`, is marked by its value (`cells.tint_quoted_times`,
`esc_times`; `facts._cut` / `_clause` on the prose side, the finding, checks, host-state and adjudication lines of a
detail block and I.1's lines on the LaTeX side), while a lateness (`delivered N ms after`) and a budget in the same
text stay plain. The legend sentence of the cell grammar (I.3, `sections.times_legend_sentence`) states the rule
and its exemptions.

The other corrections after VERIFY-COLOURS-V5 (2026-10-16): I.1 names a disagreement whose checks agree by its
outcomes ("the outcomes differ, repeat 0: no time limitation tells the outcomes apart", compare's `outcome_differs`)
instead of "a check differs"; a walk a stated resource other than time stopped (a /traverse answer's
`resource_stop`: the work or memory budget) is read from the run's raw record (`reader.Run.raw_answer`,
`blocks.resource_stops`), the detail block's adjudication names the resource, its knob and its figures (`seed 0: the
work budget (bounds.max_work_units, 100,000,000 units; used 100,000,275 units) at 3,618 bp (left), ...`), and chapter T's
prose reads `n_cap_stop_tests` / `cap_stop_ids` / `cap_stop_text` and, with the time cuts, `n_budget_stop_tests` /
`budget_stop_text` (`facts._budget_stops`), so "no walk was stopped by its budget" holds for every budget and the
"depends on the load" clause names the work cap where it stopped a walk; the reading guide's column note names the
annotated unit's own dimension and its record-key clause only where the variants carry one (`annotated_dimension`,
`annotated_record_key`), and a note whose marker sits in no cell of that matrix says "(not in this matrix)"; the
clean rule of I.3 says what the cells do (an unverified by-design pair carries the mark u, its cause in the chapter's
prose and appendix B, no detail block); the coverage gaps of R, M, N and P print "Outside <the family> every test
passed" only when that family has a footnoted test (`n_footnoted_in_gap`), else "Every test passed"; R's coverage
gap has the comma of p. 14 ("Of the four relations, ..."); M's residue claim counts the direct-remote passes beside
direct-local's; appendix F's MCP runner row states that the run did not record the service's deployed revision
(the readability review's READ-021 wording, which names no manifest version); the cover's policy sha is labelled as
the canonical JSON's, with the file's sha beside it (`reader.Run.policy_file_sha256`).

### The readability review merged, red means needs attention, the tints by the printed value (2026-10-16)

gr/suite-doc-read (the owner-requested readability review of the full v4 report, 75 comments; its record is
`~/.cache/metagraph-traverse-plans/2026-10-15/doc/full-v4-readability-review/`) is merged into gr/suite-doc; every fix
of both branches is kept (the integration notes: `.../doc/integration/MERGE.md`).

The owner's decisions of 2026-10-16. **Red means "needs attention".** (1) A test not implemented via MCP (the MCP-only
class: the tools cannot express the request) is red wherever it is printed: the ○ glyph (doc.sty's `\gni`, colour
`attnred` C62828), its counts in the standing, outcomes and agreement rows, the coverage and grouped cells, the
detail blocks' rows, adjudication, classes and observed lines, the property-level marks (`blocks.MARK_GROUP`), the
untested table's reasons, I.4's tables, appendix B's class tag, and the prose's "could not be exercised through the
tools" counts (`cells.attn` for LaTeX, `cells.attn_marked` for a phrase the prose quotes, `!!...!!` as the prose
files' inline markup; a conditional compares a value without the markers). The unit notes, the ○ cells' class words
and the rest followed after the two verifications of e3098b2e9 (below). (2) The table of
contents carries each chapter's attention tag in red under the chapter's line (`sections.toc_attention`, doc.sty's
`\tocattention`, linked to the chapter), e.g. "2 time-budget issues, 9 not implemented via MCP", and I.1 an attention
table per chapter with the totals (`sections.attention_table`). Both count TESTS per kind, the nonzero kinds only, in
`facts.ATTENTION_KINDS`' order: failures, errors, findings, disagreements, time-budget issues (a test with an
execution that hit its time budget, `facts.time_budget_issue`: a 503 deadline, a recorded `time_limited`, the
service's own deadline, or a wall printed above 55 s), not implemented via MCP (the ○ standing, or an MCP execution a
property-level not-implemented mark set aside). A test counts once per kind and may count under several. The counts
are facts (`facts.attention_facts`: `n_attention_<kind>`, `attention_text`, `n_attention_tests`, on every chapter's
values); the reading guide's legend defines each kind once (`sections.attention_legend_sentence`, from
`facts.ATTENTION_DEFINITION`). (3) Every printed execution time is tinted by its **printed** value (`cells.time_tint`
reads `fmt_t`'s figure back, `printed_ms`; a figure quoted from a text is its own printed value): a time printed
"10 s" carries no tint whatever its raw value, "56 s" is red, as gr/suite-v6's report.py marks (this branch's
report.py has no time marks). The time a walk's stop used
(`blocks.resource_stop_phrase`) is printed with its unit and tint ("used 15,593 ms" on yellow; the attempt deadline's
"requested 25,000 ms; used 22,085 ms"), the work and memory figures with theirs; the quoted texts' `client wall time
of N ms` and `elapsed N ms` are tinted like `client wall N ms`.

The readability re-verification (`.../readability/verify2/VERIFY.md`): a chapter's reading line names the tests whose
MCP standing differs from their outcome in the outcomes row, both kinds (a refusal filed under not implemented that
records service_by_design; a by-design repeat beside a host-answered one, counted mixed or inconclusive), and prints
"none in this chapter" exactly when there is none (`sections.standing_outcome_differences`); I.3 states both kinds
and defines the two rows once; the /traverse budget scope says that a label derivation reads the budget after every
k-mer and stops at it (SPEC D3); I.2's percentiles name I.3's population; the gate's pass line prints its access
paths once; a partially checked class judged over the repeats says so; a /traverse 200 finding is "a walk stopped
outside its stated rule (the seed's clock)" (`sections.finding_kinds`); the document's FINDINGS regex is pinned to
`report.verdict_text`'s own words by a test (report.py's verdict still says "without failing the request").

After the two Opus verifications of e3098b2e9 (`.../doc/integration/verify-merge/VERIFY.md`,
`.../verify-colours/VERIFY.md`), decision (1) is carried out where the text is printed and the legend says exactly
what is red. Red now as well: the unit notes' and the family prose's "Not exercised via MCP" (`!!...!!` in the
prose files), the class word of a cell, an observed line or an adjudication whose MCP access path stands ○
(`blocks.class_word_tex`: "not expressible through MCP", "refused before host execution"), appendix B's
"not implemented mcp: N (...)" evidence state, appendix G's counts, the run statements naming the class (R's
reading line "read not exercised", the edge families' "not expressible through MCP (N inputs ...)", the not-exercised
family counts `mcp_not_exercised_by_family`, the classes of `mcp_classes_text` whose members stand ○, P-modes' cell
counts, now filled from the facts). A property-level mark not implemented via MCP (manifest v5) that set a check
aside on an MCP execution of a test whose MCP access path does not stand ○ (`cells.property_not_implemented`, the one
rule of `ArmState.ni_marked` and `facts.not_implemented_issue`) is marked on its cell by a small red ○ raised after
the glyph (doc.sty's `\gnip`), counted after its glyph pair in the standing, family, grouped and coverage lines
("4 ● (3 with a property check not implemented via MCP)"), named in the attention tag and under I.1's table ("18 not
implemented via MCP (3 by a property-level mark: T-labels-0003, -0005 and -0006)"), and listed in appendix G
(`sections.property_marks_table`). I.3's legend now states the rule as applied: the glyph, its counts (a count of 0
black), the class words of what stands not implemented, the property-level mark; headings, legends, definitions,
quoted texts and the by-design table's own statements name the class in black; a footnote digit is dark grey,
not red (owner, 2026-10-16: red means only "needs attention"); it numbers the test's detail block. A walk's stop figures round half up, as every
time does. I.1's finding sentence names the seed's clock for a /traverse walk.

## Suite integration: the frozen v6 and the results document merged (2026-10-16)
gr/suite-doc (895397620: the results document generator, doc/) merged into the frozen manifest v6 (gr/suite-v6 at
a71d9177f) on gr/suite-int-v6 with a merge commit; both branches start at 07dfa4875. Two files conflicted:
- report.py (3 hunks), both sides kept: the doc branch's verdict_text (the verdict word and reason, the rules
  unchanged), EXPLORATORY_LABEL / exploratory_label / verdict_headline and gate_headline_parts, which the document
  calls, and gr/suite-v6's time_budget_issue / attention_counts / attention_line, the red marks (mark_nim_counts over
  the headline, nim_count) and the gate table's wall column. The report any run gets is the one gr/suite-v6's
  report.py prints: one exploratory label, the wording both branches had ("Exploratory verdict (not the deploy
  decision, which is the gate's above)"; the document's cover says "beside"); gate_headline_parts prints a missing
  gate.json field empty, as md() printed it. The doc branch's split of agreement_depth's not_compared (no MCP repeat
  compared the check) and partly_compared (another repeat did) stays for the document; the report's depth table lists
  both, sorted, as gr/suite-v6 printed its not_compared (skipped_on_any_repeat, a follow-up commit: the split had
  merged without a conflict and changed the table of the v5 baseline).
- README.md: the "Manifest v6" section, then the results document's section and the glossary, which stays last.

The v6 unit mapping (`doc/units/suite_v6.json`) predated the final v6: keyed to the first v6 (e416fdc5, 700
scenarios), it had no unit for the 25 /pattern warm-ups, so the generator refused a v6 run ("25 scenario(s) without
unit metadata"). It is rebuilt for the frozen v6 (content 49907dbb..., 725 scenarios) by `units_v6_build.py`
(`~/.cache/metagraph-traverse-plans/2026-10-16/int-v6/scripts/`), the successor of units_v6_rekey.py, which
refuses the final v6 by design: `python3 -I units_v6_build.py manifest/suite_v5.json manifest/suite_v6.json
doc/units/suite_v5.json doc/units/suite_v6.json` (deterministic). Every v5 scenario, unit and family keeps its v5
entry (a unit's members in v6's order, as the manifest moved them); the warm-ups P-warm-0001..0025 form one new unit, `P-warm` (family P-warm, chapter P, replicates, variant
"of <target>"), as v5's walk warm-ups form R-warm. The script checks its premises and refuses otherwise: the new ids
are warm-ups, each right before its target; the v5 ids keep their titles and, but for the /pattern scenarios
build_manifest.keep_gate_warm_repeats moves (30), their order; the relations are v5's; on v5's ids only request,
mcp, truth and tags differ; no unit label of a pattern chapter names a /pattern budget. test_mapping_complete_v6
states the same over the frozen file.

Not done here (the results document's v6 prose round): chapter P's prose has no family prose and no unit line for
P-warm, so make_doc still refuses a v6 run, naming them ("prose missing: P: families.P-warm.prose; P:
units.P-warm"); and the prose's 86 `[[manifest_version = 5 ? ... | ...]]` conditionals take their other (v4)
branch for a v6 run.

### Validation (offline; no host was contacted)
- Tests: the full suite (SUITE_TRUTH_CACHE as in round 2) at 75c190465, 879 tests OK, none skipped (1,006 s).
- The v5 baseline (runs/20261015T161813Z_refseq97-81c620ef, read only; copies in the lane's notes) re-rendered with
  a71d9177f's and the merged report.py, the clock fixed: byte-identical, ungated, gated by policy_v5_rt_mcp.json
  (GATE PASS) and gated by policy_v6_rtp_mcp.json (GATE FAIL, 116 problems: the gate's headline with problem kinds).
  The merge commit alone differed in the depth table (the follow-up above).
- Its results document built with the merged generator: 174 pages, as the doc branch's build5; no missing character, no undefined reference, the one known
  overfull box (1.1 pt); every chapter file byte-identical to build5's.
- manifest/suite_v6.json content 49907dbb95ca0a373199d40cabe888b4c136f376295a1fc526a742eae8231cd1 (computed),
  manifest/, gate/ and every suite module but report.py byte-identical to a71d9177f; the frozen v5 files to
  07dfa4875.

## Open items (the harness, 2026-10-17)

- The /pattern finalisation reserve (SD-021, the owner's decision): checks.py still reads the reserve a /pattern
  answer states (`limits.finalize_reserve_ms`, else 250) when it bounds the work time; the host's is configurable
  (`--pattern-finalize-ms`) and checks.evaluate sees no host document. The NEXT manifest states it in the check
  arguments (`reserve_ms` taken from the capabilities snapshot by the builder); checks.evaluate's API does not
  change. The frozen v5 and v6 keep the answer's reading. Done in manifest v7 (below, item 2).
- The readers of an execution's raw record (SD-093): a truncated .gz, a cut JSON, a record that is no object or
  whose answer is none reads as no body at the three sites that crashed on it (compare.direct_error_parts, which
  then keeps the line's own error codes; run_direct's previous_entry, time_limited unknown; the bench panel's
  raw_body), as the other readers already did. One shared reader in suite_common (keeping today's lookup order,
  run_dir/rel then run_dir/arms/<arm>/rel) replaces the ten copies after the open gr/suite-* lanes land. Done on
  gr/v7-harness (below, item 4).

## The v6 results document: a delta on v5 (2026-10-17)
The owner (2026-10-17): v6 is very close to v5, so its document is v5's with only what v6 changes (v6 = v5 + the
/pattern budgets of 110,000 ms, the 25 direct-local warm-ups P-warm-0001..0025, the MCP property marks and the
budget_stress tags). Every change keeps the v5 document's chapters byte-identical (the validation below).
- The version switches: the prose's `[[manifest_version = 5 ? v5 text | v4 text]]` conditionals read
  `[[manifest_version >= 5 ? ...]]` (facts.fill compares a number with `>=`; `=` stays an equality), so a v6 run
  takes every v5 branch, the text unchanged. Chapter N's version note names the run's own version
  (`{manifest_version}`) and, for v6, what v6 adds in that chapter (the budget, the eight warm-ups before its
  /pattern searches, N-nearmiss-0030's budget_stress tag, N-usecase-0005's property-level mark, the moves in the run
  order).
- Chapter P's family P-warm (the 25 direct-local warm-ups of /pattern targets) has its family prose and unit line,
  modelled on R-warm's: the targets with their warm-ups' ids, the generic checks under the target's 110,000 ms
  budget, the budget_stress tag, and the warm-ups' order, re-runs and footnoted members from the facts. make_doc
  no longer refuses a v6 run. compare classes a scenario without MCP executions `missing`; the document reads that
  as "a planned execution not recorded" only where the run's plan schedules the scenario on the MCP arm
  (`reader.Run.planned_on`, None-safe for a run without a plan), so the warm-ups, which the plan runs on
  direct-local alone, are not footnoted for their absent MCP column (it reads "not exercised: 25 –", the legend's
  "no execution: the plan had none").
- The /pattern budget in chapters P, M, N, Q, X and L: the 13 sentences that state it (every request names
  `time_budget_ms` ..., a 503 deadline at the ... ms budget) read `[[manifest_version >= 6 ? 110,000 | 60,000]]`;
  P's "(the host's default, which the service sends too)" reads, for v6, "(above the host's default of 60,000; the
  MCP arguments name it too)". Kept: P-exact's smoke run that "ran past its 60,000 ms budget" (history), the
  /resolve and /traverse budgets and caps of chapters R and T (60,000 ms, their host maximum) and the ingress
  figures. Each file's `_assumes` names 110,000 with its source; V6DeltaTest renders the six files for v5 and v6.
- The deployed service's identity, which the MCP runner does not record: `make_doc.py --identity FILE` (a JSON
  file: run, service_commit, railway_deployment_id, source, and host_release, mcp_runner_commit, manifest, which
  reader.Run checks against the run's release, the runner's git head and the manifest's content sha, refusing a
  file that names another run or contradicts the record). I.4 and appendix F's MCP runner row print it as "the
  service's deployed revision, taken from <the source's first clause> and not recorded by the runner: service
  commit ..., Railway deployment ..." (the commit shortened in I.4, in full in F); without the option the text is
  v5's ("not recorded by this run").
- P-edge-0004 (lower case): the P-edge unit's note, printed under its property line for manifest v6 only
  (`[[manifest_version = 6 ? ...]]`), says that the service's owner states the MCP's /pattern answer leaves out the
  entry's pattern text by design, that the frozen v6 manifest still marks its `entry_fields` not implemented via
  MCP, and that the next manifest revision corrects the class. sections.unit_section prints no note line for a note
  that fills to nothing, so the v5 document is unchanged.
- No local file path in the document (an internal audience beyond this machine; hosts, staging names and ids stay,
  never a secret): latex.public_paths rewrites a path inside the suite's tree from `metagraph/scripts/suite/`, a
  path under a home directory without the home and its hidden directories, a temporary one as its file name, and
  appendix F applies it to the by-design table's path, the re-evaluation's recorded path, the table warnings and
  the manifest's ground-truth sources (the index records' `~/.cache/...` report directory); I.4 to its table
  warnings, the cover to gate.py's error. The v5 document's doc.tex changes in those two appendix F rows only.
- The gates of a v6 document (the coordinator's instruction, beside the brief): the deploy gate is
  gate/policy_v6_rtp_mcp.json, passed as `--policy` (without it make_doc takes the first policy of gate/ pinning
  the manifest, policy_v6_rt_mcp.json); `--also-policy gate/policy_v6_rt_mcp.json` (repeatable) judges a further
  policy offline and prints it in I.4's table as "other gate policy", with its verdict and counts, beside the deploy
  gate, which it does not replace. The v6 document is built with

      doc/make_doc.py runs/20261017T082918Z_refseq97-81c620ef --out DIR --policy gate/policy_v6_rtp_mcp.json \
          --also-policy gate/policy_v6_rt_mcp.json --identity identity-v6-baseline.json
- A layout fix the v6 build needed (not a brief item): a detail block's "property tested?" cell with items a
  property-level mark sets aside ("yes (1 of 1); 1 marked not shown via MCP") wraps in a 30 mm top-aligned box;
  D-P-scope-0004's MCP rows, footnoted in v6, ran 44 pt past the margin. No v5 detail row has such a cell.
- A wording fix the v6 run needed (not a brief item): I.1's and I.4's inconclusive executions by reason count the
  one inconclusive by a scenario-level mark (D-P-scope-0004's MCP r0, marks_cover_repeat) as "by a scenario-level
  mark (every property check marked as not shown via MCP)", never as "the service's ... (a host state)".
- An accuracy fix (the coordinator, 2026-10-18; not a brief item): a manifest scenario-level mark
  (build_manifest.scenario_level_marks: property-level marks that covered every property check of a scenario with
  an MCP tool mapping, each part "not implemented yet via MCP: <reason> (<check>)" or "service-side by design: ...")
  is answer-side. blocks.mcp_class classes it `answer_lacks_field` ("not implemented yet via MCP (the answer lacks a
  checked field)") and blocks.short_reason words it "the request is expressible through the tools, but the MCP
  answer lacks absence_scope, the field entry_fields checks" (the field from the check's arguments), never "cannot
  be expressed through the tools", which stays for the request-side marks (support, bounds.max_query_bp, a raw body)
  and the table's not-expressible cases. v6's P-scope-0004 / -0006 read so; no v5 scenario has such a mark.
- An accuracy fix (the coordinator, 2026-10-18; not a brief item): a relation the comparator could not judge on the
  MCP arm (`not_evaluated`, its detail "a raw answer is missing") because one of its scenarios stands not
  implemented via MCP and has no MCP execution reads "not judged via MCP: P-scope-0004 is not run via MCP (not
  implemented yet via MCP)" in the detail blocks' relations line and, in the mcp column of the pairs tables and
  appendix C, "not judged: P-scope-0004 not run (not implemented yet via MCP)", with the no-execution dash
  (cells.relation_unrun / relation_unrun_text); any other not-evaluated relation keeps the comparator's words. Chapter P's P-scope sentence says, when both suffix scenarios stand not exercised (v6), that neither is
  run via MCP and the two relations are not judged there. The rule is general, so it also corrects v5's document
  where the same case printed the comparator's words: chapter C's relations line under C-caps-0003
  (REL-same_answer-0004 / -0005: C-caps-0001 and -0002 are not run via MCP) and appendix C's three C-caps rows;
  every other v5 chapter stays byte-identical.

### Validation (offline; no host was contacted; the run directories were only read)
- The v5 baseline (runs/20261015T161813Z_refseq97-81c620ef) re-rendered with this generator (4c4d0c3ad): nine
  chapter files (all but C), doc.sty, the three plot files and doc.aux / .toc / .out byte-identical to build5
  (gr/suite-doc 895397620); chapter C differs in the one relations line the relation fix corrects (C-caps-0001 /
  -0002 not run via MCP), doc.tex in appendix C's three C-caps rows (the same fix), the README pin and appendix F's
  two path rows; 174 pages, no missing character, no undefined reference, the known 1.1 pt overfull box, no local
  path in the PDF text.
- The v6 baseline (runs/20261017T082918Z_refseq97-81c620ef) built with the command above (4c4d0c3ad): 169 pages, no missing
  character, no undefined reference, no LaTeX warning, the same single 1.1 pt overfull box (T-edge-0011), no local
  path in the PDF text; the cover's GATE PASS (policy_v6_rtp_mcp: 2068 judged, 1950 passed, 0 problems) and I.4's
  policy_v6_rt_mcp row (689 judged, 635 passed, 0 problems) equal the acceptance note's offline verdicts. The builds,
  the chapter-by-chapter delta against v5 with each difference's source and the numbers checked against the run's
  data are in `~/.cache/metagraph-traverse-plans/2026-10-17/doc-v6/` (DELTA.md).
- Tests: the full suite (SUITE_TRUTH_CACHE as in round 2) at 4c4d0c3ad, 889 tests OK, none skipped (780 s).

### The delta verifier's round (2026-10-18, on gr/int-sd)
The Opus verifier of the merged document (`~/.cache/metagraph-traverse-plans/2026-10-18/doc-v6-verify/VERIFY.md`) found
every number right and three medium and nine low wording errors; each is fixed, the v5 document changing only where
the text was wrong for v5 too:
- M1: an arm the run's plan does not schedule (v6's direct-local warm-ups on direct-remote and mcp) reads "no
  execution planned" in the coverage line (cells.coverage_text, ArmState.not_planned), the agreement row ("25 with no
  MCP execution planned") and the classes line ("no MCP execution planned (the plan runs this test on direct-local
  alone)"), never "not recorded" or compare's bare "missing"; a run without a plan reads as before.
- M2: chapter P counts the answer-side scenario-level marks apart from the requests the tools cannot express (chapter
  facts n_mcp_answer_lacks, mcp_answer_lacks_ids / _fields, n_mcp_not_exercised_tools,
  mcp_not_exercised_tools_by_family); where the manifest has such a mark (sections.answer_side_mark_ids), I.3's ○ row,
  the attention definition, appendix F's vocabulary, the glossary entry and appendix G's lead add "or, by the
  manifest's scenario-level mark, the MCP answer lacks every field the test's property checks read".
- M3: I.2's /resolve row names 60,000 ms for manifest v5 and later (v5 too); L1: the client time-outs state the
  direct-remote cap of 125 s (v5 too) and v6's budget-plus-30-s rule.
- L2: N's version note keeps the warm repeats of the scenarios policy_v6_rt_mcp selects; L4: P-warm's two targets
  "whose searches all ended at the budget in the v5 baseline"; L3: M's and N's "the rest read inconclusive" only when
  a direct test did not pass (v5's chapter N too); L5: Timing 4 names the population of its unpaired count (v5 too);
  L6: counts of one in the singular (facts' :s / :v read a count carried as text); L7: a not-exercised arm whose
  reason is the manifest's own mark without a by-design case points at appendix G (blocks.limit_appendix; v5's
  T-edge-0005 too); L8: chapter C's prose gives the document relations' reason (C-caps-0001 / -0002 not run via MCP;
  v5 too); L9: a check the marks set aside reads "set aside by the marks", not "not compared".
- Validation: the v5 document (174 pages) differs from build5 only in the lines DELTA.md section 7 lists; the v6
  document (169 pages, the re-evaluated set) and the v5 one have no missing character, no undefined reference, no
  LaTeX warning, the known 1.1 pt overfull box alone and no local path. The full suite at 868134d6a: 937 tests OK,
  none skipped (744 s). The changed lines against the verifier's builds:
  `~/.cache/metagraph-traverse-plans/2026-10-17/doc-v6/CHANGED-LINES.md`.

## Manifest v7 (2026-10-18)
v7 = v6 plus exactly the owner's approved changes, each below; v6 (content 49907dbb...), v5 and every frozen policy
stay byte-identical. build_manifest.py builds v7 by default (MANIFEST_VERSION 7) and the frozen v6 again, byte for
byte, with `--manifest-version 6` (build(manifest_version=6): test_frozen_v6_is_what_the_builder_makes); the v7
steps run only for version 7.

### P-edge-0004's pattern text: service-side by design (item 1)
The service owner (2026-10-17): the MCP's /pattern answer leaves out the entry's pattern text BY DESIGN. v6 marked
P-edge-0004's entry_fields check (it reads the host's upper-cased echo of the lower-case pattern) not implemented yet
via MCP (mcp.not_implemented_checks). v7 classes it service-side by design by the suite's manifest-side mechanism for
that class at property level, `mcp.by_design_checks` [{check, case}] (lane V5; "Property-level marks" above):
`{"check": "entry_fields", "case": "pattern.lowercase"}`, the table's case the scenario already names (mcp.by_design)
for the service's handling of the pattern text (it checks the text in either case and upper-cases it before it is
stored and sent, so the host never sees the caller's text from the MCP arm). Why not T-edge-0008's route: T-edge-0008
is modelled through the table, its case traverse.duplicate_seeds stating an explains_checks entry
{check, detail_matches} that explains the not_exposed answer_fields reading (verified as a by-design pair); the same
here needs an explains_checks entry on pattern.lowercase, a change of service_by_design.json, the service's own file,
which the suite does not edit. The mark sets the not_exposed reading aside on the MCP arm exactly as the not-implemented
mark did (suite_common.checks_marked_aside reads both marks alike), the check stays required on the direct arms, and the
comparison, the report and the results document list it under "Service-side by design (property)" instead of "Not
implemented yet via MCP (property)" (compare's marked_checks mark by_design). The rule is areas.MCP_PATTERN_ENTRY_BY_DESIGN
(field -> the cases a scenario may name for it) with areas.pattern_by_design_case: an entry_fields check whose every
field the MCP's answer cannot carry is left out by design, each under a case the scenario names, is marked by design;
any other (another unexposed field beside it, no case named) keeps v6's not-implemented mark. In v7 that is
P-edge-0004 alone (the only scenario reading the pattern text). Nothing else of the scenario changes: its request,
expectations, exposes and repeats are v6's.

### The /pattern finalisation reserve in the check arguments (item 2)
The host's /pattern finalisation reserve is configurable (SPEC-pattern-search §10.2, `--pattern-finalize-ms`), and
checks.evaluate sees no host document (its API stays as it is). v5 and v6 bound a /pattern work time (SPEC §7.6:
time_budget_ms - reserve - E, checks._pattern_work_time, the bound of every stated time stop or time limitation the
stated-reason rules hold: invariants, the count checks, the entry and slot checks, deterministic) with the reserve the
judged answer states (`limits.finalize_reserve_ms`, else 250), so an answer could widen its own window. v7 states the
host's reserve in the check arguments: every POST /pattern scenario's budget_respected expectation carries
`reserve_ms`, the capabilities snapshot's `pattern_capabilities.finalize_reserve_ms` (250 in truth/staging_capabilities/,
release refseq97-6f945db4), added by areas.state_pattern_reserve after the warm-ups are placed (454 scenarios: the 429
that allow a 200 and the 25 warm-ups; the 26 refusal-only /pattern scenarios state no budget_respected); the builder
refuses a snapshot without the number. checks._scenario_pattern_reserve reads it from the scenario (its
budget_respected expectation), and _pattern_work_time takes it before the answer's statement; without the argument
(every v5 / v6 scenario) the reading is the old one, so a v5 or v6 run judges exactly as before. checks.validate_expectation
refuses a reserve_ms that is no number or below 0 (CHECK_ARGS budget_respected: rule, cap_ms, reserve_ms).
check_budget_respected itself reads no reserve (its /pattern rule is a 200's elapsed time against the budget; the late
503 window is _budget_503's). Not done (not asked): a stated limits.finalize_reserve_ms other than reserve_ms is not
flagged as a problem of its own (/resolve's invariants do that for its fixed reserve).

### The service identity, suite side (item 3)
The service lane's runner (async's suite_mcp_runner.py at 684bb327, another repository; not changed here) records the
deployed service's commit: every arms/mcp/mcp_runs.jsonl line, run.json's top level and each invocations[] entry
state `service_commit` (read from the service before the invocation: GET /health, the `deployment` block of GET
/traverse/capabilities, outside service_document_sha256), `identity_changed` (the commit read after the invocation
differs: true / false; null when either read failed) and `identity_error` (the failed read). The owner (2026-10-18):
the commit only; a Railway deployment id is neither stated nor read (a fixture that has one is ignored).
- **The merge** (merge_run.service_identity / identity_summary, over mcp_runs.jsonl, else run.json invocations[], else
  run.json alone) writes run.json `service_identity`: `recorded` (any of the three fields stated: false for a runner
  before them), `service_commits` {commit: invocations that read it; run.json's top-level commit counted 0 times when
  no invocation states it}, `service_commit` (the one commit, else null), `identity_changed` (true when any invocation
  or the top says so; false when every invocation says false; else null), `invocations`, `without_identity`,
  `without_commit`, `changed`, `unconfirmed` (identity_changed null), `identity_errors` and `source`; null without an
  MCP arm. A change during an invocation or more than one commit is a warning of the merge (the gate judges).
- **The gate** (gate.service_identity_problems, identity problems on arm mcp; gate.json `service_identity` {expected,
  recorded}) takes the deployed commit as `--expect-service-commit SHA` (7 to 40 hex digits; an abbreviation matches
  by prefix, gate.commit_matches), which run_suite.sh passes on from its own `--expect-service-commit` (the run's and
  the re-evaluation's gate; refused without a gate) and which run-gate passes from the LIVE line. Given a commit: an
  invocation that read another commit blocks, one that read none (service_commit absent or null: not confirmed) blocks,
  a run that records no identity at all (an older runner) blocks (the gate cannot confirm what it judged), and a run
  without an MCP arm blocks where the policy's arms name it. In any case identity_changed true blocks, and so does
  more than one commit over the invocations without a commit given (the run met more than one deployment; given one,
  each other commit blocks as differing). Notes: identity_changed null (one read failed; the commit read stands or
  blocks as above), without a commit given an older run ("records no service identity ... not checked"), invocations
  without a commit beside others, and the recorded commit not compared. A run resumed across a deploy blocks: the
  superseded invocations' lines stay in mcp_runs.jsonl with the earlier commit (conservative: the run met two builds).
  The commits are compared as recorded (the runner always records the full sha).
- **The results document** (sections.service_version_text) prints the recorded commit, labelled "as the MCP runner
  recorded it (read from the service before each MCP invocation)", with a change, an unread commit or a failed read
  stated; the make_doc --identity file only as a fallback, "taken from <source> and not recorded by the runner"; a
  Railway deployment id only from that file, labelled as taken from it. reader.Run refuses an identity file whose
  service_commit contradicts the recorded one.
- Older runs (v5, v6 baselines) record no identity: run.json gains `service_identity` with `recorded` false, the gate
  notes it, and the document keeps "not recorded by this run" (or the input file).

### One reader of an execution's raw record (item 4; SD-093 option 2)
The ten copies that resolved a line's `raw`, gunzipped it, swallowed the read errors and dug out the record, its
answer or its body (compare._raw_doc, load_answer and direct_error_parts; reevaluate.raw_doc and raw_answer;
merge_run.host_state_inconclusive and renormalisation; run_direct.Runner.previous_entry; bench/analyze.answer_body and
bench/make_panel.raw_body), and the results document's eleventh (doc/reader.Run.raw_answer), call one reader in
suite_common: raw_record_path(line, run_dir=None, arm_dir=None) resolves the path in today's lookup order (an absolute
path as it is; under a run directory run_dir/<raw>, then run_dir/arms/<the line's arm>/<raw>; under an arm directory
arm_dir/<raw>, then the run directory two levels up for a merged 'arms/...' path; the first that exists),
read_raw_record (a dict, else the default: no `raw`, no file, a truncated .gz, a cut JSON, a record that is no object;
`any_json` returns whatever JSON the file holds, which renormalisation counts as read, as before), raw_record_answer
(`answer`, a dict) and raw_record_body (`answer.body` as it is). Each caller keeps its own guard (no run directory, no
read) and return shape. What changes, only where today's readers would have read no file or crashed: a caller that
probed one place now also probes the other layout's (compare and the document under run_dir also try
run_dir/arms/<arm>/<raw>; renormalisation and previous_entry under the arm directory also try the run directory for an
'arms/...' path); bench/analyze probes run_dir/<raw> before run_dir/arms/<arm>/<raw> (no run directory has a top-level
raw/ in either layout, so the file is the same); a record that is no object reads as none where a copy raised
(analyze's `(doc or {}).get`, the document's AttributeError catch). A test holds that no module but suite_common opens
a raw record (tests/test_common RawRecordReader). run_direct.py and suite_common.py are snapshot files: the change
reaches new runs (a run's snapshot pairs the two); recorded runs keep theirs.

### absence_scope and strands via MCP (item 5; the owner's scope addition of 2026-10-18)
With the service bundle that carries the identity, the MCP's /pattern answer states per by_database[] entry the
host answer's `absence_scope` (verbatim, a string: "any_offset", "suffix_only", "long") and `strands` (the host's list
verbatim: ["+","-"], ["+"], ["-"], ["="] for a palindrome; some fixtures ["forward","reverse"]); for a database with
several answers (graph chunks) both are null there and each graph's values sit in by_graph[]; a key is absent when the
host's answer states none. The entry's own absence_scope (the service's merged {scope, tasks_failed}) and the predicate
block's strands ("either" / "context") are other fields.
- **The checks** (checks.MCP_PATTERN_DATABASE_FIELDS, _mcp_database_values): on the MCP arm, entry_fields reads these
  two fields from the entry's by_database[] items (each database's own non-null value), else their by_graph[] items,
  and the entry's own by_graph[]; with neither list, the entry's own value. The LIVE PATH is the last one: the MCP
  runner's normaliser (scripts/bench/suite_mcp_normalise.py, DB_ANSWER_FIELDS, from runner fcf7ae49) puts the one
  database's host string and list onto the host-shaped entry, never copies the service's merged absence_scope object,
  and leaves a chunked database's per-graph values out (listed as left out: "several answers: per graph in by_graph";
  refseq33m-experimental is one graph); the by_database[] / by_graph[] reading serves a normaliser that keeps the
  per-database view (none does today). On the entry a value of another type is a wrong value and fails, never "not
  exposed" (the v7 verification's low 1). Every stated value is held to the expectation: a present but different one
  fails ("absence_scope = 'any_offset' (by_database[0]), expected 'suffix_only'"), none stated is not exposed, with
  the same detail as before ("field not in the answer: absence_scope"). A direct arm reads the entry as before. The
  recorded v5 / v6 MCP answers carry neither field and their manifests keep the marks, so they judge exactly as
  before (with runner fcf7ae49 and a service carrying the fields, a v6 run's entry_fields would compare them too: a
  mark sets aside only a not-exposed reading; the v7 verifier injected the direct values into the 88 recorded MCP
  executions and every one passes). An older runner drops the fields, and the v7 rtp gate then blocks on them.
- **The manifest**: areas.MCP_PATTERN_ENTRY_CARRIED_V7 (absence_scope, strands) is left out of what the MCP's answer
  cannot carry for v7 (mark_pattern_mcp_fields(carried=)), so v6's not-implemented marks of entry_fields on
  P-exact-0001..0008, P-scope-0001 / -0002 / -0003 / -0005 / -0007..0012, N-usecase-0005, D-P-exact-0001 and
  D-P-scope-0004 are gone, and P-scope-0004 / -0006 (absence_scope their only property check, scenario-level marks in
  v6) get back v5's MCP mapping (pattern_search with their arguments at 110,000 ms, exposes counts_consistent,
  invariants, status). P-edge-0004 (item 1) is the one /pattern mark left. counts.mcp_not_implemented 53 -> 51,
  mcp_exposure_baseline.scenarios_with_exposes 510 -> 512.
- **The policies follow from the marks**: the rt policy is v6's but for its pins; the rtp policy requires on mcp 630
  scenarios / 1,280 executions (v6 628 / 1,275: P-scope-0004 and -0006 r0 and r1, D-P-scope-0004 r0, which
  marks_cover_repeat no longer excludes), REL-subset-0002 / -0003 on both arms (19 relations required, as before; 8
  excluded instead of 10), 91 exclusions (manifest_not_implemented 48, no marks_cover_repeat), direct-local unchanged
  (686 / 793).

### The v7 manifests and policies
- `manifest/suite_v7.json` and `manifest/suite_smoke_v7.json`, built on this tree with `SUITE_TRUTH_CACHE=<the record
  cache v4 was built from> build_manifest.py --no-fetch` from the capabilities snapshot of release refseq97-6f945db4
  (deterministic: two builds, one content): manifest_version 7, 725 scenarios, 767 executions per arm, content_sha256
  e89708cb795e098225d76ac8a69e75a116f9029549a8613c96325e8ee6a36238, relations
  39660cfa0fb98cc47ac923fad3aa33419a84a3d92d662154cf2d4164752f4a86 (v5's and v6's); smoke 14 scenarios, content_sha256
  1ef0a151aecf8e5ea58f644165561d1e41ec0328ecf115a40d3f48feb5544ced. The field-by-field diff v6 -> v7: ids, order and
  relations equal; 454 scenarios changed, 271 byte-identical: `expect[budget_respected].args.reserve_ms` 250 added on
  the 454 POST /pattern scenarios that state the check (item 2); `mcp.not_implemented_checks` removed on the 14
  absence_scope and 7 strands scenarios (item 5) and on P-edge-0004, which gains `mcp.by_design_checks`
  [{entry_fields, pattern.lowercase}] (item 1); P-scope-0004 / -0006's `mcp` back from the scenario-level mark to
  pattern_search with their arguments and exposes (item 5); no other field, no scenario of another route. counts:
  mcp_not_implemented 53 -> 51, the rest v6's; mcp_exposure_baseline.scenarios_with_exposes 510 -> 512. The smoke
  diff: reserve_ms on its 6 /pattern scenarios and P-exact-0001's mark gone.
- `gate/declared_v7_rt.json` names v7's content and declares v6's one entry (REL-subset-0007 on mcp).
- `gate/policy_v7_rt_mcp.json` (canonical sha256 2bbed2e85b71...) and `gate/policy_v7_rtp_mcp.json` (8a8c9a1c8d7e...),
  generated as v6's (`gate_policy.py manifest/suite_v7.json --name v7-rt-mcp --arms mcp,direct-local --declared
  gate/declared_v7_rt.json --out gate/policy_v7_rt_mcp.json`; the rtp one with `--name v7-rtp-mcp --routes "POST
  /resolve,POST /traverse,GET /traverse/capabilities,POST /pattern"`). The rt policy equals v6's (and so v5's) but
  for its pins: 241 selected / 36 excluded, required on mcp 205 scenarios / 416 executions, on direct-local 237 /
  273, 5 relations, and its 237 scenarios byte-identical to v6's with their indexes (test_builder_v7 FrozenV7). The rtp
  policy equals v6's but for its pins and item 5's consequences: the same 721 selected in the same order; required on
  mcp 630 / 1,280 (v6 628 / 1,275: P-scope-0004 and -0006 r0, r1; D-P-scope-0004 r0), on direct-local 686 / 793 as in
  v6; 91 exclusions (manifest_not_implemented 48, no marks_cover_repeat); 19 relations required with REL-subset-0002 /
  -0003 on both arms, 8 excluded. P-edge-0004's class (item 1) and the reserve argument (item 2) change neither policy
  (the gate applies the manifest's marks itself; the policies hold no check arguments).
- run_suite.sh keeps its defaults (v5); a v7 run names `--manifest manifest/suite_v7.json --gate
  gate/policy_v7_rtp_mcp.json` and, from the LIVE line, `--expect-service-commit SHA`. Its MCP arm needs the service
  bundle that carries the identity and absence_scope / strands per database (otherwise the rtp gate blocks on the 22
  entry_fields properties read not exposed and on the missing identity).

### Validation of v7 (offline; no host was contacted; the run directories were only read)
- Copies of the v5 baseline (20261015T161813Z_refseq97-81c620ef), the v6 baseline (20261017T082918Z_refseq97-81c620ef)
  and a v6-rtp gate run (20261017T192446Z_refseq97-81c620ef), each re-evaluated (reevaluate, merge --strict-identity,
  compare, gate) with the code of 48e04d825 and of every commit of this lane: 0 changes of any field of the
  re-evaluated lines and of the comparison; the v5 copy GATE PASS under policy_v5_rt_mcp (689 judged, 635 passed, 0
  problems), the v6 copies GATE PASS under policy_v6_rt_mcp (689 / 635) and policy_v6_rtp_mcp (2,068 / 1,950), 0
  problems; run.json gains `service_identity` (recorded false) and each verdict one note ("the MCP arm records no
  service identity ... not checked"). The recorded sets (merge, compare, gate without re-evaluation) likewise.
- The raw-record readers (item 4) on the same three runs: every reader at its call form on every line, base and head,
  77,517 calls, 0 differences.
- With the deployed commit given (--expect-service-commit 60e737f7...), the v6 baseline's copy is blocked by one
  identity problem: its runner recorded no service identity.
- 964 tests pass, none skipped (with SUITE_TRUTH_CACHE, the record cache v4 was built from; 930 at 48e04d825); each
  new rule's test fails on the code before it. manifest/suite_v6.json recomputes to 49907dbb...; every manifest,
  declared file and policy of v1-v6 and service_by_design.json are byte-identical to 48e04d825's.

## Manifest v8, lane L: area G's local library runner and the selectors (2026-10-19)
Design: DESIGN-v8 §2-§4, §7 items 2-3; the contract FORMAT-v8 §1-§3 (both in the v8 plans, outside the repository).
Nothing here touches staging; every test is offline.

### Files
- `graphlet_doc.py`: the suite's own digests of an MGT document (no library import): `body_digest` (H with its
  seed_index written `*`, J and Z left out: a server body, its standalone file and one seed at two batch positions
  share it; unmasked it is the library's parser.body_digest), `envelope_digest` (the J envelope without the top-level
  `timing` and `usage`, results[*].`timing` and, FORMAT v2 E1, results[*].annotation's lookahead counters
  `direct_reads` and `keys_mapped`; strategy settings such as output.timing and annotation.batch_kmers stay),
  `digests` (FORMAT's `doc` block). A broken frame (no final LF, no H of MGT v1 first, J anywhere but line 2, no Z counting the lines) is
  a ValueError, never a digest.
- `graphlet_select.py`: run-time selectors, pure and stdlib (async vendors it): `resolve(value, answers)` replaces
  `{ref, path}` (an RFC 6901 pointer into an earlier assembled answer), `{select, args}` and `{cursor_of}` at any
  depth; a label
  selector `{"ref": "c:3"}` (no `path`) stays a literal. Named: walk_through_merge, walk_with_continuation,
  slice_past_end, label_entered_at_merge, top_label, first_arm, second_arm_or_none, each defined in its docstring
  over the library's answer shapes. A failure is a `SelectorError` naming the argument, the op and the rule.
- `graphlet_paging.py` (FORMAT v2 A2): the paging loop's helpers (`paging_problem`, `max_pages`, `next_call`,
  `more`, `loop_error`, `error_page`, `error_record`, `paging_complete`) and `assemble(tool, pages, complete,
  volatile=())` with its tables, pure and stdlib: async's runner vendors it, so both sides follow and assemble pages
  identically by construction. `loop_error` is the `paging_runaway` rule (VERIFY-L VL6): a cursor or next_from that
  repeats one already asked for stops the loop at once (reason `repeats`), paging "all" stops at 10,000 pages
  (reason `max_pages`); a `pages` cap reached is no error.
- `run_graphlet_local.py`: the ways W1 (arms/dl-lib), W2 (arms/dr-lib), W4 (arms/mcp-lib) and drift
  (arms/drift-lib) over a run directory and a manifest, or offline over `--doc [ID=]PATH ... --plan PLAN`;
  `--export-library COMMIT --to DIR` (`git archive COMMIT metagraph/api/python | tar -x`, DIR/LIBRARY_COMMIT);
  `--drift-report` (the drift way's differences from W4, report-only). The module docstring states the record and
  file layout.
- `tests/test_graphlet_local.py` (the digests against the library; an op plan calling all 19 tools with every paging
  mode and an error case each; paging assembly against the uncut answer for every list tool and the sequence; hash
  seeds; budgeted pages resumed under local limits; a response against the same document from a file; the manifest
  mode; export and drift),
  `tests/test_graphlet_select.py`, `tests/test_graphlet_paging.py`. The library-backed classes skip below Python 3.10
  (the library's own floor).

### Decisions where FORMAT-v8 is silent
- Documents: W1/W2 take `answer.body` of the direct raw record and convert it in the way's process with
  `parser.standalone_text(results[0].graphlet, results[0], response)`; W4 and drift read the MCP's exported file
  `arms/mcp/docs/<id>.mgt.gz` (FORMAT v2 C1; `--mgt-template`); every way writes the document it processed to
  `arms/<way dir>/docs/<id>.mgt.gz`, then to its spool, and loads it by `graphlet_load`. A missing
  file, a failed seed or a refused document is the error `no_document` / `no_graphlet` / `document_refused` of every
  op of its scenario.
- Isolation: one child per way, `python -s -B [-P]`, its environment scrubbed to PATH/HOME/TMPDIR/LANG/LC_* plus
  PYTHONHASHSEED (defaults W1 1001, W2 2002, W4 4004, drift 5005; distinct per run), sys.path the library root first.
  Not `-I`: it implies -E, which ignores PYTHONHASHSEED. The child refuses a library imported from elsewhere; the
  parent refuses a child whose seed did not take effect.
- Handles in records are symbols relative to the op's scenario: `<handle>` (its own graphlet, FORMAT's form),
  `<handle:S>` (scenario S's), `<handle:OP>` (made by op OP: subtrie, load), `<handle:S/OP>`, `<handle:?>`; cursors in
  args_resolved `<cursor>`; the way's directories `<export_dir>`, `<spool_dir>`. Raw pages keep everything.
- Assembly (graphlet_paging.assemble): every paged library list is `rows` (mcp_tools._page / _page_cut; the table is
  checked against the source by a test). Rows are concatenated; every other field comes from the first page that
  states it (page 0's, and a key page 0 lacks from the first later page carrying it: nothing a page states is
  dropped, VERIFY-L VL1); the fields that say how a listing ended (total, total_at_least, complete, filtered,
  filtered_so_far, hint, local, stop) from the last page; `evidence` from the last page that states it (a page its
  budget stopped marks it `limitations.local`); compare's `fields_cut` then names only what no page carried;
  `row_truncated` / `cut_fields` merged and every row of a page that cut its row marked
  `"_cut_fields": [that page's cut_fields]` (FORMAT v2 A3); next_cursor and cursor removed. `volatile=` names keys a
  runner's pages may state differently (taken from the last page stating them): the local runner passes none, async
  the service's own per-page fields (VERIFY-L VL4). graphlet_sequence concatenates
  `sequence`, takes `to` from the last page and drops next_from (and `truncated` when the loop reached the end).
  A base field two pages state differently is `assembly_error`.
- Records add `paging_complete` (FORMAT v2 A4, graphlet_paging.paging_complete: false when the last page still
  names a next one or the answer states complete false; null on a record without an answer) and, for
  graphlet_export / graphlet_save, `file` {name, bytes, sha256, raw} (FORMAT v2 B3: the file's text kept gzipped at
  `raw`, relative to the way's directory, with handles and directories as symbols like the rest of the record;
  `bytes` the file's own size, `sha256` the kept text's); `library` adds `limits`; `error.fields` leaves out `message` and `hint`. Runner-side error codes:
  exception (the library raised), unknown_tool, paging_runaway, assembly_error and the three document codes. A plan
  names a cursor by `{"cursor_of": "<op id>"}` (graphlet_select, VERIFY-L VL5): the last next_cursor of that op's
  listing when its paging stopped early, kept by `graphlet_paging.selector_view` in the answers selectors read.

## Glossary

The full vocabulary of the suite and of its results document (`doc/make_doc.py`), one entry per term, in alphabetical order (a leading article ignored). The document's appendix D carries a one-page selection of about 25 of these terms and points here for the rest, so this section is **version-pinned**: it was written for the checks of commit c11db9dd and manifest version 4, and a results document names on its title page the checks commit and the manifest version it was built with; when the rules change, change the entry here in the same commit as the code, and the pin below.

Pinned to: checks commit `c11db9dd`, manifest version 4 (content sha256 `3a7417d60e4eb5d0…`), by-design table version 1 (service commit `22d034ba`). Definitions follow `checks.py`, `compare.py`, `report.py`, `gate.py` and the SPECs they cite (SPEC-labeled-traversal-core section 4 for /resolve). Where a figure is the run's (the warm window of 60 s, the Logan thresholds of 400 % `sra_all` CPU and load1 48, k = 31) the document states the run's own value.

- **agreement depth** — Of an agreeing MCP pair, how much was compared: fully checked (every property check the scenario expects was evaluated on both arms and agreed), partial (the MCP answer exposed only some of them; counted apart when the aggregate is not a pass), generic checks only (none of them), a refusal (both arms refused alike), on failure (both arms failed alike). The chapter lines of the results document count the fully checked pairs. *See:* MCP agreement class.
- **answer digest** — The sha256 of the canonical JSON of the normalised answer without its headers and without every key named timing (and, for a capabilities document, without its measured and observed figures); the determinism check and the direct agreement read it. *See:* normalised answer, direct agreement class.
- **answer-triggered case** — A by-design case matched by the MCP answer itself rather than named by the scenario (e.g. the host's /pattern 503 deadline relayed by the service as a failed search): the execution is by design, verified when a direct execution matched the case's host side. *See:* by-design table.
- **area, chapter** — The ten groups of scenarios by what they test, one letter each (R resolve, T traverse, P patterns up to k, L patterns longer than k, X peptides, M mismatches, Q predicates and motif, C capabilities and protocol, N negatives and near-misses, D determinism); Part II has one chapter per area. *See:* family, unit.
- **arm** — One of the three routes to the host a scenario is executed on: direct-local, direct-remote and mcp. Every table shows the three side by side. *See:* direct-local, direct-remote, mcp.
- **backend time** — The host call as the search service measured it (`mcp_phases.backend_ms`), on an MCP cell after be; the MCP overhead is the MCP wall minus each direct arm's wall on host-warm pairs, with and without the queue. *See:* mcp_phases, wall time.
- **basic graph** — The index's graph mode: it holds the forward strand of every indexed sequence only (no canonical k-mers), so the reverse complement of a window is another set of k-mers, which the record's labels do not support. *See:* k-mer.
- **by-design table** — The search service's `service_by_design.json` (131 cases at service commit 22d034ba): for each case the request it is about (`scenario_match`), the host's expected answer, the service's answer (outcome `refused_at_submit`, relayed, answered, clamped or `not_expressible`; code, text, stage) and the checks an answered case explains. Scenarios name their cases; an answer may also match a case (answer-triggered). Appendix B lists the cases of this run. *See:* service-side by design, answer-triggered case.
- **check** — A pure function of a scenario and a normalised answer, in checks.py, shared by every arm; it answers pass, fail, inconclusive, finding or not exposed with a detail. Appendix A describes every kind this run uses. *See:* generic check, property check, check outcome.
- **check outcome** — pass; fail; inconclusive (the check could not test its property on this answer); finding (a stated bound was broken on the clock the check reads: a 200 past its budget by the host's own elapsed time, a walk stopped outside its stated rule on the seed's clock, an allowed 503 deadline more than max(1 s, the budget) late on the client's clock, which carried no results); not exposed (the MCP answer lacks the field the check reads: recorded, not compared). Neither inconclusive nor finding counts as pass. *See:* execution outcome.
- **clamped** — The host lowered a request bound to its cap and stated it in limits.clamped (or strategy.clamped) with the requested and the effective value; the check of the same name reads that entry. *See:* time budget.
- **clean** — Of a test: every arm is a pass, not implemented via MCP, a verified by-design pair, or had no execution because the plan had none. A clean test costs one cell or one line; anything else carries a footnote mark and a detail block. *See:* detail block.
- **code-derived** — A tag on an expectation read from the C++ source or the fixtures that no smoke run had exercised at build time (an exact refusal text): a first failure there may be a suite error rather than the host's. *See:* expectation.
- **cold / warm** — In a cell, the executions in repeat order; an MCP pair reads host-cold then host-warm, a direct cell without a mark is host-warm, a direct time with the mark c was host-cold (a was aged). Family and chapter lines give the wall p50 and p90 for host-cold and host-warm executions apart. *See:* host cache state.
- **column label** — A taxid, one column of the annotation; a column supports a k-mer when the k-mer is in the graph and the column is in its annotation row. A column may hold several records. *See:* header label, support.
- **connect time** — The set-up of a fresh connection (TCP, and TLS on the public URL), recorded per direct execution apart from the wall time and reported per arm as part of the front cost. *See:* the front, wall time.
- **deadline (503)** — The 503 with code deadline that /resolve and /pattern answer when the answer cannot be built and written within the budget: nothing partial is sent. Where the request names a budget the SPECs allow it on any host under load, so the suite reads it as inconclusive (and as a finding when it comes more than max(1 s, the budget) late); on a /resolve without a budget it is a failure. *See:* inconclusive, finding, machine-dependent.
- **declared file** — The owner's versioned list beside a policy: exclusions (a scenario on an arm, with a reason), checks not required (a property the MCP surface never exposes) and relations excluded. The only way past the gate for a scenario the rules do not exclude; Part I prints it. *See:* policy, deploy gate.
- **deploy gate** — The strict judgement of a run against a policy (gate.py): it fails on any identity problem, any missing planned or required execution, any fail, error, finding or inconclusive on a gated arm, a mark the policy does not explain, a property check not passed, any disagreement between the gated arms, any required relation not passed, and configuration evidence not recorded live. The v4 gate runs two arms, mcp and direct-local, on the scenarios of POST /resolve, POST /traverse and GET /traverse/capabilities. *See:* policy, declared file, verdict.
- **detail block** — The explanation printed after a unit's table for every test that is not clean: the mark, the id and title, the route, timing class, tags and repeats; the expected line; per arm that is not clean its outcome, status, error, times and every check that did not pass with its detail; the agreement classes; the relations involved; and for a finding the rule it broke. *See:* clean, footnote mark.
- **direct agreement class** — How direct-local and direct-remote compare on a scenario: identical (the same answer digest and no check differing other than through the execution's own clock or repeats), same checks (another answer, the same check outcomes), timing (the outcomes differ only because one arm's answer states a time limitation the other's does not, and the stopped answer is consistent with the complete one), not comparable (no check evaluated on both), disagree. A disagreement fails the run. *See:* MCP agreement class, answer digest.
- **direct-local** — The runner on mex sending the request to the container's own port: no proxy, no network; the cleanest server timing and the reference arm of the MCP comparison. *See:* arm, reference arm.
- **direct-remote** — The runner on the laptop sending the same request through the public URL and the nginx front on dnaloc-new: the network and the front are in its wall time. The front cuts at about 120 s, so the long timing class is not run there. *See:* arm, the front.
- **discovery** — A /resolve without labels (discover: kind column or header, `max_labels`): the labels of every in-graph k-mer's row, ranked by supported k-mers, the top `max_labels` kept (ties by column id, then sequence id) and `labels_truncated` stated when more occurred than were kept (kept, total, the least support kept, the most dropped, the dropped labels at full length). *See:* support.
- **error** — An execution with no HTTP answer (a transport failure, a client time-out), an undecodable body, or a check that raised. Counted with the failures. *See:* execution outcome.
- **execution** — One request of one scenario on one arm, one repeat: its HTTP status, wall and server times, host cache state, checks and outcome, and the raw answer. One line of executions.jsonl. *See:* repeat, execution outcome.
- **execution outcome** — From the checks, by precedence: fail, finding, inconclusive, pass; inconclusive also when no check could test anything or every property check could not. Decided without the checks: error (no HTTP answer, an undecodable body, a check that raised), not implemented via MCP, service-side by design. Any outcome outside this vocabulary reads as a failure. *See:* check outcome, scenario outcome.
- **expectation** — One entry of a scenario's expect list: a check name and its arguments. The builder validates every expectation against the checks' argument table, so a misspelt argument cannot weaken a check. *See:* check, generic check, property check.
- **family** — The scenarios that share an id prefix (R-pos, T-walk, N-hard), built by one part of the builder around one question; a family has a title and one paragraph of prose and holds one or more units. Chapter D's families are all D. *See:* unit, area, chapter.
- **feature level** — What the server offers beyond the base contract, monotonic (SPEC-labeled-traversal-core section 10.3); stated in every /resolve and /traverse answer and recorded with the run's environment. *See:* release.
- **finalisation reserve** — The 250 ms a /resolve keeps of its budget to build and write the answer: the work stops at the budget less the reserve; a budget not above the reserve is refused. *See:* time budget, stated stop.
- **finding** — A stated bound broken on the clock the check reads: a 200 past its time budget (the host's own elapsed time), an allowed 503 deadline more than max(1 s, the budget) late (the client's clock: the suite's tolerance for a deadline answer, which carried no results and does not show where the time went), a walk stopped before its budget or run past it beyond its pieces (the seed's clock). Counted apart from failures; a run with findings and no failure reads FINDINGS. *See:* check outcome, verdict.
- **footnote mark** — The superscript digit on a cell that points at the test's detail block after the table; numbering restarts in every chapter and a test's cells share one mark. *See:* detail block.
- **the front** — The nginx reverse proxy on dnaloc-new that the public URL reaches; every arm but direct-local goes through it (the service too). The front cost is direct-remote minus direct-local on pairs of host-warm executions, the connection set-up apart. *See:* direct-remote, connect time.
- **generic check** — A check that holds on any answer of its route whatever the scenario asks: status, invariants, `counts_consistent` (/pattern) and `budget_respected`. The builder adds `budget_respected` to every scenario that names a time budget and allows a 200. *See:* property check.
- **ground truth** — What a scenario's right answer rests on: the record a window was cut from (its accession.version, taxid column, position in the column), the window's position and strand on the record's forward strand, the record bases around it, and for variants the edits made; for negatives the fact that the sequence is random or shuffled. Records are the 32 genomes of the search service's validation set, fetched at the indexed accession.version. *See:* record, window, last-record-verified.
- **header label** — An accession.version, one indexed sequence of a column, reached through the index's coordinates: it supports a k-mer when some coordinate of its column at the k-mer's row maps to that sequence. *See:* column label, support trace.
- **host cache state** — A recency heuristic over all arms, per execution: cold (the scenario's first execution to reach the host in this run), warm (started within 60 s of another execution of the scenario ending on the host), aged (longer after it: unknown), not reached (the host never saw it: a refusal at submit, a request the MCP cannot express), unknown (no start time). Everything timed keys on it; cache residency is never measured. *See:* cold / warm, warm after a cut, paired loop.
- **host state** — An MCP execution the service ended by one of its own limits before the host answered (a `resolve_timeout` at the service's deadline for the host call, the budget plus the reserve plus an allowance; a relayed `traversal_resolve_deadline`; a relayed /pattern deadline): inconclusive with that reason when the service's own timing shows the limit can have been reached (the host-state rule of the MCP runner, replayed by the re-evaluation), a status failure otherwise. The checks as they read on the raw answer are kept in the record as the intermediate evaluation; the results document states the adjudication first. *See:* inconclusive, relayed deadline.
- **inconclusive** — The test did not test its property: an allowed 503 deadline, a list a legitimate statement hides the source in, a walk a stated resource limit stopped first, a stop that was the property and did not happen, a determinism check with nothing to compare, an MCP answer exposing none of the property's fields, or an MCP execution the service ended by one of its own limits (a host state). Never a pass, never a failure; every inconclusive execution blocks the deploy gate. *See:* deadline (503), host state, deploy gate.
- **index fingerprint** — `index_fp`, the hash that identifies the loaded index; the manifest records the fingerprint it was built against and the deploy gate pins it. *See:* release, deploy gate.
- **intermediate evaluation** — The check lines of an MCP execution as they read on the raw answer before the host-state rule overrode the outcome; kept in the record and printed in a detail block under that label, after the adjudication. *See:* host state.
- **k-mer** — A substring of k = 31 bases; the graph's nodes. A query of L bases has L - k + 1 k-mers, numbered by their start, 0-based. On this basic graph a k-mer and its reverse complement are different nodes. *See:* basic graph, run.
- **last-record-verified** — A tag on a scenario whose record is the last of its taxid column, where the index defect of 2026-10-09 (a last record present in the sequences but unannotated) could strike; the truth records the scan's result (ok, or single-record column), and a failure there is listed apart. *See:* ground truth.
- **load telemetry** — docker stats of the staging container (CPU, memory) and of the `sra_all` container (the Logan indexing load) and the host load average, sampled every 2 s on mex and joined to each execution's window; contextual figures of the shared host, not the request's own cost. *See:* the Logan caveat.
- **the Logan caveat** — The statement that the run's timings were taken while mex carried the Logan load: given when the `sra_all` CPU mean exceeds 400 % or the load average mean exceeds 48 (half of mex's 96 cores) over the run's samples. Such timings are labelled and are not the timing baseline. *See:* load telemetry.
- **machine-dependent** — A tag on scenarios whose answer depends on the host's load (time budgets: a 200 or a 503 deadline, a stop or none). For /pattern and /resolve the checks allow the 503 deadline on every scenario that names a budget and expects a 200, tagged or not. *See:* deadline (503).
- **manifest** — The frozen, versioned set of scenarios (`manifest/suite_vN.json`), identified by its content sha256 over the canonical JSON of the scenarios; its relations are hashed apart (relations sha256). A run records both; a run on another manifest is another run. *See:* scenario, relation.
- **mcp** — The same request asked through the search service's MCP tools on its staging deployment (`traverse_resolve`, `traverse_fetch` and `graphlet_export`, `pattern_search`, the capabilities tools); the service validates and reshapes it, runs it as a job in its queue and calls the host over the public URL. Its answer is normalised to the host's shape for the fields the service exposes. The service's deadline for a resolve's host call is the budget plus the reserve plus an allowance (30.25 s for a 20,000 ms budget); the tool's inline wait (30 s with explicit labels, 10 s for a discovery) only decides whether the answer comes inline or as a job handle. *See:* arm, mcp_phases, backend time.
- **MCP agreement class** — How the MCP arm compares with the reference direct arm on a scenario, check by check and repeat by repeat on the checks both arms evaluated: agree (with its depth: fully checked, partial, generic checks only, a refusal), not comparable (no check compared), disagree, not implemented via MCP, service-side by design (verified or not), missing (no MCP line). A disagreement fails the run. *See:* reference arm, direct agreement class.
- **mcp_phases** — The decomposition of an MCP execution's wall time: submit (the tool's acknowledgement), queue (job created to started), backend (the host call), terminal lag (finished to the client seeing the terminal state) and paging (fetching the results). The closing section tabulates them. *See:* backend time.
- **mixed** — A scenario whose repeats on one arm have different outcomes (e.g. inconclusive on the cold repeat, pass on the warm one); the table shows the worst glyph with the mark m and the detail block each repeat. *See:* scenario outcome.
- **normalised answer** — The shape every check reads: the HTTP status, the body, the error (code, text) and, on a direct arm, the headers; the MCP runner maps the tool results into it for the fields the service exposes. *See:* answer digest, check.
- **not exercised** — The results document's wording for an MCP cell, aggregate or count whose executions the tools could not express (agreement class not implemented via MCP): never ok, never passed, since no property was tested there. *See:* not implemented via MCP.
- **not implemented via MCP** — The request cannot be asked through the MCP tools: the manifest marks it (`mcp.not_implemented` naming the field: support, `bounds.max_query_bp`, an unknown field, a protocol body), or a case of the by-design table with that agreement names the missing argument (e.g. a 10 kb resolve: the service always sends its own 2,000 bp cap). The results document prints such a cell as *not exercised*, never as ok; the service's own refusal, when one was recorded, is under the row or in appendix B. A known limit of the tools, not a failure and not coverage. *See:* service-side by design, by-design table.
- **paired loop** — The run's order per scenario: MCP cold, direct-local, direct-remote, MCP warm, one request in flight against the host at a time. It makes the direct executions host-warm by construction wherever the MCP cold call reached the host. *See:* host cache state, cold / warm.
- **policy** — The deploy gate's pinned selection: the manifest's content and relations hashes, the index fingerprint and the by-design table's hash; the scenarios by route; the arms; the exclusions derived from the manifest and the table (a request the MCP has no argument for, an advanced traverse override) with the owner's declared entries; the executions and relations each arm must record. *See:* deploy gate, declared file.
- **prefix** — The part of the query a stated stop resolved: its first `resolved_kmers` k-mers; a client resolves the rest from base `remainder_from_bp`. The relation `prefix_equals` compares a stopped answer with the full one on their common k-mers. *See:* stated stop, relation.
- **property check** — Every check that is not generic: the claim the scenario was built to test (`resolve_runs`, `contains_source_context`, `walk_in_record`, ...). An execution whose property checks all could not test anything is inconclusive even when its generic checks passed; a scenario with no property check (a refusal) is judged by its generic checks. *See:* generic check, inconclusive.
- **property line** — The one sentence above a unit's table that states the claim every member of the unit tests; the family prose before the units says what the family checks and why. *See:* unit.
- **re-evaluation** — Judging a recorded run again offline with the current checks, without any request (reevaluate.py): every recorded answer is re-checked, the MCP runner's by-design and host-state decisions replayed; outcomes a runner decided without the checks (error, not implemented, by design) are kept. This document reads the re-evaluated set when it exists. *See:* check.
- **record** — One indexed sequence, named by its accession.version, which is also its header label; it belongs to a taxid column and has a position in it (`seq_id`). The suite's records are the 32 genomes of the ground truth. *See:* column label, header label.
- **record key** — The three- or four-letter abbreviations of the ground truth's records in the tables (eco, mtb, pfa, ...), keyed in appendix D from the manifest's ground_truth.records (key, organism, accession.version, taxid). *See:* record.
- **reference arm** — The direct arm an MCP execution is compared with: direct-local, or, for a check item direct-local left undecided (inconclusive or not exposed) while direct-remote's item of the same repeat is decided, direct-remote's. *See:* MCP agreement class.
- **relation** — A comparison of two scenarios' answers on the same arm, evaluated by the comparator on the first repeat of each: `same_answer` (equal bodies without timing, or the same refusal), `prefix_equals` (two resolves of pieces of one sequence agree on their common k-mers), `zero_equals_absent` (a /pattern count of 0 exactly when the /resolve of the same sequence misses a k-mer), subset (every result of one entry among the results of another). Its outcome is pass, fail, inconclusive (the pair cannot test the claim: a non-200, no common k-mer) or not exposed; a failed relation fails the run. *See:* prefix.
- **relayed deadline** — The host's own 503 deadline coming back through the service under a service code (502 `search_failed` with a `pattern_deadline` text, 503 `traversal_resolve_deadline`, 504 `resolve_timeout`): the host did work, so such an execution counts as having reached the host for the cache states, while a direct 503 deadline is the host's own answer. *See:* host state, host cache state.
- **release** — The engine's release string as its capabilities state it (refseq97- and the engine's commit hash); a run that spans two releases tested no single one and never reads PASS. *See:* server instance, index fingerprint.
- **repeat** — The index of an execution among the scenario's executions on one arm: repeat 0 is the arm's first, later repeats follow within the paired loop. Every tenth scenario and every chapter-D scenario has direct repeats; an MCP scenario the service can express has two (cold, warm). *See:* execution, cold / warm.
- **run** — In a /resolve answer, a maximal interval of consecutive k-mers with one state, written as a 0-based half-open interval of k-mer starts [a, b): `graph_runs` are the runs of k-mers in the graph, a label's runs those it supports. *See:* k-mer, support.
- **run id** — The run directory's name, the UTC start time and the release's short hash (e.g. `20261014T061053Z_refseq97-6f945db4`); the arms of a run share it. *See:* release.
- **scenario** — One test of the frozen manifest: an id (the family prefix and a number, e.g. R-pos-0013), an area, a title, a route, the exact request body, the ground truth, the expectations (check name and arguments), the MCP mapping, a timing class, a repeat count and tags. The test id of the tables is the scenario id. *See:* manifest, expectation, unit.
- **scenario outcome** — Over a scenario's repeats on one arm: the common outcome, else the worst (error over fail over finding over inconclusive); repeats that differ otherwise are mixed, shown with the worst glyph and the mark m. *See:* execution outcome, mixed.
- **server instance** — The identity of the running server process as its capabilities state it; a change during a run means the server restarted. *See:* release.
- **server time** — The host's own `timing.elapsed_ms` in the answer (the top-level one, else the entries' sum); in brackets under the wall time. A refusal and a 503 carry none. *See:* wall time.
- **service-side by design** — The MCP arm legitimately differs from the direct arms because the search service validated or reshaped the request before the host saw it, or refused it in the host's place, as a case of the by-design table documents (a refusal at submit, a relayed host refusal under the service's code, an answer where the host refuses, a clamped value). A pair is verified when the direct answer matched the case's host side and the MCP answer its service side; otherwise it carries the mark u. *See:* by-design table, answer-triggered case, verified.
- **standing** — In the results document, what a cell shows: the execution outcome, or on the MCP arm the agreement class where that class states the arm's standing (not exercised via MCP, service-side by design). Family and chapter counts are by standing (by glyph); the chapter page adds the recorded outcomes as the report counts them. *See:* clean, MCP agreement class.
- **stated stop** — A /resolve answer that stopped at its budget less the reserve and says so in its stop block (reason time, `resolved_kmers`, `query_kmers`, `resolved_bp`, `remainder_from_bp`): it is exactly the resolve of the query's first `resolved_kmers` k-mers, status 200. A stop the answer states but that cannot have happened (another reason, a stop before the budget less the reserve) fails the checks that met it. *See:* time budget, prefix, deadline (503).
- **support** — Of a label at a k-mer: the k-mer is in the graph and carries the label (SPEC-labeled-traversal-core section 4.2). /resolve reports it per label as runs and as `kmers_supported`, the number of supported k-mers. *See:* run, support trace.
- **support trace** — /resolve's support: trace: a header label supports consecutive k-mers only while their coordinates in that sequence are consecutive (+1 per step); a coordinate jump ends the run. Not expressible via MCP. *See:* support.
- **time budget** — `bounds.time_budget_ms` of a /resolve (`time_budget_ms` of a /pattern, `strategy.bounds.time_budget_ms` of a /traverse): the request's deadline, lowered to the server's cap and stated in limits.clamped when lowered. For /resolve it is opt-in: without the field there is no deadline. The suite's explicit-label resolves name 20,000 ms, as the service sends it. *See:* finalisation reserve, stated stop, deadline (503), clamped.
- **timing class** — fast (under about 5 s), slow or long (over 100 s): the client time-out of an execution (120 s, 600 s, 1,800 s) and whether direct-remote runs it (long is not run there, the front cutting at about 120 s). *See:* direct-remote.
- **timing legend** — The three meanings of the slash in the results document's tables, printed under every table: repeats in order in a series or list cell (mcp: cold / warm; direct: repeat 0 / repeat 1), median / maximum in a replicates cell, p50 / p90 in a summary line. *See:* cold / warm.
- **unit** — One builder template, one property, varied over one dimension: the tests printed as one table. Its type decides the layout: series (arms as rows, variants as columns), replicates (the same test over many records, aggregated per arm and group), refusals (edge inputs, a glyph per arm), list (the refusals shape with time cells) and pairs (twins, relations and chapter D's repeats, one row each). Every scenario belongs to exactly one unit; the index (appendix E) maps each id to its unit. *See:* family, property line.
- **verdict** — The run's one word: INCOMPLETE (an arm did not finish, a planned execution not recorded), FAIL (a failed or errored execution, a disagreement, a failed relation), FINDINGS (a finding and nothing worse), INCONCLUSIVE (nothing passed, a group that tested nothing, or two deployments in one run), EMPTY, UNVERIFIED (the status-only stub judged), PASS. Inconclusive executions are counted after it, never as passes. *See:* deploy gate, finding.
- **verified** — Of a by-design pair: every side the case has was checked and held (the direct answer against the case's host side, each MCP execution against its service side; an answered case through the checks it explains). An unverified pair keeps the by-design glyph with the mark u, and the appendix says which side was not verified. *See:* service-side by design.
- **wall time** — Of a direct execution: from sending the request to the last byte of the body on an established connection, the connection set-up excluded. Of an MCP execution: the whole tool interaction, from the first submit to the answer. *See:* server time, backend time, connect time.
- **warm after a cut** — A host-warm execution whose previous host execution of the scenario ended at the host's deadline without an answer (a direct 503 deadline or a relayed one): counted apart from the warm executions in the timing statistics and the overheads, its server time being typically higher. *See:* host cache state.
- **window** — A positive control: a stretch of a record's forward strand cut at a seeded position, of A, C, G and T only, with a DUST score below 2, occurring once on either strand of the fetched sequence, at least 50 bp from any other window. The suite's checks ask for the record among the answers, never for uniqueness in the index. *See:* ground truth.
