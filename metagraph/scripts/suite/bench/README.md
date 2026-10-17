# The isolated benchmark panel (lane B, 2026-10-15)

The suite's recorded timings describe loaded runs: v3 and v4 ran next to the Logan load (the `sra_all` container at
a median ~7,700 % CPU, the host's load1 ~227). They say what those requests took then; they are no release
performance baseline (GPT's review of the suite, "Benchmark conclusions and a costed replacement protocol";
VERIFY-A confirmed the cache wording). This directory defines the replacement: a fixed panel of 48 successful
queries, a run mode that measures them in a balanced, repeated, isolated way, a load guard that refuses a verdict
from a loaded host, and an analysis with uncertainty. It is defined now and runs the first time mex is free of the
Logan load. Nothing here has run against staging yet.

## Files

| file | what |
|---|---|
| `make_panel.py` | builds a panel from a merged, compared three-arm suite run (re-runnable: v5 replaces v4) |
| `panel_v4.json` | the panel: 48 queries chosen from the v4 run, re-evaluated offline with this branch's checks |
| `run_bench.py` | runs a panel against the staging host: warm-up, 6 cycles x 3 arms, one request at a time |
| `load_guard.py` | judges each block of a run clean, loaded or unguarded from the sampler's samples |
| `analyze.py` | per-query observations, paired ratios, fixed weights, bootstrap intervals, the gated verdict; A vs B |
| `bench_lib.py` | shared: result sizes, completeness, the run manifest, the Latin-square schedule, statistics |

Tests (offline, in `../tests`): `test_bench_panel.py`, `test_bench_guard.py`, `test_bench_analyze.py`,
`test_bench_run.py` (the fake host, `fake_mcp_bench.py`, a fake `docker` for the sampler, mex in a local
directory), fixtures in `bench_fixtures.py`. `python3 -m unittest tests.test_bench_panel tests.test_bench_guard
tests.test_bench_analyze tests.test_bench_run` from the suite directory.

## The panel (`panel_v4.json`)

Built by `make_panel.py --baseline <v4 copy re-evaluated with c11db9dd> --reevaluated --manifest
../manifest/suite_v4.json --name v4 --include T-derived-0003 --accept-incomplete "T-derived-0003=<reason>"` from
run `20261014T061053Z_refseq97-6f945db4`. Every query carries its frozen scenario (and its sha256), each arm's
recorded result sizes, the recorded walls (under the Logan load: ranking and costing only) and its weight.
`panel_sha256` binds ids, scenarios, weights, each arm's recorded result sizes and any accepted incompleteness
(panel version 2); `bench_lib.read_panel` refuses a panel whose queries do not match it.

**Eligible** (460 of 681 v4 scenarios): a search route (POST /resolve, /traverse, /pattern), no protocol case
(`body_raw`, `headers`, `run_unique`), not the long class, mapped to an MCP tool; every arm recorded it, every
execution passed (decisive: no inconclusive, by-design or failing outcome), answered HTTP 200 and was not time-limited;
the direct arms identical, the MCP arm agreeing, each arm's repeats with one digest; the direct answer a complete
search (`bench_lib.answer_incomplete`: no per-pattern refusal such as `information_below_floor`, no withheld
result, no stop or cut of any reason, no cut labels, occurrences or anchors, no incomplete retrieval outside a count,
no result whose labels are not complete, no truncated or stopped traversal arm, no label lists capped at
`labels.max_labels_per_node` (label evidence a lower bound), no truncated resolve labels: capped queries are budget
stress, which the protocol runs apart); not a duplicate request (the D area repeats other areas' requests; the
original is kept). Excluded counts are in the panel (`selection.excluded`, with the ids).

**One stated exception: T-derived-0003.** Its label lists are cut at the default `labels.max_labels_per_node` 64
(label evidence `lower_bound`), on every arm and every recorded repeat alike: a deterministic cap, so its work is
the same each time. It is kept (`--accept-incomplete`, tag `capped`) because it is the panel's only large traverse
answer (96 KB) and its most label-heavy walk (70 labels); without it the traverse mode has no large answer. The
query records the incompleteness each arm stated (`accepted_incomplete`), with the cap's value: the knob, its limit
and the count observed (`labels.max_labels_per_node limit 64 observed 68` on the left, `observed 69` on the right).
The panel identity binds it (panel_v4.json: `c94f7361`), and the analysis accepts exactly that incompleteness: any
other, including the same kind of cap at another limit or count (a release lowering the default cap), is
`incomplete`, and an uncapped answer is `result_changed` (another workload).

**Selection** per mode, 16 each, deterministic: the eligible queries ranked by recorded cost (the least wall of the
direct arms: load only adds time, so the least is the best estimate on a loaded host; the MCP arm's server-stated
elapsed time can understate, e.g. Q-strands-0003's MCP answers state 12-49 ms where the direct arms state and take
9.4 s) in 16 bins of equal width in log10(cost), so the panel spans the cost range rather than following the suite's
many sub-10-ms scenarios. Tag minimums first
(`label_heavy` >= 5 labels in the answer, `long_path` >= 1000 bp of paths or a pattern extended beyond k,
`long_query` >= 900 k-mers, `large_answer` >= 64 KiB, `heavy` the top tenth by cost), then one query per bin; every
pick prefers an area of the mode not yet chosen, then a family not yet chosen, then a non-empty result (at most 2
empty results per mode). "Label-heavy" means at least 5 labels in the answer (70 at most in the panel; the v4 pool
tops out at 70 labels for traverse and 18 for resolve and pattern): as heavy as the suite's scenarios go, not a label
stress test.

| mode | areas | tags |
|---|---|---|
| resolve 16 | R 8, N 8 | 3 label-heavy (5-14 labels), 3 long queries (970-1970 k-mers), 6 heavy (0.08-5.7 s), 1 empty |
| traverse 16 | T 14, N 2 | 3 label-heavy (6-70 labels), 5 long-path (1-14 kbp), 1 large answer (the capped exception), 4 heavy (0.4-2.9 s) |
| pattern 16 | P 3, L 2, X 3, M 3, Q 2, N 3 | 3 label-heavy (5-18 labels), 4 long-path, 3 large answers (85-743 KB), 6 heavy (0.1-9.4 s) |

27 families. Weights are equal (1/48) and fixed; the analysis never re-weights silently.

| mode | query | tags | recorded cost ms | result (direct-local) |
|---|---|---|---|---|
| resolve | R-multi-0003 | - | 1.9 | labels 3, num_kmers 70, runs 2 |
| resolve | R-rc-0002 | empty_result | 2.2 | labels 1, num_kmers 70, runs 0 |
| resolve | R-discover-0002 | label_heavy | 2.4 | labels 6, num_kmers 70, runs 7 |
| resolve | N-usecase-0010 | label_heavy | 3.7 | labels 5, num_kmers 120, runs 5 |
| resolve | R-junction-0004 | - | 4.9 | labels 2, num_kmers 170, runs 2 |
| resolve | N-usecase-0012 | label_heavy | 5.4 | labels 14, num_kmers 120, runs 15 |
| resolve | N-realrep-0015 | - | 7.2 | labels 1, num_kmers 30, runs 1 |
| resolve | R-pos-0024 | long_query | 12.4 | labels 1, num_kmers 970, runs 1 |
| resolve | N-nearmiss-0034 | - | 20.5 | labels 1, num_kmers 170, runs 2 |
| resolve | R-edge-0007 | - | 38.5 | labels 1, num_kmers 70, runs 1 |
| resolve | R-pos-0025 | heavy, long_query | 75.4 | labels 1, num_kmers 1970, runs 1 |
| resolve | N-realrep-0019 | heavy | 149.5 | labels 1, num_kmers 30, runs 1 |
| resolve | N-realrep-0009 | heavy | 405.1 | labels 1, num_kmers 30, runs 1 |
| resolve | N-realrep-0021 | heavy | 1314.2 | labels 1, num_kmers 30, runs 1 |
| resolve | N-realrep-0023 | heavy | 3052.3 | labels 1, num_kmers 30, runs 1 |
| resolve | R-pos-0028 | heavy, long_query | 5651.5 | labels 1, num_kmers 1970, runs 1 |
| traverse | T-walk-0025 | - | 4.1 | paths 1, path_bp 100, labels 1 |
| traverse | T-exhaustive-0001 | - | 7.6 | paths 2, path_bp 200, labels 1 |
| traverse | T-labels-0007 | - | 10.8 | paths 2, path_bp 400, labels 1 |
| traverse | T-walk-0015 | - | 17.7 | paths 2, path_bp 200, labels 1 |
| traverse | T-walk-0008 | long_path | 26.7 | paths 1, path_bp 1000, labels 1 |
| traverse | T-derived-0006 | label_heavy, long_path | 41.4 | paths 2, path_bp 2000, labels 13 |
| traverse | T-derived-0002 | label_heavy, long_path | 56.2 | paths 17, path_bp 14198, labels 6 |
| traverse | T-derived-0003 | capped, label_heavy, large_answer, long_path | 68.3 | paths 7, path_bp 7000, labels 70 (cap 64 per node) |
| traverse | T-exhaustive-0003 | - | 101.6 | paths 2, path_bp 200, labels 1 |
| traverse | T-walk-0009 | - | 142.2 | paths 2, path_bp 200, labels 1 |
| traverse | N-realrep-0006 | - | 185.4 | paths 2, path_bp 400, labels 1 |
| traverse | N-realrep-0003 | - | 252.6 | paths 2, path_bp 200, labels 1 |
| traverse | T-graphlet-0001 | heavy | 440.7 | graphlet_lines 30, labels 1 |
| traverse | T-repeat-0001 | heavy | 788.3 | paths 2, path_bp 543, labels 1 |
| traverse | T-walk-0028 | heavy, long_path | 1609.9 | paths 1, path_bp 1000, labels 1 |
| traverse | T-repeat-0003 | heavy | 2918.1 | paths 5, path_bp 699, labels 1 |
| pattern | N-realrep-0007 | - | 1.3 | returned 2, contexts 2 |
| pattern | N-hard-0016 | long_path | 3.2 | returned 2 |
| pattern | L-verified-0001 | label_heavy, long_path | 3.2 | returned 1, labels 8 |
| pattern | P-placed-0011 | label_heavy | 4.7 | returned 7, contexts 7, labels 5 |
| pattern | Q-pred-0002 | - | 7.4 | returned 21, contexts 28 |
| pattern | P-heavy-0002 | - | 16.3 | contexts 1011 (count) |
| pattern | L-verified-0003 | label_heavy, long_path | 25.2 | returned 2, labels 18 |
| pattern | X-peptide-long-0024 | long_path | 27.5 | returned 1, labels 1 |
| pattern | M-iupac-0004 | large_answer | 28.3 | returned 381, contexts 381 |
| pattern | X-edge-0002 | - | 44.2 | returned 20, contexts 20 |
| pattern | X-peptide10-0034 | heavy, large_answer | 102.3 | returned 1554, contexts 1554 |
| pattern | M-suffix-0001 | heavy, large_answer | 142.2 | returned 3479, contexts 3479 |
| pattern | P-placed-0003 | heavy | 245.4 | returned 12, contexts 57, labels 1 |
| pattern | N-nearmiss-0026 | heavy | 414.0 | returned 2, contexts 2 |
| pattern | M-planted-0005 | heavy | 1212.8 | returned 2, contexts 2 |
| pattern | Q-strands-0003 | heavy | 9379.1 | returned 24, contexts 883, labels 1 |

**Rebuilding for v5**: `make_panel.py --baseline <a merged, compared v5 run> [--reevaluated] --include-from
panel_v4.json --name v5 --out panel_v5.json` keeps the v4 queries that are still eligible and fills the rest by the
same rules (`selection.include_not_eligible` lists the ones that dropped out). A new panel is a new workload: compare
releases only on runs of one panel (`analyze.py` refuses two panels).

## The run (`run_bench.py`)

Per release: `48 x (1 + 6 x 3) = 912` requests (`bench_lib.schedule`, the seed recorded in `schedule.json`):

- **block 0**, the warm-up: each query once through `--warmup-arm` (direct-local: the backend), in a seeded order;
  unmeasured. Repeat 0 of the old runs is not "physically cold"; this protocol does not claim cold numbers at all.
- **blocks 1-6**, the measured cycles: every query once per cycle in a fresh seeded order, its three arms back to
  back. The arm order comes from a seeded Latin square: each query runs each of the six orders of (direct-local,
  direct-remote, mcp) exactly once, each cycle runs each order 8 times, each arm comes first in 16 triples per cycle.
  The first arm of a triple pays any cache reload; the balance cancels that in the paired comparisons, and the
  analysis shows each arm's median by position.
- **one request at a time**: each step is one request through the suite's own runners, unchanged. The run's
  manifest (`bench_lib.run_manifest`) gives each query one scenario per phase (`QUERY~w0`, `QUERY~c1..c6`, `repeat`
  1), so `run_direct.py --only QUERY~c3 --warm-every 0` (direct-local on mex over the invocation's ssh master,
  direct-remote from here) and the MCP template with `{ids} QUERY~c3 {repeat} 0` each execute exactly one request.
  The MCP template is run_suite.sh's (`SUITE_MCP_CMD` or `--mcp-cmd` override it; a test keeps the default equal),
  with `{checks_dir}` the run's snapshot, as run_suite.sh passes it. A cancelled MCP job is waited for before the
  next step (`suite_common.mcp_step_result`). `--pause` (1 s) separates steps.
- **identity**: after each block one direct-local GET /capabilities (`run_direct.py --finish`) records the host's
  release and server instance (in `steps.jsonl`); the analysis reads a change as a deployment inside the run. The
  MCP runner records the service's capabilities per invocation (`mcp_runs.jsonl`).

**Isolation**: the suite lock `~/suite-dev/LOCK` on mex (run_suite.sh's protocol and document, mode `bench`), and a
refusal while any run_direct.py runs there: a suite run and a bench run never overlap. As in run_suite.sh (review N3),
a mex that does not answer the runner check never reads as "no runner": the bench refuses to start ("no answer from
..."), and within a run it stops (reason `mex_unreachable`, beside `runner_still_running`) with no closing request; a
lock it cannot release or update, and a stale lock it cannot remove, are said (SD-037). The suite's sampler runs on
mex for the whole session. The **preflight** samples 20 s and refuses to start (exit 4) unless the load guard finds
the host clean; before every block the last minute's samples are judged and a loaded host is waited for (up to 30
min, then exit 4). At each block's start and every 12 steps within a block the run stops (exit 4, resumable) when
the block turned loaded, when the sampler process is gone, or when two checks in a row found no reading of the
load: such a block can give no verdict any more (`--continue-loaded` runs on for descriptive numbers). The MCP
template's `--check-queue` pauses while the search service's queue holds foreign work (it ignores the Logan chunks,
which the load guard covers).

**Recorded**: `bench.json` (panel, seed, arms, options, the MCP template, the guard's thresholds, sessions with
mex's clock offset), `panel.json`, `manifest.json`, `schedule.json`, `steps.jsonl` (each step's times, exit, outcome
and wall, each identity probe), `scripts/` (the snapshot the arms run and judge with), `arms/<arm>/` (each runner's
own files and raw answers), `samples/`. Ctrl-C stops after the step in flight (a second one interrupts it);
`--resume RUN_DIR` continues the recorded schedule with the run's own options and snapshot.

## The load guard (`load_guard.py`)

Each block's window is the union of its executions' intervals, padded 30 s; mex's clock offset aligns samples and
direct-local executions with the orchestrator's clock. A block is **clean** when every sample in its window is at
or below both thresholds and the samples cover its executions without a gap over 15 s; **loaded** when one sample
exceeds a threshold; **unguarded** when the load is unknown (no samples, a gap, a sample without a reading, the
Logan container not reported unless `--allow-missing-sra` declares it stopped, a non-finite value: the sampler
records a NaN or inf docker figure as none, and every reader of the figures, this guard, the panel's load note, the
suite report's load caveat and compare's resource means, reads one as no reading through suite_common.finite_or_none;
SD-018). Defaults:
`sra_all` CPU <= 100 % (one core; the Logan load ran at 5,000-14,000 %) and load1 <= 16. **Neither default is
calibrated**: no quiet sample of mex exists yet, so neither the Logan container's idle CPU nor the host's idle load1
is known. The preflight prints what it saw; set `--sra-max` and `--load1-max` from it before the first release run
(the values are recorded in `bench.json`, every verdict states its thresholds, and the analysis accepts only an
override that tightens them).

**Not guarded: other people's searches on the staging container.** The guard reads the Logan container's CPU and the
host's load1 only. A foreign request to the staging container itself (a direct client, another service) runs beside
the bench unseen, unless it lifts load1 over its limit; the MCP template's `--check-queue` sees only the search
service's queue. Reserve the window with the people who use the staging host; a later guard could judge the staging
container's own `cpu_pct` in the gaps between the bench's steps.

## The analysis (`analyze.py`)

Every measured execution is classed (`success`; `missing`, `timeout`, `transport`, `refusal`, `not_pass`,
`time_limited`, `incomplete`, `unverified`, `result_changed`); latency statistics use successes only, the others
are counted and listed apart with their walls. A success must state every result size the panel recorded for its
arm, and state it alike: another size is another workload (`result_changed`), an answer not stating one cannot be
shown to be the same workload (`unverified`: a moved or wrapped body is never a success).

- per query and arm: the six walls, median, least, greatest; per query and cycle the paired differences and ratios
  mcp/direct-local, direct-remote/direct-local, mcp/direct-remote, and per query their medians;
- aggregates with the panel's fixed weights, overall and per mode: the weighted mean of the per-query median
  differences, the weighted geometric mean of the per-query median ratios, and per arm a latency index (the
  weighted geometric mean of the per-query medians);
- 95 % intervals from a seeded bootstrap (2000 resamples) that resamples queries within each mode, then each
  resampled query's cycles;
- no tail percentile (p90, p99) for any query or arm: six repeats support median comparisons, not tail claims.

**The performance verdict** is emitted only when the load guard finds every block clean with the run's recorded
thresholds (an override may only tighten them), every scheduled step is recorded, the deployment's identity is
stated and unchanged (every block's identity probe states the host's release and server instance, direct-local
states its start identity, every MCP invocation states the service's capabilities; the release, server instance,
index fingerprint and service capabilities did not change), the run's panel matches its identity, and every query
has at least 5 of 6 successes on every arm (`--min-success`). Otherwise it is refused with its reasons, and the
numbers are descriptive of an uncontrolled run. Exit 0 emitted, 3 refused.

**Release A vs B**: `analyze.py RUN_B --baseline RUN_A`: both runs must emit a verdict, on the same panel (whose
identity binds the recorded result sizes), cycles and arms, with the same harness (pause, warm-up arm, MCP template,
base URLs, checks mode, the snapshot's run_direct.py, suite_common.py and sampler.py), and their successes must
state the same result sizes for every query and arm; comparison.json records both harnesses and any threshold
override. Per arm the per-query ratio of medians (B / A), aggregated with the panel's weights, with a bootstrap interval
resampling queries within modes and each run's cycles independently. Per arm `regression` when the interval's
lower end exceeds 1 + tolerance (10 %), `improvement` when its upper end is below 1 - tolerance, `no_regression`
when its upper end is at most 1 + tolerance, else `inconclusive`. Exit 0 pass, 1 regression, 4 inconclusive, 3
refused. To limit time drift, run the two releases in alternating blocks or an ABBA schedule (A, B, B, A panels)
and compare the matched pairs; record the warm-up after each deployment transition.

## Costing

| item | value |
|---|---|
| requests per release | 912 = 48 warm-ups + 864 measured (48 x 6 cycles x 3 arms) |
| other requests | per direct arm 3 capabilities GETs at its first step and one closing GET; one direct-local GET /capabilities per block (7); per MCP step the runner's service capabilities read and queue check (service side); the release GET |
| execution estimate | 31 min of request walls + 23 min of step overhead = 0.90 h (`panel_v4.json` `costing`) |
| basis | each arm's median recorded wall per query in v4 (a loaded host's: an upper estimate) and v4's median gaps between steps: 1.37 s before a direct-local step, 1.34 s before direct-remote, 1.83 s before an MCP step (1 s pause included) |
| allowance per release | **1.5-2 host-hours**: the estimate plus 0.5-1 h of preflight, setup, quiet waits and transitions |
| two releases compared | 3-4 host-hours (two panels, alternating blocks) |

The MCP arm dominates: a pattern search through the service waits for its status poll (~3.5 s per execution in
v4), so the 16 pattern queries cost ~6 min of MCP walls alone. Six queries are of the slow class (up to 600 s
client timeout, 125 s through the front); a timed-out step costs its client limit, not its success latency: the
worst case is the sum of the steps' client limits, which the allowance does not cover. Re-estimate after one
isolated pilot run (`--max-steps` caps a session; `--dry-run` prints the schedule without any request).

## How to run

```sh
cd metagraph/scripts/suite
python3 bench/run_bench.py --panel bench/panel_v4.json --dry-run --seed 1          # the schedule, no request
METAGRAPH_ADMIN_KEY_FILE=/Users/raetsch/tmp/mgss/secrets/staging_admin_key \
  python3 bench/run_bench.py --panel bench/panel_v4.json --load1-max <from the preflight>
python3 bench/run_bench.py --resume bench/runs/<RUN>                                # after Ctrl-C / exit 3 / 4
python3 bench/analyze.py bench/runs/<RUN>                                           # again, offline
python3 bench/analyze.py bench/runs/<RUN_B> --baseline bench/runs/<RUN_A>           # release A vs B
python3 bench/load_guard.py bench/runs/<RUN>                                        # the guard alone
```

Only when mex is free of the Logan load (the preflight refuses otherwise) and no suite run holds the lock.

## Validation and limits

- Offline only: the tests drive `run_bench.py` end to end against the fake host, a fake MCP runner writing raw
  answers, a fake `docker` (clean, loaded, turning loaded mid-run, down), mex in a local directory: the schedule
  executed in order with no overlapping requests, the request counts, the identity probes, the snapshot passed as
  `{checks_dir}`, the preflight refusal, the stop on a block turning loaded, on a dead sampler and on a sampler
  without readings, `--max-steps` and `--resume` without a repeated execution, Ctrl-C after the step in flight, the
  lock, the template refusals. The analysis is tested on synthetic recorded runs (exact ratios, intervals, every
  refusal reason, regression, improvement, inconclusive). The verification of the lane (VERIFY-BENCH, D1-D8) added
  a test per defect: an answer stating no recorded size, a panel identity without the sizes, differing successful
  sizes, loosened thresholds, an unknown or changed identity (index_fp too), the label-list cap and pattern
  truncation fields, a non-finite reading, a dying sampler, a differing harness.
- Never run against staging or the real MCP runner. To check in the pilot: the real runner's per-invocation
  overhead (it reads the service's capabilities and the queue before each execution), whether its measured wall
  includes a connection set-up (each MCP step is a fresh runner process; each direct step a fresh connection, whose
  set-up run_direct.py records apart as `connect_ms`, outside `wall_ms`), the real runner's treatment of `{repeat} 0`
  (every bench MCP step uses it: the runner's "cold" index, which should change only labels such as the idempotency
  key's), the sra_all and load1 thresholds, and the costing (re-estimate from the pilot's walls and step gaps).
- The panel's recorded walls and costs come from a loaded run; they rank and cost, they are no baseline. The first
  clean run of a release is its baseline.
- Budget stress (deliberately capped queries, stop reasons, upper bounds) is a separate protocol (GPT's item 5), not
  part of this panel: capped answers are not eligible here, except T-derived-0003's deterministic default label cap,
  a stated exception (above).
