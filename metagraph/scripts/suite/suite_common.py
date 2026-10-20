#!/usr/bin/env python3
"""Shared helpers of the suite's direct runner, sampler, comparator and report.

Stdlib only: this module is copied to mex and imported by run_direct.py there (python3.12).
Everything here is pure except the small file helpers at the end, so the unit tests cover it offline.

Contract: metagraph/docs/DESIGN-test-suite.md (v2), section 8.
"""

import datetime
import gzip
import hashlib
import json
import math
import os
import re
import shlex
import sys
import zlib

ARMS = ("direct-local", "direct-remote", "mcp")
DIRECT_ARMS = ("direct-local", "direct-remote")
OUTCOMES = ("pass", "fail", "finding", "inconclusive", "not_implemented_mcp", "service_by_design", "error")
# The outcomes a check may report (DESIGN §8.2 plus the review of lane S1): `inconclusive` = the scenario could
# not test its property (an allowed 503, a withheld answer hiding the source, a stated resource stop, ...);
# `finding` = the host broke a stated bound or SPEC promise without failing the request. Anything else is read
# as a failure.
CHECK_OUTCOMES = ("pass", "fail", "not_exposed", "inconclusive", "finding")
# async's MCP runner (suite_mcp_runner.py HOST_STATE_CODES) and the design lead's rule of 2026-10-13: an MCP answer
# refused with one of these codes states an environmental host state, not a result: inconclusive (unless the scenario
# expects the code, or the plausibility rule of mcp_host_state rejects it), whatever the checks say; resolve_timeout:
# the service's resolve client gave up before the host's overrun answer under load (G23 as VERIFY-RP amends it);
# traversal_resolve_deadline: the host's /resolve 503 deadline the service relays (G22, VERIFY-2 item 3: the request's
# budget under load, SPEC-labeled-traversal-core §4.5, not a result)
MCP_HOST_STATE_CODES = ("uninterruptible_overrun_exceeds_allowance", "resolve_timeout", "traversal_resolve_deadline")
TIMING_CLASSES = ("fast", "slow", "long")

# Keys left out of the answer digest. `timing` varies by run (DESIGN §8.3); `headers` is transport metadata
# (the front proxy may change it) and not part of the answer.
DIGEST_DROP_KEYS = frozenset(["timing"])
# On the capabilities documents these prefixes name values that change with the server's life (the delivery
# calibration, the longest uninterruptible chunk seen so far), not with the document's contract.
CAPABILITY_VOLATILE_PREFIXES = ("measured_", "observed_")


def outcome_known(outcome):
    """Whether an execution line's outcome is in the vocabulary every arm and the comparator share (OUTCOMES). One
    outside it (a runner's mark not agreed in this vocabulary first) is counted apart and read as a failure by the
    comparator, the re-evaluation and the report, never as a pass."""
    return isinstance(outcome, str) and outcome in OUTCOMES


def outcome_value(outcome):
    """An outcome as a comparable and hashable value: one stored as a JSON list or object (no runner writes one) as
    its canonical JSON text, which outcome_known reads as outside the vocabulary; any other value as recorded."""
    if isinstance(outcome, (list, dict)):
        return json.dumps(outcome, sort_keys=True)
    return outcome


def recorded_outcome(line):
    """The outcome an execution line recorded, as the comparator, the report and the gate state it: merge_run reads
    an outcome outside the vocabulary as `error` and keeps the recorded one under `outcome_runner`
    (merge_run.unknown_outcomes), so for such a line (outcome error, outcome_runner present and not outcome_known)
    the recorded one, else the line's outcome; a list or object as its JSON text (outcome_value). `not
    outcome_known(recorded_outcome(line))` is the one test of an outcome outside the vocabulary, whether or not the
    set went through merge_run."""
    outcome = line.get("outcome")
    if outcome == "error" and "outcome_runner" in line and not outcome_known(line.get("outcome_runner")):
        outcome = line.get("outcome_runner")
    return outcome_value(outcome)


# ---------------------------------------------------------------------------------------------------------------
# time

def iso_utc(ts):
    """A POSIX time as an ISO-8601 UTC string with milliseconds, e.g. 2026-10-13T10:42:36.123Z."""
    dt = datetime.datetime.fromtimestamp(ts, tz=datetime.timezone.utc)
    return dt.strftime("%Y-%m-%dT%H:%M:%S.") + "%03dZ" % (dt.microsecond // 1000)


def parse_utc(text):
    """The inverse of iso_utc (also accepts no fraction and a +00:00 suffix); None for anything else."""
    if not isinstance(text, str):
        return None
    t = text.strip()
    if t.endswith("Z"):
        t = t[:-1] + "+00:00"
    try:
        dt = datetime.datetime.fromisoformat(t)
    except ValueError:
        return None
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=datetime.timezone.utc)
    return dt.timestamp()


def run_stamp(ts):
    """The UTC part of a run directory name: 20261013T104236Z."""
    return datetime.datetime.fromtimestamp(ts, tz=datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def release_short(release):
    """refseq97-b7214d29a934ff... -> refseq97-b7214d29; a name safe for a directory either way."""
    if not release:
        return "unknown"
    head, sep, tail = str(release).partition("-")
    short = head + (sep + tail[:8] if sep else "")
    return re.sub(r"[^A-Za-z0-9._-]", "_", short)


# ---------------------------------------------------------------------------------------------------------------
# manifest

def canonical_json(obj):
    """The canonical JSON of the manifest's content hash: sorted keys, no whitespace, UTF-8."""
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")


def scenarios_sha256(scenarios):
    return hashlib.sha256(canonical_json(scenarios)).hexdigest()


def relations_sha256(relations):
    """The sha256 of the canonical JSON of a manifest's `relations` (build_manifest's relations_sha256)."""
    return hashlib.sha256(canonical_json(relations if relations is not None else [])).hexdigest()


def manifest_identity(manifest):
    """The identity of a manifest's test plan, computed from its content: content_sha256 of `scenarios` and
    relations_sha256 of `relations` (of [] when it states none). The relations are evaluated across the scenarios'
    answers, so two manifests with the same scenarios but other relations test other things: both hashes must be
    equal for two arms, a resume or a gate to speak of the same manifest."""
    manifest = manifest if isinstance(manifest, dict) else {}
    return {"content_sha256": scenarios_sha256(manifest.get("scenarios", [])),
            "relations_sha256": relations_sha256(manifest.get("relations", []))}


def manifest_problems(manifest):
    """The stated hashes of a manifest that its content does not match, as sentences ([] when it verifies):
    content_sha256 always, relations_sha256 where the manifest states one."""
    ident = manifest_identity(manifest)
    out = []
    if ident["content_sha256"] != manifest.get("content_sha256"):
        out.append("content_sha256 %s does not match the canonical JSON of its scenarios (%s)" % (
            manifest.get("content_sha256"), ident["content_sha256"]))
    if "relations_sha256" in manifest and ident["relations_sha256"] != manifest.get("relations_sha256"):
        out.append("relations_sha256 %s does not match the canonical JSON of its relations (%s)" % (
            manifest.get("relations_sha256"), ident["relations_sha256"]))
    return out


def verify_manifest(manifest):
    """(ok, computed): whether content_sha256 matches the canonical JSON of `scenarios` and, where the manifest states
    relations_sha256, that one the canonical JSON of `relations` (manifest_problems); computed is the scenarios'
    sha256, the identity the arms' run.json files record as computed_sha256."""
    return not manifest_problems(manifest), scenarios_sha256(manifest.get("scenarios", []))


def parse_route(route):
    """'POST /pattern' -> ('POST', '/pattern'); a bare path is a GET."""
    parts = str(route).strip().split(None, 1)
    if len(parts) == 1:
        return "GET", parts[0]
    return parts[0].upper(), parts[1].strip()


def select_scenarios(manifest, areas=None, limit=None, ids=None):
    """The scenarios to run, in manifest order, as (global index, scenario).

    `areas` and `ids` filter (None: all); `limit` keeps the first N after filtering. The global index is the
    scenario's position in the whole manifest, so the "every tenth" warm repeat picks the same scenarios
    whatever subset runs. A warm-up (a scenario whose `warms` names another, lane V5) travels with its target: a
    selection that takes the target takes its warm-up too (the manifest keeps it right before the target), so that
    the target's first execution meets the rows it read, whatever subset runs.
    """
    area_set = set(a.strip() for a in areas) if areas else None
    id_set = set(i.strip() if isinstance(i, str) else i for i in ids) if ids else None   # 'A, B' as 'A,B' (SD-369)
    scenarios = manifest.get("scenarios", [])
    out = []
    for index, scenario in enumerate(scenarios):
        if area_set is not None and scenario.get("area") not in area_set:
            continue
        if id_set is not None and scenario.get("id") not in id_set:
            continue
        out.append((index, scenario))
        if limit is not None and limit > 0 and len(out) >= limit:
            break
    chosen = set(s.get("id") for _, s in out)
    pulled = [(i, s) for i, s in enumerate(scenarios) if s.get("warms") in chosen and s.get("id") not in chosen]
    if pulled:
        out = sorted(out + pulled, key=lambda x: x[0])
    return out


def repeat_count(index, scenario, warm_every):
    """Executions of a scenario per arm: its `repeat`, and one warm repeat for every `warm_every`-th scenario."""
    count = scenario.get("repeat", 1)
    if not isinstance(count, int) or isinstance(count, bool) or count < 1:
        count = 1
    if count == 1 and warm_every and warm_every > 0 and index % warm_every == 0:
        count = 2
    return count


def plan_executions(selected, arm, warm_every=10):
    """(planned, skipped): planned is [(index, scenario, repeat, cache)] in run order; skipped is
    [{scenario, reason}] for scenarios this arm does not run (the long class behind the front proxy)."""
    planned, skipped = [], []
    for index, scenario in selected:
        if arm == "direct-remote" and scenario.get("timing_class") == "long":
            skipped.append({"scenario": scenario.get("id"),
                            "reason": "timing_class long: over the front proxy's ~120 s cut, direct-local only"})
            continue
        for repeat in range(repeat_count(index, scenario, warm_every)):
            planned.append((index, scenario, repeat, "cold" if repeat == 0 else "warm"))
    return planned, skipped


def execution_key(line):
    return (line.get("scenario"), line.get("arm"), line.get("repeat"))


# A warm-up run again for a resume (the resume rule below): a line with the warm-up's own scenario id and repeat,
# `rerun_for_resume` true and `rerun_target` the target it ran before; the earlier execution of the same key stays.
RERUN_FIELDS = ("rerun_for_resume", "rerun_target")


def is_rerun(line):
    """Whether an execution line is a warm-up's re-run for a resume: timing evidence next to its target, never a
    judged execution of its key (merge_run, reevaluate, compare, report and the gate judge the original)."""
    return isinstance(line, dict) and line.get("rerun_for_resume") is True


def judged_lines(rows):
    """The execution lines that are executions of their key: every line but a warm-up's re-run for a resume."""
    return [r for r in rows if not is_rerun(r)]


def rerun_raw_name(scenario_id, arm, repeat, started):
    """The raw file of a warm-up's re-run: its own, beside the original execution's (raw_name), one per re-run."""
    safe = re.sub(r"[^A-Za-z0-9._-]", "_", str(scenario_id))
    return "raw/%s.%s.%d.rerun-%d.json.gz" % (safe, arm, repeat, int(round(started * 1000)))


def resume_rewarms(scope, done, target_left=True):
    """The resume rule of a runner's invocation: a warm-up (a scenario whose `warms` names its target) in the
    invocation's scope with every execution there recorded, whose target in the scope has an execution left, runs
    again right before the target's first execution left, so that execution meets the rows its warm-up read even
    after an interruption between the two. scope: [(index, scenario, repeat, cache)] (plan_executions), done: the
    recorded (scenario, arm, repeat) keys as (scenario, repeat) pairs. Returns {target id: [(index, warm-up)]} in
    scope order; a warm-up with an execution left runs as planned (it is right before its target) and is not
    listed. |target_left| False (run_direct.py --rewarm-only: the step before the target's first MCP call, which
    paired_plan's before_mcp decided) does not ask the target for an execution left on this arm."""
    by_sid, order = {}, []
    for index, scenario, repeat, _ in scope:
        sid = scenario.get("id")
        if sid not in by_sid:
            by_sid[sid] = (index, scenario, [])
            order.append(sid)
        by_sid[sid][2].append(repeat)
    out = {}
    for sid in order:
        index, scenario, repeats = by_sid[sid]
        target = scenario.get("warms")
        if not target or target not in by_sid:
            continue
        if any((sid, r) not in done for r in repeats):
            continue
        if not target_left or any((target, r) not in done for r in by_sid[target][2]):
            out.setdefault(target, []).append((index, scenario))
    return out


def feeds_previous(line):
    """Whether a recorded execution is one of the earlier repeats that checks.evaluate's `previous` holds: the host
    answered it (an HTTP status) and its runner did not end it in `error`. An execution without a status (a transport
    failure, a client timeout, an MCP call the runner gave up on) or ended in `error` (an undecodable body, a runner
    failure) has the digest of its stand-in answer, not of the host's: compared, it would fail the next repeat's
    `deterministic` although the host answered nothing. Both runners (run_direct.py, async's suite_mcp_runner.py) and
    reevaluate.py apply this rule."""
    http = line.get("http")
    return isinstance(http, int) and not isinstance(http, bool) and http > 0 and line.get("outcome") != "error"


def beyond_status(scenarios):
    """The expectation names of the scenarios that the checks stub cannot evaluate (anything but `status`)."""
    names = set()
    for scenario in scenarios:
        for exp in scenario.get("expect") or []:
            if isinstance(exp, dict) and exp.get("check") != "status":
                names.add(str(exp.get("check")))
    return sorted(names)


def single_mcp_execution(scenario):
    """Whether a scenario runs once on MCP (manifest v8, async's runner questions, 2026-10-20): an area G scenario
    (a graphlet block: the service's 300 s result cache answers a repeat from the first, so MCP repeats prove nothing,
    DESIGN-v8 §3) and a stop scenario (mcp.flow kind "stop": each repeat would stop another job). No scenario of
    v1-v7 has either."""
    if isinstance(scenario.get("graphlet"), dict):
        return True
    flow = (scenario.get("mcp") or {}).get("flow")
    return isinstance(flow, dict) and flow.get("kind") == "stop"


def mcp_repeats(index, scenario, warm_every=10):
    """The MCP executions of a scenario in paired mode, as repeat indices: repeat 0 before the direct arms
    (cold), the others after them (warm); one more than the direct arms' repeats when those run once, so
    every scenario has an MCP warm execution. A scenario the manifest marks not implemented via MCP gets
    repeat 0 only (the MCP runner records the mark), so do area G and the stop scenarios (single_mcp_execution);
    the long class none (the front cuts every arm behind it)."""
    if scenario.get("timing_class") == "long":
        return []
    if ((scenario.get("mcp") or {}).get("not_implemented")) or single_mcp_execution(scenario):
        return [0]
    return list(range(max(2, repeat_count(index, scenario, warm_every))))


def paired_plan(manifest, run_dir, areas=None, limit=None, warm_every=10, arms=("direct-local", "direct-remote"),
                ids=None):
    """The paired loop's work per scenario of the selection (select_scenarios: areas, limit, ids), in manifest order:
    [{scenario, direct-local, direct-remote, mcp_cold, mcp_warm, rewarm}] where a direct arm's value is True
    (executions left), False (all recorded) or None (the arm does not run it), mcp_cold likewise, and mcp_warm the
    warm repeat indices left. Recorded executions are read from RUN_DIR/arms/<arm>/executions.jsonl.

    rewarm {arm: [warm-up id]}: the resume rule (resume_rewarms) per arm, decided before anything of this invocation
    runs. A target with executions left on an arm whose warm-up has all of its executions there recorded (an
    interruption between the two) gets its warm-up again right before it on that arm: the paired loop names the
    warm-up with the target in the arm's step (run_direct.py --only W,T; the MCP runner --ids W,T), and the runner,
    whose scope then holds both, runs the warm-up again first (run_direct.py: repeat 0; the MCP runner: the step's
    repeat). On the MCP arm that is the target's first step left (cold, else its first warm repeat). A warm-up the
    manifest marks not implemented via MCP sends nothing there and is never run again on that arm.

    rewarm {"before_mcp": [warm-up id]}: a warm-up the MCP arm does not send (the long class, v6's P-warm: it runs on
    direct-local only and warms the host for its target's first MCP call; or one marked not implemented via MCP)
    whose direct-local executions are all recorded, while its target has an MCP step left, runs again on
    direct-local right before the target's first MCP step left (run_direct.py --rewarm-only --only W,T: the warm-up
    alone), as the paired loop's order has it without an interruption; the direct-local rule above then leaves it
    out where the target's MCP cold call comes first (that call follows the re-run already). Where the target's
    direct-local step comes before its first MCP step left (a warm repeat), that step's re-run serves (VERIFY of
    round 2, L1: the re-run inside the target's direct-local step came after its MCP cold call, which met a cold
    host again)."""
    recorded = {}
    for arm in ARMS:
        rows, _ = read_jsonl(os.path.join(run_dir, "arms", arm, "executions.jsonl"))
        recorded[arm] = set(execution_key(r) for r in judged_lines(rows) if r.get("arm") == arm)
    out = []
    by_sid = {}
    for index, scenario in select_scenarios(manifest, areas, limit, ids):
        sid = scenario.get("id")
        item = {"scenario": sid, "mcp_cold": None, "mcp_warm": [], "rewarm": {}}
        by_sid[sid] = (item, scenario)
        for arm in DIRECT_ARMS:
            if arm not in arms:
                item[arm] = None
                continue
            planned, _ = plan_executions([(index, scenario)], arm, warm_every)
            item[arm] = (any((sid, arm, p[2]) not in recorded[arm] for p in planned)) if planned else None
        if "mcp" in arms:
            reps = mcp_repeats(index, scenario, warm_every)
            if reps:
                item["mcp_cold"] = (sid, "mcp", 0) not in recorded["mcp"]
                item["mcp_warm"] = [r for r in reps[1:] if (sid, "mcp", r) not in recorded["mcp"]]
        out.append(item)
    for item in out:
        target = by_sid[item["scenario"]][1].get("warms")
        if not target or target not in by_sid:
            continue
        titem = by_sid[target][0]
        warm_up = by_sid[item["scenario"]][1]
        for arm in DIRECT_ARMS:
            if titem.get(arm) is True and item.get(arm) is False:
                titem["rewarm"].setdefault(arm, []).append(item["scenario"])
        mcp_left = item["mcp_cold"] is True or item["mcp_warm"]
        if (item["mcp_cold"] is False and not mcp_left and not (warm_up.get("mcp") or {}).get("not_implemented")
                and (titem["mcp_cold"] is True or titem["mcp_warm"])):
            titem["rewarm"].setdefault("mcp", []).append(item["scenario"])
        not_sent = item["mcp_cold"] is None or (warm_up.get("mcp") or {}).get("not_implemented")
        if ("mcp" in arms and not_sent and item.get("direct-local") is False
                and (titem["mcp_cold"] is True or titem["mcp_warm"])):
            local = titem["rewarm"].get("direct-local") or []
            if titem["mcp_cold"] is True:
                if item["scenario"] in local:
                    local.remove(item["scenario"])
                    if not local:
                        titem["rewarm"].pop("direct-local")
                titem["rewarm"].setdefault("before_mcp", []).append(item["scenario"])
            elif item["scenario"] not in local:
                titem["rewarm"].setdefault("before_mcp", []).append(item["scenario"])
    return out


# ---------------------------------------------------------------------------------------------------------------
# after an MCP step of the paired loop (async's answers to the verification of lane S2b)

# How long a request of a timing class may take (the direct runner's client timeouts); a request that states a time
# budget above it gets the budget plus CLIENT_BUDGET_MARGIN_S instead (client_timeout_s).
TIMING_CLASS_TIMEOUT_S = {"fast": 120.0, "slow": 600.0, "long": 1800.0}
# The margin a client timeout keeps above a request's stated time budget (v6: the /pattern budget of 110,000 ms
# leaves the fast class's 120 s only 10 s above the host's 503 deadline, answered at ~budget + 250 ms): the host's
# finalisation reserve, the answer's transfer and the clocks. Narrower than the window
# in which checks._budget_503 still reads a late 503 as inconclusive (max(1 s, the budget) past it): a 503 later
# than this margin is a client timeout, an error; the lateness seen so far is at most 0.3 s.
CLIENT_BUDGET_MARGIN_S = 30.0
# The host's time budget of a /pattern request that states none (SPEC-pattern-search: default_time_budget_ms).
PATTERN_DEFAULT_BUDGET_S = 60.0
# Added to a cancelled job's bound: the answer's transfer and the clocks.
CANCEL_MARGIN_S = 2.0
# jobs.jsonl events of async's runner after which the job may still run on the host.
CANCEL_EVENTS = ("canceled", "cancel_failed")


def _budget_values(obj):
    if isinstance(obj, dict):
        for key, value in obj.items():
            if key == "time_budget_ms" and _is_number(value):
                yield value
            else:
                for v in _budget_values(value):
                    yield v
    elif isinstance(obj, list):
        for item in obj:
            for v in _budget_values(item):
                yield v


def stated_request_bound_s(scenario):
    """The time budget the scenario's request states, in seconds (the largest time_budget_ms in the MCP arguments,
    else in the direct request), None when it states none."""
    stated = list(_budget_values((scenario.get("mcp") or {}).get("arguments"))) or \
        list(_budget_values(scenario.get("request")))
    return max(stated) / 1000.0 if stated else None


def client_timeout_s(scenario, timeouts=None):
    """The direct runner's client timeout of the scenario's request: its timing class's (|timeouts|, else
    TIMING_CLASS_TIMEOUT_S), or the request's stated time budget plus CLIENT_BUDGET_MARGIN_S where that is more (a
    budgeted request is answered at its budget at the latest, a /pattern 503 deadline at ~budget + 250 ms: the
    timeout must exceed it with margin, or the deadline, an inconclusive, reads as a client timeout, an error)."""
    table = timeouts or TIMING_CLASS_TIMEOUT_S
    t = table.get(scenario.get("timing_class") or "fast", table["fast"])
    bound = stated_request_bound_s(scenario)
    if bound is not None:
        t = max(t, bound + CLIENT_BUDGET_MARGIN_S)
    return t


def request_bound_s(scenario):
    """How long the host may work on the scenario's request: its time_budget_ms (the largest stated in the MCP
    arguments, else in the direct request; stated_request_bound_s), else the /pattern default budget, else the timing
    class's timeout (a traversal without a budget has no deadline of its own)."""
    stated = stated_request_bound_s(scenario)
    if stated is not None:
        return stated
    if "/pattern" in str(scenario.get("route") or ""):
        return PATTERN_DEFAULT_BUDGET_S
    return TIMING_CLASS_TIMEOUT_S.get(scenario.get("timing_class") or "fast", TIMING_CLASS_TIMEOUT_S["fast"])


def mcp_step_result(out_dir, manifest, ids, repeat, now):
    """What an MCP step (one runner invocation for `ids` and `repeat`) left: (outcomes, wait_s). outcomes: the
    recorded outcome per id, `missing` when the runner wrote none (a runner failure), `unknown` for an outcome outside
    the shared vocabulary (outcome_known; the paired loop counts both like `error`). wait_s: how long to wait
    before the next step, because an execution ended `error` with its job cancelled rather than seen terminal
    (async's runner cancels at its poll deadline; the cancelled call may still run on the host): until the
    cancellation plus the scenario's request bound (request_bound_s) and a margin, unless jobs.jsonl shows the
    job terminal (or the cancel found it already terminal)."""
    rows, _ = read_jsonl(os.path.join(out_dir, "executions.jsonl"))
    jobs, _ = read_jsonl(os.path.join(out_dir, "jobs.jsonl"))
    scen = dict((s_.get("id"), s_) for s_ in (manifest.get("scenarios") or []))
    events = {}
    for j in jobs:
        if j.get("id") is not None and j.get("event"):
            events.setdefault(str(j["id"]), []).append(j)
    outcomes, wait = [], 0.0
    for sid in ids:
        mine = [r for r in rows if r.get("arm") == "mcp" and r.get("scenario") == sid and r.get("repeat") == repeat]
        if not mine:
            outcomes.append("missing")
            continue
        line = mine[-1]
        outcomes.append(str(line.get("outcome")) if outcome_known(line.get("outcome")) else "unknown")
        if line.get("outcome") != "error":
            continue
        job = line.get("job") if isinstance(line.get("job"), dict) else {}
        evs = [e for key in (job.get("search_id"), job.get("job_id")) if key is not None
               for e in events.get(str(key), [])]
        if any(e.get("event") == "terminal" for e in evs):
            continue
        cancels = [e for e in evs if e.get("event") in CANCEL_EVENTS and e.get("detail") != "already_terminal"]
        if not cancels:
            continue
        at = max([parse_utc(e.get("at")) for e in cancels if parse_utc(e.get("at")) is not None] or
                 [parse_utc(line.get("ended_utc")) or now])
        until = at + request_bound_s(scen.get(sid) or {}) + CANCEL_MARGIN_S
        wait = max(wait, until - now)
    return outcomes, max(0.0, wait)


# ---------------------------------------------------------------------------------------------------------------
# the host's cache state per execution (review H1)

HOST_CACHE_WARM_S = 60.0
HOST_CACHE_VALUES = ("cold", "warm", "aged", "not_reached", "unknown")


# The MCP answers that relay a host deadline: the host worked on the request until its deadline (it reached the host
# and its row cache) although the service states neither the host's timing nor a backend phase. The host-state codes
# (the relayed /resolve deadline, the service's resolve client giving up before the host's overrun answer, an
# attempt's overrun) and the relayed /pattern deadline (search_failed, the text pattern_deadline: ...); the checks'
# time_limited names such a relay.
RELAYED_DEADLINE_CODES = MCP_HOST_STATE_CODES
RELAYED_PATTERN_DEADLINE = ("search_failed", "pattern_deadline")
RELAYED_TIME_LIMITED = "relayed by the service"


def backend_timed(line):
    """Whether an MCP execution records a backend call's timing: the host's own (`server.elapsed_ms`) or the
    service's backend phase (`mcp_phases.backend_ms`) is a number."""
    server = line.get("server") if isinstance(line.get("server"), dict) else {}
    phases = line.get("mcp_phases") if isinstance(line.get("mcp_phases"), dict) else {}
    return _is_number(server.get("elapsed_ms")) or _is_number(phases.get("backend_ms"))


def relayed_host_deadline(line):
    """Whether an MCP execution's answer relays a host deadline (RELAYED_DEADLINE_CODES, the relayed /pattern
    deadline, or a time_limited naming a relay): the host worked on it, without a timing on record. A resolve_timeout
    in the write phase (the line's timeout_terms, else its error's timeout_phase) is none: the service may never have
    delivered the request."""
    err = line.get("error") if isinstance(line.get("error"), dict) else {}
    code = str(err.get("code") or "")
    if code == "resolve_timeout":
        terms = line.get("timeout_terms") if isinstance(line.get("timeout_terms"), dict) else {}
        if (terms.get("timeout_phase") or err.get("timeout_phase")) == "write":
            return False
    if code in RELAYED_DEADLINE_CODES:
        return True
    if code == RELAYED_PATTERN_DEADLINE[0] and str(err.get("text") or "").startswith(RELAYED_PATTERN_DEADLINE[1]):
        return True
    tl = line.get("time_limited")
    return isinstance(tl, str) and RELAYED_TIME_LIMITED in tl


def reached_host(line):
    """Whether an execution reached the host, and so touched its row cache: a direct request that got an
    answer or timed out (not a refused connection); an MCP execution only on evidence of a backend call: a
    backend timing (backend_timed) or an answer relaying a host deadline (relayed_host_deadline: the host worked
    until its deadline, 30-60 s). Never from the MCP outcome (review N2): a refusal at submit can read `fail`
    without reaching the host, and an `error` or a service_by_design mark can follow a backend call (a relayed
    refusal does)."""
    if line.get("arm") == "mcp":
        return backend_timed(line) or relayed_host_deadline(line)
    err = line.get("error") or {}
    return not (not line.get("http") and err.get("code") == "transport")


def interrupted_on_host(line):
    """Whether an execution ended at the host's deadline without its answer: a direct 503 deadline, or an MCP answer
    relaying a host deadline (relayed_host_deadline). The next execution's host cache is warm by the recent-repeat
    heuristic, but the host's work was cut short, so its row cache may hold only part of what a complete execution
    leaves there."""
    if line.get("arm") == "mcp":
        return relayed_host_deadline(line)
    err = line.get("error") if isinstance(line.get("error"), dict) else {}
    return line.get("http") == 503 and err.get("code") == "deadline"


def line_times(line, offset_s=0.0):
    """(start, end) POSIX times of an execution on the reference clock (offset_s: how far the clock that
    stamped it is ahead); end from ended_utc, else start + wall_ms."""
    start = parse_utc(line.get("started_utc"))
    end = parse_utc(line.get("ended_utc"))
    wall = line.get("wall_ms")
    if end is None and start is not None and _is_number(wall):
        end = start + wall / 1000.0
    if start is not None:
        start -= offset_s
    if end is not None:
        end -= offset_s
    return start, end


def attach_reruns(rows, reruns, clock_offsets=None):
    """Puts each warm-up's re-run for a resume (is_rerun) next to its target, as timing evidence: the target's first
    execution on the same arm that started after the re-run ended gets it under `warm_up_reruns` (a list: scenario,
    repeat, started_utc, ended_utc, wall_ms, http, outcome, raw, and before_s, from the re-run's end to the target's
    start). Returns the re-runs in time order with `target_execution` [scenario, arm, repeat], or None when the target
    did not run after it on that arm (an interruption right after the re-run). |rows| (the judged lines, in place)
    never gain a line: a re-run is no execution of its key."""
    offsets = clock_offsets or {}
    for r in rows:
        r.pop("warm_up_reruns", None)
    out = []
    for rr in sorted(reruns, key=lambda r: parse_utc(r.get("started_utc")) or 0.0):
        arm = rr.get("arm")
        off = offsets.get(arm, 0.0) or 0.0
        _, end = line_times(rr, off)
        best = None
        for r in rows:
            if r.get("arm") != arm or r.get("scenario") != rr.get("rerun_target"):
                continue
            start, _ = line_times(r, off)
            # the stamps have milliseconds: a start in the re-run's last millisecond is after it
            if start is None or end is None or start < end - 0.001:
                continue
            if best is None or start < best[0]:
                best = (start, r)
        entry = dict((k, rr.get(k)) for k in ("scenario", "repeat", "started_utc", "ended_utc", "wall_ms", "http",
                                              "outcome", "raw"))
        summary = dict(entry, arm=arm, rerun_target=rr.get("rerun_target"), target_execution=None)
        if best is not None:
            best[1].setdefault("warm_up_reruns", []).append(dict(entry, before_s=round(max(0.0, best[0] - end), 3)))
            summary["target_execution"] = list(execution_key(best[1]))
            summary["before_s"] = round(max(0.0, best[0] - end), 3)
        out.append(summary)
    return out


def assign_host_cache(rows, clock_offsets=None, warm_s=HOST_CACHE_WARM_S):
    """Adds `host_cache`, `since_prev_s` and `host_cache_after` to every execution line (in place), from the
    order in which the run's executions of each scenario reached the host, over all arms:
    cold = the first execution of the scenario on the host in this run; warm = another execution of it ended
    at most warm_s seconds before this one started; aged = longer ago (the cache state is unknown);
    not_reached = the execution never reached the host; unknown = no start time. since_prev_s is the time from
    the end of the scenario's previous host execution to this start, host_cache_after that execution's arm, and
    `host_cache_after_interrupted` true (only then present) when that execution ended at a host deadline without
    its answer (interrupted_on_host: a relayed deadline), so a timing can tell a warm cache after a complete
    execution from one after a cut one.
    clock_offsets {arm: seconds} puts the arms on one clock (direct-local is stamped on mex)."""
    offsets = clock_offsets or {}
    arm_rank = {"direct-local": 0, "direct-remote": 1, "mcp": 2}
    timed = []
    for i, r in enumerate(rows):
        start, end = line_times(r, offsets.get(r.get("arm"), 0.0) or 0.0)
        timed.append((start if start is not None else float("inf"), arm_rank.get(r.get("arm"), 9),
                      r.get("repeat") or 0, i, start, end))
    timed.sort(key=lambda t: t[:4])
    last = {}
    for _, _, _, i, start, end in timed:
        r = rows[i]
        sid = r.get("scenario")
        prev = last.get(sid)
        since = None
        if start is None:
            label = "unknown"
        elif prev is None:
            label = "cold"
        else:
            since = round(start - prev[0], 3)
            label = "warm" if since <= warm_s else "aged"
        reached = reached_host(r)
        if not reached:
            label = "not_reached"
        r["host_cache"] = label
        r["since_prev_s"] = since
        r["host_cache_after"] = prev[1] if prev else None
        r.pop("host_cache_after_interrupted", None)
        if prev and prev[2]:
            r["host_cache_after_interrupted"] = True
        if reached and start is not None:
            stop = end if end is not None else start
            last[sid] = (max(stop, prev[0]) if prev else stop, r.get("arm"), interrupted_on_host(r))
    return rows


# the MCP runner's start refusals (async's suite_mcp_runner.py exit codes): its step recorded nothing, and the loop
# stops at once with reason mcp_exit_<rc> (run_suite.sh mcp_after, bench/run_bench.py; SD-029); 130 is an interrupt
MCP_START_REFUSALS = {2: "the MCP runner refused to run", 3: "the MCP runner refused to run",
                      4: "the MCP runner refused its by-design table", 5: "the MCP runner refused the run's snapshot"}


def mex_clock_offsets(orchestrator):
    """{arm: seconds its clock is ahead of the laptop's}: direct-local is stamped on mex, whose offset the
    orchestrator measures per invocation (orchestrator.json sessions' mex_clock_offset_s); the upper median of the
    invocations is used, {} without one. merge_run applies it and records it in run.json host_cache.clock_offsets_s,
    which the bench panel reads (SD-028)."""
    vals = []
    for sess in (orchestrator or {}).get("sessions", []) or []:
        v = sess.get("mex_clock_offset_s")
        try:
            vals.append(float(v))
        except (TypeError, ValueError):
            continue
    if not vals:
        return {}
    vals.sort()
    return {"direct-local": vals[len(vals) // 2]}


def finite_or_none(v):
    """A sampler figure (sra_all_cpu_pct, host_load1, ...) as a reading: a finite number, else None (a NaN compares
    False with every limit and an inf is no figure: no reading, never a clean one; the bench load guard's rule since
    c6306ec1a, which the report, the comparison's resources and the panel's load note read too: SD-018)."""
    return v if isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v) else None


def _number_or_none(v):
    return float(v) if isinstance(v, (int, float)) and not isinstance(v, bool) else None


# ---------------------------------------------------------------------------------------------------------------
# The service-side walls and the plausibility rules of async's MCP runner (suite_mcp_runner.py at 31a2e41f on
# feat/suite-mcp-runner-vendored: peer review R1-R3 and rule (c) as of e334715b, then 62ee4214 rule 1, 2a8acbaf rule
# 2, f94b079d rule 4, f52846fb rule 4's read floor; host_state_plausible, wall_rule_holds, service_wall, client_wall,
# client_wall_premature, timeout_statement, read_bound_below_floor, resolve_timeout_least), mirrored for the offline
# re-evaluation and the comparator: the same thresholds, the same calls, the same verdict_detail strings (lanes S5d,
# S5e). Rule (c): where no threshold or no service-side wall is known the plausibility is undecidable, and the
# execution is inconclusive for that reason (not a host state, not by design, not failed on that ground). Rule 1
# (62ee4214): a /resolve deadline on a discovery without a budget is implausible (discovery_without_budget). Rule 2
# (2a8acbaf): the client's wall is an upper bound of an unknown service-side wall, for ANONYMOUS MCP calls only
# (mcp_auth_of: the line's or raw record's `mcp_auth`; unknown on a record before 2a8acbaf unless the re-evaluation
# was told, reevaluate.py --mcp-auth).
# Rule 4 (f94b079d): a resolve_timeout is judged against the bound it states (timeout_statement), else the legacy
# bound; f52846fb: a stated READ bound below the service's own floor for a known budget (budget + reserve + 10 s -
# 1 s) is no host state whatever the wall (read_bound_below_floor). A stated DEADLINE window is what the exchange
# budget leaves after the token's acquire and the claim: deadline_window_verdict (one definition for the runner too)
# places its origin with the service's stated pre_exchange_ms, else the failing attempt's span (deadline_evidence),
# else it is undecidable below the read floor. Where the runner recorded what it judged with (raw `client_wall`,
# line `client_wall_ms`; raw / line `timeout_terms`, `deadline_evidence`), the recorded value is read and the
# recomputed one cross-checked (a mismatch is a warning).
# Every bound a rule compares with comes from the request (the scenario's MCP arguments), the host's capabilities and
# the service's settings (References: the run's recorded live documents, else the by-design table's own
# host_capabilities / service_settings blocks), never from the failing answer: a bound the answer states is checked
# against them (stated_bound_problem, tool_error_problem, the wall rule's stated reserve), not believed.

# /resolve's finalisation reserve when the refusal states none (app/traversal/resolve.py DEFAULT_FINALIZE_RESERVE_MS)
RESOLVE_RESERVE_MS = 250.0
# the least allowance the service's client timeout of a budgeted /resolve keeps over budget + reserve
# (app/traversal/resolve.py client_timeout_s, MIN_BUDGET_ALLOWANCE_S = 10 s)
RESOLVE_MIN_ALLOWANCE_MS = 10000.0
# the timeout a resolve_timeout states it applied: the first "within <N> s" of its text
STATED_TIMEOUT = re.compile(r"\bwithin ([0-9]+(?:\.[0-9]+)?(?:e[+-]?[0-9]+)?) s\b")
# the slack a resolve_timeout's least wall gets (clock and phase differences between the service and the client)
STATED_TIMEOUT_SLACK_MS = 1000.0
# the /resolve client's write-phase cap (app/traversal/backend.py _timeout: write = min(10 s, the client timeout)):
# a legacy resolve_timeout (no stated phase) is a host state from min(this, its read-timeout terms) - 1 s on (rule 4;
# the runner's RESOLVE_WRITE_CAP_MS)
RESOLVE_WRITE_TIMEOUT_MS = 10000.0
# the phases a resolve_timeout states (app/traversal/resolve.py timeout_terms: timeout_phase, bound_s, elapsed_ms)
RESOLVE_TIMEOUT_PHASES = ("read", "write", "deadline")
# the /resolve client's connect and pool phase caps (app/traversal/backend.py _timeout, TRAVERSAL_HTTP_CONNECT_TIMEOUT_S
# and TRAVERSAL_HTTP_POOL_TIMEOUT_S; with the write cap and the read timeout t they make the exchange budget,
# exchange_budget_s) and the host token's acquire timeout (app/traversal/resolve.py ACQUIRE_TIMEOUT_MS)
RESOLVE_CONNECT_TIMEOUT_MS = 10000.0
RESOLVE_POOL_TIMEOUT_MS = 5000.0
RESOLVE_ACQUIRE_TIMEOUT_MS = 2000.0
# what a deadline-phase resolve_timeout states of where its deadline began (agreed with the async session,
# 2026-10-15; each rounded to 0.1 ms, absent when unknown): pre_exchange_ms, from the absolute deadline's start (before
# the token's acquire, on the job path before the claim) to the exchange's start; acquire_ms, the acquire's own time
DEADLINE_ORIGIN_KEYS = ("pre_exchange_ms", "acquire_ms")
# bound_s (rounded to 1 ms) and pre_exchange_ms (to 0.1 ms) add up to the exchange budget within this
DEADLINE_TERMS_TOLERANCE_MS = 2.0
# the runner's record of a deadline-phase resolve_timeout's failing attempt (beside timeout_terms; null when unknown)
DEADLINE_EVIDENCE_KEYS = ("exchange_budget_s", "timeout_s", "acquire_timeout_s", "granted_at", "settled_at",
                          "span_ms", "span_source")
# rule 2: the slack the client's wall gets as the upper bound of the service-side wall
CLIENT_WALL_SLACK_MS = 1000.0
# rule 2: the tools whose execution has a client wall (the routes rules R1 and R3 judge)
CLIENT_WALL_TOOLS = ("pattern_search", "traverse_resolve")
# rule 2: the runner's `mcp_auth` (suite_mcp_client.SuiteMcp.mcp_auth): only anonymous calls bound the service-side
# wall (an authenticated identical resolve may reuse an earlier job); anything else, or none, is unknown
MCP_ANONYMOUS = "anonymous"
MCP_API_KEY = "api_key"
MCP_AUTH_VALUES = (MCP_ANONYMOUS, MCP_API_KEY)
# rule 2: what an undecidable verdict_detail adds when the bound was skipped for an API key (the runner's text)
RULE2_SKIPPED_API_KEY = ("rule 2 skipped: the MCP calls present an API key, so an identical request may reuse an "
                         "earlier job and the client wall is no upper bound")
# rule 1: the verdict_detail of a /resolve deadline on a discovery that sent no budget (the runner's NO_DEADLINE_SENT)
DISCOVERY_NO_DEADLINE = "a /resolve deadline on a discovery without a budget: no deadline was sent"
# the host's pattern finalisation reserve when the relayed text states none (SPEC-pattern §10.2: 250 in version 1)
PATTERN_RESERVE_MS = 250.0
# the reserve the host's deadline text states ("... (60000 ms, the finalisation reserve of 250 ms included)")
STATED_RESERVE = re.compile(r"finali[sz]ation reserve of ([0-9]+(?:\.[0-9]+)?) ms")
# a pattern task's held-slot line (the service's app/queue/held_span.py LINE_FORMAT)
SLOT_LINE = re.compile(r"^⧗ Slot held for .+? ([0-9.]+)s \(att ([^,]+), (\w+), acq ([0-9.]+), rel ([0-9.]+), "
                       r"http (\d+) fail (\d+) ([0-9.]+)s backoff ([0-9.]+)s, (\w+)\)")
WALL_BACKEND = "mcp_phases.backend_ms"
WALL_DEADLINE_TASK = "the deadline-ended task's one failed HTTP try (diagnostic log 'Slot held ... http 1 fail 1')"
WALL_ATTEMPT = "the settled attempt's settled_at - dispatched_at (GET /admin/traversal/jobs/{id} attempts[])"
WALL_JOB = "the job's finished_at - started_at (traverse_status)"
WALL_NONE = "no service-side timing"
# rule (c): the verdict_detail of an execution whose plausibility cannot be decided
UNDECIDABLE_NO_WALL = "no service-side timing: plausibility undecidable"
UNDECIDABLE_NO_THRESHOLD = "no threshold known: plausibility undecidable"
# a deadline window below the read floor with no span of its attempt and no stated origin (deadline_window_verdict)
UNDECIDABLE_NO_ORIGIN = "no attempt span places the deadline's origin: plausibility undecidable"
UNDECIDABLE = (UNDECIDABLE_NO_WALL, UNDECIDABLE_NO_THRESHOLD, UNDECIDABLE_NO_ORIGIN)
WALL_RULE_VALUES = ("budget_minus_reserve",)
# a stated resolve bound is rounded to 1 ms by the service (app/traversal/resolve.py timeout_terms)
STATED_BOUND_ROUNDING_S = 0.001


# ---------------------------------------------------------------------------------------------------------------
# References to the configuration: {"host": <path>} into the host's /capabilities document, {"service": <path>} into
# the service's GET /traverse/capabilities document (its root: service_settings.*, presets.<p>.budgets.<knob>,
# service_maxima.*, backends.<db>.resolve.time_budget.*). Resolution order, agreed with async's runner: the run's
# recorded live document when it states the path, else the by-design table's own block (host_capabilities /
# service_settings, shaped as those documents' roots: the versioned recorded configuration), else unknown. A path
# both state with different values is a mismatch, judged with the live value.

REFERENCE_NAMESPACES = ("host", "service")
REFERENCE_SOURCES = ("live", "table", "unknown")
# the table's block of each namespace
REFERENCE_BLOCKS = {"host": "host_capabilities", "service": "service_settings"}
# A path (a request predicate's subject, a reference): keys joined by "." (a key: any characters but . [ ] { } and
# whitespace), "[n]" (element n of a list) or "[*]" (every element of a list) after a key, and in a request path the
# whole-key placeholder "{field}"; nothing else (no empty key, no leading "[", no "a[x]").
_PATH_KEY = r"[^.\[\]{}\s]+"
_PATH_STEP = re.compile(r"(%s|\{field\})((?:\[(?:[0-9]+|\*)\])*)" % _PATH_KEY)
_PATH_INDEX = re.compile(r"\[([0-9]+|\*)\]")


def parse_path(path, field=False, star=True):
    """The steps of a dotted |path| (keys, "[n]", "[*]") by the grammar at _PATH_KEY, or None when it does not
    parse; "{field}" is a key only with |field|, "[*]" only with |star|."""
    if not isinstance(path, str) or not path:
        return None
    out = []
    for piece in path.split("."):
        m = _PATH_STEP.fullmatch(piece)
        if m is None or (m.group(1) == "{field}" and not field):
            return None
        out.append(m.group(1))
        for idx in _PATH_INDEX.findall(m.group(2)):
            if idx == "*" and not star:
                return None
            out.append("[%s]" % idx)
    return out


def path_parts(path):
    """The steps of a dotted path (parse_path with "[*]"; [] when it does not parse)."""
    return parse_path(path) or []


def lookup(doc, path):
    """(True, the value) at a dotted |path| of |doc| ("[n]" an element of a list; a path that does not parse, or
    names "[*]", is no lookup), else (False, None)."""
    parts = parse_path(path, star=False)
    if parts is None:
        return False, None
    cur = doc
    for part in parts:
        if part.startswith("["):
            if not isinstance(cur, list) or int(part[1:-1]) >= len(cur):
                return False, None
            cur = cur[int(part[1:-1])]
        elif isinstance(cur, dict) and part in cur:
            cur = cur[part]
        else:
            return False, None
    return True, cur


def override_path_listed(name, entries):
    """Whether the dotted path |name| (an MCP override, or the field a refusal names, less its prefix) is one of the
    by-design table's |entries|: the path itself, or "x.*" for a whole key, which lists the bare section x as well
    (the table's own text and the service read {"bounds": 5} as basic). The one matcher of compare's field
    predicates and the gate policy's basic overrides (SD-095)."""
    return any(name == e or (str(e).endswith(".*") and (name == e[:-2] or name.startswith(e[:-1])))
               for e in entries)


def json_equal(a, b):
    """JSON equality: a boolean equals only the same boolean (true is not 1), numbers by value (250 and 250.0 are
    one; never a numeric string), strings, null, and lists / objects element by element and key by key."""
    if isinstance(a, bool) or isinstance(b, bool):
        return isinstance(a, bool) and isinstance(b, bool) and a == b
    if _number_or_none(a) is not None or _number_or_none(b) is not None:
        return _number_or_none(a) is not None and _number_or_none(b) is not None and float(a) == float(b)
    if isinstance(a, list) or isinstance(b, list):
        return isinstance(a, list) and isinstance(b, list) and len(a) == len(b) and \
            all(json_equal(x, y) for x, y in zip(a, b))
    if isinstance(a, dict) or isinstance(b, dict):
        return isinstance(a, dict) and isinstance(b, dict) and set(a) == set(b) and \
            all(json_equal(a[k], b[k]) for k in a)
    return a == b


def same_value(a, b):
    """Whether two configuration values are the same (json_equal)."""
    return json_equal(a, b)


# The host's identity, a part of the host namespace and of a service document's identity (agreed with async's runner
# on 2026-10-15): the MCP runner's host_document states it beside its pattern and resolve parts (from the traverse
# document's backends.<db>.identity: {index_ns, index_fp, index_meta_fp, release}); a paired run's direct arms record
# release in the host's /capabilities and index_fp, index_meta_fp in its /traverse/capabilities. A host redeploy with
# unchanged limits is another deployment.
HOST_IDENTITY_KEYS = ("release", "index_fp", "index_meta_fp")
# The identity of a service configuration record (the MCP runner's run.json service_capabilities and each
# invocations[].service_capabilities: {document (GET /traverse/capabilities), host_document, sha256, read_at, base,
# database, mismatches, error}): the sha256 of the canonical JSON of its configuration view, ONE definition both the
# suite and async's runner (which vendors this module) compute (service_document_sha256). The rest of the document
# (probe times, cached flags, messages, counters) changes between the invocations of one deployment.
SERVICE_DOCUMENT_VIEW_KEYS = ("service_settings", "service_maxima", "presets")
SERVICE_BACKEND_VIEW_KEYS = ("budgets", "resolve", "effective_maxima", "limits")


def host_document_of(traverse_doc, pattern_doc, database):
    """The host's configuration as the service states it, in the shape of the host's /capabilities document (the
    {"host": path} namespace): {pattern: GET /pattern/capabilities/<db>'s capabilities.raw, resolve: {time_budget: GET
    /traverse/capabilities' backends.<db>.budgets.resolve_time_budget}, release, index_fp, index_meta_fp: from
    backends.<db>.identity (HOST_IDENTITY_KEYS)}, each part where a read states it (an identity value that is null is
    absent); None when none does."""
    out = {}
    caps = pattern_doc.get("capabilities") if isinstance(pattern_doc, dict) else None
    raw = caps.get("raw") if isinstance(caps, dict) else None
    if isinstance(raw, dict):
        out["pattern"] = raw
    backends = traverse_doc.get("backends") if isinstance(traverse_doc, dict) else None
    backend = backends.get(database) if isinstance(backends, dict) and database else None
    budgets = backend.get("budgets") if isinstance(backend, dict) else None
    terms = budgets.get("resolve_time_budget") if isinstance(budgets, dict) else None
    if isinstance(terms, dict):
        out["resolve"] = {"time_budget": terms}
    identity = backend.get("identity") if isinstance(backend, dict) else None
    for key in HOST_IDENTITY_KEYS:
        if isinstance(identity, dict) and identity.get(key) is not None:
            out[key] = identity[key]
    return out or None


def service_document_view(document, host_document=None, database=None):
    """What a service configuration record's identity covers: {service_settings, service_maxima, presets (each where
    the document states it), backends: {<database>: {budgets, resolve, effective_maxima, limits} (where stated)} of
    the record's database only, host_document (where recorded: its pattern and resolve parts and the host's identity,
    HOST_IDENTITY_KEYS)}. A key the document lacks is absent from the view, never null; another database's entry is
    no part of it."""
    doc = document if isinstance(document, dict) else {}
    view = dict((k, doc[k]) for k in SERVICE_DOCUMENT_VIEW_KEYS if k in doc)
    backends = doc.get("backends")
    backend = backends.get(database) if isinstance(backends, dict) and database else None
    if isinstance(backend, dict):
        view["backends"] = {database: dict((k, backend[k]) for k in SERVICE_BACKEND_VIEW_KEYS if k in backend)}
    if isinstance(host_document, dict):
        view["host_document"] = host_document
    return view


def service_document_sha256(record):
    """The identity of a service_capabilities record: sha256 of canonical_json(service_document_view(its document,
    its host_document, its database)); None when the view is empty (nothing was read). A read of the host document
    alone (the traverse read failed) has the identity of that document."""
    if not isinstance(record, dict):
        return None
    view = service_document_view(record.get("document"), record.get("host_document"), record.get("database"))
    return hashlib.sha256(canonical_json(view)).hexdigest() if view else None


class References(object):
    """Resolves {"host": path} / {"service": path} references: the live document of the namespace (the run's recorded
    one) when it states the path, else the table's block, else unknown. |table| and |live| map a namespace to a
    document (a dict) or None; |origins| a namespace to where its live document came from (named in reports)."""

    def __init__(self, table=None, live=None, origins=None):
        self.table = dict((ns, (table or {}).get(ns) if isinstance((table or {}).get(ns), dict) else None)
                          for ns in REFERENCE_NAMESPACES)
        self.live = dict((ns, (live or {}).get(ns) if isinstance((live or {}).get(ns), dict) else None)
                         for ns in REFERENCE_NAMESPACES)
        self.origins = dict(origins or {})

    def resolve(self, ns, path):
        """(value, source, mismatch): source "live" | "table" | "unknown"; mismatch {namespace, ref, live, table} when
        both documents state |path| with different values (the live value is returned), else None."""
        found_t, value_t = lookup(self.table.get(ns), path) if self.table.get(ns) is not None else (False, None)
        if self.live.get(ns) is not None:
            found_l, value_l = lookup(self.live[ns], path)
            if found_l:
                mismatch = None
                if found_t and not same_value(value_t, value_l):
                    mismatch = {"namespace": ns, "ref": str(path), "live": value_l, "table": value_t}
                return value_l, "live", mismatch
        if found_t:
            return value_t, "table", None
        return None, "unknown", None

    def number(self, ns, path, default=None):
        """The numeric value of a reference, else |default| (unknown, or not a number)."""
        value, source, _ = self.resolve(ns, path)
        num = _number_or_none(value) if source != "unknown" else None
        return default if num is None else num


_DEFAULT_REFERENCES = []


def default_references():
    """References to the by-design table next to this module (service_by_design.json: its host_capabilities and
    service_settings blocks), no live document; empty when the table or its blocks are missing (every rule then
    reads its built-in value)."""
    if not _DEFAULT_REFERENCES:
        table = {}
        try:
            with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "service_by_design.json")) as fh:
                data = json.load(fh)
            table = dict((ns, data.get(block)) for ns, block in REFERENCE_BLOCKS.items())
        except (OSError, ValueError):
            pass
        _DEFAULT_REFERENCES.append(References(table))
    return _DEFAULT_REFERENCES[0]



# The run's recorded configuration, the live documents the references read first (compare's and merge_run's: SD-035).
# The MCP runner's run.json record of GET /traverse/capabilities at run start (async's runner: {document,
# host_document, sha256, read_at, base, database, mismatches}, on a failed read the error instead of the document;
# each invocations[] the same); host_document has the shape of the host's /capabilities document as the service read
# it: "pattern" the /pattern/capabilities/<db> capabilities.raw, "resolve" {"time_budget": the host's block as
# received, document.backends.<db>.budgets.resolve_time_budget}
SERVICE_CAPABILITIES_KEY = "service_capabilities"
SERVICE_HOST_DOCUMENT_KEY = "host_document"
# the host's /capabilities document as the direct runners record it (run.json capabilities)
HOST_CAPABILITIES_PATH = "/capabilities"


def direct_host_document(caps, run=None):
    """The host namespace's document of a run with a direct arm: the host's recorded /capabilities, with the host's
    identity (suite_common.HOST_IDENTITY_KEYS) it does not state at its top level taken from the identity the arm
    recorded: its /traverse/capabilities (index_fp, index_meta_fp), else the run.json's own release / index_fp;
    the same paths as the MCP runner's host_document."""
    doc = dict(caps.get(HOST_CAPABILITIES_PATH) or {})
    traverse = caps.get("/traverse/capabilities") if isinstance(caps.get("/traverse/capabilities"), dict) else {}
    for key in HOST_IDENTITY_KEYS:
        if doc.get(key) is None:
            value = traverse.get(key)
            if value is None and isinstance(run, dict):
                value = run.get(key)
            if value is not None:
                doc[key] = value
    return doc


def recorded_documents(run_dir):
    """{namespace: {document, origin, read_at, sha256} or None} of a run: "service" the MCP arm's run.json
    service_capabilities.document (arms/mcp/run.json, else the merged run.json's arms.mcp; sha256 recomputed,
    suite_common.service_document_sha256, the record's own as sha256_recorded), "host" the /capabilities document a
    direct arm recorded (the merged run.json's capabilities, else arms/<direct>/run.json; with the host's identity:
    direct_host_document), else the host document the MCP runner recorded beside the service's
    (service_capabilities.host_document: an MCP-only run reads the host's references live from it); None
    where the run recorded none (a run before async's runner read it, or a failed read: then `error` says why). A
    path the live document lacks reads the table's block, else unknown (suite_common.References)."""
    out = {"service": None, "host": None}
    merged = read_json(os.path.join(run_dir, "run.json"), {}) or {}
    mcp_run = read_json(os.path.join(run_dir, "arms", "mcp", "run.json"), None)
    if not isinstance(mcp_run, dict):
        mcp_run = ((merged.get("arms") or {}).get("mcp")) if isinstance(merged.get("arms"), dict) else None
    rec = (mcp_run or {}).get(SERVICE_CAPABILITIES_KEY) if isinstance(mcp_run, dict) else None
    if isinstance(rec, dict):
        doc = rec.get("document")
        out["service"] = {"document": doc if isinstance(doc, dict) else None,
                          "origin": "arms/mcp/run.json %s" % SERVICE_CAPABILITIES_KEY,
                          "read_at": rec.get("read_at"), "sha256": service_document_sha256(rec),
                          "sha256_recorded": rec.get("sha256"),
                          "base": rec.get("base"),
                          "error": rec.get("error") if not isinstance(doc, dict) else None}
    caps = merged.get("capabilities") if isinstance(merged.get("capabilities"), dict) else None
    origin, owner = "run.json capabilities", merged
    if not (caps and isinstance(caps.get(HOST_CAPABILITIES_PATH), dict)):
        caps = None
        for name in DIRECT_ARMS:
            arm_run = read_json(os.path.join(run_dir, "arms", name, "run.json"), None)
            if isinstance(arm_run, dict) and isinstance((arm_run.get("capabilities") or {}).get(
                    HOST_CAPABILITIES_PATH), dict):
                caps, origin, owner = arm_run["capabilities"], "arms/%s/run.json capabilities" % name, arm_run
                break
    if caps and isinstance(caps.get(HOST_CAPABILITIES_PATH), dict):
        out["host"] = {"document": direct_host_document(caps, owner),
                       "origin": "%s[%s]" % (origin, HOST_CAPABILITIES_PATH), "read_at": None, "sha256": None}
    elif isinstance(rec, dict) and isinstance(rec.get(SERVICE_HOST_DOCUMENT_KEY), dict):
        out["host"] = {"document": rec[SERVICE_HOST_DOCUMENT_KEY],
                       "origin": "arms/mcp/run.json %s.%s" % (SERVICE_CAPABILITIES_KEY, SERVICE_HOST_DOCUMENT_KEY),
                       "read_at": rec.get("read_at"), "sha256": None}
    return out


def run_references(run_dir, blocks=None):
    """References of a run: its recorded live documents (recorded_documents) over |blocks| (namespace -> the table's
    host_capabilities / service_settings block; default: the table next to this module, default_references)."""
    docs = recorded_documents(run_dir)
    live = dict((ns, (d or {}).get("document")) for ns, d in docs.items())
    origins = dict((ns, (d or {}).get("origin")) for ns, d in docs.items() if d)
    return References(default_references().table if blocks is None else blocks, live, origins)


def record_answer(line, doc):
    """An MCP execution's answer as the runner's matcher reads it ({http, error, body}): its raw record's |doc|
    answer where recorded, else its line's http and error."""
    ans = (doc or {}).get("answer") if isinstance((doc or {}).get("answer"), dict) else None
    if ans is None:
        err = line.get("error") if isinstance(line.get("error"), dict) and line["error"] else None
        ans = {"http": line.get("http"), "error": err, "body": None}
    return ans


def record_error(line, doc):
    """An MCP execution's error {code, text, field}: its line's when that states one, else its raw record |doc|'s
    answer error, else the record's top-level one; None when none states one."""
    err = line.get("error")
    if isinstance(err, dict) and (err.get("code") or err.get("text")):
        return err
    if isinstance(doc, dict):
        for cand in ((doc.get("answer") or {}).get("error") if isinstance(doc.get("answer"), dict) else None,
                     doc.get("error")):
            if isinstance(cand, dict) and (cand.get("code") or cand.get("text")):
                return cand
    return None


def host_state_of_record(scenario, line, doc, refs=None):
    """host_state_verdict of an MCP execution from its raw record |doc| as compare and merge_run read a line recorded
    without the runner's inconclusive_reason (SD-035): the error of its answer (record_answer), else record_error, the
    tool error behind that answer, and |refs| (the run's References: run_references)."""
    ans = record_answer(line, doc)
    return host_state_verdict(scenario or {}, line, doc, ans.get("error") or record_error(line, doc),
                              tool_error_of(doc, ans), refs=refs)

# The terms the plausibility rules read: (namespace, path, the built-in value when the reference is unknown).
# host: the host's finalisation reserves and /resolve budget cap; service: the /resolve client's settings
# (app/traversal/resolve.py client_timeout_s, plan_budget; app/traversal/backend.py _timeout_phase).
TERMS = {
    "pattern_reserve_ms": ("host", "pattern.finalize_reserve_ms", PATTERN_RESERVE_MS),
    "resolve_reserve_ms": ("host", "resolve.time_budget.finalize_reserve_ms", RESOLVE_RESERVE_MS),
    "resolve_max_time_ms": ("host", "resolve.time_budget.max_time_ms", None),
    "timeout_s": ("service", "service_settings.resolve.timeout_s", 20.0),
    "discover_timeout_s": ("service", "service_settings.resolve.discover_timeout_s", 120.0),
    "allowance_s": ("service", "service_settings.resolve.allowance_s", RESOLVE_MIN_ALLOWANCE_MS / 1000.0),
    "default_time_budget_ms": ("service", "service_settings.resolve.default_time_budget_ms", 60000.0),
    "budget_ceiling_ms": ("service", "service_settings.resolve.budget_ceiling_ms", 240000.0),
    "write_timeout_s": ("service", "service_settings.resolve.http_timeouts_s.write", RESOLVE_WRITE_TIMEOUT_MS / 1000.0),
    "connect_timeout_s": ("service", "service_settings.resolve.http_timeouts_s.connect",
                          RESOLVE_CONNECT_TIMEOUT_MS / 1000.0),
    "pool_timeout_s": ("service", "service_settings.resolve.http_timeouts_s.pool", RESOLVE_POOL_TIMEOUT_MS / 1000.0),
    "acquire_timeout_ms": ("service", "service_settings.resolve.acquire_timeout_ms", RESOLVE_ACQUIRE_TIMEOUT_MS),
}


def term(refs, name):
    """A plausibility rule's term (TERMS) from |refs| (default_references()), its built-in value when unknown."""
    ns, path, builtin = TERMS[name]
    return (refs or default_references()).number(ns, path, builtin)


def _mcp_args(scenario):
    args = (((scenario or {}).get("mcp") or {}).get("arguments")) or {}
    return args if isinstance(args, dict) else {}


def resolve_default_budget_ms(scenario, refs=None):
    """The time budget the service sends on an explicit-labels /resolve whose request names none (app/traversal/
    resolve.py plan_budget): the service document's backends.<database>.resolve.time_budget.default_ms, else
    min(service_settings.resolve.default_time_budget_ms, the host's resolve.time_budget.max_time_ms,
    service_settings.resolve.budget_ceiling_ms)."""
    refs = refs or default_references()
    db = _mcp_args(scenario).get("database")
    if isinstance(db, str) and db:
        stated = refs.number("service", "backends.%s.resolve.time_budget.default_ms" % db)
        if stated is not None:
            return stated
    caps = [term(refs, "default_time_budget_ms"), term(refs, "resolve_max_time_ms"), term(refs, "budget_ceiling_ms")]
    return min(x for x in caps if x is not None)


def resolve_budget(scenario, refs=None):
    """(the time budget the service sent with a /resolve, whose it is) for the plausibility rules: the request's
    (the MCP arguments' time_budget_ms); an explicit-labels resolve without one the service's default
    (resolve_default_budget_ms); a discovery without one, or a request naming neither labels nor a budget, none:
    (None, why)."""
    asked = mcp_budget_ms(scenario)
    if asked is not None:
        return asked, "the request's"
    args = _mcp_args(scenario)
    if args.get("discover"):
        return None, "a discovery without a budget"
    if args.get("labels"):
        return resolve_default_budget_ms(scenario, refs), "the service's default"
    return None, "no budget"


def resolve_client_timeout_s(scenario, refs=None):
    """The service's client timeout t of the scenario's /resolve call (app/traversal/resolve.py client_timeout_s): the
    mode's base (service_settings.resolve.timeout_s, discover_timeout_s for a discovery), with a budget B (resolve_budget)
    max(base, (B + the host's reserve) / 1000 + max(service_settings.resolve.allowance_s, 10))."""
    base = term(refs, "discover_timeout_s") if _mcp_args(scenario).get("discover") else term(refs, "timeout_s")
    budget, _ = resolve_budget(scenario, refs)
    if budget is None:
        return base
    allowance = max(term(refs, "allowance_s"), RESOLVE_MIN_ALLOWANCE_MS / 1000.0)
    return max(base, (budget + term(refs, "resolve_reserve_ms")) / 1000.0 + allowance)


def tool_error_of(doc, answer):
    """The MCP tool's own error object behind an answer's error (its field, budget, reserve, ...), from the recorded
    calls' payloads of a raw MCP record, as async's runner reads it (suite_mcp_runner.py raw_error), else None."""
    code = ((answer or {}).get("error") or {}).get("code") if isinstance((answer or {}).get("error"), dict) else None
    for call in reversed((doc or {}).get("calls") or []):
        payload = call.get("payload") if isinstance(call, dict) else None
        err = payload.get("error") if isinstance(payload, dict) else None
        if isinstance(err, dict) and err.get("code") == code:
            return err
        rerr = payload.get("resolve_error") if isinstance(payload, dict) else None
        detail = rerr.get("detail") if isinstance(rerr, dict) else None
        if isinstance(detail, dict) and detail.get("code") == code:
            return detail
    return None


def mcp_tool_of(scenario):
    """The scenario's MCP tool, else the tool its route is executed with."""
    tool = ((scenario or {}).get("mcp") or {}).get("tool")
    if tool:
        return str(tool)
    return {"POST /pattern": "pattern_search", "POST /resolve": "traverse_resolve",
            "POST /traverse": "traverse_fetch"}.get(str((scenario or {}).get("route") or ""), "")


def mcp_budget_ms(scenario):
    """The time budget the MCP arm sent (mcp.arguments.time_budget_ms), or None."""
    return _number_or_none((((scenario or {}).get("mcp") or {}).get("arguments") or {}).get("time_budget_ms"))


def request_budget_ms(scenario):
    """A pattern scenario's request time_budget_ms: the MCP arguments' (what the arm sent), else the direct
    request's top-level one."""
    asked = mcp_budget_ms(scenario)
    if asked is not None:
        return asked
    req = (scenario or {}).get("request")
    return _number_or_none(req.get("time_budget_ms")) if isinstance(req, dict) else None


def _diag_tasks(diag):
    if not isinstance(diag, dict):
        return []
    tasks = [t for t in diag.get("regular_tasks") or [] if isinstance(t, dict)]
    for group in diag.get("split_task_groups") or []:
        tasks.extend(t for t in (group or {}).get("children") or [] if isinstance(t, dict))
    return tasks


def _slot_figures(task):
    found = None
    for line in task.get("diagnostic_log") or []:
        m = SLOT_LINE.match(line) if isinstance(line, str) else None
        if m:
            found = {"http_tries": int(m.group(6)), "failed_tries": int(m.group(7)),
                     "failed_post_s": float(m.group(8))}
    return found


def deadline_task_ms(doc):
    """A pattern search's deadline-ended task wall (R1's fallback): of the recorded diagnostic log, each task whose
    error starts pattern_deadline, its last held-slot line's failed POST seconds when that attempt sent one HTTP try
    and it failed; the smallest; None when no task ended so or any lacks such a line."""
    diag = next((c.get("payload") for c in reversed((doc or {}).get("calls") or []) if isinstance(c, dict)
                 and str(c.get("tool") or "").startswith("GET /searches/")
                 and str(c.get("tool")).endswith("/diagnostic-log")), None)
    ended = [t for t in _diag_tasks(diag) if str(t.get("error") or "").startswith("pattern_deadline")]
    walls = []
    for task in ended:
        fig = _slot_figures(task)
        if fig is None or fig["http_tries"] != 1 or fig["failed_tries"] != 1:
            return None
        walls.append(fig["failed_post_s"] * 1000.0)
    return min(walls) if walls else None


def settled_attempt_ms(sources):
    """The span_ms of the last settled attempt in a raw record's phase_sources.attempts, else None."""
    rows = (sources or {}).get("attempts")
    spans = [_number_or_none(r.get("span_ms")) for r in rows if isinstance(r, dict)] if isinstance(rows, list) else []
    known = [x for x in spans if x is not None]
    return known[-1] if known else None


def runner_epoch(value):
    """An ISO timestamp as epoch seconds as async's runner reads phase_sources (epoch_of: a non-empty string, every
    'Z' read as +00:00, naive = UTC), else None."""
    if not isinstance(value, str) or not value:
        return None
    try:
        t = datetime.datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        return None
    if t.tzinfo is None:
        t = t.replace(tzinfo=datetime.timezone.utc)
    return t.timestamp()


def _ms_between(a, b):
    """(b - a) in ms rounded to 0.1, None when either is unknown (the runner's ms)."""
    return None if a is None or b is None else round((b - a) * 1000.0, 1)


def service_wall(scenario, line, doc=None):
    """(the service-side wall of an MCP execution in ms, its source) for a host-state or wall rule (R3, R1; the
    runner's service_wall), in order: the line's mcp_phases.backend_ms (else the raw record's); for a pattern search
    the deadline-ended task's HTTP try (deadline_task_ms of the raw record); for a traversal job the last settled
    attempt's span, then the job's finished_at - started_at (the raw record's phase_sources, in ms to 0.1); else
    (None, WALL_NONE). Never the client's wall."""
    backend = _number_or_none(((line or {}).get("mcp_phases") or {}).get("backend_ms"))
    if backend is None:
        backend = _number_or_none(((doc or {}).get("mcp_phases") or {}).get("backend_ms"))
    if backend is not None:
        return backend, WALL_BACKEND
    sources = (doc or {}).get("phase_sources") if isinstance((doc or {}).get("phase_sources"), dict) else {}
    tool = mcp_tool_of(scenario)
    if tool == "pattern_search":
        task = deadline_task_ms(doc)
        if task is not None:
            return task, WALL_DEADLINE_TASK
    if tool in ("traverse_resolve", "traverse_fetch"):
        span = settled_attempt_ms(sources)
        if span is not None:
            return span, WALL_ATTEMPT
        job = _ms_between(runner_epoch(sources.get("started_at")), runner_epoch(sources.get("finished_at")))
        if job is not None:
            return job, WALL_JOB
    return None, WALL_NONE


def _warn(warnings, text):
    if warnings is not None and text not in warnings:
        warnings.append(text)


def mcp_auth_recorded(line, doc=None):
    """Whether the runner recorded `mcp_auth` for the execution (2a8acbaf on: the key on every line and in the raw
    record, its value the client's SuiteMcp.mcp_auth, null only from a client without one): a record before
    2a8acbaf has the key in neither."""
    return any(isinstance(src, dict) and "mcp_auth" in src for src in (line, doc))


def mcp_auth_of(line, doc=None):
    """(the MCP arm's auth as rule 2 reads it, where it came from): the `mcp_auth` the runner recorded (2a8acbaf:
    "anonymous" | "api_key" on every line and in the raw record; the line's, else the raw record's), else, for a
    record before 2a8acbaf (no `mcp_auth` key in either: mcp_auth_recorded), the re-evaluation's stated assumption
    (`reevaluation.mcp_auth_assumed` of a re-evaluated line: reevaluate.py --mcp-auth), else (None, None): unknown,
    so rule 2 decides nothing (the runner's client_wall_premature with an auth other than anonymous / api_key; a
    recorded null stays unknown, as the runner judged it)."""
    for src, name in ((line, "line"), (doc, "raw record")):
        v = src.get("mcp_auth") if isinstance(src, dict) else None
        if v is not None:
            return v, name
    if mcp_auth_recorded(line, doc):
        return None, None
    assumed = (((line or {}).get("reevaluation") or {}) if isinstance((line or {}).get("reevaluation"), dict)
               else {}).get("mcp_auth_assumed")
    if assumed is not None:
        return assumed, "the re-evaluation's --mcp-auth"
    return None, None


def client_wall(scenario, doc):
    """The execution's client wall (rule 2; the runner's client_wall, 2a8acbaf) from the raw record's calls (tool,
    t0 / t1 epoch s, payload, timed): {ms, first_submit_t0, first_submit_call, answer_t1, answer_call, answer_tool},
    ms = answer_t1 - first_submit_t0 in ms (0.1), the *_call indexes into calls. The first submit: the execution's
    first call of its submit tool (a busy or idempotency_conflict retry that joined an earlier search or job counts
    from that first send; its t0 unknown: no client wall). The call that returned the answer: traverse_resolve, the
    last timed traverse_results (the job path), else the last traverse_resolve; pattern_search, the FIRST get_status
    whose payload is terminal, else the last pattern_search. None for another tool, or when either instant is
    unknown."""
    tool = mcp_tool_of(scenario)
    if tool not in CLIENT_WALL_TOOLS:
        return None
    calls = [c if isinstance(c, dict) else {} for c in (doc or {}).get("calls") or []]
    submits = [i for i, c in enumerate(calls) if c.get("tool") == tool]
    if not submits:
        return None
    if tool == "traverse_resolve":
        reads = [i for i, c in enumerate(calls) if c.get("tool") == "traverse_results" and c.get("timed", True)]
        answer = reads[-1] if reads else submits[-1]
    else:
        terminal = [i for i, c in enumerate(calls) if c.get("tool") == "get_status"
                    and isinstance(c.get("payload"), dict) and c["payload"].get("terminal") is True]
        answer = terminal[0] if terminal else submits[-1]
    t0, t1 = _number_or_none(calls[submits[0]].get("t0")), _number_or_none(calls[answer].get("t1"))
    if t0 is None or t1 is None:
        return None
    return {"ms": round((t1 - t0) * 1000.0, 1), "first_submit_t0": t0, "first_submit_call": submits[0],
            "answer_t1": t1, "answer_call": answer, "answer_tool": calls[answer].get("tool")}


def client_wall_ms(scenario, line, doc=None, warnings=None):
    """The client wall (ms) rule 2 reads: what the runner recorded (2a8acbaf: the raw record's `client_wall` {ms,
    ...} or null, else the line's `client_wall_ms`), else client_wall of the raw record's calls; a recorded value
    the recomputation does not reproduce is kept and |warnings| (a list, optional) says so."""
    computed = client_wall(scenario, doc)
    computed_ms = None if computed is None else computed["ms"]
    if isinstance(doc, dict) and "client_wall" in doc:
        rec, where = doc.get("client_wall"), "the raw record's client_wall"
        rec = rec.get("ms") if isinstance(rec, dict) else rec
    elif isinstance(line, dict) and "client_wall_ms" in line:
        rec, where = line.get("client_wall_ms"), "the line's client_wall_ms"
    else:
        return computed_ms
    if rec is not None and _number_or_none(rec) is None:
        _warn(warnings, "%s %r is no number: recomputed from the calls (%s)" % (where, rec, computed_ms))
        return computed_ms
    rec = _number_or_none(rec)
    if (rec is None) != (computed_ms is None) or (rec is not None and abs(rec - computed_ms) > 0.05):
        _warn(warnings, "%s %s differs from the calls' client wall %s: the recorded one is read" % (
            where, rec, computed_ms))
    return rec


def client_wall_premature(scenario, line, doc, threshold_ms, warnings=None):
    """Rule 2 for an execution whose service-side wall is unknown (the runner's client_wall_premature): (True, why)
    when its MCP calls are anonymous (mcp_auth_of) and its client wall (client_wall_ms) + CLIENT_WALL_SLACK_MS <
    |threshold_ms|, why = "client wall X ms + 1 s slack < threshold Y ms: provably premature"; (False,
    RULE2_SKIPPED_API_KEY) when the calls presented an API key; (False, None) otherwise (auth unknown, no client
    wall, or not provably premature)."""
    auth = mcp_auth_of(line, doc)[0]
    if auth == MCP_API_KEY:
        return False, RULE2_SKIPPED_API_KEY
    if auth != MCP_ANONYMOUS:
        return False, None
    wall = client_wall_ms(scenario, line, doc, warnings)
    if wall is None or not wall + CLIENT_WALL_SLACK_MS < threshold_ms:
        return False, None
    return True, "client wall %.0f ms + 1 s slack < threshold %.0f ms: provably premature" % (wall, threshold_ms)


def undecidable_timing(note):
    """UNDECIDABLE_NO_WALL, with rule 2's note when it has one (the runner's _undecidable_timing)."""
    return UNDECIDABLE_NO_WALL if note is None else "%s; %s" % (UNDECIDABLE_NO_WALL, note)


def undecidable_class(why):
    """The UNDECIDABLE_* text an undecidable reason starts with (the class inconclusive_reason records; rule 2's
    API-key note may follow it), else the reason itself (the runner's undecidable_class)."""
    return next((u for u in UNDECIDABLE if str(why).startswith(u)), why)


def job_error(doc):
    """A resolve job's stored error as its traverse_results returned it (the last such call's
    resolve_error.detail, any code), else None (the runner's job_error)."""
    for call in reversed((doc or {}).get("calls") or []):
        payload = call.get("payload") if isinstance(call, dict) and call.get("tool") == "traverse_results" else None
        rerr = payload.get("resolve_error") if isinstance(payload, dict) else None
        detail = rerr.get("detail") if isinstance(rerr, dict) else None
        if isinstance(detail, dict):
            return detail
    return None


def timeout_statement(doc, error, tool_error):
    """What a resolve_timeout states fired (rule 4; the runner's timeout_statement, f94b079d): {timeout_phase,
    bound_s, elapsed_ms, source} from the FIRST of the tool error (tool_error_of: the inline tool error or a job's
    resolve_error.detail), the answer's error object (|error|), the job's stored error (job_error) that states a
    timeout_phase in RESOLVE_TIMEOUT_PHASES and a positive numeric bound_s, both in that one object, with the
    deadline's origin terms of that object where it states them (DEADLINE_ORIGIN_KEYS: pre_exchange_ms, acquire_ms,
    non-negative numbers; absent otherwise); None when none does (a legacy record)."""
    for source, err in (("tool error", tool_error), ("answer error", error), ("job error", job_error(doc))):
        if not isinstance(err, dict):
            continue
        phase, bound = err.get("timeout_phase"), _number_or_none(err.get("bound_s"))
        if phase in RESOLVE_TIMEOUT_PHASES and bound is not None and bound > 0:
            out = {"timeout_phase": phase, "bound_s": bound, "elapsed_ms": _number_or_none(err.get("elapsed_ms")),
                   "source": source}
            for key in DEADLINE_ORIGIN_KEYS:
                value = _number_or_none(err.get(key))
                if value is not None and value >= 0:
                    out[key] = value
            return out
    return None


def _valid_statement(rec):
    return rec is None or (isinstance(rec, dict) and rec.get("timeout_phase") in RESOLVE_TIMEOUT_PHASES
                           and _number_or_none(rec.get("bound_s")) is not None and rec["bound_s"] > 0)


def timeout_terms(line, doc, error, tool_error, warnings=None):
    """The statement rule 4 reads: what the runner recorded (f94b079d: the raw record's `timeout_terms`, else the
    line's; null = a legacy record), else timeout_statement of the record; a recorded statement the recomputation
    does not reproduce is kept and |warnings| (a list, optional) says so; a malformed one is replaced by the
    recomputation, with a warning. The deadline's origin terms (DEADLINE_ORIGIN_KEYS) a recorded statement lacks are
    taken from the recomputation when the rest agrees (a runner that records the statement without them)."""
    computed = timeout_statement(doc, error, tool_error)
    if isinstance(doc, dict) and "timeout_terms" in doc:
        rec, where = doc.get("timeout_terms"), "the raw record's timeout_terms"
    elif isinstance(line, dict) and "timeout_terms" in line:
        rec, where = line.get("timeout_terms"), "the line's timeout_terms"
    else:
        return computed
    if not _valid_statement(rec):
        _warn(warnings, "%s %s states no valid phase and bound: recomputed (%s)" % (
            where, json.dumps(rec, sort_keys=True), json.dumps(computed, sort_keys=True)))
        return computed
    # the deadline's origin terms a recorded statement does not state (a runner of before them) are read from the
    # record's error object; the rest is compared
    extra = dict((k, v) for k, v in (computed or {}).items() if k in DEADLINE_ORIGIN_KEYS and k not in (rec or {}))
    same = dict((k, v) for k, v in computed.items() if k not in extra) if computed is not None else None
    if json.dumps(rec, sort_keys=True) != json.dumps(same, sort_keys=True):
        _warn(warnings, "%s %s differs from the record's statement %s: the recorded one is read" % (
            where, json.dumps(rec, sort_keys=True), json.dumps(computed, sort_keys=True)))
        return rec
    return dict(rec, **extra) if extra else rec


def read_bound_below_floor(scenario, error, tool_error, doc=None, line=None, warnings=None, refs=None):
    """Rule 4's read floor (the runner's read_bound_below_floor, f52846fb): a resolve_timeout whose statement
    (timeout_terms: the runner's recorded timeout_terms, else timeout_statement) names timeout_phase "read" with a
    bound_s below the service's own read timeout less 1 s, t - 1 s (resolve_client_timeout_s: the mode's base, with a
    budget at least budget + the host's reserve + 10 s; the budget the request's, an explicit-labels resolve's
    without one the service's default, resolve_budget), is no host state: the service sized its own read timeout
    too short, whatever the wall (an unknown one included). Returns "service read bound <bound_s> s below its own
    floor <floor> s" (both %g seconds, the runner's text), else None (another phase, no stated phase, or a bound at
    or above the floor). Write and deadline phases: stated_bound_problem."""
    te = tool_error if isinstance(tool_error, dict) else {}
    stated = timeout_terms(line, doc, error, te, warnings)
    if stated is None or stated.get("timeout_phase") != "read":
        return None
    floor_ms = resolve_client_timeout_s(scenario, refs) * 1000.0 - STATED_TIMEOUT_SLACK_MS
    bound_s = float(stated["bound_s"])
    if bound_s * 1000.0 >= floor_ms:
        return None
    return "service read bound %g s below its own floor %g s" % (bound_s, floor_ms / 1000.0)


def stated_bound_problem(scenario, error, tool_error, doc=None, line=None, warnings=None, refs=None):
    """The write and deadline phases of a stated resolve_timeout (timeout_terms) against the bounds the service
    applies (app/traversal/backend.py _timeout_phase, exchange_budget_s): the write bound is min(the write timeout
    service_settings.resolve.http_timeouts_s.write, t) (t: resolve_client_timeout_s), within the 1 ms the statement
    is rounded to, else "service write bound <b> s is not its write cap <c> s"; the deadline window is the exchange
    budget less what preceded the exchange (the token's acquire, the claim), judged by deadline_window_verdict (its
    not_plausible detail; its undecidable state host_state_verdict reads through deadline_window_state). None for
    another phase, no statement, or a bound the service can have applied."""
    te = tool_error if isinstance(tool_error, dict) else {}
    stated = timeout_terms(line, doc, error, te, warnings)
    if stated is None or stated.get("timeout_phase") not in ("write", "deadline"):
        return None
    bound_s = float(stated["bound_s"])
    t_s = resolve_client_timeout_s(scenario, refs)
    if stated["timeout_phase"] == "write":
        cap_s = min(term(refs, "write_timeout_s"), t_s)
        if abs(bound_s - cap_s) <= STATED_BOUND_ROUNDING_S + 1e-9:
            return None
        return "service write bound %g s is not its write cap %g s" % (bound_s, cap_s)
    state, detail = deadline_window_state(scenario, error, te, doc, line, warnings, refs)
    return detail if state == "not_plausible" else None


def exchange_budget_s(scenario, refs=None):
    """The absolute budget B of one /resolve exchange (app/traversal/backend.py exchange_budget_s): the client timeout
    t (resolve_client_timeout_s) plus the connect, write and pool phase caps, each at most t
    (service_settings.resolve.http_timeouts_s.connect / write / pool). The service takes the exchange's deadline
    (now + B) before it acquires the host's token (exchange_deadline)."""
    t_s = resolve_client_timeout_s(scenario, refs)
    return t_s + sum(min(term(refs, name), t_s) for name in ("connect_timeout_s", "write_timeout_s", "pool_timeout_s"))


def _admin_attempts(doc):
    """The attempts of the last GET /admin/traversal/jobs/{id} the raw record's calls hold (the runner's admin read),
    [] when there is none."""
    for call in reversed((doc or {}).get("calls") or []):
        tool = str(call.get("tool") or "") if isinstance(call, dict) else ""
        if tool.startswith("GET /admin/traversal/jobs/"):
            payload = call.get("payload")
            rows = payload.get("attempts") if isinstance(payload, dict) else None
            return [a for a in rows if isinstance(a, dict)] if isinstance(rows, list) else []
    return []


def _plan_of(doc):
    """overview.plan of the last traverse_results call (the job's plan: timeout_s, max_bound_s), else {}."""
    for call in reversed((doc or {}).get("calls") or []):
        payload = call.get("payload") if isinstance(call, dict) and call.get("tool") == "traverse_results" else None
        plan = ((payload or {}).get("overview") or {}).get("plan") if isinstance(payload, dict) else None
        if isinstance(plan, dict):
            return plan
    return {}


def _evidence_number(record, key):
    return _number_or_none((record or {}).get(key)) if isinstance(record, dict) else None


def deadline_evidence_from_record(doc, refs=None):
    """deadline_evidence as the raw record of a traversal job gives it, None when it gives no figure of it: the
    failing attempt (the last settled one of the runner's GET /admin/traversal/jobs/{id} read) states
    exchange_budget_s (reserved.bound_ms / 1000, else the attempt's bound_ms), timeout_s (reserved.timeout_s),
    granted_at (its dispatched_at: the service claims the attempt at the grant of the host's token), settled_at and
    span_ms = settled_at - granted_at (span_source "attempt"); without the admin read the runner's phase_sources
    (its attempts' last span_ms, "attempt"; else the job's finished_at - started_at, "job": a one-attempt job's span,
    longer with retries) and the job plan (traverse_results overview.plan: max_bound_s, timeout_s);
    acquire_timeout_s is the service's (service_settings.resolve.acquire_timeout_ms / 1000)."""
    out = dict((k, None) for k in DEADLINE_EVIDENCE_KEYS)
    settled = [a for a in _admin_attempts(doc) if a.get("settled_at")]
    sources = (doc or {}).get("phase_sources") if isinstance((doc or {}).get("phase_sources"), dict) else {}
    plan = _plan_of(doc)
    if settled:
        a = settled[-1]
        reserved = a.get("reserved") if isinstance(a.get("reserved"), dict) else {}
        bound_ms = _evidence_number(reserved, "bound_ms")
        bound_ms = bound_ms if bound_ms is not None else _evidence_number(a, "bound_ms")
        out.update(exchange_budget_s=None if bound_ms is None else bound_ms / 1000.0,
                   timeout_s=_evidence_number(reserved, "timeout_s"),
                   granted_at=a.get("dispatched_at"), settled_at=a.get("settled_at"),
                   span_ms=_ms_between(runner_epoch(a.get("dispatched_at")), runner_epoch(a.get("settled_at"))))
        out["span_source"] = "attempt" if out["span_ms"] is not None else None
    if out["span_ms"] is None:
        span = settled_attempt_ms(sources)
        if span is not None:
            out.update(span_ms=span, span_source="attempt")
        else:
            job = _ms_between(runner_epoch(sources.get("started_at")), runner_epoch(sources.get("finished_at")))
            if job is not None:
                out.update(granted_at=sources.get("started_at"), settled_at=sources.get("finished_at"), span_ms=job,
                           span_source="job")
    if out["exchange_budget_s"] is None:
        out["exchange_budget_s"] = _evidence_number(plan, "max_bound_s")
    if out["timeout_s"] is None:
        out["timeout_s"] = _evidence_number(plan, "timeout_s")
    if all(out[k] is None for k in DEADLINE_EVIDENCE_KEYS):
        return None
    out["acquire_timeout_s"] = term(refs, "acquire_timeout_ms") / 1000.0
    return out


def deadline_evidence(line, doc, refs=None, warnings=None):
    """The record deadline_window_verdict reads: the runner's `deadline_evidence` (the raw record's, else the
    line's; DEADLINE_EVIDENCE_KEYS) when it records one, its figures cross-checked with deadline_evidence_from_record
    (a different span_ms or exchange_budget_s is a warning in |warnings|; the recorded one is read), else
    deadline_evidence_from_record; None when neither states anything."""
    computed = deadline_evidence_from_record(doc, refs)
    if isinstance(doc, dict) and "deadline_evidence" in doc:
        rec, where = doc.get("deadline_evidence"), "the raw record's deadline_evidence"
    elif isinstance(line, dict) and "deadline_evidence" in line:
        rec, where = line.get("deadline_evidence"), "the line's deadline_evidence"
    else:
        return computed
    if rec is None:
        return computed
    if not isinstance(rec, dict):
        _warn(warnings, "%s %s is no object: read from the record (%s)" % (
            where, json.dumps(rec, sort_keys=True)[:200], json.dumps(computed, sort_keys=True)))
        return computed
    for key in ("span_ms", "exchange_budget_s"):
        mine, theirs = _evidence_number(computed, key), _evidence_number(rec, key)
        if computed is not None and mine is not None and (theirs is None or abs(mine - theirs) > 0.05):
            _warn(warnings, "%s states %s %s, the record gives %s: the recorded one is read" % (
                where, key, json.dumps(rec.get(key)), mine))
    return dict((k, rec.get(k)) for k in DEADLINE_EVIDENCE_KEYS)


def deadline_window_verdict(scenario, stated, evidence=None, refs=None):
    """(state, detail) of a resolve_timeout stating timeout_phase "deadline": whether the service can have applied the
    window bound_s it states (app/traversal/resolve_exec.py takes the exchange's absolute deadline, now + B
    (exchange_budget_s), before it acquires the host's token (at most A, service_settings.resolve.acquire_timeout_ms)
    and claims the attempt; the window is what is left of B when the exchange starts, so a slow acquire or claim
    legitimately shortens it below the read floor t - 1 s). The one definition async's runner takes
    (shared_rule("deadline_window_verdict")) with its own |evidence|. |stated|: the statement (timeout_terms:
    timeout_phase, bound_s, and where the service states them pre_exchange_ms, from the deadline's start to the
    exchange's, and acquire_ms, the token acquire's own time); |evidence|: deadline_evidence's record
    (exchange_budget_s, span_ms of the failing attempt from its grant to its settling, span_source). In order:
    1. a recorded exchange_budget_s other than B (beyond the 1 ms it is rounded to): not_plausible, "the attempt's
       exchange budget X s is not t + the phase caps B s";
    2. bound_s above B: not_plausible, "service deadline window b s above its exchange budget B s";
    3. pre_exchange_ms stated (preferred: the service places the deadline's start): acquire_ms above it is
       not_plausible ("the stated acquire a ms exceeds the p ms before the exchange it is part of"); bound_s x 1000 +
       pre_exchange_ms must be B x 1000 within DEADLINE_TERMS_TOLERANCE_MS (bound_s rounded to 1 ms, pre_exchange_ms
       to 0.1 ms), else not_plausible: shorter, "service deadline window b s ended earlier than its terms allow:
       exchange budget B s - p ms before the exchange = w s"; longer, "service deadline window b s is longer than its
       terms allow: ...";
    4. bound_s at or above the read floor t - 1 s: no objection;
    5. below it, the span places the deadline's start: the attempt was granted at most max(A, acquire_ms) after it
       and settled after the deadline, so a span below B - max(A, acquire_ms) - 1 s is not_plausible ("the deadline
       began before the acquire: attempt span X ms (span_source) < exchange budget B s - acquire A s - 1 s = Y
       ms"); a span at or above it is no objection (load before the exchange); no span: no objection when step 3
       placed the start, else undecidable ("<UNDECIDABLE_NO_ORIGIN>; deadline window b s below its floor f s").
    No objection is (None, why): rule 4's wall rule follows (the stated bound - 1 s). Another phase or no
    statement: (None, None)."""
    if not isinstance(stated, dict) or stated.get("timeout_phase") != "deadline":
        return None, None
    bound_s = _number_or_none(stated.get("bound_s"))
    if bound_s is None:
        return None, None
    t_s = resolve_client_timeout_s(scenario, refs)
    budget_s = exchange_budget_s(scenario, refs)
    acquire_s = term(refs, "acquire_timeout_ms") / 1000.0
    recorded = _evidence_number(evidence, "exchange_budget_s")
    if recorded is not None and abs(recorded - budget_s) > STATED_BOUND_ROUNDING_S + 1e-9:
        return "not_plausible", "the attempt's exchange budget %g s is not t + the phase caps %g s" % (recorded,
                                                                                                    budget_s)
    if bound_s > budget_s + STATED_BOUND_ROUNDING_S + 1e-9:
        return "not_plausible", "service deadline window %g s above its exchange budget %g s" % (bound_s, budget_s)
    pre = _number_or_none(stated.get("pre_exchange_ms"))
    acq = _number_or_none(stated.get("acquire_ms"))
    if pre is not None:
        if acq is not None and acq > pre + DEADLINE_TERMS_TOLERANCE_MS:
            return "not_plausible", "the stated acquire %g ms exceeds the %g ms before the exchange it is part of" % (
                acq, pre)
        allowed_ms = budget_s * 1000.0 - pre
        terms = "exchange budget %g s - %g ms before the exchange = %g s" % (budget_s, pre, allowed_ms / 1000.0)
        if bound_s * 1000.0 < allowed_ms - DEADLINE_TERMS_TOLERANCE_MS:
            return "not_plausible", "service deadline window %g s ended earlier than its terms allow: %s" % (bound_s,
                                                                                                         terms)
        if bound_s * 1000.0 > allowed_ms + DEADLINE_TERMS_TOLERANCE_MS:
            return "not_plausible", "service deadline window %g s is longer than its terms allow: %s" % (bound_s, terms)
    floor_ms = t_s * 1000.0 - STATED_TIMEOUT_SLACK_MS
    if bound_s * 1000.0 >= floor_ms:
        return None, "deadline window %g s at or above its floor %g s" % (bound_s, floor_ms / 1000.0)
    span = _evidence_number(evidence, "span_ms")
    if span is None:
        if pre is not None:
            return None, "deadline window %g s = %s: shortened before the exchange (load)" % (bound_s, terms)
        return "undecidable", "%s; deadline window %g s below its floor %g s" % (UNDECIDABLE_NO_ORIGIN, bound_s,
                                                                                floor_ms / 1000.0)
    granted_s = max(acquire_s, acq / 1000.0) if acq is not None else acquire_s
    least_ms = (budget_s - granted_s) * 1000.0 - STATED_TIMEOUT_SLACK_MS
    source = (evidence or {}).get("span_source") or "span"
    if span < least_ms:
        return "not_plausible", "the deadline began before the acquire: attempt span %.0f ms (%s) < exchange budget " \
                                "%g s - acquire %g s - 1 s = %.0f ms" % (span, source, budget_s, granted_s, least_ms)
    return None, "deadline window %g s below its floor %g s, shortened before the exchange: attempt span %.0f ms " \
                 "(%s) >= exchange budget %g s - acquire %g s - 1 s = %.0f ms (load)" % (
                     bound_s, floor_ms / 1000.0, span, source, budget_s, granted_s, least_ms)


def deadline_window_state(scenario, error, tool_error, doc=None, line=None, warnings=None, refs=None):
    """deadline_window_verdict of a record: its statement (timeout_terms) and its deadline_evidence."""
    te = tool_error if isinstance(tool_error, dict) else {}
    stated = timeout_terms(line, doc, error, te, warnings)
    if stated is None or stated.get("timeout_phase") != "deadline":
        return None, None
    return deadline_window_verdict(scenario, stated, deadline_evidence(line, doc, refs, warnings), refs)


def tool_error_problem(scenario, tool_error, refs=None):
    """What a /resolve tool error states that the service cannot have applied (the request, the host and the service
    settings decide its bounds): a finalize_reserve_ms above the host's resolve.time_budget.finalize_reserve_ms
    ("stated finalisation reserve <r> ms above the host's <h> ms"), a time_budget_ms other than the one the service
    sent (resolve_budget: the request's, or the service's default; app/traversal/resolve.py plan_budget refuses a
    budget above the host's maximum and never lowers one: "the tool error's time_budget_ms <x> ms is not the
    request's <y> ms"); None otherwise."""
    te = tool_error if isinstance(tool_error, dict) else {}
    stated = _number_or_none(te.get("finalize_reserve_ms"))
    host = term(refs, "resolve_reserve_ms")
    if stated is not None and stated > host:
        return "stated finalisation reserve %g ms above the host's %g ms" % (stated, host)
    sent = _number_or_none(te.get("time_budget_ms"))
    budget, whose = resolve_budget(scenario, refs)
    if sent is not None and budget is not None and sent != budget:
        return "the tool error's time_budget_ms %g ms is not %s %g ms" % (sent, whose, budget)
    return None


def resolve_timeout_least(scenario, error, tool_error, doc=None, line=None, warnings=None, refs=None):
    """(the least plausible wall of a resolve_timeout, the rule as the runner's reason states it), rule 4 (the
    runner's resolve_timeout_least, f94b079d): when the record states what fired (timeout_terms: the runner's
    recorded timeout_terms, else timeout_statement), bound_s x 1000 - 1000 ms, any phase, "the stated bound <b> s
    (<phase>) - 1000 ms"; else (a legacy record) min(the write timeout, max(A, B)) - 1000 ms of whichever of A and B
    is known, "min(the write cap 10000 ms, ...) - 1000 ms" (R2's terms, capped by the write phase, which fires
    first), A = the budget (resolve_budget) + the host's reserve + 10,000 ms, without a budget the service's client
    timeout (legacy_unbudgeted_client_timeout_term), B = the timeout the refusal states ("within <N> s" of the answer
    error's text, else of the tool error's message) x 1000."""
    te = tool_error if isinstance(tool_error, dict) else {}
    stated = timeout_terms(line, doc, error, te, warnings)
    if stated is not None:
        bound = float(stated["bound_s"])
        return bound * 1000.0 - STATED_TIMEOUT_SLACK_MS, "the stated bound %g s (%s) - 1000 ms" % (
            bound, stated["timeout_phase"])
    budget, _ = resolve_budget(scenario, refs)
    reserve = term(refs, "resolve_reserve_ms")
    text = str((error or {}).get("text") or "") if isinstance(error, dict) else ""
    m = STATED_TIMEOUT.search(text) or STATED_TIMEOUT.search(str(te.get("message") or ""))
    terms = []
    if budget is not None:
        allowance = max(term(refs, "allowance_s") * 1000.0, RESOLVE_MIN_ALLOWANCE_MS)
        a = budget + reserve + allowance
        terms.append((a, "time_budget_ms %g + reserve %g + allowance %g = %g ms" % (budget, reserve, allowance, a)))
    else:
        terms.extend(legacy_unbudgeted_client_timeout_term(scenario, refs))
    if m is not None:
        terms.append((float(m.group(1)) * 1000.0, "the stated timeout %s s" % m.group(1)))
    if not terms:
        return None, "no time budget known and no stated timeout"
    top = max(t for t, _ in terms)
    named = terms[0][1] if len(terms) == 1 else "max(%s)" % ", ".join(why for _, why in terms)
    cap = term(refs, "write_timeout_s") * 1000.0
    return min(cap, top) - STATED_TIMEOUT_SLACK_MS, "min(the write cap %g ms, %s) - 1000 ms" % (cap, named)


def legacy_unbudgeted_client_timeout_term(scenario, refs=None):
    """The legacy rule's term for a resolve_timeout with no stated phase and no budget known (a discovery without
    one, or neither labels nor a budget): [(t in ms, "the client timeout <t> ms")], t the service's client timeout
    (resolve_client_timeout_s: the mode's base, 20 s / 120 s), so the threshold is min(the write cap, t) - 1 s =
    9,000 ms instead of the answer's own "within <N> s" - 1 s or rule (c). PENDING the async session's agreement
    (lane G2, after its verification): it narrows rule (c) for such records."""
    t_ms = resolve_client_timeout_s(scenario, refs) * 1000.0
    return [(t_ms, "the client timeout %g ms" % t_ms)]


def deadline_least(budget):
    """(the least plausible service-side wall of a host's 503 deadline at |budget| ms, the rule as stated): the host
    answers its deadline only once the budget has elapsed (Deadline::respond_expired; the finalisation reserve ends
    the work earlier, not the answer), less STATED_TIMEOUT_SLACK_MS of clock and phase slack. One floor for the
    relayed /resolve deadline (resolve_deadline_least) and a deadline-ended pattern task (wall_rule_verdict); the
    direct arm's client wall has its own 1 ms slack (checks._budget_503)."""
    return budget - STATED_TIMEOUT_SLACK_MS, "time_budget_ms - 1000 ms"


def resolve_deadline_least(scenario, tool_error, refs=None):
    """(the least plausible wall of a relayed /resolve deadline, the rule as stated): deadline_least of the budget
    the service sent (resolve_budget: the request's, an explicit-labels resolve's without one the service's default);
    (None, why) without a budget. The tool error's own budget and reserve are checked, not read
    (tool_error_problem)."""
    budget, _ = resolve_budget(scenario, refs)
    if budget is None:
        return None, "no time budget known"
    return deadline_least(budget)


def discovery_without_budget(scenario):
    """Whether a scenario's MCP arguments name `discover` (truthy) and no time_budget_ms (rule 1, 62ee4214; the
    runner's unbudgeted_discovery, which reads the arguments of any tool: only traverse_resolve has `discover`): the
    service sends a /resolve deadline only on request for a discovery (app/traversal/resolve.py plan_budget: a
    discovery gets none unless asked; only an explicit-labels resolve gets the service's default), and
    SPEC-labeled-traversal-core §4.5 makes the host's deadline opt-in, so a relayed deadline there is no host
    state."""
    args = (((scenario or {}).get("mcp") or {}).get("arguments") or {})
    return isinstance(args, dict) and bool(args.get("discover")) and args.get("time_budget_ms") is None


def failing_attempt_wall(scenario, line, doc=None):
    """(the failing attempt's own service-side wall in ms, its source) for an answer-triggered case's wall rule: a
    pattern search's deadline-ended task's one failed HTTP try (deadline_task_ms), a traversal job's last settled
    attempt's span (settled_attempt_ms); else (None, WALL_NONE). Never mcp_phases.backend_ms: a pattern search's
    is the sum over its tasks' calls, a traversal job's the sum over its sent attempts, so a retry would count twice."""
    sources = (doc or {}).get("phase_sources") if isinstance((doc or {}).get("phase_sources"), dict) else {}
    tool = mcp_tool_of(scenario)
    if tool == "pattern_search":
        task = deadline_task_ms(doc)
        if task is not None:
            return task, WALL_DEADLINE_TASK
    if tool in ("traverse_resolve", "traverse_fetch"):
        span = settled_attempt_ms(sources)
        if span is not None:
            return span, WALL_ATTEMPT
    return None, WALL_NONE


def host_state_verdict(scenario, line, doc=None, error=None, tool_error=None, warnings=None, refs=None):
    """(code, state, detail) of an MCP execution answered with a host-state code (MCP_HOST_STATE_CODES) the scenario
    does not expect (the runner's host_state_plausible, its reason the detail): state "host_state" (inconclusive by
    the runner's host-state rule), "undecidable" (rule (c): no threshold, UNDECIDABLE_NO_THRESHOLD, or no
    service-side wall and rule 2 decided nothing, undecidable_timing: UNDECIDABLE_NO_WALL, + RULE2_SKIPPED_API_KEY
    with an API key; a deadline window below the read floor whose origin nothing places, UNDECIDABLE_NO_ORIGIN;
    inconclusive) or "not_plausible" (the checks judge: a known wall below the threshold; rule 2, an anonymous client
    wall + 1 s below it; rule 1, a relayed /resolve deadline on a discovery that sent no budget; a tool error stating
    a reserve or budget the service cannot have applied, tool_error_problem; rule 4's read floor,
    read_bound_below_floor, and the write and deadline bounds, stated_bound_problem / deadline_window_verdict,
    whatever the wall); (None, None, None) for any other answer. traversal_resolve_deadline: discovery_without_budget,
    tool_error_problem, then resolve_deadline_least; resolve_timeout: tool_error_problem, read_bound_below_floor,
    stated_bound_problem, deadline_window_state's undecidable state, then resolve_timeout_least; the wall:
    service_wall, else rule 2 (client_wall_premature);
    uninterruptible_overrun_exceeds_allowance (the R20 lockout answers at once): no wall rule. |refs| (References;
    default_references()) the host's and the service's terms. |warnings| (a list, optional): a recorded client wall
    or statement the record does not reproduce."""
    err = error if isinstance(error, dict) else (line or {}).get("error")
    err = err if isinstance(err, dict) else {}
    code = err.get("code")
    if code is None or str(code) not in MCP_HOST_STATE_CODES:
        return None, None, None
    code = str(code)
    if code in json.dumps((scenario or {}).get("expect") or []):
        return None, None, None   # the scenario expects the code: its checks judge it
    if code not in ("resolve_timeout", "traversal_resolve_deadline"):
        return code, "host_state", ""
    if tool_error is None:
        tool_error = tool_error_of(doc, {"error": err})
    if code == "traversal_resolve_deadline" and discovery_without_budget(scenario):
        return code, "not_plausible", DISCOVERY_NO_DEADLINE
    problem = tool_error_problem(scenario, tool_error, refs)
    if problem is not None:   # the tool error states a bound the service cannot have applied: no host state
        return code, "not_plausible", problem
    if code == "traversal_resolve_deadline":
        least, rule = resolve_deadline_least(scenario, tool_error, refs)
    else:
        for rule4 in (read_bound_below_floor, stated_bound_problem):
            below = rule4(scenario, err, tool_error, doc, line, warnings, refs)
            if below is not None:   # the service's own bound was too short (or not its own), whatever the wall
                return code, "not_plausible", below
        window, why = deadline_window_state(scenario, err, tool_error, doc, line, warnings, refs)
        if window == "undecidable":   # a window below the read floor whose origin nothing places
            return code, "undecidable", why
        least, rule = resolve_timeout_least(scenario, err, tool_error, doc, line, warnings, refs)
    wall, which = service_wall(scenario, line, doc)
    if least is None:
        return code, "undecidable", UNDECIDABLE_NO_THRESHOLD
    if wall is None:   # rule 2: the anonymous client's wall is an upper bound
        premature, note = client_wall_premature(scenario, line, doc, least, warnings)
        if premature:
            return code, "not_plausible", str(note)
        return code, "undecidable", undecidable_timing(note)
    if wall < least:
        return code, "not_plausible", "%s after %.0f ms (%s), before %s = %.0f ms: not a host state" % (
            code, wall, which, rule, least)
    return code, "host_state", "%s after %.0f ms (%s) >= %s = %.0f ms" % (code, wall, which, rule, least)


def wall_rule_verdict(rule, case_id, scenario, answer, line, doc=None, warnings=None, refs=None):
    """(state, detail) of an answer-triggered case's answer_match.wall_rule for an MCP execution (R1; the runner's
    wall_rule_holds): budget_minus_reserve: the failing attempt's own wall (failing_attempt_wall: never a sum over
    a search's tasks or a job's attempts) >= the request's time_budget_ms (request_budget_ms) - 1000 ms (the engine
    answers its 503 deadline only once time_budget_ms has elapsed, Deadline::respond_expired; 1 s of clock and phase
    slack);
    a finalisation reserve the answer's text states above the host's (the host_capabilities' pattern
    .finalize_reserve_ms) is no host deadline: held back. state "holds", "held_back" (a known wall below the
    threshold; rule 2, an anonymous client wall + 1 s below it, "<case> (wall_rule <rule>): client wall X ms + 1 s
    slack < threshold Y ms: provably premature"; a stated reserve above the host's; an unknown rule: the case does
    not apply, matching continues) or "undecidable" (rule (c): no request budget, UNDECIDABLE_NO_THRESHOLD, or no
    service-side wall and rule 2 decided nothing, undecidable_timing: the case matched, the execution is
    inconclusive). |refs| and |warnings| as host_state_verdict's."""
    if rule not in WALL_RULE_VALUES:
        return "held_back", "%s: unknown wall_rule %r: the case does not apply" % (case_id, rule)
    budget = request_budget_ms(scenario)
    if budget is None:
        return "undecidable", UNDECIDABLE_NO_THRESHOLD
    err = (answer or {}).get("error")
    m = STATED_RESERVE.search(str((err or {}).get("text") or "")) if isinstance(err, dict) else None
    host = term(refs, "pattern_reserve_ms")
    if m and float(m.group(1)) > host:
        return "held_back", "%s (wall_rule %s): stated finalisation reserve %s ms above the host's %g ms: the case " \
                            "does not apply" % (case_id, rule, m.group(1), host)
    least, _ = deadline_least(budget)
    terms = "time_budget_ms %g - 1000 ms = %g ms" % (budget, least)
    wall, which = failing_attempt_wall(scenario, line, doc)
    if wall is None:   # rule 2: the anonymous client's wall is an upper bound
        premature, note = client_wall_premature(scenario, line, doc, least, warnings)
        if premature:
            return "held_back", "%s (wall_rule %s): %s" % (case_id, rule, note)
        return "undecidable", undecidable_timing(note)
    if wall < least:
        return "held_back", "%s (wall_rule %s): after %.0f ms (%s), before %s: the case does not apply" % (
            case_id, rule, wall, which, terms)
    return "holds", "%s (wall_rule %s): after %.0f ms (%s) >= %s" % (case_id, rule, wall, which, terms)


def raw_name(scenario_id, arm, repeat):
    safe = re.sub(r"[^A-Za-z0-9._-]", "_", str(scenario_id))
    return "raw/%s.%s.%d.json.gz" % (safe, arm, repeat)


# ---------------------------------------------------------------------------------------------------------------
# answers

def decode_body(raw, content_encoding):
    """The body bytes as sent before the transfer's compression."""
    enc = (content_encoding or "").strip().lower()
    if not raw or enc in ("", "identity"):
        return raw
    if enc == "gzip":
        return gzip.decompress(raw)
    if enc == "deflate":
        try:
            return zlib.decompress(raw)
        except zlib.error:
            return zlib.decompress(raw, -zlib.MAX_WBITS)
    raise ValueError("unsupported content-encoding %r" % content_encoding)


def normalise_direct(http, body_bytes, headers=None):
    """A direct answer in the shape of DESIGN §8.2: {http, body, error, headers}.

    `body` is the parsed JSON (or the text when the body is not JSON, None when empty); `error` is
    {code, text} for a non-2xx answer, from the host's {"code", "error"} body (code None when the body names
    none); `headers` keeps content-type and content-encoding for the protocol checks (left out of the digest).
    """
    text = body_bytes.decode("utf-8", "replace") if body_bytes else ""
    body = None
    if text.strip():
        try:
            body = json.loads(text)
        except ValueError:
            body = text
    error = None
    if not (isinstance(http, int) and 200 <= http < 300):
        if isinstance(body, dict):
            msg = body.get("error")
            if msg is not None and not isinstance(msg, str):
                msg = json.dumps(msg, sort_keys=True)
            error = {"code": body.get("code"), "text": msg}
        else:
            error = {"code": None, "text": text[:4000] if text else None}
    keep = {}
    for name, value in (headers or {}).items():
        lname = name.lower()
        if lname in ("content-type", "content-encoding"):
            keep[lname] = value
    return {"http": http, "body": body, "error": error, "headers": keep}


def transport_failure(code, text):
    """The normalised answer of a request that got no HTTP status (connection refused, timeout, ...)."""
    return {"http": 0, "body": None, "error": {"code": code, "text": text}, "headers": {}}


def _strip(obj, drop_keys, drop_prefixes):
    if isinstance(obj, dict):
        out = {}
        for key, value in obj.items():
            if key in drop_keys:
                continue
            if drop_prefixes and isinstance(key, str) and key.startswith(drop_prefixes):
                continue
            out[key] = _strip(value, drop_keys, drop_prefixes)
        return out
    if isinstance(obj, list):
        return [_strip(v, drop_keys, drop_prefixes) for v in obj]
    return obj


def digest_view(answer, route=None):
    """The part of a normalised answer the digest covers (see DIGEST_DROP_KEYS)."""
    prefixes = ()
    if route is not None:
        _, path = parse_route(route)
        if path.split("?", 1)[0].rstrip("/").endswith("capabilities"):
            prefixes = CAPABILITY_VOLATILE_PREFIXES
    view = {k: v for k, v in answer.items() if k != "headers"}
    return _strip(view, DIGEST_DROP_KEYS, prefixes)


def answer_digest(answer, route=None):
    return hashlib.sha256(canonical_json(digest_view(answer, route))).hexdigest()


def _is_number(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def accumulate_work(total, work, prefix=""):
    """Adds a work object's numeric counters into `total` (nested objects as dotted keys).

    Counters sum across entries (patterns, seeds, arms); `memory_bytes` and `max_*` are peaks, so they take
    the maximum.
    """
    if not isinstance(work, dict):
        return total
    for key, value in work.items():
        name = prefix + str(key)
        if isinstance(value, dict):
            accumulate_work(total, value, name + ".")
        elif _is_number(value):
            peak = str(key) == "memory_bytes" or str(key).startswith("max_")
            if name in total:
                total[name] = max(total[name], value) if peak else total[name] + value
            else:
                total[name] = value
    return total


def extract_server(body):
    """The server's own measurements of an answer: {elapsed_ms, timing, work}, or None without a JSON object.

    /resolve and /pattern carry a top-level `timing`; /pattern and /traverse carry per-entry `timing` and
    `work` (per pattern; per seed with per-arm `counters` and an `annotation` block). elapsed_ms is the
    top-level elapsed time, else the sum of the entries' elapsed times; work sums the entries' counters.
    """
    if not isinstance(body, dict):
        return None
    timing = body.get("timing") if isinstance(body.get("timing"), dict) else None
    work = {}
    entry_elapsed = []
    entry_timing = []
    entries = body.get("patterns")
    if not isinstance(entries, list):
        entries = body.get("results")
    if isinstance(entries, list):
        for entry in entries:
            if not isinstance(entry, dict):
                continue
            t = entry.get("timing")
            if isinstance(t, dict):
                entry_timing.append(t)
                if _is_number(t.get("elapsed_ms")):
                    entry_elapsed.append(t["elapsed_ms"])
            accumulate_work(work, entry.get("work"))
            arms = entry.get("arms")
            if isinstance(arms, dict):
                for arm in arms.values():
                    if isinstance(arm, dict):
                        accumulate_work(work, arm.get("counters"), "counters.")
            if isinstance(entry.get("annotation"), dict):
                accumulate_work(work, entry["annotation"], "annotation.")
    accumulate_work(work, body.get("work"))
    elapsed = None
    if timing is not None and _is_number(timing.get("elapsed_ms")):
        elapsed = timing["elapsed_ms"]
    elif entry_elapsed:
        elapsed = sum(entry_elapsed)
    if timing is None and entry_timing:
        timing = {"entries": entry_timing}
    if elapsed is None and timing is None and not work:
        return None
    return {"elapsed_ms": elapsed, "timing": timing, "work": work or None}


def find_key(obj, key):
    """The first non-null value of `key` anywhere in obj (depth-first, document order)."""
    if isinstance(obj, dict):
        if obj.get(key) is not None:
            return obj[key]
        for value in obj.values():
            found = find_key(value, key)
            if found is not None:
                return found
    elif isinstance(obj, list):
        for value in obj:
            found = find_key(value, key)
            if found is not None:
                return found
    return None


def identity_from_capabilities(docs):
    """release, server_instance and index_fp from the capabilities documents {path: body}."""
    ident = {"release": None, "server_instance": None, "index_fp": None}
    for path in ("/capabilities", "/traverse/capabilities", "/pattern/capabilities"):
        doc = docs.get(path)
        if not isinstance(doc, dict):
            continue
        for key in ident:
            if ident[key] is None:
                ident[key] = find_key(doc, key)
    return ident


# What the MCP cannot show of a scenario at PROPERTY level, by its true class (manifest v5; the owner's decision of
# 2026-10-15, lane V5 adds the fields): mcp.not_implemented_checks [{check, reason}] a missing MCP feature,
# mcp.by_design_checks [{check, case}] a deliberate service behaviour; the rest of the scenario stays tested via MCP.
# (manifest key, mark, the entry's detail key)
MCP_CHECK_MARKS = (("not_implemented_checks", "not_implemented", "reason"), ("by_design_checks", "by_design", "case"))
# how the report lists each mark
MCP_CHECK_MARK_LABELS = (("not_implemented", "Not implemented yet via MCP (property)"),
                         ("by_design", "Service-side by design (property)"))
# The checks that hold on any answer of the route (checks.GENERIC_CHECKS, for a check result that does not state its
# own `property`): every other check is a PROPERTY check, the claim the scenario was built to test
GENERIC_CHECKS = frozenset({"status", "invariants", "counts_consistent", "budget_respected"})
# The checks a repeat evaluates only when an earlier one exists (checks.evaluate: determinism has nothing to compare
# on the first repeat, so it has no result there)
LATER_REPEATS_ONLY = ("deterministic",)
# An MCP execution whose every property check a property-level mark covers tests nothing of its scenario's claim via
# MCP: that is a scenario-level mark (mcp.not_implemented / mcp.by_design), and the execution is inconclusive with
# this inconclusive_reason, never a pass on its generic checks alone (every_property_marked)
SCENARIO_LEVEL_MARK = "every property check is marked; this is a scenario-level mark"


def mcp_check_marks(scenario):
    """{check name: {"mark": "not_implemented" | "by_design", "reason" | "case": its detail}} of a scenario's
    property-level marks (MCP_CHECK_MARKS; the first entry of a name wins); {} for a manifest without them (v1-v4)."""
    mcp = (scenario or {}).get("mcp")
    out = {}
    for key, mark, detail in MCP_CHECK_MARKS:
        entries = mcp.get(key) if isinstance(mcp, dict) else None
        for e in entries if isinstance(entries, list) else []:
            name = e.get("check") if isinstance(e, dict) else None
            if isinstance(name, str) and name and name not in out:
                out[name] = {"mark": mark, detail: e.get(detail)}
    return out


def marked_aside(scenario, arm, item, marks=None):
    """The mark (mcp_check_marks' entry) that sets one check result aside: on the MCP arm, a check the scenario marks
    that reads not_exposed, what the mark declares; None otherwise. A marked check that reads anything else (pass,
    fail, finding, inconclusive) is judged as any check: a mark hides no result. |marks|: mcp_check_marks(scenario)
    when the caller has it. The readers of an execution apply checks_marked_aside, which builds on this."""
    if arm != "mcp" or not isinstance(item, dict) or item.get("outcome") != "not_exposed":
        return None
    marks = mcp_check_marks(scenario) if marks is None else marks
    return marks.get(item.get("check"))


def check_keys(checks):
    """[(check name, its occurrence among the items of that name)] of a check list, in list order: how the
    comparison, the gate and the report address one item of a check that occurs more than once."""
    seen, out = {}, []
    for c in checks if isinstance(checks, list) else []:
        name = c.get("check") if isinstance(c, dict) else None
        out.append((name, seen.get(name, 0)))
        seen[name] = seen.get(name, 0) + 1
    return out


def _property_item(item):
    """Whether a check result is a property check's: its own `property` (checks.evaluate records it), else its name
    outside GENERIC_CHECKS."""
    flag = item.get("property")
    return flag if isinstance(flag, bool) else item.get("check") not in GENERIC_CHECKS


def expected_property_checks(scenario):
    """The names of the property checks |scenario| expects (its `expect` entries outside GENERIC_CHECKS), each once,
    in order."""
    expect = scenario.get("expect") if isinstance(scenario, dict) else None
    out = []
    for e in expect if isinstance(expect, list) else []:
        name = e.get("check") if isinstance(e, dict) else None
        if isinstance(name, str) and name and name not in GENERIC_CHECKS and name not in out:
            out.append(name)
    return out


def every_property_marked(scenario, arm, checks, marks=None):
    """The sorted names of an MCP execution's property checks when the scenario's property-level marks cover every one
    of them (at least one, each recorded reading not_exposed: marked_aside) and every property check the scenario
    expects (expected_property_checks; LATER_REPEATS_ONLY aside, which a first repeat has no result of), else [].
    Such an execution tests nothing of the scenario's claim via MCP (SCENARIO_LEVEL_MARK): checks_marked_aside sets
    nothing aside for it and unmarked_checks reads its marked items inconclusive. An unmarked property check the
    scenario expects and the execution has no result of (a runner that dropped it) is no scenario-level mark: it
    stays the missing result it is. |marks|: mcp_check_marks(scenario) when the caller has it; with |scenario| None
    only the recorded checks are read."""
    if arm != "mcp" or not isinstance(checks, list):
        return []
    marks = mcp_check_marks(scenario) if marks is None else marks
    props = [c for c in checks if isinstance(c, dict) and _property_item(c)]
    if not marks or not props or any(marked_aside(None, arm, c, marks) is None for c in props):
        return []
    if any(name not in marks and name not in LATER_REPEATS_ONLY for name in expected_property_checks(scenario)):
        return []
    return sorted(set(str(c.get("check")) for c in props))


def scenario_level_mark_detail(names):
    """The detail of an execution every_property_marked names |names| of: its verdict_detail, and the detail its
    marked items read inconclusive with."""
    return "%s (%s read not_exposed via MCP)" % (SCENARIO_LEVEL_MARK, ", ".join(names))


def checks_marked_aside(scenario, arm, checks, marks=None):
    """{(check name, occurrence) (check_keys): mark} of the items of an execution's |checks| the scenario's
    property-level marks set aside: the one rule the outcome (unmarked_checks), the comparison, the gate and the
    report apply. An item is set aside only as it is recorded reading not_exposed (marked_aside); a marked check read
    otherwise, or with no result on the execution, is judged as any check (a mark hides no result, and a missing
    result stays a missing one: the runner drops nothing, it judges on unmarked_checks). Nothing is set aside when
    the marks cover every property check of the execution (every_property_marked): a scenario-level mark is no
    property-level one. {} on every arm but MCP and for a scenario without marks. |marks|: mcp_check_marks(scenario)
    when the caller has it."""
    marks = (mcp_check_marks(scenario) if marks is None else marks) if arm == "mcp" else {}
    if not marks or not isinstance(checks, list) or every_property_marked(scenario, arm, checks, marks):
        return {}
    out = {}
    for c, key in zip(checks, check_keys(checks)):
        m = marked_aside(None, arm, c, marks)
        if m is not None:
            out[key] = m
    return out


def unmarked_checks(scenario, arm, checks):
    """|checks| without the items checks_marked_aside sets aside: what the outcome rule (outcome_from_checks,
    checks.scenario_outcome) judges an MCP execution on, so a marked check never makes it inconclusive (not exposed)
    or untested; |checks| itself when nothing is set aside (every arm but MCP, every scenario without marks). When
    the marks cover every property check of the execution (every_property_marked) its marked items stay, each read
    inconclusive with scenario_level_mark_detail: either rule then gives inconclusive (fail or finding where another
    check says so), never a pass on the generic checks alone."""
    marks = mcp_check_marks(scenario) if arm == "mcp" else {}
    names = every_property_marked(scenario, arm, checks, marks)
    if names:
        detail = scenario_level_mark_detail(names)
        return [dict(c, outcome="inconclusive", detail=detail) if marked_aside(None, arm, c, marks) is not None
                else c for c in checks]
    aside = checks_marked_aside(scenario, arm, checks, marks)
    if not aside:
        return checks
    return [c for c, key in zip(checks, check_keys(checks)) if key not in aside]


def outcome_from_checks(checks):
    """An execution's outcome from its checks, the rule the MCP runner applies too: fail when a check failed (or
    reported an outcome outside CHECK_OUTCOMES); finding when a check is a finding and none failed; inconclusive
    when a check is inconclusive, or when no check could test anything (none evaluated, or every one
    not_exposed); pass otherwise. Neither finding nor inconclusive ever counts as pass."""
    outs = [c.get("outcome") for c in checks or []]
    if any(o == "fail" or o not in CHECK_OUTCOMES for o in outs):
        return "fail"
    if "finding" in outs:
        return "finding"
    if "inconclusive" in outs or all(o == "not_exposed" for o in outs):
        return "inconclusive"
    return "pass"


# ---------------------------------------------------------------------------------------------------------------
# the by-design table's explanations of check items (compare, merge_run, reevaluate, the gate and async's runner take
# these: one definition)

# the outcomes an explains_checks entry can explain: a failure, or (a {check, detail_matches} entry only) a field the
# service's path removed from the answer, which the checks read not_exposed
EXPLAINABLE_OUTCOMES = ("fail", "not_exposed")


def named_case_ids(scenario):
    """The by-design case ids a scenario names in the manifest, as async's runner reads them (suite_mcp_runner.py
    ServiceByDesign.named, less the table's scenario_ids): mcp.by_design, the scenario's service_by_design (each an id
    or a list of ids), its `sbd:<id>` tags, in that order, each once; table-free (an id the table lacks is kept, so a
    reader can name it). compare.named_cases adds the table's cases listing the scenario in scenario_ids; the gate
    and its policy read these (SD-030)."""
    scenario = scenario or {}
    ids = []
    for own in (((scenario.get("mcp") or {}).get("by_design")), scenario.get("service_by_design")):
        ids += [own] if isinstance(own, str) else [x for x in own or [] if isinstance(x, str)] \
            if isinstance(own, list) else []
    ids += [t[4:] for t in scenario.get("tags") or [] if isinstance(t, str) and t.startswith("sbd:")]
    return list(dict.fromkeys(ids))


def entry_explains(entry, check):
    """Whether one explains_checks entry names a check item: a name, every item of that check; an object {check,
    detail_matches} (peer review R4), an item of that check whose detail the regex matches completely (re.fullmatch,
    re.S)."""
    name = str(check.get("check"))
    if isinstance(entry, str):
        return entry == name
    if not isinstance(entry, dict) or str(entry.get("check")) != name:
        return False
    try:
        return re.fullmatch(str(entry.get("detail_matches")), str(check.get("detail") or ""), re.S) is not None
    except re.error:
        return False


def item_explained(entry, check):
    """Whether one explains_checks entry explains a check item of an execution (an answered / clamped case's E): an
    item of outcome fail it names (entry_explains: a name or {check, detail_matches}), or an item of outcome
    not_exposed whose detail a {check, detail_matches} entry matches completely, a field the service's path removed
    from the answer (traverse.duplicate_seeds: the merged seed's result, "field not in the answer:
    results.1.duplicate"). A bare check name never explains a not_exposed reading."""
    outcome = str(check.get("outcome"))
    if outcome == "fail":
        return entry_explains(entry, check)
    return outcome == "not_exposed" and isinstance(entry, dict) and entry_explains(entry, check)


def states_limit(entry, check):
    """Whether one explains_checks entry states a not_exposed check item as a limit of the MCP surface (the report's
    stated limits, report.stated_mcp_limit): a not_exposed item the entry names (entry_explains: an object
    {check, detail_matches} only an item whose detail it matches; a bare name every item of its check, the report's
    documented exception to item_explained, under which a bare name explains no not_exposed reading)."""
    return str(check.get("outcome")) == "not_exposed" and entry_explains(entry, check)


def explained_outcome(checks, entries, rule):
    """(the outcome, E) of an execution an answered / clamped case matched (triage-2 B-01; async's runner
    explained_outcome at e334715b): F = the failing checks (outcome fail, or one outside CHECK_OUTCOMES), E = the
    items an entry of |entries| (the matched cases' explains_checks) explains (item_explained), U = F - E. fail when
    U is non-empty; else r = |rule| over the checks without E: finding when r is finding (a by-design mark never
    erases one), service_by_design when E is non-empty, else r (a case explaining nothing leaves the checks'
    outcome). E lists each explained check's name once, in check order."""
    failing = [c for c in checks if str(c.get("outcome")) == "fail" or str(c.get("outcome")) not in CHECK_OUTCOMES]
    e_items = [c for c in checks if any(item_explained(entry, c) for entry in entries)]
    explained = list(dict.fromkeys(str(c.get("check")) for c in e_items))
    if any(all(c is not x for x in e_items) for c in failing):
        return "fail", explained
    rest = [dict(c) for c in checks if all(c is not x for x in e_items)]
    r = str(rule(rest))
    if r == "finding":
        return "finding", explained
    if e_items:
        return "service_by_design", explained
    return r, explained


def explained_items(checks, entries):
    """The indices of an execution's check items |entries| explain (item_explained), in check order: what a line
    records as explained_items (B5b-06: per ITEM, as async's runner's explained_item_indices since f386a3c0)."""
    return [i for i, c in enumerate(checks or []) if any(item_explained(entry, c) for entry in entries or [])]


# ---------------------------------------------------------------------------------------------------------------
# statistics

def percentile(values, q):
    """The q-th percentile (0-100) by linear interpolation between closest ranks; None for no values."""
    vals = sorted(v for v in values if _is_number(v))
    if not vals:
        return None
    if len(vals) == 1:
        return float(vals[0])
    pos = (len(vals) - 1) * (q / 100.0)
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return float(vals[lo])
    return vals[lo] + (vals[hi] - vals[lo]) * (pos - lo)


def summarise(values):
    """{n, p50, p90, max, mean} of the numeric values (None fields when there are none)."""
    vals = [v for v in values if _is_number(v)]
    if not vals:
        return {"n": 0, "p50": None, "p90": None, "max": None, "mean": None}
    return {"n": len(vals), "p50": percentile(vals, 50), "p90": percentile(vals, 90),
            "max": float(max(vals)), "mean": sum(vals) / float(len(vals))}


# ---------------------------------------------------------------------------------------------------------------
# files

def read_jsonl(path):
    """(rows, bad): the parsed lines of a JSON-lines file; unparsable lines (a line torn by a kill) are counted,
    not raised. A missing file is empty."""
    rows, bad = [], 0
    if not os.path.exists(path):
        return rows, bad
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                rows.append(json.loads(line))
            except ValueError:
                bad += 1
    return rows, bad


def append_jsonl(path, obj):
    """Appends one line and syncs it, so a kill loses at most the line being written; after such a kill the
    next line starts on a line of its own (the torn one stays unparsable and is skipped by read_jsonl)."""
    data = json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False) + "\n"
    try:
        with open(path, "rb") as fh:
            fh.seek(-1, os.SEEK_END)
            if fh.read(1) != b"\n":
                data = "\n" + data
    except OSError:
        pass  # missing or empty file
    with open(path, "a", encoding="utf-8") as fh:
        fh.write(data)
        fh.flush()
        os.fsync(fh.fileno())


def write_json_atomic(path, obj):
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=1, sort_keys=True, ensure_ascii=False)
        fh.write("\n")
    os.replace(tmp, path)


def read_json(path, default=None):
    if not os.path.exists(path):
        return default
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


def write_gz_json(path, obj):
    tmp = path + ".tmp"
    with gzip.open(tmp, "wt", encoding="utf-8") as fh:
        json.dump(obj, fh, sort_keys=True, ensure_ascii=False)
    os.replace(tmp, path)


def read_gz_json(path):
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        return json.load(fh)


# ---------------------------------------------------------------------------------------------------------------
# the raw record of an execution (SD-093: one reader for every caller). A runner writes it as raw/<name>.json.gz
# relative to its arm directory (run_direct, async's runner); merge_run rewrites a merged line's path to
# arms/<arm>/raw/... relative to the run directory.

def raw_record_path(line, run_dir=None, arm_dir=None):
    """The file of the raw record an execution line names (its `raw`), the first that exists of: the path itself when
    absolute; with |run_dir|, run_dir/<raw>, then run_dir/arms/<the line's arm>/<raw> (a line of an arm, its path
    relative to the arm's directory); with |arm_dir|, arm_dir/<raw>, then, for a merged line's 'arms/...' path, the
    run directory two levels up. None when the line names none or no candidate exists."""
    raw = line.get("raw") if isinstance(line, dict) else None
    if not isinstance(raw, str) or not raw:
        return None
    if os.path.isabs(raw):
        candidates = [raw]
    else:
        candidates = []
        if run_dir:
            candidates += [os.path.join(run_dir, raw), os.path.join(run_dir, "arms", line.get("arm") or "", raw)]
        if arm_dir:
            candidates.append(os.path.join(arm_dir, raw))
            if raw.startswith("arms/"):
                candidates.append(os.path.join(os.path.dirname(os.path.dirname(arm_dir)), raw))
    return next((p for p in candidates if os.path.exists(p)), None)


def read_raw_record(line, run_dir=None, arm_dir=None, any_json=False, default=None):
    """The raw record an execution line names (raw_record_path), decoded: a dict, else |default| (no `raw`, no file,
    a truncated .gz, a cut JSON, a record that is no object). |any_json|: whatever JSON the file holds, |default| only
    when it cannot be read (merge_run.renormalisation counts every record read)."""
    path = raw_record_path(line, run_dir, arm_dir)
    if path is None:
        return default
    try:
        doc = read_gz_json(path)
    except (OSError, ValueError, EOFError):
        return default
    return doc if any_json or isinstance(doc, dict) else default


def raw_record_answer(line, run_dir=None, arm_dir=None):
    """The normalised answer a raw record keeps (`answer`, a dict; both runners keep it), else None."""
    doc = read_raw_record(line, run_dir, arm_dir)
    answer = doc.get("answer") if isinstance(doc, dict) else None
    return answer if isinstance(answer, dict) else None


def raw_record_body(line, run_dir=None, arm_dir=None):
    """The answer body a raw record keeps (`answer.body`, as it is), else None."""
    answer = raw_record_answer(line, run_dir, arm_dir)
    return answer.get("body") if answer is not None else None


REEVALUATED_TAG = "reevaluated"


def reevaluated_name(name, reevaluated=True):
    """The file of a run's re-evaluated set (reevaluate.py) beside the recorded one: executions.jsonl ->
    executions.reevaluated.jsonl, report.md -> report.reevaluated.md; |name| itself when not reevaluated."""
    if not reevaluated:
        return name
    stem, ext = os.path.splitext(name)
    return "%s.%s%s" % (stem, REEVALUATED_TAG, ext)


# ---------------------------------------------------------------------------------------------------------------
# the run's snapshot of what judges its arms (P6-end-to-end-03)

SNAPSHOT_DIR = "scripts"
SNAPSHOT_FILE = "SNAPSHOT.json"
# What the arms run and judge with: the direct runner and its sampler, the checks (or the stub) and suite_common they
# import, and the by-design table the MCP runner reads beside its --checks-dir.
SNAPSHOT_FILES = ("run_direct.py", "sampler.py", "suite_common.py", "checks.py", "checks_stub.py",
                  "service_by_design.json", "service_by_design_g.json")


class SnapshotChanged(Exception):
    """A snapshot's file that is not the one its SNAPSHOT.json records: it is never executed (snapshot_refusal)."""


# the files of a checks directory load_checks_bound executes
EXECUTED_FILES = ("checks.py", "suite_common.py")


def snapshot_record(directory):
    """(SNAPSHOT.json of |directory|, an object whose `files` is an object {name: sha256}, or None, and why not: ''
    when |directory| holds no SNAPSHOT.json, else the cause: it cannot be read, is empty, is no JSON (torn), or is no
    object stating its files). Never raises."""
    path = os.path.join(directory, SNAPSHOT_FILE)
    if not os.path.lexists(path):
        return None, ""
    try:
        with open(path, "r", encoding="utf-8") as fh:
            text = fh.read()
    except (OSError, ValueError) as e:   # ValueError: bytes that are no UTF-8
        return None, "%s cannot be read (%s)" % (path, getattr(e, "strerror", None) or type(e).__name__)
    if not text.strip():
        return None, "%s is empty" % path
    try:
        doc = json.loads(text)
    except ValueError as e:
        return None, "%s is no JSON (%s)" % (path, e)
    if not isinstance(doc, dict) or not isinstance(doc.get("files"), dict):
        return None, "%s states no files" % path
    return doc, ""


def snapshot_refusal(directory, names=EXECUTED_FILES):
    """Why files |names| of |directory| may not be executed, one sentence per cause: its SNAPSHOT.json cannot be read,
    is empty, is no JSON or states no files (snapshot_record); a file of |names| it records is missing or unreadable;
    a file of |names| it holds has another sha256 than the recorded one (changed, or added after the snapshot was
    taken). [] for a directory without SNAPSHOT.json (a tree, a re-evaluation's checks) and for one whose files are the
    recorded ones; with no |names| only the record is read (a reader of its sha256, the bench's harness). Never raises:
    a refusal names its cause."""
    doc, why = snapshot_record(directory)
    if doc is None:
        return ([("%s: %s not executed" % (why, ", ".join(names))) if names else why] if why else [])
    files = doc["files"]
    out = []
    for name in names:
        path, want = os.path.join(directory, name), files.get(name)
        got = file_sha256(path)
        if got is None and want is not None:
            out.append("%s, which SNAPSHOT.json records (sha256 %s), is %s: not executed" % (
                path, str(want)[:12], "missing" if not os.path.lexists(path) else "unreadable"))
        elif got is not None and got != want:
            out.append("%s has sha256 %s, not the %s SNAPSHOT.json records: not executed" % (
                path, got[:12], str(want)[:12] if want else "none"))
    return out


def load_checks_bound(directory, name):
    """checks.py of |directory| as module |name|, its `import suite_common` resolved to |directory|'s copy when it
    has one (a run's snapshot: the digest and normaliser its arms judged with; a re-evaluation's checks directory),
    never to the caller's: the caller's suite_common module and sys.path are restored afterwards. A snapshot
    snapshot_refusal refuses (its SNAPSHOT.json unreadable, or its checks.py or suite_common.py not the one it records,
    a missing one included) is not executed: SnapshotChanged, naming the cause. Raises what loading raises."""
    import importlib.util
    refused = snapshot_refusal(directory)
    if refused:
        raise SnapshotChanged("; ".join(refused))
    own = sys.modules.get("suite_common")
    saved_path = list(sys.path)
    sys.modules.pop("suite_common", None)
    sys.path.insert(0, directory)
    try:
        spec = importlib.util.spec_from_file_location(name, os.path.join(directory, "checks.py"))
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
    finally:
        sys.path[:] = saved_path
        if own is not None:
            sys.modules["suite_common"] = own
        else:
            sys.modules.pop("suite_common", None)
    return mod


def file_sha256(path):
    """The sha256 of a file's bytes, None when it cannot be read."""
    try:
        with open(path, "rb") as fh:
            return hashlib.sha256(fh.read()).hexdigest()
    except (OSError, TypeError):
        return None


def snapshot_dir(run_dir):
    return os.path.join(run_dir, SNAPSHOT_DIR)


def take_snapshot(source, run_dir, commit=None, sessions_before=0, now=None):
    """Copies SNAPSHOT_FILES (those |source| has) into RUN_DIR/scripts and writes RUN_DIR/scripts/SNAPSHOT.json:
    {created_utc, source, scripts_commit, sessions_before (the orchestrator sessions of the run before it: a resume of
    a run older than its snapshot), files: {name: sha256}}. The arms of the run then run and judge with these copies
    (run_suite.sh), whatever happens to |source| meanwhile."""
    import shutil
    import time
    out = snapshot_dir(run_dir)
    os.makedirs(out, exist_ok=True)
    files = {}
    for name in SNAPSHOT_FILES:
        src = os.path.join(source, name)
        if not os.path.isfile(src):
            continue
        shutil.copy2(src, os.path.join(out, name))
        files[name] = file_sha256(os.path.join(out, name))
    doc = {"created_utc": iso_utc(time.time() if now is None else now), "source": os.path.abspath(source),
           "scripts_commit": commit, "sessions_before": sessions_before, "files": files}
    write_json_atomic(os.path.join(out, SNAPSHOT_FILE), doc)
    return doc


def snapshot_state(run_dir, source=None):
    """(SNAPSHOT.json of the run or None, the names whose snapshot file no longer has its recorded sha256 (a missing
    one included), the names whose file in |source| differs from the snapshot's). A SNAPSHOT.json that cannot be read,
    is empty, is no JSON or states no files (snapshot_record) gives (None, [SNAPSHOT_FILE], []): no snapshot that pins
    anything, and changed. Never raises."""
    doc, why = snapshot_record(snapshot_dir(run_dir))
    if doc is None:
        return None, ([SNAPSHOT_FILE] if why else []), []
    files = doc["files"]
    changed = sorted(n for n, sha in files.items()
                     if file_sha256(os.path.join(snapshot_dir(run_dir), n)) != sha)
    live = []
    if source:
        live = sorted(n for n in SNAPSHOT_FILES
                      if file_sha256(os.path.join(source, n)) != files.get(n))
    return doc, changed, live


def snapshot_pins(path, run_dir):
    """Whether a checks directory a runner recorded is the run's snapshot: RUN_ID/scripts at its end (an absolute path
    of the run directory as it was, which a copy of the run elsewhere keeps)."""
    if not isinstance(path, str) or not path:
        return False
    p = os.path.normpath(path)
    return os.path.basename(p) == SNAPSHOT_DIR and \
        os.path.basename(os.path.dirname(p)) == os.path.basename(os.path.normpath(run_dir))


# ---------------------------------------------------------------------------------------------------------------
# the suite lock on mex (review M4)

def parse_etime(text):
    """ps's etime ([[dd-]hh:]mm:ss) in seconds; None when unparsable."""
    t = (text or "").strip()
    if not t:
        return None
    days = 0
    if "-" in t:
        d, t = t.split("-", 1)
        try:
            days = int(d)
        except ValueError:
            return None
    try:
        parts = [int(x) for x in t.split(":")]
    except ValueError:
        return None
    secs = 0
    for x in parts:
        secs = secs * 60 + x
    return days * 86400 + secs


def host_id(env=None, system=None, scutil=None, node=None):
    """A stable id of this machine for the suite lock (review N6): SUITE_HOST_ID when set; on macOS the
    LocalHostName (`scutil --get LocalHostName`), which a DHCP lease or a VPN reconnect does not change,
    unlike `hostname` (staff-net-vpn-dhcp-<n>...: a lock left by a crashed invocation would read as another
    host's, and another machine could later get the same name); else the node name. The arguments replace
    the environment, platform.system(), the scutil call and the node name (tests)."""
    import platform
    env = os.environ if env is None else env
    if env.get("SUITE_HOST_ID"):
        return env["SUITE_HOST_ID"]
    if (system or platform.system()) == "Darwin":
        name = (scutil or _scutil_local_host_name)()
        if name:
            return name
    return node if node is not None else platform.node()


def _scutil_local_host_name():
    import subprocess
    try:
        out = subprocess.run(["scutil", "--get", "LocalHostName"], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                             universal_newlines=True, timeout=10)
    except (OSError, subprocess.SubprocessError):
        return ""
    return out.stdout.strip() if out.returncode == 0 else ""


def lock_break_command(mex_dir, held):
    """The shell command on mex that removes a stale lock only if its content is still `held` (review N6):
    the compare and the removal run inside a mkdir mutex (LOCK.break), so two invocations that judged the
    same lock stale cannot remove a lock one of them has taken since. Exit 0 removed, 3 another removal in
    progress (or a mutex left by a cut-off one), 4 the lock changed or is gone."""
    lock = "%s/LOCK" % mex_dir
    mutex = "%s/LOCK.break" % mex_dir
    return ('mkdir {m} 2>/dev/null || exit 3; '
            'if [ "$(cat {l} 2>/dev/null)" = {h} ]; then rm -f {l}; rc=0; else rc=4; fi; '
            'rmdir {m}; exit $rc').format(m=shlex.quote(mutex), l=shlex.quote(lock), h=shlex.quote(held))


def lock_verdict(lock, hostname, now, pid_alive, pid_age_s):
    """('stale' | 'live' | 'foreign', reason) for the lock document of another invocation.

    stale: the lock names a process of this host that is gone, or whose pid now belongs to a process younger
    than the lock (the pid was reused); live: that process still runs; foreign: another host's process, which
    cannot be checked from here. `hostname` is this machine's host_id(), which the lock's `host` holds (its
    `hostname`, the network name, is shown only). pid_alive / pid_age_s are the facts of this host (callables
    of the pid)."""
    if not isinstance(lock, dict):
        return "stale", "unreadable lock"
    pid = lock.get("pid")
    if lock.get("host") != hostname:
        net = lock.get("hostname")
        return "foreign", "held by pid %s on %s%s since %s" % (pid, lock.get("host"), (" (%s)" % net) if net else "",
                                                               lock.get("started_utc"))
    if not isinstance(pid, int) or not pid_alive(pid):
        return "stale", "pid %s on %s is gone" % (pid, hostname)
    age = pid_age_s(pid)
    started = lock.get("started_epoch")
    if age is not None and _is_number(started) and now - age > started + 5.0:
        return "stale", "pid %s now belongs to a process younger than the lock" % pid
    return "live", "held by pid %s on %s since %s (run %s)" % (pid, hostname, lock.get("started_utc"),
                                                              lock.get("run_id"))


def _pid_alive(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def _pid_age_s(pid):
    import subprocess
    try:
        out = subprocess.run(["ps", "-o", "etime=", "-p", str(pid)], stdout=subprocess.PIPE,
                             stderr=subprocess.DEVNULL, universal_newlines=True, timeout=10).stdout
    except (OSError, subprocess.SubprocessError):
        return None
    return parse_etime(out)


def _main(argv):
    """Small command-line helpers for run_suite.sh.

    release-short RELEASE        print the run directory's release part
    stamp                        print the run directory's UTC part for now
    manifest-sha FILE [--write]  print the manifest's content_sha256 (--write: write it, and relations_sha256 where
                                 the manifest has relations, into the file)
    beyond-status MANIFEST AREAS LIMIT [IDS]
                                 print the expected checks the stub cannot evaluate ('' for none)
    paired-plan MANIFEST RUN_DIR AREAS LIMIT WARM_EVERY ARMS [IDS]
                                 print the paired loop's work: 'SID L R MC MW WL WR WM WB' per scenario (L, R, MC:
                                 1 to do, 0 done, - not run; MW: warm MCP repeats to do, comma-separated, or -; WL,
                                 WR, WM: the warm-ups to run again right before SID on direct-local, direct-remote
                                 and the MCP arm, WB: on direct-local right before SID's first MCP step left,
                                 comma-separated, or -: paired_plan's rewarm)
    lock-check LOCK_JSON HOST_ID
                                 print 'stale|live|foreign REASON' for another invocation's lock
    host-id                      print this machine's stable id for the lock (host_id())
    mcp-step OUT_DIR MANIFEST IDS REPEAT
                                 print 'OUTCOMES WAIT_S' after an MCP step (mcp_step_result: outcomes comma-
                                 separated, the whole seconds to wait for a cancelled job)
    lock-break-cmd MEX_DIR HELD  print the shell command removing the lock only if it still holds HELD
    clock-offset T0 REMOTE T1    print REMOTE - (T0 + T1) / 2 (how far mex's clock is ahead)
    snapshot RUN_DIR SOURCE COMMIT SESSIONS_BEFORE
                                 take the run's snapshot (take_snapshot) unless it has one; print 'created' or
                                 'kept', then 'changed NAME' / 'live-differs NAME' lines (snapshot_state)
    """
    import time
    if len(argv) >= 2 and argv[0] == "release-short":
        print(release_short(argv[1]))
        return 0
    if argv and argv[0] == "stamp":
        print(run_stamp(time.time()))
        return 0
    if len(argv) >= 2 and argv[0] == "manifest-sha":
        manifest = read_json(argv[1])
        sha = scenarios_sha256(manifest.get("scenarios", []))
        if "--write" in argv[2:]:
            manifest["content_sha256"] = sha
            if "relations" in manifest or "relations_sha256" in manifest:
                manifest["relations_sha256"] = relations_sha256(manifest.get("relations", []))
            write_json_atomic(argv[1], manifest)
        print(sha)
        return 0
    if len(argv) >= 2 and argv[0] == "beyond-status":
        areas = [a for a in (argv[2] if len(argv) > 2 else "").split(",") if a] or None
        limit = int(argv[3]) if len(argv) > 3 and argv[3] else None
        ids = [i for i in (argv[4] if len(argv) > 4 else "").split(",") if i] or None
        selected = select_scenarios(read_json(argv[1]), areas, limit, ids)
        print(",".join(beyond_status([sc_ for _, sc_ in selected])))
        return 0
    if len(argv) >= 7 and argv[0] == "paired-plan":
        areas = [a for a in argv[3].split(",") if a] or None
        limit = int(argv[4]) if argv[4] else None
        arms = [a for a in argv[6].split(",") if a]
        ids = [i for i in (argv[7] if len(argv) > 7 else "").split(",") if i] or None
        flag = lambda v: "-" if v is None else ("1" if v else "0")  # noqa: E731
        for item in paired_plan(read_json(argv[1]), argv[2], areas, limit, int(argv[5]), arms, ids):
            rw = item.get("rewarm") or {}
            print("%s %s %s %s %s %s %s %s %s" % (
                item["scenario"], flag(item["direct-local"]), flag(item["direct-remote"]), flag(item["mcp_cold"]),
                ",".join(str(r) for r in item["mcp_warm"]) or "-", ",".join(rw.get("direct-local") or []) or "-",
                ",".join(rw.get("direct-remote") or []) or "-", ",".join(rw.get("mcp") or []) or "-",
                ",".join(rw.get("before_mcp") or []) or "-"))
        return 0
    if len(argv) >= 3 and argv[0] == "lock-check":
        try:
            lock = json.loads(argv[1])
        except ValueError:
            lock = None
        verdict, reason = lock_verdict(lock, argv[2], time.time(), _pid_alive, _pid_age_s)
        print("%s %s" % (verdict, reason))
        return 0
    if len(argv) >= 5 and argv[0] == "mcp-step":
        ids = [i for i in argv[3].split(",") if i]
        outcomes, wait = mcp_step_result(argv[1], read_json(argv[2]), ids, int(argv[4]), time.time())
        print("%s %d" % (",".join(outcomes) or "-", int(math.ceil(wait))))
        return 0
    if argv and argv[0] == "host-id":
        print(host_id())
        return 0
    if len(argv) >= 3 and argv[0] == "lock-break-cmd":
        print(lock_break_command(argv[1], argv[2]))
        return 0
    if len(argv) >= 5 and argv[0] == "snapshot":
        doc, changed, live = snapshot_state(argv[1], argv[2])
        if doc is None and changed:
            # a snapshot whose record is unreadable is never taken over: the run stops
            sys.stderr.write("the run's snapshot: %s\n" % snapshot_record(snapshot_dir(argv[1]))[1])
            return 2
        if doc is None:
            take_snapshot(argv[2], argv[1], argv[3] or None, int(argv[4] or 0))
            print("created")
            return 0
        print("kept")
        for name in changed:
            print("changed %s" % name)
        for name in live:
            print("live-differs %s" % name)
        return 0
    if len(argv) >= 4 and argv[0] == "clock-offset":
        try:
            t0, remote, t1 = float(argv[1]), float(argv[2]), float(argv[3])
        except ValueError:
            print("")
            return 0
        print("%.3f" % (remote - (t0 + t1) / 2.0))
        return 0
    sys.stderr.write(_main.__doc__ + "\n")
    return 2


if __name__ == "__main__":
    sys.exit(_main(sys.argv[1:]))
