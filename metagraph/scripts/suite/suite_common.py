#!/usr/bin/env python3
"""Shared helpers of the suite's direct runner, sampler, comparator and report.

Stdlib only: this module is copied to mex and imported by run_direct.py there (python3.12).
Everything here is pure except the small file helpers at the end, so the unit tests cover it offline.

Contract: DESIGN.md (v2) section 8 of ~/.cache/metagraph-traverse-plans/2026-10-13/suite/.
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


def verify_manifest(manifest):
    """(ok, computed): whether content_sha256 matches the canonical JSON of `scenarios`."""
    computed = scenarios_sha256(manifest.get("scenarios", []))
    return computed == manifest.get("content_sha256"), computed


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
    whatever subset runs.
    """
    area_set = set(a.strip() for a in areas) if areas else None
    id_set = set(ids) if ids else None
    out = []
    for index, scenario in enumerate(manifest.get("scenarios", [])):
        if area_set is not None and scenario.get("area") not in area_set:
            continue
        if id_set is not None and scenario.get("id") not in id_set:
            continue
        out.append((index, scenario))
        if limit is not None and limit > 0 and len(out) >= limit:
            break
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


def beyond_status(scenarios):
    """The expectation names of the scenarios that the checks stub cannot evaluate (anything but `status`)."""
    names = set()
    for scenario in scenarios:
        for exp in scenario.get("expect") or []:
            if isinstance(exp, dict) and exp.get("check") != "status":
                names.add(str(exp.get("check")))
    return sorted(names)


def mcp_repeats(index, scenario, warm_every=10):
    """The MCP executions of a scenario in paired mode, as repeat indices: repeat 0 before the direct arms
    (cold), the others after them (warm); one more than the direct arms' repeats when those run once, so
    every scenario has an MCP warm execution. A scenario the manifest marks not implemented via MCP gets
    repeat 0 only (the MCP runner records the mark); the long class none (the front cuts every arm behind it)."""
    if scenario.get("timing_class") == "long":
        return []
    if ((scenario.get("mcp") or {}).get("not_implemented")):
        return [0]
    return list(range(max(2, repeat_count(index, scenario, warm_every))))


def paired_plan(manifest, run_dir, areas=None, limit=None, warm_every=10, arms=("direct-local", "direct-remote")):
    """The paired loop's work per scenario, in manifest order: [{scenario, direct-local, direct-remote,
    mcp_cold, mcp_warm}] where a direct arm's value is True (executions left), False (all recorded) or None
    (the arm does not run it), mcp_cold likewise, and mcp_warm the warm repeat indices left. Recorded
    executions are read from RUN_DIR/arms/<arm>/executions.jsonl."""
    recorded = {}
    for arm in ARMS:
        rows, _ = read_jsonl(os.path.join(run_dir, "arms", arm, "executions.jsonl"))
        recorded[arm] = set(execution_key(r) for r in rows if r.get("arm") == arm)
    out = []
    for index, scenario in select_scenarios(manifest, areas, limit):
        sid = scenario.get("id")
        item = {"scenario": sid, "mcp_cold": None, "mcp_warm": []}
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
    return out


# ---------------------------------------------------------------------------------------------------------------
# after an MCP step of the paired loop (async's answers to the verification of lane S2b)

# How long a request of a timing class may take (the direct runner's client timeouts).
TIMING_CLASS_TIMEOUT_S = {"fast": 120.0, "slow": 600.0, "long": 1800.0}
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


def request_bound_s(scenario):
    """How long the host may work on the scenario's request: its time_budget_ms (the largest stated in the MCP
    arguments, else in the direct request), else the /pattern default budget, else the timing class's
    timeout (a traversal without a budget has no deadline of its own)."""
    stated = list(_budget_values((scenario.get("mcp") or {}).get("arguments"))) or \
        list(_budget_values(scenario.get("request")))
    if stated:
        return max(stated) / 1000.0
    if "/pattern" in str(scenario.get("route") or ""):
        return PATTERN_DEFAULT_BUDGET_S
    return TIMING_CLASS_TIMEOUT_S.get(scenario.get("timing_class") or "fast", TIMING_CLASS_TIMEOUT_S["fast"])


def mcp_step_result(out_dir, manifest, ids, repeat, now):
    """What an MCP step (one runner invocation for `ids` and `repeat`) left: (outcomes, wait_s). outcomes: the
    recorded outcome per id, `missing` when the runner wrote none (a runner failure). wait_s: how long to wait
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
        outcomes.append(str(line.get("outcome")))
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


def reached_host(line):
    """Whether an execution reached the host, and so touched its row cache: a direct request that got an
    answer or timed out (not a refused connection); an MCP execution only on evidence of a backend call, i.e.
    the host's own timing (`server.elapsed_ms`) or the service's backend phase (`mcp_phases.backend_ms`) is a
    number. Never from the MCP outcome (review N2): a refusal at submit can read `fail` without reaching the
    host, and an `error` or a service_by_design mark can follow a backend call (a relayed refusal does)."""
    if line.get("arm") == "mcp":
        server = line.get("server") if isinstance(line.get("server"), dict) else {}
        phases = line.get("mcp_phases") if isinstance(line.get("mcp_phases"), dict) else {}
        return _is_number(server.get("elapsed_ms")) or _is_number(phases.get("backend_ms"))
    err = line.get("error") or {}
    return not (not line.get("http") and err.get("code") == "transport")


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


def assign_host_cache(rows, clock_offsets=None, warm_s=HOST_CACHE_WARM_S):
    """Adds `host_cache`, `since_prev_s` and `host_cache_after` to every execution line (in place), from the
    order in which the run's executions of each scenario reached the host, over all arms:
    cold = the first execution of the scenario on the host in this run; warm = another execution of it ended
    at most warm_s seconds before this one started; aged = longer ago (the cache state is unknown);
    not_reached = the execution never reached the host; unknown = no start time. since_prev_s is the time from
    the end of the scenario's previous host execution to this start, host_cache_after that execution's arm.
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
        if reached and start is not None:
            stop = end if end is not None else start
            last[sid] = (max(stop, prev[0]) if prev else stop, r.get("arm"))
    return rows


def _number_or_none(v):
    return float(v) if isinstance(v, (int, float)) and not isinstance(v, bool) else None


def scenario_budget_ms(scenario):
    """The time budget a scenario's request states (bounds.time_budget_ms, else time_budget_ms), or None."""
    req = (scenario or {}).get("request")
    if not isinstance(req, dict):
        return None
    bounds = req.get("bounds") if isinstance(req.get("bounds"), dict) else {}
    budget = _number_or_none(bounds.get("time_budget_ms"))
    return budget if budget is not None else _number_or_none(req.get("time_budget_ms"))


# ---------------------------------------------------------------------------------------------------------------
# The service-side walls and the plausibility rules of async's MCP runner (suite_mcp_runner.py at 1b33afe8: peer
# review R1-R3, host_state_plausible, wall_rule_holds, service_wall, resolve_timeout_least), mirrored for the offline
# re-evaluation and the comparator, with the design lead's rule (c): where no threshold or no service-side wall is
# known the plausibility is undecidable, and the execution is inconclusive for that reason (not a host state, not by
# design, not failed on that ground). Each threshold is one function, so a changed bound changes one place.

# /resolve's finalisation reserve when the refusal states none (app/traversal/resolve.py DEFAULT_FINALIZE_RESERVE_MS)
RESOLVE_RESERVE_MS = 250.0
# the least allowance the service's client timeout of a budgeted /resolve keeps over budget + reserve
# (app/traversal/resolve.py client_timeout_s, MIN_BUDGET_ALLOWANCE_S = 10 s)
RESOLVE_MIN_ALLOWANCE_MS = 10000.0
# the timeout a resolve_timeout states it applied: the first "within <N> s" of its text
STATED_TIMEOUT = re.compile(r"\bwithin ([0-9]+(?:\.[0-9]+)?(?:e[+-]?[0-9]+)?) s\b")
# the slack a resolve_timeout's least wall gets (clock and phase differences between the service and the client)
STATED_TIMEOUT_SLACK_MS = 1000.0
# the host's pattern finalisation reserve when the relayed text states none (SPEC-pattern §10.2: 250 in version 1)
PATTERN_RESERVE_MS = 250.0
# the reserve the host's deadline text states ("... (60000 ms, the finalisation reserve of 250 ms included)")
STATED_RESERVE = re.compile(r"finali[sz]ation reserve of ([0-9]+(?:\.[0-9]+)?) ms")
# a pattern task's held-slot line (the service's app/queue/held_span.py LINE_FORMAT)
SLOT_LINE = re.compile(r"^\u29d7 Slot held for .+? ([0-9.]+)s \(att ([^,]+), (\w+), acq ([0-9.]+), rel ([0-9.]+), "
                       r"http (\d+) fail (\d+) ([0-9.]+)s backoff ([0-9.]+)s, (\w+)\)")
WALL_BACKEND = "mcp_phases.backend_ms"
WALL_DEADLINE_TASK = "the deadline-ended task's one failed HTTP try (diagnostic log 'Slot held ... http 1 fail 1')"
WALL_ATTEMPT = "the settled attempt's settled_at - dispatched_at (GET /admin/traversal/jobs/{id} attempts[])"
WALL_JOB = "the job's finished_at - started_at (traverse_status)"
WALL_NONE = "no service-side timing"
# rule (c): the verdict_detail of an execution whose plausibility cannot be decided
UNDECIDABLE_NO_WALL = "no service-side timing: plausibility undecidable"
UNDECIDABLE_NO_THRESHOLD = "no threshold known: plausibility undecidable"
UNDECIDABLE = (UNDECIDABLE_NO_WALL, UNDECIDABLE_NO_THRESHOLD)
WALL_RULE_VALUES = ("budget_minus_reserve",)


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


def service_wall(scenario, line, doc=None):
    """(the service-side wall of an MCP execution in ms, its source) for a host-state or wall rule (R3, R1), in
    order: the line's mcp_phases.backend_ms; for a pattern search the deadline-ended task's HTTP try
    (deadline_task_ms of the raw record); for a traversal job the last settled attempt's span, then the job's
    finished_at - started_at (the raw record's phase_sources); else (None, WALL_NONE). Never the client's wall."""
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
        t0, t1 = parse_utc(sources.get("started_at")), parse_utc(sources.get("finished_at"))
        if t0 is not None and t1 is not None:
            return (t1 - t0) * 1000.0, WALL_JOB
    return None, WALL_NONE


def resolve_timeout_least(scenario, error, tool_error):
    """(the least plausible wall of a resolve_timeout, the rule as stated), R2: max(A, B) - 1000 ms of whichever is
    known, A = time_budget_ms (the MCP arguments', else the tool error's) + reserve (the tool error's, else 250 ms)
    + 10,000 ms, B = the timeout the refusal states ("within <N> s" of its text, else of the tool error's message)
    x 1000; (None, why) with neither."""
    te = tool_error if isinstance(tool_error, dict) else {}
    asked = mcp_budget_ms(scenario)
    budget = asked if asked is not None else _number_or_none(te.get("time_budget_ms"))
    reserve = _number_or_none(te.get("finalize_reserve_ms"))
    reserve = RESOLVE_RESERVE_MS if reserve is None else reserve
    text = str((error or {}).get("text") or "") if isinstance(error, dict) else ""
    m = STATED_TIMEOUT.search(text) or STATED_TIMEOUT.search(str(te.get("message") or ""))
    terms = []
    if budget is not None:
        a = budget + reserve + RESOLVE_MIN_ALLOWANCE_MS
        terms.append((a, "time_budget_ms %g + reserve %g + allowance %g = %g ms" % (
            budget, reserve, RESOLVE_MIN_ALLOWANCE_MS, a)))
    if m is not None:
        terms.append((float(m.group(1)) * 1000.0, "the stated timeout %s s" % m.group(1)))
    if not terms:
        return None, "no time budget known and no stated timeout"
    top = max(t for t, _ in terms)
    named = terms[0][1] if len(terms) == 1 else "max(%s)" % ", ".join(why for _, why in terms)
    return top - STATED_TIMEOUT_SLACK_MS, "%s - 1000 ms" % named


def resolve_deadline_least(scenario, tool_error):
    """(the least plausible wall of a relayed /resolve deadline, the rule as stated): time_budget_ms (the tool
    error's, else the MCP arguments') - finalize_reserve_ms (the tool error's, else 250 ms); (None, why) without a
    budget."""
    te = tool_error if isinstance(tool_error, dict) else {}
    budget = _number_or_none(te.get("time_budget_ms"))
    budget = mcp_budget_ms(scenario) if budget is None else budget
    if budget is None:
        return None, "no time budget known"
    reserve = _number_or_none(te.get("finalize_reserve_ms"))
    return budget - (RESOLVE_RESERVE_MS if reserve is None else reserve), "time_budget_ms - finalize_reserve_ms"


def host_state_verdict(scenario, line, doc=None, error=None, tool_error=None):
    """(code, state, detail) of an MCP execution answered with a host-state code (MCP_HOST_STATE_CODES) the scenario
    does not expect: state "host_state" (inconclusive by the runner's host-state rule), "undecidable" (rule (c): no
    threshold or no service-side wall known; inconclusive, detail UNDECIDABLE_NO_THRESHOLD / UNDECIDABLE_NO_WALL) or
    "not_plausible" (a known wall below the threshold: the checks judge); (None, None, None) for any other answer.
    resolve_timeout: resolve_timeout_least; traversal_resolve_deadline: resolve_deadline_least; the wall:
    service_wall; uninterruptible_overrun_exceeds_allowance (the R20 lockout answers at once): no wall rule."""
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
    if code == "traversal_resolve_deadline":
        least, rule = resolve_deadline_least(scenario, tool_error)
    else:
        least, rule = resolve_timeout_least(scenario, err, tool_error)
    if least is None:
        return code, "undecidable", UNDECIDABLE_NO_THRESHOLD
    wall, which = service_wall(scenario, line, doc)
    if wall is None:
        return code, "undecidable", UNDECIDABLE_NO_WALL
    if wall < least:
        return code, "not_plausible", "%s after %.0f ms (%s), before %s = %.0f ms: not a host state" % (
            code, wall, which, rule, least)
    return code, "host_state", "%s after %.0f ms (%s) >= %s = %.0f ms" % (code, wall, which, rule, least)


def wall_rule_verdict(rule, case_id, scenario, answer, line, doc=None):
    """(state, detail) of an answer-triggered case's answer_match.wall_rule for an MCP execution (R1; async's
    wall_rule_holds at e334715b): budget_minus_reserve: the service-side wall (service_wall) >= the request's
    time_budget_ms (request_budget_ms) - the finalisation reserve the answer's error text states, else 250 ms. state
    "holds", "held_back" (a known wall below the threshold, or an unknown rule: the case does not apply, matching
    continues) or "undecidable" (rule (c): no request budget, UNDECIDABLE_NO_THRESHOLD, or no service-side wall,
    UNDECIDABLE_NO_WALL: the case matched, the execution is inconclusive)."""
    if rule not in WALL_RULE_VALUES:
        return "held_back", "%s: unknown wall_rule %r: the case does not apply" % (case_id, rule)
    budget = request_budget_ms(scenario)
    if budget is None:
        return "undecidable", UNDECIDABLE_NO_THRESHOLD
    err = (answer or {}).get("error")
    m = STATED_RESERVE.search(str((err or {}).get("text") or "")) if isinstance(err, dict) else None
    reserve = float(m.group(1)) if m else PATTERN_RESERVE_MS
    least = budget - reserve
    terms = "time_budget_ms %g - reserve %g (%s) = %g ms" % (budget, reserve, "stated" if m else "default", least)
    wall, which = service_wall(scenario, line, doc)
    if wall is None:
        return "undecidable", UNDECIDABLE_NO_WALL
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


REEVALUATED_TAG = "reevaluated"


def reevaluated_name(name, reevaluated=True):
    """The file of a run's re-evaluated set (reevaluate.py) beside the recorded one: executions.jsonl ->
    executions.reevaluated.jsonl, report.md -> report.reevaluated.md; |name| itself when not reevaluated."""
    if not reevaluated:
        return name
    stem, ext = os.path.splitext(name)
    return "%s.%s%s" % (stem, REEVALUATED_TAG, ext)


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
    manifest-sha FILE [--write]  print (or write into the file) the manifest's content_sha256
    beyond-status MANIFEST AREAS LIMIT
                                 print the expected checks the stub cannot evaluate ('' for none)
    paired-plan MANIFEST RUN_DIR AREAS LIMIT WARM_EVERY ARMS
                                 print the paired loop's work: 'SID L R MC MW' per scenario (L, R, MC: 1 to
                                 do, 0 done, - not run; MW: warm MCP repeats to do, comma-separated, or -)
    lock-check LOCK_JSON HOST_ID
                                 print 'stale|live|foreign REASON' for another invocation's lock
    host-id                      print this machine's stable id for the lock (host_id())
    mcp-step OUT_DIR MANIFEST IDS REPEAT
                                 print 'OUTCOMES WAIT_S' after an MCP step (mcp_step_result: outcomes comma-
                                 separated, the whole seconds to wait for a cancelled job)
    lock-break-cmd MEX_DIR HELD  print the shell command removing the lock only if it still holds HELD
    clock-offset T0 REMOTE T1    print REMOTE - (T0 + T1) / 2 (how far mex's clock is ahead)
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
            write_json_atomic(argv[1], manifest)
        print(sha)
        return 0
    if len(argv) >= 2 and argv[0] == "beyond-status":
        areas = [a for a in (argv[2] if len(argv) > 2 else "").split(",") if a] or None
        limit = int(argv[3]) if len(argv) > 3 and argv[3] else None
        selected = select_scenarios(read_json(argv[1]), areas, limit)
        print(",".join(beyond_status([sc_ for _, sc_ in selected])))
        return 0
    if len(argv) >= 7 and argv[0] == "paired-plan":
        areas = [a for a in argv[3].split(",") if a] or None
        limit = int(argv[4]) if argv[4] else None
        arms = [a for a in argv[6].split(",") if a]
        flag = lambda v: "-" if v is None else ("1" if v else "0")  # noqa: E731
        for item in paired_plan(read_json(argv[1]), argv[2], areas, limit, int(argv[5]), arms):
            print("%s %s %s %s %s" % (item["scenario"], flag(item["direct-local"]), flag(item["direct-remote"]),
                                      flag(item["mcp_cold"]),
                                      ",".join(str(r) for r in item["mcp_warm"]) or "-"))
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
