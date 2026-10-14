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


def paired_plan(manifest, run_dir, areas=None, limit=None, warm_every=10, arms=("direct-local", "direct-remote"),
                ids=None):
    """The paired loop's work per scenario of the selection (select_scenarios: areas, limit, ids), in manifest order:
    [{scenario, direct-local, direct-remote, mcp_cold, mcp_warm}] where a direct arm's value is True (executions
    left), False (all recorded) or None (the arm does not run it), mcp_cold likewise, and mcp_warm the warm repeat
    indices left. Recorded executions are read from RUN_DIR/arms/<arm>/executions.jsonl."""
    recorded = {}
    for arm in ARMS:
        rows, _ = read_jsonl(os.path.join(run_dir, "arms", arm, "executions.jsonl"))
        recorded[arm] = set(execution_key(r) for r in rows if r.get("arm") == arm)
    out = []
    for index, scenario in select_scenarios(manifest, areas, limit, ids):
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
# 1 s) is no host state whatever the wall (read_bound_below_floor). Where the runner recorded what it judged with
# (raw `client_wall`, line `client_wall_ms`; raw / line `timeout_terms`), the recorded value is read and the
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
UNDECIDABLE = (UNDECIDABLE_NO_WALL, UNDECIDABLE_NO_THRESHOLD)
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
    timeout_phase in RESOLVE_TIMEOUT_PHASES and a positive numeric bound_s, both in that one object; None when none
    does (a legacy record)."""
    for source, err in (("tool error", tool_error), ("answer error", error), ("job error", job_error(doc))):
        if not isinstance(err, dict):
            continue
        phase, bound = err.get("timeout_phase"), _number_or_none(err.get("bound_s"))
        if phase in RESOLVE_TIMEOUT_PHASES and bound is not None and bound > 0:
            return {"timeout_phase": phase, "bound_s": bound, "elapsed_ms": _number_or_none(err.get("elapsed_ms")),
                    "source": source}
    return None


def _valid_statement(rec):
    return rec is None or (isinstance(rec, dict) and rec.get("timeout_phase") in RESOLVE_TIMEOUT_PHASES
                           and _number_or_none(rec.get("bound_s")) is not None and rec["bound_s"] > 0)


def timeout_terms(line, doc, error, tool_error, warnings=None):
    """The statement rule 4 reads: what the runner recorded (f94b079d: the raw record's `timeout_terms`, else the
    line's; null = a legacy record), else timeout_statement of the record; a recorded statement the recomputation
    does not reproduce is kept and |warnings| (a list, optional) says so; a malformed one is replaced by the
    recomputation, with a warning."""
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
    if json.dumps(rec, sort_keys=True) != json.dumps(computed, sort_keys=True):
        _warn(warnings, "%s %s differs from the record's statement %s: the recorded one is read" % (
            where, json.dumps(rec, sort_keys=True), json.dumps(computed, sort_keys=True)))
    return rec


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
    deadline (t plus each phase's cap up to t) less what preceded the exchange (the token's acquire, the claim), at
    least the read floor t - 1 s, else "service deadline window <b> s below its own floor <f> s" (%g seconds). None
    for another phase, no statement, or a bound the service can have applied."""
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
    floor_ms = t_s * 1000.0 - STATED_TIMEOUT_SLACK_MS
    if bound_s * 1000.0 >= floor_ms:
        return None
    return "service deadline window %g s below its own floor %g s" % (bound_s, floor_ms / 1000.0)


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


def resolve_deadline_least(scenario, tool_error, refs=None):
    """(the least plausible wall of a relayed /resolve deadline, the rule as stated): the budget the service sent
    (resolve_budget: the request's, an explicit-labels resolve's without one the service's default) - the host's
    finalisation reserve (resolve.time_budget.finalize_reserve_ms); (None, why) without a budget. The tool error's
    own budget and reserve are checked, not read (tool_error_problem)."""
    budget, _ = resolve_budget(scenario, refs)
    if budget is None:
        return None, "no time budget known"
    return budget - term(refs, "resolve_reserve_ms"), "time_budget_ms - finalize_reserve_ms"


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
    with an API key; inconclusive) or "not_plausible" (the checks judge: a known wall below the threshold; rule 2,
    an anonymous client wall + 1 s below it; rule 1, a relayed /resolve deadline on a discovery that sent no
    budget; a tool error stating a reserve or budget the service cannot have applied, tool_error_problem; rule 4's
    read floor, read_bound_below_floor, and the write and deadline bounds, stated_bound_problem, whatever the wall);
    (None, None, None) for any other answer. traversal_resolve_deadline: discovery_without_budget, tool_error_problem,
    then resolve_deadline_least; resolve_timeout: tool_error_problem, read_bound_below_floor, stated_bound_problem,
    then resolve_timeout_least; the wall: service_wall, else rule 2 (client_wall_premature);
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
    least = budget - STATED_TIMEOUT_SLACK_MS
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
    manifest-sha FILE [--write]  print the manifest's content_sha256 (--write: write it, and relations_sha256 where
                                 the manifest has relations, into the file)
    beyond-status MANIFEST AREAS LIMIT [IDS]
                                 print the expected checks the stub cannot evaluate ('' for none)
    paired-plan MANIFEST RUN_DIR AREAS LIMIT WARM_EVERY ARMS [IDS]
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
