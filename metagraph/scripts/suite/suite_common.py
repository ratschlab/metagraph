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
# the host's cache state per execution (review H1)

HOST_CACHE_WARM_S = 60.0
HOST_CACHE_VALUES = ("cold", "warm", "aged", "not_reached", "unknown")


def reached_host(line):
    """Whether an execution reached the host, and so touched its row cache: a direct request that got an
    answer or timed out (not a refused connection), an MCP execution the host answered (outcome pass, fail,
    finding or inconclusive;
    a refusal at submit or a request the MCP cannot express never reaches it)."""
    if line.get("arm") == "mcp":
        return line.get("outcome") in ("pass", "fail", "finding", "inconclusive")
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


def lock_verdict(lock, hostname, now, pid_alive, pid_age_s):
    """('stale' | 'live' | 'foreign', reason) for the lock document of another invocation.

    stale: the lock names a process of this host that is gone, or whose pid now belongs to a process younger
    than the lock (the pid was reused); live: that process still runs; foreign: another host's process, which
    cannot be checked from here. pid_alive / pid_age_s are the facts of this host (callables of the pid)."""
    if not isinstance(lock, dict):
        return "stale", "unreadable lock"
    pid = lock.get("pid")
    if lock.get("host") != hostname:
        return "foreign", "held by pid %s on %s since %s" % (pid, lock.get("host"), lock.get("started_utc"))
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
    lock-check LOCK_JSON HOSTNAME
                                 print 'stale|live|foreign REASON' for another invocation's lock
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
