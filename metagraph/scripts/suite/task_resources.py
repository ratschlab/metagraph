#!/usr/bin/env python3
"""Per-task resource figures of a suite run: what each request cost the engine in CPU time, page faults and memory
(manifest v8; the owner, 2026-10-21: memory per task is the second critical cost next to CPU time). Report-only:
nothing here enters a digest, an identity, a check or the gate (cpu_ms and the page faults vary run to run;
memory_bytes is deterministic but stays report-only in v8).

What the engine states per request on POST /traverse, /resolve and /pattern (engine branch gr/task-resources;
SPEC-labeled-traversal-core section 7.2.1, SPEC-pattern-search section 8.4), top level and per seed / pattern entry:
  timing.cpu_ms        the CPU time of the thread that runs the request's work (the route handler's thread; a request
                       runs on one thread), from the start of the route's work to its timing block: not the body's
                       parse, the index lease or load, the answer's compression or the socket write.
  timing.page_faults   {minor, major}: the same thread's getrusage(RUSAGE_THREAD) deltas over the same span (Linux
                       only; absent elsewhere, never 0). A major fault read a page from disk: an index page the page
                       cache did not hold, the cost of an index far larger than the container's memory limit.
  work.memory_bytes    the peak of the request's memory account: what the engine models for the request (walk state,
                       results held until the answer is written, decode batches, /resolve's query and phases,
                       /pattern's retrieval account); not the process RSS, not the label / row caches, not the index
                       pages in the page cache (those are the major faults' business).

Per execution (execution_figures): each figure from the line's `server` block (suite_common.extract_server: the
top-level timing, else the sum over the entries; work's memory_bytes the maximum over the request and its entries),
else from the raw answer's body (the line's `raw` record; the MCP runner keeps some of the engine's figures only
there, e.g. /resolve's work); `sources` names where each came from: server, server_entries (the sum over the
entries a server block keeps), raw, raw_entries, or null (not stated: an older engine, a refusal, an answer without
timing).

host_memory.json (optional, in the run directory; written by hand at the run's start and end):
    {"before": SNAP, "after": SNAP}
each SNAP a flat object of the engine container's cgroup v2 files and the host's MemAvailable, e.g.
    {"label": "before", "utc": "2026-10-21T06:02:40Z", "container": "metagraph_refseq33m_experimental",
     "page_bytes": 4096, "current": 108165853184, "peak": 137438957568, "max": 137438953472,
     "stat.anon": 9465450496, "stat.file": 96841510912, "stat.pgmajfault": 233625900,
     "stat.pgscan_direct": 246266849, "stat.pgscan_kswapd": 125379262, "stat.pgsteal_direct": 220276454,
     "stat.pgsteal_kswapd": 106915779, "stat.workingset_refault_file": 63140236,
     "events.max": 83142, "events.oom": 0, "events.oom_kill": 0, "events.high": 0, "host_available_kib": 731838912}
current, peak, max: memory.current, memory.peak, memory.max in bytes (null or "max": no limit); stat.*: memory.stat
(anon and file in bytes, the pg* and workingset_* counters in pages or events); events.*: memory.events. Levels
(current, peak, stat.anon, stat.file, host_available_kib) are reported per snapshot; peak is the high-water mark
since the container started, not the run's. Counters (stat.pg*, stat.workingset_*, events.*) are cumulative since
the container started: the report states after minus before. A page is page_bytes (4096 when absent).
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import suite_common as sc  # noqa: E402

FORMAT = 1
TASK_ROUTES = ("POST /traverse", "POST /resolve", "POST /pattern")
# the figures per execution, and the families their sources are stated for
FIGURES = ("cpu_ms", "elapsed_ms", "cpu_ratio", "minor_faults", "major_faults", "memory_bytes")
SOURCE_FAMILIES = ("cpu_ms", "page_faults", "memory_bytes")
SOURCES = ("server", "server_entries", "raw", "raw_entries")
# the lists of the heaviest tests: figure, what ranks them
TOP_FIGURES = ("cpu_ms", "memory_bytes", "major_faults")
TOP_N = 20
# the host-cache labels in reading order (compare.timing_cache's; warm_after_interrupted counted apart, as C.7 does)
STATE_ORDER = ("cold", "warm", "warm_after_interrupted", "aged", "not_reached", "unknown")
HOST_MEMORY_FILE = "host_memory.json"
PAGE_BYTES = 4096
# host_memory.json: the levels (per snapshot) and the counters (after minus before), in reading order
HOST_LEVELS = ("max", "current", "peak", "stat.anon", "stat.file", "host_available_kib")
HOST_COUNTERS = ("events.max", "events.high", "events.oom", "events.oom_kill", "stat.pgsteal_direct",
                 "stat.pgsteal_kswapd", "stat.pgscan_direct", "stat.pgscan_kswapd", "stat.workingset_refault_file",
                 "stat.pgmajfault")
# the counters counted in pages (stated in GiB too)
HOST_PAGE_COUNTERS = ("stat.pgsteal_direct", "stat.pgsteal_kswapd", "stat.pgscan_direct", "stat.pgscan_kswapd",
                      "stat.workingset_refault_file")
# the graphlet tools' ways whose answers a service writes (W3: the service's graphlet tools); a figure is read only
# where such an answer states one
GRAPHLET_SERVICE_WAYS = ("W3",)
GRAPHLET_FIGURE_KEYS = ("cpu_ms", "usage")


def _num(v):
    return sc.finite_or_none(v)


# ---------------------------------------------------------------------------------------------------------------
# per execution

def _entry_list(obj, *keys):
    """The entry objects under the first of |keys| obj holds as a list (or, for by_pattern, an object by id)."""
    if not isinstance(obj, dict):
        return []
    for key in keys:
        v = obj.get(key)
        if isinstance(v, list):
            return [e for e in v if isinstance(e, dict)]
        if isinstance(v, dict):
            return [e for _, e in sorted(v.items()) if isinstance(e, dict)]
    return []


def _timing_figures(timing, entries, scope):
    """{cpu_ms, elapsed_ms, minor_faults, major_faults} and their sources from a request-level timing object and the
    entries' timing objects: the request's own figure where it states one, else the sum over the entries that state
    it (scope + '_entries')."""
    out = {"cpu_ms": None, "elapsed_ms": None, "minor_faults": None, "major_faults": None,
           "sources": {"cpu_ms": None, "page_faults": None}}
    timing = timing if isinstance(timing, dict) else {}
    if _num(timing.get("cpu_ms")) is not None:
        out["cpu_ms"] = float(timing["cpu_ms"])
        out["elapsed_ms"] = _num(timing.get("elapsed_ms"))
        out["sources"]["cpu_ms"] = scope
    else:
        with_cpu = [e for e in entries if _num(e.get("cpu_ms")) is not None]
        if with_cpu:
            out["cpu_ms"] = float(sum(e["cpu_ms"] for e in with_cpu))
            el = [_num(e.get("elapsed_ms")) for e in with_cpu]
            out["elapsed_ms"] = float(sum(el)) if all(x is not None for x in el) else None
            out["sources"]["cpu_ms"] = scope + "_entries"
    pf = timing.get("page_faults")
    if isinstance(pf, dict) and (_num(pf.get("major")) is not None or _num(pf.get("minor")) is not None):
        out["minor_faults"], out["major_faults"] = _num(pf.get("minor")), _num(pf.get("major"))
        out["sources"]["page_faults"] = scope
    else:
        pfs = [e.get("page_faults") for e in entries if isinstance(e.get("page_faults"), dict)]
        if pfs:
            for key, name in (("minor", "minor_faults"), ("major", "major_faults")):
                vals = [_num(p.get(key)) for p in pfs if _num(p.get(key)) is not None]
                out[name] = float(sum(vals)) if vals else None
            out["sources"]["page_faults"] = scope + "_entries"
    return out


def _peak_memory(work):
    """The largest memory_bytes a work object states: its own key, a flattened key ending in .memory_bytes
    (extract_server's dotted names), or one in a nested object or list (per-entry blocks such as by_pattern)."""
    best = None
    if isinstance(work, dict):
        for key, value in work.items():
            if isinstance(value, (dict, list)):
                v = _peak_memory(value)
            elif str(key) == "memory_bytes" or str(key).endswith(".memory_bytes"):
                v = _num(value)
            else:
                v = None
            if v is not None and (best is None or v > best):
                best = v
    elif isinstance(work, list):
        for value in work:
            v = _peak_memory(value)
            if v is not None and (best is None or v > best):
                best = v
    return best


def figures_of_server(server):
    """The figures a line's `server` block states (sources server / server_entries; null where it states none)."""
    server = server if isinstance(server, dict) else {}
    timing = server.get("timing") if isinstance(server.get("timing"), dict) else {}
    out = _timing_figures(timing, _entry_list(timing, "entries", "by_pattern"), "server")
    out["memory_bytes"] = _peak_memory(server.get("work"))
    out["sources"]["memory_bytes"] = "server" if out["memory_bytes"] is not None else None
    return out


def figures_of_body(body):
    """The figures an answer body states (sources raw / raw_entries): the top-level timing and work, else the
    entries' (patterns of /pattern, results of /traverse)."""
    body = body if isinstance(body, dict) else {}
    entries = _entry_list(body, "patterns", "results")
    out = _timing_figures(body.get("timing"), [e.get("timing") for e in entries if isinstance(e.get("timing"), dict)],
                          "raw")
    mem = _num((body.get("work") or {}).get("memory_bytes")) if isinstance(body.get("work"), dict) else None
    if mem is None:
        vals = [_peak_memory(e.get("work")) for e in entries]
        vals = [v for v in vals if v is not None]
        mem = max(vals) if vals else None
    out["memory_bytes"] = mem
    out["sources"]["memory_bytes"] = "raw" if mem is not None else None
    return out


def host_label(line):
    """The host-cache label the figures key on: compare.timing_cache's (cold, warm, warm after a cut, ...)."""
    import compare as cmp   # compare imports this module: read lazily
    return cmp.timing_cache(line)


def execution_figures(line, run_dir=None, read_raw=True, scenario=None):
    """One execution's figures: {scenario, area, arm, repeat, route, host, outcome, http, cpu_ms, elapsed_ms,
    cpu_ratio, minor_faults, major_faults, memory_bytes, sources}. Each figure from the line's server block, else
    (|read_raw|, a task route) from its raw answer's body; null when neither states it. The route and the area are
    the line's, else its manifest |scenario|'s."""
    scenario = scenario if isinstance(scenario, dict) else {}
    route = line.get("route") or scenario.get("route")
    f = figures_of_server(line.get("server"))
    missing = [k for k in SOURCE_FAMILIES if f["sources"][k] is None]
    if missing and read_raw and route in TASK_ROUTES:
        body = sc.raw_record_body(line, run_dir=run_dir)
        if isinstance(body, dict):
            r = figures_of_body(body)
            if "cpu_ms" in missing and r["sources"]["cpu_ms"]:
                f["cpu_ms"], f["elapsed_ms"], f["sources"]["cpu_ms"] = r["cpu_ms"], r["elapsed_ms"], r["sources"]["cpu_ms"]
            if "page_faults" in missing and r["sources"]["page_faults"]:
                f["minor_faults"], f["major_faults"] = r["minor_faults"], r["major_faults"]
                f["sources"]["page_faults"] = r["sources"]["page_faults"]
            if "memory_bytes" in missing and r["sources"]["memory_bytes"]:
                f["memory_bytes"], f["sources"]["memory_bytes"] = r["memory_bytes"], "raw"
    cpu, el = f["cpu_ms"], f["elapsed_ms"]
    area = line.get("area") or scenario.get("area") or (str(line.get("scenario") or "").split("-", 1)[0] or None)
    return {"scenario": line.get("scenario"), "area": area, "arm": line.get("arm"),
            "repeat": line.get("repeat"), "route": route, "host": host_label(line),
            "outcome": sc.recorded_outcome(line), "http": line.get("http"),
            "cpu_ms": cpu, "elapsed_ms": el,
            "cpu_ratio": (cpu / el) if cpu is not None and el is not None and el > 0 else None,
            "minor_faults": f["minor_faults"], "major_faults": f["major_faults"],
            "memory_bytes": f["memory_bytes"], "sources": f["sources"]}


def run_figures(rows, run_dir=None, scenarios=None, read_raw=True):
    """Every execution's figures (execution_figures), in the lines' order; |scenarios| {id: manifest scenario}
    fills a route or an area a line does not carry."""
    return [execution_figures(line, run_dir, read_raw, (scenarios or {}).get(line.get("scenario")))
            for line in rows if isinstance(line, dict)]


def compact(x):
    """An execution's figures as compare.json keeps them in a scenario's arm view."""
    return dict((k, x.get(k)) for k in ("repeat", "host", "cpu_ms", "elapsed_ms", "cpu_ratio", "minor_faults",
                                        "major_faults", "memory_bytes", "sources"))


# ---------------------------------------------------------------------------------------------------------------
# summaries

def stats(values):
    """{n, p10, p50, p90, max, mean} of the numeric values (sc.summarise and the 10th percentile)."""
    s = sc.summarise([v for v in values if _num(v) is not None])
    s["p10"] = sc.percentile([v for v in values if _num(v) is not None], 10)
    return s


def group_stats(execs):
    """{executions, <figure>: stats} over a group of executions."""
    out = {"executions": len(execs)}
    for k in FIGURES:
        out[k] = stats([x.get(k) for x in execs])
    return out


def _state_sort(state):
    return (STATE_ORDER.index(state) if state in STATE_ORDER else len(STATE_ORDER), str(state))


def per_scenario(execs):
    """{scenario: {arm: {executions: [compact], <host label>: group_stats}}}: per test and access path its
    executions, the host-cold one apart from the warm repeats."""
    out = {}
    for x in execs:
        a = out.setdefault(x["scenario"], {}).setdefault(x["arm"], {"executions": [], "_by": {}})
        a["executions"].append(compact(x))
        a["_by"].setdefault(x["host"], []).append(x)
    for sid in out:
        for arm, a in out[sid].items():
            by = a.pop("_by")
            for state in sorted(by, key=_state_sort):
                a[state] = group_stats(by[state])
    return out


def _sources(execs):
    out = {}
    for x in execs:
        if x.get("route") not in TASK_ROUTES:
            continue
        arm = out.setdefault(x["arm"], dict((k, {}) for k in SOURCE_FAMILIES))
        for k in SOURCE_FAMILIES:
            src = (x.get("sources") or {}).get(k) or "none"
            arm[k][src] = arm[k].get(src, 0) + 1
    return out


def top_tests(execs, figure, n=TOP_N, order=None):
    """The |n| tests with the largest |figure|: per test its execution with the largest value (ties: the first in
    the lines' order), the median over its executions that state the figure and how many do, and the host-cold
    execution's value apart from the warm repeats' median (`cold`, `warm_median`; null without one); ranked by the
    value, then by the manifest |order|, then by id."""
    pos = dict((sid, i) for i, sid in enumerate(order or []))
    best, vals, cold, warm = {}, {}, {}, {}
    for x in execs:
        v = x.get(figure)
        if _num(v) is None:
            continue
        sid = x["scenario"]
        vals.setdefault(sid, []).append(v)
        if x.get("host") == "cold":
            cold.setdefault(sid, []).append(v)
        elif x.get("host") == "warm":
            warm.setdefault(sid, []).append(v)
        if sid not in best or v > best[sid][figure]:
            best[sid] = x
    ranked = sorted(best, key=lambda s: (-best[s][figure], pos.get(s, len(pos)), str(s)))
    out = []
    for sid in ranked[:n]:
        x = best[sid]
        out.append({"scenario": sid, "area": x.get("area"), "value": x[figure], "arm": x.get("arm"),
                    "repeat": x.get("repeat"), "host": x.get("host"), "route": x.get("route"),
                    "elapsed_ms": x.get("elapsed_ms"), "cpu_ms": x.get("cpu_ms"), "cpu_ratio": x.get("cpu_ratio"),
                    "major_faults": x.get("major_faults"), "memory_bytes": x.get("memory_bytes"),
                    "median": sc.percentile(vals[sid], 50), "n": len(vals[sid]),
                    "cold": max(cold[sid]) if sid in cold else None,
                    "warm_median": sc.percentile(warm[sid], 50) if sid in warm else None})
    return out


def summarise_run(execs, order=None):
    """The run's block: {format, routes, executions {arm: n of task-route executions}, stated {arm: {figure: n}},
    sources {arm: {family: {source: n}}}, by_area {area: {arm: {host label: group_stats}}} (area "all": every
    area), top {figure: [...]}}."""
    task = [x for x in execs if x.get("route") in TASK_ROUTES]
    by_area = {}
    for x in task:
        for area in (x.get("area") or "?", "all"):
            by_area.setdefault(area, {}).setdefault(x["arm"], {}).setdefault(x["host"], []).append(x)
    by_area = dict((area, dict((arm, dict((st, group_stats(v)) for st, v in sorted(states.items(),
                                                                                   key=lambda kv: _state_sort(kv[0]))))
                               for arm, states in arms.items())) for area, arms in by_area.items())
    executions, stated = {}, {}
    for x in task:
        executions[x["arm"]] = executions.get(x["arm"], 0) + 1
        st = stated.setdefault(x["arm"], dict((k, 0) for k in ("cpu_ms", "major_faults", "memory_bytes")))
        for k in st:
            if _num(x.get(k)) is not None:
                st[k] += 1
    return {"format": FORMAT, "routes": list(TASK_ROUTES), "executions": executions, "stated": stated,
            "sources": _sources(task), "by_area": by_area,
            "top": dict((f, top_tests(task, f, TOP_N, order)) for f in TOP_FIGURES)}


def any_stated(summary, figures=("cpu_ms", "major_faults", "memory_bytes")):
    """Whether some execution of the run states one of |figures|."""
    return any(any((st or {}).get(k) for k in figures) for st in ((summary or {}).get("stated") or {}).values())


def of_run(run_dir, rows, manifest=None):
    """The run's block computed from its lines (summarise_run over run_figures), whatever its compare.json holds:
    the report and the document read the same figures however the comparison was computed."""
    order, by_id = [], {}
    scenarios = manifest.get("scenarios") if isinstance(manifest, dict) else None
    for s in scenarios if isinstance(scenarios, list) else []:
        if isinstance(s, dict) and s.get("id"):
            order.append(s["id"])
            by_id[s["id"]] = s
    return summarise_run(run_figures(rows, run_dir, by_id), order)


def engine_states_figures(rows):
    """Whether some line's server block states a figure only the per-task engine states (timing.cpu_ms or
    timing.page_faults, top level or per entry): read without the raw records (report.py's cheap test)."""
    for line in rows:
        f = figures_of_server((line or {}).get("server") if isinstance(line, dict) else None)
        if f["sources"]["cpu_ms"] or f["sources"]["page_faults"]:
            return True
    return False


# ---------------------------------------------------------------------------------------------------------------
# the engine container's memory over the run (host_memory.json)

def read_host_memory(run_dir):
    """(the parsed host_memory.json or None, why not: '' when the run directory holds none, else the problem)."""
    path = os.path.join(run_dir, HOST_MEMORY_FILE)
    if not os.path.exists(path):
        return None, ""
    try:
        doc = sc.read_json(path)
    except (OSError, ValueError) as exc:
        return None, "%s cannot be read: %s" % (HOST_MEMORY_FILE, exc)
    if not isinstance(doc, dict) or not isinstance(doc.get("before"), dict) or not isinstance(doc.get("after"), dict):
        return None, "%s is not {\"before\": SNAP, \"after\": SNAP}" % HOST_MEMORY_FILE
    return doc, ""


def _level(snap, key):
    v = snap.get(key)
    if key == "max" and (v is None or str(v).strip() == "max"):
        return None
    return _num(v)


def host_memory_figures(doc):
    """The figures of a host_memory.json: {container, page_bytes, before {label, utc}, after {...}, levels {key:
    {before, after}}, deltas {key: after - before} (null where a snapshot lacks the counter or it went down: a
    container restarted between the snapshots, `restarted` then true), pages_gib {key: GiB of the delta's pages},
    limit_bytes (after's max, else before's), notes [...]}."""
    b, a = doc.get("before") or {}, doc.get("after") or {}
    page = _num(a.get("page_bytes")) or _num(b.get("page_bytes")) or PAGE_BYTES
    levels = dict((k, {"before": _level(b, k), "after": _level(a, k)}) for k in HOST_LEVELS)
    deltas, notes, restarted = {}, [], False
    for k in HOST_COUNTERS:
        vb, va = _num(b.get(k)), _num(a.get(k))
        if vb is None or va is None:
            deltas[k] = None
        elif va < vb:
            deltas[k] = None
            restarted = True
        else:
            deltas[k] = va - vb
    if restarted:
        notes.append("a counter went down between the snapshots: the container restarted, so the deltas over the "
                     "run are not stated")
    if b.get("container") and a.get("container") and b.get("container") != a.get("container"):
        notes.append("the snapshots name different containers (%s, %s)" % (b.get("container"), a.get("container")))
    pages_gib = dict((k, (deltas[k] * page / float(1 << 30)) if deltas.get(k) is not None else None)
                     for k in HOST_PAGE_COUNTERS)
    limit = levels["max"]["after"] if levels["max"]["after"] is not None else levels["max"]["before"]
    return {"container": a.get("container") or b.get("container"), "page_bytes": page,
            "before": {"label": b.get("label"), "utc": b.get("utc")},
            "after": {"label": a.get("label"), "utc": a.get("utc")},
            "levels": levels, "deltas": deltas, "pages_gib": pages_gib, "limit_bytes": limit,
            "restarted": restarted, "notes": notes}


def host_memory_of_run(run_dir):
    """{figures or None, problem}: the run directory's host_memory.json read and turned into figures."""
    doc, why = read_host_memory(run_dir)
    return {"figures": host_memory_figures(doc) if doc is not None else None, "problem": why or None}


# ---------------------------------------------------------------------------------------------------------------
# area G: the service's graphlet tools

def graphlet_service_figures(run_dir):
    """What the service's graphlet-tool answers (ways GRAPHLET_SERVICE_WAYS, arms/mcp/graphlet_ops.jsonl) state of
    their own cost: {ops, stating, cpu_ms stats} -- `stating` counts the answers that carry a top-level `usage` or
    `cpu_ms`; None when the run holds no such records."""
    path = os.path.join(run_dir, "arms", "mcp", "graphlet_ops.jsonl")
    rows, _ = sc.read_jsonl(path)
    ops = [r for r in rows if isinstance(r, dict) and r.get("way") in GRAPHLET_SERVICE_WAYS]
    if not ops:
        return None
    stating, cpu = 0, []
    for r in ops:
        ans = r.get("answer") if isinstance(r.get("answer"), dict) else {}
        if any(k in ans for k in GRAPHLET_FIGURE_KEYS):
            stating += 1
            v = _num(ans.get("cpu_ms"))
            if v is None and isinstance(ans.get("usage"), dict):
                v = _num(ans["usage"].get("cpu_ms"))
            if v is not None:
                cpu.append(v)
    return {"ops": len(ops), "stating": stating, "cpu_ms": stats(cpu)}
