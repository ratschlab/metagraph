#!/usr/bin/env python3
"""Builds the isolated benchmark panel (lane B): a fixed set of successful, decisive queries chosen from a recorded
three-arm suite run, with the result sizes each arm recorded, fixed weights and the protocol's costing.

    make_panel.py --baseline RUN_DIR [--reevaluated] [--manifest MANIFEST] --out PANEL.json
                   [--size 48] [--quota resolve=16,traverse=16,pattern=16]
                   [--min-tags label_heavy=3,long_path=4,long_query=2,large_answer=2,heavy=2] [--max-empty 2]
                   [--include IDS] [--exclude IDS] [--include-from PANEL] [--name NAME]

RUN_DIR is a merged suite run (executions.jsonl, compare.json, run.json, the arms' raw answers; with --reevaluated
its re-evaluated set, executions.reevaluated.jsonl and compare.reevaluated.json, judged by the checks of the
re-evaluation); MANIFEST the frozen manifest it ran (default RUN_DIR/manifest.json, which must be the same: its
content_sha256 is compared with the run's). Re-runnable: v5 replaces v4 with `--baseline <a v5 run>`, and
--include-from the v4 panel keeps the v4 queries that are still eligible.

A query is ELIGIBLE (each failed rule is counted in the panel's `selection.excluded`):
  route      POST /resolve, /traverse or /pattern (a search, not a capabilities or protocol case);
  protocol   no body_raw, headers or run_unique (the request is the scenario's, the same each execution);
  long       not the long class (every arm runs it: the front cuts the long class);
  mcp        mapped to an MCP tool and not marked not implemented via MCP;
  missing    every arm recorded it;
  outcome    every recorded execution of every arm passed (decisive: no inconclusive, by-design or failure);
  http       every one answered HTTP 200 (a successful search, not a refusal);
  time       none was time-limited (the work does not depend on the host's speed);
  agreement  the direct arms identical, the MCP arm agrees, each arm's repeats with one digest;
  complete   the direct answer is a complete search (bench_lib.answer_incomplete: no per-pattern refusal, withheld
             result, stop or cut, no truncated or stopped traversal arm, no truncated labels): capped queries are
             budget stress, which the protocol runs apart;
  duplicate  not the same request as an eligible query before it (the D area repeats other areas' requests).
The answers' result sizes (bench_lib.result_sizes) are recorded per arm: a measured execution counts as the same
workload only when its sizes equal these.

SELECTION per mode (quota, default 16 each): the mode's eligible queries ranked by their recorded cost (the least
wall of the direct arms' recorded executions: ranking only, the baseline ran under the Logan load) and placed in
`quota` bins of equal width in log10(cost), so the panel spans the cost range. First the --include queries, then the
tag minimums (label_heavy, long_path, long_query, large_answer, heavy = the top tenth of the mode's cost), each
from its tagged queries at evenly spaced cost quantiles; then one query per bin holding the fewest picks, the bins
farthest from the chosen ones first, nearest the bin's centre. Every pick prefers an area of the mode not yet chosen
(the pattern route serves P, L, X, M, Q and N), then a family not yet chosen, then a non-empty result (at most
--max-empty empty results per mode). Deterministic: no random choice; ties go to the lower id.

WEIGHTS are fixed at build time: equal per query (1 / size); the analysis aggregates with them and never re-weights
silently. The COSTING states the protocol's requests (size x (1 + 6 x 3) = 912 for 48 queries) and an execution-time
estimate from the baseline's recorded walls and step gaps (a loaded host's: an upper estimate), against the planning
allowance of 1.5-2 host-hours per release. Stdlib only.
"""

import argparse
import datetime
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bench_lib as bl  # noqa: E402

sc = bl.sc

DEFAULT_SIZE = 48
DEFAULT_QUOTA = {"resolve": 16, "traverse": 16, "pattern": 16}
DEFAULT_MIN_TAGS = {"label_heavy": 3, "long_path": 4, "long_query": 2, "large_answer": 2, "heavy": 2}
DEFAULT_MAX_EMPTY = 2
HEAVY_FRACTION = 0.10
# Costs below this (ms) share the cheapest bin (a log scale needs a floor).
COST_FLOOR_MS = 0.05
TAG_ORDER = ("label_heavy", "long_path", "long_query", "large_answer", "heavy")
ELIGIBILITY_RULES = ("route", "protocol", "long", "mcp", "missing", "outcome", "http", "time", "agreement",
                     "complete", "duplicate")
# The planning allowance of GPT's protocol (host reservation hours per release, not CPU-core hours).
ALLOWANCE_HOST_HOURS = [1.5, 2.0]
# What a step adds besides its request when the baseline does not tell (seconds; v4's median gaps were 1.4 s before
# a direct step and 1.7-1.9 s before an MCP step, with a 1 s pause).
FALLBACK_STEP_GAP_S = {"direct-local": 1.4, "direct-remote": 1.4, "mcp": 1.9}


def parse_map(text, kind=int):
    out = {}
    for part in (text or "").split(","):
        part = part.strip()
        if not part:
            continue
        key, _, value = part.partition("=")
        out[key.strip()] = kind(value)
    return out


def parse_ids(text):
    return [x for x in (text or "").replace("\n", ",").replace(" ", ",").split(",") if x]


# ---------------------------------------------------------------------------------------------------------------
# the baseline run

def load_baseline(run_dir, reevaluated=False):
    """The baseline's recorded executions per scenario and arm, its comparison per scenario, its identity."""
    ex_path = os.path.join(run_dir, sc.reevaluated_name("executions.jsonl", reevaluated))
    cmp_path = os.path.join(run_dir, sc.reevaluated_name("compare.json", reevaluated))
    if not os.path.exists(ex_path) or not os.path.exists(cmp_path):
        raise SystemExit("make_panel: %s needs %s and %s (a merged and compared run%s)" % (
            run_dir, os.path.basename(ex_path), os.path.basename(cmp_path),
            "; --reevaluated reads a re-evaluation's" if reevaluated else ""))
    rows, bad = sc.read_jsonl(ex_path)
    by = {}
    for r in rows:
        by.setdefault(r.get("scenario"), {}).setdefault(r.get("arm"), []).append(r)
    for arms in by.values():
        for lines in arms.values():
            lines.sort(key=lambda r: (r.get("repeat") or 0))
    compare = sc.read_json(cmp_path) or {}
    run = sc.read_json(os.path.join(run_dir, "run.json")) or {}
    orch = sc.read_json(os.path.join(run_dir, "orchestrator.json")) or {}
    reev = sc.read_json(os.path.join(run_dir, "reevaluated.json")) if reevaluated else None
    offset = None
    for sess in orch.get("sessions") or []:
        if isinstance(sess.get("mex_clock_offset_s"), (int, float)):
            offset = sess["mex_clock_offset_s"]
            break
    return {"run_dir": os.path.abspath(run_dir), "rows": rows, "bad_lines": bad, "by": by,
            "compare": compare.get("scenarios") or {}, "compare_manifest_sha256": compare.get("manifest_sha256"),
            "run_id": run.get("run_id") or compare.get("run_id") or os.path.basename(os.path.normpath(run_dir)),
            "release": run.get("release") or compare.get("release"), "index_fp": run.get("index_fp"),
            "server_instance": run.get("server_instance") or compare.get("server_instance"),
            "executions_file": os.path.basename(ex_path), "executions_sha256": sc.file_sha256(ex_path),
            "checks": (reev or {}).get("checks"), "reevaluated": bool(reevaluated), "mex_clock_offset_s": offset,
            "load": load_context(os.path.join(run_dir, "samples.jsonl"))}


def load_context(path):
    """The baseline's load while it ran (the sampler's samples): medians of the Logan container's CPU and the host's
    load1, for the panel's honesty note (the recorded walls are a loaded host's)."""
    rows, _ = sc.read_jsonl(path)
    if not rows:
        return None
    return {"samples": len(rows),
            "sra_all_cpu_pct_median": bl.median([r.get("sra_all_cpu_pct") for r in rows]),
            "host_load1_median": bl.median([r.get("host_load1") for r in rows])}


def raw_body(run_dir, line):
    """The answer body of a recorded execution (its raw file, relative to the run directory), None when unreadable."""
    rel = line.get("raw")
    if not rel:
        return None
    for path in (os.path.join(run_dir, rel), os.path.join(run_dir, "arms", line.get("arm") or "", rel)):
        if os.path.exists(path):
            try:
                doc = sc.read_gz_json(path)
            except (OSError, ValueError, EOFError):
                return None
            return ((doc or {}).get("answer") or {}).get("body")
    return None


# ---------------------------------------------------------------------------------------------------------------
# eligibility and the candidates

def ineligible(scenario, arms, cmp_entry):
    """The eligibility rules a scenario fails (ELIGIBILITY_RULES, without `duplicate`), [] when it is eligible."""
    why = []
    if bl.mode_of(scenario.get("route")) is None:
        why.append("route")
    if scenario.get("body_raw") is not None or scenario.get("headers") or scenario.get("run_unique"):
        why.append("protocol")
    if scenario.get("timing_class") == "long":
        why.append("long")
    mcp = scenario.get("mcp") or {}
    if (not mcp.get("tool") and not mcp.get("flow")) or mcp.get("not_implemented"):
        why.append("mcp")
    lines = [r for arm in bl.ARMS for r in arms.get(arm) or []]
    if any(not arms.get(arm) for arm in bl.ARMS):
        why.append("missing")
    if any(r.get("outcome") != "pass" for r in lines):
        why.append("outcome")
    if any(r.get("http") != 200 for r in lines):
        why.append("http")
    if any(r.get("time_limited") for r in lines):
        why.append("time")
    da = (cmp_entry or {}).get("direct_agreement") or {}
    ma = (cmp_entry or {}).get("mcp_agreement") or {}
    digests = [((cmp_entry or {}).get("arms") or {}).get(arm, {}).get("digests_equal") for arm in bl.ARMS]
    # digests_equal is None for an arm that ran the scenario once (nothing to compare), False when repeats differ
    if da.get("class") != "identical" or ma.get("class") != "agree" or any(d is False for d in digests):
        why.append("agreement")
    return why


def _server_ms(line):
    server = line.get("server") or {}
    v = server.get("elapsed_ms") if isinstance(server, dict) else None
    return v if isinstance(v, (int, float)) and not isinstance(v, bool) else None


def candidate(run_dir, scenario, arms):
    """The panel record of an eligible scenario: mode, family, tags, the recorded cost and walls, and per arm the
    result sizes of its recorded answer (repeat 0; the arm's repeats share one digest)."""
    sid, route = scenario["id"], scenario["route"]
    expected, work, digests = {}, {}, {}
    for arm in bl.ARMS:
        first = arms[arm][0]
        body = raw_body(run_dir, first)
        expected[arm] = {"http": first.get("http"), "sizes": bl.result_sizes(route, body) if body is not None else None,
                         "bytes_body": (first.get("client") or {}).get("bytes_body")}
        digests[arm] = first.get("answer_digest")
        if arm == "direct-local":
            work = bl.work_counters(route, body) if body is not None else {}
            tags = bl.query_tags(route, body, expected[arm]["sizes"] or {}, expected[arm]["bytes_body"]) \
                if body is not None else []
    # the least wall of the direct arms' recorded executions: load only adds time, so the least is the best estimate
    # of the query's own cost on a loaded host (ranking only). Not the server-stated elapsed time: some answers state
    # less than their request took (a predicate's row reads), and not the MCP arm's walls (its overhead dominates).
    walls = [r.get("wall_ms") for arm in sc.DIRECT_ARMS for r in arms[arm]
             if isinstance(r.get("wall_ms"), (int, float)) and not isinstance(r.get("wall_ms"), bool)]
    cost = min(walls) if walls else None
    recorded = {}
    for arm in bl.ARMS:
        rec = {"wall_ms": [r.get("wall_ms") for r in arms[arm]], "server_elapsed_ms": [_server_ms(r) for r in arms[arm]]}
        if arm == "mcp":
            phases = {}
            for r in arms[arm]:
                for k, v in (r.get("mcp_phases") or {}).items():
                    phases.setdefault(k, []).append(v)
            rec["mcp_phases_median_ms"] = dict((k, bl.median(v)) for k, v in sorted(phases.items()))
        recorded[arm] = rec
    return {"id": sid, "mode": bl.mode_of(route), "area": scenario.get("area"), "family": bl.family(sid),
            "route": route, "timing_class": scenario.get("timing_class"), "title": scenario.get("title"),
            "tags": tags, "cost_ms": cost, "expected": expected, "work": work, "digests": digests,
            "recorded": recorded, "scenario": scenario}


def candidates(manifest, baseline):
    """(eligible candidates in manifest order, excluded {rule: count}, excluded ids {rule: [ids]})."""
    excluded, excluded_ids, out, seen = {}, {}, [], {}
    # D-area scenarios repeat other areas' requests: judge them last, so the original is the one kept
    order = sorted(manifest.get("scenarios") or [], key=lambda s: (str(s.get("area")) == "D",))
    for scenario in order:
        sid = scenario.get("id")
        arms = baseline["by"].get(sid) or {}
        why = ineligible(scenario, arms, baseline["compare"].get(sid))
        if not why:
            # the direct answer (repeat 0) a complete, successful search: no refusal, withheld result, stop or cap
            body = raw_body(baseline["run_dir"], arms["direct-local"][0])
            if body is None or bl.answer_incomplete(scenario.get("route"), body):
                why = ["complete"]
        if not why:
            key = bl.request_key(scenario)
            if key in seen:
                why = ["duplicate"]
            else:
                seen[key] = sid
        for rule in why:
            excluded[rule] = excluded.get(rule, 0) + 1
            excluded_ids.setdefault(rule, []).append(sid)
        if not why:
            out.append(candidate(baseline["run_dir"], scenario, arms))
    position = dict((s.get("id"), i) for i, s in enumerate(manifest.get("scenarios") or []))
    out.sort(key=lambda c: position[c["id"]])
    return out, excluded, excluded_ids


# ---------------------------------------------------------------------------------------------------------------
# the selection

def cost_bins(pool, n):
    """The bin (0..n-1) of each query: n equal-width bins of log10(cost) between the pool's cheapest and dearest
    query, so the bins span the cost RANGE (cheap requests measure the overhead, dear ones the engine) instead of
    following the suite's many cheap scenarios. A query without a cost falls in the last bin."""
    costs = [math.log10(max(c["cost_ms"], COST_FLOOR_MS)) for c in pool if c["cost_ms"] is not None]
    lo, hi = (min(costs), max(costs)) if costs else (0.0, 0.0)
    width = (hi - lo) / float(n) if hi > lo else 1.0
    out = {}
    for c in pool:
        if c["cost_ms"] is None:
            out[c["id"]] = n - 1
            continue
        x = math.log10(max(c["cost_ms"], COST_FLOOR_MS))
        out[c["id"]] = min(n - 1, int((x - lo) / width)) if hi > lo else 0
    return out, lo, width


def mark_heavy(pool):
    """Tags the top HEAVY_FRACTION of a mode's pool by cost `heavy` (at least one)."""
    ranked = sorted(pool, key=lambda c: (-(c["cost_ms"] or 0), c["id"]))
    n = max(1, int(round(len(ranked) * HEAVY_FRACTION))) if ranked else 0
    for c in ranked[:n]:
        if "heavy" not in c["tags"]:
            c["tags"].append("heavy")


def select_mode(pool, quota, min_tags, max_empty, include=()):
    """The quota's queries of one mode's pool, and notes on what could not be met. Each pick prefers, in order: an
    area of the mode not yet chosen (the pattern route serves P, L, X, M, Q and N), a family not yet chosen, a
    non-empty result, then the cost the step wants (a tag's quantile, a bin's centre), then the lower id. At most
    max_empty empty results per mode, the tag picks included."""
    pool = sorted(pool, key=lambda c: (c["cost_ms"] if c["cost_ms"] is not None else float("inf"), c["id"]))
    notes = []
    if len(pool) <= quota:
        if len(pool) < quota:
            notes.append("only %d eligible queries for a quota of %d" % (len(pool), quota))
        for c in pool:
            c.setdefault("selected_by", "all")
        return list(pool), notes
    picks = []
    areas, families = {}, {}

    def take(c, why):
        picks.append(c)
        areas[c["area"]] = areas.get(c["area"], 0) + 1
        families[c["family"]] = families.get(c["family"], 0) + 1
        c["selected_by"] = why

    def empties():
        return sum(1 for c in picks if "empty_result" in c["tags"])

    def allowed(c):
        return c not in picks and ("empty_result" not in c["tags"] or empties() < max_empty)

    def preference(c, distance):
        return (areas.get(c["area"], 0), families.get(c["family"], 0), "empty_result" in c["tags"], distance, c["id"])

    for c in pool:
        if c["id"] in set(include):
            take(c, "include")
    if len(picks) > quota:
        raise SystemExit("make_panel: %d included queries of mode %s exceed its quota %d" % (
            len(picks), pool[0]["mode"], quota))
    rank = dict((c["id"], i) for i, c in enumerate(pool))
    for tag in TAG_ORDER:
        want = min_tags.get(tag, 0)
        need = want - sum(1 for c in picks if tag in c["tags"])
        tagged = [c for c in pool if tag in c["tags"]]
        if need <= 0:
            continue
        if not tagged:
            continue
        free = [c for c in tagged if c not in picks]
        if len(free) < need:
            notes.append("tag %s: %d of %d wanted available" % (tag, want - need + len(free), want))
        n = min(need, len(free))
        for i in range(n):
            if len(picks) >= quota:
                notes.append("tag %s: the quota is full" % tag)
                break
            rest = [c for c in free if allowed(c)]
            if not rest:
                notes.append("tag %s: only empty results left past the cap" % tag)
                break
            # evenly spaced cost quantiles over the tagged queries
            target = rank[free[max(0, min(len(free) - 1, int(round((i + 0.5) * len(free) / float(n) - 0.5))))]["id"]]
            take(min(rest, key=lambda c: preference(c, abs(rank[c["id"]] - target))), "tag:" + tag)
    bins, lo, width = cost_bins(pool, quota)
    members = {}
    for c in pool:
        members.setdefault(bins[c["id"]], []).append(c)
    while len(picks) < quota:
        held = {}
        for c in picks:
            held[bins[c["id"]]] = held.get(bins[c["id"]], 0) + 1
        open_bins = [b for b, cs in members.items() if any(allowed(c) for c in cs)]
        if not open_bins:
            notes.append("no query left within the empty-result cap")
            break
        least = min(held.get(b, 0) for b in open_bins)
        candidates_bins = [b for b in open_bins if held.get(b, 0) == least]

        def distance(b):
            return min([abs(b - h) for h in held] or [quota])
        b = min(candidates_bins, key=lambda x: (-distance(x), x))
        centre = lo + (b + 0.5) * width
        take(min([c for c in members[b] if allowed(c)],
                 key=lambda c: preference(c, abs(math.log10(max(c["cost_ms"] or COST_FLOOR_MS, COST_FLOOR_MS))
                                                 - centre))), "bin:%d" % b)
    for c in picks:
        c["cost_bin"] = bins[c["id"]]
    return sorted(picks, key=lambda c: (c["cost_ms"] if c["cost_ms"] is not None else float("inf"), c["id"])), notes


# ---------------------------------------------------------------------------------------------------------------
# costing

def step_gaps(baseline):
    """The baseline's median gap (s) before an execution of each arm, from the previous execution's end (the
    orchestration, the runner's start-up, the pause): what a step adds besides its request."""
    events = []
    offset = baseline.get("mex_clock_offset_s") or 0.0
    for r in baseline["rows"]:
        t0, t1 = sc.parse_utc(r.get("started_utc")), sc.parse_utc(r.get("ended_utc"))
        if t0 is None or t1 is None:
            continue
        if r.get("arm") == "direct-local":
            t0, t1 = t0 - offset, t1 - offset
        events.append((t0, t1, r.get("arm")))
    events.sort()
    gaps = {}
    for a, b in zip(events, events[1:]):
        g = b[0] - a[1]
        if 0 <= g < 600:
            gaps.setdefault(b[2], []).append(g)
    return dict((arm, round(bl.median(v), 3)) for arm, v in gaps.items() if v)


def costing(chosen, baseline, cycles=bl.DEFAULT_CYCLES, warmup_arm=bl.DEFAULT_WARMUP_ARM):
    """The protocol's requests and an execution-time estimate from the baseline's recorded walls (each arm's median
    per query) and step gaps: an upper estimate, the baseline ran under the Logan load."""
    gaps = step_gaps(baseline)
    gap = dict((arm, gaps.get(arm, FALLBACK_STEP_GAP_S[arm])) for arm in bl.ARMS)
    exec_s = overhead_s = 0.0
    for c in chosen:
        med = dict((arm, (bl.median(c["recorded"][arm]["wall_ms"]) or 0.0) / 1000.0) for arm in bl.ARMS)
        exec_s += med[warmup_arm] + cycles * sum(med.values())
        overhead_s += gap[warmup_arm] + cycles * sum(gap.values())
    n = len(chosen)
    total_h = (exec_s + overhead_s) / 3600.0
    return {"requests": bl.requests_planned(n, cycles), "warmup_requests": n, "measured_requests": n * cycles * 3,
            "formula": "%d x (1 + %d x 3)" % (n, cycles),
            "other_requests": "per direct arm 3 capabilities GETs at its first step and one closing GET; one "
                              "direct-local GET /capabilities per block (the identity probe); per MCP step the "
                              "runner's service capabilities read and queue check (service-side)",
            "step_gap_s": gap, "step_gap_source": "the baseline's median gaps" if gaps else "fallback values",
            "execution_s_estimate": round(exec_s, 1), "overhead_s_estimate": round(overhead_s, 1),
            "host_hours_estimate": round(total_h, 2),
            "allowance_host_hours": ALLOWANCE_HOST_HOURS,
            "comparison_allowance_host_hours": [3.0, 4.0],
            "basis": "each arm's median recorded wall per query (the baseline's, under the Logan load: an upper "
                     "estimate) and the baseline's step gaps; the allowance adds 0.5-1 h of setup, isolation and "
                     "transitions. Re-estimate after one isolated pilot."}


# ---------------------------------------------------------------------------------------------------------------
# the panel

def build(manifest, baseline, size=DEFAULT_SIZE, quota=None, min_tags=None, max_empty=DEFAULT_MAX_EMPTY,
          include=(), exclude=(), name=None, manifest_path=None, now=None):
    quota = dict(quota or DEFAULT_QUOTA)
    if sum(quota.values()) != size:
        raise SystemExit("make_panel: the quota %s does not add up to --size %d" % (quota, size))
    min_tags = dict(DEFAULT_MIN_TAGS if min_tags is None else min_tags)
    cands, excluded, excluded_ids = candidates(manifest, baseline)
    known = set(s.get("id") for s in manifest.get("scenarios") or [])
    eligible = set(c["id"] for c in cands)
    include_missing = [i for i in include if i not in eligible]
    cands = [c for c in cands if c["id"] not in set(exclude)]
    chosen, notes = [], []
    pools = {}
    for mode in bl.MODES:
        pool = [c for c in cands if c["mode"] == mode]
        mark_heavy(pool)
        pools[mode] = len(pool)
        picks, mode_notes = select_mode(pool, quota.get(mode, 0), min_tags, max_empty, include)
        notes += ["%s: %s" % (mode, n) for n in mode_notes]
        chosen += picks
    if len(chosen) != size:
        notes.append("the panel holds %d queries, not %d" % (len(chosen), size))
    weight = 1.0 / len(chosen) if chosen else 0.0
    queries = []
    for c in chosen:
        queries.append({"id": c["id"], "mode": c["mode"], "area": c["area"], "family": c["family"],
                        "route": c["route"], "timing_class": c["timing_class"], "title": c["title"],
                        "tags": sorted(c["tags"]), "selected_by": c.get("selected_by"), "cost_bin": c.get("cost_bin"),
                        "weight": weight, "cost_ms": c["cost_ms"], "expected": c["expected"], "work": c["work"],
                        "baseline_digests": c["digests"], "recorded": c["recorded"],
                        "scenario_sha256": bl.scenario_sha256(c["scenario"]), "scenario": c["scenario"]})
    ident = sc.manifest_identity(manifest)
    by_mode = dict((m, sum(1 for q in queries if q["mode"] == m)) for m in bl.MODES)
    by_tag = {}
    for q in queries:
        for t in q["tags"]:
            by_tag[t] = by_tag.get(t, 0) + 1
    panel = {
        "panel_version": bl.PANEL_VERSION,
        "name": name or "panel",
        "created_utc": sc.iso_utc(now if now is not None else datetime.datetime.now().timestamp()),
        "builder": {"script": "bench/make_panel.py", "sha256": sc.file_sha256(os.path.abspath(__file__)),
                    "bench_lib_sha256": sc.file_sha256(os.path.join(HERE, "bench_lib.py"))},
        "source": {
            "manifest": {"path": os.path.abspath(manifest_path) if manifest_path else None,
                         "manifest_version": manifest.get("manifest_version"),
                         "content_sha256": ident["content_sha256"], "relations_sha256": ident["relations_sha256"],
                         "stated_content_sha256": manifest.get("content_sha256"), "index": manifest.get("index")},
            "baseline": {"run_id": baseline["run_id"], "run_dir": baseline["run_dir"],
                         "reevaluated": baseline["reevaluated"], "executions": baseline["executions_file"],
                         "executions_sha256": baseline["executions_sha256"], "checks": baseline["checks"],
                         "release": baseline["release"], "index_fp": baseline["index_fp"],
                         "server_instance": baseline["server_instance"], "load": baseline["load"],
                         "load_note": "the recorded walls and costs are this run's, under its load: they rank the "
                                      "queries and size the costing; they are no performance baseline"}},
        "selection": {"rules": list(ELIGIBILITY_RULES), "size": size, "quota": quota, "min_tags": min_tags,
                      "max_empty_per_mode": max_empty, "heavy_fraction": HEAVY_FRACTION,
                      "tag_thresholds": bl.TAG_THRESHOLDS, "include": list(include), "exclude": list(exclude),
                      "include_not_eligible": include_missing,
                      "include_unknown": [i for i in include if i not in known],
                      "scenarios": len(manifest.get("scenarios") or []), "eligible": len(eligible),
                      "eligible_by_mode": pools, "excluded": excluded,
                      "excluded_ids": dict((k, sorted(v)) for k, v in excluded_ids.items()),
                      "by_mode": by_mode, "by_tag": by_tag, "families": len(set(q["family"] for q in queries)),
                      "notes": notes},
        "protocol": {"warmup_per_query": 1, "warmup_arm": bl.DEFAULT_WARMUP_ARM, "cycles": bl.DEFAULT_CYCLES,
                     "arms": list(bl.ARMS), "arm_orders": bl.arm_orders(), "one_request_at_a_time": True,
                     "requests": bl.requests_planned(len(queries))},
        "weights": {"scheme": "equal", "per_query": weight,
                    "note": "fixed at build time; the analysis aggregates every arm and comparison with these "
                            "weights and states when a query had no valid observation instead of re-weighting"},
        "costing": costing(chosen, baseline),
        "queries": queries,
    }
    panel["panel_sha256"] = bl.panel_identity(queries)
    return panel


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--baseline", required=True, help="a merged, compared suite run directory")
    ap.add_argument("--reevaluated", action="store_true", help="read its re-evaluated set")
    ap.add_argument("--manifest", default=None, help="the manifest it ran (default RUN_DIR/manifest.json)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--size", type=int, default=DEFAULT_SIZE)
    ap.add_argument("--quota", default=None, help="mode=count,... (default resolve=16,traverse=16,pattern=16)")
    ap.add_argument("--min-tags", default=None,
                    help="tag=count,... (default label_heavy=3,long_path=4,long_query=2,large_answer=2,heavy=2)")
    ap.add_argument("--max-empty", type=int, default=DEFAULT_MAX_EMPTY)
    ap.add_argument("--include", default="", help="ids to include (comma-separated) when eligible")
    ap.add_argument("--include-from", default=None, help="an earlier panel whose queries to include when eligible")
    ap.add_argument("--exclude", default="", help="ids never to choose")
    ap.add_argument("--name", default=None)
    args = ap.parse_args(argv)
    manifest_path = args.manifest or os.path.join(args.baseline, "manifest.json")
    manifest = sc.read_json(manifest_path)
    if not isinstance(manifest, dict) or "scenarios" not in manifest:
        raise SystemExit("make_panel: %s is not a suite manifest" % manifest_path)
    baseline = load_baseline(args.baseline, args.reevaluated)
    computed = sc.manifest_identity(manifest)["content_sha256"]
    if baseline["compare_manifest_sha256"] and baseline["compare_manifest_sha256"] != computed:
        raise SystemExit("make_panel: %s ran another manifest (scenarios sha %s) than %s (%s)" % (
            args.baseline, baseline["compare_manifest_sha256"], manifest_path, computed))
    include = parse_ids(args.include)
    if args.include_from:
        include += [q["id"] for q in bl.read_panel(args.include_from)["queries"] if q["id"] not in include]
    quota = parse_map(args.quota) if args.quota else None
    min_tags = parse_map(args.min_tags) if args.min_tags is not None else None
    panel = build(manifest, baseline, args.size, quota, min_tags, args.max_empty, include, parse_ids(args.exclude),
                  args.name, manifest_path)
    sc.write_json_atomic(args.out, panel)
    sel = panel["selection"]
    print("panel %s: %d queries (%s) from %d eligible of %d scenarios; tags %s; families %d; %d requests, "
          "~%.2f host-hours estimated (allowance %s)" % (
              args.out, len(panel["queries"]), ", ".join("%s %d" % kv for kv in sel["by_mode"].items()),
              sel["eligible"], sel["scenarios"], sel["by_tag"], sel["families"], panel["costing"]["requests"],
              panel["costing"]["host_hours_estimate"], "-".join(str(x) for x in ALLOWANCE_HOST_HOURS)))
    for note in sel["notes"]:
        print("note: %s" % note)
    if sel["include_not_eligible"]:
        print("note: included but not eligible: %s" % ", ".join(sel["include_not_eligible"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
