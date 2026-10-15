#!/usr/bin/env python3
"""Shared helpers of the isolated benchmark panel (lane B): the panel's query modes, result sizes and tags, the run
manifest a panel derives (one scenario per query and phase), the balanced run schedule, small statistics and the
file helpers. Stdlib only; pure except the file helpers, so the offline tests cover it.

The protocol (GPT's review of the suite, 2026-10-14, "Benchmark conclusions and a costed replacement protocol"):
48 fixed successful queries spanning resolve, traverse and pattern modes; per query one unmeasured backend warm-up,
then 6 measured cycles of all three arms (48 x (1 + 6 x 3) = 912 requests); every arm order (the 6 permutations of
three arms) once per query; a recorded random seed; one request at a time. See bench/README.md.
"""

import hashlib
import itertools
import json
import math
import os
import random
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.dirname(HERE)
if SUITE not in sys.path:
    sys.path.insert(0, SUITE)
import suite_common as sc  # noqa: E402

PANEL_VERSION = 1
ARMS = tuple(sc.ARMS)  # ("direct-local", "direct-remote", "mcp")
MODES = ("resolve", "traverse", "pattern")
ROUTE_MODE = {"POST /resolve": "resolve", "POST /traverse": "traverse", "POST /pattern": "pattern"}
DEFAULT_CYCLES = 6
DEFAULT_WARMUP_ARM = "direct-local"
# A bench scenario id: the panel query's id, '~', the phase: w0 the warm-up, c1..cN the measured cycles. One scenario
# per (query, phase) lets run_direct.py --only and the MCP runner's --ids / --repeat 0 execute exactly one request
# per step, neither runner changed.
PHASE_SEP = "~"
_BENCH_ID = re.compile(r"^(?P<query>.+)~(?P<phase>w0|c(?P<cycle>[1-9][0-9]*))$")


# ---------------------------------------------------------------------------------------------------------------
# queries: mode, family, result sizes, tags

def mode_of(route):
    """resolve / traverse / pattern for the panel's routes, None for any other."""
    return ROUTE_MODE.get(str(route or "").strip())


def family(scenario_id):
    """The scenario's family: its id without the trailing number (R-pos-0003 -> R-pos, D-T-derived-0006 ->
    D-T-derived)."""
    return re.sub(r"-[0-9]+$", "", str(scenario_id))


def _num(v):
    return v if isinstance(v, (int, float)) and not isinstance(v, bool) else None


def _count_value(obj):
    """The value of a count object ({value, relation, unit}), None when it states none."""
    return _num(obj.get("value")) if isinstance(obj, dict) else None


def _sum(values):
    vals = [v for v in values if v is not None]
    return sum(vals) if vals else None


def result_sizes(route, body):
    """What an answer RETURNS (not the work it took), per mode, for the panel's "same semantic cases and result sizes":
    two executions of a query are the same workload only when these agree. A field the answer does not state is None
    (the MCP arm's normalised answer drops some): a comparison reads only the fields both sides state.

    resolve: labels, num_kmers, runs, supported_kmers, candidates; traverse: results, paths, segments, path_bp (the
    sum of the paths' length_bp), labels (label_dict, else seed.num_labels), graphlet_lines (the graphlet output);
    pattern: patterns, returned (results retrieved), contexts / occurrences / labels / anchors (the counts' values),
    by_label."""
    if not isinstance(body, dict):
        return {}
    mode = mode_of(route)
    if mode == "resolve":
        labels = body.get("labels") if isinstance(body.get("labels"), list) else None
        return {"labels": len(labels) if labels is not None else None,
                "num_kmers": _num(body.get("num_kmers")),
                "runs": sum(len(l.get("runs") or []) for l in labels if isinstance(l, dict)) if labels is not None else None,
                "supported_kmers": _sum(_num(l.get("kmers_supported")) for l in labels if isinstance(l, dict))
                if labels is not None else None,
                "candidates": len(body["candidates"]) if isinstance(body.get("candidates"), list) else None}
    if mode == "traverse":
        results = body.get("results") if isinstance(body.get("results"), list) else None
        if results is None:
            return {"results": None}
        paths = segments = path_bp = 0
        labels, graphlet = [], []
        stated_paths = False
        for r in results:
            if not isinstance(r, dict):
                continue
            for arm in (r.get("arms") or {}).values():
                if not isinstance(arm, dict):
                    continue
                if isinstance(arm.get("paths"), list):
                    stated_paths = True
                    paths += len(arm["paths"])
                    path_bp += sum(_num(p.get("length_bp")) or 0 for p in arm["paths"] if isinstance(p, dict))
                segments += len(arm.get("segments") or [])
            if _num(r.get("graphlet_lines")) is not None:
                graphlet.append(r["graphlet_lines"])
            if isinstance(r.get("label_dict"), list):
                labels.append(len(r["label_dict"]))
            elif _num((r.get("seed") or {}).get("num_labels")) is not None:
                labels.append(r["seed"]["num_labels"])
        return {"results": len(results), "paths": paths if stated_paths else None,
                "segments": segments if stated_paths else None, "path_bp": path_bp if stated_paths else None,
                "labels": sum(labels) if labels else None, "graphlet_lines": sum(graphlet) if graphlet else None}
    if mode == "pattern":
        pats = body.get("patterns") if isinstance(body.get("patterns"), list) else None
        if pats is None:
            return {"patterns": None}
        counts = [p.get("counts") or {} for p in pats if isinstance(p, dict)]
        by_label = [len(p["by_label"]) for p in pats if isinstance(p, dict) and isinstance(p.get("by_label"), list)]
        return {"patterns": len(pats),
                "returned": _sum(_num(p.get("returned")) for p in pats if isinstance(p, dict)),
                "contexts": _sum(_count_value(c.get("contexts")) for c in counts),
                "occurrences": _sum(_count_value(c.get("occurrences")) for c in counts),
                "labels": _sum(_count_value(c.get("labels")) for c in counts),
                "anchors": _sum(_count_value(c.get("anchors")) for c in counts),
                "by_label": sum(by_label) if by_label else None}
    return {}


def answer_incomplete(route, body):
    """Why an HTTP 200 answer is not a complete, successful search ([] when it is): a per-pattern or per-result
    error (a refusal inside the answer: below the information floor, a bad alphabet), a withheld result, a stop or
    cut of any reason (a time stop depends on the host's speed; a work or row cap is budget stress, which GPT's
    protocol runs apart from the panel), a traversal arm truncated or stopped, a resolve whose labels were
    truncated. The panel takes only complete answers, and the analysis counts a measured execution with such an
    answer apart from the successes."""
    if not isinstance(body, dict):
        return ["no answer body"]
    mode = mode_of(route)
    why = []
    if body.get("error"):
        why.append("error")
    if mode == "pattern":
        for p in body.get("patterns") or []:
            if not isinstance(p, dict):
                continue
            pid = p.get("id")
            if p.get("error"):
                why.append("pattern %s: error %s" % (pid, (p["error"] or {}).get("code") if isinstance(p["error"], dict)
                                                     else p["error"]))
            if p.get("withheld"):
                why.append("pattern %s: withheld %s" % (pid, (p["withheld"] or {}).get("reason")
                                                        if isinstance(p["withheld"], dict) else p["withheld"]))
            for key in ("stop", "cut"):
                if p.get(key):
                    why.append("pattern %s: %s %s" % (pid, key, (p[key] or {}).get("reason")
                                                       if isinstance(p[key], dict) else p[key]))
    elif mode == "traverse":
        for i, r in enumerate(body.get("results") or []):
            if not isinstance(r, dict):
                continue
            if r.get("error"):
                why.append("result %d: error" % i)
            # an arm `truncated` was cut by a budget (walk_domain); `complete` with the walks `partial` reached the
            # requested scope (max_extension_bp), a requested bound and no cap
            for side, arm in sorted((r.get("arms") or {}).items()):
                if not isinstance(arm, dict):
                    continue
                if arm.get("status") not in (None, "complete"):
                    why.append("result %d %s: status %s" % (i, side, arm.get("status")))
                if arm.get("stop"):
                    why.append("result %d %s: stop" % (i, side))
    elif mode == "resolve":
        if body.get("stop"):
            why.append("stop %s" % ((body["stop"] or {}).get("reason") if isinstance(body["stop"], dict)
                                    else body["stop"]))
        if body.get("labels_truncated"):
            why.append("labels truncated")
    return why


def work_counters(route, body):
    """The work an answer states it took (context only: a release may legitimately do less): resolve none; traverse
    steps, work_units, rows_requested; pattern steps, extension_edges, annotation_rows."""
    if not isinstance(body, dict):
        return {}
    mode = mode_of(route)
    if mode == "traverse":
        steps = units = rows = None
        for r in body.get("results") or []:
            if not isinstance(r, dict):
                continue
            rows = (rows or 0) + (_num((r.get("annotation") or {}).get("rows_requested")) or 0)
            for arm in (r.get("arms") or {}).values():
                c = (arm or {}).get("counters") or {}
                if _num(c.get("steps")) is not None:
                    steps = (steps or 0) + c["steps"]
                if _num(c.get("work_units")) is not None:
                    units = (units or 0) + c["work_units"]
        return {"steps": steps, "work_units": units, "rows_requested": rows}
    if mode == "pattern":
        work = [p.get("work") or {} for p in body.get("patterns") or [] if isinstance(p, dict)]
        return {"steps": _sum(_num(w.get("steps")) for w in work),
                "extension_edges": _sum(_num(w.get("extension_edges")) for w in work),
                "annotation_rows": _sum(_num(w.get("annotation_rows")) for w in work)}
    return {}


def sizes_equal(expected, observed):
    """(equal, differing fields): the fields both state must agree; a field only one side states is not compared."""
    diff = []
    for key in sorted(set(expected or {}) | set(observed or {})):
        a, b = (expected or {}).get(key), (observed or {}).get(key)
        if a is None or b is None:
            continue
        if a != b:
            diff.append(key)
    return not diff, diff


# Tag thresholds (the panel records them): label-heavy answers name at least LABEL_HEAVY labels; a long path is a
# traversal of at least LONG_PATH_BP bp of paths, or a pattern search that extends paths beyond k (scope `long`); a
# long resolve query has at least LONG_QUERY_KMERS k-mers.
LABEL_HEAVY = 5
LONG_PATH_BP = 1000
LONG_QUERY_KMERS = 900
# A large answer (the direct answer's body, decoded) costs its delivery, the MCP arm's paging above all.
LARGE_ANSWER_BYTES = 64 * 1024
TAG_THRESHOLDS = {"label_heavy_labels": LABEL_HEAVY, "long_path_bp": LONG_PATH_BP,
                  "long_query_kmers": LONG_QUERY_KMERS, "large_answer_bytes": LARGE_ANSWER_BYTES}


def query_tags(route, body, sizes=None, bytes_body=None):
    """The query's tags from its recorded direct answer: label_heavy, long_path, long_query, large_answer,
    empty_result."""
    sizes = sizes if sizes is not None else result_sizes(route, body)
    mode = mode_of(route)
    tags = []
    if _num(bytes_body) is not None and bytes_body >= LARGE_ANSWER_BYTES:
        tags.append("large_answer")
    if mode == "resolve":
        if (sizes.get("labels") or 0) >= LABEL_HEAVY:
            tags.append("label_heavy")
        if (sizes.get("num_kmers") or 0) >= LONG_QUERY_KMERS:
            tags.append("long_query")
        if not sizes.get("runs"):
            tags.append("empty_result")
    elif mode == "traverse":
        if (sizes.get("labels") or 0) >= LABEL_HEAVY:
            tags.append("label_heavy")
        if (sizes.get("path_bp") or 0) >= LONG_PATH_BP:
            tags.append("long_path")
        if sizes.get("paths") == 0 and not sizes.get("graphlet_lines"):
            tags.append("empty_result")
    elif mode == "pattern":
        if max(sizes.get("labels") or 0, sizes.get("by_label") or 0) >= LABEL_HEAVY:
            tags.append("label_heavy")
        if any(isinstance(p, dict) and p.get("scope") == "long" for p in (body or {}).get("patterns") or []):
            tags.append("long_path")
        if not sizes.get("returned") and not sizes.get("contexts") and not sizes.get("labels"):
            tags.append("empty_result")
    return tags


def scenario_sha256(scenario):
    return hashlib.sha256(sc.canonical_json(scenario)).hexdigest()


def request_key(scenario):
    """What makes two scenarios the same request: route, request (or body_raw) and the MCP arguments."""
    mcp = scenario.get("mcp") or {}
    return hashlib.sha256(sc.canonical_json([scenario.get("route"), scenario.get("request"), scenario.get("body_raw"),
                                             mcp.get("tool"), mcp.get("arguments")])).hexdigest()


def panel_identity(queries):
    """The sha256 of the canonical JSON of the panel's [(id, scenario sha256, weight)]: two runs measure the same
    panel only when this agrees."""
    return hashlib.sha256(sc.canonical_json([[q["id"], q["scenario_sha256"], q["weight"]] for q in queries])).hexdigest()


# ---------------------------------------------------------------------------------------------------------------
# the run manifest and the schedule

def bench_id(query_id, phase):
    """'w0' (warm-up) or a cycle number >= 1 -> the bench scenario id QUERY~w0 / QUERY~cN."""
    if phase == "w0":
        return "%s%sw0" % (query_id, PHASE_SEP)
    cycle = int(phase)
    if cycle < 1:
        raise ValueError("a cycle is numbered from 1, not %r" % (phase,))
    return "%s%sc%d" % (query_id, PHASE_SEP, cycle)


def parse_bench_id(text):
    """(query id, phase 'w0' or 'cN', cycle number or 0) of a bench scenario id; None for any other id."""
    m = _BENCH_ID.match(str(text or ""))
    if not m:
        return None
    return m.group("query"), m.group("phase"), int(m.group("cycle") or 0)


def run_manifest(panel, cycles=DEFAULT_CYCLES):
    """The manifest a bench run executes: per panel query one scenario per phase (QUERY~w0, QUERY~c1..cN), each the
    query's frozen scenario with that id and `repeat` 1, no relations; content_sha256 and relations_sha256 over them
    (so both runners verify it); `bench` names the panel it derives from."""
    scenarios = []
    for q in panel["queries"]:
        for phase in ["w0"] + list(range(1, cycles + 1)):
            s = json.loads(json.dumps(q["scenario"]))
            s["id"] = bench_id(q["id"], phase)
            s["repeat"] = 1
            s.pop("derived_from", None)
            scenarios.append(s)
    src = (panel.get("source") or {}).get("manifest") or {}
    out = {"manifest_version": "bench-%s" % panel.get("panel_version", PANEL_VERSION),
           "index": src.get("index"),
           "bench": {"panel_sha256": panel.get("panel_sha256"), "panel_name": panel.get("name"), "cycles": cycles,
                     "source_manifest": src.get("content_sha256")},
           "relations": [], "scenarios": scenarios}
    out["content_sha256"] = sc.scenarios_sha256(scenarios)
    out["relations_sha256"] = sc.relations_sha256([])
    return out


def arm_orders(arms=ARMS):
    """Every order of the arms (the 6 permutations of three), in a fixed order."""
    return [list(p) for p in itertools.permutations(arms)]


def schedule(query_ids, seed, cycles=DEFAULT_CYCLES, arms=ARMS, warmup_arm=DEFAULT_WARMUP_ARM):
    """The run's steps, in run order: block 0 the warm-up pass (one unmeasured request per query through
    `warmup_arm`, queries in a seeded order), then blocks 1..cycles, one measured cycle each: every query once, in a
    fresh seeded order, its arms back to back in the order the Latin square assigns. The square: the queries in a
    seeded order get offsets 0, 1, 2, ... modulo the number of orders; query q in cycle c runs order
    sigma[(offset_q + c - 1) mod P] (sigma a seeded permutation of the P orders), so each query runs every order once
    per P cycles and each cycle runs every order equally often (P divides the panel size: 48 = 8 x 6). Each step:
    {step, block, phase, cycle, query, bench_id, arm, position (0-based in the triple), order (index into
    arm_orders), measured}."""
    rng = random.Random(seed)
    orders = arm_orders(arms)
    n_orders = len(orders)
    ranked = list(query_ids)
    rng.shuffle(ranked)
    offset = dict((q, i % n_orders) for i, q in enumerate(ranked))
    sigma = list(range(n_orders))
    rng.shuffle(sigma)
    steps = []
    warm = list(query_ids)
    rng.shuffle(warm)
    for q in warm:
        steps.append({"block": 0, "phase": "w0", "cycle": 0, "query": q, "bench_id": bench_id(q, "w0"),
                      "arm": warmup_arm, "position": 0, "order": None, "measured": False})
    for c in range(1, cycles + 1):
        order_q = list(query_ids)
        rng.shuffle(order_q)
        for q in order_q:
            oi = sigma[(offset[q] + c - 1) % n_orders]
            for pos, arm in enumerate(orders[oi]):
                steps.append({"block": c, "phase": "c%d" % c, "cycle": c, "query": q, "bench_id": bench_id(q, c),
                              "arm": arm, "position": pos, "order": oi, "measured": True})
    for i, s in enumerate(steps):
        s["step"] = i
    return steps


def schedule_balance(steps, arms=ARMS):
    """{per_query: {query: {order: count}}, per_cycle: {cycle: {order: count}}, first_arm: {arm: count}} of the
    measured steps (position 0 only, one per triple): the balance the Latin square promises."""
    per_query, per_cycle, first = {}, {}, {}
    for s in steps:
        if not s.get("measured") or s.get("position") != 0:
            continue
        per_query.setdefault(s["query"], {}).setdefault(s["order"], 0)
        per_query[s["query"]][s["order"]] += 1
        per_cycle.setdefault(s["cycle"], {}).setdefault(s["order"], 0)
        per_cycle[s["cycle"]][s["order"]] += 1
        first[s["arm"]] = first.get(s["arm"], 0) + 1
    return {"per_query": per_query, "per_cycle": per_cycle, "first_arm": first}


def requests_planned(n_queries, cycles=DEFAULT_CYCLES, arms=ARMS):
    """The measured protocol's requests: one warm-up per query, then cycles x arms per query (48 x 19 = 912)."""
    return n_queries * (1 + cycles * len(arms))


# ---------------------------------------------------------------------------------------------------------------
# statistics (no tail percentiles: six repeats support a median, not a p99)

def median(values):
    vals = sorted(v for v in values if _num(v) is not None)
    if not vals:
        return None
    n = len(vals)
    mid = n // 2
    return float(vals[mid]) if n % 2 else (vals[mid - 1] + vals[mid]) / 2.0


def quantile(values, q):
    """The q-quantile (0..1) by linear interpolation; used for bootstrap intervals and quartiles only."""
    vals = sorted(v for v in values if _num(v) is not None)
    if not vals:
        return None
    pos = (len(vals) - 1) * q
    lo, hi = int(math.floor(pos)), int(math.ceil(pos))
    return float(vals[lo]) if lo == hi else vals[lo] + (vals[hi] - vals[lo]) * (pos - lo)


def geo_mean(values, weights=None):
    """The weighted geometric mean of positive values (None when there are none)."""
    pairs = [(v, (weights[i] if weights else 1.0)) for i, v in enumerate(values) if _num(v) is not None and v > 0]
    total = sum(w for _, w in pairs)
    if not pairs or total <= 0:
        return None
    return math.exp(sum(w * math.log(v) for v, w in pairs) / total)


def weighted_mean(values, weights=None):
    pairs = [(v, (weights[i] if weights else 1.0)) for i, v in enumerate(values) if _num(v) is not None]
    total = sum(w for _, w in pairs)
    if not pairs or total <= 0:
        return None
    return sum(w * v for v, w in pairs) / total


# ---------------------------------------------------------------------------------------------------------------
# files

def read_panel(path):
    """A panel document; raises ValueError when it is none, or its stated identity does not match its queries."""
    panel = sc.read_json(path)
    if not isinstance(panel, dict) or not isinstance(panel.get("queries"), list):
        raise ValueError("%s is not a bench panel" % path)
    computed = panel_identity(panel["queries"])
    if panel.get("panel_sha256") != computed:
        raise ValueError("%s: panel_sha256 %s does not match its queries (%s)" % (path, panel.get("panel_sha256"),
                                                                                computed))
    for q in panel["queries"]:
        if scenario_sha256(q.get("scenario")) != q.get("scenario_sha256"):
            raise ValueError("%s: query %s: its scenario does not match scenario_sha256" % (path, q.get("id")))
    return panel
