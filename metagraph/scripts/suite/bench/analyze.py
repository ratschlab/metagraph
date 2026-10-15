#!/usr/bin/env python3
"""Analyses a bench run (lane B), and compares two: per-query observations, paired differences and ratios between
the arms, aggregates with the panel's fixed weights, resampling intervals, the load guard, and a performance verdict
only for a clean, complete, unchanged run. Offline: no request to any host.

    analyze.py RUN_DIR [--baseline BASELINE_RUN_DIR] [--boot 2000] [--boot-seed N] [--tolerance 0.10]
               [--min-success 5] [--sra-max PCT] [--load1-max L] [--allow-missing-sra] [--no-write]

EVERY MEASURED EXECUTION (blocks 1..6; the warm-up is never measured) is classed, the first that applies:
  missing         the schedule's step has no execution line;
  transport       no HTTP status (a client timeout or a transport failure; `timeout` when it says so);
  refusal         an HTTP status other than 200 (a 4xx refusal, a 503 deadline);
  not_pass        the runner's outcome is not `pass` (fail, error, inconclusive, by design);
  time_limited    the runner or the answer states a time limitation;
  unverified      the answer could not be read;
  incomplete      the answer is no complete search (bench_lib.answer_incomplete: a refusal inside the answer, a
                  withheld result, a stop or cap);
  result_changed  a result size differs from the one the panel recorded for this arm (another workload);
  unverified      the answer does not state every result size the panel recorded for this arm (a moved or wrapped
                  body is no match), so it cannot be shown to be the same workload;
  success         otherwise.
Latency statistics use successes only; every other class is counted and listed apart with its walls (refusal and
timeout latency is never mixed into search latency).

PER QUERY and arm: the successes' walls (all six listed), median, least and greatest. PAIRED, per query and cycle,
for mcp vs direct-local, direct-remote vs direct-local and mcp vs direct-remote (both executions successes): the
difference (ms) and the ratio; per query their medians. AGGREGATES with the panel's fixed weights over the queries
that have pairs: the weighted mean of the per-query median differences and the weighted geometric mean of the
per-query median ratios, overall and per mode; per arm a latency index (the weighted geometric mean of the per-query
medians). A query without a valid observation is named, never silently re-weighted. INTERVALS: a seeded bootstrap
resampling the queries within each mode, then each resampled query's cycles (blocks); 95 % percentile intervals.
Six repeats support median comparisons: no tail percentile (p90, p99) is computed for any query or arm. The arm
order's effect is shown as each arm's median wall relative to its query median by position in the triple.

THE PERFORMANCE VERDICT is emitted only when the load guard finds every block clean (load_guard.py) with the
thresholds the run recorded (an override may only tighten them: --sra-max / --load1-max above the run's, or
--allow-missing-sra where the run did not allow it, refuse the verdict), every scheduled step is recorded, the
host's release, server instance and index fingerprint and the MCP service's capabilities are stated (every block's
identity probe, direct-local's start, every MCP invocation) and did not change during the run, and every
query has at least --min-success (5 of 6) successes on every arm. Otherwise it is refused with the reasons, and the
numbers are descriptive of an uncontrolled run.

WITH --baseline (release A vs B): both runs must emit a verdict, measure the same panel (panel_sha256), cycles and
arms; per arm and query the ratio of the medians (candidate / baseline), aggregated with the panel's weights, with a
bootstrap interval resampling queries within modes and each run's cycles independently. Per arm: `regression` when
the interval's lower end exceeds 1 + tolerance, `improvement` when its upper end is below 1 - tolerance,
`no_regression` when its upper end is at most 1 + tolerance, else `inconclusive`; the comparison passes when no arm
regresses and none is inconclusive.

Writes RUN_DIR/analysis.json and analysis.md (with --baseline: comparison.json and comparison.md in RUN_DIR).
Exit: 0 a verdict emitted (a comparison that passes), 1 a regression, 3 no verdict (refused), 4 inconclusive.
Stdlib only.
"""

import argparse
import json
import math
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bench_lib as bl  # noqa: E402
import load_guard as lg  # noqa: E402

sc = bl.sc

PAIRS = (("mcp", "direct-local"), ("direct-remote", "direct-local"), ("mcp", "direct-remote"))
CLASSES = ("success", "missing", "transport", "timeout", "refusal", "not_pass", "time_limited", "incomplete",
           "unverified", "result_changed")
DEFAULT_BOOT = 2000
DEFAULT_TOLERANCE = 0.10
DEFAULT_MIN_SUCCESS = 5
EXIT_OK, EXIT_REGRESSION, EXIT_REFUSED, EXIT_INCONCLUSIVE = 0, 1, 3, 4


def _num(v):
    return v if isinstance(v, (int, float)) and not isinstance(v, bool) else None


# ---------------------------------------------------------------------------------------------------------------
# reading a run

def load_run(run_dir):
    bench = sc.read_json(os.path.join(run_dir, "bench.json"))
    panel = sc.read_json(os.path.join(run_dir, "panel.json"))
    sched = sc.read_json(os.path.join(run_dir, "schedule.json"))
    if not isinstance(bench, dict) or not isinstance(panel, dict) or not isinstance(sched, dict):
        raise SystemExit("analyze: %s is no bench run (bench.json, panel.json, schedule.json)" % run_dir)
    panel_ok = (bl.panel_identity(panel.get("queries") or []) == panel.get("panel_sha256") ==
                (bench.get("panel") or {}).get("panel_sha256"))
    lines = {}
    for arm in bl.ARMS:
        rows, _ = sc.read_jsonl(os.path.join(run_dir, "arms", arm, "executions.jsonl"))
        for r in rows:
            if r.get("arm") == arm:
                lines[(r.get("scenario"), arm)] = r  # the last line of an execution recorded twice
    return {"dir": os.path.abspath(run_dir), "bench": bench, "panel": panel, "panel_ok": panel_ok,
            "steps": sched.get("steps") or [],
            "seed": sched.get("seed"), "lines": lines,
            "queries": dict((q["id"], q) for q in panel.get("queries") or [])}


def answer_body(run_dir, line):
    rel = line.get("raw")
    if not rel:
        return None
    arm = line.get("arm") or ""
    for path in (os.path.join(run_dir, "arms", arm, rel), os.path.join(run_dir, rel)):
        if os.path.exists(path):
            try:
                doc = sc.read_gz_json(path)
            except (OSError, ValueError, EOFError):
                return None
            answer = (doc or {}).get("answer")
            return answer.get("body") if isinstance(answer, dict) else None
    return None


def classify(run, step, line):
    """(class, detail, sizes) of a step's execution."""
    if line is None:
        return "missing", "no execution recorded", None
    http = line.get("http")
    if not isinstance(http, int) or isinstance(http, bool) or http <= 0:
        err = line.get("error") or {}
        text = " ".join(str(err.get(k) or "") for k in ("code", "text")) if isinstance(err, dict) else str(err)
        return ("timeout" if ("timeout" in text.lower() or "deadline" in text.lower()) else "transport"), \
            text.strip() or "no HTTP status", None
    if http != 200:
        return "refusal", "http %s" % http, None
    if line.get("outcome") != "pass":
        return "not_pass", "outcome %s" % line.get("outcome"), None
    if line.get("time_limited"):
        return "time_limited", str(line.get("time_limited")), None
    q = run["queries"].get(step["query"]) or {}
    body = answer_body(run["dir"], line)
    if body is None:
        return "unverified", "the answer could not be read", None
    why = bl.answer_incomplete(q.get("route"), body)
    accepted = ((q.get("accepted_incomplete") or {}).get("arms") or {}).get(step["arm"])
    if why and sorted(why) != sorted(accepted or []):
        only = " (the panel accepts only: %s)" % "; ".join(accepted) if accepted else ""
        return "incomplete", "; ".join(why[:3]) + only, None
    if accepted and not why:
        # the panel recorded a capped answer (a stated exception): an uncapped one is another workload
        return "result_changed", "not capped as the panel recorded (%s)" % "; ".join(accepted), None
    sizes = bl.result_sizes(q.get("route"), body)
    expected = ((q.get("expected") or {}).get(step["arm"]) or {}).get("sizes")
    if not expected or all(v is None for v in expected.values()):
        return "unverified", "the panel recorded no sizes for this arm", sizes
    ok, diff, unstated = bl.sizes_equal(expected, sizes)
    if diff:
        return "result_changed", "sizes differ: %s" % ", ".join(
            "%s %s -> %s" % (k, expected.get(k), sizes.get(k)) for k in diff), sizes
    if unstated:
        return "unverified", "the answer does not state %s, which the panel recorded" % ", ".join(unstated), sizes
    return "success", "", sizes


def observations(run):
    """One record per measured step: query, mode, cycle, arm, position, order, class, detail, wall and phases."""
    out = []
    for step in run["steps"]:
        if not step.get("measured"):
            continue
        line = run["lines"].get((step["bench_id"], step["arm"]))
        cls, detail, sizes = classify(run, step, line)
        q = run["queries"].get(step["query"]) or {}
        rec = {"query": step["query"], "mode": q.get("mode"), "cycle": step["cycle"], "block": step["block"],
               "arm": step["arm"], "position": step["position"], "order": step["order"], "class": cls,
               "detail": detail, "wall_ms": _num((line or {}).get("wall_ms"))}
        if cls == "success":
            rec["sizes"] = dict((k, v) for k, v in (sizes or {}).items() if v is not None)
        if line is not None:
            server = line.get("server") if isinstance(line.get("server"), dict) else {}
            rec["server_elapsed_ms"] = _num(server.get("elapsed_ms"))
            if step["arm"] == "mcp":
                rec["mcp_phases"] = dict((k, _num(v)) for k, v in (line.get("mcp_phases") or {}).items())
            else:
                client = line.get("client") or {}
                rec["ttfb_ms"] = _num(client.get("ttfb_ms"))
                rec["connect_ms"] = _num(client.get("connect_ms"))
        out.append(rec)
    return out


# ---------------------------------------------------------------------------------------------------------------
# identity: the host's release / instance and the MCP service through the run

HOST_IDENTITY = ("release", "server_instance", "index_fp")


def identity(run, arms=bl.ARMS):
    """The deployment through the run: the host's release, server instance and index fingerprint (the direct arms'
    run.json sessions, the identity probe at each block's end in steps.jsonl) and the MCP service (mcp_runs.jsonl).
    `changed`: a value that differs; `blocks_changed`: the blocks at whose end a value differs from the block before;
    `unknown`: the evidence missing (a block whose probe states no release and server instance, no start identity of
    direct-local, an MCP invocation stating no service capabilities): unknown is never read as unchanged."""
    out = {"host": [], "mcp": [], "probes": [], "changed": [], "blocks_changed": [], "unknown": []}
    starts = 0
    for arm in ("direct-local", "direct-remote"):
        doc = sc.read_json(os.path.join(run["dir"], "arms", arm, "run.json")) or {}
        for sess in doc.get("sessions") or []:
            for when in ("start", "end"):
                vals = dict((k, sess.get("%s_%s" % (k, when))) for k in HOST_IDENTITY)
                if vals["release"] or vals["server_instance"]:
                    at = sess.get("started_utc" if when == "start" else "ended_utc")
                    out["host"].append(dict(vals, arm=arm, at=at))
                if arm == "direct-local" and when == "start" and vals["release"] and vals["server_instance"]:
                    starts += 1
    rows, _ = sc.read_jsonl(os.path.join(run["dir"], "arms", "mcp", "mcp_runs.jsonl"))
    for r in rows:
        caps = r.get("service_capabilities") if isinstance(r.get("service_capabilities"), dict) else {}
        out["mcp"].append({"at": r.get("started_utc"), "service_capabilities_sha256": caps.get("sha256"),
                           "runner_git_head": r.get("git_head"), "checks_sha256": r.get("checks_sha256")})
    rows, _ = sc.read_jsonl(os.path.join(run["dir"], "steps.jsonl"))
    out["probes"] = [dict((k, r.get(k)) for k in ("block", "rc") + HOST_IDENTITY) for r in rows
                     if r.get("probe") == "identity"]
    for key in HOST_IDENTITY:
        vals = sorted(set(x[key] for x in out["host"] + out["probes"] if x.get(key)))
        if len(vals) > 1:
            out["changed"].append("%s: %s" % (key, ", ".join(vals)))
    for key in ("service_capabilities_sha256", "runner_git_head", "checks_sha256"):
        vals = sorted(set(x[key] for x in out["mcp"] if x.get(key)))
        if len(vals) > 1:
            out["changed"].append("mcp %s: %d values" % (key, len(vals)))
    # the probes at the blocks' ends: the block whose end differs from the end of the block before
    prev = {}
    for p in out["probes"]:
        if not (p.get("release") and p.get("server_instance")):
            continue
        if any(p.get(k) and prev.get(k) and p[k] != prev[k] for k in HOST_IDENTITY):
            out["blocks_changed"].append(p["block"])
        prev = dict((k, p.get(k) or prev.get(k)) for k in HOST_IDENTITY)
    blocks = sorted(set(s["block"] for s in run["steps"]))
    known = set(p["block"] for p in out["probes"] if p.get("release") and p.get("server_instance"))
    blind = [b for b in blocks if b not in known]
    if blind:
        out["unknown"].append("no identity probe states the host's release and server instance at the end of "
                              "blocks %s" % ", ".join(str(b) for b in blind))
    if not starts:
        out["unknown"].append("direct-local states no start identity (run.json)")
    if "mcp" in arms:
        if not out["mcp"]:
            out["unknown"].append("the MCP runner recorded no service identity (mcp_runs.jsonl)")
        else:
            n = sum(1 for x in out["mcp"] if not x["service_capabilities_sha256"])
            if n:
                out["unknown"].append("%d MCP invocation(s) state no service capabilities sha256" % n)
    return out


# ---------------------------------------------------------------------------------------------------------------
# statistics

def per_query(obs, arms):
    """{query: {arm: {walls by cycle (successes), classes}}}."""
    out = {}
    for o in obs:
        q = out.setdefault(o["query"], {"mode": o["mode"], "arms": {}})
        a = q["arms"].setdefault(o["arm"], {"walls": {}, "classes": {}, "others": [], "sizes": set()})
        a["classes"][o["class"]] = a["classes"].get(o["class"], 0) + 1
        if o["class"] == "success":
            a["sizes"].add(json.dumps(o.get("sizes") or {}, sort_keys=True))
        if o["class"] == "success" and o["wall_ms"] is not None:
            a["walls"][o["cycle"]] = o["wall_ms"]
        elif o["class"] != "success":
            a["others"].append({"cycle": o["cycle"], "class": o["class"], "detail": o["detail"],
                                "wall_ms": o["wall_ms"]})
    return out


def pairs_of(q, a, b):
    """{cycle: (wall a, wall b)} of a query's cycles where both arms succeeded."""
    wa, wb = (q["arms"].get(a) or {}).get("walls") or {}, (q["arms"].get(b) or {}).get("walls") or {}
    return dict((c, (wa[c], wb[c])) for c in sorted(set(wa) & set(wb)) if wa[c] > 0 and wb[c] > 0)


def pair_stats(pairs):
    if not pairs:
        return None
    diffs = [x - y for x, y in pairs.values()]
    logs = [math.log(x / y) for x, y in pairs.values()]
    return {"n": len(pairs), "median_diff_ms": bl.median(diffs), "median_ratio": math.exp(bl.median(logs)),
            "diffs_ms": [round(d, 1) for d in diffs], "ratios": [round(math.exp(v), 4) for v in logs]}


def summarise_query(q, arms):
    out = {"mode": q["mode"], "arms": {}, "pairs": {}}
    for arm in arms:
        a = q["arms"].get(arm) or {"walls": {}, "classes": {}, "others": [], "sizes": set()}
        walls = [a["walls"][c] for c in sorted(a["walls"])]
        out["arms"][arm] = {"n_success": len(walls), "walls_ms": [round(w, 1) for w in walls],
                            "median_ms": bl.median(walls), "min_ms": min(walls) if walls else None,
                            "max_ms": max(walls) if walls else None, "classes": a["classes"], "others": a["others"],
                            "sizes": [json.loads(x) for x in sorted(a.get("sizes") or ())]}
    for a, b in PAIRS:
        if a in arms and b in arms:
            out["pairs"]["%s/%s" % (a, b)] = pair_stats(pairs_of(q, a, b))
    return out


def weighted_aggregates(summaries, weights, arms, modes=None):
    """The pair aggregates (mean difference, geometric mean ratio) and the arms' latency index over the queries."""
    sel = [q for q in sorted(summaries) if modes is None or summaries[q]["mode"] in modes]
    out = {"pairs": {}, "arms": {}, "queries": len(sel)}
    for a, b in PAIRS:
        key = "%s/%s" % (a, b)
        if a not in arms or b not in arms:
            continue
        have = [q for q in sel if summaries[q]["pairs"].get(key)]
        w = [weights[q] for q in have]
        out["pairs"][key] = {
            "queries": len(have), "missing": [q for q in sel if q not in have],
            "mean_diff_ms": bl.weighted_mean([summaries[q]["pairs"][key]["median_diff_ms"] for q in have], w),
            "geo_mean_ratio": bl.geo_mean([summaries[q]["pairs"][key]["median_ratio"] for q in have], w)}
    for arm in arms:
        have = [q for q in sel if summaries[q]["arms"][arm]["median_ms"]]
        out["arms"][arm] = {"queries": len(have), "missing": [q for q in sel if q not in have],
                            "latency_index_ms": bl.geo_mean([summaries[q]["arms"][arm]["median_ms"] for q in have],
                                                            [weights[q] for q in have])}
    return out


def bootstrap(per, weights, arms, boot, seed):
    """95 % intervals of the pair aggregates and the arms' latency indices: queries resampled within each mode, then
    each query's cycles."""
    rng = random.Random(seed)
    strata = {}
    for q in sorted(per):
        strata.setdefault(per[q]["mode"], []).append(q)
    stats = dict(("%s/%s" % p, {"diff": [], "ratio": []}) for p in PAIRS if p[0] in arms and p[1] in arms)
    idx = dict((arm, []) for arm in arms)
    for _ in range(boot):
        sample = []
        for mode in sorted(strata):
            qs = strata[mode]
            sample += [qs[rng.randrange(len(qs))] for _ in qs]
        for a, b in PAIRS:
            key = "%s/%s" % (a, b)
            if key not in stats:
                continue
            diffs, logs, w = [], [], []
            for q in sample:
                pairs = list(pairs_of(per[q], a, b).values())
                if not pairs:
                    continue
                draw = [pairs[rng.randrange(len(pairs))] for _ in pairs]
                diffs.append(bl.median([x - y for x, y in draw]))
                logs.append(bl.median([math.log(x / y) for x, y in draw]))
                w.append(weights[q])
            if w:
                stats[key]["diff"].append(bl.weighted_mean(diffs, w))
                stats[key]["ratio"].append(math.exp(bl.weighted_mean(logs, w)))
        for arm in arms:
            meds, w = [], []
            for q in sample:
                walls = list(((per[q]["arms"].get(arm) or {}).get("walls") or {}).values())
                if not walls:
                    continue
                meds.append(bl.median([walls[rng.randrange(len(walls))] for _ in walls]))
                w.append(weights[q])
            if w:
                idx[arm].append(bl.geo_mean(meds, w))
    ci = lambda v: [bl.quantile(v, 0.025), bl.quantile(v, 0.975)] if v else None  # noqa: E731
    return {"boot": boot, "seed": seed,
            "pairs": dict((k, {"mean_diff_ms": ci(v["diff"]), "geo_mean_ratio": ci(v["ratio"])})
                          for k, v in stats.items()),
            "arms": dict((arm, {"latency_index_ms": ci(v)}) for arm, v in idx.items())}


def order_effect(obs, summaries):
    """Each arm's median of wall / its query's median by position in the triple (balanced orders cancel it)."""
    out = {}
    for o in obs:
        if o["class"] != "success" or o["wall_ms"] is None:
            continue
        med = summaries[o["query"]]["arms"][o["arm"]]["median_ms"]
        if med:
            out.setdefault(o["arm"], {}).setdefault(o["position"], []).append(o["wall_ms"] / med)
    return dict((arm, dict((str(p), {"n": len(v), "median_rel": bl.median(v)}) for p, v in sorted(pos.items())))
                for arm, pos in out.items())


# ---------------------------------------------------------------------------------------------------------------
# a run's analysis

def analyse(run_dir, boot=DEFAULT_BOOT, boot_seed=None, min_success=DEFAULT_MIN_SUCCESS, guard_overrides=None):
    run = load_run(run_dir)
    arms = list(run["bench"].get("arms") or bl.ARMS)
    weights = dict((q["id"], q["weight"]) for q in run["panel"]["queries"])
    obs = observations(run)
    per = per_query(obs, arms)
    summaries = dict((q, summarise_query(per[q], arms)) for q in sorted(per))
    classes = {}
    for o in obs:
        classes.setdefault(o["arm"], {}).setdefault(o["class"], 0)
        classes[o["arm"]][o["class"]] += 1
    guard = lg.judge_run(run_dir, guard_overrides)
    ident = identity(run, arms)
    scheduled = len(run["steps"])
    recorded = sum(1 for s in run["steps"] if (s["bench_id"], s["arm"]) in run["lines"])
    short = sorted("%s %s (%d of %d)" % (q, arm, summaries[q]["arms"][arm]["n_success"], run["bench"].get("cycles"))
                   for q in summaries for arm in arms if summaries[q]["arms"][arm]["n_success"] < min_success)
    absent = sorted(set(weights) - set(summaries))
    refused = []
    if not run["panel_ok"]:
        refused.append("the run's panel.json does not match its identity (panel_sha256, bench.json)")
    if guard["loosened"]:
        refused.append("thresholds loosened after the run: %s" % "; ".join(guard["loosened"]))
    if guard["verdict"] != "clean":
        refused.append("load guard: %s (loaded blocks %s, unguarded blocks %s)" % (
            guard["verdict"], ",".join(guard["loaded"]) or "-", ",".join(guard["unguarded"]) or "-"))
    if recorded < scheduled:
        refused.append("incomplete: %d of %d scheduled steps recorded" % (recorded, scheduled))
    if ident["changed"]:
        refused.append("the deployment changed during the run: %s" % "; ".join(ident["changed"]))
    if ident["unknown"]:
        refused.append("the deployment's identity is unknown: %s" % "; ".join(ident["unknown"]))
    if short:
        refused.append("%d query/arm cells with fewer than %d successes: %s" % (
            len(short), min_success, ", ".join(short[:8])))
    if absent:
        refused.append("panel queries never measured: %s" % ", ".join(absent[:8]))
    boot_seed = boot_seed if boot_seed is not None else (run["seed"] or 0)
    result = {
        "run_id": run["bench"].get("run_id"), "run_dir": run["dir"], "panel": run["bench"].get("panel"),
        "seed": run["seed"], "cycles": run["bench"].get("cycles"), "arms": arms,
        "steps": {"scheduled": scheduled, "recorded": recorded},
        "verdict": {"emitted": not refused, "refused": refused,
                    "statement": ("a clean, complete run without a deployment change: its medians and ratios are "
                                  "this release's performance baseline for the panel") if not refused else
                                 ("no performance verdict: the numbers below describe an uncontrolled run")},
        "min_success": min_success, "load_guard": guard, "identity": ident, "classes": classes,
        "aggregates": {"all": weighted_aggregates(summaries, weights, arms)},
        "intervals": bootstrap(per, weights, arms, boot, boot_seed) if boot else None,
        "order_effect": order_effect(obs, summaries),
        "queries": summaries,
        "note": "six repeats per query and arm support median comparisons; no tail percentile (p90, p99) is "
                "computed. Weights are the panel's, fixed; a query without a valid observation is listed under "
                "`missing`, never re-weighted silently. Non-successes are counted apart with their walls.",
    }
    for mode in bl.MODES:
        result["aggregates"][mode] = weighted_aggregates(summaries, weights, arms, (mode,))
    return result, per


# ---------------------------------------------------------------------------------------------------------------
# two runs

def compare(cand_dir, base_dir, boot=DEFAULT_BOOT, boot_seed=None, tolerance=DEFAULT_TOLERANCE,
            min_success=DEFAULT_MIN_SUCCESS, guard_overrides=None):
    cand, cper = analyse(cand_dir, 0, boot_seed, min_success, guard_overrides)
    base, bper = analyse(base_dir, 0, boot_seed, min_success, guard_overrides)
    refused = []
    for name, res in (("candidate", cand), ("baseline", base)):
        if not res["verdict"]["emitted"]:
            refused.append("%s: no verdict (%s)" % (name, "; ".join(res["verdict"]["refused"])))
    if (cand["panel"] or {}).get("panel_sha256") != (base["panel"] or {}).get("panel_sha256"):
        refused.append("the runs measured different panels")
    if cand["cycles"] != base["cycles"] or cand["arms"] != base["arms"]:
        refused.append("the runs differ in cycles or arms")
    # the same results: every query's successful answers state the same result sizes in both runs (the panel binds
    # the recorded ones; a size it did not record, or a release that changed a result, shows here)
    differ = []
    for q in sorted(set(cand["queries"]) & set(base["queries"])):
        for arm in [a for a in cand["arms"] if a in base["arms"]]:
            sc_, sb = cand["queries"][q]["arms"][arm]["sizes"], base["queries"][q]["arms"][arm]["sizes"]
            if sc_ and sb and sc_ != sb:
                differ.append("%s %s (%s vs %s)" % (q, arm, json.dumps(sc_, sort_keys=True),
                                                    json.dumps(sb, sort_keys=True)))
    if differ:
        refused.append("result sizes differ between the runs: %d query/arm cell(s): %s" % (
            len(differ), "; ".join(differ[:6])))
    panel = sc.read_json(os.path.join(cand_dir, "panel.json")) or {}
    weights = dict((q["id"], q["weight"]) for q in panel.get("queries") or [])
    arms = [a for a in cand["arms"] if a in base["arms"]]
    rng = random.Random(boot_seed if boot_seed is not None else 0)
    strata = {}
    for q in sorted(set(cper) & set(bper)):
        strata.setdefault(cper[q]["mode"], []).append(q)
    per_arm = {}
    for arm in arms:
        ratios, w, qs = [], [], []
        for q in sorted(set(cper) & set(bper)):
            mc = cand["queries"][q]["arms"][arm]["median_ms"]
            mb = base["queries"][q]["arms"][arm]["median_ms"]
            if mc and mb:
                ratios.append(mc / mb)
                w.append(weights.get(q, 0.0))
                qs.append(q)
        boots = []
        for _ in range(boot):
            logs, ww = [], []
            for mode in sorted(strata):
                for _ in strata[mode]:
                    q = strata[mode][rng.randrange(len(strata[mode]))]
                    wc = list(((cper[q]["arms"].get(arm) or {}).get("walls") or {}).values())
                    wb = list(((bper[q]["arms"].get(arm) or {}).get("walls") or {}).values())
                    if not wc or not wb:
                        continue
                    mc = bl.median([wc[rng.randrange(len(wc))] for _ in wc])
                    mb = bl.median([wb[rng.randrange(len(wb))] for _ in wb])
                    logs.append(math.log(mc / mb))
                    ww.append(weights.get(q, 0.0))
            if ww:
                boots.append(math.exp(bl.weighted_mean(logs, ww)))
        g = bl.geo_mean(ratios, w)
        ci = [bl.quantile(boots, 0.025), bl.quantile(boots, 0.975)] if boots else None
        if ci is None or g is None:
            verdict = "inconclusive"
        elif ci[0] > 1 + tolerance:
            verdict = "regression"
        elif ci[1] < 1 - tolerance:
            verdict = "improvement"
        elif ci[1] <= 1 + tolerance:
            verdict = "no_regression"
        else:
            verdict = "inconclusive"
        per_arm[arm] = {"geo_mean_ratio": g, "interval": ci, "verdict": verdict, "queries": len(qs),
                        "per_query": dict((q, round(r, 4)) for q, r in zip(qs, ratios))}
    if refused:
        overall = "refused"
    elif any(v["verdict"] == "regression" for v in per_arm.values()):
        overall = "regression"
    elif any(v["verdict"] == "inconclusive" for v in per_arm.values()):
        overall = "inconclusive"
    else:
        overall = "pass"
    return {"candidate": cand["run_id"], "baseline": base["run_id"], "candidate_dir": cand["run_dir"],
            "baseline_dir": base["run_dir"], "tolerance": tolerance, "boot": boot, "boot_seed": boot_seed,
            "guard": {"overrides": dict((k, v) for k, v in (guard_overrides or {}).items() if v is not None),
                      "candidate": cand["load_guard"]["limits"], "baseline": base["load_guard"]["limits"],
                      "candidate_loosened": cand["load_guard"]["loosened"],
                      "baseline_loosened": base["load_guard"]["loosened"]},
            "verdict": overall, "refused": refused, "arms": per_arm,
            "note": "per arm the ratio of the per-query medians (candidate / baseline), aggregated with the panel's "
                    "weights; the interval resamples queries within modes and each run's cycles independently"}


# ---------------------------------------------------------------------------------------------------------------
# reports

def fmt(v, digits=1):
    if v is None:
        return "-"
    if isinstance(v, float):
        return ("%%.%df" % digits) % v
    return str(v)


def fmt_ci(ci, digits=1):
    return "-" if not ci or ci[0] is None else "[%s, %s]" % (fmt(ci[0], digits), fmt(ci[1], digits))


def report_md(res):
    out = ["# Bench run %s" % res["run_id"], ""]
    v = res["verdict"]
    out.append("**Performance verdict: %s.** %s" % ("emitted" if v["emitted"] else "refused", v["statement"]))
    for r in v["refused"]:
        out.append("- %s" % r)
    out += ["", "Panel %s (%s), seed %s, %s cycles, arms %s; %d of %d scheduled steps recorded." % (
        (res["panel"] or {}).get("name"), str((res["panel"] or {}).get("panel_sha256"))[:12], res["seed"],
        res["cycles"], ", ".join(res["arms"]), res["steps"]["recorded"], res["steps"]["scheduled"]), ""]
    g = res["load_guard"]
    out += ["## Load guard: %s" % g["verdict"], "",
            "Thresholds: sra_all CPU <= %s %%, load1 <= %s, coverage gaps <= %s s, windows padded %s s%s." % (
                g["limits"]["sra_all_cpu_pct_max"], g["limits"]["host_load1_max"], g["limits"]["max_gap_s"],
                g["limits"]["pad_s"], ", a missing sra_all reading counts 0 %" if g["limits"]["allow_missing_sra"]
                else ""),
            "Overrides: %s%s." % (json.dumps(g.get("overrides") or {}, sort_keys=True),
                                  "; loosened after the run (no verdict): %s" % "; ".join(g["loosened"])
                                  if g.get("loosened") else ""),
            "", "| block | verdict | samples | sra_all max % | load1 max | reason |",
            "|---|---|---|---|---|---|"]
    for b, x in g["blocks"].items():
        out.append("| %s | %s | %d | %s | %s | %s |" % (
            b, x["verdict"], x["samples"], fmt(x["sra_all_cpu_pct_max"]), fmt(x["host_load1_max"], 2), x["reason"]))
    ident = res["identity"]
    out += ["", "## Identity", "", "Host: %s. MCP service: %s." % (
        ", ".join(sorted(set("%s/%s" % (x["release"], x["server_instance"]) for x in ident["host"]))) or "unknown",
        "%d invocation(s)" % len(ident["mcp"])) + (" Changed: %s." % "; ".join(ident["changed"])
                                                   if ident["changed"] else " Unchanged.") +
            (" Unknown: %s." % "; ".join(ident["unknown"]) if ident["unknown"] else ""), ""]
    out += ["## Executions by class (measured cycles)", "", "| arm | " + " | ".join(CLASSES) + " |",
            "|---|" + "---|" * len(CLASSES)]
    for arm in res["arms"]:
        c = res["classes"].get(arm) or {}
        out.append("| %s | %s |" % (arm, " | ".join(str(c.get(k, 0)) for k in CLASSES)))
    others = [(q, arm, o) for q, s in res["queries"].items() for arm, a in s["arms"].items() for o in a["others"]]
    if others:
        out += ["", "Non-successes (their walls are not part of any latency statistic):", ""]
        for q, arm, o in others[:60]:
            out.append("- %s %s c%s: %s (%s), wall %s ms" % (
                q, arm, o["cycle"], o["class"], o["detail"], fmt(o["wall_ms"], 0)))
        if len(others) > 60:
            out.append("- ... %d more" % (len(others) - 60))
    ints = res.get("intervals") or {"pairs": {}, "arms": {}}
    out += ["", "## Aggregates (panel weights; 95 %% bootstrap intervals, %s resamples)" % (
        (res.get("intervals") or {}).get("boot")), "", "| scope | pair | queries | mean difference ms | interval | "
            "geometric mean ratio | interval |", "|---|---|---|---|---|---|---|"]
    for scope in ("all",) + bl.MODES:
        agg = res["aggregates"][scope]
        for key, p in agg["pairs"].items():
            ci = ints["pairs"].get(key, {}) if scope == "all" else {}
            out.append("| %s | %s | %d | %s | %s | %s | %s |" % (
                scope, key, p["queries"], fmt(p["mean_diff_ms"]), fmt_ci(ci.get("mean_diff_ms")),
                fmt(p["geo_mean_ratio"], 3), fmt_ci(ci.get("geo_mean_ratio"), 3)))
    out += ["", "| scope | arm | queries | latency index ms (geometric mean of medians) | interval |",
            "|---|---|---|---|---|"]
    for scope in ("all",) + bl.MODES:
        for arm, a in res["aggregates"][scope]["arms"].items():
            ci = ints["arms"].get(arm, {}) if scope == "all" else {}
            out.append("| %s | %s | %d | %s | %s |" % (scope, arm, a["queries"], fmt(a["latency_index_ms"]),
                                                       fmt_ci(ci.get("latency_index_ms"))))
    out += ["", "## Arm order", "", "Median wall relative to the query's median, by position in the triple "
            "(the balanced orders cancel this in the aggregates):", ""]
    for arm, pos in sorted(res["order_effect"].items()):
        out.append("- %s: %s" % (arm, ", ".join(
            "position %s %s (n %d)" % (p, fmt(x["median_rel"], 3), x["n"]) for p, x in pos.items())))
    out += ["", "## Per query", "",
            "| query | mode | " + " | ".join("%s median ms [min, max] n" % a for a in res["arms"])
            + " | " + " | ".join("%s ratio" % k for k in ("mcp/direct-local", "direct-remote/direct-local")) + " |",
            "|---|---|" + "---|" * (len(res["arms"]) + 2)]
    for q, s in sorted(res["queries"].items(), key=lambda kv: (kv[1]["mode"] or "", kv[0])):
        cells = []
        for arm in res["arms"]:
            a = s["arms"][arm]
            cells.append("%s [%s, %s] %d" % (fmt(a["median_ms"]), fmt(a["min_ms"]), fmt(a["max_ms"]), a["n_success"]))
        ratios = [fmt((s["pairs"].get(k) or {}).get("median_ratio"), 3)
                  for k in ("mcp/direct-local", "direct-remote/direct-local")]
        out.append("| %s | %s | %s | %s |" % (q, s["mode"], " | ".join(cells), " | ".join(ratios)))
    out += ["", res["note"], ""]
    return "\n".join(out)


def comparison_md(cmp_):
    out = ["# Bench comparison: %s vs %s" % (cmp_["candidate"], cmp_["baseline"]), "",
           "**Verdict: %s** (tolerance %s; candidate / baseline per arm)." % (cmp_["verdict"], cmp_["tolerance"]),
           "", "Load-guard overrides: %s." % json.dumps(cmp_["guard"]["overrides"], sort_keys=True), ""]
    for r in cmp_["refused"]:
        out.append("- %s" % r)
    out += ["", "| arm | queries | geometric mean ratio | 95 % interval | verdict |", "|---|---|---|---|---|"]
    for arm, a in cmp_["arms"].items():
        out.append("| %s | %d | %s | %s | %s |" % (arm, a["queries"], fmt(a["geo_mean_ratio"], 3),
                                                   fmt_ci(a["interval"], 3), a["verdict"]))
    out += ["", cmp_["note"], ""]
    return "\n".join(out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("run_dir")
    ap.add_argument("--baseline", default=None, help="another bench run to compare with (release A vs B)")
    ap.add_argument("--boot", type=int, default=DEFAULT_BOOT)
    ap.add_argument("--boot-seed", type=int, default=None)
    ap.add_argument("--tolerance", type=float, default=DEFAULT_TOLERANCE)
    ap.add_argument("--min-success", type=int, default=DEFAULT_MIN_SUCCESS)
    ap.add_argument("--sra-max", type=float, default=None)
    ap.add_argument("--load1-max", type=float, default=None)
    ap.add_argument("--allow-missing-sra", action="store_true", default=None)
    ap.add_argument("--no-write", action="store_true")
    args = ap.parse_args(argv)
    overrides = {"sra_all_cpu_pct_max": args.sra_max, "host_load1_max": args.load1_max,
                 "allow_missing_sra": args.allow_missing_sra}
    if args.baseline:
        res = compare(args.run_dir, args.baseline, args.boot, args.boot_seed, args.tolerance, args.min_success,
                      overrides)
        if not args.no_write:
            sc.write_json_atomic(os.path.join(args.run_dir, "comparison.json"), res)
            with open(os.path.join(args.run_dir, "comparison.md"), "w") as fh:
                fh.write(comparison_md(res))
        print("comparison %s vs %s: %s" % (res["candidate"], res["baseline"], res["verdict"]))
        return {"pass": EXIT_OK, "regression": EXIT_REGRESSION, "inconclusive": EXIT_INCONCLUSIVE}.get(
            res["verdict"], EXIT_REFUSED)
    res, _ = analyse(args.run_dir, args.boot, args.boot_seed, args.min_success, overrides)
    res["guard_overrides"] = dict((k, v) for k, v in overrides.items() if v is not None)
    if not args.no_write:
        sc.write_json_atomic(os.path.join(args.run_dir, "analysis.json"), res)
        with open(os.path.join(args.run_dir, "analysis.md"), "w") as fh:
            fh.write(report_md(res))
    print("analysis %s: verdict %s%s" % (res["run_id"], "emitted" if res["verdict"]["emitted"] else "refused",
                                         "" if res["verdict"]["emitted"] else " (%s)" % "; ".join(
                                             res["verdict"]["refused"])))
    return EXIT_OK if res["verdict"]["emitted"] else EXIT_REFUSED


if __name__ == "__main__":
    sys.exit(main())
