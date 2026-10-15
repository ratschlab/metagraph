"""Fixtures of the bench's offline tests (test_bench_run.py, test_bench_analyze.py): a panel of fake-host queries,
a synthetic recorded bench run, deterministic walls and a finder of tail-percentile keys."""

import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.dirname(HERE)
sys.path.insert(0, SUITE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(SUITE, "bench"))

import bench_lib as bl  # noqa: E402
import suite_common as sc  # noqa: E402
from fake_mcp_bench import fake_answer  # noqa: E402

REQUESTS = {"POST /pattern": {"patterns": [{"id": "w", "dna": "ACGTACGTACGTACGTACGTACGTACGTACGTA"}]},
            "POST /resolve": {"sequence": "ACGT" * 10, "labels": ["1"]},
            "POST /traverse": {"seeds": [{"sequence": "ACGT" * 10}], "strategy": {}}}
MODES = {"POST /resolve": "resolve", "POST /traverse": "traverse", "POST /pattern": "pattern"}


def fake_panel(path, per_mode=2):
    """A panel of fake-host queries (2 per mode), each arm's expected sizes those of the fake host's answer."""
    queries = []
    for route, mode in sorted(MODES.items()):
        for i in range(per_mode):
            sid = "%s-q-%04d" % (mode[0].upper(), i)
            s = {"id": sid, "area": mode[0].upper(), "route": route, "timing_class": "fast", "title": sid,
                 "request": dict(REQUESTS[route], tag=sid) if route != "POST /pattern" else REQUESTS[route],
                 "repeat": 1, "expect": [], "tags": [], "truth": {}, "mcp": {"tool": "fake", "arguments": {}}}
            if route == "POST /pattern":
                s["request"] = {"patterns": [{"id": "w%d" % i, "dna": "ACGTACGTACGTACGTACGTACGTACGTACGTA"}]}
            sizes = bl.result_sizes(route, fake_answer(route, s["request"]))
            queries.append({"id": sid, "mode": mode, "area": s["area"], "family": bl.family(sid), "route": route,
                            "timing_class": "fast", "title": sid, "tags": [], "selected_by": "test",
                            "weight": 1.0 / (3 * per_mode), "cost_ms": 1.0,
                            "expected": dict((arm, {"http": 200, "sizes": sizes, "bytes_body": None})
                                             for arm in bl.ARMS),
                            "work": {}, "baseline_digests": {}, "recorded": {},
                            "scenario_sha256": bl.scenario_sha256(s), "scenario": s})
    panel = {"panel_version": 1, "name": "test", "source": {"manifest": {"index": {"database": "test"},
                                                                         "content_sha256": "x"}},
             "queries": queries, "costing": {"host_hours_estimate": 0.01, "allowance_host_hours": [1.5, 2.0]}}
    panel["panel_sha256"] = bl.panel_identity(queries)
    sc.write_json_atomic(path, panel)
    return panel


def tail_keys(obj):
    """The keys of a JSON document that name a tail percentile (p90, p95, p99, ...)."""
    out = []
    if isinstance(obj, dict):
        for k, v in obj.items():
            if re.match(r"^p(9[0-9]|99\.?[0-9]*)(_|$)", str(k)):
                out.append(k)
            out += tail_keys(v)
    elif isinstance(obj, list):
        for v in obj:
            out += tail_keys(v)
    return out


def synthetic_run(root, name, panel, wall_ms, seed=1, cycles=6, sra=1.0, release="r1", instances=("i1",),
                  sizes_override=None):
    """A recorded bench run without any runner: every scheduled step executed (wall_ms(query, arm, cycle)), raw
    answers with the panel's sizes, samples covering the run, identity probes; server instance per block from
    `instances` (cycled)."""
    d = os.path.join(root, name)
    os.makedirs(os.path.join(d, "samples"))
    steps = bl.schedule([q["id"] for q in panel["queries"]], seed, cycles)
    sc.write_json_atomic(os.path.join(d, "schedule.json"), {"seed": seed, "steps": steps})
    sc.write_json_atomic(os.path.join(d, "panel.json"), panel)
    sc.write_json_atomic(os.path.join(d, "bench.json"), {
        "run_id": name, "panel": {"panel_sha256": panel["panel_sha256"], "name": panel["name"]}, "seed": seed,
        "cycles": cycles, "arms": list(bl.ARMS), "guard": {}, "sessions": [{"clock_offset_start_s": 0.0}]})
    queries = dict((q["id"], q) for q in panel["queries"])
    t = 1791960000.0
    for s in steps:
        arm, q = s["arm"], queries[s["query"]]
        w = wall_ms(s["query"], arm, s["cycle"])
        raw = sc.raw_name(s["bench_id"], arm, 0)
        os.makedirs(os.path.join(d, "arms", arm, "raw"), exist_ok=True)
        body = fake_answer(q["route"], q["scenario"]["request"])
        if sizes_override and (s["query"], arm, s["cycle"]) in sizes_override:
            body = sizes_override[(s["query"], arm, s["cycle"])]
        sc.write_gz_json(os.path.join(d, "arms", arm, raw), {"answer": {"http": 200, "body": body}})
        sc.append_jsonl(os.path.join(d, "arms", arm, "executions.jsonl"), {
            "scenario": s["bench_id"], "arm": arm, "repeat": 0, "outcome": "pass", "http": 200, "time_limited": "",
            "wall_ms": w, "started_utc": sc.iso_utc(t), "ended_utc": sc.iso_utc(t + w / 1000.0), "raw": raw,
            "mcp_phases": {"backend_ms": 1.0} if arm == "mcp" else None})
        t += w / 1000.0 + 1.0
    sessions = [{"started_utc": sc.iso_utc(1791960000.0), "release_start": release,
                 "server_instance_start": instances[0]}]
    for b in range(cycles + 1):
        sessions.append({"finish": True, "ended_utc": sc.iso_utc(t), "release_end": release,
                         "server_instance_end": instances[b % len(instances)]})
        sc.append_jsonl(os.path.join(d, "steps.jsonl"), {
            "probe": "identity", "block": b, "release": release, "server_instance": instances[b % len(instances)]})
    sc.write_json_atomic(os.path.join(d, "arms", "direct-local", "run.json"), {"sessions": sessions})
    x = 1791960000.0 - 40
    while x < t + 40:
        sc.append_jsonl(os.path.join(d, "samples", "bench-x.jsonl"), {
            "t_utc": sc.iso_utc(x), "sra_all_cpu_pct": sra, "host_load1": 1.0, "cpu_pct": 1.0})
        x += 2.0
    return d


def base_wall(query, arm, cycle):
    """Deterministic walls: per query a level, mcp 10x direct-local, direct-remote 2x, a little cycle noise."""
    level = 10.0 * (1 + int(query[-1]))
    return level * {"direct-local": 1.0, "direct-remote": 2.0, "mcp": 10.0}[arm] * (1 + 0.01 * cycle)
