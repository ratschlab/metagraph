#!/usr/bin/env python3
"""A stand-in for the async session's MCP runner for the bench's offline tests (bench/run_bench.py): the per-id CLI
of run_suite.sh's template (--ids A --repeat N --out DIR --manifest M [--checks-dir DIR]). Per id one execution line
(outcome pass, http 200, mcp_phases, its wall FAKE_BENCH_MCP_MS ms, default 20) and its raw answer
raw/<id>.mcp.<repeat>.json.gz holding the answer the fake host (tests/fake_host.py) gives the scenario's route, so
the analysis reads the same result sizes; per invocation one mcp_runs.jsonl line with the service's capabilities
sha256 (FAKE_BENCH_MCP_CAPS, default `caps-1`). FAKE_BENCH_MCP_ERROR: ids that end `error` without an HTTP status
(the runner's poll deadline), exit 1. Skips what is recorded; logs each call to DIR/calls.log. No request is made."""

import argparse
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import suite_common as sc  # noqa: E402


def fake_answer(route, request):
    """The fake host's answer body for a route (tests/fake_host.py)."""
    if route == "POST /pattern":
        return {"patterns": [{"id": p.get("id"), "counts": {"contexts": {"value": 1, "relation": "exact"}},
                              "timing": {"elapsed_ms": 1.5}, "work": {"steps": 10, "ranges_visited": 10,
                                                                      "memory_bytes": 1000}}
                             for p in (request or {}).get("patterns", [])], "timing": {"elapsed_ms": 2.0}}
    if route == "POST /resolve":
        return {"k": 31, "num_kmers": 30, "labels": [], "timing": {"elapsed_ms": 3.0}}
    if route == "POST /traverse":
        return {"results": [{"arms": {"right": {"counters": {"steps": 5}}, "left": {"counters": {"steps": 7}}},
                             "annotation": {"rows_requested": 4}, "timing": {"walk_ms": 1.0}}]}
    return {}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--ids", required=True)
    ap.add_argument("--repeat", type=int, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--checks-dir", default=None)
    args = ap.parse_args(argv)
    os.makedirs(os.path.join(args.out, "raw"), exist_ok=True)
    manifest = sc.read_json(args.manifest)
    scen = dict((s.get("id"), s) for s in manifest.get("scenarios", []))
    errors = set(x for x in os.environ.get("FAKE_BENCH_MCP_ERROR", "").split(",") if x)
    wall = float(os.environ.get("FAKE_BENCH_MCP_MS", "20"))
    path = os.path.join(args.out, "executions.jsonl")
    done = set(sc.execution_key(r) for r in sc.read_jsonl(path)[0])
    with open(os.path.join(args.out, "calls.log"), "a") as fh:
        fh.write("%s %d checks_dir=%s\n" % (args.ids, args.repeat, args.checks_dir))
    sc.append_jsonl(os.path.join(args.out, "mcp_runs.jsonl"), {
        "started_utc": sc.iso_utc(time.time()), "git_head": "fake", "checks_dir": args.checks_dir,
        "service_capabilities": {"sha256": os.environ.get("FAKE_BENCH_MCP_CAPS", "caps-1")}})
    bad = 0
    for sid in [i for i in args.ids.split(",") if i]:
        if (sid, "mcp", args.repeat) in done:
            continue
        s = scen.get(sid) or {}
        t0 = time.time()
        time.sleep(wall / 1000.0)
        t1 = time.time()
        raw = sc.raw_name(sid, "mcp", args.repeat)
        line = {"run_id": "fake", "scenario": sid, "area": s.get("area"), "arm": "mcp", "repeat": args.repeat,
                "cache": "cold" if args.repeat == 0 else "warm", "started_utc": sc.iso_utc(t0),
                "ended_utc": sc.iso_utc(t1), "wall_ms": (t1 - t0) * 1000.0, "http": 200, "outcome": "pass",
                "checks": [{"check": "status", "outcome": "pass", "detail": "fake"}], "server": {"elapsed_ms": 1.0},
                "answer_digest": "fake", "raw": raw, "error": None, "time_limited": "",
                "mcp_phases": {"submit_ms": 1.0, "queue_ms": 2.0, "backend_ms": 3.0, "terminal_lag_ms": 4.0,
                               "paging_ms": 0.0},
                "job": {"search_id": "s-%s-%d" % (sid, args.repeat), "job_id": None}}
        if sid in errors:
            line.update(outcome="error", http=0, checks=[], raw=None,
                        error={"code": "runner_error", "text": "fake: not terminal within the deadline"})
            bad += 1
        else:
            sc.write_gz_json(os.path.join(args.out, raw), {"scenario": sid, "arm": "mcp", "repeat": args.repeat,
                                                           "answer": {"http": 200, "body": fake_answer(
                                                               s.get("route"), s.get("request")), "error": None}})
        sc.append_jsonl(path, line)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
