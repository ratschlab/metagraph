#!/usr/bin/env python3
"""A stand-in for the async session's MCP runner (scripts/bench/suite_mcp_runner.py) with its per-id CLI, for
the orchestrator's offline tests: --ids A,B --repeat N --out DIR --manifest M. It appends one execution line
per id (cache cold for repeat 0, warm otherwise; outcome not_implemented_mcp when the manifest marks the
scenario so) to DIR/executions.jsonl and its job to DIR/jobs.jsonl, skips what is recorded (resumable), and
logs each call to DIR/calls.log with whether the environment variable SUITE_TEST_SECRET reached it (the value
itself is never written). No request is made."""

import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import suite_common as sc  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--ids", required=True)
    ap.add_argument("--repeat", type=int, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--manifest", required=True)
    args = ap.parse_args(argv)
    os.makedirs(args.out, exist_ok=True)
    scen = dict((s.get("id"), s) for s in sc.read_json(args.manifest).get("scenarios", []))
    path = os.path.join(args.out, "executions.jsonl")
    done = set(sc.execution_key(r) for r in sc.read_jsonl(path)[0])
    with open(os.path.join(args.out, "calls.log"), "a") as fh:
        fh.write("%s %d env=%s\n" % (args.ids, args.repeat, "yes" if os.environ.get("SUITE_TEST_SECRET") else "no"))
    for sid in [i for i in args.ids.split(",") if i]:
        if (sid, "mcp", args.repeat) in done:
            continue
        nim = ((scen.get(sid) or {}).get("mcp") or {}).get("not_implemented")
        t0 = time.time()
        time.sleep(0.02)
        t1 = time.time()
        line = {"run_id": "fake", "scenario": sid, "area": (scen.get(sid) or {}).get("area"), "arm": "mcp",
                "repeat": args.repeat, "cache": "cold" if args.repeat == 0 else "warm",
                "started_utc": sc.iso_utc(t0), "ended_utc": sc.iso_utc(t1), "wall_ms": (t1 - t0) * 1000.0,
                "http": None if nim else 200, "outcome": "not_implemented_mcp" if nim else "pass",
                "checks": [] if nim else [{"check": "status", "outcome": "pass", "detail": "fake"}],
                "server": None, "answer_digest": None, "raw": None, "error": None,
                "mcp_phases": None if nim else {"submit_ms": 1.0, "queue_ms": 2.0, "backend_ms": 3.0,
                                                "terminal_lag_ms": 1.0, "paging_ms": 1.0},
                "job": None if nim else {"search_id": "s-%s-%d" % (sid, args.repeat), "job_id": None}}
        sc.append_jsonl(path, line)
        if line["job"]:
            sc.append_jsonl(os.path.join(args.out, "jobs.jsonl"), line["job"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
