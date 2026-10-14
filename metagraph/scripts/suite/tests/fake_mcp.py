#!/usr/bin/env python3
"""A stand-in for the async session's MCP runner (scripts/bench/suite_mcp_runner.py) with its per-id CLI, for
the orchestrator's offline tests: --ids A,B --repeat N --out DIR --manifest M [--checks-dir DIR]. It appends one execution line
per id (cache cold for repeat 0, warm otherwise; outcome not_implemented_mcp when the manifest marks the
scenario so) to DIR/executions.jsonl and its job to DIR/jobs.jsonl, skips what is recorded (resumable), and
logs each call to DIR/calls.log with whether the environment variable SUITE_TEST_SECRET reached it (the value
itself is never written). No request is made.

FAKE_MCP_MODE (as async's runner ends): `fail` every implemented execution fails its checks (exit 1);
`error_cancel` every one ends `error` at the poll deadline with its job cancelled (jobs.jsonl `created`, then
`canceled`; no backend phase; exit 1); `exit3` writes nothing and exits 3 (a start refusal); `crash` writes a traceback
to stderr and no line for any id (exit 1, a runner defect), only for a call that names FAKE_MCP_CRASH_ID when that
is set.

Like async's runner it keeps the manifest's copy in DIR/manifest.json and writes DIR/run.json with the manifest's
identity, `service_capabilities` (the live read of the service's configuration; FAKE_MCP_CAPS `missing` leaves it
out, `failed` records a failed read, `mismatch` a live / table mismatch) and the checks it judges with: `checks_dir`
(--checks-dir, the run's snapshot under run_suite.sh's {checks_dir}) and `checks_module` {file, stub}."""

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
    ap.add_argument("--checks-dir", default=None)
    args = ap.parse_args(argv)
    mode = os.environ.get("FAKE_MCP_MODE", "")
    if mode == "exit3":
        return 3
    if mode == "crash" and os.environ.get("FAKE_MCP_CRASH_ID", args.ids) in args.ids.split(","):
        sys.stderr.write("Traceback (most recent call last):\n  File \"suite_mcp_runner.py\", line 1, in verdict\n"
                         "KeyError: 'fake crash'\n")
        return 1
    os.makedirs(args.out, exist_ok=True)
    manifest = sc.read_json(args.manifest)
    scen = dict((s.get("id"), s) for s in manifest.get("scenarios", []))
    if not os.path.exists(os.path.join(args.out, "manifest.json")):
        sc.write_json_atomic(os.path.join(args.out, "manifest.json"), manifest)
    ok, computed = sc.verify_manifest(manifest)
    run = {"arm": "mcp", "manifest": {"computed_sha256": computed, "content_sha256": manifest.get("content_sha256"),
                                      "sha_verified": ok}}
    if args.checks_dir:
        run.update(checks_dir=args.checks_dir,
                   checks_module={"file": os.path.join(args.checks_dir, "checks.py"), "stub": False})
    caps = os.environ.get("FAKE_MCP_CAPS", "")
    if caps != "missing":
        doc = {"service_settings": {"advanced_enabled": False}}
        run["service_capabilities"] = {"document": None if caps == "failed" else doc, "sha256": "fake",
                                       "read_at": sc.iso_utc(time.time()), "base": "fake",
                                       "mismatches": [{"path": "service_settings.advanced_enabled", "live": True,
                                                       "table": False}] if caps == "mismatch" else []}
        if caps == "failed":
            run["service_capabilities"]["error"] = "fake: connection refused"
    sc.write_json_atomic(os.path.join(args.out, "run.json"), run)
    path = os.path.join(args.out, "executions.jsonl")
    done = set(sc.execution_key(r) for r in sc.read_jsonl(path)[0])
    with open(os.path.join(args.out, "calls.log"), "a") as fh:
        fh.write("%s %d env=%s\n" % (args.ids, args.repeat, "yes" if os.environ.get("SUITE_TEST_SECRET") else "no"))
    bad = 0
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
        if nim:
            pass
        elif mode == "fail":
            line.update(outcome="fail", checks=[{"check": "status", "outcome": "fail", "detail": "fake: http 500"}],
                        http=500)
        elif mode == "error_cancel":
            line.update(outcome="error", checks=[], http=0, server=None,
                        error={"code": "runner_error", "text": "fake: not terminal within the deadline"},
                        mcp_phases={"submit_ms": 1.0, "queue_ms": None, "backend_ms": None,
                                    "terminal_lag_ms": None, "paging_ms": None})
            jid = line["job"]["search_id"]
            sc.append_jsonl(os.path.join(args.out, "jobs.jsonl"), {"event": "created", "kind": "search", "id": jid,
                                                                   "scenario": sid, "repeat": args.repeat,
                                                                   "at": sc.iso_utc(t0)})
            sc.append_jsonl(os.path.join(args.out, "jobs.jsonl"), {"event": "canceled", "kind": "search", "id": jid,
                                                                   "detail": "ok", "at": sc.iso_utc(time.time())})
        bad += line["outcome"] in ("fail", "error")
        sc.append_jsonl(path, line)
        if line["job"] and mode != "error_cancel":
            sc.append_jsonl(os.path.join(args.out, "jobs.jsonl"), line["job"])
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
