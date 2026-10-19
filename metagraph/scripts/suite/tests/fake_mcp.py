#!/usr/bin/env python3
"""A stand-in for the async session's MCP runner (scripts/bench/suite_mcp_runner.py) with its per-id CLI, for
the orchestrator's offline tests: --ids A,B --repeat N --out DIR --manifest M [--checks-dir DIR]. It appends one
execution line per id (cache cold for repeat 0, warm otherwise; outcome not_implemented_mcp when the manifest marks the
scenario so) to DIR/executions.jsonl and its job to DIR/jobs.jsonl, skips what is recorded (resumable), and
logs each call to DIR/calls.log with whether the environment variable SUITE_TEST_SECRET reached it (the value
itself is never written). No request is made. The resume rule as async's runner applies it to the ids of one call: a
listed warm-up whose execution of --repeat is recorded while its listed target's is not runs again right before the
target (a line with rerun_for_resume true and rerun_target, the earlier one untouched).

FAKE_MCP_MODE (as async's runner ends): `fail` every implemented execution fails its checks (exit 1);
`error_cancel` every one ends `error` at the poll deadline with its job cancelled (jobs.jsonl `created`, then
`canceled`; no backend phase; exit 1); `exit3`, `exit4`, `exit5` write nothing and exit 3 (a start refusal), 4 (the
by-design table refused) or 5 (the run's snapshot refused: a checks file not SNAPSHOT.json's); `crash` writes a traceback
to stderr and no line for any id (exit 1, a runner defect), only for a call that names FAKE_MCP_CRASH_ID when that
is set; `unknown` records the outcome `skipped_by_runner`, outside the shared vocabulary (exit 0); `sigint` waits
up to 30 s on each call (only on a call that names FAKE_MCP_SIGINT_ID when that is set) and, interrupted (Ctrl-C),
prints a last line on stdout as async's runner does and notes in calls.log whether that print worked (exit 130) or met
a closed pipe (exit 1).

Like async's runner it keeps the manifest's copy in DIR/manifest.json and writes DIR/run.json with the manifest's
identity, `service_capabilities` (the live read of the service's configuration; FAKE_MCP_CAPS `missing` leaves it
out, `failed` records a failed read, `mismatch` a live / table mismatch; its sha256 service_document_sha256's) and the
checks it judges with: `checks_dir` (--checks-dir, the run's snapshot under run_suite.sh's {checks_dir}),
`checks_sha256` and `suite_common_sha256` of its files, and `checks_module` {file, stub}. FAKE_MCP_SERVICE_COMMIT set:
the deployed service's identity as async's runner 684bb327 records it, on run.json's top level and a line of
DIR/mcp_runs.jsonl per invocation: service_commit, identity_changed (false; FAKE_MCP_IDENTITY_CHANGED=1 true) and
identity_error null."""

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
    if mode in ("exit3", "exit4", "exit5"):
        return int(mode[-1])
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
        # what the runner judges with, as async's runner is to record it per invocation (one invocation here)
        run.update(checks_dir=args.checks_dir,
                   checks_sha256=sc.file_sha256(os.path.join(args.checks_dir, "checks.py")),
                   suite_common_sha256=sc.file_sha256(os.path.join(args.checks_dir, "suite_common.py")),
                   checks_module={"file": os.path.join(args.checks_dir, "checks.py"), "stub": False})
    caps = os.environ.get("FAKE_MCP_CAPS", "")
    if caps != "missing":
        doc = {"service_settings": {"advanced_enabled": False}}
        run["service_capabilities"] = {"document": None if caps == "failed" else doc,
                                       "read_at": sc.iso_utc(time.time()), "base": "fake",
                                       "mismatches": [{"path": "service_settings.advanced_enabled", "live": True,
                                                       "table": False}] if caps == "mismatch" else []}
        if caps == "failed":
            run["service_capabilities"]["error"] = "fake: connection refused"
        # the record's identity as async's runner states it (one definition: suite_common.service_document_sha256)
        run["service_capabilities"]["sha256"] = sc.service_document_sha256(run["service_capabilities"])
    if os.environ.get("FAKE_MCP_SERVICE_COMMIT"):
        run.update(service_commit=os.environ["FAKE_MCP_SERVICE_COMMIT"],
                   identity_changed=os.environ.get("FAKE_MCP_IDENTITY_CHANGED") == "1", identity_error=None)
        sc.append_jsonl(os.path.join(args.out, "mcp_runs.jsonl"), dict(run, ids=args.ids, repeat=args.repeat))
    sc.write_json_atomic(os.path.join(args.out, "run.json"), run)
    path = os.path.join(args.out, "executions.jsonl")
    done = set(sc.execution_key(r) for r in sc.read_jsonl(path)[0])
    with open(os.path.join(args.out, "calls.log"), "a") as fh:
        fh.write("%s %d env=%s\n" % (args.ids, args.repeat, "yes" if os.environ.get("SUITE_TEST_SECRET") else "no"))
    if mode == "sigint" and os.environ.get("FAKE_MCP_SIGINT_ID", args.ids) in args.ids.split(","):
        try:
            time.sleep(30)
        except KeyboardInterrupt:
            try:
                print("done: interrupted on %s" % args.ids, flush=True)
                note, rc = "final print ok", 130
            except BrokenPipeError:
                note, rc = "final print raised BrokenPipeError", 1
            with open(os.path.join(args.out, "calls.log"), "a") as fh:
                fh.write(note + "\n")
            return rc
    bad = 0
    listed = [i for i in args.ids.split(",") if i]
    steps = []
    for sid in listed:
        target = (scen.get(sid) or {}).get("warms")
        if (sid, "mcp", args.repeat) in done:
            if (target in listed and (target, "mcp", args.repeat) not in done
                    and not ((scen.get(sid) or {}).get("mcp") or {}).get("not_implemented")):
                steps.append((sid, target))
            continue
        steps.append((sid, None))
    for sid, rerun_target in steps:
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
        if rerun_target:
            line.update(rerun_for_resume=True, rerun_target=rerun_target)
            line["job"]["search_id"] += "-rerun"
        if nim:
            pass
        elif mode == "fail":
            line.update(outcome="fail", checks=[{"check": "status", "outcome": "fail", "detail": "fake: http 500"}],
                        http=500)
        elif mode == "unknown":
            line.update(outcome="skipped_by_runner")
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
