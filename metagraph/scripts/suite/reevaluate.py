#!/usr/bin/env python3
"""Re-evaluates a recorded run OFFLINE with the current checks (no request to any host).

    reevaluate.py RUN_DIR [--checks DIR] [--quiet]

For every arm under RUN_DIR/arms/ and every execution line of its executions.jsonl (the last one per scenario,
arm and repeat, as merge_run keeps it), the raw answer the runner recorded (raw/...json.gz: both runners keep the
NORMALISED answer under `answer`, the MCP runner the shape of its normaliser, so no normaliser runs again) is
judged again with checks.evaluate of DIR/checks.py (default: this directory), with what the runner passed:
`previous` = the earlier repeats of the scenario on the arm, in repeat order, as {digest, time_limited} computed
from their raw answers with the same checks (the recorded digest and time_limited where a raw answer is missing),
and `client` = {"wall_ms"}: a direct arm's recorded wall time (run_direct.py), the MCP arm's service-measured
backend time (mcp_phases.backend_ms, as async's runner passes it since its 0d127f4a). The outcome is
checks.scenario_outcome of the new checks, except where the runner decided it without them, which is kept: `error`
(no answer, a transport or decode failure, a busy deployment; unless the recorded checks are the checks module's
own exception, which a re-evaluation repeats), `not_implemented_mcp`, `service_by_design` (async's table), and the
MCP runner's host-state rule (an unexpected `uninterruptible_overrun_exceeds_allowance` is inconclusive).

Writes, keeping every original file:
  RUN_DIR/arms/<arm>/executions.reevaluated.jsonl   the re-evaluated lines (the recorded fields, with checks,
      outcome, answer_digest and time_limited replaced, and `reevaluation`: {checks_sha256, outcome_recorded,
      checks_changed: the checks whose outcome changed, added or dropped})
  RUN_DIR/reevaluated.json   the re-evaluated set's mark: the checks (directory, commit, `+dirty`, sha256 of
      checks.py and suite_common.py), the recorded checks (each arm's run.json), counts, and every execution whose
      outcome or checks changed
Then `merge_run.py --reevaluated`, `compare.py --reevaluated`, `report.py --reevaluated` (or `run_suite.sh
--reevaluate RUN_DIR`) build the re-evaluated set's executions.reevaluated.jsonl, run.reevaluated.json,
compare.reevaluated.json and report.reevaluated.md / .pdf beside the recorded ones.
"""

import argparse
import hashlib
import importlib.util
import inspect
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import merge_run  # noqa: E402
import suite_common as sc  # noqa: E402

MARK = "reevaluated.json"
# outcomes a runner decides without the checks (kept by a re-evaluation)
KEPT_OUTCOMES = ("error", "not_implemented_mcp", "service_by_design")
# async's MCP runner (suite_mcp_runner.py HOST_STATE_CODES): an environmental host state the scenario does not
# expect is inconclusive, whatever the checks say
MCP_HOST_STATE_CODES = ("uninterruptible_overrun_exceeds_allowance",)


def sha256_file(path):
    if not os.path.exists(path):
        return None
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def git_commit(directory):
    """(commit, dirty) of the tree holding |directory|: dirty when checks.py or suite_common.py differ from it;
    (None, None) outside a git tree."""
    try:
        head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=directory, stdout=subprocess.PIPE,
                              stderr=subprocess.DEVNULL, universal_newlines=True, timeout=30)
        if head.returncode != 0:
            return None, None
        st = subprocess.run(["git", "status", "--porcelain", "--", "checks.py", "suite_common.py"], cwd=directory,
                            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, universal_newlines=True, timeout=30)
        return head.stdout.strip(), bool(st.stdout.strip())
    except (OSError, subprocess.SubprocessError):
        return None, None


def load_checks(directory):
    """checks.py of |directory| as a module of its own, whose `import suite_common` resolves to that directory's
    copy (the digest and normaliser of the same checks), this script's own suite_common left in place."""
    path = os.path.join(directory, "checks.py")
    if not os.path.exists(path):
        raise SystemExit("reevaluate: no checks.py in %s" % directory)
    own = sys.modules.get("suite_common")
    saved_path = list(sys.path)
    sys.modules.pop("suite_common", None)
    sys.path.insert(0, directory)
    try:
        spec = importlib.util.spec_from_file_location("suite_checks_reevaluation", path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
    finally:
        sys.path[:] = saved_path
        if own is not None:
            sys.modules["suite_common"] = own
        else:
            sys.modules.pop("suite_common", None)
    return mod


def checks_identity(directory):
    commit, dirty = git_commit(directory)
    return {"dir": os.path.abspath(directory), "commit": commit, "dirty": dirty,
            "commit_label": (commit + ("+dirty" if dirty else "")) if commit else None,
            "checks_sha256": sha256_file(os.path.join(directory, "checks.py")),
            "suite_common_sha256": sha256_file(os.path.join(directory, "suite_common.py"))}


def arm_lines(path):
    """The arm's execution lines, the last one per scenario, arm and repeat, in first-seen order; (lines, bad)."""
    rows, bad = sc.read_jsonl(path)
    keyed, order = {}, []
    for row in rows:
        key = sc.execution_key(row)
        if key not in keyed:
            order.append(key)
        keyed[key] = row
    return [keyed[k] for k in order], bad


def raw_answer(arm_dir, line):
    """The recorded normalised answer of an execution line, or None (no raw file, or none kept in it)."""
    raw = line.get("raw")
    if not isinstance(raw, str) or not raw:
        return None
    path = raw if os.path.isabs(raw) else os.path.join(arm_dir, raw)
    if not os.path.exists(path) and raw.startswith("arms/"):
        path = os.path.join(os.path.dirname(os.path.dirname(arm_dir)), raw)
    try:
        doc = sc.read_gz_json(path)
    except (OSError, ValueError, EOFError):
        return None
    answer = (doc or {}).get("answer") if isinstance(doc, dict) else None
    return answer if isinstance(answer, dict) else None


def client_of(arm, line):
    """What the arm's runner passed as checks.evaluate's `client` (module docstring), or None."""
    if arm == "mcp":
        v = (line.get("mcp_phases") or {}).get("backend_ms")
    else:
        v = line.get("wall_ms")
    return {"wall_ms": float(v)} if isinstance(v, (int, float)) and not isinstance(v, bool) else None


def evaluate_raised(checks):
    return bool(checks) and all(c.get("check") == "evaluate" for c in checks)


def outcome_of(checks_mod, arm, scenario, answer, line, new_checks, raised):
    """The re-evaluated outcome (module docstring) and why it was kept, if it was."""
    recorded = line.get("outcome")
    if recorded in KEPT_OUTCOMES and not (recorded == "error" and evaluate_raised(line.get("checks"))):
        return recorded, "decided by the runner without the checks"
    if raised:
        return ("fail" if arm == "mcp" else "error"), None
    err = (answer or {}).get("error") or {}
    if arm == "mcp" and err and str(err.get("code")) in MCP_HOST_STATE_CODES and \
            str(err.get("code")) not in json.dumps(scenario.get("expect") or []):
        return "inconclusive", "the MCP runner's host-state rule (%s)" % err.get("code")
    rule = getattr(checks_mod, "scenario_outcome", None) or sc.outcome_from_checks
    return rule(new_checks), None


def changed_checks(old, new):
    """[{check, before, after, detail}] for every check whose outcome differs, appears or disappears (by name and
    position among the checks of that name)."""
    def keyed(rows):
        out, seen = {}, {}
        for r in rows or []:
            name = r.get("check")
            n = seen.get(name, 0)
            seen[name] = n + 1
            out[(name, n)] = r
        return out
    a, b = keyed(old), keyed(new)
    out = []
    for key in list(a) + [k for k in b if k not in a]:
        before, after = (a.get(key) or {}).get("outcome"), (b.get(key) or {}).get("outcome")
        if before != after:
            out.append({"check": key[0], "before": before, "after": after,
                        "detail": ((b.get(key) or a.get(key) or {}).get("detail") or "")[:400]})
    return out


def reevaluate(run_dir, checks_dir=None, log=None):
    """Writes the re-evaluated set (module docstring); returns the mark (reevaluated.json's content)."""
    checks_dir = os.path.abspath(checks_dir or HERE)
    checks_mod = load_checks(checks_dir)
    ident = checks_identity(checks_dir)
    takes = inspect.signature(checks_mod.evaluate).parameters
    takes_prev, takes_client = "previous" in takes, "client" in takes
    short = (ident["checks_sha256"] or "")[:16]
    top_manifest = sc.read_json(os.path.join(run_dir, "manifest.json"))
    mark = {"run_dir": os.path.abspath(run_dir), "reevaluated_utc": sc.iso_utc(time.time()), "checks": ident,
            "client_passed": {"direct-local": "wall_ms", "direct-remote": "wall_ms", "mcp": "mcp_phases.backend_ms"},
            "recorded": {}, "arms": {}, "changes": [], "warnings": []}
    for name, arm_dir in merge_run.arm_dirs(run_dir):
        run = sc.read_json(os.path.join(arm_dir, "run.json"), {}) or {}
        mark["recorded"][name] = {"checks_module": run.get("checks_module"),
                                  "scripts_commit": run.get("scripts_commit") or run.get("git_head")}
        manifest = sc.read_json(os.path.join(arm_dir, "manifest.json")) or top_manifest
        scen = dict((s.get("id"), s) for s in (manifest or {}).get("scenarios", [])) if isinstance(manifest, dict) else {}
        lines, bad = arm_lines(os.path.join(arm_dir, "executions.jsonl"))
        if bad:
            mark["warnings"].append("arms/%s/executions.jsonl: %d unparsable line(s) ignored" % (name, bad))
        counts = {"lines": len(lines), "reevaluated": 0, "kept": 0, "no_answer": 0, "no_scenario": 0,
                  "outcome_changed": 0, "checks_changed": 0}
        # each line's answer, and its {digest, time_limited} by this checks module (the recorded ones where the raw
        # answer is missing), for the later repeats' `previous`
        answers, entries, repeats = {}, {}, {}
        for line in lines:
            key = sc.execution_key(line)
            answers[key] = answer = raw_answer(arm_dir, line)
            route = (scen.get(line.get("scenario")) or {}).get("route") or line.get("route")
            try:
                if answer is not None:
                    entries[key] = {"digest": checks_mod.answer_digest(answer, route),
                                    "time_limited": checks_mod.time_limited(answer)
                                    if hasattr(checks_mod, "time_limited") else ""}
            except Exception as exc:  # the recorded values, never fatal to the re-evaluation
                mark["warnings"].append("arms/%s %s r%s: digest or time limitation raised %s: %s" % (
                    name, line.get("scenario"), line.get("repeat"), type(exc).__name__, exc))
            if key not in entries and line.get("answer_digest"):
                entries[key] = {"digest": line.get("answer_digest"), "time_limited": line.get("time_limited") or ""}
            repeats.setdefault(line.get("scenario"), []).append(line)
        out_rows = []
        for line in lines:
            sid, rep = line.get("scenario"), line.get("repeat") or 0
            scenario = scen.get(sid)
            key = sc.execution_key(line)
            answer = answers.get(key)
            row = dict(line)
            if scenario is None or answer is None:
                counts["no_scenario" if scenario is None else "no_answer"] += 1
                row["reevaluation"] = {"checks_sha256": short, "outcome_recorded": line.get("outcome"),
                                       "checks_changed": [],
                                       "kept": "no scenario in the manifest" if scenario is None else "no recorded answer"}
                out_rows.append(row)
                continue
            route = scenario.get("route")
            previous = [entries[sc.execution_key(o)] for o in sorted(repeats.get(sid, []), key=lambda o: o.get("repeat") or 0)
                        if (o.get("repeat") or 0) < rep and sc.execution_key(o) in entries]
            kwargs = {}
            if takes_prev:
                kwargs["previous"] = previous
            client = client_of(name, line)
            if takes_client and client is not None:
                kwargs["client"] = client
            raised = False
            try:
                new_checks = list(checks_mod.evaluate(scenario, answer, name, **kwargs))
            except Exception as exc:  # recorded as the runners record it, never fatal to the re-evaluation
                new_checks = [{"check": "evaluate", "outcome": "fail",
                               "detail": "exception in checks.evaluate: %s: %s" % (type(exc).__name__, exc)}]
                raised = True
            outcome, kept = outcome_of(checks_mod, name, scenario, answer, line, new_checks, raised)
            diff = changed_checks(line.get("checks"), new_checks)
            entry = entries.get(key) or {"digest": line.get("answer_digest"), "time_limited": line.get("time_limited")}
            row.update({"checks": new_checks, "outcome": outcome, "answer_digest": entry["digest"],
                        "time_limited": entry["time_limited"]})
            if line.get("answer_digest") and line.get("answer_digest") != row["answer_digest"]:
                row["answer_digest_recorded"] = line.get("answer_digest")
            row["reevaluation"] = {"checks_sha256": short, "outcome_recorded": line.get("outcome"),
                                   "checks_changed": [c["check"] for c in diff]}
            if kept:
                row["reevaluation"]["kept"] = kept
            counts["kept" if kept else "reevaluated"] += 1
            if outcome != line.get("outcome"):
                counts["outcome_changed"] += 1
            if diff:
                counts["checks_changed"] += 1
            if diff or outcome != line.get("outcome"):
                mark["changes"].append({"scenario": sid, "arm": name, "repeat": rep, "area": line.get("area"),
                                        "before": line.get("outcome"), "after": outcome, "kept": kept,
                                        "checks": diff})
            out_rows.append(row)
        merge_run.write_jsonl(os.path.join(arm_dir, sc.reevaluated_name("executions.jsonl")), out_rows)
        mark["arms"][name] = counts
        if log:
            log("%s: %d line(s), %d re-evaluated, %d kept, %d without an answer; %d outcome(s) and %d check list(s) "
                "changed" % (name, counts["lines"], counts["reevaluated"], counts["kept"], counts["no_answer"],
                             counts["outcome_changed"], counts["checks_changed"]))
    sc.write_json_atomic(os.path.join(run_dir, MARK), mark)
    return mark


def main(argv=None):
    ap = argparse.ArgumentParser(description="Re-evaluate a recorded run offline with the current checks.")
    ap.add_argument("run_dir")
    ap.add_argument("--checks", default=None, help="the directory of the checks.py to judge with (default: here)")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args(argv)
    if not os.path.isdir(os.path.join(args.run_dir, "arms")):
        sys.stderr.write("reevaluate: %s has no arms/ directory\n" % args.run_dir)
        return 2
    log = None if args.quiet else (lambda m: sys.stdout.write(m + "\n"))
    mark = reevaluate(args.run_dir, args.checks, log)
    if not args.quiet:
        sys.stdout.write("checks %s (sha256 %s)\n" % (mark["checks"]["commit_label"] or mark["checks"]["dir"],
                                                       (mark["checks"]["checks_sha256"] or "")[:16]))
        for ch in mark["changes"]:
            sys.stdout.write("  %-22s %-13s r%d  %s -> %s%s\n" % (
                ch["scenario"], ch["arm"], ch["repeat"], ch["before"], ch["after"],
                "".join("; %s %s->%s" % (c["check"], c["before"], c["after"]) for c in ch["checks"])))
        for w in mark["warnings"]:
            sys.stdout.write("warning: %s\n" % w)
    return 0


if __name__ == "__main__":
    sys.exit(main())
