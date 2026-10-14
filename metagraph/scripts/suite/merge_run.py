#!/usr/bin/env python3
"""Merges the arms of a run directory into the run's top-level files (DESIGN §6).

    merge_run.py RUN_DIR [--reevaluated]

Each arm writes RUN_DIR/arms/<arm>/ (run.json, executions.jsonl, samples.jsonl, raw/, and for the MCP arm
jobs.jsonl); the orchestrator's samplers on mex write RUN_DIR/samples/*.jsonl. This writes
RUN_DIR/executions.jsonl (raw paths rewritten relative to RUN_DIR; one line per scenario x arm x repeat, the
last recorded one when a resume repeated it; each line with `host_cache`, `since_prev_s` and
`host_cache_after` from the order in which the run's executions reached the host, see
suite_common.assign_host_cache), samples.jsonl (sorted by time, deduplicated, each sample's `phase` the arm or
loop it was taken during), jobs.jsonl, manifest.json and run.json (the identity of the run from the direct
arms, every arm's own run.json under `arms`). Re-running it rewrites the same files.

Outcomes (lane S3): every execution judged by lane S1's checks (not the stub) gets the outcome of
checks.scenario_outcome over its checks when that is stricter than the one its runner recorded (fail > finding >
inconclusive > pass; a runner's error, not_implemented_mcp or service_by_design is kept, and so is an MCP execution
the host-state rule made inconclusive: suite_common.MCP_HOST_STATE_CODES), the runner's own under
`outcome_runner`: the MCP runner of the async session applies its own rule, which does not tell property checks
from the generic ones (an execution whose property checks are all not exposed is inconclusive, not pass).

--reevaluated (lane S4): merges the re-evaluated set of reevaluate.py instead, each arm's
executions.reevaluated.jsonl into RUN_DIR/executions.reevaluated.jsonl and run.reevaluated.json (with the set's
mark, RUN_DIR/reevaluated.json, under `reevaluation`); scenario_outcome is the re-evaluation's checks'; an arm
without a re-evaluated file keeps its recorded lines (a warning says so); samples.jsonl, jobs.jsonl and
manifest.json, which no re-evaluation changes, are written only where missing, so every recorded file stays.
"""

import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import suite_common as sc  # noqa: E402

ARM_ORDER = {"direct-local": 0, "direct-remote": 1, "mcp": 2}
# the outcomes a check-derived outcome may replace, by strictness
_STRICTNESS = {"pass": 0, "inconclusive": 1, "finding": 2, "fail": 3}


def load_scenario_outcome(checks_dir=None):
    """checks.scenario_outcome of checks.py in |checks_dir| (default: next to this script), or None."""
    path = os.path.join(checks_dir or os.path.dirname(os.path.abspath(__file__)), "checks.py")
    if not os.path.exists(path):
        return None
    import importlib.util
    spec = importlib.util.spec_from_file_location("suite_checks_for_outcomes", path)
    mod = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(mod)
    except Exception:  # a checks.py that does not load leaves the runners' outcomes
        return None
    return getattr(mod, "scenario_outcome", None)


def judged_by_stub(run):
    """Whether an arm's run.json names the checks stub (run_direct: "stub"; the MCP runner: {file, stub})."""
    mod = (run or {}).get("checks_module")
    return mod == "stub" or (isinstance(mod, dict) and bool(mod.get("stub")))


def host_state_inconclusive(row):
    """An MCP execution the host-state rule made inconclusive (an answer refused with one of
    suite_common.MCP_HOST_STATE_CODES): its failing checks judge an environmental state, not a result."""
    err = row.get("error") if isinstance(row.get("error"), dict) else {}
    return row.get("arm") == "mcp" and row.get("outcome") == "inconclusive" and \
        str(err.get("code")) in sc.MCP_HOST_STATE_CODES


def stricter_outcomes(rows, arm_runs, outcome_fn):
    """Applies checks.scenario_outcome where it is stricter than the runner's outcome (module docstring); returns
    {arm: executions changed}. A re-evaluated line (`reevaluation`) was judged by checks.py whatever its runner
    ran."""
    changed = {}
    if outcome_fn is None:
        return changed
    for row in rows:
        old = row.get("outcome")
        checks = row.get("checks")
        if old not in _STRICTNESS or not checks or (
                judged_by_stub(arm_runs.get(row.get("arm"))) and "reevaluation" not in row):
            continue
        if host_state_inconclusive(row):
            continue   # the runner's host-state rule (or the re-evaluation's) decided it, not the checks
        try:
            new = outcome_fn(checks)
        except Exception:
            continue
        if _STRICTNESS.get(new, 3) > _STRICTNESS[old]:
            row.setdefault("outcome_runner", old)
            row["outcome"] = new
            changed[row.get("arm")] = changed.get(row.get("arm"), 0) + 1
    return changed


def arm_dirs(run_dir):
    root = os.path.join(run_dir, "arms")
    if not os.path.isdir(root):
        return []
    names = [n for n in os.listdir(root) if os.path.isdir(os.path.join(root, n)) and not n.startswith(".")]
    names.sort(key=lambda n: (ARM_ORDER.get(n, 9), n))
    return [(n, os.path.join(root, n)) for n in names]


def clock_offsets(orchestrator):
    """{arm: seconds its clock is ahead of the laptop's}: direct-local is stamped on mex, whose offset the
    orchestrator measures per invocation (the median of the invocations is used)."""
    vals = []
    for sess in (orchestrator or {}).get("sessions", []) or []:
        v = sess.get("mex_clock_offset_s")
        try:
            vals.append(float(v))
        except (TypeError, ValueError):
            continue
    if not vals:
        return {}
    vals.sort()
    return {"direct-local": vals[len(vals) // 2]}


def sample_files(run_dir):
    """The orchestrator's sample files (RUN_DIR/samples/*.jsonl), as (phase, path)."""
    root = os.path.join(run_dir, "samples")
    if not os.path.isdir(root):
        return []
    out = []
    for name in sorted(os.listdir(root)):
        if name.endswith(".jsonl"):
            out.append((name.split("-", 1)[0], os.path.join(root, name)))
    return out


def write_jsonl(path, rows):
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        for row in rows:
            fh.write(json.dumps(row, sort_keys=True, separators=(",", ":"), ensure_ascii=False) + "\n")
    os.replace(tmp, path)


def reevaluation_mark(run_dir):
    """RUN_DIR/reevaluated.json (reevaluate.py's mark of the re-evaluated set), or None."""
    return sc.read_json(os.path.join(run_dir, "reevaluated.json"))


def merge(run_dir, reevaluated=False):
    """Writes the merged files (of the re-evaluated set with |reevaluated|); returns the merged run.json."""
    name_of = lambda n: sc.reevaluated_name(n, reevaluated)  # noqa: E731
    mark = reevaluation_mark(run_dir) if reevaluated else None
    if reevaluated and mark is None:
        raise SystemExit("merge_run: %s has no reevaluated.json: run reevaluate.py first" % run_dir)
    warnings = []
    executions = {}
    order = []
    samples = []
    jobs = []
    arm_runs = {}
    arm_counts = {}
    manifest_copy = None
    for name, path in arm_dirs(run_dir):
        run = sc.read_json(os.path.join(path, "run.json"), {}) or {}
        arm_runs[name] = run
        exec_name = name_of("executions.jsonl")
        if reevaluated and not os.path.exists(os.path.join(path, exec_name)):
            warnings.append("arms/%s: not re-evaluated, its recorded lines are merged" % name)
            exec_name = "executions.jsonl"
        rows, bad = sc.read_jsonl(os.path.join(path, exec_name))
        if bad:
            warnings.append("arms/%s/%s: %d unparsable line(s) ignored" % (name, exec_name, bad))
        dropped = 0
        for row in rows:
            arm = row.get("arm")
            if name in ARM_ORDER and arm != name:
                dropped += 1
                continue
            row = dict(row)
            raw = row.get("raw")
            if isinstance(raw, str) and raw and not os.path.isabs(raw) and not raw.startswith("arms/"):
                row["raw"] = "arms/%s/%s" % (name, raw)
            key = sc.execution_key(row)
            if key not in executions:
                order.append(key)
            executions[key] = row
        arm_counts[name] = len(rows) - dropped
        if dropped:
            warnings.append("arms/%s: %d line(s) of another arm dropped" % (name, dropped))
        srows, sbad = sc.read_jsonl(os.path.join(path, "samples.jsonl"))
        if sbad:
            warnings.append("arms/%s/samples.jsonl: %d unparsable line(s) ignored" % (name, sbad))
        for s in srows:
            s = dict(s)
            s.setdefault("phase", name)
            samples.append(s)
        jrows, _ = sc.read_jsonl(os.path.join(path, "jobs.jsonl"))
        jobs.extend(jrows)
        if manifest_copy is None and os.path.exists(os.path.join(path, "manifest.json")):
            manifest_copy = os.path.join(path, "manifest.json")

    for phase, path in sample_files(run_dir):
        srows, sbad = sc.read_jsonl(path)
        if sbad:
            warnings.append("samples/%s: %d unparsable line(s) ignored" % (os.path.basename(path), sbad))
        for s in srows:
            s = dict(s)
            s.setdefault("phase", phase)
            samples.append(s)

    orchestrator = sc.read_json(os.path.join(run_dir, "orchestrator.json"))
    offsets = clock_offsets(orchestrator)
    merged_rows = sc.assign_host_cache([executions[k] for k in order], offsets)
    checks_dir = None
    if mark is not None:
        checks_dir = (mark.get("checks") or {}).get("dir")
        if not checks_dir or not os.path.exists(os.path.join(checks_dir, "checks.py")):
            warnings.append("the re-evaluation's checks (%s) are gone: scenario_outcome of the checks here" % checks_dir)
            checks_dir = None
    rederived = stricter_outcomes(merged_rows, arm_runs, load_scenario_outcome(checks_dir))
    for arm, n in sorted(rederived.items()):
        warnings.append("arms/%s: %d outcome(s) made stricter by checks.scenario_outcome (outcome_runner keeps the "
                        "runner's)" % (arm, n))
    seen = set()
    unique_samples = []
    for s in sorted(samples, key=lambda s: (sc.parse_utc(s.get("t_utc")) or 0.0)):
        key = (s.get("t_utc"), s.get("container"))
        if key in seen:
            continue
        seen.add(key)
        unique_samples.append(s)

    write_jsonl(os.path.join(run_dir, name_of("executions.jsonl")), merged_rows)
    # a re-evaluation changes neither the samples nor the jobs: the recorded merge's files stay
    for fname, data in (("samples.jsonl", unique_samples), ("jobs.jsonl", jobs)):
        if not reevaluated or not os.path.exists(os.path.join(run_dir, fname)):
            write_jsonl(os.path.join(run_dir, fname), data)
    top_manifest = os.path.join(run_dir, "manifest.json")
    if manifest_copy and not os.path.exists(top_manifest):
        sc.write_json_atomic(top_manifest, sc.read_json(manifest_copy))

    direct = [arm_runs[a] for a in ("direct-local", "direct-remote") if a in arm_runs]
    reference = direct[0] if direct else (next(iter(arm_runs.values())) if arm_runs else {})
    shas = set()
    for name, run in arm_runs.items():
        sha = (run.get("manifest") or {}).get("computed_sha256") or (run.get("manifest") or {}).get("content_sha256")
        if sha:
            shas.add(sha)
    if len(shas) > 1:
        warnings.append("the arms ran different manifests: %s" % ", ".join(sorted(shas)))
    instances = set()
    for run in direct:
        for sess in run.get("sessions", []):
            for key in ("server_instance_start", "server_instance_end"):
                if sess.get(key):
                    instances.add(sess[key])
    if len(instances) > 1:
        warnings.append("the server restarted during the run (server_instance %s)" % ", ".join(sorted(instances)))
    for name, run in arm_runs.items():
        for w in run.get("warnings", []) or []:
            warnings.append("arms/%s: %s" % (name, w))

    arms_summary = {}
    for name, run in arm_runs.items():
        slim = dict((k, v) for k, v in run.items() if k != "capabilities")
        slim["executions_recorded"] = arm_counts.get(name, 0)
        arms_summary[name] = slim
    merged = {
        "run_id": reference.get("run_id") or os.path.basename(os.path.abspath(run_dir)),
        "merged_utc": sc.iso_utc(time.time()),
        "manifest": reference.get("manifest"),
        "release": reference.get("release"), "server_instance": reference.get("server_instance"),
        "index_fp": reference.get("index_fp"), "capabilities": reference.get("capabilities"),
        "scripts_commit": reference.get("scripts_commit"),
        "host_load_start": (arm_runs.get("direct-local") or {}).get("host_load_start"),
        "host_load_end": (arm_runs.get("direct-local") or {}).get("host_load_end"),
        "arms": arms_summary,
        "orchestrator": orchestrator,
        "host_cache": {"warm_window_s": sc.HOST_CACHE_WARM_S, "clock_offsets_s": offsets},
        "counts": {"executions": len(merged_rows), "samples": len(unique_samples), "jobs": len(jobs)},
        "outcomes_rederived": rederived,
        "warnings": warnings,
    }
    if mark is not None:
        merged["reevaluation"] = dict((k, v) for k, v in mark.items() if k != "changes")
        merged["reevaluation"]["changes"] = len(mark.get("changes") or [])
    sc.write_json_atomic(os.path.join(run_dir, name_of("run.json")), merged)
    return merged


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    reevaluated = "--reevaluated" in argv
    argv = [a for a in argv if a != "--reevaluated"]
    if len(argv) != 1 or argv[0].startswith("-"):
        sys.stderr.write(__doc__)
        return 2
    merged = merge(argv[0], reevaluated)
    sys.stdout.write("merged %s: %d executions, %d samples, arms %s\n" % (
        argv[0], merged["counts"]["executions"], merged["counts"]["samples"], ", ".join(sorted(merged["arms"]))))
    for w in merged["warnings"]:
        sys.stdout.write("warning: %s\n" % w)
    return 0


if __name__ == "__main__":
    sys.exit(main())
