#!/usr/bin/env python3
"""Merges the arms of a run directory into the run's top-level files (DESIGN §6).

    merge_run.py RUN_DIR [--reevaluated] [--strict-identity]

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
the host-state rule made inconclusive where its plausibility rule admits it or cannot decide it (rule (c),
verdict_detail): suite_common.host_state_verdict on the raw record's service-side wall, and an MCP
execution a refused / relayed / answer-triggered by-design case decided, a finding included; the failing checks an
answered / clamped case explains, the line's explained_items (async's runner records them since f386a3c0, the
re-evaluation since lane S5c, B5b-06: the same indices), else every failing item of its explained_checks (a line
of an earlier runner), are left out: triage-2 B-01), the runner's own under
`outcome_runner`: the MCP runner of the async session applies its own rule, which does not tell property checks
from the generic ones (an execution whose property checks are all not exposed is inconclusive, not pass).

--reevaluated (lane S4): merges the re-evaluated set of reevaluate.py instead, each arm's
executions.reevaluated.jsonl into RUN_DIR/executions.reevaluated.jsonl and run.reevaluated.json (with the set's
mark, RUN_DIR/reevaluated.json, under `reevaluation`); scenario_outcome is the re-evaluation's checks'; an arm
without a re-evaluated file keeps its recorded lines (a warning says so); samples.jsonl, jobs.jsonl and
manifest.json, which no re-evaluation changes, are written only where missing, so every recorded file stays.

A re-normalised MCP arm (async's suite_mcp_renormalise.py rebuilt its answers from the recorded tool payloads after
the run): run.json `renormalisation` {arm: renormalisation()} states it (the summary's source, time, normaliser and
counts, the lines' marks, the raw records that keep the recorded answer as answer_recorded); the merge of the recorded
set warns that those lines' outcomes are stale until a re-evaluation, and the merge of a re-evaluated set copies it
into `reevaluation.renormalised`, which the report states.

The manifest identity (finding 9 of the outside review): run.json's `manifest_identity` holds the run manifest's
content_sha256 and relations_sha256 and each arm's (manifest_identity); an arm whose hashes differ, or a run manifest
whose stated hashes its content does not match, is a warning; with --strict-identity (run_suite.sh --gate) it is an
error: the files are written, then the exit status is 3.

The host's identity (release_identity): the releases and index_fps the direct arms' sessions saw at their starts and
closing GETs, in order (run.json `release_identity`, `mixed_releases`); more than one release or index_fp (a deploy
during the run, or between an interrupt and its resume) is a warning, and an error with --strict-identity.
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


def host_state_inconclusive(row, scenario=None, run_dir=None):
    """An MCP execution the host-state rule made inconclusive (an answer refused with one of
    suite_common.MCP_HOST_STATE_CODES the scenario does not expect, whose service-side wall, read from the raw record
    under |run_dir|, the plausibility rule of suite_common.host_state_verdict admits or cannot decide), or that rule
    (c) made inconclusive (its inconclusive_reason: a host-state code or plausibility undecidable; a line of an
    earlier runner without one is judged here): its failing checks judge an environmental state or an undecidable
    one, not a result."""
    if row.get("arm") != "mcp" or row.get("outcome") != "inconclusive":
        return False
    why = row.get("inconclusive_reason")   # the runner's (e334715b) or the re-evaluation's class of it
    if why is not None:
        return why in sc.UNDECIDABLE or str(why) in sc.MCP_HOST_STATE_CODES
    err = row.get("error") if isinstance(row.get("error"), dict) else {}
    if str(err.get("code")) not in sc.MCP_HOST_STATE_CODES:
        return False
    doc = None
    raw = row.get("raw")
    if run_dir and isinstance(raw, str) and raw:
        try:
            doc = sc.read_gz_json(raw if os.path.isabs(raw) else os.path.join(run_dir, raw))
        except (OSError, ValueError, EOFError):
            doc = None
    code, state, _ = sc.host_state_verdict(scenario or {}, row, doc if isinstance(doc, dict) else None)
    return code is not None and state in ("host_state", "undecidable")


def by_design_decided(row):
    """An MCP execution a refused / relayed / answer-triggered by-design case decided (service_by_design, or the
    finding a by-design mark never erases): its failing checks judge the service's refusal, not a result. An answered
    / clamped case's execution states the checks it explains (explained_checks) and is judged on the others."""
    return row.get("arm") == "mcp" and bool(row.get("service_by_design_case")) and \
        row.get("outcome") in ("service_by_design", "finding") and not row.get("explained_checks")


def stricter_outcomes(rows, arm_runs, outcome_fn, scenarios=None, run_dir=None):
    """Applies checks.scenario_outcome where it is stricter than the runner's outcome (module docstring); returns
    {arm: executions changed}. A re-evaluated line (`reevaluation`) was judged by checks.py whatever its runner
    ran. |scenarios| (id -> scenario) gives the host-state rule its budgets."""
    changed = {}
    if outcome_fn is None:
        return changed
    for row in rows:
        old = row.get("outcome")
        checks = row.get("checks")
        if old not in _STRICTNESS or not checks or (
                judged_by_stub(arm_runs.get(row.get("arm"))) and "reevaluation" not in row):
            continue
        if host_state_inconclusive(row, (scenarios or {}).get(row.get("scenario")), run_dir) or by_design_decided(row):
            continue   # the runner's host-state or by-design rule (or the re-evaluation's) decided it, not the checks
        items = row.get("explained_items")
        explained = set(row.get("explained_checks") or ())
        if isinstance(items, list) and items:
            # B5b-06: the explained check items (the runner's since f386a3c0, or the re-evaluation's), not every
            # item of an explained name
            checks = [c for i, c in enumerate(checks) if not (i in items and c.get("outcome") == "fail")]
        elif explained:
            checks = [c for c in checks if not (c.get("check") in explained and c.get("outcome") == "fail")]
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


def arm_identity(arm_dir, run):
    """The manifest identity an arm ran (suite_common.manifest_identity): computed from the manifest copy the arm
    kept (run_direct and async's runner keep one), else what its run.json states (a relations_sha256 of None: not
    stated); `source` says which."""
    copy = sc.read_json(os.path.join(arm_dir, "manifest.json"))
    if isinstance(copy, dict) and isinstance(copy.get("scenarios"), list):
        return dict(sc.manifest_identity(copy), source="manifest copy")
    info = (run or {}).get("manifest") or {}
    return {"content_sha256": info.get("computed_sha256") or info.get("content_sha256"),
            "relations_sha256": info.get("computed_relations_sha256") or info.get("relations_sha256"),
            "source": "run.json"}


def manifest_identity(run_dir, arm_runs):
    """The run's manifest identity against each arm's: {manifest: RUN_DIR/manifest.json's identity and the stated
    hashes it does not match (problems), arms: {arm: arm_identity}, differs: sentences naming each arm whose
    content_sha256 or relations_sha256 is not the run manifest's (an arm that does not state one is named
    too: what it ran cannot be shown to be the same)}. Exploratory runs read `differs` as a warning; gate.py as a
    failure."""
    top = sc.read_json(os.path.join(run_dir, "manifest.json"))
    if isinstance(top, dict):
        manifest = dict(sc.manifest_identity(top), problems=sc.manifest_problems(top))
    else:
        manifest = {"content_sha256": None, "relations_sha256": None, "problems": ["no manifest.json in the run"]}
    arms, differs = {}, []
    for name in sorted(arm_runs, key=lambda n: (ARM_ORDER.get(n, 9), n)):
        ident = arm_identity(os.path.join(run_dir, "arms", name), arm_runs[name])
        arms[name] = ident
        for key in ("content_sha256", "relations_sha256"):
            if ident.get(key) is None:
                if manifest.get(key) is not None:
                    differs.append("arms/%s states no %s (%s)" % (name, key, ident["source"]))
            elif manifest.get(key) is not None and ident[key] != manifest[key]:
                differs.append("arms/%s %s %s, the run's manifest %s (%s)" % (
                    name, key, ident[key][:12], manifest[key][:12], ident["source"]))
    if manifest.get("content_sha256") is None:
        # no top-level manifest: the arms against each other
        for key in ("content_sha256", "relations_sha256"):
            seen = sorted(set(i[key] for i in arms.values() if i.get(key)))
            if len(seen) > 1:
                differs.append("%s %s" % (key, ", ".join(x[:12] for x in seen)))
    return {"manifest": manifest, "arms": arms, "differs": differs}


def release_identity(arm_runs):
    """The host's identity over the direct arms' sessions (run_direct.py records release and index_fp at each session's
    start and its closing GET; an arm of before those fields states its first release only): {releases, index_fps: the
    distinct values in the order seen, changed: one sentence per value after the first, naming the arm, session and
    time that saw it}. Executions of a run that spans two deployments are attributed to one release that did not
    answer them all: merge_run warns, --strict-identity and the gate fail."""
    seen = []
    for arm in sc.DIRECT_ARMS:
        run = arm_runs.get(arm)
        if not isinstance(run, dict):
            continue
        t0 = sc.parse_utc(run.get("started_utc")) or 0.0
        for key in ("release", "index_fp"):
            if run.get(key):
                seen.append((t0, 0, arm, None, "first start", key, run[key]))
        for i, sess in enumerate(run.get("sessions") or []):
            for order, (where, tkey) in enumerate((("start", "started_utc"), ("end", "ended_utc"))):
                for key in ("release", "index_fp"):
                    v = sess.get("%s_%s" % (key, where))
                    if v:
                        seen.append((sc.parse_utc(sess.get(tkey)) or t0, 1 + order, arm, i, where, key, v))
    seen.sort(key=lambda x: (x[0], x[1], ARM_ORDER.get(x[2], 9)))
    out = {"releases": [], "index_fps": [], "changed": []}
    first = {}
    for t, _, arm, i, where, key, v in seen:
        values = out["releases" if key == "release" else "index_fps"]
        if v in values:
            continue
        values.append(v)
        at = "arms/%s %s" % (arm, where if i is None else "session %d %s" % (i, where))
        if key not in first:
            first[key] = (v, at)
            continue
        out["changed"].append("the host's %s changed during the run: %s (%s) -> %s (%s, %s)" % (
            key, first[key][0], first[key][1], v, at, sc.iso_utc(t) if t else "time unknown"))
    return out


# renormalised.json keys a re-evaluation states (the renormaliser's revision where it records one)
RENORMALISED_KEYS = ("source", "out", "at", "normaliser_sha256", "digest", "counts", "skipped_by_kind", "commit",
                     "git_head", "runner_sha256", "client_sha256", "input_sha256")


def renormalisation(arm_dir, rows=None):
    """What a re-normalisation of an arm's answers left in it, or None when there is no trace of one: `summary` (the
    keys of RENORMALISED_KEYS in renormalised.json, `changed` as 'SCENARIO rREPEAT'), `lines` (the executions marked
    `renormalised`: rebuilt, changed, skipped, and the changed ones) and `raw` (the raw records the lines name that
    keep the recorded answer as `answer_recorded`, of those read). The raw records are read for the MCP arm only:
    the re-normaliser rebuilds MCP answers from their tool payloads, a direct answer has nothing to rebuild from."""
    summary = sc.read_json(os.path.join(arm_dir, "renormalised.json"))
    if rows is None:
        rows, _ = sc.read_jsonl(os.path.join(arm_dir, "executions.jsonl"))
    lines = {"marked": 0, "rebuilt": 0, "changed": 0, "skipped": 0, "changed_executions": []}
    for r in rows:
        mark = r.get("renormalised")
        if not isinstance(mark, dict):
            continue
        lines["marked"] += 1
        if "skipped" in mark:
            lines["skipped"] += 1
            continue
        lines["rebuilt"] += 1
        if mark.get("changed"):
            lines["changed"] += 1
            lines["changed_executions"].append("%s r%s" % (r.get("scenario"), r.get("repeat")))
    raw = {"read": 0, "answer_recorded": 0}
    if os.path.basename(os.path.normpath(arm_dir)) == "mcp":
        for r in rows:
            rel = r.get("raw")
            if not isinstance(rel, str) or not rel:
                continue
            path = rel if os.path.isabs(rel) else os.path.join(arm_dir, rel)
            try:
                doc = sc.read_gz_json(path)
            except (OSError, ValueError, EOFError):
                continue
            raw["read"] += 1
            if isinstance(doc, dict) and "answer_recorded" in doc:
                raw["answer_recorded"] += 1
    if not isinstance(summary, dict) and not lines["marked"] and not raw["answer_recorded"]:
        return None
    out = {"lines": lines, "raw": raw}
    if isinstance(summary, dict):
        out["summary"] = dict((k, summary[k]) for k in RENORMALISED_KEYS if k in summary)
        out["summary"]["changed"] = ["%s r%s" % (c.get("scenario"), c.get("repeat"))
                                     for c in summary.get("changed") or [] if isinstance(c, dict)]
    return out


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
    manifest = sc.read_json(manifest_copy) if manifest_copy else sc.read_json(os.path.join(run_dir, "manifest.json"))
    scenarios = dict((x.get("id"), x) for x in (manifest or {}).get("scenarios", []) if isinstance(x, dict)) \
        if isinstance(manifest, dict) else {}
    rederived = stricter_outcomes(merged_rows, arm_runs, load_scenario_outcome(checks_dir), scenarios, run_dir)
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
    identity = manifest_identity(run_dir, arm_runs)
    if identity["differs"]:
        warnings.append("the arms ran different manifests: %s" % "; ".join(identity["differs"]))
    for problem in identity["manifest"]["problems"]:
        warnings.append("manifest.json: %s" % problem)
    releases = release_identity(arm_runs)
    warnings.extend(releases["changed"])
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
    renormalised = {}
    for name, path in arm_dirs(run_dir):
        info = renormalisation(path)
        if info is None:
            continue
        renormalised[name] = info
        if not reevaluated:
            warnings.append("arms/%s: its answers were re-normalised after the run (%d of %d line(s) rebuilt, %d "
                            "changed; %d raw record(s) keep the recorded answer): the outcomes of those lines are the "
                            "recorded ones, stale until a re-evaluation (run_suite.sh --reevaluate)" % (
                                name, info["lines"]["rebuilt"], arm_counts.get(name, 0), info["lines"]["changed"],
                                info["raw"]["answer_recorded"]))

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
        "manifest_identity": identity,
        "release_identity": releases,
        "mixed_releases": bool(releases["changed"]),
        "renormalisation": renormalised,
        "warnings": warnings,
    }
    if mark is not None:
        merged["reevaluation"] = dict((k, v) for k, v in mark.items() if k != "changes")
        merged["reevaluation"]["changes"] = len(mark.get("changes") or [])
        if renormalised:
            merged["reevaluation"]["renormalised"] = renormalised
    sc.write_json_atomic(os.path.join(run_dir, name_of("run.json")), merged)
    return merged


def identity_errors(merged):
    """The identity problems of a merged run.json that a gate reads as errors: the arms that ran another manifest than
    the run's, the run manifest's stated hashes its content does not match, and a host release or index_fp that changed
    during the run."""
    ident = merged.get("manifest_identity") or {}
    return list(ident.get("differs") or []) + ["manifest.json: %s" % p for p in
                                               ((ident.get("manifest") or {}).get("problems") or [])] + \
        list((merged.get("release_identity") or {}).get("changed") or [])


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    reevaluated = "--reevaluated" in argv
    strict = "--strict-identity" in argv
    argv = [a for a in argv if a not in ("--reevaluated", "--strict-identity")]
    if len(argv) != 1 or argv[0].startswith("-"):
        sys.stderr.write(__doc__)
        return 2
    merged = merge(argv[0], reevaluated)
    sys.stdout.write("merged %s: %d executions, %d samples, arms %s\n" % (
        argv[0], merged["counts"]["executions"], merged["counts"]["samples"], ", ".join(sorted(merged["arms"]))))
    for w in merged["warnings"]:
        sys.stdout.write("warning: %s\n" % w)
    errors = identity_errors(merged)
    if strict and errors:
        host = set((merged.get("release_identity") or {}).get("changed") or [])
        for e in errors:
            sys.stderr.write("merge_run: error: %s identity: %s\n" % ("host" if e in host else "manifest", e))
        return 3
    return 0


if __name__ == "__main__":
    sys.exit(main())
