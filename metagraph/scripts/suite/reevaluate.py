#!/usr/bin/env python3
"""Re-evaluates a recorded run OFFLINE with the current checks (no request to any host).

    reevaluate.py RUN_DIR [--checks DIR] [--by-design PATH] [--mcp-auth anonymous|api_key] [--quiet]

For every arm under RUN_DIR/arms/ and every execution line of its executions.jsonl (the last one per scenario,
arm and repeat, as merge_run keeps it), the raw answer the runner recorded (raw/...json.gz: both runners keep the
NORMALISED answer under `answer`, the MCP runner the shape of its normaliser, so no normaliser runs again) is
judged again with checks.evaluate of DIR/checks.py (default: this directory), with what the runner passed:
`previous` = the earlier repeats of the scenario on the arm the host answered and the runner did not end in error
(suite_common.feeds_previous of the recorded line, as both runners apply it), in repeat order, as {digest,
time_limited} computed from their raw answers with the same checks (the recorded digest and time_limited where a raw
answer is missing), and `client` = {"wall_ms"}: a direct arm's recorded wall time (run_direct.py), the MCP arm's
service-measured backend time (mcp_phases.backend_ms, as async's runner passes it since its 0d127f4a). The outcome is
checks.scenario_outcome of the new checks, except where the runner decided it without them, which is kept: an
outcome outside the vocabulary (suite_common.outcome_known; a warning, read as a failure downstream), `error`
(no answer, a transport or decode failure, a busy deployment; unless the recorded checks are the checks module's
own exception, which a re-evaluation repeats) and `not_implemented_mcp`. On the MCP arm the runner's verdict is
replayed with the current table (--by-design PATH, default compare.py's; '' none; async's runner at 31a2e41f,
`judge`, is the reference: e334715b with rules 1, 2 and 4 of 62ee4214, 2a8acbaf, f94b079d, rule 4's read floor of
f52846fb, the explained items of f386a3c0 and answer_match.mcp_outcome of 901893b7): the cases the scenario names
whose service_answer the answer matches (the field read from the recorded tool error), else the first
answer-triggered case it matches whose wall rule holds or cannot be decided, each only where its `request` holds on
the scenario's MCP arguments (the references it reads recorded as the line's request_evidence: the run's recorded
live documents, else the table's blocks; a table, or a manifest naming a case its own arguments do not satisfy, is
refused like the runner refuses it); when the first matched case states
answer_match.mcp_outcome not_implemented_mcp the execution is not_implemented_mcp (finding when a check is a
finding; the case recorded, the runner's verdict_detail; neither rule (c) nor a host-state rule applies, and a
runner's not_implemented_mcp mark that names its case is replayed, not kept); else a refused / relayed /
answer-triggered match is service_by_design, finding when a check is a finding (a by-design mark never erases
one); an answered / clamped match leaves the outcome to the checks (compare.explained_outcome, triage-2 B-01: fail
when a failing check item is not one the matched cases explain, else finding, service_by_design when some failing
item is explained, else the checks' outcome); the line records the case (`service_by_design_case`), the explained
checks (`explained_checks`, names) and items (`explained_items`, indices into checks: B5b-06) whenever a case
matched. A runner's service_by_design mark the current table does not reproduce is judged by the checks (B-06), and
says so; with no table a runner's mark, a by-design finding (B5b-10) and a rule (c) inconclusive on its case
(compare.runner_rule_c_by_design) are kept. Then the host-state rule (an
unexpected compare.MCP_HOST_STATE_CODES answer whose service-side wall the plausibility rule of
suite_common.host_state_verdict admits is inconclusive; one it rejects is judged by the checks, and so is a relayed
/resolve deadline on a discovery that sent no budget (rule 1) and an answer the client's own wall proves premature
(rule 2: ANONYMOUS MCP calls only, the line's or raw record's `mcp_auth`; a record before 2a8acbaf has no such key
and rule 2 decides nothing there unless --mcp-auth states the arm's auth, recorded in the mark and on each such line
as reevaluation.mcp_auth_assumed; a recorded null stays unknown, as the runner judged it); a resolve_timeout is
judged against the bound it states (rule 4; a stated read bound below the service's own floor for a known budget,
budget + reserve + 10 s - 1 s, is no host state whatever the wall: "service read bound B s below its own floor F
s", f52846fb); rule (c): where no threshold or wall is known, for a host state or an answer-triggered case's
wall_rule, the execution is inconclusive with verdict_detail "no threshold known / no service-side timing:
plausibility undecidable" (with an API key "; rule 2 skipped: ..."), whatever the checks say on the host-state
path, B5b-03), then a new exception of the checks (fail). The line records verdict_detail (the runner's text: why
a wall rule or a host state held or not) and inconclusive_reason as async's runner at 31a2e41f does; the runner's
recorded client wall (raw `client_wall`, line `client_wall_ms`) and timeout statement (`timeout_terms`) are read
first and cross-checked against the raw record (a mismatch is a warning: on the line's reevaluation.warnings and in
the mark).

Writes, keeping every original file:
  RUN_DIR/arms/<arm>/executions.reevaluated.jsonl   the re-evaluated lines (the recorded fields, with checks,
      outcome, answer_digest and time_limited replaced, and `reevaluation`: {checks_sha256, outcome_recorded,
      checks_changed: the checks whose outcome changed, added or dropped; mcp_auth_assumed, warnings,
      explained_items_recorded (the runner's recorded explained_items where the replay's differ) where they apply})
  RUN_DIR/reevaluated.json   the re-evaluated set's mark: the checks (directory, commit, `+dirty`, sha256 of
      checks.py and suite_common.py), the recorded checks (each arm's run.json), rule 2's auth (`mcp_auth`: the
      assumption, the lines by recorded value, the assumed lines, the values the MCP arm's run.json states; a line
      whose recorded value is none of those is a warning), counts, and every execution whose outcome or checks
      changed
Then `merge_run.py --reevaluated`, `compare.py --reevaluated`, `report.py --reevaluated` (or `run_suite.sh
--reevaluate RUN_DIR`) build the re-evaluated set's executions.reevaluated.jsonl, run.reevaluated.json,
compare.reevaluated.json and report.reevaluated.md / .pdf beside the recorded ones.
"""

import argparse
import hashlib
import inspect
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import compare  # noqa: E402
import merge_run  # noqa: E402
import suite_common as sc  # noqa: E402

MARK = "reevaluated.json"
# outcomes a runner decides without the checks (kept by a re-evaluation)
KEPT_OUTCOMES = ("error", "not_implemented_mcp", "service_by_design")
# async's MCP runner (suite_mcp_runner.py HOST_STATE_CODES): an environmental host state the scenario does not
# expect is inconclusive, whatever the checks say; with resolve_timeout (G23 as VERIFY-RP amends it: the service's
# resolve client gave up before the host's overrun answer under load), as compare.py reads them
MCP_HOST_STATE_CODES = compare.MCP_HOST_STATE_CODES


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
    if not os.path.exists(os.path.join(directory, "checks.py")):
        raise SystemExit("reevaluate: no checks.py in %s" % directory)
    try:
        return sc.load_checks_bound(directory, "suite_checks_reevaluation")
    except sc.SnapshotChanged as e:
        raise SystemExit("reevaluate: %s" % e)


def checks_identity(directory):
    commit, dirty = git_commit(directory)
    return {"dir": os.path.abspath(directory), "commit": commit, "dirty": dirty,
            "commit_label": (commit + ("+dirty" if dirty else "")) if commit else None,
            "checks_sha256": sha256_file(os.path.join(directory, "checks.py")),
            "suite_common_sha256": sha256_file(os.path.join(directory, "suite_common.py"))}


def arm_lines(path):
    """The arm's execution lines, the last one per scenario, arm and repeat, in first-seen order; (lines, bad). A
    warm-up's re-run for a resume (suite_common.is_rerun) is no execution of its key: left out here (reruns() keeps
    it), so the original execution is the one judged."""
    rows, bad = sc.read_jsonl(path)
    keyed, order = {}, []
    for row in sc.judged_lines(rows):
        key = sc.execution_key(row)
        if key not in keyed:
            order.append(key)
        keyed[key] = row
    return [keyed[k] for k in order], bad


def reruns(path):
    """The arm's warm-up re-runs for a resume, as recorded: the re-evaluated set carries them unchanged (timing
    evidence, never judged), so merge_run.py --reevaluated puts them next to their targets as the recorded merge does."""
    rows, _ = sc.read_jsonl(path)
    return [r for r in rows if sc.is_rerun(r)]


def raw_doc(arm_dir, line):
    """The recorded raw file of an execution line (a dict), or None."""
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
    return doc if isinstance(doc, dict) else None


def raw_answer(arm_dir, line):
    """The recorded normalised answer of an execution line, or None (no raw file, or none kept in it)."""
    answer = (raw_doc(arm_dir, line) or {}).get("answer")
    return answer if isinstance(answer, dict) else None


# the MCP tool's own error object behind an answer's error (compare.tool_error: the runner's raw_error)
tool_error = compare.tool_error


def client_of(arm, line):
    """What the arm's runner passed as checks.evaluate's `client` (module docstring), or None."""
    if arm == "mcp":
        v = (line.get("mcp_phases") or {}).get("backend_ms")
    else:
        v = line.get("wall_ms")
    return {"wall_ms": float(v)} if isinstance(v, (int, float)) and not isinstance(v, bool) else None


def evaluate_raised(checks):
    return bool(checks) and all(c.get("check") == "evaluate" for c in checks)


def _finding(checks):
    return any(str(c.get("outcome")) == "finding" for c in checks or [])


def outcome_of(checks_mod, arm, scenario, answer, line, new_checks, raised, cases=None, tool_err=None, info=None,
               doc=None):
    """The re-evaluated outcome (module docstring; async's judge at 31a2e41f) and why it was kept (decided without
    the checks), if it was; |info| (a dict) receives what the line records as async's runner records it:
    service_by_design_case and explained_checks (whenever a case matched), verdict_detail (why a wall rule or a host
    state held or not, or rule (c)'s text), inconclusive_reason (an inconclusive outcome's class: rule (c)'s text, the
    host-state code, suite_common.SCENARIO_LEVEL_MARK where the checks decided and the property-level marks cover
    every property check, with suite_common.scenario_level_mark_detail in verdict_detail, or "checks"); and for the
    re-evaluation's record by_design_cases, host_state and `note` (a mark the table does not reproduce, a matched case
    the checks decided) and `warnings` (a recorded client wall or timeout statement the raw record does not
    reproduce). |doc| is the execution's raw record (the service-side
    walls: suite_common.service_wall; rule 2's client wall: suite_common.client_wall_ms, its auth
    suite_common.mcp_auth_of, which reads |line|'s reevaluation.mcp_auth_assumed where the record states none). An
    undecidable host state is inconclusive whatever the checks say, as in the runner (B5b-03: lane S5b's
    addendum-11 divergence is dropped)."""
    info = {} if info is None else info
    outcome, kept = _outcome_of(checks_mod, arm, scenario, answer, line, new_checks, raised, cases, tool_err, info,
                                doc)
    if not info.get("warnings"):
        info.pop("warnings", None)
    if outcome == "inconclusive" and kept is None and "inconclusive_reason" not in info:
        whole = sc.every_property_marked(scenario, arm, new_checks)
        if whole:
            # the checks decided, and the marks cover every property check: a scenario-level mark, named as such
            info["inconclusive_reason"] = sc.SCENARIO_LEVEL_MARK
            info["verdict_detail"] = "; ".join(x for x in (info.get("verdict_detail"),
                                                           sc.scenario_level_mark_detail(whole)) if x)
    if outcome == "inconclusive":
        info.setdefault("inconclusive_reason", "checks")
    else:
        info.pop("inconclusive_reason", None)
    return outcome, kept


def _recorded_mark(line, info):
    """A runner's mark kept without a table to replay it: the line keeps the runner's own record of its case
    (service_by_design_case, explained_checks, explained_items, verdict_detail) instead of losing it."""
    for k in ("service_by_design_case", "explained_checks", "explained_items", "verdict_detail"):
        if k in line:
            info[k] = line[k]


def _outcome_of(checks_mod, arm, scenario, answer, line, new_checks, raised, cases, tool_err, info, doc):
    recorded = line.get("outcome")
    if not sc.outcome_known(recorded):
        # a runner's outcome outside the shared vocabulary is no verdict the checks may replace: kept, and read as
        # a failure by compare and the report
        info.setdefault("warnings", []).append("outcome %r outside the vocabulary (%s): kept" % (
            recorded, ", ".join(sc.OUTCOMES)))
        return recorded, "an outcome outside the vocabulary (%r): kept, read as a failure" % (recorded,)
    kept = recorded in KEPT_OUTCOMES and not (recorded == "error" and evaluate_raised(line.get("checks")))
    # a mark of a matched table case is replayed with the table: service_by_design, and (async's 901893b7) a
    # not_implemented_mcp that names its case (answer_match.mcp_outcome); the manifest's mcp.not_implemented names none
    replayed = recorded == "service_by_design" or (
        arm == "mcp" and recorded == "not_implemented_mcp" and bool(line.get("service_by_design_case")))
    if kept and not replayed:
        return recorded, "decided by the runner without the checks"
    if arm == "mcp" and not cases and recorded == "finding" and line.get("service_by_design_case"):
        # B5b-10: with no table a by-design finding (a relayed deadline's budget_respected) is the runner's mark as
        # its service_by_design is: kept
        _recorded_mark(line, info)
        return recorded, "decided by the runner's by-design case %s (a finding the mark never erases); no table " \
                         "to replay it" % line.get("service_by_design_case")
    if arm == "mcp" and not cases and compare.runner_rule_c_by_design(dict(line, arm=arm)):
        # with no table the runner's rule (c) on its case (the wall rule undecidable) is its decision, kept as its
        # service_by_design mark is: inconclusive, the case and the class recorded
        _recorded_mark(line, info)
        info["inconclusive_reason"] = line.get("inconclusive_reason")
        return recorded, "decided by the runner's by-design case %s by rule (c) (%s); no table to replay it" % (
            line.get("service_by_design_case"), line.get("inconclusive_reason"))
    rule = getattr(checks_mod, "scenario_outcome", None) or sc.outcome_from_checks
    notes = []
    warnings = info.setdefault("warnings", [])
    if arm == "mcp" and cases:
        # async's runner (judge at 31a2e41f): the named cases first, else an answer-triggered case whose request and
        # wall rule hold, or whose wall rule cannot be decided (rule (c)); a case its request or a wall holds back
        # (a known one, or rule 2) does not apply. The references the requests read: request_evidence
        walls = []
        evidence = compare.RequestEvidence()
        matched, how = compare.case_matches(cases, scenario, answer, tool_err, line, doc, walls, warnings, evidence)
        if evidence.sources:
            info["request_evidence"] = evidence.sources
        if evidence.mismatches:
            info["request_evidence_mismatches"] = evidence.mismatches
        notes = [d for state, d, _ in walls if state != "undecidable"]
        if matched:
            info["by_design_cases"] = [c.get("id") for c in matched]
            if compare.match_mark(matched[0], how) == "not_implemented_mcp":
                # async's 901893b7: the first matched case states that THIS deployment's configuration refuses the
                # request (answer_match.mcp_outcome): not implemented via MCP here, the case recorded, the runner's
                # verdict_detail (the wall-rule notes first); neither rule (c) nor a host-state rule applies, and a
                # finding is never erased
                cid = matched[0].get("id")
                info.update(service_by_design_case=cid, explained_checks=[],
                            verdict_detail=compare.not_implemented_detail(matched[0], answer, tool_err, notes))
                if _finding(new_checks):
                    return "finding", "the table's %s case %s (answer_match.mcp_outcome not_implemented_mcp); a " \
                                      "check is a finding (a mark never erases one)" % (how, cid)
                return "not_implemented_mcp", "the table's %s case %s (answer_match.mcp_outcome " \
                                              "not_implemented_mcp)" % (how, cid)
            if walls and walls[-1][0] == "undecidable":
                # rule (c): the case's answer matched, its wall rule cannot be decided: inconclusive, the case
                # recorded; a finding is never erased
                cid, why = matched[0].get("id"), walls[-1][1]
                info.update(service_by_design_case=cid, explained_checks=[], verdict_detail="; ".join(notes + [why]))
                if _finding(new_checks):
                    return "finding", "the table's %s case %s, its wall rule undecidable; a check is a finding" % (
                        how, cid)
                info["inconclusive_reason"] = sc.undecidable_class(why)
                return "inconclusive", "the table's %s case %s: %s (rule (c))" % (how, cid, why)
            if compare.judged_by_checks(matched[0]):
                # the checks a property-level mark sets aside (suite_common.unmarked_checks) decide nothing
                outcome, explained, noted = compare.answered_verdict(
                    matched, sc.unmarked_checks(scenario, arm, new_checks), rule)
                info.update(service_by_design_case=noted, explained_checks=explained,
                            explained_items=sc.explained_items(new_checks, compare.union_of(matched)))
                info["note"] = "the table's %s case %s: judged by the checks%s" % (
                    how, noted, ", %s explained" % explained if explained else "")
                return outcome, None
            cid = matched[0].get("id")
            info.update(service_by_design_case=cid, explained_checks=[])
            if notes:
                info["verdict_detail"] = "; ".join(notes)
            if _finding(new_checks):
                return "finding", "the table's %s case %s; a check is a finding (a by-design mark never erases " \
                                  "one)" % (how, cid)
            return "service_by_design", "the table's %s case %s" % (how, cid)
    if kept:
        if not cases:
            if arm == "mcp":
                _recorded_mark(line, info)
            return recorded, "decided by the runner without the checks (no table to replay its mark)"
        info["note"] = "the runner's %s mark is not reproduced by the current table: judged by the checks" % recorded
    if arm == "mcp":
        code, state, detail = sc.host_state_verdict(scenario, line, doc, (answer or {}).get("error"), tool_err,
                                                    warnings, compare.refs_of(cases))
        if detail:
            notes.append(detail)
        if notes:
            info["verdict_detail"] = "; ".join(notes)
        if state == "host_state":
            info.update(host_state=code, inconclusive_reason=code)
            return "inconclusive", "the MCP runner's host-state rule (%s)" % code
        if state == "undecidable":
            # rule (c), as the runner: inconclusive whatever the checks say (B5b-03); the class without rule 2's
            # API-key note
            info.update(host_state=code, inconclusive_reason=sc.undecidable_class(detail))
            return "inconclusive", "the MCP runner's host-state rule (%s): %s (rule (c))" % (code, detail)
    elif notes:
        info["verdict_detail"] = "; ".join(notes)
    if raised:
        return ("fail" if arm == "mcp" else "error"), None
    # an MCP check the manifest marks as not shown via MCP (mcp.not_implemented_checks / by_design_checks) that reads
    # not_exposed is set aside: it makes the execution neither inconclusive nor untested, the rest decides
    return rule(sc.unmarked_checks(scenario, arm, new_checks)), None


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


def reevaluate(run_dir, checks_dir=None, log=None, by_design=None, mcp_auth=None):
    """Writes the re-evaluated set (module docstring); returns the mark (reevaluated.json's content). |mcp_auth|
    ("anonymous" | "api_key", --mcp-auth): what rule 2 assumes for an MCP line whose record states no mcp_auth (a
    run before async's 2a8acbaf), recorded on each such line (`reevaluation.mcp_auth_assumed`) and in the mark."""
    if mcp_auth is not None and mcp_auth not in sc.MCP_AUTH_VALUES:
        raise ValueError("--mcp-auth takes %s, not %r" % (" or ".join(sc.MCP_AUTH_VALUES), mcp_auth))
    checks_dir = os.path.abspath(checks_dir or HERE)
    top_manifest = sc.read_json(os.path.join(run_dir, "manifest.json"))
    if not isinstance(top_manifest, dict):
        top_manifest = next((m for m in (sc.read_json(os.path.join(d, "manifest.json"))
                                         for _, d in merge_run.arm_dirs(run_dir)) if isinstance(m, dict)), None)
    # the table, refused like the runner refuses it, also when a scenario names a case its own arguments do not
    # satisfy; bound to the run's recorded live documents (compare.bind_run)
    cases, by_design_info = compare.load_by_design(compare.DEFAULT_BY_DESIGN if by_design is None else by_design,
                                                   (top_manifest or {}).get("scenarios") or [])
    if by_design_info is not None:
        cases = compare.bind_run(cases, run_dir)
    checks_mod = load_checks(checks_dir)
    ident = checks_identity(checks_dir)
    takes = inspect.signature(checks_mod.evaluate).parameters
    takes_prev, takes_client = "previous" in takes, "client" in takes
    # a checks module whose time_limited reads the scenario, the client and the arm (the front ingress's 504 at its
    # cut, VERIFY-V6-2 L1) gets them, as run_direct.time_limited_of passes them
    tl_takes_arm = hasattr(checks_mod, "time_limited") and \
        "arm" in inspect.signature(checks_mod.time_limited).parameters
    short = (ident["checks_sha256"] or "")[:16]
    mark = {"run_dir": os.path.abspath(run_dir), "reevaluated_utc": sc.iso_utc(time.time()), "checks": ident,
            "client_passed": {"direct-local": "wall_ms", "direct-remote": "wall_ms", "mcp": "mcp_phases.backend_ms"},
            "by_design": by_design_info, "recorded": {}, "arms": {}, "changes": [], "warnings": [],
            # where the references the table's requests and the plausibility rules read came from, and every live /
            # table mismatch (judged with the live value)
            "references": dict(compare.reference_sources(cases), mismatches=compare.reference_mismatches(cases))
            if by_design_info is not None else None,
            # rule 2's auth (async's 2a8acbaf): the MCP lines by the mcp_auth their record states, and the assumption
            # --mcp-auth made for those that state none (null: none made; they stay unknown, rule 2 decides nothing);
            # run_json: the values the MCP arm's run.json states (its invocations' included)
            "mcp_auth": {"assumed": mcp_auth, "recorded": {}, "assumed_lines": 0, "run_json": []}}
    for name, arm_dir in merge_run.arm_dirs(run_dir):
        run = sc.read_json(os.path.join(arm_dir, "run.json"), {}) or {}
        mark["recorded"][name] = {"checks_module": run.get("checks_module"),
                                  "scripts_commit": run.get("scripts_commit") or run.get("git_head")}
        # the MCP runner's run.json states the client's mcp_auth (2a8acbaf: the latest invocation's at the top, each
        # one's under invocations, as mcp_runs.jsonl): recorded in the mark, and a line's own value cross-checked
        run_auth = []
        if name == "mcp":
            for meta in [run] + [m for m in run.get("invocations") or [] if isinstance(m, dict)]:
                if meta.get("mcp_auth") is not None and meta["mcp_auth"] not in run_auth:
                    run_auth.append(meta["mcp_auth"])
            mark["mcp_auth"]["run_json"] = run_auth
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
                                    "time_limited": (checks_mod.time_limited(
                                        answer, scen.get(line.get("scenario")), client_of(name, line), name)
                                        if tl_takes_arm else checks_mod.time_limited(answer))
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
            # an earlier repeat without the host's answer (no HTTP status) or ended in error by its runner holds a
            # stand-in's digest: no earlier repeat (sc.feeds_previous, as run_direct.py and async's runner)
            previous = [entries[sc.execution_key(o)]
                        for o in sorted(repeats.get(sid, []), key=lambda o: o.get("repeat") or 0)
                        if (o.get("repeat") or 0) < rep and sc.execution_key(o) in entries and sc.feeds_previous(o)]
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
            doc = raw_doc(arm_dir, line) if name == "mcp" else None
            terr = tool_error(doc, answer) if name == "mcp" else None
            judged, assumed = line, None
            if name == "mcp":
                auth, _ = sc.mcp_auth_of(line, doc)
                stated = "unknown" if auth is None else str(auth)
                mark["mcp_auth"]["recorded"][stated] = mark["mcp_auth"]["recorded"].get(stated, 0) + 1
                if auth is not None and run_auth and auth not in run_auth:
                    mark["warnings"].append("arms/%s %s r%s: mcp_auth %r is none of run.json's %s" % (
                        name, sid, rep, auth, run_auth))
                if auth is None and mcp_auth is not None and not sc.mcp_auth_recorded(line, doc):
                    # a record before 2a8acbaf (no mcp_auth key): rule 2 reads the stated assumption
                    # (suite_common.mcp_auth_of); a recorded null stays unknown, as the runner judged it
                    assumed = mcp_auth
                    judged = dict(line, reevaluation={"mcp_auth_assumed": mcp_auth})
                    mark["mcp_auth"]["assumed_lines"] += 1
            info = {}
            outcome, kept = outcome_of(checks_mod, name, scenario, answer, judged, new_checks, raised, cases, terr,
                                       info, doc)
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
            if assumed is not None:
                row["reevaluation"]["mcp_auth_assumed"] = assumed
            for w in info.get("warnings") or []:
                # a recorded client wall or timeout statement the raw record does not reproduce (S5c-03)
                row["reevaluation"].setdefault("warnings", []).append(w)
                mark["warnings"].append("arms/%s %s r%s: %s" % (name, sid, rep, w))
            # as async's runner records them: the matched case and the failing checks it explains, whenever a case
            # matched (the replay's, not the recorded run's)
            for k in ("service_by_design_case", "explained_checks", "explained_items", "verdict_detail",
                      "inconclusive_reason", "request_evidence", "request_evidence_mismatches"):
                row.pop(k, None)
            if "service_by_design_case" in info:
                row["service_by_design_case"] = info["service_by_design_case"]
                row["explained_checks"] = info.get("explained_checks") or []
                # B5b-06: the explained check ITEMS (indices into checks), as the runner's E is keyed (its line and
                # raw record carry the same field since f386a3c0: [] when E is empty, absent when no case matched)
                row["explained_items"] = info.get("explained_items") or []
            if "explained_items" in line and line.get("explained_items") != row.get("explained_items"):
                # the runner's own list (f386a3c0) that the replay does not reproduce (another table, other checks):
                # noted, not trusted (compare notes a line's differing list the same way, S5c-05)
                row["reevaluation"]["explained_items_recorded"] = line.get("explained_items")
            for k in ("verdict_detail", "inconclusive_reason", "request_evidence", "request_evidence_mismatches"):
                if info.get(k):
                    row[k] = info[k]
            for k in ("by_design_cases", "host_state", "note"):
                if info.get(k):
                    row["reevaluation"][k] = info[k]
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
        again = reruns(os.path.join(arm_dir, "executions.jsonl"))
        if again:
            counts["reruns_for_resume"] = len(again)
        for line in again:
            row = dict(line)
            row["reevaluation"] = {"checks_sha256": short, "outcome_recorded": line.get("outcome"),
                                   "checks_changed": [],
                                   "kept": "a warm-up's re-run for a resume: timing evidence, not judged"}
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
    ap.add_argument("--by-design", default=None,
                    help="the by-design table the MCP runner's verdict is replayed with (default: compare.py's; '' none)")
    ap.add_argument("--mcp-auth", choices=sc.MCP_AUTH_VALUES, default=None,
                    help="rule 2's auth for MCP lines whose record states no mcp_auth (runs before async's 2a8acbaf): "
                         "anonymous lets the client wall bound an unknown service-side wall; recorded in the mark "
                         "and on each such line (default: none assumed, rule 2 decides nothing there)")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args(argv)
    if not os.path.isdir(os.path.join(args.run_dir, "arms")):
        sys.stderr.write("reevaluate: %s has no arms/ directory\n" % args.run_dir)
        return 2
    log = None if args.quiet else (lambda m: sys.stdout.write(m + "\n"))
    try:
        mark = reevaluate(args.run_dir, args.checks, log, args.by_design, args.mcp_auth)
    except compare.TableError as exc:
        sys.stderr.write("reevaluate: refusing: %s\n" % exc)
        return 4
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
