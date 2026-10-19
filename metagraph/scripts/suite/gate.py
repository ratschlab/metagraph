#!/usr/bin/env python3
"""The strict deploy gate of a suite run, and the run's persisted plan.

    gate.py RUN_DIR POLICY [--reevaluated] [--out FILE] [--max-lines N] [--expect-service-commit SHA]
                           [--expect-library-commit SHA] [--reuse]
    gate.py plan RUN_DIR --manifest M --arms A,B [--areas P,M] [--limit N] [--ids LIST] [--warm-every N] [--mode M]
    gate.py ids FILE --manifest M          print FILE's scenario ids, comma-separated in manifest order
    gate.py policy-ids POLICY [--arms A]  print the ids the policy's arms run (its selection less their exclusions)
    gate.py resume-ids RUN_DIR --manifest M  print the scenarios of the run's whole plan (a plain resume's selection)

Exit status of the judgement: 0 GATE PASS, 1 GATE FAIL, 2 usage (an unreadable policy or run included).

The plan. run_suite.sh writes RUN_DIR/plan.json before any arm runs: per invocation the selection (manifest identity,
areas, limit, ids, warm_every, arms, mode) and every (scenario, arm, repeat) it plans, the paired loop's plan
(suite_common.plan_executions for a direct arm, suite_common.mcp_repeats for the MCP arm; sequential mode runs the
same executions). An invocation whose plan equals the last one adds nothing; another selection on --resume is added,
and the run's plan is the union of all of them: a resume never narrows what the run has to record. A run of before
plan.json is reconstructed from orchestrator.json's sessions (the same selection fields, persisted before its arms
ran); without either, the policy's required executions are the plan.

The policy (gate_policy.py writes it from a manifest, the by-design table and the service's MCP surface): the
manifest identity (content_sha256 and relations_sha256, index_fp), the selected scenarios (by route), the arms the
gate requires, the scenarios an arm cannot express as explicit exclusions (each with its rule and reason, never
inferred from an answer), the executions required per arm, the relations required per arm, the property checks
explicitly not required (none by default) and whether the MCP arm must agree with a direct arm.

GATE PASS only when every one of these holds (each failure is listed, grouped by kind):
- identity: the run manifest's content_sha256 and relations_sha256 are the policy's, its stated hashes verify, every
  arm ran it (merge_run.manifest_identity), every plan of the run names it, every policy arm observed the index_fp
  the policy pins (and release / index_meta_fp where it pins them) and no arm observed another
  (host_identity_problems: the MCP arm's from service_capabilities.host_document), the host's release and index_fp
  did not change during the run (merge_run.release_identity), every arm of a recorded set was judged with the run's
  snapshot of its checks, every MCP invocation stating both checks_sha256 and suite_common_sha256, and a recorded
  set without a snapshot blocks (merge_run.checks_identity), the checks that judged the set are the ones their
  snapshot's SNAPSHOT.json records where they have one (suite_common.snapshot_refusal: a changed file is never
  executed), the set's comparison was judged with the checks recorded for the set (SNAPSHOT.json, a re-evaluated
  set's reevaluated.json: compared_with_checks_problems) and read the very lines the gate judges (its
  executions_sha256), and the
  comparison read the by-design table the policy's exclusions were derived from; the deployed service's identity
  (service_identity_problems, over merge_run.service_identity, the commit the MCP runner read before each
  invocation): given the deployed commit (--expect-service-commit SHA, run-gate's LIVE line), every MCP invocation
  read that commit (a differing one, or none: not confirmed, blocks; a run recording no identity at all blocks), and
  in any case the commit did not change during an invocation (identity_changed true) or between invocations (more
  than one commit); identity_changed null with the commit confirmed, and without a given commit a run recording no
  identity (an older runner), are notes;
- missing: every planned execution and every execution the policy requires is recorded;
- outcome: no execution of the policy's arms and selection is fail, error, finding, inconclusive (ANY inconclusive:
  the owner's strict gate), not_implemented_mcp (what the MCP cannot express is excluded in the policy), or
  service_by_design for a case the scenario does not name in the manifest, or a by-design pair not verified (the
  pair verified when a direct arm ran the scenario, the MCP side verified when none did);
- property: on a passing execution every property check the scenario expects (checks.is_property) passed, unless the
  policy lists it as not required for that scenario and arm, or the named answered / clamped case of a by-design
  execution explains the item (its explained_items: a failure, or a not_exposed reading whose detail an
  {check, detail_matches} entry matches, suite_common.item_explained; the pair judged below), or (the MCP arm) the
  manifest marks it not shown via MCP
  (suite_common.mcp_check_marks: mcp.not_implemented_checks / by_design_checks) and it read not_exposed
  (suite_common.checks_marked_aside; a marked check with no result is a missing result as any other); a passing MCP
  execution whose every property check is marked (suite_common.every_property_marked: a scenario-level mark) is a
  problem; with require_agreement every required property was compared with the direct reference (a marked one read
  not_exposed is neither compared nor missing);
- disagreement: no MCP disagreement (MCP a policy arm) and no direct disagreement (both direct arms policy arms) in the
  selection; with require_agreement the MCP class is agree or a verified by-design pair;
- relation: every relation the policy requires passed on each of its arms;
- configuration (the MCP arm a policy arm): configuration_problems, i.e. compare.gate_configuration_problems
  (a live document with service_settings, no live / table mismatch, no reference read from the table or unknown)
  plus a read by every invocation, one document per run (each record's identity recomputed, a recorded sha256 its
  record does not reproduce a problem), the runner's recorded mismatches, the advanced_path premise checked live and
  well-formed request_evidence; the per-line evidence is one problem per reference path;
- something passed: on each policy arm, an execution of the selection (a gate whose MCP selection is all excluded
  verified nothing via MCP, whatever a direct arm passed);
- area G (a policy with a `graphlets` block, manifest v8; gate_policy.graphlet_block), over the run's graphlet
  comparison, recomputed from the run's files every time (compare_graphlets.load_or_compute; with --reuse,
  graphlet_compare.json when every input it read is current: records, documents, twin, job and identity records,
  export files, the by-design table, the code), kind "graphlet" unless named otherwise: the comparison is of the policy's manifest (identity); every op
  call the policy requires has an entry and every way it requires a record (missing); its outcome is agree or a
  verified service_by_design pair, W3's not_exposed per the policy being no gap: a disagreement, an inconclusive or a
  precondition_failed blocks, unless the declared file sets the op call aside (graphlets.not_required, with its
  reason: a note); every pre-check of every graphlet scenario passed (library_commit not_exposed without a W3 call);
  no service check failed, and every one the policy requires passed (graphlet_service_checks.
  blocking_service_checks; a check the declared file sets aside is a note); the drift lane is a note. The deployed
  library's identity (library_identity_problems over merge_run.library_identity, the library_commit the MCP runner read per
  invocation), checked like the service's: given --expect-library-commit, every invocation read it and the policy's
  pin is it; with a graphlets block, every invocation's commit is the pin; in any case more than one commit blocks.
  Without a graphlets block and without --expect-library-commit nothing of this is judged or stated (a v4..v7
  verdict is unchanged).
"""

import argparse
import hashlib
import json
import os
import re
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import compare as cmp  # noqa: E402
import merge_run  # noqa: E402
import suite_common as sc  # noqa: E402

PLAN_FILE = "plan.json"
PLAN_VERSION = 1
POLICY_VERSIONS = (1,)
# A gate problem's kinds, in the order they are listed.
KINDS = ("identity", "configuration", "missing", "outcome", "property", "disagreement", "agreement", "relation",
         "graphlet", "nothing_passed")
# The check names of checks.GENERIC_CHECKS, for a checks.py that does not load: every other check is a property.
GENERIC_CHECKS = sc.GENERIC_CHECKS
# request_evidence sources that are no evidence of the deployment under the gate.
NOT_LIVE_EVIDENCE = ("table", "unknown")
# the outcomes a by-design or excluded execution may have besides pass
OUTCOME_BLOCKS = ("fail", "error", "finding", "inconclusive", "not_implemented_mcp")


# ---------------------------------------------------------------------------------------------------------------
# selection and plan

def read_ids(path, manifest):
    """The scenario ids of FILE (separated by commas, blanks or new lines; '#' starts a comment), in manifest order.
    An id the manifest lacks is refused (ValueError naming them): a typo must never drop a scenario silently."""
    with open(path) as fh:
        text = "\n".join(line.split("#", 1)[0] for line in fh)
    wanted = [t for t in text.replace(",", " ").split() if t]
    known = [s.get("id") for s in (manifest or {}).get("scenarios", [])]
    unknown = sorted(set(wanted) - set(known))
    if unknown:
        raise ValueError("%d id(s) not in the manifest: %s" % (len(unknown), ", ".join(unknown[:20])))
    if not wanted:
        raise ValueError("no scenario id in %s" % path)
    keep = set(wanted)
    return [sid for sid in known if sid in keep]


def planned_executions(manifest, areas=None, limit=None, ids=None, warm_every=10, arms=sc.ARMS):
    """[(scenario, arm, repeat)] a selection plans, in manifest order: the paired loop's executions (a direct arm:
    suite_common.plan_executions, the long class not on direct-remote; the MCP arm: suite_common.mcp_repeats)."""
    out = []
    for index, scenario in sc.select_scenarios(manifest, areas, limit, ids):
        sid = scenario.get("id")
        for arm in sc.ARMS:
            if arm not in arms:
                continue
            if arm == "mcp":
                reps = sc.mcp_repeats(index, scenario, warm_every)
            else:
                reps = [p[2] for p in sc.plan_executions([(index, scenario)], arm, warm_every)[0]]
            out.extend((sid, arm, r) for r in reps)
    return out


def _split(text):
    return [x.strip() for x in str(text or "").split(",") if x.strip()]


def _int_or(text, default):
    try:
        return int(str(text).strip()) if str(text if text is not None else "").strip() else default
    except ValueError:
        return default


def plan_entry(manifest, areas, limit, ids, warm_every, arms, mode, session=None):
    executions = planned_executions(manifest, areas, limit, ids, warm_every, arms)
    counts = {}
    for _, arm, _ in executions:
        counts[arm] = counts.get(arm, 0) + 1
    return {"written_utc": sc.iso_utc(time.time()), "session": session, "mode": mode,
            "manifest": sc.manifest_identity(manifest),
            "selection": {"areas": areas, "limit": limit, "ids": ids, "warm_every": warm_every, "arms": list(arms)},
            "counts": {"scenarios": len(set(e[0] for e in executions)), "executions": counts},
            "executions": [list(e) for e in executions]}


def write_plan(run_dir, entry):
    """Adds |entry| to RUN_DIR/plan.json unless the plans there, on the same manifest, already hold every execution
    it plans (a resume of the run's own selection adds nothing); returns whether it was added."""
    path = os.path.join(run_dir, PLAN_FILE)
    data = sc.read_json(path) or {"plan_version": PLAN_VERSION, "plans": []}
    held = set()
    for p in data.get("plans") or []:
        if p.get("manifest") == entry.get("manifest"):
            held.update(tuple(e) for e in p.get("executions") or [])
    if data.get("plans") and all(tuple(e) in held for e in entry.get("executions") or []):
        return False
    data["plans"].append(entry)
    sc.write_json_atomic(path, data)
    return True


def session_plan(manifest, session):
    """The executions an orchestrator.json session planned (for a run of before plan.json): its areas, limit, ids,
    warm_every and arms (the direct arms it did not skip; the MCP arm with an MCP command or a copied MCP run)."""
    arms = [arm for arm, skip in (("direct-local", "skip_local"), ("direct-remote", "skip_remote"))
            if str(session.get(skip) or "0") != "1"]
    if session.get("mcp_cmd") or session.get("mcp_run"):
        arms.append("mcp")
    ids = session.get("ids")
    ids = ids if isinstance(ids, list) else (_split(ids) or None)
    return planned_executions(manifest, _split(session.get("areas")) or None, _int_or(session.get("limit"), None),
                              ids, _int_or(session.get("warm_every"), 10), arms)


def run_plan(run_dir, manifest):
    """(keys, source, plans): the union of every plan of the run as a set of (scenario, arm, repeat), where it
    comes from ("plan.json", "orchestrator.json (reconstructed)", both, or None) and plan.json's entries. A session
    of orchestrator.json without an entry in plan.json (a run started before plan.json, resumed after it) is
    reconstructed and joins the union."""
    data = sc.read_json(os.path.join(run_dir, PLAN_FILE))
    plans = data.get("plans") or [] if isinstance(data, dict) else []
    keys, sources = set(), []
    for p in plans:
        keys.update(tuple(e) for e in p.get("executions") or [])
    if plans:
        sources.append(PLAN_FILE)
    orch = sc.read_json(os.path.join(run_dir, "orchestrator.json")) or {}
    sessions = orch.get("sessions") or []
    covered = set(p.get("session") for p in plans if isinstance(p.get("session"), int))
    rest = [s for i, s in enumerate(sessions) if i not in covered] if plans else sessions
    if plans and covered == set() and rest:
        rest = []   # entries without session numbers (written by hand or by a test): plan.json alone
    if rest and isinstance(manifest, dict):
        for s in rest:
            keys.update(session_plan(manifest, s))
        sources.append("orchestrator.json (reconstructed)")
    return keys, " and ".join(sources) or None, plans


def resume_ids(run_dir, manifest):
    """The scenarios of the run's whole plan (run_plan's union), in manifest order: the selection a resume without
    a selection of its own continues, so a plain resume after a narrower targeted one finishes the run."""
    keys, _, _ = run_plan(run_dir, manifest)
    wanted = set(k[0] for k in keys)
    return [s.get("id") for s in (manifest or {}).get("scenarios", []) if s.get("id") in wanted]


def missing_planned(run_dir, manifest, rows):
    """([(scenario, arm, repeat)] planned and not recorded, sorted in manifest order; the plan's source)."""
    keys, source, _ = run_plan(run_dir, manifest)
    recorded = set(sc.execution_key(r) for r in rows)
    pos = dict((s.get("id"), i) for i, s in enumerate((manifest or {}).get("scenarios", [])))
    arm_pos = dict((a, i) for i, a in enumerate(sc.ARMS))
    missing = sorted((k for k in keys if k not in recorded),
                     key=lambda k: (pos.get(k[0], len(pos)), arm_pos.get(k[1], 9), k[2]))
    return missing, source


# ---------------------------------------------------------------------------------------------------------------
# the policy

def load_policy(path):
    policy = sc.read_json(path)
    if not isinstance(policy, dict):
        raise ValueError("%s is not a JSON object" % path)
    if policy.get("policy_version") not in POLICY_VERSIONS:
        raise ValueError("%s: policy_version %r is not one of %s" % (path, policy.get("policy_version"),
                                                                    POLICY_VERSIONS))
    for key in ("manifest", "selection", "arms", "required"):
        if key not in policy:
            raise ValueError("%s: no %s" % (path, key))
    return policy


def excluded(policy, sid, arm, repeat=None):
    """The policy's exclusion of |sid| on |arm| (a dict), or None: of the whole scenario, or with |repeat| also one
    that excludes only some repeats (a declared entry's `repeats`)."""
    for ex in policy.get("exclusions") or []:
        if ex.get("scenario") != sid or arm not in (ex.get("arms") or sc.ARMS):
            continue
        if "repeats" not in ex or (repeat is not None and repeat in ex["repeats"]):
            return ex
    return None


def policy_ids(policy, arms=None):
    """The scenarios the policy's arms (or |arms|) run, in the policy's order: the selection less the scenarios
    excluded on every one of those arms."""
    arms = arms or policy.get("arms") or []
    return [sid for sid in policy["selection"].get("scenarios") or []
            if any(excluded(policy, sid, a) is None for a in arms)]


def required_keys(policy):
    """{(scenario, arm, repeat)} the policy requires."""
    out = set()
    for item in policy.get("required") or []:
        for r in item.get("repeats") or []:
            out.add((item.get("scenario"), item.get("arm"), r))
    return out


def not_required(policy, sid, arm, check):
    for nr in (policy.get("properties") or {}).get("not_required") or []:
        if nr.get("check") == check and nr.get("scenario") in (sid, "*") and arm in (nr.get("arms") or sc.ARMS):
            return nr
    return None


# ---------------------------------------------------------------------------------------------------------------
# the judgement

def _problem(kind, detail, scenario=None, arm=None, repeat=None):
    return {"kind": kind, "scenario": scenario, "arm": arm, "repeat": repeat, "detail": detail}


def _is_property_fn(checks_dir=None):
    mod = cmp.load_checks_module(checks_dir)
    fn = getattr(mod, "is_property", None) if mod is not None else None
    return fn if callable(fn) else (lambda name: name not in GENERIC_CHECKS)


# The checks a repeat evaluates only when an earlier one exists (checks.evaluate: determinism has nothing to compare
# on the first repeat, so it has no result there).
LATER_REPEATS_ONLY = sc.LATER_REPEATS_ONLY


def property_results(scenario, line, is_property, generic=False):
    """[(check name, occurrence, outcome or None)] of the property checks the scenario expects (every check with
    |generic|), the k-th expectation of a name answered by the line's k-th result of that name (None: no result)."""
    results = {}
    for c in line.get("checks") or []:
        results.setdefault(c.get("check"), []).append(c)
    seen, out = {}, []
    for e in (scenario or {}).get("expect") or []:
        name = e.get("check") if isinstance(e, dict) else None
        if not name or not (generic or is_property(name)):
            continue
        k = seen.get(name, 0)
        seen[name] = k + 1
        got = results.get(name) or []
        out.append((name, k, got[k].get("outcome") if k < len(got) else None))
    return out


def table_explanations(comp):
    """{(scenario, repeat): (explained item indices, the case)} of the MCP executions an answered / clamped case
    judged by the checks, as compare derived them from the table (mcp_agreement explained_items and by_design_notes,
    other_repeats included). The gate reads an execution's explanation from here, never from the line's own
    explained_items / explained_checks (a runner's E is not trusted: S5c-05)."""
    out = {}
    for sid, entry in ((comp or {}).get("scenarios") or {}).items():
        ma = (entry or {}).get("mcp_agreement") or {}
        others = ma.get("other_repeats")
        for part in [ma] + (others if isinstance(others, list) else [others]):
            if not isinstance(part, dict):
                continue
            notes = part.get("by_design_notes") or {}
            for rep, items in (part.get("explained_items") or {}).items():
                try:
                    out[(sid, int(rep))] = (list(items or []), notes.get(rep))
                except (TypeError, ValueError):
                    continue
    return out


def _explained_property(line, name, k, explanations):
    """Whether the table explains the k-th item of check |name| of an MCP line (table_explanations: indices into the
    line's checks); never on a direct arm, never by the line's own claim."""
    if line.get("arm") != "mcp":
        return False
    items, _ = explanations.get((line.get("scenario"), line.get("repeat")), ([], None))
    seen = {}
    for i, c in enumerate(line.get("checks") or []):
        n = seen.get(c.get("check"), 0)
        seen[c.get("check")] = n + 1
        if c.get("check") == name and n == k:
            return i in items
    return False


EVIDENCE_SOURCES = ("live", "table", "unknown")


def host_identity_problems(policy, arm_runs):
    """The host's identity against the policy's pins (each of suite_common.HOST_IDENTITY_KEYS its manifest block
    states: index_fp; release, index_meta_fp where pinned): every policy arm observed the key
    (merge_run.observed_host_identity: a direct arm's run.json and sessions, the MCP arm's service_capabilities
    host_document, top level and every invocation), and every value any arm observed is the pin. A policy arm that
    observed none is a problem as a conflicting value is: the pin is never taken as met unseen."""
    out = []
    want = policy.get("manifest") or {}
    arms = list(policy.get("arms") or [])
    every = sorted(set(arms) | set(arm_runs), key=lambda a: (sc.ARMS.index(a) if a in sc.ARMS else 9, a))
    for key in sc.HOST_IDENTITY_KEYS:
        pin = want.get(key)
        if not pin:
            continue
        for arm in every:
            seen = merge_run.observed_host_identity(arm, arm_runs.get(arm), key)
            if not seen:
                if arm in arms:
                    out.append(_problem("identity", "arms/%s states no observed %s; the policy pins %s (%s)" % (
                        arm, key, str(pin)[:16], "service_capabilities.host_document" if arm == "mcp"
                        else "run.json and its sessions"), arm=arm))
            elif seen != {pin}:
                out.append(_problem("identity", "arms/%s observed %s %s; the policy pins %s" % (
                    arm, key, ", ".join(sorted(x[:16] for x in seen)), str(pin)[:16]), arm=arm))
    return out


def _document_sha(rec):
    """A service_capabilities record's identity, recomputed: suite_common.service_document_sha256 (the configuration
    parts and the host document only: probe times and counters differ between the invocations of one deployment).
    The sha256 a record states is checked against it (configuration_problems), never read in its place."""
    return sc.service_document_sha256(rec)


def _snapshot_document_sha(run_dir, sha):
    """service_document_sha256 of the run's snapshot's suite_common.py, the definition its MCP runner computed the
    recorded sha256 with, when that file is still the one SNAPSHOT.json recorded (|sha|) and loads; None otherwise."""
    import importlib.util
    path = os.path.join(sc.snapshot_dir(run_dir), "suite_common.py")
    if not sha or sc.file_sha256(path) != sha:
        return None
    try:
        spec = importlib.util.spec_from_file_location("suite_common_of_the_run", path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
    except Exception:   # a snapshot that does not load: the record is compared recomputed with this definition
        return None
    fn = getattr(mod, "service_document_sha256", None)
    return fn if callable(fn) else None


def configuration_problems(run_dir, rows, policy, reevaluated=False, notes=None):
    """The MCP service's configuration evidence (agreed with async's runner on 2026-10-14), with the MCP arm gated:
    compare.gate_configuration_problems (lane G2's definition: a recorded live document stating service_settings,
    each referenced path live and equal to the table's block, no request_evidence read from the table or unknown),
    and besides it: every invocation of the MCP runner recorded a successful read and all reads are one document (the
    identity recomputed, _document_sha; a recorded sha256 the record does not reproduce is a problem, by this
    suite_common's definition when the run's snapshot's is this one, else by the snapshot's own; a note in |notes|
    when neither applies: no snapshot suite_common, or one changed since the snapshot or that does not load); the
    runner's own recorded mismatches (run.json, a line's request_evidence_mismatches); the live
    service_settings.advanced_enabled equals the policy's advanced_path where its exclusions assume it; a
    request_evidence that is not {case: {path: live | table | unknown}}. The per-line evidence not read live is one
    problem per reference path and source, naming its cases and executions."""
    if "mcp" not in (policy.get("arms") or []):
        return []
    out = []
    table_path, _ = cmp.table_for(run_dir, None, reevaluated)
    try:
        cases = cmp.load_by_design(table_path)[0] if table_path else {}
    except cmp.TableError as e:   # the comparison's refusal is the agreement problem; nothing to read the evidence by
        cases = {}
        out.append(_problem("configuration", "no by-design table to read the configuration evidence by: %s" % e,
                            arm="mcp"))
    try:
        defined = cmp.gate_configuration_problems(run_dir, cases, reevaluated)
    except (AttributeError, TypeError, ValueError) as e:
        # a malformed record compare's reading cannot take (named below line by line)
        defined = []
        out.append(_problem("configuration", "compare.gate_configuration_problems could not read the run's "
                                             "configuration records (%s: %s)" % (type(e).__name__, e), arm="mcp"))
    head = ("configuration evidence missing", "configuration mismatch")
    out.extend(_problem("configuration", x, arm="mcp") for x in defined if x.startswith(head))
    line_entries = [x for x in defined if not x.startswith(head)]
    run = sc.read_json(os.path.join(run_dir, "arms", "mcp", "run.json")) or {}
    reads = []
    if isinstance(run.get("service_capabilities"), dict):
        reads.append(("run.json", run["service_capabilities"]))
    invocations = [inv for inv in run.get("invocations") or [] if isinstance(inv, dict)]
    unread = []
    for i, inv in enumerate(invocations):
        rec = inv.get("service_capabilities")
        if not isinstance(rec, dict) or not isinstance(rec.get("document"), dict) or rec.get("error"):
            unread.append("invocation %d started %s%s" % (i, inv.get("started_utc"), (": %s" % str(
                rec.get("error"))[:120]) if isinstance(rec, dict) and rec.get("error") else ""))
            continue
        reads.append(("invocation %d" % i, rec))
    if unread:
        # one problem for the run: a runner of before the record leaves every invocation without one
        out.append(_problem("configuration", "%d of %d invocation(s) of the MCP runner recorded no read of the live "
                                             "capabilities document (%s%s)" % (
                                                 len(unread), len(invocations), "; ".join(unread[:5]),
                                                 "; ..." if len(unread) > 5 else ""), arm="mcp"))
    # a recorded sha256 its record does not reproduce by the definition the runner computed it with: this one where
    # the run's snapshot's suite_common is the pipeline's, else the snapshot's own; a note where that cannot be read
    snap = sc.snapshot_state(run_dir)[0]
    snap_common = ((snap or {}).get("files") or {}).get("suite_common.py")
    one_definition = bool(snap_common) and snap_common == sc.file_sha256(sc.__file__)
    run_definition = None if one_definition else _snapshot_document_sha(run_dir, snap_common)
    for where, rec in reads:
        got = _document_sha(rec)
        if run_definition is not None:
            try:
                want = run_definition(rec)
            except Exception as e:   # the snapshot's definition cannot read this record
                want = "%s: %s" % (type(e).__name__, e)
            if rec.get("sha256") and rec["sha256"] != want:
                out.append(_problem("configuration", "%s: the recorded capabilities sha256 %s is not its record's %s "
                                                     "(the run's snapshot's suite_common.service_document_sha256)" % (
                                                         where, str(rec["sha256"])[:12], str(want)[:12]), arm="mcp"))
            continue
        if rec.get("sha256") and got and rec["sha256"] != got:
            text = "%s: the recorded capabilities sha256 %s is not its record's %s (suite_common.service_document_" \
                   "sha256)" % (where, str(rec["sha256"])[:12], got[:12])
            if one_definition:
                out.append(_problem("configuration", text, arm="mcp"))
            elif notes is not None:
                notes.append(text + "; recorded with another suite_common than this one: compared recomputed")
    shas = sorted(set(_document_sha(rec) for _, rec in reads if isinstance(rec.get("document"), dict)))
    if len(shas) > 1:
        out.append(_problem("configuration", "%d different capabilities documents in one run (sha256 %s): the run "
                                             "spans more than one deployment" % (
                                                 len(shas), ", ".join(x[:12] for x in shas)), arm="mcp"))
    for where, rec in reads:
        for m in rec.get("mismatches") or []:
            out.append(_problem("configuration", "%s: live / table mismatch: %s" % (
                where, json.dumps(m, sort_keys=True)[:300] if not isinstance(m, str) else m[:300]), arm="mcp"))
    doc = (reads[0][1].get("document") if reads else None)
    if isinstance(doc, dict) and isinstance(doc.get("service_settings"), dict):
        premise = any("advanced_path" in (e.get("rules") or ()) for e in policy.get("exclusions") or [])
        want = bool(policy.get("advanced_path"))
        if premise or want:
            settings = doc["service_settings"]
            if "advanced_enabled" not in settings:
                out.append(_problem("configuration", "the live document does not state "
                                                     "service_settings.advanced_enabled: the policy's advanced_path "
                                                     "%s cannot be checked" % ("exclusions" if premise else "setting"),
                                    arm="mcp"))
            elif settings.get("advanced_enabled") is not want:
                out.append(_problem("configuration", "the live service_settings.advanced_enabled is %s, the policy "
                                                     "assumes %s (%s)" % (
                                                         json.dumps(settings.get("advanced_enabled")),
                                                         json.dumps(want), "its advanced_path exclusions"
                                                         if premise else "advanced_path"), arm="mcp"))
    weak, pairs_weak = {}, 0
    for r in rows:
        if r.get("arm") != "mcp":
            continue
        key = "%s r%s" % (r.get("scenario"), r.get("repeat"))
        for m in r.get("request_evidence_mismatches") or []:
            out.append(_problem("configuration", "%s recorded a live / table mismatch: %s" % (
                key, json.dumps(m, sort_keys=True)[:300]), r.get("scenario"), "mcp", r.get("repeat")))
        if "request_evidence" not in r or r.get("request_evidence") is None:
            continue
        ev = r["request_evidence"]
        bad = not isinstance(ev, dict) or any(
            not isinstance(paths, dict) or any(v not in EVIDENCE_SOURCES for v in paths.values())
            for paths in ev.values())
        if bad:
            out.append(_problem("configuration", "malformed request_evidence %s (not {case: {path: live | table | "
                                                 "unknown}})" % json.dumps(ev, sort_keys=True)[:200],
                                r.get("scenario"), "mcp", r.get("repeat")))
            continue
        for case, paths in sorted(ev.items()):
            pairs_weak += any(v != "live" for v in paths.values())
            for path, source in sorted(paths.items()):
                if source != "live":
                    weak.setdefault((path, source), {}).setdefault(case, []).append(key)
    for (path, source), by_case in sorted(weak.items()):
        n = len(set(k for keys in by_case.values() for k in keys))
        what = "was read from the table, not the live document," if source == "table" else \
            "was not found (unknown), live or in the table,"
        out.append(_problem("configuration", "the reference %s %s on %d MCP execution(s) (%s)" % (
            path, what, n, "; ".join("%s: %s" % (case, ", ".join(keys[:12]) + (" ..." if len(keys) > 12 else ""))
                                      for case, keys in sorted(by_case.items()))), arm="mcp"))
    if len(line_entries) != pairs_weak and not any(p["detail"].startswith("compare.gate_configuration_problems")
                                                   for p in out):
        # one definition: compare's per-line reading and this aggregation name the same (execution, case) pairs
        out.append(_problem("configuration", "compare.gate_configuration_problems names %d per-line evidence item(s), "
                                             "the gate's aggregation %d: %s" % (
                                                 len(line_entries), pairs_weak, "; ".join(line_entries[:5])),
                            arm="mcp"))
    return out


# the deployed commit the gate may be given (--expect-service-commit): a git sha, abbreviated to 7 hex digits or more
COMMIT_RE = re.compile(r"\A[0-9a-fA-F]{7,40}\Z")


def commit_matches(recorded, expected):
    """Whether a recorded service commit is the expected one: equal ignoring case, or one an abbreviation (a prefix of
    at least 7 hex digits) of the other."""
    a, b = str(recorded or "").strip().lower(), str(expected or "").strip().lower()
    if len(a) < 7 or len(b) < 7:
        return False
    return a.startswith(b) or b.startswith(a)


def service_identity_problems(identity, expected=None, policy_arms=()):
    """([problem detail], [note]) of the deployed service's identity the MCP arm recorded (|identity|:
    merge_run.service_identity, None without an MCP arm) against the deployed commit |expected| (None: not given).
    Given a commit: a run without an MCP arm (the policy's arms naming it), a run recording no identity at all (a
    runner before the fields), an invocation without a commit (identity not confirmed) and a commit other than the
    expected one each block. In any case the commit changing during an invocation (identity_changed true) or between
    invocations (more than one commit) blocks. identity_changed null (one of the two reads failed; the commit read
    stands or blocks as above), invocations without a commit beside others with one (no commit given), a recorded
    identity not compared (no commit given) and a run recording none (no commit given) are notes."""
    problems, notes = [], []
    want = (" the deployed commit %s (--expect-service-commit)" % expected) if expected else ""
    if identity is None:
        if expected and "mcp" in policy_arms:
            problems.append("no MCP arm recorded: the service identity cannot be confirmed against%s" % want)
        elif expected:
            notes.append("the policy's arms include no MCP arm: the service identity is not judged")
        return problems, notes
    src = identity.get("source") or "arms/mcp"
    n = identity.get("invocations") or 0
    if not identity.get("recorded"):
        if expected:
            problems.append("the MCP arm records no service identity (no service_commit or identity_changed in %s: a "
                            "runner before the identity fields): the gate cannot confirm what it judged against%s"
                            % (src, want))
        else:
            notes.append("the MCP arm records no service identity (a runner before the identity fields): not checked")
        return problems, notes
    commits = identity.get("service_commits") or {}
    errors = identity.get("identity_errors") or []
    err_text = (" (identity_error: %s)" % "; ".join(errors)) if errors else ""
    if identity.get("without_commit"):
        text = "%d of %d MCP invocation(s) record no service_commit (%s)%s" % (
            identity["without_commit"], n, src, err_text)
        if expected:
            problems.append("%s: the identity is not confirmed against%s" % (text, want))
        else:
            notes.append(text)
    if expected:
        for commit, count in commits.items():
            if not commit_matches(commit, expected):
                problems.append("the MCP arm recorded the service commit %s (%d invocation(s)), not%s" % (
                    commit, count, want))
    elif len(commits) > 1:
        problems.append("the MCP invocations recorded %d service commits (%s): the run met more than one deployment"
                        % (len(commits), ", ".join("%s on %d" % kv for kv in commits.items())))
    if identity.get("identity_changed") is True:
        problems.append("the service's commit changed during %d MCP invocation(s) (identity_changed true, %s): the run "
                        "met more than one deployment" % (max(1, identity.get("changed") or 0), src))
    if identity.get("unconfirmed"):
        notes.append("%d MCP invocation(s) record identity_changed null (one of the two reads failed%s): the commit "
                     "read before each %s" % (identity["unconfirmed"],
                                              ("; identity_error: %s" % "; ".join(errors)) if errors else "",
                                              "is judged above" if expected else "was not compared"))
    if not problems:
        if expected:
            notes.append("the MCP arm recorded the deployed service commit %s on its %d invocation(s) (%s)" % (
                expected, n, src))
        elif commits:
            notes.append("the MCP arm recorded the service commit %s (%s); no deployed commit given "
                         "(--expect-service-commit): not compared" % (
                             ", ".join("%s on %d invocation(s)" % kv for kv in commits.items()), src))
    return problems, notes


def library_identity_problems(identity, expected=None, policy_arms=(), pin=None, graphlets=False):
    """([problem detail], [note]) of the deployed graphlet library's identity (|identity|: merge_run.library_identity,
    None without an MCP arm) against the deployed commit |expected| (--expect-library-commit; None: not given) and
    the manifest's pin |pin| (the policy's graphlets.library_commit), as service_identity_problems judges the
    service's: given a commit, a run without an MCP arm (the policy's arms naming it), a run recording no
    library_commit, an invocation without one and another commit each block, and so does a pin other than the given
    commit; with a graphlets block (|graphlets|), every recorded commit is the pin, and recording none blocks; more
    than one commit blocks in any case."""
    problems, notes = [], []
    want = (" the deployed library commit %s (--expect-library-commit)" % expected) if expected else ""
    pins = [pin] if isinstance(pin, str) else list(pin or [])
    if expected:
        for p_ in pins:
            if not commit_matches(p_, expected):
                problems.append("the manifest pins the library commit %s, not%s" % (p_, want))
    if identity is None:
        if (expected or graphlets) and "mcp" in policy_arms:
            problems.append("no MCP arm recorded: the library identity cannot be confirmed")
        return problems, notes
    src = identity.get("source") or "arms/mcp"
    n = identity.get("invocations") or 0
    commits = identity.get("library_commits") or {}
    if not identity.get("recorded"):
        if expected or graphlets:
            problems.append("the MCP arm records no library_commit (%s): the gate cannot confirm the graphlet library "
                            "it judged%s" % (src, (" against" + want) if want else ""))
        return problems, notes
    if identity.get("without_commit") and (expected or graphlets):
        problems.append("%d of %d MCP invocation(s) record no library_commit (%s)" % (
            identity["without_commit"], n, src))
    if len(commits) > 1:
        problems.append("the MCP invocations recorded %d library commits (%s): the run met more than one library" % (
            len(commits), ", ".join("%s on %d" % kv for kv in commits.items())))
    for commit, count in commits.items():
        if expected and not commit_matches(commit, expected):
            problems.append("the MCP arm recorded the library commit %s (%d invocation(s)), not%s" % (
                commit, count, want))
        elif graphlets and pins and not any(commit_matches(commit, p_) for p_ in pins):
            problems.append("the MCP arm recorded the library commit %s (%d invocation(s)), not the manifest's pin "
                            "%s" % (commit, count, ", ".join(pins)))
    if not problems:
        notes.append("the MCP arm recorded the library commit %s on its %d invocation(s) (%s)%s" % (
            ", ".join(commits) or "-", n, src, "" if expected else "; no deployed library commit given "
            "(--expect-library-commit)"))
    return problems, notes


def graphlet_problems(run_dir, policy, manifest, notes, reuse=False):
    """([problem], counts) of area G (module docstring) for a policy with a graphlets block; |reuse|: --reuse."""
    import compare_graphlets as cg
    block = policy.get("graphlets") or {}
    out = []
    try:
        comp = cg.load_or_compute(run_dir, manifest, reuse)
    except Exception as e:   # never take the gate down: an uncomputable comparison blocks
        return [_problem("graphlet", "the graphlet comparison could not be computed (%s: %s)" % (
            type(e).__name__, str(e)[:200]))], {}
    want = policy.get("manifest") or {}
    if (comp.get("manifest") or {}).get("content_sha256") != want.get("content_sha256"):
        out.append(_problem("identity", "graphlet_compare.json compared the manifest content %s, the policy's is %s" % (
            str((comp.get("manifest") or {}).get("content_sha256"))[:16], str(want.get("content_sha256"))[:16])))
    if block.get("whitelist_sha256") and comp.get("whitelist_sha256") != block["whitelist_sha256"]:
        notes.append("the library's answer keys (compare_graphlets.ANSWER_KEYS) differ from the table the policy was "
                     "generated with (%s vs %s)" % (str(comp.get("whitelist_sha256"))[:12],
                                                   str(block["whitelist_sha256"])[:12]))
    for r in comp.get("run_level") or []:
        # lane O's run-level record: W1 lines no op can be found for, so area G cannot be whole (blocks)
        out.append(_problem("graphlet", "%s's run-level record: %s: %s" % (r.get("way"), r.get("declared"),
                                                                          str(r.get("note"))[:240])))
    aside = dict(((x.get("scenario"), x.get("op")), x) for x in block.get("not_required") or [])
    counts = dict((k, 0) for k in ("required", "judged", "agree", "service_by_design", "not_required"))
    scen = comp.get("scenarios") or {}
    for item in block.get("required") or []:
        sid, oid = item.get("scenario"), item.get("op")
        counts["required"] += 1
        entry = ((scen.get(sid) or {}).get("ops") or {}).get(oid)
        where = "op %s (%s)" % (oid, item.get("tool"))
        if entry is None:
            out.append(_problem("missing", "%s: no entry in the graphlet comparison" % where, sid))
            continue
        counts["judged"] += 1
        if (sid, oid) in aside:
            counts["not_required"] += 1
            notes.append("area G: %s %s %s, set aside by the declared file: %s" % (
                sid, where, entry.get("outcome"), aside[(sid, oid)].get("reason")))
            continue
        for w in item.get("ways") or []:
            if (entry.get("ways") or {}).get(w) == "missing":
                out.append(_problem("missing", "%s: no record of %s (%s)" % (where, w, gp_way_name(w)), sid))
        ne = item.get("not_exposed") or {}
        if "not_exposed" in item and not ne.get("W3"):
            out.append(_problem("graphlet", "%s: W3 left out without the owner's reason" % where, sid))
        o = entry.get("outcome")
        if o == "agree":
            counts["agree"] += 1
        elif o == "service_by_design" and all(r.get("verified") for r in entry.get("by_design") or [{}]):
            counts["service_by_design"] += 1
        elif o == "disagree":
            out.append(_problem("graphlet", "%s disagrees (%s): %s" % (
                where, ", ".join(entry.get("disagreeing") or []), cg.disagreement_text(entry)[:300]), sid))
        else:
            reason = str(entry.get("reason") or "")
            if not (o == "inconclusive" and reason.startswith("no record of") and ";" not in reason):
                # (an op inconclusive only for missing records is named by those missing problems)
                out.append(_problem("graphlet", "%s %s%s" % (where, o, (": %s" % reason[:240]) if reason else ""),
                                    sid))
    for sid in block.get("scenarios") or []:
        pre = (scen.get(sid) or {}).get("precheck")
        if not isinstance(pre, dict):
            out.append(_problem("missing", "no pre-check record", sid))
            continue
        for name in block.get("prechecks") or []:
            c = pre.get(name) or {}
            if c.get("outcome") == "pass" or (name == "library_commit" and c.get("outcome") == "not_exposed"):
                continue
            out.append(_problem("graphlet", "pre-check %s %s%s: %s" % (
                name, c.get("outcome"), " (a determinism finding)" if name in ("body", "envelope", "twin")
                and c.get("outcome") == "fail" else "", str(c.get("detail"))[:240]), sid))
    import graphlet_service_checks as gs
    aside = dict((x.get("check"), x) for x in block.get("service_checks_not_required") or [])
    for name, why in gs.blocking_service_checks(comp.get("service_checks"), block.get("service_checks") or [],
                                                aside):
        out.append(_problem("graphlet", "service check %s %s" % (name, why[:280])))
    for c in comp.get("service_checks") or []:
        if c.get("check") in aside and c.get("outcome") != "pass":
            notes.append("area G: service check %s %s, set aside by the declared file: %s" % (
                c["check"], c.get("outcome"), aside[c["check"]].get("reason")))
    checks = dict((c.get("check"), c) for c in comp.get("service_checks") or [])
    for name in block.get("service_checks") or []:
        if name not in checks:
            out.append(_problem("graphlet", "service check %s: not evaluated" % name))
    d = comp.get("drift") or {}
    if d.get("calls"):
        notes.append("area G drift lane (report-only): %d of %d compared answer(s) of the branch head's library "
                     "differ from the deployed commit's" % (len(d.get("differ") or []), d.get("compared") or 0))
    notes.append("area G: %d op call(s) required, %d agree, %d by design (verified), %d set aside by the declared file"
                 % (counts["required"], counts["agree"], counts["service_by_design"], counts["not_required"]))
    return out, counts


def gp_way_name(way):
    import graphlet_precheck as gp
    return gp.WAY_NAMES.get(way, way)


def judge(run_dir, policy, reevaluated=False, expect_service_commit=None, expect_library_commit=None, reuse=False):
    """The gate's verdict on a merged (and compared) run directory: {passed, problems, notes, counts}; with
    |expect_service_commit| (the deployed commit), the service identity the MCP arm recorded is held to it
    (service_identity_problems); with |expect_library_commit| or a policy with a graphlets block, the deployed
    library's (library_identity_problems), and with the block area G is judged (graphlet_problems)."""
    problems, notes = [], []
    name_of = lambda n: sc.reevaluated_name(n, reevaluated)  # noqa: E731
    manifest = sc.read_json(os.path.join(run_dir, "manifest.json"))
    executions_sha256 = sc.file_sha256(os.path.join(run_dir, name_of("executions.jsonl")))
    rows, bad = sc.read_jsonl(os.path.join(run_dir, name_of("executions.jsonl")))
    if not isinstance(manifest, dict):
        problems.append(_problem("identity", "the run has no manifest.json"))
        manifest = {}
    if not rows:
        problems.append(_problem("missing", "the run has no merged %s (merge_run.py)" % name_of("executions.jsonl")))
    if bad:
        problems.append(_problem("missing", "%d unparsable line(s) in %s" % (bad, name_of("executions.jsonl"))))
    scen = dict((s.get("id"), s) for s in manifest.get("scenarios", []) if isinstance(s, dict))
    want = policy["manifest"]

    # identity
    ident = sc.manifest_identity(manifest)
    for key in ("content_sha256", "relations_sha256"):
        if ident[key] != want.get(key):
            problems.append(_problem("identity", "the run manifest's %s %s is not the policy's %s" % (
                key, ident[key][:16], str(want.get(key))[:16])))
    for p in sc.manifest_problems(manifest) if manifest else []:
        problems.append(_problem("identity", "the run manifest: %s" % p))
    arm_runs = dict((name, sc.read_json(os.path.join(path, "run.json"), {}) or {})
                    for name, path in merge_run.arm_dirs(run_dir))
    for d in merge_run.manifest_identity(run_dir, arm_runs)["differs"]:
        problems.append(_problem("identity", d))
    for d in merge_run.release_identity(arm_runs)["changed"]:
        problems.append(_problem("identity", d))
    if not reevaluated:   # a re-evaluated set is judged with one checks module, the re-evaluation's
        for d in merge_run.checks_identity(run_dir, arm_runs)["errors"]:
            problems.append(_problem("identity", d))
    planned, plan_source, plans = run_plan(run_dir, manifest)
    for i, p in enumerate(plans):
        if p.get("manifest") != {"content_sha256": want.get("content_sha256"),
                                 "relations_sha256": want.get("relations_sha256")}:
            problems.append(_problem("identity", "plan %d of plan.json names another manifest (%s)" % (
                i, json.dumps(p.get("manifest"), sort_keys=True))))
    arms = list(policy.get("arms") or [])
    problems.extend(host_identity_problems(policy, arm_runs))
    service = merge_run.service_identity(run_dir, arm_runs)
    sp, sn = service_identity_problems(service, expect_service_commit, arms)
    problems.extend(_problem("identity", x, arm="mcp") for x in sp)
    notes.extend(sn)
    graphlets = bool(policy.get("graphlets"))
    library = None
    if graphlets or expect_library_commit:
        library = merge_run.library_identity(run_dir, arm_runs)
        lp, ln = library_identity_problems(library, expect_library_commit, arms,
                                           (policy.get("graphlets") or {}).get("library_commit"), graphlets)
        problems.extend(_problem("identity", x, arm="mcp") for x in lp)
        notes.extend(ln)
    selected = set(policy["selection"].get("scenarios") or [])
    problems.extend(configuration_problems(run_dir, rows, policy, reevaluated, notes))

    # completeness
    recorded = dict((sc.execution_key(r), r) for r in rows)
    need = required_keys(policy)
    if plan_source is None:
        notes.append("no plan recorded (no plan.json, no orchestrator.json): the policy's executions are the plan")
    else:
        notes.append("plan: %d execution(s) from %s" % (len(planned), plan_source))
        outside = sorted(k for k in need if k not in planned)
        if outside:
            problems.append(_problem("missing", "the run did not plan %d execution(s) the policy requires (e.g. %s)"
                                     % (len(outside), ", ".join("%s %s r%s" % k for k in outside[:5]))))
    pos = dict((sid, i) for i, sid in enumerate(policy["selection"].get("scenarios") or []))
    for k in sorted((planned | need) - set(recorded),
                    key=lambda k: (pos.get(k[0], len(pos)), sc.ARMS.index(k[1]) if k[1] in sc.ARMS else 9, k[2])):
        problems.append(_problem("missing", "planned, not recorded" if k in planned else "required, not recorded",
                                 k[0], k[1], k[2]))

    # outcomes and properties of the policy's arms over the selection
    comp = None
    try:
        comp = cmp.load_or_compute(run_dir, None, reevaluated)
    except (SystemExit, cmp.TableError) as e:
        problems.append(_problem("agreement", "compare.py refused: %s" % e))
    scomp = (comp or {}).get("scenarios") or {}
    pinned = (policy.get("by_design_table") or {}).get("sha256")
    used = ((comp or {}).get("by_design_table") or {}).get("sha256")
    if pinned and used != pinned:
        problems.append(_problem("identity", "the comparison read the by-design table sha256 %s, the policy was derived "
                                             "from %s (regenerate the policy with gate_policy.py, or compare with "
                                             "that table)" % (str(used)[:16], pinned[:16])))
    if comp is not None:
        problems.extend(_problem("identity", x) for x in compared_with_checks_problems(run_dir, comp, reevaluated))
        if comp.get("executions_sha256") != executions_sha256:
            # the lines the gate judges are not the ones the comparison read (edited after it): the comparison's
            # classes and explanations describe other lines
            problems.append(_problem("identity", "%s read %s with sha256 %s; the gate judges it with sha256 %s: "
                                                 "recompute the comparison" % (
                                                     name_of("compare.json"), name_of("executions.jsonl"),
                                                     str(comp.get("executions_sha256"))[:12]
                                                     if comp.get("executions_sha256") else "not stated",
                                                     str(executions_sha256)[:12])))
    # the property reading with the checks that judged the set: a recorded run's snapshot, a re-evaluation's checks;
    # a snapshot file that is not SNAPSHOT.json's is never executed (the reading falls back to GENERIC_CHECKS)
    judging = cmp.judging_checks(run_dir, reevaluated)
    for refused in sc.snapshot_refusal(judging["dir"]):
        problems.append(_problem("identity", "the property reading refuses the checks that judged the set (%s): %s"
                                 % (judging["source"], refused)))
    is_property = _is_property_fn(judging["dir"])
    explanations = table_explanations(comp)
    passed = judged = 0
    passed_by_arm = {}
    by_design_mcp_only = 0
    for (sid, arm, rep), r in sorted(recorded.items(), key=lambda kv: (pos.get(kv[0][0], len(pos)),
                                                                        sc.ARMS.index(kv[0][1])
                                                                        if kv[0][1] in sc.ARMS else 9, kv[0][2])):
        if sid not in selected or arm not in arms or excluded(policy, sid, arm, rep):
            continue
        judged += 1
        s = scen.get(sid) or {}
        o = r.get("outcome")
        if o == "pass":
            passed += 1
            passed_by_arm[arm] = passed_by_arm.get(arm, 0) + 1
        elif o in OUTCOME_BLOCKS:
            why = r.get("inconclusive_reason") or r.get("service_by_design_case") or ""
            if not sc.outcome_known(sc.recorded_outcome(r)):
                # the merge read an outcome outside the vocabulary as error: the problem names the recorded one
                why = "outcome %r outside the vocabulary" % (sc.recorded_outcome(r),)
            problems.append(_problem("outcome", o + (" (%s)" % why if why else ""), sid, arm, rep))
            continue
        elif o == "service_by_design":
            case = r.get("service_by_design_case")
            if case not in cmp.case_ids(s):
                problems.append(_problem("outcome", "service_by_design by case %s, which the scenario does not name "
                                                    "(the manifest names %s)" % (case, cmp.case_ids(s)), sid, arm, rep))
                continue
        else:
            problems.append(_problem("outcome", "outcome %r" % (sc.outcome_value(o),), sid, arm, rep))
            continue
        if o == "pass" or r.get("explained_checks") or r.get("explained_items") or \
                (arm == "mcp" and (sid, rep) in explanations):
            # the line's outcome is not trusted: every check the scenario expects has a result here (determinism
            # not on the first repeat), the status passed, no generic check failed or could not test, and every
            # property passed unless the policy does not require it or the named answered case explains it
            aside = sc.checks_marked_aside(s, arm, r.get("checks"))
            whole = sc.every_property_marked(s, arm, r.get("checks"))
            if whole:
                # the marks leave nothing of the scenario's claim tested via MCP: no pass, whatever the line says
                problems.append(_problem("property", "%s on a passing execution" % sc.scenario_level_mark_detail(
                    whole), sid, arm, rep))
            for name, k, out in property_results(s, r, is_property, generic=True):
                prop = is_property(name)
                label = "%s%s" % (name, (" #%d" % (k + 1)) if k else "")
                if (name, k) in aside:
                    continue   # the manifest marks it not shown via MCP (not implemented / by design): no claim
                if out is None:
                    if name in LATER_REPEATS_ONLY and rep == 0:
                        continue
                    if not (prop and not_required(policy, sid, arm, name)):
                        problems.append(_problem("property" if prop else "outcome", "the expected check %s has no "
                                                 "result on a passing execution" % label, sid, arm, rep))
                    continue
                if out == "pass":
                    continue
                if out in sc.EXPLAINABLE_OUTCOMES and _explained_property(r, name, k, explanations) and \
                        explanations[(sid, rep)][1] in cmp.case_ids(s):
                    # the named answered / clamped case explains it by the table (a failure, or a field its path
                    # removed: suite_common.item_explained): judged with the pair below
                    continue
                if prop:
                    if not not_required(policy, sid, arm, name):
                        problems.append(_problem("property", "required property %s %s" % (label, out), sid, arm, rep))
                elif name == "status" or out != "not_exposed":
                    # a generic check the arm does not expose (the MCP's budget_respected) is no claim; anything
                    # else that did not pass is
                    problems.append(_problem("outcome", "generic check %s %s on a passing execution" % (label, out),
                                             sid, arm, rep))
    for sid in policy["selection"].get("scenarios") or []:
        entry = scomp.get(sid) or {}
        ma = entry.get("mcp_agreement") or {}
        if "mcp" in arms and not excluded(policy, sid, "mcp") and ma:
            # a line whose recorded explanation is not the table's (compare's notes): the gate reads the table's,
            # and a runner that explains otherwise is named
            for key, what in (("explained_items_recorded", "explained items"),
                              ("explained_checks_recorded", "explained checks")):
                for rep, rec in sorted((ma.get(key) or {}).items(), key=lambda kv: str(kv[0])):
                    try:
                        rep_i = int(rep)
                    except (TypeError, ValueError):
                        rep_i = None
                    table = explanations.get((sid, rep_i), ([], None))[0]
                    if key == "explained_checks_recorded":
                        line = recorded.get((sid, "mcp", rep_i)) or {}
                        table = list(dict.fromkeys(str((line.get("checks") or [{}] * (i + 1))[i].get("check"))
                                                   for i in table))
                    problems.append(_problem("agreement", "the MCP line records the %s %s, the table's are %s (a "
                                                          "runner's explanation is not trusted)"
                                             % (what, json.dumps(rec), json.dumps(table)), sid, "mcp", rep_i))
            direct_ran = any((k[0] == sid and k[1] in sc.DIRECT_ARMS) for k in recorded)
            mcp_sbd = any(r.get("outcome") == "service_by_design" or r.get("explained_checks")
                          or r.get("explained_items") for (s_, a_, _), r in recorded.items() if s_ == sid and a_ == "mcp")
            cls = ma.get("class")
            if cls == "disagree":
                # compare's exact-count rule names the odd arm, the MCP's side only when it is the odd one
                what = "MCP disagrees" if ma.get("by") != "exact_counts" else \
                    "the MCP's and a direct arm's exact counts differ"
                problems.append(_problem("disagreement", "%s (%s)" % (what, _describe_agreement(ma)), sid, "mcp"))
            elif mcp_sbd or cls == "service_by_design":
                if direct_ran and ma.get("verified") is not True:
                    problems.append(_problem("outcome", "by-design pair not verified (%s)" % _describe_agreement(ma),
                                             sid, "mcp"))
                elif not direct_ran and ma.get("mcp_verified") is not True:
                    problems.append(_problem("outcome", "by-design case not verified on the MCP side (%s)"
                                             % _describe_agreement(ma), sid, "mcp"))
                elif not direct_ran:
                    by_design_mcp_only += 1
            if policy.get("require_agreement"):
                if cls not in ("agree", "service_by_design") or (cls == "service_by_design"
                                                                 and ma.get("verified") is not True):
                    if cls != "disagree":
                        problems.append(_problem("agreement", "MCP class %s: not compared with a direct arm (%s)" % (
                            cls, _describe_agreement(ma)), sid, "mcp"))
                elif cls == "agree":
                    for name in uncompared_properties(sid, recorded, ma.get("reference"), scen.get(sid), is_property):
                        if not not_required(policy, sid, "mcp", name):
                            problems.append(_problem("agreement", "required property %s not compared with %s" % (
                                name, ma.get("reference")), sid, "mcp"))
        da = entry.get("direct_agreement") or {}
        if all(a in arms for a in sc.DIRECT_ARMS) and da:
            if da.get("class") == "disagree":
                problems.append(_problem("disagreement", "direct-local and direct-remote disagree: %s" % (
                    ("exact counts differ: %s" % cmp.exact_count_text(da["exact_counts_differ"]))
                    if da.get("exact_counts_differ") else
                    da.get("outcome_differs") or "; ".join("%s %s vs %s" % (d.get("check"), d.get("a"), d.get("b"))
                                                           for d in da.get("differing_checks") or [])[:300]), sid))
            elif policy.get("require_agreement") and da.get("class") not in ("identical", "same_checks"):
                problems.append(_problem("agreement", "direct arms %s" % da.get("class"), sid))
    if by_design_mcp_only:
        notes.append("%d by-design scenario(s) verified on the MCP side only (no direct arm ran them)"
                     % by_design_mcp_only)
    judged_with = (comp or {}).get("judged_with") or {}
    if judged_with.get("live_suite_common_differs"):
        # a note, not a problem: the checks judged with the run's own suite_common, the pipeline is this tree's
        notes.append("suite_common.py of this tree (sha256 %s) differs from the one the checks judged with (%s, "
                     "sha256 %s)" % (str(judged_with.get("pipeline_suite_common_sha256"))[:12],
                                     judged_with.get("source"), str(judged_with.get("suite_common_sha256"))[:12]))

    # relations
    rels = dict((x.get("id"), x) for x in (comp or {}).get("relations") or [])
    for item in policy.get("relations") or []:
        rid = item.get("relation")
        for arm in item.get("arms") or []:
            if arm not in arms:
                continue
            got = ((rels.get(rid) or {}).get("arms") or {}).get(arm) or {}
            if got.get("outcome") != "pass":
                problems.append(_problem("relation", "%s on %s: %s%s" % (
                    rid, arm, got.get("outcome") or "not evaluated",
                    (": %s" % str(got.get("detail"))[:200]) if got.get("detail") else ""), arm=arm))
    for arm in arms:
        if not passed_by_arm.get(arm):
            problems.append(_problem("nothing_passed", "no execution of arm %s passed in the policy's selection" % arm,
                                     arm=arm))
    gcounts = None
    if graphlets:
        gprobs, gcounts = graphlet_problems(run_dir, policy, manifest, notes, reuse)
        problems.extend(gprobs)
    order = dict((k, i) for i, k in enumerate(KINDS))
    problems.sort(key=lambda p: order.get(p["kind"], len(order)))
    verdict = {"passed": not problems, "problems": problems, "notes": notes,
               "counts": {"judged": judged, "passed": passed, "planned": len(planned), "required": len(need),
                          "recorded": len(recorded)},
               "policy": {"name": policy.get("name"), "policy_version": policy.get("policy_version"),
                          "sha256": policy_sha256(policy)},
               "run_dir": os.path.abspath(run_dir), "reevaluated": reevaluated, "judged_utc": sc.iso_utc(time.time()),
               "inputs": inputs_digest(run_dir, reevaluated, graphlets),
               "service_identity": {"expected": expect_service_commit, "recorded": service}}
    if graphlets:
        verdict["counts"]["graphlet_ops"] = gcounts
    if library is not None or expect_library_commit:
        verdict["library_identity"] = {"expected": expect_library_commit, "recorded": library,
                                       "pin": (policy.get("graphlets") or {}).get("library_commit")}
    return verdict


def uncompared_properties(sid, recorded, ref_arm, scenario, is_property):
    """The property checks of an agreeing scenario not compared with the reference arm: on some MCP repeat one side
    evaluated it and the other did not (not exposed, inconclusive, no result), or no repeat compared it at all (each
    MCP repeat against the reference's same repeat, else its first; a check of later repeats only, determinism, is
    compared where both have it). A check the scenario's property-level marks set aside on an MCP repeat
    (suite_common.checks_marked_aside: read not_exposed) is left out there: the mark states why it was not compared.
    A marked check the MCP did not read is uncompared as any other."""
    if not ref_arm:
        return []
    ref = dict((k[2], r) for k, r in recorded.items() if k[0] == sid and k[1] == ref_arm)
    if not ref:
        return []
    ref0 = ref[min(ref)]
    judged = ("pass", "fail", "finding")
    compared, skipped, order = set(), set(), []
    for k, r in sorted(recorded.items(), key=lambda kv: kv[0][2]):
        if k[0] != sid or k[1] != "mcp":
            continue
        theirs = dict(((n, i), o) for n, i, o in property_results(scenario, ref.get(k[2], ref0), is_property))
        aside = sc.checks_marked_aside(scenario, "mcp", r.get("checks"))
        for name, i, o in property_results(scenario, r, is_property):
            if (name, i) in aside:
                continue   # a property-level mark of the manifest states it: neither compared nor skipped
            if (name, i) not in order:
                order.append((name, i))
            t = theirs.get((name, i))
            if o in judged and t in judged:
                compared.add((name, i))
            elif o is not None or t is not None:
                skipped.add((name, i))   # one side evaluated it, the other could not (not exposed, inconclusive)
    out = []
    for key in order:
        if (key in skipped or key not in compared) and key[0] not in out:
            out.append(key[0])
    return out


def _describe_agreement(ma):
    keys = ("case", "by", "side", "verified", "mcp_verified", "direct_verified")
    text = ", ".join("%s %s" % (k, ma.get(k)) for k in keys if ma.get(k) is not None)
    detail = ma.get("detail") or ma.get("direct_detail") or ma.get("mcp_detail")
    return (text + ("; " + str(detail)[:200] if detail else "")) or "-"


HASH_KEYS = (("checks_sha256", "checks.py"), ("suite_common_sha256", "suite_common.py"))


def compared_with_checks_problems(run_dir, comp, reevaluated=False):
    """Why a set's comparison, or the checks the gate reads beside it, are no evidence for the checks that judged the
    set: the comparison does not state that they judged it (judged_with source: the run's snapshot for a recorded set,
    the re-evaluation's checks for a re-evaluated one), or the checks.py / suite_common.py sha256 it states are not the
    ones recorded for them (SNAPSHOT.json; the checks of reevaluated.json): a comparison computed while those files
    were changed, restored before the gate, or with other checks. A re-evaluated set is evidence by its hashes, not by
    its directory: compared with this tree's checks (its re-evaluation's directory gone), it stands when their sha256
    are the recorded ones; and the re-evaluation's directory (the gate's property reading) must still hold them. The
    advice is the step that clears the problem: recompute the comparison where the checks the comparison would read
    now are the recorded ones, re-evaluate the set otherwise. [] for a recorded set without a snapshot."""
    name = sc.reevaluated_name("compare.json", reevaluated)
    if reevaluated:
        mark = (merge_run.reevaluation_mark(run_dir) or {}).get("checks") or {}
        source, whose = "reevaluation", "the re-evaluation's checks"
        where = "the re-evaluation records (reevaluated.json)"
        want = {"checks.py": mark.get("checks_sha256"), "suite_common.py": mark.get("suite_common_sha256")}
    else:
        snap, _, _ = sc.snapshot_state(run_dir)
        if snap is None:
            return []
        files = snap.get("files") or {}
        source, whose, where = "snapshot", "the run's snapshot", "SNAPSHOT.json records"
        want = {"checks.py": files.get("checks.py"), "suite_common.py": files.get("suite_common.py")}
    now = cmp.judging_checks(run_dir, reevaluated)
    out = []
    if reevaluated:
        # the checks a comparison and the gate's property reading read now (the re-evaluation's directory, or this
        # tree's where it is gone) are not the ones that judged the set
        place = "the re-evaluation's checks directory" if now.get("source") == "reevaluation" else \
            "this tree (the re-evaluation's checks directory is gone)"
        out.extend("%s (sha256 %s) of %s (%s) is not the %s %s: re-evaluate the set"
                   % (f, str(now.get(key))[:12], place, now.get("dir"), str(want[f])[:12], where)
                   for key, f in HASH_KEYS if now.get(key) != want[f])
    advice = "recompute the comparison" if not out else "re-evaluate the set"
    jw = (comp or {}).get("judged_with")
    sources = (source, "live") if reevaluated else (source,)
    if not isinstance(jw, dict) or jw.get("source") not in sources:
        return out + ["%s was not judged with %s (judged_with source %s): %s"
                      % (name, whose, jw.get("source") if isinstance(jw, dict) else "not stated", advice)]
    return out + ["%s was judged with %s sha256 %s, not the %s %s: %s"
                  % (name, f, str(jw.get(key))[:12], str(want[f])[:12], where, advice)
                  for key, f in HASH_KEYS if jw.get(key) != want[f]]


GRAPHLET_INPUTS = ("graphlet_inputs",)


def inputs_digest(run_dir, reevaluated=False, graphlets=False):
    """{file: sha256 or None} of what a verdict judged: the merged lines, the comparison, the manifest and the plan
    (a report compares it with the files it reads: a run changed after the verdict, e.g. by a later resume without
    --gate, makes the verdict stale); with |graphlets| also `graphlet_inputs`, the sha256 of every input the graphlet
    comparison reads (compare_graphlets.input_digests)."""
    out = {}
    if graphlets:
        import compare_graphlets as cg
        out["graphlet_inputs"] = hashlib.sha256(sc.canonical_json(cg.input_digests(run_dir))).hexdigest()
    for name in (sc.reevaluated_name("executions.jsonl", reevaluated), sc.reevaluated_name("compare.json", reevaluated),
                 "manifest.json", PLAN_FILE):
        path = os.path.join(run_dir, name)
        try:
            with open(path, "rb") as fh:
                out[name] = hashlib.sha256(fh.read()).hexdigest()
        except OSError:
            out[name] = None
    return out


def stale_inputs(run_dir, verdict, reevaluated=False):
    """The files whose content differs from what |verdict| judged ([] when it is current; every file for a verdict
    that records none)."""
    was = verdict.get("inputs") if isinstance(verdict.get("inputs"), dict) else {}
    now = inputs_digest(run_dir, reevaluated, any(n in was for n in GRAPHLET_INPUTS))
    return [name for name in now if was.get(name, "unrecorded") != now[name]]


def policy_sha256(policy):
    return hashlib.sha256(sc.canonical_json(policy)).hexdigest()


def render(verdict, max_lines=400):
    """The gate's verdict as text: the headline, then the problems grouped by kind, then the notes."""
    probs = verdict["problems"]
    c = verdict["counts"]
    pol = verdict["policy"]
    head = "GATE %s: %s (policy %s v%s sha256 %s; %d execution(s) judged, %d passed; %d planned, %d required)" % (
        "PASS" if verdict["passed"] else "FAIL",
        "every planned and required execution recorded, nothing blocks" if verdict["passed"]
        else "%d problem(s)" % len(probs), pol.get("name"), pol.get("policy_version"), (pol.get("sha256") or "")[:12],
        c["judged"], c["passed"], c["planned"], c["required"])
    out = [head]
    by_kind = {}
    for p in probs:
        by_kind.setdefault(p["kind"], []).append(p)
    shown = 0
    for kind in KINDS:
        items = by_kind.get(kind) or []
        if not items:
            continue
        out.append("  %s: %d" % (kind, len(items)))
        for p in items:
            if shown >= max_lines:
                break
            where = " ".join(x for x in (p.get("scenario"), p.get("arm"),
                                         ("r%s" % p["repeat"]) if p.get("repeat") is not None else None) if x)
            out.append("    %s%s" % ((where + ": ") if where else "", p["detail"]))
            shown += 1
    if shown < len(probs):
        out.append("  ... %d more (all in the gate's JSON)" % (len(probs) - shown))
    for n in verdict["notes"]:
        out.append("  note: %s" % n)
    return "\n".join(out) + "\n"


def gate_file(reevaluated):
    return sc.reevaluated_name("gate.json", reevaluated)


# ---------------------------------------------------------------------------------------------------------------
# command line

def _cmd_plan(argv):
    ap = argparse.ArgumentParser(prog="gate.py plan")
    ap.add_argument("run_dir")
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--arms", required=True)
    ap.add_argument("--areas", default="")
    ap.add_argument("--limit", default="")
    ap.add_argument("--ids", default="", help="comma-separated scenario ids")
    ap.add_argument("--warm-every", default="10")
    ap.add_argument("--mode", default="paired")
    ap.add_argument("--session", type=int, default=None)
    a = ap.parse_args(argv)
    manifest = sc.read_json(a.manifest)
    if not isinstance(manifest, dict):
        sys.stderr.write("gate.py plan: %s is not a manifest\n" % a.manifest)
        return 2
    entry = plan_entry(manifest, _split(a.areas) or None, _int_or(a.limit, None), _split(a.ids) or None,
                       _int_or(a.warm_every, 10), _split(a.arms), a.mode, a.session)
    added = write_plan(a.run_dir, entry)
    sys.stdout.write("plan %s: %d scenario(s), %s%s\n" % (
        "written" if added else "unchanged", entry["counts"]["scenarios"],
        ", ".join("%s %d" % kv for kv in sorted(entry["counts"]["executions"].items())) or "no execution",
        "" if added else " (the run's last plan)"))
    return 0


def _cmd_ids(argv):
    ap = argparse.ArgumentParser(prog="gate.py ids")
    ap.add_argument("file")
    ap.add_argument("--manifest", required=True)
    a = ap.parse_args(argv)
    try:
        ids = read_ids(a.file, sc.read_json(a.manifest))
    except (OSError, ValueError) as e:
        sys.stderr.write("gate.py ids: %s\n" % e)
        return 2
    sys.stdout.write(",".join(ids) + "\n")
    return 0


def _cmd_resume_ids(argv):
    ap = argparse.ArgumentParser(prog="gate.py resume-ids")
    ap.add_argument("run_dir")
    ap.add_argument("--manifest", required=True)
    a = ap.parse_args(argv)
    sys.stdout.write(",".join(resume_ids(a.run_dir, sc.read_json(a.manifest))) + "\n")
    return 0


def _cmd_policy_ids(argv):
    ap = argparse.ArgumentParser(prog="gate.py policy-ids")
    ap.add_argument("policy")
    ap.add_argument("--arms", default="")
    a = ap.parse_args(argv)
    try:
        policy = load_policy(a.policy)
    except (OSError, ValueError) as e:
        sys.stderr.write("gate.py policy-ids: %s\n" % e)
        return 2
    sys.stdout.write(",".join(policy_ids(policy, _split(a.arms) or None)) + "\n")
    return 0


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if argv and argv[0] == "plan":
        return _cmd_plan(argv[1:])
    if argv and argv[0] == "ids":
        return _cmd_ids(argv[1:])
    if argv and argv[0] == "policy-ids":
        return _cmd_policy_ids(argv[1:])
    if argv and argv[0] == "resume-ids":
        return _cmd_resume_ids(argv[1:])
    ap = argparse.ArgumentParser(description="The strict deploy gate of a suite run (exit 0 PASS, 1 FAIL).")
    ap.add_argument("run_dir")
    ap.add_argument("policy")
    ap.add_argument("--reevaluated", action="store_true", help="judge the run's re-evaluated set")
    ap.add_argument("--out", default=None, help="the verdict's JSON (default RUN_DIR/gate.json, "
                                                "gate.reevaluated.json with --reevaluated)")
    ap.add_argument("--max-lines", type=int, default=400)
    ap.add_argument("--expect-service-commit", metavar="SHA", default=None,
                    help="the deployed service commit (run-gate's LIVE line): every MCP invocation must have read it")
    ap.add_argument("--reuse", action="store_true",
                    help="area G: reuse graphlet_compare.json when every input it read is current (default: recompute)")
    ap.add_argument("--expect-library-commit", metavar="SHA", default=None,
                    help="the deployed graphlet library commit (the service's library_commit): every MCP invocation "
                         "must have read it, and the manifest's pin must be it")
    a = ap.parse_args(argv)
    for flag, value in (("--expect-service-commit", a.expect_service_commit),
                        ("--expect-library-commit", a.expect_library_commit)):
        if value is not None and not COMMIT_RE.match(value):
            sys.stderr.write("gate.py: %s %r is not a git commit (7 to 40 hex digits)\n" % (flag, value))
            return 2
    try:
        policy = load_policy(a.policy)
    except (OSError, ValueError) as e:
        sys.stderr.write("gate.py: the policy: %s\n" % e)
        return 2
    if not os.path.isdir(a.run_dir):
        sys.stderr.write("gate.py: %s is not a run directory\n" % a.run_dir)
        return 2
    verdict = judge(a.run_dir, policy, a.reevaluated, a.expect_service_commit, a.expect_library_commit, a.reuse)
    verdict["policy"]["path"] = os.path.abspath(a.policy)
    sc.write_json_atomic(a.out or os.path.join(a.run_dir, gate_file(a.reevaluated)), verdict)
    sys.stdout.write(render(verdict, a.max_lines))
    return 0 if verdict["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
