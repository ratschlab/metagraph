#!/usr/bin/env python3
"""The service's own behaviour in area G, checked positively (DESIGN-v8 §5, D-R6), and the declared per-way
differences verified as by-design pairs (D-R12). Two stages of compare_graphlets.py (registered on import).

Positive service checks, functions over W3's records and the scenario's facts (scenario_facts: from W1's document
by the suite's own reader, a scenario's graphlet.facts stating them where the corpus says so):
  walk_display       every graphlet_walks / graphlet_walk answer of W3 carries walk_display "first_parent";
  label_evidence     a W3 graphlet_compare answer carries label_evidence "lower_bound" exactly when one of its two
                     graphlets has cut label lists (a segment whose recorded labels are fewer than its labels_total);
  qualified_cause    present exactly with label_evidence, and consistent with the library's answer to the same call
                     (W4's, else W1's): it states each side's cap (labels.max_labels_per_node of the side's request)
                     and "no difference within the recorded evidence" exactly when the library counts no difference;
  hint               a W3 compare answer carries hint exactly when one of its graphlets has more than 64 labels on a
                     node (labels_total);
  in_service         a W3 graphlet_summary of a body over 8 MiB carries in_service false and a download; no other
                     summary carries in_service false;
  max_bytes_refused  a W3 call with max_bytes above the service's 262,144 is refused bad_argument, and the library
                     ways answer the same call (the library's ceiling is 16 MiB);
  forbidden          a W3 call the op plan sends with a wrong access_token (op probe "wrong_token") is forbidden;
  bad_handle         a W3 call with a malformed handle (op probe "malformed_handle") is bad_handle;
  library_commit     every W3 call reports the deployed library_commit and the manifest's pin (the pre-checks).
Each check: pass, fail (naming every call that breaks it), or inconclusive (exercised false: no call of the run
exercises it; or a fact it needs is unknown).

By-design pairs: an op naming a case (`by_design`: an id or a list) of the by-design table (compare.table_for: the
run's snapshot table, else service_by_design.json) whose `graphlet` block states the pair (compare.graphlet_case_
problems: tools, args predicates over args_resolved, the library ways' expectation, W3's expectation). The pair holds
when the op's tool is listed, the args predicates hold, every library way the op runs matches `library` and W3
matches `service`; then W3-W4 and W3-W5 may differ and the op is service_by_design (verified), its library pairs
still compared. A pair that does not hold is a disagreement naming the side (table, request, library, service).
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import compare as cmp  # noqa: E402
import compare_graphlets as cg  # noqa: E402
import graphlet_precheck as gp  # noqa: E402
import suite_common as sc  # noqa: E402

SERVICE_CHECKS = ("walk_display", "label_evidence", "qualified_cause", "hint", "in_service", "max_bytes_refused",
                  "forbidden", "bad_handle", "library_commit")
SERVICE_MAX_BYTES = 262144
HINT_LABELS_PER_NODE = 64
IN_SERVICE_BYTES = 8 * 1024 * 1024
NO_DIFFERENCE = "no difference within the recorded evidence"
PROBES = {"wrong_token": ("forbidden", "forbidden"), "malformed_handle": ("bad_handle", "bad_handle")}
# the pairs a by-design case lets differ (the service side against the library and the oracle)
BY_DESIGN_PAIRS = ("W3-W4", "W3-W5")


# ---------------------------------------------------------------------------------------------------------------
# facts of a scenario

def document_facts(text, checks_mod):
    """{cut_lists, max_labels_per_node, body_bytes} of an MGT document by the suite's reader (None values when it
    cannot be read)."""
    out = {"cut_lists": None, "max_labels_per_node": None, "body_bytes": None}
    if text is None:
        return out
    env = gp.envelope_of(text)
    j = len((text.split("\n", 2)[1] + "\n").encode("utf-8")) if env is not None else 0   # the J line
    out["body_bytes"] = len(text.encode("utf-8")) - j
    if checks_mod is None:
        return out
    try:
        arms = checks_mod.graphlet_arms(text)
    except Exception:
        return out
    segs = [s for a in arms.values() for s in a.get("segments") or []]
    if segs and all("labels_total" in s for s in segs):
        out["cut_lists"] = any(s.get("labels_truncated") for s in segs)
        out["max_labels_per_node"] = max(s.get("labels_total") or 0 for s in segs)
    return out


def scenario_facts(run_dir, manifest, checks_mod=None):
    """{scenario: {cut_lists, max_labels_per_node, body_bytes, over_8mib, cap}}: W1's document read by the suite,
    each fact replaced by the scenario's graphlet.facts where it states one (cut_lists, max_labels_per_node,
    body_over_8mib); cap = request.strategy.labels.max_labels_per_node."""
    checks_mod = checks_mod if checks_mod is not None else gp.load_checks(run_dir)
    out = {}
    for scen in gp.graphlet_scenarios(manifest):
        sid = scen.get("id")
        f = document_facts(gp.read_doc(run_dir, "W1", sid), checks_mod)
        stated = (scen.get("graphlet") or {}).get("facts") or {}
        for k in ("cut_lists", "max_labels_per_node"):
            if k in stated:
                f[k] = stated[k]
        f["over_8mib"] = stated["body_over_8mib"] if "body_over_8mib" in stated else (
            None if f["body_bytes"] is None else f["body_bytes"] > IN_SERVICE_BYTES)
        labels = (((scen.get("request") or {}).get("strategy") or {}).get("labels") or {})
        f["cap"] = labels.get("max_labels_per_node")
        out[sid] = f
    return out


def service_fields(rec):
    """W3's service fields of a record: the ones its runner kept, and any key of its answer outside the library's
    table (compare_graphlets.split_service_fields)."""
    out = dict(rec.get("service_fields") or {}) if isinstance(rec.get("service_fields"), dict) else {}
    out.update(cg.split_service_fields(rec.get("tool"), rec.get("answer"))[1])
    return out


def _w3_calls(recs, manifest, tools=None):
    """[(scenario, op dict, W3 record)] of the planned W3 calls (of |tools|) that answered or failed."""
    by = cg.scenario_by_id(manifest)
    out = []
    for (sid, oid), r in sorted(recs.get("W3", {}).items(), key=lambda kv: (str(kv[0][0]), str(kv[0][1]))):
        scen = by.get(sid)
        op = next((o for o in gp.scenario_ops(scen) if o["id"] == oid), None) if scen else None
        if op is None or "W3" not in gp.op_ways(scen, op) or r.get("selector_error"):
            continue
        if tools is None or r.get("tool") in tools:
            out.append((sid, op, r))
    return out


def _result(check, calls, failures, unknown=None, what=""):
    if failures:
        return {"check": check, "outcome": "fail", "exercised": True, "calls": calls, "failures": failures[:50],
                "detail": "%d of %d call(s) break it: %s" % (len(failures), calls, "; ".join(
                    "%s %s: %s" % (f["scenario"], f["op"], f["detail"]) for f in failures[:3]))}
    if not calls:
        return {"check": check, "outcome": "inconclusive", "exercised": False, "calls": 0, "failures": [],
                "detail": "not exercised: no W3 call of the run %s" % what}
    if unknown:
        return {"check": check, "outcome": "inconclusive", "exercised": True, "calls": calls, "failures": [],
                "detail": "%d call(s) without the fact they need: %s" % (len(unknown), "; ".join(unknown[:3]))}
    return {"check": check, "outcome": "pass", "exercised": True, "calls": calls, "failures": [],
            "detail": "%d call(s) %s" % (calls, what)}


def _fail(sid, op, detail):
    return {"scenario": sid, "op": op["id"], "detail": detail}


def _b_scenario(sid, op):
    b = (op.get("args") or {}).get("b")
    return b.get("graphlet") if isinstance(b, dict) and b.get("graphlet") else sid


# ---------------------------------------------------------------------------------------------------------------
# the checks

def check_walk_display(recs, manifest, facts):
    calls, bad = 0, []
    for sid, op, r in _w3_calls(recs, manifest, ("graphlet_walks", "graphlet_walk")):
        if r.get("error"):
            continue
        calls += 1
        got = service_fields(r).get("walk_display")
        if got != "first_parent":
            bad.append(_fail(sid, op, "walk_display %s" % gp.compact(got if "walk_display" in service_fields(r)
                                                                         else gp.MISSING)))
    return _result("walk_display", calls, bad, what="of graphlet_walks / graphlet_walk answered")


def _compare_calls(recs, manifest):
    return [(sid, op, r) for sid, op, r in _w3_calls(recs, manifest, ("graphlet_compare",)) if not r.get("error")]


def check_label_evidence(recs, manifest, facts):
    calls, bad, unknown = 0, [], []
    for sid, op, r in _compare_calls(recs, manifest):
        calls += 1
        sides = [facts.get(x, {}).get("cut_lists") for x in (sid, _b_scenario(sid, op))]
        if any(v is None for v in sides):
            unknown.append("%s %s: cut lists unknown" % (sid, op["id"]))
            continue
        sf = service_fields(r)
        want = any(sides)
        if want and sf.get("label_evidence") != "lower_bound":
            bad.append(_fail(sid, op, "a cut-list graphlet, label_evidence %s" % gp.compact(
                sf.get("label_evidence", gp.MISSING))))
        elif not want and "label_evidence" in sf:
            bad.append(_fail(sid, op, "no cut list, yet label_evidence %s" % gp.compact(sf["label_evidence"])))
    return _result("label_evidence", calls, bad, unknown, "of graphlet_compare answered")


def _library_answer(recs, sid, oid):
    for w in ("W4", "W1"):
        r = recs.get(w, {}).get((sid, oid))
        if r is not None and not r.get("error") and isinstance(r.get("answer"), dict):
            return r["answer"]
    return None


def check_qualified_cause(recs, manifest, facts):
    calls, bad, unknown = 0, [], []
    for sid, op, r in _compare_calls(recs, manifest):
        calls += 1
        sf = service_fields(r)
        cause = sf.get("qualified_cause")
        if ("label_evidence" in sf) != ("qualified_cause" in sf):
            bad.append(_fail(sid, op, "qualified_cause %s with label_evidence %s" % (
                "present" if "qualified_cause" in sf else "absent", "present" if "label_evidence" in sf else
                "absent")))
            continue
        if cause is None:
            continue
        if not isinstance(cause, str):
            bad.append(_fail(sid, op, "qualified_cause %s is not a text" % gp.compact(cause)))
            continue
        lib = _library_answer(recs, sid, op["id"])
        counts = (lib or {}).get("counts")
        if not isinstance(counts, dict):
            unknown.append("%s %s: no library answer with counts" % (sid, op["id"]))
            continue
        none_differ = not any((counts.get(k) or 0) for k in ("only_in_a", "only_in_b", "differ"))
        if none_differ != (NO_DIFFERENCE in cause):
            bad.append(_fail(sid, op, "the library counts %s, qualified_cause %s" % (
                "no difference" if none_differ else gp.compact(counts), gp.compact(cause))))
            continue
        for side in (sid, _b_scenario(sid, op)):
            cap = facts.get(side, {}).get("cap")
            if isinstance(cap, int) and str(cap) not in cause:
                bad.append(_fail(sid, op, "qualified_cause %s does not state the cap %d of %s" % (
                    gp.compact(cause), cap, side)))
                break
    return _result("qualified_cause", calls, bad, unknown, "of graphlet_compare answered")


def check_hint(recs, manifest, facts):
    calls, bad, unknown = 0, [], []
    for sid, op, r in _compare_calls(recs, manifest):
        calls += 1
        sides = [facts.get(x, {}).get("max_labels_per_node") for x in (sid, _b_scenario(sid, op))]
        if any(v is None for v in sides):
            unknown.append("%s %s: labels per node unknown" % (sid, op["id"]))
            continue
        want = max(sides) > HINT_LABELS_PER_NODE
        has = "hint" in service_fields(r)
        if want != has:
            bad.append(_fail(sid, op, "%d labels per node at most, hint %s" % (max(sides),
                                                                              "present" if has else "absent")))
    return _result("hint", calls, bad, unknown, "of graphlet_compare answered")


def check_in_service(recs, manifest, facts):
    calls, bad, over = 0, [], 0
    for sid, op, r in _w3_calls(recs, manifest, ("graphlet_summary",)):
        if r.get("error"):
            continue
        sf = service_fields(r)
        if facts.get(sid, {}).get("over_8mib"):
            over += 1
            calls += 1
            if sf.get("in_service") is not False or not sf.get("download"):
                bad.append(_fail(sid, op, "a body over 8 MiB: in_service %s, download %s" % (
                    gp.compact(sf.get("in_service", gp.MISSING)), "present" if sf.get("download") else "absent")))
        elif sf.get("in_service") is False:
            bad.append(_fail(sid, op, "in_service false on a body within 8 MiB"))
    res = _result("in_service", over, bad, what="of graphlet_summary on a body over 8 MiB")
    return res


def check_max_bytes_refused(recs, manifest, facts):
    calls, bad = 0, []
    for sid, op, r in _w3_calls(recs, manifest):
        mb = (r.get("args_resolved") or {}).get("max_bytes")
        if not isinstance(mb, int) or isinstance(mb, bool) or mb <= SERVICE_MAX_BYTES:
            continue
        calls += 1
        code = (r.get("error") or {}).get("code") if r.get("error") else None
        if code != "bad_argument":
            bad.append(_fail(sid, op, "max_bytes %d on W3: %s" % (mb, "error %s" % code if code else "answered")))
            continue
        refused = [w for w in gp.LIBRARY_WAYS if w in gp.op_ways(cg.scenario_by_id(manifest)[sid], op)
                   and (recs.get(w, {}).get((sid, op["id"])) or {}).get("error")]
        if refused:
            bad.append(_fail(sid, op, "max_bytes %d: the library refused it too on %s" % (mb, ", ".join(refused))))
    return _result("max_bytes_refused", calls, bad, what="with max_bytes above %d" % SERVICE_MAX_BYTES)


def _check_probe(name, probe, code, recs, manifest):
    calls, bad = 0, []
    for sid, op, r in _w3_calls(recs, manifest):
        if op.get("probe") != probe:
            continue
        calls += 1
        got = (r.get("error") or {}).get("code") if r.get("error") else None
        if got != code:
            bad.append(_fail(sid, op, "%s: %s" % (probe, "error %s" % got if got else "answered")))
    return _result(name, calls, bad, what="probing %s" % probe)


def check_library_commit(prechecks):
    calls, bad, unknown = 0, [], []
    for p in prechecks:
        lc = p.get("library_commit") or {}
        if lc.get("outcome") == "not_exposed":
            continue
        n = sum((lc.get("calls") or {}).values())
        calls += n
        if lc.get("outcome") == "fail":
            bad.append({"scenario": p["scenario"], "op": "*", "detail": lc.get("detail")})
        elif lc.get("outcome") != "pass":
            unknown.append("%s: %s" % (p["scenario"], lc.get("detail")))
    return _result("library_commit", calls, bad, unknown, "reporting the pinned and deployed library_commit")


def service_checks(recs, manifest, facts, prechecks):
    return [check_walk_display(recs, manifest, facts), check_label_evidence(recs, manifest, facts),
            check_qualified_cause(recs, manifest, facts), check_hint(recs, manifest, facts),
            check_in_service(recs, manifest, facts), check_max_bytes_refused(recs, manifest, facts),
            _check_probe("forbidden", "wrong_token", "forbidden", recs, manifest),
            _check_probe("bad_handle", "malformed_handle", "bad_handle", recs, manifest),
            check_library_commit(prechecks)]


def required_service_checks(manifest):
    """The service checks a manifest's op plan exercises (the gate requires these): walk_display with a W3 walks /
    walk op, label_evidence, qualified_cause and hint with a W3 compare op, in_service with a scenario stating
    facts.body_over_8mib, max_bytes_refused with a W3 op sending max_bytes above 262,144, forbidden / bad_handle with
    a probe op, library_commit with any W3 op."""
    out = set()
    for scen in gp.graphlet_scenarios(manifest):
        if ((scen.get("graphlet") or {}).get("facts") or {}).get("body_over_8mib"):
            out.add("in_service")
        for op in gp.scenario_ops(scen):
            if "W3" not in gp.op_ways(scen, op):
                continue
            out.add("library_commit")
            tool, args = op.get("tool"), op.get("args") or {}
            if tool in ("graphlet_walks", "graphlet_walk") and op.get("probe") is None:
                out.add("walk_display")
            if tool == "graphlet_compare" and op.get("probe") is None:
                out.update(("label_evidence", "qualified_cause", "hint"))
            mb = args.get("max_bytes")
            if isinstance(mb, int) and not isinstance(mb, bool) and mb > SERVICE_MAX_BYTES:
                out.add("max_bytes_refused")
            if op.get("probe") in PROBES:
                out.add(PROBES[op["probe"]][0])
    return [c for c in SERVICE_CHECKS if c in out]


# ---------------------------------------------------------------------------------------------------------------
# by-design pairs

def _alternatives(code):
    return [x.strip() for x in str(code).split("|") if x.strip()]


def expectation_holds(want, rec, way):
    """(holds, detail) of a way's record against a case side ({"answered": true} or {"error": {code, fields?}}, W3
    also {"service_fields": {...}})."""
    if rec is None:
        return False, "%s recorded nothing" % way
    if rec.get("selector_error"):
        return False, "%s did not run the call (selector)" % way
    err = rec.get("error") if isinstance(rec.get("error"), dict) else None
    if want.get("answered"):
        if err:
            return False, "%s refused it (%s), the case says answered" % (way, err.get("code"))
    else:
        e = want.get("error") or {}
        if not err:
            return False, "%s answered, the case says error %s" % (way, e.get("code"))
        if err.get("code") not in _alternatives(e.get("code")):
            return False, "%s's error %s, the case says %s" % (way, err.get("code"), e.get("code"))
        fields = err.get("fields") if isinstance(err.get("fields"), dict) else {}
        for k, v in sorted((e.get("fields") or {}).items()):
            if not sc.json_equal(fields.get(k), v):
                return False, "%s's error field %s %s, the case says %s" % (way, k, gp.compact(fields.get(k)),
                                                                             gp.compact(v))
    for k, v in sorted((want.get("service_fields") or {}).items()):
        got = service_fields(rec).get(k, gp.MISSING)
        if got is gp.MISSING or not sc.json_equal(got, v):
            return False, "%s's service field %s %s, the case says %s" % (way, k, gp.compact(got), gp.compact(v))
    return True, ""


def judge_by_design(entry, scen, op, recs, cases):
    """The op's by-design pair(s): entry["by_design"] and the op's outcome (module docstring)."""
    names = op.get("by_design")
    names = [names] if isinstance(names, str) else [x for x in names or [] if isinstance(x, str)]
    sid, oid = scen.get("id"), op["id"]
    ways = gp.op_ways(scen, op)
    results = []
    for cid in names:
        case = cases.get(cid)
        g = (case or {}).get("graphlet") if isinstance((case or {}).get("graphlet"), dict) else None
        res = {"case": cid, "verified": False, "side": None, "detail": ""}
        if g is None:
            res.update(side="table", detail="%s is %s" % (cid, "not an area G case (no graphlet block)" if case
                                                          else "not a case of the by-design table"))
        elif op.get("tool") not in (g.get("tools") or []):
            res.update(side="request", detail="%s does not list the tool %s" % (cid, op.get("tool")))
        else:
            ref = recs.get("W3", {}).get((sid, oid)) or recs.get("W1", {}).get((sid, oid)) or {}
            ok, why = cmp.request_holds({"id": cid, "request": g.get("args")},
                                        {"mcp": {"arguments": ref.get("args_resolved") or {}}}, None)
            if not ok:
                res.update(side="request", detail=why)
            else:
                problems = []
                for w in [x for x in gp.LIBRARY_WAYS if x in ways]:
                    held, why = expectation_holds(g.get("library") or {}, recs.get(w, {}).get((sid, oid)), w)
                    if not held:
                        problems.append(("library", why))
                if "W3" not in ways:
                    problems.append(("service", "W3 does not run the op"))
                else:
                    held, why = expectation_holds(g.get("service") or {}, recs.get("W3", {}).get((sid, oid)), "W3")
                    if not held:
                        problems.append(("service", why))
                if problems:
                    res.update(side=problems[0][0], detail="; ".join(x[1] for x in problems))
                else:
                    res.update(verified=True, detail="the library ways and W3 answer as %s states" % cid)
        results.append(res)
    entry["by_design"] = results
    if entry["outcome"] in ("precondition_failed",) or entry.get("smoke"):
        return entry
    rest = [p for p in entry.get("disagreeing") or [] if p not in BY_DESIGN_PAIRS]
    unverified = ["by-design %s: %s" % (r["case"], r["side"]) for r in results if not r["verified"]]
    if rest or unverified:
        entry["outcome"] = "disagree"
        entry["disagreeing"] = rest + unverified
    elif entry["outcome"] == "inconclusive":
        entry.pop("disagreeing", None)
    else:
        entry.pop("disagreeing", None)
        entry["outcome"] = "service_by_design"
    return entry


def by_design_stage(result, run_dir, manifest, recs, prechecks):
    """Stage of compare_graphlets: judges every op that names by-design cases, with the run's table."""
    ops = [(scen, op) for scen in gp.graphlet_scenarios(manifest) for op in gp.scenario_ops(scen)
           if op.get("by_design")]
    if not ops:
        return
    path, warnings = cmp.table_for(run_dir, None)
    result.setdefault("warnings", []).extend(w for w in warnings if w not in result.get("warnings", []))
    try:
        cases, info = cmp.load_by_design(path) if path else ({}, None)
    except cmp.TableError as e:
        cases, info = {}, None
        result["warnings"].append("by-design table refused: %s" % e)
    result["by_design_table"] = None if info is None else {
        "sha256": info.get("sha256"), "path": cmp.run_relative(info.get("path"), run_dir), "cases": info.get("cases")}
    for scen, op in ops:
        entry = result["scenarios"].get(scen.get("id"), {}).get("ops", {}).get(op["id"])
        if entry is not None:
            judge_by_design(entry, scen, op, recs, cases)


def service_checks_stage(result, run_dir, manifest, recs, prechecks):
    """Stage of compare_graphlets: the positive service checks and the scenario facts they read."""
    facts = scenario_facts(run_dir, manifest)
    result["facts"] = facts
    checks = service_checks(recs, manifest, facts, prechecks)
    result["service_checks"] = checks
    result["required_service_checks"] = required_service_checks(manifest)
    # a probe op (W3 only: nothing to compare across ways) is judged by its own probe check
    by_check = dict((c["check"], c) for c in checks)
    for scen in gp.graphlet_scenarios(manifest):
        for op in gp.scenario_ops(scen):
            probe = PROBES.get(op.get("probe"))
            entry = result["scenarios"].get(scen.get("id"), {}).get("ops", {}).get(op["id"])
            if not probe or entry is None or entry.get("outcome") in ("precondition_failed", "disagree"):
                continue
            entry["probe"] = op["probe"]
            failed = [f for f in by_check[probe[0]]["failures"] if (f["scenario"], f["op"]) == (scen.get("id"),
                                                                                                 op["id"])]
            if failed:
                entry["outcome"] = "disagree"
                entry["disagreeing"] = ["W3-probe %s" % op["probe"]]
                entry["reason"] = failed[0]["detail"]
            elif entry["ways"].get("W3") != "error":
                entry["outcome"] = "inconclusive" if entry["ways"].get("W3") == "missing" else "disagree"


for _stage in (by_design_stage, service_checks_stage):
    if _stage not in cg.STAGES:
        cg.STAGES.append(_stage)
