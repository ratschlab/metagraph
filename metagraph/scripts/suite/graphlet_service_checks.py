#!/usr/bin/env python3
"""The service's own behaviour in area G, checked positively (DESIGN-v8 §5, D-R6), and the declared per-way
differences verified as by-design pairs (D-R12). Two stages of compare_graphlets.py (registered on import).

Positive service checks, functions over W3's records and the scenario's facts (scenario_facts: from W1's document
by the suite's own reader, a scenario's graphlet.facts stating them where the corpus says so):
  walk_display       every graphlet_walks / graphlet_walk answer of W3 carries walk_display "first_parent";
  label_evidence     a W3 graphlet_compare answer carries label_evidence "lower_bound" exactly when the service adds
                     its label fields (as hint below: the verdict qualified, the notes naming label_lists, a side's
                     label lists cut), never on another compare (lane D's finding, 2026-10-20);
  qualified_cause    present exactly with label_evidence, and consistent with the library's answer to the same call
                     (W4's, else W1's): it states the cap of each side whose lists were cut (labels.max_labels_per_
                     node of the side's request; ops_runtime._label_cap names those sides only)
                     and "no difference within the recorded evidence" exactly when the library counts no difference;
  hint               a W3 compare answer carries hint exactly when its verdict is qualified (comparable), its notes
                     name label_lists and a side's label lists were cut: the side's cut_lists, or a node with more
                     labels (labels_total) than the side's cap, labels.max_labels_per_node of its scenario's strategy
                     (64, the engine's default, when it states none), not a fixed 64 (the service's rule,
                     ops_runtime._label_cap / _finish; service_settings.graphlet.compare_hint, async's limits block);
  in_service         a W3 graphlet_summary of a body over 8 MiB carries in_service false and a download; no other
                     summary carries in_service false;
  max_bytes_refused  a W3 call with max_bytes above the service's 262,144 is refused bad_argument, and the library
                     ways answer the same call (the library's ceiling is 16 MiB);
  forbidden          a W3 call the op plan sends with a wrong access_token (op probe "wrong_token") is forbidden;
  bad_handle         a W3 call with a malformed handle (op probe "bad_handle", FORMAT v2 D2) is bad_handle;
  library_commit     every W3 call reports the deployed library_commit and the manifest's pin (the pre-checks).
Each check: pass, fail (naming every call that breaks it), n/a (its precondition holds on no call of the run: no
compare call, no body over 8 MiB, no probe; the detail says which), or inconclusive (a fact it needs is unknown).
blocking_service_checks says which outcomes block (VERIFY-C VC2): any fail, a required check that did not pass, an
exercised check that is inconclusive; a check the declared file sets aside is a note.

By-design pairs: an op naming a case (`by_design`: an id or a list) of the by-design table (compare.table_for: the
run's snapshot table, else service_by_design.json) whose `graphlet` block states the pair (compare.graphlet_case_
problems: tools, args predicates over args_resolved, the library ways' expectation, W3's expectation). The pair holds
when the op's tool is listed, the args predicates hold, every library way the op runs matches `library` and W3
matches `service`; then W3-W4 and W3-W5 may differ and the op is service_by_design (verified), its library pairs
still compared. A pair that does not hold is a disagreement naming the side (table, request, library, service).
A case whose service side pins a key of the library's own answer (service_fields naming a key of compare_graphlets.
ANSWER_KEYS[tool], e.g. graphlet_compare's export: W3 rewrites the library's line, async's 2026-10-21 case) narrows
what it excuses to those keys: W3's value is read from its answer and must be the pinned one, and a W3-W4 answer
difference is excused only when the two answers agree with those keys left out of both (pinned_answer_keys).
"""

import os
import re
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
DEFAULT_LABELS_PER_NODE = 64   # the engine's labels.max_labels_per_node when a strategy states none (async's U2)
IN_SERVICE_BYTES = 8 * 1024 * 1024
NO_DIFFERENCE = "no difference within the recorded evidence"
# FORMAT v2 D2: an op's probe "wrong_token" | "bad_handle" -> (the service check, the code W3 must answer)
PROBES = {"wrong_token": ("forbidden", "forbidden"), "bad_handle": ("bad_handle", "bad_handle")}
# the pairs a by-design case lets differ (the service side against the library and the oracle)
BY_DESIGN_PAIRS = ("W3-W4", "W3-W5")
# the kinds of disagreement a verified case excuses on those pairs (an answer or an error that differs; a receipt is
# an answer): never a file (kind file: compare_graphlets compares a file by its content first, whatever the receipts)
BY_DESIGN_KINDS = ("error", "answer")


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
        out["cut_lists"] = any(s.get("labels_truncated") for s in segs) or engine_states_cut_lists(text)
        out["max_labels_per_node"] = max(s.get("labels_total") or 0 for s in segs)
    return out


def engine_states_cut_lists(text):
    """Whether the document states cut label lists as the engine writes them (F-D3, the rehearsal of 2026-10-19):
    a K record of kind label_lists ("N recorded label list(s) are cut at the cap", SPEC-core §7.5.2 K) or the O
    record's label_evidence `l` (lower_bound). A segment's entry list can be whole (labels_truncated false) while other
    recorded lists (end sets, events) are cut: on every cut-list graphlet of that rehearsal no segment showed one."""
    for line in text.split("\n"):
        f = line.split(" ")
        if f[0] == "K" and len(f) > 2 and f[2] == "label_lists":
            return True
        if f[0] == "O" and len(f) > 3 and f[3] == "l":
            return True
    return False


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
        # the check's precondition holds on no call of the run (VERIFY-C VC2: n/a, never a failure)
        return {"check": check, "outcome": "n/a", "exercised": False, "calls": 0, "failures": [],
                "detail": "n/a: no W3 call of the run %s" % what}
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
        sf = service_fields(r)
        want, why = _label_fields_expected(recs, facts, sid, op, r)
        if want is None:
            unknown.append("%s %s: %s" % (sid, op["id"], "no compare verdict" if why == "verdict" else
                                          "cut label lists unknown"))
            continue
        if want and sf.get("label_evidence") != "lower_bound":
            bad.append(_fail(sid, op, "%s, label_evidence %s" % (why, gp.compact(sf.get("label_evidence",
                                                                                         gp.MISSING)))))
        elif not want and "label_evidence" in sf:
            bad.append(_fail(sid, op, "%s, yet label_evidence %s" % (why, gp.compact(sf["label_evidence"]))))
    return _result("label_evidence", calls, bad, unknown, "of graphlet_compare answered")


def _label_fields_expected(recs, facts, sid, op, r):
    """(whether the service adds its label fields -- label_evidence, qualified_cause, hint -- to this compare answer,
    why) as it does (ops_runtime._label_cap / _finish): the verdict qualified (comparable, W3's answer, else the
    library's), the library's notes naming label_lists and a side's label lists cut (side_cut_lists: its cut_lists,
    or more labels on a node than its own cap). (None, "verdict" | "facts") when what decides it is unknown."""
    verdict, notes = _compare_verdict(recs, sid, op["id"], r)
    if verdict is None:
        return None, "verdict"
    if verdict != "qualified":
        return False, "verdict %s" % gp.compact(verdict)
    if not any(isinstance(n, str) and "label_lists" in n for n in notes):
        return False, "verdict qualified, no label_lists note"
    sides = [(x, side_cut_lists(facts.get(x, {}))) for x in (sid, _b_scenario(sid, op))]
    cut = [x for x, v in sides if v]
    if cut:
        return True, "verdict qualified, label lists cut on %s" % " and ".join(cut)
    if any(v is None for _x, v in sides):
        return None, "facts"
    return False, "verdict qualified, no cut label list"


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
            if not side_cut_lists(facts.get(side, {})):
                continue   # the service names the caps of the sides whose lists were cut only (VERIFY-D VD11)
            cap = facts.get(side, {}).get("cap")
            # the cap as a number of its own (VERIFY-C VC12): 64 is not stated by "640" or "1964"
            if isinstance(cap, int) and not re.search(r"(?<![0-9])%d(?![0-9])" % cap, cause):
                bad.append(_fail(sid, op, "qualified_cause %s does not state the cap %d of %s" % (
                    gp.compact(cause), cap, side)))
                break
    return _result("qualified_cause", calls, bad, unknown, "of graphlet_compare answered")


def side_cut_lists(f):
    """Whether a graphlet's label lists were cut (True / False / None unknown): its cut_lists (the document's
    segments, its K label_lists or O lower_bound record), or a node with more labels than its cap (labels_total
    above the scenario's labels.max_labels_per_node, DEFAULT_LABELS_PER_NODE when the strategy states none)."""
    cut = f.get("cut_lists")
    n = f.get("max_labels_per_node")
    cap = f.get("cap") if isinstance(f.get("cap"), int) and not isinstance(f.get("cap"), bool) else \
        DEFAULT_LABELS_PER_NODE
    over = n > cap if isinstance(n, int) and not isinstance(n, bool) else None
    if cut or over:
        return True
    if cut is None and over is None:
        return None
    return False


def _compare_verdict(recs, sid, oid, r):
    """(comparable, notes) of a compare call: W3's own answer, else the library's to the same call."""
    for ans in (r.get("answer"), _library_answer(recs, sid, oid)):
        if isinstance(ans, dict) and "comparable" in ans:
            notes = ans.get("notes")
            return ans["comparable"], notes if isinstance(notes, list) else []
    return None, []


def check_hint(recs, manifest, facts):
    calls, bad, unknown = 0, [], []
    for sid, op, r in _compare_calls(recs, manifest):
        calls += 1
        has = "hint" in service_fields(r)
        verdict, notes = _compare_verdict(recs, sid, op["id"], r)
        if verdict is None:
            unknown.append("%s %s: no compare verdict" % (sid, op["id"]))
            continue
        named = any(isinstance(n, str) and "label_lists" in n for n in notes)
        if verdict != "qualified" or not named:
            if has:
                bad.append(_fail(sid, op, "verdict %s%s, hint present" % (
                    gp.compact(verdict), "" if named else ", no label_lists note")))
            continue
        sides = [(x, side_cut_lists(facts.get(x, {}))) for x in (sid, _b_scenario(sid, op))]
        cut = [x for x, v in sides if v]
        if not cut and any(v is None for _x, v in sides):
            unknown.append("%s %s: cut label lists unknown" % (sid, op["id"]))
            continue
        if bool(cut) != has:
            caps = ", ".join("%s cap %s" % (x, facts.get(x, {}).get("cap") or DEFAULT_LABELS_PER_NODE)
                             for x, _v in sides)
            bad.append(_fail(sid, op, "verdict qualified, label lists %s (%s), hint %s" % (
                "cut on " + " and ".join(cut) if cut else "not cut", caps, "present" if has else "absent")))
    return _result("hint", calls, bad, unknown, "of graphlet_compare answered")


def check_in_service(recs, manifest, facts):
    calls, bad, over = 0, [], 0
    for sid, op, r in _w3_calls(recs, manifest, ("graphlet_summary",)):
        if r.get("error"):
            continue
        calls += 1
        sf = service_fields(r)
        if facts.get(sid, {}).get("over_8mib"):
            over += 1
            if sf.get("in_service") is not False or not sf.get("download"):
                bad.append(_fail(sid, op, "a body over 8 MiB: in_service %s, download %s" % (
                    gp.compact(sf.get("in_service", gp.MISSING)), "present" if sf.get("download") else "absent")))
        elif sf.get("in_service") is False:
            bad.append(_fail(sid, op, "in_service false on a body within 8 MiB"))
    # a summary within 8 MiB that says in_service false fails the check; without one over 8 MiB it is n/a otherwise
    return _result("in_service", calls if bad else over, bad, what="of graphlet_summary on a body over 8 MiB")


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
            _check_probe("bad_handle", "bad_handle", "bad_handle", recs, manifest),
            check_library_commit(prechecks)]


def blocking_service_checks(checks, required=(), set_aside=()):
    """[(check, why)] of the service checks that block the gate and the report's verdict (VERIFY-C VC2): every check
    that failed, whether or not the op plan was expected to exercise it; a required check (required_service_checks)
    that did not pass (n/a included: the plan meant to exercise it); an exercised check that is inconclusive. A check
    in |set_aside| (the declared file's graphlet_service_checks_not_required) never blocks."""
    out = []
    for c in checks or []:
        name, o = c.get("check"), c.get("outcome")
        if name in set_aside:
            continue
        if o == "fail":
            out.append((name, "fail: %s" % c.get("detail")))
        elif name in required and o != "pass":
            out.append((name, "required, %s: %s" % (o, c.get("detail"))))
        elif o == "inconclusive" and c.get("exercised"):
            out.append((name, "inconclusive: %s" % c.get("detail")))
        elif o not in ("pass", "fail", "n/a", "inconclusive"):
            out.append((name, "outcome %r outside the vocabulary" % (o,)))
    return out


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
    pinned = pinned_answer_keys(rec.get("tool"), want)
    for k, v in sorted((want.get("service_fields") or {}).items()):
        if k in pinned:   # a key of the library's own answer: W3's value is the one in its answer
            ans = rec.get("answer") if isinstance(rec.get("answer"), dict) else {}
            got = ans.get(k, gp.MISSING)
        else:
            got = service_fields(rec).get(k, gp.MISSING)
        if got is gp.MISSING or not sc.json_equal(got, v):
            return False, "%s's service field %s %s, the case says %s" % (way, k, gp.compact(got), gp.compact(v))
    return True, ""


def pinned_answer_keys(tool, side):
    """The keys a case side's service_fields pins that are keys of the library's own answer of |tool|
    (compare_graphlets.ANSWER_KEYS): W3 changes the library's value there (module docstring). -> set."""
    keys = cg.ANSWER_KEYS.get(tool) or frozenset()
    fields = (side or {}).get("service_fields") if isinstance((side or {}).get("service_fields"), dict) else {}
    return set(k for k in fields if k in keys)


def _without(rec, keys):
    """A copy of a record whose answer leaves out |keys|."""
    out = dict(rec)
    if isinstance(rec.get("answer"), dict):
        out["answer"] = dict((k, v) for k, v in rec["answer"].items() if k not in keys)
    return out


def narrow_pair_agrees(pair, sid, oid, recs, keys):
    """(agrees, detail) of a disagreeing W3-W<n> pair compared again with |keys| left out of both answers (a case
    pinning keys of the library's answer excuses those keys only): the same comparison (compare_graphlets.
    compare_records), never W5 (its answer differences are of kind oracle, which no case excuses)."""
    a, b = pair["pair"].split("-", 1)
    ra, rb = recs.get(a, {}).get((sid, oid)), recs.get(b, {}).get((sid, oid))
    if ra is None or rb is None or b not in gp.LIBRARY_WAYS:
        return False, "%s: not compared again" % pair["pair"]
    again = cg.compare_records(a, _without(ra, keys), b, _without(rb, keys), pair.get("what"))
    if again.get("outcome") == "agree":
        return True, ""
    return False, "%s still differs without %s: %s" % (pair["pair"], ", ".join(sorted(keys)),
                                                      again.get("pointer") or again.get("detail") or "?")


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
            # the case's predicates hold on every way's arguments (VERIFY-C VC11), not on W3's alone
            why = ""
            for w in [x for x in gp.WAYS if x in ways and x != "W5"]:
                ref = recs.get(w, {}).get((sid, oid))
                if ref is None:
                    continue
                ok, why = cmp.request_holds({"id": cid, "request": g.get("args")},
                                            {"mcp": {"arguments": ref.get("args_resolved") or {}}}, None)
                if not ok:
                    why = "%s on %s's args_resolved" % (why, w)
                    break
                why = ""
            if why:
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
                pinned = pinned_answer_keys(op.get("tool"), g.get("service"))
                if pinned:
                    res["pinned_answer_keys"] = sorted(pinned)
        results.append(res)
    entry["by_design"] = results
    if entry["outcome"] in ("precondition_failed",) or entry.get("smoke"):
        return entry
    # a verified case excuses only the declared difference: a W3-W4 or W3-W5 disagreement of kind error / answer,
    # never differing args_resolved, a selector (VERIFY-C VC11) or a file's content (kind file)
    kinds = dict((p.get("pair"), p.get("kind")) for p in entry.get("pairs") or [] if p.get("outcome") == "disagree")
    # a verified case pinning keys of the library's answer excuses an answer difference in those keys only; any
    # other verified case the declared pair as a whole
    broad = any(r["verified"] and not r.get("pinned_answer_keys") for r in results)
    narrow = set(k for r in results if r["verified"] for k in r.get("pinned_answer_keys") or [])
    by_pair = dict((p.get("pair"), p) for p in entry.get("pairs") or [] if p.get("outcome") == "disagree")

    def excused(p):
        if p not in BY_DESIGN_PAIRS or kinds.get(p) not in BY_DESIGN_KINDS:
            return False
        if broad:
            return True
        if narrow and kinds.get(p) == "answer":
            ok, why = narrow_pair_agrees(by_pair[p], sid, oid, recs, narrow)
            if not ok:
                entry.setdefault("by_design_narrow", []).append(why)
            return ok
        return False
    rest = [p for p in entry.get("disagreeing") or [] if not excused(p)]
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
    """Stage of compare_graphlets: judges every op that names by-design cases, with the run's tables: the
    service-by-design table and area G's own (service_by_design_g.json, VERIFY-D VD1), each the run's snapshot's,
    else the tree's; both identities recorded (by_design_table, by_design_g_table: a v8 policy pins both)."""
    ops = [(scen, op) for scen in gp.graphlet_scenarios(manifest) for op in gp.scenario_ops(scen)
           if op.get("by_design")]
    if not ops:
        return
    path, warnings = cmp.table_for(run_dir, None)
    result.setdefault("warnings", []).extend(w for w in warnings if w not in result.get("warnings", []))
    cases = {}
    for key, p_ in (("by_design_table", path), ("by_design_g_table", cmp.g_table_for(run_dir))):
        try:
            got, info = cmp.load_by_design(p_) if p_ else ({}, None)
        except cmp.TableError as e:
            got, info = {}, None
            result["warnings"].append("by-design table refused: %s" % e)
        for cid in sorted(set(got) & set(cases)):
            result["warnings"].append("by-design case %s is in both tables: %s's is read" % (cid, key))
        cases.update(got)
        result[key] = None if info is None else {
            "sha256": info.get("sha256"), "path": cmp.run_relative(info.get("path"), run_dir),
            "cases": info.get("cases")}
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
