#!/usr/bin/env python3
"""The five-way compare of area G (manifest v8; DESIGN-v8 §2, §5, FORMAT-v8 §1-§5): per op call of every graphlet
scenario, the answers of the ways that run it, after the scenario's pre-checks (graphlet_precheck.py).

    compare_graphlets.py RUN_DIR [--out FILE]     (default RUN_DIR/graphlet_compare.json; the pre-check records go to
                                                   RUN_DIR/arms/precheck/graphlet_prechecks.jsonl)

Per op call, in this order:
1. the ways it runs in (its `ways`, else the scenario's ways.default, else all five); W3 left out is `not_exposed`
   with the owner's reason (the op's not_exposed.W3, else the scenario's or the manifest's graphlets.not_exposed by
   tool or tool.argument); W3 left out without a reason is a manifest problem (inconclusive);
2. a way that recorded nothing is missing; a selector that failed in a way did not run the call there (the same
   failure in every way: inconclusive, the selector named; in some ways only: a disagreement);
3. args_resolved equal across the ways (access tokens and cursors left out): a selector that resolved differently is
   a disagreement of the pair, named as such, and that pair's answers are not compared;
4. each answer canonicalised (canonical_answer): W3's service fields (the keys outside ANSWER_KEYS, the library's
   per-tool answer keys, frozen here) stripped and kept for the service checks; handles to symbols (<handle>,
   <handle_b>, <handle:N>), cursors removed, parsed JSON values only (numbers by value); the `local` block compared
   only between library ways of one (commit, interpreter, work model) (D-R10); an export's `path` and `bytes` left
   out (the export directory, the handle and the envelope's timing make them differ by construction);
5. whole assembled answers compared by pairs: W1-W2 (transport), W1-W4 (fetch path), W3-W4 (service wrapper), and
   each of them against W5 (library vs independent reading) on W5's `checked` pointers only; W5's other pointers
   are consistency-only and listed;
6. errors by code and structured fields (never message or hint text; result_too_large's `bytes` depends on the
   handle's length);
7. a deliberate cut (the op's `cut_of`: the op id of the same call uncut) is judged in the suite's words, in each
   way against that way's uncut answer (is_cut): the same uncut answer (the uncut op's own comparison), every string
   or list of the cut answer a prefix of the uncut one, every key it lacks an optional one (OPTIONAL_KEYS, or named
   in its fields_cut / cut_fields), every key it adds a cut marker (CUT_MARKERS); the library ways' cut answers
   must still be equal (one handle length, one ceiling).
Outcomes per op: agree; disagree, naming the pairs; inconclusive (a way missing, a selector failed everywhere, a
manifest problem); precondition_failed (the scenario's preconditions failed: no op is compared). The scenario: fail
(a pre-check failed, or an op disagrees), else precondition_failed, else inconclusive, else pass. graphlet_list /
graphlet_free are a smoke outside the five-way report (D-R15): each way answered, nothing compared. The drift lane
(arms/drift-lib) is compared with W4 and listed, report-only.
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
import graphlet_precheck as gp  # noqa: E402
import suite_common as sc  # noqa: E402

COMPARE_FILE = "graphlet_compare.json"
COMPARE_VERSION = 1
WAYS = gp.WAYS

# The library's answer keys per tool (GraphletTools, metagraph/api/python/metagraph/traverse/mcp_tools.py), the keys
# a tool's answer may carry on any way: everything else a W3 answer carries is the service's own and is stripped
# (service_fields) before the cross-way comparison. Frozen here (D-R6); tests/test_compare_graphlets.py holds it to
# the library's answers on the mini_refseq fixtures. PAGE_KEYS: a list page's own keys (_page: total, rows,
# next_cursor; under local limits complete, total_at_least, filtered_so_far, stop and the local block; a single row
# larger than the page row_truncated, cut_fields).
PAGE_KEYS = frozenset(("total", "rows", "next_cursor", "complete", "total_at_least", "filtered_so_far", "stop",
                       "row_truncated", "cut_fields", "local"))
RECEIPT_KEYS = frozenset(("fields_cut", "local"))
ANSWER_KEYS = {
    "graphlet_summary": PAGE_KEYS | {"handle", "evidence", "parent", "summary", "outcome", "resource_stop"},
    "graphlet_walks": PAGE_KEYS | {"handle", "arm", "evidence", "filtered", "hint"},
    "graphlet_walk": PAGE_KEYS | {"handle", "arm", "walk", "length_bp", "evidence", "continuation"},
    "graphlet_support": PAGE_KEYS | {"handle", "arm", "walk", "evidence"},
    "graphlet_labels": PAGE_KEYS | {"handle", "evidence", "name", "ref"},
    "graphlet_splits": PAGE_KEYS | {"handle", "arm", "evidence"},
    "graphlet_claims": PAGE_KEYS | {"handle", "arm", "evidence", "filtered", "hint", "caveat"},
    "graphlet_sequence": frozenset(("handle", "arm", "orientation", "with_seed", "from", "length", "ceiling_bytes",
                                    "evidence", "to", "sequence", "truncated", "next_from", "local")),
    "graphlet_export": RECEIPT_KEYS | {"path", "bytes", "records", "parsed", "validation"},
    "graphlet_compare": PAGE_KEYS | {"a", "b", "mode", "comparable", "reason", "equal", "depth_used", "counts",
                                     "notes", "support", "scopes", "evidence", "export", "fields_cut"},
    "graphlet_subtrie": RECEIPT_KEYS | {"handle", "of", "view", "summary", "evidence"},
    "graphlet_list": PAGE_KEYS,
    "graphlet_free": frozenset(("freed",)),
    "graphlet_save": RECEIPT_KEYS | {"path", "bytes", "parsed", "validation"},
    "graphlet_load": RECEIPT_KEYS | {"handle", "summary", "summary_stop", "parsed", "validation"},
    "traverse_continue": RECEIPT_KEYS | {"request", "parent", "ceiling_bytes", "notes", "loss_budget",
                                         "branch_budget", "loss_budget_labels", "handle", "identity", "evidence",
                                         "summary", "summary_stop", "parsed", "server_summary", "hint", "stop",
                                         "index", "delivery", "graphlet_bytes", "coordinates_note"},
}
# the keys whose string values are handles (record UUIDs on the service, g_ + 12 hex in the library) at any depth
# (parent.handle, view.of, evidence.view.of, a list's rows); graphlet_compare's a and b only at the top (a compare
# row's a and b are the two sides' values)
HANDLE_KEYS = ("handle", "of", "freed")
TOP_HANDLE_KEYS = ("a", "b")
CURSOR_KEYS = ("next_cursor", "cursor")
# argument names left out of args_resolved before it is compared (per-way credentials, cursors)
ARG_SKIP = ("access_token", "access_token_b", "cursor")
# error fields that are text or advice, never compared (FORMAT-v8 §3: message text stays in the raw page)
ERROR_TEXT_FIELDS = ("message", "hint")
# per error code, structured fields that depend on the environment (result_too_large's size of what did not fit
# counts the handle's own bytes)
ERROR_ENV_FIELDS = {"result_too_large": ("bytes",)}
# per tool, answer keys whose value depends on where the way ran, left out of the comparison
ENV_KEYS = {"graphlet_export": ("path", "bytes")}
# a cut answer may add these keys (the ceiling's markers); a cut answer may lack these keys (per tool) and the keys
# its own fields_cut / cut_fields name
CUT_MARKERS = frozenset(("fields_cut", "row_truncated", "cut_fields", "caveat_effects_cut", "caveats_truncated",
                         "next_cursor", "truncated", "next_from"))
OPTIONAL_KEYS = {
    "graphlet_summary": frozenset(("effect", "message", "top_labels", "caveats")),
    "graphlet_compare": frozenset(("evidence", "export")),
}
SMOKE_TOOLS = ("graphlet_list", "graphlet_free")
# the pairs whose agreement makes W1 = W2 = W3 = W4 (a spanning tree), each with what it isolates
PAIRS = (("W1", "W2", "transport"), ("W1", "W4", "fetch path"), ("W3", "W4", "service wrapper"))
W5_WHAT = "library vs independent reading"
OP_OUTCOMES = ("agree", "disagree", "service_by_design", "inconclusive", "precondition_failed")
SCENARIO_OUTCOMES = ("fail", "precondition_failed", "inconclusive", "pass")


# ---------------------------------------------------------------------------------------------------------------
# canonical answers

def split_service_fields(tool, answer):
    """(the library's part of a W3 answer, the keys outside ANSWER_KEYS[tool]): the service's own fields are not
    compared across ways (they are checked positively, graphlet_service_checks). A tool outside the table keeps its
    whole answer."""
    keys = ANSWER_KEYS.get(tool)
    if not isinstance(answer, dict) or keys is None:
        return answer, {}
    return (dict((k, v) for k, v in answer.items() if k in keys),
            dict((k, v) for k, v in answer.items() if k not in keys))


def library_part(way, rec):
    """The answer of a record as the comparison reads it: W3's without its service fields, the others' as they
    are."""
    ans = rec.get("answer")
    return split_service_fields(rec.get("tool"), ans)[0] if way == "W3" else ans


def whitelist_gaps(tool, answer):
    """The keys of a library way's answer outside ANSWER_KEYS[tool] (the frozen table behind the library)."""
    keys = ANSWER_KEYS.get(tool)
    if keys is None or not isinstance(answer, dict):
        return []
    return sorted(k for k in answer if k not in keys)


def handle_roles(tool, answer):
    """{raw handle: symbol} of the handles an answer names in a known role: graphlet_compare's a -> <handle>, b ->
    <handle_b> (a handle that is both is <handle>); every other tool's handle -> <handle>."""
    roles = {}
    if not isinstance(answer, dict):
        return roles
    if tool == "graphlet_compare":
        for key, sym in (("a", "<handle>"), ("b", "<handle_b>")):
            if isinstance(answer.get(key), str):
                roles.setdefault(answer[key], sym)
    elif isinstance(answer.get("handle"), str):
        roles[answer["handle"]] = "<handle>"
    return roles


def _canon(obj, roles, drop_local, top=False, tool=None):
    if isinstance(obj, dict):
        out = {}
        for k in sorted(obj, key=str):
            v = obj[k]
            if k in CURSOR_KEYS or (drop_local and k == "local") or (top and k in ENV_KEYS.get(tool, ())):
                continue
            if top and tool == "graphlet_compare" and k in TOP_HANDLE_KEYS and isinstance(v, str):
                out[k] = "<handle>" if k == "a" else "<handle_b>"
            elif (k in HANDLE_KEYS or (top and k in TOP_HANDLE_KEYS)) and isinstance(v, str):
                if v not in roles:
                    roles[v] = "<handle:%d>" % (sum(1 for s in roles.values() if s.startswith("<handle:")) + 1)
                out[k] = roles[v]
            else:
                out[k] = _canon(v, roles, drop_local, tool=tool)
        if drop_local and isinstance(out.get("limitations"), dict) and "local" in out["limitations"]:
            out["limitations"] = dict((k, v) for k, v in out["limitations"].items() if k != "local")
        return out
    if isinstance(obj, list):
        return [_canon(v, roles, drop_local, tool=tool) for v in obj]
    if isinstance(obj, str):
        for raw, sym in roles.items():
            if len(raw) >= 8 and raw in obj and obj != raw:
                obj = obj.replace(raw, sym)
        return obj
    return obj


def canonical_answer(tool, answer, keep_local=False):
    """An answer as the ways compare it: handles to symbols (handle_roles, then <handle:N> in sorted-key order for
    any other handle-valued key; a raw handle inside a string too), cursors removed, the local block (and
    evidence.limitations.local) removed unless |keep_local|, the tool's ENV_KEYS removed."""
    if not isinstance(answer, dict):
        return answer
    return _canon(answer, handle_roles(tool, answer), not keep_local, top=True, tool=tool)


def canonical_error(error, keep_local=False):
    """An error as the ways compare it: its code and its structured fields without text (ERROR_TEXT_FIELDS), without
    the code's environment fields (ERROR_ENV_FIELDS), the local block unless |keep_local|, handles to symbols."""
    if not isinstance(error, dict):
        return error
    code = error.get("code")
    fields = error.get("fields") if isinstance(error.get("fields"), dict) else {}
    skip = set(ERROR_TEXT_FIELDS) | set(ERROR_ENV_FIELDS.get(code, ()))
    fields = dict((k, v) for k, v in fields.items() if k not in skip)
    return {"code": code, "fields": _canon(fields, {}, not keep_local)}


def canonical_args(args):
    """args_resolved as the ways compare it: without per-way credentials and cursors (ARG_SKIP)."""
    if not isinstance(args, dict):
        return args
    return dict((k, v) for k, v in args.items() if k not in ARG_SKIP)


def same_library_env(ra, rb):
    """Whether two library records ran one (commit, interpreter, work model): only then is `local` compared (D-R10)."""
    la, lb = ra.get("library") or {}, rb.get("library") or {}
    if not isinstance(la, dict) or not isinstance(lb, dict):
        return False
    keys = ("commit", "python", "work_model")
    return all(la.get(k) is not None for k in keys) and all(la.get(k) == lb.get(k) for k in keys)


def comparable_view(way, rec, keep_local=False):
    """{"error": canonical error} or {"answer": canonical answer} of a record."""
    if rec.get("error"):
        return {"error": canonical_error(rec["error"], keep_local)}
    return {"answer": canonical_answer(rec.get("tool"), library_part(way, rec), keep_local)}


# ---------------------------------------------------------------------------------------------------------------
# a cut of an uncut answer, in the suite's words

def is_cut(cut, uncut, optional=frozenset(), pointer=""):
    """None when |cut| is a cut of |uncut|, else (pointer, why): every string a prefix, every list a prefix whose
    items are cuts of the uncut items, every key the cut lacks optional (|optional|, or named in its fields_cut /
    cut_fields at that level), every key it adds a cut marker (CUT_MARKERS), every other value equal."""
    if isinstance(cut, dict) and isinstance(uncut, dict):
        named = set(x for k in ("fields_cut", "cut_fields") for x in (cut.get(k) or []) if isinstance(x, str))
        for k in sorted(set(cut) | set(uncut), key=str):
            p = pointer + "/" + gp.pointer_escape(k)
            if k not in uncut:
                if k in CUT_MARKERS:
                    continue
                return p, "a key only the cut answer has, not a cut marker"
            if k not in cut:
                if k in optional or k in named:
                    continue
                return p, "a key the cut answer lacks, not an optional one"
            r = is_cut(cut[k], uncut[k], optional, p)
            if r:
                return r
        return None
    if isinstance(cut, list) and isinstance(uncut, list):
        if len(cut) > len(uncut):
            return pointer or "/", "a list longer than the uncut one (%d > %d)" % (len(cut), len(uncut))
        for i, (x, y) in enumerate(zip(cut, uncut)):
            r = is_cut(x, y, optional, "%s/%d" % (pointer, i))
            if r:
                return r
        return None
    if isinstance(cut, str) and isinstance(uncut, str):
        return None if uncut.startswith(cut) else (pointer or "/", "a string that is not a prefix of the uncut one")
    if isinstance(cut, (dict, list)) or isinstance(uncut, (dict, list)):
        return pointer or "/", "a value of another type than the uncut one"
    return None if gp.leaf_equal(cut, uncut) else (pointer or "/", "a value other than the uncut one (%s vs %s)" % (
        gp.compact(cut, 60), gp.compact(uncut, 60)))


# ---------------------------------------------------------------------------------------------------------------
# the manifest's op plan

def scenario_by_id(manifest):
    return dict((s.get("id"), s) for s in gp.graphlet_scenarios(manifest))


def not_exposed_reason(scenario, op, manifest):
    """The owner's reason W3 does not run |op| (FORMAT-v8 §1, DESIGN-v8 §2's table): the op's not_exposed (an object
    {"W3": reason} or a text), else graphlet.not_exposed of the scenario, else the manifest's graphlets.not_exposed,
    each keyed by the tool or by tool.argument for an argument the op sends; None without one."""
    ne = op.get("not_exposed")
    if isinstance(ne, dict) and isinstance(ne.get("W3"), str) and ne["W3"].strip():
        return ne["W3"].strip()
    if isinstance(ne, str) and ne.strip():
        return ne.strip()
    tool = op.get("tool")
    keys = [tool] + ["%s.%s" % (tool, a) for a in sorted(op.get("args") or {})]
    for table in ((scenario.get("graphlet") or {}).get("not_exposed"),
                  ((manifest or {}).get("graphlets") or {}).get("not_exposed")
                  if isinstance((manifest or {}).get("graphlets"), dict) else None):
        if isinstance(table, dict):
            for k in keys:
                if isinstance(table.get(k), str) and table[k].strip():
                    return table[k].strip()
    return None


# ---------------------------------------------------------------------------------------------------------------
# the comparison of one op call

ANSWER_CUT = 240


def _pair(a, b, what, outcome, detail=None, diff=None, kind="answer", views=None):
    """A pair entry; a disagreeing one names the first difference (pointer, both values) and, with |views|, both
    whole canonical answers compactly (`answers`, cut at ANSWER_CUT characters)."""
    out = {"pair": "%s-%s" % (a, b), "what": what, "outcome": outcome, "kind": kind}
    if detail:
        out["detail"] = detail
    if diff:
        out["pointer"] = diff[0] or "/"
        out[a] = gp.compact(diff[1])
        out[b] = gp.compact(diff[2])
    if views is not None and outcome == "disagree":
        out["answers"] = {a: gp.compact(views[0], ANSWER_CUT), b: gp.compact(views[1], ANSWER_CUT)}
    return out


def _args_difference(ra, rb):
    return gp.first_difference(canonical_args(ra.get("args_resolved")), canonical_args(rb.get("args_resolved")))


def compare_records(a, ra, b, rb, what):
    """A pair entry for two recorded calls of one op (selector, arguments, then answers)."""
    sa, sb = ra.get("selector_error"), rb.get("selector_error")
    if sa or sb:
        if sa and sb and gp.compact(sa) == gp.compact(sb):
            return _pair(a, b, what, "not_compared", "the selector failed in both: %s" % gp.compact(sa),
                         kind="selector")
        return _pair(a, b, what, "disagree", "the selector failed in %s only: %s" % (
            a if sa else b, gp.compact(sa or sb)), kind="selector")
    d = _args_difference(ra, rb)
    if d:
        return _pair(a, b, what, "disagree", "args_resolved differ (a selector resolved differently)", d,
                     kind="args_resolved")
    keep = a in gp.LIBRARY_WAYS and b in gp.LIBRARY_WAYS and same_library_env(ra, rb)
    va, vb = comparable_view(a, ra, keep), comparable_view(b, rb, keep)
    if ("error" in va) != ("error" in vb):
        err, ok = (a, b) if "error" in va else (b, a)
        code = (va.get("error") or vb.get("error") or {}).get("code")
        return _pair(a, b, what, "disagree", "%s answered the error %s, %s an answer" % (err, code, ok),
                     ("/", va, vb), kind="error", views=(va, vb))
    d = gp.first_difference(va, vb)
    if d:
        return _pair(a, b, what, "disagree", None, d, kind="error" if "error" in va else "answer", views=(va, vb))
    return _pair(a, b, what, "agree")


def compare_w5(x, rx, r5):
    """The pair entry of way |x| against W5 on W5's checked pointers (FORMAT-v8 §4), with the pointers it checked."""
    oracle = r5.get("oracle") if isinstance(r5.get("oracle"), dict) else {}
    checked = oracle.get("checked") if isinstance(oracle.get("checked"), dict) else {}
    if rx.get("selector_error") or r5.get("selector_error"):
        return _pair(x, "W5", W5_WHAT, "not_compared", "the selector failed in %s" % (
            x if rx.get("selector_error") else "W5"), kind="selector")
    d = _args_difference(rx, r5)
    if d:
        return _pair(x, "W5", W5_WHAT, "disagree", "args_resolved differ (W5 takes W1's)", d, kind="args_resolved")
    if r5.get("error") or rx.get("error"):
        ca = (canonical_error(rx.get("error")) or {}).get("code") if rx.get("error") else None
        c5 = (r5.get("error") or {}).get("code") if r5.get("error") else None
        if ca != c5:
            return _pair(x, "W5", W5_WHAT, "disagree", "error %s vs W5's %s" % (ca, c5), ("/error/code", ca, c5),
                         kind="error")
        if r5.get("error"):
            return _pair(x, "W5", W5_WHAT, "agree", "the same error code %s" % c5)
    if not checked:
        return _pair(x, "W5", W5_WHAT, "agree", "W5 checks no pointer of this call")
    mine = canonical_answer(rx.get("tool"), library_part(x, rx)) if not rx.get("error") else {}
    theirs = canonical_answer(r5.get("tool"), r5.get("answer")) if isinstance(r5.get("answer"), dict) else {}
    for ptr in sorted(checked):
        f5, v5 = gp.pointer_get(theirs, ptr)
        fx, vx = gp.pointer_get(mine, ptr)
        if not f5:
            return _pair(x, "W5", W5_WHAT, "disagree", "W5 checks %s but its answer lacks it" % ptr,
                         (ptr, vx if fx else gp.MISSING, gp.MISSING), kind="oracle")
        if not fx:
            return _pair(x, "W5", W5_WHAT, "disagree", "%s (%s): absent in %s" % (ptr, checked[ptr], x),
                         (ptr, gp.MISSING, v5), kind="oracle")
        d = gp.first_difference(vx, v5)
        if d:
            return _pair(x, "W5", W5_WHAT, "disagree", "%s (%s)" % (ptr, checked[ptr]), (ptr + d[0], d[1], d[2]),
                         kind="oracle")
    return _pair(x, "W5", W5_WHAT, "agree", "%d checked pointer(s)" % len(checked))


def compare_cut(way, rec, ref, optional):
    """The cut check of one way: |rec| a cut of the same way's uncut |ref| (FORMAT-v8 §5 cut_of)."""
    if ref is None:
        return {"way": way, "outcome": "not_compared", "detail": "%s recorded no uncut reference" % way}
    if rec.get("error") or ref.get("error"):
        same = gp.first_difference(comparable_view(way, rec), comparable_view(way, ref)) is None
        return {"way": way, "outcome": "agree" if same else "disagree",
                "detail": "the same error" if same else "an error on one of the cut and the uncut call"}
    r = is_cut(canonical_answer(rec.get("tool"), library_part(way, rec)),
               canonical_answer(ref.get("tool"), library_part(way, ref)), optional)
    if r:
        return {"way": way, "outcome": "disagree", "pointer": r[0], "detail": r[1]}
    return {"way": way, "outcome": "agree", "detail": "a cut of %s's uncut answer" % way}


def compare_op(scenario, op, recs, precheck, manifest):
    """The entry of one op call (module docstring)."""
    sid, oid, tool = scenario.get("id"), op.get("id"), op.get("tool")
    ways = gp.op_ways(scenario, op)
    entry = {"tool": tool, "ways": {}, "pairs": [], "problems": []}
    if op.get("cut_of"):
        entry["cut_of"] = op["cut_of"]
    if "W3" not in ways:
        why = not_exposed_reason(scenario, op, manifest)
        entry["not_exposed"] = {"W3": why}
        entry["ways"]["W3"] = "not_exposed"
        if why is None:
            entry["problems"].append("W3 is left out without the owner's reason (not_exposed)")
    if ((precheck or {}).get("preconditions") or {}).get("outcome") == "precondition_failed":
        for w in ways:
            entry["ways"][w] = "not_compared"
        entry["outcome"] = "precondition_failed"
        entry["reason"] = precheck["preconditions"].get("detail")
        return entry
    present = {}
    for w in ways:
        r = recs.get(w, {}).get((sid, oid))
        if r is None:
            entry["ways"][w] = "missing"
        else:
            present[w] = r
            entry["ways"][w] = "selector_error" if r.get("selector_error") else (
                "error" if r.get("error") else "answered")
    if tool in SMOKE_TOOLS:
        bad = [w for w, r in present.items() if r.get("error") or r.get("selector_error")]
        entry["smoke"] = True
        entry["outcome"] = "inconclusive" if len(present) < len(ways) else ("disagree" if bad else "agree")
        if bad:
            entry["reason"] = "the smoke call failed in %s" % ", ".join(bad)
        return entry
    cut_of = op.get("cut_of")
    if cut_of:
        optional = OPTIONAL_KEYS.get(tool, frozenset())
        entry["cuts"] = [compare_cut(w, r, recs.get(w, {}).get((sid, cut_of)), optional)
                         for w, r in sorted(present.items()) if w != "W5"]
    for a, b, what in PAIRS:
        if a not in ways or b not in ways:
            continue
        if a not in present or b not in present:
            entry["pairs"].append(_pair(a, b, what, "not_compared", "%s recorded nothing" % (
                " and ".join(w for w in (a, b) if w not in present))))
            continue
        if cut_of and "W3" in (a, b):
            continue   # the service's ceiling cuts at its own place: W3 is judged by its own cut check
        entry["pairs"].append(compare_records(a, present[a], b, present[b], what))
    if "W5" in present:
        r5 = present["W5"]
        oracle = r5.get("oracle") if isinstance(r5.get("oracle"), dict) else {}
        entry["w5"] = {"checked": sorted((oracle.get("checked") or {}).keys())
                       if isinstance(oracle.get("checked"), dict) else [],
                       "consistency_only": sorted(oracle.get("consistency_only") or [])
                       if isinstance(oracle.get("consistency_only"), list) else [],
                       "inputs": oracle.get("inputs")}
        for x in ("W1", "W2", "W3", "W4"):
            if x in present and not (cut_of and x == "W3"):
                entry["pairs"].append(compare_w5(x, present[x], r5))
    disagree = [p["pair"] for p in entry["pairs"] if p["outcome"] == "disagree"]
    disagree += ["%s-uncut" % c["way"] for c in entry.get("cuts") or [] if c["outcome"] == "disagree"]
    missing = [w for w in ways if w not in present]
    all_selector = present and all(r.get("selector_error") for r in present.values())
    if disagree:
        entry["outcome"] = "disagree"
        entry["disagreeing"] = disagree
    elif missing or all_selector or entry["problems"] or any(
            c["outcome"] == "not_compared" for c in entry.get("cuts") or []):
        entry["outcome"] = "inconclusive"
        why = []
        if missing:
            why.append("no record of %s" % ", ".join(missing))
        if all_selector:
            why.append("the selector failed in every way: %s" % gp.compact(
                next(iter(present.values())).get("selector_error")))
        why.extend(entry["problems"])
        why.extend(c["detail"] for c in entry.get("cuts") or [] if c["outcome"] == "not_compared")
        entry["reason"] = "; ".join(why)
    else:
        entry["outcome"] = "agree"
    return entry


def pair_text(p):
    """One pair entry as text: what it isolates, its detail, and the first difference with both values."""
    a, _, b = p.get("pair", "-").partition("-")
    parts = [x for x in (p.get("detail"),) if x]
    if p.get("pointer"):
        parts.append("%s: %s %s, %s %s" % (p["pointer"], a, p.get(a), b, p.get(b)))
    return "%s (%s)%s" % (p.get("pair"), p.get("what"), (": " + "; ".join(parts)) if parts else "")


def disagreement_text(entry, limit=3):
    """What an op entry's disagreements are, as text: each disagreeing pair (pair_text), cut check, by-design pair or
    probe, the first |limit| of them."""
    out = [pair_text(p) for p in entry.get("pairs") or [] if p.get("outcome") == "disagree"]
    out += ["%s-uncut: %s%s" % (c["way"], c.get("detail"), (" at %s" % c["pointer"]) if c.get("pointer") else "")
            for c in entry.get("cuts") or [] if c.get("outcome") == "disagree"]
    out += ["by-design %s (%s side): %s" % (r.get("case"), r.get("side"), r.get("detail"))
            for r in entry.get("by_design") or [] if not r.get("verified")]
    if entry.get("probe") and entry.get("outcome") == "disagree" and entry.get("reason"):
        out.append("probe %s: %s" % (entry["probe"], entry["reason"]))
    return "; ".join(out[:limit]) + ("; ..." if len(out) > limit else "")


def scenario_outcome(precheck, ops):
    """fail (a pre-check failed or an op disagrees) > precondition_failed > inconclusive > pass."""
    pre = (precheck or {}).get("outcome")
    outs = [o.get("outcome") for o in ops.values() if not o.get("smoke")]
    if pre == "fail" or "disagree" in outs:
        return "fail"
    if pre == "precondition_failed" or "precondition_failed" in outs:
        return "precondition_failed"
    if pre != "pass" or "inconclusive" in outs:
        return "inconclusive"
    return "pass"


def _norm_pointer(ptr):
    return re.sub(r"/(0|[1-9][0-9]*)(?=/|$)", "/*", ptr)


def coverage(recs):
    """{tool: {"checked": {pointer: {"clause", "calls"}}, "consistency_only": {pointer: calls}}} over W5's records,
    list indices written *."""
    out = {}
    for (sid, oid), r in sorted(recs.get("W5", {}).items(), key=lambda kv: (str(kv[0][0]), str(kv[0][1]))):
        oracle = r.get("oracle") if isinstance(r.get("oracle"), dict) else {}
        t = out.setdefault(r.get("tool"), {"checked": {}, "consistency_only": {}})
        for ptr, clause in sorted((oracle.get("checked") or {}).items()) if isinstance(oracle.get("checked"),
                                                                                      dict) else []:
            c = t["checked"].setdefault(_norm_pointer(ptr), {"clause": clause, "calls": 0})
            c["calls"] += 1
            if clause not in c.setdefault("clauses", [c["clause"]]):
                c["clauses"].append(clause)
        for ptr in oracle.get("consistency_only") or [] if isinstance(oracle.get("consistency_only"), list) else []:
            p = _norm_pointer(ptr)
            t["consistency_only"][p] = t["consistency_only"].get(p, 0) + 1
    return out


def drift(recs):
    """The report-only drift lane: every drift-lib answer that differs from W4's of the same call."""
    out = {"calls": 0, "compared": 0, "differ": []}
    for key, rd in sorted(recs.get("drift", {}).items(), key=lambda kv: (str(kv[0][0]), str(kv[0][1]))):
        out["calls"] += 1
        r4 = recs.get("W4", {}).get(key)
        if r4 is None:
            continue
        out["compared"] += 1
        d = gp.first_difference(comparable_view("W4", r4), comparable_view("drift", rd))
        if d:
            out["differ"].append({"scenario": key[0], "op": key[1], "tool": rd.get("tool"), "pointer": d[0] or "/",
                                  "W4": gp.compact(d[1]), "drift": gp.compact(d[2]),
                                  "commit": (rd.get("library") or {}).get("commit")
                                  if isinstance(rd.get("library"), dict) else None})
    return out


def input_digests(run_dir):
    """{relative path: sha256 or None} of what a graphlet comparison read: the manifest and every way's records."""
    names = ["manifest.json"] + ["arms/%s/%s" % (gp.WAY_DIRS[w], gp.OPS_FILE) for w in WAYS + ("drift",)]
    return dict((n, sc.file_sha256(os.path.join(run_dir, n))) for n in names)


def compare_run(run_dir, manifest=None, checks_mod=None, extra=None):
    """graphlet_compare.json's content for a run directory (and the pre-check records it judged with); |extra|: the
    stages run after the per-op comparison and before the summary (STAGES)."""
    manifest = manifest if manifest is not None else sc.read_json(os.path.join(run_dir, "manifest.json"))
    recs, warnings, _ = gp.all_records(run_dir)
    prechecks = gp.prechecks(run_dir, manifest, recs, checks_mod)
    pre_by = dict((p["scenario"], p) for p in prechecks)
    scenarios, gaps, unplanned = {}, {}, []
    planned = set()
    for scen in gp.graphlet_scenarios(manifest):
        sid = scen.get("id")
        ops = {}
        for op in gp.scenario_ops(scen):
            planned.add((sid, op["id"]))
            ops[op["id"]] = compare_op(scen, op, recs, pre_by.get(sid), manifest)
            for w in gp.LIBRARY_WAYS:
                r = recs.get(w, {}).get((sid, op["id"]))
                for k in whitelist_gaps(op.get("tool"), (r or {}).get("answer")):
                    if k not in gaps.setdefault(op.get("tool"), []):
                        gaps[op.get("tool")].append(k)
        scenarios[sid] = {"title": scen.get("title"), "precheck": pre_by.get(sid), "ops": ops}
    for way in WAYS:
        for key in recs.get(way, {}):
            if key not in planned:
                unplanned.append({"way": way, "scenario": key[0], "op": key[1]})
    out = {"compare_version": COMPARE_VERSION, "generated_utc": sc.iso_utc(time.time()),
           "manifest": sc.manifest_identity(manifest), "inputs": input_digests(run_dir),
           "prechecks_sha256": hashlib.sha256(sc.canonical_json(prechecks)).hexdigest(),
           "whitelist_sha256": whitelist_sha256(), "warnings": warnings,
           "scenarios": scenarios, "unplanned_records": unplanned[:200],
           "summary": {"whitelist_gaps": gaps, "unplanned_records": len(unplanned)},
           "coverage": coverage(recs), "drift": drift(recs)}
    for fn in extra or ():
        fn(out, run_dir, manifest, recs, prechecks)
    summarise(out)
    return out, prechecks


def summarise(result):
    """Each scenario's outcome (scenario_outcome over its pre-check and ops), the op counts and the per-tool table
    of |result|, from its ops as they stand (after the stages)."""
    per_tool, counts = {}, dict((k, 0) for k in OP_OUTCOMES)
    for s in result["scenarios"].values():
        s["outcome"] = scenario_outcome(s.get("precheck"), s["ops"])
        for e in s["ops"].values():
            t = per_tool.setdefault(e.get("tool"), dict((k, 0) for k in ("calls",) + OP_OUTCOMES +
                                                        ("not_exposed", "smoke")))
            t["calls"] += 1
            t[e["outcome"]] = t.get(e["outcome"], 0) + 1
            t["not_exposed"] += 1 if e["ways"].get("W3") == "not_exposed" else 0
            t["smoke"] += 1 if e.get("smoke") else 0
            if not e.get("smoke"):
                counts[e["outcome"]] = counts.get(e["outcome"], 0) + 1
    result["summary"].update({"scenarios": dict((k, sum(1 for s in result["scenarios"].values()
                                                        if s["outcome"] == k)) for k in SCENARIO_OUTCOMES),
                              "ops": counts, "per_tool": per_tool})
    return result


def whitelist_sha256():
    return hashlib.sha256(sc.canonical_json(dict((k, sorted(v)) for k, v in sorted(ANSWER_KEYS.items())))).hexdigest()


# the stages the comparison runs after its own (graphlet_service_checks registers the service checks and the
# by-design pairs here): fn(result, run_dir, manifest, records, prechecks)
STAGES = []


def compute(run_dir, manifest=None, write=True, out_path=None):
    """Compares a run's area G (compare_run with STAGES) and writes graphlet_compare.json and the pre-check records
    (unless not |write|); returns the comparison."""
    _load_stages()
    result, prechecks = compare_run(run_dir, manifest, extra=STAGES)
    if write:
        gp.write_prechecks(run_dir, prechecks)
        sc.write_json_atomic(out_path or os.path.join(run_dir, COMPARE_FILE), result)
    return result


def stale(run_dir, result):
    """The inputs whose content differs from what |result| read ([] when current)."""
    now = input_digests(run_dir)
    was = (result or {}).get("inputs") or {}
    return [k for k in sorted(now) if was.get(k, "unrecorded") != now[k]]


def load_or_compute(run_dir, manifest=None):
    """RUN_DIR/graphlet_compare.json when it read the run's current inputs, else computed and written (the gate's
    and the report's reading, as compare.load_or_compute)."""
    data = sc.read_json(os.path.join(run_dir, COMPARE_FILE))
    if isinstance(data, dict) and not stale(run_dir, data):
        return data
    return compute(run_dir, manifest)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("run_dir")
    ap.add_argument("--out", default=None, help="default RUN_DIR/graphlet_compare.json")
    a = ap.parse_args(argv)
    manifest = sc.read_json(os.path.join(a.run_dir, "manifest.json"))
    if not isinstance(manifest, dict):
        sys.stderr.write("compare_graphlets.py: %s has no manifest.json\n" % a.run_dir)
        return 2
    result = compute(a.run_dir, manifest, out_path=a.out)
    s = result["summary"]
    sys.stdout.write("graphlet compare %s: scenarios %s; ops %s\n" % (
        a.run_dir, ", ".join("%s %d" % kv for kv in sorted(s["scenarios"].items()) if kv[1]) or "none",
        ", ".join("%s %d" % kv for kv in sorted(s["ops"].items()) if kv[1]) or "none"))
    return 0


def _load_stages():
    """Imports the modules that register comparison stages (graphlet_service_checks) when present."""
    try:
        import graphlet_service_checks  # noqa: F401
    except ImportError:
        pass


if __name__ == "__main__":
    sys.exit(main())
