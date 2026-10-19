#!/usr/bin/env python3
"""The five-way compare of area G (manifest v8; DESIGN-v8 §2, §5, FORMAT-v8 §1-§5): per op call of every graphlet
scenario, the answers of the ways that run it, after the scenario's pre-checks (graphlet_precheck.py).

    compare_graphlets.py RUN_DIR [--out FILE]     (default RUN_DIR/graphlet_compare.json; the pre-check records go to
                                                   RUN_DIR/arms/precheck/graphlet_prechecks.jsonl)

Per op call, in this order:
1. the ways it runs in (its `ways`, else the scenario's ways.default, else all five); W3 left out is `not_exposed`
   with the owner's reason (FORMAT v2 D2: the op's not_exposed, else the scenario's graphlet.not_exposed by tool or
   tool.argument); W3 left out without a reason is a manifest problem (inconclusive);
2. a way that recorded nothing is missing; a way whose runner failed (FORMAT v2 B2: no_document, exception, ...)
   did not run the call and is not compared (inconclusive); a selector that failed in a way did not run the call
   there (the same failure in every way: inconclusive, the selector named; in some ways only: a disagreement);
3. args_resolved equal across the ways (access tokens and cursors left out): a selector that resolved differently is
   a disagreement of the pair, named as such, and that pair's answers are not compared;
4. each answer canonicalised (canonical_answer): W3's service fields (the keys outside ANSWER_KEYS, the library's
   per-tool answer keys, frozen here) stripped and kept for the service checks; handles to symbols (<handle>,
   <handle_b>, <handle:N>), cursors removed, parsed JSON values only (numbers by value); the `local` block compared
   only between library ways of one (commit, interpreter, work model) (D-R10); an export's `path` and `bytes` left
   out (the export directory, the handle and the envelope's timing make them differ by construction): the file
   itself is compared by its content (FORMAT v2 B3, file_view), a file not kept leaves the pair not compared;
5. whole assembled answers compared by pairs: W1-W2 (transport), W1-W4 (fetch path), W3-W4 (service wrapper), and
   each of them against W5 (library vs independent reading) by graphlet_oracle.compare (lane O's one comparison of an
   answer against a W5 record, FORMAT v2 F2: checked pointers, unordered lists, absent keys, expected_error, cut
   rows); W5's other pointers are consistency-only and listed; without lane O's module a W5 pair is not compared
   and the op is inconclusive;
6. errors by code and structured fields (never message or hint text; result_too_large's `bytes` depends on the
   handle's length); a fitting error (FIT_ERRORS) is not comparable against W5 (FORMAT v2 F3), compared as usual
   between W1..W4;
   on a pair with W3 a row marked _cut_fields (FORMAT v2 A3), on either side, is judged by lane O's cut rule on the
   fields it names only (cut_aware_difference: a string's leading part, a list's leading whole items, never empty for
   a non-empty value; both rows cut: their common leading part), every other field equal; W3's longer handle moves
   the cut; library ways compare strictly;
7. a deliberate cut (the op's `cut_of`: the op id of the same call uncut) is judged in each way against that way's
   uncut answer (is_cut): it must show a cut (a marker, paging_complete false, fewer rows than total; else "not
   cut", a failure); a key it lacks optional (OPTIONAL_KEYS) or named in its fields_cut; a key it adds a cut marker;
   a field a row names cut by lane O's rule; the page's rows and the tool's optional lists leading parts; every other
   value equal; the library ways' cut answers must still be equal (one handle length, one ceiling).
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

try:   # lane O's W5 module: its compare() is the one comparison of an answer against a W5 record (FORMAT v2 F2)
    import graphlet_oracle  # noqa: E402
except ImportError:   # pragma: no cover - before lane O's commit lands: W5 pairs are not compared, named
    graphlet_oracle = None

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
# codes compared without their fields: a failed seed's no_graphlet is the runner's on W1/W2/W4 and the service's on
# W3 (FORMAT v2 B2), each with fields of its own
ERROR_CODE_ONLY = ("no_graphlet",)
# the runners' own failures (FORMAT v2 B2): the way did not run the call (no document, a refused document, an
# exception, an unknown tool, a runaway paging loop, an assembly error): never compared, the op inconclusive
RUNNER_ERRORS = ("no_document", "document_refused", "exception", "unknown_tool", "paging_runaway", "assembly_error")
# per tool, answer keys whose value depends on where the way ran, left out of the comparison
ENV_KEYS = {"graphlet_export": ("path", "bytes"), "graphlet_save": ("path", "bytes")}
# the tools that make a file (FORMAT v2 B3: the record's `file` {name, bytes, sha256, raw}, the text kept at
# arms/<way dir>/raw/<raw>): their files are compared by content (file_view), not only their receipts
FILE_TOOLS = ("graphlet_export", "graphlet_save")
# the per-request values of an export's json (to_json: the seed summary's timing, and the annotation's lookahead
# counters, which differ between two engine runs of one request: FORMAT v2 E1)
JSON_EXPORT_DROP = ("timing",)
ANNOTATION_LOOKAHEAD = ("direct_reads", "keys_mapped")
# a cut answer may add these keys (the ceiling's markers); a cut answer may lack these keys (per tool) and the keys
# its own fields_cut / cut_fields name
CUT_MARKERS = frozenset(("fields_cut", "row_truncated", "cut_fields", "caveat_effects_cut", "caveats_truncated",
                         "next_cursor", "truncated", "next_from", "_cut_fields"))
# FORMAT v2 A3: a row the ceiling cut carries "_cut_fields" (its page's cut_fields) in the assembled rows; the page's
# row_truncated / cut_fields stay merged at the top
ROW_CUT = "_cut_fields"
PAGE_CUT_MARKERS = ("row_truncated", "cut_fields")
OPTIONAL_KEYS = {
    "graphlet_summary": frozenset(("effect", "message", "top_labels", "caveats")),
    "graphlet_compare": frozenset(("evidence", "export")),
}
SMOKE_TOOLS = ("graphlet_list", "graphlet_free")
# the pairs whose agreement makes W1 = W2 = W3 = W4 (a spanning tree), each with what it isolates
PAIRS = (("W1", "W2", "transport"), ("W1", "W4", "fetch path"), ("W3", "W4", "service wrapper"))
W5_WHAT = "library vs independent reading"
# a pair not compared for one of these makes the op inconclusive (a missing way is named by its own reason)
UNREAD_KINDS = ("oracle_unavailable", "oracle_error", "no_file", "runner_error")
# errors of fitting an answer to a ceiling (FORMAT v2 F3, VO5): not comparable against W5, which makes no claim about
# fitting (INVENTORY §1.1: result_too_large, receipt_too_large; on W3 also op_busy, the service's worker state);
# between W1..W4 they are compared as usual (W3 calls are designed to fit, DESIGN-v8 §5)
FIT_ERRORS = ("result_too_large", "receipt_too_large")
FIT_ERRORS_W3 = FIT_ERRORS + ("op_busy",)
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
            # the local block is the answer's (or an error's) top-level key; evidence.limitations.local below
            if k in CURSOR_KEYS or (drop_local and top and k == "local") or (top and k in ENV_KEYS.get(tool, ())):
                continue
            if top and tool == "graphlet_compare" and k in TOP_HANDLE_KEYS and isinstance(v, str):
                out[k] = "<handle>" if k == "a" else "<handle_b>"
            elif (k in HANDLE_KEYS or (top and k in TOP_HANDLE_KEYS)) and isinstance(v, str) and is_symbol(v):
                out[k] = v   # a runner's symbol already (FORMAT v2 B1: <handle>, <handle:S>, <handle:OP>)
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
    if code in ERROR_CODE_ONLY:
        fields = {}
    skip = set(ERROR_TEXT_FIELDS) | set(ERROR_ENV_FIELDS.get(code, ()))
    fields = dict((k, v) for k, v in fields.items() if k not in skip)
    return {"code": code, "fields": _canon(fields, {}, not keep_local, top=True)}


def canonical_args(args):
    """args_resolved as the ways compare it: without per-way credentials and cursors (ARG_SKIP)."""
    if not isinstance(args, dict):
        return args
    return dict((k, v) for k, v in args.items() if k not in ARG_SKIP)


def file_path(run_dir, way, rec):
    """Where a record's file text is kept (FORMAT v2 B3: arms/<way dir>/raw/<file.raw>), or None."""
    f = rec.get("file") if isinstance(rec.get("file"), dict) else None
    raw = f.get("raw") if f else None
    if not isinstance(raw, str) or not raw or os.path.isabs(raw) or ".." in raw.split("/"):
        return None
    return os.path.join(run_dir, "arms", gp.WAY_DIRS[way], "raw", raw)


def file_view(rec, text):
    """A file's content as the ways compare it (INVENTORY §4.2): fasta and gfa byte for byte (their sha256); json
    parsed, a whole export without its seed summary's timing and the annotation's lookahead counters, a walk or
    compare export with its handles as symbols; mgt (an mgt export, a save) by the body digest (graphlet_doc: H's
    seed_index masked, J and Z left out, so a stored body and the library's dump with its envelope have one)."""
    args = rec.get("args_resolved") if isinstance(rec.get("args_resolved"), dict) else {}
    fmt = "mgt" if rec.get("tool") == "graphlet_save" else args.get("format")
    what = args.get("what")
    try:
        if fmt in ("fasta", "gfa"):
            return {"format": fmt, "sha256": hashlib.sha256(text.encode("utf-8")).hexdigest()}
        if fmt == "json":
            doc = json.loads(text)
            if what == "compare":
                return {"format": "json", "what": what, "json": canonical_answer("graphlet_compare", doc)}
            if what == "walk":
                return {"format": "json", "what": what, "json": canonical_answer("graphlet_walk", doc)}
            if isinstance(doc, dict):
                doc = dict((k, v) for k, v in doc.items() if k not in JSON_EXPORT_DROP)
                if isinstance(doc.get("annotation"), dict):
                    doc["annotation"] = dict((k, v) for k, v in doc["annotation"].items()
                                             if k not in ANNOTATION_LOOKAHEAD)
            return {"format": "json", "json": doc}
        if fmt == "mgt":
            digest = gp.graphlet_doc.body_digest(text) if gp.graphlet_doc is not None else \
                hashlib.sha256("\n".join(gp.body_records(text)).encode("utf-8")).hexdigest()
            return {"format": "mgt", "body_digest": digest}
    except (ValueError, TypeError) as e:
        return {"unreadable": "%s: %s" % (type(e).__name__, str(e)[:120])}
    return {"format": fmt, "sha256": hashlib.sha256(text.encode("utf-8")).hexdigest()}


def attach_files(run_dir, recs):
    """Reads every recorded file (FORMAT v2 B3) and keeps on its record `_file_text` and `_file_view` (None for a
    record of no file; {"not_kept": path} for a file whose text is not where the record says)."""
    for way, by in recs.items():
        if way not in gp.WAY_DIRS:
            continue
        for rec in by.values():
            if rec.get("tool") not in FILE_TOOLS or not isinstance(rec.get("file"), dict):
                continue
            path = file_path(run_dir, way, rec)
            text = gp._read_text(path) if path and os.path.exists(path) else None
            rec["_file_text"] = text
            rec["_file_view"] = file_view(rec, text) if text is not None else {
                "not_kept": "arms/%s/raw/%s" % (gp.WAY_DIRS[way], (rec["file"] or {}).get("raw"))}


def same_library_env(ra, rb):
    """Whether two library records ran one (commit, interpreter, work model) with one configuration (FORMAT v2 B4:
    library.limits, the GraphletTools configuration): only then is `local` compared (D-R10)."""
    la, lb = ra.get("library") or {}, rb.get("library") or {}
    if not isinstance(la, dict) or not isinstance(lb, dict):
        return False
    keys = ("commit", "python", "work_model")
    return all(la.get(k) is not None for k in keys) and all(la.get(k) == lb.get(k) for k in keys + ("limits",))


def is_symbol(value):
    """A runner's handle symbol (FORMAT v2 B1): <handle>, <handle:S>, <handle:OP>."""
    return isinstance(value, str) and value.startswith("<handle") and value.endswith(">")


def runner_error(rec):
    """The runner-side code of a record whose way did not run the call (RUNNER_ERRORS), "no_answer" for a record
    with neither an answer nor an error nor a selector error (a runner bug or a lost call, VERIFY-C VC6), or None."""
    err = rec.get("error") if isinstance(rec.get("error"), dict) else None
    if err and err.get("code") in RUNNER_ERRORS:
        return err.get("code")
    if not rec.get("error") and not rec.get("selector_error") and rec.get("answer") is None \
            and not isinstance(rec.get("oracle"), dict):   # a W5 record may declare instead of answering
        return "no_answer"
    return None


def comparable_view(way, rec, keep_local=False):
    """{"error": canonical error} or {"answer": canonical answer[, "file": file_view]} of a record."""
    if rec.get("error"):
        return {"error": canonical_error(rec["error"], keep_local)}
    out = {"answer": canonical_answer(rec.get("tool"), library_part(way, rec), keep_local)}
    if rec.get("_file_view") is not None:
        out["file"] = rec["_file_view"]
    return out


def _file_unread(way, rec):
    """Why a record of a file-making tool cannot be compared by its file (no `file`, its text not kept or not
    readable), or None."""
    if rec.get("tool") not in FILE_TOOLS or rec.get("error") or rec.get("selector_error"):
        return None
    view = rec.get("_file_view")
    if view is None:
        return "%s recorded no file (FORMAT v2 B3)" % way
    if "not_kept" in view:
        return "%s's file is not kept at %s" % (way, view["not_kept"])
    if "unreadable" in view:
        return "%s's file does not read (%s)" % (way, view["unreadable"])
    return None


# ---------------------------------------------------------------------------------------------------------------
# a cut of an uncut answer, in the suite's words

# a page's own list: a cut page holds a leading part of the whole list's rows
PAGE_LIST = "rows"
# the evidence of a cut (VERIFY-C VC4): one of these keys anywhere in the answer (next_cursor is removed by the
# runners' assembly, so it is no evidence)
CUT_EVIDENCE = frozenset(("fields_cut", "row_truncated", "cut_fields", "_cut_fields", "caveat_effects_cut",
                          "caveats_truncated", "truncated", "next_from"))


def is_cut(cut, uncut, optional=frozenset(), pointer="", fn=None):
    """None when |cut| is a cut of |uncut|, else (pointer, why) (DESIGN-v8 §5, REVIEW D-R8, with VERIFY-C VC4's
    restriction to what the answer names cut):
    - a key the cut lacks: one of |optional| (the tool's optional keys, OPTIONAL_KEYS) or named in this object's
      fields_cut;
    - a key only the cut has: a cut marker (CUT_MARKERS);
    - a field this object names cut (a row's _cut_fields, FORMAT v2 A3): lane O's rule (field_cut_problem: a
      string's leading part, a list's leading whole items, never empty for a non-empty value);
    - a page's rows and the tool's optional lists (summary top_labels, caveats): a leading part, each item a cut of
      the uncut item (so a caveat may lack its optional effect);
    - every other list: the same length, item by item; every other value, strings included, equal.
    Without lane O's rule (|fn| None and graphlet_oracle absent) a named field is uncomparable."""
    fn = fn or oracle_is_cut_of()
    if isinstance(cut, dict) and isinstance(uncut, dict):
        dropped = set(x for x in (cut.get("fields_cut") or []) if isinstance(x, str))
        shortened = set(x for x in (cut.get(ROW_CUT) or []) if isinstance(x, str))
        for k in sorted(set(cut) | set(uncut), key=str):
            p = pointer + "/" + gp.pointer_escape(k)
            if k not in uncut:
                if k in CUT_MARKERS:
                    continue
                return p, "a key only the cut answer has, not a cut marker"
            if k not in cut:
                if k in optional or k in dropped:
                    continue
                return p, "a key the cut answer lacks, neither optional nor named in its fields_cut"
            if k in shortened:
                if fn is None:
                    return p, "uncomparable: graphlet_oracle's cut rule (lane O) is not available"
                why = field_cut_problem(uncut[k], cut[k], fn)
                if why:
                    return p, "a field the row names cut: %s" % why
                continue
            if isinstance(cut[k], list) and isinstance(uncut[k], list) and (k == PAGE_LIST or k in optional):
                if len(cut[k]) > len(uncut[k]):
                    return p, "a list longer than the uncut one (%d > %d)" % (len(cut[k]), len(uncut[k]))
                if k == PAGE_LIST and uncut[k] and not cut[k]:
                    # a page holds at least its first row (a row larger than the page comes alone, cut: _page)
                    return p, "a page without a row where the uncut list has %d" % len(uncut[k])
                for i, (x, y) in enumerate(zip(cut[k], uncut[k])):
                    r = is_cut(x, y, optional, "%s/%d" % (p, i), fn)
                    if r:
                        return r
                continue
            r = is_cut(cut[k], uncut[k], optional, p, fn)
            if r:
                return r
        return None
    if isinstance(cut, list) and isinstance(uncut, list):
        if len(cut) != len(uncut):
            return pointer or "/", "a list of another length (%d vs %d) that no cut names" % (len(cut), len(uncut))
        for i, (x, y) in enumerate(zip(cut, uncut)):
            r = is_cut(x, y, optional, "%s/%d" % (pointer, i), fn)
            if r:
                return r
        return None
    if isinstance(cut, (dict, list)) or isinstance(uncut, (dict, list)):
        return pointer or "/", "a value of another type than the uncut one"
    return None if gp.leaf_equal(cut, uncut) else (pointer or "/", "a value other than the uncut one (%s vs %s)" % (
        gp.compact(cut, 60), gp.compact(uncut, 60)))


def cut_evidence(rec, answer):
    """What shows that a call was cut (VERIFY-C VC4), or None: a cut marker anywhere in its answer (CUT_EVIDENCE),
    paging_complete false (FORMAT v2 A4), or fewer rows than its total."""
    if rec.get("paging_complete") is False:
        return "paging_complete false"
    found = []

    def walk(o):
        if isinstance(o, dict):
            found.extend(k for k in o if k in CUT_EVIDENCE and o[k] not in (None, False, [], ""))
            for v in o.values():
                walk(v)
        elif isinstance(o, list):
            for v in o:
                walk(v)
    walk(answer)
    if found:
        return "the marker %s" % sorted(set(found))[0]
    if isinstance(answer, dict) and isinstance(answer.get("rows"), list) and isinstance(answer.get("total"), int) \
            and len(answer["rows"]) < answer["total"]:
        return "%d of %d rows" % (len(answer["rows"]), answer["total"])
    return None


def _rows_cut(ans):
    return isinstance(ans, dict) and isinstance(ans.get("rows"), list) and any(
        isinstance(r, dict) and r.get(ROW_CUT) for r in ans["rows"])


def oracle_is_cut_of():
    """Lane O's cut rule as fn(whole, cut) -> why |cut| is no cut of |whole| ('' when it is):
    graphlet_oracle.cut_problem(cut_value, full_value) (the suite's one cut rule, public since lane O's round 3), its
    private _is_cut_of before; None without the module. On a field's value: a string's leading part, a list's leading
    whole items (VERIFY-C VC1)."""
    if graphlet_oracle is None:
        return None
    public = getattr(graphlet_oracle, "cut_problem", None)
    if callable(public):
        return lambda whole, cut: public(cut, whole) or ""
    private = getattr(graphlet_oracle, "_is_cut_of", None)
    if callable(private):
        return lambda whole, cut: "" if private(whole, cut, "ordered") else "not a cut of the full value (ordered)"
    return None


def field_cut_problem(whole, cut, fn):
    """Why |cut| is no cut of |whole| for a field a row names cut ('' when it is): lane O's rule (|fn|: a string's
    leading part, a list's leading whole items, never a nested prefix inside an item), and never an empty value for a
    non-empty one (the coordinator's VC1 decision: empty-for-full is no cut)."""
    if isinstance(whole, (str, list)) and isinstance(cut, (str, list)) and whole and not cut:
        return "empty for a non-empty value"
    why = fn(whole, cut)
    return "not a leading part: %s" % why if why else ""


def row_cut_difference(x, y, pointer, fn):
    """The first difference of two rows of which at least one carries _cut_fields (FORMAT v2 A3), or None: a field
    named cut on one side only must be a cut of the other side's (field_cut_problem); a field named on both sides is
    compared on their common leading part (the shorter a cut of the longer); every field named on neither side
    equal; each side's own _cut_fields marker left out. -> (pointer, x's value, y's value, why) or None."""
    nx = set(x.get(ROW_CUT) or []) if isinstance(x, dict) else set()
    ny = set(y.get(ROW_CUT) or []) if isinstance(y, dict) else set()
    if not isinstance(x, dict) or not isinstance(y, dict):
        d = gp.first_difference(x, y, pointer)
        return d and d + ("",)
    for k in sorted((set(x) | set(y)) - {ROW_CUT}, key=str):
        p = pointer + "/" + gp.pointer_escape(k)
        vx, vy = x.get(k, gp.MISSING), y.get(k, gp.MISSING)
        if vx is gp.MISSING or vy is gp.MISSING:
            return p, vx, vy, "a field on one side only"
        if k in nx and k in ny:
            if isinstance(vx, (str, list)) and type(vx) is type(vy):
                longer, shorter = (vx, vy) if len(vx) >= len(vy) else (vy, vx)
                why = field_cut_problem(longer, shorter, fn)
            else:
                why = "a cut field that is no string or list"
            if why:
                return p, vx, vy, "both rows cut %s: %s" % (k, why)
        elif k in nx or k in ny:
            whole, cut = (vy, vx) if k in nx else (vx, vy)
            why = field_cut_problem(whole, cut, fn)
            if why:
                return p, vx, vy, "the row %s cut %s: %s" % ("a" if k in nx else "b", k, why)
        else:
            d = gp.first_difference(vx, vy, p)
            if d:
                return d + ("a field no row names cut",)
    return None


def cut_aware_difference(a, b):
    """The first difference of two ways' canonical answers when a row of either is marked cut (FORMAT v2 A3), else
    first_difference's: (pointer, value in a, value in b[, why]). The page markers (row_truncated, cut_fields) are
    left to the rows; rows are compared one by one, a row carrying _cut_fields by row_cut_difference: only the fields
    it names may be a cut (lane O's rule, in either direction), every other field equal. For pairs whose ceilings or
    handle lengths differ (W3 against a library way: INVENTORY §1.2, the cut moves with the handle's length), never
    between library ways, which cut identically. Without lane O's rule a cut row is not compared: ("", ..., why)
    with why starting "uncomparable"."""
    if not (_rows_cut(a) or _rows_cut(b)):
        return gp.first_difference(a, b)
    fn = oracle_is_cut_of()
    if fn is None:
        return "/rows", gp.MISSING, gp.MISSING, "uncomparable: graphlet_oracle's cut rule (lane O) is not available"
    rest_a = dict((k, v) for k, v in a.items() if k not in PAGE_CUT_MARKERS + ("rows",))
    rest_b = dict((k, v) for k, v in b.items() if k not in PAGE_CUT_MARKERS + ("rows",))
    d = gp.first_difference(rest_a, rest_b)
    if d:
        return d
    ra, rb = a.get("rows") or [], b.get("rows") or []
    if len(ra) != len(rb):
        i = min(len(ra), len(rb))
        return "/rows/%d" % i, ra[i] if i < len(ra) else gp.MISSING, rb[i] if i < len(rb) else gp.MISSING
    for i, (x, y) in enumerate(zip(ra, rb)):
        if (isinstance(x, dict) and x.get(ROW_CUT)) or (isinstance(y, dict) and y.get(ROW_CUT)):
            d = row_cut_difference(x, y, "/rows/%d" % i, fn)
        else:
            d = gp.first_difference(x, y, "/rows/%d" % i)
        if d:
            return d
    return None


# ---------------------------------------------------------------------------------------------------------------
# the manifest's op plan

def scenario_by_id(manifest):
    return dict((s.get("id"), s) for s in gp.graphlet_scenarios(manifest))


def not_exposed_reason(scenario, op, manifest=None):
    """The owner's reason W3 does not run |op| (FORMAT v2 D2: the op's not_exposed, a text or {"W3": reason}), else
    the scenario's graphlet.not_exposed keyed by the tool or by tool.argument for an argument the op sends (both
    inside content_sha256); None without one."""
    ne = op.get("not_exposed")
    if isinstance(ne, dict) and isinstance(ne.get("W3"), str) and ne["W3"].strip():
        return ne["W3"].strip()
    if isinstance(ne, str) and ne.strip():
        return ne.strip()
    tool = op.get("tool")
    keys = [tool] + ["%s.%s" % (tool, a) for a in sorted(op.get("args") or {})]
    for table in ((scenario.get("graphlet") or {}).get("not_exposed"),):
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
    ea, eb = runner_error(ra), runner_error(rb)
    if ea or eb:
        return _pair(a, b, what, "not_compared", "%s did not run the call (runner: %s)" % (
            a if ea else b, ea or eb), kind="runner_error")
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
    unread = _file_unread(a, ra) or _file_unread(b, rb)
    if unread:
        return _pair(a, b, what, "not_compared", unread, kind="no_file")
    keep = a in gp.LIBRARY_WAYS and b in gp.LIBRARY_WAYS and same_library_env(ra, rb)
    va, vb = comparable_view(a, ra, keep), comparable_view(b, rb, keep)
    if ("error" in va) != ("error" in vb):
        err, ok = (a, b) if "error" in va else (b, a)
        code = (va.get("error") or vb.get("error") or {}).get("code")
        return _pair(a, b, what, "disagree", "%s answered the error %s, %s an answer" % (err, code, ok),
                     ("/", va, vb), kind="error", views=(va, vb))
    if "answer" in va and "W3" in (a, b) and isinstance(va["answer"], dict) and isinstance(vb["answer"], dict):
        d = cut_aware_difference(va["answer"], vb["answer"])   # FORMAT v2 A3: the service's cut rows
        d = d and ("/answer" + (d[0] or ""),) + tuple(d[1:])
        if not d:   # the rest of the view (a file's content) equal as usual
            d = gp.first_difference(dict((k, v) for k, v in va.items() if k != "answer"),
                                    dict((k, v) for k, v in vb.items() if k != "answer"))
    else:
        d = gp.first_difference(va, vb)
    if d and len(d) > 3 and str(d[3]).startswith("uncomparable"):
        return _pair(a, b, what, "not_compared", d[3], kind="oracle_unavailable")
    if d:
        return _pair(a, b, what, "disagree", (d[3] or None) if len(d) > 3 else None, d[:3],
                     kind="error" if "error" in va else "answer", views=(va, vb))
    return _pair(a, b, what, "agree")


def oracle_compare_fn():
    """graphlet_oracle.compare (lane O's module; FORMAT v2 F2: the one comparison of an answer against a W5 record), or
    None when the module or the function is absent."""
    fn = getattr(graphlet_oracle, "compare", None) if graphlet_oracle is not None else None
    return fn if callable(fn) else None


def oracle_view(way, rec):
    """What a way's call gives the W5 comparison: {"error": its canonical error} for an error (a W5 record declaring
    expected_error requires one), else its assembled answer without the service fields, canonical (compare_records'
    reading)."""
    if rec.get("error"):
        return {"error": canonical_error(rec["error"])}
    if rec.get("tool") in FILE_TOOLS and callable(getattr(graphlet_oracle, "export_view", None)):
        # an export's receipt and file in the shape W5 answers (lane O's export_view, FORMAT v2 B3)
        args = rec.get("args_resolved") or {}
        fmt = "save" if rec.get("tool") == "graphlet_save" else args.get("format")
        return graphlet_oracle.export_view(fmt, library_part(way, rec), rec.get("_file_text"), args.get("what"))
    return canonical_answer(rec.get("tool"), library_part(way, rec))


# graphlet_oracle.compare's outcomes (lane O); a bare list of differences (None) is read as agree / disagree;
# inconclusive (W5 could not judge a call it should) blocks
ORACLE_OUTCOMES = ("agree", "disagree", "not_comparable", "not_checked", "inconclusive")


def _oracle_result(res):
    """(outcome or None, [difference], detail or None) of what graphlet_oracle.compare returned: lane O's
    {outcome: agree | disagree | not_comparable | not_checked, detail, checked, differences}, or a bare list of
    differences ({pointer, clause, oracle, other | missing | present})."""
    if isinstance(res, dict):
        diffs = res.get("differences") if isinstance(res.get("differences"), list) else []
        return res.get("outcome"), diffs, res.get("detail")
    return None, list(res or []), None


def fit_errors(way):
    """The fitting errors W5 makes no claim about (FORMAT v2 F3): lane O's graphlet_oracle.FIT_ERRORS (the owner's
    list) where its module loads, else result_too_large and receipt_too_large (op_busy too on W3)."""
    owner = getattr(graphlet_oracle, "FIT_ERRORS", None) if graphlet_oracle is not None else None
    if isinstance(owner, (tuple, list, set, frozenset)) and owner:
        return tuple(owner)
    return FIT_ERRORS_W3 if way == "W3" else FIT_ERRORS


def compare_w5(x, rx, r5):
    """The pair entry of way |x| against W5: graphlet_oracle.compare(W5's record, x's view) decides (FORMAT v2 F2;
    lane O's rules: checked pointers, unordered and position-ordered lists, absent keys, expected_error, cut rows,
    fitting errors); this stage adds only the selector, runner and args_resolved readings and names the first
    difference with both values. Its outcomes pass through: agree, disagree, not_comparable (a fitting error),
    not_checked (W5 makes no claim on this call)."""
    if rx.get("selector_error") or r5.get("selector_error"):
        return _pair(x, "W5", W5_WHAT, "not_compared", "the selector failed in %s" % (
            x if rx.get("selector_error") else "W5"), kind="selector")
    if runner_error(rx) or runner_error(r5):
        return _pair(x, "W5", W5_WHAT, "not_compared", "%s did not run the call (runner: %s)" % (
            x if runner_error(rx) else "W5", runner_error(rx) or runner_error(r5)), kind="runner_error")
    d = _args_difference(rx, r5)
    if d:
        return _pair(x, "W5", W5_WHAT, "disagree", "args_resolved differ (W5 takes W1's)", d, kind="args_resolved")
    code = (rx.get("error") or {}).get("code") if isinstance(rx.get("error"), dict) else None
    if code in fit_errors(x):
        return _pair(x, "W5", W5_WHAT, "not_comparable", "%s answered %s: W5 makes no claim about fitting an answer "
                                                         "to a ceiling (FORMAT v2 F3)" % (x, code), kind="fit")
    fn = oracle_compare_fn()
    if fn is None:
        return _pair(x, "W5", W5_WHAT, "not_compared", "graphlet_oracle.compare (lane O, FORMAT v2 F2) is not "
                                                       "available", kind="oracle_unavailable")
    try:
        outcome, diffs, detail = _oracle_result(fn(r5, oracle_view(x, rx)))
    except Exception as e:   # a W5 record the comparison cannot read: nothing compared, named
        return _pair(x, "W5", W5_WHAT, "not_compared", "graphlet_oracle.compare raised %s: %s" % (
            type(e).__name__, str(e)[:200]), kind="oracle_error")
    if outcome not in (None,) + ORACLE_OUTCOMES:
        # an outcome outside lane O's vocabulary is never an agreement (VERIFY-C VC8)
        return _pair(x, "W5", W5_WHAT, "not_compared", "graphlet_oracle.compare answered the outcome %r%s" % (
            outcome, (": %s" % detail) if detail else ""), kind="oracle_error")
    if outcome == "inconclusive":
        return _pair(x, "W5", W5_WHAT, "not_compared", "W5 could not judge the call (graphlet_oracle.compare): %s"
                     % (detail or "inconclusive"), kind="oracle_error")
    if outcome in ("not_comparable", "not_checked"):
        return _pair(x, "W5", W5_WHAT, outcome, detail or "graphlet_oracle.compare: %s" % outcome, kind="oracle")
    if diffs or outcome == "disagree":
        first = diffs[0] if diffs and isinstance(diffs[0], dict) else {"pointer": "/", "clause": detail}
        theirs = first["other"] if "other" in first else first["present"] if "present" in first else gp.MISSING
        # an error where W5 answers, or an answer where W5 expects an error, is a difference of kind error
        kind = "error" if (rx.get("error") or str(first.get("pointer") or "").startswith("/error")) else "oracle"
        return _pair(x, "W5", W5_WHAT, "disagree", "%d difference(s); %s" % (len(diffs), first.get("clause")),
                     (first.get("pointer") or "/", theirs, first.get("oracle", gp.MISSING)), kind=kind)
    oracle = r5.get("oracle") if isinstance(r5.get("oracle"), dict) else {}
    declared = oracle.get("declared")
    n = len(oracle.get("checked") or {}) if isinstance(oracle.get("checked"), dict) else 0
    return _pair(x, "W5", W5_WHAT, "agree", detail or ("W5 declares %s" % declared if declared else
                                                        "%d checked pointer(s)" % n))


def compare_cut(way, rec, ref, optional):
    """The cut check of one way: |rec| a cut of the same way's uncut |ref| (cut_of), which must show a cut
    (cut_evidence: the scenario was designed to cut; without one it is "not cut", a failure, VERIFY-C VC4)."""
    if ref is None:
        return {"way": way, "outcome": "not_compared", "detail": "%s recorded no uncut reference" % way}
    if rec.get("error") or ref.get("error"):
        same = gp.first_difference(comparable_view(way, rec), comparable_view(way, ref)) is None
        return {"way": way, "outcome": "agree" if same else "disagree",
                "detail": "the same error" if same else "an error on one of the cut and the uncut call"}
    mine = canonical_answer(rec.get("tool"), library_part(way, rec))
    evidence = cut_evidence(rec, mine)
    if evidence is None:
        return {"way": way, "outcome": "disagree", "detail": "not cut: the call was designed to cut (cut_of) and "
                                                             "its answer shows no cut (no marker, its paging "
                                                             "complete, every row)"}
    r = is_cut(mine, canonical_answer(ref.get("tool"), library_part(way, ref)), optional)
    if r and r[1].startswith("uncomparable"):
        return {"way": way, "outcome": "not_compared", "detail": r[1]}
    if r:
        return {"way": way, "outcome": "disagree", "pointer": r[0], "detail": r[1]}
    return {"way": way, "outcome": "agree", "detail": "a cut of %s's uncut answer (%s)" % (way, evidence)}


# the order in which a way left unconnected by PAIRS is joined to the others (VERIFY-C VC7)
ANCHOR_ORDER = ("W1", "W4", "W2", "W3")


def op_pairs(ways):
    """The pairs an op compares (W5 aside): PAIRS among its ways, then, for every way those leave unconnected (an op
    with W3 but no W4, W2 or W4 without W1), a pair with the first way of ANCHOR_ORDER outside its group, so every
    way the op runs in is compared with another (VERIFY-C VC7); a pair with W3 is "service wrapper", else "library
    ways"."""
    listed = [w for w in ANCHOR_ORDER if w in ways]
    group = dict((w, w) for w in listed)

    def find(w):
        while group[w] != w:
            w = group[w]
        return w
    out = []
    for a, b, what in PAIRS:
        if a in group and b in group:
            out.append((a, b, what))
            group[find(b)] = find(a)
    for w in listed[1:]:
        if find(w) != find(listed[0]):
            anchor = next(x for x in listed if find(x) != find(w))
            a, b = sorted((anchor, w), key=lambda x: WAYS.index(x))
            out.append((a, b, "service wrapper (no W4)" if "W3" in (a, b) else "library ways"))
            group[find(w)] = find(anchor)
    return out


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
            entry["ways"][w] = "selector_error" if r.get("selector_error") else "runner_error" if runner_error(r) \
                else "error" if r.get("error") else "answered"
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
    for a, b, what in op_pairs(ways):
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
        for key in ("declared", "unordered", "absent", "note"):
            if oracle.get(key):
                entry["w5"][key] = oracle[key]
        for x in ("W1", "W2", "W3", "W4"):
            if x in present and not (cut_of and x == "W3"):
                entry["pairs"].append(compare_w5(x, present[x], r5))
    disagree = [p["pair"] for p in entry["pairs"] if p["outcome"] == "disagree"]
    disagree += ["%s-uncut" % c["way"] for c in entry.get("cuts") or [] if c["outcome"] == "disagree"]
    missing = [w for w in ways if w not in present]
    all_selector = present and all(r.get("selector_error") for r in present.values())
    if op.get("paging", "all") == "all":
        # FORMAT v2 A4: a listing that stopped before its last page is not the whole list (VERIFY-C VC6); an op that
        # declares a page cap ({"pages": N}) or one page ("none") states its own extent
        for w, r in sorted(present.items()):
            if r.get("paging_complete") is False:
                entry["problems"].append("%s's listing stopped before its last page (paging_complete false) on a "
                                         "paging \"all\" op" % w)
    unread = [p for p in entry["pairs"] if p["outcome"] == "not_compared" and p.get("kind") in UNREAD_KINDS]
    if disagree:
        entry["outcome"] = "disagree"
        entry["disagreeing"] = disagree
    elif missing or all_selector or entry["problems"] or unread or any(
            c["outcome"] == "not_compared" for c in entry.get("cuts") or []):
        entry["outcome"] = "inconclusive"
        why = []
        if missing:
            why.append("no record of %s" % ", ".join(missing))
        if all_selector:
            why.append("the selector failed in every way: %s" % gp.compact(
                next(iter(present.values())).get("selector_error")))
        why.extend(sorted(set(str(p.get("detail")) for p in unread)))
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
    """{tool: {"checked": {pointer: {"clause", "calls", "rederived"}}, "consistency_only": {pointer: calls}}} over
    W5's records, list indices written *; `rederived` counts the calls where W5's record marks the pointer re-derived
    from DESIGN (oracle.rederived, DESIGN-v8 §2, D-R5).
    """
    out = {}
    for (sid, oid), r in sorted(recs.get("W5", {}).items(), key=lambda kv: (str(kv[0][0]), str(kv[0][1]))):
        oracle = r.get("oracle") if isinstance(r.get("oracle"), dict) else {}
        t = out.setdefault(r.get("tool"), {"checked": {}, "consistency_only": {}})
        rederived = set(oracle.get("rederived") or []) if isinstance(oracle.get("rederived"), list) else set()
        for ptr, clause in sorted((oracle.get("checked") or {}).items()) if isinstance(oracle.get("checked"),
                                                                                      dict) else []:
            c = t["checked"].setdefault(_norm_pointer(ptr), {"clause": clause, "calls": 0, "rederived": 0})
            c["calls"] += 1
            # DESIGN-v8 §2 (D-R5): a value W5 computes with a DESIGN formula, not read from an engine value
            c["rederived"] += 1 if ptr in rederived else 0
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


# the code a graphlet comparison runs (VERIFY-C VC3: a fix to any of them must be taken up)
CODE_FILES = ("compare_graphlets.py", "graphlet_precheck.py", "graphlet_service_checks.py", "graphlet_oracle.py",
              "graphlet_doc.py", "compare.py", "checks.py", "suite_common.py")


def _tree_digest(run_dir, rel_dirs):
    """sha256 over the sorted (relative path, sha256) of every file under |rel_dirs| (those that exist)."""
    h = hashlib.sha256()
    for rel in rel_dirs:
        root = os.path.join(run_dir, rel)
        if not os.path.isdir(root):
            continue
        for dirpath, dirnames, filenames in os.walk(root):
            dirnames.sort()
            for name in sorted(filenames):
                path = os.path.join(dirpath, name)
                h.update(("%s\0%s\n" % (os.path.relpath(path, run_dir), sc.file_sha256(path))).encode("utf-8"))
    return h.hexdigest()


def input_digests(run_dir, recs=None):
    """{name: sha256 or None} of everything a graphlet comparison reads (VERIFY-C VC3): the manifest; every way's
    records; the job record and the MCP runner's identity records (graphlet_jobs.jsonl, mcp_runs.jsonl, run.json);
    every kept document and twin (arms/*/docs); every export or save file the records name (arms/*/raw/<file.raw>);
    the by-design table the run reads; and the code (CODE_FILES, as this tree holds them)."""
    import compare as cmp
    names = ["manifest.json"] + ["arms/%s/%s" % (gp.WAY_DIRS[w], gp.OPS_FILE) for w in WAYS + ("drift",)] + \
        ["arms/mcp/%s" % gp.JOBS_FILE, "arms/mcp/mcp_runs.jsonl", "arms/mcp/run.json"]
    out = dict((n, sc.file_sha256(os.path.join(run_dir, n))) for n in names)
    docs = sorted(set("arms/%s/docs" % d for d in list(gp.WAY_DIRS.values()) + ["direct-local", "mcp"]))
    out["docs"] = _tree_digest(run_dir, docs)
    if recs is None:
        recs = gp.all_records(run_dir)[0]
    h = hashlib.sha256()
    for way in sorted(recs):
        for key in sorted(recs[way], key=lambda k: (str(k[0]), str(k[1]))):
            path = file_path(run_dir, way, recs[way][key]) if way in gp.WAY_DIRS else None
            if path:
                h.update(("%s\0%s\n" % (os.path.relpath(path, run_dir), sc.file_sha256(path))).encode("utf-8"))
    out["files"] = h.hexdigest()
    table, _ = cmp.table_for(run_dir, None)
    out["by_design_table"] = sc.file_sha256(table) if table else None
    out["code"] = hashlib.sha256("".join("%s\0%s\n" % (f, sc.file_sha256(os.path.join(HERE, f)))
                                         for f in CODE_FILES).encode("utf-8")).hexdigest()
    return out


def compare_run(run_dir, manifest=None, checks_mod=None, extra=None):
    """graphlet_compare.json's content for a run directory (and the pre-check records it judged with); |extra|: the
    stages run after the per-op comparison and before the summary (STAGES)."""
    manifest = manifest if manifest is not None else sc.read_json(os.path.join(run_dir, "manifest.json"))
    recs, warnings, _ = gp.all_records(run_dir)
    attach_files(run_dir, recs)
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
           "manifest": sc.manifest_identity(manifest), "inputs": input_digests(run_dir, recs),
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
    """Compares a run's area G (compare_run with STAGES), from the run's files; writes graphlet_compare.json and the
    pre-check records unless not |write| (the gate and the report compute in memory); returns the comparison."""
    _load_stages()
    result, prechecks = compare_run(run_dir, manifest, extra=STAGES)
    if write:
        gp.write_prechecks(run_dir, prechecks)
        sc.write_json_atomic(out_path or os.path.join(run_dir, COMPARE_FILE), result)
    return result


def stale(run_dir, result):
    """The inputs (input_digests) whose content differs from what |result| read ([] when every one is current)."""
    now = input_digests(run_dir)
    was = (result or {}).get("inputs") or {}
    return [k for k in sorted(now) if was.get(k, "unrecorded") != now[k]]


def load_or_compute(run_dir, manifest=None, reuse=False):
    """The gate's and the report's reading (VERIFY-C VC3): the comparison recomputed from the run's files, in memory,
    every time; with |reuse| (their --reuse flag) RUN_DIR/graphlet_compare.json when every input it read is current
    (stale: records, documents, twin, job and identity records, export files, the by-design table, the code), else
    recomputed. Never written here (compare_graphlets.py writes it)."""
    if reuse:
        data = sc.read_json(os.path.join(run_dir, COMPARE_FILE))
        if isinstance(data, dict) and data.get("compare_version") == COMPARE_VERSION and not stale(run_dir, data):
            return data
    return compute(run_dir, manifest, write=False)


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
