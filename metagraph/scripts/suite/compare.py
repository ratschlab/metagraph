#!/usr/bin/env python3
"""Comparator of the suite: joins the arms of a run per scenario and writes RUN_DIR/compare.json.

    compare.py RUN_DIR [--by-design service_by_design.json] [--reevaluated]

--reevaluated: the run's re-evaluated set (reevaluate.py; merged with merge_run --reevaluated) into
compare.reevaluated.json, the relations judged by the re-evaluation's checks.

Reads the merged run (merge_run.py; it merges first when RUN_DIR/arms exists): run.json, executions.jsonl,
samples.jsonl, manifest.json. Per scenario and arm: the outcome over the repeats, the repeats' digests (equal
or not), wall and server times. Between the direct arms: digest and per-check agreement (with the number of
checks both evaluated). For the MCP arm the outcome classes of DESIGN §1.3: agree (at least one check compared,
none differing) / not_comparable (no check both arms evaluated) / disagree / not_implemented_mcp /
service_by_design (with whether the pair was verified) / missing.
Everything timed keys on the execution's `host_cache` (merge_run: cold = the scenario's first execution on the
host in this run, warm = another one within 60 s, aged = older): overheads of §5 pair warm lines only
(direct-remote minus direct-local: network + front, without the connection set-up, which is reported apart;
MCP minus direct with and without the queue; the MCP phases), and the cold figures are reported on their own;
timing distributions and resource use (samples joined by time) per area, arm and host cache.

Service-side by design: async's table (service_by_design.json next to this script, from the search-service at
commit 22d034ba; --by-design overrides) holds cases {id, direct_expectation, service_answer, agreement}. A
scenario names its cases in `mcp.by_design` (an id or a list). agreement `not_implemented_mcp` makes the
scenario "not implemented yet via MCP"; `service_only` (no host counterpart) makes it service-side by design;
`by_design` makes it service-side by design when the pair holds: the direct answer matches
direct_expectation (http or a list, code or slot_code, the text with digits ignored, text_prefix as is,
text_contains, each where stated; nothing stated: the scenario's own status check) and the MCP answer matches service_answer
(refused_at_submit / relayed: error.code, text_contains, field with alternatives and placeholders where the MCP
execution states them; answered / clamped: no error). A pair that does not hold is a disagreement, naming the side.
Each MCP execution is matched as async's runner matches it (e334715b): the cases the scenario names first, else an
answer-triggered case (answer_match) whose wall rule holds or cannot be decided; a stated key this comparator does
not know leaves its side unverified.
"""

import argparse
import bisect
import hashlib
import json
import os
import re
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import merge_run  # noqa: E402
import suite_common as sc  # noqa: E402

OUTCOME_KEYS = ("pass", "fail", "error", "finding", "inconclusive", "not_implemented_mcp", "service_by_design",
                "mixed")
SAMPLE_BEFORE_S = 1.0
SAMPLE_AFTER_S = 2.0
SAMPLE_NEAREST_S = 5.0


# ---------------------------------------------------------------------------------------------------------------
# service-side by design

DEFAULT_BY_DESIGN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "service_by_design.json")


class TableError(ValueError):
    """The by-design table fails the load validation (table_problems): refused, as async's runner refuses it."""


# The table's vocabulary as async's runner validates it at load (suite_mcp_runner.py ServiceByDesign.problems at
# e334715b, triage-2 B-08: one semantics for the runner and this comparator): what an answer_match may state, what
# each of its sides may state, the direct_mismatch values, the service_answer outcomes the checks judge.
ANSWER_MATCH_TOP_KEYS = ("routes", "direct", "mcp", "direct_mismatch", "wall_rule")
ANSWER_MATCH_KEYS = ("http", "code", "text_contains", "field", "slot_code")
DIRECT_MISMATCH_VALUES = ("disagree", "timing")
WALL_RULE_VALUES = sc.WALL_RULE_VALUES
ANSWERED_OUTCOMES = ("answered", "clamped")
EXPLAINS_ENTRY_KEYS = ("check", "detail_matches")


def entry_problem(entry):
    """What is wrong with one explains_checks entry ('' when nothing): a non-empty check name, or an object stating
    check (a non-empty name) and detail_matches (a regex that compiles) and nothing else (peer review R4)."""
    if isinstance(entry, str):
        return "" if entry else "an empty check name"
    if not isinstance(entry, dict):
        return "%r is neither a check name nor an object" % (entry,)
    unknown = sorted(set(entry) - set(EXPLAINS_ENTRY_KEYS))
    if unknown:
        return "unknown key(s) %s (an entry states %s)" % (unknown, list(EXPLAINS_ENTRY_KEYS))
    if not isinstance(entry.get("check"), str) or not entry["check"]:
        return "an object entry without check (a check name)"
    rx = entry.get("detail_matches")
    if not isinstance(rx, str) or not rx:
        return "an object entry without detail_matches (a regex; a whole check is a plain name)"
    try:
        re.compile(rx, re.S)
    except re.error as exc:
        return "detail_matches %r is not a valid regex (%s)" % (rx, exc)
    return ""


def table_problems(cases):
    """What the load validation refuses, each naming its case (async's ServiceByDesign.problems at e334715b,
    mirrored): an answer_match without routes (a non-empty list of route strings), a key of it outside
    ANSWER_MATCH_TOP_KEYS, a side that is not an object or states a key outside ANSWER_MATCH_KEYS, an MCP side whose
    http is not one integer, that states slot_code (a direct-side key) or that states no code (peer review R5), a
    direct_mismatch outside DIRECT_MISMATCH_VALUES, a wall_rule outside WALL_RULE_VALUES; an explains_checks that is
    not a list of entries (entry_problem) or sits on a case other than answered / clamped."""
    out = []
    for case in cases:
        cid = "case %r" % (case.get("id"),)
        explains = case.get("explains_checks")
        if explains is not None:
            if not isinstance(explains, list):
                out.append("%s: explains_checks must be a list of check names" % cid)
            else:
                out.extend("%s: explains_checks[%d]: %s" % (cid, i, entry_problem(x)) for i, x in enumerate(explains)
                           if entry_problem(x))
            if not judged_by_checks(case):
                out.append("%s: explains_checks on a %s case (only answered / clamped cases leave the outcome to the "
                           "checks)" % (cid, (case.get("service_answer") or {}).get("outcome")))
        if "answer_match" not in case:
            continue
        am = case["answer_match"]
        if not isinstance(am, dict):
            out.append("%s: answer_match must be an object" % cid)
            continue
        unknown = sorted(set(am) - set(ANSWER_MATCH_TOP_KEYS))
        if unknown:
            out.append("%s: unknown answer_match key(s) %s" % (cid, unknown))
        routes = am.get("routes")
        if not isinstance(routes, list) or not routes or not all(isinstance(r, str) and r for r in routes):
            out.append("%s: answer_match without routes (a non-empty list of routes)" % cid)
        for side in ("direct", "mcp"):
            want = am.get(side)
            if not isinstance(want, dict):
                out.append("%s: answer_match.%s must be an object" % (cid, side))
                continue
            unknown = sorted(set(want) - set(ANSWER_MATCH_KEYS))
            if unknown:
                out.append("%s: unknown answer_match.%s key(s) %s" % (cid, side, unknown))
        mcp = am.get("mcp") if isinstance(am.get("mcp"), dict) else {}
        http = mcp.get("http")
        if http is not None and (isinstance(http, bool) or not isinstance(http, int)):
            out.append("%s: answer_match.mcp.http must be one integer, not %r" % (cid, http))
        if "slot_code" in mcp:
            out.append("%s: answer_match.mcp states slot_code (a direct-side key)" % cid)
        if isinstance(am.get("mcp"), dict) and not (isinstance(mcp.get("code"), str) and mcp["code"].strip()):
            out.append("%s: answer_match.mcp states no code (a side without one would make every answer on the "
                       "route that matches the rest by design)" % cid)
        if "direct_mismatch" in am and am["direct_mismatch"] not in DIRECT_MISMATCH_VALUES:
            out.append("%s: answer_match.direct_mismatch %r is not one of %s" % (
                cid, am["direct_mismatch"], list(DIRECT_MISMATCH_VALUES)))
        if "wall_rule" in am and am["wall_rule"] not in WALL_RULE_VALUES:
            out.append("%s: answer_match.wall_rule %r is not one of %s" % (cid, am["wall_rule"], list(WALL_RULE_VALUES)))
    return out


def load_by_design(path):
    """(cases by id, info) from async's table; ({}, None) without one ('' or None: no table). A table that fails the
    load validation (table_problems) raises TableError: refused, as async's runner refuses to start with it."""
    if not path or not os.path.exists(path):
        return {}, None
    with open(path, "rb") as fh:
        data_bytes = fh.read()
    data = json.loads(data_bytes.decode("utf-8"))
    listed = [c for c in data.get("cases", []) if isinstance(c, dict)]
    problems = table_problems(listed)
    if problems:
        raise TableError("the service-by-design table %s is invalid: %s" % (path, "; ".join(problems)))
    cases = dict((c["id"], c) for c in listed if c.get("id"))
    info = {"path": os.path.abspath(path), "service_commit": data.get("service_commit"),
            "version": data.get("version"), "cases": len(cases),
            "sha256": hashlib.sha256(data_bytes).hexdigest()}
    return cases, info


def table_for(run_dir, by_design=None, reevaluated=False):
    """(the table's path, warnings): |by_design| when given ('' = no table), else for a re-evaluated set the table its
    re-evaluation replayed (reevaluated.json by_design.path, B-10) where it is still there, else DEFAULT_BY_DESIGN; a
    warning when the table read differs (sha256) from the one the re-evaluation replayed."""
    warnings = []
    mark = (merge_run.reevaluation_mark(run_dir) or {}) if reevaluated else {}
    used = mark.get("by_design") if isinstance(mark.get("by_design"), dict) else None
    path = by_design
    if path is None:
        path = used.get("path") if used and used.get("path") and os.path.exists(used["path"]) else DEFAULT_BY_DESIGN
    if reevaluated and path:
        sha = None
        if os.path.exists(path):
            with open(path, "rb") as fh:
                sha = hashlib.sha256(fh.read()).hexdigest()
        if used is None:
            warnings.append("the re-evaluation replayed no by-design table; this comparison reads %s" % path)
        elif sha != used.get("sha256"):
            warnings.append("the by-design table read here (%s, sha256 %s) differs from the one the re-evaluation "
                            "replayed (%s, sha256 %s)" % (path, (sha or "missing")[:16], used.get("path"),
                                                          str(used.get("sha256"))[:16]))
    elif reevaluated and used is not None:
        warnings.append("no by-design table read here, the re-evaluation replayed %s" % used.get("path"))
    return path, warnings


def judged_by_checks(case):
    """Whether a matched case leaves the outcome to the checks: an answered / clamped case (the service answered,
    reshaping or limiting the request; triage-2 B-01), never a refused / relayed one (the host never answered)."""
    return (case.get("service_answer") or {}).get("outcome") in ANSWERED_OUTCOMES


def explains(case):
    """A case's explains_checks entries (the checks whose failure the service's path explains: a name, or {check,
    detail_matches}), [] when it states none."""
    out = []
    for x in case.get("explains_checks") or []:
        if isinstance(x, str):
            out.append(x)
        elif isinstance(x, dict) and isinstance(x.get("check"), str) and isinstance(x.get("detail_matches"), str):
            out.append({"check": x["check"], "detail_matches": x["detail_matches"]})
    return out


def entry_explains(entry, check):
    """Whether one explains_checks entry explains a check item: a name, every item of that check; an object {check,
    detail_matches} (peer review R4), an item of that check whose detail the regex matches completely (re.fullmatch,
    re.S)."""
    name = str(check.get("check"))
    if isinstance(entry, str):
        return entry == name
    if str(entry.get("check")) != name:
        return False
    try:
        return re.fullmatch(str(entry.get("detail_matches")), str(check.get("detail") or ""), re.S) is not None
    except re.error:
        return False


def union_of(cases):
    """The union of the matched cases' explains_checks entries (several named answered cases apply together)."""
    entries = [x for c in cases for x in explains(c)]
    return list(dict((json.dumps(x, sort_keys=True), x) for x in entries).values())


def _digits_free(text):
    return re.sub(r"[0-9]+", "#", text or "")


# The keys a side of a table case may state. A key the comparator matches is stated when it is present and not null;
# a descriptive key is never matched. A stated key outside both (a table newer than this comparator) is never
# ignored silently: once the known keys hold, that side is unverified (None), its detail naming the key.
DIRECT_KEYS = ("http", "code", "slot_code", "text", "text_prefix", "text_contains")
DIRECT_NOTES = ("text_source", "behaviour")
SERVICE_KEYS = ("outcome", "code", "text_contains", "field")
SERVICE_NOTES = ("tool", "http", "stage", "detail", "keys", "missing_argument", "reason")
ANSWER_DIRECT_KEYS = ("http", "code", "slot_code", "text_contains")
ANSWER_MCP_KEYS = ("http", "code", "text_contains", "field")
# the MCP runner's host-state codes (suite_common): inconclusive (reevaluate.py, merge_run.py) where the plausibility
# rule admits them, or cannot decide (suite_common.host_state_verdict), and read neither for nor against a case here
MCP_HOST_STATE_CODES = sc.MCP_HOST_STATE_CODES


def unknown_keys(side, known, notes=()):
    """The keys |side| states (present, not null) that are neither matched (|known|) nor descriptive (|notes|)."""
    if not isinstance(side, dict):
        return []
    return sorted(k for k, v in side.items() if v is not None and k not in known and k not in notes)


def _alternatives(value):
    """A stated value's alternatives: "a | b" (outside <...> placeholders), or a list."""
    if isinstance(value, list):
        return [str(v).strip() for v in value if v is not None and str(v).strip()]
    return [x.strip() for x in re.split(r"\|(?![^<]*>)", str(value)) if x.strip()]


def _http_wants(value):
    """A stated http: one status or a list of them."""
    return value if isinstance(value, list) else [value]


def direct_error_parts(line, run_dir):
    """(top-level codes, pattern slot codes, texts) of a direct answer: the error of the line and of each pattern
    slot of the raw answer."""
    top, slots, texts = set(), set(), set()
    err = line.get("error") or {}
    if err.get("code"):
        top.add(err["code"])
    if err.get("text"):
        texts.add(err["text"])
    raw = line.get("raw")
    if raw and run_dir:
        path = os.path.join(run_dir, raw)
        if os.path.exists(path):
            try:
                body = (sc.read_gz_json(path).get("answer") or {}).get("body")
            except (OSError, ValueError):
                body = None
            if isinstance(body, dict):
                for p in body.get("patterns") or []:
                    if isinstance(p, dict) and isinstance(p.get("error"), dict):
                        if p["error"].get("code"):
                            slots.add(p["error"]["code"])
                        if p["error"].get("message"):
                            texts.add(p["error"]["message"])
    return top, slots, texts


def direct_errors(line, run_dir):
    """(codes, texts) of a direct answer: the top-level error's and each pattern slot's (from the raw answer)."""
    top, slots, texts = direct_error_parts(line, run_dir)
    return top | slots, texts


def _raw_doc(line, run_dir):
    raw = (line or {}).get("raw")
    if not raw or not run_dir or not os.path.exists(os.path.join(run_dir, raw)):
        return None
    try:
        obj = sc.read_gz_json(os.path.join(run_dir, raw))
    except (OSError, ValueError, EOFError):
        return None
    return obj if isinstance(obj, dict) else None


def mcp_error(line, run_dir):
    """The MCP execution's error {code, text, field} from its line, else from its raw answer; None when it
    states none."""
    err = line.get("error")
    if isinstance(err, dict) and (err.get("code") or err.get("text")):
        return err
    obj = _raw_doc(line, run_dir)
    if obj is not None:
        for cand in ((obj.get("answer") or {}).get("error") if isinstance(obj.get("answer"), dict) else None,
                     obj.get("error")):
            if isinstance(cand, dict) and (cand.get("code") or cand.get("text")):
                return cand
    return None


# the MCP tool's own error object behind an answer's error (async's runner raw_error), suite_common's
tool_error = sc.tool_error_of


def mcp_view(line, run_dir, with_doc=False):
    """(the MCP answer as the runner's matcher reads it: {http, error, body}, the tool's error object[, the raw
    record]) of an MCP execution: its raw answer where recorded, else its line's http and error."""
    doc = _raw_doc(line, run_dir)
    ans = (doc or {}).get("answer") if isinstance((doc or {}).get("answer"), dict) else None
    if ans is None:
        err = line.get("error") if isinstance(line.get("error"), dict) and line["error"] else None
        ans = {"http": line.get("http"), "error": err, "body": None}
    return (ans, tool_error(doc, ans), doc) if with_doc else (ans, tool_error(doc, ans))


def _direct_problems(side, http, top, slots, texts):
    """(problems, number of stated keys compared) of a direct answer against a case's direct side: http (a status or
    a list), code (a top-level or a pattern slot's; alternatives "a | b"), slot_code (a pattern slot's only), text
    (whole, digits ignored), text_prefix (as is), text_contains (a string or a list)."""
    problems, compared = [], 0
    if side.get("http") is not None:
        compared += 1
        if http not in _http_wants(side["http"]):
            problems.append("http %s, the case expects %s" % (http, side["http"]))
    if side.get("code") is not None:
        compared += 1
        if not set(_alternatives(side["code"])) & (top | slots):
            problems.append("code %s, the case expects %s" % (sorted(top | slots) or None, side["code"]))
    if side.get("slot_code") is not None:
        compared += 1
        if not set(_alternatives(side["slot_code"])) & slots:
            problems.append("slot code %s, the case expects %s" % (sorted(slots) or None, side["slot_code"]))
    if side.get("text") is not None:
        compared += 1
        want = _digits_free(side["text"])
        if not any(_digits_free(t) == want for t in texts):
            problems.append("text %r, the case expects %r" % (sorted(texts)[:1], side["text"]))
    if side.get("text_prefix") is not None:
        compared += 1   # as is, no digit folding (async's table vocabulary, 2026-10-14)
        if not any(t.startswith(side["text_prefix"]) for t in texts):
            problems.append("text %r, the case expects the prefix %r" % (sorted(texts)[:1], side["text_prefix"]))
    if side.get("text_contains") is not None:
        compared += 1
        frags = side["text_contains"] if isinstance(side["text_contains"], list) else [side["text_contains"]]
        if not any(all(str(f) in t for f in frags) for t in texts):
            problems.append("text %r lacks %r" % (sorted(texts)[:1], side["text_contains"]))
    return problems, compared


def states_host_answer(expectation):
    """Whether a case's direct_expectation states anything to compare (a matched key, present and not null)."""
    return any((expectation or {}).get(k) is not None for k in DIRECT_KEYS)


def match_direct(expectation, line, run_dir):
    """(True / False / None, detail): does the direct answer match a case's direct_expectation (None: no direct
    execution, nothing stated to compare, or a stated key this comparator does not know while the known ones hold;
    a known key that does not hold is a mismatch whatever else is stated, B-05)."""
    if line is None:
        return None, "no direct execution"
    top, slots, texts = direct_error_parts(line, run_dir)
    problems, compared = _direct_problems(expectation or {}, line.get("http"), top, slots, texts)
    if problems:
        return False, "; ".join(problems)
    unknown = unknown_keys(expectation, DIRECT_KEYS, DIRECT_NOTES)
    if unknown:
        return None, "the case's direct_expectation states %s, which this comparator does not match" % unknown
    return (True if compared else None), ("matches" if compared else "nothing stated to compare")


def status_verdict(line):
    """(True / False / None, detail) of the scenario's own status check(s) on a direct execution: the direct side of a
    case that states no host answer (G11), judged by what the scenario itself expects of the host."""
    if line is None:
        return None, "no direct execution"
    st = [c for c in line.get("checks") or [] if c.get("check") == "status"]
    if not st:
        return None, "the case states no host answer, and the reference direct execution has no status check"
    outs = [c.get("outcome") for c in st]
    detail = "; ".join(str(c.get("detail") or "")[:160] for c in st)
    if "fail" in outs:
        return False, "the case states no host answer; the scenario's own status check failed on repeat %s (%s)" % (
            line.get("repeat"), detail)
    if all(o == "pass" for o in outs):
        return True, "the case states no host answer; the scenario's own status check passed on repeat %s (%s)" % (
            line.get("repeat"), detail)
    return None, "the case states no host answer; the scenario's own status check reads %s on repeat %s" % (
        outs, line.get("repeat"))


def _field_matches(observed, case_field):
    """A case's field: alternatives "a | b" (split as its code), each with `<...>` placeholders that match anything,
    dots included, as async's runner (suite_mcp_runner.py ServiceByDesign._like) matches them: "bounds.<knob>",
    "predicate<path of the first fault>" (predicate.at_least.n)."""
    if observed is None:
        return False
    for alt in _alternatives(case_field):
        rx = ".*?".join(re.escape(part) for part in re.split(r"<[^<>]+>", alt))
        if re.fullmatch(rx, str(observed), re.S) is not None:
            return True
    return False


def host_state_line(line, run_dir, scenario=None):
    """Whether an MCP execution is inconclusive by a host state (MCP_HOST_STATE_CODES) that the scenario does not
    expect and whose wall the plausibility rule admits, or cannot decide (rule (c)) (suite_common.host_state_verdict
    on the raw record's service-side wall; B-09, A-11): it reads neither for nor against a case or the direct arm.
    Only an answer refused with a host-state code: an answer-triggered case's undecidable wall rule (a relayed
    /pattern deadline, code search_failed) is no host state (B5b-09; mcp_agreement's undecidable branch)."""
    if line.get("outcome") != "inconclusive":
        return False
    code = (line.get("error") or {}).get("code") if isinstance(line.get("error"), dict) else None
    why = line.get("inconclusive_reason")
    if why is not None and code is not None:
        return str(code) in MCP_HOST_STATE_CODES and (why in sc.UNDECIDABLE or str(why) in MCP_HOST_STATE_CODES)
    ans, terr, doc = mcp_view(line, run_dir, with_doc=True)
    code, state, _ = sc.host_state_verdict(scenario or {}, line, doc, ans.get("error") or mcp_error(line, run_dir), terr)
    return code is not None and state in ("host_state", "undecidable")


def match_service(service_answer, lines, run_dir):
    """(True / False / None, detail): do the MCP executions match a case's service_answer, each decided as async's
    runner decides it (service_answers with the recorded tool error, B-02): None when no MCP execution, when none
    recorded the service's answer (answered / clamped), or when the case states a key this comparator does not know
    while the known ones hold (B-05)."""
    if not lines:
        return None, "no MCP execution"
    outcome = service_answer.get("outcome")
    for line in lines:
        ans, terr = mcp_view(line, run_dir)
        if outcome in ("answered", "clamped") and not (ans.get("error") or {}).get("code"):
            continue   # verified below only by a line that records the service's answer (review N11)
        if not service_answers(service_answer, ans, terr):
            err = ans.get("error") or {}
            if outcome in ("answered", "clamped"):
                return False, "repeat %s: http %s, error %s, the case expects the service to answer (%s)" % (
                    line.get("repeat"), ans.get("http"), err.get("code"), outcome)
            return False, "repeat %s: http %s, code %s, text %r, field %s; the case expects %s %s%s%s" % (
                line.get("repeat"), ans.get("http"), err.get("code"), str(err.get("text") or "")[:160],
                (terr or {}).get("field", err.get("field")), outcome, service_answer.get("code"),
                (" containing %r" % service_answer["text_contains"]) if service_answer.get("text_contains") else "",
                (" field %s" % service_answer["field"]) if service_answer.get("field") else "")
    unknown = unknown_keys(service_answer, SERVICE_KEYS, SERVICE_NOTES)
    if unknown:
        return None, "the case's service_answer states %s, which this comparator does not match" % unknown
    if outcome in ("answered", "clamped"):
        # Verified only when the service's answer is on record (review N11): a line without an error code says
        # nothing when it never reached the service.
        verified = any(_service_answered(x) for x in lines)
        return (True if verified else None), ("matches" if verified else "no answer of the service on record")
    return True, "matches"


def _service_answered(line):
    """Whether an MCP line records the service's answer: a 2xx status, or a backend call (sc.reached_host)."""
    http = line.get("http")
    return (isinstance(http, int) and not isinstance(http, bool) and 200 <= http < 300) or sc.reached_host(line)


def _direct_side(case, direct_line, run_dir):
    """(True / False / None, detail, by) of a case's direct side on the reference direct execution: its
    direct_expectation, or the scenario's own status check when it states no host answer (G11)."""
    d_exp = case.get("direct_expectation") or {}
    if direct_line is not None and not states_host_answer(d_exp) and not unknown_keys(d_exp, DIRECT_KEYS, DIRECT_NOTES):
        ok, detail = status_verdict(direct_line)
        return ok, detail, "scenario status check"
    ok, detail = match_direct(d_exp, direct_line, run_dir)
    return ok, detail, None


def by_design_class(case, direct_line, mcp_lines, run_dir, scenario=None):
    """The MCP outcome class a refused / relayed / service-only table case gives the MCP executions it names
    (DESIGN §1.3), with the pair's verification.

    `verified` is True only when every side the case has was checked and held: both sides for `by_design`,
    the MCP side for `service_only` (no host counterpart); a side with nothing to compare leaves it False
    (counted as an unverified by-design pair in the report). A case that states no host answer (no http, code,
    text, ... in its direct_expectation) has its direct side judged by the scenario's own status check on the
    reference direct execution (G11: pass verified, fail a disagreement, else unverified). An MCP execution that is
    inconclusive by a host state (host_state_line) is skipped: it reads neither for nor against the case; with none
    left the pair is not comparable (B-09)."""
    agreement = case.get("agreement")
    base = {"case": case.get("id"), "case_agreement": agreement,
            "service_outcome": (case.get("service_answer") or {}).get("outcome")}
    if agreement == "not_implemented_mcp":
        return dict(base, **{"class": "not_implemented_mcp"})
    skipped = [x for x in mcp_lines if host_state_line(x, run_dir, scenario)]
    judged = [x for x in mcp_lines if x not in skipped]
    if skipped:
        base["host_state_repeats"] = [x.get("repeat") for x in skipped]
    if mcp_lines and not judged:
        return dict(base, **{"class": "not_comparable", "by": "host state",
                             "detail": "every MCP execution is inconclusive by a host state (%s): none reads for or "
                                       "against the case" % sorted(set(str((mcp_error(x, run_dir) or {}).get("code"))
                                                                       for x in skipped))})
    m_ok, m_detail = match_service(case.get("service_answer") or {}, judged, run_dir)
    base.update({"mcp_verified": m_ok, "mcp_detail": m_detail})
    if agreement == "service_only":
        d_ok = True
        base["direct_detail"] = "no host counterpart"
    else:
        d_ok, d_detail, by = _direct_side(case, direct_line, run_dir)
        if by:
            base["direct_by"] = by
        base.update({"direct_verified": d_ok, "direct_detail": d_detail})
    if not mcp_lines:
        return dict(base, **{"class": "missing"})
    if d_ok is False or m_ok is False:
        base["side"] = "direct" if d_ok is False else "mcp"
        return dict(base, **{"class": "disagree"})
    base["verified"] = d_ok is True and m_ok is True
    return dict(base, **{"class": "service_by_design"})


# ---------------------------------------------------------------------------------------------------------------
# per-scenario joins

def aggregate_outcome(outcomes):
    """The outcome of a scenario on one arm over its repeats: the common one, else error > fail > finding >
    inconclusive > mixed (a repeat that could not test its property keeps the scenario from reading pass)."""
    vals = [o for o in outcomes if o]
    if not vals:
        return None
    if len(set(vals)) == 1:
        return vals[0]
    for worst in ("error", "fail", "finding", "inconclusive"):
        if worst in vals:
            return worst
    return "mixed"


def check_outcomes(checks):
    """{check name: [outcome, ...]} in the order listed (a name may appear more than once)."""
    out = {}
    for c in checks or []:
        out.setdefault(c.get("check"), []).append(c.get("outcome"))
    return out


def compare_checks(a, b):
    """(differing, unmatched, compared, not_exposed, inconclusive) of two executions' checks: names whose
    outcomes differ (pairing a name's occurrences in order), names only one side evaluated (`unmatched`, not a
    disagreement), the number of check results both sides evaluated (pass / fail / finding on both), and the
    numbers skipped because either side read not_exposed, or inconclusive (could not test its property)."""
    ca, cb = check_outcomes(a), check_outcomes(b)
    differing, unmatched = [], []
    compared = not_exposed = inconclusive = 0
    for name in sorted(set(ca) | set(cb), key=lambda n: str(n)):
        if name not in ca or name not in cb:
            unmatched.append(name)
            continue
        differs = False
        for oa, ob in zip(ca[name], cb[name]):
            if "not_exposed" in (oa, ob):
                not_exposed += 1
                continue
            if "inconclusive" in (oa, ob):
                inconclusive += 1
                continue
            compared += 1
            if oa != ob and not differs:
                differing.append({"check": name, "a": oa, "b": ob})
                differs = True
    return differing, unmatched, compared, not_exposed, inconclusive


def group_lines(rows):
    """{scenario: {arm: {repeat: line}}} and the scenario order of first appearance."""
    grouped, order = {}, []
    for r in rows:
        sid = r.get("scenario")
        if sid not in grouped:
            grouped[sid] = {}
            order.append(sid)
        grouped[sid].setdefault(r.get("arm"), {})[r.get("repeat", 0)] = r
    return grouped, order


def _num(v):
    return v if isinstance(v, (int, float)) and not isinstance(v, bool) else None


def _diff(a, b):
    a, b = _num(a), _num(b)
    return None if a is None or b is None else a - b


def _server_elapsed(line):
    return _num((line.get("server") or {}).get("elapsed_ms"))


def host_cache(line):
    """The execution's host cache state (merge_run); for a run merged before host_cache existed, the arm's
    own cold / warm mark."""
    hc = line.get("host_cache")
    if hc:
        return hc
    return "cold" if line.get("cache", "cold") == "cold" else "warm"


def wall_net(line):
    """The wall time without the connection set-up: the runner times requests on an established connection
    (client.reused present); a line of an older runner, which timed from before the connect, minus its
    connect_ms."""
    wall = _num(line.get("wall_ms"))
    client = line.get("client") or {}
    if wall is None or "reused" in client or "attempts" in client:
        return wall
    connect = _num(client.get("connect_ms"))
    return wall - connect if connect is not None else wall


def first_with(lines, label, prefer=None):
    """The arm's first execution (by repeat) with host cache `label`; `prefer` a repeat to take when it has it."""
    if prefer is not None and prefer in lines and host_cache(lines[prefer]) == label:
        return lines[prefer]
    for k in sorted(lines):
        if host_cache(lines[k]) == label:
            return lines[k]
    return None


def arm_view(lines):
    reps = [lines[k] for k in sorted(lines)]
    digests = [r.get("answer_digest") for r in reps]
    return {
        "outcome": aggregate_outcome([r.get("outcome") for r in reps]),
        "outcomes": [r.get("outcome") for r in reps],
        "http": [r.get("http") for r in reps],
        "digests_equal": (len(set(digests)) == 1) if len(digests) > 1 else None,
        "wall_ms": [r.get("wall_ms") for r in reps],
        "server_elapsed_ms": [_server_elapsed(r) for r in reps],
        "failed_checks": sorted(set(c.get("check") for r in reps for c in r.get("checks") or []
                                    if c.get("outcome") == "fail"), key=str),
    }


def _own_time_limitation(line):
    """What made a direct execution's own answer time-limited: its recorded time_limited (checks.time_limited: a
    stated time stop, a 503 deadline and its relays), else a 503 deadline of a line recorded before that; ''."""
    if line.get("time_limited"):
        return str(line["time_limited"])
    err = line.get("error") if isinstance(line.get("error"), dict) else {}
    if line.get("http") == 503 and err.get("code") == "deadline":
        return "a 503 deadline"
    return ""


# SPEC-pattern-search §7.4: the relations a count states; `tested` (§19.7: the decisions a selection pass made,
# work done) and `estimate` (not a bound) are no population counts
COUNT_RELATIONS = ("exact", "at_least", "bounds", "unknown")
NOT_POPULATION = ("tested", "estimate")


def count_range(count):
    """The true values a count object admits (SPEC-pattern-search §7.4), B5b-02: exact v -> (v, v); at_least v ->
    (v, inf); bounds -> (lower, upper) (lower defaults to its value); unknown (the phase never ran) -> None (no
    constraint); a relation outside the vocabulary, or a bound that is no number -> a str naming the contradiction
    (fail closed)."""
    rel = count.get("relation")
    if rel == "unknown":
        return None
    if rel == "exact":
        v = _num(count.get("value"))
        return (v, v) if v is not None else "exact without a numeric value"
    if rel == "at_least":
        v = _num(count.get("value"))
        return (v, float("inf")) if v is not None else "at_least without a numeric value"
    if rel == "bounds":
        lo = _num(count.get("lower", count.get("value")))
        hi = _num(count.get("upper"))
        if lo is None or hi is None:
            return "bounds without numeric lower and upper"
        return (lo, hi) if lo <= hi else "bounds [%s, %s] empty" % (lo, hi)
    return "relation %r outside %s" % (rel, "/".join(COUNT_RELATIONS))


def _describe_count(count):
    rel = count.get("relation")
    if rel == "bounds":
        return "bounds [%s, %s]" % (count.get("lower", count.get("value")), count.get("upper"))
    return "%s %s" % (rel, count.get("value"))


def _count_problems(stopped, complete, path=""):
    """The population counts of a stopped answer's entry that contradict the complete answer's count at the same
    path (B5b-02, VERIFY-REVIEW's rule): each count is the range of true values it admits (count_range), and two
    counts are consistent when their ranges overlap (a stopped exact v against a complete at_least w: v >= w; a
    stopped bounds [l, u] against a complete exact c: l <= c <= u; ...); unknown admits anything; a relation outside
    the vocabulary is a contradiction. `tested` and `estimate` are skipped (NOT_POPULATION)."""
    out = []
    if isinstance(stopped, dict) and "relation" in stopped and isinstance(complete, dict) and "relation" in complete:
        ra, rb = count_range(stopped), count_range(complete)
        where = path or "count"
        for r, side, c in ((ra, "stopped", stopped), (rb, "complete", complete)):
            if isinstance(r, str):
                out.append("%s %s: %s" % (where, side, r))
        if isinstance(ra, tuple) and isinstance(rb, tuple) and max(ra[0], rb[0]) > min(ra[1], rb[1]):
            out.append("%s %s vs %s: no common value" % (where, _describe_count(stopped), _describe_count(complete)))
    if isinstance(stopped, dict):
        for k, v in stopped.items():
            if str(k) in NOT_POPULATION:
                continue
            if isinstance(v, dict) and isinstance(complete, dict):
                out += _count_problems(v, complete.get(k), (path + "." if path else "") + str(k))
    return out


def _released_key(row):
    """A released result's identity: a context (kmer, offset, strand) or a path (sequence, strand)."""
    st = row.get("strand", row.get("orientation"))
    if row.get("kmer") is not None:
        return ("context", row.get("kmer"), row.get("offset"), st, row.get("mismatches", 0))
    return ("path", row.get("sequence"), st, row.get("mismatches", 0))


def _released_problems(stopped, complete):
    """The stopped entry's released results that the complete entry does not list although it states
    retrieval_complete (VERIFY-2 item 5's subset rule; SPEC-pattern-search §7.5: then every context is listed)."""
    if not isinstance(stopped, dict) or not isinstance(complete, dict) or \
            complete.get("retrieval_complete") is not True:
        return []
    mine = [r for r in stopped.get("results") or [] if isinstance(r, dict)]
    theirs = set(_released_key(r) for r in complete.get("results") or [] if isinstance(r, dict))
    missing = [_released_key(r) for r in mine if _released_key(r) not in theirs]
    if missing:
        return ["%d of %d released results not in the complete list (first %s)" % (len(missing), len(mine),
                                                                                    list(missing[0][1:]))]
    return []


def _walk_prefix_problems(stopped, complete, checks_mod):
    """The stopped /traverse answer's paths that are no prefix of a path of the complete one (per seed and arm, the
    flank in walking order), through the checks' reading of a result (a graphlet included)."""
    out = []
    sb, cb = stopped.get("body") or {}, complete.get("body") or {}
    for i in range(min(len(sb.get("results") or []), len(cb.get("results") or []))):
        try:
            rs, rc = checks_mod._result(stopped, i), checks_mod._result(complete, i)
        except Exception:  # an answer the checks cannot read: not shown consistent
            return ["seed %d: not readable" % i]
        for side, arm in (rs.get("arms") or {}).items():
            other = (rc.get("arms") or {}).get(side)
            if other is None:
                out.append("seed %d %s: no arm on the complete side" % (i, side))
                continue
            try:
                fs = [f if side == "right" else f[::-1] for _, f in checks_mod._arm_flanks(arm, side)]
                fc = [f if side == "right" else f[::-1] for _, f in checks_mod._arm_flanks(other, side)]
            except Exception:
                return ["seed %d %s: segment sequences not readable" % (i, side)]
            bad = [f for f in fs if not any(c.startswith(f) for c in fc)]
            if bad:
                out.append("seed %d %s: a %d bp path is no prefix of a complete path" % (i, side, len(bad[0])))
    return out


def content_consistent(stopped_line, complete_line, route, run_dir, checks_mod=None):
    """(True / False / None, detail): whether a time-limited answer's content is consistent with the other arm's
    complete one (VERIFY-2 item 5): a 503 carries no content; /pattern: every population count of each pattern entry
    consistent (_count_problems: the ranges the two counts admit overlap, B5b-02) and, where the complete entry
    states retrieval_complete, every released result among its results (_released_problems); /traverse: every
    stopped path a prefix of a complete path; None where the content cannot be read or the route has no rule (a
    /resolve stop)."""
    if stopped_line.get("http") == 503:
        return True, "a 503 carries no content"
    a, b = load_answer(stopped_line, run_dir), load_answer(complete_line, run_dir)
    if a is None or b is None or a.get("http") != 200 or b.get("http") != 200:
        return None, "no pair of answers to compare"
    if route == "POST /pattern":
        pa, pb = (a.get("body") or {}).get("patterns") or [], (b.get("body") or {}).get("patterns") or []
        if len(pa) != len(pb):
            return False, "%d vs %d pattern entries" % (len(pa), len(pb))
        probs = []
        for i, (ea, eb) in enumerate(zip(pa, pb)):
            probs += ["p%d %s" % (i, x) for x in _count_problems((ea or {}).get("counts"), (eb or {}).get("counts"))]
            probs += ["p%d %s" % (i, x) for x in _released_problems(ea, eb)]
        return (False, "; ".join(probs[:4])) if probs else (True, "every count consistent with the complete answer")
    if route == "POST /traverse" and checks_mod is not None:
        probs = _walk_prefix_problems(a, b, checks_mod)
        return (False, "; ".join(probs[:4])) if probs else (True, "every stopped path a prefix of a complete one")
    return None, "no consistency rule for %s" % route


def _determinism_only(line):
    """Whether a direct execution's only non-pass check is `deterministic` left inconclusive by an earlier repeat's
    time limitation (VERIFY-2 item 5's extension for G27)."""
    rest = [c for c in line.get("checks") or [] if c.get("outcome") != "pass"]
    return bool(rest) and all(c.get("check") == "deterministic" and c.get("outcome") == "inconclusive" and (
        "time limitation" in str(c.get("detail") or "") or "time-limited" in str(c.get("detail") or ""))
        for c in rest)


def timing_repeat(a, b, route, run_dir, checks_mod=None):
    """(True, detail) when two direct executions of one repeat whose outcomes differ differ by timing alone
    (VERIFY-2 item 5, under all its conditions): exactly one side's own answer states a time limitation, that side
    reads inconclusive (never fail or finding), the other pass or inconclusive, no check both evaluated differs, and
    the stopped answer is consistent with the complete one (content_consistent); or (the extension for G27) one side
    is inconclusive only by its determinism check, an earlier repeat having been time-limited, the other passes, and
    the two answers are the same (digest). Else (False, why)."""
    if compare_checks(a.get("checks"), b.get("checks"))[0]:
        return False, "a check differs"
    ta, tb = _own_time_limitation(a), _own_time_limitation(b)
    if bool(ta) != bool(tb):
        stopped, complete, why = (a, b, ta) if ta else (b, a, tb)
        if stopped.get("outcome") != "inconclusive" or complete.get("outcome") not in ("pass", "inconclusive"):
            return False, "outcomes %s (time-limited) and %s" % (stopped.get("outcome"), complete.get("outcome"))
        ok, detail = content_consistent(stopped, complete, route, run_dir, checks_mod)
        if ok is not True:
            return False, "the time-limited answer's content: %s" % detail
        return True, "%s on %s only (%s); %s" % (why, stopped.get("arm"), stopped.get("outcome"), detail)
    if not ta and not tb:
        for x, y in ((a, b), (b, a)):
            if x.get("outcome") == "inconclusive" and y.get("outcome") == "pass" and _determinism_only(x) and \
                    x.get("answer_digest") and x.get("answer_digest") == y.get("answer_digest"):
                return True, "%s: only its determinism check is inconclusive (an earlier repeat was time-limited); " \
                             "the same answer as %s" % (x.get("arm"), y.get("arm"))
    return False, "no time limitation tells the outcomes apart"


def planned_duplicate(a, b, scenario):
    """Why two direct executions of one repeat differ by a planned 409 duplicate, else None (A2-06 of the review of
    lane S5b): the scenario makes an id unique per run and arm (run_unique) and allows a 409 (a retained attempt's
    duplicate), one side answered 200 and the other 409, and neither failed (the 409's shape is checked:
    answer_fields when_http 409); e.g. T-attempts-0002 when one arm's attempt was retained and the other's not."""
    if not (scenario or {}).get("run_unique"):
        return None
    allowed = []
    for e in (scenario or {}).get("expect") or []:
        if e.get("check") == "status":
            w = (e.get("args") or {}).get("http", 200)
            allowed = w if isinstance(w, list) else [w]
    if 409 not in allowed or sorted([a.get("http"), b.get("http")], key=str) != [200, 409]:
        return None
    if any(x.get("outcome") not in ("pass", "inconclusive") for x in (a, b)):
        return None
    dup = a if a.get("http") == 409 else b
    return "%s met its run's retained id (409 duplicate, planned: run_unique %s)" % (
        dup.get("arm"), scenario.get("run_unique"))


def direct_agreement(local, remote, route=None, run_dir=None, checks_mod=None, scenario=None):
    """Agreement of the two direct arms on every repeat both ran: identical (equal digests), disagree (a check differs,
    or the outcomes differ other than by timing or a planned 409 duplicate), timing (the outcomes differ only on
    repeats where timing_repeat holds: counted apart, VERIFY-2 item 5), not_comparable (no check compared),
    same_checks. A repeat whose outcomes differ by a planned 409 duplicate (planned_duplicate, A2-06) is listed
    under planned_duplicates and reads neither way."""
    common = sorted(set(local) & set(remote))
    if not common:
        return None
    digest_equal = all(local[r].get("answer_digest") == remote[r].get("answer_digest") for r in common)
    differing = []
    compared = 0
    timing, untimed, planned = [], [], []
    for r in common:
        d, _, n, _, _ = compare_checks(local[r].get("checks"), remote[r].get("checks"))
        compared += n
        for x in d:
            x["repeat"] = r
        differing.extend(d)
        if local[r].get("outcome") != remote[r].get("outcome"):
            dup = planned_duplicate(local[r], remote[r], scenario)
            if dup and not d:
                planned.append({"repeat": r, "detail": dup})
                continue
            ok, why = timing_repeat(local[r], remote[r], route, run_dir, checks_mod)
            (timing if ok else untimed).append({"repeat": r, "detail": why})
    outcome_equal = all(local[r].get("outcome") == remote[r].get("outcome") for r in common)
    if digest_equal:
        cls = "identical"
    elif differing or untimed:
        cls = "disagree"
    elif timing:
        cls = "timing"
    elif not compared:
        cls = "not_comparable"
    else:
        cls = "same_checks"
    out = {"class": cls, "digest_equal": digest_equal, "outcome_equal": outcome_equal,
           "differing_checks": differing, "checks_compared": compared, "repeats": common}
    if timing:
        out["timing_repeats"] = timing
    if planned:
        out["planned_duplicates"] = planned
    if untimed and cls == "disagree":
        out["outcome_differs"] = untimed
    return out


def load_checks_module(checks_dir=None):
    """checks.py of |checks_dir| (default: next to this script) as a module, or None (direct_agreement's content
    rule reads a /traverse result through it)."""
    path = os.path.join(checks_dir or os.path.dirname(os.path.abspath(__file__)), "checks.py")
    if not os.path.exists(path):
        return None
    import importlib.util
    spec = importlib.util.spec_from_file_location("suite_checks_for_compare", path)
    mod = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(mod)
    except Exception:  # a checks.py that does not load: no /traverse content rule
        return None
    return mod


def answer_cases(cases):
    """The table cases an answer triggers (`answer_match`, whatever the scenario names): {routes, direct, mcp,
    direct_mismatch}; the direct side {http (a status or a list), code (a top-level or a pattern slot's; "a | b"),
    slot_code (a pattern slot's only), text_contains}, the mcp side {http (one status), code, text_contains, field
    (placeholders and alternatives as _field_matches)}; e.g. pattern.deadline_ends_task, the host's 503 deadline the
    service relays as a failed search."""
    return [c for c in (cases or {}).values() if isinstance(c.get("answer_match"), dict)]


def answer_side_matches(side, http, codes, texts, field=None, slot_codes=None, known=None, http_list=True):
    """Whether one side of an answer_match holds for an answer (its status, error codes, texts, the MCP error's
    field, the pattern slots' codes), as async's runner reads it (suite_mcp_runner.py ServiceByDesign.answer_matches:
    a field the answer does not name does not match; slot_code the codes of body.patterns[].error; on the MCP side
    http is one status compared by equality, |http_list| False): True or False; None when the known keys hold but the
    side states a key outside |known| (default: the union of both sides' keys), which is never ignored silently. A
    side that states nothing never matches."""
    if not isinstance(side, dict):
        return False
    known = known or tuple(sorted(set(ANSWER_DIRECT_KEYS) | set(ANSWER_MCP_KEYS)))
    if not any(side.get(k) is not None for k in known) and not unknown_keys(side, known):
        return False
    slots = set(slot_codes or ())
    if side.get("http") is not None and (http not in _http_wants(side["http"]) if http_list else http != side["http"]):
        return False
    if side.get("code") is not None and not set(_alternatives(side["code"])) & (set(codes) | slots):
        return False  # a code or alternatives "a | b", as the MCP runner's matcher reads them
    if side.get("slot_code") is not None and "slot_code" in known and not set(_alternatives(side["slot_code"])) & slots:
        return False
    if side.get("field") is not None and "field" in known and not _field_matches(field, side["field"]):
        return False
    tc = side.get("text_contains")
    if tc and not any(tc in (t or "") for t in texts):
        return False
    return None if unknown_keys(side, known) else True


def answer_case_of(cases, route, http, error, skip=()):
    """The first answer-triggered case (table order) on the scenario's route whose MCP side the MCP answer (status,
    error {code, text, field}) matches, else None (async's runner: an answer_match needs routes, B-08)."""
    err = error if isinstance(error, dict) else {}
    for case in answer_cases(cases):
        if case.get("id") in skip:
            continue
        am = case["answer_match"]
        if not am.get("routes") or route not in am["routes"]:
            continue
        if answer_side_matches(am.get("mcp"), http, set([err.get("code")]) - {None}, [err.get("text")],
                               field=err.get("field"), known=ANSWER_MCP_KEYS, http_list=False) is True:
            return case
    return None


def named_cases(cases, scenario):
    """The table cases a scenario names, as async's runner reads them (suite_mcp_runner.py ServiceByDesign.named):
    mcp.by_design (an id or a list), the scenario's service_by_design, its `sbd:<id>` tags, in that order, then every
    case listing the scenario's id in scenario_ids."""
    scenario = scenario or {}
    ids = []
    for own in (((scenario.get("mcp") or {}).get("by_design")), scenario.get("service_by_design")):
        ids += [own] if isinstance(own, str) else [str(x) for x in own or [] if isinstance(x, str)]
    ids += [t[4:] for t in scenario.get("tags") or [] if isinstance(t, str) and t.startswith("sbd:")]
    out = [cases[i] for i in dict.fromkeys(ids) if i in (cases or {})]
    sid = str(scenario.get("id"))
    out += [c for c in (cases or {}).values() if sid in (c.get("scenario_ids") or []) and c not in out]
    return out


def service_answers(want, answer, tool_error=None):
    """Whether an MCP answer matches a named case's service_answer, as async's runner reads it (ServiceByDesign.answers):
    answered / clamped: no error and http 200; refused_at_submit / relayed: an error whose code is the case's
    (`a | b`), whose text contains text_contains, and whose field (the tool's error object's, else the answer's)
    matches the case's field (alternatives, placeholders; a case field against no field does not match)."""
    want = want or {}
    answer = answer or {}
    outcome = want.get("outcome")
    err = answer.get("error")
    if outcome in ("answered", "clamped"):
        return err is None and answer.get("http") == 200
    if outcome not in ("refused_at_submit", "relayed") or not isinstance(err, dict):
        return False
    if want.get("code") is not None and str(err.get("code")) not in _alternatives(want["code"]):
        return False
    text = str(err.get("text") or "")
    if want.get("text_contains") and not any(part.strip() in text for part in str(want["text_contains"]).split("|")):
        return False   # "a | b" = either (triage-2 B-07: pattern.predicate_form_other, one text per fault)
    if want.get("field") is None:
        return True
    return _field_matches((tool_error or {}).get("field", err.get("field")), want["field"])


def case_matches(cases, scenario, answer, tool_error=None, line=None, doc=None, walls=None, notes=None):
    """(the cases an MCP answer matches, how): async's runner's ServiceByDesign.matches_why at e334715b: every case
    the scenario names whose service_answer the answer matches, in the order named ("named"; an answer matches
    either refusals or answers, never both, and several named answered cases apply together); else the first
    answer-triggered case (table order) on its route whose MCP side it matches and whose answer_match.wall_rule (if
    any) holds for the execution or cannot be decided (|line| and its raw record |doc|: suite_common.wall_rule_verdict,
    peer review R1, rule (c); a case a known wall holds back does not apply and matching continues)
    ("answer-triggered"); ([], None). |walls| (a list) receives each wall rule's (state, detail, case id): an
    undecidable one last, with its case returned; |notes| (a list) why the client's wall decided nothing
    (suite_common.wall_rule_verdict)."""
    if answer is None:
        return [], None
    named = [c for c in named_cases(cases, scenario) if service_answers(c.get("service_answer"), answer, tool_error)]
    if named:
        return named, "named"
    err = answer.get("error") if isinstance(answer.get("error"), dict) else {}
    if isinstance(tool_error, dict) and "field" in tool_error:
        err = dict(err, field=tool_error.get("field"))
    held = set()
    while True:
        case = answer_case_of(cases, (scenario or {}).get("route"), answer.get("http"), err, skip=held)
        if case is None:
            return [], None
        rule = case["answer_match"].get("wall_rule")
        if rule is None:
            return [case], "answer-triggered"
        state, detail = sc.wall_rule_verdict(str(rule), str(case.get("id")), scenario, answer, line, doc, notes)
        if walls is not None:
            walls.append((state, detail, case.get("id")))
        if state in ("holds", "undecidable"):
            return [case], "answer-triggered"
        held.add(case.get("id"))


def named_case_of(cases, scenario, answer, tool_error=None):
    """The first case the scenario names whose service_answer the MCP answer matches (the runner's matcher), else
    None."""
    for case in named_cases(cases, scenario):
        if service_answers(case.get("service_answer"), answer, tool_error):
            return case
    return None


def explained_outcome(checks, explains_names, rule):
    """(the outcome, E) of an execution an answered / clamped case matched (triage-2 B-01; async's runner
    explained_outcome at e334715b, mirrored): F = the failing checks (outcome fail, or one outside the checks'
    vocabulary), E = those of F with outcome fail that an entry of the case's explains_checks explains
    (entry_explains: a name, or {check, detail_matches}), U = F - E. fail when
    U is non-empty; else r = |rule| over the checks without E: finding when r is finding (a by-design mark never
    erases one), service_by_design when E is non-empty, else r (a case explaining nothing leaves the checks' outcome).
    E lists each explained check's name once, in check order."""
    failing = [c for c in checks if str(c.get("outcome")) == "fail" or str(c.get("outcome")) not in sc.CHECK_OUTCOMES]
    e_items = [c for c in failing if str(c.get("outcome")) == "fail"
               and any(entry_explains(entry, c) for entry in explains_names)]
    explained = list(dict.fromkeys(str(c.get("check")) for c in e_items))
    if len(e_items) < len(failing):
        return "fail", explained
    rest = [dict(c) for c in checks if all(c is not x for x in e_items)]
    r = str(rule(rest))
    if r == "finding":
        return "finding", explained
    if e_items:
        return "service_by_design", explained
    return r, explained


def answered_verdict(matched, checks, rule):
    """(the outcome, E, the recorded case's id) of an execution answered / clamped cases matched (async's judge at
    e334715b): explained_outcome against the union of their explains_checks; the case recorded is the first matched
    one that explains an item of E, else the first matched."""
    outcome, explained = explained_outcome(checks, union_of(matched), rule)
    e_items = [c for c in checks if str(c.get("outcome")) == "fail" and str(c.get("check")) in explained]
    noted = next((c for c in matched if any(entry_explains(entry, item) for entry in explains(c)
                                            for item in e_items)), matched[0])
    return outcome, explained, noted.get("id")


def explained_items(checks, entries):
    """The indices of a line's check items an answered / clamped case's explains_checks entries explain (B5b-06:
    per ITEM, as async's runner's explained_outcome keys E: an item of outcome fail that an entry explains,
    entry_explains: a name every item of that check, {check, detail_matches} an item whose detail it matches)."""
    return [i for i, c in enumerate(checks or []) if str(c.get("outcome")) == "fail"
            and any(entry_explains(entry, c) for entry in entries or [])]


def line_explained_items(line, matched):
    """The explained items of an MCP line an answered / clamped case matched: the line's recorded explained_items
    (the re-evaluation's), else explained_items against the matched cases' union."""
    rec = line.get("explained_items")
    if isinstance(rec, list) and all(isinstance(i, int) and not isinstance(i, bool) for i in rec):
        return list(rec)
    return explained_items(line.get("checks") or [], union_of(matched))


def case_ids(scenario):
    """The table cases a scenario names in mcp.by_design (an id or a list of ids)."""
    ids = ((scenario or {}).get("mcp") or {}).get("by_design")
    if isinstance(ids, str):
        return [ids]
    if isinstance(ids, list):
        return [i for i in ids if isinstance(i, str)]
    return []


def _transport_line(line):
    """Whether a direct execution has no HTTP answer (a transport error: no status, or 0)."""
    http = (line or {}).get("http")
    return not isinstance(http, int) or isinstance(http, bool) or http == 0


def answer_class(case, hit, direct_lines, run_dir):
    """The class an answer-triggered case gives the MCP executions it matched (|hit|): service-side by design (`by`
    answer), verified when a direct execution matches its direct side. Every repeat of BOTH direct arms is classified
    (B5b-04; |direct_lines| {arm: {repeat: line}}): matched, unknown (a key this
    comparator does not match) or outside. An answer outside the direct side is read by the case's
    answer_match.direct_mismatch (triage-2 B-03, async's table): "disagree": any direct HTTP answer outside it
    disagrees (side direct); "timing": a direct 200 or a direct 503 deadline is unverified (the host answered at
    another time), any other direct answer disagrees; a direct transport error or no direct answer is unverified
    whatever the value. A disagreeing repeat wins over a matching one (a match on one repeat does not excuse another
    repeat's contradiction); without direct_mismatch an outside answer is unverified and a match verifies."""
    am = case["answer_match"]
    side = am.get("direct")
    dm = am.get("direct_mismatch")
    res = {"class": "service_by_design", "by": "answer", "case": case["id"], "case_agreement": case.get("agreement"),
           "service_outcome": (case.get("service_answer") or {}).get("outcome"),
           "repeats": [r.get("repeat") for r in hit], "mcp_verified": True}
    lines = [x for arm in ("direct-local", "direct-remote")
             for _, x in sorted(((direct_lines or {}).get(arm) or {}).items(),
                                key=lambda kv: kv[0] if isinstance(kv[0], int) else -1)]
    matched, unknown, outside = [], False, []
    for x in lines:
        top, slots, texts = direct_error_parts(x, run_dir)
        ok = answer_side_matches(side, x.get("http"), top, texts, slot_codes=slots, known=ANSWER_DIRECT_KEYS)
        if ok is True:
            matched.append(x)
        elif ok is None:
            unknown = True
        else:
            outside.append(x)
    disagreeing = []
    for x in outside:
        err = x.get("error") or {}
        if _transport_line(x) or dm not in DIRECT_MISMATCH_VALUES:
            continue
        if dm == "timing" and (x.get("http") == 200 or (x.get("http") == 503 and err.get("code") == "deadline")):
            continue
        disagreeing.append(x)
    if disagreeing:
        x = disagreeing[0]
        res.update({"class": "disagree", "side": "direct", "direct_verified": False, "direct_mismatch": dm,
                    "direct_detail": "%s repeat %s: http %s (%s), outside the case's direct side %s (direct_mismatch "
                                     "%s)%s" % (x.get("arm"), x.get("repeat"), x.get("http"),
                                                (x.get("error") or {}).get("code"), side, dm,
                                                "; %d other direct repeat(s) match" % len(matched) if matched else "")})
        return res
    if matched:
        res.update(direct_verified=True, verified=True,
                   direct_detail="repeat(s) %s match" % ["%s r%s" % (x.get("arm"), x.get("repeat")) for x in matched])
        return res
    if unknown:
        res.update(direct_verified=None, verified=False,
                   direct_detail="the case's answer_match.direct states %s, which this comparator does not match" % (
                       unknown_keys(side, ANSWER_DIRECT_KEYS)))
        return res
    res.update(direct_verified=None, verified=False,
               direct_detail="no direct execution answered %s%s" % (
                   side, " (direct_mismatch %s: unverified)" % dm if dm else ""))
    return res


def _direct_ref(direct_lines):
    for arm in ("direct-local", "direct-remote"):
        if direct_lines.get(arm):
            return arm, direct_lines[arm]
    return None, None


def _split_items(checks, indices):
    """(the items at |indices|, the others, the (name, occurrence) keys of the former) of a check list."""
    keys, seen, picked, rest = set(), {}, [], []
    for i, c in enumerate(checks or []):
        name = c.get("check")
        n = seen.get(name, 0)
        seen[name] = n + 1
        if i in indices:
            picked.append(c)
            keys.add((name, n))
        else:
            rest.append(c)
    return picked, rest, keys


def _split_by_keys(checks, keys):
    """(the items whose (name, occurrence) is in |keys|, the others)."""
    seen, picked, rest = {}, [], []
    for c in checks or []:
        name = c.get("check")
        n = seen.get(name, 0)
        seen[name] = n + 1
        (picked if (name, n) in keys else rest).append(c)
    return picked, rest


def check_agreement(reps, ref, ref_arm, explained=None):
    """The per-check agreement of MCP executions with the reference direct arm (DESIGN §1.3; ref: {repeat: line},
    each MCP repeat against the reference's same repeat, else its first): `explained` {repeat: item indices} the
    check ITEMS an answered / clamped case explains on that repeat (E of explained_outcome, B5b-06), compared apart
    with the reference's items of the same check and occurrence: their differences do not count (B-01). The class:
    disagree (a differing check outside E, or one side alone errored), service_by_design (only items of E differ),
    not_comparable (no check compared), agree."""
    explained = explained or {}
    ref0 = ref[min(ref)]
    differing, by_design, unmatched = [], [], set()
    compared = not_exposed = inconclusive = 0
    for r in reps:
        theirs = (ref.get(r.get("repeat"), ref0) or {}).get("checks") or []
        mine = r.get("checks") or []
        e = set(explained.get(r.get("repeat")) or ())
        if e:
            mine_e, mine, keys = _split_items(mine, e)
            theirs_e, theirs = _split_by_keys(theirs, keys)
            d_e, _, _, _, _ = compare_checks(mine_e, theirs_e)
            for x in d_e:
                x["repeat"] = r.get("repeat")
            by_design.extend(d_e)
        d, u, n, ne, ni = compare_checks(mine, theirs)
        compared += n
        not_exposed += ne
        inconclusive += ni
        unmatched.update(str(x) for x in u)
        for x in d:
            x["repeat"] = r.get("repeat")
        differing.extend(d)
    mcp_out = aggregate_outcome([r.get("outcome") for r in reps])
    ref_out = aggregate_outcome([x.get("outcome") for x in ref.values()])
    # Agreement is on the checks both arms evaluated (DESIGN §1.3: a field the MCP does not expose is recorded,
    # not compared), so a direct failure of a check the MCP cannot see is no disagreement; an execution that
    # errored on one side only is. With no check compared at all there is nothing to agree on (review M2).
    one_sided_error = (mcp_out == "error") != (ref_out == "error")
    if differing or one_sided_error:
        cls = "disagree"
    elif by_design:
        cls = "service_by_design"
    elif not compared:
        cls = "not_comparable"
    else:
        cls = "agree"
    out = {"class": cls, "reference": ref_arm, "differing_checks": differing, "mcp_outcome": mcp_out,
           "reference_outcome": ref_out, "checks_compared": compared, "checks_not_exposed": not_exposed,
           "checks_inconclusive": inconclusive, "checks_unmatched": sorted(unmatched)}
    if by_design:
        out["explained_differences"] = by_design
    return out


def _unexplained_direct(direct_lines, sides, run_dir):
    """The direct executions no answer-triggered case of the scenario explains (a direct answer of the case's direct
    side, e.g. the host's own 503 deadline, says nothing about the other MCP executions' pair): {arm: {repeat: line}}."""
    if not sides:
        return dict((k, v) for k, v in direct_lines.items() if v)
    kept = {}
    for arm_name, lines in direct_lines.items():
        left = {}
        for k, x in (lines or {}).items():
            top, slots, texts = direct_error_parts(x, run_dir)
            if not any(answer_side_matches(s, x.get("http"), top, texts, slot_codes=slots,
                                           known=ANSWER_DIRECT_KEYS) is True for s in sides):
                left[k] = x
        if left:
            kept[arm_name] = left
    return kept


def mcp_agreement(mcp, direct_lines, scenario, cases, run_dir):
    """The MCP arm's outcome class against the reference direct arm (DESIGN §1.3).

    Order: the manifest's `mcp.not_implemented`; MCP executions inconclusive by a host state (host_state_line) are
    left out (all of them: not_comparable, B-09); every other execution is matched as async's runner matches it
    (case_matches: the cases the scenario names, else an answer-triggered case, B-08): an answer-triggered case judges
    its direct side (answer_class, direct_mismatch, B-03), a refused / relayed / service-only named case its pair
    (by_design_class), against the direct executions no answer-triggered case explains; the executions an answered /
    clamped case matched, and those no case matched, are compared check by check with the reference direct execution
    (check_agreement), the failing check items the matched answered cases explain apart (B-01, B5b-06:
    service_by_design only when nothing else differs; the case is a note); an unmatched execution the MCP runner
    marked service_by_design is, with no table, an unverified by-design part of its own beside the others' check
    agreement, and with a table judged by its checks like any other (B5b-05, B-06). A case id the scenario names that
    the table lacks, with no named case matched, is a disagreement (side table). A disagreement anywhere is the class
    (the first in execution order, the checks last); else a by-design class, verified only when every part held
    (B-04: other parts that do not agree make it unverified, shown as `other_repeats`); else the checks' class,
    naming the matched answered cases."""
    mcp = mcp or {}
    cases = cases or {}
    reps = [mcp[k] for k in sorted(mcp)]
    if ((scenario or {}).get("mcp") or {}).get("not_implemented"):
        return {"class": "not_implemented_mcp", "by": "manifest", "field": scenario["mcp"]["not_implemented"]}
    named = named_cases(cases, scenario)
    if not reps:
        if any(c.get("agreement") == "not_implemented_mcp" for c in named):
            return {"class": "not_implemented_mcp", "by": "table"}
        return {"class": "missing"}
    host = [r for r in reps if host_state_line(r, run_dir, scenario)]
    judged = [r for r in reps if not any(r is h for h in host)]
    note = {"host_state_repeats": [r.get("repeat") for r in host]} if host else {}
    if not judged:
        return dict(note, **{"class": "not_comparable", "by": "host state",
                             "detail": "every MCP execution is inconclusive by a host state (%s)" % sorted(set(
                                 str((mcp_error(x, run_dir) or {}).get("code")) for x in host))})
    # the runner's matcher, per execution
    groups, order, answered, explained, explained_names, plain = {}, [], {}, {}, {}, []
    undecided = []
    for r in judged:
        ans, terr, doc = mcp_view(r, run_dir, with_doc=True)
        walls = []
        matched, how = case_matches(cases, scenario, ans, terr, r, doc, walls)
        if matched and walls and walls[-1][0] == "undecidable":
            undecided.append(r)   # rule (c): its case's wall rule cannot be decided: neither for nor against
            continue
        if matched and judged_by_checks(matched[0]):
            _, e, noted = answered_verdict(matched, r.get("checks") or [], sc.outcome_from_checks)
            answered[r.get("repeat")] = noted
            explained[r.get("repeat")] = line_explained_items(r, matched)
            explained_names[r.get("repeat")] = e
            plain.append(r)
        elif matched:
            key = (how, matched[0].get("id"))
            if key not in groups:
                order.append(key)
                groups[key] = (matched[0], [])
            groups[key][1].append(r)
        else:
            plain.append(r)
    if undecided:
        note["undecidable_repeats"] = [r.get("repeat") for r in undecided]
        if not groups and not plain:
            return dict(note, **{"class": "not_comparable", "by": "undecidable plausibility",
                                 "detail": "every MCP execution's plausibility is undecidable (rule (c))"})
    sides = [groups[k][0]["answer_match"].get("direct") for k in order if k[0] == "answer-triggered"]
    kept = _unexplained_direct(direct_lines, sides, run_dir)
    k_arm, k_ref = _direct_ref(kept)
    k_ref0 = k_ref[min(k_ref)] if k_ref else None
    _, ref = _direct_ref(direct_lines)
    parts = []
    missing = [i for i in case_ids(scenario) if i not in cases]
    if missing and not answered and not any(k[0] == "named" for k in order):
        parts.append({"class": "disagree", "case": missing[0], "side": "table",
                      "detail": "case %s is not in the by-design table" % missing[0]})
    for key in order:
        case, hit = groups[key]
        if key[0] == "answer-triggered":
            parts.append(answer_class(case, hit, direct_lines, run_dir))
        else:
            parts.append(by_design_class(case, k_ref0, hit, run_dir, scenario))
    rest = None
    if plain:
        unmatched = [r for r in plain if r.get("repeat") not in answered]
        # B5b-05: the runner's marks no case reproduces: with no table each stays an unverified by-design part of
        # its own; with a table it is judged by its checks (B-06); every other execution is compared check by check
        marks = [r for r in unmatched if r.get("outcome") == "service_by_design"] if not cases else []
        compared_lines = [r for r in plain if not any(r is x for x in marks)]
        if any(r.get("outcome") == "not_implemented_mcp" for r in unmatched):
            rest = {"class": "not_implemented_mcp", "by": "mcp runner"}
            marks = []
        else:
            if marks:
                parts.append({"class": "service_by_design", "by": "mcp runner", "verified": False,
                              "repeats": [r.get("repeat") for r in marks]})
            if compared_lines and k_ref is None:
                rest = {"class": "no_direct_reference",
                        "mcp_outcome": aggregate_outcome([r.get("outcome") for r in compared_lines])}
            elif compared_lines:
                rest = check_agreement(compared_lines, k_ref, k_arm, explained)
    if rest is not None:
        if answered:
            rest["by_design_notes"] = dict((str(k), v) for k, v in sorted(answered.items()))
            rest["explained_checks"] = dict((str(k), v) for k, v in sorted(explained_names.items()) if v)
            rest["explained_items"] = dict((str(k), v) for k, v in sorted(explained.items()) if v)
            if rest["class"] == "service_by_design" and rest.get("explained_differences"):
                case = cases.get(answered.get(rest["explained_differences"][0].get("repeat"))) or {}
                d_ok, d_detail, by = _direct_side(case, k_ref0, run_dir)
                rest.update({"by": "explains_checks", "case": case.get("id"), "case_agreement": case.get("agreement"),
                             "service_outcome": (case.get("service_answer") or {}).get("outcome"),
                             "mcp_verified": True, "direct_verified": d_ok, "direct_detail": d_detail,
                             "verified": d_ok is True})
                if by:
                    rest["direct_by"] = by
        parts.append(rest)
    dis = [p for p in parts if p.get("class") == "disagree"]
    nim = [c for c in named if c.get("agreement") == "not_implemented_mcp"]
    if dis and nim:
        # the table states that the request is not expressible via MCP (e.g. bounds.max_query_bp): a difference
        # between the arms is that gap, as before lane S5b (the first named case not disagreeing decided)
        return dict(note, **{"class": "not_implemented_mcp", "by": "table", "case": nim[0].get("id"),
                             "disagreement": dict(dis[0])})
    if dis:
        out = dict(dis[0])
        trig = [p for p in parts if p.get("by") == "answer"]
        if trig and dis[0] is not trig[0]:
            out.update(answer_case=trig[0]["case"], answer_case_repeats=trig[0]["repeats"])
        return dict(note, **out)
    bd = [p for p in parts if p.get("class") == "service_by_design"]
    if bd:
        out = dict(bd[0])
        others = [p for p in parts if p is not bd[0]]
        out["verified"] = all(p.get("verified") is True for p in bd) and all(
            p.get("class") in ("agree", "service_by_design") for p in others)
        if others:
            out["other_repeats"] = others[0] if len(others) == 1 else others
        return dict(note, **out)
    return dict(note, **(parts[0] if parts else {"class": "missing"}))


def _client_side(line):
    """Wall time (without the connection set-up) minus the server's own elapsed time: transfer, HTTP handling
    and (remotely) the network and the front."""
    return _diff(wall_net(line), _server_elapsed(line))


def front_overhead(local, remote):
    """direct-remote minus direct-local for one scenario, on a pair of host-warm executions only (review H1):
    the raw wall difference (connection set-up excluded), the client-side difference (wall minus the server's
    elapsed time on each side) and the server elapsed difference (near 0 when both are warm). None without a
    warm execution on both arms (the scenario's cold execution is reported on its own)."""
    lw = first_with(local, "warm")
    rw = first_with(remote, "warm", prefer=lw.get("repeat") if lw else None)
    if lw is None or rw is None:
        return None
    return {"local_repeat": lw.get("repeat"), "remote_repeat": rw.get("repeat"),
            "wall_ms": _diff(wall_net(rw), wall_net(lw)),
            "net_ms": _diff(_client_side(rw), _client_side(lw)),
            "server_ms": _diff(_server_elapsed(rw), _server_elapsed(lw))}


def mcp_overhead(mcp, local, remote):
    """MCP minus direct (DESIGN §5) on host-warm executions only: the MCP warm execution against each direct
    arm's first warm one, with and without the queue wait, and the service's backend time against both; the
    MCP cold execution's wall time apart."""
    # Only executions that reached the host measure the MCP path; a refusal at submit or a request the MCP
    # cannot express says nothing about it.
    mcp = dict((k, v) for k, v in (mcp or {}).items() if sc.reached_host(v))
    if not mcp:
        return None
    m = first_with(mcp, "warm", prefer=1 if 1 in mcp else None)
    cold = first_with(mcp, "cold")
    out = {"repeat": m.get("repeat") if m else None,
           "cold_wall_ms": _num(cold.get("wall_ms")) if cold else None,
           "phases": (m.get("mcp_phases") or {}) if m else {}}
    if m is None:
        return out
    phases = out["phases"]
    queue = _num(phases.get("queue_ms"))
    backend = _num(phases.get("backend_ms"))
    wall = _num(m.get("wall_ms"))
    for name, lines in (("remote", remote), ("local", local)):
        d = first_with(lines or {}, "warm")
        if d is None:
            continue
        dwall = wall_net(d)
        out["minus_%s_ms" % name] = _diff(wall, dwall)
        out["minus_%s_no_queue_ms" % name] = _diff(wall - queue if wall is not None and queue is not None else None,
                                                   dwall)
        out["backend_minus_%s_ms" % name] = _diff(backend, dwall)
    return out


def cold_warm(lines):
    """The arm's host-cold execution and its next host-warm repeat on the same arm (only the arm that ran a
    scenario first has a cold execution)."""
    cold = first_with(lines, "cold")
    if cold is None:
        return None
    warm = None
    for k in sorted(lines):
        if k > (cold.get("repeat") or 0) and host_cache(lines[k]) == "warm":
            warm = lines[k]
            break
    if warm is None:
        return None
    c, w = wall_net(cold), wall_net(warm)
    if c is None or w is None:
        return None
    return {"cold_ms": c, "warm_ms": w, "delta_ms": w - c, "ratio": (w / c) if c > 0 else None,
            "server_cold_ms": _server_elapsed(cold), "server_warm_ms": _server_elapsed(warm)}


# ---------------------------------------------------------------------------------------------------------------
# relations between scenarios (the manifest's `relations`)
#
# A relation {id, type, a, b, title, args} states how the answers of scenarios a and b relate; it is evaluated per
# arm on the normalised answers of the first repeat of each there, from the raw files, with the outcome pass,
# fail, inconclusive (the pair cannot test the claim) or not_exposed (the MCP arm lacks a field the type reads).
# The semantics are lane S1b's (its README, "Relations"): checks.evaluate_relation of checks.py next to this
# script is the reference implementation and is used whenever checks.py is there. Without it, same_answer is
# evaluated here with the same rule (same HTTP status; not 200: the same error code and text; args.paths pairs
# of dotted paths equal, else the bodies equal without timing and without args.ignore) and the other types read
# not_evaluated.

def load_answer(line, run_dir):
    """The normalised answer of an execution from its raw file (both runners keep it under `answer`)."""
    raw = line.get("raw") if line else None
    if not raw or not run_dir:
        return None
    try:
        obj = sc.read_gz_json(os.path.join(run_dir, raw))
    except (OSError, ValueError, EOFError):
        return None
    ans = obj.get("answer") if isinstance(obj, dict) else None
    return ans if isinstance(ans, dict) else None


def path_values(obj, path):
    """The values at a dotted path; `*` steps into every element of a list (or value of an object), a number
    into that element of a list."""
    cur = [obj]
    for part in str(path).split("."):
        nxt = []
        for c in cur:
            if part == "*":
                if isinstance(c, list):
                    nxt.extend(c)
                elif isinstance(c, dict):
                    nxt.extend(c.values())
            elif isinstance(c, dict) and part in c:
                nxt.append(c[part])
            elif isinstance(c, list) and part.isdigit() and int(part) < len(c):
                nxt.append(c[int(part)])
        cur = nxt
    return cur


def drop_path(obj, path):
    """obj without the value at a dotted path (`*` as in path_values); obj itself is left unchanged."""
    obj = json.loads(json.dumps(obj))
    parts = str(path).split(".")

    def walk(node, i):
        if i == len(parts) - 1:
            if isinstance(node, dict):
                node.pop(parts[i], None)
            elif isinstance(node, list) and parts[i] == "*":
                del node[:]
            return
        nxt = []
        if parts[i] == "*":
            nxt = list(node.values()) if isinstance(node, dict) else (node if isinstance(node, list) else [])
        elif isinstance(node, dict) and parts[i] in node:
            nxt = [node[parts[i]]]
        elif isinstance(node, list) and parts[i].isdigit() and int(parts[i]) < len(node):
            nxt = [node[int(parts[i])]]
        for n in nxt:
            walk(n, i + 1)
    walk(obj, 0)
    return obj


def first_difference(a, b, path=""):
    """The first dotted path where two JSON values differ, or None."""
    if type(a) is not type(b):
        return path or "(top)"
    if isinstance(a, dict):
        for k in sorted(set(a) | set(b), key=str):
            if k not in a or k not in b:
                return (path + "." if path else "") + str(k)
            d = first_difference(a[k], b[k], (path + "." if path else "") + str(k))
            if d:
                return d
        return None
    if isinstance(a, list):
        if len(a) != len(b):
            return (path or "(top)") + " (length %d vs %d)" % (len(a), len(b))
        for i, (x, y) in enumerate(zip(a, b)):
            d = first_difference(x, y, "%s.%d" % (path, i) if path else str(i))
            if d:
                return d
        return None
    return None if a == b else (path or "(top)")


def relation_same_answer(a, b, args):
    """same_answer as lane S1b states it (for runs without checks.py)."""
    if a.get("http") != b.get("http"):
        return "fail", "http %s vs %s" % (a.get("http"), b.get("http"))
    if a.get("http") != 200:
        ea, eb = a.get("error") or {}, b.get("error") or {}
        if (ea.get("code"), ea.get("text")) != (eb.get("code"), eb.get("text")):
            return "fail", "errors differ: %s %r vs %s %r" % (ea.get("code"), (ea.get("text") or "")[:120],
                                                              eb.get("code"), (eb.get("text") or "")[:120])
        return "pass", "the same refusal (http %s, %s)" % (a.get("http"), ea.get("code"))
    ba = sc.digest_view({"body": a.get("body")})["body"]
    bb = sc.digest_view({"body": b.get("body")})["body"]
    pairs = args.get("paths")
    if pairs:
        for pa, pb in pairs:
            va = path_values(ba, pa) if pa else [ba]
            vb = path_values(bb, pb) if pb else [bb]
            if not va or not vb:
                return "fail", "path %s missing" % (pa if not va else pb)
            if va[0] != vb[0]:
                return "fail", "%s differs from %s" % (pa or "(body)", pb or "(body)")
        return "pass", "%d path pair(s) equal" % len(pairs)
    for p in args.get("ignore") or []:
        ba, bb = drop_path(ba, p), drop_path(bb, p)
    diff = first_difference(ba, bb)
    if diff is None:
        return "pass", "the same answer (timing aside)"
    return "fail", "the bodies differ at %s" % diff


BUILT_IN_RELATIONS = {"same_answer": relation_same_answer}
RELATION_TYPES = ("prefix_equals", "zero_equals_absent", "subset", "same_answer")
RELATION_OUTCOMES = ("pass", "fail", "inconclusive", "not_exposed", "not_evaluated")


def load_relation_evaluator(checks_dir=None):
    """checks.evaluate_relation of checks.py in |checks_dir| (default: next to this script; lane S1's reference),
    or None."""
    path = os.path.join(checks_dir or os.path.dirname(os.path.abspath(__file__)), "checks.py")
    if not os.path.exists(path):
        return None
    import importlib.util
    spec = importlib.util.spec_from_file_location("suite_checks_for_relations", path)
    mod = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(mod)
    except Exception:  # a checks.py that does not load leaves the built-in evaluation
        return None
    return getattr(mod, "evaluate_relation", None)


def evaluate_relations(manifest, grouped, run_dir, arms, evaluator="auto"):
    """(relations, summary): every manifest relation per arm that ran both scenarios. evaluator: lane S1's
    checks.evaluate_relation (found next to this script when "auto"), or None for the built-in evaluation."""
    if evaluator == "auto":
        evaluator = load_relation_evaluator()
    rels = (manifest.get("relations") or []) if isinstance(manifest, dict) else []
    scen = dict((x.get("id"), x) for x in manifest.get("scenarios", [])) if isinstance(manifest, dict) else {}
    out, per_arm = [], {}
    for rel in rels:
        if not isinstance(rel, dict):
            continue
        entry = {"id": rel.get("id"), "type": rel.get("type"), "a": rel.get("a"), "b": rel.get("b"),
                 "title": rel.get("title"), "args": rel.get("args") or {}, "arms": {}}
        for arm in arms:
            la_all = (grouped.get(rel.get("a")) or {}).get(arm) or {}
            lb_all = (grouped.get(rel.get("b")) or {}).get(arm) or {}
            if not la_all or not lb_all:
                continue
            la, lb = la_all[min(la_all)], lb_all[min(lb_all)]
            aa, ab = load_answer(la, run_dir), load_answer(lb, run_dir)
            if aa is None or ab is None:
                outcome, detail = "not_evaluated", "a raw answer is missing"
            elif evaluator is not None:
                try:
                    res = evaluator(rel, scen.get(rel.get("a")) or {}, aa, scen.get(rel.get("b")) or {}, ab, arm)
                    outcome, detail = res.get("outcome"), res.get("detail")
                except Exception as exc:  # never take the comparison down
                    outcome, detail = "fail", "evaluate_relation raised %s: %s" % (type(exc).__name__, exc)
            elif rel.get("type") not in RELATION_TYPES:
                outcome, detail = "fail", "unknown relation type %r" % rel.get("type")
            elif rel.get("type") in BUILT_IN_RELATIONS:
                outcome, detail = BUILT_IN_RELATIONS[rel.get("type")](aa, ab, entry["args"])
            else:
                outcome, detail = "not_evaluated", "%s needs lane S1's checks.evaluate_relation (checks.py)" % (
                    rel.get("type"))
            if outcome not in RELATION_OUTCOMES:
                detail = "outcome %r read as a failure: %s" % (outcome, detail)
                outcome = "fail"
            entry["arms"][arm] = {"outcome": outcome, "detail": detail,
                                  "repeats": [la.get("repeat"), lb.get("repeat")]}
            counts = per_arm.setdefault(arm, dict((k, 0) for k in RELATION_OUTCOMES))
            counts[outcome] += 1
        out.append(entry)
    return out, {"total": len(out), "per_arm": per_arm,
                 "implementation": "checks.evaluate_relation" if evaluator is not None else "compare.py built-in"}


# ---------------------------------------------------------------------------------------------------------------
# resources

class SampleIndex(object):
    """Samples sorted by time, for the join of executions to the resource use while they ran."""

    def __init__(self, samples):
        pairs = []
        for s in samples:
            t = sc.parse_utc(s.get("t_utc"))
            if t is not None:
                pairs.append((t, s))
        pairs.sort(key=lambda p: p[0])
        self.times = [p[0] for p in pairs]
        self.samples = [p[1] for p in pairs]

    def window(self, start, end):
        """The samples from start - 1 s to end + 2 s (docker's CPU figure averages over its own window), or
        the nearest one within 5 s when none falls inside."""
        if not self.times or start is None:
            return []
        end = start if end is None else end
        lo = bisect.bisect_left(self.times, start - SAMPLE_BEFORE_S)
        hi = bisect.bisect_right(self.times, end + SAMPLE_AFTER_S)
        if hi > lo:
            return self.samples[lo:hi]
        best, best_d = None, None
        for i in (lo - 1, lo):
            if 0 <= i < len(self.times):
                d = min(abs(self.times[i] - start), abs(self.times[i] - end))
                if best_d is None or d < best_d:
                    best, best_d = self.samples[i], d
        return [best] if best is not None and best_d <= SAMPLE_NEAREST_S else []


def line_window(line):
    start = sc.parse_utc(line.get("started_utc"))
    end = sc.parse_utc(line.get("ended_utc"))
    if end is None and start is not None and _num(line.get("wall_ms")) is not None:
        end = start + line["wall_ms"] / 1000.0
    return start, end


def resources(lines, index):
    """Peak container CPU / memory, mean sra_all CPU and host load over the executions' windows, and the peak
    server memory_bytes."""
    cpu, mem, sra, load, mem_server = [], [], [], [], []
    joined = 0
    for line in lines:
        start, end = line_window(line)
        win = index.window(start, end)
        if win:
            joined += 1
        for s in win:
            cpu.append(s.get("cpu_pct"))
            mem.append(s.get("mem_bytes"))
            sra.append(s.get("sra_all_cpu_pct"))
            load.append(s.get("host_load1"))
        work = (line.get("server") or {}).get("work") or {}
        mem_server.append(work.get("memory_bytes"))
    nums = lambda vals: [v for v in vals if _num(v) is not None]  # noqa: E731
    return {"executions": len(lines), "executions_with_samples": joined,
            "peak_cpu_pct": max(nums(cpu)) if nums(cpu) else None,
            "peak_mem_bytes": max(nums(mem)) if nums(mem) else None,
            "mean_sra_all_cpu_pct": (sum(nums(sra)) / len(nums(sra))) if nums(sra) else None,
            "mean_host_load1": (sum(nums(load)) / len(nums(load))) if nums(load) else None,
            "peak_server_memory_bytes": max(nums(mem_server)) if nums(mem_server) else None}


# ---------------------------------------------------------------------------------------------------------------
# the comparison

def inconclusive_reason(line, scenario, run_dir, cases=None):
    """An inconclusive execution's class (async's runner at e334715b, the line's inconclusive_reason): rule (c)'s
    text (no service-side timing / no threshold known: plausibility undecidable), the host-state code, or "checks"
    (the checks decided); for a line recorded without the field the same, read from the raw record: the runner's
    order, an answer-triggered case of |cases| whose wall rule is undecidable (B5b-09), else the host-state rule;
    None for any other outcome."""
    if line.get("outcome") != "inconclusive":
        return None
    if line.get("inconclusive_reason"):
        return str(line["inconclusive_reason"])
    if line.get("arm") == "mcp":
        ans, terr, doc = mcp_view(line, run_dir, with_doc=True)
        if cases:
            walls = []
            matched, _ = case_matches(cases, scenario, ans, terr, line, doc, walls)
            if matched and walls and walls[-1][0] == "undecidable":
                return walls[-1][1]
            if matched:
                return "checks"
        code, state, detail = sc.host_state_verdict(scenario or {}, line, doc, ans.get("error"), terr)
        if state == "host_state":
            return code
        if state == "undecidable":
            return detail
    return "checks"


def reevaluation_checks_dir(run_dir):
    """The directory of the checks a run's re-evaluated set was judged with (reevaluated.json), when still there."""
    mark = merge_run.reevaluation_mark(run_dir) or {}
    d = (mark.get("checks") or {}).get("dir")
    return d if d and os.path.exists(os.path.join(d, "checks.py")) else None


def compare(run_dir, by_design=None, relation_evaluator="auto", reevaluated=False, warnings=None):
    """compare.json's content for a merged run directory (of its re-evaluated set with |reevaluated|: the merged
    executions.reevaluated.jsonl and run.reevaluated.json); by_design is load_by_design's (cases, info);
    relation_evaluator as evaluate_relations takes it ("auto": checks.py next to this script when present, for
    the re-evaluated set the re-evaluation's checks.py), whose module also reads /traverse answers for the direct
    arms' timing rule; |warnings| (table_for's) are kept in the result."""
    run = sc.read_json(os.path.join(run_dir, sc.reevaluated_name("run.json", reevaluated)), {}) or {}
    rows, bad = sc.read_jsonl(os.path.join(run_dir, sc.reevaluated_name("executions.jsonl", reevaluated)))
    checks_mod = load_checks_module(reevaluation_checks_dir(run_dir) if reevaluated else None)
    if relation_evaluator == "auto":
        relation_evaluator = getattr(checks_mod, "evaluate_relation", None) if checks_mod is not None else None
    samples, _ = sc.read_jsonl(os.path.join(run_dir, "samples.jsonl"))
    manifest = sc.read_json(os.path.join(run_dir, "manifest.json"))
    cases, table_info = by_design if by_design else ({}, None)
    titles, areas, manifest_order, scen_by_id = {}, {}, [], {}
    if isinstance(manifest, dict):
        for scen in manifest.get("scenarios", []):
            titles[scen.get("id")] = scen.get("title")
            areas[scen.get("id")] = scen.get("area")
            scen_by_id[scen.get("id")] = scen
            manifest_order.append(scen.get("id"))
    grouped, first_seen = group_lines(rows)
    pos = dict((sid, i) for i, sid in enumerate(manifest_order))
    seen_pos = dict((sid, i) for i, sid in enumerate(first_seen))
    order = sorted(first_seen, key=lambda sid: (pos.get(sid, len(pos)), seen_pos[sid]))
    arms_present = sorted(set(r.get("arm") for r in rows), key=lambda a: (merge_run.ARM_ORDER.get(a, 9), a))

    scenarios = {}
    per_area = {}
    timing_vals = {}
    front_vals = {"wall_ms": [], "net_ms": [], "server_ms": []}
    front_area = {}
    front_area_net = {}
    front_unpaired = 0
    mcp_vals = {}
    mcp_cold_walls = []
    mcp_phase_vals = {}
    cw_vals = {}
    connect_vals = {}
    cold_vals = {}
    direct_summary = {"compared": 0, "identical": 0, "same_checks": 0, "not_comparable": 0, "disagree": 0,
                      "timing": 0, "disagreeing": [], "timing_scenarios": [], "checks_compared": 0}
    mcp_summary = {"compared": 0, "checks_compared": 0, "checks_not_exposed": 0, "checks_inconclusive": 0,
                   "by_design_verified": 0, "by_design_unverified": 0}
    lines_by_area_arm = {}

    for sid in order:
        by_arm = grouped[sid]
        any_line = next(iter(next(iter(by_arm.values())).values()))
        area = areas.get(sid) or any_line.get("area") or sid.split("-", 1)[0]
        entry = {"area": area, "title": titles.get(sid), "arms": {}}
        for arm, lines in by_arm.items():
            view = arm_view(lines)
            view["host_cache"] = [host_cache(lines[k]) for k in sorted(lines)]
            entry["arms"][arm] = view
            counts = per_area.setdefault(area, {}).setdefault(arm, dict((k, 0) for k in OUTCOME_KEYS + ("total",)))
            counts[view["outcome"] if view["outcome"] in counts else "mixed"] += 1
            counts["total"] += 1
            tv = timing_vals.setdefault(area, {}).setdefault(arm, {"all": [], "server": [], "by_host_cache": {}})
            for rep, line in lines.items():
                hc = host_cache(line)
                wall = wall_net(line)
                tv["all"].append(wall)
                tv["server"].append(_server_elapsed(line))
                b = tv["by_host_cache"].setdefault(hc, {"wall": [], "server": []})
                b["wall"].append(wall)
                b["server"].append(_server_elapsed(line))
                lines_by_area_arm.setdefault(area, {}).setdefault(arm, []).append(line)
                connect = _num((line.get("client") or {}).get("connect_ms"))
                if connect is not None:
                    connect_vals.setdefault(arm, []).append(connect)
                if hc == "cold":
                    cv = cold_vals.setdefault(arm, {"wall_ms": [], "server_ms": []})
                    cv["wall_ms"].append(wall)
                    cv["server_ms"].append(_server_elapsed(line))
            cw = cold_warm(lines)
            if cw:
                entry["arms"][arm]["cold_warm"] = cw
                cw_vals.setdefault(arm, {"delta_ms": [], "ratio": []})
                cw_vals[arm]["delta_ms"].append(cw["delta_ms"])
                cw_vals[arm]["ratio"].append(cw["ratio"])
        local, remote, mcp = by_arm.get("direct-local"), by_arm.get("direct-remote"), by_arm.get("mcp")
        if local and remote:
            agree = direct_agreement(local, remote, (scen_by_id.get(sid) or {}).get("route") or any_line.get("route"),
                                     run_dir, checks_mod, scen_by_id.get(sid))
            entry["direct_agreement"] = agree
            direct_summary["compared"] += 1
            direct_summary[agree["class"]] += 1
            direct_summary["checks_compared"] += agree["checks_compared"]
            if agree["class"] == "disagree":
                direct_summary["disagreeing"].append(sid)
            elif agree["class"] == "timing":
                direct_summary["timing_scenarios"].append(sid)
            fo = front_overhead(local, remote)
            entry["front_overhead"] = fo
            if fo is None:
                front_unpaired += 1
            else:
                for k in front_vals:
                    if fo.get(k) is not None:
                        front_vals[k].append(fo[k])
                if fo.get("wall_ms") is not None:
                    front_area.setdefault(area, []).append(fo["wall_ms"])
                if fo.get("net_ms") is not None:
                    front_area_net.setdefault(area, []).append(fo["net_ms"])
        if "mcp" in arms_present:
            ma = mcp_agreement(mcp, {"direct-local": local, "direct-remote": remote}, scen_by_id.get(sid), cases,
                               run_dir)
            entry["mcp_agreement"] = ma
            mcp_summary["compared"] += 1
            mcp_summary[ma["class"]] = mcp_summary.get(ma["class"], 0) + 1
            mcp_summary["checks_compared"] += ma.get("checks_compared") or 0
            mcp_summary["checks_not_exposed"] += ma.get("checks_not_exposed") or 0
            mcp_summary["checks_inconclusive"] += ma.get("checks_inconclusive") or 0
            if ma["class"] == "service_by_design" and "verified" in ma:
                mcp_summary["by_design_verified" if ma["verified"] else "by_design_unverified"] += 1
            elif ma["class"] == "service_by_design":
                mcp_summary["by_design_unverified"] += 1  # the MCP runner's own mark: nothing checked here
            if ma["class"] in ("disagree", "not_comparable"):
                mcp_summary.setdefault(ma["class"] + "_scenarios", []).append(sid)
            mo = mcp_overhead(mcp, local or {}, remote or {})
            if mo:
                entry["mcp_overhead"] = mo
                for k, v in mo.items():
                    if (k.startswith("minus_") or k.startswith("backend_minus_")) and _num(v) is not None:
                        mcp_vals.setdefault(k, []).append(v)
                if _num(mo.get("cold_wall_ms")) is not None:
                    mcp_cold_walls.append(mo["cold_wall_ms"])
                for k, v in (mo.get("phases") or {}).items():
                    if _num(v) is not None:
                        mcp_phase_vals.setdefault(k, []).append(v)
        scenarios[sid] = entry
    if "disagree_scenarios" in mcp_summary:
        mcp_summary["disagreeing"] = mcp_summary.pop("disagree_scenarios")

    relations, relations_summary = evaluate_relations(manifest, grouped, run_dir, arms_present,
                                                      evaluator=relation_evaluator)

    index = SampleIndex(samples)
    resources_out = {}
    for area, by_arm in lines_by_area_arm.items():
        for arm, lines in by_arm.items():
            resources_out.setdefault(area, {})[arm] = resources(lines, index)
    by_reason = {}
    for r in rows:
        why = inconclusive_reason(r, scen_by_id.get(r.get("scenario")), run_dir, cases)
        if why is not None:
            arm_counts = by_reason.setdefault(r.get("arm"), {})
            arm_counts[why] = arm_counts.get(why, 0) + 1
    totals = {}
    for arm in arms_present:
        arm_lines = [r for r in rows if r.get("arm") == arm]
        totals[arm] = {"executions": len(arm_lines),
                       "outcomes": dict((k, sum(1 for r in arm_lines if r.get("outcome") == k)) for k in OUTCOME_KEYS),
                       "host_cache": dict((k, sum(1 for r in arm_lines if host_cache(r) == k))
                                          for k in sc.HOST_CACHE_VALUES),
                       "wall_ms": sc.summarise([wall_net(r) for r in arm_lines]),
                       "resources": resources(arm_lines, index)}

    timing = {}
    for area, by_arm in timing_vals.items():
        for arm, tv in by_arm.items():
            timing.setdefault(area, {})[arm] = {
                "all": sc.summarise(tv["all"]), "server": sc.summarise(tv["server"]),
                "by_host_cache": dict((hc, {"wall": sc.summarise(b["wall"]), "server": sc.summarise(b["server"])})
                                      for hc, b in tv["by_host_cache"].items())}

    return {
        "run_id": run.get("run_id"), "generated_utc": sc.iso_utc(time.time()),
        "release": run.get("release"), "server_instance": run.get("server_instance"),
        "manifest_sha256": (run.get("manifest") or {}).get("computed_sha256")
        or (run.get("manifest") or {}).get("content_sha256"),
        "arms": arms_present, "unparsable_lines": bad, "by_design_table": table_info,
        "warnings": list(warnings or []),
        "host_cache": run.get("host_cache"), "reevaluation": run.get("reevaluation") if reevaluated else None,
        "order": order, "scenarios": scenarios, "relations": relations,
        "summary": {
            "relations": relations_summary,
            "per_area": per_area, "totals": totals,
            "direct_agreement": direct_summary, "mcp_agreement": mcp_summary,
            "inconclusive_by_reason": dict((a, dict(sorted(v.items(), key=lambda kv: (-kv[1], kv[0]))))
                                           for a, v in sorted(by_reason.items())),
            "timing": timing,
            "front_overhead": {"warm": dict((k, sc.summarise(v)) for k, v in front_vals.items()),
                               "per_area_wall": dict((a, sc.summarise(v)) for a, v in front_area.items()),
                               "per_area_net": dict((a, sc.summarise(v)) for a, v in front_area_net.items()),
                               "scenarios_without_warm_pair": front_unpaired,
                               "connect_ms": dict((a, sc.summarise(v)) for a, v in connect_vals.items()),
                               "cold": dict((a, dict((k, sc.summarise(v)) for k, v in cv.items()))
                                            for a, cv in cold_vals.items())},
            "mcp_overhead": dict((k, sc.summarise(v)) for k, v in mcp_vals.items()),
            "mcp_cold_wall_ms": sc.summarise(mcp_cold_walls),
            "mcp_phases": dict((k, sc.summarise(v)) for k, v in mcp_phase_vals.items()),
            "cold_warm": dict((arm, {"n": len(v["delta_ms"]), "delta_ms": sc.summarise(v["delta_ms"]),
                                     "ratio": sc.summarise(v["ratio"])}) for arm, v in cw_vals.items()),
            "resources": resources_out,
        },
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description="Join a run's arms per scenario into compare.json.")
    ap.add_argument("run_dir")
    ap.add_argument("--by-design", default=None,
                    help="async's service-by-design table (default: service_by_design.json here, for --reevaluated "
                         "the table the re-evaluation replayed; '' for none)")
    ap.add_argument("--no-merge", action="store_true", help="do not merge RUN_DIR/arms first")
    ap.add_argument("--reevaluated", action="store_true",
                    help="the run's re-evaluated set (reevaluate.py): compare.reevaluated.json")
    args = ap.parse_args(argv)
    if not args.no_merge and os.path.isdir(os.path.join(args.run_dir, "arms")):
        merge_run.merge(args.run_dir, args.reevaluated)
    path, warnings = table_for(args.run_dir, args.by_design, args.reevaluated)
    for w in warnings:
        sys.stderr.write("compare: warning: %s\n" % w)
    try:
        table = load_by_design(path)
    except TableError as exc:
        sys.stderr.write("compare: refusing: %s\n" % exc)
        return 4
    result = compare(args.run_dir, table, reevaluated=args.reevaluated, warnings=warnings)
    sc.write_json_atomic(os.path.join(args.run_dir, sc.reevaluated_name("compare.json", args.reevaluated)), result)
    s = result["summary"]
    sys.stdout.write("compare %s: %d scenarios, arms %s; direct agreement %s; mcp %s\n" % (
        args.run_dir, len(result["order"]), ", ".join(result["arms"]),
        dict((k, v) for k, v in s["direct_agreement"].items() if not isinstance(v, list)),
        dict((k, v) for k, v in s["mcp_agreement"].items() if not isinstance(v, list))))
    return 0


if __name__ == "__main__":
    sys.exit(main())
