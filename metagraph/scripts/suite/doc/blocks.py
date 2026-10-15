"""Expectations in a few words and the detail blocks of the tests that are not clean (DESIGN-DOC §4.4, §4.6).

`expect_text` turns one expectation (check kind + args) into a short phrase: one formatter per check kind, an
unknown kind printing `name(args)`. `expected_text` lines up a scenario's expectations: the status first, the
property checks, the generic checks last. `detail_block` writes the block of one test from its arm states and the
execution lines.
"""

import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

if os.path.dirname(HERE) not in sys.path:
    sys.path.insert(0, os.path.dirname(HERE))

import cells  # noqa: E402
import latex  # noqa: E402
import suite_common as sc  # noqa: E402
from latex import esc, tt, inline  # noqa: E402

GENERIC_LAST = ("invariants", "budget_respected", "counts_consistent", "deterministic", "content_encoding",
                "answer_fields", "entry_fields")


def _runs(v):
    try:
        return json.dumps(v, separators=(",", ":"))
    except (TypeError, ValueError):
        return str(v)


def _status(args):
    http = args.get("http")
    parts = []
    if isinstance(http, list):
        parts.append("http %s" % " or ".join(str(h) for h in http))
    elif http is not None:
        parts.append("http %s" % http)
    if args.get("code"):
        parts.append("code %s" % tt(args["code"]))
    if args.get("text") is not None:
        parts.append(tt(args["text"], 90))
    if args.get("text_prefix") is not None:
        parts.append("prefix %s" % tt(args["text_prefix"], 90))
    if args.get("text_contains") is not None:
        parts.append("contains %s" % tt(args["text_contains"], 90))
    return " ".join(parts) if parts else "status"


def _resolve_runs(args):
    label = tt(args.get("label", "?"))
    if args.get("full"):
        text = "full under %s" % label
    elif args.get("runs") is not None:
        text = "runs %s under %s" % (tt(_runs(args["runs"])), label)
    elif args.get("kmers") is not None:
        text = "%s k-mer%s under %s" % (args["kmers"], "" if args["kmers"] == 1 else "s", label)
    else:
        text = "runs under %s" % label
    if args.get("prefix_if_stopped"):
        text += ", prefix if stopped"
    if args.get("unless_stopped"):
        text += " unless stopped"
    return text


def _clamped(args):
    return "clamped %s %s {\\mdseries\\upshape\\glyphfont →} %s" % (tt(args.get("field", "?")), esc(args.get("requested")),
                                                   esc(args.get("effective")))


def _budget(args):
    text = "budget respected"
    if args.get("cap_ms") is not None:
        text += " (cap %s)" % cells.fmt_t(args["cap_ms"])
    return text


def _answer_fields(args):
    """`answer_fields`: the fields an answer must carry, equal, be without (dotted paths), on a status when
    `when_http` names one: `results.1.duplicate = true`, `present usage`, `absent usage (on http 409)`."""
    parts = []
    eq = args.get("equals")
    if isinstance(eq, dict) and eq:
        parts.append(", ".join("%s = %s" % (tt(k), tt(_runs(v))) for k, v in eq.items()))
    if args.get("present"):
        parts.append("present %s" % ", ".join(tt(p) for p in args["present"]))
    if args.get("absent"):
        parts.append("absent %s" % ", ".join(tt(p) for p in args["absent"]))
    text = "; ".join(parts) if parts else "fields as expected"
    if args.get("when_http") is not None:
        text += " (on http %s)" % esc(args["when_http"])
    return text


_GENERIC_BREAK = "(,=_.:/\"[{"


def _generic(name, args):
    """`name(key=value, ...)` for a check kind without a formatter, escaped, with a break opportunity after each
    ( , = _ . : / " [ { so that the text wraps inside a narrow expectation column instead of running into the next
    one (the full document, 2026-10-15); each value is cut at 40 characters."""
    text = name if not args else "%s(%s)" % (name, ", ".join(
        "%s=%s" % (k, latex.cut_text(_runs(v), 40)) for k, v in sorted(args.items())))
    out = []
    for ch in text:
        out.append(esc(ch))
        if ch in _GENERIC_BREAK:
            out.append("\\allowbreak{}")
    return "".join(out)


EXPECT_TEXT = {
    "status": _status,
    "invariants": lambda a: "invariants",
    "resolve_runs": _resolve_runs,
    "resolve_graph_runs": lambda a: "graph\\_runs %s" % tt(_runs(a.get("runs"))),
    "label_not_full": lambda a: "not full under %s" % tt(a.get("label", "?")),
    "resolve_absent": lambda a: "no k-mer in the graph",
    "resolve_discover_contains": lambda a: "discovery lists %s" % tt(a.get("label", "?")),
    "clamped": _clamped,
    "budget_respected": _budget,
    "deterministic": lambda a: "deterministic (the digest of the earlier repeats)",
    "counts_consistent": lambda a: "counts consistent",
    "answer_fields": _answer_fields,
}


def expect_text(check, args):
    """One expectation in a few words (LaTeX)."""
    args = args if isinstance(args, dict) else {}
    fn = EXPECT_TEXT.get(check)
    if fn is None:
        return _generic(check, args)
    try:
        return fn(args)
    except Exception:  # an argument shape the formatter does not know: the generic form, never a crash
        return _generic(check, args)


def ordered_expectations(scenario, is_property=None):
    """The expectations: status first, property checks, generic checks last."""
    expect = [e for e in (scenario.get("expect") or []) if isinstance(e, dict)]
    status = [e for e in expect if e.get("check") == "status"]
    rest = [e for e in expect if e.get("check") != "status"]
    generic = [e for e in rest if e.get("check") in GENERIC_LAST
               or (is_property is not None and not is_property(e.get("check")))]
    props = [e for e in rest if e not in generic]
    return status + props + generic


def expected_text(scenario, is_property=None, limit=None):
    parts = [expect_text(e.get("check"), e.get("args")) for e in ordered_expectations(scenario, is_property)]
    text = " · ".join(parts)
    return text


def refusal_expected(scenario, previous=None, is_property=None):
    """The expected column of a refusals row: the status expectation; 'the same text' when it repeats the
    previous row's; for an accepted input (a 2xx status) the property expectations the row carries follow the
    status, so that a ● there shows its claim (VERIFY-T M7: T-edge-0008's `results.1.duplicate = true`, T-edge-0012's
    clamped budget)."""
    for e in scenario.get("expect") or []:
        if isinstance(e, dict) and e.get("check") == "status":
            args = e.get("args") or {}
            text = args.get("text") or args.get("text_prefix") or args.get("text_contains")
            if previous is not None and text and text == previous:
                head = "http %s" % args.get("http") if args.get("http") is not None else ""
                return (head + " the same text").strip(), text
            head = _status(args)
            if not host_expects_refusal(scenario):
                props = [x for x in ordered_expectations(scenario, is_property)
                         if x.get("check") != "status" and (is_property is None or is_property(x.get("check")))]
                if props:
                    head += " · " + " · ".join(expect_text(x.get("check"), x.get("args")) for x in props)
            return head, text
    return "–", None


# -- the MCP arm's standing on an edge input ---------------------------------------------------------------------

MCP_CLASS_WORD = {"not_expressible": "not expressible through MCP",
                  "refused_before_host": "refused before host execution",
                  "service_answered": "answered by the service in the host's place",
                  "normalised": "answered after service normalization",
                  "relayed": "host refusal relayed through the service",
                  "answered": "answered by the host through the service",
                  "refused_unexplained": "refused by the service (no by-design case explains it)",
                  "host_state": "ended by the service's own limit (a host state)",
                  "marked": "every property check marked as not shown via MCP (a scenario-level mark)",
                  "none": "no MCP execution recorded"}

# how a detail block names a check item a property-level mark of the manifest sets aside, or one the by-design
# table explains (suite_common.MCP_CHECK_MARKS; compare's explained_items): (glyph, words)
MARK_GROUP = {"not_implemented": ("\\gni", "not implemented via MCP (the manifest's property-level mark)"),
              "by_design": ("\\gbd", "service-side by design (the manifest's property-level mark)"),
              "explained": ("\\gbd", "explained by the by-design table's case (set aside by the comparison and "
                                     "the gate)")}


def scenario_level_mark(line):
    """Whether an MCP execution is inconclusive by a scenario-level mark (suite_common.SCENARIO_LEVEL_MARK: the
    manifest's property-level marks cover every property check of the execution; the checks decided, not a host
    state)."""
    return line.get("outcome") == "inconclusive" and line.get("inconclusive_reason") == sc.SCENARIO_LEVEL_MARK


def host_expects_refusal(scenario):
    """Whether the scenario's status expectation is a refusal (an HTTP status other than 2xx)."""
    for e in scenario.get("expect") or []:
        if isinstance(e, dict) and e.get("check") == "status":
            http = (e.get("args") or {}).get("http")
            codes = http if isinstance(http, list) else [http]
            return any(isinstance(c, int) and c >= 400 for c in codes)
    return False


def mcp_class(lines, scenario, agreement, cases):
    """The class of the MCP arm's outcome on one scenario, derived from the recorded evidence, never from prose
    (the owner's review, point 5): whether the host was reached (`host_cache`), the error's source (the by-design
    case's `service_answer.outcome`: refused_at_submit, relayed, answered, not_expressible) and the status.

    Returns (class key, the evidence in a few words, LaTeX). Keys: not_expressible (no answer recorded, the manifest's
    mark or the table's not-expressible case), refused_before_host (the service answered a refusal without
    calling the host), service_answered (the service answered a 2xx from its own data without calling the host:
    the capabilities tools' own documents; VERIFY-QCND H5), relayed (the service called the host and relayed its
    refusal under its own code), normalised (the service reshaped the request and the host answered where it would
    have refused the original), answered (the host answered through the service), host_state (the service ended
    the call by its own limit), refused_unexplained (a refusal that failed: no case explains it), none (no
    execution)."""
    lines = [l for l in lines if isinstance(l, dict)]
    agreement = agreement or {}
    if not lines:
        return "none", "no MCP execution recorded"
    # a recorded not-implemented mark is no answer, whatever status a runner wrote beside it
    answered = [l for l in lines if l.get("http") is not None and l.get("outcome") != "not_implemented_mcp"]
    case_id = next((l.get("service_by_design_case") for l in lines if l.get("service_by_design_case")), None) \
        or agreement.get("case")
    case = (cases or {}).get(case_id) or {}
    sa = case.get("service_answer") or {}
    service_outcome = sa.get("outcome")
    if not answered:
        reason = (scenario.get("mcp") or {}).get("not_implemented") or sa.get("missing_argument")
        return "not_expressible", ("no MCP call made: %s" % tt(reason)) if reason else "no MCP call made"
    first = answered[0]
    http = first.get("http")
    code = (first.get("error") or {}).get("code") if isinstance(first.get("error"), dict) else None
    reached = any(cells.host_cache_of(l) in ("cold", "warm", cells.WARM_AFTER_INTERRUPTED, "aged") for l in answered)
    status_text = "http %s%s" % (esc(http), (" " + tt(code)) if code else "")
    if scenario_level_mark(first):
        return "marked", "%s; %s" % (status_text, esc(first.get("verdict_detail") or sc.SCENARIO_LEVEL_MARK, 200))
    if first.get("inconclusive_reason") and first.get("inconclusive_reason") != "checks":
        return "host_state", "%s; the service's limit %s" % (status_text, tt(first["inconclusive_reason"]))
    failed = first.get("outcome") in ("fail", "error")
    if http and http >= 400 and failed:
        # a refusal no by-design case explains: a failure of the run, never one of the four by-design classes
        return "refused_unexplained", "%s; no by-design case matched" % status_text
    if not reached:
        # the host never saw the request: the service answered at submit, in the host's place (a relayed case
        # is the one exception: the service called the host and relays its refusal under its own code)
        if service_outcome == "relayed":
            return "relayed", "%s; the host was called and its refusal relayed" % status_text
        if isinstance(http, int) and http < 400 and not code and service_outcome != "refused_at_submit":
            # a 2xx without an error code from the service's own data (the capabilities tools answer the service's
            # view of the host; the table's service outcome is `answered`): no refusal, and nothing of the host's
            # answer to compare
            return "service_answered", "%s from the service's own data, the host not reached" % status_text
        if service_outcome == "not_expressible" or agreement.get("class") == "not_implemented_mcp":
            return "refused_before_host", "%s at submit (the table: not expressible)" % status_text
        return "refused_before_host", "%s at submit" % status_text
    if http and http >= 400:
        return "relayed", "%s; the host was called and its refusal relayed" % status_text
    if host_expects_refusal(scenario) and (first.get("outcome") in ("service_by_design", "pass")
                                           or service_outcome == "answered"):
        # the host answered where the direct request is refused: the service reshaped it (de-duplicated the
        # labels, upper-cased the sequence) before calling the host, as its table documents; a failed execution
        # against an expected refusal is no normalisation but a disagreement, and reads answered
        return "normalised", "%s after the service reshaped the request (the host expectation is a refusal)" % status_text
    return "answered", status_text


# -- detail blocks ---------------------------------------------------------------------------------------------

HC_WORD = {"cold": "host-cold", "warm": "host-warm", cells.WARM_AFTER_INTERRUPTED: "warm after an interruption",
           "aged": "aged", "not_reached": "host not reached", "unknown": "host cache unknown"}
HC_SHORT = {"cold": "cold", "warm": "warm", cells.WARM_AFTER_INTERRUPTED: "warm, i", "aged": "aged",
            "not_reached": "not reached", "unknown": "unknown"}


def _err(line):
    err = line.get("error") if isinstance(line.get("error"), dict) else {}
    return err.get("code"), err.get("text")


LONG_CODE = 25


def _status_cell(line):
    """`503 deadline`; a code of more than LONG_CODE characters (`traversal_resolve_deadline`, 26) in a top-aligned
    box that wraps at its underscores, so that the block's table stays within the margin (VERIFY-QCND L7); the
    shorter codes of the frozen chapter R sample (`resolve_query_too_long`) print as before."""
    http = line.get("http")
    code, _ = _err(line)
    if http is None:
        return "–"
    if code and len(str(code)) > LONG_CODE:
        return "\\parbox[t]{27mm}{\\raggedright\\strut %s %s}" % (esc(http), tt(code))
    return "%s%s" % (esc(http), (" " + tt(code)) if code else "")


def _budget_text(run, sid):
    b = run.budget_ms(sid)
    return ("{:,} ms".format(int(b)) if isinstance(b, (int, float)) else "none")


def _property_cell(line, is_property, aside=None):
    """The 'property tested?' cell: what the execution tested of its property checks; the items a property-level
    mark of the manifest sets aside (|aside|, suite_common.checks_marked_aside) are counted apart."""
    word, passed, total = cells.property_tested(line, is_property, aside)
    marked = len(aside or ())
    suffix = ("; %d marked not shown via MCP" % marked) if marked else ""
    if word == "n/a":
        return ("n/a (status is the claim)" if total else "n/a") + suffix
    if word == "yes":
        return "yes (%d of %d)%s" % (passed, total, suffix)
    if word == "partly":
        return "partly (%d of %d)%s" % (passed, total, suffix)
    return "no" + suffix


def _finding_checks(line):
    return [c for c in (line.get("checks") or []) if isinstance(c, dict) and c.get("outcome") == "finding"]


def _split_clause(text):
    """(the evidence clause, the explanation) of a check detail: split at the first ': ' outside any parenthesis
    or bracket (a ': ' inside one belongs to the clause, as in `expected 200 (error x: ...)`); (text, '') when
    there is none."""
    text = str(text or "")
    depth = 0
    for i, ch in enumerate(text):
        if ch in "([":
            depth += 1
        elif ch in ")]":
            depth = max(0, depth - 1)
        elif ch == ":" and depth == 0 and text[i + 1:i + 2] == " ":
            return text[:i], text[i + 2:]
    return text, ""


def _balanced(text):
    """`text` with every parenthesis a cut left open closed after it (VERIFY-R2 L2)."""
    opened = text.count("(") - text.count(")")
    return text + ")" * opened if opened > 0 else text


def _first_clause(text, limit=200):
    """The evidence clause of a check detail: up to its first ': ' outside a parenthesis (the explanation
    follows), cut at `limit`, the parentheses balanced."""
    head, _ = _split_clause(text)
    return _balanced(latex.cut_text(head, limit))


def _overridden(line, is_mcp):
    """The MCP runner's host-state rule: the checks stay as they read, the outcome is overridden. A scenario-level
    mark (scenario_level_mark) is no host state: the checks decided, every property check of them marked."""
    return bool(is_mcp and line.get("outcome") == "inconclusive" and line.get("inconclusive_reason")
                and line.get("inconclusive_reason") != "checks" and not scenario_level_mark(line))


def adjudication(line, run, sid, is_mcp, state=None):
    """One phrase: the final outcome of an execution and why, the classification first (the owner's review:
    'timeout classified as a host-state outcome; the sequence property is untested')."""
    outcome = line.get("outcome")
    http = line.get("http")
    code, text = _err(line)
    wall = cells.fmt_t(line.get("wall_ms"))
    budget = _budget_text(run, sid)
    if outcome == "pass":
        return "pass: every check passed"
    if outcome == "finding":
        fc = _finding_checks(line)
        if http == 503 and code == "deadline":
            return ("\\textbf{FINDING}: an allowed 503 deadline, but delivered later than max(1 s, the budget) after "
                    "the %s budget (the evidence and the rule below); the property itself was not tested" % esc(budget))
        return "\\textbf{FINDING}: %s (the evidence and the rule below); the property itself was not tested" % esc(
            _first_clause(fc[0].get("detail")) if fc else "a stated bound broken")
    if outcome == "inconclusive":
        if _overridden(line, is_mcp):
            return ("inconclusive: the service's %s (http %s) classified as a host-state outcome by the runner's "
                    "rule (below), not as a failure; the property is untested on this repeat" % (
                        tt(line["inconclusive_reason"]), esc(http)))
        if is_mcp and scenario_level_mark(line):
            return ("inconclusive: every property check of this execution is marked in the manifest as not shown "
                    "via MCP (%s), so nothing of the claim was tested through the tools: a scenario-level mark, "
                    "never a pass on the generic checks alone" % esc(
                        _first_clause(line.get("verdict_detail") or sc.SCENARIO_LEVEL_MARK, 200)))
        if http == 503 and code == "deadline":
            return ("inconclusive: an allowed 503 deadline at %s against a %s budget, within max(1 s, the budget) "
                    "of it; nothing in the body to judge, the property not tested" % (wall, esc(budget)))
        if http is not None and http != 200:
            return "inconclusive: http %s%s, an answer the scenario allows; the property not tested" % (
                esc(http), (" " + tt(code)) if code else "")
        checks = [c for c in (line.get("checks") or []) if isinstance(c, dict) and c.get("outcome") == "inconclusive"]
        if checks:
            return "inconclusive: %s (%s)" % (esc(checks[0].get("check")), esc(_first_clause(checks[0].get("detail"), 160)))
        return "inconclusive: the property could not be tested"
    if outcome == "fail":
        failed = [c for c in (line.get("checks") or []) if isinstance(c, dict) and c.get("outcome") == "fail"]
        if failed:
            return "\\textbf{FAIL}: %s: %s" % (esc(failed[0].get("check")), esc(failed[0].get("detail"), 200))
        return "\\textbf{FAIL}"
    if outcome == "error":
        return "\\textbf{ERROR}: no HTTP answer%s" % ((": " + esc(text, 160)) if text else "")
    if outcome == "not_implemented_mcp":
        return "not exercised via MCP"
    if outcome == "service_by_design":
        return "service-side by design"
    return esc(cells.WORD.get(outcome, str(outcome)))


def _arm_adjudication(st, run, sid, arm):
    """The adjudication of one arm, or None when the execution rows already say it (the owner's review, point 1:
    one or two lines per arm for an ordinary 503): a phrase for a failure, an error or a host-state override, per
    repeat when the repeats differ; a planned execution not recorded; a disagreement with the reference arm. A
    pass, an allowed 503, a finding (its line follows) and the direct classes (the classes line) add nothing."""
    if not st.present:
        return "no execution recorded although the run planned one" if st.missing_planned else None
    phrases = []
    for l in st.lines:
        outcome = l.get("outcome")
        if outcome in ("fail", "error") or _overridden(l, arm == "mcp") or (arm == "mcp" and scenario_level_mark(l)):
            phrases.append((l.get("repeat"), adjudication(l, run, sid, arm == "mcp", st)))
    if not phrases:
        text = None
    elif len(phrases) == 1 and len(st.lines) == 1:
        text = phrases[0][1]
    elif len(set(p for _, p in phrases)) == 1 and len(phrases) == len(st.lines):
        text = "both repeats: %s" % phrases[0][1] if len(phrases) == 2 else "every repeat: %s" % phrases[0][1]
    else:
        text = "; ".join("r%s %s" % (esc(r), p) for r, p in phrases)
    if st.ne:
        text = (text + " · " if text else "") + "DISAGREES with the reference access path"
    return text


def _direct_detail_words(detail):
    """compare's `direct_detail` of an unverified by-design pair in words: 'no direct execution answered
    {'http': 503, 'code': 'deadline'} (direct_mismatch timing: unverified)' reads 'no direct execution answered
    http 503 deadline (a timing difference, not a disagreement)'; any other detail as recorded."""
    text = str(detail or "")
    m = re.match(r"no direct execution answered \{'http': (\d+)(?:, 'code': '([^']*)')?\}(?: \(direct_mismatch (\w+)[^)]*\))?",
                 text)
    if m:
        out = "no direct execution answered http %s%s" % (m.group(1), (" " + m.group(2)) if m.group(2) else "")
        if m.group(3) == "timing":
            out += " (a timing difference, not a disagreement)"
        elif m.group(3):
            out += " (direct mismatch: %s)" % m.group(3)
        return out
    return text


OUTCOME_PHRASE = {"pass": "a pass", "inconclusive": "an inconclusive execution", "fail": "a failure", "error": "an error",
                  "finding": "a finding"}


def other_repeats(st):
    """The repeats of a by-design MCP arm the case did not match, with their lines (cells.exercised_lines)."""
    return [(l.get("repeat"), l) for l in cells.exercised_lines(st)]


def unverified_reason(ma, st=None):
    """Why a by-design pair is unverified, from compare's mcp_agreement, in plain text: the direct side did not
    match the case (`direct_verified` false or None, with `direct_detail`), or the other MCP repeat outside the
    case (`other_repeats`: a separate host answer) has no direct answer to compare with or disagrees, so the pair
    as a whole is not verified although the case's own repeat is (VERIFY-PLXM H3, VERIFY-QCND H3). None when the
    pair is verified or the class is not by design."""
    ma = ma or {}
    if ma.get("class") != "service_by_design" or ma.get("verified") is not False:
        return None
    parts = []
    if ma.get("direct_verified") is not True:
        dis = ma.get("disagreement") if isinstance(ma.get("disagreement"), dict) else {}
        detail = ma.get("direct_detail") or dis.get("direct_detail")
        parts.append("the direct side did not match the case" + ((": " + _direct_detail_words(detail)) if detail else ""))
    others = ma.get("other_repeats")
    reps = [r for r, _ in other_repeats(st)] if st is not None else []
    which = ("r%s" % ", r".join(str(r) for r in reps)) if reps else "the other repeat"
    for p in ([others] if isinstance(others, dict) else (others or [])):
        if not isinstance(p, dict) or p.get("class") in ("agree", "service_by_design", None):
            continue
        outcome = p.get("mcp_outcome")
        what = OUTCOME_PHRASE.get(outcome, ("a %s" % cells.WORD.get(outcome, str(outcome))) if outcome else "an execution")
        if p.get("class") == "no_direct_reference":
            parts.append("%s outside the case (%s) has no direct answer of the same kind to compare with" % (which, what))
        elif p.get("class") == "disagree":
            parts.append("%s outside the case (%s) disagrees with the reference arm" % (which, what))
        else:
            parts.append("%s outside the case reads %s" % (which, str(p.get("class")).replace("_", " ")))
    return "; ".join(parts) if parts else "one side did not match the case"


def unverified_cause(ma):
    """'direct' (the direct side did not match the case), 'other_repeat' (a repeat outside the case without a
    direct reference or disagreeing), both joined with '+', or None for a verified or not by-design pair."""
    ma = ma or {}
    if ma.get("class") != "service_by_design" or ma.get("verified") is not False:
        return None
    causes = []
    if ma.get("direct_verified") is not True:
        causes.append("direct")
    others = ma.get("other_repeats")
    if any(isinstance(p, dict) and p.get("class") not in ("agree", "service_by_design", None)
           for p in ([others] if isinstance(others, dict) else (others or []))):
        causes.append("other_repeat")
    return "+".join(causes) or "direct"


def class_lines(st):
    """The executions an MCP arm's class (mcp_class) is read from: on an arm standing by design or not exercised
    the executions the service decided in the host's place (cells.by_design_line: the case's repeats), so that
    the class word and the status come from the same execution where the case covers one repeat only (a 200 of
    the other repeat, outside the case, is not the arm's class; VERIFY-FINAL R1); every execution otherwise."""
    return [l for l in st.lines if cells.by_design_line(st, l)] or st.lines


def _not_exercised_text(st, run, sid, short_reasons, unit_note=False):
    """The one phrase of an MCP arm the service refused or could not express (the owner's review, point 4): the
    class of the outcome with the recorded status and code (both read from the case's executions, class_lines),
    the short local reason unless the unit's note under the table carries it (`unit_note`), the other repeat
    outside the case with its recorded outcome, an unverified pair's reason, and the pointer to appendix B."""
    lines = class_lines(st)
    cls, _ = mcp_class(lines, run.scenario(sid), st.agreement, run.cases)
    case = st.agreement.get("case") or next((l.get("service_by_design_case") for l in st.lines
                                             if l.get("service_by_design_case")), None)
    glyph = "\\gni" if st.glyph_key == "not_implemented_mcp" else "\\gbd"
    head = "not exercised via MCP" if st.glyph_key == "not_implemented_mcp" else "service-side by design"
    parts = []
    if not unit_note:
        parts.append(short_reason(run, sid, case, short_reasons))
    if cls != "not_expressible":
        first = next((l for l in lines if l.get("http") is not None), None)
        code = (first.get("error") or {}).get("code") if first and isinstance(first.get("error"), dict) else None
        status = ("http %s%s" % (esc(first.get("http")), (" " + tt(code)) if code else "")) if first else None
        parts.append("%s%s" % (esc(MCP_CLASS_WORD.get(cls, cls)), (" (%s)" % status) if status else ""))
    others = other_repeats(st)
    if others:
        reps = ["r%s" % esc(r) for r, _ in others]
        parts.append("%s %s outside the case: %s (the row%s below)" % (
            "repeat" if len(others) == 1 else "repeats", esc(", ".join(reps)),
            esc("; ".join("r%s recorded %s" % (r, cells.WORD.get(l.get("outcome"), str(l.get("outcome"))))
                          for r, l in others)), "" if len(others) == 1 else "s"))
    if st.glyph_key == "service_by_design" and st.agreement.get("verified") is False:
        parts.append("UNVERIFIED pair: %s" % esc(unverified_reason(st.agreement, st), 240))
    tail = ("the note under the table and " if unit_note else "") + (("case %s in B" % tt(case)) if case else "B")
    parts.append(tail)
    return "%s\\ %s: %s" % (glyph, head, "; ".join(parts))


def short_reason(run, sid, case, short_reasons=None):
    """A short local reason for a not-exercised MCP arm: the prose's text for the case (doc/prose/cases.json
    `short`), else the manifest's not-implemented reason, else the case's missing argument or service outcome."""
    short_reasons = short_reasons or {}
    if case and short_reasons.get(case):
        return inline(short_reasons[case])
    ni = (run.scenario(sid).get("mcp") or {}).get("not_implemented")
    if ni:
        # a reason the manifest words with its own article ('the raw /pattern/capabilities document') is not
        # "the request's the ..." (VERIFY-QCND L3)
        if str(ni).lower().startswith(("the ", "a ", "an ")):
            return "%s cannot be expressed through the tools" % esc(ni, 80)
        return "the request's %s cannot be expressed through the tools" % esc(ni, 80)
    c = (run.cases or {}).get(case) or {}
    sa = c.get("service_answer") or {}
    if sa.get("missing_argument"):
        return "no MCP argument for %s" % tt(sa["missing_argument"])
    if sa.get("outcome"):
        return "the service's answer is %s" % esc(str(sa["outcome"]).replace("_", " "))
    return "the service decided the call in the host's place"


OUTCOME_WORD = {"pass": "pass", "finding": "\\textbf{FINDING}", "inconclusive": "inconclusive", "fail": "\\textbf{FAIL}",
                "error": "\\textbf{ERROR}", "not_implemented_mcp": "not exercised", "service_by_design": "by design"}


def _outcome_why(line, is_mcp):
    """The qualifier after the outcome word of an execution row, so that the row carries the adjudication of an
    ordinary case: an allowed 503, a late 503 (a finding), a host state, an inconclusive by the checks."""
    outcome = line.get("outcome")
    http = line.get("http")
    code, _ = _err(line)
    if outcome == "inconclusive":
        if _overridden(line, is_mcp):
            return "host state"
        if is_mcp and scenario_level_mark(line):
            return "scenario-level mark"
        if http == 503 and code == "deadline":
            return "allowed 503"
        if http is not None and http != 200:
            return "allowed status"
        return "by the checks"
    if outcome == "finding" and http == 503 and code == "deadline":
        return "late 503"
    return None


def _execution_rows(st, run, sid, arm, is_property):
    """One row per execution of an arm: arm / repeat, the final outcome with its qualifier, status and code, the
    wall with its host state, the requested budget, whether the property was tested (a dash for an execution the
    service decided in the host's place)."""
    rows = []
    several = len(st.lines) > 1
    for l in st.lines:
        outcome = l.get("outcome")
        if cells.by_design_line(st, l):
            # the arm's standing, as the cell shows it: a refusal the table files under not implemented reads
            # not exercised (○) in the row too, whatever outcome the runner recorded for the execution; the other
            # repeat of a by-design arm, outside the case, keeps its recorded outcome (VERIFY-PLXM H2)
            outcome = st.glyph_key
        glyph = cells.GLYPH.get(outcome, "\\gfail")
        word = OUTCOME_WORD.get(outcome, esc(cells.WORD.get(outcome, str(outcome))))
        why = _outcome_why(l, arm == "mcp")
        if why:
            word += " (%s)" % why
        hc = cells.host_cache_of(l)
        wall = cells.fmt_t(l.get("wall_ms"))
        if arm != "mcp" or hc in ("cold", "warm", cells.WARM_AFTER_INTERRUPTED):
            wall += " (%s)" % esc(HC_SHORT.get(hc, hc))
        aside = run.marked_aside(sid, l) if (arm == "mcp" and hasattr(run, "marked_aside")) else None
        tested = "–" if outcome in ("not_implemented_mcp", "service_by_design") else _property_cell(
            l, is_property, aside)
        rows.append(["%s%s" % (esc(arm), (" r%s" % esc(l.get("repeat"))) if several or arm == "mcp" else ""),
                     "%s\\ %s" % (glyph, word), _status_cell(l), wall, esc(_budget_text(run, sid)), tested])
    return rows


def _finding_items(st, arm):
    """[(arm, repeat, check, evidence, explanation)] of the finding checks of an arm's executions: the detail
    split into the evidence clause and the mechanism after it."""
    out = []
    for l in st.lines:
        if l.get("outcome") != "finding":
            continue
        for c in _finding_checks(l):
            evidence, why = _split_clause(c.get("detail"))
            out.append((arm, l.get("repeat"), c.get("check"), evidence, why))
    return out


def finding_text(items, rule):
    """The FINDING line: the evidence per execution; the mechanism once when every execution shares it (else
    with each); the rule once (VERIFY-R2 M2: never the mechanism twice)."""
    whys = []
    for _, _, _, _, why in items:
        if why and why not in whys:
            whys.append(why)
    shared = len(whys) == 1
    parts = []
    for arm, rep, check, evidence, why in items:
        text = "%s r%s, %s: %s" % (esc(arm), esc(rep), esc(check), _balanced(esc(evidence, 300)))
        if why and not shared:
            text += " — %s" % esc(why, 300)
        parts.append(text)
    text = "; ".join(parts)
    if shared:
        text += " — %s" % esc(whys[0], 400)
    if rule:
        text += " · rule: " + esc(rule, 400)
    return text


def _check_groups(line, is_mcp, marks=None, explained=None):
    """The checks of an execution that did not pass, grouped by (outcome, detail) so that a reason shared by
    several checks is said once (never 'as above', never a list of checks all untestable for one reason). An MCP
    item read not_exposed that the manifest marks at property level (|marks|: suite_common.mcp_check_marks, the
    check name's mark; set aside, or covering every property check: a scenario-level mark) is grouped under its
    mark's words with the mark's reason or case, never under not_exposed; an item the by-design table explains
    (|explained|: (indices, case), as compare derived them from the table: mcp_agreement explained_items) under
    'explained' with the case."""
    groups = []
    items, case = explained or ([], None)
    for i, c in enumerate(line.get("checks") or []):
        if not isinstance(c, dict) or c.get("outcome") == "pass":
            continue
        name = c.get("check")
        mark = (marks or {}).get(name) if c.get("outcome") == "not_exposed" else None
        if is_mcp and mark:
            key = ("marked:" + str(mark.get("mark")), str(mark.get("reason") or (
                ("case %s" % mark.get("case")) if mark.get("case") else "")))
        elif is_mcp and i in items:
            key = ("explained", ("case %s" % case) if case else "")
        else:
            key = (c.get("outcome"), c.get("detail") or "")
        for g in groups:
            if g[0] == key:
                g[1].append(c.get("check"))
                break
        else:
            groups.append((key, [c.get("check")]))
    return groups


def _group_text(key, names):
    """One group of a 'checks' line: the check names, the glyph and the reason; a marked or explained group
    (MARK_GROUP) carries the mark's words before the reason; a check read not_exposed (the MCP answer carries
    nothing the check reads: neither passed nor failed) the inconclusive glyph and the words."""
    kind, detail = key
    if kind.startswith("marked:") or kind == "explained":
        glyph, words = MARK_GROUP.get(kind.split(":", 1)[-1] if kind.startswith("marked:") else kind,
                                      ("\\gnone", kind))
        return "%s %s{} — %s%s" % (", ".join(esc(n) for n in names), glyph, words,
                                   (": " + esc(detail, 240)) if detail else "")
    if kind == "not_exposed":
        return "%s \\ginc{} — not exposed via MCP%s" % (", ".join(esc(n) for n in names),
                                                       (": " + esc(detail, 240)) if detail else "")
    # a detail cut at the column's length keeps its parentheses balanced (VERIFY-T L3)
    return "%s %s{}%s" % (", ".join(esc(n) for n in names), cells.GLYPH.get(kind, "\\gfail"),
                          (" — " + esc(_balanced(latex.cut_text(str(detail), 240)))) if detail else "")


def _texts_line(states):
    """Every distinct error text of the block, once, keyed by its status and code."""
    seen = []
    for arm in cells.ARMS:
        st = states[arm]
        for l in st.lines:
            code, text = _err(l)
            if not text or l.get("outcome") in ("not_implemented_mcp", "service_by_design"):
                continue
            key = "%s%s" % (l.get("http"), (" " + code) if code else "")
            if any(k == key and t == text for k, t in seen):
                continue
            seen.append((key, text))
    return "; ".join("%s: “%s”" % (esc(k), esc(t, 300)) for k, t in seen)


_RAW_RE = re.compile(r"^(arms/)([^/]+)(/raw/)(.+?)\.([^.]+)\.(\d+)\.json\.gz$")


def _artifacts_line(run, states):
    """The supporting artifacts, keyed once (READ-114): where every raw answer follows the one pattern
    arms/<path>/raw/<id>.<path>.<repeat>.json.gz the pattern is printed once with the (access path, repeat) keys it
    holds, else each path in full; then the executions file's physical lines."""
    raws = []
    extras = []
    numbers = []
    for arm in cells.ARMS:
        for l in states[arm].lines:
            if l.get("raw"):
                raws.append((arm, l.get("repeat"), str(l["raw"])))
            for rr in l.get("warm_up_reruns") or []:
                if isinstance(rr, dict) and rr.get("raw"):
                    extras.append("%s (the warm-up's re-run)" % tt(rr["raw"]))
            n = run.line_no(l)
            if n:
                numbers.append(str(n))
    parts = []
    patterned = [(_RAW_RE.match(raw), arm, rep) for arm, rep, raw in raws]
    sid = states[cells.ARMS[0]].lines[0].get("scenario") if states[cells.ARMS[0]].lines else None
    sid = sid or next((l.get("scenario") for a in cells.ARMS for l in states[a].lines), None)
    if raws and all(m and m.group(2) == arm and m.group(5) == arm and m.group(6) == str(rep) and m.group(4) == sid
                    for m, arm, rep in patterned):
        keys = ", ".join("%s r%s" % (cells.ARM_SHORT.get(arm, arm), esc(rep)) for arm, rep, _ in raws)
        parts.append("the raw answers %s for %s" % (
            tt("arms/<access path>/raw/%s.<access path>.<repeat>.json.gz" % sid), esc(keys)))
    else:
        parts.extend(tt(raw) for _, _, raw in raws)
    parts.extend(extras)
    text = ", ".join(parts) if parts else "no raw file recorded"
    if numbers:
        text += "; %s line%s %s" % (tt(run.executions_file), "" if len(numbers) == 1 else "s", esc(", ".join(numbers)))
    return text


DEPTH_WORD = {"full": "fully checked", "partial": "partially checked", "generic_only": "generic checks only",
              "refusal": "a refusal on both arms", "partial_non_pass": "partially checked, not a pass on both arms",
              "on_failure": "checked on a failure"}


def _repeats_text(reps):
    return facts_join(["r%s" % r for r in reps])


def facts_join(items):
    items = list(items)
    if len(items) <= 1:
        return "".join(items)
    return ", ".join(items[:-1]) + " and " + items[-1]


def depth_text(depth):
    """An agreement's depth in words over the MCP repeats compare compared, the excluded repeats named with their
    reason (the owner's second review, point 2; reader.depth_over_compared_repeats): 'fully checked on r1; r0 was
    an inconclusive host-state outcome', 'fully checked on r0 and r1', 'partially checked on r0 (not compared:
    resolve_runs)'. Plain text."""
    cls = str(depth.get("class") or "")
    text = DEPTH_WORD.get(cls, cls.replace("_", " "))
    reps = depth.get("repeats_compared") or []
    if reps:
        text += " on " + _repeats_text(reps)
    if depth.get("not_compared"):
        text += " (not compared: %s)" % latex.cut_text(", ".join(str(c) for c in depth["not_compared"]), 120)
    for e in depth.get("repeats_excluded") or []:
        if isinstance(e, dict):
            text += "; r%s was %s" % (e.get("repeat"), e.get("reason") or "left out of the comparison")
    return text


def _classes_text(run, sid, states):
    """The agreement classes in one line: the direct class with its detail where the arms differ, the MCP class
    with its depth; the by-design verification note only when the pair is unverified (else appendix B)."""
    entry = run.comp_entry(sid)
    da = entry.get("direct_agreement") or {}
    ma = entry.get("mcp_agreement") or {}
    parts = []
    cls = da.get("class")
    if cls:
        words = {"identical": "identical answers", "same_checks": "different answers, the same check outcomes",
                 "timing": "the answers differ by timing alone", "not_comparable":
                 "not comparable (no check both access paths evaluated)", "disagree": "DISAGREE"}.get(cls, cls)
        text = "direct: %s" % esc(words)
        if da.get("differing_checks"):
            text += "; differing checks: %s" % esc(", ".join(str(c) for c in da["differing_checks"]), 160)
        for t in da.get("timing_repeats") or []:
            if isinstance(t, dict) and t.get("detail"):
                text += "; r%s: %s" % (esc(t.get("repeat")), esc(t["detail"], 160))
        for d in da.get("execution_differs") or []:
            if isinstance(d, dict):
                items = ["%s local %s, remote %s" % (c.get("check"), "/".join(str(x) for x in (c.get("local") or [])),
                                                     "/".join(str(x) for x in (c.get("remote") or [])))
                         for c in (d.get("checks") or []) if isinstance(c, dict)]
                text += "; the checks differ (r%s): %s" % (esc(d.get("repeat")), esc("; ".join(items), 160))
        parts.append(text)
    mcls = ma.get("class")
    if mcls:
        text = "mcp: %s" % esc(mcls.replace("_", " ") if mcls != "not_implemented_mcp" else "not implemented via MCP")
        if ma.get("case"):
            text += ", case %s" % tt(ma["case"])
        if ma.get("checks_compared") is not None and mcls in ("agree", "disagree"):
            text += "; %s check(s) compared" % esc(ma["checks_compared"])
        depth = (run.depth or {}).get(sid)
        if depth:
            text += "; " + esc(depth_text(depth))
        if ma.get("differing_checks"):
            text += "; differing checks: %s" % esc(", ".join(str(c) for c in ma["differing_checks"]), 160)
        if mcls == "service_by_design" and ma.get("verified") is False:
            text += "; UNVERIFIED: %s" % esc(unverified_reason(ma, states.get("mcp")), 240)
        parts.append(text)
    missing = [arm for arm, st in states.items() if st.missing_planned]
    if missing:
        parts.append("planned execution(s) not recorded: %s" % esc(", ".join(missing)))
    return parts


RELATION_NOT_EVALUATED = {"not_exposed": "not evaluated: not exposed via MCP", "not_evaluated": "not evaluated"}


def _relations_text(run, sid, units):
    out = []
    comp_rel = run.comp_relations()
    for rel in run.relations:
        if sid not in (rel.get("a"), rel.get("b")):
            continue
        crel = comp_rel.get(rel.get("id")) or {}
        arms = crel.get("arms") or {}
        if all((arms.get(a) or {}).get("outcome") in (None, "pass") for a in cells.ARMS) and arms:
            continue
        per_arm = []
        details = []
        for arm in cells.ARMS:
            a = arms.get(arm) or {}
            if not a:
                continue
            o = a.get("outcome")
            if o in RELATION_NOT_EVALUATED:
                # the relation could not be evaluated on this arm (a raw answer or a field is missing via MCP):
                # neither passed nor inconclusive, the dash the sum line prints (VERIFY-QCND L5)
                per_arm.append("%s \\gnone{} (%s)" % (esc(cells.ARM_SHORT.get(arm, arm)), RELATION_NOT_EVALUATED[o]))
            else:
                per_arm.append("%s %s{}" % (esc(cells.ARM_SHORT.get(arm, arm)), cells.GLYPH.get(o, "\\ginc")))
            if a.get("detail") and o != "pass":
                details.append((arm, a["detail"]))
        other = rel.get("b") if rel.get("a") == sid else rel.get("a")
        text = "%s (%s, with %s): %s" % (latex.breakable_id(rel.get("id")), esc(rel.get("title"), 100),
                                         latex.sid(other, units.unit_of(other)), " ".join(per_arm) or "not evaluated")
        if details:
            distinct = []
            for arm, d in details:
                if d not in [x for _, x in distinct]:
                    distinct.append((arm, d))
            if len(distinct) == 1:
                text += " — %s" % esc(distinct[0][1], 120)
            else:
                text += " — " + "; ".join("%s: %s" % (esc(cells.ARM_SHORT.get(a, a)), esc(d, 100)) for a, d in distinct)
        out.append(text)
    return out


def estimate_lines(body_lines, chars_per_line=150):
    """A rough height of a block in text lines, for the short / long decision of the box (pagination)."""
    total = 0
    for text in body_lines:
        plain = re.sub(r"\\[A-Za-z]+\{[^}]*\}|\\[A-Za-z]+|[{}]", "", text)
        total += max(1, int(len(plain) / float(chars_per_line)) + 1)
    return total


SHORT_BLOCK_LINES = 28


def interpretation(run, sid, states, property_text=None):
    """The first lines of a detail block in the reader's terms (READ-114): what was requested (the unit's
    property line, else the scenario's title), what each access path did (facts.exec_phrase, one phrase per
    path) and what remains unverified (the property checks no execution of the test passed, or the statement that
    every one was shown somewhere); plain LaTeX."""
    import facts  # the facts module imports this one: a lazy import
    sc_ = run.scenario(sid)
    requested = inline(property_text) if property_text else esc(sc_.get("title"), 200)
    observed = []
    for arm in cells.ARMS:
        st = states[arm]
        if not st.present:
            if st.missing_planned:
                observed.append("%s: a planned execution not recorded" % esc(arm))
            continue
        if st.glyph_key in ("not_implemented_mcp", "service_by_design"):
            word = "not exercised via MCP" if st.glyph_key == "not_implemented_mcp" else "service-side by design"
            observed.append("%s: %s%s" % (esc(arm), word, " (pair unverified)" if "u" in st.glyph_marks else ""))
            continue
        phrases = []
        for l in st.lines:
            ph = facts.exec_phrase(run, sid, l, arm)
            if ph not in phrases:
                phrases.append(ph)
        text = phrases[0] if len(phrases) == 1 else "; ".join(
            "r%s %s" % (esc(l.get("repeat")), facts.exec_phrase(run, sid, l, arm)) for l in st.lines)
        if st.ne:
            text += " (DISAGREES with the reference access path)"
        observed.append("%s: %s" % (esc(arm), latex.join_units(esc(text))))
    props = []
    for e in sc_.get("expect") or []:
        name = e.get("check") if isinstance(e, dict) else None
        if name and run.is_property(name) and name not in props:
            props.append(name)
    unshown = [name for name in props if not any(
        c.get("check") == name and c.get("outcome") == "pass" for arm in cells.ARMS for l in states[arm].lines
        for c in (l.get("checks") or []) if isinstance(c, dict))]
    if not props:
        unverified = "none: the status, code and text are the claim, judged by the rows below"
    elif unshown:
        shown_on = []
        for arm in cells.ARMS:
            if any(c.get("check") in props and c.get("outcome") == "pass" for l in states[arm].lines
                   for c in (l.get("checks") or []) if isinstance(c, dict)):
                shown_on.append(arm)
        unverified = "%s on every access path and repeat%s" % (
            ", ".join(tt(n) for n in unshown),
            ("; the other property check%s passed on %s" % ("" if len(props) - len(unshown) == 1 else "s",
                                                             esc(", ".join(shown_on)))) if shown_on else "")
    else:
        unverified = "none of the property checks: each passed on some access path or repeat (the block explains the exception)"
    return (latex.bline("requested", requested), latex.bline("observed", "; ".join(observed) or "no execution"),
            latex.bline("unverified", unverified))


def detail_block(run, sid, states, number, units, short_reasons=None, unit_note=False, property_text=None):
    """The LaTeX of one test's detail block (the owner's review, point 1 and the first of the four further
    changes; the readability review, READ-114): the interpretation first (requested, observed, unverified), then
    the technical part: the expected line; the adjudication where the rows do not carry it (a failure, an error, a
    host-state override, a disagreement, a planned execution not recorded, the one phrase of an MCP arm the
    service decided); one row per arm and repeat (final outcome with its qualifier, HTTP and code, wall, requested
    budget, property tested?); a finding's evidence per execution with the mechanism and the rule once; the
    host-state rule and the intermediate evaluation of an overridden execution; the non-passing checks grouped by
    reason; the texts once; the classes; the relations; the supporting artifacts, keyed once. `unit_note`: the
    unit's table carries the MCP reason, so the block only points at it; `property_text`: the unit's property line
    (prose), the requested behaviour in the reader's words."""
    sc_ = run.scenario(sid)
    title = esc(sc_.get("title"), 200)
    meta = [esc(sc_.get("route") or "")]
    if sc_.get("timing_class"):
        meta.append(esc(sc_["timing_class"]))
    if sc_.get("tags"):
        meta.append(esc(", ".join(str(t) for t in sc_["tags"]), 120))
    rep = sc_.get("repeat")
    if rep:
        meta.append("%s repeat%s" % (rep, "" if rep == 1 else "s"))
    body = ["{\\color{black!60}\\small " + " · ".join(m for m in meta if m) + "}"]
    body.extend(interpretation(run, sid, states, property_text))
    body.append("{\\color{black!60}\\small\\itshape the technical record (the wire-level expectation, the rows, the "
                "checks and the artifacts):}")
    body.append(latex.bline("expected", expected_text(sc_, run.is_property)))
    adj = []
    for arm in cells.ARMS:
        st = states[arm]
        if st.present and st.glyph_key in ("not_implemented_mcp", "service_by_design"):
            text = _not_exercised_text(st, run, sid, short_reasons, unit_note)
        else:
            text = _arm_adjudication(st, run, sid, arm)
        if text:
            adj.append("\\textbf{%s}: %s" % (esc(arm), text))
    if adj:
        body.append(latex.bline("adjudication", adj[0]))
        for text in adj[1:]:
            body.append("\\bcont{%s}" % text)
    rows = []
    for arm in cells.ARMS:
        st = states[arm]
        if st.present:
            rows.extend(_execution_rows(st, run, sid, arm, run.is_property))
    if rows:
        body.append(latex.block_table(["path / repeat", "final outcome", "http · code", "wall", "requested budget",
                                       "property tested?"], rows))
    items = []
    for arm in cells.ARMS:
        items.extend(_finding_items(states[arm], arm))
    if items:
        rule = None
        for e in sc_.get("expect") or []:
            if isinstance(e, dict) and (e.get("args") or {}).get("rule"):
                rule = e["args"]["rule"]
                break
        body.append(latex.finding_line(finding_text(items, rule)))
    # the host-state rule and the intermediate evaluation of an overridden MCP execution; the non-passing checks
    # of the other executions grouped by reason
    for arm in cells.ARMS:
        st = states[arm]
        if not st.present:
            continue
        for l in st.lines:
            if cells.by_design_line(st, l):
                continue
            label = "%s r%s" % (esc(cells.ARM_SHORT.get(arm, arm)), esc(l.get("repeat")))
            for rr in l.get("warm_up_reruns") or []:
                if isinstance(rr, dict):
                    body.append(latex.bsub("warm-up re-run (%s)" % label, (
                        "%s r%s ran again for a resume %s before this execution (wall %s, http %s, %s): timing "
                        "evidence beside it, never an outcome (the warm-up's judged execution is its first)" % (
                            latex.sid(str(rr.get("scenario"))), esc(rr.get("repeat")),
                            cells.fmt_t(float(rr["before_s"]) * 1000.0) if isinstance(rr.get("before_s"), (int, float))
                            else "–", cells.fmt_t(rr.get("wall_ms")), esc(rr.get("http")),
                            esc(rr.get("outcome") or "–")))))
            if l.get("outcome") == "pass":
                continue
            if _overridden(l, arm == "mcp"):
                if l.get("verdict_detail"):
                    body.append(latex.bsub("host state (%s)" % label, esc(l["verdict_detail"])))
                groups = _check_groups(l, True, run.check_marks(sid) if hasattr(run, "check_marks") else None)
                if groups:
                    body.append(latex.bsub("intermediate (%s)" % label, "the evaluation before the host-state rule, kept as recorded: " + "; ".join(
                        _group_text((k[0], _first_clause(k[1], 160)), names) for k, names in groups)))
                continue
            marks = run.check_marks(sid) if (arm == "mcp" and hasattr(run, "check_marks")) else {}
            explained = run.explained_items(sid, l.get("repeat")) if (arm == "mcp" and hasattr(
                run, "explained_items")) else None
            groups = [g for g in _check_groups(l, arm == "mcp", marks, explained) if g[0][0] != "finding"]
            if groups and (l.get("outcome") == "fail" or (l.get("outcome") == "inconclusive" and l.get("http") == 200)
                           or l.get("outcome") == "service_by_design"):
                body.append(latex.bsub("checks (%s)" % label, "; ".join(_group_text(k, names) for k, names in groups)))
            ph = _phases_text(l) if arm == "mcp" else None
            if ph:
                body.append(latex.bsub("mcp phases (%s)" % label, esc(ph)))
    texts = _texts_line(states)
    if texts:
        body.append(latex.bline("texts", texts))
    classes = _classes_text(run, sid, states)
    if classes:
        body.append(latex.bline("classes", " · ".join(classes)))
    rels = _relations_text(run, sid, units)
    if rels:
        body.append(latex.bline("relations", "; ".join(rels)))
    body.append(latex.bline("artifacts", _artifacts_line(run, states)))
    long = estimate_lines(body) > SHORT_BLOCK_LINES
    return latex.detail_box(sid, number, title, "\n".join(body), long)


def _phases_text(line):
    ph = line.get("mcp_phases") if isinstance(line.get("mcp_phases"), dict) else None
    if not ph or ph.get("backend_ms") is None:
        return None
    keys = [("submit_ms", "submit"), ("queue_ms", "queue"), ("backend_ms", "backend"),
            ("backend_pattern_ms", "pattern"), ("terminal_lag_ms", "terminal lag"), ("paging_ms", "paging")]
    parts = ["%s %s" % (label, cells.fmt_t(ph[k])) for k, label in keys if ph.get(k) is not None]
    return ", ".join(parts)
