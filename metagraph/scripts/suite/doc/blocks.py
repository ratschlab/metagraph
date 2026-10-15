"""Expectations in a few words and the detail blocks of the tests that are not clean (DESIGN-DOC §4.4, §4.6).

`expect_text` turns one expectation (check kind + args) into a short phrase: one formatter per check kind, an
unknown kind printing `name(args)`. `expected_text` lines up a scenario's expectations: the status first, the
property checks, the generic checks last. `detail_block` writes the block of one test from its arm states and the
execution lines.
"""

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import cells  # noqa: E402
import latex  # noqa: E402
from latex import esc, tt  # noqa: E402

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
    return "clamped %s %s {\\glyphfont →} %s" % (tt(args.get("field", "?")), esc(args.get("requested")),
                                                   esc(args.get("effective")))


def _budget(args):
    text = "budget respected"
    if args.get("cap_ms") is not None:
        text += " (cap %s)" % cells.fmt_t(args["cap_ms"])
    return text


def _generic(name, args):
    if not args:
        return esc(name)
    inner = ", ".join("%s=%s" % (esc(k), esc(_runs(v), 40)) for k, v in sorted(args.items()))
    return "%s(%s)" % (esc(name), inner)


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


def refusal_expected(scenario, previous=None):
    """The expected column of a refusals row: the status expectation; 'the same text' when it repeats the
    previous row's."""
    for e in scenario.get("expect") or []:
        if isinstance(e, dict) and e.get("check") == "status":
            args = e.get("args") or {}
            text = args.get("text") or args.get("text_prefix") or args.get("text_contains")
            if previous is not None and text and text == previous:
                head = "http %s" % args.get("http") if args.get("http") is not None else ""
                return (head + " the same text").strip(), text
            return _status(args), text
    return "–", None


# -- detail blocks ---------------------------------------------------------------------------------------------

def _error_text(line):
    err = line.get("error") if isinstance(line.get("error"), dict) else None
    if not err:
        return ""
    code = err.get("code")
    text = err.get("text")
    out = ""
    if code:
        out += " %s" % tt(code)
    if text:
        out += " “%s”" % esc(text, 160)
    return out


def _phases_text(line):
    ph = line.get("mcp_phases") if isinstance(line.get("mcp_phases"), dict) else None
    if not ph or ph.get("backend_ms") is None:
        return None
    keys = [("submit_ms", "submit"), ("queue_ms", "queue"), ("backend_ms", "backend"),
            ("backend_pattern_ms", "pattern"), ("terminal_lag_ms", "terminal lag"), ("paging_ms", "paging")]
    parts = ["%s %s" % (label, cells.fmt_t(ph[k])) for k, label in keys if ph.get(k) is not None]
    return "phases: " + ", ".join(parts)


def _repeat_lines(line, run_, is_mcp):
    """The lines of one execution that is not clean: the header line and one line per check not passed."""
    out = []
    outcome = line.get("outcome")
    glyph = cells.GLYPH.get(outcome, "\\gfail")
    head = "%s %s" % (glyph, esc(cells.WORD.get(outcome, outcome)))
    if line.get("http") is not None:
        head += " · http %s" % esc(line.get("http"))
    head += _error_text(line)
    hc = cells.host_cache_of(line)
    hc_text = {"cold": "host-cold", "warm": "host-warm", cells.WARM_AFTER_INTERRUPTED:
               "host-warm, after an interrupted execution", "aged": "host cache aged", "not_reached":
               "host not reached", "unknown": "host cache unknown"}.get(hc, hc)
    head += " · wall %s (%s)" % (cells.fmt_t(line.get("wall_ms")), esc(hc_text))
    srv = (line.get("server") or {}).get("elapsed_ms") if isinstance(line.get("server"), dict) else None
    head += " · server %s" % (cells.fmt_t(srv) if srv is not None else "–")
    out.append(head)
    checks = [c for c in (line.get("checks") or []) if isinstance(c, dict)]
    passed = [c for c in checks if c.get("outcome") == "pass"]
    others = [c for c in checks if c.get("outcome") != "pass"]
    sub = []
    if passed:
        sub.append(("passed", ", ".join("%s \\gpass" % esc(c.get("check")) for c in passed)))
    prev = None
    for c in others:
        detail = c.get("detail") or ""
        shown = "as above" if (detail and detail == prev) else esc(detail, 240)
        prev = detail or prev
        sub.append(("%s %s" % (esc(c.get("check")), cells.GLYPH.get(c.get("outcome"), "\\gfail")), shown))
    if is_mcp:
        if line.get("inconclusive_reason") and line.get("inconclusive_reason") != "checks":
            sub.append(("reason", esc(line["inconclusive_reason"])))
        if line.get("verdict_detail"):
            sub.append(("host state", esc(line["verdict_detail"], 240)))
        if line.get("service_by_design_case"):
            sub.append(("by-design case", tt(line["service_by_design_case"])))
        ph = _phases_text(line)
        if ph:
            sub.append(("mcp", esc(ph)))
    if line.get("time_limited"):
        sub.append(("time-limited", esc(line["time_limited"], 160)))
    return out, sub


def _clean_arm_text(st, run, sid):
    """One line for an arm that is clean."""
    if not st.present:
        return "no execution (none planned)"
    if st.glyph_key == "not_implemented_mcp":
        case = st.agreement.get("case")
        text = "\\gni\\ not implemented via MCP"
        if case:
            c = run.cases.get(case) or {}
            text += ": table case %s" % tt(case)
            if c.get("scenario_match"):
                text += " (%s)" % esc(c["scenario_match"], 160)
        elif (run.scenario(sid).get("mcp") or {}).get("not_implemented"):
            text += ": manifest, %s" % esc((run.scenario(sid).get("mcp") or {}).get("not_implemented"), 120)
        refused = [l for l in st.lines if l.get("service_by_design_case")]
        if refused:
            l = refused[0]
            text += " · the service refused each call by design, case %s%s · %s, host not reached" % (
                tt(l["service_by_design_case"]), _error_text(l), cells.fmt_times(st.walls))
        return text
    if st.glyph_key == "service_by_design":
        case = st.agreement.get("case") or next((l.get("service_by_design_case") for l in st.lines
                                                if l.get("service_by_design_case")), None)
        text = "\\gbd\\ service-side by design"
        if case:
            text += ", case %s" % tt(case)
        text += ", verified" if st.agreement.get("verified") else (", unverified" if st.agreement.get("verified")
                                                                   is False else "")
        if st.lines:
            text += " · %s%s" % (cells.fmt_times(st.walls), _error_text(st.lines[0]))
        return text
    return "%s %s: every check passed" % (st.glyph_tex(), cells.fmt_times(st.walls, st.time_marks))


def _classes_text(run, sid, states):
    entry = run.comp_entry(sid)
    da = entry.get("direct_agreement") or {}
    ma = entry.get("mcp_agreement") or {}
    parts = []
    cls = da.get("class")
    if cls:
        words = {"identical": "identical answers", "same_checks": "different answers, the same check outcomes",
                 "timing": "the answers differ by timing alone", "not_comparable":
                 "not comparable (no check both arms evaluated)", "disagree": "DISAGREE"}.get(cls, cls)
        text = "direct-local vs direct-remote: %s" % esc(words)
        if da.get("differing_checks"):
            text += "; differing checks: %s" % esc(", ".join(str(c) for c in da["differing_checks"]), 200)
        for t in da.get("timing_repeats") or []:
            if isinstance(t, dict) and t.get("detail"):
                text += "; r%s: %s" % (esc(t.get("repeat")), esc(t["detail"], 200))
        for d in da.get("execution_differs") or []:
            if isinstance(d, dict):
                items = []
                for c in d.get("checks") or []:
                    if isinstance(c, dict):
                        items.append("%s local %s, remote %s" % (
                            c.get("check"), "/".join(str(x) for x in (c.get("local") or [])),
                            "/".join(str(x) for x in (c.get("remote") or []))))
                text += "; the execution checks differ (r%s): %s" % (esc(d.get("repeat")), esc("; ".join(items), 200))
        parts.append(text)
    mcls = ma.get("class")
    if mcls:
        text = "mcp: %s" % esc(mcls.replace("_", " ") if mcls != "not_implemented_mcp" else "not implemented via MCP")
        if ma.get("case"):
            text += ", case %s" % tt(ma["case"])
        if ma.get("by"):
            text += " by the %s" % esc(ma["by"])
        if ma.get("side"):
            text += " (%s side)" % esc(ma["side"])
        if ma.get("checks_compared") is not None and mcls in ("agree", "disagree"):
            text += "; %s check(s) compared" % esc(ma["checks_compared"])
        depth = (run.depth or {}).get(sid)
        if depth:
            text += "; depth %s" % esc(depth.get("class", "").replace("_", " "))
            if depth.get("not_compared"):
                text += " (not compared: %s)" % esc(", ".join(depth["not_compared"]), 120)
        if ma.get("differing_checks"):
            text += "; differing checks: %s" % esc(", ".join(str(c) for c in ma["differing_checks"]), 200)
        dis = ma.get("disagreement")
        if isinstance(dis, dict) and dis.get("direct_detail"):
            text += "; %s" % esc(dis["direct_detail"], 240)
        parts.append(text)
    missing = [arm for arm, st in states.items() if st.missing_planned]
    if missing:
        parts.append("planned execution(s) not recorded: %s" % esc(", ".join(missing)))
    return parts


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
        for arm in cells.ARMS:
            a = arms.get(arm) or {}
            if not a:
                continue
            per_arm.append("%s %s%s" % (esc(arm), cells.GLYPH.get(a.get("outcome"), "\\ginc"),
                                        (" " + esc(a.get("detail"), 120)) if a.get("detail") else ""))
        other = rel.get("b") if rel.get("a") == sid else rel.get("a")
        out.append("%s (%s, %s): %s" % (latex.breakable_id(rel.get("id")), esc(rel.get("title"), 120),
                                        latex.sid(other, units.unit_of(other)), " · ".join(per_arm) or "not evaluated"))
    return out


def detail_block(run, sid, states, number, units):
    """The LaTeX of one test's detail block (DESIGN-DOC §4.6)."""
    sc_ = run.scenario(sid)
    head = "\\fn{%d}\\ %s\\quad %s" % (number, latex.breakable_id(sid), esc(sc_.get("title"), 200))
    meta = [esc(sc_.get("route") or "")]
    if sc_.get("timing_class"):
        meta.append(esc(sc_["timing_class"]))
    if sc_.get("tags"):
        meta.append(esc(", ".join(str(t) for t in sc_["tags"]), 120))
    rep = sc_.get("repeat")
    if rep:
        meta.append("%s repeat%s" % (rep, "" if rep == 1 else "s"))
    body = [head + "\\par{\\color{black!60}\\small " + " · ".join(m for m in meta if m) + "}"]
    body.append(latex.bline("expected", expected_text(sc_, run.is_property)))
    for arm in cells.ARMS:
        st = states[arm]
        if st.clean:
            body.append(latex.bline(esc(arm), _clean_arm_text(st, run, sid)))
            continue
        if not st.present:
            body.append(latex.bline(esc(arm), "\\gnone\\ no execution recorded although the run planned one"))
            continue
        if st.glyph_key in ("not_implemented_mcp", "service_by_design"):
            body.append(latex.bline(esc(arm), _clean_arm_text(st, run, sid)))
            continue
        several = len(st.lines) > 1
        for line in st.lines:
            if line.get("outcome") == "pass" and several:
                body.append(latex.bsub("r%s" % esc(line.get("repeat")), "%s pass %s: every check passed" % (
                    "\\gpass", cells.fmt_t(line.get("wall_ms")))))
                continue
            heads, sub = _repeat_lines(line, run, arm == "mcp")
            key = esc(arm) if not several else "%s r%s" % (esc(arm), esc(line.get("repeat")))
            body.append(latex.bline(key, heads[0]))
            for k, v in sub:
                body.append(latex.bsub(k, v))
        if st.arm == "direct-remote" and st.agreement.get("class") in ("timing", "same_checks", "not_comparable",
                                                                       "disagree"):
            pass  # said under classes
    classes = _classes_text(run, sid, states)
    if classes:
        body.append(latex.bline("classes", " · ".join(classes)))
    rels = _relations_text(run, sid, units)
    if rels:
        body.append(latex.bline("relations", "; ".join(rels)))
    if any(c.get("outcome") == "finding" for st in states.values() for l in st.lines
           for c in (l.get("checks") or []) if isinstance(c, dict)):
        for e in sc_.get("expect") or []:
            if isinstance(e, dict) and (e.get("args") or {}).get("rule"):
                body.append(latex.bline("rule", esc(e["args"]["rule"], 400)))
                break
    return latex.detail_box(sid, "\n".join(body))
