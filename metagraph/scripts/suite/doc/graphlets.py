"""Chapter G of the results document (manifest v8): the graphlet tools compared five ways (DESIGN-v8 §1-§6), computed
from the run's records by compare_graphlets.load_or_compute, the same in-memory comparison the report's chapter G and
the deploy gate read (nothing is written): the five ways and what their agreement proves, the comparison per tool,
every disagreement with both answers, the pre-checks, the service's own behaviour, the correctness coverage (W5's
checked fields, those it re-derives from DESIGN, the consistency-only ones) and what the chapter does not cover. And
the release records of chapter S's stop scenarios. A run whose manifest has no graphlet scenario (v4..v7) has no
chapter G and nothing here is printed."""

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.dirname(HERE)
if SUITE not in sys.path:
    sys.path.insert(0, SUITE)

import cells  # noqa: E402
import facts  # noqa: E402
import latex  # noqa: E402
from latex import esc, tt, inline  # noqa: E402

WAY_ORDER = ("W1", "W2", "W3", "W4", "W5")
# the five ways (DESIGN-v8 §2): the document each processes and what processes it
WAYS = (("W1", "dl+lib", "direct-local's graphlet answer, turned into the service's document form by the library",
         "the library, locally, at the deployed library commit"),
        ("W2", "dr+lib", "the same request on direct-remote, the same conversion", "the same"),
        ("W3", "mcp", "the service's stored body", "the service's graphlet tools (its vendored library)"),
        ("W4", "mcp+lib", "the graphlet the MCP exported (byte-identical to the stored body)",
         "the library, locally, at the same commit"),
        ("W5", "oracle", "the full twin (direct-local) and the document's records",
         "the suite's own reader, without the library"))
# what an agreement of each pair proves (DESIGN-v8 §2): consistency, except W5's checked fields
PROVES = (("pre-checks", "the three library ways read one document: the bodies equal (W1, W2, W4), the envelopes "
                         "equal but for per-request values, the full twin's paths equal the graphlet's, the MCP's "
                         "strategy the one the direct requests sent"),
          ("W1 = W2", "transport: the two direct access paths answered the same graphlet, processed alike"),
          ("W1 = W4", "fetch path: the service's stored body is the engine's answer"),
          ("W3 = W4", "the service's wrapper: its tools answer as the library at the same commit"),
          ("W1 = W2 = W4", "the library's determinism: each way runs in its own process with its own hash seed"),
          ("each = W5", "correctness, on W5's checked fields only: the library's answer equals an independent "
                        "reading of the engine's values (a field W5 re-derives from a DESIGN formula catches a slip "
                        "of the implementation, not a wrong rule)"))
MAX_ITEMS = 40


def data(run):
    """The five-way comparison of |run| (compare_graphlets.load_or_compute, recomputed in memory), cached on the run;
    None for a manifest without graphlet scenarios; {"error": ...} when it cannot be computed (the chapter says so)."""
    if hasattr(run, "_graphlets"):
        return run._graphlets
    out = None
    try:
        import graphlet_precheck as gp
        if gp.graphlet_scenarios(run.manifest if isinstance(run.manifest, dict) else {}):
            import compare_graphlets as cg
            try:
                out = cg.load_or_compute(run.dir, run.manifest, False)
            except Exception as e:   # the chapter says so; the rest of the document stands
                out = {"error": "%s: %s" % (type(e).__name__, e)}
    except ImportError:
        out = None
    run._graphlets = out
    return out


def verdict(g):
    """(word, reason) of area G: report.graphlet_verdict's."""
    import report
    return report.graphlet_verdict(g)


def _json(value, limit=160):
    """A value compactly: a string the comparison already holds as text (a JSON value's text, "(absent)") as it is,
    anything else as compact JSON; cut to |limit| characters."""
    text = value if isinstance(value, str) else json.dumps(value, sort_keys=True, ensure_ascii=False,
                                                          separators=(",", ":"))
    return text if len(text) <= limit else text[:limit - 1] + "…"


def _n(x):
    return str(x) if x else "–"


def ways_section(ctx):
    rows = [[tt(w), esc(name), esc(doc), esc(by)] for w, name, doc, by in WAYS]
    out = latex.section("The five ways", "g:ways", level="subsection")
    out += latex.plain_table(latex.spec_p([12, 18, 72, 62]), ["way", "", "the document", "processed by"], rows,
                             "the ways of an op call (DESIGN-v8 §2)")
    out += latex.plain_table(latex.spec_p([26, 138]), ["agreement", "what it proves"],
                             [[esc(k), esc(v)] for k, v in PROVES],
                             "what an agreement proves: W1 to W4 agreeing is consistency (they could share a library "
                             "defect); W5 is the one correctness check")
    return out


def per_tool_section(g):
    scen = g.get("scenarios") or {}
    per_way = {}
    for sv in scen.values():
        for e in (sv.get("ops") or {}).values():
            t = per_way.setdefault(e.get("tool"), {})
            for w, st in (e.get("ways") or {}).items():
                if st in ("answered", "error", "selector_error"):
                    t[w] = t.get(w, 0) + 1
    pt = (g.get("summary") or {}).get("per_tool") or {}
    rows = []
    tot = {}
    for tool in sorted(pt):
        c = pt[tool]
        for k in ("calls", "agree", "disagree", "service_by_design", "not_exposed", "inconclusive",
                  "precondition_failed"):
            tot[k] = tot.get(k, 0) + (c.get(k) or 0)
        pw = per_way.get(tool) or {}
        rows.append([tt(tool), str(c.get("calls", 0)),
                     "{\\scriptsize %s}" % esc(", ".join("%s %d" % (w, pw[w]) for w in WAY_ORDER if pw.get(w)) or "–"),
                     _n(c.get("agree")), cells.attn(_n(c.get("disagree"))) if c.get("disagree") else "–",
                     _n(c.get("service_by_design")), _n(c.get("not_exposed")),
                     cells.attn(_n(c.get("inconclusive"))) if c.get("inconclusive") else "–",
                     cells.attn(_n(c.get("precondition_failed"))) if c.get("precondition_failed") else "–"])
    rows.append(["\\textbf{all tools}", "\\textbf{%d}" % tot.get("calls", 0), ""] + [
        "\\textbf{%s}" % _n(tot.get(k)) for k in ("agree", "disagree", "service_by_design", "not_exposed",
                                                   "inconclusive", "precondition_failed")])
    out = latex.section("The comparison per tool", "g:tools", level="subsection")
    out += latex.plain_table(latex.spec_p([30, 13, 44, 13, 13, 13, 13, 15, 13],
                                          aligns=["raggedright", "raggedleft", "raggedright"] + ["raggedleft"] * 6),
                             ["tool", "op calls", "records per way", "agree", "disagree", "by design",
                              "W3 not exposed", "incon\\-clusive", "pre\\-condition failed"], rows,
                             "per tool its op calls and their outcomes; an op call agrees when every pair of the ways "
                             "it runs on agrees; W3 not exposed counts the calls the service does not offer (the "
                             "owner's reasons, below), which are no gap of the comparison")
    return out


def _pair_line(p):
    """One disagreeing pair: its name and what it isolates, where the answers differ and both values."""
    a, _, b = str(p.get("pair") or "-").partition("-")
    text = "%s (%s)" % (esc(p.get("pair")), esc(p.get("what")))
    if p.get("pointer"):
        text += " at %s: %s %s, %s %s" % (tt(p["pointer"], 120), esc(a), tt(_json(p.get(a), 100)), esc(b),
                                           tt(_json(p.get(b), 100)))
    elif p.get("detail"):
        text += ": %s" % esc(p["detail"], 240)
    lines = [text]
    for w in (a, b):
        ans = (p.get("answers") or {}).get(w)
        if ans:
            lines.append("%s's answer: %s" % (esc(w), tt(_json(ans, 220))))
    return lines


def disagreements_section(g, max_items=MAX_ITEMS):
    """Every disagreeing op call with its pairs and both answers (compactly), then the inconclusive and
    precondition_failed calls (each blocks the gate); at most |max_items| of each, the rest counted."""
    scen = g.get("scenarios") or {}
    out = latex.section("Disagreements", "g:disagree", level="subsection")
    dis = [(sid, oid, e) for sid, sv in sorted(scen.items()) for oid, e in sorted((sv.get("ops") or {}).items())
           if e.get("outcome") == "disagree"]
    if not dis:
        out += "No op call disagreed.\n\n"
    else:
        # one entry per kind of difference (the tool, the disagreeing pairs and where they differ), its op calls
        # named, the first call's values and answers shown: a defect that strikes many calls reads once
        groups = []
        index = {}
        for sid, oid, e in dis:
            sub = []
            key_parts = [e.get("tool"), tuple(e.get("disagreeing") or [])]
            for p in e.get("pairs") or []:
                if p.get("outcome") == "disagree":
                    sub.extend(_pair_line(p))
                    key_parts.append((p.get("pair"), p.get("pointer")))
            for c in e.get("cuts") or []:
                if c.get("outcome") == "disagree":
                    sub.append("%s against its uncut answer: %s%s" % (esc(c.get("way")), latex.esc_wrap(
                        c.get("detail"), 200), (" at %s" % tt(c["pointer"])) if c.get("pointer") else ""))
                    key_parts.append(("cut", c.get("way"), c.get("pointer")))
            for r in e.get("by_design") or []:
                if not r.get("verified"):
                    sub.append("by-design %s, %s side: %s" % (tt(r.get("case")), esc(r.get("side")),
                                                             latex.esc_wrap(r.get("detail"), 200)))
                    key_parts.append(("bd", r.get("case"), r.get("side")))
            key = json.dumps(key_parts, default=str)
            if key not in index:
                index[key] = len(groups)
                groups.append({"tool": e.get("tool"), "disagreeing": e.get("disagreeing") or [], "calls": [],
                               "sub": sub})
            groups[index[key]]["calls"].append((sid, oid))
        out += ("%d op call%s disagreed, in %d kind%s of difference; each kind with its op calls, the pairs that "
                "differ, where they differ and both values of the first call, then its two answers (cut to a "
                "line):\n" % (len(dis), "" if len(dis) == 1 else "s", len(groups), "" if len(groups) == 1 else "s"))
        out += "{\\small\\begin{itemize}\\setlength\\itemsep{1pt}\n"
        for grp in groups[:max_items]:
            calls = grp["calls"]
            named = "; ".join("%s %s" % (latex.breakable_id(sid), tt(oid)) for sid, oid in calls[:8]) + (
                "; … %d more" % (len(calls) - 8) if len(calls) > 8 else "")
            out += "\\item %s, %d call%s (%s): %s" % (tt(grp["tool"]), len(calls), "" if len(calls) == 1 else "s",
                                                     named, esc(", ".join(grp["disagreeing"]) or "–"))
            if grp["sub"]:
                out += "\\begin{itemize}\\setlength\\itemsep{0pt}\n" + "".join(
                    "\\item %s\n" % x for x in grp["sub"]) + "\\end{itemize}\n"
            out += "\n"
        if len(groups) > max_items:
            out += "\\item … %d more kinds (the report's chapter G lists every call)\n" % (len(groups) - max_items)
        out += "\\end{itemize}}\n\n"
    other = [(sid, oid, e) for sid, sv in sorted(scen.items()) for oid, e in sorted((sv.get("ops") or {}).items())
             if e.get("outcome") in ("inconclusive", "precondition_failed")]
    if other:
        out += ("%d op call%s inconclusive or with failed preconditions (each blocks the gate):\n" % (
            len(other), " was" if len(other) == 1 else "s were"))
        out += "{\\small\\begin{itemize}\\setlength\\itemsep{0pt}\n"
        for sid, oid, e in other[:max_items]:
            out += "\\item %s %s (%s): %s: %s\n" % (latex.breakable_id(sid), tt(oid), tt(e.get("tool")),
                                                  esc((e.get("outcome") or "").replace("_", " ")),
                                                  latex.esc_wrap(e.get("reason"), 200))
        if len(other) > max_items:
            out += "\\item … %d more (the report's chapter G lists them all)\n" % (len(other) - max_items)
        out += "\\end{itemize}}\n\n"
    return out


PRECHECKS = ("body", "envelope", "twin", "strategy", "preconditions", "library_commit")


def _outcome_tex(word):
    """An outcome word for a narrow column: hyphenable, red unless pass, not exposed or none."""
    word = word or "–"
    text = esc(word.replace("_", " ")).replace("inconclusive", "incon\\-clusive").replace(
        "precondition", "pre\\-condition")
    return text if word in ("pass", "not_exposed", "–") else cells.attn(text)


def prechecks_section(g, max_items=MAX_ITEMS):
    scen = g.get("scenarios") or {}
    rows, bad = [], []
    for sid, sv in sorted(scen.items()):
        pre = sv.get("precheck") or {}
        rows.append([latex.breakable_id(sid)] + [_outcome_tex((pre.get(k) or {}).get("outcome")) for k in PRECHECKS]
                    + [_outcome_tex(pre.get("outcome")) + ((" (%s)" % esc(pre["finding"])) if pre.get("finding")
                                                          else "")])
        bad.extend((sid, k, pre.get(k)) for k in PRECHECKS
                   if (pre.get(k) or {}).get("outcome") not in ("pass", "not_exposed"))
    out = latex.section("Pre-checks", "g:prechecks", level="subsection")
    out += latex.plain_table(latex.spec_p([34, 15, 17, 13, 16, 21, 19, 27],
                                          aligns=["raggedright"] * 8),
                             ["scenario", "body", "envelope", "twin", "strategy", "pre\\-conditions", "library commit",
                              "scenario"], rows,
                             "the pre-checks per scenario, before any op call is compared (DESIGN-v8 §2): the bodies "
                             "of W1, W2 and W4 equal, their envelopes equal but for per-request values, the full "
                             "twin's paths equal the graphlet's, the MCP's strategy echo the direct requests' "
                             "strategy, the structural preconditions frozen in the manifest, the library commit of "
                             "every MCP call the deployed and pinned one; the last column the scenario's outcome")
    if bad:
        out += "{\\small\\begin{itemize}\\setlength\\itemsep{0pt}\n"
        for sid, k, c in bad[:max_items]:
            out += "\\item %s %s: %s: %s\n" % (latex.breakable_id(sid), esc(k.replace("_", " ")),
                                             esc(((c or {}).get("outcome") or "–").replace("_", " ")),
                                             latex.esc_wrap((c or {}).get("detail"), 240))
        if len(bad) > max_items:
            out += "\\item … %d more\n" % (len(bad) - max_items)
        out += "\\end{itemize}}\n\n"
    return out


def service_checks_section(g):
    req = set(g.get("required_service_checks") or [])
    rows = [[tt(c.get("check")), _outcome_tex(c.get("outcome")) if c.get("outcome") != "n/a" else "n/a",
             str(c.get("calls", 0)), "yes" if c.get("check") in req else "no", latex.esc_wrap(c.get("detail"), 260)]
            for c in g.get("service_checks") or []]
    out = latex.section("The service's own behaviour", "g:service", level="subsection")
    out += latex.plain_table(latex.spec_p([30, 20, 12, 15, 87], aligns=["raggedright", "raggedright", "raggedleft",
                                                                        "raggedright", "raggedright"]),
                             ["check", "outcome", "calls", "required", "detail"], rows,
                             "the service's own fields and refusals, checked positively on W3's calls (DESIGN-v8 §5): "
                             "they are stripped before the ways are compared; n/a: no call met the check's "
                             "precondition; required: the gate requires the check (the op plan exercises it)")
    return out


def coverage_section(g):
    """W5's correctness coverage per tool: the fields it checks against an engine value or a SPEC token, those it
    re-derives from a DESIGN formula, those it leaves to consistency (list indices written *)."""
    cov = g.get("coverage") or {}
    pt = (g.get("summary") or {}).get("per_tool") or {}
    rows = []
    for tool in sorted(cov):
        checked, rederived = [], []
        for ptr, c in sorted((cov[tool].get("checked") or {}).items()):
            (rederived if c.get("rederived") else checked).append(ptr)
        cons = sorted(cov[tool].get("consistency_only") or {})
        for kind, ptrs in (("checked", checked), ("checked, re-derived from DESIGN", rederived),
                           ("consistency only", cons)):
            if ptrs:
                rows.append([tt(tool), esc(kind), str(len(ptrs)),
                             "{\\scriptsize %s}" % ", ".join(tt(p) for p in ptrs[:24]) + (
                                 " … %d more" % (len(ptrs) - 24) if len(ptrs) > 24 else "")])
    out = latex.section("Correctness coverage", "g:coverage", level="subsection")
    out += ("W1 to W4 agreeing shows consistency; W5 is the only correctness check. Per tool, the fields W5 checks "
            "against an engine value or a SPEC token, the fields it re-derives from a DESIGN formula over engine "
            "values (they catch a slip of the library's implementation, not a wrong rule), and the fields it leaves "
            "to consistency. The clause each checked field rests on is in the report's chapter G.\n\n")
    out += latex.plain_table(latex.spec_p([28, 30, 10, 96], aligns=["raggedright", "raggedright", "raggedleft",
                                                                   "raggedright"]),
                             ["tool", "W5", "fields", "the fields (JSON pointers, a list index written *)"], rows,
                             "W5's coverage per tool")
    no_w5 = sorted(t for t in pt if t not in cov)
    if no_w5:
        out += "Tools W5 does not read in this run (consistency only): %s.\n\n" % ", ".join(tt(t) for t in no_w5)
    return out


def by_design_and_drift_section(g, max_items=MAX_ITEMS):
    scen = g.get("scenarios") or {}
    out = ""
    bd = [(sid, oid, r) for sid, sv in sorted(scen.items()) for oid, e in sorted((sv.get("ops") or {}).items())
          for r in e.get("by_design") or []]
    if bd:
        out += latex.section("Declared per-way differences", "g:bydesign", level="subsection")
        out += latex.plain_table(latex.spec_p([30, 26, 38, 14, 56]),
                                 ["scenario", "op", "case", "verified", "detail"],
                                 [[latex.breakable_id(sid), tt(oid), tt(r.get("case")),
                                   "yes" if r.get("verified") else cells.attn("no (%s)" % esc(r.get("side"))),
                                   esc(r.get("detail"), 200)] for sid, oid, r in bd[:max_items]],
                                 "the by-design pairs: a difference of W3 the service documents (the area G by-design "
                                 "table), verified on both sides")
    d = g.get("drift") or {}
    if d.get("calls"):
        out += latex.section("The drift lane", "g:drift", level="subsection")
        out += ("Report-only: the branch head's library over W4's documents, %d call%s, %d compared with W4, %d "
                "differing.\n\n" % (d.get("calls", 0), "" if d.get("calls") == 1 else "s", d.get("compared", 0),
                                    len(d.get("differ") or [])))
    return out


def not_covered(run):
    """{library_only: [(reason, [tools])], direct_only: [(scenario, reason)], v9: [(scenario, reason)]} of the
    manifest: the tools and arguments the service does not expose (the owner's reasons of graphlet.not_exposed, the
    store-management smoke apart), the slots the MCP cannot express, the slots left to manifest v9."""
    by_reason = {}
    direct = []
    for s in (run.manifest or {}).get("scenarios") or []:
        g = s.get("graphlet") if isinstance(s, dict) else None
        if not isinstance(g, dict):
            continue
        ne = g.get("not_exposed") or {}
        for key, why in ne.items():
            if key == "*":
                continue
            by_reason.setdefault(str(why), set()).add(key)
        if "W3" not in ((g.get("ways") or {}).get("default") or []) and ne.get("*"):
            direct.append((s["id"], str(ne["*"])))
    v9 = [(x.get("scenario"), x.get("reason")) for x in (run.manifest or {}).get("skipped") or []
          if isinstance(x, dict) and x.get("area") == "G"]
    return {"library_only": sorted((why, sorted(keys)) for why, keys in by_reason.items()), "direct_only": direct,
            "v9": v9}


def not_covered_section(run):
    nc = not_covered(run)
    out = latex.section("What this chapter does not cover", "g:notcovered", level="subsection")
    items = []
    for why, keys in nc["library_only"]:
        items.append("%s: %s" % (esc(why.rstrip(".")), ", ".join(tt(k) for k in keys)))
    out += ("Via MCP: the tools and arguments the service does not expose, by the owner's reason (they run on the "
            "library ways and W5 only, and W3's absence is no gap of the comparison):\n")
    out += "{\\small\\begin{itemize}\\setlength\\itemsep{0pt}\n" + "".join("\\item %s\n" % x for x in items) + \
        "\\end{itemize}}\n\n"
    if nc["direct_only"]:
        reasons = sorted(set(w for _, w in nc["direct_only"]))
        out += ("Direct only (W1, W2 and W5; W3 and W4 not exposed): %s (%s).\n\n" % (
            ", ".join(latex.breakable_id(s) for s, _ in nc["direct_only"]), esc("; ".join(reasons))))
    if nc["v9"]:
        out += "Left to manifest v9: %s.\n\n" % "; ".join(
            "%s (%s)" % (latex.breakable_id(s), esc(w)) for s, w in nc["v9"])
    return out


def chapter_g_sections(ctx):
    """Chapter G's comparison sections, after its background and before its families; '' without a comparison."""
    g = data(ctx.run)
    if g is None:
        return ""
    if g.get("error"):
        return "\\parhead{The five-way comparison}\nThe graphlet comparison could not be computed: %s.\n\n" % esc(
            g["error"], 300)
    out = ways_section(ctx)
    out += per_tool_section(g)
    out += disagreements_section(g)
    out += prechecks_section(g)
    out += service_checks_section(g)
    out += coverage_section(g)
    out += by_design_and_drift_section(g)
    out += not_covered_section(ctx.run)
    return out


# -- chapter S: the stops' release records --------------------------------------------------------------------

def stop_rows(run):
    """[(line, stop_flow)] of the stop scenarios' MCP executions (manifest v8: mcp.flow kind stop)."""
    out = []
    for sid in run.order:
        s = run.scenario(sid)
        if not (isinstance(s.get("mcp"), dict) and isinstance(s["mcp"].get("flow"), dict)):
            continue
        for l in run.arm_lines(sid, "mcp"):
            body = l.get("stop_flow")
            if body is None:
                ans = (l.get("answer") or {}).get("body") if isinstance(l.get("answer"), dict) else None
                body = (ans or {}).get("stop_flow") if isinstance(ans, dict) else None
            out.append((l, body if isinstance(body, dict) else None))
    return out


def chapter_s_sections(ctx):
    """Chapter S's records: per stop execution the stop's answer, the job's end, each host attempt's release record
    (reason, site, decision, detail) and the release reads; per base request its full-run times."""
    run = ctx.run
    rows = []
    for l, rec in stop_rows(run):
        if rec is None:
            rows.append([latex.breakable_id(l.get("scenario")), "r%s" % l.get("repeat"),
                         cells.attn("no stop record"), "–", "–", "–"])
            continue
        stop = rec.get("stop") or {}
        stop_text = "%s; running %s, cancel %s, bound %s s" % (
            "taken" if stop.get("took") else cells.attn("not taken"), esc(stop.get("running", "–")),
            "supported" if stop.get("cancel_supported") else "not supported", esc(stop.get("bound_seconds", "–")))
        atts = []
        for a in rec.get("attempts") or []:
            r = a.get("release") or {}
            atts.append("%s: %s, site %s, decision %s, %s%s" % (
                tt(str(a.get("attempt_id") or "–")[:8] + "…"), esc(r.get("reason") or "null"), esc(r.get("site") or "–"),
                esc(r.get("decision") if r.get("decision") is not None else "null"), esc(r.get("detail") or "–"),
                (", instance %s" % tt(a["server_instance"])) if a.get("server_instance") else ""))
        rows.append([latex.breakable_id(l.get("scenario")), "r%s" % l.get("repeat"), stop_text,
                     esc(rec.get("terminal_status") or "–"), "; ".join(atts) or cells.attn("no attempt"),
                     esc(rec.get("release_reads", "–"))])
    if not rows:
        return ""
    out = latex.section("The stops' records", "s:records", level="subsection")
    out += latex.plain_table(latex.spec_p([26, 8, 40, 16, 64, 12], aligns=["raggedright"] * 5 + ["raggedleft"]),
                             ["scenario", "", "the stop's answer", "job's end", "each host attempt's release record",
                              "release reads"], rows,
                             "per stop execution what the MCP runner recorded: the stop's answer, the job's terminal "
                             "status, each host attempt's release record as the service's admin job record states it "
                             "(reason, site, decision, detail), and how often the runner read it")
    base = []
    for sid in run.order:
        s = run.scenario(sid)
        if not str(sid).startswith("S-base-"):
            continue
        for l in run.arm_lines(sid, "direct-local"):
            base.append([latex.breakable_id(sid), "r%s" % l.get("repeat"), esc(l.get("outcome")),
                         esc(cells.fmt_t(l.get("wall_ms")) if l.get("wall_ms") is not None else "–"),
                         esc(cells.fmt_t(((l.get("server") or {}).get("elapsed_ms")))
                             if isinstance(l.get("server"), dict) and (l.get("server") or {}).get("elapsed_ms")
                             is not None else "–")])
    if base:
        out += latex.plain_table(latex.spec_p([34, 8, 22, 20, 20], aligns=["raggedright", "raggedright",
                                                                           "raggedright", "raggedleft", "raggedleft"]),
                                 ["base request", "", "outcome", "wall", "server"], base,
                                 "the stops' base requests run to completion on direct-local: repeat 1 is the warm "
                                 "full-run time a stop must beat")
    return out


# -- the chapters' prose values ------------------------------------------------------------------------------

def _count_text(counts, order):
    return ", ".join("%d %s" % (counts[k], k.replace("_", " ")) for k in order if counts.get(k)) or "none"


def chapter_values(ctx, ch):
    """The prose values of chapter G (the five-way comparison's counts and words) or S (the stops' records); {} for
    another chapter. Every figure of the run is one of these (the prose carries none)."""
    run = ctx.run
    try:
        sids = ctx.units.chapter_scenarios(ch)
        if not any(isinstance(run.scenario(s).get("graphlet"), dict)
                   or isinstance((run.scenario(s).get("mcp") or {}).get("flow"), dict) for s in sids):
            return {}
    except (AttributeError, TypeError):   # a test's stub run: no v8 chapter
        return {}
    v = {}
    if any(isinstance(run.scenario(s).get("graphlet"), dict) for s in sids):
        g = data(run) or {}
        if g.get("error"):
            v.update(g_error=str(g["error"]), g_verdict="ERROR", g_reason=str(g["error"]))
            return v
        word, reason = verdict(g)
        s = g.get("summary") or {}
        sc_, ops = s.get("scenarios") or {}, s.get("ops") or {}
        checks = g.get("service_checks") or []
        by_outcome = {}
        for c in checks:
            by_outcome[c.get("outcome")] = by_outcome.get(c.get("outcome"), 0) + 1
        cov = g.get("coverage") or {}
        n_checked = sum(1 for t in cov.values() for c in (t.get("checked") or {}).values() if not c.get("rederived"))
        n_red = sum(1 for t in cov.values() for c in (t.get("checked") or {}).values() if c.get("rederived"))
        n_cons = sum(len(t.get("consistency_only") or {}) for t in cov.values())
        dis_tools = sorted(set(e.get("tool") for sv in (g.get("scenarios") or {}).values()
                               for e in (sv.get("ops") or {}).values() if e.get("outcome") == "disagree"))
        pre_bad = sorted(set(k for sv in (g.get("scenarios") or {}).values() for k in PRECHECKS
                             if ((sv.get("precheck") or {}).get(k) or {}).get("outcome") not in ("pass", "not_exposed")))
        nc = not_covered(run)
        lib = ((run.manifest or {}).get("graphlets") or {}).get("library_commit") or ""
        v.update({
            "g_error": "", "g_verdict": word, "g_reason": reason,
            "n_g_scenarios": sum(sc_.values()), "n_g_pass": sc_.get("pass", 0), "n_g_fail": sc_.get("fail", 0),
            "n_g_inconclusive": sc_.get("inconclusive", 0), "n_g_precondition_failed": sc_.get("precondition_failed", 0),
            "g_scenarios_text": _count_text(sc_, ("pass", "fail", "precondition_failed", "inconclusive")),
            "n_g_ops": sum(ops.values()), "n_g_agree": ops.get("agree", 0), "n_g_disagree": ops.get("disagree", 0),
            "n_g_ops_inconclusive": ops.get("inconclusive", 0) + ops.get("precondition_failed", 0),
            "n_g_by_design": ops.get("service_by_design", 0),
            "g_ops_text": _count_text(ops, ("agree", "disagree", "service_by_design", "inconclusive",
                                            "precondition_failed")),
            "g_disagree_tools": ", ".join(dis_tools), "g_prechecks_failed": ", ".join(k.replace("_", " ")
                                                                                       for k in pre_bad),
            "g_service_text": _count_text(by_outcome, ("pass", "fail", "inconclusive", "n/a")),
            "n_g_service_fail": by_outcome.get("fail", 0),
            "n_g_w5_checked": n_checked, "n_g_w5_rederived": n_red, "n_g_w5_consistency": n_cons,
            "g_library_commit": lib[:12], "g_direct_only_ids": facts.ids_ranges([x for x, _ in nc["direct_only"]]),
            "n_g_direct_only": len(nc["direct_only"]), "g_v9_ids": facts.ids_ranges([x for x, _ in nc["v9"] if x]),
        })
    stops = [(l, r) for l, r in stop_rows(run) if l.get("scenario") in sids]
    if stops:
        taken = sum(1 for _, r in stops if r and (r.get("stop") or {}).get("took"))
        released = sum(1 for _, r in stops if r and (r.get("attempts") or [])
                       and all((a.get("release") or {}).get("reason") for a in r.get("attempts") or []))
        v.update({"n_stops": len(stops), "n_stops_taken": taken, "n_stops_released": released,
                  "n_stops_pass": sum(1 for l, _ in stops if l.get("outcome") == "pass")})
    return v
