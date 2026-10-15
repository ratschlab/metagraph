"""The sections of the results document (DESIGN-DOC §1, §4, §6–§8): Part I, the chapter pages, the family sections
with their unit tables and detail blocks, the closing section's tables and the appendices. Every function returns
LaTeX; the data come from a `reader.Run`, a `units.Units` and the prose dicts.
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import blocks  # noqa: E402
import cells  # noqa: E402
import latex  # noqa: E402
import report  # noqa: E402
from cells import ARMS, fmt_t, fmt_pair, fmt_n  # noqa: E402
from latex import esc, esc_wrap, tt, inline, paragraphs  # noqa: E402

RELATION_SIGN = {"prefix_equals": "{\\glyphfont ⊑}", "subset": "{\\glyphfont ⊆}", "same_answer": "=",
                 "zero_equals_absent": "{\\glyphfont ≡}"}
LOAD_CAVEAT_SRA = report.LOAD_CAVEAT_SRA_CPU_PCT
LOAD_CAVEAT_LOAD1 = report.LOAD_CAVEAT_LOAD1
# Part I's sections: key, the default title, and the id the prose file may use instead of the key
PART_ONE = [("what", "What the suite tests and how", "I.1"), ("arms", "The three arms and the paired loop", "I.2"),
            ("classes", "Outcome classes, agreement classes, the verdict rule", "I.3"),
            ("reading", "How to read a block", "I.4"), ("environment", "Versions and environment", "I.5"),
            ("verdict", "The verdict and the deploy gate", "I.6"),
            ("limits", "Known MCP limits and owner-declared exemptions", "I.7")]


class Context(object):
    """What every section reads: the run, the unit model, the prose, the chapters to build, the footnotes assigned
    (per chapter) and the per-scenario arm states (computed once)."""

    def __init__(self, run, units, prose, overview, glossary, chapters, checks_prose=None):
        self.run = run
        self.units = units
        self.prose = prose  # {chapter: prose dict}
        self.overview = overview or {}
        self.glossary = glossary or {}
        self.checks_prose = checks_prose or {}  # doc/prose/checks.json: {common: [...], checks: {kind: {...}}}
        self.chapters = chapters
        self.states = {}
        self.footnotes = {}  # chapter -> Footnotes
        self.blocks_by_chapter = {}

    def state(self, sid):
        if sid not in self.states:
            self.states[sid] = cells.scenario_states(self.run, sid)
        return self.states[sid]

    def footnote(self, sid):
        """The footnote number of a test (assigned when its chapter was built), or None."""
        ch = (self.units.by_scenario.get(sid) or {}).get("chapter")
        fn = self.footnotes.get(ch)
        return fn.get(sid) if fn else None

    def selected_scenarios(self):
        return [sid for sid in self.run.order if (self.units.by_scenario.get(sid) or {}).get("chapter")
                in self.chapters]

    def placeholders(self):
        """The values of the prose placeholders (overview.json `placeholders` names them with their sources)."""
        run = self.run
        manifest = run.manifest
        index = manifest.get("index") or {}
        counts = manifest.get("counts") or {}
        arms = run.run.get("arms") or {}
        hc = run.comp.get("host_cache") or {}
        return {"release": run.run.get("release") or "–", "run_id": run.run_id(),
                "database": index.get("database") or "–", "k": index.get("k") or "–",
                "scenarios": len(run.order), "n_scenarios": len(run.order),
                "manifest_version": manifest.get("manifest_version") or "–",
                "manifest_sha_short": (manifest.get("content_sha256") or "–")[:16],
                "n_records": len((manifest.get("ground_truth") or {}).get("records") or []),
                "n_relations": counts.get("relations", len(run.relations)),
                "server_instance": run.run.get("server_instance") or "–",
                "index_fp": (run.run.get("index_fp") or "–")[:16],
                "local_base": (arms.get("direct-local") or {}).get("base") or "–",
                "remote_base": (arms.get("direct-remote") or {}).get("base") or "–",
                "mcp_base": (arms.get("mcp") or {}).get("base") or run.mcp_run.get("base") or "–",
                "warm_window_s": fmt_n(hc.get("warm_window_s", 60), 0),
                "load_caveat_sra_cpu_pct": fmt_n(LOAD_CAVEAT_SRA, 0),
                "load_caveat_load1": fmt_n(LOAD_CAVEAT_LOAD1, 0)}

    def fill(self, text):
        """Prose placeholders `{release}` etc. filled from the run."""
        values = self.placeholders()

        class _Safe(dict):
            def __missing__(self, key):
                return "{%s}" % key

        try:
            return str(text).format_map(_Safe(values))
        except (ValueError, IndexError):
            return str(text)

    def overview_section(self, key):
        """The overview's entry of a Part I section: `sections` keyed by the section key, or a list of entries
        with an `id` (I.1 ...) or a `key`."""
        sections = self.overview.get("sections") or {}
        if isinstance(sections, dict):
            return sections.get(key) or {}
        sid = dict((k, i) for k, _, i in PART_ONE).get(key)
        for entry in sections:
            if isinstance(entry, dict) and (entry.get("key") == key or entry.get("id") == sid):
                return entry
        return {}

    def section_title(self, key, default):
        return self.overview_section(key).get("title") or default

    def prose_paragraphs(self, section_key):
        items = self.overview_section(section_key).get("paragraphs") or []
        if not items:
            return latex.placeholder("(overview prose for %s: to be written)" % section_key) + "\n\n"
        return paragraphs([self.fill(p) for p in items])


# -- small helpers --------------------------------------------------------------------------------------------

def _kv_table(rows, widths=(40, 134)):
    spec = latex.spec_p(widths)
    return latex.plain_table(spec, ["", "value"], [[k, v] for k, v in rows])


def _stat(d, key="p50"):
    return fmt_t((d or {}).get(key)) if isinstance(d, dict) else "–"


def _stats_cells(d):
    """p50 / p90 / max / n of a summarise() dict."""
    if not isinstance(d, dict) or not d.get("n"):
        return ["–", "–", "–", "0"]
    return [fmt_t(d.get("p50")), fmt_t(d.get("p90")), fmt_t(d.get("max")), str(d.get("n"))]


def _short(text, n=12):
    return esc(str(text)[:n]) if text else "–"


def _date(text):
    return esc(str(text)[:19].replace("T", " ")) if text else "–"


# -- front matter and Part I --------------------------------------------------------------------------------------

def front_matter(ctx):
    run = ctx.run
    manifest = run.manifest
    index = manifest.get("index") or {}
    arms = run.run.get("arms") or {}
    rows = [("run", tt(run.run_id())),
            ("release", tt(run.run.get("release") or "–")),
            ("server instance", tt(run.run.get("server_instance") or "–")),
            ("index", "%s, k = %s, fingerprint %s…" % (tt(index.get("database") or "–"), esc(index.get("k") or "–"),
                                                        tt(_short(run.run.get("index_fp"), 16)))),
            ("manifest", "version %s, %d scenarios, content sha256 %s…" % (
                esc(manifest.get("manifest_version")), len(run.order), tt(_short(manifest.get("content_sha256"), 16)))),
            ("arms", "; ".join("%s %s to %s" % (esc(arm), _date((info or {}).get("started_utc")),
                                               _date((info or {}).get("ended_utc")))
                               for arm, info in sorted(arms.items())) or "–"),
            ("set", "the re-evaluated set (checks commit %s)" % tt(_short(((run.reeval or {}).get("checks") or {}).get(
                "commit"), 12)) if run.reevaluated else "the recorded set")]
    return _kv_table(rows)


def verdict_line(ctx):
    v = ctx.run.verdict or "–"
    word, _, rest = v.partition(" ")
    return "\\textbf{\\large %s}\\ %s" % (esc(word), esc(rest, 400))


LEGEND_GLYPH = {"pass": "\\gpass", "fail": "\\gfail", "error": "\\gerr", "finding": "\\gfind",
                "inconclusive": "\\ginc", "not_implemented": "\\gni", "by_design": "\\gbd", "none": "\\gnone"}
LEGEND_WORD = {"pass": "pass", "fail": "fail", "error": "error", "finding": "finding", "inconclusive": "inconclusive",
               "not_implemented": "not implemented via MCP", "by_design": "service-side by design",
               "none": "no execution"}
MARK_EXAMPLE = {"superscript digit": "\\ginc\\fn{1}", "not-equal sign": "\\gpass\\gne", "t": "\\gpass\\tmark{t}",
                "s": "\\gpass\\tmark{s}", "n": "\\gpass\\tmark{n}", "u": "\\gbd\\tmark{u}", "c": "9.1 s\\tmark{c}",
                "a": "9.1 s\\tmark{a}", "m": "\\ginc\\tmark{m}"}


def legend_table(ctx=None):
    """The glyph legend, the marks and the cell grammar of I.4; the wording of the glyphs and marks comes from the
    overview prose (`legend`, `marks` of the reading section) when it has them, the symbols always from the macros
    the tables use."""
    reading = ctx.overview_section("reading") if ctx else {}
    legend = [e for e in (reading.get("legend") or []) if isinstance(e, dict) and e.get("key") in LEGEND_GLYPH]
    marks = [e for e in (reading.get("marks") or []) if isinstance(e, dict) and e.get("mark") in MARK_EXAMPLE]
    if legend and marks:
        glyph_rows = []
        for e in legend:
            meaning = str(e.get("meaning") or "")
            prefix = LEGEND_WORD[e["key"]] + ":"
            if meaning.startswith(prefix):
                meaning = meaning[len(prefix):].strip()
            glyph_rows.append([LEGEND_GLYPH[e["key"]], esc(LEGEND_WORD[e["key"]]), inline(meaning)])
        mark_rows = [[MARK_EXAMPLE[e["mark"]], esc(e.get("on") or ""), inline(e.get("meaning") or "")]
                     for e in marks]
        return _legend_tables(glyph_rows, mark_rows, with_clean=False)
    return _legend_builtin()


def _legend_tables(glyph_rows, mark_rows, with_clean=True):
    spec = latex.spec_p([8, 36, 128], aligns=["centering", "raggedright", "raggedright"])
    glyphs = latex.plain_table(spec, ["", "glyph", "meaning"], glyph_rows)
    spec2 = latex.spec_p([30, 34, 108])
    marks_t = latex.plain_table(spec2, ["mark", "on", "meaning"], mark_rows)
    grammar = ("\\begingroup\\sffamily\\footnotesize\n"
               "\\bline{line 1}{glyph, marks, wall time(s) in repeat order with the unit once: "
               "\\gpass\\ 3.1 / 0.9 s}\n"
               "\\bline{line 2}{[server elapsed time(s)]: [1.8 s / 22 ms]; [–] when the answer carries no timing "
               "(a refusal, a 503)}\n"
               "\\bline{line 3}{mcp only: the backend time the service measured for the host call: be 1.8 s / 91 ms}"
               "\n\\par\\endgroup\n\n")
    clean = ("A test is \\emph{clean} when every arm is one of: pass; not implemented via MCP; service-side by "
             "design (verified); no execution because the plan had none. Everything else gets a footnote mark on "
             "each offending cell and one detail block after the table: fail, error, finding, inconclusive, mixed "
             "repeats, a direct class disagree / not comparable / same checks / timing, an MCP class disagree / not "
             "comparable, a missing planned execution. Series matrices stack the three lines in one cell; list and "
             "pairs rows put lines 1 and 2 in the cell and line 3 in the mcp cell; replicates rows carry aggregates "
             "(ok n/n, wall p50 / max over each member's first host-warm execution, [server p50]); refusal rows "
             "carry the glyph only, their times as a maximum in the summary line.\n\n") if with_clean else ""
    return glyphs + marks_t + grammar + clean


def _legend_builtin():
    rows = [("\\gpass", "pass", "execution outcome pass on every repeat"),
            ("\\gfail", "fail", "outcome fail on any repeat"),
            ("\\gerr", "error", "no HTTP answer (an error outcome)"),
            ("\\gfind", "finding", "the host broke a stated bound or SPEC promise without failing the request"),
            ("\\ginc", "inconclusive", "the test could not test its property (an allowed 503 deadline, a withheld "
                                       "list, a service timeout); with the mark m the repeats' outcomes differ"),
            ("\\gni", "not implemented via MCP", "MCP agreement class not_implemented_mcp (the manifest's mark, "
                                                 "or a table case of that agreement)"),
            ("\\gbd", "service-side by design", "the service answered or refused in the host's place, as its table "
                                                "documents; the mark u when the pair is unverified"),
            ("\\gnone", "no execution", "the plan had none (a long scenario on direct-remote); a planned execution "
                                        "not recorded is footnoted")]
    glyph_rows = [[g, esc(w), esc(m)] for g, w, m in rows]
    marks = [("\\ginc\\fn{1}", "any cell", "a detail block follows the table (numbering restarts per chapter)"),
             ("\\gpass\\gne", "direct-remote or mcp", "disagreement with the reference arm (always footnoted)"),
             ("\\gpass\\tmark{t}", "direct-remote", "direct class timing: one arm's answer time-limited"),
             ("\\gpass\\tmark{s} / \\gpass\\tmark{n}", "direct-remote",
              "same_checks (other answer, same check outcomes) / not_comparable"),
             ("\\gbd\\tmark{u} / \\gpass\\tmark{n}", "mcp",
              "by-design pair unverified / not_comparable (no check both arms evaluated)"),
             ("9.1 s\\tmark{c} / 9.1 s\\tmark{a}", "a time", "this single execution was host-cold / aged (a direct cell "
                                                            "without a mark is host-warm)"),
             ("\\ginc\\tmark{m}", "a glyph", "the repeats' outcomes differ (the worst is shown)")]
    return _legend_tables(glyph_rows, [[a, esc(b), esc(c)] for a, b, c in marks])


def classes_tables():
    outcomes = [("pass", "every expected check passed"), ("fail", "a check failed"),
                ("error", "no HTTP answer"), ("finding", "a stated bound or SPEC promise broken without a failed "
                                                         "request (budget_respected)"),
                ("inconclusive", "a check could not test its property (an allowed 503, a withheld list, a host state)"),
                ("not_implemented_mcp", "the MCP surface cannot express the request (the manifest's mark)"),
                ("service_by_design", "the service answered or refused in the host's place, as its table says")]
    direct = [("identical", "equal answers (digest)"), ("same_checks", "different answers, the same check outcomes"),
              ("timing", "the answers differ by timing alone (one time-limited)"),
              ("not_comparable", "no check both arms evaluated"), ("disagree", "a check differs")]
    mcp = [("agree", "the checks both arms evaluated agree (the depth says how many property checks)"),
           ("not_implemented_mcp", "by the manifest or the table"),
           ("service_by_design", "a table case matched (verified when both sides match the case)"),
           ("not_comparable", "no check both arms evaluated"), ("disagree", "a check differs"),
           ("missing", "no MCP execution recorded")]
    spec = latex.spec_p([40, 134])
    out = latex.plain_table(spec, ["execution outcome", "meaning"], [[tt(k), esc(v)] for k, v in outcomes])
    out += latex.plain_table(spec, ["direct agreement class", "meaning"], [[tt(k), esc(v)] for k, v in direct])
    out += latex.plain_table(spec, ["MCP agreement class", "meaning"], [[tt(k), esc(v)] for k, v in mcp])
    return out


def environment_table(ctx):
    run = ctx.run
    r = run.run
    manifest = run.manifest
    caps = (r.get("capabilities") or {}).get("/capabilities") or {}
    reeval = run.reeval or {}
    checks = reeval.get("checks") or {}
    table = run.table_meta or {}
    comp = run.comp
    session = run.session()
    rows = [("run id", tt(run.run_id())),
            ("engine release", tt(r.get("release") or "–")),
            ("server instance", tt(r.get("server_instance") or "–")),
            ("index fingerprint", tt(r.get("index_fp") or "not stated by the host")),
            ("feature level / algorithm", "%s / %s" % (esc(caps.get("feature_level") or "–"),
                                                        esc(caps.get("algorithm_version") or "–"))),
            ("manifest", "version %s; content sha256 %s; relations sha256 %s; %d scenarios, %d relations" % (
                esc(manifest.get("manifest_version")), tt(manifest.get("content_sha256") or "–"),
                tt(manifest.get("relations_sha256") or "–"), len(run.order), len(run.relations)))]
    if run.reevaluated:
        rows.append(("checks (re-evaluation)", "commit %s%s; checks.py sha256 %s; suite\\_common sha256 %s; %s" % (
            tt(_short(checks.get("commit"), 12)), " (dirty)" if checks.get("dirty") else "",
            tt(_short(checks.get("checks_sha256"), 16)), tt(_short(checks.get("suite_common_sha256"), 16)),
            _date(reeval.get("reevaluated_utc")))))
    else:
        rows.append(("checks", "the recorded set: judged by the runners' checks (%s)" % esc(
            (run.judged or {}).get("source") or "–")))
    if table:
        rows.append(("by-design table", "version %s, %s cases, service commit %s, generated %s" % (
            esc(table.get("version")), esc(table.get("cases")), tt(_short(table.get("service_commit"), 12)),
            _date(table.get("generated_utc")))))
    mcp_run = run.mcp_run or {}
    for arm, info in sorted((r.get("arms") or {}).items()):
        info = info or {}
        ended = info.get("ended_utc") or (mcp_run.get("finished_utc") if arm == "mcp" else None)
        scripts = info.get("scripts_commit") or (mcp_run.get("git_head") if arm == "mcp" else None)
        rows.append(("arm %s" % esc(arm), "%s; state %s; %s executions recorded; runner %s, scripts %s; %s to %s" % (
            tt(info.get("base") or "–"), esc(info.get("state") or "–"), esc(info.get("executions_recorded", "–")),
            esc(info.get("runner_host") or "–"), tt(_short(scripts, 12)),
            _date(info.get("started_utc")), _date(ended))))
    if mcp_run:
        rows.append(("MCP runner", "git head %s; base %s; admin key %s; the service's deployed version is not "
                                   "recorded in a v4 run" % (
                                       tt(_short(mcp_run.get("git_head"), 12)), tt(mcp_run.get("base") or "–"),
                                       esc(mcp_run.get("admin_key") or "–"))))
    if session:
        loop = session.get("loop") or {}
        rows.append(("mode", "%s; loop %s (%s of %s scenarios); started %s; mex clock offset %s s" % (
            esc(session.get("mode") or "–"), esc(loop.get("state") or "–"), esc(loop.get("completed", "–")),
            esc(loop.get("scenarios", "–")), _date(session.get("started_utc")),
            esc(session.get("mex_clock_offset_s", "–")))))
    hl_s, hl_e = r.get("host_load_start"), r.get("host_load_end")
    rows.append(("host load (load1 / 5 / 15)", "at start %s; at end %s" % (
        esc(", ".join(fmt_n(x, 1) for x in hl_s)) if hl_s else "–",
        esc(", ".join(fmt_n(x, 1) for x in hl_e)) if hl_e else "–")))
    lf = run.load or {}
    caveat = ""
    if (lf.get("sra_mean") or 0) > LOAD_CAVEAT_SRA or (lf.get("load_mean") or 0) > LOAD_CAVEAT_LOAD1:
        caveat = (" — the Logan caveat applies: the host was busy with other work (sra\\_all above %d\\%% CPU or "
                  "load1 above %d), so the timings are upper bounds" % (LOAD_CAVEAT_SRA, LOAD_CAVEAT_LOAD1))
    rows.append(("samples", "%d; sra\\_all CPU mean %s\\%%, max %s\\%%; load1 mean %s, max %s%s" % (
        lf.get("n", 0), fmt_n(lf.get("sra_mean")), fmt_n(lf.get("sra_max")), fmt_n(lf.get("load_mean"), 1),
        fmt_n(lf.get("load_max"), 1), caveat)))
    hc = comp.get("host_cache") or {}
    rows.append(("host cache", "cold = the scenario's first execution on the host in this run; warm = another "
                               "within %s s; aged = older; the MCP pair runs cold then warm, the direct arms between "
                               "them" % esc(hc.get("warm_window_s", "60"))))
    return _kv_table(rows)


def verdict_and_gate(ctx):
    run = ctx.run
    totals = run.summary().get("totals") or {}
    out = "\\begin{quote}%s\\end{quote}\n\n" % verdict_line(ctx)
    rows = []
    for arm in ARMS:
        t = totals.get(arm) or {}
        o = t.get("outcomes") or {}
        hc = t.get("host_cache") or {}
        rows.append([esc(arm), str(t.get("executions", 0)), str(o.get("pass", 0)), str(o.get("fail", 0)),
                     str(o.get("error", 0)), str(o.get("finding", 0)), str(o.get("inconclusive", 0)),
                     str(o.get("not_implemented_mcp", 0)), str(o.get("service_by_design", 0)),
                     "%s / %s / %s" % (hc.get("cold", 0), hc.get("warm", 0), hc.get("not_reached", 0))])
    spec = "@{}l r r r r r r r r r@{}"
    out += latex.plain_table(spec, ["arm", "executions", "\\gpass", "\\gfail", "\\gerr", "\\gfind", "\\ginc", "\\gni",
                                    "\\gbd", "cold / warm / not reached"], rows)
    gate = run.gate
    if isinstance(gate, dict):
        pol = gate.get("policy") or {}
        word = "GATE PASS" if gate.get("passed") else "GATE FAIL"
        out += "\\textbf{%s}: policy %s version %s (sha256 %s…), judged %s.\n\n" % (
            word, esc(pol.get("name") or "–"), esc(pol.get("policy_version") or "–"), tt(_short(pol.get("sha256"), 12)),
            _date(gate.get("judged_utc")))
        problems = [p for p in (gate.get("problems") or []) if isinstance(p, dict)]
        if problems:
            by_kind = {}
            for p in problems:
                by_kind.setdefault(p.get("kind") or "?", []).append(p)
            rows = []
            for k, v in sorted(by_kind.items()):
                names = sorted(set("%s %s" % (p.get("scenario"), p.get("arm")) for p in v))
                rows.append([esc(k), str(len(v)), esc(", ".join(names[:20]), 400)])
            out += latex.plain_table(latex.spec_p([40, 124], extra_cols=1, extra_spec="r"),
                                     ["problem kind", "scenario, arm", "n"], [[r[0], r[2], r[1]] for r in rows])
    else:
        out += "No gate verdict is recorded for this run (no gate.json)."
        policy = run.policy
        if isinstance(policy, dict):
            counts = policy.get("counts") or {}
            out += " The policy %s (version %s) would select: %s." % (
                esc(policy.get("name") or "–"), esc(policy.get("policy_version") or "–"),
                esc(", ".join("%s %s" % (k, v) for k, v in sorted(counts.items()) if not isinstance(v, dict)), 300))
        out += "\n\n"
    return out


def mcp_limits(ctx):
    run = ctx.run
    out = ""
    # (a) the table's not-implemented cases that scenarios of the run name
    named = {}
    for sid in run.order:
        mcp = run.scenario(sid).get("mcp") or {}
        bd = mcp.get("by_design")
        for case in ([bd] if isinstance(bd, str) else (bd or [])):
            c = run.cases.get(case) or run.table.get(case) or {}
            if c.get("agreement") == "not_implemented_mcp":
                named.setdefault(case, []).append(sid)
    rows = []
    for case, sids in sorted(named.items()):
        c = run.cases.get(case) or run.table.get(case) or {}
        sa = c.get("service_answer") or {}
        rows.append([tt(case), esc_wrap(c.get("scenario_match"), 160),
                     esc_wrap(sa.get("missing_argument") or sa.get("detail"), 160), str(len(sids))])
    out += "\\paragraph{Cases of the by-design table of agreement not implemented via MCP that scenarios name}\n"
    out += latex.plain_table(latex.spec_p([50, 60, 52], extra_cols=1, extra_spec="r"),
                             ["case", "scenario match", "missing argument", "scenarios"], rows)
    # (b) the manifest's not_implemented reasons
    reasons = {}
    for sid in run.order:
        mcp = run.scenario(sid).get("mcp") or {}
        if mcp.get("not_implemented"):
            reasons.setdefault(str(mcp["not_implemented"]), []).append(sid)
    out += "\\paragraph{The manifest's reasons for not implemented via MCP}\n"
    rows = [[esc_wrap(r, 120), str(len(s)), esc(", ".join(s[:8]) + (" …" if len(s) > 8 else ""))]
            for r, s in sorted(reasons.items(), key=lambda kv: (-len(kv[1]), kv[0]))]
    out += latex.plain_table(latex.spec_p([60, 104], extra_cols=1, extra_spec="r"),
                             ["reason", "ids", "scenarios"], [[r[0], r[2], r[1]] for r in rows])
    # (c) the owner-declared file
    dec = run.declared
    out += "\\paragraph{Owner-declared exemptions}\n"
    if isinstance(dec, dict):
        out += "File %s (declared version %s).\n" % (tt(dec.get("_file") or "declared"), esc(dec.get("declared_version")))
        rows = []
        for kind in ("exclusions", "not_required", "relations_excluded"):
            for e in dec.get(kind) or []:
                if isinstance(e, dict):
                    rows.append([esc(kind.replace("_", " ")), latex.breakable_id(e.get("scenario") or e.get("relation") or "–"),
                                 esc(", ".join(e.get("arms") or [])), esc_wrap(e.get("reason"), 300)])
        out += latex.plain_table(latex.spec_p([22, 30, 22, 96]), ["kind", "scenario", "arms", "reason"], rows)
    else:
        out += "None given for this manifest.\n\n"
    return out


def part_one(ctx):
    run = ctx.run
    out = "\\part*{Part I\\qd Overview}\\addcontentsline{toc}{part}{Part I\\qd Overview}\n\\markright{Overview}\n"
    out += latex.section("I.1\\qd %s" % esc(ctx.section_title("what", PART_ONE[0][1])), "sec:I1")
    out += ctx.prose_paragraphs("what")
    manifest = run.manifest
    index = manifest.get("index") or {}
    counts = manifest.get("counts") or {}
    by_area = counts.get("by_area") or {}
    rows = [("index", "%s (namespace %s), k = %s" % (tt(index.get("database") or "–"), tt(index.get("index_ns") or "–"),
                                                      esc(index.get("k") or "–"))),
            ("scenarios", "%d: %s" % (len(run.order), esc(", ".join("%s %s" % (a, n) for a, n in sorted(
                by_area.items(), key=lambda kv: (ctx.units.chapters.index(kv[0]) if kv[0] in ctx.units.chapters
                                                 else 99, kv[0])))) or "–")),
            ("relations", "%d (%s)" % (len(run.relations), esc(", ".join("%s %s" % kv for kv in sorted(
                (counts.get("relations_by_type") or {}).items()))) or "–"))]
    gt = manifest.get("ground_truth") or {}
    for key in sorted(gt):
        value = gt[key]
        if isinstance(value, (list, dict)):
            value = "%d entries" % len(value)
        rows.append((esc(key.replace("_", " ")), esc_wrap(value, 220)))
    out += _kv_table(rows)
    out += latex.section("I.2\\qd %s" % esc(ctx.section_title("arms", PART_ONE[1][1])), "sec:I2")
    out += ctx.prose_paragraphs("arms")
    arms = run.run.get("arms") or {}
    rows = []
    for arm in ARMS:
        info = arms.get(arm) or {}
        rows.append([esc(arm), tt(info.get("base") or (run.mcp_run.get("base") if arm == "mcp" else None) or "–"),
                     esc(info.get("state") or "–"), str(info.get("executions_recorded", "–")),
                     esc(info.get("connection") or "–", 90)])
    out += latex.plain_table(latex.spec_p([24, 62, 16, 52], extra_cols=1, extra_spec="r"),
                             ["arm", "base", "state", "connection", "executions"],
                             [[r[0], r[1], r[2], r[4], r[3]] for r in rows])
    out += latex.section("I.3\\qd %s" % esc(ctx.section_title("classes", PART_ONE[2][1])), "sec:I3")
    out += ctx.prose_paragraphs("classes")
    out += classes_tables()
    out += latex.section("I.4\\qd %s" % esc(ctx.section_title("reading", PART_ONE[3][1])), "sec:I4")
    out += ctx.prose_paragraphs("reading")
    out += legend_table(ctx)
    out += latex.section("I.5\\qd %s" % esc(ctx.section_title("environment", PART_ONE[4][1])), "sec:I5")
    out += ctx.prose_paragraphs("environment")
    out += environment_table(ctx)
    out += latex.section("I.6\\qd %s" % esc(ctx.section_title("verdict", PART_ONE[5][1])), "sec:I6")
    out += ctx.prose_paragraphs("verdict")
    out += verdict_and_gate(ctx)
    out += latex.section("I.7\\qd %s" % esc(ctx.section_title("limits", PART_ONE[6][1])), "sec:I7")
    out += ctx.prose_paragraphs("limits")
    out += mcp_limits(ctx)
    return out


# -- chapters ----------------------------------------------------------------------------------------------------

def _agreement_counts(ctx, sids):
    run = ctx.run
    direct, mcp = {}, {}
    full = 0
    verified, unverified = 0, 0
    for sid in sids:
        e = run.comp_entry(sid)
        dc = (e.get("direct_agreement") or {}).get("class")
        mc = (e.get("mcp_agreement") or {}).get("class")
        if dc:
            direct[dc] = direct.get(dc, 0) + 1
        if mc:
            mcp[mc] = mcp.get(mc, 0) + 1
            if mc == "service_by_design":
                if (e.get("mcp_agreement") or {}).get("verified") is False:
                    unverified += 1
                else:
                    verified += 1
        d = (run.depth or {}).get(sid)
        if d and d.get("class") == "full":
            full += 1
    dtext = ", ".join("%d %s" % (direct[k], k.replace("_", " ")) for k in
                      ("identical", "same_checks", "timing", "not_comparable", "disagree") if direct.get(k)) or "–"
    mtext = "%d agree (fully checked %d), %d not implemented, %d by design (%d verified, %d unverified), %d not " \
            "comparable, %d disagree%s" % (
                mcp.get("agree", 0), full, mcp.get("not_implemented_mcp", 0), mcp.get("service_by_design", 0),
                verified, unverified, mcp.get("not_comparable", 0), mcp.get("disagree", 0),
                (", %d missing" % mcp["missing"]) if mcp.get("missing") else "")
    return esc(dtext), esc(mtext)


def chapter_summary(ctx, chapter, sids, n_units, n_families):
    run = ctx.run
    states = [ctx.state(sid) for sid in sids]
    out = "\\begingroup\\sffamily\\footnotesize\n"
    out += "\\bline{chapter %s}{%d tests · %d units in %d families}\n" % (esc(chapter), len(sids), n_units, n_families)
    out += "\\bline{outcomes}{%s}\n" % " · ".join("%s %s" % (esc(arm), cells.glyph_counts(states, arm)) for arm in ARMS)
    dtext, mtext = _agreement_counts(ctx, sids)
    out += "\\bline{agreement}{direct: %s · mcp: %s}\n" % (dtext, mtext)
    timing = (run.summary().get("timing") or {}).get(chapter) or {}
    cold, warm, server = [], [], []
    for arm in ARMS:
        t = (timing.get(arm) or {}).get("by_host_cache") or {}
        c = (t.get("cold") or {}).get("wall") or {}
        w = (t.get("warm") or {}).get("wall") or {}
        s = (t.get("warm") or {}).get("server") or {}
        if c.get("n"):
            cold.append("%s %s (n %d)" % (esc(arm), fmt_pair(c.get("p50"), c.get("p90")), c["n"]))
        if w.get("n"):
            warm.append("%s %s (n %d)" % (esc(arm), fmt_pair(w.get("p50"), w.get("p90")), w["n"]))
        server.append(fmt_t(s.get("p50")) if s.get("n") else "–")
    backends = [l.get("mcp_phases", {}).get("backend_ms") for sid in sids for l in run.arm_lines(sid, "mcp")
                if isinstance(l.get("mcp_phases"), dict) and cells.host_cache_of(l) == "warm"]
    be = cells.percentile(backends, 50)
    out += "\\bline{timing}{host-cold wall p50 / p90: %s; host-warm wall p50 / p90: %s; server p50 (host-warm): %s" \
           "%s}\n" % (", ".join(cold) or "–", ", ".join(warm) or "–", " · ".join(server),
                      (" (mcp backend p50 %s)" % fmt_t(be)) if be is not None else "")
    res = (run.summary().get("resources") or {}).get(chapter) or {}
    parts = []
    for arm in ARMS:
        r = res.get(arm) or {}
        if r:
            parts.append("%s: container peak %s\\%% CPU, %s; sra\\_all mean %s\\%%, load1 mean %s" % (
                esc(arm), fmt_n(r.get("peak_cpu_pct")), cells.fmt_gib(r.get("peak_mem_bytes")),
                fmt_n(r.get("mean_sra_all_cpu_pct")), fmt_n(r.get("mean_host_load1"), 1)))
    lf = run.load or {}
    caveat = " (the Logan caveat applies)" if ((lf.get("sra_mean") or 0) > LOAD_CAVEAT_SRA or
                                              (lf.get("load_mean") or 0) > LOAD_CAVEAT_LOAD1) else ""
    out += "\\bline{load}{%s%s}\n" % ("; ".join(parts) or "no samples", caveat)
    fn = ctx.footnotes.get(chapter)
    ids = fn.order if fn else []
    out += "\\bline{detail blocks}{%d%s}\n" % (len(ids), (": " + ", ".join(latex.sid(s, detail=True) for s in ids))
                                                if ids else "")
    out += "\\par\\endgroup\n\n"
    return out


def family_summary(ctx, fid, sids, relations):
    run = ctx.run
    states = [ctx.state(sid) for sid in sids]
    text = "%s: %d test%s · %s · %s" % (tt(fid), len(sids), "" if len(sids) == 1 else "s",
                                        " · ".join("%s %s" % (esc(arm), cells.glyph_counts(states, arm))
                                                   for arm in ARMS), cells.timing_text(run, sids))
    if relations:
        comp_rel = run.comp_relations()
        parts = []
        for arm in ARMS:
            counts = {}
            for rid in relations:
                o = ((comp_rel.get(rid) or {}).get("arms") or {}).get(arm, {}).get("outcome")
                counts[o] = counts.get(o, 0) + 1
            parts.append("%s %s" % (esc(arm), " ".join("%d %s" % (counts[k], cells.GLYPH.get(k, "\\ginc"))
                                                        for k in ("pass", "fail", "inconclusive") if counts.get(k))
                                    or "–"))
        text += " · relations: " + ", ".join(parts)
    return "\\sumline{%s}\n" % text


def _cell_tex(st, footnote, layout):
    l1, l2, l3 = cells.cell_parts(st, footnote)
    if layout == "series":
        return latex.cell_lines(l1, l2, l3)
    if layout in ("list", "pairs"):
        text = l1 + ((" " + l2) if l2 else "")
        if l3:
            text += "\\newline{\\scriptsize\\color{black!60}%s}" % l3
        return text
    return l1


def _footnotes_for(ctx, chapter, sid):
    states = ctx.state(sid)
    if cells.is_clean(states):
        return None
    return ctx.footnotes[chapter].assign(sid)


def _unit_blocks(ctx, chapter, uid, members):
    out = ""
    for sid in members:
        n = ctx.footnotes[chapter].get(sid)
        if n is None:
            continue
        out += blocks.detail_block(ctx.run, sid, ctx.state(sid), n, ctx.units)
        ctx.blocks_by_chapter.setdefault(chapter, []).append(sid)
    return out


def _relation_rows(ctx, relations):
    run = ctx.run
    comp_rel = run.comp_relations()
    man_rel = dict((r.get("id"), r) for r in run.relations)
    rows = []
    for rid in relations:
        rel = man_rel.get(rid) or comp_rel.get(rid) or {}
        crel = comp_rel.get(rid) or {}
        sign = RELATION_SIGN.get(rel.get("type"), esc(rel.get("type") or "?"))
        desc = "%s %s %s: %s" % (latex.sid(rel.get("a"), ctx.units.unit_of(rel.get("a"))), sign,
                                 latex.sid(rel.get("b"), ctx.units.unit_of(rel.get("b"))), esc(rel.get("title"), 90))
        cs = []
        for arm in ARMS:
            a = (crel.get("arms") or {}).get(arm) or {}
            if not a:
                cs.append("\\gnone")
                continue
            o = a.get("outcome")
            glyph = cells.GLYPH.get(o, "\\ginc") if o in cells.GLYPH else "\\ginc"
            detail = a.get("detail")
            if o == "pass" or not detail:
                cs.append(glyph)
            else:
                cs.append("%s {\\scriptsize %s}" % (glyph, esc(detail, 70)))
        rows.append([latex.breakable_id(rid), desc] + cs)
    return rows


def unit_section(ctx, chapter, uid):
    run, units = ctx.run, ctx.units
    u = units.units[uid]
    members = list(u["members"])
    prose = ctx.prose.get(chapter) or {}
    prop = (prose.get("units") or {}).get(uid)
    prop_tex = inline(prop) if prop else latex.placeholder("(unit property: to be written)")
    out = latex.unit_header(uid, esc(u.get("unit_title") or ""), u["unit_type"], prop_tex)
    layout = u["unit_type"]
    fns = dict((sid, _footnotes_for(ctx, chapter, sid)) for sid in members)
    meta = units.by_scenario
    if layout == "series":
        columns = [esc(meta[sid].get("variant") or sid) for sid in members]
        rows = []
        for arm in ARMS:
            rows.append((esc(arm), [_cell_tex(ctx.state(sid)[arm], fns[sid], "series") for sid in members]))
        out += latex.series_table(columns, rows, units_ids(ctx, members))
    elif layout == "replicates":
        groups = u.get("groups")
        rows = []
        by_group = {}
        if groups:
            for g in groups:
                gm = [sid for sid in members if meta[sid].get("group") == g]
                if gm:
                    by_group[g] = gm
        else:
            variants = [str(meta[sid].get("variant") or sid) for sid in members]
            by_group[esc(", ".join(variants), 90)] = members
        for g, gm in by_group.items():
            st = dict((sid, ctx.state(sid)) for sid in gm)
            cs = []
            for arm in ARMS:
                l1, l2 = cells.aggregate(st, arm)
                cs.append(latex.cell_lines(l1, l2))
            rows.append((esc(g) if groups else g, len(gm), cs[0], cs[1], cs[2]))
        ids = units_ids(ctx, members, by_group if groups else None)
        out += latex.replicates_table("variants" if not groups else esc(u.get("unit_title") or uid)[:30], rows, ids)
        exceptions = []
        for sid in members:
            if fns[sid] is not None:
                st = ctx.state(sid)
                arms = [a for a in ARMS if not st[a].clean]
                worst = cells.worst_of([st[a].worst for a in arms if st[a].present]) if arms else None
                exceptions.append("%s %s\\fn{%d} (%s)" % (latex.sid(sid, detail=True),
                                                          cells.GLYPH.get(worst, "\\gnone"), fns[sid],
                                                          esc(", ".join(arms))))
        out += "\\noindent{\\sffamily\\footnotesize exceptions: %s}\\par\\smallskip\n" % (
            "; ".join(exceptions) if exceptions else "none")
    elif layout == "refusals":
        rows = []
        prev = None
        for sid in members:
            st = ctx.state(sid)
            expected, prev = blocks.refusal_expected(run.scenario(sid), prev)
            cs = [st[arm].glyph_tex(fns[sid]) for arm in ARMS]
            rows.append([latex.sid(sid, detail=fns[sid] is not None) if fns[sid] else latex.breakable_id(sid),
                         esc(meta[sid].get("variant") or "", 60), expected] + cs)
        out += latex.refusals_table(rows, None)
        out += "\\noindent{\\scriptsize %s}\\par\\smallskip\n" % units_ids(ctx, members)
        states = [ctx.state(sid) for sid in members]
        out += "\\sumline{%d tests · %s · %s}\n" % (len(members), " · ".join(
            "%s %s" % (esc(arm), cells.glyph_counts(states, arm)) for arm in ARMS), cells.max_wall_text(run, members))
    elif layout == "list":
        rows = []
        for sid in members:
            st = ctx.state(sid)
            cs = [_cell_tex(st[arm], fns[sid], "list") for arm in ARMS]
            rows.append([latex.sid(sid, detail=True) if fns[sid] else latex.breakable_id(sid),
                         esc(meta[sid].get("variant") or "", 60),
                         blocks.expected_text(run.scenario(sid), run.is_property)] + cs)
        out += latex.list_table(rows, None)
    else:  # pairs
        rows = []
        for sid in members:
            st = ctx.state(sid)
            cs = []
            for arm in ARMS:
                text = _cell_tex(st[arm], fns[sid], "pairs")
                rep = run.scenario(sid).get("repeat")
                if isinstance(rep, int) and rep > 1 and st[arm].present:
                    de = run.comp_arm(sid, arm).get("digests_equal")
                    if de is True:
                        text += " ="
                    elif de is False:
                        text += " \\gne"
                cs.append(text)
            second = esc(meta[sid].get("variant") or "", 70)
            src = run.scenario(sid).get("derived_from")
            if src:
                second = "%s · %s" % (latex.sid(src, units.unit_of(src)), second)
            rows.append([latex.sid(sid, detail=True) if fns[sid] else latex.breakable_id(sid), second] + cs)
        rows += _relation_rows(ctx, u.get("relations") or [])
        out += latex.pairs_table(rows, None, "test")
    if layout in ("list", "pairs"):
        out += "\\noindent{\\scriptsize %s}\\par\\smallskip\n" % units_ids(ctx, members)
    out += _unit_blocks(ctx, chapter, uid, members)
    return out


def units_ids(ctx, members, by_group=None):
    import units as units_mod
    return units_mod.ids_text(members, None, by_group)


def family_section(ctx, chapter, fid):
    units = ctx.units
    f = units.families[fid]
    prose = ctx.prose.get(chapter) or {}
    fp = (prose.get("families") or {}).get(fid) or {}
    title = fp.get("title") or f.get("title") or fid
    out = latex.section("%s\\qd %s" % (latex.title_tt(fid), esc(title)), "fam:%s" % fid, level="subsection")
    out += inline(fp["prose"]) + "\n\n" if fp.get("prose") else latex.placeholder("(family prose: to be written)") + "\n\n"
    sids, relations = [], []
    for uid in f["units"]:
        out += unit_section(ctx, chapter, uid)
        sids.extend(units.units[uid]["members"])
        relations.extend(units.units[uid].get("relations") or [])
    out += family_summary(ctx, fid, sids, relations)
    return out


def chapter(ctx, ch):
    units = ctx.units
    ctx.footnotes[ch] = cells.Footnotes()
    prose = ctx.prose.get(ch) or {}
    title = prose.get("title") or units.titles.get(ch, ch)
    out = latex.section("%s\\qd %s" % (esc(ch), esc(title)), "ch:%s" % ch, mark="%s · %s" % (ch, title))
    out += paragraphs(prose.get("intro") or []) if prose.get("intro") else latex.placeholder(
        "(chapter prose: to be written)") + "\n\n"
    fams = units.chapter_families(ch)
    sids = units.chapter_scenarios(ch)
    body = ""
    for fid in fams:
        body += family_section(ctx, ch, fid)
    n_units = sum(len(units.families[f]["units"]) for f in fams)
    # the summary is built after the families so the footnotes are assigned
    out += chapter_summary(ctx, ch, sids, n_units, len(fams))
    return out + body


def part_two(ctx):
    out = "\\part*{Part II\\qd Results by area}\\addcontentsline{toc}{part}{Part II\\qd Results by area}\n"
    for ch in ctx.chapters:
        out += chapter(ctx, ch)
    return out


# -- the closing section -------------------------------------------------------------------------------------------

def closing(ctx):
    s = ctx.run.summary()
    out = "\\part*{Closing\\qd Cross-cutting timing}\\addcontentsline{toc}{part}{Closing\\qd Cross-cutting timing}\n"
    out += "\\markright{Cross-cutting timing}\n"
    out += latex.section("C.1\\qd MCP overhead decomposition", "sec:C1")
    out += ("Pairs of host-warm executions only: the MCP warm execution (repeat 1 preferred) minus each direct arm's "
            "first warm execution; the phases as the service measured them.\n\n")
    ph = s.get("mcp_phases") or {}
    rows = [[esc(k.replace("_ms", "").replace("_", " "))] + _stats_cells(ph.get(k)) for k in
            ("submit_ms", "queue_ms", "backend_ms", "backend_pattern_ms", "terminal_lag_ms", "paging_ms") if k in ph]
    out += latex.plain_table("@{}l r r r r@{}", ["phase", "p50", "p90", "max", "n"], rows, "MCP phases")
    mo = s.get("mcp_overhead") or {}
    rows = [[esc(k.replace("_ms", "").replace("_", " "))] + _stats_cells(mo.get(k)) for k in sorted(mo)]
    out += latex.plain_table("@{}l r r r r@{}", ["MCP minus direct", "p50", "p90", "max", "n"], rows,
                             "MCP overhead (warm MCP wall minus each direct arm; the backend time minus each direct arm)")
    if s.get("mcp_cold_wall_ms"):
        out += latex.plain_table("@{}l r r r r@{}", ["", "p50", "p90", "max", "n"],
                                 [["MCP cold wall"] + _stats_cells(s["mcp_cold_wall_ms"])], "MCP cold wall")
    out += latex.section("C.2\\qd Front cost", "sec:C2")
    fo = s.get("front_overhead") or {}
    rows = [[esc(arm)] + _stats_cells((fo.get("connect_ms") or {}).get(arm)) for arm in sorted(fo.get("connect_ms") or {})]
    out += latex.plain_table("@{}l r r r r@{}", ["connection set-up", "p50", "p90", "max", "n"], rows,
                             "connect\\_ms per arm")
    pan = fo.get("per_area_net") or {}
    rows = [[esc(a)] + _stats_cells(pan.get(a)) for a in sorted(pan, key=lambda x: (ctx.units.chapters.index(x)
                                                                                    if x in ctx.units.chapters else 99, x))]
    out += latex.plain_table("@{}l r r r r@{}", ["chapter", "p50", "p90", "max", "n"], rows,
                             "direct-remote minus direct-local, net of the server time, per chapter")
    cold = fo.get("cold") or {}
    rows = []
    for arm in sorted(cold):
        for kind in ("wall_ms", "server_ms"):
            if kind in cold[arm]:
                rows.append([esc(arm), esc(kind.replace("_ms", ""))] + _stats_cells(cold[arm][kind]))
    out += latex.plain_table("@{}l l r r r r@{}", ["arm", "", "p50", "p90", "max", "n"], rows, "host-cold executions")
    out += latex.section("C.3\\qd Cold versus warm", "sec:C3")
    cw = s.get("cold_warm") or {}
    rows = []
    for arm in sorted(cw):
        d = cw[arm] or {}
        rows.append([esc(arm), str(d.get("n", 0)), fmt_t((d.get("delta_ms") or {}).get("p50")),
                     fmt_t((d.get("delta_ms") or {}).get("p90")), fmt_n((d.get("ratio") or {}).get("p50"), 2),
                     fmt_n((d.get("ratio") or {}).get("p90"), 2)])
    out += latex.plain_table("@{}l r r r r r@{}", ["arm", "pairs", "delta p50", "delta p90", "ratio p50", "ratio p90"],
                             rows, "warm minus cold per scenario (delta), warm over cold (ratio)")
    out += latex.section("C.4\\qd Load telemetry", "sec:C4")
    lf = ctx.run.load or {}
    out += ("%d samples; sra\\_all CPU mean %s\\%% (max %s\\%%), load1 mean %s (max %s). The Logan caveat: above "
            "%d\\%% sra\\_all CPU or load1 %d the host was busy with other work and the timings are upper bounds.\n\n" % (
                lf.get("n", 0), fmt_n(lf.get("sra_mean")), fmt_n(lf.get("sra_max")), fmt_n(lf.get("load_mean"), 1),
                fmt_n(lf.get("load_max"), 1), LOAD_CAVEAT_SRA, LOAD_CAVEAT_LOAD1))
    res = s.get("resources") or {}
    totals = s.get("totals") or {}
    rows = []
    areas = sorted(res, key=lambda x: (ctx.units.chapters.index(x) if x in ctx.units.chapters else 99, x))
    for area in areas + ["all"]:
        src = res.get(area) if area != "all" else dict((a, (totals.get(a) or {}).get("resources")) for a in ARMS)
        for arm in ARMS:
            r = (src or {}).get(arm) or {}
            if not r:
                continue
            rows.append([esc(area), esc(arm), str(r.get("executions", "–")), str(r.get("executions_with_samples", "–")),
                         fmt_n(r.get("peak_cpu_pct")), cells.fmt_gib(r.get("peak_mem_bytes")),
                         fmt_n(r.get("mean_sra_all_cpu_pct")), fmt_n(r.get("mean_host_load1"), 1)])
    out += latex.plain_table("@{}l l r r r r r r@{}", ["chapter", "arm", "executions", "with samples", "peak CPU \\%",
                                                        "peak memory", "sra\\_all mean \\%", "load1 mean"], rows,
                             "the container's peak CPU and memory, the host's sra\\_all CPU and load1 over the executions")
    out += latex.section("C.5\\qd Per-chapter timing", "sec:C5")
    timing = s.get("timing") or {}
    rows = []
    wai = 0
    for area in sorted(timing, key=lambda x: (ctx.units.chapters.index(x) if x in ctx.units.chapters else 99, x)):
        for arm in ARMS:
            t = (timing[area].get(arm) or {}).get("by_host_cache") or {}
            for hc in ("cold", "warm"):
                b = t.get(hc) or {}
                w, sv = b.get("wall") or {}, b.get("server") or {}
                if not w.get("n") and not sv.get("n"):
                    continue
                rows.append([esc(area), esc(arm), esc(hc)] + _stats_cells(w) + _stats_cells(sv))
            wai += ((t.get(cells.WARM_AFTER_INTERRUPTED) or {}).get("wall") or {}).get("n") or 0
    out += latex.plain_table("@{}l l l r r r r r r r r@{}", ["chapter", "arm", "host", "wall p50", "p90", "max", "n",
                                                              "server p50", "p90", "max", "n"], rows,
                             "wall and server time per chapter, arm and host-cache state")
    if wai:
        out += ("%d execution(s) ran warm after an interrupted execution and are counted in neither bucket.\n\n" % wai)
    return out


# -- appendices --------------------------------------------------------------------------------------------------

def _catalogue_entry(ctx, name, entry, chapters_text, generic):
    """One entry of the check catalogue from checks.json: the title line, the SPEC references, then asserts /
    inconclusive / fails / args and the per-route rules as hanging lines."""
    title = inline(entry.get("title") or "")
    out = "\\Needspace{8\\baselineskip}\\noindent\\textbf{\\texttt{%s}}\\quad %s\\hfill{\\small\\itshape %s check · %s}\\par\\nopagebreak\n" % (
        esc(name), title, "generic" if generic else "property", chapters_text)
    spec = entry.get("spec")
    if isinstance(spec, list) and spec:
        out += "{\\small\\color{black!60}%s}\\par\\nopagebreak\n" % inline("; ".join(str(x) for x in spec))
    for key, label in (("asserts", "asserts"), ("inconclusive", "inconclusive"), ("fails", "fails"),
                       ("args", "arguments")):
        if entry.get(key):
            out += latex.catline(label, inline(entry[key])) + "\n"
    routes = entry.get("routes")
    if isinstance(routes, dict):
        for route, text in routes.items():
            if len(route) > 15:
                out += "\\catblock{%s}{%s}\n" % (esc(route), inline(text))
            else:
                out += latex.catline(esc(route), inline(text)) + "\n"
    return out + "\\par\\medskip\n"


def check_catalogue(ctx):
    run = ctx.run
    sids = ctx.selected_scenarios()
    kinds = {}
    for sid in sids:
        ch = ctx.units.by_scenario[sid]["chapter"]
        for e in run.scenario(sid).get("expect") or []:
            if isinstance(e, dict) and e.get("check"):
                kinds.setdefault(e["check"], set()).add(ch)
    catalogue = ctx.checks_prose.get("checks") if isinstance(ctx.checks_prose.get("checks"), dict) else {}
    prose_checks = {}
    for ch, p in ctx.prose.items():
        for name, text in (p.get("checks") or {}).items():
            if text:
                prose_checks.setdefault(name, text)
    out = ("Every check kind a scenario of the built chapters expects: what it asserts of an answer, when it is "
           "inconclusive and when it fails, its arguments and the SPEC sections it rests on (the catalogue prose of "
           "doc/prose/checks.json, else the check function's docstring in checks.py), whether it is a generic check "
           "(it holds on any answer of its route) or a property check (the scenario's own claim), and the chapters "
           "that use it.\n\n")
    common = ctx.checks_prose.get("common") or []
    if common:
        out += "\\grouphead{The rules before any check}\n" + paragraphs(common)
    generic = [k for k in sorted(kinds) if not run.is_property(k)]
    props = [k for k in sorted(kinds) if run.is_property(k)]
    if generic:
        out += "\\grouphead{Generic checks}\n"
    for name in generic + props:
        if props and name == props[0]:
            out += "\\grouphead{Property checks}\n"
        chapters_text = "chapter%s %s" % ("s" if len(kinds[name]) > 1 else "", esc(", ".join(sorted(
            kinds[name], key=lambda c: ctx.units.chapters.index(c) if c in ctx.units.chapters else 99))))
        entry = catalogue.get(name)
        if isinstance(entry, dict) and entry.get("asserts"):
            out += _catalogue_entry(ctx, name, entry, chapters_text, name in generic)
            continue
        doc = prose_checks.get(name)
        source = "prose"
        if not doc:
            doc = run.check_doc(name)
            source = "the docstring in checks.py"
        text = inline(doc) if (doc and source == "prose") else (esc_wrap(doc) if doc else latex.placeholder(
            "(no description: neither a catalogue entry nor a docstring)"))
        out += "\\Needspace{4\\baselineskip}\\noindent\\textbf{\\texttt{%s}}\\hfill{\\small\\itshape %s check · %s · %s}\\par\\nopagebreak\n{\\small %s}\\par\\medskip\n" % (
            esc(name), "generic" if name in generic else "property", chapters_text, source, text)
    return out


def by_design_table(ctx):
    run = ctx.run
    sids = ctx.selected_scenarios()
    wanted = set()
    naming = {}
    for sid in sids:
        bd = (run.scenario(sid).get("mcp") or {}).get("by_design")
        for case in ([bd] if isinstance(bd, str) else (bd or [])):
            wanted.add(case)
            naming.setdefault(case, []).append(sid)
        for l in run.arm_lines(sid, "mcp"):
            if l.get("service_by_design_case"):
                wanted.add(l["service_by_design_case"])
                naming.setdefault(l["service_by_design_case"], [])
                if sid not in naming[l["service_by_design_case"]]:
                    naming[l["service_by_design_case"]].append(sid)
    for cid, c in run.table.items():
        if set(c.get("areas") or []) & set(ctx.chapters):
            wanted.add(cid)
    rows = []
    for cid in sorted(wanted):
        c = run.table.get(cid) or run.cases.get(cid) or {}
        sa = c.get("service_answer") or {}
        de = c.get("direct_expectation") or {}
        service = []
        if sa.get("outcome"):
            service.append(esc(sa["outcome"]))
        if sa.get("stage"):
            service.append("stage %s" % esc(sa["stage"]))
        if sa.get("http") is not None:
            service.append("http %s" % esc(sa["http"]))
        if sa.get("code"):
            service.append(tt(sa["code"]))
        if sa.get("text_contains"):
            service.append("“%s”" % esc_wrap(sa["text_contains"], 100))
        if sa.get("missing_argument"):
            service.append("missing argument %s" % tt(sa["missing_argument"]))
        if c.get("explains_checks"):
            service.append("explains %s" % esc(", ".join(c["explains_checks"])))
        direct = []
        if de.get("http") is not None:
            direct.append("http %s" % esc(de["http"]))
        if de.get("text"):
            direct.append(tt(de["text"], 80))
        elif de.get("text_prefix"):
            direct.append("prefix %s" % tt(de["text_prefix"], 80))
        marks = {}
        for sid in naming.get(cid, []):
            ma = run.comp_entry(sid).get("mcp_agreement") or {}
            if ma.get("case") == cid and ma.get("class") == "service_by_design":
                mark = "verified" if ma.get("verified") else "unverified"
            elif ma.get("case") == cid:
                mark = (ma.get("class") or "").replace("_", " ")
            else:
                mark = "named"
            marks.setdefault(mark, []).append(sid)
        named = []
        for mark in sorted(marks):
            sids_ = marks[mark]
            if len(sids_) <= 6:
                named.append("%s: %s" % (esc(mark), ", ".join(latex.sid(s, ctx.units.unit_of(s)) for s in sids_)))
            else:
                named.append("%s: %d (%s, …)" % (esc(mark), len(sids_),
                                                 ", ".join(latex.sid(s, ctx.units.unit_of(s)) for s in sids_[:3])))
        rows.append([tt(cid) + "\\newline{\\scriptsize %s}" % esc(c.get("agreement") or "–"),
                     esc_wrap(c.get("scenario_match"), 220), "; ".join(service) or "–", "; ".join(direct) or "–",
                     "; ".join(named) or "–"])
    out = ("One row per case of the service-by-design table that a scenario of the built chapters names, that an "
           "MCP answer matched, or that belongs to their areas: the case and its agreement class, what it matches, "
           "the service's answer, the direct expectation, and the scenarios naming it with their verification.\n\n")
    out += latex.plain_table(latex.spec_p([34, 40, 44, 26, 30]),
                             ["case", "matches", "service", "direct", "scenarios"], rows)
    return out


def relations_appendix(ctx):
    run = ctx.run
    comp_rel = run.comp_relations()
    selected = set(ctx.selected_scenarios())
    rows = []
    for rel in run.relations:
        if rel.get("a") not in selected and rel.get("b") not in selected:
            continue
        crel = comp_rel.get(rel.get("id")) or {}
        args = rel.get("args") or {}
        arg_text = ", ".join("%s %s" % (k, v) for k, v in sorted(args.items()) if not isinstance(v, (list, dict)))
        for k, v in sorted(args.items()):
            if isinstance(v, (list, dict)):
                arg_text += ("; " if arg_text else "") + "%s %s" % (k, blocks._runs(v)[:40])
        per_arm = []
        for arm in ARMS:
            a = (crel.get("arms") or {}).get(arm) or {}
            if not a:
                per_arm.append("\\gnone")
            else:
                per_arm.append("%s {\\scriptsize %s}" % (cells.GLYPH.get(a.get("outcome"), "\\ginc"),
                                                        esc(a.get("detail"), 90)))
        rows.append([latex.breakable_id(rel.get("id")) + "\\newline{\\scriptsize %s}" % esc(rel.get("type")),
                     "%s %s %s" % (latex.sid(rel.get("a"), ctx.units.unit_of(rel.get("a"))),
                                   RELATION_SIGN.get(rel.get("type"), ""),
                                   latex.sid(rel.get("b"), ctx.units.unit_of(rel.get("b")))),
                     esc(rel.get("title"), 120) + (("\\newline{\\scriptsize %s}" % esc(arg_text, 100)) if arg_text else "")]
                    + per_arm)
    out = "Every manifest relation among the built chapters' tests, with the comparator's outcome and detail per arm.\n\n"
    out += latex.plain_table(latex.spec_p([28, 34, 40, 24, 24, 24]),
                             ["relation", "pair", "title", "direct-local", "direct-remote", "mcp"], rows)
    return out


def _sort_key(term):
    t = str(term).lower()
    for prefix in ("the ", "a "):
        if t.startswith(prefix):
            t = t[len(prefix):]
    return t


def glossary(ctx):
    """Appendix D: the entries of glossary.json ({term, definition | text, see}) in alphabetical order (a leading
    article ignored), two columns."""
    entries = [e for e in ((ctx.glossary or {}).get("entries") or []) if isinstance(e, dict) and e.get("term")]
    if not entries:
        return latex.placeholder("(glossary: to be written)") + "\n\n"
    out = "\\begin{multicols}{2}\\small\\raggedright\n"
    for e in sorted(entries, key=lambda e: _sort_key(e["term"])):
        text = inline(ctx.fill(e.get("definition") or e.get("text") or ""))
        see = e.get("see")
        see_text = inline(", ".join(str(x) for x in see)) if isinstance(see, list) and see else None
        out += latex.glossary_entry(esc(e["term"]), text, see_text)
    return out + "\\end{multicols}\n\n"


def index(ctx):
    sids = ctx.selected_scenarios()
    out = ("Every test id of the built chapters, its unit and the page of the unit's table; a footnoted id also "
           "gives the page of its detail block.\n\n")
    out += "\\begin{multicols}{3}\\fontsize{7}{8.4}\\selectfont\\ttfamily\\raggedright\n"
    for sid in sids:
        uid = ctx.units.unit_of(sid)
        entry = "%s\\ {\\color{black!60}%s}\\ \\pageref{unit:%s}" % (latex.sid(sid, uid), esc(uid),
                                                                     latex.label_key(uid))
        if ctx.footnote(sid) is not None:
            entry += "\\ d\\,\\pageref{detail:%s}" % latex.label_key(sid)
        out += "\\noindent %s\\\\\n" % entry
    out += "\\end{multicols}\n"
    return out


def appendices(ctx):
    out = "\\part*{Appendices}\\addcontentsline{toc}{part}{Appendices}\n\\markright{Appendices}\n"
    out += latex.section("A\\qd The check catalogue", "app:A")
    out += check_catalogue(ctx)
    out += latex.section("B\\qd The by-design table", "app:B")
    out += by_design_table(ctx)
    out += latex.section("C\\qd The relations", "app:C")
    out += relations_appendix(ctx)
    out += latex.section("D\\qd Glossary", "app:D")
    out += glossary(ctx)
    out += latex.section("E\\qd Index of test ids", "app:E")
    out += index(ctx)
    return out
