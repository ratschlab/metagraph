"""The sections of the results document (DESIGN-DOC §1, §4, §6–§8, the owner's review of the sample): the cover
lines, Part I in its four sections (decision summary, methods, reading guide, provenance and coverage), the
chapter pages (the conclusion first), the family sections with their unit tables and detail blocks, the closing
section (observations, plots, tables) and the appendices. Every function returns LaTeX; the data come from a
`reader.Run`, a `units.Units` and the prose dicts.
"""

import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import blocks  # noqa: E402
import cells  # noqa: E402
import compare as cmp  # noqa: E402
import facts  # noqa: E402
import gate as gate_mod  # noqa: E402
import latex  # noqa: E402
import reader  # noqa: E402
import report  # noqa: E402
from cells import ARMS, fmt_t, fmt_delta, fmt_pair, fmt_n  # noqa: E402
from latex import esc, esc_wrap, tt, inline, paragraphs  # noqa: E402

RELATION_SIGN = {"prefix_equals": "{\\mdseries\\upshape\\glyphfont ⊑}", "subset": "{\\mdseries\\upshape\\glyphfont ⊆}",
                 "same_answer": "=", "zero_equals_absent": "{\\mdseries\\upshape\\glyphfont ≡}"}
LOAD_CAVEAT_SRA = report.LOAD_CAVEAT_SRA_CPU_PCT
LOAD_CAVEAT_LOAD1 = report.LOAD_CAVEAT_LOAD1
# Part I's four sections (the owner's review, point 2): key, the default title, the id the prose may use instead
PART_ONE = [("decision", "Decision summary", "I.1"), ("methods", "Methods", "I.2"),
            ("reading", "Reading guide", "I.3"), ("provenance", "Provenance and coverage", "I.4")]
# the appendices: letter, title, label
APPENDICES = [("A", "The check catalogue"), ("B", "The by-design table"), ("C", "The relations"),
              ("D", "Record key and glossary"), ("E", "Index of test ids"),
              ("F", "Environment, vocabulary and the gate in full"),
              ("G", "Known MCP limits and owner-declared exemptions")]
README_PATH = "metagraph/scripts/suite/README.md"
README_SECTION = "Glossary"
PLOT_FILES = {"paired_dl_mcp": "plot_paired_dl_mcp.csv", "paired_dl_dr": "plot_paired_dl_dr.csv",
              "phases": "plot_phases.csv"}


class Context(object):
    """What every section reads: the run, the unit model, the prose, the chapters to build, the footnotes assigned
    (per chapter) and the per-scenario arm states (computed once)."""

    def __init__(self, run, units, prose, overview, glossary, chapters, checks_prose=None, closing=None,
                 cases=None, sample=False):
        self.run = run
        self.units = units
        self.prose = prose  # {chapter: prose dict}
        self.overview = overview or {}
        self.glossary = glossary or {}
        self.checks_prose = checks_prose or {}  # doc/prose/checks.json: {common: [...], checks: {kind: {...}}}
        self.closing = closing or {}  # doc/prose/closing.json
        self.cases_prose = cases or {}  # doc/prose/cases.json: {short: {case id: text}}
        self.sample = sample
        self.chapters = chapters
        self.states = {}
        self.footnotes = {}  # chapter -> Footnotes
        self.blocks_by_chapter = {}
        self.data_files = {}  # file name -> text (the plots' data, written beside the .tex)
        self.unresolved = []  # the prose placeholders and conditions the values lacked (facts.fill)
        self._closing_values = None
        self._untested = None
        self._facts = {}  # chapter -> facts.chapter_facts, "fam:" + family -> facts.family_facts

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

    def short_reasons(self):
        s = self.cases_prose.get("short")
        return s if isinstance(s, dict) else {}

    def placeholders(self):
        """The values of the prose placeholders (overview.json `placeholders` names them with their sources): plain
        text, never LaTeX, since the prose they fill is escaped afterwards."""
        run = self.run
        manifest = run.manifest
        index = manifest.get("index") or {}
        counts = manifest.get("counts") or {}
        arms = run.run.get("arms") or {}
        hc = run.comp.get("host_cache") or {}
        summary = run.summary()
        totals = summary.get("totals") or {}
        ma = summary.get("mcp_agreement") or {}
        gate = run.gate if isinstance(run.gate, dict) else {}
        gate_word, gate_how, gate_counts, gate_problems_text = gate_summary(self)
        gp = gate_outcome_counts(run)
        limits = mcp_limit_counts(run)
        declared = declared_counts(run)
        table = run.table_meta or {}
        return {"release": run.run.get("release") or "–", "run_id": run.run_id(),
                "database": index.get("database") or "–", "k": index.get("k") or "–",
                "scenarios": len(run.order), "n_scenarios": len(run.order),
                "manifest_version": manifest.get("manifest_version") or "–",
                "manifest_sha_short": (manifest.get("content_sha256") or "–")[:16],
                "n_records": len(run.records()),
                "n_relations": counts.get("relations", len(run.relations)),
                "server_instance": run.run.get("server_instance") or "–",
                "index_fp": (run.run.get("index_fp") or "–")[:16],
                "checks_commit": checks_commit_text(run),
                "scope_note": scope_note_text(self),
                "local_base": (arms.get("direct-local") or {}).get("base") or "–",
                "remote_base": (arms.get("direct-remote") or {}).get("base") or "–",
                "mcp_base": (arms.get("mcp") or {}).get("base") or run.mcp_run.get("base") or "–",
                "warm_window_s": fmt_n(hc.get("warm_window_s", 60), 0),
                "load_caveat_sra_cpu_pct": fmt_n(LOAD_CAVEAT_SRA, 0),
                "load_caveat_load1": fmt_n(LOAD_CAVEAT_LOAD1, 0),
                "load_caveat_line": load_caveat_line(run),
                "verdict_word": (run.verdict or "–").split(" ")[0],
                "verdict_text": run.verdict or "–",
                "n_findings": sum((totals.get(a) or {}).get("outcomes", {}).get("finding", 0) for a in ARMS),
                "n_inconclusive": sum((totals.get(a) or {}).get("outcomes", {}).get("inconclusive", 0) for a in ARMS),
                "n_failures": sum((totals.get(a) or {}).get("outcomes", {}).get(k, 0) for a in ARMS
                                  for k in ("fail", "error")),
                "n_executions": "{:,}".format(sum((totals.get(a) or {}).get("executions", 0) for a in ARMS)),
                "inconclusive_by_reason": inconclusive_by_reason_text(run),
                "n_not_implemented_mcp": (totals.get("mcp") or {}).get("outcomes", {}).get("not_implemented_mcp", 0),
                "n_service_by_design": (totals.get("mcp") or {}).get("outcomes", {}).get("service_by_design", 0),
                "n_by_design_verified": ma.get("by_design_verified", 0),
                "n_by_design_unverified": ma.get("by_design_unverified", 0),
                "n_not_expressible_scenarios": limits["not_expressible_scenarios"],
                "n_not_expressible_cases": len(limits["named"]),
                "n_manifest_not_implemented": limits["manifest_not_implemented"],
                "declared_file": declared["file"], "declared_version": declared["version"],
                "n_declared_exclusions": declared["exclusions"],
                "n_declared_not_required": declared["not_required"],
                "n_declared_relations": declared["relations_excluded"],
                "gate_word": gate_word if gate else "not evaluated",
                "gate_text": gate_text(gate_word, gate_how, gate_counts, gate_problems_text),
                "n_gate_outcome_problems": gp["outcome"], "n_gate_inconclusive_executions": gp["inconclusive"],
                "n_gate_finding_executions": gp["finding"], "n_gate_outcome_pairs": gp["pairs"],
                "gate_inconclusive_clause": gate_inconclusive_clause(run),
                "untested_scope": untested_scope_text(self),
                "by_design_cases": table.get("cases") if table.get("cases") is not None else "–",
                "by_design_sha_short": (table.get("sha256") or "–")[:12],
                "by_design_table_line": by_design_table_text(run),
                "n_chapters": len(self.units.chapters), "built_chapters": ", ".join(self.chapters),
                "readme_path": README_PATH, "readme_commit": self.readme_commit() or "not determined",
                "logan_caveat": facts.logan_caveat(run), "gate_judged": bool(gate)}

    def fill(self, text, extra=None):
        """Prose placeholders `{release}` etc. filled from the run (and `extra`, a chapter's, a family's or the
        closing section's own values), the conditionals `[[flag ? a | b]]` resolved (facts.fill); what the values
        lack stays visible and is recorded in `unresolved`."""
        values = self.placeholders()
        if extra:
            values.update(extra)
        return facts.fill(text, values, self.unresolved)

    def readme_commit(self):
        ref = self.glossary.get("readme") if isinstance(self.glossary.get("readme"), dict) else {}
        if ref.get("commit"):
            return str(ref["commit"])
        suite_dir = os.path.dirname(HERE)
        return reader.git_commit_of(os.path.join(suite_dir, "README.md"), suite_dir)

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

    def prose_paragraphs(self, section_key, sub=None, optional=False):
        """The paragraphs of a Part I section (`paragraphs`, or the slot `sub`), filled; a missing slot prints a
        visible placeholder, or nothing when the slot is `optional`."""
        sec = self.overview_section(section_key)
        items = (sec.get(sub) if sub else sec.get("paragraphs")) or []
        if isinstance(items, str):
            items = [items]
        if not items:
            if optional:
                return ""
            return latex.placeholder("(overview prose for %s%s: to be written)" % (
                section_key, ("." + sub) if sub else "")) + "\n\n"
        return paragraphs([self.fill(p) for p in items])

    def chapter_prose(self, ch):
        return self.prose.get(ch) or {}

    def chapter_values(self, ch):
        """The chapter's prose values (facts.chapter_facts), computed once per chapter after its families are
        built (the footnotes are assigned then)."""
        if ch not in self._facts:
            self._facts[ch] = facts.chapter_facts(self, ch)
        return self._facts[ch]

    def untested(self):
        if self._untested is None:
            self._untested = untested_properties(self)
        return self._untested

    def family_values(self, fid):
        """A family's own placeholders for its prose (facts.family_facts): {n_tests, n_units}, the MCP arm's
        outcome classes over its members as the tables derive them (blocks.mcp_class; the owner's review, point
        5), so that the prose's counts cannot drift from the cells (VERIFY-R2 M1): `n_mcp_<class>` per class
        (not_expressible, refused_before_host, relayed, normalised, answered, host_state, refused_unexplained,
        none), `mcp_<class>_ids` and `mcp_<class>_variants`; and the run facts of round 4 (the standing counts,
        the footnoted members' phrases, the findings, the late 503s, the relations, the declared exemptions, the
        members' facts under `x`)."""
        key = "fam:" + fid
        if key not in self._facts:
            self._facts[key] = facts.family_facts(self, fid)
        return self._facts[key]


# -- small helpers --------------------------------------------------------------------------------------------

def _variant(meta, sid=None):
    """A scenario's variant label as printed (units.variant_text: thousands separators), else its id."""
    import units as units_mod
    return units_mod.variant_text((meta or {}).get("variant") or sid or "")


def _kv_table(rows, widths=(32, 142)):
    spec = latex.spec_p(widths)
    return latex.plain_table(spec, ["", "value"], [[k, v] for k, v in rows])


def _stat(d, key="p50"):
    return fmt_t((d or {}).get(key)) if isinstance(d, dict) else "–"


def _stats_cells(d, signed=False):
    """p50 / p90 / max / n of a summarise() dict; `signed` for a difference (an overhead, a net cost, a delta),
    whose figures carry their sign."""
    if not isinstance(d, dict) or not d.get("n"):
        return ["–", "–", "–", "0"]
    f = fmt_delta if signed else fmt_t
    return [f(d.get("p50")), f(d.get("p90")), f(d.get("max")), str(d.get("n"))]


def _short(text, n=12):
    return esc(str(text)[:n]) if text else "–"


def _date(text):
    return esc(str(text)[:19].replace("T", " ")) if text else "–"


def _arm_end(run, arm, info):
    """When an arm ended: run.json's ended_utc (the direct runners), else its finished_utc (the MCP runner's
    name for it), else the MCP runner's own run.json."""
    info = info or {}
    end = info.get("ended_utc") or info.get("finished_utc")
    if not end and arm == "mcp":
        end = (run.mcp_run or {}).get("finished_utc") or (run.mcp_run or {}).get("ended_utc")
    return end


def _sessions(run):
    sessions = [x for x in ((run.run.get("orchestrator") or {}).get("sessions") or []) if isinstance(x, dict)]
    if not sessions:
        sessions = [x for x in (run.orchestrator.get("sessions") or []) if isinstance(x, dict)]
    return sessions


def _applied_offset(run):
    hc = run.comp.get("host_cache") or run.run.get("host_cache") or {}
    offsets = hc.get("clock_offsets_s") or {}
    if "direct-local" in offsets and cells._num(offsets["direct-local"]) is not None:
        return offsets["direct-local"]
    return None


def _sessions_text(run):
    """The orchestrator's sessions (run_suite.sh invocations: the first runs the loop, a later one resumes it),
    each with its start, its loop state and its direct requests, and the mex clock offset the comparison applied
    (compare.json host_cache.clock_offsets_s, the one that keyed the host-cache states)."""
    sessions = _sessions(run)
    if not sessions:
        return None
    parts = []
    for i, ses in enumerate(sessions):
        loop = ses.get("loop") or {}
        text = "session %d started %s" % (i + 1, _date(ses.get("started_utc")))
        if loop:
            text += ", loop %s" % esc(loop.get("state") or "–")
            if loop.get("reason"):
                text += " (%s)" % esc(loop["reason"], 120)
            if "completed" in loop or "scenarios" in loop:
                text += ", %s of %s scenarios" % (esc(loop.get("completed", "–")), esc(loop.get("scenarios", "–")))
            if loop.get("direct_requests") is not None:
                text += ", %s direct request%s" % (esc(loop["direct_requests"]),
                                                   "" if loop["direct_requests"] == 1 else "s")
        if ses.get("mex_clock_offset_s") is not None:
            text += " (mex clock offset measured %s s)" % esc(ses["mex_clock_offset_s"])
        parts.append(text)
    mode = esc(sessions[0].get("mode") or "–")
    applied = _applied_offset(run)
    if applied is not None:
        applied_text = "mex clock %+.3f s against the laptop, applied to direct-local by the comparison" % applied
    else:
        applied_text = "mex clock offset not applied by the comparison"
    return "%s; %s; %s" % (mode, "; ".join(parts), applied_text)


def _link_sid(ctx, sid):
    """A test id linked to its unit when its chapter is built, plain otherwise (the whole-run context)."""
    ch = (ctx.units.by_scenario.get(sid) or {}).get("chapter")
    if ch in ctx.chapters:
        return latex.sid(sid, ctx.units.unit_of(sid), detail=ctx.footnote(sid) is not None)
    return latex.breakable_id(sid)


# -- the plain-text facts the prose placeholders state ---------------------------------------------------------------

def checks_commit_text(run):
    """The commit of the checks that judged: the re-evaluation's, else the runners' scripts commit (recorded)."""
    if run.reevaluated:
        commit = ((run.reeval or {}).get("checks") or {}).get("commit")
        if commit:
            return str(commit)[:12] + (" (dirty)" if ((run.reeval or {}).get("checks") or {}).get("dirty") else "")
    arms = run.run.get("arms") or {}
    commit = next((str((info or {}).get("scripts_commit"))[:12] for info in arms.values()
                   if (info or {}).get("scripts_commit")), None)
    return ("%s (the runners' scripts, as recorded)" % commit) if commit else "not recorded"


def scope_note_text(ctx):
    """The scope sentence of a partial build (the cover's note, in prose), empty for a full build."""
    if ctx.sample:
        return "Chapter R sample; whole-run summary shown for context."
    if set(ctx.chapters) != set(ctx.units.chapters):
        return "Chapter%s %s only; whole-run summary shown for context." % (
            "s" if len(ctx.chapters) > 1 else "", ", ".join(ctx.chapters))
    return ""


def load_busy(run):
    lf = run.load or {}
    return (lf.get("sra_mean") or 0) > LOAD_CAVEAT_SRA or (lf.get("load_mean") or 0) > LOAD_CAVEAT_LOAD1


def load_caveat_line(run):
    """One sentence: whether the Logan caveat applies to the run, with the sampler's means."""
    lf = run.load or {}
    if not lf.get("n"):
        return "No load samples were recorded, so whether the Logan caveat applies is not known."
    if load_busy(run):
        return ("This run carries the Logan caveat: the sra_all CPU mean was %s %% and the load1 mean %s over %s "
                "samples, so every timing in it was taken under that load and is not the timing baseline." % (
                    fmt_n(lf.get("sra_mean")), fmt_n(lf.get("load_mean"), 1), "{:,}".format(lf.get("n", 0))))
    return "The Logan caveat does not apply to this run (sra_all CPU mean %s %%, load1 mean %s over %s samples)." % (
        fmt_n(lf.get("sra_mean")), fmt_n(lf.get("load_mean"), 1), "{:,}".format(lf.get("n", 0)))


def inconclusive_by_reason_text(run):
    """`direct-local 40 by the checks; mcp 43 by the checks, 3 by the service's traversal_resolve_deadline (a host
    state)`: compare's inconclusive_by_reason in words, per arm."""
    ibr = run.summary().get("inconclusive_by_reason") or {}
    parts = []
    for arm in ARMS:
        d = ibr.get(arm) or {}
        if not d:
            continue
        words = []
        for reason, n in sorted(d.items(), key=lambda kv: (-kv[1], kv[0])):
            if reason == "checks":
                words.append("%d by the checks" % n)
            else:
                words.append("%d by the service's %s (a host state)" % (n, reason))
        parts.append("%s %s" % (arm, ", ".join(words)))
    return "; ".join(parts) if parts else "none"


def mcp_limit_counts(run):
    """What the MCP arm could not ask: {manifest_not_implemented: scenarios with the manifest's mark, reasons:
    {reason: [sids]}, named: {case: [sids]} (the not-expressible table cases scenarios name),
    not_expressible_scenarios: how many scenarios name one, matched: the cases an MCP answer matched}."""
    reasons = {}
    named = {}
    matched = set()
    n_scen = 0
    for sid in run.order:
        mcp = run.scenario(sid).get("mcp") or {}
        if mcp.get("not_implemented"):
            reasons.setdefault(str(mcp["not_implemented"]), []).append(sid)
        bd = mcp.get("by_design")
        hit = False
        for case in ([bd] if isinstance(bd, str) else (bd or [])):
            c = run.cases.get(case) or run.table.get(case) or {}
            if c.get("agreement") == "not_implemented_mcp":
                named.setdefault(case, []).append(sid)
                hit = True
        if hit:
            n_scen += 1
        for l in run.arm_lines(sid, "mcp"):
            if l.get("service_by_design_case"):
                matched.add(l["service_by_design_case"])
    return {"manifest_not_implemented": sum(len(v) for v in reasons.values()), "reasons": reasons, "named": named,
            "not_expressible_scenarios": n_scen, "matched": matched}


def declared_counts(run):
    dec = run.declared
    if not isinstance(dec, dict):
        return {"file": "none", "version": "–", "exclusions": 0, "not_required": 0, "relations_excluded": 0}
    out = {"file": dec.get("_file") or "declared", "version": dec.get("declared_version") or "–"}
    for k in ("exclusions", "not_required", "relations_excluded"):
        out[k] = len([e for e in (dec.get(k) or []) if isinstance(e, dict)])
    return out


def gate_text(word, how, counts_text, problems_text):
    """The deploy gate in one plain line (the cover's and I.1's wording)."""
    if counts_text is None:
        return "%s (%s)" % (word, how)
    return "%s: %s; %s; %s" % (word, problems_text, counts_text, how)


def gate_outcome_counts(run):
    """The gate's outcome problems (gate.json, else the offline verdict): {judged, outcome: all of them,
    inconclusive: the inconclusive executions among them, finding: the findings, pairs: distinct (scenario, arm)}.
    These, not the run's inconclusive count, are what blocks the gate: an execution on a gated arm, in the
    policy's selection and not declared exempt (VERIFY-R2 H1). Zeros with judged False when no gate was judged."""
    gate = run.gate if isinstance(run.gate, dict) else None
    out = {"judged": gate is not None, "outcome": 0, "inconclusive": 0, "finding": 0, "pairs": 0}
    if not gate:
        return out
    pairs = set()
    for p in gate.get("problems") or []:
        if not isinstance(p, dict) or p.get("kind") != "outcome":
            continue
        out["outcome"] += 1
        detail = str(p.get("detail") or "")
        if detail.startswith("inconclusive"):
            out["inconclusive"] += 1
        elif detail.startswith("finding"):
            out["finding"] += 1
        pairs.add((p.get("scenario"), p.get("arm")))
    out["pairs"] = len(pairs)
    return out


def gate_inconclusive_clause(run):
    """Plain text for I.1: how many of the run's inconclusive executions block the gate, and why those."""
    gp = gate_outcome_counts(run)
    if not gp["judged"]:
        return "whether any of them blocks the deploy gate was not judged (no gate verdict for this run)"
    n = gp["inconclusive"]
    return ("%d of them block%s the deploy gate: %s on a gated arm, in the policy's selection and not declared "
            "exempt (the gate's outcome problems, appendix F); the others are on an ungated arm, outside the "
            "selection or exempt" % (n, "s" if n == 1 else "", "the one" if n == 1 else "those"))


def untested_scope_text(ctx):
    """`in this run` for a full build, `in the built chapter R` (chapters R, T) for a partial one (VERIFY-R2 L7)."""
    if set(ctx.chapters) == set(ctx.units.chapters):
        return "in this run"
    return "in the built chapter%s %s" % ("s" if len(ctx.chapters) > 1 else "", ", ".join(ctx.chapters))


def by_design_table_text(run):
    """The by-design table the document prints, in plain words: version, cases, service commit, sha; whether it
    is the table the re-evaluation replayed, else the comparator's warning."""
    t = run.table_meta or {}
    if not t.get("path"):
        return "no by-design table" + ((": " + "; ".join(t["warnings"])) if t.get("warnings") else "")
    text = "version %s, %s cases, service commit %s…, sha256 %s…" % (
        t.get("version"), t.get("cases"), str(t.get("service_commit") or "–")[:12], str(t.get("sha256") or "–")[:12])
    if t.get("replayed"):
        text += " (the table the re-evaluation replayed)"
    elif t.get("replayed_sha256"):
        text += " (WARNING: not the table the re-evaluation replayed, sha256 %s…, %s cases; %s)" % (
            str(t["replayed_sha256"])[:12], t.get("replayed_cases", "–"),
            "; ".join(t.get("warnings") or []) or "the replayed table was not found")
    elif t.get("warnings"):
        text += " (%s)" % "; ".join(t["warnings"])
    return text


# -- the cover ---------------------------------------------------------------------------------------------------

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
                                               _date(_arm_end(run, arm, info)))
                               for arm, info in sorted(arms.items())) or "–"),
            ("set", "the re-evaluated set (checks commit %s)" % tt(_short(((run.reeval or {}).get("checks") or {}).get(
                "commit"), 12)) if run.reevaluated else "the recorded set")]
    if ctx.sample or set(ctx.chapters) != set(ctx.units.chapters):
        rows.append(("chapters built", "%s of %s (Part II holds %s; the whole-run summary is shown for context)" % (
            esc(", ".join(ctx.chapters)), esc(", ".join(ctx.units.chapters)), esc(", ".join(ctx.chapters)))))
    return _kv_table(rows)


def verdict_line(ctx):
    v = ctx.run.verdict or "–"
    word, _, rest = v.partition(" ")
    return "\\textbf{\\large %s}\\ %s" % (esc(word), esc(rest, 400))


def gate_summary(ctx):
    """(word, how, counts text, problems text), plain text: the deploy gate's verdict from the recorded gate.json,
    else the verdict gate.py gave offline for this document, else why there is none."""
    run = ctx.run
    gate = run.gate
    if isinstance(gate, dict):
        pol = gate.get("policy") or {}
        counts = gate.get("counts") or {}
        word = "GATE PASS" if gate.get("passed") else "GATE FAIL"
        problems = [p for p in (gate.get("problems") or []) if isinstance(p, dict)]
        if run.gate_offline:
            how = ("judged offline by gate.py for this document (the run recorded no gate.json); policy %s" %
                   os.path.basename(pol.get("path") or "–"))
        else:
            how = "recorded in the run, judged %s" % str(gate.get("judged_utc") or "–")[:19].replace("T", " ")
        how += ", %s version %s (sha256 %s…), arms %s" % (
            pol.get("name") or "–", pol.get("policy_version") or "–", str(pol.get("sha256") or "–")[:12],
            ", ".join((run.policy or {}).get("arms") or []) or "–")
        counts_text = "%s executions judged, %s passed, %s planned, %s required" % (
            counts.get("judged", "–"), counts.get("passed", "–"), counts.get("planned", "–"),
            counts.get("required", "–"))
        by_kind = "; ".join("%d %s" % (sum(1 for p in problems if p.get("kind") == k), k)
                            for k in gate_mod.KINDS if any(p.get("kind") == k for p in problems))
        problems_text = "%d problem%s%s" % (len(problems), "" if len(problems) == 1 else "s",
                                             (" (%s)" % by_kind) if by_kind else "")
        return word, how, counts_text, problems_text
    how = "no gate.json in the run"
    if run.gate_error:
        how += ", and gate.py could not judge it offline: %s" % latex.cut_text(str(run.gate_error), 200)
    elif run.policy is None:
        how += ", and no policy in the suite's gate directory pins this manifest, so none was judged offline"
    return "not evaluated", how, None, None


def gate_line(ctx):
    """The deploy gate on the cover, as prominent as the run verdict (the owner's review, point 2)."""
    word, how, counts_text, problems_text = gate_summary(ctx)
    if counts_text is None:
        return "\\textbf{\\large Deploy gate: %s}\\ (%s)" % (esc(word), esc(how))
    return "\\textbf{\\large %s}: %s; %s; %s" % (esc(word), esc(problems_text), esc(counts_text), esc(how))


def sample_note(ctx):
    if ctx.sample:
        return "\\samplenote{Chapter R sample; whole-run summary shown for context}\n"
    if set(ctx.chapters) != set(ctx.units.chapters):
        return "\\samplenote{Chapter%s %s only; whole-run summary shown for context}\n" % (
            "s" if len(ctx.chapters) > 1 else "", esc(", ".join(ctx.chapters)))
    return ""


# -- Part I: decision summary -------------------------------------------------------------------------------------

def _finding_evidence(line):
    """The evidence clause of an execution's finding, without its bracketed aside (the detail block has it all)."""
    fc = [c for c in (line.get("checks") or []) if isinstance(c, dict) and c.get("outcome") == "finding"]
    text = blocks._first_clause(fc[0].get("detail"), 140) if fc else "a stated bound broken"
    return text.split(" (", 1)[0] if " (" in text else text


def principal_findings(ctx):
    """The findings, failures, errors and disagreements of the whole run, one line per scenario, linked where the
    chapter is built; empty lists when there are none."""
    run = ctx.run
    findings, failures, disagreements = {}, {}, []
    for sid in run.order:
        for arm in ARMS:
            for l in run.arm_lines(sid, arm):
                if l.get("outcome") == "finding":
                    findings.setdefault(sid, []).append("%s r%s: %s" % (esc(cells.ARM_SHORT.get(arm, arm)),
                                                                        esc(l.get("repeat")),
                                                                        latex.join_units(esc(_finding_evidence(l)))))
                elif l.get("outcome") in ("fail", "error"):
                    failed = [c for c in (l.get("checks") or []) if isinstance(c, dict) and c.get("outcome") == "fail"]
                    what = ("%s: %s" % (esc(failed[0].get("check")),
                                        latex.join_units(esc(failed[0].get("detail"), 120)))) if failed \
                        else esc(cells.WORD.get(l.get("outcome"), l.get("outcome")))
                    failures.setdefault(sid, []).append("%s r%s %s" % (esc(cells.ARM_SHORT.get(arm, arm)),
                                                                       esc(l.get("repeat")), what))
        e = run.comp_entry(sid)
        for key, who in (("direct_agreement", "direct"), ("mcp_agreement", "mcp")):
            a = e.get(key) or {}
            if a.get("class") == "disagree":
                disagreements.append("%s (%s: %s)" % (_link_sid(ctx, sid), who, esc(", ".join(
                    str(c) for c in (a.get("differing_checks") or [])) or "a check differs", 120)))
    comp_rel = run.comp_relations()
    failed_rel = []
    for rid, r in sorted(comp_rel.items()):
        arms_ = r.get("arms") or {}
        bad = [a for a in ARMS if (arms_.get(a) or {}).get("outcome") == "fail"]
        if bad:
            failed_rel.append("%s on %s" % (latex.breakable_id(rid), esc(", ".join(bad))))
    out = []
    if findings:
        items = ["%s (%s): %s" % (_link_sid(ctx, sid), esc(run.scenario(sid).get("title"), 44), "; ".join(v))
                 for sid, v in findings.items()]
        out.append(("findings (%d execution%s, %d test%s)" % (
            sum(len(v) for v in findings.values()), "" if sum(len(v) for v in findings.values()) == 1 else "s",
            len(findings), "" if len(findings) == 1 else "s"), items))
    if failures:
        items = ["%s (%s): %s" % (_link_sid(ctx, sid), esc(run.scenario(sid).get("title"), 70), "; ".join(v))
                 for sid, v in failures.items()]
        out.append(("failures and errors (%d test%s)" % (len(failures), "" if len(failures) == 1 else "s"), items))
    if disagreements:
        out.append(("disagreements between arms (%d)" % len(disagreements), disagreements))
    if failed_rel:
        out.append(("failed relations (%d)" % len(failed_rel), failed_rel))
    return out


def _arm_reason(lines, st):
    """Why an arm did not pass a property, in a few words, from its lines."""
    if not lines:
        return "no execution"
    if st.glyph_key == "not_implemented_mcp":
        return "not exercised via MCP"
    if st.glyph_key == "service_by_design":
        return "decided by the service (by design)"
    words = []
    for l in lines:
        http = l.get("http")
        code = (l.get("error") or {}).get("code") if isinstance(l.get("error"), dict) else None
        if l.get("outcome") == "pass":
            w = "pass"
        elif http == 503 and code == "deadline":
            w = "503 deadline (allowed)" if l.get("outcome") == "inconclusive" else "503 deadline, a finding"
        elif l.get("inconclusive_reason") and l.get("inconclusive_reason") != "checks":
            w = "the service's %s (host state)" % l["inconclusive_reason"]
        elif http is not None and http != 200:
            w = "http %s%s" % (http, (" " + code) if code else "")
        else:
            w = cells.WORD.get(l.get("outcome"), str(l.get("outcome")))
        if w not in words:
            words.append(w)
    return " / ".join(words)


def untested_properties(ctx):
    """The property checks that no execution of a unit's members passed on any arm (untested) or that only some
    members passed (partly), per unit of the built chapters, with the reasons per arm; and the run-wide count of
    (scenario, property check) pairs never passed, per chapter (the whole-run context)."""
    run, units = ctx.run, ctx.units
    per_unit = []
    for ch in ctx.chapters:
        for fid in units.chapter_families(ch):
            for uid in units.families[fid]["units"]:
                members = units.units[uid]["members"]
                kinds = []
                for sid in members:
                    for e in run.scenario(sid).get("expect") or []:
                        name = e.get("check") if isinstance(e, dict) else None
                        if name and run.is_property(name) and name not in kinds:
                            kinds.append(name)
                for name in kinds:
                    expecting = [sid for sid in members if any(isinstance(e, dict) and e.get("check") == name
                                                                for e in run.scenario(sid).get("expect") or [])]
                    passed = []
                    for sid in expecting:
                        if any(c.get("check") == name and c.get("outcome") == "pass"
                               for arm in ARMS for l in run.arm_lines(sid, arm)
                               for c in (l.get("checks") or []) if isinstance(c, dict)):
                            passed.append(sid)
                    if len(passed) == len(expecting):
                        continue
                    missing = [sid for sid in expecting if sid not in passed]
                    reasons = []
                    for arm in ARMS:
                        words = []
                        for sid in missing:
                            w = _arm_reason(run.arm_lines(sid, arm), ctx.state(sid)[arm])
                            if w not in words:
                                words.append(w)
                        reasons.append((arm, "; ".join(words)))
                    # the expectation in words when the missing members share it, else the check's name alone
                    texts = []
                    for sid in missing:
                        e = next((e for e in run.scenario(sid).get("expect") or []
                                  if isinstance(e, dict) and e.get("check") == name), {})
                        t = blocks.expect_text(name, e.get("args"))
                        if t not in texts:
                            texts.append(t)
                    per_unit.append({"chapter": ch, "unit": uid, "check": name,
                                     "expect": texts[0] if len(texts) == 1 else "\\texttt{%s}" % esc(name),
                                     "n": len(expecting), "passed": len(passed), "missing": missing,
                                     "reasons": reasons})
    by_chapter = {}
    for sid in run.order:
        ch = (units.by_scenario.get(sid) or {}).get("chapter")
        for e in run.scenario(sid).get("expect") or []:
            name = e.get("check") if isinstance(e, dict) else None
            if not name or not run.is_property(name):
                continue
            ok = any(c.get("check") == name and c.get("outcome") == "pass" for arm in ARMS
                     for l in run.arm_lines(sid, arm) for c in (l.get("checks") or []) if isinstance(c, dict))
            d = by_chapter.setdefault(ch, {"pairs": 0, "untested": 0, "tests": set()})
            d["pairs"] += 1
            if not ok:
                d["untested"] += 1
                d["tests"].add(sid)
    return {"per_unit": per_unit, "by_chapter": by_chapter}


def untested_lines(ctx, chapter=None):
    """The untested and partly tested properties as lines (LaTeX), of one chapter or of every built one."""
    out = []
    for u in ctx.untested()["per_unit"]:
        if chapter and u["chapter"] != chapter:
            continue
        head = "untested" if not u["passed"] else "partly tested (%d of %d)" % (u["passed"], u["n"])
        ids = ", ".join(latex.sid(s, u["unit"], detail=ctx.footnote(s) is not None) for s in u["missing"][:8])
        if len(u["missing"]) > 8:
            ids += " … (%d)" % len(u["missing"])
        reasons = "; ".join("%s: %s" % (esc(cells.ARM_SHORT.get(a, a)), esc(w)) for a, w in u["reasons"] if w)
        n_missing = len(u["missing"])
        out.append("\\textbf{%s} · \\hyperref[unit:%s]{\\texttt{%s}} · %s: no arm passed it on %s (%d test%s) — %s" % (
            head, latex.label_key(u["unit"]), esc(u["unit"]), u["expect"], ids, n_missing,
            "" if n_missing == 1 else "s", reasons))
    return out


def untested_context_line(ctx):
    """One line: the run-wide count of property claims never passed, per chapter."""
    bc = ctx.untested()["by_chapter"]
    parts = []
    total_pairs, total_untested, total_tests = 0, 0, 0
    for ch in ctx.units.chapters:
        d = bc.get(ch)
        if not d:
            continue
        total_pairs += d["pairs"]
        total_untested += d["untested"]
        total_tests += len(d["tests"])
        if d["untested"]:
            parts.append("%s %d of %d (%d test%s)" % (esc(ch), d["untested"], d["pairs"], len(d["tests"]),
                                                      "" if len(d["tests"]) == 1 else "s"))
    return ("Over the whole run %d of %d property claims (a scenario's property check) were not passed on any "
            "arm, in %d test%s%s." % (total_untested, total_pairs, total_tests, "" if total_tests == 1 else "s",
                                       (": " + "; ".join(parts)) if parts else ""))


def verdict_table(ctx):
    totals = ctx.run.summary().get("totals") or {}
    rows = []
    for arm in ARMS:
        t = totals.get(arm) or {}
        o = t.get("outcomes") or {}
        hc = t.get("host_cache") or {}
        rows.append([esc(arm), str(t.get("executions", 0)), str(o.get("pass", 0)), str(o.get("fail", 0)),
                     str(o.get("error", 0)), str(o.get("finding", 0)), str(o.get("inconclusive", 0)),
                     str(o.get("not_implemented_mcp", 0)), str(o.get("service_by_design", 0)),
                     "%s / %s / %s / %s" % (hc.get("cold", 0), hc.get("warm", 0), hc.get("aged", 0),
                                            hc.get("not_reached", 0))])
    spec = "@{}l r r r r r r r r r@{}"
    return latex.plain_table(spec, ["arm", "executions", "\\gpass", "\\gfail", "\\gerr", "\\gfind", "\\ginc", "\\gni",
                                    "\\gbd", "cold / warm / aged / not reached"], rows,
                             "the recorded outcomes per arm (the report's counts; \\gni\\ and \\gbd\\ are outcomes "
                             "here, not the tables' standing)")


def decision_summary(ctx):
    """I.1: the scope (prose), the two verdict lines and how they differ (prose `verdict`), the principal findings
    (prose `findings`, then the computed lists), the properties still untested (prose `untested`, then the computed
    list and the run-wide count) and where to look (prose `closing`); the per-arm outcome table is I.4's."""
    out = ctx.prose_paragraphs("decision")
    out += "\\begingroup\\sffamily\n\\bline{run verdict}{%s}\n\\bline{deploy gate}{%s}\n\\par\\endgroup\\medskip\n\n" % (
        verdict_line(ctx), gate_line(ctx))
    out += ctx.prose_paragraphs("decision", "verdict", optional=True)
    out += "\\parhead{Principal findings}\n"
    out += ctx.prose_paragraphs("decision", "findings", optional=True)
    groups = principal_findings(ctx)
    if not groups:
        out += "None: no finding, no failure, no error, no disagreement, no failed relation in this run.\n\n"
    for title, items in groups:
        out += "\\noindent\\textbf{%s}\\par\\nopagebreak\n\\begin{itemize}\\setlength\\itemsep{0pt}\n" % esc(title)
        for it in items[:20]:
            out += "\\item %s\n" % it
        if len(items) > 20:
            out += "\\item … %d more (Part II)\n" % (len(items) - 20)
        out += "\\end{itemize}\n"
    out += "\\parhead{Properties still untested}\n"
    if ctx.overview_section("decision").get("untested"):
        out += ctx.prose_paragraphs("decision", "untested")
    else:
        # without the prose slot the by-reason line is the generator's own
        ibr = inconclusive_by_reason_text(ctx.run)
        if ibr != "none":
            out += ("Inconclusive executions by reason: %s (the checks: a check could not test its property, "
                    "typically an allowed 503 deadline; a host state: the service's own limit ended the call).\n\n"
                    % esc(ibr))
    lines = untested_lines(ctx)
    if lines:
        out += "\\begin{itemize}\\setlength\\itemsep{0pt}\n" + "".join("\\item %s\n" % l for l in lines) + "\\end{itemize}\n"
    else:
        out += "Every property check of the built chapters was passed by at least one execution of each test.\n\n"
    out += untested_context_line(ctx) + "\n\n"
    out += ctx.prose_paragraphs("decision", "closing", optional=True)
    return out


# -- Part I: methods ----------------------------------------------------------------------------------------------

def methods_section(ctx):
    """I.2: the prose alone (the arms' bases are in it, their executions in I.4's table, their states and runners
    in appendix F)."""
    return ctx.prose_paragraphs("methods")


# -- Part I: reading guide ----------------------------------------------------------------------------------------

LEGEND_GLYPH = {"pass": "\\gpass", "fail": "\\gfail", "error": "\\gerr", "finding": "\\gfind",
                "inconclusive": "\\ginc", "not_implemented": "\\gni", "by_design": "\\gbd", "none": "\\gnone"}
LEGEND_WORD = {"pass": "pass", "fail": "fail", "error": "error", "finding": "finding", "inconclusive": "inconclusive",
               "not_implemented": "not exercised via MCP", "by_design": "service-side by design",
               "none": "no execution"}
LEGEND_PREFIXES = {"not_implemented": ("not exercised via MCP:", "not implemented via MCP:")}
MARK_EXAMPLE = {"superscript digit": "\\ginc\\fn{1}", "not-equal sign": "\\gpass\\gne", "t": "\\gpass\\tmark{t}",
                "s": "\\gpass\\tmark{s}", "n": "\\gpass\\tmark{n}", "u": "\\gbd\\tmark{u}", "c": "9.1 s\\tmark{c}",
                "a": "9.1 s\\tmark{a}", "i": "9.1 s\\tmark{i}", "m": "\\ginc\\tmark{m}"}
ANNOTATION_DEFAULT = {
    "column": "a column is one varied value of the unit's dimension; the record key in brackets names the source "
              "record (the key in appendix D)",
    "glyph": "the glyph is the arm's standing on the test: pass, fail, error, finding, inconclusive, not "
             "implemented via MCP, service-side by design (the legend below)",
    "repeats": "two numbers separated by a slash are the repeats in order (on mcp the cold call, then the warm "
               "call); the unit is written once",
    "server": "the bracketed line is the host's own server time per execution; a dash when the answer carries "
              "none (a refusal, a 503)",
    "backend": "be: the host call as the service timed it (mcp only); the gap to the wall is the service's own "
               "cost",
    "footnote": "a red superscript digit footnotes a detail block after the table; it sits on the offending cells "
                "only",
    "mark": "a small letter after a time is the host-cache state when it is not the usual one (c cold, a aged, i "
            "warm after an interrupted execution); after a glyph: m mixed repeats, t s n the direct classes, "
            "u an unverified by-design pair",
    "ids": "the ids line names the members in column order (the common prefix once); a unit id line and the index "
           "find every test",
    "summary": "the family's summary line counts its tests per arm by standing and gives the wall p50 / p90 per "
               "arm, host-cold and host-warm apart",
}


def legend_table(ctx=None):
    """The glyph legend and the marks of the reading guide; the wording of the glyphs and marks comes from the
    overview prose (`legend`, `marks` of the reading section) when it has them, the symbols always from the
    macros the tables use."""
    reading = ctx.overview_section("reading") if ctx else {}
    legend = [e for e in (reading.get("legend") or []) if isinstance(e, dict) and e.get("key") in LEGEND_GLYPH]
    marks = [e for e in (reading.get("marks") or []) if isinstance(e, dict) and e.get("mark") in MARK_EXAMPLE]
    if legend and marks:
        glyph_rows = []
        for e in legend:
            meaning = str(e.get("meaning") or "")
            for prefix in LEGEND_PREFIXES.get(e["key"], (LEGEND_WORD[e["key"]] + ":",)):
                if meaning.startswith(prefix):
                    meaning = meaning[len(prefix):].strip()
                    break
            glyph_rows.append([LEGEND_GLYPH[e["key"]], esc(LEGEND_WORD[e["key"]]), inline(meaning)])
        mark_rows = [[MARK_EXAMPLE[e["mark"]], esc(e.get("on") or ""), inline(e.get("meaning") or "")]
                     for e in marks]
        return _legend_tables(glyph_rows, mark_rows)
    return _legend_builtin()


def _legend_tables(glyph_rows, mark_rows):
    spec = latex.spec_p([8, 36, 128], aligns=["centering", "raggedright", "raggedright"])
    glyphs = latex.plain_table(spec, ["", "glyph", "meaning"], glyph_rows)
    spec2 = latex.spec_p([30, 34, 108])
    # a longtable: it may continue on the next page under its repeated header (no heading is orphaned)
    marks_t = latex.plain_table(spec2, ["mark", "on", "meaning"], mark_rows)
    return glyphs + marks_t


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
             ("9.1 s\\tmark{i}", "a time", "this execution ran host-warm after an interrupted one (a deadline cut "
                                           "without an answer); the aggregates leave it out"),
             ("\\ginc\\tmark{m}", "a glyph", "the repeats' outcomes differ (the worst is shown)")]
    return _legend_tables(glyph_rows, [[a, esc(b), esc(c)] for a, b, c in marks])


def grammar_lines(ctx):
    """The two rules of the cell grammar the reading guide's prose does not state: the replicates aggregate and
    the edge tables' MCP cell (the rest, the lines of a cell, the rounding and the slash, is in the prose and in
    the legend line under every table)."""
    out = ("\\begingroup\\sffamily\\footnotesize\n"
           "\\bline{aggregates}{a replicates cell is ok n/n · wall median / maximum over each member's first "
           "host-warm execution (else its only one) and [server median]; a wall not all host-warm carries the state's "
           "mark, alone when shared (7.0 / 7.8 s\\tmark{c}), else with its share (\\tmark{c 2/6}); members never "
           "exercised read not exercised: n, never ok 0/0; the coverage line under the table counts per arm}\n"
           "\\bline{edge tables}{the mcp cell of a refusals or list table shows the recorded status, the service's "
           "code and the class in words (refused before host execution, not expressible through MCP, answered after "
           "service normalization, host refusal relayed through the service), derived from the evidence, never from prose}"
           "\n\\par\\endgroup\n\n")
    return out


def _pick_annotated_unit(ctx):
    """The series unit of the built chapters the reading guide annotates: one with a footnoted member and a
    repeat if there is one, else the first series unit."""
    best, best_score = None, -1
    for ch in ctx.chapters:
        for fid in ctx.units.chapter_families(ch):
            for uid in ctx.units.families[fid]["units"]:
                u = ctx.units.units[uid]
                if u["unit_type"] != "series":
                    continue
                score = 0
                for sid in u["members"]:
                    st = ctx.state(sid)
                    if not cells.is_clean(st):
                        score += 2
                    if any(len(st[a].lines) > 1 for a in ("direct-local", "direct-remote")):
                        score += 1
                    if any(st[a].time_marks and any(st[a].time_marks) for a in ("direct-local", "direct-remote")):
                        score += 1
                if score > best_score:
                    best, best_score = uid, score
    return best


def annotated_matrix(ctx):
    """One real series matrix of the built chapters with numbered markers, and the list that explains each
    marker (the owner's review, point 2: one annotated matrix)."""
    uid = _pick_annotated_unit(ctx)
    if not uid:
        return ""
    u = ctx.units.units[uid]
    members = list(u["members"])
    meta = ctx.units.by_scenario
    reading = ctx.overview_section("reading")
    texts = dict(ANNOTATION_DEFAULT)
    ann = reading.get("annotations")
    if isinstance(ann, dict):
        texts.update(dict((k, str(v)) for k, v in ann.items() if v))
    used = []

    def mark(key):
        if key not in used:
            used.append(key)
        return "\\amark{%d}" % (used.index(key) + 1)

    columns = [esc(_variant(meta[sid], sid)) for sid in members]
    columns[0] = columns[0] + mark("column")
    rows = []
    repeat_marked = server_marked = backend_marked = footnote_marked = mark_marked = False
    for arm in ARMS:
        cs = []
        for i, sid in enumerate(members):
            st = ctx.state(sid)[arm]
            fn = ctx.footnote(sid)
            l1, l2, l3 = cells.cell_parts(st, fn if (fn is not None and not st.clean) else None)
            if arm == "direct-local" and i == 0:
                l1 = (l1 or "") + mark("glyph")
            if not repeat_marked and len(st.lines) > 1 and st.present and l1 and " / " in l1:
                l1 += mark("repeats")
                repeat_marked = True
            if not server_marked and l2 and arm == "direct-local" and i == len(members) - 1:
                l2 += mark("server")
                server_marked = True
            if not backend_marked and l3:
                l3 += mark("backend")
                backend_marked = True
            if not footnote_marked and fn is not None and not st.clean:
                l1 += mark("footnote")
                footnote_marked = True
            if not mark_marked and st.present and arm != "mcp" and any(st.time_marks):
                l1 += mark("mark")
                mark_marked = True
            cs.append(latex.cell_lines(l1, l2, l3))
        rows.append((esc(arm), cs))
    if not footnote_marked:
        used.append("footnote") if "footnote" not in used else None
    if not mark_marked and "mark" not in used:
        used.append("mark")
    ids = units_ids(ctx, members) + mark("ids")
    out = "\\noindent{\\small The matrix of unit \\hyperref[unit:%s]{\\texttt{%s}} (%s), annotated:}\\par\\nopagebreak\\smallskip\n" % (
        latex.label_key(uid), esc(uid), esc(u.get("unit_title") or ""))
    out += latex.series_table(columns, rows, ids, cells.TIMING_LEGEND["series"] + mark("summary"))
    out += "\\begingroup\\footnotesize\n"
    for i, key in enumerate(used):
        out += "\\annline{%d}{%s}\n" % (i + 1, inline(ctx.fill(texts.get(key, key))))
    out += "\\par\\endgroup\\medskip\n\n"
    return out


def record_key(ctx, full=False):
    """The record abbreviation key from the manifest's ground truth: key, organism, accession, taxid (compact,
    two records per row) or every field (appendix D)."""
    recs = ctx.run.records()
    if not recs:
        return ""
    if full:
        rows = []
        for r in recs:
            rows.append([tt(r.get("key")), esc(r.get("label")), tt(r.get("accession")), tt(r.get("taxid")),
                         esc(r.get("kingdom")), "{:,}".format(int(r["length"])) if isinstance(r.get("length"), (int, float))
                         else "–", esc(r.get("records_in_column", "–")), esc(r.get("transl_table", "–")),
                         esc("yes" if r.get("last_in_column") else "no")])
        # the organism column wide enough for every record label on one line: the key stays within a page
        return latex.plain_table(latex.spec_p([9, 58, 26, 12, 15], extra_cols=4, extra_spec="r r r l"),
                                 ["key", "organism (record label)", "accession", "taxid", "kingdom", "length bp",
                                  "records in column", "code", "last in column"], rows,
                                 "the %d records of the ground truth, as the manifest states them" % len(recs))
    # the compact key: two records per row, the organism cut so that a row stays one line (the full names are in
    # appendix D), so that the key fits under the reading guide's rules
    half = (len(recs) + 1) // 2
    rows = []
    for i in range(half):
        a = recs[i]
        b = recs[i + half] if i + half < len(recs) else None
        row = [tt(a.get("key")), esc(a.get("label"), 28), tt(a.get("accession")), tt(a.get("taxid"))]
        row += [tt(b.get("key")), esc(b.get("label"), 28), tt(b.get("accession")), tt(b.get("taxid"))] if b else \
            ["", "", "", ""]
        rows.append(row)
    spec = latex.spec_p([8, 42, 23, 11, 8, 42, 23, 11])
    return latex.plain_table(spec, ["key", "organism", "accession", "taxid", "key", "organism", "accession", "taxid"],
                             rows, "record key: the abbreviations the matrices use (the full names are in appendix D)",
                             needspace=min(19, len(rows) + 3))


def record_key_label(chapter=None):
    """The label of the compact record key: one per chapter (a full build prints the key on every chapter page)."""
    return "rk:%s" % chapter if chapter else "rk:compact"


def record_key_paragraph(ctx, chapter=None):
    """The compact record key, one paragraph (abbreviation and organism, in the manifest's order), printed on the
    chapter page before its first family, where the matrices need it, and labelled rk:<chapter> so that every
    matrix's legend links to it (the owner's review, point 3: the abbreviations are needed in the matrices and
    deserve easier access than the implementation terms; re-verification N-1; in I.3 it cost Part I a fifth page,
    against point 2). The full key (accession, taxid, kingdom, length) is appendix D's."""
    recs = ctx.run.records()
    if not recs:
        return ""
    items = " · ".join("%s %s" % (tt(r.get("key")), esc(r.get("label"))) for r in recs)
    return ("\\par\\smallskip\\noindent\\phantomsection\\label{%s}{\\sffamily\\footnotesize\\textbf{Record key}\\ "
            "(the abbreviations of the matrices' columns and ids; accession, taxid, kingdom and length in "
            "\\hyperref[rk:D]{appendix D}, p.\\,\\pageref{rk:D}): %s.}\\par\\smallskip\n" % (
                latex.label_key(record_key_label(chapter)), items))


def record_key_link(ctx, members, chapter=None):
    """The legend's pointer to the record key (the chapter's compact key and appendix D's full one) when a unit's
    variants or groups name a record abbreviation, else an empty string."""
    keys = set(str(r.get("key")) for r in ctx.run.records() if r.get("key"))
    if not keys:
        return ""
    words = set()
    for sid in members:
        meta = ctx.units.by_scenario.get(sid) or {}
        for field in ("variant", "group"):
            words.update(re.findall(r"[A-Za-z0-9]+", str(meta.get(field) or "")))
    if not (words & keys):
        return ""
    label = latex.label_key(record_key_label(chapter))
    return "record keys: \\hyperref[%s]{p.\\,\\pageref{%s}}, in full \\hyperref[rk:D]{appendix D}" % (label, label)


def reading_guide(ctx):
    """I.3: the prose, the annotated matrix with its marker texts, then the legend, the marks and the two grammar
    rules (longtables that break across a page without an orphaned heading). The compact record key is on each
    chapter page (record_key_paragraph), the full key appendix D's (the owner's review, point 3; point 2 holds
    Part I to about four pages: VERIFY-R2 M3, re-verification N-1)."""
    out = ctx.prose_paragraphs("reading")
    out += annotated_matrix(ctx)
    out += legend_table(ctx)
    out += grammar_lines(ctx)
    return out


# -- Part I: provenance and coverage ----------------------------------------------------------------------------------

def provenance_table(ctx):
    run = ctx.run
    r = run.run
    manifest = run.manifest
    reeval = run.reeval or {}
    checks = reeval.get("checks") or {}
    table = run.table_meta or {}
    index = manifest.get("index") or {}
    rows = [("run", "%s; release %s…" % (tt(run.run_id()), tt(_short(r.get("release"), 16)))),
            ("server", "instance %s; index %s, k = %s, fingerprint %s…" % (
                tt(r.get("server_instance") or "–"), tt(index.get("database") or "–"), esc(index.get("k") or "–"),
                tt(_short(r.get("index_fp"), 12)))),
            ("manifest", "version %s, %d scenarios, %d relations; content sha256 %s…" % (
                esc(manifest.get("manifest_version")), len(run.order), len(run.relations),
                tt(_short(manifest.get("content_sha256"), 12))))]
    if run.reevaluated:
        rows.append(("checks", "re-evaluated offline with commit %s%s, %s" % (
            tt(_short(checks.get("commit"), 12)), " (dirty)" if checks.get("dirty") else "",
            _date(reeval.get("reevaluated_utc")))))
    else:
        rows.append(("checks", "the recorded set, judged by the runners' checks"))
    if table.get("path"):
        text = "version %s, %s cases, service commit %s…, sha256 %s…" % (
            esc(table.get("version")), esc(table.get("cases")), tt(_short(table.get("service_commit"), 12)),
            tt(_short(table.get("sha256"), 12)))
        if table.get("replayed"):
            text += ", the re-evaluation's"
        elif table.get("replayed_sha256"):
            text += " (\\textbf{not} the table the re-evaluation replayed: appendix F)"
        rows.append(("by-design table", text))
    elif table.get("warnings"):
        rows.append(("by-design table", "none read: %s" % esc("; ".join(table["warnings"]), 300)))
    for arm, info in sorted((r.get("arms") or {}).items()):
        info = info or {}
        rows.append(("arm %s" % esc(arm), "%s executions, %s to %s, runner %s" % (
            esc(info.get("executions_recorded", "–")), _date(info.get("started_utc")), _date(_arm_end(run, arm, info)),
            esc(info.get("runner_host") or "–"))))
    sessions = _sessions(run)
    if sessions:
        first = sessions[0]
        loop = first.get("loop") or {}
        applied = _applied_offset(run)
        rows.append(("sessions", "%d (%s); the loop %s of %s scenarios%s; mex clock offset %s" % (
            len(sessions), esc(first.get("mode") or "–"), esc(loop.get("completed", "–")), esc(loop.get("scenarios", "–")),
            (", %s" % esc(loop.get("state"))) if loop.get("state") else "",
            ("%+.3f s applied to direct-local" % applied) if applied is not None else "not applied")))
    lf = run.load or {}
    caveat = ""
    if (lf.get("sra_mean") or 0) > LOAD_CAVEAT_SRA or (lf.get("load_mean") or 0) > LOAD_CAVEAT_LOAD1:
        caveat = ": the Logan caveat applies (timings observed under concurrent host load)"
    rows.append(("host load", "%d samples; sra\\_all CPU mean %s\\%%, load1 mean %s%s" % (
        lf.get("n", 0), fmt_n(lf.get("sra_mean")), fmt_n(lf.get("load_mean"), 1), caveat)))
    return _kv_table(rows)


def limitations_table(ctx):
    """Limitations and exemption counts with links to the full entries."""
    run = ctx.run
    limits = mcp_limit_counts(run)
    ma = run.summary().get("mcp_agreement") or {}
    rows = [["not implemented via MCP by the manifest", "%d scenarios, %d distinct reasons" % (
        limits["manifest_not_implemented"], len(limits["reasons"])), "\\hyperref[app:G]{G}"],
            ["not-expressible table cases scenarios name", "%d cases, %d scenarios" % (
                len(limits["named"]), limits["not_expressible_scenarios"]), "\\hyperref[app:G]{G}"],
            ["by-design cases an MCP answer matched", "%d cases; pairs verified %s, unverified %s" % (
                len(limits["matched"]), esc(ma.get("by_design_verified", "–")), esc(ma.get("by_design_unverified", "–"))),
             "\\hyperref[app:B]{B}"]]
    dec = declared_counts(run)
    if isinstance(run.declared, dict):
        rows.append(["owner-declared exemptions",
                     "%d exclusions, %d checks not required, %d relations excluded" % (
                         dec["exclusions"], dec["not_required"], dec["relations_excluded"]),
                     "\\hyperref[app:G]{G}"])
    else:
        rows.append(["owner-declared exemptions", "none given for this manifest", ""])
    gate = run.gate if isinstance(run.gate, dict) else None
    if gate:
        problems = [p for p in (gate.get("problems") or []) if isinstance(p, dict)]
        rows.append(["deploy-gate problems", "%d (%s)" % (len(problems), "; ".join(
            "%d %s" % (sum(1 for p in problems if p.get("kind") == k), esc(k)) for k in gate_mod.KINDS
            if any(p.get("kind") == k for p in problems)) or "none"), "\\hyperref[app:F]{F}"])
    ibr = run.summary().get("inconclusive_by_reason") or {}
    n_inc = sum(sum((ibr.get(a) or {}).values()) for a in ARMS)
    gp = gate_outcome_counts(run)
    if gp["judged"]:
        blocking = "%d block the gate (a gated arm, in the policy's selection, not declared exempt)" % gp["inconclusive"]
    else:
        blocking = "the gate was not judged"
    rows.append(["inconclusive executions", "%d; %s; each footnoted in its chapter's tables" % (n_inc, esc(blocking)),
                 "\\hyperref[sec:I1]{I.1}"])
    return latex.plain_table(latex.spec_p([62, 94, 10]), ["limitation", "count", "in full"], rows,
                             needspace=len(rows) + 3)


def provenance_section(ctx):
    """I.4: what answered and what judged (prose, the compact table), the coverage (prose `coverage`, the per-arm
    outcome table) and the limitations and exemptions (prose `exemptions`, the counts with their appendix links)."""
    out = ctx.prose_paragraphs("provenance")
    out += provenance_table(ctx)
    out += "\\parhead{Coverage}\n"
    out += ctx.prose_paragraphs("provenance", "coverage", optional=True)
    out += verdict_table(ctx)
    out += "\\parhead{Limitations and exemptions}\n"
    out += ctx.prose_paragraphs("provenance", "exemptions", optional=True)
    out += limitations_table(ctx)
    return out


def part_one(ctx):
    out = "\\part*{Part I\\qd Overview}\\addcontentsline{toc}{part}{Part I\\qd Overview}\n\\markright{Overview}\n"
    builders = {"decision": decision_summary, "methods": methods_section, "reading": reading_guide,
                "provenance": provenance_section}
    for key, default, sid in PART_ONE:
        out += latex.section("%s\\qd %s" % (sid, esc(ctx.section_title(key, default))), "sec:%s" % sid.replace(".", ""))
        out += builders[key](ctx)
    return out


# -- chapters ----------------------------------------------------------------------------------------------------

def _agreement_counts(ctx, sids):
    run = ctx.run
    direct, mcp = {}, {}
    full = 0
    fewer = []  # per fully checked pair with a repeat left out of the comparison, the repeats compared
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
            if d.get("repeats_excluded"):
                fewer.append(len(d.get("repeats_compared") or []))
    dtext = ", ".join("%d %s" % (direct[k], k.replace("_", " ")) for k in
                      ("identical", "same_checks", "timing", "not_comparable", "disagree") if direct.get(k)) or "–"
    mtext = "%d agree (fully checked %d%s), %d not implemented, %d by design (%d verified, %d unverified), %d not " \
            "comparable, %d disagree%s" % (
                mcp.get("agree", 0), full, fewer_repeats_text(fewer), mcp.get("not_implemented_mcp", 0),
                mcp.get("service_by_design", 0), verified, unverified, mcp.get("not_comparable", 0),
                mcp.get("disagree", 0), (", %d missing" % mcp["missing"]) if mcp.get("missing") else "")
    return esc(dtext), esc(mtext)


def fewer_repeats_text(fewer):
    """The qualifier of a 'fully checked' count for the pairs whose comparison left a repeat out (an MCP execution
    inconclusive by a host state; the owner's second review, point 2): `fewer` holds, per such pair, the number of
    repeats compared. ', 1 of them on one repeat only, the other repeat an inconclusive host-state outcome'; empty
    when every fully checked pair was compared on every repeat."""
    if not fewer:
        return ""
    if all(n == 1 for n in fewer):
        return ", %d of them on one repeat only, the other repeat an inconclusive host-state outcome" % len(fewer)
    return ", %d of them with a repeat left out as an inconclusive host-state outcome" % len(fewer)


def chapter_summary(ctx, chapter, sids, n_units, n_families):
    run = ctx.run
    states = [ctx.state(sid) for sid in sids]
    out = "\\begingroup\\sffamily\\footnotesize\n"
    out += "\\bline{chapter %s}{%d tests · %d units in %d families}\n" % (esc(chapter), len(sids), n_units, n_families)
    out += "\\bline{standing}{%s (by glyph, I.3)}\n" % " · ".join("%s %s" % (esc(arm), cells.glyph_counts(states, arm))
                                                                 for arm in ARMS)
    per_area = (run.summary().get("per_area") or {}).get(chapter) or {}
    if per_area:
        otext = []
        for arm in ARMS:
            o = per_area.get(arm) or {}
            items = ["%d %s" % (o[k], w) for k, w in (("pass", "pass"), ("fail", "fail"), ("error", "error"),
                                                     ("finding", "finding"), ("inconclusive", "inconclusive"),
                                                     ("mixed", "mixed"), ("not_implemented_mcp", "not implemented"),
                                                     ("service_by_design", "by design"), ("other", "other"))
                     if o.get(k)]
            otext.append("%s %s" % (esc(arm), esc(", ".join(items)) if items else "–"))
        out += "\\bline{outcomes}{%s (the recorded outcomes, as the report counts them)}\n" % " · ".join(otext)
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
                                                   for arm in ARMS), latex.join_units(cells.timing_text(run, sids)))
    if relations:
        comp_rel = run.comp_relations()
        parts = []
        for arm in ARMS:
            counts = {}
            for rid in relations:
                o = ((comp_rel.get(rid) or {}).get("arms") or {}).get(arm, {}).get("outcome")
                counts[o] = counts.get(o, 0) + 1
            parts.append("%s %s" % (esc(arm), " ".join("%d\\,%s{}" % (counts[k], cells.GLYPH.get(k, "\\ginc"))
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


def _cell_fn(st, footnote):
    """The footnote mark goes on each offending cell (I.3): none on a clean arm of a footnoted test."""
    return footnote if (footnote is not None and not st.clean) else None


def _mcp_edge_cell(ctx, sid, st, footnote, with_times):
    """The mcp cell of a refusals or list row: the glyph, the class of the outcome in words, the recorded status
    and the service's code (the owner's review, point 5); for an exercised execution with `with_times` the usual
    time cell. Returns (cell LaTeX, note under the row or None)."""
    run = ctx.run
    cls, evidence = blocks.mcp_class(st.lines, run.scenario(sid), st.agreement, run.cases)
    glyph = st.glyph_tex(_cell_fn(st, footnote))
    note = None
    if cls == "none":
        return "%s\\ no execution" % glyph, None
    if cls == "not_expressible":
        case = st.agreement.get("case")
        reason = blocks.short_reason(run, sid, case, ctx.short_reasons())
        return "%s\\ %s\\newline %s; see B" % (glyph, esc(blocks.MCP_CLASS_WORD[cls]), reason), None
    first = next((l for l in st.lines if l.get("http") is not None), None)
    code = (first.get("error") or {}).get("code") if first and isinstance(first.get("error"), dict) else None
    status = "http %s%s" % (esc(first.get("http")), (" " + tt(code)) if code else "") if first else "–"
    if cls in ("refused_before_host", "relayed", "host_state", "refused_unexplained"):
        text = "%s\\ %s\\newline %s" % (glyph, esc(blocks.MCP_CLASS_WORD[cls]), status)
        return text, None
    if cls == "normalised":
        case = next((l.get("service_by_design_case") for l in st.lines if l.get("service_by_design_case")), None) \
            or st.agreement.get("case")
        text = "%s\\ %s\\newline %s" % (glyph, esc(blocks.MCP_CLASS_WORD[cls]), status)
        if with_times and st.present and st.glyph_key not in ("not_implemented_mcp", "service_by_design"):
            text = _cell_tex(st, footnote, "list")
        if case:
            note = "{\\footnotesize mcp, %s: %s; not a service refusal (case %s, appendix B)}" % (
                latex.breakable_id(sid), blocks.short_reason(run, sid, case, ctx.short_reasons()), tt(case))
        return text, note
    # answered by the host through the service: the usual cell
    if with_times:
        return _cell_tex(st, footnote, "list"), None
    return "%s\\ %s\\newline %s" % (glyph, esc(blocks.MCP_CLASS_WORD[cls]), status), None


EDGE_LEGEND = ("mcp column: the recorded status and the service's code; the class in words follows from whether "
               "the host was reached, the error's source and the by-design case (appendix B): refused before host "
               "execution · not expressible through MCP · answered after service normalization · host refusal "
               "relayed through the service")


def _unit_blocks(ctx, chapter, uid, members, unit_note=False):
    out = ""
    for sid in members:
        n = ctx.footnotes[chapter].get(sid)
        if n is None:
            continue
        out += blocks.detail_block(ctx.run, sid, ctx.state(sid), n, ctx.units, ctx.short_reasons(), unit_note)
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


def _unit_prose(ctx, chapter, uid):
    """(property line, note) of a unit: the prose's string, or {property, note}."""
    prose = ctx.prose.get(chapter) or {}
    entry = (prose.get("units") or {}).get(uid)
    if isinstance(entry, dict):
        return entry.get("property"), entry.get("note")
    return entry, None


def mcp_unit_note(ctx, chapter, uid, members):
    """The one line under a unit's table whose MCP column is wholly not exercised (the owner's review, point 4:
    a short local reason, the explanation in appendix B): the chapter prose's `mcp_notes.units[uid]`, else a line
    built from the members' short reasons. Empty when some member was exercised via MCP."""
    run = ctx.run
    states = [ctx.state(sid)["mcp"] for sid in members]
    if not states or not all(st.present and st.glyph_key == "not_implemented_mcp" for st in states):
        return ""
    notes = ((ctx.prose.get(chapter) or {}).get("mcp_notes") or {}).get("units") or {}
    text = notes.get(uid)
    if text:
        body = inline(ctx.fill(text, facts.unit_facts(ctx, uid)))
    else:
        reasons = []
        for sid in members:
            st = ctx.state(sid)["mcp"]
            case = st.agreement.get("case") or next((l.get("service_by_design_case") for l in st.lines
                                                     if l.get("service_by_design_case")), None)
            r = blocks.short_reason(run, sid, case, ctx.short_reasons())
            if r not in reasons:
                reasons.append(r)
        body = "Not exercised via MCP: %s; see B." % "; ".join(reasons)
    return "\\noindent{\\sffamily\\footnotesize\\gni\\ %s}\\par\\smallskip\n" % body


def unit_section(ctx, chapter, uid):
    run, units = ctx.run, ctx.units
    u = units.units[uid]
    members = list(u["members"])
    prop, note = _unit_prose(ctx, chapter, uid)
    prop_tex = inline(ctx.fill(prop)) if prop else latex.placeholder("(unit property: to be written)")
    if note:
        prop_tex += "\\newline{\\color{black!60}\\itshape %s}" % inline(ctx.fill(note))
    out = latex.unit_header(uid, esc(u.get("unit_title") or ""), u["unit_type"], prop_tex)
    layout = u["unit_type"]
    fns = dict((sid, _footnotes_for(ctx, chapter, sid)) for sid in members)
    meta = units.by_scenario
    legend = cells.TIMING_LEGEND.get(layout)
    rk = record_key_link(ctx, members, chapter)
    if rk:
        legend = (legend + "; " + rk) if legend else rk
    if layout == "series":
        columns = [esc(_variant(meta[sid], sid)) for sid in members]
        rows = []
        for arm in ARMS:
            rows.append((esc(arm), [_cell_tex(ctx.state(sid)[arm], _cell_fn(ctx.state(sid)[arm], fns[sid]), "series")
                                    for sid in members]))
        out += latex.series_table(columns, rows, units_ids(ctx, members), legend)
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
            variants = [_variant(meta[sid], sid) for sid in members]
            by_group[esc(", ".join(variants), 90)] = members
        for g, gm in by_group.items():
            st = dict((sid, ctx.state(sid)) for sid in gm)
            cs = []
            for arm in ARMS:
                l1, l2 = cells.aggregate(st, arm)
                cs.append(latex.cell_lines(l1, l2))
            rows.append((esc(g) if groups else g, len(gm), cs[0], cs[1], cs[2]))
        ids = units_ids(ctx, members, by_group if groups else None)
        out += latex.replicates_table("variants" if not groups else esc(u.get("unit_title") or uid)[:30], rows, ids,
                                      legend)
        # the coverage per arm in words (exercised, passed, not exercised), then the exceptions if any: an
        # "exceptions: none" alone must not read as coverage (the owner's review, point 4)
        all_states = dict((sid, ctx.state(sid)) for sid in members)
        coverage = " · ".join("%s %s" % (esc(arm), cells.coverage_text(all_states, arm)) for arm in ARMS)
        exceptions = []
        for sid in members:
            if fns[sid] is not None:
                st = ctx.state(sid)
                arms = [a for a in ARMS if not st[a].clean]
                worst = cells.worst_of([st[a].worst for a in arms if st[a].present]) if arms else None
                exceptions.append("%s %s\\fn{%d} (%s)" % (latex.sid(sid, detail=True),
                                                          cells.GLYPH.get(worst, "\\gnone"), fns[sid],
                                                          esc(", ".join(arms))))
        out += "\\noindent{\\sffamily\\footnotesize coverage: %s%s}\\par\\smallskip\n" % (
            coverage, ("; exceptions: " + "; ".join(exceptions)) if exceptions else "")
    elif layout == "refusals":
        rows = []
        prev = None
        for sid in members:
            st = ctx.state(sid)
            expected, prev = blocks.refusal_expected(run.scenario(sid), prev)
            cs = [st[arm].glyph_tex(_cell_fn(st[arm], fns[sid])) for arm in ("direct-local", "direct-remote")]
            mcp_cell, note = _mcp_edge_cell(ctx, sid, st["mcp"], fns[sid], with_times=False)
            rows.append(([latex.sid(sid, detail=True) if fns[sid] else latex.breakable_id(sid),
                          esc_wrap(_variant(meta[sid])), expected] + cs + [mcp_cell], note))
        out += latex.refusals_table(rows, units_ids(ctx, members), EDGE_LEGEND + (("; " + rk) if rk else ""))
        states = [ctx.state(sid) for sid in members]
        out += "\\sumline{%s: %d test%s · %s · %s}\n" % (tt(uid), len(members), "" if len(members) == 1 else "s", " · ".join(
            "%s %s" % (esc(arm), cells.glyph_counts(states, arm)) for arm in ARMS),
            latex.join_units(cells.max_wall_text(run, members)))
    elif layout == "list":
        rows = []
        for sid in members:
            st = ctx.state(sid)
            cs = [_cell_tex(st[arm], _cell_fn(st[arm], fns[sid]), "list") for arm in ("direct-local", "direct-remote")]
            mcp_cell, note = _mcp_edge_cell(ctx, sid, st["mcp"], fns[sid], with_times=True)
            rows.append(([latex.sid(sid, detail=True) if fns[sid] else latex.breakable_id(sid),
                          esc_wrap(_variant(meta[sid])),
                          blocks.expected_text(run.scenario(sid), run.is_property)] + cs + [mcp_cell], note))
        out += latex.list_table(rows, units_ids(ctx, members), cells.TIMING_LEGEND["list"] + "; " + EDGE_LEGEND
                                + (("; " + rk) if rk else ""))
    else:  # pairs
        rows = []
        for sid in members:
            st = ctx.state(sid)
            cs = []
            for arm in ARMS:
                if arm == "mcp" and st[arm].present and st[arm].glyph_key in ("not_implemented_mcp", "service_by_design"):
                    text, _ = _mcp_edge_cell(ctx, sid, st[arm], fns[sid], with_times=True)
                else:
                    text = _cell_tex(st[arm], _cell_fn(st[arm], fns[sid]), "pairs")
                    rep = run.scenario(sid).get("repeat")
                    if isinstance(rep, int) and rep > 1 and st[arm].present:
                        de = run.comp_arm(sid, arm).get("digests_equal")
                        if de is True:
                            text += " ="
                        elif de is False:
                            text += " \\gne"
                cs.append(text)
            second = esc(_variant(meta[sid]))
            src = run.scenario(sid).get("derived_from")
            if src:
                second = "%s · %s" % (latex.sid(src, units.unit_of(src)), second)
            rows.append([latex.sid(sid, detail=True) if fns[sid] else latex.breakable_id(sid), second] + cs)
        rows += _relation_rows(ctx, u.get("relations") or [])
        out += latex.pairs_table(rows, units_ids(ctx, members), "test", legend)
    note = mcp_unit_note(ctx, chapter, uid, members)
    out += note
    out += _unit_blocks(ctx, chapter, uid, members, unit_note=bool(note))
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
    out = latex.section("%s\\qd %s" % (latex.title_tt(fid), esc(title)), "fam:%s" % fid, level="subsection",
                        needspace=10)
    out += inline(ctx.fill(fp["prose"], ctx.family_values(fid))) + "\n\n" if fp.get("prose") else latex.placeholder(
        "(family prose: to be written)") + "\n\n"
    sids, relations = [], []
    for uid in f["units"]:
        out += unit_section(ctx, chapter, uid)
        sids.extend(units.units[uid]["members"])
        relations.extend(units.units[uid].get("relations") or [])
    out += family_summary(ctx, fid, sids, relations)
    return out


def chapter(ctx, ch):
    """A chapter: its conclusion first (the `summary` prose and the summary block), then what it leaves untested
    (the `coverage_gap` prose and the computed lines), then the background (`intro`) and the families (the owner's
    review: lead each chapter with its conclusion and coverage gap)."""
    units = ctx.units
    ctx.footnotes[ch] = cells.Footnotes()
    prose = ctx.prose.get(ch) or {}
    title = prose.get("title") or units.titles.get(ch, ch)
    out = latex.section("%s\\qd %s" % (esc(ch), esc(title)), "ch:%s" % ch, mark="%s · %s" % (ch, title))
    fams = units.chapter_families(ch)
    sids = units.chapter_scenarios(ch)
    body = ""
    for fid in fams:
        body += family_section(ctx, ch, fid)
    n_units = sum(len(units.families[f]["units"]) for f in fams)
    values = ctx.chapter_values(ch)
    # the summary is built after the families so the footnotes are assigned
    out += "\\parhead{Conclusion}\n"
    summary = prose.get("summary")
    out += paragraphs([ctx.fill(p, values) for p in ([summary] if isinstance(summary, str) else (summary or []))]) \
        if summary else latex.placeholder("(chapter summary: to be written)") + "\n\n"
    out += chapter_summary(ctx, ch, sids, n_units, len(fams))
    out += "\\parhead{What this chapter leaves untested}\n"
    gap = prose.get("coverage_gap")
    out += paragraphs([ctx.fill(p, values) for p in ([gap] if isinstance(gap, str) else (gap or []))]) \
        if gap else latex.placeholder("(coverage gap: to be written)") + "\n\n"
    lines = untested_lines(ctx, ch)
    if lines:
        out += "\\begin{itemize}\\setlength\\itemsep{0pt}\n" + "".join("\\item %s\n" % l for l in lines) + "\\end{itemize}\n"
    else:
        out += "Every property check of this chapter was passed by at least one execution of each test.\n\n"
    out += "\\parhead{Background}\n"
    out += paragraphs([ctx.fill(p, values) for p in (prose.get("intro") or [])]) if prose.get("intro") else \
        latex.placeholder("(chapter prose: to be written)") + "\n\n"
    # the compact record key right before the first family's tables (re-verification N-1)
    out += record_key_paragraph(ctx, ch)
    return out + body


def part_two(ctx):
    out = "\\part*{Part II\\qd Results by area}\\addcontentsline{toc}{part}{Part II\\qd Results by area}\n"
    for ch in ctx.chapters:
        out += chapter(ctx, ch)
    return out


# -- the closing section -------------------------------------------------------------------------------------------

PHASES = [("submit_ms", "submit"), ("queue_ms", "queue"), ("backend_ms", "backend"),
          ("terminal_lag_ms", "terminal lag"), ("paging_ms", "paging")]


def _first_warm(run, sid, arm):
    """The first host-warm execution of a scenario on an arm that did not follow an interrupted one, or None."""
    for l in run.arm_lines(sid, arm):
        if cells.host_cache_of(l) == "warm" and cells._num(l.get("wall_ms")) is not None:
            return l
    return None


def paired_walls(ctx):
    """Per scenario the first host-warm wall per arm: [(sid, chapter, dl, dr, mcp)] with None where an arm has
    none; the plots' and the observations' pairs."""
    run = ctx.run
    out = []
    for sid in run.order:
        ch = (ctx.units.by_scenario.get(sid) or {}).get("chapter") or "?"
        walls = []
        for arm in ARMS:
            l = _first_warm(run, sid, arm)
            walls.append(l.get("wall_ms") if l else None)
        out.append((sid, ch, walls[0], walls[1], walls[2]))
    return out


def phase_population(ctx):
    """{phase key: [ms]} over the population of compare's summary.mcp_phases: per scenario the host-warm MCP
    execution with a backend timing (repeat 1 preferred), as compare.mcp_overhead picks it, so that the plot and
    C.3's table read one population (VERIFY-R2 M4). Zeros are kept here; the log-scale box leaves them out."""
    run = ctx.run
    values = dict((k, []) for k, _ in PHASES)
    for sid in run.order:
        lines = dict((l.get("repeat"), l) for l in run.arm_lines(sid, "mcp") if isinstance(l.get("repeat"), int))
        timed = dict((k, l) for k, l in lines.items() if cells.sc.backend_timed(l))
        if not timed:
            continue
        m = cmp.first_with(timed, "warm", prefer=1 if 1 in timed else None)
        ph = (m or {}).get("mcp_phases") if isinstance((m or {}).get("mcp_phases"), dict) else None
        if not ph:
            continue
        for k, _ in PHASES:
            v = ph.get(k)
            if cells._num(v) is not None:
                values[k].append(v)
    return values


def phase_quantiles(ctx):
    """[(phase, p10, p25, p50, p75, p90, n, n_zero)] over phase_population: the quantiles of the positive values
    (a log scale cannot hold a zero), n the positive ones, n_zero the zeros left out."""
    values = phase_population(ctx)
    out = []
    for k, name in PHASES:
        vals = [x for x in values[k] if x > 0]
        zeros = sum(1 for x in values[k] if x == 0)
        if len(vals) < 2:
            continue
        out.append((name, cells.percentile(vals, 10), cells.percentile(vals, 25), cells.percentile(vals, 50),
                    cells.percentile(vals, 75), cells.percentile(vals, 90), len(vals), zeros))
    return out


def plot_data(ctx):
    """The data files of the closing section's plots (deterministic text): the paired walls per arm pair and the
    phase quantiles; {file name: text}."""
    pairs = paired_walls(ctx)
    files = {}
    rows = ["scenario,chapter,dl_ms,mcp_ms"] + ["%s,%s,%.3f,%.3f" % (sid, ch, dl, mc) for sid, ch, dl, dr, mc in pairs
                                               if dl is not None and mc is not None]
    files[PLOT_FILES["paired_dl_mcp"]] = "\n".join(rows) + "\n"
    rows = ["scenario,chapter,dl_ms,dr_ms"] + ["%s,%s,%.3f,%.3f" % (sid, ch, dl, dr) for sid, ch, dl, dr, mc in pairs
                                              if dl is not None and dr is not None]
    files[PLOT_FILES["paired_dl_dr"]] = "\n".join(rows) + "\n"
    rows = ["phase,n,n_zero,p10_ms,p25_ms,p50_ms,p75_ms,p90_ms"] + [
        "%s,%d,%d,%.3f,%.3f,%.3f,%.3f,%.3f" % (name.replace(" ", "_"), n, n_zero, p10, p25, p50, p75, p90)
        for name, p10, p25, p50, p75, p90, n, n_zero in phase_quantiles(ctx)]
    files[PLOT_FILES["phases"]] = "\n".join(rows) + "\n"
    return files


# The conditions of C.1's interpretive clauses (VERIFY-R4, issue 2): the median backend overhead is small against the
# median MCP overhead below this share of it; the MCP cold-warm median ratio is close to one within this margin.
BACKEND_SHARE_SMALL = 0.5
COLD_WARM_CLOSE = 0.25


def closing_values(ctx):
    """The placeholders of the closing section's observation templates (closing.json `observations`), computed
    from the run: every value a formatted string, None when the run has no data for it; the flags of the templates'
    conditionals (`backend_share_small`, `paired_seconds_to_milliseconds`, `cold_warm_mcp_close`, ...) are bools,
    False without data, so that an interpretive clause prints only when the run's figures make it true."""
    if ctx._closing_values is not None:
        return ctx._closing_values
    run = ctx.run
    s = run.summary()
    v = {}

    def stats(prefix, d, signed=False):
        d = d if isinstance(d, dict) and d.get("n") else None
        f = fmt_delta if signed else fmt_t
        v[prefix + "_p50"] = f(d.get("p50")) if d else None
        v[prefix + "_p90"] = f(d.get("p90")) if d else None
        v[prefix + "_max"] = f(d.get("max")) if d else None
        v[prefix + "_n"] = str(d.get("n")) if d else None

    mo = s.get("mcp_overhead") or {}
    stats("mcp_over_local", mo.get("minus_local_ms"), True)
    stats("mcp_over_remote", mo.get("minus_remote_ms"), True)
    stats("mcp_over_local_no_queue", mo.get("minus_local_no_queue_ms"), True)
    stats("mcp_over_remote_no_queue", mo.get("minus_remote_no_queue_ms"), True)
    stats("backend_over_local", mo.get("backend_minus_local_ms"), True)
    stats("backend_over_remote", mo.get("backend_minus_remote_ms"), True)
    ph = s.get("mcp_phases") or {}
    for k, name in PHASES:
        stats("phase_" + name.replace(" ", "_"), ph.get(k))
    stats("phase_backend_pattern", ph.get("backend_pattern_ms"))
    stats("mcp_cold_wall", s.get("mcp_cold_wall_ms"))
    cw = s.get("cold_warm") or {}
    for arm in ARMS:
        d = cw.get(arm) or {}
        key = cells.ARM_SHORT[arm]
        v["cold_warm_%s_n" % key] = str(d.get("n")) if d.get("n") else None
        v["cold_warm_%s_delta_p50" % key] = fmt_delta((d.get("delta_ms") or {}).get("p50")) if d.get("n") else None
        v["cold_warm_%s_delta_p90" % key] = fmt_delta((d.get("delta_ms") or {}).get("p90")) if d.get("n") else None
        v["cold_warm_%s_ratio_p50" % key] = fmt_n((d.get("ratio") or {}).get("p50"), 2) if d.get("n") else None
    fo = s.get("front_overhead") or {}
    for arm in ("direct-local", "direct-remote"):
        stats("connect_" + cells.ARM_SHORT[arm], (fo.get("connect_ms") or {}).get(arm))
    nets = [(e.get("front_overhead") or {}).get("net_ms") for e in (run.comp.get("scenarios") or {}).values()
            if isinstance(e, dict) and isinstance(e.get("front_overhead"), dict)]
    nets = [x for x in nets if cells._num(x) is not None]
    v["front_net_p50"] = fmt_delta(cells.percentile(nets, 50)) if nets else None
    v["front_net_p90"] = fmt_delta(cells.percentile(nets, 90)) if nets else None
    v["front_net_n"] = str(len(nets)) if nets else None
    pairs = paired_walls(ctx)
    for key, idx in (("dl", 2), ("dr", 3), ("mcp", 4)):
        vals = [p[idx] for p in pairs if p[idx] is not None]
        v["warm_wall_%s_p50" % key] = fmt_t(cells.percentile(vals, 50)) if vals else None
        v["warm_wall_%s_p90" % key] = fmt_t(cells.percentile(vals, 90)) if vals else None
        v["warm_wall_%s_n" % key] = str(len(vals)) if vals else None
    both = [p for p in pairs if p[2] is not None and p[4] is not None]
    v["paired_dl_mcp_n"] = str(len(both)) if both else None
    ratios = [p[4] / p[2] for p in both if p[2] > 0]
    v["paired_dl_mcp_ratio_p50"] = fmt_ratio(cells.percentile(ratios, 50)) if ratios else None
    v["paired_dl_mcp_ratio_p10"] = fmt_ratio(cells.percentile(ratios, 10)) if ratios else None
    v["paired_dl_mcp_ratio_p90"] = fmt_ratio(cells.percentile(ratios, 90)) if ratios else None
    # the two sides' medians of the paired population, and whether they print in seconds against milliseconds
    dl_med = cells.percentile([p[2] for p in both], 50) if both else None
    mcp_med = cells.percentile([p[4] for p in both], 50) if both else None
    v["paired_dl_wall_p50"] = fmt_t(dl_med) if both else None
    v["paired_mcp_wall_p50"] = fmt_t(mcp_med) if both else None
    v["paired_seconds_to_milliseconds"] = bool(both) and cells._split_unit(fmt_t(mcp_med))[1] == "s" \
        and cells._split_unit(fmt_t(dl_med))[1] == "ms"
    # the backend overhead against the MCP overhead (both against direct-local, medians of separate distributions)
    over = cells._num((mo.get("minus_local_ms") or {}).get("p50")) if isinstance(mo.get("minus_local_ms"), dict) \
        and (mo.get("minus_local_ms") or {}).get("n") else None
    be = cells._num((mo.get("backend_minus_local_ms") or {}).get("p50")) \
        if isinstance(mo.get("backend_minus_local_ms"), dict) and (mo.get("backend_minus_local_ms") or {}).get("n") else None
    v["mcp_over_local_positive"] = over is not None and over > 0
    v["backend_share_small"] = over is not None and be is not None and over > 0 and be < BACKEND_SHARE_SMALL * over
    # the MCP cold-warm median ratio close to one
    cw_mcp = cw.get("mcp") or {}
    cw_ratio = cells._num((cw_mcp.get("ratio") or {}).get("p50")) if cw_mcp.get("n") else None
    v["cold_warm_mcp_close"] = cw_ratio is not None and abs(cw_ratio - 1.0) <= COLD_WARM_CLOSE
    largest = None  # (p50, phase name) of the phase with the largest median (backend_pattern left out)
    for k, name in PHASES:
        p50 = cells._num((ph.get(k) or {}).get("p50")) if isinstance(ph.get(k), dict) else None
        if p50 is not None and (largest is None or p50 > largest[0]):
            largest = (p50, name)
    v["mcp_phase_largest"] = largest[1] if largest else None
    v["mcp_phase_largest_p50"] = fmt_t(largest[0]) if largest else None
    fo_ = s.get("front_overhead") or {}
    v["scenarios_without_warm_pair"] = str(fo_["scenarios_without_warm_pair"]) \
        if cells._num(fo_.get("scenarios_without_warm_pair")) is not None else None
    v["warm_pairs_left_out_after_interrupted"] = str(fo_["warm_pairs_left_out_after_interrupted"]) \
        if cells._num(fo_.get("warm_pairs_left_out_after_interrupted")) is not None else None
    # the scenarios the MCP pairings leave out: every scenario less the pairs (VERIFY-R2 H4: the direct pairing's
    # count is C.4's, not C.3's)
    n_local = cells._num((mo.get("minus_local_ms") or {}).get("n")) if isinstance(mo.get("minus_local_ms"), dict) else None
    n_remote = cells._num((mo.get("minus_remote_ms") or {}).get("n")) if isinstance(mo.get("minus_remote_ms"), dict) else None
    v["scenarios_without_mcp_local_pair"] = str(len(run.order) - int(n_local)) if n_local is not None else None
    v["scenarios_without_mcp_remote_pair"] = str(len(run.order) - int(n_remote)) if n_remote is not None else None
    # every host-cold MCP execution (C.4's host-cold table), against the cold calls of the paired scenarios
    # (summary.mcp_cold_wall_ms, C.3); the two populations differ and the prose says so (VERIFY-R2 M4)
    stats("cold_mcp_all", ((fo_.get("cold") or {}).get("mcp") or {}).get("wall_ms"))
    n_all = cells._num((((fo_.get("cold") or {}).get("mcp") or {}).get("wall_ms") or {}).get("n"))
    n_timed = cells._num((s.get("mcp_cold_wall_ms") or {}).get("n")) if isinstance(s.get("mcp_cold_wall_ms"), dict) else None
    v["cold_mcp_untimed_n"] = str(int(n_all) - int(n_timed)) if n_all is not None and n_timed is not None else None
    v["cold_mcp_untimed_text"] = cold_mcp_untimed_text(run)
    # the phases' population: the warm execution of each paired scenario (summary.mcp_phases); paging is 0 ms
    # where the results came with the terminal status, and the log-scale box leaves those out
    phase_pop = phase_population(ctx)
    paging = phase_pop.get("paging_ms") or []
    v["phase_paging_zero_n"] = str(sum(1 for x in paging if x == 0)) if paging else None
    v["phase_paging_nonzero_n"] = str(sum(1 for x in paging if x > 0)) if paging else None
    v["load_caveat_line"] = load_caveat_line(run)
    diffs = [p[4] - p[2] for p in both]
    v["paired_dl_mcp_delta_p50"] = fmt_delta(cells.percentile(diffs, 50)) if diffs else None
    v["paired_dl_mcp_delta_p90"] = fmt_delta(cells.percentile(diffs, 90)) if diffs else None
    # the pairs whose two walls lie within 1 s of each other, in either direction (VERIFY-R2 H3: not d < 1000,
    # which counted the pairs where MCP was faster by any margin)
    v["paired_dl_mcp_under_1s"] = str(sum(1 for d in diffs if abs(d) < 1000)) if diffs else None
    v["paired_dl_mcp_faster"] = str(sum(1 for d in diffs if d < 0)) if diffs else None
    both_dr = [p for p in pairs if p[2] is not None and p[3] is not None]
    v["paired_dl_dr_n"] = str(len(both_dr)) if both_dr else None
    ddiffs = [p[3] - p[2] for p in both_dr]
    v["paired_dl_dr_delta_p50"] = fmt_delta(cells.percentile(ddiffs, 50)) if ddiffs else None
    v["paired_dl_dr_delta_p90"] = fmt_delta(cells.percentile(ddiffs, 90)) if ddiffs else None
    lf = run.load or {}
    v["load_samples"] = str(lf.get("n", 0))
    v["load_sra_mean"] = fmt_n(lf.get("sra_mean")) if lf.get("n") else None
    v["load_sra_max"] = fmt_n(lf.get("sra_max")) if lf.get("n") else None
    v["load_load1_mean"] = fmt_n(lf.get("load_mean"), 1) if lf.get("n") else None
    v["load_load1_max"] = fmt_n(lf.get("load_max"), 1) if lf.get("n") else None
    busy = (lf.get("sra_mean") or 0) > LOAD_CAVEAT_SRA or (lf.get("load_mean") or 0) > LOAD_CAVEAT_LOAD1
    v["load_caveat"] = "applies" if busy else "does not apply"
    v["load_caveat_sra_cpu_pct"] = fmt_n(LOAD_CAVEAT_SRA, 0)
    v["load_caveat_load1"] = fmt_n(LOAD_CAVEAT_LOAD1, 0)
    wai = 0
    for area, t in (s.get("timing") or {}).items():
        for arm in ARMS:
            b = ((t.get(arm) or {}).get("by_host_cache") or {}).get(cells.WARM_AFTER_INTERRUPTED) or {}
            wai += (b.get("wall") or {}).get("n") or 0
    v["warm_after_interrupted_n"] = str(wai)
    hc = run.comp.get("host_cache") or {}
    v["warm_window_s"] = fmt_n(hc.get("warm_window_s", 60), 0)
    totals = s.get("totals") or {}
    for arm in ARMS:
        hcc = (totals.get(arm) or {}).get("host_cache") or {}
        for state in ("cold", "warm", "aged", "not_reached"):
            v["executions_%s_%s" % (cells.ARM_SHORT[arm], state)] = str(hcc.get(state, 0))
    v["n_scenarios"] = str(len(run.order))
    ctx._closing_values = v
    return v


_CODE_WORDS = {"resolve_timeout": "the service's own deadline"}
_OUTCOME_WORDS = {"service_by_design": "by design", "inconclusive": "inconclusive", "not_implemented_mcp": "not exercised",
                  "finding": "a finding", "fail": "a failure", "error": "an error", "pass": "a pass"}


def cold_mcp_untimed_split(run):
    """The host-cold MCP executions without a backend timing (the calls C.4's cold column counts and C.3's cold
    wall leaves out: cold_mcp_untimed_n), grouped by what ended them, as [(n, outcome word, words)] with the
    largest group first: the words are the execution's recorded time limitation (checks.time_limited, e.g. 'the
    host's /pattern 503 deadline relayed by the service'), else its error code in words (a resolve_timeout is the
    service's own deadline), else 'no timing recorded' (re-verification of round 3, L-b: the split is computed
    from the record, not asserted)."""
    groups = {}
    for sid in run.order:
        for l in run.arm_lines(sid, "mcp"):
            if cells.host_cache_of(l) != "cold" or cells.sc.backend_timed(l) or cells._num(l.get("wall_ms")) is None:
                continue
            code = (l.get("error") or {}).get("code") if isinstance(l.get("error"), dict) else None
            tl = str(l.get("time_limited") or "")
            if tl:
                words = tl
            elif code:
                words = ("%s, `%s`" % (_CODE_WORDS[code], code)) if code in _CODE_WORDS else "`%s`" % code
            else:
                words = "no timing recorded"
            key = (_OUTCOME_WORDS.get(l.get("outcome"), str(l.get("outcome"))), words)
            groups[key] = groups.get(key, 0) + 1
    return sorted(((n, o, w) for (o, w), n in groups.items()), key=lambda t: (-t[0], t[1], t[2]))


def cold_mcp_untimed_text(run):
    """'23 by design (the host's /pattern 503 deadline relayed by the service), 3 inconclusive (the host's /resolve
    503 deadline relayed by the service) and 1 inconclusive (the service's own deadline, `resolve_timeout`)'."""
    parts = ["%d %s (%s)" % (n, o, w) for n, o, w in cold_mcp_untimed_split(run)]
    return facts._join(parts) if parts else "none"


def fmt_ratio(x):
    """A ratio with one decimal below 10, whole above (275, not 275.0)."""
    v = cells._num(x)
    if v is None:
        return "–"
    return fmt_n(v, 1 if v < 10 else 0)


def closing_section_prose(ctx, key, fallback):
    """The method / observation / limits paragraphs of a closing table section (closing.json `sections[key]`),
    filled; the generator's own line when the prose has none."""
    sec = (ctx.closing.get("sections") or {}).get(key) if isinstance(ctx.closing.get("sections"), dict) else None
    items = (sec or {}).get("paragraphs") or []
    if isinstance(items, str):
        items = [items]
    if items:
        return paragraphs([ctx.fill(p, closing_values(ctx)) for p in items])
    return fallback


def observations(ctx):
    """The observation sentences of closing.json filled from the run (the placeholders and the conditionals of
    facts.fill, the flags of closing_values); a template whose values the run lacks is replaced by a note naming
    them (never a sentence with a hole)."""
    items = ctx.closing.get("observations") or []
    values = closing_values(ctx)
    out = []
    for it in items:
        template = it.get("template") if isinstance(it, dict) else it
        if not template:
            continue
        names = re.findall(r"\{([a-z_0-9]+)\}", str(template))
        missing = [n for n in names if values.get(n) is None]
        if missing:
            out.append("{\\color{black!55}\\itshape (observation %s not computed for this run: no data for %s)}" % (
                esc((it.get("id") if isinstance(it, dict) else None) or "?"), esc(", ".join(missing))))
            continue
        out.append(inline(ctx.fill(str(template), values)))
    return out


def closing(ctx):
    cp = ctx.closing
    title = cp.get("title") or "Cross-cutting timing"
    out = "\\part*{Closing\\qd %s}\\addcontentsline{toc}{part}{Closing\\qd %s}\n" % (esc(title), esc(title))
    out += "\\markright{%s}\n" % esc(title)
    s = ctx.run.summary()
    # C.1 observations: generator-filled sentences from the prose templates, the caveats beside them
    out += latex.section("C.1\\qd Observations", "sec:C1")
    intro = cp.get("intro") or []
    out += paragraphs([ctx.fill(p, closing_values(ctx)) for p in ([intro] if isinstance(intro, str) else intro)]) \
        if intro else latex.placeholder("(closing intro: to be written)") + "\n\n"
    obs = observations(ctx)
    if obs:
        out += "\\begin{itemize}\\setlength\\itemsep{1pt}\n" + "".join("\\item %s\n" % o for o in obs) + "\\end{itemize}\n"
    else:
        out += latex.placeholder("(observations: to be written as templates in closing.json)") + "\n\n"
    caveats = cp.get("caveats") or []
    out += "\\parhead{Caveats}\n"
    out += paragraphs([ctx.fill(p, closing_values(ctx)) for p in ([caveats] if isinstance(caveats, str) else caveats)]) \
        if caveats else latex.placeholder("(caveats: the shared host load and the cold / warm heuristic: to be "
                                          "written)") + "\n\n"
    # C.2 plots
    out += latex.section("C.2\\qd Paired latencies and the MCP phases", "sec:C2")
    plots = cp.get("plots") if isinstance(cp.get("plots"), dict) else {}
    values = closing_values(ctx)
    ctx.data_files.update(plot_data(ctx))
    pairs = paired_walls(ctx)
    if any(p[2] is not None and p[4] is not None for p in pairs):
        cap = plots.get("paired_latency") or (
            "Paired host-warm walls per scenario: the warm MCP call against direct-local (left, n {paired_dl_mcp_n}) "
            "and direct-remote against direct-local (right, n {paired_dl_dr_n}); log scales, the dashed line y = x. "
            "One point per scenario; no aggregate is summed.")
        body = ("\\begin{minipage}[t]{0.49\\textwidth}\\centering\n%s\\end{minipage}\\hfill"
                "\\begin{minipage}[t]{0.49\\textwidth}\\centering\n%s\\end{minipage}" % (
                    latex.paired_plot(PLOT_FILES["paired_dl_mcp"], "dl_ms", "mcp_ms", "direct-local wall (ms)",
                                      "mcp warm wall (ms)", "fig:paired-mcp"),
                    latex.paired_plot(PLOT_FILES["paired_dl_dr"], "dl_ms", "dr_ms", "direct-local wall (ms)",
                                      "direct-remote wall (ms)", "fig:paired-dr")))
        out += latex.figure_block(body, inline(ctx.fill(cap, values)))
    boxes = phase_quantiles(ctx)
    if boxes:
        cap = plots.get("phases") or (
            "The MCP phases of the host-warm calls whose backend was reached, one box per phase: the box p25 to p75, "
            "the line p50, the whiskers p10 and p90, log scale (n per phase in the labels). The phases are drawn "
            "side by side, not summed: their medians do not add up to an end-to-end total.")
        out += latex.figure_block(latex.phase_boxplot(boxes, "fig:phases"), inline(ctx.fill(cap, values)))
    # C.3 .. C.7: the tables
    out += latex.section("C.3\\qd MCP overhead decomposition", "sec:C3")
    out += closing_section_prose(ctx, "C3", (
        "Pairs of host-warm executions only: the MCP warm execution (repeat 1 preferred) minus each direct arm's "
        "first warm execution; the phases as the service measured them. In these tables a / b never appears: "
        "p50, p90 and the maximum are columns.\n\n"))
    ph = s.get("mcp_phases") or {}
    rows = [[esc(k.replace("_ms", "").replace("_", " "))] + _stats_cells(ph.get(k)) for k in
            ("submit_ms", "queue_ms", "backend_ms", "backend_pattern_ms", "terminal_lag_ms", "paging_ms") if k in ph]
    out += latex.plain_table("@{}l r r r r@{}", ["phase", "p50", "p90", "max", "n"], rows, "MCP phases")
    mo = s.get("mcp_overhead") or {}
    rows = [[esc(k.replace("_ms", "").replace("_", " "))] + _stats_cells(mo.get(k), signed=True) for k in sorted(mo)]
    out += latex.plain_table("@{}l r r r r@{}", ["MCP minus direct", "p50", "p90", "max", "n"], rows,
                             "MCP overhead (warm MCP wall minus each direct arm; the backend time minus each direct arm); "
                             "signed")
    if s.get("mcp_cold_wall_ms"):
        out += latex.plain_table("@{}l r r r r@{}", ["", "p50", "p90", "max", "n"],
                                 [["MCP cold wall"] + _stats_cells(s["mcp_cold_wall_ms"])], "MCP cold wall")
    out += latex.section("C.4\\qd Front cost", "sec:C4")
    out += closing_section_prose(ctx, "C4", "")
    fo = s.get("front_overhead") or {}
    rows = [[esc(arm)] + _stats_cells((fo.get("connect_ms") or {}).get(arm)) for arm in sorted(fo.get("connect_ms") or {})]
    out += latex.plain_table("@{}l r r r r@{}", ["connection set-up", "p50", "p90", "max", "n"], rows,
                             "connect\\_ms per arm")
    pan = fo.get("per_area_net") or {}
    rows = [[esc(a)] + _stats_cells(pan.get(a), signed=True)
            for a in sorted(pan, key=lambda x: (ctx.units.chapters.index(x) if x in ctx.units.chapters else 99, x))]
    out += latex.plain_table("@{}l r r r r@{}", ["chapter", "p50", "p90", "max", "n"], rows,
                             "direct-remote minus direct-local, net of the server time, per chapter; signed")
    cold = fo.get("cold") or {}
    rows = []
    for arm in sorted(cold):
        for kind in ("wall_ms", "server_ms"):
            if kind in cold[arm]:
                rows.append([esc(arm), esc(kind.replace("_ms", ""))] + _stats_cells(cold[arm][kind]))
    out += latex.plain_table("@{}l l r r r r@{}", ["arm", "", "p50", "p90", "max", "n"], rows, "host-cold executions")
    out += latex.section("C.5\\qd Cold versus warm", "sec:C5")
    hc = ctx.run.comp.get("host_cache") or {}
    out += closing_section_prose(ctx, "C5", (
        "Cold and warm are the recorded heuristic of I.2 (cold: the scenario's first execution to reach the host "
        "in this run; warm: one within %s s of another; residency is never measured), so a delta is a difference "
        "between two executions of the same scenario in that order, not a measured cache effect.\n\n" % esc(
            fmt_n(hc.get("warm_window_s", 60), 0))))
    cw = s.get("cold_warm") or {}
    rows = []
    for arm in sorted(cw):
        d = cw[arm] or {}
        rows.append([esc(arm), str(d.get("n", 0)), fmt_delta((d.get("delta_ms") or {}).get("p50")),
                     fmt_delta((d.get("delta_ms") or {}).get("p90")), fmt_n((d.get("ratio") or {}).get("p50"), 2),
                     fmt_n((d.get("ratio") or {}).get("p90"), 2)])
    out += latex.plain_table("@{}l r r r r r@{}", ["arm", "pairs", "delta p50", "delta p90", "ratio p50", "ratio p90"],
                             rows, "warm minus cold per scenario (delta, signed: negative when warm was faster), "
                                   "warm over cold (ratio)")
    out += latex.section("C.6\\qd Load telemetry", "sec:C6")
    lf = ctx.run.load or {}
    out += closing_section_prose(ctx, "C6", (
        "%d samples; sra\\_all CPU mean %s\\%% (max %s\\%%), load1 mean %s (max %s). The Logan caveat: above "
        "%d\\%% sra\\_all CPU or load1 %d the host was busy with other work and the timings were "
        "observed under concurrent host load.\n\n" % (
            lf.get("n", 0), fmt_n(lf.get("sra_mean")), fmt_n(lf.get("sra_max")), fmt_n(lf.get("load_mean"), 1),
            fmt_n(lf.get("load_max"), 1), LOAD_CAVEAT_SRA, LOAD_CAVEAT_LOAD1)))
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
    out += latex.section("C.7\\qd Per-chapter timing", "sec:C7")
    out += closing_section_prose(ctx, "C7", "")
    timing = s.get("timing") or {}
    rows = []
    wai = 0
    for area in sorted(timing, key=lambda x: (ctx.units.chapters.index(x) if x in ctx.units.chapters else 99, x)):
        for arm in ARMS:
            t = (timing[area].get(arm) or {}).get("by_host_cache") or {}
            for hc_ in ("cold", "warm"):
                b = t.get(hc_) or {}
                w, sv = b.get("wall") or {}, b.get("server") or {}
                if not w.get("n") and not sv.get("n"):
                    continue
                rows.append([esc(area), esc(arm), esc(hc_)] + _stats_cells(w) + _stats_cells(sv))
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
    short = ctx.short_reasons()
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
            service.append("“%s”" % latex.join_units(esc_wrap(sa["text_contains"], 240)))
        if sa.get("missing_argument"):
            service.append("missing argument %s" % tt(sa["missing_argument"]))
        if c.get("explains_checks"):
            service.append("explains %s" % esc(", ".join(c["explains_checks"])))
        direct = []
        if de.get("http") is not None:
            direct.append("http %s" % esc(de["http"]))
        if de.get("text"):
            direct.append(latex.join_units(tt(de["text"], 240)))
        elif de.get("text_prefix"):
            direct.append("prefix %s" % latex.join_units(tt(de["text_prefix"], 240)))
        marks = {}
        unverified_detail = []
        for sid in naming.get(cid, []):
            ma = run.comp_entry(sid).get("mcp_agreement") or {}
            if ma.get("case") == cid and ma.get("class") == "service_by_design":
                mark = "verified" if ma.get("verified") else "unverified"
                if not ma.get("verified"):
                    dis = ma.get("disagreement") if isinstance(ma.get("disagreement"), dict) else {}
                    if dis.get("direct_detail"):
                        unverified_detail.append("%s: %s" % (sid, dis["direct_detail"]))
            elif ma.get("case") == cid:
                mark = (ma.get("class") or "").replace("_", " ")
                dis = ma.get("disagreement") if isinstance(ma.get("disagreement"), dict) else {}
                if dis.get("case") == cid and dis.get("direct_verified") is False and dis.get("direct_detail"):
                    unverified_detail.append("%s (direct side not verified): %s" % (sid, dis["direct_detail"]))
            else:
                # the service's refusal matched this case on an arm classed otherwise (not implemented by the table)
                matched = any(l.get("service_by_design_case") == cid for l in run.arm_lines(sid, "mcp"))
                mark = "matched by the answer" if matched else "named"
                dis = ma.get("disagreement") if isinstance(ma.get("disagreement"), dict) else {}
                if matched and dis.get("case") == cid and dis.get("direct_verified") is False and dis.get("direct_detail"):
                    unverified_detail.append("%s (direct side not verified): %s" % (sid, dis["direct_detail"]))
            marks.setdefault(mark, []).append(sid)
        named = []
        for mark in sorted(marks):
            sids_ = marks[mark]
            if len(sids_) <= 6:
                named.append("%s: %s" % (esc(mark), ", ".join(latex.sid(s, ctx.units.unit_of(s)) for s in sids_)))
            else:
                named.append("%s: %d (%s, …)" % (esc(mark), len(sids_),
                                                 ", ".join(latex.sid(s, ctx.units.unit_of(s)) for s in sids_[:3])))
        first = [tt(cid) + "\\newline{\\scriptsize %s}" % esc(c.get("agreement") or "–")]
        matches = latex.join_units(esc_wrap(c.get("scenario_match"), 600))
        if short.get(cid):
            matches += "\\newline{\\scriptsize in short: %s}" % inline(short[cid])
        rows.append(first + [matches, "; ".join(service) or "–", "; ".join(direct) or "–", "; ".join(named) or "–"])
        if unverified_detail:
            rows[-1] = (rows[-1], "{\\footnotesize not verified on the direct side: %s}" % latex.join_units(esc_wrap(
                "; ".join(unverified_detail))))
    out = ("One row per case of the service-by-design table that a scenario of the built chapters names, that an "
           "MCP answer matched, or that belongs to their areas: the case and its agreement class, what it matches "
           "(and the short reason the tables use), the service's answer, the direct expectation, and the scenarios "
           "naming it with their verification; where the direct side did not match the case, the comparison's "
           "detail is under the row.\n\n")
    out += latex.plain_table(latex.spec_p([30, 38, 40, 36, 30]),
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
    """Appendix D: the full record key first (labelled rk:D, linked from every matrix's legend and from the
    compact key of I.3: re-verification N-1), then the short glossary of glossary.json ({term, definition | text,
    see}) in alphabetical order (a leading article ignored), two columns, with the version-pinned README reference
    for the long form."""
    entries = [e for e in ((ctx.glossary or {}).get("entries") or []) if isinstance(e, dict) and e.get("term")]
    ref = ctx.glossary.get("readme") if isinstance(ctx.glossary.get("readme"), dict) else {}
    commit = ctx.readme_commit()
    out = "\\parhead{Record key}\\phantomsection\\label{rk:D}\n"
    rk = ctx.glossary.get("records_key") if isinstance(ctx.glossary.get("records_key"), dict) else {}
    if rk.get("lead"):
        out += paragraphs([ctx.fill(rk["lead"])])
    out += record_key(ctx, full=True)
    # the glossary stays one page (the owner's review, point 3): it asks for a page's space and so starts on a
    # fresh page when the record key's table ends mid-page
    out += "\\Needspace{56\\baselineskip}\n\\parhead{Glossary}\n"
    lead = ctx.glossary.get("lead")
    if lead:
        out += paragraphs([ctx.fill(lead)])
    else:
        out += ("The short form: the terms needed to read the evidence. The long form of every term, with the rules "
                "behind the classes, is the suite README's section ‘%s’ (%s), pinned to commit %s%s.\n\n" % (
                    esc(ref.get("section") or README_SECTION), tt(ref.get("path") or README_PATH),
                    tt(commit) if commit else "(not determined: no git checkout at build time)",
                    " of this repository" if commit else ""))
    if not entries:
        out += latex.placeholder("(glossary: to be written)") + "\n\n"
    else:
        # one page for the entries (the owner's review, point 3): a size between footnote and script, two columns
        out += ("\\begin{multicols}{2}\\fontsize{7.5}{8.8}\\selectfont\\raggedright\\setlength{\\columnsep}{5mm}"
                "\\setlength{\\parskip}{0pt}\n")
        for e in sorted(entries, key=lambda e: _sort_key(e["term"])):
            text = inline(ctx.fill(e.get("definition") or e.get("text") or ""))
            see = e.get("see")
            see_text = inline(", ".join(str(x) for x in see)) if isinstance(see, list) and see else None
            out += latex.glossary_entry(esc(e["term"]), text, see_text)
        out += "\\end{multicols}\n\n"
    return out


def index(ctx):
    sids = sorted(ctx.selected_scenarios())
    out = ("Every test id of the built chapters in alphabetical order, its unit and the page of the unit's table; a "
           "footnoted id also gives the page of its detail block.\n\n")
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
        rows.append(("checks (re-evaluation)", "commit %s%s; checks.py sha256 %s; suite\\_common sha256 %s; %s (the "
                                               "re-evaluation judges the recorded answers again and replaces nothing a "
                                               "runner decided without the checks)" % (
            tt(_short(checks.get("commit"), 12)), " (dirty)" if checks.get("dirty") else "",
            tt(_short(checks.get("checks_sha256"), 16)), tt(_short(checks.get("suite_common_sha256"), 16)),
            _date(reeval.get("reevaluated_utc")))))
    else:
        rows.append(("checks", "the recorded set: judged by the runners' checks (%s)" % esc(
            (run.judged or {}).get("source") or "–")))
    if table.get("path"):
        text = "version %s, %s cases, service commit %s, generated %s; sha256 %s; %s" % (
            esc(table.get("version")), esc(table.get("cases")), tt(_short(table.get("service_commit"), 12)),
            _date(table.get("generated_utc")), tt(table.get("sha256") or "–"), tt(table.get("path")))
        if table.get("replayed"):
            recorded = ((run.reeval or {}).get("by_design") or {}).get("path")
            text += "; the table the re-evaluation replayed (reevaluated.json by\\_design), found by its sha256%s" % (
                (" at this path (the re-evaluation recorded %s)" % tt(recorded))
                if recorded and recorded != table.get("path") else " at the recorded path")
        elif table.get("replayed_sha256"):
            text += "; \\textbf{WARNING}: not the table the re-evaluation replayed (sha256 %s, %s cases): %s" % (
                tt(table["replayed_sha256"]), esc(table.get("replayed_cases", "–")),
                esc("; ".join(table.get("warnings") or []) or "the replayed table was not found", 400))
        elif table.get("warnings"):
            text += "; " + esc("; ".join(table["warnings"]), 400)
        rows.append(("by-design table", text))
    elif table.get("warnings"):
        rows.append(("by-design table", "none read: %s" % esc("; ".join(table["warnings"]), 400)))
    mcp_run = run.mcp_run or {}
    for arm, info in sorted((r.get("arms") or {}).items()):
        info = info or {}
        ended = _arm_end(run, arm, info)
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
    sessions = _sessions_text(run)
    if sessions:
        rows.append(("mode", sessions))
    hl_s, hl_e = r.get("host_load_start"), r.get("host_load_end")
    rows.append(("host load (load1 / 5 / 15)", "at start %s; at end %s" % (
        esc(", ".join(fmt_n(x, 1) for x in hl_s)) if hl_s else "–",
        esc(", ".join(fmt_n(x, 1) for x in hl_e)) if hl_e else "–")))
    lf = run.load or {}
    caveat = ""
    if (lf.get("sra_mean") or 0) > LOAD_CAVEAT_SRA or (lf.get("load_mean") or 0) > LOAD_CAVEAT_LOAD1:
        caveat = (" — the Logan caveat applies: the host was busy with other work (sra\\_all above %d\\%% CPU or "
                  "load1 above %d), so the timings were observed under concurrent host load" % (
                      LOAD_CAVEAT_SRA, LOAD_CAVEAT_LOAD1))
    rows.append(("samples", "%d; sra\\_all CPU mean %s\\%%, max %s\\%%; load1 mean %s, max %s%s" % (
        lf.get("n", 0), fmt_n(lf.get("sra_mean")), fmt_n(lf.get("sra_max")), fmt_n(lf.get("load_mean"), 1),
        fmt_n(lf.get("load_max"), 1), caveat)))
    hc = comp.get("host_cache") or {}
    rows.append(("host cache", "cold = the scenario's first execution on the host in this run; warm = another "
                               "within %s s; aged = older; the MCP pair runs cold then warm, the direct arms between "
                               "them" % esc(hc.get("warm_window_s", "60"))))
    gt = manifest.get("ground_truth") or {}
    for key in sorted(gt):
        value = gt[key]
        if isinstance(value, (list, dict)):
            value = "%d entries" % len(value)
        rows.append((esc(key.replace("_", " ")), esc_wrap(value, 220)))
    return _kv_table(rows)


def gate_problems(ctx):
    """The deploy gate's problems by kind (the scenarios linked, the configuration texts each once)."""
    run = ctx.run
    gate = run.gate
    word, how, counts_text, problems_text = gate_summary(ctx)
    if not isinstance(gate, dict):
        return "\\textbf{%s}: %s.\n\n" % (esc(word), esc(how))
    out = "\\textbf{%s}: %s; %s; %s.\n\n" % (esc(word), esc(how), esc(counts_text), esc(problems_text))
    problems = [p for p in (gate.get("problems") or []) if isinstance(p, dict)]
    if problems:
        by_kind = {}
        for p in problems:
            by_kind.setdefault(p.get("kind") or "?", []).append(p)
        rows = []
        kinds = [k for k in gate_mod.KINDS if k in by_kind] + sorted(k for k in by_kind if k not in gate_mod.KINDS)
        for k in kinds:
            v = by_kind[k]
            names, details = [], []
            repeats = {}
            for p in v:
                if p.get("scenario"):
                    where = "%s %s" % (p["scenario"], p.get("arm") or "")
                    if where not in names:
                        names.append(where)
                    if p.get("repeat") is not None:
                        repeats.setdefault(where, []).append(p["repeat"])
                    continue
                text = str(p.get("detail") or "?")
                for d in details:
                    if d[0] == text:
                        d[1] += 1
                        if p.get("arm") and p["arm"] not in d[2]:
                            d[2].append(p["arm"])
                        break
                else:
                    details.append([text, 1, [p["arm"]] if p.get("arm") else []])
            parts = []
            if names:
                # every scenario-arm pair once, its repeats listed when the arm has more than one problem there
                # (r0, r1), so that the row's count (one per execution) is readable from the list
                def _name(n):
                    text = _link_sid(ctx, n.split(" ")[0])
                    arm = n.split(" ", 1)[1].strip()
                    if arm:
                        text += " " + esc(arm)
                    reps = repeats.get(n) or []
                    if len(reps) > 1:
                        text += " (%s)" % ", ".join("r%s" % esc(r) for r in sorted(reps, key=str))
                    return text
                parts.append(", ".join(_name(n) for n in names[:40]) + (
                    " … (%d in all)" % len(names) if len(names) > 40 else ""))
                n_exec = sum(1 for p in v if p.get("scenario"))
                if n_exec != len(names):
                    parts[-1] += " — %d executions in %d scenario-arm pairs" % (n_exec, len(names))
            for text, n, arms_ in details[:6]:
                parts.append("%s%s: %s" % (esc(", ".join(arms_)) + " " if arms_ else "",
                                           ("%d times" % n) if n > 1 else "once", esc_wrap(text, 240)))
            if len(details) > 6:
                parts.append("… %d more problem texts (all in gate.py's output)" % (len(details) - 6))
            rows.append([esc(k), "; ".join(parts), str(len(v))])
        out += latex.plain_table(latex.spec_p([30, 134], extra_cols=1, extra_spec="r"),
                                 ["problem kind", "scenario, arm", "n"], rows)
    for n in (gate.get("notes") or [])[:6]:
        out += "{\\small\\color{black!60}note: %s}\\par\n" % esc(n, 300)
    return out + "\n"


def mcp_limits(ctx):
    run = ctx.run
    out = ""
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
    out += "\\parhead{Cases of the by-design table of agreement not implemented via MCP that scenarios name}\n"
    out += latex.plain_table(latex.spec_p([50, 60, 52], extra_cols=1, extra_spec="r"),
                             ["case", "scenario match", "missing argument", "scenarios"], rows)
    reasons = {}
    for sid in run.order:
        mcp = run.scenario(sid).get("mcp") or {}
        if mcp.get("not_implemented"):
            reasons.setdefault(str(mcp["not_implemented"]), []).append(sid)
    out += "\\parhead{The manifest's reasons for not implemented via MCP}\n"
    rows = [[esc_wrap(r, 120), str(len(s)), esc(", ".join(s[:8]) + (" …" if len(s) > 8 else ""))]
            for r, s in sorted(reasons.items(), key=lambda kv: (-len(kv[1]), kv[0]))]
    out += latex.plain_table(latex.spec_p([60, 104], extra_cols=1, extra_spec="r"),
                             ["reason", "ids", "scenarios"], [[r[0], r[2], r[1]] for r in rows])
    dec = run.declared
    out += "\\parhead{Owner-declared exemptions}\n"
    if isinstance(dec, dict):
        out += "File %s (declared version %s).\n" % (tt(dec.get("_file") or "declared"), esc(dec.get("declared_version")))
        rows = []
        for kind in ("exclusions", "not_required", "relations_excluded"):
            for e in dec.get(kind) or []:
                if isinstance(e, dict):
                    rows.append([esc(kind.replace("_", " ")), latex.breakable_id(e.get("scenario") or e.get("relation") or "–"),
                                 esc(", ".join(e.get("arms") or [])), esc_wrap(e.get("reason"))])
        out += latex.plain_table(latex.spec_p([22, 30, 22, 96]), ["kind", "scenario", "arms", "reason"], rows)
    else:
        out += "None given for this manifest.\n\n"
    return out


def appendices(ctx):
    out = "\\part*{Appendices}\\addcontentsline{toc}{part}{Appendices}\n\\markright{Appendices}\n"
    titles = dict(APPENDICES)
    out += latex.section("A\\qd %s" % titles["A"], "app:A")
    out += check_catalogue(ctx)
    out += latex.section("B\\qd %s" % titles["B"], "app:B")
    out += by_design_table(ctx)
    out += latex.section("C\\qd %s" % titles["C"], "app:C", needspace=12)
    out += relations_appendix(ctx)
    out += latex.section("D\\qd %s" % titles["D"], "app:D")
    out += glossary(ctx)
    out += latex.section("E\\qd %s" % titles["E"], "app:E")
    out += index(ctx)
    out += latex.section("F\\qd %s" % titles["F"], "app:F")
    out += appendix_lead(ctx, "appendix_F")
    out += "\\parhead{Environment}\n" + environment_table(ctx)
    out += appendix_rules(ctx)
    out += "\\parhead{Vocabulary}\n" + classes_tables()
    out += "\\parhead{The deploy gate}\n" + gate_problems(ctx)
    out += latex.section("G\\qd %s" % titles["G"], "app:G")
    out += appendix_lead(ctx, "appendix_G")
    out += mcp_limits(ctx)
    return out


def appendix_lead(ctx, key):
    """The lead-in paragraphs of an appendix from the overview prose (`appendix_F`, `appendix_G`), if any."""
    entry = ctx.overview.get(key) if isinstance(ctx.overview.get(key), dict) else {}
    items = entry.get("paragraphs") or []
    if isinstance(items, str):
        items = [items]
    return paragraphs([ctx.fill(p) for p in items]) if items else ""


def appendix_rules(ctx):
    """The rules of judging, the verdict and the gate (overview.json `appendix_F.rules`, with its `rules_title`):
    the long form the owner's Part I leaves to an appendix (methods: the arms, the order, the ground truth and
    the timing definitions), printed before the vocabulary; nothing when the prose has none."""
    entry = ctx.overview.get("appendix_F") if isinstance(ctx.overview.get("appendix_F"), dict) else {}
    items = entry.get("rules") or []
    if isinstance(items, str):
        items = [items]
    if not items:
        return ""
    return "\\parhead{%s}\n" % esc(entry.get("rules_title") or "Judging, the verdict and the gate") + paragraphs(
        [ctx.fill(p) for p in items])
