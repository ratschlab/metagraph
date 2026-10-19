"""The sections of the results document (DESIGN-DOC §1, §4, §6–§8, the owner's review of the sample): the cover
lines, Part I in its four sections (decision summary, methods, reading guide, provenance and coverage), the
chapter pages (the conclusion first), the family sections with their unit tables and detail blocks, the closing
section (observations, plots, tables) and the appendices. Every function returns LaTeX; the data come from a
`reader.Run`, a `units.Units` and the prose dicts.
"""

import hashlib
import json
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
import merge_run  # noqa: E402
import reader  # noqa: E402
import report  # noqa: E402
import suite_common as sc  # noqa: E402
import vocabulary  # noqa: E402
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
# the closing section's displayed prefix (READ-009: chapter C is capabilities and appendix C the relations, so the
# timing sections are not C.1 .. C.7 on the page; the LaTeX labels sec:C1 .. sec:C7 stay, so that links hold)
CLOSING_PREFIX = "Timing"
# the deploy gate's problems grouped by what they ask of the owner (READ-002, READ-023): the gate's kinds stay the
# record's words; the groups are the document's reading of them. A kind decides the group except where the problem's
# detail does (an outcome problem: inconclusive is behaviour not shown, a finding or a failure is observed; a
# relation: inconclusive or not evaluated against fail).
GATE_GROUPS = (("evidence", "release evidence missing or mismatched"),
               ("unproven", "required behaviour not shown"),
               ("violation", "observed violation of a rule the suite applies"))
GATE_GROUP_ACTION = {
    "evidence": ("the problem names a record that is missing, stale or not the one the policy pins: restore or "
                 "regenerate it (re-evaluate the set with the recorded checks, recompute the comparison, regenerate "
                 "the policy from the table the comparison read) or rerun with a runner that records the live "
                 "documents. These problems concern the release evidence, not the host's or the service's behaviour."),
    "unproven": ("the property is not shown on that access path: the execution was inconclusive, a required check "
                 "had no decisive result, a by-design pair was not verified on one side, or a required relation could "
                 "not be evaluated. What clears it is a rerun on which the property is shown (a 200 within its budget, "
                 "a walk reaching its minimum) or an owner-declared exemption with its reason; nothing here says the "
                 "behaviour is wrong."),
    "violation": ("a rule the suite applies was broken in this run: a finding (a stated bound broken on the clock the "
                  "check reads; for a 503 deadline delivered late, the suite's tolerance on the client's clock, where "
                  "the time went not isolated), a failed check, a disagreement between access paths, a failed "
                  "relation, an answer the manifest does not predict. What clears it: for a late 503, a rerun within "
                  "the tolerance, or a fix or a SPEC decision once the delay is located; for a failure or a "
                  "disagreement, a fix or a SPEC decision, then a rerun; the detail block of each test holds the "
                  "evidence."),
}
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
        self.cases_prose = cases or {}  # doc/prose/cases.json: {short: {case id: text}, verification: {case id: text}}
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
        self._resolve_stops = None  # facts.resolve_stops, computed once (placeholders() is built on every fill)

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

    def current_notes(self):
        """doc/prose/cases.json `current`: a case's note on how its pinned table statement reads against the
        deployment this run met (manifest-version conditionals), printed in appendix B under its short reason."""
        s = self.cases_prose.get("current")
        return s if isinstance(s, dict) else {}

    def verification_rules(self):
        """doc/prose/cases.json `verification`: the rule by which a case's pairs are verified or not, printed in
        appendix B under the case's short reason and nowhere else (the short reason stays one clause)."""
        s = self.cases_prose.get("verification")
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
        if getattr(self, "_resolve_stops", None) is None:
            self._resolve_stops = facts.resolve_stops(run)
        stops = self._resolve_stops
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
                "verdict_text": doc_verdict_text(run),
                "verdict_headline": report.verdict_headline(doc_verdict_text(run), run.gated, "beside") if run.verdict
                else (run.verdict_headline or "–"),
                "finding_kinds_text": finding_kinds_text(run),
                "n_reruns_for_resume": len(run.reruns),
                "n_resolve_stops": stops["n_resolve_stops"], "resolve_stop_observed": stops["resolve_stop_observed"],
                "resolve_stop_ids": stops["resolve_stop_ids"],
                "n_findings": sum((totals.get(a) or {}).get("outcomes", {}).get("finding", 0) for a in ARMS),
                "n_inconclusive": sum((totals.get(a) or {}).get("outcomes", {}).get("inconclusive", 0) for a in ARMS),
                "n_failures": sum((totals.get(a) or {}).get("outcomes", {}).get(k, 0) for a in ARMS
                                  for k in ("fail", "error")),
                "n_executions": "{:,}".format(sum((totals.get(a) or {}).get("executions", 0) for a in ARMS)),
                # the two counts of I.2's classifications in words, never "execution(s)" (RV12)
                "findings_in_run": _count_phrase(sum((totals.get(a) or {}).get("outcomes", {}).get("finding", 0)
                                                     for a in ARMS), "execution", "none in this run",
                                                 " in this run: %s" % finding_kinds_text(run)),
                "failures_in_run": _count_phrase(sum((totals.get(a) or {}).get("outcomes", {}).get(k, 0) for a in ARMS
                                                     for k in ("fail", "error")), "execution", "none"),
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
                "n_declared_exclusion_pairs": declared["exclusion_pairs"],
                "declared_exclusions_clause": declared["exclusions_clause"],
                "n_declared_not_required": declared["not_required"],
                "n_declared_relations": declared["relations_excluded"],
                "gate_word": gate_word if gate else "not evaluated",
                "gate_text": gate_text(gate_word, gate_how, gate_counts, gate_problems_text, run),
                "n_gate_outcome_problems": gp["outcome"], "n_gate_inconclusive_executions": gp["inconclusive"],
                "n_gate_finding_executions": gp["finding"], "n_gate_outcome_pairs": gp["pairs"],
                "gate_inconclusive_clause": gate_inconclusive_clause(run),
                "untested_scope": untested_scope_text(self),
                "by_design_cases": table.get("cases") if table.get("cases") is not None else "–",
                "by_design_sha_short": (table.get("sha256") or "–")[:12],
                "by_design_table_line": by_design_table_text(run),
                "n_chapters": len(self.units.chapters), "built_chapters": ", ".join(self.chapters),
                "readme_path": README_PATH, "readme_commit": self.readme_commit() or "not determined",
                "logan_caveat": facts.logan_caveat(run), "gate_judged": bool(gate),
                # the readability review (2026-10-15): the gate's scope and its problems by group (READ-002,
                # READ-023), the counting units (READ-006), the standing examples (READ-007, READ-102), the
                # service's revision (READ-021), the relation examples (READ-019) and the shared vocabulary
                # (READ-022, filled through `vocab[...]`)
                "gate_routes": gate_scope(run)["routes"], "gate_arms": gate_scope(run)["arms"],
                "n_gate_selected": gate_scope(run)["selected"], "gate_scope_text": gate_scope(run)["text"],
                "gate_decision": gate_decision_short(run), "gate_groups_text": gate_groups_text(run),
                # the external review of the v6 document (V6-001): the policy's automatic exclusions and the
                # required population per gated access path
                "gate_auto_text": gate_auto_text(run) or "none (no policy)",
                "gate_required_text": gate_required_text(run),
                "n_gate_problems": len([x for x in (gate.get("problems") or []) if isinstance(x, dict)]),
                "n_gate_evidence": gate_group_counts(run)["evidence"],
                "n_gate_unproven": gate_group_counts(run)["unproven"],
                "n_gate_violations": gate_group_counts(run)["violation"],
                "gate_read_text": gate_read_text(run),
                "n_property_claims": self.untested()["by_chapter"] and sum(
                    d["pairs"] for d in self.untested()["by_chapter"].values()) or 0,
                "n_property_claims_untested": sum(d["untested"] for d in self.untested()["by_chapter"].values()),
                "counting_example": counting_example_text(self),
                "standing_example": standing_example_text(self),
                "standing_vs_outcome_example": standing_vs_outcome_example_text(self),
                "standing_vs_outcome_example_bd": standing_vs_outcome_example_text(self, kind="by_design"),
                "service_version_line": service_version_text(run),
                "relation_example_subset": relation_example_text(run, "subset"),
                "relation_example_equivalence": relation_example_text(run, "same_answer"),
                "vocab": vocabulary.glossary_texts(),
                # an answer-side scenario-level mark in the manifest (the delta verifier's M2): the glossary's and
                # appendix G's definitions of not implemented via MCP name it
                "answer_side_mark": bool(answer_side_mark_ids(run)),
                # the MCP tools the scenarios call, counted from the manifest (RV14: what the suite covers via MCP)
                "mcp_tools_text": mcp_tools_text(run),
                # the MCP runner's git head, short, for the prose that credits its normaliser (V12)
                "mcp_runner_head": str((run.mcp_run or {}).get("git_head") or "–")[:8],
                # whether this document reads a re-evaluated set: a sentence about the re-evaluation prints only then
                # (the reader review of the v7 document, RV6; the fact verifier, V12)
                "reevaluated": bool(run.reevaluated),
                # the short names of the gate policies reported, for the glossary (RV5)
                "gate_policy_names": gate_policy_names(run),
                # the kinds of property-level mark the manifest carries, for appendix G's lead-in (V10)
                "property_mark_kinds": "; ".join(w for k, w in (("not_implemented", "not implemented via MCP"),
                                                                ("by_design", "service-side by design"))
                                                 if _has_property_marks(run, k)),
                "answer_side_mark_text": answer_side_mark_text(run),
                # the annotated matrix's column note (annotated_matrix fills both from the unit it picks)
                "annotated_dimension": "the varied value", "annotated_record_key": False}

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


def _stats_cells(d, signed=False, tint=True):
    """p50 / p90 / max / n of a summarise() dict; `signed` for a difference (an overhead, a net cost, a delta),
    whose figures carry their sign and no tint; `tint=False` for a time that is no execution time (a connection
    set-up)."""
    if not isinstance(d, dict) or not d.get("n"):
        return ["–", "–", "–", "0"]
    f = fmt_delta if signed else (cells.fmt_exec if tint else fmt_t)
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
    if not commit:
        return "not recorded"
    text = "%s (the runners' scripts, as recorded" % commit
    # the MCP runner is another program at its own git head; it judged with the run's snapshot of the checks, the
    # same checks.py as the direct runners when the recorded sha256 agree (the reader review of the v7 document, RV6)
    mcp = arms.get("mcp") or {}
    head = (getattr(run, "mcp_run", None) or {}).get("git_head")
    shas = set(str((info or {}).get("checks_sha256")) for info in arms.values() if (info or {}).get("checks_sha256"))
    if head and len(shas) == 1 and mcp.get("checks_sha256"):
        text += "; the MCP runner, a program at its own git head %s%s, judged with the same checks.py, sha256 %s…" % (
            str(head)[:12], mcp_runner_repository_text(run), list(shas)[0][:12])
    return text + ")"


def mcp_tools_text(run):
    """'four tools: pattern_search (456 scenarios), traverse_resolve (118), traverse_fetch (99) and catalog (1)': the
    MCP tools the manifest's scenarios name (mcp.tool), by count; '' when none names one."""
    counts = {}
    for sid in getattr(run, "order", None) or []:
        tool = ((run.scenario(sid) or {}).get("mcp") or {}).get("tool")
        if tool:
            counts[tool] = counts.get(tool, 0) + 1
    if not counts:
        return ""
    items = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    parts = ["%s (%d scenario%s)" % (items[0][0], items[0][1], "" if items[0][1] == 1 else "s")]
    parts += ["%s (%d)" % kv for kv in items[1:]]
    return "%s tool%s: %s" % (facts.number_word(len(items)), "" if len(items) == 1 else "s", facts._join(parts))


def mcp_runner_repository_text(run):
    """', in <repository>' of the MCP runner, as the identity input states it (make_doc --identity:
    `mcp_runner_repository`; the run records the runner's git head, not its repository), '' otherwise."""
    ident = getattr(run, "identity", None)
    repo = ident.get("mcp_runner_repository") if isinstance(ident, dict) else None
    return (", in %s" % repo) if repo else ""


def mcp_runner_text(run):
    """The MCP access path's runner for the tables (RV6: never '–' beside the direct runners' hosts): the MCP client
    at its git head (and repository when given), else '–'."""
    head = (getattr(run, "mcp_run", None) or {}).get("git_head")
    return ("the MCP client, git head %s%s" % (str(head)[:12], mcp_runner_repository_text(run))) if head else "–"


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
            elif reason == getattr(sc, "SCENARIO_LEVEL_MARK", None):
                # the manifest's marks cover every property check of the execution (the v6 delta: D-P-scope-0004's
                # MCP r0, marks_cover_repeat): the mark's reading, never a host state
                words.append("%d by a scenario-level mark (every property check marked as not shown via MCP)" % n)
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
    dec = getattr(run, "declared", None)
    if not isinstance(dec, dict):
        return {"file": "none", "version": "–", "exclusions": 0, "not_required": 0, "relations_excluded": 0,
                "exclusion_pairs": 0, "exclusion_repeat_only_pairs": 0, "exclusions_clause": "no exclusion"}
    out = {"file": dec.get("_file") or "declared", "version": dec.get("declared_version") or "–"}
    for k in ("exclusions", "not_required", "relations_excluded"):
        out[k] = len([e for e in (dec.get(k) or []) if isinstance(e, dict)])
    # an exclusion entry names one scenario on one or several access paths, whole or for some repeats: the
    # scenario-arm pairs it covers are the gate's unit (READ-006; the verification of 2026-10-16, shared M6)
    pairs, repeat_only = 0, []
    for e in dec.get("exclusions") or []:
        if not isinstance(e, dict):
            continue
        n = len(e.get("arms") or []) or 1
        pairs += n
        if e.get("repeats"):
            repeat_only.append((n, ", ".join(str(r) for r in e["repeats"])))
    out["exclusion_pairs"] = pairs
    out["exclusion_repeat_only_pairs"] = sum(n for n, _ in repeat_only)
    reps = sorted(set(r for _, r in repeat_only))
    out["exclusions_clause"] = "%d exclusion%s (%d scenario-arm pair%s%s)" % (
        out["exclusions"], "" if out["exclusions"] == 1 else "s", pairs, "" if pairs == 1 else "s",
        (", %d of them for repeat%s %s only" % (out["exclusion_repeat_only_pairs"], "" if len(reps) == 1 and "," not in reps[0]
                                                 else "s", " / ".join(reps))) if repeat_only else "")
    return out


def gate_text(word, how, counts_text, problems_text, run=None):
    """The deploy gate in one plain line (the cover's and I.1's wording): the verdict word, the plain decision
    (READ-002: what the gate covers, read from the policy, and that its problem count is not a bug count), the
    problems by group and by kind, the executions judged, and how it was judged."""
    if counts_text is None:
        return "%s (%s)" % (word, how)
    scope = gate_scope(run)["text"] if run is not None else "the policy's selection"
    passed = run is not None and bool((run.gate or {}).get("passed"))
    if passed:
        return "%s: nothing blocks (%s); %s; the gate covers %s; %s" % (
            word, problems_text, counts_text, scope, _without_access_paths(how, run))
    groups = gate_groups_text(run) if run is not None else "–"
    head, kinds = _problems_head_and_kinds(problems_text)
    return ("%s: the run does not authorise a deploy. %s, by what would clear them: %s (by the gate's kinds: %s); the "
            "count is not a count of service defects; %s; the gate covers %s; %s" % (
                word, head[:1].upper() + head[1:], groups, kinds, counts_text, scope,
                _without_access_paths(how, run) if run is not None else how))


def _problems_head_and_kinds(problems_text):
    """('50 problems', '5 identity; 13 configuration; ...') of gate_summary's problems text, so that a line prints the
    kinds once (the readability review's verification of 2026-10-16, shared L5)."""
    if " (" in problems_text:
        head, _, rest = problems_text.partition(" (")
        return head, rest.rstrip(")")
    return problems_text, "–"


def _without_access_paths(how, run):
    """gate_summary's `how` without its trailing access-paths clause where the line names the gate's scope with the
    access paths already (gate_scope: a policy pins the run); unchanged when the scope cannot name them."""
    if gate_scope(run)["arms"] == "–":
        return how
    return re.sub(r", access paths [^;]+$", "", how)


def gate_outcome_counts(run):
    """The gate's outcome problems (gate.json, else the offline verdict): {judged, outcome: all of them,
    inconclusive: the inconclusive executions among them, finding: the findings, pairs: distinct (scenario, arm)}.
    These, not the run's inconclusive count, are what blocks the gate: an execution on a gated arm, in the
    policy's required population: selected, not removed by its automatic exclusions nor declared exempt (VERIFY-R2 H1, V6-001). Zeros with judged False when no gate was judged."""
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
    return ("%d of them block%s the deploy gate: %s the policy requires on that access path and repeat (the gate's "
            "outcome problems, appendix F); the others are on an ungated access path, outside the selection, removed "
            "by the policy's automatic exclusions or declared exempt" % (n, "s" if n == 1 else "",
                                                                          "the one" if n == 1 else "those"))


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


def gate_scope(run):
    """What the deploy gate covers, read from the policy, never from its file name (READ-002): {routes, arms,
    selected, text}; dashes without a policy."""
    pol = run.policy if isinstance(run.policy, dict) else {}
    sel = pol.get("selection") if isinstance(pol.get("selection"), dict) else {}
    routes = [str(r) for r in (sel.get("routes") or [])]
    arms = [str(a) for a in (pol.get("arms") or [])]
    selected = len(sel.get("scenarios") or [])
    if not pol:
        return {"routes": "–", "arms": "–", "selected": 0, "text": "no policy pins this run's manifest"}
    text = "the scenarios of %s (%d of the run's %d) on the access paths %s" % (
        facts._join(routes) if routes else "the policy's selection", selected, len(run.order),
        " and ".join(arms) if arms else "–")
    return {"routes": ", ".join(routes) or "–", "arms": ", ".join(arms) or "–", "selected": selected, "text": text}


# The policy's automatic exclusion rules (gate_policy.py), in plain words and in the ledger's order (V6-001): derived
# from the manifest, its tags and the by-design table, not decisions of the owner; `declared` (the owner's file) is
# counted apart. A rule the generator does not know prints its own name.
GATE_RULE_WORDS = (
    ("budget_stress", "tagged budget_stress: a property that depends on the host's load, which no budget makes "
                      "decisive"),
    ("manifest_not_implemented", "the manifest marks the scenario not implemented via MCP (its request, or the field "
                                 "its check reads)"),
    ("table_not_expressible", "a by-design case the scenario names states the request not expressible via MCP"),
    ("advanced_path", "traverse_fetch arguments on the advanced strategy path, which the deployment does not enable"),
    ("marks_cover_repeat", "on this MCP repeat the manifest's property-level marks set aside every property check "
                           "it can test"),
    ("needs_direct_arm", "a by-design case that only a direct execution verifies, in a policy without a direct "
                         "access path"),
)
GATE_MCP_EXPRESS_RULES = ("manifest_not_implemented", "table_not_expressible")


def gate_population(run):
    """The deploy gate's populations from the policy (V6-001): the route selection, the policy's automatic
    exclusions grouped by rules / access paths / repeats, the owner's declared exclusions, and per gated access path
    the scenarios and executions the policy requires. None without a policy."""
    pol = run.policy if isinstance(run.policy, dict) else None
    if not pol:
        return None
    rank = dict((r, i) for i, (r, _) in enumerate(GATE_RULE_WORDS))
    groups, auto, declared = {}, [], []
    for ex in pol.get("exclusions") or []:
        if not isinstance(ex, dict):
            continue
        rules = tuple(str(r) for r in (ex.get("rules") or []))
        if "declared" in rules:
            declared.append(ex)
            continue
        auto.append(ex)
        key = (rules, tuple(ex.get("arms") or []), tuple(ex.get("repeats") or []))
        groups.setdefault(key, []).append(str(ex.get("scenario")))
    order = sorted(groups, key=lambda k: (min(rank.get(r, len(rank)) for r in k[0]) if k[0] else len(rank),
                                          -len(groups[k]), k))
    need = gate_mod.required_keys(pol)
    required = []
    for arm in pol.get("arms") or []:
        keys = [k for k in need if k[1] == arm]
        required.append((arm, len(set(k[0] for k in keys)), len(keys)))
    return {"selected": len((pol.get("selection") or {}).get("scenarios") or []),
            "n_run": len(getattr(run, "order", None) or []),
            "auto": auto, "groups": [(k, groups[k]) for k in order], "declared": declared,
            "required": required, "n_required": len(need)}


def _gate_auto_kind(rules):
    """The summary's kind of an automatic exclusion: budget, mcp (the tools cannot express it), advanced, marks,
    direct, or the first rule's own name."""
    if "budget_stress" in rules:
        return "budget"
    if any(r in GATE_MCP_EXPRESS_RULES for r in rules):
        return "mcp"
    return {"advanced_path": "advanced", "marks_cover_repeat": "marks", "needs_direct_arm": "direct"}.get(
        rules[0] if rules else "", rules[0] if rules else "?")


def gate_auto_text(run):
    """The policy's automatic exclusions in one plain clause (V6-001): '94 entries: 35 budget-stress tests on
    direct-local and mcp, 57 requests the MCP tools cannot express on mcp, ...'; '' without a policy."""
    pop = gate_population(run)
    if pop is None:
        return ""
    kinds = {}
    for (rules, arms, repeats), sids in pop["groups"]:
        k = kinds.setdefault(_gate_auto_kind(rules), {"n": 0, "arms": [], "repeats": set(), "sids": []})
        k["n"] += len(sids)
        k["sids"].extend(sids)
        for a in arms:
            if a not in k["arms"]:
                k["arms"].append(a)
        k["repeats"].update(repeats)
    answer_side = set(answer_side_mark_ids(run)) if hasattr(run, "scenario") and hasattr(run, "order") else set()
    words = {"budget": ("budget-stress test", "budget-stress tests"),
             "mcp": ("request the MCP tools cannot express", "requests the MCP tools cannot express"),
             "advanced": ("request on the advanced strategy path", "requests on the advanced strategy path"),
             "marks": ("repeat whose property checks the marks all set aside",
                       "repeats whose property checks the marks all set aside"),
             "direct": ("by-design case only a direct execution verifies",
                        "by-design cases only a direct execution verifies")}
    parts = []
    for kind in ("budget", "mcp", "advanced", "marks", "direct") + tuple(sorted(k for k in kinds if k not in words)):
        if kind not in kinds:
            continue
        k = kinds[kind]
        w = words.get(kind, (kind, kind))
        arms = [a for a in sc.ARMS if a in k["arms"]] + [a for a in k["arms"] if a not in sc.ARMS]
        # the verification of the review's fixes, M1: a scenario-level answer-side mark (v6's P-scope-0004 / -0006)
        # is expressible through the tools, its answer lacking the field its check reads: the group is worded by its
        # rules and those scenarios are named apart, never counted among the requests the tools cannot express
        answer = [s for s in k["sids"] if s in answer_side] if kind == "mcp" else []
        if answer:
            n_req = k["n"] - len(answer)
            parts.append("%d that the manifest or the by-design table marks not implemented or not expressible via "
                         "MCP, on %s (%d request%s the tools cannot express; %s, whose MCP answer lacks the field "
                         "%s check%s read%s)" % (
                             k["n"], " and ".join(arms) or "every access path", n_req, "" if n_req == 1 else "s",
                             facts._ids_short(answer), "its" if len(answer) == 1 else "their",
                             "" if len(answer) == 1 else "s", "s" if len(answer) == 1 else ""))
            continue
        parts.append("%d %s%s on %s%s" % (
            k["n"], w[0] if k["n"] == 1 else w[1],
            # the term glossed where a reader of I.1 first meets it (the reader review of the v7 document, RV9)
            " (a test whose outcome depends on the host's load: excluded from the gate by rule, its outcome still "
            "reported)" if kind == "budget" and k["n"] == 1 else
            " (tests whose outcome depends on the host's load: excluded from the gate by rule, their outcomes still "
            "reported)" if kind == "budget" else "",
            " and ".join(arms) or "every access path",
            (" (repeat%s %s only)" % ("" if len(k["repeats"]) == 1 else "s", ", ".join(str(r) for r in
                                                                         sorted(k["repeats"]))))
            if k["repeats"] else ""))
    n = len(pop["auto"])
    return "%d entr%s%s%s" % (n, "y" if n == 1 else "ies", (", each one test left out on the access paths its rule "
                                                           "names" if n > 1 else ", one test left out on the access "
                                                           "paths its rule names") if n else "",
                              (": " + facts._join(parts)) if parts else "")


def gate_required_text(run):
    """'686 scenarios (793 executions) on direct-local and 628 scenarios (1275 executions) on mcp' (V6-001)."""
    pop = gate_population(run)
    if pop is None:
        return "–"
    return facts._join(["%d scenario%s (%d execution%s) on %s" % (s, "" if s == 1 else "s", e, "" if e == 1 else "s",
                                                                 arm) for arm, s, e in pop["required"]])


def gate_population_text(run):
    """The sentence that defines the gate's required population (V6-001): route selection, the policy's automatic
    exclusions (appendix F), the owner's declared exemptions (I.4, appendix G), the required scenarios and executions
    per gated access path."""
    pop = gate_population(run)
    if pop is None:
        return ""
    scope = gate_scope(run)
    dec = declared_counts(run)
    return ("Required means selected by route (%s), not removed by one of the policy's automatic exclusions, which "
            "its rules derive from the manifest, its tags and the by-design table and which are no decision of the "
            "owner (%s; appendix F lists them), and not declared exempt by the owner (%s; I.4, appendix G): %s, %d "
            "executions in all. The gate certifies those executions, not the excluded ones and not the whole "
            "service." % (scope["text"], gate_auto_text(run), ("%d exclusion%s" % (
                dec["exclusions"], "" if dec["exclusions"] == 1 else "s")) if isinstance(getattr(
                    run, "declared", None), dict)
                else "no declared file", gate_required_text(run), pop["n_required"]))


def gate_judged_outcomes(run, policy=None):
    """{outcome: executions} over what the gate judged (V6-001): the policy's selection on its access paths less
    its exclusions, as gate.judge counts `judged`; {} without a policy. |policy|: another than the deploy gate's."""
    pol = policy if isinstance(policy, dict) else (run.policy if isinstance(run.policy, dict) else None)
    if not pol:
        return {}
    selected = set((pol.get("selection") or {}).get("scenarios") or [])
    arms = set(pol.get("arms") or [])
    out = {}
    for r in getattr(run, "rows", None) or []:
        sid, arm, rep = sc.execution_key(r)
        if sid not in selected or arm not in arms or gate_mod.excluded(pol, sid, arm, rep):
            continue
        o = str(r.get("outcome"))
        out[o] = out.get(o, 0) + 1
    return out


def gate_ledger(ctx):
    """Appendix F's ledger of the policy's automatic exclusions (V6-001): one row per rule combination, access paths
    and repeats, with the scenarios; the owner's declared entries stay in appendix G."""
    pop = gate_population(ctx.run)
    if pop is None or not pop["groups"]:
        return ""
    words = dict(GATE_RULE_WORDS)
    rows = []
    for (rules, arms, repeats), sids in pop["groups"]:
        why = "; and ".join(words.get(r, r) for r in rules) or "–"
        rows.append([esc(why) + " (%s)" % ", ".join(tt(r) for r in rules),
                     esc(", ".join(arms) or "all") + ((" (repeat%s %s)" % ("" if len(repeats) == 1 else "s",
                                                                          ", ".join(str(r) for r in repeats)))
                                                      if repeats else ""),
                     ", ".join(_link_sid(ctx, s) for s in sids), str(len(sids))])
    out = "\\parhead{The policy's automatic exclusions}\n"
    out += ("%s The rules derive each exclusion from the manifest, its tags and the by-design table, never from an "
            "answer the run received; they are not the owner's decisions, which appendix G lists apart. An entry "
            "removes one scenario on the access paths named, or only the repeats named, from the gate's required "
            "executions; the scenario keeps its recorded outcome in its chapter, and nothing it did not show is "
            "shown by the exclusion. A scenario can be excluded on one access path and required on the other, so "
            "the counts below are entries, not distinct tests.\n\n" % esc(
                "The policy %s excludes %s." % (os.path.basename(ctx.run.policy_path or "") or "–",
                                                gate_auto_text(ctx.run))))
    out += latex.plain_table(latex.spec_p([52, 22, 80, 12]), ["rule", "access paths", "scenarios", "entries"],
                             rows, needspace=min(len(rows), 4) + 3)
    return out


def gate_group_of(problem):
    """The group (GATE_GROUPS) of one gate problem, from its kind and, where the kind does not decide, its detail."""
    kind = str(problem.get("kind") or "")
    detail = str(problem.get("detail") or "").strip()
    if kind in ("identity", "configuration", "missing", "agreement"):
        return "evidence"
    if kind in ("property", "nothing_passed"):
        return "unproven"
    if kind == "disagreement":
        return "violation"
    if kind == "relation":
        word = detail.split(": ", 1)[1].split(":", 1)[0].strip() if ": " in detail else detail
        return "violation" if word.startswith("fail") else "unproven"
    if kind == "outcome":
        if detail.startswith("inconclusive") or detail.startswith("by-design"):
            return "unproven"
        return "violation"
    return "violation"


def gate_group_counts(run):
    gate = run.gate if isinstance(run.gate, dict) else {}
    counts = dict((g, 0) for g, _ in GATE_GROUPS)
    for p in gate.get("problems") or []:
        if isinstance(p, dict):
            counts[gate_group_of(p)] += 1
    return counts


def gate_groups_text(run):
    """'5 release evidence problems, 44 required behaviours not shown, 1 observed violation', from the gate's
    problems (plain text); 'none' without a problem."""
    counts = gate_group_counts(run)
    words = {"evidence": ("release evidence problem", "release evidence problems"),
             "unproven": ("required behaviour not shown", "required behaviours not shown"),
             "violation": ("observed violation", "observed violations")}
    parts = ["%d %s" % (counts[g], words[g][0 if counts[g] == 1 else 1]) for g, _ in GATE_GROUPS if counts[g]]
    return ", ".join(parts) if parts else "none"


def gate_decision_short(run):
    """I.1's release decision without the population sentence, which I.1's population table carries (the second
    external review, R2-003); the long form (gate_decision_text) stays in appendix F."""
    gate = run.gate if isinstance(run.gate, dict) else None
    if gate and gate_population(run) is not None:
        if gate.get("passed"):
            return ("The gate passed: every planned execution was recorded, every execution the policy requires passed "
                    "or matched the by-design case its scenario names, and nothing blocks. It certifies those required "
                    "executions (the table below), not the excluded ones and not the whole service.")
        # the failure in the same short form, the populations in the table (the verification of round 2, R4)
        n = len([x for x in (gate.get("problems") or []) if isinstance(x, dict)])
        return ("This run does not authorise a deploy: the gate failed on %d problem%s over the executions the policy "
                "requires (the table below). The count is not a count of service defects: %s (the groups in appendix "
                "F, each with the evidence that would clear it)." % (n, "" if n == 1 else "s", gate_groups_text(run)))
    return gate_decision_text(run)


def gate_population_table(ctx):
    """I.1's table of the gate's populations (R2-003): the route selection, the policy's automatic exclusions, the
    owner's exemptions, the required scenarios and executions per gated access path and in all, with the judged
    executions by outcome; empty without a policy."""
    run = ctx.run
    pop = gate_population(run)
    gate = run.gate if isinstance(run.gate, dict) else None
    if pop is None or not gate:
        return ""
    scope = gate_scope(run)
    dec = declared_counts(run)
    rows = [["selected by route", esc("%d of the run's %d scenarios (%s)" % (pop["selected"], pop["n_run"],
                                                                         scope["routes"]))],
            ["the policy's automatic exclusions", esc("%s; rules derived from the manifest, its tags and the "
                                                      "by-design table, no decision of the owner (appendix F)"
                                                      % gate_auto_text(run))],
            ["the owner's exemptions", esc(("%d exclusion%s, %d check%s not required, %d relation%s excluded "
                                            "(appendix G)" % (dec["exclusions"], "" if dec["exclusions"] == 1 else "s",
                                                              dec["not_required"], "" if dec["not_required"] == 1
                                                              else "s", dec["relations_excluded"], "" if
                                                              dec["relations_excluded"] == 1 else "s"))
                                           if isinstance(getattr(run, "declared", None), dict) else "no declared file")]]
    for arm, s_, e_ in pop["required"]:
        rows.append(["required on %s" % esc(arm), "%d scenarios, %d executions (repeats included)" % (s_, e_)])
    counts = gate.get("counts") or {}
    rows.append(["required in all", esc("%d executions%s" % (pop["n_required"],
                                                           judged_outcomes_text(run, counts).replace(",", ":", 1)
                                                           if judged_outcomes_text(run, counts).startswith(",")
                                                           else judged_outcomes_text(run, counts)))])
    ungated = [a for a in ARMS if a not in ((run.policy or {}).get("arms") or [])]
    note = (" %s is reported but not gated by this policy." % facts._join(ungated)) if ungated else ""
    return latex.plain_table(latex.spec_p([46, 128]), ["population", "this run"], rows,
                             "the deploy gate's populations: exclusion entries are not distinct tests (an entry can "
                             "apply to one access path or one repeat), and scenario counts are not added across "
                             "access paths;%s" % esc(note), needspace=len(rows) + 4)


def gate_decision_text(run):
    """The release decision in one plain sentence (READ-002): what the gate's verdict means for a deploy, and
    what its problem count is not."""
    gate = run.gate if isinstance(run.gate, dict) else None
    if not gate:
        return ("No deploy decision was judged for this run: no gate verdict is recorded and no policy pins its "
                "manifest, so this document authorises nothing.")
    scope = gate_scope(run)
    population = gate_population_text(run)
    if gate.get("passed"):
        if population:
            return ("The gate passed: every planned execution was recorded, every execution the policy requires "
                    "passed or matched the by-design case its scenario names, and nothing blocks. %s" % population)
        return ("The gate passed: every planned and required execution of %s was recorded and nothing blocks. The "
                "gate certifies that selection on those access paths, not the whole service." % scope["text"])
    n = len([x for x in (gate.get("problems") or []) if isinstance(x, dict)])
    return ("This run does not authorise a deploy: the gate failed on %d problem%s over %s. The count is not a count "
            "of service defects: %s (the groups in appendix F, each with the evidence that would clear it).%s" % (
                n, "" if n == 1 else "s", scope["text"], gate_groups_text(run), (" " + population) if population
                else ""))


def gate_read_text(run):
    """What the gate judged (READ-023): the set it read, the comparison, the policy, and whether it was judged
    offline for this document or recorded in the run."""
    gate = run.gate if isinstance(run.gate, dict) else None
    if not gate:
        return "no gate verdict"
    pol_path = os.path.basename(run.policy_path or "") or "the recorded policy"
    where = ("judged offline by gate.py when this document was built" if run.gate_offline
             else "recorded in the run by gate.py")
    return ("%s, from the run directory's %s set (%s and %s) and the policy %s; the checks it read for the property "
            "reading are those of the checkout at build time, and an identity problem that names that checkout's "
            "checks.py or suite_common.py compares it with the checks the %s recorded" % (
                where, "re-evaluated" if run.reevaluated else "recorded", run.executions_file,
                sc.reevaluated_name("compare.json", run.reevaluated), pol_path,
                "re-evaluation" if run.reevaluated else "runners"))


def counting_example_text(ctx):
    """One scenario of the run as a worked example of the counting units (READ-006): its executions per access
    path and its property claims; plain text."""
    run = ctx.run
    best = None
    for sid in run.order:
        n = sum(len(run.arm_lines(sid, a)) for a in ARMS)
        if all(run.arm_lines(sid, a) for a in ARMS) and (best is None or n > best[1]):
            best = (sid, n)
            if n >= 5:
                break
    if not best:
        return "no scenario of this run was executed on every access path"
    sid, n = best
    per = ", ".join("%d on %s" % (len(run.arm_lines(sid, a)), a) for a in ARMS)
    props = sum(1 for e in run.scenario(sid).get("expect") or []
                if isinstance(e, dict) and e.get("check") and run.is_property(e.get("check")))
    return ("%s is one scenario; it was executed %d times (%s), so it counts %d executions, and it names %d property "
            "check%s, so it counts %d property claim%s" % (sid, n, per, n, props, "" if props == 1 else "s", props,
                                                            "" if props == 1 else "s"))


def standing_example_text(ctx):
    """A real MCP cell of the run that shows a documented refusal, not a tested property (READ-007): the first
    scenario standing service-side by design via MCP with a verified pair; plain text."""
    run = ctx.run
    for sid in run.order:
        st = ctx.state(sid)["mcp"]
        if not (st.present and st.glyph_key == "service_by_design" and "u" not in st.glyph_marks):
            continue
        lines = blocks.class_lines(st)
        first = next((l for l in lines if l.get("http") is not None), None)
        code = (first.get("error") or {}).get("code") if first and isinstance(first.get("error"), dict) else None
        status = ("http %s%s" % (first.get("http"), (" " + code) if code else "")) if first else "no status"
        case = st.agreement.get("case") or "its case"
        dl = ctx.state(sid)["direct-local"]
        dl_word = facts.STANDING_WORD.get(dl.glyph_key, "no execution") if dl.present else "no execution"
        return ("%s: via MCP %s (%s, case %s), the recorded outcome is service_by_design and the cell shows the "
                "by-design glyph without a footnote; direct-local stands %s and matched the case's host side, so the "
                "pair is verified. The cell establishes the documented service behaviour, not the property: nothing "
                "of the claim was tested through the tools" % (sid, service_decision_words(run, st), status, case, dl_word))
    return "no by-design pair in this run"


SERVICE_DECISION_WORDS = {
    "refused_at_submit": "the service refused the call at submit, in the host's place",
    "relayed": "the service relayed the host's refusal under its own code",
    "answered": "the service answered in the host's place",
    "not_expressible": "the service could not express the call",
}


def service_decision_words(run, st):
    """How the service decided an MCP call, in words from the case's `service_answer.outcome` as the comparison
    recorded it (mcp_agreement.service_outcome, else the table's case; the same evidence blocks.mcp_class reads:
    refused at submit, relayed, answered, not expressible), never a fixed phrase (the readability review's
    verification of 2026-10-16, shared M4: R-edge-0002's 502 is the host's 400 relayed)."""
    outcome = st.agreement.get("service_outcome")
    if not outcome:
        case_id = st.agreement.get("case") or next((l.get("service_by_design_case") for l in st.lines
                                                     if l.get("service_by_design_case")), None)
        case = (getattr(run, "cases", None) or {}).get(case_id) or (getattr(run, "table", None) or {}).get(case_id) or {}
        outcome = (case.get("service_answer") or {}).get("outcome")
    return SERVICE_DECISION_WORDS.get(str(outcome or ""), "the service decided the call in the host's place")


_ROW_WORD = {"service_by_design": "by design", "not_implemented_mcp": "not implemented", "mixed": "mixed"}


def standing_outcome_differences(ctx, sids=None):
    """[(sid, kind, text)] of the tests whose MCP standing (the cell's glyph: the agreement class where the service
    decided the call) differs from their outcome in the chapter's outcomes row (compare's per-test outcome over the
    repeats, compare.aggregate_outcome), of the given scenarios (a chapter's) or of the run, in their order (READ-102,
    the re-verification V2-M1). Two kinds: 'not_implemented', a refusal the by-design table files under not
    implemented that stands not exercised while its line records service_by_design; 'by_design', a test one repeat
    of which the service decided by design while the host answered the other, which stands by design while the
    outcomes row counts it by its repeats (mixed, inconclusive or the other repeat's outcome); 'other' for any
    other difference. The direct access paths cannot differ (their standing is the worst outcome, which the row's
    rule gives too). Plain text."""
    run = ctx.run
    out = []
    for sid in (run.order if sids is None else sids):
        st = ctx.state(sid)["mcp"]
        if not st.present:
            continue
        outcome = (run.comp_arm(sid, "mcp") or {}).get("outcome")
        standing = st.glyph_key if st.glyph_key != "other" else outcome
        if outcome is None or standing == outcome:
            continue
        reps = ", ".join("r%s %s" % (l.get("repeat"), cells.WORD.get(l.get("outcome"), str(l.get("outcome"))))
                         for l in st.lines)
        row = _ROW_WORD.get(outcome, str(outcome))
        if st.glyph_key == "not_implemented_mcp":
            case = st.agreement.get("case") or next((l.get("service_by_design_case") for l in st.lines
                                                      if l.get("service_by_design_case")), None)
            text = ("%s: the MCP line records service_by_design (the service refused at submit) while the table files "
                    "the case%s under %s, so the cell and the standing count read %s, the outcomes row "
                    "%s" % (sid, (" %s" % case) if case else "", cells.attn_marked("not implemented via MCP"),
                            cells.attn_marked("not exercised"), row))
            out.append((sid, "not_implemented", text))
        elif st.glyph_key == "service_by_design":
            decided = [l for l in st.lines if cells.by_design_line(st, l)]
            out.append((sid, "by_design", "%s: via MCP %s, the cell standing by design (the case covers %s) while "
                                          "the outcomes row counts it %s" % (
                                              sid, reps, facts._join("r%s" % l.get("repeat") for l in decided)
                                              or "no repeat", row)))
        else:
            out.append((sid, "other", "%s: %s via MCP, the cell standing %s, the outcomes row %s" % (
                sid, reps, facts.STANDING_WORD.get(st.glyph_key, st.glyph_key), row)))
    return out


def standing_vs_outcome_example_text(ctx, sids=None, kind="not_implemented"):
    """The first test of the given scenarios (a chapter's) or of the run whose MCP standing and outcome differ by
    `kind` (standing_outcome_differences), plain text; 'none in this run' / 'none in this chapter' when none."""
    for sid, k, text in standing_outcome_differences(ctx, sids):
        if k == kind:
            return text
    return "none in this run" if sids is None else "none in this chapter"


def standing_vs_outcome_reading(ctx, sids):
    """The reading line's clause on a chapter's differing tests (V2-M1): 'none in this chapter' exactly when no
    test's MCP standing differs from its outcome (the two rows then agree per access path), else the count and one
    named case per kind, with how many more of that kind."""
    diffs = standing_outcome_differences(ctx, sids)
    if not diffs:
        return "none in this chapter"
    parts = []
    for kind in ("not_implemented", "by_design", "other"):
        items = [t for _, k, t in diffs if k == kind]
        if items:
            parts.append(items[0] + ((" (and %d more of this kind)" % (len(items) - 1)) if len(items) > 1 else ""))
    return "in this chapter %d test%s: %s" % (len(diffs), "" if len(diffs) == 1 else "s", "; ".join(parts))


SERVICE_VERSION_KEYS = ("service_version", "service_commit", "service_git_head", "service_revision", "deployed_version")


def _commit_text(commit, full):
    commit = str(commit)
    return commit if full or len(commit) <= 12 else commit[:12] + "…"


def recorded_service_identity(run):
    """The service identity the MCP runner recorded (merge_run.service_identity: the run's `service_identity`, else
    summarised from the MCP arm's run.json), None when the run has no MCP arm."""
    ident = getattr(run, "service_identity", None)
    if isinstance(ident, dict):
        return ident
    mcp_run = getattr(run, "mcp_run", None)
    if isinstance(mcp_run, dict) and mcp_run:
        return merge_run.identity_summary(mcp_run, [mcp_run], "arms/mcp/run.json")
    return None


def service_version_text(run, full=False):
    """The service's deployed revision as the MCP runner recorded it (the commit it read from the service before each
    invocation, run.json service_identity; whether it changed or a read failed), else as the service's capabilities
    document or an identity input of this document states it (make_doc --identity, the v6 delta: taken from the
    deployment log, labelled as not recorded by the runner), or the statement that this run did not record it
    (READ-021: never a nearby commit in its place). A Railway deployment id comes only from the identity input,
    labelled as taken from it (the service states its commit only); the commit is shortened unless `full`."""
    ident = getattr(run, "identity", None)
    ident = ident if isinstance(ident, dict) and ident.get("service_commit") else None
    source = str((ident or {}).get("source") or "an input file of this document").split(";")[0].strip()
    recorded = recorded_service_identity(run) or {}
    commits = recorded.get("service_commits") or {}
    if commits:
        parts = ["service commit %s" % ", ".join(_commit_text(c, full) for c in commits)]
        if len(commits) > 1:
            parts[0] += " (%s)" % ", ".join("%s on %d invocation(s)" % (_commit_text(c, full), k)
                                           for c, k in commits.items())
        if recorded.get("identity_changed") is True:
            parts.append("changed during %d invocation(s)" % max(1, recorded.get("changed") or 0))
        if recorded.get("without_commit"):
            parts.append("not read on %d of %d invocation(s)" % (recorded["without_commit"], recorded["invocations"]))
        if recorded.get("unconfirmed"):
            parts.append("read once only on %d invocation(s) (a read failed)" % recorded["unconfirmed"])
        text = "the service's deployed revision as the MCP runner recorded it (read from the service before each MCP " \
               "invocation): %s" % ", ".join(parts)
        dep = (ident or {}).get("deployed_service_commit")
        if isinstance(dep, dict) and dep.get("commit") and len(commits) == 1:
            # held to the deployed commit of the deployment record (RV3, V3)
            same = gate_mod.commit_matches(list(commits)[0], dep["commit"])
            text += "; %s the deployed commit%s, taken from %s%s" % (
                "equal to" if same else "not", "" if same else " %s" % _commit_text(dep["commit"], full),
                str(dep.get("source") or "the identity input"),
                (" (Railway deployment %s)" % dep["railway_deployment_id"]) if dep.get("railway_deployment_id") else "")
        if ident and ident.get("railway_deployment_id"):
            text += "; Railway deployment %s, taken from %s (not recorded by the run)" % (
                ident["railway_deployment_id"], source)
        return text
    mcp_run = run.mcp_run or {}
    for key in SERVICE_VERSION_KEYS:
        if mcp_run.get(key):
            return "the service's deployed revision as the MCP runner recorded it: %s" % mcp_run[key]
    caps = mcp_run.get("service_capabilities") if isinstance(mcp_run.get("service_capabilities"), dict) else {}
    doc = caps.get("document") if isinstance(caps.get("document"), dict) else {}
    for key in ("service_version", "version", "release", "commit"):
        if isinstance(doc.get(key), (str, int)) and doc.get(key):
            return "the service's deployed revision as its capabilities document states it: %s" % doc[key]
    if ident:
        parts = ["service commit %s" % _commit_text(ident["service_commit"], full)]
        if ident.get("railway_deployment_id"):
            parts.append("Railway deployment %s" % ident["railway_deployment_id"])
        return "the service's deployed revision, taken from %s and not recorded by the runner: %s" % (
            source, ", ".join(parts))
    return ("the service's deployed revision is not recorded by this run: the MCP requests cannot be attributed to a "
            "service revision, and a later run cannot be shown to have met the same service code")


def relation_example_text(run, rtype):
    """One relation of the run of the given type, in words (READ-019): 'REL-subset-0001 (M-twin-0001 ⊆ M-twin-0002:
    the title)'; 'none in this run' when the manifest has none."""
    signs = {"subset": "⊆", "same_answer": "=", "prefix_equals": "⊑", "zero_equals_absent": "≡"}
    for rel in run.relations:
        if rel.get("type") == rtype:
            return "%s (%s %s %s: %s)" % (rel.get("id"), rel.get("a"), signs.get(rtype, rtype), rel.get("b"),
                                         str(rel.get("title") or "").rstrip("."))
    return "none in this run"


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
            ("access paths (arms)", "; ".join("%s %s to %s" % (esc(arm), _date((info or {}).get("started_utc")),
                                               _date(_arm_end(run, arm, info)))
                               for arm, info in sorted(arms.items())) or "–"),
            ("set", "the re-evaluated set (checks commit %s)" % tt(_short(((run.reeval or {}).get("checks") or {}).get(
                "commit"), 12)) if run.reevaluated else "the recorded set")]
    if ctx.sample or set(ctx.chapters) != set(ctx.units.chapters):
        rows.append(("chapters built", "%s of %s (Part II holds %s; the whole-run summary is shown for context)" % (
            esc(", ".join(ctx.chapters)), esc(", ".join(ctx.units.chapters)), esc(", ".join(ctx.chapters)))))
    return _kv_table(rows)


FINDING_KIND_WORDS = {
    "late_503": ("late 503 deadline (no results carried; beyond the suite's tolerance on the client's clock)",
                 "late 503 deadlines (no results carried; beyond the suite's tolerance on the client's clock)"),
    "over_200": ("200 past its budget (the host's own elapsed time)", "200s past their budget (the host's own elapsed time)"),
    "walk_200": ("walk stopped outside its stated rule (the seed's clock)",
                 "walks stopped outside their stated rule (the seeds' clocks)"),
    "other": ("other stated bound broken (the detail block names it)", "other stated bounds broken (the detail blocks name them)"),
}


def finding_kinds(run):
    """The run's finding executions by what the check read: {late_503: a 503 deadline delivered late (the client's
    clock), over_200: a 200 of /resolve or /pattern whose budget_respected item is a finding (the host's own elapsed
    time), walk_200: a 200 of /traverse whose budget_respected item is a finding (a time stop before the budget, or
    a walk past its pieces after the seed phase used the budget: checks.check_budget_respected, the seed's clock;
    the re-verification L-b), other: any other stated bound}."""
    counts = {"late_503": 0, "over_200": 0, "walk_200": 0, "other": 0}
    for sid in run.order:
        for arm in ARMS:
            for l in run.arm_lines(sid, arm):
                if l.get("outcome") != "finding":
                    continue
                fc = [c for c in (l.get("checks") or []) if isinstance(c, dict) and c.get("outcome") == "finding"]
                if l.get("http") == 503:
                    counts["late_503"] += 1
                elif l.get("http") == 200 and any(c.get("check") == "budget_respected" for c in fc):
                    route = str((run.scenario(sid) or {}).get("route") or "")
                    counts["walk_200" if "/traverse" in route else "over_200"] += 1
                else:
                    counts["other"] += 1
    return counts


def finding_kinds_text(run):
    """The run's findings by kind in words, plain text: '4 late 503 deadlines (no results carried; beyond the suite's
    tolerance on the client's clock)'; 'none' without a finding."""
    counts = finding_kinds(run)
    parts = ["%d %s" % (counts[k], FINDING_KIND_WORDS[k][0 if counts[k] == 1 else 1])
             for k in ("late_503", "over_200", "walk_200", "other") if counts[k]]
    return ", ".join(parts) or "none"


_FINDINGS_VERDICT = re.compile(r"^FINDINGS \((\d+) execution\(s\) where the host broke a stated bound without failing "
                               r"the request; (no failure, no disagreement)\)(.*)$", re.S)


def _count_phrase(n, noun, zero, tail=""):
    """'none in this run' for 0, '1 execution<tail>', 'N executions<tail>'."""
    return zero if not n else "%d %s%s%s" % (n, noun, "" if n == 1 else "s", tail)


_PLURAL_RE = re.compile(r"\b(\d+) ((?:[A-Za-z-]+ )*?)([a-z]+)\(s\)")


def plural_counts(text):
    """'11 inconclusive execution(s)' -> '11 inconclusive executions', '1 relation(s) failed' -> '1 relation failed',
    '0 MCP disagreement(s)' -> '0 MCP disagreements': the counts of a generated line in their proper number (RV12)."""
    return _PLURAL_RE.sub(lambda m: "%s %s%s%s" % (m.group(1), m.group(2), m.group(3),
                                                   "" if m.group(1) == "1" else "s"), text)


def doc_verdict_text(run):
    """The run verdict as this document prints it: report.verdict_text's word and figures, the reason of a FINDINGS
    verdict in the document's vocabulary (READ-003; the readability review's verification of 2026-10-16, shared
    M1): the findings by what the check read (finding_kinds_text), never 'without failing the request', since a
    late 503 deadline carried no results. Every other verdict is report.py's line as it stands."""
    # the document's convention (I.2): "access paths" in the prose, "arms" for row labels only (RV11); counts in
    # proper plurals, never "execution(s)" (RV12)
    v = (run.verdict or "–").replace("between the arms", "between the access paths")
    m = _FINDINGS_VERDICT.match(v)
    if not m:
        return plural_counts(v)
    return plural_counts("FINDINGS (%s execution(s) with a finding, a stated bound broken on the clock the check reads: "
                         "%s; %s)%s" % (m.group(1), finding_kinds_text(run), m.group(2), m.group(3)))


def verdict_line(ctx):
    """The run verdict on the cover and in I.1: the verdict word large, its reason in the document's vocabulary
    (doc_verdict_text), and for a run with a deploy gate's verdict (recorded, or judged offline for this document)
    the label report.py gives the ordinary verdict beside a gate (report.verdict_headline: the gate is the
    decision, this verdict is exploratory)."""
    v = doc_verdict_text(ctx.run)
    word, _, rest = v.partition(" ")
    out = "\\textbf{\\large %s}\\ %s" % (esc(word), esc(rest, 400))
    if ctx.run.gated:
        out += "\\\\{\\small %s}" % esc(EXPLORATORY_LABEL)
    return out


# the run verdict's label beside a deploy gate on the cover and in I.1, where the gate's box follows it (the reader
# review of the v7 document, RV10: report.py's "which is the gate's beside" is not English to a first reader)
EXPLORATORY_LABEL = "Exploratory verdict; the deploy decision is the gate's, below"


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
            how = ("judged offline by gate.py for this document (the run recorded %s); policy %s%s" % (
                "a stale gate.json: %s" % stale_gate_text(run) if run.gate_stale else "no gate.json",
                os.path.basename(pol.get("path") or "–"), deploy_label_text(run)))
            # gate.py's sha256 is that of the policy's canonical JSON; the deploy notes cite the file's bytes, so
            # both are printed (VERIFY-COLOURS-V5 item 9)
            file_sha = run.policy_file_sha256() if hasattr(run, "policy_file_sha256") else None
            how += ", %s version %s (sha256 %s… of the canonical JSON%s), access paths %s" % (
                pol.get("name") or "–", pol.get("policy_version") or "–", str(pol.get("sha256") or "–")[:12],
                (", %s… of the file" % file_sha[:12]) if file_sha else "",
                ", ".join((run.policy or {}).get("arms") or []) or "–")
        else:
            # the recorded verdict, in the words of report.py's gate section (gate_headline_parts)
            word, how, _ = report.gate_headline_parts(gate)
            how = "recorded in the run (%s; the sha256 is that of the policy's canonical JSON), access paths %s" % (
                how, ", ".join((run.policy or {}).get("arms") or gate.get("arms") or []) or "–")
        counts_text = "%s executions judged%s, %s planned, %s required" % (
            counts.get("judged", "–"), judged_outcomes_text(run, counts), counts.get("planned", "–"),
            counts.get("required", "–"))
        by_kind = "; ".join("%d %s" % (sum(1 for p in problems if p.get("kind") == k), k)
                            for k in gate_mod.KINDS if any(p.get("kind") == k for p in problems))
        problems_text = "%d problem%s%s" % (len(problems), "" if len(problems) == 1 else "s",
                                             (" (%s)" % by_kind) if by_kind else "")
        return word, how, counts_text, problems_text
    how = ("a stale gate.json in the run (%s)" % stale_gate_text(run)) if run.gate_stale else "no gate.json in the run"
    if run.gate_error:
        how += ", and gate.py could not judge it offline: %s" % latex.cut_text(latex.public_paths(str(run.gate_error)), 200)
    elif run.policy is None:
        how += ", and no policy in the suite's gate directory pins this manifest, so none was judged offline"
    return "not evaluated", how, None, None


def deploy_label_text(run):
    """' (rtp: POST /resolve, ... and POST /pattern)' after the deploy gate's policy file when another policy is
    reported beside it (RV5: the short names explained where they first appear); '' otherwise."""
    if not getattr(run, "also_gates", None) or not isinstance(getattr(run, "policy", None), dict):
        return ""
    return " (%s: %s)" % (policy_label(run.policy), facts._join(policy_routes(run.policy)) or "–")


def gate_policy_names(run):
    """The glossary's short names of the policies this document reports (RV5): 'rtp, the deploy gate (...) and rt
    (...)'; '' with one policy or none."""
    if not getattr(run, "also_gates", None) or not isinstance(getattr(run, "policy", None), dict):
        return ""
    parts = ["%s, the deploy gate (%s)" % (policy_label(run.policy), facts._join(policy_routes(run.policy)))]
    for e in run.also_gates:
        if isinstance(e.get("policy"), dict):
            parts.append("%s (%s)" % (policy_label(e["policy"]), facts._join(policy_routes(e["policy"]))))
    return facts._join(parts)


def judged_outcomes_text(run, counts, policy=None):
    """' (1950 passed, 118 service-side by design)': the gate's judged executions by outcome (V6-001), the passed
    count the gate's own; ', N passed' when the policy's population cannot be read. |policy|: another policy than
    the deploy gate's (an --also-policy one, RV5)."""
    out = gate_judged_outcomes(run, policy)
    if not out or out.get("pass", 0) != counts.get("passed"):
        return ", %s passed" % counts.get("passed", "–")
    words = [("pass", "passed"), ("service_by_design", "service-side by design")]
    parts = ["%d %s" % (out[o], w) for o, w in words if out.get(o)]
    parts += ["%d %s" % (n, o) for o, n in sorted(out.items()) if o not in dict(words)]
    if not out.get("pass"):
        parts.insert(0, "0 passed")
    return " (%s)" % ", ".join(parts)


def stale_gate_text(run):
    """Why a recorded gate.json is stale (gate.stale_inputs): the verdict of <utc> judged other contents of the
    named files than the run has now (the report's 'GATE verdict STALE')."""
    rec = run.gate_recorded or {}
    pol = rec.get("policy") or {}
    return "the verdict of %s (policy %s, %s) judged other contents of %s than the run has now" % (
        str(rec.get("judged_utc") or "–")[:19].replace("T", " "), pol.get("name") or "–",
        "PASS" if rec.get("passed") else "FAIL", ", ".join(run.gate_stale))


def gate_line(ctx):
    """The deploy gate on the cover, as prominent as the run verdict (the owner's review, point 2): the word, then
    the plain decision, the problems by group and kind, the executions judged, the scope read from the policy and
    how it was judged (READ-002)."""
    word, how, counts_text, problems_text = gate_summary(ctx)
    if counts_text is None:
        return "\\textbf{\\large Deploy gate: %s}\\ (%s)" % (esc(word), esc(how))
    run = ctx.run
    passed = bool((run.gate or {}).get("passed"))
    if passed:
        # "by design" never hyphenated in the box (RV12)
        return "\\textbf{\\large %s}: nothing blocks (%s); %s; the gate covers %s; %s" % (
            esc(word), esc(problems_text), esc(counts_text).replace("by design", "\\mbox{by design}"),
            esc(gate_scope(run)["text"]), esc(_without_access_paths(how, run)))
    head, kinds = _problems_head_and_kinds(problems_text)
    return ("\\textbf{\\large %s}: the run does not authorise a deploy. %s, by what would clear them: %s (by the "
            "gate's kinds: %s); the count is not a count of service defects; %s; the gate covers %s; %s" % (
                esc(word), esc(head[:1].upper() + head[1:]), esc(gate_groups_text(run)), esc(kinds),
                esc(counts_text).replace("by design", "\\mbox{by design}"), esc(gate_scope(run)["text"]),
                esc(_without_access_paths(how, run))))


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
                                                                        latex.join_units(cells.esc_times(_finding_evidence(l)))))
                elif l.get("outcome") in ("fail", "error"):
                    failed = [c for c in (l.get("checks") or []) if isinstance(c, dict) and c.get("outcome") == "fail"]
                    what = ("%s: %s" % (esc(failed[0].get("check")),
                                        latex.join_units(cells.esc_times(failed[0].get("detail"), 120)))) if failed \
                        else esc(cells.WORD.get(l.get("outcome"), l.get("outcome")))
                    failures.setdefault(sid, []).append("%s r%s %s" % (esc(cells.ARM_SHORT.get(arm, arm)),
                                                                       esc(l.get("repeat")), what))
        e = run.comp_entry(sid)
        for key, who in (("direct_agreement", "direct"), ("mcp_agreement", "mcp")):
            a = e.get(key) or {}
            if a.get("class") == "disagree":
                diff = [str(c) for c in (a.get("differing_checks") or [])]
                if diff:
                    what = "a check differs: " + ", ".join(diff)
                else:
                    # no check differs: the outcomes do (compare's outcome_differs: a time-limited pass against a 503
                    # deadline, say), with the repeat and compare's reason (VERIFY-COLOURS-V5 item 1)
                    od = [o for o in (a.get("outcome_differs") or []) if isinstance(o, dict)]
                    what = "the outcomes differ" + ((", repeat %s: %s" % (od[0].get("repeat"), od[0].get("detail")))
                                                    if od and od[0].get("detail") else "")
                disagreements.append("%s (%s: %s)" % (_link_sid(ctx, sid), who, esc(what, 160)))
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
        out.append(("disagreements between access paths (%d)" % len(disagreements), disagreements))
    if failed_rel:
        out.append(("failed relations (%d)" % len(failed_rel), failed_rel))
    return out


def findings_lead_text(ctx, groups=None):
    """I.1's principal findings, the result first (V6-003): what failed, broke a bound, disagreed or failed as a
    relation, else that nothing did, and how many executions were inconclusive and required by the gate."""
    run = ctx.run
    groups = principal_findings(ctx) if groups is None else groups
    n_inc = sum(1 for l in run.rows if l.get("outcome") == "inconclusive")
    gp = gate_outcome_counts(run)
    if not gp["judged"]:
        gate = ", whether any is required by a gate was not judged"
    elif n_inc == 1:
        gate = ", %s by the deploy gate" % ("required" if gp["inconclusive"] else "not required")
    elif n_inc:
        gate = ", %s of them required by the deploy gate" % (gp["inconclusive"] or "none")
    else:
        gate = ""
    tail = "%d execution%s %s inconclusive%s (Properties not shown on any access path, below)" % (
        n_inc, "" if n_inc == 1 else "s", "was" if n_inc == 1 else "were", gate)
    if not groups:
        return ("*The result.* No execution failed or errored, no resource-budget finding was recorded, no "
                "comparison between access paths disagreed and no relation failed in this run; %s." % tail)
    return ("*The result.* Listed below: %s; %s." % (facts._join(t for t, _ in groups), tail))


def _arm_reason(lines, st, name=None, run=None, sid=None):
    """Why an arm did not pass a property, in a few words, from its lines. For the MCP arm a check |name| the
    manifest marks at property level (suite_common.mcp_check_marks) that read not_exposed reads as the mark says
    (not implemented via MCP, or service-side by design), never as untested; an execution inconclusive by a
    scenario-level mark (suite_common.SCENARIO_LEVEL_MARK) says so."""
    if not lines:
        return "no execution"
    if st.glyph_key == "not_implemented_mcp":
        return cells.attn_marked("not exercised via MCP")
    words = []
    is_mcp = getattr(st, "arm", None) == "mcp"
    marks = run.check_marks(sid) if (run is not None and sid and is_mcp) else {}
    for l in lines:
        if is_mcp and cells.by_design_line(st, l):
            # the case's own repeat; the other repeat outside the case reads its recorded outcome below
            # (VERIFY-QCND H2); whether the host was reached, as the chapters' grouped cells say it (READ-208)
            w = ("by design: the host reached, the pair decided by the service" if cells.host_reached(l)
                 else "by design: decided by the service before the host")
            if w not in words:
                words.append(w)
            continue
        http = l.get("http")
        code = (l.get("error") or {}).get("code") if isinstance(l.get("error"), dict) else None
        aside = run.marked_aside(sid, l) if marks else {}
        marked = next((m for (n, _), m in aside.items() if n == name), None) if name else None
        if marked:
            w = ("%s (the manifest's property-level mark: %s)" % (cells.attn_marked("not implemented via MCP"),
                                                                  marked.get("reason"))
                 if marked.get("mark") == "not_implemented" else
                 "service-side by design (the manifest's property-level mark, case %s)" % marked.get("case"))
        elif is_mcp and blocks.scenario_level_mark(l):
            w = "every property check marked as not shown via MCP (a scenario-level mark)"
        elif l.get("outcome") == "pass":
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
                            w = _arm_reason(run.arm_lines(sid, arm), ctx.state(sid)[arm], name, run, sid)
                            if w not in words:
                                words.append(w)
                        reasons.append((arm, "; ".join(words)))
                    # the expectation in words when every member expecting the check shares it (the passed ones
                    # included: the count is over all of them), else the check's name alone
                    texts = []
                    for sid in expecting:
                        e = next((e for e in run.scenario(sid).get("expect") or []
                                  if isinstance(e, dict) and e.get("check") == name), {})
                        t = blocks.expect_text(name, e.get("args"))
                        if t not in texts:
                            texts.append(t)
                    per_unit.append({"chapter": ch, "unit": uid, "check": name,
                                     "expect": texts[0] if len(texts) == 1 else latex.tt(name),
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
        out.append("\\textbf{%s} · \\hyperref[unit:%s]{\\texttt{%s}} · %s: no access path passed it on %s (%d test%s) — %s" % (
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
    return ("Over the whole run %d of %d property claims (one claim is one property check of one scenario) were not "
            "passed on any access path, in %d test%s%s." % (total_untested, total_pairs, total_tests,
                                                              "" if total_tests == 1 else "s",
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
                     cells.attn(str(o["not_implemented_mcp"])) if o.get("not_implemented_mcp") else "0",
                     str(o.get("service_by_design", 0)),
                     "%s / %s / %s / %s" % (hc.get("cold", 0), hc.get("warm", 0), hc.get("aged", 0),
                                            hc.get("not_reached", 0))])
    spec = "@{}l r r r r r r r r r@{}"
    return latex.plain_table(spec, ["access path", "executions", "\\gpass", "\\gfail", "\\gerr", "\\gfind", "\\ginc",
                                    "\\gni", "\\gbd", "cold / warm / aged / not reached"], rows,
                             "the recorded execution outcomes per access path (the report's counts, executions not "
                             "tests; \\gni\\ and \\gbd\\ are outcomes here, not the tables' standing)")


def counting_units_table(ctx):
    """The counting units of the document with this run's figure for each (READ-006): scenario, execution, property
    claim, scenario-arm pair, gate problem; the overlaps are stated in the caption."""
    run = ctx.run
    values = ctx.placeholders()
    gp = gate_outcome_counts(run)
    rows = [["scenario (test)", "one frozen request with its ground truth and expectations; the tables' test id",
             esc(str(values["n_scenarios"]))],
            ["execution", "one request of one scenario on one access path, one repeat (the paired loop runs a "
                          "scenario on three access paths, the MCP path twice where the tools can express it)",
             esc(str(values["n_executions"]))],
            ["property claim", "one property check of one scenario (a scenario names one or several); passed when "
                               "any execution on any access path passed it", esc(str(values["n_property_claims"]))],
            ["scenario-arm pair", "one scenario on one access path, whatever its repeats: the unit of the gate's "
                                  "outcome problems and of the owner's exclusions",
             ("%s with an outcome problem" % gp["pairs"]) if gp["judged"] else "no gate judged"],
            ["gate problem", "one problem the deploy gate's policy raises against the run's records: an evidence "
                             "problem, a required behaviour not shown or an observed violation; never a count of "
                             "tests or of defects (and not a finding, which is an execution outcome)",
             esc(str(values["n_gate_problems"]))]]
    return latex.plain_table(latex.spec_p([24, 112, 38]), ["unit", "what one is", "this run"], rows,
                             "the counting units: the totals of this document count different things and are not "
                             "added to one another; the limitation counts below overlap too (a scenario may carry a "
                             "manifest mark and name a table case)", needspace=10)


def chapter_map_table(ctx):
    """The question-to-chapter map of I.1 (READ-001, READ-009): one row per chapter of the run, the reader's
    question from the overview prose (`decision.chapter_map`: [{chapter, question}]), the chapter's title and
    page; a chapter the prose has no question for prints its title alone."""
    questions = {}
    for e in ctx.overview_section("decision").get("chapter_map") or []:
        if isinstance(e, dict) and e.get("chapter"):
            questions[str(e["chapter"])] = str(e.get("question") or "")
    rows = []
    for ch in ctx.units.chapters:
        title = (ctx.prose.get(ch) or {}).get("title") or ctx.units.titles.get(ch, ch)
        where = ("\\hyperref[ch:%s]{chapter %s, p.\\,\\pageref{ch:%s}}" % (esc(ch), esc(ch), esc(ch))
                 if ch in ctx.chapters else "chapter %s (not built)" % esc(ch))
        rows.append([inline(questions.get(ch) or "–"), "%s (%s)" % (esc(title), where)])
    if not rows:
        return ""
    return latex.plain_table(latex.spec_p([96, 70]), ["the question", "where the evidence is"], rows,
                             "where to look: the chapters by the question they answer (the closing section, "
                             "Timing 1 to 7, reads the timing across chapters)", needspace=min(16, len(rows) + 3))


def untested_table(ctx, chapter=None):
    """The untested and partly tested properties as a table (READ-005): the question the unit asks (the unit
    prose's `question`, else its property line), what this run established, what remains unproven (the check in
    words), why per access path, and the tests (linked); of one chapter or of every built one. Empty when nothing is
    untested."""
    rows = []
    for u in ctx.untested()["per_unit"]:
        if chapter and u["chapter"] != chapter:
            continue
        entry = ((ctx.prose.get(u["chapter"]) or {}).get("units") or {}).get(u["unit"])
        question = property_text = None
        if isinstance(entry, dict):
            question, property_text = entry.get("question"), entry.get("property")
        elif isinstance(entry, str):
            property_text = entry
        head = inline(ctx.fill(question)) if question else (
            inline(latex.cut_text(ctx.fill(property_text), 160)) if property_text else "–")
        established = ("%d of %d tests passed it on some access path" % (u["passed"], u["n"]) if u["passed"]
                       else "no test passed it on any access path (%d test%s)" % (u["n"], "" if u["n"] == 1 else "s"))
        ids = ", ".join(latex.sid(x, u["unit"], detail=ctx.footnote(x) is not None) for x in u["missing"][:8])
        if len(u["missing"]) > 8:
            ids += " … (%d)" % len(u["missing"])
        reasons = "; ".join("%s: %s" % (esc(cells.ARM_SHORT.get(a, a)), esc(w)) for a, w in u["reasons"] if w)
        rows.append(["%s\\newline{\\scriptsize\\color{black!60}\\hyperref[unit:%s]{\\texttt{%s}}}" % (
            head, latex.label_key(u["unit"]), esc(u["unit"])), esc(established), u["expect"], reasons, ids])
    if not rows:
        return ""
    return latex.plain_table(latex.spec_p([42, 24, 34, 40, 26]),
                             ["the question (unit)", "what this run established", "unproven: the check",
                              "why, per access path (dl · dr · mcp)", "tests"], rows,
                             "properties no access path showed: a test passes a property when any execution on any "
                             "access path passed it; the reasons are the recorded outcomes, not causes")


def decision_summary(ctx):
    """I.1: the orientation and the scope (prose), the question-to-chapter map, the two verdict lines and how they
    differ (prose `verdict`), the principal findings (prose `findings`, then the computed lists), the properties
    still untested (prose `untested`, then the computed table and the run-wide count) and where to look (prose
    `closing`); the per-arm outcome table is I.4's."""
    out = ctx.prose_paragraphs("decision")
    out += chapter_map_table(ctx)
    out += "\\begingroup\\sffamily\n\\bline{run verdict}{%s}\n\\bline{deploy gate}{%s}\n\\par\\endgroup\\medskip\n\n" % (
        verdict_line(ctx), gate_line(ctx))
    out += ctx.prose_paragraphs("decision", "verdict", optional=True)
    out += gate_population_table(ctx)
    out += attention_table(ctx)
    out += "\\parhead{Principal findings}\n"
    groups = principal_findings(ctx)
    # the actual result before the definitions (the external review of the v6 document, V6-003)
    out += inline(findings_lead_text(ctx, groups)) + "\n\n"
    # the host's memory episodes, said where a reader of the summary meets the result (RV1, V9)
    lead = memory_lead_text(ctx)
    if lead:
        out += inline(lead) + "\n\n"
    out += ctx.prose_paragraphs("decision", "findings", optional=True)
    for title, items in groups:
        out += "\\noindent\\textbf{%s}\\par\\nopagebreak\n\\begin{itemize}\\setlength\\itemsep{0pt}\n" % esc(title)
        for it in items[:20]:
            out += "\\item %s\n" % it
        if len(items) > 20:
            out += "\\item … %d more (Part II)\n" % (len(items) - 20)
        out += "\\end{itemize}\n"
    out += "\\parhead{Properties not shown on any access path}\n"
    if ctx.overview_section("decision").get("untested"):
        out += ctx.prose_paragraphs("decision", "untested")
    else:
        # without the prose slot the by-reason line is the generator's own
        ibr = inconclusive_by_reason_text(ctx.run)
        if ibr != "none":
            out += ("Inconclusive executions by reason: %s (the checks: a check could not test its property, "
                    "typically an allowed 503 deadline; a host state: the service's own limit ended the call).\n\n"
                    % esc(ibr))
    table = untested_table(ctx)
    if table:
        out += table
    else:
        out += "Every property check of the built chapters was passed by at least one execution of each test.\n\n"
    out += untested_context_line(ctx) + "\n\n"
    out += ctx.prose_paragraphs("decision", "closing", optional=True)
    return out


# -- Part I: methods ----------------------------------------------------------------------------------------------

def budget_scope_table(ctx):
    """The route-by-route scope of a time budget (READ-015): the rows of the overview prose `methods.budget_scope`
    ([{route, named_by, clock, scope, not_covered, answers, verified}], constants of the SPECs), rendered after the
    methods' prose; empty when the prose has none."""
    rows_in = [r for r in (ctx.overview_section("methods").get("budget_scope") or []) if isinstance(r, dict)]
    if not rows_in:
        return ""
    # one full-width block per route, its six statements labelled (the external review of the v6 document, V6-003:
    # a seven-column table squeezed the /traverse exceptions into a word-wide column); nothing is left out
    labels = (("named_by", "The budget"), ("clock", "Its clock starts"), ("scope", "It covers"),
              ("not_covered", "It does not cover"), ("answers", "Possible answers"),
              ("verified", "What the suite verified"))
    out = ("\\Needspace{6\\baselineskip}\n{\\small\\itshape the scope of a time budget per route (the SPECs' contract; "
           "the suite's rulings in appendix A, budget\\_respected)}\\par\\nopagebreak\\smallskip\n")
    for r in rows_in:
        parts = ["\\emph{%s:} %s" % (esc(label), inline(ctx.fill(str(r.get(k) or "–")).rstrip(".")))
                 for k, label in labels]
        out += ("{\\small\\noindent\\textbf{%s.} %s\\par}\\smallskip\n" % (inline(ctx.fill(str(r.get("route") or "–"))),
                                                                        " ".join(p + "." for p in parts)))
    return out + "\n"


def methods_section(ctx):
    """I.2: the prose (the access paths' bases are in it, their executions in I.4's table, their states and runners
    in appendix F) and the budget-scope blocks (READ-015; one per route since V6-003)."""
    return ctx.prose_paragraphs("methods") + budget_scope_table(ctx)


# -- Part I: reading guide ----------------------------------------------------------------------------------------

LEGEND_GLYPH = {"pass": "\\gpass", "fail": "\\gfail", "error": "\\gerr", "finding": "\\gfind",
                "inconclusive": "\\ginc", "not_implemented": "\\gni", "by_design": "\\gbd", "none": "\\gnone"}
LEGEND_WORD = {"pass": "pass", "fail": "fail", "error": "error", "finding": "finding", "inconclusive": "inconclusive",
               "not_implemented": "not exercised via MCP", "by_design": "service-side by design",
               "none": "no execution"}
LEGEND_PREFIXES = {"not_implemented": ("not exercised via MCP:", "not implemented via MCP:")}
MARK_EXAMPLE = {"superscript digit": "\\ginc\\fn{1}", "not-equal sign": "\\gpass\\gne", "t": "\\gpass\\tmark{t}",
                "property-level mark": "\\gpass\\gnip", "by-design property mark": "\\gpass\\gbdp",
                "s": "\\gpass\\tmark{s}", "n": "\\gpass\\tmark{n}", "u": "\\gbd\\tmark{u}", "c": "9.1 s\\tmark{c}",
                "a": "9.1 s\\tmark{a}", "i": "9.1 s\\tmark{i}", "m": "\\ginc\\tmark{m}"}
ANNOTATION_DEFAULT = {
    "column": "a column is one varied value of the unit's dimension (here {annotated_dimension})"
              "[[annotated_record_key ?; the record key in brackets names the source record (the key in appendix "
              "D)]]",
    "glyph": "the glyph is the access path's standing on the test: pass, fail, error, finding, inconclusive, not "
             "implemented via MCP, service-side by design (the legend below)",
    "repeats": "two numbers separated by a slash are the repeats in order (on mcp the cold call, then the warm "
               "call); the unit is written once",
    "server": "the bracketed line is the host's own server time per execution; a dash when the answer carries "
              "none (a refusal, a 503)",
    "backend": "be: the host call as the service timed it (mcp only); the gap to the wall is the service's own "
               "cost",
    "footnote": "a dark-grey superscript digit numbers the detail block after the table (numbering restarts per "
                "chapter); it sits on the offending cells only",
    "mark": "a small letter after a time is the host-cache state when it is not the usual one (c cold, a aged, i "
            "warm after an interrupted execution); after a glyph: m mixed repeats, t s n the direct classes, "
            "u an unverified by-design pair",
    "ids": "the ids line names the members in column order (the common prefix once); a unit id line and the index "
           "find every test",
    "summary": "the family's summary line counts its tests per access path by standing and gives the wall p50 / "
               "p90 (median / 90th percentile) per access path, host-cold and host-warm apart",
}


def legend_table(ctx=None):
    """The glyph legend and the marks of the reading guide. The glyphs' wording is the shared vocabulary's
    (vocabulary.legend_rows: the same definitions appendix F's tables and the glossary print, so that the three
    cannot drift; READ-022), the symbols the macros the tables use; the marks' wording comes from the overview
    prose (`marks` of the reading section) when it has them."""
    reading = ctx.overview_section("reading") if ctx else {}
    marks = [e for e in (reading.get("marks") or []) if isinstance(e, dict) and e.get("mark") in MARK_EXAMPLE]
    # the property-level mark exists only under a manifest with property-level marks (v5): no row otherwise
    if ctx is not None and not _has_property_marks(ctx.run):
        marks = [e for e in marks if e.get("mark") != "property-level mark"]
    # the by-design kind only under a manifest with such marks (v7; RV4)
    if ctx is not None and not _has_property_marks(ctx.run, "by_design"):
        marks = [e for e in marks if e.get("mark") != "by-design property mark"]
    glyph_rows = _glyph_rows(ctx)
    if marks:
        mark_rows = [[MARK_EXAMPLE[e["mark"]], esc(e.get("on") or ""), inline(e.get("meaning") or "")]
                     for e in marks]
        return _legend_tables(glyph_rows, mark_rows)
    return _legend_builtin(glyph_rows)


def _has_property_marks(run, kind="not_implemented"):
    """Whether any scenario of the run's manifest carries a property-level mark of |kind| (not_implemented,
    by_design)."""
    if not hasattr(run, "check_marks"):
        return False
    return any(m.get("mark") == kind for sid in run.order for m in run.check_marks(sid).values())


def _legend_tables(glyph_rows, mark_rows):
    spec = latex.spec_p([8, 36, 128], aligns=["centering", "raggedright", "raggedright"])
    glyphs = latex.plain_table(spec, ["", "glyph", "meaning"], glyph_rows)
    spec2 = latex.spec_p([30, 34, 108])
    # a longtable: it may continue on the next page under its repeated header (no heading is orphaned)
    marks_t = latex.plain_table(spec2, ["mark", "on", "meaning"], mark_rows)
    return glyphs + marks_t


def answer_side_mark_ids(run):
    """The run's scenarios with an answer-side scenario-level mark (blocks.scenario_mark_parts): the request is
    expressible, the MCP answer lacks every field the property checks read (v6's P-scope-0004 / -0006)."""
    return [sid for sid in run.order if blocks.scenario_mark_parts(run.scenario(sid))]


def answer_side_mark_text(run):
    """'P-scope-0004, P-scope-0006: absence_scope': the answer-side scenario-level marks with the fields they name."""
    ids = answer_side_mark_ids(run)
    fields = []
    for sid in ids:
        for _, check, _, _ in blocks.scenario_mark_parts(run.scenario(sid)):
            for f in blocks.mark_fields(run.scenario(sid), check) or [check]:
                if f not in fields:
                    fields.append(f)
    return ("%s: %s" % (", ".join(ids), ", ".join(fields))) if ids else ""


def _glyph_rows(ctx=None):
    """The glyph legend's rows from the shared vocabulary; the not-implemented word in the attention red, as the
    glyph and its counts are (the owner's decision of 2026-10-16); the not-implemented definition widened where the
    run has an answer-side scenario-level mark."""
    side = bool(ctx is not None and answer_side_mark_ids(ctx.run))
    return [[g, cells.attn(esc(w)) if g == "\\gni" else esc(w), inline(m)] for g, w, m in vocabulary.legend_rows(side)]


def _legend_builtin(glyph_rows=None):
    glyph_rows = glyph_rows or _glyph_rows()
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


def times_legend_sentence():
    """The one sentence of the cell grammar on the tints (the owner's decision of 2026-10-15): the thresholds
    printed from cells' constants, the two words on their own tints."""
    return ("a printed execution time (a wall, a server time, a backend time or an MCP phase, alone or as a median, "
            "p90 or maximum) above %s is set on a %s tint and above %s on a %s one, decided by the value as printed "
            "(a time printed %s carries none, one printed %s is red), in the tables, the detail blocks and the prose "
            "alike, a time quoted from a check's text and the time a walk's stop used too; a difference (a delta, an "
            "overhead, how late a 503 was), a budget, a cap and a distance in time carry no tint" % (
                fmt_t(cells.TIME_TINT_YELLOW_MS), cells.TINT_MACRO["yellow"] + "{yellow}",
                fmt_t(cells.TIME_TINT_RED_MS), cells.TINT_MACRO["red"] + "{red}",
                fmt_t(cells.TIME_TINT_YELLOW_MS), fmt_t(cells.TIME_TINT_RED_MS + 1000)))


def attention_legend_sentence(answer_side=False):
    """The legend's one definition of each attention kind (the owner's decision of 2026-10-16: red means needs
    attention), from facts.ATTENTION_DEFINITION, the kinds' words in red as the contents and I.1 print them;
    `answer_side`: the run has an answer-side scenario-level mark (the not-implemented definition widened)."""
    defs = []
    for k, one, many in facts.ATTENTION_KINDS:
        d = facts.ATTENTION_DEFINITION[k]
        if k == "not_implemented" and answer_side:
            d = facts.ATTENTION_NOT_IMPLEMENTED_ANSWER_SIDE
        if k == "time_budget":
            d = d % fmt_t(cells.TIME_TINT_RED_MS)
        defs.append("%s: %s" % (cells.attn(esc(many)), esc(d)))
    return ("red means needs attention: the contents print under each chapter's line, and I.1's attention table "
            "per chapter, the tests needing attention by kind, a test once per kind and possibly under several, the "
            "kinds with none left out; %s; the \\gni\\ glyph of a test not implemented via MCP is red wherever "
            "printed, and so are its counts (a count of 0 stays black), the class words of a cell, an execution or a "
            "run statement that stands not implemented via MCP (not exercised via MCP, not expressible through MCP, "
            "refused before host execution where the arm stands \\gni) and the property-level mark (a small red "
            "\\gni\\ raised after the glyph of an MCP cell, \\gpass\\gnip: a property check set aside by the "
            "manifest's mark, the test counted under not implemented via MCP); headings, legends, definitions, quoted "
            "texts (a check's text, the policy's reasons) and the by-design table's own statements name the class in "
            "black; a footnote digit is dark grey, not red: it numbers the test's detail block, an exception to "
            "read, and marks no attention by itself" % "; ".join(defs))


def grammar_lines(ctx):
    """The rules of the cell grammar the reading guide's prose does not state: the replicates aggregate, the edge
    tables' MCP cell and the tints of the times (the rest, the lines of a cell, the rounding and the slash, is in
    the prose and in the legend line under every table)."""
    out = ("\\begingroup\\sffamily\\footnotesize\n"
           "\\bline{aggregates}{a grouped (replicates) cell is ok n/n · wall median / maximum over each member's "
           "first host-warm execution (else its only one) and [server median]; n counts the group's tests (source "
           "windows), not repeats; a wall not all host-warm carries the state's mark, alone when shared "
           "(7.0 / 7.8 s\\tmark{c}), else with its share (\\tmark{c 2/6}); members never exercised read not "
           "exercised: n, never ok 0/0; the coverage line under the table counts per access path}\n"
           "\\bline{edge tables}{the mcp cell of a refusals or list table shows the recorded status, the service's "
           "code and the class in words (refused before host execution, not expressible through MCP, answered after "
           "service normalization, host refusal relayed through the service), derived from the evidence, never from prose}"
           "\n\\bline{colour}{red marks what needs attention (I.1's table counts the kinds); a time above 10 s is set "
           "on a yellow tint, above 55 s on a red one; a footnote digit is dark grey; the full rules, with the "
           "attention kinds, are in appendix F's vocabulary}"
           "\n\\par\\endgroup\n\n")
    return out


def reading_reference(ctx):
    """The reading guide's reference rules, moved from I.3 by the second external review (R2-003): the standing
    against the outcomes row with its two examples, and the rounding of times (overview.json
    `appendix_F.reading_reference`); nothing when the prose has none."""
    entry = ctx.overview.get("appendix_F") if isinstance(ctx.overview.get("appendix_F"), dict) else {}
    items = entry.get("reading_reference") or []
    if not items:
        return ""
    return "\\parhead{Reading the tables: the reference rules}\n" + paragraphs([ctx.fill(p) for p in items])


def colour_rules(ctx):
    """Appendix F's colour rules (moved from I.3's grammar by the external review of the v6 document, V6-003): the
    tints of the printed times and the attention kinds with their red."""
    return ("\\begingroup\\sffamily\\footnotesize\n"
            "\\bline{times}{%s}" % times_legend_sentence() +
            "\n\\bline{attention}{%s}" % attention_legend_sentence(bool(ctx is not None and answer_side_mark_ids(ctx.run))) +
            "\n\\par\\endgroup\n\n")


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
    # a note whose marker sits in no cell of this matrix says so (VERIFY-COLOURS-V5 item 4)
    unmarked = []
    if not footnote_marked:
        used.append("footnote") if "footnote" not in used else None
        unmarked.append("footnote")
    if not mark_marked and "mark" not in used:
        used.append("mark")
        unmarked.append("mark")
    ids = units_ids(ctx, members) + mark("ids")
    out = "\\noindent{\\small The matrix of unit \\hyperref[unit:%s]{\\texttt{%s}} (%s), annotated:}\\par\\nopagebreak\\smallskip\n" % (
        latex.label_key(uid), esc(uid), esc(u.get("unit_title") or ""))
    out += latex.series_table(columns, rows, ids, cells.TIMING_LEGEND["series"] + mark("summary"))
    out += "\\begingroup\\footnotesize\n"
    # the column note names this unit's dimension and the record-key clause only where the variants carry a key
    # (VERIFY-COLOURS-V5 item 3)
    extra = {"annotated_dimension": str(u.get("unit_title") or "the varied value").split(" (")[0],
             "annotated_record_key": any(_variant(meta[sid], sid).endswith(")") and " (" in _variant(meta[sid], sid)
                                         for sid in members)}
    for i, key in enumerate(used):
        text = texts.get(key, key) + (" (not in this matrix)" if key in unmarked else "")
        out += "\\annline{%d}{%s}\n" % (i + 1, inline(ctx.fill(text, extra)))
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
    """The label of the compact record key: one per chapter (a full build prints the key on every chapter page),
    rk:ch-<chapter>, never rk:<chapter>: chapter D's would be appendix D's rk:D, the full key's label, and LaTeX
    would take the last of the two for every reference (VERIFY-FIX2 L5)."""
    return "rk:ch-%s" % chapter if chapter else "rk:compact"


def record_key_paragraph(ctx, chapter=None):
    """The compact record key, one paragraph (abbreviation and organism, in the manifest's order), printed on the
    chapter page before its first family, where the matrices need it, and labelled rk:ch-<chapter> so that every
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
        rows.append(("checks", "the recorded set: the outcomes as the runners' checks judged them during the run; no "
                               "re-evaluation was run afterwards, so a change to the checks since the run is not "
                               "reflected here"))
    if table.get("path"):
        text = "version %s, %s cases, service commit %s…, sha256 %s…" % (
            esc(table.get("version")), esc(table.get("cases")), tt(_short(table.get("service_commit"), 12)),
            tt(_short(table.get("sha256"), 12)))
        if table.get("replayed"):
            text += ", the re-evaluation's"
        elif table.get("replayed_sha256"):
            text += " (\\textbf{not} the table the re-evaluation replayed: appendix F)"
        tie = by_design_tie_text(run)
        if tie:
            text += "; " + esc(tie)
        rows.append(("by-design table", text))
    elif table.get("warnings"):
        rows.append(("by-design table", "none read: %s" % esc(latex.public_paths("; ".join(table["warnings"])), 300)))
    for arm, info in sorted((r.get("arms") or {}).items()):
        info = info or {}
        rows.append(("access path %s" % esc(arm), "%s executions, %s to %s, runner %s" % (
            esc(info.get("executions_recorded", "–")), _date(info.get("started_utc")), _date(_arm_end(run, arm, info)),
            esc(info.get("runner_host") or (mcp_runner_text(run) if arm == "mcp" else "–")))))
    sessions = _sessions(run)
    if sessions:
        first = sessions[0]
        loop = first.get("loop") or {}
        applied = _applied_offset(run)
        rows.append(("sessions", "%d (%s); the loop %s of %s scenarios%s; mex clock offset %s" % (
            len(sessions), esc(first.get("mode") or "–"), esc(loop.get("completed", "–")), esc(loop.get("scenarios", "–")),
            (", %s" % esc(loop.get("state"))) if loop.get("state") else "",
            ("%+.3f s applied to direct-local" % applied) if applied is not None else "not applied")))
    rows.append(("comparison", comparison_lines_text(run, short=True)))
    for entry in getattr(run, "also_gates", None) or []:
        rows.append(("other gate policy", esc(also_gate_text(entry, run))))
    if run.reruns:
        rows.append(("warm-ups run again", reruns_text(run)))
    lf = run.load or {}
    caveat = ""
    if (lf.get("sra_mean") or 0) > LOAD_CAVEAT_SRA or (lf.get("load_mean") or 0) > LOAD_CAVEAT_LOAD1:
        caveat = ": the Logan caveat applies (timings observed under concurrent host load)"
    rows.append(("host load", "%d samples; sra\\_all CPU mean %s\\%%, load1 mean %s%s" % (
        lf.get("n", 0), fmt_n(lf.get("sra_mean")), fmt_n(lf.get("load_mean"), 1), caveat)))
    return _kv_table(rows)


ROUTE_ORDER = ("POST /resolve", "POST /traverse", "GET /traverse/capabilities", "POST /pattern")


def policy_label(pol):
    """A policy's short name as this document uses it: the middle of its name ('v7-rtp-mcp' -> 'rtp'), else the
    name."""
    name = str((pol or {}).get("name") or "")
    m = re.match(r"\Av\d+-(.+)-[a-z]+\Z", name)
    return m.group(1) if m else (name or "–")


def policy_routes(pol):
    """The routes a policy selects by (its counts.by_route), in the document's order."""
    routes = list(((pol or {}).get("counts") or {}).get("by_route") or {})
    return sorted(routes, key=lambda r: (ROUTE_ORDER.index(r) if r in ROUTE_ORDER else len(ROUTE_ORDER), r))


def also_gate_text(entry, run=None):
    """A further gate policy judged offline for this document beside the deploy gate (make_doc --also-policy), plain
    text: its short name and routes against the deploy gate's (RV5), the policy, its verdict and counts with the
    judged executions by outcome, and that it is not the deploy gate."""
    name = os.path.basename(entry.get("path") or "") or "–"
    if entry.get("error") or not isinstance(entry.get("verdict"), dict):
        return "%s: not judged (%s)" % (name, latex.cut_text(latex.public_paths(str(entry.get("error") or "–")), 200))
    gate = entry["verdict"]
    pol = gate.get("policy") or {}
    counts = gate.get("counts") or {}
    problems = [p for p in (gate.get("problems") or []) if isinstance(p, dict)]
    own = entry.get("policy") or {}
    routes = policy_routes(own)
    deploy = run.policy if run is not None and isinstance(getattr(run, "policy", None), dict) else {}
    without = [r for r in policy_routes(deploy) if r not in routes]
    head = "the %spolicy %s (%s%s)" % ("narrower " if without and set(routes) < set(policy_routes(deploy)) else "",
                                      policy_label(own), facts._join(routes) or "–",
                                      (", without %s" % facts._join(without)) if without else "")
    judged = judged_outcomes_text(run, counts, own) if run is not None else ", %s passed" % counts.get("passed", "–")
    return ("%s%s, %s version %s (sha256 %s… of the canonical JSON), access paths %s: %s, %d problem%s; %s executions "
            "judged%s, %s required; judged offline by gate.py for this document beside the deploy gate, which it does "
            "not replace" % (
                (head + ": ") if own.get("name") else "", name, pol.get("name") or "–", pol.get("policy_version") or "–", str(pol.get("sha256") or "–")[:12],
                ", ".join(own.get("arms") or []) or "–",
                "GATE PASS" if gate.get("passed") else "GATE FAIL", len(problems), "" if len(problems) == 1 else "s",
                counts.get("judged", "–"), judged, counts.get("required", "–")))


def comparison_lines_text(run, short=False):
    """Whether the comparison read the very lines this document reads: compare.json's `executions_sha256` (stated
    since 2026-10-15; the gate judges only those lines with it) against the executions file's sha256 now."""
    stated = run.compare_executions_sha256
    now = run.executions_sha256
    if not stated:
        text = "%s states no executions sha256 (a comparison before 2026-10-15); this document reads %s, sha256 %s…" % (
            tt(sc.reevaluated_name("compare.json", run.reevaluated)), tt(run.executions_file), tt(_short(now, 12)))
        return text
    if stated == now:
        return "%s read %s with sha256 %s…, the lines this document reads" % (
            tt(sc.reevaluated_name("compare.json", run.reevaluated)), tt(run.executions_file), tt(_short(now, 12)))
    return ("%s read %s with sha256 %s…; the file now has sha256 %s…: \\textbf{the comparison's classes describe "
            "other lines} (recompute it)%s" % (
                tt(sc.reevaluated_name("compare.json", run.reevaluated)), tt(run.executions_file),
                tt(_short(stated, 12)), tt(_short(now, 12)), "" if short else "; the gate names this as an identity problem"))


def reruns_text(run):
    """The warm-ups run again for a resume in one line: how many, per arm, what they are (timing evidence); the
    re-runs no execution of the target followed on that arm (target_execution null) are counted apart."""
    per_arm = run.rerun_counts()
    orphans = sum(1 for r in run.reruns if not r.get("target_execution"))
    return ("%d (%s): after an interruption between a warm-up and its target the warm-up ran again right before "
            "the target%s; timing evidence beside the target's execution, never an outcome (appendix F)" % (
                len(run.reruns), ", ".join("%s %d" % (esc(a), n) for a, n in sorted(per_arm.items())) or "per access path not stated",
                ", except %d that no execution of the target followed" % orphans if orphans else ""))


def reruns_table(ctx):
    """Appendix F: every warm-up run again for a resume (run.json reruns_for_resume), as the report's section
    lists them: the warm-up, its arm and repeat, wall, http, the target and the target execution it ran before."""
    run = ctx.run
    if not run.reruns:
        return ""
    rows = []
    for r in run.reruns:
        t = r.get("target_execution")
        rows.append([latex.breakable_id(str(r.get("scenario"))), esc(r.get("arm")), esc(r.get("repeat")),
                     cells.fmt_exec(r.get("wall_ms")), esc(r.get("http")), latex.breakable_id(str(r.get("rerun_target"))),
                     ("repeat %s" % esc(t[2])) if isinstance(t, list) and len(t) > 2 else "did not run after it",
                     fmt_n(r.get("before_s"), 3) if isinstance(t, list) else "–"])
    out = "\\parhead{Warm-ups run again for a resume}\n"
    out += ("A resume after an interruption between a warm-up and its target ran the warm-up again right before the "
            "target, on every access path that ran the target (suite\\_common.resume\\_rewarms). These re-runs are timing "
            "evidence only: the warm-up's outcome is its first execution's, and no table of Part II counts them.\n\n")
    out += latex.plain_table(latex.spec_p([34, 18, 10, 14, 10, 34, 22, 14]),
                             ["warm-up", "access path", "repeat", "wall", "http", "target", "target execution",
                              "s before the target"], rows)
    return out


def limitations_table(ctx):
    """Limitations and exemption counts with links to the full entries."""
    run = ctx.run
    limits = mcp_limit_counts(run)
    ma = run.summary().get("mcp_agreement") or {}
    rows = [[cells.attn("not implemented via MCP") + " by the manifest", (cells.attn if limits["manifest_not_implemented"]
                                                                           else str)("%d scenarios" % limits[
                                                                               "manifest_not_implemented"])
             + ", %d distinct reasons" % len(limits["reasons"]), "\\hyperref[app:G]{G}"],
            ["not-expressible table cases scenarios name", "%d cases, %s" % (
                len(limits["named"]), (cells.attn if limits["not_expressible_scenarios"] else str)(
                    "%d scenarios" % limits["not_expressible_scenarios"])), "\\hyperref[app:G]{G}"],
            ["by-design cases an MCP answer matched", "%d cases; pairs verified %s, unverified %s" % (
                len(limits["matched"]), esc(ma.get("by_design_verified", "–")), esc(ma.get("by_design_unverified", "–"))),
             "\\hyperref[app:B]{B}"]]
    dec = declared_counts(run)
    if gate_population(run) is not None:
        rows.append(["the gate policy's automatic exclusions", esc(gate_auto_text(run)), "\\hyperref[app:F]{F}"])
    if isinstance(run.declared, dict):
        rows.append(["owner-declared exemptions",
                     "%d exclusion%s, %d check%s not required, %d relation%s excluded" % (
                         dec["exclusions"], "" if dec["exclusions"] == 1 else "s", dec["not_required"],
                         "" if dec["not_required"] == 1 else "s", dec["relations_excluded"],
                         "" if dec["relations_excluded"] == 1 else "s"),
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
        blocking = "%d block the gate (required by the policy on that access path and repeat)" % gp["inconclusive"]
    else:
        blocking = "the gate was not judged"
    rows.append(["inconclusive executions", "%d; %s; each footnoted in its chapter's tables" % (n_inc, esc(blocking)),
                 "\\hyperref[sec:I1]{I.1}"])
    return latex.plain_table(latex.spec_p([62, 94, 10]), ["limitation", "count", "in full"], rows,
                             needspace=len(rows) + 3)


# a changed check argument in plain words ('Since the previous edition'), with the argument's value; an argument not
# listed here is named as it is
ARGUMENT_CHANGE_WORDS = {
    ("budget_respected", "reserve_ms", "added"): "a /pattern time stop is now held to the host's documented "
                                                "finalisation reserve of %s ms, no longer to the reserve the answer "
                                                "itself reports",
}
_PREV_OUTCOME_WORDS = {"pass": "pass", "inconclusive": "inconclusive", "service_by_design": "service-side by design",
                       "not_implemented_mcp": "not exercised via MCP", "fail": "fail", "error": "error",
                       "finding": "finding", None: "no execution"}


def _mark_sets(scenario):
    """(scenario-level marks, property-level marks) of a scenario's MCP block: ({'not_implemented', 'by_design'} the
    scenario carries, {(check, mark, reason or case)})."""
    mcp = (scenario or {}).get("mcp") if isinstance((scenario or {}).get("mcp"), dict) else {}
    level = set(k for k in ("not_implemented", "by_design") if mcp.get(k))
    prop = set((name, m.get("mark"), str(m.get("reason") or m.get("case") or ""))
               for name, m in sc.mcp_check_marks(scenario).items())
    return level, prop


def _argument_changes(prev, cur):
    """{(check, argument, 'added' | 'removed' | 'changed'): (scenarios, routes, values)} of the expectations of the
    scenarios both manifests hold, item by item (check name and occurrence)."""
    out = {}
    for sid, s in cur.items():
        p = prev.get(sid)
        if not p:
            continue
        def items(sc_):
            seen, res = {}, {}
            for e in sc_.get("expect") or []:
                if isinstance(e, dict) and e.get("check"):
                    k = (e["check"], seen.get(e["check"], 0))
                    seen[e["check"]] = k[1] + 1
                    res[k] = e.get("args") if isinstance(e.get("args"), dict) else {}
            return res
        a, b = items(p), items(s)
        for key in set(a) & set(b):
            for arg in set(a[key]) | set(b[key]):
                if arg in ("rule",):   # the rule's prose is a citation, not a judged value
                    continue
                kind = ("added" if arg not in a[key] else "removed" if arg not in b[key] else
                        "changed" if a[key][arg] != b[key][arg] else None)
                if kind:
                    v = out.setdefault((key[0], arg, kind), (set(), set(), set()))
                    v[0].add(sid)
                    v[1].add(str(s.get("route") or "–"))
                    v[2].add(json.dumps(b[key].get(arg), sort_keys=True))
    return out


def since_previous(ctx):
    """I.4's 'Since the previous edition' (the reader review of the v7 document, RV7): computed from the previous
    edition's run given as an input (make_doc --previous, reader.PreviousRun): the scenarios, the manifest's MCP marks
    added and removed with this run's MCP outcomes of the tests concerned, the check arguments changed, the service
    identity recorded or not, and the executions whose outcome changed, the inconclusive ones counted; nothing
    without that input."""
    prev = getattr(ctx.run, "previous", None)
    if prev is None:
        return ""
    run = ctx.run
    pm, cm = prev.manifest, run.manifest
    items = []
    p_ids, c_ids = set(prev.scenarios), set(run.order)
    plines, clines = prev.lines(), dict(((l.get("scenario"), l.get("arm"), l.get("repeat")), l) for l in run.rows)
    if p_ids == c_ids:
        items.append("the same %d scenarios; %d executions against %d" % (len(c_ids), len(clines), len(plines)))
    else:
        items.append("%d scenarios against %d (%d added, %d removed); %d executions against %d" % (
            len(c_ids), len(p_ids), len(c_ids - p_ids), len(p_ids - c_ids), len(clines), len(plines)))
    # the manifest's MCP marks, scenario by scenario
    removed_level, removed_prop, moved, added = [], {}, [], []
    for sid in run.order:
        if sid not in prev.scenarios:
            continue
        (pl, pp), (cl, cp) = _mark_sets(prev.scenarios[sid]), _mark_sets(run.scenario(sid))
        if "not_implemented" in pl and "not_implemented" not in cl:
            removed_level.append(sid)
        p_checks = dict(((c, m), r) for c, m, r in pp)
        c_checks = dict(((c, m), r) for c, m, r in cp)
        for (check, mark), reason in p_checks.items():
            if (check, mark) in c_checks:
                continue
            other = [(m2, r2) for (c2, m2), r2 in c_checks.items() if c2 == check and m2 != mark]
            if other:
                moved.append((sid, check, mark, other[0][0], other[0][1]))
            elif mark == "not_implemented":
                removed_prop.setdefault(check, []).append(sid)
        for (check, mark), reason in c_checks.items():
            if (check, mark) not in p_checks and not any(c2 == check for c2, _ in p_checks):
                added.append((sid, check, mark, reason))

    def mcp_outcome_text(sids):
        outs = sorted(set(l.get("outcome") for k, l in clines.items() if k[0] in sids and k[1] == "mcp"),
                      key=str)
        if outs == ["pass"]:
            return "every MCP execution passed"
        return "their MCP outcomes: %s" % ", ".join(_PREV_OUTCOME_WORDS.get(o, str(o)) for o in outs)
    marks = []
    if removed_level:
        marks.append("removes %d scenario-level mark%s not implemented via MCP (%s), so %s now run%s via MCP" % (
            len(removed_level), "" if len(removed_level) == 1 else "s", facts.ids_ranges(removed_level),
            "these tests" if len(removed_level) > 1 else "the test", "" if len(removed_level) > 1 else "s"))
    for check, sids in removed_prop.items():
        marks.append("removes %d property-level mark%s not implemented via MCP on %s (%s), so that check is now "
                     "compared via MCP" % (len(sids), "" if len(sids) == 1 else "s", check, facts.ids_ranges(sids)))
    for sid, check, old, new, why in moved:
        marks.append("marks %s's %s %s instead of %s (%s)" % (
            sid, check, "service-side by design" if new == "by_design" else "not implemented via MCP",
            "not implemented via MCP" if old == "not_implemented" else "service-side by design",
            ("case %s, appendix G" % why) if new == "by_design" and why else "appendix G"))
    for sid, check, mark, why in added:
        marks.append("adds a property-level mark on %s's %s (%s)" % (
            sid, check, "service-side by design" if mark == "by_design" else "not implemented via MCP"))
    if marks:
        touched = set(removed_level) | set(s for v in removed_prop.values() for s in v)
        text = "manifest v%s %s" % (cm.get("manifest_version"), "; ".join(marks)) if pm.get("manifest_version") != \
            cm.get("manifest_version") else "the manifest %s" % "; ".join(marks)
        if touched:
            text += "; of the %d tests whose mark was removed, %s" % (len(touched), mcp_outcome_text(touched))
        items.append(text)
    for (check, arg, kind), (sids, routes, values) in sorted(_argument_changes(prev.scenarios, run.scenarios).items()):
        words = ARGUMENT_CHANGE_WORDS.get((check, arg, kind))
        value = ", ".join(sorted(values)) if len(values) <= 3 else "%d values" % len(values)
        text = "the check %s %s the argument %s on %d scenarios (%s)" % (
            check, {"added": "gains", "removed": "loses", "changed": "changes"}[kind], arg, len(sids),
            ", ".join(sorted(routes)))
        if words and len(values) == 1:
            text += ": " + (words % value if "%s" in words else words)
        items.append(text)
    p_ident = (prev.run.get("service_identity") or {}) if isinstance(prev.run.get("service_identity"), dict) else {}
    c_ident = recorded_service_identity(run) or {}
    if c_ident.get("service_commits") and not p_ident.get("service_commits"):
        items.append("the MCP runner records the service's deployed commit before each invocation, which the "
                     "previous run's runner did not (the table above and appendix F)")
    # the executions whose outcome changed, by change
    changes = {}
    for key in sorted(set(plines) | set(clines), key=lambda k: (str(k[0]), str(k[1]), str(k[2]))):
        a = (plines.get(key) or {}).get("outcome")
        b = (clines.get(key) or {}).get("outcome")
        if a != b and key[0] in c_ids:
            changes.setdefault((a, b), []).append(key)
    if changes:
        parts = []
        for (a, b), keys in sorted(changes.items(), key=lambda kv: (-len(kv[1]), str(kv[0]))):
            by_sid = {}
            for sid, arm, rep in keys:
                by_sid.setdefault(sid, []).append("%s r%s" % (arm, rep))
            parts.append("%s to %s (%s)" % (
                _PREV_OUTCOME_WORDS.get(a, str(a)) if a is not None else "no execution",
                _PREV_OUTCOME_WORDS.get(b, str(b)),
                "; ".join("%s %s" % (sid, ", ".join(v)) for sid, v in by_sid.items())))
        n = sum(len(v) for v in changes.values())
        text = "%d execution%s changed outcome, from the previous run's to this one's: %s" % (
            n, "" if n == 1 else "s", "; ".join(parts))
        # the executions that became inconclusive, placed against the host's memory episodes (Timing 6)
        became = [clines[k] for (a, b), keys in changes.items() if b == "inconclusive" for k in keys]
        eps = memory_episodes(run)
        if became and eps:
            inside = [l for l in became if any(l is x for ep in eps for x, _ in ep["rows"])]
            text += "; of the %d that became inconclusive, %s ran while the container stood at its memory limit " \
                    "(Timing 6)%s" % (len(became), facts._join(_exec_name(l) for l in inside) if inside else "none",
                                      ", the others outside that window" if inside and len(inside) < len(became)
                                      else "")
        items.append(text)
    else:
        items.append("no execution changed outcome")
    n_inc = sum(1 for l in run.rows if l.get("outcome") == "inconclusive")
    p_inc = sum(1 for l in prev.rows if l.get("outcome") == "inconclusive")
    items.append("%d inconclusive execution%s against %d" % (n_inc, "" if n_inc == 1 else "s", p_inc))
    out = "\\parhead{Since the previous edition}\n"
    out += ("Computed from the previous edition's run %s (manifest version %s, the %s set), given to this document as "
            "an input:\n" % (tt(prev.run_id()), esc(pm.get("manifest_version")),
                             "re-evaluated" if prev.reevaluated else "recorded"))
    out += "\\begin{itemize}\\setlength\\itemsep{0pt}\n" + "".join(
        "\\item %s.\n" % esc_wrap(t[0].upper() + t[1:], 2000) for t in items) + "\\end{itemize}\n"
    return out


def provenance_section(ctx):
    """I.4: what answered and what judged (prose, the compact table), what changed since the previous edition (when
    its run is given), the coverage (prose `coverage`, the per-arm outcome table) and the limitations and exemptions
    (prose `exemptions`, the counts with their appendix links)."""
    out = ctx.prose_paragraphs("provenance")
    out += provenance_table(ctx)
    out += since_previous(ctx)
    out += "\\parhead{Coverage}\n"
    out += ctx.prose_paragraphs("provenance", "coverage", optional=True)
    out += counting_units_table(ctx)
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
        if mc == "missing" and not cells.planned_on(run, sid, "mcp"):
            # no MCP execution because the plan runs the test on the direct access paths alone: not missing (M1)
            mc = "not_planned"
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
    n_ni = mcp.get("not_implemented_mcp", 0)
    mtext = "%d agree (fully checked %d%s), %s, %d by design (%d verified, %d unverified), %d not " \
            "comparable, %d disagree%s" % (
                mcp.get("agree", 0), full, fewer_repeats_text(fewer),
                cells.attn_marked("%d not implemented" % n_ni) if n_ni else "0 not implemented",
                mcp.get("service_by_design", 0), verified, unverified, mcp.get("not_comparable", 0),
                mcp.get("disagree", 0), (", %d missing" % mcp["missing"]) if mcp.get("missing") else "")
    if mcp.get("not_planned"):
        mtext += ", %d with no MCP execution planned" % mcp["not_planned"]
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
    out += "\\bline{chapter %s}{%d test%s · %d unit%s in %d famil%s}\n" % (
        esc(chapter), len(sids), "" if len(sids) == 1 else "s", n_units, "" if n_units == 1 else "s", n_families,
        "y" if n_families == 1 else "ies")
    out += "\\bline{standing}{tests · %s (one glyph per test and access path, the worst repeat; I.3)}\n" % " · ".join(
        "%s %s" % (esc(arm), cells.glyph_counts(states, arm)) for arm in ARMS)
    per_area = (run.summary().get("per_area") or {}).get(chapter) or {}
    if per_area:
        otext = []
        for arm in ARMS:
            o = per_area.get(arm) or {}
            items = [(cells.attn_marked if k == "not_implemented_mcp" else str)("%d %s" % (o[k], w))
                     for k, w in (("pass", "pass"), ("fail", "fail"), ("error", "error"),
                                  ("finding", "finding"), ("inconclusive", "inconclusive"),
                                  ("mixed", "mixed"), ("not_implemented_mcp", "not implemented"),
                                  ("service_by_design", "by design"), ("other", "other"))
                     if o.get(k)]
            otext.append("%s %s" % (esc(arm), esc(", ".join(items)) if items else "–"))
        # compare's per_area: one entry per test and access path, the recorded outcome over the test's repeats
        # (compare.aggregate_outcome: the common one, else the worst; mixed where the repeats differ otherwise), so the
        # row counts tests, as the standing row does, never executions (the verification of 2026-10-16, H1)
        out += "\\bline{outcomes}{tests · %s (each test's outcome over its repeats, as the comparison records it: the " \
               "common one, else the worst; mixed: its repeats differ otherwise; I.3)}\n" % " · ".join(otext)
    dtext, mtext = _agreement_counts(ctx, sids)
    out += "\\bline{agreement}{test pairs · direct: %s · mcp: %s}\n" % (dtext, mtext)
    out += "\\bline{reading}{%s}\n" % esc(summary_reading_text(ctx, chapter, sids))
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
        server.append(cells.fmt_exec(s.get("p50")) if s.get("n") else "–")
    backends = [l.get("mcp_phases", {}).get("backend_ms") for sid in sids for l in run.arm_lines(sid, "mcp")
                if isinstance(l.get("mcp_phases"), dict) and cells.host_cache_of(l) == "warm"]
    be = cells.percentile(backends, 50)
    out += "\\bline{timing}{host-cold wall p50 / p90: %s; host-warm wall p50 / p90: %s; server p50 (host-warm): %s" \
           "%s}\n" % (", ".join(cold) or "–", ", ".join(warm) or "–", " · ".join(server),
                      (" (mcp backend p50 %s)" % cells.fmt_exec(be)) if be is not None else "")
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


def summary_reading_text(ctx, chapter, sids=None):
    """One line under the chapter's summary rows (READ-102): both rows count tests, where they differ (the MCP
    agreement class in the standing, the recorded outcome in the outcomes row) with this chapter's own case or
    'none in this chapter', the pointer to I.3's full reading (defined once, not repeated per chapter: the
    verification of 2026-10-16, shared L5), the agreement count with the chapter's pairs short of fully checked
    (agreement_levels_text) and what the property list below counts; plain text."""
    if sids is None:
        sids = ctx.units.chapter_scenarios(chapter)
    example = standing_vs_outcome_reading(ctx, sids)
    text = ("standing and outcomes both count tests, one per access path; they differ only via MCP where the "
            "service decided the call or one repeat of it, the standing showing the agreement class and the outcomes "
            "row the recorded outcome over the repeats (%s); the full reading is I.3's. Agreement counts one "
            "comparison per test%s; "
            "the property list below counts property claims, one per check and test, over every access path" % (
                example, agreement_levels_text(ctx, chapter)))
    return text


def agreement_levels_text(ctx, chapter):
    """The agreement clause of the reading line (READ-110): the three levels an agreement can mean are I.3's; here
    the chapter's agreeing pairs short of fully checked, by what they compared (facts._mcp_agree_short: on the
    generic checks or a refusal alone, which verifies no property; with some property check compared and another
    not compared on some repeat), the first such pair named from the record (mcp_agree_partial_example)."""
    v = ctx.chapter_values(chapter)
    text = " (its three levels, identical answers, the same check outcomes, every property check compared on both: I.3)"
    if v.get("n_mcp_agree_partial"):
        text += "; %d agreeing pair%s short of fully checked: %s%s" % (
            v["n_mcp_agree_partial"], "" if v["n_mcp_agree_partial"] == 1 else "s",
            v.get("mcp_agree_short_text") or "see the detail blocks",
            ("; e.g. " + v["mcp_agree_partial_example"]) if v.get("mcp_agree_partial_example") else "")
    return text


CLAIM_STATUS_WORD = {"held": "every test", "part": "%(n)d of %(m)d tests",
                     "none": "no test passed on both (%(words)s there)",
                     "failed": "%(n)d of %(m)d tests, with a failure or a disagreement (I.1)"}


def claims_map(ctx, chapter):
    """The chapter's evidence map (READ-101): the families' claims as a compact table after the conclusion's
    headline, each claim with its family linked and what the direct access paths showed of it (the tests passing
    on both), the claims shown first; rendered when the chapter prose sets `evidence_map` true, else nothing (the
    prose then carries the claims inline, as before). The rows are the same records the conclusion's `claims_held`
    and `claims_unsupported` strings are built from (facts._claims)."""
    prose = ctx.prose.get(chapter) or {}
    if not prose.get("evidence_map"):
        return ""
    items = list((ctx.chapter_values(chapter).get("claims_items") or []))
    if not items:
        return ""
    order = {"held": 0, "part": 1, "none": 2, "failed": 3}
    items.sort(key=lambda it: order.get(it["status"], 9))
    rows = []
    for it in items:
        shown = "shown" if it["status"] in ("held", "part") else "not shown"
        rows.append(["\\hyperref[fam:%s]{\\texttt{%s}}" % (esc(it["family"]), esc(it["family"])),
                     inline(it["claim"]), "%s: %s" % (shown, esc(CLAIM_STATUS_WORD[it["status"]] % it))])
    return latex.plain_table(latex.spec_p([22, 108, 36]), ["family", "the claim (the property in words)",
                                                           "of the tests that carry the claim, those passing on both "
                                                           "direct access paths"], rows,
                             "the evidence map: what each family claims and how many of the tests that carry the claim "
                             "passed on both direct access paths (a family's other tests are not counted; the MCP access "
                             "path, the exceptions and the properties no path showed follow)")


EXCEPTION_GROUPS = (
    ("failure", "a failure, an error or a disagreement",
     "a check failed, no HTTP answer came, or two access paths disagree: a defect of the host, the service or the "
     "suite, listed in I.1 and quoted in the detail block"),
    ("finding", "a measured operational finding",
     "a stated bound was broken on the clock the check reads: a 200 past its budget (the host's own elapsed time), "
     "a walk stopped outside its stated rule (the seed's clock), "
     "or a 503 deadline delivered more than max(1 s, the budget) late on the client's clock, which carried no "
     "answer; the FINDING line of the detail block quotes the rule and the overrun"),
    ("unverified", "the property not shown on some access path",
     "an execution could not show the property there, for the recorded reason named after each test; the detail "
     "block says whether another access path or repeat showed it"),
    ("mcp_limit", "not comparable through the MCP tools",
     "no check was evaluated on both access paths (the MCP answer exposes none of the fields the checks read), so "
     "nothing is known about the property through the tools"),
    ("service_difference", "an expected service difference not fully verified",
     "the service answered in the host's place as its table documents, but one side of the pair did not match the "
     "case (appendix B states which); no property evidence through the tools"),
    ("direct_difference", "the two direct access paths answered differently",
     "direct-local and direct-remote gave different answers (a time limitation on one of them, or the same check "
     "outcomes on different bodies) without a check disagreeing"),
    ("missing", "a planned execution not recorded",
     "the plan had an execution no record holds: the standing there is unknown"),
)


def _finding_without_property(line, is_property):
    """Whether a finding execution evaluated its property checks without a pass (a late 503 deadline carried no
    answer: the property is not shown there either; the verification of 2026-10-16, chapters issue 4)."""
    props = [c for c in (line.get("checks") or []) if isinstance(c, dict) and is_property(c.get("check"))]
    return bool(props) and not any(c.get("outcome") == "pass" for c in props)


def exception_groups_of(states, is_property=None):
    """The exception groups (EXCEPTION_GROUPS keys) a footnoted test belongs to, from its arm states; with
    `is_property` (the run's), a finding execution whose property checks passed nothing counts under the property
    not shown too."""
    groups = []

    def add(key):
        if key not in groups:
            groups.append(key)
    for arm in ARMS:
        st = states[arm]
        if st.missing_planned:
            add("missing")
        if not st.present:
            continue
        if st.ne or st.worst in ("fail", "error", "other"):
            add("failure")
        if any(l.get("outcome") == "finding" for l in st.lines):
            add("finding")
            if is_property is not None and any(_finding_without_property(l, is_property) for l in st.lines
                                               if l.get("outcome") == "finding"):
                add("unverified")
        if any(l.get("outcome") == "inconclusive" for l in st.lines):
            add("unverified")
        if arm == "mcp" and "n" in st.glyph_marks:
            add("mcp_limit")
        if arm == "mcp" and "u" in st.glyph_marks:
            add("service_difference")
        if arm == "direct-remote" and any(m in st.glyph_marks for m in ("t", "s", "n")):
            add("direct_difference")
    return groups


def exception_index(ctx, chapter):
    """The chapter's footnoted tests grouped by consequence (READ-111): per group its title, the tests linked to
    their detail blocks, how many of them block the deploy gate, and the practical evidence gap in one sentence;
    empty when the chapter has no footnoted test. The timings and the raw check texts stay in the detail blocks."""
    fn = ctx.footnotes.get(chapter)
    ids = list(fn.order) if fn else []
    if not ids:
        return ""
    by_group = {}
    for sid in ids:
        for g in exception_groups_of(ctx.state(sid), ctx.run.is_property):
            by_group.setdefault(g, []).append(sid)
    blocking = set(facts._gate_blocking(ctx, ids))
    gated = isinstance(ctx.run.gate, dict)
    out = ("\\noindent{\\small The %d footnoted test%s of this chapter, by what each exception means (a test may "
           "appear under more than one heading; the detail blocks after the tables hold the timings and the raw "
           "checks):}\\par\\nopagebreak\n\\begin{itemize}\\setlength\\itemsep{0pt}\n" % (
               len(ids), "" if len(ids) == 1 else "s"))
    for key, title, gap in EXCEPTION_GROUPS:
        sids = by_group.get(key)
        if not sids:
            continue
        n_block = sum(1 for x in sids if x in blocking)
        gate_text_ = ("%d of them among the tests that block the deploy gate" % n_block
                      if gated else "the gate was not judged")
        if key == "unverified":
            # each test with its own recorded reason per access path, never a menu of possible causes (the
            # external review of the v6 document, V6-302)
            listed = "; ".join("%s: %s: %s" % (esc(arms), esc(why), ", ".join(latex.sid(x, detail=True) for x in ids))
                               for (arms, why), ids in facts.unshown_index_groups(ctx.run, sids, ctx.state))
            out += "\\item \\textbf{%s} (%d test%s; %s): %s — %s\n" % (
                esc(title), len(sids), "" if len(sids) == 1 else "s", esc(gate_text_), esc(gap), listed)
            continue
        out += "\\item \\textbf{%s} (%d test%s; %s): %s — %s\n" % (
            esc(title), len(sids), "" if len(sids) == 1 else "s", esc(gate_text_),
            ", ".join(latex.sid(x, detail=True) for x in sids), esc(gap))
    out += "\\end{itemize}\n"
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
    lines = blocks.class_lines(st)  # a by-design arm's class and status from the case's executions (VERIFY-FINAL R1)
    cls, evidence = blocks.mcp_class(lines, run.scenario(sid), st.agreement, run.cases)
    glyph = st.glyph_tex(_cell_fn(st, footnote))
    note = None
    if cls == "none":
        return "%s\\ no execution" % glyph, None
    if cls in ("not_expressible", "answer_lacks_field"):
        case = st.agreement.get("case")
        reason = blocks.short_reason(run, sid, case, ctx.short_reasons())
        return "%s\\ %s\\newline %s; see %s" % (glyph, blocks.class_word_tex(cls, st), reason,
                                                 blocks.limit_appendix(run, sid, st)), None
    first = next((l for l in lines if l.get("http") is not None), None)
    code = (first.get("error") or {}).get("code") if first and isinstance(first.get("error"), dict) else None
    status = "http %s%s" % (esc(first.get("http")), (" " + tt(code)) if code else "") if first else "–"
    if cls in ("refused_before_host", "service_answered", "relayed", "host_state", "refused_unexplained"):
        text = "%s\\ %s\\newline %s" % (glyph, blocks.class_word_tex(cls, st), status)
        return text, None
    if cls == "normalised":
        case = next((l.get("service_by_design_case") for l in st.lines if l.get("service_by_design_case")), None) \
            or st.agreement.get("case")
        text = "%s\\ %s\\newline %s" % (glyph, blocks.class_word_tex(cls, st), status)
        if with_times and st.present and st.glyph_key not in ("not_implemented_mcp", "service_by_design"):
            text = _cell_tex(st, footnote, "list")
        if case:
            note = "{\\footnotesize mcp, %s: %s; not a service refusal (case %s, appendix B)}" % (
                latex.breakable_id(sid), blocks.short_reason(run, sid, case, ctx.short_reasons()), tt(case))
        return text, note
    # answered by the host through the service: the usual cell
    if with_times:
        return _cell_tex(st, footnote, "list"), None
    return "%s\\ %s\\newline %s" % (glyph, blocks.class_word_tex(cls, st), status), None


EDGE_LEGEND = ("mcp column: the recorded status and the service's code; the class in words follows from whether "
               "the host was reached, the error's source and the by-design case (appendix B): refused before host "
               "execution · not expressible through MCP · answered by the service in the host's place · answered "
               "after service normalization · host refusal relayed through the service · answered by the host "
               "through the service")


def _edge_legend_extra(ctx, members):
    """The class words beyond EDGE_LEGEND's that a member's MCP cell prints (an answer-side scenario-level mark: not
    implemented yet via MCP, the answer lacking a checked field), so that the legend names every class shown."""
    run = ctx.run
    extra = []
    for sid in members:
        st = ctx.state(sid)["mcp"]
        cls, _ = blocks.mcp_class(blocks.class_lines(st), run.scenario(sid), st.agreement, run.cases)
        if cls == "answer_lacks_field" and blocks.MCP_CLASS_WORD[cls] not in extra:
            extra.append(blocks.MCP_CLASS_WORD[cls])
    return "".join(" · " + w for w in extra)


def _unit_blocks(ctx, chapter, uid, members, unit_note=False):
    """The detail blocks of a unit's footnoted members. The chapter prose's `details[<scenario id>]` may carry a
    reader-facing `title` in place of the manifest's and a `note` (the block's `reading` line), both filled on
    the scenario's own facts (`x[<id>]`, `manifest_version`) so that they can be conditional (READ-317, READ-319,
    READ-320)."""
    out = ""
    prop, _ = _unit_prose(ctx, chapter, uid)
    # the unit's property line states what the unit's members request together: a block's `requested` reads it only
    # for a unit of one test, else the scenario's own title and checks (the external review of the v6 document,
    # V6-302: a block never describes its neighbours' requests)
    property_text = ctx.fill(prop) if prop and len(members) == 1 else None
    details = (ctx.prose.get(chapter) or {}).get("details") or {}
    for sid in members:
        n = ctx.footnotes[chapter].get(sid)
        if n is None:
            continue
        d = details.get(sid) if isinstance(details.get(sid), dict) else {}
        values = {"x": {sid: facts.scenario_facts(ctx, sid)}, "id": sid,
                  "manifest_version": (ctx.run.manifest or {}).get("manifest_version")} if d else None
        out += blocks.detail_block(ctx.run, sid, ctx.state(sid), n, ctx.units, ctx.short_reasons(), unit_note,
                                   property_text=property_text,
                                   title=ctx.fill(d["title"], values) if d.get("title") else None,
                                   note=ctx.fill(d["note"], values) if d.get("note") else None)
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
            unrun = cells.relation_unrun(run, rel, arm, o)
            if unrun:
                # not judged: a scenario of the relation is not run via MCP (not "a raw answer is missing")
                cs.append("\\gnone {\\scriptsize %s}" % cells.relation_unrun_text(unrun, short=True))
                continue
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
        reasons, where = [], []
        for sid in members:
            st = ctx.state(sid)["mcp"]
            case = st.agreement.get("case") or next((l.get("service_by_design_case") for l in st.lines
                                                     if l.get("service_by_design_case")), None)
            r = blocks.short_reason(run, sid, case, ctx.short_reasons())
            if r not in reasons:
                reasons.append(r)
            a = blocks.limit_appendix(run, sid, st)
            if a not in where:
                where.append(a)
        body = "%s: %s; see %s." % (cells.attn("Not exercised via MCP"), "; ".join(reasons),
                                    " and ".join(sorted(where)) or "B")
    return "\\noindent{\\sffamily\\footnotesize\\gni\\ %s}\\par\\smallskip\n" % body


def unit_section(ctx, chapter, uid):
    run, units = ctx.run, ctx.units
    u = units.units[uid]
    members = list(u["members"])
    prop, note = _unit_prose(ctx, chapter, uid)
    prop_tex = inline(ctx.fill(prop)) if prop else latex.placeholder("(unit property: to be written)")
    # the note is filled on the unit's own values too (facts.unit_facts: its members under `x`), so that it can
    # state a member's gate relation from the policy (READ-318: N-usecase-defect); a note whose conditional leaves
    # nothing for this run (P-edge's v6-only note) prints nothing
    note_text = ctx.fill(note, facts.unit_facts(ctx, uid)).strip() if note else ""
    if note_text:
        prop_tex += "\\newline{\\color{black!60}\\itshape %s}" % inline(note_text)
    layout = u["unit_type"]
    key = cells.LAYOUT_KEY.get(layout)
    legend = cells.TIMING_LEGEND.get(layout)
    # a pairs table of repeated tests (chapter D) prints every repeat's wall: its key and legend name the repeats in
    # execution order and the MCP cold call as the record has it (READ-321)
    repeats = max([run.scenario(sid).get("repeat") for sid in members if isinstance(run.scenario(sid).get("repeat"), int)]
                  or [1])
    if layout == "pairs" and repeats > 1:
        key, legend = cells.repeats_texts(repeats, _mcp_cold_first(run, members))
    out = latex.unit_header(uid, esc(u.get("unit_title") or ""), u["unit_type"], prop_tex, key=key)
    fns = dict((sid, _footnotes_for(ctx, chapter, sid)) for sid in members)
    meta = units.by_scenario
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
        out += latex.replicates_table("variants" if not groups else esc(latex.cut_text(u.get("unit_title") or uid, 60)),
                                      rows, ids, legend)
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
            expected, prev = blocks.refusal_expected(run.scenario(sid), prev, run.is_property)
            cs = [st[arm].glyph_tex(_cell_fn(st[arm], fns[sid])) for arm in ("direct-local", "direct-remote")]
            mcp_cell, note = _mcp_edge_cell(ctx, sid, st["mcp"], fns[sid], with_times=False)
            rows.append(([latex.sid(sid, detail=True) if fns[sid] else latex.breakable_id(sid),
                          esc_wrap(_variant(meta[sid])), expected] + cs + [mcp_cell], note))
        out += latex.refusals_table(rows, units_ids(ctx, members), EDGE_LEGEND + _edge_legend_extra(ctx, members)
                                    + (("; " + rk) if rk else ""))
        states = [ctx.state(sid) for sid in members]
        # the unit's line (the wall max per arm) is labelled as the unit's where the unit shares its family's id,
        # so that it is not read as a repeat of the family's line with p50 / p90 (VERIFY-T L1, VERIFY-PLXM L3)
        label = tt(uid) + (" (the refusals unit)" if u.get("family") == uid else "")
        out += "\\sumline{%s: %d test%s · %s · %s}\n" % (label, len(members), "" if len(members) == 1 else "s", " · ".join(
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
                                + _edge_legend_extra(ctx, members) + (("; " + rk) if rk else ""))
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
                variant = (meta[sid] or {}).get("variant")
                second = latex.sid(src, units.unit_of(src)) if (not variant or variant == src) else \
                    "%s · %s" % (latex.sid(src, units.unit_of(src)), second)
            rows.append([latex.sid(sid, detail=True) if fns[sid] else latex.breakable_id(sid), second] + cs)
        rows += _relation_rows(ctx, u.get("relations") or [])
        out += latex.pairs_table(rows, units_ids(ctx, members), "test", legend)
    note = mcp_unit_note(ctx, chapter, uid, members)
    out += note
    out += _unit_blocks(ctx, chapter, uid, members, unit_note=bool(note))
    return out


def _mcp_cold_first(run, members):
    """Whether, over the members' MCP executions, repeat 0 is the host-cold call and every later repeat host-warm
    (the paired loop's order), so that a repeated pairs table's legend may say so; False when no MCP execution has
    a host state or the record says otherwise."""
    seen = False
    for sid in members:
        for l in run.arm_lines(sid, "mcp"):
            hc = cells.host_cache_of(l)
            if hc not in ("cold", "warm"):
                continue
            seen = True
            if (l.get("repeat") == 0) != (hc == "cold"):
                return False
    return seen


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


def intro_table(ctx, table, values):
    """One comparison table of a chapter's Background (the prose slot `intro_tables`, READ-209: the three notions
    of path support side by side): {after, caption, columns, rows, widths}; the cells are prose (inline markup,
    placeholders and conditionals), `after` the index of the Background paragraph the table follows (default: the
    last)."""
    cols = [str(c) for c in (table.get("columns") or [])]
    rows = [[inline(ctx.fill(str(c), values)) for c in r] for r in (table.get("rows") or []) if isinstance(r, list)]
    if not cols or not rows:
        return ""
    widths = table.get("widths") if isinstance(table.get("widths"), list) and len(table.get("widths")) == len(cols) \
        else [100.0 / len(cols)] * len(cols)
    caption = inline(ctx.fill(table.get("caption"), values)) if table.get("caption") else None
    return latex.plain_table(latex.spec_p(widths), [esc(c) for c in cols], rows, caption, needspace=len(rows) + 4)


def background(ctx, prose, values):
    """The chapter's Background: the `intro` items (paragraphs, or a {"table": ...} item, through `prose_blocks`:
    READ-311), each optionally followed by the `intro_tables` entries whose `after` names it (a table without
    `after` follows the last item; READ-209)."""
    intro = prose.get("intro") or []
    if isinstance(intro, str):
        intro = [intro]
    if not intro:
        return latex.placeholder("(chapter prose: to be written)") + "\n\n"
    tables = [t for t in (prose.get("intro_tables") or []) if isinstance(t, dict)]
    out = ""
    for i, p in enumerate(intro):
        out += prose_blocks(ctx, [p], values)
        for t in tables:
            after = t.get("after")
            if after == i or (after is None and i == len(intro) - 1):
                out += intro_table(ctx, t, values)
    return out


def toc_attention(ctx, ch):
    """The chapter's attention tag in the table of contents (the owner's decision of 2026-10-16): a red line under
    the chapter's own entry, linked to the chapter, with the tests needing attention per kind, the nonzero kinds
    alone in facts.ATTENTION_KINDS' order ('2 time-budget issues, 9 not implemented via MCP'); nothing when the
    chapter has none. The counts are facts.attention_facts', the ones I.1's attention table prints."""
    text = facts.attention_facts(ctx, ctx.units.chapter_scenarios(ch))["attention_text"]
    if not text:
        return ""
    return "\\addtocontents{toc}{\\protect\\tocattention{%s}{%s}}\n" % (latex.label_key("ch:%s" % ch), esc(text))


def attention_table(ctx):
    """I.1's attention table (the owner's decision of 2026-10-16): per chapter of the run the tests needing
    attention by kind (facts.attention_facts: a test counted once per kind and possibly under several), the
    nonzero counts in red, the kinds no chapter has left out, and the totals; empty when no test needs attention.
    The kinds are defined once, in appendix F's vocabulary (attention_legend_sentence; I.3 until V6-003)."""
    per, totals = [], {}
    for ch in ctx.chapters:
        v = facts.attention_facts(ctx, ctx.units.chapter_scenarios(ch))
        per.append((ch, v))
        for k, _, _ in facts.ATTENTION_KINDS:
            totals[k] = totals.get(k, 0) + v["n_attention_%s" % k]
        totals["tests"] = totals.get("tests", 0) + v["n_attention_tests"]
    kinds = [(k, many) for k, _, many in facts.ATTENTION_KINDS if totals.get(k)]
    if not kinds:
        return ""

    def cell(n):
        return cells.attn(str(n)) if n else "–"
    rows = []
    for ch, v in per:
        title = (ctx.prose.get(ch) or {}).get("title") or ctx.units.titles.get(ch, ch)
        rows.append(["\\hyperref[ch:%s]{%s\\ %s}" % (esc(ch), esc(ch), esc(title))]
                    + [cell(v["n_attention_%s" % k]) for k, _ in kinds]
                    + [cell(v["n_attention_tests"]), str(len(ctx.units.chapter_scenarios(ch)))])
    rows.append(["\\textbf{all chapters}"] + ["\\textbf{%s}" % cell(totals[k]) for k, _ in kinds]
                + ["\\textbf{%s}" % cell(totals["tests"]),
                   "\\textbf{%d}" % sum(len(ctx.units.chapter_scenarios(ch)) for ch in ctx.chapters)])
    header = ["chapter"] + [cells.attn(esc(many)) for _, many in kinds] + [cells.attn("tests needing attention"),
                                                                          "tests"]
    widths = [34] + [18] * len(kinds) + [20, 10]
    spec = latex.spec_p(widths, aligns=["raggedright"] + ["raggedleft"] * (len(kinds) + 2))
    # the tests counted under not implemented via MCP by a property-level mark alone (their MCP access path does
    # not stand ○), named so that the reader finds them (VERIFY-MERGE L1)
    prop = ["%s %s" % (esc(ch), cells.attn("%d (%s)" % (v["n_attention_property"], esc(facts.ids_ranges(
        [x.strip() for x in v["attention_property_ids"].split(",") if x.strip()])))))
            for ch, v in per if v["n_attention_property"]]
    prop_text = ("; of the not implemented via MCP, by a property-level mark alone (the MCP access path not standing "
                 "\\gni, the cell carrying the mark): %s" % "; ".join(prop)) if prop else ""
    if totals.get("time_budget"):
        # the kind glossed where a reader of I.1 first meets it in red (RV9)
        prop_text += ("; time-budget issues: tests with an execution that hit or stated a time budget, or whose wall "
                      "printed above %s; they include intended stops and slow passes (appendix F)" % esc(
                          fmt_t(cells.TIME_TINT_RED_MS)))
        # the chapters whose time-budget counts include executions in the host's memory episodes (RV1)
        mem = memory_summary(ctx)
        chs = []
        for m in mem:
            for ch in m["tb_chapters"]:
                if ch not in chs:
                    chs.append(ch)
        if chs:
            prop_text += ("; the time-budget counts of %s include executions that ran while the engine's container "
                          "stood at its memory limit, %s UTC (Timing 6)" % (
                              facts._join(chs), " and ".join(m["window"] for m in mem)))
    return latex.plain_table(spec, header, rows,
                             "needs attention: the tests of each chapter by kind, red where any (a test counts once "
                             "per kind and may count under several; the kinds are defined in appendix F's vocabulary, and the "
                             "contents print each chapter's counts under its line)%s" % prop_text,
                             needspace=min(16, len(rows) + 4))


MEMORY_LIMIT_SHARE = 0.98     # a sample at 98 % of the container's memory limit or more is "at its limit"
MEMORY_EPISODE_MIN_SAMPLES = 10
MEMORY_GAP_S = 60.0           # two samples at the limit further apart than this (the sampler paused) end an episode
MEMORY_DROP_SHARE = 0.25      # a fall of a quarter of the limit or more from one sample to the next ends an episode
                              # as a drop, reported as such
MEMORY_LATE_MS = 30000.0      # an MCP host call longer than the host's own search by more than this ...
MEMORY_LATE_SHARE = 0.90      # ... while the container was at 90 % of its limit or more: "the host answered late"
MEMORY_LIST_MAX = 12


def _utc_s(text):
    """Epoch seconds of an ISO UTC time, its fraction kept (reader.sc_time: the classification is at millisecond
    precision; the v7 baseline's M-planted-0013 started 0.8 s before the window's last sample, the same second)."""
    return reader.sc_time(text)


def _clock(t):
    import time as _time
    return _time.strftime("%H:%M:%S", _time.gmtime(t)) if t is not None else "–"


def _gib(b):
    return b / float(1 << 30)


def _mem_samples(run):
    """[(t, sample)] of the samples with a memory and a limit, in time order."""
    out = []
    for s in getattr(run, "samples", None) or []:
        if not isinstance(s, dict) or cells._num(s.get("mem_bytes")) is None or not cells._num(s.get("mem_limit_bytes")):
            continue
        t = _utc_s(s.get("t_utc"))
        if t is not None:
            out.append((t, s))
    out.sort(key=lambda x: x[0])
    return out


def _mem_at(ss, t):
    """The sample at or last before |t| (the first sample when |t| precedes them all), None without samples."""
    best = None
    for st, s in ss:
        if st > t:
            break
        best = s
    return best if best is not None else (ss[0][1] if ss else None)


def _mem_peak(ss, a, b):
    """The highest memory (bytes) over [a, b], the sample at |a| included."""
    vals = [s["mem_bytes"] for t, s in ss if a <= t <= b]
    first = _mem_at(ss, a)
    if first is not None:
        vals.append(first["mem_bytes"])
    return max(vals) if vals else None


def _span(line):
    s = _utc_s(line.get("started_utc"))
    e = _utc_s(line.get("ended_utc"))
    if s is None:
        return None, None
    return s, (e if e is not None and e >= s else s)


def _late_host(run, ss):
    """The MCP executions in which the host answered the service late: the service's host call (mcp_phases
    backend_ms) longer than the host's own search (backend_pattern_ms) by more than MEMORY_LATE_MS, the container at
    MEMORY_LATE_SHARE of its limit or more at the execution's start; each (line, share, host call start, end)."""
    out = []
    for l in run.rows:
        p = l.get("mcp_phases") if isinstance(l.get("mcp_phases"), dict) else {}
        b, h = cells._num(p.get("backend_ms")), cells._num(p.get("backend_pattern_ms"))
        s, _ = _span(l)
        if b is None or h is None or s is None or b - h <= MEMORY_LATE_MS:
            continue
        m = _mem_at(ss, s)
        share = m["mem_bytes"] / float(m["mem_limit_bytes"]) if m else None
        if share is None or share < MEMORY_LATE_SHARE:
            continue
        call = s + ((cells._num(p.get("submit_ms")) or 0) + (cells._num(p.get("queue_ms")) or 0)) / 1000.0
        out.append((l, share, call, call + b / 1000.0))
    return out


def memory_episodes(run):
    """The windows in which the engine's container stood at its memory limit (samples.jsonl: mem_bytes at
    MEMORY_LIMIT_SHARE of mem_limit_bytes or more), each a contiguous run of such samples (a sample below the share,
    or a gap of more than MEMORY_GAP_S, ends it) of at least MEMORY_EPISODE_MIN_SAMPLES samples; [] when none (every
    run before the v7 baseline's 02:41-03:05Z at 128 GiB). Each: its start and end (seconds, fraction kept), samples,
    limit and peak, the disk read meanwhile (GiB), a drop at its end (a fall of MEMORY_DROP_SHARE of the limit or more
    to the next sample) with the evidence against a restart, the memory ten minutes after, and every execution whose
    [start, end] overlaps it, placed: 'into' (started before it), 'inside', 'across' (ran past its end)."""
    ss = _mem_samples(run)
    runs_, cur = [], []
    for i, (t, s) in enumerate(ss):
        high = s["mem_bytes"] >= MEMORY_LIMIT_SHARE * s["mem_limit_bytes"]
        if high and cur and t - cur[-1][1] > MEMORY_GAP_S:
            runs_.append(cur)
            cur = []
        if high:
            cur.append((i, t, s))
        elif cur:
            runs_.append(cur)
            cur = []
    if cur:
        runs_.append(cur)
    instance = (getattr(run, "run", None) or {}).get("server_instance")
    out = []
    for high in runs_:
        if len(high) < MEMORY_EPISODE_MIN_SAMPLES:
            continue
        (i0, t0, s0), (i1, t1, s1) = high[0], high[-1]
        limit = s0["mem_limit_bytes"]
        ep = {"t0": t0, "t1": t1, "start": _clock(t0), "end": _clock(t1), "n": len(high),
              "container": s0.get("container") or "the engine's container", "limit_gib": _gib(limit),
              "peak_gib": _gib(max(s["mem_bytes"] for _, _, s in high)), "drop": None, "after_gib": None}
        blk = [cells._num(s.get("blk_read")) for _, _, s in high if cells._num(s.get("blk_read")) is not None]
        ep["read_gib"] = _gib(blk[-1] - blk[0]) if len(blk) > 1 and blk[-1] >= blk[0] else None
        after = [s for t, s in ss if t >= t1 + 600]
        ep["after_gib"] = _gib(after[0]["mem_bytes"]) if after else None
        rows = []
        for l in run.rows:
            a, b = _span(l)
            if a is None or a > t1 or b < t0:
                continue
            rows.append((l, "into" if a < t0 else "across" if b > t1 else "inside"))
        ep["rows"] = rows
        if i1 + 1 < len(ss):
            tn, sn = ss[i1 + 1]
            if s1["mem_bytes"] - sn["mem_bytes"] >= MEMORY_DROP_SHARE * limit:
                b1, bn = cells._num(s1.get("blk_read")), cells._num(sn.get("blk_read"))
                flight = [l for l, _ in rows if _span(l)[0] <= t1 and _span(l)[1] >= tn]
                # no restart: one server instance for the whole run, the disk-read counter kept counting, and every
                # request in flight across the fall completed with an answer the checks judged
                ep["drop"] = {"from_gib": _gib(s1["mem_bytes"]), "to_gib": _gib(sn["mem_bytes"]), "s": tn - t1,
                              "at": tn, "in_flight": flight,
                              "no_restart": bool(isinstance(instance, str) and instance and b1 is not None
                                                 and bn is not None and bn >= b1 and flight
                                                 and all(l.get("outcome") not in ("error",) and l.get("http") == 200
                                                         for l in flight))}
        out.append(ep)
    return out


MEMORY_SLOWER_FACTOR = 2.0    # an execution in an episode at twice its test's time elsewhere or more is "slower"


def memory_summary(ctx):
    """Per memory episode the facts the summary lines read (I.1, Timing 1, the attention caption): the window
    (HH:MM to HH:MM), the limit, the executions that ran in it, their chapters, those that did not pass, the chapters
    of the time-budget tests among them, and whether some passing search in it took MEMORY_SLOWER_FACTOR times its
    test's time elsewhere or more (this run outside every episode, or the previous edition's run when given)."""
    run = ctx.run
    eps = memory_episodes(run)
    if not eps:
        return []
    chapters = list(getattr(ctx.units, "chapters", None) or [])
    def area(l):
        return str(run.scenario(l.get("scenario")).get("area") or l.get("area") or "?")
    def ch_sort(xs):
        return sorted(set(xs), key=lambda x: (chapters.index(x) if x in chapters else 99, x))
    in_any = set(id(l) for ep in eps for l, _ in ep["rows"])
    prev = getattr(run, "previous", None)
    out = []
    for ep in eps:
        rows = [l for l, _ in ep["rows"]]
        slower = False
        for l in rows:
            sv = _server_ms(l)
            if l.get("outcome") != "pass" or sv is None:
                continue
            ref = [_server_ms(x) for x in run.rows if x.get("scenario") == l.get("scenario") and id(x) not in in_any]
            ref += [_server_ms(x) for x in (prev.rows if prev is not None else []) if x.get("scenario") == l.get("scenario")]
            ref = [x for x in ref if x is not None and x > 0]
            if ref and sv >= MEMORY_SLOWER_FACTOR * max(ref):
                slower = True
                break
        out.append({"window": "%s to %s" % (ep["start"][:5], ep["end"][:5]), "limit_gib": ep["limit_gib"],
                    "rows": rows, "chapters": ch_sort(area(l) for l in rows),
                    "bad": [l for l in rows if l.get("outcome") != "pass"],
                    "tb_chapters": ch_sort(area(l) for l in rows if facts.time_budget_issue(l)), "slower": slower})
    return out


def _memory_outcomes_text(m):
    """'30 passed, and P-warm-0004 direct-local r0 ended in a 503 deadline (inconclusive)' of an episode's
    executions."""
    bad = m["bad"]
    if not bad:
        return "every one passed"
    words = "; ".join("%s %s" % (_exec_name(l), ("ended in a 503 deadline, %s" % _stop_word(l).rsplit(", ", 1)[-1])
                                 if l.get("http") == 503 else "was %s" % _stop_word(l)) for l in bad[:MEMORY_LIST_MAX])
    return "%d passed, and %s" % (len(m["rows"]) - len(bad), words)


def memory_lead_text(ctx):
    """I.1's sentence on the host's memory episodes, after the result (the reader review of the v7 document, RV1;
    the fact verifier, V9): when, the executions that ran then and their outcomes, whether searches were slower, and
    what it is; '' without an episode."""
    mem = memory_summary(ctx)
    if not mem:
        return ""
    parts = []
    for m in mem:
        parts.append("from %s UTC the engine's container stood at its %s GiB memory limit on the host: %d execution%s "
                     "ran then (chapter%s %s)%s; %s" % (
                         m["window"], fmt_n(m["limit_gib"]), len(m["rows"]), "" if len(m["rows"]) == 1 else "s",
                         "" if len(m["chapters"]) == 1 else "s", facts._join(m["chapters"]),
                         ", and searches that ran then were slower than the same tests elsewhere" if m["slower"] else "",
                         _memory_outcomes_text(m)))
    text = "; ".join(parts)
    return ("*The host's memory.* %s%s (Timing 6). It is an operational note about the host, not a finding of the "
            "service." % (text[0].upper(), text[1:]))


def memory_timing_bullet(ctx):
    """Timing 1's observation on the host's memory episodes (RV1, V9); '' without one."""
    mem = memory_summary(ctx)
    if not mem:
        return ""
    return ("The engine's container stood at its %s GiB memory limit from %s UTC; the executions that ran then are "
            "placed in Timing 6: %s." % (fmt_n(mem[0]["limit_gib"]), " and ".join(m["window"] for m in mem), "; ".join(
                "%s: %d, %s" % (m["window"], len(m["rows"]), _memory_outcomes_text(m)) if len(mem) > 1 else
                "%d, of which %s" % (len(m["rows"]), _memory_outcomes_text(m)) for m in mem)))


def memory_chapter_note(ctx, ch):
    """A chapter's line on the host's memory for its Timing paragraph (RV1, V9): how many of its executions ran
    while the container stood at its memory limit, and its MCP executions the host answered late; '' when none."""
    run = ctx.run
    eps = memory_episodes(run)
    if not eps:
        return ""
    sids = set(ctx.units.chapter_scenarios(ch))
    parts = []
    for ep in eps:
        mine = [l for l, _ in ep["rows"] if l.get("scenario") in sids]
        if mine:
            parts.append("%d of this chapter's executions ran while the engine's container stood at its memory limit, "
                         "%s to %s UTC" % (len(mine), ep["start"][:5], ep["end"][:5]))
    for l, share, _, _ in _late_host(run, _mem_samples(run)):
        if l.get("scenario") in sids:
            parts.append("the host answered %s late, with the container at %d %% of its memory limit" % (
                _exec_name(l), int(share * 100)))
    if not parts:
        return ""
    text = "; ".join(parts)
    return "%s%s (Timing 6)" % (text[0].upper(), text[1:])


def memory_episode(run):
    """The first memory episode (memory_episodes), None without one."""
    eps = memory_episodes(run)
    return eps[0] if eps else None


def _secs(ms):
    return fmt_t(ms) if ms is not None else "–"


def _server_ms(l):
    return cells._num((l.get("server") or {}).get("elapsed_ms")) if isinstance(l.get("server"), dict) else None


def _range_text(vals, unit_word="server time"):
    vals = [v for v in vals if v is not None]
    if not vals:
        return None
    lo, hi = _secs(min(vals)), _secs(max(vals))
    return ("%s %s" % (lo, unit_word)) if lo == hi else ("%s to %s %s" % (lo, hi, unit_word))


def _exec_name(l):
    return "%s %s r%s" % (l.get("scenario"), l.get("arm"), l.get("repeat"))


def _stop_word(l):
    """How an execution that did not pass ended, in a few words: a 503 deadline (with its wall), the host's 503
    deadline the service relayed, else its outcome."""
    outcome = blocks.OUTCOME_WORD.get(l.get("outcome"), str(l.get("outcome"))) if hasattr(blocks, "OUTCOME_WORD") \
        else str(l.get("outcome"))
    outcome = re.sub(r"\\[a-z]+\{([^}]*)\}", r"\1", outcome)
    if l.get("http") == 503:
        w = cells._num(l.get("wall_ms"))
        return "a 503 deadline at %s s wall, %s" % (fmt_n(w / 1000.0, 1) if w is not None else "–", outcome)
    if "503 deadline" in str(l.get("time_limited") or ""):
        return "the host's 503 deadline relayed by the service, %s" % (
            "service-side by design" if l.get("outcome") == "service_by_design" else outcome)
    return outcome


def _kind_words(l, many):
    """How a group of executions ended, plural when |many|: 503 deadlines (their outcome), the host's 503 deadline
    the service relayed (service-side by design), else the outcome."""
    if l.get("http") == 503:
        return "%s, %s" % ("503 deadlines" if many else "a 503 deadline", _stop_word(l).rsplit(", ", 1)[-1])
    if "503 deadline" in str(l.get("time_limited") or ""):
        return "the host's 503 deadline relayed by the service, %s" % _stop_word(l).rsplit(", ", 1)[-1]
    return _stop_word(l)


def _arms_repeats(ls):
    """'direct-local r0, r1 and direct-remote r0' of executions."""
    by_arm = {}
    for l in ls:
        by_arm.setdefault(l.get("arm"), []).append("r%s" % l.get("repeat"))
    return facts._join("%s %s" % (arm, ", ".join(reps)) for arm, reps in by_arm.items())


def _count_word(n):
    return "one" if n == 1 else "%d" % n


def _not_passed(l):
    """The executions the memory reading places against the windows (RV2): every inconclusive one and every 503
    deadline, the host's relayed by the service included."""
    return (l.get("outcome") == "inconclusive" or l.get("http") == 503
            or "503 deadline" in str(l.get("time_limited") or ""))


def _previous_servers(ctx, sid, arms=None):
    """The server times of a test in the previous edition's run (make_doc --previous), on |arms| when given."""
    prev = getattr(ctx.run, "previous", None)
    if prev is None:
        return None
    return [_server_ms(l) for l in prev.rows if l.get("scenario") == sid and (arms is None or l.get("arm") in arms)
            and _server_ms(l) is not None]


def _previous_outcome(ctx, l):
    prev = getattr(ctx.run, "previous", None)
    if prev is None:
        return None
    p = prev.lines().get((l.get("scenario"), l.get("arm"), l.get("repeat")))
    return (p or {}).get("outcome")


def _service_log_text(ctx, named):
    """The service's own log as the identity input states it (make_doc --identity: `service_log` {text, times,
    source}; reader._load_identity refuses times outside the run): the text with its times, each time tied to the
    host call of a named MCP execution it falls in (|named|: [(line, call start, call end)]); '' without one."""
    ident = getattr(ctx.run, "identity", None)
    log = ident.get("service_log") if isinstance(ident, dict) and isinstance(ident.get("service_log"), dict) else None
    if not log or not log.get("text"):
        return ""
    times = [t for t in (_utc_s(x) for x in (log.get("times") or [])) if t is not None]
    parts = []
    for t in times:
        hit = [l for l, a, b in named if a <= t <= b]
        parts.append("at %s UTC%s" % (_clock(t), (", during %s's call to the host" % _exec_name(hit[0])) if hit else ""))
    return ("The service's own log, checked separately and not part of the run's files, recorded %s%s." % (
        str(log["text"]).rstrip("."), (": " + "; ".join(parts)) if parts else ""))


def memory_episode_text(ctx):
    """Timing 6's reading of the host's memory (the coordinator's addendum and the reviews of the v7 document: RV2,
    RV8, V1, V2, V4): led by what it is (an operational note about the host), then the samples (each episode, a drop
    at its end), what ran then (every execution in an episode with its outcomes, the time-budget tests among them,
    the MCP executions the host answered late, the service's log, the other time-budget tests), the run's 503
    deadlines and inconclusive executions placed against the episodes, and the limits; nothing without an episode."""
    run = ctx.run
    eps = memory_episodes(run)
    if not eps:
        return ""
    ss = _mem_samples(run)
    several = len(eps) > 1
    obs = []
    for ep in eps:
        text = ("the engine's container (%s) stood at its memory limit of %s GiB from %s to %s UTC (%d samples at %d %% "
                "of it or more, peak %s GiB)%s" % (
                    ep["container"], fmt_n(ep["limit_gib"]), ep["start"], ep["end"], ep["n"],
                    int(MEMORY_LIMIT_SHARE * 100), fmt_n(ep["peak_gib"], 1),
                    ("; it read %s GiB from disk meanwhile" % fmt_n(ep["read_gib"])) if ep["read_gib"] else ""))
        d = ep["drop"]
        if d:
            text += (". At the end of the window its memory fell from %s to %s GiB within %s s%s" % (
                fmt_n(d["from_gib"]), fmt_n(d["to_gib"]), fmt_n(d["s"]),
                " (cached index pages released, no restart: the request in flight completed)" if d["no_restart"]
                else " (whether the container restarted is not established)"))
        if ep["after_gib"] is not None:
            text += "; it was at %s GiB ten minutes after the window" % fmt_n(ep["after_gib"], 1)
        obs.append(text)
    first = _gib(ss[0][1]["mem_bytes"]) if ss else None
    observation = "*Observation.* %s.%s" % (". In a second window, ".join(o[0].upper() + o[1:] if i == 0 else o
                                                                         for i, o in enumerate(obs)),
                                           (" It started the run at %s GiB." % fmt_n(first)) if first else "")
    tags = dict((sid, (run.scenario(sid).get("tags") or [])) for sid in run.order)
    tb_rows = [l for l in run.rows if facts.time_budget_issue(l)]
    tb_tests = []
    for l in tb_rows:
        if l.get("scenario") not in tb_tests:
            tb_tests.append(l.get("scenario"))
    late = _late_host(run, ss)
    named = [(l, a, b) for l, _, a, b in late]
    ran, late_par, rest_par, place_par = [], [], [], []
    window_tests = []
    for ep in eps:
        rows = ep["rows"]
        areas = sorted(set(str(run.scenario(l.get("scenario")).get("area") or l.get("area") or "?") for l, _ in rows),
                       key=lambda x: (ctx.units.chapters.index(x) if x in ctx.units.chapters else 99, x))
        bad = [l for l, _ in rows if l.get("outcome") != "pass"]
        lead = "In the window of %s to %s, " % (ep["start"], ep["end"]) if several else ""
        text = "%s%d execution%s ran in the window (chapter%s %s)" % (
            lead, len(rows), "" if len(rows) == 1 else "s", "" if len(areas) == 1 else "s", facts._join(areas))
        into = [l for l, pos in rows if pos == "into"]
        across = [l for l, pos in rows if pos == "across"]
        if into or across:
            text = text[:-1] + "; %s)" % " and ".join(x for x in (
                ("%s started before it" % _count_word(len(into))) if into else "",
                ("%s ran past its end" % _count_word(len(across))) if across else "") if x)
        if not bad:
            text += ": every one passed"
        else:
            items = []
            for l in bad[:MEMORY_LIST_MAX]:
                prev = _previous_outcome(ctx, l)
                items.append("%s (%s%s)" % (_exec_name(l), _stop_word(l), "; it passed in the previous edition's run"
                                            if prev == "pass" else ""))
            text += ": %d passed, and %s did not: %s" % (
                len(rows) - len(bad), "%d" % len(bad) if len(bad) > 1 else "one",
                "; ".join(items) + ("; … %d more" % (len(bad) - MEMORY_LIST_MAX) if len(bad) > MEMORY_LIST_MAX else ""))
        ran.append(text + ".")
        # the time-budget tests with an execution in the window, placed by time before any tag (V1)
        tests = []
        for l, _ in rows:
            sid = l.get("scenario")
            if sid in tb_tests and sid not in tests and sid not in window_tests:
                tests.append(sid)
        chapters = list(getattr(ctx.units, "chapters", None) or [])
        tests.sort(key=lambda sid: (chapters.index(run.scenario(sid).get("area"))
                                    if run.scenario(sid).get("area") in chapters else 99, sid))
        window_tests += tests
        if tests:
            parts = []
            for sid in tests:
                mine = [(l, pos) for l, pos in rows if l.get("scenario") == sid]
                ls = [l for l, _ in mine]
                within = [l for l, pos in mine if pos != "across"]
                passed = [l for l in within if l.get("outcome") == "pass"]
                bits = []
                if "budget_stress" in tags.get(sid, []):
                    bits.append("budget_stress")
                served = _range_text([_server_ms(l) for l in passed])
                if served:
                    bits.append("%s on %d execution%s" % (served, len(passed), "" if len(passed) == 1 else "s"))
                for l in within:
                    if l.get("outcome") != "pass":
                        bits.append("%s r%s: %s" % (l.get("arm"), l.get("repeat"), "the 503 deadline above"
                                                    if l.get("http") == 503 else "named above"))
                for l, pos in mine:
                    if pos != "across":
                        continue
                    a = _span(l)[0]
                    txt = "%s r%s started %s s before the window's end and ran across %s" % (
                        l.get("arm"), l.get("repeat"), fmt_n(ep["t1"] - a, 1),
                        ("the fall to %s GiB" % fmt_n(ep["drop"]["to_gib"])) if ep["drop"] else "its end")
                    if l.get("outcome") == "pass" and _server_ms(l) is not None:
                        txt += ": %s server time" % _secs(_server_ms(l))
                    later = sorted((x for x in run.rows if x.get("scenario") == sid and _span(x)[0] is not None
                                    and _span(x)[0] > a), key=lambda x: _span(x)[0])
                    if later and _server_ms(later[0]) is not None:
                        txt += "; its next execution, %s r%s at %s UTC, took %s" % (
                            later[0].get("arm"), later[0].get("repeat"), _clock(_span(later[0])[0]),
                            _secs(_server_ms(later[0])))
                    bits.append(txt)
                outside = [_server_ms(l) for l in run.rows if l.get("scenario") == sid and not any(l is x for x in ls)]
                o_text = _range_text(outside)
                if o_text and not any(pos == "across" for _, pos in mine):
                    bits.append("outside the window %s" % o_text.replace(" server time", ""))
                prev = _previous_servers(ctx, sid, set(l.get("arm") for l in ls))
                p_text = _range_text(prev or [])
                if p_text:
                    bits.append("in the previous edition's run %s" % p_text.replace(" server time", ""))
                parts.append("%s (%s)" % (sid, "; ".join(bits)))
            ran.append("Of the %d tests with a time-budget issue (I.1), %d had an execution in %s: %s." % (
                len(tb_tests), len(tests), "this window" if several else "the window", "; ".join(parts)))
    if late:
        items = []
        for l, share, a, b in late:
            p = l.get("mcp_phases") or {}
            items.append("%s at %s UTC, the container at %d %% of its limit (a host search of %s, %s from the "
                         "service's call to the host's answer)" % (
                             _exec_name(l), _clock(_span(l)[0]), int(share * 100),
                             _secs(cells._num(p.get("backend_pattern_ms"))), _secs(cells._num(p.get("backend_ms")))))
        late_par.append("The host answered the service late (its call to the host longer than the host's own search by "
                   "more than %s s, with the container at %d %% of its memory limit or more) in %s: %s." % (
                       fmt_n(MEMORY_LATE_MS / 1000.0), int(MEMORY_LATE_SHARE * 100),
                       "one MCP execution" if len(late) == 1 else "%d MCP executions" % len(late), "; ".join(items)))
    log = _service_log_text(ctx, named)
    if log:
        late_par.append(log)
    late_tests = []
    for l, _, _, _ in late:
        sid = l.get("scenario")
        if sid in tb_tests and sid not in window_tests and sid not in late_tests:
            late_tests.append(sid)
    rest_ids = [sid for sid in tb_tests if sid not in window_tests and sid not in late_tests]
    stress = [sid for sid in rest_ids if "budget_stress" in tags.get(sid, [])]
    stops, other = [], []
    for sid in rest_ids:
        if sid in stress:
            continue
        ls = [l for l in tb_rows if l.get("scenario") == sid]
        (stops if all(str(l.get("time_limited") or "") for l in ls) else other).append(sid)
    rest = []
    if stress:
        rest.append("%d budget_stress test%s ran outside the window%s" % (
            len(stress), "" if len(stress) == 1 else "s", "s" if several else ""))
    if stops:
        ls = [l for l in tb_rows if l.get("scenario") in stops]
        rest.append("%d ended by a time budget they state (%s)%s" % (
            len(stops), facts.ids_ranges(stops), ", as their tests intend or allow: each of these executions passed"
            if all(l.get("outcome") == "pass" for l in ls) else ": %s" % ", ".join(sorted(set(
                str(l.get("outcome")) for l in ls)))))
    if other:
        rest.append("%d otherwise (%s)" % (len(other), facts.ids_ranges(other)))
    if rest:
        rest_par.append("Of the other %d tests with a time-budget issue, %s." % (len(rest_ids), "; ".join(rest)))
    # every 503 deadline and inconclusive execution placed against the windows (RV2)
    np_rows = [l for l in run.rows if _not_passed(l)]
    if np_rows:
        def in_ep(l):
            return any(l is x for ep in eps for x, _ in ep["rows"])
        inside = [l for l in np_rows if in_ep(l)]
        outside = [l for l in np_rows if not in_ep(l)]
        text = "The run's 503 deadlines and inconclusive executions, placed against the window%s: " % (
            "s" if several else "")
        text += ("inside, %s (above)" % facts._join(_exec_name(l) for l in inside)) if inside else "none inside"
        if outside:
            groups = []
            by_sid = {}
            for l in outside:
                by_sid.setdefault(l.get("scenario"), []).append(l)
            for sid, ls in by_sid.items():
                kinds = {}
                for l in ls:
                    kinds.setdefault(_kind_words(l, False), []).append(l)
                spans = [_span(l) for l in ls]
                a, b = min(x for x, _ in spans), max(y for _, y in spans)
                peak = _mem_peak(ss, a, b)
                groups.append("%s (%s; %s to %s UTC%s)" % (
                    sid, "; ".join("%s: %s" % (_arms_repeats(v), _kind_words(v[0], len(v) > 1)) for v in kinds.values()),
                    _clock(a), _clock(b), (", memory at most %s GiB" % fmt_n(_gib(peak), 1)) if peak is not None
                    else ""))
            text += "; outside, and not attributed to %s: %s" % ("them" if several else "it", "; ".join(groups))
        place_par.append(text + ".")
    has_prev = getattr(run, "previous", None) is not None
    limits = ("*Limits.* The samples show how much memory the container used, not what used it or why: the run does "
              "not isolate the cause. That a search in the window was slower is read against the same test's other "
              "executions%s, not measured." % (" and the previous edition's run" if has_prev else ""))
    lead = ("The container's memory is an operational note about the host, not a finding of the service.")
    return "\\parhead{The host's memory}\n" + paragraphs(
        [lead, observation, "*What ran then.* " + " ".join(ran)] + [" ".join(x) for x in (late_par, rest_par, place_par)
                                                                   if x] + [limits])


def chapter(ctx, ch):
    """A chapter: its conclusion first (the `summary` prose and the summary block), then what it leaves untested
    (the `coverage_gap` prose and the computed lines), then the background (`intro`) and the families (the owner's
    review: lead each chapter with its conclusion and coverage gap)."""
    units = ctx.units
    ctx.footnotes[ch] = cells.Footnotes()
    prose = ctx.prose.get(ch) or {}
    title = prose.get("title") or units.titles.get(ch, ch)
    out = latex.section("%s\\qd %s" % (esc(ch), esc(title)), "ch:%s" % ch, mark="%s · %s" % (ch, title))
    out += toc_attention(ctx, ch)
    # the chapter's purpose in the reader's terms (READ-101): what the operation lets a user ask and what kind of
    # evidence the chapter supplies, before the conclusion's inventory; optional prose (`purpose`)
    purpose = prose.get("purpose")
    fams = units.chapter_families(ch)
    sids = units.chapter_scenarios(ch)
    body = ""
    for fid in fams:
        body += family_section(ctx, ch, fid)
    n_units = sum(len(units.families[f]["units"]) for f in fams)
    values = ctx.chapter_values(ch)
    # the summary is built after the families so the footnotes are assigned
    if purpose:
        out += prose_blocks(ctx, purpose, values)
    out += "\\parhead{Conclusion}\n"
    summary = prose.get("summary")
    out += prose_blocks(ctx, summary, values) if summary else latex.placeholder("(chapter summary: to be written)") + "\n\n"
    out += claims_map(ctx, ch)
    out += chapter_summary(ctx, ch, sids, n_units, len(fams))
    out += "\\parhead{Coverage and remaining limits}\n"
    gap = prose.get("coverage_gap")
    out += prose_blocks(ctx, gap, values) if gap else latex.placeholder("(coverage gap: to be written)") + "\n\n"
    out += exception_index(ctx, ch)
    table = untested_table(ctx, ch)
    if table:
        # how the table counts, against the conclusion's 'n of m' (VERIFY-T L10: the two counts side by side)
        out += ("The properties below are counted over every access path: a test passes a property when any "
                "execution on any access path passed it, whereas the conclusion's n of m counts the tests passing on "
                "both direct access paths.\n")
        out += table
    else:
        out += "Every property check of this chapter was passed by at least one execution of each test.\n\n"
    out += "\\parhead{Background}\n"
    out += background(ctx, prose, values)
    # the compact record key right before the first family's tables (re-verification N-1)
    out += record_key_paragraph(ctx, ch)
    return out + body


def prose_blocks(ctx, items, values):
    """A chapter's prose slot (`purpose`, `summary`, `coverage_gap`, `intro`) as LaTeX: a string or a list whose
    items are paragraphs (strings, filled on the chapter's values) or one small table, written as
    {"table": {"caption", "columns", "rows", "widths"}} with every cell filled like prose (READ-311: a reader-facing
    map inside a chapter's background, where a paragraph of lists would hide what is compared)."""
    out = ""
    for item in ([items] if isinstance(items, str) else (items or [])):
        if isinstance(item, dict) and isinstance(item.get("table"), dict):
            t = item["table"]
            columns = [inline(ctx.fill(c, values)) for c in (t.get("columns") or [])]
            rows = [[inline(ctx.fill(c, values)) for c in row] for row in (t.get("rows") or [])]
            widths = t.get("widths") or [1] * len(columns)
            caption = inline(ctx.fill(t["caption"], values)) if t.get("caption") else None
            out += latex.plain_table(latex.spec_p(widths), columns, rows, caption)
        elif isinstance(item, str) and item:
            out += paragraphs([ctx.fill(item, values)])
    return out


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
    """[(phase, p10, p25, p50, p75, p90, n, n_zero, n_negative)] over phase_population: the quantiles of the
    positive values (a log scale holds neither a zero nor a negative value), n the positive ones, n_zero and
    n_negative the values left out (the external review of the v6 document, V6-006: v6's terminal lag has eight
    values below 0 ms and no zero)."""
    values = phase_population(ctx)
    out = []
    for k, name in PHASES:
        vals = [x for x in values[k] if x > 0]
        zeros = sum(1 for x in values[k] if x == 0)
        negative = sum(1 for x in values[k] if x < 0)
        if len(vals) < 2:
            continue
        out.append((name, cells.percentile(vals, 10), cells.percentile(vals, 25), cells.percentile(vals, 50),
                    cells.percentile(vals, 75), cells.percentile(vals, 90), len(vals), zeros, negative))
    return out


def phase_omitted(ctx):
    """[(phase, n, zeros, negative, min negative, max negative)] for each phase whose box leaves values out."""
    values = phase_population(ctx)
    out = []
    for k, name in PHASES:
        zeros = sum(1 for x in values[k] if x == 0)
        neg = [x for x in values[k] if x < 0]
        if zeros or neg:
            out.append((name, len(values[k]), zeros, len(neg), min(neg) if neg else None, max(neg) if neg else None))
    return out


def phase_omitted_text(ctx):
    """The caption's account of the values a log-scale box leaves out, per phase, from the plotted arrays: 'paging
    is 0 ms where the results came with the terminal status in 132 of 580 and terminal lag below 0 ms in 8 of 580
    (-10.6 to -0.8 ms)'; '' when no box leaves a value out."""
    parts = []
    for name, n, zeros, neg, lo, hi in phase_omitted(ctx):
        bits = []
        if zeros:
            bits.append("0 ms%s in %d of %d" % (" where the results came with the terminal status" if name == "paging"
                                                 else "", zeros, n))
        if neg:
            bits.append("below 0 ms in %d of %d (%s ms)" % (neg, n, fmt_n(lo, 1) if lo == hi else "%s to %s" % (
                fmt_n(lo, 1), fmt_n(hi, 1))))
        parts.append("%s is %s" % (name, " and ".join(bits)))
    return facts._join(parts)


def phase_omitted_short(ctx):
    """'paging without its 132 zeros, terminal lag without its 8 values below 0 ms' (Timing 3's method)."""
    parts = []
    for name, n, zeros, neg, lo, hi in phase_omitted(ctx):
        bits = (["its %d zero%s" % (zeros, "" if zeros == 1 else "s")] if zeros else []) + (
            ["its %d value%s below 0 ms" % (neg, "" if neg == 1 else "s")] if neg else [])
        parts.append("%s without %s" % (name, " and ".join(bits)))
    return ", ".join(parts)


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
    rows = ["phase,n,n_zero,p10_ms,p25_ms,p50_ms,p75_ms,p90_ms,n_negative"] + [
        "%s,%d,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%d" % (name.replace(" ", "_"), n, n_zero, p10, p25, p50, p75, p90, n_neg)
        for name, p10, p25, p50, p75, p90, n, n_zero, n_neg in phase_quantiles(ctx)]
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

    def stats(prefix, d, signed=False, tint=True):
        # a signed difference stays plain (fmt_delta); an execution time (a phase, a cold wall) carries its tint
        # as a marker for the prose (cells.fmt_exec_text); a connection set-up is no execution time (tint=False)
        d = d if isinstance(d, dict) and d.get("n") else None
        f = fmt_delta if signed else (cells.fmt_exec_text if tint else fmt_t)
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
        # the median paired ratio as a signed percentage off one (READ-011: '−4 %', never a ratio of two medians)
        v["cold_warm_%s_ratio_pct" % key] = fmt_pct_off_one((d.get("ratio") or {}).get("p50")) if d.get("n") else None
    fo = s.get("front_overhead") or {}
    for arm in ("direct-local", "direct-remote"):
        stats("connect_" + cells.ARM_SHORT[arm], (fo.get("connect_ms") or {}).get(arm), tint=False)
    nets = [(e.get("front_overhead") or {}).get("net_ms") for e in (run.comp.get("scenarios") or {}).values()
            if isinstance(e, dict) and isinstance(e.get("front_overhead"), dict)]
    nets = [x for x in nets if cells._num(x) is not None]
    v["front_net_p50"] = fmt_delta(cells.percentile(nets, 50)) if nets else None
    v["front_net_p90"] = fmt_delta(cells.percentile(nets, 90)) if nets else None
    v["front_net_n"] = str(len(nets)) if nets else None
    pairs = paired_walls(ctx)
    for key, idx in (("dl", 2), ("dr", 3), ("mcp", 4)):
        vals = [p[idx] for p in pairs if p[idx] is not None]
        v["warm_wall_%s_p50" % key] = cells.fmt_exec_text(cells.percentile(vals, 50)) if vals else None
        v["warm_wall_%s_p90" % key] = cells.fmt_exec_text(cells.percentile(vals, 90)) if vals else None
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
    v["paired_dl_wall_p50"] = cells.fmt_exec_text(dl_med) if both else None
    v["paired_mcp_wall_p50"] = cells.fmt_exec_text(mcp_med) if both else None
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
    v["mcp_phase_largest_p50"] = cells.fmt_exec_text(largest[0]) if largest else None
    # the flag of the observation's clause on terminal lag (READ-010: what the phase measures, its cause not isolated)
    v["mcp_phase_largest_is_terminal_lag"] = bool(largest) and largest[1] == "terminal lag"
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
    # every phase's values the log-scale box leaves out (zeros, values below 0 ms), from the plotted arrays (V6-006)
    v["phase_omitted_text"] = phase_omitted_text(ctx)   # '' (a false condition) when every value is positive
    v["phase_omitted_short"] = phase_omitted_short(ctx)
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
    ni = _OUTCOME_WORDS["not_implemented_mcp"]
    parts = ["%s (%s)" % ((cells.attn_marked if o == ni else str)("%d %s" % (n, o)), w)
             for n, o, w in cold_mcp_untimed_split(run)]
    return facts._join(parts) if parts else "none"


def fmt_pct_off_one(ratio):
    """A ratio as a signed whole percentage off one: 0.96 reads '−4 %', 1.05 '+5 %', 1.0 '0 %'; '–' for none."""
    v = cells._num(ratio)
    if v is None:
        return "–"
    pct = (v - 1.0) * 100.0
    rounded = int(abs(pct) + 0.5)
    if rounded == 0:
        return "0 %"
    return "%s%d %%" % ("−" if pct < 0 else "+", rounded)


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
    out += latex.section("%s 1\\qd Observations" % CLOSING_PREFIX, "sec:C1")
    intro = cp.get("intro") or []
    out += paragraphs([ctx.fill(p, closing_values(ctx)) for p in ([intro] if isinstance(intro, str) else intro)]) \
        if intro else latex.placeholder("(closing intro: to be written)") + "\n\n"
    obs = observations(ctx)
    # the host's memory episodes, where a reader of the closing section looks first (RV1, V9)
    mem_bullet = memory_timing_bullet(ctx)
    if mem_bullet:
        obs.append(inline(mem_bullet))
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
    out += latex.section("%s 2\\qd Paired latencies and the MCP phases" % CLOSING_PREFIX, "sec:C2")
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
    out += latex.section("%s 3\\qd MCP overhead decomposition" % CLOSING_PREFIX, "sec:C3")
    out += closing_section_prose(ctx, "C3", (
        "Pairs of host-warm executions only: the MCP warm execution (repeat 1 preferred) minus each direct access path's "
        "first warm execution; the phases as the service measured them. In these tables a / b never appears: "
        "p50, p90 and the maximum are columns.\n\n"))
    ph = s.get("mcp_phases") or {}
    rows = [[esc(k.replace("_ms", "").replace("_", " "))] + _stats_cells(ph.get(k)) for k in
            ("submit_ms", "queue_ms", "backend_ms", "backend_pattern_ms", "terminal_lag_ms", "paging_ms") if k in ph]
    out += latex.plain_table("@{}l r r r r@{}", ["phase", "p50", "p90", "max", "n"], rows, "MCP phases")
    mo = s.get("mcp_overhead") or {}
    rows = [[esc(k.replace("_ms", "").replace("_", " "))] + _stats_cells(mo.get(k), signed=True) for k in sorted(mo)]
    out += latex.plain_table("@{}l r r r r@{}", ["MCP minus direct", "p50", "p90", "max", "n"], rows,
                             "MCP overhead (warm MCP wall minus each direct access path; the backend time minus each direct access path); "
                             "signed")
    if s.get("mcp_cold_wall_ms"):
        out += latex.plain_table("@{}l r r r r@{}", ["", "p50", "p90", "max", "n"],
                                 [["MCP cold wall"] + _stats_cells(s["mcp_cold_wall_ms"])], "MCP cold wall")
    out += latex.section("%s 4\\qd Front cost" % CLOSING_PREFIX, "sec:C4")
    out += closing_section_prose(ctx, "C4", "")
    fo = s.get("front_overhead") or {}
    # a connection set-up is no execution time (wall_ms excludes it): its figures carry no tint
    rows = [[esc(arm)] + _stats_cells((fo.get("connect_ms") or {}).get(arm), tint=False)
            for arm in sorted(fo.get("connect_ms") or {})]
    out += latex.plain_table("@{}l r r r r@{}", ["connection set-up", "p50", "p90", "max", "n"], rows,
                             "connect\\_ms per access path")
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
    out += latex.plain_table("@{}l l r r r r@{}", ["access path", "", "p50", "p90", "max", "n"], rows,
                             "host-cold executions (every recorded host-cold execution with that time, whatever its outcome)")
    out += latex.section("%s 5\\qd Cold versus warm" % CLOSING_PREFIX, "sec:C5")
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
    out += latex.plain_table("@{}l r r r r r@{}", ["access path", "pairs", "delta p50", "delta p90", "ratio p50", "ratio p90"],
                             rows, "warm minus cold per scenario (delta, signed: negative when warm was faster), "
                                   "warm over cold (ratio)")
    out += latex.section("%s 6\\qd Load telemetry" % CLOSING_PREFIX, "sec:C6")
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
    # the peaks at the container's memory limit, read against the memory episodes (the reader review of the v7
    # document, RV2): which rows reached it, and when
    at_limit = {}
    eps = memory_episodes(ctx.run)
    if eps:
        limit = min(ep["limit_gib"] for ep in eps) * float(1 << 30)
        for area in areas:
            for arm in ARMS:
                peak = cells._num(((res.get(area) or {}).get(arm) or {}).get("peak_mem_bytes"))
                if peak is not None and peak >= MEMORY_LIMIT_SHARE * limit:
                    at_limit.setdefault(area, []).append(arm)
    peaks_text = ""
    if at_limit:
        peaks_text = "; the peaks at the memory limit (%s) are the window%s of %s UTC (The host's memory, below)" % (
            "; ".join("%s %s" % (esc(a), "on every access path" if len(v) == len(ARMS) else ", ".join(esc(x) for x in v))
                      for a, v in at_limit.items()), "s" if len(eps) > 1 else "",
            " and ".join("%s to %s" % (ep["start"][:5], ep["end"][:5]) for ep in eps))
    out += latex.plain_table("@{}l l r r r r r r@{}", ["chapter", "access path", "executions", "with samples", "peak CPU \\%",
                                                        "peak memory", "sra\\_all mean \\%", "load1 mean"], rows,
                             "sampled telemetry over each chapter's executions: the container's peak CPU (100\\,\\% is one "
                             "logical CPU) and peak memory (the container's memory usage as docker stats reports it), the host's sra\\_all CPU mean and load1 "
                             "mean; context of the shared host, not any request's own cost" + peaks_text)
    out += memory_episode_text(ctx)
    out += latex.section("%s 7\\qd Per-chapter timing" % CLOSING_PREFIX, "sec:C7")
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
    out += latex.plain_table("@{}l l l r r r r r r r r@{}", ["chapter", "access path", "host", "wall p50", "p90", "max",
                                                              "n", "server p50", "p90", "max", "n"], rows,
                             "wall and server time per chapter, access path and host-cache state, over every "
                             "recorded execution with that time (passes, allowed 503s, refusals and findings alike: "
                             "not the latency of successful searches); the server columns count only the answers "
                             "that carry timing.elapsed\\_ms, so their n is smaller and the columns are not paired")
    if wai:
        out += ("%d execution(s) ran warm after an interrupted execution and are counted in neither bucket.\n\n" % wai)
    return out


# -- appendices --------------------------------------------------------------------------------------------------

def _versioned(ctx, text):
    """A catalogue text with its manifest-version conditionals resolved (`[[manifest_version >= 7 ? a | b]]`: the v7
    delta), and nothing else filled: a text without one, or one carrying braces of its own, is returned as it is."""
    text = str(text)
    if "[[manifest_version" not in text or "{" in text:
        return text
    return facts.fill(text, {"manifest_version": (ctx.run.manifest or {}).get("manifest_version")}, unresolved=[])


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
            out += latex.catline(label, inline(_versioned(ctx, entry[key]))) + "\n"
    table = entry.get("decision_table")
    if isinstance(table, dict) and table.get("rows"):
        # a compact decision table (READ-014: the cases of budget_respected with the clock each uses), from the
        # catalogue prose: {caption, columns, rows}
        cols = [str(c) for c in (table.get("columns") or [])]
        rows = [[inline(str(c)) for c in r] for r in table["rows"] if isinstance(r, list)]
        widths = [100.0 / max(1, len(cols))] * len(cols)
        out += latex.plain_table(latex.spec_p(widths), [esc(c) for c in cols], rows,
                                 inline(table.get("caption") or ""), needspace=len(rows) + 4)
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
        out += "\\grouphead{The rules before any check}\n" + paragraphs([ctx.fill(p) for p in common])
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
    rules = ctx.verification_rules()
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
            # an entry is a check name (every item of it) or {check, detail_matches} (an item whose detail the
            # regex matches: suite_common.entry_explains)
            service.append("explains %s" % ", ".join(
                esc(e) if isinstance(e, str) else "%s (detail matching %s)" % (
                    esc(e.get("check")), tt(e.get("detail_matches"), 80)) for e in c["explains_checks"]
                if isinstance(e, (str, dict))))
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
                    reason = blocks.unverified_reason(ma, ctx.state(sid)["mcp"])
                    if reason:
                        unverified_detail.append("%s: %s" % (sid, reason))
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
                entry = "%s: %s" % (esc(mark), ", ".join(latex.sid(s, ctx.units.unit_of(s)) for s in sids_))
            else:
                entry = "%s: %d (%s, …)" % (esc(mark), len(sids_),
                                            ", ".join(latex.sid(s, ctx.units.unit_of(s)) for s in sids_[:3]))
            # the class word and count of the tests not implemented via MCP in red, as the case's class tag in the
            # first column (the owner's decision of 2026-10-16; VERIFY-MERGE M1 b)
            named.append(cells.attn(entry) if mark == "not implemented mcp" else entry)
        if not marks:
            # a case of the chapter's area that no scenario named and no answer matched (READ-017): documented,
            # not exercised in this run
            named.append("not exercised in this run")
        agreement = esc(c.get("agreement") or "–")
        if c.get("agreement") == "not_implemented_mcp":
            agreement = cells.attn(agreement)  # red wherever printed (the owner's decision of 2026-10-16)
        first = [tt(cid) + "\\newline{\\scriptsize %s}" % agreement]
        matches = latex.join_units(esc_wrap(c.get("scenario_match"), 600))
        if short.get(cid):
            matches += "\\newline{\\scriptsize in short: %s}" % inline(short[cid])
        # the pinned table's statement qualified for the deployment this run met (cases.json `current`, by manifest
        # version; the share review of the v7 document, V7-101): the table's own text and identity stay as they are
        notes = ctx.current_notes() if hasattr(ctx, "current_notes") else {}
        current = _versioned(ctx, notes.get(cid) or "").strip()
        if current:
            matches += "\\newline{\\scriptsize %s}" % inline(current)
        if rules.get(cid):
            matches += "\\newline{\\scriptsize verification: %s}" % inline(rules[cid])
        rows.append(first + [matches, "; ".join(service) or "– (none stated)", "; ".join(direct) or "– (none stated)",
                             "; ".join(named)])
        if unverified_detail:
            rows[-1] = (rows[-1], "{\\footnotesize %s}" % unverified_note(unverified_detail))
    out = ("One row per case of the service-by-design table that a scenario of the built chapters names, that an "
           "MCP answer matched, or that belongs to their areas: the case (its stable id, with the agreement class "
           "the table gives it), what it matches (with the short reason the tables use and, where the case has one, "
           "its verification rule), the service's documented answer and the direct expectation (both the table's "
           "statements, not observations), and this run's scenarios with their evidence.\n\n")
    out += ("\\noindent{\\small\\textbf{The evidence states of the last column.} \\emph{named}: the scenario's "
            "manifest entry names the case as its expected MCP behaviour; by itself no execution evidence. "
            "\\emph{matched by the answer}: an MCP answer of this run matched the case (observed on the MCP side). "
            "\\emph{verified}: the MCP answer matched the case's service side and the direct answer its host side "
            "(observed on both). \\emph{unverified}: observed on the MCP side, but the direct side did not match "
            "the case or another MCP repeat outside the case has no direct answer to compare with; the reason is "
            "under the row, one line per distinct reason. \\emph{not implemented mcp} and the other class words: "
            "the MCP agreement class of a scenario that names the case without the case deciding its standing. "
            "\\emph{not exercised in this run}: no scenario of the built chapters names the case and no answer "
            "matched it; the row documents the service's behaviour, which this run did not observe. A row proves "
            "the behaviour occurred only through a verified, unverified or matched state.}\\par\\medskip\n")
    out += latex.plain_table(latex.spec_p([30, 38, 40, 36, 30]),
                             ["case", "matches", "service (the table's answer)", "direct (the table's expectation)",
                              "this run's scenarios and their evidence"], rows)
    return out


def unverified_note(details):
    """The note under an unverified case's row (READ-018): one line per distinct reason with the tests it holds,
    the expectation and the observation apart and the consequence stated once; `details` are the 'sid: reason'
    (or 'sid (direct side not verified): reason') strings the table collected."""
    by_reason = {}
    order = []
    for d in details:
        sid, _, reason = str(d).partition(": ")
        reason = reason or str(d)
        if reason not in by_reason:
            by_reason[reason] = []
            order.append(reason)
        by_reason[reason].append(sid.replace(" (direct side not verified)", ""))
    lines = []
    for reason in order:
        sids = by_reason[reason]
        lines.append("%d test%s (%s): %s" % (len(sids), "" if len(sids) == 1 else "s",
                                             latex.join_units(esc_wrap(", ".join(sids))),
                                             latex.join_units(esc_wrap(reason))))
    return ("unverified pair%s, by reason — %s. Consequence: the service side matched the case, the pair as a whole "
            "is not verified; this is a gap in the comparison's evidence, not in itself a service defect." % (
                "" if len(details) == 1 else "s", "; ".join(lines)))


def relation_paths_text(paths):
    """A same_answer relation's compared paths in full and readable (the external review of the v6 document, V6-007:
    no argument is cut): '12 paths, a's against b's: release; server_instance against attempts.server_instance;
    ...', a path equal on both sides named once, the whole answer for an empty path, and a prefix every b path
    shares with its a path stated once ('34 paths, each against pattern.<the same path>: annotation, ...')."""
    pairs = [(str(a), str(b)) for a, b in paths]
    n = len(pairs)
    head = "%d path%s" % (n, "" if n == 1 else "s")
    prefix = None
    if n > 1 and all(a and b.endswith(a) and len(b) > len(a) for a, b in pairs):
        cands = set(b[:len(b) - len(a)] for a, b in pairs)
        prefix = cands.pop() if len(cands) == 1 else None
    if prefix:
        return "%s, each of a against %s<the same path> of b: %s" % (head, prefix, ", ".join(a for a, _ in pairs))
    words = []
    for a, b in pairs:
        a_, b_ = a or "the whole answer", b or "the whole answer"
        words.append(a_ if a == b else "%s against %s" % (a_, b_))
    return "%s, a's against b's: %s" % (head, "; ".join(words))


def relation_args_text(args):
    """A relation's arguments in full, as appendix C prints them under its title: the scalars, then a same_answer's
    paths (relation_paths_text) and any other list as JSON, never cut."""
    text = ", ".join("%s %s" % (k, v) for k, v in sorted(args.items()) if not isinstance(v, (list, dict)))
    for k, v in sorted(args.items()):
        if not isinstance(v, (list, dict)):
            continue
        if k == "paths" and isinstance(v, list) and all(isinstance(p, list) and len(p) == 2 for p in v):
            part = relation_paths_text(v)
        else:
            part = "%s %s" % (k, blocks._runs(v))
        text += ("; " if text else "") + part
    return text


def relations_appendix(ctx):
    run = ctx.run
    comp_rel = run.comp_relations()
    selected = set(ctx.selected_scenarios())
    rows = []
    for rel in run.relations:
        if rel.get("a") not in selected and rel.get("b") not in selected:
            continue
        crel = comp_rel.get(rel.get("id")) or {}
        arg_text = relation_args_text(rel.get("args") or {})
        per_arm = []
        for arm in ARMS:
            a = (crel.get("arms") or {}).get(arm) or {}
            unrun = cells.relation_unrun(run, rel, arm, a.get("outcome")) if a else []
            if not a:
                per_arm.append("\\gnone")
            elif unrun:
                per_arm.append("\\gnone {\\scriptsize %s}" % cells.relation_unrun_text(unrun, short=True))
            else:
                per_arm.append("%s {\\scriptsize %s}" % (cells.GLYPH.get(a.get("outcome"), "\\ginc"),
                                                        esc(a.get("detail"), 90)))
        rows.append([latex.breakable_id(rel.get("id")) + "\\newline{\\scriptsize %s}" % esc(rel.get("type")),
                     "%s %s %s" % (latex.sid(rel.get("a"), ctx.units.unit_of(rel.get("a"))),
                                   RELATION_SIGN.get(rel.get("type"), ""),
                                   latex.sid(rel.get("b"), ctx.units.unit_of(rel.get("b")))),
                     esc(rel.get("title"), 240) + (("\\newline{\\scriptsize %s}" % esc_wrap(arg_text)) if arg_text else "")]
                    + per_arm)
    out = appendix_lead(ctx, "appendix_C") or (
        "Every manifest relation among the built chapters' tests, with the comparator's outcome and detail per "
        "access path.\n\n")
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
    # the glossary at a readable size (the readability review, READ-020: two pages rather than 7.5 pt): it asks for
    # a dozen lines so that its heading is not left at a page's bottom
    out += "\\Needspace{12\\baselineskip}\n\\parhead{Glossary}\n"
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
        # two columns at footnote size (READ-020: readable at ordinary page size, two pages if need be)
        out += ("\\begin{multicols}{2}\\footnotesize\\raggedright\\setlength{\\columnsep}{5mm}"
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


def classes_tables(answer_side=False):
    """The vocabulary of appendix F, rendered from the shared definitions (doc/vocabulary.py): the three layers of
    an outcome (the check item, the execution, the displayed standing), then the two agreement classes;
    `answer_side`: the run has an answer-side scenario-level mark (the not-implemented meaning widened)."""
    out = ("\\noindent{\\small %s}\\par\\smallskip\n" % inline(vocabulary.CLASSIFICATIONS))
    out += ("\\noindent{\\small An outcome has three layers. \\textbf{A check item} (one check of one execution) "
            "answers one of five words; \\textbf{the execution's outcome} follows from its items by precedence, "
            "or is decided by the runner without them; \\textbf{the cell's standing} is what the tables show. "
            "The reading guide (I.3) and the glossary (appendix D) print these same definitions.}\\par\\smallskip\n")
    spec = latex.spec_p([30, 144])
    out += latex.plain_table(spec, ["check item", "meaning"], [[tt(k), inline(v)] for k, v in vocabulary.CHECK_OUTCOMES])
    out += latex.plain_table(latex.spec_p([30, 32, 78, 34]), ["execution outcome", "decided", "meaning", "consequence"],
                             [[tt(k), esc(d), inline(m), inline(c)] for k, d, m, c in vocabulary.outcome_table_rows(
                                 answer_side)])
    out += "\\noindent{\\small %s}\\par\\smallskip\n" % inline(vocabulary.STANDING)
    out += latex.plain_table(spec, ["direct agreement class", "meaning"],
                             [[tt(k), esc(v)] for k, v in vocabulary.DIRECT_CLASSES])
    out += latex.plain_table(spec, ["MCP agreement class", "meaning"],
                             [[tt(k), esc(v)] for k, v in vocabulary.MCP_CLASSES])
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
            _date(table.get("generated_utc")), tt(table.get("sha256") or "–"), tt(latex.public_paths(table.get("path"))))
        if table.get("replayed"):
            recorded = ((run.reeval or {}).get("by_design") or {}).get("path")
            text += "; the table the re-evaluation replayed (reevaluated.json by\\_design), found by its sha256%s" % (
                (" at this path (the re-evaluation recorded %s)" % tt(latex.public_paths(recorded)))
                if recorded and recorded != table.get("path") else " at the recorded path")
        elif table.get("replayed_sha256"):
            text += "; \\textbf{WARNING}: not the table the re-evaluation replayed (sha256 %s, %s cases): %s" % (
                tt(table["replayed_sha256"]), esc(table.get("replayed_cases", "–")),
                esc(latex.public_paths("; ".join(table.get("warnings") or []) or "the replayed table was not found"), 400))
        elif table.get("warnings"):
            text += "; " + esc(latex.public_paths("; ".join(table["warnings"])), 400)
        rows.append(("by-design table", text))
    elif table.get("warnings"):
        rows.append(("by-design table", "none read: %s" % esc(latex.public_paths("; ".join(table["warnings"])), 400)))
    mcp_run = run.mcp_run or {}
    for arm, info in sorted((r.get("arms") or {}).items()):
        info = info or {}
        ended = _arm_end(run, arm, info)
        scripts = info.get("scripts_commit") or (mcp_run.get("git_head") if arm == "mcp" else None)
        state = info.get("state")
        if not state and arm == "mcp" and info.get("exit") is not None and (info.get("finished_utc") or ended):
            # the MCP runner records its exit status, not a state (RV6)
            state = "finished, exit %s" % info.get("exit")
        rows.append(("access path %s" % esc(arm), "%s; state %s; %s executions recorded; runner %s, scripts %s; %s to %s" % (
            tt(info.get("base") or "–"), esc(state or "–"), esc(info.get("executions_recorded", "–")),
            esc(info.get("runner_host") or (mcp_runner_text(run) if arm == "mcp" else "–")), tt(_short(scripts, 12)),
            _date(info.get("started_utc")), _date(ended))))
    if mcp_run:
        # three revisions are distinct (READ-021): the runner's (the client that asked), the by-design table's service
        # commit (the documented behaviour the comparison read) and the service's deployed revision (the code that
        # served the MCP requests); the last is printed only where the run recorded it
        rows.append(("MCP runner", "git head %s (the client that asked, not the service's revision); base %s; admin "
                                   "key %s; %s" % (
                                       tt(_short(mcp_run.get("git_head"), 12)), tt(mcp_run.get("base") or "–"),
                                       esc(mcp_run.get("admin_key") or "–"), esc(service_version_text(run, full=True)))))
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
    rows.append(("comparison", comparison_lines_text(run)))
    if run.reruns:
        rows.append(("warm-ups run again for a resume", reruns_text(run)))
    hc = comp.get("host_cache") or {}
    rows.append(("host cache", "cold = the scenario's first execution on the host in this run; warm = another "
                               "within %s s; aged = older; the MCP pair runs cold then warm, the direct access paths "
                               "between them" % esc(hc.get("warm_window_s", "60"))))
    gt = manifest.get("ground_truth") or {}
    for key in sorted(gt):
        value = gt[key]
        if isinstance(value, (list, dict)):
            value = "%d entries" % len(value)
        rows.append((esc(key.replace("_", " ")), esc_wrap(latex.public_paths(value), 220)))
    return _kv_table(rows)


GATE_PROBLEMS_FILE = "gate_problems.json"
# the inputs this document read, each with its path in the suite's tree and its sha256, written beside doc.tex and
# doc.pdf (the external review of the v6 document, V6-004: a recipient of the PDF can retrieve and check every
# record behind its decision; deterministic, no clock)
EVIDENCE_FILE = "evidence.json"
EVIDENCE_RUN_FILES = ("manifest.json", "run.json", "run.reevaluated.json", "executions.jsonl",
                      "executions.reevaluated.jsonl", "compare.json", "compare.reevaluated.json", "reevaluated.json",
                      "plan.json", "orchestrator.json", "samples.jsonl", "gate.json", "gate.reevaluated.json",
                      "arms/direct-local/run.json", "arms/direct-remote/run.json", "arms/mcp/run.json")


def _file_sha256(path):
    try:
        return reader.sha256_of(path)
    except (OSError, TypeError):
        return None


def _identity_evidence(run, ident):
    """evidence.json's service identity: the recorded one where the runner recorded a commit, with an identity input
    beside it (its source and canonical sha256: e.g. the v7 baseline's service-log note); else the input's (the v6
    baseline: unchanged); else none."""
    def input_part(i):
        part = {"source": i.get("source"), "service_commit": i.get("service_commit"),
                "railway_deployment_id": i.get("railway_deployment_id")}
        if isinstance(i.get("deployed_service_commit"), dict):
            # the deployed commit the gate held the recorded identity to, with its source (RV3, V3)
            part["deployed_service_commit"] = i["deployed_service_commit"]
        if isinstance(i.get("service_log"), dict):
            # the service's log as the input states it, with its source (the document's body names neither the
            # log's file nor its internal states; RV8)
            part["service_log"] = i["service_log"]
        return dict(part, sha256_canonical=hashlib.sha256(json.dumps(i, sort_keys=True, separators=(",", ":")).encode(
            "utf-8")).hexdigest())
    rec = _recorded_identity_evidence(run)
    if rec and ident:
        rec = dict(rec, input=input_part(ident))
    return rec if rec else (input_part(ident) if ident else None)


def _recorded_identity_evidence(run):
    """The deployed service's identity as the MCP runner recorded it (manifest v7's runner: run.json
    service_identity), for evidence.json when no identity input was given; None when the run records no commit."""
    rec = recorded_service_identity(run) or {}
    if not rec.get("service_commits"):
        return None
    return {"source": "recorded by the MCP runner (run.json service_identity, read from the service before each "
                      "invocation)", "service_commit": rec.get("service_commit"),
            "service_commits": rec.get("service_commits"), "identity_changed": rec.get("identity_changed"),
            "invocations": rec.get("invocations")}


def evidence_manifest(ctx):
    """{document, run_directory, inputs, gate, by_design_table, service_identity} of this build (V6-004): every run
    file the document read (the set it read marked), the gate's record and policies, the declared file, the
    by-design table, the deployed service's identity as given; paths in the suite's tree (latex.public_paths)."""
    run = ctx.run
    inputs = []
    for name in EVIDENCE_RUN_FILES:
        path = os.path.join(run.dir, name)
        if os.path.exists(path):
            # the set this document did not read (the recorded one beside a re-evaluation, or the reverse), and a
            # recorded gate verdict the document judged again offline, are listed as present, not read
            other = (name in ("executions.jsonl", "compare.json", "run.json", "gate.json") and run.reevaluated) or \
                (name.endswith(".reevaluated.jsonl") or name.endswith(".reevaluated.json")) and not run.reevaluated
            if name.startswith("gate") and run.gate_offline:
                other = True
            inputs.append({"path": name, "sha256": _file_sha256(path), "role": "present, not read" if other
                           else "read by this document"})
    gate = run.gate if isinstance(run.gate, dict) else None
    pol = None
    if run.policy_path:
        pol = {"path": latex.public_paths(run.policy_path), "sha256_file": _file_sha256(run.policy_path),
               "sha256_canonical": ((gate or {}).get("policy") or {}).get("sha256")}
    also = [{"path": latex.public_paths(e.get("path")), "sha256_file": _file_sha256(e.get("path")),
             "verdict": None if not isinstance(e.get("verdict"), dict) else
             ("PASS" if e["verdict"].get("passed") else "FAIL")} for e in (getattr(run, "also_gates", None) or [])]
    dec = run.declared if isinstance(run.declared, dict) else None
    dec_path = os.path.join(reader.SUITE_DIR, "gate", dec["_file"]) if dec and dec.get("_file") else None
    table = run.table_meta or {}
    ident = run.identity if isinstance(getattr(run, "identity", None), dict) else None
    record = ctx.data_files.get(GATE_PROBLEMS_FILE)
    return {
        "path_bases": {"inputs": "the run directory (run_directory below; metagraph/scripts/suite/runs/ is "
                                 "git-ignored: the directory travels as a copy, its files checked by these sha256)",
                       "gate.record": "the document's build directory, beside doc.pdf and this file",
                       "other paths": "the repository root"},
        "document": {"run_id": run.run_id(), "set": "re-evaluated" if run.reevaluated else "recorded",
                     "generator": "metagraph/scripts/suite/doc/make_doc.py",
                     "generator_commit": reader.git_commit_of(HERE, os.path.dirname(HERE)),
                     "checks": checks_commit_text(run), "chapters": list(ctx.chapters)},
        "run_directory": latex.public_paths(run.dir),
        "inputs": inputs,
        "gate": {"record": GATE_PROBLEMS_FILE if gate else None,
                 "record_sha256": hashlib.sha256(record.encode("utf-8")).hexdigest() if record else None,
                 "verdict": None if not gate else ("GATE PASS" if gate.get("passed") else "GATE FAIL"),
                 "judged_offline_for_this_document": bool(run.gate_offline), "policy": pol, "also_policies": also,
                 "declared": None if not dec else {"path": latex.public_paths(dec_path) if dec_path else dec.get(
                     "_file"), "sha256_file": _file_sha256(dec_path) if dec_path else None,
                     # the policy pins the declared file by the sha256 of its canonical JSON (gate_policy.py)
                     "sha256_canonical": hashlib.sha256(sc.canonical_json(dict(
                         (k, v) for k, v in dec.items() if k != "_file"))).hexdigest()}},
        "by_design_table": {"path": latex.public_paths(table.get("path")), "sha256": table.get("sha256")}
        if table.get("path") else None,
        "service_identity": _identity_evidence(run, ident),
        # the previous edition's run, an input of 'Since the previous edition' (make_doc --previous)
        "previous_run": None if getattr(run, "previous", None) is None else {
            "run_id": run.previous.run_id(), "run_directory": latex.public_paths(run.previous.dir),
            "manifest_version": run.previous.manifest.get("manifest_version"),
            "manifest_content_sha256": run.previous.manifest.get("content_sha256"),
            "executions_file": run.previous.executions_file, "sha256": run.previous.executions_sha256},
    }


def gate_action_table(ctx, problems):
    """The gate's problems by group (GATE_GROUPS; READ-002, READ-023): per group its problems by kind, the scope they
    touch (the scenario-arm pairs, else the access paths or the run), what they imply and the evidence that would
    clear them (GATE_GROUP_ACTION)."""
    by_group = {}
    for p in problems:
        by_group.setdefault(gate_group_of(p), []).append(p)
    rows = []
    for key, title in GATE_GROUPS:
        ps = by_group.get(key) or []
        if not ps:
            continue
        kinds = {}
        for p in ps:
            kinds[p.get("kind") or "?"] = kinds.get(p.get("kind") or "?", 0) + 1
        kinds_text = ", ".join("%d %s" % (kinds[k], k) for k in gate_mod.KINDS if k in kinds) + "".join(
            ", %d %s" % (kinds[k], k) for k in sorted(k for k in kinds if k not in gate_mod.KINDS))
        pairs = []
        arms_ = []
        n_records = 0
        for p in ps:
            if p.get("scenario"):
                pair = (p["scenario"], p.get("arm") or "")
                if pair not in pairs:
                    pairs.append(pair)
            elif p.get("arm"):
                if p["arm"] not in arms_:
                    arms_.append(p["arm"])
            else:
                n_records += 1
        parts_ = []
        if pairs:
            parts_.append("%d scenario-arm pair%s: %s" % (len(pairs), "" if len(pairs) == 1 else "s", ", ".join(
                "%s %s" % (_link_sid(ctx, sid), esc(arm)) for sid, arm in pairs[:12]) + (
                " … (%d)" % len(pairs) if len(pairs) > 12 else "")))
        if arms_:
            n_arm = sum(1 for p in ps if p.get("arm") and not p.get("scenario"))
            parts_.append("%d problem%s of the access path%s %s" % (n_arm, "" if n_arm == 1 else "s",
                                                                   "" if len(arms_) == 1 else "s", esc(", ".join(arms_))))
        if n_records:
            parts_.append("%d problem%s of the run's records (the checks, the comparison, the table)" % (
                n_records, "" if n_records == 1 else "s"))
        scope = "; ".join(parts_) or "the run's records"
        rows.append([esc(title), "%d (%s)" % (len(ps), esc(kinds_text)), scope, esc(GATE_GROUP_ACTION[key])])
    return latex.plain_table(latex.spec_p([30, 26, 50, 68]),
                             ["group", "problems (by the gate's kind)", "scope", "what it implies and what clears it"],
                             rows, "the gate's problems by what they ask of the owner; the full texts follow, and every "
                                   "problem in full is in {\\upshape %s} of this document's build directory" % tt(
                                       GATE_PROBLEMS_FILE),
                             needspace=12)


def _deployed(run):
    """The identity input's deployed commit record (make_doc --identity: `deployed_service_commit` {commit, source,
    railway_deployment_id}), None without one."""
    ident = getattr(run, "identity", None)
    dep = ident.get("deployed_service_commit") if isinstance(ident, dict) else None
    return dep if isinstance(dep, dict) and dep.get("commit") else None


def _break_hex(text):
    return re.sub(r"[0-9a-f]{32,}", lambda m: "\\allowbreak{}".join(
        m.group(0)[i:i + 8] for i in range(0, len(m.group(0)), 8)), text)


def identity_note_text(ctx):
    """(appendix F's note on the service identity in plain words, the gate notes it replaces): the recorded commit
    equal to the deployed one on all invocations, or self-reported and not cross-checked when no deployed commit was
    given; ('', []) otherwise (a run recording no identity keeps the gate's note: v5, v6)."""
    run = ctx.run
    gate = run.gate if isinstance(run.gate, dict) else {}
    si = gate.get("service_identity") if isinstance(gate.get("service_identity"), dict) else {}
    rec = si.get("recorded") if isinstance(si.get("recorded"), dict) else None
    if not rec or not rec.get("recorded") or not rec.get("service_commits"):
        return "", []
    expected = si.get("expected")
    arms = (gate.get("policy") or {}).get("arms") or ["mcp"]
    problems, notes = gate_mod.service_identity_problems(rec, expected, arms)
    commits = rec["service_commits"]
    n = rec.get("invocations") or sum(commits.values())
    if problems:
        return "", []
    if expected:
        dep = _deployed(run) or {}
        src = str(dep.get("source") or "the identity input")
        rid = dep.get("railway_deployment_id")
        text = ("the service's commit, as the MCP runner read it before each invocation, equals the deployed commit %s "
                "(%s%s) on all %d MCP invocations" % (expected, src, (", Railway deployment %s" % rid) if rid else "", n))
    elif len(commits) == 1:
        text = ("the service identified itself as commit %s on every one of the %d MCP invocations, so one deployment "
                "served the whole run; the gate was given no record of the intended deployment, so this identity is "
                "self-reported, not cross-checked" % (list(commits)[0], n))
    else:
        return "", []
    return _break_hex(esc(text)), notes


def by_design_tie_text(run):
    """The by-design table tied to the service that answered (RV3): the commit and date it was generated at, and how
    many of the run's by-design pairs were verified against the deployed commit (the identity input's, equal to the
    recorded one) or the recorded one; '' for a run recording no service identity."""
    table = run.table_meta or {}
    rec = recorded_service_identity(run) or {}
    commits = list(rec.get("service_commits") or {})
    if not table.get("service_commit") or len(commits) != 1:
        return ""
    verified = unverified = 0
    for sid in run.order:
        ma = run.comp_entry(sid).get("mcp_agreement") or {}
        if ma.get("class") == "service_by_design":
            if ma.get("verified") is False:
                unverified += 1
            else:
                verified += 1
    dep = _deployed(run)
    against = ("the deployed %s" % str(dep["commit"])[:8]) if dep and gate_mod.commit_matches(commits[0], dep["commit"]) \
        else "the recorded service %s" % commits[0][:8]
    when = str(table.get("generated_utc") or "")[:10]
    pairs = ("all %d by-design pairs verified" % verified) if not unverified else (
        "%d of %d by-design pairs verified" % (verified, verified + unverified))
    return "generated at service commit %s%s; %s against %s" % (
        str(table["service_commit"])[:8], (" (%s)" % when) if when else "", pairs, against)


def gate_problems(ctx):
    """The deploy gate's problems (READ-002, READ-023): what the gate read, the problems by group with the evidence
    that clears each, then by kind (the scenarios linked, every configuration text once, none left out), and the
    full record written beside the document (GATE_PROBLEMS_FILE)."""
    run = ctx.run
    gate = run.gate
    word, how, counts_text, problems_text = gate_summary(ctx)
    if not isinstance(gate, dict):
        return "\\textbf{%s}: %s.\n\n" % (esc(word), esc(how))
    out = "\\textbf{%s}: %s; %s; %s.\n\n" % (esc(word), esc(how), esc(counts_text), esc(problems_text))
    out += "%s What the gate read: %s.\n\n" % (inline(gate_decision_text(run)), esc(gate_read_text(run)))
    out += gate_ledger(ctx)
    problems = [p for p in (gate.get("problems") or []) if isinstance(p, dict)]
    # no local path in the delivered record (the verification of the review's fixes, L5): the policy's path and any
    # path in a problem's text or a note in the suite's tree, as evidence.json and the document print them
    pol_rec = dict(gate.get("policy") or {})
    if pol_rec.get("path"):
        pol_rec["path"] = latex.public_paths(pol_rec["path"])
    ctx.data_files[GATE_PROBLEMS_FILE] = json.dumps(
        {"verdict": "GATE PASS" if gate.get("passed") else "GATE FAIL", "policy": pol_rec,
         "counts": gate.get("counts"), "judged_offline_for_this_document": bool(run.gate_offline),
         "problems": [dict(p, group=gate_group_of(p), detail=latex.public_paths(p.get("detail")))
                      if isinstance(p.get("detail"), str) else dict(p, group=gate_group_of(p)) for p in problems],
         "notes": [latex.public_paths(n) if isinstance(n, str) else n for n in (gate.get("notes") or [])]},
        indent=1, sort_keys=True) + "\n"
    if problems:
        out += gate_action_table(ctx, problems)
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
                text = latex.public_paths(str(p.get("detail") or "?"))   # no local path in F either (R3)
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
            for text, n, arms_ in details:
                parts.append("%s%s: %s" % (esc(", ".join(arms_)) + " " if arms_ else "",
                                           ("%d times" % n) if n > 1 else "once", esc_wrap(text, 320)))
            rows.append([esc(k), "; ".join(parts), str(len(v))])
        out += latex.plain_table(latex.spec_p([30, 134], extra_cols=1, extra_spec="r"),
                                 ["problem kind", "scenario and access path, or the problem's text", "n"], rows,
                                 "every problem by the gate's kind (a text cut here is complete in {\\upshape %s})" % tt(GATE_PROBLEMS_FILE))
    # the service identity in plain words (RV3, V3): the recorded commit against the deployed one the identity input
    # gives (the gate's own notes name the command-line flag and the runner's file; they stay in gate_problems.json)
    ident_text, ident_notes = identity_note_text(ctx)
    if ident_text:
        out += "{\\small\\color{black!60}note: %s}\\par\n" % ident_text
    for n in [x for x in (gate.get("notes") or []) if x not in ident_notes][:6]:
        out += gate_note_tex(n)
    return out + "\n"


def gate_note_tex(note):
    """One gate note as appendix F prints it: a long path and a full commit sha may break (the v7 baseline's identity
    note named the 40-hex commit and ran 22 pt past the margin), no local path (as the record), the text cut at 300
    characters with its end kept."""
    return "{\\small\\color{black!60}note: %s}\\par\n" % _break_hex(esc_wrap(latex.public_paths(str(note)), 300))


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
                     esc_wrap(sa.get("missing_argument") or sa.get("detail"), 160), cells.attn(str(len(sids)))])
    out += "\\parhead{Cases of the by-design table of agreement not implemented via MCP that scenarios name}\n"
    out += latex.plain_table(latex.spec_p([50, 60, 52], extra_cols=1, extra_spec="r"),
                             ["case", "scenario match", "missing argument", "scenarios"], rows)
    reasons = {}
    for sid in run.order:
        mcp = run.scenario(sid).get("mcp") or {}
        if mcp.get("not_implemented"):
            reasons.setdefault(str(mcp["not_implemented"]), []).append(sid)
    out += "\\parhead{The manifest's reasons for not implemented via MCP}\n"
    rows = [[esc_wrap(r, 120), cells.attn(str(len(s))), esc(", ".join(s[:8]) + (" …" if len(s) > 8 else ""))]
            for r, s in sorted(reasons.items(), key=lambda kv: (-len(kv[1]), kv[0]))]
    out += latex.plain_table(latex.spec_p([60, 104], extra_cols=1, extra_spec="r"),
                             ["reason", "ids", "scenarios"], [[r[0], r[2], r[1]] for r in rows])
    out += property_marks_table(ctx)
    out += property_by_design_table(ctx)
    dec = run.declared
    out += "\\parhead{Owner-declared exemptions}\n"
    out += ctx.prose_paragraphs("provenance", "exemption_kinds", optional=True) if ctx.overview_section(
        "provenance").get("exemption_kinds") else exemption_kinds_default()
    if isinstance(dec, dict):
        out += "File %s (declared version %s).\n" % (tt(dec.get("_file") or "declared"), esc(dec.get("declared_version")))
        rows = []
        for kind in ("exclusions", "not_required", "relations_excluded"):
            for e in dec.get(kind) or []:
                if isinstance(e, dict):
                    rows.append([esc(kind.replace("_", " ")), latex.breakable_id(e.get("scenario") or e.get("relation") or "–"),
                                 esc(", ".join(e.get("arms") or [])), esc(exemption_consequence(kind, e, run)),
                                 esc_wrap(e.get("reason"))])
        out += latex.plain_table(latex.spec_p([18, 26, 18, 44, 68]),
                                 ["kind", "scenario or relation", "access paths", "what it changes for this gate",
                                  "the owner's reason"], rows)
        out += ("\\noindent{\\small An exemption changes what the gate requires, not what the run observed: the test's "
                "recorded outcome stands as printed in its chapter; what its exempted executions showed stays shown, "
                "and what they did not show stays unproven there by this run (the fourth column, read from the "
                "record). A reason that names a later manifest version describes a plan for that version, not "
                "evidence from this run.}\\par\n")
    else:
        out += "None given for this manifest.\n\n"
    return out


def property_marks_table(ctx):
    """Appendix G's property-level marks not implemented via MCP (manifest v5: `mcp.not_implemented_checks`, one
    check of a scenario; VERIFY-COLOURS M1): per check and reason the scenarios that carry the mark and, in red,
    those whose MCP execution had the check set aside in this run (cells.property_not_implemented: the test counts
    under not implemented via MCP in the attention tags, its MCP cell carries the property-level mark). Empty for a
    manifest without such marks (v1 to v4)."""
    run = ctx.run
    marked = {}
    for sid in run.order:
        marks = run.check_marks(sid) if hasattr(run, "check_marks") else {}
        for check, m in sorted(marks.items()):
            if m.get("mark") == "not_implemented":
                marked.setdefault((check, str(m.get("reason") or "–")), []).append(sid)
    if not marked:
        return ""
    rows = []
    for (check, reason), sids in sorted(marked.items(), key=lambda kv: (-len(kv[1]), kv[0])):
        whole = [s for s in sids if ctx.state(s)["mcp"].glyph_key == "not_implemented_mcp"]
        aside = [s for s in sids if s not in whole
                 and cells.property_not_implemented(run, s, run.arm_lines(s, "mcp"))]
        rest = [s for s in sids if s not in aside and s not in whole]
        names = []
        if aside:
            names.append(cells.attn("set aside in this run: %s" % esc(facts._ids_short(aside))))
        if whole:
            names.append("%s: %s" % (cells.attn("the whole test not exercised via MCP (\\gni)"),
                                     esc(facts._ids_short(whole))))
        if rest:
            names.append("marked, the check read otherwise on every MCP execution (judged as any check): %s" % esc(
                facts._ids_short(rest)))
        rows.append([tt(check), esc_wrap(reason, 160), "; ".join(names), cells.attn(str(len(aside))) if aside else "0"])
    out = "\\parhead{The manifest's property-level marks not implemented via MCP}\n"
    out += ("\\noindent{\\small A mark of one check of a scenario: on the MCP access path the marked check, read not "
            "exposed, is set aside, the rest of the scenario stays tested; the test counts under not implemented via "
            "MCP in the attention tags, and its MCP cell carries the property-level mark (I.3).}\\par\\smallskip\n")
    out += latex.plain_table(latex.spec_p([34, 64, 66], extra_cols=1, extra_spec="r"),
                             ["check", "the manifest's reason", "scenarios", "set aside"], rows)
    return out


def property_by_design_table(ctx):
    """Appendix G's property-level marks service-side by design (`mcp.by_design_checks` [{check, case}]; manifest v7:
    P-edge-0004's entry_fields under pattern.lowercase): per check and case the scenarios that carry the mark. On the
    MCP access path the marked check, read not exposed, is set aside as a documented difference (appendix B gives the
    case), the rest of the test compared; it stays required on the direct access paths. Empty for a manifest without
    such marks (v5, v6)."""
    run = ctx.run
    marked = {}
    for sid in run.order:
        marks = run.check_marks(sid) if hasattr(run, "check_marks") else {}
        for check, m in sorted(marks.items()):
            if m.get("mark") == "by_design":
                marked.setdefault((check, str(m.get("case") or "–")), []).append(sid)
    if not marked:
        return ""
    rows = [[tt(check), tt(case), esc(facts._ids_short(sids)), str(len(sids))]
            for (check, case), sids in sorted(marked.items(), key=lambda kv: (-len(kv[1]), kv[0]))]
    out = "\\parhead{The manifest's property-level marks service-side by design}\n"
    out += ("\\noindent{\\small A mark of one check of a scenario as a documented service difference: on the MCP access "
            "path the marked check, read not exposed, is set aside under the by-design table's case (appendix B), the "
            "rest of the test stays tested and compared; the check stays required on the direct access paths, and the "
            "test does not count under not implemented via MCP.}\\par\\smallskip\n")
    out += latex.plain_table(latex.spec_p([34, 60, 70], extra_cols=1, extra_spec="r"),
                             ["check", "the by-design case", "scenarios", "n"], rows)
    return out


def exemption_kinds_default():
    return ("The declared file holds three kinds of entry, each a release-policy decision of the owner with its "
            "reason. An \\emph{exclusion} drops one scenario on the named access paths from the gate's required "
            "executions; a \\emph{check not required} drops one property check of one scenario on those paths, the "
            "scenario's other checks and the other paths staying required; a \\emph{relation excluded} drops one "
            "relation on those paths. None of them changes a recorded outcome or shows the omitted behaviour.\n\n")


def exempted_evidence(run, entry):
    """What the executions an exclusion exempts showed, from the record: (outcomes, shown, not_shown, partly) over
    the scenario's lines on the entry's access paths (and repeats, when it names some): the recorded outcome words
    per access path and repeat; the property checks every item of which passed on some exempted execution; those
    no item of which passed there; and {check: (passed items, items)} where some items passed and others did not
    (the readability review's verification of 2026-10-16, shared M2: an exemption shows nothing, so what passed
    stays passed and only what did not pass stays unproven)."""
    sid = entry.get("scenario")
    reps = entry.get("repeats")
    lines = []
    for arm in (entry.get("arms") or list(ARMS)):
        for l in run.arm_lines(sid, arm):
            if reps is None or l.get("repeat") in reps:
                lines.append((arm, l))
    outcomes = ["%s r%s %s" % (cells.ARM_SHORT.get(arm, arm), l.get("repeat"),
                               cells.WORD.get(l.get("outcome"), str(l.get("outcome")))) for arm, l in lines]
    tally = {}
    for _, l in lines:
        for c in l.get("checks") or []:
            if isinstance(c, dict) and run.is_property(c.get("check")):
                t = tally.setdefault(str(c.get("check")), [0, 0])
                t[1] += 1
                if c.get("outcome") == "pass":
                    t[0] += 1
    shown = sorted(n for n, (p, m) in tally.items() if p == m)
    not_shown = sorted(n for n, (p, m) in tally.items() if p == 0)
    partly = dict((n, (p, m)) for n, (p, m) in tally.items() if 0 < p < m)
    return outcomes, shown, not_shown, partly


def exemption_consequence(kind, entry, run=None):
    """The coverage consequence of one declared entry in plain words (READ-025): for an exclusion, read from the
    record when the run is given (exempted_evidence), so that a property the exempted executions passed is not
    called unproven."""
    arms_ = ", ".join(entry.get("arms") or []) or "every access path"
    if kind == "exclusions":
        reps = entry.get("repeats")
        where = (" (repeat%s %s; the other repeats stay required)" % (
            "" if len(reps) == 1 else "s", ", ".join(str(r) for r in reps))) if reps else ""
        if run is None or not entry.get("scenario"):
            return ("not required on %s%s; the exemption shows nothing: the recorded outcomes there stand, what the "
                    "exempted executions showed stays shown, and what they did not show stays unproven there by this "
                    "run" % (arms_, where))
        outcomes, shown, not_shown, partly = exempted_evidence(run, entry)
        head = "not required on %s%s; the exemption shows nothing: the recorded outcome%s there stand%s" % (
            arms_, where, "" if len(outcomes) == 1 else "s", "s" if len(outcomes) == 1 else "")
        if not outcomes:
            return head + " (no execution recorded there)"
        text = head + " (%s)" % "; ".join(outcomes)
        if not shown and not not_shown and not partly:
            scen = run.scenario(entry["scenario"]) if hasattr(run, "scenario") else None
            if not scen or blocks.host_expects_refusal(scen):
                return text + "; no property check there: the status, code and text are the claim"
            return text + "; no property check there: the test declares none of its own"
        parts = []
        if shown:
            parts.append("%s passed there and stay%s shown" % (", ".join(shown), "s" if len(shown) == 1 else ""))
        for name, (p, m) in sorted(partly.items()):
            parts.append("%s passed %d of its %d items there, the other%s not shown there by this run" % (
                name, p, m, "" if m - p == 1 else "s"))
        if not_shown:
            parts.append("%s %s not shown there by this run" % (", ".join(not_shown), "was" if len(not_shown) == 1 else "were"))
        return text + "; " + "; ".join(parts)
    if kind == "not_required":
        return ("the check %s is not required on %s; the scenario's other checks stay required there, and this "
                "check stays required on the other access paths" % (entry.get("check") or "?", arms_))
    if kind == "relations_excluded":
        return "the relation is not required on %s; the compared answers stand as recorded" % arms_
    return "–"


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
    out += reruns_table(ctx)
    out += appendix_rules(ctx)
    out += "\\parhead{Vocabulary}\n" + classes_tables(bool(answer_side_mark_ids(ctx.run)))
    out += colour_rules(ctx)
    out += reading_reference(ctx)
    out += "\\parhead{The deploy gate}\n" + gate_problems(ctx)
    # written after the gate's record, whose sha256 it carries (the second review, R2-002)
    ctx.data_files[EVIDENCE_FILE] = json.dumps(evidence_manifest(ctx), indent=1, sort_keys=True) + "\n"
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
