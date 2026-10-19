#!/usr/bin/env python3
"""Report of a suite run: RUN_DIR/report.md from run.json, executions, samples and compare.json, then
report.pdf via pandoc + xelatex (DESIGN §6).

    report.py RUN_DIR [--compare EARLIER_RUN_DIR] [--no-pdf] [--max-failures 200] [--reevaluated] [--reuse]

--reevaluated: the run's re-evaluated set (reevaluate.py) into report.reevaluated.md / .pdf, with a section on
what the re-evaluation judged with and every execution whose outcome or checks it changed.

Sections: verdict and environment (with the load caveat), outcomes per area and arm, every failure with its
expectation and the observed value, agreement of the direct arms and of the MCP arm, chapter G (manifest v8's
graphlet tools compared five ways: per tool, disagreements, pre-checks, service checks, correctness coverage; only
for a manifest with graphlet scenarios), timing per area and arm
(p50 / p90 / max; cold vs warm), overheads (front; MCP phases), resource use per area, and with --compare the
deltas against the earlier run. Tables are plain pipe tables with ASCII-safe cells. compare.json is computed
when missing (compare.py).
"""

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare as cmp  # noqa: E402
import gate as gate_mod  # noqa: E402
import suite_common as sc  # noqa: E402

ARM_LABEL = {"direct-local": "direct-local", "direct-remote": "direct-remote", "mcp": "mcp"}
# The load above which the timings carry the Logan caveat: sra_all busier than 4 cores, or load1 above half the
# 96 cores of mex.
LOAD_CAVEAT_SRA_CPU_PCT = 400.0
LOAD_CAVEAT_LOAD1 = 48.0
# The execution times the report marks (the owner, 2026-10-15): above TIME_MARK_YELLOW_MS yellow, above
# TIME_MARK_RED_MS red (fmt_time: client wall, server elapsed and MCP phase times, in the tables and the lines;
# differences between arms and connection set-ups are not execution times and stay plain). The mark is a span with
# a class and the tint inline: report_colour.lua turns the class into a \colorbox on the md -> pdf path (pandoc +
# xelatex), a Markdown viewer that renders HTML shows the inline tint, a terminal the plain number. The text stays
# black on both tints, so it stays legible.
TIME_MARK_YELLOW_MS = 10000.0
TIME_MARK_RED_MS = 55000.0
TIME_MARK_CLASS = {"yellow": "t-yellow", "red": "t-red"}
TIME_MARK_TINT = {"yellow": "#ffe866", "red": "#ff8a80"}
LUA_FILTER = os.path.join(os.path.dirname(os.path.abspath(__file__)), "report_colour.lua")

_ASCII = {
    "≤": "<=", "≥": ">=", "⊆": "subset of", "⊇": "superset of", "→": "->",
    "←": "<-", "×": "x", "–": "-", "—": "--", "‘": "'", "’": "'",
    "“": '"', "”": '"', "…": "...", "µ": "u", "μ": "u", "§": "S", "±": "+/-",
    "≠": "!=", "≈": "~", " ": " ", "−": "-", "σ": "sigma", "∑": "sum",
}
_SEQ_RE = re.compile(r"[ACGTNacgtnRYSWKMBDHVU]{41,}")
_MD_SPECIAL = re.compile(r"([\\`*_\[\]<>$~^|#{}])")


def ascii_safe(text):
    """Text with the common non-ASCII signs spelled out, anything else as '?', and long sequences cut."""
    if text is None:
        return ""
    t = str(text)
    t = "".join(_ASCII.get(ch, ch) if ord(ch) > 127 else ch for ch in t)
    t = t.encode("ascii", "replace").decode("ascii")
    t = _SEQ_RE.sub(lambda m: "%s...(%d bp)" % (m.group(0)[:16], len(m.group(0))), t)
    return t


def md(text, limit=None):
    """ASCII-safe, markdown-escaped text on one line, cut at `limit` characters."""
    t = ascii_safe(text).replace("\r", " ").replace("\n", " ")
    if limit and len(t) > limit:
        t = t[:limit - 3] + "..."
    return _MD_SPECIAL.sub(r"\\\1", t)


_NUMERIC_CELL = re.compile(r"^(-|[-+]?[0-9][0-9.,]*%?|[0-9]+/[0-9]+)$")
_SPAN_TAG = re.compile(r"</?span\b[^>]*>")


def visible(cell):
    """A cell's text as a reader sees it: the span tags of a colour mark (fmt_time) taken out."""
    return _SPAN_TAG.sub("", str(cell))


def table(header, rows, align=None):
    """A pipe table; cells are already escaped strings.

    Numeric columns are right-aligned, others left (unless `align` says). The separator's dashes follow each
    column's widest cell: pandoc sizes the columns of a wide table by them, so a short column (an area) stays
    narrow and a long one (a detail) gets the room, instead of equal widths that wrap every cell. Both read the
    cell's visible text: a colour mark's span tags take no room on the page.
    """
    if not rows:
        return "(none)\n"
    cols = len(header)
    if align is None:
        align = ["r" if all(_NUMERIC_CELL.match(visible(r[i]).strip()) for r in rows) else "l" for i in range(cols)]
    widths = [max([len(str(header[i]))] + [len(visible(r[i])) for r in rows]) for i in range(cols)]
    sep = [("-" * max(3, min(w, 60)) + ":") if a == "r" else (":" + "-" * max(3, min(w, 60)))
           for a, w in zip(align, widths)]
    out = ["| " + " | ".join(header) + " |", "|" + "|".join(sep) + "|"]
    for r in rows:
        out.append("| " + " | ".join(str(c) for c in r) + " |")
    return "\n".join(out) + "\n"


def fmt_ms(v):
    if v is None:
        return "-"
    return "%.1f" % v if abs(v) < 100 else "%.0f" % v


def time_mark(v):
    """The mark of an execution time in ms: "red" above TIME_MARK_RED_MS, "yellow" above TIME_MARK_YELLOW_MS, else
    None (a missing time too)."""
    if v is None:
        return None
    if v > TIME_MARK_RED_MS:
        return "red"
    if v > TIME_MARK_YELLOW_MS:
        return "yellow"
    return None


def marked(text, mark):
    """|text| in the span of |mark| ("yellow" / "red"; None: as it is): the class report_colour.lua colours in the
    PDF, the tint inline for a Markdown viewer that renders HTML."""
    if mark is None:
        return text
    return '<span class="%s" style="background-color:%s">%s</span>' % (TIME_MARK_CLASS[mark], TIME_MARK_TINT[mark],
                                                                        text)


def fmt_time(v):
    """An execution time in ms (client wall, server elapsed, an MCP phase) as fmt_ms shows it, marked yellow above
    TIME_MARK_YELLOW_MS and red above TIME_MARK_RED_MS (time_mark, marked) by the value as printed (fmt_ms rounds
    100 ms and above to whole milliseconds: 10,000.4 ms prints 10000 and is not above 10 s); never for a difference
    or a connection set-up, which fmt_ms prints plain."""
    text = fmt_ms(v)
    if v is None:
        return text
    try:
        shown = float(text)
    except ValueError:
        shown = v
    return marked(text, time_mark(shown))


def time_legend():
    """The one legend sentence of the colour marks."""
    return ("Execution times (client wall, server elapsed, MCP phases; ms) above %d s are marked %s and above %d s %s, "
            "in the tables and the lines, a time quoted in a check's text or a gate problem too, by the value as "
            "printed; differences between arms and connection set-ups are not marked, nor a budget, a cap or a "
            "lateness. Red means needs attention: every count and list of what is %s, and the verdict's attention "
            "line, are red too." % (
                TIME_MARK_YELLOW_MS / 1000.0, marked("yellow", "yellow"), TIME_MARK_RED_MS / 1000.0,
                marked("red", "red"), NIM_WORDS))


# The execution times a check's free-text detail quotes (a finding's or failure's check line, an observed error, a
# gate problem's detail), marked as fmt_time marks a printed one, by the value as printed (VERIFY-V6-2 M2): the figure
# right after a client wall ("client wall N ms", "client wall time of N ms", "wall time of N ms", "wall N s"), a
# server elapsed ("elapsed N ms"), a read ("a read of N ms"), "took N s" and a stop, cut or deadline "after N ms";
# a lateness ("N ms after the deadline": its figure comes before "after"), a budget, a cap, a reserve and a bound
# follow none of these words and stay plain.
_QUOTED_TIME = re.compile(r"\b(client wall time of|client wall|wall time of|wall|elapsed|took|a read of|after)(\s+)"
                          r"([0-9][0-9,]*(?:\.[0-9]+)?)(\s?)(ms|s)\b")


def mark_quoted_times(text):
    """|text| (Markdown: md's output) with every execution time it quotes (_QUOTED_TIME) marked by time_mark of the
    figure as printed (its unit read: s or ms); the unit stays outside the mark, as fmt_time's lines have it."""
    def sub(m):
        try:
            v = float(m.group(3).replace(",", ""))
        except ValueError:
            return m.group(0)
        mark = time_mark(v * 1000.0 if m.group(5) == "s" else v)
        if mark is None:
            return m.group(0)
        return "%s%s%s%s%s" % (m.group(1), m.group(2), marked(m.group(3), mark), m.group(4), m.group(5))
    return _QUOTED_TIME.sub(sub, text or "")


# What needs attention is red, as a time above TIME_MARK_RED_MS (the owner, 2026-10-16): what is not implemented yet
# via MCP (NIM_WORDS; a missing MCP feature: a whole execution, outcome not_implemented_mcp, or a property-level mark,
# mcp.not_implemented_checks) wherever the report counts or lists it with something in it (the verdict's counts and
# agreement line, the area table, the agreement section, the property marks), and the verdict's attention line
# (attention_line). A count of 0 stays plain: nothing to attend to.
NIM_WORDS = "not implemented yet via MCP"


def nim(text):
    """|text| (already Markdown) marked red: a count or a list of what is not implemented yet via MCP."""
    return marked(text, "red")


def nim_count(n, words=NIM_WORDS):
    """'N not implemented yet via MCP' (|words|), red when N is not 0."""
    text = "%d %s" % (n or 0, words)
    return nim(text) if n else text


_NIM_COUNT = re.compile(r"(?<![0-9.,])([0-9][0-9,]*) %s" % re.escape(NIM_WORDS))


def mark_nim_counts(text):
    """|text| (Markdown: a verdict composed before its marks) with every 'N not implemented yet via MCP' of a nonzero
    N marked red."""
    return _NIM_COUNT.sub(lambda m: m.group(0) if int(m.group(1).replace(",", "")) == 0 else nim(m.group(0)),
                          text or "")


def fmt_num(v, digits=1):
    if v is None:
        return "-"
    if isinstance(v, int) and not isinstance(v, bool):
        return str(v)
    return ("%%.%df" % digits) % v


def fmt_gib(v):
    return "-" if v is None else "%.2f" % (v / 1024.0 ** 3)


def fmt_mib(v):
    return "-" if v is None else "%.1f" % (v / 1024.0 ** 2)


def short(text, n=12):
    return "-" if not text else str(text)[:n]


def short_commit(commit):
    """A commit cut to 12 characters, keeping a `+dirty` marker (review L1)."""
    if not commit:
        return "-"
    text = str(commit)
    dirty = text.endswith("+dirty")
    return text.replace("+dirty", "")[:12] + ("+dirty" if dirty else "")


HOST_CACHE_ORDER = {"cold": 0, "warm": 1, cmp.WARM_AFTER_INTERRUPTED: 1.5, "aged": 2, "not_reached": 3, "unknown": 4}
HOST_CACHE_LABEL = {"cold": "cold", "warm": "warm", cmp.WARM_AFTER_INTERRUPTED: "warm after a cut", "aged": "aged",
                    "not_reached": "not reached",
                    "unknown": "unknown"}


def checks_name(info):
    """The checks module an arm ran: run_direct records "checks" / "stub", the MCP runner {file, stub}."""
    mod = (info or {}).get("checks_module")
    if isinstance(mod, dict):
        return "stub" if mod.get("stub") else (os.path.basename(str(mod.get("file") or "")) or "checks")
    return mod


def orchestrator_session(run):
    """The orchestrator's last invocation (run_suite.sh), or {}."""
    sessions = ((run.get("orchestrator") or {}).get("sessions")) or []
    return sessions[-1] if sessions else {}


def loop_text(loop):
    """The orchestrator's loop record (run_suite.sh loop_record, the last invocation's) in words: its state (reason),
    the paired loop's completed of planned scenarios (the sequential loop records no per-scenario count), the direct
    requests, the failed steps and the MCP executions the runner did not record."""
    out = ["the last invocation's loop %s%s" % (loop.get("state") or "-",
                                                (" (%s)" % loop["reason"]) if loop.get("reason") else "")]
    if "scenarios" in loop or "completed" in loop:
        out.append("%s of %s scenario(s) completed" % (loop.get("completed", "-"), loop.get("scenarios", "-")))
    if loop.get("direct_requests") is not None:
        out.append("%s direct request(s)" % loop["direct_requests"])
    if "failed_steps" in loop:
        out.append("%d failed step(s)" % len(loop.get("failed_steps") or []))
    if "mcp_missing" in loop or "mcp_missing_count" in loop:
        out.append("%s MCP execution(s) the runner did not record" % loop.get(
            "mcp_missing_count", len(loop.get("mcp_missing") or [])))
    return ", ".join(out)


def loop_problems_text(loop):
    """The loop's failed steps (ARM:SCENARIO:EXIT) and the MCP executions the runner did not record (SCENARIO:rREPEAT,
    the first 50 kept), "" when there are none."""
    steps = loop.get("failed_steps") or []
    missing = loop.get("mcp_missing") or []
    n_missing = loop.get("mcp_missing_count", len(missing))
    out = []
    if steps:
        out.append("%d failed step(s): %s" % (len(steps), ", ".join(str(x) for x in steps[:20]) +
                                               (" ..." if len(steps) > 20 else "")))
    if n_missing:
        out.append("%s MCP execution(s) the runner did not record: %s" % (
            n_missing, ", ".join(str(x) for x in missing[:20]) + (" ..." if n_missing > 20 else "")))
    return "; ".join(out)


def mixed_releases_text(run):
    """The host's release or index_fp changes during the run (merge_run.release_identity: run.json release_identity
    changed, each naming the arm, session and time that saw the new value) and the releases seen, "" when the run
    saw one deployment."""
    ident = run.get("release_identity") or {}
    changed = ident.get("changed") or []
    if not run.get("mixed_releases") and not changed:
        return ""
    return "%s (releases seen: %s)" % ("; ".join(changed) or "mixed releases", ", ".join(
        str(x) for x in ident.get("releases") or []) or "-")


def load_compare(run_dir, by_design_path=None, reevaluated=False):
    """The run's comparison (compare.load_or_compute: computed as compare.py computes it when absent; TableError
    when compare.py refuses the table)."""
    return cmp.load_or_compute(run_dir, by_design_path, reevaluated)


def area_sort_key(area):
    order = "CRTPLXMQNDS"
    return (order.index(area) if isinstance(area, str) and len(area) == 1 and area in order else len(order), str(area))


# ---------------------------------------------------------------------------------------------------------------
# sections

# the share of failed samples (the sampler's `error`: docker stats failed, e.g. a missing container or daemon) above
# which the load figures do not stand for the run
LOAD_SAMPLE_ERROR_SHARE = 0.10


def load_figures(samples):
    """The sampler's figures over the run: sra_all CPU and host load1 (mean, max, how many samples state each: a
    finite figure, suite_common.finite_or_none), the samples, and those that failed (`error`: their figures are
    null) with the first error."""
    sra = [s.get("sra_all_cpu_pct") for s in samples if sc.finite_or_none(s.get("sra_all_cpu_pct")) is not None]
    load = [s.get("host_load1") for s in samples if sc.finite_or_none(s.get("host_load1")) is not None]
    errors = [s.get("error") for s in samples if s.get("error")]
    return {"sra_mean": (sum(sra) / len(sra)) if sra else None, "sra_max": max(sra) if sra else None,
            "load_mean": (sum(load) / len(load)) if load else None, "load_max": max(load) if load else None,
            "n": len(samples), "n_sra": len(sra), "n_load": len(load), "n_error": len(errors),
            "first_error": str(errors[0]) if errors else None}


def stated_mcp_limit(line, scenario, cases, comp=None):
    """The id of the by-design case that states why an MCP execution is inconclusive, else None: the case the scenario
    names (compare.named_cases) that matched the execution (its service_by_design_case, else compare's
    by_design_notes for the repeat) and whose explains_checks state every check item of the line that did not pass
    (suite_common.states_limit), each of which reads not_exposed (a part of the answer the MCP surface does not expose
    by design, e.g. capabilities_doc of capabilities.traverse_capabilities_document; an entry {check, detail_matches}
    only an item whose detail it matches, as the verdict pipeline reads it)."""
    if line.get("arm") != "mcp" or line.get("outcome") != "inconclusive" or not cases:
        return None
    cid = line.get("service_by_design_case")
    if not cid and comp:
        ma = (((comp.get("scenarios") or {}).get(line.get("scenario")) or {}).get("mcp_agreement") or {})
        cid = (ma.get("by_design_notes") or {}).get(str(line.get("repeat")))
    named = dict((c.get("id"), c) for c in cmp.named_cases(cases, scenario))
    case = named.get(cid) if cid else None
    if not case:
        return None
    entries = case.get("explains_checks") or []
    rest = [c for c in line.get("checks") or [] if c.get("outcome") != "pass"]
    if rest and all(any(sc.states_limit(e, c) for e in entries) for c in rest):
        return cid
    return None


TABLE_ORIGIN = {"snapshot": "; the run's snapshot (scripts/service_by_design.json)",
                "live": "; this tree's", "other": "; a table given to compare"}


def pipeline_note(judged):
    """A note, not a failure: the verdict pipeline (merge_run, compare, gate, report) runs with this tree's
    suite_common.py, and it differs from the one the checks judged with (compare.judging_checks)."""
    if not judged.get("live_suite_common_differs"):
        return ""
    return "; note: suite_common.py of this tree differs from the run's snapshot (sha256 %s, the snapshot's %s): " \
        "the checks judged with the snapshot's, the verdict pipeline (merge_run, compare, gate, report) with this " \
        "tree's" % (md(short(judged.get("pipeline_suite_common_sha256"), 12)),
                    md(short(judged.get("suite_common_sha256"), 12)))


def judged_with_text(judged):
    """Which checks.py judged the comparison's relations and content rule and the property reading
    (compare.judging_checks, compare.json judged_with): the run's snapshot, a re-evaluation's checks or this tree's."""
    if not isinstance(judged, dict):
        return "not stated (a comparison of before the field)"
    sha = md(short(judged.get("checks_sha256"), 12))
    if judged.get("source") == "snapshot":
        changed = judged.get("snapshot_changed") or []
        return "checks.py of the run's snapshot (scripts/, sha256 %s, with the snapshot's suite_common.py)%s%s" % (
            sha, ("; its %s changed after it was taken (SNAPSHOT.json)" % md(", ".join(changed))) if changed else "",
            pipeline_note(judged))
    if judged.get("source") == "reevaluation":
        return "checks.py of the re-evaluation (%s, sha256 %s)%s" % (md(judged.get("dir") or "-"), sha,
                                                                     pipeline_note(judged))
    return "checks.py of this tree (%s, sha256 %s): the set has no snapshot (or re-evaluation checks) of its own, so " \
        "a change of the tree since the run reaches this comparison" % (md(judged.get("dir") or "-"), sha)


def report_cases(comp, run_dir=None):
    """The by-design table compare judged the run with (compare.json by_design_table, a path inside the run relative
    to |run_dir|: compare.run_relative), loaded; {} without one."""
    info = (comp or {}).get("by_design_table")
    path = (info or {}).get("path") if isinstance(info, dict) else None
    if path and run_dir:
        path = cmp.run_path(path, run_dir)
    if not path or not os.path.exists(path):
        return {}
    try:
        return cmp.load_by_design(path)[0]
    except cmp.TableError:
        return {}


def untested_groups(rows, manifest, cases=None, comp=None):
    """[(kind "area" | "route", name, arm, {outcome: executions})] of every area and route x arm whose executions
    include an inconclusive one and none passed: what that part of the run tested is nothing (a field the MCP no
    longer exposes leaves every execution inconclusive). Executions not implemented via MCP or by design alone are
    explicit marks, counted elsewhere, and an MCP execution inconclusive by a limit the table states for the scenario
    (stated_mcp_limit, with |cases|) is stated apart; neither makes a group untested."""
    scen = dict((s.get("id"), s) for s in (manifest or {}).get("scenarios", []) if isinstance(s, dict))
    route_of = dict((sid, s.get("route")) for sid, s in scen.items())
    groups = {}
    for r in rows:
        arm = r.get("arm")
        if stated_mcp_limit(r, scen.get(r.get("scenario")), cases, comp):
            continue   # a limit the table states for the scenario: stated in the verdict section, not untested
        for kind, name in (("area", r.get("area") or str(r.get("scenario") or "").split("-", 1)[0]),
                           ("route", route_of.get(r.get("scenario")) or r.get("route"))):
            if name:
                c = groups.setdefault((kind, name, arm), {})
                o = sc.outcome_value(r.get("outcome"))
                c[o] = c.get(o, 0) + 1
    out = []
    for (kind, name, arm), c in groups.items():
        if c.get("inconclusive") and not c.get("pass"):
            out.append((kind, name, arm, c))
    arm_pos = dict((a, i) for i, a in enumerate(sc.ARMS))
    out.sort(key=lambda g: (g[0], area_sort_key(g[1]) if g[0] == "area" else (0, str(g[1])), arm_pos.get(g[2], 9)))
    return out


def untested_text(groups, limit=12):
    return "; ".join("%s %s on %s (%s)" % (kind, name, arm, ", ".join(
        "%d %s" % (n, o) for o, n in sorted(c.items(), key=lambda kv: (-kv[1], str(kv[0]))))) for kind, name, arm, c
        in groups[:limit]) + ("; ..." if len(groups) > limit else "")


# how a verdict names the outcomes of executions that passed nothing (no_pass_text)
NO_PASS_LABELS = (("not_implemented_mcp", NIM_WORDS), ("service_by_design", "service-side by design"),
                  ("inconclusive", "inconclusive"))


# how the verdict names the property-level marks (suite_common.MCP_CHECK_MARKS)
MARK_WORDS = {"not_implemented": "not implemented yet via MCP", "by_design": "service-side by design"}


def no_pass_text(outcomes):
    """'2 not implemented yet via MCP, 1 service-side by design, 1 inconclusive' of an {outcome: executions} that
    holds no pass (every other outcome by its name)."""
    labels = dict(NO_PASS_LABELS)
    order = [k for k, _ in NO_PASS_LABELS] + sorted(k for k in outcomes if k not in labels)
    return ", ".join("%d %s" % (outcomes[k], labels.get(k, k)) for k in order if outcomes.get(k))


def time_budget_issue(line):
    """Why an execution is a time-budget issue for the attention line, or '': it hit its budget (its own time
    limitation, compare._own_time_limitation: the recorded one, a stated stop by time, a time-limited answer, the
    host's 503 deadline or the service's relay of it, the front ingress's cut; else a 503 deadline, or a 504 without a
    host code, the ingress's cut, of a line recorded without it), or it took more than TIME_MARK_RED_MS by its wall
    as printed (fmt_ms: a red mark)."""
    why = cmp._own_time_limitation(line)
    if why:
        return why
    err = line.get("error") if isinstance(line.get("error"), dict) else {}
    if line.get("http") == 504 and not err.get("code"):
        return "a 504 without a host answer (the front ingress's cut)"
    wall = line.get("wall_ms")
    if isinstance(wall, (int, float)) and not isinstance(wall, bool):
        try:
            shown = float(fmt_ms(wall))
        except ValueError:
            shown = wall
        if time_mark(shown) == "red":
            return "took %s ms" % fmt_ms(wall)
    return ""


def attention_counts(rows, manifest, comp=None):
    """{kind: the TESTS (scenario ids) of that kind} over the run's executions, for the attention line: failed (a
    fail, or an outcome outside the vocabulary), errored, with a finding, with a disagreement between the arms (the
    comparison of |comp|: the direct pair's or the MCP's class disagree, e.g. an exact count that differs while every
    execution passes), the relations that failed (their ids, on any arm), with a time-budget issue
    (time_budget_issue on any execution), not implemented yet via MCP as a whole (an MCP execution
    not_implemented_mcp, or the comparison's MCP class of |comp|, as the agreement line counts them) or by a
    property-level mark (a scenario the MCP arm ran whose manifest marks a check not implemented yet via MCP). A test
    may be of several kinds."""
    scen = dict((x.get("id"), x) for x in (manifest or {}).get("scenarios", []) if isinstance(x, dict))
    out = dict((k, set()) for k in ("failed", "errored", "finding", "disagree", "relation_failed", "time_budget",
                                    "nim_whole", "nim_marks"))
    for r in rows or []:
        sid, o = r.get("scenario"), r.get("outcome")
        if o == "fail" or not sc.outcome_known(sc.recorded_outcome(r)):
            out["failed"].add(sid)
        elif o == "error":
            out["errored"].add(sid)
        if o == "finding":
            out["finding"].add(sid)
        if time_budget_issue(r):
            out["time_budget"].add(sid)
        if r.get("arm") == "mcp":
            if o == "not_implemented_mcp":
                out["nim_whole"].add(sid)
            if any(m.get("mark") == "not_implemented" for m in sc.mcp_check_marks(scen.get(sid)).values()):
                out["nim_marks"].add(sid)
    for sid, entry in ((comp or {}).get("scenarios") or {}).items():
        if ((entry or {}).get("mcp_agreement") or {}).get("class") == "not_implemented_mcp":
            out["nim_whole"].add(sid)
        if any(((entry or {}).get(k) or {}).get("class") == "disagree" for k in ("direct_agreement", "mcp_agreement")):
            out["disagree"].add(sid)
    for rel in (comp or {}).get("relations") or []:
        if any((v or {}).get("outcome") == "fail" for v in ((rel or {}).get("arms") or {}).values()):
            out["relation_failed"].add(rel.get("id"))
    return out


def attention_line(rows, manifest, comp=None):
    """The verdict's attention line, red (the owner, 2026-10-16: red means needs attention): the TESTS per kind, only
    the kinds with any, failures, errors and findings first (with the disagreements between the arms and the failed
    relations, which fail a run whose executions all pass), then the time-budget issues, then what is not
    implemented yet via MCP (whole executions and property-level marks); '' when nothing needs attention."""
    c = attention_counts(rows, manifest, comp)
    nim_ = c["nim_whole"] | c["nim_marks"]
    # the time-budget kind counts every test the owner's definition names; those the manifest tags budget_stress (a
    # warm-up's cold run, a heavy scenario: outside the gate, slow by design) are stated beside it, never left out
    scen = dict((x.get("id"), x) for x in (manifest or {}).get("scenarios", []) if isinstance(x, dict))
    stressed = sum(1 for sid in c["time_budget"] if "budget_stress" in ((scen.get(sid) or {}).get("tags") or []))
    parts = [what % len(c[key]) for key, what in (
        ("failed", "%d failed"), ("errored", "%d errored"), ("finding", "%d with a finding"),
        ("disagree", "%d with a disagreement between the arms"), ("relation_failed", "%d relation(s) failed"),
        ("time_budget", "%%d with a time-budget issue (an execution hit its budget: a 503 deadline, a time-limited "
                        "answer, a stop by time, the ingress's cut; or took more than %d s%s)" % (
                            TIME_MARK_RED_MS / 1000.0, ("; %d of them tagged budget_stress, slow by design" % stressed)
                            if stressed else ""))) if c[key]]
    if nim_:
        parts.append("%d %s (%d as a whole, %d by a property-level mark)" % (
            len(nim_), NIM_WORDS, len(c["nim_whole"]), len(nim_ - c["nim_whole"])))
    if not parts:
        return ""
    return marked("Needs attention, in tests (scenarios): %s." % md("; ".join(parts)), "red")


def verdict_text(run, comp, rows, reevaluated=False, missing=None, manifest=None, cases=None):
    """The verdict word and its reason, as verdict_section prints it, and the facts it was decided on: {incomplete,
    loop, mixed, stub, untested, totals, no_pass_arms}. The results document (doc/make_doc.py) states the same
    verdict from the same inputs through this function."""
    s = comp["summary"]
    totals = s.get("totals", {})
    failing = sum(t["outcomes"].get("fail", 0) + t["outcomes"].get("error", 0) for t in totals.values())
    # an outcome outside the vocabulary (a runner's mark not agreed first) is read as a failure, never as a pass
    other = sum(t["outcomes"].get(cmp.OTHER_OUTCOME, 0) for t in totals.values())
    findings = sum(t["outcomes"].get("finding", 0) for t in totals.values())
    inconclusive = sum(t["outcomes"].get("inconclusive", 0) for t in totals.values())
    passed = sum(t["outcomes"].get("pass", 0) for t in totals.values())
    executions = sum(t.get("executions", 0) for t in totals.values())
    states = dict((arm, (info or {}).get("state")) for arm, info in (run.get("arms") or {}).items())
    incomplete = ["%s (%s)" % (arm, st) for arm, st in sorted(states.items())
                  if st not in (None, "done") and arm != "mcp"]
    keys, source = missing or ([], None)
    if keys:
        # the run's plan reconciled with its records: an execution planned and never recorded is never a pass
        per_arm = {}
        for _, arm, _ in keys:
            per_arm[arm] = per_arm.get(arm, 0) + 1
        incomplete.append("%d planned execution(s) not recorded (%s; plan of %s)" % (
            len(keys), ", ".join("%s %d" % kv for kv in sorted(per_arm.items())), source))
    loop = orchestrator_session(run).get("loop") or {}
    if loop and loop.get("state") not in (None, "done"):
        incomplete.append("the %s loop (%s%s)" % (loop.get("mode", "paired"), loop.get("state"),
                                                  (": %s" % loop["reason"]) if loop.get("reason") else ""))
    # a run that spans two deployments (merge_run.release_identity) tested no single one: never PASS
    mixed = mixed_releases_text(run)
    dis_direct = s["direct_agreement"].get("disagree", 0)
    dis_mcp = s["mcp_agreement"].get("disagree", 0)
    rel_failed = sum((v or {}).get("fail", 0) for v in ((s.get("relations") or {}).get("per_arm") or {}).values())
    modules = sorted(set(checks_name(info) for info in (run.get("arms") or {}).values() if checks_name(info)))
    # a re-evaluated set is judged by checks.py whatever the runners ran
    stub = "stub" in modules and not reevaluated
    untested = untested_groups(rows, manifest, cases, comp)
    arm_pos = dict((a, i) for i, a in enumerate(sc.ARMS))
    no_pass_arms = [(arm, t["outcomes"]) for arm, t in sorted(totals.items(), key=lambda kv: (
        arm_pos.get(kv[0], 9), kv[0])) if t.get("executions") and not t["outcomes"].get("pass")]
    if incomplete:
        verdict = "INCOMPLETE (not finished: %s)" % ", ".join(incomplete)
    elif failing or other or dis_direct or dis_mcp or rel_failed:
        verdict = "FAIL (%d failed or errored execution(s)%s, %d direct and %d MCP disagreement(s)%s)" % (
            failing, (", %d with an outcome outside the vocabulary" % other) if other else "", dis_direct, dis_mcp,
            (", %d failed relation(s)" % rel_failed) if rel_failed else "")
    elif findings:
        verdict = ("FINDINGS (%d execution(s) where the host broke a stated bound without failing the request; no "
                   "failure, no disagreement)" % findings)
    elif inconclusive and not passed:
        # Nothing was shown to hold: never PASS (review N5).
        verdict = ("INCONCLUSIVE (no execution passed: %d inconclusive execution(s), where no check could test its "
                   "property; no failure, no disagreement)" % inconclusive)
    elif untested:
        # a part of the run tested nothing (e.g. a field the MCP stopped exposing): never PASS
        verdict = ("INCONCLUSIVE (no execution passed in %d area / route x arm group(s): %s; no failure, no "
                   "disagreement)" % (len(untested), untested_text(untested)))
    elif mixed:
        verdict = ("INCONCLUSIVE (the run spans more than one deployment of the host: %s; no failure, no "
                   "disagreement)" % mixed)
    elif not executions:
        verdict = "EMPTY (no execution recorded)"
    elif not passed:
        # only marks (not implemented via MCP, service-side by design): nothing was shown to hold, never PASS
        verdict = ("UNVERIFIED (no execution passed: %s; nothing was shown to hold; no failure, no disagreement)"
                   % no_pass_text(dict((o, sum(t["outcomes"].get(o, 0) for t in totals.values()))
                                       for o in set(o for t in totals.values() for o in t["outcomes"]))))
    elif no_pass_arms:
        # an arm that passed nothing verified nothing through that arm, whatever another arm passed
        verdict = ("UNVERIFIED (no execution passed on %s; nothing was shown to hold there; no failure, no "
                   "disagreement)" % "; ".join("%s: %s" % (arm, no_pass_text(o)) for arm, o in no_pass_arms))
    elif stub:
        # A stub run checked the HTTP status only: never PASS (review L2).
        verdict = "UNVERIFIED (stub checks: HTTP status only): no failed status, no disagreement between the arms"
    else:
        verdict = "PASS (no failed execution, no finding, no disagreement between the arms)"
    if inconclusive and not incomplete and not verdict.startswith("INCONCLUSIVE (no execution passed"):
        verdict += "; %d inconclusive execution(s) (not counted as pass)" % inconclusive
    return verdict, {"incomplete": incomplete, "loop": loop, "mixed": mixed, "stub": stub, "untested": untested,
                     "totals": totals, "no_pass_arms": no_pass_arms}


EXPLORATORY_LABEL = "Exploratory verdict (not the deploy decision, which is the gate's %s)"


def exploratory_label(where="above"):
    """The label of the ordinary verdict beside a deploy gate's: the gate is the decision (|where|: the report
    prints the gate above the verdict, the results document's cover beside it)."""
    return EXPLORATORY_LABEL % where


def verdict_headline(verdict, gated=False, where="above"):
    """The verdict as the report's headline prints it: the verdict itself, or, for a run with a deploy gate's
    verdict (|gated|: the report leads with the gate), labelled exploratory. The results document prints the same
    headline on its cover (where="beside")."""
    return ("%s: %s" % (exploratory_label(where), verdict)) if gated else verdict


def _shown(value):
    """|value| as md() shows it before escaping: None as '' (a gate.json field that is missing)."""
    return "" if value is None else value


def gate_headline_parts(verdict):
    """(word, how, summary) of a deploy gate's verdict (gate.py's gate.json): 'GATE PASS' | 'GATE FAIL'; 'policy
    <name> v<version>, sha256 <short>, judged <utc>'; 'every planned and required execution recorded, nothing
    blocks' | '<n> problem(s): <kind> <n>, ...' in gate.KINDS order. The report's gate section and the results
    document's cover print these. A missing field (None) prints empty (_shown), as report.md() prints it."""
    probs = verdict.get("problems") or []
    pol = verdict.get("policy") or {}
    by_kind = {}
    for p in probs:
        by_kind[p.get("kind")] = by_kind.get(p.get("kind"), 0) + 1
    word = "GATE %s" % ("PASS" if verdict.get("passed") else "FAIL")
    how = "policy %s v%s, sha256 %s, judged %s" % (_shown(pol.get("name")), _shown(pol.get("policy_version")),
                                                   short(pol.get("sha256")), _shown(verdict.get("judged_utc")))
    summary = ("every planned and required execution recorded, nothing blocks" if verdict.get("passed") else
               "%d problem(s): %s" % (len(probs), ", ".join("%s %d" % (_shown(k), n) for k, n in sorted(
                   by_kind.items(), key=lambda kv: gate_mod.KINDS.index(kv[0]) if kv[0] in gate_mod.KINDS else 99))))
    return word, how, summary


def verdict_section(run, comp, rows, samples, reevaluated=False, missing=None, depth=None, manifest=None,
                    cases=None, gated=False):
    """The verdict; |missing| (gate.missing_planned: the planned executions no arm recorded, and the plan's source)
    makes it INCOMPLETE, the MCP arm's included; |depth| (agreement_depth) splits the MCP agreement; an area or route
    x arm where nothing passed and something was inconclusive (untested_groups, over |manifest|'s routes) keeps it
    from PASS; so does a run, or an arm, with executions and no pass (UNVERIFIED: its marks, not implemented via MCP
    or service-side by design, showed nothing to hold); an MCP execution inconclusive by a limit |cases| (the table)
    states for its scenario is listed apart (stated_mcp_limit). |gated|: the run has a deploy gate's verdict, which
    the report leads with; this one is labelled exploratory (the headline text: verdict_headline). Every
    nonzero count of what is not implemented yet via MCP in the headline is red (mark_nim_counts), and the
    attention line (attention_line) follows it."""
    s = comp["summary"]
    lines = ["# Verdict", ""]
    verdict, facts = verdict_text(run, comp, rows, reevaluated, missing, manifest, cases)
    totals, loop, mixed, stub, untested = (facts["totals"], facts["loop"], facts["mixed"], facts["stub"],
                                           facts["untested"])
    scen = dict((x.get("id"), x) for x in (manifest or {}).get("scenarios", []) if isinstance(x, dict))
    stated = {}
    for r in rows:
        cid = stated_mcp_limit(r, scen.get(r.get("scenario")), cases, comp)
        if cid:
            stated.setdefault((r.get("scenario"), cid), []).append(r.get("repeat"))
    lines.append("**%s**" % mark_nim_counts(md(verdict_headline(verdict, gated))))
    lines.append("")
    attention = attention_line(rows, manifest, comp)
    if attention:
        lines.append(attention)
        lines.append("")
    if mixed and "more than one deployment" not in verdict:
        lines.append("- Mixed releases: the run spans more than one deployment of the host: %s." % md(mixed))
    if loop_problems_text(loop):
        lines.append("- The %s loop: %s." % (md(loop.get("mode", "paired")), md(loop_problems_text(loop))))
    for arm in comp["arms"]:
        t = totals.get(arm, {})
        o = t.get("outcomes", {})
        lines.append("- %s: %d execution(s): %d pass, %d fail, %d error, %d finding, %d inconclusive, %s, %d "
                     "service-side by design." % (
                         md(arm), t.get("executions", 0), o.get("pass", 0), o.get("fail", 0), o.get("error", 0),
                         o.get("finding", 0), o.get("inconclusive", 0), nim_count(o.get("not_implemented_mcp", 0)),
                         o.get("service_by_design", 0)))
        if o.get(cmp.OTHER_OUTCOME):
            lines[-1] = lines[-1][:-1] + "; %d with an outcome outside the vocabulary (read as failures)." % (
                o[cmp.OTHER_OUTCOME])
    ran_mcp = set(r.get("scenario") for r in rows if r.get("arm") == "mcp")
    marked = {}
    for sid in sorted(ran_mcp, key=str):
        for name, m in sc.mcp_check_marks(scen.get(sid)).items():
            marked.setdefault(m.get("mark"), []).append((sid, name))
    if marked:
        # what the MCP cannot show at property level: set aside by the manifest's marks, listed under Agreement
        lines.append("- Property-level marks of the manifest on the MCP arm (set aside where read not_exposed, unless "
                     "they cover every property check of the execution; listed under Agreement): %s." % "; ".join(
                         (nim if mark == "not_implemented" else str)("%d check(s) of %d scenario(s) %s" % (
                             len(marked[mark]), len(set(x[0] for x in marked[mark])), MARK_WORDS.get(mark, mark)))
                         for mark, _ in sc.MCP_CHECK_MARK_LABELS if marked.get(mark)))
    outside = s.get("outcomes_outside_vocabulary") or []
    if outside:
        lines.append("- Outcomes outside the vocabulary (%s), read as failures: %s%s." % (
            md(", ".join(sc.OUTCOMES)), md("; ".join("%s %s r%s %r%s" % (
                x.get("arm"), x.get("scenario"), x.get("repeat"), x.get("outcome"),
                (" (the merge read it as %s)" % x["read_as"]) if x.get("read_as") else "") for x in outside[:20])),
            " ..." if len(outside) > 20 else ""))
    if untested and not verdict.startswith("INCONCLUSIVE (no execution passed in"):
        lines.append("- No execution passed in %d area / route x arm group(s) (only inconclusive ones, beside marks): "
                     "%s." % (len(untested), md(untested_text(untested))))
    if stated:
        lines.append("- Inconclusive via MCP by a limit the by-design table states for the scenario (not counted as "
                     "untested): %s." % md("; ".join("%s mcp r%s (%s)" % (sid, ",".join(str(x) for x in reps), cid)
                                                     for (sid, cid), reps in sorted(stated.items()))))
    ibr = s.get("inconclusive_by_reason") or {}
    if any(ibr.values()):
        # an inconclusive execution tested nothing: by its class (undecidable plausibility, a host state, the
        # checks), so a run with many such lines does not read clean
        lines.append("- Inconclusive executions by reason: %s." % "; ".join(
            "%s: %s" % (md(arm), ", ".join("%d %s" % (n, md(why)) for why, n in v.items()))
            for arm, v in ibr.items() if v))
    for arm, info in sorted((((run.get("reevaluation") or {}).get("renormalised")) or {}).items()):
        lines.append("- arms/%s: its answers were re-normalised before this re-evaluation (%d rebuilt, %d changed; "
                     "see Re-evaluation): not the answers the run recorded." % (
                         md(arm), (info.get("lines") or {}).get("rebuilt", 0), (info.get("lines") or {}).get(
                             "changed", 0)))
    ma = s["mcp_agreement"]
    if ma.get("compared"):
        dc = depth_counts(depth or {})
        missing_mcp = sum(1 for k in ((missing or ([], None))[0]) if k[1] == "mcp")
        lines.append("- MCP agreement: %d agree (%d fully checked, %d partial, %d partial with a non-pass aggregate, "
                     "%d on failure, %d generic checks only, %d refusals; see Agreement), %d not comparable (no check "
                     "both arms evaluated), %d disagree; service-side by design %d verified, %d unverified; %s; %d "
                     "scenario(s) without an MCP execution, %d planned MCP execution(s) not recorded." % (
                         ma.get("agree", 0), dc["full"], dc["partial"], dc["partial_non_pass"], dc["on_failure"],
                         dc["generic_only"], dc["refusal"], ma.get("not_comparable", 0), ma.get("disagree", 0),
                         ma.get("by_design_verified", 0), ma.get("by_design_unverified", 0),
                         nim_count(ma.get("not_implemented_mcp", 0)), ma.get("missing", 0), missing_mcp))
    if stub:
        lines.append("")
        lines.append("**Checks: the stub ran (HTTP status, error code and text only); every other expected check "
                     "reads not_exposed. This run is not a verified run.**")
    lf = load_figures(samples)
    lines.append("")
    lines.append(load_caveat(lf))
    lines.append("")
    return "\n".join(lines)


def gate_section(run_dir, reevaluated=False, max_items=200, lines=None):
    """The deploy gate's verdict (gate.py's RUN_DIR/gate.json or gate.reevaluated.json) when the run was gated:
    its headline and every problem grouped by kind, each with the wall of the execution it names (|lines|, the run's
    executions; fmt_time's mark) and the execution times its detail quotes marked (mark_quoted_times); '' when it was
    not."""
    verdict = sc.read_json(os.path.join(run_dir, gate_mod.gate_file(reevaluated)))
    if not isinstance(verdict, dict) or "problems" not in verdict:
        return ""
    probs = verdict.get("problems") or []
    pol = verdict.get("policy") or {}
    stale = gate_mod.stale_inputs(run_dir, verdict, reevaluated)
    if stale:
        # the run changed after the verdict (a later invocation without --gate): it no longer describes this run
        return "\n".join(["# Deploy gate", "",
                          "**GATE verdict STALE**: the verdict of %s (policy %s, %s) judged other contents of %s than "
                          "this run has now; run gate.py on the run again for a current verdict." % (
                              md(verdict.get("judged_utc")), md(pol.get("name")),
                              "PASS" if verdict.get("passed") else "FAIL", md(", ".join(stale))), ""])
    word, how, summary = gate_headline_parts(verdict)
    out = ["# Deploy gate", "", "**%s** (%s): %s." % (md(word), md(how), md(summary)), ""]
    by_key = dict((sc.execution_key(r), r) for r in lines or [])
    rows = []
    for p in probs[:max_items]:
        line = by_key.get((p.get("scenario"), p.get("arm"), p.get("repeat")))
        rows.append([md(p.get("kind")), md(p.get("scenario") or "-"), md(p.get("arm") or "-"),
                     md("-" if p.get("repeat") is None else p.get("repeat")),
                     fmt_time(line.get("wall_ms")) if line else "-", mark_quoted_times(md(p.get("detail"), 240))])
    if rows:
        out.append(table(["kind", "scenario", "arm", "rep", "wall ms", "detail"], rows,
                         align=["l", "l", "l", "l", "r", "l"]))
    for n in verdict.get("notes") or []:
        out.append("- %s" % md(n, 300))
    out.append("")
    return "\n".join(out)


def load_errors_text(lf):
    """'' or how many samples failed, with the first error."""
    if not lf.get("n_error"):
        return ""
    return "%d of %d samples failed (first: %s)" % (lf["n_error"], lf["n"], str(lf.get("first_error"))[:160])


def load_caveat(lf):
    if lf["n"] == 0:
        return ("Load: no resource samples in this run (the sampler did not run), so the host load during the "
                "timings is unknown; read the timings with that caveat.")
    # the sampler's figures are contextual telemetry of the host, not a control of the load the timings ran under
    busy = (lf["sra_mean"] or 0) > LOAD_CAVEAT_SRA_CPU_PCT or (lf["load_mean"] or 0) > LOAD_CAVEAT_LOAD1
    figures = "contextual telemetry: sra_all CPU mean %s %% (max %s %%), host load1 mean %s (max %s) over %d samples" % (
        fmt_num(lf["sra_mean"], 0), fmt_num(lf["sra_max"], 0), fmt_num(lf["load_mean"], 1),
        fmt_num(lf["load_max"], 1), lf["n"])
    errors = load_errors_text(lf)
    if errors:
        figures += "; " + errors
    if busy:
        return ("**Load caveat:** the timings were taken while mex carried the Logan load (%s). They are "
                "labelled as such and are not the timing baseline; the rerun after the Logan backlog drains is."
                % md(figures))
    if lf["sra_mean"] is None or lf.get("n_error", 0) > LOAD_SAMPLE_ERROR_SHARE * lf["n"]:
        # no figure of the sra_all container (or too many failed samples) says nothing about the Logan load: never
        # the positive sentence
        return ("**Load not measured:** %s, so whether mex carried the Logan load during the timings is unknown; "
                "read the timings with that caveat (%s)." % (
                    "no sample states the sra_all container's CPU" if lf["sra_mean"] is None
                    else "more than %d %% of the samples failed" % round(100 * LOAD_SAMPLE_ERROR_SHARE), md(figures)))
    return "Load: mex was not under the Logan load during this run (%s)." % md(figures)


def environment_section(run, comp, rows, samples):
    lines = ["# Environment", ""]
    man = run.get("manifest") or {}
    instance = run.get("server_instance")
    restarted = any("restarted" in w for w in run.get("warnings", []) or [])
    ident = run.get("manifest_identity") or {}
    top = ident.get("manifest") or {}
    flags = []
    if not man.get("sha_verified", True) or top.get("problems"):
        flags.append("a stated hash MISMATCH")
    if ident.get("differs"):
        flags.append("the arms ran different manifests")
    kv = [
        ("run id", md(run.get("run_id"))),
        ("release", md(run.get("release"))),
        ("server_instance", md(instance) + (" (RESTARTED during the run)" if restarted else "")),
        ("index_fp", md(short(run.get("index_fp"), 16) + "..." if run.get("index_fp") else "not stated by the host")),
        ("manifest", "version %s, sha256 %s, relations sha256 %s, %s scenario(s)%s" % (
            md(man.get("manifest_version")),
            md(short(top.get("content_sha256") or man.get("computed_sha256") or man.get("content_sha256"), 16)),
            md(short(top.get("relations_sha256") or man.get("computed_relations_sha256")
                     or man.get("relations_sha256"), 16)),
            md(man.get("scenarios_total")), (" (%s: see the warnings)" % "; ".join(flags)) if flags else "")),
        ("scripts commit", md(short_commit(run.get("scripts_commit")))),
        ("generated", md(sc.iso_utc(time.time()))),
    ]
    orch = orchestrator_session(run)
    loop = orch.get("loop") or {}
    if orch:
        mode = orch.get("mode") or "sequential"
        desc = {"paired": "paired: per scenario MCP cold, direct-local, direct-remote, MCP warm, one arm at a time",
                "sequential": "sequential: one arm after the other"}.get(mode, mode)
        if loop:
            desc += "; " + loop_text(loop)
        kv.append(("mode", md(desc)))
        if orch.get("mcp_cmd"):
            kv.append(("MCP runner", "command template `%s`" % ascii_safe(orch["mcp_cmd"]).replace("`", "'")))
    hcinfo = run.get("host_cache") or comp.get("host_cache") or {}
    offsets = hcinfo.get("clock_offsets_s") or {}
    kv.append(("host cache", "cold = a scenario's first execution on the host in this run, warm = another within "
               "%s s, aged = older (state unknown); mex clock %s" % (
                   fmt_num(hcinfo.get("warm_window_s", sc.HOST_CACHE_WARM_S), 0),
                   ("%+.3f s against the laptop (applied to direct-local)" % offsets["direct-local"])
                   if "direct-local" in offsets else "offset not measured")))
    if comp.get("by_design_table"):
        t = comp["by_design_table"]
        kv.append(("by-design table", "%s case(s), service commit %s, sha256 %s%s" % (
            md(t.get("cases")), md(short(t.get("service_commit"), 8)), md(short(t.get("sha256"), 12)),
            TABLE_ORIGIN.get((comp.get("judged_with") or {}).get("by_design"), ""))))
    kv.append(("comparison judged with", judged_with_text(comp.get("judged_with"))))
    for arm, info in sorted((run.get("arms") or {}).items(), key=lambda kv_: kv_[0]):
        info = info or {}
        load_s, load_e = info.get("host_load_start"), info.get("host_load_end")
        kv.append(("arm %s" % md(arm), "%s; state %s; %s execution(s); checks %s; %s to %s; runner %s; "
                   "runner's load1 at start %s, at end %s" % (
                       md(info.get("base") or "-"), md(info.get("state") or "-"),
                       md(info.get("executions_recorded", "-")), md(checks_name(info) or "-"),
                       md(info.get("started_utc") or "-"), md(info.get("ended_utc") or "-"),
                       md(info.get("runner_host") or "-"), fmt_num(load_s[0] if load_s else None, 1),
                       fmt_num(load_e[0] if load_e else None, 1))))
    lines.append(table(["", "value"], [[k, v] for k, v in kv], align=["l", "l"]))
    warnings = (run.get("warnings") or []) + (comp.get("warnings") or [])
    if warnings:
        lines.append("Warnings:\n")
        for w in warnings:
            lines.append("- %s" % md(w, 400))
        lines.append("")
    skipped = []
    for arm, info in (run.get("arms") or {}).items():
        for sk in (info or {}).get("skipped") or []:
            skipped.append((arm, sk))
    if skipped:
        lines.append("Not run (by design): %s." % "; ".join(
            "%s on %s (%s)" % (md(sk.get("scenario")), md(arm), md(sk.get("reason"), 80)) for arm, sk in skipped[:50]))
        if len(skipped) > 50:
            lines.append("... and %d more." % (len(skipped) - 50))
        lines.append("")
    return "\n".join(lines)


def outcomes_section(comp):
    s = comp["summary"]
    rows = []
    for area in sorted(s["per_area"], key=area_sort_key):
        for arm in comp["arms"]:
            c = s["per_area"][area].get(arm)
            if not c:
                continue
            rows.append([md(area), md(arm), c["pass"], c["fail"], c["error"], c.get("finding", 0),
                         c.get("inconclusive", 0), nim(str(c["not_implemented_mcp"])) if c["not_implemented_mcp"]
                         else c["not_implemented_mcp"], c["service_by_design"], c["mixed"], c["total"]])
    out = ["# Outcomes per area", "",
           "Scenario outcomes over their repeats (any error -> error, else any failure -> fail, else any finding "
           "-> finding, else any inconclusive -> inconclusive). A finding: the host broke a stated bound or SPEC "
           "promise without failing the request; inconclusive: the scenario could not test its property (an "
           "allowed 503, a withheld answer, a stated resource stop, or no check exposed); neither counts as pass.",
           "", table(["area", "arm", "pass", "fail", "error", "finding", "inconcl.", "not impl. MCP", "by design",
                      "mixed", "total"], rows)]
    return "\n".join(out)


def expectation_pairs(scenario, checks):
    """[(check result, its expectation args or None)]: the k-th result of a check name answers the k-th
    expectation of that name (review L3): a scenario may list a check twice with other arguments, and the checks
    may add an implicit one of another name."""
    expects = {}
    for e in (scenario or {}).get("expect") or []:
        expects.setdefault(e.get("check"), []).append(e.get("args"))
    seen, out = {}, []
    for c in checks or []:
        name = c.get("check")
        k = seen.get(name, 0)
        seen[name] = k + 1
        exp = expects.get(name) or []
        out.append((c, exp[k] if k < len(exp) else None))
    return out


def findings_section(rows, manifest, max_items):
    """Every execution with a finding: the host broke a stated bound or SPEC promise without failing."""
    scen = dict((s.get("id"), s) for s in (manifest or {}).get("scenarios", [])) if manifest else {}
    found = [r for r in rows if r.get("outcome") == "finding"]
    out = ["# Findings", ""]
    if not found:
        out.append("No finding (no check reported that the host broke a stated bound without failing).\n")
        return "\n".join(out)
    out.append("%d execution(s) with a finding: the request succeeded, but the host broke a bound or promise it "
               "states (e.g. a walk past its time budget)%s." % (
                   len(found), "" if len(found) <= max_items else "; the first %d are listed" % max_items))
    out.append("")
    for r in found[:max_items]:
        s = scen.get(r.get("scenario")) or {}
        out.append("## %s (%s, repeat %s)" % (md(r.get("scenario")), md(r.get("arm")), md(r.get("repeat"))))
        out.append("")
        if s.get("title"):
            out.append("- title: %s" % md(s.get("title"), 300))
        out.append("- route: %s; http %s; wall %s ms" % (md(r.get("route") or s.get("route") or "-"),
                                                          md(r.get("http")), fmt_time(r.get("wall_ms"))))
        for c, exp in expectation_pairs(s, r.get("checks")):
            if c.get("outcome") != "finding":
                continue
            out.append("- check %s: observed %s" % (md(c.get("check")), mark_quoted_times(md(c.get("detail"), 500))))
            if exp is not None:
                out.append("    - expected: %s" % md(json.dumps(exp, sort_keys=True), 400))
        out.append("")
    return "\n".join(out)


def failures_section(rows, manifest, max_failures):
    scen = dict((s.get("id"), s) for s in (manifest or {}).get("scenarios", [])) if manifest else {}
    # an outcome outside the vocabulary is a failure too, named as recorded (sc.recorded_outcome)
    bad = [r for r in rows if r.get("outcome") in ("fail", "error") or not sc.outcome_known(sc.recorded_outcome(r))]
    out = ["# Failures", ""]
    if not bad:
        out.append("No execution failed or errored.\n")
        return "\n".join(out)
    out.append("%d execution(s) failed or errored%s." % (
        len(bad), "" if len(bad) <= max_failures else "; the first %d are listed (all are in executions.jsonl)"
        % max_failures))
    out.append("")
    for r in bad[:max_failures]:
        s = scen.get(r.get("scenario")) or {}
        recorded = sc.recorded_outcome(r)
        out.append("## %s (%s, repeat %s): %s" % (md(r.get("scenario")), md(r.get("arm")), md(r.get("repeat")),
                                                 md(r.get("outcome")) if sc.outcome_known(recorded) else
                                                 "outcome %s outside the vocabulary (read as %s)" % (
                                                     md(repr(recorded)),
                                                     "error" if r.get("outcome") == "error" else "a failure")))
        out.append("")
        if s.get("title"):
            out.append("- title: %s" % md(s.get("title"), 300))
        out.append("- route: %s; http %s; wall %s ms" % (md(r.get("route") or s.get("route") or "-"),
                                                          md(r.get("http")), fmt_time(r.get("wall_ms"))))
        err = r.get("error") or {}
        if err.get("code") or err.get("text"):
            out.append("- observed error: %s: %s" % (md(err.get("code") or "-"),
                                                     mark_quoted_times(md(err.get("text") or "", 400))))
        for c, exp in expectation_pairs(s, r.get("checks")):
            if c.get("outcome") == "fail" or c.get("outcome") not in sc.CHECK_OUTCOMES:
                out.append("- check %s: observed %s%s" % (
                    md(c.get("check")), mark_quoted_times(md(c.get("detail"), 500)),
                    "" if c.get("outcome") == "fail" else " (outcome %s, read as a failure)" % md(c.get("outcome"))))
                if exp is not None:
                    out.append("    - expected: %s" % md(json.dumps(exp, sort_keys=True), 400))
        out.append("")
    return "\n".join(out)


DEPTH_CLASSES = (
    ("full", "fully checked: every property check compared on every MCP repeat, MCP and reference pass"),
    ("partial", "partial: a property check not compared (not exposed via MCP, or inconclusive), both pass"),
    ("partial_non_pass", "partial, and a non-pass aggregate on the MCP or the reference side"),
    ("on_failure", "on failure: every property check compared, a non-pass aggregate (e.g. both arms fail alike)"),
    ("generic_only", "generic checks only: the scenario has property checks, none compared"),
    ("refusal", "refusal: no property check, the generic checks (status, code, text) are the claim"),
)
JUDGED = ("pass", "fail", "finding")


def agreement_depth(comp, rows, is_property, manifest=None):
    """How much an MCP `agree` establishes, per scenario (VERIFY-B): compare's agree means the checks both arms
    evaluated agree, not that the answers are equal or right. Each agreeing scenario's property checks
    (checks.is_property; the generic ones, e.g. budget_respected which the MCP never sees, do not count) on each MCP
    repeat against the reference arm's same repeat (else its first; an item compare took from the other direct arm,
    reference_fallbacks, from that arm's: compare.reference_checks): {scenario: {class (DEPTH_CLASSES), not_compared
    (the property checks no repeat compared), partly_compared ({check: {compared: [repeats], skipped: [repeats]}}:
    the property checks some repeat compared and another skipped; the readability review's verification of
    2026-10-16: a check one repeat left undecided while the other repeat decided it on both arms is compared on
    that repeat, so it is not among `not_compared`), reference_fails_hidden (property checks the reference fails or
    finds where the MCP does not expose them)}}. The class reads every repeat (a skip on any repeat makes the pair
    partial). An MCP item a property-level mark of |manifest|'s scenario sets aside (suite_common.checks_marked_aside)
    is neither compared nor skipped: the mark states it, and the report lists it."""
    scen = dict((x.get("id"), x) for x in (manifest or {}).get("scenarios", []) if isinstance(x, dict))
    by = {}
    for r in rows:
        by.setdefault(r.get("scenario"), {}).setdefault(r.get("arm"), {})[r.get("repeat")] = r
    out = {}
    for sid, entry in (comp.get("scenarios") or {}).items():
        ma = entry.get("mcp_agreement") or {}
        if ma.get("class") != "agree":
            continue
        ref = (by.get(sid) or {}).get(ma.get("reference")) or {}
        mcp = (by.get(sid) or {}).get("mcp") or {}
        if not ref or not mcp:
            continue
        ref0 = ref[min(ref)]
        # the other direct arm's items where compare took a reference item from it (reference_fallbacks)
        alt = ((by.get(sid) or {}).get(ma.get("reference_alternative")) or {}) if ma.get("reference_fallbacks") else {}
        alt0 = alt[min(alt)] if alt else None
        has_prop, compared, skipped, hidden = False, set(), set(), set()
        by_rep = {}  # check -> {"compared": {repeats}, "skipped": {repeats}}
        marks = sc.mcp_check_marks(scen.get(sid))
        for rep in sorted(mcp, key=lambda k: k if isinstance(k, int) else -1):
            mine = mcp[rep].get("checks") or []
            aside = sc.checks_marked_aside(scen.get(sid), "mcp", mine, marks)
            theirs, _ = cmp.reference_checks(mine, (ref.get(rep) or ref0).get("checks") or [],
                                             (alt.get(rep) or alt0 or {}).get("checks") if alt else None)
            has_prop = has_prop or any(is_property(c.get("check")) for c in mine + theirs)
            ca, cb = {}, {}
            for c in mine:
                ca.setdefault(c.get("check"), []).append(c.get("outcome"))
            for c in theirs:
                cb.setdefault(c.get("check"), []).append(c.get("outcome"))
            for name in set(ca) & set(cb):
                if not is_property(name):
                    continue
                for k, (oa, ob) in enumerate(zip(ca[name], cb[name])):
                    if oa in JUDGED and ob in JUDGED:
                        compared.add(name)
                        by_rep.setdefault(name, {"compared": set(), "skipped": set()})["compared"].add(rep)
                    elif (name, k) in aside:
                        if ob in ("fail", "finding"):
                            hidden.add(name)   # the reference's failure stays named: a mark hides no result
                    else:
                        skipped.add(name)
                        by_rep.setdefault(name, {"compared": set(), "skipped": set()})["skipped"].add(rep)
                        if oa == "not_exposed" and ob in ("fail", "finding"):
                            hidden.add(name)
        agg_ok = ma.get("mcp_outcome") == "pass" and ma.get("reference_outcome") == "pass"
        if not has_prop:
            cls = "refusal"
        elif not compared:
            cls = "generic_only"
        elif skipped and not agg_ok:
            cls = "partial_non_pass"
        elif skipped:
            cls = "partial"
        elif not agg_ok:
            cls = "on_failure"
        else:
            cls = "full"
        def _reps(items):
            return sorted(items, key=lambda k: (not isinstance(k, int), str(k)))
        partly = dict((name, {"compared": _reps(by_rep[name]["compared"]), "skipped": _reps(by_rep[name]["skipped"])})
                      for name in sorted(skipped & compared))
        out[sid] = {"class": cls, "not_compared": sorted(skipped - compared), "partly_compared": partly,
                    "reference_fails_hidden": sorted(hidden)}
    return out


def skipped_on_any_repeat(entry):
    """The property checks of an agreement_depth |entry| some MCP repeat did not compare, sorted: those no repeat
    compared (not_compared) and those another repeat compared (partly_compared). The report's depth table lists
    these, as it did before partly_compared (gr/suite-v6); the results document tells the two apart."""
    return sorted(set(entry.get("not_compared") or []) | set(entry.get("partly_compared") or {}))


def depth_counts(depth):
    counts = dict((k, 0) for k, _ in DEPTH_CLASSES)
    for v in depth.values():
        counts[v["class"]] = counts.get(v["class"], 0) + 1
    return counts


def outcomes_text(value):
    """A check's outcomes on one execution (compare.check_outcomes: a list, one per occurrence) as text."""
    if isinstance(value, list):
        return "/".join(str(x) for x in value) or "none"
    return "none" if value is None else str(value)


def marked_checks_lines(comp, manifest, rows):
    """The report's lines of the property-level marks (suite_common.mcp_check_marks) of the scenarios the run has, per
    mark (MCP_CHECK_MARK_LABELS) and scenario: each marked check with its reason or case, the MCP executions it was
    set aside on (suite_common.checks_marked_aside), those it read not_exposed on without being set aside because the
    marks cover every property check there (suite_common.every_property_marked: a scenario-level mark) and any other
    reading (judged as usual: the mark may be stale); [] for a manifest without marks."""
    scen = dict((x.get("id"), x) for x in (manifest or {}).get("scenarios", []) if isinstance(x, dict))
    reads = {}
    for r in rows or []:
        if r.get("arm") == "mcp":
            checks = r.get("checks") or []
            aside = sc.checks_marked_aside(scen.get(r.get("scenario")), "mcp", checks)
            whole = bool(sc.every_property_marked(scen.get(r.get("scenario")), "mcp", checks))
            for c, key in zip(checks, sc.check_keys(checks)):
                how = "aside" if key in aside else "whole" if whole and c.get("outcome") == "not_exposed" else None
                reads.setdefault((r.get("scenario"), c.get("check")), []).append(
                    (r.get("repeat"), c.get("outcome"), how))
    out = []
    for mark, label in sc.MCP_CHECK_MARK_LABELS:
        items = []
        for sid in (comp or {}).get("order") or []:
            for name, m in sorted(sc.mcp_check_marks(scen.get(sid)).items()):
                if m.get("mark") != mark:
                    continue
                seen = sorted(reads.get((sid, name)) or [], key=lambda x: (str(x[0]), str(x[1])))
                aside = ["r%s" % rep for rep, _, how in seen if how == "aside"]
                whole = ["r%s" % rep for rep, _, how in seen if how == "whole"]
                other = ["r%s %s" % (rep, o) for rep, o, how in seen if how is None]
                why = m.get("reason") if mark == "not_implemented" else "case %s" % m.get("case")
                items.append("%s %s (%s; set aside on %s%s%s)" % (
                    md(sid), md(name), md(why or "-", 200), ", ".join(aside) or "no MCP execution",
                    ("; not set aside, %s: %s" % (md(sc.SCENARIO_LEVEL_MARK), ", ".join(whole))) if whole else "",
                    ("; judged as usual: %s" % ", ".join(other)) if other else ""))
        if items:
            out.append("%s: %s.\n" % (nim(label) if mark == "not_implemented" else label, "; ".join(items)))
    return out


def agreement_section(comp, depth=None, mcp_planned_missing=0, manifest=None, lines=None):
    """The agreement between the arms; |manifest| and |lines| (the executions) for the property-level marks the MCP
    part lists (marked_checks_lines)."""
    s = comp["summary"]
    da = s["direct_agreement"]
    out = ["# Agreement between the arms", "", "## direct-local vs direct-remote", ""]
    if not da.get("compared"):
        out.append("Not compared (only one direct arm ran).\n")
    else:
        out.append("%d scenario(s) compared: %d identical answers (digest), %d with the same check outcomes but "
                   "different answers, %d not comparable (different answers, no check evaluated on both), %d "
                   "differing by timing alone, %d disagreeing; %d check result(s) compared." % (
                       da["compared"], da["identical"], da["same_checks"], da.get("not_comparable", 0),
                       da.get("timing", 0), da["disagree"], da.get("checks_compared", 0)))
        out.append("")
        if da.get("timing"):
            out.append("Timing (counted apart, not a disagreement): only one arm's own answer states a time limitation "
                       "(a stop by time, a 503 deadline), that arm reads inconclusive, no check both evaluated differs "
                       "and the stopped answer is consistent with the complete one (exact counts equal, at_least <= "
                       "exact, a stopped path a prefix of a complete one), or a repeat's determinism check had only "
                       "a time-limited earlier repeat. A systematic speed gap between the arms shows here: %s." % (
                           ", ".join(md(x) for x in da.get("timing_scenarios", [])[:100])))
            out.append("")
        ed = da.get("execution_differs") or []
        if ed:
            out.append("Identical answers whose execution checks differ by arm (budget_respected reads the client's "
                       "wall, deterministic the arm's earlier repeats; each arm's own finding or failure counts in its "
                       "totals, not as a disagreement of the answers): %s.\n" % "; ".join(
                           "%s %s" % (md(sid), "; ".join("r%s %s" % (x.get("repeat"), ", ".join(
                               "%s local %s remote %s" % (md(c["check"]), md(outcomes_text(c["local"])),
                                                          md(outcomes_text(c["remote"])))
                               for c in x.get("checks") or []))
                               for x in (comp["scenarios"][sid].get("direct_agreement") or {}).get(
                                   "execution_differs") or [])) for sid in ed[:100]))
        rows = []
        for sid in da.get("disagreeing", [])[:100]:
            a = comp["scenarios"][sid].get("direct_agreement") or {}
            why = "; ".join("%s: %s vs %s" % (md(d["check"]), md(d["a"]), md(d["b"]))
                            for d in a.get("differing_checks", [])[:4])
            if a.get("exact_counts_differ"):
                # compare's exact-count rule: every check may agree (both outcomes pass), the counts do not
                why = "; ".join(x for x in (md(a.get("exact_counts_detail")
                                               or cmp.exact_count_detail(a["exact_counts_differ"]), 400), why) if x)
            rows.append([md(sid), why or "outcome differs"])
        if rows:
            out.append(table(["scenario", "differing checks or exact counts (local vs remote)"], rows,
                             align=["l", "l"]))
    out.append("## MCP")
    out.append("")
    ma = s["mcp_agreement"]
    if not ma.get("compared"):
        if mcp_planned_missing:
            # the run planned the MCP arm and it recorded nothing (a refusal at start, a crash): never "no arm"
            out.append("No MCP execution was recorded although the run planned %d (see the Verdict: INCOMPLETE).\n"
                       % mcp_planned_missing)
        else:
            out.append("No MCP arm in this run.\n")
        return "\n".join(out)
    classes = [k for k in ("agree", "not_comparable", "disagree", "not_implemented_mcp", "service_by_design",
                           "missing", "no_direct_reference") if ma.get(k)]
    out.append("%d scenario(s): %s." % (ma["compared"], ", ".join(
        (nim if k == "not_implemented_mcp" else str)("%d %s" % (ma[k], md(k))) for k in classes)))
    out.append("")
    out.append("Agreement is decided on the checks both arms evaluated: %d check result(s) compared, %d skipped "
               "as not exposed via MCP, %d as inconclusive on either side%s. `agree` needs at least one compared "
               "check; a scenario with none is "
               "`not comparable`. Service-side by design: %d pair(s) verified (each side the table case has was "
               "checked and held), %d unverified (a side stated nothing to compare, or the MCP runner's own "
               "mark).\n" % (ma.get("checks_compared", 0), ma.get("checks_not_exposed", 0),
                              ma.get("checks_inconclusive", 0),
                              (", %d set aside by a property-level mark of the manifest (listed below)"
                               % ma["checks_marked"]) if ma.get("checks_marked") else "",
                              ma.get("by_design_verified", 0), ma.get("by_design_unverified", 0)))
    if depth:
        counts = depth_counts(depth)
        out.append("What an `agree` establishes: agreement of the checks both arms evaluated, not equal or right "
                   "answers (two equally wrong answers agree; a checker defect is shared by every arm). The %d "
                   "agreeing scenario(s) by how much of their claim was compared (property checks; the generic ones, "
                   "e.g. budget_respected, which the MCP never sees, do not count):\n" % sum(counts.values()))
        rows = []
        for key, label in DEPTH_CLASSES:
            sids = [sid for sid in comp["order"] if (depth.get(sid) or {}).get("class") == key]
            if key == "full":
                shown = "-"
            else:
                shown = ", ".join("%s%s" % (md(x), (" (%s)" % md(", ".join(skipped_on_any_repeat(depth[x]))))
                                            if skipped_on_any_repeat(depth[x]) and key.startswith("partial") else "")
                                  for x in sids[:60]) + (" ..." if len(sids) > 60 else "")
            rows.append([md(label), str(len(sids)), shown or "-"])
        out.append(table(["agree, by depth", "n", "scenarios (property checks not compared)"], rows,
                         align=["l", "r", "l"]))
        hidden = [(sid, depth[sid]["reference_fails_hidden"]) for sid in comp["order"]
                  if (depth.get(sid) or {}).get("reference_fails_hidden")]
        if hidden:
            out.append("Agreeing although the reference arm fails a check the MCP does not expose (the MCP cannot "
                       "show it holds): %s.\n" % "; ".join("%s: %s" % (md(sid), md(", ".join(c))) for sid, c in hidden))
    rows = []
    for sid in comp["order"]:
        a = comp["scenarios"][sid].get("mcp_agreement") or {}
        if a.get("class") != "disagree":
            continue
        side = a.get("side")
        detail = {"direct": a.get("direct_detail"), "mcp": a.get("mcp_detail")}.get(side, a.get("detail"))
        limit = 200
        if a.get("exact_counts_differ"):
            # compare's exact-count rule: the counts it compared, whichever side it blames (a comparison of an
            # earlier compare.py states them in `detail` only)
            detail, limit = a.get("detail") or cmp.exact_count_detail(a["exact_counts_differ"]), 400
        if not detail and a.get("differing_checks"):
            detail = "; ".join("%s: mcp %s vs direct %s" % (d["check"], d["a"], d["b"])
                               for d in a["differing_checks"][:4])
        rows.append([md(sid), md(a.get("case") or "-"), md(a.get("side") or "checks"), md(detail or "-", limit)])
    if rows:
        out.append(table(["scenario", "by-design case", "side", "detail"], rows, align=["l", "l", "l", "l"]))
    nc = ma.get("not_comparable_scenarios") or []
    if nc:
        out.append("Not comparable: %s.\n" % ", ".join(md(x) for x in nc[:200]))
    unverified = [sid for sid in comp["order"]
                  if (comp["scenarios"][sid].get("mcp_agreement") or {}).get("class") == "service_by_design"
                  and not (comp["scenarios"][sid].get("mcp_agreement") or {}).get("verified")]
    if unverified:
        out.append("Service-side by design, unverified: %s.\n" % ", ".join(md(x) for x in unverified[:200]))
    by_table = [(sid, (comp["scenarios"][sid].get("mcp_agreement") or {}).get("named_by_table_only"))
                for sid in comp["order"]]
    by_table = [(sid, ids) for sid, ids in by_table if ids]
    if by_table:
        # SD-030: compare matched a case the manifest does not name (only the table's scenario_ids do); the gate
        # accepts a by-design outcome only for a case the manifest names
        out.append("By design only through the table's scenario_ids (the manifest does not name the case, and the "
                   "deploy gate does not accept such an outcome): %s.\n" % "; ".join(
                       "%s (%s)" % (md(sid), md(", ".join(ids))) for sid, ids in by_table[:200]))
    nims = [sid for sid in comp["order"]
            if (comp["scenarios"][sid].get("mcp_agreement") or {}).get("class") == "not_implemented_mcp"]
    if nims:
        out.append("%s: %s.\n" % (nim("Not implemented yet via MCP"), ", ".join(md(x) for x in nims[:200])))
    out.extend(marked_checks_lines(comp, manifest, lines))
    return "\n".join(out)


def relations_section(comp):
    rels = comp.get("relations") or []
    out = ["# Relations between scenarios", ""]
    if not rels:
        out.append("The manifest states no relation.\n")
        return "\n".join(out)
    arms = comp["arms"]
    per_arm = ((comp["summary"].get("relations") or {}).get("per_arm")) or {}
    impl = (comp["summary"].get("relations") or {}).get("implementation") or "-"
    out.append("%d relation(s), each evaluated per arm on the answers of the first repeat of both scenarios there "
               "(semantics: lane S1's README; evaluated by %s): %s." % (
                   len(rels), md(impl), "; ".join("%s %d pass, %d fail, %d inconclusive, %d not exposed, %d not "
                                                  "evaluated" % (md(arm), c.get("pass", 0), c.get("fail", 0),
                                                                 c.get("inconclusive", 0), c.get("not_exposed", 0),
                                                                 c.get("not_evaluated", 0))
                                                  for arm, c in sorted(per_arm.items())) or "none ran"))
    out.append("")
    rows = []
    for r in rels:
        rows.append([md(r.get("id")), md(r.get("type")), md(r.get("a")), md(r.get("b"))]
                    + [md((r["arms"].get(arm) or {}).get("outcome") or "-") for arm in arms])
    out.append(table(["relation", "type", "a", "b"] + [md(a) for a in arms], rows,
                     align=["l"] * (4 + len(arms))))
    bad = [(r, arm, v) for r in rels for arm, v in sorted(r["arms"].items()) if v.get("outcome") != "pass"]
    if bad:
        out.append("Not passed:\n")
        for r, arm, v in bad[:200]:
            out.append("- %s on %s: %s: %s" % (md(r.get("id")), md(arm), md(v.get("outcome")), md(v.get("detail"), 300)))
        out.append("")
    return "\n".join(out)


# ---------------------------------------------------------------------------------------------------------------
# chapter G: the graphlet tools compared five ways (manifest v8; compare_graphlets.py, graphlet_service_checks.py)

GRAPHLET_WAYS = (("W1", "dl+lib: direct-local's graphlet, the library at the deployed commit"),
                 ("W2", "dr+lib: direct-remote's graphlet, the same library"),
                 ("W3", "mcp: the service's graphlet tools (its vendored library)"),
                 ("W4", "mcp+lib: the MCP's export, the library at the same commit"),
                 ("W5", "oracle: the suite's own reading of the full twin and the MGT tokens, no library"))
GRAPHLET_PRECHECKS = ("body", "envelope", "twin", "strategy", "preconditions", "library_commit")


def graphlet_data(run_dir, manifest, reuse=False):
    """The run's graphlet comparison when its manifest has graphlet scenarios (area G), else None: a run of v4..v7
    has no chapter G and its report is unchanged. Recomputed from the run's files in memory every time
    (compare_graphlets.load_or_compute, VERIFY-C VC3); |reuse| (--reuse): graphlet_compare.json when every input it
    read is current."""
    import graphlet_precheck as gp
    if not gp.graphlet_scenarios(manifest if isinstance(manifest, dict) else {}):
        return None
    import compare_graphlets as cg
    try:
        return cg.load_or_compute(run_dir, manifest, reuse)
    except Exception as e:   # the chapter says so; the rest of the report stands
        return {"error": "%s: %s" % (type(e).__name__, e)}


def graphlet_verdict(g):
    """(word, reason) of area G: FAIL (a scenario failed: a pre-check or a disagreement; any service check failed,
    VERIFY-C VC2), INCONCLUSIVE (a scenario inconclusive or precondition_failed, a required service check not passed,
    an exercised one inconclusive), PASS; ERROR when the comparison could not be computed. The reason counts the
    service checks by outcome, n/a included (a check whose precondition held on no call).
    (The declared file's set-aside checks are the gate's: the report's verdict does not read the policy.)"""
    if g.get("error"):
        return "ERROR", "the graphlet comparison could not be computed: %s" % g["error"]
    s = g.get("summary") or {}
    sc_ = s.get("scenarios") or {}
    ops = s.get("ops") or {}
    import graphlet_service_checks as gs
    req = set(g.get("required_service_checks") or [])
    checks = g.get("service_checks") or []
    blocking = gs.blocking_service_checks(checks, req)
    failed = [name for name, why in blocking if why.startswith("fail")]
    open_ = [name for name, why in blocking if not why.startswith("fail")]
    by_outcome = {}
    for c in checks:
        by_outcome[c.get("outcome")] = by_outcome.get(c.get("outcome"), 0) + 1
    reason = "%d scenario(s): %s; op calls: %s; service checks: %s%s" % (
        sum(sc_.values()), ", ".join("%d %s" % (sc_.get(k, 0), k.replace("_", " ")) for k in
                                     ("pass", "fail", "precondition_failed", "inconclusive") if sc_.get(k)) or "none",
        ", ".join("%d %s" % (ops.get(k, 0), k.replace("_", " ")) for k in
                  ("agree", "disagree", "service_by_design", "inconclusive", "precondition_failed") if ops.get(k))
        or "none", ", ".join("%d %s" % (by_outcome[k], k) for k in ("pass", "fail", "inconclusive", "n/a")
                             if by_outcome.get(k)) or "none",
        ("; failed: %s" % ", ".join(failed)) if failed else
        (("; not passed: %s" % ", ".join(open_)) if open_ else ""))
    if sc_.get("fail") or failed:
        return "FAIL", reason
    if sc_.get("inconclusive") or sc_.get("precondition_failed") or open_:
        return "INCONCLUSIVE", reason
    if not sum(sc_.values()):
        return "EMPTY", reason
    return "PASS", reason


def graphlet_headline(g):
    word, reason = graphlet_verdict(g)
    return "**Area G (graphlet tools, five ways): %s** (%s; chapter G)." % (md(word), md(reason))


def _ways_text(per_way):
    return ", ".join("%s %d" % (w, per_way[w]) for w in ("W1", "W2", "W3", "W4", "W5") if per_way.get(w)) or "-"


def graphlet_section(g, max_items=200):
    """Chapter G: per tool its calls, ways and outcomes; every disagreement with both answers compactly; the
    pre-checks; the service checks; the correctness coverage (W5's checked fields against the consistency-only
    ones); the by-design pairs, the drift lane and what the comparison noted."""
    out = ["# G. Graphlet tools, compared five ways", ""]
    if g.get("error"):
        out.append("The graphlet comparison could not be computed: %s\n" % md(g["error"], 300))
        return "\n".join(out)
    out.append(graphlet_headline(g))
    out.append("")
    out.append("The ways: %s. Pairs: W1-W2 isolates the transport, W1-W4 the fetch path, W3-W4 the service's wrapper, "
               "each against W5 the library against an independent reading (on W5's checked fields only). Answers "
               "are compared whole, after the service's own fields are stripped, handles mapped to symbols and "
               "cursors removed; errors by code and structured fields." % "; ".join(
                   "%s %s" % (w, md(t)) for w, t in GRAPHLET_WAYS))
    out.append("")
    scen = g.get("scenarios") or {}
    # G.1 per tool
    per_tool, per_way = {}, {}
    for sid, sv in scen.items():
        for oid, e in (sv.get("ops") or {}).items():
            t = per_way.setdefault(e.get("tool"), {})
            for w, st in (e.get("ways") or {}).items():
                if st in ("answered", "error", "selector_error"):
                    t[w] = t.get(w, 0) + 1
    pt = (g.get("summary") or {}).get("per_tool") or {}
    rows = []
    for tool in sorted(pt):
        c = pt[tool]
        rows.append([md(tool), str(c.get("calls", 0)), md(_ways_text(per_way.get(tool) or {})), str(c.get("agree", 0)),
                     str(c.get("disagree", 0)), str(c.get("service_by_design", 0)), str(c.get("not_exposed", 0)),
                     str(c.get("inconclusive", 0)), str(c.get("precondition_failed", 0)),
                     str(c.get("smoke", 0))])
    out.append("## G.1 Per tool\n")
    out.append(table(["tool", "calls", "records per way", "agree", "disagree", "by design", "W3 not exposed",
                      "inconclusive", "precondition failed", "smoke"], rows,
                     align=["l", "r", "l", "r", "r", "r", "r", "r", "r", "r"]))
    # G.2 disagreements
    out.append("## G.2 Disagreements\n")
    dis = [(sid, oid, e) for sid, sv in sorted(scen.items()) for oid, e in sorted((sv.get("ops") or {}).items())
           if e.get("outcome") == "disagree"]
    if not dis:
        out.append("None.\n")
    for sid, oid, e in dis[:max_items]:
        out.append("- %s %s (%s): %s" % (md(sid), md(oid), md(e.get("tool")), md(", ".join(e.get("disagreeing")
                                                                                      or []))))
        for p in e.get("pairs") or []:
            if p.get("outcome") != "disagree":
                continue
            a, _, b = p.get("pair", "-").partition("-")
            out.append("    - %s (%s)%s%s" % (md(p.get("pair")), md(p.get("what")),
                                             (": %s" % md(p.get("detail"), 200)) if p.get("detail") else "",
                                             (" at %s: %s %s, %s %s" % (md(p["pointer"]), md(a), md(p.get(a), 100),
                                                                        md(b), md(p.get(b), 100)))
                                             if p.get("pointer") else ""))
            for w in (a, b):
                ans = (p.get("answers") or {}).get(w)
                if ans:
                    out.append("        - %s: %s" % (md(w), md(ans, 240)))
        for c in e.get("cuts") or []:
            if c.get("outcome") == "disagree":
                out.append("    - %s-uncut: %s%s" % (md(c["way"]), md(c.get("detail"), 200),
                                                    (" at %s" % md(c["pointer"])) if c.get("pointer") else ""))
        for r in e.get("by_design") or []:
            if not r.get("verified"):
                out.append("    - by-design %s, %s side: %s" % (md(r.get("case")), md(r.get("side")),
                                                                md(r.get("detail"), 240)))
        if e.get("probe"):
            out.append("    - probe %s: %s" % (md(e["probe"]), md(e.get("reason"), 200)))
    if len(dis) > max_items:
        out.append("- ... %d more (graphlet_compare.json)" % (len(dis) - max_items))
    out.append("")
    other = [(sid, oid, e) for sid, sv in sorted(scen.items()) for oid, e in sorted((sv.get("ops") or {}).items())
             if e.get("outcome") in ("inconclusive", "precondition_failed")]
    if other:
        out.append("Inconclusive and precondition_failed op calls (each blocks the gate):\n")
        for sid, oid, e in other[:max_items]:
            out.append("- %s %s (%s): %s: %s" % (md(sid), md(oid), md(e.get("tool")), md(e.get("outcome")),
                                                 md(e.get("reason"), 240)))
        out.append("")
    # G.3 pre-checks
    out.append("## G.3 Pre-checks\n")
    prow, pfail = [], []
    for sid, sv in sorted(scen.items()):
        pre = sv.get("precheck") or {}
        prow.append([md(sid)] + [md((pre.get(k) or {}).get("outcome") or "-") for k in GRAPHLET_PRECHECKS]
                    + [md(pre.get("outcome") or "-") + (" (%s)" % md(pre["finding"]) if pre.get("finding") else "")])
        pfail.extend((sid, k, pre.get(k)) for k in GRAPHLET_PRECHECKS
                     if (pre.get(k) or {}).get("outcome") not in ("pass", "not_exposed"))
    out.append(table(["scenario"] + list(GRAPHLET_PRECHECKS) + ["scenario outcome"], prow,
                     align=["l"] * (len(GRAPHLET_PRECHECKS) + 2)))
    for sid, k, c in pfail[:max_items]:
        out.append("- %s %s %s: %s" % (md(sid), md(k), md((c or {}).get("outcome")), md((c or {}).get("detail"), 300)))
    if pfail:
        out.append("")
    # G.4 service checks
    out.append("## G.4 The service's own behaviour (positive checks)\n")
    req = set(g.get("required_service_checks") or [])
    srow = [[md(c.get("check")), md(c.get("outcome")), str(c.get("calls", 0)), "yes" if c.get("check") in req else "no",
             md(c.get("detail"), 200)] for c in g.get("service_checks") or []]
    out.append(table(["check", "outcome", "calls", "required", "detail"], srow, align=["l", "l", "r", "l", "l"]))
    # G.5 coverage
    out.append("## G.5 Correctness coverage (W5)\n")
    cov = g.get("coverage") or {}
    out.append("W1 = W2 = W3 = W4 is consistency; W5 is the only correctness check. Per tool, the fields W5 checks, "
               "with the SPEC / DESIGN clause each rests on, and the fields it leaves to consistency (list indices "
               "written *). A field W5 computes with a DESIGN formula from engine values, rather than reading it from "
               "an engine value, is marked \"checked, re-derived from DESIGN\" (DESIGN-v8 §2): it catches slips of "
               "the library's implementation, not a wrong rule.\n")
    crow = []
    for tool in sorted(cov):
        for ptr, c in sorted((cov[tool].get("checked") or {}).items()):
            how = "checked" if not c.get("rederived") else "checked, re-derived from DESIGN" + (
                "" if c["rederived"] == c.get("calls") else " (%d of %d calls)" % (c["rederived"], c.get("calls", 0)))
            crow.append([md(tool), md(ptr), md(how), md("; ".join(c.get("clauses") or [c.get("clause")]), 120),
                         str(c.get("calls", 0))])
        for ptr, n in sorted((cov[tool].get("consistency_only") or {}).items()):
            crow.append([md(tool), md(ptr), "consistency only", "-", str(n)])
    out.append(table(["tool", "field", "W5", "clause", "calls"], crow, align=["l", "l", "l", "l", "r"]))
    no_w5 = sorted(t for t in pt if t not in cov)
    if no_w5:
        out.append("Tools W5 does not read in this run (consistency only): %s.\n" % md(", ".join(no_w5)))
    # G.6 by design, drift, notes
    bd = [(sid, oid, r) for sid, sv in sorted(scen.items()) for oid, e in sorted((sv.get("ops") or {}).items())
          for r in e.get("by_design") or []]
    if bd:
        out.append("## G.6 Declared per-way differences (by-design pairs)\n")
        out.append(table(["scenario", "op", "case", "verified", "detail"],
                         [[md(sid), md(oid), md(r.get("case")), "yes" if r.get("verified") else "no (%s)" % md(
                             r.get("side")), md(r.get("detail"), 200)] for sid, oid, r in bd[:max_items]],
                         align=["l"] * 5))
    d = g.get("drift") or {}
    if d.get("calls"):
        out.append("## G.7 Drift lane (report-only)\n")
        out.append("%d call(s) of the branch head's library, %d compared with W4, %d differ.\n" % (
            d.get("calls", 0), d.get("compared", 0), len(d.get("differ") or [])))
        for x in (d.get("differ") or [])[:max_items]:
            out.append("- %s %s (%s) at %s: W4 %s, drift %s" % (md(x.get("scenario")), md(x.get("op")),
                                                                md(x.get("tool")), md(x.get("pointer")),
                                                                md(x.get("W4"), 100), md(x.get("drift"), 100)))
        out.append("")
    notes = []
    gaps = (g.get("summary") or {}).get("whitelist_gaps") or {}
    if gaps:
        notes.append("keys the library ways answer outside the frozen table (compare_graphlets.ANSWER_KEYS): %s" %
                     "; ".join("%s: %s" % (t, ", ".join(k)) for t, k in sorted(gaps.items())))
    if (g.get("summary") or {}).get("unplanned_records"):
        notes.append("%d record(s) of op calls the manifest does not plan" % g["summary"]["unplanned_records"])
    notes.extend(g.get("warnings") or [])
    if notes:
        out.append("Notes:\n")
        out.extend("- %s" % md(n, 300) for n in notes)
        out.append("")
    return "\n".join(out)


def timing_section(comp):
    s = comp["summary"]
    out = ["# Timing", "",
           "Client wall time (ms; request sent to last byte on an established connection, the connection set-up "
           "not included) and the server's own elapsed time (ms; `timing.elapsed_ms`, or the sum over patterns), "
           "per area, arm and the host's cache state: cold = the scenario's first execution on the host in this "
           "run (whichever arm ran it), warm = another execution of it within 60 s, aged = older (cache state "
           "unknown), not reached = an MCP execution the host never saw (refused at submit, not expressible); warm "
           "after a cut = warm, but the scenario's previous execution ended at the host's deadline "
           "(`host_cache_after_interrupted`): counted apart from the warm ones, which the overheads and the cold vs "
           "warm pairs use.", ""]
    rows = []
    for area in sorted(s["timing"], key=area_sort_key):
        for arm in comp["arms"]:
            t = s["timing"][area].get(arm)
            if not t:
                continue
            for hc, b in sorted((t.get("by_host_cache") or {}).items(), key=lambda kv: HOST_CACHE_ORDER.get(kv[0], 9)):
                a, sv = b["wall"], b["server"]
                if not a["n"] and not sv["n"]:
                    continue  # e.g. MCP marks of scenarios it cannot express: nothing was timed
                rows.append([md(area), md(arm), md(HOST_CACHE_LABEL.get(hc, hc)), a["n"], fmt_time(a["p50"]),
                             fmt_time(a["p90"]), fmt_time(a["max"]), sv["n"], fmt_time(sv["p50"]),
                             fmt_time(sv["p90"]), fmt_time(sv["max"])])
    # The server columns cover only the answers that carry the host's timing (review N9: their own n).
    out.append(table(["area", "arm", "host cache", "n", "wall p50", "wall p90", "wall max", "server n",
                      "server p50", "server p90", "server max"], rows))
    rows = []
    for arm, cw in sorted(s["cold_warm"].items()):
        rows.append([md(arm), cw["n"], fmt_ms(cw["delta_ms"]["p50"]), fmt_ms(cw["delta_ms"]["p90"]),
                     fmt_num(cw["ratio"]["p50"], 2), fmt_num(cw["ratio"]["p90"], 2)])
    if rows:
        out.append("## Cold vs warm")
        out.append("")
        out.append("A scenario's host-cold execution and its next host-warm repeat on the same arm (warm minus "
                   "cold; warm / cold):\n")
        out.append(table(["arm", "pairs", "delta p50", "delta p90", "ratio p50", "ratio p90"], rows))
    return "\n".join(out)


def reruns_section(run):
    """The warm-ups run again for a resume (merge_run's run.json `reruns_for_resume`): timing evidence beside the
    target each one ran right before, never an outcome of its own (the warm-up's judged execution is its first)."""
    again = run.get("reruns_for_resume") or []
    if not again:
        return ""
    out = ["## Warm-ups run again for a resume", "",
           "A resume after an interruption between a warm-up and its target ran the warm-up again right before the "
           "target, on every arm that ran the target (suite_common.resume_rewarms). These re-runs are timing "
           "evidence only: the warm-up's outcome is its first execution's.", ""]
    rows = []
    for r in again:
        t = r.get("target_execution")
        rows.append([md(r.get("scenario")), md(r.get("arm")), md(r.get("repeat")), fmt_time(r.get("wall_ms")),
                     md(r.get("http")), md(r.get("rerun_target")),
                     md("repeat %s" % t[2]) if t else "did not run after it",
                     fmt_num(r.get("before_s"), 3) if t else "-"])
    out.append(table(["warm-up", "arm", "repeat", "wall", "http", "target", "target execution",
                      "s before the target"], rows))
    return "\n".join(out)


def overhead_section(comp):
    s = comp["summary"]
    fo = s["front_overhead"]
    out = ["# Overheads", "", "## Front: direct-remote minus direct-local", "",
           "The laptop's network and the front proxy on dnaloc-new, for the same request, on pairs of host-warm "
           "executions only: both executions meet the recent-repeat heuristic (another execution of the scenario "
           "ended on the host within 60 s before); cache residency and foreign load were not controlled. "
           "Each request is timed on an established connection: the TCP and TLS set-up is recorded apart "
           "(below) and is not part of the front cost. The client-side row compares wall minus the server's own "
           "elapsed time on each side; it needs the answer's timing, so error answers and the capabilities "
           "documents count in the wall row only.", ""]
    labels = {"net_ms": "client side (wall - server), warm pairs", "wall_ms": "wall, warm pairs",
              "server_ms": "server elapsed, warm pairs"}
    rows = []
    for k in ("net_ms", "wall_ms", "server_ms"):
        v = (fo.get("warm") or {}).get(k) or {}
        if v.get("n"):
            rows.append([md(labels[k]), v["n"], fmt_ms(v["p50"]), fmt_ms(v["p90"]), fmt_ms(v["max"])])
    out.append(table(["difference (ms)", "n", "p50", "p90", "max"], rows))
    if fo.get("warm_pairs_left_out_after_interrupted"):
        out.append("%d scenario(s) left out of the warm pairs, counted apart: on an arm every warm execution followed "
                   "one cut at the host's deadline (a direct 503 deadline or a deadline the service relayed: "
                   "`host_cache_after_interrupted`); the host had worked on the request, but its row cache may hold "
                   "only part of a complete execution's (the timing table shows them as warm after a cut).\n"
                   % fo["warm_pairs_left_out_after_interrupted"])
    if fo.get("scenarios_without_warm_pair"):
        out.append("%d scenario(s) ran on both direct arms without a host-warm execution on each (the first arm "
                   "had only the cold one); they are not in the front cost.\n" % fo["scenarios_without_warm_pair"])
    net = fo.get("per_area_net") or {}
    rows = []
    for a, v in sorted((fo.get("per_area_wall") or {}).items(), key=lambda kv: area_sort_key(kv[0])):
        if not v["n"]:
            continue
        n_ = net.get(a) or {}
        rows.append([md(a), n_.get("n", 0), fmt_ms(n_.get("p50")), fmt_ms(n_.get("p90")), v["n"],
                     fmt_ms(v["p50"]), fmt_ms(v["p90"])])
    if rows:
        out.append("Per area, warm pairs:\n")
        out.append(table(["area", "client-side n", "p50", "p90", "wall n", "wall p50", "wall p90"], rows))
    rows = [[md(arm), v["n"], fmt_ms(v["p50"]), fmt_ms(v["p90"]), fmt_ms(v["max"])]
            for arm, v in sorted((fo.get("connect_ms") or {}).items()) if v.get("n")]
    if rows:
        out.append("Connection set-up (TCP, and TLS on the public URL), per fresh connection; not in any wall "
                   "time above (one keep-alive connection per runner invocation, i.e. per scenario and arm in "
                   "paired mode):\n")
        out.append(table(["arm", "connections", "connect p50", "connect p90", "connect max"], rows))
    rows = []
    for arm, cv in sorted((fo.get("cold") or {}).items()):
        w, sv = cv.get("wall_ms") or {}, cv.get("server_ms") or {}
        if w.get("n"):
            rows.append([md(arm), w["n"], fmt_time(w["p50"]), fmt_time(w["p90"]), sv.get("n", 0),
                         fmt_time(sv.get("p50")), fmt_time(sv.get("p90"))])
    if rows:
        out.append("Cold executions (a scenario's first on the host; the arm that ran it first), reported apart "
                   "from the warm pairs:\n")
        out.append(table(["arm", "n", "wall p50", "wall p90", "server n", "server p50", "server p90"], rows))
    out.append("## MCP")
    out.append("")
    if not s["mcp_phases"] and not s["mcp_overhead"] and not (s.get("mcp_cold_wall_ms") or {}).get("n"):
        out.append("No MCP phases in this run.\n")
        return "\n".join(out)
    rows = [[md(k), v["n"], fmt_time(v["p50"]), fmt_time(v["p90"]), fmt_time(v["max"])]
            for k, v in sorted(s["mcp_phases"].items()) if v["n"]]
    out.append("Phases of the MCP warm executions (DESIGN section 5: submit, queue, backend as the service "
               "measured it, terminal lag, paging):\n")
    out.append(table(["phase (ms)", "n", "p50", "p90", "max"], rows))
    rows = [[md(k), v["n"], fmt_ms(v["p50"]), fmt_ms(v["p90"]), fmt_ms(v["max"])]
            for k, v in sorted(s["mcp_overhead"].items()) if v["n"]]
    out.append("MCP warm minus each direct arm's warm execution (both host-warm), with and without the queue; "
               "the backend time against each direct arm (executions that reached the host only):\n")
    out.append(table(["difference (ms)", "n", "p50", "p90", "max"], rows))
    cw = s.get("mcp_cold_wall_ms") or {}
    if cw.get("n"):
        out.append("MCP wall time of the cold execution (the scenario's first on the host), apart: n %d, p50 %s, "
                   "p90 %s, max %s ms.\n" % (cw["n"], fmt_time(cw["p50"]), fmt_time(cw["p90"]), fmt_time(cw["max"])))
    return "\n".join(out)


def resources_section(comp, lf=None):
    s = comp["summary"]
    out = ["# Resource telemetry per area (contextual)", "",
           "Contextual telemetry, not a measurement of the suite's requests and not a benchmark: samples of "
           "`docker stats` every 2 s joined to each execution's time window (start - 1 s to end + 2 s, "
           "else the nearest sample within 5 s). The CPU and memory columns are figures of the staging container "
           "while the suite ran (shared with any other traffic it served), not of the suite's request alone; a "
           "2 ms request is attributed a 3 s window. Peak server memory is the answers' own "
           "`work.memory_bytes`. The direct-remote executions join on the laptop's clock against samples stamped "
           "on mex (the offset is in the environment table).", ""]
    rows = []
    for area in sorted(s["resources"], key=area_sort_key):
        for arm in comp["arms"]:
            r = s["resources"][area].get(arm)
            if not r:
                continue
            rows.append([md(area), md(arm), "%d/%d" % (r["executions_with_samples"], r["executions"]),
                         fmt_num(r["peak_cpu_pct"], 0), fmt_gib(r["peak_mem_bytes"]),
                         fmt_num(r["mean_sra_all_cpu_pct"], 0), fmt_num(r["mean_host_load1"], 1),
                         fmt_mib(r["peak_server_memory_bytes"])])
    out.append(table(["area", "arm", "sampled", "peak CPU %", "peak mem GiB", "sra_all CPU %", "load1",
                      "server mem MiB"], rows))
    errors = load_errors_text(lf or {})
    if errors:
        out.append("Failed samples (docker stats did not answer; their figures are null and not in the columns): "
                   "%s.\n" % md(errors))
    return "\n".join(out)


def compare_runs_section(comp, old, old_dir):
    out = ["# Change against %s" % md(old.get("run_id") or os.path.basename(old_dir)), ""]
    if old.get("manifest_sha256") != comp.get("manifest_sha256"):
        out.append("**The manifests differ** (%s vs %s): scenario sets may not match.\n" % (
            md(short(old.get("manifest_sha256"))), md(short(comp.get("manifest_sha256")))))
    out.append("Release %s -> %s.\n" % (md(old.get("release")), md(comp.get("release"))))
    rows = []
    s_new, s_old = comp["summary"], old.get("summary", {})
    areas = sorted(set(s_new["per_area"]) | set(s_old.get("per_area", {})), key=area_sort_key)
    for area in areas:
        for arm in sorted(set(s_new["per_area"].get(area, {})) | set(s_old.get("per_area", {}).get(area, {}))):
            n = s_new["per_area"].get(area, {}).get(arm, {})
            o = s_old.get("per_area", {}).get(area, {}).get(arm, {})
            # Like with like (review N12): the host-warm executions of each run (an earlier run without host
            # cache states shows none).
            tn = warm_wall(s_new, area, arm)
            to = warm_wall(s_old, area, arm)
            p_new, p_old = tn.get("p50"), to.get("p50")
            delta = "-"
            if p_new is not None and p_old:
                delta = "%+.0f%%" % (100.0 * (p_new - p_old) / p_old)
            rows.append([md(area), md(arm), "%s -> %s" % (o.get("pass", "-"), n.get("pass", "-")),
                         "%s -> %s" % (o.get("fail", 0) + o.get("error", 0) if o else "-",
                                       n.get("fail", 0) + n.get("error", 0) if n else "-"),
                         "%s -> %s" % (fmt_time(p_old), fmt_time(p_new)), delta])
    out.append(table(["area", "arm", "pass", "fail+error", "warm wall p50 (ms)", "change"], rows))
    changed = []
    for sid, entry in comp["scenarios"].items():
        old_entry = old.get("scenarios", {}).get(sid)
        if not old_entry:
            continue
        for arm, view in entry["arms"].items():
            ov = old_entry.get("arms", {}).get(arm)
            if ov and ov.get("outcome") != view.get("outcome"):
                changed.append([md(sid), md(arm), md(ov.get("outcome")), md(view.get("outcome"))])
    out.append("Scenarios whose outcome changed:\n")
    out.append(table(["scenario", "arm", "before", "now"], changed[:200], align=["l", "l", "l", "l"]))
    new_only = sorted(set(comp["scenarios"]) - set(old.get("scenarios", {})))
    gone = sorted(set(old.get("scenarios", {})) - set(comp["scenarios"]))
    if new_only or gone:
        out.append("%d scenario(s) only in this run, %d only in the earlier one.\n" % (len(new_only), len(gone)))
    return "\n".join(out)


def warm_wall(summary, area, arm):
    """The wall-time summary of an area's and arm's host-warm executions in a compare.json summary ({} when the
    run has none or predates host_cache)."""
    t = ((summary.get("timing") or {}).get(area) or {}).get(arm) or {}
    return ((t.get("by_host_cache") or {}).get("warm") or {}).get("wall") or {}


def renormalised_text(arm, info):
    """One sentence on an arm whose answers were re-normalised before the re-evaluation (merge_run.renormalisation):
    by what and when where the summary records it (and that the renormaliser's own revision is missing when it
    records only the normaliser's sha256), how many answers were rebuilt and changed, which, and whether the
    re-evaluation of the recorded answers would give other totals (only when an answer changed)."""
    summ = info.get("summary") or {}
    lines = info.get("lines") or {}
    raw = info.get("raw") or {}
    counts = summ.get("counts") or {}
    rev = ", ".join("%s %s" % (md(k), md(short(summ[k], 16))) for k in ("commit", "git_head", "normaliser_sha256",
                                                                      "runner_sha256", "client_sha256", "input_sha256")
                    if summ.get(k))
    own = any(summ.get(k) for k in ("commit", "git_head", "runner_sha256", "client_sha256"))
    changed = summ.get("changed") or lines.get("changed_executions") or []
    n_changed = counts.get("changed", lines.get("changed", 0))
    return ("**arms/%s was re-normalised before this re-evaluation** (%s%s%s%s): %s of %s answer(s) rebuilt from the "
            "recorded tool payloads, %d changed%s; %s not rebuilt%s; %d of %d raw record(s) read keep the recorded "
            "answer as `answer_recorded`. %s" % (
                md(arm),
                "async's suite_mcp_renormalise.py at %s" % md(summ.get("at")) if summ.get("at") else
                "no renormalised.json: the lines' and raw records' marks only",
                ("; " + rev) if rev else "",
                "" if own or not summ else "; its own revision (commit, runner and client sha256) is not recorded",
                ("; digest %s; from %s" % (md(summ.get("digest"), 120), md(summ.get("source"), 160)))
                if summ.get("source") else "",
                counts.get("renormalised", lines.get("rebuilt", 0)), counts.get("lines", "?"), n_changed,
                (" (%s)" % md(", ".join(changed[:40]) + (" ..." if len(changed) > 40 else ""))) if changed else "",
                counts.get("skipped", lines.get("skipped", 0)),
                (" (%s)" % md(", ".join("%s %s" % kv for kv in sorted((summ.get("skipped_by_kind") or {}).items()))))
                if summ.get("skipped_by_kind") else "",
                raw.get("answer_recorded", 0), raw.get("read", 0),
                "The answers judged here are not the ones the run recorded: the recorded arm re-evaluated alone gives "
                "other totals." if n_changed else
                "As no answer changed: the re-evaluation of the recorded answers gives the same totals."))


def reevaluation_section(run_dir, max_items=200, run=None):
    """What a re-evaluation (reevaluate.py's mark, RUN_DIR/reevaluated.json) judged with, whether an arm's answers
    were re-normalised before it (run.reevaluated.json `reevaluation.renormalised`, merge_run), and every execution
    whose outcome or checks it changed: the changed records and, apart, how many changed their outcome."""
    mark = sc.read_json(os.path.join(run_dir, "reevaluated.json")) or {}
    ck = mark.get("checks") or {}
    out = ["# Re-evaluation", ""]
    out.append("The answers this run recorded, judged again offline (no request sent) on %s with the checks of %s "
               "(checks.py sha256 %s, suite_common.py %s, from %s). Outcomes a runner decided without the checks "
               "are kept: error, not implemented yet via MCP, service-side by design, the MCP runner's host-state "
               "rule. `previous` is each earlier repeat's digest and time limitation by these checks; `client` the "
               "direct arms' wall time and the MCP arm's service-measured backend time." % (
                   md(mark.get("reevaluated_utc") or "-"), md(short_commit(ck.get("commit_label")) if ck.get(
                       "commit_label") else "an unversioned directory"), md(short(ck.get("checks_sha256"), 16)),
                   md(short(ck.get("suite_common_sha256"), 16)), md(ck.get("dir") or "-")))
    out.append("")
    rec = mark.get("recorded") or {}
    rows = []
    for arm, c in sorted((mark.get("arms") or {}).items(), key=lambda kv: kv[0]):
        r = rec.get(arm) or {}
        mod = checks_name(r) or "-"
        rows.append([md(arm), md(mod), md(short_commit(r.get("scripts_commit"))), str(c.get("lines", 0)),
                     str(c.get("reevaluated", 0)), str(c.get("kept", 0)), str(c.get("no_answer", 0)),
                     str(c.get("outcome_changed", 0)), str(c.get("checks_changed", 0))])
    out.append(table(["arm", "recorded checks", "runner's commit", "lines", "re-evaluated", "kept", "no answer",
                      "outcome changed", "checks changed"], rows))
    renorm = ((run or {}).get("reevaluation") or {}).get("renormalised") or {}
    for arm, info in sorted(renorm.items()):
        out.append(renormalised_text(arm, info) + "\n")
    changes = mark.get("changes") or []
    if not changes:
        out.append("No outcome and no check changed.\n")
        return "\n".join(out)
    moved = sum(1 for ch in changes if ch.get("before") != ch.get("after"))
    out.append("%d execution record(s) changed: %d changed their outcome, %d only their checks (the outcome "
               "stayed)%s:" % (len(changes), moved, len(changes) - moved,
                               "" if len(changes) <= max_items else
                               " (the first %d listed; all are in reevaluated.json)" % max_items))
    out.append("")
    rows = []
    for ch in changes[:max_items]:
        checks = "; ".join("%s %s->%s" % (c.get("check"), c.get("before") or "-", c.get("after") or "-")
                           for c in ch.get("checks") or [])
        # the detail of a check that now gives the execution's outcome (else a failing one, else the first)
        detail = next((c.get("detail") for c in ch.get("checks") or [] if c.get("after") == ch.get("after")), None) \
            or next((c.get("detail") for c in ch.get("checks") or [] if c.get("after") in ("fail", "finding")), None) \
            or next((c.get("detail") for c in ch.get("checks") or []), "")
        rows.append([md(ch.get("scenario")), md(ch.get("arm")), str(ch.get("repeat")), md(ch.get("before")),
                     md(ch.get("after")) + (" (kept)" if ch.get("kept") else ""), md(checks, 200), md(detail, 220)])
    out.append(table(["scenario", "arm", "rep", "recorded", "re-evaluated", "checks changed", "detail"], rows))
    return "\n".join(out)


def build_report(run_dir, compare_dir=None, max_failures=200, by_design_path=None, reevaluated=False, reuse=False):
    run = sc.read_json(os.path.join(run_dir, sc.reevaluated_name("run.json", reevaluated)), {}) or {}
    rows, _ = sc.read_jsonl(os.path.join(run_dir, sc.reevaluated_name("executions.jsonl", reevaluated)))
    samples, _ = sc.read_jsonl(os.path.join(run_dir, "samples.jsonl"))
    manifest = sc.read_json(os.path.join(run_dir, "manifest.json"))
    comp = load_compare(run_dir, by_design_path, reevaluated)
    title = "MetaGraph suite: %s%s" % (ascii_safe(run.get("run_id") or os.path.basename(os.path.abspath(run_dir))),
                                       " (re-evaluated)" if reevaluated else "")
    parts = ["---", "title: \"%s\"" % title.replace('"', "'"),
             "subtitle: \"/traverse, /resolve, /pattern on %s\"" % ascii_safe(
                 ((manifest or {}).get("index") or {}).get("database") or "refseq33m-experimental"),
             "date: \"%s\"" % sc.iso_utc(time.time()), "geometry: margin=1.8cm", "fontsize: 10pt", "---", "",
             time_legend(), ""]
    is_property = gate_mod._is_property_fn(cmp.judging_checks(run_dir, reevaluated)["dir"])
    depth = agreement_depth(comp, rows, is_property, manifest)
    missing = gate_mod.missing_planned(run_dir, manifest, rows)
    # a gated run's report leads with the gate's verdict (the deploy decision); the ordinary one is exploratory
    gate_text = gate_section(run_dir, reevaluated, max_failures, rows)
    parts.append(gate_text)
    parts.append(verdict_section(run, comp, rows, samples, reevaluated, missing, depth, manifest,
                                 report_cases(comp, run_dir), gated=bool(gate_text)))
    graphlets = graphlet_data(run_dir, manifest, reuse)
    if graphlets is not None:
        parts.append(graphlet_headline(graphlets) + "\n")
    if reevaluated:
        parts.append(reevaluation_section(run_dir, max_failures, run))
    parts.append(environment_section(run, comp, rows, samples))
    parts.append(outcomes_section(comp))
    parts.append(failures_section(rows, manifest, max_failures))
    parts.append(findings_section(rows, manifest, max_failures))
    parts.append(agreement_section(comp, depth, sum(1 for k in missing[0] if k[1] == "mcp"), manifest, rows))
    parts.append(relations_section(comp))
    if graphlets is not None:
        parts.append(graphlet_section(graphlets, max_failures))
    parts.append(timing_section(comp))
    again = reruns_section(run)
    if again:
        parts.append(again)
    parts.append(overhead_section(comp))
    parts.append(resources_section(comp, load_figures(samples)))
    if compare_dir:
        # like with like: the earlier run's re-evaluated set where it has one
        old_re = reevaluated and os.path.exists(os.path.join(compare_dir, "reevaluated.json"))
        old = load_compare(compare_dir, by_design_path, old_re)
        parts.append(compare_runs_section(comp, old, compare_dir + (" (re-evaluated)" if old_re else "")))
    return "\n".join(parts) + "\n"


def to_pdf(md_path, pdf_path):
    """(ok, message): pandoc + xelatex, as DESIGN §6 says, with report_colour.lua for the colour marks (fmt_time)."""
    pandoc = shutil.which("pandoc")
    if not pandoc:
        return False, "pandoc not found: report.md only"
    cmd = [pandoc, md_path, "-o", pdf_path, "--pdf-engine=xelatex", "--lua-filter", LUA_FILTER]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    if proc.returncode != 0:
        return False, "pandoc failed (%d): %s" % (proc.returncode, proc.stderr[-2000:])
    return True, "wrote %s" % pdf_path


def main(argv=None):
    ap = argparse.ArgumentParser(description="Write report.md (and report.pdf) of a suite run.")
    ap.add_argument("run_dir")
    ap.add_argument("--compare", default=None, help="an earlier run directory to show the change against")
    ap.add_argument("--no-pdf", action="store_true")
    ap.add_argument("--max-failures", type=int, default=200)
    ap.add_argument("--by-design", default=None, help="the table compare.py uses when compare.json is missing")
    ap.add_argument("--reevaluated", action="store_true",
                    help="the run's re-evaluated set (reevaluate.py): report.reevaluated.md / .pdf")
    ap.add_argument("--reuse", action="store_true",
                    help="chapter G: reuse graphlet_compare.json when every input it read is current (default: recompute)")
    args = ap.parse_args(argv)
    try:
        text = build_report(args.run_dir, args.compare, args.max_failures, args.by_design, args.reevaluated,
                            args.reuse)
    except cmp.TableError as exc:   # a comparison to compute with a table compare.py refuses (as compare.py: 4)
        sys.stderr.write("report: refusing: %s\n" % exc)
        return 4
    md_path = os.path.join(args.run_dir, sc.reevaluated_name("report.md", args.reevaluated))
    with open(md_path, "w", encoding="ascii", errors="replace") as fh:
        fh.write(text)
    sys.stdout.write("wrote %s\n" % md_path)
    if not args.no_pdf:
        ok, msg = to_pdf(md_path, os.path.join(args.run_dir, sc.reevaluated_name("report.pdf", args.reevaluated)))
        sys.stdout.write(msg + "\n")
        if not ok:
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
