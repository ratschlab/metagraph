#!/usr/bin/env python3
"""Report of a suite run: RUN_DIR/report.md from run.json, executions, samples and compare.json, then
report.pdf via pandoc + xelatex (DESIGN §6).

    report.py RUN_DIR [--compare EARLIER_RUN_DIR] [--no-pdf] [--max-failures 200] [--reevaluated]

--reevaluated: the run's re-evaluated set (reevaluate.py) into report.reevaluated.md / .pdf, with a section on
what the re-evaluation judged with and every execution whose outcome or checks it changed.

Sections: verdict and environment (with the load caveat), outcomes per area and arm, every failure with its
expectation and the observed value, agreement of the direct arms and of the MCP arm, timing per area and arm
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


def table(header, rows, align=None):
    """A pipe table; cells are already escaped strings.

    Numeric columns are right-aligned, others left (unless `align` says). The separator's dashes follow each
    column's widest cell: pandoc sizes the columns of a wide table by them, so a short column (an area) stays
    narrow and a long one (a detail) gets the room, instead of equal widths that wrap every cell.
    """
    if not rows:
        return "(none)\n"
    cols = len(header)
    if align is None:
        align = ["r" if all(_NUMERIC_CELL.match(str(r[i]).strip()) for r in rows) else "l" for i in range(cols)]
    widths = [max([len(str(header[i]))] + [len(str(r[i])) for r in rows]) for i in range(cols)]
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


HOST_CACHE_ORDER = {"cold": 0, "warm": 1, "aged": 2, "not_reached": 3, "unknown": 4}
HOST_CACHE_LABEL = {"cold": "cold", "warm": "warm", "aged": "aged", "not_reached": "not reached",
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


def load_compare(run_dir, by_design_path=None, reevaluated=False):
    path = os.path.join(run_dir, sc.reevaluated_name("compare.json", reevaluated))
    data = sc.read_json(path)
    if data is None:
        table, warnings = cmp.table_for(run_dir, by_design_path, reevaluated)
        data = cmp.compare(run_dir, cmp.load_by_design(table), reevaluated=reevaluated, warnings=warnings)
        sc.write_json_atomic(path, data)
    return data


def area_sort_key(area):
    order = "CRTPLXMQNDS"
    return (order.index(area) if isinstance(area, str) and len(area) == 1 and area in order else len(order), str(area))


# ---------------------------------------------------------------------------------------------------------------
# sections

# the share of failed samples (the sampler's `error`: docker stats failed, e.g. a missing container or daemon) above
# which the load figures do not stand for the run
LOAD_SAMPLE_ERROR_SHARE = 0.10


def load_figures(samples):
    """The sampler's figures over the run: sra_all CPU and host load1 (mean, max, how many samples state each), the
    samples, and those that failed (`error`: their figures are null) with the first error."""
    sra = [s.get("sra_all_cpu_pct") for s in samples if cmp._num(s.get("sra_all_cpu_pct")) is not None]
    load = [s.get("host_load1") for s in samples if cmp._num(s.get("host_load1")) is not None]
    errors = [s.get("error") for s in samples if s.get("error")]
    return {"sra_mean": (sum(sra) / len(sra)) if sra else None, "sra_max": max(sra) if sra else None,
            "load_mean": (sum(load) / len(load)) if load else None, "load_max": max(load) if load else None,
            "n": len(samples), "n_sra": len(sra), "n_load": len(load), "n_error": len(errors),
            "first_error": str(errors[0]) if errors else None}


def stated_mcp_limit(line, scenario, cases, comp=None):
    """The id of the by-design case that states why an MCP execution is inconclusive, else None: the case the scenario
    names (compare.named_cases) that matched the execution (its service_by_design_case, else compare's
    by_design_notes for the repeat) and whose explains_checks name every check of the line that did not pass, each of
    which reads not_exposed (a part of the answer the MCP surface does not expose by design, e.g. capabilities_doc of
    capabilities.traverse_capabilities_document)."""
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
    names = set(e if isinstance(e, str) else (e or {}).get("check") for e in case.get("explains_checks") or [])
    rest = [c for c in line.get("checks") or [] if c.get("outcome") != "pass"]
    if rest and all(c.get("outcome") == "not_exposed" and c.get("check") in names for c in rest):
        return cid
    return None


def report_cases(comp):
    """The by-design table compare judged the run with (compare.json by_design_table), loaded; {} without one."""
    info = (comp or {}).get("by_design_table")
    path = (info or {}).get("path") if isinstance(info, dict) else None
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
                c[r.get("outcome")] = c.get(r.get("outcome"), 0) + 1
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


def verdict_section(run, comp, rows, samples, reevaluated=False, missing=None, depth=None, manifest=None,
                    cases=None):
    """The verdict; |missing| (gate.missing_planned: the planned executions no arm recorded, and the plan's source)
    makes it INCOMPLETE, the MCP arm's included; |depth| (agreement_depth) splits the MCP agreement; an area or route
    x arm where nothing passed and something was inconclusive (untested_groups, over |manifest|'s routes) keeps it
    from PASS; an MCP execution inconclusive by a limit |cases| (the table) states for its scenario is listed apart
    (stated_mcp_limit)."""
    s = comp["summary"]
    lines = ["# Verdict", ""]
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
    dis_direct = s["direct_agreement"].get("disagree", 0)
    dis_mcp = s["mcp_agreement"].get("disagree", 0)
    rel_failed = sum((v or {}).get("fail", 0) for v in ((s.get("relations") or {}).get("per_arm") or {}).values())
    modules = sorted(set(checks_name(info) for info in (run.get("arms") or {}).values() if checks_name(info)))
    # a re-evaluated set is judged by checks.py whatever the runners ran
    stub = "stub" in modules and not reevaluated
    untested = untested_groups(rows, manifest, cases, comp)
    scen = dict((x.get("id"), x) for x in (manifest or {}).get("scenarios", []) if isinstance(x, dict))
    stated = {}
    for r in rows:
        cid = stated_mcp_limit(r, scen.get(r.get("scenario")), cases, comp)
        if cid:
            stated.setdefault((r.get("scenario"), cid), []).append(r.get("repeat"))
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
    elif not executions:
        verdict = "EMPTY (no execution recorded)"
    elif stub:
        # A stub run checked the HTTP status only: never PASS (review L2).
        verdict = "UNVERIFIED (stub checks: HTTP status only): no failed status, no disagreement between the arms"
    else:
        verdict = "PASS (no failed execution, no finding, no disagreement between the arms)"
    if inconclusive and not incomplete and not verdict.startswith("INCONCLUSIVE"):
        verdict += "; %d inconclusive execution(s) (not counted as pass)" % inconclusive
    lines.append("**%s**" % md(verdict))
    lines.append("")
    for arm in comp["arms"]:
        t = totals.get(arm, {})
        o = t.get("outcomes", {})
        lines.append("- %s: %d execution(s): %d pass, %d fail, %d error, %d finding, %d inconclusive, %d not "
                     "implemented yet via MCP, %d service-side by design." % (
                         md(arm), t.get("executions", 0), o.get("pass", 0), o.get("fail", 0), o.get("error", 0),
                         o.get("finding", 0), o.get("inconclusive", 0), o.get("not_implemented_mcp", 0),
                         o.get("service_by_design", 0)))
        if o.get(cmp.OTHER_OUTCOME):
            lines[-1] = lines[-1][:-1] + "; %d with an outcome outside the vocabulary (read as failures)." % (
                o[cmp.OTHER_OUTCOME])
    outside = s.get("outcomes_outside_vocabulary") or []
    if outside:
        lines.append("- Outcomes outside the vocabulary (%s), read as failures: %s%s." % (
            md(", ".join(sc.OUTCOMES)), md("; ".join("%s %s r%s %r" % (x.get("arm"), x.get("scenario"), x.get("repeat"),
                                                                      x.get("outcome")) for x in outside[:20])),
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
                     "both arms evaluated), %d disagree; service-side by design %d verified, %d unverified; %d not "
                     "implemented via MCP; %d scenario(s) without an MCP execution, %d planned MCP execution(s) not "
                     "recorded." % (
                         ma.get("agree", 0), dc["full"], dc["partial"], dc["partial_non_pass"], dc["on_failure"],
                         dc["generic_only"], dc["refusal"], ma.get("not_comparable", 0), ma.get("disagree", 0),
                         ma.get("by_design_verified", 0), ma.get("by_design_unverified", 0),
                         ma.get("not_implemented_mcp", 0), ma.get("missing", 0), missing_mcp))
    if stub:
        lines.append("")
        lines.append("**Checks: the stub ran (HTTP status, error code and text only); every other expected check "
                     "reads not_exposed. This run is not a verified run.**")
    lf = load_figures(samples)
    lines.append("")
    lines.append(load_caveat(lf))
    lines.append("")
    return "\n".join(lines)


def gate_section(run_dir, reevaluated=False, max_items=200):
    """The deploy gate's verdict (gate.py's RUN_DIR/gate.json or gate.reevaluated.json) when the run was gated:
    its headline and every problem grouped by kind; '' when it was not."""
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
    by_kind = {}
    for p in probs:
        by_kind[p.get("kind")] = by_kind.get(p.get("kind"), 0) + 1
    out = ["# Deploy gate", "",
           "**GATE %s** (policy %s v%s, sha256 %s, judged %s): %s." % (
               "PASS" if verdict.get("passed") else "FAIL", md(pol.get("name")), md(pol.get("policy_version")),
               md(short(pol.get("sha256"))), md(verdict.get("judged_utc")),
               "every planned and required execution recorded, nothing blocks" if verdict.get("passed") else
               "%d problem(s): %s" % (len(probs), ", ".join("%s %d" % (md(k), n) for k, n in sorted(
                   by_kind.items(), key=lambda kv: gate_mod.KINDS.index(kv[0]) if kv[0] in gate_mod.KINDS else 99)))),
           ""]
    rows = []
    for p in probs[:max_items]:
        rows.append([md(p.get("kind")), md(p.get("scenario") or "-"), md(p.get("arm") or "-"),
                     md("-" if p.get("repeat") is None else p.get("repeat")), md(p.get("detail"), 240)])
    if rows:
        out.append(table(["kind", "scenario", "arm", "rep", "detail"], rows, align=["l", "l", "l", "l", "l"]))
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
            desc += "; loop %s, %s of %s scenario(s)" % (loop.get("state"), loop.get("completed", "-"),
                                                         loop.get("scenarios", "-"))
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
        kv.append(("by-design table", "%s case(s), service commit %s, sha256 %s" % (
            md(t.get("cases")), md(short(t.get("service_commit"), 8)), md(short(t.get("sha256"), 12)))))
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
                         c.get("inconclusive", 0), c["not_implemented_mcp"], c["service_by_design"], c["mixed"],
                         c["total"]])
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
                                                          md(r.get("http")), fmt_ms(r.get("wall_ms"))))
        for c, exp in expectation_pairs(s, r.get("checks")):
            if c.get("outcome") != "finding":
                continue
            out.append("- check %s: observed %s" % (md(c.get("check")), md(c.get("detail"), 500)))
            if exp is not None:
                out.append("    - expected: %s" % md(json.dumps(exp, sort_keys=True), 400))
        out.append("")
    return "\n".join(out)


def failures_section(rows, manifest, max_failures):
    scen = dict((s.get("id"), s) for s in (manifest or {}).get("scenarios", [])) if manifest else {}
    bad = [r for r in rows if r.get("outcome") in ("fail", "error")]
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
        out.append("## %s (%s, repeat %s): %s" % (md(r.get("scenario")), md(r.get("arm")), md(r.get("repeat")),
                                                 md(r.get("outcome"))))
        out.append("")
        if s.get("title"):
            out.append("- title: %s" % md(s.get("title"), 300))
        out.append("- route: %s; http %s; wall %s ms" % (md(r.get("route") or s.get("route") or "-"),
                                                          md(r.get("http")), fmt_ms(r.get("wall_ms"))))
        err = r.get("error") or {}
        if err.get("code") or err.get("text"):
            out.append("- observed error: %s: %s" % (md(err.get("code") or "-"), md(err.get("text") or "", 400)))
        for c, exp in expectation_pairs(s, r.get("checks")):
            if c.get("outcome") == "fail" or c.get("outcome") not in sc.CHECK_OUTCOMES:
                out.append("- check %s: observed %s%s" % (
                    md(c.get("check")), md(c.get("detail"), 500),
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


def agreement_depth(comp, rows, is_property):
    """How much an MCP `agree` establishes, per scenario (VERIFY-B): compare's agree means the checks both arms
    evaluated agree, not that the answers are equal or right. Each agreeing scenario's property checks
    (checks.is_property; the generic ones, e.g. budget_respected which the MCP never sees, do not count) on each MCP
    repeat against the reference arm's same repeat (else its first; an item compare took from the other direct arm,
    reference_fallbacks, from that arm's: compare.reference_checks): {scenario: {class (DEPTH_CLASSES), not_compared
    (the property checks no repeat compared), reference_fails_hidden (property checks the reference fails or finds
    where the MCP does not expose them)}}."""
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
        for rep in sorted(mcp, key=lambda k: k if isinstance(k, int) else -1):
            mine = mcp[rep].get("checks") or []
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
                for oa, ob in zip(ca[name], cb[name]):
                    if oa in JUDGED and ob in JUDGED:
                        compared.add(name)
                    else:
                        skipped.add(name)
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
        out[sid] = {"class": cls, "not_compared": sorted(skipped), "reference_fails_hidden": sorted(hidden)}
    return out


def depth_counts(depth):
    counts = dict((k, 0) for k, _ in DEPTH_CLASSES)
    for v in depth.values():
        counts[v["class"]] = counts.get(v["class"], 0) + 1
    return counts


def agreement_section(comp, depth=None, mcp_planned_missing=0):
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
                               "%s local %s remote %s" % (md(c["check"]), md(c["local"]), md(c["remote"]))
                               for c in x.get("checks") or []))
                               for x in (comp["scenarios"][sid].get("direct_agreement") or {}).get(
                                   "execution_differs") or [])) for sid in ed[:100]))
        rows = []
        for sid in da.get("disagreeing", [])[:100]:
            a = comp["scenarios"][sid].get("direct_agreement") or {}
            rows.append([md(sid), "; ".join("%s: %s vs %s" % (md(d["check"]), md(d["a"]), md(d["b"]))
                                            for d in a.get("differing_checks", [])[:4]) or "outcome differs"])
        if rows:
            out.append(table(["scenario", "differing checks (local vs remote)"], rows, align=["l", "l"]))
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
    out.append("%d scenario(s): %s." % (ma["compared"], ", ".join("%d %s" % (ma[k], md(k)) for k in classes)))
    out.append("")
    out.append("Agreement is decided on the checks both arms evaluated: %d check result(s) compared, %d skipped "
               "as not exposed via MCP, %d as inconclusive on either side. `agree` needs at least one compared "
               "check; a scenario with none is "
               "`not comparable`. Service-side by design: %d pair(s) verified (each side the table case has was "
               "checked and held), %d unverified (a side stated nothing to compare, or the MCP runner's own "
               "mark).\n" % (ma.get("checks_compared", 0), ma.get("checks_not_exposed", 0),
                              ma.get("checks_inconclusive", 0), ma.get("by_design_verified", 0),
                              ma.get("by_design_unverified", 0)))
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
                shown = ", ".join("%s%s" % (md(x), (" (%s)" % md(", ".join(depth[x]["not_compared"])))
                                            if depth[x]["not_compared"] and key.startswith("partial") else "")
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
        if not detail and a.get("differing_checks"):
            detail = "; ".join("%s: mcp %s vs direct %s" % (d["check"], d["a"], d["b"])
                               for d in a["differing_checks"][:4])
        rows.append([md(sid), md(a.get("case") or "-"), md(a.get("side") or "checks"), md(detail or "-", 200)])
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
    nim = [sid for sid in comp["order"]
           if (comp["scenarios"][sid].get("mcp_agreement") or {}).get("class") == "not_implemented_mcp"]
    if nim:
        out.append("Not implemented yet via MCP: %s.\n" % ", ".join(md(x) for x in nim[:200]))
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


def timing_section(comp):
    s = comp["summary"]
    out = ["# Timing", "",
           "Client wall time (ms; request sent to last byte on an established connection, the connection set-up "
           "not included) and the server's own elapsed time (ms; `timing.elapsed_ms`, or the sum over patterns), "
           "per area, arm and the host's cache state: cold = the scenario's first execution on the host in this "
           "run (whichever arm ran it), warm = another execution of it within 60 s, aged = older (cache state "
           "unknown), not reached = an MCP execution the host never saw (refused at submit, not expressible).", ""]
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
                rows.append([md(area), md(arm), md(HOST_CACHE_LABEL.get(hc, hc)), a["n"], fmt_ms(a["p50"]),
                             fmt_ms(a["p90"]), fmt_ms(a["max"]), sv["n"], fmt_ms(sv["p50"]), fmt_ms(sv["p90"]),
                             fmt_ms(sv["max"])])
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
    if fo.get("warm_pairs_after_interrupted"):
        out.append("%d of the warm pairs have an execution that followed one cut at the host's deadline (a direct 503 "
                   "deadline or a deadline the service relayed: host_cache_after_interrupted); the host had worked on "
                   "the request, so it reads warm, but its row cache may hold only part of a complete execution's.\n"
                   % fo["warm_pairs_after_interrupted"])
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
            rows.append([md(arm), w["n"], fmt_ms(w["p50"]), fmt_ms(w["p90"]), sv.get("n", 0), fmt_ms(sv.get("p50")),
                         fmt_ms(sv.get("p90"))])
    if rows:
        out.append("Cold executions (a scenario's first on the host; the arm that ran it first), reported apart "
                   "from the warm pairs:\n")
        out.append(table(["arm", "n", "wall p50", "wall p90", "server n", "server p50", "server p90"], rows))
    out.append("## MCP")
    out.append("")
    if not s["mcp_phases"] and not s["mcp_overhead"] and not (s.get("mcp_cold_wall_ms") or {}).get("n"):
        out.append("No MCP phases in this run.\n")
        return "\n".join(out)
    rows = [[md(k), v["n"], fmt_ms(v["p50"]), fmt_ms(v["p90"]), fmt_ms(v["max"])]
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
                   "p90 %s, max %s ms.\n" % (cw["n"], fmt_ms(cw["p50"]), fmt_ms(cw["p90"]), fmt_ms(cw["max"])))
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
                         "%s -> %s" % (fmt_ms(p_old), fmt_ms(p_new)), delta])
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


def build_report(run_dir, compare_dir=None, max_failures=200, by_design_path=None, reevaluated=False):
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
             "date: \"%s\"" % sc.iso_utc(time.time()), "geometry: margin=1.8cm", "fontsize: 10pt", "---", ""]
    is_property = gate_mod._is_property_fn(cmp.reevaluation_checks_dir(run_dir) if reevaluated else None)
    depth = agreement_depth(comp, rows, is_property)
    missing = gate_mod.missing_planned(run_dir, manifest, rows)
    parts.append(verdict_section(run, comp, rows, samples, reevaluated, missing, depth, manifest, report_cases(comp)))
    parts.append(gate_section(run_dir, reevaluated, max_failures))
    if reevaluated:
        parts.append(reevaluation_section(run_dir, max_failures, run))
    parts.append(environment_section(run, comp, rows, samples))
    parts.append(outcomes_section(comp))
    parts.append(failures_section(rows, manifest, max_failures))
    parts.append(findings_section(rows, manifest, max_failures))
    parts.append(agreement_section(comp, depth, sum(1 for k in missing[0] if k[1] == "mcp")))
    parts.append(relations_section(comp))
    parts.append(timing_section(comp))
    parts.append(overhead_section(comp))
    parts.append(resources_section(comp, load_figures(samples)))
    if compare_dir:
        # like with like: the earlier run's re-evaluated set where it has one
        old_re = reevaluated and os.path.exists(os.path.join(compare_dir, "reevaluated.json"))
        old = load_compare(compare_dir, by_design_path, old_re)
        parts.append(compare_runs_section(comp, old, compare_dir + (" (re-evaluated)" if old_re else "")))
    return "\n".join(parts) + "\n"


def to_pdf(md_path, pdf_path):
    """(ok, message): pandoc + xelatex, as DESIGN §6 says."""
    pandoc = shutil.which("pandoc")
    if not pandoc:
        return False, "pandoc not found: report.md only"
    cmd = [pandoc, md_path, "-o", pdf_path, "--pdf-engine=xelatex"]
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
    args = ap.parse_args(argv)
    text = build_report(args.run_dir, args.compare, args.max_failures, args.by_design, args.reevaluated)
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
