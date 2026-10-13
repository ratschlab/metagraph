#!/usr/bin/env python3
"""Report of a suite run: RUN_DIR/report.md from run.json, executions, samples and compare.json, then
report.pdf via pandoc + xelatex (DESIGN §6).

    report.py RUN_DIR [--compare EARLIER_RUN_DIR] [--no-pdf] [--max-failures 200]

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


def orchestrator_session(run):
    """The orchestrator's last invocation (run_suite.sh), or {}."""
    sessions = ((run.get("orchestrator") or {}).get("sessions")) or []
    return sessions[-1] if sessions else {}


def load_compare(run_dir, by_design_path=None):
    path = os.path.join(run_dir, "compare.json")
    data = sc.read_json(path)
    if data is None:
        data = cmp.compare(run_dir, cmp.load_by_design(by_design_path or cmp.DEFAULT_BY_DESIGN))
        sc.write_json_atomic(path, data)
    return data


def area_sort_key(area):
    order = "CRTPLXMQNDS"
    return (order.index(area) if isinstance(area, str) and len(area) == 1 and area in order else len(order), str(area))


# ---------------------------------------------------------------------------------------------------------------
# sections

def load_figures(samples):
    sra = [s.get("sra_all_cpu_pct") for s in samples if cmp._num(s.get("sra_all_cpu_pct")) is not None]
    load = [s.get("host_load1") for s in samples if cmp._num(s.get("host_load1")) is not None]
    return {"sra_mean": (sum(sra) / len(sra)) if sra else None, "sra_max": max(sra) if sra else None,
            "load_mean": (sum(load) / len(load)) if load else None, "load_max": max(load) if load else None,
            "n": len(samples)}


def verdict_section(run, comp, rows, samples):
    s = comp["summary"]
    lines = ["# Verdict", ""]
    totals = s.get("totals", {})
    failing = sum(t["outcomes"].get("fail", 0) + t["outcomes"].get("error", 0) for t in totals.values())
    findings = sum(t["outcomes"].get("finding", 0) for t in totals.values())
    inconclusive = sum(t["outcomes"].get("inconclusive", 0) for t in totals.values())
    states = dict((arm, (info or {}).get("state")) for arm, info in (run.get("arms") or {}).items())
    incomplete = ["%s (%s)" % (arm, st) for arm, st in sorted(states.items())
                  if st not in (None, "done") and arm != "mcp"]
    loop = orchestrator_session(run).get("loop") or {}
    if loop and loop.get("state") not in (None, "done"):
        incomplete.append("the %s loop (%s)" % (loop.get("mode", "paired"), loop.get("state")))
    dis_direct = s["direct_agreement"].get("disagree", 0)
    dis_mcp = s["mcp_agreement"].get("disagree", 0)
    rel_failed = sum((v or {}).get("fail", 0) for v in ((s.get("relations") or {}).get("per_arm") or {}).values())
    modules = sorted(set((info or {}).get("checks_module") for info in (run.get("arms") or {}).values()
                         if (info or {}).get("checks_module")))
    stub = "stub" in modules
    if incomplete:
        verdict = "INCOMPLETE (not finished: %s)" % ", ".join(incomplete)
    elif failing or dis_direct or dis_mcp or rel_failed:
        verdict = "FAIL (%d failed or errored execution(s), %d direct and %d MCP disagreement(s)%s)" % (
            failing, dis_direct, dis_mcp, (", %d failed relation(s)" % rel_failed) if rel_failed else "")
    elif findings:
        verdict = ("FINDINGS (%d execution(s) where the host broke a stated bound without failing the request; no "
                   "failure, no disagreement)" % findings)
    elif stub:
        # A stub run checked the HTTP status only: never PASS (review L2).
        verdict = "UNVERIFIED (stub checks: HTTP status only): no failed status, no disagreement between the arms"
    else:
        verdict = "PASS (no failed execution, no finding, no disagreement between the arms)"
    if inconclusive and not incomplete:
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
    ma = s["mcp_agreement"]
    if ma.get("compared"):
        lines.append("- MCP agreement: %d agree, %d not comparable (no check both arms evaluated), %d disagree; "
                     "service-side by design %d verified, %d unverified." % (
                         ma.get("agree", 0), ma.get("not_comparable", 0), ma.get("disagree", 0),
                         ma.get("by_design_verified", 0), ma.get("by_design_unverified", 0)))
    if stub:
        lines.append("")
        lines.append("**Checks: the stub ran (HTTP status, error code and text only); every other expected check "
                     "reads not_exposed. This run is not a verified run.**")
    lf = load_figures(samples)
    lines.append("")
    lines.append(load_caveat(lf))
    lines.append("")
    return "\n".join(lines)


def load_caveat(lf):
    if lf["n"] == 0:
        return ("Load: no resource samples in this run (the sampler did not run), so the host load during the "
                "timings is unknown; read the timings with that caveat.")
    busy = (lf["sra_mean"] or 0) > LOAD_CAVEAT_SRA_CPU_PCT or (lf["load_mean"] or 0) > LOAD_CAVEAT_LOAD1
    figures = "sra_all CPU mean %s %% (max %s %%), host load1 mean %s (max %s) over %d samples" % (
        fmt_num(lf["sra_mean"], 0), fmt_num(lf["sra_max"], 0), fmt_num(lf["load_mean"], 1),
        fmt_num(lf["load_max"], 1), lf["n"])
    if busy:
        return ("**Load caveat:** the timings were taken while mex carried the Logan load (%s). They are "
                "labelled as such and are not the timing baseline; the rerun after the Logan backlog drains is."
                % md(figures))
    return "Load: mex was not under the Logan load during this run (%s)." % md(figures)


def environment_section(run, comp, rows, samples):
    lines = ["# Environment", ""]
    man = run.get("manifest") or {}
    instance = run.get("server_instance")
    restarted = any("restarted" in w for w in run.get("warnings", []) or [])
    kv = [
        ("run id", md(run.get("run_id"))),
        ("release", md(run.get("release"))),
        ("server_instance", md(instance) + (" (RESTARTED during the run)" if restarted else "")),
        ("index_fp", md(short(run.get("index_fp"), 16) + "..." if run.get("index_fp") else "not stated by the host")),
        ("manifest", "version %s, sha256 %s, %s scenario(s)%s" % (
            md(man.get("manifest_version")), md(short(man.get("computed_sha256") or man.get("content_sha256"), 16)),
            md(man.get("scenarios_total")), "" if man.get("sha_verified", True) else " (content_sha256 MISMATCH)")),
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
                       md(info.get("executions_recorded", "-")), md(info.get("checks_module") or "-"),
                       md(info.get("started_utc") or "-"), md(info.get("ended_utc") or "-"),
                       md(info.get("runner_host") or "-"), fmt_num(load_s[0] if load_s else None, 1),
                       fmt_num(load_e[0] if load_e else None, 1))))
    lines.append(table(["", "value"], [[k, v] for k, v in kv], align=["l", "l"]))
    warnings = run.get("warnings") or []
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


def agreement_section(comp):
    s = comp["summary"]
    da = s["direct_agreement"]
    out = ["# Agreement between the arms", "", "## direct-local vs direct-remote", ""]
    if not da.get("compared"):
        out.append("Not compared (only one direct arm ran).\n")
    else:
        out.append("%d scenario(s) compared: %d identical answers (digest), %d with the same check outcomes but "
                   "different answers, %d not comparable (different answers, no check evaluated on both), %d "
                   "disagreeing; %d check result(s) compared." % (
                       da["compared"], da["identical"], da["same_checks"], da.get("not_comparable", 0),
                       da["disagree"], da.get("checks_compared", 0)))
        out.append("")
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
                rows.append([md(area), md(arm), md(HOST_CACHE_LABEL.get(hc, hc)), a["n"], fmt_ms(a["p50"]),
                             fmt_ms(a["p90"]), fmt_ms(a["max"]), fmt_ms(sv["p50"]), fmt_ms(sv["p90"]),
                             fmt_ms(sv["max"])])
    out.append(table(["area", "arm", "host cache", "n", "wall p50", "wall p90", "wall max", "server p50",
                      "server p90", "server max"], rows))
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
           "executions only (both arms answered from a warm row cache), so the difference is not a cache effect. "
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
            rows.append([md(arm), w["n"], fmt_ms(w["p50"]), fmt_ms(w["p90"]), fmt_ms(sv.get("p50")),
                         fmt_ms(sv.get("p90"))])
    if rows:
        out.append("Cold executions (a scenario's first on the host; the arm that ran it first), reported apart "
                   "from the warm pairs:\n")
        out.append(table(["arm", "n", "wall p50", "wall p90", "server p50", "server p90"], rows))
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


def resources_section(comp):
    s = comp["summary"]
    out = ["# Resource use per area", "",
           "Samples of `docker stats` every 2 s joined to each execution's time window (start - 1 s to end + 2 s, "
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
            tn = (s_new["timing"].get(area, {}).get(arm) or {}).get("all", {})
            to = (s_old.get("timing", {}).get(area, {}).get(arm) or {}).get("all", {})
            p_new, p_old = tn.get("p50"), to.get("p50")
            delta = "-"
            if p_new is not None and p_old:
                delta = "%+.0f%%" % (100.0 * (p_new - p_old) / p_old)
            rows.append([md(area), md(arm), "%s -> %s" % (o.get("pass", "-"), n.get("pass", "-")),
                         "%s -> %s" % (o.get("fail", 0) + o.get("error", 0) if o else "-",
                                       n.get("fail", 0) + n.get("error", 0) if n else "-"),
                         "%s -> %s" % (fmt_ms(p_old), fmt_ms(p_new)), delta])
    out.append(table(["area", "arm", "pass", "fail+error", "wall p50 (ms)", "change"], rows))
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


def build_report(run_dir, compare_dir=None, max_failures=200, by_design_path=None):
    run = sc.read_json(os.path.join(run_dir, "run.json"), {}) or {}
    rows, _ = sc.read_jsonl(os.path.join(run_dir, "executions.jsonl"))
    samples, _ = sc.read_jsonl(os.path.join(run_dir, "samples.jsonl"))
    manifest = sc.read_json(os.path.join(run_dir, "manifest.json"))
    comp = load_compare(run_dir, by_design_path)
    title = "MetaGraph suite: %s" % ascii_safe(run.get("run_id") or os.path.basename(os.path.abspath(run_dir)))
    parts = ["---", "title: \"%s\"" % title.replace('"', "'"),
             "subtitle: \"/traverse, /resolve, /pattern on %s\"" % ascii_safe(
                 ((manifest or {}).get("index") or {}).get("database") or "refseq33m-experimental"),
             "date: \"%s\"" % sc.iso_utc(time.time()), "geometry: margin=1.8cm", "fontsize: 10pt", "---", ""]
    parts.append(verdict_section(run, comp, rows, samples))
    parts.append(environment_section(run, comp, rows, samples))
    parts.append(outcomes_section(comp))
    parts.append(failures_section(rows, manifest, max_failures))
    parts.append(findings_section(rows, manifest, max_failures))
    parts.append(agreement_section(comp))
    parts.append(relations_section(comp))
    parts.append(timing_section(comp))
    parts.append(overhead_section(comp))
    parts.append(resources_section(comp))
    if compare_dir:
        old = load_compare(compare_dir, by_design_path)
        parts.append(compare_runs_section(comp, old, compare_dir))
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
    args = ap.parse_args(argv)
    text = build_report(args.run_dir, args.compare, args.max_failures, args.by_design)
    md_path = os.path.join(args.run_dir, "report.md")
    with open(md_path, "w", encoding="ascii", errors="replace") as fh:
        fh.write(text)
    sys.stdout.write("wrote %s\n" % md_path)
    if not args.no_pdf:
        ok, msg = to_pdf(md_path, os.path.join(args.run_dir, "report.pdf"))
        sys.stdout.write(msg + "\n")
        if not ok:
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
