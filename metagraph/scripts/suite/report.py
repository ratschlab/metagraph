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
    states = dict((arm, (info or {}).get("state")) for arm, info in (run.get("arms") or {}).items())
    incomplete = [arm for arm, st in states.items() if st not in (None, "done")]
    dis_direct = s["direct_agreement"].get("disagree", 0)
    dis_mcp = s["mcp_agreement"].get("disagree", 0)
    if incomplete:
        verdict = "INCOMPLETE (arm(s) not finished: %s)" % ", ".join(sorted(incomplete))
    elif failing or dis_direct or dis_mcp:
        verdict = "FAIL (%d failed or errored execution(s), %d direct and %d MCP disagreement(s))" % (
            failing, dis_direct, dis_mcp)
    else:
        verdict = "PASS (no failed execution, no disagreement between the arms)"
    lines.append("**%s**" % md(verdict))
    lines.append("")
    for arm in comp["arms"]:
        t = totals.get(arm, {})
        o = t.get("outcomes", {})
        lines.append("- %s: %d execution(s): %d pass, %d fail, %d error, %d not implemented yet via MCP, "
                     "%d service-side by design." % (md(arm), t.get("executions", 0), o.get("pass", 0),
                                                      o.get("fail", 0), o.get("error", 0),
                                                      o.get("not_implemented_mcp", 0), o.get("service_by_design", 0)))
    modules = sorted(set((info or {}).get("checks_module") for info in (run.get("arms") or {}).values()
                         if (info or {}).get("checks_module")))
    if "stub" in modules:
        lines.append("")
        lines.append("**Checks: the stub ran (HTTP status only); every other expected check reads "
                     "not_exposed. This run is not a verified run.**")
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
        ("scripts commit", md(short(run.get("scripts_commit"), 12))),
        ("generated", md(sc.iso_utc(time.time()))),
    ]
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
            rows.append([md(area), md(arm), c["pass"], c["fail"], c["error"], c["not_implemented_mcp"],
                         c["service_by_design"], c["mixed"], c["total"]])
    out = ["# Outcomes per area", "",
           "Scenario outcomes over their repeats (any error -> error, any failure -> fail).", "",
           table(["area", "arm", "pass", "fail", "error", "not impl. MCP", "by design", "mixed", "total"], rows)]
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
        expects = {}
        for e in s.get("expect") or []:
            expects.setdefault(e.get("check"), []).append(e.get("args"))
        for c in r.get("checks") or []:
            if c.get("outcome") != "fail":
                continue
            exp = expects.get(c.get("check"))
            out.append("- check %s: observed %s" % (md(c.get("check")), md(c.get("detail"), 500)))
            if exp:
                out.append("    - expected: %s" % md(json.dumps(exp[0], sort_keys=True), 400))
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
                   "different answers, %d disagreeing." % (da["compared"], da["identical"], da["same_checks"],
                                                          da["disagree"]))
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
    classes = [k for k in ("agree", "disagree", "not_implemented_mcp", "service_by_design", "missing",
                           "no_direct_reference") if ma.get(k)]
    out.append("%d scenario(s): %s." % (ma["compared"], ", ".join("%d %s" % (ma[k], md(k)) for k in classes)))
    out.append("")
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
    nim = [sid for sid in comp["order"]
           if (comp["scenarios"][sid].get("mcp_agreement") or {}).get("class") == "not_implemented_mcp"]
    if nim:
        out.append("Not implemented yet via MCP: %s.\n" % ", ".join(md(x) for x in nim[:200]))
    return "\n".join(out)


def timing_section(comp):
    s = comp["summary"]
    out = ["# Timing", "", "Client wall time (ms; request sent to last byte) and the server's own elapsed time "
           "(ms; `timing.elapsed_ms`, or the sum over patterns).", ""]
    rows = []
    for area in sorted(s["timing"], key=area_sort_key):
        for arm in comp["arms"]:
            t = s["timing"][area].get(arm)
            if not t:
                continue
            a, sv = t["all"], t["server"]
            rows.append([md(area), md(arm), a["n"], fmt_ms(a["p50"]), fmt_ms(a["p90"]), fmt_ms(a["max"]),
                         fmt_ms(sv["p50"]), fmt_ms(sv["p90"]), fmt_ms(sv["max"])])
    out.append(table(["area", "arm", "n", "wall p50", "wall p90", "wall max", "server p50", "server p90",
                      "server max"], rows))
    out.append("## Cold vs warm")
    out.append("")
    rows = []
    for area in sorted(s["timing"], key=area_sort_key):
        for arm in comp["arms"]:
            t = s["timing"][area].get(arm)
            if not t or not t["warm"]["n"]:
                continue
            rows.append([md(area), md(arm), t["cold"]["n"], fmt_ms(t["cold"]["p50"]), fmt_ms(t["cold"]["p90"]),
                         t["warm"]["n"], fmt_ms(t["warm"]["p50"]), fmt_ms(t["warm"]["p90"])])
    out.append(table(["area", "arm", "cold n", "cold p50", "cold p90", "warm n", "warm p50", "warm p90"], rows))
    rows = []
    for arm, cw in sorted(s["cold_warm"].items()):
        rows.append([md(arm), cw["n"], fmt_ms(cw["delta_ms"]["p50"]), fmt_ms(cw["delta_ms"]["p90"]),
                     fmt_num(cw["ratio"]["p50"], 2), fmt_num(cw["ratio"]["p90"], 2)])
    if rows:
        out.append("Paired repeats of the same scenario (warm minus cold; warm / cold):\n")
        out.append(table(["arm", "pairs", "delta p50", "delta p90", "ratio p50", "ratio p90"], rows))
    return "\n".join(out)


def overhead_section(comp):
    s = comp["summary"]
    out = ["# Overheads", "", "## Front: direct-remote minus direct-local", "",
           "The laptop's network and the front proxy on dnaloc-new, for the same request. The arms run one after "
           "the other, so the later arm finds the host's row cache warm: the raw wall difference mixes cache and "
           "network (and can be negative). The client-side rows compare wall minus the server's own elapsed time "
           "instead, which leaves the server's part out; they need the answer's timing, so error answers and "
           "the capabilities documents only count in the raw rows.", ""]
    fo = s["front_overhead"]
    rows = []
    labels = {"cold_net_ms": "client side (wall - server), cold", "warm_net_ms": "client side (wall - server), warm",
              "cold_ms": "raw wall, cold (repeat 0)", "warm_ms": "raw wall, warm (repeat 1)",
              "server_cold_ms": "server elapsed, cold"}
    for k in ("cold_net_ms", "warm_net_ms", "cold_ms", "warm_ms", "server_cold_ms"):
        v = fo["all"].get(k) or {}
        if v.get("n"):
            rows.append([md(labels[k]), v["n"], fmt_ms(v["p50"]), fmt_ms(v["p90"]), fmt_ms(v["max"])])
    out.append(table(["difference (ms)", "n", "p50", "p90", "max"], rows))
    net = fo.get("per_area_net") or {}
    rows = []
    for a, v in sorted(fo["per_area"].items(), key=lambda kv: area_sort_key(kv[0])):
        if not v["n"]:
            continue
        n_ = net.get(a) or {}
        rows.append([md(a), n_.get("n", 0), fmt_ms(n_.get("p50")), fmt_ms(n_.get("p90")), v["n"],
                     fmt_ms(v["p50"]), fmt_ms(v["p90"])])
    if rows:
        out.append("Per area, cold:\n")
        out.append(table(["area", "client-side n", "p50", "p90", "raw n", "raw p50", "raw p90"], rows))
    out.append("## MCP")
    out.append("")
    if not s["mcp_phases"] and not s["mcp_overhead"]:
        out.append("No MCP phases in this run.\n")
        return "\n".join(out)
    rows = [[md(k), v["n"], fmt_ms(v["p50"]), fmt_ms(v["p90"]), fmt_ms(v["max"])]
            for k, v in sorted(s["mcp_phases"].items()) if v["n"]]
    out.append("Phases of the MCP executions (DESIGN section 5: submit, queue, backend as the service measured it, "
               "terminal lag, paging):\n")
    out.append(table(["phase (ms)", "n", "p50", "p90", "max"], rows))
    rows = [[md(k), v["n"], fmt_ms(v["p50"]), fmt_ms(v["p90"]), fmt_ms(v["max"])]
            for k, v in sorted(s["mcp_overhead"].items()) if v["n"]]
    out.append("MCP (warm where repeated) minus the direct arms, with and without the queue; the backend time "
               "against each direct arm (executions that reached the host only):\n")
    out.append(table(["difference (ms)", "n", "p50", "p90", "max"], rows))
    cw = s.get("mcp_cold_wall_ms") or {}
    if cw.get("n"):
        out.append("MCP wall time of the first (cold) execution: n %d, p50 %s, p90 %s, max %s ms.\n" % (
            cw["n"], fmt_ms(cw["p50"]), fmt_ms(cw["p90"]), fmt_ms(cw["max"])))
    return "\n".join(out)


def resources_section(comp):
    s = comp["summary"]
    out = ["# Resource use per area", "",
           "Samples of `docker stats` every 2 s joined to each execution's time window (start - 1 s to end + 2 s, "
           "else the nearest sample within 5 s). Peak server memory is the answers' `work.memory_bytes`.", ""]
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
    parts.append(agreement_section(comp))
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
