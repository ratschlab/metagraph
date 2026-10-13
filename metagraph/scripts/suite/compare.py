#!/usr/bin/env python3
"""Comparator of the suite: joins the arms of a run per scenario and writes RUN_DIR/compare.json.

    compare.py RUN_DIR [--by-design service_by_design.json]

Reads the merged run (merge_run.py; it merges first when RUN_DIR/arms exists): run.json, executions.jsonl,
samples.jsonl, manifest.json. Per scenario and arm: the outcome over the repeats, the repeats' digests (equal
or not), wall and server times. Between the direct arms: digest and per-check agreement. For the MCP arm the
outcome classes of DESIGN §1.3: agree / disagree / not_implemented_mcp / service_by_design / missing.
Overheads of §5: direct-remote minus direct-local (network + front), the MCP phases and MCP minus direct with
and without the queue; cold vs warm per arm; timing distributions and resource use (samples joined by time) per
area and arm.

Service-side by design: async's table (service_by_design.json next to this script, from the search-service at
commit 22d034ba; --by-design overrides) holds cases {id, direct_expectation, service_answer, agreement}. A
scenario names its cases in `mcp.by_design` (an id or a list). agreement `not_implemented_mcp` makes the
scenario "not implemented yet via MCP"; `service_only` (no host counterpart) makes it service-side by design;
`by_design` makes it service-side by design when the pair holds: the direct answer matches
direct_expectation (http, code, and the text with digits ignored, each where stated) and the MCP answer
matches service_answer (refused_at_submit / relayed: error.code, text_contains, field where the MCP execution
states them; answered / clamped: no error). A pair that does not hold is a disagreement, naming the side.
"""

import argparse
import bisect
import hashlib
import json
import os
import re
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import merge_run  # noqa: E402
import suite_common as sc  # noqa: E402

OUTCOME_KEYS = ("pass", "fail", "error", "not_implemented_mcp", "service_by_design", "mixed")
SAMPLE_BEFORE_S = 1.0
SAMPLE_AFTER_S = 2.0
SAMPLE_NEAREST_S = 5.0


# ---------------------------------------------------------------------------------------------------------------
# service-side by design

DEFAULT_BY_DESIGN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "service_by_design.json")


def load_by_design(path):
    """(cases by id, info) from async's table; ({}, None) without one."""
    if not path or not os.path.exists(path):
        return {}, None
    with open(path, "rb") as fh:
        data_bytes = fh.read()
    data = json.loads(data_bytes.decode("utf-8"))
    cases = dict((c["id"], c) for c in data.get("cases", []) if isinstance(c, dict) and c.get("id"))
    info = {"path": os.path.abspath(path), "service_commit": data.get("service_commit"),
            "version": data.get("version"), "cases": len(cases),
            "sha256": hashlib.sha256(data_bytes).hexdigest()}
    return cases, info


def _digits_free(text):
    return re.sub(r"[0-9]+", "#", text or "")


def direct_errors(line, run_dir):
    """(codes, texts) of a direct answer: the top-level error's and each pattern slot's (from the raw answer)."""
    codes, texts = set(), set()
    err = line.get("error") or {}
    if err.get("code"):
        codes.add(err["code"])
    if err.get("text"):
        texts.add(err["text"])
    raw = line.get("raw")
    if raw and run_dir:
        path = os.path.join(run_dir, raw)
        if os.path.exists(path):
            try:
                body = (sc.read_gz_json(path).get("answer") or {}).get("body")
            except (OSError, ValueError):
                body = None
            if isinstance(body, dict):
                for p in body.get("patterns") or []:
                    if isinstance(p, dict) and isinstance(p.get("error"), dict):
                        if p["error"].get("code"):
                            codes.add(p["error"]["code"])
                        if p["error"].get("message"):
                            texts.add(p["error"]["message"])
    return codes, texts


def mcp_error(line, run_dir):
    """The MCP execution's error {code, text, field} from its line, else from its raw answer; None when it
    states none."""
    err = line.get("error")
    if isinstance(err, dict) and (err.get("code") or err.get("text")):
        return err
    raw = line.get("raw")
    if raw and run_dir and os.path.exists(os.path.join(run_dir, raw)):
        try:
            obj = sc.read_gz_json(os.path.join(run_dir, raw))
        except (OSError, ValueError):
            return None
        if isinstance(obj, dict):
            for cand in ((obj.get("answer") or {}).get("error") if isinstance(obj.get("answer"), dict) else None,
                         obj.get("error")):
                if isinstance(cand, dict) and (cand.get("code") or cand.get("text")):
                    return cand
    return None


def match_direct(expectation, line, run_dir):
    """(True / False / None, detail): does the direct answer match a case's direct_expectation (None: no
    direct execution, or nothing stated to compare)."""
    if line is None:
        return None, "no direct execution"
    problems, compared = [], 0
    if expectation.get("http") is not None:
        compared += 1
        if line.get("http") != expectation["http"]:
            problems.append("http %s, the case expects %s" % (line.get("http"), expectation["http"]))
    codes, texts = direct_errors(line, run_dir)
    if expectation.get("code") is not None:
        compared += 1
        if expectation["code"] not in codes:
            problems.append("code %s, the case expects %s" % (sorted(codes) or None, expectation["code"]))
    if expectation.get("text") is not None:
        compared += 1
        want = _digits_free(expectation["text"])
        if not any(_digits_free(t) == want for t in texts):
            problems.append("text %r, the case expects %r" % (sorted(texts)[:1], expectation["text"]))
    if problems:
        return False, "; ".join(problems)
    return (True if compared else None), ("matches" if compared else "nothing stated to compare")


def match_service(service_answer, lines, run_dir):
    """(True / False / None, detail): do the MCP executions match a case's service_answer (None: the MCP
    execution states nothing to compare, and its runner marked it service_by_design)."""
    outcome = service_answer.get("outcome")
    for line in lines:
        err = mcp_error(line, run_dir)
        if outcome in ("refused_at_submit", "relayed"):
            if err is None:
                if line.get("outcome") == "service_by_design":
                    continue
                return False, "repeat %s: no error, the case expects %s %s" % (
                    line.get("repeat"), outcome, service_answer.get("code"))
            if service_answer.get("code") and err.get("code") != service_answer["code"]:
                return False, "repeat %s: code %s, the case expects %s" % (
                    line.get("repeat"), err.get("code"), service_answer["code"])
            tc = service_answer.get("text_contains")
            if tc and err.get("text") and tc not in err["text"]:
                return False, "repeat %s: text %r lacks %r" % (line.get("repeat"), err["text"][:200], tc)
            if service_answer.get("field") and err.get("field") and err["field"] != service_answer["field"]:
                return False, "repeat %s: field %s, the case expects %s" % (
                    line.get("repeat"), err["field"], service_answer["field"])
        elif outcome in ("answered", "clamped"):
            if err is not None and err.get("code"):
                return False, "repeat %s: error %s, the case expects the service to answer (%s)" % (
                    line.get("repeat"), err.get("code"), outcome)
    if not lines:
        return None, "no MCP execution"
    verified = any(mcp_error(x, run_dir) is not None for x in lines) or outcome in ("answered", "clamped")
    return (True if verified else None), ("matches" if verified else "marked service_by_design by the MCP runner")


def by_design_class(case, direct_line, mcp_lines, run_dir):
    """The MCP outcome class a table case gives a scenario (DESIGN §1.3), with the pair's verification."""
    agreement = case.get("agreement")
    base = {"case": case.get("id"), "case_agreement": agreement,
            "service_outcome": (case.get("service_answer") or {}).get("outcome")}
    if agreement == "not_implemented_mcp":
        return dict(base, **{"class": "not_implemented_mcp"})
    if agreement == "service_only":
        return dict(base, **{"class": "service_by_design"})
    d_ok, d_detail = match_direct(case.get("direct_expectation") or {}, direct_line, run_dir)
    m_ok, m_detail = match_service(case.get("service_answer") or {}, mcp_lines, run_dir)
    base.update({"direct_verified": d_ok, "direct_detail": d_detail, "mcp_verified": m_ok, "mcp_detail": m_detail})
    if not mcp_lines:
        return dict(base, **{"class": "missing"})
    if d_ok is False or m_ok is False:
        base["side"] = "direct" if d_ok is False else "mcp"
        return dict(base, **{"class": "disagree"})
    return dict(base, **{"class": "service_by_design"})


# ---------------------------------------------------------------------------------------------------------------
# per-scenario joins

def aggregate_outcome(outcomes):
    """The outcome of a scenario on one arm over its repeats: the common one, else error > fail > mixed."""
    vals = [o for o in outcomes if o]
    if not vals:
        return None
    if len(set(vals)) == 1:
        return vals[0]
    if "error" in vals:
        return "error"
    if "fail" in vals:
        return "fail"
    return "mixed"


def check_outcomes(checks):
    """{check name: [outcome, ...]} in the order listed (a name may appear more than once)."""
    out = {}
    for c in checks or []:
        out.setdefault(c.get("check"), []).append(c.get("outcome"))
    return out


def compare_checks(a, b):
    """(differing, unmatched): check names whose outcomes differ between two executions, ignoring not_exposed
    on either side; names only one side evaluated are `unmatched` (not a disagreement)."""
    ca, cb = check_outcomes(a), check_outcomes(b)
    differing, unmatched = [], []
    for name in sorted(set(ca) | set(cb), key=lambda n: str(n)):
        if name not in ca or name not in cb:
            unmatched.append(name)
            continue
        for oa, ob in zip(ca[name], cb[name]):
            if "not_exposed" in (oa, ob):
                continue
            if oa != ob:
                differing.append({"check": name, "a": oa, "b": ob})
                break
    return differing, unmatched


def group_lines(rows):
    """{scenario: {arm: {repeat: line}}} and the scenario order of first appearance."""
    grouped, order = {}, []
    for r in rows:
        sid = r.get("scenario")
        if sid not in grouped:
            grouped[sid] = {}
            order.append(sid)
        grouped[sid].setdefault(r.get("arm"), {})[r.get("repeat", 0)] = r
    return grouped, order


def _num(v):
    return v if isinstance(v, (int, float)) and not isinstance(v, bool) else None


def _diff(a, b):
    a, b = _num(a), _num(b)
    return None if a is None or b is None else a - b


def _server_elapsed(line):
    return _num((line.get("server") or {}).get("elapsed_ms"))


def arm_view(lines):
    reps = [lines[k] for k in sorted(lines)]
    digests = [r.get("answer_digest") for r in reps]
    return {
        "outcome": aggregate_outcome([r.get("outcome") for r in reps]),
        "outcomes": [r.get("outcome") for r in reps],
        "http": [r.get("http") for r in reps],
        "digests_equal": (len(set(digests)) == 1) if len(digests) > 1 else None,
        "wall_ms": [r.get("wall_ms") for r in reps],
        "server_elapsed_ms": [_server_elapsed(r) for r in reps],
        "failed_checks": sorted(set(c.get("check") for r in reps for c in r.get("checks") or []
                                    if c.get("outcome") == "fail"), key=str),
    }


def direct_agreement(local, remote):
    """Agreement of the two direct arms on repeat 0 (and every repeat both ran)."""
    common = sorted(set(local) & set(remote))
    if not common:
        return None
    digest_equal = all(local[r].get("answer_digest") == remote[r].get("answer_digest") for r in common)
    differing = []
    for r in common:
        d, _ = compare_checks(local[r].get("checks"), remote[r].get("checks"))
        for x in d:
            x["repeat"] = r
        differing.extend(d)
    outcome_equal = all(local[r].get("outcome") == remote[r].get("outcome") for r in common)
    if digest_equal:
        cls = "identical"
    elif not differing and outcome_equal:
        cls = "same_checks"
    else:
        cls = "disagree"
    return {"class": cls, "digest_equal": digest_equal, "outcome_equal": outcome_equal,
            "differing_checks": differing, "repeats": common}


def case_ids(scenario):
    """The table cases a scenario names in mcp.by_design (an id or a list of ids)."""
    ids = ((scenario or {}).get("mcp") or {}).get("by_design")
    if isinstance(ids, str):
        return [ids]
    if isinstance(ids, list):
        return [i for i in ids if isinstance(i, str)]
    return []


def mcp_agreement(mcp, direct_lines, scenario, cases, run_dir):
    """The MCP arm's outcome class against the reference direct arm (DESIGN §1.3).

    Order: the manifest's `mcp.not_implemented`; the scenario's table cases (the first whose pair holds, else
    the first); the MCP runner's own marks; else the per-check agreement with direct-local (or direct-remote).
    """
    mcp = mcp or {}
    reps = [mcp[k] for k in sorted(mcp)]
    ref_arm, ref = None, None
    for arm in ("direct-local", "direct-remote"):
        if direct_lines.get(arm):
            ref_arm, ref = arm, direct_lines[arm]
            break
    ref0 = ref[min(ref)] if ref else None
    if ((scenario or {}).get("mcp") or {}).get("not_implemented"):
        return {"class": "not_implemented_mcp", "by": "manifest", "field": scenario["mcp"]["not_implemented"]}
    ids = case_ids(scenario)
    if ids:
        results = []
        for cid in ids:
            case = cases.get(cid)
            if case is None:
                results.append({"class": "disagree", "case": cid, "side": "table",
                                "detail": "case %s is not in the by-design table" % cid})
                continue
            results.append(by_design_class(case, ref0, reps, run_dir))
        for r in results:
            if r["class"] != "disagree":
                return r
        return results[0]
    if not reps:
        return {"class": "missing"}
    outcomes = set(r.get("outcome") for r in reps)
    if "not_implemented_mcp" in outcomes:
        return {"class": "not_implemented_mcp", "by": "mcp runner"}
    if "service_by_design" in outcomes:
        return {"class": "service_by_design", "by": "mcp runner"}
    if ref is None:
        return {"class": "no_direct_reference", "mcp_outcome": aggregate_outcome(list(outcomes))}
    differing = []
    for r in reps:
        d, _ = compare_checks(r.get("checks"), ref.get(r.get("repeat"), ref0).get("checks"))
        for x in d:
            x["repeat"] = r.get("repeat")
        differing.extend(d)
    mcp_out = aggregate_outcome([r.get("outcome") for r in reps])
    ref_out = aggregate_outcome([x.get("outcome") for x in ref.values()])
    # Agreement is on the checks both arms evaluated (DESIGN §1.3: a field the MCP does not expose is recorded,
    # not compared), so a direct failure of a check the MCP cannot see is no disagreement; an execution that
    # errored on one side only is.
    one_sided_error = (mcp_out == "error") != (ref_out == "error")
    cls = "agree" if not differing and not one_sided_error else "disagree"
    return {"class": cls, "reference": ref_arm, "differing_checks": differing, "mcp_outcome": mcp_out,
            "reference_outcome": ref_out}


def front_overhead(local, remote):
    """direct-remote minus direct-local wall time (the laptop's network + the front proxy), cold and warm."""
    out = {}
    if 0 in local and 0 in remote:
        out["cold_ms"] = _diff(remote[0].get("wall_ms"), local[0].get("wall_ms"))
        out["server_cold_ms"] = _diff(_server_elapsed(remote[0]), _server_elapsed(local[0]))
    if 1 in local and 1 in remote:
        out["warm_ms"] = _diff(remote[1].get("wall_ms"), local[1].get("wall_ms"))
    return out or None


def mcp_overhead(mcp, local, remote):
    """MCP minus direct (DESIGN §5): MCP warm (repeat 1, else 0) against each direct arm's line of the same
    repeat (else its repeat 0), with and without the queue wait, and the service's backend time against both."""
    # Only executions that reached the host measure the MCP path; a refusal at submit or a request the MCP
    # cannot express says nothing about it.
    mcp = dict((k, v) for k, v in (mcp or {}).items() if v.get("outcome") in ("pass", "fail"))
    if not mcp:
        return None
    m = mcp.get(1) or mcp.get(min(mcp))
    phases = m.get("mcp_phases") or {}
    queue = _num(phases.get("queue_ms"))
    backend = _num(phases.get("backend_ms"))
    wall = _num(m.get("wall_ms"))
    out = {"repeat": m.get("repeat"), "cold_wall_ms": _num((mcp.get(0) or {}).get("wall_ms")), "phases": phases}
    for name, lines in (("remote", remote), ("local", local)):
        if not lines:
            continue
        d = lines.get(m.get("repeat")) or lines.get(min(lines))
        dwall = _num(d.get("wall_ms"))
        out["minus_%s_ms" % name] = _diff(wall, dwall)
        out["minus_%s_no_queue_ms" % name] = _diff(wall - queue if wall is not None and queue is not None else None,
                                                   dwall)
        out["backend_minus_%s_ms" % name] = _diff(backend, dwall)
    return out


def cold_warm(lines):
    if 0 not in lines or 1 not in lines:
        return None
    cold, warm = _num(lines[0].get("wall_ms")), _num(lines[1].get("wall_ms"))
    if cold is None or warm is None:
        return None
    return {"cold_ms": cold, "warm_ms": warm, "delta_ms": warm - cold, "ratio": (warm / cold) if cold > 0 else None}


# ---------------------------------------------------------------------------------------------------------------
# resources

class SampleIndex(object):
    """Samples sorted by time, for the join of executions to the resource use while they ran."""

    def __init__(self, samples):
        pairs = []
        for s in samples:
            t = sc.parse_utc(s.get("t_utc"))
            if t is not None:
                pairs.append((t, s))
        pairs.sort(key=lambda p: p[0])
        self.times = [p[0] for p in pairs]
        self.samples = [p[1] for p in pairs]

    def window(self, start, end):
        """The samples from start - 1 s to end + 2 s (docker's CPU figure averages over its own window), or
        the nearest one within 5 s when none falls inside."""
        if not self.times or start is None:
            return []
        end = start if end is None else end
        lo = bisect.bisect_left(self.times, start - SAMPLE_BEFORE_S)
        hi = bisect.bisect_right(self.times, end + SAMPLE_AFTER_S)
        if hi > lo:
            return self.samples[lo:hi]
        best, best_d = None, None
        for i in (lo - 1, lo):
            if 0 <= i < len(self.times):
                d = min(abs(self.times[i] - start), abs(self.times[i] - end))
                if best_d is None or d < best_d:
                    best, best_d = self.samples[i], d
        return [best] if best is not None and best_d <= SAMPLE_NEAREST_S else []


def line_window(line):
    start = sc.parse_utc(line.get("started_utc"))
    end = sc.parse_utc(line.get("ended_utc"))
    if end is None and start is not None and _num(line.get("wall_ms")) is not None:
        end = start + line["wall_ms"] / 1000.0
    return start, end


def resources(lines, index):
    """Peak container CPU / memory, mean sra_all CPU and host load over the executions' windows, and the peak
    server memory_bytes."""
    cpu, mem, sra, load, mem_server = [], [], [], [], []
    joined = 0
    for line in lines:
        start, end = line_window(line)
        win = index.window(start, end)
        if win:
            joined += 1
        for s in win:
            cpu.append(s.get("cpu_pct"))
            mem.append(s.get("mem_bytes"))
            sra.append(s.get("sra_all_cpu_pct"))
            load.append(s.get("host_load1"))
        work = (line.get("server") or {}).get("work") or {}
        mem_server.append(work.get("memory_bytes"))
    nums = lambda vals: [v for v in vals if _num(v) is not None]  # noqa: E731
    return {"executions": len(lines), "executions_with_samples": joined,
            "peak_cpu_pct": max(nums(cpu)) if nums(cpu) else None,
            "peak_mem_bytes": max(nums(mem)) if nums(mem) else None,
            "mean_sra_all_cpu_pct": (sum(nums(sra)) / len(nums(sra))) if nums(sra) else None,
            "mean_host_load1": (sum(nums(load)) / len(nums(load))) if nums(load) else None,
            "peak_server_memory_bytes": max(nums(mem_server)) if nums(mem_server) else None}


# ---------------------------------------------------------------------------------------------------------------
# the comparison

def compare(run_dir, by_design=None):
    """compare.json's content for a merged run directory; by_design is load_by_design's (cases, info)."""
    run = sc.read_json(os.path.join(run_dir, "run.json"), {}) or {}
    rows, bad = sc.read_jsonl(os.path.join(run_dir, "executions.jsonl"))
    samples, _ = sc.read_jsonl(os.path.join(run_dir, "samples.jsonl"))
    manifest = sc.read_json(os.path.join(run_dir, "manifest.json"))
    cases, table_info = by_design if by_design else ({}, None)
    titles, areas, manifest_order, scen_by_id = {}, {}, [], {}
    if isinstance(manifest, dict):
        for scen in manifest.get("scenarios", []):
            titles[scen.get("id")] = scen.get("title")
            areas[scen.get("id")] = scen.get("area")
            scen_by_id[scen.get("id")] = scen
            manifest_order.append(scen.get("id"))
    grouped, first_seen = group_lines(rows)
    pos = dict((sid, i) for i, sid in enumerate(manifest_order))
    seen_pos = dict((sid, i) for i, sid in enumerate(first_seen))
    order = sorted(first_seen, key=lambda sid: (pos.get(sid, len(pos)), seen_pos[sid]))
    arms_present = sorted(set(r.get("arm") for r in rows), key=lambda a: (merge_run.ARM_ORDER.get(a, 9), a))

    scenarios = {}
    per_area = {}
    timing_vals = {}
    front_vals = {"cold_ms": [], "warm_ms": [], "server_cold_ms": []}
    front_area = {}
    mcp_vals = {}
    mcp_cold_walls = []
    mcp_phase_vals = {}
    cw_vals = {}
    direct_summary = {"compared": 0, "identical": 0, "same_checks": 0, "disagree": 0, "disagreeing": []}
    mcp_summary = {"compared": 0}
    lines_by_area_arm = {}

    for sid in order:
        by_arm = grouped[sid]
        any_line = next(iter(next(iter(by_arm.values())).values()))
        area = areas.get(sid) or any_line.get("area") or sid.split("-", 1)[0]
        entry = {"area": area, "title": titles.get(sid), "arms": {}}
        for arm, lines in by_arm.items():
            view = arm_view(lines)
            entry["arms"][arm] = view
            counts = per_area.setdefault(area, {}).setdefault(arm, dict((k, 0) for k in OUTCOME_KEYS + ("total",)))
            counts[view["outcome"] if view["outcome"] in counts else "mixed"] += 1
            counts["total"] += 1
            tv = timing_vals.setdefault(area, {}).setdefault(arm, {"cold": [], "warm": [], "all": [], "server": []})
            for rep, line in lines.items():
                tv["cold" if line.get("cache", "cold" if rep == 0 else "warm") == "cold" else "warm"].append(
                    line.get("wall_ms"))
                tv["all"].append(line.get("wall_ms"))
                tv["server"].append(_server_elapsed(line))
                lines_by_area_arm.setdefault(area, {}).setdefault(arm, []).append(line)
            cw = cold_warm(lines)
            if cw:
                entry["arms"][arm]["cold_warm"] = cw
                cw_vals.setdefault(arm, {"delta_ms": [], "ratio": []})
                cw_vals[arm]["delta_ms"].append(cw["delta_ms"])
                cw_vals[arm]["ratio"].append(cw["ratio"])
        local, remote, mcp = by_arm.get("direct-local"), by_arm.get("direct-remote"), by_arm.get("mcp")
        if local and remote:
            agree = direct_agreement(local, remote)
            entry["direct_agreement"] = agree
            direct_summary["compared"] += 1
            direct_summary[agree["class"]] += 1
            if agree["class"] == "disagree":
                direct_summary["disagreeing"].append(sid)
            fo = front_overhead(local, remote)
            if fo:
                entry["front_overhead"] = fo
                for k in front_vals:
                    if fo.get(k) is not None:
                        front_vals[k].append(fo[k])
                if fo.get("cold_ms") is not None:
                    front_area.setdefault(area, []).append(fo["cold_ms"])
        if "mcp" in arms_present:
            ma = mcp_agreement(mcp, {"direct-local": local, "direct-remote": remote}, scen_by_id.get(sid), cases,
                               run_dir)
            entry["mcp_agreement"] = ma
            mcp_summary["compared"] += 1
            mcp_summary[ma["class"]] = mcp_summary.get(ma["class"], 0) + 1
            if ma["class"] == "disagree":
                mcp_summary.setdefault("disagreeing", []).append(sid)
            mo = mcp_overhead(mcp, local or {}, remote or {})
            if mo:
                entry["mcp_overhead"] = mo
                for k, v in mo.items():
                    if (k.startswith("minus_") or k.startswith("backend_minus_")) and _num(v) is not None:
                        mcp_vals.setdefault(k, []).append(v)
                if _num(mo.get("cold_wall_ms")) is not None:
                    mcp_cold_walls.append(mo["cold_wall_ms"])
                for k, v in (mo.get("phases") or {}).items():
                    if _num(v) is not None:
                        mcp_phase_vals.setdefault(k, []).append(v)
        scenarios[sid] = entry

    index = SampleIndex(samples)
    resources_out = {}
    for area, by_arm in lines_by_area_arm.items():
        for arm, lines in by_arm.items():
            resources_out.setdefault(area, {})[arm] = resources(lines, index)
    totals = {}
    for arm in arms_present:
        arm_lines = [r for r in rows if r.get("arm") == arm]
        totals[arm] = {"executions": len(arm_lines),
                       "outcomes": dict((k, sum(1 for r in arm_lines if r.get("outcome") == k)) for k in OUTCOME_KEYS),
                       "wall_ms": sc.summarise([r.get("wall_ms") for r in arm_lines]),
                       "resources": resources(arm_lines, index)}

    timing = {}
    for area, by_arm in timing_vals.items():
        for arm, tv in by_arm.items():
            timing.setdefault(area, {})[arm] = dict((k, sc.summarise(v)) for k, v in tv.items())

    return {
        "run_id": run.get("run_id"), "generated_utc": sc.iso_utc(time.time()),
        "release": run.get("release"), "server_instance": run.get("server_instance"),
        "manifest_sha256": (run.get("manifest") or {}).get("computed_sha256")
        or (run.get("manifest") or {}).get("content_sha256"),
        "arms": arms_present, "unparsable_lines": bad, "by_design_table": table_info,
        "order": order, "scenarios": scenarios,
        "summary": {
            "per_area": per_area, "totals": totals,
            "direct_agreement": direct_summary, "mcp_agreement": mcp_summary,
            "timing": timing,
            "front_overhead": {"all": dict((k, sc.summarise(v)) for k, v in front_vals.items()),
                               "per_area": dict((a, sc.summarise(v)) for a, v in front_area.items())},
            "mcp_overhead": dict((k, sc.summarise(v)) for k, v in mcp_vals.items()),
            "mcp_cold_wall_ms": sc.summarise(mcp_cold_walls),
            "mcp_phases": dict((k, sc.summarise(v)) for k, v in mcp_phase_vals.items()),
            "cold_warm": dict((arm, {"n": len(v["delta_ms"]), "delta_ms": sc.summarise(v["delta_ms"]),
                                     "ratio": sc.summarise(v["ratio"])}) for arm, v in cw_vals.items()),
            "resources": resources_out,
        },
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description="Join a run's arms per scenario into compare.json.")
    ap.add_argument("run_dir")
    ap.add_argument("--by-design", default=DEFAULT_BY_DESIGN,
                    help="async's service-by-design table (default: service_by_design.json here; '' for none)")
    ap.add_argument("--no-merge", action="store_true", help="do not merge RUN_DIR/arms first")
    args = ap.parse_args(argv)
    if not args.no_merge and os.path.isdir(os.path.join(args.run_dir, "arms")):
        merge_run.merge(args.run_dir)
    result = compare(args.run_dir, load_by_design(args.by_design))
    sc.write_json_atomic(os.path.join(args.run_dir, "compare.json"), result)
    s = result["summary"]
    sys.stdout.write("compare %s: %d scenarios, arms %s; direct agreement %s; mcp %s\n" % (
        args.run_dir, len(result["order"]), ", ".join(result["arms"]),
        dict((k, v) for k, v in s["direct_agreement"].items() if k != "disagreeing"),
        dict((k, v) for k, v in s["mcp_agreement"].items() if k != "disagreeing")))
    return 0


if __name__ == "__main__":
    sys.exit(main())
