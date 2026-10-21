"""Offline tests of merge_run.py, compare.py and report.py on synthetic run directories: the merge of the arms,
outcome counts, direct agreement, the MCP outcome classes with async's by-design table, overheads, cold vs
warm, the resource join, the report's sections and escaping, --compare, and (when pandoc and xelatex are
installed) the PDF."""

import hashlib
import json
import os
import shutil
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.dirname(HERE)
sys.path.insert(0, os.path.dirname(HERE))

import compare as cmp  # noqa: E402
import merge_run  # noqa: E402
import report  # noqa: E402
import suite_common as sc  # noqa: E402

T0 = 1791800000.0  # a fixed epoch for the synthetic run

CASES = {
    "case.nim": {"id": "case.nim", "agreement": "not_implemented_mcp",
                 "direct_expectation": {"http": 400, "code": "invalid_request", "text": "request: unknown field 'x'"},
                 "service_answer": {"outcome": "not_expressible", "missing_argument": "x"}},
    "case.answered": {"id": "case.answered", "agreement": "by_design",
                      "direct_expectation": {"http": 400, "code": "invalid_request",
                                             "text": "request.patterns: expected a list of 1 to 16 pattern objects "
                                                     "(got 17; the server's --pattern-max-patterns)"},
                      "service_answer": {"outcome": "answered", "code": None}},
    "case.refused": {"id": "case.refused", "agreement": "by_design",
                     "direct_expectation": {"http": 400, "code": "invalid_request", "text": None},
                     "service_answer": {"outcome": "refused_at_submit", "code": "pattern_invalid",
                                        "text_contains": "PATTERN_MAX_PATTERNS_PER_SEARCH", "field": "patterns"}},
    "case.slot": {"id": "case.slot", "agreement": "by_design",
                  "direct_expectation": {"http": 200, "code": "mismatches_above_cap", "text": None},
                  "service_answer": {"outcome": "refused_at_submit", "code": "pattern_out_of_bounds"}},
    "case.only": {"id": "case.only", "agreement": "service_only",
                  "direct_expectation": {"http": None, "code": None, "text": None},
                  "service_answer": {"outcome": "refused_at_submit", "code": "rate_limited"}},
}


def scen(sid, area, by_design=None, not_implemented=None, expect_http=200):
    mcp = {"tool": "pattern_search", "arguments": {}}
    if by_design is not None:
        mcp["by_design"] = by_design
    if not_implemented:
        mcp = {"not_implemented": not_implemented}
    return {"id": sid, "area": area, "title": "title of %s with ≤ and $x_1$ | pipe" % sid,
            "route": "POST /pattern", "request": {}, "mcp": mcp, "truth": {},
            "expect": [{"check": "status", "args": {"http": expect_http}},
                       {"check": "contains_source_context", "args": {"level": 0}}],
            "timing_class": "fast", "repeat": 1, "tags": []}


# An MCP refusal at submit as async's runner records it: no host timing, the backend phase null.
SUBMIT_ONLY = {"submit_ms": 40.0, "queue_ms": None, "backend_ms": None, "terminal_lag_ms": None, "paging_ms": None}


def line(sid, arm, repeat, t, wall, outcome="pass", http=200, digest="d", checks=None, error=None,
         memory=None, phases=None, area=None, server_ms=None, reached=True):
    """An execution line; an MCP line with reached=False (or not implemented via MCP) carries no evidence of a
    backend call (server null, mcp_phases without backend_ms)."""
    checks = checks if checks is not None else [{"check": "status", "outcome": "pass", "detail": "ok"},
                                                {"check": "contains_source_context", "outcome": "pass",
                                                 "detail": "found"}]
    work = {"steps": 10}
    if memory is not None:
        work["memory_bytes"] = memory
    server = {"elapsed_ms": server_ms if server_ms is not None else wall / 2.0, "timing": None, "work": work}
    if arm == "mcp" and outcome == "not_implemented_mcp":
        server, phases = None, None
    elif arm == "mcp" and not reached:
        server, phases = None, dict(SUBMIT_ONLY)
    return {"run_id": "R", "scenario": sid, "area": area or sid.split("-")[0], "arm": arm, "repeat": repeat,
            "cache": "cold" if repeat == 0 else "warm", "started_utc": sc.iso_utc(t),
            "ended_utc": sc.iso_utc(t + wall / 1000.0), "wall_ms": wall, "http": http, "outcome": outcome,
            "checks": checks, "server": server,
            "mcp_phases": phases, "answer_digest": digest, "raw": "raw/%s.%s.%d.json.gz" % (sid, arm, repeat),
            "job": None, "error": error}


def write_arm(run_dir, arm, rows, run_extra=None, samples=None, raw_bodies=None, jobs=None):
    d = os.path.join(run_dir, "arms", arm)
    os.makedirs(os.path.join(d, "raw"), exist_ok=True)
    for r in rows:
        sc.append_jsonl(os.path.join(d, "executions.jsonl"), r)
        body = (raw_bodies or {}).get((r["scenario"], r["repeat"]))
        sc.write_gz_json(os.path.join(d, r["raw"]), {"answer": {"http": r["http"], "body": body,
                                                                "error": r.get("error")}})
    for s in samples or []:
        sc.append_jsonl(os.path.join(d, "samples.jsonl"), s)
    for j in jobs or []:
        sc.append_jsonl(os.path.join(d, "jobs.jsonl"), j)
    run = {"run_id": "R", "arm": arm, "state": "done", "checks_module": "checks", "base": "http://x",
           "manifest": {"computed_sha256": "sha1", "manifest_version": 1, "scenarios_total": 9},
           "release": "refseq97-b7214d29", "server_instance": "inst1", "index_fp": "fp",
           "sessions": [{"server_instance_start": "inst1", "server_instance_end": "inst1"}],
           "host_load_start": [200.0, 190.0, 100.0], "host_load_end": [210.0, 195.0, 110.0],
           "capabilities": {"/capabilities": {"release": "refseq97-b7214d29"}}}
    if arm == "mcp":
        run.pop("index_fp")   # async's runner states the host's index_fp in service_capabilities.host_document only
    run.update(run_extra or {})
    sc.write_json_atomic(os.path.join(d, "run.json"), run)


def sample(t, cpu, mem, sra, load, phase=None):
    s = {"t_utc": sc.iso_utc(t), "container": "metagraph_refseq33m_experimental", "cpu_pct": cpu, "mem_bytes": mem,
         "mem_limit_bytes": 1 << 37, "net_rx": 0, "net_tx": 0, "blk_read": 0, "blk_write": 0, "host_load1": load,
         "sra_all_cpu_pct": sra}
    if phase:
        s["phase"] = phase
    return s


def build_run(run_dir):
    """A run with three arms over nine scenarios exercising every comparison path."""
    scenarios = [scen("P-001", "P"), scen("P-002", "P", expect_http=400), scen("R-001", "R"),
                 scen("P-003", "P", by_design="case.answered", expect_http=400),
                 scen("P-004", "P", by_design=["case.refused"], expect_http=400),
                 scen("M-001", "M", by_design="case.slot"), scen("P-005", "P", by_design="case.nim"),
                 scen("M-002", "M", not_implemented="max_mismatches"), scen("S-001", "S", by_design="case.only")]
    manifest = {"manifest_version": 1, "content_sha256": sc.scenarios_sha256(scenarios), "scenarios": scenarios,
                "index": {"database": "refseq33m-experimental"}}
    os.makedirs(os.path.join(run_dir, "arms", "direct-local"), exist_ok=True)
    sc.write_json_atomic(os.path.join(run_dir, "arms", "direct-local", "manifest.json"), manifest)
    fail_checks = [{"check": "status", "outcome": "pass", "detail": "ok"},
                   {"check": "contains_source_context", "outcome": "fail",
                    "detail": "source context absent: offsets [3] ⊆ []"}]
    refused = {"code": "invalid_request", "text": "request.patterns: expected a list of 1 to 16 pattern objects "
                                                  "(got 20; the server's --pattern-max-patterns)"}
    t = T0
    local = [line("P-001", "direct-local", 0, t, 10.0, memory=1000),
             line("P-001", "direct-local", 1, t + 2, 6.0, memory=1500),
             line("P-002", "direct-local", 0, t + 4, 5.0, http=400, digest="e",
                  error={"code": "invalid_request", "text": "request: unknown field 'bogus'"}),
             line("R-001", "direct-local", 0, t + 6, 20.0, outcome="fail", checks=fail_checks, digest="r1"),
             line("P-003", "direct-local", 0, t + 8, 4.0, http=400, error=refused, digest="p3"),
             line("P-004", "direct-local", 0, t + 10, 4.0, http=400, error={"code": "invalid_request", "text": "x"},
                  digest="p4"),
             line("M-001", "direct-local", 0, t + 12, 4.0, digest="m1"),
             line("P-005", "direct-local", 0, t + 14, 4.0, http=400, digest="p5",
                  error={"code": "invalid_request", "text": "request: unknown field 'x'"}),
             line("M-002", "direct-local", 0, t + 16, 30.0, digest="m2"),
             line("S-001", "direct-local", 0, t + 18, 3.0)]
    # a resumed duplicate of P-001 r0: the later line wins
    local.append(line("P-001", "direct-local", 0, t + 1, 9.0, memory=1000))
    samples_local = [sample(t + 0.5, 50.0, 10 << 30, 7000.0, 230.0), sample(t + 2.5, 80.0, 11 << 30, 7500.0, 240.0),
                     sample(t + 6.5, 120.0, 12 << 30, 7400.0, 235.0), sample(t + 30.0, 10.0, 9 << 30, 100.0, 5.0)]
    write_arm(run_dir, "direct-local", local, samples=samples_local,
              raw_bodies={("M-001", 0): {"patterns": [{"id": "a", "error": {"code": "mismatches_above_cap",
                                                                             "message": "cap 1"}}]}})
    # direct-remote within 60 s of direct-local (host-warm), the MCP arm 20 s later (warm too), except R-001,
    # whose MCP execution runs first (the host-cold one, as in paired mode)
    tr = T0 + 30
    remote = [line("P-001", "direct-remote", 0, tr, 110.0, memory=1000),
              line("P-001", "direct-remote", 1, tr + 2, 100.0),
              line("P-002", "direct-remote", 0, tr + 4, 105.0, http=400, digest="e",
                   error={"code": "invalid_request", "text": "request: unknown field 'bogus'"}),
              line("R-001", "direct-remote", 0, tr + 6, 120.0, digest="r2"),  # passes remotely: disagreement
              line("P-003", "direct-remote", 0, tr + 8, 104.0, http=400, error=refused, digest="p3x"),
              line("A-999", "mcp", 0, tr, 1.0)]  # a line of another arm in this arm's file: dropped
    write_arm(run_dir, "direct-remote", remote, run_extra={"state": "interrupted"},
              samples=[sample(tr + 1, 30.0, 10 << 30, 7000.0, 220.0, phase="direct-remote")])
    tm = T0 + 50
    ph = {"submit_ms": 50.0, "queue_ms": 200.0, "backend_ms": 300.0, "terminal_lag_ms": 100.0, "paging_ms": 30.0}
    mcp_checks_fail = [{"check": "status", "outcome": "pass", "detail": "ok"},
                       {"check": "contains_source_context", "outcome": "fail", "detail": "absent"}]
    mcp = [line("P-001", "mcp", 0, tm, 900.0, phases=ph),
           line("P-001", "mcp", 1, tm + 2, 700.0, phases=ph),
           line("P-002", "mcp", 0, tm + 4, 800.0, checks=mcp_checks_fail, outcome="fail"),
           line("R-001", "mcp", 0, T0 + 5, 800.0, checks=[{"check": "status", "outcome": "pass", "detail": ""},
                                                          {"check": "contains_source_context",
                                                           "outcome": "not_exposed", "detail": ""}]),
           line("P-003", "mcp", 0, tm + 8, 800.0),  # answered by the service: by design
           line("P-004", "mcp", 0, tm + 10, 800.0, outcome="service_by_design", reached=False,
                error={"code": "pattern_invalid", "text": "more than 64 per search (PATTERN_MAX_PATTERNS_PER_SEARCH)",
                       "field": "patterns"}),
           line("M-001", "mcp", 0, tm + 12, 800.0, error={"code": "pattern_unsupported", "text": "no"}, http=400,
                checks=[{"check": "status", "outcome": "fail", "detail": "http 400, expected 200"},
                        {"check": "contains_source_context", "outcome": "fail", "detail": "no answer"}]),
           line("P-005", "mcp", 0, tm + 14, 1.0, outcome="not_implemented_mcp"),
           line("S-001", "mcp", 0, tm + 16, 50.0, outcome="service_by_design", reached=False,
                error={"code": "rate_limited", "text": "slow down"})]
    write_arm(run_dir, "mcp", mcp, jobs=[{"search_id": "s1", "token": "t1"}])
    return manifest


class MergeTest(unittest.TestCase):
    def test_merge(self):
        with tempfile.TemporaryDirectory() as d:
            build_run(d)
            merged = merge_run.merge(d)
            rows, bad = sc.read_jsonl(os.path.join(d, "executions.jsonl"))
            self.assertEqual(bad, 0)
            self.assertEqual(len(rows), 10 + 5 + 9)
            keys = [sc.execution_key(r) for r in rows]
            self.assertEqual(len(keys), len(set(keys)))
            p1 = [r for r in rows if sc.execution_key(r) == ("P-001", "direct-local", 0)][0]
            self.assertEqual(p1["wall_ms"], 9.0)
            self.assertTrue(p1["raw"].startswith("arms/direct-local/raw/"))
            self.assertTrue(os.path.exists(os.path.join(d, p1["raw"])))
            self.assertNotIn("A-999", [r["scenario"] for r in rows])
            samples, _ = sc.read_jsonl(os.path.join(d, "samples.jsonl"))
            times = [sc.parse_utc(s["t_utc"]) for s in samples]
            self.assertEqual(times, sorted(times))
            self.assertEqual(set(s["phase"] for s in samples), {"direct-local", "direct-remote"})
            self.assertEqual(sc.read_jsonl(os.path.join(d, "jobs.jsonl"))[0], [{"search_id": "s1", "token": "t1"}])
            self.assertEqual(merged["release"], "refseq97-b7214d29")
            self.assertEqual(merged["host_load_start"][0], 200.0)
            self.assertEqual(set(merged["arms"]), {"direct-local", "direct-remote", "mcp"})
            self.assertNotIn("capabilities", merged["arms"]["direct-local"])
            self.assertTrue(any("another arm" in w for w in merged["warnings"]))
            self.assertTrue(os.path.exists(os.path.join(d, "manifest.json")))
            hc = dict((sc.execution_key(r), (r["host_cache"], r["since_prev_s"], r["host_cache_after"]))
                      for r in rows)
            self.assertEqual(hc[("P-001", "direct-local", 0)], ("cold", None, None))
            self.assertEqual(hc[("P-001", "direct-local", 1)][0], "warm")
            self.assertAlmostEqual(hc[("P-001", "direct-local", 1)][1], 0.991)
            self.assertEqual(hc[("P-001", "direct-remote", 0)][::2], ("warm", "direct-local"))
            self.assertEqual(hc[("R-001", "mcp", 0)][0], "cold")
            self.assertEqual(hc[("R-001", "direct-local", 0)][::2], ("warm", "mcp"))
            self.assertEqual(hc[("S-001", "mcp", 0)][0], "not_reached")  # refused by the service
            self.assertEqual(merged["host_cache"]["warm_window_s"], 60.0)
            # idempotent
            merge_run.merge(d)
            self.assertEqual(len(sc.read_jsonl(os.path.join(d, "executions.jsonl"))[0]), len(rows))

    def test_merge_warnings(self):
        with tempfile.TemporaryDirectory() as d:
            write_arm(d, "direct-local", [line("P-1", "direct-local", 0, T0, 1.0)],
                      run_extra={"sessions": [{"server_instance_start": "a", "server_instance_end": "b"}]})
            write_arm(d, "direct-remote", [line("P-1", "direct-remote", 0, T0, 1.0)],
                      run_extra={"manifest": {"computed_sha256": "other"}})
            merged = merge_run.merge(d)
            text = " ".join(merged["warnings"])
            self.assertIn("different manifests", text)
            self.assertIn("restarted", text)


class CompareTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = self.tmp.name
        build_run(self.dir)
        merge_run.merge(self.dir)
        self.comp = cmp.compare(self.dir, (CASES, {"cases": len(CASES)}))
        self.s = self.comp["summary"]

    def tearDown(self):
        self.tmp.cleanup()

    def lines(self, sid):
        rows, _ = sc.read_jsonl(os.path.join(self.dir, "executions.jsonl"))
        out = {}
        for r in rows:
            if r["scenario"] == sid:
                out.setdefault(r["arm"], {})[r["repeat"]] = r
        return out

    def test_order_and_counts(self):
        self.assertEqual(self.comp["order"][:3], ["P-001", "P-002", "R-001"])
        self.assertEqual(self.comp["arms"], ["direct-local", "direct-remote", "mcp"])
        r_local = self.s["per_area"]["R"]["direct-local"]
        self.assertEqual((r_local["fail"], r_local["total"]), (1, 1))
        p_mcp = self.s["per_area"]["P"]["mcp"]
        self.assertEqual(p_mcp["fail"], 1)
        self.assertEqual(p_mcp["not_implemented_mcp"], 1)
        self.assertEqual(self.s["totals"]["direct-local"]["executions"], 10)

    def test_direct_agreement(self):
        sc_ = self.comp["scenarios"]
        self.assertEqual(sc_["P-001"]["direct_agreement"]["checks_compared"], 4)
        self.assertEqual(sc_["P-001"]["direct_agreement"]["class"], "identical")
        self.assertEqual(sc_["P-002"]["direct_agreement"]["class"], "identical")
        self.assertEqual(sc_["R-001"]["direct_agreement"]["class"], "disagree")
        self.assertEqual(sc_["R-001"]["direct_agreement"]["differing_checks"][0]["check"], "contains_source_context")
        self.assertEqual(sc_["P-003"]["direct_agreement"]["class"], "same_checks")
        da = self.s["direct_agreement"]
        self.assertEqual((da["compared"], da["identical"], da["same_checks"], da["disagree"]), (4, 2, 1, 1))
        self.assertEqual(da["disagreeing"], ["R-001"])

    def test_mcp_classes(self):
        cls = dict((sid, e.get("mcp_agreement", {}).get("class")) for sid, e in self.comp["scenarios"].items())
        self.assertEqual(cls["P-001"], "agree")
        self.assertEqual(cls["P-002"], "disagree")       # checks differ
        self.assertEqual(cls["R-001"], "agree")          # direct-local fails the check, MCP does not expose it
        self.assertEqual(cls["P-003"], "agree")  # answered (B-01): a note, the checks decide
        self.assertEqual(cls["P-004"], "service_by_design")  # refused_at_submit, code/text/field match
        self.assertEqual(cls["M-001"], "disagree")       # wrong service code: no case matched, its checks fail
        self.assertEqual(cls["P-005"], "not_implemented_mcp")
        self.assertEqual(cls["M-002"], "not_implemented_mcp")  # by the manifest, no MCP line
        self.assertEqual(cls["S-001"], "service_by_design")    # service_only
        m1 = self.comp["scenarios"]["M-001"]["mcp_agreement"]
        self.assertEqual([d["check"] for d in m1["differing_checks"]], ["contains_source_context", "status"])
        # the case's pair is still judged where the refusal matches: the slot code of the raw answer
        ok = cmp.by_design_class(CASES["case.slot"], self.lines("M-001")["direct-local"][0],
                                 [dict(self.lines("M-001")["mcp"][0], error={"code": "pattern_out_of_bounds"}, raw=None)],
                                 self.dir)
        self.assertEqual((ok["class"], ok["direct_verified"]), ("service_by_design", True))
        p3 = self.comp["scenarios"]["P-003"]["mcp_agreement"]
        self.assertEqual(p3["by_design_notes"], {"0": "case.answered"})
        self.assertEqual(self.s["mcp_agreement"]["disagree"], 2)
        # review M2: checks compared per comparison, by-design pairs verified on every side the case has
        self.assertEqual(self.comp["scenarios"]["R-001"]["mcp_agreement"]["checks_compared"], 1)
        self.assertEqual(self.comp["scenarios"]["R-001"]["mcp_agreement"]["checks_not_exposed"], 1)
        s1 = self.comp["scenarios"]["S-001"]["mcp_agreement"]
        self.assertTrue(s1["mcp_verified"] and s1["verified"])
        self.assertEqual((self.s["mcp_agreement"]["by_design_verified"],
                          self.s["mcp_agreement"]["by_design_unverified"]), (2, 0))

    def test_agree_needs_a_compared_check(self):
        ne = [{"check": "status", "outcome": "not_exposed", "detail": ""},
              {"check": "contains_source_context", "outcome": "not_exposed", "detail": ""}]
        d = {"direct-local": {0: line("P-1", "direct-local", 0, T0, 1.0)}}
        res = cmp.mcp_agreement({0: line("P-1", "mcp", 0, T0, 1.0, checks=ne)}, d, {"mcp": {}}, {}, None)
        self.assertEqual((res["class"], res["checks_compared"], res["checks_not_exposed"]), ("not_comparable", 0, 2))
        res = cmp.mcp_agreement({0: line("P-1", "mcp", 0, T0, 1.0, checks=ne, outcome="error")}, d, {"mcp": {}},
                                {}, None)
        self.assertEqual(res["class"], "disagree")  # an error on one side only
        res = cmp.direct_agreement({0: line("P-1", "direct-local", 0, T0, 1.0, checks=ne, digest="a")},
                                   {0: line("P-1", "direct-remote", 0, T0, 1.0, checks=ne, digest="b")})
        self.assertEqual(res["class"], "not_comparable")

    def test_service_only_verified(self):
        d = line("S-1", "direct-local", 0, T0, 1.0)
        ok = cmp.by_design_class(CASES["case.only"], d, [line("S-1", "mcp", 0, T0, 1.0, outcome="service_by_design",
                                                              error={"code": "rate_limited"})], None)
        self.assertEqual((ok["class"], ok["verified"]), ("service_by_design", True))
        # the service answered what it should have refused: no longer "by design"
        bad = cmp.by_design_class(CASES["case.only"], d, [line("S-1", "mcp", 0, T0, 1.0)], None)
        self.assertEqual((bad["class"], bad["side"]), ("disagree", "mcp"))
        alt = {"id": "q", "agreement": "service_only",
               "service_answer": {"outcome": "refused_at_submit", "code": "capacity_unavailable | traversal_queue_full"}}
        ok = cmp.by_design_class(alt, None, [line("S-1", "mcp", 0, T0, 1.0, error={"code": "traversal_queue_full"})],
                                 None)
        self.assertTrue(ok["verified"])

    def test_by_design_direct_side_and_unknown_case(self):
        sid_line = line("P-1", "direct-local", 0, T0, 1.0, http=200)
        res = cmp.by_design_class(CASES["case.refused"], sid_line, [line("P-1", "mcp", 0, T0, 1.0,
                                                                       error={"code": "pattern_invalid"})], None)
        self.assertEqual((res["class"], res["side"]), ("disagree", "direct"))
        res = cmp.mcp_agreement({0: line("P-1", "mcp", 0, T0, 1.0)}, {"direct-local": {0: sid_line}},
                                {"mcp": {"by_design": "no.such.case"}}, CASES, None)
        self.assertEqual((res["class"], res["side"]), ("disagree", "table"))
        # a list of cases: the first whose pair holds decides
        res = cmp.mcp_agreement({0: line("P-1", "mcp", 0, T0, 1.0, error={"code": "rate_limited"})},
                                {"direct-local": {0: sid_line}},
                                {"mcp": {"by_design": ["no.such.case", "case.only"]}}, CASES, None)
        self.assertEqual(res["class"], "service_by_design")
        res = cmp.mcp_agreement({}, {"direct-local": {0: sid_line}}, {"mcp": {"by_design": "case.answered"}},
                                CASES, None)
        self.assertEqual(res["class"], "missing")
        res = cmp.mcp_agreement({}, {"direct-local": {0: sid_line}}, {"mcp": {}}, CASES, None)
        self.assertEqual(res["class"], "missing")
        res = cmp.mcp_agreement({0: line("P-1", "mcp", 0, T0, 1.0)}, {}, {"mcp": {}}, CASES, None)
        self.assertEqual(res["class"], "no_direct_reference")

    def test_relayed_text_contains(self):
        case = {"id": "c", "agreement": "by_design", "direct_expectation": {"http": 400},
                "service_answer": {"outcome": "relayed", "code": "backend_fault", "text_contains": "HTTP 400: "}}
        d = line("R-1", "direct-local", 0, T0, 1.0, http=400)
        ok = cmp.by_design_class(case, d, [line("R-1", "mcp", 0, T0, 1.0,
                                                error={"code": "backend_fault", "text": "HTTP 400: too short"})], None)
        self.assertEqual(ok["class"], "service_by_design")
        bad = cmp.by_design_class(case, d, [line("R-1", "mcp", 0, T0, 1.0,
                                                 error={"code": "backend_fault", "text": "HTTP 500: boom"})], None)
        self.assertEqual(bad["class"], "disagree")
        # B-02: each line decided as the runner's matcher decides it: a runner's mark on an answer without the
        # relayed error does not hold the pair
        marked = cmp.by_design_class(case, d, [line("R-1", "mcp", 0, T0, 1.0, outcome="service_by_design")], None)
        self.assertEqual((marked["class"], marked["side"]), ("disagree", "mcp"))
        # a missing text or field no longer passes (compare read them only where the line had them)
        bounds = {"outcome": "refused_at_submit", "code": "traversal_invalid", "text_contains": "above the preset",
                  "field": "strategy.bounds.<knob>"}
        nofield = line("T-1", "mcp", 0, T0, 1.0, http=400,
                       error={"code": "traversal_invalid", "text": "time_budget_ms above the preset's row"})
        self.assertFalse(cmp.match_service(bounds, [nofield], None)[0])
        notext = line("T-1", "mcp", 0, T0, 1.0, http=400,
                      error={"code": "traversal_invalid", "field": "strategy.bounds.time_budget_ms"})
        self.assertFalse(cmp.match_service(bounds, [notext], None)[0])
        both = line("T-1", "mcp", 0, T0, 1.0, http=400, error={"code": "traversal_invalid", "text": "x above the preset",
                                                               "field": "strategy.bounds.time_budget_ms"})
        self.assertTrue(cmp.match_service(bounds, [both], None)[0])

    def test_answer_triggered_case(self):
        case = {"id": "pattern.deadline_ends_task", "agreement": "by_design",
                "direct_expectation": {"http": 503, "code": "deadline"},
                "service_answer": {"outcome": "relayed", "http": 502, "code": "search_failed",
                                   "text_contains": "pattern_deadline"},
                "answer_match": {"routes": ["POST /pattern"], "direct": {"http": 503, "code": "deadline"},
                                 "mcp": {"http": 502, "code": "search_failed", "text_contains": "pattern_deadline"}}}
        cases = dict(CASES, **{case["id"]: case})
        s = scen("N-1", "N")  # names no case
        st = lambda o, h: [{"check": "status", "outcome": o, "detail": "http %d" % h}]  # noqa: E731
        dead = line("N-1", "direct-local", 0, T0, 60000.0, outcome="fail", http=503, checks=st("fail", 503),
                    error={"code": "deadline", "text": "pattern: the answer could not be written within time_budget_ms"})
        relayed = line("N-1", "mcp", 0, T0, 61000.0, outcome="fail", http=502, checks=st("fail", 502),
                       error={"code": "search_failed", "text": "pattern_deadline: pattern: the answer could not be "
                                                                "written within time_budget_ms; retry with a larger"})
        res = cmp.mcp_agreement({0: relayed}, {"direct-local": {0: dead}}, s, cases, None)
        self.assertEqual((res["class"], res["by"], res["case"], res["verified"], res["repeats"]),
                         ("service_by_design", "answer", case["id"], True, [0]))
        # the host answered within the budget on the direct arm: by design, unverified (a timing difference)
        res = cmp.mcp_agreement({0: relayed}, {"direct-local": {0: line("N-1", "direct-local", 0, T0, 1.0)}}, s,
                                cases, None)
        self.assertEqual((res["class"], res["verified"], res["direct_verified"]), ("service_by_design", False, None))
        # another failed search, another route, or a table without the case: no trigger
        other = dict(relayed, error={"code": "search_failed", "text": "backend_down: no answer"})
        self.assertEqual(cmp.mcp_agreement({0: other}, {"direct-local": {0: dead}}, s, cases, None)["class"],
                         "agree")  # both read status fail: compared by the checks as before
        res = cmp.mcp_agreement({0: relayed}, {"direct-local": {0: line("N-1", "direct-local", 0, T0, 1.0)}},
                                dict(s, route="POST /resolve"), cases, None)
        self.assertEqual(res["class"], "disagree")
        res = cmp.mcp_agreement({0: relayed}, {"direct-local": {0: dead}}, s, CASES, None)
        self.assertNotEqual(res["class"], "service_by_design")
        # the scenario's other MCP executions are compared as usual
        ok1 = line("N-1", "direct-local", 1, T0, 1.0)
        bad1 = line("N-1", "mcp", 1, T0, 1.0, outcome="fail",
                    checks=[{"check": "status", "outcome": "pass", "detail": "ok"},
                            {"check": "contains_source_context", "outcome": "fail", "detail": "absent"}])
        res = cmp.mcp_agreement({0: relayed, 1: bad1}, {"direct-local": {0: dead, 1: ok1}}, s, cases, None)
        self.assertEqual((res["class"], res["answer_case"], res["answer_case_repeats"]),
                         ("disagree", case["id"], [0]))
        res = cmp.mcp_agreement({0: relayed, 1: line("N-1", "mcp", 1, T0, 1.0)}, {"direct-local": {0: dead, 1: ok1}},
                                s, cases, None)
        self.assertEqual((res["class"], res["verified"], res["other_repeats"]["class"]),
                         ("service_by_design", True, "agree"))
        # a code with alternatives, as the MCP runner reads it
        self.assertTrue(cmp.answer_side_matches({"code": "a | search_failed"}, 502, {"search_failed"}, []))
        self.assertFalse(cmp.answer_side_matches({}, 502, {"search_failed"}, []))
        self.assertIs(cmp.answer_case_of(cases, "POST /pattern", 502, relayed["error"]), case)
        self.assertIsNone(cmp.answer_case_of(cases, "POST /pattern", 503, relayed["error"]))

    def test_discovery_cap_cases(self):
        # async's fb66f1dd: the service's discovery cap 1000 (was 100), default 50 (was 20); the host has none. The
        # renamed cases state it with request predicates (the cap a reference into the live document, the table's
        # block its recorded value); the old ids stay as aliases the manifests v1-v4 name, without a request (v4's
        # R-edge-0023 sends 150) and with a text the refusal at any cap contains
        cases, _ = cmp.load_by_design(cmp.DEFAULT_BY_DESIGN)
        low, high = cases["resolve.discover_limit_above_cap"], cases["resolve.discover_limit_above_cap_host_answers"]
        ref = "backends.refseq33m-experimental.resolve.discover.max_labels"
        self.assertEqual(high["request"][1], {"path": "limit", "gt": {"service": ref}})
        self.assertEqual(low["request"][1], {"path": "limit", "lt": 1})
        self.assertEqual(sc.lookup(cases.blocks["service"], "backends.refseq33m-experimental.resolve.discover"),
                         (True, {"max_labels": 1000, "default_max_labels": 50}))
        for old, new in (("resolve.discover_limit_above_100", low), ("resolve.discover_limit_above_100_host_answers",
                                                                     high)):
            self.assertEqual(cases[old]["alias_of"], new["id"])
            self.assertNotIn("request", cases[old])
            self.assertEqual((cases[old]["direct_expectation"], cases[old]["agreement"]),
                             (new["direct_expectation"], new["agreement"]))
        # the manifests v1-v4 still load with the table; a v5 scenario must satisfy the renamed case's request
        for name in ("suite_v1", "suite_v2", "suite_v3", "suite_v4"):
            m = sc.read_json(os.path.join(SUITE, "manifest", name + ".json"))
            self.assertEqual(cmp.request_scenario_problems(cases, m["scenarios"]), [], name)

        def v5(case_id, limit):
            return {"id": "R-v5", "route": "POST /resolve", "mcp": {"tool": "traverse_resolve", "by_design": [case_id],
                                                                   "arguments": {"discover": True, "limit": limit}}}
        self.assertEqual(cmp.request_scenario_problems(cases, [v5(high["id"], 1001), v5(low["id"], 0)]), [])
        for case_id, limit in ((high["id"], 1000), (high["id"], 150), (low["id"], 1), (low["id"], 1500)):
            self.assertEqual(len(cmp.request_scenario_problems(cases, [v5(case_id, limit)])), 1, (case_id, limit))
        # the cap reads the run's live document first (a deployment before fb66f1dd states 100: a mismatch with the
        # table's 1000, judged with the live value)
        live = {"backends": {"refseq33m-experimental": {"resolve": {"discover": {"max_labels": 100}}}}}
        refs = sc.References(cases.blocks, live={"service": live})
        evidence = cmp.RequestEvidence()
        self.assertTrue(cmp.request_holds(high, v5(high["id"], 150), None, refs, evidence)[0])
        self.assertFalse(cmp.request_holds(high, v5(high["id"], 150), None, cases.refs)[0])
        self.assertEqual(evidence.sources, {high["id"]: {ref: "live"}})
        self.assertEqual([(m["ref"], m["live"], m["table"]) for m in evidence.mismatches], [(ref, 100, 1000)])
        # the refusal texts: the renamed cases the cap 1000's, the aliases any cap's (the recorded v3 / v4 runs)
        for text, renamed in (("`limit` must be an integer in 1..1000, got 1500", True),
                              ("`limit` must be an integer in 1..100, got 150", False)):
            ans = {"http": None, "error": {"code": "invalid_arguments", "text": text}}
            self.assertIs(cmp.service_answers(high["service_answer"], ans), renamed, text)
            self.assertTrue(cmp.service_answers(cases["resolve.discover_limit_above_100_host_answers"][
                "service_answer"], ans), text)

    def test_load_real_table(self):
        cases, info = cmp.load_by_design(cmp.DEFAULT_BY_DESIGN)
        # async's vendored table at 31a2e41f (e427d503's edit, case 131), and lane S8g's discovery cap (two renamed
        # cases, their old ids kept as aliases: 133)
        self.assertEqual(info["cases"], 133)
        # lane S7v: the vocabulary's wall_rule entry states budget_minus_reserve's threshold (B - 1 s, the name kept);
        # lane S8g: the discovery cap 1000 of async's fb66f1dd; integration stage 2: traverse.duplicate_seeds also
        # explains the merged seed's missing result
        self.assertEqual(info["sha256"], "d4b8636c748438ddbed1073a22379bcf39ae148ad90e879575dbdac38cd64f44")
        # v8 (VERIFY-D VD1): area G's cases are a table of their own, the v1-v7 table above byte-identical (async's
        # four of 2026-10-20 and the shadow run's three of 2026-10-21)
        self.assertFalse([cid for cid, c in cases.items() if "graphlet" in c])
        gcases, ginfo = cmp.load_by_design(cmp.DEFAULT_BY_DESIGN_G)
        self.assertEqual(sorted(gcases), ["graphlet.body_over_in_service_ceiling",
                                          "graphlet.body_over_in_service_ceiling.summary",
                                          "graphlet.compare_export_line",
                                          "graphlet.export_mgt_is_the_stored_body",
                                          "graphlet.max_bytes_above_service_cap",
                                          "graphlet.walks_arm_absent.one_arm",
                                          "graphlet.walks_arm_absent.two_arms"])
        self.assertTrue(ginfo["service_commit"].startswith("f994e3e2"))
        self.assertTrue(info["service_commit"].startswith("22d034ba"))
        # e427d503: traverse.duplicate_seeds explains only a merge, "N results for M seeds (N distinct)" (B5b-07),
        # and the field the merge removes from the MCP answer, results.1.duplicate (T-edge-0008, by design)
        entry, removed = cases["traverse.duplicate_seeds"]["explains_checks"]
        self.assertEqual(removed, {"check": "answer_fields",
                                   "detail_matches": r"field not in the answer: results\.1\.duplicate"})
        for outcome, detail, explained in (("not_exposed", "field not in the answer: results.1.duplicate", True),
                                           ("fail", "field not in the answer: results.1.duplicate", True),
                                           ("not_exposed", "field not in the answer: results.0.duplicate", False),
                                           ("not_exposed", "field not in the answer: results.1.duplicates", False),
                                           ("fail", "results.1.duplicate = False, expected True", False)):
            self.assertEqual(cmp.item_explained(removed, {"check": "answer_fields", "outcome": outcome,
                                                          "detail": detail}), explained, detail)
        self.assertEqual(entry, {"check": "invariants",
                                 "detail_matches": r"([1-9][0-9]*) results for [0-9]+ seeds \(\1 distinct\)"})
        for detail, explained in (("1 results for 2 seeds (1 distinct)", True),
                                  ("2 results for 3 seeds (2 distinct)", True),
                                  ("1 results for 3 seeds (2 distinct)", False),
                                  ("5 results for 2 seeds (1 distinct)", False),
                                  ("1 results for 2 seeds", False), ("1 results for 2 seeds (1 distinct); x", False)):
            self.assertEqual(cmp.entry_explains(entry, {"check": "invariants", "outcome": "fail", "detail": detail}),
                             explained, detail)
        self.assertEqual(cases["pattern.unknown_field"]["agreement"], "not_implemented_mcp")
        self.assertEqual(sorted(c["id"] for c in cmp.answer_cases(cases)),
                         ["mismatch.above_kind_cap", "pattern.caps_zero", "pattern.deadline_ends_task",
                          "traverse.advanced_path_disabled", "traverse.empty_seeds"])
        # 31a2e41f: the deployment-gated refusal reads not_implemented_mcp (answer_match.mcp_outcome); no other case
        # states the key
        self.assertEqual([c["id"] for c in cases.values() if cmp.match_mark(c) != "service_by_design"],
                         ["traverse.advanced_path_disabled"])
        adv = cases["traverse.advanced_path_disabled"]
        self.assertEqual((adv["agreement"], adv["answer_match"]["mcp_outcome"], adv["answer_match"]["direct"],
                          adv["answer_match"]["direct_mismatch"]),
                         ("not_implemented_mcp", "not_implemented_mcp", {"http": 200}, "disagree"))
        self.assertEqual(cmp.load_by_design(""), ({}, None))
        for case in cases.values():
            self.assertIn(case["agreement"], ("by_design", "not_implemented_mcp", "service_only"))
            self.assertIn(case["service_answer"]["outcome"],
                          ("refused_at_submit", "answered", "clamped", "not_expressible", "relayed"))
            # every key the table states is one the comparator matches or knows as descriptive
            self.assertEqual(cmp.unknown_keys(case["direct_expectation"], cmp.DIRECT_KEYS, cmp.DIRECT_NOTES), [],
                             case["id"])
            self.assertEqual(cmp.unknown_keys(case["service_answer"], cmp.SERVICE_KEYS, cmp.SERVICE_NOTES), [],
                             case["id"])
            am = case.get("answer_match") or {}
            self.assertEqual(cmp.unknown_keys(am.get("direct"), cmp.ANSWER_DIRECT_KEYS), [], case["id"])
            self.assertEqual(cmp.unknown_keys(am.get("mcp"), cmp.ANSWER_MCP_KEYS), [], case["id"])

    def test_overheads(self):
        # review H1: the front cost pairs host-warm executions only: P-001 local r1 (6 ms, server 3) with
        # remote r1 (100 ms, server 50); P-002's local execution is only the cold one: unpaired
        fo = self.comp["scenarios"]["P-001"]["front_overhead"]
        self.assertEqual((fo["local_repeat"], fo["remote_repeat"]), (1, 1))
        self.assertAlmostEqual(fo["wall_ms"], 94.0)
        self.assertAlmostEqual(fo["net_ms"], 47.0)
        self.assertAlmostEqual(fo["server_ms"], 47.0)
        self.assertIsNone(self.comp["scenarios"]["P-002"]["front_overhead"])
        f = self.s["front_overhead"]
        # P-001 and R-001 (whose MCP execution took the cold one, so both direct arms ran warm)
        self.assertEqual(f["warm"]["wall_ms"]["n"], 2)
        self.assertEqual(f["scenarios_without_warm_pair"], 2)
        self.assertEqual(f["per_area_net"]["P"]["n"], 1)
        self.assertEqual(f["cold"]["direct-local"]["wall_ms"]["n"], 8)  # R-001's cold execution was the MCP one
        mo = self.comp["scenarios"]["P-001"]["mcp_overhead"]
        self.assertEqual(mo["repeat"], 1)
        self.assertAlmostEqual(mo["minus_remote_ms"], 590.0)      # 700 - remote r0 (the first warm one)
        self.assertAlmostEqual(mo["minus_remote_no_queue_ms"], 390.0)
        self.assertAlmostEqual(mo["backend_minus_remote_ms"], 190.0)
        self.assertAlmostEqual(mo["minus_local_ms"], 694.0)       # 700 - local r1
        self.assertAlmostEqual(mo["backend_minus_local_ms"], 294.0)
        self.assertIsNone(mo["cold_wall_ms"])
        r1 = self.comp["scenarios"]["R-001"]["mcp_overhead"]
        self.assertEqual(r1["cold_wall_ms"], 800.0)
        self.assertNotIn("minus_local_ms", r1)                    # no MCP warm execution
        self.assertEqual(self.s["mcp_phases"]["queue_ms"]["n"], 1)
        # only MCP executions that reached the host (a backend call on record) count towards the overhead
        self.assertEqual(self.s["mcp_overhead"]["minus_local_ms"]["n"], 1)
        self.assertEqual(self.s["mcp_overhead"]["minus_remote_ms"]["n"], 3)
        self.assertEqual(self.s["mcp_cold_wall_ms"]["n"], 1)
        self.assertNotIn("mcp_overhead", self.comp["scenarios"]["P-005"])

    def test_answered_case_needs_the_answer_on_record(self):
        # review N11: an answered / clamped case is verified only by a line that records the service's answer
        want = {"outcome": "answered"}
        answered = line("P-1", "mcp", 0, T0, 800.0)
        self.assertEqual(cmp.match_service(want, [answered], None), (True, "matches"))
        mark = line("P-1", "mcp", 0, T0, 1.0, outcome="service_by_design", http=None, reached=False)
        self.assertEqual(cmp.match_service(want, [mark], None)[0], None)
        self.assertEqual(cmp.match_service({"outcome": "clamped"}, [mark, answered], None)[0], True)
        refused = line("P-1", "mcp", 0, T0, 1.0, http=400, reached=False, error={"code": "x", "text": "y"})
        self.assertEqual(cmp.match_service(want, [refused], None)[0], False)

    def test_mcp_overhead_skips_a_refusal_at_submit(self):
        # review N2: a refusal at submit recorded as `fail` is no MCP warm execution, nor does it make the
        # direct-local execution after it warm
        refused = line("T-1", "mcp", 1, T0 + 3.0, 210.0, outcome="fail", http=400, reached=False,
                       error={"code": "traversal_refused", "text": "host overrun"})
        cold = line("T-1", "mcp", 0, T0, 900.0, phases={"queue_ms": 100.0, "backend_ms": 400.0})
        local = line("T-1", "direct-local", 0, T0 + 1.5, 300.0)
        rows = [cold, local, refused]
        sc.assign_host_cache(rows)
        self.assertEqual([r["host_cache"] for r in rows], ["cold", "warm", "not_reached"])
        mo = cmp.mcp_overhead({0: cold, 1: refused}, {0: local}, {})
        self.assertIsNone(mo["repeat"])            # no MCP warm execution reached the host
        self.assertEqual(mo["cold_wall_ms"], 900.0)
        self.assertNotIn("minus_local_ms", mo)
        self.assertIsNone(cmp.mcp_overhead({1: refused}, {0: local}, {}))

    def test_wall_without_connect(self):
        new = line("P-1", "direct-remote", 0, T0, 100.0)
        new["client"] = {"connect_ms": 60.0, "reused": False, "attempts": 1}
        self.assertEqual(cmp.wall_net(new), 100.0)      # this runner times on an established connection
        old = line("P-1", "direct-remote", 0, T0, 100.0)
        old["client"] = {"connect_ms": 60.0}
        self.assertEqual(cmp.wall_net(old), 40.0)       # an older runner's wall included the connect

    def test_resources(self):
        r = self.s["resources"]["P"]["direct-local"]
        self.assertEqual(r["peak_server_memory_bytes"], 1500)
        # P-001 joins the samples at t+0.5 and t+2.5; P-002 (t+4..t+4.005) has none in its window and takes
        # the nearest, t+6.5 (CPU 120)
        self.assertEqual(r["peak_cpu_pct"], 120.0)
        # P-005 (t+14) has no sample within 5 s: not joined
        self.assertEqual((r["executions_with_samples"], r["executions"]), (5, 6))
        p1 = [x for x in sc.read_jsonl(os.path.join(self.dir, "executions.jsonl"))[0]
              if sc.execution_key(x) == ("P-001", "direct-local", 0)][0]
        win = cmp.SampleIndex(sc.read_jsonl(os.path.join(self.dir, "samples.jsonl"))[0]).window(*cmp.line_window(p1))
        self.assertEqual([s["cpu_pct"] for s in win], [50.0, 80.0])
        rr = self.s["resources"]["R"]["direct-local"]
        self.assertEqual(rr["peak_cpu_pct"], 120.0)
        self.assertEqual(self.s["resources"]["P"]["direct-remote"]["executions_with_samples"], 3)
        idx = cmp.SampleIndex([sample(T0, 1.0, 1, 1.0, 1.0)])
        self.assertEqual(len(idx.window(T0 + 4.0, T0 + 4.5)), 1)   # nearest within 5 s
        self.assertEqual(len(idx.window(T0 + 9.0, T0 + 9.5)), 0)
        self.assertEqual(idx.window(None, None), [])

    def test_resources_read_a_non_finite_sample_as_none(self):
        # SD-018: the run's resource means never carry a NaN or inf figure (compare.json would hold a NaN token)
        ln = line("P-001", "direct-local", 0, T0, 2000.0)
        samples = [dict(sample(T0 + 0.5, 50.0, 1, 1.0, 10.0)), dict(sample(T0 + 1.0, 50.0, 1, float("nan"), 30.0)),
                   dict(sample(T0 + 1.5, 50.0, 1, 3.0, float("inf")))]
        r = cmp.resources([ln], cmp.SampleIndex(samples))
        self.assertEqual(r["executions_with_samples"], 1)
        self.assertEqual((r["mean_sra_all_cpu_pct"], r["mean_host_load1"]), (2.0, 20.0))

    def test_helpers(self):
        self.assertEqual(cmp.aggregate_outcome(["pass", "pass"]), "pass")
        self.assertEqual(cmp.aggregate_outcome(["pass", "fail"]), "fail")
        self.assertEqual(cmp.aggregate_outcome(["fail", "error"]), "error")
        self.assertEqual(cmp.aggregate_outcome(["pass", "service_by_design"]), "mixed")
        self.assertIsNone(cmp.aggregate_outcome([]))
        d, u, n, ne, ni = cmp.compare_checks(
            [{"check": "a", "outcome": "pass"}, {"check": "b", "outcome": "pass"},
             {"check": "e", "outcome": "not_exposed"}, {"check": "f", "outcome": "inconclusive"},
             {"check": "g", "outcome": "finding"}],
            [{"check": "a", "outcome": "fail"}, {"check": "c", "outcome": "pass"},
             {"check": "e", "outcome": "pass"}, {"check": "f", "outcome": "pass"}, {"check": "g", "outcome": "pass"}])
        self.assertEqual(d, [{"check": "a", "a": "pass", "b": "fail"}, {"check": "g", "a": "finding", "b": "pass"}])
        self.assertEqual(u, ["b", "c"])
        self.assertEqual((n, ne, ni), (2, 1, 1))
        self.assertEqual(cmp.aggregate_outcome(["pass", "inconclusive"]), "inconclusive")
        self.assertEqual(cmp.aggregate_outcome(["inconclusive", "finding"]), "finding")
        self.assertEqual(cmp.aggregate_outcome(["finding", "fail"]), "fail")

    def test_main_writes_compare_json(self):
        self.assertEqual(cmp.main([self.dir, "--by-design", ""]), 0)
        data = sc.read_json(os.path.join(self.dir, "compare.json"))
        self.assertIsNone(data["by_design_table"])
        # without the table, case-referencing scenarios fall back to the generic comparison
        self.assertEqual(data["scenarios"]["P-005"]["mcp_agreement"]["class"], "disagree")


class TriageFixesTest(unittest.TestCase):
    """Lane S5 (the verified triage of the full run's first part): the by-design matchers (B1-B7) and the keys async
    adds to its table (answer_match.mcp field, answer_match.direct slot_code, direct_expectation http lists and
    text_prefix; a stated key the comparator does not know is never ignored)."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = self.tmp.name
        os.makedirs(os.path.join(self.dir, "raw"))

    def tearDown(self):
        self.tmp.cleanup()

    def raw_line(self, sid, arm, repeat, http=200, body=None, error=None, outcome="pass", checks=None):
        ln = line(sid, arm, repeat, T0, 10.0, outcome=outcome, http=http, error=error, checks=checks)
        sc.write_gz_json(os.path.join(self.dir, ln["raw"]), {"answer": {"http": http, "body": body, "error": error}})
        return ln

    def test_field_alternatives_and_dotted_placeholders(self):
        # B1 (G08): "predicate_strands | predicate_scope"; B2 (G07): "predicate<path of the first fault>"
        alt = {"id": "pattern.predicate_options_without_predicate", "agreement": "by_design",
               "direct_expectation": {"http": None, "code": None, "text": None, "behaviour": "not quoted"},
               "service_answer": {"outcome": "refused_at_submit", "code": "pattern_invalid",
                                  "text_contains": "it acts only with a predicate",
                                  "field": "predicate_strands | predicate_scope"}}
        m = line("Q-1", "mcp", 0, T0, 1.0, outcome="fail", http=400,
                 error={"code": "pattern_invalid", "field": "predicate_scope",
                        "text": "predicate_scope says ...: it acts only with a predicate"})
        self.assertEqual(cmp.match_service(alt["service_answer"], [m], None), (True, "matches"))
        dotted = {"outcome": "refused_at_submit", "code": "pattern_invalid", "text_contains": "predicate",
                  "field": "predicate<path of the first fault>"}
        q = line("Q-2", "mcp", 0, T0, 1.0, outcome="fail", http=400,
                 error={"code": "pattern_invalid", "field": "predicate.at_least.n",
                        "text": "predicate.at_least.n: an integer from 1 to the list's length (2), got 3"})
        self.assertEqual(cmp.match_service(dotted, [q], None)[0], True)
        self.assertEqual(cmp.match_service(dict(dotted, field="caps.<cap>"), [q], None)[0], False)

    def test_direct_http_lists_and_text_prefix(self):
        # B3 (G09, G10): http as a list (a host answering 200 or 503), text_prefix as is (no digit folding)
        d200 = self.raw_line("P-1", "direct-local", 0, http=200, body={"patterns": []})
        d503 = self.raw_line("P-1", "direct-local", 1, http=503, error={"code": "deadline", "text": "late"})
        self.assertEqual(cmp.match_direct({"http": [200, 503]}, d200, self.dir)[0], True)
        self.assertEqual(cmp.match_direct({"http": [200, 503]}, d503, self.dir)[0], True)
        self.assertEqual(cmp.match_direct({"http": [200]}, d503, self.dir)[0], False)
        nope = self.raw_line("R-1", "direct-local", 0, http=400,
                             error={"code": None, "text": "Label not found in the annotation: 'NOPE_0000.1'"})
        self.assertEqual(cmp.match_direct({"http": 400, "text_prefix": "Label not found in the annotation: "},
                                          nope, self.dir)[0], True)
        self.assertEqual(cmp.match_direct({"text_prefix": "Label not found in the annotation: 'NOPE_0"}, nope,
                                          self.dir)[0], True)
        self.assertEqual(cmp.match_direct({"text_prefix": "Label not found in the annotation: 'NOPE_9"}, nope,
                                          self.dir)[0], False)            # as is: no digit folding
        self.assertEqual(cmp.match_direct({"text_prefix": "Duplicate label"}, nope, self.dir)[0], False)
        self.assertEqual(cmp.match_direct({"text": "Label not found in the annotation: 'not-a-label-xyz'"}, nope,
                                          self.dir)[0], False)            # the whole text: label-specific
        # a stated key the comparator does not know: unverified, named; a null one is not stated
        res = cmp.match_direct({"http": 400, "text_regex": "Label .*"}, nope, self.dir)
        self.assertIsNone(res[0])
        self.assertIn("text_regex", res[1])
        self.assertEqual(cmp.match_direct({"http": 400, "text_regex": None, "behaviour": "x", "text_source": "y"},
                                          nope, self.dir)[0], True)

    def test_slot_codes(self):
        # B7: a direct side's code may be a pattern slot's; slot_code is a slot's only (not the top-level error)
        slot = self.raw_line("M-1", "direct-local", 0, http=200,
                             body={"patterns": [{"id": "w", "error": {"code": "mismatches_above_cap",
                                                                      "message": "pattern: max_mismatches 2 above"}}]})
        top = self.raw_line("M-2", "direct-local", 0, http=400,
                            error={"code": "mismatches_above_cap", "text": "pattern: max_mismatches 2 above"})
        self.assertEqual(cmp.match_direct({"http": 200, "code": "mismatches_above_cap"}, slot, self.dir)[0], True)
        self.assertEqual(cmp.match_direct({"slot_code": "mismatches_above_cap"}, slot, self.dir)[0], True)
        self.assertEqual(cmp.match_direct({"slot_code": "mismatches_above_cap"}, top, self.dir)[0], False)
        self.assertEqual(cmp.match_direct({"code": "mismatches_above_cap"}, top, self.dir)[0], True)

    def test_unstated_host_answer_judged_by_the_status_check(self):
        # B4 (G11): a case stating no host answer: the scenario's own status check on the reference direct line
        case = {"id": "pattern.duplicate_ids", "agreement": "by_design",
                "direct_expectation": {"http": None, "code": None, "text": None, "behaviour": "not quoted"},
                "service_answer": {"outcome": "refused_at_submit", "code": "pattern_invalid"}}
        m = [line("P-1", "mcp", 0, T0, 1.0, outcome="service_by_design", http=400,
                  error={"code": "pattern_invalid", "text": "pattern id 'same' is used twice"})]
        st = lambda o: [{"check": "status", "outcome": o, "detail": "http 200"}]  # noqa: E731
        res = cmp.by_design_class(case, line("P-1", "direct-local", 0, T0, 1.0, checks=st("pass")), m, None)
        self.assertEqual((res["class"], res["verified"], res["direct_by"]),
                         ("service_by_design", True, "scenario status check"))
        res = cmp.by_design_class(case, line("P-1", "direct-local", 0, T0, 1.0, checks=st("fail")), m, None)
        self.assertEqual((res["class"], res["side"]), ("disagree", "direct"))
        res = cmp.by_design_class(case, line("P-1", "direct-local", 0, T0, 1.0, checks=st("inconclusive")), m, None)
        self.assertEqual((res["class"], res["verified"], res["direct_verified"]), ("service_by_design", False, None))
        res = cmp.by_design_class(case, line("P-1", "direct-local", 0, T0, 1.0, checks=[]), m, None)
        self.assertEqual((res["class"], res["verified"]), ("service_by_design", False))
        # a stated key the comparator does not know is no "nothing stated": unverified, named
        odd = dict(case, direct_expectation={"http": None, "answer_shape": "x"})
        res = cmp.by_design_class(odd, line("P-1", "direct-local", 0, T0, 1.0, checks=st("pass")), m, None)
        self.assertEqual((res["class"], res["verified"], res["direct_verified"]), ("service_by_design", False, None))
        self.assertIn("answer_shape", res["direct_detail"])
        # the same for the service side
        odd = dict(case, service_answer={"outcome": "refused_at_submit", "code": "pattern_invalid", "text_regex": "x"})
        res = cmp.by_design_class(odd, line("P-1", "direct-local", 0, T0, 1.0, checks=st("pass")), m, None)
        self.assertEqual((res["class"], res["mcp_verified"]), ("service_by_design", None))
        self.assertIn("text_regex", res["mcp_detail"])

    def test_host_state_executions_read_neither_way(self):
        # B5 (G23): an MCP execution inconclusive by a host-state code is skipped when a named case is judged
        job_path = {"id": "resolve.job_path", "agreement": "by_design", "direct_expectation": {"http": 200},
                    "service_answer": {"outcome": "answered"}}
        d = line("R-1", "direct-local", 0, T0, 1.0)
        timeout = line("R-1", "mcp", 0, T0, 30250.0, outcome="inconclusive", http=504,
                       error={"code": "resolve_timeout", "text": "the backend did not answer within 30.25 s"})
        ok = line("R-1", "mcp", 1, T0, 4500.0)
        res = cmp.by_design_class(job_path, d, [timeout, ok], None)
        self.assertEqual((res["class"], res["verified"], res["host_state_repeats"]), ("service_by_design", True, [0]))
        res = cmp.by_design_class(job_path, d, [timeout], None)
        self.assertEqual(res["class"], "not_comparable")   # B-09: no execution reads for or against the case
        self.assertIn("resolve_timeout", res["detail"])
        # judged when it is not inconclusive (a runner of before the rule recorded it as fail)
        res = cmp.by_design_class(job_path, d, [dict(timeout, outcome="fail"), ok], None)
        self.assertEqual((res["class"], res["side"]), ("disagree", "mcp"))
        self.assertIn("resolve_timeout", cmp.MCP_HOST_STATE_CODES)
        # merge_run keeps the host-state rule's inconclusive (its failing checks judge a host state, not a result)
        st_fail = [{"check": "status", "outcome": "fail", "detail": "http 504"}]
        rows = [dict(timeout, checks=st_fail),
                dict(line("R-1", "mcp", 1, T0, 1.0, outcome="inconclusive", http=400,
                          error={"code": "other", "text": "x"}), checks=st_fail)]
        changed = merge_run.stricter_outcomes(rows, {}, lambda c: "fail")
        self.assertEqual(([r["outcome"] for r in rows], changed), (["inconclusive", "fail"], {"mcp": 1}))
        # the runner's (or the re-evaluation's) inconclusive_reason, a by-design decision with its finding, and the
        # checks an answered case explains are left as decided (lane S5b)
        rule = lambda c: "fail" if any(x["outcome"] == "fail" for x in c) else (  # noqa: E731
            "finding" if any(x["outcome"] == "finding" for x in c) else "pass")
        found = [{"check": "status", "outcome": "fail", "detail": ""},
                 {"check": "budget_respected", "outcome": "finding", "detail": ""}]
        rows = [dict(line("R-2", "mcp", 0, T0, 1.0, outcome="inconclusive", checks=st_fail),
                     inconclusive_reason="no service-side timing: plausibility undecidable"),
                dict(line("R-2", "mcp", 1, T0, 1.0, outcome="inconclusive", checks=st_fail),
                     inconclusive_reason="checks"),
                dict(line("N-2", "mcp", 0, T0, 1.0, outcome="finding", checks=found),
                     service_by_design_case="pattern.deadline_ends_task", explained_checks=[]),
                dict(line("P-2", "mcp", 0, T0, 1.0, outcome="finding", checks=found),
                     service_by_design_case="dup", explained_checks=["status"]),
                dict(line("P-3", "mcp", 0, T0, 1.0, outcome="pass", checks=found),
                     service_by_design_case="job", explained_checks=[])]
        merge_run.stricter_outcomes(rows, {}, rule)
        self.assertEqual([r["outcome"] for r in rows], ["inconclusive", "fail", "finding", "finding", "fail"])

    def deadline_cases(self, **am):
        case = {"id": "pattern.deadline_ends_task", "agreement": "by_design",
                "direct_expectation": {"http": 503, "code": "deadline"},
                "service_answer": {"outcome": "relayed", "http": 502, "code": "search_failed",
                                   "text_contains": "pattern_deadline"},
                "answer_match": {"routes": ["POST /pattern"], "direct": {"http": 503, "code": "deadline"},
                                 "mcp": {"http": 502, "code": "search_failed", "text_contains": "pattern_deadline"}}}
        for k, v in am.items():
            case["answer_match"][k] = v
        return dict(CASES, **{case["id"]: case}), case

    def test_answer_cases_before_named_cases(self):
        # B6: the answer-triggered case takes its MCP executions out, the named case is judged on the rest
        cases, case = self.deadline_cases()
        s = scen("M-1", "M", by_design="case.slot")
        relayed = line("M-1", "mcp", 0, T0, 1.0, outcome="service_by_design", http=502,
                       error={"code": "search_failed", "text": "pattern_deadline: pattern: the answer could not be"})
        refused = line("M-1", "mcp", 1, T0, 1.0, outcome="service_by_design", http=400,
                       error={"code": "pattern_out_of_bounds", "text": "max_mismatches 2 is above the maximum 1"})
        dead = self.raw_line("M-1", "direct-local", 0, http=503, error={"code": "deadline", "text": "late"})
        slot = self.raw_line("M-1", "direct-local", 1, http=200,
                             body={"patterns": [{"error": {"code": "mismatches_above_cap", "message": "above"}}]})
        res = cmp.mcp_agreement({0: relayed, 1: refused}, {"direct-local": {0: dead, 1: slot}}, s, cases, self.dir)
        self.assertEqual((res["class"], res["by"], res["case"], res["repeats"]),
                         ("service_by_design", "answer", case["id"], [0]))
        self.assertEqual((res["other_repeats"]["class"], res["other_repeats"]["case"]), ("service_by_design", "case.slot"))
        # before lane S5 the named case judged the relayed deadline too and read a disagreement
        self.assertEqual(cmp.by_design_class(CASES["case.slot"], dead, [relayed, refused], self.dir)["class"], "disagree")
        # a rest that disagrees with the named case is the class, naming the answer case
        wrong = dict(refused, error={"code": "pattern_invalid", "text": "x"}, outcome="fail",
                     checks=[{"check": "status", "outcome": "fail", "detail": "http 400"},
                             {"check": "contains_source_context", "outcome": "fail", "detail": "no answer"}])
        res = cmp.mcp_agreement({0: relayed, 1: wrong}, {"direct-local": {0: dead, 1: slot}}, s, cases, self.dir)
        self.assertEqual((res["class"], res["answer_case"], res["answer_case_repeats"],
                          [d["check"] for d in res["differing_checks"]]),
                         ("disagree", case["id"], [0], ["contains_source_context", "status"]))
        # every MCP execution taken out: no rest to judge
        res = cmp.mcp_agreement({0: relayed}, {"direct-local": {0: dead}}, s, cases, self.dir)
        self.assertEqual((res["class"], res["verified"]), ("service_by_design", True))
        self.assertNotIn("other_repeats", res)
        # the rest is judged against the direct executions the answer case does not explain: with only the host's
        # deadline on record, the named case's direct side is unverified, not a disagreement; and the whole is not
        # verified then (B-04: VERIFY-REVIEW's reading of this very assertion, which read verified True before)
        res = cmp.mcp_agreement({0: relayed, 1: refused}, {"direct-local": {0: dead}}, s, cases, self.dir)
        self.assertEqual((res["class"], res["verified"], res["other_repeats"]["class"],
                          res["other_repeats"]["verified"], res["other_repeats"]["direct_detail"]),
                         ("service_by_design", False, "service_by_design", False, "no direct execution"))

    def test_answer_match_keys(self):
        # B7 and async's addendum: the mcp side's field and http lists, the direct side's slot_code; a stated key the
        # comparator does not know: the mcp side never triggers, the direct side is unverified, named
        cases, case = self.deadline_cases(
            mcp={"http": 400, "code": "pattern_invalid | pattern_out_of_bounds", "field": "caps.<cap>"},
            direct={"http": [200, 503], "slot_code": "mismatches_above_cap"})
        s = scen("P-1", "P")
        m = line("P-1", "mcp", 0, T0, 1.0, outcome="fail", http=400,
                 error={"code": "pattern_invalid", "field": "caps.max_paths",
                        "text": "caps.max_paths: Input should be greater than or equal to 1"})
        slot = self.raw_line("P-1", "direct-local", 0, http=200,
                             body={"patterns": [{"error": {"code": "mismatches_above_cap", "message": "x"}}]})
        top = self.raw_line("P-1", "direct-local", 1, http=200, error={"code": "mismatches_above_cap", "text": "x"})
        res = cmp.mcp_agreement({0: m}, {"direct-local": {0: slot}}, s, cases, self.dir)
        self.assertEqual((res["class"], res["by"], res["verified"]), ("service_by_design", "answer", True))
        res = cmp.mcp_agreement({0: m}, {"direct-local": {1: top}}, s, cases, self.dir)
        self.assertEqual((res["class"], res["verified"]), ("service_by_design", False))  # slot_code: slots only
        other_field = dict(m, error=dict(m["error"], field="max_mismatches"))
        self.assertIsNone(cmp.answer_case_of(cases, "POST /pattern", 400, other_field["error"]))
        self.assertIs(cmp.answer_case_of(cases, "POST /pattern", 400, m["error"]), case)
        self.assertIsNone(cmp.answer_case_of(cases, "POST /pattern", 409, m["error"]))  # another http
        # unknown keys
        cases_u, _ = self.deadline_cases(mcp={"http": 400, "code": "pattern_invalid", "text_regex": "caps"})
        self.assertIsNone(cmp.answer_case_of(cases_u, "POST /pattern", 400, m["error"]))
        self.assertIsNone(cmp.answer_side_matches({"http": 400, "text_regex": "caps"}, 400, {"pattern_invalid"}, []))
        cases_d, _ = self.deadline_cases(mcp={"http": 400, "code": "pattern_invalid"},
                                         direct={"http": 200, "body_hash": "abc"})
        res = cmp.mcp_agreement({0: m}, {"direct-local": {0: slot}}, s, cases_d, self.dir)
        self.assertEqual((res["class"], res["verified"], res["direct_verified"]), ("service_by_design", False, None))
        self.assertIn("body_hash", res["direct_detail"])
        # a field on the direct side, or a slot_code on the mcp side, is a key that side does not have
        cases_f, _ = self.deadline_cases(mcp={"http": 400, "code": "pattern_invalid", "slot_code": "x"})
        self.assertIsNone(cmp.answer_case_of(cases_f, "POST /pattern", 400, m["error"]))


def damaged_raw_files(d):
    """{name: relative path} of raw records a host crash or a cut copy can leave (SD-093): a truncated .gz, a JSON list,
    a record whose answer is no object, and an intact one with a body (the control)."""
    import gzip
    out = {}
    for name, data in (("truncated", gzip.compress(b'{"answer": {"body": {"patterns": []}}}')[:-12]),
                       ("list", gzip.compress(b'[1, 2]')),
                       ("answer-text", gzip.compress(b'{"answer": "x"}')),
                       ("intact", gzip.compress(json.dumps({"answer": {"http": 200, "body": {"patterns": [
                           {"error": {"code": "slot_code", "message": "slot text"}}]}}}).encode()))):
        rel = "raw/%s.json.gz" % name
        os.makedirs(os.path.join(d, "raw"), exist_ok=True)
        with open(os.path.join(d, rel), "wb") as fh:
            fh.write(data)
        out[name] = rel
    return out


class DamagedRawRecordTest(unittest.TestCase):
    def test_direct_error_parts_of_a_damaged_raw_record(self):
        # SD-093 (guards at the crashing sites): a truncated .gz or a raw record that is no object (or whose answer
        # is none) is no body, its line's own error codes stay; it stopped compare with EOFError / AttributeError
        with tempfile.TemporaryDirectory() as d:
            files = damaged_raw_files(d)
            line = {"error": {"code": "line_code", "text": "line text"}}
            for name in ("truncated", "list", "answer-text"):
                self.assertEqual(cmp.direct_error_parts(dict(line, raw=files[name]), d),
                                 ({"line_code"}, set(), {"line text"}), name)
            self.assertEqual(cmp.direct_error_parts(dict(line, raw=files["intact"]), d),
                             ({"line_code"}, {"slot_code"}, {"line text", "slot text"}))


class ReviewFixesTest(unittest.TestCase):
    """Lane S5b: the verified review of lane S5 (B-01..B-11) and the triage's second part (item 5: the timing class),
    with async's agreed semantics (explains_checks, direct_mismatch, the table's load validation)."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = self.tmp.name
        os.makedirs(os.path.join(self.dir, "raw"))

    def tearDown(self):
        self.tmp.cleanup()

    def raw_line(self, sid, arm, repeat, http=200, body=None, error=None, outcome="pass", checks=None, **extra):
        ln = line(sid, arm, repeat, T0, 10.0, outcome=outcome, http=http, error=error, checks=checks,
                  digest=extra.pop("digest", "d"))
        ln.update(extra)
        sc.write_gz_json(os.path.join(self.dir, ln["raw"]), {"answer": {"http": http, "body": body, "error": error}})
        return ln

    CASES = {"job": {"id": "job", "agreement": "by_design", "direct_expectation": {"http": 200},
                     "service_answer": {"outcome": "answered"}, "explains_checks": []},
             "dup": {"id": "dup", "agreement": "by_design", "direct_expectation": {"http": 400},
                     "service_answer": {"outcome": "answered"}, "explains_checks": ["status"]},
             "short": {"id": "short", "agreement": "by_design", "direct_expectation": {"http": 400},
                       "service_answer": {"outcome": "refused_at_submit", "code": "sequence_too_short"}}}

    def resolve(self, names):
        return {"id": "R-1", "route": "POST /resolve", "request": {}, "mcp": {"tool": "traverse_resolve",
                                                                           "by_design": names},
                "expect": [{"check": "status", "args": {"http": 400}}]}

    def test_answered_cases_compare_the_checks(self):
        # B-01: an answered / clamped case is a note; the checks are compared, its explained checks apart
        st = lambda o: {"check": "status", "outcome": o, "detail": ""}  # noqa: E731
        runs = lambda o: {"check": "resolve_runs", "outcome": o, "detail": ""}  # noqa: E731
        d200 = self.raw_line("R-1", "direct-local", 0, checks=[st("pass"), runs("pass")])
        # P1: /resolve content failing on the MCP arm, job_path named: a disagreement, no longer by design
        m = self.raw_line("R-1", "mcp", 0, outcome="fail", checks=[st("pass"), runs("fail")])
        res = cmp.mcp_agreement({0: m}, {"direct-local": {0: d200}}, self.resolve("job"), self.CASES, self.dir)
        self.assertEqual((res["class"], res["by_design_notes"]), ("disagree", {"0": "job"}))
        ok = self.raw_line("R-1", "mcp", 0, checks=[st("pass"), runs("pass")])
        res = cmp.mcp_agreement({0: ok}, {"direct-local": {0: d200}}, self.resolve("job"), self.CASES, self.dir)
        self.assertEqual(res["class"], "agree")
        # the service answers what the host refuses: status explained by dup, nothing else differs
        d400 = self.raw_line("R-1", "direct-local", 0, http=400, error={"code": None, "text": "duplicate labels"},
                             checks=[st("pass"), runs("inconclusive")])
        m = self.raw_line("R-1", "mcp", 0, outcome="service_by_design", checks=[st("fail"), runs("pass")])
        res = cmp.mcp_agreement({0: m}, {"direct-local": {0: d400}}, self.resolve(["job", "dup"]), self.CASES,
                                self.dir)
        self.assertEqual((res["class"], res["by"], res["case"], res["verified"]),
                         ("service_by_design", "explains_checks", "dup", True))
        self.assertEqual([d["check"] for d in res["explained_differences"]], ["status"])
        # P4 (R-edge-0008): the MCP relays the host's refusal; job_path does not apply; both status checks pass
        rel = self.raw_line("R-1", "mcp", 0, http=400, outcome="pass", checks=[st("pass")],
                            error={"code": "resolve_rejected", "text": "HTTP 400: too short"})
        d = self.raw_line("R-1", "direct-local", 0, http=400, checks=[st("pass")], error={"text": "too short"})
        res = cmp.mcp_agreement({0: rel}, {"direct-local": {0: d}}, self.resolve("job"), self.CASES, self.dir)
        self.assertEqual(res["class"], "agree")
        # P2: a refusal case named after job_path; the service answers instead of refusing: a disagreement
        m = self.raw_line("R-1", "mcp", 0, outcome="fail", checks=[st("fail")])
        res = cmp.mcp_agreement({0: m}, {"direct-local": {0: d}}, self.resolve(["job", "short"]), self.CASES, self.dir)
        self.assertEqual((res["class"], [x["check"] for x in res["differing_checks"]]), ("disagree", ["status"]))
        refused = self.raw_line("R-1", "mcp", 0, http=400, outcome="service_by_design", checks=[st("fail")],
                                error={"code": "sequence_too_short", "text": "short"})
        res = cmp.mcp_agreement({0: refused}, {"direct-local": {0: d}}, self.resolve(["job", "short"]), self.CASES,
                                self.dir)
        self.assertEqual((res["class"], res["case"], res["verified"]), ("service_by_design", "short", True))

    def test_not_expressible_named_case(self):
        # a scenario naming a not_implemented_mcp case (resolve.not_expressible.bounds_max_query_bp): a difference
        # between the arms is the expressibility gap, the class not_implemented_mcp (R-budget-0001, R-twin-0002..5)
        cases = dict(self.CASES, nim={"id": "nim", "agreement": "not_implemented_mcp",
                                      "service_answer": {"outcome": "not_expressible"}})
        st = [{"check": "status", "outcome": "pass", "detail": ""}]
        refused = self.raw_line("R-1", "mcp", 0, http=400, outcome="service_by_design", checks=st,
                                error={"code": "sequence_too_short", "text": "short"})
        d200 = self.raw_line("R-1", "direct-local", 0, checks=st)
        res = cmp.mcp_agreement({0: refused}, {"direct-local": {0: d200}}, self.resolve(["job", "short", "nim"]),
                                cases, self.dir)
        self.assertEqual((res["class"], res["by"], res["disagreement"]["case"]), ("not_implemented_mcp", "table", "short"))
        d400 = self.raw_line("R-1", "direct-local", 0, http=400, checks=st, error={"text": "short"})
        res = cmp.mcp_agreement({0: refused}, {"direct-local": {0: d400}}, self.resolve(["job", "short", "nim"]),
                                cases, self.dir)
        self.assertEqual((res["class"], res["case"]), ("service_by_design", "short"))
        # only a named case's DIRECT side is the expressibility gap; a check-level disagreement (an MCP
        # backend_fault against a direct pass, as R-twin-0003 / R-edge-0013 were probed) stays a disagreement
        fault = self.raw_line("R-1", "mcp", 0, http=502, outcome="fail", error={"code": "backend_fault", "text": "x"},
                              checks=[{"check": "status", "outcome": "fail", "detail": "http 502"}])
        for named in (["job", "nim"], ["job", "short", "nim"]):
            res = cmp.mcp_agreement({0: fault}, {"direct-local": {0: d200}}, self.resolve(named), cases, self.dir)
            self.assertEqual((res["class"], res.get("side")), ("disagree", None), named)
        # every disagreement must be the gap, not the first: the named case's direct-side gap on one repeat and a
        # backend fault on another, in either order, stay a disagreement (the fault's)
        refused = self.raw_line("R-1", "mcp", 0, http=400, outcome="service_by_design", checks=st,
                                error={"code": "sequence_too_short", "text": "short"})
        f1 = self.raw_line("R-1", "mcp", 1, http=502, outcome="fail", error={"code": "backend_fault", "text": "x"},
                           checks=[{"check": "status", "outcome": "fail", "detail": "http 502"}])
        for reps in ({0: refused, 1: f1}, {0: f1, 1: refused}):
            res = cmp.mcp_agreement(reps, {"direct-local": {0: d200}}, self.resolve(["job", "short", "nim"]), cases,
                                    self.dir)
            self.assertEqual(res["class"], "disagree", res)
            self.assertNotEqual(res.get("by"), "table")
        # the gap alone, on both repeats: not_implemented_mcp by table
        r1 = self.raw_line("R-1", "mcp", 1, http=400, outcome="service_by_design", checks=st,
                           error={"code": "sequence_too_short", "text": "short"})
        res = cmp.mcp_agreement({0: refused, 1: r1}, {"direct-local": {0: d200}}, self.resolve(["job", "short", "nim"]),
                                cases, self.dir)
        self.assertEqual((res["class"], res["by"]), ("not_implemented_mcp", "table"))

    def test_named_case_reads_every_direct_repeat(self):
        # a named refusal case's direct side is judged on every repeat of both direct arms, as an answer-triggered
        # case's is (before, only the reference arm's first repeat): a contradiction on any repeat disagrees (a 503
        # deadline too: the side states the host's answer), the arm and repeat named; a transport error reads
        # neither way
        st = [{"check": "status", "outcome": "pass", "detail": ""}]
        refused = self.raw_line("R-1", "mcp", 0, http=400, outcome="service_by_design", checks=st,
                                error={"code": "sequence_too_short", "text": "short"})
        l400 = self.raw_line("R-1", "direct-local", 0, http=400, checks=st, error={"text": "short"})
        r400 = self.raw_line("R-1", "direct-remote", 0, http=400, checks=st, error={"text": "short"})
        r200 = self.raw_line("R-1", "direct-remote", 1, checks=st)
        res = cmp.mcp_agreement({0: refused}, {"direct-local": {0: l400}, "direct-remote": {0: r400, 1: r200}},
                                self.resolve(["job", "short"]), self.CASES, self.dir)
        self.assertEqual((res["class"], res["side"], res["case"]), ("disagree", "direct", "short"))
        self.assertEqual(res["direct_detail"], "direct-remote repeat 1: http 200, the case expects 400; 2 other direct "
                                               "repeat(s) match")
        dead = self.raw_line("R-1", "direct-remote", 1, http=503, outcome="inconclusive", checks=st,
                             error={"code": "deadline", "text": "late"})
        gone = dict(self.raw_line("R-1", "direct-remote", 2, http=0, outcome="error", checks=[],
                                  error={"code": "transport", "text": "reset"}))
        res = cmp.mcp_agreement({0: refused}, {"direct-local": {0: l400}, "direct-remote": {0: r400, 2: gone}},
                                self.resolve(["job", "short"]), self.CASES, self.dir)
        self.assertEqual((res["class"], res["verified"]), ("service_by_design", True))
        res = cmp.mcp_agreement({0: refused}, {"direct-local": {0: l400}, "direct-remote": {0: r400, 1: dead}},
                                self.resolve(["job", "short"]), self.CASES, self.dir)
        self.assertEqual((res["class"], res["direct_detail"][:28]), ("disagree", "direct-remote repeat 1: http"))
        # one execution as before
        self.assertEqual(cmp.by_design_class(self.CASES["short"], l400, [refused], self.dir)["class"],
                         "service_by_design")

    def test_named_200_case_sets_load_aside(self):
        # a named case whose host answers 200: a 503 (deadline or busy) or 429 on one direct execution is the host's
        # load when another direct execution holds the side (verified); it contradicts only when none holds; any other
        # status still contradicts; a 400 side keeps every 503 a contradiction (the expressibility gap of
        # R-budget-0001 / -0002, R-twin-0002 stays not_implemented_mcp by table)
        st = [{"check": "status", "outcome": "pass", "detail": ""}]
        cases = dict(self.CASES, ok200={"id": "ok200", "agreement": "by_design", "direct_expectation": {"http": 200},
                                         "service_answer": {"outcome": "refused_at_submit", "code": "caps_zero"}})
        refused = self.raw_line("R-1", "mcp", 0, http=400, outcome="service_by_design", checks=st,
                                error={"code": "caps_zero", "text": "x"})
        l200 = self.raw_line("R-1", "direct-local", 0, checks=st)
        loads = [self.raw_line("R-1", "direct-remote", 0, http=503, outcome="inconclusive", checks=st,
                               error={"code": "deadline", "text": "late"}),
                 self.raw_line("R-1", "direct-remote", 1, http=503, outcome="inconclusive", checks=st,
                               error={"code": "busy", "text": "busy"}),
                 self.raw_line("R-1", "direct-remote", 2, http=429, outcome="inconclusive", checks=st,
                               error={"code": "rate_limited", "text": "slow down"})]
        named = self.resolve(["ok200"])
        for load in loads:
            res = cmp.mcp_agreement({0: refused}, {"direct-local": {0: l200}, "direct-remote": {load["repeat"]: load}},
                                    named, cases, self.dir)
            self.assertEqual((res["class"], res["verified"]), ("service_by_design", True), load["http"])
            self.assertIn("1 direct repeat(s) answering 503 / 429 set aside (another matches)", res["direct_detail"])
        # no direct execution holds: the 503 contradicts
        res = cmp.mcp_agreement({0: refused}, {"direct-remote": {0: loads[0]}}, named, cases, self.dir)
        self.assertEqual((res["class"], res["side"]), ("disagree", "direct"))
        # another status contradicts whatever holds; it is named before a load answer
        r400 = self.raw_line("R-1", "direct-remote", 3, http=400, outcome="fail", checks=st, error={"text": "no"})
        res = cmp.mcp_agreement({0: refused}, {"direct-local": {0: l200}, "direct-remote": {0: loads[0], 3: r400}},
                                named, cases, self.dir)
        self.assertEqual((res["class"], res["direct_detail"][:28]), ("disagree", "direct-remote repeat 3: http"))
        # a 400 side (sequence_too_long's shape): every 503 contradicts, the non-load contradiction named first
        short = self.raw_line("R-1", "mcp", 0, http=400, outcome="service_by_design", checks=st,
                              error={"code": "sequence_too_short", "text": "short"})
        dl = self.raw_line("R-1", "direct-local", 0, http=503, outcome="inconclusive", checks=st,
                           error={"code": "deadline", "text": "late"})
        r200 = self.raw_line("R-1", "direct-remote", 0, checks=st)
        res = cmp.mcp_agreement({0: short}, {"direct-local": {0: dl}, "direct-remote": {0: r200}},
                                self.resolve(["job", "short"]), self.CASES, self.dir)
        self.assertEqual((res["class"], res["direct_detail"][:28]), ("disagree", "direct-remote repeat 0: http"))
        nim = dict(self.CASES, nim={"id": "nim", "agreement": "not_implemented_mcp",
                                     "service_answer": {"outcome": "not_expressible"}})
        rd = self.raw_line("R-1", "direct-remote", 0, http=503, outcome="inconclusive", checks=st,
                           error={"code": "deadline", "text": "late"})
        res = cmp.mcp_agreement({0: short}, {"direct-local": {0: dl}, "direct-remote": {0: rd}},
                                self.resolve(["job", "short", "nim"]), nim, self.dir)
        self.assertEqual((res["class"], res["by"]), ("not_implemented_mcp", "table"))

    def answer_case(self, direct, dm, http=None):
        case = {"id": "a", "agreement": "by_design", "direct_expectation": {},
                "service_answer": {"outcome": "refused_at_submit", "code": "pattern_out_of_bounds"},
                "answer_match": {"routes": ["POST /pattern"], "direct": direct,
                                 "mcp": {"http": http or 400, "code": "pattern_out_of_bounds"}}}
        if dm is not None:
            case["answer_match"]["direct_mismatch"] = dm
        return {"a": case}

    def test_direct_mismatch(self):
        # B-03 (async's a3692e2c): what a direct answer outside the case's direct side means
        s = scen("M-1", "M")
        m = self.raw_line("M-1", "mcp", 0, http=400, outcome="service_by_design",
                          error={"code": "pattern_out_of_bounds", "text": "above the cap"})
        stale = self.raw_line("M-1", "direct-local", 0, http=200, body={"patterns": [{"id": "w"}]})
        dead = self.raw_line("M-1", "direct-local", 1, http=503, error={"code": "deadline", "text": "late"})
        refusal = self.raw_line("M-1", "direct-local", 2, http=400, error={"code": "invalid_request", "text": "x"})
        transport = self.raw_line("M-1", "direct-local", 3, http=0, error={"code": "transport", "text": "reset"})
        slot = {"http": 200, "slot_code": "mismatches_above_cap"}
        cls = lambda dm, d: cmp.mcp_agreement({0: m}, {"direct-local": d}, s, self.answer_case(slot, dm), self.dir)  # noqa: E731
        # disagree: any direct HTTP answer outside the direct side (a stale cap: the host answers the slot)
        res = cls("disagree", {0: stale})
        self.assertEqual((res["class"], res["side"], res["direct_mismatch"]), ("disagree", "direct", "disagree"))
        self.assertEqual(cls("disagree", {1: dead})["class"], "disagree")
        # timing: a direct 200 or 503 deadline is unverified, a 4xx disagrees
        res = cls("timing", {0: stale})
        self.assertEqual((res["class"], res["verified"]), ("service_by_design", False))
        self.assertEqual(cls("timing", {1: dead})["class"], "service_by_design")
        self.assertEqual(cls("timing", {2: refusal})["class"], "disagree")
        # a transport error or no direct answer: unverified whatever the value; no direct_mismatch: unverified
        for dm in ("disagree", "timing", None):
            self.assertEqual(cls(dm, {3: transport})["class"], "service_by_design", dm)
            self.assertEqual(cmp.mcp_agreement({0: m}, {}, s, self.answer_case(slot, dm), self.dir)["class"],
                             "service_by_design", dm)
        self.assertEqual(cls(None, {0: stale})["class"], "service_by_design")
        # the direct side matched: verified
        good = self.raw_line("M-1", "direct-local", 4, http=200,
                             body={"patterns": [{"error": {"code": "mismatches_above_cap", "message": "cap 1"}}]})
        res = cls("disagree", {4: good})
        self.assertEqual((res["class"], res["verified"]), ("service_by_design", True))

    def test_table_validation_and_routes(self):
        # B-08 (async's 48dcbc1d): the table forms the runner refuses at load are refused here too
        base = {"id": "a", "agreement": "by_design", "service_answer": {"outcome": "refused_at_submit"}}
        bad = {"no routes": {"answer_match": {"direct": {}, "mcp": {"http": 400}}},
               "mcp http list": {"answer_match": {"routes": ["POST /pattern"], "direct": {}, "mcp": {"http": [400, 502]}}},
               "mcp slot_code": {"answer_match": {"routes": ["POST /pattern"], "direct": {}, "mcp": {"slot_code": "x"}}},
               "unknown key": {"answer_match": {"routes": ["POST /pattern"], "direct": {}, "mcp": {}, "when": 1}},
               "unknown side key": {"answer_match": {"routes": ["POST /pattern"], "direct": {"text_regex": "x"},
                                                     "mcp": {}}},
               "direct_mismatch": {"answer_match": {"routes": ["POST /pattern"], "direct": {}, "mcp": {},
                                                    "direct_mismatch": "unverified"}},
               "explains on a refusal": {"explains_checks": ["status"]},
               "explains not a list": {"service_answer": {"outcome": "answered"}, "explains_checks": "status"}}
        for what, extra in bad.items():
            self.assertTrue(cmp.table_problems([dict(base, **extra)]), what)
        self.assertEqual(cmp.table_problems([dict(base, service_answer={"outcome": "clamped"},
                                                  explains_checks=["cut"])]), [])
        path = os.path.join(self.dir, "t.json")
        sc.write_json_atomic(path, {"cases": [dict(base, **bad["mcp http list"])]})
        with self.assertRaises(cmp.TableError):
            cmp.load_by_design(path)
        self.assertEqual(cmp.main([self.dir, "--no-merge", "--by-design", path]), 4)
        real = cmp.load_by_design(cmp.DEFAULT_BY_DESIGN)[0]
        self.assertEqual(cmp.table_problems(list(real.values()), real.blocks), [])
        # an answer_match without routes never triggers (as the runner)
        case = {"a": dict(base, answer_match={"direct": {}, "mcp": {"http": 400, "code": "x"}})}
        self.assertIsNone(cmp.answer_case_of(case, "POST /pattern", 400, {"code": "x"}))
        # the MCP side's http is one status (equality); slot_code is not read there
        self.assertFalse(cmp.answer_side_matches({"http": [400]}, 400, {"x"}, [], known=cmp.ANSWER_MCP_KEYS,
                                                 http_list=False))
        self.assertIsNone(cmp.answer_side_matches({"http": 400, "slot_code": "x"}, 400, set(), [],
                                                  known=cmp.ANSWER_MCP_KEYS, http_list=False))

    def test_known_keys_first(self):
        # B-05: a known key that does not hold is a mismatch whatever unknown key is stated beside it
        d = self.raw_line("P-1", "direct-local", 0, http=500, error={"code": "internal", "text": "boom"})
        self.assertEqual(cmp.match_direct({"http": 400, "code": "invalid_request", "text_regex": "x"}, d, self.dir)[0],
                         False)
        self.assertIsNone(cmp.match_direct({"http": 500, "text_regex": "x"}, d, self.dir)[0])
        m = self.raw_line("P-1", "mcp", 0, http=400, error={"code": "totally_other", "text": "y"})
        self.assertFalse(cmp.match_service({"outcome": "refused_at_submit", "code": "pattern_invalid",
                                            "text_regex": "x"}, [m], self.dir)[0])

    def test_host_states(self):
        # B-09: all MCP executions host-state inconclusive: not comparable; a code the scenario expects is no host
        # state; A-11 / R2 / R3: judged on the service-side wall: a known wall below the threshold is judged by the
        # checks; no service-side wall: undecidable (rule (c)), read neither way
        s = {"id": "R-1", "route": "POST /resolve", "request": {"bounds": {"time_budget_ms": 20000}},
             "mcp": {"tool": "traverse_resolve", "arguments": {"time_budget_ms": 20000}},
             "expect": [{"check": "status", "args": {"http": 200}}]}
        d = self.raw_line("R-1", "direct-local", 0)
        to = self.raw_line("R-1", "mcp", 0, http=504, outcome="inconclusive",
                           error={"code": "resolve_timeout", "text": "the backend did not answer within 30.25 s"},
                           checks=[{"check": "status", "outcome": "fail", "detail": "http 504"}], wall_ms=30979.5,
                           mcp_phases={"backend_ms": 30375.0})
        self.assertTrue(cmp.host_state_line(to, self.dir, s))
        res = cmp.mcp_agreement({0: to}, {"direct-local": {0: d}}, s, {}, self.dir)
        self.assertEqual((res["class"], res["host_state_repeats"]), ("not_comparable", [0]))
        expects = dict(s, expect=[{"check": "status", "args": {"http": 504, "code": "resolve_timeout"}}])
        self.assertFalse(cmp.host_state_line(to, self.dir, expects))
        early = dict(to, mcp_phases={"backend_ms": 800.0})
        self.assertFalse(cmp.host_state_line(early, self.dir, s))
        self.assertEqual(cmp.mcp_agreement({0: early}, {"direct-local": {0: d}}, s, {}, self.dir)["class"],
                         "disagree")
        untimed = dict(to, mcp_phases={})    # the client wall is never the wall: undecidable
        self.assertTrue(cmp.host_state_line(untimed, self.dir, s))
        self.assertEqual(cmp.inconclusive_reason(untimed, s, self.dir),
                         "no service-side timing: plausibility undecidable")
        self.assertEqual(cmp.inconclusive_reason(to, s, self.dir), "resolve_timeout")
        self.assertEqual(cmp.inconclusive_reason(dict(to, inconclusive_reason="checks"), s, self.dir), "checks")
        # G22: the relayed /resolve deadline is a host state at budget - 1 s or later
        rdl = dict(to, error={"code": "traversal_resolve_deadline", "text": "within bounds.time_budget_ms (20000 ms"},
                   mcp_phases={"backend_ms": 22522.2})
        self.assertTrue(cmp.host_state_line(rdl, None, s))
        self.assertFalse(cmp.host_state_line(dict(rdl, mcp_phases={"backend_ms": 5000.0}), None, s))

    def test_wall_rule_in_the_comparison(self):
        # R1: pattern.deadline_ends_task (wall_rule budget_minus_reserve) applies only at a plausible service-side wall
        cases, _ = cmp.load_by_design(cmp.DEFAULT_BY_DESIGN)
        s = dict(scen("N-1", "N"), request={"time_budget_ms": 60000})
        st = [{"check": "status", "outcome": "fail", "detail": "http 502"}]
        err = {"code": "search_failed", "text": "pattern_deadline: pattern: the answer could not be written"}
        dead = self.raw_line("N-1", "direct-local", 0, http=503, error={"code": "deadline", "text": "late"},
                             outcome="inconclusive", checks=[{"check": "status", "outcome": "inconclusive",
                                                              "detail": ""}])
        def relay(ms):
            # the failing attempt's own wall: the deadline-ended task's one HTTP try in the raw record's diagnostic
            # log (never mcp_phases.backend_ms, a sum over the search's tasks)
            ln = self.raw_line("N-1", "mcp", 0, http=502, error=err, outcome="fail", checks=st,
                               mcp_phases={"backend_ms": 99999.0} if ms else {})
            if ms:
                slot = "\u29d7 Slot held for search s 1.0s (att 1, ok, acq 0.0, rel 0.1, http 1 fail 1 %.3fs " \
                       "backoff 0.0s, released)" % (ms / 1000.0)
                sc.write_gz_json(os.path.join(self.dir, ln["raw"]), {
                    "answer": {"http": 502, "body": None, "error": err},
                    "calls": [{"tool": "GET /searches/s/diagnostic-log", "payload": {"regular_tasks": [
                        {"error": "pattern_deadline: late", "diagnostic_log": [slot]}]}}]})
            return ln
        res = cmp.mcp_agreement({0: relay(60055.0)}, {"direct-local": {0: dead}}, s, cases, self.dir)
        self.assertEqual((res["class"], res["by"], res["verified"]), ("service_by_design", "answer", True))
        ok = self.raw_line("N-1", "direct-local", 0, checks=[{"check": "status", "outcome": "pass", "detail": ""}])
        res = cmp.mcp_agreement({0: relay(400.0)}, {"direct-local": {0: ok}}, s, cases, self.dir)
        self.assertEqual(res["class"], "disagree")   # held back: compared by its checks
        res = cmp.mcp_agreement({0: dict(relay(None), outcome="inconclusive")}, {"direct-local": {0: dead}}, s,
                                cases, self.dir)
        self.assertEqual((res["class"], res["undecidable_repeats"]), ("not_comparable", [0]))

    def test_undecidable_wall_rule_is_no_host_state(self):
        # B5b-09: a relayed /pattern deadline (code search_failed) whose wall rule is undecidable is no host state:
        # not_comparable by "undecidable plausibility", its reason replayed for a line recorded without it
        cases, _ = cmp.load_by_design(cmp.DEFAULT_BY_DESIGN)
        s = dict(scen("N-1", "N"), request={"time_budget_ms": 60000})
        st = [{"check": "status", "outcome": "fail", "detail": "http 502"}]
        err = {"code": "search_failed", "text": "pattern_deadline: pattern: the answer could not be written"}
        dead = self.raw_line("N-1", "direct-local", 0, http=503, error={"code": "deadline", "text": "late"},
                             outcome="inconclusive", checks=[{"check": "status", "outcome": "inconclusive",
                                                              "detail": ""}])
        why = "no service-side timing: plausibility undecidable"
        relay = self.raw_line("N-1", "mcp", 0, http=502, error=err, outcome="inconclusive", checks=st,
                              inconclusive_reason=why, mcp_phases={})
        self.assertFalse(cmp.host_state_line(relay, self.dir, s))
        res = cmp.mcp_agreement({0: relay}, {"direct-local": {0: dead}}, s, cases, self.dir)
        self.assertEqual((res["class"], res["by"], res["undecidable_repeats"]),
                         ("not_comparable", "undecidable plausibility", [0]))
        self.assertNotIn("host_state_repeats", res)
        bare = dict(relay)
        bare.pop("inconclusive_reason")
        self.assertEqual(cmp.inconclusive_reason(bare, s, self.dir, cases), why)
        self.assertEqual(cmp.inconclusive_reason(bare, s, self.dir), "checks")   # no table: the checks
        # a host-state line keeps its reading
        r = {"id": "R-1", "route": "POST /resolve", "request": {}, "mcp": {"tool": "traverse_resolve",
                                                                           "arguments": {"time_budget_ms": 20000}},
             "expect": [{"check": "status", "args": {"http": 200}}]}
        rdl = self.raw_line("R-1", "mcp", 0, http=503, outcome="inconclusive", inconclusive_reason=why,
                            error={"code": "traversal_resolve_deadline", "text": "late"}, mcp_phases={})
        self.assertTrue(cmp.host_state_line(rdl, self.dir, r))

    def test_client_wall_bound_in_the_comparison(self):
        # rule 2 (B5b-01, S5c-01): M-planted-0016's relayed deadline after 5 ms by the client's wall (no
        # service-side wall), ANONYMOUS calls: the answer-triggered case does not apply, the checks are compared (a
        # disagreement with a passing direct arm); a record that states no mcp_auth: undecidable
        cases, _ = cmp.load_by_design(cmp.DEFAULT_BY_DESIGN)
        s = dict(scen("M-1", "M"), request={"time_budget_ms": 60000},
                 mcp={"tool": "pattern_search", "arguments": {"time_budget_ms": 60000}})
        st = [{"check": "status", "outcome": "fail", "detail": "http 502"}]
        err = {"code": "search_failed", "text": "pattern_deadline: pattern: the answer could not be written"}
        relay = self.raw_line("M-1", "mcp", 0, http=502, error=err, outcome="fail", checks=st, mcp_phases={})
        doc = {"answer": {"http": 502, "body": None, "error": err},
               "calls": [{"tool": "pattern_search", "t0": 1000.0, "t1": 1000.002, "timed": True, "payload": {}},
                         {"tool": "get_status", "t0": 1000.002, "t1": 1000.005, "timed": True,
                          "payload": {"terminal": True, "state": "failed"}}]}
        sc.write_gz_json(os.path.join(self.dir, relay["raw"]), doc)
        ok = self.raw_line("M-1", "direct-local", 0, checks=[{"check": "status", "outcome": "pass", "detail": ""}])
        res = cmp.mcp_agreement({0: relay}, {"direct-local": {0: ok}}, s, cases, self.dir)
        self.assertEqual((res["class"], res["undecidable_repeats"]), ("not_comparable", [0]))
        for line in (dict(relay, mcp_auth="anonymous"), dict(relay, reevaluation={"mcp_auth_assumed": "anonymous"})):
            res = cmp.mcp_agreement({0: line}, {"direct-local": {0: ok}}, s, cases, self.dir)
            self.assertEqual(res["class"], "disagree")
        doc["mcp_auth"] = "anonymous"
        sc.write_gz_json(os.path.join(self.dir, relay["raw"]), doc)
        self.assertEqual(cmp.mcp_agreement({0: relay}, {"direct-local": {0: ok}}, s, cases, self.dir)["class"],
                         "disagree")
        doc["calls"][1]["t1"] = 1061.0    # 61 s: not below the threshold: undecidable
        sc.write_gz_json(os.path.join(self.dir, relay["raw"]), doc)
        res = cmp.mcp_agreement({0: dict(relay, outcome="inconclusive")}, {"direct-local": {0: ok}}, s, cases,
                                self.dir)
        self.assertEqual(res["class"], "not_comparable")

    def test_table_of_the_reevaluation(self):
        # B-10: compare --reevaluated reads the table the re-evaluation replayed and warns on another one
        other = os.path.join(self.dir, "other.json")
        with open(cmp.DEFAULT_BY_DESIGN, "rb") as fh:
            data = fh.read()
        with open(other, "wb") as fh:
            fh.write(data)
        sha = hashlib.sha256(data).hexdigest()
        sc.write_json_atomic(os.path.join(self.dir, "reevaluated.json"), {"by_design": {"path": other, "sha256": sha}})
        self.assertEqual(cmp.table_for(self.dir, None, True), (other, []))
        path, warnings = cmp.table_for(self.dir, cmp.DEFAULT_BY_DESIGN, True)
        self.assertEqual((path, warnings), (cmp.DEFAULT_BY_DESIGN, []))   # the same bytes
        sc.write_json_atomic(other, {"cases": []})
        path, warnings = cmp.table_for(self.dir, None, True)
        self.assertEqual(path, other)
        self.assertIn("differs from the one the re-evaluation replayed", warnings[0])
        self.assertEqual(cmp.table_for(self.dir, "", True)[0], "")
        self.assertEqual(cmp.table_for(self.dir, None, False), (cmp.DEFAULT_BY_DESIGN, []))

    def pattern_body(self, counts, cut=None):
        return {"patterns": [{"id": "pn", "counts": {"contexts": counts}, "cut": cut}]}

    def test_timing_class(self):
        # VERIFY-2 item 5: a stated time stop on one direct side only, consistent with the complete answer
        st = lambda o: {"check": "status", "outcome": o, "detail": ""}  # noqa: E731
        rel = lambda o: {"check": "count_relation", "outcome": o, "detail": ""}  # noqa: E731
        exact6 = {"relation": "exact", "value": 6, "by_mismatches": {"1": {"relation": "exact", "value": 6}}}
        least2 = {"relation": "at_least", "value": 2, "by_mismatches": {"1": {"relation": "at_least", "value": 2}}}
        complete = self.raw_line("N-1", "direct-remote", 0, body=self.pattern_body(exact6), checks=[st("pass"),
                                                                                                    rel("pass")],
                                 digest="b")
        stopped = self.raw_line("N-1", "direct-local", 0, body=self.pattern_body(least2, {"reason": "time"}),
                                outcome="inconclusive", checks=[st("pass"), rel("inconclusive")], digest="a",
                                time_limited="pattern 0 time-limited")
        res = cmp.direct_agreement({0: stopped}, {0: complete}, "POST /pattern", self.dir)
        self.assertEqual(res["class"], "timing", res)
        self.assertIn("pattern 0 time-limited", res["timing_repeats"][0]["detail"])
        # the stopped answer contradicts the complete one (at_least 99 > exact 6): a disagreement
        bad = self.raw_line("N-1", "direct-local", 1, body=self.pattern_body(dict(least2, value=99)),
                            outcome="inconclusive", checks=[st("pass"), rel("inconclusive")], digest="a",
                            time_limited="pattern 0 time-limited")
        complete1 = dict(complete, repeat=1, raw=complete["raw"])
        res = cmp.direct_agreement({1: bad}, {1: complete1}, "POST /pattern", self.dir)
        self.assertEqual(res["class"], "disagree")
        self.assertIn("contexts at_least 99 vs exact 6: no common value", res["outcome_differs"][0]["detail"])
        # a 503 deadline carries no content; a finding or a fail on the stopped side is no timing
        dead = self.raw_line("N-1", "direct-local", 2, http=503, error={"code": "deadline", "text": "late"},
                             outcome="inconclusive", checks=[st("inconclusive"), rel("inconclusive")], digest="c")
        complete2 = dict(complete, repeat=2)
        self.assertEqual(cmp.direct_agreement({2: dead}, {2: complete2}, "POST /pattern", self.dir)["class"], "timing")
        finding = dict(dead, outcome="finding", checks=[st("inconclusive"), {"check": "budget_respected",
                                                                             "outcome": "finding", "detail": ""}])
        self.assertEqual(cmp.direct_agreement({2: finding}, {2: complete2}, "POST /pattern", self.dir)["class"],
                         "disagree")
        # both stopped (P-heavy-0001's pair: a 503 deadline against a time-limited 200 that passed): neither answer
        # is complete, so neither shows the other wrong, timing (0605227da); a failing side stays a disagreement
        both = dict(complete2, time_limited="pattern 0 time-limited")   # outcome pass: both state a stop
        res = cmp.direct_agreement({2: dead}, {2: both}, "POST /pattern", self.dir)
        self.assertEqual(res["class"], "timing", res)
        self.assertIn("both time-limited", res["timing_repeats"][0]["detail"])
        both_fail = dict(both, outcome="fail", checks=[st("pass"), rel("fail")])
        self.assertEqual(cmp.direct_agreement({2: dead}, {2: both_fail}, "POST /pattern", self.dir)["class"],
                         "disagree")
        # the extension for G27: r1 inconclusive only by its determinism check after a time-limited r0, same answer
        det = {"check": "deterministic", "outcome": "inconclusive",
               "detail": "no earlier repeat without a time limitation (a 503 deadline): nothing to compare"}
        l1 = self.raw_line("D-1", "direct-local", 1, outcome="inconclusive", checks=[st("pass"), det], digest="x")
        r1 = self.raw_line("D-1", "direct-remote", 1, checks=[st("pass"), dict(det, outcome="pass")], digest="x")
        self.assertEqual(cmp.direct_agreement({1: l1}, {1: r1}, "POST /pattern", self.dir)["class"], "identical")
        r1b = dict(r1, answer_digest="y")
        self.assertEqual(cmp.direct_agreement({1: l1}, {1: r1b}, "POST /pattern", self.dir)["class"], "disagree")
        l1["answer_digest"], r1["answer_digest"] = "x", "x"
        l1b = dict(l1, answer_digest="x")
        r1c = dict(r1, answer_digest="x")
        both_r = {0: dict(r1c, repeat=0, answer_digest="z"), 1: r1c}
        both_l = {0: dict(l1b, repeat=0, answer_digest="z", outcome="pass", checks=[st("pass")]), 1: l1b}
        res = cmp.direct_agreement(both_l, both_r, "POST /pattern", self.dir)
        self.assertEqual(res["class"], "identical")

    def test_timing_count_ranges(self):
        # B5b-02 (VERIFY-REVIEW's rule): each population count admits a range of true values (exact [v, v],
        # at_least [v, inf), bounds [lower, upper], unknown anything); the stopped and the complete count must overlap
        ex = lambda v: {"relation": "exact", "value": v}  # noqa: E731
        al = lambda v: {"relation": "at_least", "value": v}  # noqa: E731
        bd = lambda lo, hi: {"relation": "bounds", "value": lo, "lower": lo, "upper": hi}  # noqa: E731
        unk = {"relation": "unknown", "value": None}
        bad = [(bd(10, 12), ex(6)), (bd(2, 5), ex(6)), (ex(9), al(12)), (ex(9), bd(10, 12)), (al(13), bd(10, 12)),
               (bd(1, 3), bd(10, 12)), (al(99), ex(6)), (ex(5), ex(0)), ({"relation": "at_most", "value": 9}, ex(6)),
               (ex(6), {"relation": "about", "value": 6}), (ex(None), ex(6))]
        good = [(bd(5, 7), ex(6)), (ex(12), al(12)), (ex(11), bd(10, 12)), (al(12), bd(10, 12)),
                (bd(1, 10), bd(10, 12)), (al(2), ex(6)), (al(0), al(5)), (ex(6), ex(6)), (unk, ex(6)), (ex(6), unk),
                (al(3), unk)]
        for a, b in bad:
            self.assertTrue(cmp._count_problems({"contexts": a}, {"contexts": b}), (a, b))
        for a, b in good:
            self.assertEqual(cmp._count_problems({"contexts": a}, {"contexts": b}), [], (a, b))
        # population counts only: `tested` is work done (a stopped selection pass tests fewer), estimate no bound
        stopped = {"tested": ex(50), "labels": ex(3), "occurrences": al(4), "contexts": dict(ex(9), by_offset={
            "0": bd(4, 6)}), "estimate": ex(1)}
        complete = {"tested": ex(80), "labels": ex(3), "occurrences": ex(9), "contexts": dict(ex(9), by_offset={
            "0": ex(5)}), "estimate": ex(7)}
        self.assertEqual(cmp._count_problems(stopped, complete), [])
        stopped["contexts"]["by_offset"]["0"] = bd(6, 8)
        self.assertEqual(cmp._count_problems(stopped, complete),
                         ["contexts.by_offset.0 bounds [6, 8] vs exact 5: no common value"])
        # end to end, the recorded N-nearmiss-0030 r0 pair: level 1 bounds [10, 12] against the remote's exact 6 is
        # no timing (it read timing); at_least 2 is
        st = lambda o: {"check": "status", "outcome": o, "detail": ""}  # noqa: E731
        complete = self.raw_line("N-2", "direct-remote", 0, body=self.pattern_body(
            {"relation": "exact", "value": 6, "by_mismatches": {"1": ex(6)}}), checks=[st("pass")], digest="b")
        for counts, cls in (({"relation": "at_least", "value": 2, "by_mismatches": {"1": al(2)}}, "timing"),
                            ({"relation": "at_least", "value": 10, "by_mismatches": {"1": bd(10, 12)}}, "disagree")):
            stopped = self.raw_line("N-2", "direct-local", 0, body=self.pattern_body(counts, {"reason": "time"}),
                                    outcome="inconclusive", checks=[st("pass")], digest="a",
                                    time_limited="pattern 0 time-limited")
            self.assertEqual(cmp.direct_agreement({0: stopped}, {0: complete}, "POST /pattern", self.dir)["class"],
                             cls, counts)

    def test_timing_released_subset(self):
        # VERIFY-2 item 5's subset rule: a stopped answer's released contexts are among a complete side's that
        # states retrieval_complete
        st = lambda o: {"check": "status", "outcome": o, "detail": ""}  # noqa: E731
        ctx = lambda k, o: {"kmer": k, "offset": o, "strand": "+"}  # noqa: E731
        def body(results, complete, relation="exact", n=2):
            return {"patterns": [{"id": "p", "counts": {"contexts": {"relation": relation, "value": n}},
                                  "results": results, "retrieval_complete": complete}]}
        full = self.raw_line("P-1", "direct-remote", 0, body=body([ctx("A" * 31, 0), ctx("C" * 31, 1)], True),
                             checks=[st("pass")], digest="b")
        def stopped(results):
            return self.raw_line("P-1", "direct-local", 0, body=body(results, False, "at_least", len(results)),
                                 outcome="inconclusive", checks=[st("pass")], digest="a",
                                 time_limited="pattern 0 time-limited")
        self.assertEqual(cmp.direct_agreement({0: stopped([ctx("C" * 31, 1)])}, {0: full}, "POST /pattern",
                                              self.dir)["class"], "timing")
        res = cmp.direct_agreement({0: stopped([ctx("G" * 31, 0)])}, {0: full}, "POST /pattern", self.dir)
        self.assertEqual(res["class"], "disagree")
        self.assertIn("1 of 1 released results not in the complete list", res["outcome_differs"][0]["detail"])
        # a complete side that does not state retrieval_complete licenses nothing
        partial = self.raw_line("P-1", "direct-remote", 0, body=body([ctx("A" * 31, 0)], False, "at_least", 1),
                                checks=[st("pass")], digest="b")
        self.assertEqual(cmp.direct_agreement({0: stopped([ctx("G" * 31, 0)])}, {0: partial}, "POST /pattern",
                                              self.dir)["class"], "timing")

    def test_direct_mismatch_every_repeat(self):
        # B5b-04: with direct_mismatch every repeat of both direct arms is classified and a disagreeing one wins
        # over a match (mismatch.above_kind_cap: r0 matches, r1 a 200 without the slot, r2 a 400)
        s = scen("M-1", "M")
        m = self.raw_line("M-1", "mcp", 0, http=400, outcome="service_by_design",
                          error={"code": "pattern_out_of_bounds", "text": "above the cap"})
        good = self.raw_line("M-1", "direct-local", 0, http=200,
                             body={"patterns": [{"error": {"code": "mismatches_above_cap", "message": "cap 1"}}]})
        stale = self.raw_line("M-1", "direct-local", 1, http=200, body={"patterns": [{"id": "w"}]})
        refusal = self.raw_line("M-1", "direct-local", 2, http=400, error={"code": "invalid_request", "text": "x"})
        slot = {"http": 200, "slot_code": "mismatches_above_cap"}
        cls = lambda dm, d: cmp.mcp_agreement({0: m}, d, s, self.answer_case(slot, dm), self.dir)  # noqa: E731
        res = cls("disagree", {"direct-local": {0: good, 1: stale, 2: refusal}})
        self.assertEqual((res["class"], res["side"]), ("disagree", "direct"))
        self.assertIn("1 other direct repeat(s) match", res["direct_detail"])
        self.assertEqual(cls("disagree", {"direct-local": {1: stale, 2: refusal}})["class"], "disagree")
        # timing: the 200 is unverified, the 400 disagrees, whatever matched
        self.assertEqual(cls("timing", {"direct-local": {0: good, 2: refusal}})["class"], "disagree")
        res = cls("timing", {"direct-local": {0: good, 1: stale}})
        self.assertEqual((res["class"], res["verified"]), ("service_by_design", True))
        # both direct arms count: a disagreeing remote repeat is no longer hidden behind the reference arm
        remote = dict(refusal, arm="direct-remote", repeat=0)
        self.assertEqual(cls("disagree", {"direct-local": {0: good}, "direct-remote": {0: remote}})["class"],
                         "disagree")
        res = cls("disagree", {"direct-local": {0: good}, "direct-remote": {0: dict(good, arm="direct-remote")}})
        self.assertEqual((res["class"], res["verified"]), ("service_by_design", True))
        self.assertIn("direct-remote r0", res["direct_detail"])

    def test_unmatched_runner_marks_split(self):
        # B5b-05: r0 marked service_by_design by the runner (a refusal no case of the table matches), r1 failing
        # its invariants against a passing direct arm: the disagreement of r1 decides, with or without the table
        s = dict(scen("N-1", "N"), mcp={"tool": "pattern_search", "arguments": {}})
        ok = [{"check": "status", "outcome": "pass", "detail": ""}, {"check": "invariants", "outcome": "pass",
                                                                     "detail": ""}]
        d = {0: self.raw_line("N-1", "direct-local", 0, checks=ok), 1: self.raw_line("N-1", "direct-local", 1,
                                                                                     checks=ok)}
        marked = self.raw_line("N-1", "mcp", 0, http=400, outcome="service_by_design",
                               error={"code": "unlisted_refusal", "text": "x"},
                               checks=[{"check": "status", "outcome": "fail", "detail": "http 400"}])
        broken = self.raw_line("N-1", "mcp", 1, outcome="fail", checks=[
            ok[0], {"check": "invariants", "outcome": "fail", "detail": "bad"}])
        cases, _ = cmp.load_by_design(cmp.DEFAULT_BY_DESIGN)
        for table in ({}, cases):
            res = cmp.mcp_agreement({0: marked, 1: broken}, {"direct-local": d}, s, table, self.dir)
            self.assertEqual(res["class"], "disagree", table and "table")
        # no table: the mark alone stays an unverified by-design class
        res = cmp.mcp_agreement({0: marked}, {"direct-local": d}, s, {}, self.dir)
        self.assertEqual((res["class"], res["by"], res["verified"]), ("service_by_design", "mcp runner", False))
        # with the table a mark it does not reproduce is judged by its checks (B-06): status fail vs pass
        self.assertEqual(cmp.mcp_agreement({0: marked}, {"direct-local": d}, s, cases, self.dir)["class"],
                         "disagree")

    def test_host_state_readers_of_merge_and_compare(self):
        # SD-035: merge_run.host_state_inconclusive (does a verdict rule decide the line: a host state, or rule (c)
        # by its reason) and compare.host_state_line (is it a host-state repeat) ask different questions on the
        # reason path, and read a line recorded without inconclusive_reason alike: its raw record's error and the
        # run's References (the live service document over the table's blocks)
        import merge_run
        cases, _ = cmp.load_by_design(cmp.DEFAULT_BY_DESIGN)
        res = {"id": "R-1", "route": "POST /resolve", "expect": [{"check": "status", "args": {"http": 200}}],
               "mcp": {"tool": "traverse_resolve", "arguments": {"time_budget_ms": 20000}}}
        legacy = {"code": "resolve_timeout", "text": "the backend did not answer"}
        phases = {"backend_ms": 20000.0}
        os.makedirs(os.path.join(self.dir, "arms", "mcp"))

        def both(ln):
            refs = cmp.refs_of(cmp.bind_run(cases, self.dir))
            return (merge_run.host_state_inconclusive(ln, res, self.dir),
                    cmp.host_state_line(ln, self.dir, res, refs=refs))
        # a legacy resolve_timeout at 20 s against min(the write cap, budget + reserve + allowance) - 1 s: the table's
        # 10 s write cap makes it a host state (9 s), the run's live write cap of 40 s does not (29.25 s)
        ln = self.raw_line("R-1", "mcp", 0, http=504, error=legacy, outcome="inconclusive", mcp_phases=phases)
        self.assertEqual(both(ln), (True, True))
        sc.write_json_atomic(os.path.join(self.dir, "arms", "mcp", "run.json"), {"service_capabilities": {
            "document": {"service_settings": {"resolve": {"http_timeouts_s": {"write": 40}}}}}})
        self.assertEqual(both(ln), (False, False))
        # the code only in the raw record (the line's error lost): both read the raw answer's
        bare = self.raw_line("R-1", "mcp", 1, http=504, error=legacy, outcome="inconclusive", mcp_phases=phases)
        bare["error"] = None
        os.remove(os.path.join(self.dir, "arms", "mcp", "run.json"))
        self.assertEqual(both(bare), (True, True))
        # the reason path: a host-state reason with its code decides both; rule (c)'s undecidable reason on an
        # answer-triggered case (a relayed /pattern deadline, search_failed) is the runner's decision for the merge,
        # no host-state repeat for the comparison (B5b-09)
        hs = dict(ln, inconclusive_reason="resolve_timeout")
        self.assertEqual(both(hs), (True, True))
        relayed = dict(ln, error={"code": "search_failed", "text": "pattern_deadline: x"},
                       inconclusive_reason=sc.UNDECIDABLE_NO_WALL)
        self.assertEqual(both(relayed), (True, False))

    def test_by_design_finding_without_a_table(self):
        # SD-022 (B5b-10): with no table the runner's by-design finding (a relayed deadline whose budget_respected is
        # a finding) is its mark, as reevaluate keeps it: an unverified by-design part, not a disagreement check by
        # check; with the table the case is replayed
        s = dict(scen("P-1", "P"), mcp={"tool": "pattern_search", "arguments": {}})
        ok = [{"check": "status", "outcome": "pass", "detail": ""},
              {"check": "budget_respected", "outcome": "pass", "detail": ""}]
        d = {0: self.raw_line("P-1", "direct-local", 0, checks=ok)}
        finding = self.raw_line("P-1", "mcp", 0, http=502, outcome="finding", service_by_design_case="pattern.x",
                                error={"code": "search_failed", "text": "pattern_deadline: x"},
                                checks=[{"check": "status", "outcome": "fail", "detail": "http 502"},
                                        {"check": "budget_respected", "outcome": "finding", "detail": "late"}])
        res = cmp.mcp_agreement({0: finding}, {"direct-local": d}, s, {}, self.dir)
        self.assertEqual((res["class"], res.get("by"), res.get("verified")), ("service_by_design", "mcp runner", False))
        self.assertTrue(cmp.runner_by_design_finding(finding))
        self.assertFalse(cmp.runner_by_design_finding(dict(finding, outcome="fail")))
        self.assertFalse(cmp.runner_by_design_finding(dict(finding, service_by_design_case=None)))

    def test_explained_items(self):
        # B5b-06: explains_checks entries explain check ITEMS ({check, detail_matches}); an unexplained item of the
        # same check is compared (it was hidden behind the name)
        cases = {"res": {"id": "res", "agreement": "by_design", "direct_expectation": {"http": 200},
                         "service_answer": {"outcome": "answered"},
                         "explains_checks": [{"check": "resolve_runs", "detail_matches": "merged .*"}]}}
        s = {"id": "R-1", "route": "POST /resolve", "request": {}, "mcp": {"tool": "traverse_resolve",
                                                                           "by_design": ["res"]},
             "expect": [{"check": "status", "args": {"http": 200}}]}
        rr = lambda o, d: {"check": "resolve_runs", "outcome": o, "detail": d}  # noqa: E731
        st = {"check": "status", "outcome": "pass", "detail": ""}
        mine = [st, rr("fail", "merged runs"), rr("fail", "wrong run 3")]
        self.assertEqual(cmp.explained_items(mine, cmp.union_of([cases["res"]])), [1])
        m = self.raw_line("R-1", "mcp", 0, outcome="fail", checks=mine)
        d = self.raw_line("R-1", "direct-local", 0, checks=[st, rr("pass", ""), rr("pass", "")])
        res = cmp.mcp_agreement({0: m}, {"direct-local": {0: d}}, s, cases, self.dir)
        self.assertEqual(res["class"], "disagree")
        self.assertEqual(res["differing_checks"][0]["check"], "resolve_runs")
        self.assertEqual((res["explained_items"], res["explained_checks"]), ({"0": [1]}, {"0": ["resolve_runs"]}))
        # only the explained item fails: service_by_design by explains_checks
        m2 = self.raw_line("R-1", "mcp", 0, outcome="service_by_design", checks=mine[:2] + [rr("pass", "")])
        res = cmp.mcp_agreement({0: m2}, {"direct-local": {0: d}}, s, cases, self.dir)
        self.assertEqual((res["class"], res["by"]), ("service_by_design", "explains_checks"))
        # S5c-05: the items are computed with the table in hand; a line that records others (a re-evaluation with
        # another table) is noted, not trusted
        m3 = dict(m, explained_items=[1, 2])
        res = cmp.mcp_agreement({0: m3}, {"direct-local": {0: d}}, s, cases, self.dir)
        self.assertEqual((res["class"], res["explained_items"]), ("disagree", {"0": [1]}))
        self.assertEqual(res["explained_items_recorded"], {"0": [1, 2]})
        res = cmp.mcp_agreement({0: dict(m, explained_items=[])}, {"direct-local": {0: d}}, s, cases, self.dir)
        self.assertEqual((res["explained_items"], res["explained_items_recorded"]), ({"0": [1]}, {"0": []}))
        res = cmp.mcp_agreement({0: dict(m, explained_items=[1])}, {"direct-local": {0: d}}, s, cases, self.dir)
        self.assertNotIn("explained_items_recorded", res)

    def test_mcp_outcome(self):
        # async's 901893b7 / 31a2e41f: answer_match.mcp_outcome "service_by_design" (default) | "not_implemented_mcp";
        # load validation: an answer_match key, one of the two values, never on an answered / clamped case, never a
        # top-level case key. A match of a not_implemented_mcp case (named or answer-triggered) is counted under
        # not_implemented_mcp, never service_by_design; the pair's reading is kept under `pair`
        base = {"id": "a", "agreement": "not_implemented_mcp", "direct_expectation": {"http": 200},
                "service_answer": {"outcome": "refused_at_submit", "code": "preset_unavailable"},
                "answer_match": {"routes": ["POST /traverse"], "direct": {"http": 200},
                                 "mcp": {"http": 400, "code": "preset_unavailable", "field": "overrides.<knob path>",
                                         "text_contains": "which is not enabled on this deployment"},
                                 "direct_mismatch": "disagree", "mcp_outcome": "not_implemented_mcp"},
                "request": [{"path": "{field}", "present": True}]}
        self.assertEqual(cmp.table_problems([base]), [])
        self.assertEqual(cmp.match_mark(base), "not_implemented_mcp")
        self.assertEqual(cmp.match_mark(dict(base, answer_match=dict(base["answer_match"], mcp_outcome=None))),
                         "service_by_design")
        self.assertEqual(cmp.match_mark({"id": "b"}), "service_by_design")
        bad = dict(base, answer_match=dict(base["answer_match"], mcp_outcome="maybe"))
        self.assertEqual(cmp.table_problems([bad]), [
            "case 'a': answer_match.mcp_outcome 'maybe' is not one of ['service_by_design', 'not_implemented_mcp']"])
        ans = dict(base, service_answer={"outcome": "answered"}, explains_checks=[])
        self.assertIn("case 'a': answer_match.mcp_outcome on a answered case (the checks decide an answered / clamped "
                      "case)", cmp.table_problems([ans]))
        top = dict(base, mcp_outcome="not_implemented_mcp")
        self.assertEqual(cmp.table_problems([top]), [
            "case 'a': mcp_outcome is an answer_match key (answer_match.mcp_outcome), not a case key"])
        with self.assertRaises(cmp.TableError), tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "t.json")
            sc.write_json_atomic(path, {"cases": [top]})
            cmp.load_by_design(path)
        # the class: T-labels-0004's shape (v4: MCP 400 preset_unavailable on overrides.labels.max_seed_labels, the
        # direct arms 200), the case answer-triggered (the scenario names none)
        cases = {"a": base}
        s = {"id": "T-labels-0004", "route": "POST /traverse", "request": {}, "expect": [],
             "mcp": {"tool": "traverse_fetch", "arguments": {"overrides": {"labels": {"max_seed_labels": 10000},
                                                                           "branching": {"max_paths": 2}}}}}
        err = {"code": "preset_unavailable", "field": "overrides.labels.max_seed_labels",
               "text": "overrides.labels.max_seed_labels is on the advanced strategy path, which is not enabled on "
                       "this deployment"}
        m = self.raw_line("T-labels-0004", "mcp", 0, http=400, error=err, outcome="not_implemented_mcp",
                          checks=[{"check": "status", "outcome": "fail", "detail": "http 400"}],
                          service_by_design_case="a")
        d = self.raw_line("T-labels-0004", "direct-local", 0, checks=[{"check": "status", "outcome": "pass"}])
        res = cmp.mcp_agreement({0: m}, {"direct-local": {0: d}}, s, cases, self.dir)
        self.assertEqual((res["class"], res["by"], res["case"], res["repeats"]),
                         ("not_implemented_mcp", "mcp_outcome", "a", [0]))
        self.assertEqual((res["pair"]["class"], res["pair"]["verified"]), ("service_by_design", True))
        # a direct answer outside the case's direct side (direct_mismatch "disagree") is a disagreement on the direct
        # side, the case and the direct execution named; the MCP execution stays not_implemented_mcp (the agreed
        # compare decision of 2026-10-14)
        d400 = self.raw_line("T-labels-0004", "direct-local", 0, http=400, error={"code": "x", "text": "y"},
                             checks=[{"check": "status", "outcome": "fail"}])
        res = cmp.mcp_agreement({0: m}, {"direct-local": {0: d400}}, s, cases, self.dir)
        self.assertEqual((res["class"], res["side"], res["case"], res["mcp_outcome"], res["pair"]["class"]),
                         ("disagree", "direct", "a", "not_implemented_mcp", "disagree"))
        self.assertIn("direct-local repeat 0: http 400", res["direct_detail"])
        self.assertIn("the MCP execution(s) [0] read not_implemented_mcp (case a, answer_match.mcp_outcome)",
                      res["direct_detail"])
        named_a = dict(s, mcp=dict(s["mcp"], by_design=["a"]))
        res = cmp.mcp_agreement({0: m}, {"direct-local": {0: d400}}, named_a, cases, self.dir)
        self.assertEqual((res["class"], res["side"], res["case"]), ("disagree", "direct", "a"))
        self.assertTrue(res["direct_detail"].startswith("direct-local repeat 0: "), res["direct_detail"])
        # mixed with another repeat the checks compare: a disagreement decides (here the other repeat's own: its
        # status passes where the direct one fails; each repeat is classified on its own)
        m_other = self.raw_line("T-labels-0004", "mcp", 1, checks=[{"check": "status", "outcome": "pass"}])
        res = cmp.mcp_agreement({0: m, 1: m_other}, {"direct-local": {0: d400}}, s, cases, self.dir)
        self.assertEqual((res["class"], [d["check"] for d in res["differing_checks"]]), ("disagree", ["status"]))
        self.assertEqual(res["other_repeats"]["case"], "a")
        # the T-labels-0004 pattern: the advanced-path refusal on one repeat and an MCP fault on another, the direct
        # arm inside the case's side: the fault is the scenario's class, never hidden under not_implemented_mcp
        fault = self.raw_line("T-labels-0004", "mcp", 1, http=502, outcome="fail",
                              error={"code": "backend_fault", "text": "x"},
                              checks=[{"check": "status", "outcome": "fail", "detail": "http 502"}])
        for reps in ({0: m, 1: fault}, {0: fault, 1: m}):
            res = cmp.mcp_agreement(reps, {"direct-local": {0: d}}, s, cases, self.dir)
            self.assertEqual(res["class"], "disagree", res)
            self.assertEqual(res["other_repeats"]["class"], "not_implemented_mcp")
        # the refusal alone stays not_implemented_mcp
        self.assertEqual(cmp.mcp_agreement({0: m}, {"direct-local": {0: d}}, s, cases, self.dir)["class"],
                         "not_implemented_mcp")
        # direct_mismatch "timing": a direct 200 / 503 deadline outside the side stays unverified, not a disagreement
        timing = {"a": dict(base, answer_match=dict(base["answer_match"], direct={"http": 503},
                                                   direct_mismatch="timing"))}
        res = cmp.mcp_agreement({0: m}, {"direct-local": {0: d}}, s, timing, self.dir)
        self.assertEqual((res["class"], res["pair"]["class"]), ("not_implemented_mcp", "service_by_design"))
        res = cmp.mcp_agreement({0: m}, {"direct-local": {0: d400}}, s, timing, self.dir)
        self.assertEqual((res["class"], res["pair"]["class"]), ("not_implemented_mcp", "disagree"))
        # the report lists it among the disagreements, the case and the MCP line's not_implemented_mcp stated, and
        # counts it once (disagree, not not_implemented_mcp)
        res = cmp.mcp_agreement({0: m}, {"direct-local": {0: d400}}, s, cases, self.dir)
        comp = {"order": ["T-labels-0004"], "scenarios": {"T-labels-0004": {"mcp_agreement": res}},
                "summary": {"direct_agreement": {"compared": 0},
                            "mcp_agreement": {"compared": 1, "disagree": 1, "disagreeing": ["T-labels-0004"]}}}
        text = report.agreement_section(comp)
        row = next(x for x in text.splitlines() if x.startswith("| T-labels-0004"))
        self.assertIn("| a | direct |", row)
        self.assertIn("read not_implemented_mcp (case a", row.replace("\\_", "_"))
        self.assertIn("1 disagree", text)
        self.assertNotIn("Not implemented yet via MCP", text)
        # named, the same class; with the default mark the case is service_by_design as before
        named = dict(s, mcp=dict(s["mcp"], by_design=["a"]))
        self.assertEqual(cmp.mcp_agreement({0: m}, {"direct-local": {0: d}}, named, cases, self.dir)["class"],
                         "not_implemented_mcp")
        plain = {"a": dict(base, agreement="by_design", answer_match=dict(base["answer_match"],
                                                                          mcp_outcome="service_by_design"))}
        res = cmp.mcp_agreement({0: dict(m, outcome="service_by_design")}, {"direct-local": {0: d}}, s, plain,
                                self.dir)
        self.assertEqual((res["class"], res["by"], res["verified"]), ("service_by_design", "answer", True))
        # another knob of overrides matches; a field outside overrides.* does not (the checks compare)
        other = dict(err, field="overrides.branching.max_paths")
        m2 = self.raw_line("T-labels-0006", "mcp", 0, http=400, error=other, outcome="not_implemented_mcp",
                           checks=[{"check": "status", "outcome": "fail"}])
        self.assertEqual(cmp.mcp_agreement({0: m2}, {"direct-local": {0: d}}, dict(s, id="T-labels-0006"), cases,
                                           self.dir)["class"], "not_implemented_mcp")
        m3 = self.raw_line("T-labels-0007", "mcp", 0, http=400, error=dict(err, field="preset"), outcome="fail",
                           checks=[{"check": "status", "outcome": "fail"}])
        self.assertEqual(cmp.mcp_agreement({0: m3}, {"direct-local": {0: d}}, dict(s, id="T-labels-0007"), cases,
                                           self.dir)["class"], "disagree")

    def test_planned_duplicate(self):
        # A2-06: T-attempts-0002 runs its attempt per run and arm (run_unique); a repeat that meets its retained id
        # answers 409, which the scenario plans: one arm 200, the other 409 is no disagreement (the recorded run:
        # direct-local 200 pass vs direct-remote 409, its only direct disagreement)
        s = {"id": "T-1", "route": "POST /traverse", "run_unique": ["attempt_id"],
             "expect": [{"check": "status", "args": {"http": [200, 409]}}]}
        st = lambda o: {"check": "status", "outcome": o, "detail": ""}  # noqa: E731
        af = lambda o: {"check": "answer_fields", "outcome": o, "detail": ""}  # noqa: E731
        ok = self.raw_line("T-1", "direct-local", 0, checks=[st("pass"), af("pass"), af("pass")], digest="a")
        dup = self.raw_line("T-1", "direct-remote", 0, http=409, outcome="inconclusive", digest="b",
                            error={"code": "duplicate_attempt", "text": "retained"},
                            checks=[st("pass"), af("inconclusive"), af("pass")])
        res = cmp.direct_agreement({0: ok}, {0: dup}, "POST /traverse", self.dir, None, s)
        self.assertEqual(res["class"], "same_checks", res)
        self.assertIn("direct-remote met its run's retained id", res["planned_duplicates"][0]["detail"])
        # note 3 of the review of lane S5c: with an id per run and arm a split says which arm's state was not planned
        self.assertIn("a split on repeat 0: the id was retained before this arm's first send",
                      res["planned_duplicates"][0]["detail"])
        later = cmp.planned_duplicate(dict(ok, repeat=1), dict(dup, repeat=1), s)
        self.assertIn("a split on repeat 1: direct-local's attempt of an earlier repeat was not retained", later)
        # not planned: no run_unique, no 409 in the status, or the 409's shape check failed
        self.assertEqual(cmp.direct_agreement({0: ok}, {0: dup}, "POST /traverse", self.dir, None,
                                              dict(s, run_unique=[]))["class"], "disagree")
        self.assertEqual(cmp.direct_agreement({0: ok}, {0: dup}, "POST /traverse", self.dir, None,
                                              dict(s, expect=[{"check": "status", "args": {"http": 200}}]))["class"],
                         "disagree")
        bad = dict(dup, outcome="fail", checks=[st("pass"), af("inconclusive"), af("fail")])
        self.assertEqual(cmp.direct_agreement({0: ok}, {0: bad}, "POST /traverse", self.dir, None, s)["class"],
                         "disagree")

    def test_timing_class_traverse(self):
        # a stopped walk's paths are prefixes of the complete walk's
        import checks
        st = lambda o: {"check": "status", "outcome": o, "detail": ""}  # noqa: E731
        def walk(seq, reason):
            arm = {"status": "complete", "segments": [{"id": 0, "parents": [], "from_bp": 0, "length_bp": len(seq),
                                                        "sequence": seq}],
                   "paths": [{"id": 0, "segments": [0], "length_bp": len(seq), "path_reason": reason}]}
            return {"results": [{"outcome": {"walks": "complete"}, "arms": {"right": arm}}]}
        full = "ACGT" * 250
        complete = self.raw_line("T-1", "direct-remote", 0, body=walk(full, "max_extension_bp"),
                                 checks=[st("pass")], digest="b")
        stopped = self.raw_line("T-1", "direct-local", 0, body=walk(full[:978], "time_budget"), outcome="inconclusive",
                                checks=[st("pass")], digest="a", time_limited="seed 0 right arm stopped by time_budget")
        res = cmp.direct_agreement({0: stopped}, {0: complete}, "POST /traverse", self.dir, checks)
        self.assertEqual(res["class"], "timing", res)
        off = self.raw_line("T-1", "direct-local", 1, body=walk("T" + full[1:978], "time_budget"),
                            outcome="inconclusive", checks=[st("pass")], digest="a",
                            time_limited="seed 0 right arm stopped by time_budget")
        res = cmp.direct_agreement({1: off}, {1: dict(complete, repeat=1)}, "POST /traverse", self.dir, checks)
        self.assertEqual(res["class"], "disagree")
        self.assertIn("no prefix", res["outcome_differs"][0]["detail"])

    def test_timing_class_traverse_routes(self):
        # T-labels-0003 of the v4 run: a merged walk displays each merge through ONE parent (SPEC-labeled-
        # traversal-core §7.1), and the complete walk forms merges past the stop depth: a stopped path may continue
        # in the complete answer only through a merge's other parent. The stopped walk's routes are tested against
        # the complete walk's routes through any parent, not against its displayed paths
        import checks
        st = lambda o: {"check": "status", "outcome": o, "detail": ""}  # noqa: E731

        def seg(sid, parents, seq, at):
            return {"id": sid, "parents": parents, "from_bp": at, "length_bp": len(seq), "sequence": seq}

        def walk(segments, paths, side="right"):
            arm = {"status": "complete", "segments": segments,
                   "paths": [{"id": i, "segments": p, "length_bp": 0, "path_reason": "x"} for i, p in enumerate(paths)]}
            return {"results": [{"outcome": {"walks": "complete"}, "arms": {side: arm}}]}
        # complete: ACGT, then GG or TT, merging into CCCC displayed through GG (the first parent)
        full = walk([seg(0, [], "ACGT", 0), seg(1, [0], "GG", 4), seg(2, [0], "TT", 4), seg(3, [1, 2], "CCCC", 6)],
                    [[0, 1, 3]])
        # stopped: the same split, merged after 2 bp, displayed through TT (another arrival order)
        part = walk([seg(0, [], "ACGT", 0), seg(1, [0], "GG", 4), seg(2, [0], "TT", 4), seg(3, [2, 1], "CC", 6)],
                    [[0, 2, 3]])
        tl = "seed 0 right arm stopped by time_budget"
        complete = self.raw_line("T-1", "direct-remote", 0, body=full, checks=[st("pass")], digest="b")
        stopped = self.raw_line("T-1", "direct-local", 0, body=part, outcome="inconclusive", checks=[st("pass")],
                                digest="a", time_limited=tl)
        res = cmp.direct_agreement({0: stopped}, {0: complete}, "POST /traverse", self.dir, checks)
        self.assertEqual(res["class"], "timing", res)
        self.assertIn("every stopped route a prefix of a complete route", res["timing_repeats"][0]["detail"])
        # the left arm reads leaf -> root: the same walk mirrored
        mirror = lambda body: json.loads(json.dumps(body).replace('"right"', '"left"'))  # noqa: E731
        left_full = mirror(walk([seg(0, [], "TGCA", 0), seg(1, [0], "GG", 4), seg(2, [0], "TT", 4),
                                 seg(3, [1, 2], "CCCC", 6)], [[0, 1, 3]]))
        left_part = mirror(walk([seg(0, [], "TGCA", 0), seg(1, [0], "GG", 4), seg(2, [0], "TT", 4),
                                 seg(3, [2, 1], "CC", 6)], [[0, 2, 3]]))
        c2 = self.raw_line("T-1", "direct-remote", 2, body=left_full, checks=[st("pass")], digest="b")
        s2 = self.raw_line("T-1", "direct-local", 2, body=left_part, outcome="inconclusive", checks=[st("pass")],
                           digest="a", time_limited=tl)
        self.assertEqual(cmp.direct_agreement({2: s2}, {2: c2}, "POST /traverse", self.dir, checks)["class"],
                         "timing")
        # a stopped route in no route of the complete walk still disagrees (here behind the merge's second parent)
        bad = walk([seg(0, [], "ACGT", 0), seg(1, [0], "GG", 4), seg(2, [0], "TA", 4), seg(3, [1, 2], "CC", 6)],
                   [[0, 1, 3]])
        off = self.raw_line("T-1", "direct-local", 1, body=bad, outcome="inconclusive", checks=[st("pass")],
                            digest="a", time_limited=tl)
        res = cmp.direct_agreement({1: off}, {1: dict(complete, repeat=1)}, "POST /traverse", self.dir, checks)
        self.assertEqual(res["class"], "disagree", res)
        self.assertIn("seed 0 right: a 6 bp route is no prefix of a complete route",
                      res["outcome_differs"][0]["detail"])
        # the route search's bound: not shown consistent
        self.assertIn("the route search stopped at 2 states",
                      cmp._route_prefix_problem(part["results"][0]["arms"]["right"],
                                                full["results"][0]["arms"]["right"], "right", 2))


class RequestPredicatesTest(unittest.TestCase):
    """A table case's `request` (the machine-readable scenario_match, applied to named and answer-triggered matches
    alike), the table's host_capabilities / service_settings blocks, the references' resolution order (the run's
    recorded live document, else the table's block, else unknown) and the configuration helpers lane G3's gate
    calls."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = self.tmp.name
        self.cases, self.info = cmp.load_by_design(cmp.DEFAULT_BY_DESIGN)
        with open(os.path.join(SUITE, "manifest", "suite_v4.json")) as fh:
            self.v4 = dict((x["id"], x) for x in json.load(fh)["scenarios"])

    def tearDown(self):
        self.tmp.cleanup()

    def match(self, sid, err, args=None, cases=None, evidence=None):
        scen = dict(self.v4[sid])
        if args is not None:
            scen = dict(scen, mcp=dict(scen["mcp"], arguments=args))
        walls = []
        matched, how = cmp.case_matches(cases or self.cases, scen, {"http": 400, "error": err, "body": None}, None,
                                        {"outcome": "fail"}, None, walls, None, evidence)
        return [c["id"] for c in matched], how, walls

    def test_grammar_and_load_validation(self):
        ok = [{"path": "{field}", "field_prefix": "caps.", "equals": 0}, {"path": "patterns[*].protein", "present": True},
              {"field": "overrides.", "not_in": ["bounds.*", "direction"]},
              {"service": "service_settings.advanced_enabled", "equals": False},
              {"path": "max_mismatches", "gt": {"host": "pattern.mismatches.max.protein"}},
              {"any": [[{"path": "preset", "absent": True}], [{"path": "preset", "equals": "walk"}]]}]
        for p in ok:
            self.assertEqual(cmp.request_problem(p), "", p)
        bad = [[], {"path": "x"}, {"path": "x", "equals": 1, "present": True}, {"path": "x", "service": "y",
                                                                                  "equals": 1},
               {"path": "", "present": True}, {"path": "x", "present": False}, {"path": "x", "in": []},
               {"path": "x", "gt": "big"}, {"path": "x", "equals": {"other": "p"}},
               {"path": "x", "field_prefix": "caps.", "equals": 0}, {"field": "{field}", "present": True},
               {"any": []}, {"any": [[]]}, {"any": [[{"path": "x"}]]}, {"path": "x", "when": 1, "present": True}]
        for p in bad:
            self.assertTrue(cmp.request_problem(p), p)
        blocks = self.cases.blocks
        base = {"id": "a", "agreement": "by_design", "service_answer": {"outcome": "refused_at_submit", "code": "x"}}
        trig = dict(base, answer_match={"routes": ["POST /pattern"], "direct": {"http": 200},
                                        "mcp": {"code": "x"}, "direct_mismatch": "disagree"})
        # required: an answer-triggered case without a wall_rule, a case whose field names several fields
        self.assertEqual(cmp.case_request_problems(trig, blocks),
                         ["request is required for an answer-triggered case without a wall_rule"])
        walled = dict(trig, answer_match=dict(trig["answer_match"], wall_rule="budget_minus_reserve"))
        self.assertEqual(cmp.case_request_problems(walled, blocks), [])
        multi = dict(base, service_answer=dict(base["service_answer"], field="caps.a | caps.b"))
        self.assertEqual(cmp.case_request_problems(multi, blocks), [
            "request is required for a case whose field names several request fields (caps.a | caps.b)"])
        self.assertEqual(cmp.case_request_problems(dict(multi, request=[]), blocks),
                         ["request must be a non-empty list of predicates"])
        # every reference must be stated by the table's own block (it decides without a live document)
        ref = dict(trig, request=[{"path": "x", "gt": {"service": "presets.nope.budgets.time_budget_ms"}}])
        self.assertEqual(cmp.case_request_problems(ref, blocks), [
            "request names {'service': 'presets.nope.budgets.time_budget_ms'}, which the table's service_settings "
            "block does not state"])
        self.assertTrue(cmp.case_request_problems(dict(ref, request=[{"service": "service_settings.advanced_enabled",
                                                                      "equals": False}]), {}))
        path = os.path.join(self.dir, "t.json")
        sc.write_json_atomic(path, {"cases": [trig]})
        with self.assertRaises(cmp.TableError):
            cmp.load_by_design(path)
        sc.write_json_atomic(path, {"cases": [dict(trig, request=[{"path": "x", "present": True}])],
                                    "service_settings": []})
        with self.assertRaisesRegex(cmp.TableError, "service_settings must be an object"):
            cmp.load_by_design(path)

    def test_real_table_requests(self):
        # the blocks and the cases that must state a request
        self.assertEqual(self.info["blocks"], ["host_capabilities", "service_settings"])
        b = self.cases.blocks
        self.assertEqual((sc.lookup(b["host"], "pattern.mismatches.max.protein"),
                          sc.lookup(b["host"], "resolve.time_budget.max_time_ms"),
                          sc.lookup(b["service"], "service_settings.advanced_enabled"),
                          sc.lookup(b["service"], "presets.walk.budgets.time_budget_ms"),
                          sc.lookup(b["service"], "service_settings.resolve.http_timeouts_s.phase_margin")),
                         ((True, 1), (True, 60000), (True, False), (True, 15000), (True, 25.0)))
        with_request = sorted(c["id"] for c in self.cases.values() if "request" in c)
        self.assertEqual(with_request, sorted([
            "mismatch.above_kind_cap", "pattern.caps_zero", "pattern.path_caps_without_long_search",
            "pattern.predicate_caps_without_predicate", "pattern.predicate_form_other",
            "pattern.predicate_options_noop_without_predicate", "pattern.predicate_options_without_predicate",
            "resolve.discover_limit_above_cap", "resolve.discover_limit_above_cap_host_answers",
            "resolve.time_budget_above_host_max", "resolve.time_budget_unnamed", "traverse.advanced_path_disabled",
            "traverse.bounds_budget_above_preset_row", "traverse.bounds_knob_above_service_maximum",
            "traverse.empty_seeds"]))
        self.assertEqual([c["id"] for c in self.cases.values() if cmp.request_required(c) and "request" not in c],
                         [])
        # the stale 20000 is gone: the host's maximum is 60000 on staging
        with open(cmp.DEFAULT_BY_DESIGN) as fh:
            raw = fh.read()
        self.assertNotIn("20,000 on staging", raw)
        self.assertNotIn("max_time_ms 20000 on staging", raw)
        # every frozen manifest names only cases its scenarios' own arguments satisfy
        for name in sorted(os.listdir(os.path.join(SUITE, "manifest"))):
            with open(os.path.join(SUITE, "manifest", name)) as fh:
                scenarios = json.load(fh).get("scenarios") or []
            self.assertEqual(cmp.request_scenario_problems(self.cases, scenarios), [], name)

    def test_named_scenario_must_satisfy_its_case(self):
        # a manifest scenario naming a case its own arguments do not satisfy is refused at load (a drift): P-caps-0001
        # with max_anchors 1 still naming pattern.caps_zero; T-end-0001 with a budget inside the walk row
        scenarios = [dict(self.v4["P-caps-0001"], mcp=dict(self.v4["P-caps-0001"]["mcp"], arguments=dict(
            self.v4["P-caps-0001"]["mcp"]["arguments"], max_anchors=1)))]
        problems = cmp.request_scenario_problems(self.cases, scenarios)
        self.assertEqual(len(problems), 1)
        self.assertTrue(problems[0].startswith("scenario P-caps-0001 names case pattern.caps_zero, whose request its "
                                               "own MCP arguments do not satisfy"), problems[0])
        t = self.v4["T-end-0001"]
        inside = dict(t, mcp=dict(t["mcp"], arguments=dict(t["mcp"]["arguments"], overrides={
            "bounds": {"max_extension_bp": 200, "time_budget_ms": 15000}})))
        self.assertEqual(len(cmp.request_scenario_problems(self.cases, [inside])), 1)
        with self.assertRaisesRegex(cmp.TableError, "scenario P-caps-0001 names case pattern.caps_zero"):
            cmp.load_by_design(cmp.DEFAULT_BY_DESIGN, scenarios)
        self.assertEqual(cmp.request_scenario_problems(self.cases, [self.v4[x] for x in (
            "P-caps-0001", "P-modes-0006", "L-caps-0001", "L-caps-0002", "T-end-0001", "T-edge-0012", "T-edge-0016",
            "M-protein-0010", "R-edge-0021", "Q-motif-0004", "Q-motif-0005", "Q-more-0003")]), [])
        # compare's command line refuses such a run's manifest (exit 4), as the runner refuses to start
        sc.write_json_atomic(os.path.join(self.dir, "manifest.json"), {"scenarios": scenarios})
        sc.write_json_atomic(os.path.join(self.dir, "run.json"), {"run_id": "R"})
        open(os.path.join(self.dir, "executions.jsonl"), "w").close()
        self.assertEqual(cmp.main([self.dir, "--no-merge"]), 4)

    def test_mutations_no_longer_swallowed(self):
        # VERIFY-B's finding-2 mutations through the matcher, the real table: each refusal the request did not ask
        # for no longer matches its case (the checks judge), the legitimate ones still do
        caps = {"code": "pattern_invalid", "text": "caps.max_contexts: Input should be greater than or equal to 1",
                "field": "caps.max_contexts"}
        self.assertEqual(self.match("P-exact-0001", caps)[0], [])                      # max_contexts 10000 sent
        self.assertEqual(self.match("P-exact-0001", dict(caps, field="caps.max_paths"))[0], [])   # never sent
        self.assertEqual(self.match("P-caps-0001", dict(caps, field="caps.time_budget_ms"))[0], [])
        self.assertEqual(self.match("P-caps-0001", dict(caps, field="caps.max_contexts"))[0], [])
        named, how, walls = self.match("P-caps-0001", dict(caps, field="caps.max_anchors"))
        self.assertEqual((named, how), (["pattern.caps_zero"], "named"))
        self.assertEqual(self.match("P-modes-0006", caps)[:2], (["pattern.caps_zero"], "named"))
        # the advanced path: a basic override, or a key never sent
        adv = {"code": "preset_unavailable", "field": "overrides.bounds.max_extension_bp",
               "text": "overrides.bounds.max_extension_bp is on the advanced strategy path, which is not enabled on "
                       "this deployment"}
        self.assertEqual(self.match("T-walk-0001", adv)[0], [])
        unsent = dict(adv, field="overrides.labels.max_seed_labels",
                      text=adv["text"].replace("bounds.max_extension_bp", "labels.max_seed_labels"))
        self.assertEqual(self.match("T-walk-0001", unsent)[0], [])
        ev = cmp.RequestEvidence()
        self.assertEqual(self.match("T-labels-0004", unsent, evidence=ev)[:2],
                         (["traverse.advanced_path_disabled"], "answer-triggered"))
        self.assertEqual(ev.sources, {"traverse.advanced_path_disabled": {"service_settings.advanced_enabled":
                                                                          "table"}})
        # a refused field the tool never produces for a basic-only request ("overrides." alone, "overrides.labels"
        # of a labels section): no refusal of the advanced path (T-labels-0003 sends only basic labels keys)
        for field in ("overrides.", "overrides.labels"):
            self.assertEqual(self.match("T-labels-0003", dict(unsent, field=field))[0], [], field)
        # empty seeds for a scenario that sent one; a dna pattern within its cap refused as above the kind's cap
        empty = {"code": "invalid_arguments", "text": "provide `seeds` (a list) or one `sequence` / `sequence_id`"}
        self.assertEqual(self.match("T-walk-0001", empty)[0], [])
        self.assertEqual(self.match("T-edge-0016", empty)[:2], (["traverse.empty_seeds"], "named"))
        # the tool reads an empty sequence as none
        self.assertEqual(self.match("T-edge-0016", empty, args=dict(self.v4["T-edge-0016"]["mcp"]["arguments"],
                                                                    sequence=""))[0], ["traverse.empty_seeds"])
        self.assertEqual(self.match("T-edge-0016", empty, args=dict(self.v4["T-edge-0016"]["mcp"]["arguments"],
                                                                    sequence="ACGT"))[0], [])
        kind = {"code": "pattern_out_of_bounds", "field": "max_mismatches",
                "text": "max_mismatches 2 is above the maximum 1 for protein patterns on the host of database 'x'"}
        self.assertEqual(self.match("M-protein-0010", kind)[:2], (["mismatch.above_kind_cap"], "named"))
        dna = dict(kind, text=kind["text"].replace("protein", "dna"))
        self.assertEqual(self.match("M-planted-0004", dna)[0], [])
        self.assertEqual(self.match("M-planted-0004", kind)[0], [])    # the request has no protein pattern
        # the preset row: a refusal of max_extension_bp, or of a budget the request keeps within the row
        row = {"code": "strategy_out_of_bounds", "field": "bounds.time_budget_ms",
               "text": "bounds.time_budget_ms = 60000 exceeds its maximum 15000 (the walk preset's budget row)"}
        self.assertEqual(self.match("T-end-0001", row)[:2], (["traverse.bounds_budget_above_preset_row"], "named"))
        named, _, walls = self.match("T-end-0001", dict(row, field="bounds.max_extension_bp"))
        self.assertEqual(named, [])
        self.assertEqual(self.match("T-end-0001", row, args=dict(self.v4["T-end-0001"]["mcp"]["arguments"],
                                                                 overrides={"bounds": {"time_budget_ms": 15000}}))[0],
                         [])
        named, _, walls = self.match("T-end-0001", dict(row, field="bounds.max_memory_mb"))
        self.assertEqual(named, [])
        self.assertTrue(walls[0][1].startswith("traverse.bounds_budget_above_preset_row: request predicate "
                                               '{"any": '), walls)
        # the service maximum of bounds.max_live_paths (1000; presets._check_knobs, outside SERVICE_BOUNDED)
        live = {"code": "strategy_out_of_bounds", "field": "bounds.max_live_paths",
                "text": "bounds.max_live_paths = 1001 exceeds its maximum 1000 (the service maximum)"}
        knob_args = dict(self.v4["T-end-0001"]["mcp"]["arguments"], overrides={"bounds": {"max_live_paths": 1001}})
        knob = dict(self.v4["T-end-0001"], mcp=dict(self.v4["T-end-0001"]["mcp"], arguments=knob_args,
                                                    by_design=["traverse.bounds_knob_above_service_maximum"]))
        self.v4["T-knob"] = dict(knob, id="T-knob")
        self.assertEqual(self.match("T-knob", live)[:2], (["traverse.bounds_knob_above_service_maximum"], "named"))
        self.assertEqual(self.match("T-knob", live, args=dict(knob_args, overrides={"bounds": {
            "max_live_paths": 1000}}))[0], [])
        # above the kind's cap with a host whose protein cap is another: the service text states the live cap
        two = dict(kind, text=kind["text"].replace("maximum 1", "maximum 2"))
        self.assertEqual(self.match("M-protein-0010", two)[0], ["mismatch.above_kind_cap"])
        # trie's row (30000) on a trie request
        trie = dict(self.v4["T-end-0001"]["mcp"]["arguments"], preset="trie",
                    overrides={"bounds": {"time_budget_ms": 20000}})
        self.assertEqual(self.match("T-end-0001", row, args=trie)[0], [])
        trie["overrides"] = {"bounds": {"time_budget_ms": 40000}}
        self.assertEqual(self.match("T-end-0001", row, args=trie)[0], ["traverse.bounds_budget_above_preset_row"])

    def test_request_semantics(self):
        refs = sc.References({"service": {"x": {"cap": 5}}, "host": {"k": 2}})
        args = {"a": {"b": 3}, "lst": [{"p": 1}, {"q": 2}], "zero": 0, "none": None, "flag": True}
        scen = {"mcp": {"arguments": args}}

        def holds(pred, field=None):
            return cmp.request_holds({"id": "c", "request": [pred]}, scen, field, refs)[0]
        self.assertTrue(holds({"path": "a.b", "equals": 3}))
        self.assertTrue(holds({"path": "a.b", "lt": {"service": "x.cap"}}))
        self.assertFalse(holds({"path": "a.b", "gt": {"service": "x.cap"}}))
        self.assertFalse(holds({"path": "a.b", "gt": {"service": "x.unknown"}}))   # an unknown reference
        self.assertTrue(holds({"path": "lst[*].q", "present": True}))
        self.assertTrue(holds({"path": "lst[1].q", "equals": 2}))
        self.assertTrue(holds({"path": "none", "absent": True}))                    # null is absent
        self.assertTrue(holds({"path": "missing", "not_in": [1]}))                  # not_in holds on an absent path
        self.assertFalse(holds({"path": "missing", "in": [1]}))
        self.assertFalse(holds({"path": "missing", "equals": None}))
        self.assertFalse(holds({"path": "flag", "equals": 1}))                       # true is not 1
        self.assertTrue(holds({"path": "{field}", "field_prefix": "caps.", "equals": 0}, "caps.zero"))
        self.assertFalse(holds({"path": "{field}", "field_prefix": "caps.", "equals": 0}, "other.zero"))
        self.assertFalse(holds({"path": "{field}", "present": True}, None))
        self.assertTrue(holds({"path": "{field}", "absent": True}, None))
        self.assertTrue(holds({"field": "overrides.", "not_in": ["bounds.*", "direction"]}, "overrides.labels.x"))
        self.assertFalse(holds({"field": "overrides.", "not_in": ["bounds.*"]}, "overrides.bounds.max_steps"))
        self.assertFalse(holds({"field": "overrides.", "not_in": ["bounds.*"]}, "overrides.bounds"))
        self.assertFalse(holds({"field": "overrides.", "not_in": ["bounds.*"]}, "preset"))
        self.assertTrue(holds({"field": "", "in": ["bounds.time_budget_ms"]}, "bounds.time_budget_ms"))
        self.assertTrue(holds({"host": "k", "ge": 2}))
        self.assertFalse(holds({"host": "missing", "absent": True}))                # unknown: false
        self.assertTrue(holds({"any": [[{"path": "zero", "equals": 1}], [{"path": "zero", "equals": 0}]]}))
        ok, why = cmp.request_holds({"id": "c", "request": [{"path": "zero", "equals": 1}]}, scen, "f", refs)
        self.assertEqual((ok, why), (False, 'c: request predicate {"equals": 1, "path": "zero"} does not hold '
                                            '(field f)'))
        self.assertEqual(cmp.request_holds({"id": "c"}, scen, None, refs), (True, ""))

    def test_grammar_pinned(self):
        # the readings two implementations could differ on, pinned (TABLE-DIFF section 2)
        refs = sc.References({"service": {"x": {"cap": 5}}})
        args = {"a": {"b": 3}, "lst": [{"p": 1}, {"q": 2}], "pair": [1, True], "flag": True, "caps": {"z": 0}}
        scen = {"mcp": {"arguments": args}}

        def holds(pred, field=None):
            return cmp.request_holds({"id": "c", "request": [pred]}, scen, field, refs)[0]
        # a refused field equal to its prefix names no key: not "present" (the whole arguments object)
        self.assertFalse(holds({"path": "{field}", "field_prefix": "caps.", "present": True}, "caps."))
        self.assertTrue(holds({"path": "{field}", "field_prefix": "caps.", "absent": True}, "caps."))
        self.assertFalse(holds({"field": "overrides.", "not_in": ["bounds.*"]}, "overrides."))
        # a substituted field that is no path reads absent
        self.assertFalse(holds({"path": "{field}", "present": True}, "a..b"))
        # a list at a path is ONE value; [*] steps into a list's elements, never an object's values
        self.assertTrue(holds({"path": "lst", "equals": [{"p": 1}, {"q": 2}]}))
        self.assertFalse(holds({"path": "lst", "equals": {"p": 1}}))
        self.assertFalse(holds({"path": "lst[*].p", "equals": 2}))
        self.assertTrue(holds({"path": "lst[*].q", "equals": 2}))
        self.assertFalse(holds({"path": "a[*]", "present": True}))
        # JSON equality: true is not 1, also inside lists; numbers by value
        self.assertFalse(holds({"path": "flag", "equals": 1}))
        self.assertTrue(holds({"path": "pair", "equals": [1.0, True]}))
        self.assertFalse(holds({"path": "pair", "equals": [1, 1]}))
        self.assertFalse(holds({"path": "pair", "in": [[True, 1]]}))
        self.assertTrue(sc.json_equal({"a": [1, {"b": 2.0}]}, {"a": [1.0, {"b": 2}]}))
        self.assertFalse(sc.json_equal("1", 1))
        # load: malformed paths, prefixes and references, equals null, null in a list
        for bad in ({"path": "a..b", "present": True}, {"path": "a[x]", "present": True},
                    {"path": "[0]", "present": True}, {"path": "a.{field}x", "present": True},
                    {"path": "a b", "present": True}, {"path": "{field}", "field_prefix": "caps", "equals": 0},
                    {"field": "overrides", "not_in": ["x"]}, {"service": "x[*]", "present": True},
                    {"path": "a", "gt": {"host": "p..q"}}, {"path": "a", "equals": None},
                    {"path": "a", "in": [1, None]}):
            self.assertTrue(cmp.request_problem(bad), bad)
        for good in ({"path": "patterns[*].protein", "present": True}, {"path": "a[0].b", "equals": 1},
                     {"field": "", "in": ["x"]}, {"path": "overrides.{field}", "present": True},
                     {"service": "backends.refseq33m-experimental.resolve.time_budget.default_ms", "gt": 0}):
            self.assertEqual(cmp.request_problem(good), "", good)
        # placeholders match one or more characters, never the empty string
        self.assertFalse(cmp._field_matches("overrides.", "overrides.<knob path>"))
        self.assertTrue(cmp._field_matches("overrides.x", "overrides.<knob path>"))
        self.assertTrue(cmp._field_matches("predicate.at_least.n", "predicate | predicate.<path of the first fault>"))
        self.assertFalse(cmp._field_matches("predicate.", "predicate.<path of the first fault>"))
        # a case both named and answer-triggered whose request holds it back: noted once
        t = dict(self.v4["T-walk-0001"])
        t = dict(t, mcp=dict(t["mcp"], by_design=["traverse.advanced_path_disabled"]))
        self.v4["T-named-adv"] = dict(t, id="T-named-adv")
        adv = {"code": "preset_unavailable", "field": "overrides.bounds.max_extension_bp",
               "text": "overrides.bounds.max_extension_bp is on the advanced strategy path, which is not enabled on "
                       "this deployment"}
        named, how, walls = self.match("T-named-adv", adv)
        self.assertEqual(named, [])
        self.assertEqual([w[2] for w in walls], ["traverse.advanced_path_disabled"])

    def test_reference_resolution_order(self):
        table = {"service": {"service_settings": {"advanced_enabled": False, "resolve": {"timeout_s": 20.0}}},
                 "host": {"pattern": {"finalize_reserve_ms": 250}}}
        refs = sc.References(table)
        self.assertEqual(refs.resolve("service", "service_settings.advanced_enabled"), (False, "table", None))
        self.assertEqual(refs.resolve("service", "presets.walk"), (None, "unknown", None))
        live = sc.References(table, {"service": {"service_settings": {"advanced_enabled": True},
                                                 "presets": {"walk": {"budgets": {"time_budget_ms": 15000}}}}})
        self.assertEqual(live.resolve("service", "service_settings.advanced_enabled"), (True, "live", {
            "namespace": "service", "ref": "service_settings.advanced_enabled", "live": True, "table": False}))
        self.assertEqual(live.resolve("service", "presets.walk.budgets.time_budget_ms"), (15000, "live", None))
        # a live document that lacks the path (a service before the block): the table's value
        self.assertEqual(live.resolve("service", "service_settings.resolve.timeout_s"), (20.0, "table", None))
        self.assertEqual(sc.References(table, {"service": {"service_settings": {"resolve": {"timeout_s": 20}}}}).resolve(
            "service", "service_settings.resolve.timeout_s")[2], None)          # 20 and 20.0 are one value
        # the terms the plausibility rules read, from the default table (the built-in values without it)
        self.assertEqual(sc.term(None, "resolve_reserve_ms"), 250)
        self.assertEqual(sc.term(sc.References(), "timeout_s"), 20.0)
        self.assertEqual(sc.term(sc.References(), "resolve_max_time_ms"), None)
        res = {"mcp": {"arguments": {"labels": ["562"], "database": "refseq33m-experimental"}}}
        self.assertEqual(sc.resolve_budget(res), (60000, "the service's default"))
        self.assertEqual(sc.resolve_budget({"mcp": {"arguments": {"discover": True}}}),
                         (None, "a discovery without a budget"))
        self.assertEqual(sc.resolve_client_timeout_s(res), 70.25)
        self.assertEqual(sc.resolve_client_timeout_s({"mcp": {"arguments": {"time_budget_ms": 500}}}), 20.0)
        self.assertEqual(sc.resolve_client_timeout_s({"mcp": {"arguments": {"discover": True}}}), 120.0)

    def write_run(self, service=None, host=None, lines=()):
        d = self.dir
        os.makedirs(os.path.join(d, "arms", "mcp"), exist_ok=True)
        run = {"run_id": "R", "arm": "mcp"}
        if service is not None:
            run["service_capabilities"] = service
        sc.write_json_atomic(os.path.join(d, "arms", "mcp", "run.json"), run)
        merged = {"run_id": "R"}
        if host is not None:
            merged["capabilities"] = {"/capabilities": host}
        sc.write_json_atomic(os.path.join(d, "run.json"), merged)
        with open(os.path.join(d, "executions.jsonl"), "w") as fh:
            for x in lines:
                fh.write(json.dumps(x) + "\n")

    def test_recorded_documents_and_gate_helpers(self):
        # a run before the runner read the document: the table's blocks, named; the gate's evidence is missing
        self.write_run(lines=[{"scenario": "T-labels-0004", "arm": "mcp", "repeat": 0, "request_evidence": {
            "traverse.advanced_path_disabled": {"service_settings.advanced_enabled": "table"}}}])
        bound = cmp.bind_run(self.cases, self.dir)
        self.assertEqual(cmp.reference_sources(bound)["service"],
                         {"source": "table", "origin": "the table's service_settings block"})
        self.assertEqual(cmp.reference_warnings(bound)[0],
                         "by-design references resolved from: host table (the table's host_capabilities block); "
                         "service table (the table's service_settings block)")
        problems = cmp.gate_configuration_problems(self.dir, self.cases)
        self.assertEqual(problems, [
            "configuration evidence missing: the run recorded no live service_settings (run.json "
            "service_capabilities)",
            "T-labels-0004 r0 case traverse.advanced_path_disabled read service_settings.advanced_enabled: table"])
        # the runner's record {document, sha256, read_at, base}: live, a mismatch judged with the live value
        doc = json.loads(json.dumps(self.cases.blocks["service"]))
        doc["service_settings"]["advanced_enabled"] = True
        self.write_run(service={"document": doc, "sha256": "s", "read_at": "t", "base": "b"},
                       host=self.cases.blocks["host"],
                       lines=[{"scenario": "T-labels-0004", "arm": "mcp", "repeat": 0, "request_evidence": {
                           "traverse.advanced_path_disabled": {"service_settings.advanced_enabled": "live"}}}])
        bound = cmp.bind_run(self.cases, self.dir)
        self.assertEqual(cmp.reference_sources(bound), {
            "host": {"source": "live", "origin": "run.json capabilities[/capabilities]"},
            "service": {"source": "live", "origin": "arms/mcp/run.json service_capabilities"}})
        self.assertEqual(cmp.reference_mismatches(bound), [{
            "namespace": "service", "ref": "service_settings.advanced_enabled", "live": True, "table": False}])
        self.assertIn("configuration mismatch: service service_settings.advanced_enabled is true live, false in the "
                      "table (judged with the live value)", cmp.reference_warnings(bound))
        self.assertEqual(cmp.gate_configuration_problems(self.dir, self.cases), [
            "configuration mismatch: service service_settings.advanced_enabled is true live, false in the table"])
        ev = cmp.RequestEvidence()
        unsent = {"code": "preset_unavailable", "field": "overrides.labels.max_seed_labels",
                  "text": "overrides.labels.max_seed_labels is on the advanced strategy path, which is not enabled on "
                          "this deployment"}
        self.assertEqual(self.match("T-labels-0004", unsent, cases=bound, evidence=ev)[0], [])
        self.assertEqual(ev.mismatches[0]["case"], "traverse.advanced_path_disabled")
        # the same document: no mismatch, no gate problem
        self.write_run(service={"document": self.cases.blocks["service"], "sha256": "s", "read_at": "t"},
                       host=self.cases.blocks["host"])
        self.assertEqual(cmp.gate_configuration_problems(self.dir, self.cases), [])
        # a failed read: the error instead of the document
        self.write_run(service={"error": "HTTP 503"})
        self.assertEqual(cmp.gate_configuration_problems(self.dir, self.cases), [
            "configuration evidence missing: the run recorded no live service_settings (run.json "
            "service_capabilities: HTTP 503)"])


class ReportTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = self.tmp.name
        build_run(self.dir)
        merge_run.merge(self.dir)
        sc.write_json_atomic(os.path.join(self.dir, "compare.json"), cmp.compare(self.dir, (CASES, {"cases": 5})))

    def tearDown(self):
        self.tmp.cleanup()

    def test_escaping(self):
        self.assertEqual(report.ascii_safe("a ≤ b ⊆ c µs 中"), "a <= b subset of c us ?")
        self.assertEqual(report.md("x_y|z*$1"), "x\\_y\\|z\\*\\$1")
        long_seq = "ACGT" * 30
        self.assertIn("(120 bp)", report.ascii_safe(long_seq))
        self.assertTrue(report.md("a" * 50, 10).endswith("..."))
        self.assertEqual(report.table(["a"], []), "(none)\n")
        t = report.table(["area", "n", "detail"], [["P", "12", "a long detail text"], ["R", "-", "x"]])
        sep = t.splitlines()[1]
        self.assertEqual(sep, "|:----|---:|:------------------|")

    def test_report_sections(self):
        text = report.build_report(self.dir)
        for heading in ("# Verdict", "# Environment", "# Outcomes per area", "# Failures",
                        "# Agreement between the arms", "# Timing", "## Cold vs warm", "# Overheads",
                        "# Resource telemetry per area (contextual)"):
            self.assertIn(heading, text)
        text.encode("ascii")  # ASCII-safe throughout
        self.assertIn("INCOMPLETE", text)  # direct-remote was interrupted
        self.assertIn("Load caveat", text)
        self.assertIn("R-001 (direct-local, repeat 0): fail", text)
        self.assertIn('expected: \\{"level": 0\\}', text)
        self.assertIn("- route: POST /pattern; http 200", text)
        self.assertIn("subset of", text)
        self.assertIn("M-002", text)  # listed as not implemented yet via MCP
        self.assertNotIn("checks stub", text.lower())
        self.assertIn("| P | direct-local | cold | 5 |", text)
        # review N9: the server columns state their own n (the error answers carry no host timing)
        self.assertIn("| server n |", text)
        self.assertIn("| P | direct-remote | warm | 4 |", text)
        self.assertIn("warm pairs", text)
        self.assertIn("2 pair(s) verified", text)

    def test_findings_and_inconclusive(self):
        manifest = {"scenarios": [
            {"id": "T-1", "title": "walk within its budget", "route": "POST /traverse",
             "expect": [{"check": "status", "args": {"http": 200}},
                        {"check": "budget_respected", "args": {"max_ms": 15000}}]},
            {"id": "T-2", "route": "POST /traverse", "expect": [{"check": "status", "args": {"http": 200}},
                                                                {"check": "walk_in_record", "args": {}}]}]}
        with tempfile.TemporaryDirectory() as d:
            os.makedirs(os.path.join(d, "arms", "direct-local"))
            sc.write_json_atomic(os.path.join(d, "arms", "direct-local", "manifest.json"), manifest)
            write_arm(d, "direct-local", [
                line("T-1", "direct-local", 0, T0, 34000.0, outcome="finding", checks=[
                    {"check": "status", "outcome": "pass", "detail": "http 200"},
                    {"check": "budget_respected", "outcome": "finding", "detail": "34012 ms on a 15000 ms budget"}]),
                line("T-2", "direct-local", 0, T0 + 40, 10.0, outcome="inconclusive", checks=[
                    {"check": "status", "outcome": "pass", "detail": "http 200"},
                    {"check": "walk_in_record", "outcome": "inconclusive", "detail": "stopped by max_memory"}])])
            merge_run.merge(d)
            comp = cmp.compare(d, ({}, None))
            text = report.build_report(d)
        counts = comp["summary"]["per_area"]["T"]["direct-local"]
        self.assertEqual((counts["pass"], counts["finding"], counts["inconclusive"]), (0, 1, 1))
        self.assertIn("**FINDINGS (1 execution(s)", text)
        self.assertIn("1 inconclusive execution(s) (not counted as pass)", text)
        self.assertIn("# Findings", text)
        self.assertIn("check budget\\_respected: observed 34012 ms on a 15000 ms budget", text)
        self.assertIn('expected: \\{"max\\_ms": 15000\\}', text)

    def test_merge_applies_the_property_rule(self):
        # an MCP runner whose own rule counts an execution with every property check not exposed as a pass:
        # merge_run applies checks.scenario_outcome where it is stricter, keeping the runner's outcome; a run
        # judged by the stub is left as its runner said
        nx = [{"check": "status", "outcome": "pass", "detail": ""}, {"check": "invariants", "outcome": "pass", "detail": ""},
              {"check": "count_relation", "outcome": "not_exposed", "detail": ""}]
        with tempfile.TemporaryDirectory() as d:
            write_arm(d, "mcp", [line("P-1", "mcp", 0, T0, 5.0, checks=nx),
                                 line("P-2", "mcp", 0, T0 + 1, 5.0)], run_extra={"checks_module": {"file": "x", "stub": False}})
            write_arm(d, "direct-local", [line("P-1", "direct-local", 0, T0 + 2, 5.0, checks=nx)],
                      run_extra={"checks_module": "stub"})
            merged = merge_run.merge(d)
            rows, _ = sc.read_jsonl(os.path.join(d, "executions.jsonl"))
        key = dict((sc.execution_key(r), r) for r in rows)
        self.assertEqual((key[("P-1", "mcp", 0)]["outcome"], key[("P-1", "mcp", 0)]["outcome_runner"]),
                         ("inconclusive", "pass"))
        self.assertEqual(key[("P-2", "mcp", 0)]["outcome"], "pass")
        self.assertNotIn("outcome_runner", key[("P-2", "mcp", 0)])
        self.assertEqual(key[("P-1", "direct-local", 0)]["outcome"], "pass")      # the stub's run: untouched
        self.assertEqual(merged["outcomes_rederived"], {"mcp": 1})

    def test_relations(self):
        body = lambda results, ms: {"patterns": [{"id": "w", "results": results, "timing": {"elapsed_ms": ms}}]}  # noqa: E731
        ctx = [{"kmer": "ACGT", "offset": 0, "strand": "+"}, {"kmer": "ACGT", "offset": 5, "strand": "-"}]
        manifest = {"scenarios": [{"id": i, "route": "POST /pattern", "expect": []} for i in ("P-a", "P-b", "P-c")],
                    "relations": [
                        {"id": "rel-same", "type": "same_answer", "a": "P-a", "b": "P-b"},
                        {"id": "rel-sub", "type": "subset", "a": "P-c", "b": "P-a", "args": {"key": ["offset"]}},
                        {"id": "rel-sub-bad", "type": "subset", "a": "P-a", "b": "P-c"},
                        {"id": "rel-prefix", "type": "prefix_equals", "a": "P-a", "b": "P-b"},
                        {"id": "rel-odd", "type": "bogus", "a": "P-a", "b": "P-b"},
                        {"id": "rel-gone", "type": "same_answer", "a": "P-a", "b": "P-z"}]}
        with tempfile.TemporaryDirectory() as d:
            os.makedirs(os.path.join(d, "arms", "direct-local"))
            sc.write_json_atomic(os.path.join(d, "arms", "direct-local", "manifest.json"), manifest)
            ok = [{"check": "status", "outcome": "pass", "detail": ""}]
            write_arm(d, "direct-local", [line(i, "direct-local", 0, T0 + n, 1.0, checks=ok)
                                          for n, i in enumerate(("P-a", "P-b", "P-c"))],
                      raw_bodies={("P-a", 0): body(ctx, 1.0), ("P-b", 0): body(ctx, 9.0),
                                  ("P-c", 0): body(ctx[:1], 2.0)})
            merge_run.merge(d)
            # the built-in evaluation, chosen explicitly (checks.py sits next to compare.py)
            comp = cmp.compare(d, ({}, None), relation_evaluator=None)
            sc.write_json_atomic(os.path.join(d, "compare.json"), comp)
            text = report.build_report(d)
        # the built-in evaluation: same_answer built in, the other types not evaluated
        got = dict((r["id"], (r["arms"].get("direct-local") or {}).get("outcome")) for r in comp["relations"])
        self.assertEqual(got, {"rel-same": "pass", "rel-sub": "not_evaluated", "rel-sub-bad": "not_evaluated",
                               "rel-prefix": "not_evaluated", "rel-odd": "fail", "rel-gone": None})
        self.assertEqual(comp["summary"]["relations"]["per_arm"]["direct-local"],
                         {"pass": 1, "fail": 1, "inconclusive": 0, "not_exposed": 0, "not_evaluated": 3})
        self.assertEqual(comp["summary"]["relations"]["implementation"], "compare.py built-in")
        self.assertIn("# Relations between scenarios", text)
        self.assertIn("| rel-odd | bogus | P-a | P-b | fail |", text)
        self.assertIn("1 failed relation(s)", text)
        # with lane S1's evaluate_relation: its outcome per arm, an unknown outcome read as a failure
        seen = []

        def evaluator(rel, a_sc, a_ans, b_sc, b_ans, arm):
            seen.append((rel["id"], a_sc.get("id"), b_sc.get("id"), arm, a_ans["http"]))
            return {"outcome": {"rel-sub": "inconclusive", "rel-odd": "weird"}.get(rel["id"], "pass"), "detail": "x"}
        with tempfile.TemporaryDirectory() as d2:
            os.makedirs(os.path.join(d2, "arms", "direct-local"))
            sc.write_json_atomic(os.path.join(d2, "arms", "direct-local", "manifest.json"), manifest)
            write_arm(d2, "direct-local", [line(i, "direct-local", 0, T0 + n, 1.0, checks=ok)
                                           for n, i in enumerate(("P-a", "P-b", "P-c"))],
                      raw_bodies={("P-a", 0): body(ctx, 1.0), ("P-b", 0): body(ctx, 9.0),
                                  ("P-c", 0): body(ctx[:1], 2.0)})
            merge_run.merge(d2)
            rows, _ = sc.read_jsonl(os.path.join(d2, "executions.jsonl"))
            grouped, _ = cmp.group_lines(rows)
            rels, summary = cmp.evaluate_relations(sc.read_json(os.path.join(d2, "manifest.json")), grouped, d2,
                                                   ["direct-local"], evaluator=evaluator)
        got = dict((r["id"], (r["arms"].get("direct-local") or {}).get("outcome")) for r in rels)
        self.assertEqual(got["rel-sub"], "inconclusive")
        self.assertEqual(got["rel-odd"], "fail")
        self.assertEqual(summary["implementation"], "checks.evaluate_relation")
        self.assertIn(("rel-same", "P-a", "P-b", "direct-local", 200), seen)
        res = cmp.relation_same_answer({"http": 400, "error": {"code": "c", "text": "t"}},
                                       {"http": 400, "error": {"code": "c", "text": "u"}}, {})
        self.assertEqual(res[0], "fail")
        res = cmp.relation_same_answer({"http": 200, "body": {"a": 1, "b": 2}}, {"http": 200, "body": {"a": 1, "b": 3}},
                                       {"paths": [["a", "a"]]})
        self.assertEqual(res[0], "pass")
        self.assertEqual(cmp.first_difference({"a": [1, {"b": 2}]}, {"a": [1, {"b": 3}]}), "a.1.b: 2 vs 3")
        self.assertEqual(cmp.relation_drop({"a": {"b": 1, "c": 2}, "timing": 1}, ["a.b"]), {"a": {"c": 2}})
        self.assertEqual(cmp.relation_get({"p": [{"r": [1, 2]}, {"r": [3]}]}, "p.0.r.-1"), 2)

    def test_same_answer_built_in_is_the_reference(self):
        # SD-011: compare.py's built-in same_answer (a refused or absent checks.py) reads every case as the reference
        # checks.evaluate_relation does, outcome and detail, on a direct arm and on the MCP arm
        import checks
        from unittest import mock

        def ans(body):
            return {"http": 200, "body": body, "error": None}
        lst = lambda *ts: {"p": [{"t": t, "u": 1} for t in ts]}  # noqa: E731
        cases = [
            ("equal but timing", ans({"x": 1, "timing": {"e": 1}}), ans({"x": 1, "timing": {"e": 2}}), {}),
            ("a value differs", ans({"x": [1, {"y": 2}]}), ans({"x": [1, {"y": 3}]}), {}),
            ("a key only in b", ans({"x": 1}), ans({"x": 1, "y": 2}), {}),
            ("list lengths", ans({"x": [1]}), ans({"x": [1, 2]}), {}),
            ("int vs float", ans({"x": 1}), ans({"x": 1.0}), {}),
            ("bool vs int", ans({"x": True}), ans({"x": 1}), {}),
            ("false vs 0 in a list", ans({"x": [False]}), ans({"x": [0]}), {}),
            ("bool vs bool", ans({"x": [True, False]}), ans({"x": [True, False]}), {}),
            ("bool path", ans({"x": True}), ans({"y": 1}), {"paths": [["x", "y"]]}),
            ("ignore a key", ans({"x": 1, "y": 2}), ans({"x": 1, "y": 3}), {"ignore": ["y"]}),
            ("ignore through a list", ans(lst(1, 2)), ans(lst(3, 4)), {"ignore": ["p.t"]}),
            ("ignore with '*'", ans(lst(1, 2)), ans(lst(3, 4)), {"ignore": ["p.*.t"]}),
            ("ignore with an index", ans(lst(1, 2)), ans(lst(3, 2)), {"ignore": ["p.0.t"]}),
            ("paths equal", ans({"a": {"b": 1}, "c": 2}), ans({"a": {"b": 1}, "c": 3}), {"paths": [["a.b", "a.b"]]}),
            ("paths differ", ans({"a": {"b": 1}}), ans({"a": {"b": 2}}), {"paths": [["a", "a"]]}),
            ("path with '*'", ans({"p": [1, 2]}), ans({"p": [1, 3]}), {"paths": [["p.*", "p.*"]]}),
            ("negative index", ans({"p": [1, 2]}), ans({"p": [3, 2]}), {"paths": [["p.-1", "p.-1"]]}),
            ("whole body against a path", ans({"k": 31}), ans({"pattern": {"k": 31}}), {"paths": [["", "pattern"]]}),
            ("path through timing", ans({"timing": {"a": 1, "timing": 2}}), ans({"timing": {"a": 1, "timing": 3}}),
             {"paths": [["timing", "timing"]]}),
            ("missing path", ans({"p": 1}), ans({"q": 1}), {"paths": [["p", "p"]]}),
            ("200 without a body", ans(None), ans(None), {}),
            ("malformed paths", ans({"p": 1}), ans({"p": 1}), {"paths": ["p"]}),
            ("http differs", ans({"x": 1}), {"http": 400, "error": {"code": "c", "text": "t"}}, {}),
            ("the same refusal", {"http": 400, "error": {"code": "c", "text": "t"}},
             {"http": 400, "error": {"code": "c", "text": "t"}}, {}),
            ("refusals differ", {"http": 400, "error": {"code": "c", "text": "t"}},
             {"http": 400, "error": {"code": "c", "text": "u"}}, {}),
            ("an error that is not an object", {"http": 400, "error": "boom"}, {"http": 400, "error": "boom"}, {})]
        sc_ = {"route": "POST /pattern", "request": {}, "truth": {}}
        for arm in ("direct-local", "mcp"):
            for name, a, b, args in cases:
                rel = {"id": "REL-x", "type": "same_answer", "a": "A", "b": "B", "args": args}
                ref = checks.evaluate_relation(rel, sc_, a, sc_, b, arm)
                manifest = {"scenarios": [dict(sc_, id="A"), dict(sc_, id="B")], "relations": [rel]}
                grouped = {"A": {arm: {0: {"raw": "A", "repeat": 0}}}, "B": {arm: {0: {"raw": "B", "repeat": 0}}}}
                with mock.patch.object(cmp, "load_answer", lambda line, run_dir: {"A": a, "B": b}[line["raw"]]):
                    rels, summary = cmp.evaluate_relations(manifest, grouped, None, [arm], evaluator=None)
                self.assertEqual(summary["implementation"], "compare.py built-in")
                got = rels[0]["arms"][arm]
                self.assertEqual((got["outcome"], got["detail"]), (ref["outcome"], ref["detail"]),
                                 "%s on %s" % (name, arm))
        # the owner's decision on SD-011, in both copies together: a boolean equals only the same boolean (true is
        # not 1), numbers by value (1 equals 1.0), as suite_common.json_equal
        for a, b, want in (({"x": True}, {"x": 1}, "fail"), ({"x": [False]}, {"x": [0]}, "fail"),
                           ({"x": 1}, {"x": 1.0}, "pass"), ({"x": False}, {"x": False}, "pass")):
            self.assertEqual(cmp.relation_same_answer(ans(a), ans(b), {})[0], want, (a, b))
            self.assertEqual(cmp.first_difference(a, b) == "", sc.json_equal(a, b), (a, b))

    def test_inconclusive_verdict(self):
        # review N5: a run in which no execution passed and some were inconclusive is INCONCLUSIVE, never PASS
        manifest = {"scenarios": [{"id": "T-%d" % i, "route": "POST /traverse", "expect": [
            {"check": "walk_in_record", "args": {}}]} for i in (1, 2)]}
        inc = [{"check": "walk_in_record", "outcome": "inconclusive", "detail": "stopped by max_memory"}]
        ok = [{"check": "walk_in_record", "outcome": "pass", "detail": "in record"}]

        def verdict(rows, arms=("direct-local",)):
            with tempfile.TemporaryDirectory() as d:
                for arm in arms:
                    os.makedirs(os.path.join(d, "arms", arm), exist_ok=True)
                    sc.write_json_atomic(os.path.join(d, "arms", arm, "manifest.json"), manifest)
                    write_arm(d, arm, [r for r in rows if r["arm"] == arm])
                merge_run.merge(d)
                cmp.compare(d, ({}, None))
                text = report.build_report(d)
            return [x for x in text.splitlines() if x.startswith("**")][0]

        both = [line("T-1", "direct-local", 0, T0, 10.0, outcome="inconclusive", checks=inc),
                line("T-2", "direct-local", 0, T0 + 1, 10.0, outcome="inconclusive", checks=inc)]
        v = verdict(both)
        self.assertTrue(v.startswith("**INCONCLUSIVE (no execution passed: 2 inconclusive"), v)
        self.assertNotIn("PASS", v)
        self.assertNotIn("not counted as pass", v)
        # with an MCP execution not implemented via MCP beside them: still nothing passed
        nim = line("T-1", "mcp", 0, T0 + 2, 1.0, outcome="not_implemented_mcp", checks=[])
        self.assertTrue(verdict(both + [nim], arms=("direct-local", "mcp")).startswith("**INCONCLUSIVE"))
        mixed = [both[0], line("T-2", "direct-local", 0, T0 + 1, 10.0, checks=ok)]
        v = verdict(mixed)
        self.assertTrue(v.startswith("**PASS ("), v)
        self.assertIn("1 inconclusive execution(s) (not counted as pass)", v)

    def test_unverified_verdict_and_commit(self):
        with tempfile.TemporaryDirectory() as d:
            write_arm(d, "direct-local", [line("P-1", "direct-local", 0, T0, 1.0, checks=[
                {"check": "status", "outcome": "pass", "detail": "ok"}])],
                      run_extra={"checks_module": "stub", "scripts_commit": "0123456789abcdef0123+dirty"})
            merge_run.merge(d)
            text = report.build_report(d)
        self.assertIn("**UNVERIFIED (stub checks", text)
        self.assertNotIn("PASS", text)
        self.assertIn("0123456789ab+dirty", text)
        self.assertEqual(report.short_commit("0123456789abcdef"), "0123456789ab")

    def test_failure_expectation_by_position(self):
        manifest = {"scenarios": [{"id": "Q-1", "expect": [
            {"check": "contains_source_context", "args": {"level": 0}},
            {"check": "contains_source_context", "args": {"level": 2}}]}]}
        rows = [dict(line("Q-1", "direct-local", 0, T0, 1.0, outcome="fail", checks=[
            {"check": "status", "outcome": "pass", "detail": "implicit"},
            {"check": "contains_source_context", "outcome": "pass", "detail": "ok"},
            {"check": "contains_source_context", "outcome": "fail", "detail": "absent at level 2"}]))]
        text = report.failures_section(rows, manifest, 10)
        self.assertIn('expected: \\{"level": 2\\}', text)
        self.assertNotIn('"level": 0', text)

    def test_report_stub_warning_and_quiet_host(self):
        run = sc.read_json(os.path.join(self.dir, "run.json"))
        # the MCP runner records its checks module as {file, stub}
        run["arms"]["mcp"]["checks_module"] = {"file": "/x/checks.py", "stub": True}
        sc.write_json_atomic(os.path.join(self.dir, "run.json"), run)
        os.remove(os.path.join(self.dir, "samples.jsonl"))
        text = report.build_report(self.dir)
        self.assertIn("the stub ran", text)
        self.assertIn("no resource samples", text)

    def test_compare_runs(self):
        with tempfile.TemporaryDirectory() as old:
            build_run(old)
            merge_run.merge(old)
            comp_old = cmp.compare(old, (CASES, {"cases": 5}))
            comp_old["scenarios"]["R-001"]["arms"]["direct-local"]["outcome"] = "pass"
            # review N12: the change compares host-warm executions; the mixed p50 is not used
            comp_old["summary"]["timing"]["P"]["direct-local"]["all"]["p50"] = 1000.0
            comp_old["summary"]["timing"]["P"]["direct-local"]["by_host_cache"]["warm"]["wall"]["p50"] = 4.0
            comp_old["manifest_sha256"] = "different"
            sc.write_json_atomic(os.path.join(old, "compare.json"), comp_old)
            text = report.build_report(self.dir, compare_dir=old)
        self.assertIn("# Change against", text)
        self.assertIn("The manifests differ", text)
        self.assertIn("| R-001 | direct-local | pass | fail |", text)
        change = text.split("# Change against")[1]
        self.assertIn("warm wall p50", change)
        row = [x for x in change.splitlines() if x.startswith("| P | direct-local |")][0]
        self.assertIn("+50%", row)

    def test_compare_json_computed_when_missing(self):
        os.remove(os.path.join(self.dir, "compare.json"))
        report.build_report(self.dir)
        self.assertTrue(os.path.exists(os.path.join(self.dir, "compare.json")))

    @unittest.skipUnless(shutil.which("pandoc") and shutil.which("xelatex"), "pandoc + xelatex not installed")
    def test_pdf(self):
        self.assertEqual(report.main([self.dir]), 0)
        pdf = os.path.join(self.dir, "report.pdf")
        self.assertTrue(os.path.exists(pdf))
        with open(pdf, "rb") as fh:
            self.assertEqual(fh.read(4), b"%PDF")


class PatternGateFixesTest(unittest.TestCase):
    """The /pattern gate's comparison rules (the readiness study of 2026-10-15, round 2 of v6)."""

    @staticmethod
    def body(n, time_limited=False):
        e = {"id": "w", "counts": {"contexts": {"relation": "exact", "value": n, "by_offset": {
            "19": {"relation": "exact", "value": n}}, "by_strand": {"+": {"relation": "exact", "value": n - 1},
                                                                   "-": {"relation": "exact", "value": 1}}},
            "labels": {"relation": "unknown", "value": None}}, "stop": None, "determinism": "full"}
        if time_limited:
            e.update(stop={"phase": "discovery", "reason": "time"}, determinism="time_limited")
        return {"mode": "count", "patterns": [e]}

    def compared(self, local, remote, mcp, whole=False):
        """compare's entry of one /pattern scenario whose arms answered the bodies given ({repeat: body}; None: the
        arm did not run); |whole|: the whole comparison."""
        s = dict(scen("P-010", "P"), expect=[{"check": "status", "args": {"http": 200}},
                                             {"check": "count_relation", "args": {"min": 1}}])
        manifest = {"manifest_version": 1, "content_sha256": sc.scenarios_sha256([s]), "scenarios": [s]}
        checks = [{"check": "status", "outcome": "pass", "detail": "http 200"},
                  {"check": "count_relation", "outcome": "pass", "detail": "contexts exact"}]
        with tempfile.TemporaryDirectory() as d:
            os.makedirs(os.path.join(d, "arms", "direct-local"))
            sc.write_json_atomic(os.path.join(d, "arms", "direct-local", "manifest.json"), manifest)
            for k, (arm, bodies) in enumerate((("direct-local", local), ("direct-remote", remote), ("mcp", mcp))):
                if bodies is None:
                    continue
                write_arm(d, arm, [line("P-010", arm, r, T0 + 10 * k + r, 50.0, checks=checks,
                                        digest="d%s" % bodies[r]["patterns"][0]["counts"]["contexts"]["value"])
                                   for r in sorted(bodies)],
                          raw_bodies=dict((("P-010", r), b) for r, b in bodies.items()))
            merge_run.merge(d)
            comp = cmp.compare(d, ({}, None))
            return comp if whole else comp["scenarios"]["P-010"]

    def test_the_ingress_504_tells_the_direct_outcomes_apart_by_timing(self):
        """VERIFY-V6-2 L1: a direct-local complete pass against direct-remote's front-ingress 504 at its cut (read
        inconclusive) differs by timing, as against the host's 503 deadline; it read "no time limitation tells the
        outcomes apart", a disagreement. A line that records the cut (the runners do, checks.time_limited with the
        scenario) and a line recorded before it (the judging checks' rule) alike."""
        import checks
        s = {"id": "P-011", "route": "POST /pattern", "request": {"patterns": [{"dna": "A" * 25}],
                                                                  "time_budget_ms": 110000},
             "expect": [{"check": "status", "args": {"http": 200}}]}
        ok = line("P-011", "direct-local", 0, T0, 900.0, checks=[{"check": "status", "outcome": "pass", "detail": ""}])
        for recorded in (checks.INGRESS_CUT_LIMITATION, None):
            cut = line("P-011", "direct-remote", 0, T0 + 1, 120100.0, outcome="inconclusive", http=504,
                       digest="e", error={"code": None, "text": "<html>504 Gateway Time-out</html>"},
                       checks=[{"check": "status", "outcome": "inconclusive", "detail": "http 504 from the front"}])
            cut["time_limited"] = recorded
            ok_, why = cmp.timing_repeat(ok, cut, "POST /pattern", None, checks, s)
            self.assertTrue(ok_, why)
            self.assertIn("504", why)
            self.assertEqual(cmp.direct_agreement({0: ok}, {0: cut}, "POST /pattern", None, checks, s)["class"],
                             "timing")
        # without the scenario (or the rule) an unrecorded 504 shows no time limitation, as before
        self.assertFalse(cmp.timing_repeat(ok, cut, "POST /pattern", None, checks)[0])

    def test_two_time_limited_direct_executions_differ_by_timing(self):
        """compare's timing rule when BOTH direct arms state a time limitation (P-heavy-0001 in the v5 baseline: a
        time-limited 200 at_least that passed on direct-local, the host's 503 deadline on direct-remote): neither
        answer is complete, so the pair is timing, not "no time limitation tells the outcomes apart" (a
        disagreement); two time-limited 200s whose counts contradict, or a failing side, stay a disagreement."""
        st = [{"check": "status", "outcome": "pass", "detail": "http 200"},
              {"check": "count_relation", "outcome": "pass", "detail": "contexts at_least 25"}]
        dl = [{"check": "status", "outcome": "inconclusive", "detail": "http 503 deadline"},
              {"check": "count_relation", "outcome": "inconclusive", "detail": "not testable: http 503"}]

        def limited(n, rel="at_least"):
            return {"mode": "count", "patterns": [{"id": "w", "counts": {"contexts": {"relation": rel, "value": n}},
                                                   "stop": {"phase": "discovery", "reason": "time"},
                                                   "determinism": "time_limited"}]}
        with tempfile.TemporaryDirectory() as d:
            os.makedirs(os.path.join(d, "raw"))

            def put(ln, body):
                sc.write_gz_json(os.path.join(d, ln["raw"]), {"answer": {"http": ln["http"], "body": body,
                                                                         "error": ln.get("error")}})
                return ln
            local = put(dict(line("P-012", "direct-local", 0, T0, 59900.0, checks=st),
                             time_limited="pattern 0 time-limited"), limited(25))
            remote = put(dict(line("P-012", "direct-remote", 0, T0 + 61, 60200.0, outcome="inconclusive", http=503,
                                   digest="e", error={"code": "deadline", "text": "deadline"}, checks=dl),
                              time_limited="a 503 deadline: nothing written within the time budget"), None)
            ok, why = cmp.timing_repeat(local, remote, "POST /pattern", d)
            self.assertTrue(ok, why)
            self.assertIn("both time-limited", why)
            self.assertEqual(cmp.direct_agreement({0: local}, {0: remote}, "POST /pattern", d)["class"], "timing")
            # two time-limited 200s: their counts must not contradict (at_least 25 against exact 5)
            remote2 = put(dict(line("P-012", "direct-remote", 1, T0 + 62, 59800.0, outcome="inconclusive",
                                    checks=[st[0], dict(st[1], outcome="inconclusive")]),
                               time_limited="pattern 0 time-limited"), limited(5, "exact"))
            ok, why = cmp.timing_repeat(local, remote2, "POST /pattern", d)
            self.assertFalse(ok)
            self.assertIn("no common value", why)
            # a failing side is no timing
            bad = dict(local, outcome="fail", checks=[st[0], dict(st[1], outcome="fail")])
            self.assertFalse(cmp.timing_repeat(bad, remote, "POST /pattern", d)[0])

    def test_an_exact_count_that_differs_between_access_paths_is_a_disagreement(self):
        """GPT's "doubled count" (the readiness study's DBL probe): a count inflated with all its partitions passes
        the checks on any arm, which judge one answer at a time; an exact count is the population's size on one
        index whichever arm asked, so one that differs between the access paths is a disagreement (0 false alarms
        over the v5 baseline's 3,047 pairs)."""
        same = self.compared({0: self.body(22158)}, {0: self.body(22158)}, {0: self.body(22158), 1: self.body(22158)})
        self.assertEqual(same["mcp_agreement"]["class"], "agree")
        self.assertEqual(same["direct_agreement"]["class"], "identical")
        doubled = self.compared({0: self.body(22158)}, {0: self.body(22158)},
                                {0: self.body(22158), 1: self.body(44316)})
        ma = doubled["mcp_agreement"]
        self.assertEqual((ma["class"], ma["side"], ma["by"]), ("disagree", "mcp", "exact_counts"))
        self.assertEqual(sorted(d["path"] for d in ma["exact_counts_differ"]
                                if d["a"] == ["direct-local", 0] and d["b"] == ["mcp", 1]),
                         ["p0.contexts", "p0.contexts.by_offset.19", "p0.contexts.by_strand.+"])
        self.assertIn("p0.contexts exact 22158 on direct-local r0, exact 44316 on mcp r1", ma["detail"])
        self.assertEqual(doubled["direct_agreement"]["class"], "identical")
        remote = self.compared({0: self.body(22158)}, {0: self.body(22159)}, {0: self.body(22158), 1: self.body(22158)})
        self.assertEqual(remote["direct_agreement"]["class"], "disagree")
        self.assertEqual(remote["direct_agreement"]["exact_counts_differ"][0]["path"], "p0.contexts")
        self.assertEqual(remote["mcp_agreement"]["class"], "disagree")   # it differs from direct-remote's
        # an answer a time budget limited is no population count to compare
        limited = self.compared({0: self.body(22158)}, {0: self.body(22158)},
                                {0: self.body(5, time_limited=True), 1: self.body(22158)})
        self.assertEqual(limited["mcp_agreement"]["class"], "agree")
        self.assertEqual(cmp.exact_counts({"body": self.body(3)}),
                         {"p0.contexts": 3, "p0.contexts.by_offset.19": 3, "p0.contexts.by_strand.+": 2,
                          "p0.contexts.by_strand.-": 1})

    def test_an_exact_count_disagreement_names_its_counts_and_the_odd_arm(self):
        """Round 2's verification of the checks (M1, L2): an exact-count disagreement leaves both outcomes pass, so
        the report's Agreement rows read "outcome differs" (direct) and "-" (MCP); they now state the paths and the
        values. compare blames the MCP only when its count is the odd one: a direct arm's doubled count is the direct
        side's, and with two arms alone (an rtp gate run: direct-local and the MCP) the side is undetermined. The
        scenario disagrees whichever side it is."""
        n, d2 = 22158, 44316
        mcp_odd = self.compared({0: self.body(n)}, {0: self.body(n)}, {0: self.body(n), 1: self.body(d2)}, whole=True)
        ma = mcp_odd["scenarios"]["P-010"]["mcp_agreement"]
        self.assertEqual((ma["side"], ma["exact_counts_arms"], ma["mcp_detail"]), ("mcp", ["mcp"], ma["detail"]))
        remote_odd = self.compared({0: self.body(n)}, {0: self.body(d2)}, {0: self.body(n), 1: self.body(n)},
                                   whole=True)
        ma = remote_odd["scenarios"]["P-010"]["mcp_agreement"]
        self.assertEqual((ma["class"], ma["side"], ma["by"], ma["exact_counts_arms"]),
                         ("disagree", "direct", "exact_counts", ["direct-remote"]))
        self.assertIn("(the odd arm: direct-remote; the MCP states another direct arm's value)", ma["direct_detail"])
        local_odd = self.compared({0: self.body(d2)}, {0: self.body(n)}, {0: self.body(n), 1: self.body(n)})
        self.assertEqual((local_odd["mcp_agreement"]["side"], local_odd["mcp_agreement"]["exact_counts_arms"]),
                         ("direct", ["direct-local"]))
        two = self.compared({0: self.body(d2)}, None, {0: self.body(n), 1: self.body(n)})["mcp_agreement"]
        self.assertEqual((two["class"], two["side"], two["exact_counts_arms"]),
                         ("disagree", "undetermined", ["direct-local", "mcp"]))
        self.assertIn("which arm is wrong is undetermined: direct-local, mcp", two["detail"])
        self.assertFalse(set(two) & {"direct_detail", "mcp_detail"})
        # the report's Agreement rows state the counts: the direct table and the MCP table alike
        text = report.agreement_section(remote_odd)
        rows = [x for x in text.split("\n") if x.startswith("| P-010 |")]
        self.assertEqual(len(rows), 2, text)
        self.assertIn("exact counts differ: p0.contexts exact 22158 on direct-local r0, exact 44316 on direct-remote "
                      "r0", rows[0])
        self.assertNotIn("outcome differs", rows[0])
        self.assertIn("| direct | exact counts differ (the odd arm: direct-remote; the MCP states another direct arm's "
                      "value): p0.contexts exact 44316 on direct-remote r0, exact 22158 on mcp r0", rows[1])
        # a comparison of an earlier compare.py (side mcp, the explanation in `detail` only) still shows it
        ma = mcp_odd["scenarios"]["P-010"]["mcp_agreement"]
        ma.pop("mcp_detail")
        row = [x for x in report.agreement_section(mcp_odd).split("\n") if x.startswith("| P-010 |")][-1]
        self.assertIn("| mcp | exact counts differ: p0.contexts exact 22158 on direct-local r0, exact 44316 on mcp r1",
                      row)


if __name__ == "__main__":
    unittest.main()
