"""Offline tests of merge_run.py, compare.py and report.py on synthetic run directories: the merge of the arms,
outcome counts, direct agreement, the MCP outcome classes with async's by-design table, overheads, cold vs
warm, the resource join, the report's sections and escaping, --compare, and (when pandoc and xelatex are
installed) the PDF."""

import json
import os
import shutil
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
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
           line("M-001", "mcp", 0, tm + 12, 800.0, error={"code": "pattern_unsupported", "text": "no"}),
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
        self.assertEqual(cls["P-003"], "service_by_design")  # answered; host text matched with digits ignored
        self.assertEqual(cls["P-004"], "service_by_design")  # refused_at_submit, code/text/field match
        self.assertEqual(cls["M-001"], "disagree")       # wrong service code
        self.assertEqual(cls["P-005"], "not_implemented_mcp")
        self.assertEqual(cls["M-002"], "not_implemented_mcp")  # by the manifest, no MCP line
        self.assertEqual(cls["S-001"], "service_by_design")    # service_only
        m1 = self.comp["scenarios"]["M-001"]["mcp_agreement"]
        self.assertEqual(m1["side"], "mcp")
        self.assertTrue(m1["direct_verified"])  # the slot code of the raw answer matched
        p3 = self.comp["scenarios"]["P-003"]["mcp_agreement"]
        self.assertTrue(p3["direct_verified"])
        self.assertEqual(self.s["mcp_agreement"]["disagree"], 2)
        # review M2: checks compared per comparison, by-design pairs verified on every side the case has
        self.assertEqual(self.comp["scenarios"]["R-001"]["mcp_agreement"]["checks_compared"], 1)
        self.assertEqual(self.comp["scenarios"]["R-001"]["mcp_agreement"]["checks_not_exposed"], 1)
        s1 = self.comp["scenarios"]["S-001"]["mcp_agreement"]
        self.assertTrue(s1["mcp_verified"] and s1["verified"])
        self.assertEqual((self.s["mcp_agreement"]["by_design_verified"],
                          self.s["mcp_agreement"]["by_design_unverified"]), (3, 0))

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
        unverified = cmp.by_design_class(case, d, [line("R-1", "mcp", 0, T0, 1.0, outcome="service_by_design")],
                                         None)
        self.assertEqual(unverified["class"], "service_by_design")
        self.assertIsNone(unverified["mcp_verified"])
        self.assertFalse(unverified["verified"])

    def test_load_real_table(self):
        cases, info = cmp.load_by_design(cmp.DEFAULT_BY_DESIGN)
        self.assertEqual(info["cases"], 124)
        self.assertTrue(info["service_commit"].startswith("22d034ba"))
        self.assertEqual(cases["pattern.unknown_field"]["agreement"], "not_implemented_mcp")
        self.assertEqual(cmp.load_by_design(""), ({}, None))
        for case in cases.values():
            self.assertIn(case["agreement"], ("by_design", "not_implemented_mcp", "service_only"))
            self.assertIn(case["service_answer"]["outcome"],
                          ("refused_at_submit", "answered", "clamped", "not_expressible", "relayed"))

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
                        "# Resource use per area"):
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
        self.assertIn("3 pair(s) verified", text)

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
        self.assertEqual(cmp.first_difference({"a": [1, {"b": 2}]}, {"a": [1, {"b": 3}]}), "a.1.b")
        self.assertEqual(cmp.drop_path({"a": {"b": 1, "c": 2}}, "a.b"), {"a": {"c": 2}})
        self.assertEqual(cmp.path_values({"p": [{"r": [1, 2]}, {"r": [3]}]}, "p.*.r.*"), [1, 2, 3])

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


if __name__ == "__main__":
    unittest.main()
