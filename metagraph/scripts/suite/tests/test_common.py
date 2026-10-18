"""Offline tests of suite_common (answers, digests, plans, statistics, files), checks_stub and the sampler's
parsing. Run from metagraph/scripts/suite: python3 -m unittest discover -s tests"""

import contextlib
import gzip
import io
import json
import os
import sys
import subprocess
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import checks_stub  # noqa: E402
import sampler  # noqa: E402
import suite_common as sc  # noqa: E402


class TimeAndNames(unittest.TestCase):
    def test_iso_round_trip(self):
        t = 1791800000.123
        s = sc.iso_utc(t)
        self.assertTrue(s.endswith("Z"))
        self.assertAlmostEqual(sc.parse_utc(s), t, places=3)

    def test_parse_utc_rejects(self):
        self.assertIsNone(sc.parse_utc("yesterday"))
        self.assertIsNone(sc.parse_utc(None))

    def test_release_short(self):
        self.assertEqual(sc.release_short("refseq97-b7214d29a934ff89512ad254bbb66ba793842928"), "refseq97-b7214d29")
        self.assertEqual(sc.release_short(None), "unknown")
        self.assertEqual(sc.release_short("a/b c"), "a_b_c")

    def test_parse_route(self):
        self.assertEqual(sc.parse_route("POST /pattern"), ("POST", "/pattern"))
        self.assertEqual(sc.parse_route("get /capabilities"), ("GET", "/capabilities"))
        self.assertEqual(sc.parse_route("/traverse/capabilities"), ("GET", "/traverse/capabilities"))

    def test_raw_name(self):
        self.assertEqual(sc.raw_name("M-x/1", "direct-local", 2), "raw/M-x_1.direct-local.2.json.gz")


class ManifestAndPlan(unittest.TestCase):
    def manifest(self, n=25):
        scen = [{"id": "P-%03d" % i, "area": "P" if i % 2 else "R", "repeat": 1, "timing_class": "fast"}
                for i in range(n)]
        scen[3]["repeat"] = 3
        scen[5]["timing_class"] = "long"
        return {"scenarios": scen, "content_sha256": sc.scenarios_sha256(scen)}

    def test_sha_verifies_and_detects_change(self):
        m = self.manifest()
        self.assertTrue(sc.verify_manifest(m)[0])
        m["scenarios"][0]["title"] = "changed"
        self.assertFalse(sc.verify_manifest(m)[0])

    def test_select_keeps_global_index(self):
        m = self.manifest()
        sel = sc.select_scenarios(m, areas=["P"], limit=3)
        self.assertEqual([i for i, _ in sel], [1, 3, 5])
        sel = sc.select_scenarios(m, ids=["P-010", "P-002"])
        self.assertEqual([s["id"] for _, s in sel], ["P-002", "P-010"])
        # SD-369: a hand-typed list ('P-010, P-002') selects as the clean one, ids stripped as the areas are
        sel = sc.select_scenarios(m, ids=["P-010", " P-002 "])
        self.assertEqual([s["id"] for _, s in sel], ["P-002", "P-010"])
        self.assertEqual([i for i, _ in sc.select_scenarios(m, areas=[" P"], limit=3)], [1, 3, 5])

    def test_build_manifest_areas_list(self):
        # SD-369: build_manifest --areas 'P, T' builds both areas (it dropped T without a word), empties left out
        import build_manifest as bm
        self.assertEqual(bm.area_list("P, T"), {"P", "T"})
        self.assertEqual(bm.area_list(" P,,T ,"), {"P", "T"})
        self.assertIsNone(bm.area_list(None))
        self.assertIsNone(bm.area_list(""))

    def test_repeat_count(self):
        self.assertEqual(sc.repeat_count(0, {"repeat": 1}, 10), 2)
        self.assertEqual(sc.repeat_count(10, {}, 10), 2)
        self.assertEqual(sc.repeat_count(7, {"repeat": 1}, 10), 1)
        self.assertEqual(sc.repeat_count(0, {"repeat": 1}, 0), 1)
        self.assertEqual(sc.repeat_count(4, {"repeat": 3}, 10), 3)
        self.assertEqual(sc.repeat_count(4, {"repeat": "x"}, 10), 1)

    def test_plan_local_and_remote(self):
        m = self.manifest(12)
        sel = sc.select_scenarios(m)
        planned, skipped = sc.plan_executions(sel, "direct-local", 10)
        # 12 scenarios + a warm repeat for indices 0 and 10 + 2 more repeats of index 3
        self.assertEqual(len(planned), 12 + 2 + 2)
        self.assertEqual(skipped, [])
        self.assertEqual([p[3] for p in planned[:2]], ["cold", "warm"])
        planned_r, skipped_r = sc.plan_executions(sel, "direct-remote", 10)
        self.assertEqual(len(planned_r), len(planned) - 1)
        self.assertEqual(skipped_r[0]["scenario"], "P-005")


class Answers(unittest.TestCase):
    def test_decode_body(self):
        data = b'{"a": 1}'
        self.assertEqual(sc.decode_body(gzip.compress(data), "gzip"), data)
        import zlib
        self.assertEqual(sc.decode_body(zlib.compress(data), "deflate"), data)
        raw = zlib.compressobj(wbits=-15)
        deflated = raw.compress(data) + raw.flush()
        self.assertEqual(sc.decode_body(deflated, "deflate"), data)
        self.assertEqual(sc.decode_body(data, None), data)
        with self.assertRaises(ValueError):
            sc.decode_body(data, "br")

    def test_normalise_ok(self):
        a = sc.normalise_direct(200, b'{"x": 1}', {"Content-Type": "application/json", "Server": "nginx"})
        self.assertEqual(a["body"], {"x": 1})
        self.assertIsNone(a["error"])
        self.assertEqual(a["headers"], {"content-type": "application/json"})

    def test_normalise_error_body(self):
        a = sc.normalise_direct(400, b'{"code": "invalid_request", "error": "request: unknown field \'bogus\'"}')
        self.assertEqual(a["error"], {"code": "invalid_request", "text": "request: unknown field 'bogus'"})
        a = sc.normalise_direct(503, b'{"error": "Server is currently initializing, please come back later."}')
        self.assertIsNone(a["error"]["code"])
        a = sc.normalise_direct(502, b"<html>bad gateway</html>")
        self.assertEqual(a["body"], "<html>bad gateway</html>")
        self.assertEqual(a["error"], {"code": None, "text": "<html>bad gateway</html>"})
        a = sc.normalise_direct(500, b"")
        self.assertEqual(a["error"], {"code": None, "text": None})

    def test_transport_failure(self):
        a = sc.transport_failure("timeout", "no answer")
        self.assertEqual(a["http"], 0)
        self.assertEqual(a["error"]["code"], "timeout")

    def test_digest_ignores_timing_and_headers(self):
        a = sc.normalise_direct(200, b'{"patterns": [{"id": "a", "timing": {"elapsed_ms": 1}}], "timing": {"x": 2}}',
                                {"Content-Encoding": "gzip"})
        b = sc.normalise_direct(200, b'{"timing": {"x": 9}, "patterns": [{"timing": {"elapsed_ms": 5}, "id": "a"}]}')
        self.assertEqual(sc.answer_digest(a, "POST /pattern"), sc.answer_digest(b, "POST /pattern"))
        c = sc.normalise_direct(200, b'{"patterns": [{"id": "b"}]}')
        self.assertNotEqual(sc.answer_digest(a, "POST /pattern"), sc.answer_digest(c, "POST /pattern"))

    def test_digest_capability_volatile_fields(self):
        a = sc.normalise_direct(200, b'{"deadline_check": {"observed_max_uninterruptible_ms": 1, "rule": "r"}}')
        b = sc.normalise_direct(200, b'{"deadline_check": {"observed_max_uninterruptible_ms": 7, "rule": "r"}}')
        self.assertEqual(sc.answer_digest(a, "GET /capabilities"), sc.answer_digest(b, "GET /capabilities"))
        self.assertEqual(sc.answer_digest(a, "GET /traverse/capabilities"),
                         sc.answer_digest(b, "GET /traverse/capabilities"))
        # elsewhere an observed_ field is part of the answer
        self.assertNotEqual(sc.answer_digest(a, "POST /traverse"), sc.answer_digest(b, "POST /traverse"))

    def test_extract_server_pattern(self):
        body = {"timing": {"elapsed_ms": 4.0}, "patterns": [
            {"timing": {"elapsed_ms": 1.0}, "work": {"steps": 10, "memory_bytes": 100,
                                                     "steps_by_mismatches": {"0": 4, "1": 6}}},
            {"timing": {"elapsed_ms": 2.0}, "work": {"steps": 5, "memory_bytes": 300,
                                                     "steps_by_mismatches": {"0": 5}}},
            {"error": {"code": "bad_alphabet"}}]}
        s = sc.extract_server(body)
        self.assertEqual(s["elapsed_ms"], 4.0)
        self.assertEqual(s["work"]["steps"], 15)
        self.assertEqual(s["work"]["memory_bytes"], 300)
        self.assertEqual(s["work"]["steps_by_mismatches.0"], 9)
        self.assertEqual(s["work"]["steps_by_mismatches.1"], 6)

    def test_extract_server_entries_only(self):
        body = {"patterns": [{"timing": {"elapsed_ms": 1.5}}, {"timing": {"elapsed_ms": 2.5}}]}
        s = sc.extract_server(body)
        self.assertEqual(s["elapsed_ms"], 4.0)
        self.assertEqual(len(s["timing"]["entries"]), 2)

    def test_extract_server_traverse(self):
        body = {"results": [{"arms": {"right": {"counters": {"steps": 3, "output_bp": 10}},
                                      "left": {"counters": {"steps": 4}}},
                             "annotation": {"rows_requested": 2, "access_path": "x"},
                             "timing": {"walk_ms": 1}}]}
        s = sc.extract_server(body)
        self.assertIsNone(s["elapsed_ms"])
        self.assertEqual(s["work"]["counters.steps"], 7)
        self.assertEqual(s["work"]["counters.output_bp"], 10)
        self.assertEqual(s["work"]["annotation.rows_requested"], 2)

    def test_extract_server_none(self):
        self.assertIsNone(sc.extract_server("text"))
        self.assertIsNone(sc.extract_server({"code": "x", "error": "y"}))

    def test_identity(self):
        docs = {"/capabilities": {"release": "r1", "server_instance": "s1", "pattern": {"index_fp": None}},
                "/pattern/capabilities": {"index": {"index_fp": "fp1"}}}
        self.assertEqual(sc.identity_from_capabilities(docs),
                         {"release": "r1", "server_instance": "s1", "index_fp": "fp1"})

    def test_outcome_from_checks(self):
        self.assertEqual(sc.outcome_from_checks([{"outcome": "pass"}, {"outcome": "not_exposed"}]), "pass")
        self.assertEqual(sc.outcome_from_checks([{"outcome": "pass"}, {"outcome": "fail"}]), "fail")
        self.assertEqual(sc.outcome_from_checks([]), "inconclusive")  # nothing could test anything


class Statistics(unittest.TestCase):
    def test_percentile(self):
        vals = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        self.assertEqual(sc.percentile(vals, 50), 5.5)
        self.assertAlmostEqual(sc.percentile(vals, 90), 9.1)
        self.assertEqual(sc.percentile(vals, 100), 10)
        self.assertEqual(sc.percentile(vals, 0), 1)
        self.assertEqual(sc.percentile([7], 90), 7.0)
        self.assertIsNone(sc.percentile([], 50))
        self.assertEqual(sc.percentile([None, 3, True, 1], 50), 2.0)

    def test_summarise(self):
        s = sc.summarise([3.0, 1.0, 2.0])
        self.assertEqual((s["n"], s["p50"], s["max"], s["mean"]), (3, 2.0, 3.0, 2.0))
        self.assertEqual(sc.summarise([None])["n"], 0)


class Files(unittest.TestCase):
    def test_jsonl_tolerates_torn_line(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "x.jsonl")
            sc.append_jsonl(p, {"a": 1})
            sc.append_jsonl(p, {"a": 2})
            with open(p, "a") as fh:
                fh.write('{"a": 3')
            rows, bad = sc.read_jsonl(p)
            self.assertEqual([r["a"] for r in rows], [1, 2])
            self.assertEqual(bad, 1)
            self.assertEqual(sc.read_jsonl(os.path.join(d, "missing")), ([], 0))

    def test_gz_and_atomic(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "a.json.gz")
            sc.write_gz_json(p, {"x": [1, 2]})
            self.assertEqual(sc.read_gz_json(p), {"x": [1, 2]})
            q = os.path.join(d, "b.json")
            sc.write_json_atomic(q, {"y": 1})
            self.assertEqual(sc.read_json(q), {"y": 1})
            self.assertEqual(sc.read_json(os.path.join(d, "none"), {}), {})

    def test_cli_manifest_sha(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "m.json")
            scen = [{"id": "a"}]
            with open(p, "w") as fh:
                json.dump({"scenarios": scen, "content_sha256": "old"}, fh)
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(sc._main(["manifest-sha", p, "--write"]), 0)
            self.assertEqual(sc.read_json(p)["content_sha256"], sc.scenarios_sha256(scen))


T0 = 1791800000.0


def hline(sid, arm, repeat, t, wall=100.0, outcome="pass", http=200, error=None, server=None, phases=None):
    return {"scenario": sid, "arm": arm, "repeat": repeat, "started_utc": sc.iso_utc(t),
            "ended_utc": sc.iso_utc(t + wall / 1000.0), "wall_ms": wall, "outcome": outcome, "http": http,
            "error": error, "server": server, "mcp_phases": phases}


BACKEND = {"submit_ms": 200.0, "queue_ms": 140.0, "backend_ms": 30.0, "terminal_lag_ms": 2700.0, "paging_ms": 0.0}
SUBMIT_ONLY = {"submit_ms": 200.0, "queue_ms": None, "backend_ms": None, "terminal_lag_ms": None, "paging_ms": None}


class HostCache(unittest.TestCase):
    def test_cold_warm_aged_and_not_reached(self):
        rows = [hline("A", "direct-remote", 0, T0 + 70.0),            # 69.9 s after local: aged
                hline("A", "direct-local", 0, T0),                     # first on the host: cold
                hline("A", "direct-local", 1, T0 + 1.0),               # 0.9 s later: warm
                hline("A", "mcp", 0, T0 + 70.2, outcome="service_by_design"),   # refused by the service
                hline("A", "mcp", 1, T0 + 80.0, phases=BACKEND),       # after the remote one: warm
                hline("B", "mcp", 0, T0 + 2.0, outcome="not_implemented_mcp"),  # never reached the host
                hline("B", "direct-local", 0, T0 + 3.0),               # so this is B's cold execution
                hline("C", "direct-local", 0, T0 + 4.0, http=0, error={"code": "transport", "text": "refused"}),
                hline("C", "direct-remote", 0, T0 + 5.0),
                {"scenario": "D", "arm": "direct-local", "repeat": 0}]
        sc.assign_host_cache(rows)
        got = [(r["scenario"], r["arm"], r["repeat"], r["host_cache"], r["since_prev_s"], r["host_cache_after"])
               for r in rows]
        self.assertEqual(got, [
            ("A", "direct-remote", 0, "aged", 68.9, "direct-local"),
            ("A", "direct-local", 0, "cold", None, None),
            ("A", "direct-local", 1, "warm", 0.9, "direct-local"),
            ("A", "mcp", 0, "not_reached", 0.1, "direct-remote"),
            ("A", "mcp", 1, "warm", 9.9, "direct-remote"),
            ("B", "mcp", 0, "not_reached", None, None),
            ("B", "direct-local", 0, "cold", None, None),
            ("C", "direct-local", 0, "not_reached", None, None),
            ("C", "direct-remote", 0, "cold", None, None),
            ("D", "direct-local", 0, "unknown", None, None)])

    def test_mcp_reached_only_on_backend_evidence(self):
        # review N2: an MCP line reached the host only when it carries the host's timing or the backend phase,
        # whatever its outcome
        rows = [hline("X", "mcp", 0, T0, outcome="fail", http=400, phases=SUBMIT_ONLY,
                      error={"code": "traversal_refused", "text": "host overrun"}),     # refused at submit
                hline("X", "direct-local", 0, T0 + 1.7),                                # so this one is cold
                hline("Y", "mcp", 0, T0, outcome="error", phases=BACKEND,
                      error={"code": "runner_error", "text": "not terminal in time"}),  # the backend ran
                hline("Y", "direct-local", 0, T0 + 1.7),
                hline("Z", "mcp", 0, T0, outcome="service_by_design", http=400, phases=BACKEND,
                      error={"code": "resolve_too_short", "text": "relayed"}),          # a relayed refusal
                hline("Z", "direct-local", 0, T0 + 1.7),
                hline("W", "mcp", 0, T0, outcome="pass", server={"elapsed_ms": 0.7}),   # host timing only
                hline("W", "direct-local", 0, T0 + 1.7),
                hline("V", "mcp", 0, T0, outcome="pass"),                               # no evidence at all
                hline("V", "direct-local", 0, T0 + 1.7),
                hline("U", "mcp", 0, T0, outcome="pass", server={"elapsed_ms": True}, phases={"backend_ms": "30"}),
                hline("U", "direct-local", 0, T0 + 1.7)]
        sc.assign_host_cache(rows)
        got = dict(((r["scenario"], r["arm"]), r["host_cache"]) for r in rows)
        self.assertEqual(got, {("X", "mcp"): "not_reached", ("X", "direct-local"): "cold",
                               ("Y", "mcp"): "cold", ("Y", "direct-local"): "warm",
                               ("Z", "mcp"): "cold", ("Z", "direct-local"): "warm",
                               ("W", "mcp"): "cold", ("W", "direct-local"): "warm",
                               ("V", "mcp"): "not_reached", ("V", "direct-local"): "cold",
                               ("U", "mcp"): "not_reached", ("U", "direct-local"): "cold"})
        # direct lines are unchanged: anything but a refused connection reached the host
        self.assertTrue(sc.reached_host(hline("D", "direct-remote", 0, T0, http=503)))
        self.assertFalse(sc.reached_host(hline("D", "direct-remote", 0, T0, http=0,
                                               error={"code": "transport", "text": "refused"})))

    def test_clock_offset(self):
        # mex's clock 2 s ahead: its line stamped T0 + 2.5 ran at T0 + 0.5 on the laptop, before the remote one
        rows = [hline("A", "direct-remote", 0, T0 + 1.0), hline("A", "direct-local", 0, T0 + 2.5)]
        sc.assign_host_cache(rows, {"direct-local": 2.0})
        self.assertEqual([r["host_cache"] for r in rows], ["warm", "cold"])
        sc.assign_host_cache(rows)
        self.assertEqual([r["host_cache"] for r in rows], ["cold", "warm"])


class PairedPlan(unittest.TestCase):
    def test_plan_and_resume(self):
        scen = [{"id": "C-1", "area": "C"}, {"id": "P-1", "area": "P"},
                {"id": "L-1", "area": "L", "timing_class": "long"},
                {"id": "M-1", "area": "M", "mcp": {"not_implemented": "max_mismatches"}},
                {"id": "D-1", "area": "D", "repeat": 3}]
        manifest = {"scenarios": scen}
        with tempfile.TemporaryDirectory() as d:
            plan = sc.paired_plan(manifest, d, arms=("direct-local", "direct-remote", "mcp"))
            got = [(p["scenario"], p["direct-local"], p["direct-remote"], p["mcp_cold"], p["mcp_warm"]) for p in plan]
            self.assertEqual(got, [("C-1", True, True, True, [1]), ("P-1", True, True, True, [1]),
                                   ("L-1", True, None, None, []), ("M-1", True, True, True, []),
                                   ("D-1", True, True, True, [1, 2])])
            # C-1 (index 0) has a warm repeat on the direct arms: one recorded line leaves it to do
            os.makedirs(os.path.join(d, "arms", "direct-local"))
            os.makedirs(os.path.join(d, "arms", "mcp"))
            sc.append_jsonl(os.path.join(d, "arms", "direct-local", "executions.jsonl"),
                            {"scenario": "C-1", "arm": "direct-local", "repeat": 0})
            sc.append_jsonl(os.path.join(d, "arms", "direct-local", "executions.jsonl"),
                            {"scenario": "P-1", "arm": "direct-local", "repeat": 0})
            sc.append_jsonl(os.path.join(d, "arms", "mcp", "executions.jsonl"),
                            {"scenario": "D-1", "arm": "mcp", "repeat": 1})
            plan = dict((p["scenario"], p) for p in sc.paired_plan(manifest, d, arms=("direct-local", "mcp")))
            self.assertTrue(plan["C-1"]["direct-local"])
            self.assertFalse(plan["P-1"]["direct-local"])
            self.assertIsNone(plan["P-1"]["direct-remote"])
            self.assertEqual(plan["D-1"]["mcp_warm"], [2])
            with contextlib.redirect_stdout(io.StringIO()) as buf:
                mpath = os.path.join(d, "m.json")
                with open(mpath, "w") as fh:
                    json.dump(manifest, fh)
                sc._main(["paired-plan", mpath, d, "", "", "10", "direct-local,mcp"])
            # the last four fields: the warm-ups to run again before it per arm, and on direct-local before its
            # first MCP step (none here: no warm-up)
            self.assertEqual(buf.getvalue().splitlines()[:3], ["C-1 1 - 1 1 - - - -", "P-1 0 - 1 1 - - - -",
                                                                "L-1 1 - - - - - - -"])
            with contextlib.redirect_stdout(io.StringIO()) as buf:
                sc._main(["beyond-status", mpath, "", ""])
            self.assertEqual(buf.getvalue().strip(), "")


class Lock(unittest.TestCase):
    def test_etime(self):
        self.assertEqual(sc.parse_etime("05:07"), 307)
        self.assertEqual(sc.parse_etime("01:00:00"), 3600)
        self.assertEqual(sc.parse_etime("2-01:00:01"), 2 * 86400 + 3601)
        self.assertIsNone(sc.parse_etime(""))
        self.assertIsNone(sc.parse_etime("x:y"))

    def test_verdict(self):
        lock = {"pid": 42, "host": "laptop", "started_utc": "t", "started_epoch": 1000.0, "run_id": "R"}
        alive = lambda pid: True  # noqa: E731
        dead = lambda pid: False  # noqa: E731
        self.assertEqual(sc.lock_verdict(lock, "laptop", 2000.0, alive, lambda pid: 1500.0)[0], "live")
        self.assertEqual(sc.lock_verdict(lock, "laptop", 2000.0, dead, lambda pid: None)[0], "stale")
        # pid 42 now belongs to a process started 100 s ago, long after the lock: reused
        self.assertEqual(sc.lock_verdict(lock, "laptop", 2000.0, alive, lambda pid: 100.0)[0], "stale")
        self.assertEqual(sc.lock_verdict(lock, "other", 2000.0, dead, lambda pid: None)[0], "foreign")
        self.assertEqual(sc.lock_verdict(None, "laptop", 2000.0, alive, lambda pid: 1.0)[0], "stale")
        with contextlib.redirect_stdout(io.StringIO()) as buf:
            # this process holds a lock taken now, after it started
            import time
            sc._main(["lock-check", json.dumps(dict(lock, pid=os.getpid(), host="h", started_epoch=time.time())),
                      "h"])
        self.assertTrue(buf.getvalue().startswith("live "), buf.getvalue())
        with contextlib.redirect_stdout(io.StringIO()) as buf:
            sc._main(["clock-offset", "100.0", "103.5", "101.0"])
        self.assertEqual(buf.getvalue().strip(), "3.000")


class McpStep(unittest.TestCase):
    def test_request_bound(self):
        self.assertEqual(sc.request_bound_s({"mcp": {"arguments": {"time_budget_ms": 1500}},
                                             "request": {"time_budget_ms": 9000}}), 1.5)
        self.assertEqual(sc.request_bound_s({"request": {"strategy": {"bounds": {"time_budget_ms": 5000}}}}), 5.0)
        self.assertEqual(sc.request_bound_s({"route": "POST /pattern", "request": {}}), 60.0)
        self.assertEqual(sc.request_bound_s({"route": "POST /traverse", "timing_class": "slow"}), 600.0)
        self.assertEqual(sc.request_bound_s({"route": "POST /traverse"}), 120.0)

    def test_client_timeout_exceeds_a_stated_budget(self):
        # v6: a /pattern budget of 110,000 ms leaves the fast class's 120 s 10 s above it; the client timeout is the
        # budget + the margin
        pat = {"route": "POST /pattern", "timing_class": "fast", "request": {"time_budget_ms": 110000}}
        self.assertEqual(sc.stated_request_bound_s(pat), 110.0)
        self.assertEqual(sc.client_timeout_s(pat), 110.0 + sc.CLIENT_BUDGET_MARGIN_S)
        self.assertEqual(sc.client_timeout_s(dict(pat, timing_class="slow")), 600.0)        # the class is more
        self.assertEqual(sc.client_timeout_s(dict(pat, request={"time_budget_ms": 60000})), 120.0)   # v5's budget
        self.assertEqual(sc.client_timeout_s(dict(pat, request={"time_budget_ms": 300})), 120.0)
        self.assertEqual(sc.client_timeout_s(dict(pat, mcp={"arguments": {"time_budget_ms": 200000}})), 230.0)
        self.assertIsNone(sc.stated_request_bound_s({"route": "POST /pattern", "request": {}}))
        self.assertEqual(sc.client_timeout_s({"route": "POST /pattern", "request": {}}), 120.0)   # the default budget
        self.assertEqual(sc.client_timeout_s({"route": "POST /traverse", "timing_class": "long"}), 1800.0)
        self.assertEqual(sc.client_timeout_s(pat, {"fast": 10.0, "slow": 20.0, "long": 30.0}), 140.0)
        self.assertEqual(sc.CLIENT_BUDGET_MARGIN_S, 30.0)

    def test_outcomes_and_the_wait_after_a_cancelled_job(self):
        manifest = {"scenarios": [{"id": "A", "route": "POST /traverse",
                                   "mcp": {"arguments": {"strategy": {"bounds": {"time_budget_ms": 4000}}}}},
                                  {"id": "B", "route": "POST /pattern", "mcp": {"arguments": {}}},
                                  {"id": "C", "route": "POST /pattern", "mcp": {"arguments": {}}}]}

        def ex(sid, outcome, job=None):
            return {"scenario": sid, "arm": "mcp", "repeat": 0, "outcome": outcome, "job": job,
                    "ended_utc": sc.iso_utc(T0)}
        with tempfile.TemporaryDirectory() as d:
            ej = os.path.join(d, "executions.jsonl")
            jj = os.path.join(d, "jobs.jsonl")
            sc.append_jsonl(ej, ex("A", "error", {"search_id": None, "job_id": "j1"}))
            sc.append_jsonl(ej, ex("B", "fail", {"search_id": "s2", "job_id": None}))
            sc.append_jsonl(jj, {"event": "created", "id": "j1", "at": sc.iso_utc(T0 - 300)})
            sc.append_jsonl(jj, {"event": "canceled", "id": "j1", "detail": "ok", "at": sc.iso_utc(T0)})
            outcomes, wait = sc.mcp_step_result(d, manifest, ["A", "B", "C"], 0, T0 + 1.0)
            self.assertEqual(outcomes, ["error", "fail", "missing"])
            self.assertAlmostEqual(wait, 4.0 + sc.CANCEL_MARGIN_S - 1.0)      # the cancel + A's budget + margin
            self.assertEqual(sc.mcp_step_result(d, manifest, ["A"], 0, T0 + 60.0)[1], 0.0)   # long past
            # a cancel that found the job terminal already, or a terminal event: no wait
            sc.append_jsonl(ej, ex("C", "error", {"search_id": "s3", "job_id": None}))
            sc.append_jsonl(jj, {"event": "canceled", "id": "s3", "detail": "already_terminal",
                                 "at": sc.iso_utc(T0)})
            self.assertEqual(sc.mcp_step_result(d, manifest, ["C"], 0, T0 + 1.0), (["error"], 0.0))
            sc.append_jsonl(jj, {"event": "terminal", "id": "j1", "detail": "stopped", "at": sc.iso_utc(T0 + 0.5)})
            self.assertEqual(sc.mcp_step_result(d, manifest, ["A"], 0, T0 + 1.0), (["error"], 0.0))
            # an error without a job (a transport failure at submit) has nothing to wait for
            sc.append_jsonl(ej, dict(ex("B", "error"), repeat=1))
            self.assertEqual(sc.mcp_step_result(d, manifest, ["B"], 1, T0 + 1.0), (["error"], 0.0))
            with contextlib.redirect_stdout(io.StringIO()) as buf:
                mpath = os.path.join(d, "m.json")
                with open(mpath, "w") as fh:
                    json.dump(manifest, fh)
                sc._main(["mcp-step", d, mpath, "B,C", "0"])
            self.assertEqual(buf.getvalue().strip(), "fail,error 0")


class LockIdentity(unittest.TestCase):
    def test_host_id(self):
        # review N6: the lock names this machine by a stable id, not its DHCP network name
        self.assertEqual(sc.host_id(env={"SUITE_HOST_ID": "fixed"}, system="Darwin", scutil=lambda: "x"), "fixed")
        self.assertEqual(sc.host_id(env={}, system="Darwin", scutil=lambda: "mbp-GR-3",
                                    node="staff-net-vpn-dhcp-1262.intern.ethz.ch"), "mbp-GR-3")
        self.assertEqual(sc.host_id(env={}, system="Darwin", scutil=lambda: "", node="laptop"), "laptop")
        self.assertEqual(sc.host_id(env={}, system="Linux", scutil=lambda: "never", node="box"), "box")
        self.assertTrue(sc.host_id())

    def test_foreign_reason_names_the_network_name(self):
        lock = {"pid": 7, "host": "mbp-other", "hostname": "dhcp-9.example", "started_utc": "t"}
        verdict, reason = sc.lock_verdict(lock, "mbp-GR-3", 0.0, lambda pid: True, lambda pid: None)
        self.assertEqual(verdict, "foreign")
        self.assertIn("on mbp-other (dhcp-9.example)", reason)

    def test_lock_break_command(self):
        # the stale lock is removed only while it still holds what was judged stale, inside a mutex
        held = json.dumps({"pid": 1, "host": "h", "token": "old", "note": "it's $(quoted) `x`"})
        with tempfile.TemporaryDirectory() as d:
            lock = os.path.join(d, "LOCK")

            def run():
                return subprocess.run(["bash", "-c", sc.lock_break_command(d, held)]).returncode

            with open(lock, "w") as fh:
                fh.write(held + "\n")
            self.assertEqual(run(), 0)
            self.assertFalse(os.path.exists(lock))
            self.assertFalse(os.path.exists(lock + ".break"))
            with open(lock, "w") as fh:
                fh.write(json.dumps({"pid": 2, "token": "new"}) + "\n")   # taken by another one since
            self.assertEqual(run(), 4)
            self.assertTrue(os.path.exists(lock))
            os.mkdir(lock + ".break")                                    # another removal in progress
            with open(lock, "w") as fh:
                fh.write(held + "\n")
            self.assertEqual(run(), 3)
            self.assertTrue(os.path.exists(lock))
            os.rmdir(lock + ".break")
            os.unlink(lock)
            self.assertEqual(run(), 4)                                   # gone already
            self.assertFalse(os.path.exists(lock + ".break"))


class ExecutionOutcome(unittest.TestCase):
    def test_rule(self):
        o = lambda *outs: sc.outcome_from_checks([{"check": "c%d" % i, "outcome": x} for i, x in enumerate(outs)])  # noqa: E731
        self.assertEqual(o("pass", "pass"), "pass")
        self.assertEqual(o("pass", "not_exposed"), "pass")
        self.assertEqual(o("pass", "fail", "finding"), "fail")
        self.assertEqual(o("pass", "finding", "inconclusive"), "finding")
        self.assertEqual(o("pass", "inconclusive"), "inconclusive")
        self.assertEqual(o("not_exposed", "not_exposed"), "inconclusive")
        self.assertEqual(o(), "inconclusive")
        self.assertEqual(o("pass", "weird"), "fail")  # an unknown outcome reads as a failure


class ChecksStub(unittest.TestCase):
    def test_status_pass_fail(self):
        scen = {"expect": [{"check": "status", "args": {"http": 400, "code": "invalid_request"}},
                           {"check": "contains_source_context", "args": {}}]}
        res = checks_stub.evaluate(scen, {"http": 400, "error": {"code": "invalid_request", "text": "x"}},
                                   "direct-local")
        self.assertEqual([r["outcome"] for r in res], ["pass", "not_exposed"])
        res = checks_stub.evaluate(scen, {"http": 200}, "direct-local")
        self.assertEqual(res[0]["outcome"], "fail")
        # the code is compared, as lane S1's status check does
        res = checks_stub.evaluate(scen, {"http": 400, "error": {"code": "other", "text": "x"}}, "direct-local")
        self.assertEqual(res[0]["outcome"], "fail")

    def test_implicit_status(self):
        res = checks_stub.evaluate({"expect": []}, {"http": 204}, "direct-remote")
        self.assertEqual(res, [{"check": "status", "outcome": "pass",
                                "detail": "http 204, expected 2xx (implicit: no status expectation)"}])
        self.assertEqual(checks_stub.evaluate({}, {"http": 0}, "x")[0]["outcome"], "fail")


class SamplerParsing(unittest.TestCase):
    LINES = (
        '{"BlockIO":"126GB / 0B","CPUPerc":"0.00%","Container":"metagraph_refseq33m_experimental",'
        '"ID":"a92f51b0c9ee","MemPerc":"9.27%","MemUsage":"11.87GiB / 128GiB",'
        '"Name":"metagraph_refseq33m_experimental","NetIO":"141kB / 3.85MB","PIDs":"8"}\n'
        '{"BlockIO":"440TB / 183GB","CPUPerc":"7704.04%","Container":"metagraph_sra_all","ID":"bbf00438a6b5",'
        '"MemPerc":"69.51%","MemUsage":"787.9GiB / 1.107TiB","Name":"metagraph_sra_all",'
        '"NetIO":"12.7GB / 188MB","PIDs":"513"}\n'
        'Error response from daemon: something\n')

    def test_sizes(self):
        self.assertEqual(sampler.parse_size("11.87GiB"), int(round(11.87 * 1024 ** 3)))
        self.assertEqual(sampler.parse_size("141kB"), 141000)
        self.assertEqual(sampler.parse_size("0B"), 0)
        self.assertEqual(sampler.parse_size("3.85MB"), 3850000)
        self.assertIsNone(sampler.parse_size("--"))
        self.assertIsNone(sampler.parse_size("12XB"))
        self.assertEqual(sampler.parse_pair("1kB / 2kB"), (1000, 2000))
        self.assertEqual(sampler.parse_pair("--"), (None, None))
        self.assertEqual(sampler.parse_percent("7704.04%"), 7704.04)
        self.assertIsNone(sampler.parse_percent("--"))
        # SD-018: a non-finite figure is no reading (the readers treat it as none too)
        for bad in ("NaN%", "nan", "inf%", "-inf%", "Infinity%"):
            self.assertIsNone(sampler.parse_percent(bad), bad)

    def test_stats_and_sample(self):
        stats = sampler.parse_stats_lines(self.LINES)
        self.assertEqual(set(stats), {"metagraph_refseq33m_experimental", "metagraph_sra_all"})
        s = sampler.build_sample("2026-10-13T10:00:00.000Z", stats, 240.19, sampler.CONTAINER,
                                 sampler.SRA_CONTAINER, label="direct-local")
        self.assertEqual(s["cpu_pct"], 0.0)
        self.assertEqual(s["mem_limit_bytes"], 128 * 1024 ** 3)
        self.assertEqual(s["blk_read"], 126 * 1000 ** 3)
        self.assertEqual(s["net_tx"], 3850000)
        self.assertEqual(s["sra_all_cpu_pct"], 7704.04)
        self.assertEqual(s["host_load1"], 240.19)
        self.assertEqual(s["phase"], "direct-local")
        for key in ("t_utc", "container", "cpu_pct", "mem_bytes", "mem_limit_bytes", "net_rx", "net_tx",
                    "blk_read", "blk_write", "host_load1", "sra_all_cpu_pct"):
            self.assertIn(key, s)

    def test_missing_container(self):
        s = sampler.build_sample("t", {}, None, sampler.CONTAINER, sampler.SRA_CONTAINER, error="no docker")
        self.assertIsNone(s["cpu_pct"])
        self.assertEqual(s["error"], "no docker")



class DeadlineFloor(unittest.TestCase):
    """SD-020: the host answers its 503 deadline only once the budget has elapsed (Deadline::respond_expired; the
    finalisation reserve ends the work earlier, not the answer): on the service-side wall, a relayed /resolve deadline
    (traversal_resolve_deadline, host_state_verdict) and a deadline-ended pattern task (wall_rule_verdict,
    budget_minus_reserve) are plausible from the same floor, the budget less 1 s of clock and phase slack."""

    def test_one_floor_for_the_host_deadline(self):
        scen = {"id": "R-x", "route": "POST /resolve", "expect": [],
                "mcp": {"tool": "traverse_resolve", "arguments": {"time_budget_ms": 20000}}}
        err = {"code": "traversal_resolve_deadline", "text": "the backend could not build and write its answer "
               "within bounds.time_budget_ms (20000 ms, the finalisation reserve of 250 ms included)"}
        self.assertEqual(sc.resolve_deadline_least(scen, None), (19000.0, "time_budget_ms - 1000 ms"))
        for wall, plausible in ((18900.0, False), (18999.9, False), (19000.0, True), (19100.0, True),
                                (19700.0, True), (19749.9, True), (19750.0, True), (20000.0, True)):
            doc = {"phase_sources": {"attempts": [{"span_ms": wall}]}}
            state = sc.host_state_verdict(scen, {}, doc, err, {"code": err["code"]})[1]
            rule = sc.wall_rule_verdict("budget_minus_reserve", "c", scen, {"error": err}, {}, doc)[0]
            self.assertEqual((state, rule), ("host_state" if plausible else "not_plausible",
                                             "holds" if plausible else "held_back"), wall)

if __name__ == "__main__":
    unittest.main()


class RawRecordReader(unittest.TestCase):
    """SD-093 (option 2): one reader of an execution's raw record for every caller (suite_common.raw_record_path,
    read_raw_record, raw_record_answer, raw_record_body), today's lookup order: an absolute path as it is; under a run
    directory run_dir/<raw>, then run_dir/arms/<arm>/<raw>; under an arm directory arm_dir/<raw>, then the run
    directory for a merged 'arms/...' path. An unreadable record (missing, a truncated .gz, a cut JSON) or one that is
    no object reads as none."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.run = os.path.join(self.tmp, "run")
        self.arm = os.path.join(self.run, "arms", "mcp")
        os.makedirs(os.path.join(self.arm, "raw"))
        self.rec = {"scenario": "S", "answer": {"http": 200, "body": {"x": 1}}}
        sc.write_gz_json(os.path.join(self.arm, "raw", "S.mcp.0.json.gz"), self.rec)

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_lookup(self):
        arm_line = {"arm": "mcp", "raw": "raw/S.mcp.0.json.gz"}
        merged = {"arm": "mcp", "raw": "arms/mcp/raw/S.mcp.0.json.gz"}
        absolute = {"arm": "mcp", "raw": os.path.join(self.arm, "raw", "S.mcp.0.json.gz")}
        want = os.path.join(self.arm, "raw", "S.mcp.0.json.gz")
        for line, where in ((merged, {"run_dir": self.run}), (arm_line, {"run_dir": self.run}),
                            (arm_line, {"arm_dir": self.arm}), (merged, {"arm_dir": self.arm}),
                            (absolute, {"run_dir": self.run}), (absolute, {"arm_dir": self.arm}), (absolute, {})):
            self.assertEqual(os.path.realpath(sc.raw_record_path(line, **where)), os.path.realpath(want), (line, where))
            self.assertEqual(sc.read_raw_record(line, **where), self.rec, (line, where))
            self.assertEqual(sc.raw_record_answer(line, **where), self.rec["answer"])
            self.assertEqual(sc.raw_record_body(line, **where), {"x": 1})
        # a top-level raw/ under the run directory (the comparison's test layout) comes first
        os.makedirs(os.path.join(self.run, "raw"))
        top = dict(self.rec, scenario="top")
        sc.write_gz_json(os.path.join(self.run, "raw", "S.mcp.0.json.gz"), top)
        self.assertEqual(sc.read_raw_record(arm_line, run_dir=self.run), top)
        self.assertEqual(sc.read_raw_record(arm_line, arm_dir=self.arm), self.rec)
        # nothing to read: no raw, not a string, no file, no directory given for a relative path
        for line, where in (({}, {"run_dir": self.run}), ({"raw": 5}, {"run_dir": self.run}), (None, {}),
                            ({"raw": "raw/none.json.gz", "arm": "mcp"}, {"run_dir": self.run}),
                            (arm_line, {})):
            self.assertIsNone(sc.raw_record_path(line, **where), line)
            self.assertIsNone(sc.read_raw_record(line, **where), line)
            self.assertIsNone(sc.raw_record_answer(line, **where), line)
            self.assertIsNone(sc.raw_record_body(line, **where), line)

    def test_unreadable_and_odd_records(self):
        path = os.path.join(self.arm, "raw", "S.mcp.0.json.gz")
        line = {"arm": "mcp", "raw": "raw/S.mcp.0.json.gz"}
        with open(path, "rb") as fh:
            data = fh.read()
        sentinel = object()
        for content, any_json in ((data[:len(data) // 2], sentinel),          # a truncated .gz: EOFError
                                  (gzip.compress(b'{"answer": '), sentinel),    # a cut JSON: ValueError
                                  (b"not gzip at all", sentinel),               # no gzip: OSError
                                  (gzip.compress(b"[1, 2]"), [1, 2]),           # a record that is no object
                                  (gzip.compress(b"null"), None),
                                  (gzip.compress(b'{"answer": [1]}'), {"answer": [1]})):
            with open(path, "wb") as fh:
                fh.write(content)
            self.assertIsNone(sc.raw_record_answer(line, arm_dir=self.arm), content)
            self.assertIsNone(sc.raw_record_body(line, arm_dir=self.arm), content)
            got = sc.read_raw_record(line, arm_dir=self.arm, any_json=True, default=sentinel)
            self.assertEqual(got, any_json, content)
            if not isinstance(any_json, dict):
                self.assertIsNone(sc.read_raw_record(line, arm_dir=self.arm), content)

    def test_one_reader(self):
        """Every reader of a raw record goes through suite_common (SD-093): no other module of the suite opens one
        with read_gz_json."""
        import re
        suite = os.path.dirname(HERE)
        found = []
        for sub in ("", "bench", "doc"):
            folder = os.path.join(suite, sub)
            for name in sorted(os.listdir(folder)):
                if not name.endswith(".py") or name == "suite_common.py":
                    continue
                with open(os.path.join(folder, name)) as fh:
                    if re.search(r"\bread_gz_json\(", fh.read()):
                        found.append(os.path.join(sub, name))
        self.assertEqual(found, [])
