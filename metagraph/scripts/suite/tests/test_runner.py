"""Offline end-to-end tests of run_direct.py against a local fake host (tests/fake_host.py): the execution
lines of DESIGN §8.3, raw answers, run.json identity, resume, the request cap, timeouts, transport failures,
the long class behind the front proxy, the checks fallback and a clean stop on SIGINT."""

import json
import os
import signal
import subprocess
import sys
import tempfile
import time
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.dirname(HERE)
sys.path.insert(0, SUITE)
sys.path.insert(0, HERE)

import run_direct  # noqa: E402
import suite_common as sc  # noqa: E402
from fake_host import FakeHost, RELEASE, SERVER_INSTANCE  # noqa: E402

EXEC_FIELDS = ("run_id", "scenario", "arm", "repeat", "cache", "started_utc", "wall_ms", "http", "outcome",
               "checks", "server", "mcp_phases", "answer_digest", "raw", "job")


def scenario(sid, area, route, request, http=200, timing_class="fast", repeat=1, **extra):
    s = {"id": sid, "area": area, "title": sid, "route": route, "request": request,
         "mcp": {"not_implemented": "test"}, "truth": {},
         "expect": [{"check": "status", "args": {"http": http}}], "timing_class": timing_class,
         "repeat": repeat, "tags": []}
    s.update(extra)
    return s


def make_manifest(path, scenarios):
    m = {"manifest_version": 1, "content_sha256": sc.scenarios_sha256(scenarios), "created_utc": "x",
         "index": {"database": "test"}, "ground_truth": {}, "scenarios": scenarios}
    with open(path, "w") as fh:
        json.dump(m, fh)
    return m


def basic_scenarios():
    return [
        scenario("C-1", "C", "GET /capabilities", None),
        scenario("P-1", "P", "POST /pattern", {"patterns": [{"dna": "ACGT" * 8, "id": "a"}], "mode": "count"}),
        scenario("R-1", "R", "POST /resolve", {"sequence": "ACGT" * 10, "labels": ["562"]}),
        scenario("T-1", "T", "POST /traverse", {"seeds": [{"sequence": "ACGT" * 10}]}),
        scenario("P-2", "P", "POST /pattern", {"patterns": [], "bogus": 1}, http=400),
        scenario("D-1", "D", "POST /resolve", {"sequence": "ACGT" * 10, "labels": ["562"]}, repeat=3),
    ]


def with_expect(s, *expect):
    s["expect"] = s["expect"] + list(expect)
    return s


class RunnerTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = self.tmp.name
        self.manifest_path = os.path.join(self.dir, "m.json")

    def tearDown(self):
        signal.signal(signal.SIGINT, signal.default_int_handler)
        signal.signal(signal.SIGTERM, signal.SIG_DFL)
        if hasattr(signal, "SIGHUP"):
            signal.signal(signal.SIGHUP, signal.SIG_DFL)
        self.tmp.cleanup()

    def run_direct(self, base, out, *extra):
        argv = ["--manifest", self.manifest_path, "--arm", "direct-local", "--base", base, "--out", out,
                "--pause", "0", "--quiet", "--checks", "stub"] + list(extra)
        return run_direct.main(argv)

    def test_time_limited_repeat_and_client_wall(self):
        # VERIFY-fixes.md L4: repeat 0 stopped by its time budget, repeat 1 complete: the runner passes {digest,
        # time_limited} for the earlier repeats (also for a line of before the field, from its raw answer), so
        # the determinism check is inconclusive, never a fail; and a 503 deadline long past its budget is a
        # finding by the client's wall time the runner passes to budget_respected
        req = {"patterns": [{"dna": "A" * 31, "id": "a"}], "time_budget_ms": 1000}
        late = {"sequence": "ACGT" * 10, "labels": ["562"], "bounds": {"time_budget_ms": 100}, "fake_503_after_s": 1.3}
        make_manifest(self.manifest_path, [
            with_expect(scenario("D-T", "D", "POST /timed", req, repeat=2), {"check": "deterministic", "args": {}}),
            with_expect(scenario("R-late", "R", "POST /resolve", late, http=[200, 503]),
                        {"check": "budget_respected", "args": {}})])
        out = os.path.join(self.dir, "run")
        with FakeHost() as host:
            self.assertEqual(self.run_direct(host.base, out, "--checks", "checks", "--warm-every", "0"), 0)
            rows, _ = sc.read_jsonl(os.path.join(out, "executions.jsonl"))
            by = dict(((r["scenario"], r["repeat"]), r) for r in rows)
            self.assertEqual(by[("D-T", 0)]["time_limited"], "pattern 0 time-limited")
            self.assertEqual(by[("D-T", 1)]["time_limited"], "")
            self.assertNotEqual(by[("D-T", 0)]["answer_digest"], by[("D-T", 1)]["answer_digest"])
            det = [c for c in by[("D-T", 1)]["checks"] if c["check"] == "deterministic"][0]
            self.assertEqual(det["outcome"], "inconclusive", det["detail"])
            self.assertIn("no earlier repeat without a time limitation", det["detail"])
            self.assertEqual(by[("D-T", 1)]["outcome"], "inconclusive")
            budget = [c for c in by[("R-late", 0)]["checks"] if c["check"] == "budget_respected"][0]
            self.assertEqual(budget["outcome"], "finding", budget["detail"])
            self.assertEqual(by[("R-late", 0)]["outcome"], "finding")
            # resume after a runner of before the field: repeat 1 again, repeat 0's time limitation from its raw
            kept = [dict((k, v) for k, v in r.items() if k != "time_limited") for r in rows
                    if (r["scenario"], r["repeat"]) != ("D-T", 1)]
            with open(os.path.join(out, "executions.jsonl"), "w") as fh:
                for r in kept:
                    fh.write(json.dumps(r) + "\n")
            self.assertEqual(self.run_direct(host.base, out, "--checks", "checks", "--warm-every", "0"), 0)
        rows, _ = sc.read_jsonl(os.path.join(out, "executions.jsonl"))
        again = [r for r in rows if (r["scenario"], r["repeat"]) == ("D-T", 1)][-1]
        det = [c for c in again["checks"] if c["check"] == "deterministic"][0]
        self.assertEqual(det["outcome"], "inconclusive", det["detail"])

    def test_full_run_lines_raw_and_identity(self):
        make_manifest(self.manifest_path, basic_scenarios())
        out = os.path.join(self.dir, "runs", "R1", "arms", "direct-local")
        with FakeHost() as host:
            rc = self.run_direct(host.base, out)
        self.assertEqual(rc, 0)
        rows, bad = sc.read_jsonl(os.path.join(out, "executions.jsonl"))
        self.assertEqual(bad, 0)
        # 6 scenarios; index 0 gets a warm repeat; D-1 runs 3 times
        self.assertEqual(len(rows), 6 + 1 + 2)
        for r in rows:
            for f in EXEC_FIELDS:
                self.assertIn(f, r)
            self.assertEqual(r["run_id"], "R1")
            self.assertTrue(os.path.exists(os.path.join(out, r["raw"])))
        by = dict(((r["scenario"], r["repeat"]), r) for r in rows)
        self.assertEqual(by[("C-1", 0)]["cache"], "cold")
        self.assertEqual(by[("C-1", 1)]["cache"], "warm")
        self.assertEqual(by[("C-1", 0)]["answer_digest"], by[("C-1", 1)]["answer_digest"])
        self.assertEqual(by[("P-1", 0)]["outcome"], "pass")
        self.assertEqual(by[("P-1", 0)]["client"]["content_encoding"], "gzip")
        self.assertEqual(by[("P-1", 0)]["server"]["elapsed_ms"], 2.0)
        self.assertEqual(by[("P-1", 0)]["server"]["work"]["steps"], 10)
        self.assertEqual(by[("T-1", 0)]["server"]["work"]["counters.steps"], 12)
        self.assertEqual(by[("P-2", 0)]["http"], 400)
        self.assertEqual(by[("P-2", 0)]["outcome"], "pass")
        self.assertEqual(by[("P-2", 0)]["error"]["code"], "invalid_request")
        self.assertEqual([by[("D-1", i)]["cache"] for i in range(3)], ["cold", "warm", "warm"])
        raw = sc.read_gz_json(os.path.join(out, by[("P-2", 0)]["raw"]))
        self.assertEqual(raw["answer"]["error"]["text"], "request: unknown field 'bogus'")
        self.assertEqual(raw["request"]["body"], {"patterns": [], "bogus": 1})
        run = sc.read_json(os.path.join(out, "run.json"))
        self.assertEqual(run["release"], RELEASE)
        self.assertEqual(run["server_instance"], SERVER_INSTANCE)
        self.assertEqual(run["index_fp"], "907ec986aaaa")
        self.assertEqual(set(run["capabilities"]), {"/capabilities", "/pattern/capabilities",
                                                    "/traverse/capabilities"})
        self.assertEqual(run["checks_module"], "stub")
        self.assertEqual(run["state"], "done")
        self.assertFalse(run["server_instance_changed"])
        self.assertTrue(run["manifest"]["sha_verified"])
        self.assertEqual(len(run["sessions"]), 1)
        self.assertIsNotNone(run["host_load_start"])
        self.assertTrue(os.path.exists(os.path.join(out, "manifest.json")))
        state = sc.read_json(os.path.join(out, "state.json"))
        self.assertEqual((state["state"], state["done"], state["planned"]), ("done", 9, 9))
        # 3 capabilities + 9 executions + 1 closing capabilities, all under the prefix
        self.assertEqual(len(host.log), 13)

    def test_resume_skips_recorded(self):
        make_manifest(self.manifest_path, basic_scenarios())
        out = os.path.join(self.dir, "run")
        with FakeHost() as host:
            self.assertEqual(self.run_direct(host.base, out), 0)
            n_first = len(host.log)
            self.assertEqual(self.run_direct(host.base, out), 0)
            # only the capabilities GETs again
            self.assertEqual(len(host.log) - n_first, 4)
            rows, _ = sc.read_jsonl(os.path.join(out, "executions.jsonl"))
            self.assertEqual(len(rows), 9)
            # drop the last line (and tear the file): only that execution runs again
            with open(os.path.join(out, "executions.jsonl")) as fh:
                lines = fh.readlines()
            with open(os.path.join(out, "executions.jsonl"), "w") as fh:
                fh.writelines(lines[:-1])
                fh.write('{"torn": ')
            self.assertEqual(self.run_direct(host.base, out), 0)
        rows, bad = sc.read_jsonl(os.path.join(out, "executions.jsonl"))
        self.assertEqual((len(rows), bad), (9, 1))
        self.assertTrue(rows[-1].get("after_resume"))
        run = sc.read_json(os.path.join(out, "run.json"))
        self.assertEqual(len(run["sessions"]), 3)
        self.assertEqual(run["sessions"][1]["skipped_done"], 9)

    def test_resume_refuses_other_manifest(self):
        make_manifest(self.manifest_path, basic_scenarios())
        out = os.path.join(self.dir, "run")
        with FakeHost() as host:
            self.run_direct(host.base, out, "--limit", "1")
            make_manifest(self.manifest_path, basic_scenarios()[1:])
            with self.assertRaises(SystemExit):
                self.run_direct(host.base, out)

    def test_max_requests_stops(self):
        make_manifest(self.manifest_path, basic_scenarios())
        out = os.path.join(self.dir, "run")
        with FakeHost() as host:
            rc = self.run_direct(host.base, out, "--max-requests", "5")
            # the closing GET /capabilities is inside the cap: 3 capabilities + 1 execution + 1 closing GET
            self.assertEqual(len(host.log), 5)
        self.assertEqual(rc, run_direct.EXIT_REQUEST_CAP)
        rows, _ = sc.read_jsonl(os.path.join(out, "executions.jsonl"))
        self.assertEqual(len(rows), 1)
        run = sc.read_json(os.path.join(out, "run.json"))
        self.assertEqual(run["state"], "interrupted")
        self.assertIn("max_requests", run["sessions"][0]["reason"])
        self.assertEqual(run["server_instance_end"], SERVER_INSTANCE)

    def test_areas_limit_and_ids(self):
        make_manifest(self.manifest_path, basic_scenarios())
        out = os.path.join(self.dir, "run")
        with FakeHost() as host:
            self.run_direct(host.base, out, "--areas", "P,R", "--limit", "2", "--warm-every", "0")
        rows, _ = sc.read_jsonl(os.path.join(out, "executions.jsonl"))
        self.assertEqual([r["scenario"] for r in rows], ["P-1", "R-1"])
        out2 = os.path.join(self.dir, "run2")
        with FakeHost() as host:
            self.run_direct(host.base, out2, "--ids", "T-1")
        rows, _ = sc.read_jsonl(os.path.join(out2, "executions.jsonl"))
        self.assertEqual([r["scenario"] for r in rows], ["T-1"])

    def test_timeout_and_transport_errors(self):
        scen = [scenario("X-1", "X", "POST /slow", {"sleep_s": 2.0}),
                scenario("X-2", "X", "POST /slow", {"sleep_s": 0.0})]
        make_manifest(self.manifest_path, scen)
        out = os.path.join(self.dir, "run")
        with FakeHost() as host:
            self.run_direct(host.base, out, "--timeout-fast", "0.5", "--warm-every", "0")
        rows, _ = sc.read_jsonl(os.path.join(out, "executions.jsonl"))
        self.assertEqual(rows[0]["outcome"], "error")
        self.assertEqual(rows[0]["http"], 0)
        self.assertEqual(rows[0]["error"]["code"], "timeout")
        self.assertEqual(rows[1]["outcome"], "pass")
        # nothing listens on this port: every request is a transport failure, the run still completes
        out2 = os.path.join(self.dir, "run2")
        rc = self.run_direct("http://127.0.0.1:9", out2, "--warm-every", "0")
        self.assertEqual(rc, 0)
        rows, _ = sc.read_jsonl(os.path.join(out2, "executions.jsonl"))
        self.assertEqual([r["error"]["code"] for r in rows], ["transport", "transport"])
        run = sc.read_json(os.path.join(out2, "run.json"))
        self.assertIsNone(run["release"])

    def test_long_class_skipped_remote(self):
        scen = [scenario("L-1", "L", "POST /pattern", {"patterns": []}, timing_class="long"),
                scenario("L-2", "L", "POST /pattern", {"patterns": []}, timing_class="slow")]
        make_manifest(self.manifest_path, scen)
        out = os.path.join(self.dir, "run")
        with FakeHost() as host:
            run_direct.main(["--manifest", self.manifest_path, "--arm", "direct-remote", "--base", host.base,
                             "--out", out, "--pause", "0", "--quiet", "--warm-every", "0"])
        rows, _ = sc.read_jsonl(os.path.join(out, "executions.jsonl"))
        self.assertEqual([r["scenario"] for r in rows], ["L-2"])
        self.assertEqual(rows[0]["client"]["timeout_s"], run_direct.REMOTE_TIMEOUT_CAP_S)
        run = sc.read_json(os.path.join(out, "run.json"))
        self.assertEqual(run["skipped"][0]["scenario"], "L-1")

    def test_client_timeout_follows_a_stated_budget(self):
        # v6: a /pattern request with a 110,000 ms budget (10 s below the fast class's 120 s) gets the budget + 30 s
        # on direct-local, and the front proxy's cap on direct-remote; v5's 60,000 ms keeps the class's timeout
        scen = [scenario("P-6", "P", "POST /pattern", {"patterns": [{"dna": "ACGT" * 8, "id": "a"}], "mode": "count",
                                                      "time_budget_ms": 110000}),
                scenario("P-5", "P", "POST /pattern", {"patterns": [{"dna": "ACGT" * 8, "id": "a"}], "mode": "count",
                                                      "time_budget_ms": 60000})]
        make_manifest(self.manifest_path, scen)
        for arm, want in (("direct-local", 140.0), ("direct-remote", run_direct.REMOTE_TIMEOUT_CAP_S)):
            out = os.path.join(self.dir, "run-" + arm)
            with FakeHost() as host:
                run_direct.main(["--manifest", self.manifest_path, "--arm", arm, "--base", host.base, "--out", out,
                                 "--pause", "0", "--quiet", "--warm-every", "0", "--checks", "stub"])
            rows, _ = sc.read_jsonl(os.path.join(out, "executions.jsonl"))
            self.assertEqual([(r["scenario"], r["client"]["timeout_s"]) for r in rows],
                             [("P-6", want), ("P-5", 120.0)], arm)
        self.assertEqual(run_direct.timeout_for(scen[0], "direct-local", {"fast": 120.0}), 140.0)
        self.assertEqual(run_direct.timeout_for(scen[0], "direct-local", {"fast": 200.0}), 200.0)

    def test_protocol_extensions(self):
        scen = [scenario("C-2", "C", "POST /pattern", None, http=400, body_raw="not json",
                         headers={"Accept-Encoding": None, "Content-Type": "text/plain"})]
        make_manifest(self.manifest_path, scen)
        out = os.path.join(self.dir, "run")
        with FakeHost() as host:
            self.run_direct(host.base, out, "--warm-every", "0")
        rows, _ = sc.read_jsonl(os.path.join(out, "executions.jsonl"))
        self.assertEqual(rows[0]["http"], 400)
        raw = sc.read_gz_json(os.path.join(out, rows[0]["raw"]))
        self.assertEqual(raw["request"]["body_raw"], "not json")
        self.assertNotIn("Accept-Encoding", raw["request"]["headers"])
        self.assertEqual(raw["request"]["headers"]["Content-Type"], "text/plain")

    def test_checks_exception_is_error(self):
        class Boom(object):
            @staticmethod
            def evaluate(scenario, answer, arm, previous=None):
                raise RuntimeError("bad check")
        make_manifest(self.manifest_path, basic_scenarios()[:1])
        args = run_direct.parse_args(["--manifest", self.manifest_path, "--arm", "direct-local", "--base", "x",
                                      "--out", os.path.join(self.dir, "run"), "--pause", "0", "--quiet",
                                      "--warm-every", "0"])
        with FakeHost() as host:
            args.base = host.base
            run_direct.reset_stop()
            runner = run_direct.Runner(args, sc.read_json(self.manifest_path), Boom, "boom")
            self.assertEqual(runner.run_all(), "done")
        rows, _ = sc.read_jsonl(os.path.join(self.dir, "run", "executions.jsonl"))
        self.assertEqual(rows[0]["outcome"], "error")
        self.assertIn("bad check", rows[0]["checks"][0]["detail"])

    def test_previous_digests_passed(self):
        seen = []

        class Recorder(object):
            @staticmethod
            def evaluate(scenario, answer, arm, previous=None):
                seen.append((scenario["id"], list(previous)))
                return [{"check": "status", "outcome": "pass", "detail": ""}]
        make_manifest(self.manifest_path, basic_scenarios()[-1:])  # D-1, repeat 3
        out = os.path.join(self.dir, "run")

        def run(limit_requests):
            args = run_direct.parse_args(["--manifest", self.manifest_path, "--arm", "direct-local", "--base", "x",
                                          "--out", out, "--pause", "0", "--quiet", "--warm-every", "0",
                                          "--max-requests", str(limit_requests)])
            args.base = host.base
            run_direct.reset_stop()
            run_direct.Runner(args, sc.read_json(self.manifest_path), Recorder, "recorder").run_all()
        with FakeHost() as host:
            run(5)   # 3 capabilities + repeats 0 and 1
            run(0)   # resumed: repeat 2 sees the digests recorded in the first session
        rows, _ = sc.read_jsonl(os.path.join(out, "executions.jsonl"))
        # {digest, time_limited} per earlier repeat (a checks module without time_limited: ''); an evaluate
        # without `client` is not passed one
        d = [{"digest": r["answer_digest"], "time_limited": ""} for r in rows]
        self.assertEqual(seen, [("D-1", []), ("D-1", [d[0]]), ("D-1", [d[0], d[1]])])

    def test_unanswered_repeat_is_no_earlier_repeat(self):
        # an execution the host never answered (no HTTP status: a client timeout, a refused connection) has the digest
        # of its stand-in answer: it is not passed in `previous`, so the next repeat compares with nothing instead of
        # failing `deterministic`; in one session and on a resume
        real, _, _ = run_direct.load_checks("checks")
        seen = []

        class Spy(object):
            scenario_outcome = staticmethod(real.scenario_outcome)
            time_limited = staticmethod(real.time_limited)

            @staticmethod
            def evaluate(scenario, answer, arm, previous=None, client=None):
                seen.append((scenario["id"], answer["http"], list(previous or [])))
                return real.evaluate(scenario, answer, arm, previous=previous, client=client)

        def run(out, base, *extra):
            args = run_direct.parse_args(["--manifest", self.manifest_path, "--arm", "direct-local", "--base", base,
                                          "--out", out, "--pause", "0", "--quiet", "--warm-every", "0"] + list(extra))
            run_direct.reset_stop()
            return run_direct.Runner(args, sc.read_json(self.manifest_path), Spy, "checks").run_all()

        det = {"check": "deterministic", "args": {}}
        make_manifest(self.manifest_path, [with_expect(scenario("F-1", "F", "POST /flaky", {"sleep_s": 1.5},
                                                                repeat=2), det)])
        with FakeHost() as host:
            out = os.path.join(self.dir, "one")
            self.assertEqual(run(out, host.base, "--timeout-fast", "0.5"), "done")
            out2 = os.path.join(self.dir, "resumed")
            # r0 against a port nothing listens on (3 capabilities GETs, r0, the closing GET), r1 by a resume
            self.assertEqual(run(out2, "http://127.0.0.1:9", "--max-requests", "5"), "interrupted")
            self.assertEqual(run(out2, host.base), "done")
        self.assertEqual([(s_, h, p) for s_, h, p in seen], [("F-1", 0, []), ("F-1", 200, []), ("F-1", 0, []),
                                                            ("F-1", 200, [])])
        for path in (out, out2):
            rows, _ = sc.read_jsonl(os.path.join(path, "executions.jsonl"))
            self.assertEqual([(r["repeat"], r["http"], r["outcome"]) for r in rows],
                             [(0, 0, "error"), (1, 200, "pass")], path)
            self.assertFalse([c for c in rows[1]["checks"] if c["outcome"] == "fail"], rows[1]["checks"])
        # an answered repeat still feeds `previous` (the rule: an HTTP status and no runner error)
        self.assertTrue(sc.feeds_previous({"http": 503, "outcome": "inconclusive"}))
        self.assertFalse(sc.feeds_previous({"http": 200, "outcome": "error"}))
        self.assertFalse(sc.feeds_previous({"http": None, "outcome": "not_implemented_mcp"}))

    def test_load_checks_fallback(self):
        mod, name, err = run_direct.load_checks("stub")
        self.assertEqual((name, err), ("stub", None))
        mod, name, err = run_direct.load_checks("auto")
        self.assertIn(name, ("stub", "checks"))
        self.assertTrue(hasattr(mod, "evaluate"))

    def test_stub_refused_beyond_status(self):
        scen = basic_scenarios()[:2]
        scen[1]["expect"].append({"check": "contains_source_context", "args": {"level": 0}})
        make_manifest(self.manifest_path, scen)
        out = os.path.join(self.dir, "run")
        with FakeHost() as host:
            with self.assertRaises(SystemExit) as ctx:
                self.run_direct(host.base, out)
            self.assertIn("contains_source_context", str(ctx.exception))
            self.assertEqual(host.log, [])
            # a selection without the other check runs; --allow-stub runs everything
            self.assertEqual(self.run_direct(host.base, out, "--ids", "C-1", "--warm-every", "0"), 0)
            self.assertEqual(self.run_direct(host.base, out, "--allow-stub", "--warm-every", "0"), 0)
        rows, _ = sc.read_jsonl(os.path.join(out, "executions.jsonl"))
        self.assertEqual([r["scenario"] for r in rows], ["C-1", "P-1"])
        p1 = rows[1]["checks"]
        self.assertEqual([c["outcome"] for c in p1], ["pass", "not_exposed"])

    def test_stub_status_code_and_text(self):
        import checks_stub
        answer = {"http": 400, "body": None, "error": {"code": "invalid_request", "text": "request: unknown field 'x'"}}
        ok = checks_stub.evaluate({"expect": [{"check": "status", "args": {
            "http": [400, 422], "code": "invalid_request", "text_prefix": "request:",
            "text_contains": ["unknown field"]}}]}, answer, "direct-local")
        self.assertEqual(ok[0]["outcome"], "pass")
        bad = checks_stub.evaluate({"expect": [{"check": "status", "args": {
            "http": 400, "code": "pattern_invalid"}}]}, answer, "direct-local")
        self.assertEqual(bad[0]["outcome"], "fail")
        self.assertIn("pattern_invalid", bad[0]["detail"])
        bad = checks_stub.evaluate({"expect": [{"check": "status", "args": {
            "http": 400, "text": "request: unknown field 'y'"}}]}, answer, "direct-local")
        self.assertEqual(bad[0]["outcome"], "fail")

    def test_keep_alive_connection(self):
        make_manifest(self.manifest_path, basic_scenarios())
        out = os.path.join(self.dir, "run")
        with FakeHost() as host:
            self.assertEqual(self.run_direct(host.base, out), 0)
            conns = list(host.conns)
        # one connection for the whole invocation
        self.assertEqual(len(set(conns)), 1, conns)
        rows, _ = sc.read_jsonl(os.path.join(out, "executions.jsonl"))
        # the capabilities GETs opened it: every execution reused it and has no connect time
        self.assertTrue(all(r["client"]["reused"] for r in rows))
        self.assertTrue(all(r["client"]["connect_ms"] is None for r in rows))
        self.assertTrue(all(r["client"]["attempts"] == 1 for r in rows))

    def test_dropped_connection_resent_once(self):
        make_manifest(self.manifest_path, basic_scenarios()[:3])
        out = os.path.join(self.dir, "run")
        with FakeHost(drop=True) as host:
            self.assertEqual(self.run_direct(host.base, out, "--warm-every", "0"), 0)
            n = len(host.log)
        rows, _ = sc.read_jsonl(os.path.join(out, "executions.jsonl"))
        self.assertEqual([r["outcome"] for r in rows], ["pass"] * 3)
        for r in rows:
            # either the client saw the close in time and connected afresh, or it resent on a fresh connection
            self.assertIsNotNone(r["client"]["connect_ms"])
            self.assertLessEqual(r["client"]["attempts"], 2)
        self.assertEqual(n, 3 + 3 + 1)  # a resend never reaches the server twice: the first attempt died unsent

    def test_no_resend_after_a_slow_drop(self):
        # review N10: a reused connection dropped after the server held the request for a while is not resent
        # (the server may have processed it); one dropped at once is
        import socket
        import threading

        def serve(srv, delay, accepted):
            c, _ = srv.accept()
            accepted.append(1)
            for answer in (True, False):
                data = b""
                while b"\r\n\r\n" not in data:
                    data += c.recv(65536)
                if answer:
                    c.sendall(b"HTTP/1.1 200 OK\r\nContent-Type: application/json\r\nContent-Length: 2\r\n\r\n{}")
            time.sleep(delay)
            c.close()
            srv.settimeout(1.5)
            try:
                c2, _ = srv.accept()
                accepted.append(2)
                data = b""
                while b"\r\n\r\n" not in data:
                    data += c2.recv(65536)
                c2.sendall(b"HTTP/1.1 200 OK\r\nContent-Type: application/json\r\nContent-Length: 2\r\n\r\n{}")
                c2.close()
            except socket.timeout:
                pass

        old = run_direct.STALE_RESEND_MAX_S
        run_direct.STALE_RESEND_MAX_S = 0.2
        try:
            for delay, resent in ((0.6, False), (0.0, True)):
                srv = socket.socket()
                srv.bind(("127.0.0.1", 0))
                srv.listen(2)
                accepted = []
                t = threading.Thread(target=serve, args=(srv, delay, accepted), daemon=True)
                t.start()
                client = run_direct.Client("http://127.0.0.1:%d" % srv.getsockname()[1])
                first = client.call("GET", "/a", None, {}, 10)
                second = client.call("GET", "/b", None, {}, 10)
                client.close()
                t.join(5)
                srv.close()
                self.assertEqual(first["http"], 200)
                if resent:
                    self.assertEqual((second["http"], second["attempts"], accepted), (200, 2, [1, 2]))
                else:
                    self.assertEqual((second["http"], second["attempts"], accepted), (0, 1, [1]))
                    self.assertEqual(second["error"][0], "transport")
        finally:
            run_direct.STALE_RESEND_MAX_S = old

    def test_paired_steps_only_and_finish(self):
        make_manifest(self.manifest_path, basic_scenarios())
        out = os.path.join(self.dir, "run")
        import io
        import contextlib
        with FakeHost() as host:
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                rc = run_direct.main(["--manifest", self.manifest_path, "--arm", "direct-local", "--base",
                                      host.base, "--out", out, "--pause", "0", "--quiet", "--checks", "stub",
                                      "--only", "C-1", "--status-line"])
            self.assertEqual(rc, 0)
            status = json.loads(buf.getvalue().strip().splitlines()[-1][len(run_direct.STATUS_PREFIX):])
            self.assertEqual(status, {"state": "partial", "executed": 2, "requests": 5, "reason": None})
            self.assertEqual(len(host.log), 5)  # 3 capabilities, C-1 cold and warm, no closing GET
            self.assertEqual(self.run_direct(host.base, out, "--only", "P-1,R-1"), 0)
            self.assertEqual(len(host.log), 7)  # the capabilities are not fetched again
            run = sc.read_json(os.path.join(out, "run.json"))
            self.assertEqual((run["state"], run["planned"]), ("partial", 9))
            state = sc.read_json(os.path.join(out, "state.json"))
            self.assertEqual((state["done"], state["planned"]), (4, 9))
            for sid in ("T-1", "P-2", "D-1"):
                self.assertEqual(self.run_direct(host.base, out, "--only", sid), 0)
            n = len(host.log)
            self.assertEqual(self.run_direct(host.base, out, "--finish"), 0)
            self.assertEqual(host.log[n:], [("GET", "/capabilities")])
        run = sc.read_json(os.path.join(out, "run.json"))
        self.assertEqual(run["state"], "done")
        self.assertEqual(run["server_instance_end"], SERVER_INSTANCE)
        self.assertFalse(run["server_instance_changed"])
        self.assertEqual(len(sc.read_jsonl(os.path.join(out, "executions.jsonl"))[0]), 9)

    def test_release_at_each_sessions_start_and_end(self):
        # the host's identity per session: the start's capabilities and the closing GET; a release that changes
        # between them (a deploy during the run, or before a resume) is in run.json's releases and release_changed
        import fake_host
        make_manifest(self.manifest_path, basic_scenarios()[:2])
        out = os.path.join(self.dir, "run")
        old = fake_host.RELEASE
        try:
            with FakeHost() as host:
                self.assertEqual(self.run_direct(host.base, out, "--only", "C-1"), 0)
                run = sc.read_json(os.path.join(out, "run.json"))
                self.assertEqual((run["sessions"][0]["release_start"], run["releases"], run["release_changed"]),
                                 (RELEASE, [RELEASE], False))
                self.assertEqual(run["sessions"][0]["index_fp_start"], "907ec986aaaa")
                fake_host.RELEASE = "refseq97-00000000"
                self.assertEqual(self.run_direct(host.base, out, "--finish"), 0)
        finally:
            fake_host.RELEASE = old
        run = sc.read_json(os.path.join(out, "run.json"))
        self.assertEqual(run["release"], RELEASE)            # the run's identity: its first release
        self.assertEqual((run["sessions"][1]["release_end"], run["release_end"]), ("refseq97-00000000",) * 2)
        self.assertEqual((run["releases"], run["release_changed"]), ([RELEASE, "refseq97-00000000"], True))

    def test_sigint_stops_cleanly(self):
        scen = [scenario("S-%d" % i, "S", "POST /slow", {"sleep_s": 0.3}) for i in range(20)]
        make_manifest(self.manifest_path, scen)
        out = os.path.join(self.dir, "run")
        with FakeHost() as host:
            proc = subprocess.Popen([sys.executable, os.path.join(SUITE, "run_direct.py"), "--manifest",
                                     self.manifest_path, "--arm", "direct-local", "--base", host.base, "--out", out,
                                     "--pause", "0.1", "--quiet", "--checks", "stub", "--warm-every", "0"],
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            exec_path = os.path.join(out, "executions.jsonl")
            deadline = time.time() + 20
            while time.time() < deadline and len(sc.read_jsonl(exec_path)[0]) < 2:
                time.sleep(0.05)
            proc.send_signal(signal.SIGINT)
            stdout, stderr = proc.communicate(timeout=30)
        self.assertEqual(proc.returncode, 130, stderr.decode())
        rows, bad = sc.read_jsonl(exec_path)
        self.assertEqual(bad, 0)
        self.assertLess(len(rows), 20)
        self.assertGreaterEqual(len(rows), 2)
        run = sc.read_json(os.path.join(out, "run.json"))
        self.assertEqual(run["state"], "interrupted")
        self.assertEqual(sc.read_json(os.path.join(out, "state.json"))["state"], "interrupted")

    def test_sigint_after_the_reader_of_stdout_is_gone(self):
        # a paired step's stdout is its ssh channel: Ctrl-C ends the ssh client, then the orchestrator's SIGINT stops
        # the runner after its request in flight; its next stdout write meets a closed pipe. The run is interrupted
        # (exit 130), not failed (exit 120 with a "runner failed" BrokenPipeError warning), with or without --quiet
        scen = [scenario("S-%d" % i, "S", "POST /slow", {"sleep_s": 1.0}) for i in range(3)]
        make_manifest(self.manifest_path, scen)
        for extra in ([], ["--quiet", "--status-line"]):
            out = os.path.join(self.dir, "run%d" % len(extra))
            with FakeHost() as host:
                r, w = os.pipe()
                proc = subprocess.Popen([sys.executable, os.path.join(SUITE, "run_direct.py"), "--manifest",
                                         self.manifest_path, "--arm", "direct-local", "--base", host.base, "--out", out,
                                         "--pause", "0", "--checks", "stub", "--warm-every", "0"] + extra,
                                        stdout=w, stderr=subprocess.PIPE, stdin=subprocess.DEVNULL)
                os.close(w)
                deadline = time.time() + 30
                while time.time() < deadline and ("POST", "/slow") not in host.log:
                    time.sleep(0.02)
                os.close(r)                     # the reader is gone while the request is in flight
                time.sleep(0.2)
                proc.send_signal(signal.SIGINT)
                _, stderr = proc.communicate(timeout=60)
            self.assertEqual(proc.returncode, 130, (extra, stderr.decode()))
            self.assertNotIn(b"BrokenPipeError", stderr)
            run = sc.read_json(os.path.join(out, "run.json"))
            self.assertEqual((run["state"], run["sessions"][-1]["state"]), ("interrupted", "interrupted"), extra)
            self.assertEqual(run.get("warnings", []), [], extra)
            self.assertEqual(sc.read_json(os.path.join(out, "state.json"))["state"], "interrupted")
            rows, _ = sc.read_jsonl(os.path.join(out, "executions.jsonl"))
            self.assertEqual([x["scenario"] for x in rows], ["S-0"])    # the request in flight is recorded


if __name__ == "__main__":
    unittest.main()
