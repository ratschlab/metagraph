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
        self.assertEqual(rc, run_direct.EXIT_REQUEST_CAP)
        rows, _ = sc.read_jsonl(os.path.join(out, "executions.jsonl"))
        self.assertEqual(len(rows), 2)  # 3 capabilities + 2 executions
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
        d = [r["answer_digest"] for r in rows]
        self.assertEqual(seen, [("D-1", []), ("D-1", [d[0]]), ("D-1", [d[0], d[1]])])

    def test_load_checks_fallback(self):
        mod, name = run_direct.load_checks("stub")
        self.assertEqual(name, "stub")
        mod, name = run_direct.load_checks("auto")
        self.assertIn(name, ("stub", "checks"))
        self.assertTrue(hasattr(mod, "evaluate"))

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


if __name__ == "__main__":
    unittest.main()
