"""Offline tests of run_suite.sh (the remote arm against the local fake host; no mex contact) and of the
runner's refusal of a run directory another live runner holds."""

import json
import os
import platform
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
from fake_host import FakeHost  # noqa: E402
from test_runner import basic_scenarios, make_manifest  # noqa: E402


class LiveRunnerGuard(unittest.TestCase):
    def test_refuses_directory_of_live_runner(self):
        with tempfile.TemporaryDirectory() as d:
            mpath = os.path.join(d, "m.json")
            make_manifest(mpath, basic_scenarios()[:1])
            out = os.path.join(d, "run")
            os.makedirs(out)
            sleeper = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(30)"])
            try:
                sc.write_json_atomic(os.path.join(out, "run.json"), {"arm": "direct-local",
                                                                     "runner_host": platform.node()})
                sc.write_json_atomic(os.path.join(out, "state.json"), {"state": "running", "pid": sleeper.pid})
                with self.assertRaises(SystemExit) as ctx:
                    run_direct.main(["--manifest", mpath, "--arm", "direct-local", "--base", "http://127.0.0.1:9",
                                     "--out", out, "--quiet"])
                self.assertIn("still alive", str(ctx.exception))
            finally:
                sleeper.kill()
                sleeper.wait()


class RunSuiteOffline(unittest.TestCase):
    def run_suite(self, runs, host, *args):
        env = dict(os.environ)
        env.update({"SUITE_RUNS": runs, "SUITE_REMOTE_BASE": host.base, "SUITE_REMOTE_SAMPLER": "0",
                    "SUITE_MEX": "nobody@invalid.invalid", "PYTHON": sys.executable})
        cmd = [os.path.join(SUITE, "run_suite.sh"), "--skip-local", "--no-pdf", "--pause", "0"] + list(args)
        proc = subprocess.run(cmd, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              universal_newlines=True, timeout=120)
        return proc

    def test_remote_arm_new_and_resume(self):
        with tempfile.TemporaryDirectory() as d:
            mpath = os.path.join(d, "m.json")
            make_manifest(mpath, basic_scenarios())
            runs = os.path.join(d, "runs")
            with FakeHost() as host:
                p1 = self.run_suite(runs, host, "--manifest", mpath)
                self.assertEqual(p1.returncode, 0, p1.stderr)
                run_dir = p1.stdout.strip().splitlines()[-1]
                self.assertTrue(os.path.basename(run_dir).endswith("_refseq97-b7214d29"), run_dir)
                for name in ("run.json", "executions.jsonl", "compare.json", "report.md", "orchestrator.json",
                             "manifest.json"):
                    self.assertTrue(os.path.exists(os.path.join(run_dir, name)), name)
                rows, _ = sc.read_jsonl(os.path.join(run_dir, "executions.jsonl"))
                self.assertEqual(len(rows), 9)
                self.assertTrue(all(r["arm"] == "direct-remote" for r in rows))
                self.assertTrue(rows[0]["raw"].startswith("arms/direct-remote/raw/"))
                n_requests = len(host.log)
                # --resume: nothing new runs (only the capabilities GETs), the same directory is rewritten
                p2 = self.run_suite(runs, host, "--resume", run_dir)
                self.assertEqual(p2.returncode, 0, p2.stderr)
                self.assertEqual(p2.stdout.strip().splitlines()[-1], run_dir)
                self.assertEqual(len(host.log) - n_requests, 4)
                self.assertEqual(len(sc.read_jsonl(os.path.join(run_dir, "executions.jsonl"))[0]), 9)
                orch = sc.read_json(os.path.join(run_dir, "orchestrator.json"))
                self.assertEqual(len(orch["sessions"]), 2)
                # a new invocation is a new run directory
                time.sleep(1.1)
                p3 = self.run_suite(runs, host, "--manifest", mpath, "--limit", "1", "--max-requests", "3")
                self.assertEqual(p3.returncode, 0, p3.stderr)
                self.assertNotEqual(p3.stdout.strip().splitlines()[-1], run_dir)
            with open(os.path.join(run_dir, "report.md")) as fh:
                report = fh.read()
            self.assertIn("# Verdict", report)
            self.assertIn("the stub ran", report)

    def test_bad_options(self):
        with tempfile.TemporaryDirectory() as d, FakeHost() as host:
            p = self.run_suite(d, host, "--bogus")
            self.assertEqual(p.returncode, 2)
            p = self.run_suite(d, host, "--manifest", os.path.join(d, "missing.json"))
            self.assertEqual(p.returncode, 2)
            self.assertIn("no manifest", p.stderr)


if __name__ == "__main__":
    unittest.main()
