"""Offline tests of the bench run (bench/run_bench.py): the mex side in a local directory (SUITE_MEX=local:DIR), both
direct arms against the local fake host, the MCP arm through tests/fake_mcp_bench.py, the sampler reading a fake
`docker` on PATH. Covered: the schedule executed one request at a time in its arm orders, the warm-up, the identity
probes, the {checks_dir} snapshot, the preflight refusal on a loaded host, --max-steps, Ctrl-C and --resume, the lock,
the MCP template, and the run's own analysis (classes, guard, verdict)."""

import contextlib
import io
import json
import os
import re
import shlex
import signal
import subprocess
import sys
import tempfile
import time
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.dirname(HERE)
BENCH = os.path.join(SUITE, "bench")
sys.path.insert(0, HERE)

from bench_fixtures import fake_panel, tail_keys  # noqa: E402
import analyze  # noqa: E402
import bench_lib as bl  # noqa: E402
import run_bench  # noqa: E402
import suite_common as sc  # noqa: E402
from fake_host import FakeHost  # noqa: E402

FAKE_DOCKER = r'''#!/bin/sh
# fake `docker stats --no-stream --format '{{json .}}' STAGING SRA` for the bench tests
sra="${FAKE_DOCKER_SRA:-0.50%}"
if [ -n "$FAKE_DOCKER_LOAD_FILE" ] && [ -f "$FAKE_DOCKER_LOAD_FILE" ]; then sra="7700.00%"; fi
if [ -n "$FAKE_DOCKER_DOWN_FILE" ] && [ -f "$FAKE_DOCKER_DOWN_FILE" ]; then echo "docker: daemon down" >&2; exit 1; fi
rest='"MemUsage":"1GiB / 128GiB","NetIO":"1kB / 2kB","BlockIO":"0B / 0B"'
echo "{\"Name\":\"metagraph_refseq33m_experimental\",\"CPUPerc\":\"1.00%\",$rest}"
echo "{\"Name\":\"metagraph_sra_all\",\"CPUPerc\":\"$sra\",$rest}"
'''


class BenchCase(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.d = self.tmp.name
        self.mex = os.path.join(self.d, "mex")
        os.makedirs(self.mex)
        self.runs = os.path.join(self.d, "runs")
        self.panel_path = os.path.join(self.d, "panel.json")
        self.panel = fake_panel(self.panel_path)
        self.bin = os.path.join(self.d, "bin")
        os.makedirs(self.bin)
        with open(os.path.join(self.bin, "docker"), "w") as fh:
            fh.write(FAKE_DOCKER)
        os.chmod(os.path.join(self.bin, "docker"), 0o755)

    def tearDown(self):
        self.tmp.cleanup()

    def mcp_cmd(self):
        return "%s %s --ids {ids} --repeat {repeat} --out {out} --manifest {manifest} --checks-dir {checks_dir}" % (
            shlex.quote(sys.executable), shlex.quote(os.path.join(HERE, "fake_mcp_bench.py")))

    def env(self, host, sampler=False, **extra):
        env = dict(os.environ)
        for name in ("METAGRAPH_ADMIN_KEY_FILE", "METAGRAPH_ADMIN_KEY"):
            env.pop(name, None)
        env.update({"SUITE_MCP_CMD": self.mcp_cmd(), "BENCH_RUNS": self.runs, "SUITE_REMOTE_BASE": host.base,
                    "SUITE_LOCAL_BASE": host.base, "SUITE_MEX": "local:" + self.mex, "SUITE_MEX_PY": sys.executable,
                    "PYTHON": sys.executable, "SUITE_SAMPLER": "1" if sampler else "0", "SUITE_WAIT_RUNNER_S": "4",
                    "PATH": self.bin + os.pathsep + env.get("PATH", ""),
                    # this test's runners only (they name its fake host's port), never another test's stand-in
                    "SUITE_RUNNER_PATTERN": r"[s]cripts/run_direct\.py .*--base %s" % re.escape(host.base)})
        env.update(extra)
        return env

    def bench(self, host, *args, **kw):
        cmd = [sys.executable, os.path.join(BENCH, "run_bench.py")] + list(args)
        return subprocess.run(cmd, env=self.env(host, **kw), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              universal_newlines=True, timeout=300, stdin=subprocess.DEVNULL)

    def run_dir(self, p):
        return p.stdout.strip().splitlines()[-1]


class RunBenchOffline(BenchCase):
    def test_schedule_executed_one_at_a_time(self):
        with FakeHost() as host:
            p = self.bench(host, "--panel", self.panel_path, "--seed", "5", "--cycles", "2", "--checks", "stub",
                           "--pause", "0", "--no-preflight")
            self.assertEqual(p.returncode, 0, p.stderr)
            log = list(host.log)
        run_dir = self.run_dir(p)
        sched = sc.read_json(os.path.join(run_dir, "schedule.json"))
        self.assertEqual((sched["seed"], len(sched["steps"])), (5, 6 + 6 * 2 * 3))
        rows = []
        for arm in bl.ARMS:
            rows += sc.read_jsonl(os.path.join(run_dir, "arms", arm, "executions.jsonl"))[0]
        self.assertEqual(len(rows), 42)
        by_arm = {}
        for r in rows:
            by_arm[r["arm"]] = by_arm.get(r["arm"], 0) + 1
        self.assertEqual(by_arm, {"direct-local": 18, "direct-remote": 12, "mcp": 12})
        self.assertTrue(all(r["repeat"] == 0 for r in rows))
        # in the schedule's order, one request at a time
        key = dict(((r["scenario"], r["arm"]), r) for r in rows)
        spans = [(sc.parse_utc(key[(s["bench_id"], s["arm"])]["started_utc"]),
                  sc.parse_utc(key[(s["bench_id"], s["arm"])]["ended_utc"])) for s in sched["steps"]]
        for (a0, a1), (b0, b1) in zip(spans, spans[1:]):
            self.assertLessEqual(a1, b0 + 0.005)
        # the warm-up through direct-local first
        self.assertTrue(all(s["arm"] == "direct-local" for s in sched["steps"][:6]))
        # the host saw: the release GET; direct-local 3 capabilities + 18 + 3 identity probes; direct-remote
        # 3 capabilities + 12 + its closing GET
        self.assertEqual(len(log), 1 + (3 + 18 + 3) + (3 + 12 + 1))
        probes = [r for r in sc.read_jsonl(os.path.join(run_dir, "steps.jsonl"))[0] if r.get("probe")]
        self.assertEqual([r["block"] for r in probes], [0, 1, 2])
        self.assertTrue(all(r["server_instance"] and r["release"] for r in probes))
        # the MCP runner judged with the run's snapshot
        with open(os.path.join(run_dir, "arms", "mcp", "calls.log")) as fh:
            calls = fh.read().splitlines()
        self.assertEqual(len(calls), 12)
        self.assertTrue(all(c.endswith("checks_dir=" + os.path.join(run_dir, "scripts")) for c in calls))
        bench = sc.read_json(os.path.join(run_dir, "bench.json"))
        self.assertEqual((bench["state"], bench["sessions"][-1]["state"], bench["seed"]), ("done", "done", 5))
        self.assertEqual(bench["mcp_cmd"], self.mcp_cmd())
        self.assertTrue(os.path.exists(os.path.join(run_dir, "scripts", "SNAPSHOT.json")))
        # the run's own analysis: no verdict without samples (unguarded), and none with 2 of 5 successes
        res = sc.read_json(os.path.join(run_dir, "analysis.json"))
        self.assertFalse(res["verdict"]["emitted"])
        self.assertEqual(res["load_guard"]["verdict"], "unguarded")
        self.assertEqual(res["classes"]["mcp"], {"success": 12})
        self.assertEqual(res["classes"]["direct-remote"], {"success": 12})
        self.assertFalse(tail_keys(res))  # no tail percentile anywhere
        self.assertFalse(os.path.exists(os.path.join(self.mex, "suite-dev", "LOCK")))

    def test_guarded_run_emits_a_verdict(self):
        with FakeHost() as host:
            # the laptop's own load stands in for mex's: the test lifts the load1 limit, the sra_all one stays
            p = self.bench(host, "--panel", self.panel_path, "--seed", "9", "--cycles", "2", "--checks", "stub",
                           "--pause", "0", "--preflight-s", "4", "--load1-max", "100000", sampler=True)
            self.assertEqual(p.returncode, 0, p.stderr)
        run_dir = self.run_dir(p)
        self.assertTrue(os.listdir(os.path.join(run_dir, "samples")))
        with contextlib.redirect_stdout(io.StringIO()):
            rc = analyze.main([run_dir, "--min-success", "2", "--boot", "200"])
        res = sc.read_json(os.path.join(run_dir, "analysis.json"))
        self.assertEqual(rc, 0, res["verdict"])
        self.assertTrue(res["verdict"]["emitted"])
        self.assertEqual(res["load_guard"]["verdict"], "clean")
        self.assertEqual(sorted(res["load_guard"]["blocks"]), ["0", "1", "2"])
        agg = res["aggregates"]["all"]["pairs"]["mcp/direct-local"]
        self.assertEqual(agg["queries"], 6)
        self.assertGreater(agg["geo_mean_ratio"], 0)
        self.assertEqual(len(res["intervals"]["pairs"]["mcp/direct-local"]["geo_mean_ratio"]), 2)
        with open(os.path.join(run_dir, "analysis.md")) as fh:
            self.assertIn("Performance verdict: emitted", fh.read())

    def test_preflight_refuses_a_loaded_host(self):
        with FakeHost() as host:
            p = self.bench(host, "--panel", self.panel_path, "--seed", "1", "--cycles", "2", "--checks", "stub",
                           "--pause", "0", "--preflight-s", "3", "--load1-max", "100000", sampler=True,
                           FAKE_DOCKER_SRA="7700.00%")
            self.assertEqual(p.returncode, 4, p.stderr)
            self.assertIn("preflight: the host is loaded", p.stderr)
            self.assertEqual(len(host.log), 1)  # the release GET only
        run_dir = self.run_dir(p)
        self.assertFalse(os.path.exists(os.path.join(run_dir, "arms", "direct-remote", "executions.jsonl")))
        self.assertFalse(os.path.exists(os.path.join(self.mex, "suite-dev", "LOCK")))
        bench = sc.read_json(os.path.join(run_dir, "bench.json"))
        self.assertEqual(bench["sessions"][-1]["state"], "refused")
        self.assertIn("preflight", bench["sessions"][-1]["reason"])
        self.assertFalse(os.path.exists(os.path.join(run_dir, "analysis.json")))

    def test_ctrl_c_stops_after_the_step_then_resume(self):
        with FakeHost() as host:
            cmd = [sys.executable, os.path.join(BENCH, "run_bench.py"), "--panel", self.panel_path, "--seed", "8",
                   "--cycles", "2", "--checks", "stub", "--pause", "0.3", "--no-preflight", "--no-analysis"]
            proc = subprocess.Popen(cmd, env=self.env(host), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    universal_newlines=True, stdin=subprocess.DEVNULL, start_new_session=True)
            steps = None
            deadline = time.time() + 120
            while time.time() < deadline and proc.poll() is None:
                runs = os.listdir(self.runs) if os.path.isdir(self.runs) else []
                if runs:
                    steps = os.path.join(self.runs, runs[0], "steps.jsonl")
                    if os.path.exists(steps) and len(sc.read_jsonl(steps)[0]) >= 4:
                        break
                time.sleep(0.05)
            proc.send_signal(signal.SIGINT)
            out, err = proc.communicate(timeout=120)
            self.assertEqual(proc.returncode, 130, err)
            run_dir = out.strip().splitlines()[-1]
            bench = sc.read_json(os.path.join(run_dir, "bench.json"))
            self.assertEqual(bench["sessions"][-1]["state"], "interrupted")
            done = [r for r in sc.read_jsonl(os.path.join(run_dir, "steps.jsonl"))[0] if not r.get("probe")]
            self.assertTrue(all(r["executed"] for r in done))  # the step in flight finished and was recorded
            self.assertLess(len(done), 42)
            self.assertFalse(os.path.exists(os.path.join(self.mex, "suite-dev", "LOCK")))
            q = self.bench(host, "--resume", run_dir, "--no-preflight", "--no-analysis")
            self.assertEqual(q.returncode, 0, q.stderr)
        rows = []
        for arm in bl.ARMS:
            rows += sc.read_jsonl(os.path.join(run_dir, "arms", arm, "executions.jsonl"))[0]
        self.assertEqual(len(rows), 42)

    def test_a_block_turning_loaded_stops_the_run(self):
        load = os.path.join(self.d, "logan-starts")
        with FakeHost() as host:
            cmd = [sys.executable, os.path.join(BENCH, "run_bench.py"), "--panel", self.panel_path, "--seed", "6",
                   "--cycles", "2", "--checks", "stub", "--pause", "0.2", "--preflight-s", "3", "--load1-max",
                   "100000", "--no-analysis"]
            proc = subprocess.Popen(cmd, env=self.env(host, sampler=True, FAKE_DOCKER_LOAD_FILE=load),
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
                                    stdin=subprocess.DEVNULL, start_new_session=True)
            deadline = time.time() + 120
            while time.time() < deadline and proc.poll() is None:
                runs = os.listdir(self.runs) if os.path.isdir(self.runs) else []
                steps = os.path.join(self.runs, runs[0], "steps.jsonl") if runs else None
                if steps and os.path.exists(steps) and len(sc.read_jsonl(steps)[0]) >= 8:
                    break
                time.sleep(0.05)
            with open(load, "w"):
                pass
            out, err = proc.communicate(timeout=180)
        self.assertEqual(proc.returncode, 4, err)
        self.assertIn("loaded during block", err)
        run_dir = out.strip().splitlines()[-1]
        bench = sc.read_json(os.path.join(run_dir, "bench.json"))
        self.assertEqual(bench["sessions"][-1]["state"], "stopped")
        self.assertTrue(bench["sessions"][-1]["reason"].startswith("loaded during block"))
        self.assertGreater(bench["sessions"][-1]["steps_missing"], 0)

    def start_and_wait(self, host, steps, **env):
        """run_bench in the background with the sampler; returns (proc, run_dir) once `steps` steps are logged."""
        cmd = [sys.executable, os.path.join(BENCH, "run_bench.py"), "--panel", self.panel_path, "--seed", "6",
               "--cycles", "2", "--checks", "stub", "--pause", "0.2", "--preflight-s", "3", "--load1-max",
               "100000", "--no-analysis"]
        proc = subprocess.Popen(cmd, env=self.env(host, sampler=True, **env), stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, universal_newlines=True, stdin=subprocess.DEVNULL,
                                start_new_session=True)
        deadline = time.time() + 120
        run_dir = None
        while time.time() < deadline and proc.poll() is None:
            runs = os.listdir(self.runs) if os.path.isdir(self.runs) else []
            run_dir = os.path.join(self.runs, runs[0]) if runs else None
            path = os.path.join(run_dir, "steps.jsonl") if run_dir else None
            if path and os.path.exists(path) and len(sc.read_jsonl(path)[0]) >= steps:
                break
            time.sleep(0.05)
        return proc, run_dir

    def test_a_sampler_that_dies_stops_the_run(self):
        # D7: the rest of the window would give no verdict; stop (exit 4, resumable)
        with FakeHost() as host:
            proc, run_dir = self.start_and_wait(host, 3)
            try:
                pid = int(sc.read_json(os.path.join(run_dir, "bench.json"))["sessions"][-1]["sampler_pid"])
                os.kill(pid, signal.SIGKILL)
            except (KeyError, TypeError, ValueError, OSError):
                proc.kill()
                raise
            out, err = proc.communicate(timeout=180)
        self.assertEqual(proc.returncode, 4, err)
        bench = sc.read_json(os.path.join(run_dir, "bench.json"))
        self.assertEqual(bench["sessions"][-1]["state"], "stopped")
        self.assertIn("the sampler stopped", bench["sessions"][-1]["reason"])
        self.assertGreater(bench["sessions"][-1]["steps_missing"], 0)

    def test_a_sampler_without_readings_stops_the_run(self):
        # D7: docker stops answering: two live checks in a row without a reading stop the run
        down = os.path.join(self.d, "docker-down")
        with FakeHost() as host:
            proc, run_dir = self.start_and_wait(host, 3, FAKE_DOCKER_DOWN_FILE=down)
            with open(down, "w"):
                pass
            out, err = proc.communicate(timeout=180)
        self.assertEqual(proc.returncode, 4, err)
        bench = sc.read_json(os.path.join(run_dir, "bench.json"))
        self.assertEqual(bench["sessions"][-1]["state"], "stopped")
        self.assertIn("the load is unknown", bench["sessions"][-1]["reason"])

    def test_max_steps_then_resume(self):
        with FakeHost() as host:
            p = self.bench(host, "--panel", self.panel_path, "--seed", "3", "--cycles", "2", "--checks", "stub",
                           "--pause", "0", "--no-preflight", "--max-steps", "10")
            self.assertEqual(p.returncode, 3, p.stderr)
            run_dir = self.run_dir(p)
            sched = sc.read_json(os.path.join(run_dir, "schedule.json"))
            q = self.bench(host, "--resume", run_dir, "--no-preflight")
            self.assertEqual(q.returncode, 0, q.stderr)
        self.assertEqual(sc.read_json(os.path.join(run_dir, "schedule.json")), sched)
        rows = []
        for arm in bl.ARMS:
            rows += sc.read_jsonl(os.path.join(run_dir, "arms", arm, "executions.jsonl"))[0]
        self.assertEqual(len(set((r["scenario"], r["arm"]) for r in rows)), 42)
        self.assertEqual(len(rows), 42)  # nothing ran twice
        bench = sc.read_json(os.path.join(run_dir, "bench.json"))
        self.assertEqual([s["state"] for s in bench["sessions"]], ["capped", "done"])
        steps = [r for r in sc.read_jsonl(os.path.join(run_dir, "steps.jsonl"))[0] if not r.get("probe")]
        self.assertEqual([sum(1 for r in steps if r["session"] == i) for i in (0, 1)], [10, 32])

    def test_mcp_error_is_counted_apart(self):
        with FakeHost() as host:
            bad = "R-q-0001~c2"
            p = self.bench(host, "--panel", self.panel_path, "--seed", "4", "--cycles", "2", "--checks", "stub",
                           "--pause", "0", "--no-preflight", FAKE_BENCH_MCP_ERROR=bad)
            self.assertEqual(p.returncode, 0, p.stderr)
        res = sc.read_json(os.path.join(self.run_dir(p), "analysis.json"))
        self.assertEqual(res["classes"]["mcp"], {"success": 11, "timeout": 1})
        others = res["queries"]["R-q-0001"]["arms"]["mcp"]["others"]
        self.assertEqual([(o["cycle"], o["class"]) for o in others], [(2, "timeout")])
        self.assertEqual(res["queries"]["R-q-0001"]["arms"]["mcp"]["n_success"], 1)
        self.assertEqual(res["queries"]["R-q-0001"]["pairs"]["mcp/direct-local"]["n"], 1)

    def test_the_suites_fake_mcp_runner_too(self):
        # tests/fake_mcp.py (run_suite.sh's tests' stand-in) writes no raw answer: its executions run and are
        # recorded, and the analysis cannot verify their result sizes (`unverified`, never a success)
        tmpl = "%s %s --ids {ids} --repeat {repeat} --out {out} --manifest {manifest} --checks-dir {checks_dir}" % (
            shlex.quote(sys.executable), shlex.quote(os.path.join(HERE, "fake_mcp.py")))
        with FakeHost() as host:
            p = self.bench(host, "--panel", self.panel_path, "--seed", "2", "--cycles", "1", "--checks", "stub",
                           "--pause", "0", "--no-preflight", SUITE_MCP_CMD=tmpl)
            self.assertEqual(p.returncode, 0, p.stderr)
        res = sc.read_json(os.path.join(self.run_dir(p), "analysis.json"))
        self.assertEqual(res["classes"]["mcp"], {"unverified": 6})
        self.assertEqual(res["classes"]["direct-local"], {"success": 6})
        self.assertIsNone(res["aggregates"]["all"]["pairs"]["mcp/direct-local"]["geo_mean_ratio"])

    def test_lock_held_refuses(self):
        os.makedirs(os.path.join(self.mex, "suite-dev"))
        with open(os.path.join(self.mex, "suite-dev", "LOCK"), "w") as fh:
            # a live lock of this host: this test's own process, taken after it started
            fh.write(json.dumps({"pid": os.getpid(), "host": sc.host_id(), "started_utc": "x",
                                 "started_epoch": time.time(), "run_id": "other", "mode": "paired",
                                 "token": "t"}) + "\n")
        with FakeHost() as host:
            p = self.bench(host, "--panel", self.panel_path, "--seed", "1", "--checks", "stub", "--no-preflight")
            self.assertEqual(p.returncode, 2, p.stderr)
            self.assertIn("another invocation drives the host", p.stderr)
            self.assertEqual(host.log, [])
        self.assertTrue(os.path.exists(os.path.join(self.mex, "suite-dev", "LOCK")))

    def test_mcp_template_refusals(self):
        with FakeHost() as host:
            p = self.bench(host, "--panel", self.panel_path, "--seed", "1", "--no-preflight",
                           SUITE_MCP_CMD="/no/such/python runner.py --ids {ids} --out {out}")
            self.assertEqual(p.returncode, 2)
            self.assertIn("interpreter /no/such/python is missing", p.stderr)
            p = self.bench(host, "--panel", self.panel_path, "--seed", "1", "--no-preflight",
                           SUITE_MCP_CMD=sys.executable + " runner.py --ids {ids} --out {out} --check-queue")
            self.assertEqual(p.returncode, 2)
            self.assertIn("METAGRAPH_ADMIN_KEY_FILE", p.stderr)
            self.assertEqual(host.log, [])

    def test_default_template_is_run_suite_s(self):
        with open(os.path.join(SUITE, "run_suite.sh")) as fh:
            text = fh.read()
        m = re.search(r'^ASYNC_MCP_CMD="([^"]+)"$', text, re.M)
        self.assertEqual(m.group(1), run_bench.DEFAULT_MCP_CMD)

    def test_dry_run_sends_nothing(self):
        with contextlib.redirect_stdout(io.StringIO()) as said:
            self.assertEqual(run_bench.main(["--panel", self.panel_path, "--dry-run", "--seed", "2"]), 0)
        self.assertIn("6 queries, 114 steps", said.getvalue())


if __name__ == "__main__":
    unittest.main()
