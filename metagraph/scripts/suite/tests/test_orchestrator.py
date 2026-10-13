"""Offline tests of run_suite.sh: the mex side runs in a local directory (SUITE_MEX=local:DIR), both direct arms
against the local fake host, the MCP arm through tests/fake_mcp.py. Covered: the paired order and host_cache,
the MCP command template (environment passed through, never recorded), the request cap and --resume, the lock
(held, stale, foreign), the refusal next to a running run_direct.py, the stub refusal, the samplers' files,
Ctrl-C and --resume, the sequential mode, and the runner's refusal of a directory another live runner holds.
RunSuiteOverSsh runs the ssh path through tests/fake_ssh.py: the per-invocation ssh control path."""

import json
import os
import platform
import shlex
import signal
import socket
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
from test_runner import basic_scenarios, make_manifest, scenario  # noqa: E402

SECRET = "s3cr3t-value-never-recorded"


def mcp_scenarios():
    """basic_scenarios with an MCP mapping, except P-2 (not implemented via MCP)."""
    out = []
    for s in basic_scenarios():
        if s["id"] != "P-2":
            s["mcp"] = {"tool": "pattern_search", "arguments": {}}
        out.append(s)
    return out


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


class SuiteCase(unittest.TestCase):
    """run_suite.sh with the mex side in a local directory and both direct arms against a local fake host."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.d = self.tmp.name
        self.mex = os.path.join(self.d, "mex")
        os.makedirs(self.mex)
        self.runs = os.path.join(self.d, "runs")
        self.mpath = os.path.join(self.d, "m.json")

    def tearDown(self):
        self.tmp.cleanup()

    def env(self, host, **extra):
        env = dict(os.environ)
        # never the real MCP runner nor a real key in a test
        for name in ("METAGRAPH_ADMIN_KEY_FILE", "METAGRAPH_ADMIN_KEY"):
            env.pop(name, None)
        env.update({"SUITE_MCP_CMD": "", "SUITE_RUNS": self.runs, "SUITE_REMOTE_BASE": host.base, "SUITE_LOCAL_BASE": host.base,
                    "SUITE_MEX": "local:" + self.mex, "SUITE_MEX_PY": sys.executable, "PYTHON": sys.executable,
                    "SUITE_SAMPLER": "0", "SUITE_WAIT_RUNNER_S": "4", "SUITE_POLL_S": "1",
                    "SUITE_TEST_SECRET": SECRET})
        env.update(extra)
        return env

    def run_suite(self, host, *args, **env_extra):
        cmd = [os.path.join(SUITE, "run_suite.sh"), "--no-pdf", "--pause", "0"] + list(args)
        return subprocess.run(cmd, env=self.env(host, **env_extra), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              universal_newlines=True, timeout=300, stdin=subprocess.DEVNULL)

    def mcp_cmd(self):
        return "%s %s --ids {ids} --repeat {repeat} --out {out} --manifest {manifest}" % (
            shlex.quote(sys.executable), shlex.quote(os.path.join(HERE, "fake_mcp.py")))

    def lock_path(self):
        return os.path.join(self.mex, "suite-dev", "LOCK")


class RunSuiteOffline(SuiteCase):
    def test_paired_order_host_cache_and_mcp_template(self):
        make_manifest(self.mpath, mcp_scenarios())
        with FakeHost() as host:
            p = self.run_suite(host, "--manifest", self.mpath, "--mcp-cmd", self.mcp_cmd())
            self.assertEqual(p.returncode, 0, p.stderr)
            # the release GET; per direct arm 3 capabilities + 9 executions + 1 closing GET
            self.assertEqual(len(host.log), 1 + 2 * 13)
        run_dir = p.stdout.strip().splitlines()[-1]
        rows, _ = sc.read_jsonl(os.path.join(run_dir, "executions.jsonl"))
        by_arm = {}
        for r in rows:
            by_arm.setdefault(r["arm"], []).append(r)
        self.assertEqual(dict((a, len(v)) for a, v in by_arm.items()),
                         {"direct-local": 9, "direct-remote": 9, "mcp": 12})
        key = dict((sc.execution_key(r), r) for r in rows)
        start = lambda k: sc.parse_utc(key[k]["started_utc"])  # noqa: E731
        for sid in ("C-1", "P-1", "R-1", "T-1", "D-1"):
            seq = [(sid, "mcp", 0), (sid, "direct-local", 0), (sid, "direct-remote", 0), (sid, "mcp", 1)]
            self.assertEqual(sorted(seq, key=start), seq, sid)
            self.assertEqual([key[k]["host_cache"] for k in seq], ["cold", "warm", "warm", "warm"])
        # one scenario after the other: D-1's MCP cold starts after T-1's MCP warm
        self.assertLess(start(("T-1", "mcp", 1)), start(("D-1", "mcp", 0)))
        self.assertEqual(key[("D-1", "mcp", 2)]["host_cache"], "warm")
        # P-2 is not implemented via MCP: the fake marks it, the host never saw it, direct-local was cold
        self.assertEqual(key[("P-2", "mcp", 0)]["host_cache"], "not_reached")
        self.assertEqual(key[("P-2", "direct-local", 0)]["host_cache"], "cold")
        self.assertNotIn(("P-2", "mcp", 1), key)
        orch = sc.read_json(os.path.join(run_dir, "orchestrator.json"))
        sess = orch["sessions"][-1]
        self.assertEqual((sess["mode"], sess["mcp_cmd"], sess["checks"]), ("paired", self.mcp_cmd(), "stub"))
        self.assertEqual(sess["loop"]["state"], "done")
        self.assertEqual(sess["loop"]["completed"], 6)
        self.assertEqual(sess["loop"]["direct_requests"], 26)
        self.assertIsNotNone(sess["mex_clock_offset_s"])
        run = sc.read_json(os.path.join(run_dir, "run.json"))
        self.assertEqual(run["arms"]["direct-local"]["state"], "done")
        self.assertEqual(run["arms"]["direct-remote"]["state"], "done")
        self.assertFalse(os.path.exists(self.lock_path()))
        with open(os.path.join(run_dir, "arms", "mcp", "calls.log")) as fh:
            calls = fh.read().splitlines()
        self.assertEqual(calls[:3], ["C-1 0 env=yes", "C-1 1 env=yes", "P-1 0 env=yes"])
        for name in ("orchestrator.json", "run.json", "compare.json", "report.md"):
            with open(os.path.join(run_dir, name)) as fh:
                self.assertNotIn(SECRET, fh.read(), name)
        self.assertNotIn(SECRET, p.stderr)
        with open(os.path.join(run_dir, "report.md")) as fh:
            report = fh.read()
        self.assertIn("paired: per scenario MCP cold", report)
        self.assertIn("| P | mcp | warm |", report)

    def test_cap_then_resume(self):
        make_manifest(self.mpath, basic_scenarios())
        with FakeHost() as host:
            p1 = self.run_suite(host, "--manifest", self.mpath, "--max-requests", "12")
            self.assertEqual(p1.returncode, 0, p1.stderr)
            # the release GET plus at most 12 runner requests, the two closing GETs included
            self.assertLessEqual(len(host.log), 1 + 12)
            run_dir = p1.stdout.strip().splitlines()[-1]
            sess = sc.read_json(os.path.join(run_dir, "orchestrator.json"))["sessions"][-1]
            self.assertEqual(sess["loop"]["state"], "capped")
            self.assertLessEqual(sess["loop"]["direct_requests"], 12)
            with open(os.path.join(run_dir, "report.md")) as fh:
                self.assertIn("INCOMPLETE", fh.read())
            p2 = self.run_suite(host, "--resume", run_dir)
            self.assertEqual(p2.returncode, 0, p2.stderr)
        self.assertEqual(p2.stdout.strip().splitlines()[-1], run_dir)
        rows, _ = sc.read_jsonl(os.path.join(run_dir, "executions.jsonl"))
        self.assertEqual(len(rows), 18)
        self.assertEqual(len(set(sc.execution_key(r) for r in rows)), 18)
        run = sc.read_json(os.path.join(run_dir, "run.json"))
        self.assertEqual((run["arms"]["direct-local"]["state"], run["arms"]["direct-remote"]["state"]),
                         ("done", "done"))
        orch = sc.read_json(os.path.join(run_dir, "orchestrator.json"))
        self.assertEqual([s["loop"]["state"] for s in orch["sessions"]], ["capped", "done"])

    def test_lock_held_stale_and_foreign(self):
        make_manifest(self.mpath, basic_scenarios()[:1])
        os.makedirs(os.path.join(self.mex, "suite-dev"))
        host_name = socket.gethostname()
        with FakeHost() as host:
            held = {"pid": os.getpid(), "host": host_name, "started_utc": "t", "started_epoch": time.time(),
                    "run_id": "other", "token": "x"}
            with open(self.lock_path(), "w") as fh:
                json.dump(held, fh)
            p = self.run_suite(host, "--manifest", self.mpath)
            self.assertEqual(p.returncode, 2, p.stderr)
            self.assertIn("held by pid %d" % os.getpid(), p.stderr)
            with open(self.lock_path()) as fh:
                self.assertEqual(json.load(fh)["token"], "x")  # untouched
            with open(self.lock_path(), "w") as fh:
                json.dump(dict(held, host="elsewhere"), fh)
            p = self.run_suite(host, "--manifest", self.mpath)
            self.assertEqual(p.returncode, 2, p.stderr)
            self.assertIn("on elsewhere", p.stderr)
            dead = subprocess.Popen([sys.executable, "-c", "pass"])
            dead.wait()
            with open(self.lock_path(), "w") as fh:
                json.dump(dict(held, pid=dead.pid), fh)
            p = self.run_suite(host, "--manifest", self.mpath)
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertIn("removing a stale lock", p.stderr)
        self.assertFalse(os.path.exists(self.lock_path()))

    def test_refuses_next_to_a_running_runner(self):
        make_manifest(self.mpath, basic_scenarios()[:1])
        other = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(60)", "scripts/run_direct.py"])
        try:
            with FakeHost() as host:
                p = self.run_suite(host, "--manifest", self.mpath)
                self.assertEqual(p.returncode, 2, p.stderr)
                self.assertIn("run_direct.py runs under", p.stderr)
                self.assertEqual(len(host.log), 1)  # the release GET only
        finally:
            other.kill()
            other.wait()
        self.assertFalse(os.path.exists(self.lock_path()))

    def test_stub_refused_beyond_status(self):
        scen = basic_scenarios()[:2]
        scen[1]["expect"].append({"check": "contains_source_context", "args": {"level": 0}})
        make_manifest(self.mpath, scen)
        with FakeHost() as host:
            p = self.run_suite(host, "--manifest", self.mpath)
            self.assertEqual(p.returncode, 2, p.stderr)
            self.assertIn("contains_source_context", p.stderr)
            self.assertIn("--allow-stub", p.stderr)
            p = self.run_suite(host, "--manifest", self.mpath, "--allow-stub", "--skip-local")
            self.assertEqual(p.returncode, 0, p.stderr)
        with open(os.path.join(p.stdout.strip().splitlines()[-1], "report.md")) as fh:
            self.assertIn("UNVERIFIED", fh.read())

    def test_samples_pulled_and_merged(self):
        make_manifest(self.mpath, basic_scenarios()[:2])
        with FakeHost() as host:
            p = self.run_suite(host, "--manifest", self.mpath, SUITE_SAMPLER="1")
            self.assertEqual(p.returncode, 0, p.stderr)
        run_dir = p.stdout.strip().splitlines()[-1]
        files = os.listdir(os.path.join(run_dir, "samples"))
        self.assertTrue(any(f.startswith("paired-") and f.endswith(".jsonl") for f in files), files)
        samples, _ = sc.read_jsonl(os.path.join(run_dir, "samples.jsonl"))
        self.assertTrue(samples)
        self.assertEqual(set(s["phase"] for s in samples), {"paired"})
        # the sampler was stopped
        pids = subprocess.run(["pgrep", "-f", os.path.join(self.mex, "suite-dev")], stdout=subprocess.PIPE,
                              universal_newlines=True).stdout.split()
        self.assertEqual(pids, [])

    def test_sigint_then_resume(self):
        scen = [scenario("S-%02d" % i, "S", "POST /slow", {"sleep_s": 0.3}) for i in range(12)]
        make_manifest(self.mpath, scen)
        with FakeHost() as host:
            cmd = [os.path.join(SUITE, "run_suite.sh"), "--no-pdf", "--pause", "0", "--manifest", self.mpath,
                   "--warm-every", "0"]
            proc = subprocess.Popen(cmd, env=self.env(host), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    universal_newlines=True, start_new_session=True, stdin=subprocess.DEVNULL)
            deadline = time.time() + 60
            while time.time() < deadline and sum(1 for m, p in host.log if p == "/slow") < 5:
                time.sleep(0.05)
            os.killpg(proc.pid, signal.SIGINT)  # a terminal's Ctrl-C reaches the whole group
            out, err = proc.communicate(timeout=120)
            self.assertEqual(proc.returncode, 130, err)
            self.assertIn("resume with:", err)
            self.assertFalse(os.path.exists(self.lock_path()))
            run_dir = [os.path.join(self.runs, x) for x in os.listdir(self.runs)][0]
            sess = sc.read_json(os.path.join(run_dir, "orchestrator.json"))["sessions"][-1]
            self.assertEqual(sess["loop"]["state"], "interrupted")
            p = self.run_suite(host, "--resume", run_dir, "--warm-every", "0")
            self.assertEqual(p.returncode, 0, p.stderr)
        rows, bad = sc.read_jsonl(os.path.join(run_dir, "executions.jsonl"))
        self.assertEqual(bad, 0)
        self.assertEqual(len(rows), 24)
        self.assertEqual(len(set(sc.execution_key(r) for r in rows)), 24)

    def test_sequential_remote_new_and_resume(self):
        make_manifest(self.mpath, basic_scenarios())
        with FakeHost() as host:
            p1 = self.run_suite(host, "--sequential", "--skip-local", "--manifest", self.mpath)
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
            # --resume keeps the run's mode: nothing new runs (only the capabilities GETs)
            p2 = self.run_suite(host, "--resume", run_dir, "--skip-local")
            self.assertEqual(p2.returncode, 0, p2.stderr)
            self.assertEqual(p2.stdout.strip().splitlines()[-1], run_dir)
            self.assertEqual(len(host.log) - n_requests, 4)
            orch = sc.read_json(os.path.join(run_dir, "orchestrator.json"))
            self.assertEqual([s["mode"] for s in orch["sessions"]], ["sequential", "sequential"])
            # sequential with the local arm on the (local) mex and a small cap over both arms
            time.sleep(1.1)
            p3 = self.run_suite(host, "--sequential", "--manifest", self.mpath, "--max-requests", "10",
                                "--warm-every", "0")
            self.assertEqual(p3.returncode, 0, p3.stderr)
            self.assertNotEqual(p3.stdout.strip().splitlines()[-1], run_dir)
        with open(os.path.join(run_dir, "report.md")) as fh:
            report = fh.read()
        self.assertIn("# Verdict", report)
        self.assertIn("the stub ran", report)

    def test_default_mcp_template(self):
        make_manifest(self.mpath, mcp_scenarios()[:2])
        with FakeHost() as host:
            # a full run takes SUITE_MCP_CMD's template without any flag; --resume keeps it
            p = self.run_suite(host, "--manifest", self.mpath, SUITE_MCP_CMD=self.mcp_cmd())
            self.assertEqual(p.returncode, 0, p.stderr)
            run_dir = p.stdout.strip().splitlines()[-1]
            rows, _ = sc.read_jsonl(os.path.join(run_dir, "executions.jsonl"))
            self.assertEqual(sum(1 for r in rows if r["arm"] == "mcp"), 4)
            p = self.run_suite(host, "--resume", run_dir)
            self.assertEqual(p.returncode, 0, p.stderr)
            sessions = sc.read_json(os.path.join(run_dir, "orchestrator.json"))["sessions"]
            self.assertEqual([x["mcp_cmd"] for x in sessions], [self.mcp_cmd()] * 2)
            # --no-mcp runs without it; a smoke run only with --mcp-cmd
            time.sleep(1.1)
            p = self.run_suite(host, "--manifest", self.mpath, "--no-mcp", SUITE_MCP_CMD=self.mcp_cmd())
            self.assertEqual(p.returncode, 0, p.stderr)
            rows, _ = sc.read_jsonl(os.path.join(p.stdout.strip().splitlines()[-1], "executions.jsonl"))
            self.assertFalse(any(r["arm"] == "mcp" for r in rows))
            # a template that checks the service queue refuses to start without the admin key
            p = self.run_suite(host, "--manifest", self.mpath, "--mcp-cmd", self.mcp_cmd() + " --check-queue")
            self.assertEqual(p.returncode, 2)
            self.assertIn("METAGRAPH_ADMIN_KEY_FILE", p.stderr)
            p = self.run_suite(host, "--manifest", self.mpath, "--mcp-cmd", "/no/such/python x {ids} {out}")
            self.assertIn("is missing", p.stderr)

    def test_bad_options(self):
        with FakeHost() as host:
            p = self.run_suite(host, "--bogus")
            self.assertEqual(p.returncode, 2)
            p = self.run_suite(host, "--manifest", os.path.join(self.d, "missing.json"))
            self.assertEqual(p.returncode, 2)
            self.assertIn("no manifest", p.stderr)
            p = self.run_suite(host, "--mcp-cmd", "x {ids}", "--mcp-run", self.d)
            self.assertIn("exclude each other", p.stderr)
            make_manifest(self.mpath, basic_scenarios()[:1])
            p = self.run_suite(host, "--manifest", self.mpath, "--mcp-cmd", "runner --ids {ids}")
            self.assertIn("{ids} and {out}", p.stderr)
            p = self.run_suite(host, "--manifest", self.mpath, "--mcp-run", self.d)
            self.assertIn("--sequential", p.stderr)
            self.assertEqual(p.returncode, 2)


class RunSuiteOverSsh(SuiteCase):
    """The ssh path: SUITE_MEX names a host and tests/fake_ssh.py stands in for ssh (the mex side runs in a
    local directory, its multiplexing as far as the control paths go)."""

    def setUp(self):
        super(RunSuiteOverSsh, self).setUp()
        self.bin = os.path.join(self.d, "bin")
        os.makedirs(self.bin)
        with open(os.path.join(self.bin, "ssh"), "w") as fh:
            fh.write("#!/bin/sh\nexec %s %s \"$@\"\n" % (shlex.quote(sys.executable),
                                                         shlex.quote(os.path.join(HERE, "fake_ssh.py"))))
        os.chmod(os.path.join(self.bin, "ssh"), 0o755)
        self.ssh_log = os.path.join(self.d, "ssh.log")
        self.ctl = os.path.join(self.d, "ctl")

    def env(self, host, **extra):
        env = super(RunSuiteOverSsh, self).env(host)
        env.update({"SUITE_MEX": "fake@mex", "PATH": self.bin + os.pathsep + env.get("PATH", ""),
                    "FAKE_SSH_HOME": self.mex, "FAKE_SSH_LOG": self.ssh_log, "SUITE_SSH_CONTROL_DIR": self.ctl})
        env.update(extra)
        return env

    def ssh_calls(self):
        with open(self.ssh_log) as fh:
            return [json.loads(x) for x in fh]

    def start_suite(self, host, *args, **env_extra):
        cmd = [os.path.join(SUITE, "run_suite.sh"), "--no-pdf", "--pause", "0"] + list(args)
        return subprocess.Popen(cmd, env=self.env(host, **env_extra), stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, universal_newlines=True, stdin=subprocess.DEVNULL)

    def wait_for(self, cond, timeout=60.0):
        deadline = time.time() + timeout
        while time.time() < deadline and not cond():
            time.sleep(0.05)
        self.assertTrue(cond(), "timed out")

    def test_a_refused_invocation_leaves_the_holders_connection(self):
        # review N1: the second invocation's cleanup closes its own ssh master only, never the one the
        # holder's step in flight runs through
        make_manifest(self.mpath, [scenario("S-01", "S", "POST /slow", {"sleep_s": 4})])
        with FakeHost() as host:
            holder = self.start_suite(host, "--manifest", self.mpath, "--warm-every", "0")
            self.wait_for(lambda: ("POST", "/slow") in host.log)      # the holder's direct-local step runs
            p = self.run_suite(host, "--manifest", self.mpath, "--warm-every", "0")
            self.assertEqual(p.returncode, 2, p.stderr)
            self.assertIn("held by pid %d" % holder.pid, p.stderr)
            out, err = holder.communicate(timeout=120)
            self.assertEqual(holder.returncode, 0, err)
        calls = self.ssh_calls()
        exits = [c for c in calls if c["op"] == "exit"]
        self.assertEqual([c["killed"] for c in exits], [[], []])       # no session was cut
        self.assertNotIn("direct-local: exit", err)                    # the holder's step did not fail
        own = os.path.join(self.ctl, "%%C-%d" % holder.pid)
        paths = set(c["control_path"] for c in calls if c["op"] == "run")
        self.assertIn(own, paths)
        self.assertEqual(len(paths), 2, paths)                         # one per invocation
        self.assertEqual(sorted(c["control_path"] for c in exits), sorted(paths))
        run_dir = out.strip().splitlines()[-1]
        rows, _ = sc.read_jsonl(os.path.join(run_dir, "executions.jsonl"))
        self.assertEqual(sorted((r["arm"], r["outcome"]) for r in rows),
                         [("direct-local", "pass"), ("direct-remote", "pass")])
        sess = sc.read_json(os.path.join(run_dir, "orchestrator.json"))["sessions"][-1]
        self.assertEqual(sess["loop"]["state"], "done")
        self.assertFalse(os.path.exists(self.lock_path()))

    def runner_gone(self, run_id):
        out = subprocess.run(["pgrep", "-f", "run_direct.py.*--run-id %s" % run_id], stdout=subprocess.PIPE,
                             universal_newlines=True).stdout
        return not out.strip()

    def test_mex_silent_after_a_cut_step_stops_the_loop(self):
        # review N3: the direct-local step's connection drops and mex stays silent: the wait for its runner must
        # not read the failing ssh as "no runner"; the loop stops with no further request (no direct-remote
        # step, no closing GETs), and --resume finishes the run once mex answers again
        make_manifest(self.mpath, [scenario("S-01", "S", "POST /slow", {"sleep_s": 3}),
                                   scenario("S-02", "S", "POST /slow", {"sleep_s": 0})])
        down = os.path.join(self.mex, ".ssh-down")
        with FakeHost() as host:
            proc = self.start_suite(host, "--manifest", self.mpath, "--warm-every", "0")
            self.wait_for(lambda: ("POST", "/slow") in host.log)
            open(down, "w").close()
            out, err = proc.communicate(timeout=120)
            self.assertEqual(proc.returncode, 1, err)
            self.assertIn("does not answer after the direct-local step of S-01 failed", err)
            self.assertIn("no closing steps", err)
            run_dir = out.strip().splitlines()[-1]
            run_id = os.path.basename(run_dir)
            self.wait_for(lambda: self.runner_gone(run_id))
            # the release GET, direct-local's three capabilities documents and its S-01 request: nothing else
            self.assertEqual(host.log, [("GET", "/capabilities"), ("GET", "/capabilities"),
                                        ("GET", "/pattern/capabilities"), ("GET", "/traverse/capabilities"),
                                        ("POST", "/slow")])
            sess = sc.read_json(os.path.join(run_dir, "orchestrator.json"))["sessions"][-1]
            self.assertEqual((sess["loop"]["state"], sess["loop"]["reason"]), ("failed", "mex_unreachable"))
            with open(os.path.join(run_dir, "report.md")) as fh:
                self.assertIn("INCOMPLETE", fh.read())
            self.assertTrue(os.path.exists(self.lock_path()))      # left behind: mex did not answer
            os.unlink(down)
            p = self.run_suite(host, "--resume", run_dir, "--warm-every", "0")
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertIn("removing a stale lock", p.stderr)
        rows, _ = sc.read_jsonl(os.path.join(run_dir, "executions.jsonl"))
        self.assertEqual(sorted(sc.execution_key(r) for r in rows),
                         [("S-01", "direct-local", 0), ("S-01", "direct-remote", 0),
                          ("S-02", "direct-local", 0), ("S-02", "direct-remote", 0)])
        self.assertFalse(os.path.exists(self.lock_path()))

    def test_no_closing_get_next_to_a_live_runner(self):
        # review N3: the step's connection drops, mex answers again, but the runner outlives the wait: the loop
        # stops and neither arm sends its closing GET /capabilities while the runner is on the host
        make_manifest(self.mpath, [scenario("S-01", "S", "POST /slow", {"sleep_s": 9})])
        down = os.path.join(self.mex, ".ssh-down")
        with FakeHost() as host:
            proc = self.start_suite(host, "--manifest", self.mpath, "--warm-every", "0")
            self.wait_for(lambda: ("POST", "/slow") in host.log)
            open(down, "w").close()
            self.wait_for(lambda: os.path.exists(self.ssh_log) and any(
                c["op"] == "dropped" for c in self.ssh_calls()))
            os.unlink(down)
            out, err = proc.communicate(timeout=120)
            self.assertEqual(proc.returncode, 1, err)
            self.assertIn("still runs on fake@mex", err)
            n = len(host.log)
            run_dir = out.strip().splitlines()[-1]
            self.wait_for(lambda: self.runner_gone(os.path.basename(run_dir)))
        self.assertEqual(n, 5)                                      # no direct-remote request, no closing GET
        sess = sc.read_json(os.path.join(run_dir, "orchestrator.json"))["sessions"][-1]
        self.assertEqual((sess["loop"]["state"], sess["loop"]["reason"]), ("failed", "runner_still_running"))
        self.assertFalse(os.path.exists(self.lock_path()))


if __name__ == "__main__":
    unittest.main()
