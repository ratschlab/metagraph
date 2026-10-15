"""The run's snapshot of what judges its arms (lane S7r, P6-end-to-end-03 of the two-day review): run_suite.sh copies
suite_common.SNAPSHOT_FILES into RUN_DIR/scripts once (SNAPSHOT.json with their sha256); direct-remote runs from it,
mex gets it, the MCP runner gets it as {checks_dir}; run_direct records per session the sha256 of the checks and
suite_common it imported; merge_run states every arm judged with other checks (an error with --strict-identity and for
the gate) and makes outcomes stricter with the snapshot's checks."""

import contextlib
import io
import os
import shutil
import subprocess
import sys
import tempfile
import time
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.dirname(HERE)
sys.path.insert(0, SUITE)
sys.path.insert(0, HERE)

import merge_run  # noqa: E402
import run_direct  # noqa: E402
import suite_common as sc  # noqa: E402
from fake_host import FakeHost  # noqa: E402
from test_compare_report import T0, line, write_arm  # noqa: E402
from test_orchestrator import SuiteCase, mcp_scenarios  # noqa: E402
from test_runner import basic_scenarios, make_manifest, scenario  # noqa: E402

# checks.py as a change landing in the tree during a run would leave it: every execution fails
FAILING_CHECKS = '''

def evaluate(scenario, answer, arm, previous=None, client=None):
    return [{"check": "status", "outcome": "fail", "detail": "a checks change landed during the run"}]


def scenario_outcome(results):
    return "fail"
'''


class SnapshotHelpersTest(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def test_take_and_state(self):
        src = os.path.join(self.d, "src")
        os.makedirs(src)
        for name in sc.SNAPSHOT_FILES:
            shutil.copy(os.path.join(SUITE, name), src)
        run_dir = os.path.join(self.d, "runs", "R1")
        doc = sc.take_snapshot(src, run_dir, "abc", 2)
        self.assertEqual(sorted(doc["files"]), sorted(sc.SNAPSHOT_FILES))
        self.assertEqual(doc["files"]["checks.py"], sc.file_sha256(os.path.join(SUITE, "checks.py")))
        self.assertEqual((doc["scripts_commit"], doc["sessions_before"]), ("abc", 2))
        self.assertEqual(sc.snapshot_state(run_dir, src)[1:], ([], []))
        with open(os.path.join(src, "checks.py"), "a") as fh:
            fh.write(FAILING_CHECKS)
        self.assertEqual(sc.snapshot_state(run_dir, src)[1:], ([], ["checks.py"]))
        with open(os.path.join(run_dir, "scripts", "suite_common.py"), "a") as fh:
            fh.write("\n")
        self.assertEqual(sc.snapshot_state(run_dir, src)[1], ["suite_common.py"])
        self.assertEqual(sc.snapshot_state(os.path.join(self.d, "none")), (None, [], []))
        # the CLI takes it once, then keeps it
        other = os.path.join(self.d, "runs", "R2")
        out = subprocess.run([sys.executable, os.path.join(SUITE, "suite_common.py"), "snapshot", other, src, "c", "0"],
                             stdout=subprocess.PIPE, universal_newlines=True).stdout
        self.assertEqual(out, "created\n")
        out = subprocess.run([sys.executable, os.path.join(SUITE, "suite_common.py"), "snapshot", other, SUITE, "c",
                              "0"], stdout=subprocess.PIPE, universal_newlines=True).stdout
        self.assertEqual(out, "kept\nlive-differs checks.py\n")

    def test_snapshot_pins(self):
        self.assertTrue(sc.snapshot_pins("/a/runs/R1/scripts", "/b/copy/R1"))
        self.assertTrue(sc.snapshot_pins("/a/runs/R1/scripts/", "/a/runs/R1"))
        self.assertFalse(sc.snapshot_pins("/a/metagraph/scripts/suite", "/a/runs/R1"))
        self.assertFalse(sc.snapshot_pins("/a/runs/R2/scripts", "/a/runs/R1"))
        self.assertFalse(sc.snapshot_pins(None, "/a/runs/R1"))

    def test_run_direct_records_the_checks_it_imported(self):
        mpath = os.path.join(self.d, "m.json")
        make_manifest(mpath, basic_scenarios()[:1])
        out = os.path.join(self.d, "run")
        with FakeHost() as host:
            self.assertEqual(run_direct.main(["--manifest", mpath, "--arm", "direct-local", "--base", host.base,
                                              "--out", out, "--pause", "0", "--quiet", "--checks", "checks",
                                              "--warm-every", "0"]), 0)
        run = sc.read_json(os.path.join(out, "run.json"))
        want = (sc.file_sha256(os.path.join(SUITE, "checks.py")),
                sc.file_sha256(os.path.join(SUITE, "suite_common.py")), SUITE)
        sess = run["sessions"][0]
        self.assertEqual((sess["checks_sha256"], sess["suite_common_sha256"], sess["checks_dir"]), want)
        self.assertEqual((run["checks_sha256"], run["suite_common_sha256"], run["checks_dir"]), want)


class MergeChecksIdentityTest(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.run_dir = os.path.join(self.d, "runs", "RUN1")

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def make(self, snapshot=True, direct_sha="snapshot", mcp_dir="snapshot", started=sc.iso_utc(T0 + 3600)):
        """A run whose direct-local sessions state |direct_sha| ('snapshot': the snapshot's) and whose MCP arm
        records |mcp_dir| ('snapshot': the run's)."""
        doc = sc.take_snapshot(SUITE, self.run_dir, "c", 0, now=T0) if snapshot else {"files": {}}
        files = doc["files"]
        sha = files.get("checks.py") if direct_sha == "snapshot" else direct_sha
        sessions = [{"started_utc": started, "checks_sha256": sha,
                     "suite_common_sha256": files.get("suite_common.py") if direct_sha == "snapshot" else sha}]
        write_arm(self.run_dir, "direct-local", [line("P-1", "direct-local", 0, T0 + 10, 5.0)],
                  run_extra={"sessions": sessions})
        mdir = os.path.join(self.run_dir, "scripts") if mcp_dir == "snapshot" else mcp_dir
        write_arm(self.run_dir, "mcp", [line("P-1", "mcp", 0, T0 + 20, 5.0)],
                  run_extra={"checks_dir": mdir} if mdir else {})
        return self.run_dir

    def strict(self, d):
        err = io.StringIO()
        with contextlib.redirect_stderr(err), contextlib.redirect_stdout(io.StringIO()):
            rc = merge_run.main([d, "--strict-identity"])
        return rc, err.getvalue()

    def test_a_run_judged_with_its_snapshot_is_quiet(self):
        merged = merge_run.merge(self.make())
        ident = merged["checks_identity"]
        self.assertEqual((ident["differs"], ident["unverified"], ident["errors"]), ([], [], []))
        self.assertEqual(ident["arms"]["mcp"]["checks_dirs"], [os.path.join(self.run_dir, "scripts")])
        self.assertFalse([w for w in merged["warnings"] if "checks" in w])
        self.assertNotIn("checks identity", self.strict(self.run_dir)[1])   # (the synthetic run has no manifest)

    def test_an_arm_judged_with_other_checks_is_an_error(self):
        # direct-local imported another checks.py; the MCP runner judged with the live tree, not the snapshot
        merged = merge_run.merge(self.make(direct_sha="0" * 64, mcp_dir="/somewhere/metagraph/scripts/suite"))
        errors = merged["checks_identity"]["errors"]
        self.assertTrue([e for e in errors if e.startswith("arms/direct-local: 1 of 1 session(s) judged with checks "
                                                           "sha256 000000000000 (sessions 0), the run's snapshot's")])
        self.assertTrue([e for e in errors if e.startswith("arms/mcp: 1 of 1 invocation(s) judged with other checks")])
        self.assertEqual(len([w for w in merged["warnings"] if "judged with" in w]), 3)
        rc, err = self.strict(self.run_dir)
        self.assertEqual(rc, 3)
        self.assertIn("checks identity: arms/mcp", err)

    def test_an_arm_that_states_nothing_is_unverified(self):
        # a session before the field, an MCP runner that records no checks directory: unverified, an error where
        # the run has a snapshot, nothing without one (a run of before the snapshot)
        merged = merge_run.merge(self.make(direct_sha=None, mcp_dir=None, started=sc.iso_utc(T0 - 3600)))
        ident = merged["checks_identity"]
        self.assertEqual(ident["unverified"], [
            "arms/direct-local: 1 of 1 session(s) state no checks sha256 (1 ran before the run's snapshot was taken)",
            "arms/mcp: 1 of 1 invocation(s) state no checks directory or sha256"])
        self.assertEqual(ident["errors"], ident["unverified"])
        rc, err = self.strict(self.run_dir)
        self.assertEqual(rc, 3)
        self.assertIn("checks identity: arms/mcp: 1 of 1 invocation(s) state no checks directory", err)
        other = os.path.join(self.d, "old")
        self.run_dir = other
        merged = merge_run.merge(self.make(snapshot=False, direct_sha=None, mcp_dir=None))
        self.assertEqual((merged["checks_identity"]["snapshot"], merged["checks_identity"]["errors"]), (None, []))
        self.assertNotIn("checks identity", self.strict(other)[1])

    def test_without_a_snapshot_arms_stating_different_checks_are_named(self):
        merged = merge_run.merge(self.make(snapshot=False, direct_sha="a" * 64, mcp_dir=None))
        self.assertEqual(merged["checks_identity"]["differs"], [])
        run = sc.read_json(os.path.join(self.run_dir, "arms", "mcp", "run.json"))
        run["checks_sha256"] = "b" * 64        # a runner without invocations states its one at the top
        sc.write_json_atomic(os.path.join(self.run_dir, "arms", "mcp", "run.json"), run)
        merged = merge_run.merge(self.run_dir)
        self.assertTrue(merged["checks_identity"]["differs"][0].startswith(
            "the arms were judged with different checks"))

    def test_the_mcp_runner_states_its_checks_per_invocation(self):
        # async's runner records checks_dir, checks_sha256 and suite_common_sha256 per invocation (run.json
        # invocations[]); an invocation whose suite_common is not the snapshot's is named, and a top-level
        # checks_module is not read once invocations exist
        d = self.make()
        files = sc.read_json(os.path.join(d, "scripts", "SNAPSHOT.json"))["files"]
        snap = os.path.join(d, "scripts")
        good = {"checks_dir": snap, "checks_sha256": files["checks.py"],
                "suite_common_sha256": files["suite_common.py"]}
        run = sc.read_json(os.path.join(d, "arms", "mcp", "run.json"))
        run.update(invocations=[good, dict(good, suite_common_sha256="c" * 64)],
                   checks_module={"file": "x", "stub": False, "sha256": "d" * 64})
        sc.write_json_atomic(os.path.join(d, "arms", "mcp", "run.json"), run)
        ident = merge_run.merge(d)["checks_identity"]
        self.assertEqual(ident["errors"], ["arms/mcp: 1 of 2 invocation(s) judged with other checks than the run's "
                                           "snapshot (e.g. checks_sha256 %s, suite_common_sha256 %s, checks_dir %s)"
                                           % (files["checks.py"], "c" * 64, snap)])
        self.assertEqual(ident["arms"]["mcp"]["suite_common_sha256"], sorted([files["suite_common.py"], "c" * 64]))
        run["invocations"] = [good, good]
        sc.write_json_atomic(os.path.join(d, "arms", "mcp", "run.json"), run)
        self.assertEqual(merge_run.merge(d)["checks_identity"]["errors"], [])

    def test_the_merge_judges_with_the_snapshot(self):
        # checks.scenario_outcome makes outcomes stricter: the run's snapshot's, not the live tree's
        d = self.make()
        with open(os.path.join(d, "scripts", "checks.py"), "a") as fh:
            fh.write(FAILING_CHECKS)
        merged = merge_run.merge(d)
        rows, _ = sc.read_jsonl(os.path.join(d, "executions.jsonl"))
        self.assertEqual(set(r["outcome"] for r in rows), {"fail"})
        self.assertEqual(merged["checks_identity"]["differs"],
                         ["the run's snapshot changed after it was taken: scripts/checks.py"])


class RunSuiteSnapshotTest(SuiteCase):
    def test_a_checks_change_during_the_run_does_not_reach_it(self):
        # a copy of the suite runs a paired run; a checks change lands in the copy after the first scenario: every
        # later step still judges with the run's snapshot (before it, direct-remote's later steps imported the
        # changed checks.py), each arm records the snapshot's sha256 and the merge finds nothing to state; a resume
        # says the tree differs and keeps the snapshot
        suite = os.path.join(self.d, "suite")
        os.makedirs(suite)
        for name in os.listdir(SUITE):
            if name.endswith((".py", ".sh", ".json")):
                shutil.copy2(os.path.join(SUITE, name), suite)
        make_manifest(self.mpath, [dict(s, mcp={"tool": "pattern_search", "arguments": {}}) for s in
                                   [scenario("S-%02d" % i, "S", "POST /slow", {"sleep_s": 0.2}) for i in range(4)]])
        cmd = [os.path.join(suite, "run_suite.sh"), "--no-pdf", "--pause", "0", "--manifest", self.mpath,
               "--warm-every", "0", "--checks", "checks", "--mcp-cmd", self.mcp_cmd()]
        with FakeHost() as host:
            proc = subprocess.Popen(cmd, env=self.env(host), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    universal_newlines=True, stdin=subprocess.DEVNULL)
            deadline = time.time() + 60
            while time.time() < deadline and sum(1 for _, x in host.log if x == "/slow") < 2:
                time.sleep(0.02)
            with open(os.path.join(suite, "checks.py"), "a") as fh:
                fh.write(FAILING_CHECKS)
            out, err = proc.communicate(timeout=180)
            self.assertEqual(proc.returncode, 0, err)
            run_dir = out.strip().splitlines()[-1]
            rows, _ = sc.read_jsonl(os.path.join(run_dir, "executions.jsonl"))
            self.assertEqual(len(rows), 4 * 4)
            self.assertEqual(set(r["outcome"] for r in rows), {"pass"})
            snap = sc.read_json(os.path.join(run_dir, "scripts", "SNAPSHOT.json"))
            self.assertEqual(snap["files"]["checks.py"], sc.file_sha256(os.path.join(SUITE, "checks.py")))
            merged = sc.read_json(os.path.join(run_dir, "run.json"))
            for arm in ("direct-local", "direct-remote"):
                self.assertEqual(set(x["checks_sha256"] for x in merged["arms"][arm]["sessions"]),
                                 {snap["files"]["checks.py"]}, arm)
            self.assertEqual(merged["checks_identity"]["errors"], [])
            self.assertEqual(merged["checks_identity"]["arms"]["mcp"]["checks_dirs"],
                             [os.path.join(run_dir, "scripts")])
            self.assertEqual(merged["checks_identity"]["arms"]["mcp"]["checks_sha256"], [snap["files"]["checks.py"]])
            sess = sc.read_json(os.path.join(run_dir, "orchestrator.json"))["sessions"][-1]
            self.assertEqual(sess["snapshot"], {"state": "created", "changed": [], "live_differs": [],
                                                "mcp_checks_dir": "snapshot"})
            p = subprocess.run([os.path.join(suite, "run_suite.sh"), "--no-pdf", "--resume", run_dir],
                               env=self.env(host), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               universal_newlines=True, stdin=subprocess.DEVNULL, timeout=180)
            self.assertEqual(p.returncode, 0, p.stderr)
        self.assertIn("checks.py of %s differ(s) from the run's snapshot" % suite, p.stderr)
        sess = sc.read_json(os.path.join(run_dir, "orchestrator.json"))["sessions"][-1]
        self.assertEqual((sess["snapshot"]["state"], sess["snapshot"]["live_differs"]), ("kept", ["checks.py"]))
        self.assertEqual(sc.read_json(os.path.join(run_dir, "run.json"))["checks_identity"]["errors"], [])


if __name__ == "__main__":
    unittest.main()
