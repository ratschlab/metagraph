"""The resume rule of a warm-up and its target (lane V5): after an interruption between a warm-up and its target, the
warm-up runs again right before the target on every arm that runs the target (suite_common.paired_plan's rewarm,
suite_common.resume_rewarms in run_direct.py, the MCP runner on the ids of one call), as a line with the warm-up's id
and repeat, rerun_for_resume true and rerun_target; the earlier execution stays the one judged. merge_run puts the
re-run next to its target as timing evidence (warm_up_reruns, run.json reruns_for_resume); reevaluate, compare, the
report and the gate judge the original execution and never count the re-run as a second outcome, a duplicate, a
missing or an extra execution. Offline: the fake host, tests/fake_mcp.py and the mex side in a local directory."""

import gzip
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

import gate_policy  # noqa: E402
import run_direct  # noqa: E402
import suite_common as sc  # noqa: E402
from fake_host import FakeHost  # noqa: E402
from test_orchestrator import SuiteCase  # noqa: E402
from test_runner import make_manifest, scenario  # noqa: E402

RESOLVE = {"sequence": "ACGT" * 10, "labels": ["562"]}


def scenarios(warm_up_mcp=True):
    """R-0, then the warm-up R-w of the target R-t (the manifest keeps a warm-up right before its target)."""
    out = [scenario("R-0", "R", "POST /resolve", RESOLVE),
           scenario("R-w", "R", "POST /resolve", RESOLVE, warms="R-t", tags=["warm-up"]),
           scenario("R-t", "R", "POST /resolve", RESOLVE)]
    for s in out:
        s["mcp"] = {"tool": "traverse_resolve", "arguments": {}}
    if not warm_up_mcp:
        out[1]["mcp"] = {"not_implemented": "test"}
    return out


def line(sid, arm, repeat, **extra):
    return dict({"scenario": sid, "arm": arm, "repeat": repeat, "outcome": "pass", "http": 200}, **extra)


class Rule(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.d = self.tmp.name
        self.manifest = {"scenarios": scenarios()}

    def tearDown(self):
        self.tmp.cleanup()

    def record(self, arm, rows):
        os.makedirs(os.path.join(self.d, "arms", arm), exist_ok=True)
        for r in rows:
            sc.append_jsonl(os.path.join(self.d, "arms", arm, "executions.jsonl"), r)

    def plan(self, manifest=None, arms=sc.ARMS):
        return dict((i["scenario"], i) for i in sc.paired_plan(manifest or self.manifest, self.d, warm_every=0,
                                                                arms=arms))

    def test_interrupted_between_warm_up_and_target(self):
        for arm in sc.ARMS:
            self.record(arm, [line("R-0", arm, 0), line("R-w", arm, 0)])
        self.record("mcp", [line("R-0", "mcp", 1), line("R-w", "mcp", 1)])
        plan = self.plan()
        self.assertEqual(plan["R-t"]["rewarm"], {"direct-local": ["R-w"], "direct-remote": ["R-w"], "mcp": ["R-w"]})
        self.assertEqual(plan["R-w"]["rewarm"], {})
        self.assertEqual((plan["R-w"]["direct-local"], plan["R-w"]["mcp_cold"], plan["R-w"]["mcp_warm"]),
                         (False, False, []))

    def test_only_the_arms_that_run_the_target(self):
        # direct-local: the warm-up is not recorded (it runs as planned, right before its target); direct-remote:
        # the target is recorded; MCP: the target's cold is recorded, its warm repeat is left (its first step left)
        self.record("direct-local", [line("R-0", "direct-local", 0)])
        self.record("direct-remote", [line("R-w", "direct-remote", 0), line("R-t", "direct-remote", 0)])
        self.record("mcp", [line("R-w", "mcp", 0), line("R-w", "mcp", 1), line("R-t", "mcp", 0)])
        plan = self.plan()
        self.assertEqual(plan["R-t"]["rewarm"], {"mcp": ["R-w"]})
        self.assertEqual(plan["R-t"]["mcp_warm"], [1])
        # a warm-up marked not implemented via MCP sends nothing there: never run again on that arm
        manifest = {"scenarios": scenarios(warm_up_mcp=False)}
        self.record("mcp", [line("R-t", "mcp", 1)])
        self.record("direct-local", [line("R-w", "direct-local", 0)])
        self.assertEqual(self.plan(manifest)["R-t"]["rewarm"], {"direct-local": ["R-w"]})

    def test_a_rerun_line_is_no_recorded_execution(self):
        # a re-run of a key whose own execution is missing (a torn file) does not stand in for it
        self.record("direct-local", [line("R-w", "direct-local", 0, rerun_for_resume=True, rerun_target="R-t")])
        plan = self.plan(arms=("direct-local",))
        self.assertTrue(plan["R-w"]["direct-local"])
        self.assertEqual(plan["R-t"]["rewarm"], {})

    def test_paired_plan_cli_fields(self):
        for arm in sc.ARMS:
            self.record(arm, [line("R-0", arm, 0), line("R-w", arm, 0)])
        self.record("mcp", [line("R-0", "mcp", 1), line("R-w", "mcp", 1)])
        mpath = os.path.join(self.d, "m.json")
        make_manifest(mpath, scenarios())
        out = subprocess.run([sys.executable, os.path.join(SUITE, "suite_common.py"), "paired-plan", mpath, self.d, "",
                              "", "0", "direct-local,direct-remote,mcp"], stdout=subprocess.PIPE,
                             universal_newlines=True, check=True).stdout.splitlines()
        self.assertEqual(out, ["R-0 0 0 0 - - - -", "R-w 0 0 0 - - - -", "R-t 1 1 1 1 R-w R-w R-w"])

    def test_resume_rewarms_within_a_scope(self):
        planned, _ = sc.plan_executions(list(enumerate(scenarios())), "direct-local", 0)
        self.assertEqual(sc.resume_rewarms(planned, {("R-0", 0), ("R-w", 0)}), {"R-t": [(1, scenarios()[1])]})
        self.assertEqual(sc.resume_rewarms(planned, {("R-0", 0)}), {})                  # the warm-up runs as planned
        self.assertEqual(sc.resume_rewarms(planned, {("R-w", 0), ("R-t", 0)}), {})      # the target is recorded
        # a paired step that lists the target alone (--only R-t) never runs its warm-up again
        self.assertEqual(sc.resume_rewarms([p for p in planned if p[1]["id"] == "R-t"], {("R-w", 0)}), {})


class DirectRunner(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.d = self.tmp.name
        self.mpath = os.path.join(self.d, "m.json")
        make_manifest(self.mpath, scenarios())

    def tearDown(self):
        self.tmp.cleanup()
        signal.signal(signal.SIGINT, signal.default_int_handler)
        signal.signal(signal.SIGTERM, signal.SIG_DFL)

    def run_direct(self, base, out, *extra):
        run_direct.reset_stop()
        return run_direct.main(["--manifest", self.mpath, "--arm", "direct-local", "--base", base, "--out", out,
                                "--pause", "0", "--quiet", "--checks", "stub", "--warm-every", "0"] + list(extra))

    def test_capped_between_warm_up_and_target_then_resumed(self):
        out = os.path.join(self.d, "arm")
        with FakeHost() as host:
            # 3 capabilities + R-0 + R-w + the closing GET: the cap falls between the warm-up and its target
            self.assertEqual(self.run_direct(host.base, out, "--max-requests", "6"), run_direct.EXIT_REQUEST_CAP)
            self.assertEqual([r["scenario"] for r in sc.read_jsonl(os.path.join(out, "executions.jsonl"))[0]],
                             ["R-0", "R-w"])
            n = len(host.log)
            self.assertEqual(self.run_direct(host.base, out), 0)
            # the session's capabilities, the warm-up again, the target, the closing GET
            self.assertEqual([p for _, p in host.log[n:] if p == "/resolve"], ["/resolve", "/resolve"])
            self.assertEqual(host.log[-1][1], "/capabilities")
            again = len(host.log)
            self.assertEqual(self.run_direct(host.base, out), 0)     # done: nothing runs again
            self.assertNotIn("/resolve", [p for _, p in host.log[again:]])
        rows, _ = sc.read_jsonl(os.path.join(out, "executions.jsonl"))
        self.assertEqual([(r["scenario"], r["repeat"], r.get("rerun_for_resume"), r.get("rerun_target"))
                          for r in rows],
                         [("R-0", 0, None, None), ("R-w", 0, None, None), ("R-w", 0, True, "R-t"),
                          ("R-t", 0, None, None)])
        first, rerun = rows[1], rows[2]
        self.assertNotEqual(first["raw"], rerun["raw"])
        self.assertTrue(os.path.exists(os.path.join(out, first["raw"])))
        with gzip.open(os.path.join(out, rerun["raw"]), "rt") as fh:
            doc = json.load(fh)
        self.assertEqual((doc["rerun_for_resume"], doc["rerun_target"], doc["scenario"]), (True, "R-t", "R-w"))
        self.assertLessEqual(sc.parse_utc(rerun["ended_utc"]), sc.parse_utc(rows[3]["started_utc"]))
        run = sc.read_json(os.path.join(out, "run.json"))
        self.assertEqual(run["state"], "done")
        self.assertEqual(run["sessions"][1]["reruns_for_resume"], [{"scenario": "R-w", "target": "R-t"}])
        self.assertNotIn("reruns_for_resume", run["sessions"][2])

    def test_a_step_of_the_target_alone_does_not_rerun(self):
        out = os.path.join(self.d, "arm")
        with FakeHost() as host:
            self.assertEqual(self.run_direct(host.base, out, "--only", "R-0,R-w"), 0)
            n = len(host.log)
            self.assertEqual(self.run_direct(host.base, out, "--only", "R-t"), 0)
            self.assertEqual([p for _, p in host.log[n:]], ["/resolve"])
        self.assertFalse([r for r in sc.read_jsonl(os.path.join(out, "executions.jsonl"))[0] if sc.is_rerun(r)])


class Orchestrated(SuiteCase):
    """run_suite.sh interrupted (Ctrl-C) while the MCP arm waits on the target's cold step, after the warm-up ran on
    every arm, then resumed: the warm-up runs again before the target on direct-local, direct-remote and the MCP arm;
    the merge, the comparison, the report, the re-evaluation and the gate judge the warm-up's first execution."""

    def policy(self):
        pol = gate_policy.build_policy(sc.read_json(self.mpath), routes=["POST /resolve"],
                                       arms=["direct-local", "mcp"], table_path="", warm_every=0)
        path = os.path.join(self.d, "policy.json")
        sc.write_json_atomic(path, pol)
        return path

    def calls(self, run_dir):
        with open(os.path.join(run_dir, "arms", "mcp", "calls.log")) as fh:
            return [x.split()[:2] for x in fh if x.split() and x.split()[0][:2] == "R-"]

    def check_merged(self, run_dir, suffix=""):
        rows, bad = sc.read_jsonl(os.path.join(run_dir, "executions%s.jsonl" % suffix))
        self.assertEqual(bad, 0)
        keys = [sc.execution_key(r) for r in rows]
        self.assertEqual(len(keys), len(set(keys)))                  # one line per key, no re-run among them
        self.assertFalse([r for r in rows if sc.is_rerun(r)])
        self.assertEqual(len(rows), 2 * 3 + 3 * 2)                  # direct 1 per scenario, MCP 2
        key = dict((sc.execution_key(r), r) for r in rows)
        t_mcp0 = sc.parse_utc(key[("R-t", "mcp", 0)]["started_utc"])
        for arm in sc.ARMS:
            # the warm-up's judged execution is its first, of the interrupted invocation
            self.assertLess(sc.parse_utc(key[("R-w", arm, 0)]["started_utc"]), t_mcp0 - 0.5, arm)
        run = sc.read_json(os.path.join(run_dir, "run%s.json" % suffix))
        again = run["reruns_for_resume"]
        self.assertEqual(sorted((r["arm"], r["scenario"], r["repeat"], r["rerun_target"], tuple(r["target_execution"]))
                                for r in again),
                         [("direct-local", "R-w", 0, "R-t", ("R-t", "direct-local", 0)),
                          ("direct-remote", "R-w", 0, "R-t", ("R-t", "direct-remote", 0)),
                          ("mcp", "R-w", 0, "R-t", ("R-t", "mcp", 0))])
        self.assertEqual(run["counts"]["reruns_for_resume"], 3)
        self.assertEqual(run["counts"]["executions"], 12)
        for arm in sc.ARMS:
            near = key[("R-t", arm, 0)]["warm_up_reruns"]
            self.assertEqual([(w["scenario"], w["repeat"]) for w in near], [("R-w", 0)], arm)
            self.assertGreaterEqual(near[0]["before_s"], 0.0)
            self.assertEqual(run["arms"][arm]["executions_recorded"], 3 if arm != "mcp" else 6)
            self.assertEqual(run["arms"][arm]["reruns_for_resume"], 1)
        self.assertNotIn("warm_up_reruns", key[("R-t", "mcp", 1)])
        return run

    def test_ctrl_c_between_warm_up_and_target_then_resume(self):
        make_manifest(self.mpath, scenarios())
        pol = self.policy()
        with FakeHost() as host:
            cmd = [os.path.join(SUITE, "run_suite.sh"), "--no-pdf", "--pause", "0", "--manifest", self.mpath,
                   "--mcp-cmd", self.mcp_cmd(), "--warm-every", "0", "--gate", pol, "--by-design", ""]
            proc = subprocess.Popen(cmd, env=self.env(host, FAKE_MCP_MODE="sigint", FAKE_MCP_SIGINT_ID="R-t"),
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
                                    start_new_session=True, stdin=subprocess.DEVNULL)
            deadline = time.time() + 120
            run_dir = None
            while time.time() < deadline:
                runs = os.listdir(self.runs) if os.path.isdir(self.runs) else []
                if runs:
                    run_dir = os.path.join(self.runs, runs[0])
                    if os.path.exists(os.path.join(run_dir, "arms", "mcp", "calls.log")) and \
                            ["R-t", "0"] in self.calls(run_dir):
                        break
                time.sleep(0.05)
            time.sleep(0.3)
            os.killpg(proc.pid, signal.SIGINT)
            out, err = proc.communicate(timeout=120)
            self.assertEqual(proc.returncode, 130, err)
            for arm in sc.ARMS:
                rows, _ = sc.read_jsonl(os.path.join(run_dir, "arms", arm, "executions.jsonl"))
                self.assertIn(("R-w", arm, 0), set(sc.execution_key(r) for r in rows), arm)
                self.assertNotIn("R-t", set(r["scenario"] for r in rows), arm)
            p = self.run_suite(host, "--resume", run_dir, "--by-design", "")
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertIn("GATE PASS", p.stderr)
            self.assertIn("mcp: R-w again before R-t (resume)", p.stderr)
            self.assertIn("direct-local: R-t (R-w again before it: resume)", p.stderr)
            self.assertIn("direct-remote: R-t (R-w again before it: resume)", p.stderr)
        # the MCP runner got the warm-up with the target in the target's cold step, the target alone after it
        self.assertEqual(self.calls(run_dir)[-2:], [["R-w,R-t", "0"], ["R-t", "1"]])
        for arm in sc.ARMS:
            rows, _ = sc.read_jsonl(os.path.join(run_dir, "arms", arm, "executions.jsonl"))
            seq = [(r["scenario"], r.get("rerun_target")) for r in rows if r["repeat"] == 0]
            self.assertEqual(seq[-2:], [("R-w", "R-t"), ("R-t", None)], arm)
        self.check_merged(run_dir)
        g = sc.read_json(os.path.join(run_dir, "gate.json"))
        self.assertTrue(g["passed"], g.get("problems"))
        self.assertFalse([x for x in g.get("problems") or [] if x.get("kind") == "missing"])
        comp = sc.read_json(os.path.join(run_dir, "compare.json"))
        self.assertEqual(comp["summary"]["totals"]["direct-local"]["executions"], 3)
        with open(os.path.join(run_dir, "report.md")) as fh:
            report = fh.read()
        self.assertIn("## Warm-ups run again for a resume", report)
        self.assertIn("| R-w | mcp | 0 |", report)
        self.assertNotIn("INCOMPLETE", report.split("# Environment")[0])
        # the offline path: the re-evaluated set carries the re-runs, its merge and gate judge the originals
        p = self.run_suite(host, "--reevaluate", run_dir, "--gate", pol, "--by-design", "")
        self.assertEqual(p.returncode, 0, p.stderr)
        for arm in sc.ARMS:
            rows, _ = sc.read_jsonl(os.path.join(run_dir, "arms", arm, "executions.reevaluated.jsonl"))
            self.assertEqual(len([r for r in rows if sc.is_rerun(r)]), 1, arm)
        self.check_merged(run_dir, ".reevaluated")
        self.assertTrue(sc.read_json(os.path.join(run_dir, "gate.reevaluated.json"))["passed"])

    def test_sequential_resume_after_the_warm_up(self):
        # a run of the warm-up alone (R-w's selection does not take its target), then a resume that selects the
        # target: the sequential mode's arms each run the warm-up again before it (run_direct.py over its whole
        # selection, the MCP arm's cold pass naming both)
        make_manifest(self.mpath, scenarios())
        ids = os.path.join(self.d, "w.txt")
        with open(ids, "w") as fh:
            fh.write("R-w\n")
        tid = os.path.join(self.d, "t.txt")
        with open(tid, "w") as fh:
            fh.write("R-t\n")
        with FakeHost() as host:
            p = self.run_suite(host, "--sequential", "--manifest", self.mpath, "--mcp-cmd", self.mcp_cmd(),
                               "--warm-every", "0", "--ids", ids, "--checks", "stub")
            self.assertEqual(p.returncode, 0, p.stderr)
            run_dir = p.stdout.strip().splitlines()[-1]
            p = self.run_suite(host, "--resume", run_dir, "--ids", tid)
            self.assertEqual(p.returncode, 0, p.stderr)
        self.assertEqual(self.calls(run_dir), [["R-w", "0"], ["R-w", "1"], ["R-w,R-t", "0"], ["R-t", "1"]])
        rows, _ = sc.read_jsonl(os.path.join(run_dir, "executions.jsonl"))
        self.assertEqual(sorted(sc.execution_key(r) for r in rows),
                         sorted([(s_, a, 0) for s_ in ("R-w", "R-t") for a in sc.ARMS] +
                                [("R-w", "mcp", 1), ("R-t", "mcp", 1)]))
        run = sc.read_json(os.path.join(run_dir, "run.json"))
        self.assertEqual(sorted(r["arm"] for r in run["reruns_for_resume"] if r["target_execution"]), list(sc.ARMS))


if __name__ == "__main__":
    unittest.main()
