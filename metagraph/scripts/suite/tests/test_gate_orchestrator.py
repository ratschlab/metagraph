"""run_suite.sh and the gate (lane S6c): the sequential mode records every MCP repeat the paired plan has, the run's
plan is persisted before any arm runs, --ids FILE selects scenarios, and --gate POLICY makes the gate's verdict the
exit status (offline: the mex side in a local directory, the fake host and tests/fake_mcp.py)."""

import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.dirname(HERE)
sys.path.insert(0, SUITE)
sys.path.insert(0, HERE)

import gate_policy  # noqa: E402
import suite_common as sc  # noqa: E402
from fake_host import FakeHost  # noqa: E402
from test_orchestrator import SuiteCase, mcp_scenarios  # noqa: E402
from test_runner import make_manifest  # noqa: E402


def mcp_keys(run_dir):
    rows, _ = sc.read_jsonl(os.path.join(run_dir, "executions.jsonl"))
    return sorted((r["scenario"], r["repeat"]) for r in rows if r["arm"] == "mcp")


class SequentialMcpRepeats(SuiteCase):
    def test_sequential_runs_every_warm_repeat(self):
        # D-1 has repeat 3: the paired plan gives it MCP repeats 0, 1 and 2; the sequential mode must record all
        make_manifest(self.mpath, mcp_scenarios())
        with FakeHost() as host:
            p = self.run_suite(host, "--sequential", "--skip-local", "--skip-remote", "--manifest", self.mpath,
                               "--mcp-cmd", self.mcp_cmd(), "--checks", "stub", "--warm-every", "0")
            self.assertEqual(p.returncode, 0, p.stderr)
        run_dir = p.stdout.strip().splitlines()[-1]
        keys = mcp_keys(run_dir)
        self.assertIn(("D-1", 2), keys)
        self.assertEqual([k for k in keys if k[0] == "D-1"], [("D-1", 0), ("D-1", 1), ("D-1", 2)])
        # every implemented scenario has its warm repeat, P-2 (not implemented via MCP) only its mark
        self.assertEqual(sorted(set(k[0] for k in keys if k[1] == 1)), ["C-1", "D-1", "P-1", "R-1", "T-1"])
        self.assertNotIn(("P-2", 1), keys)
        with open(os.path.join(run_dir, "arms", "mcp", "calls.log")) as fh:
            calls = [x.split()[:2] for x in fh]
        # one pass per repeat index: cold, then warm 1 of five scenarios, then warm 2 of D-1 alone
        self.assertEqual([c[1] for c in calls], ["0", "1", "2"])
        self.assertEqual(calls[2][0], "D-1")



class IdsPlanAndGate(SuiteCase):
    def write(self, name, text):
        path = os.path.join(self.d, name)
        with open(path, "w") as fh:
            fh.write(text)
        return path

    def test_ids_selection_and_persisted_plan(self):
        make_manifest(self.mpath, mcp_scenarios())
        ids = self.write("ids.txt", "# two of them\nT-1\nR-1\n")
        with FakeHost() as host:
            p = self.run_suite(host, "--manifest", self.mpath, "--mcp-cmd", self.mcp_cmd(), "--checks", "stub",
                               "--ids", ids, "--warm-every", "0")
            self.assertEqual(p.returncode, 0, p.stderr)
            run_dir = p.stdout.strip().splitlines()[-1]
            rows, _ = sc.read_jsonl(os.path.join(run_dir, "executions.jsonl"))
            self.assertEqual(sorted(set(r["scenario"] for r in rows)), ["R-1", "T-1"])
            plan = sc.read_json(os.path.join(run_dir, "plan.json"))
            self.assertEqual(len(plan["plans"]), 1)
            self.assertEqual(plan["plans"][0]["selection"]["ids"], ["R-1", "T-1"])
            self.assertEqual(plan["plans"][0]["counts"]["executions"],
                             {"direct-local": 2, "direct-remote": 2, "mcp": 4})
            self.assertEqual(sorted(tuple(e) for e in plan["plans"][0]["executions"]),
                             sorted(sc.execution_key(r) for r in rows))
            sess = sc.read_json(os.path.join(run_dir, "orchestrator.json"))["sessions"][-1]
            self.assertEqual(sess["ids"], "R-1,T-1")
            # a resume keeps the ids (no wider plan), and adds no plan entry for the same selection
            p = self.run_suite(host, "--resume", run_dir)
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertEqual(len(sc.read_json(os.path.join(run_dir, "plan.json"))["plans"]), 1)
            with open(os.path.join(run_dir, "report.md")) as fh:
                self.assertNotIn("INCOMPLETE", fh.read().split("# Environment")[0])
            # an id the manifest lacks is refused before anything runs
            bad = self.write("bad.txt", "T-1 T-9\n")
            p = self.run_suite(host, "--manifest", self.mpath, "--ids", bad, "--no-mcp", "--checks", "stub")
            self.assertEqual(p.returncode, 2)
            self.assertIn("T-9", p.stderr)

    def test_plain_resume_keeps_the_union_of_the_selections(self):
        # a targeted resume (--ids) after a wider run, then a plain resume: the plain one continues the whole run
        make_manifest(self.mpath, mcp_scenarios())
        with FakeHost() as host:
            p = self.run_suite(host, "--manifest", self.mpath, "--no-mcp", "--checks", "stub", "--skip-remote",
                               "--ids", self.write("a.txt", "T-1 R-1"), "--warm-every", "0")
            self.assertEqual(p.returncode, 0, p.stderr)
            run_dir = p.stdout.strip().splitlines()[-1]
            # (a resume keeps the run's --skip-remote and its absent MCP template)
            p = self.run_suite(host, "--resume", run_dir, "--ids", self.write("b.txt", "C-1"))
            self.assertEqual(p.returncode, 0, p.stderr)
            p = self.run_suite(host, "--resume", run_dir)
            self.assertEqual(p.returncode, 0, p.stderr)
        sessions = sc.read_json(os.path.join(run_dir, "orchestrator.json"))["sessions"]
        self.assertEqual([x["ids"] for x in sessions], ["R-1,T-1", "C-1", "C-1,R-1,T-1"])
        # the third invocation planned nothing new: plan.json holds the two selections' plans
        plan = sc.read_json(os.path.join(run_dir, "plan.json"))
        self.assertEqual([x["selection"]["ids"] for x in plan["plans"]], [["R-1", "T-1"], ["C-1"]])
        with open(os.path.join(run_dir, "report.md")) as fh:
            self.assertNotIn("INCOMPLETE", fh.read().split("# Environment")[0])

    def policy(self, scenarios):
        m = sc.read_json(self.mpath)
        pol = gate_policy.build_policy(m, routes=["POST /resolve", "POST /traverse"], arms=["mcp"], table_path="",
                                       warm_every=10)
        path = os.path.join(self.d, "policy.json")
        sc.write_json_atomic(path, pol)
        return path

    def test_gate_exit_status(self):
        make_manifest(self.mpath, mcp_scenarios())
        pol = self.policy(mcp_scenarios())
        with FakeHost() as host:
            # MCP only, the policy's selection (R-1, T-1, D-1): the gate passes, exit 0
            p = self.run_suite(host, "--manifest", self.mpath, "--mcp-cmd", self.mcp_cmd(), "--skip-local",
                               "--skip-remote", "--gate", pol, "--by-design", "")
            self.assertEqual(p.returncode, 0, p.stderr)
            run_dir = p.stdout.strip().splitlines()[-1]
            self.assertIn("GATE PASS", p.stderr)
            self.assertTrue(sc.read_json(os.path.join(run_dir, "gate.json"))["passed"])
            self.assertEqual(sorted(set(k[0] for k in mcp_keys(run_dir))), ["D-1", "R-1", "T-1"])
            with open(os.path.join(run_dir, "report.md")) as fh:
                self.assertIn("**GATE PASS**", fh.read())
            # the same with a failed read of the service's configuration: the gate fails, exit 1
            p = self.run_suite(host, "--manifest", self.mpath, "--mcp-cmd", self.mcp_cmd(), "--skip-local",
                               "--skip-remote", "--gate", pol, "--by-design", "", FAKE_MCP_CAPS="failed")
            self.assertEqual(p.returncode, 1, p.stderr)
            self.assertIn("configuration evidence missing", p.stderr)
            self.assertIn("fake: connection refused", p.stderr)
            failed_dir = p.stdout.strip().splitlines()[-1]
            # every execution failing: exit 1 with --gate, 0 without (the exploratory exit status is unchanged)
            p = self.run_suite(host, "--manifest", self.mpath, "--mcp-cmd", self.mcp_cmd(), "--skip-local",
                               "--skip-remote", "--gate", pol, "--by-design", "", FAKE_MCP_MODE="fail")
            self.assertEqual(p.returncode, 1, p.stderr)
            # a plain resume of a gated run is gated as it was (it adds nothing here): exit 1 again
            p = self.run_suite(host, "--resume", p.stdout.strip().splitlines()[-1], "--by-design", "")
            self.assertEqual(p.returncode, 1, p.stderr)
            self.assertIn("gate: FAIL", p.stderr)
            p = self.run_suite(host, "--manifest", self.mpath, "--mcp-cmd", self.mcp_cmd(), "--skip-local",
                               "--skip-remote", "--ids", self.write("ids2", "R-1"), "--by-design", "",
                               FAKE_MCP_MODE="fail")
            self.assertEqual(p.returncode, 0, p.stderr)
            # the offline path judges a recorded run again: --reevaluate --gate exits with the gate's status
            p = self.run_suite(host, "--reevaluate", run_dir, "--gate", pol, "--by-design", "")
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertTrue(os.path.exists(os.path.join(run_dir, "gate.reevaluated.json")))
            p = self.run_suite(host, "--reevaluate", failed_dir, "--gate", pol, "--by-design", "")
            self.assertEqual(p.returncode, 1, p.stderr)
            # a policy that is none
            p = self.run_suite(host, "--manifest", self.mpath, "--gate", self.write("x.json", "{}"))
            self.assertEqual(p.returncode, 2)


if __name__ == "__main__":
    unittest.main()
