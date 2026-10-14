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


if __name__ == "__main__":
    unittest.main()
