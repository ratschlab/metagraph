"""An outcome outside the shared vocabulary (suite_common.outcome_known) is never a success on the runner side (lane
S7r, P6-end-to-end-02's runner side): merge_run reads such a line as error (outcome_runner keeps the recorded one),
and the paired loop counts it like an error towards its stop."""

import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.dirname(HERE)
sys.path.insert(0, SUITE)
sys.path.insert(0, HERE)

import compare as cmp  # noqa: E402
import merge_run  # noqa: E402
import report  # noqa: E402
import suite_common as sc  # noqa: E402
from fake_host import FakeHost  # noqa: E402
from test_compare_report import T0, line, write_arm  # noqa: E402
from test_orchestrator import SuiteCase, mcp_scenarios  # noqa: E402
from test_runner import make_manifest  # noqa: E402

MANIFEST = {"scenarios": [{"id": "P-%d" % i, "route": "POST /pattern", "mcp": {"tool": "pattern_search"},
                           "expect": [{"check": "status", "args": {"http": 200}}]} for i in range(6)]}
STATUS_PASS = [{"check": "status", "outcome": "pass", "detail": "http 200"}]


class MergeUnknownOutcomeTest(unittest.TestCase):
    def test_unknown_outcome_reads_as_error(self):
        # the review's failing input: six MCP lines with outcome skipped_by_runner beside passing direct ones read
        # PASS; and a line with no outcome at all
        with tempfile.TemporaryDirectory() as d:
            os.makedirs(os.path.join(d, "arms", "direct-local"))
            sc.write_json_atomic(os.path.join(d, "arms", "direct-local", "manifest.json"), MANIFEST)
            write_arm(d, "direct-local", [line("P-%d" % i, "direct-local", 0, T0 + i, 5.0, checks=STATUS_PASS)
                                          for i in range(6)])
            mcp = [line("P-%d" % i, "mcp", 0, T0 + 10 + i, 5.0, outcome="skipped_by_runner", checks=STATUS_PASS)
                   for i in range(5)]
            mcp.append(dict(line("P-5", "mcp", 0, T0 + 20, 5.0, checks=STATUS_PASS), outcome=None))
            write_arm(d, "mcp", mcp)
            merged = merge_run.merge(d)
            rows, _ = sc.read_jsonl(os.path.join(d, "executions.jsonl"))
            mrows = [r for r in rows if r["arm"] == "mcp"]
            self.assertEqual(set(r["outcome"] for r in mrows), {"error"})
            self.assertEqual(sorted(str(r["outcome_runner"]) for r in mrows), ["None"] + ["skipped_by_runner"] * 5)
            self.assertEqual(merged["outcomes_unknown"], {"mcp": {"None": 1, "skipped_by_runner": 5}})
            self.assertIn("arms/mcp: 6 line(s) with an outcome outside the shared vocabulary (None 1, "
                          "skipped_by_runner 5) read as error (outcome_runner keeps the recorded one)",
                          merged["warnings"])
            comp = cmp.compare(d, ({}, None))
            sc.write_json_atomic(os.path.join(d, "compare.json"), comp)
            text = report.build_report(d)
        self.assertIn("**FAIL (", text)
        # the comparison and the report state the recorded outcome as outside the vocabulary (one mechanism: the
        # merge's error, read back through outcome_runner), counted once
        totals = comp["summary"]["totals"]["mcp"]["outcomes"]
        self.assertEqual((totals["other"], totals["error"]), (6, 0))
        outside = comp["summary"]["outcomes_outside_vocabulary"]
        self.assertIn({"scenario": "P-0", "arm": "mcp", "repeat": 0, "outcome": "skipped_by_runner",
                       "read_as": "error"}, outside)
        self.assertIn({"scenario": "P-5", "arm": "mcp", "repeat": 0, "outcome": None, "read_as": "error"}, outside)
        self.assertEqual(comp["scenarios"]["P-0"]["mcp_agreement"]["by"], "outcome")
        self.assertIn("6 with an outcome outside the vocabulary", text)
        self.assertIn("'skipped_by_runner' (the merge read it as error)", text.replace("\\_", "_"))
        self.assertIn("outcome 'skipped_by_runner' outside the vocabulary (read as error)", text.replace("\\_", "_"))
        # every known outcome stays as recorded
        rows = [{"arm": "mcp", "outcome": o} for o in sc.OUTCOMES]
        self.assertEqual(merge_run.unknown_outcomes(rows), {})
        self.assertEqual([r["outcome"] for r in rows], list(sc.OUTCOMES))


class PairedLoopUnknownOutcomeTest(SuiteCase):
    def test_unknown_mcp_outcomes_count_as_errors(self):
        # three MCP executions in a row with an outcome outside the vocabulary stop the loop like three errors
        make_manifest(self.mpath, mcp_scenarios())
        with FakeHost() as host:
            p = self.run_suite(host, "--manifest", self.mpath, "--mcp-cmd", self.mcp_cmd(), "--checks", "stub",
                               FAKE_MCP_MODE="unknown")
            self.assertEqual(p.returncode, 1, p.stderr)
        self.assertIn("recorded an outcome outside the shared vocabulary", p.stderr)
        run_dir = p.stdout.strip().splitlines()[-1]
        loop = sc.read_json(os.path.join(run_dir, "orchestrator.json"))["sessions"][-1]["loop"]
        self.assertEqual((loop["state"], loop["reason"]), ("failed", "mcp_errors"))
        rows, _ = sc.read_jsonl(os.path.join(run_dir, "executions.jsonl"))
        self.assertEqual(sorted((r["scenario"], r["repeat"], r["outcome"], r["outcome_runner"]) for r in rows
                                if r["arm"] == "mcp"),
                         [("C-1", 0, "error", "skipped_by_runner"), ("C-1", 1, "error", "skipped_by_runner"),
                          ("P-1", 0, "error", "skipped_by_runner")])


if __name__ == "__main__":
    unittest.main()
