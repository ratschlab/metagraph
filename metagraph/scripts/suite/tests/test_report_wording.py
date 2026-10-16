"""The report's honest wording (lane S6c; GPT review, VERIFY-A, VERIFY-B): the front cost's host-warm pairs meet a
recent-repeat heuristic, not a measured cache state; the resource figures are contextual telemetry; an MCP `agree`
is split by how much of the scenario's claim was compared (fully checked, partial, on failure, generic checks only,
refusal), and a reference failure the MCP cannot see is named."""

import contextlib
import io
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
from test_compare_report import CASES, build_run  # noqa: E402

GENERIC = ("status", "invariants", "counts_consistent", "budget_respected")


class ReportWordingTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.d = self.tmp.name
        build_run(self.d)
        merge_run.merge(self.d)
        sc.write_json_atomic(os.path.join(self.d, "compare.json"), cmp.compare(self.d, (CASES, {"cases": 5})))
        with contextlib.redirect_stdout(io.StringIO()):
            self.text = report.build_report(self.d)

    def tearDown(self):
        self.tmp.cleanup()

    def test_cache_and_telemetry_wording(self):
        self.assertIn("both executions meet the recent-repeat heuristic", self.text)
        self.assertIn("cache residency and foreign load were not controlled", self.text)
        self.assertNotIn("warm row cache), so the difference is not a cache effect", self.text)
        self.assertIn("# Resource telemetry per area (contextual)", self.text)
        self.assertIn("Contextual telemetry, not a measurement of the suite's requests", self.text)
        self.assertIn("contextual telemetry: sra\\_all CPU mean", self.text)

    def test_agreement_depth(self):
        comp = sc.read_json(os.path.join(self.d, "compare.json"))
        rows, _ = sc.read_jsonl(os.path.join(self.d, "executions.jsonl"))
        depth = report.agreement_depth(comp, rows, lambda n: n not in GENERIC)
        # P-001 compares its property on both repeats; R-001 agrees on its status only: the MCP does not expose
        # contains_source_context, which the reference (direct-local) fails
        self.assertEqual(depth["P-001"]["class"], "full")
        self.assertEqual(depth["R-001"], {"class": "generic_only", "not_compared": ["contains_source_context"],
                                          "partly_compared": {}, "reference_fails_hidden": ["contains_source_context"]})
        verdict = self.text.split("# Verdict", 1)[1].split("\n# ", 1)[0]
        self.assertIn("3 agree (2 fully checked, 0 partial, 0 partial with a non-pass aggregate, 0 on failure, "
                      "1 generic checks only, 0 refusals; see Agreement)", verdict)
        self.assertIn("generic checks only: the scenario has property checks, none compared | 1 | R-001", self.text)
        self.assertIn("Agreeing although the reference arm fails a check the MCP does not expose (the MCP cannot "
                      "show it holds): R-001: contains\\_source\\_context.", self.text)

    def test_depth_table_as_suite_v6_prints_it(self):
        # the merge of gr/suite-doc into gr/suite-v6: agreement_depth splits the property checks some MCP repeat
        # did not compare into not_compared (no repeat compared them) and partly_compared (another repeat did), for
        # the results document; the report's depth table lists both, sorted, as gr/suite-v6's report.py printed its
        # not_compared, with gr/suite-v6's header
        comp = sc.read_json(os.path.join(self.d, "compare.json"))
        split = {"P-001": {"class": "partial", "not_compared": ["zeta_check"],
                           "partly_compared": {"alpha_check": {"compared": [1], "skipped": [0]}},
                           "reference_fails_hidden": []}}
        joined = {"P-001": {"class": "partial", "not_compared": ["alpha_check", "zeta_check"],
                            "reference_fails_hidden": []}}
        self.assertEqual(report.skipped_on_any_repeat(split["P-001"]), ["alpha_check", "zeta_check"])
        text = report.agreement_section(comp, split)
        self.assertEqual(text, report.agreement_section(comp, joined))
        self.assertIn("| P-001 (alpha\\_check, zeta\\_check) |", text)
        self.assertIn("| agree, by depth | n | scenarios (property checks not compared) |", text)
        self.assertNotIn(" on r0", text)
        # agreement_depth itself keeps the doc's split
        prop = lambda o: {"check": "walk_in_record", "outcome": o}  # noqa: E731
        st = {"check": "status", "outcome": "pass"}
        c = {"scenarios": {"X": {"mcp_agreement": {"class": "agree", "reference": "direct-local",
                                                   "mcp_outcome": "pass", "reference_outcome": "pass"}}}}
        rows = [{"scenario": "X", "arm": "mcp", "repeat": 0, "checks": [st, prop("not_exposed")]},
                {"scenario": "X", "arm": "mcp", "repeat": 1, "checks": [st, prop("pass")]},
                {"scenario": "X", "arm": "direct-local", "repeat": 0, "checks": [st, prop("pass")]}]
        d = report.agreement_depth(c, rows, lambda n: n not in GENERIC)["X"]
        self.assertEqual((d["class"], d["not_compared"], d["partly_compared"]),
                         ("partial", [], {"walk_in_record": {"compared": [1], "skipped": [0]}}))
        self.assertEqual(report.skipped_on_any_repeat(d), ["walk_in_record"])

    def test_depth_classes(self):
        def run(mcp_checks, ref_checks, mcp_out="pass", ref_out="pass", expect=None):
            comp = {"scenarios": {"X": {"mcp_agreement": {"class": "agree", "reference": "direct-local",
                                                          "mcp_outcome": mcp_out, "reference_outcome": ref_out}}}}
            rows = [{"scenario": "X", "arm": "mcp", "repeat": r, "checks": c} for r, c in enumerate(mcp_checks)]
            rows += [{"scenario": "X", "arm": "direct-local", "repeat": 0, "checks": ref_checks}]
            return report.agreement_depth(comp, rows, lambda n: n not in GENERIC)["X"]["class"]

        st = {"check": "status", "outcome": "pass"}
        prop = lambda o: {"check": "walk_in_record", "outcome": o}  # noqa: E731
        bud = {"check": "budget_respected", "outcome": "not_exposed"}
        self.assertEqual(run([[st, prop("pass"), bud]] * 2, [st, prop("pass"), dict(bud, outcome="pass")]), "full")
        self.assertEqual(run([[st, prop("pass")], [st, prop("not_exposed")]], [st, prop("pass")]), "partial")
        self.assertEqual(run([[st, prop("pass")], [st, prop("inconclusive")]], [st, prop("pass")],
                             mcp_out="inconclusive"), "partial_non_pass")
        self.assertEqual(run([[st, prop("fail")]], [st, prop("fail")], "fail", "fail"), "on_failure")
        self.assertEqual(run([[st, prop("not_exposed")]], [st, prop("pass")]), "generic_only")
        self.assertEqual(run([[st]], [st]), "refusal")


if __name__ == "__main__":
    unittest.main()
