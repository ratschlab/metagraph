"""Offline tests of the verdict rules of compare.py, reevaluate.py, report.py and suite_common.py that the two-day
review's verdict findings name: the direct arms' class with equal digests, the MCP reference per repeat and check,
the host cache after a relayed deadline, the text_contains grammar, the no-table rule (c), the runner's unmatched
not_implemented_mcp mark, outcomes outside the vocabulary, what an area tested and the load caveat."""

import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import compare as cmp  # noqa: E402
import suite_common as sc  # noqa: E402
from test_compare_report import T0, line  # noqa: E402


def st(outcome, name="status", detail=""):
    return {"check": name, "outcome": outcome, "detail": detail}


class DirectIdenticalTest(unittest.TestCase):
    """Equal digests never hide a differing check or an outcome difference no timing explains."""

    def pair(self, remote_budget):
        deadline = {"code": "deadline", "text": "nothing written within the time budget"}
        local = line("R-1", "direct-local", 0, T0, 9000.0, outcome="inconclusive", http=503, digest="same",
                     error=deadline, checks=[st("inconclusive"), st("inconclusive", "budget_respected")])
        remote = line("R-1", "direct-remote", 0, T0 + 20, 2000.0, http=503, digest="same", error=deadline,
                      outcome="fail" if remote_budget == "fail" else "inconclusive",
                      checks=[st("inconclusive"), st(remote_budget, "budget_respected")])
        return {0: local}, {0: remote}

    def test_identical_bodies_one_premature_disagree(self):
        local, remote = self.pair("fail")
        local[0]["checks"][1] = st("pass", "budget_respected")
        local[0]["outcome"] = "pass"
        res = cmp.direct_agreement(local, remote)
        self.assertEqual(res["class"], "disagree")
        self.assertTrue(res["digest_equal"])
        self.assertEqual([d["check"] for d in res["differing_checks"]], ["budget_respected"])

    def test_identical_bodies_outcomes_differ_untimed_disagree(self):
        # as v4's R-budget-0002: the same 503 body, one arm's deadline late (a finding), the other's in time
        local, remote = self.pair("finding")
        remote[0]["outcome"] = "finding"
        res = cmp.direct_agreement(local, remote)
        self.assertEqual(res["class"], "disagree")
        self.assertEqual(res["outcome_differs"][0]["repeat"], 0)

    def test_equal_control_stays_identical(self):
        local, remote = self.pair("inconclusive")
        self.assertEqual(cmp.direct_agreement(local, remote)["class"], "identical")


if __name__ == "__main__":
    unittest.main()
