"""Two cases the paired three-arm smoke of suite_smoke_v3 showed: a /traverse attempt that reached its duration
bound (the MCP arm's seed error) is a stated stop, so its property checks are inconclusive, not failures; and a
by-design case may name a field with a placeholder (bounds.<knob>)."""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import checks  # noqa: E402
import compare  # noqa: E402

BOUND = ("the attempt reached the duration bound the server enforces for it (15000 ms: the seeds' time budgets "
         "plus the server's allowance, at most hard_cap_ms, the content timeout less one second; see "
         "usage.bound) while its response was built or written: nothing of it is delivered")


class AttemptBound(unittest.TestCase):
    def scenario(self):
        return {"route": "POST /traverse", "request": {"seeds": [{"sequence": "TTAGGG" * 10}]},
                "expect": [{"check": "status", "args": {"http": 200}},
                           {"check": "invariants", "args": {}},
                           {"check": "traverse_seed_outcome", "args": {"walks": ["complete", "partial"]}}]}

    def test_bound_is_inconclusive(self):
        ans = {"http": 200, "body": {"results": [{"error": BOUND}]}, "error": None}
        res = checks.evaluate(self.scenario(), ans, "mcp")
        by = {r["check"]: r["outcome"] for r in res}
        self.assertEqual("pass", by["status"])
        self.assertEqual("inconclusive", by["traverse_seed_outcome"])
        self.assertEqual("inconclusive", by["invariants"])
        self.assertEqual("inconclusive", checks.scenario_outcome(res))

    def test_other_error_still_fails(self):
        ans = {"http": 200, "body": {"results": [{"error": "Seed is not fully present in the graph"}]},
               "error": None}
        res = checks.evaluate(self.scenario(), ans, "mcp")
        self.assertEqual("fail", checks.scenario_outcome(res))


class FieldPlaceholder(unittest.TestCase):
    def test_placeholder(self):
        self.assertTrue(compare._field_matches("bounds.time_budget_ms", "bounds.<knob>"))
        self.assertFalse(compare._field_matches("strategy.time_budget_ms", "bounds.<knob>"))
        self.assertFalse(compare._field_matches("bounds.a.b", "bounds.<knob>"))
        self.assertTrue(compare._field_matches("max_mismatches", "max_mismatches"))
        self.assertFalse(compare._field_matches("max_contexts", "max_mismatches"))


if __name__ == "__main__":
    unittest.main()
