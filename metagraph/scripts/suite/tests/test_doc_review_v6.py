"""The external review of the v6 results document (2026-10-18; V6-001 .. V6-307): each agreed finding's fix in the
generator, tested where behaviour changes. The review's ids name the tests. Offline; stdlib only."""

import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.dirname(HERE)
DOC = os.path.join(SUITE, "doc")
for p in (SUITE, HERE, DOC):
    if p not in sys.path:
        sys.path.insert(0, p)

import suite_common as sc  # noqa: E402

import blocks  # noqa: E402,F401
import cells  # noqa: E402,F401
import facts  # noqa: E402
import sections  # noqa: E402

PROSE_DIR = os.path.join(DOC, "prose")


class _Stub(object):
    """A bare object with the attributes given (a Run or a Context stand-in)."""

    def __init__(self, **kw):
        for k, v in kw.items():
            setattr(self, k, v)


def prose(name):
    return sc.read_json(os.path.join(PROSE_DIR, "%s.json" % name))


class GatePopulationTest(unittest.TestCase):
    """V6-001: the gate's required population after the policy's automatic exclusions, beside the verdict."""

    def _run(self):
        policy = {
            "arms": ["direct-local", "mcp"],
            "selection": {"routes": ["POST /traverse"], "scenarios": ["T-a", "T-b", "T-c", "T-d"]},
            "exclusions": [
                {"scenario": "T-b", "arms": ["direct-local", "mcp"], "rules": ["budget_stress"], "reason": "load"},
                {"scenario": "T-c", "arms": ["mcp"], "rules": ["manifest_not_implemented", "table_not_expressible"],
                 "reason": "x"},
                {"scenario": "T-d", "arms": ["mcp"], "repeats": [0], "rules": ["marks_cover_repeat"], "reason": "y"},
                {"scenario": "T-a", "arms": ["mcp"], "repeats": [1], "rules": ["declared"], "reason": "owner"}],
            "required": [{"scenario": "T-a", "arm": "direct-local", "repeats": [0]},
                         {"scenario": "T-a", "arm": "mcp", "repeats": [0]},
                         {"scenario": "T-c", "arm": "direct-local", "repeats": [0]},
                         {"scenario": "T-d", "arm": "direct-local", "repeats": [0]},
                         {"scenario": "T-d", "arm": "mcp", "repeats": [1]}]}
        rows = [{"scenario": sid, "arm": arm, "repeat": rep, "outcome": out} for sid, arm, rep, out in (
            ("T-a", "direct-local", 0, "pass"), ("T-a", "mcp", 0, "service_by_design"), ("T-a", "mcp", 1, "pass"),
            ("T-b", "direct-local", 0, "inconclusive"), ("T-c", "direct-local", 0, "pass"),
            ("T-c", "mcp", 0, "not_implemented_mcp"), ("T-d", "direct-local", 0, "pass"),
            ("T-d", "mcp", 0, "inconclusive"), ("T-d", "mcp", 1, "pass"), ("T-e", "direct-local", 0, "fail"))]
        return _Stub(policy=policy, rows=rows, order=["T-a", "T-b", "T-c", "T-d", "T-e"], policy_path="/x/p.json",
                     declared={"exclusions": [{"scenario": "T-a", "arms": ["mcp"], "repeats": [1]}]},
                     gate={"passed": True, "problems": [], "counts": {"judged": 5, "passed": 4}})

    def test_population_and_ledger_text(self):
        run = self._run()
        pop = sections.gate_population(run)
        self.assertEqual(pop["selected"], 4)
        self.assertEqual(pop["required"], [("direct-local", 3, 3), ("mcp", 2, 2)])
        self.assertEqual(len(pop["auto"]), 3, "the owner's declared entry is not an automatic exclusion")
        self.assertEqual(sections.gate_auto_text(run),
                         "3 entries: 1 budget-stress test on direct-local and mcp, 1 request the MCP tools cannot "
                         "express on mcp and 1 repeat whose property checks the marks all set aside on mcp (repeat 0 "
                         "only)")
        self.assertEqual(sections.gate_required_text(run),
                         "3 scenarios (3 executions) on direct-local and 2 scenarios (2 executions) on mcp")
        text = sections.gate_population_text(run)
        self.assertIn("Required means selected by route (the scenarios of POST /traverse (4 of the run's 5)", text)
        self.assertIn("no decision of the owner (3 entries: ", text)
        self.assertIn("not declared exempt by the owner (1 exclusion; I.4, appendix G)", text)
        self.assertIn("5 executions in all. The gate certifies those executions, not the excluded ones", text)
        decision = sections.gate_decision_text(run)
        self.assertTrue(decision.startswith("The gate passed: every planned execution was recorded, every execution "
                                            "the policy requires passed or matched the by-design case"), decision)
        self.assertIn(text, decision)
        self.assertIn("removed by the policy's automatic exclusions", sections.gate_inconclusive_clause(run))

    def test_judged_outcomes(self):
        """The judged executions by outcome, as gate.judge counts them: T-b (budget_stress), T-c via MCP, T-d's MCP
        repeat 0, T-a's declared MCP repeat 1 and T-e (outside the selection) are not judged."""
        run = self._run()
        self.assertEqual(sections.gate_judged_outcomes(run), {"pass": 4, "service_by_design": 1})
        self.assertEqual(sections.judged_outcomes_text(run, {"passed": 4}), " (4 passed, 1 service-side by design)")
        # a count the gate's passed does not confirm falls back to the gate's own number
        self.assertEqual(sections.judged_outcomes_text(run, {"passed": 2}), ", 2 passed")
        self.assertEqual(sections.gate_population(_Stub(policy=None)), None)
        self.assertEqual(sections.gate_auto_text(_Stub(policy=None)), "")

    def test_prose_criterion(self):
        """The blocking criterion reads 'required by the policy on that access path and repeat', never 'in the
        policy's selection and not declared exempt' (which leaves out the automatic exclusions)."""
        texts = [str(prose("overview")), str(prose("glossary"))]
        for t in texts:
            self.assertNotIn("in the policy's selection and not declared exempt", t)
        self.assertIn("the policy requires on that access path and repeat",
                      prose("overview")["sections"]["decision"]["verdict"][0])
        self.assertIn("{gate_auto_text}", prose("overview")["sections"]["provenance"]["exemptions"][0])


class FullyCheckedTest(unittest.TestCase):
    """V6-201: 'fully checked' holds over the property checks no property-level mark sets aside; the chapter's MCP
    sentence names the agreeing pairs that omit a marked check, which stays unverified via MCP."""

    def test_definition(self):
        ctx = _Stub(units=_Stub(by_scenario={}))
        self.assertEqual(facts.mcp_full_definition(ctx, []), "every property check evaluated on both")
        self.assertEqual(facts.mcp_full_definition(ctx, ["T-labels-0003", "T-labels-0005", "T-labels-0006"]),
                         "every property check evaluated on both, except a check a property-level mark of the manifest "
                         "sets aside, not compared and unverified via MCP (on three of them: T-labels-0003, -0005 and "
                         "-0006)")
        many = ["P-exact-%04d" % i for i in range(1, 9)] + ["P-scope-0001", "P-edge-0004"]
        self.assertTrue(facts.mcp_full_definition(ctx, many).endswith(
            "(on ten of them: P-exact 8, P-scope 1, P-edge 1)"))

    def test_prose(self):
        for ch in ("R", "T", "P", "L", "X", "M", "Q", "N", "D"):
            text = str(prose(ch))
            self.assertEqual(text.count("{mcp_full_def}"), 2, ch)
            self.assertNotIn("(every property check evaluated on both)", text, ch)
        depth = prose("overview")["sections"]["reading"]["paragraphs"][3]
        self.assertIn("except a check a property-level mark of the manifest sets aside: that check is not compared "
                      "and stays unverified via MCP", depth)
        gloss = [e for e in prose("glossary")["entries"] if e["term"] == "agreement class, depth"][0]["definition"]
        self.assertIn("but a check a property-level mark of the manifest sets aside", gloss)


if __name__ == "__main__":
    unittest.main()
