"""Chapter G of the report (report.graphlet_section; DESIGN-v8 §1, §5): per tool its calls, ways and outcomes, every
disagreement with both answers shown compactly, the pre-check table, the service checks, the correctness-coverage
table (W5's checked fields against the consistency-only ones), the by-design pairs and the drift lane; a headline
beside the verdict; and no chapter, no new file, for a run whose manifest has no graphlet scenario."""

import copy
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import compare_graphlets as cg  # noqa: E402
import graphlet_fixtures as gf  # noqa: E402
import report  # noqa: E402
import test_gate as tg  # noqa: E402
import test_gate_graphlets as tgg  # noqa: E402
import test_graphlet_service_checks as tsc  # noqa: E402

SID = "G-merge-0001"
UUID = "0d5b9b1a-1111-2222-3333-444455556666"


class Section(tsc.ServiceBase):
    WALKS = {"id": "walks-len", "tool": "graphlet_walks", "ways": ["W1", "W2", "W3", "W4", "W5"],
             "args": {"arm": "left", "rank": "length", "n": 10, "max_bytes": 2048}}

    def comparison(self, change_w2=False, by_design=False):
        ops = [self.WALKS]
        table = None
        if by_design:
            ops.append(dict(tsc.ByDesign.OP))
            table = [tsc.CAP_CASE]
        self.make([gf.scenario(ops=ops)], table=table)
        lib = gf.walks_answer()
        w2 = dict(lib, total=3) if change_w2 else lib
        self.put(SID, self.WALKS, {"*": lib, "W2": w2, "W3": dict(lib, handle=UUID)},
                 sf={"walk_display": "first_parent"}, ways=["W1", "W2", "W3", "W4"])
        self.run_.record("W5", SID, "walks-len", "graphlet_walks", self.WALKS["args"],
                         answer={"rows": [{"length_bp": 120}]},
                         oracle={"checked": {"/rows/0/length_bp": "SPEC-core §7.1"},
                                 "consistency_only": ["/rows/0/n_alive"], "inputs": ["full", "mgt"]})
        if by_design:
            self.put(SID, tsc.ByDesign.OP, {"*": gf.walks_answer()}, errors={"W3": tsc.ByDesign.REFUSED})
        return self.result()

    def test_an_agreeing_run(self):
        g = self.comparison(by_design=True)
        self.assertEqual(report.graphlet_verdict(g)[0], "PASS")
        text = report.graphlet_section(g)
        self.assertIn("# G. Graphlet tools, compared five ways", text)
        self.assertIn("**Area G (graphlet tools, five ways): PASS**", text)
        self.assertIn("## G.1 Per tool", text)
        self.assertIn("| graphlet\\_walks | 2 | W1 2, W2 2, W3 2, W4 2, W5 1 | 1 | 0 | 1 |", text)
        self.assertIn("## G.2 Disagreements\n\nNone.", text)
        self.assertIn("| G-merge-0001 | pass | pass | pass | pass | pass | pass | pass |", text)
        self.assertIn("| walk\\_display | pass | 1 | yes |", text)
        self.assertIn("| graphlet\\_walks | /rows/\\*/length\\_bp | checked | SPEC-core S7.1 | 1 |", text)
        self.assertIn("| graphlet\\_walks | /rows/\\*/n\\_alive | consistency only | - | 1 |", text)
        self.assertIn("## G.6 Declared per-way differences", text)
        self.assertIn("graphlet.max\\_bytes\\_above\\_service\\_cap | yes |", text)

    def test_a_disagreement_shows_both_answers(self):
        g = self.comparison(change_w2=True)
        word, reason = report.graphlet_verdict(g)
        self.assertEqual(word, "FAIL")
        self.assertIn("1 fail", reason)
        text = report.graphlet_section(g)
        self.assertIn("- G-merge-0001 walks-len (graphlet\\_walks): W1-W2", text)
        self.assertIn("W1-W2 (transport) at /answer/total: W1 2, W2 3", text)
        self.assertIn('        - W1: \\{"answer":\\{"arm":"left"', text)   # markdown-escaped braces
        self.assertIn('        - W2: \\{"answer":', text)

    def test_prechecks_inconclusive_and_drift(self):
        g = self.comparison()
        os.remove(os.path.join(self.run_.dir, "arms", "direct-local", "docs", SID + ".full.json.gz"))
        self.run_.record("drift", SID, "walks-len", "graphlet_walks", self.WALKS["args"],
                         answer=dict(gf.walks_answer(), total=9))
        g = self.result()
        self.assertEqual(report.graphlet_verdict(g)[0], "INCONCLUSIVE")
        text = report.graphlet_section(g)
        self.assertIn("| G-merge-0001 | pass | pass | inconclusive |", text)
        self.assertIn("- G-merge-0001 twin inconclusive: no twin kept", text)
        self.assertIn("## G.7 Drift lane (report-only)", text)
        self.assertIn("at /answer/total: W4 2, drift 9", text)

    def test_an_uncomputable_comparison(self):
        self.assertEqual(report.graphlet_verdict({"error": "x"})[0], "ERROR")
        self.assertIn("could not be computed", report.graphlet_section({"error": "boom"}))


class BuildReport(tgg.GateG):
    def test_the_report_of_a_v8_run_has_chapter_g(self):
        self.build(change={("walks-len", "W2"): lambda kw: kw["answer"].update(total=3)})
        text = report.build_report(self.d)
        self.assertIn("**Area G (graphlet tools, five ways): FAIL**", text)
        self.assertLess(text.index("**Area G"), text.index("# G. Graphlet tools"))
        self.assertLess(text.index("# Relations between scenarios"), text.index("# G. Graphlet tools"))
        self.assertLess(text.index("# G. Graphlet tools"), text.index("# Timing"))
        self.assertFalse(os.path.exists(os.path.join(self.d, cg.COMPARE_FILE)))   # recomputed in memory (VC3)

    def test_a_run_without_graphlet_scenarios_has_none(self):
        m = tg.manifest_of()
        d = tg.make_run(os.path.join(self.tmp, "v7"), m, tg.mcp_lines(m))
        text = report.build_report(d)
        self.assertNotIn("Area G", text)
        self.assertNotIn("G. Graphlet tools", text)
        self.assertFalse(os.path.exists(os.path.join(d, cg.COMPARE_FILE)))


if __name__ == "__main__":
    unittest.main()
