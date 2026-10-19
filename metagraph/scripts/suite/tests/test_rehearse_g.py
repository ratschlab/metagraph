"""Offline tests of the area G rehearsal (rehearse_g.py, lane D): the analysis of a recorded rehearsal made of the
library's fixture graphlets (verdicts, measured facts, the shadow run of the op plans by the local library, the corpus
file and the scenarios the build makes from it), the document facts and the slot rules. No request is sent."""
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.dirname(HERE)
sys.path.insert(0, SUITE)
sys.path.insert(0, HERE)

import area_g  # noqa: E402
import suite_common as sc  # noqa: E402
from test_area_g import FIXTURES, LIB, FakeCtx, needs_library, read_json, v7_scenarios  # noqa: E402


@needs_library
class Rehearsal(unittest.TestCase):
    """rehearse_g.analyse on a recorded rehearsal made of the library's fixture graphlets: verdicts, facts, the
    shadow run, the corpus file."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="reh-g-")

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmp, ignore_errors=True)

    def record(self, scenarios, responses, twins, mcp):
        """A direct-local run dir and an MCP dir: responses[sid] the three direct answers' body, twins[sid] the
        full twin, mcp[sid] the MCP answer's body (None: not run) with the echo of the scenario's strategy."""
        dl = os.path.join(self.tmp, "dl")
        mdir = os.path.join(self.tmp, "mcp")
        for d in (os.path.join(dl, "arms", "direct-local", "raw"), os.path.join(dl, "arms", "direct-local", "docs"),
                  os.path.join(mdir, "raw")):
            os.makedirs(d, exist_ok=True)
        for s in scenarios:
            sid = s["id"]
            for rep in range(3):
                sc.write_gz_json(os.path.join(dl, "arms", "direct-local", sc.raw_name(sid, "direct-local", rep)),
                                 {"scenario": sid, "answer": {"http": 200, "body": responses[sid]}})
            sc.write_gz_json(os.path.join(dl, "arms", "direct-local", "docs", sid + ".full.json.gz"), twins[sid])
            if mcp.get(sid) is not None:
                sc.write_gz_json(os.path.join(mdir, sc.raw_name(sid, "mcp", 0)), {
                    "scenario": sid, "answer": {"http": 200, "body": mcp[sid]},
                    "calls": [{"tool": "traverse_status", "payload": {"overrides": s["request"]["strategy"]}}]})
        return dl, mdir

    def test_analyse_chooses_the_passing_candidate(self):
        import rehearse_g
        cands = dict((c["cid"], c) for c in area_g.candidates(FakeCtx(), v7_scenarios()))
        corpus = area_g.default_corpus()
        scen, meta = [], {}
        for sid, cid in (("G-arms2-c1", "two_arms.1"), ("G-arms2-c2", "two_arms.2")):
            c = cands[cid]
            ops = area_g.op_plan(sid, c["slot"], area_g.facts_of(c, corpus), True, [(sid, False)],
                                 extras={"other": sid})
            scen.append(area_g.g_scenario(sid, sid, c, corpus, None, ops, ["rehearsal"]))
            meta[sid] = {"cid": cid, "slot": c["slot"], "candidate": c}
        switch = read_json(FIXTURES["G-arms2"])
        full = read_json(FIXTURES["G-arms2"].replace(".graphlet.", ".full."))
        trace = read_json(FIXTURES["G-trace"])
        for r in (switch, trace):
            r["timing"] = {"elapsed_ms": 120.0}   # (the fixtures keep no timing; a staging answer states it)
        dl, mdir = self.record(scen, {"G-arms2-c1": switch, "G-arms2-c2": switch},
                               {"G-arms2-c1": full, "G-arms2-c2": full},
                               {"G-arms2-c1": trace, "G-arms2-c2": switch})   # c1: the MCP body differs
        manifest = {"scenarios": scen}
        res = rehearse_g.analyse(self.tmp, manifest, meta, dl, mdir, LIB, "a" * 40)
        r1, r2 = res["results"]["G-arms2-c1"], res["results"]["G-arms2-c2"]
        self.assertTrue(r1["verdicts"]["bodies"].startswith("fail"), r1["verdicts"])
        self.assertTrue(r2["passed"], r2["verdicts"])
        for k in ("echo", "bodies", "envelope", "headroom", "twin", "slot", "plan"):
            self.assertEqual(r2["verdicts"][k], "pass", (k, r2["verdicts"]))
        corpus = res["corpus"]
        self.assertEqual(corpus["slots"]["two_arms"]["candidate"], "two_arms.2")
        self.assertIn("two_arms.1: bodies fail", corpus["slots"]["two_arms"]["why"])
        e = corpus["candidates"]["two_arms.2"]
        self.assertEqual(e["preconditions"]["arms"], ["left", "right"])
        self.assertEqual(e["preconditions"]["walks_min"], {"left": 12, "right": 1})
        self.assertEqual(e["strategy"], scen[1]["request"]["strategy"])
        self.assertTrue(e["timeouts"] and all(v >= rehearse_g.TIMEOUT_FLOOR_S for v in e["timeouts"].values()))
        self.assertIn("summary_cut_bytes", e["facts"])
        for name in ("rehearsal.json", "corpus.json", "REPORT.md"):
            self.assertTrue(os.path.exists(os.path.join(self.tmp, name)), name)
        # the corpus builds: the chosen candidate with its preconditions and timeouts (seeds of v7 again)
        corp = sc.read_json(os.path.join(self.tmp, "corpus.json"))
        g = area_g.gen_G(FakeCtx(), v7_scenarios(), corp, "a" * 40)
        s = dict((x["id"], x) for x in g)["G-arms2-0001"]
        self.assertEqual(s["graphlet"]["candidate"], "two_arms.2")
        self.assertEqual(s["graphlet"]["preconditions"], e["preconditions"])

    def test_document_facts_and_slot_rules(self):
        import rehearse_g
        f = rehearse_g.document_facts(read_json(FIXTURES["G-caps"]))
        self.assertEqual(f["arms"], ["right"])
        self.assertFalse(f["bases"])
        self.assertEqual(rehearse_g.slot_problem("no_bases", f, 31), "")
        self.assertTrue(rehearse_g.slot_problem("linear", f, 31))
        fail = rehearse_g.document_facts({"results": [{"error": {"code": "x"}}]})
        self.assertTrue(fail["failed"])
        self.assertEqual(rehearse_g.slot_problem("failed_seed", fail, 31), "")
        self.assertEqual(rehearse_g.preconditions_of(fail), {"no_graphlet": True})


if __name__ == "__main__":
    unittest.main()
