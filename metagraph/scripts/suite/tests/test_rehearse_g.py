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
        for k in ("echo", "bodies", "envelope", "headroom", "twin", "slot", "plan", "fit"):
            self.assertEqual(r2["verdicts"][k], "pass", (k, r2["verdicts"]))
        corpus = res["corpus"]
        self.assertEqual(corpus["slots"]["two_arms"]["candidate"], "two_arms.2")
        self.assertIn("two_arms.1: bodies fail", corpus["slots"]["two_arms"]["why"])
        e = corpus["candidates"]["two_arms.2"]
        self.assertEqual(e["preconditions"]["arms"], ["left", "right"])
        self.assertEqual(e["preconditions"]["walks_min"], {"left": 12, "right": 1})
        self.assertEqual(e["strategy"], scen[1]["request"]["strategy"])
        # without W3's own times the ops keep their lanes; the local medians are recorded
        self.assertEqual(e["timeouts"], {})
        self.assertTrue(res["results"]["G-arms2-c2"]["op_medians_ms"])
        t, slow = rehearse_g.op_timeouts([{"sum": {"elapsed_ms": 4000.0}}, {"sum": {"elapsed_ms": 9000.0}}],
                                         [{"id": "sum", "tool": "graphlet_summary", "args": {}}])
        self.assertEqual(t, {"sum": 20.0})       # min(lane 20 s, max(5 s, 10 x 6.5 s))
        self.assertEqual(len(slow), 0)
        # a failed seed's facts, the cut lists the engine states
        self.assertFalse(rehearse_g.document_facts(switch, full).get("failed"))
        self.assertIn("summary_cut_bytes", e["facts"])
        for name in ("rehearsal.json", "corpus.json", "REPORT.md"):
            self.assertTrue(os.path.exists(os.path.join(self.tmp, name)), name)
        # the corpus builds: the chosen candidate with its preconditions and timeouts (seeds of v7 again)
        corp = sc.read_json(os.path.join(self.tmp, "corpus.json"))
        g = area_g.gen_G(FakeCtx(), v7_scenarios(), corp, "a" * 40)
        s = dict((x["id"], x) for x in g)["G-arms2-0001"]
        self.assertEqual(s["graphlet"]["candidate"], "two_arms.2")
        self.assertEqual(s["graphlet"]["preconditions"], e["preconditions"])

    def test_a_second_rehearsal_merges_into_the_first(self):
        # rehearsal 1 (G-arms2-c1 failed its twin), then rehearsal 2's twin re-check of it, direct-local only: the
        # merged result keeps the first's MCP verdicts and takes the second's twin
        import json as js
        import rehearse_g
        prior = os.path.join(self.tmp, "r1")
        os.makedirs(prior)
        cands = dict((c["cid"], c) for c in area_g.candidates(FakeCtx(), v7_scenarios()))
        corpus = area_g.default_corpus()
        scen = []
        for sid, cid in (("G-arms2-c1", "two_arms.1"), ("G-arms2-c2", "two_arms.2")):
            c = cands[cid]
            ops = area_g.op_plan(sid, c["slot"], area_g.facts_of(c, corpus), True, [(sid, False)])
            scen.append(area_g.g_scenario(sid, sid, c, corpus, None, ops, ["rehearsal"]))
        ok = dict((k, "pass") for k in ("echo", "bodies", "envelope", "headroom", "slot", "plan", "fit"))
        facts = {"arms": ["left", "right"], "walks": {"left": 3, "right": 4}, "segments": {"left": 5, "right": 22},
                 "merges": {"left": 0, "right": 5}, "walk_bp": {"left": 1000, "right": 1000}, "bases": True,
                 "label_names": ["272634"], "body_bytes": 9000, "cut_lists": False, "max_labels_per_node": 70}
        res1 = {"G-arms2-c1": {"verdicts": dict(ok, twin="fail: parents"), "passed": False, "facts": facts,
                               "times_ms": {"mcp": 2367.0}, "bodies": {"mcp": "a" * 64}, "echo": {"x": 1}},
                "G-arms2-c2": {"verdicts": dict(ok, twin="pass"), "passed": True, "facts": facts,
                               "times_ms": {"mcp": 892.0}, "bodies": {"mcp": "b" * 64}, "echo": {"x": 1}}}
        with open(os.path.join(prior, "rehearsal.json"), "w") as fh:
            js.dump({"results": res1, "meta": {"G-arms2-c1": {"cid": "two_arms.1", "slot": "two_arms"},
                                               "G-arms2-c2": {"cid": "two_arms.2", "slot": "two_arms"}}}, fh)
        with open(os.path.join(prior, "rehearsal_manifest.json"), "w") as fh:
            js.dump({"scenarios": scen}, fh)
        again = dict(scen[0], mcp={"not_implemented": rehearse_g.REHEARSAL2_DIRECT})
        res2 = {"G-arms2-c1": {"verdicts": {"twin": "pass", "plan": "pass", "fit": "pass", "slot": "pass",
                                            "bodies": "pass", "mcp": "not expressible (direct-only)"},
                               "passed": True, "facts": facts, "times_ms": {"dl0": 10.0}}}
        meta2 = {"G-arms2-c1": {"cid": "two_arms.1", "slot": "two_arms", "candidate": cands["two_arms.1"]}}
        results, meta, sc_ = rehearse_g.merge_prior(prior, res2, meta2, {"G-arms2-c1": again})
        r = results["G-arms2-c1"]
        self.assertEqual(r["verdicts"]["twin"], "pass")
        self.assertEqual(r["verdicts"]["echo"], "pass")          # the first rehearsal's MCP verdicts kept
        self.assertEqual(r["times_ms"], {"mcp": 2367.0})
        self.assertTrue(r["passed"])
        self.assertIn("rehearsal 2", r["twin_rechecked"])
        self.assertEqual(sc_["G-arms2-c1"]["mcp"].get("tool"), "traverse_fetch")   # the scenario as the build makes it
        corpus = rehearse_g.choose_corpus(results, meta, sc_, None)
        self.assertEqual(corpus["slots"]["two_arms"]["candidate"], "two_arms.1")   # first in order, now passing
        self.assertEqual(corpus["candidates"]["two_arms.1"]["preconditions"]["merges_min"], {"left": 0, "right": 5})

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
        # cut label lists as the engine states them (rehearsal 2026-10-19: every segment's entry list whole while 10
        # recorded lists were cut): a label_lists limitation, label_evidence lower_bound
        resp = read_json(FIXTURES["G-caps"])
        resp["results"][0]["limitations"] = [{"kind": "label_lists", "knob": "labels.max_labels_per_node"}]
        self.assertTrue(rehearse_g.document_facts(resp)["cut_lists"])
        resp["results"][0]["limitations"] = []
        resp["results"][0]["outcome"] = {"label_evidence": "lower_bound"}
        self.assertTrue(rehearse_g.document_facts(resp)["cut_lists"])

    def test_a_shared_seed_is_judged_by_its_cold_time(self):
        import rehearse_g
        s = {"request": {"seeds": [{"sequence": "ACGT" * 10}], "strategy": {"bounds": {"time_budget_ms": 30000}}}}
        results = {"first": {"times_ms": {"dl0": 300.0, "mcp": 14676.0}, "verdicts": {"headroom": "fail"}},
                   "later": {"times_ms": {"dl0": 170.0, "mcp": 430.0}, "verdicts": {"headroom": "pass"}}}
        seen = {}
        rehearse_g.shared_seed_headroom("first", s, {}, results, seen)
        rehearse_g.shared_seed_headroom("later", s, {}, results, seen)
        self.assertTrue(results["later"]["verdicts"]["headroom"].startswith("fail: the seed's cold engine time 14676"))
        results["pair"] = {"times_ms": {"mcp": 40.0}, "verdicts": {"headroom": "pass"}}
        rehearse_g.shared_seed_headroom("pair", s, {"pair": "radius"}, results, seen)
        self.assertEqual(results["pair"]["verdicts"]["headroom"], "pass")   # a partner runs after its graphlet anyway


if __name__ == "__main__":
    unittest.main()


def have_cache():
    import truth
    return os.path.isdir(truth.CACHE_DIR) and bool(os.listdir(truth.CACHE_DIR))


@unittest.skipUnless(have_cache(), "no record cache")
class SecondRehearsal(unittest.TestCase):
    def test_its_manifest(self):
        import rehearse_g
        ctx, v7 = rehearse_g.build_context()
        scen, meta = rehearse_g.rehearsal2_scenarios(ctx, v7)
        ids = [x["id"] for x in scen]
        self.assertEqual(ids[:9], ["G-arms2-c1", "G-pair-twoseeds-c1", "G-merge-c1", "G-pair-labelcap-c1",
                                   "G-cut-c1", "G-cut-c2", "G-cut-c3", "G-nobases-c1", "G-nobases-c2"])
        self.assertEqual(ids[9:], ["G-warm-c4-1", "G-warm-c4-2", "G-long-c4", "G-warm-c5-1", "G-warm-c5-2",
                                   "G-long-c5"])
        for x in scen[:9]:
            self.assertEqual((x["repeat"], x["mcp"]["not_implemented"]), (1, rehearse_g.REHEARSAL2_DIRECT))
            self.assertNotIn("max_memory_mb", x["graphlet"]["twin"]["strategy"]["bounds"])
        for x in (scen[9], scen[10], scen[12], scen[13]):
            self.assertEqual((x["timing_class"], x["route"]), ("long", "POST /resolve"))
        for x in (scen[11], scen[14]):
            self.assertEqual((x["repeat"], x["mcp"]["tool"]), (3, "traverse_fetch"))
            self.assertEqual(len(x["request"]["seeds"][0]["sequence"]), 6600)
        problems, _ = area_g.validate_g([x for x in scen if "graphlet" in x], set())
        self.assertEqual(problems, [])
