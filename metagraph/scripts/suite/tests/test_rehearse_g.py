"""Offline tests of the area G rehearsal (rehearse_g.py, lane D): the analysis of a recorded rehearsal made of the
library's fixture graphlets (verdicts, measured facts, the shadow run of the op plans by the local library, the corpus
file and the scenarios the build makes from it), the document facts and the slot rules. No request is sent."""
import json
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


class W3Headroom(unittest.TestCase):
    """W3's headroom rule (DESIGN-v8 §3, refined 2026-10-20): interactive ops a third of the alarm, heavy ops a third
    of the heavy timeout, an escalation allowed and recorded; the values from the capabilities' graphlet block
    (sets compared as sets; a part absent or served null: the fallback, recorded missing; never a crash)."""

    BLOCK = {"interactive": {"timeout_s": 20, "alarm_s": 3.0}, "heavy": {"timeout_s": 120},
             "heavy_ops": ["graphlet_compare", "graphlet_export:fasta|gfa|json", "graphlet_labels:name"]}

    def lanes(self, graphlet):
        import rehearse_g
        live = None if graphlet is ... else {"service_settings": {"graphlet": graphlet}}
        return rehearse_g.graphlet_lanes(sc.References(None, {"service": live}))

    def test_read_from_the_block(self):
        got = self.lanes({"lanes": dict(self.BLOCK, heavy_ops=list(reversed(self.BLOCK["heavy_ops"])))})
        self.assertEqual((got["alarm_s"], got["interactive_timeout_s"], got["heavy_timeout_s"]), (3.0, 20.0, 120.0))
        self.assertEqual(got["heavy_ops"], frozenset(self.BLOCK["heavy_ops"]))   # a set: the served order aside
        self.assertEqual(set(got["sources"].values()), {"live"})
        self.assertEqual(got["missing"], [])
        other = self.lanes({"lanes": dict(self.BLOCK, interactive={"timeout_s": 30, "alarm_s": 4.5},
                                          heavy={"timeout_s": 90})})
        self.assertEqual((other["alarm_s"], other["interactive_timeout_s"], other["heavy_timeout_s"]),
                         (4.5, 30.0, 90.0))

    def test_not_stated_is_the_fallback_recorded_missing(self):
        lanes = dict(self.BLOCK, interactive={"timeout_s": 20, "alarm_s": None})
        got = self.lanes({"lanes": lanes})
        self.assertEqual(got["alarm_s"], 3.0)
        self.assertEqual(got["sources"]["interactive.alarm_s"], "fallback")
        self.assertEqual(got["missing"], ["service_settings.graphlet.lanes.interactive.alarm_s (null; fallback 3.0)"])
        got = self.lanes({"lanes": None})
        self.assertEqual(len(got["missing"]), 4)
        self.assertTrue(all("null at service_settings.graphlet.lanes" in m for m in got["missing"]), got["missing"])
        got = self.lanes(None)
        self.assertTrue(all("null at service_settings.graphlet;" in m for m in got["missing"]), got["missing"])
        got = self.lanes(...)   # no recorded document, the table states no graphlet block: every part absent
        self.assertEqual((got["alarm_s"], got["heavy_timeout_s"]), (3.0, 120.0))
        self.assertTrue(all("(absent;" in m for m in got["missing"]), got["missing"])
        got = self.lanes({"lanes": dict(self.BLOCK, heavy_ops="graphlet_compare", heavy={"timeout_s": True})})
        self.assertEqual(got["heavy_ops"], frozenset(self.BLOCK["heavy_ops"]))
        self.assertEqual(got["heavy_timeout_s"], 120.0)
        self.assertEqual(len(got["missing"]), 2)

    def test_the_lane_of_an_op(self):
        import rehearse_g
        lanes = self.lanes({"lanes": self.BLOCK})
        for tool, args, want in (("graphlet_compare", {}, "heavy"), ("graphlet_export", {"format": "fasta"}, "heavy"),
                                 ("graphlet_export", {"format": "gfa"}, "heavy"),
                                 ("graphlet_export", {"format": "mgt"}, "interactive"),
                                 ("graphlet_labels", {"name": "x"}, "heavy"),
                                 ("graphlet_labels", {"name": None}, "interactive"),
                                 ("graphlet_walks", {"arm": "left"}, "interactive"),
                                 ("graphlet_summary", {}, "interactive")):
            self.assertEqual(rehearse_g.w3_lane({"id": "o", "tool": tool, "args": args}, lanes), want, (tool, args))

    def test_the_rule(self):
        import rehearse_g
        lanes = self.lanes({"lanes": self.BLOCK})
        ops = [{"id": i, "tool": t, "args": {}} for i, t in (
            ("fast", "graphlet_walks"), ("slow", "graphlet_walks"), ("esc", "graphlet_walks"),
            ("esc-slow", "graphlet_walks"), ("heavy", "graphlet_compare"), ("heavy-slow", "graphlet_compare"),
            ("said-no", "graphlet_summary"), ("said-yes", "graphlet_summary"), ("untimed", "graphlet_walk"))]
        w3 = {"fast": {"elapsed_ms": 800.0}, "slow": {"elapsed_ms": 1500.0}, "esc": {"elapsed_ms": 4000.0},
              "esc-slow": {"elapsed_ms": 45000.0}, "heavy": {"elapsed_ms": 30000.0},
              "heavy-slow": {"elapsed_ms": 41000.0},
              "said-no": {"elapsed_ms": 4000.0, "headers": {"X-Graphlet-Op-Escalated": "false"}},
              "said-yes": {"elapsed_ms": 2000.0, "op_meta": {"lane": "heavy", "escalated": True}},
              "untimed": {"error": {"code": "op_busy"}}}
        hr = rehearse_g.w3_headroom(ops, w3, lanes)
        self.assertEqual(hr["timed"], 8)
        # VERIFY-D VD10: no escalation inferred from the time: an unstated op over a third of the alarm fails, at
        # 1.5 s as at 4 s (monotonic); only the stated one is judged as heavy
        self.assertEqual([x.split(":")[0] for x in hr["slow"]], ["slow", "esc", "esc-slow", "heavy-slow", "said-no"])
        self.assertIn("slow: W3 1.50 s, interactive lane's third 1 s", hr["slow"])
        self.assertIn("esc: W3 4.00 s, interactive lane's third 1 s", hr["slow"])
        self.assertEqual(hr["escalations"], ["said-yes: 2.00 s, escalated to heavy (stated)"])
        stated = dict(w3, **{"esc-slow": {"elapsed_ms": 45000.0, "op_meta": {"escalated": True}}})
        hr = rehearse_g.w3_headroom(ops, stated, lanes)
        self.assertIn("esc-slow: W3 45.00 s, heavy (escalated) lane's third 40 s", hr["slow"])

    def test_lanes_of_runs(self):
        import rehearse_g
        tmp = tempfile.mkdtemp(prefix="reh-lanes-")
        try:
            got = rehearse_g.lanes_of_runs([tmp])     # nothing recorded: the table's, every part the fallback
            self.assertIsNone(got["document"])
            self.assertEqual(got["alarm_s"], 3.0)
            os.makedirs(os.path.join(tmp, "arms", "mcp"))
            doc = {"service_settings": {"graphlet": {"lanes": dict(self.BLOCK, heavy={"timeout_s": 60})}}}
            with open(os.path.join(tmp, "arms", "mcp", "run.json"), "w") as fh:
                json.dump({"service_capabilities": {"document": doc}}, fh)
            got = rehearse_g.lanes_of_runs([os.path.join(tmp, "absent"), tmp])
            self.assertEqual(got["heavy_timeout_s"], 60.0)
            self.assertIn("service_capabilities", got["document"])
            self.assertEqual(rehearse_g.lanes_record(got)["headroom_s"], {"interactive": 1.0, "heavy": 20.0})
        finally:
            import shutil
            shutil.rmtree(tmp, ignore_errors=True)


class StopRehearsal(unittest.TestCase):
    """The stop rehearsal (rehearsal 3, VERIFY-D VD6): each stop scenario three times, interleaved, MCP only; the
    constrain / annotate partner (VD8) direct-local; the entry verdict over recorded stop_flow records."""

    def test_its_manifest_and_plan(self):
        # R2-3: production's warm state first (the warmers as in a v8 run, each stop's base request once to
        # completion), then the tries; the pair (R2-5) at the end
        import rehearse_g
        scen, meta = rehearse_g.rehearsal3_scenarios(FakeCtx(), v7_scenarios())
        ids = [x["id"] for x in scen]
        self.assertEqual(ids, ["T-walk-0028", "R-warm-0005", "T-labels-0003", "S-base-c1", "S-base-c2",
                               "S-stop-c1-t1", "S-stop-c2-t1", "S-stop-c1-t2", "S-stop-c2-t2", "S-stop-c1-t3",
                               "S-stop-c2-t3", "G-pair-annot-c1", "G-annot-c1"])
        by = dict((x["id"], x) for x in scen)
        for wid in ("T-walk-0028", "R-warm-0005", "T-labels-0003"):
            self.assertTrue(meta[wid]["warmer"])
            self.assertNotIn("warms", by[wid])
        import areas
        bases = dict(areas.stop_base_requests(v7_scenarios()))
        for n, stop_id in ((1, "S-stop-0001"), (2, "S-stop-0002")):
            b_ = by["S-base-c%d" % n]
            self.assertEqual(meta[b_["id"]], {"base_of": stop_id})
            self.assertIn("not_implemented", b_["mcp"])
            self.assertEqual(b_["request"], bases[b_["route"]])
            self.assertFalse(b_.get("direct_not_run"))
        for x in scen:
            if x["id"].startswith("S-stop"):
                self.assertTrue(x["direct_not_run"])
                self.assertEqual(sc.mcp_repeats(0, x, 0), [0])
                self.assertEqual(x["mcp"]["flow"]["delay_s"], 0.0)
                self.assertEqual(meta[x["id"]]["stop_of"], "S-stop-000%s" % x["id"][8])
        partner, annot = by["G-pair-annot-c1"], by["G-annot-c1"]
        for g in (partner, annot):
            self.assertEqual(g["graphlet"]["ways"]["default"], area_g.DIRECT_ONLY_WAYS)
            self.assertEqual(g["repeat"], 3)
        self.assertFalse([o["id"] for o in partner["graphlet"]["ops"] if o["id"].startswith("cmp")])
        self.assertEqual(meta["G-pair-annot-c1"]["cid"], "annotate.1+constrain")
        self.assertEqual(meta["G-annot-c1"]["cid"], "annotate.1")
        p = rehearse_g.plan3(scen, meta)
        self.assertEqual((p["stop_tries"], p["mcp_jobs"], p["direct_local_requests"]), (6, 12, 17))

    def test_the_entry_verdict(self):
        import rehearse_g
        sys.path.insert(0, HERE)
        from test_stop_flow_checks import record
        scen, meta = rehearse_g.rehearsal3_scenarios(FakeCtx(), v7_scenarios())
        tmp = tempfile.mkdtemp(prefix="reh-s-")
        os.makedirs(os.path.join(tmp, "arms", "mcp", "raw"))
        try:
            for x in scen:
                if not x["id"].startswith("S-stop"):
                    continue
                if x["id"].startswith("S-stop-c2"):
                    from test_stop_flow_checks import ResolveStop, attempt, released
                    rec = record(tool="traverse_resolve", stop_answer=dict(ResolveStop.ANSWER),
                                 attempts=[attempt(released("grace", detail="grace_20s",
                                                            at="2026-10-20T10:00:22.6+00:00"))])
                    if x["id"].endswith("t3"):
                        rec["race"] = "terminal_before_running"
                else:
                    rec = record()
                rec["stop"]["took"] = True
                sc.write_gz_json(os.path.join(tmp, "arms", "mcp", sc.raw_name(x["id"], "mcp", 0)),
                                 {"scenario": x["id"], "answer": {"http": 200, "body": {"stop_flow": rec}}})
            got = rehearse_g.stop_try_verdicts({"scenarios": scen}, meta, tmp)
            self.assertEqual(got["entry"]["S-stop-0001"]["verdict"], "pass", got["entry"])
            self.assertEqual(got["entry"]["S-stop-0002"]["verdict"], "fail: 2 of 3 taken, 2 passed")
            self.assertEqual(got["tries"]["S-stop-c2-t3"]["race"], "terminal_before_running")
            self.assertEqual(got["tries"]["S-stop-c2-t3"]["outcome"], "inconclusive")
        finally:
            import shutil
            shutil.rmtree(tmp, ignore_errors=True)


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
        # VERIFY-D VD9: the plans the build makes from the chosen corpus, shadow-run on the rehearsal's documents
        built = res["built"]
        self.assertGreaterEqual(built["checked"], 1, built)
        self.assertNotIn("G-arms2-0001", built["not_rehearsed"])
        self.assertEqual(built["problems"], {}, built["problems"])
        self.assertIn("G-pair-radius-0001", built["not_rehearsed"])   # no partner was rehearsed here
        with open(os.path.join(self.tmp, "REPORT.md")) as fh:
            self.assertIn("## The built plans", fh.read())
        # the corpus builds: the chosen candidate with its preconditions and timeouts (seeds of v7 again)
        corp = sc.read_json(os.path.join(self.tmp, "corpus.json"))
        g = area_g.gen_G(FakeCtx(), v7_scenarios(), corp, "a" * 40)
        s = dict((x["id"], x) for x in g)["G-arms2-0001"]
        self.assertEqual(s["graphlet"]["candidate"], "two_arms.2")
        self.assertEqual(s["graphlet"]["preconditions"], e["preconditions"])
        # W3's times (the shadow run's): the lanes from the run's recorded capabilities (alarm 6 s: interactive
        # third 2 s); sum at 2.5 s misses it; walks-id-left at 7 s states its escalation: recorded, within the heavy
        # lane's third (VD10: only a stated escalation counts)
        w3dir = os.path.join(self.tmp, "w3")
        os.makedirs(os.path.join(w3dir, "arms", "mcp"))
        doc = {"service_settings": {"graphlet": {"lanes": {
            "interactive": {"timeout_s": 20, "alarm_s": 6.0}, "heavy": {"timeout_s": 120},
            "heavy_ops": ["graphlet_labels:name", "graphlet_compare", "graphlet_export:fasta|gfa|json"]}}}}
        with open(os.path.join(w3dir, "arms", "mcp", "run.json"), "w") as fh:
            json.dump({"service_capabilities": {"document": doc, "read_at": "2026-10-20T00:00:00Z"}}, fh)
        with open(os.path.join(w3dir, "arms", "mcp", "graphlet_ops.jsonl"), "w") as fh:
            for op, ms in (("sum", 2500.0), ("walks-id-left", 7000.0)):
                row = {"scenario": "G-arms2-c2", "op": op, "elapsed_ms": ms}
                if op == "walks-id-left":
                    row["headers"] = {"X-Graphlet-Op-Escalated": "true", "X-Graphlet-Op-Lane": "heavy"}
                fh.write(json.dumps(row) + "\n")
        out2 = os.path.join(self.tmp, "with-w3")
        os.makedirs(out2)
        res2 = rehearse_g.analyse(out2, manifest, meta, dl, mdir, LIB, "a" * 40, mcp_graphlet_run=w3dir)
        r2 = res2["results"]["G-arms2-c2"]
        self.assertEqual(res2["lanes"]["alarm_s"], 6.0)
        self.assertEqual(res2["lanes"]["missing"], [])
        self.assertIn("sum: W3 2.50 s, interactive lane's third 2 s", r2["verdicts"]["op_headroom"])
        self.assertNotIn("walks-id-left: W3", r2["verdicts"]["op_headroom"])
        self.assertEqual(r2["escalations"], ["walks-id-left: 7.00 s, escalated to heavy (stated)"])
        self.assertFalse(r2["passed"])
        rec = read_json(os.path.join(out2, "rehearsal.json"))
        self.assertEqual(rec["w3_lanes"]["headroom_s"], {"interactive": 2.0, "heavy": 40.0})
        self.assertEqual(rec["w3_lanes"]["sources"]["interactive.alarm_s"], "live")
        with open(os.path.join(out2, "REPORT.md")) as fh:
            report = fh.read()
        self.assertIn("## W3 headroom", report)
        self.assertIn("walks-id-left: 7.00 s, escalated to heavy", report)

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
