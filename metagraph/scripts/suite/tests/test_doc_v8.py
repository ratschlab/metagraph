"""The v8 results document as a delta on v7 (2026-10-20): manifest v8 is v7's 725 scenarios unchanged, ten warm-ups of
T-walk-0024, area G (the graphlet tools compared five ways) and the stop scenarios; the unit mapping is v7's with
chapter R's warm-up unit extended and chapters G and S added, the prose's v8 statements are `manifest_version >= 8`
branches (a v5 / v6 / v7 build stays as it was). Offline; stdlib only."""

import json
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

import units as units_mod  # noqa: E402

PROSE_DIR = os.path.join(DOC, "prose")


def prose(name):
    return sc.read_json(os.path.join(PROSE_DIR, "%s.json" % name))


class _Stub(object):
    def __init__(self, **kw):
        for k, v in kw.items():
            setattr(self, k, v)


def v8_manifest_stub():
    """A manifest shaped like v8 for the mapping: v7's scenarios as they are plus, for every scenario the v8 mapping
    adds, a stub of its area (the real v8 manifest is frozen by the owner, so the test does not depend on it)."""
    v7 = sc.read_json(os.path.join(SUITE, "manifest", "suite_v7.json"))
    raw8 = sc.read_json(os.path.join(DOC, "units", "suite_v8.json"))
    ids7 = [s["id"] for s in v7["scenarios"]]
    scen = dict((s["id"], s) for s in v7["scenarios"])
    order = [k for k in raw8 if not k.startswith("_")]
    out = []
    for sid in order:
        out.append(scen.get(sid) or {"id": sid, "area": raw8[sid]["chapter"], "title": sid})
    return {"manifest_version": 8, "scenarios": out, "relations": v7.get("relations")}, ids7, raw8


class V8MappingTest(unittest.TestCase):
    def test_mapping_is_v7_plus_the_new_scenarios(self):
        m8, ids7, raw8 = v8_manifest_stub()
        raw7 = sc.read_json(os.path.join(DOC, "units", "suite_v7.json"))
        self.assertEqual(raw8["_manifest"]["version"], 8)
        self.assertEqual(units_mod.default_mapping_path({"manifest_version": 8, "content_sha256": "not frozen"}),
                         os.path.join(DOC, "units", "suite_v8.json"))
        # v7's entries, units and families as they are; R-warm gains the ten warm-ups of T-walk-0024
        for sid in ids7:
            self.assertEqual(raw8[sid], raw7[sid], sid)
        for uid, u in raw7["_units"].items():
            if uid == "R-warm":
                continue
            self.assertEqual(raw8["_units"][uid], u, uid)
        added = ["R-warm-%04d" % i for i in range(9, 19)]
        self.assertEqual(raw8["_units"]["R-warm"]["members"], raw7["_units"]["R-warm"]["members"] + added)
        self.assertEqual(raw8["_units"]["R-warm"]["n"], 18)
        for sid in added:
            self.assertRegex(raw8[sid]["variant"], r"\Aof T-walk-0024 \(\d+ bp\)\Z")
        for fid, f in raw7["_families"].items():
            self.assertEqual(raw8["_families"][fid], f, fid)
        # chapters G and S after T, the others in v7's order
        self.assertEqual(list(raw8["_chapters"]), ["R", "T", "G", "S"] + list(raw7["_chapters"])[2:])
        self.assertEqual(sorted(f for f, v in raw8["_families"].items() if v["chapter"] in ("G", "S")),
                         ["G-corpus", "G-pair", "G-warm", "S-base", "S-stop"])
        units = units_mod.resolve(m8, units_mod.load_mapping(os.path.join(DOC, "units", "suite_v8.json")))
        self.assertEqual(units_mod.validate(units, m8), [])
        self.assertEqual(len(units.chapter_scenarios("G")), 19)
        self.assertEqual(len(units.chapter_scenarios("S")), 4)


class ChapterGTest(unittest.TestCase):
    """Chapter G's comparison sections (doc/graphlets.py) over a five-way comparison of lane C's fixture run: the ways
    and what each agreement proves, the comparison per tool, the disagreements grouped by kind with both values, the
    pre-checks, the service checks, the correctness coverage and what is not covered (computed from the manifest)."""

    def comparison(self, change_w2=False):
        import test_report_graphlets as trg
        case = trg.Section("test_an_agreeing_run")
        case.setUp()
        self.addCleanup(case.tearDown)
        g = case.comparison(change_w2=change_w2, by_design=True)
        return g, case.run_

    def test_sections(self):
        import graphlets
        g, run_ = self.comparison(change_w2=True)
        ctx = _Stub(run=_Stub(manifest=run_.manifest, dir=run_.dir))
        text = graphlets.ways_section(ctx)
        self.assertIn("\\section*{The five ways}", text.replace("\\subsection*", "\\section*"))
        self.assertIn("the suite's own reader, without the library", text)
        self.assertIn("W5 is the one correctness check", text)
        tools = graphlets.per_tool_section(g)
        self.assertIn("\\texttt{graphlet\\_\\allowbreak{}walks}", tools)
        self.assertIn("W1 2, W2 2, W3 2, W4 2, W5 1", tools)
        dis = graphlets.disagreements_section(g)
        self.assertIn("1 op call disagreed, in 1 kind of difference", dis)
        self.assertIn("W1-W2 (transport) at \\texttt{/\\allowbreak{}answer/\\allowbreak{}total}", dis)
        pre = graphlets.prechecks_section(g)
        self.assertIn("G-\\allowbreak{}merge-\\allowbreak{}0001", pre)
        service = graphlets.service_checks_section(g)
        self.assertIn("walk\\_\\allowbreak{}display", service)
        cov = graphlets.coverage_section(g)
        self.assertIn("consistency only", cov)
        nc = graphlets.not_covered(ctx.run)
        self.assertIn("v9", nc)
        word, _ = graphlets.verdict(g)
        self.assertEqual(word, "FAIL")
        g2, _ = self.comparison()
        self.assertIn("No op call disagreed.", graphlets.disagreements_section(g2))

    def test_no_graphlets_no_chapter(self):
        import graphlets
        run = _Stub(manifest={"manifest_version": 7, "scenarios": [{"id": "R-1", "area": "R"}]}, dir="/nowhere")
        self.assertIsNone(graphlets.data(run))
        self.assertEqual(graphlets.chapter_g_sections(_Stub(run=run)), "")


class ChapterSTest(unittest.TestCase):
    """Chapter S's records: per stop execution the stop's answer, the job's end and each host attempt's release
    record; the base requests' times; nothing without a stop scenario."""

    def run_stub(self, rec):
        scen = {"S-stop-0001": {"id": "S-stop-0001", "area": "S", "mcp": {"tool": "traverse_fetch",
                                                                           "flow": {"kind": "stop"}}},
                "S-base-0001": {"id": "S-base-0001", "area": "S"}}
        lines = {("S-stop-0001", "mcp"): [{"scenario": "S-stop-0001", "repeat": 0, "outcome": "pass",
                                           "stop_flow": rec}],
                 ("S-base-0001", "direct-local"): [{"scenario": "S-base-0001", "repeat": 0, "outcome": "pass",
                                                    "wall_ms": 15002.0, "server": {"elapsed_ms": 15000.0}},
                                                   {"scenario": "S-base-0001", "repeat": 1, "outcome": "pass",
                                                    "wall_ms": 9000.0, "server": {"elapsed_ms": 8990.0}}]}
        return _Stub(order=list(scen), scenario=lambda sid: scen.get(sid) or {},
                     arm_lines=lambda sid, arm: lines.get((sid, arm), []))

    def test_records(self):
        import graphlets
        rec = {"stop": {"took": True, "running": 1, "cancel_supported": True, "bound_seconds": 55},
               "terminal_status": "stopped", "release_reads": 44,
               "attempts": [{"attempt_id": "452019c1-e0c3-4460-932b-29afd4cf5a0f:1", "server_instance": "0c43fc04",
                             "release": {"reason": "completed", "site": "worker", "decision": None,
                                         "detail": "response"}}]}
        text = graphlets.chapter_s_sections(_Stub(run=self.run_stub(rec)))
        self.assertIn("taken; running 1, cancel supported, bound 55 s", text)
        self.assertIn("completed, site worker, decision null, response", text)
        self.assertIn("stopped", text)
        self.assertIn("repeat 1 is the warm full-run time", text)
        self.assertIn("& r1 & pass & 9.0 s & 9.0 s", text)
        text = graphlets.chapter_s_sections(_Stub(run=self.run_stub(None)))
        self.assertIn("no stop record", text)
        empty = _Stub(order=[], scenario=lambda sid: {}, arm_lines=lambda sid, arm: [])
        self.assertEqual(graphlets.chapter_s_sections(_Stub(run=empty)), "")


class GateV8Test(unittest.TestCase):
    """The v8 gate on the overview pages: area G's population per policy, both policies side by side, area G's op calls
    beside the executions, its problems grouped (a disagreement or a failed pre-check or service check a violation,
    an inconclusive call not shown), the graphlet library's identity in plain words."""

    def test_population_text(self):
        import sections
        pol = {"counts": {"graphlets": {"scenarios": 17, "op_calls": 1495, "w3_not_exposed": 650, "service_checks": 8,
                                        "ways": {"W1": 1493, "W2": 1493, "W3": 845, "W4": 1118, "W5": 1492}}}}
        verdict = {"counts": {"graphlet_ops": {"judged": 1495, "agree": 1490, "service_by_design": 5,
                                               "not_required": 0}}}
        self.assertEqual(sections.graphlet_population_text(pol, verdict),
                         "17 scenarios, 1495 op calls required (W1 1493, W2 1493, W3 845, W4 1118, W5 1492 records), "
                         "650 calls W3 does not expose, 8 service checks; judged 1495: 1490 agree, 5 by design")
        self.assertEqual(sections.graphlet_population_text({"counts": {}}, verdict), "")

    def test_problem_groups(self):
        import sections
        for detail, group in (("op x (graphlet_walks) disagrees (W1-W2): ...", "violation"),
                              ("op x (graphlet_walks) inconclusive: the selector failed", "unproven"),
                              ("pre-check body fail (a determinism finding): ...", "violation"),
                              ("pre-check twin inconclusive: no twin kept", "unproven"),
                              ("service check max_bytes_refused fail: 1 of 1", "violation"),
                              ("service check in_service required, not passed (n/a)", "unproven"),
                              ("the graphlet comparison could not be computed (X: y)", "evidence")):
            self.assertEqual(sections.gate_group_of({"kind": "graphlet", "detail": detail}), group, detail)

    def test_library_note(self):
        import sections
        rec = {"recorded": True, "library_commits": {"fdff625ef7322e919b28da4b81763acec65591af": 4}, "invocations": 4,
               "library_commit": "fdff625ef7322e919b28da4b81763acec65591af", "source": "arms/mcp/mcp_runs.jsonl"}
        pin = "fdff625ef7322e919b28da4b81763acec65591af"
        gate = {"library_identity": {"expected": None, "recorded": rec, "pin": pin},
                "policy": {"arms": ["mcp"], "graphlets": {}}}
        text, _ = sections.library_note_text(_Stub(run=_Stub(gate=gate, identity=None)))
        self.assertIn("on all 4 MCP invocations, the commit the manifest pins; the gate was given no record of the "
                      "deployed library", text)
        gate["library_identity"]["expected"] = pin
        ident = {"deployed_library_commit": {"commit": pin, "source": "the deployment record"}}
        text, _ = sections.library_note_text(_Stub(run=_Stub(gate=gate, identity=ident)))
        self.assertIn("and equals the deployed one (the deployment record)", text)
        self.assertTrue(sections.library_held(_Stub(gate=gate)))
        gate["library_identity"]["expected"] = "0" * 40
        self.assertFalse(sections.library_held(_Stub(gate=gate)))
        self.assertEqual(sections.library_note_text(_Stub(run=_Stub(gate={}, identity=None))), ("", []))


class AreaLineAndGlossaryTest(unittest.TestCase):
    """Area G's verdict line on the cover and in I.1 (none without a comparison), and the glossary's v8 terms, printed
    from manifest v8 on only (a v5 to v7 glossary is as it was)."""

    def test_area_line(self):
        import graphlets
        g = {"summary": {"scenarios": {"pass": 2, "fail": 1}, "ops": {"agree": 40, "disagree": 2}},
             "service_checks": [{"check": "a", "outcome": "pass"}, {"check": "b", "outcome": "n/a"}],
             "scenarios": {}}
        run = _Stub(manifest={"manifest_version": 8}, _graphlets=g)
        word, _ = graphlets.verdict(g)
        cover = graphlets.area_line(_Stub(run=run))
        self.assertIn("\\textbf{\\large %s}" % word, cover)
        self.assertIn("(3 graphlet scenarios: 2 pass, 1 fail; 42 op calls: 40 agree, 2 disagree; service checks: "
                      "1 pass, 1 n/a; chapter G)", cover)
        self.assertIn("the deploy gate requires them (I.1)", cover)
        self.assertNotIn("(I.1)", graphlets.area_line(_Stub(run=run), cover=False))
        self.assertIn("ERROR", graphlets.area_line(_Stub(run=_Stub(manifest={}, _graphlets={"error": "x"}))))
        self.assertIsNone(graphlets.area_line(_Stub(run=_Stub(manifest={}, _graphlets=None))))

    def test_glossary_since(self):
        import sections
        entries = prose("glossary")["entries"]
        new = [e["term"] for e in entries if e.get("since")]
        self.assertEqual(set(e.get("since") for e in entries if e.get("since")), {8})
        self.assertIn("way (W1 to W5)", new)
        self.assertIn("stop, release record", new)
        for v in (5, 6, 7, None):
            self.assertEqual([e["term"] for e in sections.glossary_for(entries, v)],
                             [e["term"] for e in entries if not e.get("since")])
        self.assertEqual(len(sections.glossary_for(entries, 8)), len(entries))
        terms = set(e["term"] for e in entries)
        for e in entries:
            if e.get("since"):
                for s in e.get("see") or []:
                    self.assertIn(s, terms)


class AppendixBAndSinceV8Test(unittest.TestCase):
    """Appendix B's area G part: area G's own by-design table, case by case, with the op calls that name a case and
    the comparison's evidence (none without a comparison)."""

    def test_by_design_g_section(self):
        import graphlets
        g = {"summary": {}, "scenarios": {"G-x-0001": {"ops": {"exp-mgt": {"by_design": [
            {"case": "graphlet.export_mgt_is_the_stored_body", "verified": True}]}}},
            "G-y-0001": {"ops": {"exp-mgt": {"by_design": [
                {"case": "graphlet.export_mgt_is_the_stored_body", "verified": False, "side": "W3"}]}}}}}
        manifest = {"manifest_version": 8, "scenarios": [
            {"id": "G-x-0001", "graphlet": {"ops": [{"id": "exp-mgt", "by_design":
                                                    "graphlet.export_mgt_is_the_stored_body"}]}}]}
        run = _Stub(manifest=manifest, dir="/nowhere", _graphlets=g, policy={"graphlets": {"by_design_g_table": {
            "sha256": "0" * 64}}})
        text = graphlets.by_design_g_section(_Stub(run=run))
        self.assertIn("Area G's by-design table", text)
        self.assertIn("2 op calls: 1 verified on both sides", text)
        self.assertIn("1 not verified", text)
        self.assertIn("the deploy gate's policy pins", text)
        self.assertIn("not exercised in this run", text)       # the table's other cases
        self.assertEqual(text.count("\\texttt{graphlet."), 4)
        self.assertEqual(graphlets.by_design_g_section(_Stub(run=_Stub(manifest={}, _graphlets=None))), "")


class SharedProseV8Test(unittest.TestCase):
    """The v8 branches of the shared prose: I.2's access paths name the graphlet tools and traverse_stop, its "Not
    covered" paragraph is v7's for v7 only and v8's says where the graphlet tools and the stops are tested; T-walk's
    family names the warm-ups of T-walk-0024; v5 to v7 read as before."""

    def test_branches(self):
        import facts
        pars = prose("overview")["sections"]["methods"]["paragraphs"]
        vals = {"mcp_base": "u", "local_base": "l", "remote_base": "r", "mcp_tools_text": "four tools",
                "mcp_tools_called": "TOOLS", "traversal_tools_not_called": "NOT"}
        access8 = facts.fill(pars[0], dict(vals, manifest_version=8), unresolved=[])
        self.assertIn("`catalog`, the graphlet tools (chapter G) and `traverse_stop` (chapter S))", access8)
        self.assertIn("`catalog`) on", facts.fill(pars[0], dict(vals, manifest_version=7), unresolved=[]))
        nc = [p for p in pars if "*Not covered.*" in p]
        self.assertEqual(len(nc), 2)
        for v in (6, 7):
            self.assertEqual(facts.fill(nc[1], dict(vals, manifest_version=v), unresolved=[]).strip(), "")
        self.assertIn("Not called, and so not covered: NOT;", facts.fill(nc[0], dict(vals, manifest_version=7),
                                                                         unresolved=[]))
        self.assertEqual(facts.fill(nc[0], dict(vals, manifest_version=8), unresolved=[]).strip(), "")
        v8 = facts.fill(nc[1], dict(vals, manifest_version=8), unresolved=[])
        self.assertIn("The graphlet tools are compared five ways in chapter G", v8)
        self.assertIn("so not covered here: NOT.", v8)
        walk = prose("T")["families"]["T-walk"]["prose"]
        self.assertIn("T-walk-0024 follows ten warm-ups of its region (R-warm-0009 to -0018, chapter R)",
                      facts.fill(walk, {"manifest_version": 8}, unresolved=[]))
        self.assertNotIn("ten warm-ups", facts.fill(walk, {"manifest_version": 7}, unresolved=[]))


if __name__ == "__main__":
    unittest.main()
