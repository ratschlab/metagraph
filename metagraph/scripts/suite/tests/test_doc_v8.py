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
        self.assertEqual(len(units.chapter_scenarios("G")), sum(1 for k, v in raw8.items()
                                                               if not k.startswith("_") and v["chapter"] == "G"))
        self.assertEqual(len(units.chapter_scenarios("S")), 4)
        # the frozen manifest v8 (the owner's freeze, 2026-10-21): the mapping is keyed to it and maps every scenario
        frozen = sc.read_json(os.path.join(SUITE, "manifest", "suite_v8.json"))
        self.assertEqual(raw8["_manifest"]["content_sha256"], frozen["content_sha256"])
        self.assertEqual(units_mod.default_mapping_path(frozen), os.path.join(DOC, "units", "suite_v8.json"))
        self.assertEqual(units_mod.validate(units_mod.resolve(frozen, units_mod.load_mapping(
            os.path.join(DOC, "units", "suite_v8.json"))), frozen), [])


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
                                                    "wall_ms": 15002.0, "server": {"elapsed_ms": 15000.0},
                                                    "time_limited": "seed 0 stopped by its time budget"},
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
        self.assertIn("repeat 1 is the warm time a stop must beat", text)
        self.assertIn("& r1 & pass & complete & 9.0 s & 9.0 s", text)
        # a base request its own time budget stopped (the shadow run's S-base-0001): the wall is a lower bound
        self.assertIn("& r0 & pass & seed 0 stopped by its time budget & 15 s & 15 s", text)
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
    the comparison's evidence (none without a comparison); I.4's added scenarios by area (manifest v8 on)."""

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
        import compare as cmp
        # one row per case of the table the comparison reads (here the tree's: its case count is an input, it
        # changes when a case lands or leaves, e.g. round 7's arm-absent cases)
        self.assertEqual(text.count("\\texttt{graphlet."), len(cmp.load_by_design(cmp.DEFAULT_BY_DESIGN_G)[0]))
        self.assertEqual(graphlets.by_design_g_section(_Stub(run=_Stub(manifest={}, _graphlets=None))), "")

    def test_added_scenarios(self):
        import sections
        ids = (["R-warm-%04d" % i for i in range(9, 19)] + ["G-f%d-0001" % i for i in range(19)]
               + ["S-stop-0001", "S-stop-0002", "S-base-0001", "S-base-0002"])
        run = _Stub(order=["T-walk-0001"] + ids, scenario=lambda sid: {"area": sid.split("-")[0]})
        self.assertEqual(sections.added_scenarios_text(run, set(ids)),
                         "R-warm-0009 to -0018 (chapter R); 19 scenarios (chapter G); S-stop-0001, -0002, "
                         "S-base-0001, -0002 (chapter S)")


class ChapterGReadingV8Test(unittest.TestCase):
    """Chapter G read as a whole: the store-management smoke is shown but kept out of the per-tool total (it is
    outside the five-way comparison, as the summary's op counts are), quoted messages read with typographic quotes,
    and I.1's map asks chapters G's and S's questions."""

    def test_smoke_outside_the_total(self):
        import graphlets
        g = {"scenarios": {}, "summary": {"per_tool": {
            "graphlet_walks": {"calls": 10, "agree": 9, "disagree": 1, "smoke": 0, "not_exposed": 0},
            "graphlet_list": {"calls": 1, "agree": 1, "smoke": 1, "not_exposed": 1}}}}
        text = graphlets.per_tool_section(g)
        self.assertIn("graphlet\\_\\allowbreak{}list} {\\scriptsize (smoke)}", text)
        self.assertIn("\\textbf{all tools} {\\scriptsize (without the smoke)}", text)
        self.assertIn("\\textbf{10}", text)
        self.assertIn("graphlet\\_list: a smoke outside the five-way comparison", text)

    def test_quoted(self):
        import graphlets
        self.assertEqual(graphlets._quoted('failed: "args.walk: no answer"'), "failed: \u201cargs.walk: no answer\u201d")
        self.assertIsNone(graphlets._quoted(None))

    def test_chapter_map(self):
        cm = dict((e["chapter"], e["question"]) for e in prose("overview")["sections"]["decision"]["chapter_map"])
        self.assertIn("graphlet tools answer alike", cm["G"])
        self.assertIn("traverse_stop", cm["S"])

    def test_mcp_only_pairs_counted(self):
        import sections
        entries = {"S-stop-0001": {"mcp_agreement": {"class": "no_direct_reference", "mcp_outcome": "pass"}},
                   "S-base-0001": {"mcp_agreement": {"class": "not_implemented_mcp"}}}
        run = _Stub(comp_entry=lambda sid: entries[sid], depth={})
        _, mtext = sections._agreement_counts(_Stub(run=run), ["S-stop-0001", "S-base-0001"])
        self.assertIn("1 run via MCP only (no direct answer to compare with)", mtext)
        _, mtext = sections._agreement_counts(_Stub(run=run), ["S-base-0001"])
        self.assertNotIn("via MCP only", mtext)

    def test_r_warm_v8(self):
        import facts
        r = prose("R")
        fam, unit = r["families"]["R-warm"]["prose"], r["units"]["R-warm"]["property"]
        v7 = facts.fill(fam, {"manifest_version": 7}, unresolved=[])
        v8 = facts.fill(fam, {"manifest_version": 8}, unresolved=[])
        self.assertTrue(v7.startswith("Eight resolves"))
        self.assertTrue(v8.startswith("Eighteen resolves"))
        self.assertNotIn("R-warm-0009", v7)
        self.assertIn("R-warm-0009 to -0018, ten windows of 682 to 688 bp", v8)
        self.assertIn("a taxid column for the thirteen others", v8)
        self.assertIn("T-labels-0003 and N-hard-0022)", facts.fill(unit, {"manifest_version": 7}, unresolved=[]))
        self.assertIn("N-hard-0022 and T-walk-0024)", facts.fill(unit, {"manifest_version": 8}, unresolved=[]))


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


class RealRunV8Test(unittest.TestCase):
    """What the shadow run's real records showed that the synthetic run did not (the dry build of 2026-10-21): an
    added scenario's executions are 'added', not 'changed outcome'; one execution warm after an interrupted one reads
    in the singular; a declared file given by path keeps its name and path for evidence.json; evidence.json and I.4
    name area G's own by-design table and chapter G's inputs (none before v8); the units builder keys the mapping to
    the manifest it is given."""

    def test_since_previous_added_not_changed(self):
        import sections

        def row(sid, arm, rep, outcome):
            return {"scenario": sid, "arm": arm, "repeat": rep, "outcome": outcome}
        prev_s = [{"id": "P-a-0001", "route": "POST /pattern", "area": "P"}]
        cur_s = prev_s + [{"id": "G-x-0001", "route": "POST /traverse", "area": "G"}]
        prev = _Stub(manifest={"manifest_version": 7, "scenarios": prev_s}, scenarios={"P-a-0001": prev_s[0]},
                     rows=[row("P-a-0001", "mcp", 0, "pass")], reevaluated=False, run={"run_id": "R7"})
        prev.run_id = lambda: "R7"
        prev.lines = lambda: dict(((r["scenario"], r["arm"], r["repeat"]), r) for r in prev.rows)
        cur = _Stub(manifest={"manifest_version": 8, "scenarios": cur_s}, order=[s["id"] for s in cur_s],
                    scenarios=dict((s["id"], s) for s in cur_s),
                    rows=[row("P-a-0001", "mcp", 0, "pass"), row("P-a-0001", "mcp", 1, "pass"),
                          row("G-x-0001", "direct-local", 0, "pass"), row("G-x-0001", "mcp", 0, "pass")],
                    previous=prev, service_identity={}, mcp_run={})
        cur.scenario = cur.scenarios.get
        text = sections.since_previous(_Stub(run=cur))
        self.assertIn("1 added, 0 removed", text)
        self.assertIn("added: G-x-0001 (chapter G)", text)
        self.assertIn("1 execution changed outcome", text)
        self.assertIn("no execution to pass (P-a-0001 mcp r1)", text)
        self.assertNotIn("G-x-0001 direct-local r0", text)

    def test_one_execution_warm_after_interrupted(self):
        import facts
        tpl = [t for t in prose("closing")["observations"] if t["id"] == "load"][0]["template"]
        base = {"load_sra_mean": "0", "load_load1_mean": "0.8", "load_samples": "6783", "load_caveat": "does not apply"}
        one = facts.fill(tpl, dict(base, warm_after_interrupted_n="1", warm_after_interrupted_one=True))
        self.assertIn("1 execution ran warm after an interrupted one and is in no bucket.", one)
        four = facts.fill(tpl, dict(base, warm_after_interrupted_n="4", warm_after_interrupted_one=False))
        self.assertIn("4 executions ran warm after an interrupted one and are in no bucket.", four)

    def test_declared_by_path(self):
        import tempfile
        import reader
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "declared_x.json")
            with open(path, "w") as fh:
                json.dump({"declared_version": 1, "exclusions": []}, fh)
            data = reader.Run._load_declared(_Stub(), path)
            self.assertEqual(data["_file"], "declared_x.json")
            self.assertEqual(data["_path"], os.path.abspath(path))

    def test_graphlet_evidence(self):
        import sections
        run = _Stub(manifest={"manifest_version": 7, "scenarios": []}, dir="/nowhere")
        run._graphlets = None
        self.assertEqual(sections.graphlet_evidence(run), {})
        self.assertIsNone(sections.by_design_g_table_text(run))
        run = _Stub(manifest={"manifest_version": 8, "scenarios": []}, dir="/nowhere", policy={})
        run._graphlets = {"summary": {}}
        run._graphlet_digests = {"manifest.json": "a" * 64, "arms/dl-lib/graphlet_ops.jsonl": "b" * 64,
                                 "docs": "c" * 64, "code": "d" * 64, "by_design_g_table": "e" * 64}
        ev = sections.graphlet_evidence(run)
        self.assertEqual(sorted(ev["graphlet_inputs"]), ["by_design_g_table", "code", "docs"])
        self.assertIn("by_design_table_g", ev)

    def test_by_design_g_table_text(self):
        import graphlets
        import sections
        info = {"version": 1, "cases": 5, "service_commit": "f" * 40, "sha256": "1" * 64, "path": "/x"}
        orig = graphlets.by_design_g_table
        try:
            graphlets.by_design_g_table = lambda run: ({}, info)
            run = _Stub(dir="/nowhere", policy={"graphlets": {"by_design_g_table": {"sha256": "1" * 64}}})
            run._graphlets = {"summary": {}}
            text = sections.by_design_g_table_text(run)
            self.assertIn("version 1, 5 cases", text)
            self.assertIn("the deploy gate's policy pins it", text)
            run.policy = {"graphlets": {"by_design_g_table": {"sha256": "2" * 64}}}
            self.assertIn("\\textbf{not} the table the deploy gate's policy pins", sections.by_design_g_table_text(run))
        finally:
            graphlets.by_design_g_table = orig

    def test_units_builder_keys_the_given_manifest(self):
        import subprocess
        import tempfile
        m8, ids7, raw8 = v8_manifest_stub()
        with tempfile.TemporaryDirectory() as d:
            for s in m8["scenarios"]:
                if s["id"] not in ids7:
                    s.setdefault("route", "POST /traverse")
            m8["content_sha256"] = "f" * 64
            mpath, out = os.path.join(d, "m8.json"), os.path.join(d, "units.json")
            with open(mpath, "w") as fh:
                json.dump(m8, fh)
            proc = subprocess.run([sys.executable, "-I", "-B", os.path.join(DOC, "units_v8_build.py"), mpath, out],
                                  stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            self.assertEqual(proc.returncode, 0, proc.stdout)
            built = sc.read_json(out)
            self.assertEqual(built["_manifest"], {"version": 8, "content_sha256": "f" * 64,
                                                  "scenarios": len(m8["scenarios"])})
            self.assertEqual(built["_units"]["S-base"]["unit_title"], raw8["_units"]["S-base"]["unit_title"])
            self.assertEqual(built["_families"], raw8["_families"])
            for uid, u in raw8["_units"].items():
                self.assertEqual(built["_units"][uid]["members"], u["members"], uid)

    def test_g_pair_prose_names_the_rows_pair_when_present(self):
        import facts
        text = prose("G")["families"]["G-pair"]["prose"]
        without = facts.fill(text, {"xr": {"G-pair-rows-0002": facts._Blank()}})
        self.assertIn("one seed under two label caps, and constrain against annotate;", without)
        self.assertNotIn("compare rows", without)
        withp = facts.fill(text, {"xr": {"G-pair-rows-0002": {"id": "G-pair-rows-0002"}}})
        self.assertIn("constrain against annotate, and a column walk at a small label cap against a trie of one seed, "
                      "whose compare rows page on the service's tools;", withp)

    def test_by_design_pairs_by_case(self):
        import graphlets
        ops = dict(("cmp-%02d" % i, {"by_design": [{"case": "graphlet.compare_export_line", "verified": True}]})
                   for i in range(45))
        ops["exp-mgt"] = {"by_design": [{"case": "graphlet.export_mgt_is_the_stored_body", "verified": False,
                                         "side": "W3", "detail": "W3 answered otherwise"}]}
        g = {"scenarios": {"G-x-0001": {"ops": ops}}}
        text = graphlets.by_design_and_drift_section(g)
        self.assertEqual(text.count("\\texttt{graphlet."), 3)        # two case rows, the unverified pair's row
        self.assertIn("& 45 & 45 &", text)
        self.assertIn("cmp-44", text.replace("\\_", "_").replace("\\allowbreak{}", ""))   # no cap on the names
        self.assertIn("the by-design pairs not verified", text)
        self.assertIn("W3 answered otherwise", text)

    def test_run_directory_outside_the_tree(self):
        import facts
        text = prose("overview")["appendix_F"]["paragraphs"][0]
        inside = facts.fill(text, {"run_id": "R1", "run_dir_in_tree": True})
        self.assertIn("The run directory, metagraph/scripts/suite/runs/R1, is not in the repository (runs/ is "
                      "git-ignored): it stays in the working tree that ran the suite, so", inside)
        outside = facts.fill(text, {"run_id": "R1", "run_dir_in_tree": False})
        self.assertIn("The run directory, R1, is not in the repository: it was written outside the working tree that "
                      "ran the suite (evidence.json names its folder), so", outside)
        self.assertIn("run_dir_in_tree", prose("overview")["placeholders"])

    def test_tools_not_called_with_every_export_format(self):
        import sections
        base = {"traverse_resolve": {"n": 1, "formats": {}}, "traverse_fetch": {"n": 1, "formats": {}}}
        v7 = _Stub(mcp_tool_calls=lambda: dict(base, graphlet_export={"n": 3, "formats": {"mgt": 3}}))
        self.assertIn("the graphlet processing tools (graphlet_summary", sections.traversal_tools_not_called_text(v7))
        self.assertIn("graphlet_export's other formats", sections.traversal_tools_not_called_text(v7))
        g = dict((t, {"n": 1, "formats": {}}) for t in sections.TRAVERSAL_TOOLS
                 if t.startswith("graphlet_") or t in ("traverse_status", "traverse_results", "traverse_stop"))
        g["graphlet_export"] = {"n": 4, "formats": {"fasta": 1, "gfa": 1, "json": 1, "mgt": 1}}
        v8 = _Stub(mcp_tool_calls=lambda: dict(base, **g))
        self.assertEqual(sections.traversal_tools_not_called_text(v8), "traverse_lock and traverse_capabilities")
        g["graphlet_export"] = {"n": 2, "formats": {"mgt": 1, "json": 1}}
        self.assertEqual(sections.traversal_tools_not_called_text(_Stub(mcp_tool_calls=lambda: dict(base, **g))),
                         "the graphlet processing tools (graphlet_export's other formats), traverse_lock and "
                         "traverse_capabilities")

    def test_disagreements_ragged(self):
        import graphlets
        g = {"scenarios": {"G-x-0001": {"ops": {"walks-arm-only": {"outcome": "disagree", "tool": "graphlet_walks",
                                                                   "disagreeing": ["W3-W4"], "pairs": []}}}}}
        text = graphlets.disagreements_section(g)
        self.assertIn("1 op call disagreed", text)
        self.assertIn("{\\small\\raggedright\\begin{itemize}", text)

    def test_agreement_not_planned_not_implemented(self):
        import sections
        comp = {"G-a-0001": {"mcp_agreement": {"class": "agree"}, "direct_agreement": {"class": "identical"}},
                "G-d-0001": {"mcp_agreement": {"class": "not_implemented_mcp"}},
                "G-w-0001": {"mcp_agreement": {"class": "not_implemented_mcp"}}}
        run = _Stub(comp_entry=lambda sid: comp[sid], depth={}, planned_on=lambda sid, arm: sid != "G-w-0001")
        _, mtext = sections._agreement_counts(_Stub(run=run), list(comp))
        self.assertIn("1 not implemented", mtext)
        self.assertIn("1 with no MCP execution planned", mtext)

    def test_host_memory_in_evidence(self):
        import sections
        import task_resources
        self.assertIn(task_resources.HOST_MEMORY_FILE, sections.EVIDENCE_RUN_FILES)

    def test_stale_graphlet_inputs_named(self):
        import sections
        run = _Stub(gate_recorded={"judged_utc": "2026-10-21T12:40:00Z", "policy": {"name": "v8-rtp-mcp"},
                                   "passed": True}, gate_stale=["graphlet_inputs"])
        text = sections.stale_gate_text(run)
        self.assertIn("(policy v8-rtp-mcp, PASS) judged other contents of chapter G's inputs (graphlet_inputs: one "
                      "digest over", text)
        self.assertIn("the comparison's code", text)
        run.gate_stale = ["executions.jsonl"]
        self.assertIn("judged other contents of executions.jsonl than the run has now", sections.stale_gate_text(run))


class ProblemNotesV8Test(unittest.TestCase):
    """The v8 baseline (2026-10-21): the gate's one problem, N-realrep-0023 mcp r0, got no answer although the
    service's log shows it sent: the identity input's problem_notes put that note beside the execution in I.1."""

    NOTE = {"scenario": "N-x-0001", "arm": "mcp", "repeat": 0, "source": "the operator's log check",
            "text": "the answer left the service with status 200 after 3.1 s."}

    def _ctx(self, notes):
        rows = {"mcp": [{"scenario": "N-x-0001", "arm": "mcp", "repeat": 0, "outcome": "error"}]}
        comp = {"N-x-0001": {"mcp_agreement": {"class": "disagree", "differing_checks": [], "outcome_differs": []}}}
        ident = {"service_commit": "abc"}
        if notes is not None:
            ident["problem_notes"] = notes
        run = _Stub(order=["N-x-0001"], arm_lines=lambda sid, arm: rows.get(arm, []),
                    comp_entry=lambda sid: comp.get(sid, {}), comp_relations=lambda: {},
                    scenario=lambda sid: {"title": "a title"}, identity=ident)
        return _Stub(run=run, units=_Stub(by_scenario={}), chapters=[])

    def _items(self, notes):
        import sections
        groups = dict(sections.principal_findings(self._ctx(notes)))
        return " ".join(i for t in groups for i in groups[t])

    def test_note_beside_the_execution(self):
        items = self._items([self.NOTE])
        self.assertIn("mcp r0 error (no HTTP answer): the answer left the service with status 200 after 3.1~s (the "
                      "service's own log, checked separately; not part of the run's files)", items)
        self.assertIn("the outcomes differ; the MCP execution's error above, with its note", items)
        self.assertNotIn("operator", items)

    def test_no_note_no_text(self):
        for notes in (None, [], [dict(self.NOTE, repeat=1)], [dict(self.NOTE, arm="direct-local")]):
            items = self._items(notes)
            self.assertIn("mcp r0 error (no HTTP answer)", items)
            self.assertNotIn("own log", items)
            if notes != [dict(self.NOTE, repeat=1)]:
                # the disagreement names a note on any repeat of its arm
                self.assertNotIn("with its note", items)

    def test_note_must_name_an_execution_that_did_not_pass(self):
        import tempfile
        import reader
        rows = [{"scenario": "N-x-0001", "arm": "mcp", "repeat": 0, "outcome": "error",
                 "started_utc": "2026-10-21T12:02:38Z", "ended_utc": "2026-10-21T12:04:23Z"},
                {"scenario": "N-x-0001", "arm": "mcp", "repeat": 1, "outcome": "pass",
                 "started_utc": "2026-10-21T12:04:34Z", "ended_utc": "2026-10-21T12:04:37Z"}]
        stub = _Stub(rows=rows, run={"release": None}, mcp_run={}, service_identity=None, manifest={},
                     run_id=lambda: "R")
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "identity.json")
            for note, refusal in ((self.NOTE, None),
                                  (dict(self.NOTE, repeat=1), "no recorded execution that did not pass"),
                                  (dict(self.NOTE, scenario="N-y-0001"), "no recorded execution that did not pass"),
                                  (dict(self.NOTE, text=" "), "is not {scenario, arm, repeat, text"),
                                  (dict(self.NOTE, repeat=True), "is not {scenario, arm, repeat, text")):
                with open(path, "w") as f:
                    json.dump({"run": "R", "service_commit": "abc1234", "problem_notes": [note]}, f)
                if refusal is None:
                    self.assertEqual(reader.Run._load_identity(stub, path)["problem_notes"], [note])
                else:
                    with self.assertRaises(reader.RunError) as cm:
                        reader.Run._load_identity(stub, path)
                    self.assertIn(refusal, str(cm.exception))
            with open(path, "w") as f:
                json.dump({"run": "R", "service_commit": "abc1234", "problem_notes": {"x": 1}}, f)
            with self.assertRaises(reader.RunError):
                reader.Run._load_identity(stub, path)

    def test_notes_in_the_evidence(self):
        import sections
        ev = sections._identity_evidence(_Stub(service_identity={"recorded": False}, mcp_run={}),
                                         {"source": "s", "service_commit": "abc", "problem_notes": [self.NOTE]})
        self.assertEqual(ev["problem_notes"], [self.NOTE])
