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


if __name__ == "__main__":
    unittest.main()
