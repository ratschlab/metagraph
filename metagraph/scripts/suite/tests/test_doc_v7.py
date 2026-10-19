"""The v7 results document as a delta on v6 (2026-10-19): manifest v7 (content e89708cb...) is v6 with only the
scenarios' expectations and MCP arguments changed, so the unit mapping is v6's re-keyed, the prose's v7 statements are
`manifest_version >= 7` branches (a v5 / v6 build stays as it was), and appendix G lists the property-level marks
service-side by design. Offline; stdlib only."""

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

import facts  # noqa: E402
import sections  # noqa: E402
import units as units_mod  # noqa: E402

PROSE_DIR = os.path.join(DOC, "prose")
V7_CONTENT_SHA = "e89708cb795e098225d76ac8a69e75a116f9029549a8613c96325e8ee6a36238"


def prose(name):
    return sc.read_json(os.path.join(PROSE_DIR, "%s.json" % name))


class _Stub(object):
    def __init__(self, **kw):
        for k, v in kw.items():
            setattr(self, k, v)


class V7MappingTest(unittest.TestCase):
    def test_mapping_is_v6_rekeyed(self):
        v6 = sc.read_json(os.path.join(SUITE, "manifest", "suite_v6.json"))
        v7 = sc.read_json(os.path.join(SUITE, "manifest", "suite_v7.json"))
        self.assertEqual(v7["content_sha256"], V7_CONTENT_SHA)
        self.assertEqual(sc.scenarios_sha256(v7["scenarios"]), V7_CONTENT_SHA)
        self.assertEqual([s["id"] for s in v6["scenarios"]], [s["id"] for s in v7["scenarios"]])
        m6 = units_mod.load_mapping(os.path.join(DOC, "units", "suite_v6.json"))
        path7 = os.path.join(DOC, "units", "suite_v7.json")
        self.assertEqual(units_mod.default_mapping_path({"manifest_version": 7, "content_sha256": V7_CONTENT_SHA}),
                         path7)
        m7 = units_mod.load_mapping(path7)
        raw6 = sc.read_json(os.path.join(DOC, "units", "suite_v6.json"))
        raw7 = sc.read_json(path7)
        self.assertEqual(raw7["_manifest"], {"version": 7, "content_sha256": V7_CONTENT_SHA, "scenarios": 725})
        for key in raw6:
            if key not in ("_about", "_manifest"):
                self.assertEqual(raw7[key], raw6[key], key)
        units = units_mod.resolve(v7, m7)
        self.assertEqual(sorted(units.by_scenario), sorted(s["id"] for s in v7["scenarios"]))
        self.assertTrue(m6)


class V7ProseTest(unittest.TestCase):
    """Every v7 statement is a `manifest_version >= 7` branch; v6 renders its own text."""

    def test_branches(self):
        p = prose("P")
        intro = p["intro"][4]
        self.assertIn("carries the entry's `absence_scope` and `strands`",
                      facts.fill(intro, {"manifest_version": 7}, unresolved=[]))
        self.assertIn("carries neither `absence_scope` nor `timing.elapsed_ms`",
                      facts.fill(intro, {"manifest_version": 6}, unresolved=[]))
        note = p["units"]["P-edge"]["note"]
        self.assertIn("manifest v7 marks its entry_fields check service-side by design",
                      facts.fill(note, {"manifest_version": 7}, unresolved=[]))
        self.assertNotIn("the next manifest revision corrects it", facts.fill(note, {"manifest_version": 7},
                                                                              unresolved=[]))
        self.assertIn("the next manifest revision corrects it", facts.fill(note, {"manifest_version": 6},
                                                                           unresolved=[]))
        n = prose("N")["coverage_gap"][2]
        self.assertIn("Version 7 keeps them but for that mark", facts.fill(n, {"manifest_version": 7}, unresolved=[]))
        self.assertNotIn("Version 7", facts.fill(n, {"manifest_version": 6}, unresolved=[]))
        fails = prose("checks")["checks"]["entry_fields"]["fails"]
        self.assertTrue(facts.fill(fails, {"manifest_version": 6}, unresolved=[]).endswith(
            "whose answer carries no `absence_scope` and no `timing.elapsed_ms`)."))
        self.assertNotIn("absence_scope", facts.fill(fails, {"manifest_version": 7}, unresolved=[]).split(
            "access path, whose answer")[1])

    def test_catalogue_versioned(self):
        ctx = _Stub(run=_Stub(manifest={"manifest_version": 7}))
        self.assertEqual(sections._versioned(ctx, "a [[manifest_version >= 7 ? x | y]] b"), "a x b")
        self.assertEqual(sections._versioned(ctx, "keeps {braces} [[manifest_version >= 7 ? x | y]]"),
                         "keeps {braces} [[manifest_version >= 7 ? x | y]]")


class ByDesignMarksTest(unittest.TestCase):
    def test_table(self):
        marks = {"P-edge-0004": {"entry_fields": {"mark": "by_design", "case": "pattern.lowercase"}},
                 "T-labels-0003": {"label_dict_contains": {"mark": "not_implemented", "reason": "x"}}}
        run = _Stub(order=["P-edge-0004", "T-labels-0003"], check_marks=lambda sid: marks.get(sid, {}))
        out = sections.property_by_design_table(_Stub(run=run))
        self.assertIn("The manifest's property-level marks service-side by design", out)
        plain = out.replace("\\allowbreak{}", "")
        self.assertIn("\\texttt{pattern.lowercase} & P-edge-0004 & 1", plain)
        self.assertNotIn("T-labels-0003", plain)
        self.assertEqual(sections.property_by_design_table(_Stub(run=_Stub(order=["X"], check_marks=lambda s: {}))),
                         "")



class V7RealRunTest(unittest.TestCase):
    """What the v7 baseline's build showed: P-scope-0004 / -0006 pass via MCP (the sentence says so plainly), and a
    gate note naming a 40-hex commit and a path breaks within the margin."""

    def test_p_scope_both_pass(self):
        t = prose("P")["families"]["P-scope"]["prose"]
        ok = {"mcp": {"standing": "pass"}}
        filled = facts.fill(t, {"manifest_version": 7, "x": {"P-scope-0004": ok, "P-scope-0006": ok},
                                "n_relations": 0}, unresolved=[])
        self.assertIn("kept the two off the MCP access path): via MCP both pass.", filled)

    def test_gate_note_breaks(self):
        import inspect
        src = inspect.getsource(sections.gate_problems)
        self.assertIn('re.sub(r"[0-9a-f]{32,}"', src)
        self.assertIn("esc_wrap(latex.public_paths(str(n)), 300)", src)


class MemoryEpisodeTest(unittest.TestCase):
    """The v7 baseline's host memory episode (Timing 6) beside its time-budget issues: a window at the container's
    memory limit is read from samples.jsonl; the time-budget issues are grouped as budget_stress, a stated stop by
    the test's own budget, inside the window, just before it (time between the service and the host); none for a
    run whose memory never reached its limit."""

    def _run(self, peak=128):
        gib = float(1 << 30)
        samples = []
        for i, minute in enumerate(range(0, 60)):
            mem = (peak if 20 <= minute <= 40 else 80) * gib
            samples.append({"container": "eng", "mem_bytes": mem, "mem_limit_bytes": 128 * gib,
                            "blk_read": 1e9 * minute, "t_utc": "2026-10-19T02:%02d:00.000Z" % minute})
        def line(sid, arm, rep, minute, server, backend=None, outcome="pass", tl="", wall=None):
            l = {"scenario": sid, "arm": arm, "repeat": rep, "outcome": outcome, "http": 200, "time_limited": tl,
                 "started_utc": "2026-10-19T02:%02d:00.000Z" % minute, "server": {"elapsed_ms": server},
                 "wall_ms": wall if wall is not None else server}
            if backend is not None:
                l["mcp_phases"] = {"backend_ms": backend}
            return l
        rows = [line("M-1", "direct-local", 0, 25, 96000.0), line("M-1", "direct-local", 1, 50, 500.0),
                line("X-1", "mcp", 1, 10, 34.0, backend=66000.0, wall=70000.0),
                line("R-1", "direct-local", 0, 5, 300.0, tl="resolve stopped by time"),
                line("W-1", "direct-local", 0, 30, 60000.0)]
        scen = {"M-1": {"tags": []}, "X-1": {"tags": []}, "R-1": {"tags": []}, "W-1": {"tags": ["budget_stress"]}}
        return _Stub(samples=samples, rows=rows, order=list(scen), scenario=scen.get)

    def test_episode_and_text(self):
        ep = sections.memory_episode(self._run())
        self.assertEqual((ep["start"], ep["end"], ep["n"]), ("02:20:00", "02:40:00", 21))
        self.assertEqual(sorted(ep["inside"]), ["M-1"])
        self.assertEqual(sorted(ep["before"]), ["X-1"])
        self.assertEqual(ep["stops"], ["R-1"])
        self.assertEqual(ep["stress"], ["W-1"])
        text = sections.memory_episode_text(_Stub(run=self._run()))
        self.assertIn("stood at its memory limit of 128~GiB from 02:20:00 to 02:40:00 UTC", text)
        self.assertIn("M-1 (96~s server time on 1 execution; outside the window 0.5~s)", text)
        self.assertIn("likely the same episode, not established", text)
        self.assertIn("not a finding of the service", text)
        # RV8: Observation / What ran then / Limits, every memory figure in GiB (the disk read too)
        for head in ("\\emph{Observation.}", "\\emph{What ran then.}", "\\emph{Limits.}"):
            self.assertIn(head, text)
        self.assertLess(text.index("Observation."), text.index("What ran then."))
        self.assertLess(text.index("What ran then."), text.index("Limits."))
        self.assertIn("from disk during the window", text)
        self.assertNotIn("GB from disk", text.replace("GiB", ""))

    def test_service_log_note(self):
        run = self._run()
        run.identity = {"service_commit": "abc", "service_log": {"source": "a lane's CIRCUIT_X", "text": "two slow answers"}}
        text = sections.memory_episode_text(_Stub(run=run))
        # RV8: the log in plain words; its source (a file, an internal state) only in evidence.json
        self.assertIn("The service's own log recorded two slow answers (checked separately; not part of the run's "
                      "files)", text)
        self.assertNotIn("a lane", text)
        self.assertNotIn("CIRCUIT", text)
        self.assertNotIn("likely the same episode", text)
        ev = sections._identity_evidence(_Stub(service_identity={"recorded": False}, mcp_run={}), run.identity)
        self.assertEqual(ev["service_log"]["source"], "a lane's CIRCUIT_X")

    def test_identity_evidence(self):
        rec = {"recorded": True, "service_commits": {"abc": 2}, "service_commit": "abc", "identity_changed": False,
               "invocations": 2}
        ev = sections._identity_evidence(_Stub(service_identity=rec), {"source": "s", "service_commit": "abc"})
        self.assertEqual(ev["service_commit"], "abc")
        self.assertEqual(ev["input"]["source"], "s")
        ev6 = sections._identity_evidence(_Stub(service_identity={"recorded": False}, mcp_run={}),
                                          {"source": "log", "service_commit": "def"})
        self.assertEqual(ev6["service_commit"], "def")
        self.assertNotIn("input", ev6)

    def test_no_episode(self):
        self.assertIsNone(sections.memory_episode(self._run(peak=100)))
        self.assertEqual(sections.memory_episode_text(_Stub(run=self._run(peak=100))), "")


class SincePreviousTest(unittest.TestCase):
    """I.4's 'Since the previous edition' (RV7), computed from the previous edition's run given as an input: the
    scenarios, the MCP marks removed and moved, a check argument gained, the outcomes changed and the inconclusive
    count; nothing without that input (every v5 / v6 build)."""

    def _scen(self, sid, ni=None, ni_checks=(), bd_checks=(), reserve=None):
        mcp = {"tool": "pattern_search"}
        if ni:
            mcp["not_implemented"] = ni
        if ni_checks:
            mcp["not_implemented_checks"] = [{"check": c, "reason": "r"} for c in ni_checks]
        if bd_checks:
            mcp["by_design_checks"] = [{"check": c, "case": "pattern.lowercase"} for c in bd_checks]
        args = {"rule": "x"}
        if reserve is not None:
            args["reserve_ms"] = reserve
        return {"id": sid, "route": "POST /pattern", "mcp": mcp,
                "expect": [{"check": "entry_fields", "args": {}}, {"check": "budget_respected", "args": args}]}

    def test_block(self):
        prev_s = [self._scen("P-a-0001", ni="not implemented yet"), self._scen("P-a-0002", ni_checks=["entry_fields"]),
                  self._scen("P-a-0003", ni_checks=["entry_fields"])]
        cur_s = [self._scen("P-a-0001", reserve=250), self._scen("P-a-0002", reserve=250),
                 self._scen("P-a-0003", bd_checks=["entry_fields"], reserve=250)]
        def row(sid, arm, rep, outcome):
            return {"scenario": sid, "arm": arm, "repeat": rep, "outcome": outcome}
        prev = _Stub(manifest={"manifest_version": 6, "scenarios": prev_s}, scenarios=dict((s["id"], s) for s in prev_s),
                     rows=[row("P-a-0001", "mcp", 0, "not_implemented_mcp"), row("P-a-0002", "mcp", 0, "pass"),
                           row("P-a-0002", "direct-local", 0, "inconclusive")],
                     reevaluated=True, run={"run_id": "R6"})
        prev.run_id = lambda: "R6"
        prev.lines = lambda: dict(((r["scenario"], r["arm"], r["repeat"]), r) for r in prev.rows)
        cur = _Stub(manifest={"manifest_version": 7, "scenarios": cur_s}, order=[s["id"] for s in cur_s],
                    scenarios=dict((s["id"], s) for s in cur_s),
                    rows=[row("P-a-0001", "mcp", 0, "pass"), row("P-a-0001", "mcp", 1, "pass"),
                          row("P-a-0002", "mcp", 0, "pass"), row("P-a-0002", "direct-local", 0, "pass")],
                    previous=prev, service_identity={"recorded": True, "service_commits": {"abc1234": 3}}, mcp_run={})
        cur.scenario = cur.scenarios.get
        text = sections.since_previous(_Stub(run=cur))
        self.assertIn("The same 3 scenarios; 4 executions against 3", text)
        self.assertIn("removes 1 scenario-level mark not implemented via MCP (P-a-0001), so the test now runs via MCP",
                      text)
        self.assertIn("removes 1 property-level mark not implemented via MCP on entry\\_fields (P-a-0002)", text)
        self.assertIn("marks P-a-0003's entry\\_fields service-side by design instead of not implemented via MCP", text)
        self.assertIn("every MCP execution passed", text)
        self.assertIn("gains the argument reserve\\_ms on 3 scenarios (POST /pattern): a /pattern time stop is now "
                      "held to the host's documented finalisation reserve of 250 ms", text)
        self.assertIn("records the service's deployed commit", text)
        self.assertIn("not exercised via MCP to pass (P-a-0001 mcp r0)", text)
        self.assertIn("no execution to pass (P-a-0001 mcp r1)", text)
        self.assertIn("inconclusive to pass (P-a-0002 direct-local r0)", text)
        self.assertIn("0 inconclusive executions against 1", text)
        cur.previous = None
        self.assertEqual(sections.since_previous(_Stub(run=cur)), "")


class RecordedIdentityEvidenceTest(unittest.TestCase):
    def test_recorded_identity_in_evidence(self):
        run = _Stub(service_identity={"recorded": True, "service_commits": {"abc1234": 3}, "service_commit": "abc1234",
                                      "identity_changed": False, "invocations": 3})
        ev = sections._recorded_identity_evidence(run)
        self.assertEqual(ev["service_commit"], "abc1234")
        self.assertIn("recorded by the MCP runner", ev["source"])
        self.assertIsNone(sections._recorded_identity_evidence(_Stub(service_identity={"recorded": False},
                                                                     mcp_run={})))


if __name__ == "__main__":
    unittest.main()
