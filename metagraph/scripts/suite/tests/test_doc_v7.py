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

import gate  # noqa: E402
import suite_common as sc  # noqa: E402

import facts  # noqa: E402
import sections  # noqa: E402
import units as units_mod  # noqa: E402

PROSE_DIR = os.path.join(DOC, "prose")
V7_CONTENT_SHA = "e89708cb795e098225d76ac8a69e75a116f9029549a8613c96325e8ee6a36238"


def latex_plain(text):
    return text.replace("*", "")


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
        # RV4: one reason, told plainly: the answer carries no pattern text; the case documents the upper-casing
        # no code span: the unit note is set in italics, and Fira Mono has no italic (e29e8edc0)
        self.assertIn("via MCP the entry_fields check is set aside as service-side by design (manifest v7's "
                      "mcp.by_design_checks", facts.fill(note, {"manifest_version": 7}, unresolved=[]))
        self.assertIn("the tool's answer carries no pattern text", facts.fill(note, {"manifest_version": 7},
                                                                             unresolved=[]))
        self.assertNotIn("`", facts.fill(note, {"manifest_version": 7}, unresolved=[]))
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
        v7 = facts.fill(fails, {"manifest_version": 7}, unresolved=[])
        # V8: what entry_fields reads via MCP under v7: P-edge-0004's pattern text, never timing.elapsed_ms
        self.assertNotIn("absence_scope", v7)
        self.assertNotIn("timing.elapsed_ms", v7)
        self.assertIn("only P-edge-0004's pattern text, set aside as service-side by design", v7)
        args = prose("checks")["checks"]["budget_respected"]["args"]
        self.assertIn("`reserve_ms` (manifest v7", facts.fill(args, {"manifest_version": 7}, unresolved=[]))
        self.assertNotIn("reserve_ms", facts.fill(args, {"manifest_version": 6}, unresolved=[]))

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
        # V12: the service and the MCP runner provide the field; manifest v7 only checks it now
        self.assertIn("the service states it per database and the MCP runner (git head {mcp_runner_head}) copies the "
                      "host's values onto the entry, and manifest v7 checks it via MCP", str(prose("P")))

    def test_gate_note_breaks(self):
        # V6: the behaviour, not the source: a note with a 40-hex commit breaks every 8 hex digits, a local path is
        # made public, the LaTeX special characters are escaped
        sha = "26a52cb63efd95752a109b5e09b1163bfffe60ca"
        tex = sections.gate_note_tex("the MCP arm recorded the service commit %s on 1369 invocation(s) (%s)" % (
            sha, os.path.join(SUITE, "runs", "R", "arms", "mcp", "mcp_runs.jsonl")))
        self.assertTrue(SUITE.startswith("/"))
        self.assertIn("26a52cb6\\allowbreak{}3efd9575\\allowbreak{}2a109b5e\\allowbreak{}09b1163b\\allowbreak{}"
                      "fffe60ca", tex)
        self.assertNotIn(sha, tex)
        self.assertNotIn(SUITE, tex.replace("\\allowbreak{}", ""))
        self.assertIn("mcp\\_runs", tex.replace("\\allowbreak{}", ""))
        self.assertTrue(tex.startswith("{\\small\\color{black!60}note: "))
        # W7: counts in their number
        self.assertIn("plan: 3015 executions from plan.json", sections.gate_note_tex("plan: 3015 execution(s) from "
                                                                                     "plan.json"))
        self.assertEqual(sections._break_hex("abc 0123456789abcdef0123456789abcdef"),
                         "abc 01234567\\allowbreak{}89abcdef\\allowbreak{}01234567\\allowbreak{}89abcdef")


class MemoryEpisodeTest(unittest.TestCase):
    """The host's memory episodes (Timing 6), as the reviews of the v7 document require (RV2, RV8, V1, V2, V4, V6):
    each contiguous run of samples at the container's memory limit is an episode (two print as two); an execution is
    placed by its [start, end] at millisecond precision before any tag (a budget_stress 503 inside the window is named
    as the one that did not pass); an execution running across a drop at the window's end is told apart, with the
    test's next execution; the MCP executions the host answered late are found by rule (host call against the host's
    own search, the container at 90 % of its limit or more); the service's log is printed with its times, each tied
    to the host call it falls in; every 503 and inconclusive execution is placed in or out of the window; nothing
    without an episode."""

    GIB = float(1 << 30)

    def _samples(self, two=False, peak=128):
        out = []
        t = 0
        while t <= 3600:
            hh, mm, ss_ = 2, t // 60, t % 60
            minute = t / 60.0
            if 20 * 60 <= t <= 30 * 60:
                mem = peak
            elif 30 * 60 < t:
                mem = 51 + 0.01 * ((t - 30 * 60) / 2.0)
                if two and 45 * 60 <= t <= 46 * 60:
                    mem = peak
            elif 5 <= minute < 20:
                mem = 120
            else:
                mem = 80
            out.append({"container": "eng", "mem_bytes": min(mem, 128) * self.GIB, "mem_limit_bytes": 128 * self.GIB,
                        "blk_read": 1e9 * t, "t_utc": "2026-10-19T%02d:%02d:%02d.000Z" % (hh, mm, ss_)})
            t += 2
        return out

    def _run(self, two=False, peak=128):
        def line(sid, arm, rep, start, end, server=None, outcome="pass", http=200, tl="", wall=None, phases=None,
                 error=None):
            l = {"scenario": sid, "arm": arm, "repeat": rep, "outcome": outcome, "http": http, "time_limited": tl,
                 "started_utc": "2026-10-19T%s" % start, "ended_utc": "2026-10-19T%s" % end,
                 "wall_ms": wall if wall is not None else server}
            if server is not None:
                l["server"] = {"elapsed_ms": server}
            if phases:
                l["mcp_phases"] = phases
            if error:
                l["error"] = error
            return l
        rows = [line("M-1", "direct-local", 0, "02:25:00.000Z", "02:26:36.000Z", 96000.0),
                line("M-1", "direct-local", 1, "02:50:00.000Z", "02:50:00.500Z", 500.0),
                line("W-1", "direct-local", 0, "02:27:00.000Z", "02:28:50.800Z", None, outcome="inconclusive",
                     http=503, wall=110800.0, tl="a 503 deadline", error={"code": "deadline"}),
                line("B-1", "mcp", 0, "02:29:59.500Z", "02:31:20.000Z", 80000.0, wall=80500.0),
                line("B-1", "direct-local", 0, "02:32:00.000Z", "02:32:00.600Z", 600.0),
                line("X-1", "mcp", 1, "02:10:00.000Z", "02:11:10.000Z", 34.0, wall=70000.0,
                     phases={"backend_ms": 66000.0, "backend_pattern_ms": 34.0, "submit_ms": 200.0, "queue_ms": 200.0}),
                line("R-1", "direct-local", 0, "02:03:00.000Z", "02:03:00.300Z", 300.0, tl="resolve stopped by time"),
                line("O-1", "direct-local", 0, "02:52:00.000Z", "02:52:01.000Z", 900.0, outcome="inconclusive")]
        scen = {"M-1": {"tags": [], "area": "M"}, "W-1": {"tags": ["budget_stress"], "area": "P"},
                "B-1": {"tags": [], "area": "M"}, "X-1": {"tags": [], "area": "X"}, "R-1": {"tags": [], "area": "R"},
                "O-1": {"tags": [], "area": "C"}}
        return _Stub(samples=self._samples(two, peak), rows=rows, order=list(scen), scenario=scen.get,
                     run={"server_instance": "abc"}, previous=None)

    def _text(self, run):
        return sections.memory_episode_text(_Stub(run=run, units=_Stub(chapters=["R", "P", "X", "M", "C"])))

    def test_episode_and_text(self):
        run = self._run()
        eps = sections.memory_episodes(run)
        self.assertEqual(len(eps), 1)
        ep = eps[0]
        self.assertEqual((ep["start"], ep["end"], ep["n"]), ("02:20:00", "02:30:00", 301))
        pos = dict(((l["scenario"], l["arm"], l["repeat"]), p) for l, p in ep["rows"])
        self.assertEqual(pos, {("M-1", "direct-local", 0): "inside", ("W-1", "direct-local", 0): "inside",
                               ("B-1", "mcp", 0): "across"})
        self.assertEqual(round(ep["drop"]["to_gib"]), 51)
        self.assertTrue(ep["drop"]["in_flight_completed"])
        text = self._text(run)
        self.assertIn("stood at its memory limit of 128~GiB from 02:20:00 to 02:30:00 UTC", text)
        for head in ("\\emph{Observation.}", "\\emph{What ran then.}", "\\emph{Limits.}"):
            self.assertIn(head, text)
        self.assertLess(text.index("not a finding of the service"), text.index("Observation."))
        # V7-001: observations, not a cause
        self.assertIn("At the end of the window, reported container memory fell from 128 to 51~GiB within 2~s; the "
                      "in-flight request completed. These observations do not establish the cause of the drop.", text)
        for claim in ("cached index pages", "no restart", "restarted"):
            self.assertNotIn(claim, text)
        # W1: the one observed fact beside the drop, only when run.json shows one instance throughout
        self.assertNotIn("server instance", text)
        run.run = {"server_instance": "abc", "arms": {"direct-local": {
            "server_instance": "abc", "server_instance_end": "abc", "server_instance_changed": False}}}
        self.assertIn("the in-flight request completed; the server instance was the same before and after. These "
                      "observations do not establish the cause of the drop.", self._text(run))
        run.run["arms"]["direct-local"].update(server_instance_end="def", server_instance_changed=True)
        self.assertNotIn("server instance", self._text(run))
        self.assertFalse(sections.same_server_instance(run))
        run.run = {"server_instance": "abc"}
        self.assertIn("The latency comparisons use the same test's other executions; they do not isolate the cause "
                      "of the slowdown.", text)
        self.assertNotIn("not measured", text)
        # V1: the budget_stress 503 inside the window is placed there and named as the one that did not pass
        self.assertIn("3 executions ran in the window (chapters P and M; one ran past its end): 2 passed, and one did "
                      "not: W-1 direct-local r0 (a 503 deadline at 110.8~s wall, inconclusive)", text)
        self.assertNotIn("Every one of these executions passed", text)
        self.assertIn("W-1 (budget\\_stress; direct-local r0: the 503 deadline above)", text)
        self.assertIn("M-1 (96~s server time on 1 execution; outside the window 0.5~s)", text)
        # V2: the execution that started in the window's last second and ran across the fall
        self.assertIn("B-1 (mcp r0 started 0.5~s before the window's end and ran across the fall to 51~GiB: 80~s "
                      "server time; its next execution, direct-local r0 at 02:32:00 UTC, took 0.6~s)", text)
        # V4: the late host answer by rule, and the stated stop computed from its outcome
        self.assertIn("X-1 mcp r1 at 02:10:00 UTC, the container at 93~\\% of its limit", text)
        self.assertIn("1 ended by a time budget they state (R-1), as their tests intend or allow: each of these "
                      "executions passed", text)
        # RV2: every 503 and inconclusive placed against the window
        self.assertIn("inside, W-1 direct-local r0 (above); outside, and not attributed to it: O-1 (direct-local r0: "
                      "inconclusive; 02:52:00 to 02:52:01 UTC", text)
        self.assertNotIn("likely the same episode", text)

    def test_pointers(self):
        # RV1, V9: I.1's sentence, Timing 1's observation and the chapters' Timing line, computed from the episode
        run = self._run()
        units = _Stub(chapters=["R", "P", "X", "M", "C"],
                      chapter_scenarios=lambda ch: [s for s in run.order if run.scenario(s)["area"] == ch])
        ctx = _Stub(run=run, units=units)
        lead = sections.memory_lead_text(ctx)
        self.assertIn("From 02:20 to 02:30 UTC the engine's container stood at its 128 GiB memory limit on the host: 3 "
                      "executions ran then (chapters P and M); 2 passed, and W-1 direct-local r0 ended in a 503 "
                      "deadline, inconclusive. 2 of the 3 searches in that window took at least twice as long as the "
                      "same test elsewhere (B-1, M-1, 133 to 192 times); the others did not (Timing 6). It is an "
                      "operational note about the host, not a finding of the service.", latex_plain(lead))
        # W2: a search in the window faster than the same test elsewhere is not counted as slower
        run2 = self._run()
        run2.scenario.__self__["F-1"] = {"tags": [], "area": "M"}
        run2.order.append("F-1")
        run2.rows += [dict(run2.rows[0], scenario="F-1", server={"elapsed_ms": 50.0}, wall_ms=50.0),
                      dict(run2.rows[1], scenario="F-1", server={"elapsed_ms": 5000.0}, wall_ms=5000.0)]
        m = sections.memory_summary(_Stub(run=run2, units=units))[0]
        self.assertEqual(sorted(l["scenario"] for l, _ in m["slower"]), ["B-1", "M-1"])
        self.assertIn("2 of the 4 searches in that window took at least twice as long", sections._memory_slower_text(m))
        self.assertIn("placed in Timing 6: 3, of which 2 passed", sections.memory_timing_bullet(ctx))
        self.assertEqual(sections.memory_chapter_note(ctx, "M"), "2 of this chapter's executions ran while the "
                         "engine's container stood at its memory limit, 02:20 to 02:30 UTC (Timing 6)")
        self.assertEqual(sections.memory_chapter_note(ctx, "X"), "The host answered X-1 mcp r1 late, with the "
                         "container at 93 % of its memory limit (Timing 6)")
        self.assertEqual(sections.memory_chapter_note(ctx, "R"), "")
        quiet = _Stub(run=self._run(peak=100), units=units)
        self.assertEqual((sections.memory_lead_text(quiet), sections.memory_timing_bullet(quiet),
                          sections.memory_chapter_note(quiet, "M")), ("", "", ""))

    def test_two_episodes(self):
        eps = sections.memory_episodes(self._run(two=True))
        self.assertEqual([(e["start"], e["end"]) for e in eps], [("02:20:00", "02:30:00"), ("02:45:00", "02:46:00")])
        text = self._text(self._run(two=True))
        self.assertIn("In a second window, the engine's container (eng) stood at its memory limit of 128~GiB from "
                      "02:45:00 to 02:46:00 UTC", text)
        self.assertIn("In the window of 02:20:00 to 02:30:00, 3 executions ran (chapters P and M", text)
        self.assertNotIn("ran in the window", text)

    def test_service_log_note(self):
        run = self._run()
        run.identity = {"service_commit": "abc", "service_log": {
            "source": "a lane's CIRCUIT_X", "times": ["2026-10-19T02:10:05Z", "2026-10-19T02:12:00Z"],
            "text": "two requests to the host that timed out after 5 s, at {times}, while the host's status page still "
                    "answered; the log was otherwise clean"}}
        text = self._text(run)
        # RV8, V4, W4: the log in its own words with its times; its source (a file, an internal state) only in
        # evidence.json
        self.assertIn("The service's own log recorded two requests to the host that timed out after 5~s, at 02:10:05 "
                      "and 02:12:00 UTC, while the host's status page still answered; the log was otherwise clean "
                      "(checked separately; not part of the run's files).", text)
        run.identity["service_log"]["text"] = "a slow answer"
        self.assertIn("The service's own log recorded a slow answer, at 02:10:05 and 02:12:00 UTC (checked separately",
                      self._text(run))
        self.assertNotIn("a lane", text)
        self.assertNotIn("CIRCUIT", text)
        ev = sections._identity_evidence(_Stub(service_identity={"recorded": False}, mcp_run={}), run.identity)
        self.assertEqual(ev["service_log"]["source"], "a lane's CIRCUIT_X")
        self.assertEqual(ev["service_log"]["times"], ["2026-10-19T02:10:05Z", "2026-10-19T02:12:00Z"])

    def test_service_log_time_outside_the_run_refused(self):
        import tempfile
        import reader
        run = self._run()
        stub = _Stub(rows=run.rows, run={"release": None}, mcp_run={}, service_identity=None, manifest={},
                     run_id=lambda: "R")
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "identity.json")
            for times, ok in ((["2026-10-19T02:10:05Z"], True), (["2026-10-19T05:10:05Z"], False)):
                with open(path, "w") as f:
                    json.dump({"run": "R", "service_commit": "abc1234",
                               "service_log": {"text": "x", "times": times, "source": "s"}}, f)
                if ok:
                    self.assertEqual(reader.Run._load_identity(stub, path)["service_commit"], "abc1234")
                else:
                    with self.assertRaises(reader.RunError) as cm:
                        reader.Run._load_identity(stub, path)
                    self.assertIn("lies outside the run", str(cm.exception))

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
        # a run whose memory never reached its limit (every v5 / v6 run): nothing printed, nothing pointed to
        self.assertEqual(sections.memory_episodes(self._run(peak=100)), [])
        self.assertEqual(self._text(self._run(peak=100)), "")
        self.assertEqual(sections.memory_episodes(_Stub(samples=[], rows=[], run={})), [])


class DeployedIdentityTest(unittest.TestCase):
    """The service identity held to the deployed commit (RV3, V3): the identity input's `deployed_service_commit`
    (from the deployment record, never the run's own records) is validated, passed to both gate calls, and appendix
    F's note says in plain words that the recorded commit equals it; without it the recorded identity is said to be
    self-reported; the by-design table is tied to the service that answered."""

    SHA = "26a52cb63efd95752a109b5e09b1163bfffe60ca"

    def _gate(self, expected):
        rec = {"recorded": True, "service_commits": {self.SHA: 5}, "service_commit": self.SHA, "invocations": 5,
               "identity_changed": False, "source": "arms/mcp/mcp_runs.jsonl"}
        _, notes = gate.service_identity_problems(rec, expected, ["mcp"])
        return {"service_identity": {"expected": expected, "recorded": rec}, "policy": {"arms": ["mcp"]},
                "notes": notes + ["plan: 5 execution(s)"]}

    def test_note(self):
        ident = {"deployed_service_commit": {"commit": self.SHA, "source": "the deployment record of the 18 October deploy",
                                             "railway_deployment_id": "6afb760e"}}
        run = _Stub(gate=self._gate(self.SHA), identity=ident)
        text, replaced = sections.identity_note_text(_Stub(run=run))
        self.assertIn("equals the deployed commit", text)
        self.assertIn("(the deployment record of the 18 October deploy, Railway deployment 6afb760e) on all 5 MCP invocations", text)
        self.assertIn("\\allowbreak{}", text)
        self.assertNotIn("--expect", text)
        self.assertNotIn("mcp\\_runs", text)
        self.assertEqual(len(replaced), 1)
        run = _Stub(gate=self._gate(None), identity=None)
        text, replaced = sections.identity_note_text(_Stub(run=run))
        self.assertIn("self-reported, not cross-checked", text)
        # W6: an expected commit that came with the gate's record never credits an identity input not given
        run = _Stub(gate=self._gate(self.SHA), identity=None)
        text, _ = sections.identity_note_text(_Stub(run=run))
        self.assertIn("(as the gate's record names it) on all 5 MCP invocations", text)
        self.assertNotIn("identity input", text)
        # a run recording no identity keeps the gate's own note (v5, v6)
        run = _Stub(gate={"service_identity": {"expected": None, "recorded": {"recorded": False}}}, identity=None)
        self.assertEqual(sections.identity_note_text(_Stub(run=run)), ("", []))

    def test_gate_box_clause(self):
        # W8: the gate box and F's gate paragraph say the identity was held, only when it was
        self.assertTrue(sections.identity_held(_Stub(gate=self._gate(self.SHA))))
        self.assertFalse(sections.identity_held(_Stub(gate=self._gate(None))))
        self.assertFalse(sections.identity_held(_Stub(gate=self._gate("0" * 40))))
        self.assertEqual(sections.IDENTITY_HELD_CLAUSE,
                         "; the recorded service commit was held to the deployed commit (I.4)")

    def test_identity_input_validated(self):
        import tempfile
        import reader
        stub = _Stub(rows=[], run={"release": None}, mcp_run={}, service_identity=None, manifest={},
                     run_id=lambda: "R")
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "identity.json")
            for dep, ok in (({"commit": self.SHA, "source": "LIVE"}, True), ({"commit": "not-a-sha", "source": "x"},
                                                                           False), ({"commit": self.SHA}, False)):
                with open(path, "w") as f:
                    json.dump({"run": "R", "service_commit": self.SHA, "deployed_service_commit": dep}, f)
                if ok:
                    self.assertEqual(reader.Run._load_identity(stub, path)["deployed_service_commit"]["commit"],
                                     self.SHA)
                else:
                    with self.assertRaises(reader.RunError):
                        reader.Run._load_identity(stub, path)

    def test_reader_passes_the_expected_commit_to_both_gates(self):
        import shutil
        import tempfile
        import compare as cmp
        import gate_policy
        import merge_run
        import reader
        import test_doc as td
        tmp = tempfile.mkdtemp()
        try:
            run_dir = os.path.join(tmp, "run")
            os.makedirs(run_dir)
            td.build_run(run_dir)
            merge_run.merge(run_dir)
            sc.write_json_atomic(os.path.join(run_dir, "compare.json"), cmp.compare(run_dir, (td.CASES, {"cases": 5})))
            base = reader.Run(run_dir)
            policy = gate_policy.build_policy(base.manifest, table=(td.CASES, {"cases": 5}), manifest_file="m.json",
                                              arms=("direct-local", "mcp"), warm_every=0)
            ppath, ipath = os.path.join(tmp, "policy.json"), os.path.join(tmp, "identity.json")
            with open(ppath, "w") as fh:
                json.dump(policy, fh)
            with open(ipath, "w") as fh:
                json.dump({"run": base.run_id(), "service_commit": self.SHA,
                           "deployed_service_commit": {"commit": self.SHA, "source": "LIVE"}}, fh)
            seen, orig = [], reader.gate_mod.judge

            def spy(*a, **kw):
                seen.append(kw.get("expect_service_commit"))
                return orig(*a, **kw)
            reader.gate_mod.judge = spy
            try:
                run = reader.Run(run_dir, policy=ppath, also_policies=[ppath], identity=ipath)
            finally:
                reader.gate_mod.judge = orig
            self.assertEqual(seen, [self.SHA, self.SHA])
            self.assertEqual(run.gate["service_identity"]["expected"], self.SHA)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    def test_by_design_tie(self):
        comp = {"A": {"mcp_agreement": {"class": "service_by_design", "verified": True}},
                "B": {"mcp_agreement": {"class": "service_by_design"}}, "C": {"mcp_agreement": {"class": "agree"}}}
        run = _Stub(table_meta={"service_commit": "22d034ba5c3de328", "generated_utc": "2026-10-13T11:04:20Z"},
                    service_identity={"recorded": True, "service_commits": {self.SHA: 3}, "service_commit": self.SHA},
                    order=["A", "B", "C"], comp_entry=lambda sid: comp[sid],
                    identity={"deployed_service_commit": {"commit": self.SHA, "source": "LIVE"}})
        # W6: a pair without a `verified` field counts as unverified
        self.assertEqual(sections.by_design_tie_text(run), "generated at service commit 22d034ba (2026-10-13); 1 of 2 "
                         "by-design pairs verified against the deployed 26a52cb6")
        comp["B"]["mcp_agreement"]["verified"] = True
        self.assertEqual(sections.by_design_tie_text(run), "generated at service commit 22d034ba (2026-10-13); all 2 "
                         "by-design pairs verified against the deployed 26a52cb6")
        run.identity = None
        self.assertIn("against the recorded service 26a52cb6", sections.by_design_tie_text(run))
        run.service_identity = {"recorded": False}
        run.mcp_run = {}
        self.assertEqual(sections.by_design_tie_text(run), "")


class ByDesignPropertyMarkTest(unittest.TestCase):
    """The property-level mark service-side by design is visible and told apart from the not-implemented kind (RV4,
    V10): the MCP cell carries its own raised mark (a blue diamond, \\gbdp), the legend a row of its own (only for a
    manifest with such marks), and the chapter's 'fully checked' names the check set aside by design."""

    def _run(self, kind):
        marks = {"entry_fields": {"mark": kind, "case" if kind == "by_design" else "reason": "x"}}
        return _Stub(check_marks=lambda sid: marks,
                     marked_aside=lambda sid, l: marks if l.get("aside") else {}, order=["P-edge-0004"])

    def test_cell_mark(self):
        import cells
        run = self._run("by_design")
        self.assertTrue(cells.property_by_design(run, "P-edge-0004", [{"aside": True}]))
        self.assertFalse(cells.property_by_design(run, "P-edge-0004", [{}]))
        self.assertFalse(cells.property_not_implemented(run, "P-edge-0004", [{"aside": True}]))
        st = cells.ArmState("mcp")
        st.glyph_key, st.bd_marked = "pass", True
        self.assertTrue(st.glyph_tex().endswith("\\gbdp"))
        self.assertFalse(cells.property_by_design(self._run("not_implemented"), "P-edge-0004", [{"aside": True}]))
        with open(os.path.join(DOC, "templates", "doc.sty"), encoding="utf-8") as fh:
            self.assertIn("\\newcommand\\gbdp", fh.read())

    def test_legend_and_full_definition(self):
        self.assertTrue(sections._has_property_marks(self._run("by_design"), "by_design"))
        self.assertFalse(sections._has_property_marks(self._run("by_design")))
        self.assertIn("by-design property mark", sections.MARK_EXAMPLE)
        text = facts.mcp_full_definition(_Stub(run=self._run("by_design")), ["P-edge-0004"])
        self.assertIn("except P-edge-0004's entry_fields, which appendix G lists as set aside service-side by design",
                      text)
        self.assertNotIn("not compared and unverified", text)
        text = facts.mcp_full_definition(_Stub(run=self._run("not_implemented")), ["P-edge-0004"])
        self.assertIn("except a check a property-level mark of the manifest sets aside, not compared and unverified "
                      "via MCP (on one of them: P-edge-0004)", text)


class GatePoliciesNamedTest(unittest.TestCase):
    """The two gate policies explained where they first appear (RV5): the deploy gate's short name with its routes
    after its policy file, the other policy as the narrower one with its routes and the route it leaves out, its
    judged executions split by outcome as the deploy gate's are, and the glossary's short names."""

    def test_texts(self):
        rtp = {"name": "v7-rtp-mcp", "arms": ["direct-local", "mcp"], "selection": {"scenarios": ["A", "B"]},
               "counts": {"by_route": {"POST /pattern": 1, "POST /resolve": 1, "GET /traverse/capabilities": 1,
                                       "POST /traverse": 1}}, "exclusions": []}
        rt = {"name": "v7-rt-mcp", "arms": ["mcp"], "selection": {"scenarios": ["A"]},
              "counts": {"by_route": {"POST /traverse": 1, "POST /resolve": 1}}, "exclusions": []}
        self.assertEqual((sections.policy_label(rtp), sections.policy_label(rt)), ("rtp", "rt"))
        self.assertEqual(sections.policy_routes(rtp), ["POST /resolve", "POST /traverse",
                                                       "GET /traverse/capabilities", "POST /pattern"])
        rows = [{"scenario": "A", "arm": "mcp", "repeat": 0, "outcome": "pass"},
                {"scenario": "A", "arm": "mcp", "repeat": 1, "outcome": "service_by_design"}]
        entry = {"path": "/x/policy_v7_rt_mcp.json", "policy": rt,
                 "verdict": {"passed": True, "problems": [], "policy": {"name": "v7-rt-mcp", "policy_version": 1,
                                                                        "sha256": "ab" * 32},
                             "counts": {"judged": 2, "passed": 1, "required": 2}}}
        run = _Stub(policy=rtp, also_gates=[entry], rows=rows)
        text = sections.also_gate_text(entry, run)
        self.assertTrue(text.startswith("the narrower policy rt (POST /resolve and POST /traverse, without GET "
                                        "/traverse/capabilities and POST /pattern): policy_v7_rt_mcp.json"), text)
        self.assertIn("2 executions judged (1 passed, 1 service-side by design), 2 required", text)
        self.assertEqual(sections.deploy_label_text(run), " (rtp: POST /resolve, POST /traverse, GET "
                         "/traverse/capabilities and POST /pattern)")
        self.assertIn("rtp, the deploy gate (POST /resolve", sections.gate_policy_names(run))
        self.assertEqual(sections.deploy_label_text(_Stub(policy=rtp, also_gates=[])), "")
        self.assertEqual(sections.gate_policy_names(_Stub(policy=rtp, also_gates=[])), "")


class RecordedSetTest(unittest.TestCase):
    """The recorded set defined and the MCP runner named (RV6), and no sentence about a re-evaluation in a document
    of a recorded set (RV6, V12; v5 too)."""

    def _run(self, reevaluated=False):
        arms = {"direct-local": {"scripts_commit": "1a9cee291e4b409e", "checks_sha256": "7580" * 16},
                "mcp": {"checks_sha256": "7580" * 16}}
        return _Stub(reevaluated=reevaluated, reeval={}, run={"arms": arms}, mcp_run={"git_head": "fcf7ae49465c9b61"},
                     identity={"service_commit": "x", "mcp_runner_repository": "the search service's repository"},
                     gate={"passed": True}, gate_offline=True, policy_path="/x/policy.json",
                     executions_file="executions.jsonl")

    def test_texts(self):
        run = self._run()
        self.assertEqual(sections.checks_commit_text(run), "1a9cee291e4b (the runners' scripts, as recorded; the MCP "
                         "runner, a program at its own git head fcf7ae49465c, in the search service's repository, "
                         "judged with the same checks.py, sha256 758075807580…)")
        self.assertEqual(sections.mcp_runner_text(run), "the MCP client, git head fcf7ae49465c, in the search "
                         "service's repository")
        self.assertIn("compares it with the checks the runners recorded", sections.gate_read_text(run))
        self.assertNotIn("re-evaluation", sections.gate_read_text(run))
        common = " ".join(prose("checks")["common"])
        self.assertIn("the direct runners[[reevaluated ? and the re-evaluation]] record", common)
        self.assertIn("(the MCP runner's rule[[reevaluated ?, replayed by the re-evaluation]])", common)
        for p in prose("checks")["common"]:
            filled = facts.fill(p, {"reevaluated": False}, unresolved=[])
            self.assertNotIn("re-evaluation", filled)


class NotCoveredTest(unittest.TestCase):
    """RV14, W3: I.2 says plainly, for manifest v7, which MCP tools the scenarios name (counted from the manifest),
    which the runner called (read from its raw records), how the graphlet document is read, which traversal tools
    were never called, and that the MCP's other tools are outside the suite's scope."""

    def test_tools_and_paragraph(self):
        import gzip
        import shutil
        import tempfile
        import reader
        v7 = sc.read_json(os.path.join(SUITE, "manifest", "suite_v7.json"))
        scen = dict((s["id"], s) for s in v7["scenarios"])
        run = _Stub(order=[s["id"] for s in v7["scenarios"]], scenario=scen.get)
        self.assertEqual(sections.mcp_tools_text(run), "four tools: pattern_search (456 scenarios), traverse_resolve "
                         "(118), traverse_fetch (99) and catalog (1)")
        # the runner's raw records: tools and the export's format; a download and an HTTP read are no MCP tools
        tmp = tempfile.mkdtemp()
        try:
            raw = os.path.join(tmp, "arms", "mcp", "raw")
            os.makedirs(raw)
            calls = [{"tool": "traverse_fetch"}, {"tool": "traverse_status"}, {"tool": "graphlet_export",
                                                                               "args": {"format": "mgt"}},
                     {"tool": "download"}, {"tool": "GET /admin/traversal/jobs/x"}, {"tool": "catalog"}]
            with gzip.open(os.path.join(raw, "T-1.mcp.0.json.gz"), "wt") as fh:
                json.dump({"calls": calls}, fh)
            stub = _Stub(dir=tmp, _mcp_tools=None)
            got = reader.Run.mcp_tool_calls(stub)
            self.assertEqual(sorted(got), ["catalog", "graphlet_export", "traverse_fetch", "traverse_status"])
            self.assertEqual(got["graphlet_export"]["formats"], {"mgt": 1})
            stub.mcp_tool_calls = lambda: got
            self.assertEqual(sections.mcp_tools_called_text(stub), "`traverse_fetch`, `traverse_status`, `catalog` and "
                             "`graphlet_export` in the MGT format to download the graphlet document a traversal returns")
            text = sections.traversal_tools_not_called_text(stub)
            self.assertIn("the graphlet processing tools (`graphlet_summary`", text)
            self.assertIn("`graphlet_export`'s other formats)", text)
            self.assertTrue(text.endswith("`traverse_resolve`, `traverse_results`, `traverse_stop`, `traverse_lock` and "
                                          "`traverse_capabilities`"), text)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
        par = [p for p in prose("overview")["sections"]["methods"]["paragraphs"] if "*Not covered.*" in p][0]
        self.assertEqual(facts.fill(par, {"manifest_version": 6, "mcp_tools_text": "x", "mcp_tools_called": "y",
                                          "traversal_tools_not_called": "z"}, unresolved=[]).strip(), "")
        text = facts.fill(par, {"manifest_version": 7, "mcp_tools_text": "four tools", "mcp_tools_called": "TOOLS",
                                "traversal_tools_not_called": "NOT"}, unresolved=[])
        self.assertIn("to run them the MCP runner called, as its records show, TOOLS.", text)
        self.assertIn("Not called, and so not covered: NOT;", text)
        self.assertIn("reference examples of the document format", text)
        self.assertNotIn("codec vectors", text)
        self.assertIn("The MCP's other tools (sequence search, taxonomy, reference and sequence utilities) are outside "
                      "this suite's scope.", text)
        access = prose("overview")["sections"]["methods"]["paragraphs"][0]
        self.assertIn("`pattern_search`, `catalog`)", access)
        self.assertNotIn("the capabilities tools", access)


class PluralsTest(unittest.TestCase):
    """RV12: the headline boxes print counts in their number, never report.py's "(s)"."""

    def test_plural_counts(self):
        self.assertEqual(sections.plural_counts("FAIL (0 failed or errored execution(s), 1 direct and 0 MCP "
                                                "disagreement(s)); 34 inconclusive execution(s)"),
                         "FAIL (0 failed or errored executions, 1 direct and 0 MCP disagreements); 34 inconclusive "
                         "executions")
        self.assertEqual(sections.plural_counts("1 relation(s) failed; 1 inconclusive execution(s)"),
                         "1 relation failed; 1 inconclusive execution")
        self.assertEqual(sections._count_phrase(0, "execution", "none in this run"), "none in this run")
        self.assertEqual(sections._count_phrase(2, "execution", "none", " in this run: x"), "2 executions in this run: x")


class StrandsViaMcpTest(unittest.TestCase):
    """V7-101: under manifest v7 chapter N's background says the service's answer carries each entry's strands and
    N-usecase-0005 checks them via MCP; appendix B qualifies the pinned table's pattern.strands_vocabulary statement
    without changing it; v5 / v6 read as before."""

    def test_branches(self):
        n = [p for p in prose("N")["intro"] if "pattern.strands_vocabulary" in p][0]
        v7 = facts.fill(n, {"manifest_version": 7}, unresolved=[])
        v6 = facts.fill(n, {"manifest_version": 6}, unresolved=[])
        # W9: credited as chapter P credits it, the service and the MCP runner's normaliser
        self.assertIn("(`pattern.strands_vocabulary`). In v7 the service states each entry's `strands` per database and "
                      "the MCP runner's normaliser (git head {mcp_runner_head}) copies the host's values onto the "
                      "entry, so N-usecase-0005 checks that field through MCP as well. Appendix B", v7)
        self.assertNotIn("does not echo", v7)
        self.assertIn("(`pattern.strands_vocabulary`), and the service's answer does not echo an entry's `strands`, so "
                      "that field is not exposed there. Appendix B", v6)
        # appendix B: the short reason (the pinned table's statement) unchanged, a v7 note beside it
        self.assertTrue(prose("cases")["short"]["pattern.strands_vocabulary"].endswith(
            "does not echo an entry's `strands`"))
        note = prose("cases")["current"]["pattern.strands_vocabulary"]
        ctx7 = _Stub(run=_Stub(manifest={"manifest_version": 7}))
        ctx6 = _Stub(run=_Stub(manifest={"manifest_version": 6}))
        self.assertEqual(sections._versioned(ctx7, note), "This source-table description predates v7: entry strands "
                         "are now exposed and checked.")
        self.assertEqual(sections._versioned(ctx6, note).strip(), "")


class PreviousRunRefusalTest(unittest.TestCase):
    """W6: the previous edition's run must be earlier and of the same or a lower manifest version; a service-log time
    must be UTC."""

    def _write(self, d, run_id, version):
        os.makedirs(d)
        with open(os.path.join(d, "manifest.json"), "w") as fh:
            json.dump({"manifest_version": version, "scenarios": [{"id": "A-1"}]}, fh)
        with open(os.path.join(d, "executions.jsonl"), "w") as fh:
            fh.write(json.dumps({"scenario": "A-1", "arm": "mcp", "repeat": 0, "outcome": "pass",
                                 "started_utc": "2026-10-17T08:30:00Z"}) + "\n")
        with open(os.path.join(d, "run.json"), "w") as fh:
            json.dump({"run_id": run_id}, fh)

    def test_refusals(self):
        import shutil
        import tempfile
        import reader
        tmp = tempfile.mkdtemp()
        try:
            cur = {"run_id": "20261019T005230Z_x", "manifest_version": 7,
                   "started": reader.run_start("20261019T005230Z_x", [])}
            ok = os.path.join(tmp, "v6")
            self._write(ok, "20261017T082918Z_x", 6)
            self.assertEqual(reader.PreviousRun(ok, None, cur).run_id(), "20261017T082918Z_x")
            later = os.path.join(tmp, "later")
            self._write(later, "20261020T000000Z_x", 7)
            with self.assertRaises(reader.RunError) as cm:
                reader.PreviousRun(later, None, cur)
            self.assertIn("is not earlier than this run", str(cm.exception))
            higher = os.path.join(tmp, "v8")
            self._write(higher, "20261018T000000Z_x", 8)
            with self.assertRaises(reader.RunError) as cm:
                reader.PreviousRun(higher, None, cur)
            self.assertIn("manifest version 8 is above this run's 7", str(cm.exception))
            with self.assertRaises(reader.RunError):
                reader.PreviousRun(ok, ok, cur)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    def test_service_log_time_must_be_utc(self):
        import tempfile
        import reader
        rows = [{"started_utc": "2026-10-19T02:00:00.000Z", "ended_utc": "2026-10-19T03:00:00.000Z"}]
        stub = _Stub(rows=rows, run={"release": None}, mcp_run={}, service_identity=None, manifest={},
                     run_id=lambda: "R")
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "identity.json")
            for t, ok in (("2026-10-19T02:10:05Z", True), ("2026-10-19T02:10:05+00:00", True),
                          ("2026-10-19T02:10:05+02:00", False), ("2026-10-19T02:10:05", False)):
                with open(path, "w") as fh:
                    json.dump({"run": "R", "service_commit": "abc1234", "service_log": {"text": "x", "times": [t]}}, fh)
                if ok:
                    reader.Run._load_identity(stub, path)
                else:
                    with self.assertRaises(reader.RunError) as cm:
                        reader.Run._load_identity(stub, path)
                    self.assertIn("ending in Z", str(cm.exception))


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
