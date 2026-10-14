"""Offline tests of the verdict rules of compare.py, reevaluate.py, report.py and suite_common.py that the two-day
review's verdict findings name: the direct arms' class with equal digests, the MCP reference per repeat and check,
the host cache after a relayed deadline, the text_contains grammar, the no-table rule (c), the runner's unmatched
not_implemented_mcp mark, outcomes outside the vocabulary, what an area tested and the load caveat."""

import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import compare as cmp  # noqa: E402
import suite_common as sc  # noqa: E402
from test_compare_report import T0, line, scen, write_arm  # noqa: E402


def st(outcome, name="status", detail=""):
    return {"check": name, "outcome": outcome, "detail": detail}


class DirectIdenticalTest(unittest.TestCase):
    """Equal digests never hide a differing check or an outcome difference no timing explains."""

    def pair(self, remote_budget):
        deadline = {"code": "deadline", "text": "nothing written within the time budget"}
        local = line("R-1", "direct-local", 0, T0, 9000.0, outcome="inconclusive", http=503, digest="same",
                     error=deadline, checks=[st("inconclusive"), st("inconclusive", "budget_respected")])
        remote = line("R-1", "direct-remote", 0, T0 + 20, 2000.0, http=503, digest="same", error=deadline,
                      outcome="fail" if remote_budget == "fail" else "inconclusive",
                      checks=[st("inconclusive"), st(remote_budget, "budget_respected")])
        return {0: local}, {0: remote}

    def test_identical_bodies_one_premature_disagree(self):
        local, remote = self.pair("fail")
        local[0]["checks"][1] = st("pass", "budget_respected")
        local[0]["outcome"] = "pass"
        res = cmp.direct_agreement(local, remote)
        self.assertEqual(res["class"], "disagree")
        self.assertTrue(res["digest_equal"])
        self.assertEqual([d["check"] for d in res["differing_checks"]], ["budget_respected"])

    def test_identical_bodies_outcomes_differ_untimed_disagree(self):
        # as v4's R-budget-0002: the same 503 body, one arm's deadline late (a finding), the other's in time
        local, remote = self.pair("finding")
        remote[0]["outcome"] = "finding"
        res = cmp.direct_agreement(local, remote)
        self.assertEqual(res["class"], "disagree")
        self.assertEqual(res["outcome_differs"][0]["repeat"], 0)

    def test_equal_control_stays_identical(self):
        local, remote = self.pair("inconclusive")
        self.assertEqual(cmp.direct_agreement(local, remote)["class"], "identical")


class ReferenceFallbackTest(unittest.TestCase):
    """The MCP arm is compared per repeat and check item with direct-local, an item it left undecided taken from
    direct-remote where that one decided it."""

    INC = [st("pass", detail="503 allowed"), st("inconclusive", "invariants"), st("inconclusive", "count_relation")]
    OK = [st("pass"), st("pass", "invariants"), st("pass", "count_relation")]

    def direct(self):
        dl = {0: line("S", "direct-local", 0, T0, 10.0, outcome="inconclusive", http=503, checks=list(self.INC),
                      error={"code": "deadline", "text": "deadline"})}
        dr = {0: line("S", "direct-remote", 0, T0 + 5, 10.0, checks=list(self.OK))}
        return {"direct-local": dl, "direct-remote": dr}

    def test_mcp_failure_the_other_arm_decided_disagrees(self):
        bad = [st("pass"), st("fail", "invariants", "x"), st("pass", "count_relation")]
        m = {0: line("S", "mcp", 0, T0 + 10, 10.0, outcome="fail", checks=bad)}
        ma = cmp.mcp_agreement(m, self.direct(), scen("S", "P"), {}, None)
        self.assertEqual(ma["class"], "disagree")
        self.assertEqual([d["check"] for d in ma["differing_checks"]], ["invariants"])
        self.assertEqual(ma["reference"], "direct-local")
        self.assertEqual(ma["reference_alternative"], "direct-remote")
        self.assertEqual(sorted(f["check"] for f in ma["reference_fallbacks"]), ["count_relation", "invariants"])

    def test_n_usecase_0008_shape_agrees(self):
        # v4's N-usecase-0008: direct-local r0 a 503 (every check inconclusive), direct-remote r0 200 passing, two
        # passing MCP repeats with budget_respected not exposed
        inc = [st("inconclusive"), st("inconclusive", "invariants"), st("inconclusive", "count_relation"),
               st("inconclusive", "budget_respected")]
        ok = [st("pass"), st("pass", "invariants"), st("pass", "count_relation"), st("pass", "budget_respected")]
        mine = [st("pass"), st("pass", "invariants"), st("pass", "count_relation"),
                st("not_exposed", "budget_respected")]
        direct = {"direct-local": {0: line("N", "direct-local", 0, T0, 9000.0, outcome="inconclusive", http=503,
                                           checks=inc, error={"code": "deadline", "text": "d"})},
                  "direct-remote": {0: line("N", "direct-remote", 0, T0 + 20, 900.0, checks=ok)}}
        m = dict((r, line("N", "mcp", r, T0 + 40 + r, 10.0, checks=list(mine))) for r in (0, 1))
        ma = cmp.mcp_agreement(m, direct, scen("N", "N"), {}, None)
        self.assertEqual(ma["class"], "agree")
        self.assertGreaterEqual(ma["checks_compared"], 3)
        self.assertEqual(ma["reference_outcome"], "pass")   # direct-local's 503 gave no compared item
        self.assertEqual(ma["checks_not_exposed"], 2)

    def test_reference_outcome_keeps_a_line_that_gave_items(self):
        # direct-local decides status, falls back on the rest: both lines are the reference
        m = {0: line("S", "mcp", 0, T0 + 10, 10.0, checks=list(self.OK))}
        ma = cmp.mcp_agreement(m, self.direct(), scen("S", "P"), {}, None)
        self.assertEqual(ma["class"], "agree")
        self.assertEqual(ma["reference_outcome"], "inconclusive")

    def test_no_fallback_without_the_other_arm_or_when_decided(self):
        m = {0: line("S", "mcp", 0, T0 + 10, 10.0, checks=list(self.OK))}
        only_local = {"direct-local": self.direct()["direct-local"]}
        ma = cmp.mcp_agreement(m, only_local, scen("S", "P"), {}, None)
        self.assertEqual((ma["class"], ma["checks_compared"]), ("agree", 1))
        self.assertNotIn("reference_fallbacks", ma)
        decided = self.direct()
        decided["direct-local"][0]["checks"] = list(self.OK)
        ma = cmp.mcp_agreement(m, decided, scen("S", "P"), {}, None)
        self.assertNotIn("reference_fallbacks", ma)
        self.assertEqual(ma["checks_compared"], 3)

    def test_mcp_item_undecided_takes_no_fallback(self):
        mine = [st("pass"), st("not_exposed", "invariants"), st("inconclusive", "count_relation")]
        self.assertEqual(cmp.reference_checks(mine, self.INC, self.OK), (self.INC, []))

    def test_report_depth_reads_the_fallback(self):
        import report
        direct = self.direct()
        m = {0: line("S", "mcp", 0, T0 + 10, 10.0, checks=list(self.OK))}
        rows = [x for arm in direct.values() for x in arm.values()] + list(m.values())
        ma = cmp.mcp_agreement(m, direct, scen("S", "P"), {}, None)
        depth = report.agreement_depth({"scenarios": {"S": {"mcp_agreement": ma}}}, rows,
                                       lambda n: n not in ("status", "budget_respected"))
        self.assertEqual(depth["S"]["not_compared"], [])


class RelayedDeadlineHostCacheTest(unittest.TestCase):
    """An MCP answer relaying a host deadline reached the host: the next execution is warm, not cold."""

    def rows(self, error, time_limited=""):
        return [
            {"arm": "mcp", "scenario": "S", "repeat": 0, "http": 502, "error": error,
             "mcp_phases": {"backend_ms": None}, "server": {"elapsed_ms": None},
             "started_utc": "2026-10-14T07:12:57.757Z", "ended_utc": "2026-10-14T07:14:01.239Z",
             "time_limited": time_limited},
            {"arm": "direct-local", "scenario": "S", "repeat": 0, "http": 200, "server": {"elapsed_ms": 500},
             "wall_ms": 600, "started_utc": "2026-10-14T07:14:02.739Z", "ended_utc": "2026-10-14T07:14:03.300Z"},
            {"arm": "direct-remote", "scenario": "S", "repeat": 0, "http": 200, "server": {"elapsed_ms": 200},
             "wall_ms": 400, "started_utc": "2026-10-14T07:14:21.000Z", "ended_utc": "2026-10-14T07:14:21.400Z"}]

    def test_relayed_pattern_deadline(self):
        rows = sc.assign_host_cache(self.rows({"code": "search_failed", "text": "pattern_deadline: pattern: ..."},
                                              "the host's /pattern 503 deadline relayed by the service"))
        self.assertEqual([(r["arm"], r["host_cache"], r["since_prev_s"]) for r in rows],
                         [("mcp", "cold", None), ("direct-local", "warm", 1.5), ("direct-remote", "warm", 17.7)])
        self.assertTrue(rows[1]["host_cache_after_interrupted"])
        self.assertNotIn("host_cache_after_interrupted", rows[2])
        fo = cmp.front_overhead({0: rows[1]}, {0: rows[2]})
        self.assertIsNotNone(fo)
        self.assertEqual(fo["after_interrupted"], ["direct-local"])

    def test_each_relay_code(self):
        for err in ({"code": "traversal_resolve_deadline", "text": "the backend could not build and write it"},
                    {"code": "resolve_timeout", "text": "the backend did not answer within 30.25 s"},
                    {"code": "uninterruptible_overrun_exceeds_allowance", "text": "x"}):
            rows = sc.assign_host_cache(self.rows(err))
            self.assertEqual((rows[0]["host_cache"], rows[1]["host_cache"]), ("cold", "warm"), err["code"])

    def test_refusal_at_submit_still_not_reached(self):
        for err in ({"code": "search_failed", "text": "backend unavailable"},
                    {"code": "pattern_invalid", "text": "pattern_deadline is not a field"}):
            rows = sc.assign_host_cache(self.rows(err))
            self.assertEqual((rows[0]["host_cache"], rows[1]["host_cache"]), ("not_reached", "cold"), err)
            self.assertNotIn("host_cache_after_interrupted", rows[1])

    def test_mcp_overhead_and_service_answer_need_a_backend_timing(self):
        relay = self.rows({"code": "traversal_resolve_deadline", "text": "t"})[0]
        self.assertTrue(sc.reached_host(relay))
        self.assertFalse(sc.backend_timed(relay))
        self.assertIsNone(cmp.mcp_overhead({0: relay}, {}, {}))
        self.assertFalse(cmp._service_answered(relay))

    def test_direct_deadline_marks_the_next_execution(self):
        rows = self.rows({"code": "search_failed", "text": "pattern_deadline: x"})[1:]
        rows[0].update(http=503, error={"code": "deadline", "text": "nothing written"})
        sc.assign_host_cache(rows)
        self.assertEqual(rows[1]["host_cache"], "warm")
        self.assertTrue(rows[1]["host_cache_after_interrupted"])


class TextContainsGrammarTest(unittest.TestCase):
    """An answer_match side's text_contains: a list on the direct side (every fragment within one text), one string
    on the MCP side as the runner reads it; other types refused at load."""

    TEXTS = ["request.seeds: expected a non-empty array"]

    def test_direct_list_matches_and_does_not_without_raising(self):
        side = {"http": 400, "text_contains": ["request.seeds", "non-empty"]}
        self.assertTrue(cmp.answer_side_matches(side, 400, set(), self.TEXTS, slot_codes=set(),
                                                known=cmp.ANSWER_DIRECT_KEYS))
        side = {"http": 400, "text_contains": ["request.seeds", "at least one"]}
        self.assertFalse(cmp.answer_side_matches(side, 400, set(), self.TEXTS, slot_codes=set(),
                                                 known=cmp.ANSWER_DIRECT_KEYS))
        # the fragments must be within ONE text
        self.assertFalse(cmp.answer_text_contains(["seeds", "slot"], ["seeds here", "slot there"]))

    def test_mcp_list_is_one_string(self):
        side = {"code": "invalid_arguments", "text_contains": ["request.seeds"]}
        self.assertFalse(cmp.answer_side_matches(side, 400, {"invalid_arguments"}, self.TEXTS,
                                                 known=cmp.ANSWER_MCP_KEYS, http_list=False))
        self.assertTrue(cmp.answer_text_contains(["a"], ["x ['a'] y"], mcp_side=True))

    def test_table_names_a_bad_type(self):
        case = {"id": "x.y", "agreement": "by_design", "direct_expectation": {"http": 400},
                "service_answer": {"outcome": "refused_at_submit", "code": "invalid_arguments"},
                "answer_match": {"routes": ["POST /traverse"], "wall_rule": "budget_minus_reserve",
                                 "direct": {"http": 400, "text_contains": 7, "code": ["a", "b"], "slot_code": [""]},
                                 "mcp": {"code": "invalid_arguments", "text_contains": ["provide"], "field": 3}}}
        probs = cmp.table_problems([case])
        self.assertTrue(any("answer_match.direct.text_contains must be" in p for p in probs), probs)
        self.assertTrue(any("answer_match.direct.slot_code must be" in p for p in probs), probs)
        self.assertFalse(any("answer_match.direct.code" in p for p in probs), probs)
        self.assertTrue(any("answer_match.mcp.text_contains must be one non-empty string" in p for p in probs), probs)
        self.assertTrue(any("answer_match.mcp.field must be" in p for p in probs), probs)

    def test_shipped_table_loads(self):
        # the shipped table passes the load validation (TableError otherwise), with v4's scenarios too
        manifest = sc.read_json(os.path.join(os.path.dirname(HERE), "manifest", "suite_v4.json"))
        cases, info = cmp.load_by_design(cmp.DEFAULT_BY_DESIGN, manifest["scenarios"])
        self.assertEqual(info["cases"], len(cases))
        self.assertTrue(cases)


class RuleCWithoutTableTest(unittest.TestCase):
    """With no table, a runner's rule (c) inconclusive on its by-design case is kept: inconclusive in the
    re-evaluation, neither for nor against the direct arm in the comparison."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = self.tmp.name
        os.makedirs(os.path.join(self.dir, "raw"))
        self.s = {"id": "M-1", "route": "POST /pattern", "request": {"patterns": [], "time_budget_ms": 60000},
                  "mcp": {"tool": "pattern_search", "arguments": {"time_budget_ms": 60000}},
                  "expect": [{"check": "status", "args": {"http": 200}}]}
        self.ok = [st("pass"), st("pass", "invariants"), st("pass", "count_relation")]
        self.bad = [st("fail"), st("fail", "invariants"), st("fail", "count_relation")]
        self.err = {"code": "search_failed",
                    "text": "pattern_deadline: pattern: the answer could not be written within time_budget_ms (60000 ms)"}

    def tearDown(self):
        self.tmp.cleanup()

    def raw_line(self, arm, repeat, http=200, error=None, outcome="pass", checks=None, **extra):
        ln = line("M-1", arm, repeat, T0, 10.0, outcome=outcome, http=http, error=error, checks=checks)
        ln.update(extra)
        sc.write_gz_json(os.path.join(self.dir, ln["raw"]), {"answer": {"http": http, "body": None, "error": error}})
        return ln

    def rule_c(self, **extra):
        return self.raw_line("mcp", 0, http=502, error=self.err, outcome="inconclusive", checks=list(self.bad),
                             inconclusive_reason=sc.UNDECIDABLE_NO_WALL,
                             service_by_design_case="pattern.deadline_ends_task", explained_checks=[],
                             verdict_detail=sc.UNDECIDABLE_NO_WALL, **extra)

    def test_compare_without_table_is_not_comparable(self):
        d = {"direct-local": {0: self.raw_line("direct-local", 0, checks=list(self.ok))}}
        res = cmp.mcp_agreement({0: self.rule_c()}, d, self.s, {}, self.dir)
        self.assertEqual((res["class"], res["by"]), ("not_comparable", "undecidable plausibility"))
        cases, _ = cmp.load_by_design(cmp.DEFAULT_BY_DESIGN)
        self.assertEqual(cmp.mcp_agreement({0: self.rule_c()}, d, self.s, cases, self.dir)["class"], "not_comparable")
        # another repeat failing its checks still decides
        broken = self.raw_line("mcp", 1, outcome="fail", checks=[st("pass"), st("fail", "invariants", "x"),
                                                                   st("pass", "count_relation")])
        res = cmp.mcp_agreement({0: self.rule_c(), 1: broken}, d, self.s, {}, self.dir)
        self.assertEqual(res["class"], "disagree")

    def test_reevaluate_without_table_keeps_it(self):
        import reevaluate

        class Checks(object):
            scenario_outcome = staticmethod(sc.outcome_from_checks)
        info = {}
        ans = {"http": 502, "error": self.err, "body": None}
        out, kept = reevaluate.outcome_of(Checks(), "mcp", self.s, ans, self.rule_c(), self.bad, False, {}, None,
                                          info, None)
        self.assertEqual(out, "inconclusive")
        self.assertIn("rule (c)", kept)
        self.assertEqual(info["inconclusive_reason"], sc.UNDECIDABLE_NO_WALL)
        self.assertEqual(info["service_by_design_case"], "pattern.deadline_ends_task")
        # an inconclusive line naming no case is judged by its checks as before
        plain = dict(self.rule_c(), service_by_design_case=None)
        self.assertEqual(reevaluate.outcome_of(Checks(), "mcp", self.s, ans, plain, self.bad, False, {}, None, {},
                                               None)[0], "fail")


class UnmatchedNotImplementedMarkTest(unittest.TestCase):
    """A runner's not_implemented_mcp mark on a scenario the manifest maps: with a table, unless a case reproduces it,
    the MCP side disagrees; with no table it stays the runner's class."""

    def test_mark_on_a_mapped_scenario(self):
        s = scen("P-8", "P")
        d = {"direct-local": {0: line("P-8", "direct-local", 0, T0, 10.0)}}
        m = {0: line("P-8", "mcp", 0, T0, 10.0, outcome="not_implemented_mcp", checks=[])}
        cases, _ = cmp.load_by_design(cmp.DEFAULT_BY_DESIGN)
        res = cmp.mcp_agreement(m, d, s, cases, None)
        self.assertEqual((res["class"], res["side"], res["by"]), ("disagree", "mcp", "mcp runner"))
        self.assertIn("no case of the table reproduces the mark", res["mcp_detail"])
        self.assertEqual(cmp.mcp_agreement(m, d, s, {}, None)["class"], "not_implemented_mcp")
        # the manifest's mark, a named not_implemented_mcp case and an unmapped scenario are the plan's
        nim = scen("P-9", "P", not_implemented="max_mismatches")
        self.assertEqual(cmp.mcp_agreement(m, d, nim, cases, None)["class"], "not_implemented_mcp")
        named = scen("P-10", "P", by_design="resolve.not_expressible.bounds_max_query_bp")
        self.assertEqual(cmp.mcp_agreement(m, d, named, cases, None)["class"], "not_implemented_mcp")
        unmapped = dict(s)
        unmapped.pop("mcp")
        self.assertEqual(cmp.mcp_agreement(m, d, unmapped, cases, None)["class"], "not_implemented_mcp")

    def test_other_repeats_still_compared(self):
        s = scen("P-8", "P")
        d = {"direct-local": {0: line("P-8", "direct-local", 0, T0, 10.0)}}
        ok = line("P-8", "mcp", 1, T0 + 5, 10.0)
        m = {0: line("P-8", "mcp", 0, T0, 10.0, outcome="not_implemented_mcp", checks=[]), 1: ok}
        cases, _ = cmp.load_by_design(cmp.DEFAULT_BY_DESIGN)
        res = cmp.mcp_agreement(m, d, s, cases, None)
        self.assertEqual(res["class"], "disagree")
        self.assertEqual(res["repeats"], [0])


def run_dir_with(d, scenarios, rows_by_arm, orchestrator=None):
    """A merged and compared run directory: |rows_by_arm| {arm: [lines]}, the manifest of |scenarios|."""
    import merge_run
    manifest = {"manifest_version": 1, "content_sha256": sc.scenarios_sha256(scenarios), "scenarios": scenarios}
    for arm, rows in rows_by_arm.items():
        os.makedirs(os.path.join(d, "arms", arm), exist_ok=True)
        sc.write_json_atomic(os.path.join(d, "arms", arm, "manifest.json"), manifest)
        write_arm(d, arm, rows)
    if orchestrator is not None:
        sc.write_json_atomic(os.path.join(d, "orchestrator.json"), orchestrator)
    merge_run.merge(d)
    sc.write_json_atomic(os.path.join(d, "compare.json"), cmp.compare(d, cmp.load_by_design(None)))
    return d


def verdict_of(d):
    import contextlib
    import io
    import report
    with contextlib.redirect_stdout(io.StringIO()):
        text = report.build_report(d)
    return text, text.split("# Verdict", 1)[1].strip().splitlines()[0]


class OutcomeVocabularyTest(unittest.TestCase):
    """An outcome outside the shared vocabulary is counted and read as a failure everywhere, never as a pass."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.d = self.tmp.name
        self.S = [scen("P-001", "P"), scen("P-002", "P"), scen("P-003", "P")]

    def tearDown(self):
        self.tmp.cleanup()

    def mcp_rows(self, outcome):
        return [line(x["id"], "mcp", r, T0 + i * 10 + 5 * r, 40.0, outcome=outcome)
                for i, x in enumerate(self.S) for r in (0, 1)]

    def test_paired_run_fails(self):
        orch = {"sessions": [{"mode": "paired", "mcp_cmd": "x {ids} {out}", "checks": "checks", "warm_every": "0",
                              "loop": {"mode": "paired", "state": "done", "scenarios": 3, "completed": 3}}]}
        direct = dict((arm, [line(x["id"], arm, 0, T0 + i * 10 + k, 20.0) for i, x in enumerate(self.S)])
                      for k, arm in enumerate(("direct-local", "direct-remote")))
        run_dir_with(self.d, self.S, dict(direct, mcp=self.mcp_rows("skipped_by_runner")), orch)
        comp = sc.read_json(os.path.join(self.d, "compare.json"))
        self.assertEqual(comp["summary"]["totals"]["mcp"]["outcomes"]["other"], 6)
        self.assertEqual(len(comp["summary"]["outcomes_outside_vocabulary"]), 6)
        self.assertEqual(comp["scenarios"]["P-001"]["mcp_agreement"]["class"], "disagree")
        self.assertEqual(comp["scenarios"]["P-001"]["mcp_agreement"]["by"], "outcome")
        self.assertEqual(comp["summary"]["per_area"]["P"]["mcp"]["other"], 3)
        text, v = verdict_of(self.d)
        self.assertTrue(v.startswith("**FAIL"), v)
        self.assertIn("6 with an outcome outside the vocabulary", v)
        self.assertIn("read as failures", text)

    def test_mcp_only_run_fails(self):
        orch = {"sessions": [{"mode": "paired", "mcp_cmd": "x {ids} {out}", "checks": "checks", "warm_every": "0",
                              "skip_local": "1", "skip_remote": "1"}]}
        rows = self.mcp_rows("skipped_by_runner")
        rows[0]["checks"] = [st("fail", detail="http 500")]
        run_dir_with(self.d, self.S, {"mcp": rows}, orch)
        _, v = verdict_of(self.d)
        self.assertTrue(v.startswith("**FAIL"), v)

    def test_a_known_mark_beside_it(self):
        # one unknown repeat beside a passing one: the aggregate is the unknown outcome, the class a disagreement
        self.assertEqual(cmp.aggregate_outcome(["pass", "skipped_by_runner"]), "skipped_by_runner")
        self.assertEqual(cmp.aggregate_outcome(["pass", "error"]), "error")
        self.assertTrue(sc.outcome_known("service_by_design"))
        self.assertFalse(sc.outcome_known("skipped_by_runner"))
        self.assertFalse(sc.outcome_known(None))

    def test_direct_arm_outcome_outside_disagrees(self):
        a = {0: line("R-1", "direct-local", 0, T0, 10.0, outcome="odd")}
        b = {0: line("R-1", "direct-remote", 0, T0, 10.0, outcome="odd")}
        res = cmp.direct_agreement(a, b)
        self.assertEqual(res["class"], "disagree")
        self.assertIn("outside the vocabulary", res["outcome_differs"][0]["detail"])

    def test_reevaluate_keeps_it(self):
        import reevaluate
        info = {}
        out, kept = reevaluate.outcome_of(None, "mcp", scen("P-001", "P"), {"http": 200, "body": {}, "error": None},
                                          {"outcome": "skipped_by_runner"}, [st("pass")], False, {}, None, info)
        self.assertEqual(out, "skipped_by_runner")
        self.assertIn("outside the vocabulary", kept)
        self.assertTrue(any("outside the vocabulary" in w for w in info["warnings"]))


class LoadCaveatTest(unittest.TestCase):
    """With no sra_all figure (failed samples), the report never asserts that mex was not under the Logan load."""

    def figures(self, sra, load1, error=None, n=100):
        import report
        out = []
        for i in range(n):
            x = {"t_utc": sc.iso_utc(T0 + i), "sra_all_cpu_pct": sra, "host_load1": load1}
            if error:
                x["error"] = error
            out.append(x)
        return report, report.load_figures(out)

    def test_error_only_samples(self):
        for load1 in (7.5, None):
            report, lf = self.figures(None, load1, "Error response from daemon: No such container: sra_all")
            text = report.load_caveat(lf)
            self.assertNotIn("not under the Logan load", text)
            self.assertIn("**Load not measured:**", text)
            self.assertIn("100 of 100 samples failed", text)
            self.assertIn("No such container", report.resources_section(
                {"summary": {"resources": {}}, "arms": []}, lf))

    def test_valid_low_load_still_says_so(self):
        report, lf = self.figures(50.0, 7.5)
        self.assertIn("Load: mex was not under the Logan load", report.load_caveat(lf))
        self.assertNotIn("Failed samples", report.resources_section({"summary": {"resources": {}}, "arms": []}, lf))

    def test_many_failed_samples_beside_valid_ones(self):
        import report
        good = self.figures(50.0, 7.5, n=80)[1]
        samples = [{"t_utc": sc.iso_utc(T0 + i), "sra_all_cpu_pct": 50.0, "host_load1": 7.5} for i in range(80)] + \
            [{"t_utc": sc.iso_utc(T0 + 100 + i), "sra_all_cpu_pct": None, "host_load1": None, "error": "timeout"}
             for i in range(20)]
        self.assertIn("Load: mex was not under", report.load_caveat(good))
        text = report.load_caveat(report.load_figures(samples))
        self.assertIn("**Load not measured:** more than 10 % of the samples failed", text)

    def test_busy_stays_the_caveat(self):
        report, lf = self.figures(7000.0, 220.0)
        self.assertIn("**Load caveat:**", report.load_caveat(lf))


class UntestedAndExposureTest(unittest.TestCase):
    """A field the MCP drops: with mcp.exposes a disagreement; an area or route x arm where nothing passed and
    something was inconclusive keeps the verdict from PASS."""

    NE = [st("pass"), st("not_exposed", "contains_source_context", "field not in the answer: patterns")]

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.d = self.tmp.name
        self.S = [scen("P-001", "P"), scen("P-002", "P"), scen("P-003", "P")]
        self.orch = {"sessions": [{"mode": "paired", "mcp_cmd": "x {ids} {out}", "checks": "checks",
                                   "warm_every": "0",
                                   "loop": {"mode": "paired", "state": "done", "scenarios": 3, "completed": 3}}]}

    def tearDown(self):
        self.tmp.cleanup()

    def rows(self, mcp_checks):
        direct = dict((arm, [line(x["id"], arm, 0, T0 + i * 10 + k, 20.0) for i, x in enumerate(self.S)])
                      for k, arm in enumerate(("direct-local", "direct-remote")))
        direct["mcp"] = [line(x["id"], "mcp", r, T0 + i * 10 + 5 * r, 40.0, checks=list(mcp_checks))
                         for i, x in enumerate(self.S) for r in (0, 1)]
        return direct

    def test_dropped_field_is_not_pass(self):
        # the review's case (c): every MCP property check not_exposed, the runner says pass
        run_dir_with(self.d, self.S, self.rows(self.NE), self.orch)
        text, v = verdict_of(self.d)
        self.assertTrue(v.startswith("**INCONCLUSIVE (no execution passed in"), v)
        self.assertIn("area P on mcp (6 inconclusive)", v)
        self.assertIn("route POST /pattern on mcp", v)

    def test_dropped_field_listed_as_exposed_fails(self):
        for x in self.S:
            x["mcp"]["exposes"] = ["status", "contains_source_context"]
        run_dir_with(self.d, self.S, self.rows(self.NE), self.orch)
        comp = sc.read_json(os.path.join(self.d, "compare.json"))
        ma = comp["scenarios"]["P-001"]["mcp_agreement"]
        self.assertEqual((ma["class"], ma["side"], ma["by"], ma["checks"]),
                         ("disagree", "mcp", "exposes", ["contains_source_context"]))
        _, v = verdict_of(self.d)
        self.assertTrue(v.startswith("**FAIL"), v)
        self.assertIn("3 MCP disagreement(s)", v)

    def test_exposure_of_a_baseline(self):
        rows = self.rows([st("pass"), st("pass", "contains_source_context")])["mcp"]
        rows[0]["checks"][1] = st("not_exposed", "contains_source_context")
        self.assertEqual(cmp.mcp_exposure(rows)["P-001"], ["contains_source_context", "status"])
        self.assertIsNone(cmp.exposure_drops(self.S[0], rows))   # no mcp.exposes: nothing to judge

    def test_clean_run_and_marks_still_pass(self):
        run_dir_with(self.d, self.S, self.rows([st("pass"), st("pass", "contains_source_context")]), self.orch)
        _, v = verdict_of(self.d)
        self.assertTrue(v.startswith("**PASS"), v)

    def test_one_untested_area_beside_a_tested_one(self):
        import report
        rows = [line("P-001", "mcp", 0, T0, 10.0), line("C-001", "mcp", 0, T0, 10.0, outcome="inconclusive"),
                line("X-001", "mcp", 0, T0, 10.0, outcome="not_implemented_mcp", checks=[])]
        manifest = {"scenarios": [{"id": "P-001", "route": "POST /pattern"}, {"id": "C-001", "route": "GET /c"},
                                  {"id": "X-001", "route": "GET /x"}]}
        groups = report.untested_groups(rows, manifest)
        self.assertEqual([(g[0], g[1], g[2]) for g in groups], [("area", "C", "mcp"), ("route", "GET /c", "mcp")])


if __name__ == "__main__":
    unittest.main()
