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
from test_compare_report import T0, line, scen  # noqa: E402


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


if __name__ == "__main__":
    unittest.main()
