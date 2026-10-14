"""Two cases the paired three-arm smoke of suite_smoke_v3 showed: a /traverse attempt that reached its duration
bound (the MCP arm's seed error) is a stated stop, so the property checks of that seed are inconclusive, not
failures (lane S4, after the verification of lane S3: only that seed's, only the host's exact wording, on the MCP
arm, and only once the stated bound can have been reached); and a by-design case may name a field with a
placeholder (bounds.<knob>)."""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import checks  # noqa: E402
import compare  # noqa: E402

# AttemptAtBound as the deployed host words it (run 1 of the paired smoke), and as the current source does
BOUND = ("the attempt reached the duration bound the server enforces for it (15000 ms: the seeds' time budgets "
         "plus the server's allowance, at most hard_cap_ms, the content timeout less one second; see "
         "usage.bound) while its response was built or written: nothing of it is delivered")
BOUND_IN_RAM = ("the attempt reached the duration bound the server enforces for it (15000 ms: the in_ram load (if "
                "any) plus the seeds' time budgets plus the server's allowance, at most hard_cap_ms, the content "
                "timeout less one second; see usage.bound) while its response was built or written: nothing of it "
                "is delivered")
# the MCP tool's own hint (search-service mcp_tools.py), which shares the first words
HINT = "the attempt reached the duration bound the server enforces for it; a retry needs a new attempt_id"
# a seed the attempt's walk-until stopped before it began (traverse.cpp not_started_seed_to_json), both wordings
CAUSE_OLD = ("the attempt reached the time at which the server stops walking it (the duration bound it enforces, the "
             "seeds' time budgets plus its allowance (at most hard_cap_ms), less the larger of half the allowance and "
             "the delivery reserve, kept for the delivery; see usage.bound.walk_until_ms)")
CAUSE_NEW = ("the attempt reached the time at which the server stops walking it (the duration bound it enforces, "
             "in_ram load + seeds' time budgets + allowance (at most hard_cap_ms), less the larger of half the "
             "allowance and the delivery reserve, kept for the delivery; see usage.bound.walk_until_ms)")


def not_started(cause, ms):
    return ("not started: %s before this seed began, %d ms after the request was received; no budget of the request "
            "ran out: a new attempt can traverse it" % (cause, ms))


def seed_phase(ms):
    return ("the attempt reached the time at which the server stops walking it (its duration bound less the larger of "
            "half the allowance and the delivery reserve; see usage.bound.walk_until_ms) while the seed was validated "
            "against its labels, %d ms after the request was received: no traversal was made" % ms)


class AttemptBound(unittest.TestCase):
    def scenario(self, seeds=1):
        expect = [{"check": "status", "args": {"http": 200}}, {"check": "invariants", "args": {}}]
        expect += [{"check": "traverse_seed_outcome", "args": {"walks": ["complete", "partial"], "seed": i}}
                   for i in range(seeds)]
        return {"route": "POST /traverse", "request": {"seeds": [{"sequence": "TTAGGG" * 10}] * seeds},
                "expect": expect}

    def test_bound_is_inconclusive(self):
        for text in (BOUND, BOUND_IN_RAM):
            ans = {"http": 200, "body": {"results": [{"error": text}]}, "error": None}
            res = checks.evaluate(self.scenario(), ans, "mcp")
            by = {r["check"]: r["outcome"] for r in res}
            self.assertEqual("pass", by["status"])
            self.assertEqual("inconclusive", by["traverse_seed_outcome"])
            # the invariants judge the other seeds (none here) and pass over the stopped one
            self.assertEqual("pass", by["invariants"])
            self.assertEqual("inconclusive", checks.scenario_outcome(res))
            self.assertIn("AttemptAtBound", checks.time_limited(ans))

    def test_other_error_still_fails(self):
        for text in ("Seed is not fully present in the graph", HINT, BOUND[:-1], BOUND + " (relayed)",
                     "the attempt reached its bound"):
            ans = {"http": 200, "body": {"results": [{"error": text}]}, "error": None}
            res = checks.evaluate(self.scenario(), ans, "mcp")
            self.assertEqual("fail", checks.scenario_outcome(res), text)

    def test_mcp_arm_only(self):
        # a direct arm receives AttemptAtBound as the request's 503: a 200 seed result stating it is not the host's
        ans = {"http": 200, "body": {"results": [{"error": BOUND}]}, "error": None}
        for arm in ("direct-local", "direct-remote"):
            res = checks.evaluate(self.scenario(), ans, arm)
            self.assertEqual("fail", checks.scenario_outcome(res), arm)

    def test_only_the_stopped_seed(self):
        # seed 0 stopped by the bound, seed 1 failed otherwise: seed 1 keeps its checks (and the invariants theirs)
        ans = {"http": 200, "body": {"results": [
            {"error": BOUND},
            {"error": "Seed is not fully present in the graph", "outcome": {"walks": "failed"}}]}, "error": None}
        res = checks.evaluate(self.scenario(2), ans, "mcp")
        seeds = [r for r in res if r["check"] == "traverse_seed_outcome"]
        self.assertEqual(["inconclusive", "fail"], [r["outcome"] for r in seeds])
        self.assertEqual("fail", checks.scenario_outcome(res))
        # without its outcome, seed 1 fails the invariants as before
        ans["body"]["results"][1] = {"error": "Seed is not fully present in the graph"}
        by = {r["check"]: r for r in checks.evaluate(self.scenario(2), ans, "mcp")}
        self.assertEqual("fail", by["invariants"]["outcome"])
        self.assertIn("result 1 without outcome", by["invariants"]["detail"])
        self.assertNotIn("result 0", by["invariants"]["detail"])

    def test_bound_must_have_been_reached(self):
        ans = {"http": 200, "body": {"results": [{"error": BOUND}]}, "error": None}
        # the client's wall time (the MCP runner passes the service's backend time) at or after the bound
        res = checks.evaluate(self.scenario(), ans, "mcp", client={"wall_ms": 15305.3})
        self.assertEqual("inconclusive", checks.scenario_outcome(res))
        res = checks.evaluate(self.scenario(), ans, "mcp", client={"wall_ms": 9000.0})
        by = {r["check"]: r for r in res}
        self.assertEqual("fail", by["traverse_seed_outcome"]["outcome"])
        self.assertIn("bound of 15000 ms stated as reached after 9000 ms", by["traverse_seed_outcome"]["detail"])
        self.assertEqual("fail", by["invariants"]["outcome"])
        # the answer's own usage.elapsed_ms
        ans["body"]["usage"] = {"elapsed_ms": 4000}
        res = checks.evaluate(self.scenario(), ans, "mcp")
        self.assertEqual("fail", checks.scenario_outcome(res))
        ans["body"]["usage"] = {"elapsed_ms": 15002}
        self.assertEqual("inconclusive", checks.scenario_outcome(checks.evaluate(self.scenario(), ans, "mcp")))

    def test_walk_until_stops(self):
        # a seed the attempt's walk-until stopped before it began, or in its seed phase: any arm, both wordings
        for text in (not_started(CAUSE_OLD, 14002), not_started(CAUSE_NEW, 14002), seed_phase(14002)):
            seed = {"error": text, "outcome": {"walks": "failed"},
                    "resource_stop": {"resource": "attempt_deadline", "used": 14002, "requested": 15000}}
            body = {"results": [seed], "usage": {"elapsed_ms": 14100, "bound": {"walk_until_ms": 14000}}}
            ans = {"http": 200, "body": body, "error": None}
            for arm in ("direct-local", "mcp"):
                res = checks.evaluate(self.scenario(), ans, arm)
                self.assertEqual("inconclusive", checks.scenario_outcome(res), (arm, text[:40]))
            self.assertIn("attempt's duration bound", checks.time_limited(ans))
            # before the walk-until the answer states: cannot have happened
            body["usage"]["bound"]["walk_until_ms"] = 20000
            res = checks.evaluate(self.scenario(), ans, "direct-local")
            self.assertEqual("fail", checks.scenario_outcome(res))
            body["usage"]["bound"]["walk_until_ms"] = 14000
            # a resource_stop of another resource
            seed["resource_stop"]["resource"] = "time"
            self.assertEqual("fail", checks.scenario_outcome(checks.evaluate(self.scenario(), ans, "direct-local")))
        # a paraphrase of the walk-until stop is no statement of the host
        ans = {"http": 200, "body": {"results": [{"error": "not started: the attempt ran out of time",
                                                  "outcome": {"walks": "failed"}}]}, "error": None}
        self.assertEqual("fail", checks.scenario_outcome(checks.evaluate(self.scenario(), ans, "direct-local")))


SRC = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "src"))


def literal_groups(path):
    """Every run of adjacent C++ string literals in |path| (separated by whitespace or // comments), joined."""
    with open(path, encoding="utf-8") as fh:
        text = fh.read()
    groups, cur, i = [], None, 0
    while i < len(text):
        c = text[i]
        if c == "'" and not text[i - 1].isalnum():   # a character literal ('"', '\''), not a digit separator
            j = i + 1
            while text[j] != "'":
                j += 2 if text[j] == "\\" else 1
            if cur is not None:
                groups.append(cur)
            cur = None
            i = j + 1
        elif c == '"':
            j, buf = i + 1, []
            while text[j] != '"':
                if text[j] == "\\":
                    buf.append({"n": "\n", "t": "\t"}.get(text[j + 1], text[j + 1]))
                    j += 2
                else:
                    buf.append(text[j])
                    j += 1
            cur = (cur or "") + "".join(buf)
            i = j + 1
        elif c.isspace():
            i += 1
        elif text.startswith("//", i):
            i = text.index("\n", i)
        else:
            if cur is not None:
                groups.append(cur)
            cur = None
            i += 1
    return groups


@unittest.skipUnless(os.path.isdir(SRC), "the C++ source is not beside the suite")
class HostWording(unittest.TestCase):
    """The texts the rule matches are the host's own: rebuilt from the string literals of the current source (the
    deployed host's wording, without the in_ram term, is the test constants above)."""

    def test_attempt_at_bound(self):
        fmt = [g for g in literal_groups(os.path.join(SRC, "cli", "traverse_attempts.cpp"))
               if g.startswith("the attempt reached the duration bound the server enforces for it")]
        self.assertEqual(1, len(fmt))
        text = fmt[0].replace("{}", "32798")
        self.assertTrue(checks._AT_BOUND_RE.fullmatch(text), text)
        self.assertEqual(BOUND_IN_RAM.replace("15000", "32798"), text)
        self.assertFalse(checks._AT_BOUND_RE.fullmatch(HINT))

    def test_walk_until(self):
        groups = literal_groups(os.path.join(SRC, "cli", "traverse.cpp"))
        cause = [g for g in groups if g.startswith("the attempt reached the time at which the server stops")]
        self.assertEqual([CAUSE_NEW], cause)
        self.assertIn(" before this seed began, ", groups)
        self.assertIn(" ms after the request was received; no budget of the request ran out: a new attempt can "
                      "traverse it", groups)
        self.assertTrue(checks._NOT_STARTED_RE.fullmatch(not_started(cause[0], 7)))
        walker = literal_groups(os.path.join(SRC, "graph", "traversal", "walker.cpp"))
        lead = [g for g in walker if g.startswith("the attempt reached the time at which the server stops")]
        self.assertEqual(1, len(lead))
        for part in (" while the seed was ", "read to derive its permitted set", "validated against its labels",
                     " ms after the request was received: no traversal was made"):
            self.assertIn(part, walker)
        for kind in ("read to derive its permitted set", "validated against its labels"):
            text = "%s while the seed was %s, 12 ms after the request was received: no traversal was made" % (
                lead[0], kind)
            self.assertTrue(checks._SEED_PHASE_STOP_RE.fullmatch(text), text)


class FieldPlaceholder(unittest.TestCase):
    def test_placeholder(self):
        self.assertTrue(compare._field_matches("bounds.time_budget_ms", "bounds.<knob>"))
        self.assertFalse(compare._field_matches("strategy.time_budget_ms", "bounds.<knob>"))
        # a placeholder matches anything, dots included, as async's runner (_like) reads it (lane S5's B2, G07):
        # "predicate<path of the first fault>" names the path of the first fault
        self.assertTrue(compare._field_matches("bounds.a.b", "bounds.<knob>"))
        self.assertTrue(compare._field_matches("predicate.at_least.n", "predicate<path of the first fault>"))
        self.assertTrue(compare._field_matches("predicate", "predicate<path of the first fault>"))
        self.assertFalse(compare._field_matches("max_mismatches", "predicate<path of the first fault>"))
        self.assertTrue(compare._field_matches("max_mismatches", "max_mismatches"))
        self.assertFalse(compare._field_matches("max_contexts", "max_mismatches"))
        # alternatives "a | b", split as codes are (G08)
        self.assertTrue(compare._field_matches("predicate_scope", "predicate_strands | predicate_scope"))
        self.assertTrue(compare._field_matches("predicate_strands", "predicate_strands | predicate_scope"))
        self.assertFalse(compare._field_matches("predicate", "predicate_strands | predicate_scope"))
        self.assertTrue(compare._field_matches("caps.max_paths", "caps.<cap> | max_mismatches"))
        self.assertFalse(compare._field_matches(None, "max_mismatches"))


if __name__ == "__main__":
    unittest.main()
