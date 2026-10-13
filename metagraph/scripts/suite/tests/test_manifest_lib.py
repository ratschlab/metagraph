"""Offline tests of the request builders, the MCP mapping and the frozen manifests."""
import json
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.dirname(HERE)
sys.path.insert(0, SUITE)

import checks  # noqa: E402
import manifest_lib as ml  # noqa: E402


class Builders(unittest.TestCase):
    def test_pattern_request_names_the_defaults(self):
        r = ml.pattern_request([{"dna": "ACGT" * 6}])
        self.assertEqual((r["mode"], r["max_contexts"], r["time_budget_ms"]), ("all_or_count", 10000, 60000))
        r = ml.pattern_request([{"dna": "A" * 25}], mode="partial", output={"labels": "all"})
        self.assertEqual(r["max_labels"], 1000)
        r = ml.pattern_request([{"dna": "A" * 25}], predicate={"any": ["562"]})
        self.assertEqual(r["max_predicate_contexts"], 10000)

    def test_defaults_named_where_they_act_on_both_arms(self):
        r = ml.pattern_request([{"dna": "A" * 25}], mode="count")
        self.assertEqual((r["scope"], r["strands"], r["output"]), ("any_offset", "both", {"labels": "none"}))
        self.assertNotIn("predicate_strands", r)            # a service refusal without a predicate
        r = ml.pattern_request([{"dna": "A" * 25}], output={"labels": "all"}, strands="reverse")
        self.assertEqual(r["output"], {"labels": "all", "occurrences": True})
        self.assertEqual(r["strands"], "reverse")
        r = ml.pattern_request([{"dna": "A" * 25}], output={"labels": "none", "occurrences": True})
        self.assertEqual(r["output"], {"labels": "none", "occurrences": True})   # a refusal case stays as asked
        r = ml.pattern_request([{"dna": "A" * 25}], output={"labels": "predicate_only"}, predicate={"any": ["562"]})
        self.assertEqual((r["output"]["occurrences"], r["predicate_strands"]), (True, "either"))
        a = ml.mcp_pattern(r)["arguments"]
        self.assertEqual((a["scope"], a["strands"], a["labels"], a["occurrences"], a["predicate_strands"]),
                         ("any_offset", "both", "predicate_only", True, "either"))
        a = ml.mcp_pattern(ml.pattern_request([{"dna": "A" * 25}], strands="forward"))["arguments"]
        self.assertEqual(a["strands"], "plus")
        self.assertNotIn("occurrences", a)

    def test_resolve_request_mirrors_the_service(self):
        r = ml.resolve_request("A" * 100, labels=["562"])
        self.assertEqual(r["bounds"], {"max_query_bp": 2000, "time_budget_ms": 20000})
        r = ml.resolve_request("A" * 100, discover={"max_labels": 5})
        self.assertEqual(r["bounds"], {"max_query_bp": 500})
        r = ml.resolve_request("A" * 3000, labels=["562"])
        self.assertNotIn("max_query_bp", r["bounds"])

    def test_traverse_request_is_the_preset(self):
        r = ml.traverse_request([{"sequence": "A" * 40}], {"direction": "left", "bounds": {"max_extension_bp": 100}})
        st = r["strategy"]
        self.assertEqual(st["direction"], "left")
        self.assertEqual(st["bounds"]["max_extension_bp"], 100)
        self.assertEqual(st["bounds"]["max_memory_mb"], 512)
        self.assertEqual(st["branching"]["max_label_branches"], 2)
        self.assertEqual(st["output"]["detail"], "full")


class McpMapping(unittest.TestCase):
    def test_pattern(self):
        req = ml.pattern_request([{"dna": "A" * 25, "id": "a"}], strands="forward", output={"labels": "all"})
        m = ml.mcp_pattern(req)
        self.assertEqual(m["tool"], "pattern_search")
        a = m["arguments"]
        self.assertEqual((a["strands"], a["labels"], a["mode"], a["max_contexts"]), ("plus", "all", "all_or_count", 10000))
        self.assertIn("pattern.strands_vocabulary", m["by_design"])
        self.assertNotIn("idempotency_key", a)
        m = ml.mcp_pattern(ml.pattern_request([{"dna": "A" * 25}], max_steps=1))
        self.assertEqual(m["not_implemented"], "max_steps")
        m = ml.mcp_pattern(dict(ml.pattern_request([{"dna": "A" * 25}]), max_mismatches="1"))
        self.assertIn("max_mismatches", m["not_implemented"])
        m = ml.mcp_pattern(ml.pattern_request([{"dna": "A" * 25}], max_mismatches=2))
        self.assertEqual((m["arguments"]["max_mismatches"], m["requires_service"]), (2, "I12"))
        m = ml.mcp_pattern(ml.pattern_request([{"dna": "A" * 40}]))
        self.assertIn("pattern.long_without_long_search", m["by_design"])

    def test_resolve_and_traverse(self):
        m = ml.mcp_resolve(ml.resolve_request("ACGT" * 10, labels=["562", "562"]))
        self.assertEqual(m["tool"], "traverse_resolve")
        self.assertIn("resolve.duplicate_labels", m["by_design"])
        self.assertEqual(ml.mcp_resolve(ml.resolve_request("A" * 40, labels=["1"], support="trace"))["not_implemented"], "support")
        m = ml.mcp_traverse(ml.traverse_request([{"sequence": "A" * 40, "labels": ["x"]}],
                                                {"direction": "right", "bounds": {"max_extension_bp": 100}}))
        self.assertEqual(m["arguments"]["preset"], "walk")
        self.assertEqual(m["arguments"]["overrides"], {"direction": "right", "bounds": {"max_extension_bp": 100}})
        self.assertTrue(m["arguments"]["force_fresh"])
        m = ml.mcp_traverse(ml.traverse_request([{"sequence": "A" * 40}], {}, extra={"labels": {"mode": "annotate"}}))
        self.assertEqual(m["not_implemented"], "strategy.labels.mode")
        m = ml.mcp_traverse(ml.traverse_request([{"sequence": "A" * 40}], {}, extra={"output": {"coordinates": True}}))
        self.assertEqual(m["not_implemented"], "strategy.output.coordinates")
        m = ml.mcp_traverse(ml.traverse_request([{"sequence": "A" * 40}], {"bounds": {"max_extension_bp": 50}}, preset="trie"))
        self.assertEqual(m["arguments"]["preset"], "trie")

    def test_budget_expectation(self):
        ok = {"check": "status", "args": {"http": 200}}
        s = ml.scenario("P-1", "P", "t", "POST /pattern", ml.pattern_request([{"dna": "A" * 25}]), {}, [ok])
        b = [e for e in s["expect"] if e["check"] == "budget_respected"]
        self.assertEqual(b[0]["args"]["allowance_ms"], 0)
        s = ml.scenario("T-1", "T", "t", "POST /traverse", ml.traverse_request([{"sequence": "A" * 40}]), {},
                        [{"check": "status", "args": {"http": [200, 503]}}])
        b = [e for e in s["expect"] if e["check"] == "budget_respected"]
        self.assertEqual(b[0]["args"]["allowance_ms"], 1000)    # the capabilities' delivery_reserve.stop_ms
        s = ml.scenario("R-1", "R", "t", "POST /resolve", ml.resolve_request("A" * 40, labels=["1"], time_budget_ms=250),
                        {}, [{"check": "status", "args": {"http": 400}}])
        self.assertFalse(any(e["check"] == "budget_respected" for e in s["expect"]))   # a refusal: no answer to time
        s = ml.scenario("R-2", "R", "t", "POST /resolve", ml.resolve_request("A" * 40, discover={"max_labels": 5}),
                        {}, [ok])
        self.assertFalse(any(e["check"] == "budget_respected" for e in s["expect"]))   # a discovery names none

    def test_scenario_levels(self):
        s = ml.scenario("X-1", "P", "t", "POST /pattern", ml.pattern_request([{"dna": "A" * 25}]), {},
                        [{"check": "status", "args": {"http": 200}}, {"check": "contains_source_context", "args": {}}])
        self.assertEqual(s["expect"][1]["args"]["level"], 0)


class FrozenManifests(unittest.TestCase):
    """The committed manifests: their digest, the contract's fields, every check known to checks.py."""

    def test_manifests(self):
        found = 0
        for name in ("suite_v1.json", "suite_smoke_v1.json"):
            path = os.path.join(SUITE, "manifest", name)
            if not os.path.exists(path):
                continue
            found += 1
            with open(path) as f:
                m = json.load(f)
            self.assertEqual(ml.sha256_json(m["scenarios"]), m["content_sha256"], name)
            ids = [s["id"] for s in m["scenarios"]]
            self.assertEqual(len(ids), len(set(ids)))
            for s in m["scenarios"]:
                for e in s["expect"]:
                    self.assertIn(e["check"], checks.CHECKS, s["id"])
                self.assertTrue("tool" in s["mcp"] or "not_implemented" in s["mcp"], s["id"])
        if found == 0:
            self.skipTest("no manifest built yet")


if __name__ == "__main__":
    unittest.main()
