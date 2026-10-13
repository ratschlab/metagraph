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
        self.assertIn("§7.6", b[0]["args"]["rule"])
        s = ml.scenario("T-1", "T", "t", "POST /traverse", ml.traverse_request([{"sequence": "A" * 40}]), {},
                        [{"check": "status", "args": {"http": [200, 503]}}])
        b = [e for e in s["expect"] if e["check"] == "budget_respected"]
        self.assertIn("never read at depth 0", b[0]["args"]["rule"])     # no fixed allowance: walker.cpp's pieces
        self.assertIn("first depth-1 checkpoint", b[0]["args"]["rule"])
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


class WalkExpectations(unittest.TestCase):
    """The exact walk minimum and branch depth from a record (manifest_lib.walk_expectations)."""

    def test_copies_ends_and_reentry(self):
        import random
        rng = random.Random(1)
        seq = "".join(rng.choice("ACGT") for _ in range(5000))
        seq = seq[:3000] + seq[1000:1200] + seq[3200:]          # a second copy of [1000, 1200) at 3000
        res = ml.walk_expectations(seq, [(1050, 50, ["left", "right"], 300), (100, 40, ["right"], 100),
                                         (4950, 40, ["right"], 100)])
        self.assertEqual(res[0]["left"], {"min_bp": 50, "first": "ambiguous", "branch_by_bp": 50})   # copies diverge
        self.assertEqual(res[0]["right"]["min_bp"], 100)
        self.assertEqual(res[1]["right"], {"min_bp": 100, "first": "max_extension", "branch_by_bp": None})
        self.assertEqual(res[2]["right"], {"min_bp": 10, "first": "record_end", "branch_by_bp": None})
        rep = "ACGTTG" * 20 + "".join(rng.choice("ACGT") for _ in range(300))
        r = ml.walk_expectations(rep, [(0, 40, ["right"], 100)])[0]["right"]
        self.assertLess(r["min_bp"], 100)                         # a tandem repeat re-enters its own edges
        self.assertIn(r["first"], ("reentry", "ambiguous"))


class Builder(unittest.TestCase):
    def test_tags_and_repeats(self):
        import areas
        import build_manifest as bm
        reps = areas._tandem_repeats("G" * 50 + "CA" * 30 + "TTGACCATGA" + "A" * 45 + "C")
        self.assertEqual(reps[0][1:], (60, "CA"))
        self.assertEqual(sorted(u for _, _, u in reps), ["A", "CA", "G"])      # each once, at its smallest period

        class R:
            def __init__(self, acc, last, scan, n):
                self.accession, self.last_in_column, self.last_record_scan, self.records_in_column = acc, last, scan, n
        recs = [R("NC_1.1", True, "ok", 3), R("NC_2.1", False, None, 5), R("NC_3.1", True, None, 1)]
        scs = [{"tags": [], "truth": {"record": "NC_1.1", "genome": "x"}, "expect": [{"check": "status", "args": {}}]},
               {"tags": ["smoke"], "truth": {"a": {"record": "NC_2.1", "genome": "y"}, "b": {"record": "NC_3.1", "genome": "z"}},
                "expect": [{"check": "status", "args": {"http": 400, "text": "x"}}]},
               {"tags": [], "truth": {}, "expect": [{"check": "status", "args": {"http": 400, "text": "x"}}]}]
        bm.tag_scenarios(scs, recs)
        self.assertEqual(scs[0]["tags"], ["last-record-verified"])
        self.assertEqual(scs[0]["truth"]["last_record_scan"], "ok")
        self.assertEqual(scs[1]["tags"], ["last-record-verified", "smoke"])           # exercised by the smoke run
        self.assertEqual(scs[1]["truth"]["b"]["last_record_scan"], "single-record column")
        self.assertEqual(scs[2]["tags"], ["code-derived"])


class FrozenManifests(unittest.TestCase):
    """The committed manifests: their digest, the contract's fields, every check known to checks.py."""

    def test_manifests(self):
        found = 0
        for name in ("suite_v1.json", "suite_smoke_v1.json", "suite_v2.json", "suite_smoke_v2.json"):
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
            if "relations" in m:
                self.assertEqual(ml.sha256_json(m["relations"]), m["relations_sha256"], name)
                for r in m["relations"]:
                    self.assertIn(r["type"], checks.RELATIONS, r["id"])
                    self.assertIn(r["a"], ids, r["id"])
                    self.assertIn(r["b"], ids, r["id"])
        if found == 0:
            self.skipTest("no manifest built yet")


class CapabilityRelations(unittest.TestCase):
    """The three capabilities documents of the snapshot satisfy the manifest's same_answer relations between them."""

    def test_snapshot(self):
        import areas
        caps = {}
        for n in ("capabilities", "pattern_capabilities", "traverse_capabilities"):
            with open(os.path.join(SUITE, "truth", "staging_capabilities", n + ".json")) as f:
                caps[n] = json.load(f)
        ctx = ml.Ctx([], 1)
        scen = [{"id": "c1", "area": "C", "route": "GET /traverse/capabilities"},
                {"id": "c2", "area": "C", "route": "GET /pattern/capabilities"},
                {"id": "c3", "area": "C", "route": "GET /capabilities"}]
        docs = {"c1": caps["traverse_capabilities"], "c2": caps["pattern_capabilities"], "c3": caps["capabilities"]}
        # only the capabilities part of gen_REL: its relations need nothing else
        rels = []
        cid = {x["route"]: x["id"] for x in scen}
        for r in areas.capability_relations(ctx, cid, caps):
            rels.append(r)
        self.assertEqual(len(rels), 3)
        for r in rels:
            a = checks.normalise_direct(200, docs[r["a"]])
            b = checks.normalise_direct(200, docs[r["b"]])
            res = checks.evaluate_relation(r, scen[0], a, scen[0], b, "direct-local")
            self.assertEqual(res["outcome"], "pass", res)
        bad = json.loads(json.dumps(docs["c1"]))
        bad["pattern"]["caps"]["max_contexts"] = 1
        r = rels[2]
        res = checks.evaluate_relation(r, scen[0], checks.normalise_direct(200, docs["c2"]), scen[0],
                                       checks.normalise_direct(200, bad), "direct-local")
        self.assertEqual(res["outcome"], "fail")


if __name__ == "__main__":
    unittest.main()
