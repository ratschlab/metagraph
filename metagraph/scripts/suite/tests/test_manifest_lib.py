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


class V4Builder(unittest.TestCase):
    """Lane S5: the builder changes of the verified triage (the v4 manifest; built later with a refreshed snapshot)."""

    def tearDown(self):
        ml.set_host_capabilities(None)

    def test_case_naming(self):
        p = lambda **f: ml.pattern_by_design(ml.pattern_request([{"dna": "A" * 25}], **f))  # noqa: E731
        # G02: a cap of 0 (async's pattern.caps_zero)
        self.assertIn("pattern.caps_zero", p(mode="partial", max_contexts=0))
        self.assertIn("pattern.caps_zero", p(long_search="paths", max_anchors=0))
        self.assertIn("pattern.caps_zero", p(long_search="paths", max_paths=0))
        self.assertNotIn("pattern.caps_zero", p(mode="partial", max_contexts=1))
        # G07: a predicate form fault other than a non-string name
        bad = [{"at_least": {"n": 3, "labels": ["1", "2"]}}, {"at_least": {"n": 0, "labels": ["1"]}},
               {"at_least": {"n": 1, "labels": ["1"], "x": 1}}, {"any": []}, {"any": ["1", "1"]}, {"any": [""]},
               {"any": ["1"], "none": ["2"]}, {"xor": ["1"]}, {"and": []}, {"not": {"or": []}}, {}]
        deep = {"any": ["1"]}
        for _ in range(66):
            deep = {"not": deep}
        for pred in bad + [deep]:
            self.assertIn("pattern.predicate_form_other", p(predicate=pred), pred)
        good = [{"at_least": {"n": 2, "labels": ["1", "2"]}}, {"and": [{"any": ["562"]}, {"none": ["287"]}]},
                {"or": [{"all": ["546", "615"]}, {"not": {"any": ["287"]}}]}, {"not": {"any": ["1"]}}]
        for pred in good:
            self.assertNotIn("pattern.predicate_form_other", p(predicate=pred), pred)
        self.assertNotIn("pattern.predicate_form_other", p(predicate={"any": [562]}))   # name_not_string's case
        self.assertIn("pattern.predicate_name_not_string", p(predicate={"any": [562]}))
        # the predicate options without a predicate, split by value (async's 789fb5f0)
        self.assertIn("pattern.predicate_options_without_predicate", p(predicate_scope="motif"))
        self.assertIn("pattern.predicate_options_without_predicate", p(predicate_scope="motif",
                                                                       predicate_strands="either"))
        self.assertIn("pattern.predicate_options_noop_without_predicate", p(predicate_strands="context"))
        self.assertIn("pattern.predicate_options_noop_without_predicate", p(predicate_scope="context"))
        self.assertNotIn("pattern.predicate_options_noop_without_predicate", p(predicate_scope="motif"))
        # G05: max_mismatches above a kind's cap of the snapshot
        pep = lambda m: ml.pattern_by_design(ml.pattern_request([{"protein": "MKVLAAGIVG"}], max_mismatches=m))  # noqa
        self.assertNotIn("mismatch.above_kind_cap", pep(2))       # the SPEC's cap of 3 until a snapshot says
        ml.set_host_capabilities({"mismatches": {"max": {"dna": 3, "iupac": 3, "protein": 1}}})
        self.assertIn("mismatch.above_kind_cap", pep(2))
        self.assertNotIn("mismatch.above_kind_cap", pep(1))
        self.assertNotIn("mismatch.above_kind_cap", p(max_mismatches=2))   # dna under its cap 3
        long_pep = ml.pattern_by_design(ml.pattern_request([{"protein": "MKVLAAGIVGK"}], max_mismatches=2))
        self.assertNotIn("mismatch.above_kind_cap", long_pep)           # 33 bases: the long-pattern case
        self.assertEqual(ml.mcp_pattern(ml.pattern_request([{"protein": "MKVLAAGIVG"}], max_mismatches=2))["by_design"],
                         ["mismatch.above_kind_cap"])

    def test_resolve_and_traverse_naming(self):
        r = lambda seq="ACGT" * 25, **f: ml.resolve_by_design(ml.resolve_request(seq, **f))  # noqa: E731
        # G10: the builders' absent label
        self.assertIn("resolve.unknown_label", r(labels=["562", ml.ABSENT_LABEL]))
        self.assertNotIn("resolve.unknown_label", r(labels=["562"]))
        # splits by value (async's 789fb5f0): the old ids keep what suite_v3 sends
        self.assertIn("resolve.time_budget_below_minimum", r(labels=["562"], time_budget_ms=250))
        self.assertIn("resolve.time_budget_above_reserve_below_minimum", r(labels=["562"], time_budget_ms=300))
        both = r(labels=["562"], time_budget_ms=500)
        self.assertFalse({"resolve.time_budget_below_minimum", "resolve.time_budget_above_reserve_below_minimum"} & set(both))
        self.assertIn("resolve.discover_limit_above_100", r(discover={"max_labels": 0}))
        self.assertIn("resolve.discover_limit_above_100_host_answers", r(discover={"max_labels": 101}))
        self.assertFalse([i for i in r(discover={"max_labels": 100}) if i.startswith("resolve.discover_limit")])
        # G19: an empty seed list
        req = {"seeds": [], "strategy": ml.traverse_request([])["strategy"]}
        self.assertEqual(ml.mcp_traverse(req)["by_design"], ["traverse.empty_seeds"])
        self.assertNotIn("traverse.empty_seeds", ml.traverse_by_design(ml.traverse_request([{"sequence": "A" * 40}])))

    def test_traverse_presets_and_overrides(self):
        # G15: the trace preset sends on_reconverge keep, as the service's does; the MCP still lacks coordinates
        req = ml.traverse_request([{"sequence": "A" * 40, "labels": ["x"]}], {"bounds": {"max_extension_bp": 500}},
                                  preset="trace", extra={"output": {"coordinates": True}})
        self.assertEqual((req["strategy"]["support"], req["strategy"]["branching"]),
                         ("trace", {"max_label_branches": 2, "on_reconverge": "keep"}))
        m = ml.mcp_traverse(req)
        self.assertEqual(m["not_implemented"], "strategy.output.coordinates")
        self.assertIn("traverse.trace_preset", m["by_design"])
        # G14: labels.max_seed_labels is an override traverse_fetch sends (the service bounds it to 10,000)
        req = ml.traverse_request([{"sequence": "A" * 40}], {"labels": {"seed_label_kind": "header",
                                                                         "max_seed_labels": 10000}})
        self.assertEqual(ml.mcp_traverse(req)["arguments"]["overrides"]["labels"],
                         {"seed_label_kind": "header", "max_seed_labels": 10000})

    def test_planted_window_reads_the_composition(self):
        # G03: an AT-rich record's 20-mer at m 1 meets far more contexts than the uniform model's 834
        import random
        import areas
        import truth
        rnd = random.Random(3)

        def rec(key, weights):
            seq = "".join(rnd.choices("ACGT", weights=weights, k=20000))
            return truth.Record(key=key, label=key, accession=key + ".1", kingdom="b", compartment="c",
                                gc_class="x", transl_table=11, column="1", column_index=0, seq_id=0,
                                records_in_column=1, length=len(seq), last_in_column=False, last_record_scan=None,
                                region=None, seq=seq)
        at, even = rec("at", [38, 12, 12, 38]), rec("even", [25, 25, 25, 25])
        ctx = ml.Ctx([at, even], 1)
        self.assertLess(areas._expected_ball_contexts(20, 1), areas.MM_CONTEXTS_TARGET)    # 834: uniform passes it
        t = areas._planted_window(ctx, ctx.rng("a"), at, [20, 25, 31], 1, ctx.rng("b"))
        self.assertEqual(t["length"], 25)                      # the AT-rich 20-mer is rejected by its composition
        self.assertGreater(areas._composition_factor(at, "AT" * 10), 10)
        self.assertLess(t["chance_contexts_estimate"], areas.MM_CONTEXTS_TARGET)
        lens = [areas._planted_window(ctx, ctx.rng("c", i), even, [20, 25, 31], 1, ctx.rng("d", i))["length"]
                for i in range(20)]
        self.assertIn(20, lens)                                # a balanced window keeps the uniform choice
        self.assertEqual(areas._planted_window(ctx, ctx.rng("e"), at, [31], 3, ctx.rng("f"))["length"], 31)
        # a factor below 1 is not credited: the estimate never falls below the uniform model's
        self.assertGreaterEqual(t["chance_contexts_estimate"], round(areas._expected_ball_contexts(25, 1), 1))

    def test_build_refuses_a_frozen_manifest(self):
        import build_manifest as bm
        with self.assertRaises(SystemExit):
            bm.main(["--out", os.path.join(SUITE, "manifest", "suite_v3.json"), "--smoke-out", os.devnull])

    def test_m_area_follows_the_snapshot(self):
        # G04, G05, G21 on the M area built from the record cache (skipped without it)
        import build_manifest as bm
        import tempfile
        import truth
        if not os.path.isdir(truth.CACHE_DIR) or not os.listdir(truth.CACHE_DIR):
            self.skipTest("no record cache (SUITE_TRUTH_CACHE or truth/cache)")
        with tempfile.TemporaryDirectory() as d:
            for name in ("capabilities", "pattern_capabilities", "traverse_capabilities"):
                with open(os.path.join(SUITE, "truth", "staging_capabilities", name + ".json")) as f:
                    doc = json.load(f)
                if name == "pattern_capabilities":
                    doc["mismatches"]["max"]["protein"] = 1
                with open(os.path.join(d, name + ".json"), "w") as f:
                    json.dump(doc, f)
            cap1 = bm.build(20261013, fetch=False, only_areas={"M"}, capabilities_dir=d)
        cap3 = bm.build(20261013, fetch=False, only_areas={"M"})
        s1 = dict((s["id"], s) for s in cap1["scenarios"])
        s3 = dict((s["id"], s) for s in cap3["scenarios"])
        self.assertEqual(cap1["index"]["mismatches_max"]["protein"], 1)
        # cap 1: no m = 2 peptide (ids kept, skipped with the reason), the cap slot named by design
        self.assertNotIn("M-protein-0008", s1)
        self.assertIn("M-protein-0008", " ".join(x["scenario"] for x in cap1["skipped"]))
        self.assertEqual(s1["M-protein-0010"]["mcp"]["by_design"], ["mismatch.above_kind_cap"])
        self.assertEqual(s1["M-protein-0010"]["expect"][1]["args"]["fields"], {"max_mismatches_cap": 1})
        # cap 3 (the stale snapshot): the m = 2 peptides, no cap slot
        self.assertIn("M-protein-0008", s3)
        self.assertNotIn("M-protein-0010", s3)
        # G21: every planted and protein search allows the 503 deadline and is machine-dependent
        for sid, s in s1.items():
            if sid.startswith(("M-planted", "M-protein")) and sid != "M-protein-0010":
                self.assertEqual(s["expect"][0], {"check": "status", "args": {"http": [200, 503]}}, sid)
                self.assertIn("machine-dependent", s["tags"], sid)
            if sid.startswith("M-planted"):
                self.assertLess(s["truth"]["chance_contexts_estimate"], 1000.0, sid)
        # G07, G10, G12 on the R and Q areas
        rq = dict((s["id"], s) for s in bm.build(20261013, fetch=False, only_areas={"R", "Q"})["scenarios"])
        self.assertEqual(rq["Q-more-0003"]["expect"], [{"check": "status", "args": {
            "http": 400, "code": "invalid_request", "text_contains": ["request.predicate.at_least.n"]}}])
        self.assertIn("pattern.predicate_form_other", rq["Q-more-0003"]["mcp"]["by_design"])
        self.assertIn("resolve.unknown_label", rq["R-edge-0008"]["mcp"]["by_design"])
        body = rq["R-edge-0019"]["body_raw"]
        self.assertTrue(body.startswith("{") and not body.rstrip().endswith("}"))
        with self.assertRaises(ValueError):
            json.loads(body)


class FrozenManifests(unittest.TestCase):
    """The committed manifests: their digest, the contract's fields, every check known to checks.py."""

    def test_manifests(self):
        found = 0
        for name in ("suite_v1.json", "suite_smoke_v1.json", "suite_v2.json", "suite_smoke_v2.json", "suite_v3.json",
                     "suite_smoke_v3.json"):
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
                    if m.get("manifest_version", 1) >= 3:   # built against this checks.py: every argument known
                        self.assertEqual(checks.validate_expectation(e), [], s["id"])
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
