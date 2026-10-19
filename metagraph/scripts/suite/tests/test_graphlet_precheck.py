"""The pre-checks of an area G scenario (graphlet_precheck.py; DESIGN-v8 §2-§4, FORMAT-v8 §5) on synthetic runs whose
documents are the real mini_refseq graphlets: body, envelope and twin identity (a failure is a determinism finding),
the strategy echo, the structural preconditions (precondition_failed) and library_commit against the deployed one
and the manifest's pin, each with its outcome in the suite's vocabulary and the difference named."""

import contextlib
import copy
import gzip
import io
import json
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import graphlet_fixtures as gf  # noqa: E402
import graphlet_precheck as gp  # noqa: E402

SID = "G-merge-0001"


class PrecheckBase(unittest.TestCase):
    def setUp(self):
        gf.install_graphlet_doc_stub()
        self.run_ = gf.GRun()

    def tearDown(self):
        self.run_.cleanup()

    def check(self, name=None):
        recs = gp.prechecks(self.run_.dir)
        self.assertEqual(len(recs), 1)
        return recs[0] if name is None else recs[0][name]


class AllPass(PrecheckBase):
    def test_every_check_passes_and_seed_index_is_masked(self):
        self.run_.complete_scenario()   # W2's document has seed_index 3, W1's and W4's 0
        rec = self.check()
        for name in gp.CHECKS:
            self.assertEqual(rec[name]["outcome"], "pass", (name, rec[name]))
        self.assertEqual(rec["outcome"], "pass")
        self.assertIsNone(rec["finding"])
        self.assertEqual(rec["documents"]["W1"], "arms/dl-lib/docs/%s.mgt.gz" % SID)

    def test_cli_writes_one_line_per_scenario(self):
        self.run_.complete_scenario()
        with contextlib.redirect_stdout(io.StringIO()) as out:
            self.assertEqual(gp.main([self.run_.dir]), 0)
        self.assertIn("1 scenario(s) pass 1", out.getvalue())
        with open(gp.precheck_path(self.run_.dir)) as fh:
            lines = [json.loads(x) for x in fh]
        self.assertEqual([x["scenario"] for x in lines], [SID])
        self.assertEqual(set(gp.CHECKS) - set(lines[0]), set())


class Body(PrecheckBase):
    def test_a_differing_record_fails_as_determinism_and_is_named(self):
        f = self.run_.fixture
        text = f.doc(0)
        lines = text.split("\n")
        g = next(i for i, x in enumerate(lines) if x.startswith("G "))
        lines[g] = lines[g] + "X"   # a changed G record (its bases)
        self.run_.complete_scenario(docs={"W1": f.doc(0), "W2": f.doc(0), "W4": "\n".join(lines)})
        rec = self.check()
        self.assertEqual(rec["body"]["outcome"], "fail")
        self.assertEqual(rec["body"]["pair"], "W1-W4")
        self.assertEqual(rec["body"]["difference"]["type"], "G")
        self.assertIn("record %d (G)" % g, rec["body"]["detail"])   # text line g; J is no body record
        self.assertEqual((rec["outcome"], rec["finding"]), ("fail", "determinism"))

    def test_a_way_without_digest_is_inconclusive(self):
        self.run_.complete_scenario()
        path = os.path.join(self.run_.dir, "arms", "dr-lib", "graphlet_ops.jsonl")
        with open(path) as fh:
            rows = [json.loads(x) for x in fh]
        for r in rows:
            r["doc"] = None
        with open(path, "w") as fh:
            fh.write("".join(json.dumps(r) + "\n" for r in rows))
        rec = self.check()
        self.assertEqual(rec["body"]["outcome"], "inconclusive")
        self.assertIn("W2", rec["body"]["detail"])
        self.assertEqual(rec["outcome"], "inconclusive")

    def test_two_documents_in_one_way_fail(self):
        digests = self.run_.complete_scenario()
        self.run_.record("W1", SID, "other", "graphlet_summary", {"max_bytes": 2048}, answer={},
                         doc={"body_digest": "f" * 64, "envelope_digest": digests["W1"]["envelope_digest"]})
        rec = self.check("body")
        self.assertEqual(rec["outcome"], "fail")
        self.assertIn("W1 states 2 different body_digests", rec["detail"])

    def test_a_recorded_digest_the_kept_document_does_not_reproduce_fails(self):
        self.run_.complete_scenario()
        f = self.run_.fixture
        self.run_.write_doc("W2", SID, f.doc(0).replace("\nG ", "\nG 0", 1))   # kept document changed afterwards
        rec = self.check("body")
        self.assertEqual(rec["outcome"], "fail")
        self.assertIn("W2's records state body_digest", rec["detail"])

    def test_digests_differ_without_documents_still_fail(self):
        self.run_.complete_scenario()
        for way in ("W1", "W2", "W4"):
            os.remove(os.path.join(self.run_.dir, "arms", gf.WAY_DIRS[way], "docs", SID + ".mgt.gz"))
        path = os.path.join(self.run_.dir, "arms", "mcp-lib", "graphlet_ops.jsonl")
        with open(path) as fh:
            rows = [json.loads(x) for x in fh]
        rows[0]["doc"]["body_digest"] = "0" * 64
        with open(path, "w") as fh:
            fh.write("".join(json.dumps(r) + "\n" for r in rows))
        rec = self.check()
        self.assertEqual(rec["body"]["outcome"], "fail")
        self.assertIn("not both kept", rec["body"]["detail"])
        self.assertEqual(rec["twin"]["outcome"], "inconclusive")


class Envelope(PrecheckBase):
    def test_timing_alone_differs_passes(self):
        f = self.run_.fixture
        self.run_.complete_scenario(docs={"W1": f.doc(0), "W2": f.doc(0, timing={"elapsed_ms": 99}),
                                          "W4": f.doc(0, timing={"elapsed_ms": 7})})
        self.assertEqual(self.check("envelope")["outcome"], "pass")

    def test_a_feature_level_difference_is_named(self):
        f = self.run_.fixture
        caps = copy.deepcopy(f.envelope["capabilities"])
        caps["feature_level"] = 5
        self.run_.complete_scenario(docs={"W1": f.doc(0), "W2": f.doc(0), "W4": f.doc(0, capabilities=caps)})
        rec = self.check()
        self.assertEqual(rec["envelope"]["outcome"], "fail")
        self.assertEqual(rec["envelope"]["difference"]["pointer"], "/capabilities/feature_level")
        self.assertIn("/capabilities/feature_level", rec["envelope"]["detail"])
        self.assertEqual(rec["body"]["outcome"], "pass")
        self.assertEqual((rec["outcome"], rec["finding"]), ("fail", "determinism"))

    def test_a_document_without_its_envelope_fails(self):
        f = self.run_.fixture
        self.run_.complete_scenario(docs={"W1": f.doc(0), "W2": f.doc(0), "W4": f.doc(0, envelope=False)})
        rec = self.check()
        self.assertEqual(rec["envelope"]["outcome"], "fail")
        self.assertIn("W4's document carries no J line", rec["envelope"]["detail"])
        self.assertEqual(rec["body"]["outcome"], "pass")   # the body identity leaves J out


class Twin(PrecheckBase):
    def test_a_path_length_that_differs_fails_naming_it(self):
        full = copy.deepcopy(self.run_.fixture.full)
        full["results"][0]["arms"]["right"]["paths"][0]["length_bp"] += 1
        self.run_.complete_scenario(twin=False)
        self.run_.write_twin(SID, full)
        rec = self.check()
        self.assertEqual(rec["twin"]["outcome"], "fail")
        self.assertEqual(rec["twin"]["difference"]["arm"], "right")
        self.assertEqual(rec["twin"]["difference"]["pointer"], "/0/length_bp")
        self.assertEqual((rec["outcome"], rec["finding"]), ("fail", "determinism"))

    def test_an_end_reason_or_segment_that_differs_fails(self):
        full = copy.deepcopy(self.run_.fixture.full)
        full["results"][0]["arms"]["left"]["paths"][0]["end_reasons"] = {"max_extension": 1}
        self.run_.complete_scenario(twin=False)
        self.run_.write_twin(SID, full)
        self.assertEqual(self.check("twin")["difference"]["what"], "paths")
        full = copy.deepcopy(self.run_.fixture.full)
        full["results"][0]["arms"]["left"]["segments"][0]["parents"] = [7]
        self.run_.write_twin(SID, full)
        d = self.check("twin")["difference"]
        self.assertEqual((d["what"], d["pointer"]), ("segments", "/0/parents/0"))

    def test_no_twin_or_a_refused_one_is_inconclusive(self):
        self.run_.complete_scenario(twin=False)
        self.assertEqual(self.check("twin")["outcome"], "inconclusive")
        self.run_.write_twin(SID, {"http": 503, "body": None})
        rec = self.check("twin")
        self.assertEqual(rec["outcome"], "inconclusive")
        self.assertIn("http 503", rec["detail"])

    def test_every_real_fixture_reads_as_its_twin(self):
        import checks
        for name in sorted(os.listdir(gf.MINI)):
            if not name.endswith(".graphlet.json.gz"):
                continue
            g = gf.load_fixture(name[:-len(".graphlet.json.gz")], "graphlet")
            f = gf.load_fixture(name[:-len(".graphlet.json.gz")], "full")
            self.assertIsNone(gp.twin_difference(g["results"][0]["graphlet"], f, checks), name)


class Strategy(PrecheckBase):
    def test_a_drift_of_the_echo_is_named(self):
        self.run_.complete_scenario(echo=False)
        echo = copy.deepcopy(self.run_.scenarios[0]["request"]["strategy"])
        echo["bounds"]["time_budget_ms"] = 12000
        self.run_.write_echo(SID, echo)
        rec = self.check()
        self.assertEqual(rec["strategy"]["outcome"], "fail")
        self.assertEqual(rec["strategy"]["drifts"][0]["path"], "bounds.time_budget_ms")
        self.assertIn("strategy drift at bounds.time_budget_ms: the MCP echo 12000, request 15000",
                      rec["strategy"]["detail"])
        self.assertEqual((rec["outcome"], rec["finding"]), ("fail", "strategy drift"))

    def test_numbers_by_value_and_the_echo_from_the_mcp_raw_record(self):
        self.run_.complete_scenario(echo=False)
        echo = copy.deepcopy(self.run_.scenarios[0]["request"]["strategy"])
        echo["bounds"]["time_budget_ms"] = 15000.0
        raw = {"calls": [{"tool": "traverse_fetch", "payload": {"status": "pending"}},
                         {"tool": "traverse_status", "payload": {"status": "done", "overrides": echo}}]}
        with gzip.open(self.run_.path("arms", "mcp", "raw", SID + ".mcp.0.json.gz"), "wt") as fh:
            json.dump(raw, fh)
        self.run_.append(os.path.join("arms", "mcp", "executions.jsonl"),
                         {"scenario": SID, "arm": "mcp", "repeat": 0, "raw": "raw/%s.mcp.0.json.gz" % SID})
        rec = self.check("strategy")
        self.assertEqual(rec["outcome"], "pass", rec)
        self.assertIn("raw record", rec["echo_source"])

    def test_no_echo_is_inconclusive_and_a_twin_drift_fails(self):
        self.run_.complete_scenario(echo=False)
        self.assertEqual(self.check("strategy")["outcome"], "inconclusive")
        scen = copy.deepcopy(self.run_.scenarios[0])
        scen["graphlet"]["twin"]["strategy"]["bounds"]["max_extension_bp"] = 5
        rec = gp.strategy_check(self.run_.dir, scen)
        self.assertEqual(rec["outcome"], "fail")
        self.assertEqual(rec["drifts"][0]["what"], "twin")


class Preconditions(PrecheckBase):
    def scen(self, pre):
        scen = copy.deepcopy(self.run_.scenarios[0])
        scen["graphlet"]["preconditions"] = pre
        return scen

    def test_facts_on_w1s_document(self):
        import checks
        self.run_.complete_scenario()
        docs = {"W1": self.run_.fixture.doc(0)}
        full = self.run_.fixture.full["results"][0]
        ok = {"arms": ["left", "right"], "walks_min": {"left": len(full["arms"]["left"]["paths"])},
              "segments_min": 2, "merges_min": 0, "labels_present": [full["label_dict"][0]["name"]],
              "walk_bp_min": {"right": 10}, "bases": True}
        rec = gp.preconditions_check(self.scen(ok), docs, checks)
        self.assertEqual(rec["outcome"], "pass", rec)
        bad = dict(ok, walks_min={"left": 999}, labels_present=["NO_SUCH_LABEL"])
        rec = gp.preconditions_check(self.scen(bad), docs, checks)
        self.assertEqual(rec["outcome"], "precondition_failed")
        self.assertIn("walks_min.left: want 999", rec["detail"])
        self.assertIn('labels_present: want "NO_SUCH_LABEL", seen false', rec["detail"])
        rec = gp.preconditions_check(self.scen({"arms": ["left", "right"], "colour": "blue"}), docs, checks)
        self.assertEqual(rec["outcome"], "inconclusive")
        self.assertIn("colour", rec["detail"])

    def test_a_failed_precondition_is_the_scenario_outcome(self):
        self.run_.scenarios[0]["graphlet"]["preconditions"] = {"arms": ["right"]}
        self.run_.manifest = gf.manifest(self.run_.scenarios)
        self.run_.write_json("manifest.json", self.run_.manifest)
        self.run_.complete_scenario()
        rec = self.check()
        self.assertEqual(rec["preconditions"]["outcome"], "precondition_failed")
        self.assertEqual((rec["outcome"], rec["finding"]), ("precondition_failed", None))


class LibraryCommit(PrecheckBase):
    def test_another_commit_on_a_w3_call_fails(self):
        self.run_.complete_scenario()
        self.run_.record("W3", SID, "walks-len", "graphlet_walks", {"arm": "left"}, answer={},
                         library={"commit": "b" * 40})
        rec = self.check()
        self.assertEqual(rec["library_commit"]["outcome"], "fail")
        self.assertIn("walks-len: bbbbbbbbbbbb", rec["library_commit"]["detail"])
        self.assertEqual((rec["outcome"], rec["finding"]), ("fail", "library_commit"))

    def test_deployed_other_than_the_pin_fails_and_an_abbreviation_matches(self):
        self.run_.complete_scenario(mcp_runs=False)
        self.run_.write_mcp_runs(library_commit="c" * 40)
        self.assertEqual(self.check("library_commit")["outcome"], "fail")
        os.remove(os.path.join(self.run_.dir, "arms", "mcp", "mcp_runs.jsonl"))
        self.run_.write_mcp_runs(library_commit=gf.PIN[:9])
        self.assertEqual(self.check("library_commit")["outcome"], "pass")

    def test_no_pin_or_no_deployed_record_is_inconclusive(self):
        self.run_.scenarios[0]["graphlet"].pop("library_commit")
        self.run_.manifest = gf.manifest(self.run_.scenarios)
        self.run_.write_json("manifest.json", self.run_.manifest)
        self.run_.complete_scenario()
        rec = self.check("library_commit")
        self.assertEqual(rec["outcome"], "inconclusive")
        self.assertIn("pins no library_commit", rec["detail"])
        # a top-level pin would be outside content_sha256: not read (FORMAT v2 D1)
        self.run_.manifest = gf.manifest(self.run_.scenarios, graphlets={"library_commit": gf.PIN})
        self.run_.write_json("manifest.json", self.run_.manifest)
        self.assertEqual(self.check("library_commit")["outcome"], "inconclusive")

    def test_without_a_w3_op_it_is_not_exposed(self):
        self.run_.scenarios[0]["graphlet"]["ways"]["default"] = ["W1", "W2", "W4", "W5"]
        self.run_.manifest = gf.manifest(self.run_.scenarios)
        self.run_.write_json("manifest.json", self.run_.manifest)
        self.run_.complete_scenario()
        self.assertEqual(self.check("library_commit")["outcome"], "not_exposed")


class FormatV2(PrecheckBase):
    def test_w4_reads_the_mcps_own_exported_file(self):
        self.run_.complete_scenario()
        src = os.path.join(self.run_.dir, "arms", "mcp-lib", "docs", SID + ".mgt.gz")
        os.makedirs(os.path.join(self.run_.dir, "arms", "mcp", "docs"), exist_ok=True)
        os.replace(src, os.path.join(self.run_.dir, "arms", "mcp", "docs", SID + ".mgt.gz"))
        rec = self.check()
        self.assertEqual(rec["documents"]["W4"], "arms/mcp/docs/%s.mgt.gz" % SID)
        self.assertEqual(rec["body"]["outcome"], "pass")

    def test_the_job_record(self):
        self.run_.complete_scenario(echo=False)
        echo = copy.deepcopy(self.run_.scenarios[0]["request"]["strategy"])
        self.run_.append(os.path.join("arms", "mcp", "graphlet_jobs.jsonl"),
                         {"scenario": SID, "job_id": "j1", "handle": "0d5b9b1a-1111-2222-3333-444455556666",
                          "overrides_echo": echo, "service_commit": gf.SERVICE, "library_commit": gf.PIN})
        rec = self.check()
        self.assertEqual(rec["strategy"]["outcome"], "pass", rec["strategy"])
        self.assertEqual((rec["library_commit"]["outcome"], rec["library_commit"]["job"]), ("pass", gf.PIN))
        self.run_.append(os.path.join("arms", "mcp", "graphlet_jobs.jsonl"),
                         {"scenario": SID, "job_id": "j1", "overrides": echo, "library_commit": "b" * 40})
        rec = self.check("library_commit")
        self.assertEqual(rec["outcome"], "fail")
        self.assertIn("the scenario's MCP job reports the library_commit bbbbbbbbbbbb", rec["detail"])

    def test_a_failed_seed(self):
        scen = gf.scenario(pre={"no_graphlet": True})
        self.run_.cleanup()
        self.run_ = gf.GRun(scenarios=[scen])
        self.run_.write_twin(SID, {"results": [{"seed": {"length_bp": 31}, "error": "seed_failed"}]})
        self.run_.write_echo(SID, scen["request"]["strategy"])
        self.run_.write_mcp_runs()
        err = {"code": "no_graphlet", "fields": {}}
        for way in ("W1", "W2", "W3", "W4"):
            self.run_.record(way, SID, "walks-len", "graphlet_walks", scen["graphlet"]["ops"][0]["args"], error=err)
        rec = self.check()
        for name in ("body", "envelope", "twin", "preconditions"):
            self.assertEqual(rec[name]["outcome"], "pass", (name, rec[name]))
        self.assertEqual(rec["outcome"], "pass")
        # the seed has a graphlet now: precondition_failed
        self.run_.write_doc("W1", SID, self.run_.fixture.doc(0))
        rec = self.check()
        self.assertEqual(rec["preconditions"]["outcome"], "precondition_failed")
        self.assertIn("the seed has a graphlet now", rec["preconditions"]["detail"])


class Helpers(unittest.TestCase):
    def test_pointers_and_first_difference(self):
        doc = {"a/b": [{"c~d": 1}]}
        self.assertEqual(gp.pointer_get(doc, "/a~1b/0/c~0d"), (True, 1))
        self.assertEqual(gp.pointer_get(doc, "/a~1b/1"), (False, None))
        self.assertEqual(gp.pointer_get(doc, "/a~1b/01"), (False, None))
        self.assertIsNone(gp.first_difference({"x": 1, "y": [1, 2]}, {"y": [1, 2.0], "x": 1.0}))
        self.assertEqual(gp.first_difference({"x": True}, {"x": 1})[0], "/x")
        self.assertEqual(gp.first_difference([1, 2], [1])[0], "/1")
        p, a, b = gp.first_difference({"x": 1}, {})
        self.assertEqual((p, a, b is gp.MISSING), ("/x", 1, True))

    def test_body_records_mask_and_skip(self):
        text = gf.standalone("H mgt 1 31 basic $ACGT c t k 64 1000 4 walk x\nA r c\nZ 3\n", {"k": 1})
        self.assertEqual(gp.body_records(text), ["H mgt 1 31 basic $ACGT c t k 64 1000 <seed_index> walk x", "A r c"])
        self.assertEqual(gp.envelope_of(text), {"k": 1})
        self.assertIsNone(gp.envelope_of("H x\nA r\nZ 2\n"))

    def test_scenario_outcome_order(self):
        def recs(**outs):
            return dict((k, {"outcome": outs.get(k, "pass")}) for k in gp.CHECKS)
        self.assertEqual(gp.scenario_outcome(recs()), ("pass", None))
        self.assertEqual(gp.scenario_outcome(recs(preconditions="precondition_failed", twin="inconclusive")),
                         ("precondition_failed", None))
        self.assertEqual(gp.scenario_outcome(recs(preconditions="precondition_failed", envelope="fail")),
                         ("fail", "determinism"))
        self.assertEqual(gp.scenario_outcome(recs(library_commit="not_exposed")), ("pass", None))


if __name__ == "__main__":
    unittest.main()
