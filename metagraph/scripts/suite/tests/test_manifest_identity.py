"""The manifest identity binds the relations too (GPT review finding 9, VERIFY-A): verify_manifest checks a stated
relations_sha256, run_direct records the relations' hash and a resume refuses other relations, merge_run names an
arm whose relations differ from the run's manifest, and the report shows both hashes."""

import contextlib
import io
import os
import signal
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.dirname(HERE)
sys.path.insert(0, SUITE)
sys.path.insert(0, HERE)

import merge_run  # noqa: E402
import report  # noqa: E402
import run_direct  # noqa: E402
import suite_common as sc  # noqa: E402
from fake_host import FakeHost  # noqa: E402
from test_compare_report import T0, line, write_arm  # noqa: E402
from test_runner import basic_scenarios  # noqa: E402

V4 = os.path.join(SUITE, "manifest", "suite_v4.json")


def with_relations(scenarios, relations):
    return {"manifest_version": 1, "content_sha256": sc.scenarios_sha256(scenarios), "created_utc": "x",
            "index": {"database": "test"}, "ground_truth": {}, "scenarios": scenarios, "relations": relations,
            "relations_sha256": sc.relations_sha256(relations)}


class VerifyManifestTest(unittest.TestCase):
    def setUp(self):
        self.m = sc.read_json(V4)

    def test_shipped_manifest_verifies(self):
        self.assertTrue(sc.verify_manifest(self.m)[0])
        self.assertEqual(sc.manifest_problems(self.m), [])
        ident = sc.manifest_identity(self.m)
        self.assertEqual(ident["content_sha256"], self.m["content_sha256"])
        self.assertEqual(ident["relations_sha256"], self.m["relations_sha256"])

    def test_deleted_relations_fail_verification(self):
        self.m["relations"] = []
        ok, computed = sc.verify_manifest(self.m)
        self.assertFalse(ok)
        # the scenarios' hash is still the one the arms record as computed_sha256
        self.assertEqual(computed, self.m["content_sha256"])
        self.assertIn("relations_sha256", " ".join(sc.manifest_problems(self.m)))

    def test_swapped_relation_fails_verification(self):
        rel = next(r for r in self.m["relations"] if r["type"] == "prefix_equals")
        rel["a"], rel["b"] = rel["b"], rel["a"]
        self.assertFalse(sc.verify_manifest(self.m)[0])

    def test_manifest_without_stated_relations_hash(self):
        # a manifest of before relations_sha256 verifies on its scenarios; its identity still names its relations
        m = {"scenarios": [{"id": "a"}], "content_sha256": sc.scenarios_sha256([{"id": "a"}])}
        self.assertTrue(sc.verify_manifest(m)[0])
        self.assertEqual(sc.manifest_identity(m)["relations_sha256"], sc.relations_sha256([]))

    def test_cli_writes_the_relations_hash(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "m.json")
            scen = [{"id": "a"}]
            rels = [{"id": "REL-1", "type": "same_answer", "a": "a", "b": "a", "args": {}}]
            sc.write_json_atomic(p, {"scenarios": scen, "content_sha256": "old", "relations": rels,
                                     "relations_sha256": "old"})
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(sc._main(["manifest-sha", p, "--write"]), 0)
            m = sc.read_json(p)
            self.assertEqual(m["relations_sha256"], sc.relations_sha256(rels))
            self.assertTrue(sc.verify_manifest(m)[0])


class ResumeRefusesOtherRelationsTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.d = self.tmp.name
        self.mpath = os.path.join(self.d, "m.json")

    def tearDown(self):
        signal.signal(signal.SIGINT, signal.default_int_handler)
        signal.signal(signal.SIGTERM, signal.SIG_DFL)
        if hasattr(signal, "SIGHUP"):
            signal.signal(signal.SIGHUP, signal.SIG_DFL)
        self.tmp.cleanup()

    def run_direct(self, base, out, *extra):
        return run_direct.main(["--manifest", self.mpath, "--arm", "direct-local", "--base", base, "--out", out,
                                "--pause", "0", "--quiet", "--checks", "stub"] + list(extra))

    def test_resume_refuses_changed_relations(self):
        scen = basic_scenarios()
        rels = [{"id": "REL-1", "type": "same_answer", "a": "R-1", "b": "D-1", "args": {}}]
        sc.write_json_atomic(self.mpath, with_relations(scen, rels))
        out = os.path.join(self.d, "run")
        with FakeHost() as host:
            self.run_direct(host.base, out, "--limit", "1")
            run = sc.read_json(os.path.join(out, "run.json"))
            self.assertEqual(run["manifest"]["computed_relations_sha256"], sc.relations_sha256(rels))
            self.assertEqual(run["manifest"]["relations_sha256"], sc.relations_sha256(rels))
            # the same scenarios, the relation's sides swapped: another test plan under the same scenarios' hash
            swapped = [dict(rels[0], a="D-1", b="R-1")]
            sc.write_json_atomic(self.mpath, with_relations(scen, swapped))
            with self.assertRaises(SystemExit) as ctx:
                self.run_direct(host.base, out)
            self.assertIn("other relations", str(ctx.exception))
            # an arm of before the field: the relations of the manifest copy it kept decide
            run.pop("manifest")
            run["manifest"] = {"computed_sha256": sc.scenarios_sha256(scen)}
            sc.write_json_atomic(os.path.join(out, "run.json"), run)
            with self.assertRaises(SystemExit):
                self.run_direct(host.base, out)
            # the run's own relations resume
            sc.write_json_atomic(self.mpath, with_relations(scen, rels))
            self.assertEqual(self.run_direct(host.base, out), 0)

    def test_first_run_warns_on_a_stale_relations_hash(self):
        scen = basic_scenarios()[:1]
        m = with_relations(scen, [])
        m["relations_sha256"] = "stale"
        sc.write_json_atomic(self.mpath, m)
        out = os.path.join(self.d, "run")
        with FakeHost() as host:
            self.run_direct(host.base, out)
        run = sc.read_json(os.path.join(out, "run.json"))
        self.assertFalse(run["manifest"]["sha_verified"])
        self.assertIn("relations_sha256 stale does not match", " ".join(run["warnings"]))


class MergeIdentityTest(unittest.TestCase):
    def arms(self, d, mcp_relations):
        scen = [{"id": "P-1", "area": "P"}]
        rels = [{"id": "REL-1", "type": "same_answer", "a": "P-1", "b": "P-1", "args": {}}]
        for arm, r in (("direct-local", rels), ("mcp", mcp_relations)):
            write_arm(d, arm, [line("P-1", arm, 0, T0, 1.0)])
            sc.write_json_atomic(os.path.join(d, "arms", arm, "manifest.json"), with_relations(scen, r))
        return rels

    def test_other_relations_in_an_arm_are_named(self):
        with tempfile.TemporaryDirectory() as d:
            rels = self.arms(d, [])
            merged = merge_run.merge(d)
            ident = merged["manifest_identity"]
            self.assertEqual(ident["manifest"]["relations_sha256"], sc.relations_sha256(rels))
            self.assertEqual(ident["arms"]["mcp"]["source"], "manifest copy")
            self.assertEqual(len(ident["differs"]), 1)
            self.assertIn("arms/mcp relations_sha256", ident["differs"][0])
            self.assertIn("different manifests", " ".join(merged["warnings"]))
            text = report.build_report(d)
            self.assertIn("relations sha256 %s" % sc.relations_sha256(rels)[:16], text)
            self.assertIn("the arms ran different manifests", text)

    def test_equal_identities_are_quiet(self):
        with tempfile.TemporaryDirectory() as d:
            rels = self.arms(d, None)
            for arm in ("direct-local", "mcp"):
                m = sc.read_json(os.path.join(d, "arms", arm, "manifest.json"))
                m["relations"], m["relations_sha256"] = rels, sc.relations_sha256(rels)
                sc.write_json_atomic(os.path.join(d, "arms", arm, "manifest.json"), m)
            merged = merge_run.merge(d)
            self.assertEqual(merged["manifest_identity"]["differs"], [])
            self.assertNotIn("different manifests", " ".join(merged["warnings"]))

    def test_strict_identity_is_an_error(self):
        with tempfile.TemporaryDirectory() as d:
            self.arms(d, [])
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()) as err:
                self.assertEqual(merge_run.main([d]), 0)
                self.assertEqual(merge_run.main([d, "--strict-identity"]), 3)
            self.assertIn("manifest identity: arms/mcp relations_sha256", err.getvalue())
            # the files are written all the same (the report shows the run)
            self.assertTrue(os.path.exists(os.path.join(d, "run.json")))

    def test_a_stale_stated_hash_of_the_run_manifest_is_a_warning(self):
        with tempfile.TemporaryDirectory() as d:
            self.arms(d, None)
            m = sc.read_json(os.path.join(d, "arms", "direct-local", "manifest.json"))
            m["relations"] = []   # relations_sha256 still states the deleted relation
            sc.write_json_atomic(os.path.join(d, "manifest.json"), m)
            merged = merge_run.merge(d)
            self.assertIn("manifest.json: relations_sha256", " ".join(merged["warnings"]))


if __name__ == "__main__":
    unittest.main()
