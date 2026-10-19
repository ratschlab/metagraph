"""Offline tests of run_direct.py's area G twin (manifest v8, FORMAT-v8 v2 C2, lane D) against the fake host: the
`detail: full` twin sent once on direct-local right after a graphlet scenario's first execution, kept as the bare
engine body at docs/<id>.full.json.gz (what graphlet_precheck and graphlet_oracle read) with its raw record; none on
direct-remote; a resume sends a missing twin alone; a refused twin is recorded, attempted once per session."""
import copy
import gzip
import json
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.dirname(HERE)
sys.path.insert(0, SUITE)
sys.path.insert(0, HERE)

import area_g  # noqa: E402
import suite_common as sc  # noqa: E402


class Twin(unittest.TestCase):
    """run_direct.py's area G twin (FORMAT-v8 v2 C2) against the fake host."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = self.tmp.name

    def tearDown(self):
        import signal
        signal.signal(signal.SIGINT, signal.default_int_handler)
        signal.signal(signal.SIGTERM, signal.SIG_DFL)
        if hasattr(signal, "SIGHUP"):
            signal.signal(signal.SIGHUP, signal.SIG_DFL)
        self.tmp.cleanup()

    @staticmethod
    def g(sid, route="POST /traverse", repeat=2):
        strategy = area_g.mcp_merge("walk", {"bounds": {"max_extension_bp": 100}})
        twin = copy.deepcopy(strategy)
        twin["output"]["detail"] = "full"
        return {"id": sid, "area": "G", "title": sid, "route": route,
                "request": {"seeds": [{"sequence": "ACGT" * 10}], "strategy": strategy},
                "mcp": {"not_implemented": "test"}, "truth": {}, "expect": [{"check": "status", "args": {"http": 200}}],
                "timing_class": "fast", "repeat": repeat, "tags": [],
                "graphlet": {"twin": {"strategy": twin}, "ops": []}}

    def run_arm(self, base, out, scenarios, arm="direct-local", *extra):
        import run_direct
        path = os.path.join(self.dir, "m.json")
        with open(path, "w") as fh:
            json.dump({"manifest_version": 8, "content_sha256": sc.scenarios_sha256(scenarios), "created_utc": "x",
                       "index": {"database": "test"}, "scenarios": scenarios}, fh)
        return run_direct.main(["--manifest", path, "--arm", arm, "--base", base, "--out", out, "--pause", "0",
                                "--quiet", "--checks", "stub", "--warm-every", "0"] + list(extra))

    def test_the_twin_is_sent_once_and_kept_as_the_bare_body(self):
        sys.path.insert(0, HERE)
        from fake_host import FakeHost
        out = os.path.join(self.dir, "dl")
        with FakeHost() as host:
            self.assertEqual(self.run_arm(host.base, out, [self.g("G-x-0001")]), 0)
            posts = [p for m, p in host.log if m == "POST"]
        self.assertEqual(posts, ["/traverse"] * 3)   # repeat 0, the twin, repeat 1
        with gzip.open(os.path.join(out, "docs", "G-x-0001.full.json.gz"), "rt") as fh:
            body = json.load(fh)
        self.assertIn("results", body)          # the bare engine body (graphlet_precheck.read_twin, the oracle)
        with gzip.open(os.path.join(out, "raw", "G-x-0001.direct-local.twin.json.gz"), "rt") as fh:
            raw = json.load(fh)
        self.assertEqual(raw["request"]["body"]["strategy"]["output"]["detail"], "full")
        self.assertEqual(raw["request"]["body"]["seeds"], [{"sequence": "ACGT" * 10}])
        run = sc.read_json(os.path.join(out, "run.json"))
        self.assertEqual(run["twins"]["G-x-0001"]["file"], "docs/G-x-0001.full.json.gz")
        rows, _ = sc.read_jsonl(os.path.join(out, "executions.jsonl"))
        self.assertEqual(sorted((r["scenario"], r["repeat"]) for r in rows), [("G-x-0001", 0), ("G-x-0001", 1)])
        # lane C's reader and lane O's read the kept twin
        import graphlet_precheck as gp
        os.makedirs(os.path.join(self.dir, "run", "arms"), exist_ok=True)
        os.symlink(out, os.path.join(self.dir, "run", "arms", "direct-local"))
        twin, why = gp.read_twin(os.path.join(self.dir, "run"), "G-x-0001")
        self.assertIsNotNone(twin, why)
        # a resume with the twin kept sends nothing; without it, the twin alone
        with FakeHost() as host:
            self.assertEqual(self.run_arm(host.base, out, [self.g("G-x-0001")]), 0)
            self.assertEqual([p for m, p in host.log if m == "POST"], [])
            os.remove(os.path.join(out, "docs", "G-x-0001.full.json.gz"))
            self.assertEqual(self.run_arm(host.base, out, [self.g("G-x-0001")]), 0)
            self.assertEqual([p for m, p in host.log if m == "POST"], ["/traverse"])
        self.assertTrue(os.path.exists(os.path.join(out, "docs", "G-x-0001.full.json.gz")))

    def test_no_twin_on_direct_remote_and_a_refused_twin_once(self):
        sys.path.insert(0, HERE)
        from fake_host import FakeHost
        with FakeHost() as host:
            out = os.path.join(self.dir, "dr")
            self.assertEqual(self.run_arm(host.base, out, [self.g("G-x-0001")], "direct-remote"), 0)
            self.assertEqual([p for m, p in host.log if m == "POST"], ["/traverse"] * 2)
            self.assertFalse(os.path.exists(os.path.join(out, "docs")))
            del host.log[:]
            out = os.path.join(self.dir, "dl")
            self.assertEqual(self.run_arm(host.base, out, [self.g("G-y-0001", route="POST /nowhere")]), 0)
            self.assertEqual([p for m, p in host.log if m == "POST"], ["/nowhere"] * 3)   # r0, the twin, r1
        run = sc.read_json(os.path.join(out, "run.json"))
        self.assertEqual(run["twins"]["G-y-0001"]["http"], 404)
        self.assertIsNone(run["twins"]["G-y-0001"]["file"])
        self.assertTrue(any("twin of G-y-0001" in w for w in run.get("warnings") or []))


if __name__ == "__main__":
    unittest.main()
