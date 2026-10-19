"""The deploy gate's area G (gate_policy.graphlet_block, gate.graphlet_problems, gate.library_identity_problems;
DESIGN-v8 §6) with synthetic manifests and records: every applicable way of every op call and every pre-check is
required, a disagreement, a failed pre-check, a precondition_failed and an inconclusive block, W3's not_exposed per the
manifest is no gap, the service checks the op plan exercises are required, library_commit is checked like
service_commit (--expect-library-commit) and against the manifest's pin, the owner's declared op calls are set aside
with their reason, the drift lane is report-only; a policy of a manifest without graphlet scenarios has no graphlets
block and a verdict without it states nothing new."""

import contextlib
import copy
import io
import json
import os
import shutil
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.dirname(HERE)
sys.path.insert(0, SUITE)
sys.path.insert(0, HERE)

import compare_graphlets as cg  # noqa: E402
import gate  # noqa: E402
import gate_policy  # noqa: E402
import graphlet_fixtures as gf  # noqa: E402
import merge_run  # noqa: E402
import test_gate as tg  # noqa: E402

GID = "G-merge-0001"
UUID = "0d5b9b1a-1111-2222-3333-444455556666"
WALKS = {"id": "walks-len", "tool": "graphlet_walks",
         "args": {"arm": "left", "rank": "length", "n": 10, "max_bytes": 2048}, "paging": "all"}
SPLITS = {"id": "splits", "tool": "graphlet_splits", "args": {"arm": "left", "n": 10, "max_bytes": 2048},
          "ways": ["W1", "W2", "W4", "W5"], "not_exposed": {"W3": "library-only (owner, 2026-10-19)"}}
LIST = {"id": "list", "tool": "graphlet_list", "args": {"max_bytes": 2048}, "ways": ["W1", "W2", "W4"],
        "not_exposed": {"W3": "store bookkeeping, a smoke (D-R15)"}}


def g_scenario(ops=None):
    """A graphlet scenario in test_gate's shape (its checks pass on the MCP arm) with a graphlet block."""
    base = tg.scenario(GID, "POST /traverse")
    g = gf.scenario(GID, ops=ops if ops is not None else [WALKS, SPLITS, LIST])
    base.update(request=g["request"], graphlet=g["graphlet"])
    return base


class GateG(unittest.TestCase):
    def setUp(self):
        gf.install_graphlet_doc_stub()
        gf.install_oracle_stub()
        self.tmp = tempfile.mkdtemp(prefix="lane_c_gate_")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def build(self, ops=None, records=True, change=None, name="run", **extra_records):
        self.m = tg.manifest_of(tg.scenarios() + [g_scenario(ops)])
        # a v8 run has the direct-local arm (W1's requests and the full twin); the deploy gate's arm is the MCP's
        d = tg.make_run(os.path.join(self.tmp, name), self.m, tg.mcp_lines(self.m), direct=("direct-local",))
        self.g = gf.GRun.attach(d, self.m)
        if records:
            self.g.complete_scenario(GID, ops=[])
            self.put_all(change or {})
        self.d = d
        return d

    def put_all(self, change):
        """A record per op and way (each op's ways), every answer the same; |change| {(op, way): fn(record kwargs)}."""
        mod, text = gf.install_graphlet_doc_stub(), self.g.fixture.doc(0)
        doc = {"body_digest": mod.body_digest(text), "envelope_digest": mod.envelope_digest(text)}
        scen = next(s for s in self.m["scenarios"] if s["id"] == GID)
        import graphlet_precheck as gp
        for op in gp.scenario_ops(scen):
            for way in gp.op_ways(scen, op):
                if op["tool"] == "graphlet_walks":
                    ans = gf.walks_answer()
                elif op["tool"] == "graphlet_list":
                    ans = {"total": 1, "rows": [{"handle": "g_1"}]}
                else:
                    ans = {"handle": "g_0123456789ab", "arm": "left", "total": 0, "rows": []}
                kw = {"answer": ans, "doc": doc, "service_fields": None}
                if way == "W3":
                    kw["answer"] = dict(ans, handle=UUID)
                    kw["service_fields"] = {"walk_display": "first_parent"}
                if way == "W5":
                    kw = {"answer": {}, "oracle": {"checked": {}, "consistency_only": ["/rows"]}}
                fn = change.get((op["id"], way))
                if fn == "skip":
                    continue
                if fn:
                    fn(kw)
                self.g.record(way, GID, op["id"], op["tool"], op["args"], **kw)

    def judge(self, pol=None, **kw):
        pol = pol or tg.policy_of(self.m)
        return gate.judge(self.d, pol, **kw)

    def kinds(self, v, kind):
        return [p for p in v["problems"] if p["kind"] == kind]


class Policy(GateG):
    def test_the_graphlets_block(self):
        self.build(records=False)
        pol = tg.policy_of(self.m)
        g = pol["graphlets"]
        self.assertEqual(g["scenarios"], [GID])
        self.assertEqual([(x["op"], x["ways"]) for x in g["required"]],
                         [("walks-len", ["W1", "W2", "W3", "W4", "W5"]), ("splits", ["W1", "W2", "W4", "W5"])])
        self.assertEqual(g["required"][1]["not_exposed"], {"W3": "library-only (owner, 2026-10-19)"})
        self.assertEqual(g["smoke_report_only"], [{"scenario": GID, "op": "list"}])
        self.assertEqual(g["prechecks"], ["body", "envelope", "twin", "strategy", "preconditions", "library_commit"])
        self.assertEqual(g["service_checks"], ["walk_display", "library_commit"])
        self.assertEqual((g["library_commit"], g["drift"]), (gf.PIN, "report_only"))
        self.assertEqual(pol["counts"]["graphlets"]["op_calls"], 2)
        text = gate_policy.describe(pol)
        self.assertIn("area G: 1 scenario(s), 2 op call(s) required", text)
        self.assertIn("W3 not_exposed: G-merge-0001 splits (graphlet_splits): library-only", text)

    def test_no_graphlet_scenario_no_block(self):
        m = tg.manifest_of()
        pol = tg.policy_of(m)
        self.assertNotIn("graphlets", pol)
        self.assertNotIn("graphlets", pol["counts"])
        self.assertNotIn("area G", gate_policy.describe(pol))

    def test_declared_op_calls_are_validated(self):
        self.build(records=False)
        good = {"declared_version": 1, "manifest": {"content_sha256": self.m["content_sha256"]},
                "graphlet_not_required": [{"scenario": GID, "op": "walks-len",
                                           "reason": "library defect SD-200, fixed before the freeze"}]}
        pol = tg.policy_of(self.m, declared=good)
        self.assertEqual(pol["graphlets"]["not_required"][0]["rules"], ["declared"])
        bad = dict(good, graphlet_not_required=[{"scenario": GID, "op": "nope", "reason": ""},
                                                {"scenario": "T-001", "op": "x", "reason": "r"}])
        with self.assertRaises(ValueError) as cm:
            tg.policy_of(self.m, declared=bad)
        for want in ("graphlet_not_required[0] states no reason", "has no op 'nope'",
                     "T-001 is not a selected graphlet scenario"):
            self.assertIn(want, str(cm.exception))


class Judge(GateG):
    def test_a_complete_area_g_passes(self):
        self.build()
        v = self.judge()
        self.assertTrue(v["passed"], gate.render(v))
        self.assertEqual(v["counts"]["graphlet_ops"]["agree"], 2)
        self.assertEqual(v["library_identity"]["recorded"]["library_commit"], gf.PIN)
        self.assertTrue(any(n.startswith("area G: 2 op call(s) required, 2 agree") for n in v["notes"]))
        self.assertIn("graphlet_inputs", v["inputs"])
        self.assertEqual(gate.stale_inputs(self.d, v), [])
        self.assertFalse(os.path.exists(os.path.join(self.d, cg.COMPARE_FILE)))   # computed in memory (VC3)

    def test_a_disagreement_blocks(self):
        self.build(change={("walks-len", "W2"): lambda kw: kw["answer"].update(total=3)})
        v = self.judge()
        self.assertFalse(v["passed"])
        p = self.kinds(v, "graphlet")
        self.assertEqual(len(p), 1, gate.render(v))
        self.assertEqual(p[0]["scenario"], GID)
        self.assertIn("op walks-len (graphlet_walks) disagrees (W1-W2)", p[0]["detail"])
        self.assertIn("/answer/total: W1 2, W2 3", p[0]["detail"])
        self.assertIn("graphlet: 1", gate.render(v))

    def test_a_missing_way_blocks_as_missing(self):
        self.build(change={("walks-len", "W5"): "skip"})
        v = self.judge()
        p = self.kinds(v, "missing")
        self.assertEqual([x["detail"] for x in p], ["op walks-len (graphlet_walks): no record of W5 (oracle)"])
        self.assertEqual(self.kinds(v, "graphlet"), [])

    def test_a_failed_precheck_and_an_inconclusive_block(self):
        self.build()
        self.g.write_echo(GID, dict(gf.strategy(), direction="left"))
        v = self.judge()
        self.assertTrue(any("pre-check strategy fail" in p["detail"] for p in self.kinds(v, "graphlet")),
                        gate.render(v))
        self.build(name="b", change={("splits", w): (lambda kw: kw.update(answer=None, selector_error="no split"))
                                     for w in ("W1", "W2", "W4", "W5")})
        v = self.judge()
        self.assertTrue(any("op splits (graphlet_splits) inconclusive: the selector failed in every way" in p["detail"]
                            for p in self.kinds(v, "graphlet")), gate.render(v))

    def test_a_body_mismatch_is_named_a_determinism_finding(self):
        self.build()
        f = self.g.fixture
        self.g.write_doc("W4", GID, f.doc(0).replace("\nG ", "\nG 9", 1))
        v = self.judge()
        self.assertTrue(any("pre-check body fail (a determinism finding)" in p["detail"]
                            for p in self.kinds(v, "graphlet")), gate.render(v))

    def test_precondition_failed_blocks(self):
        self.build()
        m = copy.deepcopy(self.m)
        next(s for s in m["scenarios"] if s["id"] == GID)["graphlet"]["preconditions"] = {"arms": ["right"]}
        # a manifest edit is another manifest: judged here through the comparison only
        sc_ = __import__("suite_common")
        m["content_sha256"] = sc_.scenarios_sha256(m["scenarios"])
        res, _ = cg.compare_run(self.d, m)
        self.assertEqual(res["scenarios"][GID]["ops"]["walks-len"]["outcome"], "precondition_failed")
        pol = tg.policy_of(self.m)
        orig = cg.load_or_compute
        try:
            cg.load_or_compute = lambda run_dir, manifest=None, reuse=False: dict(res, manifest=sc_.manifest_identity(self.m))
            v = self.judge(pol)
        finally:
            cg.load_or_compute = orig
        self.assertTrue(any("op walks-len (graphlet_walks) precondition_failed" in p["detail"]
                            for p in self.kinds(v, "graphlet")), gate.render(v))
        self.assertTrue(any("pre-check preconditions precondition_failed" in p["detail"]
                            for p in self.kinds(v, "graphlet")))

    def test_a_required_service_check_that_fails_blocks(self):
        self.build(change={("walks-len", "W3"): lambda kw: kw.update(service_fields={})})
        v = self.judge()
        self.assertTrue(any(p["detail"].startswith("service check walk_display fail")
                            for p in self.kinds(v, "graphlet")), gate.render(v))

    def test_declared_op_call_is_set_aside_with_its_reason(self):
        self.build(change={("walks-len", "W2"): lambda kw: kw["answer"].update(total=3)})
        declared = {"declared_version": 1, "manifest": {"content_sha256": self.m["content_sha256"]},
                    "graphlet_not_required": [{"scenario": GID, "op": "walks-len", "reason": "library defect SD-200"}]}
        v = self.judge(tg.policy_of(self.m, declared=declared))
        self.assertTrue(v["passed"], gate.render(v))
        self.assertTrue(any("set aside by the declared file: library defect SD-200" in n for n in v["notes"]))

    def test_the_drift_lane_is_a_note(self):
        self.build()
        self.g.record("drift", GID, "walks-len", "graphlet_walks", WALKS["args"], answer=dict(gf.walks_answer(),
                                                                                             total=9),
                      library={"commit": "d" * 40, "python": "3.13.1", "hashseed": 1, "work_model": 3})
        v = self.judge()
        self.assertTrue(v["passed"], gate.render(v))
        self.assertTrue(any("drift lane (report-only): 1 of 1" in n for n in v["notes"]))

    def test_a_stale_comparison_is_recomputed(self):
        self.build()
        self.assertTrue(self.judge()["passed"])
        self.g.record("W1", GID, "walks-len", "graphlet_walks", WALKS["args"], answer=dict(gf.walks_answer(), total=5))
        v = self.judge()
        self.assertFalse(v["passed"])
        self.assertTrue(any("W1-W2" in p["detail"] for p in self.kinds(v, "graphlet")))

    def test_a_comparison_of_another_manifest_is_an_identity_problem(self):
        self.build()
        cg.compute(self.d)
        path = os.path.join(self.d, cg.COMPARE_FILE)
        with open(path) as fh:
            data = json.load(fh)
        data["manifest"]["content_sha256"] = "0" * 64
        orig = cg.load_or_compute
        try:
            cg.load_or_compute = lambda run_dir, manifest=None, reuse=False: data
            v = self.judge()
        finally:
            cg.load_or_compute = orig
        self.assertTrue(any(p["detail"].startswith("graphlet_compare.json compared the manifest content 0000")
                            for p in self.kinds(v, "identity")))


class ArmDirs(GateG):
    def test_area_g_directories_are_no_arms(self):
        import graphlet_precheck as gp
        self.assertEqual(sorted(merge_run.GRAPHLET_DIRS),
                         sorted([d for w, d in gp.WAY_DIRS.items() if w != "W3"] + [gp.PRECHECK_DIR]))
        self.build()
        cg.compute(self.d)
        self.assertEqual([n for n, _ in merge_run.arm_dirs(self.d)], ["direct-local", "mcp"])


class ServiceChecksBlock(GateG):
    """VC2: a failed service check outside the policy's required list blocks; the declared file sets one aside."""
    SUMMARY = {"id": "summary", "tool": "graphlet_summary", "args": {"detail": "summary", "max_bytes": 2048},
               "ways": ["W1", "W2", "W3", "W4"]}

    def build_with_summary(self, in_service, name="run"):
        self.build(ops=[WALKS, self.SUMMARY], records=False, name=name)
        self.g.complete_scenario(GID, ops=[])
        self.put_all({("summary", "W3"): lambda kw: kw.update(service_fields={"in_service": in_service})})

    def put_all(self, change):
        import graphlet_precheck as gp
        mod, text = gf.install_graphlet_doc_stub(), self.g.fixture.doc(0)
        doc = {"body_digest": mod.body_digest(text), "envelope_digest": mod.envelope_digest(text)}
        scen = next(s for s in self.m["scenarios"] if s["id"] == GID)
        for op in gp.scenario_ops(scen):
            for way in gp.op_ways(scen, op):
                ans = gf.walks_answer() if op["tool"] == "graphlet_walks" else \
                    {"handle": "g_0123456789ab", "summary": {"seed": {}}, "evidence": {}}
                kw = {"answer": ans, "doc": doc, "service_fields": None}
                if way == "W3":
                    kw["answer"] = dict(ans, handle=UUID)
                    kw["service_fields"] = {"walk_display": "first_parent"}
                if way == "W5":
                    kw = {"answer": {}, "oracle": {"checked": {}}}
                fn = change.get((op["id"], way))
                if fn:
                    fn(kw)
                self.g.record(way, GID, op["id"], op["tool"], op["args"], **kw)

    def test_a_failed_check_the_policy_does_not_list_blocks(self):
        self.build_with_summary(False)   # in_service false on a body within 8 MiB
        pol = tg.policy_of(self.m)
        self.assertNotIn("in_service", pol["graphlets"]["service_checks"])
        v = self.judge(pol)
        self.assertTrue(any(p["detail"].startswith("service check in_service fail: 1 of 1 call(s) break it")
                            for p in self.kinds(v, "graphlet")), gate.render(v))
        # the owner sets it aside with a reason: a note
        declared = {"declared_version": 1, "manifest": {"content_sha256": self.m["content_sha256"]},
                    "graphlet_service_checks_not_required": [{"check": "in_service", "reason": "service bug SB-1"}]}
        v = self.judge(tg.policy_of(self.m, declared=declared))
        self.assertTrue(v["passed"], gate.render(v))
        self.assertTrue(any("service check in_service fail, set aside by the declared file: service bug SB-1" in n
                            for n in v["notes"]))
        bad = dict(declared, graphlet_service_checks_not_required=[{"check": "nope", "reason": ""}])
        with self.assertRaises(ValueError) as cm:
            tg.policy_of(self.m, declared=bad)
        self.assertIn("'nope' is not a service check", str(cm.exception))

    def test_a_check_whose_precondition_never_holds_is_na_and_passes(self):
        self.build_with_summary(True)
        v = self.judge()
        self.assertTrue(v["passed"], gate.render(v))
        import report
        g = cg.load_or_compute(self.d)
        word, reason = report.graphlet_verdict(g)
        self.assertEqual(word, "PASS")
        self.assertIn("service checks: 2 pass, 7 n/a", reason)


class Recompute(GateG):
    """VC3: the gate recomputes the graphlet comparison from the run's files; --reuse only with every digest current."""

    def test_deleted_twin_and_job_record_are_seen(self):
        self.build()
        self.assertTrue(self.judge()["passed"])
        cg.compute(self.d)   # a written comparison the gate must not trust by default
        os.remove(os.path.join(self.d, "arms", "direct-local", "docs", GID + ".full.json.gz"))
        os.remove(os.path.join(self.d, "arms", "mcp", "graphlet_jobs.jsonl"))
        for reuse in (False, True):
            v = self.judge(reuse=reuse)
            details = [p["detail"] for p in self.kinds(v, "graphlet")]
            self.assertTrue(any(d.startswith("pre-check twin inconclusive") for d in details), (reuse, details))
            self.assertTrue(any(d.startswith("pre-check strategy inconclusive") for d in details), (reuse, details))

    def test_reuse_takes_a_current_comparison(self):
        self.build()
        written = cg.compute(self.d)
        seen = []
        orig = cg.compute
        try:
            cg.compute = lambda *a, **k: seen.append(1) or orig(*a, **k)
            self.assertTrue(self.judge(reuse=True)["passed"])
            self.assertEqual(seen, [])   # reused
            self.assertTrue(self.judge()["passed"])
            self.assertEqual(seen, [1])  # recomputed without --reuse
        finally:
            cg.compute = orig
        self.assertEqual(written["inputs"], cg.input_digests(self.d))


class LibraryIdentity(GateG):
    def test_expect_library_commit(self):
        self.build()
        v = self.judge(expect_library_commit=gf.PIN[:10])
        self.assertTrue(v["passed"], gate.render(v))
        self.assertEqual(v["library_identity"]["expected"], gf.PIN[:10])
        v = self.judge(expect_library_commit="f" * 40)
        details = [p["detail"] for p in self.kinds(v, "identity")]
        self.assertTrue(any("the manifest pins the library commit %s, not the deployed library commit ffff" % gf.PIN
                            in x for x in details), details)
        self.assertTrue(any("recorded the library commit %s (1 invocation(s)), not the deployed" % gf.PIN in x
                            for x in details))

    def test_deployed_other_than_the_pin_and_unrecorded(self):
        self.build()
        os.remove(os.path.join(self.d, "arms", "mcp", "mcp_runs.jsonl"))
        self.g.write_mcp_runs(library_commit="c" * 40)
        v = self.judge()
        self.assertTrue(any("not the manifest's pin" in p["detail"] for p in self.kinds(v, "identity")),
                        gate.render(v))
        os.remove(os.path.join(self.d, "arms", "mcp", "mcp_runs.jsonl"))
        self.g.write_mcp_runs(library_commit=None)
        v = self.judge()
        self.assertTrue(any("1 of 1 MCP invocation(s) record no library_commit" in p["detail"]
                            for p in self.kinds(v, "identity")), gate.render(v))

    def test_unit_rules(self):
        lip = gate.library_identity_problems
        self.assertEqual(lip(None, None, ["mcp"]), ([], []))
        self.assertEqual(lip(None, "a" * 40, ["mcp"])[0], ["no MCP arm recorded: the library identity cannot be "
                                                           "confirmed"])
        unrec = {"recorded": False, "source": "s", "library_commits": {}, "invocations": 2, "without_commit": 2}
        self.assertEqual(lip(unrec, None, ["mcp"]), ([], []))   # a v7 run, no flag, no area G: nothing
        self.assertIn("records no library_commit", lip(unrec, None, ["mcp"], graphlets=True)[0][0])
        two = {"recorded": True, "source": "s", "library_commits": {"a" * 40: 1, "b" * 40: 1}, "invocations": 2,
               "without_commit": 0}
        self.assertIn("2 library commits", lip(two)[0][0])

    def test_a_v7_style_verdict_states_nothing_new(self):
        m = tg.manifest_of()
        d = tg.make_run(os.path.join(self.tmp, "v7"), m, tg.mcp_lines(m))
        v = gate.judge(d, tg.policy_of(m))
        self.assertTrue(v["passed"], gate.render(v))
        self.assertNotIn("library_identity", v)
        self.assertNotIn("graphlet_ops", v["counts"])
        self.assertNotIn("graphlet_inputs", v["inputs"])
        self.assertFalse(os.path.exists(os.path.join(d, cg.COMPARE_FILE)))

    def test_cli_flag(self):
        self.build()
        pol = os.path.join(self.tmp, "pol.json")
        with open(pol, "w") as fh:
            json.dump(tg.policy_of(self.m), fh)
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            self.assertEqual(gate.main([self.d, pol, "--expect-library-commit", "xyz"]), 2)
        self.assertIn("--expect-library-commit 'xyz' is not a git commit", err.getvalue())
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            self.assertEqual(gate.main([self.d, pol, "--expect-library-commit", gf.PIN,
                                        "--out", os.path.join(self.tmp, "v.json")]), 0)
        self.assertIn("GATE PASS", out.getvalue())


if __name__ == "__main__":
    unittest.main()
