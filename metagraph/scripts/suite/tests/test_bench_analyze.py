"""Offline tests of the bench analysis (bench/analyze.py) on synthetic recorded runs: the paired statistics and
intervals, the verdict gate (load, identity change, changed results, missing steps) and the release comparison."""

import contextlib
import copy
import io
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from bench_fixtures import base_wall, fake_panel, synthetic_run, tail_keys  # noqa: E402
from fake_mcp_bench import fake_answer  # noqa: E402
import analyze  # noqa: E402
import bench_lib as bl  # noqa: E402
import suite_common as sc  # noqa: E402


class AnalyseRuns(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.d = self.tmp.name
        self.panel = fake_panel(os.path.join(self.d, "panel.json"))

    def tearDown(self):
        self.tmp.cleanup()

    def test_clean_run_and_paired_ratios(self):
        d = synthetic_run(self.d, "a", self.panel, base_wall)
        res, _ = analyze.analyse(d, boot=300)
        self.assertTrue(res["verdict"]["emitted"], res["verdict"])
        agg = res["aggregates"]["all"]
        self.assertAlmostEqual(agg["pairs"]["mcp/direct-local"]["geo_mean_ratio"], 10.0, places=6)
        self.assertAlmostEqual(agg["pairs"]["direct-remote/direct-local"]["geo_mean_ratio"], 2.0, places=6)
        lo, hi = res["intervals"]["pairs"]["mcp/direct-local"]["geo_mean_ratio"]
        self.assertAlmostEqual(lo, 10.0, places=6)
        self.assertAlmostEqual(hi, 10.0, places=6)
        q = res["queries"]["P-q-0001"]["arms"]["direct-local"]
        self.assertEqual(q["n_success"], 6)
        self.assertEqual(len(q["walls_ms"]), 6)
        self.assertEqual(set(res["order_effect"]), set(bl.ARMS))
        md = analyze.report_md(res)
        self.assertIn("Performance verdict: emitted", md)
        self.assertFalse(tail_keys(res))

    def test_refusals(self):
        loaded = synthetic_run(self.d, "loaded", self.panel, base_wall, sra=7000.0)
        res, _ = analyze.analyse(loaded, boot=0)
        self.assertFalse(res["verdict"]["emitted"])
        self.assertIn("load guard: loaded", res["verdict"]["refused"][0])
        moved = synthetic_run(self.d, "moved", self.panel, base_wall, instances=("i1", "i1", "i2"))
        res, _ = analyze.analyse(moved, boot=0)
        self.assertTrue(any("deployment changed" in r for r in res["verdict"]["refused"]))
        # i1 i1 i2 i1 i1 i2 i1 at the ends of blocks 0..6: the instance changed within blocks 2, 3, 5 and 6
        self.assertEqual(res["identity"]["blocks_changed"], [2, 3, 5, 6])
        # a different result: another workload, never a success
        other = {"patterns": [{"id": "w0", "counts": {"contexts": {"value": 7}}}]}
        changed = synthetic_run(self.d, "changed", self.panel, base_wall,
                                sizes_override={("P-q-0000", "mcp", 3): other})
        res, _ = analyze.analyse(changed, boot=0)
        self.assertEqual(res["classes"]["mcp"].get("result_changed"), 1)
        self.assertEqual(res["queries"]["P-q-0000"]["arms"]["mcp"]["others"][0]["class"], "result_changed")
        self.assertTrue(res["verdict"]["emitted"])  # 5 of 6 still meet --min-success 5
        res, _ = analyze.analyse(changed, boot=0, min_success=6)
        self.assertFalse(res["verdict"]["emitted"])
        # a missing step: incomplete
        incomplete = synthetic_run(self.d, "incomplete", self.panel, base_wall)
        path = os.path.join(incomplete, "arms", "mcp", "executions.jsonl")
        rows = sc.read_jsonl(path)[0]
        os.replace(path, path + ".bak")
        for r in rows[1:]:
            sc.append_jsonl(path, r)
        res, _ = analyze.analyse(incomplete, boot=0)
        self.assertTrue(any(r.startswith("incomplete") for r in res["verdict"]["refused"]))
        self.assertEqual(res["classes"]["mcp"]["missing"], 1)

    def test_an_answer_stating_no_recorded_size_is_unverified(self):
        # D1: every MCP answer replaced by a body that states none of the panel's size fields
        moved = synthetic_run(self.d, "moved_body", self.panel, base_wall,
                              sizes_override=dict(((q["id"], "mcp", c), {"something": "else"})
                                                  for q in self.panel["queries"] for c in range(1, 7)))
        res, _ = analyze.analyse(moved, boot=0)
        self.assertEqual(res["classes"]["mcp"], {"unverified": 36})
        self.assertFalse(res["verdict"]["emitted"])
        others = res["queries"]["P-q-0000"]["arms"]["mcp"]["others"]
        self.assertIn("does not state", others[0]["detail"])

    def test_comparison(self):
        a = synthetic_run(self.d, "a", self.panel, base_wall)
        same = synthetic_run(self.d, "same", self.panel, base_wall, seed=2)
        slow = synthetic_run(self.d, "slow", self.panel,
                             lambda q, arm, c: base_wall(q, arm, c) * (1.5 if arm == "mcp" else 1.0), seed=3)
        cmp_ = analyze.compare(same, a, boot=200, boot_seed=1)
        self.assertEqual(cmp_["verdict"], "pass", cmp_)
        self.assertTrue(all(v["verdict"] == "no_regression" for v in cmp_["arms"].values()))
        cmp_ = analyze.compare(slow, a, boot=200, boot_seed=1)
        self.assertEqual(cmp_["verdict"], "regression")
        self.assertEqual(cmp_["arms"]["mcp"]["verdict"], "regression")
        self.assertAlmostEqual(cmp_["arms"]["mcp"]["geo_mean_ratio"], 1.5, places=6)
        self.assertEqual(cmp_["arms"]["direct-local"]["verdict"], "no_regression")
        fast = synthetic_run(self.d, "fast", self.panel,
                             lambda q, arm, c: base_wall(q, arm, c) * (0.5 if arm == "direct-local" else 1.0), seed=4)
        self.assertEqual(analyze.compare(fast, a, boot=200, boot_seed=1)["arms"]["direct-local"]["verdict"],
                         "improvement")
        loaded = synthetic_run(self.d, "loaded", self.panel, base_wall, sra=9000.0, seed=5)
        cmp_ = analyze.compare(loaded, a, boot=50, boot_seed=1)
        self.assertEqual(cmp_["verdict"], "refused")
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(analyze.main([slow, "--baseline", a, "--boot", "100"]), analyze.EXIT_REGRESSION)
        with open(os.path.join(slow, "comparison.md")) as fh:
            self.assertIn("regression", fh.read())

    def test_the_panel_identity_binds_the_recorded_sizes(self):
        # D2: the same queries with other recorded result sizes (a rebuild on another baseline) are another panel
        panel_b = copy.deepcopy(self.panel)
        for q in panel_b["queries"]:
            if q["mode"] == "pattern":
                for arm in bl.ARMS:
                    q["expected"][arm]["sizes"]["contexts"] = 7
        panel_b["panel_sha256"] = bl.panel_identity(panel_b["queries"])
        self.assertNotEqual(panel_b["panel_sha256"], self.panel["panel_sha256"])
        other = dict(((q["id"], arm, c), {"patterns": [{"id": "w", "counts": {"contexts": {"value": 7}}}]})
                     for q in panel_b["queries"] if q["mode"] == "pattern" for arm in bl.ARMS for c in range(7))
        a = synthetic_run(self.d, "a", self.panel, base_wall)
        b = synthetic_run(self.d, "b", panel_b, base_wall, seed=2, sizes_override=other)
        self.assertTrue(analyze.analyse(b, boot=0)[0]["verdict"]["emitted"])
        cmp_ = analyze.compare(b, a, boot=50, boot_seed=1)
        self.assertEqual(cmp_["verdict"], "refused")
        self.assertIn("the runs measured different panels", cmp_["refused"])

    def test_comparison_refuses_differing_successful_sizes(self):
        # D2: both runs match the panel's recorded sizes, but the candidate's answers state another result size the
        # panel did not record (returned 5): not the same result, no comparison
        extra = dict(((q["id"], "direct-local", c), {"patterns": [{"id": "w", "returned": 5,
                                                                  "counts": {"contexts": {"value": 1}}}]})
                     for q in self.panel["queries"] if q["mode"] == "pattern" for c in range(7))
        a = synthetic_run(self.d, "a", self.panel, base_wall)
        b = synthetic_run(self.d, "b", self.panel, base_wall, seed=2, sizes_override=extra)
        self.assertEqual(analyze.analyse(b, boot=0)[0]["classes"]["direct-local"], {"success": 36})
        cmp_ = analyze.compare(b, a, boot=50, boot_seed=1)
        self.assertEqual(cmp_["verdict"], "refused")
        self.assertTrue(any(r.startswith("result sizes differ between the runs") for r in cmp_["refused"]),
                        cmp_["refused"])
        self.assertIn("P-q-0000 direct-local", " ".join(cmp_["refused"]))

    def test_thresholds_may_only_tighten_after_the_run(self):
        # D3: a loaded run stays without a verdict whatever limit is passed afterwards; a tighter one may refuse more
        loaded = synthetic_run(self.d, "loaded", self.panel, base_wall, sra=7700.0)
        res, _ = analyze.analyse(loaded, boot=0, guard_overrides={"sra_all_cpu_pct_max": 10000.0})
        self.assertFalse(res["verdict"]["emitted"])
        self.assertTrue(any(r.startswith("thresholds loosened after the run") for r in res["verdict"]["refused"]),
                        res["verdict"]["refused"])
        for loosen in ({"host_load1_max": 99.0}, {"max_gap_s": 60.0}, {"pad_s": 0.0}, {"allow_missing_sra": True}):
            clean = synthetic_run(self.d, "clean-%s" % sorted(loosen)[0], self.panel, base_wall)
            res, _ = analyze.analyse(clean, boot=0, guard_overrides=loosen)
            self.assertFalse(res["verdict"]["emitted"], loosen)
        clean = synthetic_run(self.d, "tight", self.panel, base_wall)
        res, _ = analyze.analyse(clean, boot=0, guard_overrides={"host_load1_max": 5.0, "sra_all_cpu_pct_max": 50.0})
        self.assertTrue(res["verdict"]["emitted"], res["verdict"]["refused"])
        res, _ = analyze.analyse(clean, boot=0, guard_overrides={"host_load1_max": 0.5})
        self.assertEqual(res["load_guard"]["verdict"], "loaded")
        # the comparison records the overrides and refuses a loosened one
        a = synthetic_run(self.d, "a", self.panel, base_wall)
        cmp_ = analyze.compare(loaded, a, boot=50, boot_seed=1, guard_overrides={"sra_all_cpu_pct_max": 10000.0})
        self.assertEqual(cmp_["verdict"], "refused")
        self.assertEqual(cmp_["guard"]["overrides"], {"sra_all_cpu_pct_max": 10000.0})
        self.assertEqual(cmp_["guard"]["candidate"]["sra_all_cpu_pct_max"], 10000.0)
        self.assertTrue(cmp_["guard"]["candidate_loosened"])

    def test_an_unknown_identity_refuses_the_verdict(self):
        # D4: no evidence that the deployment stayed the same is no evidence of it
        blind = synthetic_run(self.d, "blind", self.panel, base_wall)
        sc.write_json_atomic(os.path.join(blind, "arms", "direct-local", "run.json"),
                             {"sessions": [{"started_utc": "x"}]})
        os.replace(os.path.join(blind, "steps.jsonl"), os.path.join(blind, "steps.jsonl.bak"))
        for blk in range(7):
            sc.append_jsonl(os.path.join(blind, "steps.jsonl"), {"probe": "identity", "block": blk, "rc": 255,
                                                                 "release": None, "server_instance": None})
        res, _ = analyze.analyse(blind, boot=0)
        self.assertFalse(res["verdict"]["emitted"])
        why = " ".join(res["verdict"]["refused"])
        self.assertIn("identity is unknown", why)
        self.assertIn("blocks 0, 1, 2, 3, 4, 5, 6", why)
        self.assertIn("direct-local states no start identity", why)
        # the MCP service's identity unknown for one invocation
        res, _ = analyze.analyse(synthetic_run(self.d, "mcp-blind", self.panel, base_wall, mcp_caps=None), boot=0)
        self.assertIn("state no service capabilities", " ".join(res["verdict"]["refused"]))
        # a block whose probe failed, the others fine
        one = synthetic_run(self.d, "one", self.panel, base_wall)
        rows = sc.read_jsonl(os.path.join(one, "steps.jsonl"))[0]
        os.replace(os.path.join(one, "steps.jsonl"), os.path.join(one, "steps.jsonl.bak"))
        for r in rows:
            if r["block"] == 4:
                r.update(release=None, server_instance=None, rc=255)
            sc.append_jsonl(os.path.join(one, "steps.jsonl"), r)
        res, _ = analyze.analyse(one, boot=0)
        self.assertIn("end of blocks 4", " ".join(res["verdict"]["refused"]))

    def test_a_changed_index_refuses_the_verdict(self):
        # D4: the index fingerprint at a block's end is compared too
        moved = synthetic_run(self.d, "fp", self.panel, base_wall, index_fps=("fp1", "fp1", "fp1", "fp2"))
        res, _ = analyze.analyse(moved, boot=0)
        self.assertFalse(res["verdict"]["emitted"])
        self.assertIn("index_fp: fp1, fp2", " ".join(res["verdict"]["refused"]))
        self.assertEqual(res["identity"]["blocks_changed"], [3, 4])  # fp1 fp1 fp1 fp2 fp1 fp1 fp1

    def test_an_accepted_cap_and_nothing_else(self):
        # D5: a panel query kept with a stated cap succeeds only with exactly that cap
        panel = copy.deepcopy(self.panel)
        why = ["result 0: label evidence lower_bound"]
        q = [x for x in panel["queries"] if x["id"] == "T-q-0000"][0]
        q["accepted_incomplete"] = {"reason": "test", "arms": dict((arm, list(why)) for arm in bl.ARMS)}
        panel["panel_sha256"] = bl.panel_identity(panel["queries"])
        capped = fake_answer(q["route"], q["scenario"]["request"])
        capped["results"][0]["outcome"] = {"label_evidence": "lower_bound"}
        stopped = copy.deepcopy(capped)
        stopped["results"][0]["arms"]["right"]["status"] = "truncated"
        bodies = dict(((q["id"], arm, c), capped) for arm in bl.ARMS for c in range(7))
        bodies[(q["id"], "mcp", 2)] = stopped
        bodies[(q["id"], "direct-remote", 3)] = fake_answer(q["route"], q["scenario"]["request"])  # not capped
        res, _ = analyze.analyse(synthetic_run(self.d, "capped", panel, base_wall, sizes_override=bodies), boot=0)
        arms = res["queries"]["T-q-0000"]["arms"]
        self.assertEqual(arms["direct-local"]["classes"], {"success": 6})
        self.assertEqual(arms["mcp"]["classes"], {"success": 5, "incomplete": 1})
        self.assertEqual(arms["direct-remote"]["classes"], {"success": 5, "result_changed": 1})
        # without the exception every capped answer is incomplete
        res, _ = analyze.analyse(synthetic_run(self.d, "capped2", self.panel, base_wall, sizes_override=bodies), boot=0)
        self.assertEqual(res["queries"]["T-q-0000"]["arms"]["direct-local"]["classes"], {"incomplete": 6})

    def test_comparison_needs_one_harness(self):
        # D8: the same harness on both sides: a pause, warm-up arm, MCP template, base URL or runner change moves
        # the walls by itself
        a = synthetic_run(self.d, "a", self.panel, base_wall)
        for i, extra in enumerate(({"pause_s": 0.5}, {"warmup_arm": "mcp"}, {"mcp_cmd": "other {ids} {out}"},
                                   {"remote_base": "https://elsewhere/api"}, {"local_base": "http://localhost:1"},
                                   {"checks": "stub"})):
            b = synthetic_run(self.d, "b%d" % i, self.panel, base_wall, seed=2, bench_extra=extra)
            cmp_ = analyze.compare(b, a, boot=20, boot_seed=1)
            self.assertEqual(cmp_["verdict"], "refused", extra)
            key = sorted(extra)[0]
            self.assertTrue(any(r.startswith("the runs' harness differs") and key in r for r in cmp_["refused"]),
                            (extra, cmp_["refused"]))
        # the snapshot's runner and sampler: their sha256 in scripts/SNAPSHOT.json
        b = synthetic_run(self.d, "snap", self.panel, base_wall, seed=3)
        for run, sha in ((a, "aaa"), (b, "bbb")):
            os.makedirs(os.path.join(run, "scripts"), exist_ok=True)
            sc.write_json_atomic(os.path.join(run, "scripts", "SNAPSHOT.json"),
                                 {"files": {"run_direct.py": sha, "suite_common.py": "s", "sampler.py": "p"}})
        cmp_ = analyze.compare(b, a, boot=20, boot_seed=1)
        self.assertTrue(any("run_direct.py" in r for r in cmp_["refused"]), cmp_["refused"])
        self.assertEqual(cmp_["harness"]["differs"], ["snapshot run_direct.py"])

    def test_heterogeneous_change_is_inconclusive(self):
        # half the queries 30 % slower, half 20 % faster: the panel's mean barely moves, its interval spans both
        a = synthetic_run(self.d, "a", self.panel, base_wall)
        mixed = synthetic_run(self.d, "mixed", self.panel,
                              lambda q, arm, c: base_wall(q, arm, c) * (1.3 if int(q[-1]) % 2 else 0.8), seed=6)
        cmp_ = analyze.compare(mixed, a, boot=400, boot_seed=2)
        self.assertEqual(cmp_["verdict"], "inconclusive", cmp_["arms"])
        lo, hi = cmp_["arms"]["direct-local"]["interval"]
        self.assertLess(lo, 1.0)
        self.assertGreater(hi, 1.1)


if __name__ == "__main__":
    unittest.main()
