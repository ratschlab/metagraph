"""Offline tests of the benchmark panel (lane B): bench_lib's result sizes, completeness, run manifest and schedule,
and make_panel.py's eligibility, selection, weights, costing and identity on a synthetic baseline run."""

import contextlib
import io
import json
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.dirname(HERE)
sys.path.insert(0, SUITE)
sys.path.insert(0, os.path.join(SUITE, "bench"))

import bench_lib as bl  # noqa: E402
import make_panel as bp  # noqa: E402
import suite_common as sc  # noqa: E402

ROUTES = {"resolve": "POST /resolve", "traverse": "POST /traverse", "pattern": "POST /pattern"}
TOOLS = {"resolve": "traverse_resolve", "traverse": "traverse", "pattern": "pattern_search"}


def resolve_body(labels=1, kmers=70, runs=1):
    return {"num_kmers": kmers, "candidates": [{"labels": ["1"]}] if runs else [],
            "labels": [{"label": str(i), "kind": "column", "runs": [[0, kmers]] * runs,
                        "kmers_supported": kmers if runs else 0} for i in range(labels)]}


def traverse_body(paths=1, bp=100, labels=1, status="complete"):
    return {"results": [{"arms": {"right": {"status": status, "paths": [{"length_bp": bp}] * paths,
                                            "segments": [{"id": i} for i in range(paths)],
                                            "counters": {"steps": bp, "work_units": 10 * bp}}},
                         "label_dict": [{"name": "L%d" % i} for i in range(labels)],
                         "outcome": {"walks": "complete"}, "annotation": {"rows_requested": bp}}]}


def pattern_body(returned=2, contexts=2, labels=None, scope="any_offset", extra=None):
    p = {"id": "w", "returned": returned, "scope": scope, "counts": {"contexts": {"value": contexts}},
         "work": {"steps": 10}}
    if labels is not None:
        p["counts"]["labels"] = {"value": labels}
        p["by_label"] = [{"column": str(i)} for i in range(labels)]
    p.update(extra or {})
    return {"patterns": [p]}


def scenario(sid, mode, area=None, request=None, **kw):
    s = {"id": sid, "area": area or sid[0], "route": ROUTES[mode], "timing_class": "fast", "title": sid,
         "request": request if request is not None else {"q": sid}, "repeat": 1, "expect": [], "tags": [],
         "truth": {}, "mcp": {"tool": TOOLS[mode], "arguments": {"q": sid}}}
    s.update(kw)
    return s


class Baseline(object):
    """A synthetic merged and compared suite run: manifest.json, executions.jsonl, compare.json, run.json and the
    raw answers, every arm with the same answer unless told otherwise."""

    def __init__(self, root):
        self.dir = os.path.join(root, "20261014T000000Z_refseq97-test")
        os.makedirs(os.path.join(self.dir, "arms"))
        self.scenarios, self.rows, self.compare = [], [], {}

    def add(self, s, body, wall=10.0, outcome="pass", http=200, time_limited="", arms=sc.ARMS, mcp_class="agree",
            direct_class="identical", digests_equal=None, mcp_body=None):
        self.scenarios.append(s)
        entry = {"arms": {}, "direct_agreement": {"class": direct_class}, "mcp_agreement": {"class": mcp_class}}
        for arm in arms:
            raw = "arms/%s/raw/%s.%s.0.json.gz" % (arm, s["id"], arm)
            os.makedirs(os.path.dirname(os.path.join(self.dir, raw)), exist_ok=True)
            b = mcp_body if (arm == "mcp" and mcp_body is not None) else body
            sc.write_gz_json(os.path.join(self.dir, raw), {"answer": {"http": http, "body": b, "error": None}})
            w = wall * (100 if arm == "mcp" else 1)
            self.rows.append({"scenario": s["id"], "arm": arm, "repeat": 0, "outcome": outcome, "http": http,
                              "time_limited": time_limited, "wall_ms": w, "raw": raw,
                              "started_utc": "2026-10-14T06:00:00.000Z", "ended_utc": "2026-10-14T06:00:01.000Z",
                              "server": {"elapsed_ms": wall / 2.0}, "answer_digest": "d-%s" % s["id"],
                              "client": {"bytes_body": len(json.dumps(b))},
                              "mcp_phases": {"backend_ms": wall, "queue_ms": 5.0} if arm == "mcp" else None})
            entry["arms"][arm] = {"digests_equal": digests_equal}
        self.compare[s["id"]] = entry

    def write(self):
        manifest = {"manifest_version": 4, "index": {"database": "test"}, "relations": [], "scenarios": self.scenarios}
        manifest["content_sha256"] = sc.scenarios_sha256(self.scenarios)
        sc.write_json_atomic(os.path.join(self.dir, "manifest.json"), manifest)
        for r in self.rows:
            sc.append_jsonl(os.path.join(self.dir, "executions.jsonl"), r)
        sc.write_json_atomic(os.path.join(self.dir, "compare.json"), {
            "manifest_sha256": manifest["content_sha256"], "run_id": os.path.basename(self.dir),
            "scenarios": self.compare})
        sc.write_json_atomic(os.path.join(self.dir, "run.json"), {"run_id": os.path.basename(self.dir),
                                                                  "release": "refseq97-test", "index_fp": "fp"})
        return manifest


def build_mixed(root, n=10):
    """A baseline with n eligible queries per mode at costs 1..n ms x 10^k, plus one of each kind of ineligible."""
    b = Baseline(root)
    for i in range(n):
        wall = 2.0 ** i
        b.add(scenario("R-pos-%04d" % i, "resolve"),
              resolve_body(labels=6 if i in (2, 7) else 1, kmers=970 if i == 5 else 70), wall=wall)
        b.add(scenario("T-walk-%04d" % i, "traverse"),
              traverse_body(bp=1500 if i % 3 == 0 else 100, labels=8 if i == 4 else 1), wall=wall)
        b.add(scenario("P-exact-%04d" % i, "pattern", area="P" if i % 2 else "M"),
              pattern_body(returned=0 if i == 1 else 3, contexts=0 if i == 1 else 3,
                           scope="long" if i in (3, 6) else "any_offset"), wall=wall)
    return b


class ResultSizes(unittest.TestCase):
    def test_sizes_per_mode(self):
        self.assertEqual(bl.result_sizes("POST /resolve", resolve_body(labels=3, kmers=70)),
                         {"labels": 3, "num_kmers": 70, "runs": 3, "supported_kmers": 210, "candidates": 1})
        t = bl.result_sizes("POST /traverse", traverse_body(paths=2, bp=500, labels=4))
        self.assertEqual((t["paths"], t["path_bp"], t["labels"], t["segments"]), (2, 1000, 4, 2))
        p = bl.result_sizes("POST /pattern", pattern_body(returned=5, contexts=9, labels=6))
        self.assertEqual((p["returned"], p["contexts"], p["labels"], p["by_label"]), (5, 9, 6, 6))
        self.assertEqual(bl.result_sizes("GET /capabilities", {"x": 1}), {})

    def test_sizes_equal_reads_only_fields_both_state(self):
        ok, diff = bl.sizes_equal({"paths": 2, "labels": 4}, {"paths": 2, "labels": None})
        self.assertTrue(ok)
        ok, diff = bl.sizes_equal({"paths": 2, "path_bp": 100}, {"paths": 3, "path_bp": 100})
        self.assertEqual((ok, diff), (False, ["paths"]))

    def test_incomplete_answers(self):
        self.assertEqual(bl.answer_incomplete("POST /pattern", pattern_body()), [])
        self.assertTrue(bl.answer_incomplete("POST /pattern", pattern_body(
            extra={"error": {"code": "information_below_floor"}})))
        self.assertTrue(bl.answer_incomplete("POST /pattern", pattern_body(extra={"withheld": {"reason": "deadline"}})))
        self.assertTrue(bl.answer_incomplete("POST /pattern", pattern_body(extra={"stop": {"reason": "max_steps"}})))
        self.assertTrue(bl.answer_incomplete("POST /pattern", pattern_body(extra={"cut": {"reason": "max_contexts"}})))
        self.assertTrue(bl.answer_incomplete("POST /traverse", traverse_body(status="truncated")))
        # a walk that reached its requested scope (walks partial, arms complete) is complete
        body = traverse_body()
        body["results"][0]["outcome"]["walks"] = "partial"
        self.assertEqual(bl.answer_incomplete("POST /traverse", body), [])
        body = resolve_body()
        body["labels_truncated"] = {"total": 99}
        self.assertTrue(bl.answer_incomplete("POST /resolve", body))
        self.assertEqual(bl.answer_incomplete("POST /resolve", resolve_body()), [])

    def test_tags(self):
        self.assertIn("label_heavy", bl.query_tags("POST /resolve", resolve_body(labels=5)))
        self.assertIn("long_query", bl.query_tags("POST /resolve", resolve_body(kmers=970)))
        self.assertIn("empty_result", bl.query_tags("POST /resolve", resolve_body(runs=0)))
        self.assertIn("long_path", bl.query_tags("POST /traverse", traverse_body(bp=1000)))
        self.assertIn("long_path", bl.query_tags("POST /pattern", pattern_body(scope="long")))
        self.assertIn("large_answer", bl.query_tags("POST /pattern", pattern_body(), bytes_body=70000))
        self.assertNotIn("empty_result", bl.query_tags("POST /traverse", {"results": [{"graphlet_lines": 30}]}))


class Schedule(unittest.TestCase):
    def test_latin_square_balance(self):
        ids = ["Q%02d" % i for i in range(48)]
        steps = bl.schedule(ids, seed=7)
        self.assertEqual(len(steps), bl.requests_planned(48))
        self.assertEqual(len(steps), 912)
        warm = [s for s in steps if s["block"] == 0]
        self.assertEqual(len(warm), 48)
        self.assertTrue(all(s["arm"] == "direct-local" and not s["measured"] for s in warm))
        self.assertEqual(steps[:48], warm)  # the warm-up pass first
        bal = bl.schedule_balance(steps)
        for q, orders in bal["per_query"].items():
            self.assertEqual(sorted(orders), list(range(6)), q)  # every order once per query
            self.assertEqual(set(orders.values()), {1})
        for c, orders in bal["per_cycle"].items():
            self.assertEqual(sorted(orders.items()), [(o, 8) for o in range(6)], c)  # 8 of each per cycle
        self.assertEqual(bal["first_arm"], {"direct-local": 96, "direct-remote": 96, "mcp": 96})
        # a triple runs back to back, in its order, within its block
        for i in range(48, len(steps), 3):
            triple = steps[i:i + 3]
            self.assertEqual(len(set((s["query"], s["cycle"], s["order"]) for s in triple)), 1)
            self.assertEqual([s["arm"] for s in triple], bl.arm_orders()[triple[0]["order"]])
            self.assertEqual([s["position"] for s in triple], [0, 1, 2])
        self.assertEqual([s["block"] for s in steps], sorted(s["block"] for s in steps))
        self.assertEqual([s["step"] for s in steps], list(range(len(steps))))

    def test_seed_reproduces(self):
        ids = ["Q%02d" % i for i in range(12)]
        self.assertEqual(bl.schedule(ids, seed=3), bl.schedule(ids, seed=3))
        self.assertNotEqual(bl.schedule(ids, seed=3), bl.schedule(ids, seed=4))

    def test_bench_ids(self):
        self.assertEqual(bl.bench_id("R-pos-0001", "w0"), "R-pos-0001~w0")
        self.assertEqual(bl.bench_id("R-pos-0001", 3), "R-pos-0001~c3")
        self.assertEqual(bl.parse_bench_id("R-pos-0001~c3"), ("R-pos-0001", "c3", 3))
        self.assertEqual(bl.parse_bench_id("R-pos-0001~w0"), ("R-pos-0001", "w0", 0))
        self.assertIsNone(bl.parse_bench_id("R-pos-0001"))
        with self.assertRaises(ValueError):
            bl.bench_id("R-pos-0001", 0)


class BuildPanel(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = self.tmp.name

    def tearDown(self):
        self.tmp.cleanup()

    def build(self, b, **kw):
        manifest = b.write()
        args = dict(size=12, quota={"resolve": 4, "traverse": 4, "pattern": 4},
                    min_tags={"label_heavy": 1, "long_path": 1, "long_query": 1, "heavy": 1})
        args.update(kw)
        return bp.build(manifest, bp.load_baseline(b.dir), now=0, **args)

    def test_eligibility_rules(self):
        b = build_mixed(self.root)
        b.add(scenario("C-caps-0001", "resolve", route="GET /capabilities"), resolve_body())
        b.add(scenario("C-protocol-0001", "pattern", body_raw="{"), pattern_body())
        b.add(scenario("R-long-0001", "resolve", timing_class="long"), resolve_body())
        b.add(scenario("R-nim-0001", "resolve", mcp={"not_implemented": "x"}), resolve_body())
        b.add(scenario("R-miss-0001", "resolve"), resolve_body(), arms=("direct-local", "direct-remote"))
        b.add(scenario("R-fail-0001", "resolve"), resolve_body(), outcome="fail")
        b.add(scenario("R-refused-0001", "resolve"), resolve_body(), http=400)
        b.add(scenario("R-timed-0001", "resolve"), resolve_body(), time_limited="stop time")
        b.add(scenario("R-disagree-0001", "resolve"), resolve_body(), mcp_class="disagree")
        b.add(scenario("R-flaky-0001", "resolve"), resolve_body(), digests_equal=False)
        b.add(scenario("P-floor-0001", "pattern"), pattern_body(extra={"error": {"code": "information_below_floor"}}))
        # the D area's copy of an eligible request is a duplicate, the original kept
        b.add(scenario("D-R-pos-0001", "resolve", area="D", request={"q": "R-pos-0003"},
                       mcp={"tool": "traverse_resolve", "arguments": {"q": "R-pos-0003"}}), resolve_body())
        panel = self.build(b)
        ex = panel["selection"]["excluded_ids"]
        self.assertEqual(ex["route"], ["C-caps-0001"])
        self.assertEqual(ex["protocol"], ["C-protocol-0001"])
        self.assertEqual(ex["long"], ["R-long-0001"])
        self.assertEqual(ex["mcp"], ["R-nim-0001"])
        self.assertEqual(ex["missing"], ["R-miss-0001"])
        self.assertEqual(ex["outcome"], ["R-fail-0001"])
        self.assertEqual(ex["http"], ["R-refused-0001"])
        self.assertEqual(ex["time"], ["R-timed-0001"])
        self.assertEqual(sorted(ex["agreement"]), ["R-disagree-0001", "R-flaky-0001"])
        self.assertEqual(ex["complete"], ["P-floor-0001"])
        self.assertEqual(ex["duplicate"], ["D-R-pos-0001"])
        self.assertEqual(panel["selection"]["eligible"], 30)

    def test_selection_quota_tags_and_identity(self):
        b = build_mixed(self.root)
        panel = self.build(b)
        q = panel["queries"]
        self.assertEqual(len(q), 12)
        self.assertEqual(panel["selection"]["by_mode"], {"resolve": 4, "traverse": 4, "pattern": 4})
        modes = dict((m, [x for x in q if x["mode"] == m]) for m in bl.MODES)
        self.assertTrue(any("label_heavy" in x["tags"] for x in modes["resolve"]))
        self.assertTrue(any("long_query" in x["tags"] for x in modes["resolve"]))
        self.assertTrue(any("long_path" in x["tags"] for x in modes["traverse"]))
        self.assertTrue(any("label_heavy" in x["tags"] for x in modes["traverse"]))
        self.assertTrue(any("long_path" in x["tags"] for x in modes["pattern"]))
        self.assertTrue(all(any("heavy" in x["tags"] for x in v) for v in modes.values()))
        # the pattern route's two areas both chosen; the cost range spanned (the cheapest and dearest bins)
        self.assertEqual(set(x["area"] for x in modes["pattern"]), {"M", "P"})
        self.assertTrue(all(x["weight"] == 1.0 / 12 for x in q))
        self.assertEqual(panel["protocol"]["requests"], 12 * 19)
        self.assertEqual(panel["costing"]["requests"], 228)
        # every query carries its frozen scenario and each arm's recorded sizes
        for x in q:
            self.assertEqual(x["scenario_sha256"], bl.scenario_sha256(x["scenario"]))
            self.assertEqual(sorted(x["expected"]), sorted(sc.ARMS))
            self.assertTrue(x["expected"]["direct-local"]["sizes"])
        path = os.path.join(self.root, "panel.json")
        sc.write_json_atomic(path, panel)
        self.assertEqual(bl.read_panel(path)["panel_sha256"], panel["panel_sha256"])
        # deterministic: a second build chooses the same queries
        again = self.build(b)
        self.assertEqual(again["panel_sha256"], panel["panel_sha256"])
        # a tampered scenario is refused
        panel["queries"][0]["scenario"]["request"] = {"q": "other"}
        sc.write_json_atomic(path, panel)
        with self.assertRaises(ValueError):
            bl.read_panel(path)

    def test_empty_cap_include_exclude(self):
        b = build_mixed(self.root)
        for i in range(6):
            b.add(scenario("N-zero-%04d" % i, "pattern", area="N"), pattern_body(returned=0, contexts=0),
                  wall=0.5 + i)
        panel = self.build(b, max_empty=1, include=["P-exact-0008"], exclude=["P-exact-0009"])
        pat = [x for x in panel["queries"] if x["mode"] == "pattern"]
        self.assertLessEqual(sum(1 for x in pat if "empty_result" in x["tags"]), 1)
        self.assertIn("P-exact-0008", [x["id"] for x in pat])
        self.assertNotIn("P-exact-0009", [x["id"] for x in panel["queries"]])
        self.assertEqual([x["selected_by"] for x in pat if x["id"] == "P-exact-0008"], ["include"])

    def test_run_manifest_verifies(self):
        b = build_mixed(self.root)
        panel = self.build(b)
        m = bl.run_manifest(panel)
        self.assertEqual(len(m["scenarios"]), 12 * 7)
        ok, _ = sc.verify_manifest(m)
        self.assertTrue(ok)
        ids = [s["id"] for s in m["scenarios"]]
        q0 = panel["queries"][0]["id"]
        self.assertEqual(ids[:7], [q0 + "~w0"] + ["%s~c%d" % (q0, c) for c in range(1, 7)])
        self.assertTrue(all(s["repeat"] == 1 for s in m["scenarios"]))
        self.assertEqual(m["bench"]["panel_sha256"], panel["panel_sha256"])
        # the bench scenario is the query's own request
        self.assertEqual(m["scenarios"][3]["request"], panel["queries"][0]["scenario"]["request"])

    def test_cli(self):
        b = build_mixed(self.root)
        b.write()
        out = os.path.join(self.root, "p.json")
        with contextlib.redirect_stdout(io.StringIO()) as said:
            rc = bp.main(["--baseline", b.dir, "--out", out, "--size", "12", "--quota",
                          "resolve=4,traverse=4,pattern=4", "--min-tags", "label_heavy=1,long_path=1"])
        self.assertEqual(rc, 0)
        self.assertIn("12 queries (resolve 4, traverse 4, pattern 4)", said.getvalue())
        self.assertEqual(len(bl.read_panel(out)["queries"]), 12)
        with self.assertRaises(SystemExit):
            bp.main(["--baseline", b.dir, "--out", out, "--size", "10", "--quota", "resolve=4,traverse=4,pattern=4"])


if __name__ == "__main__":
    unittest.main()
