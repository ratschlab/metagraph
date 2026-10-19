"""Offline tests of area G's scenarios (manifest v8, lane D): the candidates drawn from v7's T seeds, the scenarios
gen_G makes (the exact merged strategy, its MCP mapping, the full twin, the ways, the owner's not_exposed reasons),
the op plans (every argument the service takes named, library-only tools and arguments without W3, selectors that
read earlier ops, one deliberate cut per kind), a corpus file's choices, and the op plans run by the local runner on
the library's real fixture graphlets (every op answers what area_g.expected_outcome says; every cut op shows a cut).
The v8 build against the frozen v7 skips without the record cache (SUITE_TRUTH_CACHE or truth/cache)."""
import copy
import gzip
import inspect
import json
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.dirname(HERE)
LIB = os.path.normpath(os.path.join(SUITE, "..", "..", "api", "python"))
DATA = os.path.join(LIB, "tests", "data", "traverse")
sys.path.insert(0, SUITE)

import area_g  # noqa: E402
import manifest_lib as ml  # noqa: E402
import suite_common as sc  # noqa: E402
import truth  # noqa: E402

needs_library = unittest.skipIf(sys.version_info < (3, 10), "the traversal library needs Python 3.10 or later")
V7 = os.path.join(SUITE, "manifest", "suite_v7.json")


def have_cache():
    return os.path.isdir(truth.CACHE_DIR) and bool(os.listdir(truth.CACHE_DIR))


class FakeCtx(object):
    """The builder's context without the record cache: no record (the drawn long-walk candidates skip), the skips
    recorded."""

    def __init__(self):
        self.skipped = []

    def rec(self, key):
        raise KeyError(key)

    def rng(self, *parts):
        import random
        return random.Random(0)

    def skip(self, area, what, reason):
        self.skipped.append({"area": area, "scenario": what, "reason": reason})


_V7 = {}


def v7_scenarios():
    if "s" not in _V7:
        _V7["s"] = sc.read_json(V7)["scenarios"]
    return _V7["s"]


def provisional(**kw):
    ctx = FakeCtx()
    return area_g.gen_G(ctx, v7_scenarios(), **kw), ctx


class Candidates(unittest.TestCase):
    def test_every_slot_has_candidates_from_v7_seeds(self):
        cands = area_g.candidates(FakeCtx(), v7_scenarios())
        slots = set(c["slot"] for c in cands)
        # the long walk is drawn from the records (none without the cache); two databases is the rehearsal's
        self.assertEqual(slots, set(area_g.SLOT_NAMES) - {"long_walk", "two_databases"})
        ids = [c["cid"] for c in cands]
        self.assertEqual(len(ids), len(set(ids)))
        by_id = dict((s["id"], s) for s in v7_scenarios())
        for c in cands:
            src = by_id[c["source"].split(" ")[0]]
            self.assertEqual(c["seed"]["sequence"], src["request"]["seeds"][0]["sequence"], c["cid"])
            self.assertIn(c["preset"], area_g.PRESETS)
            # every candidate states v7's evidence and its expected shape
            self.assertTrue(c["evidence"] and isinstance(c["expected"], dict), c["cid"])

    def test_direct_only_candidates_are_the_inexpressible_strategies(self):
        cands = area_g.candidates(FakeCtx(), v7_scenarios())
        direct = sorted(c["cid"] for c in cands if not area_g.expressible(c))
        self.assertEqual(direct, ["annotate.1", "annotate.2", "coordinates.1", "coordinates.2", "no_bases.1",
                                  "no_bases.2"])
        for c in cands:
            if area_g.expressible(c):
                req = {"seeds": [c["seed"]], "strategy": area_g.candidate_strategy(c)}
                mcp = ml.mcp_traverse(req)
                self.assertIn("tool", mcp, c["cid"])
                self.assertEqual(area_g.mcp_merge(mcp["arguments"]["preset"], mcp["arguments"].get("overrides")),
                                 req["strategy"], c["cid"])
            else:
                self.assertIn("not_implemented", ml.mcp_traverse({"seeds": [c["seed"]],
                                                                 "strategy": area_g.candidate_strategy(c)}))

    def test_the_merged_strategy_is_the_v7_echo(self):
        # the echo of T-graphlet-0001 (v7 run 20261019T111815Z, traverse_status.overrides)
        echo = {"bounds": {"max_memory_mb": 512, "max_output_bp": 500000, "max_work_units": 100000000,
                           "time_budget_ms": 15000, "max_extension_bp": 500},
                "labels": {"mode": "constrain"}, "output": {"detail": "graphlet", "timing": True},
                "branching": {"max_label_branches": 2}, "direction": "both", "exhaustive": False}
        self.assertEqual(area_g.mcp_merge("walk", {"bounds": {"max_extension_bp": 500}}), echo)


class Scenarios(unittest.TestCase):
    def test_provisional_corpus(self):
        g, ctx = provisional()
        ids = [s["id"] for s in g]
        self.assertEqual(ids, ["G-arms2-0001", "G-pair-radius-0001", "G-pair-seeds-0001", "G-merge-0001",
                               "G-pair-cap-0001", "G-arm1-0001", "G-cut-0001", "G-stop-0001", "G-depth0-0001",
                               "G-linear-0001", "G-fail-0001", "G-trace-0001", "G-coords-0001", "G-annot-0001",
                               "G-nobases-0001"])
        skipped = dict((x["scenario"], x["reason"]) for x in ctx.skipped)
        for sid in ("G-big-0001", "G-twodb-0001", "G-long-0001"):
            self.assertIn(sid, skipped)
        self.assertTrue(skipped["G-big-0001"].startswith("slot over_8mib: v9"))
        problems, pending = area_g.validate_g(g, set())
        self.assertEqual(problems, [])
        self.assertTrue(any("preconditions" in p for p in pending))
        self.assertTrue(any("library_commit" in p for p in pending))
        self.assertTrue(any(area_g.SERVICE_CAP_CASE in p for p in pending))
        direct = ("G-coords-0001", "G-annot-0001", "G-nobases-0001")
        for s in g:
            st = s["request"]["strategy"]
            if s["id"] not in direct:
                self.assertEqual(st["output"], {"detail": "graphlet", "timing": True}, s["id"])
            self.assertEqual(st["output"]["detail"], "graphlet")
            self.assertEqual(len(s["request"]["seeds"]), 1)
            twin = copy.deepcopy(st)
            twin["output"]["detail"] = "full"
            del twin["bounds"]["max_memory_mb"]     # F-D4: the twin without the memory budget
            self.assertEqual(s["graphlet"]["twin"]["strategy"], twin)
            self.assertEqual(st["bounds"]["max_memory_mb"], 512)
            if s["id"] in direct:
                self.assertIn("not_implemented", s["mcp"])
                self.assertEqual(s["graphlet"]["ways"]["default"], area_g.DIRECT_ONLY_WAYS)
            else:
                self.assertEqual(s["mcp"]["tool"], "traverse_fetch")
                self.assertEqual(s["graphlet"]["ways"]["default"], area_g.ALL_WAYS)
            self.assertEqual(s["area"], "G")
            if s["id"] != "G-fail-0001":
                self.assertIn("rehearsal-pending", s["tags"])
        fail = dict((s["id"], s) for s in g)["G-fail-0001"]
        self.assertEqual(fail["graphlet"]["preconditions"], {"no_graphlet": True})
        self.assertEqual([o["tool"] for o in fail["graphlet"]["ops"]],
                         ["graphlet_summary", "graphlet_walks", "graphlet_walk", "graphlet_labels",
                          "graphlet_sequence", "graphlet_export", "graphlet_compare"])

    def test_the_builder_validates_area_g(self):
        import build_manifest as bm
        g, _ = provisional()
        self.assertEqual([p for p in bm.validate(g) if not p.startswith("G-")], [])
        self.assertEqual(bm.validate(g), [])

    def test_direct_only_scenarios(self):
        g, _ = provisional(include_direct_only=True)
        by_id = dict((s["id"], s) for s in g)
        for sid in ("G-annot-0001", "G-nobases-0001"):
            s = by_id[sid]
            self.assertIn("not_implemented", s["mcp"])
            self.assertIn("direct-only", s["tags"])
            self.assertEqual(s["graphlet"]["ways"]["default"], area_g.DIRECT_ONLY_WAYS)
            for o in s["graphlet"]["ops"]:
                self.assertNotIn("W3", o.get("ways") or [], (sid, o["id"]))
                self.assertNotIn("W4", o.get("ways") or [], (sid, o["id"]))
                self.assertEqual(o["not_exposed"], {"W3": area_g.REASON_DIRECT_ONLY, "W4": area_g.REASON_DIRECT_ONLY})
        self.assertEqual(area_g.REASON_DIRECT_ONLY, "not expressible via MCP on staging (advanced strategy path off)")
        coords = by_id["G-coords-0001"]
        self.assertTrue(coords["request"]["strategy"]["output"]["coordinates"])
        self.assertEqual(coords["request"]["strategy"]["support"], "trace")
        self.assertIn("labels-name-coords", [o["id"] for o in coords["graphlet"]["ops"]])
        self.assertEqual(area_g.validate_g(g, set())[0], [])
        nob = by_id["G-nobases-0001"]
        self.assertFalse(nob["request"]["strategy"]["output"]["sequences"])
        self.assertIn("seq-no-bases", [o["id"] for o in nob["graphlet"]["ops"]])

    def test_pairs_run_after_their_graphlet(self):
        g, _ = provisional()
        pos = dict((s["id"], i) for i, s in enumerate(g))
        for s in g:
            for o in s["graphlet"]["ops"]:
                b = o["args"].get("b") if o["tool"] == "graphlet_compare" else None
                if isinstance(b, dict) and "graphlet" in b:
                    self.assertLessEqual(pos[b["graphlet"]], pos[s["id"]], (s["id"], o["id"]))
        by_id = dict((s["id"], s) for s in g)
        radius = by_id["G-pair-radius-0001"]["request"]["strategy"]["bounds"]["max_extension_bp"]
        self.assertEqual(radius * 2, by_id["G-arms2-0001"]["request"]["strategy"]["bounds"]["max_extension_bp"])
        cap = by_id["G-pair-cap-0001"]["request"]["strategy"]["labels"]["max_labels_per_node"]
        self.assertEqual(cap, 2)
        seeds = by_id["G-pair-seeds-0001"]["request"]["seeds"][0]["sequence"]
        self.assertEqual(seeds[:-10], by_id["G-arms2-0001"]["request"]["seeds"][0]["sequence"][10:])
        self.assertIn("self-claims-all", [o["id"] for o in by_id["G-arms2-0001"]["graphlet"]["ops"]])

    def test_one_deliberate_cut_per_kind_and_the_single_service_ops(self):
        g, _ = provisional()
        cuts, probes, cap = {}, [], []
        for s in g:
            ops = dict((o["id"], o) for o in s["graphlet"]["ops"])
            for o in s["graphlet"]["ops"]:
                if o.get("cut_of"):
                    self.assertIn(o["cut_of"], ops)
                    uncut = ops[o["cut_of"]]
                    self.assertEqual(o["tool"], uncut["tool"])
                    cuts.setdefault(o["id"], []).append(s["id"])
                if o.get("probe"):
                    probes.append(o["probe"])
                    self.assertEqual(o["ways"], ["W3"])
                if o.get("by_design") == area_g.SERVICE_CAP_CASE:
                    cap.append(o["args"]["max_bytes"])
        # the summary cut, the ceiling's refusal, the page cut, the row cut, the compare's fields cut (the sequence
        # cut needs the long walk, drawn from the records)
        self.assertEqual(sorted(cuts), ["cmp-cut", "sum-cut", "walks-full-cut", "walks-p1", "walks-too-large"])
        self.assertTrue(all(len(v) == 1 for v in cuts.values()), cuts)
        self.assertEqual(sorted(probes), ["bad_handle", "wrong_token"])
        self.assertEqual(cap, [area_g.MAXB + 1])

    def test_the_owner_reasons_and_the_service_arguments(self):
        g, _ = provisional()
        tools = set()
        for s in g:
            for o in s["graphlet"]["ops"]:
                tools.add(o["tool"])
                ways = o.get("ways") or s["graphlet"]["ways"]["default"]
                if o["tool"] in area_g.LIBRARY_TOOL_ARGS:
                    self.assertNotIn("W3", ways)
                    self.assertIn(o["not_exposed"]["W3"], (area_g.OWNER_LIBRARY_ONLY, area_g.OWNER_NOT_SERVED,
                                                           area_g.OWNER_SMOKE, area_g.REASON_DIRECT_ONLY))
                if "W3" in ways and not o.get("probe"):
                    self.assertEqual(set(o["args"]) - {"cursor"}, set(area_g.SERVICE_ARGS[o["tool"]]),
                                     (s["id"], o["id"]))
                self.assertGreater(o["timeout_s"], 0)
        # every tool of GraphletTools but the backend ones (traverse_capabilities / _resolve / _fetch: a backend)
        self.assertEqual(tools, set(area_g.SERVICE_ARGS) | set(area_g.LIBRARY_TOOL_ARGS))

    def test_counts_within_the_design(self):
        g, _ = provisional()
        n = area_g.summary(g)
        # DESIGN-v8 §4: about 800-1,000 MCP calls and ~5,000 local calls (12 graphlets, 5 pairs)
        self.assertTrue(600 <= n["calls_per_way"]["W3"] <= 1100, n)
        self.assertTrue(sum(n["calls_per_way"][w] for w in ("W1", "W2", "W4")) <= 6000, n)


class CorpusFile(unittest.TestCase):
    def corpus(self, **cands):
        slots = dict((s, {"candidate": None, "why": "left the corpus"}) for s in area_g.SLOT_NAMES)
        slots["two_arms"] = {"candidate": "two_arms.2", "why": "passed"}
        slots["linear"] = {"candidate": "linear.3", "why": "passed"}
        return {"version": 1, "slots": slots, "candidates": cands, "rich": ["two_arms"]}

    def test_the_rehearsal_chooses(self):
        pre = {"arms": ["left", "right"], "walks_min": {"left": 4, "right": 3}, "segments_min": 15}
        c = self.corpus(**{"two_arms.2": {"facts": {"arms": ["left", "right"], "walks": {"left": 4, "right": 3},
                                                    "labels": 2, "summary_cut_bytes": 1234},
                                          "preconditions": pre, "timeouts": {"sum": 3.5},
                                          "reference": {"body_digest": "ab" * 32}}})
        g = area_g.gen_G(FakeCtx(), v7_scenarios(), c, "f" * 40)
        by_id = dict((s["id"], s) for s in g)
        self.assertEqual(sorted(by_id), ["G-arms2-0001", "G-linear-0001", "G-pair-radius-0001",
                                         "G-pair-seeds-0001"])
        s = by_id["G-arms2-0001"]
        self.assertEqual(s["graphlet"]["candidate"], "two_arms.2")
        self.assertEqual(s["graphlet"]["preconditions"], pre)
        self.assertEqual(s["graphlet"]["library_commit"], "f" * 40)
        self.assertEqual(s["graphlet"]["reference"]["body_digest"], "ab" * 32)
        ops = dict((o["id"], o) for o in s["graphlet"]["ops"])
        self.assertEqual(ops["sum"]["timeout_s"], 3.5)
        self.assertEqual(ops["sum-cut"]["args"]["max_bytes"], 1234)
        self.assertEqual(ops["walks-length-left"]["timeout_s"], area_g.LANE_S["interactive"])
        self.assertNotIn("rehearsal-pending", s["tags"])

    def test_a_direct_only_choice_waits_for_the_owner(self):
        # the owner's decision of 2026-10-19: direct-only slots are v8's by default; --no-direct-only leaves them out
        c = self.corpus()
        c["slots"]["annotate"] = {"candidate": "annotate.1", "why": "passed"}
        ids = [s["id"] for s in area_g.gen_G(FakeCtx(), v7_scenarios(), c, None)]
        self.assertIn("G-annot-0001", ids)
        ids = [s["id"] for s in area_g.gen_G(FakeCtx(), v7_scenarios(), c, None, include_direct_only=False)]
        self.assertNotIn("G-annot-0001", ids)

    def test_an_echo_that_differs_is_refused(self):
        c = self.corpus(**{"two_arms.2": {"strategy": {"bounds": {}}}})
        with self.assertRaises(SystemExit):
            area_g.gen_G(FakeCtx(), v7_scenarios(), c, None)

    def test_an_unknown_candidate_is_refused(self):
        c = self.corpus()
        c["slots"]["merges"] = {"candidate": "merges.9"}
        with self.assertRaises(SystemExit):
            area_g.gen_G(FakeCtx(), v7_scenarios(), c, None)


class OpProblems(unittest.TestCase):
    def test_validation_catches_plan_mistakes(self):
        g = {"ways": {"default": area_g.ALL_WAYS}}
        ok = area_g._op("x", "graphlet_walks", area_g.W("left"))
        ok["timeout_s"] = 20.0
        self.assertEqual(area_g.op_problems(ok, g, set(), set()), [])
        missing = dict(ok, args=dict((k, v) for k, v in ok["args"].items() if k != "tail_bp"))
        self.assertTrue(area_g.op_problems(missing, g, set(), set()))
        libarg = dict(ok, args=dict(ok["args"], coordinates=True))
        self.assertTrue(any("library-only argument" in p or "service's arguments" in p
                            for p in area_g.op_problems(libarg, g, set(), set())))
        noreason = dict(ok, ways=["W1", "W2"])
        self.assertTrue(any("reason" in p for p in area_g.op_problems(noreason, g, set(), set())))
        sel = dict(ok, args=dict(ok["args"], min_bp=area_g.ref("later", "/total")))
        self.assertTrue(any("does not run before" in p for p in area_g.op_problems(sel, g, set(), set())))
        self.assertEqual(area_g.op_problems(sel, g, set(), {"later"}), [])
        tool = dict(ok, tool="graphlet_splits", args={"arm": "left", "n": 1, "min_labels_before": 0, "max_bytes": 64})
        self.assertTrue(any("library-only tool" in p for p in area_g.op_problems(tool, g, set(), set())))


@needs_library
class LibrarySignatures(unittest.TestCase):
    def test_every_op_argument_is_the_library_tools(self):
        if LIB not in sys.path:
            sys.path.insert(0, LIB)
        from metagraph.traverse.mcp_tools import GraphletTools
        g, _ = provisional(include_direct_only=True)
        for s in g:
            for o in s["graphlet"]["ops"]:
                if o.get("probe") or "W1" not in (o.get("ways") or s["graphlet"]["ways"]["default"]):
                    continue
                params = set(inspect.signature(getattr(GraphletTools, o["tool"])).parameters) - {"self"}
                self.assertLessEqual(set(o["args"]), params | {"handle", "a"}, (s["id"], o["id"]))
        # the service's tools take a subset of the library's arguments (b for compare: the runner's handle)
        for tool, args in area_g.SERVICE_ARGS.items():
            params = set(inspect.signature(getattr(GraphletTools, tool)).parameters)
            self.assertLessEqual(set(args), params, tool)
        for tool, args in area_g.LIBRARY_TOOL_ARGS.items():
            params = set(inspect.signature(getattr(GraphletTools, tool)).parameters) - {"self", "budget", "cursor"}
            self.assertEqual(set(args) | ({"handle"} if "handle" in params else set()), params, tool)


def read_json(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as fh:
        return json.load(fh)


REAL = os.path.join(DATA, "real", "mini_refseq")
FIXTURES = {
    "G-arms2": os.path.join(REAL, "mini_ndm1_rc__switch1.graphlet.json.gz"),       # two arms, merges, cut lists
    "G-trace": os.path.join(REAL, "mini_ndm1__trace.graphlet.json.gz"),            # trace; compare a
    "G-annot": os.path.join(REAL, "mini_ndm1__annotate_exh.graphlet.json.gz"),     # annotate; compare b
    "G-coords": os.path.join(REAL, "mini_rep_00__trace_coords.graphlet.json.gz"),  # record coordinates
    "G-arm1": os.path.join(DATA, "documents", "merge.graphlet.json"),              # one arm, two merges
    "G-caps": os.path.join(DATA, "documents", "caps.graphlet.json"),               # one arm, no bases
}


@needs_library
class PlansOnRealGraphlets(unittest.TestCase):
    """The op plans run by the local runner (W1) on the library's fixture graphlets: every op answers what
    area_g.expected_outcome says (no selector failure, no unexpected error), every cut op shows a cut of its uncut
    op (compare_graphlets.compare_cut, lane C's reading)."""

    @classmethod
    def setUpClass(cls):
        import rehearse_g
        import run_graphlet_local as rgl
        cls.tmp = tempfile.mkdtemp(prefix="area-g-")
        cls.facts = {}
        for sid, path in FIXTURES.items():
            full = path.replace(".graphlet.", ".full.")
            cls.facts[sid] = rehearse_g.document_facts(read_json(path),
                                                       read_json(full) if os.path.exists(full) else None)
        cls.facts["G-coords"]["coordinates"] = True
        spec = {"G-arms2": (True, [("G-arms2", False)],
                            {"summary_cut": True, "ceiling_cut": True, "page_cut": True, "probes": True,
                             "service_cap": True, "view_and_files": True, "row_cut": True}),
                "G-trace": (True, [("G-annot", True)], {"compare_cut": "G-annot", "compare_paged": "G-annot"}),
                "G-annot": (False, [], {}), "G-coords": (False, [], {}), "G-arm1": (True, [], {}),
                "G-caps": (False, [], {})}
        cls.plan = {"scenarios": []}
        for sid, (rich, pairs, extras) in spec.items():
            ops = area_g.op_plan(sid, "test", cls.facts[sid], rich, pairs, extras=extras)
            area_g.set_timeouts(ops)
            cls.plan["scenarios"].append({"id": sid, "graphlet": {"ops": ops,
                                                                  "ways": {"default": area_g.ALL_WAYS}}})
        with open(os.path.join(cls.tmp, "plan.json"), "w") as fh:
            json.dump(cls.plan, fh)
        argv = ["--plan", os.path.join(cls.tmp, "plan.json"), "--library-root", LIB, "--out", cls.tmp, "--ways",
                "W1", "--run-id", "area-g-test", "--quiet"]
        for sid, path in sorted(FIXTURES.items()):
            argv += ["--doc", "%s=%s" % (sid, path)]
        assert rgl.main(argv) == 0
        rows, bad = sc.read_jsonl(os.path.join(cls.tmp, "arms", "dl-lib", "graphlet_ops.jsonl"))
        assert bad == 0
        cls.recs = dict(((r["scenario"], r["op"]), r) for r in rows)

    @classmethod
    def tearDownClass(cls):
        import shutil
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_every_op_answers_as_expected(self):
        wrong = []
        n = 0
        for s in self.plan["scenarios"]:
            for o in s["graphlet"]["ops"]:
                if "W1" not in (o.get("ways") or area_g.ALL_WAYS):
                    continue
                n += 1
                r = self.recs.get((s["id"], o["id"]))
                got = "missing" if r is None else "selector" if r["selector_error"] else \
                    (r["error"]["code"] if r["error"] else "ok")
                want = area_g.expected_outcome(o, self.facts[s["id"]])
                if got != want:
                    wrong.append((s["id"], o["id"], got, want))
        self.assertEqual(wrong, [])
        self.assertGreater(n, 700)

    def test_every_cut_op_shows_a_cut(self):
        import compare_graphlets as cg
        seen = 0
        for s in self.plan["scenarios"]:
            for o in s["graphlet"]["ops"]:
                if not o.get("cut_of"):
                    continue
                r, ref = self.recs[(s["id"], o["id"])], self.recs[(s["id"], o["cut_of"])]
                c = cg.compare_cut("W1", r, ref, cg.OPTIONAL_KEYS.get(o["tool"], frozenset()))
                self.assertEqual(c["outcome"], "agree", (s["id"], o["id"], c))
                seen += 1
        self.assertEqual(seen, 5)   # G-arms2: sum-cut, walks-too-large, walks-p1, walks-full-cut; G-trace: cmp-cut

    def test_whole_lists_fit_the_service_ceiling(self):
        # DESIGN-v8 §5: every W3 whole-list call fits uncut at 262,144 bytes (W3's longer handle aside: a margin)
        for (sid, oid), r in self.recs.items():
            args = r.get("args_resolved") or {}
            if args.get("max_bytes") == area_g.MAXB and r.get("answer") is not None:
                self.assertFalse(r["answer"].get("row_truncated"), (sid, oid))
                self.assertLess(len(json.dumps(r["answer"])), area_g.MAXB - 4096, (sid, oid))


@unittest.skipUnless(have_cache(), "no record cache")
class BuiltV8(unittest.TestCase):
    def test_v8_is_v7_then_area_g(self):
        import build_manifest as bm
        m = bm.build(bm.DEFAULT_SEED, fetch=False, manifest_version=8)
        v7 = sc.read_json(V7)
        n = len(v7["scenarios"])
        self.assertEqual(m["manifest_version"], 8)
        self.assertEqual(ml.canonical_json(m["scenarios"][:n]), ml.canonical_json(v7["scenarios"]))
        self.assertEqual(m["relations_sha256"], v7["relations_sha256"])
        self.assertTrue(all(s["area"] == "G" for s in m["scenarios"][n:]))
        self.assertTrue(any(s["id"] == "G-long-0001" for s in m["scenarios"][n:]))
        self.assertEqual(bm.smoke_of(m)["content_sha256"], bm.smoke_of(v7)["content_sha256"])
        self.assertTrue(m["graphlets"]["corpus"]["provisional"])


if __name__ == "__main__":
    unittest.main()
