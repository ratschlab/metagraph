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


def areas_mod():
    import areas
    return areas


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
                               "G-linear-0001", "G-fail-0001", "G-trace-0001", "G-coords-0001", "G-pair-annot-0001",
                               "G-annot-0001", "G-nobases-0001"])
        skipped = dict((x["scenario"], x["reason"]) for x in ctx.skipped)
        for sid in ("G-big-0001", "G-twodb-0001", "G-long-0001"):
            self.assertIn(sid, skipped)
        self.assertTrue(skipped["G-big-0001"].startswith("slot over_8mib: v9"))
        problems, pending = area_g.validate_g(g, set())
        self.assertEqual(problems, [])
        self.assertTrue(any("preconditions" in p for p in pending))
        self.assertTrue(any("library_commit" in p for p in pending))
        self.assertTrue(any(area_g.SERVICE_CAP_CASE in p for p in pending))
        direct = ("G-coords-0001", "G-annot-0001", "G-nobases-0001", "G-pair-annot-0001")
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


class ByDesignTable(unittest.TestCase):
    """async's three area G cases (BY-DESIGN-G-CASES.json, 2026-10-20) are in the bundled table, load and validate,
    and every case an op plan names is one of them."""

    def test_the_g_cases_load_and_the_plans_name_them(self):
        import compare as cmp
        cases, info = cmp.load_by_design(cmp.DEFAULT_BY_DESIGN_G)   # area G's own table (VERIFY-D VD1)
        for cid in (area_g.SERVICE_CAP_CASE, area_g.BIG_BODY_CASE, area_g.BIG_BODY_CASE + ".summary",
                    area_g.MGT_EXPORT_CASE):
            self.assertIn(cid, cases)
            self.assertEqual(cmp.graphlet_case_problems(cases[cid]), [], cid)
        g, _ = provisional()
        named = set(o["by_design"] for s in g for o in s["graphlet"]["ops"] if o.get("by_design"))
        self.assertTrue(named <= set(cases), named - set(cases))
        self.assertEqual(area_g.validate_g(g, set(cases))[1],
                         sorted(p for p in area_g.validate_g(g, set(cases))[1] if "by-design case" not in p))


class OpIds(unittest.TestCase):
    """The builder refuses an op named after an arm (2026-10-20, async's runner review): its raw pages
    raw/<scenario>.mcp.<page>.json.gz would overwrite the MCP execution's raw record raw/<scenario>.mcp.0.json.gz
    in arms/mcp (async's runner exits 4 on such a plan); op ids keep to a raw name's characters."""

    def test_an_op_named_mcp_is_refused(self):
        g, _ = provisional()
        self.assertEqual(area_g.validate_g(g, None)[0], [])
        for bad in ("mcp", "direct-local", "direct-remote", "a/b", "", "-x", None):
            g2 = copy.deepcopy(g)
            g2[0]["graphlet"]["ops"][0]["id"] = bad
            problems = area_g.validate_g(g2, None)[0]
            self.assertTrue(any(("op id %r" % (bad,)) in p for p in problems), (bad, problems))
        self.assertIn("collide with the execution's raw record", area_g.op_id_problems("mcp")[0])
        for ok in ("mcp-walks", "sum", "walks-id-left", "seq-0-left", "exp.mgt", "x_1"):
            self.assertEqual(area_g.op_id_problems(ok), [], ok)

    def test_every_planned_op_id_is_valid(self):
        g, _ = provisional(include_direct_only=True)
        bad = [(s["id"], o["id"]) for s in g for o in s["graphlet"]["ops"] if area_g.op_id_problems(o["id"])]
        self.assertEqual(bad, [])


class ByDesignLinks(unittest.TestCase):
    """Every op links the by-design case that applies to it (async's three, 2026-10-20) and no other; the probes stay
    within the service's 262,144 (its cap refusal comes before the token check); the summary of a body over the
    in-service ceiling sends max_bytes 262,144 (the .summary case; v9 slot, the rule set in the generator)."""

    BIG = {"over_8mib": True, "arms": ["left", "right"]}

    def test_a_body_over_the_ceiling(self):
        ops = area_g.op_plan("G-big-0001", "over_8mib", dict(self.BIG), False)
        by = dict((o["id"], o.get("by_design")) for o in ops)
        summary = next(o for o in ops if o["tool"] == "graphlet_summary")
        self.assertEqual(summary["args"]["max_bytes"], area_g.MAXB)
        self.assertEqual(summary["by_design"], area_g.BIG_BODY_CASE + ".summary")
        for o in ops:
            if o["tool"] == "graphlet_export" and o["args"].get("format") == "mgt":
                self.assertEqual(by[o["id"]], area_g.MGT_EXPORT_CASE, o["id"])   # the stored body's link
            elif o["tool"] != "graphlet_summary":
                self.assertEqual(by[o["id"]], area_g.BIG_BODY_CASE, o["id"])
        g = {"ways": {"default": area_g.ALL_WAYS}, "_facts": {"over_8mib": True}}
        for o in ops:
            o["timeout_s"] = 20.0
            self.assertFalse([m for m in area_g.op_problems(o, g, set(), set(by)) if "by_design" in m], o["id"])

    def test_the_rule(self):
        case = area_g.by_design_case
        walks = area_g._op("w", "graphlet_walks", area_g.W("left", max_bytes=area_g.MAXB + 1))
        self.assertEqual(case(walks, {}), area_g.SERVICE_CAP_CASE)
        self.assertEqual(case(walks, {"over_8mib": True}), area_g.SERVICE_CAP_CASE)   # the cap's refusal first
        self.assertIsNone(case(area_g._op("w", "graphlet_walks", area_g.W("left", max_bytes=area_g.MAXB)), {}))
        self.assertIsNone(case(area_g._op("w", "graphlet_walks", area_g.W("left", max_bytes=16 * 2 ** 20 + 1)), {}))
        seq = area_g._op("s", "graphlet_sequence", area_g.S("left", walk=0, max_bytes=area_g.MAXB + 1))
        self.assertIsNone(case(seq, {}))              # the cap case does not name graphlet_sequence
        for fmt, want in (("fasta", area_g.BIG_BODY_CASE), ("gfa", area_g.BIG_BODY_CASE),
                          ("json", area_g.BIG_BODY_CASE), ("mgt", area_g.MGT_EXPORT_CASE)):
            self.assertEqual(case(area_g._op("e", "graphlet_export", area_g.E(format=fmt)), {"over_8mib": True}),
                             want, fmt)
        cmp_ = area_g._op("c", "graphlet_compare", area_g.C("G-x-0001", mode="labels", max_bytes=2048))
        self.assertEqual(case(cmp_, {"over_8mib": True}), area_g.BIG_BODY_CASE)
        self.assertIsNone(case(cmp_, {}))
        direct = dict(walks, ways=["W1", "W2"])
        self.assertIsNone(case(direct, {}))           # no W3 call: nothing the service answers by design

    def test_validation_names_a_missing_or_a_wrong_link(self):
        g = {"ways": {"default": area_g.ALL_WAYS}, "_facts": {}}
        cap = area_g._op("w", "graphlet_walks", area_g.W("left", max_bytes=area_g.MAXB + 1), timeout_s=20.0)
        self.assertTrue(any("by_design None" in m for m in area_g.op_problems(cap, g, set(), set())))
        cap["by_design"] = area_g.SERVICE_CAP_CASE
        self.assertEqual(area_g.op_problems(cap, g, set(), set()), [])
        plain = area_g._op("w", "graphlet_walks", area_g.W("left"), timeout_s=20.0, by_design=area_g.BIG_BODY_CASE)
        self.assertTrue(any("by_design" in m for m in area_g.op_problems(plain, g, set(), set())))

    def test_probes_stay_within_the_cap(self):
        ops = area_g.op_plan("G-arms2-0001", "two_arms", {"arms": ["left", "right"]}, True, extras={"probes": True})
        probes = [o for o in ops if o.get("probe")]
        self.assertEqual(sorted(o["probe"] for o in probes), ["bad_handle", "wrong_token"])
        for o in probes:
            self.assertLessEqual(o["args"]["max_bytes"], area_g.MAXB, o["id"])
            self.assertNotIn("by_design", o)
        g = {"ways": {"default": area_g.ALL_WAYS}, "_facts": {}}
        over = dict(probes[0], args=dict(probes[0]["args"], max_bytes=area_g.MAXB + 1), timeout_s=20.0)
        self.assertTrue(any("cap's refusal" in m for m in area_g.op_problems(over, g, set(), set())))

    def test_a_failed_seed_links_no_case_and_leaves_w4_out(self):
        # VERIFY-D VD2 / VD3: no way has a graphlet; W4 not exposed (no document), and no case whose sides could hold
        ops = area_g.op_plan("G-fail-0001", "failed_seed", {"failed": True, "arms": ["left"]}, False,
                             extras={"other": "G-arms2-0001"})
        self.assertTrue(ops)
        for o in ops:
            self.assertNotIn("by_design", o, o["id"])
            self.assertEqual((o["ways"], o["not_exposed"]), (area_g.FAILED_SEED_WAYS,
                                                              {"W4": area_g.REASON_FAILED_NO_DOC}), o["id"])
        mgt = next(o for o in ops if o["id"] == "exp-mgt")
        self.assertIsNone(area_g.by_design_case(mgt, {"failed": True}))
        self.assertEqual(area_g.by_design_case(dict(mgt, ways=None), {}), area_g.MGT_EXPORT_CASE)
        g, _ = provisional()
        fail = next(s for s in g if s["id"] == "G-fail-0001")
        self.assertFalse([o["id"] for o in fail["graphlet"]["ops"] if o.get("by_design")])
        self.assertEqual(area_g.validate_g(g, None)[0], [])
        # validate_g refuses the mgt case on a failed seed
        bad = copy.deepcopy(g)
        b = next(s for s in bad if s["id"] == "G-fail-0001")
        next(o for o in b["graphlet"]["ops"] if o["id"] == "exp-mgt")["by_design"] = area_g.MGT_EXPORT_CASE
        self.assertTrue([p_ for p_ in area_g.validate_g(bad, None)[0] if "G-fail-0001/exp-mgt" in p_])

    def test_every_mgt_export_links_the_stored_body_case(self):
        g, _ = provisional(include_direct_only=True)
        seen = 0
        for s in g:
            default = s["graphlet"]["ways"]["default"]
            for o in s["graphlet"]["ops"]:
                if o["tool"] != "graphlet_export" or o["args"].get("format") != "mgt":
                    continue
                w3 = "W3" in (o.get("ways") or default) and not s["graphlet"]["preconditions"].get("no_graphlet")
                self.assertEqual(o.get("by_design"), area_g.MGT_EXPORT_CASE if w3 else None, (s["id"], o["id"]))
                seen += w3
        self.assertGreater(seen, 0)
        self.assertIsNone(area_g.by_design_case(area_g._op("e", "graphlet_export", area_g.E(format="mgt", arm="left")),
                                                {}))
        # a direct-only scenario keeps its mgt export (W1, W2, W5), without the case
        direct = [s for s in g if "direct-only" in s["tags"] and "pair" not in s["tags"]]
        self.assertTrue(direct)
        for s in direct:
            ids = [o["id"] for o in s["graphlet"]["ops"]]
            self.assertIn("exp-mgt", ids, s["id"])

    def test_the_built_plans_link_every_case(self):
        g, _ = provisional()
        self.assertFalse([p for p in area_g.validate_g(g, None)[0] if "by_design" in p or "probe" in p])
        named = [(s["id"], o["id"], o["by_design"]) for s in g for o in s["graphlet"]["ops"] if o.get("by_design")]
        self.assertIn(area_g.SERVICE_CAP_CASE, set(n[2] for n in named))
        for s in g:
            for o in s["graphlet"]["ops"]:
                if o.get("probe"):
                    self.assertLessEqual(o["args"].get("max_bytes") or 0, area_g.MAXB)


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


class McpRepeats(unittest.TestCase):
    """Area G and the stop scenarios run once on MCP (suite_common.mcp_repeats, 2026-10-20: the service's result
    cache answers a repeat from the first; a stop repeat would stop another job); v1-v7 unchanged (no scenario of
    theirs has a graphlet block or a stop flow)."""

    def test_the_rule(self):
        g = {"id": "G-x-0001", "graphlet": {}, "mcp": {"tool": "traverse_fetch"}}
        for k, we in ((0, 10), (10, 10), (3, 0)):
            self.assertEqual(sc.mcp_repeats(k, g, we), [0])
        self.assertEqual(sc.mcp_repeats(0, dict(g, repeat=3), 10), [0])
        self.assertEqual(sc.mcp_repeats(0, dict(g, timing_class="long"), 10), [])
        stop = {"id": "S-stop-0001", "mcp": {"tool": "traverse_fetch", "flow": {"kind": "stop"}}, "repeat": 3}
        self.assertEqual(sc.mcp_repeats(0, stop, 10), [0])
        other = {"id": "X", "mcp": {"tool": "traverse_fetch", "flow": {"kind": "resume"}}}
        self.assertEqual(sc.mcp_repeats(0, other, 10), [0, 1])
        self.assertEqual(sc.mcp_repeats(0, dict(other, repeat=3), 10), [0, 1, 2])
        self.assertEqual(sc.mcp_repeats(1, {"id": "Y", "mcp": {"tool": "x"}}, 10), [0, 1])

    def test_v1_to_v7_unchanged(self):
        def before(index, scenario, warm_every):   # the rule before 2026-10-20
            if scenario.get("timing_class") == "long":
                return []
            if (scenario.get("mcp") or {}).get("not_implemented"):
                return [0]
            return list(range(max(2, sc.repeat_count(index, scenario, warm_every))))
        for v in range(1, 8):
            m = sc.read_json(os.path.join(SUITE, "manifest", "suite_v%d.json" % v))
            for k, s in enumerate(m["scenarios"]):
                self.assertFalse(sc.single_mcp_execution(s), (v, s["id"]))
                for we in (10, 0):
                    self.assertEqual(sc.mcp_repeats(k, s, we), before(k, s, we), (v, s["id"], we))

    def test_the_paired_plan_of_a_graphlet_scenario(self):
        g, _ = provisional()
        m = {"scenarios": g}
        tmp = tempfile.mkdtemp(prefix="g-plan-")
        try:
            for item in sc.paired_plan(m, tmp, warm_every=10, arms=("direct-local", "mcp")):
                if item["scenario"].startswith("G-warm-"):
                    self.assertIsNone(item["mcp_cold"])
                    continue
                self.assertEqual((item["mcp_cold"], item["mcp_warm"]), (True, []), item["scenario"])
        finally:
            import shutil
            shutil.rmtree(tmp, ignore_errors=True)


class StopScenarios(unittest.TestCase):
    """The owner's decision of 2026-10-19: the stop scenarios enter only with --stop-scenarios (async's flow). VERIFY-D
    VD5 / VD6: MCP only (the direct arms skip them), the flow's expect list only, a walk CPU-bound enough to be
    running past delay_s + margin on a warm host too."""

    def test_two_flow_scenarios(self):
        import build_manifest as bm
        out = areas_mod().gen_stop(FakeCtx(), v7_scenarios())
        self.assertEqual([s["id"] for s in out], list(areas_mod().STOP_IDS))
        walk, res = out
        flow_keys = {"kind", "after", "delay_s", "observe", "poll_s", "release_wait_s"}   # async's FLOW_KEYS
        stop_checks = ["stop_answer", "job_terminal", "attempt_state", "token_release"]
        for s in out:
            self.assertLessEqual(set(s["mcp"]["flow"]), flow_keys, s["id"])
            self.assertEqual([e["check"] for e in s["expect"]], ["status"] + stop_checks, s["id"])   # VD5
            self.assertEqual(s["direct_not_run"], areas_mod().STOP_DIRECT_NOT_RUN)                    # VD6
        self.assertEqual(walk["route"], "POST /traverse")
        self.assertEqual(walk["mcp"]["tool"], "traverse_fetch")
        sw = areas_mod().STOP_WALK
        self.assertEqual(walk["request"]["seeds"], [{"labels": ["NC_007795.1"], "sequence": sw["seed"]}])
        self.assertEqual((walk["request"]["strategy"]["direction"], walk["truth"]["start"]), ("right", 1512368))
        self.assertEqual(walk["request"]["strategy"]["bounds"]["max_extension_bp"], 10000)
        for s in out:
            self.assertEqual(s["mcp"]["flow"]["delay_s"], 0.0)   # R2-1: the stop right after running
        self.assertEqual(walk["mcp"]["arguments"]["overrides"]["bounds"]["max_extension_bp"], 10000)
        sa = dict((e["check"], e["args"]) for e in walk["expect"])["stop_answer"]
        self.assertTrue(sa["cancel_supported"])
        self.assertEqual(res["route"], "POST /resolve")
        self.assertEqual(res["mcp"]["tool"], "traverse_resolve")
        sa = dict((e["check"], e["args"]) for e in res["expect"])["stop_answer"]
        self.assertEqual((sa["cancel_supported"], sa["cancel"]), (True, "disconnect"))
        self.assertEqual(res["request"]["sequence"], areas_mod().STOP_RESOLVE["query"])
        self.assertEqual(res["request"]["discover"], {"kind": "header", "max_labels": 1000})
        self.assertEqual(res["mcp"]["arguments"]["sequence"], areas_mod().STOP_RESOLVE["query"])
        self.assertNotIn("21", res["title"] + json.dumps(res["truth"]))   # R2-2: no unsourced timing claim
        self.assertEqual(bm.validate(out), [])

    @unittest.skipUnless(have_cache(), "no record cache")
    def test_the_walk_meets_no_repeat_and_the_sequences_are_the_records(self):
        # R2-1: every 30-mer from the seed to 10,260 bp past it occurs once in sau (either strand), so the header
        # walk branches nowhere within its 10 kbp; the constants are the records' bases
        recs = dict((r.key, r) for r in truth.load_with_sequences(fetch=False))
        sw, sr = areas_mod().STOP_WALK, areas_mod().STOP_RESOLVE
        sau = recs["sau"].seq
        self.assertEqual(sau[sw["start"]:sw["start"] + len(sw["seed"])], sw["seed"])
        self.assertEqual(recs["eco"].seq[sr["start"]:sr["start"] + 500], sr["query"])
        k = 30
        comp = str.maketrans("ACGT", "TGCA")

        def canon(w):
            r = w.translate(comp)[::-1]
            return min(w, r)
        window = sau[sw["start"]:sw["start"] + sw["no_repeat_bp"] + k]
        want = set(canon(window[i:i + k]) for i in range(len(window) - k + 1))
        self.assertEqual(len(want), len(window) - k + 1)   # no repeat inside the window
        hits = 0
        for i in range(len(sau) - k + 1):
            if canon(sau[i:i + k]) in want:
                hits += 1
        self.assertEqual(hits, len(want))                 # and none elsewhere in the genome

    def test_only_the_mcp_arm_runs_them(self):
        out = areas_mod().gen_stop(FakeCtx(), v7_scenarios())
        selected = list(enumerate(out))
        for arm in sc.DIRECT_ARMS:
            planned, skipped = sc.plan_executions(selected, arm, 10)
            self.assertEqual(planned, [], arm)
            self.assertEqual([x["scenario"] for x in skipped], list(areas_mod().STOP_IDS))
            self.assertTrue(all(x["reason"].startswith("direct_not_run: ") for x in skipped))
        tmp = tempfile.mkdtemp(prefix="stop-plan-")
        try:
            for item in sc.paired_plan({"scenarios": out}, tmp, warm_every=10,
                                       arms=("direct-local", "direct-remote", "mcp")):
                self.assertEqual((item["direct-local"], item["direct-remote"], item["mcp_cold"], item["mcp_warm"]),
                                 (None, None, True, []), item["scenario"])
        finally:
            import shutil
            shutil.rmtree(tmp, ignore_errors=True)


class Coverage(unittest.TestCase):
    """VERIFY-D VD8: compares that answer true (an uncut self-compare: merges.1) and unknown (depth 0, no bases),
    compare rows (D1's constrain vs annotate pair, direct-only), a limitations listing and a compare that page."""

    def test_the_gaps_are_closed(self):
        g, _ = provisional(include_direct_only=True)
        by = dict((s["id"], s) for s in g)
        ops = dict((sid, dict((o["id"], o) for o in s["graphlet"]["ops"])) for sid, s in by.items())
        for sid in ("G-merge-0001", "G-depth0-0001", "G-nobases-0001"):
            for mode in ("claims", "walks", "labels", "prefix_subset"):
                o = ops[sid]["self-%s-all" % mode]
                self.assertEqual(o["args"]["b"], {"graphlet": sid}, (sid, mode))
        # the pair: the compares on the annotate graphlet (W1, W2, W5) against its constrain partner
        partner = by["G-pair-annot-0001"]
        self.assertEqual(partner["graphlet"]["ways"]["default"], area_g.DIRECT_ONLY_WAYS)
        self.assertEqual(partner["graphlet"]["not_exposed"]["*"], area_g.REASON_ANNOTATE_PAIR)
        self.assertEqual(partner["mcp"], {"not_implemented": area_g.REASON_ANNOTATE_PAIR})
        self.assertNotEqual(partner["request"]["strategy"]["labels"].get("mode"), "annotate")
        self.assertFalse([oid for oid in ops["G-pair-annot-0001"] if oid.startswith("cmp")])
        annot = ops["G-annot-0001"]
        self.assertEqual(annot["cmp-claims-all"]["args"]["b"], {"graphlet": "G-pair-annot-0001"})
        self.assertEqual(annot["cmp-claims-all"]["not_exposed"]["W3"], area_g.REASON_DIRECT_ONLY)
        paged = annot["cmp-paged"]
        self.assertEqual((paged["args"]["mode"], paged["args"]["max_bytes"], paged["paging"]),
                         ("claims", area_g.COMPARE_PAGE, "all"))
        self.assertNotIn("cmp-paged", ops["G-pair-radius-0001"])
        self.assertIn("cmp-cut", ops["G-pair-radius-0001"])
        for sid, s in by.items():
            if "lims-paged" in ops[sid]:
                self.assertEqual(ops[sid]["lims-paged"]["args"]["max_bytes"], area_g.LIMS_PAGE, sid)
        self.assertEqual(area_g.validate_g(g, None)[0], [])
        # without the direct-only slots the pair goes and the radius partner pages its compare as before
        g2, _ = provisional(include_direct_only=False)
        ids2 = [s["id"] for s in g2]
        self.assertNotIn("G-pair-annot-0001", ids2)
        radius = next(s for s in g2 if s["id"] == "G-pair-radius-0001")
        self.assertIn("cmp-paged", [o["id"] for o in radius["graphlet"]["ops"]])


class Defaults(unittest.TestCase):
    """VERIFY-D VD13: areas.gen_G defaults to the direct-only slots as area_g.gen_G, choose and build do; the tables'
    service_commit states the commit their cases cite (the v7 table 22d034ba, area G's own f994e3e2)."""

    def test_the_defaults_agree(self):
        import inspect
        import build_manifest as bm
        import compare as cmp
        for fn in (areas_mod().gen_G, area_g.gen_G, area_g.choose, bm.build):
            self.assertIs(inspect.signature(fn).parameters["include_direct_only"].default, True, fn.__name__)
        self.assertTrue(cmp.load_by_design(cmp.DEFAULT_BY_DESIGN)[1]["service_commit"].startswith("22d034ba"))
        gcases, ginfo = cmp.load_by_design(cmp.DEFAULT_BY_DESIGN_G)
        self.assertTrue(ginfo["service_commit"].startswith("f994e3e2"))
        self.assertEqual(set(c.get("source_commit") for c in gcases.values()) - {None}, {ginfo["service_commit"]})


class Freeze(unittest.TestCase):
    """VERIFY-D VD12: --freeze refuses a provisional corpus and anything still pending."""

    def test_freeze_problems(self):
        import build_manifest as bm
        m = {"manifest_version": 8, "graphlets": {"corpus": {"file": None, "provisional": True},
                                                  "pending": ["G-x-0001: preconditions (the rehearsal freezes them)"]}}
        probs = bm.freeze_problems(m)
        self.assertTrue(any("corpus is provisional" in p for p in probs), probs)
        self.assertIn("pending: G-x-0001: preconditions (the rehearsal freezes them)", probs)
        m["graphlets"] = {"corpus": {"file": "truth/area_g/corpus.json", "provisional": False}, "pending": []}
        self.assertEqual(bm.freeze_problems(m), [])
        self.assertEqual(bm.freeze_problems({"manifest_version": 7}), [])
        self.assertTrue(os.path.exists(bm.DEFAULT_CORPUS))
        corpus = area_g.load_corpus(bm.DEFAULT_CORPUS)
        self.assertFalse(corpus.get("provisional"))


@unittest.skipUnless(have_cache(), "no record cache")
class BuiltV8(unittest.TestCase):
    def test_v8_is_v7_with_t_walk_0024s_warm_ups_then_area_g(self):
        # the owner's decisions of 2026-10-19: v8 differs from v7 only by T-walk-0024's warm-ups and area G
        import build_manifest as bm
        m = bm.build(bm.DEFAULT_SEED, fetch=False, manifest_version=8)
        v7 = sc.read_json(V7)
        n = len(v7["scenarios"])
        self.assertEqual(m["manifest_version"], 8)
        ids7 = [s["id"] for s in v7["scenarios"]]
        i = ids7.index("T-walk-0024")
        warm = list(areas_mod().T_WALK_0024_REGION["ids"])
        self.assertEqual(len(warm) % 10, 0)   # a multiple of ten keeps every later index's warm repeat
        head = m["scenarios"][:i] + m["scenarios"][i + len(warm):n + len(warm)]
        self.assertEqual(ml.canonical_json(head), ml.canonical_json(v7["scenarios"]))
        self.assertEqual([s["id"] for s in m["scenarios"][i:i + len(warm)]], warm)
        for s in m["scenarios"][i:i + len(warm)]:
            self.assertEqual((s["route"], s["warms"], s["request"]["labels"]), ("POST /resolve", "T-walk-0024",
                                                                                ["243273"]))
            self.assertEqual(s["mcp"]["tool"], "traverse_resolve")
        self.assertEqual(m["relations_sha256"], v7["relations_sha256"])
        rest = m["scenarios"][n + len(warm):]
        self.assertTrue(all(s["area"] == "G" for s in rest))
        ids = [s["id"] for s in rest]
        k = ids.index("G-long-0001")
        self.assertEqual(ids[k - 2:k], ["G-warm-0001", "G-warm-0002"])
        for s in rest[k - 2:k]:
            self.assertEqual((s["timing_class"], s["warms"], s["route"]), ("long", "G-long-0001", "POST /resolve"))
            self.assertIn("not_implemented", s["mcp"])
            self.assertLessEqual(len(s["request"]["sequence"]), 10000)
        self.assertEqual(bm.smoke_of(m)["content_sha256"], bm.smoke_of(v7)["content_sha256"])
        # the rehearsal's corpus is committed (VERIFY-D VD12): the default build reads it, the pin from it
        self.assertEqual(m["graphlets"]["corpus"]["file"], "truth/area_g/corpus.json")
        self.assertFalse(m["graphlets"]["corpus"]["provisional"])
        self.assertEqual(m["graphlets"]["library_commit"], "44583b51307ed3242e19c051aaf133c671107956")
        # the constrain / annotate pair's partner (VERIFY-D VD8) waits for its rehearsal (REHEARSAL S)
        self.assertEqual(m["graphlets"]["pending"], ["G-pair-annot-0001: preconditions (the rehearsal freezes them)"])
        self.assertEqual(bm.freeze_problems(m), ["pending: G-pair-annot-0001: preconditions (the rehearsal freezes "
                                                 "them)"])
        # the paired plan's warm repeats of every v7 scenario unchanged
        w7 = dict((s["id"], k % 10 == 0) for k, s in enumerate(v7["scenarios"]))
        w8 = dict((s["id"], k % 10 == 0) for k, s in enumerate(m["scenarios"]))
        self.assertEqual([x for x in ids7 if w7[x] != w8[x]], [])
        # the MCP plan: every v7 scenario's repeats as in v7; area G once (the warm-ups, long class, not at all)
        for we in (10, 0):
            r7 = dict((s["id"], sc.mcp_repeats(k, s, we)) for k, s in enumerate(v7["scenarios"]))
            r8 = dict((s["id"], sc.mcp_repeats(k, s, we)) for k, s in enumerate(m["scenarios"]))
            self.assertEqual([x for x in ids7 if r7[x] != r8[x]], [], we)
            for s in rest:
                self.assertEqual(r8[s["id"]], [] if s["timing_class"] == "long" else [0], s["id"])
        stop = bm.build(bm.DEFAULT_SEED, fetch=False, manifest_version=8, stop_scenarios=True)
        for k, s in enumerate(stop["scenarios"]):
            if s["id"] in areas_mod().STOP_IDS:
                self.assertEqual(sc.mcp_repeats(k, s, 10), [0], s["id"])


if __name__ == "__main__":
    unittest.main()
