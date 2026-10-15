"""Offline tests of the v5 builder (lane V5): the triage-v4 suite errors, the verified review findings on the builder
and truth, the gate's resolve / traverse scenarios decisive under load. The tests that read the records skip without
the record cache (SUITE_TRUTH_CACHE or truth/cache)."""
import json
import os
import random
import sys
import types
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.dirname(HERE)
sys.path.insert(0, SUITE)

import areas  # noqa: E402
import checks  # noqa: E402
import manifest_lib as ml  # noqa: E402
import seqtools as st  # noqa: E402
import truth  # noqa: E402

_BUILDS = {}


def have_cache():
    return os.path.isdir(truth.CACHE_DIR) and bool(os.listdir(truth.CACHE_DIR))


def built(areas_=None):
    """The manifest the builder makes from the record cache (only |areas_| when given), built once per test run."""
    import build_manifest as bm
    key = tuple(sorted(areas_)) if areas_ else ()
    if key not in _BUILDS:
        _BUILDS[key] = bm.build(bm.DEFAULT_SEED, fetch=False, only_areas=set(areas_) if areas_ else None)
    return _BUILDS[key]


def by_id(m):
    return dict((s["id"], s) for s in m["scenarios"])


def synthetic_record():
    """A record whose window w (at 40) is unique on both strands as a whole string, while two k-mers of
    revcomp(w) occur in the record (an inverted copy of w's first 31 bases at 180): the outside review's
    counterexample."""
    rng = random.Random(42)
    w = "".join(rng.choice("ACGT") for _ in range(100))
    seq = ("".join(rng.choice("ACGT") for _ in range(40)) + w + "".join(rng.choice("ACGT") for _ in range(40))
           + st.revcomp(w[:31]) + "".join(rng.choice("ACGT") for _ in range(40)))
    rec = types.SimpleNamespace(key="demo", seq=seq, offset=0, accession="SYNTHETIC.1", column="1", seq_id=0,
                                length=len(seq), region=None)
    return rec, w


class ReverseNegatives(unittest.TestCase):
    """VERIFY-D finding 8 / P4-builder-truth-02: a reverse-complement negative is checked k-mer by k-mer against the
    whole record, never only as a whole string."""

    def test_counterexample_refused(self):
        rec, w = synthetic_record()
        ctx = ml.Ctx([], 1)
        self.assertEqual(ml.reverse_kmers_in_record(rec, w), [68, 69])
        with self.assertRaises(ValueError):
            ctx.rc_negative_window(random.Random(0), rec, 100, where=(40, 40))

    def test_clean_record_accepted(self):
        rec, w = synthetic_record()
        rec.seq = rec.seq[:180] + "".join(random.Random(1).choice("ACGT") for _ in range(31)) + rec.seq[211:]
        t = ml.Ctx([], 1).rc_negative_window(random.Random(0), rec, 100, where=(40, 40))
        self.assertEqual(t["window"], w)

    def test_region_refused(self):
        rec, w = synthetic_record()
        rec.region = [1, len(rec.seq)]
        with self.assertRaises(ValueError):
            ml.reverse_kmers_in_record(rec, w)

    @unittest.skipUnless(have_cache(), "no record cache")
    def test_built_reverse_negatives(self):
        s = by_id(built({"R", "T"}))
        recs = dict((r.accession, r) for r in truth.load_with_sequences(fetch=False))
        for sid in ("R-rc-0001", "R-rc-0002", "T-strand-0001"):
            t = s[sid]["truth"]
            self.assertEqual(ml.reverse_kmers_in_record(recs[t["record"]], t["window"]), [], sid)
        # R-rc's second genome has many RefSeq records (its reverse k-mers are in the graph through them)
        self.assertEqual(s["R-rc-0002"]["truth"]["genome"], "eco")
        # T-strand: either seed refusal (both follow from the record), which one is the index's
        q = s["T-strand-0001"]["request"]["seeds"][0]["sequence"]
        self.assertEqual(s["T-strand-0001"]["expect"], [{"check": "status", "args": {
            "http": 400, "text_prefix": "seed '%s': " % q[:32]}}])
        self.assertEqual(len(s["T-strand-0001"]["truth"]["seed_refusals"]), 2)
        for text in s["T-strand-0001"]["truth"]["seed_refusals"]:
            ans = checks.normalise_direct(400, {"error": text + "o30"})
            self.assertEqual(checks.scenario_outcome(checks.evaluate(s["T-strand-0001"], ans, "direct-local")), "pass")
        walked = checks.normalise_direct(200, {"results": []})
        self.assertEqual(checks.scenario_outcome(checks.evaluate(s["T-strand-0001"], walked, "direct-local")), "fail")


class ChanceBounds(unittest.TestCase):
    """VERIFY-D: every global negative control states its chance bound; the builder refuses one above 1e-4."""

    def test_bounds(self):
        pk = areas.REFSEQ33M_KMERS / 4.0 ** 31
        q = "ACGT" * 15
        pat = {"route": "POST /pattern", "truth": {"random": q}, "request": ml.pattern_request([{"dna": q[:31]}])}
        self.assertAlmostEqual(areas.chance_bound(pat)["p"] / (2 * pk), 1.0, places=2)
        fwd = dict(pat, request=ml.pattern_request([{"dna": q[:31]}], strands="forward"))
        self.assertAlmostEqual(areas.chance_bound(fwd)["p"] / pk, 1.0, places=2)
        res = {"route": "POST /resolve", "truth": {"random": True}, "request": ml.resolve_request(q, labels=["1"])}
        self.assertAlmostEqual(areas.chance_bound(res)["p"] / (30 * pk), 1.0, places=2)
        trav = {"route": "POST /traverse", "truth": {"random": q}, "request": ml.traverse_request([{"sequence": q}])}
        self.assertAlmostEqual(areas.chance_bound(trav)["p"] / pk, 1.0, places=2)
        shuf = dict(pat, truth={"shuffled_from": {}})
        self.assertGreaterEqual(areas.chance_bound(shuf)["p_composition_worst"], areas.chance_bound(shuf)["p"])
        self.assertIsNone(areas.chance_bound(dict(pat, truth={"record": "x"})))
        # a 10-residue peptide (30 bases, below k) of 2-codon residues: 1,024 encodings, 2 strands, 2 offsets
        pep = {"route": "POST /pattern", "truth": {"random_peptide": "C" * 10},
               "request": ml.pattern_request([{"protein": "C" * 10}], genetic_code=1)}
        want = 2 * 1024 * 2 * areas.REFSEQ33M_KMERS / 4.0 ** 30
        self.assertAlmostEqual(areas.chance_bound(pep)["p"] / want, 1.0, places=2)
        self.assertGreater(areas.chance_bound(pep)["p"], areas.CHANCE_BOUND_MAX)

    def test_validate_refuses_a_loose_control(self):
        import build_manifest as bm
        sc = ml.scenario("N-x-0001", "N", "a loose control", "POST /pattern",
                         ml.pattern_request([{"id": "rp", "protein": "C" * 10}], genetic_code=1),
                         {"random_peptide": "C" * 10}, [{"check": "status", "args": {"http": 200}}])
        areas.state_chance_bounds([sc])
        self.assertTrue(any("chance bound" in p for p in bm.validate([sc])))

    @unittest.skipUnless(have_cache(), "no record cache")
    def test_built_controls(self):
        m = built({"N"})
        s = by_id(m)
        self.assertLessEqual(s["N-negative-0018"]["truth"]["chance_bound"]["p"], 1e-5)
        self.assertEqual(len(s["N-negative-0018"]["request"]["patterns"][0]["protein"]), 12)
        for sc in m["scenarios"]:
            t = sc.get("truth") or {}
            if any(k in t for k in ("random", "shuffled_from", "random_peptide")):
                self.assertLessEqual(t["chance_bound"]["p"], areas.CHANCE_BOUND_MAX, sc["id"])


class IndexEmpirical(unittest.TestCase):
    """VERIFY-D: the claims only the index's answers support are tagged with the index they were observed on."""

    @unittest.skipUnless(have_cache(), "no record cache")
    def test_tags(self):
        m = built({"R", "T", "X"})
        with open(os.path.join(SUITE, "truth", "staging_capabilities", "traverse_capabilities.json")) as f:
            fp = json.load(f)["index_fp"]
        tagged = [s for s in m["scenarios"] if "index-empirical" in s["tags"]]
        self.assertEqual(sorted(set(s["id"].rsplit("-", 1)[0] for s in tagged)), ["R-junction", "T-strand", "X-junction"])
        for s in tagged:
            self.assertEqual(s["truth"]["index_empirical"]["index_fp"], fp, s["id"])
            self.assertTrue(s["truth"]["index_empirical"]["claim"], s["id"])


class FlankCut(unittest.TestCase):
    """triage-v4 V03: a truth flank ends before a base outside ACGT (no k-mer of the graph holds it), stated as a
    structural end; the exact walk minimum stops there too."""

    def rec(self, seq):
        return types.SimpleNamespace(key="r", seq=seq, offset=1000, accession="R.1", column="1", seq_id=0,
                                     length=len(seq) + 1000, region=None)

    def test_flank_cut_both_sides(self):
        rng = random.Random(3)
        left = "".join(rng.choice("ACGT") for _ in range(50))
        right = "".join(rng.choice("ACGT") for _ in range(50))
        seq = left[:10] + "S" + left[11:] + "ACGTACGTAC" + right[:30] + "N" + right[31:]
        ctx = ml.Ctx([], 1)
        t = ctx.truth_at(self.rec(seq), 50, 10, flank=40)
        self.assertEqual(t["flank_start"], 1000 + 11)
        self.assertEqual(t["flank"], seq[11:90])
        self.assertEqual(t["flank_cut"], {"left": {"at": 1010, "base": "S"}, "right": {"at": 1090, "base": "N"}})
        w = ctx.flank_for(self.rec(seq), t, 45)
        self.assertEqual((w["flank_start"], len(w["flank"])), (1011, 79))
        clean = ctx.truth_at(self.rec(left + "ACGTACGTAC" + right), 50, 10, flank=30)
        self.assertNotIn("flank_cut", clean)
        self.assertNotIn("flank_cut", ctx.flank_for(self.rec(left + "ACGTACGTAC" + right), t, 20))

    def test_walk_minimum_stops_at_the_base(self):
        rng = random.Random(5)
        seq = "".join(rng.choice("ACGT") for _ in range(400))
        seq = seq[:300] + "S" + seq[301:]
        r = ml.walk_expectations(seq, [(100, 40, ["right"], 250)])[0]["right"]
        self.assertEqual((r["first"], r["min_bp"]), ("non_acgt", 160))
        # a k-mer holding the base is no successor: the walk's context at 200 met again elsewhere, followed by an N,
        # leaves the context one successor (counted as a second one it would read as ambiguous at 60 bp)
        ctx30 = seq[170:200]
        seq2 = seq[:300] + "A" + seq[301:] + "ACGT" + ctx30 + "N" + "TTGCA"
        r2 = ml.walk_expectations(seq2, [(100, 40, ["right"], 150)])[0]["right"]
        self.assertEqual((r2["first"], r2["min_bp"]), ("max_extension", 150))
        seq3 = seq[:300] + "A" + seq[301:] + "ACGT" + ctx30 + {"A": "C"}.get(seq[200], "A") + "TTGCA"
        r3 = ml.walk_expectations(seq3, [(100, 40, ["right"], 150)])[0]["right"]
        self.assertEqual((r3["first"], r3["min_bp"]), ("ambiguous", 60))

    @unittest.skipUnless(have_cache(), "no record cache")
    def test_t_walk_0018(self):
        t = by_id(built({"R", "T"}))["T-walk-0018"]["truth"]
        self.assertEqual(t["flank_cut"], {"left": {"at": 70345, "base": "S"}})
        self.assertEqual(t["start"] - t["flank_start"], 4898)


class DerivedLabelKinds(unittest.TestCase):
    """triage-v4 V02 (VERIFY-V4's amendment): the header-kind override only where the v4 answers counted more carriers
    than the default 1,000 (eco 6,181; cal 1), with the memory the 6,181 labels need; cal's pair walks via MCP."""

    @unittest.skipUnless(have_cache(), "no record cache")
    def test_pairs(self):
        import compare as cmp
        import gate_policy
        s = by_id(built({"R", "T"}))
        eco_h, cal_h = s["T-labels-0004"], s["T-labels-0006"]
        self.assertEqual(eco_h["request"]["strategy"]["labels"]["max_seed_labels"], ml.SERVER_MAX_SEED_LABELS)
        self.assertEqual(eco_h["request"]["strategy"]["bounds"]["max_memory_mb"], 4096)
        self.assertEqual(eco_h["truth"]["derived_labels_v4"], 6181)
        self.assertNotIn("max_seed_labels", cal_h["request"]["strategy"]["labels"])
        self.assertEqual(cal_h["request"]["strategy"]["bounds"]["max_memory_mb"], 512)
        cases, _ = cmp.load_by_design(cmp.DEFAULT_BY_DESIGN)
        self.assertIsNone(gate_policy.mcp_exclusion(cal_h, cases))                    # basic path: the MCP walks it
        self.assertIn("advanced_path", gate_policy.mcp_exclusion(eco_h, cases)["rules"])   # stays not implemented


class WorkBudget(unittest.TestCase):
    """triage-v4 G18 (VERIFY-V4: the seed phase's cost is readable only from an attempt's usage): T-budget-0003's work
    budget lies above its seed's measured seed phase, so the walk stops partway on every arm, never the seed."""

    @unittest.skipUnless(have_cache(), "no record cache")
    def test_seed_and_budget(self):
        s = by_id(built({"R", "T"}))
        b, att = s["T-budget-0003"], s["T-attempts-0002"]
        self.assertEqual(b["request"]["seeds"], [{k: v for k, v in att["request"]["seeds"][0].items()}])
        budget = b["request"]["strategy"]["bounds"]["max_work_units"]
        self.assertGreaterEqual(budget, 2 * b["truth"]["seed_phase"]["work_seed"])
        stop = [e for e in b["expect"] if e["check"] == "traverse_stop"][0]["args"]
        self.assertEqual(stop["walks"], ["partial"])
        failed = checks.normalise_direct(200, {"results": [{"outcome": {"walks": "failed"}, "resource_stop": {
            "resource": "work", "requested": budget}}]})
        res = dict((r["check"], r["outcome"]) for r in checks.evaluate(b, failed, "direct-local"))
        self.assertEqual(res["traverse_stop"], "fail")
        # the slot keeps its place and every other T-budget scenario its window
        ids = [x["id"] for x in built({"R", "T"})["scenarios"] if x["id"].startswith("T-budget")]
        self.assertEqual(ids, ["T-budget-0001", "T-budget-0002", "T-budget-0003", "T-budget-0004"])


class ResolveDecisive(unittest.TestCase):
    """Item 5 (resolve) and triage-v4 V05: the explicit-label budget is the host's cap; the budget stop follows its
    complete twin, a warm-up, and is a budget_stress scenario; eco's reference is a query this host completes."""

    def tearDown(self):
        ml.set_host_capabilities(None)

    def test_budget_and_naming(self):
        self.assertEqual(ml.RESOLVE_TIME_BUDGET_MS, 60000)
        ml.set_host_capabilities({}, {"resolve": {"time_budget": {"max_time_ms": 60000}}})
        r = lambda tb: ml.resolve_by_design(ml.resolve_request("A" * 100, labels=["1"], time_budget_ms=tb))  # noqa
        self.assertNotIn("resolve.time_budget_above_host_max", r(60000))
        self.assertIn("resolve.time_budget_above_host_max", r(70000))

    def test_builder_refuses_a_lower_cap(self):
        import build_manifest as bm
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            for name in bm.CAPABILITY_DOCUMENTS:
                with open(os.path.join(SUITE, "truth", "staging_capabilities", name + ".json")) as f:
                    doc = json.load(f)
                if name == "traverse_capabilities":
                    doc["resolve"]["time_budget"]["max_time_ms"] = 20000
                with open(os.path.join(d, name + ".json"), "w") as f:
                    json.dump(doc, f)
            with self.assertRaises(SystemExit):
                bm.build(1, fetch=False, only_areas={"C"}, capabilities_dir=d)

    def test_place_warm_ups(self):
        import build_manifest as bm
        sc = lambda i, tags=(), warms=None: {"id": i, "tags": list(tags),  # noqa: E731
                                              "truth": {"warm_up_for": warms} if warms else {}}
        got = bm.place_warm_ups([sc("A"), sc("B"), sc("C"), sc("W1", ["warm-up"], "B"), sc("W2", ["warm-up"], "B"),
                                 sc("W3", ["warm-up"], "Z"), sc("X", (), "A")])
        self.assertEqual([s["id"] for s in got], ["A", "W1", "W2", "B", "C", "W3", "X"])

    @unittest.skipUnless(have_cache(), "no record cache")
    def test_built(self):
        m = built()
        s = by_id(m)
        order = [x["id"] for x in m["scenarios"]]
        for twin, sid in (("R-twin-0002", "R-budget-0001"), ("R-twin-0004", "R-budget-0002")):
            self.assertEqual(order.index(twin) + 1, order.index(sid))                    # the warm-up right before
            self.assertIn("warm-up", s[twin]["tags"])
            self.assertEqual(s[twin]["request"]["sequence"], s[sid]["request"]["sequence"])
            self.assertEqual(s[twin]["request"]["bounds"]["time_budget_ms"], 60000)
            b = s[sid]
            self.assertIn("budget_stress", b["tags"])
            budget = b["request"]["bounds"]["time_budget_ms"]
            self.assertLess(budget - 250, b["truth"]["budget_stress"]["warm_work_ms"])  # the work deadline in the warm work
            self.assertEqual(b["repeat"], 2)
        self.assertEqual(len(s["R-budget-0001"]["request"]["sequence"]), 5000)          # eco's second half
        for x in m["scenarios"]:
            if x["route"] == "POST /resolve" and isinstance(x["request"], dict) and x["request"].get("labels") \
                    and "budget" not in x["tags"] and "edge" not in x["tags"]:
                self.assertEqual(x["request"]["bounds"]["time_budget_ms"], 60000, x["id"])
        rels = dict((r["id"], r) for r in m["relations"])
        self.assertEqual(rels["REL-prefix_equals-0002"]["args"]["a_from_kmer"], 2500)


class WalksDecisive(unittest.TestCase):
    """triage-v4 V06 and item 5 on the traverse route: less work where the cold first call under the load could not
    finish within the MCP's 15 s preset budget, the host's cap for the long column walks, warm-ups before the walks
    whose cold seed or row reads dominated; every window drawn as in v4."""

    @unittest.skipUnless(have_cache(), "no record cache")
    def test_walks(self):
        m = built()
        s = by_id(m)
        with open(os.path.join(SUITE, "manifest", "suite_v4.json")) as f:
            v4 = by_id(json.load(f))
        order = [x["id"] for x in m["scenarios"]]
        bounds = lambda sid: s[sid]["request"]["strategy"]["bounds"]  # noqa: E731
        for n in range(1, 37):
            sid = "T-walk-%04d" % n
            self.assertEqual(s[sid]["truth"]["start"], v4[sid]["truth"]["start"], sid)    # the same windows
            old = v4[sid]["request"]["strategy"]["bounds"]["max_extension_bp"]
            want = {6: 2500, 12: 1000, 18: 10000, 36: 2500}.get(n, 300 if old == 1000 else old)
            self.assertEqual(bounds(sid)["max_extension_bp"], want, sid)
            self.assertEqual(bounds(sid)["time_budget_ms"], 60000 if n in (6, 12, 18, 36) else 15000, sid)
            self.assertEqual(s[sid]["timing_class"], v4[sid]["timing_class"], sid)
        # the MCP refuses above its walk preset's budget row, by design (the long walks are the direct arms')
        self.assertIn("traverse.bounds_budget_above_preset_row", s["T-walk-0006"]["mcp"]["by_design"])
        for n in range(1, 5):
            sid = "T-repeat-%04d" % n
            self.assertEqual(bounds(sid)["max_extension_bp"], 400, sid)
            w = s[order[order.index(sid) - 1]]
            self.assertIn("warm-up", w["tags"], sid)
            self.assertEqual(w["truth"]["warm_up_for"], sid)
            self.assertIn(s[sid]["request"]["seeds"][0]["sequence"], w["request"]["sequence"])
            self.assertEqual(w["request"]["labels"], s[sid]["request"]["seeds"][0]["labels"])
        w = s[order[order.index("T-labels-0003") - 1]]
        self.assertEqual((w["truth"]["warm_up_for"], w["request"]["labels"]), ("T-labels-0003", ["511145"]))
        self.assertEqual(bounds("T-labels-0003")["max_extension_bp"], 50)
        self.assertEqual(bounds("T-labels-0005")["max_extension_bp"], 200)
        for n in range(1, 4):
            self.assertEqual(bounds("T-coords-%04d" % n)["max_extension_bp"], 250)
        w = s[order[order.index("N-hard-0022") - 1]]
        self.assertEqual(w["request"]["sequence"], "TTAGGG" * 10)
        self.assertEqual(w["request"]["labels"], s["N-hard-0022"]["request"]["seeds"][0]["labels"])
        for x in m["scenarios"]:
            if x["id"].startswith("R-warm"):     # the R-twin warm-ups are 5-10 kb, the direct arms' (by design)
                self.assertIn("warm-up", x["tags"])
                self.assertLessEqual(len(x["request"]["sequence"]), ml.RESOLVE_MAX_QUERY_BP, x["id"])
                self.assertIn("tool", x["mcp"], x["id"])
            self.assertFalse(x["id"].startswith("D-R-warm"), x["id"])


class TitlesTested(unittest.TestCase):
    """Finding 6A (VERIFY-C) and P5-v4-content-04: the record's own path on '+' for the IUPAC long pattern, the '-' path
    for the reverse one, and both primers of the forward-only pair counted."""

    @unittest.skipUnless(have_cache(), "no record cache")
    def test_checks(self):
        s = by_id(built({"L", "N"}))
        ex6 = [e for e in s["L-long-0006"]["expect"] if e["check"] == "path_spells"]
        self.assertEqual(ex6, [{"check": "path_spells", "args": {"strand": "+"}}])
        self.assertIn({"check": "path_spells", "args": {"strand": "-"}}, s["L-long-0005"]["expect"])
        t = s["L-long-0006"]["truth"]
        entry = {"id": "w", "results": [{"sequence": "A" * 60, "strand": "-", "offset": 0}], "retrieval_complete": True,
                 "counts": {"paths": {"relation": "exact", "value": 1}}}
        ans = checks.normalise_direct(200, {"patterns": [entry]})
        self.assertEqual(checks.CHECKS["path_spells"](s["L-long-0006"], ans, {"strand": "+"}, "direct-local")[0], "fail")
        entry["results"].append({"sequence": t["window"], "strand": "+", "offset": 0})
        ans = checks.normalise_direct(200, {"patterns": [entry]})
        self.assertEqual(checks.CHECKS["path_spells"](s["L-long-0006"], ans, {"strand": "+"}, "direct-local")[0], "pass")
        pair = s["N-usecase-0005"]
        got = sorted((e["check"], e["args"].get("pattern")) for e in pair["expect"] if e["args"].get("pattern") is not None)
        self.assertEqual(got, [("count_relation", 0), ("count_relation", 1), ("entry_fields", 0), ("entry_fields", 1)])


if __name__ == "__main__":
    unittest.main()
