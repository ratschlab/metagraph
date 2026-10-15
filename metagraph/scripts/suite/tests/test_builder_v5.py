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
        self.assertEqual(sorted(set(s["id"].rsplit("-", 1)[0] for s in tagged)),
                         ["R-discover", "R-junction", "T-strand", "X-junction"])    # areas R, T, X built
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
        sc = lambda i, tags=(), warms=None: dict({"id": i, "tags": list(tags), "truth": {}},  # noqa: E731
                                                 **({"warms": warms} if warms else {}))
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
            budget = b["request"]["bounds"]["time_budget_ms"]
            self.assertLess(budget - 250, b["truth"]["stop_design"]["warm_work_ms"])   # the work deadline in the warm work
            lo, hi = b["truth"]["stop_design"]["stop_window_warm_ms"]
            self.assertTrue(lo <= b["truth"]["stop_design"]["warm_work_ms"] <= hi, sid)   # the measured warm work inside
            # in the gate, required on direct-local (the owner, 2026-10-15): one repeat, no budget_stress
            self.assertEqual(b["repeat"], 1)
            self.assertNotIn("budget_stress", b["tags"])
        self.assertEqual(len(s["R-budget-0001"]["request"]["sequence"]), 5000)          # eco's second half
        # the owner, 2026-10-15: 480 ms (the batch model's stop window 0.28-2.48 s of warm work, against 1.13-1.75 s
        # at 700 ms), pae's 420 ms
        self.assertEqual([s[x]["request"]["bounds"]["time_budget_ms"] for x in ("R-budget-0001", "R-budget-0002")],
                         [480, 420])
        self.assertEqual(s["R-budget-0001"]["truth"]["stop_design"]["stop_window_warm_ms"], [284, 2480])
        self.assertIn("480 ms budget", s["R-budget-0001"]["title"])
        # one execution each: off the paired plan's warm-repeat index (the re-verification: -0002 sat at 50)
        for sid in ("R-budget-0001", "R-budget-0002"):
            self.assertNotEqual(order.index(sid) % 10, 0, sid)
        for x in m["scenarios"]:
            if x["route"] == "POST /resolve" and isinstance(x["request"], dict) and x["request"].get("labels") \
                    and "budget" not in x["tags"] and "edge" not in x["tags"]:
                self.assertEqual(x["request"]["bounds"]["time_budget_ms"], 60000, x["id"])
        rels = dict((r["id"], r) for r in m["relations"])
        self.assertEqual(rels["REL-prefix_equals-0002"]["args"]["a_from_kmer"], 2500)


class WalksDecisive(unittest.TestCase):
    """triage-v4 V06 and item 5 on the traverse route, as the owner amended them (2026-10-15): the three header walks a
    cold first call left inconclusive walk 300 bp, the other nine 1 kb; each long column walk is split into the gate's
    walk (the host's cap, a radius the slowest cold rate covers in half of it) and a long budget_stress twin right
    after it; warm-ups before the walks whose cold seed or row reads dominated; every window drawn as in v4."""

    @unittest.skipUnless(have_cache(), "no record cache")
    def test_walks(self):
        m = built()
        s = by_id(m)
        with open(os.path.join(SUITE, "manifest", "suite_v4.json")) as f:
            v4 = by_id(json.load(f))
        order = [x["id"] for x in m["scenarios"]]
        bounds = lambda sid: s[sid]["request"]["strategy"]["bounds"]  # noqa: E731
        short = {6: 1000, 12: 500, 18: 5000, 36: 1000}     # cold first-call frontiers at about 30 s (VERIFY-V5A)
        long_twin = {6: 10000, 12: 10000, 18: 10000, 36: 10000}   # the v4 test's radius (the owner, 2026-10-15)
        for n in range(1, 37):
            sid = "T-walk-%04d" % n
            self.assertEqual(s[sid]["truth"]["start"], v4[sid]["truth"]["start"], sid)    # the same windows
            old = v4[sid]["request"]["strategy"]["bounds"]["max_extension_bp"]
            want = short.get(n, 300 if n in (4, 14, 22) else old)
            self.assertEqual(bounds(sid)["max_extension_bp"], want, sid)
            self.assertEqual(bounds(sid)["time_budget_ms"], 60000 if n in short else 15000, sid)
            self.assertEqual(s[sid]["timing_class"], v4[sid]["timing_class"], sid)
            self.assertNotIn("budget_stress", s[sid]["tags"], sid)
            if n in short:
                tw = s[order[order.index(sid) + 1]]
                self.assertTrue(tw["id"].startswith("T-longwalk-"), sid)
                self.assertEqual(tw["request"]["seeds"], s[sid]["request"]["seeds"])
                self.assertEqual(tw["request"]["strategy"]["bounds"]["max_extension_bp"], long_twin[n])
                self.assertEqual(tw["request"]["strategy"]["bounds"]["max_extension_bp"],
                                 v4[sid]["request"]["strategy"]["bounds"]["max_extension_bp"])
                self.assertEqual(tw["request"]["strategy"]["bounds"]["time_budget_ms"], 60000)
                self.assertEqual(tw["truth"]["start"], v4[sid]["truth"]["start"])
                self.assertIn("budget_stress", tw["tags"])
                self.assertEqual(tw["truth"]["budget_stress"]["short_twin"], sid)
                self.assertIn(tw["id"], s[sid]["title"])
        self.assertEqual(sorted(x["id"] for x in m["scenarios"] if "budget_stress" in x["tags"]),
                         ["T-longwalk-%04d" % n for n in range(1, 5)])
        # the MCP refuses above its walk preset's budget row, by design (the long walks are the direct arms')
        self.assertIn("traverse.bounds_budget_above_preset_row", s["T-walk-0006"]["mcp"]["by_design"])
        for n in range(1, 5):
            sid = "T-repeat-%04d" % n
            self.assertEqual(bounds(sid)["max_extension_bp"], 400, sid)
            w = s[order[order.index(sid) - (2 if sid in areas.SPLIT_WARM_IDS else 1)]]
            self.assertIn("warm-up", w["tags"], sid)
            self.assertEqual(w["warms"], sid)
            self.assertIn(s[sid]["request"]["seeds"][0]["sequence"], w["request"]["sequence"])
            self.assertEqual(w["request"]["labels"], s[sid]["request"]["seeds"][0]["labels"])
        # VERIFY-V5A: a warm-up reads what its walk must reach (the seed and each arm's exact minimum plus the margin),
        # so the 6-copy B. subtilis region fits a short cold resolve; T-walk-0009's cold seed row gets one too
        for n in range(1, 5):
            sid = "T-repeat-%04d" % n
            t = s[sid]["truth"]
            radius = bounds(sid)["max_extension_bp"]
            want = [min(radius, t["walk"][side]["min_bp"] + areas.WARM_MARGIN_BP) for side in ("left", "right")]
            if sid in areas.SPLIT_WARM_IDS:
                continue
            w = s[order[order.index(sid) - 1]]
            self.assertEqual(len(w["request"]["sequence"]), want[0] + t["length"] + want[1], sid)
        # the re-verification: T-repeat-0002's 295 bp split in two warm-ups, each its own 60,000 ms resolve: the seed
        # with its left 46 bp, then the seed's last 30 bp with its right 149 bp, both right before the walk
        t = s["T-repeat-0002"]["truth"]
        seed = s["T-repeat-0002"]["request"]["seeds"][0]["sequence"]
        a, b = s["R-warm-0002"], s["R-warm-0008"]
        self.assertEqual(order[order.index("T-repeat-0002") - 2:order.index("T-repeat-0002")], ["R-warm-0002", "R-warm-0008"])
        self.assertEqual((len(a["request"]["sequence"]), len(b["request"]["sequence"])), (46 + 100, 30 + 149))
        self.assertTrue(a["request"]["sequence"].endswith(seed))
        self.assertTrue(b["request"]["sequence"].startswith(seed[-30:]))
        self.assertEqual((a["truth"]["start"], b["truth"]["start"]), (t["start"] - 46, t["start"] + 70))
        for w in (a, b):
            self.assertEqual((w["warms"], w["request"]["bounds"]["time_budget_ms"]), ("T-repeat-0002", 60000))
        w9 = s[order[order.index("T-walk-0009") - 1]]
        self.assertEqual((w9["id"], w9["request"]["labels"], len(w9["request"]["sequence"])),
                         ("R-warm-0007", ["10710"], 31 + 2 * 100))
        self.assertIn(s["T-walk-0009"]["request"]["seeds"][0]["sequence"], w9["request"]["sequence"])
        w = s[order[order.index("T-labels-0003") - 1]]
        self.assertEqual((w["warms"], w["request"]["labels"]), ("T-labels-0003", ["511145"]))
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


class GateCoverage(unittest.TestCase):
    """VERIFY-D's gaps of the resolve / traverse gate the MCP can express, and P5-v4-content-06's: forced per-arm caps,
    a cap above the service's maximum, the permit side of min_live_labels, a dropped seed label with a client seed_id."""

    def test_trie_heads(self):
        rng = random.Random(9)
        rnd = lambda n: "".join(rng.choice("ACGT") for _ in range(n))  # noqa: E731
        w = rnd(31)
        # three copies of w, each followed by 20 shared bases, then a base that differs per copy
        tail = rnd(20)
        seq = rnd(50) + w + tail + "A" + rnd(60) + w + tail + "C" + rnd(60) + w + tail + "G" + rnd(50)
        heads = ml.trie_heads(seq, 50, 31, "right", 40)
        self.assertEqual(heads[:20], [1] * 20)
        self.assertEqual(heads[20], 3)                         # the copies' continuations part at the 21st base

    @unittest.skipUnless(have_cache(), "no record cache")
    def test_scenarios(self):
        import compare as cmp
        import gate_policy
        m = built()
        s = by_id(m)
        cases, _ = cmp.load_by_design(cmp.DEFAULT_BY_DESIGN)
        recs = dict((r.key, r) for r in truth.load_with_sequences(fetch=False))
        for sid, knob in (("T-caps-0002", "max_live_paths"), ("T-caps-0003", "max_paths")):
            sc = s[sid]
            stop = [e["args"] for e in sc["expect"] if e["check"] == "traverse_stop"][0]
            self.assertEqual((stop["knob"], stop["cap_reason"], stop["walks"]), ("bounds." + knob, knob, ["partial"]))
            t = sc["truth"]
            rec = recs[t["genome"]]
            first = [x for x in (t["trie"]["left"]["first_depth_with_3_heads"],
                                 t["trie"]["right"]["first_depth_with_3_heads"]) if x is not None]
            self.assertTrue(first and min(first) < sc["request"]["strategy"]["bounds"]["max_extension_bp"])
            self.assertEqual(ml.trie_heads(rec.seq, t["start"], 31, "left", min(first))[-1] >= 3 or
                             ml.trie_heads(rec.seq, t["start"], 31, "right", min(first))[-1] >= 3, True)
            self.assertIsNone(gate_policy.mcp_exclusion(sc, cases), sid)          # the MCP can send it
        self.assertEqual(s["T-caps-0004"]["mcp"]["by_design"], ["traverse.bounds_knob_above_service_maximum"])
        self.assertEqual(s["T-seeds-0002"]["mcp"]["arguments"]["seeds"][0]["seed_id"], "suite-v5-client-seed-id")
        q = s["T-quorum-0002"]
        self.assertEqual(q["request"]["seeds"][0]["labels"], [q["truth"]["record"], q["truth"]["taxid"]])
        self.assertIn("walk", q["truth"])                                          # the exact minimum, from the record
        for sid in ("T-caps-0002", "T-caps-0003", "T-caps-0004", "T-quorum-0002", "T-seeds-0002"):
            self.assertIn("gate-coverage", s[sid]["tags"])
            self.assertNotIn("D-" + sid, s)


class RecordEndViaMcp(unittest.TestCase):
    """P4-builder-truth-01: the T-end walks run under the walk preset's budget, so the MCP walks them (v4's 60,000 ms
    made it refuse all six by design)."""

    @unittest.skipUnless(have_cache(), "no record cache")
    def test_t_end(self):
        import compare as cmp
        cases, _ = cmp.load_by_design(cmp.DEFAULT_BY_DESIGN)
        s = by_id(built({"R", "T"}))
        ends = [x for x in s.values() if x["id"].startswith("T-end")]
        self.assertEqual(len(ends), 6)
        for x in ends:
            self.assertEqual(x["request"]["strategy"]["bounds"]["time_budget_ms"], 15000, x["id"])
            self.assertEqual(ml.traverse_by_design(x["request"]), [], x["id"])
            self.assertEqual([c for c in x["mcp"].get("by_design") or []
                              if (cases.get(c) or {}).get("service_answer", {}).get("outcome") == "refused_at_submit"],
                             [], x["id"])


class DeterminismPool(unittest.TestCase):
    """P4-builder-truth-03: gen_D never twins a scenario whose answer depends on what ran before it (a run-unique,
    retained attempt id), one that repeats already or one whose status allows another answer than a 200; and V04:
    the attempt and its duplicate are two scenarios, each deciding all of its checks."""

    def sc(self, sid, area="T", **kw):
        base = {"id": sid, "area": area, "title": sid, "route": "POST /traverse", "request": {}, "truth": {},
                "expect": [{"check": "status", "args": {"http": 200}}], "timing_class": "fast", "repeat": 1,
                "tags": [], "mcp": {"not_implemented": "x"}}
        base.update(kw)
        return base

    def test_pool(self):
        pool = [self.sc("T-a", run_unique=["attempt_id"], tags=["attempts"]),
                self.sc("T-b", tags=["retention-dependent"]), self.sc("T-c", repeat=2),
                self.sc("T-d", expect=[{"check": "status", "args": {"http": [200, 409]}}]),
                self.sc("T-e"), self.sc("T-f")]
        picks = [d["derived_from"] for d in areas.gen_D(ml.Ctx([], 1), pool)]
        self.assertEqual(picks, ["T-e", "T-f"])

    @unittest.skipUnless(have_cache(), "no record cache")
    def test_attempt_pair(self):
        m = built()
        s = by_id(m)
        a, dup = s["T-attempts-0002"], s["T-attempts-0003"]
        self.assertEqual(a["request"], dup["request"])
        self.assertEqual((a["repeat"], dup["repeat"]), (1, 1))
        self.assertEqual(a["expect"][0]["args"]["http"], 200)
        self.assertEqual(dup["expect"][0]["args"]["http"], 409)
        self.assertEqual(a["run_unique"], dup["run_unique"])
        order = [x["id"] for x in m["scenarios"]]
        self.assertEqual(order.index("T-attempts-0003"), order.index("T-attempts-0002") + 1)
        body = {"error": "duplicate", "attempt": {"state": "finished", "usage": {}}}
        res = checks.evaluate(dup, checks.normalise_direct(409, body), "direct-local")
        self.assertEqual(checks.scenario_outcome(res), "pass")
        self.assertFalse([d for d in m["scenarios"] if d["area"] == "D" and d.get("run_unique")])


class RepeatReentry(unittest.TestCase):
    """P5-v4-content-02: the real-repeat walks of a short period assert the re-entry their titles claim."""

    @unittest.skipUnless(have_cache(), "no record cache")
    def test_reentry(self):
        s = by_id(built({"N"}))
        for sid in ("N-realrep-0003", "N-realrep-0006"):
            sc = s[sid]
            self.assertIn({"check": "blocked_event", "args": {"reason": "rejoined_seed", "max_at_bp": 0}},
                          sc["expect"], sid)
            arms = dict((side, {"segments": [{"id": 0, "events": []}]}) for side in ("left", "right"))
            ans = checks.normalise_direct(200, {"results": [{"arms": arms}]})
            self.assertEqual(checks.CHECKS["blocked_event"](sc, ans, {"reason": "rejoined_seed", "max_at_bp": 0},
                                                            "direct-local")[0], "fail")
            for side in arms:
                arms[side]["segments"][0]["events"] = [{"type": "blocked", "reason": "rejoined_seed", "at_bp": 0}]
            ans = checks.normalise_direct(200, {"results": [{"arms": arms}]})
            self.assertEqual(checks.CHECKS["blocked_event"](sc, ans, {"reason": "rejoined_seed", "max_at_bp": 0},
                                                            "direct-local")[0], "pass")


class HeaderExactColumnTolerant(unittest.TestCase):
    """P5-v4-content-03: exact runs on the records' headers (computed from the whole records), a tolerance on the
    taxid columns (their other records may carry the junction's or the variant's bases); and the mutated seed's
    traverse refusal as T-strand's: either seed refusal, both following from the record."""

    def test_record_runs(self):
        rec = types.SimpleNamespace(seq="C" * 10 + "ACGT" * 20 + "G" * 10)
        q = ("ACGT" * 20)[:40] + "T" * 31
        self.assertEqual(areas._record_runs(rec, q), [[0, 10]])

    @unittest.skipUnless(have_cache(), "no record cache")
    def test_junction(self):
        s = by_id(built({"R", "N"}))
        sc = s["R-junction-0001"]
        t = sc["truth"]
        args = dict((e["args"]["label"], e["args"]) for e in sc["expect"] if e["check"] == "resolve_runs")
        a_h, b_h, a_c, b_c = t["a"]["record"], t["b"]["record"], t["a"]["taxid"], t["b"]["taxid"]
        self.assertEqual(args[a_h]["runs"], [[0, t["a_run_end"]]])
        self.assertEqual(args[b_h]["runs"], [[t["b_run_start"], 170]])

        def answer(col_a_end):
            labs = [{"label": a_h, "runs": [[0, t["a_run_end"]]], "kmers_supported": t["a_run_end"]},
                    {"label": b_h, "runs": [[t["b_run_start"], 170]], "kmers_supported": 170 - t["b_run_start"]},
                    {"label": a_c, "runs": [[0, col_a_end]], "kmers_supported": col_a_end},
                    {"label": b_c, "runs": [[t["b_run_start"], 170]], "kmers_supported": 170 - t["b_run_start"]}]
            return checks.normalise_direct(200, {"num_kmers": 170, "graph_runs": [[0, 170]], "labels": labs})
        rr = checks.CHECKS["resolve_runs"]
        # another record of A's column carrying the junction's bases: the column's run reaches into the junction
        self.assertEqual(rr(sc, answer(85), args[a_c], "direct-local")[0], "pass")
        # ... but never past it into B's k-mers
        self.assertEqual(rr(sc, answer(120), args[a_c], "direct-local")[0], "fail")
        mut = s["N-nearmiss-0034"]
        self.assertEqual(sorted(mut["request"]["labels"]), sorted([mut["truth"]["taxid"], mut["truth"]["record"]]))
        self.assertFalse([e for e in mut["expect"] if e["args"].get("label") == mut["truth"]["taxid"]
                          and "max_kmers" in e["args"]])
        seed = s["N-nearmiss-0035"]
        self.assertEqual(seed["truth"]["variant_kmers_in_record"], 0)
        self.assertEqual(len(seed["truth"]["seed_refusals"]), 2)


class NothingDropsSilently(unittest.TestCase):
    """P5-v4-content-05: every planned wrong-genetic-code scenario that is not built is listed in `skipped`."""

    @unittest.skipUnless(have_cache(), "no record cache")
    def test_skipped(self):
        m = built({"X", "M", "N"})
        skipped = [x["scenario"] for x in m["skipped"]]
        for table in (2, 5, 12):
            self.assertIn("X-table1 table %d" % table, skipped)
        self.assertIn("M-protein table 12", skipped)
        self.assertEqual(len([x for x in skipped if x.startswith("N-usecase own table")]), 2)
        built_tables = set(int(x["title"].split("(table ")[1].split(")")[0]) for x in m["scenarios"]
                           if x["id"].startswith("X-table1"))
        prots = [p for p in truth.load_proteins() if p["genome"] in dict((r.key, r) for r in
                                                                      truth.load_records())]
        for table in set(p["transl_table"] for p in prots) - {1, 11}:
            self.assertTrue(table in built_tables or "X-table1 table %d" % table in skipped, table)


class Discoveries(unittest.TestCase):
    """Item 8 (lane S-checks' hand-off) and the async session's constraint: R-discover-0001/-0004 ask above v4's
    full-support tie, so the source is listed; every discovery names max_labels (at most 100 but for the by-design
    case) and kind, on both arms."""

    def test_validation(self):
        import build_manifest as bm
        mk = lambda d: ml.scenario("R-x-0001", "R", "t", "POST /resolve", ml.resolve_request("A" * 100, discover=d),  # noqa
                                   {}, [{"check": "status", "args": {"http": 200}}])
        self.assertTrue(any("max_labels and kind" in p for p in bm.validate([mk({"kind": "column"})])))
        self.assertTrue(any("max_labels and kind" in p for p in bm.validate([mk({"max_labels": 5})])))
        self.assertEqual([p for p in bm.validate([mk({"max_labels": 5, "kind": "column"})]) if "discover" in p], [])
        self.assertEqual([p for p in bm.validate([mk({"max_labels": 150, "kind": "column"})]) if "discover" in p], [])
        bad = mk({"max_labels": 1001, "kind": "column"})
        self.assertEqual([p for p in bm.validate([bad]) if "discover" in p], [])     # it names the by-design pair
        bad["mcp"]["by_design"] = ["resolve.job_path"]
        self.assertTrue(any("above the MCP's 1000" in p for p in bm.validate([bad])))
        # the service's cap of 1,000 (v5 start): 150 names no case, 1,001 the host-answers case, 0 the refusal
        self.assertEqual([c for c in mk({"max_labels": 150, "kind": "column"})["mcp"]["by_design"]
                          if "discover_limit" in c], [])
        self.assertEqual([c for c in mk({"max_labels": 1001, "kind": "column"})["mcp"]["by_design"]
                          if "discover_limit" in c], ["resolve.discover_limit_above_cap_host_answers"])
        self.assertEqual([c for c in mk({"max_labels": 0, "kind": "column"})["mcp"]["by_design"]
                          if "discover_limit" in c], ["resolve.discover_limit_above_cap"])

    @unittest.skipUnless(have_cache(), "no record cache")
    def test_limits(self):
        s = by_id(built({"R"}))
        e23, e25 = s["R-edge-0023"], s["R-edge-0025"]
        self.assertIn({"check": "resolve_discover_contains", "args": {"label": "511145"}}, e23["expect"])
        self.assertNotIn("split-case", e23["tags"])
        self.assertEqual(e25["request"]["discover"]["max_labels"], 1001)
        self.assertEqual([c for c in e25["mcp"]["by_design"] if "discover_limit" in c],
                         ["resolve.discover_limit_above_cap_host_answers"])
        self.assertEqual([c for c in s["R-edge-0022"]["mcp"]["by_design"] if "discover_limit" in c],
                         ["resolve.discover_limit_above_cap"])
        for sid in ("R-discover-0001", "R-discover-0004"):
            d = s[sid]["request"]["discover"]
            self.assertEqual(d["max_labels"], s[sid]["truth"]["discover_v4"]["full_support"] + 1, sid)
            self.assertLessEqual(d["max_labels"], 100)
            self.assertNotIn("unless truncated", s[sid]["title"])
            self.assertIn("index-empirical", s[sid]["tags"])                  # the tie is the index's (VERIFY-V5A)
            self.assertTrue(s[sid]["truth"]["index_empirical"]["claim"])


class McpExposure(unittest.TestCase):
    """Lane S-verdict's hand-off (item 9): mcp.exposes from the v4 baseline holds only the checks decided on EVERY
    recorded MCP execution of the scenario, and only where the MCP request is the baseline's."""

    def test_baseline_and_apply(self):
        import build_manifest as bm
        mcp = {"tool": "traverse_resolve", "arguments": {"labels": ["1"], "time_budget_ms": 20000}}
        man = {"scenarios": [{"id": "A", "mcp": mcp}, {"id": "B", "mcp": mcp}, {"id": "C", "mcp": mcp}]}
        line = lambda sid, rep, **o: {"arm": "mcp", "scenario": sid, "repeat": rep, "checks": [  # noqa: E731
            {"check": c, "outcome": v} for c, v in o.items()]}
        rows = [line("A", 0, status="pass", invariants="pass", resolve_runs="pass"),
                line("A", 1, status="pass", invariants="pass", resolve_runs="inconclusive", deterministic="pass"),
                line("B", 0, status="fail", resolve_runs="fail"),            # a relayed 503: nothing decided
                line("B", 1, status="pass", resolve_runs="pass"),
                line("C", 0, status="pass", resolve_runs="not_exposed", invariants="fail"),
                dict(line("A", 0, status="pass"), arm="direct-local")]
        b = bm.exposure_baseline(rows, man, {})
        got = dict((k, v["exposes"]) for k, v in b["scenarios"].items())
        self.assertEqual(got, {"A": ["invariants", "status"], "B": [], "C": ["invariants", "status"]})
        exp = [{"check": "status", "args": {}}, {"check": "invariants", "args": {}}]
        same = {"id": "A", "expect": exp, "mcp": {"tool": "traverse_resolve",
                                                  "arguments": {"labels": ["1"], "time_budget_ms": 60000}}}
        other = {"id": "C", "expect": exp[:1], "mcp": {"tool": "traverse_resolve", "arguments": {"labels": ["2"]}}}
        self.assertEqual(bm.apply_exposure([same, other], b), 1)
        self.assertEqual(same["mcp"]["exposes"], ["invariants", "status"])   # a resolve's budget aside
        self.assertNotIn("exposes", other["mcp"])                            # another request: not the baseline's
        walk = {"tool": "traverse_fetch", "arguments": {"overrides": {"bounds": {"time_budget_ms": 15000}}}}
        self.assertNotEqual(bm.mcp_request_sha256(walk), bm.mcp_request_sha256(
            {"tool": "traverse_fetch", "arguments": {"overrides": {"bounds": {"time_budget_ms": 60000}}}}))

    @unittest.skipUnless(have_cache(), "no record cache")
    def test_built(self):
        m = built()
        s = by_id(m)
        self.assertIn("resolve_runs", s["R-pos-0001"]["mcp"]["exposes"])
        self.assertNotIn("exposes", s["T-walk-0006"]["mcp"])          # its request changed
        self.assertNotIn("exposes", s["T-end-0001"]["mcp"])           # v4's MCP refused it: nothing decided
        self.assertNotIn("walk_in_record", s["T-walk-0004"]["mcp"].get("exposes", []))   # inconclusive on r0
        self.assertEqual(m["mcp_exposure_baseline"]["source"]["run"], ["20261014T061053Z_refseq97-6f945db4"])


class BudgetStressRule(unittest.TestCase):
    """gate_policy derives the budget_stress exclusions from the scenario tag (the owner approved it for v5's long
    walks on 2026-10-15), on every arm of the policy, with the manifest's reason; no declared entry is needed."""

    def test_rule(self):
        import gate_policy
        import compare as cmp
        sc = lambda sid, tags, truth=None: {"id": sid, "area": "T", "route": "POST /traverse",  # noqa: E731
                                            "request": {"seeds": [{"sequence": "A" * 40}]}, "truth": truth or {},
                                            "expect": [{"check": "status", "args": {"http": 200}}], "repeat": 1,
                                            "timing_class": "slow", "tags": tags,
                                            "mcp": {"tool": "traverse_fetch", "arguments": {"preset": "walk"}}}
        man = {"scenarios": [sc("T-a", ["walk"]), sc("T-b", ["walk", "budget_stress"], {"budget_stress": {
            "reason": "it depends on the load"}}), sc("T-c", ["budget_stress"])], "relations": [
            {"id": "REL-x", "type": "same_answer", "a": "T-a", "b": "T-b", "title": "", "args": {}}]}
        pol = gate_policy.build_policy(man, arms=("mcp", "direct-local"), table=({}, None))
        ex = dict((e["scenario"], e) for e in pol["exclusions"])
        self.assertEqual(sorted(ex), ["T-b", "T-c"])
        self.assertEqual((ex["T-b"]["rules"], ex["T-b"]["arms"], ex["T-b"]["reason"]),
                         (["budget_stress"], ["direct-local", "mcp"], "it depends on the load"))
        self.assertIn("budget_stress", ex["T-c"]["reason"])
        self.assertEqual(sorted(set(r["scenario"] for r in pol["required"])), ["T-a"])
        self.assertEqual(pol["relations"], [])
        self.assertEqual(pol["counts"]["exclusions_by_rule"], {"budget_stress": 2})
        self.assertIn("budget_stress: T-b excluded on direct-local, mcp: it depends on the load",
                      gate_policy.describe(pol))


class McpClassesInTheManifest(unittest.TestCase):
    """The owner's decision (2026-10-15): what the MCP cannot show is classed in the manifest, not declared at the
    gate: C-caps-0001 mcp.not_implemented; label_dict_contains of the T-labels pairs mcp.not_implemented_checks;
    T-edge-0008's answer_fields mcp.by_design_checks, its only property check, so the builder marks the whole scenario
    (VERIFY-G (c)); gate_policy derives the MCP exemptions; the builder validates."""

    def sc(self, **mcp):
        s = ml.scenario("T-x-0001", "T", "t", "POST /traverse", ml.traverse_request([{"sequence": "A" * 40}]), {},
                        [{"check": "status", "args": {"http": 200}}, {"check": "answer_fields", "args": {}}])
        s["mcp"].update(mcp)
        return s

    def test_validation(self):
        import build_manifest as bm
        ok = self.sc(not_implemented_checks=[{"check": "answer_fields", "reason": "r"}],
                     by_design_checks=[{"check": "answer_fields", "case": "traverse.duplicate_seeds"}])
        ok["expect"].append({"check": "walk_in_record", "args": {}})     # a property check the marks leave
        self.assertEqual([p for p in bm.validate([ok]) if "mcp." in p or "property-level" in p], [])
        for bad in (self.sc(not_implemented_checks=[{"check": "label_dict_contains", "reason": "r"}]),  # not expected
                    self.sc(not_implemented_checks=[{"check": "answer_fields"}]),                       # no reason
                    self.sc(by_design_checks=[{"check": "answer_fields", "case": "no.such_case"}])):    # no case
            self.assertTrue([p for p in bm.validate([bad]) if "mcp." in p], bad["mcp"])
        ni = self.sc(by_design_checks=[{"check": "answer_fields", "case": "traverse.duplicate_seeds"}])
        ni["mcp"] = dict(not_implemented="x", by_design_checks=ni["mcp"]["by_design_checks"])
        self.assertTrue([p for p in bm.validate([ni]) if "property-level" in p])
        # marks that cover every property check are a scenario-level mark: the builder refuses the list
        every = self.sc(by_design_checks=[{"check": "answer_fields", "case": "traverse.duplicate_seeds"}])
        self.assertTrue([p for p in bm.validate([every]) if "cover every property check" in p])

    def test_by_design_ids_in_the_table(self):
        # VERIFY-V5A re-verification item 1: every mcp.by_design id is a case of the bundled table (F-gate's
        # 422644b7: the renamed discovery cases), else the build fails
        import build_manifest as bm
        s = self.sc(by_design=["traverse.duplicate_seeds", "resolve.no_such_case"])
        self.assertEqual([p for p in bm.validate([s]) if "by-design table" in p],
                         ["T-x-0001: mcp.by_design names resolve.no_such_case, not a case of the bundled by-design "
                          "table"])
        cases = bm._table_cases()
        for cid in ("resolve.discover_limit_above_cap", "resolve.discover_limit_above_cap_host_answers"):
            self.assertIn(cid, cases)

    def test_scenario_level_marks(self):
        # VERIFY-G (c), the owner's rule: marks covering every property check become mcp.not_implemented (the
        # manifest has no whole-scenario by-design form), the generic checks (status, invariants, budget_respected)
        # never count; marks that leave a property check stay property-level
        import build_manifest as bm
        every = self.sc(by_design_checks=[{"check": "answer_fields", "case": "traverse.duplicate_seeds"}],
                        by_design=["traverse.duplicate_seeds"])
        every["expect"] += [{"check": "invariants", "args": {}}, {"check": "budget_respected", "args": {}}]
        both = self.sc(not_implemented_checks=[{"check": "answer_fields", "reason": "no field"}])
        both["id"] = "T-x-0002"
        both["expect"].append({"check": "walk_in_record", "args": {}})
        self.assertEqual(bm.scenario_level_marks([every, both]), ["T-x-0001"])
        self.assertEqual(every["mcp"], {
            "not_implemented": "service-side by design: the table's case traverse.duplicate_seeds; the MCP answer "
                               "cannot carry answer_fields",
            "by_design": ["traverse.duplicate_seeds"]})
        self.assertEqual(both["mcp"]["not_implemented_checks"], [{"check": "answer_fields", "reason": "no field"}])
        self.assertEqual([p for p in bm.validate([every, both]) if "mcp" in p or "property" in p], [])
        nim = self.sc(not_implemented_checks=[{"check": "answer_fields", "reason": "no field"}])
        bm.scenario_level_marks([nim])
        self.assertEqual(nim["mcp"]["not_implemented"], "not implemented yet via MCP: no field (answer_fields)")

    def test_policy_derives_no_property_exemption(self):
        # the property-level marks live in the manifest and gate.py applies them (suite_common.checks_marked_aside:
        # only a not_exposed reading is set aside); a policy exemption of the same check would hide its missing
        # result, failure or inconclusive on MCP, so the policy derives none and both scenarios stay required
        import gate_policy
        a = self.sc(not_implemented_checks=[{"check": "answer_fields", "reason": "no field"}])
        a["id"] = "T-a"
        b = self.sc(by_design_checks=[{"check": "answer_fields", "case": "traverse.duplicate_seeds"}])
        b["id"] = "T-b"
        pol = gate_policy.build_policy({"scenarios": [a, b], "relations": []}, arms=("mcp", "direct-local"),
                                       table=({}, None))
        self.assertEqual(pol["properties"]["not_required"], [])
        self.assertEqual(pol["counts"]["properties_not_required"], 0)
        self.assertIn("property checks not required: 0", gate_policy.describe(pol))
        for arm in ("mcp", "direct-local"):
            self.assertEqual(sorted(e["scenario"] for e in pol["required"] if e["arm"] == arm), ["T-a", "T-b"])

    @unittest.skipUnless(have_cache(), "no record cache")
    def test_built(self):
        s = by_id(built())
        self.assertIn("not_implemented", s["C-caps-0001"]["mcp"])
        for n in (3, 4, 5, 6):
            self.assertEqual(s["T-labels-%04d" % n]["mcp"]["not_implemented_checks"],
                             [{"check": "label_dict_contains", "reason": areas.LABEL_DICT_NOT_VIA_MCP}])
        # T-edge-0008's only property check is service-side by design: the whole scenario is marked
        self.assertEqual(s["T-edge-0008"]["mcp"], {
            "not_implemented": "service-side by design: the table's case traverse.duplicate_seeds; the MCP answer "
                               "cannot carry answer_fields",
            "by_design": ["traverse.duplicate_seeds"]})
        import build_manifest as bm
        self.assertFalse([x["id"] for x in s.values() if bm.covers_every_property(x)])
        cases = bm._table_cases()
        self.assertFalse([(x["id"], c) for x in s.values() for c in x["mcp"].get("by_design") or [] if c not in cases])
        with open(os.path.join(SUITE, "gate", "declared_v5_rt.json")) as f:
            d = json.load(f)
        self.assertEqual((d["exclusions"], d["not_required"], d["relations_excluded"]), ([], [], []))


class WarmRepeatGuard(unittest.TestCase):
    """VERIFY-V5A should-fix (a): the builder refuses a manifest whose run-unique scenario lands on the paired plan's
    warm repeat (every tenth scenario), which would resend its attempt id and meet the 409."""

    def test_guard(self):
        import build_manifest as bm
        import suite_common as sc_
        scs = [{"id": "X-%d" % i} for i in range(25)]
        scs[13]["run_unique"] = ["attempt_id"]
        self.assertEqual(bm.warm_repeat_problems(scs), [])
        scs[20]["run_unique"] = ["attempt_id"]
        self.assertEqual(len(bm.warm_repeat_problems(scs)), 1)
        scs[20]["expect"] = [{"check": "status", "args": {"http": 409}}]     # the duplicate meets its 409 again
        self.assertEqual(bm.warm_repeat_problems(scs), [])
        self.assertEqual(sc_.repeat_count(20, {"repeat": 1}, bm.WARM_EVERY), 2)   # the rule the guard mirrors

    def test_stop_designs_and_placement(self):
        # the re-verification: a stated-stop design decides on one execution too; keep_off_warm_repeat moves the
        # nearest movable scenario of the area from before it to right after it (warm-ups and their targets stay
        # together, run-unique scenarios stay in place)
        import build_manifest as bm

        def mk(i, area="R", **kw):
            return dict({"id": i, "area": area, "tags": []}, **kw)
        scs = [mk("R-%d" % i) for i in range(9)] + [mk("R-tw", warms="R-b", tags=["warm-up"]),
                                                     mk("R-b", truth={"stop_design": {"rule": "x"}})] + \
            [mk("R-%d" % i) for i in range(9, 12)]
        self.assertEqual(bm.warm_repeat_problems(scs), [
            "R-b: a stated-stop design at manifest index 10, which the warm repeat (every 10th scenario) would repeat"])
        self.assertEqual(bm.keep_off_warm_repeat(scs), [("R-8", "R-b")])
        self.assertEqual([x["id"] for x in scs[6:11]], ["R-6", "R-7", "R-tw", "R-b", "R-8"])
        self.assertEqual(bm.warm_repeat_problems(scs), [])
        # a run-unique pair stays together: the moved scenario goes after the pair
        pair = [mk("T-%d" % i, "T") for i in range(10)] + [mk("T-a", "T", run_unique=["a"]),
                                                           mk("T-d", "T", run_unique=["a"],
                                                              expect=[{"check": "status", "args": {"http": 409}}]),
                                                           mk("T-z", "T")]
        self.assertEqual(bm.keep_off_warm_repeat(pair), [("T-9", "T-d")])
        self.assertEqual([x["id"] for x in pair[8:13]], ["T-8", "T-a", "T-d", "T-9", "T-z"])
        # nothing movable before it in its area: the build fails
        lone = [mk("T-w%d" % i, "T", warms="T-z", tags=["warm-up"]) for i in range(10)] + \
            [mk("T-x", "T", truth={"stop_design": {"rule": "x"}})]
        with self.assertRaises(SystemExit):
            bm.keep_off_warm_repeat(lone)

    @unittest.skipUnless(have_cache(), "no record cache")
    def test_built(self):
        import build_manifest as bm
        m = built()
        self.assertEqual(bm.warm_repeat_problems(m["scenarios"]), [])
        order = [x["id"] for x in m["scenarios"]]
        for sid in ("T-attempts-0002", "R-budget-0001", "R-budget-0002"):
            self.assertNotEqual(order.index(sid) % bm.WARM_EVERY, 0, sid)
        self.assertLess(order.index("T-attempts-0002"), order.index("T-attempts-0003"))


class WarmUpsTravel(unittest.TestCase):
    """VERIFY-V5A should-fix (b): a warm-up names its target (`warms`), the builder keeps the pair adjacent, and every
    selection (ids, areas, limit; a resume's ids) that takes the target takes its warm-up too."""

    def man(self):
        mk = lambda i, area, **kw: dict({"id": i, "area": area, "route": "POST /resolve", "repeat": 1,  # noqa: E731
                                         "timing_class": "fast", "tags": [], "mcp": {"tool": "x"}}, **kw)
        return {"scenarios": [mk("R-a", "R"), mk("R-warm-0001", "R", warms="T-b", tags=["warm-up"]),
                              mk("T-b", "T"), mk("T-c", "T")]}

    def test_selections(self):
        import suite_common as sc_
        m = self.man()
        ids = lambda sel: [s["id"] for _, s in sel]  # noqa: E731
        self.assertEqual(ids(sc_.select_scenarios(m, ids=["T-b"])), ["R-warm-0001", "T-b"])
        self.assertEqual(ids(sc_.select_scenarios(m, areas=["T"])), ["R-warm-0001", "T-b", "T-c"])
        self.assertEqual(ids(sc_.select_scenarios(m, areas=["T"], limit=1)), ["R-warm-0001", "T-b"])
        self.assertEqual(ids(sc_.select_scenarios(m, ids=["T-c"])), ["T-c"])
        self.assertEqual([i for i, _ in sc_.select_scenarios(m, ids=["T-b"])], [1, 2])   # global indices kept
        import gate
        self.assertEqual(sorted(set(s for s, _, _ in gate.planned_executions(m, ids=["T-b"], arms=["direct-local"]))),
                         ["R-warm-0001", "T-b"])

    def test_policy_counts_each_scenario_once(self):
        import gate_policy
        m = dict(self.man(), relations=[])
        pol = gate_policy.build_policy(m, routes=("POST /resolve",), arms=("direct-local",), table=({}, None))
        reps = dict((e["scenario"], e["repeats"]) for e in pol["required"])
        self.assertEqual(reps["T-b"], [0])          # not its warm-up's repeats besides its own
        self.assertEqual(reps["R-warm-0001"], [0])

    @unittest.skipUnless(have_cache(), "no record cache")
    def test_built(self):
        import build_manifest as bm
        m = built()
        order = [x["id"] for x in m["scenarios"]]
        pairs = [(x["id"], x["warms"]) for x in m["scenarios"] if x.get("warms")]
        self.assertEqual(len(pairs), 10)
        for w, target in pairs:
            n = sum(1 for _, t in pairs if t == target)                    # a split pair sits together
            self.assertLessEqual(order.index(target) - order.index(w), n, w)
            self.assertGreater(order.index(target), order.index(w), w)
        smoke = [x["id"] for x in bm.smoke_of(m)["scenarios"]]
        for w, target in pairs:
            self.assertEqual(w in smoke, target in smoke, w)


if __name__ == "__main__":
    unittest.main()
