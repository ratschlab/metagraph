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


if __name__ == "__main__":
    unittest.main()
