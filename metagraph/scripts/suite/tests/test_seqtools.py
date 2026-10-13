"""Offline tests of seqtools (the values the host's SPEC states, recomputed)."""
import os
import random
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import seqtools as st  # noqa: E402


class SeqTools(unittest.TestCase):
    def test_revcomp_and_iupac(self):
        self.assertEqual(st.revcomp("AACGTN"), "NACGTT")
        self.assertEqual(st.revcomp("RYKM"), "KMRY")
        self.assertTrue(st.iupac_match("ARN", "AGT"))
        self.assertFalse(st.iupac_match("ARN", "ACT"))
        self.assertEqual(st.hamming_iupac("ARNT", "ACGA"), [1, 3])

    def test_information_bits(self):
        # SPEC §2: 2 per exact base, 1 per two-base code, 0 per N
        self.assertAlmostEqual(st.information_bits_iupac("ACGT"), 8.0)
        self.assertAlmostEqual(st.information_bits_iupac("RN"), 1.0)
        # SPEC §21.5: an exact 31-mer has |B| = 94, 4,279, 125,644 within 1, 2, 3
        p = "A" * 31
        self.assertEqual(st.mismatch_ball_size_iupac(p, 1), 94)
        self.assertEqual(st.mismatch_ball_size_iupac(p, 2), 4279)
        self.assertEqual(st.mismatch_ball_size_iupac(p, 3), 125644)
        self.assertAlmostEqual(62 - st.mismatch_information_bits_iupac(p, 1), 6.55, places=2)
        # the shortest exact DNA patterns served at 1, 2, 3 substitutions are 15, 18 and 20 (§21.5)
        for m, shortest in ((1, 15), (2, 18), (3, 20)):
            self.assertGreaterEqual(st.mismatch_information_bits_iupac("A" * shortest, m), 24)
            self.assertLess(st.mismatch_information_bits_iupac("A" * (shortest - 1), m), 24)
        # an N never substitutes
        self.assertEqual(st.mismatch_ball_size_iupac("N", 1), 4)

    def test_genetic_codes(self):
        self.assertEqual(st.translate("ATGTGATAA", 1), "M**")
        self.assertEqual(st.translate("ATGTGATAA", 4), "MW*")
        self.assertEqual(st.translate("CTG", 12), "S")
        self.assertEqual(sorted(st.codons_of("*", 1)), ["TAA", "TAG", "TGA"])
        self.assertEqual(st.codons_of("*", 27), [])
        self.assertFalse(st.has_stop(31))
        self.assertEqual(len(st.codons_of("X", 1)), 61)
        self.assertEqual(sorted(st.codons_of("B", 1)), sorted(st.codons_of("D", 1) + st.codons_of("N", 1)))
        subs = st.sense_substitutes("M", 1)
        self.assertEqual(len(subs), 60)
        self.assertNotIn("ATG", subs)
        self.assertNotIn("TAA", subs)
        # peptide bits: M and W have one codon each in table 1
        self.assertAlmostEqual(st.information_bits_peptide("MW", 1), 12.0)

    def test_generators(self):
        rng = random.Random(1)
        v, edits = st.substitute_bases("ACGTACGT", [1, 5], rng)
        self.assertEqual([e["pos"] for e in edits], [1, 5])
        for e in edits:
            self.assertNotEqual(e["from"], e["to"])
            self.assertEqual(v[e["pos"]], e["to"])
        s = st.shuffled(rng, "AACCGGTT")
        self.assertEqual(sorted(s), sorted("AACCGGTT"))
        self.assertNotEqual(s, "AACCGGTT")
        pos = st.spread_positions(rng, 31, 5)
        self.assertEqual(len(set(pos)), 5)
        c = st.clustered_positions(rng, 31, 4)
        self.assertEqual(c, list(range(c[0], c[0] + 4)))
        self.assertEqual(st.occurrences_both_strands("AACGTTTT", "ACG"), [(1, "+"), (2, "-")])
        self.assertEqual(st.iupac_code_holding("A", 4, rng), "N")


class Sdust(unittest.TestCase):
    """The sdust port against the host: every stored /pattern answer of an exact DNA pattern without a budget stop
    carries low_complexity_pattern exactly when the port flags the pattern."""

    def test_port_agrees_with_the_host(self):
        import json
        data = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "api",
                                             "python", "tests", "data", "traverse", "pattern"))
        n = flagged = 0
        for name in sorted(os.listdir(data)):
            try:
                with open(os.path.join(data, name, "request.json")) as f:
                    req = json.load(f)
                with open(os.path.join(data, name, "answer.json")) as f:
                    ans = json.load(f)
            except (OSError, ValueError):
                continue
            if not isinstance(req, dict) or not isinstance(ans, dict) or not isinstance(req.get("patterns"), list) \
                    or len(req["patterns"]) != len(ans.get("patterns") or []):
                continue
            for p, e in zip(req["patterns"], ans["patterns"]):
                if not isinstance(p, dict) or "error" in e or "counts" not in e:
                    continue
                kind = [k for k in ("dna", "iupac") if k in p]
                q = p[kind[0]].upper() if kind else ""
                if not q or any(c not in "ACGT" for c in q):
                    continue
                stop = e.get("stop") or {}
                if stop.get("reason") in ("max_steps", "time", "max_annotation_work", "max_memory") \
                        or e.get("determinism") == "time_limited":
                    continue
                want = st.sdust_flags(q)
                self.assertEqual("low_complexity_pattern" in (e.get("notes") or []), want, "%s %s" % (name, q[:40]))
                n += 1
                flagged += want
        self.assertGreater(n, 100)
        self.assertGreater(flagged, 3)

    def test_homopolymers_and_repeats(self):
        # cross-checked against the compiled sdust.c on 3,000 strings (lane S1b): these as the C code says
        self.assertTrue(st.sdust_flags("A" * 31))
        self.assertFalse(st.sdust_flags("TTAGGG" * 5))
        self.assertTrue(st.sdust_flags("TTAGGG" * 15))
        self.assertFalse(st.sdust_flags("ACGTTGCAAGGCTTACCGATCGATTGCAAGT"))


if __name__ == "__main__":
    unittest.main()
