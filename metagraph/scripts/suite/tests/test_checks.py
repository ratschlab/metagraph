"""Offline tests of checks.py: every check on stored fixture answers of the mini index (where a fixture matches the
check's shape) and on hand-made answers, a passing and a failing case each."""
import copy
import gzip
import json
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.dirname(HERE)
sys.path.insert(0, SUITE)

import checks  # noqa: E402

DATA = os.path.normpath(os.path.join(SUITE, "..", "..", "api", "python", "tests", "data", "traverse"))
PATTERN = os.path.join(DATA, "pattern")


def fixture(name):
    with open(os.path.join(PATTERN, name, "request.json")) as f:
        req = json.load(f)
    with open(os.path.join(PATTERN, name, "answer.json")) as f:
        ans = json.load(f)
    return req, ans


def ok_answer(body):
    return checks.normalise_direct(200, body)


def outcome(sc, ans, check, _arm="direct-local", previous=None, **args):
    sc = dict(sc)
    sc["expect"] = [{"check": "status", "args": {"http": ans["http"]}}, {"check": check, "args": args}]
    res = checks.evaluate(sc, ans, _arm, previous=previous)
    return res[1]["outcome"], res[1]["detail"]


# the NDM-F primer context of the fixtures: the k-mer, its offset and one of its placed records
NDM_KMER = "CCAACGGTTTGGCGATCTGGTTTTCCGCCAG"
NDM_F = "GGTTTGGCGATCTGGTTTTC"


def ndm_truth(record="NZ_LPPQ01000025.1", taxid="1296536", start=3681, seq_id=0, length=11137):
    return {"record": record, "taxid": taxid, "seq_id": seq_id, "record_length": length, "start": start,
            "length": 20, "strand": "+", "window": NDM_F, "flank": NDM_KMER, "flank_start": start - 5}


class Normalise(unittest.TestCase):
    def test_direct(self):
        a = checks.normalise_direct(400, b'{"error": "x", "code": "invalid_request"}')
        self.assertEqual(a["error"], {"code": "invalid_request", "text": "x"})
        a = checks.normalise_direct(502, b"<html>bad gateway</html>")
        self.assertEqual(a["body"], "<html>bad gateway</html>")     # suite_common keeps a non-JSON body as text
        self.assertEqual(a["error"], {"code": None, "text": "<html>bad gateway</html>"})
        import gzip as gz
        a = checks.normalise_direct(200, gz.compress(b'{"a": 1}'), {"Content-Encoding": "gzip"})
        self.assertEqual(a["body"], {"a": 1})
        self.assertEqual(a["headers"]["content-encoding"], "gzip")

    def test_digest_ignores_timing(self):
        a = ok_answer({"x": 1, "timing": {"elapsed_ms": 3}, "patterns": [{"timing": {"elapsed_ms": 1}, "v": 2}]})
        b = ok_answer({"x": 1, "timing": {"elapsed_ms": 9}, "patterns": [{"timing": {"elapsed_ms": 7}, "v": 2}]})
        c = ok_answer({"x": 2, "timing": {"elapsed_ms": 9}})
        self.assertEqual(checks.answer_digest(a), checks.answer_digest(b))
        self.assertNotEqual(checks.answer_digest(a), checks.answer_digest(c))

    def test_digest_is_the_runners(self):
        """The fallback digest (a vendored checks.py without suite_common) equals suite_common.answer_digest."""
        import suite_common
        caps = ok_answer({"a": 1, "measured_build_mbps": 3.0, "timing": {"x": 1},
                          "deadline_check": {"observed_max_uninterruptible_ms": 7}})
        caps["headers"] = {"content-encoding": "gzip"}
        saved = checks._suite_common
        try:
            for route in ("GET /traverse/capabilities", "POST /pattern", None):
                want = suite_common.answer_digest(caps, route)
                checks._suite_common = None
                got = checks.answer_digest(caps, route)
                checks._suite_common = saved
                self.assertEqual(got, want, route)
            checks._suite_common = None
            self.assertEqual(checks.normalise_direct(400, b'{"error": "x", "code": "c"}', {"Content-Type": "a/b", "X": "y"}),
                             suite_common.normalise_direct(400, b'{"error": "x", "code": "c"}', {"Content-Type": "a/b", "X": "y"}))
        finally:
            checks._suite_common = saved


class Generic(unittest.TestCase):
    def test_status(self):
        req, ans = fixture("unknown_field")
        a = checks.normalise_direct(400, ans)
        sc = {"route": "POST /pattern", "request": req, "truth": {}}
        self.assertEqual(outcome(sc, a, "status", http=400, code="invalid_request", text="request: unknown field 'bogus'")[0], "pass")
        self.assertEqual(outcome(sc, a, "status", http=400, code="invalid_request", text="request: unknown field 'x'")[0], "fail")
        self.assertEqual(outcome(sc, a, "status", http=400, code="later_increment")[0], "fail")
        self.assertEqual(outcome(sc, a, "status", http=400, text_prefix="request: unknown")[0], "pass")
        self.assertEqual(outcome(sc, a, "status", http=200)[0], "fail")
        self.assertEqual(outcome(sc, a, "status", http=[200, 400], text_contains=["bogus"])[0], "pass")

    def test_unknown_check_and_not_exposed(self):
        sc = {"route": "POST /pattern", "request": {"patterns": [{"dna": "A"}]}, "truth": {},
              "expect": [{"check": "status", "args": {"http": 200}}, {"check": "no_such_check", "args": {}},
                         {"check": "entry_fields", "args": {"pattern": 0, "equals": {"counts.contexts.value": 1}}}]}
        ans = ok_answer({"patterns": [{"id": None}]})
        res = checks.evaluate(sc, ans, "mcp")
        self.assertEqual(res[1]["outcome"], "fail")
        self.assertEqual(res[2]["outcome"], "not_exposed")
        res = checks.evaluate(sc, ans, "direct-local")
        self.assertEqual(res[2]["outcome"], "fail")

    def test_error_answer_skips_other_checks(self):
        sc = {"route": "POST /pattern", "request": {"patterns": []}, "truth": {},
              "expect": [{"check": "status", "args": {"http": [200, 503]}}, {"check": "invariants", "args": {}}]}
        res = checks.evaluate(sc, checks.normalise_direct(503, b'{"error": "d", "code": "deadline"}'), "direct-local")
        # an allowed 503 tests nothing: the body's checks are inconclusive, and so is the execution
        self.assertEqual([r["outcome"] for r in res], ["pass", "inconclusive"])
        self.assertEqual(checks.scenario_outcome(res), "inconclusive")
        sc["expect"][0]["args"]["http"] = 200
        res = checks.evaluate(sc, checks.normalise_direct(503, b'{"error": "d", "code": "deadline"}'), "direct-local")
        self.assertEqual([r["outcome"] for r in res], ["fail", "fail"])

    def test_deterministic(self):
        a = ok_answer({"x": 1, "timing": {"elapsed_ms": 3}})
        d = checks.answer_digest(a)
        sc = {"route": "POST /pattern", "request": {}, "truth": {}}
        # the first repeat has nothing to compare: evaluate leaves the check out
        sc1 = dict(sc, expect=[{"check": "status", "args": {"http": 200}}, {"check": "deterministic", "args": {}}])
        self.assertEqual([r["check"] for r in checks.evaluate(sc1, a, "direct-local", previous=None)], ["status"])
        self.assertEqual(outcome(sc, a, "deterministic", previous=[d, d])[0], "pass")
        self.assertEqual(outcome(sc, a, "deterministic", previous=[d, "0" * 64])[0], "fail")
        # a time-limited answer (this one or an earlier repeat) need not repeat (SPEC §6.8 / §7.9)
        tl = ok_answer({"x": 2, "patterns": [{"stop": {"phase": "discovery", "reason": "time"}}]})
        self.assertEqual(outcome(sc, tl, "deterministic", previous=[d])[0], "inconclusive")
        self.assertEqual(outcome(sc, a, "deterministic", previous=[{"digest": "0" * 64, "time_limited": "seed 0"}])[0],
                         "inconclusive")
        self.assertEqual(outcome(sc, a, "deterministic", previous=[{"digest": d, "time_limited": ""}])[0], "pass")
        walk = ok_answer({"results": [{"arms": {"right": {"cap_trigger": {"reason": "time_budget"}}}}]})
        self.assertEqual(outcome(sc, walk, "deterministic", previous=["0" * 64])[0], "inconclusive")

    def test_fields_encoding_capabilities_clamped(self):
        req, ans = fixture("clamped")
        a = ok_answer(ans)
        sc = {"route": "POST /pattern", "request": req, "truth": {}}
        clamp = ans["limits"]["clamped"][0]
        self.assertEqual(outcome(sc, a, "clamped", field=clamp["field"], requested=clamp["requested"], effective=clamp["effective"])[0], "pass")
        self.assertEqual(outcome(sc, a, "clamped", field=clamp["field"], requested=clamp["requested"], effective=1)[0], "fail")
        self.assertEqual(outcome(sc, a, "clamped", field="nope", requested=1, effective=1)[0], "fail")
        self.assertEqual(outcome(sc, a, "answer_fields", equals={"pattern_contract_version": 1})[0], "pass")
        self.assertEqual(outcome(sc, a, "answer_fields", equals={"pattern_contract_version": 2})[0], "fail")
        self.assertEqual(outcome(sc, a, "answer_fields", absent=["limits"])[0], "fail")
        h = checks.normalise_direct(200, ans, {"Content-Encoding": "gzip"})
        self.assertEqual(outcome(sc, h, "content_encoding", value="gzip")[0], "pass")
        self.assertEqual(outcome(sc, h, "content_encoding", value="deflate")[0], "fail")
        self.assertEqual(outcome(sc, a, "content_encoding", value="gzip")[0], "fail")      # no headers: direct fail
        caps = ok_answer(fixture("pattern_capabilities")[1])
        self.assertEqual(outcome(sc, caps, "capabilities_doc", required=["caps"], equals={"k": 31}, max_bytes=32768)[0], "pass")
        self.assertEqual(outcome(sc, caps, "capabilities_doc", equals={"k": 21})[0], "fail")
        self.assertEqual(outcome(sc, caps, "capabilities_doc", max_bytes=10)[0], "fail")


class PatternChecks(unittest.TestCase):
    def test_invariants_on_every_stored_pattern_answer(self):
        bad = []
        n = 0
        for name in sorted(os.listdir(PATTERN)):
            p = os.path.join(PATTERN, name, "answer.json")
            q = os.path.join(PATTERN, name, "request.json")
            if not (os.path.exists(p) and os.path.exists(q)):
                continue
            with open(q) as f:
                req = json.load(f)
            with open(p) as f:
                ans = json.load(f)
            if not isinstance(req, dict) or not isinstance(ans, dict) or "patterns" not in ans \
                    or not isinstance(req.get("patterns"), list) or len(req["patterns"]) != len(ans["patterns"]):
                continue
            if "graphs" in req:
                continue
            sc = {"route": "POST /pattern", "request": req, "truth": {}}
            for check in ("invariants", "counts_consistent"):
                o, d = outcome(sc, ok_answer(ans), check)
                n += 1
                if o != "pass":
                    bad.append("%s %s: %s" % (name, check, d))
        self.assertGreater(n, 150)
        self.assertEqual(bad, [])

    def test_invariants_refuse_a_mutated_answer(self):
        req, ans = fixture("mismatch_neighbours")
        sc = {"route": "POST /pattern", "request": req, "truth": {}}
        self.assertEqual(outcome(sc, ok_answer(ans), "invariants")[0], "pass")
        bad = copy.deepcopy(ans)
        bad["patterns"][0]["results"][2]["substitutions"][0]["found"] = "T"
        self.assertEqual(outcome(sc, ok_answer(bad), "invariants")[0], "fail")
        bad = copy.deepcopy(ans)
        bad["patterns"][0]["results"].reverse()
        self.assertEqual(outcome(sc, ok_answer(bad), "invariants")[0], "fail")
        bad = copy.deepcopy(ans)
        bad["patterns"][0]["counts"]["contexts"]["by_mismatches"]["1"]["value"] = 5
        self.assertEqual(outcome(sc, ok_answer(bad), "invariants")[0], "fail")
        req2, ans2 = fixture("all_or_count")
        sc2 = {"route": "POST /pattern", "request": req2, "truth": {}}
        bad = copy.deepcopy(ans2)
        e = bad["patterns"][0]
        e["results"][0]["instance"] = "A" * len(e["results"][0]["instance"])
        self.assertEqual(outcome(sc2, ok_answer(bad), "invariants")[0], "fail")
        bad = copy.deepcopy(ans2)
        bad["patterns"][0]["returned"] += 1
        self.assertEqual(outcome(sc2, ok_answer(bad), "invariants")[0], "fail")

    def test_counts_consistent_refuses(self):
        req, ans = fixture("all_or_count")
        sc = {"route": "POST /pattern", "request": req, "truth": {}}
        bad = copy.deepcopy(ans)
        e = bad["patterns"][0]
        e["results"] = e["results"][:-1]
        e["returned"] = len(e["results"])
        self.assertEqual(outcome(sc, ok_answer(bad), "counts_consistent")[0], "fail")
        req, ans = fixture("max_steps")
        sc = {"route": "POST /pattern", "request": req, "truth": {}}
        self.assertEqual(outcome(sc, ok_answer(ans), "counts_consistent")[0], "pass")
        if len(ans["patterns"]) > 1:
            bad = copy.deepcopy(ans)
            bad["patterns"][1]["stop"] = None
            self.assertEqual(outcome(sc, ok_answer(bad), "counts_consistent")[0], "fail")

    def test_source_context_label_placement(self):
        req, ans = fixture("labels_all")
        sc = {"route": "POST /pattern", "request": req, "truth": ndm_truth()}
        a = ok_answer(ans)
        self.assertEqual(outcome(sc, a, "contains_source_context", strand="+", require="all", level=0)[0], "pass")
        self.assertEqual(outcome(sc, a, "contains_source_context", strand="-", require="all")[0], "fail")
        self.assertEqual(outcome(sc, a, "contains_source_context", strand="+", level=1)[0], "fail")
        self.assertEqual(outcome(sc, a, "source_label", strand="+")[0], "pass")
        self.assertEqual(outcome(sc, a, "source_placement", strand="+")[0], "pass")
        wrong = dict(sc, truth=ndm_truth(start=3682))
        self.assertEqual(outcome(wrong, a, "source_placement", strand="+")[0], "fail")
        wrong = dict(sc, truth=ndm_truth(taxid="999999"))
        self.assertEqual(outcome(wrong, a, "source_label", strand="+")[0], "fail")
        bad = copy.deepcopy(ans)
        bad["patterns"][0]["results"] = [r for r in bad["patterns"][0]["results"] if r["kmer"] != NDM_KMER]
        self.assertEqual(outcome(sc, ok_answer(bad), "contains_source_context", strand="+")[0], "fail")
        self.assertEqual(outcome(sc, ok_answer(bad), "source_absent", strand="+")[0], "pass")   # gone from a complete list
        bad["patterns"][0]["retrieval_complete"] = False
        self.assertEqual(outcome(sc, ok_answer(bad), "source_absent", strand="+")[0], "inconclusive")   # not complete
        # a withheld list hides the source: inconclusive, never a fail of a host that withheld rightly
        hidden = copy.deepcopy(ans)
        hidden["patterns"][0].update(results=[], returned=0, retrieval_complete=False,
                                     withheld={"reason": "count_above_threshold"})
        o, d = outcome(sc, ok_answer(hidden), "contains_source_context", strand="+", require="all", level=0)
        self.assertEqual(o, "inconclusive")
        self.assertIn("count_above_threshold", d)
        self.assertEqual(outcome(sc, ok_answer(hidden), "source_label", strand="+")[0], "inconclusive")
        self.assertEqual(outcome(sc, ok_answer(hidden), "source_placement", strand="+")[0], "inconclusive")
        self.assertEqual(outcome(sc, a, "source_absent", strand="+")[0], "fail")                # present
        # the predicate projection places 562's record
        req, ans = fixture("predicate_only_record")
        sc = {"route": "POST /pattern", "request": req,
              "truth": ndm_truth(record="NZ_CP021936.1", taxid="562", start=58761, seq_id=0, length=138340)}
        self.assertEqual(outcome(sc, ok_answer(ans), "source_placement", strand="+")[0], "pass")
        self.assertEqual(outcome(sc, ok_answer(ans), "predicate_selected", strand="+", selected=True,
                                 selection_strand_in=["both"], label="562")[0], "pass")
        self.assertEqual(outcome(sc, ok_answer(ans), "predicate_selected", strand="+", selected=False)[0], "fail")
        self.assertEqual(outcome(sc, ok_answer(ans), "predicate_selected", strand="+", selected=True,
                                 selection_strand_in=["context"], label="562")[0], "fail")

    def test_mismatch_source_and_absent(self):
        req, ans = fixture("mismatch_neighbours")
        P = req["patterns"][0]["dna"]
        variant_rec = "GTTATCGTTTTTCAGCAGCCCTTTGTATCTG"   # the level-2 + context's k-mer
        t = {"record": "r", "taxid": "573", "seq_id": 0, "record_length": 1000, "start": 100, "length": 31,
             "strand": "+", "window": variant_rec, "flank": variant_rec, "flank_start": 100}
        sc = {"route": "POST /pattern", "request": req, "truth": t}
        subs = [{"position": 17, "expected": "A", "found": "G"}, {"position": 21, "expected": "C", "found": "T"}]
        self.assertEqual(checks.expected_substitutions("dna", P, variant_rec, "+"), subs)
        a = ok_answer(ans)
        self.assertEqual(outcome(sc, a, "contains_source_context", strand="+", level=2, substitutions=subs)[0], "pass")
        self.assertEqual(outcome(sc, a, "contains_source_context", strand="+", level=1)[0], "fail")
        self.assertEqual(outcome(sc, a, "contains_source_context", strand="+", level=2, substitutions=subs[:1])[0], "fail")
        self.assertEqual(outcome(sc, a, "source_absent", strand="+")[0], "fail")
        far = dict(t, window="A" * 31, flank="A" * 31)
        self.assertEqual(outcome(dict(sc, truth=far), a, "source_absent", strand="+")[0], "pass")
        self.assertEqual(outcome(sc, a, "count_relation", level=2, relation="exact", value=2)[0], "pass")
        self.assertEqual(outcome(sc, a, "count_relation", level=2, relation="exact", min=3)[0], "fail")
        # the reverse orientation reads the instance in P's orientation
        minus = "CAGATACAAAGGGCTGCTGAAAAACGATAAC"
        self.assertEqual(checks.expected_substitutions("dna", P, minus, "-"), subs)

    def test_slots_withheld_cut_stop_notes(self):
        req, ans = fixture("information_floor")
        sc = {"route": "POST /pattern", "request": req, "truth": {}}
        a = ok_answer(ans)
        txt = "pattern: 22.0 information bits, below the floor of 24.0 for scope any_offset"
        self.assertEqual(outcome(sc, a, "pattern_slot_error", pattern=0, code="information_below_floor", text=txt)[0], "pass")
        self.assertEqual(outcome(sc, a, "pattern_slot_error", pattern=0, code="bad_alphabet")[0], "fail")
        self.assertEqual(outcome(sc, a, "pattern_slot_error", pattern=3, code="information_below_floor")[0], "fail")
        req, ans = fixture("all_or_count_withheld")
        sc = {"route": "POST /pattern", "request": req, "truth": {}}
        a = ok_answer(ans)
        reason = ans["patterns"][0]["withheld"]["reason"]
        self.assertEqual(outcome(sc, a, "withheld", reason=reason)[0], "pass")
        self.assertEqual(outcome(sc, a, "withheld", reason=None)[0], "fail")
        self.assertEqual(outcome(sc, a, "retrieval_complete", value=False)[0], "pass")
        self.assertEqual(outcome(sc, a, "retrieval_complete", value=True)[0], "fail")
        req, ans = fixture("partial_cut")
        sc = {"route": "POST /pattern", "request": req, "truth": {}}
        a = ok_answer(ans)
        self.assertEqual(outcome(sc, a, "cut", reason="max_contexts")[0], "pass")
        self.assertEqual(outcome(sc, a, "cut", reason="max_steps")[0], "fail")
        req, ans = fixture("max_steps")
        sc = {"route": "POST /pattern", "request": req, "truth": {}}
        a = ok_answer(ans)
        self.assertEqual(outcome(sc, a, "stop", phase="discovery", reason="max_steps")[0], "pass")
        self.assertEqual(outcome(sc, a, "stop", null=True)[0], "fail")
        req, ans = fixture("count_low_complexity")
        sc = {"route": "POST /pattern", "request": req, "truth": {}}
        a = ok_answer(ans)
        has = [i for i, e in enumerate(ans["patterns"]) if "low_complexity_pattern" in e.get("notes", [])]
        self.assertTrue(has)
        self.assertEqual(outcome(sc, a, "notes", pattern=has[0], has=["low_complexity_pattern"])[0], "pass")
        self.assertEqual(outcome(sc, a, "notes", pattern=has[0], lacks=["low_complexity_pattern"])[0], "fail")
        self.assertEqual(outcome(sc, a, "entry_fields", pattern=has[0], equals={"scope": ans["patterns"][has[0]]["scope"]})[0], "pass")

    def test_absent(self):
        req, ans = fixture("paths")
        sc = {"route": "POST /pattern", "request": req, "truth": {}}
        a = ok_answer(ans)
        absent = [i for i, e in enumerate(ans["patterns"]) if e.get("counts", {}).get("paths", {}).get("value") == 0]
        present = [i for i, e in enumerate(ans["patterns"]) if (e.get("counts", {}).get("paths", {}).get("value") or 0) > 0]
        self.assertEqual(outcome(sc, a, "absent", pattern=absent[0], count="paths")[0], "pass")
        self.assertEqual(outcome(sc, a, "absent", pattern=present[0], count="paths")[0], "fail")

    def test_paths(self):
        req, ans = fixture("paths_labels")
        sc = {"route": "POST /pattern", "request": req, "truth": {}}
        e = ans["patterns"][0]
        r = e["results"][0]
        lab = next(l for l in r["labels"] if l["support"] == "record_verified")
        occ = lab["occurrence_list"][0]
        s, eend = [int(x) for x in occ["nt_coords"].split("-")]
        t = {"record": occ["record"], "taxid": lab["column"], "seq_id": occ["seq_id"], "record_length": occ["nt_length"],
             "start": s - 1, "length": eend - s + 1, "strand": occ["strand"], "window": r["sequence"],
             "flank": r["sequence"], "flank_start": s - 1}
        sc["truth"] = t
        a = ok_answer(ans)
        self.assertEqual(outcome(sc, a, "path_spells", strand=r["strand"])[0], "pass")
        self.assertEqual(outcome(sc, a, "path_record_verified", strand=occ["strand"])[0], "pass")
        wrong = dict(t, start=t["start"] + 1)
        self.assertEqual(outcome(dict(sc, truth=wrong), a, "path_record_verified", strand=occ["strand"])[0], "fail")
        other = dict(t, window="A" * t["length"])
        self.assertEqual(outcome(dict(sc, truth=other), a, "path_spells", strand=r["strand"])[0], "fail")
        self.assertEqual(outcome(sc, a, "path_absent", sequence=r["sequence"])[0], "fail")
        if e.get("retrieval_complete"):
            self.assertEqual(outcome(sc, a, "path_absent", sequence="C" * t["length"])[0], "pass")

    def test_motif_and_mirror(self):
        req, ans = fixture("motif_context")
        sc = {"route": "POST /pattern", "request": req, "truth": {}}
        a = ok_answer(ans)
        self.assertEqual(outcome(sc, a, "motif", selected=True, decided_by_in=["every_context"], label_present="562")[0], "pass")
        self.assertEqual(outcome(sc, a, "motif", selected=False)[0], "fail")
        self.assertEqual(outcome(sc, a, "motif", label_present="546")[0], "fail")
        req, ans = fixture("predicate_selection_strands")
        minus = next(r for r in ans["patterns"][0]["results"] if r["strand"] == "-")
        rec_kmer = checks.revcomp(minus["kmer"])
        t = {"record": "r", "taxid": "546", "seq_id": 0, "record_length": 10 ** 6, "start": 1000 + minus["offset"],
             "length": 20, "strand": "+", "window": rec_kmer[minus["offset"]:minus["offset"] + 20],
             "flank": rec_kmer, "flank_start": 1000}
        sc = {"route": "POST /pattern", "request": req, "truth": t}
        self.assertEqual(outcome(sc, ok_answer(ans), "predicate_either_mirror", label="546")[0], "pass")
        bad = copy.deepcopy(ans)
        for r in bad["patterns"][0]["results"]:
            if r["strand"] == "-":
                r["selection_strands"] = ["context"] * len(r["selection_strands"])
        self.assertEqual(outcome(sc, ok_answer(bad), "predicate_either_mirror", label="546")[0], "fail")


class MoreChecks(unittest.TestCase):
    """The checks the reviews' verification asked to make falsifiable: each fails on a plausible wrong answer."""

    def entry(self, **e):
        return ok_answer({"patterns": [e]})

    def sc(self, req=None):
        return {"route": "POST /pattern", "request": req or {"patterns": [{"dna": "A" * 31}]}, "truth": {}}

    def test_low_complexity_note(self):
        sc = self.sc()
        noted = self.entry(notes=["low_complexity_pattern"], stop=None, counts={})
        plain = self.entry(notes=[], stop=None, counts={})
        self.assertEqual(outcome(sc, noted, "low_complexity_note", expected=True)[0], "pass")
        self.assertEqual(outcome(sc, plain, "low_complexity_note", expected=True)[0], "fail")
        self.assertEqual(outcome(sc, noted, "low_complexity_note", expected=False)[0], "fail")
        self.assertEqual(outcome(sc, plain, "low_complexity_note", expected=False)[0], "pass")
        budget = self.entry(notes=[], stop={"phase": "discovery", "reason": "max_steps"}, counts={})
        self.assertEqual(outcome(sc, budget, "low_complexity_note", expected=True)[0], "inconclusive")
        both = self.entry(notes=["low_complexity_pattern"], stop={"phase": "discovery", "reason": "time"}, counts={})
        self.assertEqual(outcome(sc, both, "low_complexity_note", expected=True)[0], "fail")   # never beside it
        thr = self.entry(notes=["low_complexity_pattern"], stop={"phase": "discovery", "reason": "max_contexts"}, counts={})
        self.assertEqual(outcome(sc, thr, "low_complexity_note", expected=True)[0], "pass")    # a threshold stop

    def test_stops_stated_and_levels(self):
        sc = self.sc({"patterns": [{"dna": "A" * 20}, {"dna": "C" * 20}]})
        ok = ok_answer({"patterns": [{"stop": None, "counts": {"contexts": {"value": 3, "relation": "exact"}}},
                                     {"stop": {"phase": "discovery", "reason": "time"},
                                      "counts": {"contexts": {"value": None, "relation": "unknown"}}}]})
        self.assertEqual(outcome(sc, ok, "stops_stated", reasons=["time"])[0], "pass")
        bad = copy.deepcopy(ok)
        bad["body"]["patterns"][1]["counts"]["contexts"] = {"value": 5, "relation": "exact"}
        self.assertEqual(outcome(sc, bad, "stops_stated", reasons=["time"])[0], "fail")      # exact after a discovery stop
        bad = copy.deepcopy(ok)
        bad["body"]["patterns"][1]["stop"]["reason"] = "max_steps"
        self.assertEqual(outcome(sc, bad, "stops_stated", reasons=["time"])[0], "fail")
        bad = copy.deepcopy(ok)
        bad["body"]["patterns"][0]["counts"]["contexts"] = {"value": 3, "relation": "at_least"}
        self.assertEqual(outcome(sc, bad, "stops_stated", reasons=["time"])[0], "fail")      # inexact without a stop
        lv = self.entry(results=[{"mismatches": 0}] * 5)
        self.assertEqual(outcome(sc, lv, "results_levels", max_level=0, min_returned=5)[0], "pass")
        self.assertEqual(outcome(sc, self.entry(results=[{"mismatches": 0}, {"mismatches": 1}]), "results_levels",
                                 max_level=0)[0], "fail")
        self.assertEqual(outcome(sc, self.entry(results=[]), "results_levels", max_level=0)[0], "fail")

    def test_count_relation_after_a_budget_stop(self):
        sc = self.sc()
        stopped = self.entry(stop={"phase": "discovery", "reason": "time"},
                             counts={"contexts": {"value": 900, "relation": "at_least"}})
        self.assertEqual(outcome(sc, stopped, "count_relation", relation="exact", min=1000)[0], "inconclusive")
        unstopped = self.entry(stop=None, counts={"contexts": {"value": 900, "relation": "at_least"}})
        self.assertEqual(outcome(sc, unstopped, "count_relation", relation="exact", min=1000)[0], "fail")
        self.assertEqual(outcome(sc, self.entry(stop=None, counts={"contexts": {"value": 900, "relation": "exact"}}),
                                 "count_relation", relation="exact", min=1000)[0], "fail")
        self.assertEqual(outcome(sc, stopped, "count_relation", relation=["exact", "at_least"], min=800)[0], "pass")

    def test_mirror_needs_a_mirror(self):
        req, ans = fixture("predicate_selection_strands")
        minus = next(r for r in ans["patterns"][0]["results"] if r["strand"] == "-")
        rec_kmer = checks.revcomp(minus["kmer"])
        t = {"record": "r", "taxid": "546", "seq_id": 0, "record_length": 10 ** 6, "start": 1000 + minus["offset"],
             "length": 20, "strand": "+", "window": rec_kmer[minus["offset"]:minus["offset"] + 20],
             "flank": rec_kmer, "flank_start": 1000}
        sc = {"route": "POST /pattern", "request": req, "truth": t}
        self.assertEqual(outcome(sc, ok_answer(ans), "predicate_either_mirror", label="546")[0], "pass")
        none = copy.deepcopy(ans)
        none["patterns"][0]["results"] = [r for r in none["patterns"][0]["results"] if r["strand"] != "-"]
        self.assertEqual(outcome(sc, ok_answer(none), "predicate_either_mirror", label="546")[0], "fail")

    def test_argument_schema(self):
        self.assertEqual(checks.validate_expectation({"check": "count_relation", "args": {"relation": "exact", "min": 1}}), [])
        self.assertTrue(checks.validate_expectation({"check": "count_relation", "args": {"relation": ["exact", "at_least",
                                                                                                  "bounds"], "min": 0}}))
        self.assertTrue(checks.validate_expectation({"check": "count_relation", "args": {}}))
        self.assertTrue(checks.validate_expectation({"check": "walk_in_record", "args": {"label_scop": "record"}}))
        self.assertTrue(checks.validate_expectation({"check": "walk_end", "args": {}}))
        self.assertTrue(checks.validate_expectation({"check": "nope", "args": {}}))


class BudgetChecks(unittest.TestCase):
    """budget_respected by each route's SPEC rule and the answer's own deadline record, no fixed allowance; on
    /traverse the seed phase is outside the seed's budget (the smoke run's lambda walk: 34,201 ms under 15,000, a
    seed phase of 33,803 ms with a 20,556 ms read, stopped 19,201 ms after the deadline: within the rule)."""

    def traverse(self, elapsed, seed_phase, sad=None, piece=100.0, budget=15000, echoed=True, timing=True,
                 complete_to_bp=1, stopped=True):
        r = {"outcome": {"walks": "partial"},
             "arms": {"right": {"status": "truncated", "complete_to_bp": complete_to_bp,
                                "cap_trigger": {"reason": "time_budget", "at_bp": complete_to_bp}}}}
        if timing:
            dl = {"longest_piece": {"kind": "read", "ms": piece, "rows": 1}}
            if stopped:
                dl.update(stopped_by="time_budget", stop_after_deadline_ms=sad if sad is not None else
                          max(0.0, elapsed - budget))
            r["timing"] = {"elapsed_ms": elapsed, "seed_phase_ms": seed_phase, "deadline": dl}
        body = {"results": [r], "timing": {"elapsed_ms": elapsed + 1}}
        if echoed:
            body["strategy"] = {"bounds": {"time_budget_ms": budget}}
        return ok_answer(body)

    def sc(self, route, req):
        return {"route": route, "request": req, "truth": {}}

    def test_traverse(self):
        sc = self.sc("POST /traverse", {"seeds": [{"sequence": "A" * 40}],
                                        "strategy": {"direction": "right", "bounds": {"time_budget_ms": 15000}}})
        o, d = outcome(sc, self.traverse(34201.3, 33803.4, sad=19201.0, piece=20556.4), "budget_respected")
        self.assertEqual(o, "pass")
        self.assertIn("seed phase (33803 ms) outlasted the budget", d)
        self.assertIn("within its longest piece", d)
        # the walk saw its deadline 2,500 ms late with a longest piece of 100 ms
        o, d = outcome(sc, self.traverse(17600, 2000, sad=2500), "budget_respected")
        self.assertEqual(o, "finding")
        self.assertIn("2500 ms after max(budget 15000, seed phase 2000)", d)
        self.assertEqual(outcome(sc, self.traverse(15200, 2000, sad=90), "budget_respected")[0], "pass")
        # a budget that passed in the seed phase leaves one level
        o, d = outcome(sc, self.traverse(20500, 20000, sad=5400, piece=6000, complete_to_bp=7), "budget_respected")
        self.assertEqual(o, "finding")
        self.assertIn("one level expected", d)
        # the effective budget is the echoed strategy's (a clamp lowers it), else the request's
        self.assertEqual(outcome(sc, self.traverse(16500, 100, budget=60000, stopped=False), "budget_respected")[0],
                         "pass")
        self.assertEqual(outcome(sc, self.traverse(16600, 100, sad=1500, echoed=False), "budget_respected")[0],
                         "finding")
        # a walked seed without timing or seed phase: a direct arm fails, the MCP arm does not expose it
        self.assertEqual(outcome(sc, self.traverse(1, 0, timing=False), "budget_respected")[0], "fail")
        self.assertEqual(outcome(sc, self.traverse(1, 0, timing=False), "budget_respected", "mcp")[0], "not_exposed")

    def test_pattern_and_resolve(self):
        sc = self.sc("POST /pattern", {"patterns": [{"dna": "A" * 25}], "time_budget_ms": 300})
        over = ok_answer({"patterns": [{"stop": {"phase": "discovery", "reason": "time"}}],
                          "limits": {"time_budget_ms": 300}, "timing": {"elapsed_ms": 412.5}})
        o, d = outcome(sc, over, "budget_respected")
        self.assertEqual(o, "finding")
        self.assertIn("p0 stop discovery/time", d)
        within = ok_answer({"patterns": [{}], "limits": {"time_budget_ms": 300}, "timing": {"elapsed_ms": 299.0}})
        self.assertEqual(outcome(sc, within, "budget_respected")[0], "pass")
        sc = self.sc("POST /resolve", {"sequence": "A" * 40, "labels": ["562"], "bounds": {"time_budget_ms": 70000}})
        clamped = ok_answer({"num_kmers": 10, "limits": {"time_budget_ms": 60000}, "timing": {"elapsed_ms": 61000}})
        self.assertEqual(outcome(sc, clamped, "budget_respected")[0], "finding")
        clamped["body"]["timing"]["elapsed_ms"] = 59000
        self.assertEqual(outcome(sc, clamped, "budget_respected")[0], "pass")

    def test_every_stored_pattern_answer_within_its_budget(self):
        n = 0
        for name in sorted(os.listdir(PATTERN)):
            p, q = os.path.join(PATTERN, name, "answer.json"), os.path.join(PATTERN, name, "request.json")
            if not (os.path.exists(p) and os.path.exists(q)):
                continue
            with open(p) as f:
                ans = json.load(f)
            with open(q) as f:
                req = json.load(f)
            if not isinstance(ans, dict) or "limits" not in ans or "timing" not in ans or not isinstance(req, dict):
                continue
            o, d = outcome({"route": "POST /pattern", "request": req, "truth": {}}, ok_answer(ans),
                           "budget_respected")
            self.assertEqual(o, "pass", "%s: %s" % (name, d))
            n += 1
        self.assertGreater(n, 100)

    def test_outcome_of_an_execution(self):
        def o(*pairs):
            return checks.scenario_outcome([{"check": c, "outcome": x} for c, x in pairs])
        # precedence fail > finding > inconclusive > pass
        self.assertEqual(o(("status", "pass"), ("budget_respected", "finding")), "finding")
        self.assertEqual(o(("status", "pass"), ("budget_respected", "finding"), ("absent", "fail")), "fail")
        self.assertEqual(o(("status", "pass"), ("budget_respected", "finding"), ("absent", "inconclusive")), "finding")
        self.assertEqual(o(("status", "pass")), "pass")                                  # a refusal scenario
        self.assertEqual(o(("status", "pass"), ("invariants", "pass")), "pass")
        self.assertEqual(o(("status", "pass"), ("invariants", "pass"), ("counts_consistent", "pass"),
                           ("budget_respected", "pass"), ("contains_source_context", "inconclusive")), "inconclusive")
        self.assertEqual(o(("status", "pass"), ("contains_source_context", "pass"), ("source_label", "not_exposed")),
                         "pass")
        self.assertEqual(o(("invariants", "not_exposed"), ("source_label", "not_exposed")), "inconclusive")
        self.assertEqual(o(), "inconclusive")                                             # no check at all


class Relations(unittest.TestCase):
    """The reference semantics of the manifest's relations (README, "Relations"): each passes on a consistent pair
    and fails on a plausible wrong one."""

    def rel(self, rtype, **args):
        return {"id": "REL-x", "type": rtype, "a": "A", "b": "B", "args": args}

    def ev(self, rel, a, b, sa=None, sb=None, arm="direct-local"):
        sa = sa or {"route": "POST /pattern", "request": {}, "truth": {}}
        sb = sb or {"route": "POST /pattern", "request": {}, "truth": {}}
        r = checks.evaluate_relation(rel, sa, a, sb, b, arm)
        return r["outcome"], r["detail"]

    def resolve(self, num, graph, runs, stop=None):
        b = {"num_kmers": num, "graph_runs": graph, "labels": [{"label": "562", "runs": runs,
                                                                  "kmers_supported": sum(y - x for x, y in runs)}],
             "timing": {"elapsed_ms": 5}}
        if stop:
            b["stop"] = {"resolved_kmers": num}
        return ok_answer(b)

    def test_same_answer(self):
        a = ok_answer({"patterns": [{"id": "w", "counts": {"contexts": {"value": 7}}}], "timing": {"elapsed_ms": 1}})
        b = ok_answer({"patterns": [{"id": "w", "counts": {"contexts": {"value": 7}}}], "timing": {"elapsed_ms": 9}})
        self.assertEqual(self.ev(self.rel("same_answer"), a, b)[0], "pass")
        c = copy.deepcopy(b)
        c["body"]["patterns"][0]["counts"]["contexts"]["value"] = 8
        o, d = self.ev(self.rel("same_answer"), a, c)
        self.assertEqual(o, "fail")
        self.assertIn("patterns.0.counts.contexts.value", d)
        self.assertEqual(self.ev(self.rel("same_answer", ignore=["patterns"]), a, c)[0], "pass")
        caps = ok_answer({"release": "r1", "attempts": {"server_instance": "s"}, "pattern": {"k": 31}})
        tcaps = ok_answer({"release": "r1", "attempts": {"server_instance": "s"}, "pattern": {"k": 31, "details": "x"}})
        rel = self.rel("same_answer", paths=[["release", "release"], ["attempts.server_instance", "attempts.server_instance"],
                                             ["pattern.k", "pattern.k"]])
        self.assertEqual(self.ev(rel, caps, tcaps)[0], "pass")
        tcaps["body"]["release"] = "r2"
        self.assertEqual(self.ev(rel, caps, tcaps)[0], "fail")
        self.assertEqual(self.ev(self.rel("same_answer", paths=[["", "pattern"]]), ok_answer({"k": 31}), caps)[0], "pass")
        self.assertEqual(self.ev(self.rel("same_answer", paths=[["", "pattern"]]), ok_answer({"k": 21}), caps)[0], "fail")
        self.assertEqual(self.ev(rel, caps, ok_answer({"release": "r1"}))[0], "fail")              # a path missing
        self.assertEqual(self.ev(rel, caps, ok_answer({"release": "r1"}), arm="mcp")[0], "not_exposed")

    def test_prefix_equals(self):
        full = self.resolve(100, [[0, 100]], [[0, 40], [60, 100]])
        stopped = self.resolve(50, [[0, 50]], [[0, 40]], stop=True)
        rel = self.rel("prefix_equals", labels=["562"], a_from_kmer=0, b_from_kmer=0)
        self.assertEqual(self.ev(rel, stopped, full)[0], "pass")
        wrong = self.resolve(50, [[0, 50]], [[0, 41]], stop=True)       # a k-mer the full profile does not support
        self.assertEqual(self.ev(rel, wrong, full)[0], "fail")
        cont = self.resolve(50, [[0, 50]], [[10, 50]])                 # the resolve from k-mer 50
        rel = self.rel("prefix_equals", labels=["562"], a_from_kmer=50, b_from_kmer=0)
        self.assertEqual(self.ev(rel, cont, full)[0], "pass")
        self.assertEqual(self.ev(rel, self.resolve(50, [[0, 50]], [[0, 50]]), full)[0], "fail")
        self.assertEqual(self.ev(rel, self.resolve(0, [], [], stop=True), full)[0], "inconclusive")   # nothing common
        self.assertEqual(self.ev(rel, checks.normalise_direct(503, b'{"code": "deadline"}'), full)[0], "inconclusive")

    def test_zero_equals_absent(self):
        q = "ACGT" * 7 + "ACG"
        sa = {"route": "POST /pattern", "request": {"patterns": [{"dna": q}], "strands": "forward"}, "truth": {}}
        sb = {"route": "POST /resolve", "request": {"sequence": q}, "truth": {}}
        rel = self.rel("zero_equals_absent", pattern=0, count="contexts")
        zero = ok_answer({"patterns": [{"counts": {"contexts": {"value": 0, "relation": "exact"}}}]})
        one = ok_answer({"patterns": [{"counts": {"contexts": {"value": 1, "relation": "exact"}}}]})
        absent, present = ok_answer({"num_kmers": 1, "graph_runs": []}), ok_answer({"num_kmers": 1, "graph_runs": [[0, 1]]})
        self.assertEqual(self.ev(rel, zero, absent, sa, sb)[0], "pass")
        self.assertEqual(self.ev(rel, one, present, sa, sb)[0], "pass")
        self.assertEqual(self.ev(rel, zero, present, sa, sb)[0], "fail")
        self.assertEqual(self.ev(rel, one, absent, sa, sb)[0], "fail")
        unk = ok_answer({"patterns": [{"counts": {"contexts": {"value": 3, "relation": "at_least"}}}]})
        self.assertEqual(self.ev(rel, unk, absent, sa, sb)[0], "inconclusive")

    def test_subset(self):
        def ent(results, value, complete=True, rel_="exact"):
            return ok_answer({"patterns": [{"results": results, "retrieval_complete": complete,
                                            "counts": {"contexts": {"value": value, "relation": rel_}}}]})
        x = {"kmer": "A" * 31, "offset": 2, "strand": "+"}
        y = {"kmer": "C" * 31, "offset": 2, "strand": "+"}
        z = {"kmer": "G" * 31, "offset": 2, "strand": "-", "mismatches": 1}
        rel = self.rel("subset", a_pattern=0, b_pattern=0, key="context")
        self.assertEqual(self.ev(rel, ent([x], 1), ent([x, y], 2))[0], "pass")
        self.assertEqual(self.ev(rel, ent([x, y], 2), ent([x], 1))[0], "fail")              # b complete, y missing
        self.assertEqual(self.ev(rel, ent([x, y], 2), ent([x], 9, complete=False, rel_="at_least"))[0], "inconclusive")
        self.assertEqual(self.ev(rel, ent([x], 3), ent([x, y], 2))[0], "fail")             # exact counts out of order
        rel = self.rel("subset", a_pattern=0, b_pattern=0, key="context", b_level=0)
        self.assertEqual(self.ev(rel, ent([x, y], 2), ent([x, y, z], 3))[0], "pass")
        self.assertEqual(self.ev(rel, ent([z], 1), ent([x, y, z], 3))[0], "fail")          # z is level 1 in b


class ResolveChecks(unittest.TestCase):
    def body(self, runs, num=70, graph=None, stop=None):
        b = {"k": 31, "num_kmers": num, "graph_runs": graph if graph is not None else [[0, num]],
             "labels": [{"label": "562", "kind": "column", "kmers_supported": sum(r[1] - r[0] for r in runs), "runs": runs},
                        {"label": "9606", "kind": "column", "kmers_supported": 0, "runs": []}]}
        if stop:
            b["stop"] = stop
        return ok_answer(b)

    def test_runs(self):
        sc = {"route": "POST /resolve", "request": {"sequence": "A" * 100, "labels": ["562", "9606"]}, "truth": {}}
        full = self.body([[0, 70]])
        self.assertEqual(outcome(sc, full, "resolve_runs", label="562", full=True)[0], "pass")
        self.assertEqual(outcome(sc, full, "invariants")[0], "pass")
        self.assertEqual(outcome(sc, full, "label_not_full", label="9606")[0], "pass")
        self.assertEqual(outcome(sc, full, "label_not_full", label="562")[0], "fail")
        self.assertEqual(outcome(sc, full, "resolve_graph_runs", runs=[[0, 70]])[0], "pass")
        self.assertEqual(outcome(sc, full, "resolve_graph_runs", runs=[[0, 20]])[0], "fail")
        part = self.body([[0, 30], [40, 70]])
        self.assertEqual(outcome(sc, part, "resolve_runs", label="562", full=True)[0], "fail")
        self.assertEqual(outcome(sc, part, "resolve_runs", label="562", covers=[0, 30])[0], "pass")
        self.assertEqual(outcome(sc, part, "resolve_runs", label="562", covers=[0, 31])[0], "fail")
        self.assertEqual(outcome(sc, part, "resolve_runs", label="562", run_end_in=[25, 35])[0], "pass")
        self.assertEqual(outcome(sc, part, "resolve_runs", label="562", run_start_in=[35, 45])[0], "pass")
        self.assertEqual(outcome(sc, part, "resolve_runs", label="562", kmers=60)[0], "pass")
        self.assertEqual(outcome(sc, part, "resolve_runs", label="562", max_kmers=59)[0], "fail")
        self.assertEqual(outcome(sc, part, "resolve_runs", label="NOPE", full=True)[0], "fail")
        broken = self.body([[0, 30]], graph=[[0, 20]])
        self.assertEqual(outcome(sc, broken, "invariants")[0], "fail")
        wrongnum = self.body([[0, 69]], num=69, graph=[[0, 69]])
        self.assertEqual(outcome(sc, wrongnum, "invariants")[0], "fail")

    def test_stop_prefix_and_absent(self):
        sc = {"route": "POST /resolve", "request": {"sequence": "A" * 100, "labels": ["562", "9606"]}, "truth": {}}
        stop = {"phase": "rows", "reason": "time", "resolved_kmers": 40, "query_kmers": 70, "resolved_bp": 70,
                "query_bp": 100, "remainder_from_bp": 40}
        a = self.body([[0, 40]], num=40, graph=[[0, 40]], stop=stop)
        self.assertEqual(outcome(sc, a, "invariants")[0], "pass")
        self.assertEqual(outcome(sc, a, "resolve_runs", label="562", full=True, unless_stopped=True)[0], "pass")
        self.assertEqual(outcome(sc, a, "resolve_runs", label="562", full=True, prefix_if_stopped=True)[0], "pass")
        # a stop that claims the record's own k-mers unsupported, or absent from the graph, fails
        hole = self.body([[0, 30]], num=40, graph=[[0, 40]], stop=stop)
        self.assertEqual(outcome(sc, hole, "resolve_runs", label="562", full=True, prefix_if_stopped=True)[0], "fail")
        gone = self.body([], num=40, graph=[], stop=stop)
        self.assertEqual(outcome(sc, gone, "resolve_runs", label="562", full=True, prefix_if_stopped=True)[0], "fail")
        zero = self.body([], num=0, graph=[], stop=dict(stop, resolved_kmers=0, resolved_bp=0, remainder_from_bp=0))
        self.assertEqual(outcome(sc, zero, "resolve_runs", label="562", full=True, prefix_if_stopped=True)[0],
                         "inconclusive")
        bad = copy.deepcopy(a)
        bad["body"]["stop"]["remainder_from_bp"] = 41
        self.assertEqual(outcome(sc, bad, "invariants")[0], "fail")
        empty = ok_answer({"k": 31, "num_kmers": 70, "graph_runs": [], "labels": [
            {"label": "562", "kind": "column", "kmers_supported": 0, "runs": []}]})
        self.assertEqual(outcome(sc, empty, "resolve_absent")[0], "pass")
        self.assertEqual(outcome(sc, self.body([[0, 70]]), "resolve_absent")[0], "fail")

    def test_discover(self):
        sc = {"route": "POST /resolve", "request": {"sequence": "A" * 100, "discover": {}}, "truth": {}}
        a = self.body([[0, 70]])
        self.assertEqual(outcome(sc, a, "resolve_discover_contains", label="562")[0], "pass")
        self.assertEqual(outcome(sc, a, "resolve_discover_contains", label="1")[0], "fail")
        a["body"]["labels_truncated"] = {"kept": 2, "total": 9}
        self.assertEqual(outcome(sc, a, "resolve_discover_contains", label="1")[0], "pass")


class TraverseChecks(unittest.TestCase):
    MINI = os.path.join(DATA, "real", "mini_refseq")

    def load(self, cell, kind):
        with gzip.open(os.path.join(self.MINI, "%s.%s.json.gz" % (cell, kind))) as f:
            return json.load(f)

    def truth_from(self, result, side="right"):
        """A truth whose record is the walk itself (left flank + seed + right flank of the first paths)."""
        seed = result["seed"]
        req_seed = None
        arms = result["arms"]
        lf = checks._arm_flanks(arms["left"], "left")[0][1] if "left" in arms else ""
        rf = checks._arm_flanks(arms["right"], "right")[0][1] if "right" in arms else ""
        return lf, rf

    def test_walk_in_record_full_and_graphlet(self):
        req = self.load("mini_rep_00__trace_coords", "request")
        full = self.load("mini_rep_00__trace_coords", "full")
        r = full["results"][0]
        seed_seq = req["seeds"][0]["sequence"]
        lf, rf = self.truth_from(r)
        record = "G" * 50 + lf + seed_seq + rf + "T" * 50
        t = {"record": "x", "taxid": "1", "seq_id": 0, "record_length": len(record), "start": 50 + len(lf),
             "length": len(seed_seq), "strand": "+", "window": seed_seq, "flank": record, "flank_start": 0}
        sc = {"route": "POST /traverse", "request": req, "truth": t}
        self.assertEqual(outcome(sc, ok_answer(full), "walk_in_record", label_scope="column")[0], "pass")
        self.assertEqual(outcome(sc, ok_answer(full), "invariants")[0], "pass")
        graphlet = self.load("mini_rep_00__trace_coords", "graphlet")
        self.assertEqual(outcome(sc, ok_answer(graphlet), "walk_in_record", label_scope="column")[0], "pass")
        self.assertEqual(outcome(sc, ok_answer(graphlet), "invariants")[0], "pass")
        mutated = record[:t["start"] + len(seed_seq)] + "A" + record[t["start"] + len(seed_seq) + 1:]
        wrong = dict(t, flank=mutated) if record[t["start"] + len(seed_seq)] != "A" else dict(t, flank=record.replace(rf, "C" * len(rf)))
        self.assertEqual(outcome(dict(sc, truth=wrong), ok_answer(full), "walk_in_record", label_scope="column")[0], "fail")
        self.assertEqual(outcome(sc, ok_answer(full), "walk_in_record", label_scope="column", expect_branching=True)[0], "pass")

    def walk(self, flanks, status="complete", trig=None, reason="max_extension_bp", split_at=None, extra=None):
        """A detail-full result: per side one path whose single segment spells |flanks[side]| (natural
        orientation: a left flank ends at the seed)."""
        arms = {}
        for side, f in flanks.items():
            a = {"status": status, "segments": [{"id": 0, "parents": [], "from_bp": 0, "length_bp": len(f), "sequence": f}],
                 "paths": [{"id": 0, "segments": [0], "length_bp": len(f), "path_reason": reason,
                            "end_reasons": {reason: 1}}]}
            if trig:
                a["cap_trigger"] = {"reason": trig, "at_bp": len(f)}
            if split_at is not None:
                a["splits"] = [{"at_bp": split_at, "segment": 0}]
            a["complete_to_bp"] = 100 if status == "complete" else len(f)
            a.update(extra or {})
            arms[side] = a
        return ok_answer({"results": [{"outcome": {"walks": "complete" if status == "complete" else "partial"},
                                       "arms": arms, "label_dict": [{"name": "NC_1.1", "kind": "header"},
                                                                     {"name": "562", "kind": "column"}]}],
                          "strategy": {"bounds": {"max_extension_bp": 100}}})

    def test_walk_in_record_needs_its_minimum(self):
        import random
        rnd = random.Random(7)
        record = "".join(rnd.choice("ACGT") for _ in range(1000))
        start, L = 400, 40
        t = {"record": "x", "taxid": "1", "seq_id": 0, "record_length": 1000, "start": start, "length": L,
             "strand": "+", "window": record[start:start + L], "flank": record, "flank_start": 0}
        req = {"seeds": [{"sequence": t["window"]}], "strategy": {"direction": "both", "bounds": {"max_extension_bp": 100}}}
        sc = {"route": "POST /traverse", "request": req, "truth": t}
        right, left = record[start + L:start + L + 100], record[start - 100:start]
        self.assertEqual(outcome(sc, self.walk({"right": right, "left": left}), "walk_in_record")[0], "pass")
        # a walk that did not walk: 0 bp, or 1 bp ended without a stated reason, fails
        self.assertEqual(outcome(sc, self.walk({"right": "", "left": left}), "walk_in_record")[0], "fail")
        self.assertEqual(outcome(sc, self.walk({"right": right[:1], "left": left[-1:]}, reason="dead_end"),
                                 "walk_in_record")[0], "fail")
        # the same 1 bp stopped by the time budget: inconclusive, never a pass
        o, d = outcome(sc, self.walk({"right": right[:1], "left": left[-1:]}, status="truncated", trig="time_budget",
                                     reason="time_budget"), "walk_in_record")
        self.assertEqual(o, "inconclusive")
        self.assertIn("time_budget", d)
        # a stated branch end short of the radius passes; prefix_only asks for no minimum
        self.assertEqual(outcome(sc, self.walk({"right": right[:30], "left": left[-30:]}, reason="branch"),
                                 "walk_in_record")[0], "pass")
        self.assertEqual(outcome(sc, self.walk({"right": right[:5], "left": left[-5:]}, status="truncated",
                                               trig="max_steps", reason="max_steps"), "walk_in_record",
                                 prefix_only=True)[0], "pass")
        # a wrong base fails whatever the length
        bad = right[:50] + ("A" if right[50] != "A" else "C") + right[51:]
        self.assertEqual(outcome(sc, self.walk({"right": bad, "left": left}), "walk_in_record")[0], "fail")
        # the record's end bounds the minimum: a seed 30 bp before the end needs 30 bp on the right
        t2 = dict(t, start=930, window=record[930:970])
        sc2 = dict(sc, truth=t2, request={"seeds": [{"sequence": t2["window"]}],
                                          "strategy": {"direction": "right", "bounds": {"max_extension_bp": 100}}})
        self.assertEqual(outcome(sc2, self.walk({"right": record[970:]}, reason="dead_end"), "walk_in_record",
                                 label_scope="record")[0], "pass")
        self.assertEqual(outcome(sc2, self.walk({"right": record[970:] + "A"}, reason="dead_end"), "walk_in_record",
                                 label_scope="record")[0], "fail")                     # past the record's end
        # a seed on the reverse strand: its right arm is the reverse complement of the record's left flank
        tm = dict(t, strand="-")
        scm = dict(sc, truth=tm, request={"seeds": [{"sequence": checks.revcomp(t["window"])}],
                                          "strategy": {"direction": "both", "bounds": {"max_extension_bp": 100}}})
        self.assertEqual(outcome(scm, self.walk({"right": checks.revcomp(left), "left": checks.revcomp(right)}),
                                 "walk_in_record")[0], "pass")
        self.assertEqual(outcome(scm, self.walk({"right": right, "left": left}), "walk_in_record")[0], "fail")
        # walk_end: the end reached passes, a resource stop before it is inconclusive, a short walk without one fails
        sce = dict(sc2)
        self.assertEqual(outcome(sce, self.walk({"right": record[970:]}, reason="label_lost"), "walk_end", arm="right",
                                 reasons_in=["label_lost", "dead_end"])[0], "pass")
        self.assertEqual(outcome(sce, self.walk({"right": record[970:971]}, status="truncated", trig="time_budget",
                                                reason="time_budget"), "walk_end", arm="right")[0], "inconclusive")
        self.assertEqual(outcome(sce, self.walk({"right": record[970:980]}, reason="dead_end"), "walk_end",
                                 arm="right")[0], "fail")
        # an exact minimum from the record: a stated branch short of it fails, beyond it passes
        tx = dict(t, walk={"right": {"min_bp": 60, "first": "ambiguous"}, "left": {"min_bp": 100, "first": "max_extension"}})
        scx = dict(sc, truth=tx)
        self.assertEqual(outcome(scx, self.walk({"right": right[:30], "left": left}, reason="branch"),
                                 "walk_in_record")[0], "fail")
        self.assertEqual(outcome(scx, self.walk({"right": right[:60], "left": left}, reason="branch"),
                                 "walk_in_record")[0], "pass")
        self.assertEqual(outcome(scx, self.walk({"right": right[:30], "left": left[-30:]}, status="truncated",
                                                trig="max_steps", reason="max_steps"), "walk_in_record")[0],
                         "inconclusive")
        # expect_branching: a split at or before the depth where the record branches
        tb = dict(t, walk={"right": {"min_bp": 40, "branch_by_bp": 40}, "left": {"min_bp": 100}})
        scb = dict(sc, truth=tb)
        self.assertEqual(outcome(scb, self.walk({"right": right, "left": left}, split_at=40), "walk_in_record",
                                 expect_branching=True)[0], "pass")
        self.assertEqual(outcome(scb, self.walk({"right": right, "left": left}, split_at=41), "walk_in_record",
                                 expect_branching=True)[0], "fail")
        self.assertEqual(outcome(scb, self.walk({"right": right, "left": left}), "walk_in_record",
                                 expect_branching=True)[0], "fail")
        # the SPEC §7.0 invariant: complete exactly when complete_to_bp is the radius; walks complete only so
        self.assertEqual(outcome(sc, self.walk({"right": right, "left": left}), "invariants")[0], "pass")
        bad = self.walk({"right": right, "left": left})
        bad["body"]["results"][0]["arms"]["right"]["complete_to_bp"] = 1
        self.assertEqual(outcome(sc, bad, "invariants")[0], "fail")
        bad = self.walk({"right": right[:5], "left": left[-5:]}, status="truncated", trig="max_steps", reason="max_steps")
        bad["body"]["results"][0]["outcome"]["walks"] = "complete"
        self.assertEqual(outcome(sc, bad, "invariants")[0], "fail")
        # labels_per_node and the label dictionary
        lp = self.walk({"right": right, "left": left},
                       extra={"labels_per_node": {"cap": 1, "max_seen": 1, "nodes_truncated": 0}})
        self.assertEqual(outcome(sc, lp, "labels_per_node", cap=1)[0], "pass")
        self.assertEqual(outcome(sc, lp, "labels_per_node", cap=64)[0], "fail")
        cut = self.walk({"right": right, "left": left},
                        extra={"labels_per_node": {"cap": 1, "max_seen": 3, "nodes_truncated": 4}, "limitations": []})
        self.assertEqual(outcome(sc, cut, "labels_per_node", cap=1)[0], "fail")
        self.assertEqual(outcome(sc, lp, "label_dict_contains", name="NC_1.1", kind="header")[0], "pass")
        self.assertEqual(outcome(sc, lp, "label_dict_contains", name="562", kind="header")[0], "fail")
        self.assertEqual(outcome(sc, lp, "label_dict_contains", name="NC_1.1", kind="header", only_kind=True)[0], "fail")

        # a seed the scenario allows to fail (a budget that does not hold it)
        failed = ok_answer({"results": [{"outcome": {"walks": "failed"}, "error": "x"}]})
        self.assertEqual(outcome(sc, failed, "walk_in_record")[0], "fail")
        self.assertEqual(outcome(sc, failed, "walk_in_record", if_walked=True)[0], "inconclusive")

    def test_traverse_stop(self):
        sc = {"route": "POST /traverse", "truth": {},
              "request": {"seeds": [{"sequence": "A" * 40}],
                          "strategy": {"direction": "right", "bounds": {"max_extension_bp": 1000, "max_steps": 5,
                                                                         "max_memory_mb": 1}}}}
        def res(walks, rs=None, trig=None, lims=(), arms=True):
            r = {"outcome": {"walks": walks}, "limitations": [], "resource_stop": rs}
            if arms:
                r["arms"] = {"right": {"status": "truncated" if trig else "complete", "complete_to_bp": 5,
                                       "cap_trigger": {"reason": trig} if trig else None,
                                       "limitations": [{"kind": "walk_domain", "knob": k} for k in lims]}}
            else:
                r["limitations"] = [{"kind": "walk_domain", "knob": k} for k in lims]
            return ok_answer({"results": [r]})
        steps = res("partial", trig="max_steps", lims=["bounds.max_steps"])
        self.assertEqual(outcome(sc, steps, "traverse_stop", knob="bounds.max_steps", cap_reason="max_steps")[0], "pass")
        # a host that ignored the knob and walked on
        self.assertEqual(outcome(sc, res("complete"), "traverse_stop", knob="bounds.max_steps", cap_reason="max_steps",
                                 walks=["partial"])[0], "fail")
        self.assertEqual(outcome(sc, res("partial", trig="max_steps"), "traverse_stop", knob="bounds.max_steps",
                                 cap_reason="max_steps")[0], "fail")                     # no limitation names it
        mem = res("failed", rs={"resource": "memory", "requested": 1}, lims=["bounds.max_memory_mb"], arms=False)
        self.assertEqual(outcome(sc, mem, "traverse_stop", knob="bounds.max_memory_mb", resource="memory",
                                 cap_reason="resource_limit", walks=["partial", "failed"])[0], "pass")
        wrong = res("failed", rs={"resource": "work", "requested": 1}, lims=["bounds.max_memory_mb"], arms=False)
        self.assertEqual(outcome(sc, wrong, "traverse_stop", knob="bounds.max_memory_mb", resource="memory",
                                 walks=["partial", "failed"])[0], "fail")
        unset = res("failed", rs={"resource": "memory", "requested": 2}, lims=["bounds.max_memory_mb"], arms=False)
        self.assertEqual(outcome(sc, unset, "traverse_stop", knob="bounds.max_memory_mb", resource="memory",
                                 walks=["partial", "failed"])[0], "fail")


    def test_walk_end_and_outcome(self):
        req = self.load("mini_rep_00__trace_coords", "request")
        full = self.load("mini_rep_00__trace_coords", "full")
        r = full["results"][0]
        seed_seq = req["seeds"][0]["sequence"]
        lf, rf = self.truth_from(r)
        n = len(rf)
        record = lf + seed_seq + rf
        t = {"record": "x", "taxid": "1", "seq_id": 0, "record_length": len(record), "start": len(lf),
             "length": len(seed_seq), "strand": "+", "window": seed_seq, "flank": record, "flank_start": 0}
        sc = {"route": "POST /traverse", "request": req, "truth": t}
        self.assertEqual(outcome(sc, ok_answer(full), "walk_end", arm="right", reasons_in=["max_extension_bp"])[0], "pass")
        self.assertEqual(outcome(sc, ok_answer(full), "walk_end", arm="right", reasons_in=["dead_end"])[0], "fail")
        longer = dict(t, record_length=len(record) + 5)
        self.assertEqual(outcome(dict(sc, truth=longer), ok_answer(full), "walk_end", arm="right")[0], "fail")
        self.assertEqual(outcome(sc, ok_answer(full), "traverse_seed_outcome", walks="complete")[0], "pass")
        self.assertEqual(outcome(sc, ok_answer(full), "traverse_seed_outcome", walks="failed")[0], "fail")

    def test_trace_coordinates(self):
        req = self.load("mini_rep_00__trace_coords", "request")
        full = self.load("mini_rep_00__trace_coords", "full")
        t = {"record": "NC_019089.1", "taxid": "5", "seq_id": 6, "record_length": 10 ** 6, "start": 83467, "length": 150,
             "strand": "+", "window": req["seeds"][0]["sequence"], "flank": "", "flank_start": 0}
        sc = {"route": "POST /traverse", "request": req, "truth": t}
        self.assertEqual(outcome(sc, ok_answer(full), "trace_coordinates")[0], "pass")
        self.assertEqual(outcome(dict(sc, truth=dict(t, start=83468)), ok_answer(full), "trace_coordinates")[0], "fail")
        self.assertEqual(outcome(dict(sc, truth=dict(t, record="NOPE")), ok_answer(full), "trace_coordinates")[0], "fail")

    def test_graphlet_reader_matches_full(self):
        n = 0
        for name in sorted(os.listdir(self.MINI)):
            if not name.endswith(".graphlet.json.gz"):
                continue
            cell = name[:-len(".graphlet.json.gz")]
            g = self.load(cell, "graphlet")
            f = self.load(cell, "full")
            for rg, rf in zip(g["results"], f["results"]):
                if "graphlet" not in rg:
                    continue
                arms = checks.graphlet_arms(rg["graphlet"])
                for side, fa in rf["arms"].items():
                    self.assertEqual([x for _, x in checks._arm_flanks(fa, side)],
                                     [x for _, x in checks._arm_flanks(arms[side], side)], cell)
                    n += 1
        self.assertGreater(n, 4)


if __name__ == "__main__":
    unittest.main()
