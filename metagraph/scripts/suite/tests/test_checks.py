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
        self.assertEqual([r["outcome"] for r in res], ["pass", "pass"])
        sc["expect"][0]["args"]["http"] = 200
        res = checks.evaluate(sc, checks.normalise_direct(503, b'{"error": "d", "code": "deadline"}'), "direct-local")
        self.assertEqual([r["outcome"] for r in res], ["fail", "fail"])

    def test_deterministic(self):
        a = ok_answer({"x": 1, "timing": {"elapsed_ms": 3}})
        d = checks.answer_digest(a)
        sc = {"route": "POST /pattern", "request": {}, "truth": {}}
        self.assertEqual(outcome(sc, a, "deterministic", previous=None)[0], "pass")
        self.assertEqual(outcome(sc, a, "deterministic", previous=[d, d])[0], "pass")
        self.assertEqual(outcome(sc, a, "deterministic", previous=[d, "0" * 64])[0], "fail")

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
        self.assertEqual(outcome(sc, ok_answer(bad), "source_absent", strand="+")[0], "fail")   # inconclusive
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
