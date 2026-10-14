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


_REASON_CODE = dict((v, k) for k, v in checks._END_CODES.items())


def _ranges(ids):
    """MGT RANGES of a set of ids (SPEC §7.5.1)."""
    ids = sorted(ids)
    if not ids:
        return "."
    out, lo = [], ids[0]
    for a, b in zip(ids, ids[1:] + [None]):
        if b != a + 1:
            out.append(str(lo) if lo == a else "%d-%d" % (lo, a))
            lo = b
    return ",".join(out)


def mcp_shaped(answer, json_arms=True):
    """A `detail: full` /traverse answer in the shape the MCP's export (and `detail: graphlet`) gives: each result's
    segments and paths as an MGT graphlet (SPEC-labeled-traversal-core §7.5: A, G, T and R records; with the label
    dictionary as L records and the segments' labels as explicit G entry / end sets where the result states them),
    the arms' other fields (status, cap_trigger, limitations, counters, ...) as JSON beside it (none with json_arms
    False)."""
    out = copy.deepcopy(answer)
    for r in out["body"]["results"]:
        if "arms" not in r:
            continue
        lines = ["H mgt 1 31 basic $ACGT c k m 64 1000 0 walk * * *"]
        ld = r.pop("label_dict", None)
        prev = ""
        for x in ld or []:
            lines.append("L %s 0 * 0 %s" % ("h" if x.get("kind") == "header" else "c", x["name"]))
        for side, arm in r["arms"].items():
            trig = arm.get("cap_trigger")
            cap = "%s,%s,0,1,1,1,0" % (_REASON_CODE[trig["reason"]], trig.get("at_bp", 0)) if trig else "*"
            counters = ",".join("%s=%s" % kv for kv in (arm.get("counters") or {}).items()) or "steps=0"
            lines.append("A %s %s %s u 1 1 1 0 0 %s %s 0 * %d 0 %d 0 0 0" % (
                side[0], arm["status"][0], arm.get("complete_to_bp", 0), counters, cap, len(arm["segments"]),
                len(arm["paths"])))
            leaves = dict((p["segments"][-1], p) for p in arm["paths"])
            for seg in arm["segments"]:
                parents = ",".join(str(x) for x in seg["parents"]) or "*"
                bases = seg["sequence"] if side == "right" else seg["sequence"][::-1]
                sets = "* * *" if "labels" not in seg else "%s %s %s" % (
                    _ranges(seg["labels"]), seg.get("labels_total", "*") if seg.get("labels_truncated") else "*",
                    _ranges(seg["labels"]))
                lines.append("G %s %d %d %s * * * %s" % (parents, seg["from_bp"], seg["length_bp"], sets, bases or "."))
                for ev in seg.get("events") or []:
                    if ev.get("type") == "blocked":
                        lines.append("E %d b %s %s 1 0" % (ev["at_bp"], ev["char"], _REASON_CODE[ev["reason"]]))
                p = leaves.get(seg["id"])
                if p is not None:
                    lines.append("T %s" % (_REASON_CODE.get(p.get("path_reason"), "*")))
                    for reason in p.get("end_reasons") or {}:
                        lines.append("R %d 0 0 %d %s" % (seg["id"], p["length_bp"], _REASON_CODE[reason]))
            if json_arms:
                for key in ("segments", "paths", "splits", "branch_events"):
                    arm.pop(key, None)
        r["graphlet"] = "\n".join(lines) + "\n"
        if not json_arms:
            del r["arms"]
    return out


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
        # G21: a /pattern 503 deadline in place of a 200 is an answer SPEC-pattern §6/§7.6 allows on any host (the
        # README's machine-dependent): inconclusive, never a fail ...
        sc["expect"][0]["args"]["http"] = 200
        dl = checks.normalise_direct(503, b'{"error": "d", "code": "deadline"}')
        res = checks.evaluate(sc, dl, "direct-local")
        self.assertEqual([r["outcome"] for r in res], ["inconclusive", "inconclusive"])
        self.assertIn("503 deadline in place of 200", res[0]["detail"])
        self.assertEqual(checks.scenario_outcome(res), "inconclusive")
        # ... and budget_respected judges it by the client's wall time (a 503 long past the budget is a finding)
        sc_b = dict(sc, request={"patterns": [], "time_budget_ms": 1000},
                    expect=sc["expect"] + [{"check": "budget_respected", "args": {}}])
        self.assertEqual(checks.evaluate(sc_b, dl, "direct-local", client={"wall_ms": 1500.0})[2]["outcome"],
                         "inconclusive")
        self.assertEqual(checks.evaluate(sc_b, dl, "direct-local", client={"wall_ms": 9000.0})[2]["outcome"],
                         "finding")
        # A-01: a 503 deadline BEFORE its budget is not the SPEC's (it answers one only when the answer cannot be
        # written by the budget): v3's M-planted-0014 answered so at a wall of 5 ms against 60,000 ms fails
        sc_m = dict(sc, request={"patterns": [], "time_budget_ms": 60000}, expect=sc_b["expect"])
        res = checks.evaluate(sc_m, dl, "direct-local", client={"wall_ms": 5.0})
        self.assertEqual([r["outcome"] for r in res], ["inconclusive", "inconclusive", "fail"])
        self.assertIn("before the 60000 ms budget", res[2]["detail"])
        self.assertEqual(checks.scenario_outcome(res), "fail")
        self.assertEqual(checks.evaluate(sc_m, dl, "direct-local", client={"wall_ms": 59999.5})[2]["outcome"],
                         "inconclusive")
        self.assertEqual(checks.evaluate(sc_m, dl, "direct-local", client={"wall_ms": 60012.9})[2]["outcome"],
                         "inconclusive")
        # G22 (VERIFY-2 item 3): a /resolve 503 deadline in place of its 200 is SPEC-core §4.5's answer too, before
        # the budget it fails (the host answers it once the elapsed time reaches the budget; the finalisation reserve
        # ends the work, not the answer)
        rsc = {"route": "POST /resolve", "request": {"query": "A" * 60, "bounds": {"time_budget_ms": 20000}},
               "truth": {}, "expect": [{"check": "status", "args": {"http": 200}},
                                       {"check": "invariants", "args": {}}, {"check": "budget_respected", "args": {}}]}
        rdl = checks.normalise_direct(503, b'{"error": "resolve: the answer could not be built", "code": "deadline"}')
        res = checks.evaluate(rsc, rdl, "direct-local", client={"wall_ms": 22522.0})
        self.assertEqual([r["outcome"] for r in res], ["inconclusive", "inconclusive", "inconclusive"])
        self.assertIn("SPEC-labeled-traversal-core §4.5", res[0]["detail"])
        self.assertEqual(checks.evaluate(rsc, rdl, "direct-local", client={"wall_ms": 19800.0})[2]["outcome"],
                         "fail")
        self.assertEqual(checks.evaluate(rsc, rdl, "direct-local", client={"wall_ms": 20000.0})[2]["outcome"],
                         "inconclusive")
        res = checks.evaluate(rsc, rdl, "direct-local", client={"wall_ms": 800.0})
        self.assertEqual(res[2]["outcome"], "fail")
        self.assertIn("before the 20000 ms budget", res[2]["detail"])
        # another code, another route, or a scenario expecting a refusal: as stated
        res = checks.evaluate(sc, checks.normalise_direct(503, b'{"error": "d", "code": "loading"}'), "direct-local")
        self.assertEqual([r["outcome"] for r in res], ["fail", "fail"])
        res = checks.evaluate(dict(sc, route="POST /traverse"), dl, "direct-local")
        self.assertEqual([r["outcome"] for r in res], ["fail", "fail"])
        sc["expect"][0]["args"] = {"http": 400, "code": "invalid_request"}
        res = checks.evaluate(sc, dl, "direct-local")
        self.assertEqual([r["outcome"] for r in res], ["fail", "fail"])
        # a scenario whose status names the deadline is judged as it states
        sc["expect"][0]["args"] = {"http": 200, "text_contains": ["deadline"]}
        self.assertEqual(checks.evaluate(sc, dl, "direct-local")[0]["outcome"], "fail")
        self.assertFalse(checks.pattern_deadline_allowed(sc, dl))
        self.assertFalse(checks.deadline_allowed(sc, dl))

    def test_resolve_deadline_is_opt_in(self):
        # A2-01 (review of lane S5b, VERIFY-REVIEW: SPEC-labeled-traversal-core §4.5 l.224, capabilities
        # resolve.time_budget.default null): a /resolve without bounds.time_budget_ms has no deadline, so v4's
        # R-discover-0001 answered 503 deadline (at 5 ms, at 45 s or with no wall) fails; it was inconclusive
        disc = {"route": "POST /resolve", "request": {"sequence": "ACGT" * 25, "discover": True,
                                                       "bounds": {"max_query_bp": 2000}}, "truth": {},
                "expect": [{"check": "status", "args": {"http": 200}}, {"check": "invariants", "args": {}},
                           {"check": "resolve_discover_contains", "args": {"label": "511145"}}]}
        rdl = checks.normalise_direct(503, b'{"error": "resolve: the answer could not be built", "code": "deadline"}')
        self.assertFalse(checks.deadline_allowed(disc, rdl))
        for client in ({"wall_ms": 5.0}, {"wall_ms": 45000.0}, None):
            res = checks.evaluate(disc, rdl, "direct-local", client=client)
            self.assertEqual([r["outcome"] for r in res], ["fail", "fail", "fail"], client)
            self.assertIn("opt-in", res[0]["detail"])
            self.assertEqual(checks.scenario_outcome(res), "fail")
        # a budget in the request: the SPEC's 503 in place of the 200 (judged by budget_respected's wall rule)
        timed = dict(disc, request=dict(disc["request"], bounds={"time_budget_ms": 20000}),
                     expect=disc["expect"] + [{"check": "budget_respected", "args": {}}])
        self.assertTrue(checks.deadline_allowed(timed, rdl))
        res = checks.evaluate(timed, rdl, "direct-local", client={"wall_ms": 22522.0})
        self.assertEqual(checks.scenario_outcome(res), "inconclusive")
        self.assertEqual(checks.evaluate(timed, rdl, "direct-local", client={"wall_ms": 5.0})[3]["outcome"], "fail")
        # /pattern keeps its default budget (SPEC-pattern-search §6: default_time_budget_ms): allowed without one
        pat = {"route": "POST /pattern", "request": {"patterns": []}, "truth": {},
               "expect": [{"check": "status", "args": {"http": 200}}]}
        pdl = checks.normalise_direct(503, b'{"error": "p", "code": "deadline"}')
        self.assertTrue(checks.deadline_allowed(pat, pdl))

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
        # a time limitation that cannot have happened (a /pattern time stop after 300 ms of 60,000) explains nothing
        psc = {"route": "POST /pattern", "request": {"patterns": [{"dna": "A" * 25, "id": "w"}], "time_budget_ms": 60000},
               "truth": {}}
        early = ok_answer({"patterns": [{"id": "w", "stop": {"phase": "discovery", "reason": "time"},
                                         "determinism": "time_limited"}],
                           "limits": {"time_budget_ms": 60000, "finalize_reserve_ms": 250}, "timing": {"elapsed_ms": 300.0}})
        o, d = outcome(psc, early, "deterministic", previous=["0" * 64])
        self.assertEqual(o, "fail")
        self.assertIn("before the work time", d)
        early["body"]["timing"]["elapsed_ms"] = 59900.0
        self.assertEqual(outcome(psc, early, "deterministic", previous=["0" * 64])[0], "inconclusive")

    def test_deterministic_after_a_deadline(self):
        # G27 / A-03 (VERIFY-2 item 1, variant B): D-M-planted-0015, r0 the host's 503 deadline (or the MCP's relayed
        # 502), r1 and r2 complete and equal: a 503 deadline is time-limited, and a repeat is compared only with the
        # earlier repeats no time budget limited
        sc = {"route": "POST /pattern", "request": {}, "truth": {}}
        dl = checks.normalise_direct(503, b'{"error": "pattern: the answer could not be written", "code": "deadline"}')
        relay = {"http": 502, "body": None, "error": {"code": "search_failed",
                                                     "text": "pattern_deadline: pattern: the answer could not be written"}}
        rdl = {"http": 503, "body": None, "error": {"code": "traversal_resolve_deadline", "text": "the backend could not"}}
        for x, why in ((dl, "a 503 deadline"), (relay, "/pattern 503 deadline relayed"), (rdl, "/resolve 503 deadline")):
            self.assertIn(why, checks.time_limited(x))
        a = ok_answer({"x": 1, "timing": {"elapsed_ms": 3}})
        b = ok_answer({"x": 2, "timing": {"elapsed_ms": 3}})
        da, db, ddl = checks.answer_digest(a), checks.answer_digest(b), checks.answer_digest(dl)
        prev_dl = {"digest": ddl, "time_limited": checks.time_limited(dl)}
        # r1 after a deadline r0: nothing unlimited to compare with
        o, d = outcome(sc, a, "deterministic", previous=[prev_dl])
        self.assertEqual(o, "inconclusive")
        self.assertIn("no earlier repeat without a time limitation", d)
        # r2 equal to r1: pass; r2 differing from r1 (a real non-determinism) still fails (variant A masked it)
        self.assertEqual(outcome(sc, a, "deterministic", previous=[prev_dl, {"digest": da, "time_limited": ""}])[0],
                         "pass")
        o, d = outcome(sc, b, "deterministic", previous=[prev_dl, {"digest": da, "time_limited": ""}])
        self.assertEqual(o, "fail", d)
        # the same with the earlier time limitation in a 200 (the masking of before)
        tl = {"digest": "1" * 64, "time_limited": "pattern 0 time-limited"}
        self.assertEqual(outcome(sc, b, "deterministic", previous=[tl, {"digest": da, "time_limited": ""}])[0], "fail")
        # a deadline answer after complete ones is time-limited itself: inconclusive
        dsc = dict(sc, expect=[{"check": "status", "args": {"http": [200, 503]}},
                               {"check": "deterministic", "args": {}}])
        res = checks.evaluate(dsc, dl, "direct-local", previous=[da, da])
        self.assertEqual(res[1]["outcome"], "inconclusive", res[1]["detail"])
        # the MCP runner's plain digests count as unlimited, as before
        self.assertEqual(outcome(sc, a, "deterministic", previous=[da])[0], "pass")

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

    def test_invariants_refuse_a_result_listed_twice(self):
        # a copy of a context or path in place of a distinct one keeps the list's length and order (equal
        # keys sort as equal), so only the identity tells it
        req = {"mode": "partial", "max_contexts": 2, "patterns": [{"dna": "ACGTACGTAC", "id": "w"}]}
        sc = {"route": "POST /pattern", "request": req, "truth": {}}
        k1, k2 = "T" * 21 + "ACGTACGTAC", "G" * 21 + "ACGTACGTAC"
        r1 = {"kmer": k1, "node": 5, "offset": 21, "instance": "ACGTACGTAC", "strand": "+"}
        r2 = {"kmer": k2, "node": 9, "offset": 21, "instance": "ACGTACGTAC", "strand": "+"}

        def ans(results):
            e = {"id": "w", "mode": "partial", "returned": len(results), "results": results, "stop": None,
                 "withheld": None, "cut": {"reason": "max_contexts"}, "retrieval_complete": False,
                 "counts": {"contexts": {"relation": "exact", "value": 3, "by_offset": dict(
                     (str(o), {"relation": "exact", "value": 3 if o == 21 else 0}) for o in range(22))}}}
            return ok_answer({"patterns": [e], "limits": {"max_contexts": 2}})
        self.assertEqual(outcome(sc, ans([r1, r2]), "invariants")[0], "pass")
        self.assertEqual(outcome(sc, ans([r1, r2]), "counts_consistent")[0], "pass")
        o, d = outcome(sc, ans([r1, dict(r1)]), "invariants")
        self.assertEqual(o, "fail")
        self.assertIn("the result of r0 listed again", d)
        # a path is its nodes in an orientation (L-long-0006 with its '-' path dropped and the '+' path twice)
        P = "ACGT" * 10
        lreq = {"mode": "all_or_count", "long_search": "paths", "patterns": [{"dna": P, "id": "w"}]}
        lsc = {"route": "POST /pattern", "request": lreq, "truth": {}}
        path = {"sequence": P, "instance": P, "offset": 0, "strand": "+", "nodes": list(range(10))}

        def lans(results):
            e = {"id": "w", "mode": "all_or_count", "returned": len(results), "results": results, "stop": None,
                 "withheld": None, "retrieval_complete": True,
                 "counts": {"anchors": {"relation": "exact", "value": 2},
                            "paths": {"relation": "exact", "value": len(results)}}}
            return ok_answer({"patterns": [e], "limits": {"max_paths": 10}})
        self.assertEqual(outcome(lsc, lans([path]), "invariants")[0], "pass")
        o, d = outcome(lsc, lans([path, dict(path)]), "invariants")
        self.assertEqual(o, "fail")
        self.assertIn("listed again", d)

    def test_path_record_verified_placement(self):
        # the record's occurrence is placed with the index's seq_id and the record's length, and every
        # occurrence the label verifies lies in its record (L-verified-0001: NC_1.1 34701-34740 '+', seq_id
        # 2, 640,681 bp)
        t = {"record": "NC_1.1", "taxid": "562", "seq_id": 2, "record_length": 640681, "start": 34700, "length": 40,
             "strand": "+", "window": "ACGT" * 10}
        sc = {"route": "POST /pattern", "truth": t, "request": {"patterns": [{"dna": "ACGT" * 10, "id": "w"}],
                                                                "long_search": "paths", "output": {"labels": "all"}}}

        def ans(other=None, **occ):
            o = {"record": "NC_1.1", "seq_id": 2, "strand": "+", "nt_coords": "34701-34740", "nt_length": 640681}
            o.update(occ)
            path = {"sequence": "ACGT" * 10, "strand": "+", "labels": [
                {"column": "562", "support": "record_verified", "occurrence_list": [o] + ([other] if other else [])}]}
            return ok_answer({"patterns": [{"id": "w", "counts": {"paths": {"relation": "exact", "value": 1}},
                                            "results": [path], "returned": 1}]})
        self.assertEqual(outcome(sc, ans(), "path_record_verified")[0], "pass")
        for bad in ({"nt_length": 1}, {"seq_id": 999999999}, {"nt_length": 1, "seq_id": 999999999},
                    {"nt_length": 640681 + 1000}):
            self.assertEqual(outcome(sc, ans(**bad), "path_record_verified")[0], "fail", bad)
        other = {"record": "NC_2.1", "seq_id": 7, "strand": "-", "nt_coords": "901-940", "nt_length": 5000}
        self.assertEqual(outcome(sc, ans(other=other), "path_record_verified")[0], "pass")
        o, d = outcome(sc, ans(other=dict(other, nt_length=900)), "path_record_verified")
        self.assertEqual(o, "fail")
        self.assertIn("outside its record of nt_length 900", d)
        self.assertEqual(outcome(sc, ans(other=dict(other, nt_coords="940-901")), "path_record_verified")[0], "fail")
        # source_placement places a context's record occurrence by the same rule
        self.assertIsNone(checks._occurrence_problem({"record": "NC_1.1", "seq_id": 2, "nt_coords": "34701-34740",
                                                      "nt_length": 640681}, t))
        self.assertIn("outside its record", checks._coords_problem({"nt_coords": "0-39", "nt_length": 100}))

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
        # retrieval_complete false alone (the labels incomplete) still lists every context (SPEC §7.5: only a
        # cut or a withheld list shortens it)
        bad["patterns"][0]["retrieval_complete"] = False
        self.assertEqual(outcome(sc, ok_answer(bad), "source_absent", strand="+")[0], "pass")
        bad["patterns"][0].update(mode="partial", cut={"reason": "max_contexts"})
        cut = ok_answer(bad)
        cut["body"]["limits"]["max_contexts"] = 3                                     # a legitimate cut hides it
        self.assertEqual(outcome(sc, cut, "source_absent", strand="+")[0], "inconclusive")
        # a withheld list hides the source: inconclusive, never a fail of a host that withheld rightly
        hidden = copy.deepcopy(ans)
        hidden["patterns"][0].update(results=[], returned=0, retrieval_complete=False,
                                     withheld={"reason": "count_above_threshold"})
        hidden["limits"]["max_contexts"] = 10                                         # rightly: the count is above it
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

    def test_low_complexity_note_beside_a_later_phase(self):
        # SPEC-pattern-search §7.8: the note is the engine's, decided after its search; a stop of a later phase (the
        # labels' reads, the selection, the output) leaves it as it is, an engine phase's budget stop leaves it out
        sc = self.sc({"mode": "partial", "output": {"labels": "all"}, "patterns": [{"dna": "A" * 31}]})
        for st in ({"phase": "output", "reason": "time"}, {"phase": "output", "reason": "max_memory"},
                   {"phase": "label_discovery", "reason": "max_annotation_work"},
                   {"phase": "placement", "reason": "max_memory"},
                   {"phase": "selection", "reason": "max_predicate_work"}):
            noted = self.entry(notes=["low_complexity_pattern"], stop=st, counts={})
            plain = self.entry(notes=[], stop=st, counts={})
            self.assertEqual(outcome(sc, noted, "low_complexity_note", expected=True)[0], "pass", st)
            o, d = outcome(sc, plain, "low_complexity_note", expected=True)
            self.assertEqual(o, "fail", st)
            self.assertIn("beside a later phase's stop", d)
            self.assertEqual(outcome(sc, noted, "low_complexity_note", expected=False)[0], "fail", st)
            self.assertEqual(outcome(sc, plain, "low_complexity_note", expected=False)[0], "pass", st)
        for st in ({"phase": "discovery", "reason": "time"}, {"phase": "extension", "reason": "max_memory"},
                   {"phase": "extension", "reason": "max_annotation_work"}, {"reason": "max_steps"}):
            noted = self.entry(notes=["low_complexity_pattern"], stop=st, counts={})
            self.assertEqual(outcome(sc, noted, "low_complexity_note", expected=True)[0], "fail", st)
            self.assertEqual(outcome(sc, self.entry(notes=[], stop=st, counts={}), "low_complexity_note",
                                     expected=True)[0], "inconclusive", st)
        # a supported-path search's selection stop can be its own refusal of a held walk, inside the search
        sup = self.sc({"long_search": "supported_paths", "predicate": {"any": ["1"]}, "predicate_strands": "either",
                       "patterns": [{"dna": "A" * 40}]})
        st = {"phase": "selection", "reason": "max_memory"}
        noted = self.entry(notes=["low_complexity_pattern"], stop=st, counts={})
        plain = self.entry(notes=[], stop=st, counts={})
        self.assertEqual(outcome(sup, noted, "low_complexity_note", expected=True)[0], "pass")
        self.assertEqual(outcome(sup, plain, "low_complexity_note", expected=True)[0], "inconclusive")
        self.assertEqual(outcome(sup, noted, "low_complexity_note", expected=False)[0], "fail")
        self.assertEqual(outcome(sup, plain, "low_complexity_note", expected=False)[0], "pass")
        # the invariants: a stored answer with an output stop keeps a note; a supported-path search's own budget
        # stop does not
        req, ans = fixture("labels_all_output_budget")
        fsc = {"route": "POST /pattern", "request": req, "truth": {}}
        ans = copy.deepcopy(ans)
        ans["patterns"][0]["notes"].append("low_complexity_pattern")
        self.assertEqual(outcome(fsc, ok_answer(ans), "invariants")[0], "pass")
        ans["patterns"][0]["stop"]["reason"] = "time"
        ans["patterns"][0]["determinism"] = "time_limited"
        o, d = outcome(fsc, ok_answer(ans), "invariants")
        self.assertNotIn("low_complexity_pattern", d)
        req, ans = fixture("supported_paths_work_budget")
        fsc = {"route": "POST /pattern", "request": req, "truth": {}}
        self.assertEqual(outcome(fsc, ok_answer(ans), "invariants")[0], "pass")
        ans = copy.deepcopy(ans)
        ans["patterns"][0]["notes"].append("low_complexity_pattern")
        o, d = outcome(fsc, ok_answer(ans), "invariants")
        self.assertEqual(o, "fail")
        self.assertIn("low_complexity_pattern beside a budget stop extension/max_annotation_work", d)

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

    def test_count_relation_after_any_stop_that_leaves_it_inexact(self):
        # SPEC-pattern-search §7.4: at_least means a stop interrupted the count; every budget of §7.6, §14.4, §19.8
        # and §20.5 and every stop_at_threshold stop can be that stop (inconclusive below a minimum), no stop cannot
        def C(rel, v, **kw):
            return dict({"relation": rel, "value": v}, **kw)
        sup = self.sc({"long_search": "supported_paths", "mode": "count", "patterns": [{"dna": "ACGT" * 10}]})
        for reason in ("max_annotation_work", "max_memory"):
            e = self.entry(stop={"phase": "extension", "reason": reason},
                           counts={"anchors": C("exact", 2), "supported_paths": C("at_least", 0)})
            o, d = outcome(sup, e, "count_relation", count="supported_paths", min=1)
            self.assertEqual(o, "inconclusive", d)
            self.assertIn("after the stop extension/%s" % reason, d)
        none = self.entry(stop=None, counts={"anchors": C("exact", 2), "supported_paths": C("at_least", 0)})
        o, d = outcome(sup, none, "count_relation", count="supported_paths", min=1)
        self.assertEqual(o, "fail")
        self.assertIn("with no stop that leaves the count inexact", d)
        # a predicate's pruning states supported_paths at_least without a stop (§20.4, §20.9), any branch pruned
        # before L the paths
        pred = self.sc({"long_search": "supported_paths", "mode": "count", "predicate": {"any": ["1"]},
                        "patterns": [{"dna": "ACGT" * 10}]})
        pruned = self.entry(stop=None, counts={"anchors": C("exact", 2), "paths": C("at_least", 0),
                                               "supported_paths": C("at_least", 0, branches_pruned=2,
                                                                    branches_pruned_by_predicate=1)})
        self.assertEqual(outcome(pred, pruned, "count_relation", count="supported_paths", min=1)[0], "inconclusive")
        self.assertEqual(outcome(pred, pruned, "count_relation", count="paths", min=1)[0], "inconclusive")
        paths = self.sc({"long_search": "paths", "mode": "count", "patterns": [{"dna": "ACGT" * 10}]})
        mem = self.entry(stop={"phase": "extension", "reason": "max_memory"},
                         counts={"anchors": C("exact", 2), "paths": C("at_least", 3)})
        o, d = outcome(paths, mem, "count_relation", count="paths", relation="exact")
        self.assertEqual(o, "inconclusive", d)
        self.assertEqual(outcome(paths, self.entry(stop=None, counts={"anchors": C("exact", 2), "paths": C(
            "at_least", 3)}), "count_relation", count="paths", relation="exact")[0], "fail")
        heavy = self.sc({"mode": "count", "patterns": [{"dna": "GCTACACACACC"}]})
        o, d = outcome(heavy, self.entry(stop=None, counts={"contexts": C("at_least", 3)}), "count_relation", min=20,
                       relation=["exact", "at_least"])
        self.assertEqual(o, "fail")
        steps = self.entry(stop={"phase": "discovery", "reason": "max_steps"}, counts={"contexts": C("at_least", 3)})
        self.assertEqual(outcome(heavy, steps, "count_relation", min=20)[0], "inconclusive")
        # a stop_at_threshold stop leaves the count at_least too (§7.4)
        thr = self.sc({"mode": "count", "stop_at_threshold": True, "max_contexts": 10,
                       "patterns": [{"dna": "GCTACACACACC"}]})
        crossed = self.entry(stop={"phase": "discovery", "reason": "max_contexts"},
                             counts={"contexts": C("at_least", 11)})
        self.assertEqual(outcome(thr, crossed, "count_relation", min=20)[0], "inconclusive")
        self.assertEqual(outcome(thr, crossed, "count_relation", min=11)[0], "pass")
        # bounds decide only where they do not straddle the minimum; a value the bounds exclude fails, also after a
        # stop
        scan = self.entry(stop={"phase": "mask_scan", "reason": "max_steps"},
                          counts={"contexts": C("bounds", 5, lower=5, upper=30)})
        self.assertEqual(outcome(heavy, scan, "count_relation", min=20)[0], "inconclusive")
        self.assertEqual(outcome(heavy, scan, "count_relation", min=5)[0], "pass")
        self.assertEqual(outcome(heavy, scan, "count_relation", min=31)[0], "fail")
        self.assertEqual(outcome(heavy, scan, "count_relation", value=40)[0], "fail")
        self.assertEqual(outcome(heavy, scan, "count_relation", value=20)[0], "inconclusive")
        self.assertEqual(outcome(heavy, steps, "count_relation", value=2)[0], "fail")       # at_least 3 > 2
        self.assertEqual(outcome(heavy, steps, "count_relation", value=7)[0], "inconclusive")
        # a stop of the labels' reads or the output leaves the graph count as the search made it (§7.5): no excuse
        # for a graph count, an excuse for the labels' own counts
        labels = self.sc({"mode": "partial", "output": {"labels": "all"}, "patterns": [{"dna": "GCTACACACACC"}]})
        out = self.entry(stop={"phase": "output", "reason": "max_memory"}, counts={"contexts": C("at_least", 3),
                                                                                 "labels": C("at_least", 2)})
        self.assertEqual(outcome(labels, out, "count_relation", min=20)[0], "fail")
        self.assertEqual(outcome(labels, out, "count_relation", relation="exact")[0], "fail")
        self.assertEqual(outcome(labels, out, "count_relation", count="labels", relation="exact")[0], "inconclusive")
        # a count of relation rules of its own: at_least below a minimum is inconclusive without a stop, an
        # unexpected relation still fails
        plain = self.entry(stop=None, counts={"contexts": C("exact", 30), "labels": C("at_least", 2)})
        self.assertEqual(outcome(labels, plain, "count_relation", count="labels", min=5)[0], "inconclusive")
        self.assertEqual(outcome(labels, plain, "count_relation", count="labels", relation="exact")[0], "fail")

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
        # a list withheld at the deadline hides the mirrors (MCP Q-strands-0001 r0 of the full v4 run)
        late = copy.deepcopy(ans)
        late["patterns"][0].update(results=[], returned=0, retrieval_complete=False, withheld={"reason": "deadline"},
                                   stop={"phase": "selection", "reason": "time"}, determinism="time_limited")
        late["timing"] = dict(late.get("timing") or {}, elapsed_ms=59850.0)
        self.assertEqual(outcome(sc, ok_answer(late), "predicate_either_mirror", label="546")[0], "inconclusive")
        self.assertEqual(outcome(sc, ok_answer(late), "predicate_selected", selected=True, label="546")[0],
                         "inconclusive")
        unstated = copy.deepcopy(late)
        unstated["patterns"][0].update(stop=None, determinism="full")
        self.assertEqual(outcome(sc, ok_answer(unstated), "predicate_either_mirror", label="546")[0], "fail")
        # a hiding the request cannot have (anchor labels without labels all; a selection threshold stop
        # that left the selected count exact) fails the mirror and the selection alike
        for withheld, extra in (("anchor_labels_truncated", {"anchors_truncated": [{"cap": 2, "total": 9}]}),
                                ("threshold_crossed", {"stop": {"phase": "selection", "reason": "max_contexts"}})):
            hid = copy.deepcopy(ans)
            hid["patterns"][0].update(results=[], returned=0, retrieval_complete=False,
                                      withheld={"reason": withheld}, **extra)
            hid.setdefault("limits", {})["stop_at_threshold"] = True
            self.assertEqual(outcome(sc, ok_answer(hid), "predicate_either_mirror", label="546")[0], "fail", withheld)
            self.assertEqual(outcome(sc, ok_answer(hid), "predicate_selected", selected=True, label="546")[0],
                             "fail", withheld)
        refused = copy.deepcopy(ans)
        refused["patterns"][0] = {"id": "w", "error": {"code": "predicate_too_large"}}
        self.assertEqual(outcome(sc, ok_answer(refused), "predicate_either_mirror", label="546")[0], "fail")
        nostrand = copy.deepcopy(ans)
        for r in nostrand["patterns"][0]["results"]:
            r.pop("selection_strands", None)
        self.assertEqual(outcome(sc, ok_answer(nostrand), "predicate_either_mirror", label="546"),
                         ("fail", "field not in the answer: selection_strands"))

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

    def walk(self, elapsed, seed_phase, sad=None, piece=100.0, budget=15000, echoed=True, timing=True, at_bp=1,
             stopped=True, sides=("right",), path_bp=None, used=None, radius=1000, limitation=True):
        """A walked seed: per side one path, truncated by time_budget at at_bp with its walk_domain limitation and
        resource_stop (stopped), else complete to the radius."""
        if used is None:
            used = budget + sad if sad is not None and sad > 0 else elapsed
        arms = {}
        for side in sides:
            n = path_bp if path_bp is not None else (at_bp if stopped else radius)
            reason = "time_budget" if stopped else "max_extension_bp"
            a = {"status": "truncated" if stopped else "complete", "complete_to_bp": at_bp if stopped else radius,
                 "segments": [{"id": 0, "parents": [], "from_bp": 0, "length_bp": n, "sequence": "A" * n}],
                 "paths": [{"id": 0, "segments": [0], "length_bp": n, "path_reason": reason, "end_reasons": {reason: 1}}],
                 "limitations": []}
            if stopped:
                a["cap_trigger"] = {"reason": "time_budget", "at_bp": at_bp}
                if limitation:
                    a["limitations"] = [{"kind": "walk_domain", "knob": "bounds.time_budget_ms", "limit": budget,
                                         "observed": used}]
            arms[side] = a
        r = {"outcome": {"walks": "partial" if stopped else "complete"}, "arms": arms}
        if timing:
            dl = {"longest_piece": {"kind": "read", "ms": piece, "rows": 1}}
            if stopped:
                dl.update(stopped_by="time_budget", stop_after_deadline_ms=sad if sad is not None else
                          max(0.0, elapsed - budget))
                r["resource_stop"] = {"resource": "time", "requested": budget, "effective": budget, "used": used}
            r["timing"] = {"elapsed_ms": elapsed, "seed_phase_ms": seed_phase, "deadline": dl}
        body = {"results": [r], "timing": {"elapsed_ms": elapsed + 1}}
        if echoed:
            body["strategy"] = {"bounds": {"time_budget_ms": budget, "max_extension_bp": radius}}
        return ok_answer(body)

    def sc(self, route, req):
        return {"route": route, "request": req, "truth": {}}

    def walk_sc(self, sides=("right",), budget=15000, radius=1000):
        d = "both" if len(sides) == 2 else sides[0]
        return self.sc("POST /traverse", {"seeds": [{"sequence": "A" * 40}], "strategy": {
            "direction": d, "bounds": {"time_budget_ms": budget, "max_extension_bp": radius}}})

    def br(self, ans, sides=("right",), budget=15000, arm="direct-local"):
        return outcome(self.walk_sc(sides, budget), ans, "budget_respected", arm)

    def test_traverse_conformant_stops_are_no_finding(self):
        # VERIFY-fixes.md X1: answers SPEC §6.8 and walker.cpp allow (no reading at depth 0, the depth-0 reads
        # unpaced, a read followed by its head before the next reading)
        o, d = self.br(self.walk(34201.3, 33803.4, sad=19201.0, piece=20556.4))   # the smoke run's 34 s lambda walk
        self.assertEqual(o, "pass", d)
        self.assertIn("seed phase (33803 ms) outlasted the budget", d)
        self.assertIn("at depth 1", d)
        self.assertIn("within its pieces", d)
        # the same walk whose longest piece is the depth-0 read it made after the seed phase (397 ms): 398 ms late
        # at the first depth-1 checkpoint is within the depth-0 level's pieces
        o, d = self.br(self.walk(34201.3, 33803.4, sad=19201.0, piece=397.0))
        self.assertEqual(o, "pass", d)
        # budget 300, seed phase 250, both arms' unpaced depth-0 reads of 200 ms, stopped at depth 1 at 650 ms
        o, d = self.br(self.walk(652.0, 250.0, sad=350.0, piece=200.0, budget=300, sides=("left", "right")),
                       sides=("left", "right"), budget=300)
        self.assertEqual(o, "pass", d)
        # at a checked depth: 400 ms after the deadline with a longest piece of 300 ms (a read, then its head)
        o, d = self.br(self.walk(15400.0, 500.0, sad=400.0, piece=300.0, at_bp=500))
        self.assertEqual(o, "pass", d)
        self.assertIn("checked depth (500 bp)", d)
        # no time stop, complete to the radius 900 ms after the budget: within the pieces before the last reading
        o, d = self.br(self.walk(15900.0, 500.0, piece=300.0, stopped=False))
        self.assertEqual(o, "pass", d)
        self.assertIn("a level at the radius is never stopped", d)
        # the effective budget is the echoed strategy's (a clamp lowers it), else the request's
        self.assertEqual(outcome(self.walk_sc(budget=60000), self.walk(16500, 100, budget=60000, stopped=False),
                                 "budget_respected")[0], "pass")
        # a walked seed without timing or seed phase: a direct arm fails, the MCP arm does not expose it
        self.assertEqual(self.br(self.walk(1, 0, timing=False))[0], "fail")
        self.assertEqual(self.br(self.walk(1, 0, timing=False), arm="mcp")[0], "not_exposed")

    def test_traverse_broken_stops_are_findings(self):
        # VERIFY-fixes.md X1: answers SPEC §6.8 does not allow
        o, d = self.br(self.walk(2000.0, 500.0, sad=-13000.0, piece=300.0, at_bp=500))
        self.assertEqual(o, "finding", d)                                      # a time stop at 2 s of a 15 s budget
        self.assertIn("before its 15000 ms budget", d)
        o, d = self.br(self.walk(2000.0, 1500.0, sad=0.0, piece=300.0, at_bp=500))
        self.assertEqual(o, "finding", d)                                      # the same, clamped at 0 (walker.cpp)
        o, d = self.br(self.walk(60000.0, 500.0, piece=300.0, stopped=False))
        self.assertEqual(o, "finding", d)                                      # no time stop, to the radius at 60 s
        self.assertIn("without a time stop", d)
        o, d = self.br(self.walk(20000.0, 500.0, sad=5000.0, piece=300.0, at_bp=500))
        self.assertEqual(o, "finding", d)                                      # 5 s late at a checked depth
        self.assertIn("checked depth", d)
        o, d = self.br(self.walk(34201.3, 33803.4, sad=19201.0, piece=20556.4, at_bp=50))
        self.assertEqual(o, "finding", d)                                      # seed phase past the budget, 50 bp walked
        self.assertIn("at most 1 expected", d)
        o, d = self.br(self.walk(60000.0, 33803.4, piece=20556.4, stopped=False))
        self.assertEqual(o, "finding", d)                                      # the same without a stop, to the radius
        o, d = self.br(self.walk(40000.0, 1000.0, sad=25000.0, piece=300.0, at_bp=1))
        self.assertEqual(o, "finding", d)                                      # depth 1, far past the depth-0 pieces
        self.assertIn("at depth 1", d)
        self.assertEqual(outcome(self.walk_sc(), self.walk(16600, 100, sad=1500, echoed=False, at_bp=500),
                                 "budget_respected")[0], "finding")

    def test_503_long_past_the_budget(self):
        # VERIFY-fixes.md: a 503 deadline 8.9 s after a 750 ms budget (the smoke run's R-budget-0002) is a finding,
        # by the client's wall time the runner passes
        def sc(budget, route="POST /resolve"):
            req = {"sequence": "A" * 100, "labels": ["208964"], "bounds": {"time_budget_ms": budget}} \
                if route == "POST /resolve" else {"patterns": [{"dna": "A" * 25, "id": "w"}], "time_budget_ms": budget}
            return {"route": route, "request": req, "truth": {},
                    "expect": [{"check": "status", "args": {"http": [200, 503]}}, {"check": "invariants", "args": {}},
                               {"check": "budget_respected", "args": {}}]}
        late = checks.normalise_direct(503, b'{"code": "deadline", "error": "resolve: the answer could not be built '
                                            b'and written within bounds.time_budget_ms"}')
        res = checks.evaluate(sc(750), late, "direct-local", client={"wall_ms": 9692.9})
        self.assertEqual([r["outcome"] for r in res], ["pass", "inconclusive", "finding"])
        self.assertIn("8943 ms after the 750 ms budget", res[2]["detail"])
        self.assertEqual(checks.scenario_outcome(res), "finding")
        # 0.8 s after an 8,000 ms budget: within max(1 s, the budget), nothing to judge
        res = checks.evaluate(sc(8000), late, "direct-local", client={"wall_ms": 8803.5})
        self.assertEqual(res[2]["outcome"], "inconclusive")
        # no client wall time (the MCP runner passes none): inconclusive as before
        self.assertEqual(checks.evaluate(sc(750), late, "mcp")[2]["outcome"], "inconclusive")
        pat = checks.normalise_direct(503, b'{"code": "deadline", "error": "pattern: ..."}')
        self.assertEqual(checks.evaluate(sc(60000, "POST /pattern"), pat, "direct-remote",
                                         client={"wall_ms": 60906.0})[2]["outcome"], "inconclusive")
        self.assertEqual(checks.evaluate(sc(60000, "POST /pattern"), pat, "direct-remote",
                                         client={"wall_ms": 125000.0})[2]["outcome"], "finding")
        busy = checks.normalise_direct(503, b'{"error": "loading"}')
        self.assertEqual(checks.evaluate(sc(750), busy, "direct-local", client={"wall_ms": 99999.0})[2]["outcome"],
                         "inconclusive")

    def test_503_against_the_capped_budget(self):
        # A2-03 (review of lane S5b): R-edge-0021 asks 70,000 ms, the host's cap is 60,000 (resolve.time_budget.
        # max_time_ms; the host lowers the budget and states the clamp), so a 503 deadline at 60,100 ms is the SPEC's
        # answer: judged against min(requested, cap) when the builder states the cap (budget_respected cap_ms)
        def sc(budget, cap, route="POST /resolve"):
            req = {"sequence": "A" * 100, "labels": ["511145"], "bounds": {"time_budget_ms": budget}} \
                if route == "POST /resolve" else {"patterns": [{"dna": "A" * 25, "id": "w"}], "time_budget_ms": budget}
            args = {} if cap is None else {"cap_ms": cap}
            return {"route": route, "request": req, "truth": {},
                    "expect": [{"check": "status", "args": {"http": 200}}, {"check": "invariants", "args": {}},
                               {"check": "budget_respected", "args": args}]}
        dl = checks.normalise_direct(503, b'{"code": "deadline", "error": "resolve: the answer could not be built '
                                          b'and written within bounds.time_budget_ms (60000 ms, the finalisation '
                                          b'reserve of 250 ms included)"}')
        res = checks.evaluate(sc(70000, 60000), dl, "direct-local", client={"wall_ms": 60100.0})
        self.assertEqual(res[2]["outcome"], "inconclusive", res[2]["detail"])
        self.assertIn("60000 ms budget (the requested 70000 ms lowered to the host's cap)", res[2]["detail"])
        # without the cap (a manifest built before A2-03) the requested budget, as before: a false fail
        self.assertEqual(checks.evaluate(sc(70000, None), dl, "direct-local",
                                         client={"wall_ms": 60100.0})[2]["outcome"], "fail")
        # before the capped budget: fail; long past it: a finding; a cap above the request: the request's budget
        res = checks.evaluate(sc(70000, 60000), dl, "direct-local", client={"wall_ms": 59000.0})
        self.assertEqual(res[2]["outcome"], "fail")
        self.assertIn("before the 60000 ms budget (the requested 70000 ms lowered", res[2]["detail"])
        # the host answers /resolve's 503 only once the elapsed time reaches the budget (Deadline::respond_expired):
        # its finalisation reserve ends the work earlier, not the answer, so a 503 within it is premature too
        for wall, want in ((59900.0, "fail"), (59998.0, "fail"), (59999.0, "inconclusive"), (60000.0, "inconclusive")):
            self.assertEqual(checks.evaluate(sc(70000, 60000), dl, "direct-local",
                                             client={"wall_ms": wall})[2]["outcome"], want, wall)
        for wall, want in ((900.0, "fail"), (999.0, "inconclusive"), (1000.0, "inconclusive"), (600.0, "fail")):
            self.assertEqual(checks.evaluate(sc(1000, None), dl, "direct-local",
                                             client={"wall_ms": wall})[2]["outcome"], want, wall)
        self.assertEqual(checks.evaluate(sc(70000, 60000), dl, "direct-local",
                                         client={"wall_ms": 125000.0})[2]["outcome"], "finding")
        self.assertEqual(checks.evaluate(sc(20000, 60000), dl, "direct-local",
                                         client={"wall_ms": 19000.0})[2]["outcome"], "fail")
        pat = checks.normalise_direct(503, b'{"code": "deadline", "error": "pattern: ..."}')
        self.assertEqual(checks.evaluate(sc(700000, 600000, "POST /pattern"), pat, "direct-local",
                                         client={"wall_ms": 600100.0})[2]["outcome"], "inconclusive")
        self.assertEqual(checks.evaluate(sc(700000, 600000, "POST /pattern"), pat, "direct-local",
                                         client={"wall_ms": 599000.0})[2]["outcome"], "fail")
        self.assertEqual(checks.CHECK_ARGS["budget_respected"], ((), ("rule", "cap_ms")))

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
        # VERIFY-fixes.md H3: the property checks all not exposed (an MCP P-exact-0001 without counts), the
        # generic ones passing: inconclusive, never pass
        self.assertEqual(o(("status", "pass"), ("invariants", "pass"), ("count_relation", "not_exposed"),
                           ("entry_fields", "not_exposed"), ("budget_respected", "not_exposed")), "inconclusive")
        self.assertEqual(o(("status", "pass"), ("invariants", "pass"), ("counts_consistent", "pass"),
                           ("budget_respected", "pass")), "pass")                       # no property check: generic
        self.assertEqual(o(("status", "pass"), ("count_relation", "pass"), ("entry_fields", "not_exposed")), "pass")
        self.assertEqual(o(("status", "pass"), ("count_relation", "weird")), "fail")   # outside the vocabulary
        self.assertTrue(checks.is_property("count_relation") and checks.is_property("deterministic"))
        self.assertFalse(any(checks.is_property(n) for n in ("status", "invariants", "counts_consistent",
                                                             "budget_respected")))
        sc = {"route": "POST /pattern", "request": {"patterns": [{"dna": "A" * 12, "id": "w"}]}, "truth": {},
              "expect": [{"check": "status", "args": {"http": 200}},
                         {"check": "count_relation", "args": {"value": 1}}]}
        res = checks.evaluate(sc, ok_answer({"patterns": [{"id": "w"}]}), "mcp")
        self.assertEqual([(r["check"], r["outcome"], r["property"]) for r in res],
                         [("status", "pass", False), ("count_relation", "not_exposed", True)])
        self.assertEqual(checks.scenario_outcome(res), "inconclusive")


V2 = os.path.join(SUITE, "manifest", "suite_v2.json")


def v2_scenarios():
    with open(V2) as f:
        return dict((x["id"], x) for x in json.load(f)["scenarios"])


class StatedReasons(unittest.TestCase):
    """VERIFY-fixes.md X2: a check is inconclusive on a stated hiding, stop or time limitation only once the answer
    shows it can have happened; each reason a check accepts (README, "Stated reasons") with a case where it holds
    and one where it cannot have happened."""

    LIMITS = {"max_contexts": 100, "max_paths": 10, "max_anchors": 1000, "max_steps": 1000, "time_budget_ms": 60000,
              "finalize_reserve_ms": 250, "stop_at_threshold": False, "max_predicate_contexts": 50}

    def problems(self, counts, withheld=None, cut=None, stop=None, det="full", steps=1000, limits=None,
                 elapsed=1000.0, req=None, mode=None, index=None, predicate_echo=None, **extra):
        e = {"id": "w", "mode": mode or ("partial" if cut else "all_or_count"), "counts": counts,
             "withheld": {"reason": withheld} if withheld else None, "cut": {"reason": cut} if cut else None,
             "stop": {"phase": stop[0], "reason": stop[1]} if stop else None, "determinism": det, "results": [],
             "returned": 0, "work": {"steps": steps}}
        e.update(extra)
        lim = dict(self.LIMITS)
        lim.update(limits or {})
        body = {"patterns": [e], "limits": lim, "timing": {"elapsed_ms": elapsed}}
        if index is not None:
            body["index"] = index
        if predicate_echo is not None:
            body["predicate"] = predicate_echo
        sc = {"route": "POST /pattern", "request": req or {"patterns": [{"dna": "A" * 25, "id": "w"}]}, "truth": {}}
        return checks._pattern_reason_problems(sc, body, 0, e)

    def ok(self, *a, **k):
        p = self.problems(*a, **k)
        self.assertEqual(p, [], (a, k))

    def bad(self, needle, *a, **k):
        p = self.problems(*a, **k)
        self.assertTrue(p and needle in p[0], (needle, p, a, k))

    @staticmethod
    def c(rel, value=None, **kw):
        d = {"relation": rel, "value": value}
        d.update(kw)
        return d

    def test_withheld_reasons(self):
        ctx = lambda rel, v=None, **kw: {"contexts": self.c(rel, v, **kw)}  # noqa: E731
        self.ok(ctx("exact", 101), withheld="count_above_threshold")
        self.ok(ctx("bounds", 90, lower=90, upper=120), withheld="count_above_threshold")
        self.bad("not above the threshold 100", ctx("exact", 100), withheld="count_above_threshold")
        self.bad("not above the threshold 100", ctx("at_least", 500), withheld="count_above_threshold")
        long_req = {"patterns": [{"dna": "A" * 40, "id": "w"}], "long_search": "paths"}
        self.ok({"anchors": self.c("exact", 2), "paths": self.c("exact", 11)}, withheld="count_above_threshold",
                req=long_req)
        self.bad("not above the threshold 10", {"anchors": self.c("exact", 2), "paths": self.c("exact", 3)},
                 withheld="count_above_threshold", req=long_req)
        self.ok(ctx("at_least", 101), withheld="threshold_crossed", stop=("discovery", "max_contexts"),
                limits={"stop_at_threshold": True})
        self.bad("without stop_at_threshold", ctx("at_least", 101), withheld="threshold_crossed")
        self.ok(ctx("at_least", 9), withheld="discovery_budget", stop=("discovery", "max_steps"))
        self.bad("without a max_steps stop", ctx("at_least", 9), withheld="discovery_budget")
        self.bad("after 10 steps of max_steps 1000", ctx("at_least", 9), withheld="discovery_budget",
                 stop=("discovery", "max_steps"), steps=10)
        self.ok(ctx("unknown"), withheld="deadline", stop=("discovery", "time"), det="time_limited", elapsed=59900.0)
        self.bad("without a time stop", ctx("exact", 5), withheld="deadline")
        self.bad("before the work time", ctx("unknown"), withheld="deadline", stop=("discovery", "time"),
                 det="time_limited", elapsed=300.0)
        self.ok({"anchors": self.c("exact", 2), "paths": self.c("unknown")}, withheld="paths_later_increment",
                req={"patterns": [{"dna": "A" * 40, "id": "w"}]})
        self.bad("L <= k", ctx("exact", 5), withheld="paths_later_increment")
        self.bad("long_search 'paths'", {"anchors": self.c("exact", 2), "paths": self.c("exact", 1)},
                 withheld="paths_later_increment", req=long_req)
        self.ok({"anchors": self.c("exact", 1001), "paths": self.c("unknown")}, withheld="anchors_above_threshold",
                req=long_req)
        self.bad("not above max_anchors 1000", {"anchors": self.c("exact", 2), "paths": self.c("unknown")},
                 withheld="anchors_above_threshold", req=long_req)
        labels = {"patterns": [{"dna": "A" * 25, "id": "w"}], "output": {"labels": "all"}}
        self.ok(ctx("exact", 24), withheld="annotation_budget", stop=("label_discovery", "max_annotation_work"),
                req=labels)
        self.ok(ctx("exact", 24), withheld="annotation_budget", rows_refused=[{"row": 1}], req=labels)
        self.bad("without an annotation stop", ctx("exact", 24), withheld="annotation_budget", req=labels)
        self.ok(ctx("exact", 24), withheld="anchor_labels_truncated", anchors_truncated=[{"cap": 2, "total": 9}],
                req=labels)
        self.bad("without anchors_truncated", ctx("exact", 24), withheld="anchor_labels_truncated", req=labels)
        self.ok(ctx("exact", 24), withheld="output_budget", stop=("output", "max_memory"), req=labels)
        self.bad("without a memory stop", ctx("exact", 24), withheld="output_budget", req=labels)
        pred = {"patterns": [{"dna": "A" * 25, "id": "w"}], "predicate": {"any": ["562"]}}
        self.ok(ctx("exact", 60), withheld="predicate_above_threshold", req=pred)
        self.bad("not above max_predicate_contexts 50", ctx("exact", 10), withheld="predicate_above_threshold", req=pred)
        self.ok({"contexts": self.c("exact", 500), "selected": self.c("exact", 101)}, withheld="selected_above_threshold",
                req=pred)
        self.bad("not above 100", {"contexts": self.c("exact", 500), "selected": self.c("exact", 5)},
                 withheld="selected_above_threshold", req=pred)
        self.ok(ctx("exact", 24), withheld="predicate_budget", stop=("selection", "max_predicate_work"), req=pred)
        self.bad("without a predicate stop", ctx("exact", 24), withheld="predicate_budget", req=pred)
        self.bad("does not list", ctx("exact", 24), withheld="no_such_reason")

    def test_withheld_reasons_need_their_request(self):
        # a reason the request cannot have, whatever else the answer states; the label-free short pattern of
        # P-scope-0001 hides its source behind none of them
        ctx = lambda rel, v=None, **kw: {"contexts": self.c(rel, v, **kw)}  # noqa: E731
        refused = [{"row": 1, "phase": "label_discovery", "reason": "max_memory"}]
        self.bad("reads no annotation", ctx("exact", 19), withheld="annotation_budget", rows_refused=refused)
        self.bad("reads no annotation", ctx("exact", 19), withheld="annotation_budget",
                 stop=("label_discovery", "max_annotation_work"))
        self.bad("reads no annotation", ctx("exact", 19), rows_refused=refused)
        self.bad("without labels all", ctx("exact", 19), withheld="anchor_labels_truncated",
                 anchors_truncated=[{"cap": 2, "total": 9}])
        self.bad("without a memory account", ctx("exact", 19), withheld="output_budget",
                 stop=("output", "max_memory"))
        self.bad("without a predicate", ctx("exact", 19), withheld="predicate_budget",
                 stop=("selection", "max_predicate_work"))
        self.bad("without a predicate", ctx("exact", 60), withheld="predicate_above_threshold")
        self.bad("without a predicate", {"contexts": self.c("exact", 500), "selected": self.c("exact", 101)},
                 withheld="selected_above_threshold")
        pred = {"patterns": [{"dna": "A" * 25, "id": "w"}], "predicate": {"any": ["562"]}}
        self.bad("with a predicate on L <= k", ctx("exact", 101), withheld="count_above_threshold", req=pred)
        # a constant `true` normal form is the unfiltered request under the engine's own max_contexts (§19.6); the
        # claim needs both the constant pass and the normal form `true` (a self-contradictory answer cannot make it)
        constant = {"pass": "constant", "support": "kmer", "access": "rows"}
        both = {"contexts": self.c("exact", 101), "selected": self.c("exact", 101)}
        self.ok(both, withheld="count_above_threshold", req=pred, selection=constant,
                predicate_echo={"normal_form": True})
        for nf in (False, {"any": ["562"]}, None):
            self.bad("with a predicate on L <= k", both, withheld="count_above_threshold", req=pred, selection=constant,
                     predicate_echo={"normal_form": nf})
        self.bad("with a predicate on L <= k", both, withheld="count_above_threshold", req=pred, selection=constant)
        self.bad("with a predicate on L <= k", ctx("exact", 101), withheld="count_above_threshold", req=pred,
                 selection={"pass": "not_admitted"})
        # discovery_budget: the count the step stop interrupted is not exact (§7.4)
        self.bad("not exact", ctx("exact", 19), withheld="discovery_budget", stop=("discovery", "max_steps"))
        self.ok(ctx("bounds", 2, lower=2, upper=24), withheld="discovery_budget", stop=("mask_scan", "max_steps"))
        # only all_or_count withholds (patterns longer than k: also partial's anchors_above_threshold)
        self.bad("in mode partial", ctx("unknown"), withheld="deadline", stop=("discovery", "time"),
                 det="time_limited", elapsed=59900.0, mode="partial")
        long_req = {"patterns": [{"dna": "A" * 40, "id": "w"}], "long_search": "paths"}
        self.ok({"anchors": self.c("exact", 1001), "paths": self.c("unknown")}, withheld="anchors_above_threshold",
                req=long_req, mode="partial")
        # paths_later_increment needs anchors to extend; anchors_above_threshold a path search
        self.bad("anchors exact 0", {"anchors": self.c("exact", 0), "paths": self.c("exact", 0)},
                 withheld="paths_later_increment", req={"patterns": [{"dna": "A" * 40, "id": "w"}]})
        self.bad("without a path search", {"anchors": self.c("exact", 1001), "paths": self.c("unknown")},
                 withheld="anchors_above_threshold", req={"patterns": [{"dna": "A" * 40, "id": "w"}]})

    def test_threshold_stops(self):
        # a stop_at_threshold stop leaves the count it interrupted not exact (§7.4) and, on a masked BASIC or
        # CANONICAL graph, past its threshold: discovery stops once the running lower bound passes it, and the
        # at_least value is that bound (pattern_search.cpp threshold_crossed / running_lower / tally)
        ctx = lambda rel, v=None, **kw: {"contexts": self.c(rel, v, **kw)}  # noqa: E731
        thr = {"stop_at_threshold": True}
        masked = {"graph_mode": "basic"}
        self.ok(ctx("at_least", 101), withheld="threshold_crossed", stop=("discovery", "max_contexts"), limits=thr,
                index=masked)
        self.ok(ctx("at_least", 4), withheld="threshold_crossed", stop=("discovery", "max_contexts"),
                limits=dict(thr, max_contexts=3), index={"graph_mode": "canonical"})   # P-modes-0005
        self.bad("not exact", ctx("exact", 0), withheld="threshold_crossed", stop=("discovery", "max_contexts"),
                 limits=thr, index=masked)                                    # P-modes-0005 with its counts exact 0
        self.bad("not exact", ctx("exact", 0), withheld="threshold_crossed", stop=("discovery", "max_contexts"),
                 limits=thr)
        self.bad("not past max_contexts 100", ctx("at_least", 2), withheld="threshold_crossed",
                 stop=("discovery", "max_contexts"), limits=thr, index=masked)
        self.bad("not past max_contexts 100", ctx("at_least", 100), cut="max_contexts",
                 stop=("discovery", "max_contexts"), limits=thr, index=masked)
        # the value rule is not judged where the bound need not be the value
        self.ok(ctx("at_least", 0), withheld="threshold_crossed", stop=("discovery", "max_contexts"), limits=thr,
                index={"graph_mode": "basic", "counting": "upper_bound"})   # fixture unmasked_stop_at_threshold
        self.ok(ctx("at_least", 2), withheld="threshold_crossed", stop=("discovery", "max_contexts"), limits=thr,
                index={"graph_mode": "primary"})                            # an even-k wrapped PRIMARY union
        self.ok(ctx("at_least", 2), withheld="threshold_crossed", stop=("discovery", "max_contexts"), limits=thr)
        self.bad("without a stop_at_threshold stop", ctx("at_least", 101), withheld="threshold_crossed", limits=thr)
        self.bad("without a stop_at_threshold stop", ctx("at_least", 101), withheld="threshold_crossed",
                 stop=("discovery", "max_steps"), limits=thr)
        # with a predicate: the raw discovery's stop at max_predicate_contexts, the selection's at max_contexts
        pred = {"patterns": [{"dna": "A" * 25, "id": "w"}], "predicate": {"any": ["562"]}}
        self.ok({"contexts": self.c("at_least", 51), "selected": self.c("unknown")}, withheld="threshold_crossed",
                stop=("discovery", "max_predicate_contexts"), limits=thr, req=pred, index=masked)
        self.bad("not past max_predicate_contexts 50", {"contexts": self.c("at_least", 50),
                                                        "selected": self.c("unknown")},
                 withheld="threshold_crossed", stop=("discovery", "max_predicate_contexts"), limits=thr, req=pred,
                 index=masked)
        self.ok({"contexts": self.c("exact", 500), "selected": self.c("bounds", 101, lower=101, upper=300)},
                withheld="threshold_crossed", stop=("selection", "max_contexts"), limits=thr, req=pred, index=masked)
        self.bad("not exact", {"contexts": self.c("exact", 500), "selected": self.c("exact", 101)},
                 withheld="threshold_crossed", stop=("selection", "max_contexts"), limits=thr, req=pred)
        # a constant normal form runs the engine's own threshold in discovery: the raw count, never the selected
        self.ok({"contexts": self.c("at_least", 101), "selected": self.c("exact", 0)}, withheld="threshold_crossed",
                stop=("discovery", "max_contexts"), limits=thr, req=pred, index=masked)
        self.bad("without a predicate", ctx("at_least", 101), stop=("discovery", "max_predicate_contexts"),
                 limits=thr, mode="count")
        # max_anchors: the threshold of a pattern longer than k, with or without a path search (§7.7, §12.1)
        anchors_only = {"patterns": [{"dna": "A" * 40, "id": "w"}]}
        self.ok({"anchors": self.c("at_least", 1001), "paths": self.c("unknown")}, withheld="paths_later_increment",
                stop=("discovery", "max_anchors"), limits=thr, req=anchors_only)
        self.bad("not exact", {"anchors": self.c("exact", 2), "paths": self.c("unknown")},
                 withheld="paths_later_increment", stop=("discovery", "max_anchors"), limits=thr, req=anchors_only)
        self.bad("L <= k", ctx("at_least", 101), stop=("discovery", "max_anchors"), limits=thr, mode="count")
        # max_paths: the extension's threshold, a path search only; its paths at_least (§12.1)
        long_req = {"patterns": [{"dna": "A" * 40, "id": "w"}], "long_search": "paths"}
        self.ok({"anchors": self.c("exact", 2), "paths": self.c("at_least", 11)}, withheld="threshold_crossed",
                stop=("extension", "max_paths"), limits=thr, req=long_req)
        self.bad("not exact", {"anchors": self.c("exact", 2), "paths": self.c("exact", 11)},
                 withheld="threshold_crossed", stop=("extension", "max_paths"), limits=thr, req=long_req)
        self.bad("without a path search", {"anchors": self.c("exact", 2), "paths": self.c("at_least", 11)},
                 withheld="threshold_crossed", stop=("extension", "max_paths"), limits=thr, req=anchors_only)
        self.bad("longer than k", {"anchors": self.c("at_least", 2), "paths": self.c("unknown")},
                 stop=("discovery", "max_contexts"), limits=thr, req=long_req, mode="count")

    def test_cut_reasons(self):
        ctx = lambda rel, v=None, **kw: {"contexts": self.c(rel, v, **kw)}  # noqa: E731
        self.ok(ctx("exact", 101), cut="max_contexts")
        self.bad("within max_contexts 100", ctx("exact", 50), cut="max_contexts")
        self.ok(ctx("at_least", 9), cut="max_steps", stop=("discovery", "max_steps"))
        self.bad("without a max_steps stop", ctx("at_least", 9), cut="max_steps")
        self.ok(ctx("at_least", 0), cut="time", stop=("discovery", "time"), det="time_limited", elapsed=59900.0)
        self.ok(ctx("at_least", 9), cut="time", stop=("discovery", "max_steps"), det="time_limited", elapsed=59900.0)
        self.bad("without a time stop", ctx("exact", 9), cut="time")
        labels = {"patterns": [{"dna": "A" * 25, "id": "w"}], "output": {"labels": "all"}}
        self.ok(ctx("exact", 2000), cut="max_memory", stop=("output", "max_memory"), req=labels)
        self.bad("without a memory stop", ctx("exact", 2000), cut="max_memory", req=labels)
        self.bad("without a memory account", ctx("exact", 2000), cut="max_memory", stop=("output", "max_memory"))
        self.bad("without a predicate", ctx("exact", 1828), cut="max_predicate_contexts")
        self.bad("without a path search", ctx("exact", 11), cut="max_paths")
        self.bad("without a supported-path search", ctx("exact", 24), cut="max_annotation_work",
                 stop=("label_discovery", "max_annotation_work"), req=labels)
        long_req = {"patterns": [{"dna": "A" * 40, "id": "w"}], "long_search": "paths"}
        self.ok({"anchors": self.c("exact", 2), "paths": self.c("exact", 11)}, cut="max_paths", req=long_req)
        self.bad("within max_paths 10", {"anchors": self.c("exact", 2), "paths": self.c("exact", 3)}, cut="max_paths",
                 req=long_req)
        sp = {"patterns": [{"dna": "A" * 40, "id": "w"}], "long_search": "supported_paths"}
        self.ok({"anchors": self.c("exact", 2), "paths": self.c("exact", 2), "supported_paths": self.c("at_least", 1)},
                cut="max_annotation_work", stop=("extension", "max_annotation_work"), req=sp)
        self.bad("without that stop", {"anchors": self.c("exact", 2), "paths": self.c("exact", 2)},
                 cut="max_annotation_work", req=sp)
        pred = {"patterns": [{"dna": "A" * 25, "id": "w"}], "predicate": {"any": ["562"]}}
        self.ok({"contexts": self.c("exact", 1828), "selected": self.c("bounds", 6, lower=6, upper=1734)},
                cut="max_predicate_contexts", req=pred)
        self.bad("within 50", {"contexts": self.c("exact", 10), "selected": self.c("exact", 2)},
                 cut="max_predicate_contexts", req=pred)
        self.ok(ctx("exact", 1828), cut="max_predicate_work", stop=("selection", "max_predicate_work"), req=pred)
        self.bad("without that stop", ctx("exact", 1828), cut="max_predicate_work", req=pred)
        self.ok({"anchors": self.c("at_least", 1001)}, cut="max_anchors", stop=("discovery", "max_anchors"),
                limits={"stop_at_threshold": True}, req=long_req)
        self.bad("without that stop", {"anchors": self.c("exact", 2)}, cut="max_anchors", req=long_req)
        self.bad("does not list", ctx("exact", 2), cut="no_such_reason")

    def test_stop_reasons(self):
        ctx = {"contexts": self.c("at_least", 3)}
        self.ok(ctx, stop=("discovery", "time"), det="time_limited", elapsed=59900.0, mode="count")
        self.bad("before the work time", ctx, stop=("discovery", "time"), det="time_limited", elapsed=300.0,
                 mode="count")
        self.ok(ctx, stop=("discovery", "max_steps"), mode="count")
        self.bad("a max_steps stop after 999 steps", ctx, stop=("discovery", "max_steps"), steps=999, mode="count")
        self.ok(ctx, stop=("discovery", "max_contexts"), limits={"stop_at_threshold": True}, mode="count")
        self.bad("without stop_at_threshold", ctx, stop=("discovery", "max_contexts"), mode="count")
        labels = {"patterns": [{"dna": "A" * 25, "id": "w"}], "output": {"labels": "all"}}
        self.ok(ctx, stop=("label_discovery", "max_annotation_work"), mode="partial", req=labels)
        self.bad("reads no annotation", ctx, stop=("label_discovery", "max_annotation_work"), mode="count")
        # mode count reads no annotation whatever the projection (§14.1) and is the label-free path (§7.6)
        self.bad("reads no annotation", ctx, stop=("label_discovery", "max_annotation_work"), mode="count",
                 req=labels)
        self.bad("reads no annotation", ctx, mode="count", req=labels, rows_refused=[{"row": 1}])
        self.ok(ctx, stop=("selection", "max_predicate_work"), mode="count",
                req={"patterns": [{"dna": "A" * 25, "id": "w"}], "predicate": {"any": ["562"]}})
        self.bad("without a predicate", ctx, stop=("selection", "max_predicate_work"), mode="count")
        self.ok(ctx, stop=("output", "max_memory"), mode="partial", req=labels)
        self.bad("without a memory account", ctx, stop=("output", "max_memory"), mode="count", req=labels)
        self.ok(ctx, stop=("selection", "max_memory"), mode="count",
                req={"patterns": [{"dna": "A" * 25, "id": "w"}], "predicate": {"any": ["562"]}})
        self.ok({"anchors": self.c("exact", 2), "paths": self.c("at_least", 1)}, stop=("extension", "max_memory"),
                mode="count", req={"patterns": [{"dna": "A" * 40, "id": "w"}], "long_search": "paths"})
        self.bad("without a memory account", ctx, stop=("output", "max_memory"), mode="count")
        self.bad("does not list", ctx, stop=("discovery", "no_such_reason"), mode="count")
        # a low-complexity diagnostic the work time cut (stop null, determinism time_limited)
        self.ok({"contexts": self.c("exact", 3)}, det="time_limited", elapsed=59900.0, mode="count")
        self.bad("a time limitation after 300 ms", {"contexts": self.c("exact", 3)}, det="time_limited",
                 elapsed=300.0, mode="count")
        # the work time is a lower bound: E counts the results the answer withholds at a generous size (9,000
        # plain results: about 5.5 s; labelled ones: more than the budget, so any time stop may be one)
        self.ok({"contexts": self.c("exact", 9000)}, withheld="deadline", stop=("extraction", "time"),
                det="time_limited", elapsed=54300.0, limits={"max_contexts": 10000})
        self.bad("before the work time", {"contexts": self.c("exact", 9000)}, withheld="deadline",
                 stop=("extraction", "time"), det="time_limited", elapsed=36000.0, limits={"max_contexts": 10000})
        self.ok({"contexts": self.c("exact", 9000)}, withheld="deadline", stop=("extraction", "time"),
                det="time_limited", elapsed=36000.0, limits={"max_contexts": 10000},
                req={"patterns": [{"dna": "A" * 25, "id": "w"}], "output": {"labels": "all"}})

    def pattern_body(self, sid, entry_update, elapsed):
        """A /pattern answer of scenario |sid| (v2 manifest) with one entry of exact counts, updated."""
        sc = v2_scenarios()[sid]
        req = sc["request"]
        p = req["patterns"][0]
        L = 3 * len(p["protein"]) if "protein" in p else len(p.get("dna") or p.get("iupac"))
        offs = [str(o) for o in range(0, 31 - L + 1)]
        e = {"id": p["id"], "mode": req["mode"], "scope": req["scope"], "stop": None, "determinism": "full",
             "notes": [], "results": [], "returned": 0, "retrieval_complete": False, "cut": None,
             "withheld": None, "work": {"steps": 1000, "ranges_visited": 1000},
             "counts": {"contexts": {"relation": "exact", "value": 0,
                                     "by_offset": dict((o, {"relation": "exact", "value": 0}) for o in offs)}}}
        e.update(entry_update)
        body = {"patterns": [e], "limits": {"max_contexts": req["max_contexts"], "max_steps": 100000000,
                                            "time_budget_ms": req["time_budget_ms"], "finalize_reserve_ms": 250},
                "timing": {"elapsed_ms": elapsed}}
        return sc, ok_answer(body)

    @staticmethod
    def exact(n, offs):
        per = [n // len(offs)] * len(offs)
        per[0] += n - sum(per)
        return {"contexts": {"relation": "exact", "value": n, "by_offset": dict(
            (o, {"relation": "exact", "value": v}) for o, v in zip(offs, per))}}

    def test_verifier_examples_pattern(self):
        # X-peptide10-0009 (exact 2,810, cap 10,000) and N-realrep-0007 (exact 2) withheld as count_above_threshold
        for sid, n, offs in (("X-peptide10-0009", 2810, ["0", "1"]), ("N-realrep-0007", 2, ["0"])):
            sc, a = self.pattern_body(sid, {"counts": self.exact(n, offs),
                                            "withheld": {"reason": "count_above_threshold"}}, 100.0)
            res = dict((r["check"], r) for r in checks.evaluate(sc, a, "direct-local"))
            self.assertEqual(res["contains_source_context"]["outcome"], "fail", sid)
            self.assertIn("not above the threshold 10000", res["contains_source_context"]["detail"])
            self.assertEqual(res["invariants"]["outcome"], "fail", sid)
            self.assertEqual(checks.scenario_outcome(list(res.values())), "fail")
        # N-hard-0005: exact 2 withheld as count_above_threshold at max_contexts 100 (its v2 withheld check allows
        # the reason): the invariants fail it
        sc, a = self.pattern_body("N-hard-0005", {"counts": self.exact(2, ["0"]), "notes": ["low_complexity_pattern"],
                                                  "withheld": {"reason": "count_above_threshold"}}, 100.0)
        res = dict((r["check"], r) for r in checks.evaluate(sc, a, "direct-local"))
        self.assertEqual(res["withheld"]["outcome"], "pass")
        self.assertEqual(res["invariants"]["outcome"], "fail")
        self.assertEqual(checks.scenario_outcome(list(res.values())), "fail")
        # the 300 ms pattern stop: a discovery time stop after 300 ms of a 60,000 ms budget
        unknown = {"contexts": {"relation": "unknown", "value": None, "by_offset": {
            "0": {"relation": "unknown", "value": None}, "1": {"relation": "unknown", "value": None}}}}
        stopped = {"counts": unknown, "stop": {"phase": "discovery", "reason": "time"}, "determinism": "time_limited",
                   "withheld": {"reason": "deadline"}}
        sc, a = self.pattern_body("X-peptide10-0009", stopped, 300.0)
        res = dict((r["check"], r) for r in checks.evaluate(sc, a, "direct-local"))
        self.assertEqual(res["contains_source_context"]["outcome"], "fail")
        self.assertIn("before the work time", res["contains_source_context"]["detail"])
        self.assertEqual(res["invariants"]["outcome"], "fail")
        sc, a = self.pattern_body("X-peptide10-0009", stopped, 59800.0)      # the same stop at its work time
        res = dict((r["check"], r) for r in checks.evaluate(sc, a, "direct-local"))
        self.assertEqual(res["contains_source_context"]["outcome"], "inconclusive")
        self.assertEqual(res["invariants"]["outcome"], "pass", res["invariants"]["detail"])
        # count_relation accepts a budget stop only when it happened
        sc = {"route": "POST /pattern", "request": {"patterns": [{"dna": "A" * 25, "id": "w"}], "time_budget_ms": 60000},
              "truth": {}}
        for el, want in ((300.0, "fail"), (59900.0, "inconclusive")):
            a = ok_answer({"patterns": [{"id": "w", "counts": {"contexts": {"relation": "at_least", "value": 3}},
                                         "stop": {"phase": "discovery", "reason": "time"}, "determinism": "time_limited"}],
                           "limits": {"time_budget_ms": 60000, "finalize_reserve_ms": 250}, "timing": {"elapsed_ms": el}})
            self.assertEqual(outcome(sc, a, "count_relation", relation="exact", value=5)[0], want)

    def walk_answer(self, sc, n, trig, counters, limitation, bounds=None, timing=None, resource_stop=None):
        """A walk of scenario |sc| (v2 manifest) following its record for n bp on the requested arm, stopped by
        |trig|."""
        t = sc["truth"]
        side = sc["request"]["strategy"]["direction"]
        cont = checks._seed_frame(t, {"start": t["start"], "flank": t["flank"], "flank_start": t["flank_start"]},
                                  side, t["length"])[0]
        arm = {"status": "truncated", "complete_to_bp": n, "cap_trigger": {"reason": trig, "at_bp": n},
               "segments": [{"id": 0, "parents": [], "from_bp": 0, "length_bp": n, "sequence": cont(n)}],
               "paths": [{"id": 0, "segments": [0], "length_bp": n, "path_reason": trig, "end_reasons": {trig: 1}}],
               "counters": counters, "limitations": [limitation] if limitation else []}
        b = dict(sc["request"]["strategy"]["bounds"], max_steps=1000000)
        b.update(bounds or {})
        r = {"outcome": {"walks": "partial"}, "arms": {side: arm},
             "timing": timing or {"elapsed_ms": 2000.0, "seed_phase_ms": 500.0,
                                  "deadline": {"longest_piece": {"kind": "read", "ms": 100.0, "rows": 1}}}}
        if resource_stop:
            r["resource_stop"] = resource_stop
        return ok_answer({"results": [r], "strategy": {"bounds": b}})

    def test_verifier_examples_traverse(self):
        sc = v2_scenarios()["T-walk-0002"]                          # a left walk under the record's header, 1,000 bp
        wd = lambda knob, limit, observed: {"kind": "walk_domain", "knob": knob, "limit": limit,  # noqa: E731
                                            "observed": observed}
        # the 500-of-500,000 bp walk: "truncated by max_output_bp" at 500 bp under the request's 500,000
        a = self.walk_answer(sc, 500, "max_output_bp", {"output_bp": 500, "steps": 500},
                             wd("bounds.max_output_bp", 500000, 501))
        o, d = outcome(sc, a, "walk_in_record", label_scope="record")
        self.assertEqual(o, "fail")
        self.assertIn("observed 501 within the limit 500000", d)
        self.assertEqual(outcome(sc, a, "invariants")[0], "fail")
        a = self.walk_answer(sc, 500, "max_output_bp", {"output_bp": 500, "steps": 500}, None)
        self.assertIn("no walk_domain limitation names bounds.max_output_bp",
                      outcome(sc, a, "walk_in_record", label_scope="record")[1])
        # the same stop under a max_output_bp of 500 (the counter at the cap, the demand past it): inconclusive
        a = self.walk_answer(sc, 500, "max_output_bp", {"output_bp": 500, "steps": 500},
                             wd("bounds.max_output_bp", 500, 501), bounds={"max_output_bp": 500})
        self.assertEqual(outcome(sc, a, "walk_in_record", label_scope="record")[0], "inconclusive")
        self.assertEqual(outcome(sc, a, "invariants")[0], "pass")
        # "by max_steps" with the request's max_steps 1,000,000 far away
        a = self.walk_answer(sc, 500, "max_steps", {"output_bp": 500, "steps": 500},
                             wd("bounds.max_steps", 1000000, 1000001))
        o, d = outcome(sc, a, "walk_in_record", label_scope="record")
        self.assertEqual(o, "fail")
        self.assertIn("after 500 steps of max_steps 1000000", d)
        # a time stop 2 s into a 15 s budget, and one after it
        early = {"elapsed_ms": 2000.0, "seed_phase_ms": 500.0, "deadline": {
            "longest_piece": {"kind": "read", "ms": 100.0, "rows": 1}, "stopped_by": "time_budget",
            "stop_after_deadline_ms": 0.0}}
        a = self.walk_answer(sc, 500, "time_budget", {"output_bp": 500, "steps": 500},
                             wd("bounds.time_budget_ms", 15000, 2000.0), timing=early,
                             resource_stop={"resource": "time", "used": 2000.0})
        o, d = outcome(sc, a, "walk_in_record", label_scope="record")
        self.assertEqual(o, "fail")
        self.assertIn("before the 15000 ms budget", d)
        late = dict(early, elapsed_ms=15301.0, deadline=dict(early["deadline"], stop_after_deadline_ms=300.0))
        a = self.walk_answer(sc, 500, "time_budget", {"output_bp": 500, "steps": 500},
                             wd("bounds.time_budget_ms", 15000, 15300.0), timing=late,
                             resource_stop={"resource": "time", "used": 15300.0})
        self.assertEqual(outcome(sc, a, "walk_in_record", label_scope="record")[0], "inconclusive")
        # walk_end: a stop that cannot have happened fails, one that can is inconclusive
        sc = v2_scenarios()["T-end-0001"]
        a = self.walk_answer(sc, 30, "time_budget", {"output_bp": 30, "steps": 30},
                             wd("bounds.time_budget_ms", 60000, 2000.0), timing=early,
                             resource_stop={"resource": "time", "used": 2000.0})
        self.assertEqual(outcome(sc, a, "walk_end", arm="right", reasons_in=["dead_end"])[0], "fail")
        late = dict(late, elapsed_ms=60301.0)
        a = self.walk_answer(sc, 30, "time_budget", {"output_bp": 30, "steps": 30},
                             wd("bounds.time_budget_ms", 60000, 60300.0), timing=late,
                             resource_stop={"resource": "time", "used": 60300.0})
        self.assertEqual(outcome(sc, a, "walk_end", arm="right", reasons_in=["dead_end"])[0], "inconclusive")

    def test_mcp_shaped_stops(self):
        # G1 of the verification of lane S3: a graphlet result (the MCP's export, and `detail: graphlet`) carries its
        # arms' cap_trigger, limitations and counters as JSON beside the MGT text; they are judged there too
        sc = v2_scenarios()["T-walk-0002"]
        wd = lambda knob, limit, observed: {"kind": "walk_domain", "knob": knob, "limit": limit,  # noqa: E731
                                            "observed": observed}
        # the 500-of-500,000 bp walk "truncated by max_output_bp", with and without the walk_domain limitation
        for lim in (wd("bounds.max_output_bp", 500000, 501), None):
            full = self.walk_answer(sc, 500, "max_output_bp", {"output_bp": 500, "steps": 500}, lim)
            mcp = mcp_shaped(full)
            r = mcp["body"]["results"][0]
            self.assertIn("graphlet", r)
            self.assertNotIn("segments", r["arms"]["left"])
            for arm in ("mcp", "direct-local"):
                o, d = outcome(sc, mcp, "walk_in_record", _arm=arm, label_scope="record")
                self.assertEqual(o, "fail", (arm, lim, d))
                self.assertEqual(outcome(sc, mcp, "invariants", _arm=arm)[0], "fail", (arm, lim))
        # the same stop where it happened (max_output_bp 500): inconclusive, the invariants pass
        full = self.walk_answer(sc, 500, "max_output_bp", {"output_bp": 500, "steps": 500},
                                wd("bounds.max_output_bp", 500, 501), bounds={"max_output_bp": 500})
        for arm in ("mcp", "direct-local"):
            self.assertEqual(outcome(sc, mcp_shaped(full), "walk_in_record", _arm=arm, label_scope="record")[0],
                             "inconclusive", arm)
            self.assertEqual(outcome(sc, mcp_shaped(full), "invariants", _arm=arm)[0], "pass", arm)
        # a time stop 2 s into a 15 s budget in MCP shape (E4): walk_in_record fails, not only budget_respected
        early = {"elapsed_ms": 2000.0, "seed_phase_ms": 500.0, "deadline": {
            "longest_piece": {"kind": "read", "ms": 100.0, "rows": 1}, "stopped_by": "time_budget",
            "stop_after_deadline_ms": 0.0}}
        full = self.walk_answer(sc, 500, "time_budget", {"output_bp": 500, "steps": 500},
                                wd("bounds.time_budget_ms", 15000, 2000.0), timing=early,
                                resource_stop={"resource": "time", "used": 2000.0})
        o, d = outcome(sc, mcp_shaped(full), "walk_in_record", _arm="mcp", label_scope="record")
        self.assertEqual(o, "fail", d)
        self.assertIn("before the 15000 ms budget", d)
        # the MGT's own A record states the cap trigger where the JSON arm does not
        bare = mcp_shaped(full, json_arms=False)
        arms = checks.graphlet_arms(bare["body"]["results"][0]["graphlet"])
        self.assertEqual(arms["left"]["cap_trigger"]["reason"], "time_budget")
        self.assertEqual(arms["left"]["counters"]["output_bp"], 500)

    def test_resource_end_needs_a_cap_trigger(self):
        # G2 of the verification of lane S3: a path ended by a resource reason lies on an arm that cap stopped
        sc = v2_scenarios()["T-walk-0002"]
        wd = lambda knob, limit, observed: {"kind": "walk_domain", "knob": knob, "limit": limit,  # noqa: E731
                                            "observed": observed}
        # E7: a path ended by max_steps at 500 bp on an arm `complete` without a cap trigger
        a = self.walk_answer(sc, 500, "max_steps", {"output_bp": 500, "steps": 500}, None)
        arm = a["body"]["results"][0]["arms"]["left"]
        arm["status"], arm["complete_to_bp"] = "complete", 1000
        del arm["cap_trigger"]
        a["body"]["results"][0]["outcome"]["walks"] = "complete"
        for shape in (a, mcp_shaped(a)):
            o, d = outcome(sc, shape, "walk_in_record", label_scope="record")
            self.assertEqual(o, "fail", d)
            self.assertIn("a path ended by max_steps on an arm complete", d)
            o, d = outcome(sc, shape, "invariants")
            self.assertEqual(o, "fail", d)
        # truncated, but with no cap trigger
        arm["status"], arm["complete_to_bp"] = "truncated", 500
        a["body"]["results"][0]["outcome"]["walks"] = "partial"
        o, d = outcome(sc, a, "walk_in_record", label_scope="record")
        self.assertEqual(o, "fail", d)
        self.assertIn("without a resource cap trigger", d)
        # a path ended by max_output_bp on an arm whose cap trigger is time_budget: only as a later cap that the
        # arm's walk_domain limitations name (SPEC §7.0)
        late = {"elapsed_ms": 15301.0, "seed_phase_ms": 500.0, "deadline": {
            "longest_piece": {"kind": "read", "ms": 100.0, "rows": 1}, "stopped_by": "time_budget",
            "stop_after_deadline_ms": 300.0}}
        a = self.walk_answer(sc, 500, "max_output_bp", {"output_bp": 500, "steps": 500},
                             wd("bounds.time_budget_ms", 15000, 15300.0), timing=late,
                             resource_stop={"resource": "time", "used": 15300.0}, bounds={"max_output_bp": 500})
        a["body"]["results"][0]["arms"]["left"]["cap_trigger"]["reason"] = "time_budget"
        o, d = outcome(sc, a, "walk_in_record", label_scope="record")
        self.assertEqual(o, "fail", d)
        self.assertIn("no walk_domain limitation names a later cap", d)
        a["body"]["results"][0]["arms"]["left"]["limitations"].append(wd("bounds.max_output_bp", 500, 501))
        self.assertEqual(outcome(sc, a, "walk_in_record", label_scope="record")[0], "inconclusive")
        self.assertEqual(outcome(sc, a, "invariants")[0], "pass")
        # walk_end: a resource end on an arm without a cap trigger fails
        sc = v2_scenarios()["T-end-0001"]
        a = self.walk_answer(sc, 30, "max_steps", {"output_bp": 30, "steps": 30}, None)
        del a["body"]["results"][0]["arms"]["right"]["cap_trigger"]
        o, d = outcome(sc, a, "walk_end", arm="right", reasons_in=["dead_end"])
        self.assertEqual(o, "fail", d)
        self.assertIn("without a resource cap trigger", d)

    def test_every_stored_pattern_answer_states_possible_reasons(self):
        # the fixtures (real answers of the mini index) state only reasons that happened
        n = 0
        for name in sorted(os.listdir(PATTERN)):
            pa, pq = os.path.join(PATTERN, name, "answer.json"), os.path.join(PATTERN, name, "request.json")
            if not (os.path.exists(pa) and os.path.exists(pq)):
                continue
            with open(pa) as f:
                ans = json.load(f)
            with open(pq) as f:
                req = json.load(f)
            if not isinstance(ans, dict) or not isinstance(req, dict) or "patterns" not in ans:
                continue
            sc = {"route": "POST /pattern", "request": req, "truth": {}}
            for i, e in enumerate(ans["patterns"]):
                self.assertEqual(checks._pattern_reason_problems(sc, ans, i, e), [], name)
                n += 1
        self.assertGreater(n, 150)


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
        # a's property is its stop (R-budget, resolve_runs prefix_if_stopped): two complete profiles that
        # agree did not test it; one that disagrees still fails; the continuation keeps its comparison
        # (REL-prefix_equals-0002 / -0004)
        budget = {"route": "POST /resolve", "request": {}, "truth": {},
                  "expect": [{"check": "resolve_runs", "args": {"label": "562", "full": True, "prefix_if_stopped": True}}]}
        rel = self.rel("prefix_equals", labels=["562"], a_from_kmer=0, b_from_kmer=0)
        self.assertEqual(self.ev(rel, stopped, full, sa=budget)[0], "pass")
        o, d = self.ev(rel, full, full, sa=budget)
        self.assertEqual(o, "inconclusive")
        self.assertIn("no stop happened on a: the prefix property was not tested", d)
        self.assertEqual(self.ev(rel, self.resolve(100, [[0, 100]], [[0, 41], [60, 100]]), full, sa=budget)[0], "fail")
        self.assertEqual(self.ev(rel, full, full)[0], "pass")

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
        # exact runs (R-junction): one run more than stated fails, as a label covering the other's k-mers would
        self.assertEqual(outcome(sc, part, "resolve_runs", label="562", runs=[[0, 30], [40, 70]])[0], "pass")
        self.assertEqual(outcome(sc, part, "resolve_runs", label="562", runs=[[0, 30]])[0], "fail")
        self.assertEqual(outcome(sc, part, "resolve_runs", label="562", run_end_in=[30, 30])[0], "pass")
        self.assertEqual(outcome(sc, part, "resolve_runs", label="562", max_kmers=59)[0], "fail")
        self.assertEqual(outcome(sc, part, "resolve_runs", label="NOPE", full=True)[0], "fail")
        broken = self.body([[0, 30]], graph=[[0, 20]])
        self.assertEqual(outcome(sc, broken, "invariants")[0], "fail")
        wrongnum = self.body([[0, 69]], num=69, graph=[[0, 69]])
        self.assertEqual(outcome(sc, wrongnum, "invariants")[0], "fail")

    def test_stop_prefix_and_absent(self):
        sc = {"route": "POST /resolve", "request": {"sequence": "A" * 100, "labels": ["562", "9606"],
                                                    "bounds": {"time_budget_ms": 3000}}, "truth": {}}
        stop = {"phase": "rows", "reason": "time", "resolved_kmers": 40, "query_kmers": 70, "resolved_bp": 70,
                "query_bp": 100, "remainder_from_bp": 40}
        a = self.body([[0, 40]], num=40, graph=[[0, 40]], stop=stop)
        a["body"].update(limits={"time_budget_ms": 3000, "finalize_reserve_ms": 250}, timing={"elapsed_ms": 2800.0})
        # a stop before the budget less the reserve cannot have happened (SPEC §4.5): fail, not a prefix
        early = copy.deepcopy(a)
        early["body"]["timing"]["elapsed_ms"] = 2000.0
        o, d = outcome(sc, early, "resolve_runs", label="562", full=True, prefix_if_stopped=True)
        self.assertEqual(o, "fail")
        self.assertIn("before the budget 3000 ms less the reserve 250 ms", d)
        self.assertEqual(outcome(sc, early, "invariants")[0], "fail")
        self.assertEqual(outcome(sc, a, "invariants")[0], "pass")
        self.assertEqual(outcome(sc, a, "resolve_runs", label="562", full=True, unless_stopped=True)[0], "pass")
        self.assertEqual(outcome(sc, a, "resolve_runs", label="562", full=True, prefix_if_stopped=True)[0], "pass")
        # a stop that claims the record's own k-mers unsupported, or absent from the graph, fails
        late = {"limits": {"time_budget_ms": 3000, "finalize_reserve_ms": 250}, "timing": {"elapsed_ms": 2800.0}}
        hole = self.body([[0, 30]], num=40, graph=[[0, 40]], stop=stop)
        hole["body"].update(late)
        self.assertEqual(outcome(sc, hole, "resolve_runs", label="562", full=True, prefix_if_stopped=True)[0], "fail")
        gone = self.body([], num=40, graph=[], stop=stop)
        gone["body"].update(late)
        self.assertEqual(outcome(sc, gone, "resolve_runs", label="562", full=True, prefix_if_stopped=True)[0], "fail")
        zero = self.body([], num=0, graph=[], stop=dict(stop, resolved_kmers=0, resolved_bp=0, remainder_from_bp=0))
        zero["body"].update(late)
        self.assertEqual(outcome(sc, zero, "resolve_runs", label="562", full=True, prefix_if_stopped=True)[0],
                         "inconclusive")
        bad = copy.deepcopy(a)
        bad["body"]["stop"]["remainder_from_bp"] = 41
        self.assertEqual(outcome(sc, bad, "invariants")[0], "fail")
        # R-budget's property is the stop; a complete answer holds the profile but never tested the prefix a
        # stop answers (inconclusive), a wrong complete profile still fails; a twin's property
        # (unless_stopped) is the complete profile itself
        complete = self.body([[0, 70]])
        complete["body"].update(limits={"time_budget_ms": 3000, "finalize_reserve_ms": 250}, stop=None)
        o, d = outcome(sc, complete, "resolve_runs", label="562", full=True, prefix_if_stopped=True)
        self.assertEqual(o, "inconclusive")
        self.assertIn("no stop happened: the prefix property was not tested", d)
        holed = self.body([[0, 30]])
        self.assertEqual(outcome(sc, holed, "resolve_runs", label="562", full=True, prefix_if_stopped=True)[0], "fail")
        self.assertEqual(outcome(sc, complete, "resolve_runs", label="562", full=True, unless_stopped=True)[0], "pass")
        empty = ok_answer({"k": 31, "num_kmers": 70, "graph_runs": [], "labels": [
            {"label": "562", "kind": "column", "kmers_supported": 0, "runs": []}]})
        self.assertEqual(outcome(sc, empty, "resolve_absent")[0], "pass")
        self.assertEqual(outcome(sc, self.body([[0, 70]]), "resolve_absent")[0], "fail")

    def test_runs_on_a_stopped_answer_without_a_prefix_argument(self):
        # a stop answers the resolve of a prefix (SPEC-labeled-traversal-core §4.5); an expectation about the whole
        # query (no prefix_if_stopped / unless_stopped) is decided only where the prefix decides it
        sc = {"route": "POST /resolve", "request": {"sequence": "A" * 100, "labels": ["562", "9606"],
                                                    "bounds": {"time_budget_ms": 20000}}, "truth": {}}
        late = {"limits": {"time_budget_ms": 20000, "finalize_reserve_ms": 250}, "timing": {"elapsed_ms": 19900.0}}

        def stopped(x, runs):
            stop = {"phase": "rows", "reason": "time", "resolved_kmers": x, "query_kmers": 70,
                    "resolved_bp": x + 30 if x else 0, "query_bp": 100, "remainder_from_bp": x}
            a = self.body(runs, num=x, graph=[[0, x]] if x else [], stop=stop)
            a["body"].update(late)
            return a
        prefix = stopped(40, [[0, 40]])
        self.assertEqual(outcome(sc, prefix, "invariants")[0], "pass")
        o, d = outcome(sc, prefix, "resolve_runs", label="562", full=True)
        self.assertEqual(o, "inconclusive")
        self.assertIn("full support of the query untested", d)
        self.assertEqual(outcome(sc, stopped(40, [[0, 20], [25, 40]]), "resolve_runs", label="562", full=True)[0],
                         "fail")
        self.assertEqual(outcome(sc, stopped(0, []), "resolve_runs", label="562", full=True)[0], "inconclusive")
        holed = stopped(40, [[0, 20], [25, 40]])
        for args, want in (({"covers": [0, 20]}, "pass"), ({"covers": [10, 30]}, "fail"),
                           ({"covers": [30, 60]}, "inconclusive"), ({"covers": [10, 60]}, "fail"),
                           ({"covers": [45, 60]}, "inconclusive"),
                           ({"run_end_in": [15, 25]}, "pass"), ({"run_end_in": [30, 35]}, "fail"),
                           ({"run_start_in": [20, 30]}, "inconclusive"), ({"run_start_in": [0, 10]}, "fail"),
                           ({"runs": [[0, 20], [25, 70]]}, "inconclusive"), ({"runs": [[0, 20], [30, 70]]}, "fail"),
                           ({"kmers": 60}, "inconclusive"), ({"kmers": 30}, "fail"), ({"kmers": 70}, "fail"),
                           ({"max_kmers": 30}, "fail"), ({"max_kmers": 50}, "inconclusive"),
                           ({"min_kmers": 30}, "pass"), ({"min_kmers": 66}, "fail"),
                           ({"min_kmers": 50}, "inconclusive"),
                           ({"covers": [0, 20], "kmers": 60}, "inconclusive")):
            self.assertEqual(outcome(sc, holed, "resolve_runs", label="562", **args)[0], want, args)
        self.assertEqual(outcome(sc, prefix, "resolve_runs", label="562", run_end_in=[30, 50])[0], "inconclusive")
        self.assertEqual(outcome(sc, prefix, "resolve_runs", label="562", run_end_in=[10, 30])[0], "fail")

    def test_discover(self):
        sc = {"route": "POST /resolve", "request": {"sequence": "A" * 100, "discover": {}}, "truth": {}}
        a = self.body([[0, 70]])
        self.assertEqual(outcome(sc, a, "resolve_discover_contains", label="562")[0], "pass")
        self.assertEqual(outcome(sc, a, "resolve_discover_contains", label="1")[0], "fail")
        # SPEC-labeled-traversal-core §4.2: a truncated list keeps the top max_labels by supported k-mers; the
        # source supports every in-graph k-mer, so it may be dropped only in a tie at that support (the recorded
        # R-discover-0001: 713600 kept, the source 107806 one of 7 dropped at full length)
        disc = {"route": "POST /resolve", "request": {"sequence": "A" * 100, "discover": {"kind": "column",
                                                                                          "max_labels": 1}},
                "truth": {}}

        def listed(label, n, trunc):
            return ok_answer({"k": 31, "num_kmers": 70, "graph_runs": [[0, 70]], "labels_truncated": trunc,
                              "labels": [{"label": label, "kind": "column", "kmers_supported": n, "runs": [[0, n]]}]})
        tie = listed("713600", 70, {"kept": 1, "total": 9, "min_kept_kmers": 70, "max_dropped_kmers": 70,
                                    "dropped_full_length": 7})
        o, d = outcome(disc, tie, "resolve_discover_contains", label="107806")
        self.assertEqual(o, "pass")
        self.assertIn("ties the kept labels at all 70 in-graph k-mers", d)
        self.assertNotIn("index order", d)
        self.assertEqual(outcome(disc, tie, "invariants")[0], "pass")
        self.assertEqual(outcome(disc, tie, "resolve_discover_contains", label="713600")[0], "pass")
        # an index-order host: a label with 3 of 70 k-mers kept, one with 70 dropped
        index_order = listed("670", 3, {"kept": 1, "total": 9, "min_kept_kmers": 3, "max_dropped_kmers": 70,
                                        "dropped_full_length": 1})
        self.assertEqual(outcome(disc, index_order, "resolve_discover_contains", label="107806")[0], "fail")
        self.assertEqual(outcome(disc, index_order, "invariants")[0], "fail")
        # the kept label below full support, or no dropped label at full support: the source would be kept
        below = listed("670", 60, {"kept": 1, "total": 9, "min_kept_kmers": 60, "max_dropped_kmers": 60,
                                   "dropped_full_length": 0})
        o, d = outcome(disc, below, "resolve_discover_contains", label="107806")
        self.assertEqual(o, "fail")
        self.assertIn("supports all 70 in-graph k-mers and a kept label only 60", d)
        self.assertEqual(outcome(disc, below, "invariants")[0], "pass")
        none_full = listed("713600", 70, {"kept": 1, "total": 9, "min_kept_kmers": 70, "max_dropped_kmers": 69,
                                          "dropped_full_length": 0})
        self.assertEqual(outcome(disc, none_full, "resolve_discover_contains", label="107806")[0], "fail")
        self.assertEqual(outcome(disc, none_full, "invariants")[0], "pass")
        # labels_truncated without its fields: not exposed on the MCP arm, a failure on a direct one
        bare = listed("713600", 70, {"kept": 1, "total": 9})
        self.assertEqual(outcome(disc, bare, "resolve_discover_contains", label="107806")[0], "fail")
        self.assertEqual(outcome(disc, bare, "resolve_discover_contains", label="107806", _arm="mcp")[0],
                         "not_exposed")
        # the invariants: kept is the list's length, below total; the last kept label states min_kept_kmers; a
        # full-length drop has every in-graph k-mer; the list is ranked by supported k-mers, at most max_labels
        for trunc in ({"kept": 2, "total": 9, "min_kept_kmers": 70, "max_dropped_kmers": 70, "dropped_full_length": 7},
                      {"kept": 1, "total": 1, "min_kept_kmers": 70, "max_dropped_kmers": 70, "dropped_full_length": 0},
                      {"kept": 1, "total": 9, "min_kept_kmers": 69, "max_dropped_kmers": 69, "dropped_full_length": 0},
                      {"kept": 1, "total": 9, "min_kept_kmers": 70, "max_dropped_kmers": 60, "dropped_full_length": 2},
                      {"kept": 1, "total": 3, "min_kept_kmers": 70, "max_dropped_kmers": 70, "dropped_full_length": 3}):
            self.assertEqual(outcome(disc, listed("713600", 70, trunc), "invariants")[0], "fail", trunc)
        # only more than max_labels labels are truncated, to exactly max_labels (§4.2; default 1000): a list cut to
        # another size is no tie, whatever its counts say (R-discover-0003's shape, where a correct host lists all 18)
        def cut(max_labels, kept, total):
            sc_ = dict(disc, request=dict(disc["request"], discover={"kind": "column", "max_labels": max_labels}))
            if max_labels is None:
                sc_["request"]["discover"] = {"kind": "column"}
            ans = ok_answer({"k": 31, "num_kmers": 70, "graph_runs": [[0, 70]], "labels": [
                {"label": str(1000 + j), "kind": "column", "kmers_supported": 70, "runs": [[0, 70]]}
                for j in range(kept)], "labels_truncated": {"kept": kept, "total": total, "min_kept_kmers": 70,
                                                            "max_dropped_kmers": 70, "dropped_full_length": total - kept}})
            return sc_, ans
        for max_labels, kept, total in ((100, 1, 18), (20, 5, 18), (None, 1, 9), (20, 20, 20)):
            sc_, ans = cut(max_labels, kept, total)
            o, d = outcome(sc_, ans, "resolve_discover_contains", label="272634")
            self.assertEqual(o, "fail", (max_labels, kept, total))
            self.assertIn("only more than max_labels labels are truncated", d)
            o, d = outcome(sc_, ans, "invariants")
            self.assertEqual(o, "fail", (max_labels, kept, total))
            self.assertIn("truncates only more than max_labels labels", d)
        for max_labels, kept, total in ((20, 20, 25), (None, 1000, 1200)):
            sc_, ans = cut(max_labels, kept, total)
            self.assertEqual(outcome(sc_, ans, "resolve_discover_contains", label="272634")[0], "pass")
            self.assertEqual(outcome(sc_, ans, "invariants")[0], "pass")
        two = ok_answer({"k": 31, "num_kmers": 70, "graph_runs": [[0, 70]], "labels_truncated": None, "labels": [
            {"label": "a", "kind": "column", "kmers_supported": 10, "runs": [[0, 10]]},
            {"label": "b", "kind": "column", "kmers_supported": 70, "runs": [[0, 70]]}]})
        many = dict(disc, request=dict(disc["request"], discover={"kind": "column", "max_labels": 5}))
        o, d = outcome(many, two, "invariants")
        self.assertEqual(o, "fail")
        self.assertIn("not ranked by supported k-mers", d)
        self.assertEqual(outcome(disc, two, "invariants")[0], "fail")          # 2 labels, max_labels 1
        two["body"]["labels"].reverse()
        self.assertEqual(outcome(many, two, "invariants")[0], "pass")
        trace = dict(many, request=dict(many["request"], support="trace"))
        two["body"]["labels"].reverse()
        self.assertEqual(outcome(trace, two, "invariants")[0], "pass")         # trace states trace runs, not the rank


    def test_discover_stopped_before_an_in_graph_kmer(self):
        disc = {"route": "POST /resolve", "request": {"sequence": "A" * 100, "discover": {"kind": "column",
                                                                                          "max_labels": 1}},
                "truth": {}}
        # a discovery stopped before it read an in-graph k-mer has no profile to rank (§4.5); one that read some
        # k-mers and lists no source fails as a complete list does
        budgeted = dict(disc, request=dict(disc["request"], bounds={"time_budget_ms": 3000}))
        for resolved, graph in ((0, []), (5, [])):
            stopped = ok_answer({"k": 31, "num_kmers": resolved, "graph_runs": graph, "labels": [],
                                 "labels_truncated": None, "limits": {"time_budget_ms": 3000, "finalize_reserve_ms": 250},
                                 "timing": {"elapsed_ms": 2800.0},
                                 "stop": {"phase": "rows", "reason": "time", "resolved_kmers": resolved,
                                          "remainder_from_bp": resolved, "resolved_bp": resolved + 30 if resolved else 0,
                                          "query_kmers": 70, "query_bp": 100}})
            o, d = outcome(budgeted, stopped, "resolve_discover_contains", label="107806")
            self.assertEqual(o, "inconclusive", resolved)
            self.assertIn("no profile to rank", d)
        some = ok_answer({"k": 31, "num_kmers": 40, "graph_runs": [[0, 40]], "labels_truncated": None, "labels": [
            {"label": "713600", "kind": "column", "kmers_supported": 40, "runs": [[0, 40]]}],
            "limits": {"time_budget_ms": 3000, "finalize_reserve_ms": 250}, "timing": {"elapsed_ms": 2800.0},
            "stop": {"phase": "rows", "reason": "time", "resolved_kmers": 40, "remainder_from_bp": 40,
                     "resolved_bp": 70, "query_kmers": 70, "query_bp": 100}})
        self.assertEqual(outcome(budgeted, some, "resolve_discover_contains", label="107806")[0], "fail")


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

    def test_walk_judged_by_its_own_end_below_the_frontier(self):
        # a cap that truncated the arm stopped the record's route only where the route reaches
        # complete_to_bp; one that ended before it ended for its own stated reason
        # (SPEC-labeled-traversal-core §6.10), and where that is the record's base outside A, C, G, T
        # (T-walk-0018's 'S'), the truth's gap
        import random
        rnd = random.Random(17)
        record = "".join(rnd.choice("ACGT") for _ in range(1000))
        start, L = 400, 40

        def scen(rec, direction):
            t = {"record": "x", "taxid": "1", "seq_id": 0, "record_length": 1000, "start": start, "length": L,
                 "strand": "+", "window": rec[start:start + L], "flank": rec, "flank_start": 0}
            return {"route": "POST /traverse", "truth": t, "request": {
                "seeds": [{"sequence": t["window"]}],
                "strategy": {"direction": direction, "bounds": {"max_extension_bp": 100, "time_budget_ms": 15000}}}}
        right, left = record[start + L:start + L + 30], record[start - 30:start]
        stopped = dict(status="truncated", trig="time_budget", reason="label_lost")
        sc = scen(record, "right")
        o, d = outcome(sc, self.walk({"right": right}, **stopped), "walk_in_record")
        self.assertEqual(o, "inconclusive")      # the route reaches complete_to_bp (30): the cap may have ended it
        self.assertIn("stopped by time_budget", d)
        o, d = outcome(sc, self.walk({"right": right}, extra={"complete_to_bp": 60}, **stopped), "walk_in_record")
        self.assertEqual(o, "fail")
        self.assertIn("the record followed for only 30 of 100 bp, the path ended by ['label_lost']", d)
        # the record's next base an IUPAC code: no k-mer holds it, the route ends there (a truth gap)
        sc_s = scen(record[:start + L + 30] + "S" + record[start + L + 31:], "right")
        o, d = outcome(sc_s, self.walk({"right": right}, extra={"complete_to_bp": 60}, **stopped), "walk_in_record")
        self.assertEqual(o, "inconclusive")
        self.assertIn("before the record's base 'S', which no k-mer of the graph holds: a truth gap", d)
        o, d = outcome(sc_s, self.walk({"right": right}, reason="label_lost"), "walk_in_record")
        self.assertEqual(o, "inconclusive")      # also on a complete arm
        # a resource end stated below complete_to_bp cannot have happened (§6.10: every path that ended before it
        # ended for a semantic reason): the walk and the invariants fail, at the frontier it is a stated stop
        capped = dict(status="truncated", trig="time_budget", reason="time_budget")
        o, d = outcome(sc, self.walk({"right": right}, extra={"complete_to_bp": 60}, **capped), "walk_in_record")
        self.assertEqual(o, "fail")
        self.assertIn("below the arm's complete_to_bp 60", d)
        self.assertEqual(outcome(sc, self.walk({"right": right}, extra={"complete_to_bp": 60}, **capped),
                                 "invariants")[0], "fail")
        self.assertEqual(outcome(sc, self.walk({"right": right}, **capped), "walk_in_record")[0], "inconclusive")
        self.assertEqual(outcome(sc, self.walk({"right": right}, **capped), "invariants")[0], "pass")
        sc_l = scen(record[:start - 31] + "N" + record[start - 30:], "left")
        o, d = outcome(sc_l, self.walk({"left": left}, extra={"complete_to_bp": 60}, **stopped), "walk_in_record")
        self.assertEqual(o, "inconclusive")
        self.assertIn("'N'", d)
        self.assertEqual(outcome(scen(record, "left"), self.walk({"left": left}, extra={"complete_to_bp": 60},
                                                                 **stopped), "walk_in_record")[0], "fail")

    def merged_walk(self, record, start, L, other):
        """A right arm whose displayed path (first parents) passes a merge through another record's allele |other| at
        30 bp, while the record's own allele is the merge's second parent (SPEC §7.1, on_reconverge merge)."""
        cont = record[start + L:start + L + 100]
        segs = [{"id": 0, "parents": [], "from_bp": 0, "length_bp": 30, "sequence": cont[:30]},
                {"id": 1, "parents": [0], "from_bp": 30, "length_bp": 1, "sequence": other},
                {"id": 2, "parents": [0], "from_bp": 30, "length_bp": 1, "sequence": cont[30]},
                {"id": 3, "parents": [1, 2], "from_bp": 31, "length_bp": 69, "sequence": cont[31:],
                 "labels_via_parent": [[0], []],
                 "events": [{"at_bp": 100, "type": "label_end", "label": 0, "reason": "dead_end"}]}]
        arm = {"status": "complete", "complete_to_bp": 100, "segments": segs, "splits": [{"at_bp": 30, "segment": 0}],
               "paths": [{"id": 0, "segments": [0, 1, 3], "length_bp": 100, "path_reason": "max_extension_bp",
                          "end_reasons": {"max_extension_bp": 1}}]}
        return ok_answer({"results": [{"outcome": {"walks": "complete"}, "arms": {"right": arm}}],
                          "strategy": {"bounds": {"max_extension_bp": 100}}})

    def test_walk_in_record_takes_any_parent_at_a_merge(self):
        # G13 (as amended by VERIFY-T): the record's route through a non-first parent of a merge follows the record
        import random
        rnd = random.Random(11)
        record = "".join(rnd.choice("ACGT") for _ in range(1000))
        start, L = 400, 40
        t = {"record": "x", "taxid": "1", "seq_id": 0, "record_length": 1000, "start": start, "length": L,
             "strand": "+", "window": record[start:start + L], "flank": record, "flank_start": 0}
        req = {"seeds": [{"sequence": t["window"]}], "strategy": {"direction": "right",
                                                                   "bounds": {"max_extension_bp": 100}}}
        sc = {"route": "POST /traverse", "request": req, "truth": t}
        own = record[start + L + 30]
        other = "A" if own != "A" else "C"
        ans = self.merged_walk(record, start, L, other)
        o, d = outcome(sc, ans, "walk_in_record")
        self.assertEqual(o, "pass", d)
        self.assertIn("through a non-first parent of a merge", d)
        # the displayed path alone does not follow it (the check of before read only that)
        arm = ans["body"]["results"][0]["arms"]["right"]
        self.assertNotEqual(checks._arm_flanks(arm, "right")[0][1], record[start + L:start + L + 100])
        # with the record's allele on neither parent no route follows it: fail
        bad = self.merged_walk(record, start, L, other)
        bad["body"]["results"][0]["arms"]["right"]["segments"][2]["sequence"] = "G" if own not in "AG" else "T"
        self.assertEqual(outcome(sc, bad, "walk_in_record")[0], "fail")
        # the graphlet (the MCP's export) reads the same DAG
        o, d = outcome(sc, mcp_shaped(ans), "walk_in_record", _arm="mcp")
        self.assertEqual(o, "pass", d)
        # the route search is bounded: a search that stops before it shows the minimum is inconclusive, not fail
        frames = [checks._seed_frame(t, {"start": start, "flank": record, "flank_start": 0}, "right", L)]
        routes, beyond, capped = checks._record_routes(arm, "right", frames)
        self.assertEqual(([len(f) for _, f in routes], beyond, capped), ([100], 0, False))
        self.assertTrue(checks._record_routes(arm, "right", frames, cap=2)[2])
        saved = checks._MAX_ROUTE_STATES
        try:
            checks._MAX_ROUTE_STATES = 2
            o, d = outcome(sc, ans, "walk_in_record")
            self.assertEqual(o, "inconclusive", d)
            self.assertIn("route search stopped", d)
            self.assertEqual(outcome(sc, bad, "walk_in_record")[0], "inconclusive")
        finally:
            checks._MAX_ROUTE_STATES = saved
        # walk_end: the record's route to its end, its leaf's end reasons (the label_end event included)
        t_end = dict(t, record_length=start + L + 100)
        sc_end = dict(sc, truth=t_end)
        self.assertEqual(outcome(sc_end, ans, "walk_end", arm="right", reasons_in=["dead_end"])[0], "pass")
        self.assertEqual(outcome(sc_end, ans, "walk_end", arm="right", reasons_in=["label_lost"])[0], "fail")
        o, d = outcome(sc_end, bad, "walk_end", arm="right", reasons_in=["dead_end"])
        self.assertEqual(o, "fail")
        self.assertIn("no route through the segments follows the record", d)

    def labelled_merge(self, own_labels, other_labels, truncated=False):
        """merged_walk with segment labels: label 0 the record's header (label_dict), 1 another record's; the record's
        allele (segment 2) carries |own_labels|, the other allele (segment 1) |other_labels|."""
        import random
        rnd = random.Random(11)
        record = "".join(rnd.choice("ACGT") for _ in range(1000))
        start, L = 400, 40
        t = {"record": "NC_X.1", "taxid": "1", "seq_id": 0, "record_length": 1000, "start": start, "length": L,
             "strand": "+", "window": record[start:start + L], "flank": record, "flank_start": 0}
        req = {"seeds": [{"sequence": t["window"]}], "strategy": {"direction": "right",
                                                                   "bounds": {"max_extension_bp": 100}}}
        sc = {"route": "POST /traverse", "request": req, "truth": t}
        own = record[start + L + 30]
        ans = self.merged_walk(record, start, L, "A" if own != "A" else "C")
        r = ans["body"]["results"][0]
        r["label_dict"] = [{"name": "NC_X.1", "kind": "header"}, {"name": "NC_Y.1", "kind": "header"}]
        segs = r["arms"]["right"]["segments"]
        for sg, labs in zip(segs, ([0, 1], other_labels, own_labels, [0, 1])):
            sg["labels"] = labs
            sg["labels_total"] = len(labs)
            sg["labels_truncated"] = False
        if truncated:
            segs[2].update(labels_total=5, labels_truncated=True)
        return sc, ans

    def test_walk_routes_carry_the_record_label(self):
        # A-04: the record's sequence pieced together from another label's segment (its own allele carries only label
        # 1, the other allele label 0) does not follow the record; it passed on the sequence alone
        sc, ans = self.labelled_merge([1], [0])
        o, d = outcome(sc, ans, "walk_in_record", label_scope="derived")
        self.assertEqual(o, "fail", d)
        self.assertNotIn("sequence only", d)
        # the record's label on its own allele, or the own allele's list cut by the cap: the route follows the record
        for own, trunc in (([0, 1], False), ([1], True)):
            sc, ans = self.labelled_merge(own, [0], trunc)
            o, d = outcome(sc, ans, "walk_in_record", label_scope="derived")
            self.assertEqual(o, "pass", d)
        # the record not in the label dictionary (G14's cut), or no segment labels: the sequence alone, said so
        sc, ans = self.labelled_merge([1], [0])
        ans["body"]["results"][0]["label_dict"] = [{"name": "NC_Z.1", "kind": "header"}]
        o, d = outcome(sc, ans, "walk_in_record", label_scope="derived")
        self.assertEqual(o, "pass", d)
        self.assertIn("sequence only: the record's label is not in the seed's label dictionary", d)
        sc, ans = self.labelled_merge([1], [0])
        for sg in ans["body"]["results"][0]["arms"]["right"]["segments"]:
            del sg["labels"]
        self.assertIn("sequence only: the segments state no labels", outcome(sc, ans, "walk_in_record")[1])
        # a column label names the record too (its taxid)
        sc, ans = self.labelled_merge([1], [0])
        ans["body"]["results"][0]["label_dict"] = [{"name": "1", "kind": "column"}, {"name": "NC_Y.1", "kind": "header"}]
        self.assertEqual(outcome(sc, ans, "walk_in_record", label_scope="column")[0], "fail")
        # walk_end reads the same rule
        t_end = dict(sc["truth"], record_length=sc["truth"]["start"] + 40 + 100)
        sc, ans = self.labelled_merge([1], [0])
        self.assertEqual(outcome(dict(sc, truth=t_end), ans, "walk_end", arm="right", reasons_in=["dead_end"])[0], "fail")
        sc, ans = self.labelled_merge([0, 1], [0])
        self.assertEqual(outcome(dict(sc, truth=t_end), ans, "walk_end", arm="right", reasons_in=["dead_end"])[0], "pass")
        # the graphlet (the MCP's export): the G records' entry sets and the L records' names
        sc, ans = self.labelled_merge([1], [0])
        self.assertEqual(outcome(sc, mcp_shaped(ans), "walk_in_record", _arm="mcp", label_scope="derived")[0], "fail")
        sc, ans = self.labelled_merge([0, 1], [0])
        o, d = outcome(sc, mcp_shaped(ans), "walk_in_record", _arm="mcp", label_scope="derived")
        self.assertEqual(o, "pass", d)
        self.assertNotIn("sequence only", d)

    def test_route_search_is_memoised(self):
        # A-05: a DAG whose split children spell the same base (SPEC-violating) doubles the routes per diamond; the
        # memo on (segment, copies followed) searches each state once, so a divergent leaf fails instead of hitting the
        # bound (inconclusive before)
        import random
        rnd = random.Random(5)
        record = "".join(rnd.choice("ACGT") for _ in range(400))
        start, L = 100, 40
        t = {"record": "x", "taxid": "1", "seq_id": 0, "record_length": 400, "start": start, "length": L,
             "strand": "+", "window": record[start:start + L], "flank": record, "flank_start": 0}
        cont = record[start + L:]
        segs, last, bp = [{"id": 0, "parents": [], "from_bp": 0, "length_bp": 1, "sequence": cont[0]}], 0, 1
        for _ in range(20):
            a, b, m = len(segs), len(segs) + 1, len(segs) + 2
            segs += [{"id": a, "parents": [last], "from_bp": bp, "length_bp": 1, "sequence": cont[bp]},
                     {"id": b, "parents": [last], "from_bp": bp, "length_bp": 1, "sequence": cont[bp]},
                     {"id": m, "parents": [a, b], "from_bp": bp + 1, "length_bp": 1, "sequence": cont[bp + 1]}]
            last, bp = m, bp + 2
        wrong = "A" if cont[bp] != "A" else "C"
        segs.append({"id": len(segs), "parents": [last], "from_bp": bp, "length_bp": 1, "sequence": wrong})
        chain, cur = [len(segs) - 1], len(segs) - 1
        while segs[cur]["parents"]:
            cur = segs[cur]["parents"][0]
            chain.append(cur)
        arm = {"status": "complete", "complete_to_bp": bp + 1, "segments": segs,
               "paths": [{"id": 0, "segments": chain[::-1], "length_bp": bp + 1, "path_reason": "dead_end",
                          "end_reasons": {"dead_end": 1}}]}
        ans = ok_answer({"results": [{"outcome": {"walks": "complete"}, "arms": {"right": arm}}],
                         "strategy": {"bounds": {"max_extension_bp": 100}}})
        sc = {"route": "POST /traverse", "truth": t,
              "request": {"seeds": [{"sequence": t["window"]}], "strategy": {"direction": "right",
                                                                              "bounds": {"max_extension_bp": 100}}}}
        frames = [checks._seed_frame(t, {"start": start, "flank": record, "flank_start": 0}, "right", L)]
        routes, beyond, capped = checks._record_routes(arm, "right", frames, cap=200)
        self.assertEqual((routes, capped), ([], False))
        o, d = outcome(sc, ans, "walk_in_record")
        self.assertEqual(o, "fail", d)

    def repeat_walk(self, flank_end=None, record_scope=True):
        """N-realrep-0006's shape (G28): a seed inside a (GTTGGC)n repeat whose record's next base re-enters the seed
        (rejoined_seed at 0 bp, SPEC §6.6); the walk continues from the repeat's exit."""
        import random
        rnd = random.Random(3)
        pre = "".join(rnd.choice("ACGT") for _ in range(300))
        rep = ("GTTGGC" * 13)[:76]
        post = "".join(rnd.choice("ACT") for _ in range(1)) + "".join(rnd.choice("ACGT") for _ in range(400))
        record = pre + rep + post
        start, L = 306, 40
        flank = record if flank_end is None else record[:flank_end]
        t = {"record": "NC_R.1", "taxid": "1", "seq_id": 0, "record_length": len(record), "start": start,
             "length": L, "strand": "+", "window": record[start:start + L], "flank": flank, "flank_start": 0,
             "walk": {"right": {"min_bp": 0, "branch_by_bp": 0, "first": "ambiguous"}}}
        path = record[376:476]
        arm = {"status": "complete", "complete_to_bp": 100,
               "segments": [{"id": 0, "parents": [], "from_bp": 0, "length_bp": 100, "sequence": path,
                             "events": [{"at_bp": 0, "type": "blocked", "char": record[start + L],
                                         "reason": "rejoined_seed"},
                                        {"at_bp": 100, "type": "label_end", "label": 0, "reason": "max_extension_bp"}]}],
               "paths": [{"id": 0, "segments": [0], "length_bp": 100, "path_reason": "max_extension_bp",
                          "end_reasons": {"max_extension_bp": 1}}]}
        ans = ok_answer({"capabilities": {"regime": "basic"}, "results": [{"outcome": {"walks": "complete"},
                                                                           "arms": {"right": arm}}],
                         "strategy": {"bounds": {"max_extension_bp": 100}}})
        sc = {"route": "POST /traverse", "truth": t,
              "request": {"seeds": [{"sequence": t["window"]}], "strategy": {"direction": "right",
                                                                              "bounds": {"max_extension_bp": 100}}}}
        return sc, ans, record

    def structural_walk(self, end="branch", branches=2, text=None, blocked=True, kinds=("ambiguous", "ambiguous"),
                        strategy=None, echo=None, events=True):
        """(scenario, answer): a right arm whose route along the record passes two splits (|kinds|) at 10 and 20 bp
        and ends at 30 bp by the label end |end| (a label_end event, with |text| beside it; the lineage's |branches|
        in end_labels; a blocked event of a re-entry when |blocked|); each split's other child is one base off the
        record, ended dead_end. The request's strategy adds |strategy|, the answer's echo is |echo|."""
        import random
        rnd = random.Random(11)
        record = "".join(rnd.choice("ACGT") for _ in range(1000))
        start, L = 400, 40
        t = {"record": "x", "taxid": "1", "seq_id": 0, "record_length": 1000, "start": start, "length": L,
             "strand": "+", "window": record[start:start + L], "flank": record, "flank_start": 0}
        st = {"direction": "right", "bounds": {"max_extension_bp": 100}}
        for group, knobs in (strategy or {}).items():
            if isinstance(knobs, dict):
                st.setdefault(group, {}).update(knobs)
            else:
                st[group] = knobs
        sc = {"route": "POST /traverse", "request": {"seeds": [{"sequence": t["window"]}], "strategy": st},
              "truth": t}
        cont = record[start + L:start + L + 100]
        off = lambda b: "A" if b != "A" else "C"  # noqa: E731
        ev = {"at_bp": 30, "type": "label_end", "label": 0, "reason": end, "structural_successors": 2}
        if text:
            ev["text"] = text
        leaf_events = [ev]
        if end in ("edge_reuse", "edge_reuse_rc", "rejoined_seed") and blocked:
            leaf_events.insert(0, {"at_bp": 30, "type": "blocked", "char": cont[30], "reason": end})
        segs = [{"id": 0, "parents": [], "from_bp": 0, "length_bp": 10, "sequence": cont[:10], "events": []},
                {"id": 1, "parents": [0], "from_bp": 10, "length_bp": 10, "sequence": cont[10:20], "events": []},
                {"id": 2, "parents": [0], "from_bp": 10, "length_bp": 1, "sequence": off(cont[10]),
                 "events": [{"at_bp": 11, "type": "label_end", "label": 0, "reason": "dead_end",
                             "structural_successors": 0}]},
                {"id": 3, "parents": [1], "from_bp": 20, "length_bp": 10, "sequence": cont[20:30],
                 "events": leaf_events if events else []},
                {"id": 4, "parents": [1], "from_bp": 20, "length_bp": 1, "sequence": off(cont[20]),
                 "events": [{"at_bp": 21, "type": "label_end", "label": 0, "reason": "dead_end",
                             "structural_successors": 0}]}]
        counted = "branch" if end in checks._QUORUM_ENDS or text else end
        arm = {"status": "complete", "complete_to_bp": 100, "segments": segs,
               "splits": [{"at_bp": 10, "segment": 0, "children": [1, 2], "kind": kinds[0]},
                          {"at_bp": 20, "segment": 1, "children": [3, 4], "kind": kinds[1]}],
               "paths": [{"id": 0, "segments": [0, 1, 3], "length_bp": 30, "end_reasons": {counted: 1},
                          "end_labels": [{"label": 0, "branches": branches, "loss": 0.0}]},
                         {"id": 1, "segments": [0, 1, 4], "length_bp": 21, "end_reasons": {"dead_end": 1},
                          "end_labels": [{"label": 0, "branches": branches, "loss": 0.0}]},
                         {"id": 2, "segments": [0, 2], "length_bp": 11, "end_reasons": {"dead_end": 1},
                          "end_labels": [{"label": 0, "branches": 1, "loss": 0.0}]}]}
        body = {"results": [{"outcome": {"walks": "complete"}, "arms": {"right": arm}}],
                "strategy": echo if echo is not None else {"bounds": {"max_extension_bp": 100}}}
        return sc, ok_answer(body)

    def test_structural_ends_need_their_conditions(self):
        # SPEC-labeled-traversal-core §6.4 / §6.6: a structural end short of the radius passes walk_in_record only
        # where the request, the strategy and the arm show it can have happened; the invariants fail otherwise
        def judge(sc, ans):
            return outcome(sc, ans, "walk_in_record")[0], outcome(sc, ans, "invariants")[0]
        mlb2 = {"branching": {"max_label_branches": 2}}
        self.assertEqual(judge(*self.structural_walk(strategy=mlb2)), ("pass", "pass"))
        # the allowance unspent (an off-by-one, a lost lineage counter), or spent on no ambiguous split
        for kw in ({"branches": 0}, {"branches": 1}, {"kinds": ("ambiguous", "divergence")}):
            sc, ans = self.structural_walk(strategy=mlb2, **kw)
            o, d = outcome(sc, ans, "walk_in_record")
            self.assertEqual((o, outcome(sc, ans, "invariants")[0]), ("fail", "fail"), kw)
            self.assertIn("a stated stop that cannot have happened", d)
        # the request's allowance over the echo: a service that dropped the override is not excused by its echo
        sc, ans = self.structural_walk(strategy=mlb2, branches=0, kinds=("divergence", "divergence"),
                                       echo={"branching": {"max_label_branches": 0}})
        self.assertEqual(judge(sc, ans), ("fail", "fail"))
        sc, ans = self.structural_walk(branches=0, kinds=("divergence", "divergence"),
                                       echo={"branching": {"max_label_branches": 0}})
        self.assertEqual(judge(sc, ans), ("pass", "pass"))     # no request value: the echo's 0, nothing to spend
        sc, ans = self.structural_walk(strategy=mlb2, echo={"branching": {"max_label_branches": 0},
                                                            "clamped": [{"field": "branching.max_label_branches",
                                                                         "requested": 2, "effective": 0}]},
                                       branches=0, kinds=("divergence", "divergence"))
        self.assertEqual(judge(sc, ans), ("pass", "pass"))     # a stated clamp is the allowance
        self.assertEqual(judge(*self.structural_walk(strategy={"branching": {"max_label_branches": "unlimited"}})),
                         ("fail", "fail"))
        self.assertEqual(judge(*self.structural_walk(strategy={"exhaustive": True})), ("fail", "fail"))
        # the request's preset over the echo: an exhaustive request whose echo says max_label_branches 2 (a service
        # that rewrote the preset) is held to the preset's unlimited allowance, as is annotate mode
        for preset in ({"exhaustive": True}, {"labels": {"mode": "annotate"}}):
            self.assertEqual(judge(*self.structural_walk(strategy=preset, echo={"branching": {
                "max_label_branches": 2}})), ("fail", "fail"), preset)
        self.assertEqual(judge(*self.structural_walk(echo={"exhaustive": False, "branching": {
            "max_label_branches": 2}})), ("pass", "pass"))
        # the quorum and split-limit ends need their knob, in either form of the event
        for end, text in (("minority", None), ("branch", "minority")):
            self.assertEqual(judge(*self.structural_walk(end=end, text=text)), ("fail", "fail"))
            self.assertEqual(judge(*self.structural_walk(end=end, text=text, strategy={
                "branching": {"min_successor_labels": 2}})), ("pass", "pass"))
            self.assertEqual(judge(*self.structural_walk(end=end, text=text, strategy={
                "branching": {"min_successor_fraction": 0.2}})), ("pass", "pass"))
        self.assertEqual(judge(*self.structural_walk(end="below_min_labels")), ("fail", "fail"))
        self.assertEqual(judge(*self.structural_walk(end="below_min_labels", strategy={
            "bounds": {"min_live_labels": 2}})), ("pass", "pass"))
        self.assertEqual(judge(*self.structural_walk(end="split_limit")), ("fail", "fail"))   # 2 of 64 splits
        self.assertEqual(judge(*self.structural_walk(end="split_limit", strategy={
            "branching": {"max_splits_per_path": 2}})), ("pass", "pass"))
        self.assertEqual(judge(*self.structural_walk(end="split_limit", strategy={
            "branching": {"max_splits_per_path": 3}})), ("fail", "fail"))
        # a re-entry states its blocked event where the label ended
        for end in ("edge_reuse", "edge_reuse_rc", "rejoined_seed"):
            self.assertEqual(judge(*self.structural_walk(end=end)), ("pass", "pass"))
            self.assertEqual(judge(*self.structural_walk(end=end, blocked=False)), ("fail", "fail"))
        # a structural reason the path counts without its label_end event at the leaf
        self.assertEqual(judge(*self.structural_walk(strategy=mlb2, events=False)), ("fail", "fail"))

    def test_branch_taken_without_a_split_under_a_quorum(self):
        # under a quorum knob the walker raises a lineage's branches on every admissible successor before the
        # quorum drops the minority ones; one followed successor records no split. Label 0 is ambiguous at 4 bp
        # within its allowance (max_label_branches 1), its other successor a minority: no split; at 10 bp it is
        # ambiguous again over its allowance (branch, branches 1), label 1's only successor a minority
        def walk(msl, branches=1):
            armo = {"status": "complete", "complete_to_bp": 10,
                    "segments": [{"id": 0, "parents": [], "from_bp": 0, "length_bp": 10, "sequence": "ACGTACGTAC",
                                  "events": [{"at_bp": 4, "type": "label_end", "label": 0, "reason": "minority",
                                              "structural_successors": 2},
                                             {"at_bp": 10, "type": "label_end", "label": 0, "reason": "branch",
                                              "structural_successors": 2},
                                             {"at_bp": 10, "type": "label_end", "label": 1, "reason": "minority",
                                              "structural_successors": 2}]}],
                    "splits": [],
                    "paths": [{"id": 0, "segments": [0], "length_bp": 10, "end_reasons": {"branch": 2},
                               "end_labels": [{"label": 0, "branches": branches, "loss": 0.0},
                                              {"label": 1, "branches": 0, "loss": 0.0}]}]}
            strategy = {"direction": "right", "bounds": {"max_extension_bp": 10},
                        "branching": {"max_label_branches": 1, "min_successor_labels": msl}}
            sc = {"route": "POST /traverse", "truth": {},
                  "request": {"seeds": [{"sequence": "A" * 31}], "strategy": strategy}}
            return sc, ok_answer({"results": [{"outcome": {"walks": "complete"}, "arms": {"right": armo}}],
                                  "strategy": strategy})
        self.assertEqual(outcome(*walk(2), "invariants")[0], "pass")
        o, d = outcome(*walk(2, branches=0), "invariants")
        self.assertEqual(o, "fail")                 # the allowance is still the lineage's count at the end
        self.assertIn("a branch end after 0 branch(es) of max_label_branches 1", d)
        sc, ans = walk(1)                           # no quorum: every branch taken is a split, and no minority
        o, d = outcome(sc, ans, "invariants")
        self.assertEqual(o, "fail")
        self.assertIn("0 ambiguous split(s)", d)

    def test_structural_knobs_of_the_mcp_arm(self):
        # the MCP arm asks for its preset with the overrides of the scenario's MCP arguments; the preset mirrors the
        # scenario's direct request. Its knobs are those overrides, else the direct request, never the echo alone:
        # a service that lowers the preset's max_label_branches and echoes what it walked (T-walk-0017's MCP answer
        # with the echo and the branches lowered together) fails, named as a preset mismatch
        def mcp(branches, echo_mlb, overrides=None):
            sc, ans = self.structural_walk(strategy={"branching": {"max_label_branches": 2}}, branches=branches,
                                           echo={"branching": {"max_label_branches": echo_mlb}})
            sc["mcp"] = {"tool": "traverse_fetch", "arguments": {"preset": "walk", "overrides": dict(
                {"bounds": {"max_extension_bp": 100}, "direction": "right"}, **(overrides or {}))}}
            return sc, ans
        for check in ("invariants", "walk_in_record"):
            self.assertEqual(outcome(*mcp(2, 2), check, _arm="mcp")[0], "pass", check)
            for lowered in (1, 0):
                o, d = outcome(*mcp(lowered, lowered), check, _arm="mcp")
                self.assertEqual(o, "fail", (check, lowered))
                self.assertIn("the service walked with branching.max_label_branches %d, the scenario declares 2: a "
                              "preset mismatch" % lowered, d)
            # an override the scenario states is the MCP arm's own declaration; one the service did not apply fails
            self.assertEqual(outcome(*mcp(1, 1, {"branching": {"max_label_branches": 1}}), check, _arm="mcp")[0],
                             "pass", check)
            o, d = outcome(*mcp(1, 1, {"branching": {"max_label_branches": 2}}), check, _arm="mcp")
            self.assertEqual(o, "fail", check)
            self.assertIn("the scenario's MCP overrides state 2: an override not applied", d)
        # a direct arm is held to its request; its echo differing is judged by the end it states
        o, d = outcome(*mcp(1, 1), "invariants")
        self.assertEqual(o, "fail")
        self.assertIn("a branch end after 1 branch(es) of max_label_branches 2", d)
        # the trie preset mirrors an exhaustive request: no branch end can have happened on either arm
        sc_trie, a_trie = self.structural_walk(strategy={"exhaustive": True}, echo={"exhaustive": True})
        sc_trie["mcp"] = {"arguments": {"preset": "trie", "overrides": {"bounds": {"max_extension_bp": 100}}}}
        self.assertEqual(outcome(sc_trie, a_trie, "invariants", _arm="mcp")[0], "fail")
        self.assertEqual(outcome(sc_trie, a_trie, "invariants")[0], "fail")

    def test_structural_ends_index_each_arm_once(self):
        # the route split table is built once per arm, not per path: a wide trie (max_paths 10,000 leaves) would
        # otherwise cost time quadratic in its paths
        calls = []
        real = checks._route_splits

        def counting(armo):
            calls.append(id(armo))
            return real(armo)
        sc, ans = self.structural_walk(strategy={"branching": {"max_label_branches": 2}})
        checks._route_splits = counting
        try:
            self.assertEqual(outcome(sc, ans, "invariants")[0], "pass")
            self.assertEqual(len(calls), 1)          # one arm, three paths
            del calls[:]
            self.assertEqual(outcome(sc, ans, "walk_in_record")[0], "pass")
            self.assertEqual(len(calls), 1)
        finally:
            checks._route_splits = real

    def test_zero_extension_branch_stop(self):
        # the zero-extension mutation of T-derived-0001: both arms a 0 bp root whose every label ends
        # `branch` with branches 0 under max_label_branches 2; the dead_end control fails as well
        import random
        rnd = random.Random(5)
        record = "".join(rnd.choice("ACGT") for _ in range(3000))
        start, L = 1000, 100
        t = {"record": "x", "taxid": "1", "seq_id": 0, "record_length": 3000, "start": start, "length": L,
             "strand": "+", "window": record[start:start + L], "flank": record, "flank_start": 0}
        req = {"seeds": [{"sequence": t["window"]}],
               "strategy": {"direction": "both", "bounds": {"max_extension_bp": 1000},
                            "branching": {"max_label_branches": 2}}}
        sc = {"route": "POST /traverse", "request": req, "truth": t}

        def zero(reason, branches=0):
            arms = {}
            for side in ("left", "right"):
                labels = list(range(8))
                seg = {"id": 0, "parents": [], "from_bp": 0, "length_bp": 0, "sequence": "", "labels": labels,
                       "events": [{"at_bp": 0, "type": "label_end", "label": x, "reason": reason,
                                   "structural_successors": 2} for x in labels]}
                arms[side] = {"status": "complete", "complete_to_bp": 1000, "segments": [seg], "splits": [],
                              "paths": [{"id": 0, "segments": [0], "length_bp": 0, "end_reasons": {reason: 8},
                                         "end_labels": [{"label": x, "branches": branches, "loss": 0.0}
                                                        for x in labels]}]}
            return ok_answer({"results": [{"outcome": {"walks": "partial"}, "arms": arms}],
                              "strategy": {"bounds": {"max_extension_bp": 1000},
                                           "branching": {"max_label_branches": 2}}})
        for ans in (zero("branch"), zero("branch", branches=2), zero("dead_end")):
            self.assertEqual(outcome(sc, ans, "walk_in_record")[0], "fail")
        o, d = outcome(sc, zero("branch"), "invariants")
        self.assertEqual(o, "fail")
        self.assertIn("a branch end after 0 branch(es) of max_label_branches 2", d)
        o, d = outcome(sc, zero("branch", branches=2), "invariants")
        self.assertEqual(o, "fail")
        self.assertIn("0 ambiguous split(s)", d)
        self.assertEqual(outcome(sc, zero("dead_end"), "invariants")[0], "pass")

    def test_structural_ends_in_a_graphlet(self):
        # the graphlet's R records state the lineage's structural successors and branches (SPEC §7.5.2), its G
        # records the splits: the same conditions as the full result
        import random
        rnd = random.Random(13)
        record = "".join(rnd.choice("ACGT") for _ in range(1000))
        start, L = 400, 40
        t = {"record": "x", "taxid": "1", "seq_id": 0, "record_length": 1000, "start": start, "length": L,
             "strand": "+", "window": record[start:start + L], "flank": record, "flank_start": 0}
        req = {"seeds": [{"sequence": t["window"]}],
               "strategy": {"direction": "right", "bounds": {"max_extension_bp": 100},
                            "branching": {"max_label_branches": 2}}}
        sc = {"route": "POST /traverse", "request": req, "truth": t}
        cont = record[start + L:start + L + 100]
        off = lambda b: "A" if b != "A" else "C"  # noqa: E731

        def mgt(branches=2, split="1", successors="2"):
            return "\n".join([
                "H mgt 1 31 basic $ACGT c k m 64 1000 0 walk * * *",
                "A r c 100 u 0 0 1 1 0 steps=31 * 0 * 5 3 3 2 0 32",
                "G * 0 10 * * * * %s * %s" % (split, cont[:10]),
                "G 0 10 10 * * * * %s * %s" % (split, cont[10:20]),
                "G 0 10 1 * * * * * * %s" % off(cont[10]), "T *",
                "G 1 20 10 * * * * * * %s" % cont[20:30], "T *",
                "G 1 20 1 * * * * * * %s" % off(cont[20]), "T *",
                "R 2 0 0 11 D 0 * * 0 1 0",
                "R 3 0 0 30 R 0 * * %s %d 0" % (successors, branches),
                "R 4 0 0 21 D 0 * * 0 2 0"]) + "\n"

        def ans(**kw):
            return ok_answer({"results": [{"outcome": {"walks": "complete"}, "graphlet": mgt(**kw)}],
                              "strategy": {"bounds": {"max_extension_bp": 100},
                                           "branching": {"max_label_branches": 2}}})
        arms = checks.graphlet_arms(mgt())
        leaf = arms["right"]["paths"][1]
        self.assertEqual((leaf["length_bp"], leaf["end_reasons"], leaf["end_labels"]),
                         (30, {"branch": 1}, [{"label": 0, "branches": 2}]))
        self.assertEqual(arms["right"]["segments"][3]["events"],
                         [{"at_bp": 30, "type": "label_end", "label": 0, "reason": "branch",
                           "structural_successors": 2}])
        for arm in ("mcp", "direct-local"):
            self.assertEqual(outcome(sc, ans(), "walk_in_record", _arm=arm)[0], "pass")
            self.assertEqual(outcome(sc, ans(), "invariants", _arm=arm)[0], "pass")
            for bad in ({"branches": 0}, {"split": "0"}, {"successors": "1"}):
                self.assertEqual(outcome(sc, ans(**bad), "walk_in_record", _arm=arm)[0], "fail", bad)
                self.assertEqual(outcome(sc, ans(**bad), "invariants", _arm=arm)[0], "fail", bad)

    def test_every_stored_traverse_answer_states_possible_structural_ends(self):
        # the engine's stored answers (CLI documents of tiny indexes, full and graphlet, and the mini index's
        # cells) state only structural ends that can have happened
        n = ends = 0
        dirs = [os.path.join(DATA, x) for x in ("documents", os.path.join("documents", "budgets"), "coords",
                                                 "review3", "r21")] + [self.MINI]
        for d in dirs:
            for name in sorted(os.listdir(d)):
                if not (name.endswith(".full.json") or name.endswith(".graphlet.json") or
                        name.endswith(".full.json.gz") or name.endswith(".graphlet.json.gz")):
                    continue
                stem = name.split(".")[0]
                opener = (lambda f: gzip.open(f, "rt")) if name.endswith(".gz") else open
                with opener(os.path.join(d, name)) as f:
                    doc = json.load(f)
                req = {}
                for rq in (stem + ".request.json", stem + ".request.json.gz"):
                    if os.path.exists(os.path.join(d, rq)):
                        with ((lambda f: gzip.open(f, "rt")) if rq.endswith(".gz") else open)(os.path.join(d, rq)) as f:
                            req = json.load(f)
                sc = {"route": "POST /traverse", "request": req, "truth": {}}
                ans = ok_answer(doc)
                for i in range(len(doc.get("results") or [])):
                    r = checks._result(ans, i)
                    for side, armo in (r.get("arms") or {}).items():
                        for p in armo.get("paths") or []:
                            self.assertEqual(checks._structural_end_problems(sc, ans, armo, p), [], name)
                            ends += len(set(p.get("end_reasons") or {}) & checks._STATED_STRUCTURAL_ENDS)
                            n += 1
        self.assertGreater(n, 500)
        self.assertGreater(ends, 100)

    def test_walk_blocked_by_its_seed(self):
        # G28 (VERIFY-2 item 2, the strict rule): the record's own continuation re-enters the seed at 0 bp, the walk
        # states that block and continues from the repeat's exit: it follows the record
        self.assertEqual(checks._first_block("A" * 5 + "C" * 40, "C" * 3, "right", False), (0, "rejoined_seed"))
        import random
        rnd = random.Random(9)
        seed = "".join(rnd.choice("ACGT") for _ in range(31))
        x = "".join(rnd.choice("ACGT") for _ in range(60))
        self.assertEqual(checks._first_block(seed, x + x, "right", False), (60 + 31, "edge_reuse"))
        self.assertIsNone(checks._first_block(seed, x, "right", False))
        sc, ans, record = self.repeat_walk()
        o, d = outcome(sc, ans, "walk_in_record", label_scope="record")
        self.assertEqual(o, "pass", d)
        self.assertIn("followed for 0 bp, then blocked by rejoined_seed", d)
        self.assertEqual(outcome(sc, mcp_shaped(ans), "walk_in_record", _arm="mcp", label_scope="record")[0], "pass")
        seg = lambda a: a["body"]["results"][0]["arms"]["right"]["segments"][0]  # noqa: E731
        # the stated reason must be the computed one, the char the record's base, the block stated at all
        for change in ({"reason": "edge_reuse"}, {"char": "T" if record[346] != "T" else "A"}, {"at_bp": 1}):
            bad = copy.deepcopy(ans)
            seg(bad)["events"][0].update(change)
            self.assertEqual(outcome(sc, bad, "walk_in_record", label_scope="record")[0], "fail", change)
        bad = copy.deepcopy(ans)
        seg(bad)["events"] = seg(bad)["events"][1:]
        self.assertEqual(outcome(sc, bad, "walk_in_record", label_scope="record")[0], "fail")
        # the continuation past the block must spell the record from the repeat's exit
        bad = copy.deepcopy(ans)
        s = seg(bad)["sequence"]
        seg(bad)["sequence"] = s[:50] + ("A" if s[50] != "A" else "C") + s[51:]
        self.assertEqual(outcome(sc, bad, "walk_in_record", label_scope="record")[0], "fail")
        # a block before the truth's minimum is a truth gap: inconclusive
        gap = dict(sc, truth=dict(sc["truth"], walk={"right": {"min_bp": 10}}))
        o, d = outcome(gap, ans, "walk_in_record", label_scope="record")
        self.assertEqual(o, "inconclusive", d)
        self.assertIn("truth gap", d)
        # past the truth's flank: not verifiable (inconclusive); past the record's end with the record as the only
        # label: fail
        sc_short, ans_short, _ = self.repeat_walk(flank_end=420)
        o, d = outcome(sc_short, ans_short, "walk_in_record", label_scope="record")
        self.assertEqual(o, "inconclusive", d)
        self.assertIn("not verifiable", d)
        sc_end = dict(sc_short, truth=dict(sc_short["truth"], record_length=420))
        o, d = outcome(sc_end, ans_short, "walk_in_record", label_scope="record")
        self.assertEqual(o, "fail", d)
        self.assertIn("past the record's end", d)

    def test_walk_blocked_exit_outside_the_flank(self):
        # A2-04 (review of lane S5b, probe_w3): the truth's flank ends at 370, the repeat's exit base is 376: the
        # computed block (0 bp, rejoined_seed) has no exit in the flank; the walk that states it and continues from
        # the exit cannot be verified past the block: inconclusive (it failed: "no route ... follows the record")
        sc, ans, record = self.repeat_walk(flank_end=370)
        t = sc["truth"]
        blk = checks._computed_block(sc, ans, ans["body"]["results"][0], 0, t,
                                     {"start": t["start"], "flank": t["flank"], "flank_start": 0}, "right", 40)
        self.assertEqual((blk["d"], blk["reason"], blk["exits"]), (0, "rejoined_seed", []))
        o, d = outcome(sc, ans, "walk_in_record", label_scope="record")
        self.assertEqual(o, "inconclusive", d)
        self.assertIn("which is not in the truth's flank: not verifiable", d)
        # the block must still be stated as computed: another reason or char, another depth, or none fails
        seg = lambda a: a["body"]["results"][0]["arms"]["right"]["segments"][0]  # noqa: E731
        for change in ({"reason": "edge_reuse"}, {"char": "T" if record[346] != "T" else "A"}, {"at_bp": 1}):
            bad = copy.deepcopy(ans)
            seg(bad)["events"][0].update(change)
            self.assertEqual(outcome(sc, bad, "walk_in_record", label_scope="record")[0], "fail", change)
        bad = copy.deepcopy(ans)
        seg(bad)["events"] = seg(bad)["events"][1:]
        self.assertEqual(outcome(sc, bad, "walk_in_record", label_scope="record")[0], "fail")
        # a block before the truth's minimum: a truth gap (inconclusive)
        gap = dict(sc, truth=dict(sc["truth"], walk={"right": {"min_bp": 10}}))
        o, d = outcome(gap, ans, "walk_in_record", label_scope="record")
        self.assertEqual(o, "inconclusive", d)
        self.assertIn("truth gap", d)
        # the exit inside the flank: verified as before (pass on the full record)
        self.assertEqual(outcome(self.repeat_walk()[0], ans, "walk_in_record", label_scope="record")[0], "pass")

    def test_duplicate_seed_count(self):
        # B5b-07 (review of lane S5b): the invariants check states the distinct seeds of a request with identical
        # seeds (the service's identity: the cleaned sequence and the trimmed, de-duplicated labels), so that the
        # table explains exactly a merge (N results == N distinct), not a lost or a surplus result
        seq = "ACGT" * 10
        def sc(seeds):
            return {"route": "POST /traverse", "truth": {},
                    "request": {"seeds": seeds, "strategy": {"direction": "right"}}}
        def ans(n):
            return ok_answer({"results": [{"outcome": {"walks": "not_walked"}} for _ in range(n)]})
        one = {"sequence": seq, "labels": ["NC_1"]}
        cases = [([one, dict(one)], 1, "1 results for 2 seeds (1 distinct)"),
                 ([one, {"sequence": seq.lower() + "\n", "labels": [" NC_1 ", "NC_1"]}], 1,
                  "1 results for 2 seeds (1 distinct)"),
                 ([one, one, {"sequence": seq, "labels": ["NC_2"]}], 1, "1 results for 3 seeds (2 distinct)"),
                 ([one, dict(one)], 5, "5 results for 2 seeds (1 distinct)"),
                 ([one, {"sequence": seq, "labels": ["NC_2"]}], 1, "1 results for 2 seeds"),
                 ([one, dict(one, seed_id="b")], 1, "1 results for 2 seeds (1 distinct)")]
        tightened = r"([1-9][0-9]*) results for [0-9]+ seeds \(\1 distinct\)"
        import re
        for seeds, n, detail in cases:
            o, d = outcome(sc(seeds), ans(n), "invariants")
            self.assertEqual(o, "fail")
            self.assertEqual(d, detail)
            self.assertEqual(re.fullmatch(tightened, d, re.S) is not None, detail.startswith("1 results for 2 seeds ("),
                             detail)

    def test_paths_end(self):
        sc = {"route": "POST /traverse", "truth": {},
              "request": {"seeds": [{"sequence": "A" * 40}], "strategy": {"direction": "right"}}}
        def ans(reason, length):
            return ok_answer({"results": [{"arms": {"right": {"paths": [
                {"id": 0, "length_bp": length, "path_reason": "branch", "end_reasons": {reason: 1}}]}}}]})
        self.assertEqual(outcome(sc, ans("below_min_labels", 0), "paths_end", reasons_in=["below_min_labels"], max_bp=1)[0], "pass")
        self.assertEqual(outcome(sc, ans("max_extension_bp", 200), "paths_end", reasons_in=["below_min_labels"])[0], "fail")
        self.assertEqual(outcome(sc, ans("below_min_labels", 50), "paths_end", reasons_in=["below_min_labels"], max_bp=1)[0], "fail")
        empty = ok_answer({"results": [{"arms": {"right": {"paths": []}}}]})
        self.assertEqual(outcome(sc, empty, "paths_end", reasons_in=["below_min_labels"])[0], "fail")
        # G16: the deployed engine ends a quorum stop as `branch`, the text below_min_labels in the leaf's label_end
        # event (SPEC §13 I3 note 2): T-quorum-0001's answer
        def quorum(at_bp=0):
            seg = {"id": 0, "parents": [], "from_bp": 0, "length_bp": 0, "sequence": "",
                   "events": [{"at_bp": at_bp, "label": 0, "reason": "below_min_labels", "structural_successors": 1,
                               "type": "label_end"}]}
            return ok_answer({"results": [{"outcome": {"walks": "complete"}, "arms": {"right": {
                "status": "complete", "segments": [seg],
                "paths": [{"id": 0, "segments": [0], "length_bp": 0, "end_reasons": {"branch": 1}}]}}}]})
        self.assertEqual(outcome(sc, quorum(), "paths_end", reasons_in=["below_min_labels"], max_bp=1)[0], "pass")
        self.assertEqual(outcome(sc, quorum(), "paths_end", reasons_in=["minority"], max_bp=1)[0], "fail")
        self.assertEqual(outcome(sc, quorum(at_bp=5), "paths_end", reasons_in=["below_min_labels"])[0], "fail")
        # the graphlet states it as the R record's qualifier Rb (SPEC §7.5.2), read through _END_QUALIFIERS
        mgt = "\n".join(["A r c 0 u 0 0 1 0 0 steps=0 * 0 * 1 0 1 0 0 0", "G * 0 0 * * * * * * .", "T *",
                         "R 0 0 0 0 Rb 0 * * 1 0 0"])
        g = ok_answer({"results": [{"outcome": {"walks": "complete"}, "graphlet": mgt}]})
        self.assertEqual(outcome(sc, g, "paths_end", reasons_in=["below_min_labels"], max_bp=1, _arm="mcp")[0], "pass")
        arms = checks.graphlet_arms(mgt)
        self.assertEqual(arms["right"]["paths"][0]["end_reasons"], {"branch": 1})
        self.assertEqual(arms["right"]["segments"][0]["events"],
                         [{"at_bp": 0, "type": "label_end", "label": 0, "reason": "below_min_labels",
                           "structural_successors": 1}])
        self.assertEqual(arms["right"]["paths"][0]["end_labels"], [{"label": 0, "branches": 0}])

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
        # direction both: an arm whose walks all ended before the cap stays complete (SPEC-labeled-traversal-core
        # §6.7), for a per-arm cap and a per-seed stop alike; some arm must show the cap
        both = dict(sc, request={"seeds": [{"sequence": "A" * 40}], "strategy": {"direction": "both", "bounds": {
            "max_output_bp": 50, "max_steps": 5}}})

        def arm(status, trig, knob):
            return {"status": status, "complete_to_bp": 5, "cap_trigger": {"reason": trig} if trig else None,
                    "limitations": [{"kind": "walk_domain", "knob": knob}] if trig else []}

        def two(left, right):
            return ok_answer({"results": [{"outcome": {"walks": "partial"}, "limitations": [], "resource_stop": None,
                                           "arms": {"left": left, "right": right}}]})
        for cap, knob in (("max_output_bp", "bounds.max_output_bp"), ("max_steps", "bounds.max_steps")):
            args = dict(knob=knob, cap_reason=cap)
            self.assertEqual(outcome(both, two(arm("complete", None, knob), arm("truncated", cap, knob)),
                                     "traverse_stop", **args)[0], "pass", cap)
            self.assertEqual(outcome(both, two(arm("truncated", cap, knob), arm("truncated", cap, knob)),
                                     "traverse_stop", **args)[0], "pass", cap)
            o, d = outcome(both, two(arm("complete", None, knob), arm("complete", None, knob)), "traverse_stop", **args)
            self.assertEqual(o, "fail", cap)
            self.assertIn("no requested arm truncated by %s" % cap, d)
            self.assertEqual(outcome(both, two(arm("truncated", "max_paths", knob), arm("truncated", cap, knob)),
                                     "traverse_stop", **args)[0], "fail", cap)


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

    def test_walk_end_resource_end_below_the_frontier(self):
        # a resource end of a path shorter than its arm's complete_to_bp cannot have happened (a cap ends walks at or
        # past it, SPEC-labeled-traversal-core §6.10): the record's end not reached is no stated stop then
        import random
        rnd = random.Random(5)
        record = "".join(rnd.choice("ACGT") for _ in range(400))
        start, L = 100, 40
        t = {"record": "x", "taxid": "1", "seq_id": 0, "record_length": 400, "start": start, "length": L,
             "strand": "+", "window": record[start:start + L], "flank": record, "flank_start": 0}
        sc = {"route": "POST /traverse", "truth": t, "request": {"seeds": [{"sequence": t["window"]}], "strategy": {
            "direction": "right", "bounds": {"max_extension_bp": 1000, "time_budget_ms": 15000}}}}
        right = record[start + L:start + L + 30]
        capped = dict(status="truncated", trig="time_budget", reason="time_budget")
        o, d = outcome(sc, self.walk({"right": right}, **capped), "walk_end", arm="right")
        self.assertEqual(o, "inconclusive")
        self.assertIn("stopped by time_budget", d)
        o, d = outcome(sc, self.walk({"right": right}, extra={"complete_to_bp": 60}, **capped), "walk_end", arm="right")
        self.assertEqual(o, "fail")
        self.assertIn("below the arm's complete_to_bp 60", d)

    def test_blocked_event(self):
        # N-hard-0022: a tandem-repeat seed re-enters itself at 0 bp on each arm (SPEC §6.6)
        sc = {"route": "POST /traverse", "request": {"seeds": [{"sequence": "TTAGGG" * 10}],
                                                     "strategy": {"direction": "both"}}, "truth": {}}
        def walk(events_left, events_right):
            arms = {}
            for side, evs in (("left", events_left), ("right", events_right)):
                arms[side] = {"status": "complete", "segments": [{"id": 0, "parents": [], "from_bp": 0, "length_bp": 0,
                                                                  "sequence": "", "events": evs}], "paths": []}
            return ok_answer({"results": [{"outcome": {"walks": "complete"}, "arms": arms}]})
        rej = {"at_bp": 0, "type": "blocked", "char": "T", "reason": "rejoined_seed", "labels": [0]}
        self.assertEqual(outcome(sc, walk([rej], [rej]), "blocked_event", reason="rejoined_seed", max_at_bp=0)[0], "pass")
        o, d = outcome(sc, walk([rej], []), "blocked_event", reason="rejoined_seed", max_at_bp=0)
        self.assertEqual(o, "fail")
        self.assertIn("right: no blocked successor rejoined_seed", d)
        self.assertEqual(outcome(sc, walk([rej], [dict(rej, reason="edge_reuse")]), "blocked_event",
                                 reason="rejoined_seed")[0], "fail")
        self.assertEqual(outcome(sc, walk([rej], [dict(rej, at_bp=7)]), "blocked_event", reason="rejoined_seed",
                                 max_at_bp=0)[0], "fail")
        # the graphlet (the MCP's export) states it in its E records
        mgt = "\n".join(["A l c 0 p 0 0 1 0 0 steps=0 * 0 * 1 0 1 0 0 0", "G * 0 0 * * * * * * .", "E 0 b A J 1 0",
                         "T J .", "A r c 0 p 0 0 1 0 0 steps=0 * 0 * 1 0 1 0 0 0", "G * 0 0 * * * * * * .",
                         "E 0 b T J 1 0", "T J ."])
        g = ok_answer({"results": [{"outcome": {"walks": "complete"}, "graphlet": mgt}]})
        self.assertEqual(outcome(sc, g, "blocked_event", reason="rejoined_seed", max_at_bp=0, _arm="mcp")[0], "pass")
        arms = checks.graphlet_arms(mgt)
        self.assertEqual(arms["right"]["segments"][0]["events"],
                         [{"at_bp": 0, "type": "blocked", "char": "T", "reason": "rejoined_seed"}])

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
        n = label_ends = label_sets = 0
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
                    # the label_end events derived from the R records (SPEC §7.5.6) are the full result's
                    # the lineage's branches at each leaf (end_labels) are the full result's
                    self.assertEqual([[(x["label"], x["branches"]) for x in p.get("end_labels") or []]
                                      for p in fa["paths"]],
                                     [[(x["label"], x["branches"]) for x in p["end_labels"]]
                                      for p in arms[side]["paths"]], cell)
                    for sf, sg in zip(fa["segments"], arms[side]["segments"]):
                        ends = lambda x: sorted((e["at_bp"], e.get("label"), e["reason"],  # noqa: E731
                                                 e.get("structural_successors"))
                                                for e in x.get("events") or [] if e.get("type") == "label_end")
                        self.assertEqual(ends(sf), ends(sg), cell)
                        label_ends += len(ends(sf))
                        # A-04: the entry sets derived from the G records (SPEC §7.5.3) are the segments' labels
                        self.assertEqual((sf.get("labels"), bool(sf.get("labels_truncated"))),
                                         (sg["labels"][:len(sf.get("labels") or [])] if sf.get("labels_truncated")
                                          else sg["labels"], sg["labels_truncated"]), cell)
                        label_sets += 1
                self.assertEqual([x["name"] for x in checks.graphlet_label_dict(rg["graphlet"])],
                                 [x["name"] for x in rf.get("label_dict") or []], cell)
        self.assertGreater(n, 4)
        self.assertGreater(label_ends, 50)
        self.assertGreater(label_sets, 100)
        # a body whose sets the rules cannot read leaves them out (the route search falls back to the sequence)
        arms = checks.graphlet_arms("H mgt 1 31 basic $ACGT c k m 64 1000 0 walk * * *\nA r c 10 u 0 0 1 0 0 x *\n"
                                    "G * 0 10 !3 * * * * * ACGTACGTAC\nT X\n")
        self.assertNotIn("labels", arms["right"]["segments"][0])


if __name__ == "__main__":
    unittest.main()
