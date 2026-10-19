"""The five-way compare of area G (compare_graphlets.py; DESIGN-v8 §5, FORMAT-v8 §3-§4): the library's per-tool answer
keys frozen in the module, held to the library's answers on the mini_refseq fixtures (the tests may import the
library; the module does not); and on synthetic runs every op outcome: agree, disagree naming the pairs (W1-W2
transport, W1-W4 fetch path, W3-W4 service wrapper, W-W5 library vs independent reading), not_exposed for W3 with
the owner's reason, precondition_failed, inconclusive; args_resolved first; handles to symbols, cursors removed,
the service fields stripped; errors by code and structured fields; the cut-of-uncut rule; the drift lane."""

import copy
import gzip
import json
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import compare_graphlets as cg  # noqa: E402
import graphlet_fixtures as gf  # noqa: E402
import graphlet_precheck as gp  # noqa: E402

SID = "G-merge-0001"
UUID = "0d5b9b1a-1111-2222-3333-444455556666"
UUID_B = "9e5b9b1a-aaaa-bbbb-cccc-444455556666"


def library_tools(local_limits=False):
    """GraphletTools over a fresh store holding the five mini_refseq graphlets: (tools, [handle], tmp dir)."""
    sys.path.insert(0, gf.LIB)
    from metagraph.traverse import mcp_tools
    from metagraph.traverse.store import GraphletStore
    tmp = tempfile.mkdtemp(prefix="lane_c_lib_")
    os.makedirs(os.path.join(tmp, "exp"))
    store = GraphletStore(os.path.join(tmp, "spool"))
    tools = mcp_tools.GraphletTools(store, max_bytes=2048, export_dir=os.path.join(tmp, "exp"),
                                    local_limits=mcp_tools.ToolLimits() if local_limits else None)
    handles = []
    for name in sorted(os.listdir(gf.MINI)):
        if name.endswith(".graphlet.json.gz"):
            cell = name[:-len(".graphlet.json.gz")]
            resp, req = gf.load_fixture(cell, "graphlet"), gf.load_fixture(cell, "request")
            handles.extend(store.put(resp, req, i) for i, r in enumerate(resp["results"]) if "graphlet" in r)
    return tools, handles, tmp, mcp_tools


class Whitelist(unittest.TestCase):
    """ANSWER_KEYS against the library's answers: every key a tool answers is in its row of the table."""

    def answers(self, local_limits):
        tools, handles, tmp, mod = library_tools(local_limits)
        big = 16777216
        out = []

        def call(tool, *a, **kw):
            ans = getattr(tools, tool)(*a, **kw)
            out.append((tool, ans))
            return ans
        for h in handles:
            for det in ("summary", "limitations"):
                for mb in (big, 2048, 512):
                    call("graphlet_summary", h, detail=det, max_bytes=mb)
            for arm in ("left", "right"):
                walks = call("graphlet_walks", h, arm=arm, n=1000, spell="full", max_bytes=big)
                for rank in ("support", "length", "loss", "id"):
                    call("graphlet_walks", h, arm=arm, rank=rank, n=2, spell="tail", max_bytes=2048)
                labels = call("graphlet_labels", h, arm=arm, n=1000, max_bytes=big)
                if labels.get("rows"):
                    ref = labels["rows"][0]["ref"]
                    for rc in (True, False):
                        call("graphlet_walks", h, arm=arm, label={"ref": ref}, route_consistent=rc, max_bytes=big)
                    call("graphlet_labels", h, arm=arm, name={"ref": ref}, max_bytes=big)
                    call("graphlet_claims", h, arm=arm, labels={"ref": ref}, max_bytes=big)
                call("graphlet_claims", h, arm=arm, route_consistent=False, max_bytes=big)
                call("graphlet_splits", h, arm=arm, min_labels_before=0, n=1000, max_bytes=big)
                for row in (walks.get("rows") or [])[:3]:
                    w = row["walk"]
                    call("graphlet_walk", h, arm=arm, walk=w, max_bytes=big)
                    call("graphlet_walk", h, arm=arm, walk=w, max_bytes=300)
                    call("graphlet_support", h, arm=arm, walk=w, max_bytes=big)
                    call("graphlet_sequence", h, arm=arm, walk=w, orientation="walk", with_seed=True)
                    call("graphlet_sequence", h, arm=arm, segment=0)
                    call("traverse_continue", h, arm=arm, walk=w, execute=False)
            for fmt in ("fasta", "gfa", "json", "mgt"):
                call("graphlet_export", h, format=fmt, path="x_%s.%s" % (h, fmt))
            call("graphlet_export", h, format="json", what="walk", arm="right", walks=[0], path="w_%s.json" % h)
            lab = tools.graphlet_labels(h, max_bytes=big)
            if lab.get("rows"):
                sub = call("graphlet_subtrie", h, labels=[{"ref": lab["rows"][0]["ref"]}])
                if "handle" in sub:
                    call("graphlet_summary", sub["handle"], max_bytes=big)
            call("graphlet_save", h, path="s_%s.mgt" % h)
            call("graphlet_load", path="s_%s.mgt" % h)
        for h1 in handles:
            for h2 in handles:
                for mode in ("claims", "walks", "labels", "prefix_subset"):
                    call("graphlet_compare", h1, h2, mode=mode, max_bytes=big)
                    call("graphlet_compare", h1, h2, mode=mode, max_bytes=700)
        if local_limits:
            for h in handles[1:]:
                for wu in (1, 50, 5000):
                    call("graphlet_walks", h, arm="right", n=1000, spell="full", budget={"work_units": wu},
                         max_bytes=big)
                    call("graphlet_claims", h, arm="right", budget={"work_units": wu}, max_bytes=big)
                    call("graphlet_splits", h, arm="right", budget={"work_units": wu}, max_bytes=big)
        call("graphlet_list", max_bytes=big)
        call("graphlet_free", handles[0])
        return out, mod

    def test_every_key_the_library_answers_is_in_the_table(self):
        seen = {}
        for local in (False, True):
            answers, mod = self.answers(local)
            for tool, ans in answers:
                if "error" in ans:
                    continue
                self.assertIn(tool, cg.ANSWER_KEYS)
                self.assertEqual(sorted(set(ans) - cg.ANSWER_KEYS[tool]), [], (tool, local))
                seen.setdefault(tool, set()).update(ans)
        self.assertGreater(len(seen), 12)
        # the table has a row for every graphlet tool and traverse_continue (the backend tools answer the server)
        names = [n for n in mod.tool_names() if n.startswith("graphlet_")] + ["traverse_continue"]
        self.assertEqual(sorted(set(names) - set(cg.ANSWER_KEYS)), [])
        # the answers exercised most of each row (a key never seen is an optional one: a parent, a stop, ...)
        for tool in ("graphlet_walks", "graphlet_compare", "graphlet_summary", "graphlet_sequence"):
            self.assertGreaterEqual(len(seen[tool]) / float(len(cg.ANSWER_KEYS[tool])), 0.5, tool)

    def test_the_service_fields_are_outside_the_table(self):
        for tool, keys in (("graphlet_walks", ("walk_display",)), ("graphlet_walk", ("walk_display",)),
                           ("graphlet_compare", ("label_evidence", "qualified_cause", "hint")),
                           ("graphlet_summary", ("in_service", "in_service_reason", "download", "delivery")),
                           ("graphlet_export", ("export_id", "status", "format", "handle", "link", "gz_bytes",
                                                "in_service"))):
            for k in keys:
                self.assertNotIn(k, cg.ANSWER_KEYS[tool], (tool, k))
        self.assertIn("hint", cg.ANSWER_KEYS["graphlet_walks"])   # the library's own filter hint
        lib, svc = cg.split_service_fields("graphlet_walks", {"handle": "h", "rows": [], "walk_display": "first_parent"})
        self.assertEqual((sorted(lib), svc), (["handle", "rows"], {"walk_display": "first_parent"}))


class Canonical(unittest.TestCase):
    def test_handles_cursors_and_local(self):
        a = {"handle": "g_0123456789ab", "next_cursor": "x.y", "local": {"usage": 1},
             "evidence": {"view": {"of": "g_ffffffffffff"}, "limitations": {"local": ["local_work"], "seed": []}},
             "rows": [{"note": "see g_0123456789ab"}]}
        b = {"handle": UUID, "evidence": {"view": {"of": UUID_B}, "limitations": {"seed": []}},
             "rows": [{"note": "see %s" % UUID}]}
        self.assertEqual(cg.canonical_answer("graphlet_walks", a), cg.canonical_answer("graphlet_walks", b))
        kept = cg.canonical_answer("graphlet_walks", a, keep_local=True)
        self.assertIn("local", kept)

    def test_compare_handles_by_role_and_rows_a_b_are_values(self):
        lib = {"a": "g_0000000000aa", "b": "g_0000000000bb", "rows": [{"side": "differ", "a": "x", "b": "y"}]}
        svc = {"a": UUID, "b": UUID, "rows": [{"side": "differ", "a": "x", "b": "y"}]}   # one handle twice
        self.assertEqual(cg.canonical_answer("graphlet_compare", lib), cg.canonical_answer("graphlet_compare", svc))
        swapped = copy.deepcopy(svc)
        swapped["rows"][0].update(a="y", b="x")
        self.assertIsNotNone(gp.first_difference(cg.canonical_answer("graphlet_compare", lib),
                                                 cg.canonical_answer("graphlet_compare", swapped)))

    def test_errors_by_code_and_structured_fields(self):
        e1 = {"code": "result_too_large", "fields": {"bytes": 3000, "max_bytes": 2048, "message": "a"}}
        e2 = {"code": "result_too_large", "fields": {"bytes": 3022, "max_bytes": 2048, "message": "b"}}
        self.assertEqual(cg.canonical_error(e1), cg.canonical_error(e2))
        e3 = {"code": "unknown_handle", "fields": {"replayable": True, "hint": "replay"}}
        e4 = {"code": "unknown_handle", "fields": {"replayable": False, "reason": "job_expired"}}
        self.assertIsNotNone(gp.first_difference(cg.canonical_error(e3), cg.canonical_error(e4)))

    def test_export_path_and_bytes_are_left_out(self):
        a = {"path": "/tmp/x/g_0123456789ab.fa", "bytes": 100, "records": 4}
        b = {"path": "/srv/y.fa", "bytes": 140, "records": 4}
        self.assertEqual(cg.canonical_answer("graphlet_export", a), cg.canonical_answer("graphlet_export", b))


class Cut(unittest.TestCase):
    UNCUT = {"handle": "<handle>", "summary": {"caveats": [{"kind": "k", "effect": "long text"}, {"kind": "j"}],
                                               "top_labels": [1, 2, 3], "seed": {"seed_id": "s"}}}

    def test_a_summary_cut(self):
        cut = {"handle": "<handle>", "summary": {"caveats": [{"kind": "k"}], "caveat_effects_cut": True,
                                                 "caveats_truncated": True, "top_labels": [1], "seed": {"seed_id": "s"}}}
        self.assertIsNone(cg.is_cut(cut, self.UNCUT, cg.OPTIONAL_KEYS["graphlet_summary"]))
        self.assertIsNone(cg.is_cut(self.UNCUT, self.UNCUT))

    def test_what_is_no_cut(self):
        opt = cg.OPTIONAL_KEYS["graphlet_summary"]
        lacks = copy.deepcopy(self.UNCUT)
        del lacks["summary"]["seed"]
        self.assertEqual(cg.is_cut(lacks, self.UNCUT, opt)[0], "/summary/seed")
        extra = copy.deepcopy(self.UNCUT)
        extra["summary"]["new"] = 1
        self.assertIn("not a cut marker", cg.is_cut(extra, self.UNCUT, opt)[1])
        longer = copy.deepcopy(self.UNCUT)
        longer["summary"]["top_labels"].append(4)
        self.assertIn("longer", cg.is_cut(longer, self.UNCUT, opt)[1])
        other = copy.deepcopy(self.UNCUT)
        other["summary"]["caveats"][0]["effect"] = "long tex!"
        self.assertIn("not a prefix", cg.is_cut(other, self.UNCUT, opt)[1])
        scalar = copy.deepcopy(self.UNCUT)
        scalar["summary"]["top_labels"][0] = 9
        self.assertEqual(cg.is_cut(scalar, self.UNCUT, opt)[0], "/summary/top_labels/0")

    def test_a_row_marked_cut(self):
        uncut = {"total": 2, "rows": [{"walk": 1, "sequence": "ACGTACGT", "labels_full": [1, 2, 3]}, {"walk": 0}]}
        cut = {"total": 2, "row_truncated": True, "cut_fields": ["sequence", "labels_full"],
               "rows": [{"walk": 1, "sequence": "ACG", "labels_full": [1], "_cut_fields": ["sequence", "labels_full"]},
                        {"walk": 0}]}
        self.assertIsNone(cg.is_cut(cut, uncut))
        self.assertIsNone(cg.cut_aware_difference(cut, uncut))
        self.assertIsNone(cg.cut_aware_difference(uncut, cut))
        bad = copy.deepcopy(cut)
        bad["rows"][0]["sequence"] = "ACT"
        d = cg.cut_aware_difference(uncut, bad)
        self.assertEqual((d[0], d[1], d[2]), ("/rows/0/sequence", "ACGTACGT", "ACT"))
        self.assertIn("no cut of the other's", d[3])
        scalar = copy.deepcopy(cut)
        scalar["rows"][0]["walk"] = 2
        self.assertEqual(cg.cut_aware_difference(uncut, scalar)[0], "/rows/0/walk")
        other_row = copy.deepcopy(cut)
        other_row["rows"][1]["walk"] = 5
        self.assertEqual(cg.cut_aware_difference(uncut, other_row)[0], "/rows/1/walk")

    def test_fields_cut_and_a_truncated_row(self):
        uncut = {"a": "<handle>", "evidence": {"a": 1}, "export": "x", "total": 2,
                 "rows": [{"key": "k", "prefix": "ACGTACGT", "labels": [1, 2, 3, 4]}, {"key": "j"}]}
        cut = {"a": "<handle>", "fields_cut": ["evidence", "export"], "total": 2, "row_truncated": True,
               "cut_fields": ["prefix", "labels"], "rows": [{"key": "k", "prefix": "ACG", "labels": [1, 2]}]}
        self.assertIsNone(cg.is_cut(cut, uncut))   # evidence / export named in fields_cut


class FiveWay(unittest.TestCase):
    def setUp(self):
        gf.install_graphlet_doc_stub()
        self.oracle = gf.install_oracle_stub()
        self.op = {"id": "walks-len", "tool": "graphlet_walks",
                   "args": {"arm": "left", "rank": "length", "n": 10, "max_bytes": 2048}, "paging": "all"}
        self.ops = [self.op]
        self.run_ = None

    def tearDown(self):
        if self.run_:
            self.run_.cleanup()

    def make(self, ops=None, ways=None, pre=None, **extra):
        self.run_ = gf.GRun(scenarios=[gf.scenario(ops=ops or self.ops, ways=ways, pre=pre, **extra)])
        self.run_.complete_scenario(ops=[])
        return self.run_

    def five(self, op=None, answer=None, w3=None, w5=True, ways=("W1", "W2", "W3", "W4"), changes=None,
             errors=None, sf=None, checked=None):
        op = op or self.op
        answer = answer if answer is not None else gf.walks_answer()
        mod, text = gf.install_graphlet_doc_stub(), self.run_.fixture.doc(0)
        digests = {"body_digest": mod.body_digest(text), "envelope_digest": mod.envelope_digest(text)}
        for way in ways:
            ans = copy.deepcopy(answer)
            if way == "W3":
                ans = copy.deepcopy(w3) if w3 is not None else dict(ans, handle=UUID)
            if changes and way in changes:
                changes[way](ans)
            err = (errors or {}).get(way)
            self.run_.record(way, SID, op["id"], op["tool"], op["args"], answer=ans, error=err, doc=digests,
                             service_fields=(sf if sf is not None else {"walk_display": "first_parent"}))
        if w5:
            self.run_.record("W5", SID, op["id"], op["tool"], op["args"],
                             answer={"rows": [{"length_bp": r["length_bp"]} for r in answer.get("rows", [])]},
                             oracle={"checked": checked if checked is not None else
                                     {"/rows/0/length_bp": "SPEC-core §7.1 path length"},
                                     "consistency_only": ["/rows/0/n_alive", "/rows/1/n_alive"],
                                     "inputs": ["full", "mgt"]})

    def result(self):
        res, pre = cg.compare_run(self.run_.dir)
        return res

    def op_entry(self, res=None, oid="walks-len"):
        return (res or self.result())["scenarios"][SID]["ops"][oid]

    def test_all_five_agree(self):
        self.make()
        self.five()
        res = self.result()
        e = self.op_entry(res)
        self.assertEqual(e["outcome"], "agree", e)
        self.assertEqual(e["ways"], {"W1": "answered", "W2": "answered", "W3": "answered", "W4": "answered",
                                     "W5": "answered"})
        self.assertEqual(sorted(p["pair"] for p in e["pairs"]),
                         ["W1-W2", "W1-W4", "W1-W5", "W2-W5", "W3-W4", "W3-W5", "W4-W5"])
        self.assertEqual(res["scenarios"][SID]["outcome"], "pass")
        self.assertEqual(res["summary"]["ops"]["agree"], 1)
        self.assertEqual(res["summary"]["per_tool"]["graphlet_walks"]["agree"], 1)
        self.assertEqual(res["coverage"]["graphlet_walks"]["checked"]["/rows/*/length_bp"]["calls"], 1)
        self.assertEqual(res["coverage"]["graphlet_walks"]["consistency_only"], {"/rows/*/n_alive": 2})

    def test_transport_disagreement_is_named(self):
        self.make()
        self.five(changes={"W2": lambda a: a["rows"][0].update(length_bp=121)})
        e = self.op_entry()
        self.assertEqual(e["outcome"], "disagree")
        self.assertEqual(e["disagreeing"], ["W1-W2", "W2-W5"])
        p = next(p for p in e["pairs"] if p["pair"] == "W1-W2")
        self.assertEqual((p["what"], p["pointer"], p["W1"], p["W2"]),
                         ("transport", "/answer/rows/0/length_bp", "120", "121"))

    def test_service_wrapper_and_fetch_path(self):
        self.make()
        self.five(changes={"W3": lambda a: a["rows"][1].update(end="branch")})
        e = self.op_entry()
        self.assertEqual(e["disagreeing"], ["W3-W4"])   # W5 does not check /rows/1/end
        self.assertEqual(next(p for p in e["pairs"] if p["pair"] == "W3-W4")["what"], "service wrapper")
        self.run_.cleanup()
        self.make()
        self.five(changes={"W4": lambda a: a.update(total=3)})
        self.assertEqual(self.op_entry()["disagreeing"], ["W1-W4", "W3-W4"])

    def test_the_oracle_disagrees_on_a_checked_pointer(self):
        self.make()
        self.five()
        path = os.path.join(self.run_.dir, "arms", "oracle", "graphlet_ops.jsonl")
        with open(path) as fh:
            rec = json.loads(fh.read())
        rec["answer"]["rows"][0]["length_bp"] = 119
        with open(path, "w") as fh:
            fh.write(json.dumps(rec) + "\n")
        e = self.op_entry()
        self.assertEqual(e["disagreeing"], ["W1-W5", "W2-W5", "W3-W5", "W4-W5"])
        p = next(p for p in e["pairs"] if p["pair"] == "W1-W5")
        self.assertEqual((p["what"], p["kind"], p["pointer"]), (cg.W5_WHAT, "oracle", "/rows/0/length_bp"))
        self.assertIn("SPEC-core §7.1", p["detail"])

    def test_args_resolved_first(self):
        self.make()
        self.five()
        path = os.path.join(self.run_.dir, "arms", "mcp", "graphlet_ops.jsonl")
        with open(path) as fh:
            rec = json.loads(fh.read())
        rec["args_resolved"]["n"] = 11
        rec["args_resolved"]["access_token"] = "secret"   # never compared
        with open(path, "w") as fh:
            fh.write(json.dumps(rec) + "\n")
        e = self.op_entry()
        p = next(p for p in e["pairs"] if p["pair"] == "W3-W4")
        self.assertEqual((p["outcome"], p["kind"], p["pointer"]), ("disagree", "args_resolved", "/n"))

    def test_selectors(self):
        self.make()
        for way in ("W1", "W2", "W3", "W4", "W5"):
            self.run_.record(way, SID, "walks-len", "graphlet_walks", self.op["args"], answer=None,
                             selector_error="walk_with_continuation: no walk has a continuation")
        e = self.op_entry()
        self.assertEqual(e["outcome"], "inconclusive")
        self.assertIn("the selector failed in every way", e["reason"])
        self.run_.cleanup()
        self.make()
        self.five()
        self.run_.record("W3", SID, "walks-len", "graphlet_walks", self.op["args"], answer=None,
                         selector_error="/walks/0/id: no such path")
        e = self.op_entry()
        self.assertEqual(e["outcome"], "disagree")
        self.assertEqual(next(p for p in e["pairs"] if p["pair"] == "W3-W4")["kind"], "selector")

    def test_a_missing_way_is_inconclusive(self):
        self.make()
        self.five(ways=("W1", "W3", "W4"))
        e = self.op_entry()
        self.assertEqual((e["outcome"], e["ways"]["W2"]), ("inconclusive", "missing"))
        self.assertIn("no record of W2", e["reason"])

    def test_not_exposed_with_and_without_the_owners_reason(self):
        op = dict(self.op, id="splits", tool="graphlet_splits", args={"arm": "left", "n": 10, "max_bytes": 2048},
                  ways=["W1", "W2", "W4", "W5"], not_exposed={"W3": "library-only (owner, 2026-10-19)"})
        self.make(ops=[op])
        self.five(op=op, ways=("W1", "W2", "W4"), answer={"handle": "g_0123456789ab", "arm": "left", "total": 0,
                                                            "rows": []}, checked={})
        e = self.op_entry(oid="splits")
        self.assertEqual((e["outcome"], e["ways"]["W3"], e["not_exposed"]),
                         ("agree", "not_exposed", {"W3": "library-only (owner, 2026-10-19)"}))
        self.assertEqual(self.result()["summary"]["per_tool"]["graphlet_splits"]["not_exposed"], 1)
        # by the manifest's table (by tool, or by tool.argument)
        op2 = dict(op)
        op2.pop("not_exposed")
        scen = gf.scenario(ops=[op2])
        man = {"graphlets": {"not_exposed": {"graphlet_splits": "not exposed by the service"}}}
        self.assertEqual(cg.not_exposed_reason(scen, op2, man), "not exposed by the service")
        man = {"graphlets": {"not_exposed": {"graphlet_splits.n": "an argument the service lacks"}}}
        self.assertEqual(cg.not_exposed_reason(scen, op2, man), "an argument the service lacks")
        self.run_.cleanup()
        self.make(ops=[op2])
        self.five(op=op2, ways=("W1", "W2", "W4"), answer={"handle": "g_0123456789ab", "arm": "left", "total": 0,
                                                             "rows": []}, checked={})
        e = self.op_entry(oid="splits")
        self.assertEqual(e["outcome"], "inconclusive")
        self.assertIn("without the owner's reason", e["reason"])

    def test_precondition_failed(self):
        self.make(pre={"arms": ["right"]})
        self.five()
        res = self.result()
        self.assertEqual(self.op_entry(res)["outcome"], "precondition_failed")
        self.assertEqual(res["scenarios"][SID]["outcome"], "precondition_failed")

    def test_a_failed_precheck_fails_the_scenario_and_ops_are_still_compared(self):
        self.make()
        self.five()
        self.run_.write_echo(SID, dict(gf.strategy(), direction="left"))   # a later line wins: a drift
        res = self.result()
        self.assertEqual(self.op_entry(res)["outcome"], "agree")
        self.assertEqual(res["scenarios"][SID]["precheck"]["strategy"]["outcome"], "fail")
        self.assertEqual(res["scenarios"][SID]["outcome"], "fail")

    def test_errors(self):
        four = ["W1", "W2", "W3", "W4"]
        self.make(ways=four)
        err = {"code": "bad_argument", "fields": {"message": "rank: one of support, length"}}
        self.five(errors=dict((w, err) for w in ("W1", "W2", "W4")), w5=False)
        e = self.op_entry()
        self.assertEqual(e["outcome"], "disagree")   # W3 answered where the library refused
        p = next(p for p in e["pairs"] if p["pair"] == "W3-W4")
        self.assertIn("W4 answered the error bad_argument, W3 an answer", p["detail"])
        self.run_.cleanup()
        self.make(ways=four)
        other = {"code": "bad_argument", "fields": {"message": "another text"}}
        self.five(errors={"W1": err, "W2": err, "W3": other, "W4": err}, w5=False)
        self.assertEqual(self.op_entry()["outcome"], "agree")
        self.run_.cleanup()
        self.make(ways=four)
        self.five(errors={"W1": err, "W2": err, "W3": {"code": "bad_arm", "fields": {}}, "W4": err}, w5=False)
        p = next(p for p in self.op_entry()["pairs"] if p["pair"] == "W3-W4")
        self.assertEqual((p["outcome"], p["pointer"]), ("disagree", "/error/code"))

    def test_the_oracle_on_an_error(self):
        self.make()
        err = {"code": "bad_argument", "fields": {}}
        self.five(errors=dict((w, err) for w in ("W1", "W2", "W3", "W4")), w5=False)
        # FORMAT v2 F1: W5 states an expected error by declaring it (answer null)
        self.run_.record("W5", SID, "walks-len", "graphlet_walks", self.op["args"], answer=None,
                         oracle={"checked": {}, "consistency_only": ["/*"], "declared": "expected_error",
                                 "note": "a walk id out of range"})
        e = self.op_entry()
        self.assertEqual(e["outcome"], "agree", e)
        self.assertEqual(next(p for p in e["pairs"] if p["pair"] == "W1-W5")["detail"], "W5 declares expected_error")
        self.assertEqual(e["w5"]["declared"], "expected_error")
        # every way answered where W5 expects an error
        self.run_.cleanup()
        self.make()
        self.five(w5=False)
        self.run_.record("W5", SID, "walks-len", "graphlet_walks", self.op["args"], answer=None,
                         oracle={"checked": {}, "declared": "expected_error"})
        self.assertEqual(self.op_entry()["disagreeing"], ["W1-W5", "W2-W5", "W3-W5", "W4-W5"])

    def test_w5_is_judged_by_lane_os_compare(self):
        self.make()
        self.five()
        self.op_entry()
        # one call per way, with W5's whole record and the way's view (service fields stripped, handles symbolised)
        self.assertEqual(len(self.oracle.calls), 4)
        record, view = next((r, v) for r, v in self.oracle.calls if v.get("handle") == "<handle>")
        self.assertEqual(record["way"], "W5")
        self.assertEqual(record["oracle"]["checked"], {"/rows/0/length_bp": "SPEC-core §7.1 path length"})
        self.assertNotIn("walk_display", json.dumps([v for _, v in self.oracle.calls]))
        # what compare returns decides: a difference, not comparable, a failure to read
        cg.graphlet_oracle = gf.OracleStub([{"pointer": "/rows/0/ends", "clause": "SPEC §7.5.6 (unordered)",
                                             "oracle": ["a", "b"], "other": ["b"]}])
        e = self.op_entry()
        p = next(p for p in e["pairs"] if p["pair"] == "W1-W5")
        self.assertEqual((p["outcome"], p["pointer"], p["W1"], p["W5"]),
                         ("disagree", "/rows/0/ends", '["b"]', '["a","b"]'))
        self.assertIn("1 difference(s); SPEC §7.5.6 (unordered)", p["detail"])
        cg.graphlet_oracle = gf.OracleStub({"outcome": "not_comparable", "differences": []})
        e = self.op_entry()
        self.assertEqual((e["outcome"], next(p for p in e["pairs"] if p["pair"] == "W1-W5")["outcome"]),
                         ("agree", "not_comparable"))
        cg.graphlet_oracle = gf.OracleStub(ValueError("no oracle block"))
        e = self.op_entry()
        self.assertEqual(e["outcome"], "inconclusive")
        self.assertIn("graphlet_oracle.compare raised ValueError: no oracle block", e["reason"])

    def test_fitting_errors_are_not_comparable_against_w5_only(self):
        self.make()
        fit = {"code": "result_too_large", "fields": {"bytes": 3000, "max_bytes": 512}}
        self.five(errors=dict((w, fit) for w in ("W1", "W2", "W3", "W4")))
        e = self.op_entry()
        self.assertEqual(e["outcome"], "agree", e)
        w5 = [p for p in e["pairs"] if p["pair"].endswith("W5")]
        self.assertEqual(set(p["outcome"] for p in w5), {"not_comparable"})
        self.assertIn("W1 answered result_too_large: W5 makes no claim", w5[0]["detail"])
        self.assertEqual(self.oracle.calls, [])   # never handed to the oracle's comparison
        # between W1..W4 a fitting error is compared as usual: W3 answered where the library did not fit
        self.run_.cleanup()
        self.make()
        self.five(errors={"W1": fit, "W2": fit, "W4": fit})
        e = self.op_entry()
        self.assertEqual(e["disagreeing"], ["W3-W4"])   # W3's answer agrees with W5; the library ways' are fitting errors
        # op_busy is the service's worker state: not comparable against W5 on W3 only
        self.run_.cleanup()
        self.make()
        busy = {"code": "op_busy", "fields": {"retryable": True}}
        self.five(errors={"W3": busy})
        e = self.op_entry()
        self.assertEqual(next(p for p in e["pairs"] if p["pair"] == "W3-W5")["outcome"], "not_comparable")
        self.assertEqual(e["disagreeing"], ["W3-W4"])

    def test_without_lane_os_module_w5_is_not_compared(self):
        self.make()
        self.five()
        cg.graphlet_oracle = None
        e = self.op_entry()
        self.assertEqual(e["outcome"], "inconclusive")
        self.assertIn("graphlet_oracle.compare (lane O, FORMAT v2 F2) is not available", e["reason"])
        self.assertEqual(next(p for p in e["pairs"] if p["pair"] == "W1-W2")["outcome"], "agree")

    def test_cli(self):
        import contextlib
        import io
        self.make()
        self.five()
        with contextlib.redirect_stdout(io.StringIO()) as out:
            self.assertEqual(cg.main([self.run_.dir]), 0)
        self.assertIn("scenarios pass 1; ops agree 1", out.getvalue())

    def test_service_fields_stripped_and_whitelist_gaps_listed(self):
        self.make()
        self.five(changes={"W1": lambda a: a.update(brand_new=1), "W2": lambda a: a.update(brand_new=1),
                           "W4": lambda a: a.update(brand_new=1)})
        res = self.result()
        e = self.op_entry(res)
        # the library's new key is stripped from W3 only if W3 had it: W3 lacks it, so the service wrapper differs
        self.assertEqual(e["disagreeing"], ["W3-W4"])
        self.assertEqual(res["summary"]["whitelist_gaps"], {"graphlet_walks": ["brand_new"]})

    def test_local_compared_within_one_library_environment_only(self):
        self.make()
        self.five(changes={"W1": lambda a: a.update(local={"usage": {"work_units": 5}}),
                           "W2": lambda a: a.update(local={"usage": {"work_units": 6}}),
                           "W4": lambda a: a.update(local={"usage": {"work_units": 5}})})
        e = self.op_entry()
        self.assertEqual(e["disagreeing"], ["W1-W2"])   # one commit, interpreter and work model: local compared
        p = next(p for p in e["pairs"] if p["pair"] == "W1-W2")
        self.assertEqual(p["pointer"], "/answer/local/usage/work_units")
        self.assertEqual(next(p for p in e["pairs"] if p["pair"] == "W3-W4")["outcome"], "agree")

    def test_a_deliberate_cut(self):
        uncut_op = dict(self.op, id="walks-all", args=dict(self.op["args"], max_bytes=262144))
        cut_op = dict(self.op, id="walks-page", args=dict(self.op["args"], max_bytes=2048), paging="none",
                      cut_of="walks-all", ways=["W1", "W2", "W3", "W4"])
        self.make(ops=[uncut_op, cut_op])
        rows = gf.walks_answer()["rows"] + [dict(gf.walks_answer()["rows"][1], walk=2)]
        full = gf.walks_answer(rows=rows)
        full["total"] = 3
        self.five(op=uncut_op, answer=full)
        page = gf.walks_answer(rows=rows[:2])
        page["total"] = 3
        w3_page = dict(gf.walks_answer(rows=rows[:1]), handle=UUID, total=3)   # the service's page ends earlier
        self.five(op=cut_op, answer=page, w3=w3_page, w5=False)
        e = self.op_entry(oid="walks-page")
        self.assertEqual(e["outcome"], "agree", e)
        self.assertEqual(sorted(c["way"] for c in e["cuts"]), ["W1", "W2", "W3", "W4"])
        self.assertNotIn("W3-W4", [p["pair"] for p in e["pairs"]])
        # a W3 page that is no cut of W3's uncut answer
        path = os.path.join(self.run_.dir, "arms", "mcp", "graphlet_ops.jsonl")
        with open(path) as fh:
            lines = [json.loads(x) for x in fh]
        lines[-1]["answer"]["rows"][0]["length_bp"] = 5
        with open(path, "w") as fh:
            fh.write("".join(json.dumps(x) + "\n" for x in lines))
        e = self.op_entry(oid="walks-page")
        self.assertEqual((e["outcome"], e["disagreeing"]), ("disagree", ["W3-uncut"]))
        self.assertEqual(next(c for c in e["cuts"] if c["way"] == "W3")["pointer"], "/rows/0/length_bp")

    def test_the_service_cuts_a_row_the_library_did_not(self):
        self.make(ways=["W1", "W2", "W3", "W4"])
        lib = gf.walks_answer(rows=[dict(gf.walks_answer()["rows"][0], sequence="ACGTACGTACGT")])
        lib["total"] = 1
        w3 = copy.deepcopy(lib)
        w3.update(handle=UUID, row_truncated=True, cut_fields=["sequence"])
        w3["rows"][0].update(sequence="ACGTAC", _cut_fields=["sequence"])
        self.five(answer=lib, w3=w3, w5=False)
        e = self.op_entry()
        self.assertEqual(e["outcome"], "agree", e)
        # a library way that cut where the others did not is a disagreement (library ways cut identically)
        self.run_.cleanup()
        self.make(ways=["W1", "W2", "W3", "W4"])
        self.five(answer=lib, w3=dict(lib, handle=UUID), w5=False,
                  changes={"W2": lambda a: (a.update(row_truncated=True, cut_fields=["sequence"]),
                                            a["rows"][0].update(sequence="ACG", _cut_fields=["sequence"]))})
        self.assertEqual(self.op_entry()["disagreeing"], ["W1-W2"])
        # the service's cut row that is no prefix
        self.run_.cleanup()
        self.make(ways=["W1", "W2", "W3", "W4"])
        bad = copy.deepcopy(w3)
        bad["rows"][0]["sequence"] = "TTTT"
        self.five(answer=lib, w3=bad, w5=False)
        p = next(p for p in self.op_entry()["pairs"] if p["pair"] == "W3-W4")
        self.assertEqual((p["outcome"], p["pointer"]), ("disagree", "/answer/rows/0/sequence"))
        self.assertIn("no cut of the other's", p["detail"])

    def test_a_cut_without_its_uncut_reference_is_inconclusive(self):
        cut_op = dict(self.op, id="walks-page", cut_of="walks-all")
        self.make(ops=[cut_op])
        self.five(op=cut_op, w5=False)
        e = self.op_entry(oid="walks-page")
        self.assertEqual(e["outcome"], "inconclusive")
        self.assertIn("no uncut reference", e["reason"])

    def test_export_and_smoke(self):
        exp = {"id": "fasta", "tool": "graphlet_export", "args": {"format": "fasta", "arm": "left"},
               "ways": ["W1", "W2", "W3", "W4"]}
        lst = {"id": "list", "tool": "graphlet_list", "args": {"max_bytes": 2048}, "ways": ["W1", "W2", "W4"],
               "not_exposed": {"W3": "store bookkeeping (D-R15)"}}
        self.make(ops=[exp, lst])
        self.five(op=exp, answer={"path": "/x/g_0123456789ab.fa", "bytes": 100, "records": 4},
                  w3={"path": None, "bytes": 99, "records": 4, "export_id": "e1", "link": {"url": "u"},
                      "format": "fasta", "handle": UUID, "status": "done", "in_service": True}, w5=False)
        self.five(op=lst, answer={"total": 1, "rows": [{"handle": "g_1"}]}, ways=("W1", "W2", "W4"), w5=False)
        res = self.result()
        e = self.op_entry(res, "fasta")
        # path None vs a path: path is left out on every way; the service's receipt keys are its own
        self.assertEqual(e["outcome"], "agree", e)
        self.assertTrue(self.op_entry(res, "list")["smoke"])
        self.assertEqual(res["summary"]["ops"]["agree"], 1)   # the smoke is outside the five-way counts
        self.assertEqual(res["summary"]["per_tool"]["graphlet_list"]["smoke"], 1)

    def test_drift_lane_is_listed(self):
        self.make()
        self.five()
        ans = dict(gf.walks_answer(), total=9)
        self.run_.record("drift", SID, "walks-len", "graphlet_walks", self.op["args"], answer=ans,
                         library={"commit": "d" * 40, "python": "3.13.1", "hashseed": 1, "work_model": 3})
        res = self.result()
        self.assertEqual(res["drift"]["compared"], 1)
        self.assertEqual(res["drift"]["differ"][0]["pointer"], "/answer/total")
        self.assertEqual((res["drift"]["differ"][0]["W4"], res["drift"]["differ"][0]["drift"]), ("2", "9"))
        self.assertEqual(res["scenarios"][SID]["outcome"], "pass")   # report-only

    def test_compute_writes_and_load_or_compute_recomputes_when_stale(self):
        self.make()
        self.five()
        res = cg.compute(self.run_.dir)
        self.assertTrue(os.path.exists(os.path.join(self.run_.dir, cg.COMPARE_FILE)))
        self.assertTrue(os.path.exists(gp.precheck_path(self.run_.dir)))
        self.assertEqual(cg.stale(self.run_.dir, res), [])
        again = cg.load_or_compute(self.run_.dir)
        self.assertEqual(again["generated_utc"], res["generated_utc"])
        self.run_.record("W2", SID, "walks-len", "graphlet_walks", self.op["args"], answer=dict(gf.walks_answer(),
                                                                                                total=7))
        self.assertEqual(cg.stale(self.run_.dir, res), ["arms/dr-lib/graphlet_ops.jsonl"])
        fresh = cg.load_or_compute(self.run_.dir)
        self.assertEqual(fresh["scenarios"][SID]["ops"]["walks-len"]["outcome"], "disagree")

    def test_unplanned_records_are_counted(self):
        self.make()
        self.five()
        self.run_.record("W1", SID, "not-in-plan", "graphlet_walks", self.op["args"], answer={})
        res = self.result()
        self.assertEqual(res["summary"]["unplanned_records"], 1)


if __name__ == "__main__":
    unittest.main()
