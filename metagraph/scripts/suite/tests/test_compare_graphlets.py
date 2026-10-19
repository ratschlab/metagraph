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
        self.assertIn("a value other than the uncut one", cg.is_cut(other, self.UNCUT, opt)[1])   # no cut names it
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
        self.assertIn("not a leading part", d[3])
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
               "cut_fields": ["prefix", "labels"],
               "rows": [{"key": "k", "prefix": "ACG", "labels": [1, 2], "_cut_fields": ["prefix", "labels"]}]}
        self.assertIsNone(cg.is_cut(cut, uncut))   # evidence / export named in fields_cut, the row's fields in A3's marker
        # without the row's own marker (FORMAT v2 A3) its fields are not cut fields
        bare = copy.deepcopy(cut)
        del bare["rows"][0]["_cut_fields"]
        self.assertEqual(cg.is_cut(bare, uncut)[0], "/rows/0/labels")


class RunnerCodes(unittest.TestCase):
    """VR1: every runner-side code of FORMAT v2.3 B2 is inconclusive, never compared."""

    # FORMAT v2.3 B2 (FORMAT-v8.md, 2026-10-19): the runner-side error codes, selector_error being its own field
    FORMAT_V23_B2 = ("no_document", "no_graphlet", "document_refused", "exception", "non_finite", "unknown_tool",
                     "paging_runaway", "assembly_error", "timeout", "child_died")

    def test_the_constant_is_b2s_list(self):
        # no_graphlet is a failed seed's answer (D3): compared by its code, not a runner failure
        self.assertEqual(sorted(cg.RUNNER_ERRORS + ("no_graphlet",)), sorted(self.FORMAT_V23_B2))
        self.assertIn("no_graphlet", cg.ERROR_CODE_ONLY)

    def test_the_verifiers_gate_case(self):
        """probe_n1_gate: graphlet_support (library-only, W5 phase B) where W1, W2 and W4 hit a runner-side code, W5
        written by lane O's run step: the op is inconclusive and the gate blocks, for every code."""
        import gate
        import gate_policy
        import graphlet_oracle as go
        gf.install_graphlet_doc_stub()
        saved = cg.graphlet_oracle
        cg.graphlet_oracle = go
        walks = {"id": "walks-len", "tool": "graphlet_walks", "paging": "all", "ways": ["W1", "W2", "W3", "W4"],
                 "args": {"arm": "left", "rank": "length", "n": 1, "max_bytes": 2048}}
        sup = {"id": "sup", "tool": "graphlet_support", "args": {"arm": "left", "walk": 0, "max_bytes": 2048},
               "paging": "all", "ways": ["W1", "W2", "W4", "W5"],
               "not_exposed": {"W3": "library-only (owner, DESIGN-v8 §2)"}}
        try:
            for code in ("timeout", "child_died", "non_finite", "exception"):
                scen = gf.scenario(ops=[walks, sup])
                run = gf.GRun(scenarios=[scen])
                try:
                    run.complete_scenario(ops=[("walks-len", "graphlet_walks", walks["args"], gf.walks_answer())])
                    for w in ("W1", "W2", "W4"):
                        run.record(w, scen["id"], "sup", "graphlet_support", sup["args"],
                                   error={"code": code, "fields": {"timeout_s": 120}})
                    go.write_run(run.dir, quiet=True, force=True)
                    block = gate_policy.graphlet_block(run.manifest, {scen["id"]})
                    pol = {"graphlets": block, "manifest": {"content_sha256": run.manifest["content_sha256"]}}
                    probs, _ = gate.graphlet_problems(run.dir, pol, run.manifest, [])
                    e = cg.compute(run.dir, run.manifest, write=False)["scenarios"][scen["id"]]["ops"]["sup"]
                    self.assertEqual((e["outcome"], e["ways"]["W1"]), ("inconclusive", "runner_error"), code)
                    self.assertIn("runner: %s" % code, e["reason"])
                    self.assertTrue(any("op sup (graphlet_support) inconclusive" in p["detail"] for p in probs),
                                    (code, probs))
                finally:
                    run.cleanup()
        finally:
            cg.graphlet_oracle = saved


class CutRows(unittest.TestCase):
    """VC1: the cut-row rule on pairs with W3 is lane O's, on the fields a row names in _cut_fields only (the verifier's
    probes P1, P1b, P1c, P3 as negatives)."""

    LABELS = [{"name": "NC_000962.3", "ref": "c:0"}, {"name": "NC_000913.3", "ref": "c:1"}]

    def setUp(self):
        import graphlet_oracle
        self.saved = cg.graphlet_oracle
        cg.graphlet_oracle = graphlet_oracle

    def tearDown(self):
        cg.graphlet_oracle = self.saved

    def pair(self, w3, lib):
        ra = {"tool": "graphlet_walks", "args_resolved": {"arm": "left"}, "answer": w3}
        rb = {"tool": "graphlet_walks", "args_resolved": {"arm": "left"}, "answer": lib}
        return cg.compare_records("W3", ra, "W4", rb, "service wrapper")

    def answers(self):
        lib = gf.walks_answer()
        lib["rows"][0]["labels_full"] = copy.deepcopy(self.LABELS)
        lib["rows"][0]["sequence"] = "ACGTACGT"
        w3 = copy.deepcopy(lib)
        w3["handle"] = UUID
        return lib, w3

    def test_a_named_field_cut_to_a_leading_part(self):
        lib, w3 = self.answers()
        w3["rows"][0].update(labels_full=self.LABELS[:1], sequence="ACGT", _cut_fields=["labels_full", "sequence"])
        w3.update(row_truncated=True, cut_fields=["labels_full", "sequence"])
        self.assertEqual(self.pair(w3, lib)["outcome"], "agree")
        # the same rule in the other direction: a library-side cut against a whole W3 row
        lib2, w32 = self.answers()
        lib2["rows"][0].update(sequence="ACGTA", _cut_fields=["sequence"])
        self.assertEqual(self.pair(w32, lib2)["outcome"], "agree")
        # both rows cut: compared on their common leading part
        lib3, w33 = self.answers()
        lib3["rows"][0].update(sequence="ACGTAC", _cut_fields=["sequence"])
        w33["rows"][0].update(sequence="ACG", _cut_fields=["sequence"])
        self.assertEqual(self.pair(w33, lib3)["outcome"], "agree")

    def test_p1_an_unnamed_field_cut(self):
        lib, w3 = self.answers()
        w3["rows"][0].update(end="dead", labels_full=self.LABELS[:1], _cut_fields=["labels_full"])
        p = self.pair(w3, lib)
        self.assertEqual((p["outcome"], p["pointer"]), ("disagree", "/answer/rows/0/end"))
        self.assertIn("a field no row names cut", p["detail"])

    def test_p1b_a_nested_prefix_inside_a_cut_list(self):
        lib, w3 = self.answers()
        w3["rows"][0].update(labels_full=[{"name": "NC_0", "ref": "c:0"}], _cut_fields=["labels_full"])
        p = self.pair(w3, lib)
        self.assertEqual((p["outcome"], p["pointer"]), ("disagree", "/answer/rows/0/labels_full"))
        self.assertIn("not a leading part", p["detail"])
        self.assertIn("not a cut of the full value", p["detail"])   # lane O's cut_problem's own words

    def test_p1c_a_field_emptied_that_the_row_does_not_name(self):
        lib, w3 = self.answers()
        w3["rows"][0].update(sequence="ACGT", labels_full=[], _cut_fields=["sequence"])
        p = self.pair(w3, lib)
        self.assertEqual(p["outcome"], "disagree")
        self.assertTrue(p["pointer"].startswith("/answer/rows/0/labels_full"), p)
        self.assertIn("a field no row names cut", p["detail"])

    def test_an_empty_cut_of_a_named_field_agrees(self):
        # the library's _shrink can cut a named field to empty (a one-item list to []): lane O's one rule accepts it
        lib, w3 = self.answers()
        w3["rows"][0].update(labels_full=[], sequence="", _cut_fields=["labels_full", "sequence"])
        self.assertEqual(self.pair(w3, lib)["outcome"], "agree")
        lib2, w32 = self.answers()   # the same rule in the other direction (P3), and both rows cut (P2)
        lib2["rows"][0].update(labels_full=[], _cut_fields=["labels_full"])
        self.assertEqual(self.pair(w32, lib2)["outcome"], "agree")
        lib3, w33 = self.answers()
        lib3["rows"][0].update(labels_full=self.LABELS[:1], _cut_fields=["labels_full"])
        w33["rows"][0].update(labels_full=[], _cut_fields=["labels_full"])
        self.assertEqual(self.pair(w33, lib3)["outcome"], "agree")

    def test_an_empty_value_in_an_unnamed_field_is_refused(self):
        lib, w3 = self.answers()
        w3["rows"][0].update(labels_full=[], sequence="ACGT", _cut_fields=["sequence"])
        p = self.pair(w3, lib)
        self.assertEqual(p["outcome"], "disagree")
        self.assertTrue(p["pointer"].startswith("/answer/rows/0/labels_full"), p)
        self.assertIn("a field no row names cut", p["detail"])
        lib2, w32 = self.answers()
        w32["rows"][0].update(sequence="", _cut_fields=["labels_full"])
        p = self.pair(w32, lib2)
        self.assertEqual((p["outcome"], p["pointer"]), ("disagree", "/answer/rows/0/sequence"))

    def test_a_scalar_named_cut_and_library_pairs_strict(self):
        lib, w3 = self.answers()
        w3["rows"][0].update(length_bp=12, _cut_fields=["length_bp"])
        self.assertEqual(self.pair(w3, lib)["outcome"], "disagree")
        lib2, _ = self.answers()
        cut = copy.deepcopy(lib2)
        cut["rows"][0].update(sequence="ACGT", _cut_fields=["sequence"])
        ra = {"tool": "graphlet_walks", "args_resolved": {}, "answer": cut}
        rb = {"tool": "graphlet_walks", "args_resolved": {}, "answer": lib2}
        self.assertEqual(cg.compare_records("W1", ra, "W2", rb, "transport")["outcome"], "disagree")

    def test_without_lane_os_rule_a_cut_row_is_not_compared(self):
        lib, w3 = self.answers()
        w3["rows"][0].update(sequence="ACGT", _cut_fields=["sequence"])
        cg.graphlet_oracle = None
        p = self.pair(w3, lib)
        self.assertEqual((p["outcome"], p["kind"]), ("not_compared", "oracle_unavailable"))


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
        scen = gf.scenario(ops=[op2], not_exposed={"graphlet_splits": "not exposed by the service"})
        self.assertEqual(cg.not_exposed_reason(scen, op2), "not exposed by the service")
        scen = gf.scenario(ops=[op2], not_exposed={"graphlet_splits.n": "an argument the service lacks"})
        self.assertEqual(cg.not_exposed_reason(scen, op2), "an argument the service lacks")
        # a manifest-level table would be outside content_sha256: not read (FORMAT v2 D2)
        man = {"graphlets": {"not_exposed": {"graphlet_splits": "x"}}}
        self.assertIsNone(cg.not_exposed_reason(gf.scenario(ops=[op2]), op2, man))
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

    def test_record_conventions_of_format_v2(self):
        # B1: the runners' symbols are kept as they are
        a = {"handle": "<handle>", "evidence": {"view": {"of": "<handle:walks-sub>"}}}
        self.assertEqual(cg.canonical_answer("graphlet_walks", a), a)
        # B2: a runner's own failure is not compared, the op inconclusive
        self.make()
        self.five(errors={"W2": {"code": "exception", "fields": {"type": "KeyError"}}})
        e = self.op_entry()
        self.assertEqual((e["outcome"], e["ways"]["W2"]), ("inconclusive", "runner_error"))
        self.assertIn("W2 did not run the call (runner: exception)", e["reason"])
        # D3: a failed seed's no_graphlet is compared by its code (the runner's and the service's fields differ)
        self.run_.cleanup()
        self.make(ways=["W1", "W2", "W3", "W4"])
        lib = {"code": "no_graphlet", "fields": {"outcome": "failed"}}
        self.five(errors={"W1": lib, "W2": lib, "W4": lib, "W3": {"code": "no_graphlet", "fields": {"reason": "x"}}},
                  w5=False)
        self.assertEqual(self.op_entry()["outcome"], "agree")
        # B4: local is compared only within one library configuration (limits)
        r1 = {"library": {"commit": "c", "python": "3", "work_model": 2, "limits": {"work_units": 5}}}
        r2 = {"library": {"commit": "c", "python": "3", "work_model": 2, "limits": {"work_units": 6}}}
        self.assertFalse(cg.same_library_env(r1, r2))
        self.assertTrue(cg.same_library_env(r1, copy.deepcopy(r1)))

    def test_vc7_every_way_is_reached_by_a_pair(self):
        self.assertEqual(cg.op_pairs(["W1", "W2", "W3", "W4", "W5"]),
                         [("W1", "W2", "transport"), ("W1", "W4", "fetch path"), ("W3", "W4", "service wrapper")])
        self.assertEqual(cg.op_pairs(["W1", "W3"]), [("W1", "W3", "service wrapper (no W4)")])
        self.assertEqual(cg.op_pairs(["W2", "W3", "W4"]), [("W3", "W4", "service wrapper"), ("W2", "W4", "library ways")])
        # the verifier's P4: ways W1 and W3, W3 answers no rows against W1's two
        op = dict(self.op, ways=["W1", "W3"])
        self.make(ops=[op])
        self.five(op=op, ways=("W1", "W3"), w3=dict(gf.walks_answer(rows=[]), handle=UUID), w5=False)
        e = self.op_entry()
        self.assertEqual((e["outcome"], e["disagreeing"]), ("disagree", ["W1-W3"]))

    def test_vc6_empty_records_and_incomplete_listings(self):
        self.make(ways=["W1", "W2", "W3", "W4"])
        for way in ("W1", "W2", "W3", "W4"):
            self.run_.record(way, SID, "walks-len", "graphlet_walks", self.op["args"], answer=None)
        e = self.op_entry()
        self.assertEqual((e["outcome"], e["ways"]["W1"]), ("inconclusive", "runner_error"))
        self.assertIn("did not run the call (runner: no_answer)", e["reason"])
        self.run_.cleanup()
        self.make(ways=["W1", "W2", "W3", "W4"])
        self.five(w5=False)
        path = os.path.join(self.run_.dir, "arms", "dl-lib", "graphlet_ops.jsonl")
        with open(path) as fh:
            rec = json.loads(fh.read())
        rec["paging_complete"] = False
        with open(path, "w") as fh:
            fh.write(json.dumps(rec) + "\n")
        e = self.op_entry()
        self.assertEqual(e["outcome"], "inconclusive")
        self.assertIn("W1's listing stopped before its last page", e["reason"])
        # an op with a declared page cap states its own extent
        scen = gf.scenario(ops=[dict(self.op, paging={"pages": 1})], ways=["W1", "W2", "W3", "W4"])
        recs = gp.all_records(self.run_.dir)[0]
        pre = {"outcome": "pass", "preconditions": {"outcome": "pass"}}
        self.assertEqual(cg.compare_op(scen, scen["graphlet"]["ops"][0], recs, pre, None)["outcome"], "agree")

    def test_local_is_dropped_only_where_the_library_puts_it(self):
        a = {"handle": "h", "local": {"usage": 1}, "rows": [{"local": 1}], "evidence": {"limitations": {"local": ["x"]}}}
        self.assertEqual(cg.canonical_answer("graphlet_walks", a),
                         {"handle": "<handle>", "rows": [{"local": 1}], "evidence": {"limitations": {}}})

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
        self.assertIn("expected error", next(p for p in e["pairs"] if p["pair"] == "W1-W5")["detail"])
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
        self.oracle = cg.graphlet_oracle = gf.OracleStub()   # the calls compare_w5 makes
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
        cg.graphlet_oracle = gf.OracleStub({"outcome": "error", "detail": "W5 record unreadable", "differences": []})
        e = self.op_entry()   # VC8 (the verifier's P7): an unknown outcome is no agreement
        self.assertEqual(e["outcome"], "inconclusive")
        self.assertIn("graphlet_oracle.compare answered the outcome 'error': W5 record unreadable", e["reason"])
        cg.graphlet_oracle = gf.OracleStub({"outcome": "inconclusive", "detail": "an unresolved handle symbol",
                                            "differences": []})
        e = self.op_entry()   # lane O's inconclusive: W5 could not judge a call it should (blocks)
        self.assertEqual(e["outcome"], "inconclusive")
        self.assertIn("W5 could not judge the call (graphlet_oracle.compare): an unresolved handle symbol", e["reason"])
        cg.graphlet_oracle = gf.OracleStub(ValueError("no oracle block"))
        e = self.op_entry()
        self.assertEqual(e["outcome"], "inconclusive")
        self.assertIn("graphlet_oracle.compare raised ValueError: no oracle block", e["reason"])

    def test_fitting_errors_on_an_op_without_cut_of_are_inconclusive(self):
        # the coordinator's decision on R2-7: the plan failed to fit, the op tested nothing (blocks)
        self.make()
        fit = {"code": "result_too_large", "fields": {"bytes": 3000, "max_bytes": 512}}
        self.five(errors=dict((w, fit) for w in ("W1", "W2", "W3", "W4")))
        e = self.op_entry()
        self.assertEqual(e["outcome"], "inconclusive", e)
        self.assertIn("the plan failed to fit: W1 answered result_too_large, W2 answered result_too_large", e["reason"])
        # one way only, the others answering: inconclusive, the differing pairs named besides
        self.run_.cleanup()
        self.make()
        self.five(errors={"W1": fit, "W2": fit, "W4": fit})
        e = self.op_entry()
        self.assertEqual((e["outcome"], e["disagreeing"]), ("inconclusive", ["W3-W4"]))
        self.assertIn("the pairs W3-W4 differ besides", e["reason"])

    def test_fitting_errors_on_a_cut_op_keep_f3s_rule(self):
        uncut_op = dict(self.op, id="walks-all", args=dict(self.op["args"], max_bytes=262144))
        cut_op = dict(self.op, id="walks-page", args=dict(self.op["args"], max_bytes=64), paging="none",
                      cut_of="walks-all")
        self.make(ops=[uncut_op, cut_op])
        self.five(op=uncut_op)
        fit = {"code": "result_too_large", "fields": {"bytes": 3000, "max_bytes": 64}}
        self.five(op=cut_op, errors=dict((w, fit) for w in ("W1", "W2", "W3", "W4")))
        e = self.op_entry(oid="walks-page")
        self.assertEqual(e["outcome"], "agree", e)
        self.assertEqual(set(c["outcome"] for c in e["cuts"]), {"agree"})
        self.assertIn("the ceiling's fitting error result_too_large", e["cuts"][0]["detail"])
        w5 = [p for p in e["pairs"] if p["pair"].endswith("W5")]
        self.assertEqual(set(p["outcome"] for p in w5), {"not_comparable"})
        self.assertEqual(getattr(self.oracle, "calls", []), [])
        self.assertNotIn("W3-W5", [p["pair"] for p in e["pairs"]])   # a cut op's W3 is judged by its own cut check

    def test_vr2_only_a_ceilings_refusal_is_a_cut(self):
        uncut_op = dict(self.op, id="walks-all", args=dict(self.op["args"], max_bytes=262144))
        cut_op = dict(self.op, id="walks-page", args=dict(self.op["args"], max_bytes=64), paging="none",
                      cut_of="walks-all")
        fit = {"code": "result_too_large", "fields": {"bytes": 3000, "max_bytes": 64}}
        for code in ("op_busy", "op_timeout", "op_memory"):
            self.make(ops=[uncut_op, cut_op])
            self.five(op=uncut_op)
            self.five(op=cut_op, errors={"W1": fit, "W2": fit, "W4": fit, "W3": {"code": code, "fields": {}}})
            e = self.op_entry(oid="walks-page")
            self.assertEqual(e["outcome"], "inconclusive", (code, e))
            self.assertIn("not a cut: W3 answered %s" % code, e["reason"])
            self.run_.cleanup()
        self.make(ops=[uncut_op, cut_op])
        self.five(op=uncut_op)
        too = {"code": "body_too_large", "fields": {}}
        self.five(op=cut_op, errors={"W1": fit, "W2": fit, "W4": fit, "W3": too})
        e = self.op_entry(oid="walks-page")
        self.assertEqual(next(c for c in e["cuts"] if c["way"] == "W3")["outcome"], "agree")

    def test_a_torn_w1_line_is_w5s_inconclusive_not_differing_args(self):
        import graphlet_oracle as go
        cg.graphlet_oracle = go
        self.make()
        self.five(w5=False)
        self.run_.record("W5", SID, "walks-len", "graphlet_walks", {}, answer=None,
                         oracle={"checked": {}, "consistency_only": ["/*"], "declared": "unreadable_w1",
                                 "note": "W1's line 3 is torn: its record cannot be read"})
        path = os.path.join(self.run_.dir, "arms", "oracle", "graphlet_ops.jsonl")
        with open(path) as fh:
            rec = json.loads(fh.read())
        rec["args_resolved"] = None
        with open(path, "w") as fh:
            fh.write(json.dumps(rec) + "\n")
        e = self.op_entry()
        p = next(p for p in e["pairs"] if p["pair"] == "W1-W5")
        self.assertNotEqual(p.get("kind"), "args_resolved")
        self.assertEqual(e["outcome"], "inconclusive")
        self.assertIn("W5 could not judge the call", e["reason"])

    def test_vr4_a_failed_seed_scenario_can_pass(self):
        """A failed seed (FORMAT v2 D3): no way has a graphlet, W1/W2/W4 answer the runner's no_graphlet, W3 the
        service's, and W5 expects no_graphlet (lane O, VERIFY-C-R2 VR4: from the full twin's seed error): the op agrees
        and the scenario passes its pre-checks and its op."""
        import graphlet_oracle as go
        cg.graphlet_oracle = go
        scen = gf.scenario(pre={"no_graphlet": True})
        self.run_ = gf.GRun(scenarios=[scen])
        self.run_.write_twin(SID, {"results": [{"seed": {"length_bp": 31}, "error": "seed_failed"}]})
        self.run_.write_echo(SID, scen["request"]["strategy"])
        self.run_.write_mcp_runs()
        args = scen["graphlet"]["ops"][0]["args"]
        for way in ("W1", "W2", "W4"):
            self.run_.record(way, SID, "walks-len", "graphlet_walks", args,
                             error={"code": "no_graphlet", "fields": {"outcome": "failed"}})
        self.run_.record("W3", SID, "walks-len", "graphlet_walks", args,
                         error={"code": "no_graphlet", "fields": {"reason": "the seed's retrieval failed"}})
        self.run_.record("W5", SID, "walks-len", "graphlet_walks", args, answer=None,
                         oracle={"checked": {}, "consistency_only": ["/*"], "declared": "expected_error",
                                 "expected_error": "no_graphlet", "note": "the full twin's seed failed"})
        res = self.result()
        e = res["scenarios"][SID]["ops"]["walks-len"]
        self.assertEqual(e["outcome"], "agree", e)
        self.assertEqual(res["scenarios"][SID]["precheck"]["outcome"], "pass", res["scenarios"][SID]["precheck"])
        self.assertEqual(res["scenarios"][SID]["outcome"], "pass")

    def test_vr4_a_failed_seed_with_lane_os_run_step(self):
        """The same failed seed with W5 written by lane O's own run step (graphlet_oracle.write_run, 846eeb97d: W5
        expects no_graphlet from the full twin's seed error): the op agrees and the scenario passes."""
        import graphlet_oracle as go
        cg.graphlet_oracle = go
        scen = gf.scenario(pre={"no_graphlet": True})
        self.run_ = gf.GRun(scenarios=[scen])
        self.run_.write_twin(SID, {"results": [{"seed": {"length_bp": 31}, "error": "seed_failed"}]})
        self.run_.write_echo(SID, scen["request"]["strategy"])
        self.run_.write_mcp_runs()
        args = scen["graphlet"]["ops"][0]["args"]
        for way in ("W1", "W2", "W4"):
            self.run_.record(way, SID, "walks-len", "graphlet_walks", args,
                             error={"code": "no_graphlet", "fields": {"outcome": "failed"}})
        self.run_.record("W3", SID, "walks-len", "graphlet_walks", args,
                         error={"code": "no_graphlet", "fields": {"reason": "the seed's retrieval failed"}})
        go.write_run(self.run_.dir, quiet=True, force=True)
        res = self.result()
        e = res["scenarios"][SID]["ops"]["walks-len"]
        self.assertEqual(e["outcome"], "agree", e)
        self.assertEqual(res["scenarios"][SID]["outcome"], "pass", res["scenarios"][SID])

    def test_a_run_level_w5_record(self):
        self.make()
        self.five()
        self.run_.append(os.path.join("arms", "oracle", "graphlet_ops.jsonl"),
                         {"run_id": "test", "scenario": None, "way": "W5", "op": None, "tool": None,
                          "run_level": True, "args_resolved": None, "answer": None,
                          "error": {"code": "unreadable_w1_lines", "fields": {"lines": 2}},
                          "oracle": {"checked": {}, "declared": "unreadable_w1",
                                     "note": "2 torn W1 line(s) name no scenario and op"}})
        res = cg.compute(self.run_.dir, write=False)
        self.assertEqual(res["run_level"][0]["declared"], "unreadable_w1")
        self.assertEqual(res["summary"]["unplanned_records"], 0)
        import report
        word, reason = report.graphlet_verdict(res)
        self.assertEqual(word, "INCONCLUSIVE")
        self.assertIn("W5's run-level record: unreadable_w1: 2 torn W1 line(s)", reason)

    def test_lane_os_compare_on_the_vo4_cases(self):
        import graphlet_oracle as go
        if not callable(getattr(go, "compare", None)):
            self.skipTest("lane O's graphlet_oracle.compare is not on this tree")
        cg.graphlet_oracle = go
        op = dict(self.op, id="walk", tool="graphlet_walk", args={"arm": "left", "walk": 1, "max_bytes": 2048},
                  ways=["W1", "W2", "W3", "W4", "W5"])
        ans = {"handle": "g_0123456789ab", "arm": "left", "walk": 1, "length_bp": 120, "total": 1,
               "rows": [{"segment": 0, "ends": [{"at_bp": 5, "ref": "c:1"}, {"at_bp": 5, "ref": "c:0"}]}]}
        self.make(ops=[op])
        self.five(op=op, answer=ans, w5=False)
        # W5 lists the ends in another order and names them unordered; it also names a key that must be absent
        self.run_.record("W5", SID, "walk", "graphlet_walk", op["args"],
                         answer={"rows": [{"ends": [{"at_bp": 5, "ref": "c:0"}, {"at_bp": 5, "ref": "c:1"}]}]},
                         oracle={"checked": {"/rows/0/ends": "SPEC §7.5.6"}, "unordered": ["/rows/0/ends"],
                                 "absent": {"/continuation": "SPEC §7.1: no continuation"}})
        e = self.op_entry(oid="walk")
        self.assertEqual(e["outcome"], "agree", e)
        self.assertEqual(e["w5"]["unordered"], ["/rows/0/ends"])
        # a continuation where W5 says none
        self.run_.cleanup()
        self.make(ops=[op])
        self.five(op=op, answer=dict(ans, continuation={"length_bp": 3}), w5=False)
        self.run_.record("W5", SID, "walk", "graphlet_walk", op["args"], answer={"length_bp": 120},
                         oracle={"checked": {"/length_bp": "SPEC §7.5.6"},
                                 "absent": {"/continuation": "SPEC §7.1: no continuation"}})
        p = next(p for p in self.op_entry(oid="walk")["pairs"] if p["pair"] == "W1-W5")
        self.assertEqual((p["outcome"], p["pointer"]), ("disagree", "/continuation"))
        # an expected error of a stated code, and a call W5 makes no claim on
        self.run_.cleanup()
        self.make(ops=[op])
        err = {"code": "bad_argument", "fields": {}}
        self.five(op=op, answer=ans, errors=dict((w, err) for w in ("W1", "W2", "W3", "W4")), w5=False)
        self.run_.record("W5", SID, "walk", "graphlet_walk", op["args"], answer=None,
                         oracle={"checked": {}, "declared": "expected_error", "expected_error": "bad_arm",
                                 "note": "an arm not retrieved"})
        p = next(p for p in self.op_entry(oid="walk")["pairs"] if p["pair"] == "W1-W5")
        self.assertEqual((p["outcome"], p["pointer"], p["W1"], p["W5"]), ("disagree", "/error/code",
                                                                         '"bad_argument"', '"bad_arm"'))
        self.run_.record("W5", SID, "walk", "graphlet_walk", op["args"], answer=None,
                         oracle={"checked": {}, "declared": "not_implemented_phase_b"})
        e = self.op_entry(oid="walk")
        self.assertEqual((e["outcome"], next(p for p in e["pairs"] if p["pair"] == "W1-W5")["outcome"]),
                         ("agree", "not_checked"))

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
        self.assertIn("not a leading part", p["detail"])

    def cut_scenario(self):
        uncut_op = dict(self.op, id="walks-all", args=dict(self.op["args"], max_bytes=262144),
                        ways=["W1", "W2", "W3", "W4"])
        cut_op = dict(self.op, id="walks-page", args=dict(self.op["args"], max_bytes=512), paging="none",
                      cut_of="walks-all", ways=["W1", "W2", "W3", "W4"])
        self.make(ops=[uncut_op, cut_op])
        full = gf.walks_answer()
        self.five(op=uncut_op, answer=full, w5=False)
        return cut_op, full

    def test_vc4_a_cut_op_that_shows_no_cut_fails(self):
        # P5: cut == uncut on every way
        cut_op, full = self.cut_scenario()
        self.five(op=cut_op, answer=full, w5=False)
        e = self.op_entry(oid="walks-page")
        self.assertEqual(e["outcome"], "disagree")
        self.assertEqual(sorted(e["disagreeing"]), ["W1-uncut", "W2-uncut", "W3-uncut", "W4-uncut"])
        self.assertTrue(all(c["detail"].startswith("not cut:") for c in e["cuts"]))

    def test_vc4_p6_strings_cut_without_a_marker(self):
        cut_op, full = self.cut_scenario()
        page = copy.deepcopy(full)
        page["evidence"]["support"] = "k"
        page["rows"] = page["rows"][:1]
        page["rows"][0].update(end="", labels_full=[])
        self.five(op=cut_op, answer=page, w5=False)
        e = self.op_entry(oid="walks-page")
        self.assertEqual(e["outcome"], "disagree")
        c = next(c for c in e["cuts"] if c["way"] == "W1")
        self.assertEqual(c["pointer"], "/evidence/support")   # one row of two is evidence, the strings are not cuts

    def test_vc4_p6b_a_page_without_rows(self):
        cut_op, full = self.cut_scenario()
        self.five(op=cut_op, answer=dict(full, rows=[]), w5=False)
        c = next(c for c in self.op_entry(oid="walks-page")["cuts"] if c["way"] == "W1")
        self.assertEqual((c["outcome"], c["pointer"]), ("disagree", "/rows"))
        self.assertIn("a page without a row", c["detail"])

    def test_vc4_a_short_page_is_a_cut(self):
        cut_op, full = self.cut_scenario()
        page = copy.deepcopy(full)
        page["rows"] = page["rows"][:1]
        self.five(op=cut_op, answer=page, w5=False)   # one row of two: a page cut by its size
        e = self.op_entry(oid="walks-page")
        self.assertEqual(e["outcome"], "agree", e)
        self.assertIn("1 of 2 rows", next(c for c in e["cuts"] if c["way"] == "W1")["detail"])

    def test_vr6_a_sequence_cut_by_its_own_form(self):
        whole = {"id": "seq-all", "tool": "graphlet_sequence", "args": {"arm": "left", "walk": 0, "start": 0},
                 "paging": "all", "ways": ["W1", "W2", "W3", "W4"]}
        cut_op = {"id": "seq-page", "tool": "graphlet_sequence", "args": {"arm": "left", "walk": 0, "start": 0},
                  "paging": "none", "cut_of": "seq-all", "ways": ["W1", "W2", "W3", "W4"]}
        self.make(ops=[whole, cut_op])
        uncut = {"handle": "g_0123456789ab", "arm": "left", "orientation": "natural", "with_seed": False, "from": 0,
                 "length": 20, "ceiling_bytes": 16384, "evidence": {}, "to": 20, "sequence": "ACGTACGTACGTACGTACGT"}
        self.five(op=whole, answer=uncut, w5=False)
        page = dict(uncut, sequence="ACGTACGTAC", to=10, truncated=True, next_from=10)
        self.five(op=cut_op, answer=page, w5=False)
        e = self.op_entry(oid="seq-page")
        self.assertEqual(e["outcome"], "agree", e)
        # a slice that is no leading part, or a to past the uncut one, is no cut
        self.run_.cleanup()
        self.make(ops=[whole, cut_op])
        self.five(op=whole, answer=uncut, w5=False)
        self.five(op=cut_op, answer=dict(page, sequence="ACGTTCGTAC"), w5=False)
        c = next(c for c in self.op_entry(oid="seq-page")["cuts"] if c["way"] == "W1")
        self.assertEqual((c["outcome"], c["pointer"]), ("disagree", "/sequence"))
        self.assertEqual(cg.sequence_cut_problem(dict(page, to=25), uncut)[0], "/to")
        self.assertEqual(cg.sequence_cut_problem(dict(page, truncated=False), uncut)[0], "/truncated")

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
        self.put_files(exp, ">left_0 length_bp=80\nACGT\n",
                       answer={"path": "/x/g_0123456789ab.fa", "bytes": 100, "records": 4},
                       w3={"path": None, "bytes": 99, "records": 4, "export_id": "e1", "link": {"url": "u"},
                           "format": "fasta", "handle": UUID, "status": "done", "in_service": True})
        self.five(op=lst, answer={"total": 1, "rows": [{"handle": "g_1"}]}, ways=("W1", "W2", "W4"), w5=False)
        res = self.result()
        e = self.op_entry(res, "fasta")
        # path None vs a path: path and bytes are left out on every way; the files are compared by content
        self.assertEqual(e["outcome"], "agree", e)
        self.assertTrue(self.op_entry(res, "list")["smoke"])
        self.assertEqual(res["summary"]["ops"]["agree"], 1)   # the smoke is outside the five-way counts
        self.assertEqual(res["summary"]["per_tool"]["graphlet_list"]["smoke"], 1)

    def put_files(self, op, text, answer, w3=None, texts=None, ways=("W1", "W2", "W3", "W4"), raw_gz=False):
        """Records an export op on |ways| with its file kept (FORMAT v2 B3): |text| on every way unless |texts|."""
        for way in ways:
            raw = "%s.%s.file%s" % (SID, op["id"], ".gz" if raw_gz else "")
            f = self.run_.write_file(way, raw, (texts or {}).get(way, text))
            ans = copy.deepcopy(w3 if way == "W3" and w3 is not None else answer)
            self.run_.record(way, SID, op["id"], op["tool"], op["args"], answer=ans, file=f)

    def test_files_are_compared_by_content(self):
        exp = {"id": "fasta", "tool": "graphlet_export", "args": {"format": "fasta", "arm": "left"},
               "ways": ["W1", "W2", "W3", "W4"]}
        self.make(ops=[exp])
        rec = {"bytes": 10, "records": 1}
        self.put_files(exp, ">left_0\nACGT\n", rec, texts={"W3": ">left_0\nACGA\n"}, raw_gz=True)
        e = self.op_entry(oid="fasta")
        self.assertEqual(e["disagreeing"], ["W3-W4"])
        self.assertEqual(next(p for p in e["pairs"] if p["pair"] == "W3-W4")["pointer"], "/file/sha256")
        # a file not kept where the record says: not compared, the op inconclusive
        os.remove(os.path.join(self.run_.dir, "arms", "dr-lib", "raw", "%s.fasta.file.gz" % SID))
        e = self.op_entry(oid="fasta")
        self.assertEqual(e["outcome"], "disagree")   # W3-W4 still differs
        p = next(p for p in e["pairs"] if p["pair"] == "W1-W2")
        self.assertEqual((p["outcome"], p["kind"]), ("not_compared", "no_file"))
        self.assertIn("W2's file is not kept at arms/dr-lib/raw/", p["detail"])

    def test_json_and_mgt_files(self):
        whole = {"id": "json", "tool": "graphlet_export", "args": {"format": "json"}, "ways": ["W1", "W2", "W3", "W4"]}
        walk = {"id": "jwalk", "tool": "graphlet_export", "args": {"format": "json", "what": "walk", "arm": "left",
                                                                    "walks": [0]}, "ways": ["W1", "W2", "W3", "W4"]}
        mgt = {"id": "mgt", "tool": "graphlet_export", "args": {"format": "mgt"}, "ways": ["W1", "W2", "W3", "W4"]}
        self.make(ops=[whole, walk, mgt])
        doc = {"seed": {"length_bp": 30}, "arms": {"left": {"paths": [1]}}, "timing": {"elapsed_ms": 3},
               "annotation": {"labels": 4, "direct_reads": 10, "keys_mapped": 7}}
        other = copy.deepcopy(doc)
        other["timing"] = {"elapsed_ms": 9}
        other["annotation"].update(direct_reads=12, keys_mapped=9)   # lookahead work: per request (FORMAT v2 E1)
        self.put_files(whole, json.dumps(doc), {"records": 1}, texts={"W2": json.dumps(other)})
        wdoc = {"handle": "g_0123456789ab", "arm": "left", "walks": [{"walk": 0, "sequence": "ACGT"}],
                "evidence": {"view": None}}
        self.put_files(walk, json.dumps(wdoc), {"records": 1},
                       texts={"W3": json.dumps(dict(wdoc, handle=UUID))})
        f = self.run_.fixture
        # the library's dump carries its envelope; the service's stored body none, at another batch position
        self.put_files(mgt, f.doc(0), {"records": 9}, texts={"W3": f.doc(5, envelope=False)})
        res = self.result()
        for oid in ("json", "jwalk", "mgt"):
            self.assertEqual(self.op_entry(res, oid)["outcome"], "agree", (oid, self.op_entry(res, oid)))
        other["annotation"]["labels"] = 5
        self.put_files(whole, json.dumps(doc), {"records": 1}, texts={"W2": json.dumps(other)})
        p = next(p for p in self.op_entry(oid="json")["pairs"] if p["pair"] == "W1-W2")
        self.assertEqual(p["pointer"], "/file/json/annotation/labels")

    def test_w5_reads_the_file_through_export_view(self):
        exp = {"id": "fasta", "tool": "graphlet_export", "args": {"format": "fasta", "arm": "left"}}
        self.make(ops=[exp])
        self.oracle = cg.graphlet_oracle = gf.OracleStub()
        self.oracle.export_view = lambda fmt, receipt, text, what: {"records": receipt.get("records"),
                                                                    "fasta": [text.strip()], "fmt": fmt}
        self.put_files(exp, ">left_0\nACGT\n", {"records": 1})
        self.run_.record("W5", SID, "fasta", "graphlet_export", exp["args"], answer={"fasta": [">left_0\nACGT"]},
                         oracle={"checked": {"/fasta": "DESIGN §5 to_fasta"}})
        e = self.op_entry(oid="fasta")
        self.assertEqual(e["outcome"], "agree", e)
        views = [v for _, v in self.oracle.calls]
        self.assertEqual(views[0], {"records": 1, "fasta": [">left_0\nACGT"], "fmt": "fasta"})

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
        # VC3: without --reuse the comparison is recomputed every time, and never written
        os.remove(os.path.join(self.run_.dir, cg.COMPARE_FILE))
        cg.load_or_compute(self.run_.dir)
        self.assertFalse(os.path.exists(os.path.join(self.run_.dir, cg.COMPARE_FILE)))
        cg.compute(self.run_.dir)
        written = gp.sc.read_json(os.path.join(self.run_.dir, cg.COMPARE_FILE))
        self.assertEqual(cg.load_or_compute(self.run_.dir, reuse=True)["generated_utc"], written["generated_utc"])
        self.run_.record("W2", SID, "walks-len", "graphlet_walks", self.op["args"], answer=dict(gf.walks_answer(),
                                                                                                total=7))
        self.assertEqual(cg.stale(self.run_.dir, written), ["arms/dr-lib/graphlet_ops.jsonl"])
        fresh = cg.load_or_compute(self.run_.dir, reuse=True)
        self.assertEqual(fresh["scenarios"][SID]["ops"]["walks-len"]["outcome"], "disagree")

    def test_every_input_is_in_the_reuse_key(self):
        self.make()
        self.five()
        res = cg.compute(self.run_.dir)
        # the twin and the job record (VERIFY-C VC3's probe S1), a document, the identity records
        os.remove(os.path.join(self.run_.dir, "arms", "direct-local", "docs", SID + ".full.json.gz"))
        self.assertEqual(cg.stale(self.run_.dir, res), ["docs"])
        self.run_.write_echo(SID, {"direction": "left"})
        self.run_.write_mcp_runs(library_commit="c" * 40)
        self.assertEqual(cg.stale(self.run_.dir, res), ["arms/mcp/graphlet_jobs.jsonl", "arms/mcp/mcp_runs.jsonl",
                                                         "docs"])
        # the by-design table, the code
        changed = dict(res, inputs=dict(cg.input_digests(self.run_.dir), by_design_table="0" * 64, code="0" * 64))
        self.assertEqual(cg.stale(self.run_.dir, changed), ["by_design_table", "code"])
        # an export file the records name
        exp = {"id": "fasta", "tool": "graphlet_export", "args": {"format": "fasta"}}
        f = self.run_.write_file("W1", "x.fa", ">a\nAC\n")
        self.run_.record("W1", SID, "fasta", "graphlet_export", exp["args"], answer={"records": 1}, file=f)
        before = cg.input_digests(self.run_.dir)
        self.run_.write_file("W1", "x.fa", ">a\nAG\n")
        after = cg.input_digests(self.run_.dir)
        self.assertNotEqual(before["files"], after["files"])
        self.assertEqual(before["code"], after["code"])

    def test_unplanned_records_are_counted(self):
        self.make()
        self.five()
        self.run_.record("W1", SID, "not-in-plan", "graphlet_walks", self.op["args"], answer={})
        res = self.result()
        self.assertEqual(res["summary"]["unplanned_records"], 1)


if __name__ == "__main__":
    unittest.main()
