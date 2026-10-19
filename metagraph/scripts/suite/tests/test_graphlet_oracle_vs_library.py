"""The offline layer of manifest v8 (DESIGN-v8 §7 item 6): the library's GraphletTools and W5 (graphlet_oracle.py)
on the same engine documents with the same arguments agree on every pointer W5 checks.

Documents: every engine pair of api/python/tests/data/traverse that test_graphlet_oracle.py uses (the `detail:
graphlet` response, W5's `detail: full` twin), entered as W1 enters them (DESIGN-v8 §2): the library's own
standalone_text, then graphlet_load. The tools run as the service configures them (DESIGN-v8 §3: max_bytes 2048,
sequence pages 16 KiB, no local limits); a list is assembled over its pages (FORMAT-v8 §3: rows concatenated,
every other field from the first page), an export read from its file (graphlet_oracle.export_view). The compare
properties run on the library's compare answers: self-compare on every document, the a/b swap and the walk-set
differences on the documents that share a seed, prefix_subset on one locus at two radii.

A disagreement is a finding to analyse against the SPEC/DESIGN clause W5 names (lane O's notes, F1..F6), never a
reason to bend W5. This test imports the library; graphlet_oracle.py never does."""
import collections
import itertools
import json
import os
import re
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.dirname(HERE)
sys.path.insert(0, SUITE)
sys.path.insert(0, HERE)

import graphlet_oracle as O  # noqa: E402
from test_graphlet_oracle import DATA, NOT_TWINS, engine_pairs  # noqa: E402

LIB = os.path.normpath(os.path.join(SUITE, "..", "..", "api", "python"))
WIDE = 262144            # the service's largest page (DESIGN-v8 §5: whole-list calls fit uncut there)
LIST_TOOLS = ("graphlet_summary", "graphlet_walks", "graphlet_walk", "graphlet_labels", "graphlet_splits",
              "graphlet_claims", "graphlet_support")
SAME_SEED = [
    ["real/mini_refseq/mini_rep_00__trace_coords#0", "real/mini_refseq/mini_rep_00__trace_coords_column#0",
     "coords/rep0_column#0", "coords/rep0_cut1#0", "coords/rep0_kmer_null#0", "coords/rep0_lower_bound#0",
     "coords/rep0_mixed#0", "coords/rep0_unlimited_left#0"],
    ["real/mini_refseq/mini_ndm1__annotate_exh#0", "real/mini_refseq/mini_ndm1__trace#0"],
    ["cli_bubbles_r60", "cli_bubbles_r100", "cli_bubbles_steps150"], ["cli_quorum_tuned", "cli_quorum_exhaustive"]]
RADII = [("cli_bubbles_r60", "cli_bubbles_r100"), ("cli_bubbles_r60", "cli_bubbles_steps150")]
# library defects W5 finds, fixed on gr/labeled-traversal by a separate lane (FORMAT v2 G): each must still reproduce
# here until the fixed library is vendored, then its entry is removed (the test says so)
# (what, the calls it affects, the differences it excuses on them: nothing else on those calls is excused, R2-6)
# G1-G3 and VR4 landed with gr/v8-libfix (d39b6b38b): none is known now; a new one goes here with its call and
# differences, and the test says when the fixed library arrives
KNOWN_LIBRARY_DEFECTS = {}


def _presence_sets(d, side, p):
    """G4's reading, the library's presence profile: a run supports a base it covers where its own route (the merge
    partitions) passes the walk's segment. Built here from W5's own route machinery, for the known-defect check only
    (W5 itself follows DESIGN §5.1: _support_sets)."""
    seg_at = []
    for sid in p["segments"]:
        seg_at += [sid] * d.segment(side, sid)["length_bp"]
    sets = [set() for _ in range(p["length_bp"])]
    for i, run in enumerate(d.runs(side)):
        if run["from_bp"] == run["to_bp"]:
            continue
        route = d.route(side, i)
        if route is None:
            return None
        on = set(route)
        for b in range(run["from_bp"], min(run["to_bp"], p["length_bp"])):
            if seg_at[b] in on:
                sets[b].add(run["label"])
    return sets


# library defects of a READING (FORMAT v2 G): on the calls where the defect's reading and W5's differ, the library must
# equal the W5 record built under the defect's reading (every other value as W5 checks it) and differ from W5's own,
# so that exactly the defect's differences are excused (VERIFY-O-B VO-B1); elsewhere the library must agree with W5
KNOWN_READING_DEFECTS = {
    "G4": ("graphlet_support's displayed profile counts a run where its route passes the walk (node presence); DESIGN "
           "§5.1 derives each run's evidence at its anchor and keeps it where that chain is the walk",
           "graphlet_support", "_support_sets", _presence_sets),
}


def defect_alternative(tool, op, full, text, ctx):
    """(defect id, the W5 record under that defect's reading) for a call a reading defect affects, else (None, None)."""
    for k, (_, t, name, alt) in KNOWN_READING_DEFECTS.items():
        if t != tool:
            continue
        own = getattr(O, name)
        d = O._Doc(full, text, ctx.get("seed_index", 0))
        side = op.get("arm") if op.get("arm") in d.arms else (d.sides[0] if len(d.sides) == 1 else None)
        w = op.get("walk")
        if side is None or not isinstance(w, int) or isinstance(w, bool) or not 0 <= w < len(d.arms[side]["paths"]) \
                or d.mode != "constrain" or not d.aligned(side):
            continue
        path = d.arms[side]["paths"][w]
        if own(d, side, path) == alt(d, side, path):
            continue
        setattr(O, name, alt)
        try:
            return k, O.answer(tool, op, full, text, ctx)
        finally:
            setattr(O, name, own)
    return None, None


def excused(tool, args, diffs):
    """(the known defects that excuse some of |diffs|, the differences no known defect excuses)."""
    hit, rest = set(), []
    for d in diffs:
        ks = [k for k, (_, call, diff) in KNOWN_LIBRARY_DEFECTS.items() if call(tool, args) and diff(d)]
        if ks:
            hit.update(ks)
        else:
            rest.append(d)
    return hit, rest


def _calls(res):
    """The op plan of one document: every distinct behaviour per arm (INVENTORY §3.1), arguments named in full."""
    sides = [s for s in ("left", "right") if s in res["arms"]]
    out = [("graphlet_summary", {"max_bytes": WIDE}), ("graphlet_summary", {"detail": "limitations", "max_bytes": WIDE}),
           ("graphlet_labels", {"rank": "direct_bp", "n": 1000, "max_bytes": WIDE}),
           ("graphlet_labels", {"rank": "reach_bp", "n": 1000, "max_bytes": WIDE}),
           ("graphlet_labels", {"rank": "direct_bp", "n": 20}),           # one paging pass at the default 2048
           ("graphlet_export", {"format": "fasta"}), ("graphlet_export", {"format": "gfa"}),
           ("graphlet_export", {"format": "json"}), ("graphlet_export", {"format": "mgt"}), ("graphlet_save", {}),
           # VO1, VO14: arguments the tools refuse or answer by their own rules
           ("graphlet_labels", {"n": 0}), ("graphlet_labels", {"min_direct_bp": -1}),
           ("graphlet_labels", {"coordinates": True}), ("graphlet_labels", {"n": True}),
           ("graphlet_export", {"format": "fasta", "walks": [0]}), ("graphlet_summary", {"detail": "nope"}),
           ("graphlet_export", {"format": "fasta", "what": "walk"}), ("graphlet_subtrie", {"labels": ["no-such"]}),
           # fitting errors at small ceilings (INVENTORY §1.1 result_too_large): not_comparable against W5 (F3)
           ("graphlet_summary", {"max_bytes": 512}), ("graphlet_summary", {"detail": "limitations", "max_bytes": 300})]
    # phase B, graphlet_labels name= (INVENTORY §3.1): 3 labels by name, by ref and as a bare string, per arm and
    # both; coordinates=true with name; an unknown label
    for k, lab in enumerate(sorted({0, len(res.get("label_dict") or [1]) - 1, len(res.get("label_dict") or [1]) // 2})):
        e = (res.get("label_dict") or [{}])[lab] if res.get("label_dict") else {}
        ref = "h:%s:%s" % (e.get("column"), e.get("seq_id")) if e.get("kind") == "header" else "c:%s" % e.get("column")
        sel = [{"id": lab}, {"ref": ref}, {"name": e.get("name")}][k % 3] if e else {"id": lab}
        out.append(("graphlet_labels", {"name": sel, "max_bytes": WIDE}))
        out += [("graphlet_labels", {"name": sel, "arm": s, "max_bytes": WIDE}) for s in sides]
    if res.get("label_dict"):
        out += [("graphlet_labels", {"name": res["label_dict"][0]["name"], "max_bytes": WIDE}),
                ("graphlet_labels", {"name": {"id": 0}, "coordinates": True, "max_bytes": WIDE}),
                ("graphlet_labels", {"name": "no-such-label"}), ("graphlet_labels", {"name": {"id": 0}})]
        names = [e.get("name") for e in res["label_dict"]]
        out += [("graphlet_labels", {"name": x, "max_bytes": WIDE}) for x in sorted({x for x in names
                                                                                       if names.count(x) > 1})]
    if len(sides) == 1:   # an arm this retrieval did not walk
        other = "left" if sides == ["right"] else "right"
        out += [("graphlet_walks", {"arm": other}), ("graphlet_summary", {"arm": other})]
    nl = len(res.get("label_dict") or [])
    out.append(("graphlet_subtrie", {"labels": [{"id": 0}]}))
    if nl > 1:
        out.append(("graphlet_subtrie", {"labels": [{"id": 0}, {"id": nl - 1}], "mode": "all"}))
    for s in sides:
        a = res["arms"][s]
        n = len(a["paths"])
        lens = sorted(p["length_bp"] for p in a["paths"])
        out += [("graphlet_summary", {"arm": s, "max_bytes": WIDE}),
                ("graphlet_summary", {"arm": s, "detail": "limitations", "max_bytes": WIDE})]
        for rank in ("support", "length", "loss", "id"):
            out.append(("graphlet_walks", {"arm": s, "rank": rank, "n": 1000, "max_bytes": WIDE}))
        out += [("graphlet_walks", {"arm": s, "rank": "id", "n": 10}),   # byte-bound pages at 2048
                ("graphlet_walks", {"arm": s, "rank": "support", "n": 1000, "spell": "full", "max_bytes": WIDE}),
                ("graphlet_walks", {"arm": s, "rank": "length", "n": 1000, "spell": "tail", "tail_bp": 60,
                                    "max_bytes": WIDE}),
                ("graphlet_walks", {"arm": s, "rank": "length", "n": 1000, "spell": "tail", "tail_bp": 100000,
                                    "max_bytes": WIDE}),
                ("graphlet_walks", {"arm": s, "rank": "id", "n": 1000, "min_bp": lens[len(lens) // 2], "max_bytes": WIDE}),
                ("graphlet_walks", {"arm": s, "rank": "id", "n": 1000, "min_bp": lens[-1] + 1, "max_bytes": WIDE})]
        for lab in sorted({e["label"] for p in a["paths"] for e in p.get("end_labels") or []})[:3]:
            for rc in (True, False):
                out.append(("graphlet_walks", {"arm": s, "rank": "id", "n": 1000, "label": {"id": lab},
                                               "route_consistent": rc, "max_bytes": WIDE}))
        for w in (range(n) if n <= 12 else sorted({0, n - 1, n // 2} | {p["id"] for p in a["paths"]
                                                                          if p.get("continuation")})):
            out.append(("graphlet_walk", {"arm": s, "walk": w, "max_bytes": WIDE}))
        out += [("graphlet_walk", {"arm": s, "walk": 0}), ("graphlet_walk", {"arm": s, "walk": n}),
                ("graphlet_walk", {"arm": s, "walk": -1}),
                ("graphlet_walks", {"arm": s, "coordinates": True, "n": 1000, "max_bytes": WIDE}),
                ("graphlet_walks", {"arm": s, "rank": "id", "n": 1000, "spell": "tail", "tail_bp": 0, "max_bytes": WIDE}),
                ("graphlet_walks", {"arm": s, "rank": "id", "n": 1000, "spell": "tail", "tail_bp": 1, "max_bytes": WIDE}),
                ("graphlet_walks", {"arm": s, "n": 0}), ("graphlet_walks", {"arm": s, "min_bp": True}),
                ("graphlet_walks", {"arm": s, "route_consistent": "yes"}),
                ("graphlet_walks", {"arm": s, "label": "no-such-label"}), ("graphlet_splits", {"arm": s, "n": 0}),
                ("graphlet_export", {"format": "json", "what": "walk", "arm": s}),
                ("graphlet_sequence", {"arm": s, "walk": 0, "segment": 0}), ("graphlet_sequence", {"arm": s}),
                ("graphlet_sequence", {"arm": s, "walk": 0, "start": 5, "end": 2}),
                ("graphlet_sequence", {"arm": s, "segment": 0, "with_seed": True}),
                # R2-4: argument types as the tool checks them
                ("graphlet_sequence", {"arm": s, "walk": 0, "with_seed": "yes"}),
                ("graphlet_sequence", {"arm": s, "walk": True}), ("graphlet_sequence", {"arm": s, "segment": False}),
                ("graphlet_sequence", {"arm": s, "walk": 0, "start": True}),
                ("graphlet_sequence", {"arm": s, "walk": 0, "end": True}),
                ("graphlet_sequence", {"arm": s, "walk": -1})]
        out += [("graphlet_labels", {"arm": s, "rank": "direct_bp", "n": 1000, "max_bytes": WIDE}),
                ("graphlet_labels", {"arm": s, "rank": "reach_bp", "n": 1000, "min_direct_bp": 1, "max_bytes": WIDE})]
        for at in (0, lens[-1] // 2, lens[-1] + 5):
            out.append(("graphlet_labels", {"arm": s, "rank": "direct_bp", "n": 1000, "at_bp": at, "max_bytes": WIDE}))
        for m in (0, 2, 10 ** 6):
            out.append(("graphlet_splits", {"arm": s, "n": 1000, "min_labels_before": m, "max_bytes": WIDE}))
        for w in sorted({0, n - 1}):
            for o in ("natural", "walk"):
                for ws in (False, True):
                    out.append(("graphlet_sequence", {"arm": s, "walk": w, "orientation": o, "with_seed": ws}))
            out += [("graphlet_sequence", {"arm": s, "walk": w, "start": 5, "end": 25}),
                    ("graphlet_sequence", {"arm": s, "walk": w, "start": 3, "end": 10 ** 7})]
        for o in ("natural", "walk"):
            out.append(("graphlet_sequence", {"arm": s, "segment": len(a["segments"]) - 1, "orientation": o}))
        out += [("graphlet_export", {"format": "fasta", "arm": s}),
                ("graphlet_export", {"format": "fasta", "arm": s, "walks": [0]}),
                ("graphlet_export", {"format": "json", "what": "walk", "arm": s, "walks": list(range(min(n, 3)))}),
                ("graphlet_subtrie", {"labels": [{"id": 0}], "arm": s})]
        for p in a["paths"]:
            if p.get("continuation"):
                out.append(("traverse_continue", {"arm": s, "walk": p["id"], "execute": False, "max_bytes": WIDE}))
        conts = [p["id"] for p in a["paths"] if p.get("continuation")][:2]
        for w in conts:
            out += [("traverse_continue", {"arm": s, "walk": w, "execute": False, "max_bytes": WIDE,
                                           "overrides": {"labels": {"loss_budget": 3.0}}}),
                    ("traverse_continue", {"arm": s, "walk": w, "execute": False, "max_bytes": WIDE,
                                           "overrides": {"bounds": {"max_extension_bp": 50}}}),
                    ("traverse_continue", {"arm": s, "walk": w, "execute": False, "max_bytes": WIDE,
                                           "reset_branches": True})]
        out += [("traverse_continue", {"arm": s, "walk": p["id"], "execute": False, "max_bytes": WIDE})
                for p in a["paths"][:1] if not p.get("continuation")]
        # phase B, graphlet_claims (INVENTORY §3.1): min_bp {0, mid} x route_consistent; labels one and a list; a
        # paging pass at 2048; annotate and trace retrievals are among the documents
        runs = a.get("runs") or []
        labs = sorted({r["label"] for r in runs})[:3] or [0]
        mid = sorted(r["to_bp"] - r["from_bp"] for r in runs)[len(runs) // 2] if runs else 10
        for rc, mb in itertools.product((True, False), (0, mid)):
            out.append(("graphlet_claims", {"arm": s, "route_consistent": rc, "min_bp": mb, "max_bytes": WIDE}))
        out += [("graphlet_claims", {"arm": s, "labels": {"id": labs[0]}, "route_consistent": False, "max_bytes": WIDE}),
                ("graphlet_claims", {"arm": s, "labels": [{"id": x} for x in labs], "max_bytes": WIDE}),
                ("graphlet_claims", {"arm": s}), ("graphlet_claims", {"arm": s, "labels": []}),
                ("graphlet_claims", {"arm": s, "min_bp": -1})]
        # phase B, graphlet_support (INVENTORY §3.1): the walk sample of graphlet_walk, step None and step = k; a
        # paging pass at 2048; a walk out of range
        for w in (range(n) if n <= 12 else sorted({0, n - 1, n // 2} | {p["id"] for p in a["paths"]
                                                                          if p.get("continuation")})):
            out += [("graphlet_support", {"arm": s, "walk": w, "max_bytes": WIDE}),
                    ("graphlet_support", {"arm": s, "walk": w, "step": 20, "max_bytes": WIDE})]
        out += [("graphlet_support", {"arm": s, "walk": 0}), ("graphlet_support", {"arm": s, "walk": n}),
                ("graphlet_support", {"arm": s, "walk": 0, "step": -1})]
    return out


class OracleVsLibrary(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        sys.path.insert(0, LIB)
        try:
            from metagraph.traverse import parser
            from metagraph.traverse.mcp_tools import GraphletTools
            from metagraph.traverse.store import GraphletStore
        except ImportError as e:   # pragma: no cover - the suite's vendored copy runs without the library
            raise unittest.SkipTest("the library is not importable here: %s" % e)
        cls.tmp = tempfile.TemporaryDirectory()
        cls.exp = os.path.join(cls.tmp.name, "exports")
        os.makedirs(cls.exp)
        cls.tools = GraphletTools(GraphletStore(os.path.join(cls.tmp.name, "spool")), {}, export_dir=cls.exp)
        cls.docs = {}       # name -> (handle, full answer or None, seed index, the engine's MGT body)
        cls.pairs = engine_pairs()
        for cell, i, g, f in cls.pairs:
            # the library loads W1's document (its own standalone_text); W5 reads the engine's body (VO12), so that
            # a rewrite the conversion made would show
            text = parser.standalone_text(g["results"][i]["graphlet"], g["results"][i], g)
            name = "%s#%d" % (cell, i)
            cls.docs[name] = (cls._load(name, text), f, i, g["results"][i]["graphlet"])
        for n in ("cli_bubbles_r60", "cli_bubbles_r100", "cli_bubbles_steps150", "cli_quorum_tuned",
                  "cli_quorum_exhaustive"):
            # body-only CLI documents without a manifest: the library's own tests stamp an index_fp into H
            with open(os.path.join(DATA, "documents", "compare", n + ".mgt"), encoding="utf-8") as fh:
                h, rest = fh.read().split("\n", 1)
            f = h.split(" ")
            f[14] = "ab" * 32
            text = " ".join(f) + "\n" + rest
            cls.docs[n] = (cls._load(n, text), None, 0, text)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    @classmethod
    def _load(cls, name, text):
        path = re.sub(r"[^A-Za-z0-9_.-]", "_", name) + ".mgt"
        with open(os.path.join(cls.exp, path), "w", encoding="utf-8") as fh:
            fh.write(text)
        out = cls.tools.graphlet_load(path)
        assert "handle" in out, out
        return out["handle"]

    def assemble(self, tool, kw):
        """The assembled answer (FORMAT-v8 §3): a list tool's rows over every page, a sequence over next_from."""
        out = getattr(self.tools, tool)(**kw)
        if "error" in out:
            return out
        if tool in LIST_TOOLS or tool == "graphlet_compare":
            def page_rows(page):
                # FORMAT v2 A3: each row of a page with row_truncated carries that page's cut_fields as _cut_fields
                # (lane L's graphlet_paging.assemble; this test's stand-in until it lands on this branch)
                rows_ = [dict(r) for r in page.get("rows") or []]
                if page.get("row_truncated"):
                    for r in rows_:
                        r["_cut_fields"] = list(page.get("cut_fields") or [])
                return rows_
            rows = page_rows(out)
            cur = out.pop("next_cursor", None)
            while cur:
                nxt = getattr(self.tools, tool)(**dict(kw, cursor=cur))
                self.assertNotIn("error", nxt, (tool, kw))
                rows += page_rows(nxt)
                cur = nxt.get("next_cursor")
            if "rows" in out:
                out["rows"] = rows
        elif tool == "graphlet_sequence":
            while out.get("next_from") is not None:
                nxt = self.tools.graphlet_sequence(**dict(kw, start=out["next_from"]))
                out["sequence"] += nxt["sequence"]
                out["next_from"] = nxt.get("next_from")
        return out

    def test_every_checked_pointer_agrees(self):
        groups = collections.defaultdict(list)
        outcomes = collections.Counter()
        defects = {}
        reading, reading_gone = collections.Counter(), collections.Counter()
        n_calls = n_checked = 0
        declared = collections.Counter()
        for name, (h, f, i, text) in self.docs.items():
            if f is None:
                continue
            for tool, args in _calls(f["results"][i]):
                kw = dict(args, handle=h)
                if tool == "graphlet_save" and "path" not in args:
                    kw["path"] = "saved_%d.mgt" % n_calls
                lib = self.assemble(tool, kw)
                if tool in ("graphlet_export", "graphlet_save") and "error" not in lib:
                    with open(lib["path"], encoding="utf-8") as fh:
                        body = fh.read()
                    lib = O.export_view("save" if tool == "graphlet_save" else args.get("format", "fasta"), lib, body,
                                        args.get("what"))
                op = {k: v for k, v in args.items() if k not in ("max_bytes",)}
                rec = O.answer(tool, op, f, text, {"seed_index": i, "max_bytes": args.get("max_bytes", 2048)})
                n_calls += 1
                n_checked += len(rec["oracle"]["checked"])
                if rec["oracle"].get("declared"):
                    declared[rec["oracle"]["declared"]] += 1
                    if rec["oracle"]["declared"] == "twin_mismatch":
                        self.assertIn(name.split("#")[0], NOT_TWINS)
                verdict = O.compare(rec, lib)
                gk, galt = defect_alternative(tool, op, f, text, {"seed_index": i, "max_bytes": args.get("max_bytes", 2048)})
                if gk is not None:
                    alt_v = O.compare(galt, lib)
                    if verdict["outcome"] == "disagree" and alt_v["outcome"] == "agree":
                        reading[gk] += 1               # exactly the defect's difference: excused
                        outcomes["excused " + gk] += 1
                        continue
                    if verdict["outcome"] == "agree" and alt_v["outcome"] == "disagree":
                        reading_gone[gk] += 1          # the library follows DESIGN on this affected call
                    elif verdict["outcome"] == "disagree":
                        for d in alt_v["differences"]:   # a difference the defect does not explain
                            groups[(tool + " beyond " + gk, re.sub(r"/\d+", "/*", d["pointer"]))].append(
                                (name, args, d))
                        continue
                outcomes[verdict["outcome"]] += 1
                if verdict["outcome"] == "inconclusive":    # R2-5: only what W5 cannot judge, here the non-twins
                    self.assertIn(name.split("#")[0], NOT_TWINS, (name, tool, args, verdict["detail"]))
                known, diffs = excused(tool, op, verdict["differences"])
                for k in known:
                    defects[k] = defects.get(k, 0) + 1
                for d in diffs:
                    groups[(tool, re.sub(r"/\d+", "/*", d["pointer"]))].append((name, args, d))
        report = []
        for (tool, ptr), items in sorted(groups.items()):
            report.append("%s %s: %d in %s" % (tool, ptr, len(items), sorted({x[0] for x in items})[:4]))
            for name, args, d in items[:2]:
                report.append("    %s %s\n      W5  %s\n      lib %s" % (
                    name, json.dumps(args)[:120], json.dumps(d["oracle"])[:240],
                    "MISSING" if d.get("missing") else json.dumps(d.get("other", d.get("present")))[:240]))
        self.assertEqual(report, [], "W5 and the library disagree:\n" + "\n".join(report))
        for k, (what, _, _, _) in KNOWN_READING_DEFECTS.items():
            self.assertGreater(reading[k], 0, "known library defect %s no longer reproduces (%s): the fixed library is "
                                              "here, remove it from KNOWN_READING_DEFECTS" % (k, what))
            self.assertEqual(reading_gone[k], 0, "%s: %d affected calls where the library follows DESIGN: the fix is "
                                                 "partial" % (k, reading_gone[k]))
        print("\n[offline] known reading defects excused: %s" % dict(reading), file=sys.stderr)
        for k, (what, _, _) in KNOWN_LIBRARY_DEFECTS.items():
            self.assertGreater(defects.get(k, 0), 0, "known library defect %s no longer reproduces (%s): the fixed "
                                                    "library is here, remove it from KNOWN_LIBRARY_DEFECTS" % (k, what))
        self.assertGreater(n_calls, 3500)
        self.assertGreater(n_checked, 150000)
        self.assertEqual(declared["twin_mismatch"] > 0, True)
        self.assertGreater(outcomes["not_comparable"], 40)      # result_too_large at 512 and 300 bytes (F3)
        self.assertGreater(outcomes["agree"], 3000)

    def test_run_step_on_library_records(self):
        """FORMAT v2 F4 end to end: W1 records made from the library's answers (lane L's record shape: args_resolved
        with "<handle>", the assembled answer, paging_complete, library.limits; an error as {code, fields}; an export
        as export_view of its file), the run step over a run directory, and compare() of every W5 line with its W1
        record: no disagreement."""
        import gzip as gz
        with tempfile.TemporaryDirectory() as run:
            for d in ("arms/direct-local/raw", "arms/direct-local/docs", "arms/dl-lib"):
                os.makedirs(os.path.join(run, d))
            w1, sources, cell_text = [], {}, {}
            for cell, i, g, f in self.pairs:
                if cell in NOT_TWINS:
                    continue
                sid = "G-%s-%d" % (re.sub(r"[^A-Za-z0-9]", "_", cell), i)
                sources[sid] = (dict(f, results=[f["results"][i]]), i)
                cell_text[sid] = g["results"][i]["graphlet"]
                resp = dict(g, results=[g["results"][i]])
                with gz.open(os.path.join(run, "arms/direct-local/raw/%s.direct-local.0.json.gz" % sid), "wt",
                             encoding="utf-8") as fh:
                    json.dump({"answer": {"body": resp}}, fh)
                with gz.open(os.path.join(run, "arms/direct-local/docs/%s.full.json.gz" % sid), "wt",
                             encoding="utf-8") as fh:
                    json.dump(dict(f, results=[f["results"][i]]), fh)
                h = self.docs["%s#%d" % (cell, i)][0]
                for k, (tool, args) in enumerate(_calls(f["results"][i])):
                    kw = dict(args, handle=h)
                    if tool == "graphlet_save":
                        kw["path"] = "run_saved_%s_%d.mgt" % (sid, k)
                    lib = self.assemble(tool, kw)
                    rec = {"run_id": "R", "scenario": sid, "way": "W1", "op": "op%d" % k, "tool": tool,
                           "args_resolved": dict(args, handle="<handle>"), "selector_error": None,
                           "paging_complete": True, "answer": None, "error": None,
                           "library": {"limits": {"max_bytes": 2048, "sequence_max_bytes": 16384}}}
                    if "error" in lib:
                        rec["error"] = {"code": lib["error"], "fields": {}}
                    elif tool in ("graphlet_export", "graphlet_save"):
                        with open(lib["path"], encoding="utf-8") as fh:
                            rec["answer"] = O.export_view("save" if tool == "graphlet_save" else args["format"],
                                                          lib, fh.read(), args.get("what"))
                    else:
                        rec["answer"] = lib
                    w1.append(rec)
            with open(os.path.join(run, "arms/dl-lib/graphlet_ops.jsonl"), "w", encoding="utf-8") as fh:
                for rec in w1:
                    fh.write(json.dumps(rec) + "\n")
            path = O.write_run(run, quiet=True)
            with open(path, encoding="utf-8") as fh:
                w5 = [json.loads(x) for x in fh]
        self.assertEqual(len(w5), len(w1))
        outcomes, problems = collections.Counter(), []
        for a, b in zip(w1, w5):
            self.assertEqual((a["scenario"], a["op"], a["args_resolved"]), (b["scenario"], b["op"], b["args_resolved"]))
            v = O.compare(b, a)
            if v["outcome"] == "disagree" and a["tool"] in {t for _, t, _, _ in KNOWN_READING_DEFECTS.values()}:
                twin, _ = sources[a["scenario"]]     # the scenario's full twin, its one result (seed index 0)
                op = {x: y for x, y in a["args_resolved"].items() if x != "handle"}
                gk, galt = defect_alternative(a["tool"], op, twin, cell_text[a["scenario"]], {"seed_index": 0})
                if gk is not None and O.compare(galt, a)["outcome"] == "agree":
                    outcomes["excused " + gk] += 1
                    continue
            outcomes[v["outcome"]] += 1
            _, rest = excused(a["tool"], {x: y for x, y in a["args_resolved"].items() if x != "handle"},
                              v["differences"])
            if v["outcome"] == "disagree" and (rest or not v["differences"]):
                problems.append("%s %s %s: %s" % (a["scenario"], a["tool"], json.dumps(a["args_resolved"])[:100],
                                                  v["detail"][:200]))
        self.assertEqual(problems, [], "\n".join(problems[:10]))
        self.assertGreater(outcomes["agree"], 3000)
        self.assertGreater(outcomes["not_comparable"], 40)

    def test_known_defects_excuse_only_their_difference(self):
        """R2-6: a known defect excuses only its own difference on its own calls (the G2 entry as it stood before the
        fixed library landed, injected here so the mechanism stays tested)."""
        g2 = ("graphlet_walks spell=tail tail_bp=0 gives the whole right flank",
              lambda tool, args: tool == "graphlet_walks" and args.get("spell") == "tail"
              and args.get("tail_bp") == 0 and args.get("arm") == "right",
              lambda d: re.match(r"/rows/\d+/sequence\Z", d["pointer"]) is not None)
        KNOWN_LIBRARY_DEFECTS["G2"] = g2
        self.addCleanup(KNOWN_LIBRARY_DEFECTS.pop, "G2", None)
        args = {"arm": "right", "spell": "tail", "tail_bp": 0, "rank": "id"}
        diffs = [{"pointer": "/rows/0/sequence"}, {"pointer": "/rows/3/sequence"}, {"pointer": "/rows/0/n_full"},
                 {"pointer": "/total"}]
        hit, rest = excused("graphlet_walks", args, diffs)
        self.assertEqual((hit, [d["pointer"] for d in rest]), ({"G2"}, ["/rows/0/n_full", "/total"]))
        self.assertEqual(excused("graphlet_walks", dict(args, arm="left"), diffs[:1]), (set(), diffs[:1]))
        self.assertEqual(excused("graphlet_walks", dict(args, tail_bp=1), diffs[:1]), (set(), diffs[:1]))

    def test_g4_on_every_walk(self):
        """G4 on every walk of every twin (VERIFY-O-B VO-B1: 73 of the 399 non-empty constrain walks): where the
        library's presence profile and DESIGN §5.1's displayed profile differ, the library equals W5's record under
        the presence reading and differs from W5's own (it excuses exactly that difference); on every other walk the
        library agrees with W5."""
        counts = collections.Counter()
        problems = []
        for name, (h, f, i, text) in self.docs.items():
            if f is None or name.split("#")[0] in NOT_TWINS:
                continue
            for side, a in f["results"][i]["arms"].items():
                for p in a["paths"]:
                    op = {"arm": side, "walk": p["id"]}
                    lib = self.assemble("graphlet_support", dict(op, handle=h, max_bytes=WIDE))
                    own = O.compare(O.answer("graphlet_support", op, f, text, {"seed_index": i}), lib)["outcome"]
                    gk, galt = defect_alternative("graphlet_support", op, f, text, {"seed_index": i})
                    if gk is None:
                        counts["unaffected " + own] += 1
                        if own not in ("agree", "not_checked"):
                            problems.append("%s %s walk %d: %s on an unaffected walk" % (name, side, p["id"], own))
                        continue
                    alt = O.compare(galt, lib)["outcome"]
                    counts["affected"] += 1
                    if (own, alt) != ("disagree", "agree"):
                        problems.append("%s %s walk %d: W5 %s, under G4 %s" % (name, side, p["id"], own, alt))
        self.assertEqual(problems, [], "\n".join(problems[:10]))
        self.assertEqual(counts["affected"], 73, counts)

    def test_refused_arguments(self):
        """F-D2 against the library: paths outside the export directory, a budget without local limits, a cursor of
        another listing; W5 expects the tool's refusal (library-defined codes)."""
        name = next(k for k in self.docs if k.startswith("documents/merge#"))
        h, f, i, text = self.docs[name]
        outside = os.path.join(self.tmp.name, "outside")
        for tool, args in [("graphlet_export", {"format": "fasta", "path": outside + ".fa"}),
                           ("graphlet_export", {"format": "fasta", "path": "../escape.fa"}),
                           ("graphlet_save", {"path": outside + ".mgt"}),
                           ("graphlet_walks", {"arm": "right", "budget": {"work_units": 1000}}),
                           ("graphlet_summary", {"budget": {"work_units": 1000}})]:
            lib = getattr(self.tools, tool)(handle=h, **args)
            rec = O.answer(tool, args, f, text, {"seed_index": i})
            self.assertEqual(O.compare(rec, lib)["outcome"], "agree", (tool, args, lib))
            self.assertFalse(os.path.exists(outside + ".fa") or os.path.exists(outside + ".mgt"))
        # a cursor of one listing (walks, one row a page) sent with other arguments, or to another tool: bad_cursor
        name2 = next(k for k in self.docs if k.startswith("real/mini_refseq/mini_ndm1__trace#"))
        h2, f2, i2, text2 = self.docs[name2]
        page = self.tools.graphlet_walks(handle=h2, arm="right", rank="id", n=1, max_bytes=WIDE)
        self.assertIn("next_cursor", page)
        for lib_args in ({"rank": "length"}, {"rank": "id", "min_bp": 1}):
            lib = self.tools.graphlet_walks(handle=h2, arm="right", n=1, max_bytes=WIDE, cursor=page["next_cursor"],
                                            **lib_args)
            op = dict(lib_args, arm="right", n=1, cursor="<cursor>")
            rec = O.answer("graphlet_walks", op, f2, text2, {"seed_index": i2, "cursor_bound_elsewhere": True})
            self.assertEqual(O.compare(rec, lib)["outcome"], "agree", (lib_args, lib))
        # R2B-1, the verifier's four probes: the source op (walks, arm right, rank id, n 1) and its cursor reused;
        # the library resumes where the arguments differ only in form (arm r, spell none named, n), refuses rank length
        src = {"tool": "graphlet_walks", "args": {"handle": "<handle>", "arm": "right", "rank": "id", "n": 1}}
        for extra, want in [({"arm": "r"}, "not_checked"), ({"arm": "right", "spell": "none"}, "not_checked"),
                            ({"arm": "right", "n": 3}, "not_checked"), ({"arm": "right", "rank": "length"}, "agree")]:
            args = dict({"rank": "id", "n": 1}, **extra)
            lib = self.tools.graphlet_walks(handle=h2, max_bytes=WIDE, cursor=page["next_cursor"], **args)
            self.assertEqual("error" in lib, want == "agree", (extra, lib))
            rec = O.answer("graphlet_walks", dict(args, handle="<handle>", cursor="<cursor>"), f2, text2,
                           {"seed_index": i2, "cursor_source": src})
            self.assertEqual(O.compare(rec, lib)["outcome"], want, extra)
        lib = self.tools.graphlet_labels(handle=h2, cursor=page["next_cursor"], max_bytes=WIDE)
        rec = O.answer("graphlet_labels", {"cursor": "<cursor>"}, f2, text2, {"seed_index": i2,
                                                                            "cursor_bound_elsewhere": True})
        self.assertEqual(O.compare(rec, lib)["outcome"], "agree", lib)
        # the same listing resumed: a page W5 makes no claim about
        lib = self.tools.graphlet_walks(handle=h2, arm="right", rank="id", n=1, max_bytes=WIDE, cursor=page["next_cursor"])
        self.assertNotIn("error", lib)
        rec = O.answer("graphlet_walks", {"arm": "right", "rank": "id", "n": 1, "cursor": "<cursor>"}, f2, text2,
                       {"seed_index": i2})
        self.assertEqual(O.compare(rec, lib)["outcome"], "not_checked")

    def test_load_receipt_summary(self):
        """graphlet_load's summary, with a ceiling it fits under, is graphlet_summary's (W5 checks it only then)."""
        from metagraph.traverse.mcp_tools import GraphletTools
        wide = GraphletTools(self.tools.store, {}, export_dir=self.exp, max_bytes=16 * 1024 * 1024)
        n = 0
        for name, (h, f, i, text) in self.docs.items():
            if f is None or name.split("#")[0] in NOT_TWINS:
                continue
            path = re.sub(r"[^A-Za-z0-9_.-]", "_", name) + ".mgt"
            rec = O.answer("graphlet_load", {"path": path}, f, text, {"seed_index": i, "max_bytes": 16 * 1024 * 1024})
            self.assertEqual(O.differences(rec, wide.graphlet_load(path)), [], name)
            n += len(rec["oracle"]["checked"])
        self.assertGreater(n, 400)

    def test_sequence_pages(self):
        """A slice longer than the page (a 1 KiB ceiling here), assembled over next_from, is W5's whole slice."""
        from metagraph.traverse.mcp_tools import GraphletTools
        small = GraphletTools(self.tools.store, {}, export_dir=self.exp, sequence_max_bytes=1024)
        n = 0
        for name, (h, f, i, text) in self.docs.items():
            if f is None or name.split("#")[0] in NOT_TWINS:
                continue
            for side, a in f["results"][i]["arms"].items():
                p = max(a["paths"], key=lambda x: x["length_bp"])
                if p["length_bp"] < 1500 or "sequence" not in a["segments"][0]:
                    continue
                args = {"arm": side, "walk": p["id"], "with_seed": True}
                first = small.graphlet_sequence(handle=h, **args)
                self.assertIn("next_from", first)
                tools, self.tools = self.tools, small
                try:
                    lib = self.assemble("graphlet_sequence", dict(args, handle=h))
                finally:
                    self.tools = tools
                rec = O.answer("graphlet_sequence", args, f, text, {"seed_index": i, "sequence_max_bytes": 1024})
                self.assertEqual(O.differences(rec, lib), [], (name, side))
                n += 1
        self.assertGreater(n, 10)

    def test_compare_properties(self):
        def whole(a, b, mode, arm=None):
            kw = {"a": a, "b": b, "mode": mode, "max_bytes": WIDE}
            if arm:
                kw["arm"] = arm
            return self.assemble("graphlet_compare", kw)
        problems, seen = [], collections.Counter()
        for name, (h, f, i, text) in self.docs.items():
            fi = f if f is None else {"results": [f["results"][i]]}
            for mode in ("claims", "walks", "labels", "prefix_subset"):
                ans = whole(h, h, mode)
                seen["self %s" % ans.get("comparable")] += 1
                problems += ["self %s %s: %s" % (name, mode, p) for p in O.compare_self_properties(ans, fi, text)]
        for grp in SAME_SEED:
            for x, y in itertools.combinations(grp, 2):
                (ha, fa, ia, ta), (hb, fb, ib, tb) = self.docs[x], self.docs[y]
                fa = fa if fa is None else {"results": [fa["results"][ia]]}
                fb = fb if fb is None else {"results": [fb["results"][ib]]}
                arms = [None] + [s for s in ("left", "right") if s in O._Doc(fa, ta).arms and s in O._Doc(fb, tb).arms]
                for mode, arm in itertools.product(("claims", "walks", "labels"), arms):
                    ab, ba = whole(ha, hb, mode, arm), whole(hb, ha, mode, arm)
                    problems += ["swap %s %s %s %s: %s" % (x, y, mode, arm, p)
                                 for p in O.compare_swap_properties(ab, ba)]
                    if mode == "walks":
                        ps = O.compare_walks_properties(ab, fa, ta, fb, tb, arm=arm)
                        problems += ["walks %s %s %s: %s" % (x, y, arm, p) for p in ps]
                        if ab.get("comparable") in (True, "qualified") and ab.get("rows"):
                            seen["walk-set differences checked" + (" (one arm)" if arm else "")] += 1
        for (small, large), arm in itertools.product(RADII, (None, "right")):
            (ha, _, _, ta), (hb, _, _, tb) = self.docs[small], self.docs[large]
            ans = whole(ha, hb, "prefix_subset", arm)
            self.assertIn(ans.get("comparable"), (True, "qualified"))
            problems += ["prefix_subset %s %s %s: %s" % (small, large, arm, p)
                         for p in O.compare_prefix_subset_properties(ans, None, ta, None, tb, arm=arm)]
        self.assertEqual(problems, [], "\n".join(problems[:20]))
        self.assertGreater(seen["self True"], 30)
        self.assertGreater(seen["walk-set differences checked"], 20)
        self.assertGreater(seen["walk-set differences checked (one arm)"], 20)


if __name__ == "__main__":
    unittest.main()
