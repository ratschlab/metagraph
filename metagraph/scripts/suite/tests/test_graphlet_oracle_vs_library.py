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
LIST_TOOLS = ("graphlet_summary", "graphlet_walks", "graphlet_walk", "graphlet_labels", "graphlet_splits")
SAME_SEED = [
    ["real/mini_refseq/mini_rep_00__trace_coords#0", "real/mini_refseq/mini_rep_00__trace_coords_column#0",
     "coords/rep0_column#0", "coords/rep0_cut1#0", "coords/rep0_kmer_null#0", "coords/rep0_lower_bound#0",
     "coords/rep0_mixed#0", "coords/rep0_unlimited_left#0"],
    ["real/mini_refseq/mini_ndm1__annotate_exh#0", "real/mini_refseq/mini_ndm1__trace#0"],
    ["cli_bubbles_r60", "cli_bubbles_r100", "cli_bubbles_steps150"], ["cli_quorum_tuned", "cli_quorum_exhaustive"]]
RADII = [("cli_bubbles_r60", "cli_bubbles_r100"), ("cli_bubbles_r60", "cli_bubbles_steps150")]
# library defects W5 finds, fixed on gr/labeled-traversal by a separate lane (FORMAT v2 G): each must still reproduce
# here until the fixed library is vendored, then its entry is removed (the test says so)
KNOWN_LIBRARY_DEFECTS = {
    "G2": ("graphlet_walks spell=tail tail_bp=0 gives the whole right flank (DESIGN §6: the last 0 outward bases)",
           lambda tool, args: tool == "graphlet_walks" and args.get("spell") == "tail" and args.get("tail_bp") == 0),
}


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
           ("graphlet_export", {"format": "fasta", "what": "walk"}), ("graphlet_subtrie", {"labels": ["no-such"]})]
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
                ("graphlet_sequence", {"arm": s, "segment": 0, "with_seed": True})]
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
        cls.docs = {}       # name -> (handle, full answer or None, seed index, W1's document text)
        cls.pairs = engine_pairs()
        for cell, i, g, f in cls.pairs:
            text = parser.standalone_text(g["results"][i]["graphlet"], g["results"][i], g)
            name = "%s#%d" % (cell, i)
            cls.docs[name] = (cls._load(name, text), f, i, text)
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
                # a row cut to fit comes alone, the cut named on its page (DESIGN §6): kept on the row, where an
                # assembled answer can still say so (FORMAT question Q5)
                rows_ = [dict(r) for r in page.get("rows") or []]
                if page.get("row_truncated") and len(rows_) == 1:
                    rows_[0].update(row_truncated=True, cut_fields=page.get("cut_fields") or [])
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
        defects = {}
        n_calls = n_checked = 0
        declared = collections.Counter()
        for name, (h, f, i, text) in self.docs.items():
            if f is None:
                continue
            for tool, args in _calls(f["results"][i]):
                kw = dict(args, handle=h)
                if tool == "graphlet_save":
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
                diffs = O.differences(rec, lib)
                known = [k for k, (_, hit) in KNOWN_LIBRARY_DEFECTS.items() if hit(tool, op)]
                if known:
                    defects.setdefault(known[0], 0)
                    defects[known[0]] += bool(diffs)
                    continue
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
        for k, (what, _) in KNOWN_LIBRARY_DEFECTS.items():
            self.assertGreater(defects.get(k, 0), 0, "known library defect %s no longer reproduces (%s): the fixed "
                                                    "library is here, remove it from KNOWN_LIBRARY_DEFECTS" % (k, what))
        self.assertGreater(n_calls, 3500)
        self.assertGreater(n_checked, 150000)
        self.assertEqual(declared["twin_mismatch"] > 0, True)

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
