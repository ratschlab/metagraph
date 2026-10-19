"""Offline tests of graphlet_select (FORMAT-v8 §2): JSON pointers, op references, the named selectors on answers in
the library's shapes (mcp_tools.py), and every failure as a SelectorError. Stdlib only, no library import (the
module's own rule); test_graphlet_local runs the selectors on the library's real answers through the runner."""
import copy
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.dirname(HERE)
sys.path.insert(0, SUITE)

import graphlet_select as gs  # noqa: E402
from graphlet_select import SelectorError  # noqa: E402


def lab(i, name=None):
    return {"name": name or "label %d" % i, "ref": "c:%d" % i}


def row(segment, ins, outs, merge_of=None, cut_out=None):
    """A graphlet_walk row (mcp_tools._walk_rows) with |ins| / |outs| label ids; |cut_out|: the true out count."""
    r = {"segment": segment, "from_bp": 10 * segment, "to_bp": 10 * segment + 10, "labels_in": len(ins),
         "labels_in_total": len(ins), "labels_out": cut_out or len(outs), "labels_out_total": cut_out or len(outs),
         "exact": True, "in": [lab(i) for i in ins], "out": [lab(i) for i in outs]}
    if merge_of:
        r["merge_of"] = merge_of
    return r


def walk(arm, walk_id, rows, continuation=None, length_bp=None):
    a = {"handle": "<handle>", "arm": arm, "walk": walk_id, "length_bp": length_bp or 10 * len(rows),
         "evidence": {"complete_to_bp": {arm: 30}}, "total": len(rows), "rows": rows}
    if continuation is not None:
        a["continuation"] = continuation
    return a


def answers():
    """A two-arm graphlet's answers: walks 0 (left, no merge, a continuation), 2 and 5 (left, through merge segment 3,
    entered by c:1 and c:2 from parent 2), 1 (right, a continuation with no labels), a label table, a walks list, a
    sequence, a summary, a one-arm summary and an error."""
    return {
        "walk-l0": walk("left", 0, [row(0, [0, 1], [0, 1]), row(1, [0, 1], [0])],
                        continuation={"length_bp": 15, "loss_used": 0.0, "labels": 1}),
        "walk-l5": walk("left", 5, [row(0, [0, 1], [0, 1]), row(1, [0, 1], [0]), row(3, [0, 1, 2], [0, 1, 2],
                                                                                         merge_of=[1, 2])]),
        "walk-l2": walk("left", 2, [row(0, [0, 1], [0, 1]), row(1, [0, 1], [0]), row(3, [0, 1, 2], [0, 2],
                                                                                         merge_of=[1, 2])],
                        continuation={"length_bp": 15, "loss_used": 1.0, "labels": 2}),
        "walk-r1": walk("right", 1, [row(0, [0], [0])], continuation={"length_bp": 15, "loss_used": 0.0,
                                                                         "labels": 0}),
        "walk-cut": walk("left", 7, [row(0, [0, 1], [0, 1]), row(1, [0, 1], [0], cut_out=9),
                                     row(3, [0, 1, 2], [0], merge_of=[1, 2])]),
        "walk-same": walk("left", 8, [row(0, [0, 1], [0, 1]), row(1, [0, 1], [0, 1]),
                                      row(3, [0, 1], [0, 1], merge_of=[1, 2])]),
        "walk-wrong": walk("left", 9, [row(0, [0], [0]), row(3, [0, 1], [0, 1], merge_of=[2, 1])]),
        "labels": {"handle": "<handle>", "evidence": {}, "total": 3, "rows": [
            dict(lab(4), arm="right", direct_bp=50, reach_bp=60, reentries=0, runs=1),
            dict(lab(1), arm="left", direct_bp=40, reach_bp=40, reentries=0, runs=1),
            dict(lab(2), arm="left", direct_bp=30, reach_bp=45, reentries=1, runs=2)]},
        "label-walks": {"name": "label 1", "ref": "c:1", "handle": "<handle>", "evidence": {}, "total": 1,
                        "rows": [{"arm": "left", "from_bp": 0, "to_bp": 20, "evidence_from": 0, "end": "D",
                                  "walks_below": 1, "merged_into": None, "route": [0, 1]}]},
        "walks": {"handle": "<handle>", "arm": "left", "evidence": {}, "total": 2, "rows": [
            {"walk": 5, "length_bp": 30, "end": "dead_end"}, {"walk": 0, "length_bp": 20, "end": "dead_end"}]},
        "seq": {"handle": "<handle>", "arm": "left", "orientation": "natural", "with_seed": True, "from": 0,
                "length": 70, "ceiling_bytes": 16384, "evidence": {}, "to": 70, "sequence": "A" * 70},
        "summary": {"handle": "<handle>", "evidence": {"complete_to_bp": {"left": 30, "right": 10}},
                    "summary": {"arms": {"right": {}, "left": {}}}},
        "summary-one": {"handle": "<handle>", "evidence": {"complete_to_bp": {"right": 10}},
                        "summary": {"arms": {"right": {}}}},
        "evidence-only": {"handle": "<handle>", "evidence": {"complete_to_bp": {"left": 3, "right": 4}}},
        "failed": None,
        "error": {"error": "bad_arm", "message": "x"},
    }


def sel(name, **args):
    return {"select": name, "args": args}


class PointerTests(unittest.TestCase):
    DOC = {"a": [{"b": 1}, {"c/d": 2, "e~f": 3, "": 4}], "n": None, "s": "text"}

    def test_pointers(self):
        p = gs.pointer
        self.assertIs(p(self.DOC, ""), self.DOC)
        self.assertEqual(p(self.DOC, "/a/0/b"), 1)
        self.assertEqual(p(self.DOC, "/a/1/c~1d"), 2)
        self.assertEqual(p(self.DOC, "/a/1/e~0f"), 3)
        self.assertEqual(p(self.DOC, "/a/1/"), 4)
        self.assertIsNone(p(self.DOC, "/n"))
        self.assertEqual(p([10, 11], "/1"), 11)

    def test_pointer_failures(self):
        for path in ("a", "/x", "/a/2", "/a/01", "/a/-", "/a/b", "/s/0", "/n/x", "/a/1/e~2f", "/a/1/e~"):
            with self.subTest(path):
                with self.assertRaises(SelectorError) as cm:
                    gs.pointer(self.DOC, path)
                self.assertIn(path, str(cm.exception))


class ResolveTests(unittest.TestCase):
    def test_literals_pass(self):
        for v in (1, "left", None, True, 2.5, [1, 2], {"ref": "c:3"}, {"name": "x"}, {"id": 2},
                  {"work_units": 5}, {"ref": "c:3", "other": 1}, {"select": "x", "other": 1}):
            with self.subTest(v):
                self.assertEqual(gs.resolve(v, answers()), v)

    def test_refs_at_any_depth(self):
        ans = answers()
        args = {"arm": "left", "walk": {"ref": "walks", "path": "/rows/0/walk"},
                "walks": [{"ref": "walks", "path": "/rows/1/walk"}, 3],
                "labels": [{"ref": "c:1"}, {"ref": "labels", "path": "/rows/2/ref"}],
                "label": {"ref": "labels", "path": "/rows/0"}}
        out = gs.resolve(args, ans)
        self.assertEqual(out, {"arm": "left", "walk": 5, "walks": [0, 3], "labels": [{"ref": "c:1"}, "c:2"],
                               "label": dict(lab(4), arm="right", direct_bp=50, reach_bp=60, reentries=0, runs=1)})
        out["label"]["ref"] = "changed"   # a copy: the answers are never changed through a resolved value
        self.assertEqual(ans["labels"]["rows"][0]["ref"], "c:4")
        self.assertEqual(args["walk"], {"ref": "walks", "path": "/rows/0/walk"})

    def test_ref_failures(self):
        cases = {
            "unknown op": ({"ref": "nope", "path": "/rows"}, "no op 'nope'"),
            "failed op": ({"ref": "failed", "path": "/rows"}, "has no answer"),
            "error answer": ({"ref": "error", "path": "/rows"}, "answered the error bad_arm"),
            "missing path": ({"ref": "walks", "path": "/rows/9/walk"}, "the list has 2 items"),
            "not a string": ({"ref": 3, "path": "/rows"}, "a reference is"),
        }
        for case, (value, text) in cases.items():
            with self.subTest(case):
                with self.assertRaises(SelectorError) as cm:
                    gs.resolve({"walk": value}, answers())
                self.assertIn(text, str(cm.exception))
                self.assertTrue(str(cm.exception).startswith("args.walk: "), str(cm.exception))

    def test_select_failures(self):
        cases = {
            "unknown": ({"select": "nope"}, "no named selector 'nope'"),
            "name a list": ({"select": ["top_label"]}, "named by a string, not ['top_label']"),
            "name an object": ({"select": {"a": 1}, "args": {}}, "named by a string"),
            "name a number": ({"select": 3}, "named by a string, not 3"),
            "args not an object": ({"select": "top_label", "args": [1]}, "are an object"),
            "unexpected argument": (sel("top_label", labels="labels", colour="red"), "colour"),
            "missing argument": (sel("top_label"), "labels"),
            "nested failure": (sel("top_label", labels={"ref": "nope", "path": ""}), "args.x.top_label.args.labels"),
        }
        for case, (value, text) in cases.items():
            with self.subTest(case):
                with self.assertRaises(SelectorError) as cm:
                    gs.resolve({"x": value}, answers())
                self.assertIn(text, str(cm.exception))

    def test_cursor_of(self):
        """VL5: the last next_cursor of a listing that stopped early (graphlet_paging.selector_view keeps it)."""
        ans = answers()
        ans["walks-p1"] = dict(ans["walks"], **{gs.CURSOR_KEY: "cur-1"})
        self.assertEqual(gs.resolve({"cursor": {"cursor_of": "walks-p1"}}, ans), {"cursor": "cur-1"})
        self.assertEqual(gs.resolve([{"cursor_of": "walks-p1"}], ans), ["cur-1"])
        for op_id, text in (("walks", "followed to its end"), ("seq", "followed to its end"),
                            ("failed", "has no answer"), ("nope", "no op 'nope'"), (["x"], "not ['x']")):
            with self.subTest(op_id):
                with self.assertRaises(SelectorError) as cm:
                    gs.resolve({"cursor": {"cursor_of": op_id}}, ans)
                self.assertIn(text, str(cm.exception))
                self.assertTrue(str(cm.exception).startswith("args.cursor: "))
        self.assertEqual(gs.resolve({"cursor_of": "x", "other": 1}, ans), {"cursor_of": "x", "other": 1})

    def test_select_with_a_resolved_argument(self):
        value = sel("slice_past_end", length_of="walks", walk={"ref": "walks", "path": "/rows/1/walk"}, part="end",
                    past=5)
        self.assertEqual(gs.resolve(value, answers()), 25)

    def test_every_selector_is_documented(self):
        self.assertTrue({"walk_through_merge", "walk_with_continuation", "slice_past_end", "label_entered_at_merge",
                         "top_label", "second_arm_or_none"} <= set(gs.SELECTORS))
        for name, fn in gs.SELECTORS.items():
            with self.subTest(name):
                self.assertGreater(len((fn.__doc__ or "").split()), 25)
                self.assertEqual(fn.__name__, name)


class NamedSelectorTests(unittest.TestCase):
    def r(self, name, **args):
        return gs.resolve(sel(name, **args), answers())

    def fails(self, text, name, **args):
        with self.assertRaises(SelectorError) as cm:
            self.r(name, **args)
        self.assertIn(text, str(cm.exception))

    def test_walk_through_merge(self):
        self.assertEqual(self.r("walk_through_merge", walks=["walk-l0", "walk-l5", "walk-l2"]), 2)
        self.assertEqual(self.r("walk_through_merge", walks=["walk-l0", "walk-l5", "walk-l2"], order="given"), 5)
        self.assertEqual(self.r("walk_through_merge", walks="walk-l5"), 5)
        self.assertEqual(self.r("walk_through_merge", walks=["walk-l5", "walk-r1"], arm="left"), 5)
        self.fails("none of the walks of walk-l0, walk-r1", "walk_through_merge", walks=["walk-l0", "walk-r1"])
        self.fails("on the right arm", "walk_through_merge", walks=["walk-l5"], arm="right")
        self.fails("not a graphlet_walk answer", "walk_through_merge", walks=["labels"])
        self.fails("non-empty list of op ids", "walk_through_merge", walks=[])
        self.fails("order is one of", "walk_through_merge", walks="walk-l5", order="length")
        self.fails("has no answer", "walk_through_merge", walks=["walk-l5", "failed"])

    def test_walk_with_continuation(self):
        ops = ["walk-l5", "walk-l2", "walk-l0"]
        self.assertEqual(self.r("walk_with_continuation", walks=ops), 0)
        self.assertEqual(self.r("walk_with_continuation", walks=ops, order="given"), 2)
        self.assertEqual(self.r("walk_with_continuation", walks=ops, min_labels=2), 2)
        self.assertEqual(self.r("walk_with_continuation", walks=["walk-r1"]), 1)
        self.fails("with 1 labels or more", "walk_with_continuation", walks=["walk-r1"], min_labels=1)
        self.fails("lie on both arms", "walk_with_continuation", walks=["walk-l0", "walk-r1"])
        self.assertEqual(self.r("walk_with_continuation", walks=["walk-l0", "walk-r1"], arm="right"), 1)
        self.fails("min_labels is an integer", "walk_with_continuation", walks=ops, min_labels=-1)

    def test_slice_past_end(self):
        self.assertEqual(self.r("slice_past_end", length_of="seq", part="end"), 71)
        self.assertEqual(self.r("slice_past_end", length_of="seq", part="start", keep=20), 50)
        self.assertEqual(self.r("slice_past_end", length_of="seq", part="start"), 70)
        self.assertEqual(self.r("slice_past_end", length_of="seq", part="start", keep=500), 0)
        self.assertEqual(self.r("slice_past_end", length_of="walk-l5", part="end", past=10), 40)
        self.assertEqual(self.r("slice_past_end", length_of="walk-l5", walk=5, part="start", keep=5), 25)
        self.assertEqual(self.r("slice_past_end", length_of="walks", walk=0, part="end"), 21)
        self.fails("name the walk", "slice_past_end", length_of="walks", part="end")
        self.fails("lists no walk 4", "slice_past_end", length_of="walks", walk=4, part="end")
        self.fails("answered walk 5, not walk 2", "slice_past_end", length_of="walk-l5", walk=2, part="end")
        self.fails("past is an integer >= 1", "slice_past_end", length_of="seq", part="end", past=0)
        self.fails("part is one of", "slice_past_end", length_of="seq", part="middle")
        self.fails("not a graphlet_sequence", "slice_past_end", length_of="summary", part="end")

    def test_label_entered_at_merge(self):
        self.assertEqual(self.r("label_entered_at_merge", walk="walk-l5"), {"ref": "c:1"})
        self.assertEqual(self.r("label_entered_at_merge", walk="walk-l5", form="ref"), "c:1")
        self.assertEqual(self.r("label_entered_at_merge", walk="walk-l5", form="name"), "label 1")
        # the first walk whose merge shows a label from another parent: walk-same's merge brings none
        self.assertEqual(self.r("label_entered_at_merge", walk=["walk-same", "walk-l2"]), {"ref": "c:1"})
        self.fails("also leaves the first parent", "label_entered_at_merge", walk="walk-same")
        self.fails("no merge on the walks of walk-l0", "label_entered_at_merge", walk="walk-l0")
        self.fails("shows 1 of its 9 out labels", "label_entered_at_merge", walk="walk-cut")
        self.fails("not its first parent 2", "label_entered_at_merge", walk="walk-wrong")
        self.fails("form is one of", "label_entered_at_merge", walk="walk-l5", form="id")

    def test_top_label(self):
        self.assertEqual(self.r("top_label", labels="labels"), {"ref": "c:4"})
        self.assertEqual(self.r("top_label", labels="labels", arm="left"), {"ref": "c:1"})
        self.assertEqual(self.r("top_label", labels="labels", arm="left", index=1, form="name"), "label 2")
        self.assertEqual(self.r("top_label", labels="labels", form="ref"), "c:4")
        self.fails("lists 2 label rows on the left arm, no row 2", "top_label", labels="labels", arm="left", index=2)
        self.fails("not a graphlet_labels table answer", "top_label", labels="label-walks")
        self.fails("not a graphlet_labels table answer", "top_label", labels="walks")
        self.fails("arm is one of", "top_label", labels="labels", arm="up")

    def test_arms(self):
        self.assertEqual(self.r("second_arm_or_none", of="summary"), "right")
        self.assertIsNone(self.r("second_arm_or_none", of="summary-one"))
        self.assertEqual(self.r("second_arm_or_none", of="evidence-only"), "right")
        self.assertEqual(self.r("first_arm", of="summary"), "left")
        self.assertEqual(self.r("first_arm", of="summary-one"), "right")
        self.fails("names no arms", "second_arm_or_none", of="labels")
        self.fails("has no answer", "second_arm_or_none", of="failed")

    def test_answers_untouched(self):
        ans = answers()
        before = copy.deepcopy(ans)
        for name, args in (("walk_through_merge", {"walks": ["walk-l5", "walk-l2"]}),
                           ("label_entered_at_merge", {"walk": "walk-l5"}), ("top_label", {"labels": "labels"})):
            gs.resolve(sel(name, **args), ans)
        self.assertEqual(ans, before)


if __name__ == "__main__":
    unittest.main()
