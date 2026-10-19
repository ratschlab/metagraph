"""Offline tests of graphlet_paging (FORMAT v2 A1-A4): the paging loop's helpers, the assembly of pages into one
answer, and the list table against the library's source (read as text: no library import, like the module)."""
import os
import re
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.dirname(HERE)
LIB = os.path.normpath(os.path.join(SUITE, "..", "..", "api", "python"))
sys.path.insert(0, SUITE)

import graphlet_paging as gp  # noqa: E402


class LoopTests(unittest.TestCase):
    def test_paging_modes(self):
        for good, most in (("none", 1), ("all", gp.MAX_PAGES), ({"pages": 3}, 3)):
            self.assertIsNone(gp.paging_problem(good))
            self.assertEqual(gp.max_pages(good), most)
        for bad in ("some", {"pages": 0}, {"pages": True}, {"pages": 2, "x": 1}, None, 3):
            self.assertIn("paging is", gp.paging_problem(bad))

    def test_next_call(self):
        self.assertEqual(gp.next_call("graphlet_walks", {"rows": [], "next_cursor": "c1"}), ("cursor", "c1"))
        self.assertIsNone(gp.next_call("graphlet_walks", {"rows": []}))
        self.assertEqual(gp.next_call("graphlet_sequence", {"sequence": "A", "next_from": 0}), ("start", 0))
        self.assertIsNone(gp.next_call("graphlet_sequence", {"sequence": "A", "next_cursor": "x"}))
        self.assertIsNone(gp.next_call("graphlet_walks", None))

    def test_paging_complete(self):
        self.assertTrue(gp.paging_complete("graphlet_walks", {"rows": [1]}, {"rows": [1]}))
        self.assertFalse(gp.paging_complete("graphlet_walks", {"rows": [1], "next_cursor": "c"}, {"rows": [1]}))
        self.assertFalse(gp.paging_complete("graphlet_claims", {"rows": [1], "complete": False},
                                            {"rows": [1], "complete": False}))
        self.assertTrue(gp.paging_complete("graphlet_claims", {"rows": [1], "complete": True},
                                           {"rows": [1], "complete": True}))
        self.assertFalse(gp.paging_complete("graphlet_sequence", {"sequence": "A", "next_from": 1}, {"sequence": "A"}))
        self.assertTrue(gp.paging_complete("graphlet_summary", {"summary": {}}, {"summary": {}}))

    def test_runaway(self):
        """VL6: a continuation that repeats stops the loop at once; "all" stops at MAX_PAGES; a pages cap is none."""
        a, b = {"rows": [1], "next_cursor": "c"}, {"rows": [2], "next_cursor": "c"}
        self.assertTrue(gp.more("graphlet_walks", "all", [a]))
        self.assertIsNone(gp.loop_error("graphlet_walks", "all", [a]))
        self.assertFalse(gp.more("graphlet_walks", "all", [a, b]))
        self.assertEqual(gp.loop_error("graphlet_walks", "all", [a, b]),
                         {"code": "paging_runaway", "fields": {"pages": 2, "reason": "repeats"}})
        # the first page against the op's own cursor; paging "none" or a reached pages cap asks for no more
        self.assertEqual(gp.loop_error("graphlet_walks", "all", [a], ("cursor", "c"))["fields"]["reason"], "repeats")
        self.assertIsNone(gp.loop_error("graphlet_walks", "none", [a], ("cursor", "c")))
        self.assertIsNone(gp.loop_error("graphlet_walks", {"pages": 2}, [a, b]))
        self.assertFalse(gp.more("graphlet_walks", {"pages": 2}, [a, {"rows": [2], "next_cursor": "d"}]))
        # a sequence's next_from must go beyond the page's own start
        seq = [{"from": 5, "to": 9, "sequence": "ACGT", "next_from": 5}]
        self.assertEqual(gp.loop_error("graphlet_sequence", "all", seq)["fields"]["reason"], "repeats")
        seq[0]["next_from"] = 9
        self.assertIsNone(gp.loop_error("graphlet_sequence", "all", seq))
        self.assertTrue(gp.more("graphlet_sequence", "all", seq))
        many = [{"rows": [i], "next_cursor": "c%d" % i} for i in range(gp.MAX_PAGES)]
        self.assertEqual(gp.loop_error("graphlet_walks", "all", many),
                         {"code": "paging_runaway", "fields": {"pages": gp.MAX_PAGES, "reason": "max_pages"}})
        self.assertIsNone(gp.loop_error("graphlet_walks", "all", many[:-1] + [{"rows": [0]}]))

    def test_selector_view(self):
        """VL5: a listing stopped with a cursor left keeps it for graphlet_select's cursor_of; nothing else does."""
        answer = {"rows": [1], "total": 2}
        self.assertEqual(gp.selector_view("graphlet_walks", answer, {"rows": [1], "next_cursor": "c"}),
                         {"rows": [1], "total": 2, "_next_cursor": "c"})
        self.assertEqual(answer, {"rows": [1], "total": 2})
        self.assertIs(gp.selector_view("graphlet_walks", answer, {"rows": [1]}), answer)
        self.assertEqual(gp.selector_view("graphlet_sequence", {"sequence": "A"}, {"sequence": "A", "next_from": 1}),
                         {"sequence": "A"})
        sys.path.insert(0, SUITE)
        import graphlet_select
        self.assertEqual(gp.SELECTOR_CURSOR_KEY, graphlet_select.CURSOR_KEY)

    def test_error_pages(self):
        page = {"error": "local_budget_exceeded", "message": "m", "hint": "h", "stop": {"resource": "work"},
                "max_bytes": 64}
        self.assertTrue(gp.error_page(page))
        self.assertFalse(gp.error_page({"rows": []}))
        self.assertFalse(gp.error_page({"error": None}))
        self.assertEqual(gp.error_record(page), {"code": "local_budget_exceeded",
                                                 "fields": {"stop": {"resource": "work"}, "max_bytes": 64}})


class AssembleTests(unittest.TestCase):
    def test_rows_concatenated_and_end_fields_last(self):
        pages = [{"handle": "h", "evidence": {"x": 1}, "total": 3, "rows": [1, 2], "next_cursor": "c1"},
                 {"handle": "h", "evidence": {"x": 1}, "total": 3, "rows": [3]}]
        self.assertEqual(gp.assemble("graphlet_walks", pages, True),
                         {"handle": "h", "evidence": {"x": 1}, "total": 3, "rows": [1, 2, 3]})
        budgeted = [{"handle": "h", "total_at_least": 2, "complete": False, "filtered_so_far": {"m": 1},
                     "rows": [1, 2], "next_cursor": "c"},
                    {"handle": "h", "total": 3, "complete": True, "filtered": {"m": 2}, "rows": [3],
                     "local": {"complete": True}}]
        self.assertEqual(gp.assemble("graphlet_claims", budgeted, True),
                         {"handle": "h", "total": 3, "complete": True, "filtered": {"m": 2}, "rows": [1, 2, 3],
                          "local": {"complete": True}})

    def test_evidence_of_the_completing_page(self):
        marked = {"exact": True, "limitations": {"left": ["scope"], "local": ["local_work"]}}
        plain = {"exact": True, "limitations": {"left": ["scope"]}}
        pages = [{"handle": "h", "evidence": marked, "total_at_least": 1, "complete": False, "rows": [1],
                  "next_cursor": "c"}, {"handle": "h", "evidence": plain, "total": 2, "complete": True, "rows": [2]}]
        self.assertEqual(gp.assemble("graphlet_walk", pages, True),
                         {"handle": "h", "evidence": plain, "total": 2, "complete": True, "rows": [1, 2]})
        self.assertEqual(gp.assemble("graphlet_walk", pages[:1], False)["evidence"], marked)
        pages[1]["evidence"] = dict(plain, exact=False)
        with self.assertRaises(gp.AssemblyError):
            gp.assemble("graphlet_walk", pages, True)

    def test_keys_of_later_pages_kept(self):
        """VL1: a key page 0 lacks is taken from the first page stating it; hint and filtered from the last page (a
        budget-resumed claims listing states both only once it has filtered something)."""
        pages = [{"handle": "h", "arm": "left", "total_at_least": 1, "complete": False, "rows": [1], "next_cursor": "a"},
                 {"handle": "h", "arm": "left", "total_at_least": 2, "complete": False, "rows": [2],
                  "filtered_so_far": {"merge_entered": 1}, "hint": "x", "caveat": "c", "next_cursor": "b"},
                 {"handle": "h", "arm": "left", "total": 3, "complete": True, "rows": [3],
                  "filtered": {"merge_entered": 2}, "hint": "x", "caveat": "c"}]
        self.assertEqual(gp.assemble("graphlet_claims", pages, True),
                         {"handle": "h", "arm": "left", "total": 3, "complete": True, "rows": [1, 2, 3],
                          "filtered": {"merge_entered": 2}, "hint": "x", "caveat": "c"})
        pages[2]["caveat"] = "other"
        with self.assertRaises(gp.AssemblyError):
            gp.assemble("graphlet_claims", pages, True)
        self.assertIn("hint", gp.END_FIELDS)

    def test_cursor_removed_and_volatile_keys(self):
        """VL4: cursor goes with next_cursor; |volatile| keys may differ page by page and come from the last page
        stating them; any other key that differs is still an AssemblyError."""
        pages = [{"handle": "h", "rows": [1], "total": 2, "cursor": None, "elapsed_ms": 3, "lane": "fast",
                  "next_cursor": "c"},
                 {"handle": "h", "rows": [2], "total": 2, "cursor": "c", "elapsed_ms": 4}]
        with self.assertRaises(gp.AssemblyError):
            gp.assemble("graphlet_walks", pages, True)
        self.assertEqual(gp.assemble("graphlet_walks", pages, True, volatile=("elapsed_ms", "lane")),
                         {"handle": "h", "rows": [1, 2], "total": 2, "elapsed_ms": 4, "lane": "fast"})
        pages[1]["handle"] = "k"
        with self.assertRaises(gp.AssemblyError):
            gp.assemble("graphlet_walks", pages, True, volatile=("elapsed_ms", "lane"))
        self.assertNotIn("cursor", gp.assemble("graphlet_summary", [{"summary": {}, "cursor": "x"}], True))

    def test_compare_cuts_merged(self):
        pages = [{"a": "x", "mode": "claims", "rows": [1], "total": 3, "fields_cut": ["evidence"], "export": "e",
                  "next_cursor": "c"},
                 {"a": "x", "mode": "claims", "rows": [2], "total": 3, "evidence": {"a": 1}, "export": "e",
                  "next_cursor": "c"},
                 {"a": "x", "mode": "claims", "rows": [3], "total": 3, "fields_cut": ["evidence", "export"]}]
        self.assertEqual(gp.assemble("graphlet_compare", pages, True),
                         {"a": "x", "mode": "claims", "rows": [1, 2, 3], "total": 3, "evidence": {"a": 1},
                          "export": "e"})
        pages[1].pop("evidence")
        self.assertEqual(gp.assemble("graphlet_compare", pages, True)["fields_cut"], ["evidence"])

    def test_truncated_rows_named(self):
        pages = [{"rows": [{"s": "AC"}], "total": 2, "row_truncated": True, "cut_fields": ["s"], "next_cursor": "c"},
                 {"rows": [{"s": "G"}], "total": 2}]
        out = gp.assemble("graphlet_walks", pages, True)
        # A3: the top level merged, each row of the cut page marked with that page's own cut_fields
        self.assertEqual((out["row_truncated"], out["cut_fields"], out["rows"]),
                         (True, ["s"], [{"s": "AC", "_cut_fields": ["s"]}, {"s": "G"}]))
        self.assertEqual(pages[0]["rows"], [{"s": "AC"}])   # the pages are not changed
        pages.append({"rows": [{"s": "T", "l": [1]}], "total": 3, "row_truncated": True, "cut_fields": ["l", "s"]})
        pages[1]["next_cursor"] = "c2"
        out = gp.assemble("graphlet_walks", pages, True)
        self.assertEqual(out["cut_fields"], ["s", "l"])
        self.assertEqual([r.get("_cut_fields") for r in out["rows"]], [["s"], None, ["l", "s"]])

    def test_conflicting_pages_refused(self):
        pages = [{"handle": "h", "rows": [1], "total": 2, "next_cursor": "c"}, {"handle": "k", "rows": [2], "total": 2}]
        with self.assertRaises(gp.AssemblyError):
            gp.assemble("graphlet_walks", pages, True)
        with self.assertRaises(gp.AssemblyError):
            gp.assemble("graphlet_walks", [{"rows": [1], "next_cursor": "c"}, {"total": 1}], True)
        with self.assertRaises(gp.AssemblyError):
            gp.assemble("graphlet_summary", [{"summary": {}}, {"summary": {}}], True)

    def test_sequence(self):
        pages = [{"from": 0, "to": 3, "length": 7, "sequence": "ACG", "truncated": True, "next_from": 3},
                 {"from": 3, "to": 6, "length": 7, "sequence": "TAC", "truncated": True, "next_from": 6},
                 {"from": 6, "to": 7, "length": 7, "sequence": "G"}]
        self.assertEqual(gp.assemble("graphlet_sequence", pages, True),
                         {"from": 0, "to": 7, "length": 7, "sequence": "ACGTACG"})
        self.assertEqual(gp.assemble("graphlet_sequence", pages[:2], False),
                         {"from": 0, "to": 6, "length": 7, "sequence": "ACGTAC", "truncated": True})

    def test_single_objects_pass(self):
        self.assertEqual(gp.assemble("graphlet_summary", [{"summary": {"a": 1}}], True), {"summary": {"a": 1}})
        self.assertEqual(gp.assemble("graphlet_export", [{"path": "p", "bytes": 3}], True), {"path": "p", "bytes": 3})


class ListTableTests(unittest.TestCase):
    def test_table_matches_the_library(self):
        """Every self._page / self._page_cut call of mcp_tools.py names a tool of LIST_FIELDS and every tool there has
        one; _page builds the page's list as `rows`; graphlet_sequence states next_from."""
        with open(os.path.join(LIB, "metagraph", "traverse", "mcp_tools.py"), encoding="utf-8") as fh:
            src = fh.read()
        paged = set(re.findall(r"self\._page(?:_cut)?\(\s*'(\w+)'", src))
        self.assertEqual(paged, set(gp.LIST_FIELDS))
        self.assertEqual(set(gp.LIST_FIELDS.values()), {"rows"})
        page = src[src.index("    def _page(self"):src.index("    def _page_cut(self")]
        self.assertIn("rows=[])", page)
        self.assertIn("out['next_cursor']", page)
        self.assertIn("row_truncated=True, cut_fields=", page)
        sequence = src[src.index("    def graphlet_sequence(self"):src.index("    def graphlet_export(self")]
        self.assertIn("next_from=start + len(part)", sequence)


if __name__ == "__main__":
    unittest.main()
