"""Offline tests of the results document generator (doc/make_doc.py; DESIGN-DOC §11): the prose data complete for
every unit and family of the v4 mapping, the mapping complete over the v4 manifest, every check kind described,
the generator on the fake run of test_compare_report (the .tex, the index, the verdict, the escaping; the PDF when
xelatex is installed), the cell grammar and the expectation formatters on chapter R's real expectations."""

import json
import os
import re
import shutil
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.dirname(HERE)
DOC = os.path.join(SUITE, "doc")
for p in (SUITE, HERE, DOC):
    if p not in sys.path:
        sys.path.insert(0, p)

import compare as cmp  # noqa: E402
import merge_run  # noqa: E402
import report  # noqa: E402
import suite_common as sc  # noqa: E402
from test_compare_report import CASES, build_run  # noqa: E402

import blocks  # noqa: E402
import cells  # noqa: E402
import latex  # noqa: E402
import make_doc  # noqa: E402
import reader  # noqa: E402
import units as units_mod  # noqa: E402

V4_MANIFEST = os.path.join(SUITE, "manifest", "suite_v4.json")
V4_MAPPING = os.path.join(DOC, "units", "suite_v4.json")
PROSE_DIR = os.path.join(DOC, "prose")


def v4():
    manifest = sc.read_json(V4_MANIFEST)
    mapping = units_mod.load_mapping(V4_MAPPING)
    return manifest, mapping, units_mod.resolve(manifest, mapping)


# ---------------------------------------------------------------------------------------------------------------
# the maintenance tests on the real v4 data

class MappingTest(unittest.TestCase):
    def test_mapping_complete(self):
        """Every v4 scenario in exactly one unit, units non-empty, types in the vocabulary, series of at most 8
        distinct variants, groups consistent, every relation under exactly one unit holding a or b."""
        manifest, mapping, units = v4()
        self.assertEqual(mapping["_manifest"]["content_sha256"], manifest["content_sha256"],
                         "the mapping was built for another v4 manifest")
        self.assertEqual(units_mod.validate(units, manifest), [])
        self.assertEqual(len(units.by_scenario), len(manifest["scenarios"]))
        self.assertEqual(units.chapters, ["R", "T", "P", "L", "X", "M", "Q", "C", "N", "D"])
        self.assertEqual(len(units.chapter_scenarios("R")), 81)
        self.assertEqual(units.unit_of("R-pos-0027"), "R-pos-column-B")
        self.assertEqual(units.relation_unit["REL-prefix_equals-0001"], "R-budget")

    def test_validate_names_problems(self):
        manifest, mapping, units = v4()
        broken = json.loads(json.dumps(mapping))
        broken["_units"]["R-order"]["unit_type"] = "matrix"
        del broken["R-twin-0001"]
        with self.assertRaises(units_mod.UnitsError) as ctx:
            units_mod.resolve(manifest, broken)
        self.assertIn("R-twin-0001", str(ctx.exception))
        units2 = units_mod.resolve(manifest, mapping)
        units2.units["R-order"]["unit_type"] = "matrix"
        units2.units["R-pos-column-A"]["members"].append("R-pos-0003")
        problems = units_mod.validate(units2, manifest)
        self.assertTrue(any("R-order" in p and "matrix" in p for p in problems), problems)
        self.assertTrue(any("R-pos-0003 in 2 unit" in p for p in problems), problems)

    def test_ids_text(self):
        self.assertEqual(units_mod.ids_text(["R-pos-0001", "R-pos-0005"]), "ids: \\texttt{R-pos-}0001 · 0005")
        self.assertEqual(units_mod.ids_text(["R-pos-0001", "R-pos-0002", "R-pos-0003", "R-pos-0007"]),
                         "ids: \\texttt{R-pos-}0001–0003 · 0007")
        self.assertIn("(table 1)", units_mod.ids_text(["X-a-0001", "X-a-0002"], None, {"table 1": ["X-a-0001", "X-a-0002"]}))


class ProseTest(unittest.TestCase):
    def test_prose_complete(self):
        """The owner's rule: a suite test fails if a family or unit lacks prose (every chapter an intro, every family a
        title and prose, every unit a property line); the message names the missing keys."""
        manifest, mapping, units = v4()
        missing = []
        for ch in units.chapters:
            path = os.path.join(PROSE_DIR, "%s.json" % ch)
            if not os.path.exists(path):
                missing.append("%s: no prose file" % ch)
                continue
            data = sc.read_json(path) or {}
            if not data.get("intro"):
                missing.append("%s: intro" % ch)
            fams, us = data.get("families") or {}, data.get("units") or {}
            for fid in units.chapter_families(ch):
                entry = fams.get(fid) or {}
                if not entry.get("title"):
                    missing.append("%s: families.%s.title" % (ch, fid))
                if not entry.get("prose"):
                    missing.append("%s: families.%s.prose" % (ch, fid))
                for uid in units.families[fid]["units"]:
                    if not us.get(uid):
                        missing.append("%s: units.%s" % (ch, uid))
        self.assertEqual(missing, [], "prose missing: " + "; ".join(missing))
        for name in ("overview.json", "glossary.json"):
            self.assertTrue(os.path.exists(os.path.join(PROSE_DIR, name)), name)

    def test_check_catalogue(self):
        """Every check kind any v4 scenario expects has a docstring in checks.py or a prose entry."""
        manifest, mapping, units = v4()
        checks = cmp.load_checks_module(SUITE)
        self.assertIsNotNone(checks)
        prose_checks = set()
        for ch in units.chapters:
            data = sc.read_json(os.path.join(PROSE_DIR, "%s.json" % ch)) or {}
            prose_checks.update(k for k, v in (data.get("checks") or {}).items() if v)
        missing = set()
        for s in manifest["scenarios"]:
            for e in s.get("expect") or []:
                name = e.get("check")
                fn = checks.CHECKS.get(name)
                if fn is None:
                    missing.add("%s (not in CHECKS)" % name)
                elif not (fn.__doc__ or "").strip() and name not in prose_checks:
                    missing.add(name)
        self.assertEqual(missing, set(), "check kinds without a description: %s" % ", ".join(sorted(missing)))

    def test_scaffold_keeps_text(self):
        manifest, mapping, units = v4()
        existing = {"chapter": "R", "intro": ["written"], "families": {"R-pos": {"title": "t", "prose": "written"}},
                    "units": {"R-pos-column-A": "written"}}
        data = units_mod.prose_scaffold(units, existing, "R", ["status"])
        self.assertEqual(data["intro"], ["written"])
        self.assertEqual(data["families"]["R-pos"]["prose"], "written")
        self.assertEqual(data["units"]["R-pos-column-A"], "written")
        self.assertIn("to be written", data["families"]["R-multi"]["prose"])
        self.assertEqual(data["units"]["R-multi"], units.units["R-multi"]["property"])
        self.assertIn("to be written", data["checks"]["status"])


# ---------------------------------------------------------------------------------------------------------------
# the generator on the fake run

FAKE_MAPPING = {
    "_chapters": {"P": "Patterns up to k", "R": "Resolve", "M": "Mismatches", "S": "Service"},
    "_families": {"P-pos": {"chapter": "P", "title": "positives", "units": ["P-pos-series", "P-pos-reps"]},
                  "P-edge": {"chapter": "P", "title": "edge inputs", "units": ["P-edge-refused"]},
                  "R-one": {"chapter": "R", "title": "one resolve", "units": ["R-one"]},
                  "M-one": {"chapter": "M", "title": "mismatches", "units": ["M-one"]},
                  "S-one": {"chapter": "S", "title": "service", "units": ["S-one"]}},
    "_units": {"P-pos-series": {"family": "P-pos", "unit_type": "series", "unit_title": "length", "property": "p1",
                                "members": ["P-001", "P-003"]},
               "P-pos-reps": {"family": "P-pos", "unit_type": "replicates", "unit_title": "repeats", "property": "p2",
                              "members": ["P-004", "P-005"], "groups": ["a", "b"]},
               "P-edge-refused": {"family": "P-edge", "unit_type": "refusals", "unit_title": "refused",
                                  "property": "p3", "members": ["P-002"]},
               "R-one": {"family": "R-one", "unit_type": "pairs", "unit_title": "pair", "property": "p4",
                         "members": ["R-001"]},
               "M-one": {"family": "M-one", "unit_type": "list", "unit_title": "list", "property": "p5",
                         "members": ["M-001", "M-002"]},
               "S-one": {"family": "S-one", "unit_type": "list", "unit_title": "service only", "property": "p6",
                         "members": ["S-001"]}},
    "P-001": {"chapter": "P", "family": "P-pos", "unit_id": "P-pos-series", "unit_type": "series", "variant": "10 bp"},
    "P-003": {"chapter": "P", "family": "P-pos", "unit_id": "P-pos-series", "unit_type": "series", "variant": "20 bp"},
    "P-004": {"chapter": "P", "family": "P-pos", "unit_id": "P-pos-reps", "unit_type": "replicates", "variant": "x",
              "group": "a"},
    "P-005": {"chapter": "P", "family": "P-pos", "unit_id": "P-pos-reps", "unit_type": "replicates", "variant": "y",
              "group": "b"},
    "P-002": {"chapter": "P", "family": "P-edge", "unit_id": "P-edge-refused", "unit_type": "refusals",
              "variant": "an unknown field"},
    "R-001": {"chapter": "R", "family": "R-one", "unit_id": "R-one", "unit_type": "pairs", "variant": "the one"},
    "M-001": {"chapter": "M", "family": "M-one", "unit_id": "M-one", "unit_type": "list", "variant": "slot"},
    "M-002": {"chapter": "M", "family": "M-one", "unit_id": "M-one", "unit_type": "list", "variant": "mismatches"},
    "S-001": {"chapter": "S", "family": "S-one", "unit_id": "S-one", "unit_type": "list", "variant": "only"},
}


class FakeRunTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.run_dir = os.path.join(cls.tmp, "run")
        os.makedirs(cls.run_dir)
        build_run(cls.run_dir)
        merge_run.merge(cls.run_dir)
        sc.write_json_atomic(os.path.join(cls.run_dir, "compare.json"), cmp.compare(cls.run_dir, (CASES, {"cases": 5})))
        cls.mapping_path = os.path.join(cls.tmp, "units.json")
        with open(cls.mapping_path, "w") as fh:
            json.dump(FAKE_MAPPING, fh)
        cls.prose_dir = os.path.join(cls.tmp, "prose")
        run = reader.Run(cls.run_dir)
        units = units_mod.resolve(run.manifest, FAKE_MAPPING)
        make_doc.scaffold_prose(cls.prose_dir, units, run)
        cls.out = os.path.join(cls.tmp, "out")
        cls.rc = make_doc.main([cls.run_dir, "--out", cls.out, "--no-pdf", "--units", cls.mapping_path,
                                "--prose-dir", cls.prose_dir])
        cls.runobj = run
        cls.units = units

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def tex(self, name="doc.tex"):
        with open(os.path.join(self.out, name), encoding="utf-8") as fh:
            return fh.read()

    def test_build_tex(self):
        self.assertEqual(self.rc, 0)
        main = self.tex()
        self.assertTrue(os.path.exists(os.path.join(self.out, "doc.sty")))
        for ch in ("P", "R", "M", "S"):
            self.assertIn("\\input{chapter_%s}" % ch, main)
            self.assertTrue(os.path.exists(os.path.join(self.out, "chapter_%s.tex" % ch)))
        # every scenario id exactly once in the index
        index = main.split("Index of test ids", 1)[1].split("\\end{multicols}", 1)[0]
        for sid in self.runobj.order:
            self.assertEqual(index.count(latex.tt(sid)), 1, sid)
        # the verdict word is report.py's
        text = report.build_report(self.run_dir)
        word = re.search(r"\*\*(\S+)", text).group(1)
        self.assertEqual(self.runobj.verdict.split()[0], word)
        self.assertIn("\\textbf{\\large %s}" % word, main)
        self.assertTrue(word.startswith("INCOMPLETE"))

    def test_escaping(self):
        chapter = self.tex("chapter_R.tex")
        # R-001 fails on direct-local and the arms disagree: a detail block with the escaped title
        self.assertIn("detail:R-001", chapter)
        self.assertIn("with {\\glyphfont ≤} and \\$x\\_1\\$ | pipe", chapter)
        self.assertNotIn("$x_1$", chapter)
        self.assertEqual(latex.esc("a_b & c% $d# {e} ~f ^g \\h"),
                         "a\\_b \\& c\\% \\$d\\# \\{e\\} \\textasciitilde{}f \\textasciicircum{}g \\textbackslash{}h")
        self.assertIn("(120 bp)", latex.esc("ACGT" * 30))
        self.assertTrue(latex.esc("x" * 50, 10).endswith("…"))
        self.assertEqual(latex.inline("a `b_c` *d*"), "a \\texttt{b\\_c} \\emph{d}")
        self.assertIn("\\allowbreak{}", latex.tt("https://x.y/z"))
        self.assertIn("\\allowbreak{}", latex.tt("a" * 40))

    def test_states(self):
        st = cells.scenario_states(self.runobj, "R-001")
        self.assertEqual(st["direct-local"].glyph_key, "fail")
        self.assertFalse(st["direct-local"].clean)
        self.assertEqual(st["direct-remote"].glyph_key, "pass")
        self.assertTrue(st["direct-remote"].ne)  # the arms disagree: the mark on direct-remote
        self.assertFalse(st["direct-remote"].clean)
        self.assertFalse(cells.is_clean(st))
        st = cells.scenario_states(self.runobj, "P-001")
        self.assertTrue(cells.is_clean(st))
        self.assertEqual(st["mcp"].glyph_tex(), "\\gpass")
        # the resumed duplicate of P-001 r0 (9 ms) wins over the first line (10 ms)
        self.assertEqual(cells.fmt_times(st["direct-local"].walls), "9 / 6 ms")
        st = cells.scenario_states(self.runobj, "P-005")
        self.assertEqual(st["mcp"].glyph_key, "not_implemented_mcp")
        self.assertTrue(st["mcp"].clean)
        st = cells.scenario_states(self.runobj, "P-004")
        self.assertEqual(st["mcp"].glyph_key, "service_by_design")
        st = cells.scenario_states(self.runobj, "M-002")
        self.assertFalse(st["mcp"].present)
        st = cells.scenario_states(self.runobj, "P-002")
        self.assertEqual(st["mcp"].glyph_key, "fail")
        self.assertIn("\\fn{", self.tex("chapter_P.tex"))

    def test_build_pdf(self):
        xelatex = make_doc.find_xelatex()
        if not xelatex:
            self.skipTest("xelatex not installed")
        out = os.path.join(self.tmp, "pdf")
        rc = make_doc.main([self.run_dir, "--out", out, "--units", self.mapping_path, "--prose-dir", self.prose_dir,
                            "--chapters", "P,R"])
        self.assertEqual(rc, 0)
        self.assertTrue(os.path.exists(os.path.join(out, "doc.pdf")))

    def test_missing_prose_and_mapping_fail(self):
        out = os.path.join(self.tmp, "fail")
        empty = os.path.join(self.tmp, "empty-prose")
        os.makedirs(empty, exist_ok=True)
        rc = make_doc.main([self.run_dir, "--out", out, "--no-pdf", "--units", self.mapping_path, "--prose-dir", empty])
        self.assertEqual(rc, 2)
        broken = json.loads(json.dumps(FAKE_MAPPING))
        del broken["S-001"]
        path = os.path.join(self.tmp, "broken.json")
        with open(path, "w") as fh:
            json.dump(broken, fh)
        rc = make_doc.main([self.run_dir, "--out", out, "--no-pdf", "--units", path, "--prose-dir", self.prose_dir])
        self.assertEqual(rc, 2)


# ---------------------------------------------------------------------------------------------------------------
# the cell grammar and the expectation formatters

class CellGrammarTest(unittest.TestCase):
    def test_fmt_t(self):
        self.assertEqual(cells.fmt_t(None), "–")
        self.assertEqual(cells.fmt_t(0.3), "<1 ms")
        self.assertEqual(cells.fmt_t(26.4), "26 ms")
        self.assertEqual(cells.fmt_t(99.6), "0.1 s")
        self.assertEqual(cells.fmt_t(3140), "3.1 s")
        self.assertEqual(cells.fmt_t(18400), "18 s")

    def test_fmt_times(self):
        self.assertEqual(cells.fmt_times([3100, 900]), "3.1 / 0.9 s")
        self.assertEqual(cells.fmt_times([26, 1200]), "26 ms / 1.2 s")
        self.assertEqual(cells.fmt_times([27, 30], ["c", None]), "27\\tmark{c} / 30 ms")
        self.assertEqual(cells.fmt_times([None, 5000]), "– / 5.0 s")
        self.assertEqual(cells.fmt_times([]), "–")

    def test_worst(self):
        self.assertEqual(cells.worst_of(["pass", "inconclusive"]), "inconclusive")
        self.assertEqual(cells.worst_of(["finding", "fail"]), "fail")
        self.assertEqual(cells.worst_of(["pass", "service_by_design"]), "service_by_design")


class ExpectTextTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        manifest = sc.read_json(V4_MANIFEST)
        cls.scen = dict((s["id"], s) for s in manifest["scenarios"])
        cls.is_property = staticmethod(cmp.load_checks_module(SUITE).is_property)

    def expected(self, sid):
        return blocks.expected_text(self.scen[sid], self.is_property)

    def test_formatters(self):
        self.assertEqual(self.expected("R-pos-0001"),
                         "http 200 · full under \\texttt{511145} · graph\\_runs \\texttt{[[0,1]]} · invariants · "
                         "budget respected")
        self.assertIn("http 400 \\texttt{", self.expected("R-edge-0003"))
        self.assertIn("longer than the server limit of 10000 bp", self.expected("R-edge-0003"))
        self.assertIn("prefix \\texttt{Bad json received", self.expected("R-edge-0019"))
        self.assertIn("prefix if stopped", self.expected("R-budget-0001"))
        self.assertIn("unless stopped", self.expected("R-twin-0003"))
        self.assertIn("39 k-mers under \\texttt{511145}", self.expected("R-edge-0005"))
        self.assertIn("graph\\_runs \\texttt{[[0,20],[51,70]]}", self.expected("R-edge-0005"))
        self.assertIn("discovery lists \\texttt{107806}", self.expected("R-discover-0001"))
        self.assertIn("0 k-mers under", self.expected("R-rc-0001"))
        self.assertIn("not full under \\texttt{243232}", self.expected("R-multi-0001"))
        clamped = self.expected("R-edge-0021")
        self.assertIn("clamped", clamped)
        self.assertIn("70000", clamped)
        self.assertIn("60000", clamped)
        self.assertIn("budget respected (cap 60 s)", clamped)
        self.assertEqual(blocks.expect_text("walk_in_record", {"a": 1}), "walk\\_in\\_record(a=1)")
        self.assertEqual(blocks.expect_text("walk_in_record", {}), "walk\\_in\\_record")

    def test_refusal_expected(self):
        text, prev = blocks.refusal_expected(self.scen["R-edge-0012"])
        self.assertIn("specify exactly one", text)
        same, _ = blocks.refusal_expected(self.scen["R-edge-0013"], prev)
        self.assertEqual(same, "http 400 the same text")


if __name__ == "__main__":
    unittest.main()
