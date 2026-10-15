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
# the chapters whose prose is written: their files, the overview, the glossary and their check kinds' catalogue
# entries may hold no placeholder ("to be written"); add a chapter here when its prose lands
WRITTEN_CHAPTERS = ("R",)
PLACEHOLDER = "to be written"
# Keys whose prose is still pending may carry a placeholder (the generator ships one); empty since the prose of
# the owner's review (round 2, 2026-10-15) landed. Name a file or "<chapter>.json: <key>" here while prose is awaited.
ROUND2_PENDING = ()


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
        """Every check kind any v4 scenario expects has a catalogue entry in prose/checks.json, a chapter prose
        entry or a docstring in checks.py; every catalogue entry names a check of CHECKS, has the four texts and
        states the generic flag as checks.is_property does."""
        manifest, mapping, units = v4()
        checks = cmp.load_checks_module(SUITE)
        self.assertIsNotNone(checks)
        catalogue = (sc.read_json(os.path.join(PROSE_DIR, "checks.json")) or {}).get("checks") or {}
        self.assertTrue(catalogue, "prose/checks.json has no entries")
        for name, entry in catalogue.items():
            self.assertIn(name, checks.CHECKS, "checks.json describes %s, which CHECKS lacks" % name)
            for key in ("title", "asserts", "inconclusive", "fails", "args"):
                self.assertTrue(entry.get(key), "checks.json %s lacks %s" % (name, key))
            self.assertEqual(bool(entry.get("generic")), not checks.is_property(name),
                             "checks.json %s: generic flag disagrees with checks.is_property" % name)
        prose_checks = set(k for k, v in catalogue.items() if v.get("asserts"))
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

    def test_prose_written(self):
        """The written chapters (WRITTEN_CHAPTERS) carry no placeholder: not in their intro, summary, coverage gap,
        family prose and unit lines (a unit entry may be {property, note}), not in the overview, the glossary, the
        closing templates and the case reasons, nor in the catalogue entries of the check kinds they expect; the
        keys of ROUND2_PENDING are tolerated until their prose lands."""
        manifest, mapping, units = v4()
        problems = []
        for ch in WRITTEN_CHAPTERS:
            data = sc.read_json(os.path.join(PROSE_DIR, "%s.json" % ch)) or {}
            texts = [("intro", " ".join(data.get("intro") or []))]
            for key in ("summary", "coverage_gap"):
                v = data.get(key)
                texts.append((key, " ".join(v) if isinstance(v, list) else (v or "")))
            texts += [("families.%s" % fid, (f or {}).get("prose") or "") for fid, f in (data.get("families") or {}).items()]
            for uid, t in (data.get("units") or {}).items():
                if isinstance(t, dict):
                    texts.append(("units.%s" % uid, t.get("property") or ""))
                    if t.get("note"):
                        texts.append(("units.%s.note" % uid, t["note"]))
                else:
                    texts.append(("units.%s" % uid, t or ""))
            for uid, t in ((data.get("mcp_notes") or {}).get("units") or {}).items():
                texts.append(("mcp_notes.units.%s" % uid, t or ""))
            problems += ["%s.json: %s" % (ch, key) for key, text in texts if PLACEHOLDER in text or not text.strip()]
            kinds = set(e.get("check") for s in manifest["scenarios"] if s.get("area") == ch
                        for e in s.get("expect") or [])
            catalogue = (sc.read_json(os.path.join(PROSE_DIR, "checks.json")) or {}).get("checks") or {}
            problems += ["checks.json: %s (expected by chapter %s)" % (k, ch) for k in sorted(kinds)
                         if not (catalogue.get(k) or {}).get("asserts") or PLACEHOLDER in json.dumps(catalogue[k])]
        overview = sc.read_json(os.path.join(PROSE_DIR, "overview.json")) or {}
        for key, sec in (overview.get("sections") or {}).items():
            if not sec.get("paragraphs") or PLACEHOLDER in json.dumps(sec):
                problems.append("overview.json: %s" % key)
        for key in ("appendix_F", "appendix_G"):
            if PLACEHOLDER in json.dumps(overview.get(key) or ""):
                problems.append("overview.json: %s" % key)
        glossary = sc.read_json(os.path.join(PROSE_DIR, "glossary.json")) or {}
        for e in glossary.get("entries") or []:
            text = e.get("definition") or e.get("text") or ""
            if PLACEHOLDER in text or not text.strip():
                problems.append("glossary.json: %s" % e.get("term"))
        if PLACEHOLDER in str(glossary.get("lead") or "") or "{readme_commit}" not in str(glossary.get("lead") or ""):
            problems.append("glossary.json: lead (must name the README pin {readme_commit})")
        closing = sc.read_json(os.path.join(PROSE_DIR, "closing.json")) or {}
        for key in ("intro", "observations", "caveats", "plots", "sections"):
            if not closing.get(key) or PLACEHOLDER in json.dumps(closing.get(key) or ""):
                problems.append("closing.json: %s" % key)
        cases = sc.read_json(os.path.join(PROSE_DIR, "cases.json")) or {}
        for cid, text in (cases.get("short") or {}).items():
            if PLACEHOLDER in str(text):
                problems.append("cases.json: %s" % cid)
        pending = [p for p in problems if any(p.startswith(k) for k in ROUND2_PENDING)]
        problems = [p for p in problems if p not in pending]
        self.assertEqual(problems, [], "placeholders in written prose: " + "; ".join(problems))

    def test_placeholders_declared(self):
        """Every placeholder the prose uses is one the generator fills, and overview.json declares it."""
        run_placeholders = set()
        with open(os.path.join(DOC, "sections.py"), encoding="utf-8") as fh:
            src = fh.read()
        body = src.split("def placeholders(self):", 1)[1].split("def fill(self", 1)[0]
        run_placeholders.update(re.findall(r'"([a-z_0-9]+)":', body))
        overview = sc.read_json(os.path.join(PROSE_DIR, "overview.json")) or {}
        declared = set(overview.get("placeholders") or {})
        used = set()
        for name in ("overview.json", "glossary.json"):
            data = sc.read_json(os.path.join(PROSE_DIR, name)) or {}
            data.pop("placeholders", None)
            data.pop("_about", None)
            used.update(re.findall(r"\{([a-z_0-9]+)\}", json.dumps(data)))
        self.assertEqual(used - run_placeholders, set(), "placeholders the generator does not fill")
        self.assertEqual(used - declared, set(), "placeholders overview.json does not declare")
        self.assertEqual(declared - run_placeholders, set(), "declared placeholders the generator does not fill")

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
        self.assertIn("Index of test ids", main)
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
        self.assertIn("with {\\mdseries\\upshape\\glyphfont ≤} and \\$x\\_1\\$ | pipe", chapter)
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

    def test_replicates_cell_without_evaluated_member(self):
        """P-pos-reps on mcp: one by-design and one not-implemented member, so the cell shows the standing counts
        alone, never `ok 0/0`."""
        chapter = self.tex("chapter_P.tex")
        self.assertNotIn("ok 0/0", chapter)
        st = dict((sid, cells.scenario_states(self.runobj, sid)) for sid in ("P-004", "P-005"))
        text, line2 = cells.aggregate(st, "mcp")
        self.assertEqual(text, "not exercised: 2 (1\\,\\gni{} 1\\,\\gbd{})")
        self.assertIsNone(line2)
        # one kind alone: the count and its glyph
        one = {"P-005": st["P-005"], "X": st["P-005"], "Y": st["P-005"]}
        self.assertEqual(cells.aggregate(one, "mcp")[0], "not exercised: 3\\,\\gni{}")
        # a mixed unit: the ok count over the exercised members, the rest said apart
        mixed = {"P-004": st["P-004"], "P-005": st["P-005"], "P-001": cells.scenario_states(self.runobj, "P-001")}
        self.assertTrue(cells.aggregate(mixed, "mcp")[0].startswith("ok 1/1 · not exercised: 2 ("))
        # the coverage line under the table: words, never "exceptions: none" alone for an unexercised arm
        self.assertEqual(cells.coverage_text(st, "mcp"), "not exercised (1\\,\\gni{} 1\\,\\gbd{})")
        self.assertEqual(cells.coverage_text(st, "direct-local"), "exercised 2/2, passed 2")
        self.assertIn("coverage: direct-local exercised 2/2, passed 2", chapter)
        self.assertNotIn("exceptions: none", chapter)

    def test_footnote_on_offending_cells_only(self):
        """P-003: direct-remote is class same_checks (not clean), direct-local and mcp are clean: only the
        direct-remote cell carries the footnote mark (I.4: on each offending cell)."""
        import sections
        st = cells.scenario_states(self.runobj, "P-003")
        self.assertFalse(cells.is_clean(st))
        self.assertIsNone(sections._cell_fn(st["mcp"], 1))
        self.assertIsNone(sections._cell_fn(st["direct-local"], 1))
        self.assertEqual(sections._cell_fn(st["direct-remote"], 1), 1)
        chapter = self.tex("chapter_P.tex")
        row = [l for l in chapter.splitlines() if l.startswith("direct-local &")][0]
        self.assertNotIn("\\fn{", row)

    def test_index_sorted(self):
        main = self.tex()
        index = main.split("Index of test ids", 1)[1].split("\\end{multicols}", 1)[0]
        positions = dict((sid, index.find(latex.tt(sid))) for sid in self.runobj.order)
        self.assertTrue(all(p >= 0 for p in positions.values()), positions)
        self.assertEqual(sorted(positions, key=positions.get), sorted(self.runobj.order))
        self.assertNotEqual(list(self.runobj.order), sorted(self.runobj.order))  # the manifest is not sorted

    def test_sessions_text(self):
        """I.5's mode row names every orchestrator session and the clock offset the comparison applied (not the
        last session's measurement; VERIFY-SAMPLE H2)."""
        import sections

        class R(object):
            run = {"orchestrator": {"sessions": [
                {"mode": "paired", "started_utc": "2026-10-14T06:10:55Z", "mex_clock_offset_s": 0.12,
                 "loop": {"state": "done", "completed": 681, "scenarios": 681, "direct_requests": 1587}},
                {"mode": "paired", "started_utc": "2026-10-14T11:07:00Z", "mex_clock_offset_s": 0.025,
                 "loop": {"state": "done", "completed": 681, "scenarios": 681, "direct_requests": 3}}]}}
            comp = {"host_cache": {"clock_offsets_s": {"direct-local": 0.12}}}
            orchestrator = {}
        text = sections._sessions_text(R())
        self.assertIn("session 1 started 2026-10-14 06:10:55, loop done, 681 of 681 scenarios, 1587 direct requests", text)
        self.assertIn("session 2 started 2026-10-14 11:07:00", text)
        self.assertIn("+0.120 s against the laptop, applied to direct-local", text)
        self.assertNotIn("+0.025", text)

    def test_gate_offline(self):
        """Without a gate.json the gate is judged offline with a policy pinning the manifest and the document
        says so (VERIFY-SAMPLE M6); without any policy it says none was evaluated."""
        import gate_policy
        policy = gate_policy.build_policy(self.runobj.manifest, table=(CASES, {"cases": 5}), manifest_file="m.json",
                                          arms=("direct-local", "mcp"), warm_every=0)
        path = os.path.join(self.tmp, "policy.json")
        with open(path, "w") as fh:
            json.dump(policy, fh)
        run = reader.Run(self.run_dir, policy=path)
        self.assertTrue(run.gate_offline)
        self.assertIsInstance(run.gate, dict)
        self.assertIn("passed", run.gate)
        out = os.path.join(self.tmp, "gated")
        rc = make_doc.main([self.run_dir, "--out", out, "--no-pdf", "--units", self.mapping_path,
                            "--prose-dir", self.prose_dir, "--policy", path])
        self.assertEqual(rc, 0)
        with open(os.path.join(out, "doc.tex"), encoding="utf-8") as fh:
            main = fh.read()
        self.assertIn("judged offline by gate.py for this document", main)
        word = "GATE %s" % ("PASS" if run.gate["passed"] else "FAIL")
        self.assertIn(word, main)
        # the cover: the gate beside the verdict, as prominent (\coverbox with the large word)
        cover = main.split("\\coverbox{Deploy gate}{", 1)[1].split("\n", 1)[0]
        self.assertIn("\\textbf{\\large %s}" % word, cover)
        self.assertIn("\\coverbox{Run verdict}{\\textbf{\\large", main)
        self.assertIn("cold / warm / aged / not reached", main)
        self.assertFalse(self.runobj.gate_offline)
        # the ungated build says so with the same prominence
        ungated = self.tex()
        self.assertIn("no policy in the suite's gate directory pins this manifest", ungated)
        cover = ungated.split("\\coverbox{Deploy gate}{", 1)[1].split("\n", 1)[0]
        self.assertIn("\\textbf{\\large Deploy gate: not evaluated}", cover)

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
        self.assertEqual(cells.fmt_t(150), "0.2 s")

    def test_fmt_delta(self):
        """A difference carries its sign; a negative delta never prints as `<1 ms` (VERIFY-SAMPLE H1)."""
        self.assertEqual(cells.fmt_delta(-6.914), "−7 ms")
        self.assertEqual(cells.fmt_delta(-146.15), "−0.1 s")
        self.assertEqual(cells.fmt_delta(183.6), "+0.2 s")
        self.assertEqual(cells.fmt_delta(0.18), "0 ms")
        self.assertEqual(cells.fmt_delta(-0.3), "0 ms")
        self.assertEqual(cells.fmt_delta(-31377), "−31 s")
        self.assertEqual(cells.fmt_delta(None), "–")

    def test_host_marks(self):
        """The i mark for a warm execution after an interrupted one; the aggregate's mark alone when every member
        shares the state, with its share otherwise, none when all are host-warm."""
        self.assertEqual(cells.host_mark({"host_cache": "warm", "host_cache_after_interrupted": True}), "i")
        self.assertEqual(cells.host_mark({"host_cache": "cold"}), "c")
        self.assertIsNone(cells.host_mark({"host_cache": "warm"}))
        self.assertIsNone(cells.host_mark({"host_cache": "not_reached"}))
        self.assertEqual(cells.host_marks_text(["cold", "cold", "cold"]), "c")
        self.assertEqual(cells.host_marks_text(["warm", "cold", "cold", "warm", "warm", "warm"]), "c 2/6")
        self.assertEqual(cells.host_marks_text(["warm_after_interrupted"]), "i")
        self.assertEqual(cells.host_marks_text(["cold", "aged"]), "c 1/2 a 1/2")
        self.assertIsNone(cells.host_marks_text(["warm", "warm"]))
        self.assertIsNone(cells.host_marks_text([]))

    def test_glyph_counts_bound(self):
        """Each count is bound to its glyph with a thin non-breaking space and the pairs are set apart."""
        class S(object):
            def __init__(self, k):
                self.glyph_key = k
        states = [{"mcp": S("pass")}] * 3 + [{"mcp": S("inconclusive")}]
        self.assertEqual(cells.glyph_counts(states, "mcp"), "3\\,\\gpass{}\\enspace 1\\,\\ginc{}")

    def test_cut_at_word_boundary(self):
        self.assertEqual(latex.cut_text("below the service's 500 ms minimum", 20), "below the service's…")
        self.assertEqual(latex.cut_text("abcdefghij klmnop qrstuvwxyz", 20), "abcdefghij klmnop…")
        self.assertEqual(latex.cut_text("short", 20), "short")
        self.assertEqual(latex.cut_text("x" * 30, 10), "xxxxxxxxx…")
        self.assertTrue(latex.esc("the answer could not be built and written within bounds", 30).endswith("…"))
        self.assertNotIn("buil…", latex.esc("the answer could not be built and written within bounds", 30))

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


# ---------------------------------------------------------------------------------------------------------------
# the owner's review of the sample (round 2): the MCP outcome classes, the adjudication, the coverage wording, the
# pagination helpers, the plots' data, the closing templates

class ReviewRoundTwoTest(unittest.TestCase):
    """Fixtures in memory where the fake run has no case (a relayed host refusal, a host state), the fake run
    where it has one."""

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
        cls.runobj = reader.Run(cls.run_dir)
        cls.units = units_mod.resolve(cls.runobj.manifest, FAKE_MAPPING)
        make_doc.scaffold_prose(cls.prose_dir, cls.units, cls.runobj)
        cls.out = os.path.join(cls.tmp, "out")
        cls.rc = make_doc.main([cls.run_dir, "--out", cls.out, "--no-pdf", "--units", cls.mapping_path,
                                "--prose-dir", cls.prose_dir])

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def tex(self, name="doc.tex"):
        with open(os.path.join(self.out, name), encoding="utf-8") as fh:
            return fh.read()

    def test_mcp_classes_on_fixtures(self):
        """The four classes of the review (point 5) follow from the evidence: the host reached or not, the error's
        source (the case's service outcome) and the status; the duplicate-label case is not a refusal."""
        refusal = {"expect": [{"check": "status", "args": {"http": 400}}]}
        accepted = {"expect": [{"check": "status", "args": {"http": 200}}]}
        cases = {"c.relay": {"service_answer": {"outcome": "relayed"}},
                 "c.refuse": {"service_answer": {"outcome": "refused_at_submit"}},
                 "c.dedup": {"service_answer": {"outcome": "answered"}},
                 "c.ne": {"service_answer": {"outcome": "not_expressible", "missing_argument": "x"},
                          "agreement": "not_implemented_mcp"}}
        # not expressible: no MCP call made
        cls, ev = blocks.mcp_class([{"http": None, "outcome": "not_implemented_mcp", "host_cache": "not_reached"}],
                                   {"mcp": {"not_implemented": "support"}}, {"class": "not_implemented_mcp"}, cases)
        self.assertEqual(cls, "not_expressible")
        self.assertIn("support", ev)
        # refused before host execution: the service answered 400 at submit, the host not reached
        cls, ev = blocks.mcp_class([{"http": 400, "error": {"code": "resolve_query_too_long"}, "host_cache": "not_reached",
                                     "outcome": "service_by_design", "service_by_design_case": "c.refuse"}],
                                   refusal, {"class": "service_by_design", "case": "c.refuse"}, cases)
        self.assertEqual(cls, "refused_before_host")
        self.assertIn("\\texttt{resolve\\_\\allowbreak{}query", ev)
        self.assertIn("at submit", ev)
        # the same evidence under a not-expressible table case keeps the evidence and names the table
        cls, ev = blocks.mcp_class([{"http": 400, "error": {"code": "x"}, "host_cache": "not_reached",
                                     "outcome": "service_by_design", "service_by_design_case": "c.refuse"}],
                                   refusal, {"class": "not_implemented_mcp", "case": "c.ne"}, cases)
        self.assertEqual(cls, "refused_before_host")
        self.assertIn("not expressible", ev)
        # host refusal relayed through the service: the case says relayed (the host's 400 under the service's code)
        cls, ev = blocks.mcp_class([{"http": 502, "error": {"code": "backend_fault", "text": "HTTP 400: Query shorter than k"},
                                     "host_cache": "not_reached", "outcome": "service_by_design",
                                     "service_by_design_case": "c.relay"}],
                                   refusal, {"class": "service_by_design", "case": "c.relay"}, cases)
        self.assertEqual(cls, "relayed")
        # answered after service normalization: the host expects a refusal, the service reshaped and the host answered
        cls, ev = blocks.mcp_class([{"http": 200, "error": None, "host_cache": "cold", "outcome": "service_by_design",
                                     "service_by_design_case": "c.dedup"}],
                                   refusal, {"class": "service_by_design", "case": "c.dedup"}, cases)
        self.assertEqual(cls, "normalised")
        # answered by the host where the host expects 200: plain answered, never "normalised" because a job case says answered
        cls, ev = blocks.mcp_class([{"http": 200, "error": None, "host_cache": "warm", "outcome": "pass",
                                     "service_by_design_case": "c.dedup"}], accepted, {"class": "agree"}, cases)
        self.assertEqual(cls, "answered")
        # a host state: the service's own limit ended the call
        cls, ev = blocks.mcp_class([{"http": 504, "error": {"code": "resolve_timeout"}, "host_cache": "cold",
                                     "outcome": "inconclusive", "inconclusive_reason": "resolve_timeout"}],
                                   accepted, {"class": "agree"}, cases)
        self.assertEqual(cls, "host_state")
        self.assertEqual(blocks.mcp_class([], accepted, {}, cases)[0], "none")

    def test_mcp_classes_on_the_fake_run(self):
        """P-003 (answered by the service where the host refuses) reads normalised with a note under its row,
        P-004 (refused at submit) refused before host execution, P-005 (no call) not expressible; the edge
        tables label the expectation column 'host expectation' and carry the edge legend."""
        run = self.runobj
        for sid, expected in (("P-003", "normalised"), ("P-004", "refused_before_host"), ("P-005", "not_expressible"),
                              ("S-001", "refused_before_host")):
            st = cells.scenario_states(run, sid)["mcp"]
            cls, _ = blocks.mcp_class(st.lines, run.scenario(sid), st.agreement, run.cases)
            self.assertEqual(cls, expected, sid)
        chapter = self.tex("chapter_P.tex")
        self.assertIn("host expectation", chapter)
        self.assertIn("refused before host execution · not expressible through MCP · answered after service "
                      "normalization · host refusal relayed through the service", chapter)
        # P-002: the MCP answered 200 where the host refuses and its checks failed: a disagreement, never a
        # "normalised" answer (no by-design case explains it) and no explanation row
        row = [l for l in chapter.splitlines() if l.startswith("\\hyperref[detail:P-002]")][0]
        self.assertIn("\\gfail\\gne\\fn{2}\\ answered by the host through the service", row)
        self.assertNotIn("not a service refusal", chapter)
        st = cells.scenario_states(run, "P-002")["mcp"]
        self.assertEqual(blocks.mcp_class(st.lines, run.scenario("P-002"), st.agreement, run.cases)[0], "answered")
        # S-001 (the service's own refusal at submit) in a list table: the class and the service's code
        s_chapter = self.tex("chapter_S.tex")
        self.assertIn("refused before host execution\\newline http 200 \\texttt{rate\\_\\allowbreak{}limited}", s_chapter)
        # M-001 (a 400 the service made up, no case matched: a failure) and M-002 (no MCP line) in a list table
        m_chapter = self.tex("chapter_M.tex")
        self.assertIn("refused by the service (no by-design case explains it)\\newline http 400 "
                      "\\texttt{pattern\\_\\allowbreak{}unsupported}", m_chapter)
        self.assertNotIn("relayed through the service\\newline", m_chapter)
        self.assertIn("\\gnone\\ no execution", m_chapter)
        st = cells.scenario_states(run, "M-001")["mcp"]
        self.assertEqual(blocks.mcp_class(st.lines, run.scenario("M-001"), st.agreement, run.cases)[0], "refused_unexplained")

    def test_property_tested_and_adjudication(self):
        """One execution row reads 'property tested?' from its property checks; the adjudication names the
        classification first (a host state, an allowed 503, a finding)."""
        is_prop = lambda n: n not in ("status", "invariants", "budget_respected")
        deadline = {"http": 503, "error": {"code": "deadline", "text": "t"}, "outcome": "inconclusive", "wall_ms": 24281,
                    "checks": [{"check": "status", "outcome": "inconclusive", "property": False},
                               {"check": "resolve_runs", "outcome": "inconclusive", "property": True}]}
        self.assertEqual(cells.property_tested(deadline), ("no", 0, 1))
        passed = {"outcome": "pass", "checks": [{"check": "status", "outcome": "pass"},
                                                {"check": "resolve_runs", "outcome": "pass"},
                                                {"check": "resolve_graph_runs", "outcome": "fail"}]}
        self.assertEqual(cells.property_tested(passed, is_prop), ("partly", 1, 2))
        self.assertEqual(cells.property_tested({"checks": [{"check": "status", "outcome": "pass", "property": False}]}),
                         ("n/a", 1, 1))

        class R(object):
            def budget_ms(self, sid):
                return 20000

            def scenario(self, sid):
                return {}

        text = blocks.adjudication(deadline, R(), "x", False)
        self.assertTrue(text.startswith("inconclusive: an allowed 503 deadline at 24 s against a 20,000 ms budget"), text)
        self.assertIn("the property not tested", text)
        host_state = {"http": 504, "error": {"code": "resolve_timeout", "text": "t"}, "outcome": "inconclusive",
                      "inconclusive_reason": "resolve_timeout", "wall_ms": 31000,
                      "checks": [{"check": "status", "outcome": "fail", "detail": "http 504, expected 200"},
                                 {"check": "resolve_runs", "outcome": "fail", "detail": "no answer to check: http 504"}]}
        text = blocks.adjudication(host_state, R(), "x", True)
        self.assertTrue(text.startswith("inconclusive: the service's \\texttt{resolve"), text)
        self.assertIn("classified as a host-state outcome", text)
        self.assertIn("untested on this repeat", text)
        finding = {"http": 503, "error": {"code": "deadline", "text": "t"}, "outcome": "finding", "wall_ms": 9130,
                   "checks": [{"check": "budget_respected", "outcome": "finding",
                               "detail": "503 deadline delivered 6130 ms after the 3000 ms budget (client wall 9130 ms): the deadline is read"}]}
        text = blocks.adjudication(finding, R(), "x", False)
        self.assertTrue(text.startswith("\\textbf{FINDING}:"), text)
        self.assertIn("later than max(1 s, the budget)", text)
        # the shared reason of several checks is said once, never "as above"
        groups = blocks._check_groups(host_state, True)
        self.assertEqual(groups[0][1], ["status"])
        self.assertEqual(groups[1][1], ["resolve_runs"])

    def test_detail_block_shape(self):
        """R-001 (a failure on direct-local, a disagreement): the block leads with the adjudication, holds one row
        per arm and repeat, names its artifacts (the raw files and the executions file's lines) and never says
        'as above'; a FINDING line is set apart from the rows."""
        chapter = self.tex("chapter_R.tex")
        block = chapter.split("\\label{detail:R-001}", 1)[1].split("\\end{detail", 1)[0]
        self.assertIn("\\bline{adjudication}{\\textbf{direct-local}: \\textbf{FAIL}:", block)
        self.assertIn("\\bcont{\\textbf{direct-remote}: pass", block)
        self.assertIn("arm / repeat & final outcome & http · code & wall & requested budget & property tested?", block)
        self.assertIn("direct-local & \\gfail\\ \\textbf{FAIL} & 200 & 20 ms (warm)", block)
        self.assertIn("\\bline{artifacts}{\\texttt{arms/", block)
        self.assertIn("R-\\allowbreak{}001.\\allowbreak{}direct-\\allowbreak{}local.\\allowbreak{}0.\\allowbreak{}json.\\allowbreak{}gz}", block)
        self.assertIn("\\texttt{executions.\\allowbreak{}jsonl} lines", block)
        self.assertNotIn("as above", chapter)
        self.assertIn("DISAGREES", block)
        self.assertNotIn("\\findingline", block)
        self.assertIn("title={\\fn{1}\\ \\texttt{R-\\allowbreak{}001}", block.split("\\begin{detail", 1)[1] if "\\begin{detail" in block else chapter)
        # the artifact line numbers point at the right physical lines
        run = self.runobj
        line = run.arm_lines("R-001", "direct-local")[0]
        n = run.line_no(line)
        with open(os.path.join(self.run_dir, run.executions_file), encoding="utf-8") as fh:
            physical = fh.read().splitlines()[n - 1]
        self.assertEqual(json.loads(physical)["scenario"], "R-001")

    def test_finding_line_set_apart(self):
        """A finding keeps its evidence and rule in a \\findingline, visibly apart from an allowed 503."""
        body = latex.finding_line("evidence · rule: the rule")
        self.assertEqual(body, "\\findingline{evidence · rule: the rule}")
        with open(os.path.join(DOC, "templates", "doc.sty"), encoding="utf-8") as fh:
            sty = fh.read()
        self.assertIn("\\newcommand\\findingline", sty)
        self.assertIn("FINDING", sty.split("\\newcommand\\findingline", 1)[1].split("\n", 2)[1])

    def test_pagination_helpers(self):
        """A short block is one unbreakable box, a long one breaks and repeats the id on continuation; every
        heading asks for space below it; the estimate counts wrapped lines."""
        short = latex.detail_box("R-1", 1, "t", "body", long=False)
        self.assertIn("\\begin{detailshort}[title={\\fn{1}\\ \\texttt{R-\\allowbreak{}1}\\quad t},title after break={\\texttt{R-\\allowbreak{}1}\\ (continued)}]", short)
        long = latex.detail_box("R-1", 1, "t", "body", long=True)
        self.assertIn("\\begin{detaillong}[", long)
        self.assertTrue(latex.section("X", "sec:x").startswith("\\Needspace{8\\baselineskip}"))
        self.assertTrue(latex.section("X", "sec:x", level="subsection", needspace=10).startswith("\\Needspace{10\\baselineskip}"))
        self.assertEqual(blocks.estimate_lines(["a" * 10, "b" * 300]), 1 + 3)
        self.assertEqual(blocks.estimate_lines(["\\bline{k}{%s}" % ("x" * 140)]), 1)
        self.assertLessEqual(blocks.SHORT_BLOCK_LINES, 24)
        self.assertIn("\\Needspace{9\\baselineskip}", latex.plain_table("@{}l@{}", ["a"], [["b"]], needspace=9))
        with open(os.path.join(DOC, "templates", "doc.sty"), encoding="utf-8") as fh:
            sty = fh.read()
        self.assertIn("\\newtcolorbox{detailshort}[1][]{detailstyle,unbreakable,#1}", sty)
        self.assertIn("\\newtcolorbox{detaillong}[1][]{detailstyle,breakable,#1}", sty)
        self.assertIn("\\newcommand\\parhead[1]{\\par\\medskip\\Needspace{5\\baselineskip}", sty)
        main = self.tex()
        self.assertIn("\\parhead{Principal findings}", main)
        self.assertNotIn("\\paragraph{", main)

    def test_plot_data_files(self):
        """The plots are built from data files the generator writes (deterministic): the paired walls per arm
        pair and the phase quantiles, read by pgfplots by column name; a second build writes the same bytes."""
        for name in ("plot_paired_dl_mcp.csv", "plot_paired_dl_dr.csv", "plot_phases.csv"):
            self.assertTrue(os.path.exists(os.path.join(self.out, name)), name)
        with open(os.path.join(self.out, "plot_paired_dl_mcp.csv"), encoding="utf-8") as fh:
            lines = fh.read().splitlines()
        self.assertEqual(lines[0], "scenario,chapter,dl_ms,mcp_ms")
        self.assertTrue(len(lines) >= 2, lines)
        sid, ch, dl, mc = lines[1].split(",")
        self.assertEqual(ch, sid.split("-")[0])
        float(dl), float(mc)
        self.assertEqual(lines, sorted(lines[:1]) + sorted(lines[1:], key=lambda l: self.runobj.order.index(l.split(",")[0])))
        with open(os.path.join(self.out, "plot_phases.csv"), encoding="utf-8") as fh:
            phases = fh.read().splitlines()
        self.assertEqual(phases[0], "phase,n,p10_ms,p25_ms,p50_ms,p75_ms,p90_ms")
        main = self.tex()
        self.assertIn("table[x=dl_ms,y=mcp_ms,col sep=comma]{plot_paired_dl_mcp.csv}", main)
        self.assertNotIn("stacked", main)
        out2 = os.path.join(self.tmp, "out2")
        rc = make_doc.main([self.run_dir, "--out", out2, "--no-pdf", "--units", self.mapping_path,
                            "--prose-dir", self.prose_dir])
        self.assertEqual(rc, 0)
        for name in ("doc.tex", "chapter_R.tex", "plot_paired_dl_mcp.csv", "plot_phases.csv"):
            with open(os.path.join(self.out, name), "rb") as a, open(os.path.join(out2, name), "rb") as b:
                self.assertEqual(a.read(), b.read(), name)

    def test_closing_values_and_templates(self):
        """Every placeholder of closing.json's templates and captions is one the generator computes; a template
        whose values the run lacks prints a 'not computed' note, never a hole."""
        import sections
        closing = sc.read_json(os.path.join(PROSE_DIR, "closing.json")) or {}
        ctx = sections.Context(self.runobj, self.units, {}, {}, {}, ["P"], closing=closing)
        values = sections.closing_values(ctx)
        used = set(re.findall(r"\{([a-z_0-9]+)\}", json.dumps(closing.get("observations")) + json.dumps(closing.get("plots"))
                              + json.dumps(closing.get("intro")) + json.dumps(closing.get("caveats"))
                              + json.dumps(closing.get("sections"))))
        self.assertEqual(used - set(values), set(), "placeholders closing_values does not compute")
        self.assertIn("mcp_phase_largest", values)
        self.assertIn("load_caveat_line", values)
        self.assertEqual(sections.fmt_ratio(275.04), "275")
        self.assertEqual(sections.fmt_ratio(0.96), "1.0")
        self.assertEqual(sections.fmt_ratio(9.96), "10.0")
        # the table sections' prose is printed above the tables, filled; without it the generator's line stands
        sec_ctx = sections.Context(self.runobj, self.units, {}, {}, {}, ["P"],
                                   closing={"sections": {"C4": {"paragraphs": ["front n {front_net_n}"]}}})
        text = sections.closing_section_prose(sec_ctx, "C4", "fallback")
        self.assertTrue(text.startswith("front n "), text)
        self.assertNotIn("{", text)
        self.assertEqual(sections.closing_section_prose(sec_ctx, "C7", "fallback"), "fallback")
        ctx.closing = {"observations": [{"id": "x", "template": "n = {paired_dl_mcp_n}"},
                                        {"id": "y", "template": "{no_such_value}"}]}
        obs = sections.observations(ctx)
        self.assertEqual(obs[0], "n = %s" % values["paired_dl_mcp_n"])
        self.assertIn("not computed for this run: no data for no\\_such\\_value", obs[1])
        main = self.tex()
        self.assertIn("C.1\\qd Observations", main)
        self.assertIn("recorded heuristic", main)

    def test_untested_properties(self):
        """The untested-properties computation: a unit whose property check no member passed on any arm, with
        the reasons per arm; on the fake run every property was passed somewhere."""
        import sections
        ctx = sections.Context(self.runobj, self.units, {}, {}, {}, list(self.units.chapters))
        for ch in ctx.chapters:
            ctx.footnotes[ch] = cells.Footnotes()
        u = sections.untested_properties(ctx)
        self.assertEqual(u["per_unit"], [])
        self.assertTrue(all(d["untested"] == 0 for d in u["by_chapter"].values()))
        self.assertIn("0 of", sections.untested_context_line(ctx))
        st = cells.scenario_states(self.runobj, "P-005")["mcp"]
        self.assertEqual(sections._arm_reason(st.lines, st), "not exercised via MCP")
        line = {"http": 503, "error": {"code": "deadline"}, "outcome": "inconclusive"}

        class S(object):
            glyph_key = "inconclusive"
        self.assertEqual(sections._arm_reason([line], S()), "503 deadline (allowed)")
        main = self.tex()
        self.assertIn("\\parhead{Properties still untested}", main)
        self.assertIn("Every property check of the built chapters was passed", main)

    def test_timing_legends_and_record_key(self):
        """Every unit table carries the legend of its type (what the slash means there); the record key is built
        from the manifest's ground truth."""
        import sections
        p = self.tex("chapter_P.tex")
        self.assertIn(cells.TIMING_LEGEND["series"], p)
        self.assertIn(cells.TIMING_LEGEND["replicates"], p)
        self.assertIn(sections.EDGE_LEGEND, p)
        self.assertIn(cells.TIMING_LEGEND["list"], self.tex("chapter_M.tex"))
        self.assertIn(cells.TIMING_LEGEND["pairs"], self.tex("chapter_R.tex"))
        # a summary line names its slash itself (wall p50 / p90), the reading guide's prose states the three meanings
        self.assertIn("wall p50 / p90 host-cold", p)
        self.assertIn("\\bline{aggregates}", self.tex())
        self.assertIn("\\bline{edge tables}", self.tex())

        class R(object):
            def records(self):
                return [{"key": "eco", "label": "Escherichia coli K-12 MG1655", "accession": "NC_000913.3",
                         "taxid": "511145", "kingdom": "bacteria", "length": 4641652, "records_in_column": 93,
                         "transl_table": 11, "last_in_column": False},
                        {"key": "pfa", "label": "Plasmodium falciparum 3D7 chr 5", "accession": "NC_004326.2",
                         "taxid": "36329", "kingdom": "protists", "length": 1343557}]

        class C(object):
            run = R()
        compact = sections.record_key(C())
        self.assertIn("\\texttt{eco}", compact)
        self.assertIn("\\texttt{pfa}", compact)
        self.assertIn("record key: the abbreviations the matrices use", compact)
        full = sections.record_key(C(), full=True)
        self.assertIn("4,641,652", full)
        self.assertIn("protists", full)
        self.assertEqual(sections.record_key(type("X", (), {"run": type("Y", (), {"records": lambda self: []})()})()), "")

    def test_decision_slots_and_placeholders(self):
        """I.1 prints the prose slots between its computed parts (scope, the verdict lines, `verdict`, the findings
        head, `findings`, the lists, the untested head, `untested`, the list, `closing`); I.4 prints `coverage`
        before the per-arm outcome table and `exemptions` before the limitation counts; the review's placeholders
        are plain text filled from the run."""
        import sections
        overview = {"sections": {"decision": {"paragraphs": ["SCOPE {scope_note}"], "verdict": ["VERDICT-SLOT"],
                                              "findings": ["FINDINGS-SLOT {n_findings}"],
                                              "untested": ["UNTESTED-SLOT {inconclusive_by_reason}"],
                                              "closing": ["CLOSING-SLOT"]},
                                 "provenance": {"paragraphs": ["LEAD"], "coverage": ["COVERAGE-SLOT {n_executions}"],
                                                "exemptions": ["EXEMPTIONS-SLOT {declared_file}"]}}}
        ctx = sections.Context(self.runobj, self.units, {}, overview, {}, list(self.units.chapters), sample=True)
        for ch in ctx.chapters:
            ctx.footnotes[ch] = cells.Footnotes()
        text = sections.decision_summary(ctx)
        order = ["SCOPE Chapter R sample; whole-run summary shown for context.", "\\bline{run verdict}", "VERDICT-SLOT",
                 "\\parhead{Principal findings}", "FINDINGS-SLOT 0", "\\parhead{Properties still untested}",
                 "UNTESTED-SLOT", "CLOSING-SLOT"]
        positions = [text.find(s) for s in order]
        self.assertTrue(all(p >= 0 for p in positions), list(zip(order, positions)))
        self.assertEqual(positions, sorted(positions))
        self.assertNotIn("cold / warm / aged / not reached", text)  # the outcome table is I.4's now
        prov = sections.provenance_section(ctx)
        self.assertLess(prov.find("COVERAGE-SLOT"), prov.find("cold / warm / aged / not reached"))
        self.assertLess(prov.find("\\parhead{Limitations and exemptions}"), prov.find("EXEMPTIONS-SLOT none"))
        values = ctx.placeholders()
        self.assertTrue(values["gate_text"].startswith("not evaluated (no gate.json in the run"), values["gate_text"])
        self.assertNotIn("\\", values["gate_text"])
        self.assertEqual(values["verdict_text"], self.runobj.verdict)
        self.assertEqual(values["scope_note"], "Chapter R sample; whole-run summary shown for context.")
        self.assertEqual(sections.scope_note_text(sections.Context(self.runobj, self.units, {}, {}, {}, ["P", "R"])),
                         "Chapters P, R only; whole-run summary shown for context.")
        self.assertEqual(sections.scope_note_text(sections.Context(self.runobj, self.units, {}, {}, {},
                                                                   list(self.units.chapters))), "")
        self.assertIn("samples", values["load_caveat_line"])
        self.assertTrue(values["checks_commit"].endswith("as recorded)") or values["checks_commit"] == "not recorded",
                        values["checks_commit"])
        self.assertEqual(values["declared_file"], "none")
        self.assertEqual(values["n_executions"], "{:,}".format(len(self.runobj.rows)))
        self.assertIn("by the checks", values["inconclusive_by_reason"])
        # the unit note: nothing while some member was exercised via MCP, the short reason when none was
        self.assertEqual(sections.mcp_unit_note(ctx, "P", "P-pos-reps", ["P-004", "P-005"]), "")
        note = sections.mcp_unit_note(ctx, "P", "P-pos-reps", ["P-005"])
        self.assertTrue(note.startswith("\\noindent{\\sffamily\\footnotesize\\gni\\ Not exercised via MCP: "), note)
        self.assertTrue(note.rstrip().endswith("; see B.}\\par\\smallskip"), note)
        ctx.prose = {"P": {"mcp_notes": {"units": {"P-pos-reps": "Not exercised via MCP: *why*; see B."}}}}
        self.assertIn("\\emph{why}", sections.mcp_unit_note(ctx, "P", "P-pos-reps", ["P-005"]))
        # the appendix lead-ins
        ctx.overview = {"appendix_F": {"paragraphs": ["F-LEAD {run_id}"]}}
        self.assertIn("F-LEAD %s" % latex.esc(self.runobj.run_id()), sections.appendix_lead(ctx, "appendix_F"))
        self.assertEqual(sections.appendix_lead(ctx, "appendix_G"), "")

    def test_glossary_pin_and_sample_note(self):
        """The glossary names the README path and a commit; the sample build states its scope on the cover."""
        import sections
        main = self.tex()
        self.assertIn("pinned to commit", main)
        # the real glossary's lead carries the pin and the record key its lead
        glossary = sc.read_json(os.path.join(PROSE_DIR, "glossary.json")) or {}
        ctx = sections.Context(self.runobj, self.units, {}, {}, glossary, ["P"])
        text = sections.glossary(ctx)
        self.assertIn("as of commit", text)
        self.assertIn("\\parhead{Record key}", text)
        self.assertIn(latex.inline(glossary["records_key"]["lead"]), text)
        self.assertNotIn("{readme_commit}", text)
        self.assertIn("\\texttt{metagraph/\\allowbreak{}scripts/\\allowbreak{}suite/\\allowbreak{}README.\\allowbreak{}md}", main)
        self.assertNotIn("\\samplenote", main)  # every chapter built: no note
        ctx = sections.Context(self.runobj, self.units, {}, {}, {}, ["P"], sample=True)
        self.assertEqual(sections.sample_note(ctx), "\\samplenote{Chapter R sample; whole-run summary shown for context}\n")
        ctx = sections.Context(self.runobj, self.units, {}, {}, {}, ["P"])
        self.assertIn("Chapter P only; whole-run summary shown for context", sections.sample_note(ctx))
        self.assertIn("chapters built", sections.front_matter(ctx))


if __name__ == "__main__":
    unittest.main()
