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
        # below half a millisecond the sign and one decimal (VERIFY-R2 L1: never an unsigned "0 ms" for +0.18)
        self.assertEqual(cells.fmt_delta(0.18), "+0.2 ms")
        self.assertEqual(cells.fmt_delta(-0.3), "−0.3 ms")
        self.assertEqual(cells.fmt_delta(0), "0 ms")
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
        """R-001 (a failure on direct-local, a disagreement): the block leads with the adjudication of what the
        rows do not say (the failure, the disagreement; a pass has none), holds one row per arm and repeat with
        the outcome's qualifier, names its artifacts (the raw files and the executions file's lines) and never
        says 'as above'; a FINDING line is set apart from the rows."""
        chapter = self.tex("chapter_R.tex")
        block = chapter.split("\\label{detail:R-001}", 1)[1].split("\\end{detail", 1)[0]
        self.assertIn("\\bline{adjudication}{\\textbf{direct-local}: \\textbf{FAIL}:", block)
        self.assertIn("\\bcont{\\textbf{direct-remote}: DISAGREES with the reference arm}", block)
        self.assertNotIn("direct-remote}: pass", block)  # the row says it (VERIFY-R2 M2)
        self.assertIn("arm / repeat & final outcome & http · code & wall & requested budget & property tested?", block)
        self.assertIn("direct-local & \\gfail\\ \\textbf{FAIL} & 200 & 20 ms (warm)", block)
        self.assertIn("direct-remote & \\gpass\\ pass & 200 & 0.1 s (warm) & none & yes (1 of 1)", block)
        self.assertIn("mcp r0 & \\ginc\\ inconclusive (by the checks) & 200", block)
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
        self.assertLessEqual(blocks.SHORT_BLOCK_LINES, 30)  # a block under about 30 lines moves whole (VERIFY-R2 L8)
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
        self.assertEqual(phases[0], "phase,n,n_zero,p10_ms,p25_ms,p50_ms,p75_ms,p90_ms")
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
        self.assertIn("\\parhead{Record key}\\phantomsection\\label{rk:D}", text)
        self.assertIn(latex.inline(glossary["records_key"]["lead"]), text)
        # N-1: the record key before the glossary's entries
        self.assertLess(text.find("\\parhead{Record key}"), text.find("\\parhead{Glossary}"))
        self.assertIn("\\Needspace{56\\baselineskip}\n\\parhead{Glossary}", text)  # the glossary stays one page
        self.assertLess(text.find("\\parhead{Glossary}"), text.find("\\begin{multicols}"))
        self.assertNotIn("{readme_commit}", text)
        self.assertIn("\\texttt{metagraph/\\allowbreak{}scripts/\\allowbreak{}suite/\\allowbreak{}README.\\allowbreak{}md}", main)
        self.assertNotIn("\\samplenote", main)  # every chapter built: no note
        ctx = sections.Context(self.runobj, self.units, {}, {}, {}, ["P"], sample=True)
        self.assertEqual(sections.sample_note(ctx), "\\samplenote{Chapter R sample; whole-run summary shown for context}\n")
        ctx = sections.Context(self.runobj, self.units, {}, {}, {}, ["P"])
        self.assertIn("Chapter P only; whole-run summary shown for context", sections.sample_note(ctx))
        self.assertIn("chapters built", sections.front_matter(ctx))


# ---------------------------------------------------------------------------------------------------------------
# the verification of round 2 (VERIFY-R2): the gate's own counts, the replayed by-design table, the closing
# figures, the shorter detail blocks, the family counts from the classifier, the repeats in appendix F

class VerifyRoundTwoTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.run_dir = os.path.join(cls.tmp, "run")
        os.makedirs(cls.run_dir)
        build_run(cls.run_dir)
        merge_run.merge(cls.run_dir)
        sc.write_json_atomic(os.path.join(cls.run_dir, "compare.json"), cmp.compare(cls.run_dir, (CASES, {"cases": 5})))
        cls.runobj = reader.Run(cls.run_dir)
        cls.units = units_mod.resolve(cls.runobj.manifest, FAKE_MAPPING)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def ctx(self, **kw):
        import sections
        ctx = sections.Context(self.runobj, self.units, {}, {}, {}, list(self.units.chapters), **kw)
        for ch in ctx.chapters:
            ctx.footnotes[ch] = cells.Footnotes()
        return ctx

    def test_gate_outcome_counts(self):
        """H1: what blocks the gate is the gate's outcome problems (a gated arm, the policy's selection, no
        exemption), never the run's inconclusive count; the clause says so, and says when no gate was judged."""
        import sections

        class R(object):
            gate = {"problems": [
                {"kind": "outcome", "scenario": "A-1", "arm": "mcp", "repeat": 0, "detail": "inconclusive (checks)"},
                {"kind": "outcome", "scenario": "A-1", "arm": "mcp", "repeat": 1, "detail": "inconclusive (checks)"},
                {"kind": "outcome", "scenario": "B-1", "arm": "direct-local", "repeat": 0, "detail": "finding"},
                {"kind": "relation", "scenario": None, "arm": "direct-local", "detail": "REL-x: inconclusive"},
                {"kind": "configuration", "arm": "mcp", "detail": "no live document"}]}
        gp = sections.gate_outcome_counts(R())
        self.assertEqual(gp, {"judged": True, "outcome": 3, "inconclusive": 2, "finding": 1, "pairs": 2})
        clause = sections.gate_inconclusive_clause(R())
        self.assertTrue(clause.startswith("2 of them block the deploy gate: those on a gated arm, in the policy's selection"), clause)
        self.assertIn("the others are on an ungated arm, outside the selection or exempt", clause)

        class N(object):
            gate = None
        self.assertEqual(sections.gate_outcome_counts(N())["judged"], False)
        self.assertIn("was not judged", sections.gate_inconclusive_clause(N()))
        # the fake run has no gate: I.4's row and I.1's clause say so, never "every one blocks the gate"
        ctx = self.ctx()
        text = sections.limitations_table(ctx)
        self.assertIn("1; the gate was not judged; each footnoted in its chapter's tables", text)
        self.assertNotIn("every one blocks the gate", text)
        values = ctx.placeholders()
        self.assertEqual(values["n_gate_outcome_problems"], 0)
        self.assertEqual(values["untested_scope"], "in this run")
        self.assertEqual(sections.untested_scope_text(sections.Context(self.runobj, self.units, {}, {}, {}, ["P"])),
                         "in the built chapter P")
        self.assertEqual(sections.untested_scope_text(sections.Context(self.runobj, self.units, {}, {}, {}, ["P", "R"])),
                         "in the built chapters P, R")

    def test_by_design_table_is_the_replayed_one(self):
        """H2: with a re-evaluated run the document prints the table the re-evaluation replayed, found by its
        sha256 among the replayed path, the suite's copy and the snapshot; another table is read only with the
        warning kept, and I.4 / F print the sha and the warning."""
        import sections
        table_a = os.path.join(self.tmp, "a.json")
        with open(table_a, "w") as fh:
            json.dump({"version": 7, "service_commit": "abc", "cases": [{"id": "x.%d" % i} for i in range(5)]}, fh)
        run = reader.Run(self.run_dir)
        run.reeval = {"by_design": {"path": table_a, "sha256": reader.sha256_of(table_a), "cases": 5}}
        cases, meta = run._load_table(None, True)
        self.assertEqual(meta["path"], table_a)
        self.assertEqual(meta["cases"], 5)
        self.assertTrue(meta["replayed"])
        self.assertEqual(meta["warnings"], [])
        self.assertEqual(sorted(cases), ["x.0", "x.1", "x.2", "x.3", "x.4"])
        run.table_meta = meta
        self.assertIn("(the table the re-evaluation replayed)", sections.by_design_table_text(run))
        self.assertIn("5 cases", sections.by_design_table_text(run))
        # a replayed sha no file has: the comparator's table is read, the warning kept and printed
        run.reeval = {"by_design": {"path": table_a, "sha256": "0" * 64, "cases": 9}}
        cases, meta = run._load_table(None, True)
        self.assertFalse(meta["replayed"])
        self.assertEqual(meta["replayed_sha256"], "0" * 64)
        self.assertTrue(meta["warnings"], meta)
        run.table_meta = meta
        text = sections.by_design_table_text(run)
        self.assertIn("WARNING: not the table the re-evaluation replayed, sha256 000000000000…, 9 cases", text)
        self.assertTrue(len(meta["sha256"]) == 64)
        ctx = self.ctx()
        ctx.run = run
        self.assertIn("\\textbf{not} the table the re-evaluation replayed", sections.provenance_table(ctx))
        self.assertIn("\\textbf{WARNING}: not the table the re-evaluation replayed", sections.environment_table(ctx))
        self.assertIn("by\\_design", sections.environment_table(self.ctx()) + "by\\_design")
        self.assertIsNone(reader.sha256_of(os.path.join(self.tmp, "missing.json")))

    def test_closing_figures(self):
        """H3: the pairs within 1 s count |difference|; H4: the scenarios the MCP pairing leaves out are the
        scenarios less the MCP pairs (the direct pairing's count is another value); M4: the phase population is
        the summary's (one warm execution per paired scenario) with the zeros counted apart."""
        import sections
        ctx = self.ctx()
        v = sections.closing_values(ctx)
        pairs = [p for p in sections.paired_walls(ctx) if p[2] is not None and p[4] is not None]
        diffs = [p[4] - p[2] for p in pairs]
        self.assertEqual(v["paired_dl_mcp_under_1s"], str(sum(1 for d in diffs if abs(d) < 1000)))
        self.assertEqual(v["paired_dl_mcp_faster"], str(sum(1 for d in diffs if d < 0)))
        s = self.runobj.summary()
        n_local = (s.get("mcp_overhead") or {}).get("minus_local_ms", {}).get("n")
        if n_local:
            self.assertEqual(v["scenarios_without_mcp_local_pair"], str(len(self.runobj.order) - n_local))
        self.assertIn("cold_mcp_all_n", v)
        pop = sections.phase_population(ctx)
        ph = s.get("mcp_phases") or {}
        for key in ("submit_ms", "backend_ms", "paging_ms"):
            if key in ph:
                self.assertEqual(len(pop[key]), ph[key]["n"], key)
        for name, p10, p25, p50, p75, p90, n, n_zero in sections.phase_quantiles(ctx):
            key = dict((v_, k) for k, v_ in sections.PHASES)[name]
            self.assertEqual(n + n_zero, len(pop[key]))
            self.assertTrue(p10 > 0)
        self.assertEqual(v["phase_paging_zero_n"], str(sum(1 for x in pop["paging_ms"] if x == 0)) if pop["paging_ms"] else None)
        # N-2: a two-line tick label (the phase, then its n; `n 10 of 13` where zeros are left out), centred
        plot = latex.phase_boxplot([("paging", 1, 2, 3, 4, 5, 10, 3), ("submit", 1, 2, 3, 4, 5, 10, 0)], "fig:x")
        self.assertIn("{paging\\\\(n 10 of 13)}", plot)
        self.assertIn("{submit\\\\(n 10)}", plot)
        self.assertIn("xticklabel style={align=center", plot)
        self.assertNotIn("zeros left out", plot)

    def test_first_clause_balanced(self):
        """L2: the evidence clause splits at the first ': ' outside a parenthesis, and a cut never leaves a
        parenthesis open."""
        text = "http 504, expected 200 (error resolve_timeout: the backend did not answer within 30.25 s)"
        self.assertEqual(blocks._split_clause(text), (text, ""))
        self.assertEqual(blocks._split_clause("no answer to check: http 504 (x: y)"), ("no answer to check", "http 504 (x: y)"))
        self.assertEqual(blocks._first_clause(text, 160), text)
        cut = blocks._first_clause(text, 50)
        self.assertEqual(cut.count("("), cut.count(")"))
        self.assertTrue(cut.endswith("…)"), cut)
        self.assertEqual(blocks._balanced("a (b"), "a (b)")
        self.assertEqual(blocks._balanced("a (b)"), "a (b)")

    def test_finding_text_once(self):
        """M2: the FINDING line gives the evidence per execution, the mechanism once when shared, the rule once."""
        items = [("direct-local", 0, "budget_respected", "503 delivered 6130 ms after the 3000 ms budget", "the deadline is read where polled"),
                 ("direct-remote", 0, "budget_respected", "503 delivered 7452 ms after the 3000 ms budget", "the deadline is read where polled")]
        text = blocks.finding_text(items, "SPEC 4.5: a 503 long past the budget is a finding")
        self.assertEqual(text.count("the deadline is read where polled"), 1)
        self.assertEqual(text.count("rule:"), 1)
        self.assertIn("direct-local r0, budget\\_respected: 503 delivered 6130 ms", text)
        self.assertIn("direct-remote r0, budget\\_respected: 503 delivered 7452 ms", text)
        items[1] = items[1][:4] + ("another mechanism",)
        text = blocks.finding_text(items, None)
        self.assertIn("— the deadline is read where polled;", text)
        self.assertIn("— another mechanism", text)
        self.assertNotIn("rule:", text)

    def test_rows_carry_the_adjudication(self):
        """M2: an ordinary allowed 503 is one row (the qualifier in the outcome cell) and no adjudication line;
        a late 503 reads FINDING (late 503); a host state keeps its adjudication; a not-exercised MCP arm is one
        phrase, pointing at the unit's note when the table carries it."""
        is_prop = lambda n: n not in ("status", "invariants", "budget_respected")

        class R(object):
            cases = {}
            is_property = staticmethod(is_prop)

            def budget_ms(self, sid):
                return 20000

            def scenario(self, sid):
                return {"mcp": {}}

        def st(lines, arm, key=None, agreement=None):
            s_ = cells.ArmState(arm)
            s_.lines = lines
            s_.present = True
            s_.glyph_key = key or lines[0].get("outcome")
            s_.worst = s_.glyph_key
            s_.agreement = agreement or {}
            return s_
        allowed = {"http": 503, "error": {"code": "deadline", "text": "t"}, "outcome": "inconclusive", "wall_ms": 24281,
                   "repeat": 0, "host_cache": "warm", "checks": [{"check": "resolve_runs", "outcome": "inconclusive", "property": True}]}
        rows = blocks._execution_rows(st([allowed], "direct-local"), R(), "x", "direct-local", is_prop)
        self.assertEqual(rows[0][1], "\\ginc\\ inconclusive (allowed 503)")
        self.assertEqual(rows[0][5], "no")
        self.assertIsNone(blocks._arm_adjudication(st([allowed], "direct-local"), R(), "x", "direct-local"))
        late = dict(allowed, outcome="finding", checks=[{"check": "budget_respected", "outcome": "finding", "detail": "503 late: why"}])
        rows = blocks._execution_rows(st([late], "direct-remote"), R(), "x", "direct-remote", is_prop)
        self.assertEqual(rows[0][1], "\\gfind\\ \\textbf{FINDING} (late 503)")
        self.assertIsNone(blocks._arm_adjudication(st([late], "direct-remote"), R(), "x", "direct-remote"))
        passed = {"http": 200, "outcome": "pass", "wall_ms": 5, "repeat": 0, "host_cache": "warm",
                  "checks": [{"check": "resolve_runs", "outcome": "pass", "property": True}]}
        self.assertIsNone(blocks._arm_adjudication(st([passed], "direct-local"), R(), "x", "direct-local"))
        host = {"http": 504, "error": {"code": "resolve_timeout", "text": "t"}, "outcome": "inconclusive", "repeat": 0,
                "inconclusive_reason": "resolve_timeout", "wall_ms": 31000, "host_cache": "cold", "checks": []}
        s_ = st([host, dict(passed, repeat=1)], "mcp", "inconclusive")
        rows = blocks._execution_rows(s_, R(), "x", "mcp", is_prop)
        self.assertEqual(rows[0][1], "\\ginc\\ inconclusive (host state)")
        self.assertEqual(rows[1][0], "mcp r1")
        adj = blocks._arm_adjudication(s_, R(), "x", "mcp")
        self.assertTrue(adj.startswith("r0 inconclusive: the service's \\texttt{resolve"), adj)
        self.assertNotIn("r1 pass", adj)
        # the not-exercised MCP arm: one phrase; with the unit's note, no short reason repeated
        # as the record has it: the service's refusal at submit is a service_by_design execution whose case has
        # the agreement not_implemented_mcp, so the arm stands not exercised (○) with the class refused before host
        refused = {"http": 400, "error": {"code": "resolve_query_too_long"}, "outcome": "service_by_design",
                   "host_cache": "not_reached", "wall_ms": 58, "repeat": 0, "service_by_design_case": "c.ne"}
        R.cases = {"c.ne": {"service_answer": {"outcome": "not_expressible", "missing_argument": "bounds.max_query_bp"},
                            "agreement": "not_implemented_mcp"}}
        s_ = st([refused, dict(refused, repeat=1, wall_ms=53)], "mcp", "not_implemented_mcp",
                {"class": "not_implemented_mcp", "case": "c.ne"})
        short = {"c.ne": "the tool always sends its own bound"}
        text = blocks._not_exercised_text(s_, R(), "x", short, unit_note=True)
        self.assertEqual(text, "\\gni\\ not exercised via MCP: refused before host execution (http 400 \\texttt{resolve\\_\\allowbreak{}query\\_\\allowbreak{}too\\_\\allowbreak{}long}); the note under the table and case \\texttt{c.\\allowbreak{}ne} in B")
        text = blocks._not_exercised_text(s_, R(), "x", short, unit_note=False)
        self.assertTrue(text.startswith("\\gni\\ not exercised via MCP: the tool always sends its own bound; refused before host execution"), text)
        self.assertTrue(text.endswith("case \\texttt{c.\\allowbreak{}ne} in B"), text)
        rows = blocks._execution_rows(s_, R(), "x", "mcp", is_prop)
        self.assertEqual(rows[0][1], "\\gni\\ not exercised")
        self.assertEqual(rows[0][5], "–")
        self.assertEqual(rows[0][3], "58 ms")

    def test_family_values_from_the_classifier(self):
        """M1: a family's prose placeholders count the MCP classes as the cells derive them."""
        ctx = self.ctx()
        v = ctx.family_values("P-pos")
        self.assertEqual(v["n_tests"], 4)
        self.assertEqual(v["n_units"], 2)
        self.assertEqual(v["n_mcp_normalised"], 1)
        self.assertEqual(v["mcp_normalised_ids"], "P-003")
        self.assertEqual(v["n_mcp_refused_before_host"], 1)
        self.assertEqual(v["n_mcp_not_expressible"], 1)
        self.assertEqual(v["n_mcp_answered"], 1)
        self.assertEqual(v["mcp_relayed_ids"], "none")
        prose = {"P": {"families": {"P-pos": {"title": "t", "prose": "{n_mcp_normalised} normalised, {n_tests} tests"}},
                       "units": {}}}
        ctx.prose = prose
        import sections
        text = sections.family_section(ctx, "P", "P-pos")
        self.assertIn("1 normalised, 4 tests", text)
        # the refusals unit's summary line names its unit (L5)
        self.assertIn("\\sumline{%s: 1 test ·" % latex.tt("P-edge-refused"), sections.family_section(ctx, "P", "P-edge"))

    def test_gate_problems_list_repeats(self):
        """L3: appendix F's outcome row lists every scenario-arm pair with its repeats, and says how many
        executions the pairs hold when they differ."""
        import sections
        ctx = self.ctx()
        ctx.run.gate = {"passed": False, "policy": {"name": "p", "policy_version": 1, "sha256": "ab"}, "counts": {},
                        "problems": [
                            {"kind": "outcome", "scenario": "P-001", "arm": "mcp", "repeat": 0, "detail": "inconclusive (checks)"},
                            {"kind": "outcome", "scenario": "P-001", "arm": "mcp", "repeat": 1, "detail": "inconclusive (checks)"},
                            {"kind": "outcome", "scenario": "R-001", "arm": "direct-local", "repeat": 0, "detail": "fail"}],
                        "notes": []}
        ctx.run.gate_offline = True
        try:
            text = sections.gate_problems(ctx)
        finally:
            ctx.run.gate = None
            ctx.run.gate_offline = False
        self.assertIn("mcp (r0, r1)", text)
        self.assertIn("3 executions in 2 scenario-arm pairs", text)

    def test_appendix_rules(self):
        """M3: the judging, verdict and gate rules print in appendix F from `appendix_F.rules`; the real
        overview's methods hold the arms, the order, the ground truth and the times only."""
        import sections
        ctx = self.ctx()
        self.assertEqual(sections.appendix_rules(ctx), "")
        ctx.overview = {"appendix_F": {"rules_title": "The rules", "rules": ["JUDGING {n_scenarios}"]}}
        text = sections.appendix_rules(ctx)
        self.assertTrue(text.startswith("\\parhead{The rules}\n"), text)
        self.assertIn("JUDGING %d" % len(self.runobj.order), text)
        overview = sc.read_json(os.path.join(PROSE_DIR, "overview.json")) or {}
        heads = [p.split("*")[1] for p in overview["sections"]["methods"]["paragraphs"]]
        self.assertEqual(heads, ["The three arms.", "Execution order.", "Ground truth.", "Times and load."])
        self.assertEqual(len(overview["appendix_F"]["rules"]), 2)
        self.assertNotIn("*Coverage.*", json.dumps(overview["sections"]["provenance"]))


# ---------------------------------------------------------------------------------------------------------------
# the re-verification of round 3 (VERIFY-R2 "Re-verification"): the record key in I.3, before the glossary and
# linked from the legends (N-1), the box plot's labels (N-2, above), the prose filled from the run (N-3), the
# closing populations' wording (L-b), the uncut B note (L-d)

V4_LITERALS = ("76 of the 81", "78 on direct-remote", "79 identical", "27 of 28", "51 pairs", "14 by-design",
               "16 tests", "74 ms", "n 29", "11 s at p90", "18 s", "6.1 s", "7.5 s", "9.5 s", "25.5 s", "24 s", "8.2 s",
               "completed in 0.4 s", "four times", "five tests", "four tests", "three of the four", "kept 1 of 9",
               "7 of the 8", "70, 70, 71 and 76", "98, 99, 100 and 100", "between 1 and 25 s", "seven edge inputs",
               "six of them")
PROSE_FILES = ("R.json", "overview.json", "glossary.json", "closing.json", "cases.json", "checks.json")
PROSE_SKIP_KEYS = ("_about", "placeholders", "_assumes", "_placeholders")
# what is not a figure of a run: placeholders, code spans, ids, dates, SPEC and section references, p50 / r0
# names, ordinals, appendix letters (the conditionals' heads go first, through strip_conditional_heads)
_NOT_A_FIGURE = [r"\{[^{}]*\}", r"`[^`]*`", r"\b[A-Z]+-[a-z_]+-\d{4}\b", r"(?<=-)\d{4}\b",
                 r"\d{4}-\d{2}-\d{2}", r"\bSPEC(?:-[A-Za-z-]+)?(?:\.md)?\s*(?:§\s*)?(?:\d+(?:\.\d+)?)?",
                 r"\bsections?\s+\d+(?:\.\d+)?(?:(?:,\s*|\s+and\s+|\s+to\s+)\d+(?:\.\d+)?)*", r"§\s*\d+(?:\.\d+)?",
                 r"\b[A-Z]\.\d\b", r"\b[pr]\d+\b", r"\b\d+(?:st|nd|rd|th)\b", r"\bappendi(?:x|ces)\s+[A-G]\b"]


def prose_strings(obj, path=""):
    """(path, text) of every prose string of a prose file (its documentation keys left out)."""
    if isinstance(obj, dict):
        for k, v in obj.items():
            if k in PROSE_SKIP_KEYS:
                continue
            for x in prose_strings(v, path + "." + k):
                yield x
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            for x in prose_strings(v, "%s[%d]" % (path, i)):
                yield x
    elif isinstance(obj, str):
        yield path, obj


_COND = re.compile(r"\[\[((?:(?!\[\[|\]\]).)*)\]\]", re.S)


def strip_conditional_heads(text):
    """(the text with every conditional's head dropped and its branches kept, the heads): innermost first, as
    facts.fill reads them; a bracketed run list such as `[[0, 1]]` is prose and stays."""
    heads = []
    pos = 0
    while True:
        m = _COND.search(text, pos)
        if not m:
            return text, heads
        head, sep, rest = m.group(1).partition("?")
        if not sep:
            pos = m.end()
            continue
        heads.append(head.strip())
        text = text[:m.start()] + " " + rest.replace("|", " ") + " " + text[m.end():]
        pos = 0


def literal_numbers(text):
    """The numbers a prose string carries as text (not as a placeholder, a flag, an id or a reference)."""
    text, _ = strip_conditional_heads(text)
    for pat in _NOT_A_FIGURE:
        text = re.sub(pat, " ", text)
    return [n.rstrip(",.") for n in re.findall(r"(?<![A-Za-z_\d])\d[\d,]*(?:\.\d+)?(?:\.\.\d+)?", text)]


class VerifyRoundThreeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.run_dir = os.path.join(cls.tmp, "run")
        os.makedirs(cls.run_dir)
        build_run(cls.run_dir)
        merge_run.merge(cls.run_dir)
        sc.write_json_atomic(os.path.join(cls.run_dir, "compare.json"), cmp.compare(cls.run_dir, (CASES, {"cases": 5})))
        cls.runobj = reader.Run(cls.run_dir)
        cls.units = units_mod.resolve(cls.runobj.manifest, FAKE_MAPPING)
        cls.R = sc.read_json(os.path.join(PROSE_DIR, "R.json"))

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def ctx(self, prose=None, **kw):
        import sections
        ctx = sections.Context(self.runobj, self.units, prose or {}, {}, {}, list(self.units.chapters), **kw)
        for ch in ctx.chapters:
            ctx.footnotes[ch] = cells.Footnotes()
        return ctx

    def test_fill_conditionals_and_specs(self):
        """facts.fill: conditionals (negated, nested, compared), the plural and word specs, the index syntax, the
        spacing around a chosen or empty branch, a literal run list left alone, the unresolved recorded."""
        import facts
        v = {"a": True, "z": 0, "n": 4, "m": 1, "x": {"R-1": {"dl": {"wall": "8.2 s", "class": "relayed"}}}}
        u = []
        self.assertEqual(facts.fill("agree[[a ?, all | ({n} of)]]; x", v, u), "agree, all; x")
        self.assertEqual(facts.fill("agree[[!a ?, all | ({n} of)]]; x", v, u), "agree (4 of); x")
        self.assertEqual(facts.fill("table[[z ? (u)]]; x [[z ? y]]; z", v, u), "table; x; z")
        self.assertEqual(facts.fill("host[[a ? under load]] they", v, u), "host under load they")
        self.assertEqual(facts.fill("{n} test{n:s} follow{n:v}; {m} test{m:s} follow{m:v}; {n:word} {n:Word}", v, u),
                         "4 tests follow; 1 test follows; four Four")
        self.assertEqual(facts.fill("[[x[R-1][dl][class] = relayed ? yes | no]]/[[x[R-1][dl][class] = other ? yes | no]] {x[R-1][dl][wall]}", v, u),
                         "yes/no 8.2 s")
        self.assertEqual(facts.fill("[[a ? [[z ? inner | outer]] | c]].", v, u), "outer.")
        self.assertEqual(facts.fill("runs [[0, 1]] and [[0, 20], [51, 70]] [[a ? y]]", v, u), "runs [[0, 1]] and [[0, 20], [51, 70]] y")
        self.assertEqual(u, [])
        out = facts.fill("{missing} and [[nope ? a | b]] {x[R-9][dl][wall]}", v, u)
        self.assertEqual(out, "{missing} and (UNRESOLVED CONDITION nope) {x[R-9][dl][wall]}")
        self.assertEqual(sorted(u), ["[[nope]]", "missing", "x[R-9][dl][wall]"])
        self.assertEqual(facts.number_word(0), "zero")
        self.assertEqual(facts.number_word(13), "13")

    def test_no_v4_literal_survives(self):
        """N-3: chapter R's prose rendered on the fake run (another run entirely) prints none of the v4 sample's
        figures: every figure comes from the run, every run-dependent sentence from a computed fact."""
        import sections
        prose = {"R": dict(self.R)}
        # the fake run's one R family carries, in turn, each v4 family's prose
        texts = []
        for fid, entry in self.R["families"].items():
            prose["R"]["families"] = {"R-one": dict(entry)}
            ctx = self.ctx(prose)
            # the family's prose paragraph (the tables and blocks after it hold the fake run's own times)
            texts.append((fid, sections.family_section(ctx, "R", "R-one").split("\\unithead", 1)[0]))
        prose["R"]["families"] = {"R-one": dict(self.R["families"]["R-budget"])}
        ctx = self.ctx(prose)
        whole = sections.chapter(ctx, "R")
        conclusion = whole.split("\\parhead{Conclusion}", 1)[1].split("\\begingroup", 1)[0]
        gap = whole.split("\\parhead{What this chapter leaves untested}", 1)[1].split("\\parhead{Background}", 1)[0]
        texts.append(("chapter", conclusion + gap))
        for where, text in texts:
            for lit in V4_LITERALS:
                self.assertNotIn(lit, text, "%s: the v4 literal %r survives" % (where, lit))
            # a condition on a v4 scenario or relation the fake run lacks stays visible; every other flag resolves
            for name in re.findall(r"\(UNRESOLVED CONDITION ([^)]*)\)", text):
                self.assertTrue(name.startswith(("x[", "fam[", "rel[")), "%s: unresolved %s" % (where, name))
        # the fake run's own facts are printed in their place
        chapter = whole
        self.assertIn("0 of the 1 tests pass on direct-local and 1 on direct-remote", chapter.replace("\\", ""))
        self.assertIn("1 pair agrees", chapter)
        self.assertIn("every test was exercised through the tools", chapter)
        self.assertNotIn("Of the zero relations", whole)
        self.assertIn("the deploy gate was not judged", chapter.lower())
        # the unresolved placeholders are the v4 scenarios the fake run lacks and the family the chapter prose names
        unresolved = set(n.strip("[]").lstrip("!").split("[")[0] for n in ctx.unresolved)
        self.assertEqual(unresolved, {"fam", "x"}, sorted(set(ctx.unresolved)))

    def test_prose_placeholders_resolvable_on_v4(self):
        """Every placeholder and flag of R.json names a value the generator computes: a chapter or family key, a
        family of the v4 mapping under `fam`, a v4 scenario of chapter R under `x` with an arm and an arm key."""
        import facts
        manifest, mapping, units = v4()
        ctx = self.ctx({"R": self.R})
        chapter_keys = set(facts.chapter_facts(ctx, "R"))
        family_keys = set(facts.family_facts(ctx, "R-one"))
        unit_keys = set(facts.unit_facts(ctx, "R-one"))
        arm_keys = set(facts.arm_facts(ctx, "R-001", "mcp"))
        scenario_keys = set(facts.scenario_facts(ctx, "R-001"))
        r_sids = set(units.chapter_scenarios("R"))
        r_fams = set(units.chapter_families("R"))
        problems = []

        def check(name, scope):
            root = name.split("[", 1)[0].split(":", 1)[0]
            idx = re.findall(r"\[([^\]]+)\]", name.split(":", 1)[0])
            if root == "fam":
                if scope != "chapter" or not idx or idx[0] not in r_fams:
                    problems.append("%s: fam[...] needs a v4 family of R in the chapter prose" % name)
                elif len(idx) >= 2 and idx[1] != "x" and idx[1] not in family_keys:
                    problems.append("%s: no family value %s" % (name, idx[1]))
                elif len(idx) >= 2 and idx[1] == "x" and (len(idx) < 3 or idx[2] not in r_sids):
                    problems.append("%s: no v4 scenario of R" % name)
            elif root == "x":
                if not idx or idx[0] not in r_sids:
                    problems.append("%s: no v4 scenario of R" % name)
                elif len(idx) == 2 and idx[1] not in scenario_keys:
                    problems.append("%s: no scenario value %s" % (name, idx[1]))
                elif len(idx) >= 3 and (idx[1] not in ("dl", "dr", "mcp") or idx[2] not in arm_keys):
                    problems.append("%s: no arm value %s" % (name, idx[1:3]))
            else:
                keys = {"chapter": chapter_keys, "family": family_keys, "unit": unit_keys}[scope]
                if root not in keys and root not in ctx.placeholders():
                    problems.append("%s: no %s value %s" % (name, scope, root))

        for path, text in prose_strings(self.R):
            scope = "chapter" if path.startswith((".summary", ".coverage_gap", ".intro", ".units")) else (
                "unit" if path.startswith(".mcp_notes") else "family")
            for name in re.findall(r"\{([^{}]+)\}", text):
                check(name, scope)
            for head in strip_conditional_heads(text)[1]:
                check(head.lstrip("!").split("=", 1)[0].strip(), scope)
        self.assertEqual(problems, [], "; ".join(problems))
        self.assertEqual(self.R.get("coverage_gap_families"), ["R-budget"])
        for fid in r_fams:
            self.assertTrue(self.R["families"][fid].get("claim"), "%s: no claim" % fid)
            for sid in self.R["families"][fid].get("claim_members") or []:
                self.assertIn(sid, r_sids)

    def test_prose_numbers_filled_or_assumed(self):
        """A number a prose string carries as text is a constant the file's `_assumes` lists with its source (the
        SPEC's, the manifest's, the service's); a run's figure is a placeholder, never text."""
        problems = []
        for name in PROSE_FILES:
            data = sc.read_json(os.path.join(PROSE_DIR, name)) or {}
            assumes = data.get("_assumes") or {}
            self.assertTrue(assumes.get("_note"), "%s: _assumes without its note" % name)
            for path, text in prose_strings(data):
                for n in literal_numbers(text):
                    if n not in assumes:
                        problems.append("%s%s: %s" % (name, path, n))
        self.assertEqual(problems, [], "numbers the generator did not fill: " + "; ".join(sorted(set(problems))))
        # the scan sees what it should: a figure, not an id, a reference or a placeholder
        self.assertEqual(literal_numbers("R-pos-0027 and REL-prefix_equals-0004, -0002, I.4, p90, r0, SPEC 4.2, section 4.5, "
                                         "sections 7.4 to 7.6, `max_labels 1`, {n_tests} {x[R-pos-0001][dl][wall]} "
                                         "[[fam[R-budget][n_late503] ? a | b]] appendix D 2026-10-15 90th"), [])
        self.assertEqual(literal_numbers("76 of the 81 tests, 2,000 bp, 6.1 s, 1..100, 0.2 s"), ["76", "81", "2,000", "6.1", "1..100", "0.2"])

    def test_family_and_chapter_facts(self):
        """The facts on the fake run: the standing counts, the footnoted members' phrases, the findings, the MCP
        classes and variants, the claims sorted by what the direct arms showed, the gate flag."""
        import facts
        ctx = self.ctx({"P": {"families": {"P-pos": {"title": "t", "prose": "p", "claim": "a window resolves"},
                                           "P-edge": {"title": "e", "prose": "p", "claim": "an edge is refused"}}}})
        fv = facts.family_facts(ctx, "P-pos")
        self.assertEqual(fv["n_tests"], 4)
        self.assertEqual(fv["n_mcp_normalised"], 1)
        self.assertEqual(fv["mcp_normalised_ids"], "P-003")
        self.assertEqual(fv["mcp_normalised_variants"], "20 bp")
        self.assertIn("P-001", fv["x"])
        self.assertEqual(fv["x"]["P-001"]["dl"]["standing"], "pass")
        self.assertEqual(fv["x"]["P-001"]["dl"]["phrase"], "a pass in 9 ms / 6 ms (r0 / r1)")
        self.assertEqual(fv["x"]["P-001"]["dr"]["phrase"], "both repeats a pass in 0.1 s")
        self.assertEqual(fv["x"]["P-003"]["dl"]["phrase"], "the expected refusal (http 400) in 4 ms")
        self.assertEqual(fv["n_pass_dl"], 4)
        self.assertEqual(fv["n_relations"], 0)
        self.assertEqual(fv["n_declared_relations"], 0)
        self.assertFalse(fv["stop_observed"])
        cv = facts.chapter_facts(ctx, "P")
        self.assertEqual(cv["n_tests"], 5)
        self.assertIn("a window resolves (P-pos)", cv["claims_held"])
        self.assertIn("an edge is refused (P-edge)", cv["claims_held"])
        self.assertEqual(cv["claims_unsupported"], "")
        self.assertFalse(cv["gate_judged"])
        self.assertEqual(cv["fam"]["P-pos"]["n_tests"], 4)
        self.assertIn("direct-local", cv["direct_agreement_text"] + "direct-local")
        rv = facts.chapter_facts(self.ctx({"R": {"families": {"R-one": {"claim": "the one", "prose": "p", "title": "t"}}}}), "R")
        self.assertEqual(rv["n_direct_failed_tests"], 1)
        self.assertFalse(rv["direct_clean"])
        self.assertIn("with a failure or a disagreement", rv["claims_unsupported"])
        self.assertEqual(rv["n_footnoted"], 1)
        self.assertTrue(rv["x"]["R-001"]["dl"]["phrase"].startswith("FAIL ("), rv["x"]["R-001"]["dl"]["phrase"])
        self.assertEqual(rv["claims_held"], "")

    def test_record_key_on_the_chapter_page_and_in_legends(self):
        """N-1: the chapter page carries the compact record key before its first family (one paragraph, labelled
        rk:<chapter>), every unit table whose variants name a record links to it and to appendix D's full key in
        its legend, and appendix D's key is labelled rk:D before the glossary."""
        import sections

        class R(object):
            def records(self):
                return [{"key": "eco", "label": "Escherichia coli K-12 MG1655", "accession": "NC_000913.3", "taxid": "511145"},
                        {"key": "pfa", "label": "Plasmodium falciparum 3D7 chr 5", "accession": "NC_004326.2", "taxid": "36329"}]

        class C(object):
            run = R()
        text = sections.record_key_paragraph(C())
        self.assertIn("\\label{rk:compact}", text)
        self.assertIn("\\label{rk:P}", sections.record_key_paragraph(C(), "P"))
        self.assertIn("\\texttt{eco} Escherichia coli K-12 MG1655 · \\texttt{pfa} Plasmodium falciparum 3D7 chr 5", text)
        self.assertIn("\\hyperref[rk:D]{appendix D}", text)
        ctx = self.ctx()
        # the fake run's ground truth has no records: no key, no link
        self.assertEqual(sections.record_key_paragraph(ctx), "")
        self.assertEqual(sections.record_key_link(ctx, ["P-001"]), "")
        ctx.run.manifest["ground_truth"] = {"records": [{"key": "bp", "label": "x", "accession": "a", "taxid": "1"}]}
        try:
            self.assertIn("record keys: \\hyperref[rk:P]{p.\\,\\pageref{rk:P}}", sections.record_key_link(ctx, ["P-001"], "P"))  # variant "10 bp"
            self.assertEqual(sections.record_key_link(ctx, ["R-001"], "R"), "")  # variant "the one"
            self.assertNotIn("Record key", sections.reading_guide(ctx))  # Part I stays four pages
            chapter = sections.chapter(ctx, "P")
            self.assertLess(chapter.find("\\parhead{Background}"), chapter.find("\\label{rk:P}"))
            self.assertLess(chapter.find("\\label{rk:P}"), chapter.find("\\unithead"))
            unit = sections.unit_section(ctx, "P", "P-pos-series")
            self.assertIn("in full \\hyperref[rk:D]{appendix D}", unit)
        finally:
            ctx.run.manifest["ground_truth"] = {}

    def test_closing_untimed_cold_and_uncut_note(self):
        """L-b: the cold MCP calls C.4 counts and C.3 leaves out are a value of their own; L-d: appendix B's note
        under a row is never cut."""
        import sections
        ctx = self.ctx()
        v = sections.closing_values(ctx)
        self.assertIn("cold_mcp_untimed_n", v)
        s = ctx.run.summary()
        n_all = ((((s.get("front_overhead") or {}).get("cold") or {}).get("mcp") or {}).get("wall_ms") or {}).get("n")
        n_timed = (s.get("mcp_cold_wall_ms") or {}).get("n")
        if n_all is not None and n_timed is not None:
            self.assertEqual(v["cold_mcp_untimed_n"], str(n_all - n_timed))
        closing = sc.read_json(os.path.join(PROSE_DIR, "closing.json"))
        self.assertIn("{cold_mcp_untimed_n}", json.dumps(closing["sections"]))
        self.assertIn("cold_mcp_untimed_n", closing["placeholders"])
        self.assertNotIn("the warm execution of each paired scenario", json.dumps(closing))
        with open(os.path.join(DOC, "sections.py"), encoding="utf-8") as fh:
            src = fh.read()
        body = src.split("def by_design_table(ctx):", 1)[1].split("\ndef ", 1)[0]
        self.assertNotIn("unverified_detail), 600", body)
        self.assertIn('"; ".join(unverified_detail)))', body)
        # L-a: the limitations row no longer claims every inconclusive execution is footnoted in Part II
        self.assertIn("each footnoted in its chapter's tables", sections.limitations_table(ctx))
        glossary = sc.read_json(os.path.join(PROSE_DIR, "glossary.json"))
        inc = [e for e in glossary["entries"] if e["term"] == "inconclusive"][0]
        self.assertIn("in the policy's selection and not declared exempt", inc["definition"])
        # L-e: the by-design row of F says the table was found by its sha
        ctx.run.table_meta = {"path": "/w/service_by_design.json", "version": 1, "cases": 5, "service_commit": "abc",
                              "sha256": "d" * 64, "replayed": True, "warnings": []}
        ctx.run.reeval = {"by_design": {"path": "/main/service_by_design.json", "sha256": "d" * 64}}
        try:
            text = sections.environment_table(ctx)
        finally:
            ctx.run.table_meta = {}
            ctx.run.reeval = {}
        self.assertIn("found by its sha256 at this path (the re-evaluation recorded", text)


# ---------------------------------------------------------------------------------------------------------------
# the owner's second review (OWNER-REVIEW-SAMPLE-2.md, round 5) and the low items of the re-verification of round 3

def budget_fixture(run_dir, perturbed):
    """The fake run of test_compare_report with its chapter R read as an R-budget-like family, in two versions whose
    conclusions differ (the owner's second review, point 4): as built (R-001 fails on direct-local and the arms
    disagree; no stop, no finding), or perturbed (R-001 passes on every arm, on direct-local with a stated stop, the
    direct answers identical and the MCP pair fully checked; a second test R-002, not implemented via MCP, answers a
    503 deadline 5.8 s past its 2,500 ms budget on direct-local, a finding, and an allowed 503 deadline 0.4 s past it
    on direct-remote). Built from small data here, never from the v4 run directory."""
    from test_compare_report import line, scen, T0
    manifest = build_run(run_dir)
    arms_dir = os.path.join(run_dir, "arms")
    for s in manifest["scenarios"]:
        if s["id"] == "R-001":
            s["route"] = "POST /resolve"
            s["request"] = {"bounds": {"time_budget_ms": 2500}}
            s["expect"][1]["args"]["prefix_if_stopped"] = True
    if perturbed:
        s2 = scen("R-002", "R", not_implemented="max_query_bp")
        s2["route"] = "POST /resolve"
        s2["request"] = {"bounds": {"time_budget_ms": 2500}}
        s2["expect"][1]["args"]["prefix_if_stopped"] = True
        manifest["scenarios"].append(s2)
    manifest["content_sha256"] = sc.scenarios_sha256(manifest["scenarios"])
    sc.write_json_atomic(os.path.join(arms_dir, "direct-local", "manifest.json"), manifest)
    if perturbed:
        def rewrite(arm, edit):
            path = os.path.join(arms_dir, arm, "executions.jsonl")
            rows, _ = sc.read_jsonl(path)
            rows = edit(rows)
            os.remove(path)
            for r in rows:
                sc.append_jsonl(path, r)
                raw = os.path.join(arms_dir, arm, r["raw"])
                if not os.path.exists(raw):
                    sc.write_gz_json(raw, {"answer": {"http": r["http"], "body": None, "error": r.get("error")}})

        def local(rows):
            out = []
            for r in rows:
                if r["scenario"] == "R-001":
                    r = dict(r, outcome="pass", http=200, time_limited="resolve stopped by time",
                             checks=[{"check": "status", "outcome": "pass", "detail": "ok"},
                                     {"check": "contains_source_context", "outcome": "pass", "detail": "a prefix"}])
                out.append(r)
            late = line("R-002", "direct-local", 0, T0 + 20, 8300.0, outcome="finding", http=503, digest="r4",
                        error={"code": "deadline", "text": "resolve: nothing partial is sent"},
                        checks=[{"check": "status", "outcome": "inconclusive", "detail": "a 503 deadline"},
                                {"check": "contains_source_context", "outcome": "inconclusive", "detail": "no answer"},
                                {"check": "budget_respected", "outcome": "finding",
                                 "detail": "503 deadline delivered 5800 ms after the 2500 ms budget (client wall "
                                           "8300 ms): the deadline is read where the work polls it"}])
            late["server"] = None
            out.append(late)
            return out

        def remote(rows):
            rows = [dict(r, digest="r1", answer_digest="r1") if r["scenario"] == "R-001" else r for r in rows]
            allowed = line("R-002", "direct-remote", 0, T0 + 50, 2900.0, outcome="inconclusive", http=503, digest="r4",
                           error={"code": "deadline", "text": "resolve: nothing partial is sent"},
                           checks=[{"check": "status", "outcome": "inconclusive", "detail": "a 503 deadline"},
                                   {"check": "contains_source_context", "outcome": "inconclusive", "detail": "no answer"},
                                   {"check": "budget_respected", "outcome": "pass", "detail": "within max(1 s, the budget)"}])
            allowed["server"] = None
            return rows + [allowed]

        def mcp(rows):
            rows = [dict(r, checks=[{"check": "status", "outcome": "pass", "detail": "ok"},
                                    {"check": "contains_source_context", "outcome": "pass", "detail": "found"}])
                    if r["scenario"] == "R-001" else r for r in rows]
            return rows + [line("R-002", "mcp", 0, T0 + 60, 1.0, outcome="not_implemented_mcp")]
        rewrite("direct-local", local)
        rewrite("direct-remote", remote)
        rewrite("mcp", mcp)
    merge_run.merge(run_dir)
    sc.write_json_atomic(os.path.join(run_dir, "compare.json"), cmp.compare(run_dir, (CASES, {"cases": len(CASES)})))
    mapping = json.loads(json.dumps(FAKE_MAPPING))
    del mapping["_families"]["R-one"]
    del mapping["_units"]["R-one"]
    members = ["R-001", "R-002"] if perturbed else ["R-001"]
    mapping["_families"]["R-budget"] = {"chapter": "R", "title": "the time budget", "units": ["R-budget"]}
    mapping["_units"]["R-budget"] = {"family": "R-budget", "unit_type": "pairs", "unit_title": "the stop",
                                     "property": "a stop answers the prefix", "members": members}
    for sid, variant in (("R-001", "eco, 2500 ms stop"), ("R-002", "pae, 2500 ms stop")):
        if sid in members:
            mapping[sid] = {"chapter": "R", "family": "R-budget", "unit_id": "R-budget", "unit_type": "pairs",
                            "variant": variant}
    return mapping


class OwnerReviewTwoTest(unittest.TestCase):
    """The owner's second review (round 5): the plots' size, 'fully checked' over the repeats compared, the timing
    wordings, the conclusions of two fixtures; and the re-verification's low items: the untimed cold split, the
    number-unit space, the variant labels, the unit note's facts, the record key's wording."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.R = sc.read_json(os.path.join(PROSE_DIR, "R.json"))
        cls.dirs, cls.mappings = {}, {}
        for name, perturbed in (("plain", False), ("perturbed", True)):
            d = os.path.join(cls.tmp, name)
            os.makedirs(d)
            cls.mappings[name] = budget_fixture(d, perturbed)
            cls.dirs[name] = d

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def render(self, name):
        """(ctx, the conclusion, the coverage gap, the R-budget family prose) of chapter R on a fixture, with
        R.json's chapter prose and its R-budget family prose."""
        import sections
        runobj = reader.Run(self.dirs[name])
        units = units_mod.resolve(runobj.manifest, self.mappings[name])
        prose = {"R": dict(self.R)}
        prose["R"]["families"] = {"R-budget": dict(self.R["families"]["R-budget"])}
        ctx = sections.Context(runobj, units, prose, {}, {}, list(units.chapters))
        for ch in ctx.chapters:
            ctx.footnotes[ch] = cells.Footnotes()
        text = sections.chapter(ctx, "R")
        conclusion = text.split("\\parhead{Conclusion}", 1)[1].split("\\begingroup", 1)[0]
        gap = text.split("\\parhead{What this chapter leaves untested}", 1)[1].split("\\parhead{Background}", 1)[0]
        family = text.split("fam:R-budget", 1)[1].split("\\unithead", 1)[0]
        return ctx, conclusion, gap, family

    def test_two_fixtures_change_conclusions(self):
        """The owner's second review, point 4: chapter R rendered on two fixtures that differ in their verdicts
        changes its counts, its 'What held' branch, 'Shown' against 'Not shown', the stop sentence, 'no failure, no
        error and no disagreement', 'the run's only findings' and the family's statements, not only its numbers;
        no v4 figure and no unresolved placeholder other than the v4 scenarios the fixtures lack survive."""
        ctx_a, conc_a, gap_a, fam_a = self.render("plain")
        ctx_b, conc_b, gap_b, fam_b = self.render("perturbed")
        # the counts
        self.assertIn("0 of the 1 tests pass on direct-local and 1 on direct-remote", conc_a.replace("\\", ""))
        self.assertIn("1 of the 2 tests pass on direct-local and 1 on direct-remote", conc_b.replace("\\", ""))
        # the verdict branch of 'What held'
        self.assertIn("a direct arm failed, errored or disagreed with the other", conc_a)
        self.assertNotIn("no failure, no error and no disagreement", conc_a)
        self.assertIn("no failure, no error and no disagreement", conc_b)
        self.assertNotIn("a direct arm failed, errored or disagreed", conc_b)
        # the claims: unsupported with a disagreement, shown (in part) without
        self.assertIn("Not shown:", conc_a)
        self.assertIn("with a failure or a disagreement", conc_a)
        self.assertNotIn("Shown:", conc_a.replace("Not shown:", ""))
        self.assertIn("Shown: a 10~kb resolve under a small budget answers a stated prefix that its twins agree "
                      "with (R-budget, 1 of 2)", conc_b)
        self.assertNotIn("Not shown:", conc_b)
        # the stop sentence of the coverage gap
        self.assertIn("was not exercised successfully in this run", gap_a)
        self.assertNotIn("was answered in this run", gap_a)
        self.assertIn("was answered in this run (R-001 on direct-local (resolve stopped by time, 20~ms))", gap_b)
        self.assertNotIn("was not exercised successfully", gap_b)
        # the late 503s and the findings
        self.assertNotIn("past their budgets", gap_a)
        self.assertNotIn("only findings", gap_a)
        self.assertIn("0.4~s to 5.8~s past their budgets (R-002 5.8~s on direct-local (a finding) and 0.4~s on "
                      "direct-remote)", gap_b)
        self.assertIn("One of these 503s came more than max(1~s, the budget) after the budget and are findings, "
                      "the chapter's and the run's only findings", gap_b)
        # the footnoted tests of the family, singular
        self.assertIn("The one test of R-budget that did not pass on every arm (R-001) has its detail block", gap_a)
        self.assertIn("The one test of R-budget that did not pass on every arm (R-002) has its detail block", gap_b)
        # the family prose: the stop, the clean members, the MCP classes from the classifier
        self.assertIn("The stated stop has not been observed on this host", fam_a)
        self.assertNotIn("A stated stop was answered", fam_a)
        self.assertIn("A stated stop was answered (R-001 on direct-local (resolve stopped by time, 20~ms))", fam_b)
        self.assertIn("one of the one tests did not pass on every arm (R-001)", fam_a)
        self.assertIn("one of the two tests did not pass on every arm (R-002)", fam_b)
        self.assertIn("the one other passed on every arm exercised: R-001 (eco, 2,500~ms stop): direct-local a "
                      "pass in 20~ms; direct-remote a pass in 0.1~s; mcp a pass in 0.8~s", fam_b)
        self.assertNotIn("passed on every arm exercised", fam_a)
        for text in (fam_a, fam_b):
            self.assertIn("Via MCP: ", text)
            self.assertIn("1 answered by the host through the service", text)
            self.assertNotIn("every call was refused at submit", text)
        self.assertIn("; every test was exercised through the tools", conc_a)
        self.assertIn("1 test could not be exercised through the tools (a known limit: the time budget 1; I.4 and "
                      "appendix B)", conc_b)
        self.assertIn("Via MCP 1 pair agrees (0 fully checked)", conc_a)
        self.assertIn("Via MCP 1 pair agrees, all fully checked;", conc_b)
        self.assertIn("(2 identical pairs)", conc_b)
        # no v4 figure, nothing unresolved but the v4 scenarios the fixtures lack
        for where, text in (("a", conc_a + gap_a + fam_a), ("b", conc_b + gap_b + fam_b)):
            for lit in V4_LITERALS:
                self.assertNotIn(lit, text, "%s: the v4 literal %r survives" % (where, lit))
            for name in re.findall(r"\(UNRESOLVED CONDITION ([^)]*)\)", text):
                self.assertTrue(name.startswith("x["), "%s: unresolved %s" % (where, name))
        for ctx in (ctx_a, ctx_b):
            roots = set(n.strip("[]").lstrip("!").split("[")[0] for n in ctx.unresolved)
            self.assertEqual(roots, {"x"}, sorted(set(ctx.unresolved)))
        # the fixtures are what the test says
        self.assertIn("fam", ctx_a.chapter_values("R"))
        self.assertTrue(ctx_b.chapter_values("R")["fam"]["R-budget"]["stop_observed"])
        self.assertFalse(ctx_a.chapter_values("R")["fam"]["R-budget"]["stop_observed"])
        self.assertEqual(ctx_b.chapter_values("R")["n_finding_executions"], 1)
        self.assertTrue(ctx_b.chapter_values("R")["run_findings_only_here"])

    def test_depth_over_compared_repeats(self):
        """The owner's second review, point 2: an agreement's depth is judged over the MCP repeats compare
        compared (a repeat inconclusive by a host state left out) and stated so: the detail block's classes line,
        the chapter line's 'fully checked' count and the chapter prose name the compared repeat and the excluded
        one with its reason."""
        import sections
        d = os.path.join(self.tmp, "depth")
        shutil.copytree(self.dirs["perturbed"], d)
        comp = sc.read_json(os.path.join(d, "compare.json"))
        self.assertEqual(comp["scenarios"]["P-001"]["mcp_agreement"]["class"], "agree")
        comp["scenarios"]["P-001"]["mcp_agreement"]["host_state_repeats"] = [0]
        sc.write_json_atomic(os.path.join(d, "compare.json"), comp)
        run = reader.Run(d)
        self.assertEqual(reader.excluded_repeats(comp["scenarios"]["P-001"]), [(0, "an inconclusive host-state outcome")])
        self.assertEqual(reader.excluded_repeats(comp["scenarios"]["P-002"]), [])
        depth = run.depth["P-001"]
        self.assertEqual(depth["class"], "full")
        self.assertEqual(depth["repeats_compared"], [1])
        self.assertEqual(depth["repeats_excluded"], [{"repeat": 0, "reason": "an inconclusive host-state outcome"}])
        self.assertEqual(blocks.depth_text(depth), "fully checked on r1; r0 was an inconclusive host-state outcome")
        self.assertEqual(blocks.depth_text({"class": "full", "repeats_compared": [0, 1]}), "fully checked on r0 and r1")
        self.assertEqual(blocks.depth_text({"class": "partial", "repeats_compared": [0], "not_compared": ["a"]}),
                         "partially checked on r0 (not compared: a)")
        states = cells.scenario_states(run, "P-001")
        classes = " · ".join(blocks._classes_text(run, "P-001", states))
        self.assertIn("mcp: agree; ", classes)
        self.assertIn("fully checked on r1; r0 was an inconclusive host-state outcome", classes)
        self.assertNotIn("depth full", classes)
        units = units_mod.resolve(run.manifest, self.mappings["perturbed"])
        ctx = sections.Context(run, units, {}, {}, {}, list(units.chapters))
        for ch in ctx.chapters:
            ctx.footnotes[ch] = cells.Footnotes()
        _, mtext = sections._agreement_counts(ctx, ["P-001", "P-002", "P-003"])
        self.assertIn("fully checked 2, 1 of them on one repeat only, the other repeat an inconclusive host-state outcome)", mtext)
        _, mtext = sections._agreement_counts(ctx, ["P-003"])
        self.assertNotIn("one repeat only", mtext)
        self.assertEqual(sections.fewer_repeats_text([]), "")
        self.assertIn("with a repeat left out", sections.fewer_repeats_text([2]))
        cv = facts_mod().chapter_facts(ctx, "P")
        self.assertEqual(cv["n_mcp_full_fewer_repeats"], 1)
        self.assertEqual(cv["mcp_full_fewer_repeats_text"], "P-001 on r1, r0 an inconclusive host-state outcome")
        # the prose prints the qualifier only when a pair was checked over fewer repeats
        summary = self.R["summary"][0]
        self.assertIn("[[n_mcp_full_fewer_repeats ?", summary)
        self.assertIn("{mcp_full_fewer_repeats_text}", summary)
        # the reading guide and the glossary define the depth over the compared repeats
        overview = sc.read_json(os.path.join(PROSE_DIR, "overview.json"))
        reading = json.dumps(overview["sections"])
        self.assertIn("judged over the MCP repeats actually compared", reading)
        self.assertIn("is left out, so fully checked then holds on the other repeat only", reading)
        glossary = sc.read_json(os.path.join(PROSE_DIR, "glossary.json"))
        entry = [e for e in glossary["entries"] if e["term"] == "agreement class, depth"][0]
        self.assertIn("judged over the MCP repeats compared", entry["definition"])
        # the unperturbed run: every repeat compared, nothing excluded
        plain = reader.Run(self.dirs["plain"])
        self.assertEqual(plain.depth["P-001"]["repeats_compared"], [0, 1])
        self.assertEqual(plain.depth["P-001"]["repeats_excluded"], [])

    def test_timing_wordings(self):
        """The owner's second review, point 3: no 'upper bound', no record-versus-length cause; the owner's
        sentence in the chapter summary, 'observed under concurrent host load' in the generator's three caveats
        and the closing caveat."""
        with open(os.path.join(DOC, "sections.py"), encoding="utf-8") as fh:
            src = fh.read()
        self.assertNotIn("upper bound", src)
        self.assertEqual(src.count("observed under concurrent host load"), 3)
        for name in PROSE_FILES:
            text = json.dumps(sc.read_json(os.path.join(PROSE_DIR, name)))
            for bad in ("upper bound", "as much as the window", "record matters", "orders of magnitude",
                        "driven by", "for that reason (in this run"):
                self.assertNotIn(bad, text, "%s: %r" % (name, bad))
        summary = " ".join(self.R["summary"])
        self.assertIn("Latency varies substantially across records; these correctness controls do not isolate the "
                      "effect of query length", summary)
        self.assertIn("{x[R-pos-0027][dl][host_word]}", summary)
        self.assertNotIn("host-warm on direct-local", summary)
        self.assertIn("latency varies substantially across records", self.R["families"]["R-pos"]["prose"])
        for uid in ("R-pos-column-A", "R-pos-header-A"):
            self.assertIn("do not isolate the effect of query length", self.R["units"][uid]["note"])
        closing = sc.read_json(os.path.join(PROSE_DIR, "closing.json"))
        self.assertIn("observed under concurrent host load", closing["caveats"][0])
        # n next to the medians it counts
        self.assertIn("{wall_p50_warm_dl} at the median on direct-local (n {n_warm_dl})", summary)
        self.assertIn("(n {n_server_warm_dl} · {n_server_warm_dr} · {n_server_warm_mcp})", summary)
        # the sentences that asserted the run are conditionals on computed facts
        self.assertIn("[[mcp_all_refused_before_host ?", self.R["families"]["R-budget"]["prose"])
        self.assertIn("[[rel[REL-same_answer-0003][mcp] = pass ?", self.R["families"]["R-twin"]["prose"])
        self.assertNotIn("regularly outranked", self.R["families"]["R-discover"]["prose"])
        self.assertIn("[[mcp_all_refused_before_host ?", self.R["mcp_notes"]["units"]["R-budget"])
        self.assertNotIn("{members_text}. ", self.R["families"]["R-budget"]["prose"])
        self.assertIn("{clean_text}", self.R["families"]["R-budget"]["prose"])
        self.assertIn("Outside R-budget, {n_footnoted_elsewhere:word} test", self.R["coverage_gap"][1])

    def test_rendered_caveats_and_the_untimed_split(self):
        """The generator's caveat lines say 'observed under concurrent host load'; the host-cold MCP calls without a
        backend timing are split by what ended them from the record."""
        import sections
        run = reader.Run(self.dirs["perturbed"])
        run.load = {"n": 3, "sra_mean": 7000.0, "load_mean": 200.0, "sra_max": 1.0, "load_max": 1.0}
        units = units_mod.resolve(run.manifest, self.mappings["perturbed"])
        ctx = sections.Context(run, units, {}, {}, {}, list(units.chapters))
        self.assertIn("timings observed under concurrent host load", sections.provenance_table(ctx))
        self.assertIn("timings were observed under concurrent host load", sections.environment_table(ctx))
        self.assertIn("observed under concurrent host load", sections.closing(ctx))

        class R(object):
            order = ["a", "b", "c", "d"]
            lines = {
                "a": [{"arm": "mcp", "host_cache": "cold", "wall_ms": 5.0, "outcome": "service_by_design",
                       "error": {"code": "search_failed"}, "server": None, "mcp_phases": None,
                       "time_limited": "the host's /pattern 503 deadline relayed by the service"}],
                "b": [{"arm": "mcp", "host_cache": "cold", "wall_ms": 31000.0, "outcome": "inconclusive",
                       "error": {"code": "resolve_timeout"}, "server": None, "mcp_phases": None}],
                "c": [{"arm": "mcp", "host_cache": "cold", "wall_ms": 900.0, "outcome": "pass",
                       "server": {"elapsed_ms": 3.0}, "mcp_phases": None}],
                "d": [{"arm": "mcp", "host_cache": "warm", "wall_ms": 5.0, "outcome": "service_by_design",
                       "error": {"code": "search_failed"}, "server": None, "mcp_phases": None}]}

            def arm_lines(self, sid, arm):
                return self.lines.get(sid) or []
        self.assertEqual(sections.cold_mcp_untimed_split(R()), [
            (1, "by design", "the host's /pattern 503 deadline relayed by the service"),
            (1, "inconclusive", "the service's own deadline, `resolve_timeout`")])
        self.assertEqual(sections.cold_mcp_untimed_text(R()),
                         "1 by design (the host's /pattern 503 deadline relayed by the service) and 1 inconclusive "
                         "(the service's own deadline, `resolve_timeout`)")
        R.order = ["c", "d"]
        self.assertEqual(sections.cold_mcp_untimed_text(R()), "none")
        v = sections.closing_values(ctx)
        self.assertIn("cold_mcp_untimed_text", v)
        closing = sc.read_json(os.path.join(PROSE_DIR, "closing.json"))
        self.assertEqual(json.dumps(closing["sections"]).count("{cold_mcp_untimed_text}"), 2)
        self.assertNotIn("the calls the service's own deadline ended", json.dumps(closing))
        self.assertIn("cold_mcp_untimed_text", closing["placeholders"])

    def test_plots_at_the_owner_size(self):
        """The owner's second review, point 1: each paired scatterplot's plotting area is 62 by 55 mm (scale only
        axis), inside the half-width minipages; the phase plot keeps its two-line ticks."""
        import sections
        plot = latex.paired_plot("f.csv", "dl_ms", "mcp_ms", "x", "y", "fig:x")
        self.assertIn("scale only axis,width=62mm,height=55mm", plot)
        self.assertNotIn("0.48\\textwidth", plot)
        self.assertTrue(60 <= latex.PAIRED_PLOT_WIDTH_MM <= 70)
        run = reader.Run(self.dirs["perturbed"])
        units = units_mod.resolve(run.manifest, self.mappings["perturbed"])
        ctx = sections.Context(run, units, {}, {}, {}, list(units.chapters))
        closing = sections.closing(ctx)
        self.assertEqual(closing.count("scale only axis,width=62mm,height=55mm"), 2)
        self.assertIn("\\begin{minipage}[t]{0.49\\textwidth}", closing)
        box = latex.phase_boxplot([("submit", 1, 2, 3, 4, 5, 10, 0), ("paging", 1, 2, 3, 4, 5, 8, 2)], "fig:p")
        self.assertIn("xticklabel style={align=center", box)
        self.assertIn("{paging\\\\(n 8 of 10)}", box)

    def test_units_joined_and_variant_labels(self):
        """Low items 8 and 12: a number and its unit never break apart in the prose; a variant label prints its
        thousands separators and a scenario's label never nests parentheses."""
        self.assertEqual(latex.join_units("25 s past the 20,000 ms budget, 2,000 bp; 7762 \\% and 5 samples"),
                         "25~s past the 20,000~ms budget, 2,000~bp; 7762~\\% and 5 samples")
        self.assertEqual(latex.join_units("r0 s; 3 seconds; in 2 steps"), "r0 s; 3 seconds; in 2 steps")
        self.assertEqual(latex.inline("`max_query_bp` 2,000 bp in 0.6 s"), "\\texttt{max\\_query\\_bp} 2,000~bp in 0.6~s")
        self.assertEqual(latex.inline("*in 4 ms*"), "\\emph{in 4~ms}")
        for raw, printed in (("2000 bp (pfa)", "2,000 bp (pfa)"), ("eco, from base 5000", "eco, from base 5,000"),
                             ("eco, 3000 ms stop", "eco, 3,000 ms stop"), ("R-pos-0027", "R-pos-0027"),
                             ("10,001 bp", "10,001 bp"), ("1492R", "1492R"), ("hiv1155.e20", "hiv1155.e20"),
                             ("A^100000 on a self-loop", "A^100000 on a self-loop"), ("L-paths-0001", "L-paths-0001"),
                             ("cal.XP_720315.1 (t12)", "cal.XP_720315.1 (t12)"), ("max_work_units 300000", "max_work_units 300,000")):
            self.assertEqual(units_mod.variant_text(raw), printed, raw)
        f = facts_mod()
        self.assertEqual(f._label({"id": "R-pos-0027", "variant": "2,000 bp (pfa)"}), "R-pos-0027 (2,000 bp, pfa)")
        self.assertEqual(f._label({"id": "R-budget-0001", "variant": "eco, 3,000 ms stop"}), "R-budget-0001 (eco, 3,000 ms stop)")
        self.assertEqual(f._label({"id": "R-edge-0001", "variant": ""}), "R-edge-0001")
        # the table headers and the facts print the same label
        import sections
        run = reader.Run(self.dirs["perturbed"])
        units = units_mod.resolve(run.manifest, self.mappings["perturbed"])
        ctx = sections.Context(run, units, {}, {}, {}, list(units.chapters))
        for ch in ctx.chapters:
            ctx.footnotes[ch] = cells.Footnotes()
        self.assertIn("eco, 2,500 ms stop", sections.unit_section(ctx, "R", "R-budget"))
        self.assertEqual(f.scenario_facts(ctx, "R-001")["variant"], "eco, 2,500 ms stop")
        self.assertEqual(sections._variant({"variant": "2000 bp (pfa)"}), "2,000 bp (pfa)")
        self.assertEqual(sections._variant({}, "R-001"), "R-001")

    def test_unit_facts_and_record_key_wording(self):
        """A unit's note under its table reads the unit's own facts (the MCP classes over its members); I.2 says
        where the record key is; the README no longer places it in I.3."""
        import sections
        run = reader.Run(self.dirs["perturbed"])
        mapping = json.loads(json.dumps(self.mappings["perturbed"]))
        mapping["_units"]["R-budget"]["members"] = ["R-001"]
        mapping["_units"]["R-two"] = {"family": "R-budget", "unit_type": "pairs", "unit_title": "two",
                                      "property": "p", "members": ["R-002"]}
        mapping["_families"]["R-budget"]["units"] = ["R-budget", "R-two"]
        mapping["R-002"]["unit_id"] = "R-two"
        units = units_mod.resolve(run.manifest, mapping)
        f = facts_mod()
        prose = {"R": {"mcp_notes": {"units": {"R-two": "{n_tests} test, {mcp_classes_text}"
                                                        "[[mcp_all_refused_before_host ? (refused) | (not refused)]]"}}}}
        ctx = sections.Context(run, units, prose, {}, {}, list(units.chapters))
        for ch in ctx.chapters:
            ctx.footnotes[ch] = cells.Footnotes()
        uv = f.unit_facts(ctx, "R-two")
        self.assertEqual(uv["n_tests"], 1)
        self.assertTrue(uv["mcp_all_not_exercised"])
        self.assertFalse(uv["mcp_all_refused_before_host"])
        self.assertIn("R-002", uv["x"])
        note = sections.mcp_unit_note(ctx, "R", "R-two", ["R-002"])
        self.assertIn("1 test, 1 ", note)
        self.assertIn("(not refused)", note)
        self.assertEqual(sections.mcp_unit_note(ctx, "R", "R-budget", ["R-001"]), "")
        fv = f.family_facts(ctx, "R-budget")
        for key in ("mcp_classes_text", "mcp_all_refused_before_host", "mcp_all_not_exercised", "clean_ids",
                    "clean_text", "mcp_refused_before_host_codes"):
            self.assertIn(key, fv)
        self.assertEqual(fv["clean_ids"], "R-001")
        self.assertEqual(f.arm_facts(ctx, "R-001", "direct-local")["host_word"], "host-warm")
        self.assertEqual(f.HOST_WORD[cells.WARM_AFTER_INTERRUPTED], "warm after an interrupted execution")
        overview = sc.read_json(os.path.join(PROSE_DIR, "overview.json"))
        self.assertIn("keyed on each chapter page and in appendix D", json.dumps(overview["sections"]))
        self.assertNotIn("keyed in the glossary", json.dumps(overview))
        with open(os.path.join(SUITE, "README.md"), encoding="utf-8") as fh:
            readme = fh.read()
        self.assertNotIn("record key in I.3", readme)
        self.assertIn("test_two_fixtures_change_conclusions", readme)


def facts_mod():
    import facts
    return facts


if __name__ == "__main__":
    unittest.main()
