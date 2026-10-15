"""The report's colour marks (the owner, 2026-10-15): every printed execution time above 10 s yellow, above 55 s red,
in report.md (a span with the class and the tint inline) and on the md -> pdf path (report_colour.lua turns the class
into a \\colorbox), with one legend sentence; differences and connection set-ups stay plain."""

import contextlib
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.dirname(HERE)
sys.path.insert(0, SUITE)
sys.path.insert(0, HERE)

import compare as cmp  # noqa: E402
import merge_run  # noqa: E402
import report  # noqa: E402
import suite_common as sc  # noqa: E402
from test_compare_report import CASES, build_run  # noqa: E402

YELLOW = '<span class="t-yellow" style="background-color:#ffe866">%s</span>'
RED = '<span class="t-red" style="background-color:#ff8a80">%s</span>'


class FormatTest(unittest.TestCase):
    def test_thresholds_are_named_constants(self):
        self.assertEqual((report.TIME_MARK_YELLOW_MS, report.TIME_MARK_RED_MS), (10000.0, 55000.0))
        self.assertEqual(report.time_mark(None), None)
        self.assertEqual(report.time_mark(10000.0), None)          # above, not at
        self.assertEqual(report.time_mark(10000.1), "yellow")
        self.assertEqual(report.time_mark(55000.0), "yellow")
        self.assertEqual(report.time_mark(55000.1), "red")
        self.assertEqual(report.time_mark(60187.0), "red")
        self.assertEqual(report.time_mark(-5.0), None)

    def test_fmt_time(self):
        self.assertEqual(report.fmt_time(None), "-")
        self.assertEqual(report.fmt_time(9.96), "10.0")
        self.assertEqual(report.fmt_time(9999.0), "9999")
        self.assertEqual(report.fmt_time(12345.6), YELLOW % "12346")
        self.assertEqual(report.fmt_time(61000.0), RED % "61000")
        self.assertEqual(report.fmt_ms(61000.0), "61000")           # a difference stays plain
        self.assertEqual(report.visible(report.fmt_time(61000.0)), "61000")
        # the mark follows the value as printed: fmt_ms rounds to whole milliseconds from 100 ms, so 10,000.4 ms
        # prints 10000, which is not above 10 s (and 55,000.4 prints 55000, yellow, not red)
        self.assertEqual(report.fmt_time(10000.4), "10000")
        self.assertEqual(report.fmt_time(10000.6), YELLOW % "10001")
        self.assertEqual(report.fmt_time(55000.4), YELLOW % "55000")
        self.assertEqual(report.fmt_time(55000.6), RED % "55001")

    def test_legend(self):
        legend = report.time_legend()
        self.assertIn("above 10 s are marked %s and above 55 s %s" % (YELLOW % "yellow", RED % "red"), legend)
        self.assertIn("differences between arms and connection set-ups are not marked", legend)

    def test_table_reads_the_visible_text(self):
        t = report.table(["area", "wall"], [["P", report.fmt_time(61000.0)], ["R", report.fmt_time(5.0)]])
        lines = t.strip().split("\n")
        # right-aligned (numeric by its visible text) and the separator sized by the visible width, not the span
        self.assertEqual(lines[1], "|:----|-----:|")
        self.assertIn("| P | %s |" % (RED % "61000"), lines[2])


class ReportTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.d = self.tmp.name
        build_run(self.d)
        # P-001's warm direct-local execution takes 61 s, M-002's cold one 12 s, R-001's failing one 56 s
        path = os.path.join(self.d, "arms", "direct-local", "executions.jsonl")
        rows, _ = sc.read_jsonl(path)
        walls = {("P-001", 1): 61000.0, ("M-002", 0): 12000.0, ("R-001", 0): 56000.0}
        for r in rows:
            w = walls.get((r["scenario"], r["repeat"]))
            if w is not None:
                r["wall_ms"] = w
                r["server"]["elapsed_ms"] = w / 2.0
                r["ended_utc"] = sc.iso_utc(sc.parse_utc(r["started_utc"]) + w / 1000.0)
        os.remove(path)
        for r in rows:
            sc.append_jsonl(path, r)
        merge_run.merge(self.d)
        sc.write_json_atomic(os.path.join(self.d, "compare.json"), cmp.compare(self.d, (CASES, {"cases": 5})))
        with contextlib.redirect_stdout(io.StringIO()):
            self.text = report.build_report(self.d)

    def tearDown(self):
        self.tmp.cleanup()

    def test_marks_in_the_tables_and_the_lines(self):
        legend = report.time_legend()
        self.assertEqual(self.text.count(legend), 1)
        self.assertLess(self.text.index(legend), self.text.index("# "))      # right after the front matter
        timing = self.text.split("# Timing", 1)[1].split("\n# ", 1)[0]
        # P's warm direct-local execution: 61 s red, its server half (30.5 s) yellow; M cold: the p50 of 6 s plain,
        # the p90 and max above 10 s yellow, the server times plain; R's failing execution: 56 s red, 28 s yellow
        self.assertIn("| P | direct-local | warm | 1 | %s | %s | %s | 1 | %s | %s | %s |" % (
            (RED % "61000",) * 3 + (YELLOW % "30500",) * 3), timing)
        self.assertIn("| M | direct-local | cold | 2 | 6002 | %s | %s | 2 | 3001 | 5400 | 6000 |" % (
            YELLOW % "10800", YELLOW % "12000"), timing)
        self.assertIn("| R | direct-local | warm | 1 | %s | %s | %s | 1 | %s | %s | %s |" % (
            (RED % "56000",) * 3 + (YELLOW % "28000",) * 3), timing)
        # the failing R-001 line names its wall, marked
        failures = self.text.split("# Failures", 1)[1].split("\n# ", 1)[0]
        self.assertIn("wall %s ms" % (RED % "56000"), failures)
        # differences stay plain: the cold vs warm delta of P on direct-local is 61000 - 9 = 60991 (no span)
        self.assertIn("| direct-local | 1 | 60991 | 60991 |", timing)
        self.assertNotIn(RED % "60991", timing)
        # short times are plain everywhere
        self.assertIn("| P | direct-remote | warm | 4 | 104 | 108 | 110 | 4 | 52.2 | 54.2 | 55.0 |", timing)

    @unittest.skipUnless(shutil.which("pandoc"), "no pandoc")
    def test_latex_output_colours_the_marks(self):
        md_path = os.path.join(self.d, "report.md")
        with open(md_path, "w", encoding="ascii", errors="replace") as fh:
            fh.write(self.text)
        proc = subprocess.run([shutil.which("pandoc"), md_path, "-t", "latex", "--standalone", "--lua-filter",
                               report.LUA_FILTER], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              universal_newlines=True)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("\\definecolor{suitered}{HTML}{FF8A80}", proc.stdout)
        self.assertIn("\\definecolor{suiteyellow}{HTML}{FFE866}", proc.stdout)
        self.assertIn("{\\setlength{\\fboxsep}{1pt}\\colorbox{suitered}{61000}}", proc.stdout)
        self.assertIn("\\colorbox{suiteyellow}{30500}", proc.stdout)
        self.assertIn("\\colorbox{suiteyellow}{yellow}", proc.stdout)      # the legend's words carry their tint
        self.assertNotIn("<span", proc.stdout)
        self.assertNotIn("t-red", proc.stdout)
        self.assertIn("\\usepackage{xcolor}", proc.stdout)                # pandoc's template loads it

    @unittest.skipUnless(shutil.which("pandoc") and shutil.which("xelatex"), "no pandoc + xelatex")
    def test_pdf_is_produced(self):
        md_path = os.path.join(self.d, "report.md")
        with open(md_path, "w", encoding="ascii", errors="replace") as fh:
            fh.write(self.text)
        pdf = os.path.join(self.d, "report.pdf")
        ok, msg = report.to_pdf(md_path, pdf)
        self.assertTrue(ok, msg)
        with open(pdf, "rb") as fh:
            head = fh.read(8)
        self.assertEqual(head[:5], b"%PDF-")
        self.assertGreater(os.path.getsize(pdf), 10000)


if __name__ == "__main__":
    unittest.main()
