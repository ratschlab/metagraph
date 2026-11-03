"""Regression cases for text extraction and baseline updates."""
import contextlib
import importlib.util
import io
import json
from pathlib import Path
import sys
import unittest
from unittest import mock

SPEC = importlib.util.spec_from_file_location("check_text_rules", Path(__file__).parents[1] / "check_text_rules.py")
RULES = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RULES)


class TextExtraction(unittest.TestCase):
    def test_semicolon_spacing(self):
        for text in ("This is wrong;fix the input.", "This is wrong; fix the input.", "Stop;retry"):
            with self.subTest(text=text):
                self.assertIn("R1", [r for r, _ in RULES.check_text(text, {}, {})])

    def test_short_errors_and_terms(self):
        self.assertIn("R4", [r for r, _ in RULES.check_text("Oops", {}, {})])
        cfg = {"terms": [{"use": "k-mer", "not": ["kmer"]}]}
        self.assertIn("R8", [r for r, _ in RULES.check_text("kmer", cfg, {})])

    def test_jargon_needs_prose(self):
        cfg, scope = {"jargon": ["re:\\b[a-z]+_[a-z_]+\\b"]}, {"jargon": True}
        self.assertEqual(RULES.check_text("min_kmers", cfg, scope), [])
        self.assertIn("R5", [r for r, _ in RULES.check_text("Set min_kmers to 3 in the search form.", cfg, scope)])

    def test_clean_text_and_code(self):
        for text in ("Read the file.", "SELECT * FROM records;"):
            self.assertEqual(RULES.check_text(text, {}, {}), [])

    def test_inline_admonitions_keep_their_line(self):
        for directive in ("note", "warning", "attention"):
            units = RULES.markdown_units(f"\n.. {directive}:: Stop; retry.\n", rst=True)
            self.assertEqual(units, [(2, "Stop; retry.")])

    def test_code_directives_stay_excluded(self):
        self.assertEqual(RULES.markdown_units(".. code-block:: python\n\n   print('Oops')\n", rst=True), [])

    def test_admonition_options_are_metadata(self):
        src = ".. warning:: Stop; retry.\n   :class: kmer\n\n   Read the file.\n"
        self.assertEqual(RULES.markdown_units(src, rst=True), [(1, "Stop; retry."), (4, "Read the file.")])

    def test_include_paths_are_not_user_text(self):
        src = '#include "kmer/kmer_extractor.hpp"\nlogger->error("Oops");\n'
        units = RULES.cpp_units(src, {"cpp_all_strings": ["*.cpp"]}, "config.cpp")
        self.assertEqual(units, [(2, "Oops")])

    def test_cli_headers_are_in_scope(self):
        cfg = json.loads((Path(__file__).parents[1] / "text-rules.json").read_text())
        with mock.patch.object(RULES, "repo_files", return_value=["metagraph/src/cli/parse_sequences.hpp"]), \
                mock.patch.object(Path, "read_text", return_value='logger->warn("Oops");\n'):
            results = RULES.scan(Path("/unused"), cfg)
        self.assertEqual(results["metagraph/src/cli/parse_sequences.hpp"][0][1], "R4")


class AllowComment(unittest.TestCase):
    def scan(self, src):
        cfg = {"scopes": [{"name": "docs", "include": ["a.md"]}], "terms": [{"use": "k-mer", "not": ["kmer", "kmers"]}]}
        with mock.patch.object(RULES, "repo_files", return_value=["a.md"]), \
                mock.patch.object(Path, "read_text", return_value=src):
            return [rule for _, rule, _ in RULES.scan(Path("/unused"), cfg).get("a.md", [])]

    def test_allow_names_rule_and_term(self):
        self.assertEqual(self.scan("<!-- text-rules-allow: R8 kmer the flag name -->\nCount each kmer twice.\n"), [])

    def test_allow_without_matching_term_is_ignored(self):
        self.assertEqual(self.scan("<!-- text-rules-allow: R8 genome -->\nCount each kmer twice.\n"), ["R8"])

    def test_allow_requires_reason_before_comment_end(self):
        for comment in ("text-rules-allow: R8 kmer", "<!-- text-rules-allow: R8 kmer -->",
                        "<!-- text-rules-allow: R8 kmer   -->", "/* text-rules-allow: R8 kmer */"):
            with self.subTest(comment=comment):
                self.assertEqual(self.scan(comment + "\nCount each kmer twice.\n"), ["R8"])

    def test_allow_does_not_match_plural(self):
        self.assertEqual(self.scan("<!-- text-rules-allow: R8 kmer the flag name -->\nCount the kmers twice.\n"), ["R8"])

    def test_next_line_cannot_supply_reason(self):
        self.assertEqual(RULES.ALLOW.findall("// text-rules-allow: R8 kmer\nCount each kmer twice."), [])

    def test_allow_matches_punctuation(self):
        self.assertEqual(self.scan("<!-- text-rules-allow: R1 ; published punctuation -->\nStop; retry.\n"), [])


class BaselineCheck(unittest.TestCase):
    def test_drop_passes_with_a_note(self):
        out = io.StringIO()
        with mock.patch.object(sys, "argv", ["checker", "--root", "/unused", "--config", "config.json",
                                            "--baseline", "baseline.json"]), \
                mock.patch.object(RULES, "scan", return_value={}), \
                mock.patch.object(Path, "read_text", autospec=True,
                                  side_effect=lambda p, *a, **k: json.dumps({"a.md": {"R1": 1}} if p.name == "baseline.json" else {})), \
                mock.patch.object(Path, "exists", return_value=True), contextlib.redirect_stdout(out):
            self.assertEqual(RULES.main(), 0)
        self.assertIn("Run --update-baseline to lower it", out.getvalue())


class BaselineUpdate(unittest.TestCase):
    def update(self, base, results, exists=True):
        def read(path, *args, **kwargs):
            return json.dumps(base if path.name == "baseline.json" else {})
        with mock.patch.object(sys, "argv", ["checker", "--root", "/unused", "--config", "config.json",
                                            "--baseline", "baseline.json", "--update-baseline"]), \
                mock.patch.object(RULES, "scan", return_value=results), \
                mock.patch.object(Path, "read_text", autospec=True, side_effect=read), \
                mock.patch.object(Path, "exists", return_value=exists), \
                mock.patch.object(Path, "write_text") as write, contextlib.redirect_stdout(io.StringIO()):
            rc = RULES.main()
            return rc, write.call_args_list

    def test_increase_refused_without_writing(self):
        result = {"README.md": [(1, "R1", "A;B"), (2, "R1", "C;D")]}
        self.assertEqual(self.update({"README.md": {"R1": 1}}, result), (1, []))

    def test_new_file_refused_without_writing(self):
        self.assertEqual(self.update({}, {"new.md": [(1, "R4", "Oops")]}), (1, []))

    def test_missing_baseline_refused_for_existing_hits(self):
        self.assertEqual(self.update({}, {"README.md": [(1, "R1", "A;B")]}, exists=False), (1, []))

    def test_reduction_is_written(self):
        rc, writes = self.update({"README.md": {"R1": 2}}, {"README.md": [(1, "R1", "A;B")]})
        self.assertEqual(rc, 0)
        self.assertEqual(json.loads(writes[0].args[0]), {"README.md": {"R1": 1}})

    def test_complete_cleanup_is_written(self):
        rc, writes = self.update({"README.md": {"R1": 1}}, {})
        self.assertEqual(rc, 0)
        self.assertEqual(json.loads(writes[0].args[0]), {})


if __name__ == "__main__":
    unittest.main()
