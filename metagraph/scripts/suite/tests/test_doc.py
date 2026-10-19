"""Offline tests of the results document generator (doc/make_doc.py; DESIGN-DOC §11): the prose data complete for
every unit and family of the v4 and the v5 mapping, the mappings complete over the frozen manifests, every check kind
described, every chapter's prose free of placeholders, run figures and banned wordings and resolvable on both
mappings, the generator on the fake run of test_compare_report (the .tex, the index, the verdict, the escaping; the
PDF when xelatex is installed), the conclusions of chapters R, M and T on fixtures whose data differ, the cell
grammar and the expectation formatters on chapter R's real expectations."""

import json
import os
import inspect
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
import facts  # noqa: E402
import latex  # noqa: E402
import make_doc  # noqa: E402
import reader  # noqa: E402
import sections  # noqa: E402
import units as units_mod  # noqa: E402
import vocabulary  # noqa: E402

V4_MANIFEST = os.path.join(SUITE, "manifest", "suite_v4.json")
V4_MAPPING = os.path.join(DOC, "units", "suite_v4.json")
V5_MANIFEST = os.path.join(SUITE, "manifest", "suite_v5.json")
V5_MAPPING = os.path.join(DOC, "units", "suite_v5.json")
PROSE_DIR = os.path.join(DOC, "prose")
# the chapters whose prose is written: their files, the overview, the glossary and their check kinds' catalogue
# entries may hold no placeholder ("to be written"); every chapter since the full document (2026-10-15)
WRITTEN_CHAPTERS = ("R", "T", "P", "L", "X", "M", "Q", "C", "N", "D")
PLACEHOLDER = "to be written"
# Keys whose prose is still pending may carry a placeholder (the generator ships one); empty since the prose of
# the owner's review (round 2, 2026-10-15) landed, and again since the v5 prose round wrote R-warm (the warm-ups,
# manifest v5 only). Name a file or "<chapter>.json: <key>" here while prose is awaited.
ROUND2_PENDING = ()


# The causal or attributive timing wordings no prose file and no line of the generator may carry (the owner's
# second review, point 3; VERIFY-R4, issue 1): each was removed once and is banned since.
BANNED_TIMING_WORDINGS = ("upper bound", "as much as the window", "record matters", "orders of magnitude",
                          "driven by", "for that reason (in this run", "point to", "points to", "reflects the",
                          "reflect the", "rather than a slower host", "rather than to the cache",
                          "rather than the front", "at this load", "is read as one", "should agree",
                          "a first client sees", "host-bound part", "separates little",
                          "ratio of seconds to milliseconds", "lies in the service's own phases rather than",
                          "warmer)")


def v4():
    manifest = sc.read_json(V4_MANIFEST)
    mapping = units_mod.load_mapping(V4_MAPPING)
    return manifest, mapping, units_mod.resolve(manifest, mapping)


def v5():
    manifest = sc.read_json(V5_MANIFEST)
    mapping = units_mod.load_mapping(V5_MAPPING)
    return manifest, mapping, units_mod.resolve(manifest, mapping)


V6_MANIFEST = os.path.join(SUITE, "manifest", "suite_v6.json")
V6_MAPPING = os.path.join(DOC, "units", "suite_v6.json")
# the frozen v6 manifest's identity (gr/suite-v6 at a71d9177f, merged with this generator: v5 with every
# default-budget /pattern request at 110,000 ms and the 25 /pattern warm-ups P-warm-0001..0025, 725 scenarios; the
# gate's scenarios, their order and their required executions v5's)
V6_CONTENT_SHA = "49907dbb95ca0a373199d40cabe888b4c136f376295a1fc526a742eae8231cd1"
V6_PATTERN_BUDGET_MS = 110000
V6_WARM_UPS = ["P-warm-%04d" % i for i in range(1, 26)]


def v6_manifest():
    """The frozen v6 manifest, manifest/suite_v6.json (in the tree since the merge with gr/suite-v6; its content sha
    is the frozen one: test_mapping_complete_v6 asserts it)."""
    return sc.read_json(V6_MANIFEST)


def v6():
    manifest = v6_manifest()
    mapping = units_mod.load_mapping(V6_MAPPING)
    return manifest, mapping, units_mod.resolve(manifest, mapping)


def _prose_texts(items):
    """The prose strings of a chapter slot (`intro`, `summary`, `coverage_gap`, `purpose`): its paragraphs and the
    caption, columns and cells of a {"table": ...} item (sections.prose_blocks; chapter C's map of the documents)."""
    out = []
    for item in ([items] if isinstance(items, str) else (items or [])):
        if isinstance(item, str):
            out.append(item)
        elif isinstance(item, dict) and isinstance(item.get("table"), dict):
            t = item["table"]
            out.extend([t.get("caption") or ""] + list(t.get("columns") or [])
                       + [c for row in (t.get("rows") or []) for c in row])
    return out


# ---------------------------------------------------------------------------------------------------------------
# the maintenance tests on the real v4 data

class MappingTest(unittest.TestCase):
    def check_complete(self, manifest, mapping, units, version):
        """Every scenario in exactly one unit, units non-empty, types in the vocabulary, series of at most 8
        distinct variants, groups consistent, every relation under exactly one unit holding a or b, the mapping
        keyed to the frozen manifest's content sha."""
        self.assertEqual(mapping["_manifest"]["content_sha256"], manifest["content_sha256"],
                         "the mapping was built for another v%s manifest" % version)
        self.assertEqual(mapping["_manifest"]["version"], manifest["manifest_version"])
        self.assertEqual(units_mod.validate(units, manifest), [])
        self.assertEqual(len(units.by_scenario), len(manifest["scenarios"]))
        self.assertEqual(units.chapters, ["R", "T", "P", "L", "X", "M", "Q", "C", "N", "D"])
        self.assertEqual(units.unit_of("R-pos-0027"), "R-pos-column-B")
        self.assertEqual(units.relation_unit["REL-prefix_equals-0001"], "R-budget")
        self.assertEqual(units_mod.default_mapping_path(manifest), os.path.join(DOC, "units", "suite_v%s.json" % version))

    def test_mapping_complete(self):
        """The v4 mapping over the frozen v4 manifest (681 scenarios, 107 units)."""
        manifest, mapping, units = v4()
        self.check_complete(manifest, mapping, units, 4)
        self.assertEqual(len(units.chapter_scenarios("R")), 81)
        self.assertEqual(len(units.units), 107)
        self.assertEqual(len([u for f in units.chapter_families("R") for u in units.families[f]["units"]]), 14)

    def test_mapping_complete_v5(self):
        """The v5 mapping over the frozen v5 manifest (700 scenarios): the v4 units kept, chapter R's 14 units
        plus the warm-ups' unit, the long walks paired with the gate's walks, every v5-only scenario placed once,
        every v4 id in the unit it had (the four gate walks moved beside their long twins apart)."""
        manifest, mapping, units = v5()
        self.check_complete(manifest, mapping, units, 5)
        self.assertEqual(len(units.chapter_scenarios("R")), 90)
        self.assertEqual(len(units.units), 109)
        r_units = [u for f in units.chapter_families("R") for u in units.families[f]["units"]]
        self.assertEqual(len(r_units), 15)
        self.assertIn("R-warm", r_units)
        self.assertEqual(sorted(set(r_units) - set(u for f in v4()[2].chapter_families("R")
                                                  for u in v4()[2].families[f]["units"])), ["R-warm"])
        self.assertEqual(units.units["R-warm"]["members"], ["R-warm-%04d" % i for i in (7, 1, 2, 8, 3, 4, 5, 6)])
        self.assertEqual(units.units["T-longwalk"]["members"],
                         ["T-walk-0006", "T-longwalk-0001", "T-walk-0012", "T-longwalk-0002", "T-walk-0018",
                          "T-longwalk-0003", "T-walk-0036", "T-longwalk-0004"])
        self.assertEqual(units.units["T-longwalk"]["unit_type"], "pairs")
        self.assertEqual(units.families["T-longwalk"]["chapter"], "T")
        for sid, uid in (("R-edge-0025", "R-edge-accepted"), ("T-attempts-0003", "T-attempts"),
                         ("T-caps-0002", "T-caps"), ("T-caps-0004", "T-caps"), ("T-quorum-0002", "T-quorum"),
                         ("T-seeds-0002", "T-seeds")):
            self.assertEqual(units.unit_of(sid), uid, sid)
        v4_manifest, v4_mapping, v4_units = v4()
        v4_ids = set(s["id"] for s in v4_manifest["scenarios"])
        moved = set(sid for sid in v4_ids if v4_units.unit_of(sid) != units.unit_of(sid))
        self.assertEqual(moved, {"T-walk-0006", "T-walk-0012", "T-walk-0018", "T-walk-0036"})
        self.assertEqual(sorted(set(units.by_scenario) - v4_ids),
                         ["R-edge-0025"] + ["R-warm-%04d" % i for i in range(1, 9)] + ["T-attempts-0003"]
                         + ["T-caps-%04d" % i for i in (2, 3, 4)] + ["T-longwalk-%04d" % i for i in range(1, 5)]
                         + ["T-quorum-0002", "T-seeds-0002"])
        # a v4 id whose title did not change keeps its v4 label; a changed title gets a label from the new one
        for sid in ("R-pos-0001", "R-edge-0004", "T-walk-0001"):
            self.assertEqual(units.by_scenario[sid]["variant"], v4_units.by_scenario[sid]["variant"], sid)
        self.assertEqual(units.by_scenario["R-budget-0001"]["variant"], "eco, 480 ms stop")
        self.assertEqual(units.by_scenario["R-twin-0002"]["variant"], "eco, complete, 60 s (warm-up)")
        self.assertEqual(units.by_scenario["T-walk-0004"]["group"], "right, 300 bp")

    def test_mapping_complete_v6(self):
        """The v6 mapping over the frozen v6 manifest (725 scenarios: v5's 700, every default-budget /pattern request
        at 110,000 ms, and the 25 /pattern warm-ups P-warm-0001..0025, each right before its target; the free /pattern
        scenarios build_manifest.keep_gate_warm_repeats moves change places): every v5 scenario keeps its v5 entry,
        every v5 unit and family its entry (a unit's members in v6's order); the warm-ups form one unit, P-warm (family P-warm, chapter P, replicates,
        variant 'of <target>'), as v5's walk warm-ups form R-warm; the generator picks the mapping by the sha; no
        label of a pattern chapter names a /pattern budget (60 s, 110,000 or 120,000 ms; the budgets those chapters
        name are P-caps-0004's 300 ms; '60 bp' and 'A^60' are lengths)."""
        manifest, mapping, units = v6()
        self.assertEqual(manifest["content_sha256"], V6_CONTENT_SHA)
        self.assertEqual(sc.scenarios_sha256(manifest["scenarios"]), V6_CONTENT_SHA)
        self.assertEqual(len(manifest["scenarios"]), 725)
        self.check_complete(manifest, mapping, units, 6)
        self.assertEqual(units_mod.default_mapping_path({"manifest_version": 6, "content_sha256": V6_CONTENT_SHA}),
                         V6_MAPPING)
        v5_manifest, v5_mapping, v5_units = v5()
        ids5 = [s["id"] for s in v5_manifest["scenarios"]]
        ids6 = [s["id"] for s in manifest["scenarios"]]
        warm = [s for s in manifest["scenarios"] if s["id"] not in set(ids5)]
        self.assertEqual(sorted(s["id"] for s in warm), V6_WARM_UPS)
        titles5 = dict((s["id"], s.get("title")) for s in v5_manifest["scenarios"])
        self.assertEqual([s.get("title") for s in manifest["scenarios"] if s["id"] in titles5],
                         [titles5[x] for x in ids6 if x in titles5])
        self.assertEqual([r.get("id") for r in manifest["relations"]], [r.get("id") for r in v5_manifest["relations"]])
        # every v5 entry is v5's; the warm-ups' entries, unit and family are the only new ones
        for key in v5_mapping:
            if key in ("_about", "_manifest", "_units", "_families"):
                continue
            self.assertEqual(mapping[key], v5_mapping[key], key)
        for part in ("_units", "_families"):
            for key in v5_mapping[part]:
                got, want = dict(mapping[part][key]), dict(v5_mapping[part][key])
                if part == "_units" and isinstance(want.get("members"), list):
                    # a unit's members are v5's, in v6's order (keep_gate_warm_repeats moved some)
                    self.assertEqual(got.pop("members"), sorted(want.pop("members"), key=ids6.index), key)
                self.assertEqual(got, want, key)
            self.assertEqual(sorted(set(mapping[part]) - set(v5_mapping[part])), ["P-warm"])
        self.assertEqual(sorted(k for k in mapping if not k.startswith("_") and k not in v5_mapping), V6_WARM_UPS)
        self.assertEqual([k for k in mapping if not k.startswith("_")], ids6)
        self.assertEqual(mapping["_manifest"], {"version": 6, "content_sha256": V6_CONTENT_SHA, "scenarios": 725})
        self.assertIn("manifest v6", mapping["_about"])
        self.assertIn(V6_CONTENT_SHA, mapping["_about"])
        self.assertNotIn(v5_manifest["content_sha256"], mapping["_about"])
        for sid in ids5:
            self.assertEqual(units.by_scenario[sid], v5_units.by_scenario[sid], sid)
        for uid, u in v5_units.units.items():
            self.assertEqual(dict(units.units[uid], members=None), dict(u, members=None), uid)
            self.assertEqual(sorted(units.units[uid]["members"]), sorted(u["members"]), uid)
        for fid, f in v5_units.families.items():
            self.assertEqual(units.families[fid], f, fid)
        self.assertEqual(len(units.units), 110)
        # the warm-ups' unit: in chapter P, in the manifest's order, each named by its target
        pw = units.units["P-warm"]
        self.assertEqual((pw["family"], pw["unit_type"]), ("P-warm", "replicates"))
        self.assertEqual(units.families["P-warm"]["chapter"], "P")
        self.assertEqual(units.families["P-warm"]["units"], ["P-warm"])
        self.assertEqual(pw["members"], [s["id"] for s in warm])
        for s in warm:
            self.assertEqual(units.unit_of(s["id"]), "P-warm")
            self.assertEqual(units.by_scenario[s["id"]]["chapter"], "P")
            self.assertEqual(units.by_scenario[s["id"]]["variant"], "of %s" % s["warms"])
            self.assertEqual(ids6[ids6.index(s["id"]) + 1], s["warms"])
            self.assertIn("warm-up", s.get("tags") or [])
        for s in warm:
            self.assertIn(s["warms"].rsplit("-", 1)[0], pw["property"])  # every target's family is named
        budget = re.compile(r"\b60 s\b|\b60,000 ms\b|\b60000\b|\b110,000\b|\b110000\b|\b120,000\b|\b120000\b")
        for uid, u in units.units.items():
            if u["family"][0] in "PLXMQ":
                self.assertIsNone(budget.search(json.dumps(u, ensure_ascii=False)), uid)
        for sid, m in units.by_scenario.items():
            if m["chapter"] in "PLXMQ":
                self.assertIsNone(budget.search("%s %s" % (m.get("variant"), m.get("group"))), sid)
        self.assertEqual(units.by_scenario["P-caps-0004"]["variant"], "16 patterns under a 300 ms budget")
        # the pattern scenarios of v6 state the 110,000 ms budget: v5's 449 default-budget requests and the warm-ups
        budgets = [s["request"].get("time_budget_ms") for s in manifest["scenarios"]
                   if s.get("route") == "POST /pattern" and isinstance(s.get("request"), dict)]
        self.assertEqual(budgets.count(V6_PATTERN_BUDGET_MS), 449 + len(V6_WARM_UPS))
        self.assertNotIn(60000, budgets)
        self.assertNotIn(120000, budgets)

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
        title and prose, every unit a property line); the message names the missing keys. Over the v4, the v5 and
        the v6 mapping (the v5-only family and units have at least the mapping's draft; v6's P-warm its prose, the
        v6 delta of 2026-10-17)."""
        for version in (v4, v5, v6):
            self.check_prose_complete(version())

    def check_prose_complete(self, loaded):
        manifest, mapping, units = loaded
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
        v5_manifest = v5()[0]
        problems = []
        for ch in WRITTEN_CHAPTERS:
            data = sc.read_json(os.path.join(PROSE_DIR, "%s.json" % ch)) or {}
            # a prose slot's list may hold a {"table": ...} item beside its paragraphs (sections.prose_blocks; chapter C's
            # map of the capabilities documents): its cells are prose too
            texts = [("intro", " ".join(_prose_texts(data.get("intro"))))]
            for key in ("summary", "coverage_gap"):
                v = data.get(key)
                texts.append((key, " ".join(_prose_texts(v)) if isinstance(v, list) else (v or "")))
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
            kinds = set(e.get("check") for m in (manifest, v5_manifest) for s in m["scenarios"] if s.get("area") == ch
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

    def test_no_local_paths(self):
        """The v6 delta, item 6: no local file path in the document (the run directory, the suite's checkout, a
        home or temporary directory); the by-design table's path prints from metagraph/scripts/suite/ or as a
        file name."""
        local = re.compile(r"/Users/|/home/|~/|\\textasciitilde\{\}/|/var/folders/|/private/|/tmp/")
        for name in sorted(os.listdir(self.out)):
            if not name.endswith(".tex"):
                continue
            text = self.tex(name).replace("\\allowbreak{}", "").replace("\\_", "_")
            self.assertNotIn(self.tmp, text, name)
            self.assertNotIn(os.path.realpath(self.tmp), text, name)
            self.assertNotIn(SUITE, text, name)
            self.assertIsNone(local.search(text), "%s: %s" % (name, (local.search(text) or [None])[0]))

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
        # a member the service decided by design is said apart from one never exercised, with whether the host
        # was reached (READ-208; P-004's refusal came at submit)
        # the not-implemented count in the attention red (the owner's decision of 2026-10-16)
        self.assertEqual(text, "\\attn{not exercised: 1\\,\\gni{}} · by design: 1\\,\\gbd{} (decided before the host)")
        self.assertIsNone(line2)
        # one kind alone: the count and its glyph
        one = {"P-005": st["P-005"], "X": st["P-005"], "Y": st["P-005"]}
        self.assertEqual(cells.aggregate(one, "mcp")[0], "\\attn{not exercised: 3\\,\\gni{}}")
        # a mixed unit: the ok count over the exercised members, the rest said apart
        mixed = {"P-004": st["P-004"], "P-005": st["P-005"], "P-001": cells.scenario_states(self.runobj, "P-001")}
        self.assertTrue(cells.aggregate(mixed, "mcp")[0].startswith(
            "ok 1/1 · \\attn{not exercised: 1\\,\\gni{}} · by design: 1\\,\\gbd{} (decided before the host)"))
        # the coverage line under the table: words, never "exceptions: none" alone for an unexercised arm
        self.assertEqual(cells.coverage_text(st, "mcp"),
                         "exercised 0/2; \\attn{not exercised 1\\,\\gni{}}; by design 1\\,\\gbd{} (decided before the host)")
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

    def test_a_refused_table_without_compare_json(self):
        # L5 of the verification of SD-017: the document reads the comparison as report.py does (compare's
        # load_or_compute); a run without compare.json whose table compare.py refuses makes make_doc refuse with
        # compare.py's and report.py's exit 4, naming the refusal, never a traceback, and caches no comparison
        import contextlib
        import io
        d = os.path.join(self.tmp, "refused-run")
        shutil.copytree(self.run_dir, d)
        os.remove(os.path.join(d, "compare.json"))
        table = os.path.join(self.tmp, "refused-table.json")
        cases = [dict(c) for c in CASES.values()]
        for c in cases:
            if c["id"] == "case.refused":
                c["request"] = [{"path": "patterns", "present": True}]   # a request that is no object
        with open(table, "w") as fh:
            json.dump({"version": 1, "cases": cases}, fh)
        self.assertEqual(cmp.main([d, "--no-merge", "--by-design", table]), 4)
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            rc = make_doc.main([d, "--out", os.path.join(self.tmp, "refused-out"), "--no-pdf", "--units",
                                self.mapping_path, "--prose-dir", self.prose_dir, "--by-design", table])
        self.assertEqual(rc, 4)
        self.assertIn("make_doc: refusing: ", err.getvalue())
        self.assertIn("the service-by-design table %s is invalid" % table, err.getvalue())
        self.assertFalse(os.path.exists(os.path.join(d, "compare.json")))
        with self.assertRaises(reader.TableRefused):
            reader.Run(d, by_design=table)
        self.assertTrue(issubclass(reader.TableRefused, reader.RunError))


# ---------------------------------------------------------------------------------------------------------------
# the cell grammar and the expectation formatters

class CellGrammarTest(unittest.TestCase):
    def test_mcp_missing_only_where_planned(self):
        """The v6 delta (2026-10-17): compare's MCP class `missing` (no MCP execution) is a planned execution not
        recorded only where the run's plan schedules the scenario on the MCP arm; the direct-local warm-ups P-warm
        (timing class long) are clean there, not footnoted. A run without a plan keeps the comparator's reading."""
        class Stub(object):
            def __init__(self, planned):
                self.planned = planned

            def arm_lines(self, sid, arm):
                return []

            def comp_entry(self, sid):
                return {"mcp_agreement": {"class": "missing"}}

            def comp_arm(self, sid, arm):
                return {}

            def planned_missing(self, sid, arm):
                return False

            def planned_on(self, sid, arm):
                return self.planned

        st = cells.arm_state(Stub(False), "P-warm-0001", "mcp")
        self.assertFalse(st.present)
        self.assertFalse(st.missing_planned)
        self.assertTrue(st.clean)
        st = cells.arm_state(Stub(True), "P-warm-0001", "mcp")
        self.assertTrue(st.missing_planned)
        self.assertFalse(st.clean)
        stub = Stub(False)
        del Stub.planned_on
        self.assertTrue(cells.arm_state(stub, "P-warm-0001", "mcp").missing_planned)

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
        # a kind without a formatter: name(args), breakable after ( , = _ so that it wraps in a narrow column
        self.assertEqual(blocks.expect_text("walk_in_record", {"a": 1}),
                         "walk\\_\\allowbreak{}in\\_\\allowbreak{}record(\\allowbreak{}a=\\allowbreak{}1)")
        self.assertEqual(blocks.expect_text("walk_in_record", {}), "walk\\_\\allowbreak{}in\\_\\allowbreak{}record")
        self.assertEqual(blocks.expect_text("cut", {"reason": "max_paths"}),
                         "cut(\\allowbreak{}reason=\\allowbreak{}\"\\allowbreak{}max\\_\\allowbreak{}paths\"\\allowbreak{})")

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
        self.assertIn("refused before host execution · not expressible through MCP · answered by the service in the "
                      "host's place · answered after service normalization · host refusal relayed through the "
                      "service · answered by the host through the service", chapter)
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
        # the wall is a printed execution time above 10 s: on its tint (the budget plain)
        self.assertTrue(text.startswith("inconclusive: an allowed 503 deadline at \\tyellow{24 s} against a 20,000 ms "
                                        "budget"), text)
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
        self.assertIn("\\bcont{\\textbf{direct-remote}: DISAGREES with the reference access path}", block)
        self.assertNotIn("direct-remote}: pass", block)  # the row says it (VERIFY-R2 M2)
        self.assertIn("path / repeat & final outcome & http · code & wall & requested budget & property tested?", block)
        self.assertIn("direct-local & \\gfail\\ \\textbf{FAIL} & 200 & 20 ms (warm)", block)
        self.assertIn("direct-remote & \\gpass\\ pass & 200 & 0.1 s (warm) & none & yes (1 of 1)", block)
        self.assertIn("mcp r0 & \\ginc\\ inconclusive (by the checks) & 200", block)
        # the artifacts keyed once (READ-114): the one pattern with the (access path, repeat) keys it holds
        self.assertIn("\\bline{artifacts}{the raw answers \\texttt{arms/", block)
        self.assertIn("R-\\allowbreak{}001.\\allowbreak{}<access path>.\\allowbreak{}<repeat>.\\allowbreak{}json.\\allowbreak{}gz} for dl r0, dr r0, mcp r0", block)
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
        self.assertEqual(phases[0], "phase,n,n_zero,p10_ms,p25_ms,p50_ms,p75_ms,p90_ms,n_negative")
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
        flags = set(n.lstrip("!").split("=")[0].strip()
                    for n in re.findall(r"\[\[([^?\[\]]+)\?", json.dumps(closing.get("observations"))))
        self.assertEqual(flags - set(values), set(), "flags closing_values does not compute")
        for flag in flags:
            self.assertIsInstance(values[flag], bool, flag)
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
        self.assertIn("Timing 1\\qd Observations", main)
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
        self.assertEqual(sections._arm_reason(st.lines, st), "\ue003not exercised via MCP\ue002")
        self.assertEqual(latex.esc(sections._arm_reason(st.lines, st)), "\\attn{not exercised via MCP}")
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
        self.assertTrue(note.startswith("\\noindent{\\sffamily\\footnotesize\\gni\\ \\attn{Not exercised via MCP}: "), note)
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
        self.assertIn("\\Needspace{12\\baselineskip}\n\\parhead{Glossary}", text)  # readable size, two pages if need be (READ-020)
        self.assertIn("\\begin{multicols}{2}\\footnotesize", text)
        self.assertNotIn("fontsize{7.5}", text)
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
        # V6-001: the criterion names the policy's required population, its automatic exclusions included
        self.assertTrue(clause.startswith("2 of them block the deploy gate: those the policy requires on that access "
                                          "path and repeat"), clause)
        self.assertIn("the others are on an ungated access path, outside the selection, removed by the policy's "
                      "automatic exclusions or declared exempt", clause)

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
        # V6-006 (the external review of the v6 document): the values below 0 ms are counted apart too
        for name, p10, p25, p50, p75, p90, n, n_zero, n_negative in sections.phase_quantiles(ctx):
            key = dict((v_, k) for k, v_ in sections.PHASES)[name]
            self.assertEqual(n + n_zero + n_negative, len(pop[key]))
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
        self.assertEqual(text, "\\gni\\ \\attn{not exercised via MCP}: \\attn{refused before host execution} (http 400 \\texttt{resolve\\_\\allowbreak{}query\\_\\allowbreak{}too\\_\\allowbreak{}long}); the note under the table and case \\texttt{c.\\allowbreak{}ne} in B")
        text = blocks._not_exercised_text(s_, R(), "x", short, unit_note=False)
        self.assertTrue(text.startswith("\\gni\\ \\attn{not exercised via MCP}: the tool always sends its own bound; \\attn{refused before host execution}"), text)
        self.assertTrue(text.endswith("case \\texttt{c.\\allowbreak{}ne} in B"), text)
        rows = blocks._execution_rows(s_, R(), "x", "mcp", is_prop)
        self.assertEqual(rows[0][1], "\\gni\\ \\attn{not exercised}")
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
        self.assertEqual(heads, ["The three access paths.", "Execution order.", "Ground truth.", "Times and load."])
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
PROSE_FILES = tuple("%s.json" % ch for ch in WRITTEN_CHAPTERS) + ("overview.json", "glossary.json", "closing.json",
                                                                   "cases.json", "checks.json")
PROSE_SKIP_KEYS = ("_about", "placeholders", "_assumes", "_placeholders")
# what is not a figure of a run: placeholders, code spans, ids, dates, SPEC and section references, p50 / r0
# names, ordinals, appendix letters (the conditionals' heads go first, through strip_conditional_heads)
_NOT_A_FIGURE = [r"\{[^{}]*\}", r"`[^`]*`", r"\b[A-Z]+-[a-z_]+-\d{4}\b", r"(?<=-)\d{4}\b",
                 r"\d{4}-\d{2}-\d{2}", r"\bSPEC(?:-[A-Za-z-]+)?(?:\.md)?\s*(?:§\s*)?(?:\d+(?:\.\d+)?)?",
                 r"\bsections?\s+\d+(?:\.\d+)?(?:(?:,\s*|\s+and\s+|\s+to\s+)\d+(?:\.\d+)?)*", r"§\s*\d+(?:\.\d+)?",
                 r"\b[A-Z]\.\d\b", r"\b[pr]\d+\b", r"\b\d+(?:st|nd|rd|th)\b", r"\bappendi(?:x|ces)\s+[A-G]\b",
                 # the closing section's displayed names (Timing 1 .. 7; READ-009) are references, not figures
                 r"\bTiming\s+\d(?:\s*(?:to|and|,)\s*\d)*\b"]


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
        gap = whole.split("\\parhead{Coverage and remaining limits}", 1)[1].split("\\parhead{Background}", 1)[0]
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

    def test_prose_placeholders_resolvable(self):
        """Every placeholder and flag of a written chapter's prose names a value the generator computes: a chapter,
        family or unit key, a family of the chapter under `fam`, a unit of the chapter under `unit`, a scenario of
        the chapter under `x` with an arm and an arm key; a unit's property line reads the run-wide placeholders
        alone, a unit note under its table the unit's values. Over the v4 and the v5 mapping: a scenario named in
        the chapter prose or in the prose of a family both mappings hold is one both manifests hold, so that the
        same prose renders on a run of either (a sentence true in one manifest alone sits under a condition on
        `manifest_version`); a family one mapping alone holds (R-warm, T-longwalk) may name its own members."""
        import facts
        ctx = self.ctx({"R": self.R})
        chapter_keys = set(facts.chapter_facts(ctx, "R"))
        family_keys = set(facts.family_facts(ctx, "R-one"))
        unit_keys = set(facts.unit_facts(ctx, "R-one"))
        arm_keys = set(facts.arm_facts(ctx, "R-001", "mcp"))
        scenario_keys = set(facts.scenario_facts(ctx, "R-001"))
        global_keys = set(ctx.placeholders())
        loaded = {4: v4()[2], 5: v5()[2]}
        problems = []
        for ch in WRITTEN_CHAPTERS:
            data = sc.read_json(os.path.join(PROSE_DIR, "%s.json" % ch))
            sids = dict((v, set(u.chapter_scenarios(ch))) for v, u in loaded.items())
            fams = dict((v, set(u.chapter_families(ch))) for v, u in loaded.items())
            uids = dict((v, set(x for f in u.chapter_families(ch) for x in u.families[f]["units"])) for v, u in loaded.items())
            both_sids = sids[4] & sids[5]
            any_fams, any_uids = fams[4] | fams[5], uids[4] | uids[5]

            def members(fid):
                """The ids a family's prose may name: both manifests' when both hold the family, else that one's."""
                holders = [v for v in loaded if fid in fams[v]]
                if len(holders) == 2:
                    return both_sids
                return set(s for v in holders for u in loaded[v].families[fid]["units"] for s in loaded[v].units[u]["members"])

            def check(name, scope, allowed):
                root = name.split("[", 1)[0].split(":", 1)[0]
                idx = re.findall(r"\[([^\]]+)\]", name.split(":", 1)[0])
                where = "%s.json %s" % (ch, name)
                if root == "fam":
                    if scope != "chapter" or not idx or idx[0] not in any_fams:
                        problems.append("%s: fam[...] needs a family of %s in the chapter prose" % (where, ch))
                    elif len(idx) >= 2 and idx[1] != "x" and idx[1] not in family_keys:
                        problems.append("%s: no family value %s" % (where, idx[1]))
                    elif len(idx) >= 2 and idx[1] == "x" and (len(idx) < 3 or idx[2] not in members(idx[0])):
                        problems.append("%s: no scenario of %s in both manifests" % (where, idx[0]))
                elif root == "unit":
                    if scope != "chapter" or not idx or idx[0] not in any_uids:
                        problems.append("%s: unit[...] needs a unit of %s in the chapter prose" % (where, ch))
                    elif len(idx) >= 2 and idx[1] != "x" and idx[1] not in unit_keys:
                        problems.append("%s: no unit value %s" % (where, idx[1]))
                elif root == "x":
                    if not idx or idx[0] not in allowed:
                        problems.append("%s: no scenario of %s in both manifests" % (where, ch))
                    elif len(idx) == 2 and idx[1] not in scenario_keys:
                        problems.append("%s: no scenario value %s" % (where, idx[1]))
                    elif len(idx) >= 3 and (idx[1] not in ("dl", "dr", "mcp") or idx[2] not in arm_keys):
                        problems.append("%s: no arm value %s" % (where, idx[1:3]))
                elif root == "rel":
                    if scope != "family":
                        problems.append("%s: rel[...] is a family value" % where)
                else:
                    keys = {"chapter": chapter_keys, "family": family_keys, "unit": unit_keys, "global": set()}[scope]
                    if root not in keys and root not in global_keys:
                        problems.append("%s: no %s value %s" % (where, scope, root))

            for path, text in prose_strings(data):
                if path.startswith((".summary", ".coverage_gap", ".intro")):
                    scope, allowed = "chapter", both_sids
                elif path.startswith(".units"):
                    scope, allowed = "global", both_sids
                elif path.startswith(".mcp_notes"):
                    scope, allowed = "unit", both_sids
                elif path.startswith(".families."):
                    fid = path.split(".")[2]
                    scope, allowed = "family", members(fid)
                else:
                    continue
                for name in re.findall(r"\{([^{}]+)\}", text):
                    check(name, scope, allowed)
                for head in strip_conditional_heads(text)[1]:
                    check(re.split(r">?=", head.lstrip("!"), maxsplit=1)[0].strip(), scope, allowed)
            for fid in any_fams:
                entry = data["families"].get(fid) or {}
                if not entry.get("claim") and not any(p.startswith("%s.json: families.%s" % (ch, fid)) for p in ROUND2_PENDING):
                    problems.append("%s.json: %s has no claim" % (ch, fid))
                for _, wanted in facts.claim_items(entry):
                    for sid in wanted or []:
                        if sid not in members(fid):
                            problems.append("%s.json: claim member %s is not a member of %s" % (ch, sid, fid))
            for fid in data.get("coverage_gap_families") or []:
                if fid not in any_fams:
                    problems.append("%s.json: coverage_gap_families names %s" % (ch, fid))
        self.assertEqual(problems, [], "; ".join(problems))
        self.assertEqual(self.R.get("coverage_gap_families"), ["R-budget"])

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
                                         "[[fam[R-budget][n_late503] ? a | b]] appendix D 2026-10-15 90th Timing 3 to 7, "
                                         "Timing 1"), [])
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
        # the host state of a direct execution is said where it is not host-warm (VERIFY-PLXM M5)
        self.assertEqual(fv["x"]["P-001"]["dl"]["phrase"], "a pass in 9 ms (host-cold) / 6 ms (r0 / r1)")
        self.assertEqual(fv["x"]["P-001"]["dr"]["phrase"], "both repeats a pass in 0.1 s")
        self.assertEqual(fv["x"]["P-003"]["dl"]["phrase"], "the expected refusal (http 400) in 4 ms (host-cold)")
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
        rk:ch-<chapter>), every unit table whose variants name a record links to it and to appendix D's full key in
        its legend, and appendix D's key is labelled rk:D before the glossary. The chapter label is rk:ch-D for
        chapter D, not appendix D's rk:D (VERIFY-FIX2 L5: a label defined twice)."""
        import sections

        class R(object):
            def records(self):
                return [{"key": "eco", "label": "Escherichia coli K-12 MG1655", "accession": "NC_000913.3", "taxid": "511145"},
                        {"key": "pfa", "label": "Plasmodium falciparum 3D7 chr 5", "accession": "NC_004326.2", "taxid": "36329"}]

        class C(object):
            run = R()
        text = sections.record_key_paragraph(C())
        self.assertIn("\\label{rk:compact}", text)
        self.assertIn("\\label{rk:ch-P}", sections.record_key_paragraph(C(), "P"))
        self.assertEqual(sections.record_key_label("D"), "rk:ch-D")
        self.assertNotIn("\\label{rk:D}", sections.record_key_paragraph(C(), "D"))
        self.assertIn("\\texttt{eco} Escherichia coli K-12 MG1655 · \\texttt{pfa} Plasmodium falciparum 3D7 chr 5", text)
        self.assertIn("\\hyperref[rk:D]{appendix D}", text)
        ctx = self.ctx()
        # the fake run's ground truth has no records: no key, no link
        self.assertEqual(sections.record_key_paragraph(ctx), "")
        self.assertEqual(sections.record_key_link(ctx, ["P-001"]), "")
        ctx.run.manifest["ground_truth"] = {"records": [{"key": "bp", "label": "x", "accession": "a", "taxid": "1"}]}
        try:
            self.assertIn("record keys: \\hyperref[rk:ch-P]{p.\\,\\pageref{rk:ch-P}}", sections.record_key_link(ctx, ["P-001"], "P"))  # variant "10 bp"
            self.assertEqual(sections.record_key_link(ctx, ["R-001"], "R"), "")  # variant "the one"
            self.assertNotIn("Record key", sections.reading_guide(ctx))  # Part I stays four pages
            chapter = sections.chapter(ctx, "P")
            self.assertLess(chapter.find("\\parhead{Background}"), chapter.find("\\label{rk:ch-P}"))
            self.assertLess(chapter.find("\\label{rk:ch-P}"), chapter.find("\\unithead"))
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
        self.assertIn("unverified_note(unverified_detail)", body)  # grouped by reason, never cut (READ-018)
        self.assertNotIn("cut_text", src.split("def unverified_note(details):", 1)[1].split("\ndef ", 1)[0])
        # L-a: the limitations row no longer claims every inconclusive execution is footnoted in Part II
        self.assertIn("each footnoted in its chapter's tables", sections.limitations_table(ctx))
        glossary = sc.read_json(os.path.join(PROSE_DIR, "glossary.json"))
        inc = [e for e in glossary["entries"] if e["term"] == "inconclusive"][0]
        self.assertIn("not removed by the policy's automatic exclusions nor declared exempt", inc["definition"])
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
    number-unit space, the variant labels, the unit note's facts, the record key's wording. The final patch
    (VERIFY-R4): the closing section's conditionals on two fixtures, the banned timing wordings, the count grammar
    and the joined units outside the prose."""

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
        self.last_whole = text  # the whole chapter, for the evidence map after the conclusion (READ-101)
        conclusion = text.split("\\parhead{Conclusion}", 1)[1].split("\\begingroup", 1)[0]
        gap = text.split("\\parhead{Coverage and remaining limits}", 1)[1].split("\\parhead{Background}", 1)[0]
        family = text.split("fam:R-budget", 1)[1].split("\\unithead", 1)[0]
        return ctx, conclusion, gap, family

    def test_v4_prefix_relations_wording(self):
        """VERIFY-R5-2 issue 2 (the v4 wording error): the v4 branch of R-budget's prose, unit line and claim states the
        prefix_equals relations as the v4 manifest has them: REL-prefix_equals-0001 and -0003 compare the stop with its
        full twin (R-twin-0002 and -0004), -0002 and -0004 the continuation from base 5,000 (R-twin-0003 and -0005) with
        the full twin over the k-mers from that base; no sentence compares the stop against both twins."""
        ctx, conclusion, gap, family = self.render("perturbed")
        plain = family.replace("\\_", "_").replace("~", " ")
        self.assertIn("The full twin of each stop, the same 10 kb under 20,000 ms (R-twin-0002 and -0004, full unless stopped), "
                      "is what the \\texttt{prefix_equals} relations compare the stop against (REL-prefix_equals-0001 and -0003), "
                      "and the query's second half from base 5,000 (R-twin-0003 and -0005) is compared with the full twin over "
                      "the k-mers from that base (REL-prefix_equals-0002 and -0004): on their common k-mers", plain)
        for text in ("are what the", "its twins", "The two twins"):
            self.assertNotIn(text, plain, text)
        # the claim is in the evidence map after the conclusion's headline (READ-101)
        self.assertIn("its full twin agrees with", self.last_whole)
        self.assertNotIn("its twins agree with", self.last_whole)
        # the unit line and the claim of the real v4 prose, filled as the generator fills them (manifest_version 4)
        entry = self.R["units"]["R-budget"]  # {property, question} since the readability review (READ-005's slot)
        unit = facts.fill(entry["property"] if isinstance(entry, dict) else entry, {"manifest_version": 4})
        self.assertIn("its full 20 s twin agrees with it on their common k-mers, and the continuation from base 5,000 agrees "
                      "with the full twin on theirs (relations `prefix_equals`)", unit)
        self.assertNotIn("its twins", unit)
        self.assertEqual(facts.fill(self.R["families"]["R-budget"]["claim"], {"manifest_version": 4}),
                         "a 10 kb resolve under a small budget answers a stated prefix that its full twin agrees with")

    def test_two_fixtures_change_conclusions(self):
        """The owner's second review, point 4: chapter R rendered on two fixtures that differ in their verdicts
        changes its counts, its 'What held' branch, 'Shown' against 'Not shown', the stop sentence, 'no failure, no
        error and no disagreement', 'the run's only findings' and the family's statements, not only its numbers;
        no v4 figure and no unresolved placeholder other than the v4 scenarios the fixtures lack survive."""
        ctx_a, conc_a, gap_a, fam_a = self.render("plain")
        whole_a = self.last_whole
        ctx_b, conc_b, gap_b, fam_b = self.render("perturbed")
        whole_b = self.last_whole
        # the counts
        self.assertIn("0 of the 1 tests pass on direct-local and 1 on direct-remote", conc_a.replace("\\", ""))
        self.assertIn("1 of the 2 tests pass on direct-local and 1 on direct-remote", conc_b.replace("\\", ""))
        # the verdict branch of 'What held'
        self.assertIn("a direct access path failed, errored or disagreed with the other", conc_a)
        self.assertNotIn("no failure, no error and no disagreement", conc_a)
        self.assertIn("no failure, no error and no disagreement", conc_b)
        self.assertNotIn("a direct access path failed, errored or disagreed", conc_b)
        # the claims: unsupported with a disagreement, shown (in part) without
        self.assertIn("Not shown:", conc_a)
        self.assertIn("with a failure or a disagreement", conc_a)
        self.assertNotIn("Shown:", conc_a.replace("Not shown:", ""))
        self.assertIn("Every family's claim was shown, in full or in part, on the direct access paths.", conc_b)
        self.assertNotIn("Not shown:", conc_b)
        # the evidence map after the headline (READ-101): the claim with its family linked and the tests passing on
        # both direct access paths; a failing family reads not shown with the failure named
        self.assertIn("\\hyperref[fam:R-budget]{\\texttt{R-budget}} & a 10~kb resolve under a small budget answers a "
                      "stated prefix that its full twin agrees with & shown: 1 of 2 tests", whole_b)
        self.assertIn("\\hyperref[fam:R-budget]{\\texttt{R-budget}} & a 10~kb resolve under a small budget answers a "
                      "stated prefix that its full twin agrees with & not shown: 0 of 1 tests, with a failure or a "
                      "disagreement (I.1)", whole_a)
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
        self.assertIn("One of these 503s came more than max(1~s, the budget) after the budget and is a finding, "
                      "the chapter's and the run's only finding;", gap_b)
        self.assertNotIn("are findings", gap_b)
        # the footnoted tests of the family, singular
        self.assertIn("The one test of R-budget that did not pass on every access path (R-001) has its detail block", gap_a)
        self.assertIn("The one test of R-budget that did not pass on every access path (R-002) has its detail block", gap_b)
        # the family prose: the stop, the clean members, the MCP classes from the classifier
        self.assertIn("The stated stop has not been observed on this host", fam_a)
        self.assertNotIn("A stated stop was answered", fam_a)
        self.assertIn("A stated stop was answered (R-001 on direct-local (resolve stopped by time, 20~ms))", fam_b)
        # R-001 fails on direct-local without any 503 deadline: the 503 clause is not printed (VERIFY-R5 M1)
        self.assertIn("the only test did not pass on every access path (R-001); its detail block follows the table", fam_a)
        self.assertNotIn("503 deadlines past the budget", fam_a)
        self.assertNotIn("of the one tests", fam_a)
        self.assertIn("one of the two tests did not pass on every access path (R-002): the coverage gap on the chapter page "
                      "reads its 503 deadlines past the budget and its detail block follows the table", fam_b)
        # the clean members by id alone: the chronology stays in the detail blocks (READ-111)
        self.assertIn("the one other passed on every access path exercised: R-001", fam_b)
        self.assertNotIn("direct-local a pass in 20~ms; direct-remote a pass in 0.1~s", fam_b)
        self.assertNotIn("passed on every access path exercised", fam_a)
        for text in (fam_a, fam_b):
            self.assertIn("Via MCP: ", text)
            self.assertIn("1 answered by the host through the service", text)
            self.assertNotIn("every call was refused at submit", text)
        self.assertIn("; every test was exercised through the tools", conc_a)
        self.assertIn("\\attn{1 test could not be exercised through the tools} (a known limit: \\attn{the time budget 1}; "
                      "I.4 and appendix B)", conc_b)
        # the agreement's levels (READ-110): fully checked apart from the checks that ran alone
        self.assertIn("Via MCP 1 pair agrees on the checks both access paths evaluated (0 fully checked, every property "
                      "check evaluated on both; 1 on the generic checks or a refusal alone, which verifies no property)", conc_a)
        self.assertIn("Via MCP 1 pair agrees on the checks both access paths evaluated, every one fully checked (every "
                      "property check evaluated on both);", conc_b)
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

    def test_closing_two_fixtures_change_conclusions(self):
        """VERIFY-R4, issue 2: the closing section's interpretive clauses follow the computed values, shown on two
        fixtures whose conditions differ. On the perturbed fixture as built the median backend overhead is below
        half of the median MCP overhead, the paired medians print in seconds against milliseconds and the MCP arm
        has no cold-warm pairs; on a copy whose summary and walls are edited each clause flips: the host call
        carries half or more, both medians print in milliseconds, the cold-warm ratio is far from one and then
        close to one, and a negative MCP overhead prints neither backend clause."""
        import sections
        closing = sc.read_json(os.path.join(PROSE_DIR, "closing.json"))
        for key in ("mcp_over_local_positive", "backend_share_small", "paired_seconds_to_milliseconds",
                    "cold_warm_mcp_close"):
            self.assertIn("[[%s ?" % key, json.dumps(closing["observations"]))

        def render(d):
            run = reader.Run(d)
            run.load = {"n": 3, "sra_mean": 7000.0, "load_mean": 200.0, "sra_max": 1.0, "load_max": 1.0}
            units = units_mod.resolve(run.manifest, self.mappings["perturbed"])
            ctx = sections.Context(run, units, {}, {}, {}, list(units.chapters), closing=closing)
            obs = " ".join(sections.observations(ctx))
            self.assertNotIn("UNRESOLVED", obs)
            self.assertEqual(ctx.unresolved, [])
            return ctx, obs, sections.closing_values(ctx)

        ctx_a, obs_a, v_a = render(self.dirs["perturbed"])
        self.assertTrue(v_a["mcp_over_local_positive"])
        self.assertTrue(v_a["backend_share_small"])
        self.assertIn("that median is small beside the median MCP-minus-local difference of", obs_a)
        self.assertIn("two medians of separate per-scenario distributions, not a share of any request's overhead", obs_a)
        self.assertNotIn("at least half", obs_a)
        self.assertTrue(v_a["paired_seconds_to_milliseconds"])
        self.assertIn("the two sides' separate medians are 0.9~s against 6~ms, seconds against milliseconds, and their "
                      "quotient is not that ratio", obs_a)
        self.assertIn("the median of the per-scenario ratios", obs_a)
        self.assertFalse(v_a["cold_warm_mcp_close"])
        self.assertIn("observation cold-warm not computed for this run", obs_a)
        self.assertNotIn("quarter", obs_a)
        # the copy: the summary's backend overhead near the MCP overhead, an MCP cold-warm entry far from one,
        # the paired MCP wall in milliseconds
        d = os.path.join(self.tmp, "closing-flipped")
        shutil.copytree(self.dirs["perturbed"], d)
        comp = sc.read_json(os.path.join(d, "compare.json"))
        mo = comp["summary"]["mcp_overhead"]
        self.assertLess(mo["backend_minus_local_ms"]["p50"], 0.5 * mo["minus_local_ms"]["p50"])
        mo["backend_minus_local_ms"]["p50"] = 0.9 * mo["minus_local_ms"]["p50"]
        comp["summary"]["cold_warm"]["mcp"] = {"n": 3, "delta_ms": {"p50": 2000.0, "p90": 3000.0},
                                               "ratio": {"p50": 2.5, "p90": 3.0}}
        sc.write_json_atomic(os.path.join(d, "compare.json"), comp)
        path = os.path.join(d, "executions.jsonl")
        rows, _ = sc.read_jsonl(path)
        os.remove(path)
        for r in rows:
            if r["arm"] == "mcp" and r["scenario"] == "P-001":
                r["wall_ms"] = 50.0
            sc.append_jsonl(path, r)
        ctx_b, obs_b, v_b = render(d)
        self.assertFalse(v_b["backend_share_small"])
        self.assertIn("that median is at least half of the median MCP-minus-local difference of", obs_b)
        self.assertNotIn("is small beside", obs_b)
        self.assertFalse(v_b["paired_seconds_to_milliseconds"])
        self.assertIn("the two sides' separate medians are 50~ms against 6~ms, and their quotient", obs_b)
        self.assertNotIn("seconds against milliseconds", obs_b)
        self.assertFalse(v_b["cold_warm_mcp_close"])
        self.assertIn("(the median paired ratio warm over cold 2.50, that is +150~\\% off one; n 3)", obs_b)
        self.assertIn("lies outside 0.75 to 1.25 in this run, whatever the cause", obs_b)
        self.assertNotIn("differ little", obs_b)
        # the ratio close to one
        comp["summary"]["cold_warm"]["mcp"]["ratio"]["p50"] = 1.05
        sc.write_json_atomic(os.path.join(d, "compare.json"), comp)
        _, obs_c, v_c = render(d)
        self.assertTrue(v_c["cold_warm_mcp_close"])
        self.assertIn("lies between 0.75 and 1.25, so the recorded cold and warm executions differ little in this run, "
                      "whatever the cause", obs_c)
        self.assertNotIn("quarter of one", obs_c)
        self.assertIn("that is +5~\\% off one", obs_c)
        self.assertNotIn("lies outside 0.75 to 1.25", obs_c)
        # a negative MCP overhead prints neither backend clause
        mo["minus_local_ms"]["p50"] = -5.0
        sc.write_json_atomic(os.path.join(d, "compare.json"), comp)
        _, obs_d, v_d = render(d)
        self.assertFalse(v_d["mcp_over_local_positive"])
        self.assertFalse(v_d["backend_share_small"])
        self.assertNotIn("MCP-minus-local difference", obs_d)
        self.assertIn("took −5~ms longer", obs_d)
        # the whole closing section renders on both without an unresolved placeholder or a banned wording
        for ctx in (ctx_a, ctx_b):
            text = sections.closing(ctx)
            self.assertNotIn("UNRESOLVED", text)
            for bad in BANNED_TIMING_WORDINGS:
                self.assertNotIn(bad, text)

    def test_units_joined_outside_the_prose(self):
        """VERIFY-R4, issue 4: a number is joined to its unit in the findings list of I.1, in the summary lines
        under the tables and in appendix B's quoted texts, not only in the prose paragraphs."""
        import sections
        run = reader.Run(self.dirs["perturbed"])
        units = units_mod.resolve(run.manifest, self.mappings["perturbed"])
        ctx = sections.Context(run, units, {}, {}, {}, list(units.chapters))
        for ch in ctx.chapters:
            ctx.footnotes[ch] = cells.Footnotes()
        groups = dict(sections.principal_findings(ctx))
        items = " ".join(item for title in groups for item in groups[title])
        self.assertIn("5800~ms after the 2500~ms budget", items)
        self.assertNotRegex(items, r"\d (ms|s|bp)\b")
        line = sections.family_summary(ctx, "R-budget", ["R-001", "R-002"], [])
        self.assertRegex(line, r"\d~(ms|s)\b")
        self.assertNotRegex(line, r"\d (ms|s)\b")
        run.table["x-units"] = {"areas": ["R"], "agreement": "refused_before_host",
                                "scenario_match": "a query of 10001 bp",
                                "service_answer": {"text_contains": "longer than the server limit of 10000 bp"},
                                "direct_expectation": {"http": 400, "text": "answered within 3000 ms"}}
        table = sections.by_design_table(ctx)
        for joined in ("10001~bp", "10000~bp", "3000~ms"):
            self.assertIn(joined, table)
        self.assertNotIn("10000 bp", table)

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
        and the closing caveat. VERIFY-R4, issue 1: no causal timing reading in the closing section either (the
        banned wordings), its limits name what is measured and what is not isolated, and its interpretive
        clauses are conditionals on computed flags (issue 2)."""
        with open(os.path.join(DOC, "sections.py"), encoding="utf-8") as fh:
            src = fh.read()
        for bad in BANNED_TIMING_WORDINGS:
            self.assertNotIn(bad, src, "sections.py: %r" % bad)
        self.assertEqual(src.count("observed under concurrent host load"), 3)
        for name in PROSE_FILES:
            text = json.dumps(sc.read_json(os.path.join(PROSE_DIR, name)))
            for bad in BANNED_TIMING_WORDINGS:
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
        obs = json.dumps(closing["observations"])
        for key in ("mcp_over_local_positive", "backend_share_small", "paired_seconds_to_milliseconds",
                    "cold_warm_mcp_close"):
            self.assertIn("[[%s ?" % key, obs)
        secs = json.dumps(closing["sections"])
        self.assertIn("the record does not isolate the service's queue, the host's load of that minute and the cache "
                      "from one another", secs)
        self.assertIn("the run does not isolate the loop's order from the host's load", secs)
        self.assertIn("records each answer's body size but does not isolate its effect", secs)
        self.assertIn("the run does not isolate the host's load from the service's queue in either", secs)
        self.assertIn("what the host had cached, the loop's order and the load are not isolated from one another", secs)
        self.assertIn("Nothing here attributes a slow execution to the load", secs)
        # n next to the medians it counts
        self.assertIn("{wall_p50_warm_dl} at the median on direct-local (n {n_warm_dl})", summary)
        self.assertIn("(n {n_server_warm_dl} · {n_server_warm_dr} · {n_server_warm_mcp})", summary)
        # the sentences that asserted the run are conditionals on computed facts
        self.assertIn("[[mcp_all_refused_before_host ?", self.R["families"]["R-budget"]["prose"])
        self.assertIn("[[rel[REL-same_answer-0003][mcp] = pass ?", self.R["families"]["R-twin"]["prose"])
        self.assertNotIn("regularly outranked", self.R["families"]["R-discover"]["prose"])
        self.assertIn("[[mcp_all_refused_before_host ?", self.R["mcp_notes"]["units"]["R-budget"])
        self.assertNotIn("{members_text}. ", self.R["families"]["R-budget"]["prose"])
        self.assertIn("{clean_ids}", self.R["families"]["R-budget"]["prose"])
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
        self.assertIn("1 test, \\attn{1 ", note)  # the class of a member standing ○ is a red count
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

    def test_glossary_stated_stop_is_conditional(self):
        """The v5 merge (2026-10-15): appendix D's entry 'stated stop, prefix' states from the record whether the run
        answered one (facts.resolve_stops over every /resolve scenario on every arm, the run-wide placeholders
        n_resolve_stops, resolve_stop_observed and resolve_stop_ids), never as a fixed sentence: the plain fixture
        has none, the perturbed one R-001's stop on direct-local; a /traverse seed stopped by its time budget and a
        503 deadline are no stated stop of /resolve."""
        import sections
        glossary = sc.read_json(os.path.join(PROSE_DIR, "glossary.json"))
        entry = [e for e in glossary["entries"] if e["term"] == "stated stop, prefix"][0]
        self.assertIn("[[resolve_stop_observed ? Observed in this run ({resolve_stop_ids}). | Not observed in this run.]]",
                      entry["definition"])
        self.assertNotIn("met it. Not observed in this run.", json.dumps(glossary))
        rendered = {}
        for name in ("plain", "perturbed"):
            run = reader.Run(self.dirs[name])
            units = units_mod.resolve(run.manifest, self.mappings[name])
            ctx = sections.Context(run, units, {}, {}, glossary, list(units.chapters))
            values = ctx.placeholders()
            text = sections.glossary(ctx)
            self.assertNotIn("UNRESOLVED", text)
            self.assertEqual(ctx.unresolved, [])
            rendered[name] = (values, text)
        values, text = rendered["plain"]
        self.assertEqual((values["n_resolve_stops"], values["resolve_stop_observed"], values["resolve_stop_ids"]),
                         (0, False, "none"))
        self.assertIn("fails the checks that met it. Not observed in this run.", text)
        self.assertNotIn("Observed in this run (", text)
        values, text = rendered["perturbed"]
        self.assertEqual((values["n_resolve_stops"], values["resolve_stop_observed"], values["resolve_stop_ids"]),
                         (1, True, "R-001"))
        self.assertIn("fails the checks that met it. Observed in this run (R-001).", text)
        self.assertNotIn("Not observed in this run", text)
        # the fact itself: /resolve only, every arm, the executions counted and the scenarios named once
        lines = {("t", "direct-local"): [{"http": 200, "time_limited": "seed 0 stopped by its time budget"}],
                 ("r", "direct-local"): [{"http": 200, "time_limited": "resolve stopped by time"}],
                 ("r", "mcp"): [{"http": 200, "time_limited": "resolve stopped by time"}],
                 ("s", "direct-local"): [{"http": 503, "time_limited": "a 503 deadline: nothing written within the time budget"}],
                 ("u", "direct-local"): [{"http": 200, "time_limited": "resolve stopped by time"}, {"http": 200, "time_limited": ""}]}
        routes = {"t": "POST /traverse", "r": "POST /resolve", "s": "POST /resolve", "u": "POST /resolve"}
        stub = _Stub(order=["t", "r", "s", "u"], scenario=lambda sid: {"route": routes[sid]},
                     arm_lines=lambda sid, arm: lines.get((sid, arm), []))
        self.assertEqual(facts_mod().resolve_stops(stub),
                         {"n_resolve_stops": 3, "resolve_stop_observed": True, "resolve_stop_ids": "r, u"})
        stub.order = ["t", "s"]
        self.assertEqual(facts_mod().resolve_stops(stub),
                         {"n_resolve_stops": 0, "resolve_stop_observed": False, "resolve_stop_ids": "none"})

    def test_closing_load_clauses_are_conditional(self):
        """The v5 merge (2026-10-15): C.1's caveat 'observed under concurrent host load' and C.3's 'the MetaGraph
        host carried the Logan load' print only when the sampler's means make the Logan caveat true (logan_caveat),
        as the closing section's other interpretive clauses do; with the caveat false the sentences still read."""
        import sections
        closing = sc.read_json(os.path.join(PROSE_DIR, "closing.json"))
        self.assertIn("is the shared host's[[logan_caveat ?, observed under concurrent host load]]: not a baseline",
                      closing["caveats"][0])
        self.assertIn("[[logan_caveat ? The service ran on its own host and the MetaGraph host carried the Logan load "
                      "(Timing 6), so the overhead is this run's, not the service's floor. | The service ran on its own host, "
                      "and the MetaGraph host ran under the load Timing 6 records, so", closing["sections"]["C3"]["paragraphs"][2])

        def render(load):
            run = reader.Run(self.dirs["perturbed"])
            run.load = load
            units = units_mod.resolve(run.manifest, self.mappings["perturbed"])
            ctx = sections.Context(run, units, {}, {}, {}, list(units.chapters), closing=closing)
            values = sections.closing_values(ctx)
            caveat = ctx.fill(closing["caveats"][0], values)
            c3 = sections.closing_section_prose(ctx, "C3", "")
            self.assertNotIn("UNRESOLVED", caveat + c3)
            self.assertEqual(ctx.unresolved, [])
            return ctx.placeholders()["logan_caveat"], caveat, c3

        busy, caveat, c3 = render({"n": 3, "sra_mean": 7000.0, "load_mean": 200.0, "sra_max": 1.0, "load_max": 1.0})
        self.assertTrue(busy)
        self.assertIn("This run carries the Logan caveat", caveat)
        self.assertIn("Every timing in this section is the shared host's, observed under concurrent host load: not a "
                      "baseline for the release.", caveat)
        self.assertIn("The service ran on its own host and the MetaGraph host carried the Logan load (Timing 6), so the "
                      "overhead is this run", c3)
        self.assertNotIn("under the load Timing 6 records", c3)
        idle, caveat, c3 = render({"n": 3, "sra_mean": 10.0, "load_mean": 2.0, "sra_max": 1.0, "load_max": 1.0})
        self.assertFalse(idle)
        self.assertIn("The Logan caveat does not apply to this run", caveat)
        self.assertIn("Every timing in this section is the shared host's: not a baseline for the release.", caveat)
        self.assertNotIn("concurrent host load", caveat)
        self.assertIn("The service ran on its own host, and the MetaGraph host ran under the load Timing 6 records, so the "
                      "overhead is this run", c3)
        self.assertNotIn("carried the Logan load", c3)


def facts_mod():
    import facts
    return facts


# ---------------------------------------------------------------------------------------------------------------
# the rebase onto the suite at 07dfa487 (F-gate, the frozen v5, the oracle): the verdict and the gate line by
# report.py's code, the newer run fields read without breaking v4

class RebasedFieldsTest(unittest.TestCase):
    """The document states the report's verdict and gate line by the same code as report.py (verdict_text,
    verdict_headline, gate_headline_parts): a run with a gate verdict labels the run verdict exploratory; a stale
    recorded gate.json is named and the gate judged offline. The newer run fields of merge_run / compare / the v5
    manifest are read: warm-ups run again for a resume (warm_up_reruns, reruns_for_resume), the comparison's
    executions_sha256, the table's explained items, the property-level marks (a marked check reads as not
    implemented via MCP for that property, never as untested) and the scenario-level mark. v4 runs unchanged."""

    @classmethod
    def setUpClass(cls):
        import gate_policy
        from test_compare_report import T0, line, scen, write_arm
        cls.tmp = tempfile.mkdtemp()
        cls.run_dir = os.path.join(cls.tmp, "run")
        os.makedirs(cls.run_dir)
        manifest = build_run(cls.run_dir)
        # v5-shaped additions: P-001 marks contains_source_context not implemented via MCP; a warm-up W-001 of
        # P-001 whose re-run for a resume precedes P-001's direct-local execution; S-002 reads every property
        # check not_exposed via MCP under the marks (a scenario-level mark)
        scenarios = manifest["scenarios"]
        p1 = next(s for s in scenarios if s["id"] == "P-001")
        p1["expect"].append({"check": "source_label", "args": {}})
        p1["mcp"]["not_implemented_checks"] = [{"check": "contains_source_context",
                                                "reason": "the MCP's records carry no context offsets"}]
        warm = scen("W-001", "P")
        warm["warms"] = "R-001"
        warm["tags"] = ["warm-up"]
        whole = scen("S-002", "S")
        whole["mcp"]["not_implemented_checks"] = [{"check": "contains_source_context", "reason": "no contexts"}]
        scenarios.extend([warm, whole])
        manifest["content_sha256"] = sc.scenarios_sha256(scenarios)
        for arm in ("direct-local", "direct-remote", "mcp"):
            sc.write_json_atomic(os.path.join(cls.run_dir, "arms", arm, "manifest.json"), manifest)
        not_exposed = [{"check": "status", "outcome": "pass", "detail": "ok"},
                       {"check": "contains_source_context", "outcome": "not_exposed", "detail": "no field"}]
        t = T0 - 20
        extra_local = [line("W-001", "direct-local", 0, t, 5.0),
                       dict(line("W-001", "direct-local", 0, T0 + 4.5, 4.0), rerun_for_resume=True,
                            rerun_target="R-001"),
                       line("S-002", "direct-local", 0, T0 + 20, 5.0)]
        extra_mcp = [line("W-001", "mcp", 0, t + 5, 50.0), line("W-001", "mcp", 1, t + 7, 40.0),
                     line("S-002", "mcp", 0, T0 + 52, 60.0, checks=not_exposed),
                     line("S-002", "mcp", 1, T0 + 54, 60.0, checks=not_exposed)]
        # P-001's MCP lines read the marked check not_exposed beside a passing status and a passing second property
        p1_checks = not_exposed + [{"check": "source_label", "outcome": "pass", "detail": "ok"}]
        write_arm(cls.run_dir, "direct-local", extra_local)
        write_arm(cls.run_dir, "mcp", extra_mcp + [line("P-001", "mcp", 0, T0 + 50, 900.0, checks=p1_checks),
                                                   line("P-001", "mcp", 1, T0 + 52, 700.0, checks=p1_checks)])
        merge_run.merge(cls.run_dir)
        sc.write_json_atomic(os.path.join(cls.run_dir, "compare.json"), cmp.compare(cls.run_dir, (CASES, {"cases": 5})))
        cls.manifest = manifest
        mapping = json.loads(json.dumps(FAKE_MAPPING))
        mapping["_units"]["S-one"]["members"].append("S-002")
        mapping["_units"]["P-pos-reps"]["members"].append("W-001")
        mapping["S-002"] = {"chapter": "S", "family": "S-one", "unit_id": "S-one", "unit_type": "list", "variant": "whole"}
        mapping["W-001"] = {"chapter": "P", "family": "P-pos", "unit_id": "P-pos-reps", "unit_type": "replicates",
                            "variant": "w", "group": "a"}
        cls.mapping_path = os.path.join(cls.tmp, "units.json")
        with open(cls.mapping_path, "w") as fh:
            json.dump(mapping, fh)
        cls.prose_dir = os.path.join(cls.tmp, "prose")
        run = reader.Run(cls.run_dir)
        make_doc.scaffold_prose(cls.prose_dir, units_mod.resolve(run.manifest, mapping), run)
        cls.policy_path = os.path.join(cls.tmp, "policy.json")
        with open(cls.policy_path, "w") as fh:
            json.dump(gate_policy.build_policy(manifest, table=(CASES, {"cases": 5}), manifest_file="m.json",
                                               arms=("direct-local", "mcp"), warm_every=0), fh)
        cls.runobj = run

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def build(self, name, run_dir=None, extra=()):
        out = os.path.join(self.tmp, name)
        rc = make_doc.main([run_dir or self.run_dir, "--out", out, "--no-pdf", "--units", self.mapping_path,
                            "--prose-dir", self.prose_dir] + list(extra))
        self.assertEqual(rc, 0)
        with open(os.path.join(out, "doc.tex"), encoding="utf-8") as fh:
            main = fh.read()
        chapters = {}
        for ch in ("P", "R", "M", "S"):
            with open(os.path.join(out, "chapter_%s.tex" % ch), encoding="utf-8") as fh:
                chapters[ch] = fh.read()
        return main, chapters

    def test_verdict_and_gate_by_report_code(self):
        import gate
        run = self.runobj
        # ungated: the headline is the verdict; the report's word is the document's
        self.assertFalse(run.gated)
        self.assertEqual(run.verdict_headline, run.verdict)
        self.assertEqual(run.verdict, report.verdict_text(run.run, run.comp, run.rows, run.reevaluated,
                                                          (run.missing, run.plan_source), run.manifest, run.cases)[0])
        main, _ = self.build("ungated")
        self.assertNotIn("Exploratory verdict", main)
        self.assertIn("\\textbf{\\large Deploy gate: not evaluated}", main)
        # a recorded gate.json: its words are report.gate_headline_parts', the run verdict labelled exploratory
        gated = os.path.join(self.tmp, "gated-run")
        shutil.copytree(self.run_dir, gated)
        policy = gate.load_policy(self.policy_path)
        verdict = gate.judge(gated, policy)
        sc.write_json_atomic(os.path.join(gated, "gate.json"), verdict)
        run2 = reader.Run(gated)
        self.assertTrue(run2.gated)
        self.assertFalse(run2.gate_offline)
        self.assertEqual(run2.gate_stale, [])
        self.assertEqual(run2.verdict_headline, report.verdict_headline(run2.verdict, True, "beside"))
        self.assertTrue(run2.verdict_headline.startswith(report.exploratory_label("beside") + ": "))
        word, how, summary = report.gate_headline_parts(verdict)
        main, _ = self.build("gated", gated)
        cover = main.split("\\coverbox{Run verdict}{", 1)[1].split("\n", 1)[0]
        self.assertIn(latex.esc(sections.EXPLORATORY_LABEL), cover)
        self.assertNotIn("which is the gate's beside", cover)
        gate_box = main.split("\\coverbox{Deploy gate}{", 1)[1].split("\n", 1)[0]
        self.assertIn("\\textbf{\\large %s}" % word, gate_box)
        self.assertIn(latex.esc(how), gate_box)
        self.assertIn("recorded in the run", gate_box)
        # the kinds and the access paths once each on the cover and in I.1's gate line
        self.assertEqual(gate_box.count("by the gate's kinds"), 1)
        self.assertEqual(gate_box.count("access paths"), 1, gate_box)
        word2, how2, counts2, problems2 = sections.gate_summary(_Stub(run=run2))
        line = sections.gate_text(word2, how2, counts2, problems2, run2)
        self.assertEqual(line.count("access paths"), 1, line)
        kinds = line.split("by the gate's kinds: ", 1)[1].split(")", 1)[0]
        self.assertEqual(line.count(kinds), 1, line)
        self.assertNotIn("problems (", line)
        # the report's headline uses the same parts
        text = report.build_report(gated)
        self.assertIn("**%s** (%s): %s." % (report.md(word), report.md(how), report.md(summary)), text)
        self.assertIn(report.exploratory_label("above"), text)
        # a stale recorded verdict (the run changed after it): named, and the gate judged offline
        stale = os.path.join(self.tmp, "stale-run")
        shutil.copytree(gated, stale)
        rows, _ = sc.read_jsonl(os.path.join(stale, "executions.jsonl"))
        os.remove(os.path.join(stale, "executions.jsonl"))
        for r in rows[1:]:
            sc.append_jsonl(os.path.join(stale, "executions.jsonl"), r)
        run3 = reader.Run(stale, policy=self.policy_path)
        self.assertEqual(run3.gate_stale, ["executions.jsonl"])
        self.assertTrue(run3.gate_offline)
        self.assertTrue(run3.gated)
        main, _ = self.build("stale", stale, ["--policy", self.policy_path])
        self.assertIn("a stale gate.json", main)
        self.assertIn("judged other contents of executions.jsonl than the run has now", main)
        self.assertIn("judged offline by gate.py for this document", main)

    def test_reruns_and_the_comparison_lines(self):
        run = self.runobj
        rows, _ = sc.read_jsonl(os.path.join(self.run_dir, "executions.jsonl"))
        self.assertFalse(any(sc.is_rerun(r) for r in rows), "a re-run is never a line of the merged set")
        self.assertEqual(len(run.reruns), 1)
        self.assertEqual(run.reruns[0]["target_execution"], ["R-001", "direct-local", 0])
        self.assertEqual(run.rerun_counts(), {"direct-local": 1})
        target = run.arm_lines("R-001", "direct-local")[0]
        self.assertEqual([r["scenario"] for r in target["warm_up_reruns"]], ["W-001"])
        self.assertEqual(run.compare_executions_sha256, run.executions_sha256)
        main, chapters = self.build("fields")
        self.assertIn("warm-ups run again", main)
        self.assertIn("Warm-ups run again for a resume", main)
        self.assertIn("the lines this document reads", main)
        # the rerun is timing evidence in its target's block, never an execution of its own
        index = main.split("Index of test ids", 1)[1].split("\\end{multicols}", 1)[0]
        self.assertEqual(index.count(latex.tt("W-001")), 1)
        self.assertIn("warm-up re-run (dl", chapters["R"])
        self.assertIn("ran again for a resume", chapters["R"])
        # a v4 run (no such fields) reads as before

    def test_property_level_marks_read_as_the_mark(self):
        run = self.runobj
        # P-001: the marked check read not_exposed is set aside; the execution passes (suite_common.unmarked_checks)
        mcp0 = run.arm_lines("P-001", "mcp")[0]
        self.assertEqual(mcp0["outcome"], "pass")
        aside = run.marked_aside("P-001", mcp0)
        self.assertEqual(aside, {("contains_source_context", 0): {"mark": "not_implemented",
                                                                 "reason": "the MCP's records carry no context offsets"}})
        self.assertEqual(run.marked_aside("P-001", run.arm_lines("P-001", "direct-local")[0]), {})
        self.assertEqual(cells.property_tested(mcp0, run.is_property, aside), ("yes", 1, 1))
        self.assertEqual(cells.property_tested(mcp0, run.is_property), ("partly", 1, 2))
        self.assertEqual(blocks._property_cell(mcp0, run.is_property, aside), "yes (1 of 1); 1 marked not shown via MCP")
        # in the detail block's rows that long cell wraps in a top-aligned box (the v6 delta: D-P-scope-0004's MCP
        # rows ran 44 pt past the margin); a row without a mark keeps its plain cell
        rows = blocks._execution_rows(cells.arm_state(run, "P-001", "mcp"), run, "P-001", "mcp", run.is_property)
        self.assertEqual(rows[0][-1], "\\parbox[t]{30mm}{\\raggedright\\strut yes (1 of 1); 1 marked not shown via MCP}")
        rows = blocks._execution_rows(cells.arm_state(run, "P-001", "direct-local"), run, "P-001", "direct-local",
                                      run.is_property)
        self.assertNotIn("parbox", rows[0][-1])
        groups = blocks._check_groups(mcp0, True, run.check_marks("P-001"))
        self.assertEqual(groups, [(("marked:not_implemented", "the MCP's records carry no context offsets"),
                                   ["contains_source_context"])])
        self.assertIn("\\attn{not implemented via MCP} (the manifest's property-level mark)", blocks._group_text(*groups[0]))
        st = cells.arm_state(run, "P-001", "mcp")
        self.assertEqual(sections._arm_reason(mcp0 and [mcp0], st, "contains_source_context", run, "P-001"),
                         "\ue003not implemented via MCP\ue002 (the manifest's property-level mark: the MCP's records carry no "
                         "context offsets)")
        # the untested-properties lines never read the marked check as untested on the MCP arm: the direct arms
        # passed it, so P-001 is not listed
        main, _ = self.build("marks")
        self.assertNotIn("P-pos-reps} \\texttt{contains", main)
        # S-002: every property check marked and read not_exposed: a scenario-level mark, inconclusive, named
        s0 = run.arm_lines("S-002", "mcp")[0]
        self.assertEqual((s0["outcome"], s0.get("inconclusive_reason")), ("inconclusive", sc.SCENARIO_LEVEL_MARK))
        self.assertTrue(blocks.scenario_level_mark(s0))
        self.assertFalse(blocks._overridden(s0, True))
        self.assertEqual(blocks.mcp_class([s0], run.scenario("S-002"), {}, run.cases)[0], "marked")
        self.assertIn("scenario-level mark", blocks.adjudication(s0, run, "S-002", True))
        self.assertIn("scenario-level mark", facts.exec_phrase(run, "S-002", s0, "mcp"))
        self.assertEqual(sections._arm_reason([s0], cells.arm_state(run, "S-002", "mcp"), "contains_source_context",
                                              run, "S-002"),
                         "every property check marked as not shown via MCP (a scenario-level mark)")
        self.assertEqual(run.marked_aside("S-002", s0), {})
        self.assertEqual(blocks._check_groups(s0, True, run.check_marks("S-002")),
                         [(("marked:not_implemented", "no contexts"), ["contains_source_context"])])
        _, chapters = self.build("whole")
        self.assertIn("scenario-level mark", chapters["S"])
        self.assertIn("\\attn{not implemented via MCP} (the manifest's property-level mark): no contexts", chapters["S"])

    def test_explained_items_from_the_table(self):
        run = self.runobj
        # P-003: answered by the service, the table's case.answered judged by the checks: compare's explained items
        items, case = run.explained_items("P-003", 0)
        ma = run.comp_entry("P-003").get("mcp_agreement") or {}
        self.assertEqual(items, list((ma.get("explained_items") or {}).get("0") or []))
        self.assertEqual(run.explained_items("P-001", 0), ([], None))

    def test_mapping_by_manifest_sha(self):
        manifest = sc.read_json(V4_MANIFEST)
        self.assertEqual(units_mod.default_mapping_path(manifest), V4_MAPPING)
        other = {"manifest_version": 9, "content_sha256": manifest["content_sha256"]}
        self.assertEqual(units_mod.default_mapping_path(other), V4_MAPPING)
        unknown = {"manifest_version": 7, "content_sha256": "0" * 64}
        self.assertEqual(os.path.basename(units_mod.default_mapping_path(unknown)), "suite_v7.json")


class V5ScaffoldTest(unittest.TestCase):
    """The generator on a FAKE v5-shaped run (tests/fake_doc_v5.py: every scenario of the frozen v5 manifest,
    synthetic answers): the mapping is picked by the manifest's content sha, every v5 id is in the index once, the
    v5-only units print, the manifest's property-level marks and the warm-up's re-run are read."""

    @classmethod
    def setUpClass(cls):
        from fake_doc_v5 import build_v5_run
        cls.tmp = tempfile.mkdtemp()
        cls.run_dir = os.path.join(cls.tmp, "run")
        cls.manifest = build_v5_run(cls.run_dir)
        cls.out = os.path.join(cls.tmp, "doc")
        cls.rc = make_doc.main([cls.run_dir, "--out", cls.out, "--no-pdf"])

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def read(self, name):
        with open(os.path.join(self.out, name), encoding="utf-8") as fh:
            return fh.read()

    def test_scaffold(self):
        self.assertEqual(self.rc, 0)
        run = reader.Run(self.run_dir)
        self.assertEqual(run.manifest["content_sha256"], sc.read_json(V5_MANIFEST)["content_sha256"])
        units = units_mod.resolve(run.manifest, units_mod.load_mapping(units_mod.default_mapping_path(run.manifest)))
        self.assertEqual(units.source, "mapping")
        self.assertEqual(len(units.units), 109)
        main = self.read("doc.tex")
        for ch in ("R", "T", "P", "L", "X", "M", "Q", "C", "N", "D"):
            self.assertIn("\\input{chapter_%s}" % ch, main)
        index = main.split("Index of test ids", 1)[1].split("\\end{multicols}", 1)[0]
        for sid in run.order:
            self.assertEqual(index.count(latex.tt(sid)), 1, sid)
        self.assertEqual(len(run.order), 700)
        chapter_r, chapter_t = self.read("chapter_R.tex"), self.read("chapter_T.tex")
        self.assertIn("\\label{unit:R-warm}", chapter_r)
        self.assertIn("\\label{unit:T-longwalk}", chapter_t)
        self.assertLess(chapter_t.index("T-\\allowbreak{}walk-\\allowbreak{}0006"),
                        chapter_t.index("T-\\allowbreak{}longwalk-\\allowbreak{}0001"))
        # the marked checks of the T-labels pairs read as set aside (compare: checks_marked), never as untested
        self.assertEqual(run.summary()["mcp_agreement"].get("checks_marked"), 8)
        self.assertEqual(run.marked_aside("T-labels-0003", run.arm_lines("T-labels-0003", "mcp")[0]),
                         {("label_dict_contains", 0): {"mark": "not_implemented",
                                                      "reason": "the MCP traverse answers carry no label dictionary"}})
        self.assertNotIn("label\\_dict\\_contains", main.split("Properties still untested", 1)[1].split("\\parhead", 1)[0])
        # the warm-up's re-run is timing evidence in the provenance and appendix F, not an execution
        self.assertEqual(len(run.reruns), 1)
        self.assertIn("Warm-ups run again for a resume", main)
        self.assertEqual(index.count(latex.tt("R-warm-0001")), 1)

    def test_chapter_r_v5_prose(self):
        """The v5 prose round of chapter R: the warm-ups' family and unit are written (no placeholder), the v5
        branches print (R-budget's budgets and twins, the service's 60,000 ms, the discovery cap, the junctions'
        headers, eco's reverse complement) and nothing of the v4 wording for the same facts; the warm-ups' order
        facts (facts._warm_ups) read the fake run: eight warm-ups before their targets on every arm, the two complete
        twins before their stops, one warm-up run again for a resume."""
        import sections
        run = reader.Run(self.run_dir)
        units = units_mod.resolve(run.manifest, units_mod.load_mapping(units_mod.default_mapping_path(run.manifest)))
        chapter_r = self.read("chapter_R.tex")
        self.assertNotIn(PLACEHOLDER, chapter_r)
        self.assertNotIn("UNRESOLVED CONDITION", chapter_r)
        for text in ("5,000~bp of E. coli under 480~ms and 10,000~bp of P. aeruginosa under 420~ms",
                     "The complete twin of each stop, the same query under 60,000~ms (R-twin-0002 and -0004, full unless "
                     "stopped, run before the stop), is what the",
                     "is compared with the complete twin over the k-mers from that base",
                     "its complete twin agrees with it on their common k-mers, and the continuation from the middle base "
                     "(2,500 and 5,000) agrees with the complete twin on theirs",
                     "answers a stated prefix that its complete twin agrees with",
                     "Each complete twin's first execution in this run started before its stop's on each direct access path that ran both (4 pairs of twin and access path compared)",
                     "neither stop (R-budget-0001 under 480~ms, R-budget-0002 under 420~ms, each run after its complete twin) nor the complete 60~s twins answered one",
                     "Each request names the 60,000~ms budget", "(70.25~s here)", "a discovery limit outside 1 to 1,000",
                     "with \\texttt{max\\_labels} 9, 3, 20 or 100",
                     "resolved under both records' headers and both taxid columns",
                     "the three discovery edges, R-edge-0022, -0023 and -0025, take its first 60~bp",
                     "(S. aureus, E. coli; eco takes the place of v4's M. jannaschii",
                     "In this run every warm-up's first execution started before its target's on each access path that ran both (24 pairs of warm-up and access path compared)",
                     "A resume after an interruption ran one warm-up again (R-warm-0001 on direct-local before T-repeat-0001)",
                     "under a label its target walks under (a record's header for the five repeat warm-ups, a taxid column for the three others)",
                     "a warm-up's query resolves fully within the budget under a label its target walks under",
                     "A warm-up of a traverse scenario of chapters T and N",
                     "within the service's discovery cap of 1,000 (100 at the v4 run)",
                     "that none is the other genome's is the index's answer"):
            self.assertIn(text, chapter_r, text)
        for text in ("10~kb of a record (E. coli under 3,000~ms", "the same 10~kb under 20,000~ms", "neither 10~kb stop",
                     "Each request names the 20,000~ms budget", "(30.25~s here)", "M. jannaschii), resolved",
                     "the two discovery edges, R-edge-0022 and -0023, take",
                     # VERIFY-R5: the stop is compared with its complete twin, not with the continuation; the 503 clause
                     # (no 503 on the fake) and the v4 wordings
                     "its twins agree with", "its twins, the complete resolve", "503 deadlines past the budget",
                     "since manifest v5", "never the other genome's", "or search scenario", "within its radius, plus a margin",
                     "then doubling to at most", "again right before the target (",
                     # VERIFY-R5-2 issue 5: T-labels-0003 walks under derived column labels, 511145 among them
                     "the label its target walks under", "the label the target walks under", "its target's label"):
            self.assertNotIn(text, chapter_r, text)
        prose = dict((ch, sc.read_json(os.path.join(PROSE_DIR, "%s.json" % ch))) for ch in ("R",))
        ctx = sections.Context(run, units, prose, {}, {}, ["R"])
        ctx.footnotes["R"] = cells.Footnotes()
        fw = facts.family_facts(ctx, "R-warm")
        self.assertEqual((fw["n_warm_ups"], fw["n_warm_up_arm_pairs"], fw["n_warm_ups_after_target"]), (8, 24, 0))
        self.assertTrue(fw["warm_ups_all_before_target"])
        self.assertEqual(fw["warm_ups_after_target_text"], "")
        self.assertEqual(fw["n_warm_up_reruns"], 1)
        self.assertEqual(fw["warm_up_reruns_text"], "R-warm-0001 on direct-local before T-repeat-0001")
        self.assertIn("R-warm-0007 (T-walk-0009); R-warm-0001 (T-repeat-0001)", fw["warm_up_targets_text"])
        fb = facts.family_facts(ctx, "R-budget")
        self.assertEqual((fb["n_warm_ups"], fb["n_warm_up_arm_pairs"], fb["n_warm_ups_after_target"]), (2, 6, 0))
        self.assertEqual(fb["warm_up_targets_text"], "R-twin-0002 (R-budget-0001); R-twin-0004 (R-budget-0002)")
        # the direct arms alone (VERIFY-R5 L4): the MCP arm's pair is left out of R-budget's order sentence
        self.assertEqual((fb["n_warm_up_direct_pairs"], fb["n_warm_ups_after_target_direct"]), (4, 0))
        self.assertTrue(fb["warm_ups_all_before_target_direct"])
        self.assertEqual((fw["n_warm_up_direct_pairs"], fw["n_warm_ups_after_target_direct"]), (16, 0))
        self.assertEqual(fw["n_warm_up_reruns_without_target"], 0)
        # no 503 on the fake: the footnoted-with-a-503 flag is false, nothing footnoted
        self.assertEqual((fb["n_footnoted"], fb["n_late503"], fb["n_footnoted_late503"]), (0, 0, 0))
        self.assertFalse(fb["footnoted_all_late503"])
        fp = facts.family_facts(ctx, "R-pos")
        self.assertEqual((fp["n_warm_ups"], fp["n_warm_up_arm_pairs"]), (0, 0))
        self.assertFalse(fp["warm_ups_all_before_target"])
        self.assertEqual(facts.unit_facts(ctx, "R-warm")["n_warm_ups"], 8)
        # a run whose warm-up started after its target: the fact says so, with the arm and the delay
        rows, _ = sc.read_jsonl(os.path.join(self.run_dir, "arms", "direct-remote", "executions.jsonl"))
        late = dict((r["scenario"], r) for r in rows)["R-warm-0003"]
        target = dict((r["scenario"], r) for r in rows)["T-repeat-0003"]
        moved = dict(late, started_utc=sc.iso_utc(sc.parse_utc(target["started_utc"]) + 12.0))
        orig = run.arm_lines
        run.arm_lines = lambda sid, arm: ([moved] if (sid, arm) == ("R-warm-0003", "direct-remote") else orig(sid, arm))
        fw2 = facts.family_facts(ctx, "R-warm")
        self.assertEqual((fw2["n_warm_up_arm_pairs"], fw2["n_warm_ups_after_target"]), (24, 1))
        self.assertFalse(fw2["warm_ups_all_before_target"])
        self.assertEqual(fw2["warm_ups_after_target_text"],
                         "R-warm-0003 on direct-remote (started 12~s after its target T-repeat-0003)".replace("~", " "))
        self.assertEqual((fw2["n_warm_up_direct_pairs"], fw2["n_warm_ups_after_target_direct"]), (16, 1))
        self.assertEqual(fw2["warm_ups_after_target_direct_text"], fw2["warm_ups_after_target_text"])
        self.assertEqual(fw2["warm_ups_after_target_direct_arms"], "direct-remote")
        # the complete twin after its stop on both direct arms, or on one (VERIFY-R5-2 issue 4): the sentence names them
        budget_prose = prose["R"]["families"]["R-budget"]["prose"]

        def late_twin(arm):
            stop = orig("R-budget-0002", arm)[0]
            return dict(orig("R-twin-0004", arm)[0], started_utc=sc.iso_utc(sc.parse_utc(stop["started_utc"]) + 12.0))
        run.arm_lines = lambda sid, arm: ([late_twin(arm)] if sid == "R-twin-0004" and arm != "mcp" else orig(sid, arm))
        fb2 = facts.family_facts(ctx, "R-budget")
        self.assertEqual((fb2["n_warm_up_direct_pairs"], fb2["n_warm_ups_after_target_direct"]), (4, 2))
        self.assertEqual(fb2["warm_ups_after_target_direct_arms"], "both direct access paths")
        self.assertIn("In this run two executions of the complete twins did not start before their stops' on both direct access paths "
                      "(R-twin-0004 on direct-local (started 12 s after its target R-budget-0002); R-twin-0004 on direct-remote "
                      "(started 12 s after its target R-budget-0002)), against the design.", facts.fill(budget_prose, fb2))
        run.arm_lines = lambda sid, arm: ([late_twin(arm)] if (sid, arm) == ("R-twin-0004", "direct-local") else orig(sid, arm))
        fb3 = facts.family_facts(ctx, "R-budget")
        self.assertEqual((fb3["n_warm_ups_after_target_direct"], fb3["warm_ups_after_target_direct_arms"]), (1, "direct-local"))
        self.assertIn("In this run one execution of a complete twin did not start before its stop's on a direct access path (R-twin-0004 "
                      "on direct-local (started 12 s after its target R-budget-0002)), against the design.",
                      facts.fill(budget_prose, fb3))
        run.arm_lines = orig
        # a tie is "started with its target", not "<1 ms after" (VERIFY-R5 L10)
        tied = dict(late, started_utc=target["started_utc"])
        run.arm_lines = lambda sid, arm: ([tied] if (sid, arm) == ("R-warm-0003", "direct-remote") else orig(sid, arm))
        fw3 = facts.family_facts(ctx, "R-warm")
        self.assertEqual(fw3["n_warm_ups_after_target"], 1)
        self.assertEqual(fw3["warm_ups_after_target_text"], "R-warm-0003 on direct-remote (started with its target T-repeat-0003)")
        # a re-run no execution of the target followed is said so, not "before the target"
        run.reruns = [dict(run.reruns[0], target_execution=None)]
        fw4 = facts.family_facts(ctx, "R-warm")
        self.assertEqual((fw4["n_warm_up_reruns"], fw4["n_warm_up_reruns_without_target"]), (1, 1))
        self.assertEqual(fw4["warm_up_reruns_text"], "R-warm-0001 on direct-local (no execution of its target T-repeat-0001 followed)")
        self.assertIn(", except 1 that no execution of the target followed", sections.reruns_text(run))

    def _budget_fixture(self, name, mode):
        """The fake v5 run with the stops' direct executions rewritten (re-merged, compared, chapter R rendered):
        `complete`, both stops complete within their budgets on the direct arms (resolve_runs inconclusive, no stop
        happened); `mixed`, R-budget-0001 completes and R-budget-0002 answers a 503 deadline 0.8 s past its budget on
        direct-local (inconclusive: under max(1 s, the budget)) and completes on direct-remote; `all_late`, both stops
        answer a 503 deadline on both direct arms, R-budget-0001 0.9 s past its 480 ms budget (allowed: under max(1 s,
        the budget), VERIFY-R5-3 issue 4) and R-budget-0002 2.5 s past its 420 ms budget (a finding), while the
        complete twins keep the fake's passes (VERIFY-R5-2 issue 1: v5's likely outcome, the stops late and the twins
        clean); `stop_error`, R-budget-0001 completes and R-budget-0002 answers an http 500 on both direct arms (an
        error, no 503); `over_budget`, R-budget-0001 completes and R-budget-0002 answers a 200 past its budget on
        both direct arms (budget_respected a finding, no 503); `late_and_over`, R-budget-0001 answers a 503 deadline
        1.3 s past its budget on both direct arms (findings) and R-budget-0002 a 200 past its budget on both (findings
        too; VERIFY-R5-3 issue 1). Returns (chapter R's text with ~ and allowbreak removed, the coverage gap's
        paragraph, R-budget's family facts)."""
        p = os.path.join(self.tmp, name)
        shutil.copytree(self.run_dir, p)
        complete = [{"check": "status", "outcome": "pass", "detail": "ok"},
                    {"check": "invariants", "outcome": "pass", "detail": "ok"},
                    {"check": "resolve_runs", "outcome": "inconclusive",
                     "detail": "no stop happened: the prefix property was not tested (label 511145: 4970 of 4970 k-mers, 1 run(s))"},
                    {"check": "budget_respected", "outcome": "pass", "detail": "ok"}]

        def late(r, budget_ms, ms, finding):
            err = {"code": "deadline", "text": "resolve: the answer could not be built and written within bounds.time_budget_ms "
                                               "(%d ms, the finalisation reserve of 250 ms included): nothing partial is sent" % budget_ms}
            checks = [{"check": "status", "outcome": "inconclusive", "detail": "http 503 deadline in place of 200: allowed", "property": False},
                      {"check": "invariants", "outcome": "inconclusive", "detail": "not testable: http 503", "property": False},
                      {"check": "resolve_runs", "outcome": "inconclusive", "detail": "not testable: http 503", "property": False},
                      {"check": "budget_respected", "outcome": "finding" if finding else "pass",
                       "detail": "503 deadline %d ms past the budget%s" % (ms, "" if finding else ": within max(1 s, the budget)")}]
            return dict(r, outcome="finding" if finding else "inconclusive", http=503, error=err,
                        wall_ms=float(budget_ms + ms), checks=checks)

        def errored(r):
            return dict(r, outcome="error", http=500, error={"code": "internal", "text": "internal error"}, wall_ms=300.0,
                        checks=[{"check": "status", "outcome": "fail", "detail": "http 500, expected 200", "property": False}])

        def over_budget(r, budget_ms, ms):
            return dict(r, outcome="finding", http=200, wall_ms=float(budget_ms + ms),
                        checks=[c if c["check"] != "budget_respected" else
                                {"check": "budget_respected", "outcome": "finding",
                                 "detail": "a 200 at a client wall time of %d ms, %d ms past the %d ms budget" % (budget_ms + ms, ms, budget_ms)}
                                for c in complete])

        for arm in ("direct-local", "direct-remote"):
            path = os.path.join(p, "arms", arm, "executions.jsonl")
            rows, _ = sc.read_jsonl(path)
            os.remove(path)
            for r in rows:
                if r["scenario"] == "R-budget-0001" and mode == "all_late":
                    r = late(r, 480, 900, False)
                elif r["scenario"] == "R-budget-0002" and mode == "all_late":
                    r = late(r, 420, 2500, True)
                elif r["scenario"] == "R-budget-0002" and mode == "stop_error":
                    r = errored(r)
                elif r["scenario"] == "R-budget-0002" and mode == "over_budget":
                    r = over_budget(r, 420, 900)
                elif r["scenario"] == "R-budget-0001" and mode == "late_and_over":
                    r = late(r, 480, 1300, True)
                elif r["scenario"] == "R-budget-0002" and mode == "late_and_over":
                    r = over_budget(r, 420, 900)
                elif r["scenario"] == "R-budget-0002" and mode == "mixed" and arm == "direct-local":
                    r = late(r, 420, 800, False)
                elif r["scenario"] in ("R-budget-0001", "R-budget-0002"):
                    r = dict(r, outcome="inconclusive", checks=complete)
                sc.append_jsonl(path, r)
        merge_run.merge(p)
        sc.write_json_atomic(os.path.join(p, "compare.json"), cmp.compare(p, cmp.load_by_design(cmp.DEFAULT_BY_DESIGN)))
        out = os.path.join(self.tmp, name + "-doc")
        self.assertEqual(make_doc.main([p, "--out", out, "--chapters", "R", "--no-pdf"]), 0)
        with open(os.path.join(out, "chapter_R.tex"), encoding="utf-8") as fh:
            tex = fh.read().replace("\\allowbreak{}", "").replace("~", " ")
        self.assertNotIn(PLACEHOLDER, tex)
        self.assertNotIn("UNRESOLVED CONDITION", tex)
        gap = tex.split("What remains unsupported.", 1)[1].split("\n", 1)[0]
        import sections
        run = reader.Run(p)
        units = units_mod.resolve(run.manifest, units_mod.load_mapping(units_mod.default_mapping_path(run.manifest)))
        prose = dict((ch, sc.read_json(os.path.join(PROSE_DIR, "%s.json" % ch))) for ch in ("R",))
        ctx = sections.Context(run, units, prose, {}, {}, ["R"])
        ctx.footnotes["R"] = cells.Footnotes()
        return tex, gap, facts.family_facts(ctx, "R-budget")

    def test_chapter_r_budget_without_503(self):
        """VERIFY-R5 M1: a v5 stop that completes within its budget is footnoted (inconclusive) without a 503 deadline,
        so R-budget's run sentence must not say the coverage gap reads 503 deadlines; the coverage gap says the
        stops completed. With one footnoted member answering a 503 and the other not, the run sentence names it (the
        noun by the count of 503s, a comma before 'and their detail blocks'; VERIFY-R5-2 issue 3) and the coverage
        gap names the test that answered the 503, prints one value when the earliest and the latest lateness print
        alike, and says that the other stop completed (VERIFY-R5-2 issue 1)."""
        tex, gap, fb = self._budget_fixture("budget-complete", "complete")
        self.assertEqual((fb["n_footnoted"], fb["n_late503"], fb["n_footnoted_late503"]), (2, 0, 0))
        self.assertFalse(fb["footnoted_all_late503"])
        self.assertEqual((fb["late503_id_list"], fb["late503_ids_text"], fb["late503_range_text"]), ([], "", "–"))
        self.assertEqual((fb["n_stops_complete"], fb["stops_complete_text"]), (2, "R-budget-0001 and R-budget-0002"))
        self.assertIn("The stops completed within their budgets without a stop, so the property was not exercised.", gap)
        self.assertIn("two of the six tests did not pass on every access path (R-budget-0001, R-budget-0002); their detail "
                      "blocks follow the table; the four others passed on every access path exercised", tex)
        self.assertNotIn("503 deadlines past the budget", tex)
        self.assertNotIn("answered 503 deadline", tex)
        tex, gap, fb = self._budget_fixture("budget-mixed", "mixed")
        self.assertEqual((fb["n_footnoted"], fb["n_late503"], fb["n_footnoted_late503"]), (2, 1, 1))
        self.assertFalse(fb["footnoted_all_late503"])
        self.assertEqual((fb["late503_ids"], fb["late503_id_list"], fb["late503_ids_text"]),
                         ("R-budget-0002", ["R-budget-0002"], "R-budget-0002"))
        self.assertEqual((fb["late503_min"], fb["late503_max"], fb["late503_range_text"]), ("0.8 s", "0.8 s", "0.8 s"))
        self.assertEqual((fb["n_stops_complete"], fb["stops_complete_text"]), (1, "R-budget-0001"))
        self.assertIn("two of the six tests did not pass on every access path (R-budget-0001, R-budget-0002): the coverage gap on the "
                      "chapter page reads the 503 deadline past the budget of R-budget-0002, and their detail blocks follow "
                      "the table", tex)
        self.assertIn("nor the complete 60 s twins answered one. On this host", gap)
        self.assertIn("R-budget-0002 answered 503 deadline, 0.8 s past its budget (R-budget-0002 0.8 s on direct-local). "
                      "R-budget-0001 completed within its budget without a stop. So the prefix property is inconclusive "
                      "wherever it was asked", gap)
        for text in ("they answered 503 deadline", "0.8 s to 0.8 s", "their 503 deadlines past the budget",
                     "503 deadlines past the budget of", "The stops completed within their budgets",
                     "R-budget-0002 and their detail blocks"):
            self.assertNotIn(text, tex, text)

    def test_chapter_r_budget_all_late_twins_clean(self):
        """VERIFY-R5-2 issue 1: both stops answer 503 deadlines past their budgets on both direct arms while the complete
        twins pass everywhere (v5's likely outcome). The coverage gap names the two stops as the tests that answered
        503, with the range of their lateness, and claims no 503 of a twin; the run sentence reads 'their 503
        deadlines' for the two footnoted stops."""
        tex, gap, fb = self._budget_fixture("budget-all-late", "all_late")
        self.assertEqual((fb["n_footnoted"], fb["n_late503"], fb["n_late503_tests"], fb["n_footnoted_late503"]), (2, 4, 2, 2))
        self.assertTrue(fb["footnoted_all_late503"])
        self.assertEqual(fb["late503_id_list"], ["R-budget-0001", "R-budget-0002"])
        self.assertEqual((fb["late503_ids_text"], fb["late503_range_text"]), ("R-budget-0001 and R-budget-0002", "0.9 s to 2.5 s"))
        self.assertEqual((fb["n_stops_complete"], fb["stops_complete_text"], fb["n_finding_executions"]), (0, "", 2))
        self.assertEqual((fb["n_stops"], fb["stops_all_complete"]), (2, False))
        self.assertEqual((fb["n_late503_findings"], fb["late503_findings_only_in_run"]), (2, True))
        self.assertIn("nor the complete 60 s twins answered one. On this host", gap)
        self.assertIn("R-budget-0001 and R-budget-0002 answered 503 deadline, 0.9 s to 2.5 s past their budgets (R-budget-0001 "
                      "0.9 s on direct-local and 0.9 s on direct-remote; R-budget-0002 2.5 s on direct-local (a finding) and "
                      "2.5 s on direct-remote (a finding)). Two of these 503s came more than max(1 s, the budget) after the "
                      "budget and are findings, the chapter's and the run's only findings;", gap)
        # no twin is said to have answered a 503: the twins are named only where the gap says they answered no stop
        sentence_503 = gap.split("answered one.", 1)[1].split("So the prefix property", 1)[0]
        self.assertNotIn("twin", sentence_503)
        self.assertNotIn("they answered", gap)
        self.assertNotIn("completed within", gap)
        self.assertIn("two of the six tests did not pass on every access path (R-budget-0001, R-budget-0002): the coverage gap on the "
                      "chapter page reads their 503 deadlines past the budget and their detail blocks follow the table; the four "
                      "others passed on every access path exercised: R-twin-0002", tex)

    def test_chapter_r_budget_no_503_other_outcomes(self):
        """VERIFY-R5-3 issue 1: without a 503 the coverage gap says 'The stops completed within their budgets' only when
        every stop did; a stop that errored or answered a 200 past its budget is neither, so the gap names the one that
        completed. 'Of these 503s' counts the 503 deadlines that are findings, so a 200 past its budget (a finding
        too) is not called a 503, and 'the run's only findings' is claimed only when the run's findings are exactly
        the 503 ones."""
        tex, gap, fb = self._budget_fixture("budget-stop-error", "stop_error")
        self.assertEqual((fb["n_stops"], fb["n_stops_complete"], fb["stops_all_complete"], fb["n_late503"]), (2, 1, False, 0))
        self.assertEqual((fb["n_finding_executions"], fb["n_late503_findings"], fb["late503_findings_only_in_run"]), (0, 0, False))
        self.assertIn("nor the complete 60 s twins answered one. R-budget-0001 completed within its budget without a stop. "
                      "So the prefix property is inconclusive wherever it was asked", gap)
        for text in ("The stops completed within their budgets", "of these 503s", "answered 503 deadline"):
            self.assertNotIn(text, gap, text)
        tex, gap, fb = self._budget_fixture("budget-over-budget", "over_budget")
        self.assertEqual((fb["n_stops"], fb["n_stops_complete"], fb["stops_all_complete"], fb["n_late503"]), (2, 1, False, 0))
        self.assertEqual((fb["n_finding_executions"], fb["n_late503_findings"], fb["late503_findings_only_in_run"]), (2, 0, False))
        self.assertIn("nor the complete 60 s twins answered one. R-budget-0001 completed within its budget without a stop. "
                      "So the prefix property is inconclusive wherever it was asked", gap)
        for text in ("The stops completed within their budgets", "of these 503s", "answered 503 deadline", "only finding"):
            self.assertNotIn(text, gap, text)
        tex, gap, fb = self._budget_fixture("budget-late-and-over", "late_and_over")
        self.assertEqual((fb["n_stops"], fb["n_stops_complete"], fb["stops_all_complete"], fb["n_late503"]), (2, 0, False, 2))
        self.assertEqual((fb["n_finding_executions"], fb["n_late503_findings"], fb["late503_findings_only_in_run"]), (4, 2, False))
        self.assertIn("R-budget-0001 answered 503 deadline, 1.3 s past its budget (R-budget-0001 1.3 s on direct-local (a finding) "
                      "and 1.3 s on direct-remote (a finding)). Two of these 503s came more than max(1 s, the budget) after the "
                      "budget and are findings; SPEC section 4.5 names the pieces", gap)
        for text in ("Four of these 503s", "only finding", "completed within"):
            self.assertNotIn(text, gap, text)
        # the complete fixture: every stop completed, the plural sentence (unchanged)
        tex, gap, fb = self._budget_fixture("budget-complete-2", "complete")
        self.assertEqual((fb["n_stops"], fb["n_stops_complete"], fb["stops_all_complete"]), (2, 2, True))
        self.assertIn("The stops completed within their budgets without a stop, so the property was not exercised.", gap)


# ---------------------------------------------------------------------------------------------------------------
# the full document (2026-10-15): the chapters T..D written, the facts their prose reads, two more fixture pairs

def mismatch_fixture(run_dir, perturbed):
    """The fake run with its chapter M read as M.json's M-planted and M-protein families (M-001 the planted window,
    M-002 the substituted residue), in two versions: as built (both pass on direct-local), or perturbed (both
    answer a 200 on direct-local whose list the deadline withheld: `time_limited` names the pattern, the property
    check inconclusive; P-001's second direct-local repeat is time-limited too but passes, so it is not counted)."""
    from test_compare_report import line
    manifest = build_run(run_dir)
    if perturbed:
        arms_dir = os.path.join(run_dir, "arms")
        path = os.path.join(arms_dir, "direct-local", "executions.jsonl")
        rows, _ = sc.read_jsonl(path)
        os.remove(path)
        withheld = [{"check": "status", "outcome": "pass", "detail": "ok"},
                    {"check": "contains_source_context", "outcome": "inconclusive",
                     "detail": "0 of 3 source contexts listed (withheld deadline): the list is not complete"}]
        for r in rows:
            if r["scenario"] in ("M-001", "M-002"):
                r = dict(r, outcome="inconclusive", wall_ms=60000.0, time_limited="pattern 0 time-limited", checks=withheld)
            elif r["scenario"] == "P-001" and r["repeat"] == 1:
                r = dict(r, time_limited="pattern 0 time-limited")
            sc.append_jsonl(path, r)
    merge_run.merge(run_dir)
    sc.write_json_atomic(os.path.join(run_dir, "compare.json"), cmp.compare(run_dir, (CASES, {"cases": len(CASES)})))
    mapping = json.loads(json.dumps(FAKE_MAPPING))
    del mapping["_families"]["M-one"]
    del mapping["_units"]["M-one"]
    mapping["_families"]["M-planted"] = {"chapter": "M", "title": "planted substitutions", "units": ["M-planted"]}
    mapping["_families"]["M-protein"] = {"chapter": "M", "title": "protein mismatches", "units": ["M-protein-level1"]}
    mapping["_units"]["M-planted"] = {"family": "M-planted", "unit_type": "list", "unit_title": "planted",
                                      "property": "p", "members": ["M-001"]}
    mapping["_units"]["M-protein-level1"] = {"family": "M-protein", "unit_type": "list", "unit_title": "residues",
                                             "property": "p", "members": ["M-002"]}
    mapping["M-001"] = {"chapter": "M", "family": "M-planted", "unit_id": "M-planted", "unit_type": "list", "variant": "slot"}
    mapping["M-002"] = {"chapter": "M", "family": "M-protein", "unit_id": "M-protein-level1", "unit_type": "list",
                        "variant": "residue"}
    return mapping


def walk_fixture(run_dir, perturbed):
    """The fake run with its chapter R read as T.json's T-walk family (R-001 the walk), in two versions: as built
    (R-001 fails on direct-local, passes on direct-remote), or perturbed (every arm answers a 200 whose seed the
    time budget stopped before the record's minimum: `walk_in_record` inconclusive on each, so no arm showed the
    property)."""
    manifest = build_run(run_dir)
    if perturbed:
        arms_dir = os.path.join(run_dir, "arms")
        cut = [{"check": "status", "outcome": "pass", "detail": "ok"},
               {"check": "walk_in_record", "outcome": "inconclusive",
                "detail": "right: the record followed for 656 of 1000 bp, stopped by time_budget: a stated stop the "
                          "budget allows"}]
        for arm in ("direct-local", "direct-remote", "mcp"):
            path = os.path.join(arms_dir, arm, "executions.jsonl")
            rows, _ = sc.read_jsonl(path)
            os.remove(path)
            for r in rows:
                if r["scenario"] == "R-001":
                    r = dict(r, outcome="inconclusive", http=200, time_limited="seed 0 stopped by its time budget",
                             checks=cut, digest="walk")
                sc.append_jsonl(path, r)
    merge_run.merge(run_dir)
    sc.write_json_atomic(os.path.join(run_dir, "compare.json"), cmp.compare(run_dir, (CASES, {"cases": len(CASES)})))
    mapping = json.loads(json.dumps(FAKE_MAPPING))
    del mapping["_families"]["R-one"]
    del mapping["_units"]["R-one"]
    mapping["_families"]["T-walk"] = {"chapter": "R", "title": "walks", "units": ["T-walk"]}
    mapping["_units"]["T-walk"] = {"family": "T-walk", "unit_type": "pairs", "unit_title": "the walk", "property": "p",
                                   "members": ["R-001"]}
    mapping["R-001"] = {"chapter": "R", "family": "T-walk", "unit_id": "T-walk", "unit_type": "pairs",
                        "variant": "cje, 1,000 bp"}
    return mapping


class FullDocumentTest(unittest.TestCase):
    """The chapters T..D written (2026-10-15): every chapter's prose rendered on the fake run prints no figure of
    its own and leaves nothing unresolved but the scenarios, families, units and relations the fake run lacks; the
    conclusions of chapters M and T change with the data on two fixtures each, through the facts their prose reads
    (the time limitations that left a 200 inconclusive, in M's and T's readings, and the unit root); the facts."""

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
        cls.prose = dict((ch, sc.read_json(os.path.join(PROSE_DIR, "%s.json" % ch))) for ch in WRITTEN_CHAPTERS)
        cls.dirs, cls.mappings = {}, {}
        for name, builder, perturbed in (("m-plain", mismatch_fixture, False), ("m-perturbed", mismatch_fixture, True),
                                         ("t-plain", walk_fixture, False), ("t-perturbed", walk_fixture, True)):
            d = os.path.join(cls.tmp, name)
            os.makedirs(d)
            cls.mappings[name] = builder(d, perturbed)
            cls.dirs[name] = d

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def ctx(self, prose, run=None, units=None):
        import sections
        run = run or self.runobj
        units = units or self.units
        ctx = sections.Context(run, units, prose, {}, {}, list(units.chapters))
        for ch in ctx.chapters:
            ctx.footnotes[ch] = cells.Footnotes()
        return ctx

    @staticmethod
    def parts(text):
        conclusion = text.split("\\parhead{Conclusion}", 1)[1].split("\\begingroup", 1)[0]
        gap = text.split("\\parhead{Coverage and remaining limits}", 1)[1].split("\\parhead{Background}", 1)[0]
        return conclusion, gap

    def test_every_chapter_renders_from_computed_facts(self):
        """Each written chapter's page and every family's prose, rendered on the fake run through its chapter P
        (another run entirely): every flag resolves except those naming a scenario, family, unit or relation the
        fake run lacks, no run figure of the v4 sample survives on the chapter page, and the fake run's own counts
        are printed (its four-member family prints 'four tests', which is the point)."""
        import sections
        for ch in WRITTEN_CHAPTERS:
            data = self.prose[ch]
            for fid, entry in data["families"].items():
                prose = {"P": dict(data)}
                prose["P"]["families"] = {"P-pos": dict(entry)}
                ctx = self.ctx(prose)
                sections.family_section(ctx, "P", "P-pos")  # the fake run's own four-member family prints its counts
                roots = set(n.strip("[]").lstrip("!").split("[")[0] for n in ctx.unresolved)
                self.assertTrue(roots <= {"x", "fam", "rel", "unit"}, "%s %s: unresolved %s" % (ch, fid, sorted(set(ctx.unresolved))))
            prose = {"P": dict(data)}
            prose["P"]["families"] = {}
            ctx = self.ctx(prose)
            whole = sections.chapter(ctx, "P")
            conclusion, gap = self.parts(whole)
            for lit in V4_LITERALS:
                self.assertNotIn(lit, conclusion + gap, "%s: the v4 literal %r survives" % (ch, lit))
            for name in re.findall(r"\(UNRESOLVED CONDITION ([^)]*)\)", conclusion + gap):
                self.assertTrue(name.startswith(("x[", "fam[", "rel[", "unit[")), "%s: unresolved %s" % (ch, name))
            roots = set(n.strip("[]").lstrip("!").split("[")[0] for n in ctx.unresolved)
            self.assertTrue(roots <= {"x", "fam", "rel", "unit"}, "%s: unresolved %s" % (ch, sorted(set(ctx.unresolved))))
            self.assertIn("5 of the 5 tests pass on direct-local", conclusion.replace("\\", ""), ch)
            # the gate consequence sits in the coverage gap, or in the conclusion where the chapter leads with its
            # evidence gaps and the gate consequence (M and N after READ-301, READ-312)
            self.assertIn("the deploy gate was not judged", (conclusion + gap).lower(), ch)

    def render(self, name, ch, data, families):
        """(ctx, the conclusion, the coverage gap, {family: its prose paragraph}) of the fake run's chapter `ch`
        rendered with a written chapter's prose and the named families' entries."""
        import sections
        runobj = reader.Run(self.dirs[name])
        units = units_mod.resolve(runobj.manifest, self.mappings[name])
        prose = {ch: dict(data)}
        prose[ch]["families"] = dict((fid, dict(data["families"][fid])) for fid in families)
        ctx = self.ctx(prose, runobj, units)
        text = sections.chapter(ctx, ch)
        conclusion, gap = self.parts(text)
        fams = {}
        for fid in families:
            fams[fid] = text.split("fam:%s" % fid, 1)[1].split("\\unithead", 1)[0]
        return ctx, conclusion, gap, fams

    def test_mismatch_fixtures_change_conclusions(self):
        """Chapter M on two fixtures: the withheld lists the deadline left (a 200 with a time limitation naming the
        pattern, the property inconclusive) are counted in the conclusion, the coverage gap (the family's and the
        unit's counts, through `fam[M-planted]` and `unit[M-protein-level1]`) and the family prose; a time-limited
        200 that passed is not counted."""
        M = self.prose["M"]
        ctx_a, conc_a, gap_a, fam_a = self.render("m-plain", "M", M, ["M-planted", "M-protein"])
        ctx_b, conc_b, gap_b, fam_b = self.render("m-perturbed", "M", M, ["M-planted", "M-protein"])
        self.assertNotIn("whose list the deadline withheld", conc_a)
        self.assertIn("every planted window completed on both direct access paths", conc_a)
        self.assertIn("In this run every M-planted test was answered within its budget on the direct access paths.", gap_a)
        self.assertIn("Every substituted residue was found within its budget on the direct access paths.", gap_a)
        # the conclusion counts the exceptions and points at the exception index; the ids are the index's and the
        # detail blocks' (READ-301, READ-312)
        self.assertIn("Two tests read inconclusive on direct-local and zero on direct-remote, and two answered a 200 "
                      "whose list the deadline withheld on some access path; the exception index on this page groups "
                      "them by what the exception means", conc_b)
        self.assertNotIn("(M-001, M-002)", conc_b)
        self.assertIn("the planted searches have completion gaps where the budget bound (1 of the 1 planted windows read "
                      "inconclusive on direct-local and 0 on direct-remote", conc_b)
        self.assertIn("In this run 1 of the 1 M-planted tests read inconclusive on direct-local and 0 on direct-remote: "
                      "one answered a 200 whose list the deadline withheld on some access path (1 on direct-local, 0 on "
                      "direct-remote, 0 via MCP); the exception index below names each test", gap_b)
        self.assertIn("Of the 1 substituted residues 1 read inconclusive on direct-local and 0 on direct-remote; no "
                      "substituted residue was found at level 1 on direct-local, so the residue claim rests on no "
                      "evidence from this host", gap_b)
        self.assertIn("In this run every test passed on every access path exercised.", fam_a["M-protein"])
        self.assertIn("In this run one of the 1 tests did not pass on every access path (M-002); the detail blocks follow "
                      "the tables.", fam_b["M-protein"])
        self.assertIn("1 inconclusive on direct-local, 0 on direct-remote; the detail blocks follow the table.", fam_b["M-planted"])
        for ctx in (ctx_a, ctx_b):
            roots = set(n.strip("[]").lstrip("!").split("[")[0] for n in ctx.unresolved)
            self.assertTrue(roots <= {"x"}, sorted(set(ctx.unresolved)))
        v = ctx_b.chapter_values("M")
        self.assertEqual((v["n_time_limited_200"], v["n_time_limited_200_tests"], v["time_limited_200_ids"]), (2, 2, "M-001, M-002"))
        self.assertEqual(v["time_limited_200_text"], "M-001 \ue00160 s\ue002 on direct-local; M-002 \ue00160 s\ue002 on direct-local")
        self.assertEqual(v["time_limited_text"], "M-001 on direct-local (pattern 0 time-limited); M-002 on direct-local (pattern 0 time-limited)")
        self.assertEqual(v["n_time_cut_tests"], 0)
        self.assertEqual(v["unit"]["M-protein-level1"]["n_inconclusive_dl"], 1)
        self.assertEqual(v["fam"]["M-planted"]["n_time_limited_200_tests"], 1)
        self.assertFalse(v["direct_all_pass"])
        self.assertEqual(v["manifest_version"], 1)
        self.assertEqual(v["x"]["M-001"]["dl"]["time_limited"], "pattern 0 time-limited")
        self.assertEqual(v["x"]["M-001"]["dl"]["checks"]["contains_source_context"]["outcome"], "inconclusive")
        self.assertEqual(v["x"]["M-001"]["dr"]["checks"], {})
        # a time-limited 200 that passed (P-001 r1) is not a limitation that left the claim untested
        p = ctx_b.chapter_values("P")
        self.assertEqual((p["n_time_limited"], p["n_time_limited_200"]), (0, 0))
        self.assertTrue(ctx_a.chapter_values("M")["direct_all_pass"])

    def test_walk_fixtures_change_conclusions(self):
        """Chapter T on two fixtures: a walk every arm's time budget stopped before the record's minimum is a time
        cut the coverage gap counts, named among the tests no arm showed the property on; T-walk's and T-longwalk's
        prose state it; as built (a failure, no cut) the other branches print."""
        import sections
        T = self.prose["T"]
        ctx_a, conc_a, gap_a, fam_a = self.render("t-plain", "R", T, ["T-walk"])
        ctx_b, conc_b, gap_b, fam_b = self.render("t-perturbed", "R", T, ["T-walk"])
        self.assertIn("In this run no walk was cut by its time budget or the attempt's bound before it showed its "
                      "property.", gap_a)
        self.assertIn("In this run 1 test had such an execution (R-001). On one of them (R-001) no access path and no repeat "
                      "showed the property, so the walk at that radius rests on no evidence from this host", gap_b)
        self.assertIn("a direct access path failed, errored or disagreed with the other", conc_a)
        self.assertIn("no failure, no error and no disagreement", conc_b)
        self.assertIn("In this run 1 test carries a detail block after the table (an access path did not pass, or the "
                      "paths answered differently: R-001).", fam_a["T-walk"])
        self.assertIn("(an access path did not pass, or the paths answered differently: R-001); 1 of the family's tests "
                      "was cut by a time stop before the minimum (R-001).", fam_b["T-walk"])
        # the long walks' prose reads the cut's evidence
        prose = {"R": dict(T)}
        prose["R"]["families"] = {"T-walk": dict(T["families"]["T-longwalk"])}
        runobj = reader.Run(self.dirs["t-perturbed"])
        ctx = self.ctx(prose, runobj, units_mod.resolve(runobj.manifest, self.mappings["t-perturbed"]))
        text = sections.family_section(ctx, "R", "T-walk").split("\\unithead", 1)[0]
        cut = "right: the record followed for 656 of 1000~bp, stopped by time\\_budget: a stated stop the budget allows"
        self.assertIn("one of the one walks was stopped by a stated budget before the record's minimum (R-001 on "
                      "direct-local r0: %s, on direct-remote r0: %s and on mcp r0: %s); none passed on both direct "
                      "access paths." % (cut, cut, cut), text)
        self.assertNotIn("the walk preset's work cap", text)   # no resource stop in this fixture
        self.assertIn("since whether one 60,000~ms call completes it depends on the load. Both are direct-path tests", text)
        v_t = ctx.chapter_values("R")
        self.assertEqual((v_t["n_cap_stop_tests"], v_t["cap_stop_ids"], v_t["n_budget_stop_tests"], v_t["budget_stop_ids"]),
                         (0, "", 1, "R-001"))
        self.assertEqual(v_t["budget_stop_text"], v_t["time_cut_text"])
        v = ctx_b.chapter_values("R")
        self.assertEqual((v["n_time_cut_tests"], v["time_cut_ids"], v["n_time_cut_untested_tests"], v["time_cut_untested_ids"]),
                         (1, "R-001", 1, "R-001"))
        self.assertEqual(v["n_time_limited"], 3)
        self.assertEqual(v["n_time_limited_200"], 0)
        self.assertEqual(v["x"]["R-001"]["dl"]["time_limited"], "seed 0 stopped by its time budget")
        # the first arm's clause with its name, the explanation after the second ': ' dropped (VERIFY-T M6)
        self.assertEqual(v["x"]["R-001"]["dl"]["checks"]["walk_in_record"]["evidence"],
                         "right: the record followed for 656 of 1000 bp, stopped by time_budget")
        self.assertEqual(v["x"]["R-001"]["dl"]["checks"]["walk_in_record"]["detail"],
                         "right: the record followed for 656 of 1000 bp, stopped by time_budget: a stated stop the budget allows")
        self.assertEqual(v["fam"]["T-walk"]["n_time_cut_tests"], 1)
        self.assertEqual(v["unit"]["T-walk"]["n_time_cut_tests"], 1)
        a = ctx_a.chapter_values("R")
        self.assertEqual((a["n_time_cut_tests"], a["n_time_limited"]), (0, 0))
        for ctx in (ctx_a, ctx_b):
            roots = set(n.strip("[]").lstrip("!").split("[")[0] for n in ctx.unresolved)
            self.assertTrue(roots <= {"x", "fam"}, sorted(set(ctx.unresolved)))

    def test_catalogue_covers_every_check_kind(self):
        """Every check kind the v4 and the v5 manifest expect has a catalogue entry with its texts, and one entry
        serves a kind two chapters use (the routes of `invariants` list /resolve, /pattern and /traverse)."""
        catalogue = (sc.read_json(os.path.join(PROSE_DIR, "checks.json")) or {}).get("checks") or {}
        kinds = set(e.get("check") for m in (v4()[0], v5()[0]) for s in m["scenarios"] for e in s.get("expect") or [])
        self.assertEqual(sorted(kinds - set(catalogue)), [])
        self.assertEqual(sorted(catalogue["invariants"]["routes"]), ["POST /pattern", "POST /resolve", "POST /traverse"])
        self.assertNotIn("this round lists", catalogue["invariants"]["asserts"])
        cases = (sc.read_json(os.path.join(PROSE_DIR, "cases.json")) or {}).get("short") or {}
        named = set()
        for m in (v4()[0], v5()[0]):
            for s in m["scenarios"]:
                bd = (s.get("mcp") or {}).get("by_design")
                for cid in (bd if isinstance(bd, list) else [bd]):
                    if cid:
                        named.add(cid)
        self.assertEqual(sorted(named - set(cases)), [])


class _Stub(object):
    """A bare object with the attributes given (a Run or a Context stand-in for the unit tests below)."""

    def __init__(self, **kw):
        for k, v in kw.items():
            setattr(self, k, v)


class VerifiedFullDocumentTest(unittest.TestCase):
    """The fixes after the Opus verifications of the full v4 document (VERIFY-T, VERIFY-PLXM, VERIFY-QCND,
    2026-10-15): the /traverse budget read from strategy.bounds; a by-design arm's repeat outside the case shown
    with its recorded outcome; why a pair is unverified, from the comparison's fields; a 2xx the service answered
    in the host's place; a check that did not apply left out of 'property tested'; an accepted input's property
    expectations in a refusals row; a walk check's clause; the host state on a direct phrase; several claims per
    family, each over its members, with a conditional on the manifest; the coverage line's exercised repeats."""

    def test_budget_ms_reads_every_route(self):
        scenarios = {"t": {"route": "POST /traverse", "request": {"strategy": {"bounds": {"time_budget_ms": 15000}}}},
                     "r": {"route": "POST /resolve", "request": {"bounds": {"time_budget_ms": 2500}}},
                     "p": {"route": "POST /pattern", "request": {"time_budget_ms": 60000}},
                     "n": {"route": "POST /pattern", "request": {}}}
        run = _Stub(scenario=lambda sid: scenarios.get(sid) or {})
        run.budget_ms = lambda sid: reader.Run.budget_ms(run, sid)
        self.assertEqual(reader.Run.budget_ms(run, "t"), 15000)
        self.assertEqual(reader.Run.budget_ms(run, "r"), 2500)
        self.assertEqual(reader.Run.budget_ms(run, "p"), 60000)
        self.assertIsNone(reader.Run.budget_ms(run, "n"))
        self.assertEqual(blocks._budget_text(run, "t"), "15,000 ms")
        self.assertEqual(blocks._budget_text(run, "n"), "none")

    @staticmethod
    def _by_design_arm(other_outcome="pass", other_class="no_direct_reference", direct_verified=True, case_repeat=0):
        """An MCP arm standing by design on r0 (the host's deadline relayed) with r1 outside the case; with
        `case_repeat` 1 the case covers r1 and r0 is the repeat outside it (M-planted-0015 on the v4 run)."""
        st = cells.ArmState("mcp")
        st.present = True
        st.glyph_key = "service_by_design"
        st.glyph_marks = ["u"]
        covered = {"scenario": "N-1", "arm": "mcp", "repeat": 0, "outcome": "service_by_design", "http": 502,
                   "wall_ms": 63000.0, "host_cache": "cold", "service_by_design_case": "pattern.deadline_ends_task",
                   "error": {"code": "search_failed", "text": "pattern_deadline"},
                   "checks": [{"check": "status", "outcome": "fail", "detail": "http 502"},
                              {"check": "source_absent", "outcome": "fail", "detail": "no answer", "property": True}]}
        other = {"scenario": "N-1", "arm": "mcp", "repeat": 1, "outcome": other_outcome, "http": 200, "wall_ms": 22000.0,
                 "host_cache": "warm", "inconclusive_reason": "checks" if other_outcome == "inconclusive" else None,
                 "time_limited": "pattern 0 time-limited" if other_outcome == "inconclusive" else "",
                 "checks": [{"check": "status", "outcome": "pass", "detail": "http 200"},
                            {"check": "source_absent", "outcome": other_outcome, "property": True,
                             "detail": "absent" if other_outcome == "pass" else "list cut (cut time)"}]}
        st.lines = [covered, other]
        if case_repeat == 1:
            covered["repeat"], other["repeat"] = 1, 0
            covered["host_cache"], other["host_cache"] = "warm", "cold"
            st.lines = [other, covered]
        st.agreement = {"class": "service_by_design", "by": "answer", "case": "pattern.deadline_ends_task",
                        "repeats": [case_repeat], "verified": False, "direct_verified": direct_verified,
                        "direct_detail": "repeat(s) ['direct-local r0', 'direct-remote r0'] match" if direct_verified
                        else "no direct execution answered {'http': 503, 'code': 'deadline'} (direct_mismatch timing: unverified)",
                        "mcp_verified": True, "other_repeats": {"class": other_class, "mcp_outcome": other_outcome}}
        return st

    def test_by_design_repeat_outside_the_case(self):
        """VERIFY-PLXM H2, VERIFY-QCND H2: the case's repeat reads by design, the other repeat its recorded outcome
        with its property cell; the coverage line counts the arm as exercised on that repeat."""
        st = self._by_design_arm("pass")
        self.assertTrue(cells.by_design_line(st, st.lines[0]))
        self.assertFalse(cells.by_design_line(st, st.lines[1]))
        self.assertEqual([l["repeat"] for l in cells.exercised_lines(st)], [1])
        run = _Stub(scenario=lambda sid: {"route": "POST /pattern", "request": {"time_budget_ms": 60000}},
                    marked_aside=lambda sid, line: {}, cases={}, budget_ms=lambda sid: 60000)
        rows = blocks._execution_rows(st, run, "N-1", "mcp", lambda n: n == "source_absent")
        self.assertEqual(rows[0][1], "\\gbd\\ by design")
        self.assertEqual(rows[0][5], "–")
        self.assertEqual(rows[1][1], "\\gpass\\ pass")
        self.assertEqual(rows[1][5], "yes (1 of 1)")
        self.assertEqual(rows[1][4], "60,000 ms")
        st2 = self._by_design_arm("inconclusive")
        rows = blocks._execution_rows(st2, run, "N-1", "mcp", lambda n: n == "source_absent")
        self.assertEqual(rows[1][1], "\\ginc\\ inconclusive (by the checks)")
        self.assertEqual(rows[1][5], "no")
        # the one phrase of the arm names the other repeat and why the pair is unverified
        run.cases = {"pattern.deadline_ends_task": {"service_answer": {"outcome": "relayed"}}}
        text = blocks._not_exercised_text(st, run, "N-1", {"pattern.deadline_ends_task": "the deadline relayed"})
        self.assertIn("repeat r1 outside the case: r1 recorded pass (the row below)", text)
        self.assertIn("UNVERIFIED pair: r1 outside the case (a pass) has no direct answer of the same kind to compare with", text)
        # the coverage line: exercised on one repeat, passed there
        states = {"N-1": {"mcp": st}, "N-2": {"mcp": st2}}
        self.assertEqual(cells.coverage_text(states, "mcp"), "exercised 2/2 (2 on one repeat only, the other by design), passed 1")
        # a not-exercised arm's refusal rows stay not exercised, whatever the runner wrote
        ni = cells.ArmState("mcp")
        ni.present, ni.glyph_key = True, "not_implemented_mcp"
        ni.lines = [{"repeat": 0, "outcome": "service_by_design", "http": 400}]
        self.assertTrue(cells.by_design_line(ni, ni.lines[0]))
        self.assertEqual(cells.coverage_text({"R-1": {"mcp": ni}}, "mcp"), "exercised 0/1; \\attn{not exercised 1\\,\\gni{}}")

    def test_class_from_the_cases_execution(self):
        """VERIFY-FINAL R1: where the case covers r1 only (r0 a 200 outside the case), the class word and the
        status of the arm's phrase, the arm facts and the edge cell all come from r1, the relayed 502: never
        'answered by the host' (r0's class) beside 'http 502' (r1's status)."""
        st = self._by_design_arm("inconclusive", case_repeat=1)
        self.assertEqual([l["repeat"] for l in st.lines], [0, 1])
        self.assertEqual([l["repeat"] for l in blocks.class_lines(st)], [1])
        scenario = {"route": "POST /pattern", "request": {"time_budget_ms": 60000},
                    "expect": [{"check": "status", "args": {"http": 200}}]}
        cases = {"pattern.deadline_ends_task": {"service_answer": {"outcome": "relayed"}}}
        run = _Stub(scenario=lambda sid: scenario, cases=cases, budget_ms=lambda sid: 60000,
                    is_property=lambda name: name == "source_absent", marked_aside=lambda sid, line: {})
        # the first answered execution alone (r0, the 200) would read answered: the regression
        self.assertEqual(blocks.mcp_class(st.lines, scenario, st.agreement, cases)[0], "answered")
        self.assertEqual(blocks.mcp_class(blocks.class_lines(st), scenario, st.agreement, cases)[0], "relayed")
        short = (sc.read_json(os.path.join(PROSE_DIR, "cases.json")) or {}).get("short") or {}
        text = blocks._not_exercised_text(st, run, "N-1", short)
        self.assertIn("host refusal relayed through the service (http 502 \\texttt{search\\_\\allowbreak{}failed})", text)
        self.assertNotIn("answered by the host", text)
        self.assertIn("repeat r0 outside the case: r0 recorded inconclusive (the row below)", text)
        # the short reason is the one clause; the verification rule is appendix B's (R2)
        self.assertIn("comes back as 502", text)
        self.assertNotIn("the pair is verified where", text)
        # the arm facts (the class word the prose may print) read the same execution
        f = facts.arm_facts(_Stub(run=run, state=lambda sid: {"mcp": st}), "N-1", "mcp")
        self.assertEqual((f["class"], f["class_word"]), ("relayed", "host refusal relayed through the service"))
        # and so do its status and code: r1's 502, not r0's 200 (VERIFY-FIX2 L3)
        self.assertEqual((f["http"], f["code"]), ("502", "search_failed"))
        self.assertTrue(f["phrase"].startswith("by design (http 502 search_failed), the pair unverified"), f["phrase"])
        self.assertIn("r0 outside the case inconclusive", f["phrase"])
        # the case on r0 (the 19 other blocks of the v4 run) reads as before
        st0 = self._by_design_arm("inconclusive")
        self.assertEqual([l["repeat"] for l in blocks.class_lines(st0)], [0])
        self.assertIn("host refusal relayed through the service (http 502", blocks._not_exercised_text(st0, run, "N-1", short))
        # an arm the service did not decide reads every execution
        plain = cells.ArmState("mcp")
        plain.present, plain.glyph_key, plain.lines = True, "pass", [{"repeat": 0, "http": 200}, {"repeat": 1, "http": 200}]
        self.assertEqual(len(blocks.class_lines(plain)), 2)

    def test_short_reasons_are_short_sentences(self):
        """VERIFY-FIX2 L1: every short local reason of cases.json is one short sentence (the owner's review, point 4:
        a short local reason, the long form in appendix B): at most 50 words and at most two clauses (one ';', the
        second where the host refuses too, with its own text), never the 77 words R2 removed. The reasons the
        owner's chapter R reviews passed run to about 45 words."""
        short = (sc.read_json(os.path.join(PROSE_DIR, "cases.json")) or {}).get("short") or {}
        self.assertGreater(len(short), 50)
        for case, text in sorted(short.items()):
            self.assertLessEqual(len(text.split()), 50, "%s: %d words" % (case, len(text.split())))
            self.assertLessEqual(text.count(";"), 1, case)
            self.assertFalse(text.rstrip().endswith("."), case)

    def test_short_reason_one_clause_and_the_verification_rule_in_b(self):
        """VERIFY-FINAL R2: the deadline case's short reason is one clause (the owner's rule for cases.json), and
        its verification rule is printed in appendix B under the short reason, nowhere else."""
        cases = sc.read_json(os.path.join(PROSE_DIR, "cases.json")) or {}
        short = cases["short"]["pattern.deadline_ends_task"]
        self.assertNotIn(";", short)
        self.assertLessEqual(len(short.split()), 25, short)
        rule = cases["verification"]["pattern.deadline_ends_task"]
        self.assertIn("verified where a direct access path answered the 503 too", rule)
        self.assertIn("unverified where no direct access path answered the 503", rule)
        self.assertIn("no direct 200 to compare with", rule)
        self.assertEqual(sorted(set(cases["verification"]) - set(cases["short"])), [])
        table = {"pattern.deadline_ends_task": {"areas": ["M"], "agreement": "by_design", "scenario_match": "a 503",
                                                "service_answer": {"outcome": "relayed", "http": 502, "code": "search_failed"},
                                                "direct_expectation": {"http": 503}},
                 "pattern.caps_zero": {"areas": ["M"], "agreement": "by_design", "scenario_match": "a cap of 0",
                                       "service_answer": {"outcome": "refused_at_submit"}}}
        ctx = _Stub(run=_Stub(order=[], table=table, cases={}), units=_Stub(by_scenario={}, unit_of=lambda s: None),
                    chapters=["M"], selected_scenarios=lambda: [], short_reasons=lambda: cases["short"],
                    verification_rules=lambda: cases["verification"])
        text = sections.by_design_table(ctx)
        self.assertIn("in short: the host's /pattern 503 deadline ends the service's task at once", text)
        self.assertEqual(text.count("verification: "), 1)
        self.assertLess(text.index("in short: the host's /pattern 503"), text.index("verification: the pair is verified where"))
        self.assertIn("where the case has one, its verification rule", text)
        # a case without a rule prints its short reason alone
        self.assertIn("in short: the service's request model takes every cap", text)
        self.assertEqual(sections.Context.verification_rules(_Stub(cases_prose={})), {})

    def test_unverified_reason_from_the_comparison(self):
        """VERIFY-PLXM H3, VERIFY-QCND H3: the reason is the comparison's, not a fixed sentence."""
        st = self._by_design_arm("pass")
        self.assertEqual(blocks.unverified_reason(st.agreement, st),
                         "r1 outside the case (a pass) has no direct answer of the same kind to compare with")
        self.assertEqual(blocks.unverified_cause(st.agreement), "other_repeat")
        st = self._by_design_arm("pass", other_class="agree", direct_verified=False)
        self.assertEqual(blocks.unverified_reason(st.agreement, st),
                         "the direct side did not match the case: no direct execution answered http 503 deadline "
                         "(a timing difference, not a disagreement)")
        self.assertEqual(blocks.unverified_cause(st.agreement), "direct")
        st = self._by_design_arm("inconclusive", direct_verified=None)
        self.assertEqual(blocks.unverified_cause(st.agreement), "direct+other_repeat")
        self.assertIn("; r1 outside the case (an inconclusive execution) has no direct answer", blocks.unverified_reason(st.agreement, st))
        self.assertIsNone(blocks.unverified_reason({"class": "service_by_design", "verified": True}))
        self.assertIsNone(blocks.unverified_cause({"class": "agree"}))
        # the classes line and the facts read the same reason
        self.assertIn("UNVERIFIED: r1 outside the case", blocks._classes_text(
            _Stub(comp_entry=lambda sid: {"mcp_agreement": self._by_design_arm("pass").agreement}, depth={}),
            "N-1", {"mcp": self._by_design_arm("pass")})[0])

    def test_service_answered_in_the_hosts_place(self):
        """VERIFY-QCND H5: a 2xx the service answered from its own data (the host not reached) is no refusal."""
        cases = {"capabilities.traverse_capabilities_document": {"service_answer": {"outcome": "answered"}},
                 "case.refused": {"service_answer": {"outcome": "refused_at_submit", "code": "pattern_invalid"}}}
        answered = [{"http": 200, "outcome": "inconclusive", "inconclusive_reason": "checks", "host_cache": "not_reached",
                     "service_by_design_case": "capabilities.traverse_capabilities_document"}]
        cls, words = blocks.mcp_class(answered, {}, {"case": "capabilities.traverse_capabilities_document"}, cases)
        self.assertEqual(cls, "service_answered")
        self.assertIn("from the service's own data, the host not reached", words)
        self.assertEqual(blocks.MCP_CLASS_WORD[cls], "answered by the service in the host's place")
        refused = [{"http": 200, "outcome": "service_by_design", "host_cache": "not_reached",
                    "error": {"code": "pattern_invalid", "text": "x"}}]
        self.assertEqual(blocks.mcp_class(refused, {}, {"case": "case.refused"}, cases)[0], "refused_before_host")
        self.assertIn("answered by the service in the host's place", sections.EDGE_LEGEND)
        self.assertIn("answered by the host through the service", sections.EDGE_LEGEND)

    def test_check_that_did_not_apply(self):
        """VERIFY-T M9: answer_fields with when_http on another status tested nothing."""
        line = {"checks": [{"check": "status", "outcome": "pass", "detail": "http 200"},
                           {"check": "answer_fields", "outcome": "pass", "detail": "fields as expected", "property": True},
                           {"check": "answer_fields", "outcome": "pass", "property": True,
                            "detail": "http 200: the fields are checked on an http 409 only"}]}
        self.assertTrue(cells.check_not_applicable(line["checks"][2]))
        self.assertFalse(cells.check_not_applicable(line["checks"][1]))
        self.assertEqual(cells.property_tested(line), ("yes", 1, 1))
        self.assertEqual(blocks._property_cell(line, None), "yes (1 of 1)")

    def test_accepted_input_in_a_refusals_row(self):
        """VERIFY-T M7: an accepted input's property expectations follow its status."""
        is_prop = lambda n: n in ("answer_fields", "clamped")
        dup = {"expect": [{"check": "status", "args": {"http": 200}}, {"check": "invariants", "args": {}},
                          {"check": "answer_fields", "args": {"equals": {"results.1.duplicate": True}}}]}
        text, _ = blocks.refusal_expected(dup, None, is_prop)
        self.assertEqual(text, "http 200 · \\texttt{results.\\allowbreak{}1.\\allowbreak{}duplicate} = \\texttt{true}")
        clamp = {"expect": [{"check": "status", "args": {"http": 200}},
                            {"check": "clamped", "args": {"field": "bounds.time_budget_ms", "requested": 70000,
                                                         "effective": 60000, "where": "strategy"}}]}
        self.assertTrue(blocks.refusal_expected(clamp, None, is_prop)[0].startswith("http 200 · clamped "))
        refused = {"expect": [{"check": "status", "args": {"http": 400, "text": "no"}},
                              {"check": "answer_fields", "args": {"present": ["usage"]}}]}
        self.assertEqual(blocks.refusal_expected(refused, None, is_prop)[0], "http 400 \\texttt{no}")
        self.assertEqual(blocks.expect_text("answer_fields", {"present": ["usage"], "absent": ["x"], "when_http": 409}),
                         "present \\texttt{usage}; absent \\texttt{x} (on http 409)")

    def test_walk_clause_and_host_state(self):
        """VERIFY-T M6 and VERIFY-PLXM M5: the walk check's first arm with its clause; a direct execution's host
        state beside its time."""
        self.assertEqual(facts._clause("left: the record followed for 793 of 1000 bp, stopped by time_budget; right: 1 route(s)"),
                         "left: the record followed for 793 of 1000 bp, stopped by time_budget")
        self.assertEqual(facts._clause("not testable: seed 0 stopped by its attempt's duration bound (stated by the host: x)"),
                         "not testable: seed 0 stopped by its attempt's duration bound")
        self.assertEqual(facts._clause("seed 0 not walked (outcome 'failed'): nothing to follow"), "seed 0 not walked (outcome 'failed')")
        self.assertEqual(facts._wall_text({"wall_ms": 1900.0, "host_cache": "cold"}, "direct-local"), "1.9 s (host-cold)")
        self.assertEqual(facts._wall_text({"wall_ms": 99.0, "host_cache": "warm"}, "direct-remote"), "99 ms")
        self.assertEqual(facts._wall_text({"wall_ms": 1900.0, "host_cache": "cold"}, "mcp"), "1.9 s")
        run = _Stub(scenario=lambda sid: {}, budget_ms=lambda sid: None)
        self.assertEqual(facts.exec_phrase(run, "x", {"outcome": "pass", "http": 200, "wall_ms": 70.0, "host_cache": "cold"},
                                           "direct-local"), "a pass in 70 ms (host-cold)")

    def test_claim_items(self):
        """VERIFY-PLXM H1 and VERIFY-T H4: a family's claims as a list, each over its members."""
        self.assertEqual(facts.claim_items({"claim": "a", "claim_members": ["M-1"]}), [("a", ["M-1"])])
        self.assertEqual(facts.claim_items({"claim": [{"claim": "a", "members": ["M-1"]}, "b"]}), [("a", ["M-1"]), ("b", None)])
        self.assertEqual(facts.claim_items({}), [])
        self.assertEqual(facts.fill("a[[manifest_version = 5 ?, and b]]", {"manifest_version": 4}), "a")
        self.assertEqual(facts.fill("a[[manifest_version = 5 ?, and b]]", {"manifest_version": 5}), "a, and b")
        # the v6 delta (2026-10-17): `>=` compares numbers, so that manifest v6 takes the v5 branch
        self.assertEqual(facts.fill("a[[manifest_version >= 5 ?, and b | , or c]]", {"manifest_version": 4}), "a, or c")
        self.assertEqual(facts.fill("a[[manifest_version >= 5 ?, and b | , or c]]", {"manifest_version": 5}), "a, and b")
        self.assertEqual(facts.fill("a[[manifest_version >= 5 ?, and b | , or c]]", {"manifest_version": 6}), "a, and b")
        self.assertEqual(facts.fill("a[[manifest_version = 5 ?, and b]]", {"manifest_version": 6}), "a")


def tint_fixture(run_dir):
    """The fake run of test_compare_report with slow executions (the owner's decision of 2026-10-15, the tints):
    P-001 r0 on direct-local at 61 s with a 12 s server time (a series cell), P-004 at 30 s on direct-local (a
    replicates member), P-001's warm MCP call at 14 s with a 12 s backend and a 56 s terminal lag (a cell's be
    line), P-002's failing MCP call with a 12 s backend phase (a detail block's phases line) and R-001's failing
    direct-local execution at 61 s (a detail block's row). Built from small data here, never from a run directory."""
    build_run(run_dir)
    arms_dir = os.path.join(run_dir, "arms")

    def rewrite(arm, edit):
        path = os.path.join(arms_dir, arm, "executions.jsonl")
        rows, _ = sc.read_jsonl(path)
        os.remove(path)
        for r in rows:
            sc.append_jsonl(path, edit(r))

    def local(r):
        if r["scenario"] == "P-001" and r["repeat"] == 0:
            r = dict(r, wall_ms=61000.0, server=dict(r["server"], elapsed_ms=12000.0))
        elif r["scenario"] == "P-004":
            r = dict(r, wall_ms=30000.0, server=dict(r["server"], elapsed_ms=200.0))
        elif r["scenario"] == "R-001":
            r = dict(r, wall_ms=61000.0, server=dict(r["server"], elapsed_ms=30500.0))
        return r

    def mcp(r):
        if r["scenario"] == "P-001" and r["repeat"] == 1:
            r = dict(r, wall_ms=14000.0, mcp_phases=dict(r["mcp_phases"], backend_ms=12000.0, terminal_lag_ms=56000.0))
        elif r["scenario"] == "P-002":
            r = dict(r, mcp_phases={"submit_ms": 50.0, "queue_ms": 200.0, "backend_ms": 12000.0,
                                    "terminal_lag_ms": 100.0, "paging_ms": 30.0})
        return r
    rewrite("direct-local", local)
    rewrite("mcp", mcp)
    merge_run.merge(run_dir)
    sc.write_json_atomic(os.path.join(run_dir, "compare.json"), cmp.compare(run_dir, (CASES, {"cases": len(CASES)})))


class TimeTintTest(unittest.TestCase):
    """The tints of a printed execution time (the owner's decision of 2026-10-15): above 10 s yellow, above 55 s
    red, strictly, by cells' two constants and nothing else; on the cells' walls, server and backend times, the
    detail blocks' rows and phases, the replicates aggregates and the timing tables; never on a difference, a
    bound, a connection set-up or a figure quoted in the prose; the colour tokens in doc.sty; one legend sentence
    added to the cell grammar's lines; a tinted time as wide as an untinted one (the pagination does not move)."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.run_dir = os.path.join(cls.tmp, "run")
        os.makedirs(cls.run_dir)
        tint_fixture(cls.run_dir)
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

    def test_thresholds_and_formatters(self):
        """Strictly above 10 s yellow, strictly above 55 s red; fmt_exec is fmt_t on the tint; fmt_times tints each
        figure alone with the mark outside; a difference (fmt_delta) and the plain fmt_t never carry a tint."""
        self.assertEqual((cells.TIME_TINT_YELLOW_MS, cells.TIME_TINT_RED_MS), (10000.0, 55000.0))
        self.assertIsNone(cells.time_tint(None))
        self.assertIsNone(cells.time_tint(10000))
        # decided by the printed value (the owner's decision of 2026-10-16): 10,400 ms prints 10 s, no tint; 10,500
        # prints 11 s, yellow; 55,400 prints 55 s, yellow; 55,500 prints 56 s, red
        self.assertIsNone(cells.time_tint(10000.5))
        self.assertIsNone(cells.time_tint(10400))
        self.assertEqual(cells.time_tint(10500), "yellow")
        self.assertEqual(cells.time_tint(55000), "yellow")
        self.assertEqual(cells.time_tint(55000.5), "yellow")
        self.assertEqual(cells.time_tint(55400), "yellow")
        self.assertEqual(cells.time_tint(55500), "red")
        self.assertEqual((cells.printed_ms("10 s"), cells.printed_ms("15,593 ms"), cells.printed_ms("3.1 s"),
                          cells.printed_ms("<1 ms"), cells.printed_ms("–")), (10000.0, 15593.0, 3100.0, 0.0, None))
        # a quoted figure is its own printed value: 10452 ms is yellow, 10000 ms is not
        self.assertIn("\ue000", cells.tint_quoted_times("client wall 10452 ms"))
        self.assertNotIn("\ue000", cells.tint_quoted_times("client wall 10000 ms"))
        self.assertEqual(cells.fmt_exec(61000), "\\tred{61 s}")
        self.assertEqual(cells.fmt_exec(12000), "\\tyellow{12 s}")
        self.assertEqual(cells.fmt_exec(3140), "3.1 s")
        self.assertEqual(cells.fmt_exec(None), "–")
        self.assertEqual(cells.fmt_t(61000), "61 s")
        self.assertEqual(cells.fmt_delta(61000), "+61 s")
        self.assertEqual(cells.fmt_delta(-61000), "−61 s")
        self.assertEqual(cells.fmt_times([61000, 900]), "\\tred{61} / 0.9 s")
        self.assertEqual(cells.fmt_times([12000, 30], ["c", None]), "\\tyellow{12 s}\\tmark{c} / 30 ms")
        self.assertEqual(cells.fmt_times([900, 61000], [None, "a"]), "0.9 / \\tred{61 s}\\tmark{a}")
        self.assertEqual(cells.fmt_times([61000, 900], tint=False), "61 / 0.9 s")
        self.assertEqual(cells.fmt_pair(12000, 61000), "\\tyellow{12} / \\tred{61 s}")
        # the closing tables: a signed difference and a connection set-up stay plain, an execution time is tinted
        import sections
        d = {"p50": 12000.0, "p90": 61000.0, "max": 61000.0, "n": 3}
        self.assertEqual(sections._stats_cells(d), ["\\tyellow{12 s}", "\\tred{61 s}", "\\tred{61 s}", "3"])
        self.assertEqual(sections._stats_cells(d, signed=True), ["+12 s", "+61 s", "+61 s", "3"])
        self.assertEqual(sections._stats_cells(d, tint=False), ["12 s", "61 s", "61 s", "3"])
        # the constants are the document's own: no other module of the suite defines them
        for name in ("suite_common.py", "report.py", "compare.py", "checks.py"):
            path = os.path.join(SUITE, name)
            if os.path.exists(path):
                with open(path, encoding="utf-8") as fh:
                    self.assertNotIn("TIME_TINT_", fh.read(), name)

    def test_cells_and_aggregates(self):
        """A cell's wall, server and backend lines and a replicates aggregate on synthetic arm states."""
        st = cells.ArmState("direct-local")
        st.present, st.worst, st.glyph_key = True, "pass", "pass"
        st.lines = [{"wall_ms": 61000.0, "host_cache": "cold", "server": {"elapsed_ms": 12000.0}},
                    {"wall_ms": 900.0, "host_cache": "warm", "server": {"elapsed_ms": 400.0}}]
        l1, l2, l3 = cells.cell_parts(st)
        # the unit once: the first figure bare on its tint, the last with the unit inside the tint
        self.assertEqual(l1, "\\gpass\\ \\tred{61}\\tmark{c} / 0.9 s")
        self.assertEqual(l2, "[\\tyellow{12} / 0.4 s]")
        self.assertIsNone(l3)
        st.lines[1]["wall_ms"] = 14000.0
        self.assertEqual(cells.cell_parts(st)[0], "\\gpass\\ \\tred{61}\\tmark{c} / \\tyellow{14 s}")
        mcp = cells.ArmState("mcp")
        mcp.present, mcp.worst, mcp.glyph_key = True, "pass", "pass"
        mcp.lines = [{"wall_ms": 14000.0, "host_cache": "warm", "server": {"elapsed_ms": 11000.0},
                      "mcp_phases": {"backend_ms": 12000.0}}]
        l1, l2, l3 = cells.cell_parts(mcp)
        self.assertEqual((l1, l2, l3), ("\\gpass\\ \\tyellow{14 s}", "[\\tyellow{11 s}]", "be \\tyellow{12 s}"))
        warm = cells.ArmState("direct-local")
        warm.present, warm.worst, warm.glyph_key = True, "pass", "pass"
        warm.lines = [{"wall_ms": 30000.0, "host_cache": "warm", "server": {"elapsed_ms": 12000.0}}]
        cold = cells.ArmState("direct-local")
        cold.present, cold.worst, cold.glyph_key = True, "pass", "pass"
        cold.lines = [{"wall_ms": 30000.0, "host_cache": "cold", "server": {"elapsed_ms": 12000.0}}]
        text, line2 = cells.aggregate({"a": {"direct-local": warm}, "b": {"direct-local": cold}}, "direct-local")
        self.assertEqual(text, "ok 2/2 · \\tyellow{30} / \\tyellow{30 s}\\tmark{c 1/2}")
        self.assertEqual(line2, "[\\tyellow{12 s}]")
        self.assertEqual(cells.fmt_times([3100, 900]), "3.1 / 0.9 s")

    def test_rendered_document(self):
        """On the fixture: the series cell (wall and server), the mcp cell's be line, the replicates aggregate,
        the family summary's p50 / p90, the detail blocks' rows and phases line, the closing tables; the legend
        sentence once in the reading guide with the thresholds from the constants; a block's size estimate counts
        a tinted figure as itself."""
        self.assertEqual(self.rc, 0)
        p = self.tex("chapter_P.tex")
        self.assertIn("\\gpass\\ \\tred{61 s}\\tmark{c} / 6 ms", p)  # P-001 on direct-local (the units differ)
        self.assertIn("[\\tyellow{12 s} / 3 ms]", p)  # its server times
        self.assertIn("be 0.3 / \\tyellow{12 s}", p)  # P-001's mcp cell: the backend line, r0 then r1
        self.assertIn("\\gpass\\ 0.9 / \\tyellow{14 s}", p)  # its walls
        self.assertIn("ok 1/1 · \\tyellow{30} / \\tyellow{30 s}\\tmark{c}", p)  # the replicates aggregate (P-004)
        self.assertIn("backend \\tyellow{12 s}, terminal lag 0.1 s", p)  # P-002's block: the phases line
        self.assertIn("host-cold wall p50 / p90: direct-local 5 ms / \\tyellow{49 s}", p)  # the chapter's timing line
        self.assertIn("wall p50 / p90 host-cold: \\tyellow{15} / \\tyellow{52~s}", p)  # P-pos's summary line
        self.assertNotIn("\\tyellow{\\tyellow", p)
        self.assertNotIn("\\tred{\\t", p)
        # a budget, a time distance and the 'requested budget' column never carry a tint
        self.assertNotIn("\\tyellow{60,000", p)
        self.assertNotIn("none\\tyellow", p)
        r = self.tex("chapter_R.tex")
        self.assertIn("& \\tred{61 s} (warm) & none &", r)  # R-001's block: the direct-local row, the budget plain
        self.assertIn("\\gfail\\fn{1}\\ \\tred{61 s} [\\tyellow{31 s}]", r)  # R-001's pairs cell: wall, server 30.5 s
        main = self.tex()
        self.assertIn("\\tred{56 s}", main)  # C.3: the terminal lag's maximum
        self.assertIn("\\tred{61 s}", main)  # C.7: chapter P's host-cold wall maximum on direct-local
        legend = main.split("\\bline{times}{", 1)
        self.assertEqual(len(legend), 2)
        sentence = legend[1].split("}\n", 1)[0]
        self.assertIn("above 10 s is set on a \\tyellow{yellow} tint and above 55 s on a \\tred{red} one", sentence)
        self.assertIn("decided by the value as printed (a time printed 10 s carries none, one printed 56 s is red), in "
                      "the tables, the detail blocks and the prose alike, a time quoted from a check's text and the time "
                      "a walk's stop used too; a difference (a delta, an overhead, how late a 503 was), a budget, a cap "
                      "and a distance in time carry no tint", sentence)
        self.assertEqual(main.count("\\bline{times}"), 1)
        # the two reviewed lines of the cell grammar are as they were, the new one after them
        import sections
        g = sections.grammar_lines(None)
        self.assertLess(g.index("\\bline{aggregates}"), g.index("\\bline{edge tables}"))
        # V6-003: the tint and attention rules are appendix F's; I.3 keeps a short colour key
        self.assertLess(g.index("\\bline{edge tables}"), g.index("\\bline{colour}"))
        self.assertNotIn("\\bline{times}", g)
        self.assertIn("\\bline{times}", sections.colour_rules(None))
        self.assertIn("members never exercised read not exercised: n, never ok 0/0; the coverage line under the "
                      "table counts per access path}", g)
        self.assertIn("derived from the evidence, never from prose}", g)
        # the connection set-up table and the signed tables carry no tint
        c4 = main.split("connect\\_ms per arm", 1)[0].rsplit("\\begin{longtable}", 1)[1]
        self.assertNotIn("\\tyellow", c4)
        self.assertNotIn("\\tred", c4)
        # the prose's figures are plain text (escaped afterwards): a scenario fact prints the time with its tint as
        # a marker, never as a macro, and latex.esc sets the macro last
        f = facts.scenario_facts(sections.Context(self.runobj, self.units, {}, {}, {}, list(self.units.chapters)),
                                 "R-001")
        self.assertEqual(f["dl"]["wall"], "\ue00161 s\ue002")
        self.assertEqual(f["dl"]["walls"], "\ue00161 s\ue002")
        self.assertIn(f["dl"]["server"], ("\ue00031 s\ue002", "–"))
        self.assertNotIn("\\t", json.dumps(f["dl"]))
        self.assertEqual(latex.esc("a pass in %s (host-cold)" % f["dl"]["wall"]), "a pass in \\tred{61 s} (host-cold)")
        p1 = facts.scenario_facts(sections.Context(self.runobj, self.units, {}, {}, {}, list(self.units.chapters)),
                                  "P-001")
        self.assertIn("a pass in \ue00161 s\ue002", p1["dl"]["phrase"])  # P-001's phrase, as the prose quotes it
        self.assertIn("a pass in \\tred{61 s}", latex.esc(p1["dl"]["phrase"]))
        self.assertEqual(blocks.estimate_lines(["\\bline{a}{\\tred{61 s} and \\tyellow{12} s}"]),
                         blocks.estimate_lines(["\\bline{a}{61 s and 12 s}"]))
        self.assertEqual(blocks.estimate_lines(["x" * 149 + " \\tred{61 s}"]), blocks.estimate_lines(["x" * 149 + " 61 s"]))

    def test_prose_markers_and_quoted_times(self):
        """A figure quoted in prose carries its tint as a marker (fmt_exec_text) that latex.esc sets as the macro
        after escaping; an execution time quoted inside a check's or the runner's text (`client wall N ms`, `...
        after N ms`) is marked by its value (tint_quoted_times, esc_times) and a lateness, a budget or any other
        figure of the text stays; a marker a cut left open is closed; a close without an open is dropped."""
        self.assertEqual(cells.fmt_exec_text(900), "0.9 s")
        self.assertEqual(cells.fmt_exec_text(None), "–")
        self.assertEqual(cells.fmt_exec_text(12000), "\ue00012 s\ue002")
        self.assertEqual(cells.fmt_exec_text(61000), "\ue00161 s\ue002")
        self.assertEqual(latex.esc("in \ue00012 s\ue002, 5 ms & \ue00161 s\ue002"), "in \\tyellow{12 s}, 5 ms \\& \\tred{61 s}")
        self.assertEqual(latex.esc("in \ue00012 s\ue002", limit=10), "in \\tyellow{12 s}")
        self.assertEqual(latex.esc("in \ue00012 s\ue002 and more", limit=9), "in \\tyellow{12…}")   # cut inside: closed
        self.assertEqual(latex.set_tints("x\ue002y"), "xy")
        self.assertEqual(latex.join_units(latex.esc("a pass in \ue00161 s\ue002")), "a pass in \\tred{61~s}")
        quoted = "503 deadline delivered 7452 ms after the 3000 ms budget (client wall 10452 ms)"
        self.assertEqual(cells.tint_quoted_times(quoted),
                         "503 deadline delivered 7452 ms after the 3000 ms budget (client wall \ue00010452 ms\ue002)")
        self.assertEqual(cells.esc_times(quoted),
                         "503 deadline delivered 7452 ms after the 3000 ms budget (client wall \\tyellow{10452 ms})")
        self.assertEqual(cells.esc_times("client wall 9000 ms"), "client wall 9000 ms")
        self.assertEqual(cells.esc_times("resolve_timeout after 30395 ms (the settled attempt's settled_at - dispatched_at)"),
                         "resolve\\_timeout after \\tyellow{30395 ms} (the settled attempt's settled\\_at - dispatched\\_at)")
        self.assertEqual(cells.esc_times("traversal_resolve_deadline after 61218 ms"),
                         "traversal\\_resolve\\_deadline after \\tred{61218 ms}")
        self.assertEqual(cells.esc_times("delivered 25496 ms after the 20000 ms budget"),
                         "delivered 25496 ms after the 20000 ms budget")
        self.assertEqual(facts._cut("x (client wall 45496 ms) y", 400), "x (client wall \ue00045496 ms\ue002) y")
        self.assertEqual(facts._clause("503 deadline delivered 25496 ms after the 20000 ms budget (client wall 45496 ms): "
                                       "the deadline is read where the work polls it"),
                         "503 deadline delivered 25496 ms after the 20000 ms budget (client wall \ue00045496 ms\ue002)")

    def test_resource_stop_phrase(self):
        """A walk a stated resource stopped: the resource, its knob and its figures from the answer's resource_stop
        (blocks.resource_stops, reader.Run.raw_answer), '' without a raw record or a stop."""
        import gzip
        run_dir = os.path.join(self.tmp, "stop-run")
        os.makedirs(os.path.join(run_dir, "arms", "direct-local", "raw"))
        raw = "arms/direct-local/raw/T-001.direct-local.0.json.gz"
        body = {"results": [{"resource_stop": {"resource": "work", "requested": 100000000, "effective": 100000000,
                                               "used": 100000275, "phase": "traversal"},
                             "arms": {"left": {"cap_trigger": {"at_bp": 3618, "reason": "resource_limit"}, "complete_to_bp": 3618},
                                      "right": {"complete_to_bp": 3617}}},
                            {"outcome": {}}]}
        with gzip.open(os.path.join(run_dir, raw), "wt", encoding="utf-8") as fh:
            json.dump({"answer": {"http": 200, "body": body}}, fh)

        class FakeRun(object):
            dir = run_dir
            raw_answer = reader.Run.raw_answer
            policy_path = None
            policy_file_sha256 = reader.Run.policy_file_sha256

            def budget_ms(self, sid):
                return None
        run = FakeRun()
        line = {"raw": raw, "checks": [{"check": "walk_in_record", "outcome": "inconclusive",
                                        "detail": "left: the record followed for 3618 of 10000 bp, stopped by resource_limit"}]}
        self.assertEqual(blocks.resource_stops(run, line), [
            {"seed": 0, "resource": "work", "knob": "bounds.max_work_units", "requested": 100000000, "used": 100000275,
             "at": "3,618 bp (left), 3,617 bp (right)"}])
        self.assertEqual(blocks.resource_stop_phrase(run, line),
                         "seed 0: the work budget (bounds.max_work_units, 100,000,000 units; used 100,000,275 units) at "
                         "3,618 bp (left), 3,617 bp (right)")
        self.assertEqual(blocks.resource_stop_phrase(run, {"raw": "arms/none.json.gz"}), "")
        # a resource outside the three knobs (the MCP runner's attempt_deadline) is named as the answer names it
        raw2 = "arms/mcp/raw/T-001.mcp.0.json.gz"
        os.makedirs(os.path.join(run_dir, "arms", "mcp", "raw"))
        with gzip.open(os.path.join(run_dir, raw2), "wt", encoding="utf-8") as fh:
            json.dump({"answer": {"http": 200, "body": {"results": [{"resource_stop": {
                "resource": "attempt_deadline", "requested": 25000, "used": 22085}}]}}}, fh)
        # its figures are milliseconds: the time used is an execution time, with its unit and on its tint by the
        # printed value (VERIFY-V6-2 M1); the requested bound stays plain
        self.assertEqual(blocks.resource_stop_phrase(run, {"raw": raw2}),
                         "seed 0: a stop by attempt_deadline (requested 25,000 ms; used \ue00022,085 ms\ue002)")
        self.assertEqual(latex.esc(blocks.resource_stop_phrase(run, {"raw": raw2})),
                         "seed 0: a stop by attempt\\_deadline (requested 25,000 ms; used \\tyellow{22,085 ms})")
        self.assertEqual(blocks.resource_stop_phrase(run, {}), "")
        self.assertIsNone(run.policy_file_sha256())
        text = blocks.adjudication(dict(line, outcome="inconclusive", http=200), run, "T-001", False)
        self.assertIn("inconclusive: walk\\_in\\_record (left: the record followed for 3618 of 10000 bp, stopped by "
                      "resource\\_limit) — the stop, as the answer states it: seed 0: the work budget "
                      "(bounds.max\\_work\\_units, 100,000,000 units; used 100,000,275 units) at 3,618 bp (left), 3,617 bp "
                      "(right)",
                      text)

    def test_sty_tokens_and_width(self):
        """doc.sty defines the two colour tokens and the macros; with xelatex, a tinted time sets exactly as wide
        as an untinted one and the PDF's text is the figure alone."""
        with open(os.path.join(DOC, "templates", "doc.sty"), encoding="utf-8") as fh:
            sty = fh.read()
        self.assertIn("\\definecolor{timeyellow}{HTML}{FFE866}", sty)
        self.assertIn("\\definecolor{timered}{HTML}{FF8A80}", sty)
        self.assertIn("\\newcommand\\tyellow[1]{\\timetint{timeyellow}{#1}}", sty)
        self.assertIn("\\newcommand\\tred[1]{\\timetint{timered}{#1}}", sty)
        self.assertEqual(cells.TINT_MACRO, {"yellow": "\\tyellow", "red": "\\tred"})
        xelatex = make_doc.find_xelatex()
        if not xelatex:
            self.skipTest("xelatex not installed")
        d = os.path.join(self.tmp, "sty")
        os.makedirs(d)
        shutil.copyfile(os.path.join(DOC, "templates", "doc.sty"), os.path.join(d, "doc.sty"))
        with open(os.path.join(d, "t.tex"), "w", encoding="utf-8") as fh:
            fh.write("\\documentclass[10pt]{article}\\usepackage{doc}\\newcommand\\docrunid{tint}\n"
                     "\\newlength\\wplain\\newlength\\wtint\\begin{document}\n"
                     "\\settowidth\\wplain{[61 s / 0.9 s] be 12 s}\\settowidth\\wtint{[\\tred{61 s} / 0.9 s] be \\tyellow{12 s}}\n"
                     "\\typeout{TINTWIDTHS \\the\\wplain\\space\\the\\wtint}\n"
                     "\\tred{61 s} \\tyellow{12 s} 3.1 s\n\\end{document}\n")
        ok, message = make_doc.build_pdf(d, os.path.join(d, "t.tex"), xelatex)
        self.assertTrue(ok, message)
        with open(os.path.join(d, "t.log"), encoding="utf-8", errors="replace") as fh:
            log = fh.read()
        m = re.search(r"TINTWIDTHS (\S+) (\S+)", log)
        self.assertIsNotNone(m, "no width line in the log")
        self.assertEqual(m.group(1), m.group(2))
        pdftotext = shutil.which("pdftotext")
        if pdftotext:
            import subprocess
            text = subprocess.run([pdftotext, os.path.join(d, "t.pdf"), "-"], stdout=subprocess.PIPE).stdout.decode()
            self.assertIn("61 s 12 s 3.1 s", text)


# ---------------------------------------------------------------------------------------------------------------
# the readability review of the full v4 report (READ-001 .. READ-025 and the shared templates, 2026-10-15)

class ReadabilityReviewTest(unittest.TestCase):
    """The readability review (REVIEW-INDEX.md, 2026-10-15): one vocabulary rendered in the legend, appendix F and the
    glossary (READ-022); the gate's problems by what clears them, the decision sentence and the scope read from the
    policy (READ-002, READ-023); the counting units (READ-006); the untested table and the exception index
    (READ-005, READ-111); the detail block's interpretation and keyed artifacts (READ-114); the budget-scope table
    and the ruling's decision table (READ-014, READ-015); the Timing names and the paired ratios (READ-009 to
    READ-011); appendix B's evidence states and grouped notes (READ-017, READ-018); the exemption kinds (READ-025);
    the chapter purpose slot, the unit heading's timing key and the summary's denominators (READ-101, READ-102,
    READ-105); the access-path wording (READ-008)."""

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
        cls.overview = sc.read_json(os.path.join(PROSE_DIR, "overview.json"))
        cls.glossary = sc.read_json(os.path.join(PROSE_DIR, "glossary.json"))
        cls.closing = sc.read_json(os.path.join(PROSE_DIR, "closing.json"))
        cls.checks = sc.read_json(os.path.join(PROSE_DIR, "checks.json"))
        cls.budget_dirs, cls.budget_mappings = {}, {}
        for name, perturbed in (("plain", False), ("perturbed", True)):
            d = os.path.join(cls.tmp, "budget-" + name)
            os.makedirs(d)
            cls.budget_mappings[name] = budget_fixture(d, perturbed)
            cls.budget_dirs[name] = d

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def ctx(self, prose=None, overview=None, glossary=None, closing=None, run=None, units=None):
        import sections
        run = run or self.runobj
        units = units or self.units
        ctx = sections.Context(run, units, prose or {}, overview or {}, glossary or {}, list(units.chapters),
                               closing=closing)
        for ch in ctx.chapters:
            ctx.footnotes[ch] = cells.Footnotes()
        return ctx

    def test_one_vocabulary(self):
        """READ-022: the execution outcomes of the vocabulary are suite_common's, the check items its CHECK_OUTCOMES;
        the reading guide's legend, appendix F's tables and the glossary's core entries print the same texts."""
        import sections
        import vocabulary
        self.assertEqual(set(e["key"] for e in vocabulary.EXECUTION_OUTCOMES), set(sc.OUTCOMES))
        self.assertEqual(set(k for k, _ in vocabulary.CHECK_OUTCOMES), set(sc.CHECK_OUTCOMES))
        rows = vocabulary.legend_rows()
        self.assertEqual([w for _, w, _ in rows][-1], "no execution")
        self.assertEqual(len(rows), len(sc.OUTCOMES) + 1)
        legend = sections.legend_table(self.ctx(overview=self.overview))
        tables = sections.classes_tables()
        for e in vocabulary.EXECUTION_OUTCOMES:
            self.assertIn(latex.inline(e["legend"]), legend, e["key"])
            self.assertIn(latex.tt(e["key"]), tables)
            self.assertIn(latex.inline(e["meaning"]), tables, e["key"])
        self.assertIn("check item & meaning", tables)
        self.assertIn("execution outcome & decided & meaning & consequence", tables)
        self.assertIn(latex.inline(vocabulary.STANDING), tables)
        self.assertNotIn("legend", self.overview["sections"]["reading"])  # the glyph legend is the vocabulary's
        # the error outcome names the raised-check difference between the layers, as the code has it
        err = vocabulary.execution_outcome("error")["meaning"]
        self.assertIn("the evaluation itself (checks.evaluate) raised on a direct access path", err)
        self.assertIn("via MCP the same line reads fail", err)
        self.assertIn("a single check that raises inside the evaluation is a failing item, so that execution reads fail "
                      "on every access path", err)
        self.assertIn("the evaluation itself (checks.evaluate) raised", vocabulary.execution_outcome("error")["legend"])
        self.assertNotIn("a check that raised on a direct access path", vocabulary.execution_outcome("error")["legend"])
        self.assertIn("a check that raised is recorded as a failing item", vocabulary.CHECK_OUTCOMES[1][1])
        # the finding names the clock the check reads, never 'without failing the request' (READ-003)
        for text in (vocabulary.execution_outcome("finding")["legend"], vocabulary.execution_outcome("finding")["meaning"],
                     dict(vocabulary.CHECK_OUTCOMES)["finding"]):
            self.assertIn("on the clock the check reads", text)
            self.assertIn("client's clock", text)
            self.assertNotIn("without failing the request", text)
            self.assertNotIn("while the request was answered", text)
        # the standing text says what each count counts: the chapter's outcomes row counts tests (READ-102)
        self.assertIn("a chapter's outcomes row counts tests too", vocabulary.STANDING)
        self.assertIn("I.4's outcome table counts executions", vocabulary.STANDING)
        self.assertIn("mixed where the repeats differ otherwise", vocabulary.STANDING)
        # a joined consequence starts its sentence with a capital (READ-020)
        texts = vocabulary.glossary_texts()
        self.assertTrue(texts["error"].endswith(". Counted with the failures."), texts["error"])
        self.assertIn(". The verdict reads FINDINGS", texts["finding"])
        self.assertIn(". Never a pass and never a failure", texts["inconclusive"])
        # the glossary reads the same texts through the vocab placeholders and renders them
        terms = dict((e["term"], e["definition"]) for e in self.glossary["entries"])
        for term, key in (("execution outcome", "execution_outcome"), ("check outcome", "check_outcome"),
                          ("error", "error"), ("finding", "finding"), ("standing (the glyph)", "standing"),
                          ("access path (arm)", "access_path")):
            self.assertIn("{vocab[%s]}" % key, terms[term], term)
        ctx = self.ctx(glossary=self.glossary)
        text = sections.glossary(ctx)
        self.assertNotIn("UNRESOLVED", text)
        self.assertNotIn("{vocab", text)
        self.assertIn(latex.inline(vocabulary.glossary_texts()["error"]), text)
        self.assertIn("\\textbf{traversal arm (direction)}", text)
        self.assertIn("\\textbf{access path (arm)}", text)
        self.assertEqual(ctx.placeholders()["vocab"], vocabulary.glossary_texts())

    def test_gate_groups_decision_and_scope(self):
        """READ-002, READ-023: a problem's group follows its kind and, for an outcome or a relation, its detail; the
        decision sentence, the groups text, the scope read from the policy (never its file name), the action table,
        the full record beside the document, and the ungated run's statements."""
        import gate_policy
        import sections
        g = sections.gate_group_of
        self.assertEqual(g({"kind": "identity", "detail": "x"}), "evidence")
        self.assertEqual(g({"kind": "configuration", "detail": "x"}), "evidence")
        self.assertEqual(g({"kind": "missing", "detail": "planned, not recorded"}), "evidence")
        self.assertEqual(g({"kind": "agreement", "detail": "x"}), "evidence")
        self.assertEqual(g({"kind": "outcome", "detail": "inconclusive (checks)"}), "unproven")
        self.assertEqual(g({"kind": "outcome", "detail": "by-design pair not verified (x)"}), "unproven")
        self.assertEqual(g({"kind": "outcome", "detail": "finding"}), "violation")
        self.assertEqual(g({"kind": "outcome", "detail": "fail"}), "violation")
        self.assertEqual(g({"kind": "outcome", "detail": "service_by_design by case x, which the scenario does not name"}), "violation")
        self.assertEqual(g({"kind": "property", "detail": "required property x inconclusive"}), "unproven")
        self.assertEqual(g({"kind": "nothing_passed", "detail": "x"}), "unproven")
        self.assertEqual(g({"kind": "relation", "detail": "REL-x on direct-local: inconclusive: no two profiles"}), "unproven")
        self.assertEqual(g({"kind": "relation", "detail": "REL-x on direct-local: fail: differs"}), "violation")
        self.assertEqual(g({"kind": "disagreement", "detail": "x"}), "violation")
        problems = [{"kind": "identity", "detail": "a"}, {"kind": "configuration", "arm": "mcp", "detail": "b"},
                    {"kind": "outcome", "scenario": "P-001", "arm": "mcp", "repeat": 0, "detail": "inconclusive (checks)"},
                    {"kind": "outcome", "scenario": "P-001", "arm": "mcp", "repeat": 1, "detail": "inconclusive (checks)"},
                    {"kind": "outcome", "scenario": "R-001", "arm": "direct-local", "repeat": 0, "detail": "finding"},
                    {"kind": "relation", "arm": "direct-local", "detail": "REL-x on direct-local: inconclusive: y"}]
        stub = _Stub(gate={"problems": problems, "passed": False}, policy=None, policy_path=None, order=list(self.runobj.order),
                     gate_offline=True, reevaluated=False, executions_file="executions.jsonl")
        self.assertEqual(sections.gate_group_counts(stub), {"evidence": 2, "unproven": 3, "violation": 1})
        self.assertEqual(sections.gate_groups_text(stub),
                         "2 release evidence problems, 3 required behaviours not shown, 1 observed violation")
        self.assertEqual(sections.gate_scope(stub)["text"], "no policy pins this run's manifest")
        stub.policy = {"selection": {"routes": ["POST /resolve", "POST /traverse"], "scenarios": ["P-001", "R-001"]},
                       "arms": ["direct-local", "mcp"]}
        scope = sections.gate_scope(stub)
        self.assertEqual(scope["routes"], "POST /resolve, POST /traverse")
        self.assertEqual(scope["arms"], "direct-local, mcp")
        self.assertEqual(scope["selected"], 2)
        self.assertIn("the scenarios of POST /resolve and POST /traverse (2 of the run's %d) on the access paths "
                      "direct-local and mcp" % len(self.runobj.order), scope["text"])
        decision = sections.gate_decision_text(stub)
        self.assertTrue(decision.startswith("This run does not authorise a deploy: the gate failed on 6 problems over the "
                                            "scenarios of POST /resolve"), decision)
        self.assertIn("The count is not a count of service defects: 2 release evidence problems", decision)
        self.assertIn("judged offline by gate.py when this document was built", sections.gate_read_text(stub))
        stub.gate = {"problems": [], "passed": True}
        # V6-001 (the external review of the v6 document): with a policy, the passed decision names the required
        # population after the policy's automatic exclusions
        self.assertTrue(sections.gate_decision_text(stub).startswith(
            "The gate passed: every planned execution was recorded, every execution the policy requires passed"))
        self.assertIn("not the whole service", sections.gate_decision_text(stub))
        # the ungated fake run: no decision judged, and the cover says so
        ctx = self.ctx(overview=self.overview)
        values = ctx.placeholders()
        self.assertTrue(values["gate_decision"].startswith("No deploy decision was judged for this run"))
        self.assertEqual(values["gate_groups_text"], "none")
        self.assertEqual(values["n_gate_problems"], 0)
        self.assertEqual(values["gate_read_text"], "no gate verdict")
        self.assertIn("\\textbf{\\large Deploy gate: not evaluated}", sections.gate_line(ctx))
        # the fake run gated by a policy built from its own manifest: the cover line carries the decision and the scope
        policy = gate_policy.build_policy(self.runobj.manifest, table=(CASES, {"cases": 5}), manifest_file="m.json",
                                          arms=("direct-local", "mcp"), warm_every=0)
        path = os.path.join(self.tmp, "policy.json")
        with open(path, "w") as fh:
            json.dump(policy, fh)
        run = reader.Run(self.run_dir, policy=path)
        self.assertTrue(run.gate_offline)
        ctx = self.ctx(overview=self.overview, run=run)
        line = sections.gate_line(ctx)
        self.assertIn("the gate covers the scenarios of", line)
        self.assertIn("access paths direct-local and mcp", line)
        if run.gate.get("passed"):
            self.assertIn("nothing blocks", line)
        else:
            self.assertIn("the run does not authorise a deploy", line)
            self.assertIn("by what would clear them", line)
            self.assertIn("the count is not a count of service defects", line)
        text = sections.gate_problems(ctx)
        self.assertIn("What the gate read:", text)
        self.assertIn(sections.GATE_PROBLEMS_FILE, ctx.data_files)
        record = json.loads(ctx.data_files[sections.GATE_PROBLEMS_FILE])
        self.assertEqual(len(record["problems"]), len(run.gate["problems"]))
        for p in record["problems"]:
            self.assertIn(p["group"], ("evidence", "unproven", "violation"))
        if run.gate["problems"]:
            self.assertIn("the gate's problems by what they ask of the owner", text)
            self.assertNotIn("more problem texts", text)  # every text is printed (READ-023)
        self.assertNotIn("more problem texts", sections.gate_action_table(ctx, problems))
        action = sections.gate_action_table(ctx, problems)
        self.assertIn("release evidence missing or mismatched", action)
        self.assertIn("1 problem of the access path mcp", action)
        self.assertIn("1 problem of the run's records", action)
        self.assertEqual(action.count("1 scenario-arm pair: "), 2)  # P-001 mcp (two repeats) under unproven, R-001 under violation
        self.assertIn("1 problem of the access path direct-local", action)  # the relation, without a scenario
        self.assertIn("observed violation of a rule the suite applies", action)

    def test_counting_units_and_examples(self):
        """READ-006, READ-007, READ-102, READ-019, READ-021: the counting units table with the run's figures and the
        worked example; the standing examples; the relation examples; the service revision statement."""
        import sections
        ctx = self.ctx(overview=self.overview)
        table = sections.counting_units_table(ctx)
        for unit in ("scenario (test)", "execution", "property claim", "scenario-arm pair", "gate problem"):
            self.assertIn(unit + " &", table)
        self.assertIn("& %d\\\\" % len(self.runobj.order), table)
        self.assertIn("no gate judged", table)
        values = ctx.placeholders()
        self.assertEqual(values["n_property_claims"], sum(d["pairs"] for d in ctx.untested()["by_chapter"].values()))
        example = values["counting_example"]
        self.assertRegex(example, r"^[A-Z]-\d{3} is one scenario; it was executed \d+ times \(")
        self.assertIn("property claim", example)
        standing = values["standing_example"]
        self.assertTrue(standing.startswith("P-004: via MCP the service refused the call at submit, in the host's place"), standing)
        self.assertIn("nothing of the claim was tested through the tools", standing)
        self.assertNotIn("decided the call in the host's place", standing)
        # the words follow the case's service_answer.outcome (READ-007: R-edge-0002's 502 is the host's 400 relayed)
        relayed = _Stub(cases={"c.relay": {"service_answer": {"outcome": "relayed"}}}, table={})
        st = cells.ArmState("mcp")
        st.agreement = {"case": "c.relay"}
        self.assertEqual(sections.service_decision_words(relayed, st), "the service relayed the host's refusal under its own code")
        st.agreement = {"case": "nowhere"}
        self.assertEqual(sections.service_decision_words(relayed, st), "the service decided the call in the host's place")
        answered = _Stub(cases={"c": {"service_answer": {"outcome": "answered"}}}, table={})
        st.agreement = {"case": "c"}
        self.assertEqual(sections.service_decision_words(answered, st), "the service answered in the host's place")
        self.assertIn(values["standing_vs_outcome_example"], ("none in this run",) if not any(
            self.runobj.arm_lines(sid, "mcp") and cells.scenario_states(self.runobj, sid)["mcp"].glyph_key == "not_implemented_mcp"
            and any(l.get("outcome") == "service_by_design" for l in self.runobj.arm_lines(sid, "mcp"))
            for sid in self.runobj.order) else (values["standing_vs_outcome_example"],))
        reading = sections.summary_reading_text(ctx, "P")
        self.assertIn("standing and outcomes both count tests, one per access path", reading)
        self.assertNotIn("count executions", reading)
        self.assertNotIn("counts twice", reading)
        # 'none' once, and only when no test's MCP standing differs from its outcome (V2-M1)
        self.assertIn("the outcomes row the recorded outcome over the repeats (none in this chapter)", reading)
        self.assertNotIn("none in this chapter)", reading.replace("(none in this chapter)", ""))
        self.assertIn("the full reading is I.3's", reading)
        self.assertIn("Agreement counts one comparison per test (its three levels, identical answers, the same check "
                      "outcomes, every property check compared on both: I.3)", reading)
        self.assertEqual(sections.standing_vs_outcome_example_text(ctx, ["P-001"]), "none in this chapter")
        # the outcomes row of every chapter counts tests: per access path it sums to the chapter's tests with an
        # execution there (compare's per_area, one entry per test), and is labelled so (READ-102)
        per_area = self.runobj.summary()["per_area"]
        for ch in self.units.chapters:
            sids = self.units.chapter_scenarios(ch)
            text = sections.chapter_summary(ctx, ch, sids, 1, 1)
            self.assertIn("\\bline{outcomes}{tests ·", text)
            self.assertIn("each test's outcome over its repeats, as the comparison records it", text)
            for arm in cells.ARMS:
                o = (per_area.get(ch) or {}).get(arm) or {}
                if not o:
                    continue
                n_tests = sum(1 for sid in sids if self.runobj.arm_lines(sid, arm))
                self.assertEqual(sum(v for k, v in o.items() if k != "total"), n_tests, (ch, arm, o))
                self.assertEqual(o["total"], n_tests, (ch, arm))
        # the declared exclusions in both units (READ-006)
        self.assertEqual(values["declared_exclusions_clause"], "no exclusion")
        dec = _Stub(declared={"_file": "d.json", "declared_version": 1,
                              "exclusions": [{"scenario": "A", "arms": ["mcp"]}, {"scenario": "B", "arms": ["mcp", "direct-local"]},
                                             {"scenario": "C", "arms": ["direct-local"], "repeats": [1]}],
                              "not_required": [], "relations_excluded": []})
        counts = sections.declared_counts(dec)
        self.assertEqual((counts["exclusions"], counts["exclusion_pairs"], counts["exclusion_repeat_only_pairs"]), (3, 4, 1))
        self.assertEqual(counts["exclusions_clause"], "3 exclusions (4 scenario-arm pairs, 1 of them for repeat 1 only)")
        self.assertIn("{declared_exclusions_clause}", " ".join(self.overview["sections"]["provenance"]["exemptions"]))
        self.assertNotIn("{n_declared_exclusions} scenario-arm pairs", " ".join(self.overview["sections"]["provenance"]["exemptions"]))
        # the findings by what the check read, and the verdict's reason in the document's words (READ-003)
        self.assertEqual(values["finding_kinds_text"], "none")
        self.assertEqual(values["verdict_text"], self.runobj.verdict)  # no FINDINGS verdict here: report.py's line
        perturbed = reader.Run(self.budget_dirs["perturbed"])
        self.assertEqual(sections.finding_kinds(perturbed), {"late_503": 1, "over_200": 0, "walk_200": 0, "other": 0})
        self.assertEqual(sections.finding_kinds_text(perturbed),
                         "1 late 503 deadline (no results carried; beyond the suite's tolerance on the client's clock)")
        # the perturbed fixture's own verdict is INCOMPLETE (an interrupted access path): report.py's line as it is
        self.assertEqual(sections.doc_verdict_text(perturbed), perturbed.verdict)
        findings_run = _Stub(verdict="FINDINGS (1 execution(s) where the host broke a stated bound without failing the "
                                     "request; no failure, no disagreement); 2 inconclusive execution(s) (not counted as pass)",
                             order=perturbed.order, arm_lines=perturbed.arm_lines)
        doc_verdict = sections.doc_verdict_text(findings_run)
        self.assertEqual(doc_verdict, "FINDINGS (1 execution(s) with a finding, a stated bound broken on the clock the "
                                      "check reads: 1 late 503 deadline (no results carried; beyond the suite's tolerance "
                                      "on the client's clock); no failure, no disagreement); 2 inconclusive execution(s) "
                                      "(not counted as pass)")
        self.assertNotIn("without failing the request", doc_verdict)
        pctx = self.ctx(run=perturbed, units=units_mod.resolve(perturbed.manifest, self.budget_mappings["perturbed"]))
        self.assertIn(latex.esc(perturbed.verdict.split(" ", 1)[1], 400), sections.verdict_line(pctx))
        self.assertEqual(pctx.placeholders()["verdict_text"], perturbed.verdict)
        self.assertEqual(pctx.placeholders()["finding_kinds_text"], sections.finding_kinds_text(perturbed))
        self.assertNotIn("without failing the request", json.dumps(self.overview) + json.dumps(self.glossary) + json.dumps(self.closing))
        self.assertEqual(sections.relation_example_text(self.runobj, "subset"), "none in this run")
        rel_run = _Stub(relations=[{"id": "REL-subset-0001", "type": "subset", "a": "A", "b": "B", "title": "every A among B."}])
        self.assertEqual(sections.relation_example_text(rel_run, "subset"), "REL-subset-0001 (A ⊆ B: every A among B)")
        self.assertTrue(values["service_version_line"].startswith("the service's deployed revision is not recorded by this run"))
        self.assertIn("cannot be attributed to a service revision", values["service_version_line"])
        with_version = _Stub(mcp_run={"service_commit": "abc123"})
        self.assertEqual(sections.service_version_text(with_version),
                         "the service's deployed revision as the MCP runner recorded it (read from the service before "
                         "each MCP invocation): service commit abc123")
        with_legacy = _Stub(mcp_run={"service_version": "abc123"})   # a key of an earlier runner's
        self.assertEqual(sections.service_version_text(with_legacy),
                         "the service's deployed revision as the MCP runner recorded it: abc123")
        env = sections.environment_table(ctx)
        self.assertIn("the client that asked, not the service's revision", env)
        self.assertIn("is not recorded by this run", env)

    def test_untested_table_exception_index_and_interpretation(self):
        """READ-005, READ-111, READ-114: the untested properties as a table with the unit's question, the chapter's
        exceptions grouped by what they mean, and the detail block led by the interpretation with the artifacts
        keyed once."""
        import sections
        # the fake run passes every property somewhere: no table
        self.assertEqual(sections.untested_table(self.ctx()), "")
        # the perturbed budget fixture: R-002's 503s leave a property untested; the table names the unit's question
        run = reader.Run(self.budget_dirs["perturbed"])
        units = units_mod.resolve(run.manifest, self.budget_mappings["perturbed"])
        prose = {"R": {"units": {"R-budget": {"property": "the long property line", "question": "Does a stop answer the prefix?"}}}}
        ctx = self.ctx(prose=prose, run=run, units=units)
        text = sections.chapter(ctx, "R")
        table = sections.untested_table(ctx, "R")
        if table:
            self.assertIn("the question (unit) & what this run established & unproven: the check", table)
            self.assertIn("Does a stop answer the prefix?", table)
            self.assertIn("passed it on some access path", table)
            self.assertIn("no access path showed", table)
        self.assertIn("\\bline{standing}{tests ·", text)
        self.assertIn("\\bline{outcomes}{tests ·", text)
        self.assertIn("\\bline{agreement}{test pairs ·", text)
        self.assertIn("\\bline{reading}{standing and outcomes both count tests", text)
        index = sections.exception_index(ctx, "R")
        self.assertIn("by what each exception means", index)
        self.assertIn("\\textbf{a measured operational finding}", index)
        self.assertIn("the gate was not judged", index)
        self.assertIn("\\hyperref[detail:R-002]", index)
        self.assertEqual(sections.exception_index(self.ctx(), "S"), "")  # no footnoted test, no index
        # the groups of a test from its arm states
        st = cells.scenario_states(self.runobj, "R-001")  # fails on direct-local, the arms disagree, mcp inconclusive
        self.assertEqual(sections.exception_groups_of(st), ["failure", "unverified"])
        self.assertEqual(sections.exception_groups_of(cells.scenario_states(self.runobj, "P-001")), [])
        st2 = cells.scenario_states(run, "R-002")
        self.assertIn("finding", sections.exception_groups_of(st2))
        self.assertIn("unverified", sections.exception_groups_of(st2))
        # a finding execution alone: under the property not shown too when its property checks passed nothing
        # (READ-111: a late 503 carried no answer), with the run's is_property; not without it, not when one passed
        only = dict((arm, cells.ArmState(arm)) for arm in cells.ARMS)
        only["direct-local"].present = True
        only["direct-local"].worst = "finding"
        only["direct-local"].lines = [{"outcome": "finding", "checks": [{"check": "budget_respected", "outcome": "finding"},
                                                                      {"check": "resolve_runs", "outcome": "inconclusive"}]}]
        is_prop = lambda name: name == "resolve_runs"
        self.assertEqual(sections.exception_groups_of(only), ["finding"])
        self.assertEqual(sections.exception_groups_of(only, is_prop), ["finding", "unverified"])
        only["direct-local"].lines[0]["checks"][1]["outcome"] = "pass"
        self.assertEqual(sections.exception_groups_of(only, is_prop), ["finding"])
        finding_gap = dict((k, g) for k, _, g in sections.EXCEPTION_GROUPS)["finding"]
        self.assertIn("which carried no answer", finding_gap)
        self.assertNotIn("while the request was answered", finding_gap)
        # V6-302 (the external review of the v6 document): the entry names each test's recorded reason, never a menu
        self.assertIn("for the recorded reason named after each test",
                      dict((k, g) for k, _, g in sections.EXCEPTION_GROUPS)["unverified"])
        # the detail block: requested, observed, unverified first, then the technical record; the artifacts keyed once
        block = blocks.detail_block(self.runobj, "R-001", st, 1, self.units, property_text="the claim in words")
        self.assertIn("\\bline{requested}{the claim in words}", block)
        self.assertIn("\\bline{observed}{direct-local: FAIL (", block)
        self.assertIn("DISAGREES with the reference access path", block)
        self.assertIn("\\bline{unverified}{", block)
        self.assertLess(block.find("\\bline{unverified}"), block.find("the technical record"))
        self.assertIn("\n\\par\\noindent{\\color{black!60}\\small\\itshape the technical record", block)
        self.assertLess(block.find("the technical record"), block.find("\\bline{expected}"))
        self.assertIn("path / repeat & final outcome", block)
        self.assertIn("arms/\\allowbreak{}<access path>/\\allowbreak{}raw/\\allowbreak{}R-\\allowbreak{}001.\\allowbreak{}<access path>.\\allowbreak{}<repeat>.\\allowbreak{}json.\\allowbreak{}gz", block)
        self.assertIn("for dl r0, dr r0, mcp r0", block)
        self.assertIn("\\texttt{executions.\\allowbreak{}jsonl} lines", block)
        # a refusal scenario: the status is the claim
        st_p2 = cells.scenario_states(self.runobj, "P-002")
        block2 = blocks.detail_block(self.runobj, "P-002", st_p2, 2, self.units)
        self.assertIn("\\bline{requested}{", block2)
        self.assertIn("\\bline{unverified}{", block2)
        self.assertIn("\\bline{observed}{direct-local: the expected refusal (http 400)", block2)

    def test_budget_scope_and_decision_table(self):
        """READ-014, READ-015: the route-by-route budget-scope table after the methods, and the ruling's decision
        table in the catalogue's budget_respected entry."""
        import sections
        rows = self.overview["sections"]["methods"]["budget_scope"]
        self.assertEqual([r["route"] for r in rows], ["POST /resolve", "POST /pattern", "POST /traverse"])
        for r in rows:
            for key in ("named_by", "clock", "scope", "not_covered", "answers", "verified"):
                self.assertTrue(r.get(key), (r["route"], key))
        self.assertIn("per seed", rows[2]["named_by"])
        self.assertIn("the other seeds of the request", rows[2]["not_covered"])
        ctx = self.ctx(overview=self.overview)
        text = sections.methods_section(ctx)
        self.assertIn("the scope of a time budget per route", text)
        # V6-003: one full-width block per route, its six statements labelled (no seven-column table)
        self.assertIn("\\textbf{POST /resolve.} \\emph{The budget:}", text)
        for label in ("Its clock starts", "It covers", "It does not cover", "Possible answers", "What the suite verified"):
            self.assertEqual(text.count("\\emph{%s:}" % label), len(rows), label)
        self.assertLess(text.find("\\emph{Times and load.}"), text.find("the scope of a time budget per route"))
        self.assertEqual(sections.budget_scope_table(self.ctx()), "")
        entry = self.checks["checks"]["budget_respected"]
        table = entry["decision_table"]
        self.assertEqual(table["columns"], ["the answer", "the clock read", "the threshold", "the outcome"])
        outcomes = " ".join(r[3] for r in table["rows"])
        for word in ("pass when within it; finding when past it", "fail: the one failure of the check's own rule",
                     "inconclusive (allowed", "finding: a late deadline answer", "a seed phase longer than the budget is a performance fact"):
            self.assertIn(word, outcomes)
        self.assertIn("not a proof that the delay was the server's computation", entry["asserts"])
        self.assertIn("the one failure of the check's own rule (a missing `timing.elapsed_ms` or `limits` field on a direct "
                      "access path fails by the common rule", entry["fails"])
        self.assertNotIn("the one failure the check can give", entry["fails"])
        # the budget-scope table states the clock's start and keeps outside the clock apart from on it (READ-015)
        self.assertIn("which starts when the seed's work starts (its seed phase counts against the budget)", rows[2]["clock"])
        # a label derivation reads the budget after every k-mer and stops at it (SPEC D3; the re-verification L-a)
        self.assertIn("read in a label derivation (after every k-mer) and from depth 1 of the walk", rows[2]["clock"])
        self.assertIn("a label derivation stops at the budget (SPEC D3", rows[2]["not_covered"])
        self.assertTrue(rows[2]["not_covered"].startswith("outside its clock: "), rows[2]["not_covered"])
        self.assertIn("On the clock but never stopped by it: the seed's k-mer mapping and the resolution of named "
                      "labels", rows[2]["not_covered"])
        self.assertTrue(rows[0]["not_covered"].startswith("outside its clock: "), rows[0]["not_covered"])
        self.assertIn("On the clock but never interrupted by it", rows[0]["not_covered"])
        rendered = sections._catalogue_entry(ctx, "budget_respected", entry, "chapters R", True)
        self.assertIn("the answer & the clock read & the threshold & the outcome", rendered)
        self.assertIn("the ruling, case by case", rendered)
        self.assertIn("whether an absence may be read from the returned list", self.checks["checks"]["retrieval_complete"]["title"])
        self.assertIn("an exact count of 0 in the searched scope", self.checks["checks"]["retrieval_complete"]["asserts"])

    def test_timing_names_and_paired_ratios(self):
        """READ-009 to READ-011: the closing section prints Timing 1 .. 7 (its labels unchanged), no prose names a
        C.n section, the percentage off one, and the Timing references are no figures to the number scan."""
        import sections
        ctx = self.ctx(closing=self.closing)
        text = sections.closing(ctx)
        for i, title in ((1, "Observations"), (3, "MCP overhead decomposition"), (7, "Per-chapter timing")):
            self.assertIn("\\section*{Timing %d\\qd %s}" % (i, title), text)
            self.assertIn("\\label{sec:C%d}" % i, text)
        self.assertNotIn("\\section*{C.", text)
        for name in ("closing.json", "overview.json", "glossary.json") + tuple("%s.json" % ch for ch in WRITTEN_CHAPTERS):
            data = json.dumps(sc.read_json(os.path.join(PROSE_DIR, name)))
            self.assertNotRegex(data, r"\bC\.[1-7]\b", name)
        self.assertEqual(sections.fmt_pct_off_one(0.96), "−4 %")
        self.assertEqual(sections.fmt_pct_off_one(1.05), "+5 %")
        self.assertEqual(sections.fmt_pct_off_one(1.0), "0 %")
        self.assertEqual(sections.fmt_pct_off_one(None), "–")
        values = sections.closing_values(ctx)
        self.assertIn("cold_warm_mcp_ratio_pct", values)
        self.assertIn("cold_warm_mcp_ratio_pct", self.closing["placeholders"])
        obs = json.dumps(self.closing["observations"])
        self.assertIn("the median of the per-scenario ratios", obs)
        self.assertIn("their quotient is not that ratio", obs)
        self.assertNotIn("less than half of the MCP overhead", obs)
        self.assertIn("[[mcp_phase_largest_is_terminal_lag ?", obs)
        c7 = " ".join(self.closing["sections"]["C7"]["paragraphs"])
        self.assertIn("no column is the latency of successful searches", c7)
        self.assertIn("their n is smaller than the wall's", c7)
        self.assertIn("they are not paired", c7)
        c6 = " ".join(self.closing["sections"]["C6"]["paragraphs"])
        self.assertIn("100 % is one logical CPU", c6)
        self.assertNotIn("ran threads in parallel", c6)
        self.assertNotIn("the loaded index's residency", c6)
        self.assertIn("Nothing here attributes a slow execution to the load", c6)
        self.assertEqual(literal_numbers("Timing 3 and Timing 1 to 7; see Timing 6"), [])

    def test_by_design_states_and_grouped_note(self):
        """READ-017, READ-018: appendix B's evidence states defined before the table, a case nobody exercised says
        so, and an unverified case's note is grouped by reason with the consequence once."""
        import sections
        ctx = self.ctx()
        ctx.run.table["x-nobody"] = {"areas": ["P"], "agreement": "by_design", "scenario_match": "a request",
                                     "service_answer": {"outcome": "refused_at_submit"}}
        try:
            text = sections.by_design_table(ctx)
        finally:
            del ctx.run.table["x-nobody"]
        self.assertIn("The evidence states of the last column.", text)
        self.assertIn("\\emph{not exercised in this run}", text)
        self.assertIn(latex.tt("x-nobody"), text)
        self.assertIn("not exercised in this run\\\\", text)
        self.assertIn("this run's scenarios and their evidence", text)
        self.assertIn("service (the table's answer)", text)
        note = sections.unverified_note(["M-1: the direct side did not match the case: no direct execution answered http 503",
                                         "M-2: the direct side did not match the case: no direct execution answered http 503",
                                         "N-1 (direct side not verified): r1 outside the case (a pass) has no direct answer"])
        self.assertIn("2 tests (M-1, M-2): the direct side did not match the case", note)
        self.assertIn("1 test (N-1): r1 outside the case (a pass) has no direct answer", note)
        self.assertEqual(note.count("Consequence:"), 1)
        self.assertIn("not in itself a service defect", note)

    def test_exemptions_and_relations_lead(self):
        """READ-019, READ-024, READ-025: the exemption kinds and each entry's consequence, the later-version note, and
        appendix C's lead with the relation examples."""
        import sections
        self.assertIn("not required on mcp; the exemption shows nothing: the recorded outcomes there stand",
                      sections.exemption_consequence("exclusions", {"scenario": "A", "arms": ["mcp"]}))
        self.assertIn("(repeat 1; the other repeats stay required)", sections.exemption_consequence(
            "exclusions", {"scenario": "A", "arms": ["direct-local"], "repeats": [1]}))
        # read from the record (READ-025): what the exempted executions passed stays shown, what they did not stays
        # unproven; a refusal has no property check; the perturbed fixture's R-002 (a late 503 on direct-local)
        run = reader.Run(self.budget_dirs["perturbed"])
        text = sections.exemption_consequence("exclusions", {"scenario": "R-002", "arms": ["direct-local"]}, run)
        self.assertTrue(text.startswith("not required on direct-local; the exemption shows nothing: the recorded outcome "
                                        "there stands (dl r0 finding); "), text)
        self.assertIn("not shown there by this run", text)
        self.assertNotIn("passed there", text)
        props = [c["check"] for l in run.arm_lines("R-002", "direct-local") for c in l["checks"] if run.is_property(c["check"])]
        self.assertTrue(props)
        for name in props:
            self.assertIn(name, text)
        passed = sections.exemption_consequence("exclusions", {"scenario": "P-001", "arms": ["mcp"]}, self.runobj)
        self.assertIn("passed there and stay", passed)
        self.assertNotIn("not shown there", passed)
        outcomes, shown, not_shown, partly = sections.exempted_evidence(self.runobj, {"scenario": "P-001", "arms": ["mcp"]})
        self.assertEqual(len(outcomes), len(self.runobj.arm_lines("P-001", "mcp")))
        self.assertTrue(shown and not not_shown and not partly, (shown, not_shown, partly))
        status_only = _Stub(arm_lines=lambda sid, arm: [{"repeat": 0, "outcome": "pass",
                                                         "checks": [{"check": "status", "outcome": "pass"}]}],
                            is_property=lambda name: name != "status")
        refusal = sections.exemption_consequence("exclusions", {"scenario": "E-1", "arms": ["direct-local"]}, status_only)
        self.assertIn("(dl r0 pass); no property check there: the status, code and text are the claim", refusal)
        self.assertIn("(no execution recorded there)", sections.exemption_consequence(
            "exclusions", {"scenario": "P-002", "arms": ["nowhere"]}, self.runobj))
        self.assertIn("the check walk_in_record is not required on mcp; the scenario's other checks stay required there",
                      sections.exemption_consequence("not_required", {"scenario": "A", "arms": ["mcp"], "check": "walk_in_record"}))
        self.assertIn("the relation is not required on direct-local", sections.exemption_consequence(
            "relations_excluded", {"relation": "REL-x", "arms": ["direct-local"]}))
        ctx = self.ctx(overview=self.overview)
        ctx.run.declared = {"_file": "declared_x.json", "declared_version": 1,
                            "exclusions": [{"scenario": "P-001", "arms": ["mcp"], "reason": "r1 (manifest v5)"}],
                            "not_required": [{"scenario": "P-003", "arms": ["mcp"], "check": "c", "reason": "r2"}],
                            "relations_excluded": [{"relation": "REL-1", "arms": ["direct-local"], "reason": "r3"}]}
        try:
            text = sections.mcp_limits(ctx)
        finally:
            ctx.run.declared = None
        self.assertIn("kind & scenario or relation & access paths & what it changes for this gate & the owner's reason", text)
        self.assertIn("The declared file holds three kinds of entry", text)
        self.assertIn("An exemption changes what the gate requires, not what the run observed", text)
        self.assertIn("what its exempted executions showed stays shown", text)
        self.assertNotIn("remains unproven on that access path by this run", text)
        self.assertIn("A reason that names a later manifest version describes a plan", text)
        self.assertIn("exclusions & \\texttt{P-\\allowbreak{}001} & mcp & not required on mcp", text)
        lead = sections.appendix_lead(ctx, "appendix_C")
        self.assertIn("two answers can share an omission", lead)
        self.assertIn("none in this run", lead)  # the fake run's manifest has no relation
        rel_ctx = self.ctx(overview=self.overview)
        rel_ctx.run.relations = [{"id": "REL-subset-0001", "type": "subset", "a": "P-001", "b": "P-003", "title": "t"}]
        try:
            self.assertIn("REL-subset-0001 (P-001", sections.relations_appendix(rel_ctx))
        finally:
            rel_ctx.run.relations = []
        # the glossary labels the later manifest's mark as not exercised by this run
        entry = [e for e in self.glossary["entries"] if e["term"] == "property-level mark (manifest v5)"][0]
        self.assertIn("[[manifest_version >= 5 ? A mark of this run's manifest. | Not part of this run's manifest", entry["definition"])
        self.assertIn("supplies no evidence about it on that path", entry["definition"])
        rendered = facts.fill(entry["definition"], {"manifest_version": 4})
        self.assertIn("Not part of this run's manifest (version 4)", rendered)
        self.assertNotIn("A mark of this run's manifest", rendered)
        self.assertIn("A mark of this run's manifest.", facts.fill(entry["definition"], {"manifest_version": 5}))

    def test_chapter_purpose_heading_keys_and_access_paths(self):
        """READ-101, READ-105, READ-008, READ-001: the chapter's purpose slot before the Conclusion, the unit
        heading's timing key and the grouped table's n of tests, the access-path wording with the one definition of
        arm, the question-to-chapter map and the orientation before the run's identifiers."""
        import sections
        prose = {"P": {"purpose": ["PURPOSE-SLOT {n_tests} tests"], "intro": ["i"], "summary": ["s"], "coverage_gap": ["g"],
                       "families": {"P-pos": {"title": "t", "prose": "p"}, "P-edge": {"title": "t", "prose": "p"}},
                       "units": {"P-pos-series": "u1", "P-pos-reps": "u2", "P-edge-refused": "u3"}}}
        ctx = self.ctx(prose=prose)
        text = sections.chapter(ctx, "P")
        self.assertIn("PURPOSE-SLOT 5 tests", text)
        self.assertLess(text.find("PURPOSE-SLOT"), text.find("\\parhead{Conclusion}"))
        self.assertIn("\\unithead{P-pos-series}{length}{series · times: wall per repeat, r0 / r1 (mcp: cold / warm)}", text)
        self.assertIn("{grouped · times: wall median / maximum over the group's tests; n = tests}", text)
        self.assertIn("& n (tests) & direct-local & direct-remote & mcp", text)
        self.assertIn("{refusals · no times in the rows}", text)
        self.assertIn(cells.TIMING_LEGEND["replicates"], text)
        self.assertIn("never cold / warm; n counts the group's tests", cells.TIMING_LEGEND["replicates"])
        self.assertIn("the median and the 90th percentile", cells.TIMING_LEGEND["summary"])
        self.assertEqual(set(cells.LAYOUT_KEY), {"series", "list", "pairs", "replicates", "refusals"})
        # no purpose: nothing printed before the Conclusion but the heading
        text2 = sections.chapter(self.ctx(prose={"P": dict(prose["P"], purpose=None)}), "P")
        head = text2.split("\\parhead{Conclusion}", 1)[0]
        self.assertNotIn("PURPOSE", head)
        # R.json carries the example purpose
        R = sc.read_json(os.path.join(PROSE_DIR, "R.json"))
        self.assertTrue(R.get("purpose"))
        self.assertIn("does it know this sequence", " ".join(R["purpose"]))
        # the access paths: one definition of arm in the methods, the traversal arm qualified in T's background
        methods = " ".join(self.overview["sections"]["methods"]["paragraphs"])
        self.assertIn("The suite's records call them arms", methods)
        self.assertIn("a \\*traversal arm\\*, a direction of a walk, is always written with that qualifier".replace("\\*", "*"), methods)
        self.assertIn("mex, the machine that runs the engine's container", methods)
        self.assertIn("the front, the nginx reverse proxy", methods)
        T = sc.read_json(os.path.join(PROSE_DIR, "T.json"))
        self.assertIn("From each requested traversal arm, that is, direction", " ".join(T["intro"]))
        self.assertNotIn("access path", " ".join(T["intro"]).split("traversal arm, that is")[0])
        # the orientation and the map in I.1
        decision = self.overview["sections"]["decision"]
        self.assertTrue(decision["paragraphs"][0].startswith("*What is tested.*"))
        self.assertIn("*The result units.*", decision["paragraphs"][1])
        self.assertIn("*Two ways in.*", decision["paragraphs"][2])
        self.assertTrue(decision["paragraphs"][3].startswith("*Scope.*"))
        self.assertEqual([e["chapter"] for e in decision["chapter_map"]], ["R", "T", "P", "L", "X", "M", "Q", "C", "N", "D"])
        ctx = self.ctx(overview=self.overview)
        table = sections.chapter_map_table(ctx)
        self.assertIn("where to look: the chapters by the question they answer", table)
        for ch in self.units.chapters:
            self.assertIn("\\hyperref[ch:%s]{chapter %s" % (ch, ch), table)
        self.assertIn("Service (\\hyperref[ch:S]", table)  # a chapter without a question prints its title alone
        summary = sections.decision_summary(ctx)
        self.assertLess(summary.find("\\emph{What is tested.}"), summary.find("\\emph{Scope.}"))
        self.assertLess(summary.find("\\emph{Scope.}"), summary.find("where to look: the chapters"))
        self.assertLess(summary.find("where to look: the chapters"), summary.find("\\bline{run verdict}"))
        self.assertIn("MCP, the Model Context Protocol", summary)
        self.assertNotIn("UNRESOLVED", summary)
        # V6-003 (the external review of the v6 document): the result first; the three classifications' full text is
        # appendix F's vocabulary, I.1 names them and points there
        self.assertIn("three classifications, never read for one another (appendix F's vocabulary)", summary)
        self.assertLess(summary.find("\\emph{The result.}"), summary.find("A \\emph{finding} is the suite's category"))
        self.assertIn("one claim is one property check of one scenario", summary)


class ChapterRTReadabilityTest(unittest.TestCase):
    """The readability review's chapter comments of R and T (READ-101 to READ-117, lane RT, 2026-10-16): the evidence
    map after the conclusion's headline (READ-101); the agreement's three levels in the reading line and the classes
    line led by the uncompared property (READ-110); the FINDING line's labelled check note and the document's reading
    of a late 503 (READ-108); no bare 'arm' in the chapters' prose or the generated phrases (READ-008 / READ-112); a
    reader-level question per unit (READ-005); the inventories trimmed to ids (READ-111); and the prose statements of
    READ-103, READ-104, READ-106, READ-107, READ-109, READ-113, READ-115, READ-116 and READ-117."""

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
        cls.T = sc.read_json(os.path.join(PROSE_DIR, "T.json"))

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def ctx(self, prose):
        import sections
        ctx = sections.Context(self.runobj, self.units, prose, {}, {}, list(self.units.chapters))
        for ch in ctx.chapters:
            ctx.footnotes[ch] = cells.Footnotes()
        return ctx

    def test_evidence_map(self):
        """READ-101: with `evidence_map` true the families' claims render as a table after the conclusion's prose
        and before the summary block, each with its family linked and what the direct access paths showed; the rows
        are the records the claim strings are built from; without the opt-in nothing is added."""
        import sections
        base = {"intro": ["i"], "summary": ["s [[claims_unsupported ? Not shown: {claims_unsupported}. | shown]]"],
                "coverage_gap": ["g"],
                "families": {"P-pos": {"title": "t", "prose": "p", "claim": "a window resolves"},
                             "P-edge": {"title": "e", "prose": "p", "claim": "an edge is refused"}},
                "units": {"P-pos-series": "u1", "P-pos-reps": "u2", "P-edge-refused": "u3"}}
        ctx = self.ctx({"P": dict(base, evidence_map=True)})
        text = sections.chapter(ctx, "P")
        self.assertIn("the evidence map: what each family claims", text)
        self.assertIn("\\hyperref[fam:P-pos]{\\texttt{P-pos}} & a window resolves & shown: every test", text)
        self.assertIn("\\hyperref[fam:P-edge]{\\texttt{P-edge}} & an edge is refused & shown: every test", text)
        self.assertLess(text.find("\\parhead{Conclusion}"), text.find("the evidence map"))
        self.assertLess(text.find("the evidence map"), text.find("\\bline{standing}"))
        items = ctx.chapter_values("P")["claims_items"]
        self.assertEqual([(it["family"], it["status"], it["n"], it["m"]) for it in items],
                         [("P-pos", "held", 4, 4), ("P-edge", "held", 1, 1)])
        self.assertNotIn("the evidence map", sections.chapter(self.ctx({"P": dict(base)}), "P"))
        ctx3 = self.ctx({"R": {"evidence_map": True, "intro": ["i"], "summary": ["s"], "coverage_gap": ["g"],
                               "families": {"R-one": {"title": "t", "prose": "p", "claim": "the one"}},
                               "units": {"R-one": "u"}}})
        text3 = sections.chapter(ctx3, "R")
        self.assertIn("\\hyperref[fam:R-one]{\\texttt{R-one}} & the one & not shown: 0 of 1 tests, with a failure or a "
                      "disagreement (I.1)", text3)
        v = ctx3.chapter_values("R")
        self.assertEqual(v["claims_items"][0]["status"], "failed")
        self.assertIn("with a failure or a disagreement", v["claims_unsupported"])
        for key in ("n_mcp_agree_partial", "mcp_agree_partial_example", "footnoted_elsewhere_ids"):
            self.assertIn(key, v)
        self.assertEqual(v["footnoted_elsewhere_ids"], "R-001 (the one)")

    def test_finding_and_classes_lines(self):
        """READ-108, READ-110: the FINDING line labels the check's mechanism as the check's note and adds the
        document's reading of a late 503 (an explanation the contract permits, not a measured cause); the classes
        line of an agreement short of fully checked leads with the property it did not compare."""
        note = ("the deadline is read where the work and the delivery poll it, so an unpolled piece of about that length "
                "ran (the wall also holds the network and any wait for a request thread)")
        items = [("direct-local", 0, "budget_respected",
                  "503 deadline delivered 6130 ms after the 3000 ms budget (client wall 9130 ms)", note)]
        text = blocks.finding_text(items, "the rule")
        self.assertIn("— the check's note: the deadline is read where the work and the delivery poll it", text)
        self.assertIn("the document's reading: what is established is a 503 deadline answer delivered late on the "
                      "client's clock", text)
        self.assertIn("isolates neither the piece nor its length · rule: the rule", text)
        self.assertEqual(text.count("the document's reading"), 1)
        other = blocks.finding_text([("direct-local", 0, "budget_respected", "elapsed 2 s over the 1 s budget",
                                      "another mechanism")], None)
        self.assertIn("— the check's note: another mechanism", other)
        self.assertNotIn("the document's reading", other)
        run = _Stub(comp_entry=lambda sid: {"direct_agreement": {"class": "identical"},
                                            "mcp_agreement": {"class": "agree", "checks_compared": 8}},
                    depth={"T-1": {"class": "partial_non_pass", "repeats_compared": [0, 1],
                                   "not_compared": ["walk_in_record"]}})
        classes = blocks._classes_text(run, "T-1", {})
        self.assertEqual(classes[1], "mcp: the two access paths agree on the checks both evaluated, but walk\\_in\\_record "
                                     "was not compared on any repeat (not decided on both access paths): the agreement does "
                                     "not establish that property (class agree, 8 check(s) compared; partially checked over r0 "
                                     "and r1 (the aggregate not a pass on both access paths))")
        # a check one repeat left undecided while the other decided it on both: established on that repeat alone
        # (the verification of 2026-10-16, chapters issue 3: T-walk-0004's MCP r1 passed walk_in_record)
        run.depth = {"T-1": {"class": "partial_non_pass", "repeats_compared": [0, 1], "not_compared": [],
                             "partly_compared": {"walk_in_record": {"compared": [1], "skipped": [0]}}}}
        self.assertEqual(blocks.agreement_gap_text(run.depth["T-1"]),
                         "the two access paths agree on the checks both evaluated, but walk_in_record was not compared on r0 "
                         "(left undecided there) and decided alike on both access paths on r1: the agreement establishes "
                         "that property on r1 alone")
        # the class is the aggregate over the repeats, said so (the re-verification L-f); the not-compared clause
        # stays where no gap sentence precedes it, and is not repeated after one
        self.assertEqual(blocks.depth_text(run.depth["T-1"]), "partially checked over r0 and r1 (the aggregate not a pass "
                                                              "on both access paths) (not compared on r0: walk_in_record "
                                                              "(compared on r1))")
        self.assertEqual(blocks.depth_text(run.depth["T-1"], gap_said=True),
                         "partially checked over r0 and r1 (the aggregate not a pass on both access paths)")
        self.assertNotIn("(compared on r1))", blocks._classes_text(run, "T-1", {})[1])
        # a generic-only agreement's gap sentence names no check, so its depth keeps the not-compared clause
        run.depth = {"T-1": {"class": "generic_only", "repeats_compared": [0, 1], "not_compared": ["entry_fields"]}}
        self.assertIn("generic checks only on r0 and r1 (not compared: entry\\_fields))", blocks._classes_text(run, "T-1", {})[1])
        run.depth = {"T-1": {"class": "partial_non_pass", "repeats_compared": [0, 1], "not_compared": [],
                             "partly_compared": {"walk_in_record": {"compared": [1], "skipped": [0]}}}}
        self.assertIn("establishes that property on r1 alone", blocks._classes_text(run, "T-1", {})[1])
        # report.agreement_depth tells the two apart: no repeat compared against some repeat compared
        is_prop = lambda name: name == "walk_in_record"
        comp = {"scenarios": {"T-1": {"mcp_agreement": {"class": "agree", "reference": "direct-local",
                                                        "mcp_outcome": "inconclusive", "reference_outcome": "pass"}}}}
        rows = [{"scenario": "T-1", "arm": "direct-local", "repeat": 0, "checks": [{"check": "walk_in_record", "outcome": "pass"}]},
                {"scenario": "T-1", "arm": "mcp", "repeat": 0, "checks": [{"check": "walk_in_record", "outcome": "inconclusive"}]},
                {"scenario": "T-1", "arm": "mcp", "repeat": 1, "checks": [{"check": "walk_in_record", "outcome": "pass"}]}]
        depth = report.agreement_depth(comp, rows, is_prop)["T-1"]
        self.assertEqual(depth["class"], "partial_non_pass")
        self.assertEqual(depth["not_compared"], [])
        self.assertEqual(depth["partly_compared"], {"walk_in_record": {"compared": [1], "skipped": [0]}})
        rows[2]["checks"][0]["outcome"] = "inconclusive"
        depth = report.agreement_depth(comp, rows, is_prop)["T-1"]
        self.assertEqual((depth["class"], depth["not_compared"], depth["partly_compared"]), ("generic_only", ["walk_in_record"], {}))
        run.depth = {"T-1": {"class": "full", "repeats_compared": [0, 1]}}
        self.assertEqual(blocks._classes_text(run, "T-1", {})[1], "mcp: agree; 8 check(s) compared; fully checked on r0 and r1")
        self.assertEqual(blocks.agreement_gap_text({"class": "generic_only"}),
                         "the two access paths agree on the generic checks alone; no property check was compared, so the "
                         "agreement establishes no property")
        self.assertIn("agree on a refusal", blocks.agreement_gap_text({"class": "refusal"}))
        self.assertEqual(blocks.DEPTH_WORD["partial_non_pass"],
                         "partially checked (the aggregate over the repeats not a pass on both access paths)")

    def test_reading_levels_and_prose_words(self):
        """READ-110: the reading line states the three levels of an agreement. READ-008 / READ-112: the chapters' prose
        says access path for a route to the host and traversal arm for a side of a walk, never a bare arm, and so do
        the generated phrases. READ-005: every unit of R and T carries a question. READ-111: the inventories are ids.
        The statements the chapter comments asked for are in the prose data."""
        import facts
        import sections
        ctx = self.ctx({"P": {"intro": ["i"], "summary": ["s"], "coverage_gap": ["g"],
                              "families": {"P-pos": {"title": "t", "prose": "p", "claim": "c"},
                                           "P-edge": {"title": "t", "prose": "p", "claim": "c"}}, "units": {}}})
        levels = sections.agreement_levels_text(ctx, "P")
        self.assertIn("(its three levels, identical answers, the same check outcomes, every property check compared on "
                      "both: I.3)", levels)
        self.assertNotIn("which verifies no property", levels)  # the clause is the count's, printed only with a count
        self.assertIn("its three levels", sections.summary_reading_text(ctx, "P"))
        # I.3 carries the one explanation of the three levels and of what an agreement short of fully checked means
        overview = sc.read_json(os.path.join(PROSE_DIR, "overview.json"))
        depth_par = overview["sections"]["reading"]["paragraphs"][3]
        self.assertIn("an agreement on the generic checks or a refusal alone verifies no property", depth_par)
        self.assertIn("establishes that check on the repeats where both access paths decided it alike", depth_par)
        self.assertIn("one of three things: identical answers (the direct pair's digests)", depth_par)
        # the chapter's agreements short of fully checked by what they compared (facts._mcp_agree_short)
        short = facts._mcp_agree_short(_Stub(comp_entry=lambda sid: {"mcp_agreement": {"class": "agree"}},
                                             depth={"a": {"class": "generic_only"}, "b": {"class": "refusal"},
                                                    "c": {"class": "partial"}, "d": {"class": "partial_non_pass"},
                                                    "e": {"class": "full"}, "f": {"class": "on_failure"}}),
                                       ["a", "b", "c", "d", "e", "f"])
        self.assertEqual((short["n_mcp_agree_generic"], short["n_mcp_agree_partial_checked"], short["n_mcp_agree_on_failure"]), (2, 2, 1))
        self.assertEqual(short["mcp_agree_short_text"],
                         "2 on the generic checks or a refusal alone, which verifies no property; 2 with some property check "
                         "compared and agreeing on both access paths and another not compared on some repeat (the detail "
                         "block names which); 1 with every property check compared on a failing answer")
        self.assertEqual(facts._mcp_agree_short(_Stub(comp_entry=lambda sid: {"mcp_agreement": {"class": "agree"}},
                                                      depth={"e": {"class": "full"}}), ["e"])["mcp_agree_short_text"], "")
        for name, data in (("R.json", self.R), ("T.json", self.T)):
            for path, text in prose_strings(data):
                for m in re.finditer(r"\\barms?\\b", text):
                    before = text[max(0, m.start() - 10):m.start()]
                    self.assertTrue(before.endswith("traversal "), "%s%s: a bare arm in %r" % (
                        name, path, text[max(0, m.start() - 40):m.end() + 20]))
        with open(os.path.join(DOC, "facts.py"), encoding="utf-8") as fh:
            src = fh.read()
        for bad in ('"no exercised arm"', 'on every arm"', '"both direct arms"'):
            self.assertNotIn(bad, src, bad)
        self.assertIn("supplied seed", " ".join(self.T["purpose"]))
        self.assertIn("label support", " ".join(self.R["purpose"]))
        for name, data in (("R", self.R), ("T", self.T)):
            self.assertTrue(data.get("evidence_map"), name)
            for uid, entry in data["units"].items():
                self.assertTrue(isinstance(entry, dict) and entry.get("question"), "%s.json units.%s: no question" % (name, uid))
            self.assertNotIn("{footnoted_text}", json.dumps(data["families"]))
            self.assertNotIn("{footnoted_elsewhere_text}", json.dumps(data["coverage_gap"]))
            self.assertNotIn("{clean_text}", json.dumps(data["families"]))
            self.assertIn("{footnoted_elsewhere_ids}", json.dumps(data["coverage_gap"]))
        R, T = self.R, self.T
        self.assertIn("not established by the builder", R["families"]["R-rc"]["prose"])
        self.assertIn("holding under one stated assumption about the columns' content", R["families"]["R-junction"]["prose"])
        self.assertIn("a mismatch to this oracle", R["families"]["R-junction"]["prose"])
        self.assertIn("accepts two kinds of evidence", R["families"]["R-discover"]["prose"])
        self.assertIn("a 100 bp query at k = 31 has 70 k-mer starts", R["intro"][0])
        self.assertIn("label support, I.1", R["intro"][0])
        self.assertIn("isolated by that measurement", R["families"]["R-budget"]["prose"])
        self.assertIn("no comparison with the full detail", T["families"]["T-graphlet"]["claim"])
        self.assertIn("lifecycle in three stages", T["families"]["T-attempts"]["prose"])
        self.assertIn("A 1 ms budget is no 1 ms response", T["families"]["T-budget"]["prose"])
        self.assertIn("{x[T-budget-0001][dl][checks][budget_respected][detail]}", T["families"]["T-budget"]["prose"])
        self.assertIn("What a pass proves", T["intro"][2])
        self.assertIn("witness leaf", T["families"]["T-exhaustive"]["prose"])
        self.assertIn("one witness route per direction", T["families"]["T-walk"]["claim"])
        self.assertIn("three statements are kept apart".capitalize(), T["intro"][1])


class ChapterPLXTest(unittest.TestCase):
    """The readability review's chapter comments of P, L and X (lane PLX, READ-201 .. READ-215): a grouped cell
    and its coverage line say a by-design member apart from one never exercised, with whether the host was
    reached, and the detail block's observed line names the MCP class (READ-208); the Background's comparison-table
    slot `intro_tables` (READ-209); the prose facts the comments asked for, on both mappings."""

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

    def ctx(self, prose):
        ctx = sections.Context(self.runobj, self.units, prose, {}, {}, list(self.units.chapters))
        for ch in ctx.chapters:
            ctx.footnotes[ch] = cells.Footnotes()
        return ctx

    @staticmethod
    def by_design_state(host_cache, http):
        st = cells.ArmState("mcp")
        st.present = True
        st.glyph_key = "service_by_design"
        st.worst = "service_by_design"
        st.lines = [{"outcome": "service_by_design", "host_cache": host_cache, "http": http, "wall_ms": 50.0}]
        return st

    def test_by_design_apart_from_not_exercised(self):
        """READ-208: a by-design member is said apart from a never-exercised one, with whether the host was reached
        (a relayed host refusal: the host executed the request; a refusal at submit: it did not), in the grouped
        cell and in the coverage line; the detail block's observed line carries the class and that the property was
        not evaluated there."""
        reached = self.by_design_state("warm", 502)
        before = self.by_design_state("not_reached", 400)
        self.assertEqual(cells.by_design_reached_text([reached]), " (the host reached, the pair decided by the service)")
        self.assertEqual(cells.by_design_reached_text([before]), " (decided before the host)")
        self.assertEqual(cells.by_design_reached_text([reached, before]), " (1 the host reached, the rest decided before it)")
        self.assertEqual(cells.by_design_reached_text([]), "")
        passing = cells.scenario_states(self.runobj, "P-001")["mcp"]
        self.assertEqual(passing.worst, "pass")
        states = {"a": {"mcp": passing}, "b": {"mcp": reached}}
        text, _ = cells.aggregate(states, "mcp")
        self.assertTrue(text.startswith("ok 1/1 · by design: 1\\,\\gbd{} (the host reached, the pair decided by the service)"), text)
        self.assertNotIn("not exercised", text)
        self.assertEqual(cells.coverage_text(states, "mcp"),
                         "exercised 1/2, passed 1; by design 1\\,\\gbd{} (the host reached, the pair decided by the service)")
        # a member with no execution recorded stays apart from both
        absent = cells.ArmState("mcp")
        self.assertEqual(cells.coverage_text({"a": {"mcp": passing}, "z": {"mcp": absent}}, "mcp"),
                         "exercised 1/2, passed 1; not recorded 1\\,\\gnone{}")
        st = cells.scenario_states(self.runobj, "P-004")
        block = blocks.detail_block(self.runobj, "P-004", st, 1, self.units)
        self.assertIn("mcp: service-side by design, %s (" % blocks.MCP_CLASS_WORD["refused_before_host"], block)
        self.assertIn("; the property not evaluated there", block)
        st5 = cells.scenario_states(self.runobj, "P-005")
        block5 = blocks.detail_block(self.runobj, "P-005", st5, 2, self.units)
        self.assertIn("mcp: \\attn{not exercised via MCP}, \\attn{%s} (" % blocks.MCP_CLASS_WORD["not_expressible"], block5)

    def test_intro_tables_slot(self):
        """READ-209: a chapter's `intro_tables` entry renders as a plain table after the Background paragraph it
        names (default: the last), its cells prose with placeholders and inline markup."""
        base = {"summary": ["s"], "coverage_gap": ["g"],
                "families": {"P-pos": {"title": "t", "prose": "p"}, "P-edge": {"title": "t", "prose": "p"}},
                "units": {"P-pos-series": "u1", "P-pos-reps": "u2", "P-edge-refused": "u3"}}
        table = {"after": 0, "caption": "cap {n_tests}", "columns": ["a", "b"], "rows": [["`x`", "*y* {n_tests}"]]}
        prose = {"P": dict(base, intro=["FIRST-PARAGRAPH", "SECOND-PARAGRAPH"], intro_tables=[table])}
        bg = sections.chapter(self.ctx(prose), "P").split("\\parhead{Background}", 1)[1]
        self.assertLess(bg.find("FIRST-PARAGRAPH"), bg.find("\\begin{longtable}"))
        self.assertLess(bg.find("\\begin{longtable}"), bg.find("SECOND-PARAGRAPH"))
        self.assertIn("{\\itshape cap 5}", bg)
        self.assertIn("\\texttt{x} & \\emph{y} 5\\\\", bg)
        prose2 = {"P": dict(base, intro=["FIRST-PARAGRAPH", "SECOND-PARAGRAPH"], intro_tables=[dict(table, after=None)])}
        bg2 = sections.chapter(self.ctx(prose2), "P").split("\\parhead{Background}", 1)[1]
        self.assertLess(bg2.find("SECOND-PARAGRAPH"), bg2.find("\\begin{longtable}"))
        prose3 = {"P": dict(base, intro=["ONLY"])}
        self.assertNotIn("longtable", sections.chapter(self.ctx(prose3), "P").split("\\parhead{Background}", 1)[1]
                         .split("\\subsection", 1)[0])

    def test_plx_prose_facts(self):
        """The prose of P, L and X carries what the comments asked for: the chapter purposes (READ-101), the worked
        example before the mode catalogue (READ-201), the scope of a green count (READ-202), the absence unit per
        query class (READ-205), N and the leading-run shortcut (READ-206), the accepted set apart from the observed
        branch (READ-207), the host reached via MCP (READ-208), the three notions of support with the comparison
        table and the hypothetical mosaic (READ-209), the memory account's scope (READ-210), the v4 proof limit of
        L-long-other beside the v5 check (READ-211), the search model before the automaton (READ-212), the
        unconditional qualifier in the unit title of both mappings (READ-213), the three absences of the junction
        (READ-214), the floor as an admission rule (READ-215); every unit of the three chapters has its question."""
        P = sc.read_json(os.path.join(PROSE_DIR, "P.json"))
        L = sc.read_json(os.path.join(PROSE_DIR, "L.json"))
        X = sc.read_json(os.path.join(PROSE_DIR, "X.json"))
        for d in (P, L, X):
            self.assertTrue(d.get("purpose"))
            for uid, e in d["units"].items():
                self.assertTrue(isinstance(e, dict) and e.get("question"), "%s.%s has no question" % (d["chapter"], uid))
        intro = " ".join(P["intro"])
        self.assertIn("k - L + 1 = 7 contexts", intro)
        self.assertIn("not seven organisms or seven biological occurrences", intro)
        self.assertLess(intro.find("A worked example"), intro.find("Mode `count` returns no list"))
        self.assertIn("no check of this chapter verifies an index-wide total independently", P["summary"][0])
        self.assertIn("What they do not establish: the total number of k-mers", P["families"]["P-exact"]["prose"])
        neg = P["families"]["P-negative"]["prose"]
        self.assertIn("`exact` 0 paths says that no walk of the index's retained k-mers spells it", neg)
        self.assertIn("a single k-mer cannot contain a 40-base pattern", neg)
        self.assertIn("not a shorter query", P["families"]["P-iupac"]["prose"])
        caps = P["families"]["P-caps"]["prose"]
        self.assertIn("a set of allowed outcomes, not one branch", caps)
        self.assertIn("{x[P-caps-0004][dl][evidence]}", caps)
        self.assertIn("{x[P-caps-0004][mcp][phrase]}", caps)
        self.assertIn("[[x[P-heavy-0001][mcp][class] = relayed ?", P["families"]["P-heavy"]["prose"])
        self.assertIn("an admission rule, not a search result", P["families"]["P-floor"]["prose"])
        self.assertIn("A refusal therefore proves no absence", P["families"]["P-floor"]["prose"])
        self.assertEqual(len(L["intro_tables"]), 1)
        self.assertEqual(L["intro_tables"][0]["after"], 1)
        self.assertEqual([r[0] for r in L["intro_tables"][0]["rows"]],
                         ["a path (`long_search: paths`)", "label support (`label_intersection`)",
                          "record-verified support (`record_verified`)"])
        self.assertIn("A hypothetical case", L["intro"][1])
        self.assertIn("described, not discriminated", L["intro"][2])
        self.assertIn("not the graph, the process or the container's resident memory", L["families"]["L-memory"]["prose"])
        self.assertIn("{x[L-memory-0001][dl][checks][stop][evidence]}", L["families"]["L-memory"]["prose"])
        self.assertIn("not part of this run's evidence", L["families"]["L-long"]["prose"])
        line = L["units"]["L-long-other"]["property"]
        self.assertIn("no v4 check reads the spelled path", facts.fill(line, {"manifest_version": 4}))
        self.assertNotIn("v4", facts.fill(line, {"manifest_version": 5}).split(";")[0])
        self.assertIn("path_spells on -", facts.fill(line, {"manifest_version": 5}))
        self.assertTrue(X["intro"][0].startswith("A `protein` pattern searches amino acids by their DNA"))
        self.assertIn("no claim about remote homology", X["intro"][0])
        self.assertLess(X["intro"][0].find("exact at the residue level"), X["intro"][0].find("codon automaton"))
        self.assertIn("a property of the query model", X["families"]["X-stop"]["prose"])
        for path in (V4_MAPPING, V5_MAPPING):
            self.assertEqual(sc.read_json(path)["_units"]["X-stop-none"]["unit_title"],
                             "tables without an unconditional stop codon")
        self.assertIn("Three absences are kept apart", X["families"]["X-junction"]["prose"])
        self.assertIn("not an exon-aware spliced alignment", X["families"]["X-junction"]["prose"])

    def test_access_paths_and_agreement_levels(self):
        """READ-008, the chapter-prose part for P, L and X (the merger's reconciliation of the three chapter lanes): a
        route to the host is an access path, never a bare arm, in the three prose files, as in the other chapters'.
        READ-110 (READ-203 merged into it): every chapter's 'Via MCP' sentence says that the agreement is on the checks
        both access paths evaluated and counts the agreements that verify no property (n_mcp_agree_partial), as R's
        and T's do."""
        for ch in ("P", "L", "X"):
            data = sc.read_json(os.path.join(PROSE_DIR, "%s.json" % ch))
            for path, text in prose_strings(data):
                for m in re.finditer(r"\barms?\b", text):
                    before = text[max(0, m.start() - 10):m.start()]
                    self.assertTrue(before.endswith("traversal "), "%s.json%s: a bare arm in %r" % (
                        ch, path, text[max(0, m.start() - 40):m.end() + 20]))
        for ch in ("R", "T", "P", "L", "X", "M", "Q", "N", "D"):
            data = sc.read_json(os.path.join(PROSE_DIR, "%s.json" % ch))
            summary = " ".join(_prose_texts(data.get("summary")))
            self.assertIn("on the checks both access paths evaluated[[mcp_all_full ?", summary, ch)
            self.assertIn("{mcp_agree_short_text}", summary, ch)
            self.assertNotIn("which verifies no property", summary, ch)
            self.assertNotIn("for instance a walk its budget stopped on both sides", summary, ch)
        self.assertIn("(SPEC-pattern-search section 23.2)", " ".join(_prose_texts(sc.read_json(
            os.path.join(PROSE_DIR, "C.json"))["intro"])))
        # the second round's local corrections (READ-212, READ-106, READ-206): synonymous codons, the junction
        # oracle's qualifier at the chapter's opening, the 0-based N positions said so
        X = sc.read_json(os.path.join(PROSE_DIR, "X.json"))
        self.assertIn("differ only by synonymous codons match the same peptide", X["intro"][0])
        self.assertNotIn("third codon positions", X["intro"][0])
        R = sc.read_json(os.path.join(PROSE_DIR, "R.json"))
        self.assertIn("junctions must end where the reference records say (under the one stated assumption about the "
                      "column's other records", " ".join(R["purpose"]))
        P = sc.read_json(os.path.join(PROSE_DIR, "P.json"))
        self.assertIn("positions 10 to 12 (counted from 0", P["families"]["P-iupac"]["prose"])
        # appendix A's traversal checks write the direction as traversal arm (READ-008): a bare arm is the request's
        # own field in a code span, or the SPEC's section title quoted
        checks = sc.read_json(os.path.join(PROSE_DIR, "checks.json"))
        allowed_before = ("traversal ", "access-path ", "access ", "`")
        for path, text in prose_strings(checks):
            for m in re.finditer(r"\b[Aa]n? arms?\b|\barms?\b", text):
                before = text[max(0, m.start() - 12):m.start()]
                word = m.group(0)
                ok = (any(before.endswith(a) for a in allowed_before) or word.startswith(("a ", "an ", "A ", "An ")) and False
                      or text[m.start():m.end() + 10].startswith("arms and steps") or text[m.start() - 1:m.start()] == "`")
                if word in ("arm", "arms", "An arm", "A arm", "an arm", "a arm") and not ok:
                    # a bare direction word: only the SPEC's title "(arms and steps" may stay
                    self.assertTrue(text[m.start():].startswith("arms and steps"), "checks.json%s: a bare arm in %r" % (
                        path, text[max(0, m.start() - 50):m.end() + 30]))



def _v4_branch(text):
    """A prose string as a v4 run reads it: every `[[manifest_version >= 5 ? a | b]]` resolved to b, every other
    conditional's two branches kept with its head dropped, innermost first (READ-313: no v4 legend names v5)."""
    pos = 0
    while True:
        m = _COND.search(text, pos)
        if not m:
            return text
        head, sep, rest = m.group(1).partition("?")
        if not sep:
            pos = m.end()
            continue
        yes, _, no = rest.partition("|")
        chosen = no if head.strip() == "manifest_version >= 5" else (yes + " " + no)
        text = text[:m.start()] + " " + chosen + " " + text[m.end():]
        pos = 0


class ChapterMQCNDTest(unittest.TestCase):
    """The readability review's chapter comments on M, Q, C, N and D (READ-301 to READ-321, 2026-10-16) and the small
    generator parts they needed: a chapter prose slot may hold a table item (sections.prose_blocks; C's map of the
    capabilities documents, READ-311); a detail block takes a reader title and a `reading` note from the chapter
    prose's `details`, filled on the scenario's facts (READ-317, READ-319, READ-320); a pairs table of repeated
    tests prints a three-repeat key and legend (cells.repeats_texts, READ-321); every v4 legend of N and D names no
    v5 check (READ-313); every unit of the five chapters carries a `question`, every chapter a `purpose`; the prose
    says access path, not arm (READ-008)."""

    CHAPTERS = ("M", "Q", "C", "N", "D")

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
        cls.prose = dict((ch, sc.read_json(os.path.join(PROSE_DIR, "%s.json" % ch))) for ch in cls.CHAPTERS)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def ctx(self, prose=None):
        import sections
        ctx = sections.Context(self.runobj, self.units, prose or {}, {}, {}, list(self.units.chapters))
        for ch in ctx.chapters:
            ctx.footnotes[ch] = cells.Footnotes()
        return ctx

    def test_prose_table_item(self):
        """READ-311: a {"table": ...} item of a chapter slot renders as a small table with every cell filled like
        prose; C's background carries the three documents' map; the prose scan reads the table's texts."""
        import sections
        ctx = self.ctx()
        text = sections.prose_blocks(ctx, ["para {n_tests}", {"table": {"caption": "cap {n_tests}", "columns": ["a", "b"],
                                                                        "rows": [["x", "y {n_tests}"]]}}], {"n_tests": 5})
        self.assertIn("para 5", text)
        self.assertIn("\\itshape cap 5}", text)
        self.assertIn("a & b\\\\", text)
        self.assertIn("x & y 5\\\\", text)
        self.assertEqual(sections.prose_blocks(ctx, None, {}), "")
        tables = [p for p in self.prose["C"]["intro"] if isinstance(p, dict)]
        self.assertEqual(len(tables), 1)
        rows = tables[0]["table"]["rows"]
        self.assertEqual([r[0] for r in rows], ["`GET /capabilities`", "`GET /pattern/capabilities`", "`GET /traverse/capabilities`"])
        self.assertEqual(tables[0]["table"]["columns"][0], "route")
        for r in rows:
            self.assertIn("REL-same_answer-000", r[3])
        self.assertIn("never the whole documents", tables[0]["table"]["caption"])
        self.assertIn("REL-same_answer-0004", " ".join(_prose_texts(self.prose["C"]["intro"])))
        # the chapter renders the table between its paragraphs
        prose = {"C": dict(self.prose["C"])}
        prose["C"]["families"] = {}
        self.assertIn("client feature-availability check", " ".join(_prose_texts(prose["C"]["intro"])))
        self.assertIn("None of this chapter's host-document or protocol claims was validated through MCP",
                      self.prose["C"]["coverage_gap"][0])

    def test_detail_title_and_note(self):
        """READ-317, READ-319, READ-320: a detail block's reader title in place of the manifest's and its `reading`
        line after the interpretation, both from the chapter prose's `details` and filled on the scenario's facts."""
        import sections
        st = cells.scenario_states(self.runobj, "R-001")
        block = blocks.detail_block(self.runobj, "R-001", st, 1, self.units, title="TITLE-X", note="NOTE-Y")
        self.assertIn("\\quad TITLE-X}", block)
        self.assertIn("\\bline{reading}{NOTE-Y}", block)
        self.assertLess(block.find("\\bline{unverified}"), block.find("\\bline{reading}"))
        self.assertLess(block.find("\\bline{reading}"), block.find("the technical record"))
        plain = blocks.detail_block(self.runobj, "R-001", st, 1, self.units)
        self.assertNotIn("\\bline{reading}", plain)
        self.assertIn(latex.esc(self.runobj.scenario("R-001").get("title"), 200), plain)
        prose = {"R": {"intro": ["i"], "summary": ["s"], "coverage_gap": ["g"],
                       "families": {"R-one": {"title": "t", "prose": "p"}}, "units": {"R-one": "u"},
                       "details": {"R-001": {"title": "[[x[R-001][dl][standing] = fail ? FAILED-TITLE | OTHER-TITLE]]",
                                             "note": "note on {x[R-001][dl][standing]} ({manifest_version})"}}}}
        ctx = self.ctx(prose)
        text = sections.chapter(ctx, "R")
        self.assertIn("\\quad FAILED-TITLE}", text)
        self.assertNotIn("OTHER-TITLE", text)
        self.assertIn("\\bline{reading}{note on fail (", text)
        self.assertEqual([u for u in ctx.unresolved if "R-001" in u], [])
        # the prose files: N-usecase-0005's title states the current manifest's assertion, N-hard-0022 and
        # D-P-scope-0004 carry a reading
        n_details = self.prose["N"]["details"]
        self.assertEqual(set(n_details), {"N-usecase-0005", "N-hard-0022"})
        title4 = facts.fill(n_details["N-usecase-0005"]["title"], {"manifest_version": 4})
        self.assertIn("the first entry reports forward-only strands", title4)
        self.assertNotIn("both", title4)
        self.assertIn("both entries report + and both are counted", facts.fill(n_details["N-usecase-0005"]["title"], {"manifest_version": 5}))
        self.assertIn("I.2's budget scope", n_details["N-hard-0022"]["note"])
        self.assertIn("not about end-to-end latency", n_details["N-hard-0022"]["note"])
        self.assertEqual(set(self.prose["D"]["details"]), {"D-P-scope-0004"})
        self.assertIn("repeatability and the source property", self.prose["D"]["details"]["D-P-scope-0004"]["note"])

    def test_repeats_key_and_legend(self):
        """READ-321: a pairs table of repeated tests names the repeats in execution order, the MCP cold call as the
        record has it, the server time and the be line per execution, and the digest marks."""
        import sections
        key, legend = cells.repeats_texts(3, True)
        self.assertEqual(key, "pairs · times: wall per repeat, r0 / r1 / r2 (mcp: r0 cold, r1 and r2 warm)")
        self.assertIn("in execution order, r0 / r1 / r2 (via MCP r0 is the host-cold call, r1 and r2 the warm ones)", legend)
        self.assertIn("be the host call as the service timed it, per execution (MCP only)", legend)
        self.assertIn("= and ≠: the normalised answer digests of the repeats on that access path are equal / differ", legend)
        key2, legend2 = cells.repeats_texts(2, False)
        self.assertIn("r0 / r1 (mcp: see the detail blocks)", key2)
        self.assertIn("via MCP the detail blocks state each call's host state", legend2)
        self.assertIn("r1 to r3", cells.repeats_texts(4, True)[0])
        lines = {"A": {"mcp": [{"repeat": 0, "host_cache": "cold"}, {"repeat": 1, "host_cache": "warm"},
                               {"repeat": 2, "host_cache": "warm"}]}}
        stub = _Stub(arm_lines=lambda sid, arm: lines.get(sid, {}).get(arm, []))
        self.assertTrue(sections._mcp_cold_first(stub, ["A"]))
        lines["A"]["mcp"][1]["host_cache"] = "cold"
        self.assertFalse(sections._mcp_cold_first(stub, ["A"]))
        self.assertFalse(sections._mcp_cold_first(stub, ["B"]))
        # a pairs unit of single-repeat tests keeps the plain key and legend
        prose = {"R": {"intro": ["i"], "summary": ["s"], "coverage_gap": ["g"],
                       "families": {"R-one": {"title": "t", "prose": "p"}}, "units": {"R-one": "u"}}}
        text = sections.chapter(self.ctx(prose), "R")
        self.assertIn("\\unithead{R-one}{pair}{pairs · times: wall per repeat, r0 / r1 (mcp: cold / warm)}", text)
        self.assertIn(cells.TIMING_LEGEND["pairs"], text)

    def test_slots_questions_and_access_paths(self):
        """READ-101, READ-005, READ-008 for M, Q, C, N and D: a purpose per chapter, a question per unit, and the
        prose says access path (a traversal arm, a direction, keeps its qualifier)."""
        _, _, units4 = v4()
        for ch in self.CHAPTERS:
            data = self.prose[ch]
            self.assertTrue(data.get("purpose"), ch)
            self.assertIn("I.2's", " ".join(data["purpose"]), ch)
            for uid in [u for f in units4.chapter_families(ch) for u in units4.families[f]["units"]]:
                entry = data["units"].get(uid)
                self.assertIsInstance(entry, dict, (ch, uid))
                self.assertTrue(entry.get("property"), (ch, uid))
                self.assertTrue(entry.get("question", "").rstrip().endswith("?") or "A pass confirms" in entry.get("question", ""), (ch, uid))
            for path, text in prose_strings(data):
                for m in re.finditer(r"\b(\w+) (arms?)\b", text):
                    self.assertEqual(m.group(1), "traversal", "%s%s: %s" % (ch, path, m.group(0)))
        self.assertNotIn("on the direct arms)", inspect.getsource(facts._claim))

    def test_v4_legends_name_no_v5(self):
        """READ-313: read as a v4 run reads them, the legends and prose of N and D name no v5 check; the version
        note names them as not exercised by this run, and the v5 branch reads in the present tense."""
        for ch in ("N", "D"):
            data = self.prose[ch]
            for path, text in prose_strings(data):
                if path == ".coverage_gap[2]" and ch == "N":
                    continue
                v4_text = _v4_branch(text)
                self.assertNotRegex(v4_text, r"\bv5\b|\bversion 5\b|\bmanifest 5\b", "%s%s" % (ch, path))
                self.assertNotIn("from v5", v4_text, "%s%s" % (ch, path))
        note = self.prose["N"]["coverage_gap"][2]
        v4_note = facts.fill(note, {"manifest_version": 4})
        self.assertIn("This run executed manifest version 4", v4_note)
        self.assertIn("none of it is exercised by this run and none of it is evidence about this run", v4_note)
        self.assertIn("both primers of the forward-only pair asserted and counted", v4_note)
        v5_note = facts.fill(note, {"manifest_version": 5})
        self.assertIn("This run executed manifest version 5", v5_note)
        self.assertNotIn("not exercised by this run", v5_note)
        self.assertTrue(self.prose["N"]["intro"][0].startswith("This run executed manifest version {manifest_version}"))
        d_intro = facts.fill(self.prose["D"]["intro"][2], {"manifest_version": 4})
        self.assertIn("a later manifest version marks these per check, which this run does not exercise", d_intro)
        # the unit legends read the current manifest's assertion
        strand = facts.fill(self.prose["N"]["units"]["N-usecase-strand"]["property"], {"manifest_version": 4})
        self.assertIn("its first entry on + only (the second primer is asserted by no check of this manifest)", strand)
        self.assertIn("Neither tests an amplicon", strand)
        peptide = facts.fill(self.prose["N"]["units"]["N-negative-peptide"]["note"], {"manifest_version": 4})
        self.assertIn("The shortest is 10 residues", peptide)

    def test_chapter_sentences(self):
        """The local corrections in words: M's substitution and floor definitions and its batch sentence (READ-302,
        READ-304, READ-305), Q's assumption and strand axes (READ-308, READ-309), N's absence scopes, coordinates,
        repeat counts and defect pin (READ-314 to READ-316, READ-318), D's two dimensions (READ-320)."""
        M = self.prose["M"]
        intro = " ".join(M["intro"])
        self.assertIn("N allows any base, so a position holding N matches every base without spending a substitution", intro)
        self.assertIn("The information floor is an admission rule", intro)
        self.assertIn("not the confidence or significance of any match", intro)
        self.assertIn("What a context counts", intro)
        self.assertIn("as an illustration, not a measurement", intro)
        refusal = M["families"]["M-refusal"]["prose"]
        self.assertIn("[[x[M-refusal-0010][mcp][class] = refused_before_host ?", refusal)
        self.assertIn("a batch holding one unsupported long mismatch pattern returns nothing through the tools", refusal)
        self.assertIn("not a promise of when it will be", refusal)
        self.assertIn("Positive protein recovery", M["summary"][0])
        self.assertIn("admission tests, not recoveries", M["summary"][0])
        Q = self.prose["Q"]
        self.assertIn("An assumption of these fixtures, not a proved fact", Q["families"]["Q-pred"]["prose"])
        self.assertIn("unproved assumption", Q["families"]["Q-pred"]["claim"])
        self.assertIn("Two axes are in play", Q["families"]["Q-strands"]["prose"])
        self.assertIn("A toy example, not an executed test", Q["families"]["Q-motif"]["prose"])
        self.assertIn("What a predicate selects", Q["intro"][0])
        N = self.prose["N"]
        self.assertIn("What an empty answer rules out differs by route", N["intro"][0])
        self.assertIn("excludes that source alone, not every similar sequence", N["intro"][0])
        self.assertIn("two coordinate systems meet here", N["intro"][1])
        self.assertIn("starts 0 to 169", N["intro"][1])
        self.assertIn("What these tests do not measure", N["families"]["N-hard"]["prose"])
        self.assertIn("a count of two is not two copies", N["families"]["N-hard"]["prose"])
        self.assertIn("test no amplicon", N["families"]["N-usecase"]["prose"])
        self.assertTrue(N["units"]["N-usecase-defect"]["property"].startswith("Diagnostic pin of a known index defect"))
        # V6-304 (the external review of the v6 document): a failure is consistent with the repair only under the
        # stated observation, never by itself
        self.assertIn("a failed test alone does not establish the repair", N["units"]["N-usecase-defect"]["note"])
        # the gate relation stated from the policy, never left as a condition for the reader (READ-318)
        note = N["units"]["N-usecase-defect"]["note"]
        self.assertNotIn("Where the gate's selection includes this scenario", note)
        selected = facts.fill(note, {"gate_judged": True, "x": {"N-usecase-0014": {"gate_selected": True}}})
        self.assertIn("The gate's policy selects this scenario", selected)
        self.assertIn("would count as an observed violation until the expectation is revised", selected)
        self.assertIn("The gate's policy does not select this scenario", facts.fill(note, {"gate_judged": True, "x": {"N-usecase-0014": {"gate_selected": False}}}))
        self.assertIn("The deploy gate was not judged for this run", facts.fill(note, {"gate_judged": False, "x": {"N-usecase-0014": {"gate_selected": False}}}))
        # N's version note pointer reads right on either manifest (READ-313)
        self.assertIn("says how this version differs from version 4", facts.fill(N["intro"][0], {"manifest_version": 5}))
        self.assertIn("says what a later version changes in this chapter", facts.fill(N["intro"][0], {"manifest_version": 4}))
        self.assertNotIn("later version", facts.fill(N["intro"][0], {"manifest_version": 5}))
        # M cites the SPEC's step counts without a causal paraphrase (READ-301)
        self.assertNotIn("where the budget binds first", M["summary"][0])
        self.assertIn("section 21.10 states the work per allowance (about 700, 14,000 and 200,000 steps at 1, 2 and 3 "
                      "substitutions on refseq33m", M["summary"][0])
        self.assertIn("does not isolate the allowance's share", M["summary"][0])
        self.assertIn("The evidence gaps.", N["summary"][0])
        self.assertNotIn("{time_limited_text}", N["summary"][0])  # counts, not the id lists (READ-312)
        self.assertNotIn("{late503_text}", N["summary"][0])
        self.assertNotIn("{footnoted_elsewhere_text}", N["coverage_gap"][1])
        D = self.prose["D"]
        self.assertIn("two things, to be read apart", D["intro"][1])
        self.assertIn("validates nothing about the omitted field", D["intro"][1])
        self.assertIn("[[x[D-P-scope-0004][mcp][standing] = inconclusive ?", D["summary"][0])
        for uid in D["units"]:
            self.assertIn("where that path exposes the fields they read", D["units"][uid]["property"])


# ---------------------------------------------------------------------------------------------------------------
# the integration of the readability review with the tints (2026-10-16): the owner's decisions on red (needs
# attention), the printed-value tints, the readability re-verification's V2-M1 and lows

class IntegrationRoundTest(unittest.TestCase):
    """The owner's decisions of 2026-10-16: red means needs attention; a test not implemented via MCP is red wherever
    printed; the contents carry each chapter's attention tag and I.1 an attention table, counted by tests per kind
    (facts.attention_of), the kinds defined once in I.3's legend; every printed execution time tinted by its printed
    value. The readability re-verification: the reading line names a chapter's standing-outcome differences, or none
    exactly when there are none (V2-M1); a /traverse 200 finding has its own kind (L-b); the verdict regex is pinned
    to report.verdict_text (L-i)."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.run_dir = os.path.join(cls.tmp, "run")
        os.makedirs(cls.run_dir)
        tint_fixture(cls.run_dir)
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

    def ctx(self):
        ctx = sections.Context(self.runobj, self.units, {}, {}, {}, list(self.units.chapters))
        for ch in ctx.chapters:
            ctx.footnotes[ch] = cells.Footnotes()
        return ctx

    def test_attention_counts_by_test(self):
        """Each test once per kind, possibly under several kinds, the kinds in their order: a failure, a
        disagreement, a time-budget issue (a wall printed above 55 s here: P-001 r0 and R-001 at 61 s), not
        implemented via MCP (P-005's ○); a wall printed at 30 s (P-004) is no time-budget issue."""
        self.assertEqual(self.rc, 0)
        ctx = self.ctx()
        by = facts.attention_of(ctx, self.units.chapter_scenarios("P"))
        self.assertEqual(by["failure"], ["P-002"])
        self.assertEqual(by["time_budget"], ["P-001"])
        self.assertEqual(by["not_implemented"], ["P-005"])
        self.assertEqual((by["error"], by["finding"]), ([], []))
        v = facts.attention_facts(ctx, self.units.chapter_scenarios("P"))
        self.assertEqual(v["attention_text"], "1 failure, 1 disagreement, 1 time-budget issue, 1 not implemented via MCP")
        self.assertEqual(v["n_attention_tests"], 3)
        self.assertEqual(facts.attention_facts(ctx, self.units.chapter_scenarios("S"))["attention_text"], "")
        self.assertEqual([k for k, _, _ in facts.ATTENTION_KINDS],
                         ["failure", "error", "finding", "disagreement", "time_budget", "not_implemented"])
        self.assertEqual(facts.attention_text({"not_implemented": 9, "time_budget": 2}),
                         "2 time-budget issues, 9 not implemented via MCP")
        # the execution rule: a 503 deadline, a recorded time limitation, the service's deadline, a red wall
        self.assertTrue(facts.time_budget_issue({"http": 503, "error": {"code": "deadline"}, "wall_ms": 100}))
        self.assertTrue(facts.time_budget_issue({"http": 200, "time_limited": "pattern 0 time-limited", "wall_ms": 100}))
        self.assertTrue(facts.time_budget_issue({"http": 504, "inconclusive_reason": "resolve_timeout", "wall_ms": 100}))
        self.assertTrue(facts.time_budget_issue({"http": 200, "wall_ms": 55500}))
        self.assertFalse(facts.time_budget_issue({"http": 200, "wall_ms": 55400}))  # prints 55 s

    def test_contents_tags_table_and_legend(self):
        """The contents: a red tag under each chapter's line with tests needing attention, none for a chapter
        without; I.1: the attention table with the kinds some chapter has, the totals; I.3: one legend line that
        defines every kind."""
        tags = {}
        for ch in self.units.chapters:
            text = self.tex("chapter_%s.tex" % ch)
            m = re.search(r"\\addtocontents\{toc\}\{\\protect\\tocattention\{ch:%s\}\{([^}]*)\}\}" % ch, text)
            tags[ch] = m.group(1) if m else None
        self.assertEqual(tags, {"R": "1 failure, 1 disagreement, 1 time-budget issue",
                                "P": "1 failure, 1 disagreement, 1 time-budget issue, 1 not implemented via MCP",
                                "M": "1 failure, 1 disagreement", "S": None})
        main = self.tex()
        self.assertIn("chapter & \\attn{failures} & \\attn{disagreements} & \\attn{time-budget issues} & "
                      "\\attn{not implemented via MCP} & \\attn{tests needing attention} & tests\\\\", main)
        table = main.split("needs attention: the tests of each chapter", 1)[1].split("\\end{longtable}", 1)[0]
        self.assertNotIn("\\attn{errors}", table)  # a kind no chapter has is left out
        self.assertIn("\\textbf{all chapters} & \\textbf{\\attn{3}} & \\textbf{\\attn{3}} & \\textbf{\\attn{2}} & "
                      "\\textbf{\\attn{1}} & \\textbf{\\attn{5}} & \\textbf{9}\\\\", main)
        self.assertLess(main.index("\\bline{deploy gate}"), main.index("needs attention: the tests of each chapter"))
        self.assertLess(main.index("needs attention: the tests of each chapter"), main.index("\\parhead{Principal findings}"))
        self.assertEqual(main.count("\\bline{attention}{red means needs attention"), 1)
        legend = main.split("\\bline{attention}{", 1)[1].split("}\n", 1)[0]
        for k, _, many in facts.ATTENTION_KINDS:
            self.assertEqual(legend.count("\\attn{%s}:" % many), 1, many)
        self.assertIn("or whose wall prints above 55 s", legend)
        with open(os.path.join(DOC, "templates", "doc.sty"), encoding="utf-8") as fh:
            sty = fh.read()
        self.assertIn("\\definecolor{attnred}{HTML}{C62828}", sty)
        self.assertIn("\\newcommand\\gni{\\leavevmode{\\mdseries\\upshape\\glyphfont\\color{attnred}", sty)
        self.assertIn("\\newcommand\\tocattention[2]", sty)

    def test_not_implemented_red_wherever_printed(self):
        """The ○ count of a standing row, the outcomes and agreement rows, the coverage and aggregate cells, the
        detail block's rows, adjudication and classes, the property marks and the prose's phrases: red."""
        st = dict((sid, cells.scenario_states(self.runobj, sid)) for sid in ("P-001", "P-005"))
        self.assertIn("\\attn{1\\,\\gni{}}", cells.glyph_counts(list(st.values()), "mcp"))
        self.assertNotIn("\\attn{1\\,\\gpass", cells.glyph_counts(list(st.values()), "mcp"))
        chapter = self.tex("chapter_P.tex")
        self.assertIn("\\attn{1 not implemented}", chapter)  # the outcomes row and the agreement line
        self.assertEqual(latex.esc(cells.attn_marked("not exercised via MCP")), "\\attn{not exercised via MCP}")
        self.assertEqual(latex.inline("a !!3 tests could not be exercised!! b"), "a \\attn{3 tests could not be exercised} b")
        self.assertIn("\\attn{not implemented via MCP}", blocks.MARK_GROUP["not_implemented"][1])
        # a conditional compares the printed value without the markers
        self.assertEqual(facts.fill("[[p = not exercised via MCP ? yes | no]]", {"p": cells.attn_marked("not exercised via MCP")}), "yes")

    def test_standing_outcome_differences(self):
        """V2-M1: the reading line names a chapter's tests whose MCP standing differs from its outcome, by kind; it
        prints 'none in this chapter' exactly when there is none (here: the fake run's chapters have none)."""
        ctx = self.ctx()
        for ch in self.units.chapters:
            diffs = sections.standing_outcome_differences(ctx, self.units.chapter_scenarios(ch))
            rows_agree = all(ctx.state(sid)["mcp"].glyph_key == (self.runobj.comp_arm(sid, "mcp") or {}).get("outcome")
                             for sid in self.units.chapter_scenarios(ch) if ctx.state(sid)["mcp"].present)
            self.assertEqual(diffs == [], rows_agree, ch)
            self.assertEqual("(none in this chapter)" in sections.summary_reading_text(ctx, ch), rows_agree, ch)

        # a by-design repeat beside a host-answered one (M-planted-0009 in the v4 run): stands by design, counts mixed
        def state(key, outcomes, marks=()):
            s = cells.ArmState("mcp")
            s.present, s.glyph_key, s.glyph_marks = True, key, list(marks)
            s.lines = [{"repeat": i, "outcome": o} for i, o in enumerate(outcomes)]
            s.agreement = {"case": "pattern.deadline_ends_task"}
            return s
        states = {"M-1": state("service_by_design", ["service_by_design", "pass"]),
                  "M-2": state("service_by_design", ["service_by_design", "inconclusive"]),
                  "M-3": state("not_implemented_mcp", ["service_by_design"]),
                  "M-4": state("pass", ["pass", "pass"])}
        outcomes = {"M-1": "mixed", "M-2": "inconclusive", "M-3": "service_by_design", "M-4": "pass"}
        stub = _Stub(run=_Stub(order=list(states), comp_arm=lambda sid, arm: {"outcome": outcomes[sid]}),
                     state=lambda sid: {"mcp": states[sid]})
        diffs = sections.standing_outcome_differences(stub, list(states))
        self.assertEqual([(sid, kind) for sid, kind, _ in diffs],
                         [("M-1", "by_design"), ("M-2", "by_design"), ("M-3", "not_implemented")])
        self.assertEqual(diffs[0][2], "M-1: via MCP r0 service-side by design, r1 pass, the cell standing by design (the "
                                      "case covers r0) while the outcomes row counts it mixed")
        reading = sections.standing_vs_outcome_reading(stub, list(states))
        self.assertTrue(reading.startswith("in this chapter 3 tests: M-3: the MCP line records service_by_design"), reading)
        self.assertIn("; M-1: via MCP r0 service-side by design, r1 pass", reading)
        self.assertIn("(and 1 more of this kind)", reading)
        self.assertEqual(sections.standing_vs_outcome_reading(stub, ["M-4"]), "none in this chapter")
        self.assertEqual(sections.standing_vs_outcome_example_text(stub, ["M-4"], kind="by_design"), "none in this chapter")

    def test_finding_kind_of_a_walk_and_the_verdict_regex(self):
        """L-b: a /traverse 200 whose budget_respected is a finding is a walk stopped outside its stated rule, not a
        200 past its budget. L-i: the document's FINDINGS regex matches report.verdict_text's own words."""
        lines = {("T-1", "direct-local"): [{"outcome": "finding", "http": 200,
                                            "checks": [{"check": "budget_respected", "outcome": "finding"}]}],
                 ("P-1", "direct-local"): [{"outcome": "finding", "http": 200,
                                            "checks": [{"check": "budget_respected", "outcome": "finding"}]}],
                 ("R-1", "direct-remote"): [{"outcome": "finding", "http": 503, "checks": []}]}
        routes = {"T-1": "POST /traverse", "P-1": "POST /pattern", "R-1": "POST /resolve"}
        run = _Stub(order=["T-1", "P-1", "R-1"], scenario=lambda sid: {"route": routes[sid]},
                    arm_lines=lambda sid, arm: lines.get((sid, arm), []))
        self.assertEqual(sections.finding_kinds(run), {"late_503": 1, "over_200": 1, "walk_200": 1, "other": 0})
        self.assertIn("1 walk stopped outside its stated rule (the seed's clock)", sections.finding_kinds_text(run))
        comp = {"summary": {"totals": {"direct-local": {"executions": 2, "outcomes": {"pass": 1, "finding": 1}}},
                            "direct_agreement": {}, "mcp_agreement": {}, "relations": {}}}
        verdict = report.verdict_text({"arms": {"direct-local": {"state": "done"}}}, comp, [])[0]
        self.assertTrue(sections._FINDINGS_VERDICT.match(verdict), verdict)
        doc = sections.doc_verdict_text(_Stub(verdict=verdict, order=[]))
        self.assertTrue(doc.startswith("FINDINGS (1 execution(s) with a finding, a stated bound broken"), doc)


class ColourRoundTwoTest(unittest.TestCase):
    """The two Opus verifications of e3098b2e9 (integration/verify-merge, integration/verify-colours): decision (1)
    carried out where the text is printed (the class word of what stands ○, the unit notes, appendix B's evidence
    state, the prose's class phrases), the property-level mark not implemented via MCP marked, counted and named, the
    legend stating exactly what is red, a walk's stop figures rounded half up."""

    @staticmethod
    def state(key, marks=(), ni_marked=False):
        st = cells.ArmState("mcp")
        st.present, st.glyph_key, st.glyph_marks, st.ni_marked = True, key, list(marks), ni_marked
        return st

    def test_class_word_red_where_the_arm_stands_not_implemented(self):
        """VERIFY-MERGE M1 (c): the ○ cell's, the observed line's and the adjudication's class word in red; a
        by-design arm's class word black."""
        ni, bd = self.state("not_implemented_mcp"), self.state("service_by_design")
        self.assertEqual(blocks.class_word_tex("not_expressible", ni), "\\attn{not expressible through MCP}")
        self.assertEqual(blocks.class_word_tex("refused_before_host", ni), "\\attn{refused before host execution}")
        self.assertEqual(blocks.class_word_tex("refused_before_host", bd), "refused before host execution")
        self.assertEqual(blocks.class_word_tex("answered"), "answered by the host through the service")
        # the refusals table's mcp cell of a ○ arm: the glyph and its class word red, the reason black
        ni.lines = [{"repeat": 0, "outcome": "not_implemented_mcp", "http": None}]
        ni.agreement = {"class": "not_implemented_mcp"}
        run = _Stub(scenario=lambda sid: {"mcp": {"not_implemented": "graph"}}, cases={})
        ctx = _Stub(run=run, short_reasons=lambda: {})
        cell, note = sections._mcp_edge_cell(ctx, "R-9", ni, None, False)
        # the manifest's own mark without a by-design case: its reason is in appendix G (the delta verifier's L7)
        self.assertEqual(cell, "\\gni\\ \\attn{not expressible through MCP}\\newline the request's graph cannot be "
                               "expressed through the tools; see G")
        self.assertIsNone(note)
        run.scenario = lambda sid: {"mcp": {"not_implemented": "graph", "by_design": ["resolve.not_expressible.graph"]}}
        self.assertTrue(sections._mcp_edge_cell(ctx, "R-9", ni, None, False)[0].endswith("; see B"))
        # the detail block's adjudication phrase
        text = blocks._not_exercised_text(ni, run, "R-9", {})
        self.assertTrue(text.startswith("\\gni\\ \\attn{not exercised via MCP}: the request's graph"), text)
        # the reading line's 'read not exercised' of a refusal filed under not implemented
        m3 = self.state("not_implemented_mcp")
        m3.lines = [{"repeat": 0, "outcome": "service_by_design"}]
        m3.agreement = {"case": "resolve.not_expressible.bounds_max_query_bp"}
        stub = _Stub(run=_Stub(order=["R-1"], comp_arm=lambda sid, arm: {"outcome": "service_by_design"}),
                     state=lambda sid: {"mcp": m3})
        text = latex.esc(sections.standing_outcome_differences(stub, ["R-1"])[0][2])
        self.assertIn("so the cell and the standing count read \\attn{not exercised}, the outcomes row by design", text)

    def test_unit_notes_and_prose_class_phrases_red(self):
        """VERIFY-MERGE M1 (a), (c) and VERIFY-COLOURS L1: every unit note of the prose files and every family prose's
        'Not exercised via MCP' carries the attention markup; inline sets it red."""
        notes = 0
        for ch in ("C", "D", "L", "M", "N", "P", "Q", "R", "T", "X"):
            with open(os.path.join(DOC, "prose", "%s.json" % ch), encoding="utf-8") as fh:
                data = json.load(fh)
            for uid, text in (((data.get("mcp_notes") or {}).get("units")) or {}).items():
                notes += 1
                self.assertTrue(text.startswith("!!Not exercised via MCP!!: "), (ch, uid))
                self.assertTrue(latex.inline(text).startswith("\\attn{Not exercised via MCP}: "), (ch, uid))
            raw = json.dumps(data, ensure_ascii=False)
            self.assertNotIn(" Not exercised via MCP:", raw.replace("!!Not exercised via MCP!!:", ""), ch)
        self.assertEqual(notes, 6)  # C-protocol, L-memory, R-trace, R-budget, T-coords, T-attempts
        # P-modes' cell counts are facts now (VERIFY-COLOURS L3), never the literal 'three cells'
        with open(os.path.join(DOC, "prose", "P.json"), encoding="utf-8") as fh:
            p = fh.read()
        self.assertNotIn("three cells by design or not exercised", p)
        self.assertIn("!!{n_mcp_not_exercised:word} not exercised!!", p)

    def test_appendix_b_count_and_class_counts_red(self):
        """VERIFY-MERGE M1 (b), VERIFY-COLOURS L2: appendix B's 'not implemented mcp: N (...)' evidence state in red;
        a class whose members stand not exercised is a red count in mcp_classes_text, a mixed class names its
        not-exercised members apart."""
        x = {"A": {"mcp": {"class": "not_expressible", "standing": "not exercised", "code": None}, "variant": "a"},
             "B": {"mcp": {"class": "refused_before_host", "standing": "by design", "code": "x"}, "variant": "b"},
             "C": {"mcp": {"class": "refused_before_host", "standing": "not exercised", "code": "x"}, "variant": "c"}}
        v = facts._mcp_class_values(x, ["A", "B", "C"])
        self.assertEqual(latex.esc(v["mcp_classes_text"]), "\\attn{1 not expressible through MCP}, 2 refused before "
                                                           "host execution (\\attn{1 of them not exercised})")
        src = inspect.getsource(sections.by_design_table)
        self.assertIn('named.append(cells.attn(entry) if mark == "not implemented mcp" else entry)', src)

    def test_property_level_mark(self):
        """VERIFY-COLOURS M1, VERIFY-MERGE L1: an MCP execution whose check the manifest's property-level mark not
        implemented via MCP set aside (read not exposed) marks the cell with the small red ○, counts after its glyph
        pair, and is named in the attention tag; a marked check read otherwise sets nothing aside."""
        marks = {"label_dict_contains": {"mark": "not_implemented", "reason": "no dictionary via MCP"}}

        class Run(object):
            def check_marks(self, sid):
                return marks if sid.startswith("T-labels") else {}

            def marked_aside(self, sid, line):
                return dict(((c["check"], 0), self.check_marks(sid)[c["check"]]) for c in line.get("checks") or []
                            if c.get("outcome") == "not_exposed" and c["check"] in self.check_marks(sid))
        run = Run()
        aside = [{"checks": [{"check": "walk_in_record", "outcome": "pass"},
                             {"check": "label_dict_contains", "outcome": "not_exposed"}]}]
        judged = [{"checks": [{"check": "walk_in_record", "outcome": "pass"},
                              {"check": "label_dict_contains", "outcome": "pass"}]}]
        self.assertTrue(cells.property_not_implemented(run, "T-labels-0003", aside))
        self.assertFalse(cells.property_not_implemented(run, "T-labels-0003", judged))
        self.assertFalse(cells.property_not_implemented(run, "T-walk-0001", aside))
        self.assertFalse(cells.property_not_implemented(_Stub(), "T-labels-0003", aside))  # no marks (v1 to v4)
        st = self.state("pass", ni_marked=True)
        self.assertEqual(st.glyph_tex(), "\\gpass\\gnip")
        states = [{"mcp": self.state("pass", ni_marked=True)}, {"mcp": self.state("pass", ni_marked=True)},
                  {"mcp": self.state("pass")}, {"mcp": self.state("not_implemented_mcp")}]
        self.assertEqual(cells.glyph_counts(states, "mcp"),
                         "3\\,\\gpass{}\\ \\attn{(2 with a property check not implemented via MCP)}\\enspace "
                         "\\attn{1\\,\\gni{}}")
        self.assertEqual(facts.attention_text({"not_implemented": 18, "time_budget": 1},
                                              ["T-labels-0003", "T-labels-0005", "T-labels-0006"]),
                         "1 time-budget issue, 18 not implemented via MCP (3 by a property-level mark: T-labels-0003, "
                         "-0005 and -0006)")
        self.assertEqual(facts.attention_text({"not_implemented": 2}), "2 not implemented via MCP")
        # the attention count applies the same rule: a ○ standing, or the mark that set a check aside
        st = self.state("pass")
        st.lines = aside
        self.assertTrue(facts.not_implemented_issue(run, "T-labels-0003", st))
        st.lines = judged
        self.assertFalse(facts.not_implemented_issue(run, "T-labels-0003", st))
        with open(os.path.join(DOC, "templates", "doc.sty"), encoding="utf-8") as fh:
            self.assertIn("\\newcommand\\gnip{\\textsuperscript{", fh.read())
        self.assertIn("property-level mark", sections.MARK_EXAMPLE)

    def test_legend_states_what_is_red(self):
        """VERIFY-MERGE M1 (the legend narrowed to what is red), VERIFY-COLOURS L4 (the footnote digit stated)."""
        legend = sections.attention_legend_sentence()
        self.assertIn("the \\gni\\ glyph of a test not implemented via MCP is red wherever printed", legend)
        self.assertIn("(a count of 0 stays black)", legend)
        self.assertIn("headings, legends, definitions, quoted texts", legend)
        self.assertIn("a footnote digit is dark grey, not red", legend)
        self.assertNotIn("counts, marks and class words are red wherever printed", legend)

    def test_stop_figures_round_half_up(self):
        """VERIFY-COLOURS L5: a walk's stop figures round half up, as every printed time does (15,003.67 ms prints
        15,004 ms, on its tint by that printed value)."""
        line = {"raw": "x"}
        run = _Stub(raw_answer=lambda l: {"body": {"results": [{"resource_stop": {
            "resource": "time", "requested": 15000, "used": 15003.67}}]}})
        self.assertEqual(latex.esc(blocks.resource_stop_phrase(run, line)),
                         "seed 0: the time budget (bounds.time\\_budget\\_ms, 15,000 ms; used \\tyellow{15,004 ms})")
        run = _Stub(raw_answer=lambda l: {"body": {"results": [{"resource_stop": {
            "resource": "time", "requested": 10000, "used": 10000.4}}]}})
        self.assertIn("used 10,000 ms", latex.esc(blocks.resource_stop_phrase(run, line)))



class V6DeltaTest(unittest.TestCase):
    """The v6 results document as a delta on v5 (the owner, 2026-10-17: v6 is very close to v5): what v6 changes
    is a branch on `manifest_version`, so a v5 build stays byte-identical."""

    # the /pattern budget mentions a v6 run keeps at 60,000 ms: the v5 baseline's history, not the budget in force
    V6_KEEPS_60 = ("at the 60,000 ms budget then in force", "ran past its 60,000 ms budget on the staging host",
                   "above the host's default of 60,000", "the /pattern budget of 110,000 ms (60,000 before)")

    def test_pattern_budget_by_version(self):
        """Item 3: chapters P, M, N, Q, X and L state the /pattern budget as 110,000 ms for v6 and 60,000 ms for v5
        (v4 too); their other 60 s figures stay."""
        for ch in ("P", "M", "N", "Q", "X", "L"):
            data = sc.read_json(os.path.join(PROSE_DIR, "%s.json" % ch))
            v5_text, v6_text = [], []
            for path, text in prose_strings(data):
                v6_text.append(facts.fill(text, {"manifest_version": 6}))
                if ".P-warm" not in path:  # v6's own family (the warm-ups at 110,000 ms) is in no v5 document
                    v5_text.append(facts.fill(text, {"manifest_version": 5}))
            v5_text, v6_text = "\n".join(v5_text), "\n".join(v6_text)
            self.assertNotIn("110,000", v5_text, ch)
            self.assertIn("110,000", v6_text, ch)
            rest = v6_text
            for keep in self.V6_KEEPS_60:
                rest = rest.replace(keep, "")
            self.assertNotIn("60,000", rest, "%s: a 60,000 ms /pattern budget left in the v6 text" % ch)
            self.assertIn("60,000", v5_text, ch)

    def test_service_identity_input(self):
        """Item 4: the deployed service's identity from an input file (make_doc --identity), labelled as taken from
        its source (the deployment log) and not recorded by the runner: short in I.4, in full in appendix F; a run
        whose runner recorded none and no input keeps the statement that it is not recorded; a file naming another
        run, or contradicting what the run recorded, is refused."""
        ident = {"run": "R", "service_commit": "60e737f75ff696c170841c983c0464d1157425c9",
                 "railway_deployment_id": "fa013db1-70f0-4dc9-bca1-275c7eefa053",
                 "source": "the deployment log (deploy 2); the MCP runner records no service commit"}
        run = _Stub(mcp_run={"git_head": "9d76035e"}, identity=ident)
        short = sections.service_version_text(run)
        self.assertEqual(short, "the service's deployed revision, taken from the deployment log (deploy 2) and not "
                                "recorded by the runner: service commit 60e737f75ff6…, Railway deployment "
                                "fa013db1-70f0-4dc9-bca1-275c7eefa053")
        self.assertIn("service commit 60e737f75ff696c170841c983c0464d1157425c9,",
                      sections.service_version_text(run, full=True))
        self.assertIn("is not recorded by this run", sections.service_version_text(_Stub(mcp_run={}, identity=None)))
        tmp = tempfile.mkdtemp()
        try:
            run_dir = os.path.join(tmp, "run")
            os.makedirs(run_dir)
            build_run(run_dir)
            merge_run.merge(run_dir)
            sc.write_json_atomic(os.path.join(run_dir, "compare.json"), cmp.compare(run_dir, (CASES, {"cases": 5})))
            run_id = reader.Run(run_dir).run_id()
            path = os.path.join(tmp, "identity.json")
            sc.write_json_atomic(path, dict(ident, run=run_id))
            self.assertEqual(reader.Run(run_dir, identity=path).identity["service_commit"], ident["service_commit"])
            sc.write_json_atomic(path, dict(ident, run="another-run"))
            with self.assertRaises(reader.RunError):
                reader.Run(run_dir, identity=path)
            sc.write_json_atomic(path, dict(ident, run=run_id, manifest="v6, content " + "0" * 64))
            with self.assertRaises(reader.RunError):
                reader.Run(run_dir, identity=path)
        finally:
            shutil.rmtree(tmp)

    def test_service_identity_recorded(self):
        """Manifest v7, item 3: the commit the MCP runner recorded (run.json service_identity, read from the service
        before each invocation) is printed, labelled as recorded; the identity input only as a fallback, its Railway
        deployment id (the service states its commit only) labelled as taken from it; an input whose commit
        contradicts the recorded one is refused."""
        sha = "60e737f75ff696c170841c983c0464d1157425c9"
        recorded = merge_run.identity_summary({}, [{"service_commit": sha, "identity_changed": False,
                                                    "identity_error": None}] * 2, "arms/mcp/mcp_runs.jsonl")
        ident = {"run": "R", "service_commit": sha, "railway_deployment_id": "fa013db1-70f0-4dc9-bca1-275c7eefa053",
                 "source": "the deployment log (deploy 3)"}
        run = _Stub(mcp_run={"git_head": "9d76035e"}, service_identity=recorded, identity=None)
        self.assertEqual(sections.service_version_text(run),
                         "the service's deployed revision as the MCP runner recorded it (read from the service before "
                         "each MCP invocation): service commit 60e737f75ff6…")
        self.assertIn("service commit %s" % sha, sections.service_version_text(run, full=True))
        run.identity = ident
        self.assertEqual(sections.service_version_text(run),
                         "the service's deployed revision as the MCP runner recorded it (read from the service before "
                         "each MCP invocation): service commit 60e737f75ff6…; Railway deployment "
                         "fa013db1-70f0-4dc9-bca1-275c7eefa053, taken from the deployment log (deploy 3) (not recorded "
                         "by the run)")
        # a change, a read that failed, a second commit: stated
        odd = merge_run.identity_summary({}, [{"service_commit": sha, "identity_changed": True},
                                              {"service_commit": "1009f77c7aa5", "identity_changed": None,
                                               "identity_error": "after: timeout"},
                                              {"service_commit": None, "identity_changed": None}], "x")
        text = sections.service_version_text(_Stub(mcp_run={}, service_identity=odd, identity=None))
        for part in ("service commit 60e737f75ff6…, 1009f77c7aa5 (60e737f75ff6… on 1 invocation(s), 1009f77c7aa5 on 1 "
                     "invocation(s))", "changed during 1 invocation(s)", "not read on 1 of 3 invocation(s)",
                     "read once only on 2 invocation(s) (a read failed)"):
            self.assertIn(part, text)
        # nothing recorded: the input file is the fallback, saying so (test_service_identity_input)
        older = merge_run.identity_summary({}, [{}], "arms/mcp/run.json")
        self.assertIn("taken from the deployment log (deploy 3) and not recorded by the runner",
                      sections.service_version_text(_Stub(mcp_run={}, service_identity=older, identity=ident)))
        tmp = tempfile.mkdtemp()
        try:
            run_dir = os.path.join(tmp, "run")
            os.makedirs(run_dir)
            build_run(run_dir)
            path = os.path.join(run_dir, "arms", "mcp", "run.json")
            mcp = sc.read_json(path)
            mcp.update(service_commit=sha, identity_changed=False, identity_error=None)
            sc.write_json_atomic(path, mcp)
            sc.append_jsonl(os.path.join(run_dir, "arms", "mcp", "mcp_runs.jsonl"),
                            {"service_commit": sha, "identity_changed": False, "identity_error": None})
            merge_run.merge(run_dir)
            sc.write_json_atomic(os.path.join(run_dir, "compare.json"), cmp.compare(run_dir, (CASES, {"cases": 5})))
            r = reader.Run(run_dir)
            self.assertEqual(r.service_identity["service_commits"], {sha: 1})
            self.assertTrue(sections.service_version_text(r).startswith(
                "the service's deployed revision as the MCP runner recorded it (read from the service"))
            ipath = os.path.join(tmp, "identity.json")
            sc.write_json_atomic(ipath, dict(ident, run=r.run_id()))
            self.assertIn("Railway deployment fa013db1", sections.service_version_text(reader.Run(run_dir,
                                                                                                  identity=ipath)))
            sc.write_json_atomic(ipath, dict(ident, run=r.run_id(), service_commit="1009f77c7aa5"))
            with self.assertRaises(reader.RunError) as cm:
                reader.Run(run_dir, identity=ipath)
            self.assertIn("service_commit 1009f77c7aa5, the MCP runner recorded %s" % sha, str(cm.exception))
        finally:
            shutil.rmtree(tmp)

    def test_public_paths(self):
        """Item 6: a local path is rewritten (inside the suite's tree from metagraph/scripts/suite/, under a home
        without the home and its hidden directories, a temporary one as its file name); routes, URLs and a lone ~
        stay."""
        pp = latex.public_paths
        self.assertEqual(pp("/Users/u/git/services/metagraph/metagraph/scripts/suite/runs/R/scripts/service_by_design"
                            ".json"), "metagraph/scripts/suite/runs/R/scripts/service_by_design.json")
        self.assertEqual(pp("decoded 2026-10-09 (~/.cache/metagraph-index-defects/r/data/records.tsv, last.tsv)"),
                         "decoded 2026-10-09 (metagraph-index-defects/r/data/records.tsv, last.tsv)")
        self.assertEqual(pp("/private/var/folders/k9/x/T/tmp1/run/compare.json"), "compare.json")
        self.assertEqual(pp("/home/bob/.local/share/f.txt; /tmp/a/b.json"), "share/f.txt; b.json")
        for keep in ("GET /traverse/capabilities", "POST /pattern", "https://metagraph.ethz.ch:8080/x", "named with `~`"):
            self.assertEqual(pp(keep), keep)
        self.assertIsNone(pp(None))

    def test_also_policy_row(self):
        """The coordinator's instruction for the v6 document: a further gate policy (make_doc --also-policy: the
        resolve/traverse policy_v6_rt_mcp beside the deploy gate policy_v6_rtp_mcp) is judged offline and named in
        I.4 with its verdict, never in the deploy gate's place; a policy that cannot be judged says why."""
        entry = {"path": "/x/gate/policy_v6_rt_mcp.json", "policy": {"arms": ["mcp", "direct-local"]}, "error": None,
                 "verdict": {"passed": True, "problems": [], "counts": {"judged": 689, "passed": 635, "required": 689},
                             "policy": {"name": "v6-rt-mcp", "policy_version": 1, "sha256": "1ed3f31501c5" + "0" * 52}}}
        self.assertEqual(sections.also_gate_text(entry),
                         "policy_v6_rt_mcp.json, v6-rt-mcp version 1 (sha256 1ed3f31501c5… of the canonical JSON), "
                         "access paths mcp, direct-local: GATE PASS, 0 problems; 689 executions judged, 635 passed, "
                         "689 required; judged offline by gate.py for this document beside the deploy gate, which it "
                         "does not replace")
        self.assertEqual(sections.also_gate_text({"path": "/x/p.json", "error": "ValueError: bad", "verdict": None}),
                         "p.json: not judged (ValueError: bad)")
        self.assertIn("other gate policy", inspect.getsource(sections.provenance_table)
                      if hasattr(sections, "provenance_table") else inspect.getsource(sections))

    def test_scenario_level_mark_inconclusive_wording(self):
        """The v6 run's D-P-scope-0004 MCP r0 is inconclusive by a scenario-level mark (marks_cover_repeat): I.1 and
        I.4 count it by the mark, never as the service's limit or a host state."""
        run = _Stub(summary=lambda: {"inconclusive_by_reason": {
            "direct-local": {"checks": 4}, "mcp": {"checks": 3, sc.SCENARIO_LEVEL_MARK: 1,
                                                   "traversal_resolve_deadline": 2}}})
        self.assertEqual(sections.inconclusive_by_reason_text(run),
                         "direct-local 4 by the checks; mcp 3 by the checks, 2 by the service's "
                         "traversal_resolve_deadline (a host state), 1 by a scenario-level mark (every property "
                         "check marked as not shown via MCP)")

    def test_answer_side_scenario_mark(self):
        """The coordinator's accuracy fix (2026-10-18): a manifest scenario-level mark (build_manifest.
        scenario_level_marks: property-level marks covering every property check, answer-side) reads as an
        expressible request whose MCP answer lacks the checked field, not implemented yet via MCP (the field named
        from the mark's check), never "cannot be expressed through the tools"; a request-side mark (a field the
        tools have no argument for) keeps the not-expressible wording."""
        mark = ("not implemented yet via MCP: the MCP's /pattern answer does not carry the entry's absence_scope (the "
                "service states the request's scope instead) (entry_fields)")
        scen = {"id": "P-scope-0004", "mcp": {"not_implemented": mark},
                "expect": [{"check": "status", "args": {"http": 200}},
                           {"check": "entry_fields", "args": {"equals": {"absence_scope": "suffix_only"}}}]}
        self.assertEqual(blocks.scenario_mark_parts(scen), [("not_implemented", "entry_fields", mark[29:-15], None)])
        self.assertEqual(blocks.answer_mark_text(scen), "the request is expressible through the tools, but the MCP "
                                                        "answer lacks absence_scope, the field entry_fields checks")
        line = {"scenario": "P-scope-0004", "arm": "mcp", "repeat": 0, "outcome": "not_implemented_mcp", "http": None}
        cls, evidence = blocks.mcp_class([line], scen, {"class": "not_implemented_mcp", "by": "manifest"}, {})
        self.assertEqual(cls, "answer_lacks_field")
        self.assertIn("no MCP call made: the request is expressible", evidence)
        self.assertEqual(blocks.MCP_CLASS_WORD[cls], "not implemented yet via MCP (the answer lacks a checked field)")
        run = _Stub(scenario=lambda sid: scen, cases={})
        reason = blocks.short_reason(run, "P-scope-0004", None)
        self.assertNotIn("cannot be expressed", reason)
        self.assertIn("absence\\_scope", reason)
        # two parts, one of them service-side by design, in the mark's order
        both = dict(scen, mcp={"not_implemented": "service-side by design: the table's case c.x; the MCP answer "
                                                  "cannot carry deterministic; " + mark})
        self.assertEqual([p[:2] for p in blocks.scenario_mark_parts(both)],
                         [("by_design", "deterministic"), ("not_implemented", "entry_fields")])
        self.assertIn("what deterministic checks (service-side by design: the table's case c.x)",
                      blocks.answer_mark_text(both))
        # a request-side mark keeps its wording and class
        req = {"id": "R-x", "mcp": {"not_implemented": "support"}, "expect": []}
        self.assertEqual(blocks.scenario_mark_parts(req), [])
        self.assertIsNone(blocks.answer_mark_text(req))
        self.assertEqual(blocks.mcp_class([dict(line, scenario="R-x")], req, {}, {})[0], "not_expressible")
        self.assertEqual(blocks.short_reason(_Stub(scenario=lambda sid: req, cases={}), "R-x", None),
                         "the request's support cannot be expressed through the tools")

    def test_relation_not_run_via_mcp(self):
        """The coordinator's accuracy fix (2026-10-18): a relation the comparator could not judge on the MCP arm
        (`not_evaluated`, "a raw answer is missing") because one of its scenarios stands not implemented via MCP
        (no MCP execution) reads "not judged via MCP: <id> is not run via MCP (not implemented yet via MCP)" with
        the no-execution dash, in the pairs tables, appendix C and the detail blocks; a relation not evaluated for
        another reason keeps the comparator's words."""
        class Stub(object):
            relations = [{"id": "REL-subset-0002", "type": "subset", "a": "A", "b": "B", "title": "t"}]
            depth = {}

            def property_not_implemented(self, *a):
                return False

            def arm_lines(self, sid, arm):
                if sid == "A" and arm == "mcp":
                    return [{"scenario": "A", "arm": "mcp", "repeat": 0, "outcome": "not_implemented_mcp", "http": None}]
                return [{"scenario": sid, "arm": arm, "repeat": 0, "outcome": "pass", "http": 200}]

            def comp_entry(self, sid):
                return {"mcp_agreement": {"class": "not_implemented_mcp" if sid == "A" else "agree"}}

            def comp_arm(self, sid, arm):
                return {}

            def planned_missing(self, sid, arm):
                return False

            def planned_on(self, sid, arm):
                return True

            def comp_relations(self):
                return {"REL-subset-0002": {"arms": {
                    "direct-local": {"outcome": "pass", "detail": "a's 60 results within b's 640"},
                    "direct-remote": {"outcome": "pass", "detail": "a's 60 results within b's 640"},
                    "mcp": {"outcome": "not_evaluated", "detail": "a raw answer is missing"}}}}

        run = Stub()
        rel = run.relations[0]
        self.assertEqual(cells.relation_unrun(run, rel, "mcp", "not_evaluated"), ["A"])
        self.assertEqual(cells.relation_unrun(run, rel, "mcp", "inconclusive"), [])
        self.assertEqual(cells.relation_unrun(run, rel, "direct-local", "not_evaluated"), [])
        self.assertEqual(cells.relation_unrun_text(["A"]),
                         "not judged via MCP: A is not run via MCP (\\attn{not implemented yet via MCP})")
        self.assertIn(" are not run via MCP", cells.relation_unrun_text(["A", "B"]))
        units = _Stub(unit_of=lambda sid: "U")
        text = blocks._relations_text(run, "B", units)[0]
        self.assertIn("mcp \\gnone{} (not judged via MCP: A is not run via MCP (\\attn{not implemented yet via MCP}))",
                      text)
        self.assertNotIn("a raw answer is missing", text)
        # the tables' mcp column says via MCP itself: the short form there
        self.assertEqual(cells.relation_unrun_text(["A"], short=True),
                         "not judged: A not run (\\attn{not implemented yet via MCP})")
        for source in (inspect.getsource(sections._relation_rows), inspect.getsource(sections.relations_appendix)):
            self.assertIn("cells.relation_unrun_text(unrun, short=True)", source)
        # chapter P's P-scope sentence: both suffix scenarios not run via MCP (v6) -> the relations not judged there;
        # v5's inconclusive branch unchanged
        prose = sc.read_json(os.path.join(PROSE_DIR, "P.json"))["families"]["P-scope"]["prose"]
        ne = {"mcp": {"standing": "not exercised"}}
        filled = facts.fill(prose, {"x": {"P-scope-0004": ne, "P-scope-0006": ne}, "n_relations": 0,
                                    "manifest_version": 6})
        self.assertIn("so via MCP neither is run (!!not implemented yet via MCP!!: the MCP answer lacks that field), and "
                      "the two relations are not judged there.", filled)
        inc = {"mcp": {"standing": "inconclusive"}}
        self.assertIn("they read inconclusive (the check not exposed",
                      facts.fill(prose, {"x": {"P-scope-0004": inc, "P-scope-0006": inc}, "n_relations": 0,
                                         "manifest_version": 5}))

    def test_verify_round_plurals_and_limits(self):
        """The delta verifier's L5 / L6 (2026-10-18): a count carried as text takes its plural from the number; the
        closing section's Timing 4 limits name the population of the unpaired count and word one pair in the
        singular; I.4's declared counts follow their numbers."""
        self.assertEqual(facts.fill("the {n} call{n:s}", {"n": "1"}), "the 1 call")
        self.assertEqual(facts.fill("the {n} call{n:s}", {"n": "12"}), "the 12 calls")
        closing = sc.read_json(os.path.join(PROSE_DIR, "closing.json"))
        text = " ".join(closing["sections"]["C4"]["paragraphs"])
        self.assertNotIn("the long scenarios not run on direct-remote among them", text)
        filled = facts.fill(text, {"scenarios_without_warm_pair": "104", "warm_pairs_left_out_after_interrupted": "1"})
        self.assertIn("Of the scenarios run on both direct access paths, 104 have no host-warm pair of the two, and 1 "
                      "pair warm after an interrupted execution is left out.", filled)
        self.assertIn('"" if dec["relations_excluded"] == 1 else "s"', inspect.getsource(sections))

    def test_unplanned_absence_is_no_execution(self):
        """The delta verifier's M1 (2026-10-18): an arm the run's plan does not schedule for a test (v6's direct-local
        warm-ups on direct-remote and mcp) reads "no execution planned" in the coverage line, the agreement row and
        the classes line, never "not recorded" or compare's bare "missing"; a planned execution not recorded keeps
        its words."""
        def st(present, planned_missing=False):
            a = cells.ArmState("mcp")
            a.present = present
            a.missing_planned = planned_missing
            a.not_planned = not present and not planned_missing
            if present:
                a.glyph_key, a.worst, a.lines = "pass", "pass", [{"outcome": "pass"}]
            return a
        states = {"W-1": {"mcp": st(False)}, "W-2": {"mcp": st(False)}, "X-1": {"mcp": st(False, True)}}
        text = cells.coverage_text(states, "mcp")
        self.assertIn("no execution planned 2\\,\\gnone{}", text)
        self.assertIn("not recorded 1\\,", text)

        class Run(object):
            depth = {}

            def comp_entry(self, sid):
                return {"mcp_agreement": {"class": "missing"}}

            def planned_on(self, sid, arm):
                return arm == "direct-local"

        self.assertEqual(blocks._classes_text(Run(), "W-1", {"mcp": st(False)}),
                         ["mcp: no MCP execution planned (the plan runs this test on direct-local alone)"])
        ctx = _Stub(run=Run())
        self.assertIn(", 2 with no MCP execution planned", sections._agreement_counts(ctx, ["W-1", "W-2"])[1])
        self.assertNotIn("missing", sections._agreement_counts(ctx, ["W-1", "W-2"])[1])

    def test_answer_side_definitions(self):
        """The delta verifier's M2 (2026-10-18): where the manifest has an answer-side scenario-level mark, chapter
        P's coverage sentence counts those tests apart from the ones the tools cannot express, and the definitions
        of not implemented via MCP (I.3's legend, the attention kind, appendix F's vocabulary, the glossary,
        appendix G's lead) name the case; a run without one prints the definitions as before."""
        gap = sc.read_json(os.path.join(PROSE_DIR, "P.json"))["coverage_gap"][0]
        vals = {"n_mcp_not_exercised": 9, "n_mcp_answer_lacks": 2, "n_mcp_not_exercised_tools": 7,
                "mcp_not_exercised_tools_by_family": "modes 1, caps 2, edge inputs 4",
                "mcp_not_exercised_by_family": "scope and strands 2, modes 1, caps 2, edge inputs 4",
                "mcp_answer_lacks_ids": "P-scope-0004, P-scope-0006", "mcp_answer_lacks_fields": "absence_scope"}
        v6 = facts.fill(gap, vals)
        self.assertIn("!!Nine tests are not run via MCP!!: !!seven!! the tools cannot express, having no argument for "
                      "what they send (modes 1, caps 2, edge inputs 4; appendix B)", v6)
        self.assertIn("and !!two!! whose MCP answer lacks the field their property check reads (P-scope-0004, "
                      "P-scope-0006: absence_scope; not implemented yet via MCP", v6)
        v5 = facts.fill(gap, dict(vals, n_mcp_answer_lacks=0, n_mcp_not_exercised=7))
        self.assertIn("!!Seven tests are carried by the direct access paths alone!!, the tools having no argument", v5)
        legend = dict((w, m) for g, w, m in vocabulary.legend_rows(True))["not implemented via MCP"]
        self.assertIn("or, by the manifest's scenario-level mark, the MCP answer lacks every field", legend)
        self.assertNotIn("scenario-level", dict((w, m) for g, w, m in vocabulary.legend_rows())["not implemented via MCP"])
        self.assertIn("by the manifest's scenario-level mark, the MCP answer lacks every field its property checks read",
                      sections.attention_legend_sentence(True))
        self.assertNotIn("scenario-level mark, the MCP answer", sections.attention_legend_sentence())
        meaning = dict((k, m) for k, d, m, c in vocabulary.outcome_table_rows(True))["not_implemented_mcp"]
        self.assertIn("the request is expressible but the MCP answer lacks every field", meaning)
        entry = [e for e in sc.read_json(os.path.join(PROSE_DIR, "glossary.json"))["entries"]
                 if e["term"] == "not exercised via MCP (not implemented)"][0]["definition"]
        self.assertIn("Or, by the manifest's scenario-level mark (P-scope-0004: absence_scope), the request is expressible",
                      facts.fill(entry, {"answer_side_mark": True, "answer_side_mark_text": "P-scope-0004: absence_scope"}))
        self.assertIn("a protocol body). The cell says", facts.fill(entry, {"answer_side_mark": False}))

    def test_limit_appendix(self):
        """The delta verifier's L7 (2026-10-18): a not-exercised MCP arm points at appendix B when a by-design case
        is named, at appendix G when the manifest's own mark is the reason (G lists the manifest's reasons)."""
        st = cells.ArmState("mcp")
        st.agreement = {"class": "not_implemented_mcp", "by": "manifest"}
        run = _Stub(scenario=lambda sid: {"mcp": {"not_implemented": "graphs", "by_design": []}})
        self.assertEqual(blocks.limit_appendix(run, "T-edge-0005", st), "G")
        run = _Stub(scenario=lambda sid: {"mcp": {"not_implemented": "graphs", "by_design": ["resolve.x"]}})
        self.assertEqual(blocks.limit_appendix(run, "R-edge-0017", st), "B")
        st.agreement = {"class": "service_by_design", "case": "pattern.lowercase"}
        run = _Stub(scenario=lambda sid: {"mcp": {}})
        self.assertEqual(blocks.limit_appendix(run, "P-edge-0004", st), "B")

    def test_marked_checks_set_aside_not_uncompared(self):
        """The delta verifier's L9 (2026-10-18): a check the manifest's property-level marks set aside on the MCP arm
        reads "set aside by the marks", neither compared nor skipped, in the reading line's example and the classes
        line; an unmarked check not compared keeps "not compared"."""
        depth = {"class": "partial", "not_compared": ["entry_fields", "answer_fields"]}
        self.assertEqual(blocks.not_compared_text(depth), "not compared: entry_fields, answer_fields")
        self.assertEqual(blocks.not_compared_text(depth, ["entry_fields"]),
                         "not compared: answer_fields; set aside by the marks: entry_fields")
        gap = blocks.agreement_gap_text({"class": "partial", "not_compared": ["entry_fields"]}, ["entry_fields"])
        self.assertIn("entry_fields is set aside by the manifest's marks on the MCP access path", gap)
        self.assertNotIn("not compared", gap)

    def test_p_edge_0004_note(self):
        """Item 5: the P-edge unit's note on P-edge-0004 prints for the frozen v6 manifest only (the service's owner
        states the omission is by design; the next manifest revision corrects the class), so a v5 build prints no
        note line under P-edge (unit_section skips a note that fills to nothing)."""
        note = sc.read_json(os.path.join(PROSE_DIR, "P.json"))["units"]["P-edge"]["note"]
        v6 = facts.fill(note, {"manifest_version": 6})
        self.assertTrue(v6.startswith("P-edge-0004 (lower case): the service's owner states"), v6)
        self.assertIn("by design", v6)
        self.assertIn("the frozen v6 manifest", v6)
        self.assertIn("the next manifest revision corrects it", v6)
        for version in (4, 5):
            self.assertEqual(facts.fill(note, {"manifest_version": version}).strip(), "")
        # the v7 delta: manifest v7 marks the check service-side by design, which the note says instead
        v7 = facts.fill(note, {"manifest_version": 7})
        self.assertIn("via MCP the `entry_fields` check is set aside as service-side by design (manifest v7's", v7)
        self.assertNotIn("the next manifest revision corrects it", v7)
        self.assertIn('note_text = ctx.fill(note, facts.unit_facts(ctx, uid)).strip() if note else ""',
                      inspect.getsource(sections.unit_section))

if __name__ == "__main__":
    unittest.main()
