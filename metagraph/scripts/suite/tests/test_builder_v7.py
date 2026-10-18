"""Offline tests of manifest v7 (the owner, 2026-10-18): v6 with exactly the approved changes. Item 1: P-edge-0004's
entry_fields check, which reads the pattern text the MCP's /pattern answer leaves out by design (the service owner's
statement, 2026-10-17), is classed service-side by design at property level (mcp.by_design_checks under the table's
pattern.lowercase) instead of not implemented yet via MCP. The tests that build skip without the record cache
(SUITE_TRUTH_CACHE or truth/cache)."""
import copy
import json
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.dirname(HERE)
sys.path.insert(0, SUITE)
sys.path.insert(0, HERE)

import areas  # noqa: E402
import build_manifest as bm  # noqa: E402
import suite_common as sc  # noqa: E402
from test_builder_v5 import built, have_cache  # noqa: E402
from test_builder_v6 import by_id, leaf_diffs  # noqa: E402

PATTERN_TEXT_BY_DESIGN = {"P-edge-0004": "pattern.lowercase"}


def read(name):
    return sc.read_json(os.path.join(SUITE, "manifest", name))


def pattern_scenario(sid, equals, by_design=("pattern.lowercase",)):
    """A /pattern scenario the MCP runs whose entry_fields reads |equals| (and a property the MCP carries)."""
    return {"id": sid, "route": "POST /pattern", "area": "P",
            "mcp": dict({"tool": "pattern_search", "arguments": {}}, **({"by_design": list(by_design)}
                                                                        if by_design else {})),
            "expect": [{"check": "status", "args": {"http": 200}}, {"check": "invariants", "args": {}},
                       {"check": "entry_fields", "args": {"equals": equals}},
                       {"check": "contains_source_context", "args": {}}]}


class PatternTextByDesign(unittest.TestCase):
    """The rule (areas.MCP_PATTERN_ENTRY_BY_DESIGN, pattern_by_design_case): an entry_fields check whose every field the
    MCP's answer cannot carry is left out by design, with a case of the field that the scenario names, is marked
    service-side by design (mcp.by_design_checks); otherwise not implemented yet via MCP, as in v6."""

    def test_the_rule(self):
        self.assertEqual(areas.MCP_PATTERN_ENTRY_BY_DESIGN, {"pattern": ("pattern.lowercase",)})
        text = pattern_scenario("A", {"pattern": "ACGT"})
        unnamed = pattern_scenario("B", {"pattern": "ACGT"}, by_design=())
        mixed = pattern_scenario("C", {"pattern": "ACGT", "strands": "both"})
        other = pattern_scenario("D", {"strands": "both"})
        scen = [copy.deepcopy(x) for x in (text, unnamed, mixed, other)]
        self.assertEqual(areas.mark_pattern_mcp_fields(scen, by_design=True),
                         {"A": ["entry_fields"], "B": ["entry_fields"], "C": ["entry_fields"], "D": ["entry_fields"]})
        a, b, c, d = scen
        self.assertEqual(a["mcp"]["by_design_checks"], [{"check": "entry_fields", "case": "pattern.lowercase"}])
        self.assertNotIn("not_implemented_checks", a["mcp"])
        self.assertEqual(sc.mcp_check_marks(a), {"entry_fields": {"mark": "by_design", "case": "pattern.lowercase"}})
        for x, what in ((b, "the entry's pattern text"), (c, "the entry's pattern text; the entry's strands"),
                        (d, "the entry's strands")):
            self.assertNotIn("by_design_checks", x["mcp"], x["id"])
            self.assertEqual(x["mcp"]["not_implemented_checks"], [
                {"check": "entry_fields", "reason": "the MCP's /pattern answer does not carry %s" % what}], x["id"])
        # v6's rule (by_design False) marks every one not implemented, as the frozen v6 has it
        v6 = [copy.deepcopy(x) for x in (text, unnamed, mixed, other)]
        areas.mark_pattern_mcp_fields(v6)
        self.assertTrue(all("by_design_checks" not in x["mcp"] and x["mcp"]["not_implemented_checks"] for x in v6))
        # a check marked already is left as it is (either mark)
        again = copy.deepcopy(a)
        self.assertEqual(areas.mark_pattern_mcp_fields([again], by_design=True), {})
        self.assertEqual(again["mcp"], a["mcp"])
        # the builder accepts the mark (a check it expects, a case of the table) and refuses one outside the table
        bad = copy.deepcopy(a)
        bad["mcp"]["by_design_checks"] = [{"check": "entry_fields", "case": "no.such_case"}]
        self.assertTrue(any("by_design_checks" in p for p in bm.validate([self.full(bad)])))
        self.assertEqual([p for p in bm.validate([self.full(a)]) if "by_design_checks" in p], [])

    @staticmethod
    def full(s):
        """|s| with the fields the builder's validation reads besides the marks."""
        return dict(s, title="t", request={"patterns": [{"dna": "ACGT", "id": "w"}], "mode": "all_or_count",
                                           "max_contexts": 10, "time_budget_ms": 110000},
                    truth={}, timing_class="fast", repeat=1, tags=[])

    def test_the_versions_the_builder_makes(self):
        self.assertEqual(bm.MANIFEST_VERSION, 7)
        self.assertEqual(bm.BUILDABLE_VERSIONS, (6, 7))
        self.assertEqual(bm.SMOKE_TAG, "smoke-v5")   # v7 keeps v5's smoke selection (the gate's scenarios)
        with self.assertRaises(SystemExit):
            bm.build(bm.DEFAULT_SEED, fetch=False, manifest_version=5)

    @unittest.skipUnless(have_cache(), "no record cache")
    def test_built_v7_against_the_frozen_v6(self):
        """The v7 the builder makes differs from the frozen v6 in the approved changes only: P-edge-0004's MCP mark."""
        m, v6 = built(version=7), read("suite_v6.json")
        self.assertEqual(m["manifest_version"], 7)
        s6, s7 = by_id(v6), by_id(m)
        self.assertEqual([x["id"] for x in m["scenarios"]], [x["id"] for x in v6["scenarios"]])
        self.assertEqual(m["relations"], v6["relations"])
        changed = {}
        for sid, s in s6.items():
            for path, old, new in leaf_diffs(s, s7[sid]):
                changed.setdefault(sid, []).append((path, old, new))
        self.assertEqual(sorted(changed), sorted(PATTERN_TEXT_BY_DESIGN))
        for sid, case in PATTERN_TEXT_BY_DESIGN.items():
            self.assertEqual(changed[sid], [
                (".mcp.by_design_checks", None, [{"check": "entry_fields", "case": case}]),
                (".mcp.not_implemented_checks", [{"check": "entry_fields", "reason": "the MCP's /pattern answer does "
                                                                           "not carry the entry's pattern text"}],
                 None)])
            self.assertIn(case, s7[sid]["mcp"]["by_design"])
        c6, c7 = json.loads(json.dumps(v6["counts"])), json.loads(json.dumps(m["counts"]))
        self.assertEqual(c7, c6)
