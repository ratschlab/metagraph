"""Offline tests of manifest v7 (the owner, 2026-10-18): v6 with exactly the approved changes. Item 1: P-edge-0004's
entry_fields check, which reads the pattern text the MCP's /pattern answer leaves out by design (the service owner's
statement, 2026-10-17), is classed service-side by design at property level (mcp.by_design_checks under the table's
pattern.lowercase) instead of not implemented yet via MCP. Item 2: every /pattern scenario's budget_respected states
the host's finalisation reserve (reserve_ms, the capabilities snapshot's pattern finalize_reserve_ms). The tests that
build skip without the record cache (SUITE_TRUTH_CACHE or truth/cache). Item 5 (the owner's scope addition of
2026-10-18): the MCP's answer carries absence_scope and strands per database, so the checks that read them are no
longer marked not implemented via MCP and P-scope-0004 / -0006 run on the MCP arm again."""
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
# the /pattern scenarios of v6 with a budget_respected expectation (every one that allows a 200: 429 + the 25 warm-ups)
N_PATTERN_BUDGET = 454
# the snapshot's pattern_capabilities finalize_reserve_ms (truth/staging_capabilities/, release refseq97-6f945db4)
SNAPSHOT_RESERVE_MS = 250
# item 5: v6's marks of what v7's MCP answer carries per database (test_builder_v6's ABSENCE_SCOPE and STRANDS)
ABSENCE_SCOPE = ["P-exact-%04d" % i for i in range(1, 9)] + ["P-scope-%04d" % i for i in (1, 2, 3, 5)] + \
    ["D-P-exact-0001", "D-P-scope-0004"]
STRANDS = ["P-scope-%04d" % i for i in range(7, 13)] + ["N-usecase-0005"]
UNMARKED_SCENARIOS = ["P-scope-0004", "P-scope-0006"]   # v6's scenario-level marks (absence_scope their only check)


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
        """The v7 the builder makes differs from the frozen v6 in the approved changes only: P-edge-0004's MCP mark
        and budget_respected's reserve_ms on every /pattern scenario that states the check."""
        m, v6 = built(version=7), read("suite_v6.json")
        self.assertEqual(m["manifest_version"], 7)
        s6, s7 = by_id(v6), by_id(m)
        self.assertEqual([x["id"] for x in m["scenarios"]], [x["id"] for x in v6["scenarios"]])
        self.assertEqual(m["relations"], v6["relations"])
        changed, reserve = {}, []
        for sid, s in s6.items():
            for path, old, new in leaf_diffs(s, s7[sid]):
                k = next((i for i, e in enumerate(s["expect"]) if e["check"] == "budget_respected"), None)
                if (path, old, new) == (".expect[%s].args.reserve_ms" % k, None, SNAPSHOT_RESERVE_MS):
                    reserve.append(sid)
                    continue
                changed.setdefault(sid, []).append((path, old, new))
        self.assertEqual(sorted(reserve), sorted(x["id"] for x in v6["scenarios"] if x["route"] == "POST /pattern"
                                                 and any(e["check"] == "budget_respected" for e in x["expect"])))
        self.assertEqual(len(reserve), N_PATTERN_BUDGET)
        what = {"absence_scope": "the entry's absence_scope (the service states the request's scope instead)",
                "strands": "the entry's strands"}
        for ids, field in ((ABSENCE_SCOPE, "absence_scope"), (STRANDS, "strands")):
            for sid in ids:
                self.assertEqual(changed.pop(sid), [(".mcp.not_implemented_checks", [
                    {"check": "entry_fields", "reason": "the MCP's /pattern answer does not carry %s" % what[field]}],
                    None)], sid)
        v5 = by_id(read("suite_v5.json"))
        for sid in UNMARKED_SCENARIOS:
            # the MCP mapping v5 had (its tool, arguments and exposes), with v6's budget
            self.assertNotIn("not_implemented", s7[sid]["mcp"])
            want = dict(v5[sid]["mcp"], arguments=dict(v5[sid]["mcp"]["arguments"], time_budget_ms=110000))
            self.assertEqual(s7[sid]["mcp"], want, sid)
            self.assertEqual(sorted(p for p, _, _ in changed.pop(sid)),
                             [".mcp.arguments", ".mcp.exposes", ".mcp.not_implemented", ".mcp.tool"], sid)
            self.assertEqual(s7[sid]["request"]["time_budget_ms"], s7[sid]["mcp"]["arguments"]["time_budget_ms"])
        self.assertEqual(sorted(changed), sorted(PATTERN_TEXT_BY_DESIGN))
        for sid, case in PATTERN_TEXT_BY_DESIGN.items():
            self.assertEqual(changed[sid], [
                (".mcp.by_design_checks", None, [{"check": "entry_fields", "case": case}]),
                (".mcp.not_implemented_checks", [{"check": "entry_fields", "reason": "the MCP's /pattern answer does "
                                                                           "not carry the entry's pattern text"}],
                 None)])
            self.assertIn(case, s7[sid]["mcp"]["by_design"])
        c6, c7 = json.loads(json.dumps(v6["counts"])), json.loads(json.dumps(m["counts"]))
        self.assertEqual(c6.pop("mcp_not_implemented") - c7.pop("mcp_not_implemented"), len(UNMARKED_SCENARIOS))
        self.assertEqual(c7, c6)
        self.assertEqual(m["mcp_exposure_baseline"]["scenarios_with_exposes"],
                         v6["mcp_exposure_baseline"]["scenarios_with_exposes"] + len(UNMARKED_SCENARIOS))
        # no /pattern check is marked for absence_scope or strands any more; P-edge-0004's by-design mark stays
        marks = [(x["id"], e.get("reason") or e.get("case")) for x in m["scenarios"]
                 for key in ("not_implemented_checks", "by_design_checks") for e in x["mcp"].get(key) or []
                 if x["route"] == "POST /pattern"]
        self.assertEqual(marks, [("P-edge-0004", "pattern.lowercase")])


class PatternReserve(unittest.TestCase):
    """Item 2 (the owner's decision on SD-021): the builder states the host's /pattern finalisation reserve in the
    check arguments (budget_respected's reserve_ms), from the capabilities snapshot it reads."""

    def test_state_pattern_reserve(self):
        shared = {"rule": "r"}
        pat = {"id": "P", "route": "POST /pattern", "expect": [{"check": "status", "args": {"http": 200}},
                                                             {"check": "budget_respected", "args": shared}]}
        res = {"id": "R", "route": "POST /resolve", "expect": [{"check": "budget_respected", "args": shared}]}
        refusal = {"id": "Q", "route": "POST /pattern", "expect": [{"check": "status", "args": {"http": 400}}]}
        scen = [pat, res, refusal]
        self.assertEqual(areas.state_pattern_reserve(scen, {"finalize_reserve_ms": 250}), ["P"])
        self.assertEqual(pat["expect"][1], {"check": "budget_respected", "args": {"rule": "r", "reserve_ms": 250}})
        self.assertEqual(shared, {"rule": "r"})   # the builder's shared args are not mutated
        self.assertEqual(res["expect"], [{"check": "budget_respected", "args": {"rule": "r"}}])
        self.assertEqual(refusal["expect"], [{"check": "status", "args": {"http": 400}}])
        for caps in ({}, {"finalize_reserve_ms": None}, {"finalize_reserve_ms": "250"}, {"finalize_reserve_ms": True},
                     {"finalize_reserve_ms": -1}):
            with self.assertRaises(SystemExit):
                areas.state_pattern_reserve([dict(pat)], caps)

    def test_snapshot_reserve(self):
        caps = bm.load_capabilities(bm.DEFAULT_CAPABILITIES)
        self.assertEqual(caps["pattern_capabilities"]["finalize_reserve_ms"], SNAPSHOT_RESERVE_MS)


class PerDatabaseFields(unittest.TestCase):
    """Item 5: the fields the MCP's answer carries per database from v7 (areas.MCP_PATTERN_ENTRY_CARRIED_V7) are not
    among what it cannot carry; v6's rule still marks them."""

    def test_carried(self):
        self.assertEqual(areas.MCP_PATTERN_ENTRY_CARRIED_V7, ("absence_scope", "strands"))
        import checks
        self.assertEqual(sorted(checks.MCP_PATTERN_DATABASE_FIELDS), sorted(areas.MCP_PATTERN_ENTRY_CARRIED_V7))
        both = pattern_scenario("A", {"absence_scope": "suffix_only", "strands": ["+"]})
        self.assertEqual(areas.pattern_mcp_unexposed(both, areas.MCP_PATTERN_ENTRY_CARRIED_V7), {})
        self.assertEqual(areas.pattern_mcp_unexposed(both), {"entry_fields": [
            "the entry's absence_scope (the service states the request's scope instead)", "the entry's strands"]})
        mixed = pattern_scenario("B", {"strands": ["+"], "limits.max_contexts": 5})
        self.assertEqual(areas.pattern_mcp_unexposed(mixed, areas.MCP_PATTERN_ENTRY_CARRIED_V7),
                         {"entry_fields": ["the entry's limits"]})
        v7, v6 = [copy.deepcopy(both)], [copy.deepcopy(both)]
        self.assertEqual(areas.mark_pattern_mcp_fields(v7, by_design=True, carried=areas.MCP_PATTERN_ENTRY_CARRIED_V7),
                         {})
        self.assertNotIn("not_implemented_checks", v7[0]["mcp"])
        self.assertEqual(areas.mark_pattern_mcp_fields(v6), {"A": ["entry_fields"]})
