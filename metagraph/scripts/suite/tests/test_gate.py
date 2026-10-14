"""The strict deploy gate (gate.py, gate_policy.py) and the plan reconciliation of the report (GPT review finding 1,
verified by VERIFY-A; the owner's strict gate of 2026-10-14): a run fulfilling its declared plan passes, an
intentionally MCP-only one without a direct arm; a missing execution, a scenario absent from every arm, any
inconclusive, a not-implemented or undeclared by-design outcome, a required property not exposed, a disagreement, a
failed relation, a manifest identity mismatch and missing configuration evidence fail."""

import contextlib
import io
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
import gate  # noqa: E402
import gate_policy  # noqa: E402
import merge_run  # noqa: E402
import report  # noqa: E402
import suite_common as sc  # noqa: E402
from test_compare_report import CASES, T0, line, write_arm  # noqa: E402

# an answered case explaining a check (the service answers what the host refuses, as resolve.duplicate_labels)
ANSWERED = {"id": "case.answered_explains", "agreement": "by_design",
            "direct_expectation": {"http": 400, "code": "invalid_request", "text": None},
            "service_answer": {"outcome": "answered", "code": None}, "explains_checks": ["status"]}
TABLE = (dict(list(CASES.items()) + [(ANSWERED["id"], ANSWERED)]),
         {"sha256": "table-sha", "cases": len(CASES) + 1, "version": 1, "service_commit": "x"})
REFUSED = {"code": "pattern_invalid", "text": "more than 64 per search (PATTERN_MAX_PATTERNS_PER_SEARCH)",
           "field": "patterns"}
CAPS = {"document": {"service_settings": {"advanced_enabled": False}}, "sha256": "c", "read_at": "t", "base": "b",
        "mismatches": []}


def scenario(sid, route, repeat=1, by_design=None, not_implemented=None, arguments=None, tool="traverse_fetch",
             refusal=False):
    mcp = {"tool": tool, "arguments": arguments or {}}
    if by_design:
        mcp["by_design"] = by_design
    if not_implemented:
        mcp = {"not_implemented": not_implemented}
    expect = [{"check": "status", "args": {"http": 200}}, {"check": "contains_source_context", "args": {}},
              {"check": "answer_fields", "args": {}}]
    if refusal:
        # a refusal: its status is the claim, no property check
        expect = [{"check": "status", "args": {"http": 400}}]
    if repeat > 1:
        expect.append({"check": "deterministic", "args": {}})
    return {"id": sid, "area": sid.split("-")[0], "title": sid, "route": route, "request": {}, "mcp": mcp,
            "truth": {}, "expect": expect, "timing_class": "fast", "repeat": repeat, "tags": []}


def scenarios():
    return [scenario("T-001", "POST /traverse"), scenario("R-002", "POST /resolve", tool="traverse_resolve"),
            scenario("R-003", "POST /resolve", by_design=["case.refused"], tool="traverse_resolve", refusal=True),
            scenario("T-004", "POST /traverse", not_implemented="strategy.bogus"),
            scenario("T-005", "POST /traverse", arguments={"overrides": {"labels": {"max_seed_labels": 9}}}),
            scenario("P-006", "POST /pattern", tool="pattern_search"),
            scenario("D-007", "POST /resolve", repeat=3, tool="traverse_resolve")]


def relations():
    return [{"id": "REL-1", "type": "same_answer", "a": "T-001", "b": "R-002", "args": {}},
            {"id": "REL-2", "type": "same_answer", "a": "T-001", "b": "P-006", "args": {}}]


def manifest_of(scen=None, rels=None):
    scen = scen if scen is not None else scenarios()
    rels = rels if rels is not None else relations()
    return {"manifest_version": 9, "content_sha256": sc.scenarios_sha256(scen), "scenarios": scen,
            "relations": rels, "relations_sha256": sc.relations_sha256(rels),
            "index": {"database": "test", "index_fp": "fp"}}


def checks_for(sid, rep, prop="pass"):
    out = [{"check": "status", "outcome": "pass", "detail": "ok"},
           {"check": "contains_source_context", "outcome": "pass", "detail": "x"},
           {"check": "answer_fields", "outcome": prop, "detail": "y"}]
    if sid.startswith("D-") and rep > 0:
        out.append({"check": "deterministic", "outcome": "pass", "detail": "same"})
    return out


def mcp_lines(m, skip=(), change=None):
    """Every MCP execution the paired plan has (warm_every 0), passing, R-003 refused by design; |change| maps a
    key (scenario, repeat) to the line's fields to replace."""
    rows, t = [], T0
    for index, s in sc.select_scenarios(m):
        for rep in sc.mcp_repeats(index, s, 0):
            key = (s["id"], rep)
            t += 1
            if key in skip:
                continue
            if (s.get("mcp") or {}).get("not_implemented"):
                r = line(s["id"], "mcp", rep, t, 1.0, outcome="not_implemented_mcp", checks=[])
            elif s["id"] == "R-008":
                # the service answers what the host refuses: the status fails, explained by the named answered case
                r = line(s["id"], "mcp", rep, t, 5.0, outcome="service_by_design",
                         checks=[{"check": "status", "outcome": "fail", "detail": "http 200, expected 400"}])
                r.update(service_by_design_case=ANSWERED["id"], explained_checks=["status"], explained_items=[0])
            elif s["id"] == "R-003":
                r = line(s["id"], "mcp", rep, t, 5.0, outcome="service_by_design", reached=False, http=400,
                         error=dict(REFUSED), checks=[{"check": "status", "outcome": "fail", "detail": "400"}])
                r["service_by_design_case"] = "case.refused"
            else:
                r = line(s["id"], "mcp", rep, t, 10.0, checks=checks_for(s["id"], rep))
            r.update((change or {}).get(key) or {})
            rows.append(r)
    return rows


def direct_lines(m, arm, change=None):
    rows, t = [], T0 + 500 + (0 if arm == "direct-local" else 100)
    for index, s in sc.select_scenarios(m):
        for p in sc.plan_executions([(index, s)], arm, 0)[0]:
            t += 1
            if s["id"] in ("R-003", "R-008"):
                r = line(s["id"], arm, p[2], t, 3.0, http=400, error={"code": "invalid_request", "text": "x"},
                         checks=[{"check": "status", "outcome": "pass", "detail": "400"}])
            else:
                r = line(s["id"], arm, p[2], t, 3.0, checks=checks_for(s["id"], p[2]))
            r.update((change or {}).get((s["id"], p[2])) or {})
            rows.append(r)
    return rows


def bodies(rows):
    """The same answer body for every 200 answer (the relations compare them)."""
    return dict(((r["scenario"], r["repeat"]), {"x": 1}) for r in rows if r.get("http") == 200)


def make_run(d, m, mcp_rows, direct=False, caps=CAPS, plan=True, arm_manifest=None, direct_change=None):
    """A gated run directory as run_suite.sh leaves it: the manifest, the orchestrator's session, plan.json, the
    arms (the MCP arm with its manifest copy and service_capabilities; |direct| True both direct arms, or a tuple of
    them), merged and compared with TABLE."""
    os.makedirs(d, exist_ok=True)
    sc.write_json_atomic(os.path.join(d, "manifest.json"), m)
    direct = ("direct-local", "direct-remote") if direct is True else tuple(direct or ())
    sc.write_json_atomic(os.path.join(d, "orchestrator.json"), {"sessions": [{
        "areas": "", "limit": "", "skip_local": "0" if "direct-local" in direct else "1",
        "skip_remote": "0" if "direct-remote" in direct else "1", "mode": "paired", "warm_every": "0",
        "mcp_cmd": "runner {ids} {out}", "mcp_run": "",
        "loop": {"mode": "paired", "state": "done", "reason": ""}}]})
    arms = ["mcp"] + list(direct)
    if plan:
        gate.write_plan(d, gate.plan_entry(m, None, None, None, 0, arms, "paired", 0))
    run_extra = {"manifest": {"computed_sha256": m["content_sha256"]}}
    if caps is not None:
        run_extra["service_capabilities"] = caps
    write_arm(d, "mcp", mcp_rows, run_extra=run_extra, raw_bodies=bodies(mcp_rows))
    sc.write_json_atomic(os.path.join(d, "arms", "mcp", "manifest.json"), arm_manifest or m)
    if direct:
        for arm in direct:
            rows = direct_lines(m, arm, direct_change)
            write_arm(d, arm, rows, raw_bodies=bodies(rows))
            sc.write_json_atomic(os.path.join(d, "arms", arm, "manifest.json"), m)
    merge_run.merge(d)
    sc.write_json_atomic(os.path.join(d, "compare.json"), cmp.compare(d, TABLE))
    return d


def policy_of(m, **kw):
    kw.setdefault("warm_every", 0)
    return gate_policy.build_policy(m, table=TABLE, manifest_file="m.json", **kw)


class GateCase(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.m = manifest_of()

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def run_dir(self, name="run"):
        return os.path.join(self.tmp, name)

    def problems(self, verdict, kind=None):
        return [p for p in verdict["problems"] if kind is None or p["kind"] == kind]


class PolicyTest(GateCase):
    def test_selection_by_route_and_exclusions(self):
        pol = policy_of(self.m)
        self.assertEqual(pol["selection"]["scenarios"], ["T-001", "R-002", "R-003", "T-004", "T-005", "D-007"])
        ex = dict((e["scenario"], e["rules"]) for e in pol["exclusions"])
        self.assertEqual(ex, {"T-004": ["manifest_not_implemented"], "T-005": ["advanced_path"]})
        req = dict((i["scenario"], i["repeats"]) for i in pol["required"])
        self.assertEqual(req, {"T-001": [0, 1], "R-002": [0, 1], "R-003": [0, 1], "D-007": [0, 1, 2]})
        self.assertEqual([r["relation"] for r in pol["relations"]], ["REL-1"])
        self.assertEqual(pol["relations_excluded"][0]["relation"], "REL-2")
        self.assertEqual(pol["manifest"]["relations_sha256"], self.m["relations_sha256"])
        self.assertEqual(pol["manifest"]["index_fp"], "fp")
        # a deployment with the advanced path: T-005 is required; a named not-expressible case excludes too
        self.assertNotIn("T-005", [e["scenario"] for e in policy_of(self.m, advanced_path=True)["exclusions"]])
        cases = dict(TABLE[0])
        m2 = manifest_of([scenario("R-009", "POST /resolve", by_design=["case.nim"], tool="traverse_resolve")], [])
        pol2 = gate_policy.build_policy(m2, table=(cases, TABLE[1]), warm_every=0)
        self.assertEqual(pol2["exclusions"][0]["rules"], ["table_not_expressible"])
        self.assertEqual(gate.policy_ids(pol), ["T-001", "R-002", "R-003", "D-007"])

    def test_v4_route_policy_counts(self):
        # the frozen v4: 222 scenarios on the resolve / traverse routes (221 POST + GET /traverse/capabilities), 31 the
        # MCP cannot express (22 manifest, 7 more by a named not-expressible case, 2 on the advanced path); an
        # MCP-only policy also leaves out the 4 explaining answered pairs it cannot verify
        m = sc.read_json(os.path.join(SUITE, "manifest", "suite_v4.json"))
        pol = gate_policy.build_policy(m, manifest_file="manifest/suite_v4.json")
        c = pol["counts"]
        self.assertEqual(c["selected"], 222)
        self.assertEqual(c["by_area"], {"C": 1, "D": 6, "N": 41, "R": 81, "T": 93})
        self.assertEqual(c["excluded"], 35)
        self.assertEqual(c["exclusions_by_rule"], {"manifest_not_implemented": 22, "table_not_expressible": 26,
                                                   "advanced_path": 2, "needs_direct_arm": 4})
        self.assertEqual(sorted(e["scenario"] for e in pol["exclusions"] if "needs_direct_arm" in e["rules"]),
                         ["C-caps-0001", "R-edge-0009", "T-edge-0007", "T-edge-0008"])
        self.assertEqual(c["required"], {"mcp": {"scenarios": 187, "executions": 380}})
        self.assertEqual([r["relation"] for r in pol["relations"]], ["REL-same_answer-0003"])
        # the deploy gate's two arms with the owner's declared entries
        two = gate_policy.build_policy(m, arms=("mcp", "direct-local"), manifest_file="manifest/suite_v4.json",
                                       declared=sc.read_json(os.path.join(SUITE, "gate", "declared_v4_rt.json")))
        c = two["counts"]
        self.assertEqual(c["exclusions_by_rule"], {"manifest_not_implemented": 22, "table_not_expressible": 26,
                                                   "advanced_path": 2, "declared": 4})
        self.assertEqual(c["required"], {"mcp": {"scenarios": 189, "executions": 384},
                                         "direct-local": {"scenarios": 220, "executions": 253}})
        self.assertEqual([r["relation"] for r in two["relations"]], ["REL-same_answer-0003"])
        self.assertEqual(sorted(r["relation"] for r in two["relations_excluded"] if "declared" in (r.get("rules") or ())),
                         ["REL-prefix_equals-%04d" % i for i in (1, 2, 3, 4)])
        self.assertEqual(sorted((n["scenario"], n["check"]) for n in two["properties"]["not_required"]),
                         [("T-edge-0008", "answer_fields"), ("T-labels-0003", "label_dict_contains"),
                          ("T-labels-0005", "label_dict_contains")])

    def test_query_cap_not_expressible_is_moot_below_the_cap(self):
        # the service always sends its own bounds.max_query_bp (2,000 bp with labels, 500 for a discovery): a query
        # no longer than that cap, whose direct request states no bound, is the same request through the MCP
        m = sc.read_json(os.path.join(SUITE, "manifest", "suite_v4.json"))
        pol = gate_policy.build_policy(m, manifest_file="manifest/suite_v4.json")
        ex = dict((e["scenario"], e) for e in pol["exclusions"])
        for sid in ("R-edge-0012", "R-edge-0013"):
            self.assertNotIn(sid, ex)
        for sid in ("R-budget-0001", "R-budget-0002", "R-edge-0003", "R-twin-0002", "R-twin-0005"):
            self.assertEqual(ex[sid]["rules"], ["table_not_expressible"], sid)
        # R-edge-0004 states the bound itself (and the manifest marks it)
        self.assertIn("table_not_expressible", ex["R-edge-0004"]["rules"])
        case = {"id": "resolve.not_expressible.bounds_max_query_bp", "agreement": "not_implemented_mcp",
                "service_answer": {"outcome": "not_expressible", "missing_argument": "bounds.max_query_bp"}}
        s = scenario("R-010", "POST /resolve", by_design=[case["id"]], tool="traverse_resolve",
                     arguments={"sequence": "A" * 501, "discover": True})
        s["request"] = {"sequence": "A" * 501, "discover": {}}
        self.assertIsNotNone(gate_policy.mcp_exclusion(s, {case["id"]: case}))   # above the discovery cap
        s["mcp"]["arguments"]["sequence"] = "A" * 500
        self.assertIsNone(gate_policy.mcp_exclusion(s, {case["id"]: case}))
        s["request"]["bounds"] = {"max_query_bp": 100}
        self.assertIsNotNone(gate_policy.mcp_exclusion(s, {case["id"]: case}))   # the direct request states a bound

    def test_basic_override_list_comes_from_the_table(self):
        # traverse.advanced_path_disabled's request predicate (overrides. not_in [...]) is the one list; the module's
        # copy of the service's list is only the fallback of a table without the case
        args = {"overrides": {"direction": "right", "bounds": {"max_extension_bp": 9}}}
        self.assertEqual(gate_policy.advanced_keys(args), [])
        case = {"id": "traverse.advanced_path_disabled", "agreement": "not_implemented_mcp",
                "request": [{"service": "service_settings.advanced_enabled", "equals": False},
                            {"field": "overrides.", "not_in": ["", "bounds.*", "exhaustive", "labels"]}]}
        allowed = gate_policy.basic_overrides({case["id"]: case})
        self.assertEqual(allowed, ["", "bounds.*", "exhaustive", "labels"])
        self.assertEqual(gate_policy.advanced_keys(args, allowed), ["overrides.direction"])
        s = scenario("T-010", "POST /traverse", arguments=args)
        ex = gate_policy.mcp_exclusion(s, {case["id"]: case})
        self.assertEqual(ex["rules"], ["advanced_path"])
        self.assertIn("overrides.direction", ex["reason"])
        real = gate_policy.basic_overrides(cmp.load_by_design(cmp.DEFAULT_BY_DESIGN)[0])
        self.assertIn("labels.seed_label_kind", real)
        self.assertEqual(gate_policy.advanced_keys({"overrides": {"labels": {"max_seed_labels": 1}}}, real),
                         ["overrides.labels.max_seed_labels"])

    def test_shipped_v4_policy_is_regenerated_exactly(self):
        shipped = sc.read_json(os.path.join(SUITE, "gate", "policy_v4_rt_mcp.json"))
        m = sc.read_json(os.path.join(SUITE, "manifest", "suite_v4.json"))
        fresh = gate_policy.build_policy(m, arms=("mcp", "direct-local"), manifest_file="manifest/suite_v4.json",
                                         name="v4-rt-mcp", declared_file="gate/declared_v4_rt.json",
                                         declared=sc.read_json(os.path.join(SUITE, "gate", "declared_v4_rt.json")))
        for p in (shipped, fresh):
            p.pop("generated_utc")
        self.assertEqual(shipped, fresh)

    def test_cli_writes_policy_and_ids(self):
        mpath, pol, ids = (os.path.join(self.tmp, n) for n in ("m.json", "p.json", "ids.txt"))
        sc.write_json_atomic(mpath, self.m)
        with contextlib.redirect_stdout(io.StringIO()) as out:
            self.assertEqual(gate_policy.main([mpath, "--by-design", "", "--warm-every", "0", "--out", pol,
                                               "--ids-out", ids]), 0)
        self.assertIn("excluded: 2 (advanced_path 1, manifest_not_implemented 1)", out.getvalue())
        self.assertIn("required on mcp: 4 scenario(s), 9 execution(s)", out.getvalue())
        self.assertEqual(gate.read_ids(ids, self.m), ["T-001", "R-002", "R-003", "D-007"])
        self.assertEqual(gate.load_policy(pol)["counts"]["selected"], 6)
        # a manifest whose stated hash does not verify is refused
        bad = dict(self.m, relations_sha256="stale")
        sc.write_json_atomic(mpath, bad)
        with contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(gate_policy.main([mpath, "--by-design", ""]), 2)

    def test_read_ids(self):
        path = os.path.join(self.tmp, "ids")
        with open(path, "w") as fh:
            fh.write("# the gate\nD-007, T-001\nR-002  # resolve\n")
        self.assertEqual(gate.read_ids(path, self.m), ["T-001", "R-002", "D-007"])
        with open(path, "w") as fh:
            fh.write("T-001 T-999\n")
        with self.assertRaises(ValueError) as ctx:
            gate.read_ids(path, self.m)
        self.assertIn("T-999", str(ctx.exception))


def with_answered(m):
    """|m| with R-008, a refusal of the host that the service answers by design (an explaining answered case)."""
    scen = m["scenarios"] + [scenario("R-008", "POST /resolve", by_design=[ANSWERED["id"]], tool="traverse_resolve",
                                      refusal=True)]
    return manifest_of(scen, m["relations"])


class TwoArmAndDeclaredTest(GateCase):
    """The answered by-design pairs need a direct arm: an MCP-only policy excludes them (needs_direct_arm, from the
    manifest and the table), the two-arm policy (mcp + direct-local) verifies them; owner-declared exclusions and
    MCP-not-required properties come from a versioned file with a reason each."""

    def setUp(self):
        super(TwoArmAndDeclaredTest, self).setUp()
        self.m = with_answered(self.m)

    def declared(self, **kw):
        d = {"declared_version": 1, "manifest": {"content_sha256": self.m["content_sha256"]},
             "exclusions": [{"scenario": "T-001", "arms": ["mcp"], "reason": "its capabilities document lacks release"}],
             "not_required": [{"scenario": "R-002", "arms": ["mcp"], "check": "answer_fields",
                               "reason": "the MCP never exposes the field"}]}
        d.update(kw)
        return d

    def test_mcp_only_policy_excludes_answered_explaining_pairs(self):
        pol = policy_of(self.m)
        ex = dict((e["scenario"], e) for e in pol["exclusions"])
        self.assertEqual(ex["R-008"]["rules"], ["needs_direct_arm"])
        self.assertIn(ANSWERED["id"], ex["R-008"]["reason"])
        v = gate.judge(make_run(self.run_dir(), self.m, mcp_lines(self.m)), pol)
        self.assertTrue(v["passed"], gate.render(v))

    def test_mcp_only_run_fails_an_answered_pair_it_cannot_verify(self):
        # the same run judged by a policy that requires R-008 on the MCP arm without a direct arm
        pol = policy_of(self.m)
        pol["exclusions"] = [e for e in pol["exclusions"] if e["scenario"] != "R-008"]
        pol["required"].append({"scenario": "R-008", "arm": "mcp", "repeats": [0, 1]})
        v = gate.judge(make_run(self.run_dir(), self.m, mcp_lines(self.m)), pol)
        self.assertIn("R-008", [p["scenario"] for p in v["problems"]], gate.render(v))

    def test_two_arm_policy_verifies_answered_pairs(self):
        pol = policy_of(self.m, arms=("mcp", "direct-local"))
        self.assertNotIn("R-008", [e["scenario"] for e in pol["exclusions"]])
        self.assertIn({"scenario": "R-008", "arm": "direct-local", "repeats": [0]}, pol["required"])
        v = gate.judge(make_run(self.run_dir(), self.m, mcp_lines(self.m), direct=("direct-local",)), pol)
        self.assertTrue(v["passed"], gate.render(v))
        # without the direct arm the two-arm policy's direct-local executions are missing
        v = gate.judge(make_run(self.run_dir("b"), self.m, mcp_lines(self.m)), pol)
        self.assertIn("missing", [p["kind"] for p in v["problems"]])

    def test_declared_entries(self):
        pol = policy_of(self.m, arms=("mcp", "direct-local"), declared=self.declared())
        ex = dict((e["scenario"], e) for e in pol["exclusions"])
        self.assertEqual(ex["T-001"]["rules"], ["declared"])
        self.assertEqual(ex["T-001"]["arms"], ["mcp"])
        self.assertEqual(pol["properties"]["not_required"][0]["check"], "answer_fields")
        self.assertEqual(pol["declared"]["declared_version"], 1)
        text = gate_policy.describe(pol)
        self.assertIn("declared: T-001 excluded on mcp: its capabilities document lacks release", text)
        self.assertIn("declared: R-002 answer_fields not required on mcp: the MCP never exposes the field", text)
        # T-001 still required on direct-local; its MCP lines are not judged
        self.assertIn({"scenario": "T-001", "arm": "direct-local", "repeats": [0]}, pol["required"])
        change = {("T-001", 0): {"outcome": "inconclusive"},
                  ("R-002", 0): {"checks": checks_for("R-002", 0, prop="not_exposed")}}
        v = gate.judge(make_run(self.run_dir(), self.m, mcp_lines(self.m, change=change), direct=("direct-local",)),
                       pol)
        self.assertTrue(v["passed"], gate.render(v))

    def test_declared_repeat_and_relation_entries(self):
        # a planned duplicate of one repeat (inconclusive by the scenario's design) and a relation the host cannot
        # produce, both declared on direct-local: neither is required there, everything else stays
        decl = self.declared(exclusions=[{"scenario": "D-007", "arms": ["direct-local"], "repeats": [2],
                                          "reason": "a planned duplicate"}],
                             relations_excluded=[{"relation": "REL-1", "arms": ["direct-local"],
                                                  "reason": "never testable on this host"}])
        pol = policy_of(self.m, arms=("mcp", "direct-local"), declared=decl)
        req = dict(((i["scenario"], i["arm"]), i["repeats"]) for i in pol["required"])
        self.assertEqual(req[("D-007", "direct-local")], [0, 1])
        self.assertEqual(req[("D-007", "mcp")], [0, 1, 2])
        self.assertEqual([(r["relation"], r["arms"]) for r in pol["relations"]], [("REL-1", ["mcp"])])
        self.assertIn({"relation": "REL-1", "arms": ["direct-local"], "rules": ["declared"],
                       "reason": "never testable on this host"}, pol["relations_excluded"])
        self.assertIn("D-007", gate.policy_ids(pol))
        text = gate_policy.describe(pol)
        self.assertIn("declared: D-007 r2 excluded on direct-local: a planned duplicate", text)
        self.assertIn("declared: REL-1 excluded on direct-local: never testable on this host", text)
        base = make_run(self.run_dir(), self.m, mcp_lines(self.m), direct=("direct-local",),
                        direct_change={("D-007", 2): {"outcome": "inconclusive"}})
        comp = sc.read_json(os.path.join(base, "compare.json"))
        rel = [r for r in comp["relations"] if r["id"] == "REL-1"][0]
        rel["arms"]["direct-local"] = {"outcome": "inconclusive", "detail": "no two profiles"}
        sc.write_json_atomic(os.path.join(base, "compare.json"), comp)
        v = gate.judge(base, pol)
        self.assertTrue(v["passed"], gate.render(v))
        v = gate.judge(make_run(self.run_dir("b"), self.m, mcp_lines(self.m), direct=("direct-local",),
                                direct_change={("D-007", 1): {"outcome": "inconclusive"}}), pol)
        self.assertEqual([(p["scenario"], p["arm"], p["repeat"]) for p in self.problems(v, "outcome")],
                         [("D-007", "direct-local", 1)])

    def test_declared_file_is_validated(self):
        for bad, why in ((self.declared(declared_version=2), "declared_version"),
                         (self.declared(exclusions=[{"scenario": "T-001", "arms": ["mcp"]}]), "reason"),
                         (self.declared(exclusions=[{"scenario": "T-999", "arms": ["mcp"], "reason": "x"}]), "T-999"),
                         (self.declared(not_required=[{"scenario": "R-002", "arms": ["mcp"], "check": "bogus",
                                                       "reason": "x"}]), "bogus"),
                         (self.declared(manifest={"content_sha256": "other"}), "content_sha256"),
                         (self.declared(exclusions=[{"scenario": "T-001", "arms": ["mcp"], "repeats": ["1"],
                                                     "reason": "x"}]), "repeats"),
                         (self.declared(relations_excluded=[{"relation": "REL-9", "arms": ["mcp"], "reason": "x"}]),
                          "REL-9")):
            with self.assertRaises(ValueError) as ctx:
                policy_of(self.m, declared=bad)
            self.assertIn(why, str(ctx.exception))

    def test_cli_declared_option(self):
        mpath, dpath, pol = (os.path.join(self.tmp, n) for n in ("m.json", "declared.json", "p.json"))
        sc.write_json_atomic(mpath, self.m)
        sc.write_json_atomic(dpath, self.declared())
        with contextlib.redirect_stdout(io.StringIO()) as out:
            self.assertEqual(gate_policy.main([mpath, "--by-design", "", "--arms", "mcp,direct-local",
                                               "--declared", dpath, "--out", pol]), 0)
        self.assertIn("declared: T-001 excluded on mcp", out.getvalue())
        self.assertEqual(gate.load_policy(pol)["declared"]["file"], dpath)


class PlanTest(GateCase):
    def test_plan_union_and_reconstruction(self):
        d = self.run_dir()
        os.makedirs(d)
        first = gate.plan_entry(self.m, None, None, ["T-001"], 0, ["mcp"], "paired", 0)
        self.assertTrue(gate.write_plan(d, first))
        self.assertFalse(gate.write_plan(d, dict(first, written_utc="later")))   # the same plan again: nothing
        second = gate.plan_entry(self.m, None, None, ["R-002"], 0, ["mcp"], "paired", 1)
        self.assertTrue(gate.write_plan(d, second))
        keys, source, plans = gate.run_plan(d, self.m)
        self.assertEqual(source, "plan.json")
        self.assertEqual(len(plans), 2)
        # a resume never narrows: the union of both selections
        self.assertEqual(keys, {("T-001", "mcp", 0), ("T-001", "mcp", 1), ("R-002", "mcp", 0), ("R-002", "mcp", 1)})
        os.remove(os.path.join(d, "plan.json"))
        sc.write_json_atomic(os.path.join(d, "orchestrator.json"), {"sessions": [{
            "areas": "", "limit": "", "skip_local": "1", "skip_remote": "0", "warm_every": "0", "mcp_cmd": "x",
            "ids": "D-007"}]})
        keys, source, _ = gate.run_plan(d, self.m)
        self.assertEqual(source, "orchestrator.json (reconstructed)")
        self.assertEqual(keys, {("D-007", a, r) for a in ("direct-remote", "mcp") for r in (0, 1, 2)})


class PlanSessionsTest(GateCase):
    def test_sessions_before_plan_json_join_the_union(self):
        # session 0 ran before plan.json existed (orchestrator.json only), session 1 wrote its plan entry
        d = self.run_dir()
        os.makedirs(d)
        sc.write_json_atomic(os.path.join(d, "orchestrator.json"), {"sessions": [
            {"areas": "", "limit": "", "skip_local": "1", "skip_remote": "1", "warm_every": "0", "mcp_cmd": "x",
             "ids": "D-007"},
            {"areas": "", "limit": "", "skip_local": "1", "skip_remote": "1", "warm_every": "0", "mcp_cmd": "x",
             "ids": "T-001"}]})
        gate.write_plan(d, gate.plan_entry(self.m, None, None, ["T-001"], 0, ["mcp"], "paired", 1))
        keys, source, _ = gate.run_plan(d, self.m)
        self.assertEqual(source, "plan.json and orchestrator.json (reconstructed)")
        self.assertEqual(sorted(set(k[0] for k in keys)), ["D-007", "T-001"])
        self.assertEqual(gate.resume_ids(d, self.m), ["T-001", "D-007"])
        # a later plan adding nothing new is not written; one adding an execution is
        self.assertFalse(gate.write_plan(d, gate.plan_entry(self.m, None, None, ["T-001"], 0, ["mcp"], "paired", 2)))
        self.assertTrue(gate.write_plan(d, gate.plan_entry(self.m, None, None, ["R-002"], 0, ["mcp"], "paired", 3)))


class JudgeTest(GateCase):
    def judge(self, d, **kw):
        return gate.judge(d, policy_of(self.m, **kw))

    def test_complete_mcp_only_run_passes(self):
        v = self.judge(make_run(self.run_dir(), self.m, mcp_lines(self.m)))
        self.assertTrue(v["passed"], gate.render(v))
        # R-003 refused by its named case, verified on the MCP side only (no direct arm in the plan)
        self.assertIn("1 by-design scenario(s) verified on the MCP side only", " ".join(v["notes"]))
        self.assertEqual(v["counts"]["judged"], 9)

    def test_three_arm_run_with_agreement_passes(self):
        d = make_run(self.run_dir(), self.m, mcp_lines(self.m), direct=True)
        v = gate.judge(d, policy_of(self.m, arms=("direct-local", "direct-remote", "mcp"), require_agreement=True))
        self.assertTrue(v["passed"], gate.render(v))

    def test_missing_mcp_execution_fails(self):
        v = self.judge(make_run(self.run_dir(), self.m, mcp_lines(self.m, skip=[("R-002", 1)])))
        self.assertFalse(v["passed"])
        self.assertEqual([(p["scenario"], p["arm"], p["repeat"]) for p in self.problems(v, "missing")],
                         [("R-002", "mcp", 1)])

    def test_scenario_absent_from_every_arm_fails(self):
        v = self.judge(make_run(self.run_dir(), self.m, mcp_lines(self.m, skip=[("T-001", 0), ("T-001", 1)]),
                                plan=False))
        self.assertEqual(len(self.problems(v, "missing")), 2)
        self.assertIn("orchestrator.json (reconstructed)", " ".join(v["notes"]))

    def test_any_inconclusive_fails(self):
        v = self.judge(make_run(self.run_dir(), self.m, mcp_lines(self.m, change={
            ("D-007", 2): {"outcome": "inconclusive", "inconclusive_reason": "traversal_resolve_deadline"}})))
        self.assertEqual([p["detail"] for p in self.problems(v, "outcome")],
                         ["inconclusive (traversal_resolve_deadline)"])

    def test_not_implemented_on_a_selected_scenario_fails(self):
        # a flipped configuration refusing a basic request: not an exclusion of the policy, so it blocks
        v = self.judge(make_run(self.run_dir(), self.m, mcp_lines(self.m, change={
            ("T-001", 1): {"outcome": "not_implemented_mcp", "service_by_design_case": "traverse.advanced_path_disabled"}})))
        self.assertIn("not_implemented_mcp", self.problems(v, "outcome")[0]["detail"])

    def test_excluded_scenarios_are_not_judged(self):
        v = self.judge(make_run(self.run_dir(), self.m, mcp_lines(self.m, change={
            ("T-004", 0): {"outcome": "fail"}})))
        self.assertTrue(v["passed"], gate.render(v))

    def test_undeclared_by_design_case_fails(self):
        v = self.judge(make_run(self.run_dir(), self.m, mcp_lines(self.m, change={
            ("R-002", 0): {"outcome": "service_by_design", "service_by_design_case": "case.refused", "http": 400,
                           "error": dict(REFUSED)}})))
        self.assertIn("which the scenario does not name", self.problems(v, "outcome")[0]["detail"])

    def test_by_design_not_verified_fails(self):
        # the MCP answer is not the case's service answer: neither side verifies the pair
        bad = dict(REFUSED, code="other_code")
        v = self.judge(make_run(self.run_dir(), self.m, mcp_lines(self.m, change={
            ("R-003", 0): {"error": bad}, ("R-003", 1): {"error": bad}})))
        self.assertFalse(v["passed"])
        self.assertTrue(self.problems(v, "disagreement") or self.problems(v, "outcome"), gate.render(v))

    def test_required_property_not_exposed_fails_unless_not_required(self):
        change = {("T-001", 1): {"checks": checks_for("T-001", 1, prop="not_exposed")}}
        d = make_run(self.run_dir(), self.m, mcp_lines(self.m, change=change))
        v = self.judge(d)
        self.assertEqual([(p["scenario"], p["repeat"], p["detail"]) for p in self.problems(v, "property")],
                         [("T-001", 1, "required property answer_fields not_exposed")])
        v = self.judge(d, not_required=[{"scenario": "T-001", "arms": ["mcp"], "check": "answer_fields",
                                         "reason": "not exposed via MCP"}])
        self.assertTrue(v["passed"], gate.render(v))

    def test_property_evaluated_on_a_later_repeat_only(self):
        # determinism exists from repeat 1 on: absent on repeat 0 is no problem, absent on every repeat is
        v = self.judge(make_run(self.run_dir(), self.m, mcp_lines(self.m)))
        self.assertEqual(self.problems(v, "property"), [])
        change = dict(((("D-007", r), {"checks": checks_for("T-001", r)}) for r in (1, 2)))
        v = self.judge(make_run(self.run_dir("b"), self.m, mcp_lines(self.m, change=change)))
        self.assertEqual([p["detail"] for p in self.problems(v, "property")],
                         ["the expected check deterministic has no result on a passing execution"] * 2)

    def test_generic_checks_are_not_trusted_from_the_outcome(self):
        # a pass line without its status result, and a refusal passing with no check at all, fail the gate
        no_status = [c for c in checks_for("T-001", 0) if c["check"] != "status"]
        v = self.judge(make_run(self.run_dir(), self.m, mcp_lines(self.m, change={
            ("T-001", 0): {"checks": no_status}})))
        self.assertEqual([(p["scenario"], p["repeat"], p["detail"]) for p in self.problems(v, "outcome")],
                         [("T-001", 0, "the expected check status has no result on a passing execution")])
        d = make_run(self.run_dir("b"), self.m, mcp_lines(self.m), direct=("direct-local",),
                     direct_change={("R-003", 0): {"checks": []}})
        v = gate.judge(d, policy_of(self.m, arms=("mcp", "direct-local")))
        self.assertIn(("R-003", "direct-local", "the expected check status has no result on a passing execution"),
                      [(p["scenario"], p["arm"], p["detail"]) for p in v["problems"]])
        # a generic check not passing on a passing line (status not exposed) blocks too
        bad = [dict(c, outcome="not_exposed") if c["check"] == "status" else c for c in checks_for("R-002", 1)]
        v = self.judge(make_run(self.run_dir("c"), self.m, mcp_lines(self.m, change={("R-002", 1): {"checks": bad}})))
        self.assertIn("status not_exposed", " ".join(p["detail"] for p in self.problems(v, "outcome")))

    def test_property_absent_on_one_repeat_fails(self):
        # only determinism may be absent, and only on the first repeat (nothing to compare)
        v = self.judge(make_run(self.run_dir(), self.m, mcp_lines(self.m, change={
            ("T-001", 0): {"checks": [c for c in checks_for("T-001", 0) if c["check"] != "answer_fields"]}})))
        self.assertEqual([(p["scenario"], p["repeat"], p["detail"]) for p in self.problems(v, "property")],
                         [("T-001", 0, "the expected check answer_fields has no result on a passing execution")])
        v = self.judge(make_run(self.run_dir("b"), self.m, mcp_lines(self.m, change={
            ("D-007", 1): {"checks": checks_for("D-007", 0)}})))
        self.assertEqual([(p["scenario"], p["repeat"]) for p in self.problems(v, "property")], [("D-007", 1)])

    def test_mcp_disagreement_fails(self):
        change = {("T-001", 0): {"checks": checks_for("T-001", 0, prop="fail"), "outcome": "fail"}}
        d = make_run(self.run_dir(), self.m, mcp_lines(self.m), direct=True, direct_change=change)
        v = self.judge(d)
        self.assertEqual([p["scenario"] for p in self.problems(v, "disagreement")], ["T-001"])

    def test_require_agreement_needs_every_property_compared(self):
        change = {("T-001", 0): {"checks": checks_for("T-001", 0, prop="not_exposed")}}
        d = make_run(self.run_dir(), self.m, mcp_lines(self.m, change=change), direct=True)
        v = gate.judge(d, policy_of(self.m, arms=("direct-local", "direct-remote", "mcp"), require_agreement=True,
                                    not_required=[{"scenario": "T-001", "arms": ["mcp"],
                                                   "check": "answer_fields", "reason": "x"}]))
        self.assertTrue(v["passed"], gate.render(v))   # not required: neither exposed nor compared
        v = gate.judge(d, policy_of(self.m, arms=("direct-local", "direct-remote", "mcp"), require_agreement=True))
        self.assertIn("required property answer_fields not compared with direct-local",
                      [p["detail"] for p in self.problems(v, "agreement")])

    def test_failed_relation_fails(self):
        d = make_run(self.run_dir(), self.m, mcp_lines(self.m))
        comp = sc.read_json(os.path.join(d, "compare.json"))
        rel = [r for r in comp["relations"] if r["id"] == "REL-1"][0]
        rel["arms"]["mcp"] = {"outcome": "not_evaluated", "detail": "no answer"}
        sc.write_json_atomic(os.path.join(d, "compare.json"), comp)
        v = self.judge(d)
        self.assertEqual([p["detail"][:30] for p in self.problems(v, "relation")], ["REL-1 on mcp: not_evaluated: n"])

    def test_identity_mismatches_fail(self):
        other = manifest_of(rels=[dict(relations()[0], a="R-002", b="T-001")])
        pol = policy_of(self.m)
        # the run's manifest: same scenarios, the relation swapped
        v = gate.judge(make_run(self.run_dir(), other, mcp_lines(other)), pol)
        self.assertIn("relations_sha256", " ".join(p["detail"] for p in self.problems(v, "identity")))
        # an arm that ran other relations than the run's manifest
        v = gate.judge(make_run(self.run_dir("b"), self.m, mcp_lines(self.m), arm_manifest=other), pol)
        self.assertIn("arms/mcp relations_sha256", " ".join(p["detail"] for p in self.problems(v, "identity")))
        # the by-design table the policy was derived from
        v = gate.judge(make_run(self.run_dir("c"), self.m, mcp_lines(self.m)),
                       dict(pol, by_design_table={"sha256": "another"}))
        self.assertIn("by-design table", " ".join(p["detail"] for p in self.problems(v, "identity")))

    def test_configuration_evidence(self):
        v = self.judge(make_run(self.run_dir(), self.m, mcp_lines(self.m), caps=None))
        self.assertIn("records no service_capabilities", self.problems(v, "configuration")[0]["detail"])
        v = self.judge(make_run(self.run_dir("b"), self.m, mcp_lines(self.m),
                                caps=dict(CAPS, document=None, error="refused")))
        self.assertIn("read of the live capabilities failed", self.problems(v, "configuration")[0]["detail"])
        v = self.judge(make_run(self.run_dir("c"), self.m, mcp_lines(self.m),
                                caps=dict(CAPS, mismatches=[{"path": "service_settings.advanced_enabled"}])))
        self.assertIn("live / table mismatch", self.problems(v, "configuration")[0]["detail"])
        ev = {"case.refused": {"service_settings.advanced_enabled": "table"}}
        v = self.judge(make_run(self.run_dir("d"), self.m, mcp_lines(self.m, change={
            ("R-003", 0): {"request_evidence": ev}})))
        self.assertEqual([(p["scenario"], p["repeat"]) for p in self.problems(v, "configuration")], [("R-003", 0)])
        ev = {"case.refused": {"service_settings.advanced_enabled": "live"}}
        v = self.judge(make_run(self.run_dir("e"), self.m, mcp_lines(self.m, change={
            ("R-003", 0): {"request_evidence": ev}})))
        self.assertTrue(v["passed"], gate.render(v))

    def test_nothing_passed_fails(self):
        rows = [r for r in mcp_lines(self.m) if r["scenario"] in ("R-003", "T-004", "T-005")]
        m = manifest_of([s for s in scenarios() if s["id"] in ("R-003", "T-004", "T-005")], [])
        v = gate.judge(make_run(self.run_dir(), m, rows), policy_of(m))
        self.assertIn("nothing_passed", [p["kind"] for p in v["problems"]])

    def test_cli_exit_status_and_verdict_file(self):
        d = make_run(self.run_dir(), self.m, mcp_lines(self.m))
        pol = os.path.join(self.tmp, "policy.json")
        sc.write_json_atomic(pol, policy_of(self.m))
        cmd = [sys.executable, "-I", os.path.join(SUITE, "gate.py")]
        p = subprocess.run(cmd + [d, pol], capture_output=True, text=True)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertTrue(p.stdout.startswith("GATE PASS"))
        self.assertTrue(sc.read_json(os.path.join(d, "gate.json"))["passed"])
        d2 = make_run(self.run_dir("b"), self.m, mcp_lines(self.m, skip=[("T-001", 1)]))
        p = subprocess.run(cmd + [d2, pol], capture_output=True, text=True)
        self.assertEqual(p.returncode, 1)
        self.assertIn("GATE FAIL", p.stdout)
        self.assertIn("T-001 mcp r1: planned, not recorded", p.stdout)
        sc.write_json_atomic(pol, {"policy_version": 99})
        self.assertEqual(subprocess.run(cmd + [d, pol], capture_output=True, text=True).returncode, 2)


class StaleGateTest(GateCase):
    """A gate verdict no longer describes a run that changed after it (a later resume without --gate): the report
    shows it as stale, never as GATE PASS."""

    def test_report_marks_a_stale_verdict(self):
        d = make_run(self.run_dir(), self.m, mcp_lines(self.m))
        verdict = gate.judge(d, policy_of(self.m))
        self.assertTrue(verdict["passed"])
        self.assertEqual(set(verdict["inputs"]), {"executions.jsonl", "compare.json", "manifest.json", "plan.json"})
        sc.write_json_atomic(os.path.join(d, "gate.json"), verdict)
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertIn("**GATE PASS**", report.build_report(d))
        # a later invocation changed the run: one MCP line fails now, merged and compared again
        rows, _ = sc.read_jsonl(os.path.join(d, "arms", "mcp", "executions.jsonl"))
        rows[0]["outcome"] = "fail"
        os.remove(os.path.join(d, "arms", "mcp", "executions.jsonl"))
        for r in rows:
            sc.append_jsonl(os.path.join(d, "arms", "mcp", "executions.jsonl"), r)
        merge_run.merge(d)
        sc.write_json_atomic(os.path.join(d, "compare.json"), cmp.compare(d, TABLE))
        with contextlib.redirect_stdout(io.StringIO()):
            text = report.build_report(d)
        self.assertNotIn("**GATE PASS**", text)
        self.assertIn("**GATE verdict STALE**", text)
        self.assertIn("executions.jsonl", text.split("# Deploy gate", 1)[1].split("\n# ", 1)[0])


class ReportPlanTest(GateCase):
    """report.py's verdict reconciles the run's plan: a planned execution no arm recorded is never a pass."""

    def verdict(self, d):
        with contextlib.redirect_stdout(io.StringIO()):
            text = report.build_report(d)
        return text, text.split("# Verdict", 1)[1].strip().splitlines()[0]

    def test_missing_mcp_execution_is_incomplete(self):
        d = make_run(self.run_dir(), self.m, mcp_lines(self.m, skip=[("R-002", 0), ("R-002", 1)]))
        _, v = self.verdict(d)
        self.assertIn("INCOMPLETE", v)
        self.assertIn("2 planned execution(s) not recorded (mcp 2; plan of plan.json)", v)

    def test_scenario_absent_everywhere_is_incomplete_without_plan_file(self):
        d = make_run(self.run_dir(), self.m, mcp_lines(self.m, skip=[("D-007", 2)]), plan=False)
        _, v = self.verdict(d)
        self.assertIn("1 planned execution(s) not recorded (mcp 1; plan of orchestrator.json (reconstructed))", v)

    def test_complete_selective_mcp_only_run_still_passes(self):
        d = make_run(self.run_dir(), self.m, mcp_lines(self.m))
        text, v = self.verdict(d)
        self.assertTrue(v.startswith("**PASS"), v)
        self.assertNotIn("# Deploy gate", text)

    def test_gate_section(self):
        d = make_run(self.run_dir(), self.m, mcp_lines(self.m, skip=[("T-001", 1)]))
        verdict = gate.judge(d, policy_of(self.m))
        sc.write_json_atomic(os.path.join(d, "gate.json"), verdict)
        text, _ = self.verdict(d)
        self.assertIn("# Deploy gate", text)
        self.assertIn("**GATE FAIL**", text)
        self.assertIn("missing 1", text)


if __name__ == "__main__":
    unittest.main()
