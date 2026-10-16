"""The strict deploy gate (gate.py, gate_policy.py) and the plan reconciliation of the report (GPT review finding 1,
verified by VERIFY-A; the owner's strict gate of 2026-10-14): a run fulfilling its declared plan passes, an
intentionally MCP-only one without a direct arm; a missing execution, a scenario absent from every arm, any
inconclusive, a not-implemented or undeclared by-design outcome, a required property not exposed, a disagreement, a
failed relation, a manifest identity mismatch and missing configuration evidence fail."""

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
import gate  # noqa: E402
import gate_policy  # noqa: E402
import merge_run  # noqa: E402
import reevaluate  # noqa: E402
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
# a read of the service's configuration as async's runner records it: the host's identity in host_document (the
# service's backends.<db>.identity), the record's sha256 suite_common.service_document_sha256
CAPS = {"document": {"service_settings": {"advanced_enabled": False}}, "host_document": {"index_fp": "fp"},
        "database": "test", "read_at": "t", "base": "b", "mismatches": []}
CAPS["sha256"] = sc.service_document_sha256(CAPS)


def rehashed(caps, **fields):
    """|caps| with |fields| replaced and its sha256 recomputed (another document, recorded correctly)."""
    out = dict(caps, **fields)
    out["sha256"] = sc.service_document_sha256(out)
    return out


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


def make_run(d, m, mcp_rows, direct=False, caps=CAPS, plan=True, arm_manifest=None, direct_change=None,
             snapshot=True):
    """A gated run directory as run_suite.sh leaves it: the manifest, the orchestrator's session, plan.json, the
    run's snapshot of its checks (|snapshot|), the arms (the MCP arm with its manifest copy, service_capabilities and
    the checks it judged with, stated as a runner without invocations does at run.json's top; each direct session
    with the sha256 of its checks; |direct| True both direct arms, or a tuple of them), merged and compared with
    TABLE."""
    os.makedirs(d, exist_ok=True)
    files = sc.take_snapshot(SUITE, d, "test", 0, now=T0 - 60)["files"] if snapshot else {}
    judged = {"checks_sha256": files.get("checks.py"), "suite_common_sha256": files.get("suite_common.py")}
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
    if snapshot:
        run_extra.update(judged, checks_dir=os.path.join(d, "scripts"))
    if caps is not None:
        run_extra["service_capabilities"] = caps
    write_arm(d, "mcp", mcp_rows, run_extra=run_extra, raw_bodies=bodies(mcp_rows))
    sc.write_json_atomic(os.path.join(d, "arms", "mcp", "manifest.json"), arm_manifest or m)
    if direct:
        for arm in direct:
            rows = direct_lines(m, arm, direct_change)
            sessions = [{"server_instance_start": "inst1", "server_instance_end": "inst1",
                         "started_utc": sc.iso_utc(T0), **(judged if snapshot else {})}]
            write_arm(d, arm, rows, raw_bodies=bodies(rows), run_extra={"sessions": sessions})
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
        # the reason fits each of them: C-caps-0001's service answers its own view, it refuses nothing
        caps = [e for e in pol["exclusions"] if e["scenario"] == "C-caps-0001"][0]
        self.assertNotIn("refuses", caps["reason"])
        self.assertIn("capabilities.traverse_capabilities_document", caps["reason"])
        self.assertIn("only against a direct execution", caps["reason"])
        # the MCP-only policy with the declared file (the shipped declared entries on the MCP arm)
        only = gate_policy.build_policy(m, manifest_file="manifest/suite_v4.json",
                                        declared=sc.read_json(os.path.join(SUITE, "gate", "declared_v4_rt.json")))
        self.assertEqual(only["counts"]["excluded"], 36)
        self.assertEqual(only["counts"]["required"], {"mcp": {"scenarios": 186, "executions": 378}})
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
        # REL-prefix_equals-0004 (twin half vs full twin, no stated stop) stays required on direct-local; -0002 is
        # declared with its own reason, -0001 / -0003 with V05's
        self.assertEqual([(r["relation"], r["arms"]) for r in two["relations"]],
                         [("REL-same_answer-0003", ["direct-local", "mcp"]),
                          ("REL-prefix_equals-0004", ["direct-local"])])
        declared_rel = dict((r["relation"], r["reason"]) for r in two["relations_excluded"]
                            if "declared" in (r.get("rules") or ()))
        self.assertEqual(sorted(declared_rel), ["REL-prefix_equals-0001", "REL-prefix_equals-0002",
                                                "REL-prefix_equals-0003"])
        self.assertIn("a stated stop answers the prefix", declared_rel["REL-prefix_equals-0001"])
        self.assertIn("a stated stop answers the prefix", declared_rel["REL-prefix_equals-0003"])
        self.assertNotIn("stated stop", declared_rel["REL-prefix_equals-0002"])
        self.assertIn("never completed on this host", declared_rel["REL-prefix_equals-0002"])
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

    def test_shipped_v5_policy_is_regenerated_exactly(self):
        # the frozen v5 (the owner's approval, 2026-10-15) and its deploy-gate policy, regenerated with the command
        # DECLARED-v5 states: the figures the owner approved, no property exemption (gate.py applies the marks)
        shipped = sc.read_json(os.path.join(SUITE, "gate", "policy_v5_rt_mcp.json"))
        m = sc.read_json(os.path.join(SUITE, "manifest", "suite_v5.json"))
        self.assertEqual((m["content_sha256"], m["relations_sha256"]),
                         ("72181c4a275e6bb94a905e6d43eaad4e179a6b436dd5e00175e3936ea3aa1875",
                          "39660cfa0fb98cc47ac923fad3aa33419a84a3d92d662154cf2d4164752f4a86"))
        smoke = sc.read_json(os.path.join(SUITE, "manifest", "suite_smoke_v5.json"))
        self.assertEqual(smoke["content_sha256"], "abc0596fef7977dc0fa7d51bf7b8d0194e6427ed726f0a9acf12605116546e8d")
        declared = sc.read_json(os.path.join(SUITE, "gate", "declared_v5_rt.json"))
        self.assertEqual(declared["manifest"]["content_sha256"], m["content_sha256"])
        fresh = gate_policy.build_policy(m, arms=("mcp", "direct-local"), manifest_file="manifest/suite_v5.json",
                                         name="v5-rt-mcp", declared_file="gate/declared_v5_rt.json",
                                         declared=declared)
        for p in (shipped, fresh):
            p.pop("generated_utc")
        self.assertEqual(shipped, fresh)
        c = shipped["counts"]
        self.assertEqual((c["selected"], c["excluded"], c["relations_required"], c["properties_not_required"]),
                         (241, 36, 5, 0))
        self.assertEqual(c["exclusions_by_rule"], {"advanced_path": 1, "budget_stress": 4,
                                                   "manifest_not_implemented": 24, "table_not_expressible": 27})
        self.assertEqual(c["required"], {"mcp": {"scenarios": 205, "executions": 416},
                                         "direct-local": {"scenarios": 237, "executions": 273}})
        self.assertEqual([r["relation"] for r in shipped["relations"]],
                         ["REL-same_answer-0003"] + ["REL-prefix_equals-%04d" % i for i in range(1, 5)])

    def test_shipped_v6_policy_is_regenerated_exactly_and_is_v5s(self):
        # v6 (the owner, 2026-10-15 / -16): v5 with every default-budget /pattern request at 110,000 ms; its policy is
        # generated as v5's was (0 declared entries, 0 property exemptions, the same by-design table) and must equal
        # v5's but for the pins (the manifest's file, content and version, the name, the declared file and its sha256,
        # generated_utc), so that the v5 acceptance carries over: the same selection, required executions, relations,
        # exclusions, properties and counts
        shipped = sc.read_json(os.path.join(SUITE, "gate", "policy_v6_rt_mcp.json"))
        m = sc.read_json(os.path.join(SUITE, "manifest", "suite_v6.json"))
        self.assertEqual((m["manifest_version"], m["content_sha256"], m["relations_sha256"]),
                         (6, "935eec18a7612404871d345938df1dfbf415bb6d959991c4f1cb254d43bf1e85",
                          "39660cfa0fb98cc47ac923fad3aa33419a84a3d92d662154cf2d4164752f4a86"))
        smoke = sc.read_json(os.path.join(SUITE, "manifest", "suite_smoke_v6.json"))
        self.assertEqual(smoke["content_sha256"], "ee731e491d7e1168665acb439de6a705b66c0bf3c3cbacae7007d2611277f7dd")
        declared = sc.read_json(os.path.join(SUITE, "gate", "declared_v6_rt.json"))
        self.assertEqual(declared["manifest"], {"file": "manifest/suite_v6.json", "content_sha256": m["content_sha256"]})
        self.assertEqual((declared["exclusions"], declared["not_required"], declared["relations_excluded"]), ([], [], []))
        fresh = gate_policy.build_policy(m, arms=("mcp", "direct-local"), manifest_file="manifest/suite_v6.json",
                                         name="v6-rt-mcp", declared_file="gate/declared_v6_rt.json",
                                         declared=declared)
        for p in (shipped, fresh):
            p.pop("generated_utc")
        self.assertEqual(shipped, fresh)
        v5 = sc.read_json(os.path.join(SUITE, "gate", "policy_v5_rt_mcp.json"))
        self.assertEqual(shipped["by_design_table"]["sha256"], v5["by_design_table"]["sha256"])
        self.assertEqual(shipped["by_design_table"]["sha256"][:8], "d4b8636c")

        def pinned(p):
            q = json.loads(json.dumps(p))
            q.pop("generated_utc", None)
            q.pop("name")
            q["manifest"] = dict((k, v) for k, v in q["manifest"].items()
                                 if k not in ("file", "content_sha256", "manifest_version"))
            q["declared"] = dict((k, v) for k, v in q["declared"].items() if k not in ("file", "sha256"))
            return q
        self.assertEqual(pinned(shipped), pinned(v5))
        self.assertEqual((shipped["name"], shipped["manifest"]["file"], shipped["manifest"]["manifest_version"],
                          shipped["declared"]["file"]),
                         ("v6-rt-mcp", "manifest/suite_v6.json", 6, "gate/declared_v6_rt.json"))
        self.assertEqual(shipped["manifest"]["index_fp"], v5["manifest"]["index_fp"])
        c = shipped["counts"]
        self.assertEqual((c["selected"], c["excluded"], c["relations_required"], c["properties_not_required"]),
                         (241, 36, 5, 0))
        self.assertEqual(c["required"], {"mcp": {"scenarios": 205, "executions": 416},
                                         "direct-local": {"scenarios": 237, "executions": 273}})
        self.assertNotIn("POST /pattern", c["by_route"])

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


# an answered case whose path removes a field: the duplicate's result (as traverse.duplicate_seeds explains on MCP)
MERGED = {"id": "case.answered_merge", "agreement": "by_design",
          "direct_expectation": {"http": 200, "code": None, "text": None},
          "service_answer": {"outcome": "answered", "code": None},
          "explains_checks": [{"check": "contains_source_context", "detail_matches": "1 results for 2 seeds"},
                              {"check": "answer_fields",
                               "detail_matches": r"field not in the answer: results\.1\.duplicate"}]}


class ExplainedNotExposedTest(GateCase):
    """An answered case that explains a field its path removes ({check, detail_matches} on a not_exposed reading,
    compare.item_explained): the MCP execution reads service_by_design with the item in E, and the gate accepts the
    verified pair as for any answered case (gate.judge's property loop skips an explained item, the pair is judged
    in the agreement step); without the item in E the property blocks as not_exposed."""

    @staticmethod
    def merged_checks(field="results.1.duplicate"):
        return [{"check": "status", "outcome": "pass", "detail": "ok"},
                {"check": "contains_source_context", "outcome": "fail", "detail": "1 results for 2 seeds"},
                {"check": "answer_fields", "outcome": "not_exposed", "detail": "field not in the answer: " + field}]

    def run_with(self, name, explained_items, checks_by_rep=None, names_only=False, direct_change=None):
        """The T-010 run whose MCP lines (r0, r1) record |explained_items| (and the names of those items; with
        |names_only| the names alone) on |checks_by_rep| (default: the merged reading on both)."""
        scen = [s for s in scenarios() if s["id"] in ("T-001", "R-002")]
        scen.append(scenario("T-010", "POST /traverse", by_design=[MERGED["id"]]))
        m = manifest_of(scen, [])
        change = {}
        for r in (0, 1):
            checks = (checks_by_rep or {}).get(r) or self.merged_checks()
            change[("T-010", r)] = {"outcome": "service_by_design", "checks": checks,
                                    "service_by_design_case": MERGED["id"],
                                    "explained_checks": [checks[i]["check"] for i in explained_items]}
            if not names_only:
                change[("T-010", r)]["explained_items"] = list(explained_items)
        global TABLE
        saved = TABLE
        TABLE = (dict(list(saved[0].items()) + [(MERGED["id"], MERGED)]), dict(saved[1], cases=len(saved[0]) + 1))
        try:
            d = make_run(self.run_dir(name), m, mcp_lines(m, change=change), direct=("direct-local",),
                         direct_change=direct_change)
            pol = policy_of(m, arms=("mcp", "direct-local"))
        finally:
            TABLE = saved
        return d, pol

    def test_the_rule(self):
        checks = [{"check": "status", "outcome": "pass", "detail": "ok"},
                  {"check": "contains_source_context", "outcome": "fail", "detail": "1 results for 2 seeds"},
                  {"check": "answer_fields", "outcome": "not_exposed",
                   "detail": "field not in the answer: results.1.duplicate"}]
        entries = cmp.explains(MERGED)
        rule = sc.outcome_from_checks
        self.assertEqual(cmp.explained_outcome(checks, entries, rule),
                         ("service_by_design", ["contains_source_context", "answer_fields"]))
        self.assertEqual(cmp.explained_items(checks, entries), [1, 2])
        # another field missing, or the field present with another value: not explained
        other = [checks[0], checks[1], dict(checks[2], detail="field not in the answer: results.0.walk")]
        self.assertEqual(cmp.explained_items(other, entries), [1])
        wrong = [checks[0], checks[1], dict(checks[2], outcome="fail", detail="results.1.duplicate = False, "
                                                                              "expected True")]
        self.assertEqual(cmp.explained_outcome(wrong, entries, rule), ("fail", ["contains_source_context"]))
        # any unexplained failure still fails; a name entry never explains a not_exposed reading
        more = checks + [{"check": "invariants", "outcome": "fail", "detail": "result 0 without outcome"}]
        self.assertEqual(cmp.explained_outcome(more, entries, rule)[0], "fail")
        self.assertFalse(cmp.item_explained("answer_fields", checks[2]))
        self.assertTrue(cmp.item_explained("contains_source_context", checks[1]))

    def test_one_definition_for_every_reader(self):
        # the explanation of a check item is suite_common's, the shared rule async's runner calls (shared_rule takes
        # it from the vendored suite_common): compare uses those functions, reevaluate records explained_items with
        # them, merge_run drops a recorded explained item of either explainable outcome, the gate skips one
        import importlib.util
        for name in ("entry_explains", "item_explained", "explained_outcome", "explained_items"):
            self.assertIs(getattr(cmp, name), getattr(sc, name), name)
        self.assertEqual(sc.EXPLAINABLE_OUTCOMES, ("fail", "not_exposed"))
        spec = importlib.util.spec_from_file_location("suite_common_as_vendored", os.path.join(SUITE, "suite_common.py"))
        vendored = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(vendored)
        checks = [{"check": "status", "outcome": "pass", "detail": "ok"},
                  {"check": "contains_source_context", "outcome": "fail", "detail": "1 results for 2 seeds"},
                  {"check": "answer_fields", "outcome": "not_exposed",
                   "detail": "field not in the answer: results.1.duplicate"}]
        entries = cmp.explains(MERGED)
        self.assertEqual(vendored.explained_outcome(checks, entries, vendored.outcome_from_checks),
                         ("service_by_design", ["contains_source_context", "answer_fields"]))
        self.assertEqual(vendored.explained_items(checks, entries), [1, 2])
        self.assertFalse(vendored.item_explained("answer_fields", checks[2]))
        # merge_run: a recorded line's explained items are taken out before the checks' rule, the not_exposed one too
        seen = []
        row = {"scenario": "T-010", "arm": "mcp", "repeat": 0, "outcome": "finding", "checks": checks,
               "explained_checks": ["contains_source_context", "answer_fields"], "explained_items": [1, 2]}
        merge_run.stricter_outcomes([row], {}, lambda c: seen.append([x["check"] for x in c]) or "pass")
        self.assertEqual(seen, [["status"]])

    def test_the_gate_accepts_the_verified_pair(self):
        d, pol = self.run_with("merged", (1, 2))
        v = gate.judge(d, pol)
        self.assertTrue(v["passed"], gate.render(v))
        ma = sc.read_json(os.path.join(d, "compare.json"))["scenarios"]["T-010"]["mcp_agreement"]
        self.assertEqual((ma["class"], ma["verified"], ma["explained_items"]),
                         ("service_by_design", True, {"0": [1, 2], "1": [1, 2]}))
        # a line that records only the failure in E (a runner without the shared rule): the table explains the
        # removed field whatever the line says, so no property problem; the line's other explanation is named
        d, pol = self.run_with("failure-only", (1,))
        v = gate.judge(d, pol)
        self.assertEqual([(p["kind"], p["scenario"], p["repeat"], p["detail"]) for p in v["problems"]],
                         [("agreement", "T-010", r, "the MCP line records the explained items [1], the table's are "
                                                    "[1, 2] (a runner's explanation is not trusted)") for r in (0, 1)])

    def test_the_gate_takes_the_explanation_from_the_table_not_the_line(self):
        # VERIFY-G stage 2: a line's own claim never explains a reading at the gate: answer_fields not_exposed on
        # another path, claimed explained by items or by names, or on r1 beside a legitimate r0, blocks as the
        # property it is; a direct-arm line carrying the case's fields is no MCP explanation at all
        other = self.merged_checks("results.0.labels")
        missing = ("property", "T-010", "mcp", "required property answer_fields not_exposed")
        for name, kw, want in (
                ("claimed-items", dict(checks_by_rep={0: other, 1: other}), [(0, missing), (1, missing)]),
                ("claimed-names", dict(checks_by_rep={0: other, 1: other}, names_only=True),
                 [(0, missing), (1, missing)]),
                ("r1-claimed", dict(checks_by_rep={1: other}), [(1, missing)])):
            d, pol = self.run_with(name, (1, 2), **kw)
            v = gate.judge(d, pol)
            self.assertFalse(v["passed"], name)
            got = [(p["repeat"], (p["kind"], p["scenario"], p["arm"], p["detail"])) for p in v["problems"]
                   if p["kind"] == "property"]
            self.assertEqual(got, want, gate.render(v))
            # the claim that is not the table's is named
            agreement = [p for p in v["problems"] if p["kind"] == "agreement"]
            self.assertEqual(sorted(p["repeat"] for p in agreement), [r for r, _ in want], gate.render(v))
            self.assertTrue(all("(a runner's explanation is not trusted)" in p["detail"] for p in agreement))
        direct = self.merged_checks()
        direct[1] = dict(direct[1], outcome="pass", detail="2 results")
        d, pol = self.run_with("direct-claim", (1, 2), direct_change={("T-010", 0): {
            "checks": direct, "service_by_design_case": MERGED["id"], "explained_checks": ["answer_fields"],
            "explained_items": [2]}})
        v = gate.judge(d, pol)
        self.assertIn(("property", "direct-local", 0, "required property answer_fields not_exposed"),
                      [(p["kind"], p["arm"], p["repeat"], p["detail"]) for p in v["problems"]])


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

    def test_outcome_outside_the_vocabulary_is_named(self):
        # the merge reads it as error; the gate's problem names the outcome the runner recorded
        v = self.judge(make_run(self.run_dir(), self.m, mcp_lines(self.m, change={
            ("D-007", 2): {"outcome": "skipped_by_runner"}})))
        self.assertEqual([p["detail"] for p in self.problems(v, "outcome")],
                         ["error (outcome 'skipped_by_runner' outside the vocabulary)"])

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

    def test_release_changed_during_the_run_fails(self):
        # a deploy during the run (or between an interrupt and its resume): a session's closing GET states another
        # release than the run's first; the merge warns and states mixed_releases, --strict-identity exits 3 and the
        # gate fails on identity
        d = make_run(self.run_dir(), self.m, mcp_lines(self.m), direct=True)
        pol = policy_of(self.m, arms=("direct-local", "direct-remote", "mcp"), require_agreement=True)
        self.assertTrue(gate.judge(d, pol)["passed"])
        path = os.path.join(d, "arms", "direct-local", "run.json")
        run = sc.read_json(path)
        run["sessions"].append({"started_utc": "2026-10-14T10:00:00Z", "ended_utc": "2026-10-14T10:05:00Z",
                                "release_start": run["release"], "release_end": "refseq97-00000000"})
        sc.write_json_atomic(path, run)
        merged = merge_run.merge(d)
        self.assertTrue(merged["mixed_releases"])
        self.assertEqual(merged["release_identity"]["releases"], ["refseq97-b7214d29", "refseq97-00000000"])
        self.assertTrue([w for w in merged["warnings"] if "release changed during the run" in w])
        err = io.StringIO()
        with contextlib.redirect_stderr(err), contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(merge_run.main([d, "--strict-identity"]), 3)
        self.assertIn("host identity: the host's release changed during the run: refseq97-b7214d29 (arms/direct-local "
                      "first start) -> refseq97-00000000 (arms/direct-local session 1 end", err.getvalue())
        v = gate.judge(d, pol)
        self.assertFalse(v["passed"])
        self.assertIn("release changed during the run", " ".join(p["detail"] for p in self.problems(v, "identity")))

    def test_configuration_evidence(self):
        v = self.judge(make_run(self.run_dir(), self.m, mcp_lines(self.m), caps=None))
        self.assertIn("configuration evidence missing", self.problems(v, "configuration")[0]["detail"])
        v = self.judge(make_run(self.run_dir("b"), self.m, mcp_lines(self.m),
                                caps=dict(CAPS, document=None, error="refused")))
        self.assertIn("configuration evidence missing", self.problems(v, "configuration")[0]["detail"])
        self.assertIn("refused", self.problems(v, "configuration")[0]["detail"])
        v = self.judge(make_run(self.run_dir("c"), self.m, mcp_lines(self.m),
                                caps=dict(CAPS, mismatches=[{"path": "service_settings.advanced_enabled"}])))
        self.assertIn("live / table mismatch", self.problems(v, "configuration")[0]["detail"])
        ev = {"case.refused": {"service_settings.advanced_enabled": "table"}}
        v = self.judge(make_run(self.run_dir("d"), self.m, mcp_lines(self.m, change={
            ("R-003", 0): {"request_evidence": ev}})))
        self.assertEqual([p["detail"] for p in self.problems(v, "configuration")],
                         ["the reference service_settings.advanced_enabled was read from the table, not the live "
                          "document, on 1 MCP execution(s) (case.refused: R-003 r0)"])
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


class ConfigurationTest(GateCase):
    """Configuration evidence through compare.gate_configuration_problems (lane G2's definition) and what the gate
    adds: every invocation read one document, the advanced_path exclusions' premise holds live, request_evidence is
    well formed, and the per-line table evidence is one problem per reference path."""

    def judge(self, d, **kw):
        return gate.judge(d, policy_of(self.m, **kw))

    def conf(self, v):
        return [p["detail"] for p in v["problems"] if p["kind"] == "configuration"]

    def set_mcp_run(self, d, **fields):
        path = os.path.join(d, "arms", "mcp", "run.json")
        run = sc.read_json(path)
        run.update(fields)
        sc.write_json_atomic(path, run)

    def test_a_document_without_service_settings(self):
        v = self.judge(make_run(self.run_dir(), self.m, mcp_lines(self.m), caps=dict(CAPS, document={})))
        self.assertIn("configuration evidence missing", " ".join(self.conf(v)))

    def test_every_invocation_reads_one_document(self):
        d = make_run(self.run_dir(), self.m, mcp_lines(self.m))
        self.set_mcp_run(d, invocations=[{"started_utc": "a", "service_capabilities": CAPS},
                                         {"started_utc": "b"}])
        self.assertIn("1 of 2 invocation(s) of the MCP runner recorded no read of the live capabilities document "
                      "(invocation 1 started b)", " ".join(self.conf(self.judge(d))))
        other = dict(CAPS, sha256="d", document={"service_settings": {"advanced_enabled": True}})
        self.set_mcp_run(d, invocations=[{"started_utc": "a", "service_capabilities": CAPS},
                                         {"started_utc": "b", "service_capabilities": other}])
        self.assertIn("2 different capabilities documents in one run", " ".join(self.conf(self.judge(d))))
        self.set_mcp_run(d, invocations=[{"started_utc": "a", "service_capabilities": CAPS},
                                         {"started_utc": "b", "service_capabilities": CAPS}])
        self.assertEqual(self.conf(self.judge(d)), [])

    def test_advanced_path_premise_checked_live(self):
        on = dict(CAPS, document={"service_settings": {"advanced_enabled": True}})
        v = self.judge(make_run(self.run_dir(), self.m, mcp_lines(self.m), caps=on))
        self.assertIn("the live service_settings.advanced_enabled is true", " ".join(self.conf(v)))
        # a policy for a deployment with the advanced path: required, and the live document must say so
        v = self.judge(make_run(self.run_dir("b"), self.m, mcp_lines(self.m)), advanced_path=True)
        self.assertIn("the live service_settings.advanced_enabled is false", " ".join(self.conf(v)))
        unstated = dict(CAPS, document={"service_settings": {}})
        v = self.judge(make_run(self.run_dir("c"), self.m, mcp_lines(self.m), caps=unstated))
        self.assertIn("does not state service_settings.advanced_enabled", " ".join(self.conf(v)))

    def test_malformed_request_evidence(self):
        for bad in ({"case.refused": "table"}, ["live"], {"case.refused": {"p": "maybe"}}):
            v = self.judge(make_run(self.run_dir(str(len(os.listdir(self.tmp)))), self.m, mcp_lines(
                self.m, change={("R-003", 1): {"request_evidence": bad}})))
            self.assertIn("malformed request_evidence", " ".join(self.conf(v)), bad)

    def test_table_evidence_is_one_problem_per_reference_path(self):
        ev = {"case.refused": {"service_settings.advanced_enabled": "table",
                               "presets.walk.budgets.max_extension_bp": "unknown"}}
        change = dict(((sid, r), {"request_evidence": ev}) for sid in ("R-003", "T-001") for r in (0, 1))
        v = self.judge(make_run(self.run_dir(), self.m, mcp_lines(self.m, change=change)))
        self.assertEqual(self.conf(v), [
            "the reference presets.walk.budgets.max_extension_bp was not found (unknown), live or in the table, on 4 "
            "MCP execution(s) (case.refused: T-001 r0, T-001 r1, R-003 r0, R-003 r1)",
            "the reference service_settings.advanced_enabled was read from the table, not the live document, on 4 MCP "
            "execution(s) (case.refused: T-001 r0, T-001 r1, R-003 r0, R-003 r1)"])

    def test_line_mismatches_are_read(self):
        v = self.judge(make_run(self.run_dir(), self.m, mcp_lines(self.m, change={("R-003", 0): {
            "request_evidence_mismatches": [{"namespace": "service", "ref": "x", "live": 1, "table": 2}]}})))
        self.assertIn("R-003 r0 recorded a live / table mismatch", " ".join(self.conf(v)))


class IdentityEvidenceTest(GateCase):
    """The identity a GATE PASS asserts needs its evidence (GPT review 2's B4, verified by VERIFY-G): with a snapshot
    every MCP invocation states both checks_sha256 and suite_common_sha256 (a checks directory stands in for
    neither); every policy arm observed the index_fp the policy pins, the MCP arm in its service_capabilities'
    host_document (the real runner writes no top-level index_fp); each capabilities record's identity is recomputed
    and a recorded sha256 its record does not reproduce blocks; a recorded set without a snapshot blocks."""

    def real_run(self, name="run", arms=("direct-local", "mcp"), invocations=2):
        """A snapshot-backed run in the real runner's shape: the MCP arm states its checks and its configuration per
        invocation (run.json invocations[], the latest at the top), no top-level index_fp."""
        direct = tuple(a for a in arms if a != "mcp")
        d = make_run(self.run_dir(name), self.m, mcp_lines(self.m), direct=direct)
        files = sc.read_json(os.path.join(d, "scripts", "SNAPSHOT.json"))["files"]
        path = os.path.join(d, "arms", "mcp", "run.json")
        run = sc.read_json(path)
        for k in ("checks_sha256", "suite_common_sha256", "checks_dir"):
            run.pop(k, None)
        run["invocations"] = [{"started_utc": sc.iso_utc(T0 + i), "checks_dir": os.path.join(d, "scripts"),
                               "checks_sha256": files["checks.py"], "suite_common_sha256": files["suite_common.py"],
                               "service_capabilities": dict(CAPS, read_at="t%d" % i)} for i in range(invocations)]
        sc.write_json_atomic(path, run)
        self.assertNotIn("index_fp", run)
        return d

    def edit(self, d, arm, fn):
        path = os.path.join(d, "arms", arm, "run.json")
        run = sc.read_json(path)
        fn(run)
        sc.write_json_atomic(path, run)
        merge_run.merge(d)
        sc.write_json_atomic(os.path.join(d, "compare.json"), cmp.compare(d, TABLE))
        return d

    def judge(self, d, arms=("direct-local", "mcp")):
        return gate.judge(d, policy_of(self.m, arms=arms))

    def blocked(self, v, kind, text):
        self.assertFalse(v["passed"])
        details = [p["detail"] for p in v["problems"] if p["kind"] == kind]
        self.assertTrue([x for x in details if text in x], (text, gate.render(v)))

    def test_the_real_shape_passes(self):
        for arms in (("direct-local", "mcp"), ("mcp",)):
            v = self.judge(self.real_run("-".join(arms), arms), arms)
            self.assertTrue(v["passed"], gate.render(v))

    def test_every_mcp_invocation_states_both_sha256(self):
        for i, drop in enumerate((("suite_common_sha256",), ("checks_sha256",),
                                  ("checks_sha256", "suite_common_sha256"))):
            d = self.edit(self.real_run(str(i)), "mcp",
                          lambda r: [r["invocations"][1].pop(k) for k in drop])
            self.blocked(self.judge(d), "identity", "arms/mcp: 1 of 2 invocation(s) do not state both checks_sha256 "
                                                    "and suite_common_sha256 (%s missing on 1" % drop[0])

    def test_every_policy_arm_observes_the_pinned_index(self):
        def no_host_index(r):
            for rec in [r["service_capabilities"]] + [x["service_capabilities"] for x in r["invocations"]]:
                rec["host_document"] = {}
                rec["sha256"] = sc.service_document_sha256(rec)
        for arms in (("direct-local", "mcp"), ("mcp",)):
            d = self.edit(self.real_run("mcp-" + "-".join(arms), arms), "mcp", no_host_index)
            self.blocked(self.judge(d, arms), "identity", "arms/mcp states no observed index_fp; the policy pins fp "
                                                          "(service_capabilities.host_document)")
        d = self.edit(self.real_run("direct"), "direct-local", lambda r: r.pop("index_fp"))
        self.blocked(self.judge(d), "identity", "arms/direct-local states no observed index_fp")
        # a direct arm outside the policy need not observe it, but one it observed must be the pin
        d = self.edit(self.real_run("other", ("direct-local", "mcp")), "direct-local", lambda r: r.pop("index_fp"))
        self.assertTrue(self.judge(d, ("mcp",))["passed"])
        d = self.edit(self.real_run("other2", ("direct-local", "mcp")), "direct-local",
                      lambda r: r["sessions"][0].update(index_fp_end="other"))
        self.blocked(self.judge(d, ("mcp",)), "identity", "arms/direct-local observed index_fp fp, other; the policy "
                                                          "pins fp")

    def test_a_conflicting_mcp_host_index_blocks(self):
        def other(r):
            for rec in [r["service_capabilities"]] + [x["service_capabilities"] for x in r["invocations"]]:
                rec["host_document"] = {"index_fp": "other"}
                rec["sha256"] = sc.service_document_sha256(rec)
        for arms in (("direct-local", "mcp"), ("mcp",)):
            d = self.edit(self.real_run("-".join(arms), arms), "mcp", other)
            self.blocked(self.judge(d, arms), "identity", "arms/mcp observed index_fp other; the policy pins fp")

    def test_a_stale_recorded_document_hash_blocks(self):
        def stale_invocation(r):
            rec = r["invocations"][1]["service_capabilities"]
            rec["document"] = {"service_settings": {"advanced_enabled": True}}   # its sha256 kept
        v = self.judge(self.edit(self.real_run("a"), "mcp", stale_invocation))
        self.blocked(v, "configuration", "invocation 1: the recorded capabilities sha256 %s is not its record's"
                     % CAPS["sha256"][:12])
        self.blocked(v, "configuration", "2 different capabilities documents in one run")

        def stale_top(r):
            r["service_capabilities"]["host_document"] = {"index_fp": "other"}   # its sha256 kept
        v = self.judge(self.edit(self.real_run("b"), "mcp", stale_top))
        self.blocked(v, "configuration", "run.json: the recorded capabilities sha256")
        self.blocked(v, "identity", "arms/mcp observed index_fp fp, other")
        # rehashed, the changed document is another deployment's
        v = self.judge(self.edit(self.real_run("c"), "mcp", lambda r: r["invocations"][1].update(
            service_capabilities=rehashed(CAPS, document={"service_settings": {"advanced_enabled": True}}))))
        self.blocked(v, "configuration", "2 different capabilities documents in one run")
        self.assertFalse([p for p in v["problems"] if "is not its record's" in p["detail"]])

    def test_a_hash_recorded_with_another_suite_common_is_checked_with_the_snapshots(self):
        # the run's snapshot's suite_common is not this one: the recorded sha256 is checked with the definition the
        # runner computed it with, the snapshot's own; a note only where that cannot be read
        def run_with(name, tail, top_sha=None, every_sha=None):
            d = self.real_run(name)
            with open(os.path.join(d, "scripts", "suite_common.py"), "a") as fh:
                fh.write(tail)
            snap = sc.read_json(os.path.join(d, "scripts", "SNAPSHOT.json"))
            snap["files"]["suite_common.py"] = sc.file_sha256(os.path.join(d, "scripts", "suite_common.py"))
            sc.write_json_atomic(os.path.join(d, "scripts", "SNAPSHOT.json"), snap)

            def restate(r):
                for x in r.get("invocations") or r.get("sessions"):
                    x["suite_common_sha256"] = snap["files"]["suite_common.py"]
                    if every_sha and "service_capabilities" in x:
                        x["service_capabilities"]["sha256"] = every_sha
                if "service_capabilities" in r:
                    r["service_capabilities"]["sha256"] = top_sha or every_sha or r["service_capabilities"]["sha256"]
            self.edit(d, "direct-local", restate)
            return self.judge(self.edit(d, "mcp", restate))

        # the same definition in another file: a stale hash blocks as it does with this suite_common
        v = run_with("same", "\n", top_sha="0" * 64)
        self.blocked(v, "configuration", "run.json: the recorded capabilities sha256 000000000000 is not its record's "
                                         "%s (the run's snapshot's suite_common.service_document_sha256)"
                     % CAPS["sha256"][:12])
        self.assertEqual(len(self.problems(v)), 1, gate.render(v))
        # another definition, every record computed with it: passes, nothing noted
        other = "\n\ndef service_document_sha256(rec):\n    return '1' * 64\n"
        v = run_with("other", other, every_sha="1" * 64)
        self.assertTrue(v["passed"], gate.render(v))
        self.assertFalse([x for x in v["notes"] if "recorded capabilities sha256" in x])
        # ... and a record that its definition does not reproduce blocks
        self.blocked(run_with("other-stale", other, top_sha="0" * 64, every_sha="1" * 64), "configuration",
                     "run.json: the recorded capabilities sha256 000000000000 is not its record's 111111111111 (the "
                     "run's snapshot's suite_common.service_document_sha256)")
        # a snapshot whose definition cannot be read: compared recomputed with this one, a note
        v = run_with("unread", "\ndel service_document_sha256\n", top_sha="0" * 64)
        self.assertTrue(v["passed"], gate.render(v))
        self.assertIn("run.json: the recorded capabilities sha256 000000000000 is not its record's %s "
                      "(suite_common.service_document_sha256); recorded with another suite_common than this one: "
                      "compared recomputed" % CAPS["sha256"][:12], v["notes"])

    def test_invocations_of_one_deployment_differ_in_read_times_only(self):
        def probes(r):
            for i, x in enumerate(r["invocations"]):
                rec = x["service_capabilities"]
                rec["document"] = dict(rec["document"], probed_at="2026-10-15T0%d:00:00Z" % i, cached=bool(i))
                rec["sha256"] = sc.service_document_sha256(rec)
        v = self.judge(self.edit(self.real_run(), "mcp", probes))
        self.assertTrue(v["passed"], gate.render(v))

    def test_a_recorded_set_without_a_snapshot_blocks(self):
        d = make_run(self.run_dir(), self.m, mcp_lines(self.m), snapshot=False)
        v = gate.judge(d, policy_of(self.m))
        self.assertEqual([p["detail"] for p in v["problems"]], [merge_run.UNPINNED])
        # its re-evaluated set is judged with the re-evaluation's checks (reevaluated.json), not on that ground
        for arm in ("mcp",):
            shutil.copy(os.path.join(d, "arms", arm, "executions.jsonl"),
                        os.path.join(d, "arms", arm, "executions.reevaluated.jsonl"))
        sc.write_json_atomic(os.path.join(d, "reevaluated.json"), {"checks": {
            "dir": SUITE, "checks_sha256": sc.file_sha256(os.path.join(SUITE, "checks.py")),
            "suite_common_sha256": sc.file_sha256(os.path.join(SUITE, "suite_common.py"))}})
        merge_run.merge(d, reevaluated=True)
        sc.write_json_atomic(os.path.join(d, "compare.reevaluated.json"), cmp.compare(d, TABLE, reevaluated=True))
        v = gate.judge(d, policy_of(self.m), reevaluated=True)
        self.assertNotIn(merge_run.UNPINNED, [p["detail"] for p in v["problems"]])
        self.assertEqual([p for p in v["problems"] if p["kind"] == "identity"], [])


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



class GatedHeadlineTest(GateCase):
    """The report's headline beside the gate (GPT review 2's B8, VERIFY-G's addition): a gated report leads with the
    gate's verdict and labels the ordinary one exploratory; a run, or an arm, with executions and no pass reads
    UNVERIFIED, never PASS; the gate's nothing_passed applies per policy arm."""

    def bold(self, d):
        with contextlib.redirect_stdout(io.StringIO()):
            text = report.build_report(d)
        return [x for x in text.splitlines() if x.startswith("**")], text

    def gated(self, d, pol):
        sc.write_json_atomic(os.path.join(d, "gate.json"), gate.judge(d, pol))
        return self.bold(d)

    def marks_only(self):
        # R-003 refused by design, T-004 not implemented via MCP: a legitimate MCP-only selection with no pass
        m = manifest_of([s for s in scenarios() if s["id"] in ("R-003", "T-004")], [])
        return m, make_run(self.run_dir(), m, mcp_lines(m))

    def test_marks_only_run_is_unverified(self):
        m, d = self.marks_only()
        lines, _ = self.bold(d)
        self.assertEqual(lines[0], "**UNVERIFIED (no execution passed: 1 not implemented via MCP, 2 service-side by "
                                   "design; nothing was shown to hold; no failure, no disagreement)**")
        lines, text = self.gated(d, policy_of(m))
        self.assertTrue(lines[0].startswith("**GATE FAIL**"), lines)
        self.assertIn("nothing\\_passed 1", lines[0])
        self.assertEqual(lines[1], "**Exploratory verdict (not the deploy decision, which is the gate's above): "
                                   "UNVERIFIED (no execution passed: 1 not implemented via MCP, 2 service-side by "
                                   "design; nothing was shown to hold; no failure, no disagreement)**")
        self.assertLess(text.index("# Deploy gate"), text.index("# Verdict"))

    def test_an_arm_that_passed_nothing(self):
        # direct-local passes T-004, which the MCP cannot express: the run reads UNVERIFIED on mcp, and the two-arm
        # gate fails nothing_passed for mcp (it passed on direct-local alone before)
        m = manifest_of([s for s in scenarios() if s["id"] == "T-004"], [])
        d = make_run(self.run_dir(), m, mcp_lines(m), direct=("direct-local",))
        lines, _ = self.bold(d)
        self.assertEqual(lines[0], "**UNVERIFIED (no execution passed on mcp: 1 not implemented via MCP; nothing was "
                                   "shown to hold there; no failure, no disagreement)**")
        v = gate.judge(d, policy_of(m, arms=("direct-local", "mcp")))
        self.assertEqual([(p["kind"], p["arm"], p["detail"]) for p in v["problems"]],
                         [("nothing_passed", "mcp", "no execution of arm mcp passed in the policy's selection")])
        self.assertTrue(gate.judge(d, policy_of(m, arms=("direct-local",)))["passed"])

    def test_mixed_and_clean_gated_runs_lead_with_the_gate(self):
        d = make_run(self.run_dir(), self.m, mcp_lines(self.m, change={
            ("D-007", 2): {"outcome": "inconclusive", "inconclusive_reason": "traversal_resolve_deadline"}}))
        lines, _ = self.gated(d, policy_of(self.m))
        self.assertTrue(lines[0].startswith("**GATE FAIL**"), lines)
        self.assertTrue(lines[1].startswith("**Exploratory verdict (not the deploy decision, which is the gate's "
                                            "above): PASS (no failed execution, no finding, no disagreement between "
                                            "the arms); 1 inconclusive execution(s)"), lines)
        d = make_run(self.run_dir("clean"), self.m, mcp_lines(self.m))
        lines, _ = self.gated(d, policy_of(self.m))
        self.assertTrue(lines[0].startswith("**GATE PASS**"), lines)
        self.assertTrue(lines[1].startswith("**Exploratory verdict (not the deploy decision, which is the gate's "
                                            "above): PASS"), lines)
        # a stale verdict leads too, and the ordinary verdict stays exploratory
        rows, _ = sc.read_jsonl(os.path.join(d, "executions.jsonl"))
        os.remove(os.path.join(d, "executions.jsonl"))
        for r in rows[1:]:
            sc.append_jsonl(os.path.join(d, "executions.jsonl"), r)
        lines, _ = self.bold(d)
        self.assertTrue(lines[0].startswith("**GATE verdict STALE**"), lines)
        self.assertTrue(lines[1].startswith("**Exploratory verdict"), lines)



class PropertyMarksTest(GateCase):
    """The v5 manifest's property-level marks of what the MCP cannot show (owner decision 2026-10-15):
    mcp.not_implemented_checks [{check, reason}] and mcp.by_design_checks [{check, case}]. A marked check the MCP
    reads not_exposed is set aside: the execution is neither inconclusive nor untested for it, it is no required-
    property gate problem and no agreement skip; the report lists it per scenario; any other reading of it is judged
    as usual. A manifest without the fields (v1-v4) is judged exactly as before."""

    NI = {"check": "answer_fields", "reason": "the MCP's records carry no such field"}
    BD = {"check": "answer_fields", "case": "case.answered_explains"}

    def manifest(self, marked=True, both=False):
        scen = [s for s in scenarios() if s["id"] in ("T-001", "R-002", "D-007")]
        if marked:
            scen[0]["mcp"]["not_implemented_checks"] = [dict(self.NI)] + (
                [{"check": "contains_source_context", "reason": "x"}] if both else [])
            scen[1]["mcp"]["by_design_checks"] = [dict(self.BD)]
        return manifest_of(scen, [])

    def marked_run(self, m, name="run", prop="not_exposed", direct=("direct-local",), change=None):
        """The MCP arm reads answer_fields |prop| on T-001 and R-002 (and contains_source_context too where
        |change| says so)."""
        hide = dict(((sid, r), {"checks": checks_for(sid, r, prop=prop)}) for sid in ("T-001", "R-002")
                    for r in (0, 1))
        hide.update(change or {})
        return make_run(self.run_dir(name), m, mcp_lines(m, change=hide), direct=direct)

    def rows(self, d):
        return dict((sc.execution_key(r), r) for r in sc.read_jsonl(os.path.join(d, "executions.jsonl"))[0])

    def report(self, d):
        with contextlib.redirect_stdout(io.StringIO()):
            return report.build_report(d)

    def test_helpers(self):
        m = self.manifest(both=True)
        s = m["scenarios"][0]
        self.assertEqual(sc.mcp_check_marks(s), {
            "answer_fields": {"mark": "not_implemented", "reason": self.NI["reason"]},
            "contains_source_context": {"mark": "not_implemented", "reason": "x"}})
        self.assertEqual(sc.mcp_check_marks(m["scenarios"][1]), {"answer_fields": {"mark": "by_design",
                                                                                   "case": self.BD["case"]}})
        checks = checks_for("T-001", 0, prop="not_exposed")
        self.assertEqual(sc.unmarked_checks(s, "mcp", checks), checks[:2])
        self.assertIs(sc.unmarked_checks(s, "direct-local", checks), checks)
        self.assertIs(sc.unmarked_checks(scenarios()[0], "mcp", checks), checks)   # no marks: unchanged
        failing = checks_for("T-001", 0, prop="fail")
        self.assertIs(sc.unmarked_checks(s, "mcp", failing), failing)   # a mark hides no result
        self.assertEqual(sc.mcp_check_marks({"mcp": {"not_implemented_checks": "x", "by_design_checks": [{}]}}), {})
        # marks covering every property check of the execution: nothing set aside, the marked items read
        # inconclusive (a scenario-level mark); the item's own `property` decides what a property check is
        whole = [dict(c, outcome="not_exposed") if c["check"] != "status" else c for c in checks]
        detail = sc.scenario_level_mark_detail(["answer_fields", "contains_source_context"])
        self.assertEqual(detail, "every property check is marked; this is a scenario-level mark (answer_fields, "
                                 "contains_source_context read not_exposed via MCP)")
        self.assertEqual(sc.every_property_marked(s, "mcp", whole), ["answer_fields", "contains_source_context"])
        self.assertEqual(sc.checks_marked_aside(s, "mcp", whole), {})
        self.assertEqual(sc.unmarked_checks(s, "mcp", whole), [whole[0]] + [
            dict(c, outcome="inconclusive", detail=detail) for c in whole[1:]])
        self.assertEqual([c["outcome"] for c in whole], ["pass", "not_exposed", "not_exposed"])   # not changed
        self.assertEqual(sc.every_property_marked(s, "direct-local", whole), [])
        self.assertEqual(sc.every_property_marked(self.manifest()["scenarios"][0], "mcp", whole), [])   # one marked
        flagged = [dict(whole[0], property=False), dict(whole[1], property=True), dict(whole[2], property=False)]
        self.assertEqual(sc.every_property_marked(s, "mcp", flagged), ["contains_source_context"])
        self.assertEqual(sc.every_property_marked(s, "mcp", whole[:1]), [])   # no property check: a refusal
        self.assertEqual(sc.GENERIC_CHECKS, cmp.load_checks_module(SUITE).GENERIC_CHECKS)

    def test_a_marked_check_makes_no_execution_inconclusive(self):
        # T-001's answer_fields reads not_exposed via MCP beside a passing contains_source_context: set aside, the
        # execution passes in the merge and the re-evaluation as it does without the mark
        mod = cmp.load_checks_module(SUITE)
        for marked in (False, True):
            m = self.manifest(marked)
            line = self.rows(self.marked_run(m, "one%s" % marked))[("T-001", "mcp", 0)]
            self.assertEqual(line["outcome"], "pass")
            self.assertEqual(reevaluate.outcome_of(mod, "mcp", m["scenarios"][0], {"http": 200, "body": {}}, line,
                                                   line["checks"], False)[0], "pass")

    def test_marks_covering_every_property_check_are_a_scenario_level_mark(self):
        # T-001's only property checks both read not_exposed via MCP and both are marked: nothing of its claim is
        # tested via MCP, which a scenario-level mark states; inconclusive with that reason in the merge (the runner
        # recorded pass), the re-evaluation and the comparison's fallback, and gate-blocking, never a pass on the
        # generic checks alone. Without the marks: inconclusive by the checks, as before.
        change = dict((("T-001", r), {"checks": [dict(c, outcome="not_exposed") if c["check"] != "status" else c
                                                 for c in checks_for("T-001", r)]}) for r in (0, 1))
        detail = sc.scenario_level_mark_detail(["answer_fields", "contains_source_context"])
        plain_dir = self.marked_run(self.manifest(False), "plain", change=change, direct=False)
        marked_dir = self.marked_run(self.manifest(both=True), "marked", change=change, direct=False)
        plain, marked = self.rows(plain_dir), self.rows(marked_dir)
        self.assertEqual(plain[("T-001", "mcp", 0)]["outcome"], "inconclusive")
        self.assertNotIn("inconclusive_reason", plain[("T-001", "mcp", 0)])
        for r in (0, 1):
            row = marked[("T-001", "mcp", r)]
            self.assertEqual((row["outcome"], row["outcome_runner"], row["inconclusive_reason"], row["verdict_detail"]),
                             ("inconclusive", "pass", sc.SCENARIO_LEVEL_MARK, detail))
        s = self.manifest(both=True)["scenarios"][0]
        line = dict(marked[("T-001", "mcp", 0)], outcome="pass")
        for k in ("outcome_runner", "inconclusive_reason", "verdict_detail"):
            line.pop(k)
        mod = cmp.load_checks_module(SUITE)
        info = {}
        self.assertEqual(reevaluate.outcome_of(mod, "mcp", s, {"http": 200, "body": {}}, line, line["checks"],
                                               False, info=info)[0], "inconclusive")
        self.assertEqual((info["inconclusive_reason"], info["verdict_detail"]), (sc.SCENARIO_LEVEL_MARK, detail))
        info = {}
        self.assertEqual(reevaluate.outcome_of(mod, "mcp", scenarios()[0], {"http": 200, "body": {}}, line,
                                               line["checks"], False, info=info)[0], "inconclusive")
        self.assertEqual(info["inconclusive_reason"], "checks")
        fallback = dict(line, outcome="inconclusive", arm="mcp")
        self.assertEqual(cmp.inconclusive_reason(fallback, s, marked_dir), sc.SCENARIO_LEVEL_MARK)
        self.assertEqual(cmp.inconclusive_reason(fallback, scenarios()[0], plain_dir), "checks")
        # the gate: an outcome problem naming the reason; the report counts it by reason and says it was not set aside
        m = self.manifest(both=True)
        v = gate.judge(marked_dir, policy_of(m))
        self.assertEqual([(p["kind"], p["scenario"], p["repeat"], p["detail"]) for p in v["problems"]],
                         [("outcome", "T-001", r, "inconclusive (%s)" % sc.SCENARIO_LEVEL_MARK) for r in (0, 1)])
        text = self.report(marked_dir)
        self.assertIn("2 %s" % report.md(sc.SCENARIO_LEVEL_MARK), text)
        self.assertIn("T-001 answer\\_fields (the MCP's records carry no such field; set aside on no MCP execution; "
                      "not set aside, %s: r0, r1)" % report.md(sc.SCENARIO_LEVEL_MARK), text)
        # a line that reads pass all the same (a merge without the checks' rule) is a property problem at the gate
        path = os.path.join(marked_dir, "executions.jsonl")
        rows = sc.read_jsonl(path)[0]
        for row in rows:
            if row["scenario"] == "T-001" and row["arm"] == "mcp":
                row["outcome"] = "pass"
                for k in ("outcome_runner", "inconclusive_reason", "verdict_detail"):
                    row.pop(k, None)
        with open(path, "w") as fh:
            fh.write("".join(json.dumps(row) + "\n" for row in rows))
        v = gate.judge(marked_dir, policy_of(m))
        self.assertFalse(v["passed"])
        self.assertEqual([(p["scenario"], p["repeat"], p["detail"]) for p in self.problems(v, "property")
                          if sc.SCENARIO_LEVEL_MARK in p["detail"]],
                         [("T-001", r, detail + " on a passing execution") for r in (0, 1)])

    def test_the_gate_reads_a_mark_as_no_required_property(self):
        for marked in (False, True):
            m = self.manifest(marked)
            d = self.marked_run(m, str(marked), direct=False)
            v = gate.judge(d, policy_of(m))
            props = [(p["scenario"], p["repeat"], p["detail"]) for p in self.problems(v, "property")]
            if marked:
                self.assertTrue(v["passed"], gate.render(v))
            else:
                self.assertEqual(props, [(sid, r, "required property answer_fields not_exposed")
                                         for sid in ("T-001", "R-002") for r in (0, 1)])
        # with require_agreement a marked check is no property left uncompared
        for marked in (False, True):
            m = self.manifest(marked)
            d = self.marked_run(m, "agree%s" % marked, direct=True)
            v = gate.judge(d, policy_of(m, arms=("direct-local", "direct-remote", "mcp"), require_agreement=True))
            agreement = [p["detail"] for p in self.problems(v, "agreement")]
            self.assertEqual(agreement, [] if marked else
                             ["required property answer_fields not compared with direct-local"] * 2)
        # a marked check the MCP reads fail is judged: an outcome problem, never set aside
        m = self.manifest()
        v = gate.judge(self.marked_run(m, "fail", prop="fail", direct=False), policy_of(m))
        self.assertFalse(v["passed"])
        self.assertIn(("T-001", "mcp", 0), [(p["scenario"], p["arm"], p["repeat"])
                                            for p in self.problems(v, "outcome")])

    def test_a_marked_check_without_a_result_is_a_missing_one(self):
        # the mark declares a not_exposed reading: a passing MCP line that has no result for the marked check (a
        # runner that dropped it) is the missing expected result it is without the mark, and with require_agreement
        # an uncompared property
        m = self.manifest()
        dropped = dict((("T-001", r), {"checks": checks_for("T-001", r)[:2], "outcome": "pass"}) for r in (0, 1))
        self.assertEqual(sc.checks_marked_aside(m["scenarios"][0], "mcp", dropped[("T-001", 0)]["checks"]), {})
        v = gate.judge(self.marked_run(m, "dropped", direct=False, change=dropped), policy_of(m))
        self.assertFalse(v["passed"])
        self.assertEqual([(p["scenario"], p["repeat"], p["detail"]) for p in self.problems(v, "property")],
                         [("T-001", r, "the expected check answer_fields has no result on a passing execution")
                          for r in (0, 1)])
        d = self.marked_run(m, "dropped-agree", direct=True, change=dropped)
        v = gate.judge(d, policy_of(m, arms=("direct-local", "direct-remote", "mcp"), require_agreement=True))
        self.assertIn("required property answer_fields not compared with direct-local",
                      [p["detail"] for p in self.problems(v, "agreement") if p["scenario"] == "T-001"])
        # the reading the mark declares stays set aside (R-002's by-design mark, read not_exposed, in the same runs)
        self.assertFalse([p for p in v["problems"] if p["scenario"] == "R-002"], gate.render(v))
        rows = self.rows(d)
        self.assertEqual(sc.checks_marked_aside(m["scenarios"][1], "mcp", rows[("R-002", "mcp", 0)]["checks"]),
                         {("answer_fields", 0): {"mark": "by_design", "case": self.BD["case"]}})

    def test_a_dropped_unmarked_property_is_no_scenario_level_mark(self):
        # T-001 marks answer_fields only; a runner that dropped the unmarked contains_source_context leaves the
        # marked answer_fields, read not_exposed, its only property item. The marks do not cover the property checks
        # the scenario expects: no scenario-level mark (the merge keeps the recorded pass, the comparison's fallback
        # names the checks), the marked reading stays set aside, and the gate names the unmarked check's missing
        # result
        m = self.manifest()
        s = m["scenarios"][0]
        dropped = dict((("T-001", r), {"checks": [c for c in checks_for("T-001", r, prop="not_exposed")
                                                  if c["check"] != "contains_source_context"], "outcome": "pass"})
                       for r in (0, 1))
        checks = dropped[("T-001", 0)]["checks"]
        self.assertEqual(sc.expected_property_checks(s), ["contains_source_context", "answer_fields"])
        self.assertEqual(sc.every_property_marked(s, "mcp", checks), [])
        self.assertEqual(sc.checks_marked_aside(s, "mcp", checks),
                         {("answer_fields", 0): {"mark": "not_implemented", "reason": self.NI["reason"]}})
        self.assertEqual(sc.unmarked_checks(s, "mcp", checks), checks[:1])
        d = self.marked_run(m, "dropped-unmarked", direct=False, change=dropped)
        rows = self.rows(d)
        for r in (0, 1):
            self.assertEqual(rows[("T-001", "mcp", r)]["outcome"], "pass")
            self.assertNotIn("inconclusive_reason", rows[("T-001", "mcp", r)])
        self.assertEqual(cmp.inconclusive_reason(dict(rows[("T-001", "mcp", 0)], outcome="inconclusive"), s, d),
                         "checks")
        v = gate.judge(d, policy_of(m))
        self.assertFalse(v["passed"])
        self.assertFalse([p for p in v["problems"] if sc.SCENARIO_LEVEL_MARK in p["detail"]], gate.render(v))
        self.assertEqual([(p["kind"], p["scenario"], p["repeat"], p["detail"]) for p in v["problems"]],
                         [("property", "T-001", r, "the expected check contains_source_context has no result on a "
                                                   "passing execution") for r in (0, 1)])
        self.assertIn("T-001 answer\\_fields (the MCP's records carry no such field; set aside on r0, r1)",
                      self.report(d))
        # the marks covering every property check the scenario expects stay a scenario-level mark: with a marked
        # check dropped (its missing result named beside it), and on a first repeat, where determinism has no result
        both = self.manifest(both=True)["scenarios"][0]
        self.assertEqual(sc.every_property_marked(both, "mcp", checks), ["answer_fields"])
        d7 = dict(scenarios()[6], mcp=dict(scenarios()[6]["mcp"], not_implemented_checks=[
            dict(self.NI), {"check": "contains_source_context", "reason": "x"}]))
        self.assertIn("deterministic", sc.expected_property_checks(d7))
        first = [dict(c, outcome="not_exposed") if c["check"] != "status" else c for c in checks_for("D-007", 0)]
        later = [dict(c, outcome="not_exposed") if c["check"] not in ("status", "deterministic") else c
                 for c in checks_for("D-007", 1)]
        self.assertEqual(sc.every_property_marked(d7, "mcp", first), ["answer_fields", "contains_source_context"])
        self.assertEqual(sc.every_property_marked(d7, "mcp", later), [])

    def test_compare_and_report_set_it_aside_and_list_it(self):
        texts = {}
        for marked in (False, True):
            d = self.marked_run(self.manifest(marked), str(marked))
            comp = sc.read_json(os.path.join(d, "compare.json"))
            ma = comp["scenarios"]["T-001"]["mcp_agreement"]
            texts[marked] = self.report(d)
            if marked:
                self.assertEqual((ma["class"], ma["checks_not_exposed"], ma["checks_marked"]), ("agree", 0, 2))
                self.assertEqual(ma["marked_checks"], [{"repeat": 0, "check": "answer_fields",
                                                        "mark": "not_implemented"},
                                                       {"repeat": 1, "check": "answer_fields",
                                                        "mark": "not_implemented"}])
                self.assertEqual(comp["summary"]["mcp_agreement"]["checks_marked"], 4)
            else:
                self.assertEqual((ma["class"], ma["checks_not_exposed"]), ("agree", 2))
                self.assertNotIn("checks_marked", ma)
                self.assertNotIn("checks_marked", comp["summary"]["mcp_agreement"])
        # the agreement's depth: partial without the marks, full with them
        self.assertIn("T-001 (answer\\_fields)", texts[False])
        depth_rows = [x for x in texts[True].splitlines() if x.startswith("| partial: ")]
        self.assertNotIn("T-001", depth_rows[0])
        self.assertIn("Not implemented yet via MCP (property): T-001 answer\\_fields (the MCP's records carry no "
                      "such field; set aside on r0, r1).", texts[True])
        self.assertIn("Service-side by design (property): R-002 answer\\_fields (case case.answered\\_explains; set "
                      "aside on r0, r1).", texts[True])
        self.assertIn("- Property-level marks of the manifest on the MCP arm (set aside where read not_exposed, "
                      "unless they cover every property check of the execution; listed under Agreement): 1 check(s) "
                      "of 1 scenario(s) not implemented yet via MCP; 1 check(s) of 1 scenario(s) service-side by "
                      "design.", texts[True])
        self.assertIn("4 set aside by a property-level mark of the manifest", texts[True])
        for x in ("(property)", "Property-level marks", "set aside by a property-level mark"):
            self.assertNotIn(x, texts[False])
        # a marked check read otherwise is judged and named so
        d = self.marked_run(self.manifest(), "pass", prop="pass")
        self.assertIn("T-001 answer\\_fields (the MCP's records carry no such field; set aside on no MCP execution; "
                      "judged as usual: r0 pass, r1 pass)", self.report(d))


if __name__ == "__main__":
    unittest.main()
