"""Offline tests of manifest v6 (the owner, 2026-10-15 / -16): v5 with every default-budget /pattern request at
110,000 ms instead of 60,000, on both arms, and nothing else; the deploy gate's scenarios and relations byte-identical
to v5's;
the MCP exposure baseline rebased to the budget-free request digest, its lists unchanged. The test that rebuilds v6
skips without the record cache (SUITE_TRUTH_CACHE or truth/cache)."""
import collections
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
import manifest_lib as ml  # noqa: E402
import suite_common as sc  # noqa: E402
from test_builder_v5 import built, have_cache  # noqa: E402

V5_CONTENT = "72181c4a275e6bb94a905e6d43eaad4e179a6b436dd5e00175e3936ea3aa1875"
V6_CONTENT = "49907dbb95ca0a373199d40cabe888b4c136f376295a1fc526a742eae8231cd1"
RELATIONS = "39660cfa0fb98cc47ac923fad3aa33419a84a3d92d662154cf2d4164752f4a86"
V6_SMOKE_CONTENT = "ee731e491d7e1168665acb439de6a705b66c0bf3c3cbacae7007d2611277f7dd"
# manifest_lib.PATTERN_TIME_BUDGET_MS of v6 (the owner, 2026-10-16: 110,000, not 120,000)
V6_BUDGET = 110000
# round 2 (the owner, 2026-10-16): the /pattern checks the MCP's answer cannot carry (areas.pattern_mcp_unexposed), by
# what they read: exactly the 27 scenarios whose MCP executions read them not_exposed in the v5 baseline 81c620ef,
# less the 3 whose field the MCP answer does carry per entry (MCP_PATTERN_PER_ENTRY: Q-pred-0011 / -0012's
# predicate block, P-edge-0009's load_ms; the checks read them there, checks.MCP_PATTERN_ENTRY_FIELDS). Two are
# scenario-level (their only property check).
ABSENCE_SCOPE = ["P-exact-%04d" % i for i in range(1, 9)] + ["P-scope-%04d" % i for i in range(1, 7)] + \
    ["D-P-exact-0001", "D-P-scope-0004"]
STRANDS = ["P-scope-%04d" % i for i in range(7, 13)] + ["N-usecase-0005"]
PATTERN_TEXT = ["P-edge-0004"]
SCENARIO_MARKED = {"P-scope-0004", "P-scope-0006"}
# round 2: the /pattern scenarios whose MCP cold call (r0) timed out or was time-limited in the v5 baseline, all 25 (the
# owner's decision 3; round2/verify-manifest/warm_targets.out), each with its warm-up: in 23 a later execution (the
# warm MCP repeat or a direct arm) completed (round2/pg/warm_up_targets.txt), in M-protein-0002 and P-heavy-0001
# (P-warm-0024 / -0025, round 2's verification, L4) none did
WARM_TARGETS = ["M-planted-%04d" % i for i in (7, 8, 9, 13, 14, 15, 16, 17, 18)] + \
    ["M-protein-%04d" % i for i in (1, 3, 4, 5, 6, 7)] + ["N-nearmiss-%04d" % i for i in (6, 13, 20, 27, 29, 30, 32)] + \
    ["N-usecase-0008", "M-protein-0002", "P-heavy-0001"]
N_WARM = 25
# round 2: the heavy /pattern scenarios the owner tags budget_stress (outside the gate, still run and reported)
STRESSED = {"M-protein-0001", "M-protein-0002", "M-protein-0003", "M-protein-0007", "P-heavy-0001", "N-nearmiss-0030"}
PER_ENTRY_UNMARKED = {"Q-pred-0011", "Q-pred-0012", "P-edge-0009"}
# the /pattern scenarios of v5 whose request names another budget than the default, or no JSON object at all
DELIBERATE = {"P-caps-0004": 300, "M-protein-0010": 30000}
NO_OBJECT = {"P-edge-0012", "C-protocol-0004", "C-protocol-0005", "C-protocol-0006"}


def read(name):
    return sc.read_json(os.path.join(SUITE, "manifest", name))


def by_id(m):
    return dict((s["id"], s) for s in m["scenarios"])


def leaf_diffs(a, b, path=""):
    """Every leaf difference between a and b as (path, old, new)."""
    out = []
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            if k not in a or k not in b:
                out.append((path + "." + k, a.get(k), b.get(k)))
            else:
                out.extend(leaf_diffs(a[k], b[k], path + "." + k))
    elif isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
        for i, (x, y) in enumerate(zip(a, b)):
            out.extend(leaf_diffs(x, y, "%s[%d]" % (path, i)))
    elif a != b:
        out.append((path, a, b))
    return out


def checks_is_property(name):
    import checks
    return checks.is_property(name)


class ManifestV6(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.v5, cls.v6 = read("suite_v5.json"), read("suite_v6.json")
        cls.s5, cls.s6 = by_id(cls.v5), by_id(cls.v6)

    def test_versions_and_digests(self):
        self.assertEqual(bm.MANIFEST_VERSION, 6)
        self.assertEqual(bm.SMOKE_TAG, "smoke-v5")   # v6 keeps v5's smoke selection and its tag (the gate's scenarios)
        self.assertEqual((self.v5["manifest_version"], self.v5["content_sha256"]), (5, V5_CONTENT))
        self.assertEqual((self.v6["manifest_version"], self.v6["content_sha256"]), (6, V6_CONTENT))
        self.assertEqual((self.v5["relations_sha256"], self.v6["relations_sha256"]), (RELATIONS, RELATIONS))
        self.assertEqual(sc.manifest_problems(self.v6), [])
        c5, c6 = json.loads(json.dumps(self.v5["counts"])), json.loads(json.dumps(self.v6["counts"]))
        self.assertEqual(c6.pop("mcp_not_implemented") - c5.pop("mcp_not_implemented"), len(SCENARIO_MARKED))
        n = len(WARM_TARGETS)
        for key, more in (("scenarios", n), ("executions_per_arm", n)):
            self.assertEqual(c6.pop(key) - c5.pop(key), more, key)
        self.assertEqual(c6["by_area"].pop("P") - c5["by_area"].pop("P"), n)
        for tag, more in (("budget_stress", n + len(STRESSED)), ("warm-up", n), ("machine-dependent", n)):
            self.assertEqual(c6["by_tag"].pop(tag) - c5["by_tag"].pop(tag), more, tag)
        self.assertEqual(c6, c5)
        self.assertEqual(self.v6["counts"]["by_tag"]["smoke-v5"], 14)

    def test_the_warm_ups_and_the_budget_stress_tags(self):
        """Round 2 (the owner, 2026-10-16): a warm-up right before each /pattern scenario whose MCP cold call timed
        out in v5 (its request, direct-local only: timing_class long; generic checks; budget_stress, as a cold run by
        design), so that the MCP arm's first call meets a warm host; budget_stress for the six heavy scenarios no
        budget below the ingress's cut makes decisive. The gate's executions do not move (keep_gate_warm_repeats)."""
        order = [s["id"] for s in self.v6["scenarios"]]
        warm = dict((s["warms"], s) for s in self.v6["scenarios"] if s["id"].startswith("P-warm-"))
        self.assertEqual(sorted(warm), sorted(WARM_TARGETS))
        self.assertEqual(sorted(w["id"] for w in warm.values()), ["P-warm-%04d" % i for i in range(1, N_WARM + 1)])
        for target, w in warm.items():
            t = self.s6[target]
            self.assertEqual(order.index(w["id"]) + 1, order.index(target), target)   # right before it
            self.assertEqual(w["request"], t["request"], target)                       # its own question
            self.assertEqual((w["route"], w["area"], w["timing_class"]), ("POST /pattern", "P", "long"), target)
            self.assertEqual(sc.mcp_repeats(order.index(w["id"]), w), [], target)       # the MCP arm sends nothing
            self.assertEqual(sorted(w["tags"]), ["budget_stress", "machine-dependent", "warm-up"], target)
            self.assertEqual([e["check"] for e in w["expect"] if checks_is_property(e["check"])], [], target)
            self.assertIn(areas.PATTERN_WARM_UPS[target][1], w["truth"]["warm_up"])
            self.assertEqual(w["truth"]["budget_stress"]["reason"], areas.PATTERN_WARM_STRESS)
        for sid in STRESSED:
            s = self.s6[sid]
            self.assertIn("budget_stress", s["tags"], sid)
            self.assertEqual(s["truth"]["budget_stress"]["reason"], areas.PATTERN_BUDGET_STRESS[sid], sid)
        self.assertEqual(sorted(s["id"] for s in self.v6["scenarios"] if "budget_stress" in s["tags"]
                                and not s["id"].startswith("P-warm-")),
                         sorted(STRESSED | set("T-longwalk-%04d" % i for i in range(1, 5))))
        # the deploy gate's executions: every scenario on its routes keeps v5's warm repeat
        i5 = dict((s["id"], i) for i, s in enumerate(self.v5["scenarios"]))
        i6 = dict((s["id"], i) for i, s in enumerate(self.v6["scenarios"]))
        for sid, s in self.s5.items():
            if s["route"] in bm.GATE_ROUTES:
                self.assertEqual(sc.repeat_count(i5[sid], s, 10), sc.repeat_count(i6[sid], self.s6[sid], 10), sid)
        import gate_policy
        self.assertEqual(bm.GATE_ROUTES, gate_policy.DEFAULT_ROUTES)
        # v5's scenarios keep their order outside area N, where keep_gate_warm_repeats moved some of its /pattern ones,
        # but for P-twin-0001: with 25 warm-ups R-twin-0005 (POST /resolve) would gain a warm repeat at index 690,
        # and no free scenario of area R precedes its block, so the nearest free one, P-twin-0001 (a /pattern twin of
        # P-edge-0004, another sequence than R-twin-0001's), goes right after the block (round 2's verification, L4)
        kept = [x for x in order if x in self.s5 and self.s5[x]["area"] != "N" and x != "P-twin-0001"]
        self.assertEqual(kept, [s["id"] for s in self.v5["scenarios"] if s["area"] != "N" and s["id"] != "P-twin-0001"])
        self.assertEqual(order[order.index("R-twin-0005") + 1], "P-twin-0001")
        self.assertEqual(order[order.index("R-twin-0001") - 1], "M-twin-0002")

    def test_keep_gate_warm_repeats(self):
        def s(sid, route="POST /pattern", area="N", **kw):
            return dict({"id": sid, "route": route, "area": area, "repeat": 1, "tags": [], "expect": []}, **kw)
        base = [s("A"), s("B"), s("C"), s("G1", "POST /resolve"), s("D"), s("E"), s("F"), s("G"), s("H"),
                s("I"), s("G2", "POST /resolve"), s("J")]
        ref = bm.gate_warm_status(base, 10)
        self.assertEqual(ref, {"G1": False, "G2": True})
        grown = [s("W", tags=["warm-up"], warms="A")] + base   # a warm-up before everything shifts both
        moved = bm.keep_gate_warm_repeats(grown, ref, 10)
        self.assertEqual(bm.gate_warm_problems(grown, ref, 10), [])
        self.assertEqual(bm.gate_warm_status(grown, 10), ref)
        self.assertEqual([x["id"] for x in grown if x["route"] == "POST /resolve"], ["G1", "G2"])
        self.assertTrue(moved)
        with self.assertRaises(SystemExit):   # its warm repeat would change, and nothing before it may move
            bm.keep_gate_warm_repeats([s("W", tags=["warm-up"], warms="X"), s("X"), s("G3", "POST /resolve",
                                                                                     area="R")], {"G3": False}, 2)

    def test_the_mcp_marks(self):
        """What the MCP's /pattern answer cannot carry is marked not implemented yet via MCP at property level (the
        owner's rule of 2026-10-15): entry_fields reading absence_scope, strands or the pattern text; scenario-level
        where that is the scenario's only property check (P-scope-0004, -0006). The predicate block and load_ms the
        MCP's answer carries per entry stay unmarked: the checks read them there (checks.MCP_PATTERN_ENTRY_FIELDS)."""
        marked = {}
        for s in self.v6["scenarios"]:
            for e in (s["mcp"].get("not_implemented_checks") or []):
                if s["route"] == "POST /pattern":
                    marked.setdefault(s["id"], []).append((e["check"], e["reason"]))
        self.assertEqual(sorted(marked), sorted(set(ABSENCE_SCOPE + STRANDS + PATTERN_TEXT) - SCENARIO_MARKED))
        for ids, what in ((ABSENCE_SCOPE, "the entry's absence_scope"), (STRANDS, "the entry's strands"),
                          (PATTERN_TEXT, "the entry's pattern text")):
            for sid in ids:
                if sid in SCENARIO_MARKED:
                    self.assertNotIn("tool", self.s6[sid]["mcp"], sid)
                    self.assertIn("not implemented yet via MCP: the MCP's /pattern answer does not carry %s" % what,
                                  self.s6[sid]["mcp"]["not_implemented"], sid)
                    continue
                self.assertEqual([c for c, _ in marked[sid]], ["entry_fields"], sid)
                self.assertIn(what, marked[sid][0][1], sid)
                self.assertFalse(bm.covers_every_property(self.s6[sid]), sid)
        for sid in PER_ENTRY_UNMARKED:
            self.assertNotIn("not_implemented_checks", self.s6[sid]["mcp"], sid)
            self.assertEqual(areas.pattern_mcp_unexposed(self.s6[sid]), {}, sid)
        self.assertEqual(areas.MCP_PATTERN_PER_ENTRY, ("predicate", "timing.load_ms"))

    def test_only_the_pattern_budgets_changed(self):
        """The diff v5 -> v6: the budgets, and the MCP marks of round 2 (test_the_mcp_marks), nothing else."""
        self.assertEqual(sorted(self.s5), sorted(set(self.s6) - set("P-warm-%04d" % i for i in range(1, N_WARM + 1))))
        kinds = collections.Counter()
        changed_routes = set()
        for sid, s in self.s5.items():
            s6 = json.loads(json.dumps(self.s6[sid]))
            if sid in STRESSED:   # test_the_warm_ups_and_the_budget_stress_tags
                s6["tags"] = sorted(set(s6["tags"]) - {"budget_stress"})
                s6["truth"].pop("budget_stress")
            if sid in SCENARIO_MARKED:
                s6["mcp"] = dict(s["mcp"], arguments=dict(s["mcp"]["arguments"], time_budget_ms=V6_BUDGET))
            elif s["route"] == "POST /pattern":
                s6["mcp"].pop("not_implemented_checks", None)
            diffs = leaf_diffs(s, s6)
            if diffs:
                changed_routes.add(s["route"])
            for path, old, new in diffs:
                kinds[(path, old, new)] += 1
        self.assertEqual(changed_routes, {"POST /pattern"})
        self.assertEqual(ml.PATTERN_TIME_BUDGET_MS, V6_BUDGET)
        self.assertEqual(kinds, {(".request.time_budget_ms", 60000, V6_BUDGET): 449,
                                 (".mcp.arguments.time_budget_ms", 60000, V6_BUDGET): 429})
        for sid, s in self.s6.items():
            if s["route"] != "POST /pattern":
                continue
            req = s.get("request")
            if sid in NO_OBJECT:
                self.assertFalse(isinstance(req, dict), sid)
                continue
            self.assertEqual(req["time_budget_ms"], DELIBERATE.get(sid, V6_BUDGET), sid)
            args = (s.get("mcp") or {}).get("arguments")
            if isinstance(args, dict) and "time_budget_ms" in args:
                self.assertEqual(args["time_budget_ms"], req["time_budget_ms"], sid)   # one question on every arm
        self.assertEqual(self.v6["relations"], self.v5["relations"])

    def test_gate_scenarios_and_relations_are_v5s(self):
        policy = sc.read_json(os.path.join(SUITE, "gate", "policy_v5_rt_mcp.json"))
        selected = sorted(set(r["scenario"] for r in policy["required"]))
        self.assertEqual(len(selected), 237)
        for sid in selected:
            self.assertNotEqual(self.s5[sid]["route"], "POST /pattern", sid)
            self.assertEqual(ml.canonical_json(self.s5[sid]), ml.canonical_json(self.s6[sid]), sid)
        rel5 = dict((r["id"], r) for r in self.v5["relations"])
        rel6 = dict((r["id"], r) for r in self.v6["relations"])
        required = [r["relation"] for r in policy["relations"]]
        self.assertEqual(len(required), 5)
        for rid in required:
            self.assertEqual(rel5[rid], rel6[rid], rid)
        # the v6 policy selects the same
        p6 = sc.read_json(os.path.join(SUITE, "gate", "policy_v6_rt_mcp.json"))
        self.assertEqual(sorted(set(r["scenario"] for r in p6["required"])), selected)
        self.assertEqual([r["relation"] for r in p6["relations"]], required)

    def test_no_gate_scenario_reaches_the_front_ingress_cut(self):
        """checks.ingress_cut_allowed (the front ingress's 504 at the cut reads inconclusive) covers a request whose
        stated budget reaches FRONT_INGRESS_CUT_MS - the tolerance: v6's 110,000 ms /pattern budget (inclusively, and
        its late window reaches past the cut measured a tolerance late) and nothing the deploy gate selects (every
        gated budget is at most 70,000 ms), so a gate run is judged as under v5."""
        import checks
        edge = checks.FRONT_INGRESS_CUT_MS - checks.FRONT_INGRESS_CUT_TOLERANCE_MS
        self.assertGreaterEqual(V6_BUDGET, edge)
        self.assertGreaterEqual(V6_BUDGET + max(checks._LATE_503_MIN_MS, V6_BUDGET),
                                checks.FRONT_INGRESS_CUT_MS + checks.FRONT_INGRESS_CUT_TOLERANCE_MS)
        a504 = checks.normalise_direct(504, b"<html>504 Gateway Time-out</html>")
        for name in ("policy_v5_rt_mcp.json", "policy_v6_rt_mcp.json"):
            policy = sc.read_json(os.path.join(SUITE, "gate", name))
            for sid in sorted(set(r["scenario"] for r in policy["required"])):
                s = self.s6[sid]
                bound = sc.stated_request_bound_s(s)
                self.assertTrue(bound is None or bound * 1000.0 < edge, (name, sid, bound))
                self.assertFalse(checks.ingress_cut_allowed(s, a504, {"wall_ms": 120000.0}), (name, sid))
        covered = [s["id"] for s in self.v6["scenarios"] if checks.ingress_cut_allowed(s, a504, {"wall_ms": 120000.0})]
        # the 449 scenarios at the 110,000 ms budget less the 22 that expect a refusal (a 400, answered at once), and
        # the 25 warm-ups of round 2 (their targets' requests; direct-local only, where the rule never applies)
        self.assertEqual(len(covered), 427 + N_WARM)
        self.assertEqual(sorted(x for x in covered if x.startswith("P-warm-")),
                         ["P-warm-%04d" % i for i in range(1, N_WARM + 1)])
        self.assertTrue(all(self.s6[x]["timing_class"] == "long" for x in covered if x.startswith("P-warm-")))
        self.assertTrue(all(self.s6[sid]["route"] == "POST /pattern" and
                            self.s6[sid]["request"]["time_budget_ms"] == V6_BUDGET for sid in covered))
        refused = [s["id"] for s in self.v6["scenarios"] if s["route"] == "POST /pattern" and isinstance(
            s.get("request"), dict) and s["request"].get("time_budget_ms") == V6_BUDGET and s["id"] not in covered]
        self.assertEqual(len(refused), 22)
        for sid in refused:
            st = [e for e in self.s6[sid]["expect"] if e["check"] == "status"]
            self.assertEqual(st[0]["args"].get("http"), 400, sid)
        self.assertEqual([s["id"] for s in self.v5["scenarios"] if checks.ingress_cut_allowed(s, a504, {"wall_ms": 120000.0})], [])

    def test_smoke_v6_is_smoke_v5_with_the_budget(self):
        sm5, sm6 = read("suite_smoke_v5.json"), read("suite_smoke_v6.json")
        self.assertEqual(sm6["content_sha256"], V6_SMOKE_CONTENT)
        self.assertEqual(sc.manifest_problems(sm6), [])
        self.assertEqual([s["id"] for s in sm5["scenarios"]], [s["id"] for s in sm6["scenarios"]])
        self.assertEqual(sm6["relations"], sm5["relations"])
        self.assertIn("smoke-v5 (full content_sha256 %s)" % V6_CONTENT, sm6["selection"])
        for a, b in zip(sm5["scenarios"], sm6["scenarios"]):
            diffs = leaf_diffs(a, b)
            if a["route"] == "POST /pattern":
                marked = [".mcp.not_implemented_checks"] if a["id"] in ABSENCE_SCOPE + STRANDS + PATTERN_TEXT else []
                self.assertEqual(sorted(d[0] for d in diffs),
                                 [".mcp.arguments.time_budget_ms"] + marked + [".request.time_budget_ms"], a["id"])
            else:
                self.assertEqual(diffs, [], a["id"])

    def test_exposure_baseline_rebased(self):
        """The baseline v6 reads is v4's with its digests recomputed by the budget-free rule from the v4 manifest:
        the exposes lists byte for byte v4's, and every /pattern scenario keeps its list in v6."""
        old_path = os.path.join(SUITE, "truth", "mcp_exposure", "v4-20261014T061053Z.json")
        new_path = os.path.join(SUITE, "truth", "mcp_exposure", "v4-20261014T061053Z-budget-free.json")
        old, new = sc.read_json(old_path), sc.read_json(new_path)
        v4 = read("suite_v4.json")
        self.assertEqual(bm.rebase_exposure(old, v4, old_path), new)
        self.assertEqual(sorted(new["scenarios"]), sorted(old["scenarios"]))
        for sid, e in old["scenarios"].items():
            self.assertEqual(new["scenarios"][sid]["exposes"], e["exposes"], sid)
        self.assertEqual(new["source"]["rebased_from"]["file"], "truth/mcp_exposure/v4-20261014T061053Z.json")
        self.assertEqual(new["source"]["manifest_content_sha256"], v4["content_sha256"])
        self.assertEqual(new["rule"], bm.EXPOSURE_RULE)
        with self.assertRaises(ValueError):
            bm.rebase_exposure(old, self.v6, old_path)   # not the manifest it was recorded against
        self.assertEqual(os.path.relpath(bm.DEFAULT_EXPOSURE, SUITE),
                         os.path.join("truth", "mcp_exposure", "v4-20261014T061053Z-budget-free.json"))
        self.assertEqual(self.v6["mcp_exposure_baseline"]["file"], "truth/mcp_exposure/v4-20261014T061053Z-budget-free.json")
        self.assertEqual(self.v6["mcp_exposure_baseline"]["sha256"], ml.sha256_json(new))
        # the lists of v6 are v5's, scenario for scenario, but where the MCP no longer runs a scenario (a
        # scenario-level mark of round 2)
        for sid, s in self.s5.items():
            want = None if sid in SCENARIO_MARKED else (s.get("mcp") or {}).get("exposes")
            self.assertEqual(want, (self.s6[sid].get("mcp") or {}).get("exposes"), sid)
        self.assertEqual(sum(1 for s in self.v6["scenarios"] if (s.get("mcp") or {}).get("exposes")),
                         self.v5["mcp_exposure_baseline"]["scenarios_with_exposes"] -
                         sum(1 for sid in SCENARIO_MARKED if (self.s5[sid].get("mcp") or {}).get("exposes")))

    def test_client_timeouts_exceed_every_stated_budget(self):
        """Every v6 request's client timeout exceeds its stated time budget by the margin on direct-local, and the
        front proxy's cap on direct-remote still exceeds the budget's deadline answer (~budget + 250 ms)."""
        import run_direct
        over = collections.Counter()
        for s in self.v6["scenarios"]:
            bound = sc.stated_request_bound_s(s)
            if bound is None:
                continue
            local = run_direct.timeout_for(s, "direct-local", sc.TIMING_CLASS_TIMEOUT_S)
            self.assertGreaterEqual(local, bound + sc.CLIENT_BUDGET_MARGIN_S, s["id"])
            if s.get("timing_class") != "long":
                remote = run_direct.timeout_for(s, "direct-remote", sc.TIMING_CLASS_TIMEOUT_S)
                self.assertGreater(remote, bound + 1.0, s["id"])
            over[(s["route"], bound, local)] += 1
        self.assertEqual(over[("POST /pattern", 110.0, 140.0)], 397)   # the fast ones; the 52 slow ones keep 600
        self.assertEqual(over[("POST /pattern", 110.0, 600.0)], 52)
        self.assertEqual(over[("POST /resolve", 60.0, 120.0)] + over[("POST /resolve", 60.0, 600.0)],
                         sum(1 for s in self.v6["scenarios"] if s["route"] == "POST /resolve"
                             and sc.stated_request_bound_s(s) == 60.0))

    def test_request_digest_leaves_the_pattern_budget_out(self):
        a = {"tool": "pattern_search", "arguments": {"patterns": [{"dna": "A" * 25}], "time_budget_ms": 60000}}
        b = {"tool": "pattern_search", "arguments": {"patterns": [{"dna": "A" * 25}], "time_budget_ms": 120000}}
        c = {"tool": "pattern_search", "arguments": {"patterns": [{"dna": "C" * 25}], "time_budget_ms": 120000}}
        self.assertEqual(bm.mcp_request_sha256(a), bm.mcp_request_sha256(b))
        self.assertNotEqual(bm.mcp_request_sha256(b), bm.mcp_request_sha256(c))
        fetch = {"tool": "traverse_fetch", "arguments": {"overrides": {"bounds": {"time_budget_ms": 60000}}}}
        other = {"tool": "traverse_fetch", "arguments": {"overrides": {"bounds": {"time_budget_ms": 120000}}}}
        self.assertNotEqual(bm.mcp_request_sha256(fetch), bm.mcp_request_sha256(other))   # a fetch's budget shapes the walk
        self.assertEqual(bm.BUDGET_FREE_DIGEST_TOOLS, ("traverse_resolve", "pattern_search"))

    @unittest.skipUnless(have_cache(), "no record cache")
    def test_frozen_v6_is_what_the_builder_makes(self):
        m = built()
        self.assertEqual(m["content_sha256"], V6_CONTENT)
        self.assertEqual(m["relations_sha256"], RELATIONS)
        self.assertEqual(bm.smoke_of(m)["content_sha256"], V6_SMOKE_CONTENT)
        self.assertEqual(ml.canonical_json(m["scenarios"]), ml.canonical_json(self.v6["scenarios"]))



class BuilderRules(unittest.TestCase):
    """The builder's validation rules of the harness round of 2026-10-17 (the owner's decisions on the
    semantic-duplication findings), held on hand-made scenarios and on every manifest under manifest/."""

    @staticmethod
    def manifests():
        d = os.path.join(SUITE, "manifest")
        return [(n, sc.read_json(os.path.join(d, n))) for n in sorted(os.listdir(d)) if n.endswith(".json")]

    def test_same_answer_ignore_paths_name_object_keys(self):
        # SD-011: args.ignore names object keys, dot-separated, looked through lists (checks._drop and compare's
        # relation_drop alike); a step that is a list index or '*' matches nothing in either reader, so the builder
        # refuses it, as it refuses a list that is no list of such paths
        scen = [{"id": "A", "route": "POST /pattern"}, {"id": "B", "route": "POST /pattern"}]
        rel = lambda ignore: {"id": "REL-x", "type": "same_answer", "a": "A", "b": "B",  # noqa: E731
                              "args": {"ignore": ignore}}
        self.assertEqual(bm.validate_relations([rel(["timing", "patterns.stats.load_ms"])], scen), [])
        for bad in (["p.0.t"], ["p.*.t"], ["*"], ["p.-1"], ["p..t"], [""], "p.t", [5], {"p": 1}):
            probs = bm.validate_relations([rel(bad)], scen)
            self.assertTrue(any("args.ignore" in p for p in probs), (bad, probs))
        # no manifest uses args.ignore (the rule changes no build)
        for name, m in self.manifests():
            for r in m.get("relations") or []:
                self.assertNotIn("ignore", r.get("args") or {}, (name, r.get("id")))


    def test_mcp_by_design_is_a_list_of_case_ids(self):
        # SD-030: a string mcp.by_design is refused as such (validate iterated it character by character and named
        # each letter as an unknown case); every manifest states a list
        table = set(["pattern.caps_zero"])
        base = {"id": "P-x", "area": "P", "title": "t", "route": "POST /pattern", "truth": {}, "timing_class": "fast",
                "repeat": 1, "tags": [], "expect": [{"check": "status", "args": {"http": 200}}],
                "request": {"patterns": [], "mode": "count", "max_contexts": 1, "time_budget_ms": 1000}}
        ok = dict(base, mcp={"tool": "pattern_search", "arguments": {}, "by_design": ["pattern.caps_zero"]})
        self.assertEqual([p for p in bm.validate([ok], table) if "by_design" in p], [])
        bad = dict(base, mcp={"tool": "pattern_search", "arguments": {}, "by_design": "pattern.caps_zero"})
        probs = [p for p in bm.validate([bad], table) if "by_design" in p]
        self.assertEqual(probs, ["P-x: mcp.by_design \"pattern.caps_zero\" is not a list of case ids"])
        for name, m in self.manifests():
            for s in m.get("scenarios") or []:
                self.assertIsInstance((s.get("mcp") or {}).get("by_design", []), list, (name, s.get("id")))


    def test_one_status_expectation_and_one_reader(self):
        # SD-094 (the decision brief's option 2): validate refuses a scenario with more than one status expectation
        # or an expectation whose args is no object, so that the readers taking the first, the last or every entry
        # agree; the builder side reads the allowed statuses through manifest_lib.allowed_statuses, and the M filter
        # keeps above a kind cap only a scenario whose status admits no 200 (it read [200, 400] as a refusal)
        base = {"id": "P-x", "area": "P", "title": "t", "route": "POST /pattern", "truth": {}, "timing_class": "fast",
                "repeat": 1, "tags": [], "mcp": {"tool": "pattern_search", "arguments": {}},
                "request": {"patterns": [{"dna": "A" * 25}], "mode": "count", "max_contexts": 1,
                            "time_budget_ms": 1000, "max_mismatches": 9}}
        st = lambda http: {"check": "status", "args": {"http": http}}  # noqa: E731
        self.assertEqual([p for p in bm.validate([dict(base, expect=[st(200)])], set()) if "status" in p], [])
        two = dict(base, expect=[st(200), st(400)])
        self.assertIn("P-x: 2 status expectations, not one", bm.validate([two], set()))
        listed = dict(base, expect=[{"check": "status", "args": [200]}])
        self.assertIn("P-x: the expectation status's args [200] is not an object", bm.validate([listed], set()))
        self.assertEqual([ml.allowed_statuses(e) for e in ([], [st(404)], [st([200, 400])], [{"check": "status"}])],
                         [[200], [404], [200, 400], [200]])

        class Ctx(object):
            skipped = []

            def skip(self, area, what, why):
                self.skipped.append(what)
        kept = areas._within_kind_caps(Ctx(), "M", [dict(base, id="M-%d" % i, expect=[st(h)])
                                                   for i, h in enumerate((200, [200, 400], 400, [400, 413]))])
        self.assertEqual([x["id"] for x in kept], ["M-2", "M-3"])
        self.assertEqual([areas._status_200_only(dict(base, expect=[st(h)])) for h in (200, [200], [200, 400], 400)],
                         [True, True, False, False])
        # every manifest: at most one status expectation per scenario, every args an object
        for name, m in self.manifests():
            for sc_ in m.get("scenarios") or []:
                exp = sc_.get("expect") or []
                self.assertLessEqual(len([e for e in exp if e.get("check") == "status"]), 1, (name, sc_["id"]))
                self.assertTrue(all(isinstance(e.get("args", {}), dict) for e in exp), (name, sc_["id"]))


if __name__ == "__main__":
    unittest.main()
