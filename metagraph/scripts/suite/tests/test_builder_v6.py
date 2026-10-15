"""Offline tests of manifest v6 (the owner, 2026-10-15): v5 with every default-budget /pattern request at 120,000 ms
instead of 60,000, on both arms, and nothing else; the deploy gate's scenarios and relations byte-identical to v5's;
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

import build_manifest as bm  # noqa: E402
import manifest_lib as ml  # noqa: E402
import suite_common as sc  # noqa: E402
from test_builder_v5 import built, have_cache  # noqa: E402

V5_CONTENT = "72181c4a275e6bb94a905e6d43eaad4e179a6b436dd5e00175e3936ea3aa1875"
V6_CONTENT = "e416fdc53241fd6e8da5114a0bba97c9a12253d949390495b60d232c14f8ec4e"
RELATIONS = "39660cfa0fb98cc47ac923fad3aa33419a84a3d92d662154cf2d4164752f4a86"
V6_SMOKE_CONTENT = "183295a285c6da48aa44e3a47b580d574ab8678376a5242ff36eff7be085587e"
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
        self.assertEqual(self.v6["counts"], self.v5["counts"])
        self.assertEqual(self.v6["counts"]["by_tag"]["smoke-v5"], 14)

    def test_only_the_pattern_budgets_changed(self):
        self.assertEqual([s["id"] for s in self.v5["scenarios"]], [s["id"] for s in self.v6["scenarios"]])
        kinds = collections.Counter()
        changed_routes = set()
        for sid, s in self.s5.items():
            diffs = leaf_diffs(s, self.s6[sid])
            if diffs:
                changed_routes.add(s["route"])
            for path, old, new in diffs:
                kinds[(path, old, new)] += 1
        self.assertEqual(changed_routes, {"POST /pattern"})
        self.assertEqual(kinds, {(".request.time_budget_ms", 60000, 120000): 449,
                                 (".mcp.arguments.time_budget_ms", 60000, 120000): 429})
        for sid, s in self.s6.items():
            if s["route"] != "POST /pattern":
                continue
            req = s.get("request")
            if sid in NO_OBJECT:
                self.assertFalse(isinstance(req, dict), sid)
                continue
            self.assertEqual(req["time_budget_ms"], DELIBERATE.get(sid, 120000), sid)
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
        stated budget reaches FRONT_INGRESS_CUT_MS - the tolerance: v6's 120,000 ms /pattern budget and nothing the
        deploy gate selects (every gated budget is at most 70,000 ms), so a gate run is judged as under v5."""
        import checks
        edge = checks.FRONT_INGRESS_CUT_MS - checks.FRONT_INGRESS_CUT_TOLERANCE_MS
        a504 = checks.normalise_direct(504, b"<html>504 Gateway Time-out</html>")
        for name in ("policy_v5_rt_mcp.json", "policy_v6_rt_mcp.json"):
            policy = sc.read_json(os.path.join(SUITE, "gate", name))
            for sid in sorted(set(r["scenario"] for r in policy["required"])):
                s = self.s6[sid]
                bound = sc.stated_request_bound_s(s)
                self.assertTrue(bound is None or bound * 1000.0 < edge, (name, sid, bound))
                self.assertFalse(checks.ingress_cut_allowed(s, a504, {"wall_ms": 120000.0}), (name, sid))
        covered = [s["id"] for s in self.v6["scenarios"] if checks.ingress_cut_allowed(s, a504, {"wall_ms": 120000.0})]
        # the 449 scenarios at the 120,000 ms budget less the 22 that expect a refusal (a 400, answered at once)
        self.assertEqual(len(covered), 427)
        self.assertTrue(all(self.s6[sid]["route"] == "POST /pattern" and
                            self.s6[sid]["request"]["time_budget_ms"] == 120000 for sid in covered))
        refused = [s["id"] for s in self.v6["scenarios"] if s["route"] == "POST /pattern" and isinstance(
            s.get("request"), dict) and s["request"].get("time_budget_ms") == 120000 and s["id"] not in covered]
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
                self.assertEqual(sorted(d[0] for d in diffs),
                                 [".mcp.arguments.time_budget_ms", ".request.time_budget_ms"], a["id"])
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
        # the lists of v6 are v5's, scenario for scenario
        for sid, s in self.s5.items():
            self.assertEqual((s.get("mcp") or {}).get("exposes"), (self.s6[sid].get("mcp") or {}).get("exposes"), sid)
        self.assertEqual(sum(1 for s in self.v6["scenarios"] if (s.get("mcp") or {}).get("exposes")),
                         self.v5["mcp_exposure_baseline"]["scenarios_with_exposes"])

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
        self.assertEqual(over[("POST /pattern", 120.0, 150.0)], 397)   # the fast ones; the 52 slow ones keep 600
        self.assertEqual(over[("POST /pattern", 120.0, 600.0)], 52)
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


if __name__ == "__main__":
    unittest.main()
