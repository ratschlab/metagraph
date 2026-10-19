"""The service's own behaviour in area G checked positively (graphlet_service_checks.py; DESIGN-v8 §5, D-R6) and the
declared per-way differences verified as by-design pairs against the table (D-R12): walk_display, label_evidence and
qualified_cause on cut-list graphlets, hint above 64 labels per node, in_service false plus download over 8 MiB,
max_bytes 262,145 refused on W3 and answered locally, a wrong token forbidden, a malformed handle bad_handle,
library_commit; the area G case schema of the by-design table (compare.graphlet_case_problems)."""

import copy
import json
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import compare as cmp  # noqa: E402
import compare_graphlets as cg  # noqa: E402
import graphlet_fixtures as gf  # noqa: E402
import graphlet_precheck as gp  # noqa: E402
import graphlet_service_checks as gs  # noqa: E402

SID, SID_B = "G-merge-0001", "G-cut-0002"
UUID = "0d5b9b1a-1111-2222-3333-444455556666"
UUID_B = "9e5b9b1a-aaaa-bbbb-cccc-444455556666"
FOUR = ["W1", "W2", "W3", "W4"]

CAP_CASE = {"id": "graphlet.max_bytes_above_service_cap", "areas": ["G"], "agreement": "by_design",
            "scenario_match": "a graphlet_* call with max_bytes above 262,144",
            "direct_expectation": {"behaviour": "the library's ceiling is 16 MiB"},
            "service_answer": {"tool": "graphlet_*", "outcome": "refused_at_submit", "code": "bad_argument",
                               "field": "max_bytes"},
            "graphlet": {"tools": ["graphlet_walks", "graphlet_summary"],
                         "args": [{"path": "max_bytes", "gt": 262144}],
                         "library": {"answered": True},
                         "service": {"error": {"code": "bad_argument"}}},
            "source": ["test"]}


def compare_answer(handle_a="g_0000000000aa", handle_b="g_0000000000bb", counts=None, comparable="qualified"):
    return {"a": handle_a, "b": handle_b, "mode": "claims", "comparable": comparable, "reason": "lower bound",
            "equal": None, "depth_used": {"left": 100}, "counts": counts or {"only_in_a": 0, "only_in_b": 0,
                                                                             "differ": 0},
            "notes": [], "support": ["kmer", "kmer"], "scopes": [{}, {}], "total": 0, "rows": []}


class ServiceBase(unittest.TestCase):
    def setUp(self):
        gf.install_graphlet_doc_stub()
        gf.install_oracle_stub()
        self.run_ = None

    def tearDown(self):
        if self.run_:
            self.run_.cleanup()

    def make(self, scenarios, table=None):
        self.run_ = gf.GRun(scenarios=scenarios)
        for s in scenarios:
            self.run_.complete_scenario(sid=s["id"], ops=[], mcp_runs=(s is scenarios[0]))
        if table is not None:
            data = {"version": 1, "service_commit": gf.SERVICE, "cases": table}
            with open(self.run_.path("scripts", "service_by_design.json"), "w") as fh:
                json.dump(data, fh)
        return self.run_

    def put(self, sid, op, answers, errors=None, sf=None, ways=FOUR):
        """Records op |op| of |sid| on |ways|: answers {way: answer} (a default for the others), errors {way: error},
        sf: W3's service fields."""
        mod, text = gf.install_graphlet_doc_stub(), self.run_.fixture.doc(0)
        doc = {"body_digest": mod.body_digest(text), "envelope_digest": mod.envelope_digest(text)}
        for way in ways:
            ans = copy.deepcopy(answers.get(way, answers.get("*")))
            self.run_.record(way, sid, op["id"], op["tool"], op["args"], answer=ans, error=(errors or {}).get(way),
                             doc=doc, service_fields=sf if way == "W3" else None)

    def result(self):
        res, _ = cg.compare_run(self.run_.dir, extra=cg.STAGES)
        return res

    def check(self, res, name):
        return next(c for c in res["service_checks"] if c["check"] == name)


class WalkDisplay(ServiceBase):
    OP = {"id": "walks", "tool": "graphlet_walks", "args": {"arm": "left", "rank": "id", "n": 10, "max_bytes": 2048},
          "ways": FOUR}

    def test_pass_fail_and_not_exercised(self):
        self.make([gf.scenario(ops=[self.OP])])
        lib = gf.walks_answer()
        self.put(SID, self.OP, {"*": lib, "W3": dict(lib, handle=UUID)}, sf={"walk_display": "first_parent"})
        res = self.result()
        c = self.check(res, "walk_display")
        self.assertEqual((c["outcome"], c["calls"]), ("pass", 1))
        self.assertEqual(self.check(res, "hint")["outcome"], "n/a")   # no compare call: its precondition (VC2)
        self.assertFalse(self.check(res, "hint")["exercised"])
        self.assertIn("n/a: no W3 call of the run of graphlet_compare answered", self.check(res, "hint")["detail"])
        # the runner left it in the answer: read as a service field all the same
        self.put(SID, self.OP, {"*": lib, "W3": dict(lib, handle=UUID, walk_display="last_parent")}, sf={})
        c = self.check(self.result(), "walk_display")
        self.assertEqual(c["outcome"], "fail")
        self.assertIn('walk_display "last_parent"', c["detail"])
        self.put(SID, self.OP, {"*": lib, "W3": dict(lib, handle=UUID)}, sf={})
        self.assertIn("walk_display (absent)", self.check(self.result(), "walk_display")["detail"])


class CompareFields(ServiceBase):
    def scenarios(self, cut_b=True, labels_b=70):
        op = {"id": "cmp", "tool": "graphlet_compare",
              "args": {"b": {"graphlet": SID_B}, "mode": "claims", "max_bytes": 262144}, "ways": FOUR}
        a = gf.scenario(sid=SID, ops=[op], facts={"cut_lists": False, "max_labels_per_node": 10})
        b = gf.scenario(sid=SID_B, ops=[], facts={"cut_lists": cut_b, "max_labels_per_node": labels_b})
        b["request"]["strategy"]["labels"]["max_labels_per_node"] = 32
        return [a, b], op

    def run_with(self, sf, counts=None, cut_b=True, labels_b=70):
        scens, op = self.scenarios(cut_b, labels_b)
        self.make(scens)
        lib = compare_answer(counts=counts)
        self.put(SID, op, {"*": lib, "W3": compare_answer(UUID, UUID_B, counts=counts)}, sf=sf)
        return self.result()

    def test_a_cut_list_graphlet_with_every_field(self):
        sf = {"label_evidence": "lower_bound", "qualified_cause": "caps 64 and 32: %s" % gs.NO_DIFFERENCE,
              "hint": "re-fetch both with overrides.labels.max_labels_per_node >= the label count"}
        res = self.run_with(sf)
        for name in ("label_evidence", "qualified_cause", "hint"):
            self.assertEqual(self.check(res, name)["outcome"], "pass", self.check(res, name))
        self.assertEqual(res["scenarios"][SID]["ops"]["cmp"]["outcome"], "agree")   # the fields are stripped

    def test_label_evidence_missing_or_out_of_place(self):
        res = self.run_with({"qualified_cause": "x"})
        self.assertEqual(self.check(res, "label_evidence")["outcome"], "fail")
        self.assertIn("qualified_cause present with label_evidence absent",
                      self.check(res, "qualified_cause")["detail"])
        self.run_.cleanup()
        res = self.run_with({"label_evidence": "lower_bound", "qualified_cause": "caps 64 and 32: %s" %
                             gs.NO_DIFFERENCE}, cut_b=False, labels_b=10)
        self.assertIn("no cut list, yet label_evidence", self.check(res, "label_evidence")["detail"])

    def test_qualified_cause_against_the_library(self):
        res = self.run_with({"label_evidence": "lower_bound", "qualified_cause": "caps 64 and 32: %s" %
                             gs.NO_DIFFERENCE, "hint": "h"}, counts={"only_in_a": 2, "only_in_b": 0, "differ": 0})
        c = self.check(res, "qualified_cause")
        self.assertEqual(c["outcome"], "fail")
        self.assertIn("the library counts", c["detail"])
        self.run_.cleanup()
        res = self.run_with({"label_evidence": "lower_bound", "qualified_cause": "cap 64: %s" % gs.NO_DIFFERENCE,
                             "hint": "h"})
        self.assertIn("does not state the cap 32 of %s" % SID_B, self.check(res, "qualified_cause")["detail"])

    def test_hint_iff_more_than_64_labels_per_node(self):
        res = self.run_with({"label_evidence": "lower_bound", "qualified_cause": "caps 64 and 32: %s" %
                             gs.NO_DIFFERENCE})
        self.assertIn("70 labels per node at most, hint absent", self.check(res, "hint")["detail"])
        self.run_.cleanup()
        res = self.run_with({"label_evidence": "lower_bound", "qualified_cause": "caps 64 and 32: %s" %
                             gs.NO_DIFFERENCE, "hint": "h"}, labels_b=64)
        self.assertIn("64 labels per node at most, hint present", self.check(res, "hint")["detail"])


class InService(ServiceBase):
    OP = {"id": "summary", "tool": "graphlet_summary", "args": {"detail": "summary", "max_bytes": 2048},
          "ways": FOUR}

    def test_over_8_mib(self):
        self.make([gf.scenario(ops=[self.OP], facts={"body_over_8mib": True})])
        lib = {"handle": "g_0123456789ab", "summary": {"seed": {}}, "evidence": {}}
        w3 = {"handle": UUID, "summary": {"seed": {}}, "evidence": {}}
        self.put(SID, self.OP, {"*": lib, "W3": w3}, sf={"in_service": False, "in_service_reason": "body > 8 MiB",
                                                         "download": {"url": "u"}, "delivery": "inline"})
        self.assertEqual(self.check(self.result(), "in_service")["outcome"], "pass")
        self.put(SID, self.OP, {"*": lib, "W3": w3}, sf={"in_service": False})
        self.assertIn("download absent", self.check(self.result(), "in_service")["detail"])

    def test_within_8_mib(self):
        self.make([gf.scenario(ops=[self.OP])])
        lib = {"handle": "g_0123456789ab", "summary": {"seed": {}}, "evidence": {}}
        self.put(SID, self.OP, {"*": lib, "W3": dict(lib, handle=UUID)}, sf={})
        res = self.result()
        self.assertEqual(res["facts"][SID]["over_8mib"], False)   # read from W1's document
        c = self.check(res, "in_service")
        self.assertEqual((c["outcome"], c["exercised"]), ("n/a", False))
        self.put(SID, self.OP, {"*": lib, "W3": dict(lib, handle=UUID)}, sf={"in_service": False})
        self.assertIn("in_service false on a body within 8 MiB", self.check(self.result(), "in_service")["detail"])


class MaxBytesAndProbes(ServiceBase):
    def test_max_bytes_262145(self):
        op = {"id": "walks-cap", "tool": "graphlet_walks", "ways": FOUR,
              "args": {"arm": "left", "rank": "id", "n": 10, "max_bytes": 262145}}
        self.make([gf.scenario(ops=[op])])
        refused = {"code": "bad_argument", "fields": {"field": "max_bytes"}}
        self.put(SID, op, {"*": gf.walks_answer()}, errors={"W3": refused})
        self.assertEqual(self.check(self.result(), "max_bytes_refused")["outcome"], "pass")
        self.put(SID, op, {"*": gf.walks_answer(), "W3": dict(gf.walks_answer(), handle=UUID)})
        self.assertIn("max_bytes 262145 on W3: answered", self.check(self.result(), "max_bytes_refused")["detail"])
        self.put(SID, op, {"*": None}, errors=dict((w, refused) for w in FOUR))
        self.assertIn("the library refused it too", self.check(self.result(), "max_bytes_refused")["detail"])

    def test_wrong_token_and_malformed_handle(self):
        ops = [{"id": "tok", "tool": "graphlet_summary", "args": {"max_bytes": 2048}, "ways": ["W3"],
                "probe": "wrong_token"},
               {"id": "bad", "tool": "graphlet_summary", "args": {"max_bytes": 2048}, "ways": ["W3"],
                "probe": "bad_handle"}]
        self.make([gf.scenario(ops=ops)])
        self.put(SID, ops[0], {}, errors={"W3": {"code": "forbidden", "fields": {}}}, ways=["W3"])
        self.put(SID, ops[1], {}, errors={"W3": {"code": "unknown_handle", "fields": {}}}, ways=["W3"])
        res = self.result()
        self.assertEqual(self.check(res, "forbidden")["outcome"], "pass")
        c = self.check(res, "bad_handle")
        self.assertEqual(c["outcome"], "fail")
        self.assertIn("bad_handle: error unknown_handle", c["detail"])
        self.assertEqual(res["required_service_checks"], ["forbidden", "bad_handle", "library_commit"])
        ops_ = res["scenarios"][SID]["ops"]
        self.assertEqual((ops_["tok"]["outcome"], ops_["tok"]["probe"]), ("agree", "wrong_token"))
        self.assertEqual((ops_["bad"]["outcome"], ops_["bad"]["disagreeing"]),
                         ("disagree", ["W3-probe bad_handle"]))

    def test_library_commit_from_the_prechecks(self):
        op = WalkDisplay.OP
        self.make([gf.scenario(ops=[op])])
        self.put(SID, op, {"*": gf.walks_answer(), "W3": dict(gf.walks_answer(), handle=UUID)},
                 sf={"walk_display": "first_parent"})
        self.assertEqual(self.check(self.result(), "library_commit")["outcome"], "pass")
        self.run_.record("W3", SID, "walks", "graphlet_walks", op["args"], answer=dict(gf.walks_answer(), handle=UUID),
                         service_fields={"walk_display": "first_parent"}, library={"commit": "e" * 40})
        self.assertEqual(self.check(self.result(), "library_commit")["outcome"], "fail")


class Blocking(unittest.TestCase):
    """VC2: any failed service check blocks; a required one must pass; n/a blocks only when required."""

    def test_rules(self):
        checks = [{"check": "in_service", "outcome": "fail", "detail": "in_service false on a body within 8 MiB"},
                  {"check": "hint", "outcome": "n/a", "exercised": False, "detail": "n/a: no compare call"},
                  {"check": "walk_display", "outcome": "n/a", "exercised": False, "detail": "n/a: no walks call"},
                  {"check": "label_evidence", "outcome": "inconclusive", "exercised": True, "detail": "facts"},
                  {"check": "library_commit", "outcome": "pass"}]
        out = dict(gs.blocking_service_checks(checks, required=["walk_display", "library_commit"]))
        self.assertEqual(sorted(out), ["in_service", "label_evidence", "walk_display"])
        self.assertTrue(out["in_service"].startswith("fail"))
        self.assertTrue(out["walk_display"].startswith("required, n/a"))
        self.assertNotIn("hint", out)
        self.assertEqual([n for n, _ in gs.blocking_service_checks(checks, (), set_aside=["in_service",
                                                                                         "label_evidence"])], [])
        self.assertEqual(gs.blocking_service_checks([{"check": "hint", "outcome": "error"}]),
                         [("hint", "outcome 'error' outside the vocabulary")])


class Required(unittest.TestCase):
    def test_from_the_op_plan(self):
        ops = [WalkDisplay.OP,
               {"id": "cmp", "tool": "graphlet_compare", "args": {"b": {"graphlet": SID}, "max_bytes": 2048}},
               {"id": "cap", "tool": "graphlet_summary", "args": {"max_bytes": 262145}},
               {"id": "splits", "tool": "graphlet_splits", "args": {"max_bytes": 262145}, "ways": ["W1", "W2"]}]
        man = gf.manifest([gf.scenario(ops=ops, facts={"body_over_8mib": True})])
        self.assertEqual(gs.required_service_checks(man),
                         ["walk_display", "label_evidence", "qualified_cause", "hint", "in_service",
                          "max_bytes_refused", "library_commit"])
        self.assertEqual(gs.required_service_checks(gf.manifest([gf.scenario(ways=["W1", "W2", "W4"])])), [])


class ByDesign(ServiceBase):
    OP = {"id": "walks-cap", "tool": "graphlet_walks", "ways": FOUR, "by_design": CAP_CASE["id"],
          "args": {"arm": "left", "rank": "id", "n": 10, "max_bytes": 262145}}
    REFUSED = {"code": "bad_argument", "fields": {"field": "max_bytes"}}

    def test_a_verified_pair(self):
        self.make([gf.scenario(ops=[self.OP])], table=[CAP_CASE])
        self.put(SID, self.OP, {"*": gf.walks_answer()}, errors={"W3": self.REFUSED})
        res = self.result()
        e = res["scenarios"][SID]["ops"]["walks-cap"]
        self.assertEqual(e["outcome"], "service_by_design", e)
        self.assertEqual(e["by_design"][0]["verified"], True)
        self.assertNotIn("disagreeing", e)
        self.assertEqual(res["summary"]["ops"]["service_by_design"], 1)
        self.assertEqual(res["scenarios"][SID]["outcome"], "pass")
        self.assertEqual(res["by_design_table"]["path"], "scripts/service_by_design.json")

    def test_a_pair_that_does_not_hold_names_its_side(self):
        self.make([gf.scenario(ops=[self.OP])], table=[CAP_CASE])
        ans = gf.walks_answer()
        self.put(SID, self.OP, {"*": ans, "W3": dict(ans, handle=UUID)})   # W3 answered: no difference at all
        e = self.result()["scenarios"][SID]["ops"]["walks-cap"]
        self.assertEqual(e["outcome"], "disagree")
        self.assertEqual(e["disagreeing"], ["by-design %s: service" % CAP_CASE["id"]])
        self.assertIn("W3 answered, the case says error bad_argument", e["by_design"][0]["detail"])
        # the library ways refused too: the library side
        self.put(SID, self.OP, {"*": None}, errors=dict((w, self.REFUSED) for w in FOUR))
        e = self.result()["scenarios"][SID]["ops"]["walks-cap"]
        self.assertEqual(e["by_design"][0]["side"], "library")

    def test_library_pairs_still_compared_and_request_and_table_sides(self):
        self.make([gf.scenario(ops=[self.OP])], table=[CAP_CASE])
        ans = gf.walks_answer()
        self.put(SID, self.OP, {"*": ans, "W2": dict(ans, total=5)}, errors={"W3": self.REFUSED})
        e = self.result()["scenarios"][SID]["ops"]["walks-cap"]
        self.assertEqual((e["outcome"], e["disagreeing"]), ("disagree", ["W1-W2"]))
        self.run_.cleanup()
        op = dict(self.OP, args=dict(self.OP["args"], max_bytes=2048))   # the args predicate does not hold
        self.make([gf.scenario(ops=[op])], table=[CAP_CASE])
        self.put(SID, op, {"*": ans}, errors={"W3": self.REFUSED})
        e = self.result()["scenarios"][SID]["ops"]["walks-cap"]
        self.assertEqual(e["by_design"][0]["side"], "request")
        self.run_.cleanup()
        op = dict(self.OP, by_design="graphlet.no_such_case")
        self.make([gf.scenario(ops=[op])], table=[CAP_CASE])
        self.put(SID, op, {"*": ans}, errors={"W3": self.REFUSED})
        e = self.result()["scenarios"][SID]["ops"]["walks-cap"]
        self.assertEqual((e["by_design"][0]["side"], e["outcome"]), ("table", "disagree"))

    def test_the_case_schema(self):
        self.assertEqual(cmp.table_problems([CAP_CASE]), [])
        bad = copy.deepcopy(CAP_CASE)
        bad["agreement"] = "service_only"
        bad["graphlet"].update(tools=[], args=[{"field": "x.", "equals": 1}], zz=1,
                               library={"answered": True, "error": {"code": "a"}},
                               service={"error": {"code": ""}, "service_fields": []})
        text = "; ".join(cmp.table_problems([bad]))
        for want in ("unknown key(s) ['zz']", "agreement by_design", "non-empty list of graphlet tool names",
                     "path predicates over args_resolved only", "graphlet.library states exactly one",
                     "graphlet.service.error must be an object with a code", "service_fields: an object"):
            self.assertIn(want, text)
        lib_sf = copy.deepcopy(CAP_CASE)
        lib_sf["graphlet"]["library"]["service_fields"] = {"x": 1}
        self.assertIn("on the service side only", "; ".join(cmp.table_problems([lib_sf])))
        # the real table loads unchanged (it has no area G case)
        cases, info = cmp.load_by_design(cmp.DEFAULT_BY_DESIGN)
        self.assertFalse(any("graphlet" in c for c in cases.values()))

    def test_a_refused_table_is_a_warning_and_every_pair_unverified(self):
        bad = copy.deepcopy(CAP_CASE)
        bad["graphlet"]["tools"] = []
        self.make([gf.scenario(ops=[self.OP])], table=[bad])
        self.put(SID, self.OP, {"*": gf.walks_answer()}, errors={"W3": self.REFUSED})
        res = self.result()
        self.assertTrue(any("by-design table refused" in w for w in res["warnings"]))
        self.assertEqual(res["scenarios"][SID]["ops"]["walks-cap"]["by_design"][0]["side"], "table")


class Facts(unittest.TestCase):
    def test_document_facts_of_the_fixtures(self):
        import checks
        f = gf.Fixture()
        facts = gs.document_facts(f.doc(0), checks)
        self.assertIs(facts["cut_lists"], False)
        self.assertGreater(facts["max_labels_per_node"], 0)
        # the body without its J line (the Z count may differ by one digit)
        self.assertLessEqual(abs(facts["body_bytes"] - len(f.doc(0, envelope=False).encode("utf-8"))), 1)
        self.assertEqual(gs.document_facts(None, checks)["cut_lists"], None)


if __name__ == "__main__":
    unittest.main()
