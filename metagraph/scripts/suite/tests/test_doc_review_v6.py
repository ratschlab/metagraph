"""The external review of the v6 results document (2026-10-18; V6-001 .. V6-307): each agreed finding's fix in the
generator, tested where behaviour changes. The review's ids name the tests. Offline; stdlib only."""

import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.dirname(HERE)
DOC = os.path.join(SUITE, "doc")
for p in (SUITE, HERE, DOC):
    if p not in sys.path:
        sys.path.insert(0, p)

import suite_common as sc  # noqa: E402

import blocks  # noqa: E402,F401
import cells  # noqa: E402,F401
import facts  # noqa: E402
import latex  # noqa: E402
import sections  # noqa: E402

PROSE_DIR = os.path.join(DOC, "prose")


class _Stub(object):
    """A bare object with the attributes given (a Run or a Context stand-in)."""

    def __init__(self, **kw):
        for k, v in kw.items():
            setattr(self, k, v)


def prose(name):
    return sc.read_json(os.path.join(PROSE_DIR, "%s.json" % name))


class GatePopulationTest(unittest.TestCase):
    """V6-001: the gate's required population after the policy's automatic exclusions, beside the verdict."""

    def _run(self):
        policy = {
            "arms": ["direct-local", "mcp"],
            "selection": {"routes": ["POST /traverse"], "scenarios": ["T-a", "T-b", "T-c", "T-d"]},
            "exclusions": [
                {"scenario": "T-b", "arms": ["direct-local", "mcp"], "rules": ["budget_stress"], "reason": "load"},
                {"scenario": "T-c", "arms": ["mcp"], "rules": ["manifest_not_implemented", "table_not_expressible"],
                 "reason": "x"},
                {"scenario": "T-d", "arms": ["mcp"], "repeats": [0], "rules": ["marks_cover_repeat"], "reason": "y"},
                {"scenario": "T-a", "arms": ["mcp"], "repeats": [1], "rules": ["declared"], "reason": "owner"}],
            "required": [{"scenario": "T-a", "arm": "direct-local", "repeats": [0]},
                         {"scenario": "T-a", "arm": "mcp", "repeats": [0]},
                         {"scenario": "T-c", "arm": "direct-local", "repeats": [0]},
                         {"scenario": "T-d", "arm": "direct-local", "repeats": [0]},
                         {"scenario": "T-d", "arm": "mcp", "repeats": [1]}]}
        rows = [{"scenario": sid, "arm": arm, "repeat": rep, "outcome": out} for sid, arm, rep, out in (
            ("T-a", "direct-local", 0, "pass"), ("T-a", "mcp", 0, "service_by_design"), ("T-a", "mcp", 1, "pass"),
            ("T-b", "direct-local", 0, "inconclusive"), ("T-c", "direct-local", 0, "pass"),
            ("T-c", "mcp", 0, "not_implemented_mcp"), ("T-d", "direct-local", 0, "pass"),
            ("T-d", "mcp", 0, "inconclusive"), ("T-d", "mcp", 1, "pass"), ("T-e", "direct-local", 0, "fail"))]
        return _Stub(policy=policy, rows=rows, order=["T-a", "T-b", "T-c", "T-d", "T-e"], policy_path="/x/p.json",
                     declared={"exclusions": [{"scenario": "T-a", "arms": ["mcp"], "repeats": [1]}]},
                     gate={"passed": True, "problems": [], "counts": {"judged": 5, "passed": 4}})

    def test_population_and_ledger_text(self):
        run = self._run()
        pop = sections.gate_population(run)
        self.assertEqual(pop["selected"], 4)
        self.assertEqual(pop["required"], [("direct-local", 3, 3), ("mcp", 2, 2)])
        self.assertEqual(len(pop["auto"]), 3, "the owner's declared entry is not an automatic exclusion")
        self.assertEqual(sections.gate_auto_text(run),
                         "3 entries: 1 budget-stress test on direct-local and mcp, 1 request the MCP tools cannot "
                         "express on mcp and 1 repeat whose property checks the marks all set aside on mcp (repeat 0 "
                         "only)")
        self.assertEqual(sections.gate_required_text(run),
                         "3 scenarios (3 executions) on direct-local and 2 scenarios (2 executions) on mcp")
        text = sections.gate_population_text(run)
        self.assertIn("Required means selected by route (the scenarios of POST /traverse (4 of the run's 5)", text)
        self.assertIn("no decision of the owner (3 entries: ", text)
        self.assertIn("not declared exempt by the owner (1 exclusion; I.4, appendix G)", text)
        self.assertIn("5 executions in all. The gate certifies those executions, not the excluded ones", text)
        decision = sections.gate_decision_text(run)
        self.assertTrue(decision.startswith("The gate passed: every planned execution was recorded, every execution "
                                            "the policy requires passed or matched the by-design case"), decision)
        self.assertIn(text, decision)
        self.assertIn("removed by the policy's automatic exclusions", sections.gate_inconclusive_clause(run))

    def test_judged_outcomes(self):
        """The judged executions by outcome, as gate.judge counts them: T-b (budget_stress), T-c via MCP, T-d's MCP
        repeat 0, T-a's declared MCP repeat 1 and T-e (outside the selection) are not judged."""
        run = self._run()
        self.assertEqual(sections.gate_judged_outcomes(run), {"pass": 4, "service_by_design": 1})
        self.assertEqual(sections.judged_outcomes_text(run, {"passed": 4}), " (4 passed, 1 service-side by design)")
        # a count the gate's passed does not confirm falls back to the gate's own number
        self.assertEqual(sections.judged_outcomes_text(run, {"passed": 2}), ", 2 passed")
        self.assertEqual(sections.gate_population(_Stub(policy=None)), None)
        self.assertEqual(sections.gate_auto_text(_Stub(policy=None)), "")

    def test_prose_criterion(self):
        """The blocking criterion reads 'required by the policy on that access path and repeat', never 'in the
        policy's selection and not declared exempt' (which leaves out the automatic exclusions)."""
        texts = [str(prose("overview")), str(prose("glossary"))]
        for t in texts:
            self.assertNotIn("in the policy's selection and not declared exempt", t)
        self.assertIn("the policy requires on that access path and repeat",
                      prose("overview")["sections"]["decision"]["verdict"][0])
        self.assertIn("{gate_auto_text}", prose("overview")["sections"]["provenance"]["exemptions"][0])


class FullyCheckedTest(unittest.TestCase):
    """V6-201: 'fully checked' holds over the property checks no property-level mark sets aside; the chapter's MCP
    sentence names the agreeing pairs that omit a marked check, which stays unverified via MCP."""

    def test_definition(self):
        ctx = _Stub(units=_Stub(by_scenario={}))
        self.assertEqual(facts.mcp_full_definition(ctx, []), "every property check evaluated on both")
        self.assertEqual(facts.mcp_full_definition(ctx, ["T-labels-0003", "T-labels-0005", "T-labels-0006"]),
                         "every property check evaluated on both, except a check a property-level mark of the manifest "
                         "sets aside, not compared and unverified via MCP (on three of them: T-labels-0003, -0005 and "
                         "-0006)")
        many = ["P-exact-%04d" % i for i in range(1, 9)] + ["P-scope-0001", "P-edge-0004"]
        self.assertTrue(facts.mcp_full_definition(ctx, many).endswith(
            "(on ten of them: P-exact 8, P-scope 1, P-edge 1)"))

    def test_prose(self):
        for ch in ("R", "T", "P", "L", "X", "M", "Q", "N", "D"):
            text = str(prose(ch))
            self.assertEqual(text.count("{mcp_full_def}"), 2, ch)
            self.assertNotIn("(every property check evaluated on both)", text, ch)
        depth = prose("overview")["sections"]["reading"]["paragraphs"][3]
        self.assertIn("except a check a property-level mark of the manifest sets aside: that check is not compared "
                      "and stays unverified via MCP", depth)
        gloss = [e for e in prose("glossary")["entries"] if e["term"] == "agreement class, depth"][0]["definition"]
        self.assertIn("but a check a property-level mark of the manifest sets aside", gloss)


class _Arm(object):
    def __init__(self, lines, present=True):
        self.lines, self.present = lines, present


class ExceptionSummaryTest(unittest.TestCase):
    """V6-302: an exception summary names the actual request, the unverified property per access path and the
    recorded cause; never a menu of possible causes, never the whole unit's request."""

    def setUp(self):
        self.run = _Stub(is_property=lambda name: name not in ("status", "invariants", "budget_respected"),
                         budget_ms=lambda sid: 110000)

    def _states(self, **arms):
        return dict((arm, _Arm(arms.get(arm.replace("-", "_"), []), bool(arms.get(arm.replace("-", "_")))))
                    for arm in ("direct-local", "direct-remote", "mcp"))

    def test_reasons(self):
        hidden = {"outcome": "inconclusive", "http": 200, "repeat": 0, "checks": [
            {"check": "capabilities_doc", "outcome": "not_exposed", "detail": "field not in the answer: features"}]}
        late503 = {"outcome": "inconclusive", "http": 503, "error": {"code": "deadline"}, "repeat": 0,
                   "wall_ms": 110200, "checks": []}
        walk = {"outcome": "inconclusive", "http": 200, "repeat": 0, "checks": [
            {"check": "walk_in_record", "outcome": "inconclusive",
             "detail": "left: the record followed for 3618 of 10000 bp, stopped by resource_limit"}]}
        self.assertEqual(facts.unshown_reason(self.run, "C", hidden, "mcp", short=True),
                         "capabilities_doc not exposed via MCP (the answer lacks features)")
        self.assertEqual(facts.unshown_reason(self.run, "C", hidden, "mcp"),
                         "capabilities_doc not exposed via MCP: field not in the answer: features")
        self.assertEqual(facts.unshown_reason(self.run, "P", dict(late503), "direct-local"),
                         "an allowed 503 deadline (no answer to check)")
        self.assertEqual(facts.unshown_reason(self.run, "T", walk, "direct-local", short=True),
                         "walk_in_record could not decide")
        self.assertEqual(facts.unshown_reason(self.run, "T", {"outcome": "pass", "checks": []}, "mcp"), "")
        self.assertEqual(facts.exec_phrase(self.run, "C", hidden, "mcp", wall="0.3 s"),
                         "inconclusive (capabilities_doc not exposed via MCP: field not in the answer: features)")
        states = {"C": self._states(mcp=[hidden, dict(hidden, repeat=1)]),
                  "T": self._states(direct_local=[walk], direct_remote=[walk]),
                  "T2": self._states(direct_local=[walk], direct_remote=[walk])}
        self.assertEqual(facts.unshown_index_groups(self.run, ["C", "T", "T2"], states.get),
                         [(("mcp", "capabilities_doc not exposed via MCP (the answer lacks features)"), ["C"]),
                          (("direct-local and direct-remote", "walk_in_record could not decide"), ["T", "T2"])])
        self.assertEqual(facts.unshown_reasons_text(self.run, "C", states["C"]),
                         "mcp r0, r1: capabilities_doc not exposed via MCP: field not in the answer: features")

    def test_requested_is_the_scenario(self):
        """The unit's property line is a block's `requested` only for a unit of one test."""
        import inspect
        src = inspect.getsource(sections._unit_blocks)
        self.assertIn("property_text = ctx.fill(prop) if prop and len(members) == 1 else None", src)
        src = inspect.getsource(blocks.interpretation)
        self.assertIn('requested += " (property check%s %s)"', src)
        self.assertIn("on every access path: none, each property check passed on some access path or repeat; ", src)


class WarmUpTimeoutTest(unittest.TestCase):
    """V6-204: a test without a property check that expects a 200 (a warm-up) is not 'the status is the claim':
    its detail block says it declares no property of its own and that the allowed 503 left its generic checks
    undecided; a refusal keeps 'the status, code and text are the claim'."""

    def test_property_cell(self):
        line = {"outcome": "inconclusive", "checks": [{"check": "status", "outcome": "inconclusive"}]}
        generic = lambda name: False  # noqa: E731
        self.assertEqual(blocks._property_cell(line, generic, refusal=False), "n/a (no property declared)")
        self.assertEqual(blocks._property_cell(line, generic), "n/a (status is the claim)")

    def test_interpretation(self):
        warm = {"id": "P-warm-0025", "title": "warm-up of P-heavy-0001: its request", "warms": "P-heavy-0001",
                "expect": [{"check": "status", "args": {"http": 200}}, {"check": "invariants", "args": {}}]}
        refusal = {"id": "R-edge-0002", "title": "k-1 bases: refused",
                   "expect": [{"check": "status", "args": {"http": 400, "text": "Query shorter than k"}}]}
        line = {"outcome": "inconclusive", "http": 503, "error": {"code": "deadline"}, "repeat": 0, "wall_ms": 1000,
                "checks": [{"check": "status", "outcome": "inconclusive"}]}
        scen = {"P-warm-0025": warm, "R-edge-0002": refusal}
        run = _Stub(scenario=scen.get, is_property=lambda name: False, budget_ms=lambda sid: 110000)
        states = {"direct-local": _Arm([line]), "direct-remote": _Arm([], False), "mcp": _Arm([], False)}
        for st in states.values():
            st.missing_planned, st.glyph_key, st.ne = False, "inconclusive", False
        req, _, unv = blocks.interpretation(run, "P-warm-0025", states)
        self.assertIn("(no property check: the generic checks alone; the property is its target's, P-heavy-0001)", req)
        self.assertIn("no property of its own declared", unv)
        self.assertIn("left undecided on direct-local: an allowed 503 deadline (no answer to check)", unv)
        self.assertNotIn("status, code and text are the claim", unv)
        line_r = dict(line, http=400, error={"code": "invalid_request"}, outcome="pass")
        states["direct-local"] = _Arm([line_r])
        states["direct-local"].missing_planned, states["direct-local"].glyph_key = False, "pass"
        states["direct-local"].ne = False
        req, _, unv = blocks.interpretation(run, "R-edge-0002", states)
        self.assertIn("the status, code and text are the claim", unv)
        self.assertIn("(no property check: the refusal's status, code and text are the claim)", req)


class WarmUpCausalTest(unittest.TestCase):
    """V6-203: the P-warm prose states the warm-up's intent and the gate exclusion, and that the run does not
    isolate what decides completion; it never rules out the deployment."""

    def test_prose(self):
        text = prose("P")["families"]["P-warm"]["prose"]
        self.assertNotIn("not on the deployment", text)
        self.assertIn("this run does not isolate what decides whether it completes within the budget", text)
        self.assertIn("host-cold is the suite's execution-order label, not a measured cache state", text)
        self.assertIn("Searches meant to read, before a /pattern scenario", text)


class ProteinAdmissionTest(unittest.TestCase):
    """V6-202: the admission-only sentence names the cap and refusal tests, never 'the passing protein tests' as a
    whole, and points at M-protein-level1 for the recovery where it passed."""

    def test_sentence(self):
        text = prose("M")["summary"][0]
        self.assertNotIn("The passing protein tests are admission tests", text)
        self.assertIn("The passing protein-cap and refusal tests are admission tests, not recoveries", text)
        filled = facts.fill(text, {"unit": {"M-protein-level1": {"n_pass_dl": 7, "direct_all_pass": True}}})
        self.assertIn("which only M-protein-level1's passes show", filled)
        filled = facts.fill(text, {"unit": {"M-protein-level1": {"n_pass_dl": 0, "direct_all_pass": False}}})
        self.assertNotIn("which only M-protein-level1's passes show", filled)


class DefectPinTest(unittest.TestCase):
    """V6-304: the known-defect pin's note names the observation that would support the repair; any other failure
    is diagnosed, and a failure alone never licenses revising the expectation."""

    def test_note(self):
        note = prose("N")["units"]["N-usecase-defect"]["note"]
        self.assertNotIn("a failure signals the repair", note)
        self.assertIn("consistent with the repair only when its answer shows the k-mer still in the graph and its "
                      "header now supporting it", note)
        self.assertIn("is diagnosed like any failure", note)


class SeedPhaseBudgetTest(unittest.TestCase):
    """V6-002: the catalogue's budget_respected never places the seed phase outside the /traverse budget; a seed-only
    answer is a limit of the check (a coverage gap), and the 503 contract is /pattern's and /resolve's."""

    def test_catalogue(self):
        b = prose("checks")["checks"]["budget_respected"]
        text = str(b)
        self.assertNotIn("which the budget does not govern", text)
        self.assertIn("a coverage gap of the suite, not an exemption of the server", b["inconclusive"])
        self.assertIn("on /pattern and /resolve past the budget the answer is a 503", b["asserts"])
        row = [r for r in b["decision_table"]["rows"] if r[0] == "a /traverse answer with no walked seed"][0]
        self.assertIn("no check of this suite judges its termination", row[3])
        note = prose("N")["details"]["N-hard-0022"]["note"]
        self.assertNotIn("phases the budget does not cover", note)
        self.assertIn("read on the seed's own clock, which starts with the seed phase", note)


class JunctionPremiseTest(unittest.TestCase):
    """V6-101: the junction's exact header runs rest on the builder's whole-record look-up (areas._record_runs),
    which the prose states; the local continuations explain the runs, they do not define them."""

    def test_prose_and_builder(self):
        import inspect
        import areas
        text = str(prose("R"))
        self.assertIn("looks up every k-mer of the 170 bp query in the whole of each fetched record", text)
        self.assertIn("no junction k-mer occurs elsewhere in either record", text)
        self.assertIn("every query k-mer looked up in the whole record", prose("R")["units"]["R-junction"]["property"])
        src = inspect.getsource(areas._record_runs)
        self.assertIn("if q[i:i + K] in rec.seq:", src)   # the whole record, not the window's neighbourhood
        for version in (5, 6):
            m = sc.read_json(os.path.join(SUITE, "manifest", "suite_v%d.json" % version))
            for s in m["scenarios"]:
                if not s["id"].startswith("R-junction"):
                    continue
                t = s.get("truth") or {}
                runs = [e["args"]["runs"] for e in s["expect"] if e["check"] == "resolve_runs" and "runs" in e["args"]]
                self.assertEqual(runs, [[[0, t["a_run_end"]]], [[t["b_run_start"], 170]]], s["id"])


class AnnotateWitnessTest(unittest.TestCase):
    """V6-102: annotate mode's all-successor behaviour is the SPEC's contract; the test checks a source-route witness
    and no claim, property line or question says every successor was shown."""

    def test_prose(self):
        t = prose("T")
        fam = t["families"]["T-labels"]
        self.assertTrue(fam["claim"].startswith("in annotate mode a returned route follows the source record"))
        self.assertIn("not an enumeration of every successor", fam["prose"])
        self.assertNotIn("annotate mode lists every structural successor", t["units"]["T-labels"]["property"])
        self.assertNotIn("Does annotate mode follow every structural successor", t["units"]["T-labels"]["question"])


class StrandRefusalStageTest(unittest.TestCase):
    """V6-103: T-strand names the seed refusal the run observed (facts' seed_refusal, from the answer's text) and
    bounds the claim: a graph-absence refusal does not show strand-specific label support."""

    def test_prose_branches(self):
        text = prose("T")["families"]["T-strand"]["prose"]
        graph = facts.fill(text, {"x": {"T-strand-0001": {"dl": {"seed_refusal": "graph"}, "mcp": {"class": ""}}},
                                  "manifest_version": 6, "n_footnoted": 0})
        self.assertIn("a plain absence, so the test shows that this input is refused and not that the header tells "
                      "the two orientations apart", graph)
        self.assertNotIn("the traverse-side statement of the strand fact", graph)
        label = facts.fill(text, {"x": {"T-strand-0001": {"dl": {"seed_refusal": "label"}, "mcp": {"class": ""}}},
                                  "manifest_version": 6, "n_footnoted": 0})
        self.assertIn("the traverse-side statement of the strand fact", label)

    def test_fact_from_text(self):
        import inspect
        src = inspect.getsource(facts.arm_facts)
        self.assertIn('"graph" if "not fully present in the graph" in text', src)


class ZeroEqualsAbsentTest(unittest.TestCase):
    """V6-305: N-zero reads exact 0 per query length (a context for k bases, a path for a longer query) and a
    shortened graph_runs as an absence from the graph, never from the records or from biology."""

    def test_prose(self):
        text = prose("N")["families"]["N-zero"]["prose"]
        self.assertNotIn("no retained k-mer holds the pattern", text)
        self.assertIn("no path of the graph spells the whole query on that strand", text)
        self.assertIn("some k-mer of the query is absent from the graph", text)
        self.assertNotIn("in no record of the index", text)


class GroundTruthWordingTest(unittest.TestCase):
    """V6-008: the overview keeps the server's completeness statement apart from an independent enumeration.
    V6-301: chapter Q names exactly the expectations that depend on the assumed human non-carrier."""

    def test_overview_completeness(self):
        text = prose("overview")["sections"]["methods"]["paragraphs"][2]
        self.assertNotIn("which no check of this suite establishes other than through a complete list", text)
        self.assertIn("enumerates no matching result of the index independently", text)
        self.assertIn("lower bounds on recovery", text)

    def test_q_assumption(self):
        q = prose("Q")
        fam = q["families"]["Q-pred"]["prose"]
        self.assertIn("The assumption decides the expectations of `all(A, human)`, `at_least` 2 of (A, human) and "
                      "`and(any A, not any human)`; `none(A)` is expected false from A's known presence alone", fam)
        self.assertNotIn("the human column (hsa) does not;", str(q))
        self.assertIn("none(A) rests on A's presence alone", q["units"]["Q-pred-forms"]["property"])


class PrimerCountTest(unittest.TestCase):
    """V6-303: the forward-only primer counts are counts of matching contexts, not of records or gene copies."""

    def test_prose(self):
        text = prose("N")["families"]["N-usecase"]["prose"]
        self.assertNotIn("the records holding the gene on the other strand", text)
        self.assertIn("counts the matching contexts on the strand asked, not records and not verified copies of the "
                      "gene", text)


class LongPeptideCoverageTest(unittest.TestCase):
    """V6-205: the coverage gap concerns positive recovery of long peptides; X-junction's paths searches are the
    negative control, named apart."""

    def test_prose(self):
        text = str(prose("X"))
        self.assertNotIn("Peptides longer than k are tested as supported paths with a predicate alone, never as", text)
        self.assertIn("Positive recovery of a peptide whose coding DNA (three bases per residue) is longer than k", text)
        self.assertIn("X-junction searches its 20-residue windows as `paths`, but for the absence of the spliced source "
                      "DNA", text)


class RecordBoundaryTest(unittest.TestCase):
    """V6-206: k - L + 1 is an interior-window rule; both record boundaries limit the contexts, the unit tests
    starts."""

    def test_prose(self):
        text = prose("P")["families"]["P-scope"]["prose"]
        self.assertNotIn("A record start is the one place", text)
        self.assertIn("a record boundary leaves fewer (a window at a record's end too), and this unit exercises record "
                      "starts alone", text)


class ExpectedStopTest(unittest.TestCase):
    """V6-209: chapter M's purpose separates a budget stop that keeps a witness check from deciding (inconclusive)
    from a test whose property is the stated stop (a pass)."""

    def test_prose(self):
        text = prose("M")["purpose"][0]
        self.assertNotIn("a search the budget ended counts as inconclusive, never as a pass", text)
        self.assertIn("a test whose property is a correctly stated stop", text)


class DeterminismStopsTest(unittest.TestCase):
    """V6-306: chapter D and the catalogue keep time-limited stops (not repeated, not promised) apart from a work-limit
    stop (max_steps), which D-M-order-0003 repeats."""

    def test_prose(self):
        gap = prose("D")["coverage_gap"]
        self.assertIn("the repeat of a stop by elapsed time or of a 503 deadline is not asked", gap[0])
        self.assertIn("a stop by a work limit, which no clock decides, is repeated like any answer", gap[0])
        self.assertIn("(a stop or a cut by elapsed time, a 503 deadline)", gap[1])
        det = prose("checks")["checks"]["deterministic"]["inconclusive"]
        self.assertIn("A stop by a work limit (`max_steps`), which no clock decides, is compared like any answer", det)


class DemonstratedStopTest(unittest.TestCase):
    """V6-106: a demonstrated stop is not filed under 'unsupported': the coverage heading is neutral and R's lead
    depends on whether the stop was observed; the MCP limit is named apart."""

    def test_heading_and_lead(self):
        import inspect
        self.assertIn("\\parhead{Coverage and remaining limits}", inspect.getsource(sections.chapter))
        gap = prose("R")["coverage_gap"][0]
        shown = facts.fill(gap, {"fam": {"R-budget": {"stop_observed": True, "stop_observed_text": "R-budget-0001 on "
                                                      "direct-local"}},
                                 "x": {"R-budget-0001": {"mcp": {"standing": "not exercised"}}}, "manifest_version": 6})
        self.assertTrue(shown.startswith("*Shown in this run.* The stated stop"), shown[:80])
        self.assertIn("it remains untested through MCP, which cannot express these requests (see B)", shown)
        not_shown = facts.fill(gap, {"fam": {"R-budget": {"stop_observed": False}}, "manifest_version": 6})
        self.assertTrue(not_shown.startswith("*What remains unsupported.*"), not_shown[:80])


class WarmedTwinMarkerTest(unittest.TestCase):
    """V6-107: R-budget explains why a stop run right after its warming twin can carry the host-cold marker."""

    def test_prose(self):
        text = prose("R")["families"]["R-budget"]["prose"]
        v6 = facts.fill(text, {"manifest_version": 6, "n_warm_up_direct_pairs": 0})
        self.assertIn("the marker tracks the earlier executions of the same scenario (I.2), not every read of those "
                      "rows", v6)
        self.assertNotIn("the marker tracks", facts.fill(text, {"manifest_version": 4}))


class PlantedTitleTest(unittest.TestCase):
    """V6-208: the planted-substitution unit names t (planted) and m (max_mismatches) as its two dimensions."""

    def test_mappings(self):
        for v in (4, 5, 6):
            text = open(os.path.join(DOC, "units", "suite_v%d.json" % v), encoding="utf-8").read()
            self.assertNotIn("planted substitutions by m and max_mismatches", text, v)
            self.assertIn('"unit_title": "planted substitutions t and max_mismatches m"', text, v)


class FootnoteColourTest(unittest.TestCase):
    """V6-005: every legend describes the footnote digit as rendered (doc.sty's fnink, dark grey), never as red."""

    def test_one_description(self):
        import inspect
        sty = open(os.path.join(DOC, "templates", "doc.sty"), encoding="utf-8").read()
        self.assertIn("\\definecolor{fnink}{HTML}{37474F}", sty)
        reading = prose("overview")["sections"]["reading"]
        self.assertIn("a dark-grey superscript digit", reading["annotations"]["footnote"])
        mark = [m for m in reading["marks"] if m.get("mark") == "superscript digit"][0]
        self.assertIn("dark grey, not red", mark["meaning"])
        src = inspect.getsource(sections)
        self.assertNotIn("a red superscript digit", src)
        self.assertNotIn("a red superscript digit", str(prose("overview")))


class PairsLegendTest(unittest.TestCase):
    """V6-307: an ordinary pairs table's legend defines every timing row its cells print (be included) and no digest
    mark it cannot print; the repeated (chapter D) legend keeps the marks."""

    def test_legends(self):
        legend = cells.TIMING_LEGEND["pairs"]
        self.assertIn("be the host call as the service timed it, in the same repeat order (MCP only)", legend)
        self.assertNotIn("≠", legend)
        _, repeated = cells.repeats_texts(3)
        self.assertIn("= and ≠", repeated)
        self.assertIn("be the host call as the service timed it", repeated)


class RepeatsMetaTest(unittest.TestCase):
    """V6-105: a detail block's metadata calls the manifest's count the manifest's and names the recorded executions
    where they differ from it."""

    def test_source(self):
        import inspect
        src = inspect.getsource(blocks.detail_block)
        self.assertIn('text = "%s repeat%s in the manifest"', src)
        self.assertIn('text += "; recorded executions %s (the rows below)"', src)


class UnitHeadingTest(unittest.TestCase):
    """V6-108: a unit heading's key never runs into its title: a key that does not fit moves whole to the next line,
    flush right, with a minimum gap (the TeXbook's \\signed construction)."""

    def test_macro(self):
        sty = open(os.path.join(DOC, "templates", "doc.sty"), encoding="utf-8").read()
        line = [l for l in sty.splitlines() if l.startswith("\\newcommand\\unithead")][0]
        self.assertIn("\\unskip\\nobreak\\hfil\\penalty50\\hskip1.5em\\hbox{}\\nobreak\\hfil\\mbox{", line)
        self.assertIn("\\parfillskip=0pt", line)
        self.assertNotIn("#2\\hfill{", line)


class PhasePopulationTest(unittest.TestCase):
    """V6-006: every phase box's left-out values (zeros, values below 0 ms) are counted from the plotted arrays and
    stated, and the tick label's n reconciles with the population."""

    def _ctx(self, values):
        ctx = _Stub()
        orig = sections.phase_population
        self.addCleanup(setattr, sections, "phase_population", orig)
        sections.phase_population = lambda c: values
        return ctx

    def test_omitted(self):
        values = {"submit_ms": [1.0, 2.0, 3.0], "queue_ms": [1.0, 2.0], "backend_ms": [5.0, 6.0],
                  "terminal_lag_ms": [-10.6, -0.8, 100.0, 200.0, 300.0], "paging_ms": [0, 0, 10.0, 20.0]}
        ctx = self._ctx(values)
        self.assertEqual(sections.phase_omitted_text(ctx),
                         "terminal lag is below 0 ms in 2 of 5 (-10.6 to -0.8 ms) and paging is 0 ms where the results "
                         "came with the terminal status in 2 of 4")
        self.assertEqual(sections.phase_omitted_short(ctx),
                         "terminal lag without its 2 values below 0 ms, paging without its 2 zeros")
        boxes = dict((b[0], b) for b in sections.phase_quantiles(ctx))
        self.assertEqual(boxes["terminal lag"][6:], (3, 0, 2))
        self.assertIn("(n 3 of 5)", latex.phase_boxplot([boxes["terminal lag"]], "fig:x"))
        values["terminal_lag_ms"] = [-8.7, 100.0, 200.0]
        self.assertIn("below 0 ms in 1 of 3 (-8.7 ms)", sections.phase_omitted_text(ctx))


class RelationArgumentsTest(unittest.TestCase):
    """V6-007: appendix C prints every relation argument in full, the compared paths readable and never cut, so the
    lead's 'printed in full' holds."""

    def test_paths(self):
        paths = [["release", "release"], ["server_instance", "attempts.server_instance"]] + [
            ["field_%02d" % i, "field_%02d" % i] for i in range(40)]
        text = sections.relation_args_text({"paths": paths})
        self.assertTrue(text.startswith("42 paths, a's against b's: release; server_instance against "
                                        "attempts.server_instance; field_00"), text[:120])
        self.assertTrue(text.endswith("field_39"))
        self.assertNotIn("…", text)
        self.assertEqual(sections.relation_args_text({"paths": [["", "pattern"]]}),
                         "1 path, a's against b's: the whole answer against pattern")
        self.assertEqual(sections.relation_args_text({"paths": [["k", "pattern.k"], ["mask", "pattern.mask"]]}),
                         "2 paths, each of a against pattern.<the same path> of b: k, mask")
        self.assertEqual(sections.relation_args_text({"a_from_kmer": 0, "b_from_kmer": 0, "labels": ["511145"]}),
                         'a_from_kmer 0, b_from_kmer 0; labels ["511145"]')
        import inspect
        self.assertIn("esc_wrap(arg_text))", inspect.getsource(sections.relations_appendix))
        self.assertIn("The arguments of each relation are printed in full under its title",
                      prose("overview")["appendix_C"]["paragraphs"][0])


if __name__ == "__main__":
    unittest.main()
