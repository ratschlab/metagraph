"""The run's snapshot of what judges its arms (lane S7r, P6-end-to-end-03 of the two-day review): run_suite.sh copies
suite_common.SNAPSHOT_FILES into RUN_DIR/scripts once (SNAPSHOT.json with their sha256); direct-remote runs from it,
mex gets it, the MCP runner gets it as {checks_dir}; run_direct records per session the sha256 of the checks and
suite_common it imported; merge_run states every arm judged with other checks (an error with --strict-identity and for
the gate) and makes outcomes stricter with the snapshot's checks."""

import contextlib
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.dirname(HERE)
sys.path.insert(0, SUITE)
sys.path.insert(0, HERE)

import compare as cmp  # noqa: E402
import gate  # noqa: E402
import merge_run  # noqa: E402
import report  # noqa: E402
import run_direct  # noqa: E402
import suite_common as sc  # noqa: E402
from fake_host import FakeHost  # noqa: E402
from test_compare_report import T0, line, write_arm  # noqa: E402
from test_gate import make_run, manifest_of, mcp_lines, policy_of  # noqa: E402
from test_orchestrator import SuiteCase, mcp_scenarios  # noqa: E402
from test_runner import basic_scenarios, make_manifest, scenario  # noqa: E402

# checks.py as a change landing in the tree during a run would leave it: every execution fails
FAILING_CHECKS = '''

def evaluate(scenario, answer, arm, previous=None, client=None):
    return [{"check": "status", "outcome": "fail", "detail": "a checks change landed during the run"}]


def scenario_outcome(results):
    return "fail"
'''
# checks.py as a change landing in the tree during a run would leave it: every relation fails
FAILING_RELATIONS = '''

def evaluate_relation(rel, a_scenario, a_answer, b_scenario, b_answer, arm):
    return {"relation": rel.get("id"), "type": rel.get("type"), "a": rel.get("a"), "b": rel.get("b"),
            "outcome": "fail", "detail": "a checks change that landed during the run"}
'''


class SnapshotHelpersTest(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def test_take_and_state(self):
        src = os.path.join(self.d, "src")
        os.makedirs(src)
        for name in sc.SNAPSHOT_FILES:
            shutil.copy(os.path.join(SUITE, name), src)
        run_dir = os.path.join(self.d, "runs", "R1")
        doc = sc.take_snapshot(src, run_dir, "abc", 2)
        self.assertEqual(sorted(doc["files"]), sorted(sc.SNAPSHOT_FILES))
        self.assertEqual(doc["files"]["checks.py"], sc.file_sha256(os.path.join(SUITE, "checks.py")))
        self.assertEqual((doc["scripts_commit"], doc["sessions_before"]), ("abc", 2))
        self.assertEqual(sc.snapshot_state(run_dir, src)[1:], ([], []))
        with open(os.path.join(src, "checks.py"), "a") as fh:
            fh.write(FAILING_CHECKS)
        self.assertEqual(sc.snapshot_state(run_dir, src)[1:], ([], ["checks.py"]))
        with open(os.path.join(run_dir, "scripts", "suite_common.py"), "a") as fh:
            fh.write("\n")
        self.assertEqual(sc.snapshot_state(run_dir, src)[1], ["suite_common.py"])
        self.assertEqual(sc.snapshot_state(os.path.join(self.d, "none")), (None, [], []))
        # the CLI takes it once, then keeps it
        other = os.path.join(self.d, "runs", "R2")
        out = subprocess.run([sys.executable, os.path.join(SUITE, "suite_common.py"), "snapshot", other, src, "c", "0"],
                             stdout=subprocess.PIPE, universal_newlines=True).stdout
        self.assertEqual(out, "created\n")
        out = subprocess.run([sys.executable, os.path.join(SUITE, "suite_common.py"), "snapshot", other, SUITE, "c",
                              "0"], stdout=subprocess.PIPE, universal_newlines=True).stdout
        self.assertEqual(out, "kept\nlive-differs checks.py\n")

    def test_snapshot_pins(self):
        self.assertTrue(sc.snapshot_pins("/a/runs/R1/scripts", "/b/copy/R1"))
        self.assertTrue(sc.snapshot_pins("/a/runs/R1/scripts/", "/a/runs/R1"))
        self.assertFalse(sc.snapshot_pins("/a/metagraph/scripts/suite", "/a/runs/R1"))
        self.assertFalse(sc.snapshot_pins("/a/runs/R2/scripts", "/a/runs/R1"))
        self.assertFalse(sc.snapshot_pins(None, "/a/runs/R1"))

    def test_run_direct_records_the_checks_it_imported(self):
        mpath = os.path.join(self.d, "m.json")
        make_manifest(mpath, basic_scenarios()[:1])
        out = os.path.join(self.d, "run")
        with FakeHost() as host:
            self.assertEqual(run_direct.main(["--manifest", mpath, "--arm", "direct-local", "--base", host.base,
                                              "--out", out, "--pause", "0", "--quiet", "--checks", "checks",
                                              "--warm-every", "0"]), 0)
        run = sc.read_json(os.path.join(out, "run.json"))
        want = (sc.file_sha256(os.path.join(SUITE, "checks.py")),
                sc.file_sha256(os.path.join(SUITE, "suite_common.py")), SUITE)
        sess = run["sessions"][0]
        self.assertEqual((sess["checks_sha256"], sess["suite_common_sha256"], sess["checks_dir"]), want)
        self.assertEqual((run["checks_sha256"], run["suite_common_sha256"], run["checks_dir"]), want)


class MergeChecksIdentityTest(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.run_dir = os.path.join(self.d, "runs", "RUN1")

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def make(self, snapshot=True, direct_sha="snapshot", mcp_dir="snapshot", started=sc.iso_utc(T0 + 3600),
             mcp_sha="snapshot"):
        """A run whose direct-local sessions state |direct_sha| ('snapshot': the snapshot's) and whose MCP arm
        records |mcp_dir| ('snapshot': the run's) and the checks_sha256 and suite_common_sha256 |mcp_sha|
        ('snapshot': the snapshot's; None: neither)."""
        doc = sc.take_snapshot(SUITE, self.run_dir, "c", 0, now=T0) if snapshot else {"files": {}}
        files = doc["files"]
        sha = files.get("checks.py") if direct_sha == "snapshot" else direct_sha
        sessions = [{"started_utc": started, "checks_sha256": sha,
                     "suite_common_sha256": files.get("suite_common.py") if direct_sha == "snapshot" else sha}]
        write_arm(self.run_dir, "direct-local", [line("P-1", "direct-local", 0, T0 + 10, 5.0)],
                  run_extra={"sessions": sessions})
        mdir = os.path.join(self.run_dir, "scripts") if mcp_dir == "snapshot" else mcp_dir
        extra = {"checks_dir": mdir} if mdir else {}
        if mcp_sha == "snapshot":
            extra.update(checks_sha256=files.get("checks.py"), suite_common_sha256=files.get("suite_common.py"))
        elif mcp_sha:
            extra.update(checks_sha256=mcp_sha, suite_common_sha256=mcp_sha)
        write_arm(self.run_dir, "mcp", [line("P-1", "mcp", 0, T0 + 20, 5.0)], run_extra=extra)
        return self.run_dir

    def strict(self, d):
        err = io.StringIO()
        with contextlib.redirect_stderr(err), contextlib.redirect_stdout(io.StringIO()):
            rc = merge_run.main([d, "--strict-identity"])
        return rc, err.getvalue()

    def test_a_run_judged_with_its_snapshot_is_quiet(self):
        merged = merge_run.merge(self.make())
        ident = merged["checks_identity"]
        self.assertEqual((ident["differs"], ident["unverified"], ident["errors"]), ([], [], []))
        self.assertEqual(ident["arms"]["mcp"]["checks_dirs"], [os.path.join(self.run_dir, "scripts")])
        self.assertFalse([w for w in merged["warnings"] if "checks" in w])
        self.assertNotIn("checks identity", self.strict(self.run_dir)[1])   # (the synthetic run has no manifest)

    def test_an_arm_judged_with_other_checks_is_an_error(self):
        # direct-local imported another checks.py; the MCP runner judged with the live tree, not the snapshot
        merged = merge_run.merge(self.make(direct_sha="0" * 64, mcp_dir="/somewhere/metagraph/scripts/suite",
                                           mcp_sha=None))
        errors = merged["checks_identity"]["errors"]
        self.assertTrue([e for e in errors if e.startswith("arms/direct-local: 1 of 1 session(s) judged with checks "
                                                           "sha256 000000000000 (sessions 0), the run's snapshot's")])
        self.assertTrue([e for e in errors if e.startswith("arms/mcp: 1 of 1 invocation(s) judged with other checks")])
        self.assertEqual(len([w for w in merged["warnings"] if "judged with" in w]), 3)
        rc, err = self.strict(self.run_dir)
        self.assertEqual(rc, 3)
        self.assertIn("checks identity: arms/mcp", err)

    def test_an_arm_that_states_nothing_is_unverified(self):
        # a session before the field, an MCP runner that records no checks: unverified, an error where the run has
        # a snapshot; without one (a run of before the snapshot) the unverified arms are named as before and the
        # run itself is unpinned, an error for the strict reading
        merged = merge_run.merge(self.make(direct_sha=None, mcp_dir=None, started=sc.iso_utc(T0 - 3600),
                                           mcp_sha=None))
        ident = merged["checks_identity"]
        self.assertEqual(ident["unverified"], [
            "arms/direct-local: 1 of 1 session(s) state no checks sha256 (1 ran before the run's snapshot was taken)",
            "arms/mcp: 1 of 1 invocation(s) do not state both checks_sha256 and suite_common_sha256 (checks_sha256 "
            "missing on 1, suite_common_sha256 missing on 1; a checks directory does not stand in for them)"])
        self.assertEqual(ident["errors"], ident["unverified"])
        rc, err = self.strict(self.run_dir)
        self.assertEqual(rc, 3)
        self.assertIn("checks identity: arms/mcp: 1 of 1 invocation(s) do not state both", err)
        other = os.path.join(self.d, "old")
        self.run_dir = other
        merged = merge_run.merge(self.make(snapshot=False, direct_sha=None, mcp_dir=None, mcp_sha=None))
        ident = merged["checks_identity"]
        self.assertEqual((ident["snapshot"], ident["unverified"][1:], ident["errors"]), (
            None, ["arms/mcp: 1 of 1 invocation(s) state no checks directory or sha256"], [merge_run.UNPINNED]))
        self.assertFalse([w for w in merged["warnings"] if "snapshot" in w])   # readable as before: no warning
        rc, err = self.strict(other)
        self.assertEqual(rc, 3)
        self.assertIn("checks identity: the run has no snapshot (scripts/SNAPSHOT.json)", err)

    def test_with_a_snapshot_every_mcp_invocation_states_both_sha256(self):
        # the checks directory (the run's snapshot) never stands in for a sha256: an invocation missing either one is
        # unverified, an error; a stated one other than the snapshot's differs; the stub's for a stub-judged arm
        files = sc.read_json(os.path.join(self.make(), "scripts", "SNAPSHOT.json"))["files"]
        snap = os.path.join(self.run_dir, "scripts")
        full = {"checks_dir": snap, "checks_sha256": files["checks.py"],
                "suite_common_sha256": files["suite_common.py"]}
        path = os.path.join(self.run_dir, "arms", "mcp", "run.json")
        run = sc.read_json(path)

        def errors(invocations, **top):
            sc.write_json_atomic(path, dict(run, invocations=invocations, **top))
            return merge_run.merge(self.run_dir)["checks_identity"]["errors"]
        self.assertEqual(errors([full, full]), [])
        for drop in (["suite_common_sha256"], ["checks_sha256"], ["checks_sha256", "suite_common_sha256"]):
            got = errors([full, dict((k, v) for k, v in full.items() if k not in drop)])
            self.assertEqual(len(got), 1, drop)
            self.assertTrue(got[0].startswith("arms/mcp: 1 of 2 invocation(s) do not state both checks_sha256 and "
                                              "suite_common_sha256 (%s missing on 1" % drop[0]), got)
        got = errors([full, dict(full, suite_common_sha256="e" * 64)])
        self.assertTrue(got[0].startswith("arms/mcp: 1 of 2 invocation(s) judged with other checks"), got)
        stub = dict(full, checks_sha256=files["checks_stub.py"])
        self.assertTrue(errors([stub])[0].startswith("arms/mcp: 1 of 1 invocation(s) judged with other checks"))
        self.assertEqual(errors([stub], checks_module={"file": "checks_stub.py", "stub": True}), [])

    def test_without_a_snapshot_arms_stating_different_checks_are_named(self):
        merged = merge_run.merge(self.make(snapshot=False, direct_sha="a" * 64, mcp_dir=None))
        self.assertEqual(merged["checks_identity"]["differs"], [])
        run = sc.read_json(os.path.join(self.run_dir, "arms", "mcp", "run.json"))
        run["checks_sha256"] = "b" * 64        # a runner without invocations states its one at the top
        sc.write_json_atomic(os.path.join(self.run_dir, "arms", "mcp", "run.json"), run)
        merged = merge_run.merge(self.run_dir)
        self.assertTrue(merged["checks_identity"]["differs"][0].startswith(
            "the arms were judged with different checks"))

    def test_the_mcp_runner_states_its_checks_per_invocation(self):
        # async's runner records checks_dir, checks_sha256 and suite_common_sha256 per invocation (run.json
        # invocations[]); an invocation whose suite_common is not the snapshot's is named, and a top-level
        # checks_module is not read once invocations exist
        d = self.make()
        files = sc.read_json(os.path.join(d, "scripts", "SNAPSHOT.json"))["files"]
        snap = os.path.join(d, "scripts")
        good = {"checks_dir": snap, "checks_sha256": files["checks.py"],
                "suite_common_sha256": files["suite_common.py"]}
        run = sc.read_json(os.path.join(d, "arms", "mcp", "run.json"))
        run.update(invocations=[good, dict(good, suite_common_sha256="c" * 64)],
                   checks_module={"file": "x", "stub": False, "sha256": "d" * 64})
        sc.write_json_atomic(os.path.join(d, "arms", "mcp", "run.json"), run)
        ident = merge_run.merge(d)["checks_identity"]
        self.assertEqual(ident["errors"], ["arms/mcp: 1 of 2 invocation(s) judged with other checks than the run's "
                                           "snapshot (e.g. checks_sha256 %s, suite_common_sha256 %s, checks_dir %s)"
                                           % (files["checks.py"], "c" * 64, snap)])
        self.assertEqual(ident["arms"]["mcp"]["suite_common_sha256"], sorted([files["suite_common.py"], "c" * 64]))
        run["invocations"] = [good, good]
        sc.write_json_atomic(os.path.join(d, "arms", "mcp", "run.json"), run)
        self.assertEqual(merge_run.merge(d)["checks_identity"]["errors"], [])

    def test_the_merge_judges_with_the_snapshot(self):
        # checks.scenario_outcome makes outcomes stricter: the run's snapshot's, not the live tree's
        d = self.make()
        with open(os.path.join(d, "scripts", "checks.py"), "a") as fh:
            fh.write(FAILING_CHECKS)
        merged = merge_run.merge(d)
        rows, _ = sc.read_jsonl(os.path.join(d, "executions.jsonl"))
        self.assertEqual(set(r["outcome"] for r in rows), {"fail"})
        self.assertEqual(merged["checks_identity"]["differs"],
                         ["the run's snapshot changed after it was taken: scripts/checks.py"])


class ComparisonSnapshotTest(unittest.TestCase):
    """compare, the gate's and the report's property reading judge a recorded set with the run's snapshot (its
    checks.py and table) when it has one, else with this tree's, and say which."""

    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.run_dir = os.path.join(self.d, "runs", "RUN1")

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def test_the_checks_and_the_table_of_a_recorded_set(self):
        os.makedirs(self.run_dir)
        live = cmp.judging_checks(self.run_dir)
        self.assertEqual((live["source"], live["dir"]), ("live", SUITE))
        self.assertEqual(cmp.table_for(self.run_dir)[0], cmp.DEFAULT_BY_DESIGN)
        sc.take_snapshot(SUITE, self.run_dir, "c", 0)
        snap = os.path.join(self.run_dir, "scripts")
        sc_sha = sc.file_sha256(os.path.join(SUITE, "suite_common.py"))
        self.assertEqual(cmp.judging_checks(self.run_dir), {
            "source": "snapshot", "dir": snap, "checks_sha256": sc.file_sha256(os.path.join(SUITE, "checks.py")),
            "snapshot_changed": [], "suite_common_sha256": sc_sha, "pipeline_suite_common_sha256": sc_sha,
            "live_suite_common_differs": False})
        self.assertEqual(cmp.table_for(self.run_dir)[0], os.path.join(snap, "service_by_design.json"))
        self.assertEqual(cmp.table_for(self.run_dir, "")[0], "")             # a table given (none) still wins
        # a re-evaluated set is judged with the re-evaluation's checks and table, never the snapshot's
        self.assertEqual(cmp.judging_checks(self.run_dir, reevaluated=True)["source"], "live")
        self.assertEqual(cmp.table_for(self.run_dir, None, True)[0], cmp.DEFAULT_BY_DESIGN)
        with open(os.path.join(snap, "checks.py"), "a") as fh:
            fh.write(FAILING_RELATIONS)
        self.assertEqual(cmp.judging_checks(self.run_dir)["snapshot_changed"], ["checks.py"])

    def test_the_gate_and_the_report_read_properties_with_the_snapshot(self):
        m = manifest_of()
        d = make_run(self.run_dir, m, mcp_lines(m))
        sc.take_snapshot(SUITE, d, "c", 0)
        seen = []
        real = gate._is_property_fn

        def spy(checks_dir=None):
            seen.append(checks_dir)
            return real(checks_dir)
        gate._is_property_fn = spy
        try:
            gate.judge(d, policy_of(m))
            os.remove(os.path.join(d, "compare.json"))
            text = report.build_report(d)
        finally:
            gate._is_property_fn = real
        self.assertEqual(seen, [os.path.join(d, "scripts")] * 2)
        self.assertIn("checks.py of the run's snapshot (scripts/, sha256 %s, with the snapshot's suite_common.py) |"
                      % sc.file_sha256(os.path.join(SUITE, "checks.py"))[:12], text)
        self.assertIn("the run's snapshot (scripts/service_by_design.json)", text)
        comp = sc.read_json(os.path.join(d, "compare.json"))
        self.assertEqual((comp["judged_with"]["source"], comp["judged_with"]["by_design"]), ("snapshot", "snapshot"))
        # a run without a snapshot says the tree's checks judged it
        shutil.rmtree(os.path.join(d, "scripts"))
        os.remove(os.path.join(d, "compare.json"))
        self.assertIn("checks.py of this tree (%s" % SUITE.replace("_", "\\_"), report.build_report(d))


# a checks.py whose relation reads a suite_common helper (answer_digest): a relation rule a later checks.py may have
DIGEST_RELATION = '''

def evaluate_relation(rel, a_scenario, a_answer, b_scenario, b_answer, arm):
    a, b = answer_digest(a_answer), answer_digest(b_answer)
    return {"relation": rel.get("id"), "type": rel.get("type"), "a": rel.get("a"), "b": rel.get("b"),
            "outcome": "pass" if a == b else "fail", "detail": "digests %s / %s" % (a[:24], b[:24])}
'''
# suite_common.py as a change landing in the tree after the snapshot would leave it: every digest a new sentinel
SENTINEL_DIGEST = '''

_SENTINELS = []


def answer_digest(answer, route=None):
    _SENTINELS.append(1)
    return "SENTINEL-LIVE-SUITE-COMMON-%d" % len(_SENTINELS)
'''
# the verdict pipeline of the live tree LIVE judging the recorded run RUN_DIR (a fresh interpreter): the module the
# snapshot's checks.py bound, the relations, compare.json's judged_with and the report's verdict and note
JUDGE = r'''
import json, os, sys
live, run_dir = sys.argv[1], sys.argv[2]
sys.path.insert(0, live)
import compare as cmp, report, suite_common as sc
j = cmp.judging_checks(run_dir)
mod = cmp.load_checks_module(j["dir"])
table, _ = cmp.table_for(run_dir)
comp = cmp.compare(run_dir, cmp.load_by_design(table))
sc.write_json_atomic(os.path.join(run_dir, "compare.json"), comp)
text = report.build_report(run_dir)
print(json.dumps({
    "bound": os.path.relpath(mod._suite_common.__file__, run_dir),
    "pipeline": os.path.relpath(sys.modules["suite_common"].__file__, live),
    "relations": dict((r["id"], r["arms"]["mcp"]["outcome"]) for r in comp["relations"]),
    "verdict": text.split("# Verdict", 1)[1].strip().splitlines()[0],
    "judged_with": comp["judged_with"],
    "note": [x for x in text.splitlines() if "suite_common.py of this tree" in x]}))
'''


class MovedRunTest(unittest.TestCase):
    """compare.json names the run's snapshot (its table, its checks) relative to RUN_DIR: a report regenerated after
    the run directory moved reads the same table and states the same."""

    CASE = "capabilities.traverse_capabilities_document"

    def setUp(self):
        self.d = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def test_a_moved_copy_re_reports_identically(self):
        # C-caps-0001 inconclusive via MCP by a limit the table states (the snapshot's table, read through compare.json)
        caps = {"id": "C-caps-0001", "area": "C", "route": "GET /traverse/capabilities", "request": None,
                "mcp": {"tool": "traverse_capabilities", "arguments": {}, "by_design": [self.CASE]},
                "expect": [{"check": "status", "args": {"http": 200}}, {"check": "capabilities_doc", "args": {}}]}
        S = [dict(caps, id="P-001", area="P", route="POST /pattern", mcp={"tool": "pattern_search", "arguments": {}},
                  expect=[{"check": "status", "args": {"http": 200}}], title="p"), caps]
        manifest = {"manifest_version": 1, "content_sha256": sc.scenarios_sha256(S), "scenarios": S}
        run_dir = os.path.join(self.d, "runs", "RUN1")
        status = {"check": "status", "outcome": "pass", "detail": "http 200"}
        mcp = [line("P-001", "mcp", r, T0 + 10 + r, 20.0, checks=[status]) for r in (0, 1)]
        for r in (0, 1):
            ln = line("C-caps-0001", "mcp", r, T0 + 20 + r, 20.0, area="C", outcome="inconclusive",
                      checks=[status, {"check": "capabilities_doc", "outcome": "not_exposed",
                                       "detail": "field not in the answer: release"}])
            ln.update(service_by_design_case=self.CASE, explained_checks=[], explained_items=[])
            mcp.append(ln)
        rows = {"direct-local": [line(x["id"], "direct-local", 0, T0 + i, 20.0, area=x["area"],
                                      checks=[dict(status, check=e["check"]) for e in x["expect"]])
                                 for i, x in enumerate(S)], "mcp": mcp}
        for arm, lines in rows.items():
            os.makedirs(os.path.join(run_dir, "arms", arm), exist_ok=True)
            sc.write_json_atomic(os.path.join(run_dir, "arms", arm, "manifest.json"), manifest)
            write_arm(run_dir, arm, lines)
        sc.write_json_atomic(os.path.join(run_dir, "orchestrator.json"), {"sessions": [{
            "mode": "paired", "mcp_cmd": "x {ids} {out}", "checks": "checks", "warm_every": "0", "skip_remote": "1",
            "loop": {"mode": "paired", "state": "done", "scenarios": 2, "completed": 2}}]})
        sc.take_snapshot(SUITE, run_dir, "c", 0)
        merge_run.merge(run_dir)
        comp = cmp.compare(run_dir, cmp.load_by_design(cmp.table_for(run_dir)[0]))
        sc.write_json_atomic(os.path.join(run_dir, "compare.json"), comp)

        def text_of(d):
            return [x for x in report.build_report(d).splitlines() if not x.startswith(("date:", "| generated"))]
        before = text_of(run_dir)
        self.assertIn("C-caps-0001 mcp r0,1 (capabilities.traverse\\_capabilities\\_document)", "\n".join(before))
        moved = os.path.join(self.d, "moved", "RUN1")
        os.makedirs(os.path.dirname(moved))
        os.rename(run_dir, moved)
        self.assertEqual(text_of(moved), before)
        self.assertEqual((comp["by_design_table"]["path"], comp["judged_with"]["dir"]),
                         (os.path.join("scripts", "service_by_design.json"), "scripts"))


class SnapshotSuiteCommonTest(unittest.TestCase):
    """The snapshot's checks.py imports the snapshot's suite_common in the verdict pipeline (compare, the gate's and
    the report's property reading, merge_run's stricter outcomes), never this tree's: a change of a suite_common helper
    the checks call, landing in the tree after the run, cannot move the recorded run's verdict; that the pipeline's
    own suite_common differs from the snapshot's is recorded and noted."""

    def setUp(self):
        self.d = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def judge(self, live, run_dir):
        p = subprocess.run([sys.executable, "-I", "-B", "-c", JUDGE, live, run_dir], stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE, universal_newlines=True, timeout=120)
        self.assertEqual(p.returncode, 0, p.stderr)
        return json.loads(p.stdout.strip().splitlines()[-1])

    def test_a_live_suite_common_change_cannot_move_the_verdict(self):
        live = os.path.join(self.d, "live")
        os.makedirs(live)
        for name in os.listdir(SUITE):
            if name.endswith((".py", ".json")):
                shutil.copy2(os.path.join(SUITE, name), live)
        with open(os.path.join(live, "checks.py"), "a") as fh:
            fh.write(DIGEST_RELATION)
        m = manifest_of()
        run_dir = make_run(os.path.join(self.d, "runs", "RUN1"), m, mcp_lines(m))
        sc.take_snapshot(live, run_dir, "c", 0)
        before = self.judge(live, run_dir)
        self.assertEqual(before["relations"], {"REL-1": "pass", "REL-2": "pass"})
        with open(os.path.join(live, "suite_common.py"), "a") as fh:
            fh.write(SENTINEL_DIGEST)
        after = self.judge(live, run_dir)
        self.assertEqual((after["relations"], after["verdict"]), (before["relations"], before["verdict"]))
        # the snapshot's checks bound the snapshot's suite_common; the pipeline keeps its own
        self.assertEqual(before["bound"], os.path.join("scripts", "suite_common.py"))
        self.assertEqual((after["bound"], after["pipeline"]), (os.path.join("scripts", "suite_common.py"),
                                                               "suite_common.py"))
        self.assertFalse(before["judged_with"]["live_suite_common_differs"])
        self.assertEqual(before["note"], [])
        jw = after["judged_with"]
        self.assertTrue(jw["live_suite_common_differs"])
        self.assertEqual(jw["suite_common_sha256"], sc.file_sha256(os.path.join(run_dir, "scripts", "suite_common.py")))
        self.assertEqual(jw["pipeline_suite_common_sha256"], sc.file_sha256(os.path.join(live, "suite_common.py")))
        self.assertEqual(len(after["note"]), 1, after["note"])
        self.assertIn("suite_common.py of this tree differs from the run's snapshot", after["note"][0])

    def test_the_gate_notes_a_differing_pipeline_suite_common(self):
        # a note of the gate, never a problem: the checks judged with the snapshot's suite_common
        m = manifest_of()
        d = make_run(os.path.join(self.d, "runs", "RUN1"), m, mcp_lines(m))
        sc.take_snapshot(SUITE, d, "c", 0)
        comp = cmp.compare(d, cmp.load_by_design(cmp.table_for(d)[0]))
        comp["judged_with"].update(live_suite_common_differs=True, pipeline_suite_common_sha256="a" * 64)
        sc.write_json_atomic(os.path.join(d, "compare.json"), comp)
        v = gate.judge(d, policy_of(m))
        self.assertIn("suite_common.py of this tree (sha256 aaaaaaaaaaaa) differs from the one the checks judged with "
                      "(snapshot, sha256 %s)" % sc.file_sha256(os.path.join(SUITE, "suite_common.py"))[:12], v["notes"])
        self.assertFalse([p for p in v["problems"] if "suite_common" in p["detail"]])

    def test_the_bound_loader_leaves_the_callers_suite_common(self):
        snap = os.path.join(self.d, "snap")
        os.makedirs(snap)
        for name in ("checks.py", "suite_common.py"):
            shutil.copy2(os.path.join(SUITE, name), snap)
        with open(os.path.join(snap, "suite_common.py"), "a") as fh:
            fh.write(SENTINEL_DIGEST)
        path = list(sys.path)
        mod = sc.load_checks_bound(snap, "suite_checks_probe")
        self.assertTrue(mod.answer_digest({"http": 200, "body": {}}).startswith("SENTINEL-LIVE-SUITE-COMMON"))
        self.assertIs(sys.modules["suite_common"], sc)
        self.assertEqual(sys.path, path)
        for fn in (cmp.load_checks_module(snap), merge_run.load_scenario_outcome(snap),
                   cmp.load_relation_evaluator(snap)):
            self.assertIsNotNone(fn)
        self.assertTrue(cmp.load_checks_module(snap).answer_digest({"http": 200}).startswith("SENTINEL"))
        self.assertIs(sys.modules["suite_common"], sc)


class RunSuiteSnapshotTest(SuiteCase):
    def test_a_checks_change_during_the_run_does_not_reach_the_comparison(self):
        # a copy of the suite runs a paired run with one same_answer relation; after the first request the copy's
        # checks.py gets an evaluate_relation failing every relation: the comparison judges the relations with the
        # run's snapshot and its table, so the verdict stays PASS (before, compare loaded the tree's checks.py and
        # read FAIL, "1 failed relation" per arm), and the report says what judged it
        suite = os.path.join(self.d, "suite")
        os.makedirs(suite)
        for name in os.listdir(SUITE):
            if name.endswith((".py", ".sh", ".json")):
                shutil.copy2(os.path.join(SUITE, name), suite)
        m = make_manifest(self.mpath, [dict(s, mcp={"tool": "pattern_search", "arguments": {}}) for s in
                                       [scenario("S-%02d" % i, "S", "POST /slow", {"sleep_s": 0.2}) for i in range(4)]])
        m["relations"] = [{"id": "REL-1", "type": "same_answer", "a": "S-00", "b": "S-01"}]
        with open(self.mpath, "w") as fh:
            json.dump(m, fh)
        cmd = [os.path.join(suite, "run_suite.sh"), "--no-pdf", "--pause", "0", "--manifest", self.mpath,
               "--warm-every", "0", "--checks", "checks", "--mcp-cmd", self.mcp_cmd()]
        with FakeHost() as host:
            proc = subprocess.Popen(cmd, env=self.env(host), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    universal_newlines=True, stdin=subprocess.DEVNULL)
            deadline = time.time() + 60
            while time.time() < deadline and sum(1 for _, x in host.log if x == "/slow") < 1:
                time.sleep(0.02)
            with open(os.path.join(suite, "checks.py"), "a") as fh:
                fh.write(FAILING_RELATIONS)
            out, err = proc.communicate(timeout=180)
        self.assertEqual(proc.returncode, 0, err)
        run_dir = out.strip().splitlines()[-1]
        comp = sc.read_json(os.path.join(run_dir, "compare.json"))
        self.assertEqual([(r["id"], dict((arm, v["outcome"]) for arm, v in r["arms"].items()))
                          for r in comp["relations"]],
                         [("REL-1", {"direct-local": "pass", "direct-remote": "pass", "mcp": "not_evaluated"})])
        self.assertEqual((comp["judged_with"]["source"], comp["judged_with"]["dir"], comp["judged_with"]["by_design"]),
                         ("snapshot", "scripts", "snapshot"))
        self.assertEqual(comp["by_design_table"]["path"], os.path.join("scripts", "service_by_design.json"))
        with open(os.path.join(run_dir, "report.md")) as fh:
            text = fh.read()
        self.assertIn("**PASS", text)
        self.assertIn("checks.py of the run's snapshot", text)

    def test_a_checks_change_during_the_run_does_not_reach_it(self):
        # a copy of the suite runs a paired run; a checks change lands in the copy after the first scenario: every
        # later step still judges with the run's snapshot (before it, direct-remote's later steps imported the
        # changed checks.py), each arm records the snapshot's sha256 and the merge finds nothing to state; a resume
        # says the tree differs and keeps the snapshot
        suite = os.path.join(self.d, "suite")
        os.makedirs(suite)
        for name in os.listdir(SUITE):
            if name.endswith((".py", ".sh", ".json")):
                shutil.copy2(os.path.join(SUITE, name), suite)
        make_manifest(self.mpath, [dict(s, mcp={"tool": "pattern_search", "arguments": {}}) for s in
                                   [scenario("S-%02d" % i, "S", "POST /slow", {"sleep_s": 0.2}) for i in range(4)]])
        cmd = [os.path.join(suite, "run_suite.sh"), "--no-pdf", "--pause", "0", "--manifest", self.mpath,
               "--warm-every", "0", "--checks", "checks", "--mcp-cmd", self.mcp_cmd()]
        with FakeHost() as host:
            proc = subprocess.Popen(cmd, env=self.env(host), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    universal_newlines=True, stdin=subprocess.DEVNULL)
            deadline = time.time() + 60
            while time.time() < deadline and sum(1 for _, x in host.log if x == "/slow") < 2:
                time.sleep(0.02)
            with open(os.path.join(suite, "checks.py"), "a") as fh:
                fh.write(FAILING_CHECKS)
            out, err = proc.communicate(timeout=180)
            self.assertEqual(proc.returncode, 0, err)
            run_dir = out.strip().splitlines()[-1]
            rows, _ = sc.read_jsonl(os.path.join(run_dir, "executions.jsonl"))
            self.assertEqual(len(rows), 4 * 4)
            self.assertEqual(set(r["outcome"] for r in rows), {"pass"})
            snap = sc.read_json(os.path.join(run_dir, "scripts", "SNAPSHOT.json"))
            self.assertEqual(snap["files"]["checks.py"], sc.file_sha256(os.path.join(SUITE, "checks.py")))
            merged = sc.read_json(os.path.join(run_dir, "run.json"))
            for arm in ("direct-local", "direct-remote"):
                self.assertEqual(set(x["checks_sha256"] for x in merged["arms"][arm]["sessions"]),
                                 {snap["files"]["checks.py"]}, arm)
            self.assertEqual(merged["checks_identity"]["errors"], [])
            self.assertEqual(merged["checks_identity"]["arms"]["mcp"]["checks_dirs"],
                             [os.path.join(run_dir, "scripts")])
            self.assertEqual(merged["checks_identity"]["arms"]["mcp"]["checks_sha256"], [snap["files"]["checks.py"]])
            sess = sc.read_json(os.path.join(run_dir, "orchestrator.json"))["sessions"][-1]
            self.assertEqual(sess["snapshot"], {"state": "created", "changed": [], "live_differs": [],
                                                "mcp_checks_dir": "snapshot"})
            p = subprocess.run([os.path.join(suite, "run_suite.sh"), "--no-pdf", "--resume", run_dir],
                               env=self.env(host), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               universal_newlines=True, stdin=subprocess.DEVNULL, timeout=180)
            self.assertEqual(p.returncode, 0, p.stderr)
        self.assertIn("checks.py of %s differ(s) from the run's snapshot" % suite, p.stderr)
        sess = sc.read_json(os.path.join(run_dir, "orchestrator.json"))["sessions"][-1]
        self.assertEqual((sess["snapshot"]["state"], sess["snapshot"]["live_differs"]), ("kept", ["checks.py"]))
        self.assertEqual(sc.read_json(os.path.join(run_dir, "run.json"))["checks_identity"]["errors"], [])


if __name__ == "__main__":
    unittest.main()
