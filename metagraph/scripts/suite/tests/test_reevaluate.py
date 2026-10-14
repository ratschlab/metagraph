"""Offline tests of reevaluate.py and the --reevaluated set of merge_run.py, compare.py, report.py and run_suite.sh
on a synthetic recorded run: what is judged again and what is kept, what the checks are passed (previous, client),
the mark of the set, and that every recorded file stays as it was."""

import hashlib
import json
import os
import subprocess
import sys
import tempfile
import textwrap
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.dirname(HERE)
sys.path.insert(0, SUITE)

import compare as cmp  # noqa: E402
import merge_run  # noqa: E402
import reevaluate  # noqa: E402
import report  # noqa: E402
import suite_common as sc  # noqa: E402

T0 = 1791800000.0
BOUND = ("the attempt reached the duration bound the server enforces for it (15000 ms: the seeds' time budgets "
         "plus the server's allowance, at most hard_cap_ms, the content timeout less one second; see "
         "usage.bound) while its response was built or written: nothing of it is delivered")
WALKED = {"results": [{"outcome": {"walks": "complete"},
                       "arms": {"left": {"status": "complete", "segments": [], "paths": []},
                                "right": {"status": "complete", "segments": [], "paths": []}}}]}


def scenario(sid, route, expect, mcp=None, request=None):
    return {"id": sid, "area": sid.split("-")[0], "title": "title of %s" % sid, "route": route,
            "request": request if request is not None else ({"seeds": [{"sequence": "ACGT" * 10}]}
                                                            if route == "POST /traverse" else {}),
            "mcp": mcp or {"tool": "x", "arguments": {}}, "truth": {}, "expect": expect, "timing_class": "fast",
            "repeat": 1, "tags": []}


STATUS = [{"check": "status", "args": {"http": 200}}]
TRAVERSE = STATUS + [{"check": "invariants", "args": {}},
                     {"check": "traverse_seed_outcome", "args": {"walks": ["complete", "partial"]}}]
SCENARIOS = [
    scenario("T-0001", "POST /traverse", TRAVERSE),
    scenario("P-0001", "POST /pattern", STATUS),
    scenario("P-0002", "POST /pattern", STATUS, mcp={"not_implemented": "max_steps"}),
    scenario("P-0003", "POST /pattern", STATUS, mcp={"tool": "x", "arguments": {}, "by_design": "case.x"}),
    scenario("P-0004", "POST /pattern", STATUS),
    scenario("P-0005", "POST /pattern", STATUS),
    scenario("P-0006", "POST /pattern", STATUS),
    scenario("R-0001", "POST /resolve", [{"check": "status", "args": {"http": [200, 503]}},
                                         {"check": "budget_respected", "args": {}}],
             request={"sequence": "ACGT" * 10, "bounds": {"time_budget_ms": 3000}}),
    scenario("D-0001", "POST /pattern", STATUS + [{"check": "deterministic", "args": {}}]),
]
TIME_LIMITED = {"patterns": [{"id": "w", "determinism": "time_limited", "stop": {"phase": "discovery", "reason": "time"}}],
                "limits": {"time_budget_ms": 60000, "finalize_reserve_ms": 250}, "timing": {"elapsed_ms": 59900.0}}
COMPLETE = {"patterns": [{"id": "w", "determinism": "full"}], "limits": {"time_budget_ms": 60000}}
DEADLINE = {"http": 503, "body": {"code": "deadline", "error": "late"}, "error": {"code": "deadline", "text": "late"}}

# (scenario, arm, repeat): (recorded outcome, recorded checks, answer or None, line extras)
OLD_FAIL = [{"check": "status", "outcome": "fail", "detail": "old checks"}]
EXECS = {
    ("T-0001", "direct-local", 0): ("pass", [{"check": "status", "outcome": "pass", "detail": "http 200"},
                                            {"check": "invariants", "outcome": "pass", "detail": "ok"},
                                            {"check": "traverse_seed_outcome", "outcome": "pass", "detail": "ok"}],
                                    {"http": 200, "body": WALKED, "error": None}, {}),
    ("T-0001", "mcp", 0): ("fail", OLD_FAIL, {"http": 200, "body": {"results": [{"error": BOUND}]}, "error": None},
                           {"mcp_phases": {"backend_ms": 15305.3}}),
    ("T-0001", "mcp", 1): ("fail", OLD_FAIL, {"http": 200, "body": {"results": [{"error": BOUND}]}, "error": None},
                           {"mcp_phases": {"backend_ms": 9000.0}}),
    ("P-0001", "direct-local", 0): ("fail", OLD_FAIL, {"http": 200, "body": COMPLETE, "error": None}, {}),
    ("P-0002", "mcp", 0): ("not_implemented_mcp", [], None, {}),
    ("P-0003", "mcp", 0): ("service_by_design", [], {"http": 400, "body": None,
                                                     "error": {"code": "pattern_invalid", "text": "x"}}, {}),
    ("P-0004", "direct-remote", 0): ("error", [], {"http": 0, "body": None,
                                                   "error": {"code": "transport", "text": "refused"}}, {}),
    ("P-0005", "direct-local", 0): ("error", [{"check": "evaluate", "outcome": "fail", "detail": "exception in "
                                               "checks.evaluate: KeyError"}],
                                    {"http": 200, "body": COMPLETE, "error": None}, {}),
    ("P-0006", "mcp", 0): ("inconclusive", [], {"http": 503, "body": None, "error": {
        "code": "uninterruptible_overrun_exceeds_allowance", "text": "R20 lockout"}}, {"mcp_phases": {}}),
    ("R-0001", "direct-local", 0): ("inconclusive", [], DEADLINE, {"wall_ms": 9257.0}),
    ("R-0001", "mcp", 0): ("inconclusive", [], DEADLINE, {"mcp_phases": {"backend_ms": 4000.0}}),
    ("D-0001", "mcp", 0): ("pass", None, {"http": 200, "body": TIME_LIMITED, "error": None}, {}),
    ("D-0001", "mcp", 1): ("fail", [{"check": "status", "outcome": "pass", "detail": "http 200"},
                                    {"check": "deterministic", "outcome": "fail", "detail": "digest differs"}],
                           {"http": 200, "body": COMPLETE, "error": None}, {}),
}


def build_run(run_dir):
    manifest = {"manifest_version": 1, "scenarios": SCENARIOS, "relations": []}
    os.makedirs(run_dir)
    sc.write_json_atomic(os.path.join(run_dir, "manifest.json"), manifest)
    t = T0
    for (sid, arm, rep), (outcome, checks, answer, extra) in sorted(EXECS.items(), key=lambda kv: (kv[0][1], kv[0][0], kv[0][2])):
        d = os.path.join(run_dir, "arms", arm)
        os.makedirs(os.path.join(d, "raw"), exist_ok=True)
        raw = sc.raw_name(sid, arm, rep)
        line = {"run_id": "R", "scenario": sid, "area": sid.split("-")[0], "arm": arm, "repeat": rep,
                "cache": "cold" if rep == 0 else "warm", "started_utc": sc.iso_utc(t), "ended_utc": sc.iso_utc(t + 1),
                "wall_ms": 1000.0, "http": (answer or {}).get("http"), "outcome": outcome,
                "checks": checks if checks is not None else [{"check": "status", "outcome": "pass", "detail": "ok"}],
                "server": None, "mcp_phases": None, "answer_digest": "recorded-%s-%d" % (sid, rep), "raw": raw,
                "job": None, "error": (answer or {}).get("error"),
                "route": next(s["route"] for s in SCENARIOS if s["id"] == sid)}
        line.update(extra)
        t += 2
        sc.append_jsonl(os.path.join(d, "executions.jsonl"), line)
        doc = {"scenario": sid, "arm": arm, "repeat": rep}
        if answer is not None:
            doc["answer"] = answer
        else:
            doc["not_implemented"] = "max_steps"
        sc.write_gz_json(os.path.join(d, raw), doc)
    for arm in ("direct-local", "direct-remote", "mcp"):
        d = os.path.join(run_dir, "arms", arm)
        run = {"run_id": "R", "arm": arm, "state": "done", "release": "refseq97-x", "scripts_commit": "6f945db4",
               "checks_module": {"file": "/x/checks.py", "stub": False} if arm == "mcp" else "checks",
               "manifest": {"computed_sha256": "sha"}}
        sc.write_json_atomic(os.path.join(d, "run.json"), run)
        sc.write_json_atomic(os.path.join(d, "manifest.json"), manifest)


def digests(root):
    out = {}
    for d, _, files in os.walk(root):
        for f in files:
            p = os.path.join(d, f)
            with open(p, "rb") as fh:
                out[os.path.relpath(p, root)] = hashlib.sha256(fh.read()).hexdigest()
    return out


class Reevaluate(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.run_dir = os.path.join(self.tmp.name, "run")
        build_run(self.run_dir)

    def tearDown(self):
        self.tmp.cleanup()

    def lines(self, arm):
        rows, bad = sc.read_jsonl(os.path.join(self.run_dir, "arms", arm, "executions.reevaluated.jsonl"))
        self.assertEqual(bad, 0)
        return dict(((r["scenario"], r["repeat"]), r) for r in rows)

    def test_reevaluated_outcomes(self):
        mark = reevaluate.reevaluate(self.run_dir)
        mcp, local, remote = self.lines("mcp"), self.lines("direct-local"), self.lines("direct-remote")
        # the attempt bound of the MCP arm: reached (backend 15,305 ms) -> inconclusive; not reached -> fail
        self.assertEqual(mcp[("T-0001", 0)]["outcome"], "inconclusive")
        self.assertEqual(mcp[("T-0001", 1)]["outcome"], "fail")
        self.assertIn("stated as reached after 9000 ms",
                      " ".join(c["detail"] for c in mcp[("T-0001", 1)]["checks"]))
        # a stale recorded check outcome is judged again; an evaluate exception is judged again
        self.assertEqual(local[("P-0001", 0)]["outcome"], "pass")
        self.assertEqual(local[("P-0005", 0)]["outcome"], "pass")
        # kept: not implemented, a transport error, the MCP runner's host-state rule; a runner's by-design mark the
        # current table does not reproduce (case.x is not in it) is judged by the checks (B-06), and says so
        self.assertEqual(mcp[("P-0002", 0)]["outcome"], "not_implemented_mcp")
        self.assertEqual(mcp[("P-0002", 0)]["reevaluation"]["kept"], "no recorded answer")
        self.assertEqual(mcp[("P-0003", 0)]["outcome"], "fail")
        self.assertIn("not reproduced by the current table", mcp[("P-0003", 0)]["reevaluation"]["note"])
        self.assertNotIn("service_by_design_case", mcp[("P-0003", 0)])
        self.assertEqual(remote[("P-0004", 0)]["outcome"], "error")
        self.assertEqual(mcp[("P-0006", 0)]["outcome"], "inconclusive")
        self.assertIn("host-state", mcp[("P-0006", 0)]["reevaluation"]["kept"])
        # the client: a direct arm's wall (9,257 ms on a 3,000 ms budget: a finding), the MCP arm's backend time
        self.assertEqual(local[("R-0001", 0)]["outcome"], "finding")
        self.assertEqual(mcp[("R-0001", 0)]["outcome"], "inconclusive")
        # previous as {digest, time_limited}: the time-limited cold repeat explains the warm one's difference
        self.assertEqual(mcp[("D-0001", 1)]["outcome"], "inconclusive")
        det = [c for c in mcp[("D-0001", 1)]["checks"] if c["check"] == "deterministic"][0]
        self.assertIn("no earlier repeat without a time limitation", det["detail"])
        self.assertEqual(mcp[("D-0001", 0)]["time_limited"], "pattern 0 time-limited")
        # the digest by these checks, the recorded one kept beside it
        self.assertEqual(mcp[("D-0001", 0)]["answer_digest"],
                         sc.answer_digest({"http": 200, "body": TIME_LIMITED, "error": None}, "POST /pattern"))
        self.assertEqual(mcp[("D-0001", 0)]["answer_digest_recorded"], "recorded-D-0001-0")
        # the mark
        sha = mark["checks"]["checks_sha256"]
        with open(os.path.join(SUITE, "checks.py"), "rb") as fh:
            self.assertEqual(sha, hashlib.sha256(fh.read()).hexdigest())
        self.assertEqual(mcp[("T-0001", 0)]["reevaluation"]["checks_sha256"], sha[:16])
        self.assertEqual(mcp[("T-0001", 0)]["reevaluation"]["outcome_recorded"], "fail")
        self.assertEqual(sc.read_json(os.path.join(self.run_dir, "reevaluated.json"))["checks"]["checks_sha256"], sha)
        changed = sorted((c["scenario"], c["arm"], c["repeat"], c["before"], c["after"]) for c in mark["changes"])
        self.assertIn(("T-0001", "mcp", 0, "fail", "inconclusive"), changed)
        self.assertIn(("P-0001", "direct-local", 0, "fail", "pass"), changed)
        self.assertIn(("D-0001", "mcp", 1, "fail", "inconclusive"), changed)
        self.assertNotIn("T-0001", [c[0] for c in changed if c[1] == "direct-local"])
        self.assertEqual(mark["arms"]["mcp"]["no_answer"], 1)
        self.assertEqual(mark["recorded"]["direct-local"]["scripts_commit"], "6f945db4")

    def test_merge_compare_report_keep_the_recorded_files(self):
        merge_run.merge(self.run_dir)
        sc.write_json_atomic(os.path.join(self.run_dir, "compare.json"),
                             cmp.compare(self.run_dir, cmp.load_by_design("")))
        report.main([self.run_dir, "--no-pdf"])
        before = digests(self.run_dir)
        reevaluate.main([self.run_dir, "--quiet"])
        merge_run.main([self.run_dir, "--reevaluated"])
        cmp.main([self.run_dir, "--by-design", "", "--reevaluated"])
        report.main([self.run_dir, "--no-pdf", "--reevaluated"])
        after = digests(self.run_dir)
        for path, digest in before.items():
            self.assertEqual(after.get(path), digest, path)
        new = sorted(set(after) - set(before))
        self.assertEqual(new, sorted(["arms/direct-local/executions.reevaluated.jsonl",
                                      "arms/direct-remote/executions.reevaluated.jsonl",
                                      "arms/mcp/executions.reevaluated.jsonl", "reevaluated.json",
                                      "executions.reevaluated.jsonl", "run.reevaluated.json",
                                      "compare.reevaluated.json", "report.reevaluated.md"]))
        run = sc.read_json(os.path.join(self.run_dir, "run.reevaluated.json"))
        self.assertEqual(run["reevaluation"]["checks"]["dir"], SUITE)
        self.assertGreater(run["reevaluation"]["changes"], 0)
        rows, _ = sc.read_jsonl(os.path.join(self.run_dir, "executions.reevaluated.jsonl"))
        by = dict(((r["scenario"], r["arm"], r["repeat"]), r["outcome"]) for r in rows)
        self.assertEqual(by[("T-0001", "mcp", 0)], "inconclusive")
        comp = sc.read_json(os.path.join(self.run_dir, "compare.reevaluated.json"))
        self.assertEqual(comp["scenarios"]["T-0001"]["arms"]["mcp"]["outcomes"], ["inconclusive", "fail"])
        self.assertEqual(comp["summary"]["relations"]["implementation"], "checks.evaluate_relation")
        self.assertIsNotNone(comp["reevaluation"])
        with open(os.path.join(self.run_dir, "report.reevaluated.md")) as fh:
            text = fh.read()
        self.assertIn("(re-evaluated)", text)
        self.assertIn("# Re-evaluation", text)
        self.assertIn("T-0001", text.split("# Re-evaluation", 1)[1].split("\n# ", 1)[0])
        with open(os.path.join(self.run_dir, "report.md")) as fh:
            self.assertNotIn("# Re-evaluation", fh.read())

    def test_checks_of_another_directory(self):
        # a checks.py of its own directory: it is passed what the runners pass, and the relations use it too
        d = os.path.join(self.tmp.name, "other")
        os.makedirs(d)
        with open(os.path.join(d, "checks.py"), "w") as fh:
            fh.write(textwrap.dedent('''
                import json
                def evaluate(scenario, answer, arm, previous=None, client=None):
                    return [{"check": "probe", "outcome": "pass", "property": True,
                             "detail": json.dumps({"previous": previous, "client": client}, sort_keys=True)}]
                def scenario_outcome(results):
                    return "pass"
                def answer_digest(answer, route=None):
                    return "D%s" % answer.get("http")
                def time_limited(answer):
                    return "TL" if "determinism" in json.dumps(answer.get("body")) else ""
            '''))
        mark = reevaluate.reevaluate(self.run_dir, d)
        self.assertEqual(mark["checks"]["dir"], d)
        self.assertIsNone(mark["checks"]["commit"])
        mcp, local = self.lines("mcp"), self.lines("direct-local")
        probe = json.loads(mcp[("D-0001", 1)]["checks"][0]["detail"])
        self.assertEqual(probe["previous"], [{"digest": "D200", "time_limited": "TL"}])
        self.assertIsNone(probe["client"])                                     # no backend time recorded
        self.assertEqual(json.loads(mcp[("T-0001", 1)]["checks"][0]["detail"])["client"], {"wall_ms": 9000.0})
        self.assertEqual(json.loads(local[("R-0001", 0)]["checks"][0]["detail"]),
                         {"client": {"wall_ms": 9257.0}, "previous": []})
        self.assertEqual(local[("P-0005", 0)]["outcome"], "pass")
        merge_run.merge(self.run_dir, reevaluated=True)
        comp = cmp.compare(self.run_dir, cmp.load_by_design(""), reevaluated=True)
        self.assertEqual(comp["summary"]["relations"]["implementation"], "compare.py built-in")

    def test_run_suite_reevaluate_is_offline(self):
        env = dict(os.environ)
        # any request or ssh would fail: nothing is reachable
        env.update({"SUITE_MEX": "local:" + os.path.join(self.tmp.name, "no-mex"), "SUITE_REMOTE_BASE": "http://127.0.0.1:9",
                    "SUITE_LOCAL_BASE": "http://127.0.0.1:9", "SUITE_MCP_CMD": "", "PYTHON": sys.executable})
        p = subprocess.run([os.path.join(SUITE, "run_suite.sh"), "--reevaluate", self.run_dir, "--no-pdf"], env=env,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, timeout=300,
                           stdin=subprocess.DEVNULL)
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertTrue(os.path.samefile(p.stdout.strip().splitlines()[-1], self.run_dir))
        for name in ("reevaluated.json", "executions.reevaluated.jsonl", "compare.reevaluated.json",
                     "report.reevaluated.md"):
            self.assertTrue(os.path.exists(os.path.join(self.run_dir, name)), name)
        self.assertFalse(os.path.exists(os.path.join(self.tmp.name, "no-mex")))
        self.assertFalse(os.path.exists(os.path.join(self.run_dir, "orchestrator.json")))
        p = subprocess.run([os.path.join(SUITE, "run_suite.sh"), "--reevaluate", self.tmp.name], env=env,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, timeout=60,
                           stdin=subprocess.DEVNULL)
        self.assertEqual(p.returncode, 2)
        self.assertIn("no arms/ directory", p.stderr)

    def test_answer_triggered_outcome(self):
        cases, _ = cmp.load_by_design(cmp.DEFAULT_BY_DESIGN)  # the real table: pattern.deadline_ends_task
        sc_ = scenario("N-0001", "POST /pattern", STATUS, request={"patterns": [], "time_budget_ms": 60000},
                       mcp={"tool": "pattern_search", "arguments": {}})
        relayed = {"http": 502, "body": None,
                   "error": {"code": "search_failed", "text": "pattern_deadline: pattern: the answer could not be "
                                                             "written within time_budget_ms (60000 ms, the "
                                                             "finalisation reserve of 250 ms included); retry"}}
        fail = [{"check": "status", "outcome": "fail", "detail": "http 502"}]
        at = lambda ms: {"outcome": "fail", "wall_ms": 99999.0, "mcp_phases": {"backend_ms": ms}}  # noqa: E731
        info = {}
        out, kept = reevaluate.outcome_of(None, "mcp", sc_, relayed, at(60055.0), fail, False, cases, None, info)
        self.assertEqual((out, kept), ("service_by_design",
                                       "the table's answer-triggered case pattern.deadline_ends_task"))
        self.assertIn(">= time_budget_ms 60000 - reserve 250 (stated) = 59750 ms", info["verdict_detail"])
        # R1 (wall_rule budget_minus_reserve): a known wall below budget - reserve holds the case back (the checks
        # judge); rule (c): no service-side wall, or no request budget: inconclusive, the case recorded
        self.assertEqual(reevaluate.outcome_of(None, "mcp", sc_, relayed, at(400.0), fail, False, cases)[0], "fail")
        info = {}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", sc_, relayed, {"outcome": "fail", "wall_ms": 61000.0},
                                               fail, False, cases, None, info)[0], "inconclusive")
        self.assertEqual((info["inconclusive_reason"], info["service_by_design_case"]),
                         ("no service-side timing: plausibility undecidable", "pattern.deadline_ends_task"))
        info = {}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", dict(sc_, request={}), relayed, at(61000.0), fail, False,
                                               cases, None, info)[0], "inconclusive")
        self.assertEqual(info["inconclusive_reason"], "no threshold known: plausibility undecidable")
        finding = fail + [{"check": "budget_respected", "outcome": "finding", "detail": ""}]
        self.assertEqual(reevaluate.outcome_of(None, "mcp", sc_, relayed, {"outcome": "fail"}, finding, False,
                                               cases)[0], "finding")
        # the wall of the deadline-ended task's one HTTP try (the diagnostic log's held-slot line)
        slot = "\u29d7 Slot held for search s 61.0s (att 1, ok, acq 0.1, rel 61.1, http 1 fail 1 60.10s backoff 0.0s, " \
               "released)"
        doc = {"calls": [{"tool": "GET /searches/s/diagnostic-log", "payload": {"regular_tasks": [
            {"error": "pattern_deadline: late", "diagnostic_log": ["x", slot]}]}}]}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", sc_, relayed, {"outcome": "fail"}, fail, False, cases,
                                               None, None, doc)[0], "service_by_design")
        two = dict(doc, calls=[{"tool": "GET /searches/s/diagnostic-log", "payload": {"regular_tasks": [
            {"error": "pattern_deadline: late", "diagnostic_log": [slot.replace("http 1 fail 1", "http 2 fail 1")]}]}}])
        self.assertEqual(reevaluate.outcome_of(None, "mcp", sc_, relayed, {"outcome": "fail"}, fail, False, cases,
                                               None, None, two)[0], "inconclusive")
        # another failed search, a direct arm, another route or no table: judged by the checks
        other = dict(relayed, error={"code": "search_failed", "text": "backend_down"})
        self.assertEqual(reevaluate.outcome_of(None, "mcp", sc_, other, at(61000.0), fail, False, cases)[0], "fail")
        self.assertEqual(reevaluate.outcome_of(None, "direct-local", sc_, relayed, at(61000.0), fail, False,
                                               cases)[0], "fail")
        res = scenario("R-0009", "POST /resolve", STATUS)
        self.assertEqual(reevaluate.outcome_of(None, "mcp", res, relayed, at(61000.0), fail, False, cases)[0], "fail")
        self.assertEqual(reevaluate.outcome_of(None, "mcp", sc_, relayed, at(61000.0), fail, False)[0], "fail")

    def test_host_state_and_answer_case_order(self):
        cases, _ = cmp.load_by_design(cmp.DEFAULT_BY_DESIGN)
        fail = [{"check": "status", "outcome": "fail", "detail": "http 504"}]
        finding = fail + [{"check": "budget_respected", "outcome": "finding", "detail": ""}]
        res = scenario("R-0009", "POST /resolve", STATUS, mcp={"tool": "traverse_resolve",
                                                               "arguments": {"time_budget_ms": 20000}})
        to = {"http": 504, "body": None, "error": {"code": "resolve_timeout",
                                                   "text": "the backend did not answer within 30.25 s"}}
        at = lambda ms: {"outcome": "fail", "wall_ms": 99999.0, "mcp_phases": {"backend_ms": ms}}  # noqa: E731
        # R2: least = max(20000 + 250 + 10000, 30.25 s) - 1000 ms = 29,250 ms on the service-side wall (R3)
        info = {}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", res, to, at(30375.2), fail, False, cases, None, info),
                         ("inconclusive", "the MCP runner's host-state rule (resolve_timeout)"))
        self.assertEqual(info["inconclusive_reason"], "resolve_timeout")
        info = {}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", res, to, at(800.0), fail, False, cases, None, info),
                         ("fail", None))
        self.assertIn("before max(time_budget_ms 20000 + reserve 250 + allowance 10000", info["verdict_detail"])
        # R3: the settled attempt's span, then the job's finished_at - started_at; never the client's wall
        line = {"outcome": "fail", "wall_ms": 30979.5}
        job = {"phase_sources": {"started_at": "2026-10-13T17:02:44.567542+00:00",
                                 "finished_at": "2026-10-13T17:03:14.942774+00:00"}}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", res, to, line, fail, False, cases, None, None, job)[0],
                         "inconclusive")
        span = {"phase_sources": dict(job["phase_sources"], attempts=[{"span_ms": 900.0}, {"span_ms": 1200.0}])}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", res, to, line, fail, False, cases, None, None, span)[0],
                         "fail")
        # rule (c): no service-side wall, or no threshold (no budget, no stated timeout): inconclusive, the reason
        info = {}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", res, to, line, fail, False, cases, None, info)[0],
                         "inconclusive")
        self.assertEqual((info["inconclusive_reason"], info["host_state"]),
                         ("no service-side timing: plausibility undecidable", "resolve_timeout"))
        bare = dict(to, error={"code": "resolve_timeout", "text": "the backend did not answer"})
        info = {}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", scenario("R-0011", "POST /resolve", STATUS), bare,
                                               at(30000.0), fail, False, cases, None, info)[0], "inconclusive")
        self.assertEqual(info["inconclusive_reason"], "no threshold known: plausibility undecidable")
        # an unbudgeted one judged against the timeout it states (N s - 1 s)
        self.assertEqual(reevaluate.outcome_of(None, "mcp", scenario("R-0011", "POST /resolve", STATUS), to,
                                               at(29300.0), fail, False, cases)[0], "inconclusive")
        # addendum 11 (a divergence from async's runner, which returns inconclusive first): an undecidable host
        # state never swallows a finding
        self.assertEqual(reevaluate.outcome_of(None, "mcp", res, to, line, finding, False, cases)[0], "finding")
        # G22: the relayed /resolve deadline: budget (the tool error's, else the request's) - reserve
        rdl = {"http": 503, "body": None, "error": {"code": "traversal_resolve_deadline", "text": "the backend could "
                                                    "not build and write its answer within bounds.time_budget_ms "
                                                    "(20000 ms, the finalisation reserve of 250 ms included)"}}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", res, rdl, at(22522.2), fail, False, cases)[0],
                         "inconclusive")
        self.assertEqual(reevaluate.outcome_of(None, "mcp", res, rdl, at(19700.0), fail, False, cases)[0], "fail")
        te = {"code": "traversal_resolve_deadline", "time_budget_ms": 5000, "finalize_reserve_ms": 250}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", res, rdl, at(5000.0), fail, False, cases, te)[0],
                         "inconclusive")
        info = {}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", res, rdl, {"outcome": "fail", "wall_ms": 22522.2}, fail,
                                               False, cases, None, info)[0], "inconclusive")   # N-realrep-0002 r0
        self.assertEqual(info["inconclusive_reason"], "no service-side timing: plausibility undecidable")
        named = scenario("R-0010", "POST /resolve", [{"check": "status", "args": {"http": 504,
                                                                                 "code": "resolve_timeout"}}])
        self.assertEqual(reevaluate.outcome_of(None, "mcp", named, to, at(30375.2), fail, False, cases)[0], "fail")
        # the R20 lockout answers at once: no wall rule
        lock = {"http": 503, "body": None, "error": {"code": "uninterruptible_overrun_exceeds_allowance", "text": "x"}}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", res, lock, {"outcome": "fail"}, fail, False, cases)[0],
                         "inconclusive")
        # B-06 and a runner error
        sc_ = scenario("N-0001", "POST /pattern", STATUS)
        info = {}
        out, kept = reevaluate.outcome_of(None, "mcp", sc_, {"http": 400, "error": {"code": "x"}},
                                          {"outcome": "service_by_design"}, fail, False, cases, None, info)
        self.assertEqual((out, kept), ("fail", None))   # B-06: judged by the checks
        self.assertIn("not reproduced by the current table", info["note"])
        # without a table to replay it the runner's mark is kept
        self.assertEqual(reevaluate.outcome_of(None, "mcp", sc_, {"http": 400, "error": {"code": "x"}},
                                               {"outcome": "service_by_design"}, fail, False, {})[0],
                         "service_by_design")
        self.assertEqual(reevaluate.outcome_of(None, "mcp", sc_, to, {"outcome": "error"}, fail, False, cases)[0],
                         "error")

    def test_named_cases_replayed(self):
        # the runner's matcher with the current table: G08's field alternation (Q-motif-0004 named the case, the
        # runner of the run could not read "predicate_strands | predicate_scope"), the field from the tool error
        cases, _ = cmp.load_by_design(cmp.DEFAULT_BY_DESIGN)
        q = scenario("Q-0004", "POST /pattern", [{"check": "status", "args": {"http": 400, "code": "invalid_request"}}],
                     mcp={"tool": "pattern_search", "arguments": {},
                          "by_design": ["pattern.predicate_options_without_predicate"]})
        ans = {"http": 400, "body": None, "error": {"code": "pattern_invalid", "field": "predicate_scope",
                                                    "text": "predicate_scope says ...: it acts only with a predicate"}}
        fail = [{"check": "status", "outcome": "fail", "detail": "code 'pattern_invalid'"}]
        self.assertEqual(reevaluate.outcome_of(None, "mcp", q, ans, {"outcome": "fail"}, fail, False, cases),
                         ("service_by_design", "the table's named case pattern.predicate_options_without_predicate"))
        # the runner's reading: a case field against an answer that names none does not match
        nofield = dict(ans, error={k: v for k, v in ans["error"].items() if k != "field"})
        self.assertEqual(reevaluate.outcome_of(None, "mcp", q, nofield, {"outcome": "fail"}, fail, False, cases)[0],
                         "fail")
        self.assertEqual(reevaluate.outcome_of(None, "mcp", q, nofield, {"outcome": "fail"}, fail, False, cases,
                                               {"code": "pattern_invalid", "field": "predicate_scope"})[0],
                         "service_by_design")
        doc = {"calls": [{"payload": {"error": {"code": "pattern_invalid", "field": "predicate_scope"}}}]}
        self.assertEqual(reevaluate.tool_error(doc, ans), {"code": "pattern_invalid", "field": "predicate_scope"})
        self.assertIsNone(reevaluate.tool_error(doc, {"error": {"code": "other"}}))
        # an answered case: no error and a 200; a direct arm is never matched
        self.assertTrue(cmp.service_answers({"outcome": "answered"}, {"http": 200, "error": None}))
        self.assertFalse(cmp.service_answers({"outcome": "answered"}, {"http": 400, "error": {"code": "x"}}))
        self.assertEqual(reevaluate.outcome_of(None, "direct-local", q, ans, {"outcome": "fail"}, fail, False,
                                               cases)[0], "fail")

    def test_answered_cases_judged_by_the_checks(self):
        # triage-2 B-01 as async's runner judges it (judge / explained_outcome at 48dcbc1d): an answered / clamped
        # case never decides before the checks; F failing, E = F & explains_checks, U = F - E
        cases = {"job": {"id": "job", "agreement": "by_design", "direct_expectation": {"http": 200},
                         "service_answer": {"outcome": "answered"}, "explains_checks": []},
                 "dup": {"id": "dup", "agreement": "by_design", "direct_expectation": {"http": 400},
                         "service_answer": {"outcome": "answered"}, "explains_checks": ["status"]},
                 "short": {"id": "short", "agreement": "by_design", "direct_expectation": {"http": 400},
                           "service_answer": {"outcome": "refused_at_submit", "code": "sequence_too_short"}}}
        ok = {"http": 200, "body": {}, "error": None}
        st = lambda o: {"check": "status", "outcome": o, "detail": ""}  # noqa: E731
        runs = {"check": "resolve_runs", "outcome": "pass", "detail": ""}
        by = lambda ids: scenario("R-0100", "POST /resolve", STATUS, mcp={"tool": "x", "by_design": ids})  # noqa: E731
        info = {}
        # nothing explained, nothing failing: the checks' outcome, the case a note
        self.assertEqual(reevaluate.outcome_of(None, "mcp", by("job"), ok, {"outcome": "service_by_design"},
                                               [st("pass"), runs], False, cases, None, info), ("pass", None))
        self.assertEqual((info["service_by_design_case"], info["explained_checks"]), ("job", []))
        # the content of every /resolve answer can fail again (P1 of VERIFY-REVIEW)
        bad = dict(runs, outcome="fail")
        self.assertEqual(reevaluate.outcome_of(None, "mcp", by("job"), ok, {"outcome": "pass"},
                                               [st("pass"), bad], False, cases)[0], "fail")
        # several named answered cases apply together: status explained by dup (R-edge-0009's naming)
        info = {}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", by(["job", "dup"]), ok, {"outcome": "fail"},
                                               [st("fail"), runs], False, cases, None, info)[0], "service_by_design")
        self.assertEqual((info["service_by_design_case"], info["explained_checks"]), ("dup", ["status"]))
        # U non-empty: fail; E with a finding left: finding (a mark never erases one)
        self.assertEqual(reevaluate.outcome_of(None, "mcp", by(["job", "dup"]), ok, {"outcome": "fail"},
                                               [st("fail"), bad], False, cases)[0], "fail")
        finding = {"check": "budget_respected", "outcome": "finding", "detail": ""}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", by(["job", "dup"]), ok, {"outcome": "fail"},
                                               [st("fail"), runs, finding], False, cases)[0], "finding")
        # a refusal case named after job_path still fails a 200 (P2): job_path matches, status fails unexplained
        self.assertEqual(reevaluate.outcome_of(None, "mcp", by(["job", "short"]), ok, {"outcome": "pass"},
                                               [st("fail"), runs], False, cases)[0], "fail")
        refused = {"http": 400, "body": None, "error": {"code": "sequence_too_short", "text": "short"}}
        info = {}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", by(["job", "short"]), refused, {"outcome": "fail"},
                                               [st("fail")], False, cases, None, info),
                         ("service_by_design", "the table's named case short"))
        self.assertEqual(info["explained_checks"], [])
        # a refused / relayed / answer-triggered match with a finding: finding, the case recorded
        info = {}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", by(["job", "short"]), refused, {"outcome": "fail"},
                                               [st("fail"), finding], False, cases, None, info)[0], "finding")
        self.assertEqual(info["service_by_design_case"], "short")

    def test_named_cases_before_answer_cases(self):
        # B-08 h: one case order, the runner's: a named case the answer matches wins over an answer-triggered one
        cases, _ = cmp.load_by_design(cmp.DEFAULT_BY_DESIGN)
        named = {"id": "x.relayed", "agreement": "by_design", "direct_expectation": {"http": 503},
                 "service_answer": {"outcome": "relayed", "code": "search_failed", "text_contains": "pattern_deadline"}}
        cases = dict(cases, **{"x.relayed": named})
        sc_ = scenario("N-0002", "POST /pattern", STATUS, mcp={"tool": "x", "by_design": ["x.relayed"]})
        relayed = {"http": 502, "body": None, "error": {"code": "search_failed", "text": "pattern_deadline: late"}}
        fail = [{"check": "status", "outcome": "fail", "detail": "http 502"}]
        at = {"outcome": "fail", "mcp_phases": {"backend_ms": 60100.0}}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", sc_, relayed, at, fail, False, cases),
                         ("service_by_design", "the table's named case x.relayed"))
        unnamed = dict(sc_, mcp={}, request={"time_budget_ms": 60000})
        self.assertEqual(reevaluate.outcome_of(None, "mcp", unnamed, relayed, at, fail, False, cases),
                         ("service_by_design", "the table's answer-triggered case pattern.deadline_ends_task"))
        # a relayed deadline whose budget_respected is a finding keeps it (Q1 of VERIFY-REVIEW)
        self.assertEqual(reevaluate.outcome_of(None, "mcp", unnamed, relayed, at,
                                               fail + [{"check": "budget_respected", "outcome": "finding"}], False,
                                               cases)[0], "finding")

    def test_rows_record_the_case_and_no_table(self):
        # the re-evaluated lines record the matched case and the explained checks (the replay's); '' = no table
        reevaluate.reevaluate(self.run_dir, by_design="")
        mark = sc.read_json(os.path.join(self.run_dir, "reevaluated.json"))
        self.assertIsNone(mark["by_design"])
        self.assertEqual(self.lines("mcp")[("P-0003", 0)]["outcome"], "service_by_design")   # kept: no table
        bad = os.path.join(self.tmp.name, "bad.json")
        sc.write_json_atomic(bad, {"cases": [{"id": "c", "agreement": "by_design", "service_answer": {"outcome": "relayed"},
                                              "answer_match": {"direct": {}, "mcp": {"http": [400, 502]}}}]})
        self.assertEqual(reevaluate.main([self.run_dir, "--by-design", bad, "--quiet"]), 4)
        with self.assertRaises(cmp.TableError) as ctx:
            cmp.load_by_design(bad)
        self.assertIn("answer_match without routes", str(ctx.exception))
        self.assertIn("must be one integer", str(ctx.exception))

    def test_merge_without_a_reevaluation(self):
        with self.assertRaises(SystemExit):
            merge_run.merge(self.run_dir, reevaluated=True)


if __name__ == "__main__":
    unittest.main()
