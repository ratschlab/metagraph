"""Offline tests of reevaluate.py and the --reevaluated set of merge_run.py, compare.py, report.py and run_suite.sh
on a synthetic recorded run: what is judged again and what is kept, what the checks are passed (previous, client),
the mark of the set, and that every recorded file stays as it was."""

import contextlib
import hashlib
import io
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
        # R2 as addendum 1 item 4 caps it for a refusal that states no phase: least = min(10,000 ms (the write
        # phase), max(20000 + 250 + 10000, 30.25 s)) - 1000 ms = 9,000 ms on the service-side wall (R3)
        info = {}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", res, to, at(30375.2), fail, False, cases, None, info),
                         ("inconclusive", "the MCP runner's host-state rule (resolve_timeout)"))
        self.assertEqual(info["inconclusive_reason"], "resolve_timeout")
        info = {}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", res, to, at(800.0), fail, False, cases, None, info),
                         ("fail", None))
        self.assertEqual(info["verdict_detail"], "resolve_timeout after 800 ms (mcp_phases.backend_ms), before min(the "
                                                 "write cap 10000 ms, max(time_budget_ms 20000 + reserve 250 + "
                                                 "allowance 10000 = 30250 ms, the stated timeout 30.25 s)) - 1000 ms = "
                                                 "9000 ms: not a host state")
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
        # B5b-03: an undecidable host state is inconclusive whatever the checks say, as in async's runner (lane S5b's
        # addendum-11 divergence, a finding there, is dropped)
        info = {}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", res, to, line, finding, False, cases, None, info)[0],
                         "inconclusive")
        self.assertEqual(info["inconclusive_reason"], "no service-side timing: plausibility undecidable")
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

    @staticmethod
    def calls(*spec):
        """A raw record's calls: (tool, t0, t1[, payload])."""
        return {"calls": [{"tool": x[0], "t0": x[1], "t1": x[2], "timed": True,
                           "payload": x[3] if len(x) > 3 else {}} for x in spec]}

    def test_client_wall_rule2(self):
        # rule 2 (async's 2a8acbaf; B5b-01, S5c-01, S5c-03): the service-side wall unknown, the client's wall from
        # the FIRST submit call's t0 to the answering call's t1 bounds it from above for ANONYMOUS MCP calls only;
        # + 1000 ms below the threshold the answer is provably premature: R1 the case does not apply, R3
        # implausible (the checks judge); the runner's verdict_detail texts
        cases, _ = cmp.load_by_design(cmp.DEFAULT_BY_DESIGN)
        fail = [{"check": "status", "outcome": "fail", "detail": "http 502"}]
        pat = scenario("M-planted-0016", "POST /pattern", STATUS, request={"patterns": [], "time_budget_ms": 60000},
                       mcp={"tool": "pattern_search", "arguments": {"time_budget_ms": 60000}})
        relayed = {"http": 502, "body": None,
                   "error": {"code": "search_failed", "text": "pattern_deadline: pattern: the answer could not be "
                                                             "written within time_budget_ms (60000 ms, the "
                                                             "finalisation reserve of 250 ms included)"}}
        line = {"outcome": "fail", "wall_ms": 5.0, "mcp_phases": {"submit_ms": 2.0}}
        anon = dict(line, mcp_auth="anonymous")
        done = {"state": "failed", "terminal": True}
        fast = self.calls(("pattern_search", 1000.0, 1000.002), ("get_status", 1000.002, 1000.005, done),
                          ("get_results", 1000.005, 1000.4))
        premature = ("pattern.deadline_ends_task (wall_rule budget_minus_reserve): client wall 5 ms + 1 s slack < "
                     "threshold 59750 ms: provably premature")
        # S5c-01: a record without mcp_auth (every run before 2a8acbaf) is unknown: rule 2 decides nothing
        info = {}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", pat, relayed, line, fail, False, cases, None, info,
                                               fast)[0], "inconclusive")
        self.assertEqual((info["inconclusive_reason"], info["verdict_detail"], info["service_by_design_case"]),
                         (sc.UNDECIDABLE_NO_WALL, sc.UNDECIDABLE_NO_WALL, "pattern.deadline_ends_task"))
        # anonymous (the line's, the raw record's, or the re-evaluation's stated assumption): premature
        for ln, doc in ((anon, fast), (line, dict(fast, mcp_auth="anonymous")),
                        (dict(line, reevaluation={"mcp_auth_assumed": "anonymous"}), fast)):
            info = {}
            self.assertEqual(reevaluate.outcome_of(None, "mcp", pat, relayed, ln, fail, False, cases, None, info,
                                                   doc)[0], "fail")
            self.assertEqual(info["verdict_detail"], premature)
        # a recorded mcp_auth wins over the assumption; an unknown value is no auth
        for ln in (dict(anon, mcp_auth="oauth"), dict(line, mcp_auth="api_key",
                                                     reevaluation={"mcp_auth_assumed": "anonymous"})):
            self.assertEqual(reevaluate.outcome_of(None, "mcp", pat, relayed, ln, fail, False, cases, None, None,
                                                   fast)[0], "inconclusive")
        # an API key: skipped, the runner's note after the undecidable text, the class without it
        info = {}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", pat, relayed, dict(line, mcp_auth="api_key"), fail,
                                               False, cases, None, info, fast)[0], "inconclusive")
        self.assertEqual(info["verdict_detail"], sc.UNDECIDABLE_NO_WALL + "; " + sc.RULE2_SKIPPED_API_KEY)
        self.assertEqual(info["inconclusive_reason"], sc.UNDECIDABLE_NO_WALL)
        # a busy retry: the first submit (idempotency_conflict) 100 s before the accepted one: no bound below it
        busy = self.calls(("pattern_search", 900.0, 900.1, {"error": {"code": "idempotency_conflict"}}),
                          ("pattern_search", 1000.0, 1000.002), ("get_status", 1000.002, 1000.005, done))
        info = {}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", pat, relayed, anon, fail, False, cases, None, info,
                                               busy)[0], "inconclusive")
        self.assertEqual(info["verdict_detail"], sc.UNDECIDABLE_NO_WALL)
        self.assertEqual(sc.client_wall(pat, busy)["ms"], 100005.0)
        # S5c-03 (a): the first submit's t0 unknown: no client wall (not the first submit that has one)
        no_t0 = self.calls(("pattern_search", 1000.0, 1000.002), ("pattern_search", 1100.0, 1100.002),
                           ("get_status", 1100.002, 1100.005, done))
        del no_t0["calls"][0]["t0"]
        self.assertIsNone(sc.client_wall(pat, no_t0))
        self.assertEqual(reevaluate.outcome_of(None, "mcp", pat, relayed, anon, fail, False, cases, None, None,
                                               no_t0)[0], "inconclusive")
        # S5c-03 (b): the FIRST terminal get_status answered (5 ms), not a later one (100 s)
        two = self.calls(("pattern_search", 1000.0, 1000.002), ("get_status", 1000.002, 1000.005, done),
                         ("get_status", 1000.005, 1100.0, done))
        self.assertEqual(sc.client_wall(pat, two)["answer_call"], 1)
        self.assertEqual(reevaluate.outcome_of(None, "mcp", pat, relayed, anon, fail, False, cases, None, None,
                                               two)[0], "fail")
        # within the slack of the threshold: undecidable as before (no raw calls: undecidable too)
        near = self.calls(("pattern_search", 1000.0, 1000.002), ("get_status", 1000.002, 1059.0, done))
        self.assertEqual(reevaluate.outcome_of(None, "mcp", pat, relayed, anon, fail, False, cases, None, None,
                                               near)[0], "inconclusive")
        # R3: a relayed /resolve deadline 30 ms and a resolve_timeout 40 ms against a 20,000 ms budget
        res = scenario("R-pos-0001", "POST /resolve", STATUS, mcp={"tool": "traverse_resolve",
                                                                   "arguments": {"time_budget_ms": 20000}})
        rdl = {"http": 503, "body": None, "error": {"code": "traversal_resolve_deadline", "text": "within "
                                                    "bounds.time_budget_ms (20000 ms, the finalisation reserve of "
                                                    "250 ms included)"}}
        to = {"http": 504, "body": None, "error": {"code": "resolve_timeout",
                                                   "text": "the backend did not answer within 30.25 s"}}
        for ans, ms, least in ((rdl, 0.030, 19750), (to, 0.040, 9000)):
            info = {}
            doc = self.calls(("traverse_resolve", 2000.0, 2000.0 + ms, {"error": ans["error"]}))
            self.assertEqual(reevaluate.outcome_of(None, "mcp", res, ans, {"outcome": "fail"}, fail, False, cases,
                                                   None, info, doc)[0], "inconclusive", ans["error"]["code"])
            self.assertEqual(info["inconclusive_reason"], sc.UNDECIDABLE_NO_WALL)   # S5c-01: no mcp_auth
            info = {}
            self.assertEqual(reevaluate.outcome_of(None, "mcp", res, ans, {"outcome": "fail", "mcp_auth": "anonymous"},
                                                   fail, False, cases, None, info, doc)[0], "fail")
            self.assertEqual(info["verdict_detail"], "client wall %.0f ms + 1 s slack < threshold %d ms: provably "
                                                     "premature" % (ms * 1000, least))
            # N-realrep-0002 r0's shape (22,522 ms >= 19,750 ms): undecidable as recorded
            doc = self.calls(("traverse_resolve", 2000.0, 2022.522, {"error": ans["error"]}))
            self.assertEqual(reevaluate.outcome_of(None, "mcp", res, ans, {"outcome": "fail", "mcp_auth": "anonymous"},
                                                   fail, False, cases, None, None, doc)[0], "inconclusive")
        # the job path: the answering call is the last TIMED traverse_results
        job = self.calls(("traverse_resolve", 2000.0, 2000.01, {"job_id": "j"}), ("traverse_status", 2000.01, 2030.0),
                         ("traverse_results", 2030.0, 2030.2), ("traverse_results", 2030.2, 2090.0))
        job["calls"][-1]["timed"] = False
        wall = sc.client_wall(res, job)
        self.assertEqual((wall["ms"], wall["answer_call"], wall["answer_tool"]), (30200.0, 2, "traverse_results"))
        self.assertIsNone(sc.client_wall(res, None))
        self.assertIsNone(sc.client_wall(scenario("T-1", "POST /traverse", STATUS,
                                                  mcp={"tool": "traverse_fetch", "arguments": {}}), job))

    def test_recorded_client_wall_and_terms(self):
        # S5c-03: the runner's recorded client wall (raw client_wall, else the line's client_wall_ms) and timeout
        # statement (raw / line timeout_terms) are read first; a recomputation that differs is a warning
        cases, _ = cmp.load_by_design(cmp.DEFAULT_BY_DESIGN)
        fail = [{"check": "status", "outcome": "fail", "detail": "http 504"}]
        res = scenario("R-pos-0001", "POST /resolve", STATUS, mcp={"tool": "traverse_resolve",
                                                                   "arguments": {"time_budget_ms": 20000}})
        to = {"http": 504, "body": None, "error": {"code": "resolve_timeout",
                                                   "text": "the backend did not answer within 30.25 s"}}
        doc = dict(self.calls(("traverse_resolve", 2000.0, 2000.04, {"error": to["error"]})), mcp_auth="anonymous")
        line = {"outcome": "fail"}
        self.assertEqual(sc.client_wall_ms(res, line, doc), 40.0)
        warns = []
        self.assertEqual(sc.client_wall_ms(res, line, dict(doc, client_wall={"ms": 40.0}), warns), 40.0)
        self.assertEqual(warns, [])
        info = {}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", res, to, line, fail, False, cases, None, info,
                                               dict(doc, client_wall={"ms": 9500.0}))[0], "inconclusive")
        self.assertIn("the raw record's client_wall 9500.0 differs from the calls' client wall 40.0",
                      info["warnings"][0])
        self.assertEqual(sc.client_wall_ms(res, dict(line, client_wall_ms=None), doc, warns), None)
        self.assertIn("the line's client_wall_ms None differs", warns[0])
        # timeout_terms: recorded (a stated write 10 s) against a record that states none; null = legacy
        stated = {"timeout_phase": "write", "bound_s": 10.0, "elapsed_ms": 10012.0, "source": "tool error"}
        at = {"outcome": "fail", "mcp_phases": {"backend_ms": 9500.0}}
        info = {}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", res, to, at, fail, False, cases, None, info,
                                               {"timeout_terms": stated})[0], "inconclusive")
        self.assertIn("the stated bound 10 s (write) - 1000 ms = 9000 ms", info["verdict_detail"])
        self.assertIn("the raw record's timeout_terms", info["warnings"][0])
        info = {}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", res, to, dict(at, timeout_terms=None), fail, False, cases,
                                               {"code": "resolve_timeout", "timeout_phase": "write", "bound_s": 5},
                                               info)[0], "inconclusive")
        self.assertIn("min(the write cap 10000 ms,", info["verdict_detail"])
        self.assertIn("the line's timeout_terms null differs", info["warnings"][0])
        # a malformed recorded statement is replaced by the record's own
        info = {}
        reevaluate.outcome_of(None, "mcp", res, to, dict(at, timeout_terms={"timeout_phase": "connect"}), fail, False,
                              cases, None, info)
        self.assertIn("states no valid phase and bound", info["warnings"][0])

    def test_mcp_auth_option(self):
        # S5c-01: a run before async's 2a8acbaf records no mcp_auth: rule 2 decides nothing unless the re-evaluation
        # is told (--mcp-auth anonymous), which the mark and each line it applied to record; compare on the
        # re-evaluated set reads the same assumption
        run = os.path.join(self.tmp.name, "auth")
        pat = scenario("M-planted-0016", "POST /pattern", STATUS, request={"patterns": [], "time_budget_ms": 60000},
                       mcp={"tool": "pattern_search", "arguments": {"time_budget_ms": 60000}})
        manifest = {"manifest_version": 1, "scenarios": [pat], "relations": []}
        err = {"code": "search_failed", "text": "pattern_deadline: pattern: the answer could not be written within "
                                                "time_budget_ms (60000 ms, the finalisation reserve of 250 ms "
                                                "included)"}
        os.makedirs(run)
        sc.write_json_atomic(os.path.join(run, "manifest.json"), manifest)
        for arm, http, error in (("mcp", 502, err), ("direct-local", 200, None)):
            d = os.path.join(run, "arms", arm)
            os.makedirs(os.path.join(d, "raw"))
            raw = sc.raw_name(pat["id"], arm, 0)
            sc.append_jsonl(os.path.join(d, "executions.jsonl"), {
                "run_id": "R", "scenario": pat["id"], "arm": arm, "repeat": 0, "cache": "cold", "http": http,
                "outcome": "fail" if arm == "mcp" else "pass", "checks": [], "mcp_phases": {} if arm == "mcp" else None,
                "raw": raw, "error": error, "route": "POST /pattern", "started_utc": sc.iso_utc(T0), "wall_ms": 5.0})
            doc = {"answer": {"http": http, "body": None if error else {"patterns": []}, "error": error}}
            if arm == "mcp":
                doc["calls"] = self.calls(("pattern_search", 1000.0, 1000.002),
                                          ("get_status", 1000.002, 1000.005, {"terminal": True}))["calls"]
            sc.write_gz_json(os.path.join(d, raw), doc)
            sc.write_json_atomic(os.path.join(d, "run.json"), {"run_id": "R", "arm": arm, "checks_module": "checks"})
            sc.write_json_atomic(os.path.join(d, "manifest.json"), manifest)
        mark = reevaluate.reevaluate(run)
        row = sc.read_jsonl(os.path.join(run, "arms", "mcp", "executions.reevaluated.jsonl"))[0][0]
        self.assertEqual((row["outcome"], row["inconclusive_reason"]), ("inconclusive", sc.UNDECIDABLE_NO_WALL))
        self.assertEqual(mark["mcp_auth"], {"assumed": None, "recorded": {"unknown": 1}, "assumed_lines": 0})
        self.assertNotIn("mcp_auth_assumed", row["reevaluation"])
        self.assertEqual(reevaluate.main([run, "--quiet", "--mcp-auth", "anonymous"]), 0)
        mark = sc.read_json(os.path.join(run, "reevaluated.json"))
        row = sc.read_jsonl(os.path.join(run, "arms", "mcp", "executions.reevaluated.jsonl"))[0][0]
        self.assertEqual(mark["mcp_auth"], {"assumed": "anonymous", "recorded": {"unknown": 1}, "assumed_lines": 1})
        self.assertEqual((row["outcome"], row["reevaluation"]["mcp_auth_assumed"]), ("fail", "anonymous"))
        self.assertIn("client wall 5 ms + 1 s slack < threshold 59750 ms", row["verdict_detail"])
        direct = sc.read_jsonl(os.path.join(run, "arms", "direct-local", "executions.reevaluated.jsonl"))[0][0]
        cases, _ = cmp.load_by_design(cmp.DEFAULT_BY_DESIGN)
        row["raw"] = os.path.join("arms", "mcp", row["raw"])
        res = cmp.mcp_agreement({0: row}, {"direct-local": {0: direct}}, pat, cases, run)
        self.assertEqual(res["class"], "disagree")      # the case held back: the checks compared
        bare = dict(row, reevaluation={})
        self.assertEqual(cmp.mcp_agreement({0: dict(bare, outcome="inconclusive")}, {"direct-local": {0: direct}},
                                           pat, cases, run)["class"], "not_comparable")
        # the option takes anonymous or api_key only
        with self.assertRaises(SystemExit), contextlib.redirect_stderr(io.StringIO()):
            reevaluate.main([run, "--quiet", "--mcp-auth", "maybe"])

    def test_discovery_without_budget(self):
        # A2-01 MCP side with addendum 1 item 1: a relayed /resolve deadline on a discovery that names no
        # time_budget_ms is implausible (the service sends no deadline for a discovery; SPEC-core §4.5 opt-in):
        # the checks judge (status fails); an explicit-labels resolve without a caller budget is "budget unknown"
        cases, _ = cmp.load_by_design(cmp.DEFAULT_BY_DESIGN)
        fail = [{"check": "status", "outcome": "fail", "detail": "http 503"}]
        disc = scenario("R-discover-0001", "POST /resolve", STATUS,
                        mcp={"tool": "traverse_resolve", "arguments": {"database": "db", "discover": True,
                                                                       "kind": "taxid", "sequence": "ACGT"}})
        rdl = {"http": 503, "body": None, "error": {"code": "traversal_resolve_deadline", "text": "late"}}
        te = {"code": "traversal_resolve_deadline", "time_budget_ms": None}
        for line in ({"outcome": "fail", "mcp_phases": {"backend_ms": 4.0}}, {"outcome": "fail"}):
            info = {}
            self.assertEqual(reevaluate.outcome_of(None, "mcp", disc, rdl, line, fail, False, cases, te, info)[0],
                             "fail")
            self.assertEqual(info["verdict_detail"], sc.DISCOVERY_NO_DEADLINE)
        self.assertTrue(sc.discovery_without_budget(disc))
        # the runner's unbudgeted_discovery reads the arguments whatever the tool (only traverse_resolve has them)
        self.assertTrue(sc.discovery_without_budget(dict(disc, mcp={"tool": "x", "arguments": {"discover": True}})))
        self.assertFalse(sc.discovery_without_budget(dict(disc, mcp={"tool": "x", "arguments": ["discover"]})))
        timed = dict(disc, mcp={"tool": "traverse_resolve", "arguments": dict(disc["mcp"]["arguments"],
                                                                                time_budget_ms=20000)})
        self.assertFalse(sc.discovery_without_budget(timed))
        self.assertEqual(reevaluate.outcome_of(None, "mcp", timed, rdl, {"outcome": "fail", "mcp_phases": {
            "backend_ms": 20500.0}}, fail, False, cases)[0], "inconclusive")
        # explicit labels, no caller budget: the stated budget's threshold, else undecidable
        labels = scenario("R-0001", "POST /resolve", STATUS, mcp={"tool": "traverse_resolve",
                                                                  "arguments": {"labels": ["562"], "sequence": "AC"}})
        self.assertFalse(sc.discovery_without_budget(labels))
        at = {"outcome": "fail", "mcp_phases": {"backend_ms": 20000.0}}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", labels, rdl, at, fail, False, cases,
                                               {"code": "traversal_resolve_deadline", "time_budget_ms": 20000})[0],
                         "inconclusive")
        info = {}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", labels, rdl, at, fail, False, cases, None, info)[0],
                         "inconclusive")
        self.assertEqual(info["inconclusive_reason"], "no threshold known: plausibility undecidable")

    def test_resolve_timeout_stated_phase(self):
        # rule 4 (async's f94b079d): a resolve_timeout that states its phase and bound is a host state at bound_s x
        # 1000 - 1000 ms or later, any phase (the phase named in the detail); a legacy one at min(10 s, max(A, B)) -
        # 1 s; the runner's verdict_detail texts
        cases, _ = cmp.load_by_design(cmp.DEFAULT_BY_DESIGN)
        fail = [{"check": "status", "outcome": "fail", "detail": "http 504"}]
        res = scenario("R-0009", "POST /resolve", STATUS, mcp={"tool": "traverse_resolve",
                                                               "arguments": {"time_budget_ms": 20000}})
        to = {"http": 504, "body": None, "error": {"code": "resolve_timeout",
                                                   "text": "the backend did not answer within 10 s"}}
        at = lambda ms: {"outcome": "fail", "mcp_phases": {"backend_ms": ms}}  # noqa: E731
        write = {"code": "resolve_timeout", "timeout_phase": "write", "bound_s": 10, "elapsed_ms": 10012}
        info = {}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", res, to, at(9500.0), fail, False, cases, write, info)[0],
                         "inconclusive")
        self.assertEqual(info["verdict_detail"], "resolve_timeout after 9500 ms (mcp_phases.backend_ms) >= the stated "
                                                 "bound 10 s (write) - 1000 ms = 9000 ms")
        self.assertEqual(reevaluate.outcome_of(None, "mcp", res, to, at(8000.0), fail, False, cases, write)[0],
                         "fail")
        read = dict(write, timeout_phase="read", bound_s=30.25)
        self.assertEqual(reevaluate.outcome_of(None, "mcp", res, to, at(30375.0), fail, False, cases, read)[0],
                         "inconclusive")
        self.assertEqual(reevaluate.outcome_of(None, "mcp", res, to, at(20000.0), fail, False, cases, read)[0],
                         "fail")
        # legacy: no stated phase; the recorded ~30.4 s lines stay host states, a 9.5 s write timeout is one now
        self.assertEqual(sc.resolve_timeout_least(res, to["error"], None), (9000.0, (
            "min(the write cap 10000 ms, max(time_budget_ms 20000 + reserve 250 + allowance 10000 = 30250 ms, the "
            "stated timeout 10 s)) - 1000 ms")))
        self.assertEqual(reevaluate.outcome_of(None, "mcp", res, to, at(30375.2), fail, False, cases)[0],
                         "inconclusive")
        self.assertEqual(reevaluate.outcome_of(None, "mcp", res, to, at(9500.0), fail, False, cases)[0],
                         "inconclusive")
        self.assertEqual(reevaluate.outcome_of(None, "mcp", res, to, at(8000.0), fail, False, cases)[0], "fail")
        bare = scenario("R-0011", "POST /resolve", STATUS)
        short = dict(to, error={"code": "resolve_timeout", "text": "the backend did not answer within 5 s"})
        self.assertEqual(sc.resolve_timeout_least(bare, short["error"], None)[0], 4000.0)
        # S5c-02: the statement is validated as the runner's timeout_statement: a phase in read | write | deadline
        # AND a positive numeric bound_s, both in ONE object (the tool error, else the answer error, else the
        # job's stored error); anything else is legacy (backend_ms 8 against the legacy 9,000 ms: a fault)
        for bad, answer_error in ((dict(write, bound_s=0), {}), (dict(write, bound_s=-5), {}),
                                  (dict(write, timeout_phase="connect", bound_s=0.5), {}),
                                  ({"code": "resolve_timeout", "timeout_phase": "write"}, {"bound_s": 0.5}),
                                  (dict(write, bound_s="10"), {})):
            ans = dict(to, error=dict(to["error"], **answer_error))
            self.assertIsNone(sc.timeout_statement(None, ans["error"], bad), bad)
            self.assertEqual(reevaluate.outcome_of(None, "mcp", res, ans, at(8.0), fail, False, cases, bad)[0],
                             "fail", bad)
        # the answer error states it when the tool error does not; the job's stored error (any code) third
        self.assertEqual(sc.timeout_statement(None, dict(write, code="resolve_timeout"), {"code": "resolve_timeout"})[
            "source"], "answer error")
        job = {"calls": [{"tool": "traverse_resolve", "payload": {"job_id": "j"}},
                         {"tool": "traverse_results", "payload": {"resolve_error": {
                             "status": 504, "detail": {"code": "other", "timeout_phase": "write", "bound_s": 5}}}}]}
        self.assertEqual(sc.timeout_statement(job, to["error"], None)["source"], "job error")
        info = {}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", res, to, at(6000.0), fail, False, cases, None, info,
                                               job)[0], "inconclusive")
        self.assertIn("the stated bound 5 s (write) - 1000 ms = 4000 ms", info["verdict_detail"])

    def test_resolve_timeout_read_floor(self):
        # rule 4's read floor (async's f52846fb): a stated READ bound with a known budget (the MCP arguments', else the
        # tool error's time_budget_ms) below budget + reserve (the tool error's finalize_reserve_ms, else 250) +
        # 10,000 - 1000 ms is no host state whatever the wall (an unknown one included): the checks judge,
        # verdict_detail "service read bound {b} s below its own floor {f} s" (%g); write / deadline phases keep
        # bound_s - 1 s; no budget known and a legacy record are unchanged (the runner's test, mirrored)
        cases, _ = cmp.load_by_design(cmp.DEFAULT_BY_DESIGN)
        fail = [{"check": "status", "outcome": "fail", "detail": "http 504"}]

        def res(budget):
            args = {"labels": ["562"]} if budget is None else {"time_budget_ms": budget}
            return scenario("R-floor", "POST /resolve", STATUS, mcp={"tool": "traverse_resolve", "arguments": args})
        to = {"http": 504, "body": None, "error": {"code": "resolve_timeout",
                                                   "text": "the backend did not answer within 20.25 s"}}

        def stated(phase, bound, **extra):
            return dict({"code": "resolve_timeout", "timeout_phase": phase, "bound_s": bound,
                         "elapsed_ms": bound * 1000 + 3.0}, **extra)

        def judge(budget, te, wall_ms, doc=None, line=None):
            info = {}
            line = dict(line or {"outcome": "fail"}, mcp_phases={"backend_ms": wall_ms})
            out = reevaluate.outcome_of(None, "mcp", res(budget), to, line, fail, False, cases, te, info, doc)[0]
            return out, info.get("verdict_detail"), info.get("inconclusive_reason")
        # below the floor (20 + 0.25 + 10 - 1 = 29.25 s): the checks judge at a known, a long and an unknown wall
        below = "service read bound 20.25 s below its own floor 29.25 s"
        for wall in (20300.0, 45000.0):
            self.assertEqual(judge(20000, stated("read", 20.25), wall), ("fail", below, None))
        anon = dict(self.calls(("traverse_resolve", 2000.0, 2020.3, {"error": stated("read", 20.25)})),
                    mcp_auth="anonymous")
        self.assertEqual(judge(20000, stated("read", 20.25), None, anon), ("fail", below, None))
        self.assertEqual(judge(20000, stated("read", 20.25), None), ("fail", below, None))
        self.assertEqual(sc.host_state_verdict(res(20000), {"mcp_phases": {}}, None, to["error"],
                                               stated("read", 20.25)), ("resolve_timeout", "not_plausible", below))
        self.assertEqual(sc.read_bound_below_floor(res(20000), to["error"], stated("read", 29.249)),
                         "service read bound 29.249 s below its own floor 29.25 s")
        # at or above the floor the stated bound - 1 s decides
        for bound in (29.25, 30.25):
            self.assertIsNone(sc.read_bound_below_floor(res(20000), to["error"], stated("read", bound)))
            out, detail, why = judge(20000, stated("read", bound), bound * 1000)
            self.assertEqual((out, why), ("inconclusive", "resolve_timeout"), bound)
            self.assertIn("the stated bound %g s (read)" % bound, detail)
        # the reserve the tool error states raises the floor (20 + 2 + 10 - 1 = 31 s)
        self.assertEqual(judge(20000, stated("read", 30.25, finalize_reserve_ms=2000.0), 30250.0),
                         ("fail", "service read bound 30.25 s below its own floor 31 s", None))
        # the tool error's budget when the arguments name none; no budget known: no floor
        self.assertEqual(judge(None, stated("read", 20.25, time_budget_ms=20000), 20300.0), ("fail", below, None))
        self.assertEqual(judge(None, stated("read", 20.0), 19000.0)[::2], ("inconclusive", "resolve_timeout"))
        # write and deadline phases keep bound_s - 1 s on a budgeted call
        for phase, bound in (("write", 10.0), ("deadline", 15.0)):
            self.assertIsNone(sc.read_bound_below_floor(res(20000), to["error"], stated(phase, bound)))
            self.assertEqual(judge(20000, stated(phase, bound), bound * 1000 - 1000)[::2],
                             ("inconclusive", "resolve_timeout"), phase)
        # a legacy record (no stated phase) keeps min(10 s, ...) - 1 s
        self.assertEqual(judge(20000, {"code": "resolve_timeout"}, 9000.0)[::2], ("inconclusive", "resolve_timeout"))
        # the runner's recorded statement is what the floor reads (the raw record's timeout_terms)
        rec = {"timeout_phase": "read", "bound_s": 20.25, "elapsed_ms": 20253.0, "source": "tool error"}
        self.assertEqual(judge(20000, stated("read", 20.25), 20300.0, {"timeout_terms": rec}), ("fail", below, None))
        self.assertEqual(judge(20000, stated("read", 20.25), 20300.0, None, {"outcome": "fail", "timeout_terms": rec}),
                         ("fail", below, None))
        # merge_run (and compare) read the same rule for a line recorded without inconclusive_reason: the answer
        # error's statement below the floor is no host state, one at the floor is
        for bound, host in ((20.25, False), (29.25, True)):
            row = {"arm": "mcp", "outcome": "inconclusive", "mcp_phases": {"backend_ms": 45000.0},
                   "error": dict(to["error"], timeout_phase="read", bound_s=bound)}
            self.assertIs(merge_run.host_state_inconclusive(row, res(20000)), host, bound)

    def test_by_design_finding_without_a_table(self):
        # B5b-10: with no table a recorded by-design finding (a relayed deadline's budget_respected) is kept, as
        # the runner's service_by_design mark is; with the table it is replayed (finding: a mark never erases one)
        sc_ = scenario("N-0001", "POST /pattern", STATUS, request={"patterns": [], "time_budget_ms": 60000},
                       mcp={"tool": "pattern_search", "arguments": {}})
        relayed = {"http": 502, "body": None, "error": {"code": "search_failed", "text": "pattern_deadline: late"}}
        checks = [{"check": "status", "outcome": "fail", "detail": "http 502"},
                  {"check": "budget_respected", "outcome": "finding", "detail": "late"}]
        line = {"outcome": "finding", "service_by_design_case": "pattern.deadline_ends_task",
                "mcp_phases": {"backend_ms": 125000.0}}
        out, kept = reevaluate.outcome_of(None, "mcp", sc_, relayed, line, checks, False, {})
        self.assertEqual(out, "finding")
        self.assertIn("pattern.deadline_ends_task", kept)
        cases, _ = cmp.load_by_design(cmp.DEFAULT_BY_DESIGN)
        self.assertEqual(reevaluate.outcome_of(None, "mcp", sc_, relayed, line, checks, False, cases),
                         ("finding", "the table's answer-triggered case pattern.deadline_ends_task; a check is a "
                                     "finding (a by-design mark never erases one)"))
        # a finding without a case is the checks' (re-evaluated), with or without the table
        self.assertEqual(reevaluate.outcome_of(None, "mcp", sc_, relayed, {"outcome": "finding"}, checks[:1], False,
                                               {})[0], "fail")

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
        self.assertEqual((info["service_by_design_case"], info["explained_checks"], info["explained_items"]),
                         ("dup", ["status"], [0]))
        # B5b-06: E per item: {check, detail_matches} explains the matching item only; the line records its index,
        # and merge_run leaves out that item alone (not every item of the name)
        item_cases = dict(cases, res={"id": "res", "agreement": "by_design", "direct_expectation": {"http": 200},
                                      "service_answer": {"outcome": "answered"},
                                      "explains_checks": [{"check": "resolve_runs", "detail_matches": "merged .*"}]})
        two = [st("pass"), dict(bad, detail="merged runs"), dict(bad, detail="wrong run 3")]
        info = {}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", by(["res"]), ok, {"outcome": "fail"}, two, False,
                                               item_cases, None, info)[0], "fail")
        self.assertEqual((info["explained_checks"], info["explained_items"]), (["resolve_runs"], [1]))
        row = {"arm": "mcp", "outcome": "pass", "checks": two, "explained_checks": ["resolve_runs"],
               "explained_items": [1], "reevaluation": {}}
        merge_run.stricter_outcomes([row], {}, sc.outcome_from_checks)
        self.assertEqual(row["outcome"], "fail")   # the unexplained item still fails
        row = dict(row, outcome="pass", checks=two[:2] + [dict(runs)])
        merge_run.stricter_outcomes([row], {}, sc.outcome_from_checks)
        self.assertEqual(row["outcome"], "pass")
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
