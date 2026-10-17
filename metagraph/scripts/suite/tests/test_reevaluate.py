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


def slot_doc(seconds, tries="http 1 fail 1"):
    """A pattern search's raw record whose deadline-ended task's held-slot line states one failed HTTP try of
    |seconds| (the failing attempt's own wall the wall rule reads)."""
    slot = "\u29d7 Slot held for search s %.1fs (att 1, ok, acq 0.1, rel 0.1, %s %.3fs backoff 0.0s, released)" % (
        seconds, tries, seconds)
    return {"calls": [{"tool": "GET /searches/s/diagnostic-log", "payload": {"regular_tasks": [
        {"error": "pattern_deadline: late", "diagnostic_log": [slot]}]}}]}


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

    def probe_checks(self):
        """A checks.py of its own directory whose one check states what it was passed (previous, client)."""
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
        return d

    def test_checks_of_another_directory(self):
        # a checks.py of its own directory: it is passed what the runners pass, and the relations use it too
        d = self.probe_checks()
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

    def test_unanswered_repeat_feeds_no_previous(self):
        # an earlier repeat the host never answered (no HTTP status) or its runner ended in error is no earlier repeat
        # in `previous` (suite_common.feeds_previous, as run_direct.py applies it): its digest is a stand-in's, and
        # the next repeat's deterministic would fail although the host answered nothing
        ok = {"http": 200, "body": COMPLETE, "error": None}
        added = {("direct-local", 0): ("error", {"http": 0, "body": None,
                                                 "error": {"code": "transport", "text": "timed out"}}),
                 ("direct-local", 1): ("pass", ok), ("direct-remote", 0): ("inconclusive", DEADLINE),
                 ("direct-remote", 1): ("pass", ok)}
        for (arm, rep), (outcome, answer) in sorted(added.items()):
            raw = sc.raw_name("D-0001", arm, rep)
            sc.append_jsonl(os.path.join(self.run_dir, "arms", arm, "executions.jsonl"), {
                "run_id": "R", "scenario": "D-0001", "area": "D", "arm": arm, "repeat": rep, "route": "POST /pattern",
                "started_utc": sc.iso_utc(T0 + 100 + rep), "ended_utc": sc.iso_utc(T0 + 101 + rep), "wall_ms": 50.0,
                "http": answer["http"], "outcome": outcome, "checks": [], "answer_digest": "stand-in", "raw": raw,
                "error": answer["error"]})
            sc.write_gz_json(os.path.join(self.run_dir, "arms", arm, raw),
                             {"scenario": "D-0001", "arm": arm, "repeat": rep, "answer": answer})
        # the MCP arm's r0: answered (http 200), ended in error by its runner
        path = os.path.join(self.run_dir, "arms", "mcp", "executions.jsonl")
        rows, _ = sc.read_jsonl(path)
        with open(path, "w") as fh:
            for r in rows:
                if (r["scenario"], r["repeat"]) == ("D-0001", 0):
                    r["outcome"] = "error"
                fh.write(json.dumps(r) + "\n")
        reevaluate.reevaluate(self.run_dir, self.probe_checks())
        previous = dict((arm, json.loads(self.lines(arm)[("D-0001", 1)]["checks"][0]["detail"])["previous"])
                        for arm in ("direct-local", "direct-remote", "mcp"))
        self.assertEqual(previous, {"direct-local": [], "mcp": [],
                                    "direct-remote": [{"digest": "D503", "time_limited": ""}]})

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

        def task(seconds, tries="http 1 fail 1"):
            # the deadline-ended task's held-slot line: its one failed HTTP try is the failing attempt's own wall
            slot = "\u29d7 Slot held for search s %.1fs (att 1, ok, acq 0.1, rel 61.1, %s %.3fs backoff 0.0s, " \
                   "released)" % (seconds, tries, seconds)
            return {"calls": [{"tool": "GET /searches/s/diagnostic-log", "payload": {"regular_tasks": [
                {"error": "pattern_deadline: late", "diagnostic_log": ["x", slot]}]}}]}
        info = {}
        out, kept = reevaluate.outcome_of(None, "mcp", sc_, relayed, {"outcome": "fail"}, fail, False, cases, None,
                                          info, task(60.055))
        self.assertEqual((out, kept), ("service_by_design",
                                       "the table's answer-triggered case pattern.deadline_ends_task"))
        self.assertIn(">= time_budget_ms 60000 - 1000 ms = 59000 ms", info["verdict_detail"])
        # R1 (wall_rule budget_minus_reserve): a known wall below the request's budget - 1 s holds the case back (the
        # checks judge): the engine answers its deadline only once the budget has elapsed
        self.assertEqual(reevaluate.outcome_of(None, "mcp", sc_, relayed, {"outcome": "fail"}, fail, False, cases,
                                               None, None, task(0.4))[0], "fail")
        self.assertEqual(reevaluate.outcome_of(None, "mcp", sc_, relayed, {"outcome": "fail"}, fail, False, cases,
                                               None, None, task(59.5))[0], "service_by_design")
        info = {}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", sc_, relayed, {"outcome": "fail"}, fail, False, cases,
                                               None, info, task(58.9))[0], "fail")
        self.assertIn("before time_budget_ms 60000 - 1000 ms = 59000 ms: the case does not apply",
                      info["verdict_detail"])
        # a reserve the answer states above the host's (host_capabilities pattern.finalize_reserve_ms 250) is no
        # host deadline, however long the wall (a 60,000 ms reserve would have made any wall plausible)
        big = dict(relayed, error=dict(relayed["error"], text=relayed["error"]["text"].replace(
            "reserve of 250 ms", "reserve of 60000 ms")))
        info = {}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", sc_, big, {"outcome": "fail"}, fail, False, cases, None,
                                               info, task(60.055))[0], "fail")
        self.assertIn("stated finalisation reserve 60000 ms above the host's 250 ms: the case does not apply",
                      info["verdict_detail"])
        # the failing attempt's own wall, never mcp_phases.backend_ms (a pattern search's is the sum over its tasks'
        # calls, a retry's included): with no usable slot line rule (c) holds, whatever backend_ms says
        info = {}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", sc_, relayed, at(61000.0), fail, False, cases, None,
                                               info)[0], "inconclusive")
        self.assertEqual((info["inconclusive_reason"], info["service_by_design_case"]),
                         ("no service-side timing: plausibility undecidable", "pattern.deadline_ends_task"))
        self.assertEqual(reevaluate.outcome_of(None, "mcp", sc_, relayed, at(400.0), fail, False, cases)[0],
                         "inconclusive")
        # rule (c): no service-side wall, or no request budget: inconclusive, the case recorded
        info = {}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", sc_, relayed, {"outcome": "fail", "wall_ms": 61000.0},
                                               fail, False, cases, None, info)[0], "inconclusive")
        self.assertEqual((info["inconclusive_reason"], info["service_by_design_case"]),
                         ("no service-side timing: plausibility undecidable", "pattern.deadline_ends_task"))
        info = {}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", dict(sc_, request={}), relayed, {"outcome": "fail"},
                                               fail, False, cases, None, info, task(61.0))[0], "inconclusive")
        self.assertEqual(info["inconclusive_reason"], "no threshold known: plausibility undecidable")
        finding = fail + [{"check": "budget_respected", "outcome": "finding", "detail": ""}]
        self.assertEqual(reevaluate.outcome_of(None, "mcp", sc_, relayed, {"outcome": "fail"}, finding, False,
                                               cases)[0], "finding")
        # a task whose slot shows more than one try gives no wall (a sum over a retry)
        self.assertEqual(reevaluate.outcome_of(None, "mcp", sc_, relayed, {"outcome": "fail"}, fail, False, cases,
                                               None, None, task(61.0, "http 2 fail 1"))[0], "inconclusive")
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
        # no budget and no stated timeout: the service's client timeout is the term (20 s, the mode's base:
        # legacy_unbudgeted_client_timeout_term, pending the async session's agreement), so the threshold is
        # min(10 s, 20 s) - 1 s = 9,000 ms, not rule (c)
        bare = dict(to, error={"code": "resolve_timeout", "text": "the backend did not answer"})
        info = {}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", scenario("R-0011", "POST /resolve", STATUS), bare,
                                               at(30000.0), fail, False, cases, None, info)[0], "inconclusive")
        self.assertEqual(info["inconclusive_reason"], "resolve_timeout")
        self.assertIn(">= min(the write cap 10000 ms, the client timeout 20000 ms) - 1000 ms = 9000 ms",
                      info["verdict_detail"])
        self.assertEqual(reevaluate.outcome_of(None, "mcp", scenario("R-0011", "POST /resolve", STATUS), bare,
                                               at(8000.0), fail, False, cases)[0], "fail")
        # an unbudgeted one stating its timeout: the larger term, capped by the write phase
        self.assertEqual(reevaluate.outcome_of(None, "mcp", scenario("R-0011", "POST /resolve", STATUS), to,
                                               at(29300.0), fail, False, cases)[0], "inconclusive")
        # B5b-03: an undecidable host state is inconclusive whatever the checks say, as in async's runner (lane S5b's
        # addendum-11 divergence, a finding there, is dropped)
        info = {}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", res, to, line, finding, False, cases, None, info)[0],
                         "inconclusive")
        self.assertEqual(info["inconclusive_reason"], "no service-side timing: plausibility undecidable")
        # G22: the relayed /resolve deadline: the budget the service sent (the request's) - 1 s of slack (SD-020: the
        # host answers its deadline only once the budget has elapsed; the reserve ends the work, not the answer)
        rdl = {"http": 503, "body": None, "error": {"code": "traversal_resolve_deadline", "text": "the backend could "
                                                    "not build and write its answer within bounds.time_budget_ms "
                                                    "(20000 ms, the finalisation reserve of 250 ms included)"}}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", res, rdl, at(22522.2), fail, False, cases)[0],
                         "inconclusive")
        self.assertEqual(reevaluate.outcome_of(None, "mcp", res, rdl, at(19700.0), fail, False, cases)[0],
                         "inconclusive")
        self.assertEqual(reevaluate.outcome_of(None, "mcp", res, rdl, at(18900.0), fail, False, cases)[0], "fail")
        # the tool error's budget is checked against the request's, never read: the service refuses a budget above
        # the host's maximum and never lowers one, so a tool error stating 5000 ms for a request of 20000 ms is a
        # service fault, not a host state (bounds only from the request, the host and the service settings)
        te = {"code": "traversal_resolve_deadline", "time_budget_ms": 5000, "finalize_reserve_ms": 250}
        info = {}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", res, rdl, at(5000.0), fail, False, cases, te, info)[0],
                         "fail")
        self.assertEqual(info["verdict_detail"], "the tool error's time_budget_ms 5000 ms is not the request's 20000 ms")
        te = {"code": "traversal_resolve_deadline", "time_budget_ms": 20000, "finalize_reserve_ms": 19999}
        info = {}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", res, rdl, at(5.0), fail, False, cases, te, info)[0],
                         "fail")
        self.assertEqual(info["verdict_detail"], "stated finalisation reserve 19999 ms above the host's 250 ms")
        te = {"code": "traversal_resolve_deadline", "time_budget_ms": 20000, "finalize_reserve_ms": 250}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", res, rdl, at(19800.0), fail, False, cases, te)[0],
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
                     "threshold 59000 ms: provably premature")
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
        for ans, ms, least in ((rdl, 0.030, 19000), (to, 0.040, 9000)):
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
            # N-realrep-0002 r0's shape (22,522 ms >= 19,000 ms): undecidable as recorded
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
        self.assertEqual(mark["mcp_auth"], {"assumed": None, "recorded": {"unknown": 1}, "assumed_lines": 0,
                                            "run_json": []})
        self.assertNotIn("mcp_auth_assumed", row["reevaluation"])
        self.assertEqual(reevaluate.main([run, "--quiet", "--mcp-auth", "anonymous"]), 0)
        mark = sc.read_json(os.path.join(run, "reevaluated.json"))
        row = sc.read_jsonl(os.path.join(run, "arms", "mcp", "executions.reevaluated.jsonl"))[0][0]
        self.assertEqual(mark["mcp_auth"], {"assumed": "anonymous", "recorded": {"unknown": 1}, "assumed_lines": 1,
                                            "run_json": []})
        self.assertEqual((row["outcome"], row["reevaluation"]["mcp_auth_assumed"]), ("fail", "anonymous"))
        self.assertIn("client wall 5 ms + 1 s slack < threshold 59000 ms", row["verdict_detail"])
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
        # lane S5e, where the runner records it (2a8acbaf on): mcp_auth on every line and in the raw record, the
        # client's in run.json (each invocation's). A recorded null (a client without mcp_auth) is unknown to the
        # runner: no assumption applies; run.json's values are in the mark and a line's other value is a warning
        d = os.path.join(run, "arms", "mcp")
        rows = sc.read_jsonl(os.path.join(d, "executions.jsonl"))[0]
        with open(os.path.join(d, "executions.jsonl"), "w") as fh:
            fh.write(json.dumps(dict(rows[0], mcp_auth=None)) + "\n")
        sc.write_json_atomic(os.path.join(d, "run.json"), {"run_id": "R", "arm": "mcp", "checks_module": "checks",
                                                           "mcp_auth": "anonymous",
                                                           "invocations": [{"mcp_auth": "api_key"}]})
        mark = reevaluate.reevaluate(run, mcp_auth="anonymous")
        row = sc.read_jsonl(os.path.join(run, "arms", "mcp", "executions.reevaluated.jsonl"))[0][0]
        self.assertEqual(mark["mcp_auth"], {"assumed": "anonymous", "recorded": {"unknown": 1}, "assumed_lines": 0,
                                            "run_json": ["anonymous", "api_key"]})
        self.assertEqual((row["outcome"], row["inconclusive_reason"]), ("inconclusive", sc.UNDECIDABLE_NO_WALL))
        self.assertNotIn("mcp_auth_assumed", row["reevaluation"])
        self.assertEqual(sc.mcp_auth_of(dict(row, reevaluation={"mcp_auth_assumed": "anonymous"})), (None, None))
        self.assertTrue(sc.mcp_auth_recorded({}, {"mcp_auth": None}))
        self.assertFalse(sc.mcp_auth_recorded({"outcome": "fail"}, None))
        with open(os.path.join(d, "executions.jsonl"), "w") as fh:
            fh.write(json.dumps(dict(rows[0], mcp_auth="oauth")) + "\n")
        mark = reevaluate.reevaluate(run)
        self.assertEqual(mark["mcp_auth"]["recorded"], {"oauth": 1})
        self.assertIn("mcp_auth 'oauth' is none of run.json's ['anonymous', 'api_key']", mark["warnings"][-1])

    def test_discovery_without_budget(self):
        # A2-01 MCP side with addendum 1 item 1: a relayed /resolve deadline on a discovery that names no
        # time_budget_ms is implausible (the service sends no deadline for a discovery; SPEC-core §4.5 opt-in):
        # the checks judge (status fails); an explicit-labels resolve without a caller budget has the service's default
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
        # explicit labels, no caller budget: the service sends its default (service_settings.resolve
        # .default_time_budget_ms held to the host's resolve.time_budget.max_time_ms: 60000 ms), so the threshold is
        # known; a tool error stating another budget is no host state (a relayed no-budget deadline must state the
        # service's default)
        labels = scenario("R-0001", "POST /resolve", STATUS, mcp={"tool": "traverse_resolve",
                                                                  "arguments": {"labels": ["562"], "sequence": "AC"}})
        self.assertFalse(sc.discovery_without_budget(labels))
        self.assertEqual(sc.resolve_budget(labels), (60000.0, "the service's default"))
        at = {"outcome": "fail", "mcp_phases": {"backend_ms": 20000.0}}
        info = {}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", labels, rdl, at, fail, False, cases,
                                               {"code": "traversal_resolve_deadline", "time_budget_ms": 20000},
                                               info)[0], "fail")
        self.assertEqual(info["verdict_detail"],
                         "the tool error's time_budget_ms 20000 ms is not the service's default 60000 ms")
        info = {}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", labels, rdl, at, fail, False, cases, None, info)[0],
                         "fail")
        self.assertIn("before time_budget_ms - 1000 ms = 59000 ms: not a host state", info["verdict_detail"])
        self.assertEqual(reevaluate.outcome_of(None, "mcp", labels, rdl, {"outcome": "fail", "mcp_phases": {
            "backend_ms": 60010.0}}, fail, False, cases, {"code": "traversal_resolve_deadline",
                                                          "time_budget_ms": 60000})[0], "inconclusive")

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
        # no budget: the client timeout (20 s) is a term beside the stated 5 s, so a stated "within 0.001 s" or "5 s"
        # no longer lowers the threshold below min(10 s, t) - 1 s (pending the async session's agreement)
        self.assertEqual(sc.resolve_timeout_least(bare, short["error"], None), (9000.0, (
            "min(the write cap 10000 ms, max(the client timeout 20000 ms, the stated timeout 5 s)) - 1000 ms")))
        tiny = {"code": "resolve_timeout", "text": "the backend did not answer within 0.001 s"}
        self.assertEqual(sc.resolve_timeout_least(scenario("R-d", "POST /resolve", STATUS, mcp={
            "tool": "traverse_resolve", "arguments": {"discover": True}}), tiny, None)[0], 9000.0)
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
                             "status": 504, "detail": {"code": "other", "timeout_phase": "write", "bound_s": 10}}}}]}
        self.assertEqual(sc.timeout_statement(job, to["error"], None)["source"], "job error")
        info = {}
        self.assertEqual(reevaluate.outcome_of(None, "mcp", res, to, at(9500.0), fail, False, cases, None, info,
                                               job)[0], "inconclusive")
        self.assertIn("the stated bound 10 s (write) - 1000 ms = 9000 ms", info["verdict_detail"])
        # the write bound is the client's write phase, min(service_settings' write timeout 10 s, t): every budgeted
        # call has t >= 10 s, so a stated write bound of 5 s (or 0.001 s) is no bound the service applies, whatever
        # the wall: the checks judge
        for bound in (5, 0.001):
            short_write = dict(write, bound_s=bound)
            info = {}
            self.assertEqual(reevaluate.outcome_of(None, "mcp", res, to, at(9500.0), fail, False, cases, short_write,
                                                   info)[0], "fail", bound)
            self.assertEqual(info["verdict_detail"], "service write bound %g s is not its write cap 10 s" % bound)
        # within the 1 ms the service rounds the statement to
        self.assertEqual(reevaluate.outcome_of(None, "mcp", res, to, at(9500.0), fail, False, cases,
                                               dict(write, bound_s=10.001))[0], "inconclusive")

    def test_resolve_timeout_read_floor(self):
        # rule 4's read floor (async's f52846fb): a stated READ bound below the service's own read timeout t - 1 s
        # (t = max(the mode's base, budget + the host's reserve + 10 s), the budget the request's, an explicit-labels
        # resolve's without one the service's default) is no host state whatever the wall (an unknown one
        # included): the checks judge, verdict_detail "service read bound {b} s below its own floor {f} s" (%g);
        # the bounds come from the request, the host's capabilities and the service's settings, never from the
        # tool error, whose budget and reserve are checked; a legacy record is unchanged
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
        # a reserve the tool error states is checked against the host's (resolve.time_budget.finalize_reserve_ms
        # 250, host_capabilities), never read: one above it is no bound the service applied
        self.assertEqual(judge(20000, stated("read", 30.25, finalize_reserve_ms=2000.0), 30250.0),
                         ("fail", "stated finalisation reserve 2000 ms above the host's 250 ms", None))
        # no budget in the arguments: an explicit-labels resolve gets the service's default (60000 ms: t = 70.25 s),
        # a tool error stating another budget is no host state, and the floor holds whatever the tool error says
        self.assertEqual(judge(None, stated("read", 20.25, time_budget_ms=20000), 20300.0),
                         ("fail", "the tool error's time_budget_ms 20000 ms is not the service's default 60000 ms",
                          None))
        self.assertEqual(judge(None, stated("read", 20.0), 19000.0),
                         ("fail", "service read bound 20 s below its own floor 69.25 s", None))
        self.assertEqual(judge(None, stated("read", 70.25), 70300.0)[::2], ("inconclusive", "resolve_timeout"))
        # a discovery without a budget: the discovery timeout (service_settings.resolve.discover_timeout_s 120 s)
        disc = scenario("R-floor", "POST /resolve", STATUS, mcp={"tool": "traverse_resolve",
                                                                 "arguments": {"discover": True}})
        self.assertEqual(sc.read_bound_below_floor(disc, to["error"], stated("read", 20.25)),
                         "service read bound 20.25 s below its own floor 119 s")
        self.assertIsNone(sc.read_bound_below_floor(disc, to["error"], stated("read", 120.0)))
        # the write phase: min(the write timeout, t) = 10 s on a budgeted call, then bound_s - 1 s
        self.assertIsNone(sc.read_bound_below_floor(res(20000), to["error"], stated("write", 10.0)))
        self.assertEqual(judge(20000, stated("write", 10.0), 9000.0)[::2], ("inconclusive", "resolve_timeout"))
        # the deadline window: the exchange budget (t plus each phase's cap) less what preceded the exchange
        # (acquire, claim); below the read floor t - 1 s the attempt's span places its origin, and without one
        # nothing does (undecidable: test_resolve_timeout_deadline_window)
        self.assertEqual(judge(20000, stated("deadline", 15.0), 45000.0),
                         ("inconclusive", sc.UNDECIDABLE_NO_ORIGIN + "; deadline window 15 s below its floor 29.25 s",
                          sc.UNDECIDABLE_NO_ORIGIN))
        self.assertEqual(judge(20000, stated("deadline", 30.0), 29000.0)[::2], ("inconclusive", "resolve_timeout"))
        self.assertEqual(judge(20000, stated("deadline", 30.0), 28000.0)[0], "fail")
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

    def test_resolve_timeout_deadline_window(self):
        # GPT review 2's B7 (VERIFY-G): the service takes the exchange's deadline (now + B, B = t + the connect,
        # write and pool caps = 55.25 s for t = 30.25 s) before it acquires the host's token (at most A = 2 s) and
        # claims the attempt, so a slow acquire or claim legitimately leaves a deadline window below the read floor
        # t - 1 s = 29.25 s. Below the floor the failing attempt's span (from its grant to its settling) places the
        # deadline's origin: shorter than B - A - 1 s = 52.25 s, the deadline began before the acquire (fail);
        # otherwise load (inconclusive); no span: undecidable. The service's stated pre_exchange_ms / acquire_ms
        # (agreed 2026-10-15) are preferred where present: bound_s + pre_exchange_ms = B within 2 ms
        cases, _ = cmp.load_by_design(cmp.DEFAULT_BY_DESIGN)
        fail = [{"check": "status", "outcome": "fail", "detail": "http 504"}]
        res = scenario("R-pos-0001", "POST /resolve", STATUS, mcp={"tool": "traverse_resolve",
                                                                   "arguments": {"time_budget_ms": 20000}})
        to = {"http": 504, "body": None, "error": {"code": "resolve_timeout", "text": "the backend's answer did not "
                                                   "complete within the exchange window (deadline)"}}
        self.assertEqual((sc.resolve_client_timeout_s(res), sc.exchange_budget_s(res)), (30.25, 55.25))

        def te(bound, **extra):
            return dict({"code": "resolve_timeout", "timeout_phase": "deadline", "bound_s": bound,
                         "elapsed_ms": bound * 1000.0}, **extra)

        def doc(span_ms=None, job_ms=None, admin=None, evidence=None):
            d = {"phase_sources": {"path": "job", "attempts": [{"attempt_id": "x:1", "sent": False, "span_ms": span_ms,
                                                                 "elapsed_ms": None, "http_status": None}]}}
            if job_ms is not None:
                d["phase_sources"].update(started_at="2026-10-15T00:00:00+00:00",
                                          finished_at=sc.iso_utc(sc.parse_utc("2026-10-15T00:00:00Z") + job_ms / 1e3))
            if admin is not None:
                d["calls"] = [{"tool": "GET /admin/traversal/jobs/j", "payload": {"attempts": admin}}]
            if evidence is not None:
                d["deadline_evidence"] = evidence
            return d

        def judge(tool_error, raw=None, wall_ms=None):
            info = {}
            line = {"outcome": "fail", "mcp_phases": {"backend_ms": wall_ms}}
            out = reevaluate.outcome_of(None, "mcp", res, to, line, fail, False, cases, tool_error, info, raw)[0]
            return out, info.get("verdict_detail"), info.get("inconclusive_reason")
        # the full window: the floor holds, the stated bound - 1 s decides (unchanged)
        self.assertEqual(judge(te(55.14), doc(55290.0))[::2], ("inconclusive", "resolve_timeout"))
        # legitimate: a claim 27.5 s after the grant, an acquire of 2 s and a claim of 26 s (span from the attempt),
        # a one-attempt job's span only: load
        for raw in (doc(54800.0), doc(53300.0), doc(None, job_ms=54800.0)):
            self.assertEqual(judge(te(27.25), raw)[::2], ("inconclusive", "resolve_timeout"), raw)
        self.assertEqual(sc.deadline_window_verdict(res, te(27.25), {"span_ms": 54800.0, "span_source": "attempt"}), (
            None, "deadline window 27.25 s below its floor 29.25 s, shortened before the exchange: attempt span 54800 "
                  "ms (attempt) >= exchange budget 55.25 s - acquire 2 s - 1 s = 52250 ms (load)"))
        # no span: undecidable, inconclusive by rule (c)'s class
        self.assertEqual(judge(te(27.25), doc(None)), (
            "inconclusive", sc.UNDECIDABLE_NO_ORIGIN + "; deadline window 27.25 s below its floor 29.25 s",
            sc.UNDECIDABLE_NO_ORIGIN))
        row = {"arm": "mcp", "outcome": "inconclusive", "inconclusive_reason": sc.UNDECIDABLE_NO_ORIGIN}
        self.assertTrue(merge_run.host_state_inconclusive(row, res))
        # defects: the origin at the job's creation 28 s before the grant, an exchange budget of 20 s, a tiny window
        began = "the deadline began before the acquire: attempt span %d ms (attempt) < exchange budget 55.25 s - " \
                "acquire 2 s - 1 s = 52250 ms"
        for bound, span in ((27.14, 27290), (19.89, 20040), (0.001, 150)):
            self.assertEqual(judge(te(bound), doc(float(span)), wall_ms=bound * 1000.0 + 5.0),
                             ("fail", began % span, None), bound)
        # the attempt's exchange budget the record states (its reserved.bound_ms) is checked, never believed
        admin = [{"attempt_id": "x:1", "dispatched_at": "2026-10-15T00:00:00.000+00:00",
                  "settled_at": "2026-10-15T00:00:54.800+00:00", "reserved": {"bound_ms": 55250, "timeout_s": 30.25}}]
        self.assertEqual(judge(te(27.25), doc(None, admin=admin), 27300.0)[::2], ("inconclusive", "resolve_timeout"))
        self.assertEqual(sc.deadline_evidence_from_record(doc(None, admin=admin)), {
            "exchange_budget_s": 55.25, "timeout_s": 30.25, "acquire_timeout_s": 2.0,
            "granted_at": "2026-10-15T00:00:00.000+00:00", "settled_at": "2026-10-15T00:00:54.800+00:00",
            "span_ms": 54800.0, "span_source": "attempt"})
        admin[0]["reserved"]["bound_ms"] = 50000
        self.assertEqual(judge(te(27.25), doc(None, admin=admin)), (
            "fail", "the attempt's exchange budget 50 s is not t + the phase caps 55.25 s", None))
        # a window above the exchange budget is none the service applied
        self.assertEqual(judge(te(56.0), doc(56100.0))[:2],
                         ("fail", "service deadline window 56 s above its exchange budget 55.25 s"))
        # the runner's deadline_evidence is read (cross-checked with the record: a different span is a warning)
        ev = {"exchange_budget_s": 55.25, "timeout_s": 30.25, "acquire_timeout_s": 2.0, "granted_at": "g",
              "settled_at": "s", "span_ms": 20000.0, "span_source": "attempt"}
        self.assertEqual(judge(te(27.25), doc(None, evidence=ev))[0], "fail")
        self.assertEqual(judge(te(27.25), doc(None, evidence=dict(ev, span_ms=54800.0)))[0], "inconclusive")
        warnings = []
        self.assertEqual(sc.deadline_evidence({}, doc(54800.0, evidence=ev), None, warnings)["span_ms"], 20000.0)
        self.assertEqual(warnings, ["the raw record's deadline_evidence states span_ms 20000.0, the record gives "
                                    "54800.0: the recorded one is read"])
        self.assertEqual(sc.deadline_evidence({"deadline_evidence": dict(ev, span_ms=54800.0)}, None)["span_ms"],
                         54800.0)
        # the stated origin (pre_exchange_ms, acquire_ms) is preferred: consistent, load without any span
        self.assertEqual(judge(te(27.25, pre_exchange_ms=28000, acquire_ms=1500), doc(None), 27300.0)[::2],
                         ("inconclusive", "resolve_timeout"))
        self.assertEqual(sc.deadline_window_verdict(res, te(27.25, pre_exchange_ms=28000)), (
            None, "deadline window 27.25 s = exchange budget 55.25 s - 28000 ms before the exchange = 27.25 s: "
                  "shortened before the exchange (load)"))
        self.assertEqual(judge(te(27.251, pre_exchange_ms=27999), doc(None), 27300.0)[0], "inconclusive")  # rounding
        # ... and it catches what the floor does not: a window that ended earlier than its terms allow, above the
        # floor too (the origin 5 s early), one longer than they allow, an acquire longer than the time before
        self.assertEqual(judge(te(27.25, pre_exchange_ms=20000), doc(54800.0), 27300.0), (
            "fail", "service deadline window 27.25 s ended earlier than its terms allow: exchange budget 55.25 s - "
                    "20000 ms before the exchange = 35.25 s", None))
        self.assertEqual(judge(te(50.14, pre_exchange_ms=110), doc(50290.0), 50200.0)[0], "fail")
        self.assertEqual(judge(te(55.14, pre_exchange_ms=110), doc(55290.0), 55200.0)[::2],
                         ("inconclusive", "resolve_timeout"))
        self.assertEqual(judge(te(27.25, pre_exchange_ms=30000), doc(None), 27300.0)[1],
                         "service deadline window 27.25 s is longer than its terms allow: exchange budget 55.25 s - "
                         "30000 ms before the exchange = 25.25 s")
        self.assertEqual(judge(te(27.25, pre_exchange_ms=28000, acquire_ms=28500), doc(None), 27300.0)[1],
                         "the stated acquire 28500 ms exceeds the 28000 ms before the exchange it is part of")
        # a self-consistent statement whose origin precedes the acquire: the span still places it; a measured acquire
        # longer than the acquire timeout widens what the span may lack
        self.assertEqual(judge(te(27.14, pre_exchange_ms=28110), doc(27290.0), 27200.0)[1], began % 27290)
        self.assertEqual(judge(te(27.25, pre_exchange_ms=28000, acquire_ms=3500), doc(51800.0), 27300.0)[::2],
                         ("inconclusive", "resolve_timeout"))
        self.assertEqual(judge(te(27.25, pre_exchange_ms=28000), doc(51800.0), 27300.0)[0], "fail")
        # a runner's recorded timeout_terms without the origin terms: they are read from the record's error
        rec = {"timeout_phase": "deadline", "bound_s": 27.25, "elapsed_ms": 27250.0, "source": "tool error"}
        warnings = []
        got = sc.timeout_terms({"timeout_terms": rec}, None, to["error"], te(27.25, pre_exchange_ms=20000), warnings)
        self.assertEqual((got, warnings), (dict(rec, pre_exchange_ms=20000.0), []))
        self.assertEqual(sc.timeout_terms({"timeout_terms": rec}, None, to["error"], te(27.25)), rec)
        self.assertEqual(sc.timeout_statement(None, to["error"], te(27.25, pre_exchange_ms=-1, acquire_ms="5")), rec)

    def test_explained_items_as_the_runner(self):
        # async's f386a3c0: the line (and raw record) records explained_items beside explained_checks, the indices
        # into checks of the results with outcome fail that an entry of the matched cases' union explains; [] when E
        # is empty (a refused / relayed case, an answered case explaining nothing); absent when no case matched
        cases, _ = cmp.load_by_design(cmp.DEFAULT_BY_DESIGN)
        ok = {"http": 200, "body": {}, "error": None}
        dup = scenario("T-edge-0008", "POST /traverse", TRAVERSE, mcp={"tool": "traverse_fetch", "arguments": {},
                                                                     "by_design": ["traverse.duplicate_seeds"]})
        merged = {"check": "invariants", "outcome": "fail", "detail": "1 results for 2 seeds (1 distinct)"}
        checks = [{"check": "status", "outcome": "pass", "detail": ""}, merged,
                  {"check": "traverse_seed_outcome", "outcome": "pass", "detail": ""}]

        def judge(chk, sc_=dup, ans=ok):
            info = {}
            out = reevaluate.outcome_of(None, "mcp", sc_, ans, {"outcome": "fail"}, chk, False, cases, None, info)[0]
            return out, info.get("explained_checks"), info.get("explained_items")
        self.assertEqual(judge(checks), ("service_by_design", ["invariants"], [1]))
        # one invariants failure explained, the other not: fail, one item
        two = checks + [{"check": "invariants", "outcome": "fail", "detail": "result 0 without outcome"}]
        self.assertEqual(judge(two), ("fail", ["invariants"], [1]))
        # the older wording (no distinct count) is no merge for the table's regex: nothing explained, [] recorded
        self.assertEqual(judge([checks[0], dict(merged, detail="1 results for 2 seeds"), checks[2]]),
                         ("fail", [], []))
        whole = [{"check": "x", "outcome": "pass"}, {"check": "status", "outcome": "fail"},
                 {"check": "status", "outcome": "fail"}]
        self.assertEqual(cmp.explained_items(whole, ["status"]), [1, 2])
        self.assertEqual(cmp.explained_items(whole, []), [])
        # no case: nothing recorded
        self.assertEqual(judge(checks, scenario("P-0001", "POST /pattern", STATUS)), ("fail", None, None))
        # end to end: resolve.duplicate_labels (explains status) answered 200 where the host refuses; the runner's
        # recorded list is kept when the replay reproduces it and noted when not
        run = os.path.join(self.tmp.name, "items")
        res = scenario("R-dup", "POST /resolve", [{"check": "status", "args": {"http": 400}}],
                       mcp={"tool": "traverse_resolve", "by_design": ["resolve.duplicate_labels"],
                            "arguments": {"sequence": "A" * 40, "labels": ["562", "562"]}})
        manifest = {"manifest_version": 1, "scenarios": [res], "relations": []}
        d = os.path.join(run, "arms", "mcp")
        os.makedirs(os.path.join(d, "raw"))
        sc.write_json_atomic(os.path.join(run, "manifest.json"), manifest)
        sc.write_json_atomic(os.path.join(d, "manifest.json"), manifest)
        sc.write_json_atomic(os.path.join(d, "run.json"), {"run_id": "R", "arm": "mcp", "checks_module": "checks"})
        for rep, items in ((0, [0]), (1, [1]), (2, None)):
            raw = sc.raw_name(res["id"], "mcp", rep)
            line = {"run_id": "R", "scenario": res["id"], "arm": "mcp", "repeat": rep, "cache": "cold", "http": 200,
                    "outcome": "service_by_design", "checks": [], "mcp_phases": {}, "raw": raw, "error": None,
                    "route": "POST /resolve", "started_utc": sc.iso_utc(T0), "wall_ms": 5.0, "mcp_auth": "anonymous",
                    "service_by_design_case": "resolve.duplicate_labels", "explained_checks": ["status"]}
            if items is not None:
                line["explained_items"] = items
            sc.append_jsonl(os.path.join(d, "executions.jsonl"), line)
            sc.write_gz_json(os.path.join(d, raw), {"answer": {"http": 200, "body": {"results": []}, "error": None},
                                                    "explained_items": items, "mcp_auth": "anonymous"})
        reevaluate.reevaluate(run)
        rows = sc.read_jsonl(os.path.join(d, "executions.reevaluated.jsonl"))[0]
        for row in rows:
            self.assertEqual((row["outcome"], row["service_by_design_case"], row["explained_checks"],
                              row["explained_items"]), ("service_by_design", "resolve.duplicate_labels", ["status"],
                                                        [0]), row["repeat"])
            self.assertEqual([c["outcome"] for c in row["checks"]], ["fail"])
            self.assertEqual(list(row).index("explained_items"), list(row).index("explained_checks") + 1)
        self.assertNotIn("explained_items_recorded", rows[0]["reevaluation"])
        self.assertEqual(rows[1]["reevaluation"]["explained_items_recorded"], [1])
        self.assertNotIn("explained_items_recorded", rows[2]["reevaluation"])   # a runner before f386a3c0

    def test_mcp_outcome_not_implemented(self):
        # async's 901893b7 (the runner's tests, mirrored): a matched case stating answer_match.mcp_outcome
        # not_implemented_mcp makes the execution not_implemented_mcp, the case recorded, explained checks / items [],
        # verdict_detail naming the deployment condition (the service's own text); service_by_design stays the
        # default; neither rule (c) nor a host-state code makes it inconclusive; a finding is never erased; a named
        # match applies it too
        text = ("overrides.labels.max_seed_labels is on the advanced strategy path, which is not enabled on this "
                "deployment")

        def case(cid="test.deployment_gate", route="POST /traverse", **am_extra):
            am = {"routes": [route], "direct": {"http": 200},
                  "mcp": {"http": 400, "code": "preset_unavailable", "field": "overrides.<knob path>",
                          "text_contains": "advanced strategy path, which is not enabled on this deployment"},
                  "direct_mismatch": "disagree", "mcp_outcome": "not_implemented_mcp"}
            am.update(am_extra)
            return {"id": cid, "agreement": "not_implemented_mcp", "direct_expectation": {"http": 200},
                    "service_answer": {"outcome": "refused_at_submit", "code": am["mcp"]["code"],
                                       "text_contains": am["mcp"]["text_contains"], "field": am["mcp"].get("field")},
                    "answer_match": dict((k, v) for k, v in am.items() if v is not None)}
        table = {"test.deployment_gate": case()}
        t4 = scenario("T-labels-0004", "POST /traverse", STATUS, mcp={"tool": "traverse_fetch", "arguments": {}})
        ans = {"http": 400, "body": None, "error": {"code": "preset_unavailable", "text": text,
                                                    "field": "overrides.labels.max_seed_labels"}}
        failing = [{"check": "status", "outcome": "fail", "detail": "http 400, expected 200"},
                   {"check": "invariants", "outcome": "fail", "detail": "no answer to check: http 400"}]

        def judge(sc_, ans_, cases, checks=failing, line=None, doc=None, terr=None):
            info = {}
            out = reevaluate.outcome_of(None, "mcp", sc_, ans_, line or {"outcome": "fail"}, checks, False, cases,
                                        terr, info, doc)[0]
            return out, info
        out, info = judge(t4, ans, table)
        self.assertEqual((out, info["service_by_design_case"], info["explained_checks"], info.get("explained_items"),
                          info.get("inconclusive_reason")),
                         ("not_implemented_mcp", "test.deployment_gate", [], None, None))
        self.assertEqual(info["verdict_detail"], cmp.MCP_NOT_IMPLEMENTED_DETAIL.format(
            case="test.deployment_gate", field="overrides.labels.max_seed_labels", text=text))
        self.assertEqual(info["verdict_detail"], (
            "not implemented via MCP on this deployment (case test.deployment_gate, answer_match.mcp_outcome "
            "not_implemented_mcp): the service refused overrides.labels.max_seed_labels at submit: %r; on a "
            "deployment that enables it the request is accepted and the case does not apply" % text))
        # the tool error's field wins when it states the key (None: "the request"), as the runner's
        self.assertIn("refused overrides.x at submit", judge(t4, ans, table, terr={"code": "preset_unavailable",
                                                                                    "field": "overrides.x"})[1][
            "verdict_detail"])
        nofield = {"n": case("n", mcp={"http": 400, "code": "preset_unavailable", "text_contains": "advanced"})}
        bare = dict(ans, error={"code": "preset_unavailable", "text": text})
        self.assertIn("refused the request at submit", judge(t4, bare, nofield)[1]["verdict_detail"])
        self.assertIn("refused the request at submit", judge(t4, ans, nofield, terr={"code": "preset_unavailable",
                                                                                     "field": None})[1][
            "verdict_detail"])
        # the default class: absent or stated service_by_design
        for default in (None, "service_by_design"):
            self.assertEqual(judge(t4, ans, {"c": case("c", mcp_outcome=default)})[0], "service_by_design")
        # a finding is never erased
        finding = failing + [{"check": "budget_respected", "outcome": "finding", "detail": "34 s on a 15 s budget"}]
        out, info = judge(t4, ans, table, finding)
        self.assertEqual((out, info["service_by_design_case"]), ("finding", "test.deployment_gate"))
        self.assertIn("not implemented via MCP", info["verdict_detail"])
        # rule (c) does not apply: a wall-ruled case whose wall is undecidable still gives its class
        deadline = {"http": 502, "body": None, "error": {"code": "search_failed", "text": "pattern_deadline: pattern: "
                    "the answer could not be written within time_budget_ms (60000 ms, the finalisation reserve of "
                    "250 ms included)"}}
        walled = {"test.walled": case("test.walled", "POST /pattern", wall_rule="budget_minus_reserve",
                                      mcp={"http": 502, "code": "search_failed", "text_contains": "pattern_deadline"})}
        pat = scenario("N-1", "POST /pattern", STATUS, mcp={"tool": "pattern_search",
                                                           "arguments": {"time_budget_ms": 60000}})
        out, info = judge(pat, deadline, walled)
        self.assertEqual((out, info["service_by_design_case"], info.get("inconclusive_reason")),
                         ("not_implemented_mcp", "test.walled", None))
        self.assertTrue(info["verdict_detail"].startswith("not implemented via MCP on this deployment (case "
                                                          "test.walled"), info["verdict_detail"])
        # a known wall below the threshold holds the case back in matching, as before; its note leads the detail of
        # a later case that matches (the failing attempt's own wall: the deadline-ended task's one HTTP try)
        slot = "\u29d7 Slot held for search s 0.1s (att 1, ok, acq 0.0, rel 0.1, http 1 fail 1 0.005s backoff 0.0s, " \
               "released)"
        quick = {"calls": [{"tool": "GET /searches/s/diagnostic-log", "payload": {"regular_tasks": [
            {"error": "pattern_deadline: late", "diagnostic_log": [slot]}]}}]}
        out, info = judge(pat, deadline, walled, doc=quick)
        self.assertEqual((out, info.get("service_by_design_case")), ("fail", None))
        both = dict(walled, **{"test.gate2": case("test.gate2", "POST /pattern",
                                                  mcp={"http": 502, "code": "search_failed",
                                                       "text_contains": "pattern_deadline"})})
        out, info = judge(pat, deadline, both, doc=quick)
        self.assertEqual((out, info["service_by_design_case"]), ("not_implemented_mcp", "test.gate2"))
        self.assertTrue(info["verdict_detail"].startswith("test.walled (wall_rule budget_minus_reserve): after 5 ms"))
        self.assertIn("; not implemented via MCP on this deployment (case test.gate2", info["verdict_detail"])
        # the host-state rules do not apply
        timed_out = {"http": 504, "body": None, "error": {"code": "resolve_timeout", "text": "timed out after 30 s"}}
        gated = {"test.gated_resolve": case("test.gated_resolve", "POST /resolve",
                                            mcp={"http": 504, "code": "resolve_timeout", "text_contains": "timed"})}
        res = scenario("R-1", "POST /resolve", STATUS, mcp={"tool": "traverse_resolve", "arguments": {}})
        job = {"outcome": "fail", "mcp_phases": {"backend_ms": 31000.0}}
        self.assertEqual(judge(res, timed_out, gated, line=job)[0], "not_implemented_mcp")
        self.assertEqual(judge(res, timed_out, {"x": case("x")}, line=job)[0], "inconclusive")   # the host state
        # a named match applies it too (another route)
        named = scenario("X-1", "GET /elsewhere", STATUS, mcp={"tool": "traverse_fetch",
                                                               "by_design": ["test.deployment_gate"]})
        self.assertEqual(judge(named, ans, table)[0], "not_implemented_mcp")
        # the real table's case 131 on v4's T-labels-0004 / -0006 (MCP 400 preset_unavailable on
        # overrides.labels.max_seed_labels, which the scenario sends): not_implemented_mcp; a field outside
        # overrides.* is not the case
        cases, _ = cmp.load_by_design(cmp.DEFAULT_BY_DESIGN)
        t4 = scenario("T-labels-0004", "POST /traverse", STATUS, mcp={"tool": "traverse_fetch", "arguments": {
            "overrides": {"bounds": {"max_extension_bp": 200},
                          "labels": {"max_seed_labels": 10000, "seed_label_kind": "header"}}}})
        for sid in ("T-labels-0004", "T-labels-0006"):
            out, info = judge(dict(t4, id=sid), ans, cases)
            self.assertEqual((out, info["service_by_design_case"]), ("not_implemented_mcp",
                                                                    "traverse.advanced_path_disabled"), sid)
            # the deployment setting it read: the table's service_settings block (the run recorded no live document)
            self.assertEqual(info["request_evidence"],
                             {"traverse.advanced_path_disabled": {"service_settings.advanced_enabled": "table"}})
        self.assertEqual(judge(t4, dict(ans, error=dict(ans["error"], field="preset")), cases)[0], "fail")
        self.assertEqual(judge(t4, {"http": 200, "body": {}, "error": None}, cases)[0], "fail")
        self.assertEqual(judge(dict(t4, route="POST /resolve"), ans, cases)[0], "fail")
        # the case's request (VERIFY-B's mutations): a refusal of a key the scenario never sent, or of a basic
        # override, is no refusal of the advanced path: the checks judge, the request's note in the detail
        unsent = dict(t4, mcp={"tool": "traverse_fetch", "arguments": {}})
        out, info = judge(unsent, ans, cases)
        self.assertEqual((out, info.get("service_by_design_case")), ("fail", None))
        self.assertEqual(info["verdict_detail"], 'traverse.advanced_path_disabled: request predicate {"path": '
                                                 '"{field}", "present": true} does not hold (field '
                                                 'overrides.labels.max_seed_labels)')
        basic = dict(ans, error=dict(ans["error"], field="overrides.bounds.max_extension_bp",
                                     text=text.replace("labels.max_seed_labels", "bounds.max_extension_bp")))
        out, info = judge(t4, basic, cases)
        self.assertEqual(out, "fail")
        self.assertIn('request predicate {"field": "overrides.", "not_in": ["", "bounds.*"', info["verdict_detail"])
        # a deployment whose recorded live document enables the advanced path: the case does not apply (judged with
        # the live value, the mismatch with the table's block recorded)
        live = cmp.ByDesignCases(cases, cases.blocks, sc.References(cases.blocks, {"service": {
            "service_settings": {"advanced_enabled": True}}}))
        out, info = judge(t4, ans, live)
        self.assertEqual(out, "fail")
        self.assertEqual(info["request_evidence"],
                         {"traverse.advanced_path_disabled": {"service_settings.advanced_enabled": "live"}})
        self.assertEqual(info["request_evidence_mismatches"], [{
            "namespace": "service", "ref": "service_settings.advanced_enabled", "live": True, "table": False,
            "case": "traverse.advanced_path_disabled"}])
        # a runner's not_implemented_mcp mark that names its case is replayed with the table (not kept unread): the
        # case reproduces it; a table without the case judges the checks (noted); no table keeps it with its case
        mark = {"outcome": "not_implemented_mcp", "service_by_design_case": "traverse.advanced_path_disabled",
                "explained_checks": [], "explained_items": [], "verdict_detail": "recorded detail"}
        out, info = judge(t4, ans, cases, line=mark)
        self.assertEqual((out, info["service_by_design_case"]), ("not_implemented_mcp",
                                                                "traverse.advanced_path_disabled"))
        out, info = judge(t4, ans, {"x": case("x", "POST /resolve")}, line=mark)
        self.assertEqual(out, "fail")
        self.assertIn("the runner's not_implemented_mcp mark is not reproduced", info["note"])
        info = {}
        out, kept = reevaluate.outcome_of(None, "mcp", t4, ans, mark, failing, False, {}, None, info)
        self.assertEqual((out, info["service_by_design_case"], info["verdict_detail"]),
                         ("not_implemented_mcp", "traverse.advanced_path_disabled", "recorded detail"))
        self.assertIn("no table", kept)
        # the manifest's mcp.not_implemented (no case named) stays the runner's decision
        self.assertEqual(reevaluate.outcome_of(None, "mcp", t4, ans, {"outcome": "not_implemented_mcp"}, failing,
                                               False, cases), ("not_implemented_mcp",
                                                               "decided by the runner without the checks"))

    def test_request_evidence_on_the_line(self):
        # a re-evaluated MCP line records the references its matched cases' requests read (request_evidence, the
        # runner's field), the mark where they came from; a run whose manifest names a case its own arguments do not
        # satisfy is refused (exit 4), as the runner refuses to start
        run = os.path.join(self.tmp.name, "req")
        args = {"database": "refseq33m-experimental", "preset": "walk", "seeds": [{"sequence": "ACGT" * 10}],
                "overrides": {"labels": {"max_seed_labels": 10000}}}
        t4 = scenario("T-labels-0004", "POST /traverse", STATUS, mcp={"tool": "traverse_fetch", "arguments": args})
        err = {"code": "preset_unavailable", "field": "overrides.labels.max_seed_labels",
               "text": "overrides.labels.max_seed_labels is on the advanced strategy path, which is not enabled on "
                       "this deployment"}
        d = os.path.join(run, "arms", "mcp")
        os.makedirs(os.path.join(d, "raw"))
        manifest = {"manifest_version": 4, "scenarios": [t4], "relations": []}
        sc.write_json_atomic(os.path.join(run, "manifest.json"), manifest)
        sc.write_json_atomic(os.path.join(d, "manifest.json"), manifest)
        sc.write_json_atomic(os.path.join(d, "run.json"), {"run_id": "R", "arm": "mcp", "state": "done"})
        raw = sc.raw_name("T-labels-0004", "mcp", 0)
        sc.write_gz_json(os.path.join(d, raw), {"answer": {"http": 400, "body": None, "error": err}})
        sc.append_jsonl(os.path.join(d, "executions.jsonl"), {
            "run_id": "R", "scenario": "T-labels-0004", "area": "T", "arm": "mcp", "repeat": 0, "http": 400,
            "outcome": "fail", "checks": [], "raw": raw, "error": err, "route": "POST /traverse"})
        mark = reevaluate.reevaluate(run)
        row = sc.read_jsonl(os.path.join(d, "executions.reevaluated.jsonl"))[0][0]
        self.assertEqual((row["outcome"], row["service_by_design_case"]),
                         ("not_implemented_mcp", "traverse.advanced_path_disabled"))
        self.assertEqual(row["request_evidence"],
                         {"traverse.advanced_path_disabled": {"service_settings.advanced_enabled": "table"}})
        self.assertNotIn("request_evidence_mismatches", row)
        self.assertEqual(mark["references"]["service"],
                         {"source": "table", "origin": "the table's service_settings block"})
        # the runner's live document (run.json service_capabilities) enabling the advanced path: judged with it
        cases, _ = cmp.load_by_design(cmp.DEFAULT_BY_DESIGN)
        doc = json.loads(json.dumps(cases.blocks["service"]))
        doc["service_settings"]["advanced_enabled"] = True
        sc.write_json_atomic(os.path.join(d, "run.json"), {"run_id": "R", "arm": "mcp", "state": "done",
                                                           "service_capabilities": {"document": doc, "sha256": "x",
                                                                                    "read_at": "t", "base": "b"}})
        mark = reevaluate.reevaluate(run)
        row = sc.read_jsonl(os.path.join(d, "executions.reevaluated.jsonl"))[0][0]
        self.assertEqual((row["outcome"], row.get("service_by_design_case")), ("fail", None))
        self.assertEqual(row["request_evidence"],
                         {"traverse.advanced_path_disabled": {"service_settings.advanced_enabled": "live"}})
        self.assertEqual(row["request_evidence_mismatches"][0]["live"], True)
        self.assertEqual(mark["references"]["mismatches"], [{"namespace": "service", "ref":
                                                             "service_settings.advanced_enabled", "live": True,
                                                             "table": False}])
        # a manifest naming traverse.empty_seeds on a scenario that sends a seed: refused at load
        bad = dict(t4, mcp=dict(t4["mcp"], by_design=["traverse.empty_seeds"]))
        sc.write_json_atomic(os.path.join(run, "manifest.json"), dict(manifest, scenarios=[bad]))
        with contextlib.redirect_stderr(io.StringIO()) as err_out:
            self.assertEqual(reevaluate.main([run, "--quiet"]), 4)
        self.assertIn("scenario T-labels-0004 names case traverse.empty_seeds", err_out.getvalue())

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
        self.assertEqual(reevaluate.outcome_of(None, "mcp", sc_, relayed, line, checks, False, cases, None, None,
                                               slot_doc(125.0)),
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
                     mcp={"tool": "pattern_search", "arguments": {"predicate_scope": "motif"},
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
        # the case's request: predicate_scope "motif" and no predicate (Q-motif-0004's arguments); a scenario sending
        # predicate_scope "context" asks the other case (pattern.predicate_options_noop_without_predicate)
        info = {}
        other = dict(q, mcp=dict(q["mcp"], arguments={"predicate_scope": "context"}))
        self.assertEqual(reevaluate.outcome_of(None, "mcp", other, ans, {"outcome": "fail"}, fail, False, cases, None,
                                               info)[0], "fail")
        self.assertEqual(info["verdict_detail"], 'pattern.predicate_options_without_predicate: request predicate '
                                                 '{"equals": "motif", "path": "predicate_scope"} does not hold (field '
                                                 'predicate_scope)')
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
        self.assertEqual(reevaluate.outcome_of(None, "mcp", unnamed, relayed, at, fail, False, cases, None, None,
                                               slot_doc(60.1)),
                         ("service_by_design", "the table's answer-triggered case pattern.deadline_ends_task"))
        # a relayed deadline whose budget_respected is a finding keeps it (Q1 of VERIFY-REVIEW)
        self.assertEqual(reevaluate.outcome_of(None, "mcp", unnamed, relayed, at,
                                               fail + [{"check": "budget_respected", "outcome": "finding"}], False,
                                               cases, None, None, slot_doc(60.1))[0], "finding")

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
