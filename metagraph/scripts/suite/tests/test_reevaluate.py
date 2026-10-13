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
        # kept: not implemented, by design, a transport error, the MCP runner's host-state rule
        self.assertEqual(mcp[("P-0002", 0)]["outcome"], "not_implemented_mcp")
        self.assertEqual(mcp[("P-0002", 0)]["reevaluation"]["kept"], "no recorded answer")
        self.assertEqual(mcp[("P-0003", 0)]["outcome"], "service_by_design")
        self.assertEqual(remote[("P-0004", 0)]["outcome"], "error")
        self.assertEqual(mcp[("P-0006", 0)]["outcome"], "inconclusive")
        self.assertIn("host-state", mcp[("P-0006", 0)]["reevaluation"]["kept"])
        # the client: a direct arm's wall (9,257 ms on a 3,000 ms budget: a finding), the MCP arm's backend time
        self.assertEqual(local[("R-0001", 0)]["outcome"], "finding")
        self.assertEqual(mcp[("R-0001", 0)]["outcome"], "inconclusive")
        # previous as {digest, time_limited}: the time-limited cold repeat explains the warm one's difference
        self.assertEqual(mcp[("D-0001", 1)]["outcome"], "inconclusive")
        det = [c for c in mcp[("D-0001", 1)]["checks"] if c["check"] == "deterministic"][0]
        self.assertIn("time budget limited it", det["detail"])
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

    def test_merge_without_a_reevaluation(self):
        with self.assertRaises(SystemExit):
            merge_run.merge(self.run_dir, reevaluated=True)


if __name__ == "__main__":
    unittest.main()
