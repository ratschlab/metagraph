"""Offline tests of shadow_run.py (the coordinator's SHADOW RUN pipeline, round 6): the steps and their arguments
(dry run), their order and the run directory the post-run steps read (a recording runner on a synthetic run),
the stop on a failing step with its resume hint, the gate's exit passed through, the declared file's re-pin, and
the expected load of a manifest. Nothing is sent."""
import contextlib
import io
import json
import os
import shutil
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.dirname(HERE)
sys.path.insert(0, SUITE)

import shadow_run  # noqa: E402
import suite_common as sc  # noqa: E402

PIN = "fdff625ef7322e919b28da4b81763acec65591af"


class Recorder(object):
    def __init__(self, fail=None, gate_rc=0):
        self.calls, self.fail, self.gate_rc = [], fail, gate_rc

    def __call__(self, argv, env=None, stdin=None):
        self.calls.append((argv, env))
        name = os.path.basename(argv[1]) if len(argv) > 1 else argv[0]
        if self.fail and self.fail in " ".join(argv):
            return 3
        return self.gate_rc if name == "gate.py" else 0


class ShadowRun(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="shadow-")
        self.out = os.path.join(self.tmp, "out")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_dry_run_names_every_step(self):
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rc = shadow_run.main(["--out", self.out, "--dry-run", "--library-commit", PIN])
        text = buf.getvalue()
        self.assertEqual(rc, 0)
        names = [ln.split()[0] for ln in text.splitlines()[1:] if ln.startswith("  ")]
        self.assertEqual(names, list(shadow_run.STEPS))
        self.assertIn("--stop-scenarios", text)
        self.assertIn("--paired", text)
        self.assertNotIn("--skip-remote", text)                  # all three arms
        self.assertIn("--ways W1,W2,W4,drift", text)
        self.assertIn("--library-commit %s" % PIN, text)
        self.assertIn("declared_v8_shadow.json", text)
        self.assertIn("METAGRAPH_ADMIN_KEY_FILE=", text)          # the run step's admin reads (the stop flows)
        self.assertFalse(os.path.exists(self.out))                # a dry run writes nothing

    def test_the_post_run_steps_read_the_run_and_the_gate_rc_passes_through(self):
        run = os.path.join(self.out, "runs", "20261021T000000Z_refseq97-x")
        os.makedirs(run)
        rec = Recorder(gate_rc=1)
        rc = shadow_run.main(["--out", self.out, "--from", "local", "--library-commit", PIN], runner=rec)
        self.assertEqual(rc, 1)                                   # GATE FAIL passes through
        tools = [os.path.basename(argv[1]) for argv, _ in rec.calls]
        self.assertEqual(tools, ["run_graphlet_local.py", "graphlet_oracle.py", "compare_graphlets.py", "report.py",
                                 "gate.py"])
        for argv, _env in rec.calls:
            self.assertIn(run, argv)                              # the newest run under DIR/runs
        local = rec.calls[0][0]
        self.assertEqual(local[local.index("--manifest") + 1], os.path.join(run, "manifest.json"))
        self.assertEqual(local[local.index("--library-root") + 1], os.path.join(self.out, "library-fdff625ef732"))
        gate = rec.calls[-1][0]
        self.assertEqual(gate[-1], os.path.join(self.out, "policy_v8_shadow.json"))
        with open(os.path.join(self.out, "shadow.log")) as fh:
            log = fh.read()
        self.assertIn("== gate exit 1", log)

    def test_a_failing_step_stops_with_its_resume_hint(self):
        run = os.path.join(self.out, "runs", "r1")
        os.makedirs(run)
        rec = Recorder(fail="graphlet_oracle.py")
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            rc = shadow_run.main(["--out", self.out, "--from", "local", "--run-dir", run, "--library-commit", PIN],
                                 runner=rec)
        self.assertEqual(rc, 3)
        self.assertEqual([os.path.basename(a[1]) for a, _ in rec.calls], ["run_graphlet_local.py",
                                                                          "graphlet_oracle.py"])
        self.assertIn("resume with --from oracle", err.getvalue())
        # no run directory for the post-run steps: refused before anything runs
        rc = shadow_run.main(["--out", os.path.join(self.tmp, "empty"), "--from", "compare", "--library-commit",
                              PIN], runner=Recorder())
        self.assertEqual(rc, 2)

    def test_the_run_step_sends_all_three_arms_with_async_s_runner(self):
        rec = Recorder()
        shadow_run.main(["--out", self.out, "--from", "run", "--to", "run", "--library-commit", PIN], runner=rec)
        argv, env = rec.calls[0]
        self.assertTrue(argv[0].endswith("run_suite.sh"))
        self.assertIn("--paired", argv)
        self.assertNotIn("--skip-remote", argv)
        self.assertEqual(argv[argv.index("--mcp-cmd") + 1], "default")   # run_suite.sh's ASYNC_MCP_CMD
        self.assertEqual(env["SUITE_RUNS"], os.path.join(self.out, "runs"))

    def test_resume_sets_earlier_records_aside(self):
        # R4-3: --from local / --from oracle after a partial run: the earlier records move aside, logged
        run = os.path.join(self.out, "runs", "r1")
        for d in ("dl-lib", "drift-lib", "oracle", "mcp"):
            os.makedirs(os.path.join(run, "arms", d))
        rec = Recorder()
        shadow_run.main(["--out", self.out, "--from", "local", "--to", "oracle", "--library-commit", PIN], runner=rec)
        arms = sorted(os.listdir(os.path.join(run, "arms")))
        self.assertIn("mcp", arms)                                   # what the run step wrote stays
        for d in ("dl-lib", "drift-lib", "oracle"):
            self.assertNotIn(d, arms)
            self.assertTrue(any(x.startswith(d + ".prev-") for x in arms), arms)
        with open(os.path.join(self.out, "shadow.log")) as fh:
            log = fh.read()
        self.assertIn("local: the earlier records set aside", log)
        self.assertIn("oracle: the earlier records set aside", log)
        self.assertRegex(log, r"local \(code [0-9a-f]{40}(\+dirty)?\)")   # R4-4: the code's commit per step

    def test_the_library_is_reused_only_when_whole(self):
        lib = os.path.join(self.out, "library-" + PIN[:12])
        os.makedirs(os.path.join(lib, "metagraph", "api", "python"))
        rec = Recorder()
        shadow_run.main(["--out", self.out, "--from", "library", "--to", "library", "--library-commit", PIN],
                        runner=rec)
        self.assertEqual(len(rec.calls), 1)                          # no mark: extracted again
        with open(os.path.join(lib, shadow_run.LIBRARY_MARK), "w") as fh:
            fh.write(PIN + "\n")
        rec = Recorder()
        shadow_run.main(["--out", self.out, "--from", "library", "--to", "library", "--library-commit", PIN],
                        runner=rec)
        self.assertEqual(rec.calls, [])
        with open(os.path.join(self.out, "shadow.log")) as fh:
            self.assertIn("library: reused", fh.read())

    def test_derive_declared(self):
        m = {"manifest_version": 8, "scenarios": [{"id": "T-x-0001", "route": "POST /traverse"}], "relations": []}
        ident = sc.manifest_identity(m)
        m.update(content_sha256=ident["content_sha256"], relations_sha256=ident["relations_sha256"])
        mp, outp = os.path.join(self.tmp, "m.json"), os.path.join(self.tmp, "d.json")
        with open(mp, "w") as fh:
            json.dump(m, fh)
        d = shadow_run.derive_declared(mp, outp)
        with open(shadow_run.DECLARED) as fh:
            v7 = json.load(fh)
        self.assertEqual(d["manifest"]["content_sha256"], ident["content_sha256"])
        for k in ("exclusions", "not_required", "relations_excluded"):
            self.assertEqual(d.get(k), v7.get(k), k)
        self.assertIn("re-pin", d["about"])
        with open(outp) as fh:
            self.assertEqual(json.load(fh), d)

    def test_the_expected_load(self):
        scen = [{"id": "T-x-%04d" % i, "route": "POST /traverse", "mcp": {"tool": "traverse_fetch"}, "repeat": 1}
                for i in range(20)]
        scen.append({"id": "S-stop-0001", "route": "POST /traverse", "direct_not_run": "mcp only",
                     "mcp": {"tool": "traverse_fetch", "flow": {"kind": "stop"}}})
        got = shadow_run.load({"scenarios": scen})
        self.assertEqual((got["scenarios"], got["direct_local"], got["direct_remote"]), (21, 22, 22))
        self.assertEqual(got["mcp"], 20 * 2 + 1)
        self.assertGreater(got["expected_hours"], 0)


class Baseline(unittest.TestCase):
    """Round 8 item 5: --baseline judges the run by the committed freeze package (gate/declared_v8_rt.json and both
    v8 policies): the build a dry freeze, the gate files copied and checked against it, both policies judging with the
    deployed service and library commits."""

    SERVICE = "9abe8517"

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="baseline-")
        self.out = os.path.join(self.tmp, "out")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_dry_run_names_the_package(self):
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rc = shadow_run.main(["--baseline", "--out", self.out, "--dry-run", "--library-commit", PIN,
                                  "--expect-service-commit", self.SERVICE])
        self.assertEqual(rc, 0)
        text = buf.getvalue()
        self.assertTrue(text.startswith("baseline run: %s" % self.out), text)
        lines = dict((ln.split()[0], ln) for ln in text.splitlines()[1:] if ln.startswith("  "))
        self.assertEqual(list(lines), list(shadow_run.STEPS))
        self.assertIn("--stop-scenarios", lines["build"])
        self.assertIn("--freeze", lines["build"])                                  # a dry freeze
        self.assertIn(os.path.join(self.out, "suite_v8.json"), lines["build"])
        self.assertIn("--copy-gate", lines["declared"])
        self.assertIn("declared_v8_rt.json", lines["declared"])
        self.assertIn("policy_v8_rtp_mcp.json policy_v8_rt_mcp.json", lines["policy"])
        self.assertNotIn("gate_policy.py", lines["policy"])
        self.assertIn("--judge", lines["gate"])
        self.assertIn("--expect-service-commit %s" % self.SERVICE, lines["gate"])
        self.assertIn("--expect-library-commit %s" % PIN, lines["gate"])   # the library pin by default
        self.assertIn("--paired", lines["run"])
        self.assertNotIn("--skip-remote", lines["run"])
        self.assertFalse(os.path.exists(self.out))

    def test_the_service_commit_is_required(self):
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            rc = shadow_run.main(["--baseline", "--out", self.out, "--dry-run", "--library-commit", PIN])
        self.assertEqual(rc, 2)
        self.assertIn("--expect-service-commit", err.getvalue())

    def _manifest(self, content=None):
        m = {"manifest_version": 8, "scenarios": [{"id": "T-x-0001"}], "relations": []}
        ident = sc.manifest_identity(m)
        m.update(content_sha256=content or ident["content_sha256"], relations_sha256=ident["relations_sha256"])
        path = os.path.join(self.tmp, "suite_v8.json")
        with open(path, "w") as fh:
            json.dump(m, fh)
        return path, ident

    def test_copy_gate_checks_every_pin(self):
        mpath, ident = self._manifest()
        gate = os.path.join(self.tmp, "gate")
        os.makedirs(gate)
        sc.write_json_atomic(os.path.join(gate, "a.json"), {"manifest": {"file": "manifest/suite_v8.json",
                                                                          "content_sha256": ident["content_sha256"]}})
        sc.write_json_atomic(os.path.join(gate, "b.json"), {"manifest": dict(ident, file="manifest/suite_v8.json")})
        sc.write_json_atomic(os.path.join(gate, "c.json"), {"manifest": dict(ident, content_sha256="0" * 64)})
        dst = os.path.join(self.tmp, "pkg")
        self.assertEqual(shadow_run.copy_gate(mpath, dst, ["a.json", "b.json"], gate), 0)
        self.assertEqual(sorted(os.listdir(dst)), ["a.json", "b.json"])
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            self.assertEqual(shadow_run.copy_gate(mpath, os.path.join(self.tmp, "pkg2"), ["a.json", "c.json"],
                                                  gate), 1)
        self.assertIn("gate/c.json pins content_sha256 0000000000000000", err.getvalue())
        self.assertFalse(os.path.exists(os.path.join(self.tmp, "pkg2")))     # nothing copied

    def test_the_committed_package_is_complete(self):
        for name in (shadow_run.BASELINE_DECLARED,) + shadow_run.BASELINE_POLICIES:
            with open(os.path.join(shadow_run.GATE_DIR, name)) as fh:
                self.assertEqual(json.load(fh)["manifest"]["file"], "manifest/suite_v8.json", name)

    def test_judge_all_runs_both_policies(self):
        run = os.path.join(self.tmp, "run")
        os.makedirs(run)
        calls = []

        def runner(argv, env=None, stdin=None):
            calls.append(argv)
            return 1 if os.path.basename(argv[3]) == "policy_v8_rt_mcp.json" else 0
        pols = [os.path.join(self.out, x) for x in shadow_run.BASELINE_POLICIES]
        with contextlib.redirect_stdout(io.StringIO()):
            rc = shadow_run.judge_all(run, pols, ["--expect-service-commit", self.SERVICE], runner)
        self.assertEqual(rc, 1)                                       # one policy failing fails the baseline
        self.assertEqual([os.path.basename(a[1]) for a in calls], ["gate.py", "gate.py"])
        self.assertEqual([a[a.index("--out") + 1] for a in calls],
                         [os.path.join(run, "gate.json"), os.path.join(run, "gate_policy_v8_rt_mcp.json")])
        for a in calls:
            self.assertEqual(a[-2:], ["--expect-service-commit", self.SERVICE])
        # the pipeline's gate step: shadow_run --judge with both copies
        rec = Recorder()
        shadow_run.main(["--baseline", "--out", self.out, "--from", "gate", "--run-dir", run, "--library-commit",
                         PIN, "--expect-service-commit", self.SERVICE], runner=rec)
        argv = rec.calls[0][0]
        self.assertEqual(argv[2:5], ["--judge", run, pols[0]])
        self.assertEqual(argv[5], pols[1])
        with open(os.path.join(self.out, "baseline.log")) as fh:
            self.assertIn("== gate exit 0", fh.read())


if __name__ == "__main__":
    unittest.main()
