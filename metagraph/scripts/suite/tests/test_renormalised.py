"""A re-evaluation of re-normalised answers says so (GPT review finding 10, VERIFY-A): merge_run records what async's
suite_mcp_renormalise.py left in the arm (its summary, the lines' marks, the raw records keeping answer_recorded),
a plain merge warns that those lines' outcomes are stale, and the re-evaluated report states the re-normalisation
and the changed records apart from the changed outcomes."""

import contextlib
import io
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.dirname(HERE)
sys.path.insert(0, SUITE)
sys.path.insert(0, HERE)

import merge_run  # noqa: E402
import reevaluate  # noqa: E402
import report  # noqa: E402
import suite_common as sc  # noqa: E402
from test_reevaluate import build_run  # noqa: E402


def renormalise(run_dir, changed_key):
    """The traces async's re-normaliser leaves in a copy of the MCP arm: every line marked (one changed), each raw
    record with the recorded answer kept, renormalised.json."""
    arm = os.path.join(run_dir, "arms", "mcp")
    rows, _ = sc.read_jsonl(os.path.join(arm, "executions.jsonl"))
    os.remove(os.path.join(arm, "executions.jsonl"))
    changed = []
    for r in rows:
        key = (r["scenario"], r["repeat"])
        r["renormalised"] = {"changed": key == changed_key, "answer_digest_recorded": r.get("answer_digest")}
        if key == changed_key:
            changed.append({"scenario": key[0], "repeat": key[1], "raw": r.get("raw")})
        sc.append_jsonl(os.path.join(arm, "executions.jsonl"), r)
        raw = os.path.join(arm, r["raw"])
        doc = sc.read_gz_json(raw)
        doc["answer_recorded"] = doc.get("answer")
        sc.write_gz_json(raw, doc)
    sc.write_json_atomic(os.path.join(arm, "renormalised.json"), {
        "source": "/x/arms/mcp", "out": "/x/mcp-renormalised", "at": "2026-10-14T02:07:28.559Z",
        "normaliser_sha256": "db1bd7b6" * 8, "digest": "suite_common.answer_digest of /x",
        "counts": {"lines": len(rows), "renormalised": len(rows), "changed": len(changed),
                   "unchanged": len(rows) - len(changed), "skipped": 0}, "skipped_by_kind": {}, "changed": changed, "skipped": []})
    return rows


class RenormalisedProvenanceTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.run_dir = os.path.join(self.tmp.name, "run")
        build_run(self.run_dir)

    def tearDown(self):
        self.tmp.cleanup()

    def test_traces_counted(self):
        rows = renormalise(self.run_dir, ("T-0001", 1))
        info = merge_run.renormalisation(os.path.join(self.run_dir, "arms", "mcp"))
        self.assertEqual(info["lines"]["rebuilt"], len(rows))
        self.assertEqual(info["lines"]["changed_executions"], ["T-0001 r1"])
        self.assertEqual(info["raw"], {"read": len(rows), "answer_recorded": len(rows)})
        self.assertEqual(info["summary"]["normaliser_sha256"], "db1bd7b6" * 8)
        self.assertEqual(info["summary"]["changed"], ["T-0001 r1"])
        # an arm without any trace: None
        self.assertIsNone(merge_run.renormalisation(os.path.join(self.run_dir, "arms", "direct-local")))

    def test_plain_merge_warns_stale_outcomes(self):
        renormalise(self.run_dir, ("T-0001", 1))
        merged = merge_run.merge(self.run_dir)
        self.assertIn("mcp", merged["renormalisation"])
        self.assertIn("stale until a re-evaluation", " ".join(merged["warnings"]))

    def test_reevaluated_report_states_it(self):
        renormalise(self.run_dir, ("T-0001", 1))
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            mark = reevaluate.reevaluate(self.run_dir)
            merged = merge_run.merge(self.run_dir, reevaluated=True)
            text = report.build_report(self.run_dir, reevaluated=True)
        self.assertEqual(merged["reevaluation"]["renormalised"]["mcp"]["lines"]["changed"], 1)
        self.assertNotIn("stale until a re-evaluation", " ".join(merged["warnings"]))
        self.assertIn("arms/mcp was re-normalised before this re-evaluation", text)
        self.assertIn("normaliser\\_sha256 db1bd7b6db1bd7b6", text)
        self.assertIn("1 changed (T-0001 r1)", text)
        self.assertIn("arms/mcp: its answers were re-normalised before this re-evaluation", text.split("# Re-ev")[0])
        changes = mark.get("changes") or []
        moved = sum(1 for c in changes if c["before"] != c["after"])
        self.assertIn("%d execution record(s) changed: %d changed their outcome, %d only their checks" % (
            len(changes), moved, len(changes) - moved), text)


    def test_wording_without_a_changed_answer(self):
        renormalise(self.run_dir, None)
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            reevaluate.reevaluate(self.run_dir)
            merge_run.merge(self.run_dir, reevaluated=True)
            text = report.build_report(self.run_dir, reevaluated=True)
        sentence = [x for x in text.splitlines() if x.startswith("**arms/mcp was re-normalised")][0]
        self.assertIn("0 changed", sentence)
        self.assertNotIn("other totals", sentence)
        self.assertIn("no answer changed: the re-evaluation of the recorded answers gives the same totals", sentence)
        # the summary records the normaliser's sha256 only: the renormaliser's own revision is said to be missing
        self.assertIn("its own revision (commit, runner and client sha256) is not recorded", sentence)

    def test_wording_with_the_revision_recorded(self):
        renormalise(self.run_dir, ("T-0001", 1))
        path = os.path.join(self.run_dir, "arms", "mcp", "renormalised.json")
        summary = sc.read_json(path)
        summary.update(commit="7a65e913aaaa", runner_sha256="ab" * 32)
        sc.write_json_atomic(path, summary)
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            reevaluate.reevaluate(self.run_dir)
            merge_run.merge(self.run_dir, reevaluated=True)
            text = report.build_report(self.run_dir, reevaluated=True)
        sentence = [x for x in text.splitlines() if x.startswith("**arms/mcp was re-normalised")][0]
        self.assertIn("commit 7a65e913aaaa", sentence)
        self.assertNotIn("is not recorded", sentence)
        self.assertIn("the recorded arm re-evaluated alone gives other totals", sentence)


if __name__ == "__main__":
    unittest.main()
