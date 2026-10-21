"""The results document's Timing 8, "Resources per task" (doc/resources.py, closing.json `resources`): only from
manifest v8 on, the figures not stated when no answer states them, the tables when they do, the container's memory
from host_memory.json, the service's graphlet-op figures, the prose's placeholders, and the overview's v8 branch
naming Timing 8. On test_doc's budget fixture (test_compare_report's synthetic run). Offline; stdlib only."""

import json
import os
import re
import shutil
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.dirname(HERE)
DOC = os.path.join(SUITE, "doc")
for p in (SUITE, HERE, DOC):
    if p not in sys.path:
        sys.path.insert(0, p)

import facts  # noqa: E402
import reader  # noqa: E402
import resources  # noqa: E402
import sections  # noqa: E402
import suite_common as sc  # noqa: E402
import task_resources as tr  # noqa: E402
import units as units_mod  # noqa: E402
from test_task_resources import SNAP_BEFORE, snap_after  # noqa: E402

PROSE_DIR = os.path.join(DOC, "prose")


class ResourcesSectionTest(unittest.TestCase):

    def setUp(self):
        import test_doc
        self.tmp = tempfile.mkdtemp()
        self.d = os.path.join(self.tmp, "run")
        os.makedirs(self.d)
        self.mapping = test_doc.budget_fixture(self.d, False)
        self.closing = sc.read_json(os.path.join(PROSE_DIR, "closing.json"))

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def ctx(self, version=8):
        run = reader.Run(self.d)
        run.manifest = dict(run.manifest, manifest_version=version)
        units = units_mod.resolve(run.manifest, self.mapping)
        return sections.Context(run, units, {}, {}, {}, list(units.chapters), closing=self.closing)

    def add_figures(self):
        path = os.path.join(self.d, "executions.jsonl")
        rows, _ = sc.read_jsonl(path)
        with open(path, "w") as fh:
            for i, r in enumerate(rows):
                if isinstance(r.get("server"), dict):
                    el = r["server"].get("elapsed_ms") or 1.0
                    r["server"]["timing"] = {"elapsed_ms": el, "cpu_ms": el * 0.4,
                                             "page_faults": {"minor": 10, "major": 3 + i}}
                fh.write(json.dumps(r) + "\n")

    def test_only_from_v8(self):
        self.assertNotIn("Resources per task", sections.closing(self.ctx(version=7)))
        text = sections.closing(self.ctx())
        self.assertIn("\\section*{Timing 8\\qd Resources per task}", text)
        self.assertIn("\\label{sec:C8}", text)
        self.assertNotIn("UNRESOLVED", text)

    def test_not_stated(self):
        # the fixture's answers state a memory account (work.memory_bytes on /pattern) but no CPU time or faults
        text = resources.section(self.ctx())
        self.assertIn("not stated", text)
        self.assertIn("states the engine's CPU time (timing.cpu\\_ms): not stated", text)
        self.assertIn("states its page faults", text)
        self.assertIn("the largest memory account per task", text)
        self.assertNotIn("the most CPU time per task", text)
        self.assertIn("holds no \\texttt{host\\_\\allowbreak{}memory", text)
        self.assertNotIn("UNRESOLVED", text)
        # nothing stated at all: one sentence instead of the tables
        path = os.path.join(self.d, "executions.jsonl")
        rows, _ = sc.read_jsonl(path)
        with open(path, "w") as fh:
            for r in rows:
                if isinstance(r.get("server"), dict):
                    r["server"]["work"] = {"steps": 1}
                fh.write(json.dumps(r) + "\n")
        text = resources.section(self.ctx())
        self.assertIn("No task execution of this run states the engine's CPU time, page faults or memory account",
                      text)
        self.assertNotIn("per chapter, access path and host-cache", text)

    def test_figures_and_host_memory(self):
        self.add_figures()
        with open(os.path.join(self.d, tr.HOST_MEMORY_FILE), "w") as fh:
            json.dump({"before": SNAP_BEFORE, "after": snap_after(**{
                "events.max": 83150, "stat.pgsteal_direct": SNAP_BEFORE["stat.pgsteal_direct"] + 262144})}, fh)
        text = resources.section(self.ctx())
        self.assertIn("the engine's CPU time per task (timing.cpu\\_ms)", text)
        self.assertIn("the 20 tests with the most CPU time per task", text)
        self.assertIn("the 20 tests with the most major page faults per task", text)
        self.assertIn("server block", text)
        self.assertIn("limit hits (events max) & ~ & ~ & 8", text)
        self.assertIn("262,144 pages (1.0 GiB)", text)
        self.assertIn("the index (648~GiB on disk) is about 5 times the container's memory limit of 128~GiB", text)
        self.assertIn("(2026-10-21T06:02:40Z)", text)
        self.assertNotIn("UNRESOLVED", text)
        self.assertNotIn("No task execution of this run states", text)

    def test_graphlet_ops(self):
        path = os.path.join(self.d, "arms", "mcp", "graphlet_ops.jsonl")
        sc.append_jsonl(path, {"way": "W3", "tool": "graphlet_summary", "answer": {"summary": {}}})
        self.assertIn("state no cost of their own in this run: not stated", resources.section(self.ctx()))
        sc.append_jsonl(path, {"way": "W3", "tool": "graphlet_walk", "answer": {"usage": {"cpu_ms": 3.0}}})
        self.assertIn("state their cost on 1 of 2 calls", resources.section(self.ctx()))

    def test_prose_placeholders(self):
        ctx = self.ctx()
        s = tr.of_run(self.d, ctx.run.rows, ctx.run.manifest)
        vals = resources.values(ctx, s, tr.host_memory_of_run(self.d), None)
        text = json.dumps(self.closing["resources"])
        used = set(re.findall(r"\{([a-z_0-9]+)\}", text))
        flags = set(n.strip() for n in re.findall(r"\[\[\s*([a-z_0-9]+)\s*\?", text))
        self.assertEqual(used - set(vals), set())
        self.assertEqual(flags - set(vals), set())
        for f in flags:
            self.assertIsInstance(vals[f], bool, f)


class OverviewBranchTest(unittest.TestCase):

    def test_timing_eight_named_from_v8(self):
        ov = sc.read_json(os.path.join(PROSE_DIR, "overview.json"))
        texts = [t for t in json.dumps(ov).split('", "') if "Timing 1 to" in t]
        self.assertEqual(len(texts), 2)
        for t in texts:
            v7 = facts.fill(t, {"manifest_version": 7})
            v8 = facts.fill(t, {"manifest_version": 8})
            self.assertIn("(Timing 1 to 7)", v7)
            self.assertNotIn("Timing 8", v7)
            self.assertIn("(Timing 1 to 8)", v8)


if __name__ == "__main__":
    unittest.main()
