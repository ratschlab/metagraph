"""Offline tests of the bench run's load guard (bench/load_guard.py) on synthetic samples and a synthetic run
directory: clean, loaded and unguarded blocks, the coverage rule, the clock offset and the missing Logan container."""

import contextlib
import io
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.dirname(HERE)
sys.path.insert(0, SUITE)
sys.path.insert(0, os.path.join(SUITE, "bench"))

import bench_lib as bl  # noqa: E402
import load_guard as lg  # noqa: E402
import suite_common as sc  # noqa: E402

T0 = 1791960000.0  # 2026-10-14


def samples(t_from, t_to, step=2.0, sra=3.0, load1=2.0, cpu=50.0, **overrides):
    """Samples every `step` s on mex's clock; overrides {offset_s: {field: value}} change single ones."""
    out, t = [], t_from
    while t <= t_to:
        s = {"t_utc": sc.iso_utc(t), "sra_all_cpu_pct": sra, "host_load1": load1, "cpu_pct": cpu, "phase": "bench"}
        out.append(s)
        t += step
    return out


class JudgeSamples(unittest.TestCase):
    def test_clean_loaded_unguarded(self):
        self.assertEqual(lg.judge_samples(samples(T0, T0 + 20), None)["verdict"], "clean")
        s = samples(T0, T0 + 20)
        s[3]["sra_all_cpu_pct"] = 7700.0
        res = lg.judge_samples(s, None)
        self.assertEqual((res["verdict"], res["breach_count"]), ("loaded", 1))
        self.assertEqual(res["breaches"][0]["over"], ["sra_all_cpu_pct"])
        s = samples(T0, T0 + 20, load1=30.0)
        self.assertEqual(lg.judge_samples(s, None)["verdict"], "loaded")
        self.assertEqual(lg.judge_samples(s, {"host_load1_max": 40.0})["verdict"], "clean")
        self.assertEqual(lg.judge_samples([], None)["verdict"], "unguarded")

    def test_missing_logan_reading(self):
        s = samples(T0, T0 + 10, sra=None)
        self.assertEqual(lg.judge_samples(s, None)["verdict"], "unguarded")
        # declared stopped: a sample whose staging container was read counts 0 %; one docker did not answer stays unread
        self.assertEqual(lg.judge_samples(s, {"allow_missing_sra": True})["verdict"], "clean")
        s[0]["cpu_pct"] = None
        self.assertEqual(lg.judge_samples(s, {"allow_missing_sra": True})["verdict"], "unguarded")


class JudgeBlocks(unittest.TestCase):
    def test_windows_coverage_and_offset(self):
        s = samples(T0 - 60, T0 + 600)
        # block 1 ran 100-200 s, block 2 300-400 s; a spike at 350 s loads block 2 only
        for x in s:
            if abs(sc.parse_utc(x["t_utc"]) - (T0 + 350)) < 1:
                x["sra_all_cpu_pct"] = 5000.0
        res = lg.judge_blocks({1: [[T0 + 100, T0 + 200]], 2: [[T0 + 300, T0 + 400]]}, s, None)
        self.assertEqual(res["blocks"]["1"]["verdict"], "clean")
        self.assertEqual(res["blocks"]["2"]["verdict"], "loaded")
        self.assertEqual((res["verdict"], res["loaded"]), ("loaded", ["2"]))
        # the pad (30 s) counts the load just before a block
        res = lg.judge_blocks({2: [[T0 + 375, T0 + 400]]}, s, None)
        self.assertEqual(res["verdict"], "loaded")
        res = lg.judge_blocks({2: [[T0 + 385, T0 + 400]]}, s, None)
        self.assertEqual(res["verdict"], "clean")
        # mex's clock 100 s ahead: the same samples now cover 250-300 s of the orchestrator's clock as 350-400 s
        res = lg.judge_blocks({2: [[T0 + 250, T0 + 260]]}, s, None, offset_s=100.0)
        self.assertEqual(res["verdict"], "loaded")

    def test_gap_and_no_execution(self):
        s = samples(T0, T0 + 100) + samples(T0 + 160, T0 + 300)
        res = lg.judge_blocks({1: [[T0 + 50, T0 + 200]]}, s, None)
        self.assertEqual(res["blocks"]["1"]["verdict"], "unguarded")
        self.assertTrue(res["blocks"]["1"]["gaps"])
        res = lg.judge_blocks({1: [[T0 + 50, T0 + 200]]}, s, {"max_gap_s": 120.0})
        self.assertEqual(res["verdict"], "clean")
        self.assertEqual(lg.judge_blocks({}, s, None)["verdict"], "unguarded")
        # samples that stop before the block ends leave its end uncovered
        res = lg.judge_blocks({1: [[T0 + 250, T0 + 400]]}, s, None)
        self.assertEqual(res["verdict"], "unguarded")

    def test_resumed_block_has_two_windows(self):
        w = lg.merge_windows([[T0, T0 + 10], [T0 + 20, T0 + 30], [T0 + 1000, T0 + 1010]], 30.0)
        self.assertEqual(len(w), 2)
        self.assertEqual(w[0], [T0 - 30, T0 + 60])


class RunDirectory(unittest.TestCase):
    def test_judge_run_reads_schedule_executions_samples(self):
        with tempfile.TemporaryDirectory() as d:
            steps = bl.schedule(["Q1", "Q2"], seed=1, cycles=2)
            sc.write_json_atomic(os.path.join(d, "schedule.json"), {"seed": 1, "steps": steps})
            offset = 5.0
            sc.write_json_atomic(os.path.join(d, "bench.json"), {
                "guard": {"host_load1_max": 8.0}, "sessions": [{"clock_offset_start_s": offset,
                                                                "clock_offset_end_s": offset}]})
            t = T0
            for s in steps:
                arm = s["arm"]
                mex = offset if arm == "direct-local" else 0.0
                os.makedirs(os.path.join(d, "arms", arm), exist_ok=True)
                sc.append_jsonl(os.path.join(d, "arms", arm, "executions.jsonl"), {
                    "scenario": s["bench_id"], "arm": arm, "repeat": 0, "started_utc": sc.iso_utc(t + mex),
                    "ended_utc": sc.iso_utc(t + 1 + mex), "wall_ms": 1000.0})
                t += 3
            os.makedirs(os.path.join(d, "samples"))
            rows = samples(T0 - 40 + offset, t + 40 + offset, load1=4.0)
            for r in rows:
                sc.append_jsonl(os.path.join(d, "samples", "bench-x.jsonl"), r)
            res = lg.judge_run(d)
            self.assertEqual(sorted(res["blocks"]), ["0", "1", "2"])
            self.assertEqual(res["verdict"], "clean")
            self.assertEqual(res["limits"]["host_load1_max"], 8.0)  # the run's own threshold
            self.assertEqual(res["clock_offset_s"], offset)
            self.assertEqual(lg.judge_run(d, {"host_load1_max": 3.0})["verdict"], "loaded")
            with contextlib.redirect_stdout(io.StringIO()) as said:
                self.assertEqual(lg.main([d]), 0)
                self.assertEqual(lg.main([d, "--load1-max", "3"]), 1)
            self.assertIn("run: clean", said.getvalue())
            self.assertIn("run: loaded", said.getvalue())


if __name__ == "__main__":
    unittest.main()
