"""The per-task resource figures (task_resources.py) on synthetic answers: each figure from the server block, else
from the raw answer, per entry where the answer states it per entry only; null where nothing states it (an older
engine, a refusal); the summaries per area, arm and host cache, per test with the cold execution apart, the heaviest
tests; host_memory.json's levels and deltas; the service's graphlet-op figures; compare.json's block and the report's
section, which stays away from a v1-v7 report. Offline; stdlib only."""

import copy
import json
import os
import shutil
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.dirname(HERE)
for p in (SUITE, HERE):
    if p not in sys.path:
        sys.path.insert(0, p)

import compare as cmp  # noqa: E402
import merge_run  # noqa: E402
import report  # noqa: E402
import suite_common as sc  # noqa: E402
import task_resources as tr  # noqa: E402
import test_compare_report as tcr  # noqa: E402

SNAP_BEFORE = {"label": "before", "utc": "2026-10-21T06:02:40Z", "container": "metagraph_refseq33m_experimental",
               "page_bytes": 4096, "current": 108165853184, "peak": 137438957568, "max": 137438953472,
               "stat.anon": 9465450496, "stat.file": 96841510912, "stat.pgmajfault": 233625900,
               "stat.pgscan_direct": 246266849, "stat.pgscan_kswapd": 125379262, "stat.pgsteal_direct": 220276454,
               "stat.pgsteal_kswapd": 106915779, "stat.workingset_refault_file": 63140236, "events.max": 83142,
               "events.oom": 0, "events.oom_kill": 0, "events.high": 0, "host_available_kib": 731838912}


def snap_after(**changes):
    s = dict(SNAP_BEFORE, label="after", utc="2026-10-21T12:35:00Z")
    s.update(changes)
    return s


def traverse_body(cpu=12.5, el=20.0, major=7, minor=300, mem=4096, seed_mem=(1024, 2048)):
    """A /traverse answer as the per-task engine writes it: top-level timing and work, per seed the same blocks."""
    return {"timing": {"elapsed_ms": el, "serialize_ms": 0.1, "cpu_ms": cpu,
                       "page_faults": {"minor": minor, "major": major}},
            "work": {"memory_bytes": mem},
            "results": [{"timing": {"elapsed_ms": el / 2, "cpu_ms": cpu / 2, "page_faults": {"minor": 1, "major": 1}},
                         "work": {"memory_bytes": m}} for m in seed_mem]}


def pattern_body_entries_only():
    """A /pattern answer the service passes through without its top-level timing (load_ms only)."""
    return {"timing": {"load_ms": 0.0},
            "patterns": [{"id": "a", "timing": {"elapsed_ms": 4.0, "cpu_ms": 3.0, "page_faults": {"minor": 2, "major": 5}},
                          "work": {"steps": 1, "memory_bytes": 700}},
                         {"id": "b", "timing": {"elapsed_ms": 6.0, "cpu_ms": 1.0, "page_faults": {"minor": 3, "major": 0}},
                          "work": {"steps": 1, "memory_bytes": 900}}]}


def exec_line(sid, arm, repeat, server, route="POST /traverse", host="warm", raw=None, http=200, area=None):
    return {"scenario": sid, "arm": arm, "repeat": repeat, "route": route, "host_cache": host, "http": http,
            "outcome": "pass", "server": server, "raw": raw, "area": area or sid.split("-")[0]}


class FiguresTest(unittest.TestCase):

    def test_server_block_top_level(self):
        srv = sc.extract_server(traverse_body())
        f = tr.figures_of_server(srv)
        self.assertEqual((f["cpu_ms"], f["elapsed_ms"], f["minor_faults"], f["major_faults"]), (12.5, 20.0, 300, 7))
        # the request's account is the maximum over the request and its seeds (extract_server's max)
        self.assertEqual(f["memory_bytes"], 4096)
        self.assertEqual(f["sources"], {"cpu_ms": "server", "page_faults": "server", "memory_bytes": "server"})

    def test_server_block_entries_summed(self):
        srv = sc.extract_server(pattern_body_entries_only())
        # extract_server keeps the top-level timing whole (load_ms only): the entries come from the body ...
        self.assertIsNone(tr.figures_of_server(srv)["cpu_ms"])
        # ... while a runner's server block with the entries (entries / by_pattern) sums them
        for timing in ({"entries": [{"elapsed_ms": 4.0, "cpu_ms": 3.0, "page_faults": {"minor": 2, "major": 5}},
                                    {"elapsed_ms": 6.0, "cpu_ms": 1.0, "page_faults": {"minor": 3, "major": 0}}]},
                       {"by_pattern": {"a": {"elapsed_ms": 4.0, "cpu_ms": 3.0, "page_faults": {"minor": 2, "major": 5}},
                                       "b": {"elapsed_ms": 6.0, "cpu_ms": 1.0,
                                             "page_faults": {"minor": 3, "major": 0}}}}):
            f = tr.figures_of_server({"elapsed_ms": 10.0, "timing": timing,
                                      "work": {"by_pattern": {"a": {"memory_bytes": 700}, "b": {"memory_bytes": 900}}}})
            self.assertEqual((f["cpu_ms"], f["elapsed_ms"], f["minor_faults"], f["major_faults"]), (4.0, 10.0, 5, 5))
            self.assertEqual(f["memory_bytes"], 900)
            self.assertEqual(f["sources"]["cpu_ms"], "server_entries")
            self.assertEqual(f["sources"]["page_faults"], "server_entries")
        # a flattened key (extract_server's dotted names) counts too
        self.assertEqual(tr.figures_of_server({"work": {"by_pattern.a.memory_bytes": 5, "steps": 9}})["memory_bytes"], 5)

    def test_absent_figures_are_null_never_zero(self):
        # an older engine: timing without cpu_ms, work without memory_bytes; a macOS engine: no page_faults
        f = tr.figures_of_server({"elapsed_ms": 3.0, "timing": {"elapsed_ms": 3.0, "cpu_ms": 2.0}, "work": {"steps": 4}})
        self.assertEqual(f["cpu_ms"], 2.0)
        self.assertIsNone(f["major_faults"])
        self.assertIsNone(f["memory_bytes"])
        self.assertEqual(f["sources"], {"cpu_ms": "server", "page_faults": None, "memory_bytes": None})
        for server in (None, {}, {"elapsed_ms": 1.0, "timing": None, "work": None}):
            f = tr.figures_of_server(server)
            self.assertEqual([f[k] for k in ("cpu_ms", "elapsed_ms", "minor_faults", "major_faults", "memory_bytes")],
                             [None] * 5)

    def test_body_figures(self):
        f = tr.figures_of_body(traverse_body())
        self.assertEqual((f["cpu_ms"], f["major_faults"], f["memory_bytes"]), (12.5, 7, 4096))
        self.assertEqual(f["sources"], {"cpu_ms": "raw", "page_faults": "raw", "memory_bytes": "raw"})
        # a body without the request's work: the seeds' maximum
        body = traverse_body()
        body.pop("work")
        self.assertEqual(tr.figures_of_body(body)["memory_bytes"], 2048)
        # /pattern without its top-level timing: the entries summed
        f = tr.figures_of_body(pattern_body_entries_only())
        self.assertEqual((f["cpu_ms"], f["elapsed_ms"], f["major_faults"], f["memory_bytes"]), (4.0, 10.0, 5, 900))
        self.assertEqual(f["sources"]["cpu_ms"], "raw_entries")
        # /resolve: top level only
        f = tr.figures_of_body({"timing": {"elapsed_ms": 9.0, "cpu_ms": 8.0}, "work": {"memory_bytes": 77}})
        self.assertEqual((f["cpu_ms"], f["memory_bytes"], f["major_faults"]), (8.0, 77, None))

    def test_execution_falls_back_to_the_raw_answer(self):
        """The MCP shape (async's runner a52bd3d6): server.timing is the body's, server.work the service's usage
        without memory_bytes, the engine's account in the raw answer only."""
        d = tempfile.mkdtemp()
        try:
            os.makedirs(os.path.join(d, "raw"))
            sc.write_gz_json(os.path.join(d, "raw", "T-1.mcp.0.json.gz"),
                             {"answer": {"http": 200, "body": traverse_body(mem=5000)}})
            server = {"elapsed_ms": 20.0, "timing": traverse_body()["timing"],
                      "work": {"attempts_total": 1, "work_units_known": 10}}
            line = exec_line("T-1", "mcp", 0, server, raw="raw/T-1.mcp.0.json.gz", host="cold")
            x = tr.execution_figures(line, d)
            self.assertEqual((x["cpu_ms"], x["major_faults"], x["memory_bytes"]), (12.5, 7, 5000))
            self.assertEqual(x["sources"], {"cpu_ms": "server", "page_faults": "server", "memory_bytes": "raw"})
            self.assertEqual(x["cpu_ratio"], 12.5 / 20.0)
            self.assertEqual(x["host"], "cold")
            # without reading the raw records: the memory account not stated
            self.assertIsNone(tr.execution_figures(line, d, read_raw=False)["memory_bytes"])
            # a refusal: neither the server block nor the body states anything
            sc.write_gz_json(os.path.join(d, "raw", "T-2.direct-local.0.json.gz"),
                             {"answer": {"http": 400, "body": {"error": {"code": "invalid_request"}}}})
            x = tr.execution_figures(exec_line("T-2", "direct-local", 0, None, raw="raw/T-2.direct-local.0.json.gz",
                                               http=400), d)
            self.assertEqual([x[k] for k in ("cpu_ms", "major_faults", "memory_bytes", "cpu_ratio")], [None] * 4)
            self.assertEqual(x["sources"], {"cpu_ms": None, "page_faults": None, "memory_bytes": None})
            # a line without a route takes its scenario's; another route is never read from the raw record
            line = exec_line("T-1", "mcp", 0, None, route=None, raw="raw/T-1.mcp.0.json.gz")
            self.assertEqual(tr.execution_figures(line, d, scenario={"route": "POST /traverse", "area": "T"})
                             ["memory_bytes"], 5000)
            line["route"] = "GET /traverse/capabilities"
            self.assertIsNone(tr.execution_figures(line, d)["memory_bytes"])
        finally:
            shutil.rmtree(d)


def _srv(cpu, el, major, mem):
    return {"elapsed_ms": el, "timing": {"elapsed_ms": el, "cpu_ms": cpu, "page_faults": {"minor": 1, "major": major}},
            "work": {"memory_bytes": mem}}


class SummaryTest(unittest.TestCase):

    def lines(self):
        out = []
        for i, (sid, area) in enumerate((("T-1", "T"), ("T-2", "T"), ("R-1", "R"))):
            base = 10.0 * (i + 1)
            out.append(exec_line(sid, "mcp", 0, _srv(base * 0.2, base * 2, 500 + i, 1000 * (i + 1)), host="cold",
                                 route="POST /resolve" if area == "R" else "POST /traverse"))
            out.append(exec_line(sid, "mcp", 1, _srv(base * 0.9, base, i, 1000 * (i + 1)),
                                 route="POST /resolve" if area == "R" else "POST /traverse"))
            out.append(exec_line(sid, "direct-local", 0, _srv(base, base * 1.1, 0, 1000 * (i + 1)),
                                 route="POST /resolve" if area == "R" else "POST /traverse"))
        out.append(exec_line("T-3", "direct-local", 0, None, http=400))     # a refusal: counted, states nothing
        out.append(exec_line("C-1", "direct-local", 0, None, route="GET /capabilities"))   # not a task route
        return out

    def test_by_area_arm_and_host(self):
        execs = tr.run_figures(self.lines())
        s = tr.summarise_run(execs, order=["T-1", "T-2", "R-1", "T-3"])
        self.assertEqual(s["format"], tr.FORMAT)
        self.assertEqual(s["executions"], {"mcp": 6, "direct-local": 4})
        self.assertEqual(s["stated"]["direct-local"], {"cpu_ms": 3, "major_faults": 3, "memory_bytes": 3})
        t = s["by_area"]["T"]["mcp"]
        self.assertEqual(sorted(t), ["cold", "warm"])
        self.assertEqual(t["cold"]["cpu_ms"]["n"], 2)
        self.assertEqual(t["cold"]["major_faults"]["max"], 501)
        self.assertEqual(t["warm"]["cpu_ratio"]["p50"], 0.9)
        self.assertEqual(s["by_area"]["all"]["mcp"]["cold"]["executions"], 3)
        self.assertEqual(s["by_area"]["T"]["direct-local"]["warm"]["executions"], 3)   # the refusal counted ...
        self.assertEqual(s["by_area"]["T"]["direct-local"]["warm"]["cpu_ms"]["n"], 2)  # ... stating nothing
        self.assertNotIn("C", s["by_area"])
        self.assertEqual(s["sources"]["direct-local"]["cpu_ms"], {"server": 3, "none": 1})

    def test_per_scenario_cold_apart(self):
        ps = tr.per_scenario(tr.run_figures(self.lines()))
        m = ps["T-2"]["mcp"]
        self.assertEqual([e["repeat"] for e in m["executions"]], [0, 1])
        self.assertEqual(m["cold"]["cpu_ms"]["max"], 4.0)
        self.assertEqual(m["warm"]["cpu_ms"]["max"], 18.0)
        self.assertNotIn("cold", ps["T-2"]["direct-local"])

    def test_top_tests(self):
        execs = tr.run_figures(self.lines())
        top = tr.top_tests(execs, "cpu_ms", n=2, order=["T-1", "T-2", "R-1"])
        self.assertEqual([e["scenario"] for e in top], ["R-1", "T-2"])
        self.assertEqual((top[0]["value"], top[0]["arm"], top[0]["repeat"], top[0]["n"]), (30.0, "direct-local", 0, 3))
        self.assertEqual(top[0]["median"], 27.0)
        # the host-cold execution apart from the warm repeats: R-1's mcp r0 (cold) 6 ms, warm 27 and 30 ms
        self.assertEqual((top[0]["cold"], top[0]["warm_median"]), (6.0, 28.5))
        # ties: the manifest's order
        mem = tr.top_tests(execs, "memory_bytes", n=3, order=["T-1", "T-2", "R-1"])
        self.assertEqual([e["scenario"] for e in mem], ["R-1", "T-2", "T-1"])
        faults = tr.top_tests(execs, "major_faults", n=1)
        self.assertEqual((faults[0]["scenario"], faults[0]["value"], faults[0]["host"]), ("R-1", 502, "cold"))

    def test_nothing_stated(self):
        lines = [exec_line("T-1", "mcp", 0, {"elapsed_ms": 3.0, "timing": {"elapsed_ms": 3.0}, "work": None})]
        s = tr.summarise_run(tr.run_figures(lines, read_raw=False))
        self.assertFalse(tr.any_stated(s))
        self.assertEqual(s["top"], {"cpu_ms": [], "memory_bytes": [], "major_faults": []})
        self.assertFalse(tr.engine_states_figures(lines))
        # memory_bytes alone (an older /pattern answer) is stated, but is no per-task engine's figure
        lines[0]["server"]["work"] = {"memory_bytes": 5}
        self.assertTrue(tr.any_stated(tr.summarise_run(tr.run_figures(lines, read_raw=False))))
        self.assertFalse(tr.engine_states_figures(lines))
        lines[0]["server"]["timing"]["cpu_ms"] = 1.0
        self.assertTrue(tr.engine_states_figures(lines))


class HostMemoryTest(unittest.TestCase):

    def setUp(self):
        self.d = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.d)

    def write(self, doc):
        with open(os.path.join(self.d, tr.HOST_MEMORY_FILE), "w") as fh:
            fh.write(doc if isinstance(doc, str) else json.dumps(doc))

    def test_absent(self):
        self.assertEqual(tr.read_host_memory(self.d), (None, ""))
        self.assertEqual(tr.host_memory_of_run(self.d), {"figures": None, "problem": None})

    def test_malformed(self):
        self.write("{not json")
        doc, why = tr.read_host_memory(self.d)
        self.assertIsNone(doc)
        self.assertIn("cannot be read", why)
        self.write({"before": SNAP_BEFORE})
        self.assertIn("is not", tr.read_host_memory(self.d)[1])

    def test_levels_and_deltas(self):
        self.write({"before": SNAP_BEFORE, "after": snap_after(**{
            "current": 120 << 30, "stat.pgsteal_direct": SNAP_BEFORE["stat.pgsteal_direct"] + 262144,
            "stat.pgsteal_kswapd": SNAP_BEFORE["stat.pgsteal_kswapd"] + 1024,
            "stat.workingset_refault_file": SNAP_BEFORE["stat.workingset_refault_file"] + 524288,
            "events.max": 83150, "stat.pgmajfault": SNAP_BEFORE["stat.pgmajfault"] + 99})})
        f = tr.host_memory_of_run(self.d)["figures"]
        self.assertEqual(f["limit_bytes"], 137438953472)
        self.assertEqual(f["levels"]["current"], {"before": 108165853184, "after": 120 << 30})
        self.assertEqual(f["levels"]["peak"]["before"], 137438957568)
        self.assertEqual(f["deltas"]["events.max"], 8)
        self.assertEqual(f["deltas"]["events.oom_kill"], 0)
        self.assertEqual(f["deltas"]["stat.pgmajfault"], 99)
        self.assertEqual(f["pages_gib"]["stat.pgsteal_direct"], 1.0)        # 262144 pages of 4 KiB
        self.assertEqual(f["pages_gib"]["stat.workingset_refault_file"], 2.0)
        self.assertNotIn("stat.pgmajfault", f["pages_gib"])                # a count of faults, not of pages
        self.assertFalse(f["restarted"])
        self.assertEqual(f["notes"], [])
        self.assertEqual((f["before"]["utc"], f["after"]["utc"]), ("2026-10-21T06:02:40Z", "2026-10-21T12:35:00Z"))

    def test_restart_no_limit_other_container(self):
        self.write({"before": dict(SNAP_BEFORE, max="max"),
                    "after": snap_after(**{"max": "max", "events.max": 3, "container": "other"})})
        f = tr.host_memory_of_run(self.d)["figures"]
        self.assertIsNone(f["limit_bytes"])
        self.assertIsNone(f["deltas"]["events.max"])
        self.assertTrue(f["restarted"])
        self.assertEqual(len(f["notes"]), 2)

    def test_page_size_and_missing_counter(self):
        before = dict(SNAP_BEFORE, page_bytes=16384)
        before.pop("stat.pgscan_direct")
        self.write({"before": before, "after": snap_after(**{
            "page_bytes": 16384, "stat.pgsteal_kswapd": SNAP_BEFORE["stat.pgsteal_kswapd"] + 65536})})
        f = tr.host_memory_of_run(self.d)["figures"]
        self.assertEqual(f["pages_gib"]["stat.pgsteal_kswapd"], 1.0)
        self.assertIsNone(f["deltas"]["stat.pgscan_direct"])


class GraphletOpsTest(unittest.TestCase):

    def test_stated_or_not(self):
        d = tempfile.mkdtemp()
        try:
            self.assertIsNone(tr.graphlet_service_figures(d))
            path = os.path.join(d, "arms", "mcp", "graphlet_ops.jsonl")
            os.makedirs(os.path.dirname(path))
            for r in ({"way": "W3", "tool": "graphlet_summary", "answer": {"summary": {}}},
                      {"way": "W4", "tool": "graphlet_summary", "answer": {"cpu_ms": 99.0}}):
                sc.append_jsonl(path, r)
            g = tr.graphlet_service_figures(d)
            self.assertEqual((g["ops"], g["stating"], g["cpu_ms"]["n"]), (1, 0, 0))
            sc.append_jsonl(path, {"way": "W3", "tool": "graphlet_walk", "answer": {"usage": {"cpu_ms": 4.0}}})
            sc.append_jsonl(path, {"way": "W3", "tool": "graphlet_walk", "answer": {"cpu_ms": 2.0}})
            g = tr.graphlet_service_figures(d)
            self.assertEqual((g["ops"], g["stating"], g["cpu_ms"]["max"]), (3, 2, 4.0))
        finally:
            shutil.rmtree(d)


class CompareReportTest(unittest.TestCase):
    """compare.json's block and the report's section on test_compare_report's synthetic run."""

    def setUp(self):
        self.d = tempfile.mkdtemp()
        tcr.build_run(self.d)
        merge_run.merge(self.d)

    def tearDown(self):
        shutil.rmtree(self.d)

    def rows(self):
        return sc.read_jsonl(os.path.join(self.d, "executions.jsonl"))[0]

    def write_rows(self, rows):
        with open(os.path.join(self.d, "executions.jsonl"), "w") as fh:
            for r in rows:
                fh.write(json.dumps(r) + "\n")

    def with_figures(self):
        rows = self.rows()
        for r in rows:
            if isinstance(r.get("server"), dict):
                el = r["server"].get("elapsed_ms") or 1.0
                r["server"]["timing"] = {"elapsed_ms": el, "cpu_ms": el / 2.0, "page_faults": {"minor": 3, "major": 2}}
        self.write_rows(rows)
        return rows

    def test_compare_block_and_unchanged_verdicts(self):
        table = (tcr.CASES, {"cases": len(tcr.CASES)})
        plain = cmp.compare(self.d, table)
        self.assertIn("task_resources", plain["summary"])
        self.assertFalse(tr.engine_states_figures(self.rows()))
        self.with_figures()
        figured = cmp.compare(self.d, table)
        tblock = figured["summary"]["task_resources"]
        self.assertEqual(tblock["format"], tr.FORMAT)
        self.assertTrue(tblock["stated"]["direct-local"]["cpu_ms"])
        view = figured["scenarios"]["P-001"]["arms"]["direct-local"]["task_resources"]
        self.assertEqual([e["repeat"] for e in view["executions"]], [0, 1])
        self.assertEqual(view["executions"][0]["sources"]["cpu_ms"], "server")
        # the figures are report-only: every verdict of the comparison is as without them
        for key in ("direct_agreement", "mcp_agreement", "per_area", "relations", "inconclusive_by_reason"):
            self.assertEqual(plain["summary"][key], figured["summary"][key], key)
        for sid, entry in plain["scenarios"].items():
            other = figured["scenarios"][sid]
            for k in ("direct_agreement", "mcp_agreement"):
                self.assertEqual(entry.get(k), other.get(k), (sid, k))
            for arm, v in entry["arms"].items():
                a = dict(v, task_resources=None)
                b = dict(other["arms"][arm], task_resources=None)
                self.assertEqual(a, b, (sid, arm))

    def test_report_section(self):
        manifest = sc.read_json(os.path.join(self.d, "manifest.json"))
        rows = self.rows()
        hm = tr.host_memory_of_run(self.d)
        # a v1-v7 run whose answers state no per-task figure: no section, the report as it was
        self.assertFalse(report.task_resources_wanted(manifest, rows, hm))
        self.assertNotIn("# Resources per task", report.build_report(self.d))
        # manifest v8: the section, the figures not stated
        v8 = dict(manifest, manifest_version=8)
        self.assertTrue(report.task_resources_wanted(v8, rows, hm))
        text = report.task_resources_section(self.d, rows, v8)
        self.assertIn("# Resources per task", text)
        self.assertIn("No `host_memory.json` in the run directory", text)
        # the figures on record (the engine states them) and a host_memory.json
        rows = self.with_figures()
        with open(os.path.join(self.d, tr.HOST_MEMORY_FILE), "w") as fh:
            json.dump({"before": SNAP_BEFORE, "after": snap_after(**{"events.max": 83150})}, fh)
        self.assertTrue(report.task_resources_wanted(manifest, rows, tr.host_memory_of_run(self.d)))
        os.remove(os.path.join(self.d, "compare.json")) if os.path.exists(os.path.join(self.d, "compare.json")) else None
        full = report.build_report(self.d)
        self.assertIn("# Resources per task", full)
        self.assertIn("## CPU time per area, arm and host cache", full)
        self.assertIn("## The 20 tests with the most CPU time", full)
        self.assertIn("| events.max | count |  |  | 8 |", full)
        self.assertIn("| direct-local | ", full)
        self.assertIn("server ", full)


if __name__ == "__main__":
    unittest.main()
