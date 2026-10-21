"""The closing section's "Resources per task" (Timing 8, manifest v8 on): what each task cost the engine, read from the
engine's own per-request figures (task_resources.py: timing.cpu_ms, timing.page_faults, work.memory_bytes) and the
engine container's cgroup memory over the run (host_memory.json); the prose is closing.json's `resources`.

Report-only: none of these figures is read by a check, a digest, an identity or the gate.
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
SUITE = os.path.dirname(HERE)
if SUITE not in sys.path:
    sys.path.insert(0, SUITE)

import cells  # noqa: E402
import latex  # noqa: E402
import task_resources as tr  # noqa: E402
from cells import ARMS, fmt_t, fmt_n  # noqa: E402
from latex import esc, paragraphs  # noqa: E402

LABEL = "sec:C8"
TITLE = "Resources per task"
# the indexes' size on disk, by the manifest's index database (not in any run record): refseq33m-experimental's
# files on mex, measured 2026-10-21 (ls -laL: graph 362.0 GB, annotation 332.9 GB, .seqs 0.64 GB, edge mask
# 0.97 GB)
INDEX_ON_DISK_BYTES = {"refseq33m-experimental": 648 * (1 << 30)}
HOST_WORDS = {"cold": "cold", "warm": "warm", "warm_after_interrupted": "warm after a cut", "aged": "aged",
              "not_reached": "not reached", "unknown": "unknown"}
SOURCE_WORDS = (("server", "server block"), ("server_entries", "server block (entries summed)"),
                ("raw", "raw answer"), ("raw_entries", "raw answer (entries summed)"), ("none", "not stated"))
HOST_FIGURE_WORDS = {
    "max": "memory limit (memory.max)", "current": "usage (memory.current)",
    "peak": "peak since the container started (memory.peak)", "stat.anon": "anonymous memory (stat anon)",
    "stat.file": "page cache (stat file)", "host_available_kib": "host's available memory (MemAvailable)",
    "events.max": "limit hits (events max)", "events.high": "high-boundary hits (events high)",
    "events.oom": "out-of-memory events (events oom)", "events.oom_kill": "OOM kills (events oom\\_kill)",
    "stat.pgsteal_direct": "evicted under the container's limit (pgsteal\\_direct)",
    "stat.pgsteal_kswapd": "evicted by host reclaim (pgsteal\\_kswapd)",
    "stat.pgscan_direct": "scanned, direct reclaim (pgscan\\_direct)",
    "stat.pgscan_kswapd": "scanned, host reclaim (pgscan\\_kswapd)",
    "stat.workingset_refault_file": "refaulted file pages (workingset\\_refault\\_file)",
    "stat.pgmajfault": "major faults, whole container (pgmajfault)"}


def fmt_bytes(v):
    """A byte count with its unit (B, KiB, MiB, GiB); `–` for none."""
    x = cells._num(v)
    if x is None:
        return "–"
    for unit, size in (("GiB", 1 << 30), ("MiB", 1 << 20), ("KiB", 1 << 10)):
        if abs(x) >= size:
            return "%.1f %s" % (x / float(size), unit)
    return "%d B" % x


def fmt_count(v):
    x = cells._num(v)
    return "–" if x is None else "{:,}".format(int(round(x)))


def fmt_ratio(v):
    return fmt_n(v, 2)


def _chapter_key(ctx, area):
    chapters = list(ctx.units.chapters)
    return (area == "all", chapters.index(area) if area in chapters else len(chapters), str(area))


def _host_key(state):
    return tr.STATE_ORDER.index(state) if state in tr.STATE_ORDER else len(tr.STATE_ORDER)


def limit_bytes(run, host):
    """The container's memory limit: host_memory.json's memory.max, else the sampler's mem_limit_bytes."""
    f = (host or {}).get("figures") or {}
    if cells._num(f.get("limit_bytes")) is not None:
        return f["limit_bytes"]
    lims = [cells._num(s.get("mem_limit_bytes")) for s in (run.samples or []) if isinstance(s, dict)]
    lims = [x for x in lims if x is not None and x > 0]
    return max(lims) if lims else None


def values(ctx, s, host, graphlet):
    """The placeholders and flags of closing.json's `resources` paragraphs."""
    run = ctx.run
    db = (((run.manifest or {}).get("index") or {}).get("database")) if isinstance(run.manifest, dict) else None
    index = INDEX_ON_DISK_BYTES.get(db)
    limit = limit_bytes(run, host)
    f = (host or {}).get("figures") or {}
    stated = tr.any_stated(s)
    return {
        "task_limit_known": limit is not None,
        "task_limit_gib": cells.fmt_gib(limit),
        "task_index_known": index is not None and limit is not None,
        "task_index_gib": cells.fmt_gib(index),
        "task_index_ratio": ("%.0f" % (index / float(limit))) if index and limit else "–",
        "task_top_n": str(tr.TOP_N),
        "task_any_stated": stated,
        "task_host_memory": bool(f),
        "task_host_before": ((f.get("before") or {}).get("utc")) or "–",
        "task_host_after": ((f.get("after") or {}).get("utc")) or "–",
        "task_page_bytes": fmt_count(f.get("page_bytes")),
        "task_graphlet_stated": bool(graphlet and graphlet.get("stating")),
        "task_graphlet_ops": fmt_count((graphlet or {}).get("ops")),
        "task_graphlet_stating": fmt_count((graphlet or {}).get("stating")),
    }


def _sources_text(counts):
    return "; ".join("%s %s" % (word, fmt_count(counts[key])) for key, word in SOURCE_WORDS if counts.get(key)) or "–"


def sources_table(s):
    rows = []
    for arm in [a for a in ARMS if a in s["executions"]] + sorted(a for a in s["executions"] if a not in ARMS):
        src = s["sources"].get(arm) or {}
        rows.append([esc(arm), fmt_count(s["executions"][arm])] +
                    [esc(_sources_text(src.get(k) or {})) for k in tr.SOURCE_FAMILIES])
    spec = latex.spec_p((30, 18, 42, 42, 42), aligns=["raggedright", "raggedleft", "raggedright", "raggedright",
                                                       "raggedright"])
    return latex.plain_table(spec, ["access path", "task executions", "CPU time from", "page faults from",
                                    "memory account from"], rows,
                             "where each access path's figures come from: the execution's recorded server block, "
                             "else its raw answer (an answer that states the figure per entry only: their sum, "
                             "their maximum for the memory account); not stated: no answer states it")


def _group_rows(ctx, s):
    """[(chapter, arm, host label, group_stats)] in reading order, the groups that state some figure."""
    out = []
    for area in sorted(s["by_area"], key=lambda a: _chapter_key(ctx, a)):
        arms = s["by_area"][area]
        for arm in [a for a in ARMS if a in arms] + sorted(a for a in arms if a not in ARMS):
            for state in sorted(arms[arm], key=_host_key):
                g = arms[arm][state]
                if g["cpu_ms"]["n"] or g["major_faults"]["n"] or g["memory_bytes"]["n"]:
                    out.append((area, arm, state, g))
    return out


def cpu_table(ctx, s):
    rows = []
    for area, arm, state, g in _group_rows(ctx, s):
        c, e, r = g["cpu_ms"], g["elapsed_ms"], g["cpu_ratio"]
        if not c["n"]:
            continue
        rows.append([esc(area), esc(arm), esc(HOST_WORDS.get(state, state)), "%d/%d" % (c["n"], g["executions"]),
                     fmt_t(c["p50"]), fmt_t(c["p90"]), fmt_t(c["max"]), fmt_t(e["p50"]), fmt_ratio(r["p10"]),
                     fmt_ratio(r["p50"])])
    return latex.plain_table("@{}l l l r r r r r r r@{}", ["chapter", "access path", "host", "n", "CPU p50", "p90",
                                                          "max", "elapsed p50", "CPU/el. p10", "p50"], rows,
                             "the engine's CPU time per task (timing.cpu\\_ms) per chapter, access path and host-cache "
                             "state, beside the engine's elapsed time of the same answers; CPU/el.: CPU time over "
                             "elapsed time per execution (well below 1: the task waited); n: the executions stating "
                             "it of the task executions; chapter all: every chapter")


def memory_table(ctx, s):
    rows = []
    for area, arm, state, g in _group_rows(ctx, s):
        mj, mb = g["major_faults"], g["memory_bytes"]
        rows.append([esc(area), esc(arm), esc(HOST_WORDS.get(state, state)), "%d/%d" % (mj["n"], g["executions"]),
                     fmt_count(mj["p50"]), fmt_count(mj["p90"]), fmt_count(mj["max"]),
                     "%d/%d" % (mb["n"], g["executions"]), fmt_bytes(mb["p50"]), fmt_bytes(mb["p90"]),
                     fmt_bytes(mb["max"])])
    return latex.plain_table("@{}l l l r r r r r r r r@{}", ["chapter", "access path", "host", "n", "major p50", "p90",
                                                            "max", "n", "memory p50", "p90", "max"], rows,
                             "major page faults per task (timing.page\\_faults.major: pages read from disk) and the "
                             "request's memory account (work.memory\\_bytes, its peak) per chapter, access path and "
                             "host-cache state; n: the executions stating each of the task executions")


def top_table(ctx, entries, figure):
    import sections   # sections imports this module: read lazily
    rows = []
    for i, e in enumerate(entries, 1):
        at = "%s r%s %s" % (e.get("arm"), e.get("repeat"), HOST_WORDS.get(e.get("host"), e.get("host")))
        head = [str(i), sections._link_sid(ctx, e["scenario"]), esc(e.get("area")), esc(at)]
        if figure == "cpu_ms":
            body = [fmt_t(e["value"]), fmt_t(e.get("elapsed_ms")), fmt_ratio(e.get("cpu_ratio")), fmt_t(e.get("cold")),
                    fmt_t(e.get("warm_median")), fmt_t(e.get("median"))]
        elif figure == "memory_bytes":
            body = [fmt_bytes(e["value"]), fmt_bytes(e.get("median"))]
        else:
            body = [fmt_count(e["value"]), fmt_count(e.get("cold")), fmt_count(e.get("warm_median")),
                    fmt_count(e.get("median"))]
        rows.append(head + body + [str(e.get("n"))])
    if figure == "cpu_ms":
        header, spec = ["CPU", "elapsed", "CPU/el.", "cold", "warm p50", "median"], "r r r r r r"
        cap = "the %d tests with the most CPU time per task" % tr.TOP_N
    elif figure == "memory_bytes":
        header, spec = ["memory", "median"], "r r"
        cap = "the %d tests with the largest memory account per task" % tr.TOP_N
    else:
        header, spec = ["major faults", "cold", "warm p50", "median"], "r r r r"
        cap = "the %d tests with the most major page faults per task" % tr.TOP_N
    cap += (": per test its execution with the largest figure (access path, repeat, host-cache state)%s, the median "
            "over its executions that state the figure and how many do" % (
                "" if figure == "memory_bytes" else ", the host-cold execution's figure apart from the warm "
                                                    "repeats' median"))
    return latex.plain_table("@{}r l l l %s r@{}" % spec, ["\\#", "test", "chapter", "execution"] + header + ["n"],
                             rows, cap)


def host_table(host):
    f = host["figures"]
    rows = []
    for k in tr.HOST_LEVELS:
        lv = f["levels"][k]
        if k == "host_available_kib":
            fmt = lambda v: fmt_gib1(v * 1024.0) if cells._num(v) is not None else "–"  # noqa: E731
        else:
            fmt = fmt_gib1
        rows.append([HOST_FIGURE_WORDS.get(k, esc(k)), fmt(lv["before"]), fmt(lv["after"]), ""])
    for k in tr.HOST_COUNTERS:
        d = f["deltas"].get(k)
        if k in tr.HOST_PAGE_COUNTERS:
            gib = f["pages_gib"].get(k)
            cell = "–" if d is None else "%s pages (%.1f GiB)" % (fmt_count(d), gib)
        else:
            cell = fmt_count(d)
        rows.append([HOST_FIGURE_WORDS.get(k, esc(k)), "", "", cell])
    spec = latex.spec_p((64, 22, 22, 50), aligns=["raggedright", "raggedleft", "raggedleft", "raggedleft"])
    out = latex.plain_table(spec, ["cgroup figure", "before", "after", "over the run"], rows,
                            "the engine container's memory (%s): levels at the run's start and end, counters as their "
                            "growth over the run" % esc(f.get("container") or "container not named"))
    for n in f.get("notes") or []:
        out += "%s.\n\n" % esc(n[:1].upper() + n[1:])
    return out


def fmt_gib1(v):
    x = cells._num(v)
    return "–" if x is None else "%.1f GiB" % (x / float(1 << 30))


def section(ctx):
    """Timing 8: the per-task figures of a manifest v8 run (LaTeX)."""
    run = ctx.run
    prose = ctx.closing.get("resources") if isinstance(ctx.closing.get("resources"), dict) else {}
    s = tr.of_run(run.dir, run.rows, run.manifest)
    host = tr.host_memory_of_run(run.dir)
    graphlet = tr.graphlet_service_figures(run.dir)
    vals = values(ctx, s, host, graphlet)
    out = latex.section("%s 8\\qd %s" % ("Timing", esc(prose.get("title") or TITLE)), LABEL)
    items = prose.get("paragraphs") or []
    if isinstance(items, str):
        items = [items]
    if items:
        out += paragraphs([ctx.fill(p, vals) for p in items])
    else:
        out += latex.placeholder("(resources per task: to be written in closing.json `resources`)") + "\n\n"
    out += sources_table(s)
    if vals["task_any_stated"]:
        if tr.any_stated(s, ("cpu_ms",)):
            out += cpu_table(ctx, s)
        else:
            out += ("No task execution of this run states the engine's CPU time (timing.cpu\\_ms): not stated.\n\n")
        out += memory_table(ctx, s)
        if not tr.any_stated(s, ("major_faults",)):
            out += ("No task execution of this run states its page faults (timing.page\\_faults, stated on Linux): "
                    "not stated.\n\n")
        for figure in tr.TOP_FIGURES:
            if s["top"].get(figure):
                out += top_table(ctx, s["top"][figure], figure)
    else:
        out += ("No task execution of this run states the engine's CPU time, page faults or memory account: the "
                "figures per task are not stated (an engine before these fields).\n\n")
    if host.get("figures"):
        out += host_table(host)
    elif host.get("problem"):
        out += "The run's host\\_memory.json is not read: %s.\n\n" % esc(host["problem"])
    return out
