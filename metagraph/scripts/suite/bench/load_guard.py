#!/usr/bin/env python3
"""The bench run's load guard (lane B): no performance verdict from a block that ran on a loaded host.

    load_guard.py RUN_DIR [--sra-max PCT] [--load1-max L] [--allow-missing-sra] [--max-gap-s S] [--pad-s S] [--json]

It reads the resource samples the suite's sampler took on mex during the run (sampler.py: one `docker stats` read
of the staging container and the Logan container sra_all, and the host's load1, every ~2 s) and judges every block
of the run's schedule (block 0 the warm-up pass, blocks 1..6 the measured cycles) over the time its requests ran:

  clean      every sample in the block's window is at or below both thresholds (the Logan container's CPU,
             default 100 %, one core: the Logan load ran at ~5,000-14,000 %; the host's load1, default 16), and the
             samples cover the window without a gap longer than --max-gap-s (15 s);
  loaded     at least one sample in the window exceeds a threshold (its breaches are listed);
  unguarded  the load during the window is unknown: no samples, a gap, or a sample without a reading (docker did not
             answer; the Logan container not reported, unless --allow-missing-sra declares it stopped on purpose,
             which counts a sample whose staging container was read as 0 % for sra_all).

A block's window is the union of its executions' [start, end] intervals, each widened by --pad-s (30 s), so the load
between a block's requests and just around them counts too; the samples must cover the window less that padding (the
executions themselves). A block a resume finished later has two windows. Executions on mex
(direct-local) are stamped with mex's clock, the others with the orchestrator's; samples carry mex's clock: the run's
recorded clock offset (bench.json) puts all on the orchestrator's. The run's verdict is clean only when every block
is clean; analyze.py emits a performance verdict only then. Thresholds default to the run's own (bench.json `guard`),
and what was used is part of every verdict. Stdlib only.
"""

import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bench_lib as bl  # noqa: E402

sc = bl.sc

DEFAULTS = {"sra_all_cpu_pct_max": 100.0, "host_load1_max": 16.0, "allow_missing_sra": False, "max_gap_s": 15.0,
            "pad_s": 30.0}
MAX_BREACHES_LISTED = 10


def thresholds(overrides=None, base=None):
    """DEFAULTS, then |base| (a run's bench.json `guard`), then the non-None |overrides|."""
    out = dict(DEFAULTS)
    for src in (base or {}, overrides or {}):
        for k, v in src.items():
            if k in DEFAULTS and v is not None:
                out[k] = v
    return out


def sample_time(sample, offset_s=0.0):
    """A sample's time on the orchestrator's clock (epoch s): its mex time less mex's offset; None when unstamped."""
    t = sc.parse_utc(sample.get("t_utc")) if isinstance(sample, dict) else None
    return None if t is None else t - (offset_s or 0.0)


def _num(v):
    return v if isinstance(v, (int, float)) and not isinstance(v, bool) else None


def reading(sample, limits):
    """(sra_all CPU %, load1) of a sample, None for what it does not state (an unreported Logan container counts as
    0 % when the limits allow a missing one and the staging container was read)."""
    sra, load1 = _num(sample.get("sra_all_cpu_pct")), _num(sample.get("host_load1"))
    if sra is None and limits.get("allow_missing_sra") and _num(sample.get("cpu_pct")) is not None:
        sra = 0.0
    return sra, load1


def judge_samples(samples, limits):
    """The verdict on a set of samples (no coverage rule): {verdict, reason, samples, sra_all_cpu_pct_max,
    host_load1_max, breaches (the first MAX_BREACHES_LISTED), breach_count, unread}."""
    limits = thresholds(limits)
    out = {"samples": len(samples), "sra_all_cpu_pct_max": None, "host_load1_max": None, "breaches": [],
           "breach_count": 0, "unread": 0}
    for s in samples:
        sra, load1 = reading(s, limits)
        if sra is None or load1 is None:
            out["unread"] += 1
        if sra is not None:
            out["sra_all_cpu_pct_max"] = max(sra, out["sra_all_cpu_pct_max"] or sra)
        if load1 is not None:
            out["host_load1_max"] = max(load1, out["host_load1_max"] or load1)
        over = []
        if sra is not None and sra > limits["sra_all_cpu_pct_max"]:
            over.append("sra_all_cpu_pct")
        if load1 is not None and load1 > limits["host_load1_max"]:
            over.append("host_load1")
        if over:
            out["breach_count"] += 1
            if len(out["breaches"]) < MAX_BREACHES_LISTED:
                out["breaches"].append({"t_utc": s.get("t_utc"), "sra_all_cpu_pct": sra, "host_load1": load1,
                                        "over": over})
    if not samples:
        out.update(verdict="unguarded", reason="no samples")
    elif out["breach_count"]:
        out.update(verdict="loaded", reason="%d sample(s) above the thresholds (sra_all %s %%, load1 %s)" % (
            out["breach_count"], limits["sra_all_cpu_pct_max"], limits["host_load1_max"]))
    elif out["unread"]:
        out.update(verdict="unguarded", reason="%d sample(s) without a reading" % out["unread"])
    else:
        out.update(verdict="clean", reason="every sample at or below the thresholds")
    return out


def merge_windows(intervals, pad_s):
    """The union of the [t0, t1] intervals each widened by pad_s, as sorted disjoint [t0, t1]."""
    spans = sorted([float(a) - pad_s, float(b) + pad_s] for a, b in intervals if a is not None and b is not None)
    out = []
    for a, b in spans:
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return out


def coverage_gaps(times, windows, max_gap_s, pad_s=0.0):
    """The gaps (s) longer than max_gap_s between a window's start, its samples' times and its end, the window less
    its padding (the samples must cover the executions, not the margin around them)."""
    gaps = []
    for w0, w1 in windows:
        w0, w1 = w0 + pad_s, w1 - pad_s
        inside = sorted(t for t in times if w0 <= t <= w1)
        marks = [w0] + inside + [w1]
        for a, b in zip(marks, marks[1:]):
            if b - a > max_gap_s:
                gaps.append({"from": sc.iso_utc(a), "to": sc.iso_utc(b), "s": round(b - a, 1)})
    return gaps


def judge_block(samples, intervals, limits, offset_s=0.0):
    """The verdict on one block: judge_samples over the samples inside its windows, `unguarded` on a coverage gap
    (unless already `loaded`)."""
    limits = thresholds(limits)
    windows = merge_windows(intervals, limits["pad_s"])
    timed = [(sample_time(s, offset_s), s) for s in samples]
    inside = [s for t, s in timed if t is not None and any(w0 <= t <= w1 for w0, w1 in windows)]
    out = judge_samples(inside, limits)
    out["windows"] = [[sc.iso_utc(a), sc.iso_utc(b)] for a, b in windows]
    out["gaps"] = coverage_gaps([t for t, _ in timed if t is not None], windows, limits["max_gap_s"], limits["pad_s"])
    if not windows:
        out.update(verdict="unguarded", reason="no execution in the block")
    elif out["verdict"] == "clean" and out["gaps"]:
        out.update(verdict="unguarded", reason="%d gap(s) in the samples longer than %s s" % (
            len(out["gaps"]), limits["max_gap_s"]))
    return out


def judge_blocks(blocks, samples, limits, offset_s=0.0):
    """{blocks: {block: verdict}, verdict, loaded, unguarded, limits} for {block: [[t0, t1], ...]}."""
    limits = thresholds(limits)
    per = dict((str(b), judge_block(samples, iv, limits, offset_s)) for b, iv in sorted(blocks.items()))
    loaded = sorted((b for b, v in per.items() if v["verdict"] == "loaded"), key=int)
    unguarded = sorted((b for b, v in per.items() if v["verdict"] == "unguarded"), key=int)
    verdict = "loaded" if loaded else ("unguarded" if unguarded or not per else "clean")
    return {"verdict": verdict, "blocks": per, "loaded": loaded, "unguarded": unguarded, "limits": limits,
            "clock_offset_s": offset_s}


# ---------------------------------------------------------------------------------------------------------------
# a run directory

def read_samples(run_dir):
    """Every sample of the run (RUN_DIR/samples/*.jsonl), deduplicated by (t_utc, phase), in time order."""
    out, seen = [], set()
    d = os.path.join(run_dir, "samples")
    if os.path.isdir(d):
        for name in sorted(os.listdir(d)):
            if not name.endswith(".jsonl"):
                continue
            rows, _ = sc.read_jsonl(os.path.join(d, name))
            for r in rows:
                key = (r.get("t_utc"), r.get("phase"))
                if key not in seen:
                    seen.add(key)
                    out.append(r)
    out.sort(key=lambda r: r.get("t_utc") or "")
    return out


def clock_offset(bench):
    """mex's clock offset (s, mex ahead) of the run: the mean of its sessions' recorded offsets, 0 when none."""
    vals = []
    for sess in (bench or {}).get("sessions") or []:
        for key in ("clock_offset_start_s", "clock_offset_end_s"):
            if _num(sess.get(key)) is not None:
                vals.append(sess[key])
    return sum(vals) / len(vals) if vals else 0.0


def execution_interval(line, offset_s):
    """An execution's [start, end] on the orchestrator's clock (direct-local lines carry mex's clock)."""
    t0, t1 = sc.parse_utc(line.get("started_utc")), sc.parse_utc(line.get("ended_utc"))
    if t0 is None:
        return None
    if t1 is None:
        t1 = t0 + (_num(line.get("wall_ms")) or 0.0) / 1000.0
    if line.get("arm") == "direct-local":
        t0, t1 = t0 - offset_s, t1 - offset_s
    return [t0, t1]


def run_blocks(run_dir, offset_s):
    """{block: [[t0, t1], ...]} of the run's recorded executions, by the schedule's (bench id, arm) -> block."""
    sched = sc.read_json(os.path.join(run_dir, "schedule.json")) or {}
    block_of = dict(((s["bench_id"], s["arm"]), s["block"]) for s in sched.get("steps") or [])
    blocks = {}
    for arm in bl.ARMS:
        rows, _ = sc.read_jsonl(os.path.join(run_dir, "arms", arm, "executions.jsonl"))
        for r in rows:
            b = block_of.get((r.get("scenario"), r.get("arm")))
            iv = execution_interval(r, offset_s)
            if b is not None and iv is not None:
                blocks.setdefault(b, []).append(iv)
    return blocks


def judge_run(run_dir, overrides=None):
    """The load guard's verdict on a bench run directory."""
    bench = sc.read_json(os.path.join(run_dir, "bench.json")) or {}
    limits = thresholds(overrides, bench.get("guard"))
    offset = clock_offset(bench)
    out = judge_blocks(run_blocks(run_dir, offset), read_samples(run_dir), limits, offset)
    out["overrides"] = dict((k, v) for k, v in (overrides or {}).items() if v is not None)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("run_dir")
    ap.add_argument("--sra-max", type=float, default=None, help="the Logan container's CPU %% limit")
    ap.add_argument("--load1-max", type=float, default=None, help="the host's load1 limit")
    ap.add_argument("--allow-missing-sra", action="store_true", default=None,
                    help="a sample that does not report sra_all (the container stopped on purpose) counts as 0 %%")
    ap.add_argument("--max-gap-s", type=float, default=None)
    ap.add_argument("--pad-s", type=float, default=None)
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args(argv)
    res = judge_run(args.run_dir, {"sra_all_cpu_pct_max": args.sra_max, "host_load1_max": args.load1_max,
                                   "allow_missing_sra": args.allow_missing_sra, "max_gap_s": args.max_gap_s,
                                   "pad_s": args.pad_s})
    if args.json:
        print(json.dumps(res, indent=1, sort_keys=True))
    else:
        for b, v in res["blocks"].items():
            print("block %s: %s (%s; %d samples, sra_all max %s %%, load1 max %s)" % (
                b, v["verdict"], v["reason"], v["samples"], v["sra_all_cpu_pct_max"], v["host_load1_max"]))
        print("run: %s" % res["verdict"])
    return 0 if res["verdict"] == "clean" else 1


if __name__ == "__main__":
    sys.exit(main())
