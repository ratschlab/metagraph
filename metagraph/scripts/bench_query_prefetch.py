#!/usr/bin/env python3
"""Benchmark `metagraph server_query` query latency, broken down per stage.

Spawns a `metagraph server_query` process against a multi-graph CSV (or a
single -i/-a pair), sends one or more `/search` requests over HTTP, and
parses the server's trace logs to summarize:

  * "[Query graph construction] Contigs mapped to the full graph ..."
        -> k-mer mapping latency per batch
  * "RD query [...] -- traversal: ... call_rd_rows: ..."
        -> annotation (row-diff) traversal + decoding latency per batch

Designed for A/B testing optimizations such as MADV_WILLNEED prefetch:
run the same workload before and after a change, optionally with the OS
page cache dropped between runs, and diff the resulting summaries (or
the `--json` output).

Stdlib only -- no `requests` / `pandas` dependency.

Example:

    scripts/bench_query_prefetch.py \\
        --metagraph build/metagraph/metagraph \\
        --graphs graphs.csv \\
        --query reads.fa \\
        --mmap --madv-random \\
        -p 4 --threads-each 8 \\
        --warmup 1 --repeats 3 \\
        --json results_after.json

Add `--drop-cache` to evict the OS page cache between repeats (requires
root: writes `/proc/sys/vm/drop_caches` on Linux, runs `purge` on macOS).
Without root the flag is a no-op with a warning -- the rest of the script
runs unprivileged.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import re
import shlex
import signal
import statistics
import subprocess
import sys
import tempfile
import threading
import time
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Trace-line regexes. Anchored to the substantive part of the log line so
# changes to the "[YYYY-MM-DD HH:MM:SS.mmm] [trace] " prefix don't break us.
# ---------------------------------------------------------------------------

KMER_MAP_RE = re.compile(
    r"\[Query graph construction\] Contigs mapped to the full graph "
    r"\[threads: (?P<threads>\d+), contigs: (?P<contigs>\d+), "
    r"chunk_size: (?P<chunk_size>\d+)\] "
    r"\(found (?P<found>\d+) / (?P<total>\d+) k-mers\) "
    r"in (?P<sec>[\d.]+) sec"
)

RD_QUERY_RE = re.compile(
    r"RD query \[threads: (?P<threads>\d+), "
    r"rows: (?P<rows_in>\d+) -> (?P<rows_out>\d+) "
    r"\((?P<expansion>[\d.]+)x\)\] -- "
    r"traversal: (?P<traversal>[\d.]+) sec, "
    r"call_rd_rows: (?P<call_rd_rows>[\d.]+) sec "
    r"\(set bits: (?P<set_bits>\d+), capacity: (?P<capacity>\d+)\), "
    r"decoding: (?P<decoding>[\d.]+) sec, "
    r"reconstruction: (?P<reconstruction>[\d.]+) sec"
)

# server.cpp emits this once all graphs have been loaded.
SERVER_READY_RE = re.compile(r"All graphs were loaded.*Ready to serve queries")

# server.cpp also emits this single-graph variant when -i/-a is used (no CSV).
SINGLE_GRAPH_READY_RE = re.compile(r"\[Server\] Will listen on")


# ---------------------------------------------------------------------------
# Stats helpers
# ---------------------------------------------------------------------------


def summarize(xs: List[float]) -> Dict[str, float]:
    """Return count / mean / median / p50 / p90 / p99 / max / total of `xs`."""
    if not xs:
        return {"count": 0}
    s = sorted(xs)
    return {
        "count": len(xs),
        "mean": statistics.fmean(xs),
        "median": statistics.median(s),
        "p50": _percentile(s, 50),
        "p90": _percentile(s, 90),
        "p99": _percentile(s, 99),
        "max": s[-1],
        "total": sum(xs),
    }


def _percentile(sorted_xs: List[float], pct: float) -> float:
    if not sorted_xs:
        return 0.0
    if len(sorted_xs) == 1:
        return sorted_xs[0]
    # linear interpolation between closest ranks (numpy's default)
    k = (len(sorted_xs) - 1) * (pct / 100.0)
    lo = int(k)
    hi = min(lo + 1, len(sorted_xs) - 1)
    frac = k - lo
    return sorted_xs[lo] * (1 - frac) + sorted_xs[hi] * frac


def fmt_summary(label: str, units: str, s: Dict[str, float]) -> str:
    if not s.get("count"):
        return f"  {label:24s}  (no events)"
    return (
        f"  {label:24s} "
        f"n={s['count']:<5d}  "
        f"mean={s['mean']:.3f}{units}  "
        f"med={s['median']:.3f}{units}  "
        f"p99={s['p99']:.3f}{units}  "
        f"max={s['max']:.3f}{units}  "
        f"total={s['total']:.2f}{units}"
    )


# ---------------------------------------------------------------------------
# Server lifecycle
# ---------------------------------------------------------------------------


@dataclass
class ServerProcess:
    proc: subprocess.Popen
    log_path: Path
    drain_thread: threading.Thread
    log_buffer: List[str] = field(default_factory=list)
    log_lock: threading.Lock = field(default_factory=threading.Lock)

    def snapshot(self) -> str:
        """Return everything captured from the server's stderr so far."""
        with self.log_lock:
            return "".join(self.log_buffer)

    def stop(self, timeout: float = 10.0) -> None:
        if self.proc.poll() is None:
            self.proc.send_signal(signal.SIGTERM)
            try:
                self.proc.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                self.proc.kill()
                self.proc.wait()
        self.drain_thread.join(timeout=2.0)


def spawn_server(
    metagraph: str,
    graphs_csv: Optional[str],
    single_index: Optional[Tuple[str, str]],
    port: int,
    host: str,
    mmap: bool,
    madv_random: bool,
    parallel: int,
    threads_each: int,
    extra_args: List[str],
    log_dir: Path,
    ready_timeout: float,
) -> ServerProcess:
    """Spawn `metagraph server_query`, return once it's accepting connections."""
    cmd: List[str] = [metagraph, "server_query"]
    if graphs_csv is not None:
        cmd.append(graphs_csv)
    elif single_index is not None:
        cmd += ["-i", single_index[0], "-a", single_index[1]]
    else:
        raise ValueError("Either graphs_csv or single_index must be provided.")

    cmd += ["--port", str(port)]
    if host:
        cmd += ["--address", host]
    cmd += ["-p", str(parallel), "--threads-each", str(threads_each)]
    if mmap:
        cmd.append("--mmap")
    if madv_random:
        cmd.append("--madv-random")
    cmd.append("-v")  # required for trace lines to be emitted
    cmd += extra_args

    log_path = log_dir / "server.log"
    print(f"[bench] spawning: {' '.join(shlex.quote(c) for c in cmd)}", file=sys.stderr)
    print(f"[bench] server log -> {log_path}", file=sys.stderr)

    log_file = open(log_path, "w", buffering=1)  # line-buffered
    # Note: server logs to stderr; merge stdout+stderr for completeness.
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        bufsize=1,
        text=True,
        # New process group so SIGTERM only hits the server, not us.
        preexec_fn=os.setsid if hasattr(os, "setsid") else None,
    )

    server = ServerProcess(proc=proc, log_path=log_path, drain_thread=None)  # type: ignore[arg-type]

    def drain() -> None:
        # Stream the server's combined stderr+stdout into both the in-memory
        # buffer (for ready-detection and parsing) and the on-disk log.
        assert proc.stdout is not None
        for line in proc.stdout:
            with server.log_lock:
                server.log_buffer.append(line)
            log_file.write(line)
        log_file.close()

    drain_thread = threading.Thread(target=drain, daemon=True, name="server-drain")
    drain_thread.start()
    server.drain_thread = drain_thread

    # Multi-graph servers print "Ready to serve queries"; single-graph ones load
    # asynchronously and only emit "Will listen on ...". For the latter we also
    # need to retry the actual HTTP request because /search returns 503-ish
    # while the lazy load is in flight; we handle that in send_query().
    ready_re = SERVER_READY_RE if graphs_csv is not None else SINGLE_GRAPH_READY_RE
    deadline = time.monotonic() + ready_timeout
    while time.monotonic() < deadline:
        if proc.poll() is not None:
            tail = server.snapshot()[-2000:]
            raise RuntimeError(
                f"Server exited (code {proc.returncode}) before ready.\n"
                f"--- last log ---\n{tail}"
            )
        if ready_re.search(server.snapshot()):
            print(f"[bench] server ready on {host or '0.0.0.0'}:{port}", file=sys.stderr)
            return server
        time.sleep(0.2)

    server.stop()
    raise TimeoutError(
        f"Server did not become ready within {ready_timeout:.1f}s. "
        f"See {log_path} for details."
    )


# ---------------------------------------------------------------------------
# Cache drop (best-effort; per-platform)
# ---------------------------------------------------------------------------


def drop_page_cache() -> Optional[str]:
    """Try to drop the OS page cache. Returns None on success, else an error
    message. Best-effort: silently no-ops if the script lacks privileges."""
    sys_name = platform.system()
    try:
        if sys_name == "Linux":
            subprocess.run(["sync"], check=True)
            with open("/proc/sys/vm/drop_caches", "w") as f:
                f.write("3\n")
            return None
        elif sys_name == "Darwin":
            r = subprocess.run(["purge"], check=False, capture_output=True)
            if r.returncode != 0:
                return f"`purge` failed: {r.stderr.decode().strip()}"
            return None
        else:
            return f"unsupported OS: {sys_name}"
    except PermissionError:
        return "permission denied (run as root)"
    except FileNotFoundError as e:
        return f"missing tool: {e}"
    except Exception as e:  # pragma: no cover
        return repr(e)


# ---------------------------------------------------------------------------
# Query
# ---------------------------------------------------------------------------


def fasta_payload(query_path: Path, max_seqs: Optional[int]) -> str:
    """Read `query_path` (FASTA/FASTQ) and return its contents as FASTA text.
    Optionally cap to first `max_seqs` records."""
    text = query_path.read_text()
    if not max_seqs:
        return text
    out: List[str] = []
    seen = 0
    in_seq = False
    for line in text.splitlines():
        if line.startswith(">"):
            seen += 1
            if seen > max_seqs:
                break
            in_seq = True
            out.append(line)
        elif in_seq:
            out.append(line)
    return "\n".join(out) + "\n"


def send_query(
    host: str,
    port: int,
    fasta: str,
    discovery_fraction: float,
    top_labels: int,
    graph_names: Optional[List[str]],
    timeout: float,
    initial_retries: int,
) -> Tuple[float, int]:
    """POST /search with `fasta`. Returns (wall_seconds, http_status)."""
    payload: Dict[str, Any] = {
        "FASTA": fasta,
        "discovery_fraction": discovery_fraction,
        "top_labels": top_labels,
        "count_labels": True,
    }
    if graph_names is not None:
        payload["graphs"] = graph_names

    body = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(
        f"http://{host}:{port}/search",
        data=body,
        method="POST",
        headers={"Content-Type": "application/json"},
    )

    last_err: Optional[Exception] = None
    for attempt in range(initial_retries + 1):
        t0 = time.monotonic()
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                resp.read()
                return time.monotonic() - t0, resp.status
        except urllib.error.HTTPError as e:
            # 503 etc. while server is still warming up -- retry briefly.
            last_err = e
            if attempt < initial_retries:
                time.sleep(0.5)
                continue
            raise
        except (urllib.error.URLError, ConnectionResetError) as e:
            last_err = e
            if attempt < initial_retries:
                time.sleep(0.5)
                continue
            raise
    raise RuntimeError(f"unreachable: {last_err!r}")


# ---------------------------------------------------------------------------
# Log parsing
# ---------------------------------------------------------------------------


@dataclass
class KmerMapEvent:
    threads: int
    contigs: int
    chunk_size: int
    found: int
    total: int
    sec: float

    def to_dict(self) -> Dict[str, Any]:
        return self.__dict__.copy()


@dataclass
class RdQueryEvent:
    threads: int
    rows_in: int
    rows_out: int
    expansion: float
    traversal: float
    call_rd_rows: float
    set_bits: int
    capacity: int
    decoding: float
    reconstruction: float

    def to_dict(self) -> Dict[str, Any]:
        return self.__dict__.copy()


def parse_log(text: str) -> Tuple[List[KmerMapEvent], List[RdQueryEvent]]:
    kmer_events: List[KmerMapEvent] = []
    rd_events: List[RdQueryEvent] = []
    for m in KMER_MAP_RE.finditer(text):
        d = m.groupdict()
        kmer_events.append(KmerMapEvent(
            threads=int(d["threads"]),
            contigs=int(d["contigs"]),
            chunk_size=int(d["chunk_size"]),
            found=int(d["found"]),
            total=int(d["total"]),
            sec=float(d["sec"]),
        ))
    for m in RD_QUERY_RE.finditer(text):
        d = m.groupdict()
        rd_events.append(RdQueryEvent(
            threads=int(d["threads"]),
            rows_in=int(d["rows_in"]),
            rows_out=int(d["rows_out"]),
            expansion=float(d["expansion"]),
            traversal=float(d["traversal"]),
            call_rd_rows=float(d["call_rd_rows"]),
            set_bits=int(d["set_bits"]),
            capacity=int(d["capacity"]),
            decoding=float(d["decoding"]),
            reconstruction=float(d["reconstruction"]),
        ))
    return kmer_events, rd_events


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    src = p.add_argument_group("server input (one required)")
    src.add_argument("--graphs", type=Path,
                     help="CSV with rows '<name>,<graph>,<annotation>' for "
                          "multi-graph mode")
    src.add_argument("--single-index", nargs=2, metavar=("GRAPH", "ANNO"),
                     help="single (graph, annotation) pair instead of a CSV")

    srv = p.add_argument_group("server")
    srv.add_argument("--metagraph", required=True,
                     help="path to the metagraph binary")
    srv.add_argument("--port", type=int, default=5555)
    srv.add_argument("--host", default="127.0.0.1",
                     help="server bind address (default: 127.0.0.1)")
    srv.add_argument("--mmap", action="store_true",
                     help="pass --mmap to server_query")
    srv.add_argument("--madv-random", action="store_true",
                     help="pass --madv-random to server_query")
    srv.add_argument("-p", "--parallel", type=int, default=1,
                     help="--parallel: max parallel connections")
    srv.add_argument("--threads-each", type=int, default=1,
                     help="--threads-each: threads per graph")
    srv.add_argument("--server-arg", action="append", default=[],
                     help="extra arg forwarded to server_query (repeatable)")
    srv.add_argument("--ready-timeout", type=float, default=600.0,
                     help="seconds to wait for the server to become ready")

    qry = p.add_argument_group("query")
    qry.add_argument("--query", required=True, type=Path,
                     help="FASTA file to send as the query body")
    qry.add_argument("--graph-name", action="append",
                     help="restrict query to these graph names "
                          "(multi-graph mode; repeatable). Default: all.")
    qry.add_argument("--discovery-fraction", type=float, default=0.0)
    qry.add_argument("--top-labels", type=int, default=10000)
    qry.add_argument("--max-seqs", type=int, default=0,
                     help="cap query to first N sequences (0 = all)")
    qry.add_argument("--query-timeout", type=float, default=3600.0,
                     help="HTTP timeout per request")

    bnch = p.add_argument_group("benchmark")
    bnch.add_argument("--warmup", type=int, default=0,
                      help="number of warmup queries (not counted)")
    bnch.add_argument("--repeats", type=int, default=1,
                      help="number of timed queries")
    bnch.add_argument("--drop-cache", action="store_true",
                      help="best-effort drop of OS page cache before each run "
                           "(needs root; warns if it fails)")

    out = p.add_argument_group("output")
    out.add_argument("--json", type=Path,
                     help="write summary JSON to this path")
    out.add_argument("--log-dir", type=Path,
                     help="dir to keep server logs (default: temp, deleted)")

    args = p.parse_args()
    if (args.graphs is None) == (args.single_index is None):
        p.error("exactly one of --graphs or --single-index is required")
    return args


def main() -> int:
    args = parse_args()

    if args.drop_cache:
        err = drop_page_cache()
        if err:
            print(f"[bench] WARN: cache drop unavailable: {err}", file=sys.stderr)

    # --log-dir = persistent; otherwise a TemporaryDirectory we clean up.
    if args.log_dir is not None:
        args.log_dir.mkdir(parents=True, exist_ok=True)
        log_dir_ctx: Any = _NoopCtx(args.log_dir)
    else:
        log_dir_ctx = tempfile.TemporaryDirectory(prefix="bench_metagraph_")

    with log_dir_ctx as log_dir_str:
        log_dir = Path(log_dir_str)

        server = spawn_server(
            metagraph=args.metagraph,
            graphs_csv=str(args.graphs) if args.graphs else None,
            single_index=tuple(args.single_index) if args.single_index else None,
            port=args.port,
            host=args.host,
            mmap=args.mmap,
            madv_random=args.madv_random,
            parallel=args.parallel,
            threads_each=args.threads_each,
            extra_args=args.server_arg,
            log_dir=log_dir,
            ready_timeout=args.ready_timeout,
        )

        try:
            fasta = fasta_payload(args.query, args.max_seqs)
            print(f"[bench] query payload: {len(fasta)} bytes", file=sys.stderr)

            wall_times: List[float] = []
            log_marks: List[int] = []  # log offset before each TIMED request

            # Warmup
            for i in range(args.warmup):
                wall, status = send_query(
                    host=args.host, port=args.port, fasta=fasta,
                    discovery_fraction=args.discovery_fraction,
                    top_labels=args.top_labels,
                    graph_names=args.graph_name,
                    timeout=args.query_timeout,
                    initial_retries=20,  # retry while async load completes
                )
                print(f"[bench] warmup {i+1}/{args.warmup}: "
                      f"{wall:.2f}s (HTTP {status})", file=sys.stderr)

            # Timed runs
            for i in range(args.repeats):
                if args.drop_cache and i > 0:  # cache already (just) dropped pre-spawn
                    err = drop_page_cache()
                    if err:
                        print(f"[bench] WARN: cache drop failed: {err}",
                              file=sys.stderr)
                # Mark log offset so we can attribute trace lines to this run.
                log_marks.append(len(server.snapshot()))
                wall, status = send_query(
                    host=args.host, port=args.port, fasta=fasta,
                    discovery_fraction=args.discovery_fraction,
                    top_labels=args.top_labels,
                    graph_names=args.graph_name,
                    timeout=args.query_timeout,
                    initial_retries=2,
                )
                wall_times.append(wall)
                print(f"[bench] run {i+1}/{args.repeats}: "
                      f"{wall:.2f}s (HTTP {status})", file=sys.stderr)

            log_marks.append(len(server.snapshot()))

        finally:
            # Give the server a moment to flush trailing trace lines, then stop.
            time.sleep(0.5)
            server.stop()

        # --- parse + report --------------------------------------------------
        full_log = server.snapshot()

        # Per-run aggregates: only count events that landed between consecutive marks.
        per_run_kmer: List[List[KmerMapEvent]] = []
        per_run_rd: List[List[RdQueryEvent]] = []
        for i in range(args.repeats):
            chunk = full_log[log_marks[i]:log_marks[i + 1]]
            k_ev, r_ev = parse_log(chunk)
            per_run_kmer.append(k_ev)
            per_run_rd.append(r_ev)

        # Pooled across all timed runs.
        all_kmer = [e for run in per_run_kmer for e in run]
        all_rd = [e for run in per_run_rd for e in run]

        result = build_result(args, wall_times, all_kmer, all_rd, per_run_kmer, per_run_rd)
        print_human_report(result)

        if args.json is not None:
            args.json.write_text(json.dumps(result, indent=2, default=str) + "\n")
            print(f"[bench] wrote {args.json}", file=sys.stderr)

    return 0


def build_result(
    args: argparse.Namespace,
    wall_times: List[float],
    all_kmer: List[KmerMapEvent],
    all_rd: List[RdQueryEvent],
    per_run_kmer: List[List[KmerMapEvent]],
    per_run_rd: List[List[RdQueryEvent]],
) -> Dict[str, Any]:
    return {
        "config": {
            "graphs": str(args.graphs) if args.graphs else None,
            "single_index": list(args.single_index) if args.single_index else None,
            "query": str(args.query),
            "graph_names": args.graph_name,
            "mmap": args.mmap,
            "madv_random": args.madv_random,
            "parallel": args.parallel,
            "threads_each": args.threads_each,
            "warmup": args.warmup,
            "repeats": args.repeats,
            "drop_cache": args.drop_cache,
            "discovery_fraction": args.discovery_fraction,
            "top_labels": args.top_labels,
        },
        "client_wall_sec_per_run": wall_times,
        "client_wall_sec": summarize(wall_times),
        "kmer_mapping": {
            "events_per_run": [len(r) for r in per_run_kmer],
            "sec": summarize([e.sec for e in all_kmer]),
            "contigs": summarize([e.contigs for e in all_kmer]),
            "found_ratio": summarize([
                e.found / e.total if e.total else 0.0 for e in all_kmer
            ]),
        },
        "rd_query": {
            "events_per_run": [len(r) for r in per_run_rd],
            "rows_in": summarize([e.rows_in for e in all_rd]),
            "rows_out": summarize([e.rows_out for e in all_rd]),
            "expansion": summarize([e.expansion for e in all_rd]),
            "traversal_sec": summarize([e.traversal for e in all_rd]),
            "call_rd_rows_sec": summarize([e.call_rd_rows for e in all_rd]),
            "decoding_sec": summarize([e.decoding for e in all_rd]),
            "reconstruction_sec": summarize([e.reconstruction for e in all_rd]),
            "set_bits": summarize([e.set_bits for e in all_rd]),
            "capacity": summarize([e.capacity for e in all_rd]),
        },
        "raw": {
            "kmer_mapping_events": [e.to_dict() for e in all_kmer],
            "rd_query_events": [e.to_dict() for e in all_rd],
        },
    }


def print_human_report(r: Dict[str, Any]) -> None:
    cfg = r["config"]
    print()
    print("=" * 72)
    print("metagraph server_query benchmark")
    print("=" * 72)
    print(f"  query:       {cfg['query']}")
    if cfg["graphs"]:
        print(f"  graphs csv:  {cfg['graphs']}")
    elif cfg["single_index"]:
        print(f"  index:       {cfg['single_index'][0]}")
        print(f"  annotation:  {cfg['single_index'][1]}")
    if cfg["graph_names"]:
        print(f"  graph names: {', '.join(cfg['graph_names'])}")
    print(f"  flags:       --parallel={cfg['parallel']} "
          f"--threads-each={cfg['threads_each']} "
          f"--mmap={cfg['mmap']} --madv-random={cfg['madv_random']}")
    print(f"  runs:        warmup={cfg['warmup']}  repeats={cfg['repeats']}  "
          f"drop_cache={cfg['drop_cache']}")
    print()

    print("== Client wall time per request ==")
    print(fmt_summary("wall", "s", r["client_wall_sec"]))
    print()

    print("== K-mer mapping  (Query graph construction / Contigs mapped) ==")
    print(fmt_summary("time per event",       "s", r["kmer_mapping"]["sec"]))
    print(fmt_summary("contigs per event",    "",  r["kmer_mapping"]["contigs"]))
    print(fmt_summary("found-ratio",          "",  r["kmer_mapping"]["found_ratio"]))
    print(f"  events per run: {r['kmer_mapping']['events_per_run']}")
    print()

    print("== Annotation queries  (RD query) ==")
    print(fmt_summary("traversal time",       "s", r["rd_query"]["traversal_sec"]))
    print(fmt_summary("call_rd_rows time",    "s", r["rd_query"]["call_rd_rows_sec"]))
    print(fmt_summary("decoding time",        "s", r["rd_query"]["decoding_sec"]))
    print(fmt_summary("reconstruction time",  "s", r["rd_query"]["reconstruction_sec"]))
    print(fmt_summary("rows in",              "",  r["rd_query"]["rows_in"]))
    print(fmt_summary("rows out",             "",  r["rd_query"]["rows_out"]))
    print(fmt_summary("expansion",            "x", r["rd_query"]["expansion"]))
    print(f"  events per run: {r['rd_query']['events_per_run']}")
    print()


class _NoopCtx:
    """Mimics tempfile.TemporaryDirectory's context but for a fixed path."""
    def __init__(self, p: Path) -> None:
        self._p = str(p)
    def __enter__(self) -> str:
        return self._p
    def __exit__(self, *exc: Any) -> None:
        return None


if __name__ == "__main__":
    sys.exit(main())
