#!/usr/bin/env python3
"""Resource sampler of the suite (DESIGN §5, samples schema §8.4): runs on mex next to the direct runner.

Every --interval seconds (start to start) one READ-ONLY `docker stats --no-stream` call reads the staging
container (CPU, memory, net and block IO) and the sra_all container (CPU only: the Logan load), and
/proc/loadavg the host load; one JSON line per sample goes to --out. Stops on SIGTERM / SIGINT / SIGHUP, when
--stop-file exists, after --max-seconds, when the --parent-pid process is gone, or when the --heartbeat file
was not touched for --heartbeat-max-age seconds (the orchestrator on the laptop touches it per scenario, so a
sampler outlives a dead orchestrator by at most that long). Stdlib only (python3.12 on mex).

    sampler.py --out RUN_DIR/samples.jsonl [--interval 2] [--label direct-local]
"""

import argparse
import json
import os
import re
import signal
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import suite_common as sc  # noqa: E402

CONTAINER = "metagraph_refseq33m_experimental"
SRA_CONTAINER = "metagraph_sra_all"

# docker prints memory in binary units and IO in decimal ones; both spellings are accepted for either.
_UNITS = {
    "b": 1, "kb": 1000, "mb": 1000 ** 2, "gb": 1000 ** 3, "tb": 1000 ** 4, "pb": 1000 ** 5,
    "kib": 1024, "mib": 1024 ** 2, "gib": 1024 ** 3, "tib": 1024 ** 4, "pib": 1024 ** 5,
}
_SIZE_RE = re.compile(r"^\s*([0-9]*\.?[0-9]+)\s*([A-Za-z]*)\s*$")


def parse_size(text):
    """'11.87GiB' -> bytes (int); '141kB' -> 141000; None for '--' or anything unparsable."""
    if text is None:
        return None
    m = _SIZE_RE.match(str(text))
    if not m:
        return None
    unit = m.group(2).lower() or "b"
    if unit not in _UNITS:
        return None
    return int(round(float(m.group(1)) * _UNITS[unit]))


def parse_pair(text):
    """'11.87GiB / 128GiB' -> (bytes, bytes); (None, None) when unparsable."""
    if not text or "/" not in str(text):
        return None, None
    left, right = str(text).split("/", 1)
    return parse_size(left), parse_size(right)


def parse_percent(text):
    """'7704.04%' -> 7704.04; None for '--' and for a figure that is not finite ('NaN%', 'inf%': no reading,
    suite_common.finite_or_none)."""
    if text is None:
        return None
    t = str(text).strip().rstrip("%").strip()
    try:
        return sc.finite_or_none(float(t))
    except ValueError:
        return None


def parse_stats_lines(output):
    """{container name: parsed stats} from the `docker stats --format '{{json .}}'` output."""
    out = {}
    for line in (output or "").splitlines():
        line = line.strip()
        if not line.startswith("{"):
            continue
        try:
            row = json.loads(line)
        except ValueError:
            continue
        name = row.get("Name") or row.get("Container")
        mem, mem_limit = parse_pair(row.get("MemUsage"))
        net_rx, net_tx = parse_pair(row.get("NetIO"))
        blk_read, blk_write = parse_pair(row.get("BlockIO"))
        out[name] = {"cpu_pct": parse_percent(row.get("CPUPerc")), "mem_bytes": mem, "mem_limit_bytes": mem_limit,
                     "net_rx": net_rx, "net_tx": net_tx, "blk_read": blk_read, "blk_write": blk_write}
    return out


def build_sample(t_utc, stats, load1, container, sra_container, label=None, error=None):
    """One samples.jsonl line (§8.4), null fields for a container docker did not report."""
    c = stats.get(container) or {}
    sra = stats.get(sra_container) or {}
    sample = {"t_utc": t_utc, "container": container,
              "cpu_pct": c.get("cpu_pct"), "mem_bytes": c.get("mem_bytes"),
              "mem_limit_bytes": c.get("mem_limit_bytes"), "net_rx": c.get("net_rx"), "net_tx": c.get("net_tx"),
              "blk_read": c.get("blk_read"), "blk_write": c.get("blk_write"),
              "host_load1": load1, "sra_all_cpu_pct": sra.get("cpu_pct")}
    if label:
        sample["phase"] = label
    if error:
        sample["error"] = error
    return sample


def read_load1():
    try:
        with open("/proc/loadavg") as fh:
            return float(fh.read().split()[0])
    except (OSError, ValueError, IndexError):
        try:
            return os.getloadavg()[0]
        except (OSError, AttributeError):
            return None


def docker_stats(containers, call_timeout):
    """(stdout, error) of one read-only `docker stats --no-stream` call."""
    cmd = ["docker", "stats", "--no-stream", "--format", "{{json .}}"] + list(containers)
    try:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=call_timeout,
                              universal_newlines=True)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return "", "%s: %s" % (type(exc).__name__, exc)
    err = proc.stderr.strip() if proc.returncode != 0 else None
    return proc.stdout, err


def process_alive(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def heartbeat_stale(path, max_age_s, now=None):
    """Whether the heartbeat file is missing or was last touched more than max_age_s ago."""
    try:
        mtime = os.stat(path).st_mtime
    except OSError:
        return True
    return ((now if now is not None else time.time()) - mtime) > max_age_s


_stop = {"flag": False}


def _on_signal(signum, frame):
    _stop["flag"] = True


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", required=True, help="samples.jsonl to append to")
    ap.add_argument("--interval", type=float, default=2.0, help="seconds between sample starts (2)")
    ap.add_argument("--container", default=CONTAINER)
    ap.add_argument("--sra-container", default=SRA_CONTAINER)
    ap.add_argument("--label", default=None, help="a `phase` field on every sample (e.g. the arm)")
    ap.add_argument("--stop-file", default=None, help="stop when this file exists")
    ap.add_argument("--max-seconds", type=float, default=0, help="stop after this long (0: no limit)")
    ap.add_argument("--pid-file", default=None, help="write the sampler's pid here")
    ap.add_argument("--parent-pid", type=int, default=0,
                    help="stop when this process is gone (the runner that started the sampler)")
    ap.add_argument("--heartbeat", default=None, help="stop when this file is older than --heartbeat-max-age")
    ap.add_argument("--heartbeat-max-age", type=float, default=3600.0)
    ap.add_argument("--once", action="store_true", help="take one sample and exit")
    args = ap.parse_args(argv)

    for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
        signal.signal(sig, _on_signal)
    if args.pid_file:
        with open(args.pid_file, "w") as fh:
            fh.write("%d\n" % os.getpid())
    out_dir = os.path.dirname(os.path.abspath(args.out))
    os.makedirs(out_dir, exist_ok=True)

    t_begin = time.time()
    taken = 0
    while not _stop["flag"]:
        t0 = time.time()
        stdout, err = docker_stats([args.container, args.sra_container], call_timeout=max(10.0, 5 * args.interval))
        t1 = time.time()
        stats = parse_stats_lines(stdout)
        # docker averages CPU over its own measuring window inside the call: the midpoint stamps it best.
        sample = build_sample(sc.iso_utc((t0 + t1) / 2.0), stats, read_load1(), args.container,
                              args.sra_container, args.label, err)
        sc.append_jsonl(args.out, sample)
        taken += 1
        if args.once:
            break
        if args.max_seconds and time.time() - t_begin >= args.max_seconds:
            break
        if args.stop_file and os.path.exists(args.stop_file):
            break
        if args.parent_pid and not process_alive(args.parent_pid):
            break
        if args.heartbeat and heartbeat_stale(args.heartbeat, args.heartbeat_max_age):
            break
        while not _stop["flag"] and time.time() < t0 + args.interval:
            time.sleep(min(0.2, max(0.0, t0 + args.interval - time.time())))
    if args.pid_file:
        try:
            os.remove(args.pid_file)
        except OSError:
            pass
    return 0


if __name__ == "__main__":
    sys.exit(main())
