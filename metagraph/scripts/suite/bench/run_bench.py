#!/usr/bin/env python3
"""Runs the isolated benchmark panel against the staging host (lane B): one unmeasured warm-up per query, then
six measured cycles of all three arms, every arm order once per query, a recorded seed, one request at a time.

    run_bench.py --panel PANEL [--seed N] [--arms direct-local,direct-remote,mcp] [--cycles 6]
                 [--warmup-arm direct-local] [--mcp-cmd TEMPLATE] [--pause 1] [--checks auto|checks|stub]
                 [--sra-max PCT] [--load1-max L] [--allow-missing-sra] [--preflight-s 20] [--no-preflight]
                 [--quiet-wait-s 1800] [--continue-loaded] [--max-steps N] [--no-analysis]
    run_bench.py --resume RUN_DIR [--max-steps N] [--no-preflight] [--quiet-wait-s S] [--continue-loaded]
                 [--no-analysis]
    run_bench.py --panel PANEL --dry-run [--seed N]

THE RUN (bench_lib.schedule): block 0 is the warm-up pass, one request per query through --warmup-arm (direct-local:
the backend), queries in a seeded order; blocks 1..6 are the measured cycles: every query once, in a fresh seeded
order, its three arms back to back in the order a seeded Latin square assigns (each query runs each of the six arm
orders once; each cycle runs each order 8 times for 48 queries). 48 x (1 + 6 x 3) = 912 requests. The run's manifest
(bench_lib.run_manifest) gives each query one scenario per phase (QUERY~w0, QUERY~c1..c6), so every step executes
exactly one request with the suite's own runners, unchanged:
  direct-local   run_direct.py --only BENCH_ID on mex over this invocation's ssh connection (the run's snapshot);
  direct-remote  run_direct.py --only BENCH_ID from here through the public URL (the run's snapshot);
  mcp            the MCP runner template (run_suite.sh's, SUITE_MCP_CMD or --mcp-cmd) with {ids} BENCH_ID,
                 {repeat} 0, {repeats} 1, {out} RUN_DIR/arms/mcp, {manifest} RUN_DIR/manifest.json and {checks_dir}
                 RUN_DIR/scripts (the run's snapshot, as run_suite.sh passes it), run with `bash -c`; only the
                 template is recorded. A cancelled MCP job is waited for (suite_common.mcp_step_result) before the
                 next step, so one request is in flight.
--pause seconds separate two steps. After each block one direct-local GET /capabilities (run_direct.py --finish)
records the host's release and server instance at the block's end (analyze.py reads a change as a new deployment
inside the run); the run ends with each direct arm's closing GET.

ISOLATION. Like run_suite.sh it first takes the suite lock ~/suite-dev/LOCK on mex and refuses to start next to a
running run_direct.py (a suite run and a bench run never overlap), then names the run from one GET /capabilities.
The suite's sampler runs on mex for the whole session (SUITE_SAMPLER=0: none; the run is then unguarded and
analyze.py emits no verdict). The preflight samples --preflight-s seconds and refuses to start (exit 4) unless the
load guard (load_guard.py) finds the host clean; before every later block the last minute's samples are judged and a
loaded host is waited for, up to --quiet-wait-s, then the run stops (exit 4, resumable); at each block's start and
every 12 steps within a block the run stops (exit 4) when the block turned loaded, the sampler stopped, or two
checks in a row found no reading of the load: such a block can give no verdict (--continue-loaded runs on for
descriptive numbers). The thresholds
(--sra-max 100 %, --load1-max 16, --allow-missing-sra) are recorded in bench.json, which the analysis uses.

RESUMABLE: --resume RUN_DIR continues the recorded schedule (never a new seed) with the run's snapshot and options;
a step whose execution is recorded is skipped. Ctrl-C (or SIGTERM) stops after the step in flight; a second one
interrupts the step (its runner records its request in flight). --max-steps N stops this invocation after N steps
(exit 3). The run directory, BENCH_RUNS (default bench/runs)/<UTC>_<release-short>-bench: bench.json (identity,
options, sessions), panel.json, manifest.json, schedule.json, steps.jsonl (one line per step: times, exit, outcome),
scripts/ (the snapshot), arms/<arm>/ (each runner's own files), samples/; at the end analyze.py writes analysis.json
and analysis.md (offline; --no-analysis skips it). Exit: 0 every step recorded, 1 a failed or incomplete run,
2 a setup error, 3 --max-steps, 4 refused or stopped on load, 130 interrupted.

Environment as run_suite.sh: SUITE_MEX (srv-metagraph@mex.inf.ethz.ch; local:DIR runs the mex side in DIR, for
tests), SUITE_MEX_DIR (suite-dev), SUITE_MEX_PY (python3), SUITE_LOCAL_BASE, SUITE_REMOTE_BASE, PYTHON,
SUITE_SAMPLER, SUITE_WAIT_RUNNER_S (150), SUITE_SSH_CONTROL_DIR, SUITE_HOST_ID, SUITE_RUNNER_PATTERN,
SUITE_MCP_CMD; BENCH_RUNS. The MCP template with --check-queue needs METAGRAPH_ADMIN_KEY_FILE (or
METAGRAPH_ADMIN_KEY). Stdlib only.
"""

import argparse
import json
import os
import random
import re
import shlex
import signal
import subprocess
import sys
import threading
import time
import urllib.request
import uuid

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import bench_lib as bl  # noqa: E402
import load_guard as lg  # noqa: E402

sc = bl.sc

# run_suite.sh's ASYNC_MCP_CMD (the tests keep the two equal)
DEFAULT_MCP_CMD = ("/Users/raetsch/git/services/metagraph-search-service/.venv/bin/python "
                   "/Users/raetsch/tmp/mgss/wt-suite/scripts/bench/suite_mcp_runner.py --manifest {manifest} "
                   "--ids {ids} --repeat {repeat} --out {out} --check-queue --queue-ignore sra-logan-chunks "
                   "--checks-dir {checks_dir}")
STATUS_PREFIX = "RUN_DIRECT_STATUS "
PROGRESS = re.compile(r"^\[(?P<arm>[\w-]+)\] \d+/\d+ (?P<sid>\S+) r(?P<rep>\d+) (?P<outcome>\S+) http=(?P<http>\S+) "
                      r"(?P<ms>[0-9.]+) ms")
EXIT_DONE, EXIT_FAILED, EXIT_SETUP, EXIT_CAPPED, EXIT_LOADED, EXIT_INTERRUPTED = 0, 1, 2, 3, 4, 130
RECENT_S = 60.0          # the samples judged before a block and during one
CHECK_EVERY_STEPS = 12   # how often a running block's samples are judged
QUIET_POLL_S = 30.0      # how often a loaded host is looked at again
HEARTBEAT_EVERY_S = 300  # the sampler's heartbeat (it stops 3600 s after the last)
SAMPLER_MAX_S = 6 * 3600


class Setup(Exception):
    """A setup refusal: exit 2 (or the code it names) with its message."""

    def __init__(self, text, code=EXIT_SETUP):
        Exception.__init__(self, text)
        self.code = code


def env(name, default):
    v = os.environ.get(name)
    return default if v is None else v


class Config(object):
    def __init__(self):
        self.mex = env("SUITE_MEX", "srv-metagraph@mex.inf.ethz.ch")
        self.mex_dir = env("SUITE_MEX_DIR", "suite-dev")
        self.mex_py = env("SUITE_MEX_PY", "python3")
        self.local_base = env("SUITE_LOCAL_BASE", "http://localhost:42010")
        self.remote_base = env("SUITE_REMOTE_BASE", "https://metagraph.ethz.ch:8080/api/refseq33m-experimental")
        self.py = env("PYTHON", sys.executable or "python3")
        self.sampler = env("SUITE_SAMPLER", "1") == "1"
        self.wait_runner_s = float(env("SUITE_WAIT_RUNNER_S", "150"))
        self.ctl_dir = env("SUITE_SSH_CONTROL_DIR", "/tmp/suite-ssh-%d" % os.getuid())
        self.runner_pattern = env("SUITE_RUNNER_PATTERN", r"[s]cripts/run_direct\.py")
        self.mcp_cmd = env("SUITE_MCP_CMD", DEFAULT_MCP_CMD)
        self.runs = env("BENCH_RUNS", os.path.join(HERE, "runs"))


# ---------------------------------------------------------------------------------------------------------------
# mex: ssh (one control master of this invocation) or, for tests, a local directory

class Mex(object):
    def __init__(self, cfg):
        self.target = cfg.mex
        self.local = cfg.mex[len("local:"):] if cfg.mex.startswith("local:") else None
        self.opts = ["-o", "BatchMode=yes", "-o", "ConnectTimeout=20", "-o", "ServerAliveInterval=30",
                     "-o", "LogLevel=ERROR"]
        self.ctl = None
        if not self.local:
            os.makedirs(cfg.ctl_dir, exist_ok=True)
            os.chmod(cfg.ctl_dir, 0o700)
            self.ctl = os.path.join(cfg.ctl_dir, "%C-" + str(os.getpid()))
            self.opts += ["-o", "ControlMaster=auto", "-o", "ControlPath=" + self.ctl, "-o", "ControlPersist=300"]

    def argv(self, cmd):
        if self.local:
            return ["bash", "-c", cmd]
        return ["ssh"] + self.opts + [self.target, cmd]

    def run(self, cmd):
        """(exit, stdout) of a shell command in mex's home directory (stderr dropped, no stdin)."""
        try:
            p = subprocess.run(self.argv(cmd), cwd=self.local, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                               stderr=subprocess.DEVNULL, universal_newlines=True)
        except OSError as exc:
            return 255, str(exc)
        return p.returncode, p.stdout

    def popen(self, cmd, stderr):
        return subprocess.Popen(self.argv(cmd), cwd=self.local, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                                stderr=stderr, universal_newlines=True, start_new_session=True)

    def path(self, rel):
        return os.path.join(self.local, rel) if self.local else "%s:%s" % (self.target, rel)

    def rsync(self, src, dst, delete=False, excludes=()):
        cmd = ["rsync", "-a"] + (["--delete"] if delete else []) + ["--exclude=%s" % e for e in excludes]
        if not self.local:
            cmd += ["-e", "ssh " + " ".join(shlex.quote(o) for o in self.opts)]
        try:
            return subprocess.run(cmd + [src, dst], stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                                  stderr=subprocess.DEVNULL).returncode
        except OSError:
            return 255

    def close(self):
        if self.ctl:
            subprocess.run(["ssh", "-o", "ControlPath=" + self.ctl, "-O", "exit", self.target],
                           stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


# ---------------------------------------------------------------------------------------------------------------
# signals: the first stops after the step in flight, the second interrupts it

class Stop(object):
    def __init__(self):
        self.count = 0
        self.proc = None

    def handler(self, signum, frame):
        self.count += 1
        if self.count >= 2 and self.proc is not None and self.proc.poll() is None:
            try:
                os.killpg(self.proc.pid, signal.SIGINT)
            except OSError:
                pass

    @property
    def requested(self):
        return self.count > 0

    @property
    def hard(self):
        return self.count > 1


STOP = Stop()


def log(text):
    sys.stderr.write("[run_bench %s] %s\n" % (time.strftime("%H:%M:%S", time.gmtime()), text))
    sys.stderr.flush()


def scripts_commit():
    try:
        head = subprocess.run(["git", "-C", SUITE, "rev-parse", "HEAD"], stdout=subprocess.PIPE,
                              stderr=subprocess.DEVNULL, universal_newlines=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", SUITE, "status", "--porcelain", "--", "."], stdout=subprocess.PIPE,
                               stderr=subprocess.DEVNULL, universal_newlines=True).stdout.splitlines()
    except OSError:
        return "unknown"
    dirty = [x for x in dirty if not re.match(r"^\?\? (metagraph/scripts/suite/)?(bench/)?runs/", x)]
    return (head or "unknown") + ("+dirty" if dirty else "")


def expand_template(template, values):
    """The MCP template with each {name} replaced by its shell-quoted value (run_suite.sh's mcp_step)."""
    for key, value in values.items():
        template = template.replace("{%s}" % key, shlex.quote(str(value)))
    return template


def check_mcp(template, arms):
    """Refuses an MCP arm whose template cannot run (before any request): run_suite.sh's rules."""
    if "mcp" not in arms:
        return
    tmpl = template or ""
    if "{ids}" not in tmpl or "{out}" not in tmpl:
        raise Setup("the MCP template needs {ids} and {out} (SUITE_MCP_CMD / --mcp-cmd; --arms without mcp runs "
                    "without the MCP arm)")
    exe = tmpl.split(" ", 1)[0]
    if exe.startswith("/") and not os.access(exe, os.X_OK):
        raise Setup("the MCP runner's interpreter %s is missing (--mcp-cmd names another)" % exe)
    if "--check-queue" in tmpl and not (os.environ.get("METAGRAPH_ADMIN_KEY_FILE")
                                        or os.environ.get("METAGRAPH_ADMIN_KEY")):
        raise Setup("the MCP runner checks the service queue and needs METAGRAPH_ADMIN_KEY_FILE (or "
                    "METAGRAPH_ADMIN_KEY) in the environment")
    if "{checks_dir}" not in tmpl:
        log("warning: the MCP template names no {checks_dir}: its runner judges with other checks than the run's "
            "snapshot")


class Bench(object):
    def __init__(self, args, cfg):
        self.args, self.cfg = args, cfg
        self.mex = Mex(cfg)
        self.lock_held = False
        self.lock_token = uuid.uuid4().hex
        self.lock_host = sc.host_id()
        self.lock_started = time.time()
        self.run_id = None
        self.run_dir = None
        self.sampler_pid = None
        self.sampler_rel = None
        self.hb_stop = threading.Event()
        self.hb_thread = None
        self.session = None
        self.used_requests = 0
        self.ran_arms = set()
        self.schedule = []
        self.unguarded_checks = 0

    # -- the lock (run_suite.sh's protocol and document) ------------------------------------------------------
    def lock_json(self):
        return json.dumps({"pid": os.getpid(), "host": self.lock_host, "hostname": os.uname()[1],
                           "started_utc": sc.iso_utc(self.lock_started), "started_epoch": round(self.lock_started, 3),
                           "run_id": self.run_id or "", "mode": "bench", "token": self.lock_token})

    def take_lock(self):
        d = self.cfg.mex_dir
        rc, _ = self.mex.run("mkdir -p %s" % d)
        if rc != 0:
            raise Setup("ssh %s failed" % self.cfg.mex)
        for _ in range(3):
            rc, _ = self.mex.run("(set -o noclobber; printf '%%s\\n' %s > %s/LOCK)" % (
                shlex.quote(self.lock_json()), d))
            if rc == 0:
                self.lock_held = True
                log("lock ~/%s/LOCK taken on %s" % (d, self.cfg.mex))
                return
            _, held = self.mex.run("cat %s/LOCK" % d)
            held = held.strip()
            if not held:
                continue
            try:
                doc = json.loads(held)
            except ValueError:
                doc = None
            verdict, reason = sc.lock_verdict(doc, self.lock_host, time.time(), sc._pid_alive, sc._pid_age_s)
            if verdict != "stale":
                raise Setup("the suite lock ~/%s/LOCK on %s is %s %s; another invocation drives the host" % (
                    d, self.cfg.mex, verdict, reason))
            log("removing a stale lock (%s)" % reason)
            if self.mex.run(sc.lock_break_command(d, held))[0] != 0:
                time.sleep(1)
        raise Setup("could not take the lock ~/%s/LOCK on %s" % (d, self.cfg.mex))

    def lock_update(self):
        d = self.cfg.mex_dir
        self.mex.run("grep -q %s %s/LOCK && printf '%%s\\n' %s > %s/LOCK.%s && mv -f %s/LOCK.%s %s/LOCK" % (
            self.lock_token, d, shlex.quote(self.lock_json()), d, self.lock_token, d, self.lock_token, d))

    def release_lock(self):
        if not self.lock_held:
            return
        d = self.cfg.mex_dir
        self.mex.run("if grep -q %s %s/LOCK 2>/dev/null; then rm -f %s/LOCK; fi" % (self.lock_token, d, d))
        self.lock_held = False

    def runner_pids(self):
        """The suite's run_direct.py processes on mex ('' none); None when mex did not answer."""
        rc, out = self.mex.run("printf 'runners:'; pgrep -u $(id -u) -f %s | tr '\\n' ' '; echo"
                               % shlex.quote(self.cfg.runner_pattern))
        if "runners:" not in out:
            return None
        return out.split("runners:", 1)[1].strip()

    def stop_mex_runner(self):
        """SIGINT to this run's runner on mex (it stops after its request in flight), then up to 120 s for it."""
        pattern = "%s .*--run-id %s" % (self.cfg.runner_pattern, self.run_id)
        self.mex.run("pkill -INT -u $(id -u) -f %s" % shlex.quote(pattern))
        if not self.wait_no_runner(120):
            log("the direct-local runner still runs on %s; --resume waits for it" % self.cfg.mex)

    def wait_no_runner(self, max_s):
        waited = 0.0
        while True:
            pids = self.runner_pids()
            if pids == "":
                return True
            if waited >= max_s:
                return False
            time.sleep(5)
            waited += 5

    # -- the run directory --------------------------------------------------------------------------------------
    def release_of_host(self):
        try:
            with urllib.request.urlopen(self.cfg.remote_base.rstrip("/") + "/capabilities", timeout=30) as r:
                return json.loads(r.read().decode("utf-8")).get("release") or ""
        except (OSError, ValueError):
            return ""

    def new_run(self):
        panel = bl.read_panel(self.args.panel)
        arms = [a for a in self.args.arms.split(",") if a]
        if any(a not in bl.ARMS for a in arms) or len(set(arms)) != len(arms) or not arms:
            raise Setup("--arms takes distinct arms of %s" % ",".join(bl.ARMS))
        if self.args.warmup_arm not in bl.ARMS:
            raise Setup("--warmup-arm is one of %s" % ",".join(bl.ARMS))
        seed = self.args.seed if self.args.seed is not None else random.SystemRandom().randrange(1, 2 ** 31)
        self.take_lock()
        if not self.wait_no_runner(self.cfg.wait_runner_s):
            raise Setup("a run_direct.py runs under ~/%s on %s; wait for it" % (self.cfg.mex_dir, self.cfg.mex))
        log("release: one GET /capabilities on %s (names the run directory)" % self.cfg.remote_base)
        release = self.release_of_host()
        self.run_id = "%s_%s-bench" % (sc.run_stamp(time.time()), sc.release_short(release))
        self.run_dir = os.path.join(self.cfg.runs, self.run_id)
        if os.path.exists(self.run_dir):
            raise Setup("%s exists already" % self.run_dir)
        self.lock_update()
        os.makedirs(os.path.join(self.run_dir, "arms"))
        sc.write_json_atomic(os.path.join(self.run_dir, "panel.json"), panel)
        manifest = bl.run_manifest(panel, self.args.cycles)
        sc.write_json_atomic(os.path.join(self.run_dir, "manifest.json"), manifest)
        steps = bl.schedule([q["id"] for q in panel["queries"]], seed, self.args.cycles, arms, self.args.warmup_arm)
        sc.write_json_atomic(os.path.join(self.run_dir, "schedule.json"), {
            "seed": seed, "cycles": self.args.cycles, "arms": arms, "warmup_arm": self.args.warmup_arm,
            "arm_orders": bl.arm_orders(arms), "balance": bl.schedule_balance(steps, arms), "steps": steps})
        commit = scripts_commit()
        sc.take_snapshot(SUITE, self.run_dir, commit)
        bench = {"run_id": self.run_id, "created_utc": sc.iso_utc(time.time()), "release_at_start": release or None,
                 "panel": {"path": os.path.abspath(self.args.panel), "name": panel.get("name"),
                           "panel_sha256": panel["panel_sha256"], "queries": len(panel["queries"])},
                 "manifest_sha256": manifest["content_sha256"], "seed": seed, "cycles": self.args.cycles,
                 "arms": arms, "warmup_arm": self.args.warmup_arm, "pause_s": self.args.pause,
                 "checks": self.args.checks, "allow_stub": self.args.allow_stub,
                 "mcp_cmd": self.cfg.mcp_cmd if "mcp" in arms else None,
                 "local_base": self.cfg.local_base, "remote_base": self.cfg.remote_base, "mex": self.cfg.mex,
                 "scripts_commit": commit, "requests_planned": len(steps),
                 "guard": lg.thresholds({"sra_all_cpu_pct_max": self.args.sra_max,
                                         "host_load1_max": self.args.load1_max,
                                         "allow_missing_sra": self.args.allow_missing_sra or None}),
                 "sampler": self.cfg.sampler, "sessions": [], "state": "new"}
        sc.write_json_atomic(self.bench_path(), bench)
        log("run directory %s (seed %d, %d steps)" % (self.run_dir, seed, len(steps)))
        return bench

    def resume_run(self):
        self.run_dir = os.path.abspath(self.args.resume)
        bench = sc.read_json(self.bench_path())
        if not isinstance(bench, dict):
            raise Setup("--resume %s holds no bench.json" % self.run_dir)
        self.run_id = bench["run_id"]
        self.take_lock()
        self.lock_update()
        if not self.wait_no_runner(self.cfg.wait_runner_s):
            raise Setup("a run_direct.py still runs under ~/%s on %s (this run's cut-off step?); wait for it, then "
                        "--resume again" % (self.cfg.mex_dir, self.cfg.mex))
        doc, changed, live = sc.snapshot_state(self.run_dir, SUITE)
        if doc is None:
            unreadable = sc.snapshot_refusal(sc.snapshot_dir(self.run_dir), ())
            raise Setup("%s: the run's snapshot record is unreadable: %s" % (self.run_dir, unreadable[0]) if unreadable
                        else "%s has no snapshot of its scripts" % self.run_dir)
        if changed:
            log("warning: the run's snapshot changed after it was taken (%s)" % ", ".join(changed))
        if live:
            log("scripts: %s differ from the run's snapshot: the run keeps judging with its snapshot" % ", ".join(live))
        log("resuming %s (seed %s)" % (self.run_id, bench.get("seed")))
        return bench

    def bench_path(self):
        return os.path.join(self.run_dir, "bench.json")

    def save_session(self, bench, **fields):
        self.session.update(fields)
        doc = sc.read_json(self.bench_path()) or bench
        sessions = doc.get("sessions") or []
        if sessions and sessions[-1].get("id") == self.session["id"]:
            sessions[-1] = self.session
        else:
            sessions.append(self.session)
        doc["sessions"] = sessions
        if "state" in fields:
            doc["state"] = fields["state"]
        sc.write_json_atomic(self.bench_path(), doc)

    # -- mex: scripts, clock, sampler -----------------------------------------------------------------------------
    def arm_rel(self):
        return "runs/%s/arms/direct-local" % self.run_id

    def prepare_mex(self):
        d = self.cfg.mex_dir
        if self.mex.run("mkdir -p %s/scripts %s/manifest %s/%s %s/runs/%s/samples" % (
                d, d, d, self.arm_rel(), d, self.run_id))[0] != 0:
            raise Setup("ssh %s failed" % self.cfg.mex)
        if self.mex.rsync(os.path.join(self.run_dir, "scripts") + "/", self.mex.path(d + "/scripts/"), delete=True,
                          excludes=("__pycache__/", "*.pyc")) != 0:
            raise Setup("rsync of the run's snapshot to %s failed" % self.cfg.mex)
        if self.mex.rsync(os.path.join(self.run_dir, "manifest.json"),
                          self.mex.path("%s/manifest/%s.json" % (d, self.run_id))) != 0:
            raise Setup("rsync of the run's manifest to %s failed" % self.cfg.mex)
        local = os.path.join(self.run_dir, "arms", "direct-local")
        if os.path.isdir(local) and self.mex.run("test -f %s/%s/run.json" % (d, self.arm_rel()))[0] != 0:
            log("the mex copy of direct-local is gone: pushing the local one back")
            self.mex.rsync(local + "/", self.mex.path("%s/%s/" % (d, self.arm_rel())))
        else:
            self.pull_local_arm()

    def pull_local_arm(self):
        d = self.cfg.mex_dir
        if self.mex.run("test -d %s/%s" % (d, self.arm_rel()))[0] == 0:
            os.makedirs(os.path.join(self.run_dir, "arms", "direct-local"), exist_ok=True)
            if self.mex.rsync(self.mex.path("%s/%s/" % (d, self.arm_rel())),
                              os.path.join(self.run_dir, "arms", "direct-local") + "/") != 0:
                log("pulling direct-local from %s failed" % self.cfg.mex)

    def pull_samples(self):
        d = self.cfg.mex_dir
        if self.mex.run("test -d %s/runs/%s/samples" % (d, self.run_id))[0] == 0:
            os.makedirs(os.path.join(self.run_dir, "samples"), exist_ok=True)
            self.mex.rsync(self.mex.path("%s/runs/%s/samples/" % (d, self.run_id)),
                           os.path.join(self.run_dir, "samples") + "/", excludes=(".heartbeat",))

    def clock_offset(self):
        t_a = time.time()
        _, out = self.mex.run("t=$(date +%%s.%%N 2>/dev/null); case \"$t\" in ''|*[!0-9.]*) %s -c "
                              "'import time; print(\"%%.3f\" %% time.time())' ;; *) echo \"$t\" ;; esac"
                              % shlex.quote(self.cfg.mex_py))
        t_b = time.time()
        try:
            return round(float(out.strip()) - (t_a + t_b) / 2.0, 3)
        except ValueError:
            return None

    def start_sampler(self):
        if not self.cfg.sampler:
            log("no sampler (SUITE_SAMPLER=0): the run is unguarded and analyze.py emits no performance verdict")
            return
        d = self.cfg.mex_dir
        self.sampler_rel = "runs/%s/samples/bench-%s.jsonl" % (self.run_id, sc.run_stamp(time.time()))
        hb = "runs/%s/samples/.heartbeat" % self.run_id
        cmd = " ".join(shlex.quote(x) for x in [self.cfg.mex_py, "scripts/sampler.py", "--out", self.sampler_rel,
                                                "--label", "bench", "--max-seconds", str(SAMPLER_MAX_S),
                                                "--heartbeat", hb, "--heartbeat-max-age", "3600"])
        _, out = self.mex.run("cd %s && touch %s && (nohup %s > /dev/null 2>&1 < /dev/null & echo $!)" % (d, hb, cmd))
        pid = out.strip()
        if re.match(r"^[0-9]+$", pid):
            self.sampler_pid = pid
            log("sampler on %s (pid %s)" % (self.cfg.mex, pid))
            self.hb_stop.clear()
            self.hb_thread = threading.Thread(target=self.heartbeat_loop, args=(d + "/" + hb,), daemon=True)
            self.hb_thread.start()
        else:
            log("no sampler on %s (got %r): this session has no resource samples" % (self.cfg.mex, pid))

    def heartbeat_loop(self, path):
        while not self.hb_stop.wait(HEARTBEAT_EVERY_S):
            self.mex.run("touch %s" % path)

    def stop_sampler(self):
        self.hb_stop.set()
        if self.sampler_pid:
            self.mex.run("kill -TERM %s" % self.sampler_pid)
            for _ in range(6):
                if self.mex.run("kill -0 %s" % self.sampler_pid)[0] != 0:
                    break
                time.sleep(1)
            self.sampler_pid = None

    def recent_samples(self, seconds):
        """The sampler's samples of the last |seconds| on mex's clock (read from its file on mex); none when its newest
        is older than that."""
        if not self.sampler_rel:
            return []
        _, out = self.mex.run("tail -n %d %s/%s 2>/dev/null" % (int(seconds) + 10, self.cfg.mex_dir, self.sampler_rel))
        rows = []
        for line in out.splitlines():
            try:
                rows.append(json.loads(line))
            except ValueError:
                pass
        times = [sc.parse_utc(r.get("t_utc")) for r in rows]
        newest = max([t for t in times if t is not None] or [0])
        # samples older than that are no reading of the host now (a sampler that stopped): none
        if newest < time.time() + (self.session.get("clock_offset_start_s") or 0.0) - seconds:
            return []
        return [r for r, t in zip(rows, times) if t is not None and t >= newest - seconds]

    def judge_now(self, bench, seconds):
        return lg.judge_samples(self.recent_samples(seconds), bench.get("guard"))

    def sampler_alive(self):
        return bool(self.sampler_pid) and self.mex.run("kill -0 %s" % self.sampler_pid)[0] == 0

    def live_check(self, bench, block):
        """Within the run (every CHECK_EVERY_STEPS steps and at each block's start): None to go on, else why the run
        stops. A block that turned loaded, a sampler that stopped, or two checks in a row without a reading of the
        host's load can give no performance verdict: stop rather than spend the reserved window (--continue-loaded
        runs on, its numbers descriptive)."""
        if not self.sampler_alive():
            why = "the sampler stopped (pid %s): the load of the rest of the run is unknown" % self.sampler_pid
        else:
            res = self.judge_now(bench, RECENT_S)
            if res["verdict"] == "unguarded":
                self.unguarded_checks += 1
            else:
                self.unguarded_checks = 0
            if res["verdict"] == "loaded":
                why = "loaded during block %d (%s)" % (block, res["reason"])
            elif self.unguarded_checks >= 2:
                why = "the load is unknown during block %d (%s, twice in a row)" % (block, res["reason"])
            else:
                return None
        if self.args.continue_loaded:
            log("%s: no performance verdict from this run (--continue-loaded runs on)" % why)
            return None
        log("%s: stopping (no performance verdict can come from it; --continue-loaded runs on)" % why)
        return why

    def preflight(self, bench):
        if self.args.no_preflight:
            log("preflight skipped (--no-preflight)")
            return {"skipped": "--no-preflight"}
        if not self.sampler_pid:
            log("preflight skipped: no sampler")
            return {"skipped": "no sampler"}
        log("preflight: sampling the host for %s s" % self.args.preflight_s)
        end = time.time() + self.args.preflight_s
        while time.time() < end and not STOP.requested:
            time.sleep(0.5)
        res = self.judge_now(bench, self.args.preflight_s + 5)
        log("preflight: %s (%s; %d samples, sra_all max %s %%, load1 max %s)" % (
            res["verdict"], res["reason"], res["samples"], res["sra_all_cpu_pct_max"], res["host_load1_max"]))
        if res["verdict"] != "clean":
            raise Setup("preflight: the host is %s (%s): no bench run now (--no-preflight overrides; the analysis "
                        "would emit no verdict)" % (res["verdict"], res["reason"]), EXIT_LOADED)
        return res

    def quiet_before(self, bench, block):
        """Waits while the last minute's samples are loaded, up to --quiet-wait-s; False when still loaded."""
        if not self.sampler_pid:
            return True
        waited = 0.0
        while True:
            res = self.judge_now(bench, RECENT_S)
            if res["verdict"] != "loaded":
                if res["verdict"] == "unguarded":
                    log("before block %d: the load is unknown (%s)" % (block, res["reason"]))
                return True
            if STOP.requested:
                return True  # the loop stops as interrupted
            if waited >= self.args.quiet_wait_s:
                log("before block %d: the host is still loaded after %d s (%s): stopping" % (
                    block, waited, res["reason"]))
                return False
            if waited == 0:
                log("before block %d: the host is loaded (%s): waiting up to %s s" % (
                    block, res["reason"], self.args.quiet_wait_s))
            time.sleep(QUIET_POLL_S)
            waited += QUIET_POLL_S

    # -- steps --------------------------------------------------------------------------------------------------
    def recorded(self):
        """{(bench id, arm)} of the executions every arm recorded (direct-local's local copy, pulled)."""
        out = set()
        for arm in bl.ARMS:
            rows, _ = sc.read_jsonl(os.path.join(self.run_dir, "arms", arm, "executions.jsonl"))
            out |= set((r.get("scenario"), r.get("arm")) for r in rows if r.get("arm") == arm)
        return out

    def direct_args(self, bench, arm, ids):
        args = ["--arm", arm, "--run-id", self.run_id, "--checks", bench.get("checks") or "auto", "--warm-every", "0",
                "--pause", "0", "--scripts-commit", bench.get("scripts_commit") or "unknown", "--ids", ",".join(ids)]
        if bench.get("allow_stub"):
            args.append("--allow-stub")
        return args

    def run_proc(self, proc):
        STOP.proc = proc
        try:
            out, _ = proc.communicate()
        finally:
            STOP.proc = None
        return proc.returncode, out or ""

    def direct_step(self, bench, step, ids, finish=False):
        arm = step["arm"] if step else "direct-local"
        what = ["--finish"] if finish else ["--only", step["bench_id"]]
        if arm == "direct-local":
            d = self.cfg.mex_dir
            argv = [self.cfg.mex_py, "scripts/run_direct.py", "--base", self.cfg.local_base, "--out", self.arm_rel(),
                    "--manifest", "manifest/%s.json" % self.run_id] + self.direct_args(bench, arm, ids) + what + \
                ["--status-line"]
            cmd = "cd %s && %s 2>> %s/runner.log" % (d, " ".join(shlex.quote(x) for x in argv), self.arm_rel())
            err_path = os.path.join(self.run_dir, ".mex-step.err")
            errf = open(err_path, "w")
            try:
                rc, out = self.run_proc(self.mex.popen(cmd, errf))
            finally:
                errf.close()
            if rc != 0:
                with open(err_path) as fh:
                    for line in fh.read().splitlines()[-5:]:
                        log("    ssh: %s" % line)
            os.remove(err_path)
        else:
            out_dir = os.path.join(self.run_dir, "arms", arm)
            os.makedirs(out_dir, exist_ok=True)
            argv = [self.cfg.py, os.path.join(self.run_dir, "scripts", "run_direct.py"), "--base",
                    self.cfg.remote_base if arm == "direct-remote" else self.cfg.local_base, "--out", out_dir,
                    "--manifest", os.path.join(self.run_dir, "manifest.json")] + \
                self.direct_args(bench, arm, ids) + what + ["--status-line"]
            with open(os.path.join(out_dir, "runner.log"), "a") as errf:
                proc = subprocess.Popen(argv, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=errf,
                                        universal_newlines=True, start_new_session=True)
                rc, out = self.run_proc(proc)
        status, outcome, ms = {}, None, None
        for line in out.splitlines():
            if line.startswith(STATUS_PREFIX):
                try:
                    status = json.loads(line[len(STATUS_PREFIX):])
                except ValueError:
                    status = {}
            m = PROGRESS.match(line)
            if m:
                outcome, ms = m.group("outcome"), float(m.group("ms"))
        self.used_requests += int(status.get("requests") or 0)
        return {"rc": rc, "executed": bool(status.get("executed")), "outcome": outcome, "wall_ms": ms,
                "requests": status.get("requests"), "state": status.get("state")}

    def mcp_step(self, bench, step):
        out_dir = os.path.join(self.run_dir, "arms", "mcp")
        os.makedirs(out_dir, exist_ok=True)
        cmd = expand_template(bench.get("mcp_cmd") or "", {
            "ids": step["bench_id"], "repeat": 0, "repeats": 1, "out": out_dir,
            "manifest": os.path.join(self.run_dir, "manifest.json"),
            "checks_dir": os.path.join(self.run_dir, "scripts")})
        with open(os.path.join(out_dir, "runner.log"), "a") as logf:
            logf.write("[run_bench %s] mcp: %s\n" % (sc.iso_utc(time.time()), step["bench_id"]))
            logf.flush()
            proc = subprocess.Popen(["bash", "-c", cmd], stdin=subprocess.DEVNULL, stdout=logf,
                                    stderr=subprocess.STDOUT, start_new_session=True)
            STOP.proc = proc
            try:
                rc = proc.wait()
            finally:
                STOP.proc = None
            logf.write("[run_bench %s] mcp: exit %s\n" % (sc.iso_utc(time.time()), rc))
        manifest = sc.read_json(os.path.join(self.run_dir, "manifest.json"))
        outcomes, wait_s = sc.mcp_step_result(out_dir, manifest, [step["bench_id"]], 0, time.time())
        if wait_s > 0:
            log("mcp: the cancelled job of %s may still run on the host: waiting up to %d s" % (
                step["bench_id"], int(wait_s)))
            while wait_s > 0 and not STOP.hard:
                time.sleep(min(10.0, wait_s))
                outcomes, wait_s = sc.mcp_step_result(out_dir, manifest, [step["bench_id"]], 0, time.time())
        outcome = outcomes[0] if outcomes else "missing"
        wall = None
        rows, _ = sc.read_jsonl(os.path.join(out_dir, "executions.jsonl"))
        for r in rows:
            if r.get("scenario") == step["bench_id"] and r.get("repeat") == 0:
                wall = r.get("wall_ms")
        return {"rc": rc, "executed": outcome != "missing", "outcome": outcome, "wall_ms": wall}

    def identity_probe(self, bench, block):
        """direct-local's GET /capabilities at a block's end (run_direct.py --finish): the host's release and server
        instance then, recorded in steps.jsonl (analyze.py reads a change as a deployment inside the run)."""
        res = self.direct_step(bench, None, self.arm_ids(bench, "direct-local"), finish=True)
        self.pull_local_arm()
        doc = sc.read_json(os.path.join(self.run_dir, "arms", "direct-local", "run.json")) or {}
        last = (doc.get("sessions") or [{}])[-1]
        self.append_step({"probe": "identity", "block": block, "arm": "direct-local", "rc": res["rc"],
                          "requests": res["requests"], "t_end": round(time.time(), 3),
                          "release": last.get("release_end"), "server_instance": last.get("server_instance_end"),
                          "index_fp": last.get("index_fp_end")})
        if res["rc"] != 0:
            log("block %d: the identity probe failed (exit %s)" % (block, res["rc"]))

    def arm_ids(self, bench, arm):
        sched = self.schedule
        return [s["bench_id"] for s in sched if s["arm"] == arm]

    def append_step(self, line):
        line = dict(line)
        line.setdefault("session", self.session["id"])
        sc.append_jsonl(os.path.join(self.run_dir, "steps.jsonl"), line)

    # -- the session ----------------------------------------------------------------------------------------------
    def run(self):
        if not self.args.resume:
            check_mcp(self.cfg.mcp_cmd, [a for a in self.args.arms.split(",") if a])
        bench = self.resume_run() if self.args.resume else self.new_run()
        self.schedule = (sc.read_json(os.path.join(self.run_dir, "schedule.json")) or {}).get("steps") or []
        self.session = {"id": len(bench.get("sessions") or []), "started_utc": sc.iso_utc(time.time()),
                        "host": os.uname()[1], "argv": [a for a in sys.argv[1:]], "state": "running"}
        arms = bench.get("arms") or list(bl.ARMS)
        if self.args.resume:
            check_mcp(bench.get("mcp_cmd"), arms)
        self.prepare_mex()
        self.save_session(bench, clock_offset_start_s=self.clock_offset())
        self.start_sampler()
        self.save_session(bench, sampler_file=self.sampler_rel, sampler_pid=self.sampler_pid)
        self.session["preflight"] = self.preflight(bench)
        state, reason = self.loop(bench, arms)
        return self.finish(bench, arms, state, reason)

    def loop(self, bench, arms):
        done = self.recorded()
        todo = [s for s in self.schedule if (s["bench_id"], s["arm"]) not in done]
        log("%d of %d steps to run (%d recorded)" % (len(todo), len(self.schedule), len(self.schedule) - len(todo)))
        fails = dict((a, 0) for a in bl.ARMS)
        mcp_bad = 0
        ids = dict((a, self.arm_ids(bench, a)) for a in bl.ARMS)
        block, probed, ran, state, reason = None, None, 0, "done", None
        total = len(self.schedule)
        for step in todo:
            if STOP.requested:
                state, reason = "interrupted", "signal"
                break
            if self.args.max_steps is not None and ran >= self.args.max_steps:
                state, reason = "capped", "max_steps %d" % self.args.max_steps
                break
            if step["block"] != block:
                if block is not None:
                    self.identity_probe(bench, block)
                    probed = block
                if not (block is None and self.session.get("preflight", {}).get("verdict") == "clean"):
                    if not self.quiet_before(bench, step["block"]):
                        state, reason = "stopped", "loaded before block %d" % step["block"]
                        break
                    why = self.live_check(bench, step["block"]) if self.sampler_pid and block is not None else None
                    if why:
                        state, reason = "stopped", why
                        break
                block = step["block"]
            t0 = time.time()
            if step["arm"] == "mcp":
                res = self.mcp_step(bench, step)
            else:
                res = self.direct_step(bench, step, ids[step["arm"]])
            t1 = time.time()
            ran += 1
            self.ran_arms.add(step["arm"])
            line = dict((k, step[k]) for k in ("step", "block", "phase", "cycle", "query", "bench_id", "arm",
                                               "position", "order", "measured"))
            line.update({"t_start": round(t0, 3), "t_end": round(t1, 3), "started_utc": sc.iso_utc(t0),
                         "ended_utc": sc.iso_utc(t1)}, **res)
            self.append_step(line)
            log("%d/%d block %d %s %s %s: %s%s" % (
                step["step"] + 1, total, step["block"], step["bench_id"], step["arm"],
                "" if step["order"] is None else "(order %d, pos %d)" % (step["order"], step["position"]),
                res.get("outcome") or ("exit %s" % res["rc"]),
                "" if res.get("wall_ms") is None else " %.0f ms" % res["wall_ms"]))
            if step["arm"] == "mcp":
                if res["rc"] in sc.MCP_START_REFUSALS:
                    log("error: %s (exit %d): stopping the bench" % (sc.MCP_START_REFUSALS[res["rc"]], res["rc"]))
                    state, reason = "failed", "mcp_exit_%d" % res["rc"]
                    break
                if res["rc"] == 130:
                    state, reason = "interrupted", "mcp_interrupted"
                    break
                mcp_bad = mcp_bad + 1 if res["outcome"] in ("error", "missing", "unknown") else 0
                if mcp_bad >= 3:
                    state, reason = "failed", "mcp_errors"
                    break
            else:
                if res["rc"] != 0 or not res["executed"]:
                    fails[step["arm"]] += 1
                    if step["arm"] == "direct-local":
                        if STOP.hard:
                            self.stop_mex_runner()
                        if not self.wait_no_runner(self.cfg.wait_runner_s):
                            state, reason = "failed", "runner_still_running"
                            break
                else:
                    fails[step["arm"]] = 0
                if max(fails.values()) >= 3:
                    state, reason = "failed", "failing_steps"
                    break
            if STOP.hard and step["arm"] == "direct-local":
                self.stop_mex_runner()
            if STOP.requested:
                state, reason = "interrupted", "signal"
                break
            if self.sampler_pid and ran % CHECK_EVERY_STEPS == 0:
                why = self.live_check(bench, block)
                if why:
                    state, reason = "stopped", why
                    break
            pause_end = time.time() + float(bench.get("pause_s") or 0)
            while time.time() < pause_end and not STOP.requested:
                time.sleep(min(0.1, pause_end - time.time()))
        if block is not None and block != probed and reason != "runner_still_running":
            self.identity_probe(bench, block)
        return state, reason

    def finish(self, bench, arms, state, reason):
        # the closing GET of direct-remote (direct-local's is the last block's identity probe); none while a
        # direct-local runner may still run on the host
        if "direct-remote" in self.ran_arms and reason != "runner_still_running":
            res = self.direct_step(bench, {"arm": "direct-remote"}, self.arm_ids(bench, "direct-remote"), finish=True)
            if res["rc"] != 0:
                log("direct-remote: the closing GET failed (exit %s)" % res["rc"])
        self.pull_local_arm()
        self.stop_sampler()
        self.pull_samples()
        recorded = self.recorded()
        missing = sum(1 for s in self.schedule if (s["bench_id"], s["arm"]) not in recorded)
        if state == "done" and missing:
            state, reason = "incomplete", "%d step(s) not recorded" % missing
        self.save_session(bench, ended_utc=sc.iso_utc(time.time()), clock_offset_end_s=self.clock_offset(),
                          state=state, reason=reason, steps_missing=missing, direct_requests=self.used_requests)
        log("%s%s: %d of %d steps recorded" % (
            state, "" if not reason else " (%s)" % reason, len(self.schedule) - missing, len(self.schedule)))
        return {"done": EXIT_DONE, "capped": EXIT_CAPPED, "stopped": EXIT_LOADED,
                "interrupted": EXIT_INTERRUPTED}.get(state, EXIT_FAILED)

    def cleanup(self):
        self.stop_sampler()
        self.release_lock()
        self.mex.close()


def dry_run(args):
    panel = bl.read_panel(args.panel)
    arms = [a for a in args.arms.split(",") if a]
    seed = args.seed if args.seed is not None else 1
    steps = bl.schedule([q["id"] for q in panel["queries"]], seed, args.cycles, arms, args.warmup_arm)
    bal = bl.schedule_balance(steps, arms)
    print("panel %s (%s): %d queries, %d steps (%s), seed %d; first arms %s; costing ~%s host-hours (allowance %s)"
          % (args.panel, panel.get("name"), len(panel["queries"]), len(steps),
             bl.requests_planned(len(panel["queries"]), args.cycles, arms), seed, bal["first_arm"],
             panel.get("costing", {}).get("host_hours_estimate"), panel.get("costing", {}).get("allowance_host_hours")))
    return 0


def parse_args(argv):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    what = ap.add_mutually_exclusive_group(required=True)
    what.add_argument("--panel")
    what.add_argument("--resume")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--arms", default=",".join(bl.ARMS))
    ap.add_argument("--cycles", type=int, default=bl.DEFAULT_CYCLES)
    ap.add_argument("--warmup-arm", default=bl.DEFAULT_WARMUP_ARM)
    ap.add_argument("--mcp-cmd", default=None)
    ap.add_argument("--pause", type=float, default=1.0)
    ap.add_argument("--checks", choices=("auto", "checks", "stub"), default="auto")
    ap.add_argument("--allow-stub", action="store_true")
    ap.add_argument("--sra-max", type=float, default=None)
    ap.add_argument("--load1-max", type=float, default=None)
    ap.add_argument("--allow-missing-sra", action="store_true")
    ap.add_argument("--preflight-s", type=float, default=20.0)
    ap.add_argument("--no-preflight", action="store_true")
    ap.add_argument("--quiet-wait-s", type=float, default=1800.0)
    ap.add_argument("--continue-loaded", action="store_true",
                    help="run on when a block turns loaded (its numbers are then descriptive only)")
    ap.add_argument("--max-steps", type=int, default=None)
    ap.add_argument("--no-analysis", action="store_true")
    args = ap.parse_args(argv)
    if args.cycles < 1:
        ap.error("--cycles must be 1 or more")
    return args


def main(argv=None):
    args = parse_args(argv)
    if args.dry_run:
        if not args.panel:
            raise SystemExit("--dry-run needs --panel")
        return dry_run(args)
    cfg = Config()
    if args.mcp_cmd is not None:
        cfg.mcp_cmd = args.mcp_cmd
    STOP.count, STOP.proc = 0, None
    for sig in (signal.SIGINT, signal.SIGTERM):
        signal.signal(sig, STOP.handler)
    bench = Bench(args, cfg)
    rc = EXIT_SETUP
    try:
        rc = bench.run()
    except Setup as exc:
        log("error: %s" % exc)
        rc = exc.code
        if bench.session is not None:
            bench.save_session({}, ended_utc=sc.iso_utc(time.time()), state="refused", reason=str(exc))
    finally:
        bench.cleanup()
    if bench.run_dir and os.path.isdir(bench.run_dir):
        executed = any(os.path.exists(os.path.join(bench.run_dir, "arms", arm, "executions.jsonl"))
                       for arm in bl.ARMS)
        if not args.no_analysis and executed:
            subprocess.run([cfg.py, os.path.join(HERE, "analyze.py"), bench.run_dir], stdin=subprocess.DEVNULL)
        if rc not in (EXIT_DONE,):
            log("resume with: %s --resume %s" % (sys.argv[0], bench.run_dir))
        print(bench.run_dir)
    return rc


if __name__ == "__main__":
    sys.exit(main())
