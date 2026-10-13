#!/usr/bin/env python3
"""Direct runner of the suite: executes a frozen scenario manifest against one host arm.

    run_direct.py --manifest M --arm direct-local --base http://localhost:42010 --out RUN_DIR
                  [--areas M,P] [--limit N] [--ids ID,...] [--sample-docker] [--pause 1]

Arms (DESIGN §1.2): `direct-local` runs on mex against the container's port (no front proxy);
`direct-remote` runs from the laptop through the public URL. Scenarios run in manifest order, one request in
flight, --pause seconds between requests; each scenario runs `repeat` times (repeat 0 cold, later warm; every
--warm-every-th scenario of the manifest gets one warm repeat). Per execution one line of executions.jsonl
(DESIGN §8.3) and the gzipped answer under raw/; run.json holds the run's identity (release, server_instance,
index_fp, the three capabilities documents, host load at start and end). The answers are judged by
`checks.evaluate` (checks.py of lane S1, or checks_stub.py until it exists).

Resumable: a second invocation with the same --out skips the executions already recorded (same manifest
required). SIGINT / SIGTERM stop cleanly after the request in flight (a second SIGINT aborts it); run.json and
state.json then say `interrupted` (exit 130; exit 3 when --max-requests stopped the run). Stdlib only: runs with python3.12 on mex.
"""

import argparse
import http.client
import json
import os
import platform
import signal
import socket
import ssl
import subprocess
import sys
import time
import traceback
import urllib.parse

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import suite_common as sc  # noqa: E402

DEFAULT_TIMEOUTS = {"fast": 120.0, "slow": 600.0, "long": 1800.0}
# Exit codes: 0 done, 130 stopped by a signal, EXIT_REQUEST_CAP stopped at --max-requests (a planned stop: the
# orchestrator goes on to the report), 1 failed.
EXIT_REQUEST_CAP = 3
# The front proxy cuts at ~120 s; a little more lets its own answer (a 504) arrive instead of a client timeout.
REMOTE_TIMEOUT_CAP_S = 125.0
CAPABILITY_PATHS = ("/capabilities", "/pattern/capabilities", "/traverse/capabilities")
CAPABILITY_TIMEOUT_S = 60.0
DEFAULT_HEADERS = {
    "Content-Type": "application/json",
    "Accept": "application/json",
    "Accept-Encoding": "gzip, deflate",
    "Connection": "close",
    "User-Agent": "metagraph-suite-run-direct/1",
}


class StopState(object):
    def __init__(self):
        self.flag = False
        self.reason = None
        self.signals = 0


STOP = StopState()


def reset_stop():
    """A fresh stop state (the tests run several runs in one process)."""
    STOP.flag, STOP.reason, STOP.signals = False, None, 0


def _on_signal(signum, frame):
    STOP.signals += 1
    STOP.flag = True
    if STOP.reason is None:
        STOP.reason = "signal %d" % signum
    if STOP.signals >= 2 and signum == signal.SIGINT:
        raise KeyboardInterrupt()


def install_signals():
    signal.signal(signal.SIGINT, _on_signal)
    signal.signal(signal.SIGTERM, _on_signal)
    # Under nohup SIGHUP is ignored on purpose (the run outlives the ssh session): keep it ignored then.
    if hasattr(signal, "SIGHUP") and signal.getsignal(signal.SIGHUP) != signal.SIG_IGN:
        signal.signal(signal.SIGHUP, _on_signal)


def load_checks(mode):
    """(module, name): lane S1's checks.py unless --checks stub, else (or when it is absent) the stub."""
    if mode != "stub":
        try:
            import checks  # noqa: F401
            return checks, "checks"
        except ImportError:
            if mode == "checks":
                raise
    import checks_stub
    return checks_stub, "stub"


def host_load():
    try:
        with open("/proc/loadavg") as fh:
            return [float(x) for x in fh.read().split()[:3]]
    except (OSError, ValueError):
        try:
            return list(os.getloadavg())
        except (OSError, AttributeError):
            return None


def pid_alive(pid):
    try:
        os.kill(int(pid), 0)
    except (ProcessLookupError, ValueError, TypeError):
        return False
    except PermissionError:
        return True
    return True


def scripts_commit(explicit):
    if explicit:
        return explicit
    try:
        out = subprocess.run(["git", "-C", HERE, "rev-parse", "HEAD"], stdout=subprocess.PIPE,
                             stderr=subprocess.DEVNULL, universal_newlines=True, timeout=10)
        return out.stdout.strip() or None
    except (OSError, subprocess.SubprocessError):
        return None


def http_call(base, method, path, body, headers, timeout):
    """One request on a fresh connection.

    Returns {http, raw (decoded body bytes), headers, wall_ms, connect_ms, ttfb_ms, bytes_wire, error}; http is
    0 and error (code, text) set when no HTTP status arrived. wall_ms runs from before the connect to the last
    body byte; ttfb_ms to the status line.
    """
    u = urllib.parse.urlsplit(base)
    prefix = u.path.rstrip("/")
    full_path = prefix + path
    if u.scheme == "https":
        conn = http.client.HTTPSConnection(u.hostname, u.port or 443, timeout=timeout,
                                           context=ssl.create_default_context())
    else:
        conn = http.client.HTTPConnection(u.hostname, u.port or 80, timeout=timeout)
    result = {"http": 0, "raw": b"", "headers": {}, "wall_ms": None, "connect_ms": None, "ttfb_ms": None,
              "bytes_wire": 0, "error": None}
    t0 = time.perf_counter()
    try:
        conn.connect()
        result["connect_ms"] = (time.perf_counter() - t0) * 1000.0
        conn.request(method, full_path, body=body, headers=headers)
        resp = conn.getresponse()
        result["ttfb_ms"] = (time.perf_counter() - t0) * 1000.0
        wire = resp.read()
        result["wall_ms"] = (time.perf_counter() - t0) * 1000.0
        result["http"] = resp.status
        result["headers"] = dict((k.lower(), v) for k, v in resp.getheaders())
        result["bytes_wire"] = len(wire)
        try:
            result["raw"] = sc.decode_body(wire, result["headers"].get("content-encoding"))
        except (ValueError, OSError, EOFError) as exc:
            result["raw"] = wire
            result["error"] = ("decode", "%s: %s" % (type(exc).__name__, exc))
    except socket.timeout as exc:
        result["wall_ms"] = (time.perf_counter() - t0) * 1000.0
        result["error"] = ("timeout", "no answer within %.0f s: %s" % (timeout, exc))
    except (OSError, http.client.HTTPException) as exc:
        result["wall_ms"] = (time.perf_counter() - t0) * 1000.0
        result["error"] = ("transport", "%s: %s" % (type(exc).__name__, exc))
    finally:
        conn.close()
    return result


def build_request(scenario):
    """(method, path, body bytes or None, headers, the request as recorded) of a scenario.

    The body is the compact JSON of `request` (the exact host body, DESIGN §8.1), or `body_raw` verbatim for
    protocol cases; `headers` in the scenario override the defaults (a null value removes one). A GET sends no
    body unless `body_raw` says so.
    """
    method, path = sc.parse_route(scenario.get("route", "GET /capabilities"))
    headers = dict(DEFAULT_HEADERS)
    for name, value in (scenario.get("headers") or {}).items():
        for existing in list(headers):
            if existing.lower() == name.lower():
                del headers[existing]
        if value is not None:
            headers[name] = value
    if scenario.get("body_raw") is not None:
        raw = scenario["body_raw"]
        body = raw.encode("utf-8") if isinstance(raw, str) else bytes(raw)
        recorded = {"body_raw": raw if isinstance(raw, str) else None}
    elif method == "GET" and scenario.get("request") is None:
        body = None
        recorded = {"body": None}
    else:
        body = json.dumps(scenario.get("request"), separators=(",", ":"), ensure_ascii=False).encode("utf-8")
        recorded = {"body": scenario.get("request")}
    if body is None:
        headers.pop("Content-Type", None)
    recorded.update({"method": method, "path": path, "headers": headers})
    return method, path, body, headers, recorded


def timeout_for(scenario, arm, timeouts):
    t = timeouts.get(scenario.get("timing_class") or "fast", timeouts["fast"])
    if arm == "direct-remote":
        t = min(t, REMOTE_TIMEOUT_CAP_S)
    return t


class Runner(object):
    def __init__(self, args, manifest, checks, checks_name):
        self.args = args
        self.manifest = manifest
        self.checks = checks
        self.checks_name = checks_name
        self.out = os.path.abspath(args.out)
        self.raw_dir = os.path.join(self.out, "raw")
        self.exec_path = os.path.join(self.out, "executions.jsonl")
        self.run_path = os.path.join(self.out, "run.json")
        self.state_path = os.path.join(self.out, "state.json")
        self.requests = 0
        self.executed = 0
        self.sampler = None
        self.timeouts = {"fast": args.timeout_fast, "slow": args.timeout_slow, "long": args.timeout_long}

    # -- bookkeeping ------------------------------------------------------------------------------------------
    def default_run_id(self):
        parent = os.path.basename(os.path.dirname(self.out))
        if parent == "arms":
            return os.path.basename(os.path.dirname(os.path.dirname(self.out)))
        return os.path.basename(self.out)

    def write_state(self, state, planned, done, current=None):
        sc.write_json_atomic(self.state_path, {
            "run_id": self.run["run_id"], "arm": self.args.arm, "pid": os.getpid(), "state": state,
            "reason": STOP.reason, "planned": planned, "done": done, "executed_session": self.executed,
            "requests_session": self.requests, "current": current, "updated_utc": sc.iso_utc(time.time())})

    def open_run(self, planned_count, skipped):
        ok, computed = sc.verify_manifest(self.manifest)
        existing = sc.read_json(self.run_path)
        manifest_info = {"path": os.path.abspath(self.args.manifest),
                         "manifest_version": self.manifest.get("manifest_version"),
                         "content_sha256": self.manifest.get("content_sha256"), "computed_sha256": computed,
                         "sha_verified": ok, "index": self.manifest.get("index"),
                         "scenarios_total": len(self.manifest.get("scenarios", []))}
        state = sc.read_json(self.state_path) or {}
        if (state.get("state") == "running" and state.get("pid") not in (None, os.getpid())
                and existing and existing.get("runner_host") == platform.node() and pid_alive(state["pid"])):
            raise SystemExit("run_direct: %s is in use by a runner still alive (pid %s); stop it or wait"
                             % (self.out, state["pid"]))
        if existing:
            if existing.get("arm") != self.args.arm:
                raise SystemExit("run_direct: %s holds the arm %s, not %s; use another --out"
                                 % (self.out, existing.get("arm"), self.args.arm))
            if existing.get("manifest", {}).get("computed_sha256") != computed:
                raise SystemExit("run_direct: %s was run with another manifest (sha %s); use another --out"
                                 % (self.out, existing.get("manifest", {}).get("computed_sha256")))
            self.run = existing
        else:
            self.run = {"run_id": self.args.run_id or self.default_run_id(), "arm": self.args.arm,
                        "base": self.args.base, "manifest": manifest_info, "sessions": [],
                        "started_utc": sc.iso_utc(time.time())}
        self.run.update({
            "scripts_commit": scripts_commit(self.args.scripts_commit), "checks_module": self.checks_name,
            "runner_host": platform.node(), "python": platform.python_version(), "planned": planned_count,
            "skipped": skipped,
            "options": {"areas": self.args.areas, "limit": self.args.limit, "ids": self.args.ids,
                        "pause_s": self.args.pause, "warm_every": self.args.warm_every,
                        "timeouts_s": self.timeouts, "remote_timeout_cap_s": REMOTE_TIMEOUT_CAP_S,
                        "max_requests": self.args.max_requests, "sample_docker": self.args.sample_docker}})
        if not ok:
            self.run.setdefault("warnings", []).append(
                "manifest content_sha256 %s does not match the canonical JSON of its scenarios (%s)"
                % (self.manifest.get("content_sha256"), computed))
        self.session = {"started_utc": sc.iso_utc(time.time()), "host_load_start": host_load(), "state": "running"}
        self.run["sessions"].append(self.session)
        if "host_load_start" not in self.run:
            self.run["host_load_start"] = self.session["host_load_start"]
        # The run directory keeps the manifest it ran, so a later report or --compare needs nothing else.
        copy = os.path.join(self.out, "manifest.json")
        if not os.path.exists(copy):
            sc.write_json_atomic(copy, self.manifest)

    def save_run(self):
        sc.write_json_atomic(self.run_path, self.run)

    # -- requests ---------------------------------------------------------------------------------------------
    def may_request(self, force=False):
        if force:
            return True
        if STOP.flag:
            return False
        if self.args.max_requests and self.requests >= self.args.max_requests:
            STOP.flag = True
            STOP.reason = "max_requests %d reached" % self.args.max_requests
            return False
        return True

    def pause(self):
        end = time.time() + self.args.pause
        while not STOP.flag and time.time() < end:
            time.sleep(min(0.1, max(0.0, end - time.time())))

    def snapshot_capabilities(self, at_end=False):
        """The three capabilities documents at the start; at the end GET /capabilities only, whatever stopped
        the run (it states whether the server restarted during the run)."""
        docs, statuses = {}, {}
        paths = CAPABILITY_PATHS[:1] if at_end else CAPABILITY_PATHS
        for path in paths:
            if not self.may_request(force=at_end):
                break
            res = http_call(self.args.base, "GET", path, None, {"Accept-Encoding": "gzip, deflate",
                                                                  "Connection": "close"}, CAPABILITY_TIMEOUT_S)
            self.requests += 1
            statuses[path] = res["http"]
            if res["http"]:
                try:
                    docs[path] = json.loads(res["raw"].decode("utf-8", "replace"))
                except ValueError:
                    docs[path] = res["raw"].decode("utf-8", "replace")[:2000]
            else:
                docs[path] = {"transport_error": list(res["error"])}
            if not at_end:
                self.pause()
        ident = sc.identity_from_capabilities(docs)
        return docs, statuses, ident

    def execute(self, index, scenario, repeat, cache, after_resume, previous):
        sid = scenario.get("id")
        method, path, body, headers, recorded = build_request(scenario)
        timeout = timeout_for(scenario, self.args.arm, self.timeouts)
        started = time.time()
        res = http_call(self.args.base, method, path, body, headers, timeout)
        ended = time.time()
        self.requests += 1
        if res["http"]:
            answer = sc.normalise_direct(res["http"], res["raw"], res["headers"])
        else:
            answer = sc.transport_failure(res["error"][0], res["error"][1])
        route = scenario.get("route", "%s %s" % (method, path))
        try:
            checks = list(self.checks.evaluate(scenario, answer, self.args.arm, previous=previous))
            outcome = sc.outcome_from_checks(checks)
        except Exception as exc:  # a check that crashes is an error of this execution, never of the run
            checks = [{"check": "evaluate", "outcome": "fail",
                       "detail": "exception in checks.evaluate: %s: %s | %s"
                                 % (type(exc).__name__, exc, traceback.format_exc(limit=3).replace("\n", " "))}]
            outcome = "error"
        if not res["http"] or (res["error"] and res["error"][0] == "decode"):
            outcome = "error"
        raw_rel = sc.raw_name(sid, self.args.arm, repeat)
        sc.write_gz_json(os.path.join(self.out, raw_rel), {
            "scenario": sid, "arm": self.args.arm, "repeat": repeat, "started_utc": sc.iso_utc(started),
            "request": recorded, "response": {"http": res["http"], "headers": res["headers"],
                                              "bytes_wire": res["bytes_wire"], "bytes_body": len(res["raw"])},
            "answer": answer})
        line = {
            "run_id": self.run["run_id"], "scenario": sid, "area": scenario.get("area"), "route": route,
            "timing_class": scenario.get("timing_class"), "arm": self.args.arm, "repeat": repeat, "cache": cache,
            "started_utc": sc.iso_utc(started), "ended_utc": sc.iso_utc(ended), "wall_ms": res["wall_ms"],
            "http": answer["http"], "outcome": outcome, "checks": checks,
            "server": sc.extract_server(answer["body"]), "mcp_phases": None,
            "answer_digest": sc.answer_digest(answer, route), "raw": raw_rel, "job": None,
            "error": answer["error"],
            "client": {"connect_ms": res["connect_ms"], "ttfb_ms": res["ttfb_ms"], "bytes_wire": res["bytes_wire"],
                       "bytes_body": len(res["raw"]), "content_encoding": res["headers"].get("content-encoding"),
                       "timeout_s": timeout},
        }
        if after_resume:
            line["after_resume"] = True
        sc.append_jsonl(self.exec_path, line)
        self.executed += 1
        return line

    # -- sampler ----------------------------------------------------------------------------------------------
    def start_sampler(self):
        cmd = [sys.executable, os.path.join(HERE, "sampler.py"), "--out", os.path.join(self.out, "samples.jsonl"),
               "--interval", str(self.args.sample_interval), "--label", self.args.arm,
               "--parent-pid", str(os.getpid())]
        log = open(os.path.join(self.out, "sampler.log"), "a")
        # A session of its own: a terminal's Ctrl-C reaches the runner, which then stops the sampler itself.
        try:
            self.sampler = subprocess.Popen(cmd, stdout=log, stderr=log, stdin=subprocess.DEVNULL,
                                            start_new_session=True)
        finally:
            log.close()  # the sampler holds its own descriptor

    def stop_sampler(self):
        if self.sampler is None:
            return
        try:
            self.sampler.terminate()
            self.sampler.wait(timeout=15)
        except subprocess.TimeoutExpired:
            self.sampler.kill()
        except OSError:
            pass
        self.sampler = None

    # -- the run ----------------------------------------------------------------------------------------------
    def run_all(self):
        os.makedirs(self.raw_dir, exist_ok=True)
        selected = sc.select_scenarios(self.manifest, self.args.areas, self.args.limit, self.args.ids)
        planned, skipped = sc.plan_executions(selected, self.args.arm, self.args.warm_every)
        self.open_run(len(planned), skipped)
        rows, bad = sc.read_jsonl(self.exec_path)
        done = set(sc.execution_key(r) for r in rows)
        # The answer digests of each scenario's repeats so far, in repeat order: checks.evaluate's `previous`
        # (lane S1's determinism check compares a repeat with the earlier ones).
        digests = {}
        for r in sorted(rows, key=lambda r: (r.get("repeat") or 0)):
            if r.get("arm") == self.args.arm and r.get("answer_digest"):
                digests.setdefault(r.get("scenario"), []).append(r["answer_digest"])
        if bad:
            self.run.setdefault("warnings", []).append("%d unparsable line(s) in executions.jsonl ignored" % bad)
        todo = [p for p in planned if (p[1].get("id"), self.args.arm, p[2]) not in done]
        self.session["skipped_done"] = len(planned) - len(todo)
        self.save_run()
        self.write_state("running", len(planned), len(planned) - len(todo))
        if self.args.sample_docker:
            self.start_sampler()
        state = "done"
        try:
            docs, statuses, ident = self.snapshot_capabilities()
            self.session["capabilities_http"] = statuses
            self.session["server_instance_start"] = ident["server_instance"]
            if "capabilities" not in self.run or not self.run.get("release"):
                self.run["capabilities"] = docs
                self.run["capabilities_http"] = statuses
                self.run.update(ident)
            self.save_run()
            ran_in_session = set()
            ndone = len(planned) - len(todo)
            for index, scenario, repeat, cache in todo:
                if not self.may_request():
                    break
                sid = scenario.get("id")
                self.write_state("running", len(planned), ndone, sid)
                after_resume = repeat > 0 and (sid, repeat - 1) not in ran_in_session
                line = self.execute(index, scenario, repeat, cache, after_resume, list(digests.get(sid, [])))
                digests.setdefault(sid, []).append(line["answer_digest"])
                ran_in_session.add((sid, repeat))
                ndone += 1
                if not self.args.quiet:
                    sys.stdout.write("[%s] %d/%d %s r%d %s http=%s %.0f ms\n" % (
                        self.args.arm, ndone, len(planned), sid, repeat, line["outcome"], line["http"],
                        line["wall_ms"] or 0.0))
                    sys.stdout.flush()
                self.pause()
            if STOP.flag:
                state = "interrupted"
        except KeyboardInterrupt:
            state = "interrupted"
            STOP.reason = STOP.reason or "KeyboardInterrupt"
        except Exception:
            state = "failed"
            self.run.setdefault("warnings", []).append("runner failed: " + traceback.format_exc(limit=5))
            raise
        finally:
            end_ident = {"server_instance": None}
            if state != "failed" and STOP.signals < 2:
                try:
                    _, _, end_ident = self.snapshot_capabilities(at_end=True)
                except KeyboardInterrupt:
                    pass
            self.session.update({"ended_utc": sc.iso_utc(time.time()), "state": state, "reason": STOP.reason,
                                 "executed": self.executed, "requests": self.requests,
                                 "host_load_end": host_load(),
                                 "server_instance_end": end_ident.get("server_instance")})
            self.run["ended_utc"] = self.session["ended_utc"]
            self.run["host_load_end"] = self.session["host_load_end"]
            self.run["state"] = state
            starts = set(s.get("server_instance_start") for s in self.run["sessions"]
                         if s.get("server_instance_start"))
            ends = set(s.get("server_instance_end") for s in self.run["sessions"] if s.get("server_instance_end"))
            self.run["server_instance_end"] = end_ident.get("server_instance")
            self.run["server_instance_changed"] = len(starts | ends) > 1
            self.save_run()
            rows, _ = sc.read_jsonl(self.exec_path)
            done_now = len(set(sc.execution_key(r) for r in rows if r.get("arm") == self.args.arm))
            self.write_state(state, len(planned), min(done_now, len(planned)))
            self.stop_sampler()
        return state


def parse_args(argv):
    ap = argparse.ArgumentParser(description="Run a suite manifest against one direct arm (DESIGN §8.5).")
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--arm", required=True, choices=sc.DIRECT_ARMS)
    ap.add_argument("--base", required=True,
                    help="host base URL, e.g. http://localhost:42010 or "
                         "https://metagraph.ethz.ch:8080/api/refseq33m-experimental")
    ap.add_argument("--out", required=True, help="the arm's run directory (created; resumed when it exists)")
    ap.add_argument("--areas", default=None, help="comma-separated areas (default: all)")
    ap.add_argument("--ids", default=None, help="comma-separated scenario ids (default: all)")
    ap.add_argument("--limit", type=int, default=None, help="the first N selected scenarios only")
    ap.add_argument("--pause", type=float, default=1.0, help="seconds between requests (1)")
    ap.add_argument("--warm-every", type=int, default=10,
                    help="one warm repeat for every Nth scenario of the manifest (10; 0: none)")
    ap.add_argument("--timeout-fast", type=float, default=DEFAULT_TIMEOUTS["fast"])
    ap.add_argument("--timeout-slow", type=float, default=DEFAULT_TIMEOUTS["slow"])
    ap.add_argument("--timeout-long", type=float, default=DEFAULT_TIMEOUTS["long"])
    ap.add_argument("--max-requests", type=int, default=0,
                    help="stop after this many requests in this session, capabilities included (0: no cap)")
    ap.add_argument("--sample-docker", action="store_true", help="run sampler.py alongside (on mex)")
    ap.add_argument("--sample-interval", type=float, default=2.0)
    ap.add_argument("--checks", choices=("auto", "checks", "stub"), default="auto",
                    help="auto: checks.py when present, else the stub")
    ap.add_argument("--run-id", default=None)
    ap.add_argument("--scripts-commit", default=None, help="the scripts' commit (the mex copy has no git)")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args(argv)
    args.areas = [a for a in args.areas.split(",") if a] if args.areas else None
    args.ids = [i for i in args.ids.split(",") if i] if args.ids else None
    return args


def main(argv=None):
    args = parse_args(argv)
    manifest = sc.read_json(args.manifest)
    if not isinstance(manifest, dict) or "scenarios" not in manifest:
        raise SystemExit("run_direct: %s is not a suite manifest" % args.manifest)
    checks, checks_name = load_checks(args.checks)
    reset_stop()
    install_signals()
    runner = Runner(args, manifest, checks, checks_name)
    state = runner.run_all()
    if not args.quiet:
        sys.stdout.write("[%s] %s: %d executed, %d requests, run dir %s\n"
                         % (args.arm, state, runner.executed, runner.requests, runner.out))
    if state == "interrupted" and (STOP.reason or "").startswith("max_requests"):
        return EXIT_REQUEST_CAP
    return {"done": 0, "interrupted": 130}.get(state, 1)


if __name__ == "__main__":
    sys.exit(main())
