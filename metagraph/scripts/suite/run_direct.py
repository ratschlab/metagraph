#!/usr/bin/env python3
"""Direct runner of the suite: executes a frozen scenario manifest against one host arm.

    run_direct.py --manifest M --arm direct-local --base http://localhost:42010 --out RUN_DIR
                  [--areas M,P] [--limit N] [--ids ID,...] [--only ID,...] [--finish] [--sample-docker]
                  [--pause 1] [--checks auto|checks|stub] [--allow-stub] [--status-line]

Arms (DESIGN §1.2): `direct-local` runs on mex against the container's port (no front proxy);
`direct-remote` runs from the laptop through the public URL. Scenarios run in manifest order, one request in
flight, --pause seconds between requests; each scenario runs `repeat` times (repeat 0 cold, later warm; every
--warm-every-th scenario of the manifest gets one warm repeat). Per execution one line of executions.jsonl
(DESIGN §8.3) and the gzipped answer under raw/; run.json holds the run's identity (release, server_instance,
index_fp, the three capabilities documents, host load at start and end). The answers are judged by
`checks.evaluate` (checks.py of lane S1, or checks_stub.py, which is refused for a selection whose expectations
go beyond the HTTP status unless --allow-stub).

Timing: one keep-alive connection per invocation, reused while the server keeps it open; `wall_ms` runs from
sending the request to the last body byte on an established connection, and the connection set-up (TCP, TLS)
is recorded apart as `client.connect_ms` (null when the request reused the connection).

Paired steps (run_suite.sh --paired): `--only IDS` executes just these scenarios of the selection (the three
capabilities documents only when the arm has none yet, no closing GET); `--finish` executes nothing and takes the
closing GET /capabilities. Either way run.json's `planned` / `state` cover the whole selection: `done` when every
planned execution is recorded, else `partial`.

Resumable: a second invocation with the same --out skips the executions already recorded (same manifest
required). SIGINT / SIGTERM stop cleanly after the request in flight (a second SIGINT aborts it); run.json and
state.json then say `interrupted` (exit 130; exit 3 when --max-requests stopped the run). --max-requests counts
every request of the invocation, the closing GET /capabilities included. Stdlib only: runs with python3.12 on
mex.
"""

import argparse
import inspect
import http.client
import json
import os
import platform
import re
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

DEFAULT_TIMEOUTS = dict(sc.TIMING_CLASS_TIMEOUT_S)
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
    "User-Agent": "metagraph-suite-run-direct/1",
}
STATUS_PREFIX = "RUN_DIRECT_STATUS "


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
    """(module, name, import error): lane S1's checks.py unless --checks stub, else (or, with auto, when it
    cannot be imported) the stub; the import error says why auto fell back."""
    error = None
    if mode != "stub":
        try:
            import checks  # noqa: F401
            return checks, "checks", None
        except ImportError as exc:
            if mode == "checks":
                raise
            error = "%s: %s" % (type(exc).__name__, exc)
    import checks_stub
    return checks_stub, "stub", error


def accepts_keyword(fn, name):
    """Whether |fn| takes the keyword |name| (a vendored checks.evaluate of before may not take `client`)."""
    try:
        params = inspect.signature(fn).parameters
    except (TypeError, ValueError):
        return False
    return name in params or any(x.kind is inspect.Parameter.VAR_KEYWORD for x in params.values())


def time_limited_of(checks, answer):
    """What made an answer time-limited (checks.time_limited), '' when nothing did, None when the checks module
    cannot tell (the stub)."""
    fn = getattr(checks, "time_limited", None)
    if fn is None:
        return None
    try:
        return fn(answer) or ""
    except Exception:  # never take the run down over a note
        return None


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


# Errors of a request sent on a kept-alive connection that the server had closed in the meantime: the request
# never reached the application, so it is sent once more on a fresh connection.
STALE_CONNECTION_ERRORS = (http.client.RemoteDisconnected, BrokenPipeError, ConnectionResetError,
                           ConnectionAbortedError, http.client.CannotSendRequest)
# A closed idle connection fails at once (about a round trip); one that fails later may have carried a request
# the server was already processing, which is not sent again (review N10).
STALE_RESEND_MAX_S = 1.0


class Client(object):
    """One keep-alive HTTP(S) connection to the base URL, opened on demand and reused while the server keeps it
    open (it closes it after an answer with `Connection: close`)."""

    def __init__(self, base):
        u = urllib.parse.urlsplit(base)
        self.https = u.scheme == "https"
        self.host = u.hostname
        self.port = u.port or (443 if self.https else 80)
        self.prefix = u.path.rstrip("/")
        self.context = ssl.create_default_context() if self.https else None
        self.conn = None

    def close(self):
        if self.conn is not None:
            try:
                self.conn.close()
            except OSError:
                pass
        self.conn = None

    def _open(self, timeout):
        self.close()
        if self.https:
            self.conn = http.client.HTTPSConnection(self.host, self.port, timeout=timeout, context=self.context)
        else:
            self.conn = http.client.HTTPConnection(self.host, self.port, timeout=timeout)
        self.conn.connect()

    def call(self, method, path, body, headers, timeout):
        """One request. Returns {http, raw (decoded body bytes), headers, wall_ms, connect_ms, ttfb_ms,
        bytes_wire, reused, attempts, error}; http is 0 and error (code, text) set when no HTTP status arrived.

        wall_ms runs from sending the request to the last body byte, ttfb_ms to the status line, both on an
        established connection; connect_ms is the set-up of a fresh connection (null when reused). A request on
        a reused connection the server had closed meanwhile is sent once more on a fresh one (attempts 2), when
        it failed within STALE_RESEND_MAX_S of being sent; a later failure is recorded as a transport error.
        """
        result = {"http": 0, "raw": b"", "headers": {}, "wall_ms": None, "connect_ms": None, "ttfb_ms": None,
                  "bytes_wire": 0, "reused": False, "attempts": 0, "error": None}
        for attempt in (0, 1):
            result["attempts"] = attempt + 1
            reused = self.conn is not None and self.conn.sock is not None
            result["reused"] = reused
            t0 = time.perf_counter()
            try:
                if reused:
                    self.conn.timeout = timeout
                    self.conn.sock.settimeout(timeout)
                else:
                    self._open(timeout)
                    result["connect_ms"] = (time.perf_counter() - t0) * 1000.0
                t0 = time.perf_counter()
                self.conn.request(method, self.prefix + path, body=body, headers=headers)
                resp = self.conn.getresponse()
                result["ttfb_ms"] = (time.perf_counter() - t0) * 1000.0
                wire = resp.read()
                result["wall_ms"] = (time.perf_counter() - t0) * 1000.0
            except STALE_CONNECTION_ERRORS as exc:
                self.close()
                if reused and attempt == 0 and time.perf_counter() - t0 <= STALE_RESEND_MAX_S:
                    continue
                result["wall_ms"] = (time.perf_counter() - t0) * 1000.0
                result["error"] = ("transport", "%s: %s" % (type(exc).__name__, exc))
                return result
            except socket.timeout as exc:
                self.close()
                result["wall_ms"] = (time.perf_counter() - t0) * 1000.0
                result["error"] = ("timeout", "no answer within %.0f s: %s" % (timeout, exc))
                return result
            except (OSError, http.client.HTTPException) as exc:
                self.close()
                result["wall_ms"] = (time.perf_counter() - t0) * 1000.0
                result["error"] = ("transport", "%s: %s" % (type(exc).__name__, exc))
                return result
            result["http"] = resp.status
            result["headers"] = dict((k.lower(), v) for k, v in resp.getheaders())
            result["bytes_wire"] = len(wire)
            if resp.will_close:
                self.close()
            try:
                result["raw"] = sc.decode_body(wire, result["headers"].get("content-encoding"))
            except (ValueError, OSError, EOFError) as exc:
                result["raw"] = wire
                result["error"] = ("decode", "%s: %s" % (type(exc).__name__, exc))
            return result
        return result


def http_call(base, method, path, body, headers, timeout):
    """One request on a fresh connection (closed afterwards); see Client.call."""
    client = Client(base)
    try:
        return client.call(method, path, body, headers, timeout)
    finally:
        client.close()


# an id field's pattern (SPEC-labeled-traversal-core §5: attempt_id, budget_id, locus_id)
_ID_UNSAFE = re.compile(r"[^A-Za-z0-9._:-]")


def run_unique_request(scenario, run_id):
    """The scenario's request with each top-level field its `run_unique` lists made unique per run (A-10): '-' and
    the run id appended (the run directory's name, which the arms of a run share: within a run the first arm's
    request runs, the later ones meet its retained id's 409), cut to the id pattern's 128 characters. The manifest
    keeps the base value, so it stays the same for every run; without a run id or a run_unique the request as is."""
    request = scenario.get("request")
    fields = scenario.get("run_unique") or []
    if not run_id or not fields or not isinstance(request, dict):
        return request
    out = dict(request)
    suffix = _ID_UNSAFE.sub("_", str(run_id))
    for f in fields:
        if isinstance(out.get(f), str):
            out[f] = (out[f] + "-" + suffix)[:128]
    return out


def build_request(scenario, run_id=None):
    """(method, path, body bytes or None, headers, the request as recorded) of a scenario.

    The body is the compact JSON of `request` (the exact host body, DESIGN §8.1; its `run_unique` fields made unique
    per run with |run_id|: run_unique_request), or `body_raw` verbatim for protocol cases; `headers` in the scenario
    override the defaults (a null value removes one). A GET sends no body unless `body_raw` says so.
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
        request = run_unique_request(scenario, run_id)
        body = json.dumps(request, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
        recorded = {"body": request}
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
    def __init__(self, args, manifest, checks, checks_name, checks_error=None):
        self.args = args
        self.manifest = manifest
        self.checks = checks
        self.checks_name = checks_name
        self.evaluate_takes_client = accepts_keyword(getattr(checks, "evaluate", None), "client")
        self.checks_error = checks_error
        self.out = os.path.abspath(args.out)
        self.raw_dir = os.path.join(self.out, "raw")
        self.exec_path = os.path.join(self.out, "executions.jsonl")
        self.run_path = os.path.join(self.out, "run.json")
        self.state_path = os.path.join(self.out, "state.json")
        self.requests = 0
        self.executed = 0
        self.sampler = None
        self.client = Client(args.base)
        self.timeouts = {"fast": args.timeout_fast, "slow": args.timeout_slow, "long": args.timeout_long}
        # Which capabilities GETs this invocation makes: a paired step (--only) takes the three documents only
        # when the arm has none yet and no closing GET; --finish takes only the closing GET.
        self.closing_get = not args.only
        self.start_gets = not args.finish

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
            "skipped": skipped, "connection": "keep-alive; wall_ms excludes the connection set-up (connect_ms)",
            "options": {"areas": self.args.areas, "limit": self.args.limit, "ids": self.args.ids,
                        "pause_s": self.args.pause, "warm_every": self.args.warm_every,
                        "timeouts_s": self.timeouts, "remote_timeout_cap_s": REMOTE_TIMEOUT_CAP_S,
                        "max_requests": self.args.max_requests, "sample_docker": self.args.sample_docker,
                        "allow_stub": self.args.allow_stub}})
        if self.checks_error:
            self.run["checks_import_error"] = self.checks_error
        if not ok:
            warning = ("manifest content_sha256 %s does not match the canonical JSON of its scenarios (%s)"
                       % (self.manifest.get("content_sha256"), computed))
            if warning not in self.run.get("warnings", []):
                self.run.setdefault("warnings", []).append(warning)
        self.session = {"started_utc": sc.iso_utc(time.time()), "host_load_start": host_load(), "state": "running"}
        if self.args.only:
            self.session["only"] = self.args.only
        if self.args.finish:
            self.session["finish"] = True
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
        """Whether one more request may go out: no stop asked and --max-requests not reached, keeping one
        request back for the closing GET /capabilities when this invocation makes one."""
        if force:
            return True
        if STOP.flag:
            return False
        reserve = 1 if self.closing_get else 0
        if self.args.max_requests and self.requests + reserve >= self.args.max_requests:
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
            res = self.client.call("GET", path, None, {"Accept-Encoding": "gzip, deflate"}, CAPABILITY_TIMEOUT_S)
            self.requests += res["attempts"]
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
        method, path, body, headers, recorded = build_request(scenario, self.run["run_id"])
        timeout = timeout_for(scenario, self.args.arm, self.timeouts)
        started = time.time()
        res = self.client.call(method, path, body, headers, timeout)
        ended = time.time()
        self.requests += res["attempts"]
        if res["http"]:
            answer = sc.normalise_direct(res["http"], res["raw"], res["headers"])
        else:
            answer = sc.transport_failure(res["error"][0], res["error"][1])
        route = scenario.get("route", "%s %s" % (method, path))
        try:
            kw = {"previous": previous}
            if self.evaluate_takes_client:
                # the client's wall time: budget_respected judges an allowed 503 deadline by it
                kw["client"] = {"wall_ms": res["wall_ms"]}
            checks = list(self.checks.evaluate(scenario, answer, self.args.arm, **kw))
            # lane S1's rule (checks.scenario_outcome: property checks apart from the generic ones), else the
            # stub's own (suite_common.outcome_from_checks)
            outcome = getattr(self.checks, "scenario_outcome", sc.outcome_from_checks)(checks)
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
            "time_limited": time_limited_of(self.checks, answer),
            "error": answer["error"],
            "client": {"connect_ms": res["connect_ms"], "reused": res["reused"], "attempts": res["attempts"],
                       "ttfb_ms": res["ttfb_ms"], "bytes_wire": res["bytes_wire"],
                       "bytes_body": len(res["raw"]), "content_encoding": res["headers"].get("content-encoding"),
                       "timeout_s": timeout},
        }
        if after_resume:
            line["after_resume"] = True
        sc.append_jsonl(self.exec_path, line)
        self.executed += 1
        return line

    def previous_entry(self, line):
        """{digest, time_limited} of a recorded execution, for checks.evaluate's `previous`."""
        tl = line.get("time_limited")
        if tl is None and line.get("raw"):
            try:
                tl = time_limited_of(self.checks, (sc.read_gz_json(os.path.join(self.out, line["raw"])) or {}).get("answer"))
            except (OSError, ValueError):
                tl = None
        return {"digest": line.get("answer_digest"), "time_limited": tl or ""}

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
    def recorded_keys(self):
        rows, _ = sc.read_jsonl(self.exec_path)
        return set(sc.execution_key(r) for r in rows if r.get("arm") == self.args.arm)

    def final_state(self, planned):
        """done when every planned execution of the selection is recorded, else partial."""
        done = self.recorded_keys()
        missing = [p for p in planned if (p[1].get("id"), self.args.arm, p[2]) not in done]
        return ("done" if not missing else "partial"), len(planned) - len(missing)

    def run_all(self):
        os.makedirs(self.raw_dir, exist_ok=True)
        selected = sc.select_scenarios(self.manifest, self.args.areas, self.args.limit, self.args.ids)
        planned, skipped = sc.plan_executions(selected, self.args.arm, self.args.warm_every)
        if self.checks_name == "stub" and not self.args.allow_stub:
            names = sc.beyond_status([s for _, s in selected])
            if names:
                raise SystemExit("run_direct: the checks stub cannot evaluate %s (expected by the selected "
                                 "scenarios)%s; install lane S1's checks.py, or pass --allow-stub for a run that "
                                 "checks the HTTP status only"
                                 % (", ".join(names[:8]), (" and checks.py did not import: %s" % self.checks_error)
                                    if self.checks_error else ""))
        self.open_run(len(planned), skipped)
        rows, bad = sc.read_jsonl(self.exec_path)
        done = set(sc.execution_key(r) for r in rows)
        # The earlier repeats of each scenario so far, in repeat order, as {digest, time_limited}: checks.evaluate's
        # `previous` (lane S1's determinism check compares a repeat with the earlier ones, and a time-limited earlier
        # repeat need not repeat: SPEC-labeled-traversal-core §6.8, SPEC-pattern-search §7.9). A line of a runner
        # before the field gets it from its raw answer.
        digests = {}
        for r in sorted(rows, key=lambda r: (r.get("repeat") or 0)):
            if r.get("arm") == self.args.arm and r.get("answer_digest"):
                digests.setdefault(r.get("scenario"), []).append(self.previous_entry(r))
        if bad:
            self.run.setdefault("warnings", []).append("%d unparsable line(s) in executions.jsonl ignored" % bad)
        # The executions this invocation is for: the whole selection, a paired step's scenarios (--only), or
        # none (--finish); those already recorded are skipped.
        scope = planned
        if self.args.finish:
            scope = []
        elif self.args.only:
            only = set(self.args.only)
            scope = [p for p in planned if p[1].get("id") in only]
        todo = [p for p in scope if (p[1].get("id"), self.args.arm, p[2]) not in done]
        ndone = sum(1 for p in planned if (p[1].get("id"), self.args.arm, p[2]) in done)
        self.session["skipped_done"] = len(scope) - len(todo)
        self.save_run()
        self.write_state("running", len(planned), ndone)
        if self.args.sample_docker:
            self.start_sampler()
        state = "done"
        try:
            if self.start_gets and (not self.args.only or "capabilities" not in self.run):
                docs, statuses, ident = self.snapshot_capabilities()
                self.session["capabilities_http"] = statuses
                self.session["server_instance_start"] = ident["server_instance"]
                if "capabilities" not in self.run or not self.run.get("release"):
                    self.run["capabilities"] = docs
                    self.run["capabilities_http"] = statuses
                    self.run.update(ident)
                self.save_run()
            ran_in_session = set()
            for index, scenario, repeat, cache in todo:
                if not self.may_request():
                    break
                sid = scenario.get("id")
                self.write_state("running", len(planned), ndone, sid)
                after_resume = repeat > 0 and (sid, repeat - 1) not in ran_in_session
                line = self.execute(index, scenario, repeat, cache, after_resume, list(digests.get(sid, [])))
                digests.setdefault(sid, []).append(self.previous_entry(line))
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
            if self.closing_get and state != "failed" and STOP.signals < 2:
                try:
                    _, _, end_ident = self.snapshot_capabilities(at_end=True)
                except KeyboardInterrupt:
                    pass
            self.client.close()
            if state == "done":
                state, ndone = self.final_state(planned)
            self.session.update({"ended_utc": sc.iso_utc(time.time()), "state": state, "reason": STOP.reason,
                                 "executed": self.executed, "requests": self.requests,
                                 "host_load_end": host_load()})
            if self.closing_get:
                self.session["server_instance_end"] = end_ident.get("server_instance")
                self.run["server_instance_end"] = end_ident.get("server_instance")
            self.run["ended_utc"] = self.session["ended_utc"]
            self.run["host_load_end"] = self.session["host_load_end"]
            self.run["state"] = state
            starts = set(s.get("server_instance_start") for s in self.run["sessions"]
                         if s.get("server_instance_start"))
            ends = set(s.get("server_instance_end") for s in self.run["sessions"] if s.get("server_instance_end"))
            self.run["server_instance_changed"] = len(starts | ends) > 1
            self.save_run()
            done_now = len(self.recorded_keys())
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
    ap.add_argument("--only", default=None,
                    help="execute only these scenario ids of the selection (a paired step: the capabilities "
                         "documents only when the arm has none yet, no closing GET)")
    ap.add_argument("--finish", action="store_true",
                    help="execute nothing: the closing GET /capabilities and the selection's state (paired mode)")
    ap.add_argument("--pause", type=float, default=1.0, help="seconds between requests (1)")
    ap.add_argument("--warm-every", type=int, default=10,
                    help="one warm repeat for every Nth scenario of the manifest (10; 0: none)")
    ap.add_argument("--timeout-fast", type=float, default=DEFAULT_TIMEOUTS["fast"])
    ap.add_argument("--timeout-slow", type=float, default=DEFAULT_TIMEOUTS["slow"])
    ap.add_argument("--timeout-long", type=float, default=DEFAULT_TIMEOUTS["long"])
    ap.add_argument("--max-requests", type=int, default=0,
                    help="stop after this many requests in this invocation, the capabilities GETs (the closing "
                         "one included) and any resend on a fresh connection counted (0: no cap)")
    ap.add_argument("--sample-docker", action="store_true", help="run sampler.py alongside (on mex)")
    ap.add_argument("--sample-interval", type=float, default=2.0)
    ap.add_argument("--checks", choices=("auto", "checks", "stub"), default="auto",
                    help="auto: checks.py when it imports, else the stub")
    ap.add_argument("--allow-stub", action="store_true",
                    help="run with the stub even when the selection expects checks beyond the HTTP status")
    ap.add_argument("--run-id", default=None)
    ap.add_argument("--scripts-commit", default=None, help="the scripts' commit (the mex copy has no git)")
    ap.add_argument("--status-line", action="store_true",
                    help="print a last line '%s{json}' (state, executed, requests) for the orchestrator"
                         % STATUS_PREFIX)
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args(argv)
    args.areas = [a for a in args.areas.split(",") if a] if args.areas else None
    args.ids = [i for i in args.ids.split(",") if i] if args.ids else None
    args.only = [i for i in args.only.split(",") if i] if args.only else None
    if args.only and args.finish:
        ap.error("--only and --finish exclude each other")
    return args


def main(argv=None):
    args = parse_args(argv)
    manifest = sc.read_json(args.manifest)
    if not isinstance(manifest, dict) or "scenarios" not in manifest:
        raise SystemExit("run_direct: %s is not a suite manifest" % args.manifest)
    checks, checks_name, checks_error = load_checks(args.checks)
    reset_stop()
    install_signals()
    runner = Runner(args, manifest, checks, checks_name, checks_error)
    state = runner.run_all()
    if not args.quiet:
        sys.stdout.write("[%s] %s: %d executed, %d requests, run dir %s\n"
                         % (args.arm, state, runner.executed, runner.requests, runner.out))
    if args.status_line:
        sys.stdout.write(STATUS_PREFIX + json.dumps({"state": state, "executed": runner.executed,
                                                     "requests": runner.requests, "reason": STOP.reason},
                                                    sort_keys=True) + "\n")
    sys.stdout.flush()
    if state == "interrupted" and (STOP.reason or "").startswith("max_requests"):
        return EXIT_REQUEST_CAP
    return {"done": 0, "partial": 0, "interrupted": 130}.get(state, 1)


if __name__ == "__main__":
    sys.exit(main())
